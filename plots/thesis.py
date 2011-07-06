# encoding: utf-8

import os
import math
from glob import glob
from array import array

import sqlite3 as sqlite
import ROOT
import analysis

from analysis.asym import AsymmetryGenerator, ScalarCounts, Polarizations
from analysis.histos2 import HistogramManager
from analysis.histos import shifted
from analysis.plots import pid_calibration
from analysis.runlists import long2_run6 as runlist
from analysis.runlists import final_runlist_run5 as runlist05
from analysis.runlists import transverse_run6 as transverse_runlist
from analysis.util import getRun, fillList, hadd_interactive
from analysis.plots import graphics

histDir         = '/Users/kocolosk/data/run6/spin2008/hist'
transHistDir    = '/Users/kocolosk/data/run6/spin2008/hist-transverse'
mcasymFile      = '/Users/kocolosk/data/run6/simu/mcasym_run6_thesis.root'
# mcasymFile      = '/Users/kocolosk/minbias_reweighted.root'
datamcFile      = '/Users/kocolosk/data/run6/spin2008/merged.cphist.root'

histDir05       = '/Users/kocolosk/data/run5/hist'

zbins = [
    # 0.075, 0.125, ## these two are biased b/c of the track pT cut
    0.2, 0.3, 0.45, 0.65, 1.0
]

ptbins = [2.0, 3.18, 4.56, 6.32, 8.8, 12.84]

def final_result_run5():
    ROOT.gStyle.SetErrorX(0)
    ROOT.gStyle.SetOptDate(0)
    
    f = ROOT.TFile('/Users/kocolosk/data/run5/final_result.root')
    hm = f.Get('final_minus')
    hp = f.Get('final_plus')
    
    ## systematic uncertainties x 10^-3 
    
    # ## trigger bias dropping MAX
    # trig_m_l = [12.00, 13.00, 16.00, 21.00,  9.00]
    # trig_m_h = [ 8.00, 10.00, 12.00, 13.00, 18.00]
    # trig_p_l = [ 8.00, 35.00, 27.00, 41.00, 10.00]
    # trig_p_h = [ 8.00, 14.00, 12.00, 31.00, 10.00]
    # 
    # ## trigger bias dropping MIN, MAX
    # trig_m_l = [12.00, 13.00, 16.00, 21.00,  9.00]
    # trig_m_h = [ 5.00,  8.00,  8.00, 10.00, 11.00]
    # trig_p_l = [ 5.00, 35.00, 27.00, 41.00, 10.00]
    # trig_p_h = [ 5.00, 14.00, 11.00, 12.00, 10.00]
    # 
    # ## trigger bias dropping MIN, MAX, P*
    # trig_m_l = [5.00,  5.00,  7.00,  5.00,  4.00]
    # trig_m_h = [5.00,  8.00,  7.00,  3.00, 11.00]
    # trig_p_l = [4.00, 13.00, 11.00, 21.00,  7.00]
    # trig_p_h = [4.00, 14.00, 11.00, 12.00,  4.00]
    # 
    # ## trig+reco dropping MIN, MAX, P*
    # trig_m_l = [5.30,  4.90,  6.70,  4.50,  3.40]
    # trig_m_h = [5.30,  9.10,  6.70,  3.00,  6.50]
    # trig_p_l = [4.20, 13.20,  9.30, 20.30, 13.70]
    # trig_p_h = [4.20, 14.30,  2.40,  9.10,  4.00]
    # 
    # ## reco dropping MIN, MAX, P*
    # # trig_m_l = [0.30,  1.30,  1.20,  4.70,  6.40]
    # # trig_m_h = [0.30,  1.50,  1.20,  4.70,  6.40]
    # # trig_p_l = [0.30,  1.10,  8.60,  3.80,  6.90]
    # # trig_p_h = [0.30,  1.10,  5.40,  3.80,  4.50]
    # 
    # ## trig+reco after rescaling, dropping MAX,MIN,P*
    # trig_m_l = [5.90,  8.30,  9.30,  7.20,  4.80]
    # trig_m_h = [2.70,  7.20,  3.40,  3.60,  5.70]
    # trig_p_l = [6.10, 12.80, 17.60, 20.90, 15.20]
    # trig_p_h = [6.10,  6.60, 10.10, 18.60,  6.20]
    # 
    # # asigma_m = [1.80, 1.80, 1.80, 1.80, 1.80]
    # # asigma_p = [1.80, 1.80, 1.80, 1.80, 1.80]
    # 
    # asigma_m = [2.91, 2.30, 3.74, 5.98, 12.07]
    # asigma_p = [1.34, 2.56, 5.96, 5.72, 15.11]
    # 
    # pid_m = [0.00, 0.00, 0.00, 0.00, 0.00]
    # pid_p = [0.00, 0.00, 0.00, 0.00, 0.00]
    # 
    # rellumi = [0.98, 0.98, 0.98, 0.98, 0.98]
    
    ## trig+reco after rescaling, M030-P030
    trig_m_h = [1.50,  2.30,  2.40,  4.10,  5.20]
    trig_m_l = [7.60,  10.1,  10.5,  9.70,  7.0]
    trig_p_h = [1.50,  2.1, 2.5, 4.1, 5.1]
    trig_p_l = [2.10,  14.90, 20.20, 22.80,  16.20]
    
    asigma_m = [2.91, 2.30, 3.74, 5.98, 12.07]
    asigma_p = [1.34, 2.56, 5.96, 5.72, 15.11]
    
    rellumi = [0.98, 0.98, 0.98, 0.98, 0.98]
    
    sum_quadrature = lambda *args: math.sqrt( sum([a**2 for a in args]) )/1E3
    
    syst_m = ROOT.TGraphAsymmErrors(5, 
        array('d', [hm.GetBinCenter(bin) for bin in range(1,6)]),
        array('d', [hm.GetBinContent(bin) for bin in range(1,6)]),
        array('d', [0.2 for i in range(5)]),      ## exl
        array('d', [0.2 for i in range(5)]),      ## exh
        array('d', map(sum_quadrature, trig_m_l, asigma_m, rellumi)),
        array('d', map(sum_quadrature, trig_m_h, asigma_m, rellumi))
    )
    
    syst_p = ROOT.TGraphAsymmErrors(5,
        array('d', [hp.GetBinCenter(bin) for bin in range(1,6)]),
        array('d', [hp.GetBinContent(bin) for bin in range(1,6)]),
        array('d', [0.2 for i in range(5)]),      ## exl
        array('d', [0.2 for i in range(5)]),      ## exh
        array('d', map(sum_quadrature, trig_p_l, asigma_p, rellumi)),
        array('d', map(sum_quadrature, trig_p_h, asigma_p, rellumi))
    )
    
    ## auto-generate LaTeX tables
    for i in range(5):
        print '%.2f - %.2f & %.1f & $\pm$ %.1f & -%.1f +%.1f &  %.1f & $\pm$ %.1f & -%.1f +%.1f\\\\'\
            % (ptbins[i], ptbins[i+1], 
            hm.GetBinContent(i+1)*1000, hm.GetBinError(i+1)*1000,
            syst_m.GetErrorYlow(i)*1000, syst_m.GetErrorYhigh(i)*1000,
            hp.GetBinContent(i+1)*1000, hp.GetBinError(i+1)*1000,
            syst_p.GetErrorYlow(i)*1000, syst_p.GetErrorYhigh(i)*1000)
    
    ## NLO curves
    from analysis.asym import theoryCurves
    nlo_m = {
        'STD': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_std, 
            analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
        'ZERO': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_zero, 
            analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
        # 'MAX': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_max, 
        #     analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
        'MIN': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_min, 
            analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),    
        'GSC': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_gsc, 
            analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph()
    }
    nlo_p = {
        'STD': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_std, 
            analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
        'ZERO': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_zero, 
            analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
        # 'MAX': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_max, 
        #     analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
        'MIN': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_min, 
            analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
        'GSC': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_gsc, 
            analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph()
    }
    
    
    nlo_m['DSSV'] = ROOT.TGraphErrors(15,
        array('d', [float(i) for i in range(1,16)]),
        array('d', [
            -2.685E-05, -4.188E-05, -8.716E-05, -1.562E-04, -1.841E-04,
            -1.420E-04, -5.366E-05,  6.018E-05,  1.274E-04,  1.235E-04,
             7.536E-06, -2.645E-04, -6.844E-04, -1.261E-03, -1.999E-03
        ]),
        array('d', [0 for i in range(15)]),
        array('d', [
            3.011E-05, 3.331E-05, 1.549E-04, 3.664E-04, 5.990E-04,
            8.074E-04, 9.963E-04, 1.168E-03, 1.345E-03, 1.513E-03,
            1.684E-03, 1.851E-03, 1.990E-03, 2.117E-03, 2.223E-03
        ]))
    
    nlo_p['DSSV'] = ROOT.TGraphErrors(15,
        array('d', [float(i) for i in range(1,16)]),
        array('d', [
            -3.800E-05, -3.658E-05,  7.296E-07,  1.634E-04,  6.061E-04,
             1.404E-03,  2.588E-03,  4.144E-03,  5.984E-03,  8.075E-03,
             1.031E-02,  1.271E-02,  1.528E-02,  1.779E-02,  2.039E-02
        ]),
        array('d', [0 for i in range(15)]),
        array('d', [
            1.821E-05, 8.788E-05, 3.054E-04, 6.487E-04, 1.038E-03, 
            1.406E-03, 1.758E-03, 2.090E-03, 2.417E-03, 2.745E-03, 
            3.080E-03, 3.430E-03, 3.750E-03, 4.090E-03, 4.390E-03
        ]))
    
    for nlo in (nlo_m, nlo_p):
        ## line style from pi0 paper
        # nlo['MAX'].SetLineStyle(2)
        nlo['MIN'].SetLineStyle(3)
        nlo['STD'].SetLineStyle(4)
        nlo['ZERO'].SetLineStyle(8)
        nlo['GSC'].SetLineStyle(7)
        
        nlo['ZERO'].SetLineColor(ROOT.kBlue)
        # nlo['MAX'].SetLineStyle(4)
        # nlo['MAX'].SetLineColor(ROOT.kRed)
        # nlo['MIN'].SetLineStyle(2)
        nlo['MIN'].SetLineColor(ROOT.kGreen)
        # nlo['GSC'].SetLineStyle(4)
        nlo['GSC'].SetLineColor(ROOT.kMagenta)
        nlo['DSSV'].SetFillColor(ROOT.kYellow)
        # nlo['DSSV'].SetLineStyle(5)
        # nlo['DSSV'].SetLineStyle(1)
        [gr.SetLineWidth(2) for gr in nlo.values()]
        nlo['DSSV'].SetLineWidth(3)
        # nlo['DSSV'].SetLineWidth(2)
    
    ## graphics setup
    for h in (hm,hp):
        # h.GetYaxis().SetRangeUser(-0.14, 0.09)
        h.GetYaxis().SetRangeUser(-0.12, 0.04)
        h.GetYaxis().SetNdivisions(507)
        h.GetYaxis().SetTitle('A_{LL}  ')
        h.GetYaxis().SetTitleSize(0.06)
        h.GetYaxis().SetTitleOffset(0.7)
        h.GetYaxis().SetLabelOffset(0.015)
        h.GetXaxis().SetTitle('p_{T}  [GeV/c]  ')
        h.GetXaxis().SetTitleSize(0.045)
        h.GetXaxis().SetTitleOffset(1.16)
    
    for g in (syst_m,syst_p):
        g.SetLineColor(ROOT.kRed)
        g.SetFillColor(15)
    
    line = ROOT.TLine(2.00, 0.0, 12.84, 0.0)
    line.SetLineStyle(2)
    
    # leg = ROOT.TLegend(0.17, 0.17, 0.52, 0.47)
    # leg = ROOT.TLegend(0.17, 0.17, 0.7, 0.3)
    leg = ROOT.TLegend(0.17, 0.17, 0.85, 0.35)
    leg.SetNColumns(2)
    leg.SetColumnSeparation(-0.2)
    leg.SetBorderSize(0)
    leg.AddEntry(nlo_m['STD'], 'GRSV STD', 'l')
    leg.AddEntry(nlo_m['DSSV'], 'DSSV', 'l')
    leg.AddEntry(nlo_m['MIN'], 'GRSV MIN', 'l')
    leg.AddEntry(nlo_m['GSC'], 'GS Set C', 'l')
    leg.AddEntry(nlo_m['ZERO'], 'GRSV ZERO', 'l')
    # leg.AddEntry(nlo_m['MAX'], 'MAX', 'l')
    # leg.AddEntry(nlo_m['DSSV'], 'DSSV', 'lf')
    
    latex = ROOT.TLatex()
    latex.SetTextSize(0.24)
    
    poltext = ROOT.TLatex()
    poltext.SetTextSize(0.0415)
    
    c = graphics.canvas2()
    
    for i in (1, 2):
        pad = c.cd(i)
        pad.SetTopMargin(0.03)
        pad.SetBottomMargin(0.12)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.02)
        
    pad1 = c.cd(1)
    hm.Draw('e1')
    # nlo_m['DSSV'].Draw('3')
    nlo_m['DSSV'].Draw('lx')
    [ g.Draw('lx') for g in nlo_m.values() ]
    line.Draw()
    # latex.DrawLatex(2.5, 0.05, '#pi^{-}')
    latex.DrawLatex(3.0, -0.06, '#pi^{-}')
    text = '#pm 9.4% polarization scale uncertainty not shown'
    # poltext.DrawLatex(2.5, -0.125, text)
    poltext.DrawLatex(2.9, -0.105, text)
    syst_m.Draw('2p')
    hm.Draw('e1 same')
    
    pad2 = c.cd(2)
    hp.Draw('e1')
    # nlo_p['DSSV'].Draw('3')
    nlo_p['DSSV'].Draw('lx')
    [ g.Draw('lx') for g in nlo_p.values() ]
    line.Draw()
    leg.Draw()
    # latex.DrawLatex(2.5, 0.05, '#pi^{+}')
    latex.DrawLatex(3.0, -0.06, '#pi^{+}')
    syst_p.Draw('2p')
    hp.Draw('e1 same')
    
    graphics.maybe_save(c)


def final_result_run6():
    """
    result plot showing asymmetries versus z
    """    
    asym_p = AsymmetryGenerator('asym_p', bins=zbins, key='z_away2')
    asym_m = AsymmetryGenerator('asym_m', bins=zbins, key='z_away2')
    
    db = sqlite.connect('/Users/kocolosk/data/analysis.db')
    dbc = db.cursor()
    polarizations = Polarizations.Final
    
    allFiles = glob(histDir + '/chargedPions_*.hist.root')
    for fname in allFiles[:]:
        run = getRun(fname)
        if run in runlist:
            print fname, run
            mgr = HistogramManager(ROOT.TFile(fname), keys=['z_away2'])
            uu,ud,du,dd = analysis.util.scaler_counts(dbc, run)
            pol = polarizations[analysis.util.fill(run)]
            asym_p.FillFromHistogramManager(mgr, 'jetpatch', 1, uu,ud,du,dd, \
                pol.py,pol.pb)
            asym_m.FillFromHistogramManager(mgr, 'jetpatch', -1, uu,ud,du,dd, \
                pol.py,pol.pb)
    
    ## fun with graphics
    ROOT.gStyle.SetErrorX(0.0)
    ROOT.gStyle.SetOptDate(0)
    
    c = graphics.canvas2()
    c.SetLogy(0)
    
    for i in (1, 2):
        pad = c.cd(i)
        pad.SetTopMargin(0.03)
        pad.SetBottomMargin(0.1)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.02)
    
    line = ROOT.TLine(zbins[0], 0.0, zbins[-1], 0.0)
    line.SetLineStyle(2)
    
    prelim = ROOT.TText()
    prelim.SetTextColor(44)
    prelim.SetTextAlign(21)
    
    latex = ROOT.TLatex()
    latex.SetTextSize(0.25)
    
    syst = systematic_uncertainties_run6()
    syst_m = ROOT.TGraphAsymmErrors(len(zbins))
    syst_p = ROOT.TGraphAsymmErrors(len(zbins))
    
    for g in (syst_m,syst_p):
        g.SetMarkerColor(15)
        g.SetFillColor(15)
        g.SetPoint(len(zbins), 1.0, 0)
        g.SetPointError(len(zbins), 0., 0., 0., 0.)
        g.GetXaxis().SetRangeUser(zbins[0], zbins[-1])
    
    hm = asym_m.GetAsymmetry('ll')
    hp = asym_p.GetAsymmetry('ll')
    for i in range(len(zbins)-1):
        z = (zbins[i] + zbins[i+1])/2
        syst_m.SetPoint(i, z, hm.GetBinContent(i+1))
        syst_m.SetPointError(i, 0.01, 0.01, syst['minus']['low'][i], syst['minus']['high'][i])
        syst_p.SetPoint(i, z, hp.GetBinContent(i+1))
        syst_p.SetPointError(i, 0.01, 0.01, syst['plus']['low'][i], syst['plus']['high'][i])
    
    [ g.GetXaxis().SetRangeUser(zbins[0], zbins[-1]) for g in (syst_m,syst_p) ]
    
    for h in (hm,hp, syst_m, syst_p):
        h.SetTitle('')
        h.GetXaxis().SetTitle('z')
        h.GetYaxis().SetRangeUser(-0.05, 0.08)
        h.GetYaxis().SetTitle('A_{LL}  ')
        h.GetYaxis().SetTitleSize(0.06)
        h.GetYaxis().SetTitleOffset(1.0)
        h.GetXaxis().SetTitle('z  ')
        h.GetXaxis().SetTitleSize(0.065)
        h.GetXaxis().SetTitleOffset(0.7)
    
    title = {
        'STD': 'GRSV-STD',
        'GSC': 'GS Set C',
        'DSSV': 'DSSV'
    }
    
    def make_graph(num, denom):
        g = ROOT.TGraphErrors(len(num))
        for i,n in enumerate(num):
            d = denom[i]
            val = n.y / d.y
            err = val*math.sqrt((n.sys/n.y)**2 + (d.sys/d.y)**2)
            g.SetPoint(i, n.x, val)
            g.SetPointError(i, 0, 0)
        return g
    
    nlo_m = {
        'STD': make_graph(analysis.deflorian.minus.std, 
            analysis.deflorian.minus.mrst),
        'DSSV': make_graph(analysis.deflorian.minus.dssv, 
            analysis.deflorian.minus.mrst),
        'GSC': make_graph(analysis.deflorian.minus.gsc, 
            analysis.deflorian.minus.mrst),
    }
    nlo_p = {
        'STD': make_graph(analysis.deflorian.plus.std, 
            analysis.deflorian.plus.mrst),
        'DSSV': make_graph(analysis.deflorian.plus.dssv, 
            analysis.deflorian.plus.mrst),
        'GSC': make_graph(analysis.deflorian.plus.gsc, 
            analysis.deflorian.plus.mrst),
    }
    
    for nlo in (nlo_m, nlo_p):
        nlo['GSC'].SetLineColor(ROOT.kMagenta)
        nlo['GSC'].SetLineStyle(7)
        nlo['STD'].SetLineStyle(4)
        [gr.SetLineWidth(3) for gr in nlo.values()]
        [gr.SetFillColor(ROOT.kOrange) for gr in nlo.values()]
    
    leg = ROOT.TLegend(0.69, 0.15, 0.97, 0.36)
    leg.SetBorderSize(0)
    leg.SetNColumns(1)
    leg.AddEntry(nlo_m['STD'], 'GRSV-STD', 'l')
    leg.AddEntry(nlo_m['DSSV'], 'DSSV', 'l')
    leg.AddEntry(nlo_m['GSC'], 'GS-C', 'l')
    
    latex2 = ROOT.TLatex()
    
    c.cd(1)
    hm.SetMarkerStyle(20)
    syst_m.Draw('a2p')
    # prelim.DrawText(0.42,-0.041,"STAR Preliminary")
    latex.DrawLatex(0.05 + zbins[0], 0.044, '#pi^{-}')
    latex2.DrawLatex(0.45, 0.06, '10 < jet p_{T} < 30')
    latex2.DrawLatex(0.5, 0.045, '#pi p_{T} > 2')
    [ g.Draw('3') for g in nlo_m.values() ]
    [ g.Draw('lx') for g in nlo_m.values() ]
    line.Draw()
    syst_m.Draw('2p')
    hm.Draw('e1 same')
    
    c.cd(2)
    hp.SetMarkerStyle(21)
    syst_p.Draw('a2p')
    # prelim.DrawText(0.42,-0.041,"STAR Preliminary")
    latex.DrawLatex(0.05 + zbins[0], 0.044, '#pi^{+}')
    [ g.Draw('3') for g in nlo_p.values() ]
    [ g.Draw('lx') for g in nlo_p.values() ]
    line.Draw()
    leg.Draw()
    syst_p.Draw('2p')
    hp.Draw('e1 same')
    
    ## auto-generate LaTeX tables
    for i in range(4):
        print '%.2f - %.2f & %.1f & $\pm$ %.1f & -%.1f +%.1f &  %.1f & $\pm$ %.1f & -%.1f +%.1f\\\\'\
            % (zbins[i], zbins[i+1], 
            hm.GetBinContent(i+1)*1000, hm.GetBinError(i+1)*1000,
            syst_m.GetErrorYlow(i)*1000, syst_m.GetErrorYhigh(i)*1000,
            hp.GetBinContent(i+1)*1000, hp.GetBinError(i+1)*1000,
            syst_p.GetErrorYlow(i)*1000, syst_p.GetErrorYhigh(i)*1000)
    
    graphics.maybe_save(c)
    return (hm.Clone(), hp.Clone())

def systematic_uncertainties_run5():
    ## systematic uncertainties x 10^-3 
    
    ## trig+reco after rescaling, dropping MAX,MIN,P*
    # trig_m_l = [5.90,  8.30,  9.30,  7.20,  4.80]
    # trig_m_h = [2.70,  7.20,  3.40,  3.60,  5.70]
    # trig_p_l = [6.10, 12.80, 17.60, 20.90, 15.20]
    # trig_p_h = [6.10,  6.60, 10.10, 18.60,  6.20]
    
    ## trig+reco after rescaling, M030-P030
    trig_m_h = [1.50,  2.30,  2.40,  4.10,  5.20]
    trig_m_l = [7.60,  10.1,  10.5,  9.70,  7.0]
    trig_p_h = [1.50,  2.1, 2.5, 4.1, 5.1]
    trig_p_l = [2.10,  14.90, 20.20, 22.80,  16.20]
    
    asigma_m = [2.91, 2.30, 3.74, 5.98, 12.07]
    asigma_p = [1.34, 2.56, 5.96, 5.72, 15.11]
    
    rellumi = [0.98, 0.98, 0.98, 0.98, 0.98]
    
    sum_quadrature = lambda *args: math.sqrt( sum([a**2 for a in args]) )/1E3
    
    return {
        'minus': {
            'low':  map(sum_quadrature, trig_m_l, asigma_m, rellumi),
            'high': map(sum_quadrature, trig_m_h, asigma_m, rellumi)
        },
        'plus': {
            'low':  map(sum_quadrature, trig_p_l, asigma_p, rellumi),
            'high': map(sum_quadrature, trig_p_h, asigma_p, rellumi)
        }
    }

def systematic_uncertainties_run6():
    ## systematic uncertainties x 10^-3 
    
    # 0.20 - 0.30  & -0.0005 +0.0052 & -0.0009 +0.0134 \\
    # 0.30 - 0.45  & -0.0010 +0.0113 & -0.0023 +0.0165 \\
    # 0.45 - 0.65  & -0.0017 +0.0017 & -0.0040 +0.0032 \\
    # 0.65 - 1.00  & -0.0035 +0.0108 & -0.0076 +0.0076 \\
    
    ## trig+reco after rescaling, M030 - P030
    trig_m_h = [0.5, 1.0, 1.7, 3.5]
    trig_m_l = [5.2, 11.3, 1.7, 10.8]
    trig_p_h = [0.9, 2.3, 4.0, 7.6]
    trig_p_l = [13.4, 16.5, 3.2, 7.6]
    
    asigma_m = [0.1, 0.16, 0.4, 0.44]
    asigma_p = [0.09, 0.15, 0.23, 0.37]
    
    rellumi = [0.98, 0.98, 0.98, 0.98]
    
    jetpt_m = [3, 5, 16, 10]
    jetpt_p = [4, 7, 16, 16]
    
    sum_quadrature = lambda *args: math.sqrt( sum([a**2 for a in args]) )/1E3
    
    d = {
        'minus': {
            'low':  map(sum_quadrature, trig_m_l, asigma_m, rellumi, jetpt_m),
            'high': map(sum_quadrature, trig_m_h, asigma_m, rellumi, jetpt_m)
        },
        'plus': {
            'low':  map(sum_quadrature, trig_p_l, asigma_p, rellumi, jetpt_p),
            'high': map(sum_quadrature, trig_p_h, asigma_p, rellumi, jetpt_p)
        }
    }
    
    for i in range(4):
        print '-%.1f +%.1f & -%.1f +%.1f' % (
            1000*d['minus']['low'][i], 1000*d['minus']['high'][i],
            1000*d['plus']['low'][i], 1000*d['plus']['high'][i])
            
    return d

def systematic_uncertainties():
    """tabulates sources of uncertainty and sums them in quadrature"""
    
    result_m = [
         0.066,     # [0.07-0.12]   0.066 ± 0.019
         0.019,     # [0.12-0.20]   0.019 ± 0.009
         0.002,     # [0.20-0.30]   0.002 ± 0.009
        -0.006,     # [0.30-0.45]  -0.006 ± 0.014
         0.007,     # [0.45-0.65]   0.007 ± 0.023
         0.012      # [0.65-1.00]   0.012 ± 0.040
    ]
    result_p = [
         0.026,     # [0.07-0.12]   0.026 ± 0.019
         0.021,     # [0.12-0.20]   0.021 ± 0.008
         0.002,     # [0.20-0.30]   0.002 ± 0.009
        -0.014,     # [0.30-0.45]  -0.014 ± 0.013
         0.024,     # [0.45-0.65]   0.024 ± 0.022
         0.046      # [0.65-1.00]   0.046 ± 0.037
    ]
    
    pid_contamination = 0.10
    pid_asym_m = [
        ( 0.051 , 0.038),      # [0.07-0.12]   0.051 ± 0.038
        (-0.017 , 0.016),      # [0.12-0.20]  -0.017 ± 0.016
        (-0.032 , 0.016),      # [0.20-0.30]  -0.032 ± 0.016
        (-0.006 , 0.023),      # [0.30-0.45]  -0.006 ± 0.023
        (-0.031 , 0.042),      # [0.45-0.65]  -0.031 ± 0.042
        ( 0.089 , 0.085)       # [0.65-1.00]   0.089 ± 0.085
    ]
    pid_asym_p = [
        ( 0.005 , 0.036),      # [0.07-0.12]   0.005 ± 0.036
        ( 0.006 , 0.015),      # [0.12-0.20]   0.006 ± 0.015
        (-0.006 , 0.015),      # [0.20-0.30]  -0.006 ± 0.015
        ( 0.018 , 0.020),      # [0.30-0.45]   0.018 ± 0.020
        (-0.038 , 0.032),      # [0.45-0.65]  -0.038 ± 0.032
        ( 0.142 , 0.059)       # [0.65-1.00]   0.142 ± 0.059
    ]
    
    for i in range(len(pid_asym_m)):
        val, err = pid_asym_m[i]
        pid_asym_m[i] = max( val-result_m[i], err)
    for i in range(len(pid_asym_p)):
        val, err = pid_asym_p[i]
        pid_asym_p[i] = max( val-result_p[i], err)
    
    beam_vector = 0.0102
    asigma_m = [
        0.035,      # [0.07-0.12]   0.005 ± 0.035
        0.015,      # [0.12-0.20]  -0.012 ± 0.015
        0.016,      # [0.20-0.30]  -0.014 ± 0.016
        0.027,      # [0.30-0.45]  -0.027 ± 0.023
        0.066,      # [0.45-0.65]  -0.066 ± 0.040
        0.073       # [0.65-1.00]  -0.072 ± 0.073
    ]
    asigma_p = [
        0.034,      # [0.07-0.12]  -0.001 ± 0.034
        0.014,      # [0.12-0.20]  -0.007 ± 0.014
        0.015,      # [0.20-0.30]   0.007 ± 0.015
        0.025,      # [0.30-0.45]  -0.025 ± 0.022
        0.039,      # [0.45-0.65]  -0.039 ± 0.037
        0.061       # [0.65-1.00]   0.033 ± 0.061
    ]
    
    mcasym_m = [
        0.0066,     # [0.07-0.12]   0.0012 ± 0.0066
        0.0057,     # [0.12-0.20]   0.0057 ± 0.0025
        0.0089,     # [0.20-0.30]   0.0089 ± 0.0020
        0.0077,     # [0.30-0.45]   0.0077 ± 0.0026
        0.0042,     # [0.45-0.65]   0.0038 ± 0.0042
        0.0070      # [0.65-1.00]   0.0053 ± 0.0070
    ]
    mcasym_p = [
        0.0047,     # [0.07-0.12]  -0.0014 ± 0.0047
        0.0077,     # [0.12-0.20]   0.0077 ± 0.0024
        0.0147,     # [0.20-0.30]   0.0147 ± 0.0023
        0.0105,     # [0.30-0.45]   0.0105 ± 0.0024
        0.0057,     # [0.45-0.65]   0.0057 ± 0.0044
        0.0112      # [0.65-1.00]   0.0112 ± 0.0081
    ]
    
    pt_shift_m = [ 0, 0,
        0.003, # [0.20-0.30]   0.006 low, 0.001 high, 0.003 avg
        0.005, # [0.30-0.45]   0.007 low, 0.003 high, 0.005 avg
        0.016, # [0.45-0.65]   0.020 low, 0.012 high, 0.016 avg
        0.010  # [0.65-1.00]   0.011 low, 0.008 high, 0.010 avg
    ]
    
    pt_shift_p = [ 0, 0,
        0.004, # [0.20-0.30]   0.005 low, 0.003 high, 0.004 avg
        0.007, # [0.30-0.45]   0.008 low, 0.006 high, 0.007 avg
        0.016, # [0.45-0.65]   0.023 low, 0.008 high, 0.016 avg
        0.016  # [0.65-1.00]   0.012 low, 0.020 high, 0.016 avg
    ]
    
    relative_luminosity = 9.4e-4
    
    minus = [0.0 for bin in zbins[:-1]]
    plus = [0.0 for bin in zbins[:-1]]
    
    start = len(zbins) == 5 and 2 or 0
    for i in range(start, start+len(zbins)-1):
        minus[i-start] = math.sqrt(
            pow(relative_luminosity, 2) + 
            # pow(pid_contamination*pid_asym_m[i], 2) +
            pow(beam_vector*asigma_m[i], 2) +
            pow(mcasym_m[i], 2) + 
            pow(pt_shift_m[i], 2)
        )
        plus[i-start] = math.sqrt(
            pow(relative_luminosity, 2) +
            # pow(pid_contamination*pid_asym_p[i], 2) + 
            pow(beam_vector*asigma_p[i], 2) + 
            pow(mcasym_p[i], 2) + 
            pow(pt_shift_p[i], 2)
        )
    
    return {'minus':minus, 'plus':plus}

def ch1_2005_predictions():
    ROOT.gStyle.SetErrorX(0)
    ROOT.gStyle.SetOptDate(0)
    
    h = ROOT.TH1D('h', '', 1, 2.0, 12.74)
    
    ## NLO curves
    from analysis.asym import theoryCurves
    nlo_m = {
        'STD': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_std, 
            analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
        'ZERO': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_zero, 
            analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
        'MIN': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_min, 
            analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),    
        'GSC': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_gsc, 
            analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph()
    }
    nlo_p = {
        'STD': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_std, 
            analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
        'ZERO': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_zero, 
            analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
        'MIN': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_min, 
            analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
        'GSC': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_gsc, 
            analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph()
    }
    
    nlo_m['DSSV'] = ROOT.TGraphErrors(15,
        array('d', [float(i) for i in range(1,16)]),
        array('d', [
            -2.685E-05, -4.188E-05, -8.716E-05, -1.562E-04, -1.841E-04,
            -1.420E-04, -5.366E-05,  6.018E-05,  1.274E-04,  1.235E-04,
             7.536E-06, -2.645E-04, -6.844E-04, -1.261E-03, -1.999E-03
        ]),
        array('d', [0 for i in range(15)]),
        array('d', [
            3.011E-05, 3.331E-05, 1.549E-04, 3.664E-04, 5.990E-04,
            8.074E-04, 9.963E-04, 1.168E-03, 1.345E-03, 1.513E-03,
            1.684E-03, 1.851E-03, 1.990E-03, 2.117E-03, 2.223E-03
        ]))
    
    nlo_p['DSSV'] = ROOT.TGraphErrors(15,
        array('d', [float(i) for i in range(1,16)]),
        array('d', [
            -3.800E-05, -3.658E-05,  7.296E-07,  1.634E-04,  6.061E-04,
             1.404E-03,  2.588E-03,  4.144E-03,  5.984E-03,  8.075E-03,
             1.031E-02,  1.271E-02,  1.528E-02,  1.779E-02,  2.039E-02
        ]),
        array('d', [0 for i in range(15)]),
        array('d', [
            1.821E-05, 8.788E-05, 3.054E-04, 6.487E-04, 1.038E-03, 
            1.406E-03, 1.758E-03, 2.090E-03, 2.417E-03, 2.745E-03, 
            3.080E-03, 3.430E-03, 3.750E-03, 4.090E-03, 4.390E-03
        ]))
    
    for nlo in (nlo_m, nlo_p):
        nlo['MIN'].SetLineStyle(3)
        nlo['STD'].SetLineStyle(4)
        nlo['ZERO'].SetLineStyle(8)
        nlo['GSC'].SetLineStyle(7)
        
        nlo['ZERO'].SetLineColor(ROOT.kBlue)
        nlo['MIN'].SetLineColor(ROOT.kGreen)
        nlo['GSC'].SetLineColor(ROOT.kMagenta)
        nlo['DSSV'].SetFillColor(ROOT.kYellow)
        [gr.SetLineWidth(3) for gr in nlo.values()]
        nlo['DSSV'].SetLineWidth(3)
    
    ## graphics setup
    h.GetYaxis().SetRangeUser(-0.042, 0.042)
    h.GetYaxis().SetNdivisions(507)
    h.GetYaxis().SetTitle('A_{LL}    ')
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.7)
    h.GetYaxis().SetLabelOffset(0.015)
    h.GetXaxis().SetTitle('p_{T}  [GeV/c]  ')
    h.GetXaxis().SetTitleSize(0.045)
    h.GetXaxis().SetTitleOffset(1.16)
    
    line = ROOT.TLine(2.00, 0.0, 12.84, 0.0)
    line.SetLineStyle(2)
    
    leg = ROOT.TLegend(0.27, 0.2, 0.95, 0.38)
    leg.SetNColumns(2)
    leg.SetColumnSeparation(-0.2)
    leg.SetBorderSize(0)
    leg.AddEntry(nlo_m['STD'], 'GRSV STD', 'l')
    leg.AddEntry(nlo_m['DSSV'], 'DSSV', 'l')
    leg.AddEntry(nlo_m['MIN'], 'GRSV MIN', 'l')
    leg.AddEntry(nlo_m['GSC'], 'GS Set C', 'l')
    leg.AddEntry(nlo_m['ZERO'], 'GRSV ZERO', 'l')
    
    latex = ROOT.TLatex()
    latex.SetTextSize(0.23)
    
    poltext = ROOT.TLatex()
    poltext.SetTextSize(0.0415)
    
    c = graphics.canvas2()
    
    for i in (1, 2):
        pad = c.cd(i)
        pad.SetTopMargin(0.03)
        pad.SetBottomMargin(0.12)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.02)
        
    pad1 = c.cd(1)
    h.DrawCopy()
    nlo_m['DSSV'].Draw('lx')
    [ g.Draw('lx') for g in nlo_m.values() ]
    line.Draw()
    latex.DrawLatex(2.8, -0.033, '#pi^{-}')
    
    pad2 = c.cd(2)
    h.Draw()
    nlo_p['DSSV'].Draw('lx')
    [ g.Draw('lx') for g in nlo_p.values() ]
    line.Draw()
    latex.DrawLatex(2.8, -0.033, '#pi^{+}')
    
    graphics.maybe_save(c)


def ch1_2006_predictions():
    ## fun with graphics
    ROOT.gStyle.SetErrorX(0.0)
    ROOT.gStyle.SetOptDate(0)
    
    c = graphics.canvas2()
    c.SetLogy(0)
    
    for i in (1, 2):
        pad = c.cd(i)
        pad.SetTopMargin(0.03)
        pad.SetBottomMargin(0.12)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.02)
        
    
    line = ROOT.TLine(zbins[0], 0.0, zbins[-1], 0.0)
    line.SetLineStyle(2)
    line.SetLineColor(10)
    
    latex = ROOT.TLatex()
    latex.SetTextSize(0.23)
    # latex.SetTextAlign(21)
    
    h = ROOT.TH1D('h', '', 1, 0.2, 1.0)
    h.GetYaxis().SetNdivisions(507)
    h.GetYaxis().SetTitle('A_{LL}    ')
    h.GetXaxis().SetTitle('z = p_{T}(#pi)/p_{T}(jet)')
    h.GetYaxis().SetRangeUser(-0.042, 0.042)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(0.7)
    h.GetYaxis().SetLabelOffset(0.015)
    h.GetXaxis().SetTitleSize(0.045)
    h.GetXaxis().SetTitleOffset(1.16)
    
    c.cd(1)
    line.Draw()
    latex.DrawLatex(0.07 + zbins[0], -0.085, '#pi^{-}')
    
    c.cd(2)
    line.Draw()
    latex.DrawLatex(0.07 + zbins[0], -0.085, '#pi^{+}')
    
    title = {
        'STD': 'GRSV-STD',
        'GSC': 'GS Set C',
        'DSSV': 'DSSV'
    }
    
    def make_graph(num, denom):
        g = ROOT.TGraphErrors(len(num))
        for i,n in enumerate(num):
            d = denom[i]
            val = n.y / d.y
            err = val*math.sqrt((n.sys/n.y)**2 + (d.sys/d.y)**2)
            g.SetPoint(i, n.x, val)
            g.SetPointError(i, 0, 0) # no error bars here
        return g
    
    nlo_m = {
        'STD': make_graph(analysis.deflorian.minus.std, 
            analysis.deflorian.minus.mrst),
        'DSSV': make_graph(analysis.deflorian.minus.dssv, 
            analysis.deflorian.minus.mrst),
        'GSC': make_graph(analysis.deflorian.minus.gsc, 
            analysis.deflorian.minus.mrst),
    }
    nlo_p = {
        'STD': make_graph(analysis.deflorian.plus.std, 
            analysis.deflorian.plus.mrst),
        'DSSV': make_graph(analysis.deflorian.plus.dssv, 
            analysis.deflorian.plus.mrst),
        'GSC': make_graph(analysis.deflorian.plus.gsc, 
            analysis.deflorian.plus.mrst),
    }
    
    for nlo in (nlo_m, nlo_p):
        # nlo['MIN'].SetLineStyle(3)
        nlo['STD'].SetLineStyle(4)
        # nlo['ZERO'].SetLineStyle(8)
        nlo['GSC'].SetLineStyle(7)
        # nlo['DSSV'].SetLineStyle(2)
        nlo['GSC'].SetLineColor(ROOT.kMagenta)
        # nlo['DSSV'].SetLineColor(ROOT.kGreen+2)
        [gr.SetLineWidth(3) for gr in nlo.values()]
        [gr.SetFillColor(ROOT.kOrange) for gr in nlo.values()]
    
    legmin = ROOT.TGraph()
    legmin.SetLineStyle(3)
    legmin.SetLineWidth(3)
    legmin.SetLineColor(ROOT.kGreen)
    
    legzero = ROOT.TGraph()
    legzero.SetLineStyle(8)
    legzero.SetLineWidth(3)
    legzero.SetLineColor(ROOT.kBlue)
    
    leg = ROOT.TLegend(0.39, 0.18, 1.05, 0.38)
    leg.SetFillStyle(0)
    leg.SetNColumns(2)
    leg.SetColumnSeparation(-0.2)
    leg.SetBorderSize(0)
    leg.AddEntry(nlo_m['STD'], 'GRSV STD', 'l')
    leg.AddEntry(nlo_m['DSSV'], 'DSSV', 'l')
    leg.AddEntry(legmin, 'GRSV MIN', 'l')
    leg.AddEntry(nlo_m['GSC'], 'GS Set C', 'l')
    leg.AddEntry(legzero, 'GRSV ZERO', 'l')
    # leg.AddEntry(nlo_m['MIN'], 'GRSV MIN', 'l')
    # leg.AddEntry(nlo_m['ZERO'], 'GRSV ZERO', 'l')
    # leg = ROOT.TLegend(0.62, 0.15, 0.9, 0.36)
    # leg.SetBorderSize(0)
    # leg.SetNColumns(1)
    # leg.AddEntry(nlo_m['STD'], 'GRSV-STD', 'l')
    # leg.AddEntry(nlo_m['DSSV'], 'DSSV', 'l')
    # leg.AddEntry(nlo_m['GSC'], 'GS-C', 'l')
    
    latex2 = ROOT.TLatex()
    
    c.cd(1)
    h.DrawCopy()
    # hm.SetMarkerStyle(20)
    # syst_m.Draw('a2p')
    # prelim.DrawText(0.42,-0.041,"STAR Preliminary")
    latex.DrawLatex(0.05 + zbins[0], -0.033, '#pi^{-}')
    latex2.DrawLatex(0.27, 0.032, '10 < jet p_{T} < 30')
    latex2.DrawLatex(0.275,  0.025, '#pi p_{T} > 2')
    [ g.Draw('3') for g in nlo_m.values() ]
    [ g.Draw('lx') for g in nlo_m.values() ]
    line.Draw()
    # syst_m.Draw('2p')
    # hm.Draw('e1 same')
    
    c.cd(2)
    h.DrawCopy()
    latex2.DrawLatex(0.27, 0.032, '10 < jet p_{T} < 30')
    latex2.DrawLatex(0.275,  0.025, '#pi p_{T} > 2')
    # hp.SetMarkerStyle(21)
    # syst_p.Draw('a2p')
    # prelim.DrawText(0.42,-0.041,"STAR Preliminary")
    latex.DrawLatex(0.05 + zbins[0], -0.033, '#pi^{+}')
    [ g.Draw('3') for g in nlo_p.values() ]
    [ g.Draw('lx') for g in nlo_p.values() ]
    line.Draw()
    leg.Draw()
    h.DrawCopy('same')
    # syst_p.Draw('2p')
    # hp.Draw('e1 same')
    graphics.maybe_save(c)

def beam_polarizations():
    ROOT.gStyle.SetOptDate(0)
    pol = Polarizations.Final
    fills05 = set([analysis.util.fill(r) for r in runlist05])
    fills06 = set([analysis.util.fill(r) for r in runlist])
    hblue05 = ROOT.TH1D('hblue05', '', len(fills05), 1, len(fills05)+1)
    hblue06 = ROOT.TH1D('hblue06', '', len(fills06), 1, len(fills06)+1)
    hyellow05 = ROOT.TH1D('hyellow05', '', len(fills05), 1, len(fills05)+1)
    hyellow06 = ROOT.TH1D('hyellow06', '', len(fills06), 1, len(fills06)+1)
    
    [h.SetMarkerColor(ROOT.kBlue) for h in (hblue05,hblue06)]
    [h.SetMarkerColor(ROOT.kYellow) for h in (hyellow05,hyellow06)]
    [h.SetMarkerStyle(25) for h in (hblue05,hblue06)]
    [h.SetMarkerStyle(20) for h in (hyellow05,hyellow06)]
    
    for i,fill in enumerate(fills05):
        hblue05.SetBinContent(i+1, pol[fill].pb)
        hblue05.SetBinError(i+1, pol[fill].eb)
        if i%10 == 2:
            hblue05.GetXaxis().SetBinLabel(i+1, str(fill))

        hyellow05.SetBinContent(i+1, pol[fill].py)
        hyellow05.SetBinError(i+1, pol[fill].ey)

    for i,fill in enumerate(fills06):
        hblue06.SetBinContent(i+1, pol[fill].pb)
        hblue06.SetBinError(i+1, pol[fill].eb)
        if i%5 == 1:
            hblue06.GetXaxis().SetBinLabel(i+1, str(fill))

        hyellow06.SetBinContent(i+1, pol[fill].py)
        hyellow06.SetBinError(i+1, pol[fill].ey)

    for h in hblue05,hblue06:
        h.GetXaxis().LabelsOption('d')
        h.GetYaxis().SetRangeUser(0.3,0.7)
        h.GetYaxis().SetNdivisions(5,2,0)
        h.GetYaxis().SetLabelSize(0.045)
        h.GetXaxis().SetTitleOffset(1.65)
        h.GetXaxis().SetTitleSize(0.05)
        h.GetXaxis().SetLabelSize(0.06)
        h.GetXaxis().SetLabelOffset(0.01)
    hblue06.SetXTitle('RHIC fill number')
    
    c05 = graphics.canvas1(None, 750, 350)
    hblue05.Draw('e1 x0')
    hyellow05.Draw('e1 x0 same')

    c06 = graphics.canvas1(None, 750, 350)
    hblue06.Draw('e1 x0')
    hyellow06.Draw('e1 x0 same')

    leg = ROOT.TLegend(0.5, 0.25, 0.95, 0.4)
    leg.AddEntry(hblue06, 'Blue Beam', 'p')
    leg.AddEntry(hyellow06, 'Yellow Beam', 'p')
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.Draw()
    
    for pad in c05,c06:
        pad.SetRightMargin(0.005)
        pad.SetTopMargin(0.06)
        pad.SetLeftMargin(0.035)
        pad.SetBottomMargin(0.17)
    
    # calculate final beam polarization
    poly05 = sum([pol[fill].py / (pol[fill].ey**2) for fill in fills05]) / \
        sum([1.0/(pol[fill].ey**2) for fill in fills05])
    polb05 = sum([pol[fill].pb / (pol[fill].eb**2) for fill in fills05]) / \
        sum([1.0/(pol[fill].eb**2) for fill in fills05])
    poly06 = sum([pol[fill].py / (pol[fill].ey**2) for fill in fills06]) / \
        sum([1.0/(pol[fill].ey**2) for fill in fills06])
    polb06 = sum([pol[fill].pb / (pol[fill].eb**2) for fill in fills06]) / \
        sum([1.0/(pol[fill].eb**2) for fill in fills06])
    print 'Yellow 05', poly05
    print 'Blue   05', polb05
    print 'Yellow 06', poly06
    print 'Blue   06', polb06

    graphics.maybe_save(c05)
    graphics.maybe_save(c06)

def r3():
    db = sqlite.connect('/Users/kocolosk/data/analysis.db')
    dbc = db.cursor()
    
    h05 = ROOT.TH1D('h05', '2005', 70, 0.65, 1.35)
    h06 = ROOT.TH1D('h06', '2006', 70, 0.65, 1.35)
    
    for run in runlist05:
        uu,ud,du,dd = analysis.util.scaler_counts(dbc, run)
        h05.Fill(float(uu+dd)/(ud+du))
    
    for run in runlist:
        uu,ud,du,dd = analysis.util.scaler_counts(dbc, run)
        h06.Fill(float(uu+dd)/(ud+du))
    
    for h in h05, h06:
        h.SetXTitle('R')
        h.SetFillColor(14)
    c = graphics.canvas2()
    c.cd(1)
    ROOT.gPad.SetRightMargin(0.005)
    ROOT.gPad.SetLeftMargin(0.06)
    ROOT.gPad.SetBottomMargin(0.08)
    h05.Draw()
    c.cd(2)
    ROOT.gPad.SetRightMargin(0.005)
    ROOT.gPad.SetLeftMargin(0.06)
    ROOT.gPad.SetBottomMargin(0.08)
    h06.Draw()
    graphics.maybe_save(c)

def vertex_stats():
    h = {}
    for trig in ('96011', '96221', '96233', '117001', 'jetpatch'):
        hdir = trig.startswith('9') and histDir05 or histDir
        rd = trig.startswith('9') and runlist05 or runlist
        h[trig] = hadd_interactive(hdir, rd, trig, 'anyspin', None, 'nVertices')
    for histo in h.values():
        graphics.canvas1()
        h.Draw()
    raw_input('hold')

def vertex_stats(trig):
    hdir = trig.startswith('9') and histDir05 or histDir
    rd = trig.startswith('9') and runlist05 or runlist
    h = hadd_interactive(hdir, rd, trig, 'anyspin', None, 'nVertices')
    total_events = h.GetEntries()
    with_vertex = h.Integral(2, h.GetNbinsX())
    print int(total_events), '&', int(with_vertex), '&', \
        100*with_vertex/total_events

def asigma_2005():
    asigma('2005', analysis.runlists.transverse_run5)

def asigma_2006():
    asigma('2006', analysis.runlists.transverse_run6)

def asigma(year, runlist):
    if year == '2005':
        asym_p = AsymmetryGenerator('asym_p', bins=ptbins, key='pt')
        asym_m = AsymmetryGenerator('asym_m', bins=ptbins, key='pt')
    else:
        asym_p = AsymmetryGenerator('asym_p', bins=zbins, key='z_away2')
        asym_m = AsymmetryGenerator('asym_m', bins=zbins, key='z_away2')
        
    db = sqlite.connect('/Users/kocolosk/data/analysis.db')
    dbc = db.cursor()
    
    polarizations = Polarizations.Final
    
    ## generate the asymmetries
    if year == '2005':
        transHistDir = '/Users/kocolosk/work/2010-04-11-transverse-pt-histos'
    else:
        transHistDir = '/Users/kocolosk/data/run6/spin2008/hist-transverse'
    allFiles = glob(transHistDir + '/chargedPions_*.hist.root')
    for fname in allFiles[:]:
        run = getRun(fname)
        if run in runlist:
            print fname, run
            if year == '2005':
                mgr = HistogramManager(ROOT.TFile(fname), ['pt'])
            else:
                mgr = HistogramManager(ROOT.TFile(fname), ['z_away2'])
            uu,ud,du,dd = analysis.util.scaler_counts(dbc, run)
            pol = polarizations[analysis.util.fill(run)]
            asym_p.FillFromHistogramManager(mgr, 'jetpatch', 1, uu,ud,du,dd, \
                pol.py,pol.pb)
            asym_m.FillFromHistogramManager(mgr, 'jetpatch', -1, uu,ud,du,dd, \
                pol.py,pol.pb)
    
    c = graphics.canvas2()
    c.SetLogy(0)
    ROOT.gStyle.SetErrorX(0.0)
    
    if year == '2005':
        line = ROOT.TLine(ptbins[0], 0.0, ptbins[-1], 0.0)
    else:
        line = ROOT.TLine(zbins[0], 0.0, zbins[-1], 0.0)        
    line.SetLineStyle(2)
    
    c.cd(1)
    ROOT.gPad.SetRightMargin(0.005)
    hm = asym_m.GetAsymmetry('ll')
    hm.SetMarkerStyle(20)
    hm.SetTitle(year + ' A_{#Sigma} #pi^{-}')
    # hm.Fit('pol0', 'q')
    hm.Draw('e1')
    line.Draw('same')
    
    c.cd(2)
    ROOT.gPad.SetRightMargin(0.005)
    hp = asym_p.GetAsymmetry('ll')
    hp.SetTitle(year + ' A_{#Sigma} #pi^{+}')
    hp.SetMarkerStyle(21)
    # hp.Fit('pol0', 'q')
    hp.Draw('e1')
    line.Draw('same')
    
    for h in (hm,hp):
        if year == '2005':
            h.GetXaxis().SetTitle('p_{T}')
            h.GetYaxis().SetRangeUser(-1.0, 1.0)
        else:
            h.GetXaxis().SetTitle('p_{T}(#pi)/p_{T}(jet)')            
            h.GetYaxis().SetRangeUser(-0.1, 0.1)
    
    raw_input('wait here:')
    c.Print('asigma_%s.png' % (year,))
    
    print '\nπ- A_{σ}'
    for bin in range(1, hm.GetNbinsX()+1):
        print '[%.2f-%.2f]  % .3f ± %.3f' % (hm.GetBinLowEdge(bin), 
            hm.GetBinLowEdge(bin+1), hm.GetBinContent(bin), hm.GetBinError(bin))
    
    print '\nπ+ A_{σ}'
    for bin in range(1, hp.GetNbinsX()+1):
        print '[%.2f-%.2f]  % .3f ± %.3f' % (hp.GetBinLowEdge(bin), 
            hp.GetBinLowEdge(bin+1), hp.GetBinContent(bin), hp.GetBinError(bin))

def pt_shift_uncertainty():
    ## this is my guess based on inverting Dave's fit
    measured_pt = [10.31, 12.91, 15.58, 19.06, 23.20, 28.45]
    
    ## http://cyclotron.tamu.edu/star/2005n06Jets/PRDweb/
    ## "numbers are from the preliminary analysis"
    total_uncertainty = [ 0.72, 0.90, 0.99, 1.09, 1.27, 1.52 ]
    
    c = graphics.canvas1()
    
    low = ROOT.TGraph(len(measured_pt))
    for i,(pt,error) in enumerate(zip(measured_pt, total_uncertainty)):
        low.SetPoint(i, pt, shifted(pt)-error)
    lowfit = ROOT.TF1('lowfit', 'pol2', 10, 30)
    lowfit.SetLineStyle(2)
    low.Fit(lowfit)
    low.GetXaxis().SetTitle('measured jet p_{T}')
    low.GetXaxis().SetRangeUser(10,30)
    low.GetYaxis().SetTitle('corrected jet p_{T}')
    low.GetYaxis().SetRangeUser(9,26)
    low.SetTitle('Uncertainty on p_{T} shift')
    low.SetLineStyle(2)
    low.Draw('ap2')
    
    high = ROOT.TGraph(len(measured_pt))
    for i,(pt,error) in enumerate(zip(measured_pt, total_uncertainty)):
        high.SetPoint(i, pt, shifted(pt)+error)
    highfit = ROOT.TF1('lowfit', 'pol2', measured_pt[0], measured_pt[-1])
    highfit.SetLineStyle(2)
    high.Fit(highfit)
    high.SetLineStyle(2)
    high.Draw('same')
    
    mid = ROOT.TF1('mid', lambda x: shifted(x[0]), 10.31, 28.45)
    mid.Draw('same')
    
    raw_input('wait here:')

def pt_shift_uncertainty_asymmetry():
    """
    calculate A_{LL} assuming different jet pt shifts
    """
    asym_pl = AsymmetryGenerator('asym_pl', bins=zbins, key='z_away2_low')
    asym_pm = AsymmetryGenerator('asym_pm', bins=zbins, key='z_away2')
    asym_ph = AsymmetryGenerator('asym_ph', bins=zbins, key='z_away2_high')
    asym_ml = AsymmetryGenerator('asym_ml', bins=zbins, key='z_away2_low')
    asym_mm = AsymmetryGenerator('asym_mm', bins=zbins, key='z_away2')
    asym_mh = AsymmetryGenerator('asym_mh', bins=zbins, key='z_away2_high')
    
    scalars = ScalarCounts(os.environ['STAR'] + 
        '/StRoot/StSpinPool/StTamuRelLum/inputs/run6.txt')
    
    polarizations = Polarizations.Final
    
    ## generate the asymmetries
    allFiles = glob(histDir + '/chargedPions_*.hist.root')
    for fname in allFiles[:]:
        run = getRun(fname)
        if run in runlist:
            print fname, run
            mgr = HistogramManager(ROOT.TFile(fname), ['z_away2', 
                'z_away2_low', 'z_away2_high'])
            
            try:
                bin6 = scalars[str(run) + '-5-6']
                bin7 = scalars[str(run) + '-5-7']
                bin8 = scalars[str(run) + '-5-8']
                bin9 = scalars[str(run) + '-5-9']
            except KeyError:
                bin6 = scalars[str(run) + '-6-6']
                bin7 = scalars[str(run) + '-6-7']
                bin8 = scalars[str(run) + '-6-8']
                bin9 = scalars[str(run) + '-6-9']
            uu = bin6.uu + bin7.uu + bin8.uu + bin9.uu
            ud = bin6.ud + bin7.ud + bin8.ud + bin9.ud
            du = bin6.du + bin7.du + bin8.du + bin9.du
            dd = bin6.dd + bin7.dd + bin8.dd + bin9.dd
            
            pol = polarizations[bin7.fill]
            
            asym_ml.FillFromHistogramManager(mgr, 'jetpatch', -1, uu,ud,du,dd, \
                pol.py,pol.pb)
            asym_mm.FillFromHistogramManager(mgr, 'jetpatch', -1, uu,ud,du,dd, \
                pol.py,pol.pb)
            asym_mh.FillFromHistogramManager(mgr, 'jetpatch', -1, uu,ud,du,dd, \
                pol.py,pol.pb)
            asym_pl.FillFromHistogramManager(mgr, 'jetpatch', 1, uu,ud,du,dd, \
                pol.py,pol.pb)
            asym_pm.FillFromHistogramManager(mgr, 'jetpatch', 1, uu,ud,du,dd, \
                pol.py,pol.pb)
            asym_ph.FillFromHistogramManager(mgr, 'jetpatch', 1, uu,ud,du,dd, \
                pol.py,pol.pb)
    
    
    c = graphics.canvas2()
    c.SetLogy(0)
    ROOT.gStyle.SetErrorX(0.0)
    
    c.cd(1)
    hml = asym_ml.GetAsymmetry('ll')
    hmm = asym_mm.GetAsymmetry('ll')
    hmh = asym_mh.GetAsymmetry('ll')
    hmm.SetMarkerStyle(20)
    hmm.SetTitle('p_{T} shift uncertainty on A_{LL} #pi^{-}')
    # hml.Fit('pol0', 'q')
    hmm.Draw('e1')
    gml = ROOT.TGraphErrors(hml)
    gmh = ROOT.TGraphErrors(hmh)
    for i in range(hmm.GetNbinsX()):
        gml.SetPoint(i, hml.GetBinCenter(i+1)-0.02, hml.GetBinContent(i+1))
        gmh.SetPoint(i, hmh.GetBinCenter(i+1)+0.02, hmh.GetBinContent(i+1))
    gml.Draw('p')
    gmh.Draw('p')
    
    c.cd(2)
    hpl = asym_pl.GetAsymmetry('ll')
    hpm = asym_pm.GetAsymmetry('ll')
    hph = asym_ph.GetAsymmetry('ll')
    hpm.SetTitle('p_{T} shift uncertainty on A_{LL} #pi^{+}')
    hpm.SetMarkerStyle(21)
    # hpl.Fit('pol0', 'q')
    hpm.Draw('e1')
    gpl = ROOT.TGraphErrors(hpl)
    gph = ROOT.TGraphErrors(hph)
    for i in range(hpm.GetNbinsX()):
        gpl.SetPoint(i, hpl.GetBinCenter(i+1)-0.02, hpl.GetBinContent(i+1))
        gph.SetPoint(i, hph.GetBinCenter(i+1)+0.02, hph.GetBinContent(i+1))
    gpl.Draw('p')
    gph.Draw('p')
    
    for g in (gml, gmh):
        g.SetMarkerStyle(24)
    for g in (gpl, gph):
        g.SetMarkerStyle(25)
    
    for h in (hmm,hpm):
        h.GetXaxis().SetTitle('p_{T}(#pi)/p_{T}(jet)')
        h.GetYaxis().SetRangeUser(-0.06, 0.08)
    
    raw_input('wait here:')
    c.Print('pt_shift_uncertainty_asymmetry.png')
    
    print '\nπ- shift asymmetry'
    for bin in range(1, hmm.GetNbinsX()+1):
        low = hml.GetBinContent(bin)
        mid = hmm.GetBinContent(bin)
        high = hmh.GetBinContent(bin)
        diffl = abs(mid-low)
        diffh = abs(mid-high)
        diff = (diffl+diffh)/2
        print '[%.2f-%.2f]  % .3f low, %.3f high, %.3f avg' % \
            (hmm.GetBinLowEdge(bin), hmm.GetBinLowEdge(bin+1),diffl,diffh,diff)
    
    print '\nπ+ shift asymmetry'
    for bin in range(1, hpm.GetNbinsX()+1):
        low = hpl.GetBinContent(bin)
        mid = hpm.GetBinContent(bin)
        high = hph.GetBinContent(bin)
        diffl = abs(mid-low)
        diffh = abs(mid-high)
        diff = (diffl+diffh)/2
        print '[%.2f-%.2f]  % .3f low, %.3f high, %.3f avg' % \
            (hpm.GetBinLowEdge(bin), hpm.GetBinLowEdge(bin+1),diffl,diffh,diff)


def pt_shift_uncertainty():
    # these values are actually obtained from inverting the staszak fit params
    measured = [10.31, 12.91, 15.58, 19.06, 23.20, 28.45]
    
    # http://cyclotron.tamu.edu/star/2006Jets/sep28_2007/
    measured = [10.37, 12.76, 15.70, 19.31, 23.75, 29.21]
    
    # these values are from staszak
    corrected = [10.0589, 12.1506, 14.277, 17.0064, 20.204, 24.18]
    
    # these values are from trent
    # corrected = [9.53, 11.50, 13.91, 16.78, 20.25, 24.49]
    
    shift = map(lambda c,m: c-m, corrected, measured)
    
    ## uncertainties
    shift_stats = [0.17, 0.14, 0.14, 0.06, 0.09, 0.09]
    detsim = [0.40, 0.46, 0.55, 0.66, 0.80, 0.97]
    had_ue = [0.57, 0.76, 0.81, 0.86, 0.98, 1.17]
    
    ## yes, these are the same
    sum_quadrature = lambda *args: math.sqrt( sum([a**2 for a in args]) )
    total_uncertainty = map(sum_quadrature, shift_stats, detsim, had_ue)
    total_uncertainty = map(sum_quadrature, detsim, had_ue)
    # total_uncertainty = [0.72, 0.9, 0.99, 1.09, 1.27, 1.52]
    
    graph = ROOT.TGraphErrors(len(shift),
        array('d', measured),
        array('d', shift),
        array('d', [0 for i in range(len(shift))]),
        array('d', shift_stats)
    )
    
    graph_syst = ROOT.TGraphErrors(len(shift),
        array('d', measured),
        array('d', shift),
        array('d', [0 for i in range(len(shift))]),
        array('d', total_uncertainty)
    )
    
    c = graphics.canvas1()
    ROOT.gStyle.SetOptDate(0)
    ROOT.gPad.SetTopMargin(0.02)
    ROOT.gPad.SetRightMargin(0.02)
    
    graph.SetMarkerStyle(20)
    graph.GetYaxis().SetTitle('(corrected - measured) p_{T}')
    graph.GetXaxis().SetTitle('measured p_{T} [GeV/c]')
    graph.SetTitle('')
    graph.Draw('ap')
    
    # the true range in measured pT used to obtain the corrected values
    # graph.GetXaxis().SetRangeUser(9.30, 32.22)
    
    # the range of accepted measured jet pT used in the analysis
    graph.GetXaxis().SetRangeUser(10.0, 30.0)
    graph.GetYaxis().SetRangeUser(-7.0, 1.0)
    
    graph_syst.SetFillColor(16)
    graph_syst.SetFillStyle(3007)
    graph_syst.Draw('[] same')
    
    ## shall i plot these?
    shift_small = map(lambda s,e: s+e, shift, total_uncertainty)
    shift_large = map(lambda s,e: s-e, shift, total_uncertainty)
    
    # shifted = lambda pt: 1.538 + (0.8439-1)*pt - 0.001691*pt*pt
    # fit = ROOT.TGraph(200,
    #     array('d', [10+0.1*i for i in range(200)]),
    #     array('d', [shifted(10+0.1*i) for i in range(200)])
    # )
    # fit.Draw('same')
    
    graph.Fit('pol2')
    graphics.maybe_save()

def nlo_graphs(year, charge):
    from analysis.asym import theoryCurves
    def make_graph(num, denom):
        g = ROOT.TGraphErrors(len(num))
        for i,n in enumerate(num):
            d = denom[i]
            val = n.y / d.y
            err = val*math.sqrt((n.sys/n.y)**2 + (d.sys/d.y)**2)
            g.SetPoint(i, n.x, val)
            g.SetPointError(i, 0, 0)
        return g
    
    if year == 2005 and charge in (-1, 'minus'):
        graphs = {
            'STD': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_std, 
                analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
            'ZERO': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_zero, 
                analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
            'MAX': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_max, 
                analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
            'MIN': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_min, 
                analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),    
            'GSC': theoryCurves(analysis.asym.werner_minus_dss_cteqm5_gsc, 
                analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
            'DSSV':ROOT.TGraphErrors(15,
                array('d', [float(i) for i in range(1,16)]),
                array('d', [
                    -2.685E-05, -4.188E-05, -8.716E-05, -1.562E-04, -1.841E-04,
                    -1.420E-04, -5.366E-05,  6.018E-05,  1.274E-04,  1.235E-04,
                     7.536E-06, -2.645E-04, -6.844E-04, -1.261E-03, -1.999E-03
                ]),
                array('d', [0 for i in range(15)]),
                array('d', [
                    3.011E-05, 3.331E-05, 1.549E-04, 3.664E-04, 5.990E-04,
                    8.074E-04, 9.963E-04, 1.168E-03, 1.345E-03, 1.513E-03,
                    1.684E-03, 1.851E-03, 1.990E-03, 2.117E-03, 2.223E-03
                ]))
        }
    elif year == 2005 and charge in (1, 'plus'):
        graphs = {
            'STD': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_std, 
                analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
            'ZERO': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_zero, 
                analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
            'MAX': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_max, 
                analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
            'MIN': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_min, 
                analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
            'GSC': theoryCurves(analysis.asym.werner_plus_dss_cteqm5_gsc, 
                analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
            'DSSV': ROOT.TGraphErrors(15,
                array('d', [float(i) for i in range(1,16)]),
                array('d', [
                    -3.800E-05, -3.658E-05,  7.296E-07,  1.634E-04,  6.061E-04,
                     1.404E-03,  2.588E-03,  4.144E-03,  5.984E-03,  8.075E-03,
                     1.031E-02,  1.271E-02,  1.528E-02,  1.779E-02,  2.039E-02
                ]),
                array('d', [0 for i in range(15)]),
                array('d', [
                    1.821E-05, 8.788E-05, 3.054E-04, 6.487E-04, 1.038E-03, 
                    1.406E-03, 1.758E-03, 2.090E-03, 2.417E-03, 2.745E-03, 
                    3.080E-03, 3.430E-03, 3.750E-03, 4.090E-03, 4.390E-03
                ]))
        }
    elif year == 2006 and charge in (-1, 'minus'):
        graphs = {
            'STD': make_graph(analysis.deflorian.minus.std, 
                analysis.deflorian.minus.mrst),
            'DSSV': make_graph(analysis.deflorian.minus.dssv, 
                analysis.deflorian.minus.mrst),
            'GSC': make_graph(analysis.deflorian.minus.gsc, 
                analysis.deflorian.minus.mrst),
        }
    elif year == 2006 and charge in (1, 'plus'):
        graphs = {
            'STD': make_graph(analysis.deflorian.plus.std, 
                analysis.deflorian.plus.mrst),
            'DSSV': make_graph(analysis.deflorian.plus.dssv, 
                analysis.deflorian.plus.mrst),
            'GSC': make_graph(analysis.deflorian.plus.gsc, 
                analysis.deflorian.plus.mrst),
        }
    return graphs


def make_fun(graph):
    return lambda x,p: graph.Eval(x[0])

def confidence_levels():
    args = {'m05':(2005,-1), 'p05':(2005,1), 'm06':(2006,-1), 'p06':(2006,1)}
    
    # TF1 versions of all theory curves
    fit = {}
    for prefix,arg in args.items():
        b = (arg[0] == 2005) and (ptbins[0],ptbins[-1]) or (zbins[0],zbins[-1])
        for name,graph in nlo_graphs(*arg).items():
            tf1 = ROOT.TF1(prefix+name, make_fun(graph), b[0], b[1], 1)
            fit.setdefault(prefix,{})[name] = tf1
    
    # graphs of final results
    f05 = ROOT.TFile('/Users/kocolosk/data/run5/final_result.root')
    asym_p = AsymmetryGenerator('asym_p', bins=zbins, key='z_away2')
    asym_m = AsymmetryGenerator('asym_m', bins=zbins, key='z_away2')
    db = sqlite.connect('/Users/kocolosk/data/analysis.db')
    dbc = db.cursor()
    polarizations = Polarizations.Final
    allFiles = glob(histDir + '/chargedPions_*.hist.root')
    for fname in allFiles[:]:
        run = getRun(fname)
        if run in runlist:
            print fname, run
            mgr = HistogramManager(ROOT.TFile(fname), keys=['z_away2'])
            uu,ud,du,dd = analysis.util.scaler_counts(dbc, run)
            pol = polarizations[analysis.util.fill(run)]
            asym_p.FillFromHistogramManager(mgr, 'jetpatch', 1, uu,ud,du,dd, \
                pol.py,pol.pb)
            asym_m.FillFromHistogramManager(mgr, 'jetpatch', -1, uu,ud,du,dd, \
                pol.py,pol.pb)
    
    # including the point-to-point systematic uncertainties
    sum_quad = lambda *args: math.sqrt( sum([a**2 for a in args]) )
    syst05 = systematic_uncertainties_run5()
    syst06 = systematic_uncertainties_run6()
    hm05 = f05.Get('final_minus')
    hp05 = f05.Get('final_plus')
    hm06 = asym_m.GetAsymmetry('ll')
    hp06 = asym_p.GetAsymmetry('ll')
    # for bin in range(1,5):
    #     hm06.SetBinError(bin, sum_quad(hm06.GetBinError(bin), syst06['minus'][bin-1]))
    #     hp06.SetBinError(bin, sum_quad(hp06.GetBinError(bin), syst06['plus'][bin-1]))
    
    result_graph = {
        'm05': ROOT.TGraphAsymmErrors(hm05),
        'p05': ROOT.TGraphAsymmErrors(hp05),
        'm06': ROOT.TGraphAsymmErrors(hm06),
        'p06': ROOT.TGraphAsymmErrors(hp06)
    }
    
    for i in range(4):
        ml = sum_quad(hm06.GetBinError(i+1), syst06['minus']['low'][i])
        mh = sum_quad(hm06.GetBinError(i+1), syst06['minus']['high'][i])
        pl = sum_quad(hp06.GetBinError(i+1), syst06['plus']['low'][i])
        ph = sum_quad(hp06.GetBinError(i+1), syst06['plus']['high'][i])
        result_graph['m06'].SetPointError(i, 0., 0., ml, mh)
        result_graph['p06'].SetPointError(i, 0., 0., pl, ph)
    
    for i in range(5):
        ml = sum_quad(hm05.GetBinError(i+1), syst05['minus']['low'][i])
        mh = sum_quad(hm05.GetBinError(i+1), syst05['minus']['high'][i])
        pl = sum_quad(hp05.GetBinError(i+1), syst05['plus']['low'][i])
        ph = sum_quad(hp05.GetBinError(i+1), syst05['plus']['high'][i])
        result_graph['m05'].SetPointError(i, 0., 0., ml, mh)
        result_graph['p05'].SetPointError(i, 0., 0., pl, ph)
    
    chisquare = {}
    prob = {}
    for measurement, graph in result_graph.items():
        chisquare[measurement] = dict()
        prob[measurement] = dict()
        for scenario, f in fit[measurement].items():
            # graph.Fit(f)
            chi2 = graph.Chisquare(f)
            chisquare[measurement][scenario] = chi2
            prob[measurement][scenario] = ROOT.TMath.Prob(chi2, graph.GetN())
    
    def name(short):
        if short == 'MAX':
            return 'GRSV $\\Delta g$ = g'
        elif short == 'MIN':
            return 'GRSV $\\Delta g$ = -g'
        elif short == 'STD':
            return 'GRSV STD'
        elif short == 'ZERO':
            return 'GRSV $\\Delta g$ = 0'
        elif short == 'GSC':
            return 'GS Set C'
        else:
            return short
    
    for x in ('MAX', 'MIN', 'ZERO', 'STD', 'GSC', 'DSSV'):
        print '%s & %.2g & %.2g & %.2g & %.2g \\\\' % (name(x), prob['m05'][x],
            prob['p05'][x], prob['m06'].get(x, -99), prob['p06'].get(x, -99))

def mcasym(spin = 'anyspin'):
    """
    comparison of MC asymmetries for minbias and 137222
    """
    if spin == 'anyspin':
        stitle = ''
    else:
        stitle = '_%(spin)s' % locals()
    
    f = ROOT.TFile(mcasymFile)
    keys = ['z_M030', 'z_P030']
    # keys = ['MIN']
    mgr = HistogramManager(f, keys=keys)
    
    line = ROOT.TLine(zbins[0], 0.0, zbins[-1], 0.0)
    line.SetLineStyle(2)
    ROOT.gStyle.SetErrorX()
    
    smooth_factor = 1
    
    c = graphics.canvas2('Reweighted Asymmetry Differences')
    diffs = {}
    for i,key in enumerate(keys):
        mb_m = mgr[spin]['117001'].tracks_minus[key].Clone()
        mb_p = mgr[spin]['117001'].tracks_plus[key].Clone()
        jp_m = mgr[spin]['137222'].tracks_minus[key].Clone()
        jp_p = mgr[spin]['137222'].tracks_plus[key].Clone()            
        
        opt = i>0 and 'e2 same' or 'e2'
        
        for h in (mb_m, mb_p, jp_m, jp_p):
            # h.Smooth(smooth_factor)
            h.GetXaxis().SetTitle('p_{T}(#pi)/p_{T}(jet)')
            h.GetXaxis().SetRangeUser(zbins[0], zbins[-1])
            h.GetYaxis().SetRangeUser(-0.06,0.06)
        
        diff_m = jp_m.Clone()
        diff_m.Add(mb_m, -1)
        
        diff_p = jp_p.Clone()
        diff_p.Add(mb_p, -1)
        
        for h in (diff_m, diff_p):
            h.GetYaxis().SetRangeUser(-0.06, 0.06)
            h.GetXaxis().SetRangeUser(zbins[0], zbins[-1])
            h.SetMarkerStyle(20)
        
        cdiff = graphics.canvas2(key)
        
        leg = ROOT.TLegend(0.15, 0.15, 0.6, 0.35)
        leg.SetHeader('Compare triggered asym to:')
        leg.AddEntry(diff_m, 'p_{T} reweighted MB', 'p')
        
        cdiff.cd(1)
        if spin == 'anyspin':
            mtitle = '#pi^{-} JP1 - MB for %(key)s' % locals()
            ptitle = '#pi^{+} JP1 - MB for %(key)s' % locals()            
        else:
            mtitle = '#pi^{-} JP1 - MB for %(key)s -- %(spin)s processes only' \
                % locals()
            ptitle = '#pi^{+} JP1 - MB for %(key)s -- %(spin)s processes only' \
                % locals()
        
        mb_m.SetTitle(mtitle)
        mb_m.Draw('e1')
        jp_m.Draw('e1 same')
        line.Draw()
        leg.Draw()
        
        cdiff.cd(2)
        mb_p.SetTitle(ptitle)
        mb_p.Draw('e1')
        jp_p.Draw('e1 same')
        line.Draw()
        
        opt2 = i>0 and 'e1 same' or 'e1'
        
        c.cd(1)
        line.Draw()
        diff_m.Draw(opt2)
        
        c.cd(2)
        line.Draw()
        diff_p.Draw(opt2)
        
        diffs[key] = (diff_m, diff_p)
        raw_input('pause')
        del cdiff
    
    # now report the asymmetry difference for each scenario
    for key in keys:
        diff_m, diff_p = diffs[key]
        print '\nπ- asymmetry differences for %(key)s' % locals()
        for bin in range(1, diff_m.GetNbinsX()+1):
            print '[%.2f-%.2f]  % .4f ± %.4f' % (
                diff_m.GetBinLowEdge(bin), diff_m.GetBinLowEdge(bin+1), 
                diff_m.GetBinContent(bin), diff_m.GetBinError(bin) )
                
        print '\nπ+ asymmetry differences for %(key)s' % locals()
        for bin in range(1, diff_m.GetNbinsX()+1):
            print '[%.2f-%.2f]  % .4f ± %.4f' % (
                diff_p.GetBinLowEdge(bin), diff_p.GetBinLowEdge(bin+1), 
                diff_p.GetBinContent(bin), diff_p.GetBinError(bin) )
    
    c.Update()
    raw_input('wait here:')
    
    c.Print('mcasym_diff_reweight%(stitle)s.png' % locals())


def reweighted_mb_mcasym():
    ROOT.gStyle.SetOptDate(0)
    
    f1 = ROOT.TFile(analysis.plots.spin2008.mcasymFile)
    f2 = ROOT.TFile(mcasymFile)
    
    mgr1 = HistogramManager(f1, keys=['STDw', 'DSSVw'])
    mgr2 = HistogramManager(f2, keys=['z_M030', 'z_P030'])
    
    color = {
        'STDw': 15,
        'DSSVw': ROOT.kBlack,
        'z_M030': ROOT.kMagenta,
        'z_P030': ROOT.kMagenta
    }
    
    name = {
        'STDw': 'GRSV STD',
        'DSSVw': 'DSSV',
        'z_M030': '\Delta G = -0.3',
        'z_P030': '\Delta G = \pm 0.3'
    }
    
    c = graphics.canvas2()
    line = ROOT.TLine(zbins[0], 0.0, zbins[-1], 0.0)
    line.SetLineStyle(2)
    leg = ROOT.TLegend(0.4, 0.4, 0.88, 0.7)
    
    begin = True
    for key in ('STDw', 'DSSVw', 'z_M030', 'z_P030'):
        mgr = key in ('STDw', 'DSSVw') and mgr1 or mgr2
        hm = mgr['anyspin']['117001'].tracks_minus[key]
        hp = mgr['anyspin']['117001'].tracks_plus[key]
        
        # take a look at JP instead
        # try:
        #     hm = mgr['anyspin']['jetpatch'].tracks_minus[key]
        #     hp = mgr['anyspin']['jetpatch'].tracks_plus[key]
        # except KeyError:
        #     hm = mgr['anyspin']['137222'].tracks_minus[key]
        #     hp = mgr['anyspin']['137222'].tracks_plus[key]
        
        for h in (hm,hp):
            h.SetLineColor(color[key])
            h.SetMarkerColor(color[key])
            h.SetMarkerStyle(20)
            h.GetXaxis().SetRangeUser(zbins[0], zbins[-1])
            h.GetYaxis().SetRangeUser(-0.005, 0.075)
            h.SetXTitle('z')
            h.SetYTitle('A_{LL}')
            h.GetYaxis().SetTitleSize(0.06)
            h.GetYaxis().SetTitleOffset(1.0)
            h.SetTitle('')
        
        if key in ('z_M030', 'z_P030'):
            h.SetMarkerStyle(21)
        
        if key != 'z_M030':
            print key
            leg.AddEntry(hm, name[key], 'lp')
        
        latex = ROOT.TLatex()
        latex.SetTextSize(0.23)
        
        opt = begin and 'e1' or 'e1 same'
        c.cd(1)
        hm.DrawCopy(opt)
        line.Draw()
        leg.Draw()
        latex.DrawLatex(0.05 + zbins[0], 0.054, '#pi^{-}')
        c.cd(2)
        hp.DrawCopy(opt)
        line.Draw()
        latex.DrawLatex(0.05 + zbins[0], 0.054, '#pi^{+}')
        begin = False
    
    for i in (1, 2):
        pad = c.cd(i)
        pad.SetTopMargin(0.03)
        pad.SetBottomMargin(0.1)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.02)
    raw_input('wait')
    c.Print('mb_mcasym_2006.png')

def syst_mcasym_2006():
    f = ROOT.TFile(mcasymFile)
    mgr = HistogramManager(f, keys=['z_M030', 'z_P030'])
    
    mb_m_m = mgr['anyspin']['117001'].tracks_minus['z_M030']
    mb_m_p = mgr['anyspin']['117001'].tracks_minus['z_P030']
    mb_p_m = mgr['anyspin']['117001'].tracks_plus['z_M030']
    mb_p_p = mgr['anyspin']['117001'].tracks_plus['z_P030']

    jp_m_m = mgr['anyspin']['137222'].tracks_minus['z_M030']
    jp_m_p = mgr['anyspin']['137222'].tracks_minus['z_P030']
    jp_p_m = mgr['anyspin']['137222'].tracks_plus['z_M030']
    jp_p_p = mgr['anyspin']['137222'].tracks_plus['z_P030']
    
    diff_m_m = jp_m_m
    diff_m_p = jp_m_p
    diff_p_m = jp_p_m
    diff_p_p = jp_p_p
    
    diff_m_m.Add(mb_m_m, -1)
    diff_m_p.Add(mb_m_p, -1)
    diff_p_m.Add(mb_p_m, -1)
    diff_p_p.Add(mb_p_p, -1)
    
    # for bin in range(4,8):
    #     print diff_p_p.GetBinCenter(bin), jp_p_p.GetBinContent(bin), mb_p_p.GetBinContent(bin), diff_p_p.GetBinContent(bin)
    #     diff_m_m.SetBinContent(bin, jp_m_m.GetBinContent(bin) - mb_m_m.GetBinContent(bin))
    #     diff_m_p.SetBinContent(bin, jp_m_p.GetBinContent(bin) - mb_m_p.GetBinContent(bin))
    #     diff_p_m.SetBinContent(bin, jp_p_m.GetBinContent(bin) - mb_p_m.GetBinContent(bin))
    #     diff_p_p.SetBinContent(bin, jp_p_p.GetBinContent(bin) - mb_p_p.GetBinContent(bin))
    
    gm, gem = make_syst_trig_graphs(diff_m_m, diff_m_p, 4, 4)
    gp, gep = make_syst_trig_graphs(diff_p_m, diff_p_p, 4, 4)
    
    h = diff_p_p
    for i in range(4):
        print '%.2f - %.2f  & -%.4f +%.4f & -%.4f +%.4f \\\\' % \
            (h.GetBinLowEdge(i+4), h.GetBinLowEdge(i+5),
                max(gm.GetErrorYlow(i), gem.GetErrorY(i)),
                max(gm.GetErrorYhigh(i), gem.GetErrorY(i)),
                max(gp.GetErrorYlow(i), gep.GetErrorY(i)),
                max(gp.GetErrorYhigh(i), gep.GetErrorY(i)),
            )
    
    c = graphics.canvas2()
    shell = ROOT.TH2D('shell', '', 1, 0.2, 1.0, 1, -0.008, 0.024)
    shell.SetXTitle('z')
    shell.SetYTitle('A_{LL}^{JP} - A_{LL}^{MB}')
    shell.GetYaxis().SetTitleSize(0.055)
    shell.GetYaxis().SetTitleOffset(1.2)
    line = ROOT.TLine(zbins[0], 0.0, zbins[-1], 0.0)
    line.SetLineStyle(2)
    
    for g in (gem, gep):
        g.SetMarkerColor(15)
        g.SetFillColor(15)
    
    for g in (gm, gp):
        g.SetMarkerStyle(21)
    
    c.cd(1)
    shell.Draw()
    line.Draw()
    gem.Draw('2')
    gm.Draw('[]')
    
    c.cd(2)
    shell.Draw()
    line.Draw()
    gep.Draw('2')
    gp.Draw('[]')
    
    for i in (1, 2):
        pad = c.cd(i)
        pad.SetTopMargin(0.03)
        pad.SetBottomMargin(0.1)
        pad.SetLeftMargin(0.15)
        pad.SetRightMargin(0.02)
    raw_input('wait')
    c.Print('syst_mcasym_2006.png')


def syst_mcasym_2005():
    # making graphs by eye
    gm = ROOT.TGraphAsymmErrors(5)
    gp = ROOT.TGraphAsymmErrors(5)
    gem = ROOT.TGraphErrors(5)
    gep = ROOT.TGraphErrors(5)
    
    for g in (gm,gp,gem,gep):
        for i in range(5):
            g.SetPoint(i, (ptbins[i+1] + ptbins[i])/2, 0.0)
    
    
    gm.SetPointError(0, 0, 0, 0.0015, 0.0076)
    gm.SetPointError(1, 0, 0, 0.0014, 0.0101)
    gm.SetPointError(2, 0, 0, 0.0011, 0.0105)
    gm.SetPointError(3, 0, 0, 0.0008, 0.0097)
    gm.SetPointError(4, 0, 0, 0.0020, 0.0070)
    
    gp.SetPointError(0, 0, 0, 0.0008, 0.0021)
    gp.SetPointError(1, 0, 0, 0.0008, 0.0149)
    gp.SetPointError(2, 0, 0, 0.0008, 0.0202)
    gp.SetPointError(3, 0, 0, 0.0008, 0.0228)
    gp.SetPointError(4, 0, 0, 0.0008, 0.0162)
    
    gem.SetPointError(0, 0.15, 0.0010)
    gem.SetPointError(1, 0.15, 0.0023)
    gem.SetPointError(2, 0.15, 0.0024)
    gem.SetPointError(3, 0.15, 0.0041)
    gem.SetPointError(4, 0.15, 0.0052)
    
    gep.SetPointError(0, 0.15, 0.0015)
    gep.SetPointError(1, 0.15, 0.0021)
    gep.SetPointError(2, 0.15, 0.0025)
    gep.SetPointError(3, 0.15, 0.0041)
    gep.SetPointError(4, 0.15, 0.0051)
    
    print '2.00 - 3.18  & -%.4f +%.4f & -%.4f +%.4f \\\\' % (0.0015, 0.0076, 0.0015, 0.0021)
    print '3.18 - 4.56  & -%.4f +%.4f & -%.4f +%.4f \\\\' % (0.0023, 0.0101, 0.0021, 0.0149)
    print '4.56 - 6.32  & -%.4f +%.4f & -%.4f +%.4f \\\\' % (0.0024, 0.0105, 0.0025, 0.0202)
    print '6.32 - 8.80  & -%.4f +%.4f & -%.4f +%.4f \\\\' % (0.0041, 0.0097, 0.0041, 0.0228)
    print '8.80 - 12.84 & -%.4f +%.4f & -%.4f +%.4f \\\\' % (0.0052, 0.0070, 0.0051, 0.0162)
    
    c = graphics.canvas2()
    shell = ROOT.TH2D('shell', '', 1, ptbins[0], ptbins[-1], 1, -0.008, 0.024)
    shell.SetXTitle('p_{T}')
    shell.SetYTitle('A_{LL}^{JP} - A_{LL}^{MB}')
    shell.GetYaxis().SetTitleSize(0.055)
    shell.GetYaxis().SetTitleOffset(1.2)
    line = ROOT.TLine(ptbins[0], 0.0, ptbins[-1], 0.0)
    line.SetLineStyle(2)
    
    for g in (gem, gep):
        g.SetMarkerColor(15)
        g.SetFillColor(15)
    
    for g in (gm, gp):
        g.SetMarkerStyle(21)
    
    c.cd(1)
    shell.Draw()
    line.Draw()
    gem.Draw('2')
    gm.Draw('[]')
    
    c.cd(2)
    shell.Draw()
    line.Draw()
    gep.Draw('2')
    gp.Draw('[]')
    
    for i in (1, 2):
        pad = c.cd(i)
        pad.SetTopMargin(0.03)
        pad.SetBottomMargin(0.1)
        pad.SetLeftMargin(0.15)
        pad.SetRightMargin(0.02)
    
    raw_input('wait')
    c.Print('syst_mcasym_2005.png')

def make_syst_trig_graphs(low, high, nbins, offset):
    g1 = ROOT.TGraphAsymmErrors(nbins)
    g2 = ROOT.TGraphErrors(nbins)
    
    i = 0
    for bin in range(offset, nbins+offset):
        g1.SetPoint(i, low.GetBinCenter(bin), 0.0)
        g2.SetPoint(i, low.GetBinCenter(bin), 0.0)
        g1.SetPointError(i, 0., 0., max(0.0008, low.GetBinContent(bin)), max(0.0008, high.GetBinContent(bin)))
        g2.SetPointError(i, 0.015, low.GetBinError(bin))
        i += 1
    return (g1,g2)

