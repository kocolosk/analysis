# encoding: utf-8

import os
import math
from glob import glob

import ROOT

import analysis.deflorian

from analysis.asym import AsymmetryGenerator, ScalarCounts, Polarizations
from analysis.histos import HistogramManager, shifted
from analysis.plots import pid_calibration
from analysis.runlists import long2_run6 as runlist
from analysis.runlists import transverse_run6 as transverse_runlist
from analysis.util import getRun, fillList, hadd_interactive
from analysis.plots import graphics
from analysis.plots.spin2008 import systematic_uncertainties

histDir         = '/Users/kocolosk/data/run6/spin2008/hist'
transHistDir    = '/Users/kocolosk/data/run6/spin2008/hist-transverse'
mcasymFile      = '/Users/kocolosk/data/run6/simu/mcasym_10.cphist.root'
datamcFile      = '/Users/kocolosk/data/run6/spin2008/merged.cphist.root'

zbins = [
    # 0.075, 0.125, ## these two are biased b/c of the track pT cut
    0.2, 0.3, 0.45, 0.65, 1.0
]

def run6_result_theory_curves():
    """
    result plot showing asymmetries versus z
    """    
    asym_p = AsymmetryGenerator('asym_p', bins=zbins, key='z_away2')
    asym_m = AsymmetryGenerator('asym_m', bins=zbins, key='z_away2')
    
    scalars = ScalarCounts(os.environ['STAR'] + 
        '/StRoot/StSpinPool/StTamuRelLum/inputs/run6.txt')
    
    polarizations = Polarizations.Final
    
    allFiles = glob(histDir + '/chargedPions_*.hist.root')
    for fname in allFiles[:]:
        run = getRun(fname)
        if run in runlist:
            print fname, run
            mgr = HistogramManager(ROOT.TFile(fname), ['z_away2'])
            
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
    
    syst = systematic_uncertainties()
    syst_m = ROOT.TGraphErrors(len(zbins))
    syst_p = ROOT.TGraphErrors(len(zbins))
    
    for g in (syst_m,syst_p):
        g.SetMarkerColor(15)
        g.SetFillColor(15)
        g.SetPoint(len(zbins), 1.0, 0)
        g.SetPointError(len(zbins), 0., 0.)
        g.GetXaxis().SetRangeUser(zbins[0], zbins[-1])
    
    hm = asym_m.GetAsymmetry('ll')
    hp = asym_p.GetAsymmetry('ll')
    for i in range(len(zbins)-1):
        z = (zbins[i] + zbins[i+1])/2
        syst_m.SetPoint(i, z, hm.GetBinContent(i+1))
        syst_m.SetPointError(i, 0.01, syst['minus'][i])
        syst_p.SetPoint(i, z, hp.GetBinContent(i+1))
        syst_p.SetPointError(i, 0.01, syst['plus'][i])
    
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
            g.SetPointError(i, 0, err)
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
        nlo['DSSV'].SetLineStyle(2)
        nlo['GSC'].SetLineColor(ROOT.kMagenta)
        nlo['DSSV'].SetLineColor(ROOT.kGreen+2)
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
    prelim.DrawText(0.42,-0.041,"STAR Preliminary")
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
    prelim.DrawText(0.42,-0.041,"STAR Preliminary")
    latex.DrawLatex(0.05 + zbins[0], 0.044, '#pi^{+}')
    [ g.Draw('3') for g in nlo_p.values() ]
    [ g.Draw('lx') for g in nlo_p.values() ]
    line.Draw()
    leg.Draw()
    syst_p.Draw('2p')
    hp.Draw('e1 same')
    
    graphics.maybe_save(c)
