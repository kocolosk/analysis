# encoding: utf-8
# i need a better header template :-/
# this module is supposed to collect all the plots i'll be showing at SPIN 2008

import os
import math
from glob import glob

import ROOT

from .asym import AsymmetryGenerator, ScalarCounts, Polarizations
from .histos import HistogramManager
from .plots import pid_calibration
from .runlists import long2_run6 as runlist
from .util import getRun, fillList, hadd_interactive
from . import graphics

histDir = '/Users/kocolosk/data/run6/spin2008/hist'
transHistDir = '/Users/kocolosk/data/run6/spin2008/hist-transverse'
zbins = [
    0.075, 0.125, ## these two are biased b/c of the track pT cut
    0.2, 0.3, 0.45, 0.65, 1.0
]

def result():
    """
    result plot showing asymmetries versus z
    """    
    asym_p = AsymmetryGenerator('asym_p', bins=zbins, key='z_away2')
    asym_m = AsymmetryGenerator('asym_m', bins=zbins, key='z_away2')
    
    scalars = ScalarCounts(os.environ['STAR'] + 
        '/StRoot/StSpinPool/StTamuRelLum/inputs/run6.txt')
    
    polarizations = Polarizations.Final
    
    ## generate the asymmetries
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
    
    def makeGraph(h, xshift):
        g = ROOT.TGraphErrors(h)
        for bin in range(1,h.GetNbinsX()+1):
            g.SetPoint(bin-1, h.GetBinCenter(bin)+xshift, h.GetBinContent(bin))
        return g
    
    
    hp = asym_p.GetAsymmetry('ll')
    gp = makeGraph(hp, 0.01)
    gp.SetMarkerStyle(21)
    
    hm = asym_m.GetAsymmetry('ll')
    gm = makeGraph(hm, -0.01)
    gm.SetMarkerColor(ROOT.kRed)
    gm.SetMarkerStyle(21)
    gm.SetLineColor(ROOT.kBlack)
    
    ### I HATE TGraphErrors
    syst = ROOT.TGraphErrors(2*hp.GetNbinsX())
    for i in range(hp.GetNbinsX()):
        syst.SetPoint(i, hp.GetBinCenter(i+1)+0.01, hp.GetBinContent(i+1))
        syst.SetPointError(i, 0.004, 0.01)
    for i in range(hp.GetNbinsX(), hp.GetNbinsX()+hm.GetNbinsX()):
        bin = i - hp.GetNbinsX() + 1
        syst.SetPoint(i, hm.GetBinCenter(bin)-0.01, hm.GetBinContent(bin))
        syst.SetPointError(i, 0.004, 0.01)
    
    syst.GetXaxis().Set(1,0.0,1.0)
    syst.GetXaxis().SetTitle('p_{T}(#pi) / p_{T}(jet)')
    syst.GetYaxis().SetRangeUser(-0.06, 0.12)
    syst.GetYaxis().SetTitleOffset(1.2)
    syst.GetYaxis().SetTitle('A_{LL}')
    syst.SetTitle('')
    syst.SetFillColor(16)
    syst.SetLineColor(16)
    
    line = ROOT.TLine(0.0, 0.0, 1.0, 0.0)
    line.SetLineStyle(2)
    
    leg = ROOT.TLegend(0.2,0.7,0.3,0.87)
    leg.AddEntry(gm, '#pi^{-}', 'p')
    leg.AddEntry(gp, '#pi^{+}', 'p')
    leg.SetBorderSize(0)
    
    c1 = graphics.canvas1('A_{LL}')
    syst.Draw('a2')
    syst.Draw('p')
    gp.Draw('p')
    gm.Draw('p')
    leg.Draw()
    
    latex = ROOT.TLatex()
    latex.DrawLatex(0.03,0.125,
    " #vec{p} + #vec{p} #rightarrow jet + #pi + X at #sqrt{s}=200 GeV \
             -1< #eta^{#pi}< 1 ")
    
    prelim = ROOT.TText()
    prelim.SetTextColor(44)
    prelim.DrawText(0.3,0.065,"2006 STAR Preliminary");
    
    raw_input('wait here:')
    


def systematic_uncertainties():
    """tabulates sources of uncertainty and sums them in quadrature"""
    
    pid_contamination = 0.10
    pid_asym_m = [
        0.051,      # [0.07-0.12]   0.051 ± 0.038
        0.017,      # [0.12-0.20]  -0.017 ± 0.016
        0.032,      # [0.20-0.30]  -0.032 ± 0.016
        0.023,      # [0.30-0.45]  -0.006 ± 0.023
        0.042,      # [0.45-0.65]  -0.031 ± 0.042
        0.089       # [0.65-1.00]   0.089 ± 0.085
    ]
    pid_asym_p = [
        0.036,      # [0.07-0.12]   0.005 ± 0.036
        0.015,      # [0.12-0.20]   0.006 ± 0.015
        0.015,      # [0.20-0.30]  -0.006 ± 0.015
        0.020,      # [0.30-0.45]   0.018 ± 0.020
        0.038,      # [0.45-0.65]  -0.038 ± 0.032
        0.142       # [0.65-1.00]   0.142 ± 0.059
    ]
    
    beam_vector = 0.0102
    asigma_m = [
        0.107,      # [0.07-0.12]   0.107 ± 0.099
        0.056,      # [0.12-0.20]   0.039 ± 0.056
        0.067,      # [0.20-0.30]   0.050 ± 0.067
        0.102,      # [0.30-0.45]  -0.059 ± 0.102
        0.266,      # [0.45-0.65]  -0.266 ± 0.179
        0.313       # [0.65-1.00]   0.030 ± 0.313
    ]
    asigma_p = [
        0.099,      # [0.07-0.12]  -0.094 ± 0.099
        0.055,      # [0.12-0.20]   0.037 ± 0.055
        0.065,      # [0.20-0.30]   0.007 ± 0.065
        0.096,      # [0.30-0.45]  -0.033 ± 0.096
        0.161,      # [0.45-0.65]   0.096 ± 0.161
        0.287       # [0.65-1.00]  -0.168 ± 0.287
    ]
    
    relative_luminosity = 9.4e-4
    
    minus = [0.0 for bin in zbins[:-1]]
    plus = [0.0 for bin in zbins[:-1]]
    
    for i in range(len(zbins)-1):
        minus[i] = math.sqrt(
            pow(relative_luminosity, 2) + 
            pow(pid_contamination*pid_asym_m[i], 2) +
            pow(beam_vector*asigma_m[i], 2)
        )
        plus[i] = math.sqrt(
            pow(relative_luminosity, 2) +
            pow(pid_contamination*pid_asym_p[i], 2) + 
            pow(beam_vector*asigma_p[i], 2)
        )
    
    return {'minus':minus, 'plus':plus}


def ssa():
    """
    plots single-spin asymmetries versus fill and z, saves them in PNG format
    """
    asym_zp = AsymmetryGenerator('asym_zp', bins=zbins, key='z_away2')
    asym_zm = AsymmetryGenerator('asym_zm', bins=zbins, key='z_away2')
    asym_fp = {}
    asym_fm = {}
    
    fills = [int(f) for f in fillList(runlist)]
    for fill in fills:
        asym_fp[fill] = AsymmetryGenerator('p%d' % fill, bins=[1,-0.5,0.5], 
            key='one')
        asym_fm[fill] = AsymmetryGenerator('m%d' % fill, bins=[1,-0.5,0.5], 
            key='one')
    
    scalars = ScalarCounts(os.environ['STAR'] + 
        '/StRoot/StSpinPool/StTamuRelLum/inputs/run6.txt')
    
    polarizations = Polarizations.Final
    
    allFiles = glob(histDir + '/chargedPions_*.hist.root')
    for fname in allFiles[:]:
        run = getRun(fname)
        if run in runlist:
            print fname, run
            mgr = HistogramManager(ROOT.TFile(fname), ['z_away2','one'])
            
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
            
            asym_zp.FillFromHistogramManager(mgr, 'jetpatch', 1, uu,ud,du,dd,
                pol.py,pol.pb)
            asym_zm.FillFromHistogramManager(mgr, 'jetpatch', -1, uu,ud,du,dd,
                pol.py,pol.pb)
            asym_fp[bin7.fill].FillFromHistogramManager(mgr, 'jetpatch', 1, 
                uu,ud,du,dd, pol.py,pol.pb)
            asym_fm[bin7.fill].FillFromHistogramManager(mgr, 'jetpatch', -1, 
                uu,ud,du,dd, pol.py,pol.pb)
            
    title = {
        'ly':'Yellow Beam', 
        'lb':'Blue Beam', 
        'ls':'Like-Sign', 
        'us':'Unlike-Sign'
    }
    marker_color = {
        'ly':ROOT.kYellow, 
        'lb':ROOT.kBlue, 
        'ls':ROOT.kRed, 
        'us':ROOT.kBlack
    }
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetErrorX(0.0)
    ROOT.gStyle.SetOptFit(111)
    
    canvases = []
    for key in ('ly','lb','ls','us'):
        final_plus = ROOT.TH1D('final_plus_%s' % key, '', len(fills), 0.5, 
            len(fills)+0.5)
        final_minus = ROOT.TH1D('final_minus_%s' % key, '', len(fills), 0.5, 
            len(fills)+0.5)
        canvases.append([final_minus, final_plus])
        for i,f in enumerate(fills):
            hplus = asym_fp[f].GetAsymmetry(key)
            final_plus.SetBinContent( i+1, hplus.GetBinContent(1) )
            final_plus.SetBinError( i+1, hplus.GetBinError(1) )
            hplus.Delete()
            
            hminus = asym_fm[f].GetAsymmetry(key)
            final_minus.SetBinContent( i+1, hminus.GetBinContent(1) )
            final_minus.SetBinError( i+1, hminus.GetBinError(1) )
            hminus.Delete()
            
        for var,p,m in [('z',asym_zp,asym_zm), ('fill',asym_fp,asym_fm)]:
            c = graphics.canvas2()
            canvases.append(c)
        
            hp = (var == 'z') and p.GetAsymmetry(key) or final_plus
            hp.SetTitle(title[key] + ' SSA for #pi^{+}')
            hp.SetMarkerStyle(21)
    
            hm = (var == 'z') and m.GetAsymmetry(key) or final_minus
            hm.SetTitle(title[key] + ' SSA for #pi^{-}')
            hm.SetMarkerStyle(20)
    
            for h in (hp,hm):
                h.SetMarkerColor(marker_color[key])
                if var == 'fill':
                    h.GetYaxis().SetRangeUser(-0.2, 0.2)
                    h.SetXTitle('fill index')
                elif var == 'z':
                    h.SetXTitle('p_{T}(#pi) / p_{T}(jet)')
                    if key in ('ly', 'lb'):
                        h.GetYaxis().SetRangeUser(-0.05, 0.05)
                    else:
                        h.GetYaxis().SetRangeUser(-0.1, 0.1)
        
            c.cd(1)
            hm.Fit('pol0')
            hm.Draw('e1')
        
            c.cd(2)
            hp.Fit('pol0')
            hp.Draw('e1')
        
            c.Print('%s_ssa_%s.png' % (key, var))
    
    raw_input('wait here:')


def pid_background_asymmetry():
    """
    calculate A_{LL} for particles in the p/K sideband of the dE/dx window
    """
    asym_p = AsymmetryGenerator('asym_p', bins=zbins, key='z_away2_bg')
    asym_m = AsymmetryGenerator('asym_m', bins=zbins, key='z_away2_bg')
    
    scalars = ScalarCounts(os.environ['STAR'] + 
        '/StRoot/StSpinPool/StTamuRelLum/inputs/run6.txt')
    
    polarizations = Polarizations.Final
    
    ## generate the asymmetries
    allFiles = glob(histDir + '/chargedPions_*.hist.root')
    for fname in allFiles[:]:
        run = getRun(fname)
        if run in runlist:
            print fname, run
            mgr = HistogramManager(ROOT.TFile(fname), ['z_away2_bg'])
            
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
    
    c = graphics.canvas2()
    c.SetLogy(0)
    ROOT.gStyle.SetErrorX(0.0)
    
    c.cd(1)
    hm = asym_m.GetAsymmetry('ll')
    hm.SetMarkerStyle(20)
    hm.SetTitle('Background A_{LL} #pi^{-}')
    hm.Fit('pol0', 'q')
    hm.Draw('e1')
    
    c.cd(2)
    hp = asym_p.GetAsymmetry('ll')
    hp.SetTitle('Background A_{LL} #pi^{+}')
    hp.SetMarkerStyle(21)
    hp.Fit('pol0', 'q')
    hp.Draw('e1')
    
    for h in (hm,hp):
        h.GetXaxis().SetTitle('p_{T}(#pi)/p_{T}(jet)')
        h.GetYaxis().SetRangeUser(-0.1, 0.1)
    
    raw_input('wait here:')
    c.Print('pid_background_asymmetry.png')
    
    print '\nπ- Background Asymmetry'
    for bin in range(1, hm.GetNbinsX()+1):
        print '[%.2f-%.2f]  % .3f ± %.3f' % (hm.GetBinLowEdge(bin), 
            hm.GetBinLowEdge(bin+1), hm.GetBinContent(bin), hm.GetBinError(bin))
    
    print '\nπ+ Background Asymmetry'
    for bin in range(1, hp.GetNbinsX()+1):
        print '[%.2f-%.2f]  % .3f ± %.3f' % (hp.GetBinLowEdge(bin), 
            hp.GetBinLowEdge(bin+1), hp.GetBinContent(bin), hp.GetBinError(bin))
    


def pid_background_fraction():
    """simplistic calculation of background contamination fraction"""
    nsig_p = hadd_interactive(histDir, runlist, 'jetpatch', 'anyspin', 'plus', 
        'away2_nSigmaPion')
    nsig_m = hadd_interactive(histDir, runlist, 'jetpatch', 'anyspin', 'minus', 
        'away2_nSigmaPion')
    
    nsig_p.SetTitle('n#sigma(#pi) for #pi^{+}')
    nsig_m.SetTitle('n#sigma(#pi) for #pi^{-}')
    
    print 'Fit summary for π+'
    pfits = pid_calibration(nsig_p)
    print 'Fit summary for π-'
    mfits = pid_calibration(nsig_m)
    c = graphics.canvas2()
    
    c.cd(1)
    nsig_p.Draw()
    [fit.Draw('same') for fit in pfits[1:]]
    
    c.cd(2)
    nsig_m.Draw()
    [fit.Draw('same') for fit in mfits[1:]]
    
    raw_input('wait here:')
    c.Print('pid_background_fraction.png')


def asigma():
    asym_p = AsymmetryGenerator('asym_p', bins=zbins, key='z_away2')
    asym_m = AsymmetryGenerator('asym_m', bins=zbins, key='z_away2')
    
    scalars = ScalarCounts(os.environ['STAR'] + 
        '/StRoot/StSpinPool/StTamuRelLum/inputs/run6.txt')
    
    polarizations = Polarizations.Final
    
    ## generate the asymmetries
    allFiles = glob(transHistDir + '/chargedPions_*.hist.root')
    for fname in allFiles[:10]:
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
    
    c = graphics.canvas2()
    c.SetLogy(0)
    ROOT.gStyle.SetErrorX(0.0)
    
    c.cd(1)
    hm = asym_m.GetAsymmetry('ll')
    hm.SetMarkerStyle(20)
    hm.SetTitle('A_{#sigma} #pi^{-}')
    hm.Fit('pol0', 'q')
    hm.Draw('e1')
    
    c.cd(2)
    hp = asym_p.GetAsymmetry('ll')
    hp.SetTitle('A_{#sigma} #pi^{+}')
    hp.SetMarkerStyle(21)
    hp.Fit('pol0', 'q')
    hp.Draw('e1')
    
    for h in (hm,hp):
        h.GetXaxis().SetTitle('p_{T}(#pi)/p_{T}(jet)')
        h.GetYaxis().SetRangeUser(-0.1, 0.1)
    
    raw_input('wait here:')
    c.Print('asigma.png')
    
    print '\nπ- A_{σ}'
    for bin in range(1, hm.GetNbinsX()+1):
        print '[%.2f-%.2f]  % .3f ± %.3f' % (hm.GetBinLowEdge(bin), 
            hm.GetBinLowEdge(bin+1), hm.GetBinContent(bin), hm.GetBinError(bin))
    
    print '\nπ+ A_{σ}'
    for bin in range(1, hp.GetNbinsX()+1):
        print '[%.2f-%.2f]  % .3f ± %.3f' % (hp.GetBinLowEdge(bin), 
            hp.GetBinLowEdge(bin+1), hp.GetBinContent(bin), hp.GetBinError(bin))
    
    ## conservative estimate of |tan(θY)tan(θB)cos(ϕY-ϕB)| = 0.0102 ± 0.0002


def ffcomp():
    """
    this function still needs some work ... seems very buggy
    """
    from .util import tf1
    from .ff import dss
    from math import sqrt
    
    f = ROOT.TFile('pythia_ff.hist.root')
    ROOT.gStyle.SetOptLogy()
    ROOT.gStyle.SetErrorX()
    charge = 'minus'
    
    line = ROOT.TLine(0,1,1,1)
    line.SetLineStyle(2)
    def dostuff(c, h, f):
        c.cd(1)
        ROOT.gPad.SetLogy()
        scalefactor = f.Integral(0.2, 0.8) / h.Integral(40, 160, 'width')
        h.Scale(scalefactor)
        h.GetXaxis().SetTitle('z')
        ratio = h.Clone()
        h.Draw()
        f.SetLineColor(ROOT.kRed)
        f.Draw('same')
        c.cd(2)
        ROOT.gPad.SetLogy(0)
        ratio.SetTitle('')
        ratio.SetXTitle('')
        ratio.Sumw2()
        ratio.Divide(f)
        for bin in range(1, ratio.GetNbinsX()+1):
            ratio.SetBinError(bin, ratio.GetBinError(bin)*sqrt(scalefactor))
        ratio.GetYaxis().SetRangeUser(0,2)
        ratio.GetYaxis().SetTitle('PYTHIA/DSS')
        ratio.Draw()
        line.Draw('same')
    
    cu = graphics.canvas1e('u+ubar')
    uubar = f.Get('u_%(charge)s' % locals())
    uubar.Add(f.Get('ubar_%(charge)s' % locals()))
    uubar.SetTitle('u+ubar => #pi^{-}')
    dssu = tf1(dss, 0, 1, flavor='u+ubar', charge=charge)
    dostuff(cu, uubar, dssu)
    
    leg = ROOT.TLegend(0.7, 0.7, 0.85, 0.85)
    leg.AddEntry(uubar, 'PYTHIA')
    leg.AddEntry(dssu, 'DSS')
    
    cu.cd(1)
    leg.Draw()
    cu.Print('uubar.png')
    
    cd = graphics.canvas1e('d+dbar')
    ddbar = f.Get('d_%(charge)s' % locals())
    ddbar.Add(f.Get('dbar_%(charge)s' % locals()))
    ddbar.SetTitle('d+dbar => #pi^{-}')
    dssd = tf1(dss, 0, 1, flavor='d+dbar', charge=charge)
    dostuff(cd, ddbar, dssd)
    
    cd.cd(1)
    leg.Draw()
    cd.Print('ddbar.png')
    
    cg = graphics.canvas1e('g')
    g = f.Get('g_%(charge)s' % locals())
    g.SetTitle('g => #pi^{-}')
    dssg = tf1(dss, 0, 1, flavor='g', charge=charge)
    dostuff(cg, g, dssg)
    
    cg.cd(1)
    leg.Draw()
    cg.Print('g.png')
    
    raw_input('wait here')

