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
from .runlists import transverse_run6 as transverse_runlist
from .util import getRun, fillList, hadd_interactive
from . import graphics

histDir         = '/Users/kocolosk/data/run6/spin2008/hist'
transHistDir    = '/Users/kocolosk/data/run6/spin2008/hist-transverse'
mcasymFile      = '/Users/kocolosk/data/run6/spin2008/mcasym.cphist.root'
datamcFile      = '/Users/kocolosk/data/run6/spin2008/merged.cphist.root'

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
    
    line = ROOT.TLine(0.075, 0.0, 1.0, 0.0)
    line.SetLineStyle(2)
    
    prelim = ROOT.TText()
    prelim.SetTextColor(44)
    prelim.SetTextAlign(21)
    
    latex = ROOT.TLatex()
    latex.SetTextSize(0.25)
    # latex.SetTextAlign(21)
    
    c.cd(1)
    hm = asym_m.GetAsymmetry('ll')
    hm.SetMarkerStyle(20)
    hm.SetTitle('')
    hm.Draw('e1')
    line.Draw()
    prelim.DrawText(0.5,0.065,"2006 STAR Preliminary")
    latex.DrawLatex(0.15, -0.085, '#pi^{-}')
    
    c.cd(2)
    hp = asym_p.GetAsymmetry('ll')
    hp.SetTitle('')
    hp.SetMarkerStyle(21)
    hp.Draw('e1')
    line.Draw()
    prelim.DrawText(0.5,0.065,"2006 STAR Preliminary")
    latex.DrawLatex(0.15, -0.085, '#pi^{+}')
    
    for h in (hm,hp):
        h.GetXaxis().SetTitle('p_{T}(#pi)/p_{T}(jet)')
        h.GetYaxis().SetRangeUser(-0.1, 0.1)
    
    latex = ROOT.TLatex()
    latex.DrawLatex(0.03,0.125,
    " #vec{p} + #vec{p} #rightarrow jet + #pi + X at #sqrt{s}=200 GeV \
             -1< #eta^{#pi}< 1 ")
    
    raw_input('wait here:')
    
    print '\nπ- A_{LL}'
    for bin in range(1, hm.GetNbinsX()+1):
        print '[%.2f-%.2f]  % .3f ± %.3f' % (hm.GetBinLowEdge(bin), 
            hm.GetBinLowEdge(bin+1), hm.GetBinContent(bin), hm.GetBinError(bin))
    
    print '\nπ+ A_{LL}'
    for bin in range(1, hp.GetNbinsX()+1):
        print '[%.2f-%.2f]  % .3f ± %.3f' % (hp.GetBinLowEdge(bin), 
            hp.GetBinLowEdge(bin+1), hp.GetBinContent(bin), hp.GetBinError(bin))


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
    
    for val in pid_asym_m:
        print val 
    for val in pid_asym_p:
        print val
    
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
        h.GetYaxis().SetRangeUser(-0.2, 0.2)
    
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
    for fname in allFiles[:]:
        run = getRun(fname)
        if run in transverse_runlist:
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


def trigger_bias():
    """
    comparison of MC asymmetries for minbias and 137222
    """
    f = ROOT.TFile(mcasymFile)
    # keys = ['STD','MAX','MIN','ZERO','GS_NLOC']
    keys = ['STD']
    mgr = HistogramManager(f, keys=keys + ['denom'])
    
    line = ROOT.TLine(0.1, 0.0, 0.9, 0.0)
    line.SetLineStyle(2)
    ROOT.gStyle.SetErrorX()
    
    color = {
        'STD': ROOT.kBlack,
        'MAX': ROOT.kRed,
        'MIN': ROOT.kGreen,
        'ZERO': ROOT.kBlue,
        'GS_NLOC': ROOT.kMagenta
    }
    
    smooth_factor = 2
    
    cmb = graphics.canvas2()
    cjp = graphics.canvas2()
    keepme = []
    diffs = {}
    for i,key in enumerate(keys):
        mb_m = mgr.anyspin['117001'].tracks_minus[key].Clone()
        mb_p = mgr.anyspin['117001'].tracks_plus[key].Clone()
        jp_m = mgr.anyspin['jetpatch'].tracks_minus[key].Clone()
        jp_p = mgr.anyspin['jetpatch'].tracks_plus[key].Clone()
        
        opt = i>0 and 'e2 same' or 'e2'
        
        for h in (mb_m, mb_p, jp_m, jp_p):
            h.Smooth(smooth_factor)
            h.GetXaxis().SetTitle('p_{T}(#pi)/p_{T}(jet)')
            h.GetXaxis().SetRangeUser(0.1, 0.9)
            h.GetYaxis().SetRangeUser(-0.08,0.08)
            h.SetLineColor(color[key])
            h.SetFillColor(color[key])
          
        for h in (mb_m, mb_p):
            h.SetLineColor(ROOT.kRed)
            h.SetFillColor(ROOT.kRed)
        
        cmb.cd(1)
        line.Draw()
        mb_m.SetTitle('MB MC Asymmetries for #pi^{-}')
        mb_m.Draw(opt)
        
        cmb.cd(2)
        line.Draw()
        mb_p.SetTitle('MB MC Asymmetries for #pi^{+}')
        mb_p.Draw(opt)
        
        cjp.cd(1)
        line.Draw()
        jp_m.SetTitle('JP1 MC Asymmetries for #pi^{-}')
        jp_m.Draw(opt)
        
        cjp.cd(2)
        line.Draw()
        jp_p.SetTitle('JP1 MC Asymmetries for #pi^{+}')
        jp_p.Draw(opt)
        
        diff_m = jp_m.Clone()
        diff_m.Add(mb_m, -1)
        
        diff_p = jp_p.Clone()
        diff_p.Add(mb_p, -1)
        
        for h in (diff_m, diff_p):
            h.GetYaxis().SetRangeUser(-0.03, 0.03)
            h.SetMarkerStyle(20)
        
        cdiff = graphics.canvas2(key)
        
        cdiff.cd(1)
        line.Draw()
        diff_m.SetTitle('#pi^{-} JP1 - MC for %(key)s' % locals())
        diff_m.Draw('e1')
        
        cdiff.cd(2)
        line.Draw()
        diff_p.SetTitle('#pi^{+} JP1 - MC for %(key)s' % locals())
        diff_p.Draw('e1')
        
        diffs[key] = (diff_m, diff_p)
        
        keepme.extend([jp_m, jp_p])
    
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
    
    raw_input('wait here:')


def subprocess_shift():
    import histos
    histos.simu = True
    f = ROOT.TFile(datamcFile)
    mgr = HistogramManager(f, ['hardP'])
    mb = {
        'all': mgr.anyspin['117001']['hardP'],
        'gg': mgr.gg['117001']['hardP'],
        'qg': mgr.qg['117001']['hardP'],
        'qq': mgr.qq['117001']['hardP']
    }
    jp = {
        'all': mgr.anyspin['jetpatch']['hardP'],
        'gg': mgr.gg['jetpatch']['hardP'],
        'qg': mgr.qg['jetpatch']['hardP'],
        'qq': mgr.qq['jetpatch']['hardP']
    }
    print mb
    print jp
    
    for key in ('gg','qg','qq'):
        mb[key].Divide(mb['all'])
        jp[key].Divide(jp['all'])
    
    for h in (mb,jp):
        h['gg'].SetLineColor(ROOT.kRed)
        h['qg'].SetLineColor(ROOT.kBlue)
        h['qq'].SetLineColor(ROOT.kGreen)
    
    c = graphics.canvas2()
    c.cd(1)
    mb['gg'].Draw()
    mb['gg'].SetTitle('MB')
    mb['gg'].GetYaxis().SetRangeUser(0., 0.8)
    mb['gg'].GetXaxis().SetRangeUser(0, 30)
    mb['qg'].Draw('same')
    mb['qq'].Draw('same')
    
    c.cd(2)
    jp['gg'].Draw()
    jp['gg'].SetTitle('JP')
    jp['gg'].SetXTitle('hardP')
    jp['gg'].GetYaxis().SetRangeUser(0., 0.8)
    jp['gg'].GetXaxis().SetRangeUser(0, 30)
    jp['qg'].Draw('same')
    jp['qq'].Draw('same')
    
    raw_input('wait here:')
    


def subprocess_fraction_z():
    import histos
    histos.simu = True
    ROOT.gStyle.SetErrorX()
    f = ROOT.TFile(datamcFile)
    charge = -1
    mgr = HistogramManager(f, ['z_away2'])
    mb = {
        'all': mgr.anyspin['117001'].trackHistograms(charge)['z_away2'],
        'gg': mgr.gg['117001'].trackHistograms(charge)['z_away2'],
        'qg': mgr.qg['117001'].trackHistograms(charge)['z_away2'],
        'qq': mgr.qq['117001'].trackHistograms(charge)['z_away2']
    }
    jp = {
        'all': mgr.anyspin['jetpatch'].trackHistograms(charge)['z_away2'],
        'gg': mgr.gg['jetpatch'].trackHistograms(charge)['z_away2'],
        'qg': mgr.qg['jetpatch'].trackHistograms(charge)['z_away2'],
        'qq': mgr.qq['jetpatch'].trackHistograms(charge)['z_away2']
    }
    print mb
    print jp
    
    for key in ('gg','qg','qq'):
        mb[key].Divide(mb['all'])
        jp[key].Divide(jp['all'])
    
    for h in (mb,jp):
        h['gg'].SetLineColor(ROOT.kRed)
        h['qg'].SetLineColor(ROOT.kBlue)
        h['qq'].SetLineColor(ROOT.kGreen)
    
    c = graphics.canvas2()
    c.cd(1)
    mb['gg'].Draw()
    mb['gg'].SetTitle('MB')
    mb['gg'].GetYaxis().SetRangeUser(0., 0.8)
    mb['gg'].GetXaxis().SetRangeUser(0, 30)
    mb['qg'].Draw('same')
    mb['qq'].Draw('same')
    
    c.cd(2)
    jp['gg'].Draw()
    jp['gg'].SetTitle('JP')
    jp['gg'].GetYaxis().SetRangeUser(0., 0.8)
    jp['gg'].GetXaxis().SetRangeUser(0, 30)
    jp['qg'].Draw('same')
    jp['qq'].Draw('same')
    
    raw_input('wait here:')


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


def meanpt():
    """plots mean pT in each z bin for for pions and tracks"""
    meanpt_p = hadd_interactive(histDir, runlist, 'jetpatch', 'anyspin', 
        'plus', 'meanpt')
    meanpt_m = hadd_interactive(histDir, runlist, 'jetpatch', 'anyspin', 
        'minus', 'meanpt')
    
    meanjetpt_p = hadd_interactive(histDir, runlist, 'jetpatch', 'anyspin', 
        'plus', 'meanjetpt')
    meanjetpt_m = hadd_interactive(histDir, runlist, 'jetpatch', 'anyspin', 
        'minus', 'meanjetpt')
    
    ## now get some results from simulations
    f = ROOT.TFile(datamcFile)
    mgr = HistogramManager(f, keys=['meanpt', 'meanjetpt'])
    meansim_p = mgr.anyspin['jetpatch'].tracks_plus['meanpt']
    meansim_m = mgr.anyspin['jetpatch'].tracks_minus['meanpt']
    meanjetsim_p = mgr.anyspin['jetpatch'].tracks_plus['meanjetpt']
    meanjetsim_m = mgr.anyspin['jetpatch'].tracks_minus['meanjetpt']
    
    mbsim_p = mgr.anyspin['117001'].tracks_plus['meanpt']
    mbsim_m = mgr.anyspin['117001'].tracks_minus['meanpt']
    mbjetsim_p = mgr.anyspin['117001'].tracks_plus['meanjetpt']
    mbjetsim_m = mgr.anyspin['117001'].tracks_minus['meanjetpt']
    
    meanjetpt_m.SetTitle('JP data, #pi -')
    meanjetpt_p.SetTitle('JP data, #pi +')
    
    meanjetsim_m.SetTitle('Monte Carlo #pi -')
    meanjetsim_p.SetTitle('Monte Carlo #pi +')
    
    for h in (meanjetpt_m, meanjetpt_p, meanjetsim_m, meanjetsim_p):
        h.SetMarkerStyle(21)
        h.SetXTitle('z')
        h.SetYTitle('< pT >')
    
    [h.SetMarkerStyle(20) for h in (mbjetsim_m, mbjetsim_p)]
    [h.SetMarkerStyle(24) for h in (mbsim_m, mbsim_p)]
    
    for h in (meanpt_m, meanpt_p, meansim_m, meansim_p):
        h.SetMarkerStyle(25)
    
    for h in (meansim_m, meansim_p, meanjetsim_m, meanjetsim_p):
        h.SetMarkerColor(ROOT.kRed)
        h.SetLineColor(ROOT.kRed)
    
    leg = ROOT.TLegend(.6, .7, .8, .88)
    leg.AddEntry(meanjetpt_p, 'jet')
    leg.AddEntry(meanpt_p, '#pi')
    
    ## first compare JP data/MC
    c = graphics.canvas2('Data / Monte Carlo comparison for JP')
    c.cd(1)
    meanjetpt_m.Draw()
    meanjetsim_m.Draw('][ hist same')
    meanpt_m.Draw('same')
    meansim_m.Draw('][ hist same')
    
    c.cd(2)
    meanjetpt_p.Draw()
    meanjetsim_p.Draw('][ hist same')
    meanpt_p.Draw('same')
    meansim_p.Draw('][ hist same')
    leg.Draw()
    
    ## now compare JP and MB simulations
    c2 = graphics.canvas2('Comparison of MB and JP simulations')
    
    c2.cd(1)
    meanjetsim_m.Draw('hist p')
    mbjetsim_m.Draw('hist p same')
    mbsim_m.Draw('hist p same')
    meansim_m.Draw('hist p same')
    
    c2.cd(2)
    meanjetsim_p.Draw('hist p')
    mbjetsim_p.Draw('hist p same')
    mbsim_p.Draw('hist p same')
    meansim_p.Draw('hist p same')
    
    leg2 = ROOT.TLegend(.68, .7, .88, .88)
    leg2.AddEntry(meanjetsim_p, 'JP')
    leg2.AddEntry(mbjetsim_p, 'MB')
    leg2.Draw()
    
    raw_input('wait here:')
    
    c.Print('jp-means.png')
    c2.Print('simu-means.png')

