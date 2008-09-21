# i need a better header template :-/
# this module is supposed to collect all the plots i'll be showing at SPIN 2008

import os
from glob import glob

import ROOT

from .asym import AsymmetryGenerator, ScalarCounts, Polarizations
from .histos import HistogramManager
from .runlists import long2_run6 as runlist
from .util import getRun
from . import graphics

histDir = '/Users/kocolosk/data/run6/hist'
zbins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]

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
    

def ssa():
    asym_zp = AsymmetryGenerator('asym_zp', bins=zbins, key='z_away2')
    asym_zm = AsymmetryGenerator('asym_zm', bins=zbins, key='z_away2')
    
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
            
            asym_zp.FillFromHistogramManager(mgr, 'jetpatch', 1, uu,ud,du,dd, \
                pol.py,pol.pb)
            asym_zm.FillFromHistogramManager(mgr, 'jetpatch', -1, uu,ud,du,dd, \
                pol.py,pol.pb)
            
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
        c = graphics.canvas2()
        canvases.append(c)
        
        hp = asym_zp.GetAsymmetry(key)
        hp.SetTitle(title[key] + ' SSA for #pi^{+}')
        hp.SetMarkerStyle(21)
    
        hm = asym_zm.GetAsymmetry(key)
        hm.SetTitle(title[key] + ' SSA for #pi^{-}')
        hm.SetMarkerStyle(20)
    
        for h in (hp,hm):
            h.SetMarkerColor(marker_color[key])
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
        
        c.Print('%s_ssa_%s.png' % (key, 'z'))
    
    raw_input('wait here:')

