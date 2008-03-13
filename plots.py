import os
import math
from glob import glob

import ROOT
import analysis

run5_tree_dir = '/Users/kocolosk/data/run5/tree'
run5_hist_dir = '/Users/kocolosk/data/run5/hist'

run6_tree_dir = '/Users/kocolosk/data/run6/tree'
run6_hist_dir = '/Users/kocolosk/data/run6/hist'

run5_scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
run6_scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run6.txt'

def spin2006_asymmetries():
    """asymmetries for charged pion production shown at SPIN 2006"""
    asym_plus = analysis.AsymmetryGenerator('asym_plus')
    asym_minus = analysis.AsymmetryGenerator('asym_minus')
    
    runlist = analysis.asym.golden_runlist_c
    
    scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
    scalars = analysis.ScalarCounts(scalar_path)
    
    polarizations = analysis.Polarizations.Online
    
    theory = analysis.asym.theoryCurves()
    plusGraphs = [ theory.getGraph('plus',key) for key in ('std','zero','max','min') ]
    minusGraphs = [ theory.getGraph('minus',key) for key in ('std','zero','max','min') ]
    
    ## systematic uncertainties
    baseline = -0.1
    syst_x = [3.0, 5.0, 7.0, 9.0]
    syst = {'plus': [7.3, 8.4, 7.5, 5.1], 'minus': [5.7, 6.0, 5.7, 7.1] }
    systGraph = {'plus': ROOT.TGraph(len(syst_x)+3), 'minus': ROOT.TGraph(len(syst_x)+3) }
    for charge in ('plus','minus'):
        systGraph[charge].SetPoint(0, 3.0, baseline)
        systGraph[charge].SetPoint(5, 9.0, baseline)
        systGraph[charge].SetPoint(6, 3.0, baseline)
        for i,val in enumerate(syst[charge]):
            systGraph[charge].SetPoint(i+1, syst_x[i], baseline + (val/1000.0))
    
    ## generate the asymmetries
    allFiles = glob(run5_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            print fname, run
            tfile = ROOT.TFile(fname)
            mgr = analysis.HistogramManager(tfile,['pt'])
            
            try:
                bin7 = scalars[str(run) + '-5-7']
                bin8 = scalars[str(run) + '-5-8']
                bin9 = scalars[str(run) + '-5-9']
            except KeyError:
                print run, 'is not in the scalars database'
                continue
            uu = bin7.uu + bin8.uu + bin9.uu
            ud = bin7.ud + bin8.ud + bin9.ud
            du = bin7.du + bin8.du + bin9.du
            dd = bin7.dd + bin8.dd + bin9.dd
            
            try:
                pol = polarizations[bin7.fill]
            except KeyError:
                print fill, 'has no online polarization values'
            
            asym_plus.FillFromHistogramManager(mgr, 'jetpatch', 1, uu, ud, du, dd, pol.py, pol.pb)
            asym_minus.FillFromHistogramManager(mgr, 'jetpatch', -1, uu, ud, du, dd, pol.py, pol.pb)
            tfile.Close()
    
    ## fun with graphics
    h1 = asym_plus.GetAsymmetry('ll')
    g1 = ROOT.TGraphErrors(h1)
    h2 = asym_minus.GetAsymmetry('ll')
    g2 = ROOT.TGraphErrors(h2)
    
    ## ignore bin width errors
    for gr in (g1,g2):
        for point in range(gr.GetN()):
            gr.SetPointError(point, 0.0, gr.GetErrorY(point))
    
    line = ROOT.TLine(2.0, 0.0, 10.0, 0.0)
    line.SetLineStyle(2)
    
    latex = ROOT.TLatex()
    
    leg = ROOT.TLegend(0.13, 0.65, 0.35, 0.88)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(plusGraphs[0],' GRSV-STD', 'l')
    leg.AddEntry(plusGraphs[1],' #Delta G =  0', 'l')
    leg.AddEntry(plusGraphs[2],' #Delta G =  G', 'l')
    leg.AddEntry(plusGraphs[3],' #Delta G = -G', 'l')
    
    bg = ROOT.TH1D(h1)
    bg.Reset()
    bg.SetYTitle(' A_{LL}')
    bg.GetYaxis().SetRangeUser(-0.11, 0.11)
    
    ## pi-plus
    c1 = ROOT.TCanvas('c1','A_{LL} for #pi^{+}')
    bg.SetXTitle('#pi^{+} P_{T} (GeV/c)')
    bg.DrawCopy()
    g1.SetMarkerSize(0.9);
    g1.SetMarkerStyle(21)
    g1.Draw('p')
    [ g.Draw('l') for g in plusGraphs ]
    systGraph['plus'].SetLineColor(1)
    systGraph['plus'].SetFillColor(15)
    systGraph['plus'].Draw('fl')
    line.Draw('same')
    leg.Draw('p')
    latex.DrawLatex(2.3,0.12," #vec{p} + #vec{p} #rightarrow #pi^{+} + X at #sqrt{s}=200 GeV \
                        -1< #eta^{#pi}< 1 ")
    latex.DrawLatex(2.6,-0.07,"2005 STAR Preliminary");
    
    ## pi-minus
    c2 = ROOT.TCanvas('c2','A_{LL} for #pi^{-}')
    bg.SetXTitle('#pi^{-} P_{T} (GeV/c)')
    bg.DrawCopy()
    g2.SetMarkerSize(0.9);
    g2.SetMarkerStyle(20)
    g2.Draw('p')
    [ g.Draw('l') for g in minusGraphs ]
    systGraph['minus'].SetLineColor(1)
    systGraph['minus'].SetFillColor(15)
    systGraph['minus'].Draw('fl')
    line.Draw('same')
    leg.Draw('p')
    latex.DrawLatex(2.3,0.12," #vec{p} + #vec{p} #rightarrow #pi^{-} + X at #sqrt{s}=200 GeV \
                        -1< #eta^{#pi}< 1 ")
    latex.DrawLatex(2.6,-0.07,"2005 STAR Preliminary");
    
    raw_input('wait here:')


def asymmetries_for_publication_run5(runlist=None):
    """final results for inclusive asymmetries from Run 5"""
    asym_plus = analysis.AsymmetryGenerator('asym_plus', key='pt')
    asym_minus = analysis.AsymmetryGenerator('asym_minus', key='pt')
    
    scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
    scalars = analysis.ScalarCounts(scalar_path)
    
    polarizations = analysis.Polarizations.Final
    
    #theory = analysis.asym.theoryCurves()
    #plusGraphs = [ theory.getGraph('plus',key) for key in ('std','zero','max','min') ]
    #minusGraphs = [ theory.getGraph('minus',key) for key in ('std','zero','max','min') ]
    from analysis.asym import theoryCurves
    plusGraphs = [
    theoryCurves(analysis.asym.werner_plus_dss_cteqm5_std, analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
    theoryCurves(analysis.asym.werner_plus_dss_cteqm5_zero, analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
    theoryCurves(analysis.asym.werner_plus_dss_cteqm5_max, analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),
    theoryCurves(analysis.asym.werner_plus_dss_cteqm5_min, analysis.xsec.werner_plus_dss_cteqm5_pt).getGraph(),   
    ]
    minusGraphs = [
    theoryCurves(analysis.asym.werner_minus_dss_cteqm5_std, analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
    theoryCurves(analysis.asym.werner_minus_dss_cteqm5_zero, analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
    theoryCurves(analysis.asym.werner_minus_dss_cteqm5_max, analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),
    theoryCurves(analysis.asym.werner_minus_dss_cteqm5_min, analysis.xsec.werner_minus_dss_cteqm5_pt).getGraph(),    
    ]
    
    for grList in (plusGraphs, minusGraphs):
        grList[1].SetLineStyle(3)
        grList[1].SetLineColor(ROOT.kBlue)
        grList[2].SetLineStyle(4)
        grList[2].SetLineColor(ROOT.kRed)
        grList[3].SetLineStyle(2)
        grList[3].SetLineColor(ROOT.kGreen)
        for gr in grList:
            gr.SetLineWidth(3)
    
    ## systematic uncertainties
    baseline = -0.1
    syst_x = [3.0, 5.0, 7.0, 9.0]
    syst = {'plus': [7.3, 8.4, 7.5, 5.1], 'minus': [5.7, 6.0, 5.7, 7.1] }
    systGraph = {'plus': ROOT.TGraph(len(syst_x)+3), 'minus': ROOT.TGraph(len(syst_x)+3) }
    for charge in ('plus','minus'):
        systGraph[charge].SetPoint(0, 3.0, baseline)
        systGraph[charge].SetPoint(5, 9.0, baseline)
        systGraph[charge].SetPoint(6, 3.0, baseline)
        for i,val in enumerate(syst[charge]):
            systGraph[charge].SetPoint(i+1, syst_x[i], baseline + (val/1000.0))
    
    ## generate the asymmetries
    allFiles = glob(run5_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            print fname, run
            tfile = ROOT.TFile(fname)
            mgr = analysis.HistogramManager(tfile,['pt', 'pt_near', 'pt_away'])
            
            try:
                bin7 = scalars[str(run) + '-5-7']
                bin8 = scalars[str(run) + '-5-8']
                bin9 = scalars[str(run) + '-5-9']
            except KeyError:
                print run, 'is not in the scalars database'
                continue
            uu = bin7.uu + bin8.uu + bin9.uu
            ud = bin7.ud + bin8.ud + bin9.ud
            du = bin7.du + bin8.du + bin9.du
            dd = bin7.dd + bin8.dd + bin9.dd
            
            try:
                pol = polarizations[bin7.fill]
            except KeyError:
                print bin7.fill, 'has no final polarization values'
                continue
            
            asym_plus.FillFromHistogramManager(mgr, 'alltrigs', 1, uu, ud, du, dd, pol.py, pol.pb)
            asym_minus.FillFromHistogramManager(mgr, 'alltrigs', -1, uu, ud, du, dd, pol.py, pol.pb)
            #asym_plus.FillFromHistogramManager(mgr, 'jetpatch', 1, uu, ud, du, dd, pol.py, pol.pb)
            #asym_minus.FillFromHistogramManager(mgr, 'jetpatch', -1, uu, ud, du, dd, pol.py, pol.pb)
            tfile.Close()
            
    ## fun with graphics
    h1 = asym_plus.GetAsymmetry('ll')
    g1 = ROOT.TGraphErrors(h1)
    h2 = asym_minus.GetAsymmetry('ll')
    g2 = ROOT.TGraphErrors(h2)
    
    ## ignore bin width errors
    for gr in (g1,g2):
        for point in range(gr.GetN()):
            gr.SetPointError(point, 0.0, gr.GetErrorY(point))
            
    line = ROOT.TLine(2.0, 0.0, 10.0, 0.0)
    line.SetLineStyle(2)
    
    latex = ROOT.TLatex()
    
    leg = ROOT.TLegend(0.13, 0.67, 0.40, 0.89)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(plusGraphs[0],' GRSV-std', 'l')
    leg.AddEntry(plusGraphs[2],' GRSV #Delta g =     g input', 'l')
    leg.AddEntry(plusGraphs[1],' GRSV #Delta g =     0 input', 'l')
    leg.AddEntry(plusGraphs[3],' GRSV #Delta g = -g input', 'l')
    
    bg = ROOT.TH1D(h1)
    bg.Reset()
    bg.SetYTitle(' A_{LL}')
    bg.GetYaxis().SetRangeUser(-0.11, 0.11)
    
    ## pi-plus
    c1 = ROOT.TCanvas('c1','A_{LL} for #pi^{+}')
    bg.SetXTitle('#pi^{+} P_{T} [GeV/c]')
    bg.DrawCopy()
    g1.SetMarkerSize(0.9);
    g1.SetMarkerStyle(21)
    g1.Draw('p')
    [ g.Draw('l') for g in plusGraphs ]
    systGraph['plus'].SetLineColor(1)
    systGraph['plus'].SetFillColor(12)
    systGraph['plus'].Draw('fl')
    line.Draw('same')
    leg.Draw('p')
    latex.DrawLatex(2.2,0.12,"STAR #vec{p} + #vec{p} #rightarrow #pi^{+} + X at #sqrt{s}=200 GeV \
                -1< #eta^{#pi}< 1 ")
    
    ## pi-minus
    c2 = ROOT.TCanvas('c2','A_{LL} for #pi^{-}')
    bg.SetXTitle('#pi^{-} P_{T} [GeV/c]')
    bg.DrawCopy()
    g2.SetMarkerSize(0.9);
    g2.SetMarkerStyle(20)
    g2.Draw('p')
    [ g.Draw('l') for g in minusGraphs ]
    systGraph['minus'].SetLineColor(1)
    systGraph['minus'].SetFillColor(12)
    systGraph['minus'].Draw('fl')
    line.Draw('same')
    leg.Draw('p')
    latex.DrawLatex(2.2,0.12,"STAR #vec{p} + #vec{p} #rightarrow #pi^{-} + X at #sqrt{s}=200 GeV \
                -1< #eta^{#pi}< 1 ")
    
    
    ## combo plot
    c3 = ROOT.TCanvas('c3', 'A_{LL} combined', 0, 0, 800, 400)
    c3.Draw()
    titlepad = ROOT.TPad('titlepad', '', 0.0, 0.9, 1.0, 1.0)
    
    leftpad  = ROOT.TPad('leftpad','', 0.0, 0.0, 0.5, 1.0)
    leftpad.SetLeftMargin(0.14)
    leftpad.SetRightMargin(0.02)
    
    rightpad = ROOT.TPad('rightpad','', 0.5, 0.0, 1.0, 1.0)
    rightpad.SetLeftMargin(0.11)
    rightpad.SetRightMargin(0.05)
    
    for pad in (titlepad, leftpad, rightpad):
        pad.Draw()
        pad.SetFillColor(10)
        pad.SetBorderMode(0)
        pad.SetFillStyle(4000) ## make it transparent
        
    leg2 = ROOT.TLegend(0.16, 0.67, 0.54, 0.90)
    leg2.SetFillStyle(0)
    leg2.SetBorderSize(0)
    leg2.AddEntry(plusGraphs[0],' GRSV-std', 'l')
    leg2.AddEntry(plusGraphs[2],' GRSV #Deltag =  g input', 'l')
    leg2.AddEntry(plusGraphs[1],' GRSV #Deltag =  0 input', 'l')
    leg2.AddEntry(plusGraphs[3],' GRSV #Deltag = -g input', 'l')
    
    titlepad.cd()
    latex.SetTextSize(0.7)
    latex.SetTextAlign(21)
    latex.DrawLatex(0.5,0.2,"STAR #vec{p} + #vec{p} #rightarrow #pi + X at #sqrt{s}=200 GeV \
         |#eta^{#pi}|<1.0")
    
    leftpad.cd()
    bg.SetXTitle('')
    bg.SetYTitle('A_{LL}')
    bg.GetYaxis().SetTitleSize(0.05)
    bg.GetYaxis().SetTitleOffset(1.22)
    bg.DrawCopy()
    g2.SetMarkerSize(0.9);
    g2.SetMarkerStyle(20)
    g2.Draw('p')
    [ g.Draw('l') for g in minusGraphs ]
    systGraph['minus'].SetLineColor(1)
    systGraph['minus'].SetFillColor(12)
    systGraph['minus'].Draw('fl')
    line.Draw('same')
    leg2.Draw('p')
    latex.SetTextSize(0.2)
    latex.SetTextAlign(21)
    latex.DrawLatex(4.0,-0.075,'#pi^{-}')
    
    rightpad.cd()
    bg.SetXTitle('#pi P_{T} [GeV/c]')
    bg.SetYTitle('')
    bg.DrawCopy()
    g1.SetMarkerSize(0.9);
    g1.SetMarkerStyle(21)
    g1.Draw('p')
    [ g.Draw('l') for g in plusGraphs ]
    systGraph['plus'].SetLineColor(1)
    systGraph['plus'].SetFillColor(12)
    systGraph['plus'].Draw('fl')
    line.Draw('same')
    latex.DrawLatex(4.0,-0.075,'#pi^{+}')
    
    'pi-plus fit to pol0'
    h1.Fit('pol0', 'N', '', 2.0, 10.0)
    'pi-minus fit to pol0'
    h2.Fit('pol0', 'N', '', 2.0, 10.0)
    
    raw_input('wait here:')


def jet_correlations_run5():
    """3-D deta-dphi plot -- possible paper plot at one time. \
    Also plots the uncorrected pion momentum fraction for near-side
    and away-side
    """
    style = ROOT.TStyle(ROOT.gStyle)
    style.SetOptStat(0)   
    style.SetLabelOffset(-0.01,'xy')
    style.SetLabelSize(0.035,'xy')
    style.SetTitleOffset(1.2,'y')
    style.cd()
    
    #fig3 = ROOT.TH2D('figure3','',50,-1.5*math.pi,0.5*math.pi,50,-2.0,0.4)
    #fig3.SetYTitle('#eta pion - #eta jet')
    #fig3.SetXTitle('#phi pion - #phi jet')
    
    #fig3b = ROOT.TH2D('figure3b','',50,-math.pi,math.pi,50,-1.5,1.5)
    
    #xbins = [2.0, 3.0, 4.0, 5.5, 7.0, 10.0]
    #ar = array('d',xbins)
    #fig3c = ROOT.TH2D('figure3c','Uncorrected Pion Momentum Fraction',len(xbins)-1,ar,50,0.,1.)
    #fig3c_away = ROOT.TH2D('figure3c_away','',len(xbins)-1,ar,50,0.,1.)
    
    runlist = None
    
    h1 = None
    h2 = None
    h3 = None
    
    ## silly hack
    keepMeOpen = []
    
    allFiles = glob(run5_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            print fname, run
            tfile = ROOT.TFile(fname)
            mgr = analysis.HistogramManager(tfile,['dphi_deta', 'z', 'z_away'])
            
            if h1 is None:
                h1 = mgr['anyspin']['jetpatch'].tracks_sum['dphi_deta'].Clone()
                h2 = mgr['anyspin']['jetpatch'].tracks_sum['z'].Clone()
                h3 = mgr['anyspin']['jetpatch'].tracks_sum['z_away'].Clone()
                keepMeOpen.append(tfile)
            else:
                h1.Add(mgr['anyspin']['jetpatch'].tracks_sum['dphi_deta'])
                h2.Add(mgr['anyspin']['jetpatch'].tracks_sum['z'])
                h3.Add(mgr['anyspin']['jetpatch'].tracks_sum['z_away'])
    
    c1 = ROOT.TCanvas('c1')
    c1.SetLogz()
    h1.SetXTitle('#phi pion - #phi jet')
    h1.SetYTitle('#eta pion - #eta jet')
    h1.GetYaxis().SetRangeUser(-2.0, 0.4)
    h1.DrawCopy('lego2')
    
    #reset some styles
    style.SetLabelOffset(0.005,'xy')
    style.SetLabelSize(0.04,'xy')
    style.SetTitleOffset(1,'y')
    
    c2 = ROOT.TCanvas('c2')
    h2.FitSlicesY()
    fig3c_mean = ROOT.gDirectory.Get('%s_1' % (h2.GetName(),))
    fig3c_mean.SetTitle('Uncorrected pion momentum fraction')
    fig3c_mean.SetXTitle('#pi p_{T} [GeV/c]')
    fig3c_mean.SetYTitle('< p_{T,#pi} / p_{T,jet} >')
    fig3c_mean.SetAxisRange(0,1,'y')
    fig3c_mean.SetMarkerStyle(21)
    
    h3.FitSlicesY()
    fig3c_away_mean = ROOT.gDirectory.Get('%s_1' % (h3.GetName(),))
    fig3c_away_mean.SetMarkerStyle(25)
    fig3c_away_mean.SetMarkerColor(ROOT.kRed)
    fig3c_away_mean.SetLineColor(ROOT.kRed)
    
    leg = ROOT.TLegend(0.75,0.8,0.89,0.89)
    leg.AddEntry(fig3c_mean,'dR < 0.4','p')
    leg.AddEntry(fig3c_away_mean,'dR > 1.5','p')
    
    fig3c_mean.DrawCopy()
    fig3c_away_mean.Draw('same')
    leg.Draw('same')
    
    c3 = ROOT.TCanvas('c3')
    h3.Draw()
    
    c4 = ROOT.TCanvas('combo plot')
    mainpad = ROOT.TPad('mainpad', '', 0.0, 0.0, 1.0, 1.0)
    mainpad.SetLeftMargin(0.1)
    mainpad.SetRightMargin(0.03)
    mainpad.SetTopMargin(0.05)
    mainpad.SetBottomMargin(0.1)
    #mainpad.SetTicky()
    
    insetpad = ROOT.TPad('insetpad', '', 0.52, 0.13, 0.97, 0.6)
    insetpad.SetLogz()
    #insetpad.SetLeftMargin(0.)
    #insetpad.SetRightMargin(0.)
    #insetpad.SetTopMargin(0.)
    #insetpad.SetBottomMargin(0.)
    
    for pad in (mainpad, insetpad):
        pad.Draw()
        pad.SetFillColor(10)
        pad.SetBorderMode(0)
        pad.SetFillStyle(4000) ## make it transparent
    
    mainpad.cd()
    fig3c_mean.SetTitle('')
    fig3c_mean.GetYaxis().SetRangeUser(0., 0.8)
    fig3c_mean.Draw()
    fig3c_away_mean.Draw('same')
    
    insetpad.cd()
    h1.GetXaxis().SetTitle('#Delta #phi')
    h1.GetXaxis().SetTitleSize(0.1)
    h1.GetXaxis().SetLabelSize(0.06)
    h1.GetYaxis().SetTitle('#Delta #eta')
    h1.GetYaxis().SetTitleSize(0.1)
    h1.GetYaxis().SetLabelSize(0.06)
    
    ## this is temporary till we get the cuts right
    h1.GetZaxis().SetRangeUser(10, 600000)
    #h1.GetYaxis().SetRangeUser(-1.5, 0.4)
    
    h1.DrawCopy('lego2')
    
    mainpad.cd()
    leg2 = ROOT.TLegend(0.15, 0.78, 0.45, 0.92)
    leg2.AddEntry(fig3c_mean,'< z > trigger jet','p')
    leg2.AddEntry(fig3c_away_mean,'< z > away-side jet','p')
    leg2.Draw()
    
    raw_input('wait here:')


def jet_correlations_run6():
    """3-D deta-dphi plot -- possible paper plot at one time. \
    Also plots the uncorrected pion momentum fraction for near-side
    and away-side
    """
    style = ROOT.TStyle(ROOT.gStyle)
    style.SetOptStat(0)   
    style.SetLabelOffset(-0.01,'xy')
    style.SetLabelSize(0.035,'xy')
    style.SetTitleOffset(1.2,'y')
    style.cd()
    
    runlist = None
    
    h1 = None
    h2 = None
    h3 = None
    
    ## silly hack
    keepMeOpen = []
    
    allFiles = glob(run6_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            print fname, run
            tfile = ROOT.TFile(fname)
            mgr = analysis.HistogramManager(tfile,['dphi_deta', 'z', 'z_away'])
            
            if h1 is None:
                h1 = mgr['anyspin']['alltrigs'].tracks_sum['dphi_deta'].Clone()
                h2 = mgr['anyspin']['alltrigs'].tracks_sum['z'].Clone()
                h3 = mgr['anyspin']['alltrigs'].tracks_sum['z_away'].Clone()
                keepMeOpen.append(tfile)
            else:
                h1.Add(mgr['anyspin']['alltrigs'].tracks_sum['dphi_deta'])
                h2.Add(mgr['anyspin']['alltrigs'].tracks_sum['z'])
                h3.Add(mgr['anyspin']['alltrigs'].tracks_sum['z_away'])
    
    c1 = ROOT.TCanvas('c1')
    c1.SetLogz()
    h1.SetXTitle('#phi pion - #phi jet')
    h1.SetYTitle('#eta pion - #eta jet')
    h1.DrawCopy('lego2')
    
    #reset some styles
    style.SetLabelOffset(0.005,'xy')
    style.SetLabelSize(0.04,'xy')
    style.SetTitleOffset(1,'y')
    
    c2 = ROOT.TCanvas('c2')
    h2.FitSlicesY()
    fig3c_mean = ROOT.gDirectory.Get('%s_1' % (h2.GetName(),))
    fig3c_mean.SetTitle('Uncorrected pion momentum fraction')
    fig3c_mean.SetXTitle('#pi p_{T}')
    fig3c_mean.SetYTitle('< p_{T,#pi} / p_{T,jet} >')
    fig3c_mean.SetAxisRange(0,1,'y')
    fig3c_mean.SetMarkerStyle(21)
    
    h3.FitSlicesY()
    fig3c_away_mean = ROOT.gDirectory.Get('%s_1' % (h3.GetName(),))
    fig3c_away_mean.SetMarkerStyle(25)
    fig3c_away_mean.SetMarkerColor(ROOT.kRed)
    fig3c_away_mean.SetLineColor(ROOT.kRed)
    
    leg = ROOT.TLegend(0.75,0.8,0.89,0.89)
    leg.AddEntry(fig3c_mean,'dR < 0.4','p')
    leg.AddEntry(fig3c_away_mean,'dR > 1.5','p')
    
    fig3c_mean.Draw()
    fig3c_away_mean.Draw('same')
    leg.Draw('same')
    
    c3 = ROOT.TCanvas('c3')
    h3.Draw()
    
    raw_input('wait here:')


def asymmetry_statistics_comparison():
    """plots statistical prescision of A_{LL} for Run 5 (prelim + final)
    and Run 5 + Run 6 combined away-side measurement
    """
    prelim = analysis.AsymmetryGenerator('prelim')
    final = analysis.AsymmetryGenerator('final')
    combo = analysis.AsymmetryGenerator('combo', key='pt_away')
    
    scalars_run5 = analysis.ScalarCounts(run5_scalar_path)
    scalars_run6 = analysis.ScalarCounts(run6_scalar_path)
    
    polarizations = analysis.Polarizations.Final
    polarizations_prelim = analysis.Polarizations.Online
    
    ## generate the asymmetries
    allFiles = glob(run5_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        print fname, run
        tfile = ROOT.TFile(fname)
        mgr = analysis.HistogramManager(tfile,['pt', 'pt_away'])
            
        try:
            bin7 = scalars_run5[str(run) + '-5-7']
            bin8 = scalars_run5[str(run) + '-5-8']
            bin9 = scalars_run5[str(run) + '-5-9']
        except KeyError:
            print run, 'is not in the scalars database'
            continue
        uu = bin7.uu + bin8.uu + bin9.uu
        ud = bin7.ud + bin8.ud + bin9.ud
        du = bin7.du + bin8.du + bin9.du
        dd = bin7.dd + bin8.dd + bin9.dd
        
        if run in analysis.golden_runlist_c:
            try:
                pol = polarizations_prelim[bin7.fill]
                prelim.FillFromHistogramManager(mgr, 'jetpatch', 1, uu, ud, du, dd, pol.py, pol.pb)
            except KeyError:
                print bin7.fill, 'has no preliminary polarization values'
            
        try:
            pol = polarizations[bin7.fill]
            final.FillFromHistogramManager(mgr, 'jetpatch', 1, uu, ud, du, dd, pol.py, pol.pb)
            combo.FillFromHistogramManager(mgr, 'alltrigs', 1, uu, ud, du, dd, pol.py, pol.pb)
        except KeyError:
            print bin7.fill, 'has no final polarization values'
        
        tfile.Close()
    
    allFiles = glob(run6_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        print fname, run
        tfile = ROOT.TFile(fname)
        mgr = analysis.HistogramManager(tfile,['pt', 'pt_away'])
            
        try:
            bin6 = scalars_run6[str(run) + '-5-6']
            bin7 = scalars_run6[str(run) + '-5-7']
            bin8 = scalars_run6[str(run) + '-5-8']
            bin9 = scalars_run6[str(run) + '-5-9']
        except KeyError:
            print run, 'is not in the scalars database'
            continue
        uu = bin6.uu + bin7.uu + bin8.uu + bin9.uu
        ud = bin6.ud + bin7.ud + bin8.ud + bin9.ud
        du = bin6.du + bin7.du + bin8.du + bin9.du
        dd = bin6.dd + bin7.dd + bin8.dd + bin9.dd
        
        try:
            pol = polarizations[bin6.fill]
            combo.FillFromHistogramManager(mgr, 'alltrigs', 1, uu, ud, du, dd, pol.py, pol.pb)
        except KeyError:
            print bin6.fill, 'has no final polarization values'
        
        tfile.Close()
        
    h1 = prelim.GetAsymmetry('ll')
    h2 = final.GetAsymmetry('ll')
    h3 = combo.GetAsymmetry('ll')
    
    g1 = ROOT.TGraphErrors(h1)
    g2 = ROOT.TGraphErrors(h2)
    g3 = ROOT.TGraphErrors(h3)
    
    for point in range(g1.GetN()):
        x = h1.GetBinCenter(point+1)
        g1.SetPoint(point, x-0.2, 0)
        g2.SetPoint(point, x, 0)
        g3.SetPoint(point, x+0.2, 0)
        for gr in (g1,g2,g3):
            gr.SetMarkerStyle(21)
            gr.SetPointError(point, 0.0, gr.GetErrorY(point))
            
    g2.SetMarkerColor(ROOT.kRed)
    g2.SetLineColor(ROOT.kRed)
    g3.SetMarkerColor(ROOT.kGreen)
    g3.SetLineColor(ROOT.kGreen)
    
    leg = ROOT.TLegend(0.13, 0.67, 0.40, 0.89)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(g1, '2005 prelim', 'p')
    leg.AddEntry(g2, '2005 final', 'p')
    leg.AddEntry(g3, '2005/6 away-side', 'p')
    
    bg = ROOT.TH1D(h1)
    bg.Reset()
    bg.SetTitle('Statistical Precision of Various A_{LL} Measurements')
    bg.SetXTitle('p_{T}')
    bg.GetYaxis().SetRangeUser(-0.06, 0.06)
    
    c = ROOT.TCanvas()
    bg.DrawCopy()
    g1.Draw('p')
    g2.Draw('p')
    g3.Draw('p')
    leg.Draw()
    
    raw_input('press enter:')
    


def trigger_bias_using_away_side(runlist=None, trgname='alltrigs'):
    """plot asym vs. pt for near + away, fit with pol0 and compare"""
    generator = { 
    'near_plus'  : analysis.AsymmetryGenerator('near_plus', key='pt_near'), 
    'near_minus' : analysis.AsymmetryGenerator('near_minus', key='pt_near'),
    'away_plus'  : analysis.AsymmetryGenerator('away_plus', key='pt_away'), 
    'away_minus' : analysis.AsymmetryGenerator('away_minus', key='pt_away') 
    }
    
    scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
    scalars = analysis.ScalarCounts(scalar_path)
    
    polarizations = analysis.Polarizations.Final
    
    ## generate the asymmetries
    allFiles = glob(run5_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            print fname, run
            tfile = ROOT.TFile(fname)
            mgr = analysis.HistogramManager(tfile,['pt_near', 'pt_away'])
            
            try:
                bin7 = scalars[str(run) + '-5-7']
                bin8 = scalars[str(run) + '-5-8']
                bin9 = scalars[str(run) + '-5-9']
            except KeyError:
                print run, 'is not in the scalars database'
                continue
            uu = bin7.uu + bin8.uu + bin9.uu
            ud = bin7.ud + bin8.ud + bin9.ud
            du = bin7.du + bin8.du + bin9.du
            dd = bin7.dd + bin8.dd + bin9.dd
            
            try:
                pol = polarizations[bin7.fill]
            except KeyError:
                print bin7.fill, 'has no final polarization values'
                continue
            
            generator['near_plus'].FillFromHistogramManager(mgr, trgname, 1, uu, ud, du, dd, pol.py, pol.pb)
            generator['near_minus'].FillFromHistogramManager(mgr, trgname, -1, uu, ud, du, dd, pol.py, pol.pb)
            generator['away_plus'].FillFromHistogramManager(mgr, trgname, 1, uu, ud, du, dd, pol.py, pol.pb)
            generator['away_minus'].FillFromHistogramManager(mgr, trgname, -1, uu, ud, du, dd, pol.py, pol.pb)
            tfile.Close()
    
    #ROOT.gStyle.SetOptStat('n')
    #ROOT.gStyle.SetOptFit(111)
    fit = {
    'near_plus' : ROOT.TF1('near_plus', 'pol0'),
    'near_minus' : ROOT.TF1('near_minus', 'pol0'),
    'away_plus' : ROOT.TF1('away_plus', 'pol0'),
    'away_minus' : ROOT.TF1('away_minus', 'pol0'),
    }
    fit['away_plus'].SetLineColor(ROOT.kRed)
    fit['away_minus'].SetLineColor(ROOT.kRed)
    
    c1 = ROOT.TCanvas('plus','Comparison of near and away side for #pi^{+}')
    h1_near = generator['near_plus'].GetAsymmetry('ll')
    h1_near.GetYaxis().SetRangeUser(-0.11, 0.11)
    h1_away = generator['away_plus'].GetAsymmetry('ll')
    h1_away.SetLineColor(ROOT.kRed)
    h1_near.Draw()
    h1_near.Fit(fit['near_plus'],'','same')
    h1_away.Draw('same')
    h1_away.Fit(fit['away_plus'],'','same')
    
    leg1 = ROOT.TLegend(0.13,0.7,0.43,0.89)
    leg1.AddEntry(fit['near_plus'],'%f +/- %f' % 
        (fit['near_plus'].GetParameter(0), fit['near_plus'].GetParError(0)),'l')
    leg1.AddEntry(fit['away_plus'],'%f +/- %f' % 
        (fit['away_plus'].GetParameter(0), fit['away_plus'].GetParError(0)),'l')
    leg1.Draw()
    
    c2 = ROOT.TCanvas('minus','Comparison of near and away side for #pi^{-}')
    h2_near = generator['near_minus'].GetAsymmetry('ll')
    h2_near.GetYaxis().SetRangeUser(-0.11, 0.11)
    h2_away = generator['away_minus'].GetAsymmetry('ll')
    h2_away.SetLineColor(ROOT.kRed)
    h2_near.Draw()
    h2_near.Fit(fit['near_minus'],'','same')
    h2_away.Draw('same')
    h2_away.Fit(fit['away_minus'],'','same')
    
    leg2 = ROOT.TLegend(0.13,0.7,0.43,0.89)
    leg2.AddEntry(fit['near_minus'],'%f +/- %f' % 
        (fit['near_minus'].GetParameter(0), fit['near_minus'].GetParError(0)),'l')
    leg2.AddEntry(fit['away_minus'],'%f +/- %f' % 
        (fit['away_minus'].GetParameter(0), fit['away_minus'].GetParError(0)),'l')
    leg2.Draw()
    
    print 'Size of systematic assigned if we take the difference btw the fits with errors:'
    val = math.fabs( fit['near_plus'].GetParameter(0) - fit['away_plus'].GetParameter(0) )
    err = math.sqrt(fit['near_plus'].GetParError(0) ** 2 + fit['away_plus'].GetParError(0) ** 2)
    print 'plus  : %f' % (val+err,)
    val = math.fabs( fit['near_minus'].GetParameter(0) - fit['away_minus'].GetParameter(0) )
    err = math.sqrt(fit['near_minus'].GetParError(0) ** 2 + fit['away_minus'].GetParError(0) ** 2)
    print 'minus : %f' % (val+err,)
    
    
    raw_input('press enter:')


def qa_pid():
    allFiles = glob(run5_hist_dir + '/chargedPions_*.hist.root')
    fill_runlists = {}
    reverse_dict = {}
    
    for fname in allFiles:
        run = analysis.getRun(fname)
        reverse_dict[run] = 0.0
    answer = analysis.getAllFills(reverse_dict.keys())
    for run,fill in answer:
        reverse_dict[run] = int(fill)
        try:
            fill_runlists[fill].append(run)
        except KeyError:
            fill_runlists[fill] = [run]    
    
    nSigmaRun = ROOT.TH1D('nSigmaRun', 'Mean nSigmaPion per run', len(allFiles), 0.5, len(allFiles)+0.5)
    nSigmaFill = ROOT.TH1D('nSigmaFill', 'blerg2', len(fill_runlists), 0.5, len(fill_runlists)+0.5)
    
    ps = ROOT.TPostScript('blerg.ps')
    c = ROOT.TCanvas('c','',100,100,600,800)
    
    pad = 1
    ROOT.gStyle.SetOptStat('m')
    for row,fname in enumerate(allFiles):
        if row % 15 == 0:
            c.Update()
            ps.NewPage()
            c.Clear()
            c.Divide(3,5)
            pad = 1 
        c.cd(pad)    
        print fname
        tfile = ROOT.TFile(fname)
        mgr = analysis.HistogramManager(tfile,['nSigmaPion'])
        h = mgr.anyspin['alltrigs'].tracks_sum['nSigmaPion']
        run = analysis.getRun(fname)
        h.SetTitle('%3d - %d - %d' % (row+1, run, reverse_dict[run]))
        h.SetStats(True)
        if h.GetMean() > 0: 
            h.SetLineColor(ROOT.kRed)
            print 'bad = fill, run'
        h.DrawCopy()
        nSigmaRun.SetBinContent(row+1, h.GetMean())
        nSigmaRun.SetBinError(row+1, h.GetMeanError())
        pad += 1
    ps.Close()
    c = ROOT.TCanvas()
    nSigmaRun.Draw()
    raw_input('press enter:')


def jetpatch_phi_correlation(tree, patchNumber):
    patchPhi = analysis.histos.JetCuts.patchPhi2006
    h = ROOT.TH1D('h','',720,-360,360)
    for entry in tree:
        for i in range(12):
            adc = tree.event.jetPatchAdc(i)
            if adc > analysis.histos.JetCuts.triggerThresholds[137221]:
                for jet in tree.event.jets():
                    diff = math.degrees(jet.Phi()) - patchPhi[i]
                    h.Fill(diff)
    h.Draw()
    raw_input('press enter:')


def pid_calibration(h, useElectrons=True):
    """takes an nSigmaPion histo and does the triple-Gaussian fit"""
    h.SetStats(True)
    #h.SetTitle('n#sigma(#pi) calibration for F7305')
    h.SetXTitle('n#sigma(#pi)')
    ROOT.gStyle.SetOptFit(111)
    
    if useElectrons:
        fit = ROOT.TF1('fit','gaus(0)+gaus(3)+gaus(6)', -6.0, 6.0)
    else:
        fit = ROOT.TF1('fit','gaus(0)+gaus(3)', -6.0, 6.0)
    fit.SetParameter(0, h.GetMaximum() * 0.9)
    fit.SetParameter(1, 0.0)
    fit.SetParameter(2, 1.0)
    fit.SetParameter(3, h.GetMaximum() * 0.5)
    fit.SetParameter(4, -2.0)
    fit.SetParameter(5, 1.0)
    if useElectrons:
        fit.SetParameter(6, h.GetMaximum() * 0.05)
        #fit.SetParLimits(6, 0.0, h.GetMaximum() * 0.1)
        fit.SetParameter(7, 3.)
        #fit.SetParLimits(7, 1.5, 5.0)  ## this one drops b/g by 10% !
        fit.SetParameter(8, 1.0)
    
    fit.SetParName(0, '#pi magnitude')
    fit.SetParName(1, '#pi mean')
    fit.SetParName(2, '#pi width')
    fit.SetParName(3, 'p/K magnitude')
    fit.SetParName(4, 'p/K mean')
    fit.SetParName(5, 'p/K width')
    if useElectrons:
        fit.SetParName(6, 'ele magnitude')
        fit.SetParName(7, 'ele mean')
        fit.SetParName(8, 'ele width')
    
    h.Fit(fit, 'rq')
    h.DrawCopy()
    
    pifit = ROOT.TF1('pifit', 'gaus', -6.0, 6.0)
    pifit.FixParameter(0, fit.GetParameter(0))
    pifit.FixParameter(1, fit.GetParameter(1))
    pifit.FixParameter(2, fit.GetParameter(2))
    pifit.SetLineColor(ROOT.kRed)
    
    pkfit = ROOT.TF1('pkfit', 'gaus', -6.0, 6.0)
    pkfit.FixParameter(0, fit.GetParameter(3))
    pkfit.FixParameter(1, fit.GetParameter(4))
    pkfit.FixParameter(2, fit.GetParameter(5))
    pkfit.SetLineColor(ROOT.kGreen)
    
    if useElectrons:
        elefit = ROOT.TF1('elefit', 'gaus', -6.0, 6.0)
        elefit.FixParameter(0, fit.GetParameter(6))
        elefit.FixParameter(1, fit.GetParameter(7))
        elefit.FixParameter(2, fit.GetParameter(8))
        elefit.SetLineColor(ROOT.kBlue)
    
    pifit.DrawCopy('same')
    pkfit.DrawCopy('same')
    if useElectrons:
        elefit.DrawCopy('same')
    
    lowBound = pifit.GetParameter(1) - 1.0 * pifit.GetParameter(2)
    highBound = pifit.GetParameter(1) + 2.0 * pifit.GetParameter(2)
    totalPions = pifit.Integral(-6.0, 6.0)
    oldPions = pifit.Integral(-1.0, 2.0)
    newPions = pifit.Integral(lowBound, highBound)
    if useElectrons:
        oldBg = pkfit.Integral(-1.0, 2.0) + elefit.Integral(-1.0, 2.0)
        newBg = pkfit.Integral(lowBound, highBound) + elefit.Integral(lowBound, highBound)
    else:
        oldBg = pkfit.Integral(-1.0, 2.0)
        newBg = pkfit.Integral(lowBound, highBound)
    print 'tot=%7.2f old eff=%.2f new eff=%.2f old bg=%.2f new bg=%.2f' % \
        (totalPions, oldPions/totalPions, newPions/totalPions, oldBg/fit.Integral(-1.0, 2.0), newBg/fit.Integral(lowBound, highBound))
    #print 'old', pifit.Integral(-1.0, 2.0), 'new', pifit.Integral(lowBound, highBound)
    #print 'old', pkfit.Integral(-1.0, 2.0), 'new', pkfit.Integral(lowBound, highBound)
    #print 'old', elefit.Integral(-1.0, 2.0), 'new', elefit.Integral(lowBound, highBound)
    
    
    #print h.GetTitle()[-4:], pifit.GetParameter(1)
    #raw_input('press enter to continue:')
    if useElectrons:
        return fit, pifit, pkfit, elefit
    else:
        return fit, pifit, pkfit


def pid_calibration_allfills(mydir='/Users/kocolosk/data/run5/hist-by-fill'):
    """generates a PDF of triple-Gaussian fits for all fills, plus a histogram of pion means"""
    allFiles = os.listdir(mydir)
    hfill = ROOT.TH1D('hfill','mean of pion Gaussian by RHIC Fill', len(allFiles), 0.5, len(allFiles)+0.5)
    ps = ROOT.TPostScript('pid.ps')
    c = ROOT.TCanvas('c','',100,100,600,800)
    pad = 1
    
    ## some cumulative stats
    nEntries = 0
    nTotalPions = 0
    nOldPions = 0
    nNewPions = 0
    nOldBg = 0
    nNewBg = 0
    nOldCounts = 0
    nNewCounts = 0
    
    myrecords = [] ## fname, pi mean, pi sigma
    counter = 0
    for fname in allFiles:
        if not fname.endswith('.root'): continue
        if counter % 15 == 0:
            c.Update()
            ps.NewPage()
            c.Clear()
            c.Divide(3,5)
            pad = 1
        counter += 1 
        c.cd(pad)    
        print fname
        tfile = ROOT.TFile(os.path.join(mydir, fname))
        mgr = analysis.HistogramManager(tfile, ['nSigmaPion'])
        h = mgr.anyspin['alltrigs'].tracks_sum['nSigmaPion']
        h.SetTitle('n#sigma(#pi) calibration for F%s' % (fname[-14:-10],))
        fit, pifit, pkfit, elefit = pid_calibration(h)
        mean = fit.GetParameter(1)
        error = fit.GetParError(1)
        sigma = fit.GetParameter(2)
        hfill.SetBinContent(counter+1, mean)
        hfill.SetBinError(counter+1, error)
        myrecords.append((fname, mean, sigma))
        pad += 1
        
        ## stats
        lowBound = pifit.GetParameter(1) - 1.0 * pifit.GetParameter(2)
        highBound = pifit.GetParameter(1) + 2.0 * pifit.GetParameter(2)
        binWidth = h.GetBinWidth(1)
        nEntries += h.GetEntries()
        nTotalPions += pifit.Integral(-6.0, 6.0) / binWidth
        nOldPions += pifit.Integral(-1.0, 2.0) / binWidth
        nNewPions += pifit.Integral(lowBound, highBound) / binWidth
        nOldBg += (pkfit.Integral(-1.0, 2.0) + elefit.Integral(-1.0, 2.0)) / binWidth
        nNewBg += (pkfit.Integral(lowBound, highBound) + elefit.Integral(lowBound, highBound)) / binWidth
        nOldCounts += fit.Integral(-1.0, 2.0) / binWidth
        nNewCounts += fit.Integral(lowBound, highBound) / binWidth
    ps.Close()
    c = ROOT.TCanvas()
    hfill.GetYaxis().SetRangeUser(-0.5, 0.8)
    hfill.SetXTitle('fill index')
    hfill.Draw('e')
    
    for r in myrecords:
        print '%d : (% 1.6f, %1.6f),' % (int(r[0][13:17]), r[1], r[2])
    
    print 'Old Total Efficiency    = %.4f' % (nOldPions/nTotalPions,)
    print 'New Total Efficiency    = %.4f' % (nNewPions/nTotalPions,)
    print 'Old Background Fraction = %.4f' % (nOldBg/nOldCounts,)
    print 'New Background Fraction = %.4f' % (nNewBg/nNewCounts,)
    print 'Total Statistics Old = %.0f New %.0f' % (nOldCounts, nNewCounts)
    raw_input('press enter:')


def print_statistics(runlist, keytype='event'):
    """Run 5 only: prints detailed (event/track) statistics for all runs in runlist"""
    a = [0, 0, 0, 0, 0, 0, 0, 0]
    print ' ---------------------------------------------------------------------------------- '
    print '|                  96011   96201   96211   96221   96233 |      HT      JP     ALL |'
    for row,run in enumerate(runlist):
        if row % 10 == 0: 
            print ' ---------------------------------------------------------------------------------- '
        f = ROOT.TFile('~/data/run5/hist/chargedPions_%d.hist.root' % run)
        if keytype=='event':
            h = analysis.HistogramManager(f,'bx7').anyspin
            b = [
                h['96011']['bx7'].GetEntries(),
                h['96201']['bx7'].GetEntries(),
                h['96211']['bx7'].GetEntries(),
                h['96221']['bx7'].GetEntries(),
                h['96233']['bx7'].GetEntries(),
                h['hightower']['bx7'].GetEntries(),
                h['jetpatch']['bx7'].GetEntries(),
                h['alltrigs']['bx7'].GetEntries() 
            ]
        else:
            h = analysis.HistogramManager(f,'pt').anyspin
            b = [
                h['96011'].tracks_sum['pt'].GetEntries(),
                h['96201'].tracks_sum['pt'].GetEntries(),
                h['96211'].tracks_sum['pt'].GetEntries(),
                h['96221'].tracks_sum['pt'].GetEntries(),
                h['96233'].tracks_sum['pt'].GetEntries(),
                h['hightower'].tracks_sum['pt'].GetEntries(),
                h['jetpatch'].tracks_sum['pt'].GetEntries(),
                h['alltrigs'].tracks_sum['pt'].GetEntries() 
            ]
        print '| %d %d : %7d %7d %7d %7d %7d | %7d %7d %7d |' % (analysis.getFill(run), run, \
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7] )
        for i in range(len(b)):
            a[i] += b[i]
    print ' ---------------------------------------------------------------------------------- '
    print '|     sum      : %7d %7d %7d %7d %7d | %7d %7d %7d |' % \
        (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7])


def runlist_luminosity(runlist):
    """prints integrated luminosity seen by minbias trigger in runs of this list"""
    lumiSum = 0
    for run in runlist:
        if run > 7000000:
            f = ROOT.TFile('~/data/run6/tree/chargedPions_%d.tree.root' % (run,))
            minBiasId = 117001
        else:
            path = '/Users/kocolosk/data/run5/tree/chargedPions_%d.tree.root' % run
            if os.path.isfile(path):
                f = ROOT.TFile(path)
            else:
                path = '/Users/kocolosk/data/run5/tree-minbias/chargedPions_%d.tree.root' % run
                f = ROOT.TFile(path)
            minBiasId = 96011
        try:
            tree = f.tree
            lumi = analysis.tree.integratedLuminosity(tree, minBiasId)
            print '%d: %3.6f nb^-1' % (run, lumi)
            lumiSum += lumi
        except AttributeError:
            pass
    print 'Integrated Recorded Luminosity for Runlist: %3.6f pb^-1' % (lumiSum/1000,)


def spinInfoForFrank():
    """not actually a plot"""
    chain = ROOT.TChain('tree')
    chain.Add('~/data/run5/tree/chargedPions_*')
    
    ## only adding these while the other ones spin
    chain.Add('~/data/run5/tree/backup-2008-01-08-trigger-prescales/chargedPions_*')
    
    chain.SetBranchStatus('*',0)
    chain.SetBranchStatus('mRunId',1)
    chain.SetBranchStatus('mEventId',1)
    chain.SetBranchStatus('mSpinBit',1)
    chain.SetBranchStatus('mBx7',1)
    chain.SetBranchStatus('mSpinQA',1)
    
    f = ROOT.TFile('spin_info.root','recreate')
    nt = ROOT.TNtuple('nt','spin info','run:event:spinbit:bx7:qa')
        
    for entry in chain:
        ev = entry.event
        nt.Fill(ev.runId(), ev.eventId(), ev.spinBit(), ev.bx7(), ev.isSpinValid())
    
    nt.Write()
    f.Close()


def ssa_by_fill(asym_key='ly', runlist=analysis.final_runlist_run5):
    """generates canvas for asym_key (ly,lb,ls,us). Each data point is the SSA for a fill"""
    tuples = analysis.getAllFills(runlist)
    fills = []
    
    for run, fill in tuples:
        ## temporary hacks
        ## http://www.star.bnl.gov/HyperNews-star/protected/get/starspin/3324.html
        if 6144002 <= run <= 6144029: fill = 7128
        if 6144041 <= run <= 6144042: fill = 7129        
        if 6145067 <= run <= 6145068: fill = 7136
        if 6146001 <= run <= 6146026: fill = 7138
        if run in runlist:
            fills.append(int(fill))
    fills = analysis.uniqify(fills)
    
    asym_plus = {}
    asym_minus = {}
    for f in fills:
        asym_plus[f] = analysis.AsymmetryGenerator(name='F%s_plus' % f, bins=[1,2.0,10.0])
        asym_minus[f] = analysis.AsymmetryGenerator(name='F%s_minus' % f, bins=[1,2.0,10.0])
    
    scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
    scalars = analysis.ScalarCounts(scalar_path)
    
    polarizations = analysis.Polarizations.Final
    
    allFiles = glob(run5_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            print fname, run
            tfile = ROOT.TFile(fname)
            mgr = analysis.HistogramManager(tfile,('pt',))
            
            try:
                bin7 = scalars[str(run) + '-5-7']
                bin8 = scalars[str(run) + '-5-8']
                bin9 = scalars[str(run) + '-5-9']
            except KeyError:
                print run, 'is not in the scalars database'
                continue
            uu = bin7.uu + bin8.uu + bin9.uu
            ud = bin7.ud + bin8.ud + bin9.ud
            du = bin7.du + bin8.du + bin9.du
            dd = bin7.dd + bin8.dd + bin9.dd
            
            try:
                pol = polarizations[bin7.fill]
            except KeyError:
                print bin7.fill, 'has no final polarization values'
                continue
            
            asym_plus[bin7.fill].FillFromHistogramManager(mgr, 'jetpatch', 1, uu, ud, du, dd, pol.py, pol.pb)
            asym_minus[bin7.fill].FillFromHistogramManager(mgr, 'jetpatch', -1, uu, ud, du, dd, pol.py, pol.pb)
            tfile.Close()
    
    title = {'ly':'Yellow Beam', 'lb':'Blue Beam', 'ls':'Like-Sign', 'us':'Unlike-Sign'}
    final_hist_plus = ROOT.TH1D('final_hist_plus','#pi^{+} %s SSA' % title[asym_key], len(fills), 0.5, len(fills)+0.5)
    final_hist_minus = ROOT.TH1D('final_hist_minus','#pi^{-} %s SSA' % title[asym_key], len(fills), 0.5, len(fills)+0.5)
    marker_color = {'ly':ROOT.kYellow, 'lb':ROOT.kBlue, 'ls':ROOT.kRed, 'us':ROOT.kBlack}
    for h in (final_hist_plus, final_hist_minus):
        h.SetMarkerColor( marker_color[asym_key] )
        h.GetYaxis().SetRangeUser(-0.2, 0.2)
        h.SetXTitle('fill index')
    final_hist_plus.SetMarkerStyle(21)
    final_hist_minus.SetMarkerStyle(20)
    for i,f in enumerate(fills):
        hplus = asym_plus[f].GetAsymmetry(asym_key)
        final_hist_plus.SetBinContent( i+1, hplus.GetBinContent(1) )
        final_hist_plus.SetBinError( i+1, hplus.GetBinError(1) )
        print '%d % .4f % .4f % .4f' % (f, hplus.GetBinContent(1), hplus.GetBinError(1), hplus.GetBinContent(1)/hplus.GetBinError(1))
        hplus.Delete()
        
        hminus = asym_minus[f].GetAsymmetry(asym_key)
        final_hist_minus.SetBinContent( i+1, hminus.GetBinContent(1) )
        final_hist_minus.SetBinError( i+1, hminus.GetBinError(1) )
        hminus.Delete()
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetErrorX(0.0)
    ROOT.gStyle.SetOptFit(111)
    
    cp = ROOT.TCanvas('%s_ssa_fill_plus' % asym_key)
    final_hist_plus.Draw('e1')
    final_hist_plus.Fit('pol0')
    cm = ROOT.TCanvas('%s_ssa_fill_minus' % asym_key)
    final_hist_minus.Draw('e1')
    final_hist_minus.Fit('pol0')
    raw_input('wait here')
    cp.Print('.png')
    cm.Print('.png')


def ssa(asym_key='ly', runlist=analysis.final_runlist_run5, variable='pt'):
    """plots a ssa against the given variable"""
    ## note: need to redefine binning if variable != pt
    asym_plus = analysis.AsymmetryGenerator('asym_plus', key=variable)
    asym_minus = analysis.AsymmetryGenerator('asym_minus', key=variable)
    scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
    scalars = analysis.ScalarCounts(scalar_path)
    polarizations = analysis.Polarizations.Final
    allFiles = glob(run5_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            print fname, run
            tfile = ROOT.TFile(fname)
            mgr = analysis.HistogramManager(tfile,('pt',))
            
            try:
                bin7 = scalars[str(run) + '-5-7']
                bin8 = scalars[str(run) + '-5-8']
                bin9 = scalars[str(run) + '-5-9']
            except KeyError:
                print run, 'is not in the scalars database'
                continue
            uu = bin7.uu + bin8.uu + bin9.uu
            ud = bin7.ud + bin8.ud + bin9.ud
            du = bin7.du + bin8.du + bin9.du
            dd = bin7.dd + bin8.dd + bin9.dd
            
            try:
                pol = polarizations[bin7.fill]
            except KeyError:
                print bin7.fill, 'has no final polarization values'
                continue
            
            asym_plus.FillFromHistogramManager(mgr, 'jetpatch', 1, uu, ud, du, dd, pol.py, pol.pb)
            asym_minus.FillFromHistogramManager(mgr, 'jetpatch', -1, uu, ud, du, dd, pol.py, pol.pb)
            tfile.Close()
            
    title = {'ly':'Yellow Beam', 'lb':'Blue Beam', 'ls':'Like-Sign', 'us':'Unlike-Sign'}
    marker_color = {'ly':ROOT.kYellow, 'lb':ROOT.kBlue, 'ls':ROOT.kRed, 'us':ROOT.kBlack}
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetErrorX(0.0)
    ROOT.gStyle.SetOptFit(111)
    
    hp = asym_plus.GetAsymmetry(asym_key)
    hp.SetTitle(title[asym_key] + ' SSA for #pi^{+}')
    hp.SetMarkerStyle(21)
    
    hm = asym_minus.GetAsymmetry(asym_key)
    hm.SetTitle(title[asym_key] + ' SSA for #pi^{-}')
    hm.SetMarkerStyle(20)
    
    for h in (hp,hm):
        h.SetMarkerColor(marker_color[asym_key])
        h.SetXTitle(variable)
        if asym_key in ('ly', 'lb'):
            h.GetYaxis().SetRangeUser(-0.05, 0.05)
        else:
            h.GetYaxis().SetRangeUser(-0.1, 0.1)
    
    cp = ROOT.TCanvas('cp')
    hp.Fit('pol0')
    hp.Draw('e1')
    
    cm = ROOT.TCanvas('cm')
    hm.Fit('pol0')
    hm.Draw('e1')
    
    raw_input('wait here')
    cp.Print('%s_ssa_%s_plus.png' % (asym_key, variable))
    cm.Print('%s_ssa_%s_minus.png' % (asym_key, variable))
    


def pid_background_asymmetry(runlist=analysis.final_runlist_run5):
    """plots A_{LL} for charged tracks outside PID window and fits with a pol0"""
    asym_plus = analysis.AsymmetryGenerator('asym_plus', key='pt_bg')
    asym_minus = analysis.AsymmetryGenerator('asym_minus', key='pt_bg')
    scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
    scalars = analysis.ScalarCounts(scalar_path)
    polarizations = analysis.Polarizations.Final
    allFiles = glob(run5_hist_dir + '/chargedPions_*.hist.root')
    for fname in allFiles:
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            print fname, run
            tfile = ROOT.TFile(fname)
            mgr = analysis.HistogramManager(tfile,('pt_bg',))
            
            try:
                bin7 = scalars[str(run) + '-5-7']
                bin8 = scalars[str(run) + '-5-8']
                bin9 = scalars[str(run) + '-5-9']
            except KeyError:
                print run, 'is not in the scalars database'
                continue
            uu = bin7.uu + bin8.uu + bin9.uu
            ud = bin7.ud + bin8.ud + bin9.ud
            du = bin7.du + bin8.du + bin9.du
            dd = bin7.dd + bin8.dd + bin9.dd
            
            try:
                pol = polarizations[bin7.fill]
            except KeyError:
                print bin7.fill, 'has no final polarization values'
                continue
            
            asym_plus.FillFromHistogramManager(mgr, 'jetpatch', 1, uu, ud, du, dd, pol.py, pol.pb)
            asym_minus.FillFromHistogramManager(mgr, 'jetpatch', -1, uu, ud, du, dd, pol.py, pol.pb)
            tfile.Close()
            
    #title = {'ly':'Yellow Beam', 'lb':'Blue Beam', 'ls':'Like-Sign', 'us':'Unlike-Sign'}
    #marker_color = {'ly':ROOT.kYellow, 'lb':ROOT.kBlue, 'ls':ROOT.kRed, 'us':ROOT.kBlack}
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetErrorX(0.0)
    ROOT.gStyle.SetOptFit(111)
    
    hp = asym_plus.GetAsymmetry('ll')
    hp.SetTitle('PID Background A_{LL} for #pi^{+}')
    hp.SetMarkerStyle(21)
    
    hm = asym_minus.GetAsymmetry('ll')
    hm.SetTitle('PID Background A_{LL} for #pi^{-}')
    hm.SetMarkerStyle(20)
    
    for h in (hp,hm):
        #h.SetMarkerColor(marker_color[asym_key])
        h.SetXTitle('p_{T}')
       # if asym_key in ('ly', 'lb'):
        h.GetYaxis().SetRangeUser(-0.05, 0.05)
        #else:
        #    h.GetYaxis().SetRangeUser(-0.1, 0.1)
    
    cp = ROOT.TCanvas('cp')
    hp.Fit('pol0')
    hp.Draw('e1')
    
    cm = ROOT.TCanvas('cm')
    hm.Fit('pol0')
    hm.Draw('e1')
    
    raw_input('wait here')
    cp.Print('pid_background_asymmetry_plus.png')
    cm.Print('pid_background_asymmetry_minus.png')


def print_ssa(runlist=analysis.final_runlist_run5, charge=1):
    """prints (val +/- err) => (nSig) for y,b,ls,us for each run in list"""
    
    scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
    scalars = analysis.ScalarCounts(scalar_path)
    polarizations = analysis.Polarizations.Final
    
    print '------------------------------------------------------------------------------------\
---------------------------------------------------------------'
    print ' Fill |   Run    |          Yellow Beam           |           Blue Beam            |\
           Like Sign            |          Unlike Sign          '
    print '------------------------------------------------------------------------------------\
---------------------------------------------------------------'
    for run in runlist:
        gen = analysis.AsymmetryGenerator(name=str(run), bins=[1,2.0,10.0])
        f = ROOT.TFile(run5_hist_dir + '/chargedPions_%d.hist.root' % run)
        mgr = analysis.HistogramManager(f,('pt',))
        
        try:
            bin7 = scalars[str(run) + '-5-7']
            bin8 = scalars[str(run) + '-5-8']
            bin9 = scalars[str(run) + '-5-9']
        except KeyError:
            print run, 'is not in the scalars database'
            continue
        uu = bin7.uu + bin8.uu + bin9.uu
        ud = bin7.ud + bin8.ud + bin9.ud
        du = bin7.du + bin8.du + bin9.du
        dd = bin7.dd + bin8.dd + bin9.dd
        
        try:
            pol = polarizations[bin7.fill]
        except KeyError:
            print bin7.fill, 'has no final polarization values'
            continue
        
        gen.FillFromHistogramManager(mgr, 'jetpatch', charge, uu, ud, du, dd, pol.py, pol.pb)
        f.Close()
        y = gen.GetAsymmetry('ly')
        line = 'F%d | R%d | (% .5f +- %.5f) => % .2f |' % (bin7.fill, run, y.GetBinContent(1), \
            y.GetBinError(1), (y.GetBinContent(1)/y.GetBinError(1)) )
        b = gen.GetAsymmetry('lb')
        line += ' (% .5f +- %.5f) => % .2f |' % (b.GetBinContent(1), \
            b.GetBinError(1), (b.GetBinContent(1)/b.GetBinError(1)))
        ls = gen.GetAsymmetry('ls')
        line += ' (% .5f +- %.5f) => % .2f |' % (ls.GetBinContent(1), \
            ls.GetBinError(1), (ls.GetBinContent(1)/ls.GetBinError(1)))
        us = gen.GetAsymmetry('us')
        line += ' (% .5f +- %.5f) => % .2f' % (us.GetBinContent(1), \
            us.GetBinError(1), (us.GetBinContent(1)/us.GetBinError(1)))
        print line
        [ h.Delete() for h in (y,b,ls,us) ]


def asigma(runlist=analysis.transverse_run5):
    """plots asigma -- duh"""
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetErrorX(0.0)
    ROOT.gStyle.SetOptFit(111)
    asym_plus = analysis.AsymmetryGenerator('asym_plus', key='pt')
    asym_minus = analysis.AsymmetryGenerator('asym_minus', key='pt')
    #scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
    scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run6.txt'
    scalars = analysis.ScalarCounts(scalar_path)
    polarizations = analysis.Polarizations.Final
    allFiles = glob(run5_hist_dir + '-transverse/chargedPions_*.hist.root')
    allFiles += glob(run6_hist_dir + '-transverse/chargedPions_*.hist.root')
    print allFiles
    for fname in allFiles:
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            print fname, run
            tfile = ROOT.TFile(fname)
            mgr = analysis.HistogramManager(tfile,('pt',))
            
            try:
                bin7 = scalars[str(run) + '-5-7']
                bin8 = scalars[str(run) + '-5-8']
                bin9 = scalars[str(run) + '-5-9']
            except KeyError:
                print run, 'is not in the scalars database'
                continue
            uu = bin7.uu + bin8.uu + bin9.uu
            ud = bin7.ud + bin8.ud + bin9.ud
            du = bin7.du + bin8.du + bin9.du
            dd = bin7.dd + bin8.dd + bin9.dd
            
            try:
                pol = polarizations[bin7.fill]
            except KeyError:
                print bin7.fill, 'has no final polarization values'
                continue
            
            asym_plus.FillFromHistogramManager(mgr, 'jetpatch', 1, uu, ud, du, dd, pol.py, pol.pb)
            asym_minus.FillFromHistogramManager(mgr, 'jetpatch', -1, uu, ud, du, dd, pol.py, pol.pb)
            tfile.Close()
    
    hp = asym_plus.GetAsymmetry('ll')
    hp.SetTitle('A_{#Sigma} for #pi^{+} using Run 6 data')
    hp.SetMarkerStyle(21)
    
    hm = asym_minus.GetAsymmetry('ll')
    hm.SetTitle('A_{#Sigma} for #pi^{-} using Run 6 data')
    hm.SetMarkerStyle(20)
    
    for h in (hp,hm):
        h.SetXTitle('p_{T}')
        h.GetYaxis().SetRangeUser(-1.0, 1.0)
    
    cp = ROOT.TCanvas('cp')
    hp.Fit('pol0')
    hp.Draw('e1')
    
    cm = ROOT.TCanvas('cm')
    hm.Fit('pol0')
    hm.Draw('e1')
    
    raw_input('wait here')
    cp.Print('asigma_plus.png')
    cm.Print('asigma_minus.png')


def pid_pt_dependence(runlist=analysis.final_runlist_run5, tfile=None):
    """plots recalibrated nSigmaPion for charge-summed pions in each pT bin"""
    if tfile is None:
        nsig = [ ROOT.TH1D('nsig_%d' % i,'',100,-6.0,6.0) for i in range(5)]
        nsig[0].SetTitle('Recalibrated n#sigma(#pi), 2<p_{T}<10')
        nsig[1].SetTitle('Recalibrated n#sigma(#pi), 8<p_{T}<10')
        nsig[2].SetTitle('Recalibrated n#sigma(#pi), 6<p_{T}<8')
        nsig[3].SetTitle('Recalibrated n#sigma(#pi), 4<p_{T}<6')
        nsig[4].SetTitle('Recalibrated n#sigma(#pi), 2<p_{T}<4')
        ptbins = [2.0, 4.0, 6.0, 8.0]
        ptbins.reverse()
        ecuts = analysis.histos.EventCuts()
        for run in runlist:
            fill = analysis.getFill(run)
            tcuts = analysis.histos.TrackCuts(fill)
            pidFit = analysis.histos.pidCalibration[fill]
            if run > 7000000:
                fname = run6_tree_dir + '/chargedPions_%d.tree.root' % run
            else:
                fname = run5_tree_dir + '/chargedPions_%d.tree.root' % run
                if not os.path.isfile(fname):
                    fname = run5_tree_dir + '-minbias/chargedPions_%d.tree.root' % run
            print fname
            f = ROOT.TFile(fname)
            for t in f.tree:
                ecuts.set(t.event)
                if not ecuts.all: continue
                for track in t.event.tracks():
                    tcuts.set(track)
                    if tcuts.eta and tcuts.dca and tcuts.fit:
                        ## recalibrate nSigmaPion
                        nSigmaPion = (track.nSigmaPion() - pidFit[0]) / pidFit[1]
                        nsig[0].Fill(nSigmaPion)
                        for i,ptcut in enumerate(ptbins):
                            if track.pt() > ptcut:
                                nsig[i+1].Fill(nSigmaPion)
                                break
    else:
        nsig = [ tfile.Get('nsig_%d' % i) for i in range(5) ]
    
    
    c = []
    outfile = ROOT.TFile('/tmp/pid_pt_dependence.root','recreate')
    for h in nsig:
        fits = pid_calibration(h)
        #fits = pid_calibration(h, useElectrons=(h.GetName() != 'nsig_1'))
        [ h.GetListOfFunctions().Add(f) for f in fits[1:] ]
        c.append(ROOT.TCanvas())
        h.Draw()
        h.Write()
        c[-1].Print('.png')
    
    raw_input('wait here:')
    outfile.Close()

