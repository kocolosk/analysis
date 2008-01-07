import os
import math
from glob import glob

import ROOT
import analysis

run5_tree_dir = '/Users/kocolosk/data/run5/tree'
run5_hist_dir = '/Users/kocolosk/data/run5/hist'

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


def asymmetries_for_publication_run5():
   """final results for inclusive asymmetries from Run 5"""
   asym_plus = analysis.AsymmetryGenerator('asym_plus')
   asym_minus = analysis.AsymmetryGenerator('asym_minus')
   
   runlist = None
   
   scalar_path = os.environ['STAR'] + '/StRoot/StSpinPool/StTamuRelLum/inputs/run5.txt'
   scalars = analysis.ScalarCounts(scalar_path)
   
   polarizations = analysis.Polarizations.Final
   
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
            print bin7.fill, 'has no final polarization values'
            continue
         
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
   
   raw_input('wait here:')


def eta_phi_jet_correlation_run5():
   """3-D deta-dphi plot -- possible paper plot at one time"""
   """this is broken -- need to add deta-dphi correlations to mgr"""
   style = ROOT.TStyle()
   style.SetOptStat(0)    
   style.SetLabelOffset(-0.01,'xy')
   style.SetLabelSize(0.035,'xy')
   style.SetTitleOffset(1.2,'y')
   style.cd()
   
   fig3 = ROOT.TH2D('figure3','',50,-1.5*math.pi,0.5*math.pi,50,-2.0,0.4)
   fig3.SetYTitle('#eta pion - #eta jet')
   fig3.SetXTitle('#phi pion - #phi jet')
   
   fig3b = ROOT.TH2D('figure3b','',50,-math.pi,math.pi,50,-1.5,1.5)
   
   #reset some styles
   style.SetLabelOffset(0.005,'xy')
   #ROOT.gStyle.SetLabelSize(0.04,'xy')
   style.SetTitleOffset(1,'y')
   
   xbins = [2.0, 3.0, 4.0, 5.5, 7.0, 10.0]
   ar = array('d',xbins)
   fig3c = ROOT.TH2D('figure3c','Uncorrected Pion Momentum Fraction',len(xbins)-1,ar,50,0.,1.)
   fig3c_away = ROOT.TH2D('figure3c_away','',len(xbins)-1,ar,50,0.,1.)
   
   pions = 0
   
   timer = ROOT.TStopwatch()
   for i in xrange(reader.GetEntries()):
       if i == 0:
           timer.Stop()
           print 'generated a TEventList after %f CPU seconds and %f real seconds' % (timer.CpuTime(),timer.RealTime())
           timer.Start(True)
       
       if i % 200000 == 0: print 'processed %d of %d entries' % (i,reader.GetEntries())
       reader.GetEntry(i)
       
       t1 = reader.event().trigger(96221)
       t2 = reader.event().trigger(96233)
       if (t1 is not None and t1.didFire()) or (t2 is not None and t2.didFire()):
           for j in range(reader.nJets()):
               for k in range(reader.nChargedPions()):
                   deta = -1*(reader.jet(j).Eta() - reader.chargedPion(k).Eta())
                   dphi = -1*(reader.jet(j).Phi() - reader.chargedPion(k).Phi())
           
                   #normalize dphi
                   if dphi < -1.5*math.pi:
                       dphi = dphi + 2*math.pi
                   if dphi > 0.5*math.pi:
                       dphi = dphi - 2*math.pi
           
                   fig3.Fill(dphi,deta)
           
                   #renormalize
                   if dphi < -math.pi:
                       dphi = dphi + 2*math.pi
                   if dphi > math.pi:
                       dphi = dphi - 2*math.pi
               
                   if j == 0: pions = pions + 1
           
                   fig3b.Fill(dphi,deta)
                   #figure3.Fill(reader.chargedPion(k).Phi(),reader.chargedPion(k).Eta())
                   
                   dR = math.sqrt(deta*deta + dphi*dphi)
                   if dR < 0.4: fig3c.Fill(reader.chargedPion(k).Pt(),(reader.chargedPion(k).Pt()/reader.jet(j).Pt()))
                   elif dR > 1.5: fig3c_away.Fill(reader.chargedPion(k).Pt(),(reader.chargedPion(k).Pt()/reader.jet(j).Pt()))
   
   print 'processed the events after %f CPU seconds and %f real seconds' % (timer.CpuTime(),timer.RealTime())
   
   print 'counted %d pions' % (pions,)
   
   c2 = ROOT.TCanvas('c2')
   c2.SetLogz()
   fig3.Draw('lego2')
   
   c3 = ROOT.TCanvas('c3')
   #c3.SetLogz()
   fig3c.Draw('col')
   fig3c.FitSlicesY()
   fig3c_mean = ROOT.gDirectory.Get('figure3c_1')
   fig3c_mean.SetTitle(fig3c.GetTitle())
   fig3c_mean.SetXTitle('#pi p_{T}')
   fig3c_mean.SetYTitle('< p_{T,#pi} / p_{T,jet} >')
   fig3c_mean.SetAxisRange(0,1,'y')
   fig3c_mean.SetMarkerStyle(21)
   
   fig3c_away.FitSlicesY()
   fig3c_away_mean = ROOT.gDirectory.Get('figure3c_away_1')
   fig3c_away_mean.SetMarkerStyle(25)
   fig3c_away_mean.SetMarkerColor(ROOT.kRed)
   fig3c_away_mean.SetLineColor(ROOT.kRed)
   
   leg = ROOT.TLegend(0.75,0.8,0.89,0.89)
   leg.AddEntry(fig3c_mean,'dR < 0.4','p')
   leg.AddEntry(fig3c_away_mean,'dR > 1.5','p')
   
   c4 = ROOT.TCanvas('c4')
   #reset some styles
   fig3c_mean.Draw()
   fig3c_away_mean.Draw('same')
   
   leg.Draw('same')    
   
   raw_input('press enter to continue: ')