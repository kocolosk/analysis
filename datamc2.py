import ROOT
import minimc
import histos

def driver(trigId):
   dFile = ROOT.TFile('~/data/run5/allruns.hist.root')
   mFile = ROOT.TFile('~/data/run5/minimc/combined.plus.hist.root')
   
   mcTrigId = str(trigId)
   if trigId == 96011: mcTrigId = 'minbias'
   
   mgr = histos.HistogramManager(dFile)
   mc  = minimc.MiniMcHistos(mcTrigId,'',mFile)
   data = mgr['anyspin'][str(trigId)]
   compare(data,mc)
   

def compare(data, mc, charge = 1):
   """plot var for (Track)HistogramCollection and MiniMcHistos on same canvas"""
   c = ROOT.TCanvas()
   c.Divide(3,2)
   c.cd(1)
   
   ratioc = ROOT.TCanvas('ratioc')
   ratioc.Divide(3,2)
   
   c2 = ROOT.TCanvas('c2','p_{T} Spectra for Charged Pions',800,400)
   c2.Divide(2,1)
   
   c.cd(1)
   
   dvz = data['vz'].Clone()
   dvz.Rebin()
   dvz.Rebin()
   dvz.Sumw2()
   dvz.DrawNormalized()
   
   mvz = mc['vz'].Clone()
   mvz.Rebin()
   mvz.Rebin()
   mvz.SetLineColor(ROOT.kRed)
   mvz.DrawNormalized('same')
   
   ratioc.cd(1)
   ratiovz = dvz.Clone()
   ratiovz.Divide(mvz)
   ratiovz.Scale(mvz.Integral() / dvz.Integral())
   ratiovz.GetYaxis().SetRangeUser(0.,3.0)
   ratiovz.Draw()
   
   for counter,var in enumerate( ('pt','eta','phi','dcaG','nHitsFit') ):
      c.cd(counter+2)
      h = data.trackHistograms(charge)[var]
      h.Sumw2()
      h.SetXTitle(var)
      h.DrawNormalized()
      mc[var].SetLineColor(ROOT.kRed)
      mc[var].DrawNormalized('same')
      
      if var in ('pt'):
         ROOT.gPad.SetLogy()
         
      ratioc.cd(counter+2)
      ratio = h.Clone()
      ratio.Divide(mc[var])
      ratio.Scale(mc[var].Integral() / h.Integral())
      ratio.GetYaxis().SetRangeUser(0.,3.0)
      ratio.Draw()
      
      ## special treatment for pt normalization over smaller range
      if var == 'pt':
         c.cd(counter+2)
         dataBin1 = h.FindBin(2.01)
         dataBin2 = h.FindBin(8.01)
         
         mcBin1   = mc['ptRescaled'].FindBin(2.01)
         mcBin2   = mc['ptRescaled'].FindBin(8.01)
         
         dataIntegral = h.Integral(dataBin1, dataBin2)
         mcIntegral   = mc['ptRescaled'].Integral(mcBin1, mcBin2)
         h.Scale(1.0/dataIntegral)
         mc['ptRescaled'].Scale(1.0/mcIntegral)
         
         h.Draw()
         mc['ptRescaled'].Draw('same')
         
         pad = ratioc.cd(counter+2)
         #ratio.Scale(mcIntegral / (30*dataIntegral))
         ratio = h.Clone()
         ratio.Divide(mc['ptRescaled'])
         ratio.GetYaxis().SetRangeUser(0.,3.0)
         ratio.Draw()
         
         c2.cd(1)
         ROOT.gPad.SetLogy()
         h.Draw()
         mc['ptRescaled'].Draw('same')
         c2.cd(2)
         ratio.Draw()
         
   
   
   raw_input('blah')
   c.Print('comp.eps')
   ratioc.Print('ratio.eps')
   c2.Print('pt.eps')

