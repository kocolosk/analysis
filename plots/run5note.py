import ROOT
from analysis.plots import graphics
from analysis.histos2 import HistogramManager

simuFile = '/Users/kocolosk/data/run5-simu/hist/merged.cphist.root'

def mcjet_ptprofile(fname='ptprofiles.cphist.root'):
    ROOT.gStyle.SetErrorX(False)
    
    mgr = HistogramManager(ROOT.TFile(fname))
    minbias = mgr['anyspin']['96011']['sum']['ptprofiles_mcjetpt']
    jp1     = mgr['anyspin']['96221']['sum']['ptprofiles_mcjetpt']
    jp2     = mgr['anyspin']['96233']['sum']['ptprofiles_mcjetpt']
    
    jp1.SetMarkerStyle(27)
    jp1.SetMarkerColor(ROOT.kBlue)
    jp1.SetLineColor(ROOT.kBlue)
    
    jp2.SetMarkerStyle(25)
    jp2.SetMarkerColor(ROOT.kRed)
    jp2.SetLineColor(ROOT.kRed)
    
    minbias.SetTitle('')
    minbias.SetMarkerStyle(21)
    minbias.GetXaxis().SetRangeUser(1, 12)
    minbias.GetYaxis().SetRangeUser(3, 25)
    
    c = graphics.canvas1()
    minbias.Draw('e1')
    jp2.Draw('same e1')
    jp1.Draw('same e1')
    
    leg = ROOT.TLegend(0.15, 0.65, 0.3, 0.85)
    leg.AddEntry(minbias.obj, 'MB')
    leg.AddEntry(jp1.obj, 'JP1')
    leg.AddEntry(jp2.obj, 'JP2')
    leg.Draw()
    
    raw_input('wait here:')

def jet_ptprofile(fname='ptprofiles.cphist.root'):
    ROOT.gStyle.SetErrorX(False)
    
    mgr = HistogramManager(ROOT.TFile(fname))
    minbias = mgr['anyspin']['96011']['sum']['ptprofiles_jetpt']
    jp1     = mgr['anyspin']['96221']['sum']['ptprofiles_jetpt']
    jp2     = mgr['anyspin']['96233']['sum']['ptprofiles_jetpt']
    
    jp1.SetMarkerStyle(27)
    jp1.SetMarkerColor(ROOT.kBlue)
    jp1.SetLineColor(ROOT.kBlue)
    
    jp2.SetMarkerStyle(25)
    jp2.SetMarkerColor(ROOT.kRed)
    jp2.SetLineColor(ROOT.kRed)
    
    minbias.SetTitle('')
    minbias.SetMarkerStyle(21)
    minbias.GetXaxis().SetRangeUser(1, 12)
    minbias.GetYaxis().SetRangeUser(3, 25)
    
    c = graphics.canvas1()
    minbias.Draw('e1')
    jp2.Draw('same e1')
    jp1.Draw('same e1')
    
    leg = ROOT.TLegend(0.15, 0.65, 0.3, 0.85)
    leg.AddEntry(minbias, 'MB')
    leg.AddEntry(jp1, 'JP1')
    leg.AddEntry(jp2, 'JP2')
    leg.Draw()
    
    raw_input('wait here:')

def x_ptprofile(fname='ptprofiles.cphist.root'):
    ROOT.gStyle.SetErrorX(False)
    
    mgr = HistogramManager(ROOT.TFile(fname))
    minbias = mgr['anyspin']['96011']['sum']['ptprofiles_x']
    jp1     = mgr['anyspin']['96221']['sum']['ptprofiles_x']
    jp2     = mgr['anyspin']['96233']['sum']['ptprofiles_x']
    
    jp1.SetMarkerStyle(27)
    jp1.SetMarkerColor(ROOT.kBlue)
    jp1.SetLineColor(ROOT.kBlue)
    
    jp2.SetMarkerStyle(25)
    jp2.SetMarkerColor(ROOT.kRed)
    jp2.SetLineColor(ROOT.kRed)
    
    minbias.SetTitle('')
    minbias.SetMarkerStyle(21)
    minbias.GetXaxis().SetRangeUser(1, 12)
    minbias.GetYaxis().SetRangeUser(0.03, 0.25)
    
    c = graphics.canvas1()
    minbias.Draw('e1')
    jp2.Draw('same e1')
    jp1.Draw('same e1')
    
    leg = ROOT.TLegend(0.15, 0.65, 0.3, 0.85)
    leg.AddEntry(minbias, 'MB')
    leg.AddEntry(jp1, 'JP1')
    leg.AddEntry(jp2, 'JP2')
    leg.Draw()
    
    raw_input('wait here:')

def subprocess_shift(fname = simuFile):
    mgr = HistogramManager(ROOT.TFile(fname))
    mb = {
        'all': mgr.anyspin['96011']['hardP'],
        'gg': mgr.gg['96011']['hardP'],
        'qg': mgr.qg['96011']['hardP'],
        'qq': mgr.qq['96011']['hardP']
    }
    jp = {
        'all': mgr.anyspin['96233']['hardP'],
        'gg': mgr.gg['96233']['hardP'],
        'qg': mgr.qg['96233']['hardP'],
        'qq': mgr.qq['96233']['hardP']
    }
    
    for key in ('gg','qg','qq'):
        mb[key].Divide(mb['all'])
        jp[key].Divide(jp['all'])
    
    for h in mb.values():
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.8)
    
    for h in jp.values():
        h.SetMarkerStyle(24)
        h.SetMarkerSize(0.8)
    
    for h in (mb,jp):
        h['gg'].SetLineColor(ROOT.kRed)
        h['gg'].SetMarkerColor(ROOT.kRed)
        h['qg'].SetLineColor(ROOT.kBlue)
        h['qg'].SetMarkerColor(ROOT.kBlue)
        h['qq'].SetLineColor(ROOT.kGreen)
        h['qq'].SetMarkerColor(ROOT.kGreen)
    
    mb['gg'].SetTitle('Subprocess Fraction')
    mb['gg'].SetXTitle('partonic p_{T}')
    mb['gg'].GetYaxis().SetRangeUser(0., 0.7)
    mb['gg'].GetXaxis().SetRangeUser(3, 30)
    
    c = graphics.canvas1()
    mb['gg'].Draw('e1')
    mb['qg'].Draw('e1 same')
    mb['qq'].Draw('e1 same')
    
    jp['gg'].Draw('e1 same')
    jp['qg'].Draw('e1 same')
    jp['qq'].Draw('e1 same')
    
    leg = ROOT.TLegend(0.6, 0.8, 0.88, 0.88, 'filled = MB, open = JP2')
    leg.SetNColumns(3)
    leg.AddEntry(mb['gg'], 'gg')
    leg.AddEntry(mb['qg'], 'qg')
    leg.AddEntry(mb['qq'], 'qq')
    leg.Draw()
    
    raw_input('wait here:')
