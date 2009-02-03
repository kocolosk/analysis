# encoding: utf-8

import ROOT
from analysis.plots import graphics
from analysis.histos2 import HistogramManager
from analysis.util import hadd_interactive
from analysis.runlists import final_runlist_run5 as runlist

simuFile = '/Users/kocolosk/data/run5-simu/hist/merged.cphist.root'
mcasymFile = '/Users/kocolosk/data/run6-simu/hist/mcasym_10.cphist.root'
histDir = '/Users/kocolosk/data/run5/hist'

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
    c.Print('jetpt-profile.png')

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
    c.Print('mcjet-profile.png')

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
    c.Print('x-profile.png')

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
    c.Print('subprocess-shift.png')

def track_cuts_summary():
    canvas = [ graphics.canvas1() for i in range(4) ]
    
    def dostuff(name, xmin, xmax):
        h = hadd_interactive(histDir, runlist, 'jetpatch', 'anyspin', 'sum', name)
        integral = h.Integral()
        bin1 = h.FindBin(xmin)
        bin2 = h.FindBin(xmax)
        h2 = h.Clone()
        for bin in range(1, bin1):
            h2.SetBinContent(bin, 0)
        for bin in range(bin2+1, h.GetNbinsX()):
            h2.SetBinContent(bin, 0)
        h2.SetFillColor(ROOT.kGreen)
        efficiency = 100*h2.Integral()/integral
        h.Draw()
        h.SetTitle(h.GetTitle() + ' -- %.0f%% pass' % (efficiency,))
        h2.Draw('same')
    
    canvas[0].cd()
    dostuff('eta', -0.99, 0.99)
    
    canvas[1].cd()
    dostuff('dcaG', 0., 0.99)
    
    canvas[2].cd()
    dostuff('nHitsFit', 25, 100)
    
    canvas[3].cd()
    dostuff('nSigmaPion', -0.99, 1.99)
    
    raw_input('wait here:')
    canvas[0].Print('eta-cut.png')
    canvas[1].Print('nHitsFit-cut.png')
    canvas[2].Print('dcaG-cut.png')
    canvas[3].Print('nSigmaPion-cut.png')

def z_profile():
    ROOT.gStyle.SetErrorX(False)
    
    mgr = HistogramManager(ROOT.TFile('zprofiles.cphist.root'))
    minbias = mgr['anyspin']['96011']['sum']['ptprofiles_z']
    jp1     = mgr['anyspin']['96221']['sum']['ptprofiles_z']
    jp2     = mgr['anyspin']['96233']['sum']['ptprofiles_z']
    
    jp1.SetMarkerStyle(27)
    jp1.SetMarkerColor(ROOT.kBlue)
    jp1.SetLineColor(ROOT.kBlue)
    
    jp2.SetMarkerStyle(25)
    jp2.SetMarkerColor(ROOT.kRed)
    jp2.SetLineColor(ROOT.kRed)
    
    minbias.SetTitle('')
    minbias.SetMarkerStyle(21)
    minbias.GetXaxis().SetRangeUser(1, 12)
    minbias.GetYaxis().SetRangeUser(0.1, 1.0)
    
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
    c.Print('z-profile.png')

def mcz_profile():
    ROOT.gStyle.SetErrorX(False)
    
    mgr = HistogramManager(ROOT.TFile('zprofiles.cphist.root'))
    minbias = mgr['anyspin']['96011']['sum']['ptprofiles_mcz']
    jp1     = mgr['anyspin']['96221']['sum']['ptprofiles_mcz']
    jp2     = mgr['anyspin']['96233']['sum']['ptprofiles_mcz']
    
    jp1.SetMarkerStyle(27)
    jp1.SetMarkerColor(ROOT.kBlue)
    jp1.SetLineColor(ROOT.kBlue)
    
    jp2.SetMarkerStyle(25)
    jp2.SetMarkerColor(ROOT.kRed)
    jp2.SetLineColor(ROOT.kRed)
    
    minbias.SetTitle('')
    minbias.SetMarkerStyle(21)
    minbias.GetXaxis().SetRangeUser(1, 12)
    minbias.GetYaxis().SetRangeUser(0.1, 1.0)
    
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
    c.Print('mcz-profile.png')

def mcasym(spin = 'anyspin'):
    """
    comparison of MC asymmetries for minbias and 137222
    """
    if spin == 'anyspin':
        stitle = ''
    else:
        stitle = '_%(spin)s' % locals()
    
    keys = ['STD','MIN','ZERO','MAX']
    mgr = HistogramManager(ROOT.TFile(mcasymFile))
    
    line = ROOT.TLine(0, 0.0, 12, 0.0)
    line.SetLineStyle(2)
    ROOT.gStyle.SetErrorX()
    
    color = {
        'STD': ROOT.kBlack,
        'MAX': ROOT.kRed,
        'MIN': ROOT.kGreen,
        'ZERO': ROOT.kBlue,
        'GS_NLOC': ROOT.kMagenta
    }
    
    smooth_factor = 1
    
    cmb = graphics.canvas2()
    # cmbw = graphics.canvas2()
    cjp = graphics.canvas2()
    alldiffs = graphics.canvas2('Asymmetry Differences')
    # alldiffsw = graphics.canvas2('Reweighted Asymmetry Differences')
    keepme = []
    diffs = {}
    for i,key in enumerate(keys):
        mb_m = mgr[spin]['117001'].tracks_minus[key].Clone()
        mb_p = mgr[spin]['117001'].tracks_plus[key].Clone()
        jp_m = mgr[spin]['jetpatch'].tracks_minus[key].Clone()
        jp_p = mgr[spin]['jetpatch'].tracks_plus[key].Clone()
        # mb_mw = mgr[spin]['96011'].tracks_minus[key+'w'].Clone()
        # mb_pw = mgr[spin]['96011'].tracks_plus[key+'w'].Clone()
        
        opt = i>0 and 'e2 same' or 'e2'
        
        # for h in (mb_m, mb_p, jp_m, jp_p, mb_mw, mb_pw):
        for h in (mb_m, mb_p, jp_m, jp_p):
            h.Smooth(smooth_factor)
            h.GetXaxis().SetTitle('p_{T}(#pi)/p_{T}(jet)')
            h.GetXaxis().SetRangeUser(0, 12)
            h.GetYaxis().SetRangeUser(-0.08,0.08)
            h.SetLineColor(color[key])
            h.SetFillColor(color[key])
        
        cmb.cd(1)
        line.Draw()
        mb_m.SetTitle('MB MC Asymmetries for #pi^{-} %(stitle)s' % locals())
        mb_m.Draw(opt)
        
        cmb.cd(2)
        line.Draw()
        mb_p.SetTitle('MB MC Asymmetries for #pi^{+} %(stitle)s' % locals())
        mb_p.Draw(opt)
        
        # cmbw.cd(1)
        # line.Draw()
        # mb_mw.SetTitle('Reweighted MB MC Asymmetries for #pi^{-} %(stitle)s' \
        #     % locals())
        # mb_mw.Draw(opt)
        # 
        # cmbw.cd(2)
        # line.Draw()
        # mb_pw.SetTitle('Reweighted MB MC Asymmetries for #pi^{+} %(stitle)s' \
        #     % locals())
        # mb_pw.Draw(opt)
        
        cjp.cd(1)
        line.Draw()
        jp_m.SetTitle('JP1 MC Asymmetries for #pi^{-} %(stitle)s' % locals())
        jp_m.Draw(opt)
        
        cjp.cd(2)
        line.Draw()
        jp_p.SetTitle('JP1 MC Asymmetries for #pi^{+} %(stitle)s' % locals())
        jp_p.Draw(opt)
        
        diff_m = jp_m.Clone()
        diff_m.Add(mb_m, -1)
        
        # diff_mw = jp_m.Clone()
        # diff_mw.Add(mb_mw, -1)
        
        diff_p = jp_p.Clone()
        diff_p.Add(mb_p, -1)
        
        # diff_pw = jp_p.Clone()
        # diff_pw.Add(mb_pw, -1)
        
        for h in (diff_m, diff_p):
            h.GetYaxis().SetRangeUser(-0.03, 0.03)
            h.GetXaxis().SetRangeUser(0, 12)
            h.SetMarkerStyle(20)
        
        
        # for h in (diff_mw, diff_pw):
        #     h.GetYaxis().SetRangeUser(-0.03, 0.03)
        #     h.GetXaxis().SetRangeUser(zbins[0], zbins[-1])
        #     h.SetMarkerStyle(24)
        
        cdiff = graphics.canvas2(key)
        
        leg = ROOT.TLegend(0.15, 0.15, 0.6, 0.35)
        leg.SetHeader('Compare triggered asym to:')
        leg.AddEntry(diff_m, 'MB', 'p')
        # leg.AddEntry(diff_mw, 'p_{T} reweighted MB', 'p')
        
        cdiff.cd(1)
        if spin == 'anyspin':
            mtitle = '#pi^{-} JP1 - MB for %(key)s' % locals()
            ptitle = '#pi^{+} JP1 - MB for %(key)s' % locals()            
        else:
            mtitle = '#pi^{-} JP1 - MB for %(key)s -- %(spin)s processes only' \
                % locals()
            ptitle = '#pi^{+} JP1 - MB for %(key)s -- %(spin)s processes only' \
                % locals()
        
        diff_m.SetTitle(mtitle)
        # diff_mw.SetTitle(mtitle)
        diff_m.Draw('e1')
        # diff_mw.Draw('e1 same')
        line.Draw()
        leg.Draw()
        
        cdiff.cd(2)
        diff_p.SetTitle(ptitle)
        # diff_pw.SetTitle(ptitle)
        diff_p.Draw('e1')
        # diff_pw.Draw('e1 same')
        line.Draw()
        
        cdiff.Print('mcasym_%(key)s_diffs%(stitle)s.png' % locals())
        
        opt2 = i>0 and 'e1 same' or 'e1'
        
        alldiffs.cd(1)
        line.Draw()
        diff_m.Draw(opt2)
        
        alldiffs.cd(2)
        line.Draw()
        diff_p.Draw(opt2)
        
        # alldiffsw.cd(1)
        # line.Draw()
        # diff_mw.Draw(opt2)
        
        # alldiffsw.cd(2)
        # line.Draw()
        # diff_pw.Draw(opt2)
        
        # diffs[key] = (diff_mw, diff_pw)
        
        keepme.extend([jp_m, jp_p, cdiff, leg])
    
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
    
    cmb.Print('mcasym_minbias%(stitle)s.png' % locals())
    # cmbw.Print('mcasym_minbias_reweight%(stitle)s.png' % locals())
    cjp.Print('mcasym_jetpatch%(stitle)s.png' % locals())
    alldiffs.Print('mcasym_diff%(stitle)s.png' % locals())
    # alldiffsw.Print('mcasym_diff_reweight%(stitle)s.png' % locals())
    
