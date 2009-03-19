import ROOT
import math
import os.path

xsec = { 
'2_3'       : 8.150E+06,
'3_4'       : 1.302E+06,
'4_5'       : 3.158E+05,
'5_7'       : 1.372E+05,
'7_9'       : 2.290E+04,
'9_11'      : 5.495E+03,
'11_15'     : 2.220E+03,
'15_25'     : 3.907E+02,
'25_35'     : 1.074E+01,
'above_35'  : 5.300E-01,
'35_45'     : 5.000E-01,
'45_55'     : 2.857E-02,
'55_65'     : 1.451E-03
}

samples = {
### 2005
'rcf1224'   : '5_7',
'rcf1225'   : '7_9',
'rcf1226'   : '9_11',
'rcf1227'   : '11_15',
'rcf1228'   : '15_25',
'rcf1229'   : '25_35',
'rcf1230'   : 'above_35',
'rcf1231'   : '2_3',
'rcf1232'   : '3_4',
'rcf1233'   : '4_5',
'rcf1235'   : 'minbias', ## ?
'rcf1270'   : '45_55',
'rcf1271'   : '55_65',
'rcf1273'   : '0_2',
### 2006
'rcf1302'   : '45_55',
'rcf1303'   : '35_45',
'rcf1304'   : '55_65',
'rcf1306'   : '25_35',
'rcf1307'   : '15_25',
'rcf1308'   : '11_15',
'rcf1309'   : '9_11',
'rcf1310'   : '7_9',
'rcf1311'   : '5_7',
'rcf1317'   : '4_5',
'rcf1318'   : '3_4',
'rcf1319'   : 'minbias'
}

def mergeHistos(histos, eventCounts, sampleIds):
    out = histos[0].Clone()
    out.SetDirectory(None)
    out.Reset('ice')
    isProfile = out.ClassName().startswith('TProfile')
    nbins = out.GetBin(out.GetNbinsX(), out.GetNbinsY(), out.GetNbinsZ())
    for h,nevents,sample in zip(histos, eventCounts, sampleIds):
        if sample == 'minbias': continue
        wp = xsec[sample] / nevents
        if isProfile:
            out.Add(h, wp)
        else:
            # print sample, nevents, wp
            for bin in range(1, nbins+1):
                content = out.GetBinContent(bin)
                error = out.GetBinError(bin)**2
                nparticles = h.GetBinContent(bin)
                content += nparticles * wp
                error   += wp * wp * nparticles * (1+nparticles/nevents)
                out.SetBinContent(bin,content)
                out.SetBinError(bin, math.sqrt(error))
    return out

def mergeSamples(outName, inputFileNames):
    outFile = ROOT.TFile(outName, 'recreate')
    hlist = merge_samples(inputFileNames)
    outFile.cd()
    [h.Write() for h in hlist]
    outFile.Close()

def merge_samples(inputFileNames):
    merged = []
    inputFiles = [ROOT.TFile(n) for n in inputFileNames]
    eventCounts = [f.Get('eventCounter').GetEntries() for f in inputFiles]
    sampleIds = [samples[os.path.basename(n)[:7]] for n in inputFileNames]
    keys = inputFiles[0].GetListOfKeys()
    for key in keys:
        if key.GetName() == 'eventCounter': continue
        histos = [f.Get(key.GetName()) for f in inputFiles]
        print key
        merged.append( mergeHistos(histos, eventCounts, sampleIds) )
        [h.Delete() for h in histos]
    return merged

def mcasym(outName, inputFileNames, triggers=('jetpatch','117001'), keys=None):
    """
    caveats:
    select different nparticles histos for 05(pt) and 06(z_away2.Rebin())
    """
    outFile = ROOT.TFile(outName, 'recreate')
    inputFiles = [ROOT.TFile(n) for n in inputFileNames]
    oldkeys = keys or ['STD','MAX','MIN','ZERO','GS_NLOC','DSSV']
    keys = []
    [keys.extend([key, key+'w']) for key in oldkeys]
    for sub in ('anyspin', 'gg', 'qg', 'qq'):
        for trigger in triggers:
            for key in keys:
                minus_inputs = map(lambda f: \
                {
                    'id': samples[os.path.basename(f.GetName())[:7]],
                    'xsec': xsec[samples[os.path.basename(f.GetName())[:7]]],
                    'nevents': f.Get('eventCounter').GetEntries(),
                    'num': f.Get('_%s_%s_minus_%s' % (trigger, sub, key)),
                    'denom': f.Get('_%s_%s_minus_denom%s' % (trigger, sub,
                        key.endswith('w') and 'w' or '')),
                    # 'nparticles': f.Get('_%s_%s_minus_pt' % (trigger, sub))
                    'nparticles': f.Get('_%s_%s_minus_z_away2' % (trigger, sub))
                }, inputFiles)
                plus_inputs = map(lambda f: \
                {
                    'id': samples[os.path.basename(f.GetName())[:7]],
                    'xsec': xsec[samples[os.path.basename(f.GetName())[:7]]],
                    'nevents': f.Get('eventCounter').GetEntries(),
                    'num': f.Get('_%s_%s_plus_%s' % (trigger, sub, key)),
                    'denom': f.Get('_%s_%s_plus_denom%s' % (trigger, sub,
                        key.endswith('w') and 'w' or '')),
                    # 'nparticles': f.Get('_%s_%s_plus_pt' % (trigger, sub))
                    'nparticles': f.Get('_%s_%s_plus_z_away2' % (trigger, sub))
                }, inputFiles)
                outFile.cd()
                _mcasym_merge(minus_inputs).Write()
                _mcasym_merge(plus_inputs).Write()
    outFile.Close()
        
def _mcasym_merge(inputs, minParticlesToAccept=10):
    out = inputs[0]['num'].Clone()
    out.Reset('ice')
    
    ## don't include a bin from an individual sample in content or error
    ## if it has fewer than this minParticlesToAccept # of particles
    
    ## first do the bin contents
    bottom = []
    for bin in range(1, out.GetNbinsX()+1):
        top = 0.0
        bot = 0.0
        error = 0.0
        for sample in inputs:
            nparticles = sample['nparticles'].GetBinContent(bin)
            if nparticles < minParticlesToAccept: continue
            wp = sample['xsec']/sample['nevents']
            top += sample['num'].GetBinContent(bin) * wp
            bot += sample['denom'].GetBinContent(bin) * wp
        
        bottom.append(bot)
        content = bot > 0 and (top/bot) or 0.0
        out.SetBinContent(bin, content)
    
    ## now do the errors
    for bin in range(1, out.GetNbinsX()+1):
        error = 0.0
        
        for sample in inputs:
            nparticles = sample['nparticles'].GetBinContent(bin)
            if nparticles < minParticlesToAccept: continue
            wp = sample['xsec']/sample['nevents']
            nevents = sample['nevents']
            term0 = sample['num'].GetBinContent(bin)/nparticles
            term1 = term0 - out.GetBinContent(bin)
            term2 = sample['num'].GetBinError(bin)**2/nparticles - term0**2
            error += wp*wp*nparticles*(term1*term1*(1+nparticles/nevents)+term2)
        
        ferror = bottom[bin-1]>0 and math.sqrt(error/(bottom[bin-1]**2)) or 0.0
        out.SetBinError(bin, ferror)
        
    return out
    
    
def partonicCrossSection(sample='2_3', nevents=1000, sqrts=200):
    """runs standalone Pythia to determine xsec for weighting purposes.
    returns xsec"""
    import ROOT
    pythia = ROOT.TPythia6()
    
    ckMin = sample.split('_')[0]
    if ckMin == 'above': 
        ckMin = 35
    else:
        ckMin = float(ckMin)
    
    ckMax = sample.split('_')[1]
    if ckMax == '35':
        ckMax = 1000
    else:
        ckMax = float(ckMax)
    
    #pythia.SetMSTP(51, 8)
    #pythia.SetMSEL(1)
    
    # CDF Tune A for STAR
    pythia.SetMSEL(1) # could also be 2, but gives VERY different results
    pythia.SetMSTP(51, 7)
    pythia.SetMSTP(81, 1)
    pythia.SetMSTP(82, 4)
    pythia.SetPARP(67, 4.0)
    pythia.SetPARP(83, 0.5)
    pythia.SetPARP(84, 0.4)
    pythia.SetPARP(85, 0.9)
    pythia.SetPARP(86, 0.95)
    pythia.SetPARP(89, 1800)
    pythia.SetPARP(90, 0.25)
    pythia.SetPARP(91, 1.0)
    
    pythia.SetCKIN(3, ckMin)
    pythia.SetCKIN(4, ckMax)
    
    pythia.Initialize('CMS', 'p', 'p', sqrts)
    
    for i in range(nevents):
        if i % 1000 == 0: print 'generating event', i
        pythia.GenerateEvent()
        
    pythia.Pystat(1)
    
    return pythia.GetPARI(1)


def nlo_pion_fractions():
    simuFile = '/Users/kocolosk/data/run5-simu/hist/merged.cphist.root'
    
    import operator
    from array import array
    
    from analysis.plots import graphics
    from analysis.histos2 import HistogramManager
    
    bin_centers = array('d', [2.59, 3.87, 5.44, 7.56, 10.82])
    
    gg_gg = [4.439015E+06, 1.996755E+05, 1.263479E+04, 7.499839E+02, 2.810099E+01]
    qg_qg = [3.388824E+06, 1.976496E+05, 1.670780E+04, 1.382847E+03, 7.863113E+01]
    qq_qq = [2.527259E+05, 2.286667E+04, 2.908361E+03, 3.641042E+02, 3.269253E+01]
    qb_qb = [1.441496E+05, 1.107764E+04, 1.164664E+03, 1.152362E+02, 7.462816E+00]
    qb_gg = [2.679754E+04, 1.085830E+03, 7.030866E+01, 4.713206E+00, 2.204383E-01]
    
    total = [8.301383E+06, 4.332231E+05, 3.347704E+04, 2.609205E+03, 1.463661E+02]
    
    pythia_id = {
        11: qq_qq,
        12: qb_qb,
        13: qb_gg,
        28: qg_qg,
        68: gg_gg
    }
    
    
    graph = {
        11: ROOT.TGraph(5, bin_centers, array('d', 
            map(operator.div, pythia_id[11], total))),
        12: ROOT.TGraph(5, bin_centers, array('d', 
            map(operator.div, pythia_id[12], total))),
        # 13: ROOT.TGraph(5, bin_centers, array('d',
        #     map(operator.div, pythia_id[13], total))),
        28: ROOT.TGraph(5, bin_centers, array('d', 
            map(operator.div, pythia_id[28], total))),
        68: ROOT.TGraph(5, bin_centers, array('d', 
            map(operator.div, pythia_id[68], total)))
    }
    
    mgr = HistogramManager(ROOT.TFile(simuFile))
    hist = {
        11: mgr.qq['96011']['plus']['pt'],
        28: mgr.qg['96011']['plus']['pt'],
        68: mgr.gg['96011']['plus']['pt']
    }
    htotal = mgr.anyspin['96011']['plus']['pt']
    [h.Divide(htotal) for h in hist.values()]
    
    graph[68].SetLineColor(ROOT.kRed)
    graph[28].SetLineColor(ROOT.kGreen)
    graph[11].SetLineColor(ROOT.kBlue)
    graph[12].SetLineColor(ROOT.kMagenta)
    [g.SetLineWidth(2) for g in graph.values()]
    
    hist[68].SetMarkerColor(ROOT.kRed)
    hist[28].SetMarkerColor(ROOT.kGreen)
    hist[11].SetMarkerColor(ROOT.kBlue)
    [h.SetMarkerStyle(24) for h in hist.values()]
    
    bg = ROOT.TH2D('bg', '', 1, 2.00, 12.84, 1, 0., 0.6)
    bg.SetTitle('Subprocess Fractions (solid=NLO+DSS, points=Pythia)')
    bg.SetXTitle('#pi p_{T}')
    
    leg = ROOT.TLegend(0.8, 0.2, 0.88, 0.8)
    leg.AddEntry(graph[68], 'gg', 'l')
    leg.AddEntry(graph[28], 'qg', 'l')
    leg.AddEntry(graph[11], 'qq', 'l')
    leg.AddEntry(graph[12], 'q#bar{q}', 'l')
    
    c = graphics.canvas1()
    bg.Draw()
    [ g.Draw() for g in graph.values() ]
    [ h.Draw('same') for h in hist.values() ]
    leg.Draw()
    
    graphics.maybe_save()

