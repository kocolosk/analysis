import ROOT
import math
# // $Id$

xsec = { 
'2_3'       : 8.150,
'3_4'       : 1.302,
'4_5'       : 3.158E-01,
'5_7'       : 1.372E-01,
'7_9'       : 2.290E-02,
'9_11'      : 5.495E-03,
'11_15'     : 2.220E-03,
'15_25'     : 3.907E-04,
'25_35'     : 1.074E-05,
'above_35'  : 5.300E-07,
'35_45'     : 5.000E-07,
'45_55'     : 2.857E-08,
'55_65'     : 1.451E-09
}

run6_samples = {
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
    out.Reset('ice')
    nbins = out.GetBin(out.GetNbinsX(), out.GetNbinsY(), out.GetNbinsZ())
    for h,nevents,sample in zip(histos, eventCounts, sampleIds):
        if sample == 'minbias': continue
        wp = xsec[sample] / nevents
        # print sample, nevents, wp
        for bin in range(1, nbins+1):
            content = out.GetBinContent(bin)
            error = out.GetBinError(bin)**2
            if h is None: 
                print 'wtf! why is this histo None?'
                continue
            nparticles = h.GetBinContent(bin)
            # if nparticles >= minCounts:
            content += nparticles * wp
            error   += wp * wp * nparticles * (1+nparticles/nevents)
            out.SetBinContent(bin,content)
            out.SetBinError(bin, math.sqrt(error))
    return out


def mergeSamples(outName, inputFileNames):
    outFile = ROOT.TFile(outName, 'recreate')
    inputFiles = [ROOT.TFile(n) for n in inputFileNames]
    eventCounts = [f.Get('eventCounter').GetEntries() for f in inputFiles]
    ## this will throw an exception on run5 samples
    sampleIds = [run6_samples[n[:7]] for n in inputFileNames]
    
    keys = inputFiles[0].GetListOfKeys()
    for key in keys:
        print key
        if key.GetName() == 'eventCounter': continue
        histos = [f.Get(key.GetName()) for f in inputFiles]
        
        # skip 2D histos for now b/c they are SLOW
        if histos[0].ClassName().startswith('TH2'): 
            [h.Delete() for h in histos]
            continue
        
        out = mergeHistos(histos, eventCounts, sampleIds)
        [h.Delete() for h in histos]
        outFile.cd()
        out.Write()
    
    outFile.Close()
    
def partonicCrossSection(sample='2_3', nevents=1000, sqrts=200):
    """runs standalone Pythia to determine xsec for weighting purposes, returns xsec"""
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


# /*****************************************************************************
#  * $Log$
#  *****************************************************************************/
