#!/usr/bin/env python
# encoding: utf-8
"""
myTest.py

Created by Adam Kocoloski on 2007-04-12.
Copyright (c) 2007 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import math
import sets
import getopt
import re
import paths, eventlist

import ROOT
#ROOT.gSystem.Load('StarSpinAnalyses')
#ROOT.gSystem.Load('StJetMaker')
#ROOT.gSystem.Load('StChargedPionAnalysisMaker')
#mystyle.use(1)
ROOT.gStyle.SetOptTitle(1)

help_message = '''
Here's how to use this script:

qa.py --update file1 file2 ...
creates a runnumber.qa.root file in $PWD for each input file specified

qa.py --
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class QAHistograms(list):
    QAHistogramTypes = ['pt','eta','phi','vz','nhitsfit','dedx','nsigmapion','dcaG']
    LogyValue = [1,0,0,0,0,0,0,0]
    
    def __init__(self,name,charge=None,trigIds=None,tfile=None):
        self.charge = charge
        self.trigIds = trigIds
        self.namePrefix = name
        self.names = ['%s_%s' % (name,elem) for elem in self.QAHistogramTypes]
        if tfile is not None:
            print 'reading from tfile'
            self.Read(tfile)
        else:
            self.append(ROOT.TH1D(self.names[0],'',10,2.,12.))
            self.append(ROOT.TH1D(self.names[1],'',60,-2.,2.))
            self.append(ROOT.TH1D(self.names[2],'',15,-math.pi,math.pi))
            self.append(ROOT.TH1D(self.names[3],'',100,-100.,100.))
            self.append(ROOT.TH1D(self.names[4],'',40,9.5,49.5))
            self.append(ROOT.TH1D(self.names[5],'',60,2.,5.))
            self.append(ROOT.TH1D(self.names[6],'',60,-3.,3.))
            self.append(ROOT.TH1D(self.names[7],'',50,0.,5.))
    
    def Fill(self,name,val):
        self[self.QAHistogramTypes.index(name)].Fill(val)
    
    def Write(self,name='',option=None):
        [elem.Write(name,option) for elem in self]
    
    def Read(self,tfile):
        self = [tfile.Get(elem) for elem in self.names]
        print self
        
    def Reset(self):
        [h.Reset() for h in self]
        
    def Print(self,title='',path='./',format='.eps'):
        c = ROOT.TCanvas()
        index = 0
        print self
        for elem in self:
            elem.SetXTitle(self.QAHistogramTypes[index])
            elem.SetTitle(title)
            elem.Draw()
            ROOT.gPad.SetLogy(self.LogyValue[index])
            c.SaveAs(path + elem.GetName() + format)
            index += 1


class QAHistogramsSummary(QAHistograms):
    SummaryYTitle = ['< p_{T} >', '< #eta >', '< #phi >', '< v_{z} >', '< nHitsFit >', 
                     '< dEdx >', '< nSigmaPion >', '< DCA_{global} >']
                     
    def __init__(self,name,numberOfRuns,year):
        self.currentRunIndex = 1
        self.numberOfRuns = numberOfRuns
        self.namePrefix = name
        self.names = ['%s_%s' % (name,elem) for elem in self.QAHistogramTypes]
        
        #self = []
        for i in xrange(len(self.QAHistogramTypes)):
            self.append(ROOT.TH1F(self.names[i],'',numberOfRuns,0.5,numberOfRuns+0.5))
            self[i].SetMarkerStyle(24)
            self[i].SetMarkerSize(0.7)
                    
        self.bg = []
        if year == 2006:
            self.bg.append(ROOT.TH2D(self.names[0] + '_bg','',1,0.5,numberOfRuns+0.5, 1,3.,4.))
            self.bg.append(ROOT.TH2D(self.names[1] + '_bg','',1,0.5,numberOfRuns+0.5, 1,-0.3,0.3))
            self.bg.append(ROOT.TH2D(self.names[2] + '_bg','',1,0.5,numberOfRuns+0.5, 1,-0.5,0.5))
            self.bg.append(ROOT.TH2D(self.names[3] + '_bg','',1,0.5,numberOfRuns+0.5, 1,-20.,20.))
            self.bg.append(ROOT.TH2D(self.names[4] + '_bg','',1,0.5,numberOfRuns+0.5, 1,30.,40.))
            self.bg.append(ROOT.TH2D(self.names[5] + '_bg','',1,0.5,numberOfRuns+0.5, 1,2.5,3.5))
            self.bg.append(ROOT.TH2D(self.names[6] + '_bg','',1,0.5,numberOfRuns+0.5, 1,-1.0,1.0))
            self.bg.append(ROOT.TH2D(self.names[7] + '_bg','',1,0.5,numberOfRuns+0.5, 1,0.,1.5))
        if year == 2005:
            self.bg.append(ROOT.TH2D(self.names[0] + '_bg','',1,0.5,numberOfRuns+0.5, 1,2.5,3.5))
            self.bg.append(ROOT.TH2D(self.names[1] + '_bg','',1,0.5,numberOfRuns+0.5, 1,-0.3,0.7))
            self.bg.append(ROOT.TH2D(self.names[2] + '_bg','',1,0.5,numberOfRuns+0.5, 1,-0.5,0.5))
            self.bg.append(ROOT.TH2D(self.names[3] + '_bg','',1,0.5,numberOfRuns+0.5, 1,-30.,10.))
            self.bg.append(ROOT.TH2D(self.names[4] + '_bg','',1,0.5,numberOfRuns+0.5, 1,30.,40.))
            self.bg.append(ROOT.TH2D(self.names[5] + '_bg','',1,0.5,numberOfRuns+0.5, 1,2.5,3.5))
            self.bg.append(ROOT.TH2D(self.names[6] + '_bg','',1,0.5,numberOfRuns+0.5, 1,-1.0,1.0))
            self.bg.append(ROOT.TH2D(self.names[7] + '_bg','',1,0.5,numberOfRuns+0.5, 1,0.,1.5))
        
        for i in xrange(len(self.bg)):
            self.bg[i].SetXTitle('run index')
            self.bg[i].SetYTitle(self.SummaryYTitle[i])
    
    def AddHistogramsFromFile(self,inputFileName):
        inputFile = ROOT.TFile(inputFileName,'read')
        for elem in self:
            tmp = inputFile.Get(elem.GetName())
            if tmp is not None and tmp.GetEntries() > 0:
                error = 1.0/math.sqrt(tmp.GetEntries())
                elem.SetBinContent(self.currentRunIndex,tmp.GetMean())
                elem.SetBinError(self.currentRunIndex,error)
                tmp.Delete()
        self.currentRunIndex = self.currentRunIndex + 1
        inputFile.Close()
    
    def Print(self,title='',path='./',format='.eps'):
        c = ROOT.TCanvas()
        for i in xrange(len(self)):
            self.bg[i].SetTitle(self.namePrefix)
            self.bg[i].Draw()
            self[i].Draw('same')
            c.SaveAs(path + self[i].GetName() + format)
    


def printMap(mymap):
    #example loop over STL map  
    counter = 0
    try:
        myiter = iter(mymap)
        while True:
            mypair = myiter.next() #note mypair is *not* an iterator, more like *(it++)
            print 'mapsize=%d, imap=%d, (id,adc)=(%d,%d)' % (mymap.size(),counter,mypair.first,mypair.second)
            counter = counter + 1
    except StopIteration: pass


def initQaList(run):
    if run == 2006:
        return [    QAHistograms('bjp1_plus',     charge=(1,),    trigIds=(137221,137222)     ),
                    QAHistograms('bjp1_minus',    charge=(-1,),   trigIds=(137221,137222)     ),
                    QAHistograms('bjp1_sum',      charge=(1,-1),  trigIds=(137221,137222)     ),
                    QAHistograms('l2jetB_plus',   charge=(1,),    trigIds=(137622,)           ),
                    QAHistograms('l2jetB_minus',  charge=(-1,),   trigIds=(137622,)           ),
                    QAHistograms('l2jetB_sum',    charge=(1,-1),  trigIds=(137622,)           ),
                    QAHistograms('l2gammaB_plus', charge=(1,),    trigIds=(5,137611)          ),
                    QAHistograms('l2gammaB_minus',charge=(-1,),   trigIds=(5,137611)          ),
                    QAHistograms('l2gammaB_sum',  charge=(1,-1),  trigIds=(5,137611)          ) 
        ]
    if run == 2005:
        return [    QAHistograms('mb_plus',     charge=(1,),    trigIds=(96011,)),
                    QAHistograms('mb_minus',    charge=(-1,),   trigIds=(96011,)),
                    QAHistograms('mb_sum',      charge=(1,-1),  trigIds=(96011,)),
                    QAHistograms('bht1_plus',   charge=(1,),    trigIds=(96201,)),
                    QAHistograms('bht1_minus',  charge=(-1,),   trigIds=(96201,)),
                    QAHistograms('bht1_sum',    charge=(1,-1),  trigIds=(96201,)),
                    QAHistograms('bht2_plus',   charge=(1,),    trigIds=(96211,)),
                    QAHistograms('bht2_minus',  charge=(-1,),   trigIds=(96211,)),
                    QAHistograms('bht2_sum',    charge=(1,-1),  trigIds=(96211,)),
                    QAHistograms('bjp1_plus',   charge=(1,),    trigIds=(96221,)),
                    QAHistograms('bjp1_minus',  charge=(-1,),   trigIds=(96221,)),
                    QAHistograms('bjp1_sum',    charge=(1,-1),  trigIds=(96221,)),
                    QAHistograms('bjp2_plus',   charge=(1,),    trigIds=(96233,)),
                    QAHistograms('bjp2_minus',  charge=(-1,),   trigIds=(96233,)),
                    QAHistograms('bjp2_sum',    charge=(1,-1),  trigIds=(96233,))
        ]
    return None


def initQaSummaryList(run,length):
    if run == 2006:
        return [    QAHistogramsSummary('bjp1_plus',     length, run),
                    QAHistogramsSummary('bjp1_minus',    length, run),
                    QAHistogramsSummary('bjp1_sum',      length, run),
                    QAHistogramsSummary('l2jetB_plus',   length, run),
                    QAHistogramsSummary('l2jetB_minus',  length, run),
                    QAHistogramsSummary('l2jetB_sum',    length, run),
                    QAHistogramsSummary('l2gammaB_plus', length, run),
                    QAHistogramsSummary('l2gammaB_minus',length, run),
                    QAHistogramsSummary('l2gammaB_sum',  length, run) 
        ]
    if run == 2005:
        return [    QAHistogramsSummary('mb_plus',     length, run),
                    QAHistogramsSummary('mb_minus',    length, run),
                    QAHistogramsSummary('mb_sum',      length, run),
                    QAHistogramsSummary('bht1_plus',   length, run),
                    QAHistogramsSummary('bht1_minus',  length, run),
                    QAHistogramsSummary('bht1_sum',    length, run),
                    QAHistogramsSummary('bht2_plus',   length, run),
                    QAHistogramsSummary('bht2_minus',  length, run),
                    QAHistogramsSummary('bht2_sum',    length, run),
                    QAHistogramsSummary('bjp1_plus',   length, run),
                    QAHistogramsSummary('bjp1_minus',  length, run),
                    QAHistogramsSummary('bjp1_sum',    length, run),
                    QAHistogramsSummary('bjp2_plus',   length, run),
                    QAHistogramsSummary('bjp2_minus',  length, run),
                    QAHistogramsSummary('bjp2_sum',    length, run)
        ]
    return None


def updateHistograms(runnumber,qalist):
    [qa.Reset() for qa in qalist]
    
    year = 2005
    if runnumber > 7000000: year = 2006
    
    #open jetSkim file
    jetSkimFileName = '%s%d.tree.root' % (paths.path(site='MIT',fileType='jetSkim',run=year),runnumber)
    jetSkimFile = ROOT.TFile(jetSkimFileName)
    jetSkimTree = jetSkimFile.Get('jetSkimTree')
    ev = ROOT.StJetSkimEvent()
    jetSkimTree.SetBranchAddress('skimEventBranch',ev)
    
    #create event list -- really just the trigger selection in case of QA
    eventList = eventlist.eventlist()
    if year == 2006:    eventList.addTriggers([5, 137611, 137622, 137221, 137222])
    else:               eventList.addTriggers([96011, 96201, 96211, 96221, 96233])
    eventList.requireDidFire = True
    jetSkimTree.Draw('>>elist',str(eventList))
    elist = ROOT.gDirectory.Get('elist')
    jetSkimTree.SetEventList(elist)
    
    #open chargedPion file and create index
    chargedPionFileName = '%s%d.tree.root' % (paths.path(site='MIT',fileType='chargedPions',run=year),runnumber)
    chargedPionFile = ROOT.TFile(chargedPionFileName)
    chargedPionTree = chargedPionFile.Get('chargedPionTree')
    tracks  = ROOT.TClonesArray('StChargedPionTrack',100)
    try:
        chargedPionTree.SetBranchAddress('primaries',tracks)
    except AttributeError: return
    chargedPionTree.BuildIndex('run','event')
            
    for i in xrange(elist.GetN()):
        jetSkimTree.GetEntry(i)
        chargedPionTree.GetEntryWithIndex(ev.runId(),ev.eventId())
        
        #make a set of all fired triggers
        triggersDidFire = [elem.trigId() for elem in list(ev.triggers()) if elem.didFire() > 0]
        triggersDidFire = sets.Set(triggersDidFire)
                
        #event level quantities == vz
        if ev.bestVert() is not None:
            vz = ev.bestVert().position()[2]
            [qa.Fill('vz',vz) for qa in qalist if len(triggersDidFire & sets.Set(qa.trigIds)) > 0]
        
        #now loop through tracks
        for track in list(tracks):
            if track.nHitsFit() < 11: continue
            q = int(track.charge())            
            for qa in qalist:
                if q in qa.charge and len(triggersDidFire & sets.Set(qa.trigIds)) > 0:                    
                    qa.Fill('pt',track.pt())
                    qa.Fill('eta',track.eta())
                    qa.Fill('phi',track.phi())
                    qa.Fill('nhitsfit',track.nHitsFit())
                    qa.Fill('dedx',track.dEdx()*10**6)
                    qa.Fill('nsigmapion',track.nSigmaPion())
                    qa.Fill('dcaG',track.globalDca().mag())
    
    #ok, that's it -- time to go
    qaFileName = '%d.qa.root' % (runnumber,)
    qaFile = ROOT.TFile(qaFileName,'recreate')
    [qa.Write('',ROOT.TObject.kOverwrite) for qa in qalist]
    qaFile.Close()
    print 'wrote file ' + qaFileName


def updateHistogramsInFile(inputFileName):
    #open the file and load the tree
    inputFile = ROOT.TFile(inputFileName,'update')
    inputTree = ROOT.gDirectory.Get('spinAnalysisTree')
    
    #create the objects needed for the branch addresses
    theEvent    = ROOT.StJetSkimEvent()
    primTracks  = ROOT.TClonesArray('StChargedPionTrack',100)
    globTracks  = ROOT.TClonesArray('StChargedPionTrack',100)
    
    #set the branch addresses
    inputTree.SetBranchAddress('skimEventBranch',theEvent)
    inputTree.SetBranchAddress('chargedPionPrimaries',primTracks)
    inputTree.SetBranchAddress('chargedPionGlobals',globTracks)
    
    #histograms
    qalist = []
    qalist.append(QAHistograms('bjp1_plus',     charge=(1,),    trigIds=(137221,137222)     ))
    qalist.append(QAHistograms('bjp1_minus',    charge=(-1,),   trigIds=(137221,137222)     ))
    qalist.append(QAHistograms('bjp1_sum',      charge=(1,-1),  trigIds=(137221,137222)     ))
    qalist.append(QAHistograms('l2jetB_plus',   charge=(1,),    trigIds=(137622,)           ))
    qalist.append(QAHistograms('l2jetB_minus',  charge=(-1,),   trigIds=(137622,)           ))
    qalist.append(QAHistograms('l2jetB_sum',    charge=(1,-1),  trigIds=(137622,)           ))
    qalist.append(QAHistograms('l2gammaB_plus', charge=(1,),    trigIds=(5,137611)          ))
    qalist.append(QAHistograms('l2gammaB_minus',charge=(-1,),   trigIds=(5,137611)          ))
    qalist.append(QAHistograms('l2gammaB_sum',  charge=(1,-1),  trigIds=(5,137611)          ))
    
    entries = inputTree.GetEntries()
    print 'update QA histograms using %6d spinAnalysisTree entries in file %s' % (entries,inputFileName)
    for i in xrange(entries):
        inputTree.GetEntry(i)
        
        #make a set of all fired triggers
        triggersDidFire = []
        for k in xrange(theEvent.triggers().GetEntries()):
            trig = theEvent.triggers().At(k)
            if trig.didFire() > 0: triggersDidFire.append(trig.trigId())
        triggersDidFire = sets.Set(triggersDidFire)
        
        #event level quantities
        vz = theEvent.bestVert().position()[2]
        for elem in qalist:
            triggersOfInterest = sets.Set(elem.trigIds)
            if len(triggersOfInterest & triggersDidFire) > 0: 
                elem.Fill('vz',vz)
        
        #now loop through tracks
        for j in xrange(primTracks.GetEntries()):
            track = primTracks.At(j)
            globalTrack = globTracks.At(j)
            
            for elem in qalist:
                q = int(track.charge())
                triggersOfInterest = sets.Set(elem.trigIds)
                if q in elem.charge and len(triggersOfInterest & triggersDidFire) > 0:
                    elem.Fill('pt',track.pt())
                    elem.Fill('eta',track.eta())
                    elem.Fill('phi',track.phi())
                    elem.Fill('nhits',track.nHits())
                    elem.Fill('nhitsratio',float(track.nHits())/track.nHitsPoss())
                    elem.Fill('nhitsfit',track.nHitsFit())
                    elem.Fill('nhitsfitratio',float(track.nHitsFit())/track.nHitsPoss())
                    elem.Fill('nsigmapion',track.nSigmaPion())
                    elem.Fill('dcaG',globalTrack.dca().mag())
    
    if inputFile.cd('qa') is 0:
        inputFile.mkdir('qa')
        inputFile.cd('qa')
    [elem.Write('',ROOT.TObject.kOverwrite) for elem in qalist]
    inputFile.Close()


def printHistogramsInFile(inputFileName,path='./',format='.eps'):
    inputFile = ROOT.TFile(inputFileName,'read')
    
    runnumber = inputFileName.split('.')[0]
    path = path + runnumber + '/'
    
    try:
        os.mkdir(path)
    except OSError: pass
    
    qalist = []
    qalist.append(QAHistograms('bjp1_plus',     tfile=inputFile))
    qalist.append(QAHistograms('bjp1_minus',    tfile=inputFile))
    qalist.append(QAHistograms('bjp1_sum',      tfile=inputFile))
    qalist.append(QAHistograms('l2jetB_plus',   tfile=inputFile))
    qalist.append(QAHistograms('l2jetB_minus',  tfile=inputFile))
    qalist.append(QAHistograms('l2jetB_sum',    tfile=inputFile))
    qalist.append(QAHistograms('l2gammaB_plus', tfile=inputFile))
    qalist.append(QAHistograms('l2gammaB_minus',tfile=inputFile))
    qalist.append(QAHistograms('l2gammaB_sum',  tfile=inputFile))
    
    [elem.Print(runnumber,path,format) for elem in qalist]
    inputFile.Close()


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "bhups", ["batch", "help", 'update', 'print','summary'])
        except getopt.error, msg:
            raise Usage(msg)
    
        ROOT.gStyle.SetCanvasColor(10)
        ROOT.gStyle.SetCanvasBorderMode(0)
    
        #magic!
        throwitaway = ROOT.std.map(long,float)()
        
        # option processing
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
            if option in ('-u','--update','-s','--summary'):
                runs = [re.search('(6|7)\d{6}',inputPath).group() for inputPath in args]
                if int(runs[0]) > 7000000:  year = 2006
                else:                       year = 2005
                for i in range(len(runs)):
                    print '%3d     %s' % (i+1, runs[i])
            if option in ('-u','--update'):
                qalist = initQaList(year)
                [updateHistograms(int(run),qalist) for run in runs]
            if option in ('-p','--print'):
                [printHistogramsInFile(elem,'./','.gif') for elem in args]
            if option in ('-s','--summary'):
                ROOT.gStyle.SetOptStat(0)
                sumlist = initQaSummaryList(year, len(args))
                for elem in sumlist:
                    [elem.AddHistogramsFromFile(inputFile) for inputFile in args]
                    print 'printing QAHistogramsSummary for %s' % (elem.namePrefix,)
                    path = "./%s/" % (elem.namePrefix,)
                    try:
                        os.mkdir(path)
                    except OSError: pass
                    elem.Print('',path,'.gif')
                
                
        return 0
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
        
if __name__ == '__main__':
    sys.exit(main())
