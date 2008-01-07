#!/usr/bin/env python
# encoding: utf-8
"""
maintainSpinTrees.py

Created by Adam Kocoloski on 2007-05-09.
Copyright (c) 2007 __MyCompanyName__. All rights reserved.
"""
import os, sys, getopt, re, math, time
from array import array
from glob import glob
print 'pyroot -- import ROOT'
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kFatal

help_message = '''
This is a set of scripts used to manage the creation and updating of common Spin PWG analysis trees.  Usage:  

maintainSpinTrees.py -c jetSkimFile1 jetSkimFile2 ...
creates new spinTree files in the current working directory starting from the jetSkimEvent files.  No particle branches are added

maintainSpinTrees.py -u { ConeJets | ConeJets12 | ConeJetsEMC | chargedPions | bemcPions | bemcElectrons | gammas } file1 file2 ...
copies each file into $PWD/new and adds / updates the specified branch.
'''

allPaths = { '2005':{}, '2006':{} }

allPaths['2005']['MIT'] = {  'jets':'/Volumes/scratch/common/run5/jets/jets_',
                    'jetSkim':'/Volumes/scratch/common/run5/jetSkim/jetSkim_',
                    'chargedPions':'/Volumes/scratch/common/run5/chargedPions/chargedPions_',
                    'bemcPions':'/Volumes/scratch/common/run5/bemcPions/NeutralPionTree.root'
}

allPaths['2005']['RCF'] = {  'jets':'/star/institutions/mit/common/run5/jets/jets_',
                    'jetSkim':'/star/institutions/mit/common/run5/jetSkim/jetSkim_',
                    'chargedPions':'/star/institutions/mit/common/run5/chargedPions/chargedPions_'
}

allPaths['2005']['PDSF'] = { 'jets':'/dante3/starspin/msar/2005/jets/jets_',
                    'jetSkim':'/home/kocolosk/analysis/run5/may04/jetSkim_merged/jetSkim_',
                    'chargedPions':'/home/kocolosk/analysis/run5/may04/tracks_merged/chargedPions_'
}

allPaths['2006']['MIT'] = {  'jets':'/Volumes/scratch/common/run6/jets/jets_',
                    'jetSkim':'/Volumes/scratch/common/run6/jetSkim/jetSkim_',
                    'chargedPions':'/Volumes/scratch/common/run6/chargedPions/chargedPions_',
                    'bemcPions':'/Volumes/scratch/common/run6/bemcPions/NeutralPionTree.root',
                    #'bemcElectrons':'/Volumes/scratch/common/run6/bemcElectrons/',
                    'bemcElectrons':'root://deltag5.lns.mit.edu//Volumes/data01/run6/bemcElectrons/',
                    'gammaCandidates':'/Volumes/scratch/common/run6/gammas/gamma_'
}

allPaths['2006']['Xgrid'] = {  'jets':'/Volumes/deltag5.lns.mit.edu/common/run6/jets/jets_',
                    'jetSkim':'/Volumes/deltag5.lns.mit.edu/common/run6/jetSkim/jetSkim_',
                    'chargedPions':'/Volumes/deltag5.lns.mit.edu/common/run6/chargedPions/chargedPions_',
                    'bemcPions':'/Volumes/deltag5.lns.mit.edu/common/run6/bemcPions/NeutralPionTree.root',
                    'bemcElectrons':'/Volumes/deltag5.lns.mit.edu/common/run6/bemcElectrons/',
                    'gammaCandidates':'/Volumes/deltag5.lns.mit.edu/common/run6/gammas/gamma_'
}

allPaths['2006']['RCF'] = {  'jets':'/star/institutions/mit/common/run6/jets/jets_',
                    'jetSkim':'/star/institutions/mit/common/run6/jetSkim/jetSkim_',
                    'chargedPions':'/star/institutions/mit/common/run6/chargedPions/chargedPions_',
                    'bemcPions':'/star/institutions/mit/common/run6/bemcPions/NeutralPionTree.root',
                    'gammaCandidates':'/star/institutions/mit/betan/ROOT_Files/Gamma/gamma_'
}

allPaths['2006']['PDSF'] = { 'jets':'/dante3/starspin/msar/2006/TTrees/jets_',
                    'jetSkim':'/dante3/starspin/msar/2006/TTrees/skimTrees.new/jetSkim_',
                    'chargedPions':'/???/'
}

year = '2006'
site = 'MIT'
paths = allPaths[year][site]


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def acceptEvent(ev):
    if ev.isValid() == 0:               return False
    if ev.isPolLong() == 0:             return False
    if ev.isMaskedUsingBx48() != 0:     return False
    if ev.offsetBx48minusBX7() != 0:    return False
    
    #require a found vertex
    #if ev.vertices().GetEntries() == 0: return False
    
    #trigger selection:  production + 5
    nTrigs = ev.triggers().GetEntries()
    for i in xrange(nTrigs):
        currentTrig = ev.triggers().At(i)
        if currentTrig.trigId() < 1000 and currentTrig.trigId() != 5: continue
        return True
    
    return False


def acceptJet(jet,bestVert):
    if bestVert is None:                return False
    
    if year == '2005':
        if jet.detEta(bestVert.position()[2]) < 0.2:    return False
        if jet.detEta(bestVert.position()[2]) > 0.8:    return False
    
    if year == '2006':
        pass
    
    return True


def acceptChargedPion(pion):
    if math.fabs(pion.eta()) > 1:       return False
    if pion.nHitsFit() < 25:            return False
    
    if pion.globalDca().mag() > 1:      return False
    
    if pion.nSigmaPion() < -1:          return False
    if pion.nSigmaPion() > 2:           return False
    
    return True


def acceptBemcPion(pion):
    return True


def acceptGammaCandidate(gamma):
    return True


def getRunNumber(inputPath):
    regex = re.search('(6|7)\d{6}',inputPath)
    return regex.group()



def loadLibraries():
    print 'pyroot -- loading libraries'
    ROOT.gSystem.Load('libPhysics')
    ROOT.gSystem.Load('libTable')
    ROOT.gSystem.Load('StarRoot')
    ROOT.gSystem.Load('St_base')
    ROOT.gSystem.Load('StarClassLibrary')
    ROOT.gSystem.Load('StChain')
    ROOT.gSystem.Load('St_Tables')
    ROOT.gSystem.Load('StDetectorDbMaker')
    ROOT.gSystem.Load('StEvent')
    ROOT.gSystem.Load('StEmcUtil')
    ROOT.gSystem.Load('StStrangeMuDstMaker')
    ROOT.gSystem.Load('StMuDSTMaker')
    ROOT.gSystem.Load('StSpinDbMaker')
    ROOT.gSystem.Load('StEmcTriggerMaker')
    ROOT.gSystem.Load('StMCAsymMaker')
    ROOT.gSystem.Load('StJetFinder')
    
    ROOT.gSystem.Load('StJetMaker')
    
    ROOT.gSystem.Load('StMcEvent')
    ROOT.gSystem.Load('StChargedPionAnalysisMaker')
    
    ROOT.gSystem.Load('StEEmcUtil')
    ROOT.gSystem.Load('StEEmcA2EMaker')
    ROOT.gSystem.Load('StGammaMaker')
    
    ROOT.gSystem.Load('StSpinTree')
    print 'pyroot -- loading libraries complete'


def loadLibrariesMIT():
    ROOT.gSystem.Load('StarSpinAnalyses')
    ROOT.gSystem.Load('StJets')
    ROOT.gSystem.Load('StEmcUtilities')


def loadLibrariesXgrid():
    ROOT.gSystem.Load('/Volumes/star1.lns.mit.edu/lib/StarSpinAnalyses.so')
    ROOT.gSystem.Load('/Volumes/star1.lns.mit.edu/lib/StJets.so')
    ROOT.gSystem.Load('/Volumes/star1.lns.mit.edu/lib/StEmcUtilities.so')


def fixSkimTree(inputPath, outDir):    
    jetSkimFile = ROOT.TFile.Open(inputPath, 'read')
    jetSkimTree = jetSkimFile.Get('jetSkimTree')
    oldEvent = ROOT.StJetSkimEvent()
    jetSkimTree.SetBranchAddress('skimEventBranch', oldEvent)
    
    runnumber = re.search('(6|7)\d{6}',inputPath).group()
    outputName = '%sjetSkim_%s.tree.root' % (outDir,runnumber)
    outputFile = ROOT.TFile(outputName, 'recreate')
    newEvent = ROOT.StJetSkimEvent()
    ROOT.StJetSkimTrig.Class().IgnoreTObjectStreamer()
    ROOT.StPythiaEvent.Class().IgnoreTObjectStreamer()
    outputTree = ROOT.TTree('jetSkimTree','StJetSkimEvent tree')
    outputTree.Branch ("skimEventBranch", "StJetSkimEvent", newEvent, 64000, 99);
    
    newHeaders = ROOT.TClonesArray('StJetSkimTrigHeader',50)
    oldHeaders = list(jetSkimTree.GetUserInfo().At(0))
    newList = [h for h in oldHeaders]
    for j in range(len(newList)):
        newHeaders[j] = newList[j]
    outputTree.GetUserInfo().Add(newHeaders)
    outputTree.GetUserInfo().At(0).SetUniqueID(jetSkimTree.GetUserInfo().At(0).GetUniqueID())
        
    for i in xrange(jetSkimTree.GetEntries()):
        jetSkimTree.GetEntry(i)
        
        newEvent.setFill(oldEvent.fill())
        newEvent.setRunId(oldEvent.runId())
        newEvent.setEventId(oldEvent.eventId())
        newEvent.setMudstFileName(oldEvent.mudstFileName())
        
        newEvent.setBx7(oldEvent.bx7())
        newEvent.setBx48(oldEvent.bx48())
        newEvent.setSpinBits(oldEvent.spinBits())
        newEvent.setEbbc(oldEvent.eBbc())
        newEvent.setWbbc(oldEvent.wBbc())
        newEvent.setBbcTimeBin(oldEvent.bbcTimeBin())
        
        newEvent.setIsValid(oldEvent.isValid())
        newEvent.setIsPolLong(oldEvent.isPolLong())
        newEvent.setIsPolTrans(oldEvent.isPolTrans())
        newEvent.setIsMaskedUsingBx48(oldEvent.isMaskedUsingBx48())
        newEvent.setOffsetBx48minusBX7(oldEvent.offsetBx48minusBX7())
        newEvent.setSpin4UsingBx48(oldEvent.spin4usingBx48())
        
        newEvent.setL2Result(oldEvent.l2Result())
        
        newEvent.setMcEvent(oldEvent.mcEvent())
        
        for trig in oldEvent.triggers():
            newEvent.setTrig(trig)
            
        for vert in oldEvent.vertices():
            newEvent.setVert(vert)
        
        newEvent.setBestVert(0)
        
        newEvent.setTrigHeaderArray(outputTree.GetUserInfo().At(0))
        
        outputTree.Fill()
        newEvent.Clear()
    
    jetSkimFile.cd()
    jetSkimFile.Close()    
    
    outputFile.cd()
    outputTree.Write()
    outputFile.Close()
    
    
    print 'wrote', outputName


def makeSpinTree(inputPath,outDir):
    jetSkimFile = ROOT.TFile.Open(inputPath,'read')
    jetSkimTree = jetSkimFile.Get('jetSkimTree')
    event = ROOT.StJetSkimEvent()
    jetSkimTree.SetBranchAddress('skimEventBranch',event)
    
    runnumber = re.search('(6|7)\d{6}',inputPath).group()
    
    outputName = '%sspinAnalyses_%s.tree.root' % (outDir,runnumber)
    outputFile = ROOT.TFile(outputName,'recreate')
    outputTree = jetSkimTree.CloneTree(0)
    
    #slightly ugly hack to keep TRefs working for trigHeaders
    headers = outputTree.GetUserInfo()
    for i in xrange(headers.GetEntries()):
        headers.At(i).SetUniqueID(jetSkimTree.GetUserInfo().At(i).GetUniqueID())
    
    oldentries = jetSkimTree.GetEntries()
    for i in xrange(oldentries):
        jetSkimTree.GetEntry(i)
        if acceptEvent(event) == True: outputTree.Fill()
        
    jetSkimFile.Close()
    
    outputTree.SetName('spinTree')
    outputTree.SetTitle('updated %s' % (time.ctime(),))
    
    outputTree.Write()
    outputFile.Close()
    
    print 'built a spin tree at %s' % (outputName,)


def updateChargedPions(inputPath):
    runnumber = getRunNumber(inputPath)
    
    spinFile = ROOT.TFile.Open(inputPath,'read')
    if spinFile is None: raise Usage('could not open spin file at %s' % (inputPath,))
    spinTree = spinFile.Get('spinTree')
    if spinTree is None: raise Usage('could not get spinTree in file at %s' % (inputPath,))
    
    event = ROOT.StJetSkimEvent()
    spinTree.SetBranchAddress('skimEventBranch',event)
    
    #create the output file now so friendTree is disk-resident
    try: os.mkdir('new')
    except OSError: pass
    outPath = 'new/spinAnalyses_%s.tree.root' % (runnumber,)
    outFile = ROOT.TFile(outPath,'recreate')
    copySpinTree(outFile,spinTree)
    
    friendTree = ROOT.TTree('chargedPions','updated %s' % (time.ctime(),))
    goodTracks = ROOT.TClonesArray('StChargedPionTrack',100)
    nParticles = array('b',[-1])
    friendTree.Branch('chargedPions',goodTracks)
    friendTree.Branch('nChargedPions',nParticles,'nChargedPions/B')
    
    isRealFile = True
    
    inputName = '%s%s.tree.root' % (paths['chargedPions'],runnumber)
    inputFile = ROOT.TFile.Open(inputName,'read')
    if inputFile is None or inputFile.IsOpen() == False: isRealFile = False
    
    if isRealFile:
        inputTree = inputFile.Get('chargedPionTree')
        if inputTree is None: raise Usage('could not connect tree in file %s' % (inputName,))
        oldTracks = ROOT.TClonesArray('StChargedPionTrack',100)
        inputTree.SetBranchAddress('primaries',oldTracks)
        inputTree.BuildIndex('run','event')
    
    for i in range(spinTree.GetEntries()):
        spinTree.GetEntry(i)
        if isRealFile:
            nBytesRead = inputTree.GetEntryWithIndex(event.runId(),event.eventId())
            if nBytesRead > 0:
                oldList = list(oldTracks)
                goodList = [track for track in oldList if acceptChargedPion(track)]
                nParticles[0] = len(goodList)
                for j in range(len(goodList)): goodTracks[j] = goodList[j]
        friendTree.Fill()
        goodTracks.Clear()
        nParticles[0] = -1
    
    inputFile.Close()
    
    saveSpinTree(outFile,spinTree,friendTree,[spinFile.Get('bemcPions'),
                                              spinFile.Get('ConeJets'),
                                              spinFile.Get('ConeJetsEMC'),
                                              spinFile.Get('bemcElectrons'),
                                              spinFile.Get('gammaCandidates')])
    
    spinFile.Close()
    
    print 'updated the chargedPions tree at %s' % (outPath,)


def updateBemcPions(inputPath):
    runnumber = getRunNumber(inputPath)
    
    spinFile = ROOT.TFile.Open(inputPath,'read')
    if spinFile is None: raise Usage('could not open spin file at %s' % (inputPath,))
    
    spinTree = spinFile.Get('spinTree')
    if spinTree is None: raise Usage('could not get spinTree in file at %s' % (inputPath,))
    
    event = ROOT.StJetSkimEvent()
    spinTree.SetBranchAddress('skimEventBranch',event)
    
    #create the output file now so friendTree is disk-resident
    try: os.mkdir('new')
    except OSError: pass
    outPath = 'new/spinAnalyses_%s.tree.root' % (runnumber,)
    outFile = ROOT.TFile(outPath,'recreate')
    copySpinTree(outFile,spinTree)
    
    friendTree = ROOT.TTree('bemcPions','updated %s' % (time.ctime(),))
    goodTracks = ROOT.TClonesArray('TPi0Candidate',100)
    nParticles = array('b',[-1])
    friendTree.Branch('bemcPions',goodTracks)
    friendTree.Branch('nBemcPions',nParticles,'nBemcPions/B')
    
    inputFile = ROOT.TFile.Open(paths['bemcPions'],'read')
    if inputFile is None: raise Usage('could not open file at %s' % (paths['bemcPions'],))
        
    inputTree = inputFile.Get('NeutralPionTree')
    if inputTree is None: raise Usage('could not connect tree in file %s' % (paths['bemcPions'],))
    
    oldTracks = ROOT.TClonesArray('TPi0Candidate',100)
    inputTree.SetBranchAddress('CandidateBranch',oldTracks)
    inputTree.BuildIndex('runNumber','eventNumber')
    
    for i in range(spinTree.GetEntries()):
        spinTree.GetEntry(i)
        nBytesRead = inputTree.GetEntryWithIndex(event.runId(),event.eventId())
        if nBytesRead > 0:
            oldList = list(oldTracks)
            goodList = [track for track in oldList if acceptBemcPion(track)]
            nParticles[0] = len(goodList)
            for j in range(len(goodList)): goodTracks[j] = goodList[j]
        friendTree.Fill()
        goodTracks.Clear()
        nParticles[0] = -1
    
    inputFile.Close()
    
    saveSpinTree(outFile,spinTree,friendTree,[spinFile.Get('chargedPions'),
                                              spinFile.Get('bemcElectrons'),
                                              spinFile.Get('ConeJets'),
                                              spinFile.Get('ConeJetsEMC'),
                                              spinFile.Get('gammaCandidates')])
    
    spinFile.Close()
    
    print 'updated the bemcPions tree at %s' % (outPath,)


def updateBemcElectrons(inputPath):
    runnumber = getRunNumber(inputPath)
    
    spinFile = ROOT.TFile.Open(inputPath,'read')
    if spinFile is None: raise Usage('could not open spin file at %s' % (inputPath,))
    
    spinTree = spinFile.Get('spinTree')
    if spinTree is None: raise Usage('could not get spinTree in file at %s' % (inputPath,))
    
    event = ROOT.StJetSkimEvent()
    spinTree.SetBranchAddress('skimEventBranch',event)
    
    #create the output file now so friendTree is disk-resident
    try: os.mkdir('new')
    except OSError: pass
    outPath = 'new/spinAnalyses_%s.tree.root' % (runnumber,)
    outFile = ROOT.TFile(outPath,'recreate')
    copySpinTree(outFile,spinTree)
    
    friendTree = ROOT.TTree('bemcElectrons','updated %s' % (time.ctime(),))
    goodPrimaries   = ROOT.TClonesArray('StPrimaryElectron',100)
    goodGlobals     = ROOT.TClonesArray('StGlobalElectron',500)
    nParticles = array('b',[-1])
    friendTree.Branch('PrimaryElectrons',goodPrimaries)
    friendTree.Branch('GlobalElectrons',goodGlobals)
    friendTree.Branch('nBemcElectrons',nParticles,'nBemcElectrons/B')
    
    isRealFile = True
    
    inputName = '%s%s.e.pMuDst.root' % (paths['bemcElectrons'], runnumber)
    inputFile = ROOT.TFile.Open(inputName,'read')
    if inputFile is None or inputFile.IsOpen() == False:
        isRealFile = False
    
    if isRealFile:
        inputTree = inputFile.Get('electronTree')
        if inputTree is None: raise Usage('could not connect tree in file %s' % (inputName,))
        
        oldPrimaries   = ROOT.TClonesArray('StPrimaryElectron',100)
        oldGlobals     = ROOT.TClonesArray('StGlobalElectron',500)
        inputTree.SetBranchAddress('PrimaryElectrons',oldPrimaries)
        inputTree.SetBranchAddress('GlobalElectrons',oldGlobals)
        inputTree.BuildIndex('Run','Event')
    
    for i in range(spinTree.GetEntries()):
        spinTree.GetEntry(i)
        
        trigs = [ event.trigger(117001), event.trigger(137213), event.trigger(137221), event.trigger(5),
                  event.trigger(137222), event.trigger(137585), event.trigger(137611), event.trigger(137622) ]
        
        if isRealFile:
            nBytesRead = inputTree.GetEntryWithIndex(event.runId(),event.eventId())
            if nBytesRead > 0:
                nParticles[0] = 0 # found matching entry
                
                accept = False
                for trig in trigs:
                    if trig is not None:
                        if trig.didFire() or trig.shouldFire() > 0: accept = True
                if accept == False: 
                    friendTree.Fill()
                    continue
                
                liPrimaries = list(oldPrimaries)
                liGlobals   = list(oldGlobals)
                
                nParticles[0] = len(liPrimaries)
                
                copy1 = [track for track in liPrimaries]
                copy2 = [track for track in liGlobals]
                
                for j in range(len(copy1)): goodPrimaries[j]    = copy1[j]
                for j in range(len(copy2)): goodGlobals[j]      = copy2[j]
            
        friendTree.Fill()
        goodPrimaries.Clear()
        goodGlobals.Clear()
        nParticles[0] = -1
        
    if isRealFile: inputFile.Close()
    
    saveSpinTree(outFile,spinTree,friendTree, [ spinFile.Get('chargedPions'),
                                                spinFile.Get('bemcPions'),
                                                spinFile.Get('ConeJets'),
                                                spinFile.Get('ConeJetsEMC'),
                                                spinFile.Get('gammaCandidates') ] )
    
    spinFile.Close()
    
    print 'updated the bemcElectrons tree at %s' % (outPath,)


def setTriggerBits(jet, triggers, particles):
    if year == '2005': 
        patchPhi = [90., 30., -30., -90., -150., 150.]
        minPhi = 40.
        maxPhi = 320.
    if year == '2006':
        patchPhi = [150., 90., 30., -30., -90., -150., 150., 90., 30., -30., -90., -150.]
        minPhi = 36.
        maxPhi = 324.
    
    for trig in list(triggers):
        if trig.shouldFire()>0:
            #HT
            if trig.trigId() in (96201,96211):
                towers = trig.towersAboveThreshold(0)
                ROOT.SetOwnership(towers,True)
                try:
                    myiter = iter(towers)
                    while True:
                        mypair = myiter.next()
                        towerId = mypair.first
                        for t2j in particles:
                            if t2j.detectorId() == 9 and t2j.trackIndex() == towerId:
                                jet.addGeomTrigger(trig.trigId())
                except StopIteration: pass
                                
            #HTTP
            if trig.trigId() in (137822,):
                towers = trig.towersAboveThreshold(0)
                trigPatches = trig.triggerPatchesAboveThreshold(0)
                ROOT.SetOwnership(towers,True)
                ROOT.SetOwnership(trigPatches,True)
                try:
                    titer = iter(towers)
                    while True:
                        tpair = titer.next()
                        towerId = tpair.first
                        for t2j in particles:
                            if t2j.detectorId() == 9 and t2j.trackIndex() == towerId:
                                ROOT.gROOT.ProcessLine("decoder.GetTriggerPatchFromTowerId(%d,patchId);" % (towerId,))
                                try:
                                    piter = iter(trigPatches)
                                    while True:
                                        ppair = piter.next()
                                        if ROOT.patchId == ppair.first:
                                            jet.addGeomTrigger(trig.trigId())
                                except StopIteration: pass
                except StopIteration: pass
            
            #JP
            if trig.trigId() in (96221,96233,137221,137222):
                jp = trig.jetPatchesAboveThreshold(0)
                ROOT.SetOwnership(jp,True)
                try:
                    myiter = iter(jp)
                    while True:
                        patchId = myiter.next().first
                        dPhi = math.degrees(jet.Phi()) - patchPhi[patchId]
                        if math.fabs(dPhi) < minPhi or math.fabs(dPhi) > maxPhi:
                            if year == '2006' and jet.zVertex > -900:
                                if patchId <= 5 and jet.detEta() > -0.1 and jet.detEta() < 1.1:
                                    jet.addGeomTrigger(trig.trigId())
                                if patchId >= 6 and jet.detEta() > -1.1 and jet.detEta() < 0.1:
                                    jet.addGeomTrigger(trig.trigId())
                            if year == '2005':
                                jet.addGeomTrigger(trig.trigId())
                except StopIteration: pass


def updateJets(inputPath,jetBranchName):
    runnumber = getRunNumber(inputPath)
        
    spinFile = ROOT.TFile.Open(inputPath,'read')
    if spinFile is None: raise Usage('could not open spin file at %s' % (inputPath,))
    spinTree = spinFile.Get('spinTree')
    if spinTree is None: raise Usage('could not get spinTree in file at %s' % (inputPath,))
    
    event = ROOT.StJetSkimEvent()
    spinTree.SetBranchAddress('skimEventBranch',event)
    
    #create the output file now so friendTree is disk-resident
    try: os.mkdir('new')
    except OSError: pass
    outPath = 'new/spinAnalyses_%s.tree.root' % (runnumber,)
    outFile = ROOT.TFile(outPath,'recreate')
    copySpinTree(outFile,spinTree)
    
    newBranchName = jetBranchName
    if newBranchName == 'ConeJets12': newBranchName = 'ConeJets'
    jetFriend = ROOT.TTree(newBranchName,'updated %s' % (time.ctime(),))
    goodJets = ROOT.TClonesArray('StJet',100)
    nParticles = array('b',[-1])
    jetFriend.Branch(newBranchName,goodJets)
    jetFriend.Branch('n%s' % (newBranchName,),nParticles,'n%s/B' % (newBranchName,))
    
    jetName = '%s%s.tree.root' % (paths['jets'],runnumber)
    jetFile = ROOT.TFile.Open(jetName,'read')
    if jetFile is None:
        raise Usage('could not open file at %s' % (jetName,))
    jetTree = jetFile.Get('jet')
    if jetTree is None:
        raise Usage('could not connect tree in file %s' % (jetName,))
    
    bigJets = ROOT.StJets()
    jetTree.SetBranchAddress(jetBranchName,bigJets)
    jetTree.BuildIndex('mRunNumber','mEventNumber')
    
    for i in xrange(spinTree.GetEntries()):
        spinTree.GetEntry(i)
        nBytesRead = jetTree.GetEntryWithIndex(event.runId(),event.eventId())
        if nBytesRead > 0:
            theJets = bigJets.jets()
            
            #this messiness is temporary until the next pass over jet trees
            copyList = [ROOT.StJet(jet) for jet in list(theJets)]
            counter = 0
            for jet in copyList:
                if event.bestVert() is not None: jet.zVertex = event.bestVert().position()[2]
                else: jet.zVertex = -999
                setTriggerBits(jet, event.triggers(), bigJets.particles(counter))
                counter += 1
            goodList = [jet for jet in copyList if acceptJet(jet,event.bestVert())]
            
            nParticles[0] = len(goodList)
            
            for j in range(len(goodList)): goodJets[j] = goodList[j]
        jetFriend.Fill()
        goodJets.Clear()
        nParticles[0] = -1
    
    jetFile.Close()
    
    if newBranchName == 'ConeJets': otherJet = spinFile.Get('ConeJetsEMC')
    else: otherJet = spinFile.Get('ConeJets')
    saveSpinTree(outFile,spinTree,jetFriend,[ spinFile.Get('chargedPions'),
                                              spinFile.Get('bemcPions'),
                                              otherJet,
                                              spinFile.Get('bemcElectrons'),
                                              spinFile.Get('gammaCandidates')])
    
    spinFile.Close()
    
    print 'updated the %s branch at %s' % (newBranchName,outPath)


def updateGammaCandidates(inputPath):
    runnumber = getRunNumber(inputPath)
        
    spinFile = ROOT.TFile.Open(inputPath,'read')
    if spinFile is None: raise Usage('could not open spin file at %s' % (inputPath,))
    spinTree = spinFile.Get('spinTree')
    if spinTree is None: raise Usage('could not get spinTree in file at %s' % (inputPath,))
    
    event = ROOT.StJetSkimEvent()
    spinTree.SetBranchAddress('skimEventBranch',event)
    
    #create the output file now so friendTree is disk-resident
    try: os.mkdir('new')
    except OSError: pass
    outPath = 'new/spinAnalyses_%s.tree.root' % (runnumber,)
    outFile = ROOT.TFile(outPath,'recreate')
    copySpinTree(outFile,spinTree)
    
    friendTree = ROOT.TTree('gammas','updated %s' % (time.ctime(),))
    goodTracks = ROOT.TClonesArray('StGammaCandidate',100)
    nParticles = array('b',[-1])
    friendTree.Branch('gammas',goodTracks)
    friendTree.Branch('nGammas',nParticles,'nGammas/B')
    
    isRealFile = True
    
    #glob for matching input files then chain together
    filelist = glob('%s%s.root' % (paths['gammaCandidates'],runnumber))
    if len(filelist) == 0: 
        isRealFile = False
    else:
        inputFile = ROOT.TFile.Open(filelist[0],'read')
        inputTree = inputFile.Get('gammas')
    #inputTree = ROOT.TChain('gammas')
    #[inputTree.AddFile(f) for f in filelist]
    #inputTree.AddFile(filelist[0])
    
    if isRealFile:
        gammaEvent = ROOT.StGammaEvent()
        inputTree.SetBranchAddress('GammaEvent',gammaEvent)
        inputTree.BuildIndex('mRunNumber','mEventNumber')
    
    for i in range(spinTree.GetEntries()):
        spinTree.GetEntry(i)
        if isRealFile:
            nBytesRead = inputTree.GetEntryWithIndex(event.runId(),event.eventId())
            if nBytesRead > 0:
                #oldList = list(gammaEvent.candidates())
                #goodList = [track for track in oldList if acceptGammaCandidate(track)]
                goodList = []
                for i in range(gammaEvent.numberOfCandidates()):
                    g = gammaEvent.candidate(i)
                    if acceptGammaCandidate(g): goodList.append(g)
                nParticles[0] = len(goodList)
                for j in range(len(goodList)): goodTracks[j] = goodList[j]
        friendTree.Fill()
        goodTracks.Clear()
        nParticles[0] = -1
    
    #inputFile.Close()
        
    saveSpinTree(outFile,spinTree,friendTree,[spinFile.Get('bemcPions'),
                                              spinFile.Get('ConeJets'),
                                              spinFile.Get('ConeJetsEMC'),
                                              spinFile.Get('bemcElectrons'),
                                              spinFile.Get('chargedPions')])
    
    spinFile.Close()
    
    print 'updated the gammaCandidates tree at %s' % (outPath,)
    


def copySpinTree(f,spinTree):
    spinTree2 = spinTree.CloneTree()
    spinTree2.SetBranchStatus('*')
    ROOT.SetOwnership(spinTree2,True)
    
    #slightly ugly hack to keep TRefs working for trigHeaders
    headers = spinTree2.GetUserInfo()
    for i in xrange(headers.GetEntries()):
        headers.At(i).SetUniqueID(spinTree.GetUserInfo().At(i).GetUniqueID())
    
    f.cd()
    spinTree2.Write()
    del spinTree


def saveSpinTree(f,spinTree,friendTree,otherTrees):
    f.cd()
    del spinTree
    
    if friendTree is not None: friendTree.Write()
    
    #copy all particle trees into new file (saves space vs. overwriting in old file)
    for t in otherTrees:
        if t is not None:
            t2 = t.CloneTree()
            t2.SetBranchStatus('*')
            ROOT.SetOwnership(t2,True)
            t2.Write()
    f.Close()

def slimGammaEvent(oldEvent):
    ev = ROOT.StGammaEvent()
    ev.SetRunNumber(oldEvent.runNumber())
    ev.SetEventNumber(oldEvent.eventNumber())
    ev.SetMudstFileName(oldEvent.muDstFileName())
    ev.SetVertex(oldEvent.vertex())
    ev.SetMagneticField(oldEvent.magneticField())
    ev.SetPythia(oldEvent.pythia())
    
    for i in range(oldEvent.numberOfCandidates()):
        old = oldEvent.candidate(i)
        c = ev.newCandidate()
        c.SetMomentum(old.momentum())
        c.SetPosition(old.position())
        c.SetEnergy(old.energy())
        c.SetSeedEnergy(old.seedEnergy())
        c.SetPre1Energy(old.pre1Energy())
        c.SetPre2Energy(old.pre2Energy())
        c.SetPostEnergy(old.postEnergy())
        c.SetSmduEnergy(old.smduEnergy())
        c.SetSmdvEnergy(old.smdvEnergy())
        
        for j in range(old.numberOfMyTracks()):
            ev.newTrack(old.mytrack(j))
            c.addMyTrack(ev.track(ev.numberOfTracks()-1))
        
        for j in range(old.numberOfMyTowers()):
            ev.newTower(old.mytower(j))
            c.addMyTower(ev.tower(ev.numberOfTowers()-1))
        
        for j in range(old.numberOfMyPreshower1()):
            ev.newPreshower1(old.mypreshower1(j))
            c.addMyPreshower1(ev.preshower1(ev.numberOfPreshower1()-1))
        
        for j in range(old.numberOfMyPreshower2()):
            ev.newPreshower2(old.mypreshower2(j))
            c.addMyPreshower2(ev.preshower2(ev.numberOfPreshower2()-1))
            
        for j in range(old.numberOfMyPostshower()):
            ev.newPostshower(old.mypostshower(j))
            c.addMyPostshower(ev.postshower(ev.numberOfPostshower()-1))

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hcu:", ['help',"create", "update=", 'site='])
        except getopt.error, msg: raise Usage(msg)
        
        global year
        global site
        global paths
        
        #magic!
        throwitaway = ROOT.std.map(long,float)()
        
        #default site = MIT
        libs = loadLibrariesMIT
        
        if len(args) == 0: raise Usage('no files specified')
        
        if args[0].endswith('.list'):
            filelist = file(args[0])
            inputFiles = filelist.readlines()
            inputFiles = [f.strip() for f in inputFiles] #get rid of weird line endings
        else:
            inputFiles = args
        
        print inputFiles
        
        #paths
        firstRun = getRunNumber(inputFiles[0])
        if firstRun.startswith('6'): year = '2005'
        if firstRun.startswith('7'): year = '2006'
        paths = allPaths[year][site]
        
        # option processing
        for option, value in opts:
            if option in ('-h','--help'): raise Usage(help_message)
            
            if option == '--site' and value in ('RCF','PDSF'): 
                libs = loadLibraries
                site = value
                paths = allPaths[year][site]
            
            if option == '--site' and value == 'Xgrid': 
                libs = loadLibrariesXgrid
                site = value
                paths = allPaths[year][site]
            
            if option in ('-c','-u','--create','--update'): libs()
            
            if option in ("-c", "--create"): [makeSpinTree(elem,'./') for elem in inputFiles]
            
            if option in ("-u", "--update"):
                if value in ('ConeJets','ConeJets12','ConeJetsEMC'):
                    ROOT.gROOT.ProcessLine("int patchId; StEmcDecoder decoder;")
                    [updateJets(elem,value) for elem in inputFiles]
                if value == 'chargedPions':
                    [updateChargedPions(elem) for elem in inputFiles]
                if value == 'bemcPions':
                    [updateBemcPions(elem) for elem in inputFiles]
                if value == 'bemcElectrons':
                    [updateBemcElectrons(elem) for elem in inputFiles]
                if value == 'gammas':
                    [updateGammaCandidates(elem) for elem in inputFiles]
                if value == 'skim':
                    [fixSkimTree(elem,'./') for elem in inputFiles]
    
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
