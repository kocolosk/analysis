import ROOT, minimc, math, os, uuid
from datetime import datetime

spin4 = {'uu':5, 'du':6, 'ud':9, 'dd':10}

def loadSpinTree():
    ROOT.gSystem.Load('libTable')
    ROOT.gSystem.Load('StarRoot')
    ROOT.gSystem.Load('StarClassLibrary')
    ROOT.gROOT.Macro('loadMuDst.C')
    ROOT.gSystem.Load('StDaqLib')
    ROOT.gSystem.Load('StDetectorDbMaker')
    ROOT.gSystem.Load('StEmcTriggerMaker')
    ROOT.gSystem.Load('StMCAsymMaker')
    ROOT.gSystem.Load('StSpinDbMaker')
    ROOT.gSystem.Load('StJetFinder')
    ROOT.gSystem.Load('StJetMaker')
    ROOT.gSystem.Load('StChargedPionAnalysisMaker')
    ROOT.gSystem.Load('StSpinTree')

class EventInfo:
    """simple class to store event-level stuff we care about"""
    def __init__(self,):
        self.clear()
    
    
    def clear(self):
        """docstring for clear"""
        self.run        = -1
        self.event      = -1
        self.bbcTimeBin = -1
        self.spinBit    = -1
        self.nVertices  = -1
        self.vertices   = []
        self.trigs      = []
        self.simuTrigs  = []
    
    
    def fill(self,skimEvent):
        self.clear()
        
        self.run        = skimEvent.runId()
        self.event      = skimEvent.eventId()
        self.bbcTimeBin = skimEvent.bbcTimeBin() / 32
        self.nVertices  = skimEvent.vertices().GetEntries()
        
        if skimEvent.isValid() and skimEvent.isPolLong() and not skimEvent.isPolTrans() and \
            not skimEvent.isMaskedUsingBx48() and skimEvent.offsetBx48minusBX7() == 0:
            self.spinBit = skimEvent.spin4usingBx48()
        
        for v in skimEvent.vertices():
            #self.vertices.append( ROOT.TVector3(v.position()[0], v.position()[1], v.position()[2]) )
            self.vertices.append( v.position() )
        
        for t in skimEvent.triggers():
            if t.didFire() > 0:     
                self.trigs.append(t.trigId())
                
            if t.shouldFire() > 0:
                self.simuTrigs.append(t.trigId())
                
            if t.trigId() == 96011 and skimEvent.eBbc() > 0 and skimEvent.wBbc() > 0:
                self.simuTrigs.append(t.trigId())
    
    def fill2(self,myEvent):
        self.clear()
        
        self.run        = myEvent.runId()
        self.event      = myEvent.eventId()
        self.bbcTimeBin = myEvent.bbcTimeBin() / 32
        self.nVertices  = myEvent.nVertices()
        self.spinBit    = myEvent.isSpinValid() and myEvent.spinBit() or 0
        
        self.vertices   = myEvent.vertices()
        
        for trigId in (96011, 96201, 96211, 96221, 96233):
            if myEvent.isTrigger(trigId): 
                self.trigs.append(trigId)
            if myEvent.isSimuTrigger(trigId):
                self.simuTrigs.append(trigId)
        
    def __str__(self):
        return 'R%d E%d BBC:%2d Spin4:%2d nVerts:%2d vz:' % \
            (self.run, self.event, self.bbcTimeBin, self.spinBit, self.nVertices) + \
            str([v.z() for v in self.vertices]) + ' trigs:' + str(self.trigs) + ' simuTrigs:' + str(self.simuTrigs)
    

class EventCuts:
    def __init__(self, triggers):
        self.spin    = False
        self.vertex  = False
        self.bbc     = False
        self.trigger = False
        self.all     = False
        
        self.__triggers = triggers
        
    def set(self, event):
        self.spin   = event.isSpinValid()
        self.vertex = event.nVertices() > 0
        
        timeDiff = event.bbcTimeBin()
        if timeDiff % 32 != 0:
            self.bbc = False
        else:
            bin = timeDiff / 32
            self.bbc = bin in (7,8,9) or (event.runId() > 7000000 and bin == 6)
            
        self.trig = False
        for trig in self.__triggers:
            if event.isTrigger(trig) and event.isSimuTrigger(trig):
                self.trigger = True
                
        self.all = self.spin and self.vertex and self.bbc and self.trigger
        
    


class TrackCuts:
    def __init__(self):
        self.eta = False
        self.dca = False
        self.fit = False
        self.pid = False
        self.all = False
    
    def set(self, track):
        self.eta = math.fabs( track.eta() ) < 1.0
        self.dca = math.fabs( track.globalDca().mag() ) < 1.0
        self.pid = track.nSigmaPion() > -1.0 and track.nSigmaPion() < 2.0
        self.fit = track.nHitsFit() > 25
        self.all = self.eta and self.dca and self.pid and self.fit

    
class DataHistos2(dict):
    """generates and manages raw histograms from data.  All cuts are contained here.
    This version uses 2-D histograms to save spin-sorted histograms.  A 4bit-2bit converter
    is given as a dict class member."""
    spin2 = {0:0, 5:1, 6:2, 9:3, 10:4, 'all':0, 'any':0, 'uu':1, 'du':2, 'ud':3, 'dd':4}
    
    mcPtBins        = minimc.MiniMcHistos.ptBins
    mcVzBins        = minimc.MiniMcHistos.vzBins
    
    def __init__(self, name, eventCuts, tfile=None):
        super(DataHistos2, self).__init__()
        
        self.name       = name
        self.spinIndex  = 0
        self.ecuts      = eventCuts
        
        self['pt'] = ROOT.TH2D('_%s_pt' % (name,), '', 4, 0.5, 4.5, self.mcPtBins[0], self.mcPtBins[1], self.mcPtBins[2])
        self['vz'] = ROOT.TH2D('_%s_vz' % (name,), '', 4, 0.5, 4.5, self.mcVzBins[0], self.mcVzBins[1], self.mcVzBins[2])
    
    def setSpin(self, spin4):
        """return spin-summed histos by default.  Alternatively, one can select a spinbit here"""
        try:
            self.spinbit = self.spin2[spin4]
        except KeyError:
            print spin4, 'is not a valid 4-bit spin'
    
    def fillEvent(self, ev):
        """fill all appropriate histos from an StChargedPionEvent"""
        self.spinIndex = self.ecuts.spin and DataHistos2.spin2[ev.spinBit()] or 0
        if self.ecuts.vertex:
            self['vz'].Fill(self.spinIndex, ev.vertex(0).z())
            
    def fillTrack(self, track, tcuts):
        if self.ecuts.all:
            if tcuts.all:
                self['pt'].Fill(self.spinIndex, track.pt())
    
    def get(self, key, spin4=None):
        """return a given histogram.  Prefer this to dict lookups b/c of spin-sorting"""
        localspin = spin4 and self.spin2[spin4] or self.spinbit
        newname = self[key].GetName() + str(uuid.uuid1())
        if localspin > 0:
            return self[key].ProjectionY((newname), localspin, localspin)
        else:
            return self[key].ProjectionY((newname), 1, 4)
        
    


def testDataHistos2():
    chain = ROOT.TChain('tree')
    chain.Add('~/work/charged-pion-event-2/chargedPions_6119032.tree.root')
    chain.GetEntry(0)
    event = chain.event
    
    elistFile = ROOT.TFile('~/work/charged-pion-event-2/eventLists.root')
    elist = elistFile.Get('has_trigger')
    #elist.Intersect( elistFile.Get('has_track') )
    
    ecuts = (
        EventCuts( triggers=(96011,)                             ),
        EventCuts( triggers=(96201,)                             ),
        EventCuts( triggers=(96211,)                             ),
        EventCuts( triggers=(96221,)                             ),
        EventCuts( triggers=(96233,)                             ),
        EventCuts( triggers=(96011, 96201, 96211)               ),
        EventCuts( triggers=(96011, 96221, 96233)               ),
        EventCuts( triggers=(96011, 96201, 96211, 96221, 96233) ),
        )
    
    plusFile = ROOT.TFile('test.plus.hist.root', 'recreate')
    hPlus = (
        DataHistos2('96011',        ecuts[0]),
        DataHistos2('96201',        ecuts[1]),
        DataHistos2('96211',        ecuts[2]),
        DataHistos2('96221',        ecuts[3]),
        DataHistos2('96233',        ecuts[4]),
        DataHistos2('hightower',    ecuts[5]),
        DataHistos2('jetpatch',     ecuts[6]),
        DataHistos2('alltrigs',     ecuts[7]),
        )
        
    minusFile = ROOT.TFile('test.minus.hist.root', 'recreate')
    hMinus = (
        DataHistos2('96011',        ecuts[0]),
        DataHistos2('96201',        ecuts[1]),
        DataHistos2('96211',        ecuts[2]),
        DataHistos2('96221',        ecuts[3]),
        DataHistos2('96233',        ecuts[4]),
        DataHistos2('hightower',    ecuts[5]),
        DataHistos2('jetpatch',     ecuts[6]),
        DataHistos2('alltrigs',     ecuts[7]),
        )
            
    allHists = hPlus + hMinus
    
    tcuts = TrackCuts()
    
    fname = ''
    for i in range(elist.GetN()):
        if chain.GetEntry(elist.GetEntry(i)) == 0: break
        if fname != chain.GetCurrentFile().GetName():
           fname = chain.GetCurrentFile().GetName()
           print fname
        
        [ e.set(event) for e in ecuts ]
        
        [ h.fillEvent(event) for h in allHists ]
        
        #hists.fillEvent(event, ecuts)        
        for track in event.tracks():
            tcuts.set(track)
            #hists.fillTrack(track, tcuts)
            if track.charge() == 1:
                [ h.fillTrack(track, tcuts) for h in hPlus ]
            else:
                [ h.fillTrack(track, tcuts) for h in hMinus ]
    
    print 'Goodbye'
    elistFile.Close()
    
    plusFile.cd()
    for dh in hPlus:
        [ h.Write() for h in dh.values() ]
    plusFile.Close()
    
    minusFile.cd()
    for dh in hMinus:
        [ h.Write() for h in dh.values() ]
    minusFile.Close()
    
    
    
class DataHistos(dict):
    """generates and manages raw histograms from data.  All cuts are contained here"""
    mcPtBins        = minimc.MiniMcHistos.ptBins
    mcEtaBins       = minimc.MiniMcHistos.etaBins
    mcPhiBins       = minimc.MiniMcHistos.phiBins
    mcVzBins        = minimc.MiniMcHistos.vzBins
    mcNFitBins      = minimc.MiniMcHistos.nFitBins
    mcPBins         = minimc.MiniMcHistos.pBins
    mcDEdxBins      = minimc.MiniMcHistos.dEdxBins
    mcDcaGBins      = minimc.MiniMcHistos.dcaGBins
    
    def __init__(self, trigIdString, spinString='anyspin', tfile=None):
        super(DataHistos, self).__init__()
        
        self.clear()
        
        self.trigIdString = trigIdString        
        if 'all' in self.trigIdString:
            self.triggers = [96011, 96201, 96211, 96221, 96233]
        elif 'mb' or 'minbias' in self.trigIdString:
            self.triggers = [96011]
        elif 'ht' or 'hightower' in self.trigIdString:
            self.triggers = [96011, 96201, 96211]
        elif 'jp' or 'jetpatch' in self.trigIdString:
            self.triggers = [96011, 96221, 96233]
        else:
            self.triggers = [int(self.trigIdString)]
            
        self.spinString = spinString
        if 'any' in self.spinString:
            self.spinBits = [5,6,9,10]
        else:
            self.spinBits = [spin4[spinString]]
        
        #for trigIdString in ('96011', '96201', '96211', '96221', '96233', 'all', 'ht', 'jp'):
        #    for spinString in ('anyspin', 'uu', 'ud', 'du', 'dd'):
        self['pt'] = ROOT.TH1D('_%s_pt_%s' % (trigIdString,spinString), '', self.mcPtBins[0], self.mcPtBins[1], self.mcPtBins[2])
        self['eta'] = ROOT.TH1D('_%s_eta_%s' % (trigIdString,spinString), '', self.mcEtaBins[0], self.mcEtaBins[1], self.mcEtaBins[2])
        self['phi'] = ROOT.TH1D('_%s_phi_%s' % (trigIdString,spinString), '', self.mcPhiBins[0], self.mcPhiBins[1], self.mcPhiBins[2])
        self['nHitsFit'] = ROOT.TH1D('_%s_nHitsFit_%s' % (trigIdString,spinString), '', self.mcNFitBins[0], self.mcNFitBins[1], self.mcNFitBins[2])
        self['dEdx'] = ROOT.TH1D('_%s_dEdx_%s' % (trigIdString,spinString), '', self.mcDEdxBins[0], self.mcDEdxBins[1], self.mcDEdxBins[2])
        self['dcaG'] = ROOT.TH1D('_%s_dcaG_%s' % (trigIdString,spinString), '', self.mcDcaGBins[0], self.mcDcaGBins[1], self.mcDcaGBins[2])
        
        self['nVertices'] = ROOT.TH1D('_%s_nVertices_%s' % (trigIdString,spinString), '', 15, -0.5, 14.5)
        self['vz'] = ROOT.TH1D('_%s_vz_%s' % (trigIdString,spinString), '', self.mcVzBins[0], self.mcVzBins[1], self.mcVzBins[2])        
        
        ## vz after requiring nVertices > 0 and BBC in [(6),7,8,9]
        self['nVerticesBBC'] = ROOT.TH1D('_%s_nVerticesBBC_%s' % (trigIdString,spinString), '', 15, -0.5, 14.5)
        self['vzBBC'] = ROOT.TH1D('_%s_vzBBC_%s' % (trigIdString,spinString), '', self.mcVzBins[0], self.mcVzBins[1], self.mcVzBins[2])
        
        if tfile is not None:
            for key in self.keys():
                self[key].Delete()
                self[key] = tfile.Get('_%s_%s_%s' % (trigIdString,key,spinString))
        
    
    
    def clear(self):
        """reset any stored data from the previous event"""        
        ## general event cuts
        self.vertCut = False
        self.trigCut = False
        self.spinCut = False
    
    
    def checkCuts(self, ev):
        self.clear()
        
        ## vertex cut (TPC vertex exists and BBC timebin in specified range)
        if ev.run > 7000000:
            self.vertCut = ev.nVertices > 0 and ev.bbcTimeBin in [6,7,8,9]
        else:
            self.vertCut = ev.nVertices > 0 and ev.bbcTimeBin in [7,8,9]
        
        ## trigger cut (both hardware and software required)
        for trigId in self.triggers:
            if trigId in ev.trigs and trigId in ev.simuTrigs:
                self.trigCut = True
                break
        
        ## spinCut (spinDb good and matches spinbit we care about)
        self.spinCut = ev.spinBit in self.spinBits
    
    def checkCuts2(self, ev):
        self.clear()
        
        ## vertex cut (TPC vertex exists and BBC timebin in specified range)
        if ev.runId() > 7000000:
            self.vertCut = ev.nVertices() > 0 and ev.bbcTimeBin() % 32 == 0 and ev.bbcTimeBin() / 32 in (6,7,8,9)
        else:
            self.vertCut = ev.nVertices() > 0 and ev.bbcTimeBin() % 32 == 0 and ev.bbcTimeBin() / 32 in (7,8,9)
        
        ## trigger cut (both hardware and software required)
        for trigId in self.triggers:
            if ev.isTrigger(trigId) and ev.isSimuTrigger(trigId):
                self.trigCut = True
                break
            
        
    
    def fillEvent(self, ev):
        """"check cuts and fill event-level histograms"""
        self.checkCuts2(ev)
        
        ## event-wise histograms
        if self.trigCut and self.spinCut:
            self['nVertices'].Fill(ev.nVertices())
            for v in ev.vertices(): self['vz'].Fill(v.z())
            
        if self.trigCut and self.spinCut and self.vertCut:
            self['vzBBC'].Fill(ev.vertex(0).z())
            
        
    
    
    def fillTrack(self, track, cuts):
        if self.trigCut and self.spinCut and self.vertCut:
            if cuts.eta and cuts.dca and cuts.fit and cuts.pid:
                self['pt'].Fill(track.pt())
                self['phi'].Fill(track.phi())
                self['dEdx'].Fill(track.dEdx() * 1e7)
            
            if cuts.dca and cuts.fit and cuts.pid:            
                self['eta'].Fill(track.eta())
            
            if cuts.eta and cuts.fit and cuts.pid:
                self['dcaG'].Fill(track.globalDca().mag())
            
            if cuts.eta and cuts.dca and cuts.pid:
                self['nHitsFit'].Fill(track.nHitsFit())
            
    


def makeDataHistograms(charges=[-1,1],runnumber=None):
    #reader = ROOT.StSpinTreeReader()
    #
    #reader.selectDataset('/Users/kocolosk/data/run5/dataset.txt')
    #
    #reader.selectTrigger(96011)
    #reader.selectTrigger(96201)
    #reader.selectTrigger(96211)
    #reader.selectTrigger(96221)
    #reader.selectTrigger(96233)
    #
    #reader.requireDidFire   = True
    #reader.requireShouldFire = True
    #
    #reader.connectJets          = False
    #reader.connectNeutralJets   = False
    #reader.connectBemcPions     = False
    #reader.connectBemcElectrons = False
    #reader.connectEemcPions     = False
    #
    event = EventInfo()
    
    reader = ROOT.TChain('tree')
    reader.Add('~/work/charged-pion-event-2/chargedPions_611903*')
    reader.SetBranchStatus('mJets*',0)
    
    elistFile = ROOT.TFile('~/work/charged-pion-event-2/eventLists.root')
    elist = elistFile.Get('has_trigger')
    
    cuts  = TrackCuts()
    
    if runnumber is None: 
        name = 'allruns'
        #reader.selectRunlist('/Library/STAR/10.5/star/DEV/StRoot/StSpinPool/StSpinTree/filters/run5_chargedPions.runlist')
    else: 
        name = str(runnumber)
        #reader.selectRun(runnumber)
    
    outPlus  = ROOT.TFile('%s.plus.hist.root' % (name,), 'recreate')
    dataPlus = []
    for trig in ('all', 'mb', 'ht', 'jp', '96201', '96211', '96221', '96233'):
        for spin in ('any', 'uu', 'ud', 'du', 'dd'):
            dataPlus.append( DataHistos(trig, spin) )
            
    outMinus  = ROOT.TFile('%s.minus.hist.root' % (name,), 'recreate')
    dataMinus = []
    for trig in ('all', 'mb', 'ht', 'jp', '96201', '96211', '96221', '96233'):
        for spin in ('any', 'uu', 'ud', 'du', 'dd'):
            dataMinus.append( DataHistos(trig, spin) )
    
    outSum  = ROOT.TFile('%s.sum.hist.root' % (name,), 'recreate')
    dataSum = []
    for trig in ('all', 'mb', 'ht', 'jp', '96201', '96211', '96221', '96233'):
        for spin in ('any', 'uu', 'ud', 'du', 'dd'):
            dataSum.append( DataHistos(trig, spin) )
    
    entries = reader.GetEntries()
    #entries = 100000
    print 'Histogramming %d events ...' % (entries,)
    for i in range(entries):
        if i % 5000 == 1:
            print 'processing event %-10d [ %5.2f%% ]' % (i, float(100*i)/entries)
        
        if reader.GetEntry(elist.GetEntry(i)) == 0: break
        event.fill2(reader.event)
        
        ## let's try checking the cuts here
        #ev = reader.event
        #hasVert = ev.nVertices() > 0
        #bbcOk = ev.bbcTimeBin() % 32 == 0 and ev.bbcTimeBin() / 32 in (7,8,9)
        
        
        [ d.fillEvent(reader.event) for d in dataSum  ]
        [ d.fillEvent(reader.event) for d in dataPlus  ]
        [ d.fillEvent(reader.event) for d in dataMinus ]
        
        for track in reader.event.tracks():
            cuts.set(track)
            
            [ d.fillTrack(track, cuts) for d in dataSum ]
            
            if track.charge() == 1:
                [ d.fillTrack(track, cuts) for d in dataPlus ]
                
            if track.charge() == -1:
                [ d.fillTrack(track, cuts) for d in dataMinus ]
    
    outPlus.cd()
    for d in dataPlus: [ h.Write() for h in d.values() ]
    outPlus.Close()
    
    outMinus.cd()
    for d in dataMinus: [ h.Write() for h in d.values() ]
    outMinus.Close()
    
    outSum.cd()
    for d in dataSum: [ h.Write() for h in d.values() ]
    outSum.Close()


def awesomeTree(runnumber=6119032):
    """generator of charged pion event trees"""
    skimFile = ROOT.TFile('~/data/run5/jetSkim/jetSkim_%d.tree.root' % (runnumber,))
    
    trackFile = ROOT.TFile('~/data/run5/tracks/chargedPions_%d.tree.root' % (runnumber,))
    trackFile.chargedPionTree.BuildIndex('run','event')
    
    #jetFile = ROOT.TFile('~/data/run5/jets/jets_%d.tree.root' % (runnumber,))
    jetFile = ROOT.TFile('/Volumes/scratch/common/run5/jets/jets_%d.tree.root' % (runnumber,))
    jetFile.jet.BuildIndex('mRunId','mEventId')
    jetFile.jet.GetEntry(0)
    ConeJets = jetFile.jet.ConeJets
    
    ## keep this tree on a diet ... no need for extra 8 bytes / object
    ## this code should be in the library now
    #ROOT.StChargedPionEvent.Class().IgnoreTObjectStreamer()    
    #ROOT.StChargedPionVertex.Class().IgnoreTObjectStreamer()
    #ROOT.StChargedPionTrack.Class().IgnoreTObjectStreamer()
    #ROOT.StChargedPionJet.Class().IgnoreTObjectStreamer()
    #ROOT.StChargedPionJetParticle.Class().IgnoreTObjectStreamer()
    
    outFile = ROOT.TFile('chargedPions_%d.tree.root' % (runnumber,), 'recreate')
    outTree = ROOT.TTree('tree', 'created %s' % (str(datetime.now()),))
    ev = ROOT.StChargedPionEvent()
    outTree.Branch('event','StChargedPionEvent',ev)
    
    myJet = ROOT.StChargedPionJet()
    
    for entry in skimFile.jetSkimTree:
        # event info
        sk = skimFile.jetSkimTree.skimEventBranch
        ROOT.StChargedPionMaker.translateEvent(sk, ev)
        
        # charged pion info
        trackFile.chargedPionTree.GetEntryWithIndex(sk.runId(), sk.eventId())
        assert(trackFile.chargedPionTree.event == ev.eventId())
        for track in trackFile.chargedPionTree.primaries:
            ev.addTrack(track)
            
        # jet info using maker.translateJet
        jetFile.jet.GetEntryWithIndex(sk.runId(), sk.eventId())
        assert(ConeJets.eventId() == ev.eventId())
        ROOT.StChargedPionMaker.translateJets(ConeJets, ev)
        
        outTree.Fill()
        ev.Clear()
    
    outFile.cd()
    outTree.Write()
    outFile.Close()


def veryAwesome():
    """driver for awesome tree function"""
    dir = '/Users/kocolosk/data/run5/jetSkim/'
    files =  os.listdir(dir)
    for f in files[360:]:
       if f.endswith('.root'):
          print f
          awesomeTree(int(f[8:15]))


def measureReadSpeed():
    import time
    chain = ROOT.TChain('tree')
    
    ROOT.StChargedPionJetParticle.Class().IgnoreTObjectStreamer()
    chain.Add('~/work/charged-pion-event-2/chargedPions_61*')
    
    chain.SetBranchStatus('mJets*',0)
    #chain.SetBranchStatus('mJets.mParticles*',0)
    
    #elistFile = ROOT.TFile('~/work/charged-pion-event/eventLists.root')    
    elistFile = ROOT.TFile('~/work/charged-pion-event-2/eventLists.root')
    elist = elistFile.Get('all_cuts')
    
    counter = 0
    beginTime = time.clock()
    for i in range(elist.GetN()):
        if chain.GetEntry(elist.GetEntry(i)) == 0: break
        #for entry in chain:
        ev = chain.event
        if counter % 1000 == 1: print ev.runId(), ev.eventId()
        counter += 1
        mb = ev.isTrigger(96011) and ev.isSimuTrigger(96011)
        nVerts = ev.nVertices()
        for jet in ev.jets():
            pt = jet.Pt()
            for particle in jet.particles():
                det = particle.detectorId()
        for track in ev.tracks():
            eta = track.Eta()
    endTime = time.clock()
    print 'Processed', counter, ' events in ', endTime-beginTime, ' CPU seconds'
        
    
    print 'Goodbye'


def test2():
    f = ROOT.TFile('~/work/charged-pion-event/chargedPions_6119063.tree.root')
    f.tree.GetEntry(0)
    ev = f.tree.event
    
    for entry in f.tree:
        print ev.runId(), ev.eventId()
        for track in ev.tracks():
            print '...', track.pt(), track.eta(), track.phi()
    
    print f.tree.GetEntries()

loadSpinTree()
#makeDataHistograms([-1],6119063)
#ROOT.StChargedPionJetParticle.Class().IgnoreTObjectStreamer()
#measureReadSpeed()

def elist1():
    f = ROOT.TFile('~/work/charged-pion-event/chargedPions_6119063.tree.root')
    f.tree.Draw('>>elist1','mTracks.fE>0')
    elist = ROOT.gDirectory.Get('elist1')
    
    for i in range(elist.GetN()):
        f.tree.GetEntry(elist.GetEntry(i))
        print f.tree.event.eventId(), f.tree.event.nTracks()
    #f.tree.GetEn

def particleTuple(runnumber):
    """silly attempt at schema evolution"""
    f = ROOT.TFile('~/work/charged-pion-event/chargedPions_%d.tree.root' % (runnumber,))
    
    outFile = ROOT.TFile('particles_%d.nt.root' % (runnumber,),'recreate')
    nt = ROOT.TNtuple('nt','','event:pt:eta:phi:e:index:det:jetIndex')
    
    for entry in f.tree:
        ev = f.tree.event
        for row,jet in enumerate(ev.jets()):
            for p in jet.particles():
                nt.Fill(ev.eventId(), p.Pt(), p.Eta(), p.Phi(), p.E(), p.index(), p.detectorId(), row)
    
    nt.Write()
    outFile.Close()

def makeParticles():
    """driver for particle tuple function"""
    ROOT.StChargedPionJetParticle.Class().IgnoreTObjectStreamer()
    dir = '/Users/kocolosk/work/charged-pion-event/'
    files =  os.listdir(dir)
    for f in files[697:]:
        if 'chargedPions_' in f:
            print f
            particleTuple(int(f[13:20]))
    
def recombinate(runnumber):
    f1 = ROOT.TFile('~/work/charged-pion-event/chargedPions_%d.tree.root' % (runnumber,))
    f1.tree.SetBranchStatus('mJets.mParticles', 0)
    
    f2 = ROOT.TFile('~/work/charged-pion-event/particles_%d.nt.root' % (runnumber,))
    
    outFile = ROOT.TFile('chargedPions_%d.tree.root' % (runnumber,), 'recreate')
    outTree = ROOT.TTree('tree', 'created %s' % (str(datetime.now()),))
    ev = ROOT.StChargedPionEvent()
    outTree.Branch('event','StChargedPionEvent',ev)
    
    jet = ROOT.StChargedPionJet()
    particle = ROOT.StChargedPionJetParticle()
    
    ntentry = 0
    maxentries = f2.nt.GetEntries()
    
    for entry in f1.tree:
        old = f1.tree.event
        ev.copy(old)
        #print old.eventId(), ev.eventId()
        #ev.setRunId( old.runId() )
        #ev.setEventId( old.eventId() )
        #ev.setBx7( old.bx7() )
        #ev.setBbcTimeBin( old.bbcTimeBin() )
        #
        #ev.setSpinBit( old.spinBit() )
        #ev.setPolValid( old.isPolValid() )
        #ev.setPolLong( old.isPolLong() )
        #ev.setPolTrans( old.isPolTrans() )
        #ev.setBxingMasked( old.isBxingMasked() )
        #ev.setBxingOffset( not old.isNullOffset() )
        #
        #for trig in (96011,96201,96211,96221,96233):
        #    if old.isTrigger(trig): ev.addTrigger(trig)
        #    if old.isSimuTrigger(trig): ev.addSimuTrigger(trig)
        #
        #for v in old.vertices():
        #    ev.addVertex(v)
        #    
        #for t in old.tracks():
        #    ev.addTrack(t)
        
        for row,jet in enumerate(ev.jets()):
            #jet.SetPx( j.Px() )
            #jet.SetPy( j.Py() )
            #jet.SetPz( j.Pz() )
            #jet.SetE( j.E() )
            #jet.setCharge( j.charge() )
            #jet.setNTpcTracks( j.nTpcTracks() )
            #jet.setNBarrelTowers( j.nBarrelTowers() )
            #jet.setNEndcapTowers( j.nEndcapTowers() )
            #jet.setTpcEtSum( j.tpcEtSum() )
            #jet.setBarrelEtSum( j.barrelEtSum() )
            #jet.setEndcapEtSum( j.endcapEtSum() )
            #jet.setVertexZ( j.vertexZ() )
            
            if ntentry<maxentries:
                f2.nt.GetEntry(ntentry)
                while int(f2.nt.event) == ev.eventId() and int(f2.nt.jetIndex) == row:
                    particle.SetPt( f2.nt.pt )
                    particle.SetEta( f2.nt.eta )
                    particle.SetPhi( f2.nt.phi )
                    particle.SetE( f2.nt.e )
                    particle.setIndex( int(f2.nt.index) )
                    particle.setDetectorId( int(f2.nt.det) )
                    jet.addParticle( particle )
                    ntentry += 1
                    if ntentry == maxentries: break
                    f2.nt.GetEntry(ntentry)
                
            #ev.addJet( jet )
            #jet.Clear()
            
        outTree.Fill()
        ev.Clear()
        
    outFile.cd()
    outTree.Write()
    outFile.Close()
    


def recombineAll():
    """driver for particle tuple function"""
    dir = '/Users/kocolosk/work/charged-pion-event/'
    files =  os.listdir(dir)
    for f in files[42:]:
        if 'chargedPions_' in f:
            print f
            recombinate(int(f[13:20]))


def makeEventLists(fname='~/work/charged-pion-event-2/chargedPions_6119063.tree.root'):
    #ROOT.StChargedPionJetParticle.Class().IgnoreTObjectStreamer()
    #f = ROOT.TFile(fname, 'update')
    #ch = f.tree
    ch = ROOT.TChain('tree')
    ch.Add('~/data/run5/tree/chargedPions_*.tree.root')
    
    elists = []
    
    outFile = ROOT.TFile('eventLists.root','recreate')
    
    print 'generating has_trigger list'
    trigIds = [96011,96201,96211,96221,96233]
    trigString = ''
    for trig in trigIds[:-1]:
        trigString += 'event.isTrigger(%d) || ' % (trig,)
    trigString += 'event.isTrigger(%d)' % (trigIds[-1])
    ch.Draw('>>has_trigger',trigString)
    elists.append( ROOT.gDirectory.Get('has_trigger') )
    
    print 'generating has_vertex list'
    ch.Draw('>>has_vertex','mVertices.mRanking>0')
    elists.append( ROOT.gDirectory.Get('has_vertex') )
    
    print 'generating bbc_789 list'
    #ch.Draw('>>bbc_789','mBbcTimeBin%32==0 && mBbcTimeBin/32>6 && mBbcTimeBin/32<9')
    ## dropping the discrete timebin cut -- APK 2008-01-09
    ## also fix a BUG!!! 9 vs. 10
    ch.Draw('>>bbc_789','mBbcTimeBin%32==0 && mBbcTimeBin/32>6 && mBbcTimeBin/32<10')
    elists.append( ROOT.gDirectory.Get('bbc_789') )
    
    print 'generating spinDbOk list'
    ch.Draw('>>spinDbOk','event.isSpinValid()')
    elists.append( ROOT.gDirectory.Get('spinDbOk') )
    
    print 'generating has_track list'
    ch.Draw('>>has_track','mTracks.fE>0')
    elists.append( ROOT.gDirectory.Get('has_track') )
    
    print 'generating has_good_pion list'
    ch.Draw('>>has_good_pion','abs(mTracks.eta())<1.0 && mTracks.globalDca().mag()<1.0 && mTracks.nHitsFit()>25 && mTracks.nSigmaPion()>-1 && mTracks.nSigmaPion()<2')
    elists.append( ROOT.gDirectory.Get('has_good_pion') )
    
    print 'now get intersection'
    all_cuts = ROOT.TEventList('all_cuts','intersection of above lists')
    all_cuts.Add( elists[0] )
    print elists
    for e in elists: 
        print e.GetName(), e.GetN()
        all_cuts.Intersect(e)
    print all_cuts.GetName(), all_cuts.GetN()
    
    outFile.cd()
    [e.Write() for e in elists]
    all_cuts.Write()
    outFile.Close()


def test():
    """docstring for test"""
    datahists = [ DataHistos('alltrigs'), DataHistos('minbias'), DataHistos('hightower'),
                  DataHistos('jetpatch'), DataHistos('96201'), DataHistos('96211'), 
                  DataHistos('96221'), DataHistos('96233') ]
    
    reader = ROOT.StSpinTreeReader()
    reader.selectDataset('/Users/kocolosk/data/run5/dataset.txt')
    reader.connectNeutralJets = False
    reader.connectBemcElectrons = False
    for i in range(100):
        reader.GetEntry(i)
        [ h.fillEvent(reader.event()) for h in datahists ]
        for track in reader.chargedPions():
            [ h.fillTrack(track) for h in datahists ]
    
    outfile = ROOT.TFile('out.root','recreate')
    for d in datahists:
        [ h.Write() for h in d.values() ]
    outfile.Close()
    raw_input('The End')


if __name__ == '__main__':
    #makeDataHistograms()
    testDataHistos2()