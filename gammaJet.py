import math
import ROOT
ROOT.gSystem.Load('StarSpinAnalyses')

def test():
    reader = ROOT.StSpinTreeReader()
    
    #reader.selectDataset('/Users/kocolosk/data/run6/dataset.txt')
    #reader.selectFile('/Users/kocolosk/data/run6/spinTree/spinAnalyses_7132007.tree.root')
    reader.selectDataset('/Users/kocolosk/Xcode/StarSpinLibraries/StRoot/StSpinPool/StSpinTree/datasets/run6_mit.dataset')
    reader.selectRun(7132007)
    
    reader.requireJet = True
    
    neutralPionJet = 0
    diJet = 0
    neutralDiJet = 0
    gammaJet = 0
    
    for i in xrange(reader.GetEntries()):
        reader.GetEntry(i)
        
        #jets
        if reader.nJets() != 2: continue
        jet1 = reader.jet(0)
        jet2 = reader.jet(1)
        dPhi = math.fabs(jet1.Phi() - jet2.Phi())
        if dPhi < (math.pi / 3): continue
        diJet += 1
        if jet1.tpcEtSum > 0 and jet2.tpcEtSum > 0: continue
        neutralDiJet += 1
        if jet1.tpcEtSum == 0: 
            neutralJet  = jet1
            otherJet    = jet2
        else:
            neutralJet  = jet2
            otherJet    = jet1
        
        #pi0s
        for i in range(reader.nBemcPions()):
            dEta = reader.bemcPion(i).eta() - neutralJet.Eta()
            dPhi = math.fabs(reader.bemcPion(i).phi() - neutralJet.Phi())
            if dPhi > 2 * math.pi: dPhi -= 2 * math.pi
            dR = math.sqrt(dEta*dEta + dPhi*dPhi)
            if dR < 0.1:
                neutralPionJet += 1
                continue
            
        #maybe it's a photon-gamma!
        for i in range(reader.nGammas()):
            dEta = reader.gamma(i).Eta() - neutralJet.Eta()
            dPhi = math.fabs(reader.gamma(i).Phi() - neutralJet.Phi())
            if dPhi > 2 * math.pi: dPhi -= 2 * math.pi
            dR = math.sqrt(dEta*dEta + dPhi*dPhi)
            if dR < 0.1:
                gammaJet += 1
    
    print 'total =', reader.GetEntries()
    print 'diJet =', diJet
    print 'neutralDiJet =', neutralDiJet
    print 'neutralPionJet =', neutralPionJet
    print 'gammaJet =', gammaJet