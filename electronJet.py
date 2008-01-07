#/usr/bin/env python

import math,sys
import ROOT
ROOT.gSystem.Load('StarSpinAnalyses')
#import mystyle; mystyle.use(1)
ROOT.gStyle.SetOptLogy(1)
ROOT.gStyle.SetOptStat(0)

def initReaderRun6(reader):
    #reader.selectDataset('/Users/kocolosk/data/run6/dataset.txt')
    #reader.selectDataset('/Users/kocolosk/data/run6/dataset_dg5.txt')
    reader.selectFile('/Users/kocolosk/new/spinAnalyses_7132007.tree.root')
    
    reader.connectJets          = True
    reader.connectNeutralJets   = False
    reader.connectChargedPions  = False
    reader.connectBemcPions     = False
    reader.connectEemcPions     = False
    reader.connectElectrons     = True
    
    #reader.selectRunlist('/Users/kocolosk/data/run6/runlists/713.runlist')
    #reader.selectRunlist('/Users/kocolosk/Xcode/StarSpinLibraries/StRoot/StSpinPool/StSpinTree/filters/run6_jets.runlist')
    #reader.selectRunlist('/Users/kocolosk/Xcode/StarSpinLibraries/StRoot/StSpinPool/StSpinTree/filters/run6_electrons.runlist')
    #reader.selectRun(7132007)
    
    #reader.selectTrigger(137221)
    #reader.selectTrigger(137222)
    #reader.selectTrigger(137611)
    
    #reader.requireDidFire       = True
    reader.requireShouldFire    = False

def sanity():
    reader = ROOT.StSpinTreeReader()
    reader.selectDataset('/Users/kocolosk/data/run6/dataset_dg5.txt')
    
    reader.connectJets          = False
    reader.connectNeutralJets   = False
    reader.connectChargedPions  = False
    reader.connectBemcPions     = False
    reader.connectEemcPions     = False
    reader.connectElectrons     = True
    
    reader.selectRun(7132068)
    
    for i in range(reader.GetEntries()):
        reader.GetEntry(i)
        
        for electron in reader.bemcElectrons():
            meth = getattr(electron,'global')
            gl = meth()
            if gl is None: 
                print 'global TRef broken'
                continue
            print electron.pt(), electron.eta(), electron.phi()
            print gl.pt(), gl.eta(), gl.phi()
            print '---------------------------------------------'


def make():
    reader = ROOT.StSpinTreeReader()
    initReaderRun6(reader)
    
    h_inclusive = ROOT.TH1D('h_inclusive','',100,0.,1.)
    h_electron  = ROOT.TH1D('h_electron','',100,0.,1.)
    
    track_mult_inc = ROOT.TH1D('track_mult_inc','',20,-0.5,19.5)
    track_mult_ele = ROOT.TH1D('track_mult_ele','',20,-0.5,19.5)
    
    ptComp = ROOT.TH2D('ptComp','',100,0,50,100,0,20)
    
    for i in range(reader.GetEntries()):
        reader.GetEntry(i)
        
        for jet in reader.jets():
            R = (jet.btowEtSum + jet.etowEtSum) / jet.Et()
            if R == 0. or R > 0.99: continue
            #if jet.geomTrigger(137221) == False and jet.geomTrigger(137222) == False: continue
            if jet.detEta() > 0.9 or jet.detEta() < -0.7: continue
            
            h_inclusive.Fill(R)
            track_mult_inc.Fill(jet.nTracks)
            for electron in reader.bemcElectrons():
                deta = math.fabs(electron.eta() - jet.Eta())
                dphi = math.fabs(electron.phi() - jet.Phi())
                if dphi > math.pi: dphi = math.fabs(dphi - 2*math.pi)
                dR = math.sqrt(deta**2 + dphi**2)
                if dR < 0.7:
                    h_electron.Fill(R)
                    track_mult_ele.Fill(jet.nTracks)
                    ptComp.Fill(jet.Pt(),electron.pt())
                    
    h_inclusive.SetXTitle('R_{T}')
    h_inclusive.SetTitle('L2GammaB')
    h_electron.SetLineColor(ROOT.kRed)
    
    c = ROOT.TCanvas()    
    h_inclusive.Draw()
    h_electron.Draw('same')
    
    leg = ROOT.TLegend(0.15,0.73,0.3,0.87)
    leg.AddEntry(h_inclusive,'inclusive jets','l')
    leg.AddEntry(h_electron,'electron jets','l')
    leg.Draw()
    
    c2 = ROOT.TCanvas('c2')
    track_mult_inc.SetTitle('L2GammaB')
    track_mult_inc.SetXTitle('track multiplicity')
    track_mult_inc.Draw()
    track_mult_ele.SetLineColor(ROOT.kRed)
    track_mult_ele.Draw('same')
    leg.Draw()
    
    c3 = ROOT.TCanvas('c3')
    ptComp.SetTitle('Electron jets in L2GammaB')
    ptComp.SetYTitle('electron p_{T}')
    ptComp.SetXTitle('jet p_{T}')
    ptComp.Draw('col')
    c3.SetLogy(0)
    
    print 'inclusive jets = ', int(h_inclusive.GetEntries())
    print 'electron jets  = ', int(h_electron.GetEntries())
    raw_input('press enter to continue:')
    
    
if __name__ == '__main__':
    sys.exit(make())