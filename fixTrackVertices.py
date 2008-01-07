#!/usr/bin/env python
# encoding: utf-8
"""
fixTrackVertices.py

Created by Adam Kocoloski on 2007-05-23.
Copyright (c) 2007 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import glob
import re
import paths
import ROOT

ROOT.gSystem.Load('StarSpinAnalyses')

def fixTrackVertices(runnumber):
    jetSkimFileName = '%s%d.tree.root' % (paths.path('jetSkim','MIT',5),runnumber)
    jetSkimFile = ROOT.TFile(jetSkimFileName)
    jetSkimTree = jetSkimFile.Get('jetSkimTree')
    
    chargedPionFileName = '%s%d.tree.root' % (paths.path('chargedPions','MIT',5),runnumber)
    chargedPionFile = ROOT.TFile(chargedPionFileName)
    chargedPionTree = chargedPionFile.Get('chargedPionTree')
    
    ev = ROOT.StJetSkimEvent()
    jetSkimTree.SetBranchAddress('skimEventBranch',ev)
    
    tracks = ROOT.TClonesArray('StChargedPionTrack',100)
    chargedPionTree.SetBranchAddress('primaries',tracks)
    
    outputFileName = 'chargedPions_%d.tree.root' % (runnumber,)
    outputFile = ROOT.TFile(outputFileName,'recreate')
    outputTree = chargedPionTree.CloneTree(0)
    
    for i in xrange(jetSkimTree.GetEntries()):
        jetSkimTree.GetEntry(i)
        chargedPionTree.GetEntry(i)
        
        trackVertices = [ROOT.StThreeVectorF(v.position()[0],v.position()[1],v.position()[2]) for v in list(ev.vertices())]
        for track in list(tracks):
            track.setVertex(trackVertices[track.vertexIndex()])
        outputTree.Fill()
    
    badTrackHisto = chargedPionFile.Get('badTracks')
    
    outputFile.cd()
    outputTree.Write()
    if badTrackHisto is not None: badTrackHisto.Write()
    outputFile.Close()
    print 'wrote file ' + outputFileName

def main():
    paths = glob.glob('/Volumes/scratch/common/run5/chargedPions/chargedPions_6*')
    for path in paths:
        runnumber = re.search('(6|7)\d{6}',path).group()
        fixTrackVertices(int(runnumber))
    


if __name__ == '__main__':
    main()

