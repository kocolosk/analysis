import ROOT

libs_to_load = [ 'libPhysics', 'libTable', 'StarRoot', 'StarClassLibrary', 'St_base', 
'StChain', 'St_Tables', 'StUtilities', 'StTreeMaker', 'StIOMaker', 'StTriggerDataMaker', 
'StBichsel', 'StEvent', 'StEventUtilities', 'StDbLib', 'StEmcUtil', 'StTofUtil', 'StPmdUtil', 
'StStrangeMuDstMaker', 'StMuDSTMaker', 'StDaqLib', 'StDetectorDbMaker', 'StEmcTriggerMaker', 
'StMCAsymMaker', 'StSpinDbMaker', 'St_db_Maker', 'StTriggerUtilities', 'StJetFinder', 'StJetMaker', 
'StChargedPionAnalysisMaker', 'StSpinTree', 'StMiniMcEvent']

print 'analysis : loading shared libraries ...'
libs_already_loaded = ROOT.gSystem.GetLibraries()
for library in libs_to_load:
   if library not in libs_already_loaded:
      ROOT.gSystem.Load(library)
      #print 'analysis : loaded', library
print 'analysis : loading complete'

## style stuff
ROOT.gStyle.SetCanvasColor(10)
ROOT.gStyle.SetFillColor(10)
ROOT.gStyle.SetStatColor(0)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetCanvasBorderMode(0)
ROOT.gStyle.SetOptDate(1)

import asym
import minimc
import histos
import tree
import datamc2 as datamc

## classes
from asym   import AsymmetryGenerator, ScalarCounts, Polarizations
from histos import HistogramManager
from minimc import MiniMcHistos

## utility methods
from util   import *

## runlists
from asym   import golden_runlist_c, minbias_runs, final_runlist_run5
from runlists import *

__all__ = ['asym','datamc2','histos','minimc']
