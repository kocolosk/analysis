import ROOT

import socket
if socket.gethostname().endswith('bnl.gov'):
    libs_to_load = [ 
    'libPhysics', 'libTable', 'StarRoot', 'StarClassLibrary', 'St_base',
    'StChain', 'St_Tables', 'StUtilities', 'StTreeMaker', 'StIOMaker',
    'StTriggerDataMaker', 'StBichsel', 'StEvent', 'StEventUtilities', 'StDbLib',
    'StEmcUtil', 'StTofUtil', 'StPmdUtil', 'StStrangeMuDstMaker',
    'StMuDSTMaker', 'StDaqLib', 'StDetectorDbMaker', 'StEmcTriggerMaker',
    'StJetSkimEvent', 'StJets', 'StMCAsymMaker', 'StSpinDbMaker', 'St_db_Maker',
    'StTriggerUtilities', 'StEEmcUtil', 'StEmcRawMaker', 'StEmcADCtoEMaker',
    'StJetFinder', 'StJetMaker', 'StMiniMcEvent', 'StChargedPionAnalysisMaker',
    'StSpinTree'
    ]
else:
    libs_to_load = [ 'StChargedPionEvent' ]
del socket

print 'analysis : loading shared libraries ...'
libs_already_loaded = ROOT.gSystem.GetLibraries()
for library in libs_to_load:
   if library not in libs_already_loaded:
      ROOT.gSystem.Load(library)
      #print 'analysis : loaded', library
print 'analysis : loading complete'

## this is nifty ... can extend classes on-the-fly
ROOT.StTinyMcTrack.charge = ROOT.StTinyMcTrack.chargeMc
ROOT.StMiniMcPair.pt = ROOT.StMiniMcPair.ptPr
ROOT.StMiniMcPair.eta = ROOT.StMiniMcPair.etaPr
ROOT.StMiniMcPair.nHitsFit = ROOT.StMiniMcPair.fitPts
def _minimc_globalDca(self):
    try:
        self.dcaVec.set(0, 0, self.dcaGl())
    except AttributeError:
        self.dcaVec = ROOT.StThreeVectorF(0, 0, self.dcaGl())
    return self.dcaVec
ROOT.StMiniMcPair.globalDca = _minimc_globalDca
def _minimc_DeltaPhi(self, jet):
    try:
        self.lv.SetPtEtaPhiM(self.ptPr(), self.etaPr(), self.phiPr(), 0.135)
    except AttributeError:
        self.lv = ROOT.TLorentzVector()
        self.lv.SetPtEtaPhiM(self.ptPr(), self.etaPr(), self.phiPr(), 0.135)
    return self.lv.DeltaPhi(jet)
ROOT.StMiniMcPair.DeltaPhi = _minimc_DeltaPhi

## style stuff
ROOT.gStyle.SetCanvasColor(10)
ROOT.gStyle.SetFillColor(10)
ROOT.gStyle.SetStatColor(0)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetCanvasBorderMode(0)
ROOT.gStyle.SetOptDate(1)

## runlists
from runlists import *

import asym
import minimc
import histos
import tree
import simu
import ff
import histos2
import config

## spin2008 module uses relative imports
from sys import version_info
if version_info[1] >= 4:
    import spin2008

## classes
from asym   import AsymmetryGenerator, ScalarCounts, Polarizations
from histos import HistogramManager
from minimc import MiniMcHistos
from xsec   import datapoint as DataPoint

## utility methods
from util   import *

import plots
from plots import spin2008
from plots import graphics

__all__ = ['asym','datamc2','histos','minimc', 'runlists', 'util']
