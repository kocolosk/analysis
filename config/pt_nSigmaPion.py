name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

import ROOT
from array import array
from analysis import pid

class_ = ROOT.TH2D

binning = {
    'nbinsx': 5,
    'xbins': array('d', [2.0, 3.18, 4.56, 6.32, 8.8, 12.84]),
    'nbinsy': 240,
    'ybins': (-6.0, 6.0)
}

props = {
    'SetXTitle': ('p_{T}',),
    'SetYTitle': ('n#sigma(#pi)',)
}

branches = ('mVertices*', 'mBbcTimeBin', 'mRunId', 'mTracks*')

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    simu = isinstance(event, ROOT.StChargedPionMcEvent)
    bin = event.bbcTimeBin()/32
    bbc_cut = simu or bin in (7,8,9) or (event.runId() > 7000000 and bin==6)
    return vertex_cut and bbc_cut

def accept_track(event, track):
    eta_cut = abs( track.eta() ) < 1.0
    dca_cut = abs( track.globalDca().mag() ) < 1.0
    fit_cut = track.nHitsFit() > 25
    return eta_cut and dca_cut and fit_cut

def analyze(event, jet_trigger_filter, **kw):
    for track in event.tracks():
        if event.charge_filter(track) and accept_track(event, track):
            yield (track.Pt(), track.nSigmaPion())

