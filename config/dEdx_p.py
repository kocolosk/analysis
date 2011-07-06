name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

import ROOT
from analysis import pid

class_ = ROOT.TH2D

binning = {
    'nbinsx': 200,
    'xbins': (0.0, 20.0),
    'nbinsy': 200,
    'ybins': (1E-6, 6E-6)
}

props = {
    'SetXTitle': (name,)
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
    if isinstance(event, ROOT.StChargedPionMcEvent):
        pid_cut = True
    else:
        nsigpi = pid.shift(event.runId(), track.nSigmaPion())
        pt = track.Pt()
        if pt < 3.18:
            pid_cut = -1.1 < nsigpi < 2.3
        elif pt < 4.56:
            pid_cut = -1.4 < nsigpi < 2.1
        elif pt < 6.32:
            pid_cut = -1.4 < nsigpi < 1.8
        elif pt < 8.80:
            pid_cut = -1.4 < nsigpi < 1.8
        else:
            pid_cut = -1.3 < nsigpi < 1.4
    return eta_cut and dca_cut and fit_cut and pid_cut

def analyze(event, **kw):
    for track in event.tracks():
        if event.charge_filter(track) and accept_track(event, track):
            yield (track.P(), track.dEdx())

