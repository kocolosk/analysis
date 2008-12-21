name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

import ROOT
from .. import pid

class_ = ROOT.TH1D

binning = {
    'nbinsx': 45,
    'xbins': (0.5, 45.5)
}

props = {
    'SetXTitle': (name,)
}

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    simu = isinstance(event, ROOT.StChargedPionMcEvent)
    bin = event.bbcTimeBin()/32
    bbc_cut = simu or bin in (7,8,9) or (event.runId() > 7000000 and bin==6)
    return vertex_cut and bbc_cut

def accept_track(event, track):
    eta_cut = abs( track.eta() ) < 1.0
    dca_cut = abs( track.globalDca().mag() ) < 1.0
    if isinstance(event, ROOT.StChargedPionMcEvent):
        pid_cut = True
    else:
        pid_min = pid.min(event.runId())
        pid_max = pid.max(event.runId())
        pid_cut = pid_min < track.nSigmaPion() < pid_max
    return eta_cut and dca_cut and pid_cut

def analyze(event, **kw):
    for track in event.tracks():
        if event.charge_filter(track) and accept_track(event, track):
            yield (track.nHitsFit(),)

