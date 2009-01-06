name    = '_'.join(__name__.split('.')[-2:])
VERSION = '$Id$'[5:-2]

import ROOT
import mcasym

class_ = ROOT.TH1D

binning = {
    'nbinsx': 40,
    'xbins': (0.0, 20.0)
}

props = {
    'SetXTitle': ('p_{T}',),
    'SetYTitle': (name,)
}

branches = ('mVertices*', 'mMatchedPairs*')

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    simu_cut = isinstance(event, ROOT.StChargedPionMcEvent)
    return vertex_cut and simu_cut

def accept_track(track):
    eta_cut = abs( track.eta() ) < 1.0
    dca_cut = abs( track.globalDca().mag() ) < 1.0
    fit_cut = track.nHitsFit() > 25
    pid_cut = track.geantId() in (8,9)
    return eta_cut and dca_cut and fit_cut and pid_cut

def analyze(event, **kw):
    for track in event.matchedPairs():
        if event.charge_filter(t) and accept_track(t):
            yield (track.pt(), mcasym.num('STD', event))
