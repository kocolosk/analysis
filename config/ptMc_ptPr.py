name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

from array import array
import ROOT

class_ = ROOT.TH2D

binning = {
    'nbinsx': 200,
    'xbins': (0.0, 20.0),
    'nbinsy': 5,
    'ybins': array('d', [2.0, 3.18, 4.56, 6.32, 8.8, 12.84])
}

props = {
    'SetXTitle': ('ptMc',),
    'SetYTitle': ('ptPr',)
}

branches = ('mVertices*', 'mMatchedPairs*')

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    simu_cut = isinstance(event, ROOT.StChargedPionMcEvent)
    return vertex_cut and simu_cut

def accept_track(track):
    pid_cut = track.geantId() in (8,9)
    return pid_cut

def analyze(event, **kw):
    for track in event.matchedPairs():
        if event.charge_filter(track) and accept_track(track):
            yield (track.ptMc(), track.ptPr())

