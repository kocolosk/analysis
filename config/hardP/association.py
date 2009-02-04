name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

from array import array
import ROOT

class_ = ROOT.TH2D

binning = {
    'nbinsx': 5,
    'xbins': array('d', [2.0, 3.18, 4.56, 6.32, 8.8, 12.84]),
    'nbinsy': 200,
    'ybins': (0.0, 50.0)
}

props = {
    'SetXTitle': ('#pi p_{T}',),
    'SetYTitle': ('partonic p_{T}',)
}

branches = ('mMcVertex*', 'mVertices*', 'mMatchedPairs*', 'mHardP')

def accept_event(event):
    simu_cut = isinstance(event, ROOT.StChargedPionMcEvent)
    vz_cut = abs(event.mcVertex().z()) < 60
    vertex_cut = event.nVertices() > 0
    return simu_cut and vz_cut and vertex_cut

def accept_track(track):
    return track.geantId() in (8,9)

def analyze(event, **kw):
    for track in event.matchedPairs():
        if event.charge_filter(track) and accept_track(track):
            yield (track.ptMc(), event.hardP())

