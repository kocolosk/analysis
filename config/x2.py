name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

import ROOT

class_ = ROOT.TH1D

binning = {
    'nbinsx': 200,
    'xbins': (0.0, 1.0)
}

props = {
    'SetXTitle': ('partonic p_{T}',)
}

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    simu_cut = isinstance(event, ROOT.StChargedPionMcEvent)
    return vertex_cut and simu_cut

def analyze(event):
    yield (event.x2(),)

