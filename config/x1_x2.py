name    = __name__.split('.')[-1]
VERSION = '$Id$'.strip('$Id: ')

import ROOT

class_ = ROOT.TH2D

binning = {
    'nbinsx': 200,
    'xbins': (0.0, 1.0),
    'nbinsy': 200,
    'ybins': (0.0, 1.0)
}

props = {
    'SetXTitle': ('x1',),
    'SetYTitle': ('x2',)
}

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    simu_cut = isinstance(event, ROOT.StChargedPionMcEvent)
    return vertex_cut and simu_cut

def analyze(event):
    yield (event.x1(), event.x2())

