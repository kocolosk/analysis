name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

import ROOT

class_ = ROOT.TH1D

binning = {
    'nbinsx': 500,
    'xbins': (-150., 150.)
}

props = {
    'SetXTitle': (name,)
}

branches = ('mVertices*', )

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    return vertex_cut

def analyze(event, **kw):
    yield (event.vertex(0).z(),)

