name    = __name__.split('.')[-1]
VERSION = '$Id$'.strip('$Id: ')

import ROOT

class_ = ROOT.TH1D

binning = {
    'nbinsx': 500,
    'xbins': (-150., 150.)
}

props = {
    'SetXTitle': (name,)
}

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    return vertex_cut

def analyze(event):
    yield (event.vertex(0).z(),)

