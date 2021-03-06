name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

import ROOT

class_ = ROOT.TH1D

binning = {
    'nbinsx': 15,
    'xbins': (-0.5, 14.5)
}

props = {
    'SetXTitle': (name,)
}

branches = ('mVertices*', )

def accept_event(event):
    return True

def analyze(event, **kw):
    yield (event.nVertices(),)

