name    = __name__.split('.')[-1]
VERSION = '$Id$'.strip('$Id: ')

import ROOT

class_ = ROOT.TH1D

binning = {
    'nbinsx': 15,
    'xbins': (-0.5, 14.5)
}

props = {
    'SetXTitle': (name,)
}

def accept_event(event):
    return True

def analyze(event):
    yield (event.nVertices(),)

