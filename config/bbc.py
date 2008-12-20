name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

import ROOT

class_ = ROOT.TH1D

binning = {
    'nbinsx': 400,
    'xbins': (-0.5, 399.5)
}

props = {
    'SetXTitle': (name,)
}

def accept_event(event):
    return True

def analyze(event):
    yield (event.bbcTimeBin(),)
