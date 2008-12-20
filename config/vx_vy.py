name    = __name__.split('.')[-1]
VERSION = '$Id$'[5:-2]

import ROOT

class_ = ROOT.TH2D

binning = {
    'nbinsx': 400,
    'xbins': (-2.0, 2.0),
    'nbinsy': 400,
    'ybins': (-2.0, 2.0)
}

props = {
    'SetXTitle': ('vx',),
    'SetYTitle': ('vy',)
}

def accept_event(event):
    vertex_cut = event.nVertices() > 0
    simu = isinstance(event, ROOT.StChargedPionMcEvent)
    bin = event.bbcTimeBin()/32
    bbc_cut = simu or bin in (7,8,9) or (event.runId() > 7000000 and bin==6)
    return vertex_cut and bbc_cut

def analyze(event, **kw):
    yield (event.vertex(0).x(), event.vertex(0).y())

