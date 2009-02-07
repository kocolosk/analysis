name    = '_'.join(__name__.split('.')[-2:])
VERSION = '$Id$'[5:-2]

from array import array
import ROOT
import mcasym

class_ = ROOT.TH1D

binning = {
    'nbinsx': 5,
    'xbins': array('d', [2.0, 3.18, 4.56, 6.32, 8.8, 12.84])
}

props = {
    'SetXTitle': ('p_{T}',),
    'SetYTitle': (name,)
}

branches = ('mFlavor*', 'mX1', 'mParton1*', 'mParton2*', 'mParton3*',  
    'mProcessId', 'mMcVertex*', 'mMcTracks*')

def accept_event(event):
    simu_cut = isinstance(event, ROOT.StChargedPionMcEvent)
    vz_cut = abs(event.mcVertex().z()) < 60
    return simu_cut and vz_cut

def accept_track(track):
    etamc_cut = abs(track.etaMc()) < 1.0
    pid_cut = track.geantId() in (8,9)
    return etamc_cut and pid_cut

def analyze(event, **kw):
    for track in event.mcTracks():
        if event.charge_filter(track) and accept_track(track):
            yield (track.ptMc(), mcasym.num('MAX', event))

