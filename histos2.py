import ROOT
from array import array

class Histo(object):
    """
    wrapper around ROOT histograms.  It uses a config module to determine how
    to instantiate and fill the histogram, and it overloads the histogram's
    Fill method to treat multi-particle statistics correctly.
    
    spin and charge are needed only to generate a unique name for ROOT.  The 
    actual filtering is done elsewhere.
    """
    def __init__(self, trig, spin, mod, tfile=None, charge=None):
        self.mod = mod
        self.trigger_filter = trigger_filter(trig)
        if charge:
            name = '_%s_%s_%s_%s' % (trig, spin, charge, mod.name)
        else:
            name = '_%s_%s_%s' % (trig, spin, mod.name)
        title = '%s@%s' % (mod.name, mod.VERSION)
        if tfile:
            self.h = tfile.Get(name)
        else:
            self.__construct(name, title, **mod.binning)
        [ getattr(self.h, key)(*val) for key,val in mod.props.items() ]
        self.profile = (mod.class_ == ROOT.TProfile)
        if not self.profile:
            self.h.Sumw2()
        self.vals = []
    
    def __getattr__(self, name):
        """
        fall back to ROOT method if I didn't define a replacement
        """
        return getattr(self.h, name)
    
    def __repr__(self):
        return '<Histo(%s) "%s" at %#x>' % (self.mod.class_.__name__, \
            self.h.GetName(), id(self))
    
    @staticmethod
    def __checkBins(nbins, bins):
        if isinstance(bins, array):
            assert nbins == len(bins)-1
        elif bins:
            assert len(bins) == 2
        else:
            assert nbins == 0
    
    def __construct(self, name, title, nbinsx, xbins, nbinsy=0, nbinsz=0, **kw):
        ybins = kw.get('ybins')
        zbins = kw.get('zbins')
        self.__checkBins(nbinsx, xbins)
        self.__checkBins(nbinsy, ybins)
        self.__checkBins(nbinsz, zbins)
        class_ = self.mod.class_
        if isinstance(zbins, array):
            binning = (nbinsx, xbins, nbinsy, ybins, nbinsz, zbins)
        elif zbins:
            binning = (nbinsx, xbins[0], xbins[1], nbinsy, ybins[0], ybins[1], \
                nbinsz, zbins[0], zbins[1])
        elif isinstance(xbins, array) and isinstance(ybins, array):
            binning = (nbinsx, xbins, nbinsy, ybins)
        elif isinstance(xbins, array):
            binning = (nbinsx, xbins, nbinsy, ybins[0], ybins[1])
        elif isinstance(ybins, array):
            binning = (nbinsx, xbins[0], xbins[1], nbinsy, ybins)
        elif ybins:
            binning = (nbinsx, xbins[0], xbins[1], nbinsy, ybins[0], ybins[1])
        elif isinstance(xbins, array):
            binning = (nbinsx, xbins)
        else:
            binning = (nbinsx, xbins[0], xbins[1])
        self.h = class_(name, title, *binning)
    
    def Fill(self, x, y=None, z=None):
        if self.profile:
            self.h.Fill(x,y)
        else:
            self.vals.append((x,y,z))
    
    def Flush(self):
        """
        Ends the event and actually fills the ROOT histogram, taking multi-
        particle statistics into account.
        """
        if self.profile: return
        weight = {}
        keep = {}
        for x,y,z in self.vals:
            bin = self.h.FindBin(x, y or 0, z or 0)
            weight[bin] = weight.get(bin,0) + 1
            keep[bin] = (x,y,z)
        for bin,w in weight.items():
            x,y,z = keep[bin]
            if z:
                self.h.Fill(x,y,z,w)
            elif y:
                try:
                    self.h.Fill(x,y,w)
                except TypeError:
                    # assume this is an mcasym histo which has its own weight
                    self.h.Fill(x,y)
            else:
                self.h.Fill(x,w)
        self.vals = []
    
    def analyze(self, event):
        if self.trigger_filter(event) and self.mod.accept_event(event):
            [self.Fill(*vals) for vals in self.mod.analyze(event)]
        self.Flush()
    


trigger_cache = {}
def passed(ev, trigId):
    try:
        return trigger_cache[trigId]
    except KeyError:
        val = ev.isSimuTrigger(trigId) and \
            (ev.isTrigger(trigId) or isinstance(ev, ROOT.StChargedPionMcEvent))
        trigger_cache[trigId] = val
        return val


def trigger_filter(trigId):
    fun = lambda ev: False
    try:
        itrig = int(trigId)
        fun = lambda ev: passed(ev,itrig)
    except ValueError:
        if trigId == 'jetpatch':
            fun = lambda ev: passed(ev,137222) or passed(ev,137221) or \
                passed(ev,96233) or passed(ev,96221)
    return fun


def update(modlist, triggers, tree, tfile=None):
    spinKeys = {5:'uu', 6:'du', 9:'ud', 10:'dd'}
    subProcessKeys = {68:'gg', 28:'qg', 11:'qq'}
    trackHistos = ['pt','eta','phi','nHitsFit','dcaG','dEdx','nSigmaPion']
    global trigger_cache
    
    tree.GetEntry(0)
    simu = isinstance(tree.event, ROOT.StChargedPionMcEvent)
    spinlist = ['other']
    spinlist.extend(simu and ('gg','qg','qq') or ('uu','ud','du','dd'))
    
    hevent  = dict.fromkeys(spinlist, [])
    hsum    = dict.fromkeys(spinlist, [])
    hplus   = dict.fromkeys(spinlist, [])
    hminus  = dict.fromkeys(spinlist, [])
    for spin in spinlist:
        for trig in triggers:
            for mod in filter(lambda m: m.name not in trackHistos, modlist):
                hevent[spin].append(Histo(trig, spin, mod))
            for mod in filter(lambda m: m.name in trackHistos, modlist):
                hsum[spin].append(Histo(trig, spin, mod, charge='sum'))                
                hplus[spin].append(Histo(trig, spin, mod, charge='plus'))
                hminus[spin].append(Histo(trig, spin, mod, charge='minus'))
    
    for entry in tree:
        ev = tree.event
        
        spin = 'other'
        if simu:
            spin = subProcessKeys.get(ev.processId(), 'other')
        elif ev.isSpinValid():
            spin = spinKeys.get(ev.spinBit(), 'other')
        
        [ h.analyze(ev) for h in hevent[spin] ]
        
        ev.charge_filter = lambda t: True
        [ h.analyze(ev) for h in hsum[spin] ]
        
        ev.charge_filter = lambda t: t.charge() == 1
        [ h.analyze(ev) for h in hplus[spin] ]
        
        ev.charge_filter = lambda t: t.charge() == -1
        [ h.analyze(ev) for h in hminus[spin] ]
        
        trigger_cache = {}
    
    ## don't forget post-processing for 'anyspin' histos
    for coll in (hevent, hsum, hplus, hminus):
        coll['anyspin'] = [ h.Clone() for h in coll['other'] ]
        for i in range(len(coll['anyspin'])):
            for spin in spinlist[1:]:
                coll['anyspin'][i].Add(coll[spin][i].h)
    
    if tfile:
        tfile.cd()
        for coll in (hevent, hsum, hplus, hminus):
            for hlist in coll.values():
                [ h.Write() for h in hlist ]
    return {'event':hevent, 'plus':hplus, 'minus':hminus, 'sum':hsum}


def write_histograms(treeDir='~/data/run5/tree', globber='*', **kw):
    ## FIXME: this won't work with simulations b/c it makes one .root/file
    from os.path import basename
    from glob import glob
    from . import config
    
    modlist = kw.get('modlist') or [getattr(config, m) for m in config.modules]
    
    triggers = kw.get('trigList') or \
        ('96011','96221','96233','117001','137221','137222','jetpatch')
    
    files = glob(treeDir + '/' + globber + '.root')
    for fname in files:
        print fname
        f = ROOT.TFile(fname)
        f.tree.GetEntry(0)
        simu = isinstance(f.tree.event, ROOT.StChargedPionMcEvent)
        if simu:
            outname = basename(fname)[:9] + '.cphist.root'
        else:
            outname = basename(fname).replace('.tree.','.hist.')
        
        outFile = ROOT.TFile(outname, 'update')
        update(modlist, triggers, f.tree, outFile)
        outFile.Close()

