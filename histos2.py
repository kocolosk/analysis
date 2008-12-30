import ROOT
import math
from array import array
from sets import Set

class Histo(object):
    """
    wrapper to write ROOT histograms.  It uses a config module to determine how
    to instantiate and fill the histogram, and it overloads the histogram's
    Fill method to treat multi-particle statistics correctly.
    
    spin and charge are needed only to generate a unique name for ROOT.  The 
    actual filtering is done elsewhere.
    """
    def __init__(self, trig, spin, mod, charge=None):
        self.mod = mod
        self.trigger_filter = trigger_filter(trig)
        self.jet_filter = jet_trigger_filter(trig)
        if charge:
            name = '_%s_%s_%s_%s' % (trig, spin, charge, mod.name)
        else:
            name = '_%s_%s_%s' % (trig, spin, mod.name)
        title = '%s@%s' % (mod.name, mod.VERSION)
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
    
    def __checkBins(nbins, bins):
        if isinstance(bins, array):
            assert nbins == len(bins)-1
        elif bins:
            assert len(bins) == 2
        else:
            assert nbins == 0
    __checkBins = staticmethod(__checkBins)
    
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
        elif ybins and isinstance(xbins, array):
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
            [self.Fill(*vals) for vals in self.mod.analyze(event, \
                jet_trigger_filter=self.jet_filter)]
        self.Flush()
    


class HColl(dict):
    """
    Just a dict whose keys are accessible as attributes.  Of limited utility
    since integers can't be attributes, so all specific triggers are 
    inaccessible.
    """
    def __getattr__(self, name):
        try:
            val = self[name]
        except KeyError:
            raise AttributeError
        setattr(self, name, val)
        return val
    
    def trackHistograms(self, charge):
        """
        will raise AttributeError if called on a non-track-level collection
        """
        if charge == 1:  return self.tracks_plus
        if charge == -1: return self.tracks_minus
        if charge == 0:  return self.tracks_sum


class HKey(object):
    """
    allows for lazy loading of objects from a TFile in HistogramManager
    """
    def __init__(self, key):
        self.key = key
        self.obj = None
    
    def __getattr__(self, name):
        if not self.obj:
            self.obj = self.key.ReadObj()
        return getattr(self.obj, name)
    


class HistogramManager(HColl):
    """
    loads TKeys from the TFile into a deep dictionary
    """
    def __init__(self, tfile=None, **kw):
        super(HistogramManager, self).__init__()
        self.tfile = tfile
        
        keyparts = [ k.GetName().split('_')[1:] for k in tfile.GetListOfKeys() ]
        triggers = Set([k[0] for k in keyparts])
        spins = Set([k[1] for k in keyparts])
        
        for s in spins:
            self[s] = HColl.fromkeys(triggers, HColl())
            for t in triggers:
                for charge in ('plus', 'minus', 'sum'):
                    self[s][t][charge] = HColl()
                ## for backwards compatibility
                self[s][t].tracks_plus = self[s][t]['plus']
                self[s][t].tracks_minus = self[s][t]['minus']
                self[s][t].tracks_sum = self[s][t]['sum']
        
        for key in tfile.GetListOfKeys():
            k = key.GetName().split('_')[1:]
            if len(k) == 3:
                self[k[1]][k[0]][k[2]] = HKey(key)
            elif k[2] in ('plus', 'minus', 'sum'):
                self[k[1]][k[0]][k[2]]['_'.join(k[3:])] = HKey(key)
            else:
                self[k[1]][k[0]]['_'.join(k[2:])] = HKey(key)
    


trigger_cache = {}
def passed(ev, trigId):
    try:
        return trigger_cache[trigId]
    except KeyError:
        val = ev.isSimuTrigger(trigId) and \
            (isinstance(ev, ROOT.StChargedPionMcEvent) or ev.isTrigger(trigId))
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


def jet_passed(event, jet, trigId):
    if trigId in (96011, 117001): 
        return True
    elif trigId == 96201:
        for particle in jet.particles():
            if particle.detectorId() == ROOT.kBarrelEmcTowerId:
                if event.highTowerAdc(particle.index()) > 13:
                    return True
        return False
    elif trigId == 96211:
        for particle in jet.particles():
            if particle.detectorId() == ROOT.kBarrelEmcTowerId:
                if event.highTowerAdc(particle.index()) > 17:
                    return True
        return False
    elif trigId == 96221:
        patchPhi = (90.,30.,-30.,-90.,-150.,150.)
        for patchId in range(6):
            if event.jetPatchAdc(patchId) > 66:
                dPhi = abs(math.degrees(jet.Phi()) - patchPhi[patchId])
                if dPhi < 40 or dPhi > 320:
                    return True
        return False
    elif trigId == 96233:
        patchPhi = (90.,30.,-30.,-90.,-150.,150.)
        for patchId in range(6):
            if event.jetPatchAdc(patchId) > 83:
                dPhi = abs(math.degrees(jet.Phi()) - patchPhi[patchId])
                if dPhi < 40 or dPhi > 320:
                    return True
        return False
    elif trigId == 137221:
        patchPhi = (150.,90.,30.,-30.,-90.,-150.,150.,90.,30.,-30.,-90.,-150.)
        for patchId in range(12):
            if event.jetPatchAdc(patchId) > 58:
                dPhi = abs(math.degrees(jet.Phi()) - patchPhi[patchId])
                if dPhi < 36 or dPhi > 324:
                    return True
        return False
    elif trigId == 137222:
        patchPhi = (150.,90.,30.,-30.,-90.,-150.,150.,90.,30.,-30.,-90.,-150.)
        for patchId in range(12):
            if event.jetPatchAdc(patchId) > 60:
                dPhi = abs(math.degrees(jet.Phi()) - patchPhi[patchId])
                if dPhi < 36 or dPhi > 324:
                    return True


def jet_trigger_filter(trigId):
    fun = lambda ev,jet: False
    try:
        itrig = int(trigId)
        fun = lambda ev,jet: jet_passed(ev,jet,itrig)
    except ValueError:
        if trigId == 'jetpatch':
            fun = lambda ev,jet: \
                (passed(ev, 96221)  and jet_passed(ev, jet, 96221)  ) or \
                (passed(ev, 96233)  and jet_passed(ev, jet, 96233)  ) or \
                (passed(ev, 137221) and jet_passed(ev, jet, 137221) ) or \
                (passed(ev, 137222) and jet_passed(ev, jet, 137222) )
    return fun


def update(modlist, triggers, tree, tfile=None):
    spinKeys = {5:'uu', 6:'du', 9:'ud', 10:'dd'}
    subProcessKeys = {68:'gg', 28:'qg', 11:'qq'}
    trackHistos = [ m.name for m in 
        filter(lambda m: hasattr(m, 'accept_track'), modlist) ]
    global trigger_cache
    
    tree.GetEntry(0)
    simu = isinstance(tree.event, ROOT.StChargedPionMcEvent)
    
    if not simu:
        year = tree.event.runId() > 7000000 and 2006 or 2005
    else:
        ## muDstName looks like '/rcf1302_01_2000evts.MuDst.root'
        dataset = int(tree.event.muDstName()[4:8])
        year = (dataset > 1300) and 2006 or 2005
    
    spinlist = ['other']
    spinlist.extend(simu and ('gg','qg','qq') or ('uu','ud','du','dd'))
    
    hevent  = dict.fromkeys(spinlist)
    hsum    = dict.fromkeys(spinlist)
    hplus   = dict.fromkeys(spinlist)
    hminus  = dict.fromkeys(spinlist)
    for spin in spinlist:
        hevent[spin]    = list()
        hsum[spin]      = list()
        hplus[spin]     = list()
        hminus[spin]    = list()
        for trig in triggers:
            for mod in filter(lambda m: m.name not in trackHistos, modlist):
                hevent[spin].append(Histo(trig, spin, mod))
            for mod in filter(lambda m: m.name in trackHistos, modlist):
                hsum[spin].append(Histo(trig, spin, mod, charge='sum'))
                hplus[spin].append(Histo(trig, spin, mod, charge='plus'))
                hminus[spin].append(Histo(trig, spin, mod, charge='minus'))
    
    ## only default branches are the ones used by the trigger filter
    tree.SetBranchStatus('*', 0)
    tree.SetBranchStatus('mSimuTriggerBits', 1)
    if tree.GetBranch('mTriggerBits'):
        tree.SetBranchStatus('mTriggerBits', 1)
    
    for mod in modlist:
        for branchname in mod.branches:
            tree.SetBranchStatus(branchname, 1)
    
    fname = ''
    for entry in tree:
        if fname != tree.GetCurrentFile().GetName():
            fname = tree.GetCurrentFile().GetName()
            print 'starting', fname
        
        ev = tree.event
        ev.year = year
        
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
        [ h.SetName( h.GetName().replace('other','anyspin') ) for h in \
            coll['anyspin']]
        for i in range(len(coll['anyspin'])):
            for spin in spinlist[1:]:
                coll['anyspin'][i].Add(coll[spin][i].h)
    
    if tfile:
        tfile.cd()
        for coll in (hevent, hsum, hplus, hminus):
            for hlist in coll.values():
                [h.Write() for h in filter(lambda h: h.GetEntries() > 0, hlist)]
    return {'event':hevent, 'plus':hplus, 'minus':hminus, 'sum':hsum}


def write_histograms(treeDir='~/data/run5/tree', globber='*', **kw):
    from os.path import basename, join
    from glob import glob
    from analysis import config
    
    def all_modules(mod):
        """
        return a list of all histogram modules in config
        """
        ret = list()
        if hasattr(mod, 'modules'):
            [ ret.extend(all_modules(getattr(mod, m))) for m in mod.modules ]
        else:
            ret.append(mod)
        return ret
    
    modlist = kw.get('modlist') or all_modules(config)
    
    triggers = kw.get('trigList') or \
        ('96011','96221','96233','117001','137221','137222','jetpatch')
    
    histdir = kw.get('histdir') or './'
    
    files = glob(treeDir + '/' + globber + '.root')
    for fname in files:
        f = ROOT.TFile(fname)
        f.tree.GetEntry(0)
        simu = isinstance(f.tree.event, ROOT.StChargedPionMcEvent)
        if simu:
            outname = join(histdir, basename(fname)[:7] + '.cphist.root')
            tree = ROOT.TChain('tree')
            tree.Add(treeDir + '/' + globber + '.root')
        else:
            outname = join(histdir, basename(fname).replace('.tree.','.hist.'))
            tree = f.tree
        
        outFile = ROOT.TFile(outname, 'update')
        update(modlist, triggers, tree, outFile)
        outFile.Close()
        if simu: break

def condor_simu(treeDir, trigList=None, modlist=None, histdir='./'):
    """
    submits a single writeHistograms job to Condor for each *sample* in treeDir.
    combines condor() and hadd_simu() into one step to (hopefully) save time
    """
    import os
    import sets
    allfiles = os.listdir(treeDir)
    try:
        os.mkdir('out')
        os.mkdir('err')
        os.mkdir('log')
    except OSError: 
        pass
        
    files = os.listdir(treeDir)
    prefixes = [os.path.basename(f)[:7] for f in files if f.endswith('.root')]
    uniq = sets.Set(prefixes)
    
    ## need to stringify modlist
    _modlist = [mod.__name__ for mod in modlist]
    
    ## build the script that we will run -- note trick with sys.argv
    f = open('job.py', 'w')
    f.write('import sys\n')
    f.write('import analysis\n')
    f.write('analysis.histos2.writeHistograms(\'%s\', globber=\'*%%s*\' %% \
        sys.argv[1], trigList=%s, modlist=%s, histdir=\'%s\')\n' % (treeDir, 
        trigList, _modlist, histdir))
    
    ## build the submit.condor file
    f = open('submit.condor', 'w')
    f.write('executable = /usr/bin/python\n')
    f.write('getenv = True\n')
    f.write('notification = Error\n')
    f.write('notify_user = kocolosk@rcf.rhic.bnl.gov\n')
    f.write('universe = vanilla\n')
    f.write('stream_output = True\n')
    f.write('stream_error = True\n')
    f.write('transfer_executable = False\n\n')
    
    ## add jobs for each simulation sample
    for sample in uniq:
        f.write('output = out/%s.out\n' % sample)
        f.write('error = err/%s.err\n' % sample)
        f.write('log = log/%s.condor.log\n' % sample)
        f.write('arguments = %s/job.py %s\n' % (os.getcwd(), sample))
        f.write('queue\n\n')
    
    f.close()
    
    ## and off we go
    os.system('condor_submit submit.condor')

