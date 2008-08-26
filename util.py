def getRun(path):
    """searches for an integer runnumber in the supplied path.
    if not found, return the filename after stripping path and suffix info"""
    import re
    import os
    regex = re.search('(6|7)\d{6}',path)
    try:
        return int(regex.group())
    except AttributeError:
        return os.path.basename(path).split('.')[0]


def getFill(runnumber):
    """queries DB to get fill for this runnumber"""
    import MySQLdb
    db = MySQLdb.connect(host='star1.lns.mit.edu', port=3316, db='RunLog_onl')
    dbc = db.cursor()
    nrows = dbc.execute('select blueFillNumber from beamInfo where runNumber=%d and \
        deactive=0 order by beginTime desc limit 1' % (runnumber,))
    if nrows ==0:
        print runnumber, 'has no beamInfo record in the RunLog_onl database'
        return 0
    else:
        fill = int(dbc.fetchone()[0])
        ## temporary hacks
        ## http://www.star.bnl.gov/HyperNews-star/protected/get/starspin/3324.html
        if 6144002 <= runnumber <= 6144029: fill = 7128
        if 6144041 <= runnumber <= 6144042: fill = 7129
        if 6145067 <= runnumber <= 6145068: fill = 7136
        if 6146001 <= runnumber <= 6146026: fill = 7138
        return fill

def getAllFills(runlist):
    """returns tuple of (run,fill) tuples"""
    import MySQLdb
    db = MySQLdb.connect(host='star1.lns.mit.edu', port=3316, db='RunLog_onl')
    dbc = db.cursor()
    nrows = dbc.execute('select runNumber,blueFillNumber from beamInfo where runNumber>=%d \
        and runNumber<=%d and deactive=0 order by runNumber desc' % (min(runlist),max(runlist)))
    return dbc.fetchall()


def getAllFills2(runlist):
    """this one returns a run,fill dict containing only entries from specified runs"""
    tuples = getAllFills(runlist)
    d = {}
    for run,fill in tuples:
        if run in runlist:
            d[run] = int(fill)
    return d


def uniqify(seq, idfun=None):
    # order-preserving method to return unique elements in list
    if idfun is None: 
        def idfun(x): return x 
    seen = {} 
    result = [] 
    for item in seq: 
        marker = idfun(item) 
        # in old Python versions: 
        # if seen.has_key(marker) 
        # but in new ones: 
        if marker in seen: continue 
        seen[marker] = 1 
        result.append(item) 
    return result

def daqEventCount(run, trigId):
    """queries RunLog Browser to get event count for this (run,trigId) pair"""
    """basically the world's stupidest HTML parser"""
    import urllib2
    sock = urllib2.urlopen('http://online.star.bnl.gov/RunLogRun6/Summary.php?run=%d' % run)
    for line in sock:
        items = line.split()
        for row,item in enumerate(items):
            if item.find(str(trigId)) > 0:
                eventCountItem = items[row+4]
                return int(eventCountItem.lstrip('align="right">').rstrip('</TD>R'))
    
    

def lafferty_wyatt_point(lowedge, highedge, expo_slope):
    """calculates the l-w point for a bin where the true distribution is an
    exponential characterized by expo_slope.
    """
    import math
    rhs = (math.exp(expo_slope*highedge) - math.exp(expo_slope*lowedge)) / expo_slope
    rhs /= (highedge - lowedge)
    return math.log(rhs) / expo_slope
    
def hadd_interactive(histDir, runlist, trig, spin, charge, key):
    import ROOT
    from glob import glob
    import os.path
    allFiles = glob( os.path.join(histDir, "*.root") )
    if charge == None:
        keystring = "_%s_%s_%s" % (trig,spin,key)
    else:
        keystring = "_%s_%s_%s_%s" % (trig,spin,charge,key)
    ## should check that the first file is in the runlist
    f = ROOT.TFile(allFiles[0])
    h = f.Get(keystring).Clone()
    h.SetDirectory(0)
    f.Close()
    for fname in allFiles[1:]:
        run = getRun(fname)
        if runlist is None or run in runlist:
            print fname, keystring
            tmp = ROOT.TFile(fname)
            h.Add( tmp.Get(keystring) )
            tmp.Close()
    return h
