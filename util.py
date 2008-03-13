def getRun(path):
    """searches for an integer runnumber in the supplied path"""
    import re
    regex = re.search('(6|7)\d{6}',path)
    return int(regex.group())


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
