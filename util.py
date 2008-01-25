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
        return int(dbc.fetchone()[0])

def getAllFills(runlist):
    """returns tuple of (run,fill) tuples"""
    import MySQLdb
    db = MySQLdb.connect(host='star1.lns.mit.edu', port=3316, db='RunLog_onl')
    dbc = db.cursor()
    nrows = dbc.execute('select runNumber,blueFillNumber from beamInfo where runNumber>=%d \
        and runNumber<=%d and deactive=0 order by runNumber desc' % (min(runlist),max(runlist)))
    return dbc.fetchall()
