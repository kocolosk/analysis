class taimask:
    def __init__(self):
        f = open('/Users/kocolosk/data/run6/fillPatternBxRun6_070412.txt')
        self.d = {}
        for line in f.readlines():
            l2 = line.split()
            l3 = [int(elem) for elem in l2]
            self.d[l3[0]] = l3[1:]
            
    def useBxing(self, fill, bx):
        try:
            record = self.d[fill]
            return record[bx]    
        except KeyError:
            print 'no record found for fill %d' % (fill,)
            return 0
