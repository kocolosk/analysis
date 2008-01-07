import sys
import getopt
import re
import math
print 'pyroot -- import ROOT'
import ROOT

paths = {  'jets':'/Volumes/scratch/common/run6/jets/jets_',
                    'jetSkim':'/Volumes/scratch/common/run6/jetSkim/jetSkim_',
                    'chargedPions':'/Volumes/scratch/common/run6/chargedPions/chargedPions_',
                    'bemcPions':'/Volumes/scratch/common/run6/neutralPions/NeutralPionTree.root'
}

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg
        
def removeFriend(inputPath):
    print inputPath
    spinFile = ROOT.TFile.Open(inputPath,'update')
    if spinFile is None:
        raise Usage('could not open spin file at %s' % (inputPath,))
    spinTree = spinFile.Get('spinTree')
    if spinTree is None:
        raise Usage('could not get spinTree in file at %s' % (inputPath,))
    
    #friends = [spinFile.Get('ConeJets12'),spinFile.Get('ConeJetsEMC'),spinFile.Get('chargedPions'),spinFile.Get('bemcPions')]
    friends = [spinFile.Get('ConeJets'),spinFile.Get('chargedPions')]
    [elem.SetTitle('this can be a friend') for elem in friends]
    
    #remove the previous jet friend tree and create a new one
    listOfFriends = spinTree.GetListOfFriends()
    if listOfFriends is not None:
        for i in xrange(listOfFriends.GetEntries()):
            f = listOfFriends.At(0)
            print f.GetName()
            spinTree.RemoveFriend(f.GetTree())
    
    spinFile.cd()
    spinTree.Write('',ROOT.TObject.kOverwrite)
    [elem.Write('',ROOT.TObject.kOverwrite) for elem in friends]
    spinFile.Close()

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hcu:", ['help',"create", "update=", 'site='])
        except getopt.error, msg:
            raise Usage(msg)

        #magic!
        throwitaway = ROOT.std.map(long,float)()

        ROOT.gSystem.Load('StarSpinAnalyses')
        [removeFriend(elem) for elem in args]
            

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
