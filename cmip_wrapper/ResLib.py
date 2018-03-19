#
# Residue library management
# Attypes and partial charges
# data from CMIP reslib formatted file
#
import sys
class ResiduesDataLib():
    def __init__(self,fname):
        self.RData={}
        try:
            fh = open(fname,"r")
        except OSError:
            print ("#ERROR while loading library file (" ,fname,")" )
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            r = Residue(data)
            self.RData[r.id]=r
        self.nres = len(self.RData)

    def getParams (self,resid,atid):
        id = resid+':'+atid
        if id in self.RData:
            return self.RData[id]
        else:
            #print ('WARNING: atom not found in library ('+id+')')
            return Residue([resid,atid,'X',0.])

class Residue():
    def __init__(self,data):
        self.id     = data[0]+':'+data[1]
        self.atType = data[2]
        self.charg  = float(data[3])

