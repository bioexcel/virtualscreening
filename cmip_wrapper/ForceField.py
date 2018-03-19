#
# Manage Forcefield parameters (vdw, Srf)
# uses modified CMIP vdwprm file
#
import sys

class VdwParamset():
    def __init__ (self, fname):
        self.atTypes = {}
        try:
            fh = open(fname,"r")
        except OSError:
            print ("#ERROR while loading parameter file (" ,fname,")" )
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            self.atTypes[data[0]]=AtType(data)
        self.ntypes = len(self.atTypes)
        fh.close()

class AtType():
    def __init__(self,data):
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612



