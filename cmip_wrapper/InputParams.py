#!/usr/bin/python3

#CMIP.InputParams - Object to run prepare input parameters for CMIP jobs
# inputP = CMIP.InputParams(title, grid, [grid0])
#
# hashref=inputP.keywords
#
# string=inputP.asParamString
#
# inputP.addKeyword ({'key1':'value1', 'key2':'value2', ...})
#
# inputP.delKeyword (key)
# 
# inputP.grid(grid)
#
# inputP.grid0(grid);

import re
import sys

class InputParams():
    
    def __init__ (self, title='', gr='', gr0=''):
        self.title = title
        self.grid = gr
        self.keywords = {}
        if gr0:
            self.grid0 = gr0
            self.keywords['PBFOCUS'] = 1
        self.keywords['WRITELOG'] = 1
        
    def loadFromInput(self, fname):  ## provisional
        try:
            INPUT = open (fname, 'r')
        except OSError:
            print ("#ERROR while loading file (", fname, ")")
            sys.exit(2)
        KEYS = {}
        title = INPUT.readline()
        for line in INPUT:
            if re.search("&", line):
                continue
            for prm in  re.split(" *, *", line):
                [k, v] = re.split (" *= *", prm)
                k = k.replace(' ', '')
                KEYS[k.uppercase] = v
        INPUT.close()
        grid = CMIP.Grid()
        grid0 = CMIP.Grid()
        if 'READGRID' in KEYS:
            grid.setreadGrid(KEYS['READGRID'])
        grid.int(KEYS['INTX'], KEYS['INTY'], KEYS['INTZ'])
        grid.cen(KEYS['CENX'], KEYS['CENY'], KEYS['CENZ'])
        grid.dim(KEYS['DIMX'], KEYS['DIMY'], KEYS['DIMZ'])
        if 'PBFOCUS' in KEYS:
            if 'READGRID0' in KEYS:
                grid0.readGrid(KEYS['READGRID0'])
            grid0.int(KEYS['INTX0'], KEYS['INTY0'], KEYS['INTZ0'])
            grid0.cen(KEYS['CENX0'], KEYS['CENY0'], KEYS['CENZ0'])	
            grid0.dim(KEYS['DIMX0'], KEYS['DIMY0'], KEYS['DIMZ0'])
        inp = CMIP.InputParams(title, grid, grid0)
        for k in KEYS.keys():
            inp.addKeyword(k, KEYS[k])
        return inp
    
    def setgrid(self, gr):
        if gr:
            self.grid = gr
        return self['grid']


    def setgrid0(self, gr):
        if gr:
            self.grid0 = gr
            self.keywords['PBFOCUS'] = 1
        return self.grid0

    def __str__(self): 
        out = self.title + "\n&cntrl\n"
        for k in self.keywords.keys():
            out = out + " {}={}\n".format(k,str(self.keywords[k]))
        if 'PBFOCUS' in self.keywords:
            out = out + self.grid0.__str__()
        out = out + self.grid.__str__()
        out = out + "&end\n"
        return out

    def addKeyword (self, KEYW):
        for k in KEYW.keys():
            self.keywords[k] = KEYW[k]
        return self

    def delKeyword (self, key):
        self.keywords.pop(key, None)
        return self


