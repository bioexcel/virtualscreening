#
#
import re
import sys

DEFAULT_FILES = {
    'i' : 'param',
    'vdw' : 'vdwprm',
    'hs' : 'hostpdb',
    'epr' : 'energies',
    'grdin' : 'gridin',
    'o' :'output',
    'grdout' : 'gridout',
    'outpdb' : 'outpdb',
    'grdb' : 'grdB',
    'rst' : 'restart',
    'pr' : 'probpdb',
    'grad' : 'gradient',
    'gds' :  'gdesolv',
    'srf' : 'surface',
    'eps' : 'outeps',
    'ats' : 'ats',
    'eps2' : 'outeps2',
    'fr' : 'outfrac',
    'byat' : 'byatom',
    'atsgrd' : 'atsgrd',
    'l' : 'log',
    'cube' : 'cube'
 }


class Result():
    
    def __init__(self, workdir='', files={}):
        self.stdout=''
        self.stderr=''
        self.workdir=workdir
        self.logData={}
        self.files=files
        self.errstr=''
        if 'o' in files:
            self.output = _getFile(files['o'])
            if self.output == 0:
                self.errstr = self.output
        if 'stdout' in files:
            self.stdout = _getFile(files['stdout'])
            if self.output == 0:
                self.errstr = self.output
        if 'stderr' in files:
            self.stderr = _getFile(files['stderr'])
            if self.output == 0:
                self.errstr = self.output
        if self.stderr == '' and files['l']:
            try:
                LOG  = open (files['l'],"r")
                runtitle=''
                strtitle=''
                for line in LOG:
                    line = line.rstrip()
                    if runtitle:
                        if re.search(runtitle,line):
                            self.logData[runtitle.replace(' ','')] = strtitle
                            runtitle=''
                            strtitle=''
                        else:
                            strtitle = strtitle + line
                    else:
                        if re.search("title",line):
                            [runtitle,v] = line.split('=')
                            strtitle=v
                        else:
                            line= line.replace(' ','')
                            [k,v] = line.split('=')
                            k = k.replace(' ','')
                            v = v.replace(' ','')
                            self.logData[k]=v
                LOG.close()
            except IOError:
                self.errstr = "CMIP Log file not found"
        else:
            self.errstr = subprocess.run ("tail -1 "+ files['o'])

    def getFileName(self,key):
        if key in self.files:
            return self.files[key]
        else:
            return self.workdir+'/'+DEFAULT_FILES[key]

    def getFileContents(self, key):
        fn = self.getFileName(key)
        return _getFile(fn)

    def asString(self):
        txt=""
        for k in sorted(self.logData.keys()):
            txt = txt + k +" = " + self.logData[k]+"\n"
        return txt

def _getFile(fn):
    try:
        with open(fn,"r") as content:
            file = content.read()
        return file
    except IOError:
        return "file "+ fn + " not found"
        
