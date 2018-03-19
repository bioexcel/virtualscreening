#!/uisr/bin/python3
#
#package CMIP::Run;
#
#=head1 NAME
#
#CMIP::Run - Object to run CMIP jobs
#
#=head1 SYNOPSIS
#
#use CMIP::Run;
#use CMIP::InputParameters;
#
#$cmipRun = CMIP::Run->new($inputParameters);
#
#=head1 METHODS
#
#$value=$cmipRun->type('direct'|'batch');
#
#$cmipRun->addFile(id1=>file1, id2=>file2,...);
#
#$cmipRun->delFile(id); 
#
#$hashref = $cmipRun->files;
#
#$resultref = $cmipRun->execute($queue);
#
#$cmipRun->rmTmpDir;
#
#=cut

import os
import shutil
import tempfile
import sys
import subprocess

class Run():
    def __init__(self,inputP, type='direct'):
        workdir=os.getcwd()
        self.workdir = workdir
        self.tmpdir = tempfile.mkdtemp(prefix="CMIP")
        self.param = inputP
        self.type = type
        #self.files = {'vdw': Local.VDWPRM}
        self.files = {}
        self.queue = ''
        self.exefile = ''

    def addFile(self, FILES):
        for k in FILES.keys():
            self.files[k] = FILES[k]
        return self

    def delFile(self,file):
        self.files.pop(file, None)
        return self

    def setup(self):
        if 'i' not in self.files:
            try:
                TMPPAR = open (self.tmpdir+"/param", "w")
            except OSError:
                print ("Error: open "+ self.tmpdir + "/param")
                sys.exit(1)
            
            TMPPAR.write(str(self.param))
            TMPPAR.close()
            if self.exefile:
                cmd = self.exefile + " -i " + self.tmpdir +"/param"
            else:
                cmd = CMIP.Local.CMIP + " -i " + self.tmpdir + "/param"
        if 'o' not in self.files:
            self.addFile({'o' : self.tmpdir + "/cmip.out"})
        if 'l' not in self.files:
            self.addFile({'l' : self.tmpdir + "/cmip.log"})
        if  self.type ==  'direct':
            self.addFile(
                {
                    "stderr": self.tmpdir + "/stderr.log",
                    "stdout": self.tmpdir + "/stdout.log"
                }
            )
        for k in self.files.keys():
            cmd = cmd + " -{} {}".format(k,self.files[k])

        runscript = self.tmpdir + "/run.csh"
        try:
            CSH = open (runscript,"w")
        except OSError as e:
            print ("Error: {} {} ({})".format(e.errno, e.strerror,runscript))
            sys.exit(1)
        #print (("#!/bin/tcsh -f\nhostname\ncd "+ self.tmpdir + "\n" + cmd + "\n"))
        CSH.write ("#!/bin/tcsh -f\nhostname\ncd "+ self.tmpdir + "\n" + cmd + "\n")
        CSH.close()
        return "/bin/tcsh ",runscript 
    
    def execute(self):
        cmdline = self.setup()        
        cproc = subprocess.getoutput (cmdline)            
        return CMIP.Result(self.files)
    		
    def rmTmpDir(self):
        shutil.rmtree(self.tmpdir)
