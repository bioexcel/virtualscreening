#! /usr/bin/python
"""Python wrapper module for CMIP
"""
__author__ = "gelpi"
__date__ = "$24-nov-2017 12:30:59$"

import re
import sys
from Grid import Grid
from InputParams import InputParams
from Run import Run
from Result import Result
import tools.file_utils as fu
from command_wrapper import cmd_wrapper



class CMIPWrapper():
    """Wrapper class for the 2.7 version of CMIP
        Adapted from CMIP standard python wrapper
    """

    def __init__(self, paths, props):
        self.cmipPaths=[]
        self.paths={}
        
        for k in paths.keys():
            if k == 'cmip':
                self.cmipPaths=paths[k]
            else:
                self.paths[k]=paths[k]
        self.cmipPaths.append({'hs':paths['input_pdb_path']})
        
        self.properties={}
        self.properties['step'] =  props['step']
        for k in props.keys():
            if k == 'cmip':
                self.properties['cmipkwds']=props[k]
            else:
                self.properties[k]=props[k]        
        
        
        gr = self._prepGrid()
        if 'pbfocus' in self.properties['cmipkwds'] and\
            self.properties['cmipkwds']['pbfocus'] == 1:
            gr0 = _prepGrid(1)
            inpP = InputParams(self.properties['step'], gr, gr0)
        else:
            inpP = InputParams(self.properties['step'], gr)

        for prm in self.properties['cmipkwds'].keys():
            inpP.addKeyword ({prm.upper():self.properties['cmipkwds'][prm]})
        
        self.run = Run(inpP,'wrapper')        
        
        if 'titration' in self.properties['cmipkwds']:
            self.run.exefile = self.properties['CMIP_root']+"/src/titration"
        else:
            self.run.exefile = self.properties['CMIP_root']+"/src/cmip"
        
        for fn in self.cmipPaths:
            self.run.addFile(fn)
        
        if 'vdw' not in self.run.files:
            self.run.addFile({'vdw':self.properties['CMIP_root']+'/dat/vdwprm'})    
    
    def launch(self,deleteTmpDir=True):
    # using cmd_wrapper
        self.out_log, self.err_log = fu.get_logs(path=self.properties['path'],step=self.properties['step'])
        command = cmd_wrapper.CmdWrapper(self.run.setup(), self.out_log, self.err_log)
        command.launch()
    # using original CMIP execute
        # self.result = self.run.execute
        self.result = Result(self.run.tmpdir,self.run.files)
        if deleteTmpDir:
            self.run.rmTmpDir()
        return self.result
        
    def _prepGrid(self,outgrid=0):        
        if outgrid:
            lb='0'
        else:
            lb=''
        gr = Grid(outgrid)
        if 'readgrid'+lb in self.properties['cmipkwds']:
            gr.setreadGrid(self.properties['cmipkwds']['readgrid'+lb])
            del self.properties['cmipkwds']['readgrid'+lb]
        if 'perfill'+lb in self.properties['cmipkwds']:
            gr.setperfill (self.properties['cmipkwds']['perfill'+lb])
            del self.properties['cmipkwds']['perfill'+lb]
        if 'cen'+lb in self.properties['cmipkwds']:
            gr.cen = self.properties['cmipkwds']['cen'+lb]
            del self.properties['cmipkwds']['cen'+lb]
        if 'dim'+lb in self.properties['cmipkwds']:
            gr.dim = self.properties['cmipkwds']['dim'+lb]
            del self.properties['cmipkwds']['dim'+lb]
        if 'int'+lb in self.properties['cmipkwds']:
            gr.int = self.properties['cmipkwds']['int'+lb]
            del self.properties['cmipkwds']['int'+lb]
        gr.update()
        return gr
        
        

def main():
    step=sys.argv[1]
    propfn=sys.argv[2]
    props = settings.YamlReader(propfn, 'linux').get_prop_dic()[step]
    props['step']=step
    paths = settings.YamlReader(propfn, 'linux').get_paths_dic()[step]

    cw = CMIPWrapper(paths,props)
    cw.launch(True)

if __name__ == "__main__":
    main()
