#! /usr/bin/python
"""Python wrapper module for preppdb.py
"""
__author__ = "gelpi"
__date__ = "$24-nov-2017 12:30:59$"

import re
import sys
import tools.file_utils as fu
from command_wrapper import cmd_wrapper


class prepPDBWrapper():
    """Wrapper class for preppdb script        
    """

    def __init__(self, paths, props):
        self.properties={}
        for k in props:
            self.properties[k]=props[k]
        if 'res_lib' not in paths:
            self.res_lib = props['CMIP_root']+"/dat/res.lib"
        else:
            self.res_lib = paths['res_lib']
        self.pdb_in = paths['input_pdb_path']
        self.pdb_out = paths['output_pdb_path']        
        self.script = props['CMIP_root']+"/python/preppdb.py"
    
    def launch(self):
    # using cmd_wrapper
        out_log, err_log = fu.get_logs(path=self.properties['path'],step=self.properties['step'])
        cmd =["python ",self.script, self.res_lib, self.pdb_in, self.pdb_out]
        command = cmd_wrapper.CmdWrapper(cmd, out_log, err_log)
        command.launch()

def main():
    step=sys.argv[1]
    propfn=sys.argv[2]
    props = settings.YamlReader(propfn, 'linux').get_prop_dic()[step]
    props['step']=step
    paths = settings.YamlReader(propfn, 'linux').get_paths_dic()[step]
    prepPDBWrapper(paths,props).launch()

if __name__ == "__main__":
    main()
