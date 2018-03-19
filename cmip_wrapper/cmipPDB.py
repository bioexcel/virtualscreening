#!/usr/bin/env python

"""Python wrapper module for SCWRL
"""
import sys
import re
import json
from os.path import join as opj
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import configuration.settings as settings
from command_wrapper import cmd_wrapper
from tools import file_utils as fu

class CMIPPDB(object):
    """Class to set up PDB files to be used by CMIP with COORFMT=2.
    Args:
        input_pdb_path (str): Path to the input PDB file.
        output_pdb_path (srt): Path to the output prepared PDB file.
        properties (dic): All properties and system path
    """
    def __init__(self, input_pdb_path, output_pdb_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        self.input_pdb_path = input_pdb_path
        self.output_pdb_path = output_pdb_path
        self.CMIP_root_path = properties.get('CMIP_root_path', None)
        self.path = properties.get('path','')
        self.step = properties.get('step','')
        self.mutation = properties.get('mutation', None)

    def launch(self):
        """Launches the execution of the SCWRL binary.
        """
        out_log, err_log = fu.get_logs(path=self.path, mutation=self.mutation, step=self.step)

        aaLib = ResiduesDataLib(opj(self.CMIP_root_path, 'dat', 'res.lib')
        out_log.info("Residue or atom pairs loaded from: "+)




        scrwl = 'Scwrl4' if self.scwrl4_path is None else self.scwrl4_path
        cmd = [scrwl, '-i', prepared_file_path, '-o', self.output_pdb_path]

        command = cmd_wrapper.CmdWrapper(cmd, out_log, err_log)
        return command.launch()

#Creating a main function to be compatible with CWL
def main():
    system=sys.argv[1]
    step=sys.argv[2]
    properties_file=sys.argv[3]
    prop = settings.YamlReader(properties_file, system).get_prop_dic()[step]
    Scwrl4(input_pdb_path=sys.argv[4], output_pdb_path=sys.argv[5], properties=prop).launch()

if __name__ == '__main__':
    main()
