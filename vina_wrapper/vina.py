#!/usr/bin/env python

"""Python wrapper module for AutoDock Vina
"""
import sys
import re
import json
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import configuration.settings as settings
from command_wrapper import cmd_wrapper
from tools import file_utils as fu


class AutoDockVina(object):
    """Wrapper class for AutoDock Vina software"""
    def __init__(self, ligand_pdb_path, ligand_name, receptor_pdb_path,
                 box_path, output_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        self.ligand_pdb_path = ligand_pdb_path
        self.ligand_name = ligand_name
        self.receptor_pdb_path = receptor_pdb_path
        self.box_path = box_path
        self.output_path = output_path
        self.vina_path = properties.get('vina_path',None)
        self.path = properties.get('path','')
        self.step = properties.get('step','')


    def launch(self):
        """Launches the execution of the AutoDock Vina docking"""
        out_log, err_log = fu.get_logs(path=self.path, step=self.step)

        vina = 'vina' if self.vina_path is None else self.vina_path
        x0, y0, z0 = (0, 0, 0)
        sidex, sidey, sidez = (5, 5, 5)
        cmd = [vina, '--ligand', self.ligand_pdb_path,  '--receptor', self.receptor_pdb_path,
               '--center_x=' + x0, '--center_y=' + y0, '--center_z=' + z0,
               '--size_x=' + sidex, '--size_y=' + sidey, '--size_z=' + sidez,
               "--out ", out_log + " --log " + logfile', prepared_file_path, '-o', self.output_pdb_path, '-h', '-t']
        if self.mutation:
            cmd.append('-s')
            cmd.append(sequence_file_path)

        command = cmd_wrapper.CmdWrapper(cmd, out_log, err_log)
        return command.launch()

#Creating a main function to be compatible with CWL
def main():
    system=sys.argv[1]
    step=sys.argv[2]
    properties_file=sys.argv[3]
    prop = settings.YamlReader(properties_file, system).get_prop_dic()[step]
    AutoDockVina(ligand_pdb_path=sys.argv[4], ligand_name=sys.argv[5],
                 receptor_pdb_path=sys.argv[6], box_path=sys.argv[7],
                 output_path=sys.argv[8], properties=prop).launch()

if __name__ == '__main__':
    main()
