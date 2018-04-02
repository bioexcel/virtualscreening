#!/usr/bin/env python

"""Python wrapper module for AutoDock Vina
"""
import sys
import re
import os
import json
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import configuration.settings as settings
from command_wrapper import cmd_wrapper
from tools import file_utils as fu


class AutoDockVina(object):
    """Wrapper class for AutoDock Vina software"""
    def __init__(self, ligand_pdbqt_path, receptor_pdbqt_path,
                 box_path, log_file, output_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        self.ligand_pdbqt_path = ligand_pdbqt_path
        self.receptor_pdbqt_path = receptor_pdbqt_path
        self.box_path = box_path
        self.log_file = log_file
        self.output_path = output_path
        self.vina_path = properties.get('vina_path',None)
        self.path = properties.get('path','')
        self.step = properties.get('step','')

    @staticmethod
    def _calculate_box(box_file_path):
        with open(box_file_path, 'r') as box_file:
            for line in box_file:
                line = line.rstrip(os.linesep)
                if line.startswith("REMARK BOX CENTER"):
                    fields = line.split()
                    center = fields[3:6]
                    size = fields[-3:]
                    return center[0], center[1], center[2], size[0], size[1], size[2]
            return 0, 0 ,0, 0, 0, 0

    def launch(self):
        """Launches the execution of the AutoDock Vina docking"""
        out_log, err_log = fu.get_logs(path=self.path, step=self.step)

        vina = 'vina' if self.vina_path is None else self.vina_path
        x0, y0, z0, sidex, sidey, sidez = AutoDockVina._calculate_box(self.box_path)
        cmd = [vina,
               '--ligand', self.ligand_pdbqt_path,
               '--receptor', self.receptor_pdbqt_path,
               '--center_x=' + x0, '--center_y=' + y0, '--center_z=' + z0,
               '--size_x=' + sidex, '--size_y=' + sidey, '--size_z=' + sidez,
               '--out', out_log,
               '--log', log_file]

        command = cmd_wrapper.CmdWrapper(cmd, out_log, err_log)
        return command.launch()

#Creating a main function to be compatible with CWL
def main():
    system=sys.argv[1]
    step=sys.argv[2]
    properties_file=sys.argv[3]
    prop = settings.YamlReader(properties_file, system).get_prop_dic()[step]
    AutoDockVina(ligand_pdbqt_path=sys.argv[4], receptor_pdbqt_path=sys.argv[5],
                 box_path=sys.argv[5], log_file=sys.argv[6], properties=prop).launch()

if __name__ == '__main__':
    main()
