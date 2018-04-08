#!/usr/bin/env python

"""Python wrapper module for the AutoDockVina sdf2pdb module"""
import os
import sys
import json
import pybel
import configuration.settings as settings
from tools import file_utils as fu


class SDF2PDB(object):
    """Encapsulates logic to convert an SDF file to PDB structure"""
    def __init__(self, input_sdf_path, properties, **kwargs):
        self.mutation = properties.get('mutation',None)
        self.step = properties.get('step',None)
        self.path = properties.get('path','')
        self.max_pdbs = int(properties.get('max_pdbs', 0))
        self.input_sdf_path = input_sdf_path

    def launch(self):
        """Launches the execution of the sdf2pdb module."""
        out_log, err_log = fu.get_logs(path=self.path, mutation=self.mutation, step=self.step)
        pdb_list = []
        for mol_index, mol in enumerate(pybel.readfile("sdf", self.input_sdf_path)):
            if mol_index >= self.max_pdbs and self.max_pdbs !=0 : break
            fname=str(mol.title)+'.pdb'
            mol.write("pdb", fname, overwrite=True)
            pdb_list.append(os.path.abspath(fname))

        return pdb_list


def main():
    """Main function to be compatible with CWL"""
    system=sys.argv[1]
    step=sys.argv[2]
    properties_file=sys.argv[3]
    prop = settings.YamlReader(properties_file, system).get_prop_dic()[step]
    SDF2PDB(input_sdf_path=sys.argv[4],
            output_pdb_path=sys.argv[5],
            properties=prop).launch()


if __name__ == '__main__':
    main()
