#!/usr/bin/env python

"""Python wrapper module for the AutoDockVina sdf2pdb module"""
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
        self.input_sdf_path = input_sdf_path
        self.output_pdb_path = output_pdb_path

    def launch(self):
        """Launches the execution of the sdf2pdb module."""
        out_log, err_log = fu.get_logs(path=self.path, mutation=self.mutation, step=self.step)
        pdb_list = []
        for mol in pybel.readfile("sdf", self.input_sdf_path):
            fname=str(mol.title)+'.pdb'
            mol.write("pdb", fname)
            pdblist.append(os.path.abspath(fname))

        return pdb_list


def main():
    """Main function to be compatible with CWL"""
    system=sys.argv[1]
    step=sys.argv[2]
    properties_file=sys.argv[3]
    prop = settings.YamlReader(properties_file, system).get_prop_dic()[step]
    SDF2PDB(input_sdf_path=sys.argv[4],
            output_pdb_path=sys.argv[5],
            properties=prop.launch()


if __name__ == '__main__':
    main()
