#!/usr/bin/env python

"""Python wrapper module for the AutoDockVina sdf2pdb module"""
import sys
import json
import openbabel
import configuration.settings as settings
from tools import file_utils as fu


class SDF2PDB(object):
    """Encapsulates logic to convert an SDF file to PDB structure"""
    def __init__(self, input_sdf_path, output_pdb_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        self.input_sdf_path = input_sdf_path
        self.output_pdb_path = output_pdb_path

    def launch(self):
        """Launches the execution of the sdf2pdb module."""
        out_log, err_log = fu.get_logs(path=self.path)

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("sdf", "pdb")

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, self.input_sdf_path)
        obConversion.WriteFile(mol, self.output_pdb_path)


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
