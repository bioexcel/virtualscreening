#!/usr/bin/env python

"""Python wrapper module for the AutoDockVina prepare_receptor module"""
import sys
import json
import os
from command_wrapper import cmd_wrapper
import configuration.settings as settings
from tools import file_utils as fu
from vina_errors import MGLToolsConfigurationError


class VinaPrepareReceptor(object):
    """Encapsulates logic to prepare a receptor PDB structure for docking with AutoDockVina"""
    def __init__(self, input_receptor_pdb_path, output_receptor_pdbqt_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        self.input_receptor_pdb_path = input_receptor_pdb_path
        self.output_receptor_pdbqt_path = output_receptor_pdbqt_path
        self.mgltools_path = properties.get('mgltools_path', None)
        if not self.mgltools_path:
            raise MGLToolsConfigurationError("Path not found")
        self.executor_path = os.path.join(self.mgltools_path, 'bin', 'pythonsh')
        self.script_path = os.path.join(self.mgltools_path, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_receptor4.py')
        self.path = properties.get('path','')

    def launch(self):
        """Launches the execution of the MGLTools prepare_receptor module."""
        out_log, err_log = fu.get_logs(path=self.path)

        cmd = [self.executor_path, self.script_path,
               '-r', self.input_receptor_pdb_path,
               '-o', self.output_receptor_pdbqt_path]

        command = cmd_wrapper.CmdWrapper(cmd, out_log, err_log)
        return command.launch()


def main():
    """Main function to be compatible with CWL"""
    system=sys.argv[1]
    step=sys.argv[2]
    properties_file=sys.argv[3]
    prop = settings.YamlReader(properties_file, system).get_prop_dic()[step]
    VinaPrepareReceptor(input_receptor_pdb_path=sys.argv[4],
                      output_receptor_pdbqt_path=sys.argv[5],
                      properties=prop).launch()


if __name__ == '__main__':
    main()
