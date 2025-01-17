#!/usr/bin/env python

"""Python wrapper for the GROMACS cluster module
"""
import os
import sys
import json
import ntpath
import numpy as np
import configuration.settings as settings
from command_wrapper import cmd_wrapper
from tools import file_utils as fu

class Cluster(object):
    """Wrapper for the 5.1.2 version of the rms module
    Args:
        input_gro_path (str): Path to the original (before launching the trajectory) GROMACS structure file GRO.
        input_xtc_path (str): Path to the GROMACS trajectory.
        output_pdb_path (str): Path to the output cluster PDB file.
        properties (dic):
            gmx_path (str): Path to the GROMACS executable binary.
            dista (bool): Use RMSD of distances instead of RMS deviation
            method (str): Method for cluster determination: linkage, jarvis-patrick, monte-carlo, diagonalization, gromos
            cutoff (float): RMSD cut-off (nm) for two structures to be neighbor
    """

    def __init__(self, input_gro_path, input_xtc_path, output_pdb_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        self.input_gro_path = input_gro_path
        self.input_xtc_path = input_xtc_path
        self.output_pdb_path = output_pdb_path
        self.gmx_path = properties.get('gmx_path',None)
        self.mutation = properties.get('mutation',None)
        self.step = properties.get('step',None)
        self.path = properties.get('path','')
        self.mpirun = properties.get('mpirun', False)
        self.mpirun_np = properties.get('mpirun_np', None)
        self.cutoff = properties.get('cutoff', 0.1)


    def launch(self):
        """Launches the execution of the GROMACS rms module.
        """
        out_log, err_log = fu.get_logs(path=self.path, mutation=self.mutation, step=self.step)
        gmx = 'gmx' if self.gmx_path is 'None' else self.gmx_path

        cmd = [gmx, 'cluster',
               '-s', self.input_gro_path,
               '-f', self.input_xtc_path,
               '-cl', self.output_pdb_path,
               '-cutoff', str(self.cutoff)]

        if self.mpirun_np is not None:
            cmd.insert(0, str(self.mpirun_np))
            cmd.insert(0, '-np')
        if self.mpirun:
            cmd.insert(0, 'mpirun')
            cmd.append('<<<')
            cmd.append('\"'+"$'0\n0\n'"+'\"')
        else:
            cmd.insert(0, '|')
            cmd.insert(0, '\"'+'1 1'+'\"')
            cmd.insert(0, 'echo')
        command = cmd_wrapper.CmdWrapper(cmd, out_log, err_log)
        return command.launch()

#Creating a main function to be compatible with CWL
def main():
    system=sys.argv[1]
    step=sys.argv[2]
    properties_file=sys.argv[3]
    prop = settings.YamlReader(properties_file, system).get_prop_dic()[step]
    Rms(input_gro_path=sys.argv[4],
        input_trr_path=sys.argv[5],
        output_xvg_path=sys.argv[6],
        properties=prop).launch()

if __name__ == '__main__':
    main()
