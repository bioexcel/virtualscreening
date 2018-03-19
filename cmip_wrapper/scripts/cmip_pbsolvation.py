"""Python wrapper module for SCWRL
"""
import sys
import re
import json
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import configuration.settings as settings
from command_wrapper import cmd_wrapper
from tools import file_utils as fu

class Scwrl4(object):
    """Wrapper class for the 4.0 version of SCWRL.
    Args:
        input_pdb_path (str): Path to the input PDB file.
        output_pdb_path (srt): Path to the output mutated PDB file.
        properties (dic): All properties and system path
    """
    def __init__(self, input_pdb_path, output_pdb_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        self.input_pdb_path = input_pdb_path
        self.output_pdb_path = output_pdb_path
        pattern = re.compile(("(?P<chain>[a-zA-Z*]{1}).(?P<wt>[a-zA-Z]{3})(?P<resnum>\d+)(?P<mt>[a-zA-Z]{3})"))
        self.mutation = pattern.match(properties['mutation']).groupdict()
        self.scwrl4_path = properties.get('scwrl4_path',None)
        self.path = properties.get('path','')
        self.mut = properties.get('mutation','')
        self.step = properties.get('step','')


    def launch(self):
        """Launches the execution of the SCWRL binary.
        """
        out_log, err_log = fu.get_logs(path=self.path, mutation=self.mut, step=self.step)
        if self.mutation is not None:
            # Read structure with Biopython
            parser = PDBParser(PERMISSIVE=1)
            st = parser.get_structure('s', self.input_pdb_path)  # s random id never used

            # Remove the side chain of the AA to be mutated
            if self.mutation['chain']!='*':
                chains = [self.mutation['chain']]
            else:
                chains = [chain.id for chain in st[0]]

            resnum = int(self.mutation['resnum'])

            for chain in chains:
                residue = st[0][chain][(' ', resnum, ' ')]
                backbone_atoms = ['N', 'CA', 'C', 'O', 'CB']
                not_backbone_atoms = []

                # The following formula does not work. Biopython bug?
                # for atom in residue:
                #     if atom.id not in backbone_atoms:
                #         residue.detach_child(atom.id)

                for atom in residue:
                    if atom.id not in backbone_atoms:
                        not_backbone_atoms.append(atom)
                for atom in not_backbone_atoms:
                    residue.detach_child(atom.id)

                # Change residue name
                residue.resname = self.mutation['mt'].upper()

            # Write resultant structure
            w = PDBIO()
            w.set_structure(st)
            prepared_file_path = self.mut+self.step+'prepared.pdb'
            w.save(prepared_file_path)
        else:
            prepared_file_path = self.input_pdb_path

        scrwl = 'Scwrl4' if self.scwrl4_path is None else self.scwrl4_path
        cmd = [scrwl, '-i', prepared_file_path, '-o', self.output_pdb_path]

        command = cmd_wrapper.CmdWrapper(cmd, out_log, err_log)
        return command.launch()

#Creating a main function to be compatible with CWL
def main():
    step=sys.argv[3]
    prop=sys.argv[4]
    step, system, mut = step.split(':')
    prop = settings.YamlReader(prop, system).get_prop_dic(mut)[step]
    prop['path']=''
    Scwrl4(input_pdb_path=sys.argv[1],
           output_pdb_path=sys.argv[2],
           step=step,
           properties=prop).launch()

if __name__ == '__main__':
    main()
