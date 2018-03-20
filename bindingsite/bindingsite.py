#!/usr/bin/env python

"""Python wrapper module for the GROMACS solvate module
"""
import sys
import json
import configuration.settings as settings
from command_wrapper import cmd_wrapper
from tools import file_utils as fu

import os
import warnings
import Bio.PDB
import Bio.pairwise2
import Bio.SubsMat.MatrixInfo
from Bio.Data.SCOPData import protein_letters_3to1 as prot_one_letter


class BindingSite(object):
    """Finds the binding site of a PDB file based on the ligands' location of the PDB cluster90 members
    Args:
        structure_path (str): Path to the PDB protein structure.
        clusterPDBs_zip_path (str): Path to the TAR file containing the PDBs of cluster90 members.
        bindingsite_path (str): Path to the PDB containig the residues belonging to the binding site
        properties (dic):
            clusterPDBs_path (str): Path the uncompressed file containing the PDBs of cluster90 members
	    matrix (str): Substitution matrices for use in alignments. Available matrices are: 'benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65', 'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180', 'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure'
	    gap_open (float): Gap open penalty
	    gap_extend (float): Gap extend penalty
	    trim_ends (boolean): Cut unaligned sequence ends
    """

    def __init__(self, structure_path, clusterPDBs_zip_path, bindingsite_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
	#inputs
        self.structure_path       = structure_path
        self.clusterPDBs_zip_path = clusterPDBs_zip_path
	#outputs
        self.bindingsite_path     = bindingsite_path
	# path tmp files
        self.clusterPDBs_path     = properties.get('clusterPDBs_path','tmp/cluster')
	# options pairwise align
	self.matrix               = properties.get('matrix', Bio.SubsMat.MatrixInfo.blosum62)
	self.gap_open             = properties.get('gap_open', -10.0)
	self.gap_extend           = properties.get('gap_extend', -0.5)
	self.trim_ends            = properties.get('trim_ends', True)



    def launch(self):
        """Launches the pipeline to find the binding site of the given structure
        """
	#out_log, err_log = fu.get_logs(path=self.path, mutation=self.mutation, step=self.step)

	##	
	## Loading and parsing reference PDB structure
	print "\nParsing input PDB structure %s..." % self.structure_path

	structPDB      = {}

	if os.path.isfile(self.structure_path):
	    parser         = Bio.PDB.PDBParser()
	    structure_name = os.path.basename(self.structure_path.split('.')[0])
	    structPDB      = parser.get_structure(structure_name,self.structure_path)[0]
	else:
	    raise Exception('Input file "structure_path" {1} is not a file'.format(self.structure_path))

	# Use only the fist chain
	n_chains = structPDB.get_list()
	if len(n_chains) != 1:
	    warnings.warn('More than one chain found in "structure_path". Using only first chain to find the binding site')

	for struct_chain in structPDB.get_chains():
	    structPDB = struct_chain
	
	# Get AA sequence
	structPDB_seq   = self.__get_pdb_sequence(structPDB)


	##
	## Selecting CA atoms from first model, first chain of reference PDB structure

	#struct_atoms = []
  	#for struct_res in structPDB:
	#    try:
	#        struct_atoms.append(struct_res['CA'])
	#    except KeyError:
	#        #do not append hets
	#	pass

	#if len(struct_atoms)==0:	
	#    raise Exception('Cannot find CA atoms (first chain, first model) for "structure_path" {1}. Cannot find binding site in a empty selection.'.format(self.structure_path))
	#else:
	#   print "    Found %s protein residues" % len(struct_atoms)


	##
        ## Untar clusterPDBs

	print "\nExtracting cluster PDB structures from %s" % self.clusterPDBs_zip_path
        self.untar_file(tar_file=self.clusterPDBs_zip_path, dest_dir=self.clusterPDBs_path)

	##
	## Loop each clusterPDB

	for (dirpath, dirname, filenames) in os.walk(self.clusterPDBs_path):
	    if not filenames:
		next
	    for filename in filenames:

		##
		## Loading and parsing cluster PDB structure
		print "\nReading sequence cluster member  %s..." % filename

		clusterPDB = {}
		if not filename.endswith(('pdb', 'ent')):
		    raise Exception('clusterPDB TAR file contains file {0}. Format not supported. Must be .pdb/.ent'.format(filename))

		clusterPDB_name = os.path.basename(filename.split('.')[0])
		clusterPDB      = parser.get_structure(clusterPDB_name, dirpath+"/"+filename)[0]

		# Use only the fist chain
		for cluster_chain in clusterPDB.get_chains():
		    clusterPDB = cluster_chain
	
		# Get AA sequence
		clusterPDB_seq = self.__get_pdb_sequence(clusterPDB)

		##
		## Mapping residues by sequence alignment to match structPDB-clusterPDB paired residues

    		aln, residue_map = self.__align_sequences(structPDB_seq,clusterPDB_seq)
		print "    Aligned sequence is:\n    %s" % aln[1];

	        # Calculate (gapless) sequence identity
	        seq_identity, gap_seq_identity = self.__calculate_alignment_identity(aln[0], aln[1])
		print "    Sequence identity (%%): %s" % seq_identity
		print "    Gap less identity (%%): %s" % gap_seq_identity


		##
		## Selecting aligned CA atoms from first model, first chain

		print residue_map

		struct_atoms  = []
		cluster_atoms = []
    		for struct_res in residue_map:
		    struct_atoms.append(structPDB[struct_res]['CA'])
		    if (filename == "4WRG_A.pdb"):
			    print "clusterPDB[%s] = " % residue_map[struct_res]
			    print clusterPDB[residue_map[struct_res]]

		    cluster_atoms.append(clusterPDB[residue_map[struct_res]]['CA'])

		if len(cluster_atoms)==0:
		    raise Exception('Cannot find CA atoms (1st model, 1st chain) matching "structure_path" {1} with the cluster member {2}. Ignoring this member.'.format(self.structure_path,filename))
		else:
		    print "    Found %s aligned protein residues" % len(cluster_atoms)

		#cluster_atoms = []
		#for cluster_res in clusterPDB:
		#    try:
		#        cluster_atoms.append(cluster_res['CA'])
		#    except KeyError:
		#        #do not append hets
		#    	pass


		##
		## Align against input structure

		si = Bio.PDB.Superimposer()
		si.set_atoms(struct_atoms, cluster_atoms)
		si.apply(clusterPDB.get_atoms())


		# Print RMSD:
		print "    Cluster PDB superimposed. RMSD=%s" %si.rms

		# Save the aligned version of 1UBQ.pdb
		io = Bio.PDB.PDBIO()
		io.set_structure(clusterPDB)
		clusterPDB_aligned_path = dirpath+"/"+clusterPDB_name+"_aligned.pdb"
		print "    Writting transformed cluster member into %s" % clusterPDB_aligned_path
		io.save(clusterPDB_aligned_path)

	
    def __get_pdb_sequence(seft, structure):
        """
        Retrieves the AA sequence from a PDB structure.
        """

        aa = lambda r: (r.id[1], prot_one_letter.get(r.resname, 'X'))
	seq=[]
        for r in structure.get_residues():
	    if Bio.PDB.Polypeptide.is_aa(r):
		seq.append(aa(r))
        return seq

    def __align_sequences(self, seqA, seqB):
	"""
	Performs a global pairwise alignment between two sequences using the Needleman-Wunsch algorithm as implemented in Biopython.
	Returns the alignment and the residue mapping between both original sequences.
	"""

	#seq list to seq string
	sequence_A = ''.join([i[1] for i in seqA])
	sequence_B = ''.join([i[1] for i in seqB])

	# Do pairwaise alignment
        alns = Bio.pairwise2.align.globalds(sequence_A, sequence_B,self.matrix, self.gap_open, self.gap_extend, penalize_end_gaps=(False, False) )


	best_aln = alns[0]
	aligned_A, aligned_B, score, begin, end = best_aln

	# Equivalent residue numbering. Relative to reference
	mapping = {}
	aa_i_A, aa_i_B = 0, 0
	for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
	    if aa_aln_A == '-':
	        if aa_aln_B != '-':
		    aa_i_B += 1
	    elif aa_aln_B == '-':
	        if aa_aln_A != '-':
		    aa_i_A += 1
	    else:
	        assert seqA[aa_i_A][1] == aa_aln_A
	        assert seqB[aa_i_B][1] == aa_aln_B
	        mapping[seqA[aa_i_A][0]] = seqB[aa_i_B][0]
	        aa_i_A += 1
	        aa_i_B += 1

	return ((aligned_A, aligned_B), mapping)


    def __calculate_alignment_identity(self,alignedA,alignedB):
        """
        Returns the percentage of identical characters between two sequences
        """
        matches = [alignedA[i] == alignedB[i] for i in xrange(len(alignedA))]
        seq_id = (100 * sum(matches)) / len(alignedA)

        gapless_sl = sum([1 for i in xrange(len(alignedA)) if (alignedA[i] != '-' and alignedB[i] != '-')])
        gap_id = (100 * sum(matches)) / gapless_sl
        return (seq_id, gap_id)




    def untar_file(selft, tar_file=None, dest_dir=None):
	"""Untar (and uncompress - gzip,bz2) files into destionation folder
	"""
	import tarfile
	import errno

	# open TAR file
	tar_path = os.path.abspath(tar_file)
	tar = tarfile.open(tar_path)

	# create destination folder
	try:
	    os.makedirs(dest_dir)
	except OSError as exc:
	    if exc.errno == 17 and os.path.isdir(dest_dir):
                pass
            else:
                raise

	# extract all TAR files
	tar.extractall(path=dest_dir)
	tar.close() 


#Creating a main function to be compatible with CWL
def main():
    system         = sys.argv[1]
    step           = sys.argv[2]
    properties_file= sys.argv[3]
    prop           = settings.YamlReader(properties_file, system).get_prop_dic()[step]

    BindingSite(structure_path       = sys.argv[4],
                clusterPDBs_zip_path = sys.argv[5],
                bindingsite_path     = sys.argv[6],
                properties           = prop        ).launch()


if __name__ == '__main__':
    main()
