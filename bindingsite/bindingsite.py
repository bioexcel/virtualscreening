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
import re
from mmb_api import pdb
import Bio.PDB
import Bio.pairwise2
import Bio.SubsMat.MatrixInfo
from Bio.Data.SCOPData import protein_letters_3to1 as prot_one_letter

import numpy
#from Bio.Struct.Geometry import center_of_mass
#from Bio import Struct


class BindingSite(object):
    """Finds the binding site of the input_pdb file based on the ligands' location of similar structures (members of the sequence identity cluster)
    Args:
        pdb_code (str): PDB code for the protein structure where to binding site is to be find
        pdb_chain (str): Chain id for the pdb_code where to binding site is to be find
        output_pdb_path (str): Path to the PDB containig the residues belonging to the binding site
        properties (dic):
            radius: Cut-off distance(Amstrongs) around ligand atoms to consider a protein atom as a binding site atom. Default: 5
            identity_cluster: Minimal sequence identity (%) of the cluster members shared with the input_pdb. Options: '95', '90', '75' ,'50'. Default: 90
            max_num_ligands: Total number of superimposed ligands to be extracted from the identity cluster. For populated clusters, the restriction avoids to superimpose redundant structures. If 0, all ligands extracted will be considered. Default: 15.
            matrix (str): Substitution matrices for use in alignments. Available matrices are: 'benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65', 'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180', 'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure'. Default: blosum62
            gap_open (float): Gap open penalty. Default: -10.0
            gap_extend (float): Gap extend penalty. Default: -0.5
            trim_ends (True|False): Cut unaligned sequence ends. Default: True
    """

    def __init__(self, output_pdb_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        #i/o
        self.pdb_code             = properties.get('pdb_code','').lower()
        self.pdb_chain            = properties.get('pdb_chain','').upper()
        self.output_pdb_path      = output_pdb_path
        # options
        self.radius               = properties.get('radius', 5)
        self.identity_cluster     = str(properties.get('identity_cluster', '90'))
        self.max_num_ligands      = properties.get('max_num_ligands',15)
        matrix_name               = properties.get('matrix', "blosum62")
        self.gap_open             = properties.get('gap_open', -10.0)
        self.gap_extend           = properties.get('gap_extend', -0.5)
        trim_ends                 = properties.get('trim_ends', "True")
        self.path                 = properties.get('path', '')

        # eval booleans
        if trim_ends == "True": self.trim_ends=True; self.trim_ends=False
        # eval obj attributes
        self.matrix = getattr(Bio.SubsMat.MatrixInfo, matrix_name)

        # load modres and ions list
        self.modres   = self.__modres()
        self.ions     = self.__ions()


    def launch(self):
        """Launches the pipeline to find the binding site of the given structure
        """
        #out_log, err_log = fu.get_logs(path=self.path, mutation=self.mutation, step=self.step)

        ##
        ## Download reference PDB structure
        structure_name  = self.pdb_code + "_" + self.pdb_chain
        input_pdb_path  = self.path + "/" + structure_name + ".pdb"

        print "\nDownloading input PDB %s to %s" % (structure_name,input_pdb_path)
        pdb.MmbPdb().get_pdb(self.pdb_code,input_pdb_path,"filter=:" + self.pdb_chain)

        if not os.path.isfile(input_pdb_path):
            raise Exception('Failed to download input PDB file % from MMB PDB server' % input_pdb_path )

        ##
        ## Loading and parsing reference PDB structure
        print "\nParsing input PDB structure %s..." % input_pdb_path

        # Parse structure
        parser      = Bio.PDB.PDBParser()
        structPDB   = parser.get_structure(structure_name,input_pdb_path)
        if len(structPDB):
            structPDB = structPDB[0]

        # Use only one chain
        n_chains = structPDB.get_list()
        if len(n_chains) != 1:
            warnings.warn('More than one chain found in the input PDB structure. Using only the first chain to find the binding site')

        for struct_chain in structPDB.get_chains():
            if self.pdb_chain  != struct_chain.get_id():
                raise Exception('The input PDB %s does not contain the given chain %s. Found %s' % (input_pdb_path,self.pdb_chain,struct_chain.get_id()) )
            structPDB = struct_chain

        # Get AA sequence
        structPDB_seq   = self.__get_pdb_sequence(structPDB)
        if len(structPDB_seq) == 0:
                raise Exception('Cannot extract AA sequence from the input PDB structure {1}. Wrong format?'.format(input_pdb_path))

        else:
                print "    Found %s residues (chain %s)" % (len(structPDB_seq),self.pdb_chain)

        ##
        ## Getting PDBs for members of the identity cluster

        cluster_name = "cl-" + self.identity_cluster
        cluster_list = pdb.MmbPdb().get_cluster_pdb_code_chains(self.pdb_code+"_"+self.pdb_chain,cluster_name)

        if len(cluster_list) == 0:
                raise Exception('Found no similar structures into identity cluster %s for input PDB structure %s. Use a less restrictive identity cluster? Check available options for argument "identity_cluster"' % (cluster_name,structure_name))
        else:
                print "Found %s similar structures into the identity cluster %s%%" % (len(cluster_list),cluster_name)


        ##
        ## Loop each clusterPDB

        print "\nAligning cluster members ligands to input PDB structure, one by one..."

        clusterPDB_ligands_aligned = []
        clusterPDB_ligands_num     = 0


        for cluster_name in cluster_list:

                print "\n- Cluster member: %s" % cluster_name

                clusterPDB_path         = self.path + "/cluster_" + cluster_name + ".pdb"
                clusterPDB_aligned_path = self.path + "/cluster_" + cluster_name + "_aligned.pdb"

                clusterPDB_code, clusterPDB_chain = self.__get_pdb_code_from_filename(cluster_name)

                ##
                ## Donwloading and parsing cluster PDB structure

                # Download PDB from MMB server
                if os.path.exists(clusterPDB_path) == False:
                    print "    Downloading to  %s" % clusterPDB_path
                    pdb.MmbPdb().get_pdb(clusterPDB_code,clusterPDB_path,"filter=:" + clusterPDB_chain.upper())

                if os.path.exists(clusterPDB_path) == False:
                    raise Exception('Cluster member with pdbCode=%s and pdbChain=%s not found at the MMB PDB server. Ignoring this member.' % clusterPDB_code,clusterPDB_chain.upper())
                    continue;

                # Load and Parse PDB
                clusterPDB = {}
                clusterPDB = parser.get_structure(cluster_name, clusterPDB_path)[0]

                # Use only the fist chain
                for cluster_chain in clusterPDB.get_chains():
                    clusterPDB = cluster_chain

                # Looking for ligands
                clusterPDB_ligands = self.__get_ligand_residues(clusterPDB)
                if (len(clusterPDB_ligands)) == 0:
                    print "    No ligands found that could guide the binding site search. Ignoring this member.";
                    continue;

                # Selecting the largest ligand, if more than one
                lig_atoms_num     = 0
                clusterPDB_ligand = {}
                if len(clusterPDB_ligands) > 1:
                    for lig_res in clusterPDB_ligands:
                        lig_res_atoms_num = len(lig_res.get_list())
                        print "    Ligand found: %s  (%s atoms) " % (lig_res.get_resname(), lig_res_atoms_num);
                        if lig_res_atoms_num > lig_atoms_num:
                            clusterPDB_ligand = lig_res
                            lig_atoms_num     = lig_res_atoms_num
                else:
                    clusterPDB_ligand = clusterPDB_ligands[0]
                    lig_atoms_num = len(clusterPDB_ligands[0].get_list())

                print "    Member accepted. Valid ligand found: %s (%s atoms)" % (clusterPDB_ligand.get_resname(), lig_atoms_num)

                ##
                ## Mapping residues by sequence alignment to match structPDB-clusterPDB paired residues

                # Get AA sequence
                clusterPDB_seq = self.__get_pdb_sequence(clusterPDB)

                # Pairwise align
                aln, residue_map = self.__align_sequences(structPDB_seq,clusterPDB_seq)
                print "    Matching residues to input PDB structure. Alignment is:\n    %s" % aln[1];

                # Calculate (gapless) sequence identity
                seq_identity, gap_seq_identity = self.__calculate_alignment_identity(aln[0], aln[1])
                print "    Sequence identity (%%): %s" % seq_identity
                print "    Gap less identity (%%): %s" % gap_seq_identity


                ##
                ## Selecting aligned CA atoms from first model, first chain

                struct_atoms  = []
                cluster_atoms = []

                for struct_res in residue_map:
                    try:
                        cluster_atoms.append(clusterPDB[residue_map[struct_res]]['CA'])
                        struct_atoms.append(structPDB[struct_res]['CA'])
                    except KeyError:
                        print "Cannot find CA atom for residue %s  (input PDB  %s)" % (structPDB[struct_res],struct_res)
                        pass

                if len(cluster_atoms)==0:
                    raise Exception('Cannot find CA atoms (1st model, 1st chain) in cluster member {1} when aligning against {2}. Ignoring this member.'.format(clusterPDB_path,structure_name))
                else:
                    print "    Superimposing %s aligned protein residues" % len(cluster_atoms)

                #cluster_atoms = []
                #for cluster_res in clusterPDB:
                #    try:
                #        cluster_atoms.append(cluster_res['CA'])
                #    except KeyError:
                #        #do not append hets
                #       pass


                ##
                ## Align against input structure

                si = Bio.PDB.Superimposer()
                si.set_atoms(struct_atoms, cluster_atoms)
                si.apply(clusterPDB.get_atoms())
                print "    RMSD: %s" %si.rms

                # Save transformed structure (and ligand)
                clusterPDB_ligand_aligned = clusterPDB[clusterPDB_ligand.get_id()]
                print "    Saving transformed ligand coordinates"

                clusterPDB_ligands_aligned.append(clusterPDB_ligand_aligned)

                #iioo = Bio.PDB.PDBIO()
                #iioo.set_structure(clusterPDB_ligand_aligned)
                #iioo.save(self.path + "/cluster_" + cluster_name +  "_ligand.pdb")

                ##
                ##  Stop after n accepted cluster members

                clusterPDB_ligands_num += 1

                if clusterPDB_ligands_num > self.max_num_ligands:
                    break


        ##
        ## Select binding site atoms as those around cluster superimposed ligands

        print "\nDefining binding site residues as those %sA around the %s cluster superimposed ligands" % (self.radius,clusterPDB_ligands_num)

        # select Atoms from aligned ligands
        clusterPDB_ligands_aligned2      = [res for res in clusterPDB_ligands_aligned]
        clusterPDB_ligands_aligned_atoms = Bio.PDB.Selection.unfold_entities(clusterPDB_ligands_aligned2, 'A')

        # select Atoms from input PDB structure
        structPDB_atoms = [ atom for atom in structPDB.get_atoms()]

        # compute neighbors for aligned ligands in the input PDB structure
        structPDB_bs_residues_raw = {}
        structPDB_neighbors = Bio.PDB.NeighborSearch(structPDB_atoms)
        for ligand_atom in clusterPDB_ligands_aligned_atoms:
                # look for PDB atoms 5A around each ligand atom
                k_l = structPDB_neighbors.search(ligand_atom.coord, self.radius, 'R')
                for k in k_l:
                    structPDB_bs_residues_raw[k.get_id()] = k.get_full_id()

        ##
        ## Save binding site to PDB

        io = Bio.PDB.PDBIO()
        print "Writing binding site residues into %s" % self.output_pdb_path

        # unselect input PDB atoms not in binding site
        structPDB_bs_atoms = 0
        for res in list(structPDB):
            if res.get_id() not in structPDB_bs_residues_raw.keys():
                structPDB.detach_child(res.get_id())
            else:
                structPDB_bs_atoms += len(res.get_list())

        # write PDB file
        io.set_structure(structPDB)
        io.save(self.output_pdb_path)

        return 0


    #############
    #############


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


    def __get_ligand_residues(self,PDBchain,ignore_wats=True,ignore_small_molec=True, ignore_ions=True, ignore_modres=True):
        """
        Returns heteroatoms residues.
        Args:
            PDBchain (Bio.PDB.PDBParser chain object): PDB selection of the chain to be scanned
            ignore_wats (boolean): If True, water residues will be skipped, and not returned as ligand residues
            ignore_small_molec (boolean): If True, small ligands (< 5 atoms) will be skipped, and not returned as ligand residues
            ignore_ions (boolean): If True, ion residues will be skipped, and not returned as ligand residues
            ignore_modres (boolean): If True, modified aminoa acid residues will be skipped, and not returned as ligand residues
        """

        small_molec_atoms_min = 5
        ligands = []

        for res in PDBchain.get_residues():
            res_entity  = res.get_full_id()
            res_hetflag = res_entity[3][0]

            # skip aminoacids
            if res_hetflag == " ":
                continue

            # skip waters, if defined
            if res_hetflag == "W":
                if not ignore_wats:
                    ligands.append(res)
                continue
            # skip small_molec (< small_molec_atoms_min)
            #if ignore_small_molec:
            #    if len(res.get_list()) < small_molec_atoms_min:
            #            continue
            # skip ions
            if ignore_ions:
                if res.get_resname().strip() in self.ions:
                    continue
            # skip modres
            if ignore_modres:
                if res.get_resname().strip() in self.modres.keys():
                    continue

            # add as ligand
            ligands.append(res)

        return ligands

    def __get_binding_site_residues(self, struct_chain, ligands, distance):
        """
        Prints out residues located at a given distance from ligand
        """
        chains = model.child_dict
        for group in ligands.values():
            for ligand_res in group:
                for protein_res in struct_chain.child_list:
                    if protein_res not in group:
                        dist = residue_dist_to_ligand(protein_res, ligand_res)
                        if dist and dist < distance :
                            print protein_res.resname, protein_res.id[1], dist


    def __get_pdb_code_from_filename(self,filename):
        """
        Parses pdb_code (and pdb_chain if there) from a filename. RE *[pdbCode](-|.|_)([pdbChain])*
        """
        p = re.compile('([0-9][A-Za-z0-9]{3})[\W_]*([A-Za-z])*')
        m = p.match(filename)
        return (m.group(1),m.group(2))


    def __ions(self):
        return {
        'UNX',  #UNKNOWN ATOM OR ION
        'LI',  #LITHIUM ION
        'OH',  #HYDROXIDE ION
        'NH4',  #AMMONIUM ION
        'F',  #FLUORIDE ION
        'ND4',  #AMMONIUM CATION WITH D
        'NA',  #SODIUM ION
        'MG',  #MAGNESIUM ION
        'CYN',  #CYANIDE ION
        'AL',  #ALUMINUM ION
        '2FK',  #SUPEROXO ION
        'PER',  #PEROXIDE ION
        '3P8',  #methylammonium ion
        'CL',  #CHLORIDE ION
        'K',  #POTASSIUM ION
        'CA',  #CALCIUM ION
        'AZI',  #AZIDE ION
        'NO2',  #NITRITE ION
        '4TI',  #TITANIUM ION
        'V',  #VANADIUM ION
        'CR',  #CHROMIUM ION
        'MN',  #MANGANESE (II) ION
        'MN3',  #MANGANESE (III) ION
        'FE',  #FE (III) ION
        'FE2',  #FE (II) ION
        'SCN',  #THIOCYANATE ION
        '3NI',  #NICKEL (III) ION
        'NI',  #NICKEL (II) ION
        '3CO',  #COBALT (III) ION
        'CO',  #COBALT (II) ION
        'ACT',  #ACETATE ION
        'CO3',  #CARBONATE ION
        'BCT',  #BICARBONATE ION
        'NO3',  #NITRATE ION
        'CU',  #COPPER (II) ION
        'CU1',  #COPPER (I) ION
        'CU3',  #COPPER (III) ION
        'ZN',  #ZINC ION
        'BEF',  #BERYLLIUM TRIFLUORIDE ION
        'GA',  #GALLIUM (III) ION
        'MH2',  #MANGANESE ION
        'TMA',  #TETRAMETHYLAMMONIUM ION
        'BO4',  #BORATE ION
        'PO3',  #PHOSPHITE ION
        'BR',  #BROMIDE ION
        'SO3',  #SULFITE ION
        'LCO',  #CHLORATE ION
        'BF4',  #BERYLLIUM TETRAFLUORIDE ION
        'RB',  #RUBIDIUM ION
        'SR',  #STRONTIUM ION
        'OXL',  #OXALATE ION
        'Y1',  #YTTRIUM ION
        'YT3',  #YTTRIUM (III) ION
        'ZR',  #ZIRCONIUM ION
        'PO4',  #PHOSPHATE ION
        '4MO',  #MOLYBDENUM(IV) ION
        '6MO',  #MOLYBDENUM(VI) ION
        'PI',  #HYDROGENPHOSPHATE ION
        'SO4',  #SULFATE ION
        '2HP',  #DIHYDROGENPHOSPHATE ION
        'DMI',  #2
        'FPO',  #FLUORO-PHOSPHITE ION
        'VN3',  #VANADATE ION
        'LCP',  #PERCHLORATE ION
        '3MT',  #3-METHYLTHIAZOLIUM ION
        'HAI',  #CYCLOHEXYLAMMONIUM ION
        'RU',  #RUTHENIUM ION
        'MLI',  #MALONATE ION
        'TEA',  #TRIETHYLAMMONIUM ION
        'RH3',  #RHODIUM(III) ION
        'ALF',  #TETRAFLUOROALUMINATE ION
        'CHT',  #CHOLINE ION
        'SEK',  #SELENOCYANATE ION
        'PD',  #PALLADIUM ION
        'AG',  #SILVER ION
        'CD',  #CADMIUM ION
        'DTI',  #3
        'IN',  #INDIUM (III) ION
        'VO4',  #VANADATE ION
        'SB',  #ANTIMONY (III) ION
        'IOD',  #IODIDE ION
        'CON',  #COBALT TETRAAMMINE ION
        'CUA',  #DINUCLEAR COPPER ION
        'BSY',  #BISELENITE ION
        'NET',  #TETRAETHYLAMMONIUM ION
        'OAA',  #OXALOACETATE ION
        'CS',  #CESIUM ION
        'THE',  #THREONATE ION
        'CAC',  #CACODYLATE ION
        'BA',  #BARIUM ION
        'LA',  #LANTHANUM (III) ION
        'CE',  #CERIUM (III) ION
        'PR',  #PRASEODYMIUM ION
        'SE4',  #SELENATE ION
        'MOW',  #Oxo(sulfanyl)molybdenum(IV) ION
        'SM',  #SAMARIUM (III) ION
        'EU',  #EUROPIUM ION
        'EU3',  #EUROPIUM (III) ION
        'GD3',  #GADOLINIUM ION
        'TB',  #TERBIUM(III) ION
        'MOO',  #MOLYBDATE ION
        'SMO',  #DIOXOSULFIDOMOLYBDENUM(VI) ION
        'MOS',  #DIOXOTHIOMOLYBDENUM(VI) ION
        'DY',  #DYSPROSIUM ION
        'TCN',  #TETRACYANONICKELATE ION
        'EDR',  #EDROPHONIUM ION
        'ER3',  #ERBIUM (III) ION
        'TRA',  #ACONITATE ION
        'YB',  #YTTERBIUM (III) ION
        'YB2',  #YTTERBIUM (II) ION
        'LU',  #LUTETIUM (III) ION
        '1AL',  #ALLANTOATE ION
        'W',  #TUNGSTEN ION
        'ATH',  #4-HYDROXY-ACONITATE ION
        'FLC',  #CITRATE ANION
        'OS',  #OSMIUM ION
        'OS4',  #OSMIUM 4+ ION
        'T1A',  #TETRAETHYLARSONIUM ION
        'IR',  #IRIDIUM ION
        'IR3',  #IRIDIUM (III) ION
        'PT',  #PLATINUM (II) ION
        'PT4',  #PLATINUM (IV) ION
        'AU',  #GOLD ION
        'AU3',  #GOLD 3+ ION
        'HG',  #MERCURY (II) ION
        'NRU',  #RUTHENIUM (III) HEXAAMINE ION
        'TL',  #THALLIUM (I) ION
        'RHD',  #RHODIUM HEXAMINE ION
        'EUD',  #EUDESMANE CATION
        'PB',  #LEAD (II) ION
        'BS3',  #Bismuth(III) ION
        'PDV',  #Divanadate ion
        'MMC',  #METHYL MERCURY ION
        'EMC',  #ETHYL MERCURY ION
        'TH',  #THORIUM ION
        'TBA',  #TETRABUTYLAMMONIUM ION
        'AM',  #AMERICIUM ION
        '4PU',  #PLUTONIUM ION
        'PTN',  #PLATINUM TRIAMINE ION
        'ZCM',  #CURIUM ION
        'AUC',  #GOLD (I) CYANIDE ION
        'DSC',  #DODECANESULFONATE ION
        'CF',  #CALIFORNIUM ION
        'PBM',  #TRIMETHYL LEAD ION
        'DME',  #DECAMETHONIUM ION
        'MAC',  #MERCURY ACETATE ION
        'WO5',  #TUNGSTATE(VI) ION
        'IUM',  #URANYL (VI) ION
        'CUZ',  #(MU-4-SULFIDO)-TETRA-NUCLEAR COPPER ION
        'I3M'  #Tri-iodode Anion
        }

    def __modres(self):
        return {
        '0CS':'ALA', ##  0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
        '1AB':'PRO', ##  1AB PRO  1,4-DIDEOXY-1,4-IMINO-D-ARABINITOL
        '1LU':'LEU', ##  1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
        '1PA':'PHE', ##  1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
        '1TQ':'TRP', ##  1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
        '1TY':'TYR', ##  1TY TYR
        '23F':'PHE', ##  23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
        '23S':'TRP', ##  23S TRP  MODIFIED TRYPTOPHAN
        '2BU':'ALA', ##  2BU ADE
        '2ML':'LEU', ##  2ML LEU  2-METHYLLEUCINE
        '2MR':'ARG', ##  2MR ARG  N3, N4-DIMETHYLARGININE
        '2MT':'PRO', ##  2MT PRO
        '2OP':'ALA', ##  2OP (2S  2-HYDROXYPROPANAL
        '2TY':'TYR', ##  2TY TYR
        '32S':'TRP', ##  32S TRP  MODIFIED TRYPTOPHAN
        '32T':'TRP', ##  32T TRP  MODIFIED TRYPTOPHAN
        '3AH':'HIS', ##  3AH HIS
        '3MD':'ASP', ##  3MD ASP  2S,3S-3-METHYLASPARTIC ACID
        '3TY':'TYR', ##  3TY TYR  MODIFIED TYROSINE
        '4DP':'TRP', ##  4DP TRP
        '4F3':'ALA', ##  4F3 ALA  CYCLIZED
        '4FB':'PRO', ##  4FB PRO  (4S)-4-FLUORO-L-PROLINE
        '4FW':'TRP', ##  4FW TRP  4-FLUOROTRYPTOPHANE
        '4HT':'TRP', ##  4HT TRP  4-HYDROXYTRYPTOPHAN
        '4IN':'TRP', ##  4IN TRP  4-AMINO-L-TRYPTOPHAN
        '4PH':'PHE', ##  4PH PHE  4-METHYL-L-PHENYLALANINE
        '5CS':'CYS', ##  5CS CYS
        '6CL':'LYS', ##  6CL LYS  6-CARBOXYLYSINE
        '6CW':'TRP', ##  6CW TRP  6-CHLORO-L-TRYPTOPHAN
        'A0A':'ASP', ##  A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
        'AA4':'ALA', ##  AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
        'AAR':'ARG', ##  AAR ARG  ARGININEAMIDE
        'AB7':'GLU', ##  AB7 GLU  ALPHA-AMINOBUTYRIC ACID
        'ABA':'ALA', ##  ABA ALA  ALPHA-AMINOBUTYRIC ACID
        'ACB':'ASP', ##  ACB ASP  3-METHYL-ASPARTIC ACID
        'ACL':'ARG', ##  ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
        'ACY':'GLY', ##  ACY GLY  POST-TRANSLATIONAL MODIFICATION
        'AEI':'THR', ##  AEI THR  ACYLATED THR
        'AFA':'ASN', ##  AFA ASN  N-[7-METHYL-OCT-2,4-DIENOYL]ASPARAGINE
        'AGM':'ARG', ##  AGM ARG  4-METHYL-ARGININE
        'AGT':'CYS', ##  AGT CYS  AGMATINE-CYSTEINE ADDUCT
        'AHB':'ASN', ##  AHB ASN  BETA-HYDROXYASPARAGINE
        'AHO':'ALA', ##  AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
        'AHP':'ALA', ##  AHP ALA  2-AMINO-HEPTANOIC ACID
        'AIB':'ALA', ##  AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
        'AKL':'ASP', ##  AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
        'ALA':'ALA', ##  ALA ALA
        'ALC':'ALA', ##  ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
        'ALG':'ARG', ##  ALG ARG  GUANIDINOBUTYRYL GROUP
        'ALM':'ALA', ##  ALM ALA  1-METHYL-ALANINAL
        'ALN':'ALA', ##  ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
        'ALO':'THR', ##  ALO THR  ALLO-THREONINE
        'ALS':'ALA', ##  ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
        'ALT':'ALA', ##  ALT ALA  THIOALANINE
        'ALY':'LYS', ##  ALY LYS  N(6)-ACETYLLYSINE
        'AME':'MET', ##  AME MET  ACETYLATED METHIONINE
        'AP7':'ALA', ##  AP7 ADE
        'APH':'ALA', ##  APH ALA  P-AMIDINOPHENYL-3-ALANINE
        'API':'LYS', ##  API LYS  2,6-DIAMINOPIMELIC ACID
        'APK':'LYS', ##  APK LYS
        'AR2':'ARG', ##  AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
        'AR4':'GLU', ##  AR4 GLU
        'ARG':'ARG', ##  ARG ARG
        'ARM':'ARG', ##  ARM ARG  DEOXY-METHYL-ARGININE
        'ARO':'ARG', ##  ARO ARG  C-GAMMA-HYDROXY ARGININE
        'ASA':'ASP', ##  ASA ASP  ASPARTIC ALDEHYDE
        'ASB':'ASP', ##  ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
        'ASI':'ASP', ##  ASI ASP  L-ISO-ASPARTATE
        'ASK':'ASP', ##  ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
        'ASL':'ASP', ##  ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
        'ASN':'ASN', ##  ASN ASN
        'ASP':'ASP', ##  ASP ASP
        'AYA':'ALA', ##  AYA ALA  N-ACETYLALANINE
        'AYG':'ALA', ##  AYG ALA
        'AZK':'LYS', ##  AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
        'B2A':'ALA', ##  B2A ALA  ALANINE BORONIC ACID
        'B2F':'PHE', ##  B2F PHE  PHENYLALANINE BORONIC ACID
        'B2I':'ILE', ##  B2I ILE  ISOLEUCINE BORONIC ACID
        'B2V':'VAL', ##  B2V VAL  VALINE BORONIC ACID
        'B3A':'ALA', ##  B3A ALA  (3S)-3-AMINOBUTANOIC ACID
        'B3D':'ASP', ##  B3D ASP  3-AMINOPENTANEDIOIC ACID
        'B3E':'GLU', ##  B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
        'B3K':'LYS', ##  B3K LYS  (3S)-3,7-DIAMINOHEPTANOIC ACID
        'B3S':'SER', ##  B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
        'B3X':'ASN', ##  B3X ASN  (3S)-3,5-DIAMINO-5-OXOPENTANOIC ACID
        'B3Y':'TYR', ##  B3Y TYR
        'BAL':'ALA', ##  BAL ALA  BETA-ALANINE
        'BBC':'CYS', ##  BBC CYS
        'BCS':'CYS', ##  BCS CYS  BENZYLCYSTEINE
        'BCX':'CYS', ##  BCX CYS  BETA-3-CYSTEINE
        'BFD':'ASP', ##  BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
        'BG1':'SER', ##  BG1 SER
        'BHD':'ASP', ##  BHD ASP  BETA-HYDROXYASPARTIC ACID
        'BIF':'PHE', ##  BIF PHE
        'BLE':'LEU', ##  BLE LEU  LEUCINE BORONIC ACID
        'BLY':'LYS', ##  BLY LYS  LYSINE BORONIC ACID
        'BMT':'THR', ##  BMT THR
        'BNN':'ALA', ##  BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
        'BOR':'ARG', ##  BOR ARG
        'BPE':'CYS', ##  BPE CYS
        'BTR':'TRP', ##  BTR TRP  6-BROMO-TRYPTOPHAN
        'BUC':'CYS', ##  BUC CYS  S,S-BUTYLTHIOCYSTEINE
        'BUG':'LEU', ##  BUG LEU  TERT-LEUCYL AMINE
        'C12':'ALA', ##  C12 ALA
        'C1X':'LYS', ##  C1X LYS  MODIFIED LYSINE
        'C3Y':'CYS', ##  C3Y CYS  MODIFIED CYSTEINE
        'C5C':'CYS', ##  C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
        'C6C':'CYS', ##  C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
        'C99':'ALA', ##  C99 ALA
        'CAB':'ALA', ##  CAB ALA  4-CARBOXY-4-AMINOBUTANAL
        'CAF':'CYS', ##  CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
        'CAS':'CYS', ##  CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
        'CCS':'CYS', ##  CCS CYS  CARBOXYMETHYLATED CYSTEINE
        'CGU':'GLU', ##  CGU GLU  CARBOXYLATION OF THE CG ATOM
        'CH6':'ALA', ##  CH6 ALA
        'CH7':'ALA', ##  CH7 ALA
        'CHG':'GLY', ##  CHG GLY  CYCLOHEXYL GLYCINE
        'CHP':'GLY', ##  CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
        'CHS':'PHE', ##  CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
        'CIR':'ARG', ##  CIR ARG  CITRULLINE
        'CLB':'ALA', ##  CLB ALA
        'CLD':'ALA', ##  CLD ALA
        'CLE':'LEU', ##  CLE LEU  LEUCINE AMIDE
        'CLG':'LYS', ##  CLG LYS
        'CLH':'LYS', ##  CLH LYS
        'CLV':'ALA', ##  CLV ALA
        'CME':'CYS', ##  CME CYS  MODIFIED CYSTEINE
        'CML':'CYS', ##  CML CYS
        'CMT':'CYS', ##  CMT CYS  O-METHYLCYSTEINE
        'CQR':'ALA', ##  CQR ALA
        'CR2':'ALA', ##  CR2 ALA  POST-TRANSLATIONAL MODIFICATION
        'CR5':'ALA', ##  CR5 ALA
        'CR7':'ALA', ##  CR7 ALA
        'CR8':'ALA', ##  CR8 ALA
        'CRK':'ALA', ##  CRK ALA
        'CRO':'THR', ##  CRO THR  CYCLIZED
        'CRQ':'TYR', ##  CRQ TYR
        'CRW':'ALA', ##  CRW ALA
        'CRX':'ALA', ##  CRX ALA
        'CS1':'CYS', ##  CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
        'CS3':'CYS', ##  CS3 CYS
        'CS4':'CYS', ##  CS4 CYS
        'CSA':'CYS', ##  CSA CYS  S-ACETONYLCYSTEIN
        'CSB':'CYS', ##  CSB CYS  CYS BOUND TO LEAD ION
        'CSD':'CYS', ##  CSD CYS  3-SULFINOALANINE
        'CSE':'CYS', ##  CSE CYS  SELENOCYSTEINE
        'CSI':'ALA', ##  CSI ALA
        'CSO':'CYS', ##  CSO CYS  INE S-HYDROXYCYSTEINE
        'CSR':'CYS', ##  CSR CYS  S-ARSONOCYSTEINE
        'CSS':'CYS', ##  CSS CYS  1,3-THIAZOLE-4-CARBOXYLIC ACID
        'CSU':'CYS', ##  CSU CYS  CYSTEINE-S-SULFONIC ACID
        'CSW':'CYS', ##  CSW CYS  CYSTEINE-S-DIOXIDE
        'CSX':'CYS', ##  CSX CYS  OXOCYSTEINE
        'CSY':'ALA', ##  CSY ALA  MODIFIED TYROSINE COMPLEX
        'CSZ':'CYS', ##  CSZ CYS  S-SELANYL CYSTEINE
        'CTH':'THR', ##  CTH THR  4-CHLOROTHREONINE
        'CWR':'ALA', ##  CWR ALA
        'CXM':'MET', ##  CXM MET  N-CARBOXYMETHIONINE
        'CY0':'CYS', ##  CY0 CYS  MODIFIED CYSTEINE
        'CY1':'CYS', ##  CY1 CYS  ACETAMIDOMETHYLCYSTEINE
        'CY3':'CYS', ##  CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
        'CY4':'CYS', ##  CY4 CYS  S-BUTYRYL-CYSTEIN
        'CY7':'CYS', ##  CY7 CYS  MODIFIED CYSTEINE
        'CYD':'CYS', ##  CYD CYS
        'CYF':'CYS', ##  CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
        'CYG':'CYS', ##  CYG CYS
        'CYJ':'LYS', ##  CYJ LYS  MODIFIED LYSINE
        'CYQ':'CYS', ##  CYQ CYS
        'CYR':'CYS', ##  CYR CYS
        'CYS':'CYS', ##  CYS CYS
        'CZ2':'CYS', ##  CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
        'CZZ':'CYS', ##  CZZ CYS  THIARSAHYDROXY-CYSTEINE
        'DA2':'ARG', ##  DA2 ARG  MODIFIED ARGININE
        'DAB':'ALA', ##  DAB ALA  2,4-DIAMINOBUTYRIC ACID
        'DAH':'PHE', ##  DAH PHE  3,4-DIHYDROXYDAHNYLALANINE
        'DAL':'ALA', ##  DAL ALA  D-ALANINE
        'DAM':'ALA', ##  DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
        'DAR':'ARG', ##  DAR ARG  D-ARGININE
        'DAS':'ASP', ##  DAS ASP  D-ASPARTIC ACID
        'DBU':'ALA', ##  DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
        'DBY':'TYR', ##  DBY TYR  3,5 DIBROMOTYROSINE
        'DBZ':'ALA', ##  DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
        'DCL':'LEU', ##  DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
        'DCY':'CYS', ##  DCY CYS  D-CYSTEINE
        'DDE':'HIS', ##  DDE HIS
        'DGL':'GLU', ##  DGL GLU  D-GLU
        'DGN':'GLN', ##  DGN GLN  D-GLUTAMINE
        'DHA':'ALA', ##  DHA ALA  2-AMINO-ACRYLIC ACID
        'DHI':'HIS', ##  DHI HIS  D-HISTIDINE
        'DHL':'SER', ##  DHL SER  POST-TRANSLATIONAL MODIFICATION
        'DIL':'ILE', ##  DIL ILE  D-ISOLEUCINE
        'DIV':'VAL', ##  DIV VAL  D-ISOVALINE
        'DLE':'LEU', ##  DLE LEU  D-LEUCINE
        'DLS':'LYS', ##  DLS LYS  DI-ACETYL-LYSINE
        'DLY':'LYS', ##  DLY LYS  D-LYSINE
        'DMH':'ASN', ##  DMH ASN  N4,N4-DIMETHYL-ASPARAGINE
        'DMK':'ASP', ##  DMK ASP  DIMETHYL ASPARTIC ACID
        'DNE':'LEU', ##  DNE LEU  D-NORLEUCINE
        'DNG':'LEU', ##  DNG LEU  N-FORMYL-D-NORLEUCINE
        'DNL':'LYS', ##  DNL LYS  6-AMINO-HEXANAL
        'DNM':'LEU', ##  DNM LEU  D-N-METHYL NORLEUCINE
        'DPH':'PHE', ##  DPH PHE  DEAMINO-METHYL-PHENYLALANINE
        'DPL':'PRO', ##  DPL PRO  4-OXOPROLINE
        'DPN':'PHE', ##  DPN PHE  D-CONFIGURATION
        'DPP':'ALA', ##  DPP ALA  DIAMMINOPROPANOIC ACID
        'DPQ':'TYR', ##  DPQ TYR  TYROSINE DERIVATIVE
        'DPR':'PRO', ##  DPR PRO  D-PROLINE
        'DSE':'SER', ##  DSE SER  D-SERINE N-METHYLATED
        'DSG':'ASN', ##  DSG ASN  D-ASPARAGINE
        'DSN':'SER', ##  DSN SER  D-SERINE
        'DTH':'THR', ##  DTH THR  D-THREONINE
        'DTR':'TRP', ##  DTR TRP  D-TRYPTOPHAN
        'DTY':'TYR', ##  DTY TYR  D-TYROSINE
        'DVA':'VAL', ##  DVA VAL  D-VALINE
        'DYG':'ALA', ##  DYG ALA
        'DYS':'CYS', ##  DYS CYS
        'EFC':'CYS', ##  EFC CYS  S,S-(2-FLUOROETHYL)THIOCYSTEINE
        'ESB':'TYR', ##  ESB TYR
        'ESC':'MET', ##  ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
        'FCL':'PHE', ##  FCL PHE  3-CHLORO-L-PHENYLALANINE
        'FGL':'ALA', ##  FGL ALA  2-AMINOPROPANEDIOIC ACID
        'FGP':'SER', ##  FGP SER
        'FHL':'LYS', ##  FHL LYS  MODIFIED LYSINE
        'FLE':'LEU', ##  FLE LEU  FUROYL-LEUCINE
        'FLT':'TYR', ##  FLT TYR  FLUOROMALONYL TYROSINE
        'FME':'MET', ##  FME MET  FORMYL-METHIONINE
        'FOE':'CYS', ##  FOE CYS
        'FOG':'PHE', ##  FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
        'FOR':'MET', ##  FOR MET
        'FRF':'PHE', ##  FRF PHE  PHE FOLLOWED BY REDUCED PHE
        'FTR':'TRP', ##  FTR TRP  FLUOROTRYPTOPHANE
        'FTY':'TYR', ##  FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
        'GHG':'GLN', ##  GHG GLN  GAMMA-HYDROXY-GLUTAMINE
        'GHP':'GLY', ##  GHP GLY  4-HYDROXYPHENYLGLYCINE
        'GL3':'GLY', ##  GL3 GLY  POST-TRANSLATIONAL MODIFICATION
        'GLH':'GLN', ##  GLH GLN
        'GLN':'GLN', ##  GLN GLN
        'GLU':'GLU', ##  GLU GLU
        'GLY':'GLY', ##  GLY GLY
        'GLZ':'GLY', ##  GLZ GLY  AMINO-ACETALDEHYDE
        'GMA':'GLU', ##  GMA GLU  1-AMIDO-GLUTAMIC ACID
        'GMU':'ALA', ##  GMU 5MU
        'GPL':'LYS', ##  GPL LYS  LYSINE GUANOSINE-5'-MONOPHOSPHATE
        'GT9':'CYS', ##  GT9 CYS  SG ALKYLATED
        'GVL':'SER', ##  GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
        'GYC':'CYS', ##  GYC CYS
        'GYS':'GLY', ##  GYS GLY
        'H5M':'PRO', ##  H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
        'HHK':'ALA', ##  HHK ALA  (2S)-2,8-DIAMINOOCTANOIC ACID
        'HIA':'HIS', ##  HIA HIS  L-HISTIDINE AMIDE
        'HIC':'HIS', ##  HIC HIS  4-METHYL-HISTIDINE
        'HIP':'HIS', ##  HIP HIS  ND1-PHOSPHONOHISTIDINE
        'HIQ':'HIS', ##  HIQ HIS  MODIFIED HISTIDINE
        'HIS':'HIS', ##  HIS HIS
        'HLU':'LEU', ##  HLU LEU  BETA-HYDROXYLEUCINE
        'HMF':'ALA', ##  HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
        'HMR':'ARG', ##  HMR ARG  BETA-HOMOARGININE
        'HPE':'PHE', ##  HPE PHE  HOMOPHENYLALANINE
        'HPH':'PHE', ##  HPH PHE  PHENYLALANINOL GROUP
        'HPQ':'PHE', ##  HPQ PHE  HOMOPHENYLALANINYLMETHANE
        'HRG':'ARG', ##  HRG ARG  L-HOMOARGININE
        'HSE':'SER', ##  HSE SER  L-HOMOSERINE
        'HSL':'SER', ##  HSL SER  HOMOSERINE LACTONE
        'HSO':'HIS', ##  HSO HIS  HISTIDINOL
        'HTI':'CYS', ##  HTI CYS
        'HTR':'TRP', ##  HTR TRP  BETA-HYDROXYTRYPTOPHANE
        'HY3':'PRO', ##  HY3 PRO  3-HYDROXYPROLINE
        'HYP':'PRO', ##  HYP PRO  4-HYDROXYPROLINE
        'IAM':'ALA', ##  IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
        'IAS':'ASP', ##  IAS ASP  ASPARTYL GROUP
        'IGL':'ALA', ##  IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
        'IIL':'ILE', ##  IIL ILE  ISO-ISOLEUCINE
        'ILE':'ILE', ##  ILE ILE
        'ILG':'GLU', ##  ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
        'ILX':'ILE', ##  ILX ILE  4,5-DIHYDROXYISOLEUCINE
        'IML':'ILE', ##  IML ILE  N-METHYLATED
        'IPG':'GLY', ##  IPG GLY  N-ISOPROPYL GLYCINE
        'IT1':'LYS', ##  IT1 LYS
        'IYR':'TYR', ##  IYR TYR  3-IODO-TYROSINE
        'KCX':'LYS', ##  KCX LYS  CARBAMOYLATED LYSINE
        'KGC':'LYS', ##  KGC LYS
        'KOR':'CYS', ##  KOR CYS  MODIFIED CYSTEINE
        'KST':'LYS', ##  KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
        'KYN':'ALA', ##  KYN ALA  KYNURENINE
        'LA2':'LYS', ##  LA2 LYS
        'LAL':'ALA', ##  LAL ALA  N,N-DIMETHYL-L-ALANINE
        'LCK':'LYS', ##  LCK LYS
        'LCX':'LYS', ##  LCX LYS  CARBAMYLATED LYSINE
        'LDH':'LYS', ##  LDH LYS  N~6~-ETHYL-L-LYSINE
        'LED':'LEU', ##  LED LEU  POST-TRANSLATIONAL MODIFICATION
        'LEF':'LEU', ##  LEF LEU  2-5-FLUOROLEUCINE
        'LET':'LYS', ##  LET LYS  ODIFIED LYSINE
        'LEU':'LEU', ##  LEU LEU
        'LLP':'LYS', ##  LLP LYS
        'LLY':'LYS', ##  LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
        'LME':'GLU', ##  LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
        'LNT':'LEU', ##  LNT LEU
        'LPD':'PRO', ##  LPD PRO  L-PROLINAMIDE
        'LSO':'LYS', ##  LSO LYS  MODIFIED LYSINE
        'LYM':'LYS', ##  LYM LYS  DEOXY-METHYL-LYSINE
        'LYN':'LYS', ##  LYN LYS  2,6-DIAMINO-HEXANOIC ACID AMIDE
        'LYP':'LYS', ##  LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
        'LYR':'LYS', ##  LYR LYS  MODIFIED LYSINE
        'LYS':'LYS', ##  LYS LYS
        'LYX':'LYS', ##  LYX LYS  N''-(2-COENZYME A)-PROPANOYL-LYSINE
        'LYZ':'LYS', ##  LYZ LYS  5-HYDROXYLYSINE
        'M0H':'CYS', ##  M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
        'M2L':'LYS', ##  M2L LYS
        'M3L':'LYS', ##  M3L LYS  N-TRIMETHYLLYSINE
        'MAA':'ALA', ##  MAA ALA  N-METHYLALANINE
        'MAI':'ARG', ##  MAI ARG  DEOXO-METHYLARGININE
        'MBQ':'TYR', ##  MBQ TYR
        'MC1':'SER', ##  MC1 SER  METHICILLIN ACYL-SERINE
        'MCL':'LYS', ##  MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
        'MCS':'CYS', ##  MCS CYS  MALONYLCYSTEINE
        'MDO':'ALA', ##  MDO ALA
        'MEA':'PHE', ##  MEA PHE  N-METHYLPHENYLALANINE
        'MEG':'GLU', ##  MEG GLU  (2S,3R)-3-METHYL-GLUTAMIC ACID
        'MEN':'ASN', ##  MEN ASN  GAMMA METHYL ASPARAGINE
        'MET':'MET', ##  MET MET
        'MEU':'GLY', ##  MEU GLY  O-METHYL-GLYCINE
        'MFC':'ALA', ##  MFC ALA  CYCLIZED
        'MGG':'ARG', ##  MGG ARG  MODIFIED D-ARGININE
        'MGN':'GLN', ##  MGN GLN  2-METHYL-GLUTAMINE
        'MHL':'LEU', ##  MHL LEU  N-METHYLATED, HYDROXY
        'MHO':'MET', ##  MHO MET  POST-TRANSLATIONAL MODIFICATION
        'MHS':'HIS', ##  MHS HIS  1-N-METHYLHISTIDINE
        'MIS':'SER', ##  MIS SER  MODIFIED SERINE
        'MLE':'LEU', ##  MLE LEU  N-METHYLATED
        'MLL':'LEU', ##  MLL LEU  METHYL L-LEUCINATE
        'MLY':'LYS', ##  MLY LYS  METHYLATED LYSINE
        'MLZ':'LYS', ##  MLZ LYS  N-METHYL-LYSINE
        'MME':'MET', ##  MME MET  N-METHYL METHIONINE
        'MNL':'LEU', ##  MNL LEU  4,N-DIMETHYLNORLEUCINE
        'MNV':'VAL', ##  MNV VAL  N-METHYL-C-AMINO VALINE
        'MPQ':'GLY', ##  MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
        'MSA':'GLY', ##  MSA GLY  (2-S-METHYL) SARCOSINE
        'MSE':'MET', ##  MSE MET  ELENOMETHIONINE
        'MSO':'MET', ##  MSO MET  METHIONINE SULFOXIDE
        'MTY':'PHE', ##  MTY PHE  3-HYDROXYPHENYLALANINE
        'MVA':'VAL', ##  MVA VAL  N-METHYLATED
        'N10':'SER', ##  N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
        'NAL':'ALA', ##  NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
        'NAM':'ALA', ##  NAM ALA  NAM NAPTHYLAMINOALANINE
        'NBQ':'TYR', ##  NBQ TYR
        'NC1':'SER', ##  NC1 SER  NITROCEFIN ACYL-SERINE
        'NCB':'ALA', ##  NCB ALA  CHEMICAL MODIFICATION
        'NEP':'HIS', ##  NEP HIS  N1-PHOSPHONOHISTIDINE
        'NFA':'PHE', ##  NFA PHE  MODIFIED PHENYLALANINE
        'NIY':'TYR', ##  NIY TYR  META-NITRO-TYROSINE
        'NLE':'LEU', ##  NLE LEU  NORLEUCINE
        'NLN':'LEU', ##  NLN LEU  NORLEUCINE AMIDE
        'NLO':'LEU', ##  NLO LEU  O-METHYL-L-NORLEUCINE
        'NMC':'GLY', ##  NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
        'NMM':'ARG', ##  NMM ARG  MODIFIED ARGININE
        'NPH':'CYS', ##  NPH CYS
        'NRQ':'ALA', ##  NRQ ALA
        'NVA':'VAL', ##  NVA VAL  NORVALINE
        'NYC':'ALA', ##  NYC ALA
        'NYS':'CYS', ##  NYS CYS
        'NZH':'HIS', ##  NZH HIS
        'OAS':'SER', ##  OAS SER  O-ACETYLSERINE
        'OBS':'LYS', ##  OBS LYS  MODIFIED LYSINE
        'OCS':'CYS', ##  OCS CYS  CYSTEINE SULFONIC ACID
        'OCY':'CYS', ##  OCY CYS  HYDROXYETHYLCYSTEINE
        'OHI':'HIS', ##  OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
        'OHS':'ASP', ##  OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
        'OLT':'THR', ##  OLT THR  O-METHYL-L-THREONINE
        'OMT':'MET', ##  OMT MET  METHIONINE SULFONE
        'OPR':'ARG', ##  OPR ARG  C-(3-OXOPROPYL)ARGININE
        'ORN':'ALA', ##  ORN ALA  ORNITHINE
        'ORQ':'ARG', ##  ORQ ARG  N~5~-ACETYL-L-ORNITHINE
        'OSE':'SER', ##  OSE SER  O-SULFO-L-SERINE
        'OTY':'TYR', ##  OTY TYR
        'OXX':'ASP', ##  OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
        'P1L':'CYS', ##  P1L CYS  S-PALMITOYL CYSTEINE
        'P2Y':'PRO', ##  P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
        'PAQ':'TYR', ##  PAQ TYR  SEE REMARK 999
        'PAT':'TRP', ##  PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
        'PBB':'CYS', ##  PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
        'PBF':'PHE', ##  PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
        'PCA':'PRO', ##  PCA PRO  5-OXOPROLINE
        'PCS':'PHE', ##  PCS PHE  PHENYLALANYLMETHYLCHLORIDE
        'PEC':'CYS', ##  PEC CYS  S,S-PENTYLTHIOCYSTEINE
        'PF5':'PHE', ##  PF5 PHE  2,3,4,5,6-PENTAFLUORO-L-PHENYLALANINE
        'PFF':'PHE', ##  PFF PHE  4-FLUORO-L-PHENYLALANINE
        'PG1':'SER', ##  PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
        'PG9':'GLY', ##  PG9 GLY  D-PHENYLGLYCINE
        'PHA':'PHE', ##  PHA PHE  PHENYLALANINAL
        'PHD':'ASP', ##  PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
        'PHE':'PHE', ##  PHE PHE
        'PHI':'PHE', ##  PHI PHE  IODO-PHENYLALANINE
        'PHL':'PHE', ##  PHL PHE  L-PHENYLALANINOL
        'PHM':'PHE', ##  PHM PHE  PHENYLALANYLMETHANE
        'PIA':'ALA', ##  PIA ALA  FUSION OF ALA 65, TYR 66, GLY 67
        'PLE':'LEU', ##  PLE LEU  LEUCINE PHOSPHINIC ACID
        'PM3':'PHE', ##  PM3 PHE
        'POM':'PRO', ##  POM PRO  CIS-5-METHYL-4-OXOPROLINE
        'PPH':'LEU', ##  PPH LEU  PHENYLALANINE PHOSPHINIC ACID
        'PPN':'PHE', ##  PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
        'PR3':'CYS', ##  PR3 CYS  INE DTT-CYSTEINE
        'PRO':'PRO', ##  PRO PRO
        'PRQ':'PHE', ##  PRQ PHE  PHENYLALANINE
        'PRR':'ALA', ##  PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
        'PRS':'PRO', ##  PRS PRO  THIOPROLINE
        'PSA':'PHE', ##  PSA PHE
        'PSH':'HIS', ##  PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
        'PTH':'TYR', ##  PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
        'PTM':'TYR', ##  PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
        'PTR':'TYR', ##  PTR TYR  O-PHOSPHOTYROSINE
        'PYA':'ALA', ##  PYA ALA  3-(1,10-PHENANTHROL-2-YL)-L-ALANINE
        'PYC':'ALA', ##  PYC ALA  PYRROLE-2-CARBOXYLATE
        'PYR':'SER', ##  PYR SER  CHEMICALLY MODIFIED
        'PYT':'ALA', ##  PYT ALA  MODIFIED ALANINE
        'PYX':'CYS', ##  PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
        'R1A':'CYS', ##  R1A CYS
        'R1B':'CYS', ##  R1B CYS
        'R1F':'CYS', ##  R1F CYS
        'R7A':'CYS', ##  R7A CYS
        'RC7':'ALA', ##  RC7 ALA
        'RCY':'CYS', ##  RCY CYS
        'S1H':'SER', ##  S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
        'SAC':'SER', ##  SAC SER  N-ACETYL-SERINE
        'SAH':'CYS', ##  SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
        'SAR':'GLY', ##  SAR GLY  SARCOSINE
        'SBD':'SER', ##  SBD SER
        'SBG':'SER', ##  SBG SER  MODIFIED SERINE
        'SBL':'SER', ##  SBL SER
        'SC2':'CYS', ##  SC2 CYS  N-ACETYL-L-CYSTEINE
        'SCH':'CYS', ##  SCH CYS  S-METHYL THIOCYSTEINE GROUP
        'SCS':'CYS', ##  SCS CYS  MODIFIED CYSTEINE
        'SCY':'CYS', ##  SCY CYS  CETYLATED CYSTEINE
        'SDP':'SER', ##  SDP SER
        'SEB':'SER', ##  SEB SER  O-BENZYLSULFONYL-SERINE
        'SEC':'ALA', ##  SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
        'SEL':'SER', ##  SEL SER  2-AMINO-1,3-PROPANEDIOL
        'SEP':'SER', ##  SEP SER  E PHOSPHOSERINE
        'SER':'SER', ##  SER SER
        'SET':'SER', ##  SET SER  AMINOSERINE
        'SGB':'SER', ##  SGB SER  MODIFIED SERINE
        'SGR':'SER', ##  SGR SER  MODIFIED SERINE
        'SHC':'CYS', ##  SHC CYS  S-HEXYLCYSTEINE
        'SHP':'GLY', ##  SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
        'SIC':'ALA', ##  SIC ALA
        'SLZ':'LYS', ##  SLZ LYS  L-THIALYSINE
        'SMC':'CYS', ##  SMC CYS  POST-TRANSLATIONAL MODIFICATION
        'SME':'MET', ##  SME MET  METHIONINE SULFOXIDE
        'SMF':'PHE', ##  SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
        'SNC':'CYS', ##  SNC CYS  S-NITROSO CYSTEINE
        'SNN':'ASP', ##  SNN ASP  POST-TRANSLATIONAL MODIFICATION
        'SOC':'CYS', ##  SOC CYS  DIOXYSELENOCYSTEINE
        'SOY':'SER', ##  SOY SER  OXACILLOYL-ACYLATED SERINE
        'SUI':'ALA', ##  SUI ALA
        'SUN':'SER', ##  SUN SER  TABUN CONJUGATED SERINE
        'SVA':'SER', ##  SVA SER  SERINE VANADATE
        'SVV':'SER', ##  SVV SER  MODIFIED SERINE
        'SVX':'SER', ##  SVX SER  MODIFIED SERINE
        'SVY':'SER', ##  SVY SER  MODIFIED SERINE
        'SVZ':'SER', ##  SVZ SER  MODIFIED SERINE
        'SXE':'SER', ##  SXE SER  MODIFIED SERINE
        'TBG':'GLY', ##  TBG GLY  T-BUTYL GLYCINE
        'TBM':'THR', ##  TBM THR
        'TCQ':'TYR', ##  TCQ TYR  MODIFIED TYROSINE
        'TEE':'CYS', ##  TEE CYS  POST-TRANSLATIONAL MODIFICATION
        'TH5':'THR', ##  TH5 THR  O-ACETYL-L-THREONINE
        'THC':'THR', ##  THC THR  N-METHYLCARBONYLTHREONINE
        'THR':'THR', ##  THR THR
        'TIH':'ALA', ##  TIH ALA  BETA(2-THIENYL)ALANINE
        'TMD':'THR', ##  TMD THR  N-METHYLATED, EPSILON C ALKYLATED
        'TNB':'CYS', ##  TNB CYS  S-(2,3,6-TRINITROPHENYL)CYSTEINE
        'TOX':'TRP', ##  TOX TRP
        'TPL':'TRP', ##  TPL TRP  TRYTOPHANOL
        'TPO':'THR', ##  TPO THR  HOSPHOTHREONINE
        'TPQ':'ALA', ##  TPQ ALA  2,4,5-TRIHYDROXYPHENYLALANINE
        'TQQ':'TRP', ##  TQQ TRP
        'TRF':'TRP', ##  TRF TRP  N1-FORMYL-TRYPTOPHAN
        'TRN':'TRP', ##  TRN TRP  AZA-TRYPTOPHAN
        'TRO':'TRP', ##  TRO TRP  2-HYDROXY-TRYPTOPHAN
        'TRP':'TRP', ##  TRP TRP
        'TRQ':'TRP', ##  TRQ TRP
        'TRW':'TRP', ##  TRW TRP
        'TRX':'TRP', ##  TRX TRP  6-HYDROXYTRYPTOPHAN
        'TTQ':'TRP', ##  TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
        'TTS':'TYR', ##  TTS TYR
        'TY2':'TYR', ##  TY2 TYR  3-AMINO-L-TYROSINE
        'TY3':'TYR', ##  TY3 TYR  3-HYDROXY-L-TYROSINE
        'TYB':'TYR', ##  TYB TYR  TYROSINAL
        'TYC':'TYR', ##  TYC TYR  L-TYROSINAMIDE
        'TYI':'TYR', ##  TYI TYR  3,5-DIIODOTYROSINE
        'TYN':'TYR', ##  TYN TYR  ADDUCT AT HYDROXY GROUP
        'TYO':'TYR', ##  TYO TYR
        'TYQ':'TYR', ##  TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
        'TYR':'TYR', ##  TYR TYR
        'TYS':'TYR', ##  TYS TYR  INE SULPHONATED TYROSINE
        'TYT':'TYR', ##  TYT TYR
        'TYX':'CYS', ##  TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
        'TYY':'TYR', ##  TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
        'TYZ':'ARG', ##  TYZ ARG  PARA ACETAMIDO BENZOIC ACID
        'UMA':'ALA', ##  UMA ALA
        'VAD':'VAL', ##  VAD VAL  DEAMINOHYDROXYVALINE
        'VAF':'VAL', ##  VAF VAL  METHYLVALINE
        'VAL':'VAL', ##  VAL VAL
        'VDL':'VAL', ##  VDL VAL  (2R,3R)-2,3-DIAMINOBUTANOIC ACID
        'VLL':'VAL', ##  VLL VAL  (2S)-2,3-DIAMINOBUTANOIC ACID
        'VME':'VAL', ##  VME VAL  O- METHYLVALINE
        'X9Q':'ALA', ##  X9Q ALA
        'XX1':'LYS', ##  XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
        'XXY':'ALA', ##  XXY ALA
        'XYG':'ALA', ##  XYG ALA
        'YCM':'CYS', ##  YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
        'YOF':'TYR'
        }



#Creating a main function to be compatible with CWL
def main():
    system         = sys.argv[1]
    step           = sys.argv[2]
    properties_file= sys.argv[3]
    prop           = settings.YamlReader(properties_file, system).get_prop_dic()[step]

    BindingSite(pdb_code        = sys.argv[4],
                pdb_chain       = sys.argv[5],
                output_pdb_path = sys.argv[6],
                properties           = prop    ).launch()


if __name__ == '__main__':
    main()
