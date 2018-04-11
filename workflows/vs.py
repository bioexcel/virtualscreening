import os
import sys
import time
import tools.file_utils as fu
import configuration.settings as settings
import gromacs_wrapper.pdb2gmx as pdb2gmx
import gromacs_wrapper.grompp as grompp
import scwrl_wrapper.scwrl as scwrl
import gromacs_wrapper.solvate as solvate
import gromacs_wrapper.editconf as editconf
import gromacs_wrapper.genion as genion
import gromacs_wrapper.mdrun as mdrun
import gromacs_wrapper.make_ndx as make_ndx
import gromacs_wrapper.genrestr as genrestr
import mmb_api.pdb as pdb
import mmb_api.uniprot as uniprot
import gromacs_wrapper.rms as rms
import gnuplot_wrapper.gnuplot as gnuplot
import gromacs_extra.ndx2resttop as ndx2resttop
import vina_wrapper.sdf2pdb as sdf2pdb
import vina_wrapper.prepare_ligand as prepare_ligand
import vina_wrapper.prepare_receptor as prepare_receptor
import vina_wrapper.vina as vina
import bindingsite.bindingsite as bindingsite
import bindingsite.box as box
import md_cluster
import Bio.PDB
import dude_api.dude as dude
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


def main():
    start_time = time.time()
    yaml_path=sys.argv[1]
    system=sys.argv[2]

    conf = settings.YamlReader(yaml_path, system)
    workflow_path = conf.properties[system]['workflow_path']
    fu.create_dir(os.path.abspath(workflow_path))
    out_log, _ = fu.get_logs(path=workflow_path, console=True, level='DEBUG')
    paths = conf.get_paths_dic()
    prop = conf.get_prop_dic(global_log=out_log)

    out_log.info('')
    out_log.info('_______GROMACS FULL WORKFLOW_______')
    out_log.info('')
    out_log.info("Command Executed:")
    out_log.info(" ".join(sys.argv))
    out_log.info('Workflow_path: '+workflow_path)
    out_log.info('Config File: '+yaml_path)
    out_log.info('System: '+system)
    out_log.info('')

    #Download initial structure
    out_log.info( 'step1:  mmbpdb -- Get PDB')
    structure = conf.properties[system].get('initial_structure_pdb_path', None)
    if structure is None or not os.path.isfile(structure):
        out_log.info( '     Selected PDB code: ' + prop['step1_mmbpdb']['pdb_code'])
        fu.create_dir(prop['step1_mmbpdb']['path'])
        pdb.MmbPdb().get_pdb(prop['step1_mmbpdb']['pdb_code'], paths['step1_mmbpdb']['output_pdb_path'], "filter=:" + prop['step1_mmbpdb']['pdb_chain'])
        structure = paths['step1_mmbpdb']['output_pdb_path']

    decoys_sdf_path = conf.properties[system].get('decoys_sdf_path', None)
    if not decoys_sdf_path:
        #Dowload decoys
        out_log.info('Downloading decoys from DUDe')
        fu.create_dir(prop['dude']['path'])
        out_log.debug('\nPaths:\n'+str(paths['dude'])+'\nProperties:\n'+str(prop['dude'])+'\n')
        decoys_sdf_path = dude.Dude(properties=prop['dude'], **paths['dude']).launch()
    #SDF to PDB list
    out_log.info('Converting sdf to pdb')
    fu.create_dir(prop['sdf2pdb']['path'])
    out_log.debug('\nPaths:\n'+str(paths['sdf2pdb'])+'\nProperties:\n'+str(prop['sdf2pdb'])+'\n')
    decoys_pdb_list = sdf2pdb.SDF2PDB(input_sdf_path=decoys_sdf_path, properties=prop['sdf2pdb']).launch()


    #Get Actives
    actives_dir_path = conf.properties[system].get('actives_dir_path', None)
    actives_pdb_list = []
    if not actives_dir_path:
        pass
    for dirpath,_,filenames in os.walk(actives_dir_path):
       for f in filenames:
           actives_pdb_list.append(os.path.abspath(os.path.join(dirpath, f)))

    #Merge actives and decoys
    small_molecules_list = decoys_pdb_list + actives_pdb_list

    #Get bindingsite
    out_log.info('bindingsite')
    fu.create_dir(prop['bindingsite']['path'])
    out_log.debug('\nPaths:\n'+str(paths['bindingsite'])+'\nProperties:\n'+str(prop['bindingsite'])+'\n')
    bindingsite.BindingSite(properties=prop['bindingsite'], **paths['bindingsite']).launch()

    # Get Receptor pdb cluster
    receptors_pdb_path = conf.properties[system].get('receptors_pdb_path', None)
    if not receptors_pdb_path:
        out_log.info('md_clustering workflow')
        fu.create_dir(prop['md_clustering']['path'])
        out_log.debug('\nPaths:\n'+str(paths['md_clustering'])+'\nProperties:\n'+str(prop['md_clustering'])+'\n')
        receptors_pdb_path = md_cluster.MDCluster(yaml_path=prop['md_clustering']['yaml_path'], system=prop['md_clustering']['system'], workflow_path=paths['md_clustering']['workflow_path'], structure_pdb_path=structure).launch()

    # Get BindingSite BOX
    out_log.info('box')
    fu.create_dir(prop['box']['path'])
    out_log.debug('\nPaths:\n'+str(paths['box'])+'\nProperties:\n'+str(prop['box'])+'\n')
    box.Box(properties=prop['box'], input_pdb_path=receptors_pdb_path, **paths['box']).launch()

    # Split receptor PDB into multiple PDBs
    parser=Bio.PDB.PDBParser()
    st = parser.get_structure('st', receptors_pdb)
    for i in range(len(st)):
        structure_tmp = st.copy()
        for j in range(len(st)):
            if j != i:
                structure_tmp.detach_child(j)
        io = Bio.PDB.PDBIO()
        io.set_structure(structure_tmp)
        receptor_name='receptor_'+str(i)+'.pdb'
        io.save(receptor_name)
        receptors_list.append(receptor_name)

    max_receptors = conf.properties[system].get('max_receptors', len(receptors_list))
    max_ligands = conf.properties[system].get('max_ligands', len(small_molecules_list))
    for receptor_index, receptor_path in enumerate(receptor_list):
        if receptor_index >= max_receptors: break
        #Prepare receptor
        out_log.info('Prepare receptor')
        fu.create_dir(prop['prepare_receptor']['path'])
        out_log.debug('\nPaths:\n'+str(paths['prepare_receptor'])+'\nProperties:\n'+str(prop['prepare_receptor'])+'\n')
        prepare_receptor.VinaPrepareReceptor(properties=prop['prepare_receptor'], input_receptor_pdb_path=receptor_path, **paths['prepare_receptor']).launch()

        for ligand_index, ligand_path in enumerate(small_molecules_list):
            if ligand_index >= max_ligands: break
            #Prepare Ligand
            out_log.info('Prepare ligand')
            fu.create_dir(prop['prepare_ligand']['path'])
            out_log.debug('\nPaths:\n'+str(paths['prepare_ligand'])+'\nProperties:\n'+str(prop['prepare_ligand'])+'\n')
            prepare_ligand.VinaPrepareLigand(properties=prop['prepare_ligand'], input_ligand_pdb_path=ligand_path, **paths['prepare_ligand']).launch()
            #Launch vina docking
            out_log.info('vina')
            fu.create_dir(prop['vina']['path'])
            out_log.debug('\nPaths:\n'+str(paths['vina'])+'\nProperties:\n'+str(prop['vina'])+'\n')
            vina.Vina(properties=prop['vina'], **paths['vina']).launch()


    elapsed_time = time.time() - start_time
    out_log.info('')
    out_log.info('')
    out_log.info('Execution sucessful: ')
    out_log.info('  Workflow_path: '+workflow_path)
    out_log.info('  Config File: '+yaml_path)
    out_log.info('  System: '+system)
    out_log.info('')
    out_log.info('Elapsed time: '+str(elapsed_time)+' seconds')
    out_log.info('')

if __name__ == '__main__':
    main()
