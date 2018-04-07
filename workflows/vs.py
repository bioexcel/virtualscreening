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
import workflows.egfr_md as egfr_md






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

    out_log.info('md_clustering workflow')
    fu.create_dir(prop['md_clustering']['path'])
    out_log.debug('\nPaths:\n'+str(paths['md_clustering'])+'\nProperties:\n'+str(prop['md_clustering'])+'\n')
    receptors_pdb = egfr_md.Workflow(yaml_path=prop['md_clustering']['yaml_path'], system=prop['md_clustering']['system'], workflow=paths['md_clustering']['workflow_path']).launch()

    out_log.info('Receptors pdb file: '+receptors_pdb)

    # out_log.info('SDF2PDB: Convert SDF2PDB')
    # fu.create_dir(prop['SDF2PDB']['path'])
    # paths['SDF2PDB']['input_sdf_path']=conf.properties[system].get('initial_structure_pdb_path', None)
    # out_log.debug('\nPaths:\n'+str(paths['prepare_ligand'])+'\nProperties:\n'+str(prop['prepare_ligand'])+'\n')
    # scwrl.Scwrl4(properties=prop['prepare_ligand'], **paths['prepare_ligand']).launch())
    #
    #
    # out_log.info('Prepare Ligands:')
    # fu.create_dir(prop['prepare_ligand']['path'])
    # paths['prepare_ligand']['input_ligand_pdb_path']=conf.properties[system].get('initial_structure_pdb_path', None)
    # out_log.debug('\nPaths:\n'+str(paths['prepare_ligand'])+'\nProperties:\n'+str(prop['prepare_ligand'])+'\n')
    # scwrl.Scwrl4(properties=prop['prepare_ligand'], **paths['prepare_ligand']).launch()



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
