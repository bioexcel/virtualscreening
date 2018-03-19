import os
import sys
import time
import shutil
import tools.file_utils as fu
import configuration.settings as settings
import gromacs_wrapper.pdb2gmx as pdb2gmx
import gromacs_wrapper.grompp as grompp
import scwrl_wrapper.scwrl as scwrl
import gromacs_wrapper.solvate as solvate
import gromacs_wrapper.editconf as editconf
import gromacs_wrapper.genion as genion
import gromacs_wrapper.mdrun as mdrun
import mmb_api.pdb as pdb
import mmb_api.uniprot as uniprot
import gromacs_wrapper.rms as rms
import gnuplot_wrapper.gnuplot as gnuplot


def main():
    start_time = time.time()
    yaml_path=sys.argv[1]
    system=sys.argv[2]
    conf = settings.YamlReader(yaml_path, system)
    workflow_path = conf.properties[system]['workflow_path']
    fu.create_dir(os.path.abspath(workflow_path))
    out_log, _ = fu.get_logs(path=workflow_path, console=True)
    paths = conf.get_paths_dic()
    prop = conf.get_prop_dic(global_log=out_log)

    out_log.info('')
    out_log.info('_______MD LAUNCH_______')
    out_log.info('')

    out_log.info('step1: gppeq ---- Preprocessing: Molecular dynamics')
    fu.create_dir(prop['step1_gppeq']['path'])
    grompp.Grompp(properties=prop['step1_gppeq'], **paths['step1_gppeq']).launch()

    out_log.info('step2: mdeq ----- Running: Molecular dynamics')
    fu.create_dir(prop['step2_mdeq']['path'])
    mdrun.Mdrun(properties=prop['step2_mdeq'], **paths['step2_mdeq']).launch()

    #Create setupfiles dir and copy files
    mdfiles_path = os.path.join(workflow_path,'mdfiles')
    fu.create_dir(mdfiles_path)
    shutil.copy(paths['step2_mdeq']['input_tpr_path'], os.path.join(mdfiles_path, 'md.tpr'))
    shutil.copy(paths['step2_mdeq']['output_cpt_path'], os.path.join(mdfiles_path, 'md.cpt'))
    shutil.copy(paths['step2_mdeq']['output_gro_path'], os.path.join(mdfiles_path, 'md.gro'))
    shutil.copy('_step2_mdeq_mdeq.xtc', os.path.join(mdfiles_path, 'md.xtc'))

    removed_list = fu.remove_temp_files(['#', '.top', '.plotscript', '.edr', '.xtc', '.itp', '.top', '.log', '.pdb', '.cpt', '.mdp', '.trr'])
    out_log.info('')
    out_log.info('Removing unwanted files: ')
    for removed_file in removed_list:
        out_log.info('    X    ' + removed_file)

    elapsed_time = time.time() - start_time
    out_log.info('')
    out_log.info('')
    out_log.info('Execution sucessful: ')
    out_log.info('  Workflow_path: '+workflow_path)
    out_log.info('  Config File: '+yaml_path)
    out_log.info('  System: '+system)
    if len(sys.argv) >= 4:
        out_log.info('  Nodes: '+sys.argv[3])
    out_log.info('')
    out_log.info('Elapsed time: '+str(elapsed_time)+' seconds')
    out_log.info('')

if __name__ == '__main__':
    main()
