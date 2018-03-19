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
    out_log.info('_______MD SETUP FULL WORKFLOW_______')
    out_log.info('')

    out_log.info('step1:  mmbpdb --- Get PDB')
    structure = conf.properties[system].get('initial_structure_pdb_path', None)
    if structure is None or not os.path.isfile(structure):
        out_log.info( '                   Selected PDB code: ' + prop['step1_mmbpdb']['pdb_code'])
        fu.create_dir(prop['step1_mmbpdb']['path'])
        pdb.MmbPdb().get_pdb(prop['step1_mmbpdb']['pdb_code'], paths['step1_mmbpdb']['output_pdb_path'])
        structure = paths['step1_mmbpdb']['output_pdb_path']

    out_log.info('step2:  p2g ------ Create gromacs topology')
    fu.create_dir(prop['step2_p2g']['path'])
    pdb2gmx.Pdb2gmx(properties=prop['step2_p2g'], **paths['step2_p2g']).launch()

    out_log.info('step3:  ec ------- Define box dimensions')
    fu.create_dir(prop['step3_ec']['path'])
    editconf.Editconf(properties=prop['step3_ec'], **paths['step3_ec']).launch()

    out_log.info('step4:  sol ------ Fill the box with water molecules')
    fu.create_dir(prop['step4_sol']['path'])
    solvate.Solvate(properties=prop['step4_sol'], **paths['step4_sol']).launch()

    out_log.info('step5:  gppions -- Preprocessing: Adding monoatomic ions')
    fu.create_dir(prop['step5_gppions']['path'])
    grompp.Grompp(properties=prop['step5_gppions'], **paths['step5_gppions']).launch()

    out_log.info('step6:  gio ------ Running: Adding monoatomic ions')
    fu.create_dir(prop['step6_gio']['path'])
    genion.Genion(properties=prop['step6_gio'], **paths['step6_gio']).launch()

    out_log.info('step7:  gppmin --- Preprocessing: Energy minimization')
    fu.create_dir(prop['step7_gppmin']['path'])
    grompp.Grompp(properties=prop['step7_gppmin'], **paths['step7_gppmin']).launch()

    out_log.info('step8:  mdmin ---- Running: Energy minimization')
    fu.create_dir(prop['step8_mdmin']['path'])
    mdrun.Mdrun(properties=prop['step8_mdmin'], **paths['step8_mdmin']).launch()

    out_log.info('step9:  gppnvt --- Preprocessing: nvt constant number of molecules, volume and temp')
    fu.create_dir(prop['step9_gppnvt']['path'])
    grompp.Grompp(properties=prop['step9_gppnvt'], **paths['step9_gppnvt']).launch()

    out_log.info('step10: mdnvt ---- Running: nvt constant number of molecules, volume and temp')
    fu.create_dir(prop['step10_mdnvt']['path'])
    mdrun.Mdrun(properties=prop['step10_mdnvt'], **paths['step10_mdnvt']).launch()

    out_log.info('step11: gppnpt --- Preprocessing: npt constant number of molecules, pressure and temp')
    fu.create_dir(prop['step11_gppnpt']['path'])
    grompp.Grompp(properties=prop['step11_gppnpt'], **paths['step11_gppnpt']).launch()

    out_log.info('step12: mdnpt ---- Running: npt constant number of molecules, pressure and temp')
    fu.create_dir(prop['step12_mdnpt']['path'])
    mdrun.Mdrun(properties=prop['step12_mdnpt'], **paths['step12_mdnpt']).launch()

    out_log.info('step13: gppeq ---- Preprocessing: Molecular dynamics Equilibration')
    fu.create_dir(prop['step13_gppeq']['path'])
    grompp.Grompp(properties=prop['step13_gppeq'], **paths['step13_gppeq']).launch()

    out_log.info('step14: mdeq ----- Running: Molecular dynamics Equilibration')
    fu.create_dir(prop['step14_mdeq']['path'])
    mdrun.Mdrun(properties=prop['step14_mdeq'], **paths['step14_mdeq']).launch()

    #Create setupfiles dir and copy files
    setupfiles_path = os.path.join(workflow_path,'setupfiles')
    fu.create_dir(setupfiles_path)
    shutil.copy(paths['step14_mdeq']['output_gro_path'], os.path.join(setupfiles_path, 'md_setup.gro'))
    shutil.copy(paths['step14_mdeq']['output_cpt_path'], os.path.join(setupfiles_path, 'md_setup.cpt'))
    shutil.copy(paths['step6_gio']['output_top_zip_path'], os.path.join(setupfiles_path, 'md_setup.zip'))

    removed_list = fu.remove_temp_files(['#', '.top', '.plotscript', '.edr', '.xtc', '.itp', '.top', '.log', '.pdb', '.cpt', '.mdp'])
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
