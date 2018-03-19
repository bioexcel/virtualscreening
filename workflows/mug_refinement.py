"""
Gromacs full setup from a PDB code
"""


import os
import sys
import time
import shutil
from os.path import join as opj
from command_wrapper import cmd_wrapper
import tools.file_utils as fu
import configuration.settings as settings
import gromacs_wrapper.pdb2gmx as pdb2gmx
import gromacs_wrapper.grompp as grompp
import gromacs_wrapper.solvate as solvate
import gromacs_wrapper.editconf as editconf
import gromacs_wrapper.genion as genion
import gromacs_wrapper.mdrun as mdrun

def main():
    start_time = time.time()
    structure_pdb_path_in=os.path.abspath(sys.argv[1])
    structure_pdb_path_out=os.path.abspath(sys.argv[2])
    yaml_path='/home/user/pymdsetup/workflows/conf/conf_mug_refinement.yaml'
    system='mug_pymdsetup'
    conf = settings.YamlReader(yaml_path, system)
    workflow_path = conf.properties[system]['workflow_path']
    fu.create_dir(os.path.abspath(workflow_path))
    out_log, _ = fu.get_logs(path=workflow_path, console=True)
    paths = conf.get_paths_dic()
    prop = conf.get_prop_dic(global_log=out_log)
    #TODO: Source of problems
    #Change directories always creates problems
    os.chdir(workflow_path)

    out_log.info('\n\n_______MUG REFINEMENT_______\n\n')

    out_log.info('in ----------- Get PDB structure')
    fu.create_dir(prop['step1_mmbpdb']['path'])
    shutil.copy(structure_pdb_path_in, paths['step1_mmbpdb']['output_pdb_path'])

    out_log.info('sed ---------- Replacing atom names')
    sed_path = opj(workflow_path, 'step2_sed')
    fu.create_dir(sed_path)
    sed_pdb_path = opj(sed_path,'sed_replaced.pdb')
    shutil.copy(paths['step1_mmbpdb']['output_pdb_path'], sed_pdb_path)
    cmd = ['sed', '-i', "'s/  O1P  /  OP1  /g'", sed_pdb_path]
    sed1_out_log, sed1_err_log = fu.get_logs(path=sed_path, step='step2_sed1')
    command = cmd_wrapper.CmdWrapper(cmd, sed1_out_log, sed1_err_log)
    command.launch()
    cmd = ['sed', '-i', "'s/  O2P  /  OP2  /g'", sed_pdb_path]
    sed2_out_log, sed2_err_log = fu.get_logs(path=sed_path, step='step2_sed2')
    command = cmd_wrapper.CmdWrapper(cmd, sed2_out_log, sed2_err_log)
    command.launch()

    out_log.info('pdb2gmx ------ Create gromacs topology')
    fu.create_dir(prop['step4_p2g']['path'])
    pdb2gmx.Pdb2gmx(properties=prop['step4_p2g'],
                    input_structure_pdb_path=sed_pdb_path,
                    output_gro_path=paths['step4_p2g']['output_gro_path'],
                    output_top_zip_path=paths['step4_p2g']['output_top_zip_path']).launch()

    out_log.info('editconf ----- Define box dimensions')
    fu.create_dir(prop['step5_ec']['path'])
    editconf.Editconf(properties=prop['step5_ec'], **paths['step5_ec']).launch()

    out_log.info('solvate ------ Fill the box with water molecules')
    fu.create_dir(prop['step6_sol']['path'])
    solvate.Solvate(properties=prop['step6_sol'], **paths['step6_sol']).launch()

    out_log.info('grompp_ions -- Preprocessing: Adding monoatomic ions')
    fu.create_dir(prop['step7_gppions']['path'])
    grompp.Grompp(properties=prop['step7_gppions'], **paths['step7_gppions']).launch()

    out_log.info('genion ------- Running: Adding monoatomic ions')
    fu.create_dir(prop['step8_gio']['path'])
    genion.Genion(properties=prop['step8_gio'], **paths['step8_gio']).launch()

    out_log.info('grompp_min --- Preprocessing: Energy minimization')
    fu.create_dir(prop['step9_gppmin']['path'])
    grompp.Grompp(properties=prop['step9_gppmin'], **paths['step9_gppmin']).launch()

    out_log.info('mdrun_min ---- Running: Energy minimization')
    fu.create_dir(prop['step10_mdmin']['path'])
    mdrun.Mdrun(properties=prop['step10_mdmin'], **paths['step10_mdmin']).launch()

    out_log.info('grompp_nvt --- Preprocessing: nvt constant number of molecules, volume and temp')
    fu.create_dir(prop['step11_gppnvt']['path'])
    grompp.Grompp(properties=prop['step11_gppnvt'], **paths['step11_gppnvt']).launch()

    out_log.info('mdrun_nvt ---- Running: nvt constant number of molecules, volume and temp')
    fu.create_dir(prop['step12_mdnvt']['path'])
    mdrun.Mdrun(properties=prop['step12_mdnvt'], **paths['step12_mdnvt']).launch()

    out_log.info('grompp_npt --- Preprocessing: npt constant number of molecules, pressure and temp')
    fu.create_dir(prop['step13_gppnpt']['path'])
    grompp.Grompp(properties=prop['step13_gppnpt'], **paths['step13_gppnpt']).launch()

    out_log.info('mdrun_npt ---- Running: npt constant number of molecules, pressure and temp')
    fu.create_dir(prop['step14_mdnpt']['path'])
    mdrun.Mdrun(properties=prop['step14_mdnpt'], **paths['step14_mdnpt']).launch()

    out_log.info('grompp_eq ---- Preprocessing: 100ps Molecular dynamics Equilibration')
    fu.create_dir(prop['step15_gppeq']['path'])
    grompp.Grompp(properties=prop['step15_gppeq'], **paths['step15_gppeq']).launch()

    out_log.info('mdrun_eq ----- Running: 100ps Molecular dynamics Equilibration')
    fu.create_dir(prop['step16_mdeq']['path'])
    mdrun.Mdrun(properties=prop['step16_mdeq'], **paths['step16_mdeq']).launch()

    out_log.info('trjconv ------ Extract last snapshot')
    step17_path=opj(workflow_path,'step17_trjconv')
    step17_index=opj(step17_path, 'step17_make_ndx.ndx')
    step17_pdb=opj(step17_path, 'step17_trjconv.pdb')
    fu.create_dir(step17_path)

    #TODO: source of problems
    # should create a wrapper for the make_ndx tool and call it using subprocess
    os.system('printf "! \\"Water_and_ions\\" \nq\n" | '+prop['step16_mdeq']['gmx_path']+' make_ndx -f '+paths['step16_mdeq']['output_gro_path']+' -o '+step17_index+' > '+opj(step17_path, 'make_ndx.out')+' 2> '+ opj(step17_path, 'make_ndx.err'))
    cmd = ['echo', 'Protein_DNA', '|',
           prop['step16_mdeq']['gmx_path'], "trjconv",
           "-s", paths['step15_gppeq']['output_tpr_path'],
           "-f", paths['step16_mdeq']['output_trr_path'],
           "-o", step17_pdb,
           "-n", step17_index,
           "-dump", '1']
    step17_out_log, step17_err_log = fu.get_logs(path=step17_path, step='step17_trjconv')
    command = cmd_wrapper.CmdWrapper(cmd, step17_out_log, step17_err_log)
    command.launch()

    elapsed_time = time.time() - start_time
    removed_list = fu.remove_temp_files(['#', '.top', '.plotscript', '.edr', '.xtc', '.itp', '.top', '.log', '.pdb', '.cpt', '.mdp', '.ndx'])
    shutil.copy(step17_pdb, structure_pdb_path_out)
    out_log.info('')
    out_log.info('Removing unwanted files: ')
    for removed_file in removed_list:
        out_log.info('    X    ' + removed_file)

    out_log.info('\n\nExecution sucessful: ')
    out_log.info('  Output_pdb_path: ' + structure_pdb_path_out)
    out_log.info('  Workflow_path: ' + workflow_path)
    out_log.info('  Config File: ' + yaml_path)
    out_log.info('  System: ' + system)
    out_log.info('  Elapsed time: ' + str(elapsed_time) + 'seconds')

if __name__ == '__main__':
    main()
