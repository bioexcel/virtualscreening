import os
import sys
import time
from os.path import join as opj
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
from pycompss.api.parameter import *
from pycompss.api.task import task
from pycompss.api.constraint import constraint

def main():
    from pycompss.api.api import compss_open
    start_time = time.time()
    yaml_path=sys.argv[1]
    system=sys.argv[2]
    conf = settings.YamlReader(yaml_path, system)
    workflow_path = conf.properties[system]['workflow_path']
    fu.create_dir(os.path.abspath(workflow_path))
    out_log, _ = fu.get_logs(path=workflow_path, console=True)
    paths_glob = conf.get_paths_dic()
    prop_glob = conf.get_prop_dic()

    out_log.info('')
    out_log.info('_______GROMACS FULL WORKFLOW_______')
    out_log.info('')

    out_log.info( 'step1:  mmbpdb -- Get PDB')
    structure = conf.properties[system].get('initial_structure_pdb_path', None)
    if structure is None or not os.path.isfile(structure):
        out_log.info( '     Selected PDB code: ' + prop_glob['step1_mmbpdb']['pdb_code'])
        fu.create_dir(prop_glob['step1_mmbpdb']['path'])
        pdb.MmbPdb().get_pdb(prop_glob['step1_mmbpdb']['pdb_code'], paths_glob['step1_mmbpdb']['output_pdb_path'])
        structure = paths_glob['step1_mmbpdb']['output_pdb_path']

    out_log.info( 'step2:  mmbuniprot -- Get mutations')
    mutations = conf.properties.get('input_mapped_mutations_list', None)
    if mutations is None or len(mutations) < 7:
        mmbuniprot = uniprot.MmbVariants(prop_glob['step1_mmbpdb']['pdb_code'])
        mutations = mmbuniprot.get_pdb_variants()
        if mutations is None or len(mutations) == 0: return
    else:
        mutations = [m.strip() for m in conf.properties.get('input_mapped_mutations_list').split(',')]

    mutations_limit = min(len(mutations), int(prop_glob.get('mutations_limit', len(mutations))))
    out_log.info('')
    out_log.info('Number of mutations to be modelled: ' + str(mutations_limit))

    rms_list = []
    mutations_counter = 0
    for mut in mutations:
        if mutations_counter == mutations_limit: break
        mutations_counter += 1
        paths = conf.get_paths_dic(mut)
        prop = conf.get_prop_dic(mut)


        out_log.info('')
        out_log.info('-------------------------')
        out_log.info(str(mutations_counter) + '/' + str(mutations_limit) + ' ' + mut)
        out_log.info('-------------------------')
        out_log.info('')

        out_log.info('step3:  scw ------ Model mutation')
        fu.create_dir(prop['step3_scw']['path'])
        paths['step3_scw']['input_pdb_path']=structure
        scwrl_pc(properties=prop['step3_scw'], **paths['step3_scw'])

        out_log.info('step4:  p2g ------ Create gromacs topology')
        fu.create_dir(prop['step4_p2g']['path'])
        pdb2gmx_pc(properties=prop['step4_p2g'], **paths['step4_p2g'])

        out_log.info('step5:  ec ------- Define box dimensions')
        fu.create_dir(prop['step5_ec']['path'])
        editconf_pc(properties=prop['step5_ec'], **paths['step5_ec'])

        out_log.info('step6:  sol ------ Fill the box with water molecules')
        fu.create_dir(prop['step6_sol']['path'])
        solvate_pc(properties=prop['step6_sol'], **paths['step6_sol'])

        out_log.info('step7:  gppions -- Preprocessing: Adding monoatomic ions')
        fu.create_dir(prop['step7_gppions']['path'])
        grompp_pc(properties=prop['step7_gppions'], **paths['step7_gppions'])

        out_log.info('step8:  gio ------ Running: Adding monoatomic ions')
        fu.create_dir(prop['step8_gio']['path'])
        genion_pc(properties=prop['step8_gio'], **paths['step8_gio'])

        out_log.info('step9:  gppmin --- Preprocessing: Energy minimization')
        fu.create_dir(prop['step9_gppmin']['path'])
        grompp_pc(properties=prop['step9_gppmin'], **paths['step9_gppmin'])

        out_log.info('step10: mdmin ---- Running: Energy minimization')
        fu.create_dir(prop['step10_mdmin']['path'])
        mdrun_pc(properties=prop['step10_mdmin'], **paths['step10_mdmin'])

        out_log.info('step11: gppnvt --- Preprocessing: nvt constant number of molecules, volume and temp')
        fu.create_dir(prop['step11_gppnvt']['path'])
        grompp_pc(properties=prop['step11_gppnvt'], **paths['step11_gppnvt'])

        out_log.info('step12: mdnvt ---- Running: nvt constant number of molecules, volume and temp')
        fu.create_dir(prop['step12_mdnvt']['path'])
        mdrun_pc_cpt(properties=prop['step12_mdnvt'], **paths['step12_mdnvt'])

        out_log.info('step13: gppnpt --- Preprocessing: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step13_gppnpt']['path'])
        grompp_pc_cpt(properties=prop['step13_gppnpt'], **paths['step13_gppnpt'])

        out_log.info('step14: mdnpt ---- Running: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step14_mdnpt']['path'])
        mdrun_pc_cpt(properties=prop['step14_mdnpt'], **paths['step14_mdnpt'])

        out_log.info('step15: gppeq ---- Preprocessing: 1ns Molecular dynamics Equilibration')
        fu.create_dir(prop['step15_gppeq']['path'])
        grompp_pc_cpt(properties=prop['step15_gppeq'], **paths['step15_gppeq'])

        out_log.info('step16: mdeq ----- Running: 1ns Molecular dynamics Equilibration')
        fu.create_dir(prop['step16_mdeq']['path'])
        mdrun_pc(properties=prop['step16_mdeq'], **paths['step16_mdeq'])

        out_log.info('step17: rmsd ----- Computing RMSD')
        fu.create_dir(prop['step17_rmsd']['path'])
        rms_list.append(rms_pc(properties=prop['step17_rmsd'], **paths['step17_rmsd']))

    xvg_dict = reduce(merge_dictionaries, rms_list)
    out_log.info('step18: gnuplot ----- Creating RMSD plot')
    fu.create_dir(prop_glob['step18_gnuplot']['path'])
    output_png_path = paths_glob['step18_gnuplot']['output_png_path']
    properties = prop_glob['step18_gnuplot']
    gnuplot_pc(xvg_dict, output_png_path, properties)
    png = compss_open(output_png_path)
    elapsed_time = time.time() - start_time

    with open(opj(workflow_path, 'time.txt'), 'a') as time_file:
        time_file.write('Elapsed time: ')
        time_file.write(str(elapsed_time))
        time_file.write('\n')
        time_file.write('Config File: ')
        time_file.write(sys.argv[1])
        time_file.write('\n')
        time_file.write('Sytem: ')
        time_file.write(sys.argv[2])
        time_file.write('\n')
        if len(sys.argv) >= 4:
            time_file.write('Nodes: ')
            time_file.write(sys.argv[3])
            time_file.write('\n')


############################## PyCOMPSs functions #############################
#Minotauro
#computing_units = "0"
#MareNostrum4
computing_units = "48"


@task(returns=dict)
def merge_dictionaries(a, b):
    return dict(a, **b)

@task(input_pdb_path=FILE_IN, output_pdb_path=FILE_OUT)
def scwrl_pc(input_pdb_path, output_pdb_path, properties, **kwargs):
    scwrl.Scwrl4(input_pdb_path=input_pdb_path, output_pdb_path=output_pdb_path, properties=properties, **kwargs).launch()

@task(input_structure_pdb_path=FILE_IN, output_gro_path=FILE_OUT, output_top_zip_path=FILE_OUT)
def pdb2gmx_pc(input_structure_pdb_path, output_gro_path, output_top_zip_path, properties, **kwargs):
    pdb2gmx.Pdb2gmx(input_structure_pdb_path=input_structure_pdb_path, output_gro_path=output_gro_path, output_top_zip_path=output_top_zip_path, properties=properties, **kwargs).launch()

@task(input_gro_path=FILE_IN, output_gro_path=FILE_OUT)
def editconf_pc(input_gro_path, output_gro_path, properties, **kwargs):
    editconf.Editconf(input_gro_path=input_gro_path, output_gro_path=output_gro_path, properties=properties, **kwargs).launch()

@task(input_solute_gro_path=FILE_IN, output_gro_path=FILE_OUT, input_top_zip_path=FILE_IN, output_top_zip_path=FILE_OUT)
def solvate_pc(input_solute_gro_path, output_gro_path, input_top_zip_path, output_top_zip_path, properties, **kwargs):
    solvate.Solvate(input_solute_gro_path=input_solute_gro_path, output_gro_path=output_gro_path, input_top_zip_path=input_top_zip_path, output_top_zip_path=output_top_zip_path, properties=properties, **kwargs).launch()

@task(input_gro_path=FILE_IN, input_top_zip_path=FILE_IN, input_mdp_path=FILE_IN, output_tpr_path=FILE_OUT,  input_cpt_path=FILE_IN)
def grompp_pc_cpt(input_gro_path, input_top_zip_path, input_mdp_path, output_tpr_path, properties, input_cpt_path, **kwargs):
    grompp.Grompp(input_gro_path=input_gro_path, input_top_zip_path=input_top_zip_path, input_mdp_path=input_mdp_path, output_tpr_path=output_tpr_path, properties=properties, input_cpt_path=input_cpt_path, **kwargs).launch()

@task(input_gro_path=FILE_IN, input_top_zip_path=FILE_IN, input_mdp_path=FILE_IN, output_tpr_path=FILE_OUT)
def grompp_pc(input_gro_path, input_top_zip_path, input_mdp_path, output_tpr_path, properties, **kwargs):
    grompp.Grompp(input_gro_path=input_gro_path, input_top_zip_path=input_top_zip_path, input_mdp_path=input_mdp_path, output_tpr_path=output_tpr_path, properties=properties, **kwargs).launch()

@task(input_tpr_path=FILE_IN, output_gro_path=FILE_OUT, input_top_zip_path=FILE_IN, output_top_zip_path=FILE_OUT)
def genion_pc(input_tpr_path, output_gro_path, input_top_zip_path, output_top_zip_path, properties, **kwargs):
    genion.Genion(input_tpr_path=input_tpr_path, output_gro_path=output_gro_path, input_top_zip_path=input_top_zip_path, output_top_zip_path=output_top_zip_path, properties=properties, **kwargs).launch()

@constraint(ComputingUnits=computing_units)
@task(input_tpr_path=FILE_IN, output_trr_path=FILE_OUT, output_gro_path=FILE_OUT, output_cpt_path=FILE_OUT)
def mdrun_pc_cpt(input_tpr_path, output_trr_path, output_gro_path, properties, output_cpt_path, **kwargs):
    mdrun.Mdrun(input_tpr_path=input_tpr_path, output_trr_path=output_trr_path, properties=properties, output_gro_path=output_gro_path, output_cpt_path=output_cpt_path, **kwargs).launch()

@constraint(ComputingUnits=computing_units)
@task(input_tpr_path=FILE_IN, output_trr_path=FILE_OUT, output_gro_path=FILE_OUT)
def mdrun_pc(input_tpr_path, output_trr_path, output_gro_path, properties, **kwargs):
    mdrun.Mdrun(input_tpr_path=input_tpr_path, output_trr_path=output_trr_path, properties=properties, output_gro_path=output_gro_path, **kwargs).launch()

@task(input_gro_path=FILE_IN, input_trr_path=FILE_IN, output_xvg_path=FILE_OUT, returns=dict)
def rms_pc(input_gro_path, input_trr_path, output_xvg_path, properties, **kwargs):
    return rms.Rms(input_gro_path=input_gro_path, input_trr_path=input_trr_path, output_xvg_path=output_xvg_path, properties=properties, **kwargs).launch()

@task(output_png_path=FILE_OUT)
def gnuplot_pc(input_xvg_path_dict, output_png_path, properties):
    gnuplot.Gnuplot(input_xvg_path_dict=input_xvg_path_dict, output_png_path=output_png_path, properties=properties).launch()

if __name__ == '__main__':
    main()
