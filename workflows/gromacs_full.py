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
import mmb_api.pdb as pdb
import mmb_api.uniprot as uniprot
import gromacs_wrapper.rms as rms
import gnuplot_wrapper.gnuplot as gnuplot


def main():
    start_time = time.time()
    yaml_path=sys.argv[1]
    system=sys.argv[2]
    n_mutations=sys.argv[3]

    conf = settings.YamlReader(yaml_path, system)
    workflow_path = conf.properties[system]['workflow_path']
    fu.create_dir(os.path.abspath(workflow_path))
    out_log, _ = fu.get_logs(path=workflow_path, console=True)
    paths_glob = conf.get_paths_dic()
    prop_glob = conf.get_prop_dic()

    out_log.info('')
    out_log.info('_______GROMACS FULL WORKFLOW_______')
    out_log.info('')

    out_log.info("Command Executed:")
    out_log.info(" ".join(sys.argv))
    out_log.info('Workflow_path: '+workflow_path)
    out_log.info('Config File: '+yaml_path)
    out_log.info('System: '+system)
    out_log.info('Mutations limit: '+n_mutations)
    if len(sys.argv) >= 5:
        out_log.info('Nodes: '+sys.argv[4])

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

    mutations_limit = min(len(mutations), int(n_mutations))
    out_log.info('')
    out_log.info('Number of mutations to be modelled: ' + str(mutations_limit))

    rms_list = []
    mutations_counter = 0
    for mut in mutations:
        if mutations_counter == mutations_limit: break
        mutations_counter += 1
        paths = conf.get_paths_dic(mut)
        prop = conf.get_prop_dic(mut, global_log=out_log)

        out_log.info('')
        out_log.info('-------------------------')
        out_log.info(str(mutations_counter) + '/' + str(mutations_limit) + ' ' + mut)
        out_log.info('-------------------------')
        out_log.info('')

        out_log.info('step3:  scw ------ Model mutation')
        fu.create_dir(prop['step3_scw']['path'])
        paths['step3_scw']['input_pdb_path']=structure
        scwrl.Scwrl4(properties=prop['step3_scw'], **paths['step3_scw']).launch()

        out_log.info('step4:  p2g ------ Create gromacs topology')
        fu.create_dir(prop['step4_p2g']['path'])
        pdb2gmx.Pdb2gmx(properties=prop['step4_p2g'], **paths['step4_p2g']).launch()

        out_log.info('step5:  ec ------- Define box dimensions')
        fu.create_dir(prop['step5_ec']['path'])
        editconf.Editconf(properties=prop['step5_ec'], **paths['step5_ec']).launch()

        out_log.info('step6:  sol ------ Fill the box with water molecules')
        fu.create_dir(prop['step6_sol']['path'])
        solvate.Solvate(properties=prop['step6_sol'], **paths['step6_sol']).launch()

        out_log.info('step7:  gppions -- Preprocessing: Adding monoatomic ions')
        fu.create_dir(prop['step7_gppions']['path'])
        grompp.Grompp(properties=prop['step7_gppions'], **paths['step7_gppions']).launch()

        out_log.info('step8:  gio ------ Running: Adding monoatomic ions')
        fu.create_dir(prop['step8_gio']['path'])
        genion.Genion(properties=prop['step8_gio'], **paths['step8_gio']).launch()

        out_log.info('step9:  gppmin --- Preprocessing: Energy minimization')
        fu.create_dir(prop['step9_gppmin']['path'])
        grompp.Grompp(properties=prop['step9_gppmin'], **paths['step9_gppmin']).launch()

        out_log.info('step10: mdmin ---- Running: Energy minimization')
        fu.create_dir(prop['step10_mdmin']['path'])
        mdrun.Mdrun(properties=prop['step10_mdmin'], **paths['step10_mdmin']).launch()

        out_log.info('step11: gppnvt --- Preprocessing: nvt constant number of molecules, volume and temp')
        fu.create_dir(prop['step11_gppnvt']['path'])
        grompp.Grompp(properties=prop['step11_gppnvt'], **paths['step11_gppnvt']).launch()

        out_log.info('step12: mdnvt ---- Running: nvt constant number of molecules, volume and temp')
        fu.create_dir(prop['step12_mdnvt']['path'])
        mdrun.Mdrun(properties=prop['step12_mdnvt'], **paths['step12_mdnvt']).launch()

        out_log.info('step13: gppnpt --- Preprocessing: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step13_gppnpt']['path'])
        grompp.Grompp(properties=prop['step13_gppnpt'], **paths['step13_gppnpt']).launch()

        out_log.info('step14: mdnpt ---- Running: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step14_mdnpt']['path'])
        mdrun.Mdrun(properties=prop['step14_mdnpt'], **paths['step14_mdnpt']).launch()

        out_log.info('step15: gppeq ---- Preprocessing: 1ns Molecular dynamics Equilibration')
        fu.create_dir(prop['step15_gppeq']['path'])
        grompp.Grompp(properties=prop['step15_gppeq'], **paths['step15_gppeq']).launch()

        out_log.info('step16: mdeq ----- Running: Free Molecular dynamics Equilibration')
        fu.create_dir(prop['step16_mdeq']['path'])
        mdrun.Mdrun(properties=prop['step16_mdeq'], **paths['step16_mdeq']).launch()

        out_log.info('step17: rmsd ----- Computing RMSD')
        fu.create_dir(prop['step17_rmsd']['path'])
        rms_list.append(rms.Rms(properties=prop['step17_rmsd'], **paths['step17_rmsd']).launch())
        removed_list = fu.remove_temp_files(['#', '.top', '.plotscript', '.edr', '.xtc', '.itp', '.top', '.log', '.pdb', '.cpt', '.mdp', '.xvg', '.grp', '.seq'])

    xvg_dict=reduce(lambda a, b: dict(a, **b), rms_list)
    out_log.info('step18: gnuplot ----- Creating RMSD plot')
    fu.create_dir(prop_glob['step18_gnuplot']['path'])
    gnuplot.Gnuplot(input_xvg_path_dict=xvg_dict, properties=prop_glob['step18_gnuplot'], **paths_glob['step18_gnuplot']).launch()
    elapsed_time = time.time() - start_time


    out_log.info('')
    out_log.info('')
    out_log.info('Execution sucessful: ')
    out_log.info('  Workflow_path: '+workflow_path)
    out_log.info('  Config File: '+yaml_path)
    out_log.info('  System: '+system)
    if len(sys.argv) >= 5:
        out_log.info('  Nodes: '+sys.argv[4])
    out_log.info('')
    out_log.info('Elapsed time: '+str(elapsed_time)+' seconds')
    out_log.info('')

if __name__ == '__main__':
    main()
