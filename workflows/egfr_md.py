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


class Workflow(object):
    def __init__(self, yaml_path, system, workflow_path):
        self.yaml_path=yaml_path
        self.system=system
        self.workflow_path=workflow_path

    def launch():
        yaml_path=self.yaml_path
        system=self.system
        workflow_path=self.workflow_path
        start_time = time.time()
        conf = settings.YamlReader(yaml_path, system)
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

        out_log.info( 'step1:  mmbpdb ------ Get PDB')
        structure = conf.properties[system].get('initial_structure_pdb_path', None)
        out_log.info( 22*' '+'Selected PDB structure: ' + structure)

        out_log.info('step2 and 3:  scw --- Check and Repair PDB')
        fu.create_dir(prop['step3_scw']['path'])
        paths['step3_scw']['input_pdb_path']=structure
        out_log.debug('\nPaths:\n'+str(paths['step3_scw'])+'\nProperties:\n'+str(prop['step3_scw'])+'\n')
        scwrl.Scwrl4(properties=prop['step3_scw'], **paths['step3_scw']).launch()

        out_log.info('step4:  p2g --------- Create gromacs topology')
        fu.create_dir(prop['step4_p2g']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step4_p2g'])+'\nProperties:\n'+str(prop['step4_p2g'])+'\n')
        pdb2gmx.Pdb2gmx(properties=prop['step4_p2g'], **paths['step4_p2g']).launch()

        out_log.info('step5:  ec ---------- Define box dimensions')
        fu.create_dir(prop['step5_ec']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step5_ec'])+'\nProperties:\n'+str(prop['step5_ec'])+'\n')
        editconf.Editconf(properties=prop['step5_ec'], **paths['step5_ec']).launch()

        out_log.info('step6:  sol --------- Fill the box with water molecules')
        fu.create_dir(prop['step6_sol']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step6_sol'])+'\nProperties:\n'+str(prop['step6_sol'])+'\n')
        solvate.Solvate(properties=prop['step6_sol'], **paths['step6_sol']).launch()

        out_log.info('step7:  gppions ----- Preprocessing: Adding monoatomic ions')
        fu.create_dir(prop['step7_gppions']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step7_gppions'])+'\nProperties:\n'+str(prop['step7_gppions'])+'\n')
        grompp.Grompp(properties=prop['step7_gppions'], **paths['step7_gppions']).launch()

        out_log.info('step8:  gio --------- Running: Adding monoatomic ions')
        fu.create_dir(prop['step8_gio']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step8_gio'])+'\nProperties:\n'+str(prop['step8_gio'])+'\n')
        genion.Genion(properties=prop['step8_gio'], **paths['step8_gio']).launch()

        out_log.info('Step9: gppndx ------- Preprocessing index creation')
        fu.create_dir(prop['step9_gppndx']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step9_gppndx'])+'\nProperties:\n'+str(prop['step9_gppndx'])+'\n')
        grompp.Grompp(properties=prop['step9_gppndx'], **paths['step9_gppndx']).launch()

        out_log.info('Step10: make_ndx ---- Create restrain index')
        fu.create_dir(prop['step10_make_ndx']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step10_make_ndx'])+'\nProperties:\n'+str(prop['step10_make_ndx'])+'\n')
        make_ndx.MakeNdx(properties=prop['step10_make_ndx'], **paths['step10_make_ndx']).launch()

        out_log.info('Step11: genrest - Create restrain topology')
        fu.create_dir(prop['step11_genrestr']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step11_genrestr'])+'\nProperties:\n'+str(prop['step11_genrestr'])+'\n')
        genrestr.Genrestr(properties=prop['step11_genrestr'], **paths['step11_genrestr']).launch()

        out_log.info('step15: gppmin ------ Preprocessing: minimization')
        fu.create_dir(prop['step15_gppmin']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step15_gppmin'])+'\nProperties:\n'+str(prop['step15_gppmin'])+'\n')
        grompp.Grompp(properties=prop['step15_gppmin'], **paths['step15_gppmin']).launch()

        out_log.info('step16: mdmin ------- Running: minimization')
        fu.create_dir(prop['step16_mdmin']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step16_mdmin'])+'\nProperties:\n'+str(prop['step16_mdmin'])+'\n')
        mdrun.Mdrun(properties=prop['step16_mdmin'], **paths['step16_mdmin']).launch()

        out_log.info('Step17: genrestr - Create restrain topology')
        fu.create_dir(prop['step17_genrestr']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step17_genrestr'])+'\nProperties:\n'+str(prop['step17_genrestr'])+'\n')
        genrestr.Genrestr(properties=prop['step17_genrestr'], **paths['step17_genrestr']).launch()

        out_log.info('step18: gppsa ------- Preprocessing: simulated annealing')
        fu.create_dir(prop['step18_gppsa']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step18_gppsa'])+'\nProperties:\n'+str(prop['step18_gppsa'])+'\n')
        grompp.Grompp(properties=prop['step18_gppsa'], **paths['step18_gppsa']).launch()

        out_log.info('step19: mdsa -------- Running: simulated annealing')
        fu.create_dir(prop['step19_mdsa']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step19_mdsa'])+'\nProperties:\n'+str(prop['step19_mdsa'])+'\n')
        mdrun.Mdrun(properties=prop['step19_mdsa'], **paths['step19_mdsa']).launch()

        out_log.info('step20: gppnvt_1000 - Preprocessing: nvt constant number of molecules, volume and temp')
        fu.create_dir(prop['step20_gppnvt_1000']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step20_gppnvt_1000'])+'\nProperties:\n'+str(prop['step20_gppnvt_1000'])+'\n')
        grompp.Grompp(properties=prop['step20_gppnvt_1000'], **paths['step20_gppnvt_1000']).launch()

        out_log.info('step21: mdnvt_1000 -- Running: nvt constant number of molecules, volume and temp')
        fu.create_dir(prop['step21_mdnvt_1000']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step21_mdnvt_1000'])+'\nProperties:\n'+str(prop['step21_mdnvt_1000'])+'\n')
        mdrun.Mdrun(properties=prop['step21_mdnvt_1000'], **paths['step21_mdnvt_1000']).launch()

        out_log.info('Step22: genrestr - Create restrain topology')
        fu.create_dir(prop['step22_genrestr']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step22_genrestr'])+'\nProperties:\n'+str(prop['step22_genrestr'])+'\n')
        genrestr.Genrestr(properties=prop['step22_genrestr'], **paths['step22_genrestr']).launch()

        out_log.info('step23: gppnvt_800 -- Preprocessing: nvt constant number of molecules, volume and temp')
        fu.create_dir(prop['step23_gppnvt_800']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step23_gppnvt_800'])+'\nProperties:\n'+str(prop['step23_gppnvt_800'])+'\n')
        grompp.Grompp(properties=prop['step23_gppnvt_800'], **paths['step23_gppnvt_800']).launch()

        out_log.info('step24: mdnvt_800 --- Running: nvt constant number of molecules, volume and temp')
        fu.create_dir(prop['step24_mdnvt_800']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step24_mdnvt_800'])+'\nProperties:\n'+str(prop['step24_mdnvt_800'])+'\n')
        mdrun.Mdrun(properties=prop['step24_mdnvt_800'], **paths['step24_mdnvt_800']).launch()

        out_log.info('Step25: genrestr - Create restrain topology')
        fu.create_dir(prop['step25_genrestr']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step25_genrestr'])+'\nProperties:\n'+str(prop['step25_genrestr'])+'\n')
        genrestr.Genrestr(properties=prop['step25_genrestr'], **paths['step25_genrestr']).launch()

        out_log.info('step26: gppnpt_500 -- Preprocessing: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step26_gppnpt_500']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step26_gppnpt_500'])+'\nProperties:\n'+str(prop['step26_gppnpt_500'])+'\n')
        grompp.Grompp(properties=prop['step26_gppnpt_500'], **paths['step26_gppnpt_500']).launch()

        out_log.info('step27: mdnpt_500 --- Running: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step27_mdnpt_500']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step27_mdnpt_500'])+'\nProperties:\n'+str(prop['step27_mdnpt_500'])+'\n')
        mdrun.Mdrun(properties=prop['step27_mdnpt_500'], **paths['step27_mdnpt_500']).launch()

        out_log.info('Step28: genrestr - Create restrain topology')
        fu.create_dir(prop['step28_genrestr']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step28_genrestr'])+'\nProperties:\n'+str(prop['step28_genrestr'])+'\n')
        genrestr.Genrestr(properties=prop['step28_genrestr'], **paths['step28_genrestr']).launch()

        out_log.info('step29: gppnpt_300 -- Preprocessing: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step29_gppnpt_300']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step29_gppnpt_300'])+'\nProperties:\n'+str(prop['step29_gppnpt_300'])+'\n')
        grompp.Grompp(properties=prop['step29_gppnpt_300'], **paths['step29_gppnpt_300']).launch()

        out_log.info('step30: mdnpt_300 --- Running: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step30_mdnpt_300']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step30_mdnpt_300'])+'\nProperties:\n'+str(prop['step30_mdnpt_300'])+'\n')
        mdrun.Mdrun(properties=prop['step30_mdnpt_300'], **paths['step30_mdnpt_300']).launch()

        out_log.info('Step31: genrestr - Create restrain topology')
        fu.create_dir(prop['step31_genrestr']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step31_genrestr'])+'\nProperties:\n'+str(prop['step31_genrestr'])+'\n')
        genrestr.Genrestr(properties=prop['step31_genrestr'], **paths['step31_genrestr']).launch()

        out_log.info('step32: gppnpt_200 -- Preprocessing: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step32_gppnpt_200']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step32_gppnpt_200'])+'\nProperties:\n'+str(prop['step32_gppnpt_200'])+'\n')
        grompp.Grompp(properties=prop['step32_gppnpt_200'], **paths['step32_gppnpt_200']).launch()

        out_log.info('step33: mdnpt_200 --- Running: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step33_mdnpt_200']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step33_mdnpt_200'])+'\nProperties:\n'+str(prop['step33_mdnpt_200'])+'\n')
        mdrun.Mdrun(properties=prop['step33_mdnpt_200'], **paths['step33_mdnpt_200']).launch()

        out_log.info('Step34: genrestr - Create restrain topology')
        fu.create_dir(prop['step34_genrestr']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step34_genrestr'])+'\nProperties:\n'+str(prop['step34_genrestr'])+'\n')
        genrestr.Genrestr(properties=prop['step34_genrestr'], **paths['step34_genrestr']).launch()

        out_log.info('step35: gppnpt_100 -- Preprocessing: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step35_gppnpt_100']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step35_gppnpt_100'])+'\nProperties:\n'+str(prop['step35_gppnpt_100'])+'\n')
        grompp.Grompp(properties=prop['step35_gppnpt_100'], **paths['step35_gppnpt_100']).launch()

        out_log.info('step36: mdnpt_100 --- Running: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step36_mdnpt_100']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step36_mdnpt_100'])+'\nProperties:\n'+str(prop['step36_mdnpt_100'])+'\n')
        mdrun.Mdrun(properties=prop['step36_mdnpt_100'], **paths['step36_mdnpt_100']).launch()

        out_log.info('Step37: genrestr - Create restrain topology')
        fu.create_dir(prop['step37_genrestr']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step37_genrestr'])+'\nProperties:\n'+str(prop['step37_genrestr'])+'\n')
        genrestr.Genrestr(properties=prop['step37_genrestr'], **paths['step37_genrestr']).launch()

        out_log.info('step38: gppnpt ------ Preprocessing: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step38_gppnpt']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step38_gppnpt'])+'\nProperties:\n'+str(prop['step38_gppnpt'])+'\n')
        grompp.Grompp(properties=prop['step38_gppnpt'], **paths['step38_gppnpt']).launch()

        out_log.info('step39: mdnpt ------- Running: npt constant number of molecules, pressure and temp')
        fu.create_dir(prop['step39_mdnpt']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step39_mdnpt'])+'\nProperties:\n'+str(prop['step39_mdnpt'])+'\n')
        mdrun.Mdrun(properties=prop['step39_mdnpt'], **paths['step39_mdnpt']).launch()

        out_log.info('step40: gppmd ------- Preprocessing: Free Molecular dynamics')
        fu.create_dir(prop['step40_gppmd']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step40_gppmd'])+'\nProperties:\n'+str(prop['step40_gppmd'])+'\n')
        grompp.Grompp(properties=prop['step40_gppmd'], **paths['step40_gppmd']).launch()

        out_log.info('step41: md ---------- Running: Free Molecular dynamics')
        fu.create_dir(prop['step41_md']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step41_md'])+'\nProperties:\n'+str(prop['step41_md'])+'\n')
        mdrun.Mdrun(properties=prop['step41_md'], **paths['step41_md']).launch()

        out_log.info('step42: cluster ---------- Clustering MD')
        fu.create_dir(prop['step42_cluster']['path'])
        out_log.debug('\nPaths:\n'+str(paths['step42_cluster'])+'\nProperties:\n'+str(prop['step42_cluster'])+'\n')
        mdrun.Mdrun(properties=prop['step42_cluster'], **paths['step42_cluster']).launch()

        #fu.remove_temp_files(['#', '.top', '.plotscript', '.edr', '.xtc', '.itp', '.top', '.log', '.pdb', '.cpt', '.mdp', '.xvg', '.seq'])

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

        return paths['step42_cluster']['output_pdb_path']

def main():
    yaml_path=sys.argv[1]
    system=sys.argv[2]
    workflow_path=sys.arv[3]
    Workflow(yaml_path, system, workflow_path).launch()

if __name__ == '__main__':
    main()
