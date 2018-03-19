#! /usr/bin/python
"""Python wrapper module for CMIP
"""
__author__ = "gelpi"
__date__ = "$24-nov-2017 12:30:59$"

import time
import sys
import os
import re
import configuration.settings as settings
import tools.file_utils as fu
import cmip_wrapper.CMIPWrapper as CW
import cmip_wrapper.prepPDBWrapper as PPW
from shutil import copyfile

def main():
    start_time = time.time()
    yaml_path=sys.argv[1]
    system=sys.argv[2]
    conf = settings.YamlReader(yaml_path, system)
    workflow_path = conf.properties[system]['workflow_path']
    fu.create_dir(os.path.abspath(workflow_path))
    out_log, _ = fu.get_logs(path=workflow_path, console=True)
    paths = conf.get_paths_dic()
    props = conf.get_prop_dic(global_log=out_log)

    out_log.info('')
    out_log.info('_______TEST CMIP TITRATION WORKFLOW_______')
    out_log.info('')

    stepid = 'step1_mmbpdb'
    out_log.info(stepid)
    props['step'] = stepid

    # For testing only
    fu.create_dir(props[stepid]['path'])

    copyfile(paths['step1_mmbpdb']['input_test_pdb_path'], paths['step1_mmbpdb']['output_pdb_path'])
    out_log.info('Test input file copy from: '+paths['step1_mmbpdb']['input_test_pdb_path']+' to: '+paths['step1_mmbpdb']['output_pdb_path'])

    stepid = 'step2_preppdbCMIP'
    out_log.info(stepid)
    props['step'] = stepid
    fu.create_dir(props[stepid]['path'])
    out_log.info('Step 2 prep pdb CMIP create dir: '+props[stepid]['path'])
    PPW.prepPDBWrapper(paths[stepid],props[stepid]).launch()

    # stepid = 'step3_CMIPTitration'
    # out_log.info(stepid)
    # props['step'] = stepid
    #
    # fu.create_dir(props[stepid]['path'])
    #
    # cw = CW.CMIPWrapper(paths[stepid],props[stepid])
    # #Delaying removal of temporary dir until output is processed
    # cw.launch(False)
    # out_log.debug ("Temp directory: "+ cw.run.tmpdir)
    # out_log.info(stepid + " CMIP Output")
    #
    # if cw.result.errstr:
    #     out_log.error(cw.result.errstr)
    #     sys.exit(1)
    # else:
    #     out_log.info(cw.result.output)
    # out_log.info('')
    #
    # outpdbFn = cw.result.getFileName('outpdb')+".pdb"
    # out_log.info('Building output pdb file at')
    # out_log.info(paths[stepid]['output_pdb_path'])
    #
    # out_log.info('Base File: '+ paths[stepid]['input_pdb_path'])
    # out_log.info('Titration File: '+ outpdbFn)
    #
    # outFh = open(paths[stepid]['output_pdb_path'],"w")
    # baseFh = open(paths[stepid]['input_pdb_path'],"r")
    # addFh = open(outpdbFn,"r")
    #
    # outFh.write(baseFh.read())
    # baseFh.close()
    # for line in addFh:
    #     if re.match('^ATOM',line):
    #         outFh.write(line)
    # outFh.close()
    # addFh.close()
    #
    # out_log.info('')
    #
    # cw.run.rmTmpDir()
    #
    elapsed_time = time.time() - start_time
    out_log.info('')
    out_log.info('')
    out_log.info('Execution successful: ')
    out_log.info('  Workflow_path: '+workflow_path)
    out_log.info('  Config File: '+yaml_path)
    out_log.info('  System: '+system)
    if len(sys.argv) >= 4:
        out_log.info('  Nodes: '+sys.argv[3])
    out_log.info('')
    out_log.info('Elapsed time: '+str(elapsed_time)+' seconds')
    out_log.info('')

if __name__ == "__main__":
    main()
