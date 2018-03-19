"""Boiler plate functions for testsys
"""
import os
import sys
import shutil
from os.path import join as opj
from configuration import settings
from tools import file_utils as fu
def test_setup(test_object, dict_key):
    test_object.data_dir = opj(os.path.dirname(os.path.abspath(sys.modules[__name__].__file__)),'data')
    test_object.yaml_path= opj(test_object.data_dir, 'conf.yaml')
    test_object.system=os.getenv('testsys')
    if test_object.system is None:
        print 'WARNING: "testsys" env variable should be set, "linux" will be used by default value.'
        print '     Please, try: "export testsys=linux"'
        test_object.system='linux'
    conf = settings.YamlReader(test_object.yaml_path, test_object.system)
    test_object.properties = conf.get_prop_dic()[dict_key]
    fu.create_dir(test_object.properties['path'])

def test_teardown(test_object):
    shutil.rmtree(test_object.properties['workflow_path'])
    fu.remove_temp_files(['#', '.top', '.plotscript', '.edr', '.xtc', '.itp', '.top', '.log', '.pdb', '.cpt', '.mdp', '.xvg', '.grp', '.ndx'])

def exe_success(return_code):
    return return_code == 0

def not_empty(file_path):
    return ( os.path.isfile(file_path) and os.path.getsize(file_path) > 0 )
