from os.path import join as opj
from test import fixtures as fx
from gromacs_wrapper.editconf import Editconf


class TestEditconf(object):
    def setUp(self):
        fx.test_setup(self,'editconf')

    def tearDown(self):
        fx.test_teardown(self)

    def test_launch(self):
        output_gro_path = opj(self.properties['path'], self.properties['output_gro_path'])
        returncode = Editconf(input_gro_path=opj(self.data_dir, self.properties['input_gro_path']),
                              output_gro_path=output_gro_path,
                              properties=self.properties).launch()
        assert fx.exe_success(returncode)
        assert fx.not_empty(output_gro_path)
