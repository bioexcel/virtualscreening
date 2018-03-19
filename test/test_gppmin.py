from os.path import join as opj
from test import fixtures as fx
from gromacs_wrapper.grompp import Grompp


class TestGrompp(object):
    def setUp(self):
        fx.test_setup(self,'gppmin')

    def tearDown(self):
        fx.test_teardown(self)

    def test_launch(self):
        output_tpr_path = opj(self.properties['path'], self.properties['output_tpr_path'])
        returncode = Grompp(input_gro_path=opj(self.data_dir, self.properties['input_gro_path']),
                              input_top_zip_path=opj(self.data_dir, self.properties['input_top_zip_path']),
                              output_tpr_path=output_tpr_path,
                              properties=self.properties).launch()

        assert fx.exe_success(returncode)
        assert fx.not_empty(output_tpr_path)
