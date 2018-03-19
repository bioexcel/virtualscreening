from os.path import join as opj
from test import fixtures as fx
from gromacs_wrapper.mdrun import Mdrun


class TestMdrun(object):
    def setUp(self):
        fx.test_setup(self,'mdrun')

    def tearDown(self):
        fx.test_teardown(self)

    def test_launch(self):
        output_trr_path = opj(self.properties['path'], self.properties['output_trr_path'])
        output_gro_path = opj(self.properties['path'], self.properties['output_gro_path'])
        returncode = Mdrun(input_tpr_path=opj(self.data_dir, self.properties['input_tpr_path']),
                              output_trr_path=output_trr_path,
                              output_gro_path=output_gro_path,
                              properties=self.properties).launch()
        assert fx.exe_success(returncode)
        assert fx.not_empty(output_trr_path)
        assert fx.not_empty(output_gro_path)
