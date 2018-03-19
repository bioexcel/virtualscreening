from os.path import join as opj
from test import fixtures as fx
from gromacs_wrapper.rms import Rms


class TestRms(object):
    def setUp(self):
        fx.test_setup(self,'rms')

    def tearDown(self):
        fx.test_teardown(self)

    def test_launch(self):
        output_xvg_path = opj(self.properties['path'], self.properties['output_xvg_path'])
        Rms(input_gro_path=opj(self.data_dir, self.properties['input_gro_path']),
            input_trr_path=opj(self.data_dir, self.properties['input_trr_path']),
            output_xvg_path=output_xvg_path, properties=self.properties).launch()
        assert fx.not_empty(output_xvg_path)
