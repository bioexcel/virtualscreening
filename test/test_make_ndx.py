from os.path import join as opj
from test import fixtures as fx
from gromacs_wrapper.make_ndx import MakeNdx


class TestMakeNdx(object):
    def setUp(self):
        fx.test_setup(self,'make_ndx')

    def tearDown(self):
        fx.test_teardown(self)

    def test_launch(self):
        output_ndx_path = opj(self.properties['path'], self.properties['output_ndx_path'])
        returncode = MakeNdx(input_structure_path=opj(self.data_dir, self.properties['input_structure_path']),
                            output_ndx_path=output_ndx_path,
                            properties=self.properties).launch()
        assert fx.exe_success(returncode)
        assert fx.not_empty(output_ndx_path)
