from os.path import join as opj
from test import fixtures as fx
from gromacs_wrapper.genrestr import Genrestr


class TestMakeNdx(object):
    def setUp(self):
        fx.test_setup(self,'genrestr')

    def tearDown(self):
        fx.test_teardown(self)

    def test_launch(self):
        output_top_zip_path = opj(self.properties['path'], self.properties['output_top_zip_path'])
        returncode = Genrestr(input_structure_path=opj(self.data_dir, self.properties['input_structure_path']),
                              input_ndx_path=opj(self.data_dir, self.properties['input_ndx_path']),
                              input_top_zip_path=opj(self.data_dir, self.properties['input_top_zip_path']),
                              output_top_zip_path=output_top_zip_path,
                              properties=self.properties).launch()
        assert fx.exe_success(returncode)
        assert fx.not_empty(output_top_zip_path)
