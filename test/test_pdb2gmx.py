from os.path import join as opj
from test import fixtures as fx
from gromacs_wrapper.pdb2gmx import Pdb2gmx


class TestPdb2gmx(object):
    def setUp(self):
        fx.test_setup(self,'pdb2gmx')

    def tearDown(self):
        fx.test_teardown(self)

    def test_launch(self):
        output_gro_path=opj(self.properties['path'], self.properties['output_gro_path'])
        output_top_zip_path=opj(self.properties['path'], self.properties['output_top_zip_path'])
        returncode = Pdb2gmx(input_structure_pdb_path=opj(self.data_dir, self.properties['input_structure_pdb_path']),
                             output_gro_path=output_gro_path,
                             output_top_zip_path=output_top_zip_path,
                             properties=self.properties).launch()

        assert fx.exe_success(returncode)
        assert fx.not_empty(output_gro_path)
        assert fx.not_empty(output_top_zip_path)
