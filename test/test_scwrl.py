from os.path import join as opj
from test import fixtures as fx
from scwrl_wrapper.scwrl import Scwrl4


class TestScwrl(object):
    def setUp(self):
        fx.test_setup(self,'scwrl')

    def tearDown(self):
        fx.test_teardown(self)

    def test_launch(self):
        output_pdb_path=opj(self.properties['path'], self.properties['output_pdb_path'])
        returncode = Scwrl4(input_pdb_path=opj(self.data_dir, self.properties['input_pdb_path']),
                            output_pdb_path=output_pdb_path,
                            properties=self.properties).launch()

        assert fx.exe_success(returncode)
        assert fx.not_empty(output_pdb_path)
