from os.path import join as opj
from test import fixtures as fx
from bindingsite.box import Box


class TestBindingSite(object):
    def setUp(self):
        fx.test_setup(self,'box')

    def tearDown(self):
        pass
        #fx.test_teardown(self)

    def test_launch(self):
        output_pdb_path=opj(self.properties['path'], self.properties['output_pdb_path'])

        Box(input_pdb_path  = opj(self.data_dir, self.properties['input_pdb_path']),
		    resid_pdb_path  = opj(self.data_dir, self.properties['resid_pdb_path']),
                    output_pdb_path = opj(self.properties['path'], self.properties['output_pdb_path']),
                    properties      = self.properties).launch()

        assert fx.not_empty(output_pdb_path)
