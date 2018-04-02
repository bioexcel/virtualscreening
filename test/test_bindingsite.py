from os.path import join as opj
from test import fixtures as fx
from bindingsite.bindingsite import BindingSite


class TestBindingSite(object):
    def setUp(self):
        fx.test_setup(self,'bindingsite')

    def tearDown(self):
        pass
        #fx.test_teardown(self)

    def test_launch(self):
        output_pdb_path=opj(self.properties['path'], self.properties['output_pdb_path'])

        BindingSite(pdb_code        = self.properties['pdb_code'],
                    pdb_chain       = self.properties['pdb_chain'],
                    output_pdb_path = opj(self.properties['path'], self.properties['output_pdb_path']),
                    properties      = self.properties).launch()

        assert fx.not_empty(output_pdb_path)
