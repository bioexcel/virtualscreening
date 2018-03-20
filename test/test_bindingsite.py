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
        
        print opj(self.data_dir, self.properties['input_pdb_path']), opj(self.data_dir,self.properties['clusterPDBs_zip_path']), output_pdb_path, self.properties
        
        bs =  BindingSite(input_pdb_path       = opj(self.data_dir, self.properties['input_pdb_path']),
                          clusterPDBs_zip_path = opj(self.data_dir,self.properties['clusterPDBs_zip_path']),
                          output_pdb_path      = output_pdb_path,
                          properties           = self.properties)
        bs.launch()

        assert fx.not_empty(output_pdb_path)
