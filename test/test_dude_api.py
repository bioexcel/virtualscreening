from os.path import join as opj
from test import fixtures as fx
from dude_api.dude import Dude


class TestDudeApi(object):
    def setUp(self):
        fx.test_setup(self,'dude')

    def tearDown(self):
        pass
        #fx.test_teardown(self)

    def test_launch(self):
        output_sdf_path=opj(self.properties['path'], self.properties['output_sdf_path'])
        Dude(output_sdf_path=output_sdf_path, properties=self.properties).get_decoys_from_pdb()
        assert fx.not_empty(output_sdf_path)
