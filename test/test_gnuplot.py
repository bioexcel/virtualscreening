from os.path import join as opj
from test import fixtures as fx
import numpy as np
from gnuplot_wrapper.gnuplot import Gnuplot


class TestEditconf(object):
    def setUp(self):
        fx.test_setup(self,'gnuplot')

    def tearDown(self):
        fx.test_teardown(self)

    def test_launch(self):
        output_png_path = opj(self.properties['path'], self.properties['output_png_path'])
        returncode = Gnuplot(input_xvg_path_dict={u'A.Lys58Glu': np.array([[ 0.,2.3958132],[ 1.,2.4249675]]),
                                                  u'A.Thr74Ala': np.array([[ 0.,2.4072435],[ 1.,2.4521089]])},
                             output_png_path=output_png_path,
                             properties=self.properties).launch()
        assert fx.exe_success(returncode)
        assert fx.not_empty(output_png_path)
