from __future__ import division, print_function, unicode_literals

__author__ = 'setten'

import numpy
from pymatgen.util.testing import PymatgenTest
from pymatgen.util.convergence import determine_convergence


class ConvergenceTest(PymatgenTest):

    def test_determine_convergence(self):
        self.maxDiff = None
        xs = [1, 2, 3, 4, 5, 6]
        # a converging example:
        ys = [4, 5, 6, 6, 6, 6]
        self.assertEqual(determine_convergence(xs, ys, name='name', tol=0.1, plots=False),
                         [True, 4, 6, 3, 6.0, 0.091269841269839724])
        #self.assertTrue(os.path.isfile('name.fitdat'))
        #self.assertTrue(os.path.isfile('plot-fits'))
        # another converging example
        ys = [1/1, 1/2, 1/3, 1/4, 1/5, 1/6]
        self.assertEqual(determine_convergence(xs, ys, name='name', tol=0.3, plots=False),
                         [True, 3, 0.3333333333333333, 2, 0.0, -0.12813051146384496])
        # a non converging example
        ys = [4, 5, 6, 7, 8, 9]
        self.assertEqual(determine_convergence(xs, ys, name='name', tol=0.01, plots=False),
                         [False, numpy.inf, None, None, 14.607906815185412, None])
        # another non converging example
        ys = [4, 5, 4, 5, 4, 5]
        self.assertEqual(determine_convergence(xs, ys, name='name', tol=0.01, plots=False),
                         [False, numpy.inf, None, None, 11.368169147574115, None])
        # os.remove('name.fitdat')
        # os.remove('plot-fits')
        # files = os.listdir('.')
        # for f in files:
        #    if f.startswith("convdat.") or f.endswith(".gif"):
        #        os.remove(f)
