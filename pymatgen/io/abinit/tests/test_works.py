# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

import os.path
import collections

from tempfile import mkdtemp
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinit import *

_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files', "abinit")

def ref_file(filename):
    return os.path.join(_test_dir, filename)


#class WorkTestCase(PymatgenTest):
#
#    def test_pseudoconvergence(self):
#        workdir = mkdtemp(prefix="test_pseudoconvergence")
#        manager = TaskManager.sequential()
#        pseudo = ref_file("14si.pspnc")
#        ecut_list = range(10, 40, 2)
#
#        pptest_wf = PseudoConvergence(pseudo, ecut_list, atols_mev=(10, 1, 0.1), 
#                                      workdir=workdir, manager=manager)
#        pptest_wf.allocate()
#
#        # Test pickle
#        # FIXME: protocol 2 does not work due to __new__
#        self.serialize_with_pickle(pptest_wf, protocols=[0,1], test_eq=True)
#
#        print(repr(pptest_wf))
#        print(pptest_wf)
#
#        self.assertTrue(isinstance(pptest_wf, collections.Iterable))
#        self.assertTrue(pptest_wf.isnc)
#
#        pptest_wf.build()
#        pptest_wf.rmtree()


if __name__ == "__main__":
    import unittest2 as unittest
    unittest.main()
