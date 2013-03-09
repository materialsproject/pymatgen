#!/usr/bin/env python
from __future__ import division, print_function

import unittest
import os.path
import collections

from pymatgen.io.abinitio import *

test_dir = os.path.join(os.path.dirname(__file__))

def filepath(basename):
    return os.path.join(test_dir, basename)


class WorkflowTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_pseudoecuttest(self):

        pptest_wf = PseudoEcutTest_Workflow(filepath("14si.pspnc"), range(10,40,2))

        print(repr(pptest_wf))
        print(pptest_wf)

        pptest_wf.show_inputs()

        self.assertTrue(isinstance(pptest_wf, collections.Iterable))
        self.assertTrue(pptest_wf.isnc)

