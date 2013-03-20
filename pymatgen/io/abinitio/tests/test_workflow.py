#!/usr/bin/env python
from __future__ import division, print_function

import unittest
import os.path
import collections
from tempfile import mkdtemp

from pymatgen.io.abinitio import *
from pymatgen.io.abinitio.workflow import Work

test_dir = os.path.join(os.path.dirname(__file__))

def filepath(basename):
    return os.path.join(test_dir, basename)


class WorkflowTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_pseudoecutconvergence(self):

        workdir = "test_pseudoecutconvergence"
        #workdir = mkdtemp()

        pptest_wf = PseudoEcutConvergence(workdir, filepath("14si.pspnc"), range(10,40,2))

        #with self.assertRaises(Workflow.Error):
        #    cannot_have_another_wf_in_same_workdir = PseudoEcutTest_Workflow(workdir, filepath("14si.pspnc"), range(10,40,2))

        print(repr(pptest_wf))
        print(pptest_wf)

        pptest_wf.show_inputs()

        self.assertTrue(isinstance(pptest_wf, collections.Iterable))
        self.assertTrue(pptest_wf.isnc)

        pptest_wf.build()

        pptest_wf.destroy()
