# coding: utf-8

from __future__ import unicode_literals, division, print_function

import os
import tempfile
import shutil

from pymatgen.util.testing import PymatgenTest
from monty.functools import lazy_property
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.abinitio import *
from pymatgen.io.abinitio.flows import *
from pymatgen.io.abinitio.tasks import *

_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')


def ref_file(filename):
    return os.path.join(_test_dir, filename)


class FakeAbinitInput(object):
    """Emulate an Abinit input."""
    @lazy_property
    def pseudos(self):
        return ref_file("14si.pspnc")

    @lazy_property
    def structure(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                          [1.9200989668, 3.3257101909, 0.00],
                          [0.00, -2.2171384943, 3.1355090603]])
        return Structure(lattice, ["Si", "Si"], coords)


class FlowUnitTest(PymatgenTest):
    """Ppovides helper function for testing Abinit flows."""
    def setUp(self):
        """Initialization phase."""
        super(FlowUnitTest, self).setUp()

        # Temporary directory for the flow.
        self.workdir = tempfile.mkdtemp()

        # Create the TaskManager.
        self.manager = TaskManager.from_file(os.path.join(_test_dir, "taskmanager.yml"))

        # Fake input file
        self.fake_input = FakeAbinitInput()

    def tearDown(self):
        """Delete workdir"""
        shutil.rmtree(self.workdir)


class AbinitFlowTest(FlowUnitTest):

    def test_base(self):
        """Testing AbinitFlow..."""
        flow = AbinitFlow(workdir=self.workdir, manager=self.manager)

        # Build a workflow with a task
        task0_w0 = flow.register_task(self.fake_input)
        self.assertTrue(len(flow) == 1)
        self.assertEqual(flow.num_tasks, 1)

        # Build a workflow containing two tasks depending on task0_w0
        work = Work()
        work.register(self.fake_input)
        work.register(self.fake_input)
        self.assertTrue(len(work) == 2)

        flow.register_work(work, deps={task0_w0: "WFK"})
        self.assertTrue(len(flow) == 2)


        # Add another work without dependencies.
        task0_w2 = flow.register_task(self.fake_input)
        self.assertTrue(len(flow) == 3)

        # Allocate internal tables
        flow.allocate()

        # Check dependecies.
        self.assertTrue(flow[1].depends_on(task0_w0))
        self.assertTrue(flow[1][0].depends_on(task0_w0))
        self.assertFalse(flow[2][0].depends_on(task0_w0))
        self.assertEqual(flow[1].pos, 1)
        self.assertEqual(flow[1][0].pos, (1, 0))
        self.assertEqual(flow[2][0].pos, (2, 0))

        self.assertFalse(flow.all_ok)
        self.assertEqual(flow.num_tasks, 4)
        self.assertEqual(flow.ncores_inuse, 0)

        # Check for deadlocks
        flow.check_dependencies()

        # TODO: Fix pickle for flow. Test is temporarily disabled for now by the Hulk.
        # Save the flow in pickle format.
        flow.build_and_pickle_dump()

        # Find the pickle file in workdir and recreate the flow.
        same_flow = AbinitFlow.pickle_load(self.workdir)

        self.assertEqual(same_flow, flow)

        # Test show_status
        flow.show_status()


#class BandStructureFlowTest(FlowUnitTest):
#    def test_base(self):
#        """Testing bandstructure flow..."""
#        flow = bandstructure_flow(self.workdir, self.manager, self.fake_input, self.fake_input)


if __name__ == '__main__':
    import unittest
    unittest.main()
