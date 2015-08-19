# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

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
from pymatgen.io.abinitio.pseudos import Pseudo

_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')


def ref_file(filename):
    return os.path.join(_test_dir, filename)


class FakeAbinitInput(object):
    """Emulate an Abinit input."""
    @lazy_property
    def pseudos(self):
        return [Pseudo.as_pseudo(ref_file("14si.pspnc"))]

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
    """Provides helper function for testing Abinit flows."""
    MANAGER = """\
policy:
    autoparal: 1
qadapters:
    - priority: 1
      queue:
        qtype: slurm
        qname: Oban
        qparams:
            mail_user: nobody@nowhere
      limits:
        timelimit: 0:20:00
        min_cores: 4
        max_cores: 12
        #condition: {"$eq": {omp_threads: 2}}
      hardware:
        num_nodes: 10
        sockets_per_node: 1
        cores_per_socket: 2
        mem_per_node: 4 Gb
      job:
        modules:
            - intel/compilerpro/13.0.1.117
            - fftw3/intel/3.3
        shell_env:
            PATH: /home/user/tmp_intel13/src/98_main/:/home/user//NAPS/intel13/bin:$PATH
            LD_LIBRARY_PATH: /home/user/NAPS/intel13/lib:$LD_LIBRARY_PATH
        mpi_runner: mpirun

# Connection to the MongoDb database (optional) 
db_connector:
    database: abinit
    collection: test
    #host: 0.0.0.0 
    #port: 8080 
    #user: gmatteo
    #password: helloworld
"""
    def setUp(self):
        """Initialization phase."""
        super(FlowUnitTest, self).setUp()

        # Temporary directory for the flow.
        self.workdir = tempfile.mkdtemp()

        # Create the TaskManager.
        self.manager = TaskManager.from_string(self.MANAGER)

        # Fake input file
        self.fake_input = FakeAbinitInput()

    def tearDown(self):
        """Delete workdir"""
        shutil.rmtree(self.workdir)


class FlowTest(FlowUnitTest):

    def test_base(self):
        """Testing Flow..."""
        aequal, atrue, afalse = self.assertEqual, self.assertTrue, self.assertFalse
        flow = Flow(workdir=self.workdir, manager=self.manager)

        # Build a work with a task
        work = flow.register_task(self.fake_input)
        assert work.is_work
        task0_w0 = work[0]
        atrue(task0_w0.is_task)
        afalse(task0_w0.has_subnodes)
        print(task0_w0.status.colored)
        atrue(len(flow) == 1) 
        aequal(flow.num_tasks, 1)
        atrue(flow.has_db) 

        #print(task0_w0.input_structure)
        print(task0_w0.make_input)

        assert flow.select_tasks(nids=task0_w0.node_id)[0] == task0_w0
        assert flow.select_tasks(wslice=slice(0,1,1)) == [task0_w0]

        # Build a workflow containing two tasks depending on task0_w0
        work = Work()
        atrue(work.is_work)
        atrue(work.has_subnodes)
        work.register(self.fake_input)
        work.register(self.fake_input)
        aequal(len(work), 2)

        flow.register_work(work, deps={task0_w0: "WFK"})
        atrue(flow.is_flow)
        atrue(flow.has_subnodes)
        aequal(len(flow), 2)

        # Add another work without dependencies.
        task0_w2 = flow.register_task(self.fake_input)[0]
        atrue(len(flow) == 3)
        afalse(flow.is_work)

        # Allocate internal tables
        flow.allocate()

        # Check dependecies.
        atrue(flow[1].depends_on(task0_w0))
        atrue(flow[1][0].depends_on(task0_w0))
        atrue(flow[1][0] in task0_w0.get_children())
        atrue(task0_w0 in flow[1][0].get_parents())
        afalse(flow[2][0].depends_on(task0_w0))
        afalse(flow[2][0] in task0_w0.get_children())
        afalse(task0_w0 in flow[2][0].get_parents())
        aequal(flow[1].pos, 1)
        aequal(flow[1][0].pos, (1, 0))
        aequal(flow[2][0].pos, (2, 0))

        afalse(flow.all_ok)
        aequal(flow.num_tasks, 4)
        aequal(flow.ncores_inuse, 0)

        # API for iterations
        aequal(len(list(flow.iflat_tasks(status="Initialized"))), sum(len(work) for work in flow))
        aequal(list(flow.iflat_tasks(nids=task0_w0.node_id)), [task0_w0])

        aequal([task0_w0], flow.tasks_from_nids(task0_w0.node_id))
        aequal([(0, 0)], flow.wti_from_nids(task0_w0.node_id))
        aequal([task0_w2], flow.tasks_from_nids([task0_w2.node_id]))
        aequal([(2, 0)], flow.wti_from_nids([task0_w2.node_id]))

        # Check for deadlocks
        flow.check_dependencies()

        # Save the flow in pickle format.
        flow.build_and_pickle_dump()

        # Find the pickle file in workdir and recreate the flow.
        same_flow = Flow.pickle_load(self.workdir)

        same_flow = Flow.pickle_load(self.workdir)
        aequal(same_flow, flow)

        # Test show_status
        flow.show_status()


#class BandStructureFlowTest(FlowUnitTest):
#    def test_base(self):
#        """Testing bandstructure flow..."""
#        flow = bandstructure_flow(self.workdir, self.manager, self.fake_input, self.fake_input)


if __name__ == '__main__':
    import unittest
    unittest.main()
