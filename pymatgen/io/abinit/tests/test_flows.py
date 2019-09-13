# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import tempfile
import shutil

from pymatgen.util.testing import PymatgenTest
from monty.functools import lazy_property
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.abinit import *
from pymatgen.io.abinit.flows import *
from pymatgen.io.abinit.works import *
from pymatgen.io.abinit.tasks import *
from pymatgen.io.abinit.pseudos import Pseudo

_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                         'test_files', "abinit")


def ref_file(filename):
    return os.path.join(_test_dir, filename)


class FakeAbinitInput:
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

    def get(self, key, default=None):
        """The real AbinitInput is a dict-like object."""
        if default is not None: return default
        return key


class FlowUnitTest(PymatgenTest):
    """Provides helper function for testing Abinit flows."""
    MANAGER = """\
policy:
    autoparal: 1
qadapters:
    - &batch
      priority: 1
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

batch_adapter: *batch
"""
    def setUp(self):
        """Initialization phase."""
        super().setUp()

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
        assert flow.isinstance(Flow)
        assert not flow.isinstance(None)
        assert not flow.has_scheduler

        # Build a work with a task
        work = flow.register_task(self.fake_input)
        atrue(work.is_work)
        atrue(len(work.color_hex) == 7)
        atrue(work.color_hex.startswith("#"))
        task0_w0 = work[0]
        atrue(task0_w0.is_task)
        print(task0_w0.status.colored)
        atrue(len(flow) == 1)
        aequal(flow.num_tasks, 1)
        atrue(flow.has_db)

        #print(task0_w0.input_structure)
        print(task0_w0.make_input)

        # Task history
        assert len(task0_w0.history) == 0
        task0_w0.history.info("Hello %s", "world")
        assert len(task0_w0.history) == 1
        print(task0_w0.history)
        record = task0_w0.history.pop()
        print(record, repr(record))
        assert record.get_message(asctime=False) == "Hello world"
        assert len(task0_w0.history) == 0
        assert flow.select_tasks(nids=task0_w0.node_id)[0] == task0_w0
        assert flow.select_tasks(wslice=slice(0,1,1)) == [task0_w0]
        assert flow.select_tasks(task_class="DfptTask") == []
        assert flow.get_task_scfcycles() == []

        # Build a workflow containing two tasks depending on task0_w0
        work = Work()
        atrue(work.is_work)
        work.register(self.fake_input)
        work.register(self.fake_input)
        aequal(len(work), 2)

        flow.register_work(work, deps={task0_w0: "WFK"})
        atrue(flow.is_flow)
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
        aequal(flow.ncores_used, 0)

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
        aequal(same_flow, flow)

        # to/from string
        # FIXME This does not work with py3k
        #s = flow.pickle_dumps(protocol=0)
        #same_flow = Flow.pickle_loads(s)
        #aequal(same_flow, flow)

        self.assertMSONable(flow)

        flow.show_info()
        flow.show_summary()
        flow.show_inputs()
        flow.show_inputs(varnames="znucl")

        df_vars = flow.get_vars_dataframe("ecut", "acell")
        atrue("ecut" in df_vars)

        # Test show_status
        flow.show_status()
        flow.show_tricky_tasks()
        flow.show_event_handlers()

    def test_workdir(self):
        """Testing if one can use workdir=None in flow.__init__ and then flow.allocate(workdir)."""
        flow = Flow(workdir=None, manager=self.manager)
        flow.register_task(self.fake_input)
        #flow.register_work(work)
        work = Work()
        work.register_scf_task(self.fake_input)
        flow.register_work(work)

        # If flow.workdir is None, we should used flow.allocate(workdir)
        with self.assertRaises(RuntimeError): flow.allocate()

        tmpdir = tempfile.mkdtemp()
        flow.allocate(workdir=tmpdir)

        print(flow)
        assert len(flow) == 2

        flow.build()

        for i, work in enumerate(flow):
            assert work.workdir == os.path.join(tmpdir, "w%d" % i)
            for t, task in enumerate(work):
                assert task.workdir == os.path.join(work.workdir, "t%d" % t)


class TestFlowInSpectatorMode(FlowUnitTest):

    def test_spectator(self):
        flow = Flow(workdir=self.workdir, manager=self.manager)

        work0 = Work()
        gs_task = work0.register_scf_task(self.fake_input)
        assert gs_task.isinstance(ScfTask)
        assert gs_task.isinstance("ScfTask")
        task = work0.register_scf_task(self.fake_input)
        assert task.is_abinit_task
        assert not task.is_optic_task
        assert not task.is_anaddb_task

        work1 = Work()
        work1.register_scf_task(self.fake_input)

        flow.register_work(work0)
        flow.register_work(work1)

        flow.disconnect_signals()
        flow.disconnect_signals()

        flow.connect_signals()
        flow.connect_signals()

        for mode in [False, True]:
            flow.set_spectator_mode(mode=mode)
            assert flow.in_spectator_mode == mode
            for node in flow.iflat_nodes():
                assert node.in_spectator_mode == mode

        assert len(list(flow.iflat_nodes())) == 1 + len(flow.works) + sum(len(work) for work in flow)
        assert flow.node_from_nid(flow.node_id) == flow

        flow.set_spectator_mode(mode=False)
        flow.build_and_pickle_dump()

        # picke load always returns a flow in spectator mode.
        flow = Flow.pickle_load(flow.workdir)
        assert flow.in_spectator_mode

        #with self.assertRaises(flow.SpectatorError): flow.pickle_dump()
        #with self.assertRaises(flow.SpectatorError): flow.make_scheduler().start()

        work = flow[0]
        assert work.send_signal(work.S_OK) is None
        #with self.assertRaises(work.SpectatorError): work.on_ok()
        #with self.assertRaises(work.SpectatorError): work.on_all_ok()

        task = work[0]
        assert task.send_signal(task.S_OK) is None
        #with self.assertRaises(task.SpectatorError): task._on_done()
        #with self.assertRaises(task.SpectatorError): task.on_ok()
        #with self.assertRaises(task.SpectatorError): task._on_ok()


class TestBatchLauncher(FlowUnitTest):

    def test_batchlauncher(self):
        """Testing BatchLauncher methods."""
        # Create the TaskManager.
        manager = TaskManager.from_string(self.MANAGER)
        print("batch_adapter", manager.batch_adapter)
        assert manager.batch_adapter is not None

        def build_flow_with_name(name):
            """Build a flow with workdir None and the given name."""
            flow = Flow(workdir=None, manager=self.manager)
            flow.set_name(name)

            flow.register_task(self.fake_input)
            work = Work()
            work.register_scf_task(self.fake_input)
            flow.register_work(work)

            return flow

        from pymatgen.io.abinit.launcher import BatchLauncher
        tmpdir = tempfile.mkdtemp()
        batch = BatchLauncher(workdir=tmpdir, manager=manager)
        print(batch)

        flow0 = build_flow_with_name("flow0")
        flow1 = build_flow_with_name("flow1")
        flow2_same_name = build_flow_with_name("flow1")

        batch.add_flow(flow0)

        # Cannot add the same flow twice.
        with self.assertRaises(batch.Error):
            batch.add_flow(flow0)

        batch.add_flow(flow1)

        # Cannot add two flows with the same name.
        with self.assertRaises(batch.Error):
            batch.add_flow(flow2_same_name)

        batch.submit(dry_run=True)

        for i, flow in enumerate([flow0, flow1]):
            assert flow.workdir == os.path.join(batch.workdir, "flow%d" % i)

        batch.pickle_dump()
        batch_from_pickle = BatchLauncher.pickle_load(batch.workdir)
        assert all(f1 == f2 for f1, f2 in zip(batch.flows, batch_from_pickle.flows))


if __name__ == '__main__':
    import unittest
    unittest.main()
