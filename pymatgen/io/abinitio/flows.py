"""
Abinit Flows
"""
from __future__ import division, print_function

import os
import collections
import itertools
import cPickle as pickle

try:
    from pydispatch import dispatcher
except ImportError:
    pass

#from pymatgen.util.string_utils import list_strings, pprint_table
#from pymatgen.io.abinitio.tasks import (AbinitTask, Dependency, Node, ScfTask, NscfTask, HaydockBseTask)
from pymatgen.io.abinitio.utils import Directory
from pymatgen.io.abinitio.workflows import Workflow

import logging
logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

__all__ = [
    "AbinitFlow",
    "bandstructure_flow",
    "g0w0_flow",
]


class AbinitFlow(collections.Iterable):
    """
    This object is a container of workflows. Its main task is managing the 
    possible inter-depencies among the workflows and the creation of
    dynamic worflows that are generates by callabacks registered by the user.
    """
    VERSION = "0.1"

    PICKLE_FNAME = "__AbinitFlow__.pickle"

    def __init__(self, workdir, manager):
        """
        Args:
            workdir:
                String specifying the directory where the workflows will be produced.
            manager:
                `TaskManager` object responsible for the submission of the jobs.
        """
        self.workdir = os.path.abspath(workdir)

        self.manager = manager.deepcopy()

        # List of workflows.
        self._works = []

        # List of callbacks that must be executed when the dependencies reach S_OK
        self._callbacks = []

        # Directories with (input|output|temporary) data.
        self.indir = Directory(os.path.join(self.workdir, "indata"))
        self.outdir = Directory(os.path.join(self.workdir, "outdata"))
        self.tmpdir = Directory(os.path.join(self.workdir, "tmpdata"))

    def __repr__(self):
        return "<%s at %s, workdir=%s>" % (self.__class__.__name__, id(self), self.workdir)

    def __str__(self):
        return "<%s, workdir=%s>" % (self.__class__.__name__, self.workdir)

    def __len__(self):
        return len(self.works)

    def __iter__(self):
        return self.works.__iter__()

    def __getitem__(self, slice):
        return self.works[slice]

    @property
    def works(self):
        """List of `Workflow` objects contained in self.."""
        return self._works

    @property
    def ncpus_reserved(self):
        """
        Returns the number of CPUs reserved in this moment.
        A CPUS is reserved if it's still not running but 
        we have submitted the task to the queue manager.
        """
        return sum(work.ncpus_reverved for work in self)

    @property
    def ncpus_allocated(self):
        """
        Returns the number of CPUs allocated in this moment.
        A CPU is allocated if it's running a task or if we have
        submitted a task to the queue manager but the job is still pending.
        """
        return sum(work.ncpus_allocated for work in self)

    @property
    def ncpus_inuse(self):
        """
        Returns the number of CPUs used in this moment.
        A CPU is used if there's a job that is running on it.
        """
        return sum(work.ncpus_inuse for work in self)

    def used_ids(self):
        """
        Returns a set with all the ids used so far to identify `Task` and `Workflow`.
        """
        ids = []
        for work in self:
            ids.append(work.node_id)
            for task in work:
                ids.append(task.node_id)

        used_ids = set(ids)
        assert len(ids_set) == len(ids)

        return used_ids

    def generate_new_nodeid(self):
        """Returns an unused node identifier."""
        used_ids = self.used_ids()

        for nid in itertools.count():
            if nid not in used_ids:
                return nid

    def check_status(self):
        """Check the status of the workflows in self."""
        for work in self:
            work.check_status()

    def build(self, *args, **kwargs):
        self.indir.makedirs()
        self.outdir.makedirs()
        self.tmpdir.makedirs()

        for work in self:
            work.build(*args, **kwargs)

    def build_and_pickle_dump(self, protocol=-1):
        self.build()
        self.pickle_dump(protocol=protocol)

    def pickle_dump(self, protocol=-1):
        """Save the status of the object in pickle format."""
        filepath = os.path.join(self.workdir, self.PICKLE_FNAME)

        with open(filepath, mode="w" if protocol == 0 else "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)

        # Atomic transaction.
        #filepath_new = filepath + ".new"
        #filepath_save = filepath + ".save"
        #shutil.copyfile(filepath, filepath_save)

        #try:
        #    with open(filepath_new, mode="w" if protocol == 0 else "wb") as fh:
        #        pickle.dump(self, fh, protocol=protocol)

        #    os.rename(filepath_new, filepath)
        #except IOError:
        #    os.rename(filepath_save, filepath)
        #finally:
        #    try
        #        os.remove(filepath_save)
        #    except:
        #        pass
        #    try
        #        os.remove(filepath_new)
        #    except:
        #        pass

    @staticmethod
    def pickle_load(filepath):
        with open(filepath, "rb") as fh:
            flow = pickle.load(fh)

        flow.connect_signals()
        return flow

    def register_task(self, input, deps=None, manager=None, task_class=None):
        """
        Utility function that generates a `Workflow` made of a single task

        Args:
            input:
                Abinit Input file or `Strategy` object.
            deps:
                List of `Dependency` objects specifying the dependency of this node.
                An empy list of deps implies that this node has no dependencies.
            manager:
                The `TaskManager` responsible for the submission of the task. 
                If manager is None, we use the `TaskManager` specified during the creation of the workflow.
            task_class:
                Task subclass to instantiate. Default: `AbinitTask` 

        Returns:   
            The generated `Task`.
        """
        work = Workflow(manager=manager)
        task = work.register(input, deps=deps, task_class=task_class)

        self.register_work(work)
        return task

    def register_work(self, work, deps=None, manager=None):
        """
        Register a new `Workflow` and add it to the internal list, 
        taking into account possible dependencies.

        Args:
            work:
                `Workflow` object.
            deps:
                List of `Dependency` objects specifying the dependency of this node.
                An empy list of deps implies that this node has no dependencies.
            manager:
                The `TaskManager` responsible for the submission of the task. 
                If manager is None, we use the `TaskManager` specified during the creation of the workflow.

        Returns:   
            The workflow.
        """
        # Directory of the workflow.
        work_workdir = os.path.join(self.workdir, "work_" + str(len(self)))

        # Make a deepcopy since manager is mutable and we might change it at run-time.
        manager = self.manager.deepcopy() if manager is None else manager.deepcopy()
        #print("flow manager", manager)

        work.set_workdir(work_workdir)
        work.set_manager(manager)

        self.works.append(work)

        if deps:
            deps = [Dependency(node, exts) for node, exts in deps.items()]
            work.add_deps(deps)

        return work

    def register_cbk(self, cbk, cbk_data, deps, work_class, manager=None):
        """
        Registers a callback function that will generate the Task of the Workflow.

        Args:
            cbk:
                Callback function.
            cbk_data
                Additional data pased to the callback function.
            deps:
                List of `Dependency` objects specifying the dependency of the workflow.
            work_class:
                `Workflow` class to instantiate.
            manager:
                The `TaskManager` responsible for the submission of the task. 
                If manager is None, we use the `TaskManager` specified during the creation of the `Flow`.
                                                                                                            
        Returns:   
            The `Workflow` that will be finalized by the callback.
        """
        # TODO: pass a workflow factory instead of a class
        # Directory of the workflow.
        work_workdir = os.path.join(self.workdir, "work_" + str(len(self)))
                                                                                                            
        # Make a deepcopy since manager is mutable and we might change it at run-time.
        manager = self.manager.deepcopy() if manager is None else manager.deepcopy()
                                                                                                            
        # Create an empty workflow and register the callback
        work = work_class(workdir=work_workdir, manager=manager)
        
        self._works.append(work)
                                                                                                            
        deps = [Dependency(node, exts) for node, exts in deps.items()]
        if not deps:
            raise ValueError("A callback must have deps!")

        work.add_deps(deps)

        # Wrap the callable in a Callback object and save 
        # useful info such as the index of the workflow and the callback data.
        cbk = Callback(cbk, work, deps=deps, cbk_data=cbk_data)
                                                                                                            
        self._callbacks.append(cbk)
                                                                                                            
        return work

    def allocate(self):
        for work in self:
            work.allocate()

        return self

    def show_dependencies(self):
        for work in self:
            work.show_intrawork_deps()

    def on_dep_ok(self, signal, sender):
        # TODO
        # Replace this callback with dynamic dispatch
        # on_all_S_OK for workflow
        # on_S_OK for task
        print("on_dep_ok with sender %s, signal %s" % (str(sender), signal))

        for i, cbk in enumerate(self._callbacks):

            if not cbk.handle_sender(sender):
                print("Do not handle")
                continue

            if not cbk.can_execute():
                print("cannot execute")
                continue 

            # Execute the callback to generate the workflow.
            print("about to build new workflow")
            #empty_work = self._works[cbk.w_idx]

            # TODO better treatment of ids
            # Make sure the new workflow has the same id as the previous one.
            #new_work_idx = cbk.w_idx
            work = cbk(flow=self)
            work.add_deps(cbk.deps)

            # Disable the callback.
            cbk.disable()

            # Update the database.
            self.pickle_dump()

    #def finalize(self):
    #    """This method is called when the flow is completed."""

    def connect_signals(self):
        """
        Connect the signals within the workflow.
        self is responsible for catching the important signals raised from 
        its task and raise new signals when some particular condition occurs.
        """
        # Connect the signals inside each Workflow.
        for work in self:
            work.connect_signals()

        # Observe the nodes that must reach S_OK in order to call the callbacks.
        for cbk in self._callbacks:
            for dep in cbk.deps:
                print("connecting %s \nwith sender %s, signal %s" % (str(cbk), dep.node, dep.node.S_OK))
                dispatcher.connect(self.on_dep_ok, signal=dep.node.S_OK, sender=dep.node, weak=False)

        self.show_receivers()

    def show_receivers(self, sender=dispatcher.Any, signal=dispatcher.Any):
        print("*** live receivers ***")
        for rec in dispatcher.liveReceivers(dispatcher.getReceivers(sender, signal)):
            print("receiver -->", rec)
        print("*** end live receivers ***")


class Callback(object):

    def __init__(self, func, work, deps, cbk_data):
        """
        Initialize the callback.

        Args:
            func:
                The function to execute. Must have signature .... TODO
            work:
                Reference to the `Workflow` that will be filled.
            deps:
                List of dependencies associated to the callback
            cbk_data:
                Additional data to pass to the callback.
        """
        self.func  = func
        self.work = work
        self.deps = deps
        self.cbk_data = cbk_data or {}
        self._disabled = False

    def __call__(self, flow):
        """Execute the callback."""
        if self.can_execute():
            print("in callback")
            #print("in callback with sender %s, signal %s" % (sender, signal))
            cbk_data = self.cbk_data.copy()
            #cbk_data["_w_idx"] = self.w_idx
            return self.func(flow=flow, work=self.work, cbk_data=cbk_data)
        else:
            raise Exception("Cannot execute")

    def can_execute(self):
        """True if we can execut the callback."""
        return not self._disabled and [dep.status == dep.node.S_OK  for dep in self.deps]

    def disable(self):
        """
        True if the callback has been disabled.
        This usually happens when the callback has been executed.
        """
        self._disabled = True

    def handle_sender(self, sender):
        """
        True if the callback is associated to the sender
        i.e. if the node who sent the signal appears in the 
        dependencies of the callback.
        """
        return sender in [d.node for d in self.deps]


# Factory functions.
def bandstructure_flow(workdir, manager, scf_input, nscf_input):
    """
    Build an `AbinitFlow` for band structure calculations.

    Args:
        workdir:
            Working directory.
        manager:
            `TaskManager` object used to submit the jobs
        scf_input:
            Input for the GS SCF run.
        nscf_input:
            Input for the NSCF run (band structure run).

    Returns:
        `AbinitFlow`
    """
    flow = AbinitFlow(workdir, manager)
    work = BandStructureWorkflow(scf_input, nscf_input)
    flow.register_work(work)
    return flow.allocate()


def g0w0_flow(workdir, manager, scf_input, nscf_input, scr_input, sigma_input):
    """
    Build an `AbinitFlow` for one-shot G0W0 calculations.

    Args:
        workdir:
            Working directory.
        manager:
            `TaskManager` object used to submit the jobs
        scf_input:
            Input for the GS SCF run.
        nscf_input:
            Input for the NSCF run (band structure run).
        scr_input:
            Input for the SCR run.
        sigma_input:
            Input for the SIGMA run.

    Returns:
        `AbinitFlow`
    """
    flow = AbinitFlow(workdir, manager)
    work = G0W0_Workflow(scf_input, scf_input, nscf_input, scr_input, sigma_input)
    flow.register_work(work)
    return flow.allocate()
