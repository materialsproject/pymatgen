"""
Abinit Flows
"""
from __future__ import division, print_function

import os
import sys
import time
import collections
import itertools
import cPickle as pickle
#import pickle as pickle

try:
    from pydispatch import dispatcher
except ImportError:
    pass

from pymatgen.util.io_utils import FileLock
from pymatgen.util.string_utils import pprint_table
from pymatgen.io.abinitio.tasks import Dependency, Node, Task, ScfTask, PhononTask 
from pymatgen.io.abinitio.utils import Directory
from pymatgen.io.abinitio.abiinspect import yaml_read_irred_perts
from pymatgen.io.abinitio.workflows import Workflow, BandStructureWorkflow, PhononWorkflow, G0W0_Workflow

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
    "phonon_flow",
]


class AbinitFlow(Node):
    """
    This object is a container of workflows. Its main task is managing the 
    possible inter-depedencies among the workflows and the creation of
    dynamic worfflows that are generates by callbacks registered by the user.

    .. attributes:

        creation_date:
            String with the creation_date

        pickle_protocol: 
            Protocol for Pickle database (default: -1 i.e. latest protocol)
    """
    VERSION = "0.1"

    PICKLE_FNAME = "__AbinitFlow__.pickle"

    def __init__(self, workdir, manager, auto_restart=False, pickle_protocol=-1):
        """
        Args:
            workdir:
                String specifying the directory where the workflows will be produced.
            manager:
                `TaskManager` object responsible for the submission of the jobs.
            auto_restart:
                True if unconverged calculations should be restarted automatically.
            pickle_procol:
                Pickle protocol version used for saving the status of the object.
                -1 denotes the latest version supported by the python interpreter.
        """
        super(AbinitFlow, self).__init__()

        self.workdir = os.path.abspath(workdir)
        self.creation_date = time.asctime()

        self.manager = manager.deepcopy()
        self.auto_restart = auto_restart

        # List of workflows.
        self._works = []

        # List of callbacks that must be executed when the dependencies reach S_OK
        self._callbacks = []

        # Directories with (input|output|temporary) data.
        self.indir = Directory(os.path.join(self.workdir, "indata"))
        self.outdir = Directory(os.path.join(self.workdir, "outdata"))
        self.tmpdir = Directory(os.path.join(self.workdir, "tmpdata"))

        self.pickle_protocol = int(pickle_protocol)

        # Signal slots: a dictionary with the list 
        # of callbacks indexed by node_id and SIGNAL_TYPE.
        # When the node changes its status, it broadcast a signal.
        # The flow is listening to all the nodes of the calculation
        # [node_id][SIGNAL] = list_of_signal_handlers
        #self._sig_slots =  slots = {}
        #for work in self:
        #    slots[work] = {s: [] for s in work.S_ALL}

        #for task in self.iflat_tasks():
        #    slots[task] = {s: [] for s in work.S_ALL}

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
    def all_ok(self):
        """True if all the tasks in workflows have reached S_OK."""
        return all(work.all_ok for work in self)

    @property
    def num_tasks_with_error(self):
        """The number of tasks whose status is `S_ERROR`."""
        return len(list(self.iflat_tasks(status=self.S_ERROR, op="=")))

    #@property
    #def completed(self):
    #    """True if all the tasks of the flow have reached S_OK."""
    #    return all(task.status == task.S_OK for task in self.iflat_tasks())

    def iflat_tasks_wti(self, status=None, op="="):
        """
        Returns:
            (task, work_index, task_index)
        """
        return self._iflat_tasks_wti(status=status, op=op, with_wti=True)

    def iflat_tasks(self, status=None, op="="):
        """
        Returns:
            task
        """
        return self._iflat_tasks_wti(status=status, op=op, with_wti=False)

    def _iflat_tasks_wti(self, status=None, op="=", with_wti=True):
        """
        Generators that produces a flat sequence of task.
        if status is not None, only the tasks with the specified status are selected.

        Returns:
            (task, work_index, task_index) if with_wti is True else task
        """
        if status is None:
            for wi, work in enumerate(self):
                for ti, task in enumerate(work):
                    if with_wti:
                        yield task, wi, ti
                    else:
                        yield task

        else:
            # Get the operator from the string.
            import operator
            op = {
                "=": operator.eq,
                "!=": operator.ne,
                ">": operator.gt,
                ">=": operator.ge,
                "<": operator.lt,
                "<=": operator.le,
            }[op]

            for wi, work in enumerate(self):
                for ti, task in enumerate(work):
                    if op(task.status, status):
                        if with_wti:
                            yield task, wi, ti
                        else:
                            yield task

    @property
    def status_counter(self):
        """
        Returns a `Counter` object that counts the number of task with 
        given status (use the string representation of the status as key).
        """
        # Count the number of tasks with given status in each workflow.
        counter = self[0].status_counter
        for work in self[1:]:
            counter += work.status_counter

        return counter

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

    def check_status(self):
        """Check the status of the workflows in self."""
        for work in self:
            work.check_status()

        # Test whether some task should be restarted.
        if self.auto_restart:
            num_restarts = 0
            for task, wt in self.iflat_tasks(status=Task.S_UNCONVERGED):
                msg = "AbinitFlow will try restart task %s" % task
                print(msg)
                logger.info(msg)
                retcode = task.restart_if_needed()
                if retcode == 0: 
                    num_restarts += 1

            if num_restarts:
                print("num_restarts done successfully: ", num_restarts)
                self.pickle_dump()

    def show_status(self, stream=sys.stdout):
        """
        Report the status of the workflows and the status 
        of the different tasks on the specified stream.
        """
        for i, work in enumerate(self):
            print(80*"=")
            print("Workflow #%d: %s, Finalized=%s\n" % (i, work, work.finalized) )

            table = [[
                     "Task", "Status", "Queue_id", 
                     "Errors", "Warnings", "Comments", 
                     "MPI", "OMP", 
                     "num_restarts", "max_restarts", "Task Class"
                     ]]

            for task in work:
                task_name = os.path.basename(task.name)

                # Parse the events in the main output.
                report = task.get_event_report()

                events = map(str, 3*["N/A"])
                if report is not None: 
                    events = map(str, [report.num_errors, report.num_warnings, report.num_comments])

                cpu_info = map(str, [task.mpi_ncpus, task.omp_ncpus])
                task_info = map(str, [task.num_restarts, task.max_num_restarts, task.__class__.__name__])

                table.append(
                    [task_name, str(task.status), str(task.queue_id)] + 
                    events + 
                    cpu_info + 
                    task_info
                    )

            pprint_table(table, out=stream)

    def build(self, *args, **kwargs):
        """Make directories and files of the `Flow`."""
        self.indir.makedirs()
        self.outdir.makedirs()
        self.tmpdir.makedirs()

        for work in self:
            work.build(*args, **kwargs)

    def build_and_pickle_dump(self):
        """
        Build dirs and file of the `Flow` and save the object in pickle format.

        Returns:
            0 if success
        """
        self.build()
        return self.pickle_dump()

    def pickle_dump(self):
        """
        Save the status of the object in pickle format.

        Returns:
            0 if success
        """
        protocol = self.pickle_protocol
        filepath = os.path.join(self.workdir, self.PICKLE_FNAME)

        with FileLock(filepath) as lock:
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
        return 0

    @classmethod
    def pickle_load(cls, filepath, disable_signals=False):
        """
        Loads the object from a pickle file and performs initial setup.

        Args:
            filepath:
                Filename or directory name. It filepath is a directory, we 
                scan the directory tree starting from filepath and we 
                read the first pickle database.
            disable_signals:
                If True, the nodes of the flow are not connected by signals.
                This option is usually used when we want to read a flow 
                in read-only mode and we want to avoid any possible side effect.
        """
        if os.path.isdir(filepath):
            # Walk through each directory inside path and find the pickle database.
            for dirpath, dirnames, filenames in os.walk(filepath):
                fnames = [f for f in filenames if f == cls.PICKLE_FNAME]
                if fnames:
                    assert len(fnames) == 1
                    filepath = os.path.join(dirpath, fnames[0])
                    break
            else:
                err_msg = "Cannot find %s inside directory %s" % (cls.PICKLE_FNAME, path)
                raise ValueError(err_msg)

        with FileLock(filepath) as lock:
            with open(filepath, "rb") as fh:
                flow = pickle.load(fh)

        if not disable_signals:
            flow.connect_signals()

        # Recompute the status of each task since tasks that
        # have been submitted previously might be completed.
        flow.check_status()
        return flow

    def register_task(self, input, deps=None, manager=None, task_class=None):
        """
        Utility function that generates a `Workflow` made of a single task

        Args:
            input:
                Abinit Input file or `Strategy` object of `Task` object.
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

    def register_work(self, work, deps=None, manager=None, workdir=None):
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
            workdir:
                The name of the directory used for the `Workflow`.

        Returns:   
            The registered `Workflow`.
        """
        # Directory of the workflow.
        if workdir is None:
            work_workdir = os.path.join(self.workdir, "work_" + str(len(self)))
        else:
            work_workdir = os.path.join(self.workdir, os.path.basenane(workdir))

        work.set_workdir(work_workdir)

        if manager is not None:
            work.set_manager(manager)

        self.works.append(work)

        if deps:
            deps = [Dependency(node, exts) for node, exts in deps.items()]
            work.add_deps(deps)

        return work

    def register_cbk(self, cbk, cbk_data, deps, work_class, manager=None):
        """
        Registers a callback function that will generate the `Task` of the `Workflow`.

        Args:
            cbk:
                Callback function.
            cbk_data
                Additional data passed to the callback function.
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

    def allocate(self, manager=None):
        """
        Allocate the `AbinitFlow` i.e. assign the `workdir` and (optionally) 
        the `TaskManager` to the different tasks in the Flow.
        """
        for work in self:
            work.allocate(manager=self.manager)
            work.set_flow(self)

        for task in self.iflat_tasks():
            task.set_flow(self)

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

        # Associate to each signal the callback _on_signal
        # (bound method of the node that will be called by `AbinitFlow`
        # Each node will set its attribute _done_signal to True to tell
        # the flow that this callback should be disabled.

        # Register the callbacks for the Workflows.
        #for work in self:
        #    slot = self._sig_slots[work]
        #    for signal in S_ALL:
        #        done_signal = getattr(work, "_done_ " + signal, False)
        #        if not done_sig:
        #            cbk_name = "_on_" + str(signal)
        #            cbk = getattr(work, cbk_name, None)
        #            if cbk is None: continue
        #            slot[work][signal].append(cbk)
        #            print("connecting %s\nwith sender %s, signal %s" % (str(cbk), dep.node, dep.node.S_OK))
        #            dispatcher.connect(self.on_dep_ok, signal=signal, sender=dep.node, weak=False)

        # Register the callbacks for the Tasks.

        #self.show_receivers()

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
        self.func = func
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
def bandstructure_flow(workdir, manager, scf_input, nscf_input, dos_inputs=None):
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
        dos_inputs:
            Input(s) for the NSCF run (dos run).

    Returns:
        `AbinitFlow`
    """
    flow = AbinitFlow(workdir, manager)
    work = BandStructureWorkflow(scf_input, nscf_input, dos_inputs=dos_inputs)
    flow.register_work(work)
    return flow.allocate()


def g0w0_flow(workdir, manager, scf_input, nscf_input, scr_input, sigma_inputs):
    """
    Build an `AbinitFlow` for one-shot $G_0W_0$ calculations.

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
        sigma_inputs:
            List of inputs for the SIGMA run.

    Returns:
        `AbinitFlow`
    """
    flow = AbinitFlow(workdir, manager)
    work = G0W0_Workflow(scf_input, nscf_input, scr_input, sigma_inputs)
    flow.register_work(work)
    return flow.allocate()


def phonon_flow(workdir, manager, scf_input, ph_inputs):
    """
    Build an `AbinitFlow` for phonon calculations.

    Args:
        workdir:
            Working directory.
        manager:
            `TaskManager` used to submit the jobs
        scf_input:
            Input for the GS SCF run.
        ph_inputs:
            List of Inputs for the phonon runs.

    Returns:
        `AbinitFlow`
    """
    natom = len(scf_input.structure)

    # Create the container that will manage the different workflows.
    flow = AbinitFlow(workdir, manager)

    # Register the first workflow (GS calculation)
    scf_task = flow.register_task(scf_input, task_class=ScfTask)

    # Build a temporary workflow with a shell manager just to run 
    # ABINIT to get the list of irreducible pertubations for this q-point.
    shell_manager = manager.to_shell_manager(mpi_ncpus=1)
    print(shell_manager)

    if not isinstance(ph_inputs, (list, tuple)):
        ph_inputs = [ph_inputs]

    for i, ph_input in enumerate(ph_inputs):
        fake_input = ph_input.deepcopy()

        # Run abinit on the front-end to get the list of irreducible pertubations.
        tmp_dir = os.path.join(workdir, "__ph_run" + str(i) + "__")
        w = Workflow(workdir=tmp_dir, manager=shell_manager)
        fake_task = w.register(fake_input)
        #print(w)
        #print("manager",w.manager)

        # Use the magic value paral_rf = -1 
        # to get the list of irreducible perturbations for this q-point.
        vars = dict(paral_rf=-1,
                    rfatpol=[1, natom],  # Set of atoms to displace.
                    rfdir=[1, 1, 1],     # Along this set of reduced coordinate axis.
                   )

        fake_task.strategy.add_extra_abivars(vars)
        w.allocate()
        w.start(wait=True)

        # Parse the file to get the perturbations.
        irred_perts = yaml_read_irred_perts(fake_task.log_file.path)
        print(irred_perts)
        w.rmtree()

        # Now we can build the final list of workflows:
        # One workflow per q-point, each workflow computes all 
        # the irreducible perturbations for a singe q-point.
        work_qpt = PhononWorkflow()

        for irred_pert in irred_perts:
            print(irred_pert)
            new_input = ph_input.deepcopy()

            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            qpt = irred_pert["qpt"]
            idir = irred_pert["idir"]
            ipert = irred_pert["ipert"]

            # TODO this will work for phonons, but not for the other types of perturbations.
            rfdir = 3 * [0]
            rfdir[idir -1] = 1
            rfatpol = [ipert, ipert]

            new_input.set_variables(
                #rfpert=1,
                qpt=qpt,
                rfdir=rfdir,
                rfatpol=rfatpol,
            )

            work_qpt.register(new_input, deps={scf_task: "WFK"}, task_class=PhononTask)

        flow.register_work(work_qpt)
                                            
    return flow.allocate()
