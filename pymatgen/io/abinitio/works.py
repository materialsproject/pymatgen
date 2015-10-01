# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Works for Abinit:
"""
from __future__ import unicode_literals, division, print_function

import os
import shutil
import time
import abc
import collections
import numpy as np
import six
import copy

from six.moves import filter
from monty.collections import AttrDict
from monty.itertools import chunks
from monty.functools import lazy_property
from monty.fnmatch import WildCard
from pydispatch import dispatcher
from pymatgen.core.units import EnergyArray
from . import wrappers
from .nodes import Dependency, Node, NodeError, NodeResults, check_spectator
from .tasks import (Task, AbinitTask, ScfTask, NscfTask, PhononTask, DdkTask, 
                    BseTask, RelaxTask, DdeTask, BecTask, ScrTask, SigmaTask)

from .utils import Directory
from .netcdf import ETSF_Reader, NetcdfReader
from .abitimer import AbinitTimerParser

import logging
logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


__all__ = [
    "Work",
    "BandStructureWork",
    "RelaxWork",
    "G0W0Work",
    "QptdmWork",
    "SigmaConvWork",
    "BseMdfWork",
    "PhononWork",
]


class WorkResults(NodeResults):
    JSON_SCHEMA = NodeResults.JSON_SCHEMA.copy() 

    @classmethod
    def from_node(cls, work):
        """Initialize an instance from a :class:`Work` instance."""
        new = super(WorkResults, cls).from_node(work)

        # Will put all files found in outdir in GridFs 
        # Warning: assuming binary files.
        d = {os.path.basename(f): f for f in work.outdir.list_filepaths()}
        new.register_gridfs_files(**d)

        return new


class WorkError(NodeError):
    """Base class for the exceptions raised by Work objects."""


class BaseWork(six.with_metaclass(abc.ABCMeta, Node)):
    Error = WorkError

    Results = WorkResults

    # interface modeled after subprocess.Popen
    @abc.abstractproperty
    def processes(self):
        """Return a list of objects that support the `subprocess.Popen` protocol."""

    def poll(self):
        """
        Check if all child processes have terminated. Set and return returncode attribute.
        """
        return [task.poll() for task in self]

    def wait(self):
        """
        Wait for child processed to terminate. Set and return returncode attribute.
        """
        return [task.wait() for task in self]

    def communicate(self, input=None):
        """
        Interact with processes: Send data to stdin. Read data from stdout and
        stderr, until end-of-file is reached.
        Wait for process to terminate. The optional input argument should be a
        string to be sent to the child processed, or None, if no data should be
        sent to the children.

        communicate() returns a list of tuples (stdoutdata, stderrdata).
        """
        return [task.communicate(input) for task in self]

    @property
    def returncodes(self):
        """
        The children return codes, set by poll() and wait() (and indirectly by communicate()).
        A None value indicates that the process hasn't terminated yet.
        A negative value -N indicates that the child was terminated by signal N (Unix only).
        """
        return [task.returncode for task in self]

    @property
    def ncores_reserved(self):
        """
        Returns the number of cores reserved in this moment.
        A core is reserved if it's still not running but 
        we have submitted the task to the queue manager.
        """
        return sum(task.manager.num_cores for task in self if task.status == task.S_SUB)

    @property
    def ncores_allocated(self):
        """
        Returns the number of CPUs allocated in this moment.
        A core is allocated if it's running a task or if we have
        submitted a task to the queue manager but the job is still pending.
        """
        return sum(task.manager.num_cores for task in self if task.status in [task.S_SUB, task.S_RUN])

    @property
    def ncores_used(self):
        """
        Returns the number of cores used in this moment.
        A core is used if there's a job that is running on it.
        """
        return sum(task.manager.num_cores for task in self if task.status == task.S_RUN)

    def fetch_task_to_run(self):
        """
        Returns the first task that is ready to run or 
        None if no task can be submitted at present"

        Raises:
            `StopIteration` if all tasks are done.
        """
        # All the tasks are done so raise an exception 
        # that will be handled by the client code.
        if all(task.is_completed for task in self):
            raise StopIteration("All tasks completed.")

        for task in self:
            if task.can_run:
                return task

        # No task found, this usually happens when we have dependencies. 
        # Beware of possible deadlocks here!
        logger.warning("Possible deadlock in fetch_task_to_run!")
        return None

    def fetch_alltasks_to_run(self):
        """
        Returns a list with all the tasks that can be submitted.
        Empty list if not task has been found.
        """
        return [task for task in self if task.can_run]

    @abc.abstractmethod
    def setup(self, *args, **kwargs):
        """Method called before submitting the calculations."""

    def _setup(self, *args, **kwargs):
        self.setup(*args, **kwargs)

    def connect_signals(self):
        """
        Connect the signals within the work.
        The :class:`Work` is responsible for catching the important signals raised from 
        its task and raise new signals when some particular condition occurs.
        """
        for task in self:
            dispatcher.connect(self.on_ok, signal=task.S_OK, sender=task)

    def disconnect_signals(self):
        """
        Disable the signals within the work. This function reverses the process of `connect_signals`
        """
        for task in self:
            try:
                dispatcher.disconnect(self.on_ok, signal=task.S_OK, sender=task)
            except dispatcher.errors.DispatcherKeyError as exc:
                logger.debug(str(exc))

    @property
    def all_ok(self):
        return all(task.status == task.S_OK for task in self)

    #@check_spectator
    def on_ok(self, sender):
        """
        This callback is called when one task reaches status `S_OK`.
        It executes on_all_ok when all task in self have reached `S_OK`.
        """
        logger.debug("in on_ok with sender %s" % sender)

        if self.all_ok: 
            if self.finalized:
                return AttrDict(returncode=0, message="Work has been already finalized")
            else:
                # Set finalized here, because on_all_ok might change it (e.g. Relax + EOS in a single work)
                self.finalized = True
                try:
                    results = AttrDict(**self.on_all_ok())
                except Exception as exc:
                    self.history.critical("on_all_ok raises %s" % str(exc))
                    self.finalized = False
                    raise

                # Signal to possible observers that the `Work` reached S_OK
                self.history.info("Work %s is finalized and broadcasts signal S_OK" % str(self))

                if self._finalized:
                    self.send_signal(self.S_OK)

                return results

        return AttrDict(returncode=1, message="Not all tasks are OK!")

    #@check_spectator
    def on_all_ok(self):
        """
        This method is called once the `Work` is completed i.e. when all the tasks
        have reached status S_OK. Subclasses should provide their own implementation

        Returns:
            Dictionary that must contain at least the following entries:
                returncode:
                    0 on success. 
                message: 
                    a string that should provide a human-readable description of what has been performed.
        """
        return dict(returncode=0, message="Calling on_all_ok of the base class!")

    def get_results(self, **kwargs):
        """
        Method called once the calculations are completed.
        The base version returns a dictionary task_name: TaskResults for each task in self.
        """
        results = self.Results.from_node(self)
        return results


class NodeContainer(six.with_metaclass(abc.ABCMeta)):
    """
    Mixin classes for `Work` and `Flow` objects providing helperf functions
    to register tasks in the container. The helperfunctios call the
    `register` method of the container.
    """
    # TODO: Abstract protocol for containers

    @abc.abstractmethod
    def register_task(self, *args, **kwargs):
        """
        Register a task in the container.
        """
        # TODO: shall flow.register_task return a Task or a Work?

    # Helper functions
    def register_scf_task(self, *args, **kwargs):
        """Register a Scf task."""
        kwargs["task_class"] = ScfTask
        return self.register_task(*args, **kwargs)
                                                    
    def register_nscf_task(self, *args, **kwargs):
        """Register a nscf task."""
        kwargs["task_class"] = NscfTask
        return self.register_task(*args, **kwargs)
                                                    
    def register_relax_task(self, *args, **kwargs):
        """Register a task for structural optimization."""
        kwargs["task_class"] = RelaxTask
        return self.register_task(*args, **kwargs)

    def register_phonon_task(self, *args, **kwargs):
        """Register a phonon task."""
        kwargs["task_class"] = PhononTask
        return self.register_task(*args, **kwargs)

    def register_ddk_task(self, *args, **kwargs):
        """Register a ddk task."""
        kwargs["task_class"] = DdkTask
        return self.register_task(*args, **kwargs)

    def register_scr_task(self, *args, **kwargs):
        """Register a screening task."""
        kwargs["task_class"] = ScrTask
        return self.register_task(*args, **kwargs)

    def register_sigma_task(self, *args, **kwargs):
        """Register a sigma task."""
        kwargs["task_class"] = SigmaTask
        return self.register_task(*args, **kwargs)

    # TODO: Remove
    def register_dde_task(self, *args, **kwargs):
        """Register a Dde task."""
        kwargs["task_class"] = DdeTask
        return self.register_task(*args, **kwargs)

    def register_bec_task(self, *args, **kwargs):
        """Register a BEC task."""
        kwargs["task_class"] = BecTask
        return self.register_task(*args, **kwargs)

    def register_bse_task(self, *args, **kwargs):
        """Register a nscf task."""
        kwargs["task_class"] = BseTask
        return self.register_task(*args, **kwargs)


class Work(BaseWork, NodeContainer):
    """
    A Work is a list of (possibly connected) tasks.
    """
    def __init__(self, workdir=None, manager=None):
        """
        Args:
            workdir: Path to the working directory.
            manager: :class:`TaskManager` object.
        """
        super(Work, self).__init__()

        self._tasks = []

        if workdir is not None:
            self.set_workdir(workdir)

        if manager is not None:
            self.set_manager(manager)

    def set_manager(self, manager):
        """Set the :class:`TaskManager` to use to launch the :class:`Task`."""
        self.manager = manager.deepcopy()
        for task in self:
            task.set_manager(manager)

    @property
    def flow(self):
        """The flow containing this :class:`Work`."""
        return self._flow

    def set_flow(self, flow):
        """Set the flow associated to this :class:`Work`."""
        if not hasattr(self, "_flow"):
            self._flow = flow
        else: 
            if self._flow != flow:
                raise ValueError("self._flow != flow")

    @lazy_property
    def pos(self):
        """The position of self in the :class:`Flow`"""
        for i, work in enumerate(self.flow):
            if self == work: 
                return i
        raise ValueError("Cannot find the position of %s in flow %s" % (self, self.flow))

    @property
    def pos_str(self):
        """String representation of self.pos"""
        return "w" + str(self.pos) 

    def set_workdir(self, workdir, chroot=False):
        """Set the working directory. Cannot be set more than once unless chroot is True"""
        if not chroot and hasattr(self, "workdir") and self.workdir != workdir:
            raise ValueError("self.workdir != workdir: %s, %s" % (self.workdir,  workdir))

        self.workdir = os.path.abspath(workdir)
                                                                       
        # Directories with (input|output|temporary) data.
        # The work will use these directories to connect 
        # itself to other works and/or to produce new data 
        # that will be used by its children.
        self.indir = Directory(os.path.join(self.workdir, "indata"))
        self.outdir = Directory(os.path.join(self.workdir, "outdata"))
        self.tmpdir = Directory(os.path.join(self.workdir, "tmpdata"))

    def chroot(self, new_workdir):
        self.set_workdir(new_workdir, chroot=True)

        for i, task in enumerate(self):
            new_tdir = os.path.join(self.workdir, "t" + str(i))
            task.set_workdir(new_tdir, chroot=True)

    def __len__(self):
        return len(self._tasks)

    def __iter__(self):
        return self._tasks.__iter__()

    def __getitem__(self, slice):
        return self._tasks[slice]

    def chunks(self, chunk_size):
        """Yield successive chunks of tasks of lenght chunk_size."""
        for tasks in chunks(self, chunk_size):
            yield tasks

    def opath_from_ext(self, ext):
        """
        Returns the path of the output file with extension ext.
        Use it when the file does not exist yet.
        """
        return self.indir.path_in("in_" + ext)

    def opath_from_ext(self, ext):
        """
        Returns the path of the output file with extension ext.
        Use it when the file does not exist yet.
        """
        return self.outdir.path_in("out_" + ext)

    @property
    def processes(self):
        return [task.process for task in self]

    @property
    def all_done(self):
        """True if all the :class:`Task` objects in the :class:`Work` are done."""
        return all(task.status >= task.S_DONE for task in self)

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return all(task.isnc for task in self)

    @property
    def ispaw(self):
        """True if PAW calculation."""
        return all(task.ispaw for task in self)

    @property
    def status_counter(self):
        """
        Returns a `Counter` object that counts the number of task with 
        given status (use the string representation of the status as key).
        """
        counter = collections.Counter() 

        for task in self:
            counter[str(task.status)] += 1

        return counter

    def allocate(self, manager=None):
        """
        This function is called once we have completed the initialization 
        of the :class:`Work`. It sets the manager of each task (if not already done)
        and defines the working directories of the tasks.

        Args:
            manager: :class:`TaskManager` object or None
        """
        for i, task in enumerate(self):

            if not hasattr(task, "manager"):
                # Set the manager
                # Use the one provided in input else the one of the work.
                task.set_manager(manager) if manager is not None else task.set_manager(self.manager)

            task_workdir = os.path.join(self.workdir, "t" + str(i))

            if not hasattr(task, "workdir"):
                task.set_workdir(task_workdir)
            else:
                if task.workdir != task_workdir:
                    raise ValueError("task.workdir != task_workdir: %s, %s" % (task.workdir, task_workdir))

    def register(self, obj, deps=None, required_files=None, manager=None, task_class=None):
        """
        Registers a new :class:`Task` and add it to the internal list, taking into account possible dependencies.

        Args:
            obj: :class:`AbinitInput` instance.
            deps: Dictionary specifying the dependency of this node.
                  None means that this obj has no dependency.
            required_files: List of strings with the path of the files used by the task.
                Note that the files must exist when the task is registered.
                Use the standard approach based on Works, Tasks and deps 
                if the files will be produced in the future.
            manager:
                The :class:`TaskManager` responsible for the submission of the task. If manager is None, we use
                the `TaskManager` specified during the creation of the :class:`Work`.
            task_class: Task subclass to instantiate. Default: :class:`AbinitTask`

        Returns:   
            :class:`Task` object
        """
        task_workdir = None
        if hasattr(self, "workdir"):
            task_workdir = os.path.join(self.workdir, "t" + str(len(self)))

        if isinstance(obj, Task):
            task = obj

        else:
            # Set the class
            if task_class is None:
                task_class = AbinitTask

            #from .strategies import HtcStrategy
            #if isinstance(obj, HtcStrategy):
            #    # Create the new task (note the factory so that we create subclasses easily).
            #    raise NotImplementedError("HtcStrategy")
            #    task = task_class(obj, task_workdir, manager)
            #
            #else:
            task = task_class.from_input(obj, task_workdir, manager)

        self._tasks.append(task)

        # Handle possible dependencies.
        if deps is not None:
            deps = [Dependency(node, exts) for node, exts in deps.items()]
            task.add_deps(deps)

        # Handle possible dependencies.
        if required_files is not None:
            task.add_required_files(required_files)

        return task

    # Needed by NodeContainer
    register_task = register

    def path_in_workdir(self, filename):
        """Create the absolute path of filename in the working directory."""
        return os.path.join(self.workdir, filename)

    def setup(self, *args, **kwargs):
        """
        Method called before running the calculations.
        The default implementation is empty.
        """

    def build(self, *args, **kwargs):
        """Creates the top level directory."""
        # Create the directories of the work.
        self.indir.makedirs()
        self.outdir.makedirs()
        self.tmpdir.makedirs()

        # Build dirs and files of each task.
        for task in self:
            task.build(*args, **kwargs)

        # Connect signals within the work.
        self.connect_signals()

    @property
    def status(self):
        """
        Returns the status of the work i.e. the minimum of the status of the tasks.
        """
        return self.get_all_status(only_min=True)

    def get_all_status(self, only_min=False):
        """
        Returns a list with the status of the tasks in self.

        Args:
            only_min: If True, the minimum of the status is returned.
        """
        if len(self) == 0:
            # The work will be created in the future.
            if only_min:
                return self.S_INIT
            else:
                return [self.S_INIT]

        self.check_status()
        status_list = [task.status for task in self]

        if only_min:
            return min(status_list)
        else:
            return status_list

    def check_status(self):
        """Check the status of the tasks."""
        # Recompute the status of the tasks
        for task in self:
            if task.status == task.S_LOCKED: continue
            task.check_status()

        # Take into account possible dependencies. Use a list instead of generators 
        for task in self:
            if task.status == task.S_LOCKED: continue
            if task.status < task.S_SUB and all([status == task.S_OK for status in task.deps_status]):
                task.set_status(task.S_READY, "Status set to Ready")

    def rmtree(self, exclude_wildcard=""):
        """
        Remove all files and directories in the working directory

        Args:
            exclude_wildcard: Optional string with regular expressions separated by `|`.
                Files matching one of the regular expressions will be preserved.
                example: exclude_wildard="*.nc|*.txt" preserves all the files
                whose extension is in ["nc", "txt"].
        """
        if not exclude_wildcard:
            shutil.rmtree(self.workdir)

        else:
            w = WildCard(exclude_wildcard)

            for dirpath, dirnames, filenames in os.walk(self.workdir):
                for fname in filenames:
                    path = os.path.join(dirpath, fname)
                    if not w.match(fname):
                        os.remove(path)

    def rm_indatadir(self):
        """Remove all the indata directories."""
        for task in self:
            task.rm_indatadir()

    def rm_outdatadir(self):
        """Remove all the indata directories."""
        for task in self:
            task.rm_outatadir()

    def rm_tmpdatadir(self):
        """Remove all the tmpdata directories."""
        for task in self:
            task.rm_tmpdatadir()

    def move(self, dest, isabspath=False):
        """
        Recursively move self.workdir to another location. This is similar to the Unix "mv" command.
        The destination path must not already exist. If the destination already exists
        but is not a directory, it may be overwritten depending on os.rename() semantics.

        Be default, dest is located in the parent directory of self.workdir, use isabspath=True
        to specify an absolute path.
        """
        if not isabspath:
            dest = os.path.join(os.path.dirname(self.workdir), dest)

        shutil.move(self.workdir, dest)

    def submit_tasks(self, wait=False):
        """
        Submits the task in self and wait.
        TODO: change name.
        """
        for task in self:
            task.start()
    
        if wait:
            for task in self: task.wait()

    def start(self, *args, **kwargs):
        """
        Start the work. Calls build and _setup first, then submit the tasks.
        Non-blocking call unless wait is set to True
        """
        wait = kwargs.pop("wait", False)

        # Initial setup
        self._setup(*args, **kwargs)

        # Build dirs and files.
        self.build(*args, **kwargs)

        # Submit tasks (does not block)
        self.submit_tasks(wait=wait)

    def read_etotals(self, unit="Ha"):
        """
        Reads the total energy from the GSR file produced by the task.

        Return a numpy array with the total energies in Hartree
        The array element is set to np.inf if an exception is raised while reading the GSR file.
        """
        if not self.all_done:
            raise self.Error("Some task is still in running/submitted state")

        etotals = []
        for task in self:
            # Open the GSR file and read etotal (Hartree)
            gsr_path = task.outdir.has_abiext("GSR")
            etot = np.inf
            if gsr_path:
                with ETSF_Reader(gsr_path) as r:
                    etot = r.read_value("etotal")
                
            etotals.append(etot)

        return EnergyArray(etotals, "Ha").to(unit)

    def parse_timers(self):
        """
        Parse the TIMER section reported in the ABINIT output files.

        Returns:
            :class:`AbinitTimerParser` object
        """
        filenames = list(filter(os.path.exists, [task.output_file.path for task in self]))
                                                                           
        parser = AbinitTimerParser()
        parser.parse(filenames)
                                                                           
        return parser


class BandStructureWork(Work):
    """Work for band structure calculations."""

    def __init__(self, scf_input, nscf_input, dos_inputs=None, workdir=None, manager=None):
        """
        Args:
            scf_input: Input for the SCF run 
            nscf_input: Input for the NSCF run defining the band structure calculation.
            dos_inputs: Input(s) for the DOS. DOS is computed only if dos_inputs is not None.
            workdir: Working directory.
            manager: :class:`TaskManager` object.
        """
        super(BandStructureWork, self).__init__(workdir=workdir, manager=manager)

        # Register the GS-SCF run.
        self.scf_task = self.register_scf_task(scf_input)

        # Register the NSCF run and its dependency.
        self.nscf_task = self.register_nscf_task(nscf_input, deps={self.scf_task: "DEN"})

        # Add DOS computation(s) if requested.
        self.dos_tasks = []
        if dos_inputs is not None:
            if not isinstance(dos_inputs, (list, tuple)):
                dos_inputs = [dos_inputs]

            for dos_input in dos_inputs:
                dos_task = self.register_nscf_task(dos_input, deps={self.scf_task: "DEN"})
                self.dos_tasks.append(dos_task)

    def plot_ebands(self, **kwargs):
        """
        Plot the band structure. kwargs are passed to the plot method of :class:`ElectronBands`.

        Returns:
            `matplotlib` figure
        """
        with self.nscf_task.open_gsr() as gsr: 
            return gsr.ebands.plot(**kwargs)

    def plot_ebands_with_edos(self, dos_pos=0, method="gaussian", step=0.01, width=0.1, **kwargs):
        """
        Plot the band structure and the DOS.

        Args:
            dos_pos: Index of the task from which the DOS should be obtained (note: 0 refers to the first DOS task).
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            kwargs: Keyword arguments passed to `plot_with_edos` method to customize the plot.

        Returns:
            `matplotlib` figure.
        """
        with self.nscf_task.open_gsr() as gsr: 
            gs_ebands = gsr.ebands

        with self.dos_tasks[dos_pos].open_gsr() as gsr: 
            dos_ebands = gsr.ebands

        edos = dos_ebands.get_edos(method=method, step=step, width=width)
        return gs_ebands.plot_with_edos(edos, **kwargs)

    def plot_edoses(self, dos_pos=None, method="gaussian", step=0.01, width=0.1, **kwargs):
        """
        Plot the band structure and the DOS.

        Args:
            dos_pos: Index of the task from which the DOS should be obtained. 
                     None is all DOSes should be displayed. Accepts integer or list of integers.
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            kwargs: Keyword arguments passed to `plot` method to customize the plot.

        Returns:
            `matplotlib` figure.
        """
        if dos_pos is not None and not isinstance(dos_pos, (list, tuple)): dos_pos = [dos_pos]

        from abipy.electrons.ebands import ElectronDosPlotter
        plotter = ElectronDosPlotter()
        for i, task in enumerate(self.dos_tasks):
            if dos_pos is not None and i not in dos_pos: continue
            with task.open_gsr() as gsr:
                edos = gsr.ebands.get_edos(method=method, step=step, width=width)
                ngkpt = task.get_inpvar("ngkpt")
                plotter.add_edos("ngkpt %s" % str(ngkpt), edos)

        return plotter.plot(**kwargs)


class RelaxWork(Work):
    """
    Work for structural relaxations. The first task relaxes the atomic position
    while keeping the unit cell parameters fixed. The second task uses the final 
    structure to perform a structural relaxation in which both the atomic positions
    and the lattice parameters are optimized.
    """
    def __init__(self, ion_input, ioncell_input, workdir=None, manager=None, target_dilatmx=None):
        """
        Args:
            ion_input: Input for the relaxation of the ions (cell is fixed)
            ioncell_input: Input for the relaxation of the ions and the unit cell.
            workdir: Working directory.
            manager: :class:`TaskManager` object.
        """
        super(RelaxWork, self).__init__(workdir=workdir, manager=manager)

        self.ion_task = self.register_relax_task(ion_input)

        # Note:
        #   1) It would be nice to restart from the WFK file but ABINIT crashes due to the
        #      different unit cell parameters if paral_kgb == 1
        #paral_kgb = ion_input[0]["paral_kgb"]
        #if paral_kgb == 1:

        #deps = {self.ion_task: "WFK"}  # --> FIXME: Problem in rwwf
        #deps = {self.ion_task: "DEN"}
        deps = None

        self.ioncell_task = self.register_relax_task(ioncell_input, deps=deps)

        # Lock ioncell_task as ion_task should communicate to ioncell_task that 
        # the calculation is OK and pass the final structure.
        self.ioncell_task.lock(source_node=self)
        self.transfer_done = False

        self.target_dilatmx = target_dilatmx

    #@check_spectator
    def on_ok(self, sender):
        """
        This callback is called when one task reaches status S_OK.
        If sender == self.ion_task, we update the initial structure
        used by self.ioncell_task and we unlock it so that the job can be submitted.
        """
        logger.debug("in on_ok with sender %s" % sender)

        if sender == self.ion_task and not self.transfer_done:
            # Get the relaxed structure from ion_task
            ion_structure = self.ion_task.get_final_structure()

            # Transfer it to the ioncell task (we do it only once).
            self.ioncell_task._change_structure(ion_structure)
            self.transfer_done = True

            # Unlock ioncell_task so that we can submit it.
            self.ioncell_task.unlock(source_node=self)

        elif sender == self.ioncell_task and self.target_dilatmx:
            actual_dilatmx = self.ioncell_task.get_inpvar('dilatmx', 1.)
            if self.target_dilatmx < actual_dilatmx:
                self.ioncell_task.reduce_dilatmx(target=self.target_dilatmx)
                logger.info('Converging dilatmx. Value reduce from {} to {}.'
                            .format(actual_dilatmx, self.ioncell_task.get_inpvar('dilatmx')))
                self.ioncell_task.reset_from_scratch()

        return super(RelaxWork, self).on_ok(sender)

    def plot_ion_relaxation(self, **kwargs):
        """
        Plot the history of the ion-cell relaxation. 
        kwargs are passed to the plot method of :class:`HistFile`

        Return `matplotlib` figure or None if hist file is not found.
        """
        with self.ion_task.open_hist() as hist:
            return hist.plot(**kwargs) if hist else None

    def plot_ioncell_relaxation(self, **kwargs):
        """
        Plot the history of the ion-cell relaxation.
        kwargs are passed to the plot method of :class:`HistFile`

        Return `matplotlib` figure or None if hist file is not found.
        """
        with self.ioncell_task.open_hist() as hist:
            return hist.plot(**kwargs) if hist else None


class G0W0Work(Work):
    """
    Work for G0W0 calculations.
    """
    def __init__(self, scf_input, nscf_input, scr_input, sigma_inputs,
                 workdir=None, manager=None, spread_scr=False, nksmall=None):
        """
        Args:
            scf_input: Input for the SCF run
            nscf_input: Input for the NSCF run
            scr_input: Input for the screening run
            sigma_inputs: List of :class:AbinitInput`for the self-energy run.
            workdir: Working directory of the calculation.
            manager: :class:`TaskManager` object.
            spread_scr: Attach a screening task to every sigma task
                if false only one screening task with the max ecuteps and nbands for all sigma tasks
            nksmall: if not None add a dos and bands calculation to the Work
        """
        super(G0W0Work, self).__init__(workdir=workdir, manager=manager)

        # Register the GS-SCF run.
        # register all scf_inputs but link the nscf only the last scf in the list
        #MG: FIXME Why this?
        if isinstance(scf_input, (list, tuple)):
            for single_scf_input in scf_input:
                self.scf_task = self.register_scf_task(single_scf_input)
        else:
            self.scf_task = self.register_scf_task(scf_input)



        nogw = False

        if nksmall:
            raise NotImplementedError("with nksmall but strategies have been removed")
            # if nksmall add bandstructure and dos calculations as well

            from abiobjects import KSampling
            if nksmall < 0:
                nksmall = -nksmall
                nogw = True
            scf_in = scf_input[-1] if isinstance(scf_input, (list, tuple)) else scf_input
            logger.info('added band structure calculation')
            bands_input = NscfStrategy(scf_strategy=scf_in,
                                       ksampling=KSampling.path_from_structure(ndivsm=nksmall, structure=scf_in.structure),
                                       nscf_nband=scf_in.electrons.nband, ecut=scf_in.ecut, chksymbreak=0, tolwfr=1e-18)
            self.bands_task = self.register_nscf_task(bands_input, deps={self.scf_task: "DEN"})
            # note we don not let abinit print the dos, since this is inconpatible with parakgb
            # the dos will be evaluated later using abipy
            dos_input = NscfStrategy(scf_strategy=scf_in,
                                     ksampling=KSampling.automatic_density(kppa=nksmall**3, structure=scf_in.structure,
                                                                           shifts=(0.0, 0.0, 0.0)),
                                     nscf_nband=scf_in.electrons.nband, ecut=scf_in.ecut, chksymbreak=0)

            self.dos_task = self.register_nscf_task(dos_input, deps={self.scf_task: "DEN"})

            #from abiobjects import KSampling
            #if nksmall < 0:
            #    nksmall = -nksmall
            #    nogw = True
            #scf_in = scf_input[-1] if isinstance(scf_input, (list, tuple)) else scf_input
            #logger.info('added band structure calculation')
            #bands_input = NscfStrategy(scf_strategy=scf_in,
            #                           ksampling=KSampling.path_from_structure(ndivsm=nksmall, structure=scf_in.structure),
            #                           nscf_nband=scf_in.electrons.nband, ecut=scf_in.ecut, chksymbreak=0)

            #self.bands_task = self.register_nscf_task(bands_input, deps={self.scf_task: "DEN"})
            ## note we don not let abinit print the dos, since this is inconpatible with parakgb
            ## the dos will be evaluated later using abipy
            #dos_input = NscfStrategy(scf_strategy=scf_in,
            #                         ksampling=KSampling.automatic_density(kppa=nksmall**3, structure=scf_in.structure,
            #                                                               shifts=(0.0, 0.0, 0.0)),
            #                         nscf_nband=scf_in.electrons.nband, ecut=scf_in.ecut, chksymbreak=0)

            #self.dos_task = self.register_nscf_task(dos_input, deps={self.scf_task: "DEN"})

        # Register the SIGMA runs.
        if not nogw:
            # Construct the input for the NSCF run.
            self.nscf_task = nscf_task = self.register_nscf_task(nscf_input, deps={self.scf_task: "DEN"})

            # Register the SCREENING run.
            if not spread_scr:
                self.scr_task = scr_task = self.register_scr_task(scr_input, deps={nscf_task: "WFK"})
            else:
                self.scr_tasks = []

            if not isinstance(sigma_inputs, (list, tuple)):
                sigma_inputs = [sigma_inputs]

            self.sigma_tasks = []
            for sigma_input in sigma_inputs:
                if spread_scr:
                    new_scr_input = copy.deepcopy(scr_input)
                    new_scr_input.screening.ecuteps = sigma_input.sigma.ecuteps
                    new_scr_input.screening.nband = sigma_input.sigma.nband
                    new_scr_input.electrons.nband = sigma_input.sigma.nband
                    scr_task = self.register_scr_task(new_scr_input, deps={nscf_task: "WFK"})

                task = self.register_sigma_task(sigma_input, deps={nscf_task: "WFK", scr_task: "SCR"})
                self.sigma_tasks.append(task)


class SigmaConvWork(Work):
    """
    Work for self-energy convergence studies.
    """
    def __init__(self, wfk_node, scr_node, sigma_inputs, workdir=None, manager=None):
        """
        Args:
            wfk_node: The node who has produced the WFK file or filepath pointing to the WFK file.
            scr_node: The node who has produced the SCR file or filepath pointing to the SCR file.
            sigma_inputs: List of :class:`AbinitInput` for the self-energy runs.
            workdir: Working directory of the calculation.
            manager: :class:`TaskManager` object.
        """
        # Cast to node instances.
        wfk_node, scr_node = Node.as_node(wfk_node), Node.as_node(scr_node)

        super(SigmaConvWork, self).__init__(workdir=workdir, manager=manager)

        # Register the SIGMA runs.
        if not isinstance(sigma_inputs, (list, tuple)): 
            sigma_inputs = [sigma_inputs]

        for sigma_input in sigma_inputs:
            self.register_sigma_task(sigma_input, deps={wfk_node: "WFK", scr_node: "SCR"})


class BseMdfWork(Work):
    """
    Work for simple BSE calculations in which the self-energy corrections
    are approximated by the scissors operator and the screening is modeled
    with the model dielectric function.
    """
    def __init__(self, scf_input, nscf_input, bse_inputs, workdir=None, manager=None):
        """
        Args:
            scf_input: Input for the SCF run.
            nscf_input: Input for the NSCF run.
            bse_inputs: List of Inputs for the BSE run.
            workdir: Working directory of the calculation.
            manager: :class:`TaskManager`.
        """
        super(BseMdfWork, self).__init__(workdir=workdir, manager=manager)

        # Register the GS-SCF run.
        self.scf_task = self.register_scf_task(scf_input)

        # Construct the input for the NSCF run.
        self.nscf_task = self.register_nscf_task(nscf_input, deps={self.scf_task: "DEN"})

        # Construct the input(s) for the BSE run.
        if not isinstance(bse_inputs, (list, tuple)):
            bse_inputs = [bse_inputs]

        for bse_input in bse_inputs:
            self.register_bse_task(bse_input, deps={self.nscf_task: "WFK"})

    def get_mdf_robot(self):
        """Builds and returns a :class:`MdfRobot` for analyzing the results in the MDF files."""
        from abilab.robots import MdfRobot
        robot = MdfRobot()
        for task in self[2:]:
            mdf_path = task.outdir.has_abiext(robot.EXT)
            if mdf_path:
                robot.add_file(str(task), mdf_path)
        return robot

    #def plot_conv_mdf(self, **kwargs)
    #    with self.get_mdf_robot() as robot:
    #        robot.get_mdf_plooter()
    #    plotter.plot(**kwargs)


class QptdmWork(Work):
    """
    This work parallelizes the calculation of the q-points of the screening.
    It also provides the callback `on_all_ok` that calls mrgscr to merge
    all the partial screening files produced.
    """
    def create_tasks(self, wfk_file, scr_input):
        """
        Create the SCR tasks and register them in self.

        Args:
            wfk_file: Path to the ABINIT WFK file to use for the computation of the screening.
            scr_input: Input for the screening calculation.
        """
        assert len(self) == 0
        wfk_file = self.wfk_file = os.path.abspath(wfk_file)

        # Build a temporary work in the tmpdir that will use a shell manager
        # to run ABINIT in order to get the list of q-points for the screening.
        shell_manager = self.manager.to_shell_manager(mpi_procs=1)

        w = Work(workdir=self.tmpdir.path_join("_qptdm_run"), manager=shell_manager)

        fake_input = scr_input.deepcopy()
        fake_task = w.register(fake_input)
        w.allocate()
        w.build()

        # Create the symbolic link and add the magic value
        # nqpdm = -1 to the input to get the list of q-points.
        fake_task.inlink_file(wfk_file)
        fake_task._set_inpvars({"nqptdm": -1})
        fake_task.start_and_wait()

        # Parse the section with the q-points
        with NetcdfReader(fake_task.outdir.has_abiext("qptdms.nc")) as reader:
            qpoints = reader.read_value("reduced_coordinates_of_kpoints")
        #print("qpoints)
        #w.rmtree()

        # Now we can register the task for the different q-points
        for qpoint in qpoints:
            qptdm_input = scr_input.deepcopy()
            qptdm_input.set_vars(nqptdm=1, qptdm=qpoint)
            new_task = self.register_scr_task(qptdm_input, manager=self.manager)
            # Add the garbage collector.
            if self.flow.gc is not None: 
                new_task.set_gc(self.flow.gc)

        self.allocate()

    def merge_scrfiles(self, remove_scrfiles=True):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Work`.
        If remove_scrfiles is True, the partial SCR files are removed after the merge.
        """
        scr_files = list(filter(None, [task.outdir.has_abiext("SCR") for task in self]))

        logger.debug("will call mrgscr to merge %s:\n" % str(scr_files))
        assert len(scr_files) == len(self)

        mrgscr = wrappers.Mrgscr(manager=self[0].manager, verbose=1)
        final_scr = mrgscr.merge_qpoints(self.outdir.path, scr_files, out_prefix="out")

        if remove_scrfiles:
            for scr_file in scr_files:
                try:
                    os.remove(scr_file)
                except IOError:
                    pass

        return final_scr

    #@check_spectator
    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Work`.
        """
        final_scr = self.merge_scrfiles()
        return self.Results(node=self, returncode=0, message="mrgscr done", final_scr=final_scr)


def build_oneshot_phononwork(scf_input, ph_inputs, workdir=None, manager=None, work_class=None):
    """
    Returns a work for the computation of phonon frequencies
    ph_inputs is a list of input for Phonon calculation in which all the independent perturbations 
    are explicitly computed i.e. 

        * rfdir 1 1 1
        * rfatpol 1 natom

    .. warning::
        This work is mainly used for simple calculations, e.g. convergence studies.
        Use :class:`PhononWork` for better efficiency.
    """
    work_class = OneShotPhononWork if work_class is None else work_class
    work = work_class(workdir=workdir, manager=manager)
    scf_task = work.register_scf_task(scf_input)
    ph_inputs = [ph_inputs] if not isinstance(ph_inputs, (list, tuple)) else ph_inputs

    for phinp in ph_inputs:
        # Check rfdir and rfatpol.
        rfdir = np.array(phinp.get("rfdir", [0, 0, 0]))
        if len(rfdir) != 3 or any(rfdir != (1, 1, 1)):
            raise ValueError("Expecting rfdir == (1, 1, 1), got %s" % rfdir)

        rfatpol = np.array(phinp.get("rfatpol", [1, 1]))
        if len(rfatpol) != 2 or any(rfatpol != (1, len(phinp.structure))):
            raise ValueError("Expecting rfatpol == (1, natom), got %s" % rfatpol)

        # cannot use PhononTaks here because the Task is not able to deal with multiple phonon calculations
        ph_task = work.register(phinp, deps={scf_task: "WFK"})

    return work


class OneShotPhononWork(Work):
    """
    Simple and very inefficient work for the computation of the phonon frequencies
    It consists of a GS task and a DFPT calculations for all the independent perturbations.
    The main advantage is that one has direct access to the phonon frequencies that
    can be computed at the end of the second task without having to call anaddb.

    Use ``build_oneshot_phononwork`` to construct this work from the input files.
    """
    def read_phonons(self):
        """
        Read phonon frequencies from the output file.

        Return:
            List of namedtuples. Each `namedtuple` has the following attributes:
                
                - qpt: ndarray with the q-point in reduced coordinates.
                - freqs: ndarray with 3 x Natom phonon frequencies in meV
        """
        # 
        #   Phonon wavevector (reduced coordinates) :  0.00000  0.00000  0.00000
        #  Phonon energies in Hartree :
        #    1.089934E-04  4.990512E-04  1.239177E-03  1.572715E-03  1.576801E-03
        #    1.579326E-03
        #  Phonon frequencies in cm-1    :
        # -  2.392128E+01  1.095291E+02  2.719679E+02  3.451711E+02  3.460677E+02
        # -  3.466221E+02
        BEGIN = "  Phonon wavevector (reduced coordinates) :"
        END = " Phonon frequencies in cm-1    :"

        ph_tasks, qpts, phfreqs = self[1:], [], []
        for task in ph_tasks:

            # Parse output file.
            with open(task.output_file.path, "r") as fh:
                qpt, inside = None, 0 
                for line in fh:
                    if line.startswith(BEGIN):
                        qpts.append([float(s) for s in line[len(BEGIN):].split()])
                        inside, omegas = 1, []
                    elif line.startswith(END):
                        break
                    elif inside:
                        inside += 1
                        if inside > 2:
                            omegas.extend((float(s) for s in line.split()))
                else:
                    raise ValueError("Cannot find %s in file %s" % (END, task.output_file.path))

                phfreqs.append(omegas)

        # Use namedtuple to store q-point and frequencies in meV
        phonon = collections.namedtuple("phonon", "qpt freqs")
        return [phonon(qpt=qpt, freqs=freqs_meV) for qpt, freqs_meV in zip(qpts, EnergyArray(phfreqs, "Ha").to("meV") )]

    def get_results(self, **kwargs):
        results = super(OneShotPhononWork, self).get_results()
        phonons = self.read_phonons()
        results.update(phonons=phonons)
        return results


class MergeDdb(object):
    """Mixin classes for Works that have to merge the DDB files produced by the tasks."""

    def merge_ddb_files(self):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.

        Returns:
            path to the output DDB file
        """
        ddb_files = list(filter(None, [task.outdir.has_abiext("DDB") for task in self]))

        self.history.info("Will call mrgddb to merge %s:\n" % str(ddb_files))
        # assert len(ddb_files) == len(self)

        #if len(ddb_files) == 1:
        # Avoid the merge. Just move the DDB file to the outdir of the work

        # Final DDB file will be produced in the outdir of the work.
        out_ddb = self.outdir.path_in("out_DDB")
        desc = "DDB file merged by %s on %s" % (self.__class__.__name__, time.asctime())

        mrgddb = wrappers.Mrgddb(manager=self[0].manager, verbose=0)
        mrgddb.merge(self.outdir.path, ddb_files, out_ddb=out_ddb, description=desc)

        return out_ddb


class PhononWork(Work, MergeDdb):
    """
    This work usually consists of nirred Phonon tasks where nirred is 
    the number of irreducible perturbations for a given q-point.
    It provides the callback method (on_all_ok) that calls mrgddb to merge the partial DDB files produced 
    """
    @classmethod
    def from_scf_task(cls, scf_task, qpt, tolerance=None):
        """
        Construct a `PhononWork` from a :class:`ScfTask` object.
        The input file for phonons is automatically generated from the input of the ScfTask.
        Each phonon task depends on the WFK file produced by scf_task.

        Args:
            scf_task: ScfTask object. 
            qpt: q-point for phonons in reduced coordinates.
        """
        if not isinstance(scf_task, ScfTask):
            raise TypeError("task %s does not inherit from ScfTask" % scf_task)

        new = cls() #manager=scf_task.manager)

        multi = scf_task.input.make_ph_inputs_qpoint(qpt, tolerance=tolerance)

        for ph_inp in multi:
            new.register_phonon_task(ph_inp, deps={scf_task: "WFK"})

        return new

    # TODO
    #def compute_phonons(self)
    #    """
    #    Call anaddb to compute the phonon frequencies for this q-point and
    #    store the results in the outdir of the work.
    #    """
    #    #atask = AnaddbTask(anaddb_input, ddb_node,
    #    #         gkk_node=None, md_node=None, ddk_node=None, workdir=None, manager=None)
    #    #atask.start_and_wait()
    #    return phonons

    #@check_spectator 
    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.
        """
        # Merge DDB files.
        out_ddb = self.merge_ddb_files()

        results = self.Results(node=self, returncode=0, message="DDB merge done")
        results.register_gridfs_files(DDB=(out_ddb, "t"))

        return results


class BecWork(Work, MergeDdb):
    """
    Work for the computation of the Born effective charges.

    This work consists of DDK tasks and phonon + electric fiel perturbation
    It provides the callback method (on_all_ok) that calls mrgddb to merge the partial DDB files produced 
    """
    @classmethod
    def from_scf_task(cls, scf_task, ddk_tolerance=None):
        """Build a BecWork from a ground-state task."""
        if not isinstance(scf_task, ScfTask):
            raise TypeError("task %s does not inherit from GsTask" % scf_task)

        new = cls() #manager=scf_task.manager)

        # DDK calculations
        multi_ddk = scf_task.input.make_ddk_inputs(tolerance=ddk_tolerance)

        ddk_tasks = []
        for ddk_inp in multi_ddk:
            ddk_task = new.register_ddk_task(ddk_inp, deps={scf_task: "WFK"})
            ddk_tasks.append(ddk_task)

        # Build the list of inputs for electric field perturbation and phonons
        # Each bec task is connected to all the previous DDK task and to the scf_task.
        bec_deps = {ddk_task: "DDK" for ddk_task in ddk_tasks}
        bec_deps.update({scf_task: "WFK"})

        bec_inputs = scf_task.input.make_bec_inputs() #tolerance=efile
        for bec_inp in bec_inputs:
             new.register_bec_task(bec_inp, deps=bec_deps)

        return new

    def on_all_ok(self):
        """
        This method is called when all the task reach S_OK
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.
        """
        # Merge DDB files.
        out_ddb = self.merge_ddb_files()

        results = self.Results(node=self, returncode=0, message="DDB merge done")
        results.register_gridfs_files(DDB=(out_ddb, "t"))

        return results
