# coding: utf-8
"""Abinit Workflows"""
from __future__ import unicode_literals, division, print_function

import os
import shutil
import time
import abc
import collections
import numpy as np
import six

from six.moves import filter
from prettytable import PrettyTable
from monty.collections import AttrDict
from monty.itertools import chunks
from monty.pprint import pprint_table
from monty.functools import lazy_property
from pymatgen.core.units import ArrayWithUnit
from pymatgen.serializers.json_coders import PMGSONable, json_pretty_dump
from pymatgen.util.string_utils import WildCard
from . import wrappers
from .tasks import (Task, AbinitTask, Dependency, Node, NodeResults, ScfTask, NscfTask, PhononTask, DdkTask, BseTask, RelaxTask)
from .strategies import HtcStrategy # ScfStrategy, RelaxStrategy
from .utils import Directory
from .netcdf import ETSF_Reader
from .abitimer import AbinitTimerParser
from .abiinspect import yaml_read_kpoints

try:
    from pydispatch import dispatcher
except ImportError:
    pass

import logging
logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


__all__ = [
    "Workflow",
    "BandStructureWorkflow",
    "RelaxWorkflow",
    "G0W0_Workflow",
    "QptdmWorkflow",
    "SigmaConvWorkflow",
    "BSEMDF_Workflow",
    "PhononWorkflow",
]


class WorkResults(NodeResults):
    JSON_SCHEMA = NodeResults.JSON_SCHEMA.copy() 

    @classmethod
    def from_node(cls, work):
        """Initialize an instance from a WorkFlow instance."""
        new = super(WorkResults, cls).from_node(work)

        #new.update(
        #    #input=work.strategy
        #)

        # Will put all files found in outdir in GridFs 
        # Warning: assuming binary files.
        d = {os.path.basename(f): f for f in work.outdir.list_filepaths()}
        new.add_gridfs_files(**d)

        return new


class WorkflowError(Exception):
    """Base class for the exceptions raised by Workflow objects."""


class BaseWorkflow(six.with_metaclass(abc.ABCMeta, Node)):
    Error = WorkflowError

    Results = WorkResults

    # interface modeled after subprocess.Popen
    @abc.abstractproperty
    def processes(self):
        """Return a list of objects that support the subprocess.Popen protocol."""

    def poll(self):
        """
        Check if all child processes have terminated. Set and return
        returncode attribute.
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

    def show_intrawork_deps(self):
        """Show the dependencies within the `Workflow`."""
        table = PrettyTable(["Task #"] + [str(i) for i in range(len(self))])

        for ii, task1 in enumerate(self):
            line = (1 + len(self)) * [""]
            line[0] = str(ii)
            for jj, task2 in enumerate(self):
                if task1.depends_on(task2):
                    line[jj+1] = "^"
            table.add_row(line)

        print(table)

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
        return sum(task.tot_cores for task in self if task.status == task.S_SUB)

    @property
    def ncores_allocated(self):
        """
        Returns the number of CPUs allocated in this moment.
        A core is allocated if it's running a task or if we have
        submitted a task to the queue manager but the job is still pending.
        """
        return sum(task.tot_cores for task in self if task.status in [task.S_SUB, task.S_RUN])

    @property
    def ncores_inuse(self):
        """
        Returns the number of cores used in this moment.
        A core is used if there's a job that is running on it.
        """
        return sum(task.tot_cores for task in self if task.status == task.S_RUN)

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
        #if all(task.is_completed for task in self) return []
        return [task for task in self if task.can_run]

    @abc.abstractmethod
    def setup(self, *args, **kwargs):
        """Method called before submitting the calculations."""

    def _setup(self, *args, **kwargs):
        self.setup(*args, **kwargs)

    def connect_signals(self):
        """
        Connect the signals within the workflow.
        self is responsible for catching the important signals raised from 
        its task and raise new signals when some particular condition occurs.
        """
        for task in self:
            dispatcher.connect(self.on_ok, signal=task.S_OK, sender=task)

    @property
    def all_ok(self):
        return all(task.status == task.S_OK for task in self)

    def on_ok(self, sender):
        """
        This callback is called when one task reaches status S_OK.
        It executes on_all_ok when all task in self have reached S_OK.
        """
        logger.debug("in on_ok with sender %s" % sender)

        if self.all_ok: 
            if self.finalized:
                return AttrDict(returncode=0, message="Workflow has been already finalized")

            else:
                # Set finalized here, because on_all_ok might change it (e.g. Relax + EOS in a single workflow)
                self.finalized = True
                try:
                    results = AttrDict(**self.on_all_ok())
                except:
                    self.finalized = False
                    raise

                # Signal to possible observers that the `Workflow` reached S_OK
                logger.info("Workflow %s is finalized and broadcasts signal S_OK" % str(self))
                logger.info("Workflow %s status = %s" % (str(self), self.status))

                if self._finalized:
                    dispatcher.send(signal=self.S_OK, sender=self)

                return results

        return AttrDict(returncode=1, message="Not all tasks are OK!")

    def on_all_ok(self):
        """
        This method is called once the `workflow` is completed i.e. when all the tasks 
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


class Workflow(BaseWorkflow):
    """
    A Workflow is a list of (possibly connected) tasks.
    """
    Error = WorkflowError

    def __init__(self, workdir=None, manager=None):
        """
        Args:
            workdir:
                Path to the working directory.
            manager:
                `TaskManager` object.
        """
        super(Workflow, self).__init__()

        self._tasks = []

        if workdir is not None:
            self.set_workdir(workdir)

        if manager is not None:
            self.set_manager(manager)

    def set_manager(self, manager):
        """Set the `TaskManager` to use to launch the Task."""
        self.manager = manager.deepcopy()
        for task in self:
            task.set_manager(manager)

    @property
    def flow(self):
        """The flow containing this `Workflow`."""
        return self._flow

    def set_flow(self, flow):
        """Set the flow associated to this `Workflow`."""
        if not hasattr(self, "_flow"):
            self._flow = flow
        else: 
            if self._flow != flow:
                raise ValueError("self._flow != flow")

    @lazy_property
    def pos(self):
        """The position of self in the Flow"""
        for i, work in enumerate(self.flow):
            if self == work: 
                return i
        raise ValueError("Cannot find the position of %s in flow %s" % (self, self.flow))

    def set_workdir(self, workdir, chroot=False):
        """Set the working directory. Cannot be set more than once unless chroot is True"""
        if not chroot and hasattr(self, "workdir") and self.workdir != workdir:
            raise ValueError("self.workdir != workdir: %s, %s" % (self.workdir,  workdir))

        self.workdir = os.path.abspath(workdir)
                                                                       
        # Directories with (input|output|temporary) data.
        # The workflow will use these directories to connect 
        # itself to other workflows and/or to produce new data 
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
        """True if all the `Task` in the `Workflow` are done."""
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
        of the `Workflow`. It sets the manager of each task (if not already done)
        and defines the working directories of the tasks.

        Args:
            manager:
                `TaskManager` object or None
        """
        for i, task in enumerate(self):

            if not hasattr(task, "manager"):
                # Set the manager
                # Use the one provided in input else the one of the workflow.
                task.set_manager(manager) if manager is not None else task.set_manager(self.manager)

            task_workdir = os.path.join(self.workdir, "t" + str(i))

            if not hasattr(task, "workdir"):
                task.set_workdir(task_workdir)
            else:
                if task.workdir != task_workdir:
                    raise ValueError("task.workdir != task_workdir: %s, %s" % (task.workdir, task_workdir))

    def register(self, obj, deps=None, required_files=None, manager=None, task_class=None):
        """
        Registers a new `Task` and add it to the internal list, taking into account possible dependencies.

        Args:
            obj:
                `Strategy` object or `AbinitInput` instance.
                if Strategy object, we create a new `AbinitTask` from the input strategy and add it to the list.
            deps:
                Dictionary specifying the dependency of this node.
                None means that this obj has no dependency.
            required_files:
                List of strings with the path of the files used by the task.
                Note that the files must exist when the task is registered.
                Use the standard approach based on Workflows, Tasks and deps 
                if the files will be produced in the future.
            manager:
                The `TaskManager` responsible for the submission of the task. If manager is None, we use 
                the `TaskManager` specified during the creation of the `Workflow`.
            task_class:
                Task subclass to instantiate. Default: `AbinitTask` 

        Returns:   
            `Task` object
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

            if isinstance(obj, HtcStrategy):
                # Create the new task (note the factory so that we create subclasses easily).
                task = task_class(obj, task_workdir, manager)

            else:
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

    # Helper functions
    def register_scf_task(self, *args, **kwargs):
        """Register a Scf task."""
        kwargs["task_class"] = ScfTask
        return self.register(*args, **kwargs)
                                                    
    def register_nscf_task(self, *args, **kwargs):
        """Register a nscf task."""
        kwargs["task_class"] = NscfTask
        return self.register(*args, **kwargs)
                                                    
    def register_relax_task(self, *args, **kwargs):
        """Register a task for structural optimization."""
        kwargs["task_class"] = RelaxTask
        return self.register(*args, **kwargs)

    def register_phonon_task(self, *args, **kwargs):
        """Register a phonon task."""
        kwargs["task_class"] = PhononTask
        return self.register(*args, **kwargs)

    def register_ddk_task(self, *args, **kwargs):
        """Register a nscf task."""
        kwargs["task_class"] = DdkTask
        return self.register(*args, **kwargs)

    def register_bse_task(self, *args, **kwargs):
        """Register a nscf task."""
        kwargs["task_class"] = BseTask

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
        # Create the directories of the workflow.
        self.indir.makedirs()
        self.outdir.makedirs()
        self.tmpdir.makedirs()

        # Build dirs and files of each task.
        for task in self:
            task.build(*args, **kwargs)

        # Connect signals within the workflow.
        self.connect_signals()

    @property
    def status(self):
        """
        Returns the status of the workflow i.e. the minimum of the status of the tasks.
        """
        return self.get_all_status(only_min=True)

    #def set_status(self, status):

    def get_all_status(self, only_min=False):
        """
        Returns a list with the status of the tasks in self.

        Args:
            only_min:
                If True, the minimum of the status is returned.
        """
        if len(self) == 0:
            # The workflow will be created in the future.
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
            task.check_status()

        # Take into account possible dependencies. Use a list instead of generators 
        for task in self:
            # changed <= to <
            # todo should this not be < ? a task that is already submitted should not be put to ready
            # it does no harm because of the lock file but logically it seems wrong also gives the wrong infromation
            if task.status < task.S_SUB and all([status == task.S_OK for status in task.deps_status]):
                task.set_status(task.S_READY)

    def rmtree(self, exclude_wildcard=""):
        """
        Remove all files and directories in the working directory

        Args:
            exclude_wildcard:
                Optional string with regular expressions separated by `|`.
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

        return ArrayWithUnit(etotals, "Ha").to(unit)

    def parse_timers(self):
        """
        Parse the TIMER section reported in the ABINIT output files.

        Returns:
            `AbinitTimerParser` object
        """
        filenames = list(filter(os.path.exists, [task.output_file.path for task in self]))
                                                                           
        parser = AbinitTimerParser()
        parser.parse(filenames)
                                                                           
        return parser


class BandStructureWorkflow(Workflow):
    """Workflow for band structure calculations."""
    def __init__(self, scf_input, nscf_input, dos_inputs=None, workdir=None, manager=None):
        """
        Args:
            scf_input:
                Input for the SCF run or `SCFStrategy` object.
            nscf_input:
                Input for the NSCF run or `NSCFStrategy` object defining the band structure calculation.
            dos_inputs:
                Input(s) for the DOS. DOS is computed only if dos_inputs is not None.
            workdir:
                Working directory.
            manager:
                `TaskManager` object.
        """
        super(BandStructureWorkflow, self).__init__(workdir=workdir, manager=manager)

        # Register the GS-SCF run.
        self.scf_task = self.register(scf_input, task_class=ScfTask)

        # Register the NSCF run and its dependency.
        self.nscf_task = self.register(nscf_input, deps={self.scf_task: "DEN"}, task_class=NscfTask)

        # Add DOS computation(s) if requested.
        if dos_inputs is not None:
            if not isinstance(dos_inputs, (list, tuple)):
                dos_inputs = [dos_inputs]

            for dos_input in dos_inputs:
                self.register(dos_input, deps={self.scf_task: "DEN"}, task_class=NscfTask)


class RelaxWorkflow(Workflow):
    """
    Workflow for structural relaxations. The first task relaxes the atomic position
    while keeping the unit cell parameters fixed. The second task uses the final 
    structure to perform a structural relaxation in which both the atomic positions
    and the lattice parameters are optimized.
    """
    def __init__(self, ion_input, ioncell_input, workdir=None, manager=None):
        """
        Args:
            ion_input:
                Input for the relaxation of the ions (cell is fixed)
            ioncell_input:
                Input for the relaxation of the ions and the unit cell.
            workdir:
                Working directory.
            manager:
                `TaskManager` object.
        """
        super(RelaxWorkflow, self).__init__(workdir=workdir, manager=manager)

        self.ion_task = self.register(ion_input, task_class=RelaxTask)

        # Use WFK for the time being since I don't know why Abinit produces all these _TIM?_DEN files.
        #self.ioncell_task = self.register(ioncell_input, deps={self.ion_task: "DEN"}, task_class=RelaxTask)
        self.ioncell_task = self.register(ioncell_input, deps={self.ion_task: "WFK"}, task_class=RelaxTask)

        # Lock ioncell_task as ion_task should communicate to ioncell_task that 
        # the calculation is OK and pass the final structure.
        self.ioncell_task.set_status(self.S_LOCKED)

        self.transfer_done = False

    def on_ok(self, sender):
        """
        This callback is called when one task reaches status S_OK.
        If sender == self.ion_task, we update the initial structure
        used by self.ioncell_task and we unlock it so that the job can be submitted.
        """
        logger.debug("in on_ok with sender %s" % sender)

        if sender == self.ion_task and not self.transfer_done:
            # Get the relaxed structure from ion_task
            ion_structure = self.ion_task.read_final_structure()
            print("Got relaxed ion_structure", ion_structure)

            # Transfer it to the ioncell task (we do it only once).
            self.ioncell_task.change_structure(ion_structure)
            self.transfer_done = True

            # Unlock ioncell_task so that we can submit it.
            self.ioncell_task.set_status(self.S_READY)

        return super(RelaxWorkflow, self).on_ok(sender)


class G0W0_Workflow(Workflow):
    """
    Workflow for G0W0 calculations.
    """
    def __init__(self, scf_input, nscf_input, scr_input, sigma_inputs,
                 workdir=None, manager=None):
        """
        Args:
            scf_input:
                Input for the SCF run or `SCFStrategy` object.
            nscf_input:
                Input for the NSCF run or `NSCFStrategy` object.
            scr_input:
                Input for the screening run or `ScrStrategy` object 
            sigma_inputs:
                List of Strategies for the self-energy run.
            workdir:
                Working directory of the calculation.
            manager:
                `TaskManager` object.
        """
        super(G0W0_Workflow, self).__init__(workdir=workdir, manager=manager)

        # Register the GS-SCF run.
        # register all scf_inputs but link the nscf only the last scf in the list
        if isinstance(scf_input, (list, tuple)):
            for single_scf_input in scf_input:
                self.scf_task = self.register(single_scf_input, task_class=ScfTask)
        else:
            self.scf_task = self.register(scf_input, task_class=ScfTask)

        # Construct the input for the NSCF run.
        self.nscf_task = nscf_task = self.register(nscf_input, deps={self.scf_task: "DEN"}, task_class=NscfTask)

        # Register the SCREENING run.
        self.scr_task = scr_task = self.register(scr_input, deps={nscf_task: "WFK"})

        # Register the SIGMA runs.
        if not isinstance(sigma_inputs, (list, tuple)): 
            sigma_inputs = [sigma_inputs]

        self.sigma_tasks = []
        for sigma_input in sigma_inputs:
            task = self.register(sigma_input, deps={nscf_task: "WFK", scr_task: "SCR"})
            self.sigma_tasks.append(task)


class SigmaConvWorkflow(Workflow):
    """
    Workflow for self-energy convergence studies.
    """
    def __init__(self, wfk_node, scr_node, sigma_inputs, workdir=None, manager=None):
        """
        Args:
            wfk_node:
                The node who has produced the WFK file or filepath pointing to the WFK file.
            scr_node:
                The node who has produced the SCR file or filepath pointing to the SCR file.
            sigma_inputs:
                List of Strategies for the self-energy run.
            workdir:
                Working directory of the calculation.
            manager:
                `TaskManager` object.
        """
        # Cast to node instances.
        wfk_node = Node.as_node(wfk_node)
        scr_node = Node.as_node(scr_node)

        super(SigmaConvWorkflow, self).__init__(workdir=workdir, manager=manager)

        # Register the SIGMA runs.
        if not isinstance(sigma_inputs, (list, tuple)): 
            sigma_inputs = [sigma_inputs]

        for sigma_input in sigma_inputs:
            self.register(sigma_input, deps={wfk_node: "WFK", scr_node: "SCR"})


class BSEMDF_Workflow(Workflow):
    """
    Workflow for simple BSE calculations in which the self-energy corrections
    are approximated by the scissors operator and the screening in modeled
    with the model dielectric function.
    """
    def __init__(self, scf_input, nscf_input, bse_inputs, workdir=None, manager=None):
        """
        Args:
            scf_input:
                Input for the SCF run or `ScfStrategy` object.
            nscf_input:
                Input for the NSCF run or `NscfStrategy` object.
            bse_inputs:
                List of Inputs for the BSE run or `BSEStrategy` object.
            workdir:
                Working directory of the calculation.
            manager:
                `TaskManager`.
        """
        super(BSEMDF_Workflow, self).__init__(workdir=workdir, manager=manager)

        # Register the GS-SCF run.
        self.scf_task = self.register(scf_input, task_class=ScfTask)

        # Construct the input for the NSCF run.
        self.nscf_task = self.register(nscf_input, deps={self.scf_task: "DEN"}, task_class=NscfTask)

        # Construct the input(s) for the BSE run.
        if not isinstance(bse_inputs, (list, tuple)):
            bse_inputs = [bse_inputs]

        for bse_input in bse_inputs:
            self.register(bse_input, deps={self.nscf_task: "WFK"}, task_class=BseTask)


class QptdmWorkflow(Workflow):
    """
    This workflow parallelizes the calculation of the q-points of the screening.
    It also provides the callback `on_all_ok` that calls mrgscr to merge
    all the partial screening files produced.
    """
    def create_tasks(self, wfk_file, scr_input):
        """
        Create the SCR tasks and register them in self.

        Args:
            wfk_file:
                Path to the ABINIT WFK file to use for the computation of the screening.
            scr_input:
                Input for the screening calculation.
        """
        assert len(self) == 0
        wfk_file = self.wfk_file = os.path.abspath(wfk_file)

        # Build a temporary workflow in the tmpdir that will use a shell manager
        # to run ABINIT in order to get the list of q-points for the screening.
        shell_manager = self.manager.to_shell_manager(mpi_procs=1)

        w = Workflow(workdir=self.tmpdir.path_join("_qptdm_run"), manager=shell_manager)

        fake_input = scr_input.deepcopy()
        fake_task = w.register(fake_input)
        w.allocate()
        w.build()

        # Create the symbolic link and add the magic value
        # nqpdm = -1 to the input to get the list of q-points.
        fake_task.inlink_file(wfk_file)
        fake_task.strategy.add_extra_abivars({"nqptdm": -1})
        fake_task.start_and_wait()

        # Parse the section with the q-points
        try:
            qpoints = yaml_read_kpoints(fake_task.log_file.path, doc_tag="!Qptdms")
            #print(qpoints)
        finally:
            w.rmtree()

        # Now we can register the task for the different q-points
        for qpoint in qpoints:
            qptdm_input = scr_input.deepcopy()
            qptdm_input.set_variables(nqptdm=1, qptdm=qpoint)

            self.register(qptdm_input, manager=self.manager)

        self.allocate()

    def merge_scrfiles(self, remove_scrfiles=True):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Workflow`.
        If remove_scrfiles is True, the partial SCR files are removed after the merge.
        """
        scr_files = list(filter(None, [task.outdir.has_abiext("SCR") for task in self]))

        logger.debug("will call mrgscr to merge %s:\n" % str(scr_files))
        assert len(scr_files) == len(self)

        # TODO: Propapagate the manager to the wrappers
        mrgscr = wrappers.Mrgscr(verbose=1)
        mrgscr.set_mpi_runner("mpirun")
        final_scr = mrgscr.merge_qpoints(scr_files, out_prefix="out", cwd=self.outdir.path)

        if remove_scrfiles:
            for scr_file in scr_files:
                try:
                    os.remove(scr_file)
                except IOError:
                    pass

        return final_scr

    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgscr` in sequential on the local machine to produce
        the final SCR file in the outdir of the `Workflow`.
        """
        final_scr = self.merge_scrfiles()
        return self.Results(node=self,returncode=0, message="mrgscr done", final_scr=final_scr)


def build_oneshot_phononwork(workdir, manager, scf_input, ph_inputs, work_class=None):
    """Assuming ph_input compute all the perturbations"""
    work_class = OneShotPhononWorkflow if work_class is None else work_class
    work = work_class(workdir=workdir, manager=manager)
    scf_task = work.register_scf_task(scf_input)
    ph_task = work.register_phonon_task(ph_input, deps={scf_task: "WFK"})
    return work


class OneShotPhononWorkflow(Workflow):
    def read_phonons(self):
        ph_task = self[-1]
        ph_task.output_file.path
        # call the parser here
        return phonons

    def on_all_ok(self):
        """
        """
        # Read phonon frequencies from the output file.
        results = super(OneShotPhononWorkflow, self).get_results()
        results.update(phonons=self.read_phonons())
        return results


class PhononWorkflow(Workflow):
    """
    This workflow usually consists of nirred Phonon tasks where nirred is 
    the number of irreducible perturbations for a given q-point.
    It provides the callback method (on_all_ok) that calls mrgddb to merge 
    the partial DDB files and mrgggkk to merge the GKK files.
    """
    def merge_ddb_files(self):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Workflow`.

        Returns:
            path to the output DDB file
        """
        ddb_files = list(filter(None, [task.outdir.has_abiext("DDB") for task in self]))

        logger.debug("will call mrgddb to merge %s:\n" % str(ddb_files))
        assert len(ddb_files) == len(self)

        #if len(ddb_files) == 1:
        # Avoid the merge. Just move the DDB file to the outdir of the workflow

        # Final DDB file will be produced in the outdir of the workflow.
        out_ddb = self.outdir.path_in("out_DDB")
        desc = "DDB file merged by %s on %s" % (self.__class__.__name__, time.asctime())

        # TODO: propagate the taskmanager
        mrgddb = wrappers.Mrgddb(verbose=1)
        mrgddb.set_mpi_runner("mpirun")
        mrgddb.merge(ddb_files, out_ddb=out_ddb, description=desc, cwd=self.outdir.path)
        return out_ddb

    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Workflow`.
        """
        # Merge DDB files.
        out_ddb = self.merge_ddb_files()

        results = self.Results(node=self,returncode=0, message="DDB merge done")
        results.add_gridfs_files(DDB=(out_ddb, "t"))

        return results

