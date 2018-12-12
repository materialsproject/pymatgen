# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Works for Abinit
"""

import os
import shutil
import time
import abc
import collections
import numpy as np
import copy

from monty.collections import AttrDict
from monty.itertools import chunks
from monty.functools import lazy_property
from monty.fnmatch import WildCard
from pydispatch import dispatcher
from pymatgen.core.units import EnergyArray
from . import wrappers
from .nodes import Dependency, Node, NodeError, NodeResults, FileNode, check_spectator
from .tasks import (Task, AbinitTask, ScfTask, NscfTask, DfptTask, PhononTask, ElasticTask, DdkTask,
                    BseTask, RelaxTask, DdeTask, BecTask, ScrTask, SigmaTask, TaskManager,
                    DteTask, EphTask, CollinearThenNonCollinearScfTask)

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
    "PhononWfkqWork",
    "GKKPWork",
    "BecWork",
    "DteWork",
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


class BaseWork(Node, metaclass=abc.ABCMeta):
    Error = WorkError

    Results = WorkResults

    # interface modeled after subprocess.Popen
    @property
    @abc.abstractmethod
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
        It executes on_all_ok when all tasks in self have reached `S_OK`.
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
        This method is called once the `Work` is completed i.e. when all tasks
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

    def get_graphviz(self, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate task graph in the DOT language (only parents and children of this work).

        Args:
            engine: Layout command used. ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']
            graph_attr: Mapping of (attribute, value) pairs for the graph.
            node_attr: Mapping of (attribute, value) pairs set for all nodes.
            edge_attr: Mapping of (attribute, value) pairs set for all edges.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        from graphviz import Digraph
        fg = Digraph("work", #filename="work_%s.gv" % os.path.basename(self.workdir),
            engine="fdp" if engine == "automatic" else engine)

        # Set graph attributes.
        # https://www.graphviz.org/doc/info/
        #fg.attr(label="%s@%s" % (self.__class__.__name__, self.relworkdir))
        fg.attr(label=repr(self))
        #fg.attr(fontcolor="white", bgcolor='purple:pink')
        fg.attr(rankdir="LR", pagedir="BL")
        #fg.attr(constraint="false", pack="true", packMode="clust")
        fg.node_attr.update(color='lightblue2', style='filled')
        #fg.node_attr.update(ranksep='equally')

        # Add input attributes.
        if graph_attr is not None:
            fg.graph_attr.update(**graph_attr)
        if node_attr is not None:
            fg.node_attr.update(**node_attr)
        if edge_attr is not None:
            fg.edge_attr.update(**edge_attr)

        def node_kwargs(node):
            return dict(
                #shape="circle",
                color=node.color_hex,
                label=(str(node) if not hasattr(node, "pos_str") else
                    node.pos_str + "\n" + node.__class__.__name__),
            )

        edge_kwargs = dict(arrowType="vee", style="solid")
        cluster_kwargs = dict(rankdir="LR", pagedir="BL", style="rounded", bgcolor="azure2")

        # Build cluster with tasks in *this* work
        cluster_name = "cluster%s" % self.name
        with fg.subgraph(name=cluster_name) as wg:
            wg.attr(**cluster_kwargs)
            wg.attr(label="%s (%s)" % (self.__class__.__name__, self.name))
            for task in self:
                wg.node(task.name, **node_kwargs(task))
                # Connect task to children
                for child in task.get_children():
                    # Test if child is in this cluster (self).
                    myg = wg if child in self else fg
                    myg.node(child.name, **node_kwargs(child))
                    # Find file extensions required by this task
                    i = [dep.node for dep in child.deps].index(task)
                    edge_label = "+".join(child.deps[i].exts)
                    myg.edge(task.name, child.name, label=edge_label, color=task.color_hex,
                             **edge_kwargs)

                # Connect task to parents
                for parent in task.get_parents():
                    # Test if parent is in this cluster (self).
                    myg = wg if parent in self else fg
                    myg.node(parent.name, **node_kwargs(parent))
                    # Find file extensions required by this task
                    i = [dep.node for dep in task.deps].index(parent)
                    edge_label = "+".join(task.deps[i].exts)
                    myg.edge(parent.name, task.name, label=edge_label, color=parent.color_hex,
                             **edge_kwargs)

        # Treat the case in which we have a work producing output for tasks in *this* work.
        #for work in self.flow:
        #    children = work.get_children()
        #    if not children or all(child not in self for child in children):
        #        continue
        #    cluster_name = "cluster%s" % work.name
        #    seen = set()
        #    for child in children:
        #        if child not in self: continue
        #        # This is not needed, too much confusing
        #        #fg.edge(cluster_name, child.name, color=work.color_hex, **edge_kwargs)
        #        # Find file extensions required by work
        #        i = [dep.node for dep in child.deps].index(work)
        #        for ext in child.deps[i].exts:
        #            out = "%s (%s)" % (ext, work.name)
        #            fg.node(out)
        #            fg.edge(out, child.name, **edge_kwargs)
        #            key = (cluster_name, out)
        #            if key not in seen:
        #                fg.edge(cluster_name, out, color=work.color_hex, **edge_kwargs)
        #                seen.add(key)

        return fg


class NodeContainer(metaclass=abc.ABCMeta):
    """
    Mixin classes for `Work` and `Flow` objects providing helper functions
    to register tasks in the container. The helper function call the
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

    def register_collinear_then_noncollinear_scf_task(self, *args, **kwargs):
        """Register a Scf task that perform a SCF run first with nsppol = 2 and then nspinor = 2"""
        kwargs["task_class"] = CollinearThenNonCollinearScfTask
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

    def register_elastic_task(self, *args, **kwargs):
        """Register an elastic task."""
        kwargs["task_class"] = ElasticTask
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

    def register_dde_task(self, *args, **kwargs):
        """Register a Dde task."""
        kwargs["task_class"] = DdeTask
        return self.register_task(*args, **kwargs)

    def register_dte_task(self, *args, **kwargs):
        """Register a Dte task."""
        kwargs["task_class"] = DteTask
        return self.register_task(*args, **kwargs)

    def register_bec_task(self, *args, **kwargs):
        """Register a BEC task."""
        kwargs["task_class"] = BecTask
        return self.register_task(*args, **kwargs)

    def register_bse_task(self, *args, **kwargs):
        """Register a Bethe-Salpeter task."""
        kwargs["task_class"] = BseTask
        return self.register_task(*args, **kwargs)

    def register_eph_task(self, *args, **kwargs):
        """Register an electron-phonon task."""
        kwargs["task_class"] = EphTask
        return self.register_task(*args, **kwargs)

    def walknset_vars(self, task_class=None, *args, **kwargs):
        """
        Set the values of the ABINIT variables in the input files of the nodes

        Args:
            task_class: If not None, only the input files of the tasks belonging
                to class `task_class` are modified.

        Example:

            flow.walknset_vars(ecut=10, kptopt=4)
        """
        def change_task(task):
            if task_class is not None and task.__class__ is not task_class: return False
            return True

        if self.is_work:
            for task in self:
                if not change_task(task): continue
                task.set_vars(*args, **kwargs)

        elif self.is_flow:
            for task in self.iflat_tasks():
                if not change_task(task): continue
                task.set_vars(*args, **kwargs)

        else:
            raise TypeError("Don't know how to set variables for object class %s"  % self.__class__.__name__)


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
        self.wdir = Directory(self.workdir)

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
                # Use the one provided in input else the one of the work/flow.
                if manager is not None:
                    task.set_manager(manager)
                else:
                    # Look first in work and then in the flow.
                    if hasattr(self, "manager"):
                        task.set_manager(self.manager)
                    else:
                        task.set_manager(self.flow.manager)

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
            obj: :class:`AbinitInput` instance or `Task` object.
            deps: Dictionary specifying the dependency of this node or list of dependencies
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

            task = task_class.from_input(obj, task_workdir, manager)

        self._tasks.append(task)

        # Handle possible dependencies given either as dict or list.
        if deps is not None:
            if hasattr(deps, "items"):
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
        # Ignore OK and LOCKED tasks.
        for task in self:
            if task.status in (task.S_OK, task.S_LOCKED): continue
            task.check_status()

        # Take into account possible dependencies. Use a list instead of generators
        for task in self:
            if task.status == task.S_LOCKED: continue
            if task.status < task.S_SUB and all(status == task.S_OK for status in task.deps_status):
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

        return plotter.combiplot(**kwargs)


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
                self.history.info('Converging dilatmx. Value reduce from {} to {}.'
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
    Work for general G0W0 calculations.
    All input can be either single inputs or lists of inputs
    """
    def __init__(self, scf_inputs, nscf_inputs, scr_inputs, sigma_inputs,
                 workdir=None, manager=None):
        """
        Args:
            scf_inputs: Input(s) for the SCF run, if it is a list add all but only link
                to the last input (used for convergence studies on the KS band gap)
            nscf_inputs: Input(s) for the NSCF run, if it is a list add all but only
                link to the last (i.e. addditiona DOS and BANDS)
            scr_inputs: Input for the screening run
            sigma_inputs: List of :class:AbinitInput`for the self-energy run.
                if scr and sigma are lists of the same length, every sigma gets its own screening.
                if there is only one screening all sigma inputs are linked to this one
            workdir: Working directory of the calculation.
            manager: :class:`TaskManager` object.
        """
        super(G0W0Work, self).__init__(workdir=workdir, manager=manager)

        spread_scr = (isinstance(sigma_inputs, (list, tuple)) and
                      isinstance(scr_inputs, (list, tuple)) and
                      len(sigma_inputs) == len(scr_inputs))
        #print("spread_scr", spread_scr)

        self.sigma_tasks = []

        # Register the GS-SCF run.
        # register all scf_inputs but link the nscf only the last scf in the list
        # multiple scf_inputs can be provided to perform convergence studies
        if isinstance(scf_inputs, (list, tuple)):
            for scf_input in scf_inputs:
                self.scf_task = self.register_scf_task(scf_input)
        else:
            self.scf_task = self.register_scf_task(scf_inputs)

        # Register the NSCF run (s).
        if isinstance(nscf_inputs, (list, tuple)):
            for nscf_input in nscf_inputs:
                self.nscf_task = nscf_task = self.register_nscf_task(nscf_input, deps={self.scf_task: "DEN"})
        else:
            self.nscf_task = nscf_task = self.register_nscf_task(nscf_inputs, deps={self.scf_task: "DEN"})

        # Register the SCR and SIGMA run(s).
        if spread_scr:
            for scr_input, sigma_input in zip(scr_inputs, sigma_inputs):
                scr_task = self.register_scr_task(scr_input, deps={nscf_task: "WFK"})
                sigma_task = self.register_sigma_task(sigma_input, deps={nscf_task: "WFK", scr_task: "SCR"})
                self.sigma_tasks.append(sigma_task)
        else:
            # Sigma work(s) connected to the same screening.
            scr_task = self.register_scr_task(scr_inputs, deps={nscf_task: "WFK"})
            if isinstance(sigma_inputs, (list, tuple)):
                for inp in sigma_inputs:
                    task = self.register_sigma_task(inp, deps={nscf_task: "WFK", scr_task: "SCR"})
                    self.sigma_tasks.append(task)
            else:
                task = self.register_sigma_task(sigma_inputs, deps={nscf_task: "WFK", scr_task: "SCR"})
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
        fake_task.set_vars({"nqptdm": -1})
        fake_task.start_and_wait()

        # Parse the section with the q-points
        with NetcdfReader(fake_task.outdir.has_abiext("qptdms.nc")) as reader:
            qpoints = reader.read_value("reduced_coordinates_of_kpoints")
        #print("qpoints)

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

        self.history.info("Will call mrgscr to merge %s SCR files:\n" % len(scr_files))
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

# TODO: MergeDdb --> DfptWork(Work) postpone it because it may break pickle.

class MergeDdb:
    """Mixin class for Works that have to merge the DDB files produced by the tasks."""

    def add_becs_from_scf_task(self, scf_task, ddk_tolerance, ph_tolerance):
        """
        Build tasks for the computation of Born effective charges and add them to the work.

        Args:
            scf_task: ScfTask object.
            ddk_tolerance: dict {"varname": value} with the tolerance used in the DDK run.
                None to use AbiPy default.
            ph_tolerance: dict {"varname": value} with the tolerance used in the phonon run.
                None to use AbiPy default.

	Return:
	    (ddk_tasks, bec_tasks)
	"""
        if not isinstance(scf_task, ScfTask):
            raise TypeError("task `%s` does not inherit from ScfTask" % scf_task)

	# DDK calculations (self-consistent to get electric field).
        multi_ddk = scf_task.input.make_ddk_inputs(tolerance=ddk_tolerance)

        ddk_tasks = []
        for ddk_inp in multi_ddk:
            ddk_task = self.register_ddk_task(ddk_inp, deps={scf_task: "WFK"})
            ddk_tasks.append(ddk_task)

        # Build the list of inputs for electric field perturbation and phonons
        # Each BEC task is connected to all the previous DDK task and to the scf_task.
        bec_deps = {ddk_task: "DDK" for ddk_task in ddk_tasks}
        bec_deps.update({scf_task: "WFK"})

        bec_inputs = scf_task.input.make_bec_inputs(tolerance=ph_tolerance)
        bec_tasks = []
        for bec_inp in bec_inputs:
             bec_task = self.register_bec_task(bec_inp, deps=bec_deps)
             bec_tasks.append(bec_task)

        return ddk_tasks, bec_tasks

    def merge_ddb_files(self, delete_source_ddbs=True, only_dfpt_tasks=True,
                        exclude_tasks=None, include_tasks=None):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.

        Args:
            delete_source_ddbs: True if input DDB should be removed once final DDB is created.
            only_dfpt_tasks: False to merge all DDB files produced by the tasks of the work
                Useful e.g. for finite stress corrections in which the stress in the
                initial configuration should be merged in the final DDB.
            exclude_tasks: List of tasks that should be excluded when merging the partial DDB files.
            include_tasks: List of tasks that should be included when merging the partial DDB files.
                Mutually exclusive with exclude_tasks.

        Returns:
            path to the output DDB file
        """
        if exclude_tasks:
            my_tasks = [task for task in self if task not in exclude_tasks]
        elif include_tasks:
            my_tasks = [task for task in self if task in include_tasks]
        else:
            my_tasks = [task for task in self]

        if only_dfpt_tasks:
            ddb_files = list(filter(None, [task.outdir.has_abiext("DDB") for task in my_tasks \
                                       if isinstance(task, DfptTask)]))
        else:
            ddb_files = list(filter(None, [task.outdir.has_abiext("DDB") for task in my_tasks]))

        self.history.info("Will call mrgddb to merge %s DDB files:" % len(ddb_files))
        # DDB files are always produces so this should never happen!
        if not ddb_files:
            raise RuntimeError("Cannot find any DDB file to merge by the task of " % self)

        # Final DDB file will be produced in the outdir of the work.
        out_ddb = self.outdir.path_in("out_DDB")

        if len(ddb_files) == 1:
            # Avoid the merge. Just copy the DDB file to the outdir of the work.
            shutil.copy(ddb_files[0], out_ddb)
        else:
            # Call mrgddb
            desc = "DDB file merged by %s on %s" % (self.__class__.__name__, time.asctime())
            mrgddb = wrappers.Mrgddb(manager=self[0].manager, verbose=0)
            mrgddb.merge(self.outdir.path, ddb_files, out_ddb=out_ddb, description=desc,
                         delete_source_ddbs=delete_source_ddbs)

        return out_ddb

    def merge_pot1_files(self, delete_source=True):
        """
        This method is called when all the q-points have been computed.
        It runs `mrgdvdb` in sequential on the local machine to produce
        the final DVDB file in the outdir of the `Work`.

        Args:
            delete_source: True if POT1 files should be removed after (successful) merge.

        Returns:
            path to the output DVDB file. None if not DFPT POT file is found.
        """
        natom = len(self[0].input.structure)
        max_pertcase = 3 * natom

        pot1_files = []
        for task in self:
            if not isinstance(task, DfptTask): continue
            paths = task.outdir.list_filepaths(wildcard="*_POT*")
            for path in paths:
                # Include only atomic perturbations i.e. files whose ext <= 3 * natom
                i = path.rindex("_POT")
                pertcase = int(path[i+4:].replace(".nc", ""))
                if pertcase <= max_pertcase:
                    pot1_files.append(path)

        # prtpot = 0 disables the output of the DFPT POT files so an empty list is not fatal here.
        if not pot1_files: return None

        self.history.info("Will call mrgdvdb to merge %s files:" % len(pot1_files))

        # Final DDB file will be produced in the outdir of the work.
        out_dvdb = self.outdir.path_in("out_DVDB")

        if len(pot1_files) == 1:
            # Avoid the merge. Just move the DDB file to the outdir of the work
            shutil.copy(pot1_files[0], out_dvdb)
        else:
            # FIXME: The merge may require a non-negligible amount of memory if lots of qpts.
            # Besides there are machines such as lemaitre3 that are problematic when
            # running MPI applications on the front-end
            mrgdvdb = wrappers.Mrgdvdb(manager=self[0].manager, verbose=0)
            mrgdvdb.merge(self.outdir.path, pot1_files, out_dvdb, delete_source=delete_source)

        return out_dvdb


class PhononWork(Work, MergeDdb):
    """
    This work consists of nirred Phonon tasks where nirred is
    the number of irreducible atomic perturbations for a given set of q-points.
    It provides the callback method (on_all_ok) that calls mrgddb (mrgdv) to merge
    all the partial DDB (POT) files produced. The two files are available in the
    output directory of the Work.
    """

    @classmethod
    def from_scf_task(cls, scf_task, qpoints, is_ngqpt=False, tolerance=None, with_becs=False,
                      ddk_tolerance=None, manager=None):
        """
        Construct a `PhononWork` from a :class:`ScfTask` object.
        The input file for phonons is automatically generated from the input of the ScfTask.
        Each phonon task depends on the WFK file produced by the `scf_task`.

        Args:
            scf_task: ScfTask object.
            qpoints: q-points in reduced coordinates. Accepts single q-point, list of q-points
                or three integers defining the q-mesh if `is_ngqpt`.
            is_ngqpt: True if `qpoints` should be interpreted as divisions instead of q-points.
            tolerance: dict {"varname": value} with the tolerance to be used in the phonon run.
                None to use AbiPy default.
            with_becs: Activate calculation of Electric field and Born effective charges.
            ddk_tolerance: dict {"varname": value} with the tolerance used in the DDK run if with_becs.
                None to use AbiPy default.
            manager: :class:`TaskManager` object.
        """
        if not isinstance(scf_task, ScfTask):
            raise TypeError("task `%s` does not inherit from ScfTask" % scf_task)

        if is_ngqpt:
            qpoints = scf_task.input.abiget_ibz(ngkpt=qpoints, shiftk=[0, 0, 0], kptopt=1).points
        qpoints = np.reshape(qpoints, (-1, 3))

        new = cls(manager=manager)
        if with_becs:
            new.add_becs_from_scf_task(scf_task, ddk_tolerance, ph_tolerance=tolerance)

        for qpt in qpoints:
            if with_becs and np.sum(qpt ** 2) < 1e-12: continue
            multi = scf_task.input.make_ph_inputs_qpoint(qpt, tolerance=tolerance)
            for ph_inp in multi:
                new.register_phonon_task(ph_inp, deps={scf_task: "WFK"})

        return new

    @classmethod
    def from_scf_input(cls, scf_input, qpoints, is_ngqpt=False, tolerance=None,
                       with_becs=False, ddk_tolerance=None, manager=None):
        """
        Similar to `from_scf_task`, the difference is that this method requires
        an input for SCF calculation. A new ScfTask is created and added to the Work.
        This API should be used if the DDB of the GS task should be merged.
        """
        if is_ngqpt:
            qpoints = scf_input.abiget_ibz(ngkpt=qpoints, shiftk=[0, 0, 0], kptopt=1).points

        qpoints = np.reshape(qpoints, (-1, 3))

        new = cls(manager=manager)
        # Create ScfTask
        scf_task = new.register_scf_task(scf_input)

        if with_becs:
            new.add_becs_from_scf_task(scf_task, ddk_tolerance, ph_tolerance=tolerance)

        for qpt in qpoints:
            if with_becs and np.sum(qpt ** 2) < 1e-12: continue
            multi = scf_task.input.make_ph_inputs_qpoint(qpt, tolerance=tolerance)
            for ph_inp in multi:
                new.register_phonon_task(ph_inp, deps={scf_task: "WFK"})

        return new

    #@check_spectator
    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.
        """
        # Merge DDB files.
        out_ddb = self.merge_ddb_files()

        # Merge DVDB files.
        out_dvdb = self.merge_pot1_files()

        return self.Results(node=self, returncode=0, message="DDB merge done")


class PhononWfkqWork(Work, MergeDdb):
    """
    This work computes phonons with DFPT on an arbitrary q-mesh (usually denser than the k-mesh for electrons)
    by computing WKQ files for each q-point.
    The number of irreducible atomic perturbations for each q-point are taken into account.
    It provides the callback method (on_all_ok) that calls mrgddb (mrgdv) to merge
    all the partial DDB (POT) files produced. The two files are available in the
    output directory of the Work. The WKQ files are removed at runtime.
    """

    @classmethod
    def from_scf_task(cls, scf_task, ngqpt, ph_tolerance=None, tolwfr=1.0e-22, nband=None,
                      with_becs=False, ddk_tolerance=None, shiftq=(0, 0, 0), is_ngqpt=True, remove_wfkq=True,
                      manager=None):
        """
        Construct a `PhononWfkqWork` from a :class:`ScfTask` object.
        The input files for WFQ and phonons are automatically generated from the input of the ScfTask.
        Each phonon task depends on the WFK file produced by scf_task and the associated WFQ file.

        Args:
            scf_task: ScfTask object.
            ngqpt: three integers defining the q-mesh
            with_becs: Activate calculation of Electric field and Born effective charges.
            ph_tolerance: dict {"varname": value} with the tolerance for the phonon run.
                None to use AbiPy default.
            tolwfr: tolerance used to compute WFQ.
            ddk_tolerance: dict {"varname": value} with the tolerance used in the DDK run if with_becs.
                None to use AbiPy default.
            shiftq: Q-mesh shift. Multiple shifts are not supported.
            is_ngqpt: the ngqpt is interpreted as a set of integers defining the q-mesh, otherwise
                      is an explicit list of q-points
            remove_wfkq: Remove WKQ files when the children are completed.
            manager: :class:`TaskManager` object.

        .. note:

            Use k-meshes with one shift and q-meshes that are multiple of ngkpt
            to decrease the number of WFQ files to be computed.
        """
        if not isinstance(scf_task, ScfTask):
            raise TypeError("task `%s` does not inherit from ScfTask" % scf_task)

        shiftq = np.reshape(shiftq, (3,))
        if is_ngqpt:
            qpoints = scf_task.input.abiget_ibz(ngkpt=ngqpt, shiftk=shiftq, kptopt=1).points
        else:
            qpoints = ngqpt

        new = cls(manager=manager)
        new.remove_wfkq = remove_wfkq
        new.wfkq_tasks = []
        new.wfkq_task_children = collections.defaultdict(list)

        if with_becs:
            # Add DDK and BECS.
            new.add_becs_from_scf_task(scf_task, ddk_tolerance, ph_tolerance)

        # Get ngkpt, shift for electrons from input.
        # Won't try to skip WFQ if multiple shifts or off-diagonal kptrlatt
        ngkpt, shiftk = scf_task.input.get_ngkpt_shiftk()
        try_to_skip_wfkq = True
        if ngkpt is None or len(shiftk) > 1 and is_ngqpt:
            try_to_skip_wfkq = True

        # TODO: One could avoid kptopt 3 by computing WFK in the IBZ and then rotating.
        # but this has to be done inside Abinit.
        for qpt in qpoints:
            is_gamma = np.sum(qpt ** 2) < 1e-12
            if with_becs and is_gamma: continue

            # Avoid WFQ if k + q = k (requires ngkpt, multiple shifts are not supported)
            need_wfkq = True
            if is_gamma:
                need_wfkq = False
            elif try_to_skip_wfkq:
                # k = (i + shiftk) / ngkpt
                qinds = np.rint(qpt * ngqpt - shiftq)
                f = (qinds * ngkpt) % ngqpt
                need_wfkq = np.any(f != 0)

            if need_wfkq:
                nscf_inp = scf_task.input.new_with_vars(qpt=qpt, nqpt=1, iscf=-2, kptopt=3, tolwfr=tolwfr)
                if nband:
                    nbdbuf = max(2,nband*0.1)
                    nscf_inp.set_vars(nband=nband+nbdbuf, nbdbuf=nbdbuf)
                wfkq_task = new.register_nscf_task(nscf_inp, deps={scf_task: ["DEN", "WFK"]})
                new.wfkq_tasks.append(wfkq_task)

            multi = scf_task.input.make_ph_inputs_qpoint(qpt, tolerance=ph_tolerance)
            for ph_inp in multi:
                deps = {scf_task: "WFK", wfkq_task: "WFQ"} if need_wfkq else {scf_task: "WFK"}
                #ph_inp["prtwf"] = -1
                t = new.register_phonon_task(ph_inp, deps=deps)
                if need_wfkq:
                    new.wfkq_task_children[wfkq_task].append(t)

        return new

    def on_ok(self, sender):
        """
        This callback is called when one task reaches status `S_OK`.
        It removes the WFKQ file if all its children have reached `S_OK`.
        """
        if self.remove_wfkq:
            for task in self.wfkq_tasks:
                if task.status != task.S_OK: continue
                children = self.wfkq_task_children[task]
                if all(child.status == child.S_OK for child in children):
                   path = task.outdir.has_abiext("WFQ")
                   if path:
                       self.history.info("Removing WFQ: %s" % path)
                       os.remove(path)

        return super(PhononWfkqWork, self).on_ok(sender)

    #@check_spectator
    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.
        """
        # Merge DDB files.
        out_ddb = self.merge_ddb_files()

        # Merge DVDB files.
        out_dvdb = self.merge_pot1_files()

        return self.Results(node=self, returncode=0, message="DDB merge done")

class GKKPWork(Work):
    """
    This work computes electron-phonon matrix elements for all the q-points
    present in a DVDB and DDB file
    """
    @classmethod
    def from_den_ddb_dvdb(cls, inp, den_path, ddb_path, dvdb_path, mpiprocs=1, remove_wfkq=True,
                          qpath=None, with_ddk=True, expand=True, manager=None):
        """
        Construct a `PhononWfkqWork` from a DDB and DVDB file.
        For each q found, a WFQ task and an EPH task computing the matrix elements are created.
        """
        import abipy.abilab as abilab

        # Create file nodes
        den_file = FileNode(den_path)
        ddb_file = FileNode(ddb_path)
        dvdb_file = FileNode(dvdb_path)

        # Create new work
        new = cls(manager=manager)
        new.remove_wfkq = remove_wfkq
        new.wfkq_tasks = []
        new.wfkq_task_children = collections.defaultdict(list)
        if manager is None: manager = TaskManager.from_user_config()
        tm = manager.new_with_fixed_mpi_omp(mpiprocs, 1)

        # Create a WFK task
        kptopt = 1 if expand else 3
        nscf_inp = inp.new_with_vars(iscf=-2, kptopt=kptopt)
        wfk_task = new.register_nscf_task(nscf_inp, deps={den_file: "DEN"},manager=tm)
        new.wfkq_tasks.append(wfk_task)
        new.wfk_task = wfk_task

        # Read path and regular grid from DDB file
        with abilab.abiopen(ddb_path) as ddb:
            q_frac_coords = np.array([k.frac_coords for k in ddb.qpoints])
            ddb_ngqpt = ddb.guessed_ngqpt

        # If qpath is set, we read the list of q-points to be used to interpolate the DVDB file.
        # The DVDB and DDB file have to correspond to a regular grid.
        dvdb = dvdb_file
        if qpath is None:
            qpath = q_frac_coords
        else:
            interp_inp = inp.new_with_vars(optdriver=7, eph_task=-5, ddb_ngqpt=ddb_ngqpt,
                                           ph_nqpath=len(qpath), ph_qpath=qpath, prtphdos=0)
            dvdb = new.register_eph_task(interp_inp, deps={wfk_task: "WFK", ddb_file: "DDB", dvdb_file: "DVDB"},
                                          manager=tm)

        # Create a WFK expansion task
        if expand:
            fbz_nscf_inp = inp.new_with_vars(optdriver=8)
            fbz_nscf_inp.set_spell_check(False)
            fbz_nscf_inp.set_vars(wfk_task="wfk_fullbz")
            tm_serial = manager.new_with_fixed_mpi_omp(1,1)
            wfk_task = new.register_nscf_task(fbz_nscf_inp, deps={wfk_task: "WFK", den_file: "DEN"},
                                              manager=tm_serial)
            new.wfkq_tasks.append(wfk_task)
            new.wfk_task = wfk_task

        if with_ddk:
            kptopt = 3 if expand else 1
            ddk_inp = inp.new_with_vars(optdriver=8,kptopt=kptopt)
            ddk_inp.set_spell_check(False)
            ddk_inp.set_vars(wfk_task="wfk_ddk")
            ddk_task = new.register_nscf_task(ddk_inp, deps={wfk_task: "WFK", den_file: "DEN"}, manager=tm)
            new.wfkq_tasks.append(ddk_task)

        # For each qpoint
        for qpt in qpath:
            is_gamma = np.sum(qpt ** 2) < 1e-12
            if is_gamma:
                # Create a link from WFK to WFQ on_ok
                wfkq_task = wfk_task
                deps = {wfk_task: ["WFK","WFQ"], ddb_file: "DDB", dvdb: "DVDB" }
            else:
                # Create a WFQ task
                nscf_inp = nscf_inp.new_with_vars(kptopt=3, qpt=qpt, nqpt=1)
                wfkq_task = new.register_nscf_task(nscf_inp, deps={den_file: "DEN"}, manager=tm)
                new.wfkq_tasks.append(wfkq_task)
                deps = {wfk_task: "WFK", wfkq_task: "WFQ", ddb_file: "DDB", dvdb: "DVDB" }

            # Create a EPH task
            eph_inp = inp.new_with_vars(optdriver=7, prtphdos=0, eph_task=-2, kptopt=3,
                                        ddb_ngqpt=[1,1,1], nqpt=1, qpt=qpt)
            t = new.register_eph_task(eph_inp, deps=deps, manager=tm)
            new.wfkq_task_children[wfkq_task].append(t)

        return new

    @classmethod
    def from_phononwfkq_work(cls, phononwfkq_work, nscf_vars={}, remove_wfkq=True, with_ddk=True, manager=None):
        """
        Construct a `GKKPWork` from a `PhononWfkqWork` object.
        The WFQ are the ones used for PhononWfkqWork so in principle have only valence bands
        """
        # Get list of qpoints from the the phonon tasks in this work
        qpoints = []
        qpoints_deps = []
        for task in phononwfkq_work:
            if isinstance(task,PhononTask):
                # Store qpoints
                qpt = task.input.get("qpt", [0,0,0])
                qpoints.append(qpt)
                # Store dependencies
                qpoints_deps.append(task.deps)

        # Create file nodes
        ddb_path  = phononwfkq_work.outdir.has_abiext("DDB")
        dvdb_path = phononwfkq_work.outdir.has_abiext("DVDB")
        ddb_file = FileNode(ddb_path)
        dvdb_file = FileNode(dvdb_path)

        # Get scf_task from first q-point
        for dep in qpoints_deps[0]:
            if isinstance(dep.node,ScfTask) and dep.exts[0] == 'WFK':
                scf_task = dep.node

        # Create new work
        new = cls(manager=manager)
        new.remove_wfkq = remove_wfkq
        new.wfkq_tasks = []
        new.wfk_task = []

        # Add one eph task per qpoint
        for qpt,qpoint_deps in zip(qpoints,qpoints_deps):
            # Create eph task
            eph_input = scf_task.input.new_with_vars(optdriver=7, prtphdos=0, eph_task=-2,
                                                     ddb_ngqpt=[1,1,1], nqpt=1, qpt=qpt)
            deps = {ddb_file: "DDB", dvdb_file: "DVDB" }
            for dep in qpoint_deps:
                deps[dep.node] = dep.exts[0]
            # If no WFQ in deps link the WFK with WFQ extension
            if 'WFQ' not in deps.values():
                inv_deps = dict((v, k) for k, v in deps.items())
                wfk_task = inv_deps['WFK']
                wfk_path = wfk_task.outdir.has_abiext("WFK")
                # Check if netcdf
                filename, extension = os.path.splitext(wfk_path)
                infile = 'out_WFQ' + extension
                wfq_path = os.path.join(os.path.dirname(wfk_path), infile)
                if not os.path.isfile(wfq_path): os.symlink(wfk_path, wfq_path)
                deps[FileNode(wfq_path)] = 'WFQ'
            new.register_eph_task(eph_input, deps=deps)

        return new

    def on_ok(self, sender):
        """
        This callback is called when one task reaches status `S_OK`.
        It removes the WFKQ file if all its children have reached `S_OK`.
        """
        if self.remove_wfkq:
            for task in self.wfkq_tasks:
                if task.status != task.S_OK: continue
                children = self.wfkq_task_children[task]
                if all(child.status == child.S_OK for child in children):
                   path = task.outdir.has_abiext("WFQ")
                   if path:
                       self.history.info("Removing WFQ: %s" % path)
                       os.remove(path)

        # If wfk task we create a link to a wfq file so abinit is happy
        if sender == self.wfk_task:
            wfk_path = self.wfk_task.outdir.has_abiext("WFK")
            # Check if netcdf
            filename, extension = os.path.splitext(wfk_path)
            infile = 'out_WFQ' + extension
            infile = os.path.join(os.path.dirname(wfk_path), infile)
            os.symlink(wfk_path, infile)

        return super(GKKPWork, self).on_ok(sender)


class BecWork(Work, MergeDdb):
    """
    Work for the computation of the Born effective charges.

    This work consists of DDK tasks and phonon + electric field perturbation
    It provides the callback method (on_all_ok) that calls mrgddb to merge the
    partial DDB files produced by the work.
    """

    @classmethod
    def from_scf_task(cls, scf_task, ddk_tolerance=None, ph_tolerance=None, manager=None):
        """
        Build tasks for the computation of Born effective charges from a ground-state task.

        Args:
            scf_task: ScfTask object.
            ddk_tolerance: tolerance used in the DDK run if with_becs. None to use AbiPy default.
            ph_tolerance: dict {"varname": value} with the tolerance used in the phonon run.
                None to use AbiPy default.
            manager: :class:`TaskManager` object.
	"""
        new = cls(manager=manager)
        new.add_becs_from_scf_task(scf_task, ddk_tolerance, ph_tolerance)
        return new

    def on_all_ok(self):
        """
        This method is called when all tasks reach S_OK
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.
        """
        # Merge DDB files.
        out_ddb = self.merge_ddb_files()
        return self.Results(node=self, returncode=0, message="DDB merge done")


class DteWork(Work, MergeDdb):
    """
    Work for the computation of the third derivative of the energy.

    This work consists of DDK tasks and electric field perturbation.
    It provides the callback method (on_all_ok) that calls mrgddb to merge the partial DDB files produced
    """
    @classmethod
    def from_scf_task(cls, scf_task, ddk_tolerance=None, manager=None):
        """
	Build a DteWork from a ground-state task.

        Args:
            scf_task: ScfTask object.
            ddk_tolerance: tolerance used in the DDK run if with_becs. None to use AbiPy default.
            manager: :class:`TaskManager` object.
	"""
        if not isinstance(scf_task, ScfTask):
            raise TypeError("task `%s` does not inherit from ScfTask" % scf_task)

        new = cls(manager=manager)

        # DDK calculations
        multi_ddk = scf_task.input.make_ddk_inputs(tolerance=ddk_tolerance)

        ddk_tasks = []
        for ddk_inp in multi_ddk:
            ddk_task = new.register_ddk_task(ddk_inp, deps={scf_task: "WFK"})
            ddk_tasks.append(ddk_task)

        # Build the list of inputs for electric field perturbation
        # Each task is connected to all the previous DDK, DDE task and to the scf_task.
        multi_dde = scf_task.input.make_dde_inputs(use_symmetries=False)

        # To compute the nonlinear coefficients all the directions of the perturbation
        # have to be taken in consideration
        # DDE calculations
        dde_tasks = []
        dde_deps = {ddk_task: "DDK" for ddk_task in ddk_tasks}
        dde_deps.update({scf_task: "WFK"})
        for dde_inp in multi_dde:
            dde_task = new.register_dde_task(dde_inp, deps=dde_deps)
            dde_tasks.append(dde_task)

        # DTE calculations
        dte_deps = {scf_task: "WFK DEN"}
        dte_deps.update({dde_task: "1WF 1DEN" for dde_task in dde_tasks})

        multi_dte = scf_task.input.make_dte_inputs()
        dte_tasks = []
        for dte_inp in multi_dte:
             dte_task = new.register_dte_task(dte_inp, deps=dte_deps)
             dte_tasks.append(dte_task)

        return new

    def on_all_ok(self):
        """
        This method is called when all tasks reach S_OK
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Work`.
        """
        # Merge DDB files.
        out_ddb = self.merge_ddb_files()
        return self.Results(node=self, returncode=0, message="DDB merge done")
