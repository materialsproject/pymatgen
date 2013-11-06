"""
Abinit Workflows
"""
from __future__ import division, print_function

import sys
import os
import shutil
import time
import abc
import collections
import numpy as np

try:
    from pydispatch import dispatcher
except ImportError:
    pass

from pymatgen.core.units import ArrayWithUnit, Ha_to_eV
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.design_patterns import Enum, AttrDict
from pymatgen.serializers.json_coders import MSONable, json_pretty_dump
from pymatgen.io.smartio import read_structure
from pymatgen.util.num_utils import iterator_from_slice, chunks, monotonic
from pymatgen.util.string_utils import list_strings, pprint_table, WildCard
from pymatgen.io.abinitio import wrappers
from pymatgen.io.abinitio.tasks import (Task, AbinitTask, Dependency, Node, ScfTask, NscfTask, HaydockBseTask, RelaxTask)
from pymatgen.io.abinitio.strategies import Strategy
from pymatgen.io.abinitio.utils import File, Directory
from pymatgen.io.abinitio.netcdf import ETSF_Reader
from pymatgen.io.abinitio.abiobjects import Smearing, AbiStructure, KSampling, Electrons
from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.io.abinitio.strategies import ScfStrategy
from pymatgen.io.abinitio.eos import EOS
from pymatgen.io.abinitio.abitimer import AbinitTimerParser

import logging
logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

__all__ = [
    "Workflow",
    "IterativeWorkflow",
    "BandStructureWorkflow",
    "RelaxWorkflow",
    "DeltaFactorWorkflow",
    "G0W0_Workflow",
    "BSEMDF_Workflow",
    "PhononWorkflow",
]


class WorkflowError(Exception):
    """Base class for the exceptions raised by Workflow objects."""


class BaseWorkflow(Node):
    __metaclass__ = abc.ABCMeta

    Error = WorkflowError

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
        table = [["Task #"] + [str(i) for i in range(len(self))]]

        for ii, task1 in enumerate(self):
            line = (1 + len(self)) * [""]
            line[0] = str(ii)
            for jj, task2 in enumerate(self):
                if task1.depends_on(task2):
                    line[jj+1] = "^"

            table.append(line)

        pprint_table(table)

    @property
    def returncodes(self):
        """
        The children return codes, set by poll() and wait() (and indirectly by communicate()).
        A None value indicates that the process hasn't terminated yet.
        A negative value -N indicates that the child was terminated by signal N (Unix only).
        """
        return [task.returncode for task in self]

    @property
    def ncpus_reserved(self):
        """
        Returns the number of CPUs reserved in this moment.
        A CPUS is reserved if it's still not running but 
        we have submitted the task to the queue manager.
        """
        return sum(task.tot_ncpus for task in self if task.status == task.S_SUB)

    @property
    def ncpus_allocated(self):
        """
        Returns the number of CPUs allocated in this moment.
        A CPU is allocated if it's running a task or if we have
        submitted a task to the queue manager but the job is still pending.
        """
        return sum(task.tot_ncpus for task in self if task.status in [task.S_SUB, task.S_RUN])

    @property
    def ncpus_inuse(self):
        """
        Returns the number of CPUs used in this moment.
        A CPU is used if there's a job that is running on it.
        """
        return sum(task.tot_ncpus for task in self if task.status == task.S_RUN)

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
                #print(task, str(task.status), [task.deps_status])
                return task

        # No task found, this usually happens when we have dependencies. 
        # Beware of possible deadlocks here!
        logger.warning("Possible deadlock in fetch_task_to_run!")
        return None

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
        """
        logger.debug("in on_ok with sender %s" % sender)

        if self.all_ok: 
            if self.finalized:
                return AttrDict(returncode=0, message="Workflow has been already finalized")

            else:
                results = AttrDict(**self.on_all_ok())
                self._finalized = True
                # Signal to possible observers that the `Workflow` reached S_OK
                print("Workflow %s is finalized and broadcasts signal S_OK" % str(self))
                print("Workflow %s status = %s" % (str(self), self.status))
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
        return dict(returncode=0, 
                    message="Calling on_all_ok of the base class!",
                    )

    def get_results(self):
        """
        Method called once the calculations are completed.

        The base version returns a dictionary task_name : TaskResults for each task in self.
        """
        return WorkflowResults(task_results={task.name: task.results for task in self})


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

    def set_workdir(self, workdir):
        """Set the working directory. Cannot be set more than once."""

        if hasattr(self, "workdir") and self.workdir != workdir:
            raise ValueError("self.workdir != workdir: %s, %s" % (self.workdir,  workdir))

        self.workdir = os.path.abspath(workdir)
                                                                       
        # Directories with (input|output|temporary) data.
        # The workflow will use these directories to connect 
        # itself to other workflows and/or to produce new data 
        # that will be used by its children.
        self.indir = Directory(os.path.join(self.workdir, "indata"))
        self.outdir = Directory(os.path.join(self.workdir, "outdata"))
        self.tmpdir = Directory(os.path.join(self.workdir, "tmpdata"))

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
            #print(hasattr(task, "manager"))
            if not hasattr(task, "manager"):
                # Set the manager

                if manager is not None:
                    # Use the one provided in input.
                    task.set_manager(manager)
                else:
                    # Use the one of the workflow.
                    task.set_manager(self.manager)

            task_workdir = os.path.join(self.workdir, "task_" + str(i))

            if not hasattr(task, "workdir"):
                task.set_workdir(task_workdir)
            else:
                if task.workdir != task_workdir:
                    raise ValueError("task.workdir != task_workdir: %s, %s" % (task.workdir, task_workdir))

    def register(self, obj, deps=None, manager=None, task_class=None):
        """
        Registers a new `Task` and add it to the internal list, taking into account possible dependencies.

        Args:
            obj:
                `Strategy` object or `AbinitInput` instance.
                if Strategy object, we create a new `AbinitTask` from the input strategy and add it to the list.
            deps:
                Dictionary specifying the dependency of this node.
                None means that this obj has no dependency.
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
            task_workdir = os.path.join(self.workdir, "task_" + str(len(self)))

        if isinstance(obj, Task):
            task = obj

        else:
            # Set the class
            if task_class is None:
                task_class = AbinitTask

            if isinstance(obj, Strategy):
                # Create the new task (note the factory so that we create subclasses easily).
                task = task_class(obj, task_workdir, manager)

            else:
                task = task_class.from_input(obj, task_workdir, manager)

        self._tasks.append(task)

        # Handle possible dependencies.
        if deps is not None:
            deps = [Dependency(node, exts) for (node, exts) in deps.items()]
            task.add_deps(deps)

        return task

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
        #print("status_list", status_list)

        if only_min:
            return min(status_list)
        else:
            return status_list

    def check_status(self):
        """Check the status of the tasks."""
        # Recompute the status of the tasks
        for task in self:
            task.check_status()

        # Take into account possible dependencies.Use a list instead of generators 
        for task in self:
            if task.status <= task.S_SUB and all([status == task.S_OK for status in task.deps_status]): 
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

        # Build dirs and files.
        self.build(*args, **kwargs)

        # Initial setup
        self._setup(*args, **kwargs)

        # Submit tasks (does not block)
        self.submit_tasks(wait=wait)

    def read_etotal(self):
        """
        Reads the total energy from the GSR file produced by the task.

        Return a numpy array with the total energies in Hartree
        The array element is set to np.inf if an exception is raised while reading the GSR file.
        """
        if not self.all_done:
            raise self.Error("Some task is still in running/submitted state")

        etotal = []
        for task in self:
            # Open the GSR file and read etotal (Hartree)
            gsr_path = task.outdir.has_abiext("GSR")
            etot = np.inf
            if gsr_path:
                with ETSF_Reader(gsr_path) as r:
                    etot = r.read_value("etotal")
                
            etotal.append(etot)

        return etotal

    def parse_timers(self):
        """
        Parse the TIMER section reported in the ABINIT output files.

        Returns:
            `AbinitTimerParser` object
        """
        filenames = filter(os.path.exists, [task.output_file.path for task in self])
                                                                           
        parser = AbinitTimerParser()
        parser.parse(filenames)
                                                                           
        return parser


class IterativeWorkflow(Workflow):
    """
    This object defines a `Workflow` that produces `Tasks` until a particular 
    condition is satisfied (mainly used for convergence studies or iterative algorithms.)
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, strategy_generator, max_niter=25, workdir=None, manager=None):
        """
        Args:
            strategy_generator:
                Generator object that produces `Strategy` objects.
            max_niter:
                Maximum number of iterations. A negative value or zero value
                is equivalent to having an infinite number of iterations.
            workdir:
                Working directory.
            manager:
                `TaskManager` class.
        """
        super(IterativeWorkflow, self).__init__(workdir, manager)

        self.strategy_generator = strategy_generator

        self._max_niter = max_niter
        self.niter = 0

    @property
    def max_niter(self):
        return self._max_niter

    #def set_max_niter(self, max_niter):
    #    self._max_niter = max_niter

    #def set_inputs(self, inputs):
    #    self.strategy_generator = list(inputs)

    def next_task(self):
        """
        Generate and register a new `Task`.

        Returns: 
            New `Task` object
        """
        try:
            next_strategy = next(self.strategy_generator)

        except StopIteration:
            raise

        self.register(next_strategy)
        assert len(self) == self.niter

        return self[-1]

    def submit_tasks(self, *args, **kwargs):
        """
        Run the tasks till self.exit_iteration says to exit 
        or the number of iterations exceeds self.max_niter

        Returns: 
            dictionary with the final results
        """
        self.niter = 1

        while True:
            if self.niter > self.max_niter > 0:
                logger.debug("niter %d > max_niter %d" % (self.niter, self.max_niter))
                break

            try:
                task = self.next_task()
            except StopIteration:
                break

            # Start the task and block till completion.
            task.start(*args, **kwargs)
            task.wait()

            data = self.exit_iteration(*args, **kwargs)

            if data["exit"]:
                break

            self.niter += 1

    @abc.abstractmethod
    def exit_iteration(self, *args, **kwargs):
        """
        Return a dictionary with the results produced at the given iteration.
        The dictionary must contains an entry "converged" that evaluates to
        True if the iteration should be stopped.
        """


def check_conv(values, tol, min_numpts=1, mode="abs", vinf=None):
    """
    Given a list of values and a tolerance tol, returns the leftmost index for which

        abs(value[i] - vinf) < tol if mode == "abs"

    or

        abs(value[i] - vinf) / vinf < tol if mode == "rel"

    returns -1 if convergence is not achieved. By default, vinf = values[-1]

    Args:
        tol:
            Tolerance
        min_numpts:
            Minimum number of points that must be converged.
        mode:
            "abs" for absolute convergence, "rel" for relative convergence.
        vinf:
            Used to specify an alternative value instead of values[-1].
    """
    vinf = values[-1] if vinf is None else vinf

    if mode == "abs":
        vdiff = [abs(v - vinf) for v in values]
    elif mode == "rel":
        vdiff = [abs(v - vinf) / vinf for v in values]
    else:
        raise ValueError("Wrong mode %s" % mode)

    numpts = len(vdiff)
    i = -2

    if (numpts > min_numpts) and vdiff[-2] < tol:
        for i in range(numpts-1, -1, -1):
            if vdiff[i] > tol:
                break
        if (numpts - i -1) < min_numpts: i = -2

    return i + 1


def compute_hints(ecut_list, etotal, atols_mev, pseudo, min_numpts=1, stream=sys.stdout):
    de_low, de_normal, de_high = [a / (1000 * Ha_to_eV) for a in atols_mev]

    num_ene = len(etotal)
    etotal_inf = etotal[-1]

    ihigh   = check_conv(etotal, de_high, min_numpts=min_numpts)
    inormal = check_conv(etotal, de_normal)
    ilow    = check_conv(etotal, de_low)

    accidx = {"H": ihigh, "N": inormal, "L": ilow}

    table = []; app = table.append

    app(["iter", "ecut", "etotal", "et-e_inf [meV]", "accuracy",])
    for idx, (ec, et) in enumerate(zip(ecut_list, etotal)):
        line = "%d %.1f %.7f %.3f" % (idx, ec, et, (et-etotal_inf) * Ha_to_eV * 1.e+3)
        row = line.split() + ["".join(c for c,v in accidx.items() if v == idx)]
        app(row)

    if stream is not None:
        stream.write("pseudo: %s\n" % pseudo.name)
        pprint_table(table, out=stream)

    ecut_high, ecut_normal, ecut_low = 3 * (None,)
    exit = (ihigh != -1)

    if exit:
        ecut_low    = ecut_list[ilow]
        ecut_normal = ecut_list[inormal]
        ecut_high   = ecut_list[ihigh]

    aug_ratios = [1,]
    aug_ratio_low, aug_ratio_normal, aug_ratio_high = 3 * (1,)

    data = {
        "exit"       : ihigh != -1,
        "etotal"     : list(etotal),
        "ecut_list"  : ecut_list,
        "aug_ratios" : aug_ratios,
        "low"        : {"ecut": ecut_low, "aug_ratio": aug_ratio_low},
        "normal"     : {"ecut": ecut_normal, "aug_ratio": aug_ratio_normal},
        "high"       : {"ecut": ecut_high, "aug_ratio": aug_ratio_high},
        "pseudo_name": pseudo.name,
        "pseudo_path": pseudo.path,
        "atols_mev"  : atols_mev,
        "dojo_level" : 0,
    }

    return data


def plot_etotal(ecut_list, etotals, aug_ratios, **kwargs):
    """
    Uses Matplotlib to plot the energy curve as function of ecut

    Args:
        ecut_list:
            List of cutoff energies
        etotals:
            Total energies in Hartree, see aug_ratios
        aug_ratios:
            List augmentation rations. [1,] for norm-conserving, [4, ...] for PAW
            The number of elements in aug_ration must equal the number of (sub)lists
            in etotals. Example:

                - NC: etotals = [3.4, 4,5 ...], aug_ratios = [1,]
                - PAW: etotals = [[3.4, ...], [3.6, ...]], aug_ratios = [4,6]

        =========     ==============================================================
        kwargs        description
        =========     ==============================================================
        show          True to show the figure
        savefig       'abc.png' or 'abc.eps'* to save the figure to a file.
        =========     ==============================================================

    Returns:
        `matplotlib` figure.
    """
    show = kwargs.pop("show", True)
    savefig = kwargs.pop("savefig", None)

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    npts = len(ecut_list)

    if len(aug_ratios) != 1 and len(aug_ratios) != len(etotals):
        raise ValueError("The number of sublists in etotal must equal the number of aug_ratios")

    if len(aug_ratios) == 1:
        etotals = [etotals,]

    lines, legends = [], []

    emax = -np.inf
    for (aratio, etot) in zip(aug_ratios, etotals):
        emev = np.array(etot) * Ha_to_eV * 1000
        emev_inf = npts * [emev[-1]]
        yy = emev - emev_inf

        emax = np.max(emax, np.max(yy))

        line, = ax.plot(ecut_list, yy, "-->", linewidth=3.0, markersize=10)

        lines.append(line)
        legends.append("aug_ratio = %s" % aratio)

    ax.legend(lines, legends, 'upper right', shadow=True)

    # Set xticks and labels.
    ax.grid(True)
    ax.set_xlabel("Ecut [Ha]")
    ax.set_ylabel("$\Delta$ Etotal [meV]")
    ax.set_xticks(ecut_list)

    #ax.yaxis.set_view_interval(-10, emax + 0.01 * abs(emax))
    #ax.xaxis.set_view_interval(-10, 20)
    ax.yaxis.set_view_interval(-10, 20)

    ax.set_title("$\Delta$ Etotal Vs Ecut")

    if show:
        plt.show()

    if savefig is not None:
        fig.savefig(savefig)

    return fig


class PseudoConvergence(Workflow):

    def __init__(self, workdir, manager, pseudo, ecut_list, atols_mev,
                 toldfe=1.e-8, spin_mode="polarized", 
                 acell=(8, 9, 10), smearing="fermi_dirac:0.1 eV"):

        super(PseudoConvergence, self).__init__(workdir, manager)

        # Temporary object used to build the strategy.
        generator = PseudoIterativeConvergence(workdir, manager, pseudo, ecut_list, atols_mev,
                                               toldfe    = toldfe,
                                               spin_mode = spin_mode,
                                               acell     = acell,
                                               smearing  = smearing,
                                               max_niter = len(ecut_list),
                                              )
        self.atols_mev = atols_mev
        self.pseudo = Pseudo.aspseudo(pseudo)

        self.ecut_list = []
        for ecut in ecut_list:
            strategy = generator.strategy_with_ecut(ecut)
            self.ecut_list.append(ecut)
            self.register(strategy)

    def get_results(self):

        # Get the results of the tasks.
        wf_results = super(PseudoConvergence, self).get_results()

        etotal = self.read_etotal()
        data = compute_hints(self.ecut_list, etotal, self.atols_mev, self.pseudo)

        plot_etotal(data["ecut_list"], data["etotal"], data["aug_ratios"],
            show=False, savefig=self.path_in_workdir("etotal.pdf"))

        wf_results.update(data)

        if not monotonic(etotal, mode="<", atol=1.0e-5):
            logger.warning("E(ecut) is not decreasing")
            wf_results.push_exceptions("E(ecut) is not decreasing:\n" + str(etotal))

        #if kwargs.get("json_dump", True):
        #    wf_results.json_dump(self.path_in_workdir("results.json"))

        return wf_results


class PseudoIterativeConvergence(IterativeWorkflow):

    def __init__(self, workdir, manager, pseudo, ecut_list_or_slice, atols_mev,
                 toldfe=1.e-8, spin_mode="polarized", 
                 acell=(8, 9, 10), smearing="fermi_dirac:0.1 eV", max_niter=50,):
        """
        Args:
            workdir:
                Working directory.
            pseudo:
                string or Pseudo instance
            ecut_list_or_slice:
                List of cutoff energies or slice object (mainly used for infinite iterations).
            atols_mev:
                List of absolute tolerances in meV (3 entries corresponding to accuracy ["low", "normal", "high"]
            manager:
                `TaskManager` object.
            spin_mode:
                Defined how the electronic spin will be treated.
            acell:
                Lengths of the periodic box in Bohr.
            smearing:
                Smearing instance or string in the form "mode:tsmear". Default: FemiDirac with T=0.1 eV
        """
        self.pseudo = Pseudo.aspseudo(pseudo)

        self.atols_mev = atols_mev
        self.toldfe = toldfe
        self.spin_mode = spin_mode
        self.smearing = Smearing.assmearing(smearing)
        self.acell = acell

        if isinstance(ecut_list_or_slice, slice):
            self.ecut_iterator = iterator_from_slice(ecut_list_or_slice)
        else:
            self.ecut_iterator = iter(ecut_list_or_slice)

        # Construct a generator that returns strategy objects.
        def strategy_generator():
            for ecut in self.ecut_iterator:
                yield self.strategy_with_ecut(ecut)

        super(PseudoIterativeConvergence, self).__init__(strategy_generator(), 
              max_niter=max_niter, workdir=workdir, manager=manager, )

        if not self.isnc:
            raise NotImplementedError("PAW convergence tests are not supported yet")

    def strategy_with_ecut(self, ecut):
        """Return a Strategy instance with given cutoff energy ecut."""

        # Define the system: one atom in a box of lenghts acell.
        boxed_atom = AbiStructure.boxed_atom(self.pseudo, acell=self.acell)

        # Gamma-only sampling.
        gamma_only = KSampling.gamma_only()

        # Setup electrons.
        electrons = Electrons(spin_mode=self.spin_mode, smearing=self.smearing)

        # Don't write WFK files.
        extra_abivars = {
            "ecut" : ecut,
            "prtwf": 0,
            "toldfe": self.toldfe,
        }

        strategy = ScfStrategy(boxed_atom, self.pseudo, gamma_only,
                               spin_mode=self.spin_mode, smearing=self.smearing,
                               charge=0.0, scf_algorithm=None,
                               use_symmetries=True, **extra_abivars)

        return strategy

    @property
    def ecut_list(self):
        """The list of cutoff energies computed so far"""
        return [float(task.strategy.ecut) for task in self]

    def check_etotal_convergence(self, *args, **kwargs):
        return compute_hints(self.ecut_list, self.read_etotal(), self.atols_mev,
                             self.pseudo)

    def exit_iteration(self, *args, **kwargs):
        return self.check_etotal_convergence(self, *args, **kwargs)

    def get_results(self):
        """Return the results of the tasks."""
        wf_results = super(PseudoIterativeConvergence, self).get_results()

        data = self.check_etotal_convergence()

        ecut_list, etotal, aug_ratios = data["ecut_list"],  data["etotal"], data["aug_ratios"]

        plot_etotal(ecut_list, etotal, aug_ratios,
            show=False, savefig=self.path_in_workdir("etotal.pdf"))

        wf_results.update(data)

        if not monotonic(data["etotal"], mode="<", atol=1.0e-5):
            logger.warning("E(ecut) is not decreasing")
            wf_results.push_exceptions("E(ecut) is not decreasing\n" + str(etotal))

        #if kwargs.get("json_dump", True):
        #    wf_results.json_dump(self.path_in_workdir("results.json"))

        return wf_results


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
        """
        logger.debug("in on_ok with sender %s" % sender)

        if sender == self.ion_task and not self.transfer_done:
            # Get the relaxed structure.
            ion_structure = self.ion_task.read_final_structure()
            print("ion_structure", ion_structure)

            # Transfer it to the ioncell task (do it only once).
            self.ioncell_task.change_structure(ion_structure)
            self.transfer_done = True

            # Finally unlock ioncell_task so that we can submit it.
            self.ioncell_task.set_status(self.S_READY)

        base_results = super(RelaxWorkflow, self).on_ok(sender)
        return base_results


class DeltaFactorWorkflow(Workflow):

    def __init__(self, structure_or_cif, pseudo, kppa,
                 spin_mode="polarized", toldfe=1.e-8, smearing="fermi_dirac:0.1 eV",
                 accuracy="normal", ecut=None, pawecutdg=None, ecutsm=0.05, chksymbreak=0, workdir=None, manager=None): 
                 # FIXME Hack in chksymbreak
        """
        Build a `Workflow` for the computation of the deltafactor.

        Args:   
            structure_or_cif:
                Structure objec or string with the path of the CIF file.
            pseudo:
                String with the name of the pseudopotential file or `Pseudo` object.` object.` object.` 
            kppa:
                Number of k-points per atom.
            spin_mode:
                Spin polarization mode.
            toldfe:
                Tolerance on the energy (Ha)
            smearing:
                Smearing technique.
            workdir:
                String specifing the working directory.
            manager:
                `TaskManager` responsible for the submission of the tasks.
        """
        super(DeltaFactorWorkflow, self).__init__(workdir=workdir, manager=manager)

        if isinstance(structure_or_cif, Structure):
            structure = structure_or_cif
        else:
            # Assume CIF file
            structure = read_structure(structure_or_cif)

        self.pseudo = Pseudo.aspseudo(pseudo)

        structure = AbiStructure.asabistructure(structure)

        smearing = Smearing.assmearing(smearing)

        self._input_structure = structure

        v0 = structure.volume

        # From 94% to 106% of the equilibrium volume.
        self.volumes = v0 * np.arange(94, 108, 2) / 100.

        for vol in self.volumes:
            new_lattice = structure.lattice.scale(vol)

            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)
            new_structure = AbiStructure.asabistructure(new_structure)

            extra_abivars = dict(
                pawecutdg=pawecutdg,
                ecutsm=ecutsm,
                toldfe=toldfe,
                prtwf=0,
                paral_kgb=0,
            )

            if ecut is not None:
                extra_abivars.update({"ecut": ecut})

            ksampling = KSampling.automatic_density(new_structure, kppa,
                                                    chksymbreak=chksymbreak)

            scf_input = ScfStrategy(new_structure, self.pseudo, ksampling,
                                    accuracy=accuracy, spin_mode=spin_mode,
                                    smearing=smearing, **extra_abivars)

            self.register(scf_input, task_class=ScfTask)

    def get_results(self):
        num_sites = self._input_structure.num_sites

        etotal = ArrayWithUnit(self.read_etotal(), "Ha").to("eV")

        wf_results = super(DeltaFactorWorkflow, self).get_results()

        wf_results.update({
            "etotal"    : list(etotal),
            "volumes"   : list(self.volumes),
            "natom"     : num_sites,
            "dojo_level": 1,
        })

        try:
            #eos_fit = EOS.Murnaghan().fit(self.volumes/num_sites, etotal/num_sites)
            #print("murn",eos_fit)
            #eos_fit.plot(show=False, savefig=self.path_in_workdir("murn_eos.pdf"))

            # Use same fit as the one employed for the deltafactor.
            eos_fit = EOS.DeltaFactor().fit(self.volumes/num_sites, etotal/num_sites)

            eos_fit.plot(show=False, savefig=self.outdir.path_in("eos.pdf"))

            # FIXME: This object should be moved to pseudo_dojo.
            # Get reference results (Wien2K).
            from pseudo_dojo.refdata.deltafactor import df_database, df_compute
            wien2k = df_database().get_entry(self.pseudo.symbol)
                                                                                                 
            # Compute deltafactor estimator.
            dfact = df_compute(wien2k.v0, wien2k.b0_GPa, wien2k.b1, eos_fit.v0, eos_fit.b0_GPa, eos_fit.b1, b0_GPa=True)

            print("delta",eos_fit)
            print("Deltafactor = %.3f meV" % dfact)

            wf_results.update({
                "v0": eos_fit.v0,
                "b0": eos_fit.b0,
                "b0_GPa": eos_fit.b0_GPa,
                "b1": eos_fit.b1,
            })

        except EOS.Error as exc:
            wf_results.push_exceptions(exc)

        #if kwargs.get("json_dump", True):
        #    wf_results.json_dump(self.path_in_workdir("results.json"))

        # Write data for the computation of the delta factor
        with open(self.outdir.path_in("deltadata.txt"), "w") as fh:
            fh.write("# Deltafactor = %s meV\n" % dfact)
            fh.write("# Volume/natom [Ang^3] Etotal/natom [eV]\n")
            for (v, e) in zip(self.volumes, etotal):
                fh.write("%s %s\n" % (v/num_sites, e/num_sites))

        return wf_results

    def on_all_ok(self):
        return self.get_results()

    #def make_report(self, results, **kwargs):
    #    d = dict(v0=v0,
    #             b0_GPa=b0_GPa,
    #             b1=b1,
    #             dfact=dfact
    #            )
    #    if results.exceptions:
    #        d["_exceptions"] = str(results.exceptions)
    #                                                                                         
    #    d = {self.accuracy: d}


class G0W0_Workflow(Workflow):

    def __init__(self, scf_input, nscf_input, scr_input, sigma_inputs,
                 workdir=None, manager=None):
        """
        Workflow for G0W0 calculations.

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
        self.scf_task = scf_task = self.register(scf_input, task_class=ScfTask)

        # Construct the input for the NSCF run.
        self.nscf_task = nscf_task = self.register(nscf_input, deps={scf_task: "DEN"}, task_class=NscfTask)

        # Register the SCREENING run.
        self.scr_task = scr_task = self.register(scr_input, deps={nscf_task: "WFK"})

        # Register the SIGMA runs.
        if not isinstance(sigma_inputs, (list, tuple)): 
            sigma_inputs = [sigma_inputs]

        self.sigma_tasks = []
        for sigma_input in sigma_inputs:
            self.sigma_tasks.append(self.register(sigma_input, deps={nscf_task: "WFK", scr_task: "SCR"}))


#class SCGW_Workflow(Workflow):
#
#    def __init__(self, scr_input, sigma_input, workdir=None, manager=None):
#        """
#        Workflow for G0W0 calculations.
#
#        Args:
#            scr_input:
#                Input for the screening run or `ScrStrategy` object 
#            sigma_input:
#                Strategy for the self-energy run.
#            workdir:
#                Working directory of the calculation.
#            manager:
#                `TaskManager` object.
#        """
#        super(SCGW_Workflow, self).__init__(workdir=workdir, manager=manager)
#
#        # Register the SCREENING run.
#        self.scr_task = self.register(scr_input, deps={nscf_task: "WFK"})
#
#        # Register the SIGMA run.
#        self.sigma_task = self.register(sigma_input, deps={self.nscf_task: "WFK", self.scr_task: "SCR"})
#
#    def not_converged(self):
#       return self.sigma_task.not_converged()
#
#    def restart(self):
#        ext = "QPS"
#        qps_file = self.sigma_task.outdir.has_abiext(ext)
#        irdvars = irdvars_for_ext(ext)
#
#        if not qps_file:
#            raise TaskRestartError("Cannot find the QPS file to restart from.")
#
#        # Move the QPS file produced by the SIGMA task to 
#        # the indir of the SCR task and the indir of the SIGMA task.
#        scr_infile = self.scr_task.indir.path_in(os.path.basename(qps_file)
#        sigma_infile = self.sigma_task.indir.path_in(os.path.basename(qps_file)
#        shutil.copy(qps_file, scr_infile)
#        shutil.move(qps_file, sigma_infile)
#
#        # Add the appropriate variable for reading the QPS file.
#        self.scr_task.strategy.add_extra_abivars(irdvars)
#        self.sigma_task.strategy.add_extra_abivars(irdvars)
#
#        # Now we can resubmit the job.
#        #for task in self.
#        #    task.reset()
#        self._restart()


class BSEMDF_Workflow(Workflow):

    def __init__(self, scf_input, nscf_input, bse_input, workdir=None, manager=None):
        """
        Workflow for simple BSE calculations in which the self-energy corrections 
        are approximated by the scissors operator and the screening in modeled 
        with the model dielectric function.

        Args:
            scf_input:
                Input for the SCF run or `ScfStrategy` object.
            nscf_input:
                Input for the NSCF run or `NscfStrategy` object.
            bse_input:
                Input for the BSE run or `BSEStrategy` object.
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

        # Construct the input for the BSE run.
        self.bse_task = self.register(bse_input, deps={self.nscf_task: "WFK"}, task_class=HaydockBseTask)


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
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Workflow`.
        """
        ddb_files = filter(None, [task.outdir.has_abiext("DDB") for task in self])

        logger.debug("will call mrgddb to merge %s:\n" % str(ddb_files))
        assert len(ddb_files) == len(self)

        #if len(ddb_files) == 1:
        # Avoid the merge. Just move the DDB file to the outdir of the workflow

        # Final DDB file will be produced in the outdir of the workflow.
        out_ddb = self.outdir.path_in("out_DDB")
        desc = "DDB file merged by %s on %s" % (self.__class__.__name__, time.asctime())

        mrgddb = wrappers.Mrgddb(verbose=1)
        mrgddb.merge(ddb_files, out_ddb=out_ddb, description=desc, cwd=self.outdir.path)

    def merge_gkk_files(self):
        """
        This method is called when all the q-points have been computed.
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Workflow`.
        """
        gkk_files = filter(None, [task.outdir.has_abiext("GKK") for task in self])
                                                                                         
        logger.debug("Will call mrggkk to merge %s:\n" % str(gkk_files))
        assert len(gkk) == len(self)

        #if len(gkk) == 1:
        # Avoid the merge. Just move the GKK file to the outdir of the workflow
                                                                                         
        # Final GKK file will be produced in the outdir of the workflow.
        out_ggk = self.outdir.path_in("out_GKK")

        mrggkk = wrappers.Mrggkk(verbose=1)
        raise NotImplementedError("Have to check mrggkk")
        #mrggkk.merge(gswfk_file, dfpt_files, gkk_files, out_fname, binascii=0, cwd=self.outdir.path)

    def on_all_ok(self):
        """
        This method is called when all the q-points have been computed.
        Ir runs `mrgddb` in sequential on the local machine to produce
        the final DDB file in the outdir of the `Workflow`.
        """
        # Merge DDB files.
        self.merge_ddb_files()

        # Merge GKK files.
        #self.merge_gkk_files()

        results = dict(
            returncode=0,
            message="DDB merge done",
        )

        return results


class WorkflowResults(dict, MSONable):
    """
    Dictionary used to store some of the results produce by a Task object
    """
    _MANDATORY_KEYS = [
        "task_results",
    ]

    _EXC_KEY = "_exceptions"

    def __init__(self, *args, **kwargs):
        super(WorkflowResults, self).__init__(*args, **kwargs)

        if self._EXC_KEY not in self:
            self[self._EXC_KEY] = []

    @property
    def exceptions(self):
        return self[self._EXC_KEY]

    def push_exceptions(self, *exceptions):
        for exc in exceptions:
            newstr = str(exc)
            if newstr not in self.exceptions:
                self[self._EXC_KEY] += [newstr]

    def assert_valid(self):
        """
        Returns empty string if results seem valid.

        The try assert except trick allows one to get a string with info on the exception.
        We use the += operator so that sub-classes can add their own message.
        """
        # Validate tasks.
        for tres in self.task_results:
            self[self._EXC_KEY] += tres.assert_valid()

        return self[self._EXC_KEY]

    @property
    def to_dict(self):
        d = {k: v for k,v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        mydict = {k: v for k, v in d.items() if k not in ["@module", "@class"]}
        return cls(mydict)
