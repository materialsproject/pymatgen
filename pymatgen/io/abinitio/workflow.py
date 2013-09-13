"""
Abinit workflows
"""
from __future__ import division, print_function

import sys
import os
import shutil
import abc
import collections
import functools
import numpy as np
import cPickle as pickle

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.design_patterns import Enum, AttrDict
import pymatgen.core.physical_constants as const 
from pymatgen.serializers.json_coders import MSONable, json_pretty_dump
from pymatgen.io.smartio import read_structure
from pymatgen.util.num_utils import iterator_from_slice, chunks
from pymatgen.io.abinitio.task import task_factory, Task, AbinitTask
from pymatgen.io.abinitio.strategies import Strategy

from .utils import File
from .netcdf import ETSF_Reader
from .abiobjects import Smearing, AbiStructure, KSampling, Electrons
from .pseudos import Pseudo
from .strategies import ScfStrategy
from .task import RunMode
from .eos import EOS

#import logging
#logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

__all__ = [
    "Workflow",
]

################################################################################


def map_method(method):
    """Decorator that calls item.method for all items in a iterable object."""
    @functools.wraps(method)
    def wrapped(iter_obj, *args, **kwargs):
        return [getattr(item, method.__name__)(*args, **kwargs)
                for item in iter_obj]
    return wrapped



class Product(object):
    """
    A product represents a file produced by an AbinitTask instance, file
    that is needed by another task in order to start the calculation.
    """
    # TODO
    # It would be nice to pass absolute paths to abinit with getden_path
    # so that I can avoid creating symbolic links before running but
    # the presence of the C-bindings complicates the implementation
    # (gfortran SIGFAULTs if I add strings to dataset_type!
    _EXT2ABIVARS = {
        "_DEN": {"irdden": 1},
        "_WFK": {"irdwfk": 1},
        "_SCR": {"irdscr": 1},
        "_QPS": {"irdqps": 1},
    }

    def __init__(self, ext, path):
        self.ext = ext
        self.file = File(path)

    def __str__(self):
        return "ext = %s, file = %s" % (self.ext, self.file)

    def get_filepath(self):
        return self.file.path

    def get_abivars(self):
        return self._EXT2ABIVARS[self.ext].copy()


class WorkLink(object):
    """
    This object describes the dependencies among the tasks contained in a Work instance.

    A WorkLink is a task that produces a list of products (files) that are
    reused by the other tasks belonging to a Work instance.
    One usually instantiates the object by calling work.register and produces_exts.
    Example:

        # Register the SCF task in work and get the link.
        scf_link = work.register(scf_strategy)

        # Register the NSCF calculation and its dependency on the SCF run.
        nscf_link = work.register(nscf_strategy, links=scf_link.produces_exts("_DEN"))
    """
    def __init__(self, task, exts=None):
        """
        Args:
            task:
                The task associated to the link.
            exts:
                Extensions of the output files that are needed for running the other tasks.
        """
        self._task = task

        self._products = []

        if exts is not None:
            if isinstance(exts, str):
                exts = [exts,]

            for ext in exts:
                prod = Product(ext, task.odata_path_from_ext(ext))
                self._products.append(prod)

    def __str__(self):
        s = "%s: task %s with products\n %s" % (
                self.__class__.__name__, repr(self._task), "\n".join(str(p) for p in self.products))
        return s

    @property
    def products(self):
        return self._products

    def produces_exts(self, exts):
        return WorkLink(self._task, exts=exts)

    def get_abivars(self):
        """
        Returns a dictionary with the abinit variables that must
        be added to the input file in order to connect the two tasks.
        """
        abivars = {}
        for prod in self._products:
            abivars.update(prod.get_abivars())
        return abivars

    def get_filepaths_and_exts(self):
        """Returns the paths of the output files produced by self and its extensions"""
        filepaths = [prod.get_filepath() for prod in self._products]
        exts = [prod.ext for prod in self._products]
        return filepaths, exts

    @property
    def status(self):
        """The status of the link, equivalent to the task status"""
        return self._task.status


class WorkflowError(Exception):
    """Base class for the exceptions raised by Workflow objects."""


class BaseWorkflow(object):
    __metaclass__ = abc.ABCMeta

    Error = WorkflowError

    # Basename of the pickle database.
    _PICKLE_FNAME = "__workflow__.pickle"

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
        Wait for child processed to terminate. Set and return returncode
        attributes.
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
    def ncpus_reserved(self):
        """Returns the number of CPUs reserved in this moment."""
        ncpus = 0
        for task in self:
            if task.is_allocated:
                ncpus += task.tot_ncpus

        return ncpus

    def fetch_task_to_run(self):
        """
        Returns the first task that is ready to run or None if no task can be submitted at present"

        Raises `StopIteration` if all tasks are done.
        """
        for task in self:
            if task.can_run:
                return task

        # All the tasks are done so raise an exception 
        # that will be handled by the client code.
        if all([task.is_completed for task in self]):
            raise StopIteration("All tasks completed.")

        # No task found, this usually happens when we have dependencies. 
        # Beware of possible deadlocks here!
        return None

    @abc.abstractmethod
    def setup(self, *args, **kwargs):
        """Method called before submitting the calculations."""

    def _setup(self, *args, **kwargs):
        self.setup(*args, **kwargs)

    def get_results(self, *args, **kwargs):
        """
        Method called once the calculations completes.

        The base version returns a dictionary task_name : TaskResults for each task in self.
        """
        return WorkFlowResults(task_results={task.name: task.results for task in self})

    def pickle_dump(self, protocol=0):
        """Save the status of the object in pickle format."""
        filepath = os.path.join(self.workdir, self._PICKLE_FNAME)

        with open(filepath, "w") as fh:
            pickle.dump(self, fh, protocol=protocol)

    @classmethod
    def pickle_load(cls, path):
        """
        Load the object from a pickle file. 

        Args:
            path:
                filepath or directory name.
        """
        if os.path.isdir(path):
            path = os.path.join(path, cls._PICKLE_FNAME)

        # TODO atomic transaction.
        with open(path, "rb") as fh:
            new = pickle.load(fh)

        #new.__class__ = cls
        return new


class WorkFlowResults(dict, MSONable):
    """
    Dictionary used to store some of the results produce by a Task object
    """
    _MANDATORY_KEYS = [
        "task_results",
    ]
    _EXC_KEY = "_exceptions"

    def __init__(self, *args, **kwargs):
        super(WorkFlowResults, self).__init__(*args, **kwargs)

        if self._EXC_KEY not in self:
            self[self._EXC_KEY] = []

    @property
    def exceptions(self):
        return self[self._EXC_KEY]

    def push_exceptions(self, *exceptions):
        for exc in exceptions:
            newstr = str(exc)
            if newstr not in self.exceptions:
                self[self._EXC_KEY] += [newstr,]

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
        mydict = {k: v for k,v in d.items() if k not in ["@module", "@class",]}
        return cls(mydict)
    

class Workflow(BaseWorkflow):
    """
    A Workflow is a list of (possibly connected) tasks.
    """
    Error = WorkflowError

    def __init__(self, workdir, runmode):
        """
        Args:
            workdir:
                Path to the working directory.
            runmode:
                RunMode instance or string "sequential"
        """
        self.workdir = os.path.abspath(workdir)

        self.runmode = RunMode.asrunmode(runmode)

        self._tasks = []

        # Dict with the dependencies of each task, indexed by task.id
        self._links_dict = collections.defaultdict(list)

    def __len__(self):
        return len(self._tasks)

    def __iter__(self):
        return self._tasks.__iter__()

    def chunks(self, chunk_size):
        """Yield successive chunks of tasks of lenght chunk_size."""
        for tasks in chunks(self, chunk_size):
            yield tasks

    def __getitem__(self, slice):
        return self._tasks[slice]

    def __repr__(self):
        return "<%s at %s, workdir = %s>" % (self.__class__.__name__, id(self), str(self.workdir))

    def __str__(self):
        return self.__repr__()

    @property
    def processes(self):
        return [task.process for task in self]

    @property
    def alldone(self):
        """True if all the Task in the `Workflow` are completed."""
        return all([task.status == Task.S_DONE for task in self])

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return all(task.isnc for task in self)

    @property
    def ispaw(self):
        """True if PAW calculation."""
        return all(task.ispaw for task in self)

    #@property
    #def to_dict(self):
    #    d = dict(
    #        workdir=self.workdir,
    #        runmode=self.runmode.to_dict,
    #        kwargs=self._kwargs
    #    )
    #    d["@module"] = self.__class__.__module__
    #    d["@class"] = self.__class__.__name__
    #    return d

    #@classmethod
    #def from_dict(cls, d):
    #    return cls(d["workdir"], d["runmode"])

    def register(self, obj, links=()):
        """
        Registers a new task and add it to the internal list, taking into account possible dependencies.

        Args:
            obj:
                `Strategy` object or `AbinitInput` instance.
                if Strategy object, we create a new `AbinitTask` from the input strategy and add it to the list.

        Returns:   
            `WorkLink` object
        """
        # Handle possible dependencies.
        if links and not isinstance(links, collections.Iterable):
            links = [links,]

        task_id = len(self) + 1
        task_workdir = os.path.join(self.workdir, "task_" + str(task_id))

        if isinstance(obj, Strategy):
            # Create the new task (note the factory so that we create subclasses easily).
            task = task_factory(obj, task_workdir, self.runmode, task_id=task_id, links=links)

        else:
            # Create the new task from the input. Note that no subclasses are instanciated here.
            try:
                task = AbinitTask.from_input(obj, task_workdir, self.runmode, task_id=task_id, links=links)

            except:
                raise TypeError("Don't know how to create a Task from object type %s" % str(type(obj)))

        self._tasks.append(task)

        if links:
            self._links_dict[task_id].extend(links)
            print("task_id %s needs\n %s" % (task_id, [str(l) for l in links]))

        return WorkLink(task)

    def path_in_workdir(self, filename):
        """Create the absolute path of filename in the workind directory."""
        return os.path.join(self.workdir, filename)

    def setup(self, *args, **kwargs):
        """
        Method called before running the calculations.
        The default implementation is empty.
        """

    def build(self, *args, **kwargs):
        """Creates the top level directory."""
        # Create top level directory.
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # Build dirs and files of each task.
        for task in self:
            task.build(*args, **kwargs)

        #for task in self:
        #    task.connect()

    def get_status(self, only_highest_rank=False):
        """Get the status of the tasks in self."""
        status_list = [task.status for task in self]

        if only_highest_rank:
            return max(status_list)
        else:
            return status_list

    def recheck_status(self):
        for task in self:
            task.recheck_status()

    def show_inputs(self, stream=sys.stdout):
        lines = []
        app = lines.append

        width = 120
        for task in self:
            app("\n")
            app(repr(task)+"\n")
            app("Input file: %s\n" % task.input_file.path)
            app(str(task.make_input()))
            app(width*"=" + "\n")

        stream.write("\n".join(lines))

    def rmtree(self, exclude_wildcard=""):
        """
        Remove all files and directories in the working directory

        Args:
            exclude_wildcard:
                Optional string with regular expressions separated by |.
                Files matching one of the regular expressions will be preserved.
                example: exclude_wildard="*.nc|*.txt" preserves all the files
                whose extension is in ["nc", "txt"].

        """
        if not exclude_wildcard:
            shutil.rmtree(self.workdir)

        else:
            pats = exclude_wildcard.split("|")

            import fnmatch

            def keep(fname):
                for pat in pats:
                    if fnmatch.fnmatch(fname, pat):
                        return True
                return False

            for dirpath, dirnames, filenames in os.walk(self.workdir):
                for fname in filenames:
                    path = os.path.join(dirpath, fname)
                    if not keep(fname):
                        os.remove(path)

    def rm_tmpdatadir(self):
        """Remove all the tmpdata directories."""
        for task in self:
            task.rm_tmpdatadir()

    def rm_indatadir(self):
        """Remove all the indata directories."""
        for task in self:
            task.rm_indatadir()

    def move(self, dst, isabspath=False):
        """
        Recursively move self.workdir to another location. This is similar to the Unix "mv" command.
        The destination path must not already exist. If the destination already exists
        but is not a directory, it may be overwritten depending on os.rename() semantics.

        Be default, dst is located in the parent directory of self.workdir, use isabspath=True
        to specify an absolute path.
        """
        if not isabspath:
            dst = os.path.join(os.path.dirname(self.workdir), dst)

        shutil.move(self.workdir, dst)

    def submit_tasks(self, *args, **kwargs):
        """
        Submits the task in self.
        """
        for task in self:
            task.start(*args, **kwargs)
            # FIXME
            task.wait()

    def start(self, *args, **kwargs):
        """
        Start the work. Calls build and _setup first, then the tasks are submitted.
        Non-blocking call
        """
        # Build dirs and files.
        self.build(*args, **kwargs)

        # Initial setup
        self._setup(*args, **kwargs)

        # Submit tasks (does not block)
        self.submit_tasks(*args, **kwargs)

    def read_etotal(self):
        """
        Reads the total energy from the GSR file produced by the task.

        Return a numpy array with the total energies in Hartree
        The array element is set to np.inf if an exception is raised while reading the GSR file.
        """
        if not self.alldone:
            raise self.Error("Some task is still in running/submitted state")

        etotal = []
        for task in self:
            # Open the GSR file and read etotal (Hartree)
            with ETSF_Reader(task.odata_path_from_ext("_GSR")) as ncdata:
                etotal.append(ncdata.read_value("etotal"))

        return etotal

    def json_dump(self, filename):
        json_pretty_dump(self.to_dict, filename)
                                                  
    @classmethod
    def json_load(cls, filename):
        return cls.from_dict(json_load(filename))


class IterativeWork(Workflow):
    """
    This object defined a workflow that produces tasks until a particular 
    condition is satisfied (mainly used for simple convergence studies).
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, workdir, runmode, strategy_generator, max_niter=25):
        """
        Args:
            workdir:
                Working directory.
            runmode:
                `RunMode` object.
            strategy_generator:
                Strategy generator.
            max_niter:
                Maximum number of iterations. A negative value or zero value
                is equivalent to having an infinite number of iterations.
        """
        super(IterativeWork, self).__init__(workdir, runmode)

        self.strategy_generator = strategy_generator

        self.max_niter = max_niter

    def next_task(self):
        """
        Generate and register a new task

        Returns: 
            New `Task` object
        """
        try:
            next_strategy = next(self.strategy_generator)
        except StopIteration:
            raise StopIteration

        self.register_task(next_strategy)
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
                print("niter %d > max_niter %d" % (self.niter, self.max_niter))
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

##########################################################################################

def strictly_increasing(values):
    return all(x < y for x, y in zip(values, values[1:]))


def strictly_decreasing(values):
    return all(x > y for x, y in zip(values, values[1:]))


def non_increasing(values):
    return all(x >= y for x, y in zip(values, values[1:]))


def non_decreasing(values):
    return all(x <= y for x, y in zip(values, values[1:]))


def monotonic(values, mode="<", atol=1.e-8):
    """
    Returns False if values are not monotonic (decreasing|increasing).
    mode is "<" for a decreasing sequence, ">" for an increasing sequence.
    Two numbers are considered equal if they differ less that atol.

    .. warning:
        Not very efficient for large data sets.

    >>> values = [1.2, 1.3, 1.4]
    >>> monotonic(values, mode="<")
    False
    >>> monotonic(values, mode=">")
    True
    """
    if len(values) == 1:
        return True

    if mode == ">":
        for i in range(len(values)-1):
            v, vp = values[i], values[i+1]
            if abs(vp - v) > atol and vp <= v:
                return False

    elif mode == "<":
        for i in range(len(values)-1):
            v, vp = values[i], values[i+1]
            if abs(vp - v) > atol and vp >= v:
                return False

    else:
        raise ValueError("Wrong mode %s" % str(mode))

    return True


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
    de_low, de_normal, de_high = [a / (1000 * const.HA_TO_EV) for a in atols_mev]

    num_ene = len(etotal)
    etotal_inf = etotal[-1]

    ihigh   = check_conv(etotal, de_high, min_numpts=min_numpts)
    inormal = check_conv(etotal, de_normal)
    ilow    = check_conv(etotal, de_low)

    accidx = {"H": ihigh, "N": inormal, "L": ilow}

    table = []; app = table.append

    app(["iter", "ecut", "etotal", "et-e_inf [meV]", "accuracy",])
    for idx, (ec, et) in enumerate(zip(ecut_list, etotal)):
        line = "%d %.1f %.7f %.3f" % (idx, ec, et, (et-etotal_inf) * const.HA_TO_EV * 1.e+3)
        row = line.split() + ["".join(c for c,v in accidx.items() if v == idx)]
        app(row)

    if stream is not None:
        from pymatgen.util.string_utils import pprint_table
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

##########################################################################################

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
        emev = const.Ha2meV(etot)
        emev_inf = npts * [emev[-1]]
        yy = emev - emev_inf

        emax = max(emax, np.max(yy))

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

##########################################################################################

class PseudoConvergence(Workflow):

    def __init__(self, workdir, pseudo, ecut_list, atols_mev,
                 runmode="sequential", toldfe=1.e-8, spin_mode="polarized", 
                 acell=(8, 9, 10), smearing="fermi_dirac:0.1 eV",):

        super(PseudoConvergence, self).__init__(workdir, runmode)

        # Temporary object used to build the strategy.
        generator = PseudoIterativeConvergence(workdir, pseudo, ecut_list, atols_mev,
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

    def get_results(self, *args, **kwargs):

        # Get the results of the tasks.
        wf_results = super(PseudoConvergence, self).get_results()

        etotal = self.read_etotal()
        data = compute_hints(self.ecut_list, etotal, self.atols_mev, self.pseudo)

        plot_etotal(data["ecut_list"], data["etotal"], data["aug_ratios"],
            show=False, savefig=self.path_in_workdir("etotal.pdf"))

        wf_results.update(data)

        if not monotonic(etotal, mode="<", atol=1.0e-5):
            print("E(ecut) is not decreasing")
            wf_results.push_exceptions("E(ecut) is not decreasing:\n" + str(etotal))

        if kwargs.get("json_dump", True):
            wf_results.json_dump(self.path_in_workdir("results.json"))

        return wf_results


class PseudoIterativeConvergence(IterativeWork):

    def __init__(self, workdir, pseudo, ecut_list_or_slice, atols_mev,
                 runmode="sequential", toldfe=1.e-8, spin_mode="polarized", 
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

        super(PseudoIterativeConvergence, self).__init__(
            workdir, runmode, strategy_generator(), max_niter=max_niter)

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

    def get_results(self, *args, **kwargs):
        # Get the results of the tasks.
        wf_results = super(PseudoIterativeConvergence, self).get_results()

        data = self.check_etotal_convergence()

        ecut_list, etotal, aug_ratios = data["ecut_list"],  data["etotal"], data["aug_ratios"]

        plot_etotal(ecut_list, etotal, aug_ratios,
            show=False, savefig=self.path_in_workdir("etotal.pdf"))

        wf_results.update(data)

        if not monotonic(data["etotal"], mode="<", atol=1.0e-5):
            print("E(ecut) is not decreasing")
            wf_results.push_exceptions("E(ecut) is not decreasing\n" + str(etotal))

        if kwargs.get("json_dump", True):
            wf_results.json_dump(self.path_in_workdir("results.json"))

        return wf_results

################################################################################

class BandStructure(Workflow):
    """Workflow for band structure calculations."""
    def __init__(self, workdir, runmode, scf_strategy, nscf_strategy,
                 dos_strategy=None):
        """
        Args:
            workdir:
                Working directory.
            runmode:
                `RunMode` instance.
            scf_strategy:
                `SCFStrategy` instance
            nscf_strategy:
                `NSCFStrategy` instance defining the band structure calculation.
            dos_strategy:
                `NSCFStrategy` instance defining the DOS calculation. 
                DOS is computed only if dos_strategy is not None.
        """
        super(BandStructure, self).__init__(workdir, runmode)

        # Register the GS-SCF run.
        scf_link = self.register(scf_strategy)

        # Register the NSCF run and its dependency
        self.register(nscf_strategy, links=scf_link.produces_exts("_DEN"))

        # Add DOS computation
        if dos_strategy is not None:
            self.register(dos_strategy, links=scf_link.produces_exts("_DEN"))

################################################################################

class Relaxation(Workflow):

    def __init__(self, workdir, runmode, relax_strategy):
        """
        Args:
            workdir:
                Working directory.
            runmode:
                `RunMode` instance.
            relax_strategy:
                `RelaxStrategy` instance
        """
        super(Relaxation, self).__init__(workdir, runmode)

        link = self.register(relax_strategy)

################################################################################

class DeltaTest(Workflow):

    def __init__(self, workdir, runmode, structure_or_cif, pseudos, kppa,
                 spin_mode="polarized", toldfe=1.e-8, smearing="fermi_dirac:0.1 eV",
                 accuracy="normal", ecut=None, ecutsm=0.05, chksymbreak=0): # FIXME Hack

        super(DeltaTest, self).__init__(workdir, runmode)

        if isinstance(structure_or_cif, Structure):
            structure = structure_or_cif
        else:
            # Assume CIF file
            structure = read_structure(structure_or_cif)

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

            extra_abivars = {
                "ecutsm": ecutsm,
                "toldfe": toldfe,
                "prtwf" : 0,
                "paral_kgb": 0,
            }

            if ecut is not None:
                extra_abivars.update({"ecut": ecut})

            ksampling = KSampling.automatic_density(new_structure, kppa,
                                                    chksymbreak=chksymbreak)

            scf_strategy = ScfStrategy(new_structure, pseudos, ksampling,
                                       accuracy=accuracy, spin_mode=spin_mode,
                                       smearing=smearing, **extra_abivars)

            self.register(scf_strategy)

    def get_results(self, *args, **kwargs):

        num_sites = self._input_structure.num_sites

        etotal = const.Ha2eV(self.read_etotal())

        wf_results = super(DeltaTest, self).get_results()

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
            print("delta",eos_fit)
            eos_fit.plot(show=False, savefig=self.path_in_workdir("eos.pdf"))

            wf_results.update({
                "v0": eos_fit.v0,
                "b0": eos_fit.b0,
                "b0_GPa": eos_fit.b0_GPa,
                "b1": eos_fit.b1,
            })

        except EOS.Error as exc:
            wf_results.push_exceptions(exc)

        if kwargs.get("json_dump", True):
            wf_results.json_dump(self.path_in_workdir("results.json"))

        # Write data for the computation of the delta factor
        with open(self.path_in_workdir("deltadata.txt"), "w") as fh:
            fh.write("# Volume/natom [Ang^3] Etotal/natom [eV]\n")
            for (v, e) in zip(self.volumes, etotal):
                fh.write("%s %s\n" % (v/num_sites, e/num_sites))

        return wf_results

################################################################################

class GW_Workflow(Workflow):

    def __init__(self, workdir, runmode, scf_strategy, nscf_strategy,
                 scr_strategy, sigma_strategy):
        """
        Workflow for GW calculations.

        Args:
            workdir:
                Working directory of the calculation.
            runmode:
                `RunMode` instance.
            scf_strategy:
                `SCFStrategy` instance
            nscf_strategy:
                `NSCFStrategy` instance
            scr_strategy:
                Strategy for the screening run.
            sigma_strategy:
                Strategy for the self-energy run.
        """
        super(GW_Workflow, self).__init__(workdir, runmode)

        # Register the GS-SCF run.
        scf_link = self.register(scf_strategy)

        # Construct the input for the NSCF run.
        nscf_link = self.register(nscf_strategy, links=scf_link.produces_exts("_DEN"))

        # Register the SCREENING run.
        screen_link = self.register(scr_strategy, links=nscf_link.produces_exts("_WFK"))

        # Register the SIGMA run.
        self.register(sigma_strategy, links=[nscf_link.produces_exts("_WFK"),
                                             screen_link.produces_exts("_SCR")])


class BSEMDF_Workflow(Workflow):

    def __init__(self, workdir, runmode, scf_strategy, nscf_strategy, bse_strategy): 
        """
        Workflow for simple BSE calculations in which the self-energy corrections 
        are approximated by the scissors operator and the screening in modeled 
        with the model dielectric function.

        Args:
            workdir:
                Working directory of the calculation.
            runmode:
                Run mode.
            scf_strategy:
                ScfStrategy instance
            nscf_strategy:
                NscfStrategy instance
            bse_strategy:
                BSEStrategy instance.
        """
        super(BSEMDF_Workflow, self).__init__(workdir, runmode)

        # Register the GS-SCF run.
        scf_link = self.register(scf_strategy)

        # Construct the input for the NSCF run.
        nscf_link = self.register(nscf_strategy, links=scf_link.produces_exts("_DEN"))

        # Construct the input for the BSE run.
        bse_link = self.register(bse_strategy, links=nscf_link.produces_exts("_WFK"))

