"""
Classes defining Abinit calculations and workflows
"""
from __future__ import division, print_function

import sys
import os
import os.path
import shutil
import abc
import collections
import functools
import numpy as np

from pprint import pprint

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.design_patterns import Enum
from pymatgen.core.physical_constants import Bohr2Ang, Ang2Bohr, Ha2eV, Ha_eV, Ha2meV
from pymatgen.serializers.json_coders import MSONable, json_pretty_dump #, PMGJSONDecoder
from pymatgen.io.smartio import read_structure
from pymatgen.util.num_utils import iterator_from_slice, chunks
from pymatgen.io.abinitio.task import task_factory, Task

from .utils import abinit_output_iscomplete, File
from .netcdf import GSR_Reader
from .abiobjects import Smearing
from .pseudos import Pseudo, PseudoDatabase, PseudoTable, get_abinit_psp_dir
from .input import Input, ElectronsCard, SystemCard, ControlCard, KpointsCard 
from .task import RunMode

#import logging
#logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

#__all__ = [
#]

##########################################################################################

def map_method(method):
    "Decorator that calls item.method for all items in a iterable object."
    @functools.wraps(method)
    def wrapped(iter_obj, *args, **kwargs):
        return [getattr(item, method.__name__)(*args, **kwargs) for item in iter_obj]
    return wrapped

##########################################################################################

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
    _ext2abivars = {
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
        return self._ext2abivars[self.ext].copy()

class WorkLink(object):
    """
    This object describes the dependencies among the tasks contained in a Work instance.

    A  WorkLink is a task that produces a list of products (files) that are 
    reused by the tasks belonging to a Work instance. 
    One usually instantiates the object by calling work.register_input and produces_exts.
    Example:

        # Generate the input file for the SCF calculation
        scf_input = ...

        # Register the SCF input in work and get the link.
        scf_link = work.register_input(scf_input)

        # Generate the input file for NSCF calculation
        nscf_input = ...

        # Register the NSCF input and its dependency on the SCF run.
        nscf_link = work.register_input(nscf_input, links=scf_link.produces_exts("_DEN"))
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
        abivars =  {}
        for prod in self._products:
            abivars.update(prod.get_abivars())
        return abivars

    def get_filepaths_and_exts(self):
        "Returns the paths of the output files produced by self and its extensions"
        filepaths = [prod.get_filepath() for prod in self._products]
        exts = [prod.ext for prod in self._products]
        return filepaths, exts

    @property
    def status(self):
        "The status of the link, equivalent to the task status"
        return self._task.status

##########################################################################################

class WorkError(Exception):
    "Base class for the exceptions raised by Work objects"

class BaseWork(object):
    __metaclass__ = abc.ABCMeta

    Error = WorkError

    # interface modeled after subprocess.Popen
    @abc.abstractproperty
    def processes(self):
        "Return a list of objects that support the subprocess.Popen protocol."

    def poll(self):
        "Check if all child processes have terminated. Set and return returncode attribute."
        return [task.poll() for task in self]

    def wait(self):
        "Wait for child processed to terminate. Set and return returncode attributes."
        return [task.wait() for task in self]

    def communicate(self, input=None):
        """
        Interact with processes: Send data to stdin. Read data from stdout and stderr, until end-of-file is reached. 
        Wait for process to terminate. The optional input argument should be a string to be sent to the 
        child processed, or None, if no data should be sent to the children

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
        "Returns the number of CPUs reserved in this moment."
        ncpus = 0
        for task in self:
            if task.status in [task.S_SUB, task.S_RUN]:
                ncpus += task.tot_ncpus
        return ncpus

    def fetch_task_to_run(self):
        """
        Returns the first task that is ready to run or None if no task can be submitted at present" 

        Raises StopIteration if all tasks are done. 
        """
        for task in self:
            # The task is ready to run if its status is S_READY and all the other links (if any) are done!
            if (task.status == task.S_READY) and all([link_stat==task.S_DONE for link_stat in task.links_status]):
                return task

        # All the tasks are done so raise an exception that will be handled by the client code.
        if all([task.status == task.S_DONE for task in self]): 
            raise StopIteration

        # No task found, this usually happens when we have dependencies. Beware of possible deadlocks here!
        return None

    @abc.abstractmethod
    def setup(self, *args, **kwargs):
        """
        Method called before submitting the calculations. 
        """
                                                         
    def _setup(self, *args, **kwargs):
        self.setup(*args, **kwargs)

    def get_results(self, *args, **kwargs):
        """
        Method called once the calculations completes.
        The base version returns a dictionary task_name : TaskResults for each task in self.
        """
        work_results = collections.OrderedDict()
        for task in self:
            work_results[task.name] = task.results
        return work_results
                                                         
##########################################################################################

class Work(BaseWork, MSONable):
    """
    A work is a list of (possibly connected) tasks.
    """
    Error = WorkError

    @classmethod
    def from_task(cls, task):
        "Build a Work instance from a task object"
        workdir, tail = os.path.dirname(task.workdir)
        new = cls(workdir, taks.runmode)
        new.register_input(task.input)
        return new

    def __init__(self, workdir, runmode, **kwargs):
        """
        Args:
            workdir:
                Path to the working directory.
            runmode:
                RunMode instance or string "sequential"
        """
        self.workdir = os.path.abspath(workdir)

        self.runmode = RunMode.asrunmode(runmode)

        self._kwargs = kwargs

        self._tasks = []

        # Dict with the dependencies of each task, indexed by task.id
        self._links_dict = collections.defaultdict(list)

    def __len__(self):
        return len(self._tasks)

    def __iter__(self):
        return self._tasks.__iter__()

    def chunks(self, chunk_size):
        "Yield successive chunks of tasks of lenght chunk_size."
        for tasks in chunks(self, chunk_size):
            yield tasks

    def __getitem__(self, slice):
        return self._tasks[slice]

    def __repr__(self):
        return "<%s at %s, workdir = %s>" % (self.__class__.__name__, id(self), str(self.workdir))

    @property
    def to_dict(self):
        d = {"workdir": self.workdir,
             "runmode": self.runmode.to_dict,
             "kwargs" : self._kwargs,
            }
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        return Work(d["workdir"], d["runmode"], **d["kwargs"])

    @property
    def alldone(self):
        return all([task.status == Task.S_DONE for task in self])

    @property
    def isnc(self):
        "True if norm-conserving calculation"
        return all(task.isnc for task in self)
                                                
    @property
    def ispaw(self):
        "True if PAW calculation"
        return all(task.ispaw for task in self)

    def path_in_workdir(self, filename):
        "Create the absolute path of filename in the workind directory."
        return os.path.join(self.workdir, filename)

    def setup(self, *args, **kwargs):
        """
        Method called before running the calculations. 
        The default implementation is empty.
        """

    def show_inputs(self, stream=sys.stdout):
        lines = []
        app = lines.append

        width = 120
        for task in self:
            app("\n")
            app(repr(task))
            app("\ninput: %s" % task.input_file.path)
            app("\n")
            app(str(task.input))
            app(width*"=" + "\n")

        stream.write("\n".join(lines))

    def register_input(self, input, links=()):
        """
        Registers a new Input:

            - creates a new AbinitTask from input.
            - adds the new task to the internal list, taking into account possible dependencies.

        Returns: WorkLink object
        """
        task_id = len(self) + 1

        task_workdir = os.path.join(self.workdir, "task_" + str(task_id))
        
        # Handle possible dependencies.
        if links:
            if not isinstance(links, collections.Iterable): 
                links = [links,]

        # Create the new task (note the factory so that we create subclasses easily).
        task = task_factory(input, task_workdir, self.runmode, task_id=task_id, links=links)

        # Add it to the internal list.
        #if task.id in self._links_dict:
        #    raise ValueError("task is already registered!")

        self._tasks.append(task)
                                         
        if links:
            self._links_dict[task_id].extend(links)
            print("task_id %s neeeds\n %s" % (task_id, [str(l) for l in links]))

        return WorkLink(task)

    def build(self, *args, **kwargs):
        "Creates the top level directory"
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

    def get_status(self, only_highest_rank=False):
        "Get the status of the tasks in self."
        status_list = [task.status for task in self]

        if only_highest_rank:
            return max(status_list)
        else:
            return status_list

    @property
    def processes(self):
        return [task.process for task in self]

    def rmtree(self, *args, **kwargs):
        """
        Remove all calculation files and directories.
                                                                                   
        Keyword arguments:
            force: (False)
                Do not ask confirmation.
            verbose: (0)
                Print message if verbose is not zero.
        """
        if kwargs.pop('verbose', 0):
            print('Removing directory tree: %s' % self.workdir)
                                                                                    
        shutil.rmtree(self.workdir)

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
            with GSR_Reader(task.odata_path_from_ext("_GSR")) as ncdata:
                etotal.append(ncdata.get_value("etotal"))
                                                            
        return etotal

##########################################################################################

class IterativeWork(Work):
    """
    TODO
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, workdir, runmode, input_generator, max_niter=25):
        """
        Args:
            workdir:
            input_generator:
            max_niter:
                Maximum number of iterations. A negative value or zero value 
                is equivalent to having an infinite number of iterations.
        """
        super(IterativeWork, self).__init__(workdir, runmode)

        self.input_generator = input_generator

        self.max_niter = max_niter

    def next_task(self):
        """
        Generate and register a new task

        Return: task object
        """
        try:
            next_input = next(self.input_generator)
        except StopIteration:
            raise StopIteration

        self.register_input(next_input)
        assert len(self) == self.niter

        return self[-1]

    def submit_tasks(self, *args, **kwargs):
        """
        Run the tasks till self.exit_iteration says to exit or the number of iterations exceeds self.max_niter

        Return dictionary with the final results
        """
        self.niter = 1

        while True:
            if self.max_niter > 0 and self.niter > self.max_niter: 
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
    de_low, de_normal, de_high = [a / (1000 * Ha_eV) for a in atols_mev]

    num_ene = len(etotal)
    etotal_inf = etotal[-1]

    ihigh   = check_conv(etotal, de_high, min_numpts=min_numpts)
    inormal = check_conv(etotal, de_normal)
    ilow    = check_conv(etotal, de_low)

    accidx = {"H": ihigh, "N": inormal, "L": ilow}

    table = []
    app = table.append
                                                                                 
    app(["step", "ecut", "etotal", "et-e_inf [meV]", "accuracy",])
    for idx, (ec, et) in enumerate(zip(ecut_list, etotal)):
        line = "%d %.1f %g %.3f" % (idx, ec, et, (et-etotal_inf)* Ha_eV * 1.e+3)
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
        "exit"        : ihigh != -1,
        "etotal"      : list(etotal),
        "ecut_list"   : ecut_list,
        "aug_ratios"  : aug_ratios,
        "low"         : {"ecut": ecut_low, "aug_ratio": aug_ratio_low},
        "normal"      : {"ecut": ecut_normal, "aug_ratio": aug_ratio_normal},
        "high"        : {"ecut": ecut_high, "aug_ratio": aug_ratio_high},
        "pseudo_name" : pseudo.name,
        "pseudo_path" : pseudo.path,
        "atols_mev"   : atols_mev,
        "dojo_level"  : 0,
    }

    return data

def plot_etotal(ecut_list, etotals, aug_ratios, show=True, savefig=None, *args, **kwargs):
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

        show:
            True to show the figure
        savefig:
            'abc.png' or 'abc.eps'* to save the figure to a file.
    """
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
        emev = Ha2meV(etot)
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
    ax.yaxis.set_view_interval(-10, 20)

    ax.set_title("$\Delta$ Etotal Vs Ecut")

    if show:
        plt.show()
                             
    if savefig is not None:
        fig.savefig(savefig)

class PseudoConvergence(Work):

    def __init__(self, workdir, pseudo, ecut_list, atols_mev,
                 runmode       = "sequential",
                 spin_mode     = "polarized", 
                 acell         = (8, 9, 10), 
                 smearing      = "fermi_dirac:0.1 eV",
                ):

        super(PseudoConvergence, self).__init__(workdir, runmode)

        # Temporary object used to build the input files
        generator = PseudoIterativeConvergence(workdir, pseudo, ecut_list, atols_mev, 
                                               spin_mode     = spin_mode, 
                                               acell         = acell, 
                                               smearing      = smearing,
                                               max_niter     = len(ecut_list),
                                              )
        self.atols_mev = atols_mev
        self.pseudo = Pseudo.aspseudo(pseudo)

        self.ecut_list = []
        for ecut in ecut_list:
            input = generator.input_with_ecut(ecut)
            self.ecut_list.append(ecut)
            self.register_input(input)

    def get_results(self, *args, **kwargs):

        # Get the results of the tasks.
        work_results = super(PseudoConvergence, self).get_results()

        data = compute_hints(self.ecut_list, self.read_etotal(), self.atols_mev, self.pseudo)

        plot_etotal(data["ecut_list"], data["etotal"], data["aug_ratios"], 
            show=False, savefig=self.path_in_workdir("etotal.pdf"))

        work_results.update(data)

        if kwargs.get("json_dump", True):
            json_pretty_dump(work_results, self.path_in_workdir("results.json"))

        return work_results

class PseudoIterativeConvergence(IterativeWork):

    def __init__(self, workdir, pseudo, ecut_list_or_slice, atols_mev,
                 runmode       = "sequential",
                 spin_mode     = "polarized", 
                 acell         = (8, 9, 10), 
                 smearing      = "fermi_dirac:0.1 eV",
                 max_niter     = 50,
                ):
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
        self.spin_mode = spin_mode
        self.smearing  = Smearing.assmearing(smearing)
        self.acell     = acell

        if isinstance(ecut_list_or_slice, slice):
            self.ecut_iterator = iterator_from_slice(ecut_list_or_slice)
        else:
            self.ecut_iterator = iter(ecut_list_or_slice)

        # Construct a generator that returns AbinitInput objects.
        def input_generator():
            for ecut in self.ecut_iterator:
                yield self.input_with_ecut(ecut)

        super(PseudoIterativeConvergence, self).__init__(workdir, runmode, input_generator(), max_niter=max_niter)

        if not self.isnc:
            raise NotImplementedError("PAW convergence tests are not supported yet")

    def input_with_ecut(self, ecut):
        "Return an Input instance with given cutoff energy ecut"

        # Define the system: one atom in a box of lenghts acell. 
        boxed_atom = SystemCard.boxed_atom(self.pseudo, acell=self.acell)

        # Gamma-only sampling.
        gamma_only = KpointsCard.gamma_only()

        # Setup electrons.
        electrons = ElectronsCard(spin_mode=self.spin_mode, smearing=self.smearing)

        # Generate inputs with different values of ecut and register the task.
        control = ControlCard(boxed_atom, electrons, gamma_only,
                              ecut        = ecut, 
                              prtwf       = 0,
                              want_forces = False,
                              want_stress = False,
                             )

        return Input(control, boxed_atom, gamma_only, electrons)

    @property
    def ecut_list(self):
        "The list of cutoff energies computed so far"
        return [float(task.input.get_variable("ecut")) for task in self]

    def check_etotal_convergence(self, *args, **kwargs):
        return compute_hints(self.ecut_list, self.read_etotal(), self.atols_mev, self.pseudo)

    def exit_iteration(self, *args, **kwargs):
        return self.check_etotal_convergence(self, *args, **kwargs)
                                                                                       
    def get_results(self, *args, **kwargs):
        # Get the results of the tasks.
        work_results = super(PseudoIterativeConvergence, self).get_results()

        data = self.check_etotal_convergence()

        plot_etotal(data["ecut_list"], data["etotal"], data["aug_ratios"], 
            show=False, savefig=self.path_in_workdir("etotal.pdf"))

        work_results.update(data)

        if kwargs.get("json_dump", True):
            json_pretty_dump(work_results, self.path_in_workdir("results.json"))

        return work_results

##########################################################################################

class BandStructure(Work):

    def __init__(self, workdir, runmode, structure, pseudos, scf_ngkpt, nscf_nband,
                 spin_mode    = "polarized",
                 scf_smearing = "fermi_dirac:0.1eV",
                 kpath_bounds = None,
                 ndivsm       = 15, 
                 dos_ngkpt    = None
                ):

        super(BandStructure, self).__init__(workdir, runmode)

        scf_smearing = Smearing.assmearing(scf_smearing)

        # Construct the input for the GS-SCF run.
        #scf_input = gs_strategy.make_input()
        scf_input = Input.SCF_groundstate(structure, pseudos, 
                                          ngkpt     = scf_ngkpt, 
                                          spin_mode = spin_mode,
                                          smearing  = scf_smearing,
                                         )

        scf_link = self.register_input(scf_input)

        # Construct the input for the NSCF run.
        #nscf_strategy.learn(scf_input=scf_input)
        #nscf_input = nscf_strategy.make_input()
        nscf_input = Input.NSCF_kpath_from_SCF(scf_input, nscf_nband, ndivsm=ndivsm, kpath_bounds=kpath_bounds)

        self.register_input(nscf_input, links=scf_link.produces_exts("_DEN"))

        # Add DOS computation
        if dos_ngkpt is not None:
            #dos_strategy.learn(scf_input=scf_input)
            #dos_input = dos_strategy.make_input()
            dos_input = Input.NSCF_kmesh_from_SCF(scf_input, nscf_nband, dos_ngkpt)

            self.register_input(dos_input, links=scf_link.produces_exts("_DEN"))

##########################################################################################

class Relaxation(Work):

    def __init__(self, workdir, runmode, structure, pseudos, strategy,
                 ngkpt     = None,
                 kppa      = None,
                 spin_mode = "polarized",
                 smearing  = "fermi_dirac:0.1 eV",
                 **kwargs
                ):
                                                                                                   
        super(Relaxation, self).__init__(workdir, runmode)

        smearing = Smearing.assmearing(smearing)

        #relax_input = strategy.make_input()

        relax_input = Input.Relax(structure, pseudos, strategy,
                                  ngkpt     = ngkpt,
                                  kppa      = kppa,
                                  spin_mode = "polarized", 
                                  smearing  = smearing,
                                  **kwargs
                                 )

        link = self.register_input(relax_input)

##########################################################################################

class DeltaTest(Work):

    def __init__(self, workdir, runmode, structure_or_cif, pseudos, kppa,
                 spin_mode = "polarized",
                 smearing  = "fermi_dirac:0.1 eV",
                 accuracy  = "normal",
                 ecutsm    = 0.05,
                ):

        super(DeltaTest, self).__init__(workdir, runmode)

        if isinstance(structure_or_cif, Structure):
            structure = structure_or_cif
        else:
            # Assume CIF file
            structure = read_structure(structure_or_cif)

        smearing = Smearing.assmearing(smearing)

        self._input_structure = structure

        v0 = structure.volume

        self.volumes = v0 * np.arange(90, 112, 2) / 100.

        for vol in self.volumes:

            new_lattice = structure.lattice.scale(vol)

            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)

            scf_input = Input.SCF_groundstate(new_structure, pseudos, 
                                              kppa      = kppa,
                                              spin_mode = spin_mode,
                                              smearing  = smearing,
                                              # **kwargs
                                              accuracy  = accuracy,
                                              ecutsm    = ecutsm,
                                              prtwf     = 0,
                                             )
            self.register_input(scf_input)

    def get_results(self, *args, **kwargs):

        num_sites = self._input_structure.num_sites

        etotal = Ha2eV(self.read_etotal())

        from .eos import EOS
        eos_fit = EOS.Murnaghan().fit(self.volumes, etotal)

        print(eos_fit)

        eos_fit.plot(show=False, savefig=self.path_in_workdir("eos.pdf"))

        results = {
            "etotal"     : list(etotal),
            "volumes"    : list(self.volumes),
            "natom"      : num_sites,
            "v0"         : eos_fit.v0,
            "b"          : eos_fit.b,
            "bp"         : eos_fit.bp,
            "dojo_level" : 1,
        }

        if kwargs.get("json_dump", True):
            json_pretty_dump(results, self.path_in_workdir("results.json"))

        # Write data for the computation of the delta factor
        with open(self.path_in_workdir("deltadata.txt"), "w") as fh:
            fh.write("# Volume/natom [Ang^3] Etotal/natom [eV]\n")
            for (v, e) in zip(self.volumes, etotal):
                fh.write("%s %s\n" % (v/num_sites, e/num_sites))

##########################################################################################

class G0W0(Work):

    def __init__(self, workdir, runmode, structure, pseudos, scf_ngkpt, nscf_ngkpt, scr_strategy, sigma_strategy,
                 spin_mode = "polarized", smearing  = "fermi_dirac:0.1 eV"):
        """
            Args:
                workdir:
                    Working directory of the calculation.
                runmode:
                structure: 
                    pymatgen structure.
                pseudos: 
                    List of pseudopotentials
                # FIXME
                scf_ngkpt: 
                nscf_ngkpt: 

                scr_strategy: 
                    Strategy for the screening run.
                sigma_strategy:
                    Strategy for the self-energy run.
                spin_mode: 
                    Spin polarization.
                smearing:
                    Smearing technique
        """

        super(G0W0, self).__init__(workdir, runmode)

        smearing = Smearing.assmearing(smearing)

        # Construct the input for the GS-SCF run.
        #scf_input = gs_strategy.make_input()
        scf_input = Input.SCF_groundstate(structure, pseudos, 
                                          ngkpt     = scf_ngkpt, 
                                          spin_mode = spin_mode, 
                                          smearing  = smearing,
                                         )

        scf_link = self.register_input(scf_input)

        scr_nband = scr_strategy.get_varvalue("nband")

        sigma_nband = sigma_strategy.get_varvalue("nband")

        max_nband = max(scr_nband, sigma_nband) 

        nscf_nband = int(max_nband + 0.05 * max_nband)

        istwfk = "*1" # FIXME

        # Construct the input for the NSCF run.
        #nscf_strategy.learn(scf_input=scf_input)
        #nscf_input = nscf_strategy.make_input()
        nscf_input = Input.NSCF_kmesh_from_SCF(scf_input, nscf_nband, nscf_ngkpt, istwfk=istwfk)

        nscf_link = self.register_input(nscf_input, links=scf_link.produces_exts("_DEN"))

        # Construct the input for the SCR run.
        #scr_strategy.learn(scf_input=scf_input, nscf_input=nscf_input)
        #screen_input = scr_strategy.make_input()
        screen_input = Input.SCR_from_NSCF(nscf_input, scr_strategy, smearing=smearing, istwfk=istwfk)

        screen_link = self.register_input(screen_input, links=nscf_link.produces_exts("_WFK"))

        # Construct the input for the SIGMA run.
        #sigma_strategy.learn(scf_input=scf_input, nscf_input=nscf_input, scr_input=scr_input)
        #sigma_input = sigma_strategy.make_input()
        sigma_input = Input.SIGMA_from_SCR(screen_input, sigma_strategy, smearing=smearing, istwfk=istwfk)

        sigma_links = [nscf_link.produces_exts("_WFK"), screen_link.produces_exts("_SCR"),]

        self.register_input(sigma_input, links=sigma_links)

##########################################################################################

class PPConvergenceFactory(object):
    "Factory object"

    def work_for_pseudo(self, workdir, pseudo, ecut_range, 
                        runmode   = "sequential",
                        atols_mev = (10, 1, 0.1),
                        spin_mode = "polarized", 
                        acell     = (8, 9, 10),
                        smearing  = "fermi_dirac:0.1 eV",
                       ):
        """
        Return a Work object given the pseudopotential pseudo

        Args:
            workdir: 
                Working directory.
            pseudo: 
                Pseudo object.
            ecut_range:
                range of cutoff energies in Ha units.
            runmode: 
            atols_mev:
                Tolerances in meV for accuracy in ["low", "normal", "high"]
            spin_mode:
                Spin polarization.
            acell: 
                Length of the real space lattice (Bohr units)
            smearing: 
                Defines the smearing technique.
        """
        workdir = os.path.abspath(workdir)

        smearing = Smearing.assmearing(smearing)

        if isinstance(ecut_range, slice):

            work = PseudoIterativeConvergence(workdir, pseudo, ecut_range, atols_mev,
                                              runmode    = runmode,
                                              spin_mode  = spin_mode,
                                              acell      = acell, 
                                              smearing   = smearing,
                                             )

        else:

            work = PseudoConvergence(workdir, pseudo, ecut_range, atols_mev,
                                     runmode    = runmode,
                                     spin_mode  = spin_mode,
                                     acell      = acell, 
                                     smearing   = smearing,
                                    )

        return work

##########################################################################################
