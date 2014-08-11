"""
Classes defining Abinit calculations and workflows
"""
from __future__ import division, print_function

import os
import time
import shutil
import collections
import abc
import warnings
import copy
import yaml
import numpy as np
from pymatgen.io.abinitio import myaml

from pymatgen.io.abinitio import abiinspect
from pymatgen.io.abinitio import events 

try:
    from pydispatch import dispatcher
except ImportError:
    pass

from monty.json import loadf 
from pymatgen.core.design_patterns import Enum, AttrDict
from pymatgen.util.io_utils import FileLock
from pymatgen.util.string_utils import stream_has_colours, is_string, list_strings, WildCard
from pymatgen.serializers.json_coders import MSONable, json_pretty_dump
from pymatgen.io.abinitio.utils import File, Directory, irdvars_for_ext, abi_splitext, abi_extensions, FilepathFixer, Condition

from pymatgen.io.abinitio.qadapters import qadapter_class
from pymatgen.io.abinitio.netcdf import ETSF_Reader
from pymatgen.io.abinitio.strategies import StrategyWithInput, OpticInput, AnaddbInput, order_pseudos

import logging
logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

__all__ = [
    "TaskManager",
    "ParalHintsParser",
    "ScfTask",
    "NscfTask",
    "RelaxTask",
    "DDK_Task",
    "PhononTask",
    "G_Task",
    "HaydockBseTask",
    "OpticTask",
    "AnaddbTask",
]

# Tools and helper functions.
def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()

class TaskResults(dict, MSONable):
    """
    Dictionary used to store the most important results produced by a Task.
    """
    _MANDATORY_KEYS = [
        "task_name",
        "task_returncode",
        "task_status",
        #"task_events",
    ]

    EXC_KEY = "_exceptions"

    def __init__(self, *args, **kwargs):
        super(TaskResults, self).__init__(*args, **kwargs)
                                                               
        if self.EXC_KEY not in self:
            self[self.EXC_KEY] = []

    @property
    def exceptions(self):
        return self[self.EXC_KEY]

    def push_exceptions(self, *exceptions):
        for exc in exceptions:
            newstr = str(exc)
            if newstr not in self.exceptions:
                self[self.EXC_KEY] += [newstr,]

    def assert_valid(self):
        """
        Returns an empty list if results seem valid. 

        The try assert except trick allows one to get a string with info on the exception.
        We use the += operator so that sub-classes can add their own message.
        """
        # TODO Better treatment of events.
        try:
            assert (self["task_returncode"] == 0 and self["task_status"] == self.S_OK)

        except AssertionError as exc:
            self.push_exceptions(exc)

        return self.exceptions

    @property
    def to_dict(self):
        d = {k: v for k,v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d
                                                                                
    @classmethod
    def from_dict(cls, d):
        return cls({k: v for k,v in d.items() if k not in ["@module", "@class",]})

    def json_dump(self, filename):
        json_pretty_dump(self.to_dict, filename) 

    @classmethod
    def json_load(cls, filename):
        return cls.from_dict(loadf(filename))


class ParalHintsError(Exception):
    """Base error class for `ParalHints`."""


class ParalConf(AttrDict):
    """
    This object store the parameters associated to one 
    of the possible parallel configurations reported by ABINIT.
    Essentially it is a dictionary whose values can also be accessed 
    as attributes. It also provides default values for selected keys
    that might not be present in the ABINIT dictionary.

    Example:

        --- !Autoparal
        info: 
            version: 1
            autoparal: 1
            max_ncpus: 108
        configurations:
            -   tot_ncpus: 2         # Total number of CPUs
                mpi_ncpus: 2         # Number of MPI processes.
                omp_ncpus: 1         # Number of OMP threads (1 if not present)
                mem_per_cpu: 10     # Estimated memory requirement per MPI processor in Megabytes.
                efficiency: 0.4      # 1.0 corresponds to an "expected" optimal efficiency (strong scaling).
                vars: {              # Dictionary with the variables that should be added to the input.
                      varname1: varvalue1
                      varname2: varvalue2
                      }
            -
        ...

    For paral_kgb we have:
    nproc     npkpt  npspinor    npband     npfft    bandpp    weight   
       108       1         1        12         9         2        0.25
       108       1         1       108         1         2       27.00
        96       1         1        24         4         1        1.50
        84       1         1        12         7         2        0.25
    """
    _DEFAULTS = {
        "omp_ncpus": 1,     
        "mem_per_cpu": 0.0, 
        "vars": {}       
    }

    def __init__(self, *args, **kwargs):
        super(ParalConf, self).__init__(*args, **kwargs)
        
        # Add default values if not already in self.
        for k, v in self._DEFAULTS.items():
            if k not in self:
                self[k] = v

    @property
    def speedup(self):
        """Estimated speedup reported by ABINIT."""
        return self.efficiency * self.tot_ncpus

    @property
    def tot_mem(self):
        """Estimated total memory in Mbs (computed from mem_per_cpu)"""
        return self.mem_per_cpu * self.tot_ncpus



class ParalHintsParser(object):

    Error = ParalHintsError

    def parse(self, filename):
        """
        Read the AutoParal section (YAML forma) from filename.
        Assumes the file contains only one section.
        """
        with abiinspect.YamlTokenizer(filename) as r:
            doc = r.next_doc_with_tag("!Autoparal")
            try:
                d = myaml.load(doc.text_notag)
                return ParalHints(info=d["info"], confs=d["configurations"])

            except:
                import traceback
                print("traceback", traceback.format_exc())
                raise self.Error("Wrong YAML doc:\n %s" % doc.text)


class ParalHints(collections.Iterable):
    """
    Iterable with the hints for the parallel execution reported by ABINIT.
    """
    Error = ParalHintsError

    def __init__(self, info, confs):
        self.info = info
        self._confs = [ParalConf(**d) for d in confs]

    def __getitem__(self, key):
        return self._confs[key]

    def __iter__(self):
        return self._confs.__iter__()

    def __len__(self):
        return self._confs.__len__()

    def __str__(self):
        return "\n".join(str(conf) for conf in self)

    def copy(self):
        return copy.copy(self)

    def select_with_condition(self, condition):
        """
        Remove all the configurations that do not satisfy the given condition.

            Args:
                `Condition` object with operators expressed with a Mongodb-like syntax
        """
        new_confs = []

        for conf in self:
            add_it = condition.apply(obj=conf)
            #print("add_it", add_it, "conf", conf)
            if add_it:
                new_confs.append(conf)

        self._confs = new_confs

    def sort_by_efficiency(self, reverse=False):
        """
        Sort the configurations in place so that conf with lowest efficieny 
        appears in the first positions.
        """
        self._confs.sort(key=lambda c: c.efficiency, reverse=reverse)

    def sort_by_speedup(self, reverse=False):
        """
        Sort the configurations in place so that conf with lowest speedup 
        appears in the first positions.
        """
        self._confs.sort(key=lambda c: c.speedup, reverse=reverse)

    def sort_by_mem_cpu(self, reverse=False):
        """
        Sort the configurations in place so that conf with lowest memory per CPU
        appears in the first positions.
        """
        # Avoid sorting if mem_per_cpu is not available.
        has_mem_info = any(c.mem_per_cpu > 0.0 for c in self)

        if has_mem_info:
            self._confs.sort(key=lambda c: c.mem_per_cpu, reverse=reverse)

    def select_optimal_conf(self, policy):
        """
        Find the optimal configuration according to the `TaskPolicy` policy.
        """
        # Make a copy since we are gonna change the object in place.
        #hints = self.copy()

        hints = ParalHints(self.info, confs=[c for c in self if c.tot_ncpus <= policy.max_ncpus])
        #print(hints)

        # First select the configurations satisfying the 
        # condition specified by the user (if any)
        if policy.condition:
            print("condition",policy.condition)
            hints.select_with_condition(policy.condition)
            print("after condition", hints)

            # If no configuration fullfills the requirements, 
            # we return the one with the highest speedup.
            if not hints:
                #print("no configuration")
                hints = self.copy()
                hints.sort_by_speedup()
                return hints[-1].copy()

        hints.sort_by_speedup()

        # Find the optimal configuration according to policy.mode.
        #mode = policy.mode
        #if mode in ["default", "aggressive"]:
        #    hints.sort_by_spedup()

        #elif mode == "conservative":
        #    hints.sort_by_efficiency()
        #    # Remove tot_ncpus == 1
        #    hints.pop(tot_ncpus==1)
        #    if not hints:

        #else:
        #    raise ValueError("Wrong value for mode: %s" % str(mode))

        # Return a copy of the configuration.
        optimal = hints[-1].copy()
        logger.debug("Will relaunch the job with optimized parameters:\n %s" % optimal)

        return optimal


class TaskPolicy(object):
    """
    This object stores the parameters used by the `TaskManager` to 
    create the submission script and/or to modify the ABINIT variables 
    governing the parallel execution. A `TaskPolicy` object contains 
    a set of variables that specify the launcher, as well as the options
    and the condition used to select the optimal configuration for the parallel run 
    """

    def __init__(self, autoparal=0, automemory=0, mode="default", max_ncpus=None, use_fw=False, condition=None): 
        """
        Args:
            autoparal: 
                Value of ABINIT autoparal input variable. None to disable the autoparal feature.
            automemory:
                int defining the memory policy. 
                If > 0 the memory requirements will be computed at run-time from the autoparal section
                produced by ABINIT. In this case, the job script will report the autoparal memory
                instead of the one specified by the user.
            mode:
                Select the algorith to select the optimal configuration for the parallel execution.
                Possible values: ["default", "aggressive", "conservative"]
            max_ncpus:
                Max number of CPUs that can be used (must be specified if autoparal > 0).
            use_fw: 
                True if we are using fireworks.
            condition: 
                condition used to filter the autoparal configuration (Mongodb syntax)
        """
        self.autoparal = autoparal
        self.automemory = automemory
        self.mode = mode 
        self.max_ncpus = max_ncpus
        self.use_fw = use_fw 
        self.condition = Condition(condition) if condition is not None else condition

        if self.autoparal and self.max_ncpus is None:
            raise ValueError("When autoparal is not zero, max_ncpus must be specified.")

    def __str__(self):
        lines = [self.__class__.__name__ + ":"]
        app = lines.append
        for k, v in self.__dict__.items():
            if k.startswith("_"): continue
            app("%s: %s" % (k, v))

        return "\n".join(lines)


class TaskManager(object):
    """
    A `TaskManager` is responsible for the generation of the job script and the submission 
    of the task, as well as of the specification of the parameters passed to the resource manager
    (e.g. Slurm, PBS ...) and/or the run-time specification of the ABINIT variables governing the 
    parallel execution. A `TaskManager` delegates the generation of the submission
    script and the submission of the task to the `QueueAdapter`. 
    A `TaskManager` has a `TaskPolicy` that governs the specification of the 
    parameters for the parallel executions.
    """
    YAML_FILE = "taskmanager.yml"

    def __init__(self, qtype, qparams=None, setup=None, modules=None, shell_env=None, omp_env=None, 
                 pre_run=None, post_run=None, mpi_runner=None, policy=None):

        qad_class = qadapter_class(qtype)

        self.qadapter = qad_class(qparams=qparams, setup=setup, modules=modules, shell_env=shell_env, omp_env=omp_env, 
                                  pre_run=pre_run, post_run=post_run, mpi_runner=mpi_runner)

        if policy is None:
            # Use default policy.
            self.policy = TaskPolicy()
        else:
            if isinstance(policy, TaskPolicy):
                self.policy = policy
            else:
                # Assume dict-like object.
                self.policy = TaskPolicy(**policy) 

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        app("tot_ncpus %d, mpi_ncpus %d, omp_ncpus %s" % (self.tot_ncpus, self.mpi_ncpus, self.omp_ncpus))
        app("MPI_RUNNER %s" % str(self.qadapter.mpi_runner))
        app("policy: %s" % str(self.policy))

        return "\n".join(lines)

    @classmethod
    def from_dict(cls, d):
        """Create an instance from a dictionary."""
        return cls(**d)

    @classmethod
    def from_file(cls, filename):
        """Read the configuration parameters from a Yaml file."""
        with open(filename, "r") as fh:
            return cls.from_dict(myaml.load(fh))

    @classmethod
    def from_user_config(cls):
        """
        Initialize the `TaskManager` from the YAML file 'taskmanager.yaml'.
        Search first in the working directory and then in the configuration
        directory of abipy.

        Raises:
            RuntimeError if file is not found.
        """
        # Try in the current directory.
        path = os.path.join(os.getcwd(), cls.YAML_FILE)

        if os.path.exists(path):
            return cls.from_file(path)

        # Try in the configuration directory.
        dirpath = os.path.join(os.getenv("HOME"), ".abinit", "abipy")
        path = os.path.join(dirpath, cls.YAML_FILE)

        if os.path.exists(path):
            return cls.from_file(path)
    
        raise RuntimeError("Cannot locate %s neither in current directory nor in %s" % (cls.YAML_FILE, dirpath))

    @classmethod 
    def sequential(cls):
        """
        Build a simple `TaskManager` that submits jobs via a simple shell script.
        Assume the shell environment has been already initialized.
        """
        return cls(qtype="shell")

    @classmethod 
    def simple_mpi(cls, mpi_runner="mpirun", mpi_ncpus=1, policy=None):
        """
        Build a `TaskManager` that submits jobs with a simple shell script and mpirun.
        Assume the shell environment is already properly initialized.
        """
        return cls(qtype="shell", qparams=dict(MPI_NCPUS=mpi_ncpus), mpi_runner=mpi_runner, policy=policy)

    @property
    def tot_ncpus(self):
        """Total number of CPUs used to run the task."""
        return self.qadapter.tot_ncpus

    @property
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""
        return self.qadapter.mpi_ncpus

    @property
    def omp_ncpus(self):
        """Number of CPUs used for OpenMP."""
        return self.qadapter.omp_ncpus

    def to_shell_manager(self, mpi_ncpus=1, policy=None):
        """
        Returns a new `TaskManager` with the same parameters as self but replace the `QueueAdapter` 
        with a `ShellAdapter` with mpi_ncpus so that we can submit the job without passing through the queue.
        Replace self.policy with a `TaskPolicy` with autoparal==0.
        """
        cls = self.__class__
        qad = self.qadapter.deepcopy()

        policy = TaskPolicy(autoparal=0) if policy is None else policy

        new = cls("shell", qparams={"MPI_NCPUS": mpi_ncpus}, setup=qad.setup, modules=qad.modules, 
                  shell_env=qad.shell_env, omp_env=qad.omp_env, pre_run=qad.pre_run, 
                  post_run=qad.post_run, mpi_runner=qad.mpi_runner, policy=policy)

        return new

    def new_with_policy(self, policy):
        """
        Returns a new `TaskManager` with same parameters as self except for policy.
        """
        new = self.deepcopy()
        new.policy = policy
        return new

    def copy(self):
        """Shallow copy of self."""
        return copy.copy(self)

    def deepcopy(self):
        """Deep copy of self."""
        return copy.deepcopy(self)

    def set_mpi_ncpus(self, mpi_ncpus):
        """Set the number of MPI nodes to use."""
        self.qadapter.set_mpi_ncpus(mpi_ncpus)

    def set_omp_ncpus(self, omp_ncpus):
        """Set the number of OpenMp threads to use."""
        self.qadapter.set_omp_ncpus(omp_ncpus)

    def set_mem_per_cpu(self, mem_mb):
        """Set the memory (in Megabytes) per CPU."""
        self.qadapter.set_mem_per_cpu(mem_mb)

    def write_jobfile(self, task):
        """
        Write the submission script.

        Args:
            task:
                `AbinitTask` object.

        Returns:
            The path of the script file.
        """
        # Construct the submission script.
        script = self.qadapter.get_script_str(
            job_name=task.name, 
            launch_dir=task.workdir, 
            executable=task.executable,
            qout_path=task.qout_file.path,
            qerr_path=task.qerr_file.path,
            stdin=task.files_file.path, 
            stdout=task.log_file.path,
            stderr=task.stderr_file.path,
            )

        # Write the script.
        script_file = task.job_file.path

        with open(script_file, "w") as fh:
            fh.write(script)

        return script_file

    def launch(self, task):
        """
        Build the input files and submit the task via the `Qadapter` 

        Args:
            task:
                `TaskObject`
        
        Returns:
            Process object.
        """
        task.build()

        script_file = self.write_jobfile(task)

        # Submit the task.
        task.set_status(task.S_SUB)

        # FIXME: CD to script file dir?
        process, queue_id = self.qadapter.submit_to_queue(script_file)

        # Save the queue id.
        task.set_queue_id(queue_id)

        return process


# The code below initializes a counter from a file when the module is imported 
# and save the counter's updated value automatically when the program terminates 
# without relying on the application making an explicit call into this module at termination.
conf_dir = os.path.join(os.getenv("HOME"), ".abinit", "abipy")

if not os.path.exists(conf_dir):
    os.makedirs(conf_dir)

_COUNTER_FILE = os.path.join(conf_dir, "nodecounter")
del conf_dir

try:
    with open(_COUNTER_FILE, "r") as fh:
        _COUNTER = int(fh.read())

except IOError:
    _COUNTER = -1

def get_newnode_id():
    """
    Returns a new node identifier used both for `Task` and `Workflow` objects.

    .. warnings:
        The id is unique inside the same python process so be careful when 
        Workflows and Task are constructed at run-time or when threads are used.
    """
    global _COUNTER
    _COUNTER += 1
    return _COUNTER


def save_lastnode_id():
    """Save the id of the last node created."""
    with FileLock(_COUNTER_FILE) as lock:
        with open(_COUNTER_FILE, "w") as fh:
            fh.write("%d" % _COUNTER)

import atexit
atexit.register(save_lastnode_id)


class FakeProcess(object):
    """
    This object is attached to a Task instance if the task has not been submitted
    This trick allows us to simulate a process that is still running so that 
    we can safely poll task.process.
    """
    def poll(self):
        return None

    def wait(self):
        raise RuntimeError("Cannot wait a FakeProcess")

    def communicate(self, input=None):
        raise RuntimeError("Cannot communicate with a FakeProcess")

    def kill(self):
        raise RuntimeError("Cannot kill a FakeProcess")

    @property
    def returncode(self):
        return None


class Product(object):
    """
    A product represents an output file produced by ABINIT instance.
    This file is needed to start another `Task` or another `Workflow`.
    """
    def __init__(self, ext, path):
        """
        Args:
            ext:
                ABINIT file extension
            path:
                (asbolute) filepath
        """
        if ext not in abi_extensions():
            raise ValueError("Extension %s has not been registered in the internal database" % str(ext))

        self.ext = ext
        self.file = File(path)

    @classmethod
    def from_file(cls, filepath):
        """Build a `Product` instance from a filepath."""
        # Find the abinit extension.
        for i in range(len(filepath)):
            if filepath[i:] in abi_extensions():
                ext = filepath[i:]
                break
        else:
            raise ValueError("Cannot detect abinit extension in %s" % filepath)
        
        return cls(ext, filepath)

    def __str__(self):
        return "File=%s, Extension=%s, " % (self.file.path, self.ext)

    @property
    def filepath(self):
        """Absolute path of the file."""
        return self.file.path

    def connecting_vars(self):
        """
        Returns a dictionary with the ABINIT variables that 
        must be used to make the code use this file.
        """
        return irdvars_for_ext(self.ext)


class Dependency(object):
    """
    This object describes the dependencies among the nodes of a calculation.

    A `Dependency` consists of a `Node` that produces a list of products (files) 
    that are used by the other nodes (`Task` or `Workflow`) to start the calculation.
    One usually creates the object by calling work.register 

    Example:

        # Register the SCF task in work.
        scf_task = work.register(scf_strategy)

        # Register the NSCF calculation and its dependency on the SCF run via deps.
        nscf_task = work.register(nscf_strategy, deps={scf_task: "DEN"})
    """
    def __init__(self, node, exts=None):
        """
        Args:
            node:
                The task or the worfklow associated to the dependency.
            exts:
                Extensions of the output files that are needed for running the other tasks.
        """
        self._node = node

        if exts and is_string(exts):
            exts = exts.split()

        self.exts = exts or []

    def __hash__(self):
        return hash(self._node)

    def __repr__(self):
        return "Node %s will produce: %s " % (repr(self.node), repr(self.exts))

    def __str__(self):
        return "Node %s will produce: %s " % (str(self.node), str(self.exts))

    @property
    def info(self):
        return str(self.node)

    @property
    def node(self):
        """The node associated to the dependency."""
        return self._node

    @property
    def status(self):
        """The status of the dependency, i.e. the status of the node."""
        return self.node.status

    @property
    def products(self):
        """List of output files produces by self."""
        try:
            return self._products
        except:
            self._products = []
            for ext in self.exts:
                prod = Product(ext, self.node.opath_from_ext(ext))
                self._products.append(prod)

            return self._products

    def connecting_vars(self):
        """
        Returns a dictionary with the variables that must be added to the 
        input file in order to connect this `Node` to its dependencies.
        """
        vars = {}
        for prod in self.products:
            vars.update(prod.connecting_vars())

        return vars

    def get_filepaths_and_exts(self):
        """Returns the paths of the output files produced by self and its extensions"""
        filepaths = [prod.filepath for prod in self.products]
        exts = [prod.ext for prod in self.products]

        return filepaths, exts


# Possible status of the node.
STATUS2STR = collections.OrderedDict([
    (1, "Initialized"),   # Node has been initialized
    (2, "Locked"),        # Task is locked an must be explicitly unlocked by en external subject (Workflow).
    (3, "Ready"),         # Node is ready i.e. all the depencies of the node have status S_OK
    (4, "Submitted"),     # Node has been submitted (The `Task` is running or we have started to finalize the Workflow)
    (5, "Running"),       # Node is running.
    (6, "Done"),          # Node done, This does not imply that results are ok or that the calculation completed successfully
    (7, "Error"),         # Node raised some kind of Error (the submission process, the queue manager or ABINIT ...).
    (8, "Unconverged"),   # This usually means that an iterative algorithm didn't converge.
    (9, "Completed"),     # Execution completed successfully.
])

class Status(int):
    """This object is an integer representing the status of the `Node`."""
    def __repr__(self):
        return "<%s: %s, at %s>" % (self.__class__.__name__, str(self), id(self))

    def __str__(self):
        """String representation."""
        return STATUS2STR[self]

    @classmethod
    def from_string(cls, s):
        """Return a `Status` instance from its string representation."""
        for num, text in STATUS2STR.items():
            if text == s:
                return cls(num)
        else:
            raise ValueError("Wrong string %s" % s)


class Node(object):
    """
    Abstract base class defining the interface that must be 
    implemented by the nodes of the calculation.

    Nodes are hashable and can be tested for equality
    (hash uses the node identifier, while eq uses workdir).
    """
    __metaclass__ = abc.ABCMeta

    # Possible status of the node.
    S_INIT = Status(1)
    S_LOCKED = Status(2)
    S_READY = Status(3)
    S_SUB = Status(4)
    S_RUN = Status(5)
    S_DONE = Status(6)
    S_ERROR = Status(7)
    S_UNCONVERGED = Status(8)
    S_OK = Status(9)

    ALL_STATUS = [
        S_INIT,
        S_LOCKED,
        S_READY,
        S_SUB,
        S_RUN,
        S_DONE,
        S_ERROR,
        S_UNCONVERGED,
        S_OK,
    ]

    def __init__(self):
        # Node identifier.
        self._node_id = get_newnode_id()

        # List of dependencies
        self._deps = []

        # List of files (products) needed by this node.
        self._required_files = []

        # Used to push additional info during the execution. 
        self.history = collections.deque(maxlen=100)

        # Set to true if the node has been finalized.
        self._finalized = False

    def __eq__(self, other):
        if not isinstance(other, Node):
            return False

        return (self.__class__ == other.__class__ and 
                self.workdir == other.workdir)
               #self.node_id == other.node_id and 
                                                       
    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.node_id)

    def __repr__(self):
        try:
            return "<%s, node_id %s, workdir=%s>" % (
                self.__class__.__name__, self.node_id, os.path.relpath(self.workdir))

        except AttributeError:
            # this usually happens when workdir has not been initialized
            return "<%s, node_id %s, workdir=None>" % (self.__class__.__name__, self.node_id)
                                                                                            
    def __str__(self):
        try:
            return "<%s, workdir=%s>" % (self.__class__.__name__, os.path.relpath(self.workdir))
        except AttributeError:
            # this usually happens when workdir has not been initialized
            return "<%s, workdir=None>" % self.__class__.__name__

    @property
    def name(self):
        """
        The name of the node 
        (only used for facilitating its identification in the user interface).
        """
        try:
            return self._name
        except AttributeError:
            return os.path.relpath(self.workdir)

    def set_name(self, name):
        """Set the name of the Node."""
        self._name = name

    @property
    def node_id(self):
        """Node identifier."""
        return self._node_id
                                                         
    def set_node_id(self, node_id):
        """Set the node identifier. Use it carefully!"""
        self._node_id = node_id

    @property
    def finalized(self):
        """True if the `Workflow` has been finalized."""
        return self._finalized

    @finalized.setter
    def finalized(self, boolean):
        self._finalized = boolean
        self.history.append("Finalized on %s" % time.asctime())

    @property
    def str_history(self):
        """String representation of history."""
        return "\n".join(self.history)
                                                         
    #@abc.abstractproperty
    #def workdir(self):

    @property
    def has_subnodes(self):
        """True if self contains sub-nodes e.g. `Workflow` object."""
        return isinstance(self, collections.Iterable)

    @property
    def deps(self):
        """
        List of `Dependency` objects defining the dependencies 
        of this `Node`. Empty list if this `Node` does not have dependencies.
        """
        return self._deps

    def add_deps(self, deps):
        """
        Add a list of dependencies to the `Node`.

        Args:
            deps:
                List of `Dependency` objects specifying the 
                dependencies of the node.
        """
        # We want a list
        if not isinstance(deps, (list, tuple)):
            deps = [deps]

        assert all(isinstance(d, Dependency) for d in deps)

        # Add the dependencies to the node
        self._deps.extend(deps)

        if self.has_subnodes:
            # This means that the node contains sub-nodes 
            # that should inherit the same dependency.
            for task in self:
                task.add_deps(deps)

    def remove_deps(self, deps):
        """
        Remove a list of dependencies from the `Node`.

        Args:
            deps:
                List of `Dependency` objects specifying the 
                dependencies of the node.
        """
        if not isinstance(deps, (list, tuple)):
            deps = [deps]
                                                                                      
        assert all(isinstance(d, Dependency) for d in deps)

        self._deps = [d for d in self._deps if d not in deps]
                                                                                      
        if self.has_subnodes:
            # This means that the node consists of sub-nodes 
            # that should remove the same list of dependencies.
            for task in self:
                task.remove_deps(deps)                                                                                                                                        

    @property
    def deps_status(self):
        """Returns a list with the status of the dependencies."""
        if not self.deps:
            return [self.S_OK]
                                                                  
        return [d.status for d in self.deps]

    def depends_on(self, other):
        """True if this node depends on the other node."""
        return other in [d.node for d in self.deps]

    def str_deps(self):
        """Return the string representation of the dependecies of the node."""
        lines = []
        app = lines.append

        app("Dependencies of node %s:" % str(self))
        for i, dep in enumerate(self.deps):
            app("%d) %s, status=%s" % (i, dep.info, str(dep.status)))

        return "\n".join(lines)

    @property
    def required_files(self):
        """
        List of `Product` objects with info on the files needed by the `Node`.
        """
        return self._required_files

    def add_required_files(self, files):
        """
        Add a list of path to the list of files required by the `Node`.

        Args:
            files:
                string or list of strings with the path of the files
        """
        # We want a list of absolute paths.
        files = map(os.path.abspath, list_strings(files))

        # Convert to list of products.
        files = [Product.from_file(path) for path in files]

        # Add the dependencies to the node.
        self._required_files.extend(files)

    #@abc.abstractmethod
    #def set_status(self, status):
    #    """Set the status of the `Node`."""
    #
    #def status(self):
    #    """Return the status of the `Node`."""
    #   return self._status

    #@abc.abstractmethod
    #def check_status(self, status):
    #    """Check the status of the `Node`."""

    #@abc.abstractmethod
    #def connect_signals():
    #    """Connect the signals."""


class TaskError(Exception):
    """Base Exception for `Task` methods"""


class TaskRestartError(TaskError):
    """Exception raised while trying to restart the `Task`."""


class Task(Node):
    __metaclass__ = abc.ABCMeta

    Error = TaskError

    # List of `AbinitEvent` subclasses that are tested in the not_converged method. 
    # Subclasses should provide their own list if they need to check the converge status.
    CRITICAL_EVENTS = [
    ]

    # Prefixes for Abinit (input, output, temporary) files.
    Prefix = collections.namedtuple("Prefix", "idata odata tdata")
    pj = os.path.join

    prefix = Prefix(pj("indata", "in"), pj("outdata", "out"), pj("tmpdata", "tmp"))
    del Prefix, pj

    def __init__(self, strategy, workdir=None, manager=None, deps=None, required_files=None):
        """
        Args:
            strategy: 
                Input file or `Strategy` instance defining the calculation.
            workdir:
                Path to the working directory.
            manager:
                `TaskManager` object.
            deps:
                Dictionary specifying the dependency of this node.
                None means that this obj has no dependency.
            required_files:
                List of strings with the path of the files used by the task.
        """
        # Init the node
        super(Task, self).__init__()

        # Save the strategy to use to generate the input file.
        # FIXME
        #self.strategy = strategy.deepcopy()
        self.strategy = strategy
                                                               
        if workdir is not None:
            self.set_workdir(workdir)
                                                               
        if manager is not None:
            self.set_manager(manager)

        # Handle possible dependencies.
        if deps:
            deps = [Dependency(node, exts) for (node, exts) in deps.items()]
            self.add_deps(deps)

        if required_files:
            self.add_required_files(required_files)

        # Set the initial status.
        self.set_status(self.S_INIT)

        # Number of restarts effectuated.
        self.num_restarts = 0

    def __getstate__(self):
        """
        Return state is pickled as the contents for the instance.
                                                                                      
        In this case we just remove the process since Subprocess objects cannot be pickled.
        This is the reason why we have to store the returncode in self._returncode instead
        of using self.process.returncode.
        """
        return {k:v for k,v in self.__dict__.items() if k not in ["_process",]}

    def set_workdir(self, workdir, chroot=False):
        """Set the working directory. Cannot be set more than once unless chroot is True"""
        if not chroot and hasattr(self, "workdir") and self.workdir != workdir:
                raise ValueError("self.workdir != workdir: %s, %s" % (self.workdir,  workdir))

        self.workdir = os.path.abspath(workdir)

        # Files required for the execution.
        self.input_file = File(os.path.join(self.workdir, "run.abi"))
        self.output_file = File(os.path.join(self.workdir, "run.abo"))
        self.files_file = File(os.path.join(self.workdir, "run.files"))
        self.job_file = File(os.path.join(self.workdir, "job.sh"))
        self.log_file = File(os.path.join(self.workdir, "run.log"))
        self.stderr_file = File(os.path.join(self.workdir, "run.err"))
        self.start_lockfile = File(os.path.join(self.workdir, "__startlock__"))

        # Directories with input|output|temporary data.
        self.indir = Directory(os.path.join(self.workdir, "indata"))
        self.outdir = Directory(os.path.join(self.workdir, "outdata"))
        self.tmpdir = Directory(os.path.join(self.workdir, "tmpdata"))

        # stderr and output file of the queue manager.
        self.qerr_file = File(os.path.join(self.workdir, "queue.err"))
        self.qout_file = File(os.path.join(self.workdir, "queue.out"))

    def set_manager(self, manager):
        """Set the `TaskManager` to use to launch the Task."""
        self.manager = manager.deepcopy()

    @property
    def flow(self):
        """The flow containing this `Task`."""
        return self._flow

    def set_flow(self, flow):
        """Set the flow associated to this `Task`."""
        if not hasattr(self, "_flow"):
            self._flow = flow
        else: 
            if self._flow != flow:
                raise ValueError("self._flow != flow")

    def make_input(self):
        """Construct and write the input file of the calculation."""
        return self.strategy.make_input()

    def ipath_from_ext(self, ext):
        """
        Returns the path of the input file with extension ext.
        Use it when the file does not exist yet.
        """
        return os.path.join(self.workdir, self.prefix.idata + "_" + ext)

    def opath_from_ext(self, ext):
        """
        Returns the path of the output file with extension ext.
        Use it when the file does not exist yet.
        """
        return os.path.join(self.workdir, self.prefix.odata + "_" + ext)

    @abc.abstractproperty
    def executable(self):
        """
        Path to the executable associated to the task (internally stored in self._executable).
        """

    def set_executable(self, executable):
        """Set the executable associate to this task."""
        self._executable = executable

    @property
    def process(self):
        try:
            return self._process
        except AttributeError:
            # Attach a fake process so that we can poll it.
            return FakeProcess()

    #@property
    #def is_allocated(self):
    #    """
    #    True if the task has been allocated, 
    #    i.e. if it has been submitted or if it's running.
    #    """
    #    return self.status in [self.S_SUB, self.S_RUN]

    @property
    def is_completed(self):
        """True if the task has been executed."""
        return self.status >= self.S_DONE

    @property
    def can_run(self):
        """The task can run if its status is < S_SUB and all the other depencies (if any) are done!"""
        all_ok = all([stat == self.S_OK for stat in self.deps_status])
        return self.status < self.S_SUB and all_ok

    def not_converged(self):
        """Return True if the calculation is not converged."""
        logger.debug("not_converged method of the base class will always return False")
        report = self.get_event_report()
        return report.filter_types(self.CRITICAL_EVENTS)

    def cancel(self):
        """
        Cancel the job. Returns 1 if job was cancelled.
        """
        if self.queue_id is None: return 0 
        if self.status >= self.S_DONE: return 0 

        exit_status = self.manager.qadapter.cancel(self.queue_id)
        if exit_status != 0: return 0

        # Remove output files and reset the status.
        self.reset()
        return 1

    def _on_done(self):
        self.fix_ofiles()

    def _on_ok(self):
        self.fix_ofiles()
        results = self.on_ok()
        self._finalized = True
        return results

    def on_ok(self):
        """
        This method is called once the `Task` has reached status S_OK. 
        Subclasses should provide their own implementation

        Returns:
            Dictionary that must contain at least the following entries:
                returncode:
                    0 on success. 
                message: 
                    a string that should provide a human-readable description of what has been performed.
        """
        return dict(returncode=0, 
                    message="Calling on_all_ok of the base class!")

    def fix_ofiles(self):
        """
        This method is called when the task reaches S_OK.
        It changes the extension of particular output files
        produced by Abinit so that the 'official' extension
        is preserved e.g. out_1WF14 --> out_1WF
        """
        filepaths = self.outdir.list_filepaths()
        logger.info("in fix_ofiles with filepaths %s" % filepaths) 

        old2new = FilepathFixer().fix_paths(filepaths)

        for old, new in old2new.items():
            logger.debug("will rename old %s to new %s" % (old, new))
            os.rename(old, new)

    def _restart(self):
        """
        Called by restart once we have finished preparing the task for restarting.
        """
        self.set_status(self.S_READY, info_msg="Restarted on %s" % time.asctime())

        # Increase the counter.
        self.num_restarts += 1
        self.history.append("Restarted on %s, num_restarts %d" % (time.asctime(), self.num_restarts))

        # Remove the lock file
        self.start_lockfile.remove()
 
        # Relaunch the task.
        fired = self.start()

        if not fired:
            self.history.append("[%s], restart failed" % time.asctime())

        return fired

    def restart(self):
        """
        Restart the calculation.  Subclasses should provide a concrete version that 
        performs all the actions needed for preparing the restart and then calls self._restart
        to restart the task. The default implementation is empty.

        Returns:
            1 if job was restarted, 0 otherwise.
        """
        logger.debug("Calling the **empty** restart method of the base class")
        return 0

    def poll(self):
        """Check if child process has terminated. Set and return returncode attribute."""
        self._returncode = self.process.poll()

        if self._returncode is not None:
            self.set_status(self.S_DONE)

        return self._returncode

    def wait(self):
        """Wait for child process to terminate. Set and return returncode attribute."""
        self._returncode = self.process.wait()
        self.set_status(self.S_DONE)

        return self._returncode

    def communicate(self, input=None):
        """
        Interact with process: Send data to stdin. Read data from stdout and stderr, until end-of-file is reached. 
        Wait for process to terminate. The optional input argument should be a string to be sent to the 
        child process, or None, if no data should be sent to the child.

        communicate() returns a tuple (stdoutdata, stderrdata).
        """
        stdoutdata, stderrdata = self.process.communicate(input=input)
        self._returncode = self.process.returncode
        self.set_status(self.S_DONE)

        return stdoutdata, stderrdata 

    def kill(self):
        """Kill the child."""
        self.process.kill()
        self._returncode = self.process.returncode

    @property
    def returncode(self):
        """
        The child return code, set by poll() and wait() (and indirectly by communicate()). 
        A None value indicates that the process hasn't terminated yet.
        A negative value -N indicates that the child was terminated by signal N (Unix only).
        """
        try: 
            return self._returncode
        except AttributeError:
            return 0

    def reset(self):
        """
        Reset the task status. Mainly used if we made a silly mistake in the initial
        setup of the queue manager and we want to fix it and rerun the task.

        Returns:
            0 on success, 1 if reset failed.
        """
        # Can only reset tasks that are done.
        if self.status < self.S_DONE: return 1

        self.set_status(self.S_INIT, info_msg="Reset on %s" % time.asctime())
        self.set_queue_id(None)

        # Remove output files otherwise the EventParser will think the job is still running
        self.output_file.remove()
        self.log_file.remove()
        self.stderr_file.remove()
        self.start_lockfile.remove()
        self.qerr_file.remove()
        self.qout_file.remove()

        # TODO send a signal to the flow 
        #self.workflow.check_status()
        return 0

    @property
    def queue_id(self):
        """Queue identifier returned by the Queue manager. None if not set"""
        try:
            return self._queue_id
        except AttributeError:
            return None

    def set_queue_id(self, queue_id):
        """Set the task identifier."""
        self._queue_id = queue_id

    @property
    def has_queue_manager(self):
        """True if we are submitting jobs via a queue manager."""
        return self.manager.qadapter.QTYPE.lower() != "shell"

    @property
    def tot_ncpus(self):
        """Total number of CPUs used to run the task."""
        return self.manager.tot_ncpus
                                                         
    @property
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""
        return self.manager.mpi_ncpus
                                                         
    @property
    def omp_ncpus(self):
        """Number of CPUs used for OpenMP."""
        return self.manager.omp_ncpus

    @property
    def status(self):
        """Gives the status of the task."""
        return self._status

    def set_status(self, status, info_msg=None):
        """Set the status of the task."""
        # Accepts strings as well.
        if not isinstance(status, Status):
            try:
                status = getattr(Node, status)
            except AttributeError:
                status = Status.from_string(status)

        assert status in STATUS2STR

        changed = True
        if hasattr(self, "_status"):
            changed = (status != self._status)

        self._status = status

        # Add new entry to history only if the status has changed.
        if changed:
            if status == self.S_SUB: 
                self._submission_time = time.time()
                self.history.append("Submitted on %s" % time.asctime())

            if status == self.S_OK:
                self.history.append("Completed on %s" % time.asctime())

            if status == self.S_ERROR:
                self.history.append("Error info:\n %s" % str(info_msg))

        if status == self.S_DONE:
            self._on_done()
                                                                                
        if status == self.S_OK:
            #if status == self.S_UNCONVERGED:
            #    logger.debug("Task %s broadcasts signal S_UNCONVERGED" % self)
            #    dispatcher.send(signal=self.S_UNCONVERGED, sender=self)
                                                                                
            # Finalize the task.
            if not self.finalized:
                self._on_ok()
                                                                                
            logger.debug("Task %s broadcasts signal S_OK" % self)
            dispatcher.send(signal=self.S_OK, sender=self)

        return status

    def check_status(self):
        """
        This function check the status of the task by inspecting the output and the 
        error files produced by the application and by the queue manager.
        """
        # A locked task can only be unlocked by calling set_status explicitly.
        black_list = [self.S_LOCKED]
        if self.status in black_list: return

        # Check the returncode of the process first.
        if self.returncode != 0:
            return self.set_status(self.S_ERROR, info_msg="return code %s" % self.returncode)

        # Start to check when the output file has been created.
        if not self.output_file.exists:
            logger.debug("output_file does not exists")

            if not self.stderr_file.exists and not self.qerr_file.exists:
                # The job is still in the queue.
                return self.status

            else:
                # Analyze the standard error of the executable:
                if self.stderr_file.exists:
                    err_msg = self.stderr_file.read()
                    if err_msg:
                        logger.critical("%s: executable stderr:\n %s" % (self, err_msg))
                        return self.set_status(self.S_ERROR, info_msg=err_msg)

                # Analyze the error file of the resource manager.
                if self.qerr_file.exists:
                    err_info = self.qerr_file.read()
                    if err_info:
                        logger.critical("%s: queue stderr:\n %s" % (self, err_msg))
                        return self.set_status(self.S_ERROR, info_msg=err_info)

                return self.status

        # Check if the run completed successfully.
        report = self.get_event_report()

        if report.run_completed:
            # Check if the calculation converged.
            not_ok = self.not_converged()

            if not_ok:
                return self.set_status(self.S_UNCONVERGED)
            else:
                return self.set_status(self.S_OK)

        # This is the delicate part since we have to discern among different possibilities:
        #
        # 1) Calculation stopped due to an Abinit Error or Bug.
        #
        # 2) Segmentation fault that (by definition) was not handled by ABINIT.
        #    In this case we check if the ABINIT standard error is not empty.
        #    hoping that nobody has written to sdterr (e.g. libraries in debug mode)
        #
        # 3) Problem with the resource manager and/or the OS (walltime error, resource error, phase of the moon ...)
        #    In this case we check if the error file of the queue manager is not empty.
        #    Also in this case we *assume* that there's something wrong if the stderr of the queue manager is not empty
        # 
        # 4) Calculation is still running!
        #
        # Point 2) and 3) are the most complicated since there's no standard!

        # 1) Search for possible errors or bugs in the ABINIT **output** file.
        if report.errors or report.bugs:
            logger.critical("%s: Found Errors or Bugs in ABINIT main output!" % self)
            return self.set_status(self.S_ERROR, info_msg=str(report.errors) + str(report.bugs))

        # 2) Analyze the stderr file for Fortran runtime errors.
        if self.stderr_file.exists:
            err_info = self.stderr_file.read()
            if err_info:
                return self.set_status(self.S_ERROR, info_msg=err_info)

        # 3) Analyze the error file of the resource manager.
        if self.qerr_file.exists:
            err_info = self.qerr_file.read()
            if err_info:
                return self.set_status(self.S_ERROR, info_msg=err_info)

        # 4) Assume the job is still running.
        return self.set_status(self.S_RUN)

    def out_to_in(self, out_file):
        """
        Move an output file to the output data directory of the `Task` 
        and rename the file so that ABINIT will read it as an input data file.

        Returns:
            The absolute path of the new file in the indata directory.
        """
        in_file = os.path.basename(out_file).replace("out", "in", 1)
        dest = os.path.join(self.indir.path, in_file)
                                                                           
        if os.path.exists(dest) and not os.path.islink(dest):
           logger.warning("Will overwrite %s with %s" % (dest, out_file))
                                                                           
        os.rename(out_file, dest)
        return dest

    def inlink_file(self, filepath):
        """
        Create a symbolic link to the specified file in the 
        directory containing the input files of the task.
        """
        if not os.path.exists(filepath): 
            logger.debug("Creating symbolic link to not existent file %s" % filepath)

        # Extract the Abinit extension and add the prefix for input files.
        root, abiext = abi_splitext(filepath)

        infile = "in_" + abiext
        infile = self.indir.path_in(infile)

        # Link path to dest if dest link does not exist.
        # else check that it points to the expected file.
        logger.debug("Linking path %s --> %s" % (filepath, infile))
        print("Linking path %s --> %s" % (filepath, infile))
                                                             
        if not os.path.exists(infile):
            os.symlink(filepath, infile)
        else:
            if os.path.realpath(infile) != filepath:
                raise self.Error("infile %s does not point to filepath %s" % (infile, filepath))

    def make_links(self):
        """
        Create symbolic links to the output files produced by the other tasks.

        ..warning:
            
            This method should be called only when the calculation is READY because
            it uses a heuristic approach to find the file to link.
        """
        for dep in self.deps:
            filepaths, exts = dep.get_filepaths_and_exts()

            for (path, ext) in zip(filepaths, exts):
                logger.info("Need path %s with ext %s" % (path, ext))
                dest = self.ipath_from_ext(ext)

                if not os.path.exists(path): 
                    # Try netcdf file. TODO: this case should be treated in a cleaner way.
                    path += "-etsf.nc"
                    if os.path.exists(path): dest += "-etsf.nc"

                if not os.path.exists(path):
                    err_msg = "%s: %s is needed by this task but it does not exist" % (self, path)
                    logger.critical(err_msg)
                    raise self.Error(err_msg)

                # Link path to dest if dest link does not exist.
                # else check that it points to the expected file.
                logger.debug("Linking path %s --> %s" % (path, dest))

                if not os.path.exists(dest):
                    os.symlink(path, dest)
                else:
                    if os.path.realpath(dest) != path:
                        raise self.Error("dest %s does not point to path %s" % (dest, path))

        for f in self.required_files:
            path, dest = f.filepath, self.ipath_from_ext(f.ext)
      
            # Link path to dest if dest link does not exist.
            # else check that it points to the expected file.
            print("Linking path %s --> %s" % (path, dest))
                                                                                         
            if not os.path.exists(dest):
                os.symlink(path, dest)
            else:
                if os.path.realpath(dest) != path:
                    raise self.Error("dest %s does not point to path %s" % (dest, path))


    @abc.abstractmethod
    def setup(self):
        """Public method called before submitting the task."""

    def _setup(self):
        """
        This method calls self.setup after having performed additional operations
        such as the creation of the symbolic links needed to connect different tasks.
        """
        self.make_links()
        self.setup()

    # TODO: For the time being, we inspect the log file,
    # We will start to use the output file when the migration to YAML is completed
    def get_event_report(self, source="log"):
        """
        Analyzes the main output file for possible Errors or Warnings.

        Args:
            source:
                "output" for the main output file.
                "log" for the log file.

        Returns:
            `EventReport` instance or None if the main output file does not exist.
        """
        file = {
            "output": self.output_file,
            "log": self.log_file,
        }[source]

        if not file.exists:
            return None

        parser = events.EventsParser()
        try:
            return parser.parse(file.path)

        except parser.Error as exc:
            # Return a report with an error entry with info on the exception.
            logger.critical("%s: Exception while parsing ABINIT events:\n %s" % (file, str(exc)))
            self.set_status(self.S_ERROR, info_msg=str(exc))
            return parser.report_exception(file.path, exc)

    @property
    def results(self):
        """The results produced by the task. Set by get_results"""
        try:
            return self._results

        except AttributeError:
            self._results = self.get_results()
            return self._results 

    def get_results(self, *args, **kwargs):
        """
        Method called once the calculation is completed, 
        Updates self._results and returns TaskResults instance.
        Subclasses should extend this method (if needed) by adding 
        specialized code that performs some kind of post-processing.
        """
        # Check whether the process completed.
        if self.returncode is None:
            raise self.Error("return code is None, you should call wait, communitate or poll")

        if self.status is None or self.status < self.S_DONE:
            raise self.Error("Task is not completed")

        return TaskResults({
            "task_name"      : self.name,
            "task_returncode": self.returncode,
            "task_status"    : self.status,
            #"task_events"    : self.events.to_dict
        })

    def move(self, dest, is_abspath=False):
        """
        Recursively move self.workdir to another location. This is similar to the Unix "mv" command.
        The destination path must not already exist. If the destination already exists
        but is not a directory, it may be overwritten depending on os.rename() semantics.

        Be default, dest is located in the parent directory of self.workdir.
        Use is_abspath=True to specify an absolute path.
        """
        if not is_abspath:
            dest = os.path.join(os.path.dirname(self.workdir), dest)

        shutil.move(self.workdir, dest)

    def in_files(self):
        """Return all the input data files used."""
        return self.indir.list_filepaths()

    def out_files(self):
        """Return all the output data files produced."""
        return self.outdir.list_filepaths()

    def tmp_files(self):
        """Return all the input data files produced."""
        return self.tmpdir.list_filepaths()

    def path_in_workdir(self, filename):
        """Create the absolute path of filename in the top-level working directory."""
        return os.path.join(self.workdir, filename)

    def rename(self, src_basename, dest_basename, datadir="outdir"):
        """
        Rename a file located in datadir.

        src_basename and dest_basename are the basename of the source file
        and of the destination file, respectively.
        """
        directory = {
            "indir": self.indir,
            "outdir": self.outdir,
            "tmpdir": self.tmpdir,
        }[datadir]

        src = directory.path_in(src_basename)
        dest = directory.path_in(dest_basename)

        os.rename(src, dest)

    def build(self, *args, **kwargs):
        """
        Creates the working directory and the input files of the `Task`.
        It does not overwrite files if they already exist.
        """
        # Create dirs for input, output and tmp data.
        self.indir.makedirs()
        self.outdir.makedirs()
        self.tmpdir.makedirs()

        # Write files file and input file.
        if not self.files_file.exists:
            self.files_file.write(self.filesfile_string)

        self.input_file.write(self.make_input())

        self.manager.write_jobfile(self)

    def rmtree(self, exclude_wildcard=""):
        """
        Remove all files and directories in the working directory

        Args:
            exclude_wildcard:
                Optional string with regular expressions separated by |.
                Files matching one of the regular expressions will be preserved.
                example: exclude_wildcard="*.nc|*.txt" preserves all the files
                whose extension is in ["nc", "txt"].
        """
        if not exclude_wildcard:
            shutil.rmtree(self.workdir)

        else:
            w = WildCard(exclude_wildcards)

            for dirpath, dirnames, filenames in os.walk(self.workdir):
                for fname in filenames:
                    filepath = os.path.join(dirpath, fname)
                    if not w.match(fname):
                        os.remove(filepath)

    def remove_files(self, *filenames):
        """Remove all the files listed in filenames."""
        filenames = list_strings(filenames)

        for dirpath, dirnames, fnames in os.walk(self.workdir):
            for fname in fnames:
                if fname in filenames:
                    filepath = os.path.join(dirpath, fname)
                    os.remove(filepath)

    def setup(self):
        """Base class does not provide any hook."""

    def start(self):
        """
        Starts the calculation by performing the following steps:

            - build dirs and files
            - call the _setup method
            - execute the job file by executing/submitting the job script.

        Returns:
            1 if task was started, 0 otherwise.
            
        """
        if self.status >= self.S_SUB:
            raise self.Error("Task status: %s" % str(self.status))

        if self.start_lockfile.exists:
            logger.critical("Found lock file: %s" % self.start_lockfile.relpath)
            return 0

        self.start_lockfile.write("Started on %s" % time.asctime())

        self.build()
        self._setup()

        # Add the variables needed to connect the node.
        for d in self.deps:
            vars = d.connecting_vars()
            logger.debug("Adding connecting vars %s " % vars)
            self.strategy.add_extra_abivars(vars)

        # Add the variables needed to read the required files
        for f in self.required_files:
            #raise NotImplementedError("")
            vars = irdvars_for_ext("DEN")
            print("Adding connecting vars %s " % vars)
            self.strategy.add_extra_abivars(vars)

        # Automatic parallelization
        if hasattr(self, "autoparal_fake_run"):
            try:
                self.autoparal_fake_run()
            except:
                # Log the exception and continue with the parameters specified by the user.
                logger.critical("autoparal_fake_run raised:\n%s" % straceback())
                self.set_status(self.S_ERROR)
                return 0

        # Start the calculation in a subprocess and return.
        self._process = self.manager.launch(self)

        return 1

    #def start_and_wait(self, *args, **kwargs):
    #    """
    #    Helper method to start the task and wait.

    #    Mainly used when we are submitting the task via the shell
    #    without passing through a queue manager.
    #    """
    #    self.start(*args, **kwargs)
    #    return self.wait()


class AbinitTask(Task):
    """
    Base class defining an ABINIT calculation
    """
    #IN = "in"
    #OUT = "out"
    #TMP = "tmp"

    @classmethod
    def from_input(cls, abinit_input, workdir=None, manager=None):
        """
        Create an instance of `AbinitTask` from an ABINIT input.
    
        Args:
            abinit_input:
                `AbinitInput` object.
            workdir:
                Path to the working directory.
            manager:
                `TaskManager` object.
        """
        # TODO: Find a better way to do this. I will likely need to refactor the Strategy object
        strategy = StrategyWithInput(abinit_input)

        return cls(strategy, workdir=workdir, manager=manager)

    def setup(self):
        """
        Abinit has the very *bad* habit of changing the file extension by appending the characters in [A,B ..., Z] 
        to the output file, and this breaks a lot of code that relies of the use of a unique file extension.
        Here we fix this issue by renaming run.abo to run.abo_[number] if the output file "run.abo" already
        exists. A few lines of code in python, a lot of problems if you try to implement this trick in Fortran90. 
        """
        # God rot the FORTRAN committee who was not able to give Fortran a decent standard library as well as $Windows$ OS!
        # I don't really care if Fortran2003 provides support for OOP programming. 
        # What I need is a standardized interface to communicate with the OS!
        if self.output_file.exists:
            # Find the index of the last file (if any) and push.
            # TODO: Maybe it's better to use run.abo --> run(1).abo
            fnames = [f for f in os.listdir(self.workdir) if f.startswith(self.output_file.basename)]
            nums = [int(f) for f in [f.split("_")[-1] for f in fnames] if f.isdigit()]
            last = max(nums) if nums else 0
            new_path = self.output_file.path + "_" + str(last+1)

            # Call os.rename and are done! It's really amazing the that I can write Fortran code that runs on 10**3 processors 
            # whereas a simple mv in F90 requires a lot of unportable tricks (where unportable means "not supported by $windows$").
            print("Will rename %s to %s" % (self.output_file.path, new_path))
            os.rename(self.output_file.path, new_path)

    @property
    def executable(self):
        """Path to the executable required for running the Task."""
        try:
            return self._executable
        except AttributeError:
            return "abinit"

    @property
    def pseudos(self):
        """List of pseudos used in the calculation."""
        return self.strategy.pseudos

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return all(p.isnc for p in self.pseudos)

    @property
    def ispaw(self):
        """True if PAW calculation"""
        return all(p.ispaw for p in self.pseudos)

    @property
    def filesfile_string(self):
        """String with the list of files and prefixex needed to execute ABINIT."""
        lines = []
        app = lines.append
        pj = os.path.join

        app(self.input_file.path)                 # Path to the input file
        app(self.output_file.path)                # Path to the output file
        app(pj(self.workdir, self.prefix.idata))  # Prefix for input data
        app(pj(self.workdir, self.prefix.odata))  # Prefix for output data
        app(pj(self.workdir, self.prefix.tdata))  # Prefix for temporary data

        # Paths to the pseudopotential files.
        # Note that here the pseudos **must** be sorted according to znucl.
        for pseudo in self.pseudos:
            app(pseudo.path)

        return "\n".join(lines)

    def autoparal_fake_run(self):
        """
        Find an optimal set of parameters for the execution of the task 
        using the options specified in `TaskPolicy`.
        This method can change the ABINIT input variables and/or the 
        parameters passed to the `TaskManager` e.g. the number of CPUs for MPI and OpenMp.

        Returns:
           confs, optimal 
           where confs is a `ParalHints` object with the configuration reported by 
           autoparal and optimal is the optimal configuration selected.
           Returns (None, None) if some problem occurred.
        """
        logger.info("in autoparal_fake_run")
        manager = self.manager
        policy = manager.policy

        if policy.autoparal == 0 or policy.max_ncpus in [None, 1]: 
            msg = "Nothing to do in autoparal, returning (None, None)"
            logger.info(msg)
            print(msg)
            return None, None

        if policy.autoparal != 1:
            raise NotImplementedError("autoparal != 1")

        # 1) Run ABINIT in sequential to get the possible configurations with max_ncpus

        # Set the variables for automatic parallelization
        autoparal_vars = dict(
            autoparal=policy.autoparal,
            max_ncpus=policy.max_ncpus,
        )

        self.strategy.add_extra_abivars(autoparal_vars)

        # Build a simple manager to run the job in a shell subprocess on the frontend
        # we don't want to make a request to the queue manager for this simple job!
        seq_manager = manager.to_shell_manager(mpi_ncpus=1)

        # Return code is always != 0 
        process = seq_manager.launch(self)
        #print("launched")
        process.wait()  

        # Remove the variables added for the automatic parallelization
        self.strategy.remove_extra_abivars(autoparal_vars.keys())

        # 2) Parse the autoparal configurations from the main output file.
        #print("parsing")
        parser = ParalHintsParser()

        try:
            confs = parser.parse(self.output_file.path)
            #print("confs", confs)
            #self.all_autoparal_confs = confs 

        except parser.Error:
            logger.critical("Error while parsing Autoparal section:\n%s" % straceback())
            return None, None

        # 3) Select the optimal configuration according to policy
        optimal = confs.select_optimal_conf(policy)
        print("optimal Autoparal conf:\n %s" % optimal)

        # 4) Change the input file and/or the submission script
        self.strategy.add_extra_abivars(optimal.vars)
                                                                  
        # Change the number of MPI nodes.
        manager.set_mpi_ncpus(optimal.mpi_ncpus)

        # Change the number of OpenMP threads.
        #if optimal.omp_ncpus > 1:
        #    manager.set_omp_ncpus(optimal.omp_ncpus)
        #else:
        #    manager.disable_omp()

        # Change the memory per node if automemory evaluates to True.
        mem_per_cpu = optimal.mem_per_cpu

        if policy.automemory and mem_per_cpu:
            # mem_per_cpu = max(mem_per_cpu, policy.automemory)
            manager.set_mem_per_cpu(mem_per_cpu)

        # Reset the status, remove garbage files ...
        self.set_status(self.S_INIT)

        # Remove the output file since Abinit likes to create new files 
        # with extension .outA, .outB if the file already exists.
        os.remove(self.output_file.path)
        os.remove(self.log_file.path)
        os.remove(self.stderr_file.path)

        return confs, optimal


# TODO
# Enable restarting capabilites:
# Before doing so I need:
#   1) Preliminary standardization of the ABINT events and critical WARNINGS (YAML)
#   2) Change the parser so that we can use strings in the input file.
#      We need this change for restarting structural relaxations so that we can read 
#      the initial structure from file.

class ScfTask(AbinitTask):
    """
    Self-consistent ground-state calculations.
    Provide support for in-place restart via (WFK|DEN) files
    """
    CRITICAL_EVENTS = [
        events.ScfConvergenceWarning,
    ]

    def restart(self):
        """SCF calculations can be restarted if we have either the WFK file or the DEN file."""
        # Prefer WFK over DEN files since we can reuse the wavefunctions.
        for ext in ["WFK", "DEN"]:
            restart_file = self.outdir.has_abiext(ext)
            irdvars = irdvars_for_ext(ext)
            if restart_file:
                break

        if not restart_file:
            raise TaskRestartError("Cannot find WFK or DEN file to restart from.")

        # Move out --> in.
        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

        # Now we can resubmit the job.
        return self._restart()

    def inspect(self, **kwargs):
        """
        Plot the SCF cycle results with matplotlib.

        Returns
            `matplotlib` figure, None is some error occurred. 
        """
        scf_cycle = abiinspect.GroundStateScfCycle.from_file(self.output_file.path)
        if scf_cycle is not None:
            return scf_cycle.plot(**kwargs)


class NscfTask(AbinitTask):
    """
    Non-Self-consistent GS calculation.
    Provide in-place restart via WFK files
    """
    CRITICAL_EVENTS = [
        events.NscfConvergenceWarning,
    ]

    def restart(self):
        """NSCF calculations can be restarted only if we have the WFK file."""
        ext = "WFK"
        restart_file = self.outdir.has_abiext(ext)
        irdvars = irdvars_for_ext(ext)

        if not restart_file:
            raise TaskRestartError("Cannot find the WFK file to restart from.")

        # Move out --> in.
        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

        # Now we can resubmit the job.
        return self._restart()


class RelaxTask(AbinitTask):
    """
    Structural optimization.

    .. attributes:

        initial_structure:
        final_structure:
    """
    # What about a possible ScfConvergenceWarning?
    CRITICAL_EVENTS = [
        events.RelaxConvergenceWarning,
    ]

    def change_structure(self, structure):
        """Change the input structure."""
        print("changing structure")
        self.strategy.abinit_input.set_structure(structure)

    def read_final_structure(self, save=False):
        """Read the final structure from the GSR file and save it in self.final_structure."""
        # We already have it in memory.
        if hasattr(self, "final_structure"):
            return self.final_structure
        
        # Read it from file and save it if save is True.
        gsr_file = self.outdir.has_abiext("GSR")
        if not gsr_file:
            raise TaskRestartError("Cannot find the GSR file with the final structure to restart from.")

        with ETSF_Reader(gsr_file) as r:
            final_structure = r.read_structure()

        if save:
            self.final_structure = final_structure

        return final_structure

    def restart(self):
        # Structure relaxations can be restarted only if we have the WFK file or the DEN or the GSR file.
        # from which we can read the last structure (mandatory) and the wavefunctions (not mandatory but useful).
        # Prefer WFK over other files since we can reuse the wavefunctions.
        for ext in ["WFK", "DEN"]:
            ofile = self.outdir.has_abiext(ext)
            if ofile:

                irdvars = irdvars_for_ext(ext)
                infile = self.out_to_in(ofile)
                break

        if not ofile:
            raise TaskRestartError("Cannot find the WFK|DEN file to restart from.")

        # Read the relaxed structure from the GSR file.
        structure = self.read_final_structure()
                                                           
        # Change the structure.
        #self.change_structure(structure)

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

        # Now we can resubmit the job.
        return self._restart()

    def inspect(self, **kwargs):
        """
        Plot the evolution of the structural relaxation with matplotlib.

        Returns
            `matplotlib` figure, None is some error occurred. 
        """
        relaxation = abiinspect.Relaxation.from_file(self.output_file.path)
        if relaxation is not None:
            return relaxation.plot(**kwargs)


class DDK_Task(AbinitTask):
    """Task for DDK calculations."""


class PhononTask(AbinitTask):
    """
    DFPT calculations for a single atomic perturbation.
    Provide support for in-place restart via (1WF|1DEN) files
    """
    # TODO: 
    # for the time being we don't discern between GS and PhononCalculations.
    # Restarting Phonon calculation is more difficult due to the crazy rules employed in ABINIT 
    CRITICAL_EVENTS = [
        events.ScfConvergenceWarning,
    ]

    def restart(self):
        # Phonon calculations can be restarted only if we have the 1WF file or the 1DEN file.
        # from which we can read the first-order wavefunctions or the first order density.
        # Prefer 1WF over 1DEN since we can reuse the wavefunctions.
        #self.fix_ofiles()
        for ext in ["1WF", "1DEN"]:
            restart_file = self.outdir.has_abiext(ext)
            irdvars = irdvars_for_ext(ext)
            if restart_file:
                break

        if not restart_file:
            raise TaskRestartError("Cannot find the 1WF|1DEN|file to restart from.")

        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

        # Now we can resubmit the job.
        return self._restart()

    def inspect(self, **kwargs):
        """
        Plot the Phonon SCF cycle results with matplotlib.

        Returns
            `matplotlib` figure, None is some error occurred. 
        """
        scf_cycle = abiinspect.PhononScfCycle.from_file(self.output_file.path)
        if scf_cycle is not None:
            return scf_cycle.plot(**kwargs)


class G_Task(AbinitTask):
    """
    Tasks for SIGMA calculations employing the self-consistent G approximation 
    Provide support for in-place restart via QPS files
    """
    CRITICAL_EVENTS = [
        events.QPSConvergenceWarning,
    ]

    def restart(self):
        # G calculations can be restarted only if we have the QPS file 
        # from which we can read the results of the previous step.
        ext = "QPS"
        restart_file = self.outdir.has_abiext(ext)
        irdvars = irdvars_for_ext(ext)

        if not restart_file:
            raise TaskRestartError("Cannot find the QPS file to restart from.")

        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

        # Now we can resubmit the job.
        return self._restart()


class BseTask(AbinitTask):
    """
    Task for Bethe-Salpeter calculations.

    .. note:

        The BSE codes provides both iterative and direct schemes
        for the computation of the dielectric function. 
        The direct diagonalization cannot be restarted whereas 
        Haydock and CG support restarting.
    """

#class CgBseTask(BseTask):
#    """Bethe-Salpeter calculations with the conjugate-gradient method."""


class HaydockBseTask(BseTask):
    """
    Bethe-Salpeter calculations with Haydock iterative scheme.
    Provide in-place restart via (BSR|BSC) files
    """
    CRITICAL_EVENTS = [
        events.HaydockConvergenceWarning,
    ]

    def restart(self):
        # BSE calculations with Haydock can be restarted only if we have the 
        # excitonic Hamiltonian and the HAYDR_SAVE file.
        # TODO: This version seems to work but the main output file is truncated
        # the log file is complete though.
        irdvars = {}

        # Move the BSE blocks to indata.
        # This is done only once at the end of the first run.
        # Successive restarts will use the BSR|BSC files in the indir directory
        # to initialize the excitonic Hamiltonian
        count = 0
        for ext in ["BSR", "BSC"]:
            ofile = self.outdir.has_abiext(ext)
            if ofile:
                count += 1
                irdvars.update(irdvars_for_ext(ext))
                self.out_to_in(ofile)

        if not count:
            # outdir does not contain the BSR|BSC file.
            # This means that num_restart > 1 and the files should be in task.indir
            count = 0
            for ext in ["BSR", "BSC"]:
                ifile = self.indir.has_abiext(ext)
                if ifile:
                    count += 1

            if not count:
                raise TaskRestartError("Cannot find BSR|BSC files in %s" % self.indir)

        # Rename HAYDR_SAVE files
        count = 0
        for ext in ["HAYDR_SAVE", "HAYDC_SAVE"]:
            ofile = self.outdir.has_abiext(ext)
            if ofile:
                count += 1
                irdvars.update(irdvars_for_ext(ext))
                self.out_to_in(ofile)

        if not count:
            raise TaskRestartError("Cannot find the HAYDR_SAVE file to restart from.")

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

        # Now we can resubmit the job.
        return self._restart()


class OpticTask(Task):
    # TODO
    # FIx the problem with report.is_completed
    # all the executables in Abinit should signal the successful completion with the same format.
    # possibly with YAML

    def __init__(self, optic_input, nscf_node, ddk_nodes, workdir=None, manager=None):
        """
        Create an instance of `OpticTask` from an string containing the input.
    
        Args:
            optic_input:
                string with the optic variables (filepaths will be added at run time).
            nscf_node:
                The NSCF task that will produce thw WFK file.
            ddk_nodes:
                List of `DDK_Task` nodes that will produce the DDK files.
            workdir:
                Path to the working directory.
            manager:
                `TaskManager` object.
        """
        if is_string(nscf_node):
            assert all(is_string(f) for f in ddk_nodes)
            nscf_node = FileNode(nscf_node)
            ddk_nodes = [FileNode(fname) for fname in ddk_nodes]

        deps = {task: "1WF" for task in ddk_nodes}
        deps.update({nscf_node: "WFK"})
        print("deps",deps)

        strategy = OpticInput(optic_input)
        super(OpticTask, self).__init__(strategy=strategy, workdir=workdir, manager=manager, deps=deps)

        # Keep a reference to the nscf_task and the ddk tasks
        assert len(ddk_nodes) == 3
        self.nscf_node, self.ddk_nodes = nscf_node, ddk_nodes

    def set_workdir(self, workdir):
        super(OpticTask, self).set_workdir(workdir)
        # Small hack: the log file of optics is actually the main output file. 
        self.output_file = self.log_file

    @property
    def executable(self):
        """Path to the executable required for running the `OpticTask`."""
        try:
            return self._executable
        except AttributeError:
            return "optic"

    @property
    def filesfile_string(self):
        """String with the list of files and prefixes needed to execute ABINIT."""
        lines = []
        app = lines.append
        pj = os.path.join

        #optic.in     ! Name of input file
        #optic.out    ! Unused
        #optic        ! Root name for all files that will be produced
        app(self.input_file.path)                 # Path to the input file
        app(pj(self.workdir, "unused"))           # Path to the output file
        app(pj(self.workdir, self.prefix.odata))  # Prefix for output data

        return "\n".join(lines)

    @property
    def wfk_filepath(self):
        """Returns (at runtime) the absolute path of the WFK file produced by the NSCF run."""
        return self.nscf_node.outdir.has_abiext("WFK")

    @property
    def ddk_filepaths(self):
        """Returns (at runtime) the absolute path of the DDK files produced by the DDK runs."""
        return [ddk_task.outdir.has_abiext("1WF") for ddk_task in self.ddk_nodes]

    def make_input(self):
        """Construct and write the input file of the calculation."""
        # Set the file paths.
        files = "\n".join(self.ddk_filepaths + [self.wfk_filepath]) + "\n"

        # Get the input specified by the user
        user_inp = self.strategy.make_input()

        # Join them.
        return files + user_inp

    def setup(self):
        """Public method called before submitting the task."""

    def make_links(self):
        """
        Optic allows the user to specify the paths of the input file.
        hence we don't need to create symbolic links.
        """

class AnaddbTask(Task):
    # TODO
    # FIx the problem with report.is_completed
    # all the executables in Abinit should signal the successful completion with the same format.
    # possibly with YAML
    def __init__(self, anaddb_input, ddb_node, 
                 gkk_node=None, md_node=None, ddk_node=None, workdir=None, manager=None):
        """
        Create an instance of `AnaddbTask` from an string containing the input.

        Args:
            anaddb_input:
                string with the anaddb variables.
            ddb_node:
                The node that will produce the DDB file (can be either `Task` or `Workflow` object)
            gkk_node:
                The node that will produce the GKK file (can be either `Task` or `Workflow` object)
                optional.
            md_node:
                The node that will produce the MD file (can be either `Task` or `Workflow` object)
                optional.
            gkk_node:
                The node that will produce the GKK file (can be either `Task` or `Workflow` object)
                optional.
            workdir:
                Path to the working directory.
            manager:
                `TaskManager` object.
        """
        # Keep a reference to the nodes.
        if is_string(ddb_node): ddb_node = FileNode(ddb_node)
        deps = {ddb_node: "DDB"}
        self.ddb_node = ddb_node

        if is_string(gkk_node): gkk_node = FileNode(gkk_node)
        if gkk_node is not None: deps.update({gkk_node: "GKK"})
        self.gkk_node = gkk_node

        # TODO: I never used it!
        if is_string(md_node): md_node = FileNode(md_node)
        if md_node is not None: deps.update({md_node: "MD"})
        self.md_node = md_node

        if is_string(ddk_node): ddk_node = FileNode(ddk_node)
        if ddk_node is not None: deps.update({ddk_node: "DDK"})
        self.ddk_node = ddk_node
                                                                                                        
        # TODO Refactor this code.
        strategy = AnaddbInput(anaddb_input)
        super(AnaddbTask, self).__init__(strategy=strategy, workdir=workdir, manager=manager, deps=deps)

    @property
    def executable(self):
        """Path to the executable required for running the `AnaddbTask`."""
        try:
            return self._executable
        except AttributeError:
            return "anaddb"

    @property
    def filesfile_string(self):
        """String with the list of files and prefixes needed to execute ABINIT."""
        lines = []
        app = lines.append

        app(self.input_file.path)          # 1) Path of the input file
        app(self.output_file.path)         # 2) Path of the output file
        app(self.ddb_filepath)             # 3) Input derivative database e.g. t13.ddb.in
        app(self.md_filepath)              # 4) Output molecular dynamics e.g. t13.md
        app(self.gkk_filepath)             # 5) Input elphon matrix elements  (GKK file)
        # FIXME check this one
        app(self.outdir.path_join("out"))  # 6) Base name for elphon output files e.g. t13
        app(self.ddk_filepath)             # 7) File containing ddk filenames for elphon/transport.

        return "\n".join(lines)

    @property
    def ddb_filepath(self):
        """Returns (at runtime) the absolute path of the input DDB file."""
        path = self.ddb_node.outdir.has_abiext("DDB")
        return path if path else "DDB_FILE_DOES_NOT_EXIST"

    @property
    def md_filepath(self):
        """Returns (at runtime) the absolute path of the input MD file."""
        if self.md_node is None:
            return "MD_FILE_DOES_NOT_EXIST"

        path = self.md_node.outdir.has_abiext("MD")
        return path if path else "MD_FILE_DOES_NOT_EXIST"

    @property
    def gkk_filepath(self):
        """Returns (at runtime) the absolute path of the input GKK file."""
        if self.gkk_node is None:
            return "GKK_FILE_DOES_NOT_EXIST"

        path = self.gkk_node.outdir.has_abiext("GKK")
        return path if path else "GKK_FILE_DOES_NOT_EXIST"

    @property
    def ddk_filepath(self):
        """Returns (at runtime) the absolute path of the input DKK file."""
        if self.ddk_node is None:
            return "DDK_FILE_DOES_NOT_EXIST"

        path = self.ddk_node.outdir.has_abiext("DDK")
        return path if path else "DDK_FILE_DOES_NOT_EXIST"

    def setup(self):
        """Public method called before submitting the task."""

    def make_links(self):
        """
        Anaddb allows the user to specify the paths of the input file.
        hence we don't need to create symbolic links.
        """
