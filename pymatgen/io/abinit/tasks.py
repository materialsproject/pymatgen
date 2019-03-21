# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""This module provides functions and classes related to Task objects."""

import os
import time
import datetime
import shutil
import collections
import abc
import copy
import ruamel.yaml as yaml
from io import StringIO
import numpy as np

from pprint import pprint
from itertools import product

from monty.string import is_string, list_strings
from monty.termcolor import colored, cprint
from monty.collections import AttrDict
from monty.functools import lazy_property, return_none_if_raise
from monty.json import MSONable
from monty.fnmatch import WildCard
from pymatgen.core.units import Memory
from pymatgen.util.serialization import json_pretty_dump, pmg_serialize
from .utils import File, Directory, irdvars_for_ext, abi_splitext, FilepathFixer, Condition, SparseHistogram
from .qadapters import make_qadapter, QueueAdapter, QueueAdapterError
from . import qutils as qu
from .db import DBConnector
from .nodes import Status, Node, NodeError, NodeResults, NodeCorrections, FileNode, check_spectator
from . import abiinspect
from . import events


__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

__all__ = [
    "TaskManager",
    "AbinitBuild",
    "ParalHintsParser",
    "ParalHints",
    "AbinitTask",
    "ScfTask",
    "NscfTask",
    "RelaxTask",
    "DdkTask",
    "PhononTask",
    "ElasticTask",
    "SigmaTask",
    "EphTask",
    "OpticTask",
    "AnaddbTask",
    "set_user_config_taskmanager",
]

import logging
logger = logging.getLogger(__name__)

# Tools and helper functions.

def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def lennone(PropperOrNone):
    if PropperOrNone is None:
        return 0
    else:
        return len(PropperOrNone)


def nmltostring(nml):
    """Convert a dictionary representing a Fortran namelist into a string."""
    if not isinstance(nml,dict):
      raise ValueError("nml should be a dict !")

    curstr = ""
    for key,group in nml.items():
       namelist = ["&" + key]
       for k, v in group.items():
         if isinstance(v, list) or isinstance(v, tuple):
           namelist.append(k + " = " + ",".join(map(str, v)) + ",")
         elif is_string(v):
           namelist.append(k + " = '" + str(v) + "',")
         else:
           namelist.append(k + " = " + str(v) + ",")
       namelist.append("/")

       curstr = curstr + "\n".join(namelist) + "\n"

    return curstr


class TaskResults(NodeResults):

    JSON_SCHEMA = NodeResults.JSON_SCHEMA.copy()
    JSON_SCHEMA["properties"] = {
        "executable": {"type": "string", "required": True},
    }

    @classmethod
    def from_node(cls, task):
        """Initialize an instance from an :class:`AbinitTask` instance."""
        new = super(TaskResults, cls).from_node(task)

        new.update(
            executable=task.executable,
            #executable_version:
            #task_events=
            pseudos=[p.as_dict() for p in task.input.pseudos],
            #input=task.input
        )

        new.register_gridfs_files(
            run_abi=(task.input_file.path, "t"),
            run_abo=(task.output_file.path, "t"),
        )

        return new


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
                mem_per_cpu: 10      # Estimated memory requirement per MPI processor in Megabytes.
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

    def __str__(self):
        stream = StringIO()
        pprint(self, stream=stream)
        return stream.getvalue()

    @property
    def num_cores(self):
        return self.mpi_procs * self.omp_threads

    @property
    def mem_per_proc(self):
        return self.mem_per_cpu

    @property
    def mpi_procs(self):
        return self.mpi_ncpus

    @property
    def omp_threads(self):
        return self.omp_ncpus

    @property
    def speedup(self):
        """Estimated speedup reported by ABINIT."""
        return self.efficiency * self.num_cores

    @property
    def tot_mem(self):
        """Estimated total memory in Mbs (computed from mem_per_proc)"""
        return self.mem_per_proc * self.mpi_procs


class ParalHintsError(Exception):
    """Base error class for `ParalHints`."""


class ParalHintsParser:

    Error = ParalHintsError

    def __init__(self):
        # Used to push error strings.
        self._errors = collections.deque(maxlen=100)

    def add_error(self, errmsg):
        self._errors.append(errmsg)

    def parse(self, filename):
        """
        Read the `AutoParal` section (YAML format) from filename.
        Assumes the file contains only one section.
        """
        with abiinspect.YamlTokenizer(filename) as r:
            doc = r.next_doc_with_tag("!Autoparal")
            try:
                d = yaml.safe_load(doc.text_notag)
                return ParalHints(info=d["info"], confs=d["configurations"])
            except:
                import traceback
                sexc = traceback.format_exc()
                err_msg = "Wrong YAML doc:\n%s\n\nException:\n%s" % (doc.text, sexc)
                self.add_error(err_msg)
                logger.critical(err_msg)
                raise self.Error(err_msg)


class ParalHints(collections.abc.Iterable):
    """
    Iterable with the hints for the parallel execution reported by ABINIT.
    """
    Error = ParalHintsError

    def __init__(self, info, confs):
        self.info = info
        self._confs = [ParalConf(**d) for d in confs]

    @classmethod
    def from_mpi_omp_lists(cls, mpi_procs, omp_threads):
        """
        Build a list of Parallel configurations from two lists
        containing the number of MPI processes and the number of OpenMP threads
        i.e. product(mpi_procs, omp_threads).
        The configuration have parallel efficiency set to 1.0 and no input variables.
        Mainly used for preparing benchmarks.
        """
        info = {}
        confs = [ParalConf(mpi_ncpus=p, omp_ncpus=p, efficiency=1.0)
                 for p, t in product(mpi_procs, omp_threads)]

        return cls(info, confs)

    def __getitem__(self, key):
        return self._confs[key]

    def __iter__(self):
        return self._confs.__iter__()

    def __len__(self):
        return self._confs.__len__()

    def __repr__(self):
        return "\n".join(str(conf) for conf in self)

    def __str__(self):
        return repr(self)

    @lazy_property
    def max_cores(self):
        """Maximum number of cores."""
        return max(c.mpi_procs * c.omp_threads for c in self)

    @lazy_property
    def max_mem_per_proc(self):
        """Maximum memory per MPI process."""
        return max(c.mem_per_proc for c in self)

    @lazy_property
    def max_speedup(self):
        """Maximum speedup."""
        return max(c.speedup for c in self)

    @lazy_property
    def max_efficiency(self):
        """Maximum parallel efficiency."""
        return max(c.efficiency for c in self)

    @pmg_serialize
    def as_dict(self, **kwargs):
        return {"info": self.info, "confs": self._confs}

    @classmethod
    def from_dict(cls, d):
        return cls(info=d["info"], confs=d["confs"])

    def copy(self):
        """Shallow copy of self."""
        return copy.copy(self)

    def select_with_condition(self, condition, key=None):
        """
        Remove all the configurations that do not satisfy the given condition.

            Args:
                condition: dict or :class:`Condition` object with operators expressed with a Mongodb-like syntax
                key: Selects the sub-dictionary on which condition is applied, e.g. key="vars"
                    if we have to filter the configurations depending on the values in vars
        """
        condition = Condition.as_condition(condition)
        new_confs = []

        for conf in self:
            # Select the object on which condition is applied
            obj = conf if key is None else AttrDict(conf[key])
            add_it = condition(obj=obj)
            #if key is "vars": print("conf", conf, "added:", add_it)
            if add_it: new_confs.append(conf)

        self._confs = new_confs

    def sort_by_efficiency(self, reverse=True):
        """Sort the configurations in place. items with highest efficiency come first"""
        self._confs.sort(key=lambda c: c.efficiency, reverse=reverse)
        return self

    def sort_by_speedup(self, reverse=True):
        """Sort the configurations in place. items with highest speedup come first"""
        self._confs.sort(key=lambda c: c.speedup, reverse=reverse)
        return self

    def sort_by_mem_per_proc(self, reverse=False):
        """Sort the configurations in place. items with lowest memory per proc come first."""
        # Avoid sorting if mem_per_cpu is not available.
        if any(c.mem_per_proc > 0.0 for c in self):
            self._confs.sort(key=lambda c: c.mem_per_proc, reverse=reverse)
        return self

    def multidimensional_optimization(self, priorities=("speedup", "efficiency")):
        # Mapping property --> options passed to sparse_histogram
        opts = dict(speedup=dict(step=1.0), efficiency=dict(step=0.1), mem_per_proc=dict(memory=1024))
        #opts = dict(zip(priorities, bin_widths))

        opt_confs = self._confs
        for priority in priorities:
            histogram = SparseHistogram(opt_confs, key=lambda c: getattr(c, priority), **opts[priority])
            pos = 0 if priority == "mem_per_proc" else -1
            opt_confs = histogram.values[pos]

        #histogram.plot(show=True, savefig="hello.pdf")
        return self.__class__(info=self.info, confs=opt_confs)

    #def histogram_efficiency(self, step=0.1):
    #    """Returns a :class:`SparseHistogram` with configuration grouped by parallel efficiency."""
    #    return SparseHistogram(self._confs, key=lambda c: c.efficiency, step=step)

    #def histogram_speedup(self, step=1.0):
    #    """Returns a :class:`SparseHistogram` with configuration grouped by parallel speedup."""
    #    return SparseHistogram(self._confs, key=lambda c: c.speedup, step=step)

    #def histogram_memory(self, step=1024):
    #    """Returns a :class:`SparseHistogram` with configuration grouped by memory."""
    #    return SparseHistogram(self._confs, key=lambda c: c.speedup, step=step)

    #def filter(self, qadapter):
    #    """Return a new list of configurations that can be executed on the `QueueAdapter` qadapter."""
    #    new_confs = [pconf for pconf in self if qadapter.can_run_pconf(pconf)]
    #    return self.__class__(info=self.info, confs=new_confs)

    def get_ordered_with_policy(self, policy, max_ncpus):
        """
        Sort and return a new list of configurations ordered according to the :class:`TaskPolicy` policy.
        """
        # Build new list since we are gonna change the object in place.
        hints = self.__class__(self.info, confs=[c for c in self if c.num_cores <= max_ncpus])

        # First select the configurations satisfying the condition specified by the user (if any)
        bkp_hints = hints.copy()
        if policy.condition:
            logger.info("Applying condition %s" % str(policy.condition))
            hints.select_with_condition(policy.condition)

            # Undo change if no configuration fullfills the requirements.
            if not hints:
                hints = bkp_hints
                logger.warning("Empty list of configurations after policy.condition")

        # Now filter the configurations depending on the values in vars
        bkp_hints = hints.copy()
        if policy.vars_condition:
            logger.info("Applying vars_condition %s" % str(policy.vars_condition))
            hints.select_with_condition(policy.vars_condition, key="vars")

            # Undo change if no configuration fullfills the requirements.
            if not hints:
                hints = bkp_hints
                logger.warning("Empty list of configurations after policy.vars_condition")

        if len(policy.autoparal_priorities) == 1:
            # Example: hints.sort_by_speedup()
            if policy.autoparal_priorities[0] in ['efficiency', 'speedup', 'mem_per_proc']:
                getattr(hints, "sort_by_" + policy.autoparal_priorities[0])()
            elif isinstance(policy.autoparal_priorities[0], collections.Mapping):
                if policy.autoparal_priorities[0]['meta_priority'] == 'highest_speedup_minimum_efficiency_cutoff':
                    min_efficiency = policy.autoparal_priorities[0].get('minimum_efficiency', 1.0)
                    hints.select_with_condition({'efficiency': {'$gte': min_efficiency}})
                    hints.sort_by_speedup()
        else:
            hints = hints.multidimensional_optimization(priorities=policy.autoparal_priorities)
            if len(hints) == 0: raise ValueError("len(hints) == 0")

        #TODO: make sure that num_cores == 1 is never selected when we have more than one configuration
        #if len(hints) > 1:
        #    hints.select_with_condition(dict(num_cores={"$eq": 1)))

        # Return final (orderded ) list of configurations (best first).
        return hints


class TaskPolicy:
    """
    This object stores the parameters used by the :class:`TaskManager` to
    create the submission script and/or to modify the ABINIT variables
    governing the parallel execution. A `TaskPolicy` object contains
    a set of variables that specify the launcher, as well as the options
    and the conditions used to select the optimal configuration for the parallel run
    """
    @classmethod
    def as_policy(cls, obj):
        """
        Converts an object obj into a `:class:`TaskPolicy. Accepts:

            * None
            * TaskPolicy
            * dict-like object
        """
        if obj is None:
            # Use default policy.
            return TaskPolicy()
        else:
            if isinstance(obj, cls):
                return obj
            elif isinstance(obj, collections.abc.Mapping):
                return cls(**obj)
            else:
                raise TypeError("Don't know how to convert type %s to %s" % (type(obj), cls))

    @classmethod
    def autodoc(cls):
        return """
    autoparal:                # (integer). 0 to disable the autoparal feature (DEFAULT: 1 i.e. autoparal is on)
    condition:                # condition used to filter the autoparal configurations (Mongodb-like syntax).
                              # DEFAULT: empty i.e. ignored.
    vars_condition:           # Condition used to filter the list of ABINIT variables reported by autoparal
                              # (Mongodb-like syntax). DEFAULT: empty i.e. ignored.
    frozen_timeout:           # A job is considered frozen and its status is set to ERROR if no change to
                              # the output file has been done for `frozen_timeout` seconds. Accepts int with seconds or
                              # string in slurm form i.e. days-hours:minutes:seconds. DEFAULT: 1 hour.
    precedence:               # Under development.
    autoparal_priorities:     # Under development.
"""

    def __init__(self, **kwargs):
        """
        See autodoc
        """
        self.autoparal = kwargs.pop("autoparal", 1)
        self.condition = Condition(kwargs.pop("condition", {}))
        self.vars_condition = Condition(kwargs.pop("vars_condition", {}))
        self.precedence = kwargs.pop("precedence", "autoparal_conf")
        self.autoparal_priorities = kwargs.pop("autoparal_priorities", ["speedup"])
        #self.autoparal_priorities = kwargs.pop("autoparal_priorities", ["speedup", "efficiecy", "memory"]
        # TODO frozen_timeout could be computed as a fraction of the timelimit of the qadapter!
        self.frozen_timeout = qu.slurm_parse_timestr(kwargs.pop("frozen_timeout", "0-1:00:00"))

        if kwargs:
            raise ValueError("Found invalid keywords in policy section:\n %s" % str(kwargs.keys()))

        # Consistency check.
        if self.precedence not in ("qadapter", "autoparal_conf"):
            raise ValueError("Wrong value for policy.precedence, should be qadapter or autoparal_conf")

    def __str__(self):
        lines = []
        app = lines.append
        for k, v in self.__dict__.items():
            if k.startswith("_"): continue
            app("%s: %s" % (k, v))
        return "\n".join(lines)


class ManagerIncreaseError(Exception):
    """
    Exception raised by the manager if the increase request failed
    """


class FixQueueCriticalError(Exception):
    """
    error raised when an error could not be fixed at the task level
    """


# Global variable used to store the task manager returned by `from_user_config`.
_USER_CONFIG_TASKMANAGER = None


def set_user_config_taskmanager(manager):
    """Change the default manager returned by TaskManager.from_user_config."""
    global _USER_CONFIG_TASKMANAGER
    _USER_CONFIG_TASKMANAGER  = manager


class TaskManager(MSONable):
    """
    A `TaskManager` is responsible for the generation of the job script and the submission
    of the task, as well as for the specification of the parameters passed to the resource manager
    (e.g. Slurm, PBS ...) and/or the run-time specification of the ABINIT variables governing the parallel execution.
    A `TaskManager` delegates the generation of the submission script and the submission of the task to the :class:`QueueAdapter`.
    A `TaskManager` has a :class:`TaskPolicy` that governs the specification of the parameters for the parallel executions.
    Ideally, the TaskManager should be the **main entry point** used by the task to deal with job submission/optimization
    """
    YAML_FILE = "manager.yml"
    USER_CONFIG_DIR = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")

    ENTRIES = {"policy", "qadapters", "db_connector", "batch_adapter"}

    @classmethod
    def autodoc(cls):
        from .db import DBConnector
        s = """
# TaskManager configuration file (YAML Format)

policy:
    # Dictionary with options used to control the execution of the tasks.

qadapters:
    # List of qadapters objects (mandatory)
    -  # qadapter_1
    -  # qadapter_2

db_connector:
    # Connection to MongoDB database (optional)

batch_adapter:
    # Adapter used to submit flows with batch script. (optional)

##########################################
# Individual entries are documented below:
##########################################

"""
        s += "policy: " + TaskPolicy.autodoc() + "\n"
        s += "qadapter: " + QueueAdapter.autodoc() + "\n"
        #s += "db_connector: " + DBConnector.autodoc()
        return s

    @classmethod
    def from_user_config(cls):
        """
        Initialize the :class:`TaskManager` from the YAML file 'manager.yaml'.
        Search first in the working directory and then in the AbiPy configuration directory.

        Raises:
            RuntimeError if file is not found.
        """
        global _USER_CONFIG_TASKMANAGER
        if _USER_CONFIG_TASKMANAGER is not None:
            return _USER_CONFIG_TASKMANAGER

        # Try in the current directory then in user configuration directory.
        path = os.path.join(os.getcwd(), cls.YAML_FILE)
        if not os.path.exists(path):
            path = os.path.join(cls.USER_CONFIG_DIR, cls.YAML_FILE)

        if not os.path.exists(path):
            raise RuntimeError(colored(
		"\nCannot locate %s neither in current directory nor in %s\n"
                "!!! PLEASE READ THIS: !!!\n"
                "To use AbiPy to run jobs this file must be present\n"
                "It provides a description of the cluster/computer you are running on\n"
                "Examples are provided in abipy/data/managers." % (cls.YAML_FILE, path), color="red"))

        _USER_CONFIG_TASKMANAGER = cls.from_file(path)
        return _USER_CONFIG_TASKMANAGER

    @classmethod
    def from_file(cls, filename):
        """Read the configuration parameters from the Yaml file filename."""
        try:
            with open(filename, "r") as fh:
                return cls.from_dict(yaml.safe_load(fh))
        except Exception as exc:
            print("Error while reading TaskManager parameters from %s\n" % filename)
            raise

    @classmethod
    def from_string(cls, s):
        """Create an instance from string s containing a YAML dictionary."""
        return cls.from_dict(yaml.safe_load(s))

    @classmethod
    def as_manager(cls, obj):
        """
        Convert obj into TaskManager instance. Accepts string, filepath, dictionary, `TaskManager` object.
        If obj is None, the manager is initialized from the user config file.
        """
        if isinstance(obj, cls): return obj
        if obj is None: return cls.from_user_config()

        if is_string(obj):
            if os.path.exists(obj):
                return cls.from_file(obj)
            else:
                return cls.from_string(obj)

        elif isinstance(obj, collections.abc.Mapping):
            return cls.from_dict(obj)
        else:
            raise TypeError("Don't know how to convert type %s to TaskManager" % type(obj))

    @classmethod
    def from_dict(cls, d):
        """Create an instance from a dictionary."""
        return cls(**{k: v for k, v in d.items() if k in cls.ENTRIES})

    @pmg_serialize
    def as_dict(self):
        return copy.deepcopy(self._kwargs)

    def __init__(self, **kwargs):
        """
        Args:
            policy:None
            qadapters:List of qadapters in YAML format
            db_connector:Dictionary with data used to connect to the database (optional)
        """
        # Keep a copy of kwargs
        self._kwargs = copy.deepcopy(kwargs)

        self.policy = TaskPolicy.as_policy(kwargs.pop("policy", None))

        # Initialize database connector (if specified)
        self.db_connector = DBConnector(**kwargs.pop("db_connector", {}))

        # Build list of QAdapters. Neglect entry if priority == 0 or `enabled: no"
        qads = []
        for d in kwargs.pop("qadapters"):
            if d.get("enabled", False): continue
            qad = make_qadapter(**d)
            if qad.priority > 0:
                qads.append(qad)
            elif qad.priority < 0:
                raise ValueError("qadapter cannot have negative priority:\n %s" % qad)

        if not qads:
            raise ValueError("Received emtpy list of qadapters")
        #if len(qads) != 1:
        #    raise NotImplementedError("For the time being multiple qadapters are not supported! Please use one adapter")

        # Order qdapters according to priority.
        qads = sorted(qads, key=lambda q: q.priority)
        priorities = [q.priority for q in qads]
        if len(priorities) != len(set(priorities)):
            raise ValueError("Two or more qadapters have same priority. This is not allowed. Check taskmanager.yml")

        self._qads, self._qadpos = tuple(qads), 0

        # Initialize the qadapter for batch script submission.
        d = kwargs.pop("batch_adapter", None)
        self.batch_adapter = None
        if d: self.batch_adapter = make_qadapter(**d)
        #print("batch_adapter", self.batch_adapter)

        if kwargs:
            raise ValueError("Found invalid keywords in the taskmanager file:\n %s" % str(list(kwargs.keys())))

    @lazy_property
    def abinit_build(self):
        """:class:`AbinitBuild` object with Abinit version and options used to build the code"""
        return AbinitBuild(manager=self)

    def to_shell_manager(self, mpi_procs=1):
        """
        Returns a new `TaskManager` with the same parameters as self but replace the :class:`QueueAdapter`
        with a :class:`ShellAdapter` with mpi_procs so that we can submit the job without passing through the queue.
        """
        my_kwargs = copy.deepcopy(self._kwargs)
        my_kwargs["policy"] = TaskPolicy(autoparal=0)

        # On BlueGene we need at least two qadapters.
        # One for running jobs on the computing nodes and another one
        # for running small jobs on the fronted. These two qadapters
        # will have different enviroments and different executables.
        # If None of the q-adapters has qtype==shell, we change qtype to shell
        # and we return a new Manager for sequential jobs with the same parameters as self.
        # If the list contains a qadapter with qtype == shell, we ignore the remaining qadapters
        # when we build the new Manager.
        has_shell_qad = False
        for d in my_kwargs["qadapters"]:
            if d["queue"]["qtype"] == "shell": has_shell_qad = True
        if has_shell_qad:
            my_kwargs["qadapters"] = [d for d in my_kwargs["qadapters"] if d["queue"]["qtype"] == "shell"]

        for d in my_kwargs["qadapters"]:
            d["queue"]["qtype"] = "shell"
            d["limits"]["min_cores"] = mpi_procs
            d["limits"]["max_cores"] = mpi_procs

            # If shell_runner is specified, replace mpi_runner with shell_runner
            # in the script used to run jobs on the frontend.
            # On same machines based on Slurm, indeed, mpirun/mpiexec is not available
            # and jobs should be executed with `srun -n4 exec` when running on the computing nodes
            # or with `exec` when running in sequential on the frontend.
            if "job" in d and "shell_runner" in d["job"]:
                shell_runner = d["job"]["shell_runner"]
                #print("shell_runner:", shell_runner, type(shell_runner))
                if not shell_runner or shell_runner == "None": shell_runner = ""
                d["job"]["mpi_runner"] = shell_runner
                #print("shell_runner:", shell_runner)

        #print(my_kwargs)
        new = self.__class__(**my_kwargs)
        new.set_mpi_procs(mpi_procs)

        return new

    def new_with_fixed_mpi_omp(self, mpi_procs, omp_threads):
        """
        Return a new `TaskManager` in which autoparal has been disabled.
        The jobs will be executed with `mpi_procs` MPI processes and `omp_threads` OpenMP threads.
        Useful for generating input files for benchmarks.
        """
        new = self.deepcopy()
        new.policy.autoparal = 0
        new.set_mpi_procs(mpi_procs)
        new.set_omp_threads(omp_threads)
        return new

    @property
    def has_queue(self):
        """True if we are submitting jobs via a queue manager."""
        return self.qadapter.QTYPE.lower() != "shell"

    @property
    def qads(self):
        """List of :class:`QueueAdapter` objects sorted according to priorities (highest comes first)"""
        return self._qads

    @property
    def qadapter(self):
        """The qadapter used to submit jobs."""
        return self._qads[self._qadpos]

    def select_qadapter(self, pconfs):
        """
        Given a list of parallel configurations, pconfs, this method select an `optimal` configuration
        according to some criterion as well as the :class:`QueueAdapter` to use.

        Args:
            pconfs: :class:`ParalHints` object with the list of parallel configurations

        Returns:
            :class:`ParallelConf` object with the `optimal` configuration.
        """
        # Order the list of configurations according to policy.
        policy, max_ncpus = self.policy, self.max_cores
        pconfs = pconfs.get_ordered_with_policy(policy, max_ncpus)

        if policy.precedence == "qadapter":

            # Try to run on the qadapter with the highest priority.
            for qadpos, qad in enumerate(self.qads):
                possible_pconfs = [pc for pc in pconfs if qad.can_run_pconf(pc)]

                if qad.allocation == "nodes":
                #if qad.allocation in ["nodes", "force_nodes"]:
                    # Select the configuration divisible by nodes if possible.
                    for pconf in possible_pconfs:
                        if pconf.num_cores % qad.hw.cores_per_node == 0:
                            return self._use_qadpos_pconf(qadpos, pconf)

                # Here we select the first one.
                if possible_pconfs:
                    return self._use_qadpos_pconf(qadpos, possible_pconfs[0])

        elif policy.precedence == "autoparal_conf":
            # Try to run on the first pconf irrespectively of the priority of the qadapter.
            for pconf in pconfs:
                for qadpos, qad in enumerate(self.qads):

                    if qad.allocation == "nodes" and not pconf.num_cores % qad.hw.cores_per_node == 0:
                        continue # Ignore it. not very clean

                    if qad.can_run_pconf(pconf):
                        return self._use_qadpos_pconf(qadpos, pconf)

        else:
            raise ValueError("Wrong value of policy.precedence = %s" % policy.precedence)

        # No qadapter could be found
        raise RuntimeError("Cannot find qadapter for this run!")

    def _use_qadpos_pconf(self, qadpos, pconf):
        """
        This function is called when we have accepted the :class:`ParalConf` pconf.
        Returns pconf
        """
        self._qadpos = qadpos

        # Change the number of MPI/OMP cores.
        self.set_mpi_procs(pconf.mpi_procs)
        if self.has_omp: self.set_omp_threads(pconf.omp_threads)

        # Set memory per proc.
        #FIXME: Fixer may have changed the memory per proc and should not be resetted by ParalConf
        #self.set_mem_per_proc(pconf.mem_per_proc)
        return pconf

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        #app("[Task policy]\n%s" % str(self.policy))

        for i, qad in enumerate(self.qads):
            app("[Qadapter %d]\n%s" % (i, str(qad)))
        app("Qadapter selected: %d" % self._qadpos)

        if self.has_db:
            app("[MongoDB database]:")
            app(str(self.db_connector))

        return "\n".join(lines)

    @property
    def has_db(self):
        """True if we are using MongoDB database"""
        return bool(self.db_connector)

    @property
    def has_omp(self):
        """True if we are using OpenMP parallelization."""
        return self.qadapter.has_omp

    @property
    def num_cores(self):
        """Total number of CPUs used to run the task."""
        return self.qadapter.num_cores

    @property
    def mpi_procs(self):
        """Number of MPI processes."""
        return self.qadapter.mpi_procs

    @property
    def mem_per_proc(self):
        """Memory per MPI process."""
        return self.qadapter.mem_per_proc

    @property
    def omp_threads(self):
        """Number of OpenMP threads"""
        return self.qadapter.omp_threads

    def deepcopy(self):
        """Deep copy of self."""
        return copy.deepcopy(self)

    def set_mpi_procs(self, mpi_procs):
        """Set the number of MPI processes to use."""
        self.qadapter.set_mpi_procs(mpi_procs)

    def set_omp_threads(self, omp_threads):
        """Set the number of OpenMp threads to use."""
        self.qadapter.set_omp_threads(omp_threads)

    def set_mem_per_proc(self, mem_mb):
        """Set the memory (in Megabytes) per CPU."""
        self.qadapter.set_mem_per_proc(mem_mb)

    @property
    def max_cores(self):
        """
        Maximum number of cores that can be used.
        This value is mainly used in the autoparal part to get the list of possible configurations.
        """
        return max(q.hint_cores for q in self.qads)

    def get_njobs_in_queue(self, username=None):
        """
        returns the number of jobs in the queue,
        returns None when the number of jobs cannot be determined.

        Args:
            username: (str) the username of the jobs to count (default is to autodetect)
        """
        return self.qadapter.get_njobs_in_queue(username=username)

    def cancel(self, job_id):
        """Cancel the job. Returns exit status."""
        return self.qadapter.cancel(job_id)

    def write_jobfile(self, task, **kwargs):
        """
        Write the submission script. Return the path of the script

        ================  ============================================
        kwargs            Meaning
        ================  ============================================
        exec_args         List of arguments passed to task.executable.
                          Default: no arguments.

        ================  ============================================
        """
        script = self.qadapter.get_script_str(
            job_name=task.name,
            launch_dir=task.workdir,
            executable=task.executable,
            qout_path=task.qout_file.path,
            qerr_path=task.qerr_file.path,
            stdin=task.files_file.path,
            stdout=task.log_file.path,
            stderr=task.stderr_file.path,
            exec_args=kwargs.pop("exec_args", []),
        )

        # Write the script.
        with open(task.job_file.path, "w") as fh:
            fh.write(script)
            task.job_file.chmod(0o740)
            return task.job_file.path

    def launch(self, task, **kwargs):
        """
        Build the input files and submit the task via the :class:`Qadapter`

        Args:
            task: :class:`TaskObject`

        Returns:
            Process object.
        """
        if task.status == task.S_LOCKED:
            raise ValueError("You shall not submit a locked task!")

        # Build the task
        task.build()

        # Pass information on the time limit to Abinit (we always assume ndtset == 1)
        if isinstance(task, AbinitTask):
            args = kwargs.get("exec_args", [])
            if args is None: args = []
            args = args[:]
            args.append("--timelimit %s" % qu.time2slurm(self.qadapter.timelimit))
            kwargs["exec_args"] = args

        # Write the submission script
        script_file = self.write_jobfile(task, **kwargs)

        # Submit the task and save the queue id.
        try:
            qjob, process = self.qadapter.submit_to_queue(script_file)
            task.set_status(task.S_SUB, msg='Submitted to queue')
            task.set_qjob(qjob)
            return process

        except self.qadapter.MaxNumLaunchesError as exc:
            # TODO: Here we should try to switch to another qadapter
            # 1) Find a new parallel configuration in those stored in task.pconfs
            # 2) Change the input file.
            # 3) Regenerate the submission script
            # 4) Relaunch
            task.set_status(task.S_ERROR, msg="max_num_launches reached: %s" % str(exc))
            raise

    def get_collection(self, **kwargs):
        """Return the MongoDB collection used to store the results."""
        return self.db_connector.get_collection(**kwargs)

    def increase_mem(self):
        # OLD
        # with GW calculations in mind with GW mem = 10,
        # the response fuction is in memory and not distributed
        # we need to increase memory if jobs fail ...
        # return self.qadapter.more_mem_per_proc()
        try:
            self.qadapter.more_mem_per_proc()
        except QueueAdapterError:
            # here we should try to switch to an other qadapter
            raise ManagerIncreaseError('manager failed to increase mem')

    def increase_ncpus(self):
        """
        increase the number of cpus, first ask the current qadapter, if that one raises a QadapterIncreaseError
        switch to the next qadapter. If all fail raise an ManagerIncreaseError
        """
        try:
            self.qadapter.more_cores()
        except QueueAdapterError:
            # here we should try to switch to an other qadapter
            raise ManagerIncreaseError('manager failed to increase ncpu')

    def increase_resources(self):
        try:
            self.qadapter.more_cores()
            return
        except QueueAdapterError:
            pass

        try:
            self.qadapter.more_mem_per_proc()
        except QueueAdapterError:
            # here we should try to switch to an other qadapter
            raise ManagerIncreaseError('manager failed to increase resources')

    def exclude_nodes(self, nodes):
        try:
            self.qadapter.exclude_nodes(nodes=nodes)
        except QueueAdapterError:
            # here we should try to switch to an other qadapter
            raise ManagerIncreaseError('manager failed to exclude nodes')

    def increase_time(self):
        try:
            self.qadapter.more_time()
        except QueueAdapterError:
            # here we should try to switch to an other qadapter
            raise ManagerIncreaseError('manager failed to increase time')


class AbinitBuild:
    """
    This object stores information on the options used to build Abinit

        .. attribute:: info
            String with build information as produced by `abinit -b`

        .. attribute:: version
            Abinit version number e.g 8.0.1 (string)

        .. attribute:: has_netcdf
            True if netcdf is enabled.

        .. attribute:: has_omp
            True if OpenMP is enabled.

        .. attribute:: has_mpi
            True if MPI is enabled.

        .. attribute:: has_mpiio
            True if MPI-IO is supported.
    """
    def __init__(self, workdir=None, manager=None):
        manager = TaskManager.as_manager(manager).to_shell_manager(mpi_procs=1)

        # Build a simple manager to run the job in a shell subprocess
        import tempfile
        workdir = tempfile.mkdtemp() if workdir is None else workdir

        # Generate a shell script to execute `abinit -b`
        stdout = os.path.join(workdir, "run.abo")
        script = manager.qadapter.get_script_str(
            job_name="abinit_b",
            launch_dir=workdir,
            executable="abinit",
            qout_path=os.path.join(workdir, "queue.qout"),
            qerr_path=os.path.join(workdir, "queue.qerr"),
            #stdin=os.path.join(workdir, "run.files"),
            stdout=stdout,
            stderr=os.path.join(workdir, "run.err"),
            exec_args=["-b"],
        )

        # Execute the script.
        script_file = os.path.join(workdir, "job.sh")
        with open(script_file, "wt") as fh:
            fh.write(script)
        qjob, process = manager.qadapter.submit_to_queue(script_file)
        process.wait()

        if process.returncode != 0:
            logger.critical("Error while executing %s" % script_file)
            print("stderr:\n", process.stderr.read())
            #print("stdout:", process.stdout.read())

        # To avoid: ResourceWarning: unclosed file <_io.BufferedReader name=87> in py3k
        process.stderr.close()

        with open(stdout, "rt") as fh:
            self.info = fh.read()

        # info string has the following format.
        """
        === Build Information ===
         Version       : 8.0.1
         Build target  : x86_64_darwin15.0.0_gnu5.3
         Build date    : 20160122

        === Compiler Suite ===
         C compiler       : gnu
         C++ compiler     : gnuApple
         Fortran compiler : gnu5.3
         CFLAGS           :  -g -O2 -mtune=native -march=native
         CXXFLAGS         :  -g -O2 -mtune=native -march=native
         FCFLAGS          :  -g -ffree-line-length-none
         FC_LDFLAGS       :

        === Optimizations ===
         Debug level        : basic
         Optimization level : standard
         Architecture       : unknown_unknown

        === Multicore ===
         Parallel build : yes
         Parallel I/O   : yes
         openMP support : no
         GPU support    : no

        === Connectors / Fallbacks ===
         Connectors on : yes
         Fallbacks on  : yes
         DFT flavor    : libxc-fallback+atompaw-fallback+wannier90-fallback
         FFT flavor    : none
         LINALG flavor : netlib
         MATH flavor   : none
         TIMER flavor  : abinit
         TRIO flavor   : netcdf+etsf_io-fallback

        === Experimental features ===
         Bindings            : @enable_bindings@
         Exports             : no
         GW double-precision : yes

        === Bazaar branch information ===
         Branch ID : gmatteo@gmac-20160112110440-lf6exhneqim9082h
         Revision  : 1226
         Committed : 0
        """
        self.version = "0.0.0"
        self.has_netcdf = False
        self.has_omp = False
        self.has_mpi, self.has_mpiio = False, False

        def yesno2bool(line):
            ans = line.split()[-1].lower()
            return dict(yes=True, no=False, auto=True)[ans]

        # Parse info.
        for line in self.info.splitlines():
            if "Version" in line: self.version = line.split()[-1]
            if "TRIO flavor" in line:
                self.has_netcdf = "netcdf" in line
            if "openMP support" in line: self.has_omp = yesno2bool(line)
            if "Parallel build" in line: self.has_mpi = yesno2bool(line)
            if "Parallel I/O" in line: self.has_mpiio = yesno2bool(line)

    def __str__(self):
        lines = []
        app = lines.append
        app("Abinit Build Information:")
        app("    Abinit version: %s" % self.version)
        app("    MPI: %s, MPI-IO: %s, OpenMP: %s" % (self.has_mpi, self.has_mpiio, self.has_omp))
        app("    Netcdf: %s" % self.has_netcdf)
        return "\n".join(lines)

    def version_ge(self, version_string):
        """True is Abinit version is >= version_string"""
        return self.compare_version(version_string, ">=")

    def compare_version(self, version_string, op):
        """Compare Abinit version to `version_string` with operator `op`"""
        from pkg_resources import parse_version
        from monty.operator import operator_from_str
        op = operator_from_str(op)
        return op(parse_version(self.version), parse_version(version_string))


class FakeProcess:
    """
    This object is attached to a :class:`Task` instance if the task has not been submitted
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


class MyTimedelta(datetime.timedelta):
    """A customized version of timedelta whose __str__ method doesn't print microseconds."""
    def __new__(cls, days, seconds, microseconds):
        return datetime.timedelta.__new__(cls, days, seconds, microseconds)

    def __str__(self):
        """Remove microseconds from timedelta default __str__"""
        s = super(MyTimedelta, self).__str__()
        microsec = s.find(".")
        if microsec != -1: s = s[:microsec]
        return s

    @classmethod
    def as_timedelta(cls, delta):
        """Convert delta into a MyTimedelta object."""
        # Cannot monkey patch the __class__ and must pass through __new__ as the object is immutable.
        if isinstance(delta, cls): return delta
        return cls(delta.days, delta.seconds, delta.microseconds)


class TaskDateTimes:
    """
    Small object containing useful :class:`datetime.datatime` objects associated to important events.

    .. attributes:

        init: initialization datetime
        submission: submission datetime
        start: Begin of execution.
        end: End of execution.
    """
    def __init__(self):
        self.init = datetime.datetime.now()
        self.submission, self.start, self.end = None, None, None

    def __str__(self):
        lines = []
        app = lines.append

        app("Initialization done on: %s" % self.init)
        if self.submission is not None: app("Submitted on: %s" % self.submission)
        if self.start is not None: app("Started on: %s" % self.start)
        if self.end is not None: app("Completed on: %s" % self.end)

        return "\n".join(lines)

    def reset(self):
        """Reinitialize the counters."""
        self = self.__class__()

    def get_runtime(self):
        """:class:`timedelta` with the run-time, None if the Task is not running"""
        if self.start is None: return None

        if self.end is None:
            delta = datetime.datetime.now() - self.start
        else:
            delta = self.end - self.start

        return MyTimedelta.as_timedelta(delta)

    def get_time_inqueue(self):
        """
        :class:`timedelta` with the time spent in the Queue, None if the Task is not running

        .. note:

            This value is always greater than the real value computed by the resource manager
            as we start to count only when check_status sets the `Task` status to S_RUN.
        """
        if self.submission is None: return None

        if self.start is None:
            delta = datetime.datetime.now() - self.submission
        else:
            delta = self.start - self.submission
            # This happens when we read the exact start datetime from the ABINIT log file.
            if delta.total_seconds() < 0: delta = datetime.timedelta(seconds=0)

        return MyTimedelta.as_timedelta(delta)


class TaskError(NodeError):
    """Base Exception for :class:`Task` methods"""


class TaskRestartError(TaskError):
    """Exception raised while trying to restart the :class:`Task`."""


class Task(Node, metaclass=abc.ABCMeta):
    """
    A Task is a node that performs some kind of calculation.
    This is base class providing low-level methods.
    """
    # Use class attributes for TaskErrors so that we don't have to import them.
    Error = TaskError
    RestartError = TaskRestartError

    # List of `AbinitEvent` subclasses that are tested in the check_status method.
    # Subclasses should provide their own list if they need to check the converge status.
    CRITICAL_EVENTS = []

    # Prefixes for Abinit (input, output, temporary) files.
    Prefix = collections.namedtuple("Prefix", "idata odata tdata")
    pj = os.path.join

    prefix = Prefix(pj("indata", "in"), pj("outdata", "out"), pj("tmpdata", "tmp"))
    del Prefix, pj

    def __init__(self, input, workdir=None, manager=None, deps=None):
        """
        Args:
            input: :class:`AbinitInput` object.
            workdir: Path to the working directory.
            manager: :class:`TaskManager` object.
            deps: Dictionary specifying the dependency of this node.
                  None means that this Task has no dependency.
        """
        # Init the node
        super(Task, self).__init__()

        self._input = input

        if workdir is not None:
            self.set_workdir(workdir)

        if manager is not None:
            self.set_manager(manager)

        # Handle possible dependencies.
        if deps:
            self.add_deps(deps)

        # Date-time associated to submission, start and end.
        self.datetimes = TaskDateTimes()

        # Count the number of restarts.
        self.num_restarts = 0

        self._qjob = None
        self.queue_errors = []
        self.abi_errors = []

        # two flags that provide, dynamically, information on the scaling behavious of a task. If any process of fixing
        # finds none scaling behaviour, they should be switched. If a task type is clearly not scaling they should be
        # swiched.
        self.mem_scales = True
        self.load_scales = True

    def __getstate__(self):
        """
        Return state is pickled as the contents for the instance.

        In this case we just remove the process since Subprocess objects cannot be pickled.
        This is the reason why we have to store the returncode in self._returncode instead
        of using self.process.returncode.
        """
        return {k: v for k, v in self.__dict__.items() if k not in ["_process"]}

    #@check_spectator
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
        # This file is produced by Abinit if nprocs > 1 and MPI_ABORT.
        self.mpiabort_file = File(os.path.join(self.workdir, "__ABI_MPIABORTFILE__"))

        # Directories with input|output|temporary data.
        self.wdir = Directory(self.workdir)
        self.indir = Directory(os.path.join(self.workdir, "indata"))
        self.outdir = Directory(os.path.join(self.workdir, "outdata"))
        self.tmpdir = Directory(os.path.join(self.workdir, "tmpdata"))

        # stderr and output file of the queue manager. Note extensions.
        self.qerr_file = File(os.path.join(self.workdir, "queue.qerr"))
        self.qout_file = File(os.path.join(self.workdir, "queue.qout"))

    def set_manager(self, manager):
        """Set the :class:`TaskManager` used to launch the Task."""
        self.manager = manager.deepcopy()

    @property
    def work(self):
        """The :class:`Work` containing this `Task`."""
        return self._work

    def set_work(self, work):
        """Set the :class:`Work` associated to this `Task`."""
        if not hasattr(self, "_work"):
            self._work = work
        else:
            if self._work != work:
                raise ValueError("self._work != work")

    @property
    def flow(self):
        """The :class:`Flow` containing this `Task`."""
        return self.work.flow

    @lazy_property
    def pos(self):
        """The position of the task in the :class:`Flow`"""
        for i, task in enumerate(self.work):
            if self == task:
                return self.work.pos, i
        raise ValueError("Cannot find the position of %s in flow %s" % (self, self.flow))

    @property
    def pos_str(self):
        """String representation of self.pos"""
        return "w" + str(self.pos[0]) + "_t" + str(self.pos[1])

    @property
    def num_launches(self):
        """
        Number of launches performed. This number includes both possible ABINIT restarts
        as well as possible launches done due to errors encountered with the resource manager
        or the hardware/software."""
        return sum(q.num_launches for q in self.manager.qads)

    @property
    def input(self):
        """AbinitInput object."""
        return self._input

    def get_inpvar(self, varname, default=None):
        """Return the value of the ABINIT variable varname, None if not present."""
        return self.input.get(varname, default)

    def set_vars(self, *args, **kwargs):
        """
        Set the values of the ABINIT variables in the input file. Return dict with old values.
        """
        kwargs.update(dict(*args))
        old_values = {vname: self.input.get(vname) for vname in kwargs}
        self.input.set_vars(**kwargs)
        if kwargs or old_values:
            self.history.info("Setting input variables: %s" % str(kwargs))
            self.history.info("Old values: %s" % str(old_values))

        return old_values

    @property
    def initial_structure(self):
        """Initial structure of the task."""
        return self.input.structure

    def make_input(self, with_header=False):
        """Construct the input file of the calculation."""
        s = str(self.input)
        if with_header: s = str(self) + "\n" + s
        return s

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

    @property
    @abc.abstractmethod
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

    @property
    def is_abinit_task(self):
        """True if this task is a subclass of AbinitTask."""
        return isinstance(self, AbinitTask)

    @property
    def is_anaddb_task(self):
        """True if this task is a subclass of OpticTask."""
        return isinstance(self, AnaddbTask)

    @property
    def is_optic_task(self):
        """True if this task is a subclass of OpticTask."""
        return isinstance(self, OpticTask)

    @property
    def is_completed(self):
        """True if the task has been executed."""
        return self.status >= self.S_DONE

    @property
    def can_run(self):
        """The task can run if its status is < S_SUB and all the other dependencies (if any) are done!"""
        all_ok = all(stat == self.S_OK for stat in self.deps_status)
        return self.status < self.S_SUB and self.status != self.S_LOCKED and all_ok

    #@check_spectator
    def cancel(self):
        """Cancel the job. Returns 1 if job was cancelled."""
        if self.queue_id is None: return 0
        if self.status >= self.S_DONE: return 0

        exit_status = self.manager.cancel(self.queue_id)
        if exit_status != 0:
            logger.warning("manager.cancel returned exit_status: %s" % exit_status)
            return 0

        # Remove output files and reset the status.
        self.history.info("Job %s cancelled by user" % self.queue_id)
        self.reset()
        return 1

    def with_fixed_mpi_omp(self, mpi_procs, omp_threads):
        """
        Disable autoparal and force execution with `mpi_procs` MPI processes
        and `omp_threads` OpenMP threads. Useful for generating benchmarks.
        """
        manager = self.manager if hasattr(self, "manager") else self.flow.manager
        self.manager = manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)

    #def set_max_ncores(self, max_ncores):
    #    """
    #    """
    #    manager = self.manager if hasattr(self, "manager") else self.flow.manager
    #    self.manager = manager.new_with_max_ncores(mpi_procs, omp_threads)

    #@check_spectator
    def _on_done(self):
        self.fix_ofiles()

    #@check_spectator
    def _on_ok(self):
        # Fix output file names.
        self.fix_ofiles()

        # Get results
        results = self.on_ok()

        self.finalized = True

        return results

    #@check_spectator
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
        return dict(returncode=0, message="Calling on_all_ok of the base class!")

    #@check_spectator
    def fix_ofiles(self):
        """
        This method is called when the task reaches S_OK.
        It changes the extension of particular output files
        produced by Abinit so that the 'official' extension
        is preserved e.g. out_1WF14 --> out_1WF
        """
        filepaths = self.outdir.list_filepaths()
        #logger.info("in fix_ofiles with filepaths %s" % list(filepaths))

        old2new = FilepathFixer().fix_paths(filepaths)

        for old, new in old2new.items():
            self.history.info("will rename old %s to new %s" % (old, new))
            os.rename(old, new)

    #@check_spectator
    def _restart(self, submit=True):
        """
        Called by restart once we have finished preparing the task for restarting.

        Return:
            True if task has been restarted
        """
        self.set_status(self.S_READY, msg="Restarted on %s" % time.asctime())

        # Increase the counter.
        self.num_restarts += 1
        self.history.info("Restarted, num_restarts %d" % self.num_restarts)

        # Reset datetimes
        self.datetimes.reset()

        # Remove the lock file
        self.start_lockfile.remove()

        if submit:
            # Relaunch the task.
            fired = self.start()
            if not fired: self.history.warning("Restart failed")
        else:
            fired = False

        return fired

    #@check_spectator
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
            self.set_status(self.S_DONE, "status set to Done")

        return self._returncode

    def wait(self):
        """Wait for child process to terminate. Set and return returncode attribute."""
        self._returncode = self.process.wait()
        try:
            self.process.stderr.close()
        except:
            pass
        self.set_status(self.S_DONE, "status set to Done")

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
        self.set_status(self.S_DONE, "status set to Done")

        return stdoutdata, stderrdata

    def kill(self):
        """Kill the child."""
        self.process.kill()
        self.set_status(self.S_ERROR, "status set to Error by task.kill")
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
        # One should be able to reset 'Submitted' tasks (sometimes, they are not in the queue
        # and we want to restart them)
        #if self.status != self.S_SUB and self.status < self.S_DONE: return 1

        # Remove output files otherwise the EventParser will think the job is still running
        self.output_file.remove()
        self.log_file.remove()
        self.stderr_file.remove()
        self.start_lockfile.remove()
        self.qerr_file.remove()
        self.qout_file.remove()
        if self.mpiabort_file.exists:
            self.mpiabort_file.remove()

        self.set_status(self.S_INIT, msg="Reset on %s" % time.asctime())
        self.num_restarts = 0
        self.set_qjob(None)

        # Reset finalized flags.
        self.work.finalized = False
        self.flow.finalized = False

        return 0

    @property
    @return_none_if_raise(AttributeError)
    def queue_id(self):
        """Queue identifier returned by the Queue manager. None if not set"""
        return self.qjob.qid

    @property
    @return_none_if_raise(AttributeError)
    def qname(self):
        """Queue name identifier returned by the Queue manager. None if not set"""
        return self.qjob.qname

    @property
    def qjob(self):
        return self._qjob

    def set_qjob(self, qjob):
        """Set info on queue after submission."""
        self._qjob = qjob

    @property
    def has_queue(self):
        """True if we are submitting jobs via a queue manager."""
        return self.manager.qadapter.QTYPE.lower() != "shell"

    @property
    def num_cores(self):
        """Total number of CPUs used to run the task."""
        return self.manager.num_cores

    @property
    def mpi_procs(self):
        """Number of CPUs used for MPI."""
        return self.manager.mpi_procs

    @property
    def omp_threads(self):
        """Number of CPUs used for OpenMP."""
        return self.manager.omp_threads

    @property
    def mem_per_proc(self):
        """Memory per MPI process."""
        return Memory(self.manager.mem_per_proc, "Mb")

    @property
    def status(self):
        """Gives the status of the task."""
        return self._status

    def lock(self, source_node):
        """Lock the task, source is the :class:`Node` that applies the lock."""
        if self.status != self.S_INIT:
            raise ValueError("Trying to lock a task with status %s" % self.status)

        self._status = self.S_LOCKED
        self.history.info("Locked by node %s", source_node)

    def unlock(self, source_node, check_status=True):
        """
        Unlock the task, set its status to `S_READY` so that the scheduler can submit it.
        source_node is the :class:`Node` that removed the lock
        Call task.check_status if check_status is True.
        """
        if self.status != self.S_LOCKED:
            raise RuntimeError("Trying to unlock a task with status %s" % self.status)

        self._status = self.S_READY
        if check_status: self.check_status()
        self.history.info("Unlocked by %s", source_node)

    #@check_spectator
    def set_status(self, status, msg):
        """
        Set and return the status of the task.

        Args:
            status: Status object or string representation of the status
            msg: string with human-readable message used in the case of errors.
        """
        # truncate string if it's long. msg will be logged in the object and we don't want to waste memory.
        if len(msg) > 2000:
            msg = msg[:2000]
            msg += "\n... snip ...\n"

        # Locked files must be explicitly unlocked
        if self.status == self.S_LOCKED or status == self.S_LOCKED:
            err_msg = (
                 "Locked files must be explicitly unlocked before calling set_status but\n"
                 "task.status = %s, input status = %s" % (self.status, status))
            raise RuntimeError(err_msg)

        status = Status.as_status(status)

        changed = True
        if hasattr(self, "_status"):
            changed = (status != self._status)

        self._status = status

        if status == self.S_RUN:
            # Set datetimes.start when the task enters S_RUN
            if self.datetimes.start is None:
                self.datetimes.start = datetime.datetime.now()

        # Add new entry to history only if the status has changed.
        if changed:
            if status == self.S_SUB:
                self.datetimes.submission = datetime.datetime.now()
                self.history.info("Submitted with MPI=%s, Omp=%s, Memproc=%.1f [Gb] %s " % (
                    self.mpi_procs, self.omp_threads, self.mem_per_proc.to("Gb"), msg))

            elif status == self.S_OK:
                self.history.info("Task completed %s", msg)

            elif status == self.S_ABICRITICAL:
                self.history.info("Status set to S_ABI_CRITICAL due to: %s", msg)

            else:
                self.history.info("Status changed to %s. msg: %s", status, msg)

        #######################################################
        # The section belows contains callbacks that should not
        # be executed if we are in spectator_mode
        #######################################################
        if status == self.S_DONE:
            # Execute the callback
            self._on_done()

        if status == self.S_OK:
            # Finalize the task.
            if not self.finalized:
                self._on_ok()

                # here we remove the output files of the task and of its parents.
                if self.gc is not None and self.gc.policy == "task":
                    self.clean_output_files()

            if self.status == self.S_OK:
                # Because _on_ok might have changed the status.
                self.send_signal(self.S_OK)

        return status

    def check_status(self):
        """
        This function checks the status of the task by inspecting the output and the
        error files produced by the application and by the queue manager.
        """
        # 1) see it the job is blocked
        # 2) see if an error occured at submitting the job the job was submitted, TODO these problems can be solved
        # 3) see if there is output
        # 4) see if abinit reports problems
        # 5) see if both err files exist and are empty
        # 6) no output and no err files, the job must still be running
        # 7) try to find out what caused the problems
        # 8) there is a problem but we did not figure out what ...
        # 9) the only way of landing here is if there is a output file but no err files...

        # 1) A locked task can only be unlocked by calling set_status explicitly.
        # an errored task, should not end up here but just to be sure
        black_list = (self.S_LOCKED, self.S_ERROR)
        #if self.status in black_list: return self.status

        # 2) Check the returncode of the job script
        if self.returncode != 0:
            msg = "job.sh return code: %s\nPerhaps the job was not submitted properly?" % self.returncode
            return self.set_status(self.S_QCRITICAL, msg=msg)

        # If we have an abort file produced by Abinit
        if self.mpiabort_file.exists:
            return self.set_status(self.S_ABICRITICAL, msg="Found ABINIT abort file")

        # Analyze the stderr file for Fortran runtime errors.
        # getsize is 0 if the file is empty or it does not exist.
        err_msg = None
        if self.stderr_file.getsize() != 0:
            err_msg = self.stderr_file.read()

        # Analyze the stderr file of the resource manager runtime errors.
        # TODO: Why are we looking for errors in queue.qerr?
        qerr_info = None
        if self.qerr_file.getsize() != 0:
            qerr_info = self.qerr_file.read()

        # Analyze the stdout file of the resource manager (needed for PBS !)
        qout_info = None
        if self.qout_file.getsize():
            qout_info = self.qout_file.read()

        # Start to check ABINIT status if the output file has been created.
        #if self.output_file.getsize() != 0:
        if self.output_file.exists:
            try:
                report = self.get_event_report()
            except Exception as exc:
                msg = "%s exception while parsing event_report:\n%s" % (self, exc)
                return self.set_status(self.S_ABICRITICAL, msg=msg)

            if report is None:
                return self.set_status(self.S_ERROR, msg="got None report!")

            if report.run_completed:
                # Here we  set the correct timing data reported by Abinit
                self.datetimes.start = report.start_datetime
                self.datetimes.end = report.end_datetime

                # Check if the calculation converged.
                not_ok = report.filter_types(self.CRITICAL_EVENTS)
                if not_ok:
                    return self.set_status(self.S_UNCONVERGED, msg='status set to unconverged based on abiout')
                else:
                    return self.set_status(self.S_OK, msg="status set to ok based on abiout")

            # Calculation still running or errors?
            if report.errors:
                # Abinit reported problems
                logger.debug('Found errors in report')
                for error in report.errors:
                    logger.debug(str(error))
                    try:
                        self.abi_errors.append(error)
                    except AttributeError:
                        self.abi_errors = [error]

                # The job is unfixable due to ABINIT errors
                logger.debug("%s: Found Errors or Bugs in ABINIT main output!" % self)
                msg = "\n".join(map(repr, report.errors))
                return self.set_status(self.S_ABICRITICAL, msg=msg)

            # 5)
            if self.stderr_file.exists and not err_msg:
                if self.qerr_file.exists and not qerr_info:
                    # there is output and no errors
                    # The job still seems to be running
                    return self.set_status(self.S_RUN, msg='there is output and no errors: job still seems to be running')

        # 6)
        if not self.output_file.exists:
            logger.debug("output_file does not exists")
            if not self.stderr_file.exists and not self.qerr_file.exists:
                # No output at allThe job is still in the queue.
                return self.status

        # 7) Analyze the files of the resource manager and abinit and execution err (mvs)
        # MG: This section has been disabled: several portability issues
        # Need more robust logic in error_parser, perhaps logic provided by users via callbacks.
        if False and (qerr_info or qout_info):
            from pymatgen.io.abinit.scheduler_error_parsers import get_parser
            scheduler_parser = get_parser(self.manager.qadapter.QTYPE, err_file=self.qerr_file.path,
                                          out_file=self.qout_file.path, run_err_file=self.stderr_file.path)

            if scheduler_parser is None:
                return self.set_status(self.S_QCRITICAL,
                                       msg="Cannot find scheduler_parser for qtype %s" % self.manager.qadapter.QTYPE)

            scheduler_parser.parse()

            if scheduler_parser.errors:
                # Store the queue errors in the task
                self.queue_errors = scheduler_parser.errors
                # The job is killed or crashed and we know what happened
                msg = "scheduler errors found:\n%s" % str(scheduler_parser.errors)
                return self.set_status(self.S_QCRITICAL, msg=msg)

            elif lennone(qerr_info) > 0:
                # if only qout_info, we are not necessarily in QCRITICAL state,
                # since there will always be info in the qout file
                self.history.info('Found unknown message in the queue qerr file: %s' % str(qerr_info))
                #try:
                #    rt = self.datetimes.get_runtime().seconds
                #except:
                #    rt = -1.0
                #tl = self.manager.qadapter.timelimit
                #if rt > tl:
                #    msg += 'set to error : runtime (%s) exceded walltime (%s)' % (rt, tl)
                #    print(msg)
                #    return self.set_status(self.S_ERROR, msg=msg)
                # The job may be killed or crashed but we don't know what happened
                # It may also be that an innocent message was written to qerr, so we wait for a while
                # it is set to QCritical, we will attempt to fix it by running on more resources

        # 8) analyzing the err files and abinit output did not identify a problem
        # but if the files are not empty we do have a problem but no way of solving it:
        # The job is killed or crashed but we don't know what happend
        # it is set to QCritical, we will attempt to fix it by running on more resources
        if err_msg:
            msg = 'Found error message:\n %s' % str(err_msg)
            self.history.warning(msg)
            #return self.set_status(self.S_QCRITICAL, msg=msg)

        # 9) if we still haven't returned there is no indication of any error and the job can only still be running
        # but we should actually never land here, or we have delays in the file system ....
        # print('the job still seems to be running maybe it is hanging without producing output... ')

        # Check time of last modification.
        if self.output_file.exists and \
           (time.time() - self.output_file.get_stat().st_mtime > self.manager.policy.frozen_timeout):
            msg = "Task seems to be frozen, last change more than %s [s] ago" % self.manager.policy.frozen_timeout
            return self.set_status(self.S_ERROR, msg=msg)

        # Handle weird case in which either run.abo, or run.log have not been produced
        #if self.status not in (self.S_INIT, self.S_READY) and (not self.output.file.exists or not self.log_file.exits):
        #    msg = "Task have been submitted but cannot find the log file or the output file"
        #    return self.set_status(self.S_ERROR, msg)

        return self.set_status(self.S_RUN, msg='final option: nothing seems to be wrong, the job must still be running')

    def reduce_memory_demand(self):
        """
        Method that can be called by the scheduler to decrease the memory demand of a specific task.
        Returns True in case of success, False in case of Failure.
        Should be overwritten by specific tasks.
        """
        return False

    def speed_up(self):
        """
        Method that can be called by the flow to decrease the time needed for a specific task.
        Returns True in case of success, False in case of Failure
        Should be overwritten by specific tasks.
        """
        return False

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
        self.history.info("Linking path %s --> %s" % (filepath, infile))

        if not os.path.exists(infile):
            os.symlink(filepath, infile)
        else:
            if os.path.realpath(infile) != filepath:
                raise self.Error("infile %s does not point to filepath %s" % (infile, filepath))

    def make_links(self):
        """
        Create symbolic links to the output files produced by the other tasks.

        .. warning::

            This method should be called only when the calculation is READY because
            it uses a heuristic approach to find the file to link.
        """
        for dep in self.deps:
            filepaths, exts = dep.get_filepaths_and_exts()

            for path, ext in zip(filepaths, exts):
                logger.info("Need path %s with ext %s" % (path, ext))
                dest = self.ipath_from_ext(ext)

                if not os.path.exists(path):
                    # Try netcdf file.
                    # TODO: this case should be treated in a cleaner way.
                    path += ".nc"
                    if os.path.exists(path): dest += ".nc"

                if not os.path.exists(path):
                    raise self.Error("%s: %s is needed by this task but it does not exist" % (self, path))

                if path.endswith(".nc") and not dest.endswith(".nc"): # NC --> NC file
                    dest += ".nc"

                # Link path to dest if dest link does not exist.
                # else check that it points to the expected file.
                logger.debug("Linking path %s --> %s" % (path, dest))
                if not os.path.exists(dest):
                    os.symlink(path, dest)
                else:
                    # check links but only if we haven't performed the restart.
                    # in this case, indeed we may have replaced the file pointer with the
                    # previous output file of the present task.
                    if os.path.realpath(dest) != path and self.num_restarts == 0:
                        raise self.Error("\nDestination:\n %s\ndoes not point to path:\n %s" % (dest, path))

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

    def get_event_report(self, source="log"):
        """
        Analyzes the main logfile of the calculation for possible Errors or Warnings.
        If the ABINIT abort file is found, the error found in this file are added to
        the output report.

        Args:
            source: "output" for the main output file,"log" for the log file.

        Returns:
            :class:`EventReport` instance or None if the source file file does not exist.
        """
        # By default, we inspect the main log file.
        ofile = {
            "output": self.output_file,
            "log": self.log_file}[source]

        parser = events.EventsParser()

        if not ofile.exists:
            if not self.mpiabort_file.exists:
                return None
            else:
                # ABINIT abort file without log!
                abort_report = parser.parse(self.mpiabort_file.path)
                return abort_report

        try:
            report = parser.parse(ofile.path)
            #self._prev_reports[source] = report

            # Add events found in the ABI_MPIABORTFILE.
            if self.mpiabort_file.exists:
                logger.critical("Found ABI_MPIABORTFILE!!!!!")
                abort_report = parser.parse(self.mpiabort_file.path)
                if len(abort_report) != 1:
                    logger.critical("Found more than one event in ABI_MPIABORTFILE")

                # Weird case: empty abort file, let's skip the part
                # below and hope that the log file contains the error message.
                #if not len(abort_report): return report

                # Add it to the initial report only if it differs
                # from the last one found in the main log file.
                last_abort_event = abort_report[-1]
                if report and last_abort_event != report[-1]:
                    report.append(last_abort_event)
                else:
                    report.append(last_abort_event)

            return report

        #except parser.Error as exc:
        except Exception as exc:
            # Return a report with an error entry with info on the exception.
            msg = "%s: Exception while parsing ABINIT events:\n %s" % (ofile, str(exc))
            self.set_status(self.S_ABICRITICAL, msg=msg)
            return parser.report_exception(ofile.path, exc)

    def get_results(self, **kwargs):
        """
        Returns :class:`NodeResults` instance.
        Subclasses should extend this method (if needed) by adding
        specialized code that performs some kind of post-processing.
        """
        # Check whether the process completed.
        if self.returncode is None:
            raise self.Error("return code is None, you should call wait, communicate or poll")

        if self.status is None or self.status < self.S_DONE:
            raise self.Error("Task is not completed")

        return self.Results.from_node(self)

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

    #@check_spectator
    def build(self, *args, **kwargs):
        """
        Creates the working directory and the input files of the :class:`Task`.
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

    #@check_spectator
    def rmtree(self, exclude_wildcard=""):
        """
        Remove all files and directories in the working directory

        Args:
            exclude_wildcard: Optional string with regular expressions separated by |.
                Files matching one of the regular expressions will be preserved.
                example: exclude_wildcard="*.nc|*.txt" preserves all the files whose extension is in ["nc", "txt"].
        """
        if not exclude_wildcard:
            shutil.rmtree(self.workdir)

        else:
            w = WildCard(exclude_wildcard)

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

    def clean_output_files(self, follow_parents=True):
        """
        This method is called when the task reaches S_OK. It removes all the output files
        produced by the task that are not needed by its children as well as the output files
        produced by its parents if no other node needs them.

        Args:
            follow_parents: If true, the output files of the parents nodes will be removed if possible.

        Return:
            list with the absolute paths of the files that have been removed.
        """
        paths = []
        if self.status != self.S_OK:
            logger.warning("Calling task.clean_output_files on a task whose status != S_OK")

        # Remove all files in tmpdir.
        self.tmpdir.clean()

        # Find the file extensions that should be preserved since these files are still
        # needed by the children who haven't reached S_OK
        except_exts = set()
        for child in self.get_children():
            if child.status == self.S_OK: continue
            # Find the position of self in child.deps and add the extensions.
            i = [dep.node for dep in child.deps].index(self)
            except_exts.update(child.deps[i].exts)

        # Remove the files in the outdir of the task but keep except_exts.
        exts = self.gc.exts.difference(except_exts)
        #print("Will remove its extensions: ", exts)
        paths += self.outdir.remove_exts(exts)
        if not follow_parents: return paths

        # Remove the files in the outdir of my parents if all the possible dependencies have been fulfilled.
        for parent in self.get_parents():

            # Here we build a dictionary file extension --> list of child nodes requiring this file from parent
            # e.g {"WFK": [node1, node2]}
            ext2nodes = collections.defaultdict(list)
            for child in parent.get_children():
                if child.status == child.S_OK: continue
                i = [d.node for d in child.deps].index(parent)
                for ext in child.deps[i].exts:
                    ext2nodes[ext].append(child)

            # Remove extension only if no node depends on it!
            except_exts = [k for k, lst in ext2nodes.items() if lst]
            exts = self.gc.exts.difference(except_exts)
            #print("%s removes extensions %s from parent node %s" % (self, exts, parent))
            paths += parent.outdir.remove_exts(exts)

        self.history.info("Removed files: %s" % paths)
        return paths

    def setup(self):
        """Base class does not provide any hook."""

    #@check_spectator
    def start(self, **kwargs):
        """
        Starts the calculation by performing the following steps:

            - build dirs and files
            - call the _setup method
            - execute the job file by executing/submitting the job script.

        Main entry point for the `Launcher`.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        autoparal       False to skip the autoparal step (default True)
        exec_args       List of arguments passed to executable.
        ==============  ==============================================================

        Returns:
            1 if task was started, 0 otherwise.

        """
        if self.status >= self.S_SUB:
            raise self.Error("Task status: %s" % str(self.status))

        if self.start_lockfile.exists:
            self.history.warning("Found lock file: %s" % self.start_lockfile.path)
            return 0

        self.start_lockfile.write("Started on %s" % time.asctime())

        self.build()
        self._setup()

        # Add the variables needed to connect the node.
        for d in self.deps:
            cvars = d.connecting_vars()
            self.history.info("Adding connecting vars %s" % cvars)
            self.set_vars(cvars)

            # Get (python) data from other nodes
            d.apply_getters(self)

        # Automatic parallelization
        if kwargs.pop("autoparal", True) and hasattr(self, "autoparal_run"):
            try:
                self.autoparal_run()
            #except QueueAdapterError as exc:
            #    # If autoparal cannot find a qadapter to run the calculation raises an Exception
            #    self.history.critical(exc)
            #    msg = "Error while trying to run autoparal in task:%s\n%s" % (repr(task), straceback())
            #    cprint(msg, "yellow")
            #    self.set_status(self.S_QCRITICAL, msg=msg)
            #    return 0
            except Exception as exc:
                # Sometimes autoparal_run fails because Abinit aborts
                # at the level of the parser e.g. cannot find the spacegroup
                # due to some numerical noise in the structure.
                # In this case we call fix_abicritical and then we try to run autoparal again.
                self.history.critical("First call to autoparal failed with `%s`. Will try fix_abicritical" % exc)
                msg = "autoparal_fake_run raised:\n%s" % straceback()
                logger.critical(msg)

                fixed = self.fix_abicritical()
                if not fixed:
                    self.set_status(self.S_ABICRITICAL, msg="fix_abicritical could not solve the problem")
                    return 0

                try:
                    self.autoparal_run()
                    self.history.info("Second call to autoparal succeeded!")
                    #cprint("Second call to autoparal succeeded!", "green")

                except Exception as exc:
                    self.history.critical("Second call to autoparal failed with %s. Cannot recover!", exc)
                    msg = "Tried autoparal again but got:\n%s" % straceback()
                    cprint(msg, "red")
                    self.set_status(self.S_ABICRITICAL, msg=msg)
                    return 0

        # Start the calculation in a subprocess and return.
        self._process = self.manager.launch(self, **kwargs)
        return 1

    def start_and_wait(self, *args, **kwargs):
        """
        Helper method to start the task and wait for completion.

        Mainly used when we are submitting the task via the shell without passing through a queue manager.
        """
        self.start(*args, **kwargs)
        retcode = self.wait()
        return retcode

    def get_graphviz(self, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate task graph in the DOT language (only parents and children of this task).

        Args:
            engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']
            graph_attr: Mapping of (attribute, value) pairs for the graph.
            node_attr: Mapping of (attribute, value) pairs set for all nodes.
            edge_attr: Mapping of (attribute, value) pairs set for all edges.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph
        fg = Digraph("task", # filename="task_%s.gv" % os.path.basename(self.workdir),
            engine="dot" if engine == "automatic" else engine)

        # Set graph attributes.
        #fg.attr(label="%s@%s" % (self.__class__.__name__, self.relworkdir))
        fg.attr(label=repr(self))
        #fg.attr(fontcolor="white", bgcolor='purple:pink')
        #fg.attr(rankdir="LR", pagedir="BL")
        #fg.attr(constraint="false", pack="true", packMode="clust")
        fg.node_attr.update(color='lightblue2', style='filled')

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

        # Build cluster with tasks.
        cluster_name = "cluster%s" % self.work.name
        with fg.subgraph(name=cluster_name) as wg:
            wg.attr(**cluster_kwargs)
            wg.attr(label="%s (%s)" % (self.__class__.__name__, self.name))
            wg.node(self.name, **node_kwargs(self))

            # Connect task to children.
            for child in self.get_children():
                # Test if child is in the same work.
                myg = wg if child in self.work else fg
                myg.node(child.name, **node_kwargs(child))
                # Find file extensions required by this task
                i = [dep.node for dep in child.deps].index(self)
                edge_label = "+".join(child.deps[i].exts)
                myg.edge(self.name, child.name, label=edge_label, color=self.color_hex,
                         **edge_kwargs)

            # Connect task to parents
            for parent in self.get_parents():
                # Test if parent is in the same work.
                myg = wg if parent in self.work else fg
                myg.node(parent.name, **node_kwargs(parent))
                # Find file extensions required by self (task)
                i = [dep.node for dep in self.deps].index(parent)
                edge_label = "+".join(self.deps[i].exts)
                myg.edge(parent.name, self.name, label=edge_label, color=parent.color_hex,
                         **edge_kwargs)

        # Treat the case in which we have a work producing output for other tasks.
        #for work in self:
        #    children = work.get_children()
        #    if not children: continue
        #    cluster_name = "cluster%s" % work.name
        #    seen = set()
        #    for child in children:
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


class DecreaseDemandsError(Exception):
    """
    exception to be raised by a task if the request to decrease some demand, load or memory, could not be performed
    """


class AbinitTask(Task):
    """
    Base class defining an ABINIT calculation
    """
    Results = TaskResults

    @classmethod
    def from_input(cls, input, workdir=None, manager=None):
        """
        Create an instance of `AbinitTask` from an ABINIT input.

        Args:
            ainput: `AbinitInput` object.
            workdir: Path to the working directory.
            manager: :class:`TaskManager` object.
        """
        return cls(input, workdir=workdir, manager=manager)

    @classmethod
    def temp_shell_task(cls, inp, mpi_procs=1, workdir=None, manager=None):
        """
        Build a Task with a temporary workdir. The task is executed via the shell with 1 MPI proc.
        Mainly used for invoking Abinit to get important parameters needed to prepare the real task.

        Args:
            mpi_procs: Number of MPI processes to use.
        """
        # Build a simple manager to run the job in a shell subprocess
        import tempfile
        workdir = tempfile.mkdtemp() if workdir is None else workdir
        if manager is None: manager = TaskManager.from_user_config()

        # Construct the task and run it
        task = cls.from_input(inp, workdir=workdir, manager=manager.to_shell_manager(mpi_procs=mpi_procs))
        task.set_name('temp_shell_task')
        return task

    def setup(self):
        """
        Abinit has the very *bad* habit of changing the file extension by appending the characters in [A,B ..., Z]
        to the output file, and this breaks a lot of code that relies of the use of a unique file extension.
        Here we fix this issue by renaming run.abo to run.abo_[number] if the output file "run.abo" already
        exists. A few lines of code in python, a lot of problems if you try to implement this trick in Fortran90.
        """
        def rename_file(afile):
            """Helper function to rename :class:`File` objects. Return string for logging purpose."""
            # Find the index of the last file (if any).
            # TODO: Maybe it's better to use run.abo --> run(1).abo
            fnames = [f for f in os.listdir(self.workdir) if f.startswith(afile.basename)]
            nums = [int(f) for f in [f.split("_")[-1] for f in fnames] if f.isdigit()]
            last = max(nums) if nums else 0
            new_path = afile.path + "_" + str(last+1)

            os.rename(afile.path, new_path)
            return "Will rename %s to %s" % (afile.path, new_path)

        logs = []
        if self.output_file.exists: logs.append(rename_file(self.output_file))
        if self.log_file.exists: logs.append(rename_file(self.log_file))

        if logs:
            self.history.info("\n".join(logs))

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
        return self.input.pseudos

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return self.input.isnc

    @property
    def ispaw(self):
        """True if PAW calculation"""
        return self.input.ispaw

    @property
    def is_gs_task(self):
        """True if task is GsTask subclass."""
        return isinstance(self, GsTask)

    @property
    def is_dfpt_task(self):
        """True if task is a DftpTask subclass."""
        return isinstance(self, DfptTask)

    @lazy_property
    def cycle_class(self):
        """
        Return the subclass of ScfCycle associated to the task or
        None if no SCF algorithm if associated to the task.
        """
        if isinstance(self, RelaxTask):
            return abiinspect.Relaxation
        elif isinstance(self, GsTask):
            return abiinspect.GroundStateScfCycle
        elif self.is_dfpt_task:
            return abiinspect.D2DEScfCycle

        return None

    @property
    def filesfile_string(self):
        """String with the list of files and prefixes needed to execute ABINIT."""
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
        # Here we reorder the pseudos if the order is wrong.
        ord_pseudos = []

        znucl = [specie.number for specie in
                 self.input.structure.types_of_specie]

        for z in znucl:
            for p in self.pseudos:
                if p.Z == z:
                    ord_pseudos.append(p)
                    break
            else:
                raise ValueError("Cannot find pseudo with znucl %s in pseudos:\n%s" % (z, self.pseudos))

        for pseudo in ord_pseudos:
            app(pseudo.path)

        return "\n".join(lines)

    def set_pconfs(self, pconfs):
        """Set the list of autoparal configurations."""
        self._pconfs = pconfs

    @property
    def pconfs(self):
        """List of autoparal configurations."""
        try:
            return self._pconfs
        except AttributeError:
            return None

    def uses_paral_kgb(self, value=1):
        """True if the task is a GS Task and uses paral_kgb with the given value."""
        paral_kgb = self.get_inpvar("paral_kgb", 0)
        # paral_kgb is used only in the GS part.
        return paral_kgb == value and isinstance(self, GsTask)

    def _change_structure(self, new_structure):
        """Change the input structure."""
        # Compare new and old structure for logging purpose.
        # TODO: Write method of structure to compare self and other and return a dictionary
        old_structure = self.input.structure
        old_lattice = old_structure.lattice

        abc_diff = np.array(new_structure.lattice.abc) - np.array(old_lattice.abc)
        angles_diff = np.array(new_structure.lattice.angles) - np.array(old_lattice.angles)
        cart_diff = new_structure.cart_coords - old_structure.cart_coords
        displs = np.array([np.sqrt(np.dot(v, v)) for v in cart_diff])

        recs, tol_angle, tol_length = [], 10**-2, 10**-5

        if np.any(np.abs(angles_diff) > tol_angle):
            recs.append("new_agles - old_angles = %s" % angles_diff)

        if np.any(np.abs(abc_diff) > tol_length):
            recs.append("new_abc - old_abc = %s" % abc_diff)

        if np.any(np.abs(displs) > tol_length):
            min_pos, max_pos = displs.argmin(), displs.argmax()
            recs.append("Mean displ: %.2E, Max_displ: %.2E (site %d), min_displ: %.2E (site %d)" %
                (displs.mean(), displs[max_pos], max_pos, displs[min_pos], min_pos))

        self.history.info("Changing structure (only significant diffs are shown):")
        if not recs:
            self.history.info("Input and output structure seems to be equal within the given tolerances")
        else:
            for rec in recs:
                self.history.info(rec)

        self.input.set_structure(new_structure)
        #assert self.input.structure == new_structure

    def autoparal_run(self):
        """
        Find an optimal set of parameters for the execution of the task
        This method can change the ABINIT input variables and/or the
        submission parameters e.g. the number of CPUs for MPI and OpenMp.

        Set:
           self.pconfs where pconfs is a :class:`ParalHints` object with the configuration reported by
           autoparal and optimal is the optimal configuration selected.
           Returns 0 if success
        """
        policy = self.manager.policy

        if policy.autoparal == 0: # or policy.max_ncpus in [None, 1]:
            logger.info("Nothing to do in autoparal, returning (None, None)")
            return 0

        if policy.autoparal != 1:
            raise NotImplementedError("autoparal != 1")

        ############################################################################
        # Run ABINIT in sequential to get the possible configurations with max_ncpus
        ############################################################################

        # Set the variables for automatic parallelization
        # Will get all the possible configurations up to max_ncpus
        # Return immediately if max_ncpus == 1
        max_ncpus = self.manager.max_cores
        if max_ncpus == 1: return 0

        autoparal_vars = dict(autoparal=policy.autoparal, max_ncpus=max_ncpus, mem_test=0)
        self.set_vars(autoparal_vars)

        # Run the job in a shell subprocess with mpi_procs = 1
        # we don't want to make a request to the queue manager for this simple job!
        # Return code is always != 0
        process = self.manager.to_shell_manager(mpi_procs=1).launch(self)
        self.history.pop()
        retcode = process.wait()
        # To avoid: ResourceWarning: unclosed file <_io.BufferedReader name=87> in py3k
        process.stderr.close()
        #process.stdout.close()

        # Remove the variables added for the automatic parallelization
        self.input.remove_vars(list(autoparal_vars.keys()))

        ##############################################################
        # Parse the autoparal configurations from the main output file
        ##############################################################
        parser = ParalHintsParser()
        try:
            pconfs = parser.parse(self.output_file.path)
        except parser.Error:
            # In principle Abinit should have written a complete log file
            # because we called .wait() but sometimes the Yaml doc is incomplete and
            # the parser raises. Let's wait 5 secs and then try again.
            time.sleep(5)
            try:
                pconfs = parser.parse(self.output_file.path)
            except parser.Error:
                logger.critical("Error while parsing Autoparal section:\n%s" % straceback())
                return 2

        ######################################################
        # Select the optimal configuration according to policy
        ######################################################
        optconf = self.find_optconf(pconfs)

        ####################################################
        # Change the input file and/or the submission script
        ####################################################
        self.set_vars(optconf.vars)

        # Write autoparal configurations to JSON file.
        d = pconfs.as_dict()
        d["optimal_conf"] = optconf
        json_pretty_dump(d, os.path.join(self.workdir, "autoparal.json"))

        ##############
        # Finalization
        ##############
        # Reset the status, remove garbage files ...
        self.set_status(self.S_INIT, msg='finished autoparal run')

        # Remove the output file since Abinit likes to create new files
        # with extension .outA, .outB if the file already exists.
        os.remove(self.output_file.path)
        os.remove(self.log_file.path)
        os.remove(self.stderr_file.path)

        return 0

    def find_optconf(self, pconfs):
        """Find the optimal Parallel configuration."""
        # Save pconfs for future reference.
        self.set_pconfs(pconfs)

        # Select the partition on which we'll be running and set MPI/OMP cores.
        optconf = self.manager.select_qadapter(pconfs)
        return optconf

    def select_files(self, what="o"):
        """
        Helper function used to select the files of a task.

        Args:
            what: string with the list of characters selecting the file type
                  Possible choices:
                    i ==> input_file,
                    o ==> output_file,
                    f ==> files_file,
                    j ==> job_file,
                    l ==> log_file,
                    e ==> stderr_file,
                    q ==> qout_file,
                    all ==> all files.
        """
        choices = collections.OrderedDict([
            ("i", self.input_file),
            ("o", self.output_file),
            ("f", self.files_file),
            ("j", self.job_file),
            ("l", self.log_file),
            ("e", self.stderr_file),
            ("q", self.qout_file),
        ])

        if what == "all":
            return [getattr(v, "path") for v in choices.values()]

        selected = []
        for c in what:
            try:
                selected.append(getattr(choices[c], "path"))
            except KeyError:
                logger.warning("Wrong keyword %s" % c)

        return selected

    def restart(self):
        """
        general restart used when scheduler problems have been taken care of
        """
        return self._restart()

    #@check_spectator
    def reset_from_scratch(self):
        """
        Restart from scratch, this is to be used if a job is restarted with more resources after a crash

        Move output files produced in workdir to _reset otherwise check_status continues
        to see the task as crashed even if the job did not run
        """
        # Create reset directory if not already done.
        reset_dir = os.path.join(self.workdir, "_reset")
        reset_file = os.path.join(reset_dir, "_counter")
        if not os.path.exists(reset_dir):
            os.mkdir(reset_dir)
            num_reset = 1
        else:
            with open(reset_file, "rt") as fh:
                num_reset = 1 + int(fh.read())

        # Move files to reset and append digit with reset index.
        def move_file(f):
            if not f.exists: return
            try:
                f.move(os.path.join(reset_dir, f.basename + "_" + str(num_reset)))
            except OSError as exc:
                logger.warning("Couldn't move file {}. exc: {}".format(f, str(exc)))

        for fname in ("output_file", "log_file", "stderr_file", "qout_file", "qerr_file"):
            move_file(getattr(self, fname))

        with open(reset_file, "wt") as fh:
            fh.write(str(num_reset))

        self.start_lockfile.remove()

        # Reset datetimes
        self.datetimes.reset()

        return self._restart(submit=False)

    #@check_spectator
    def fix_abicritical(self):
        """
        method to fix crashes/error caused by abinit

        Returns:
            1 if task has been fixed else 0.
        """
        event_handlers = self.event_handlers
        if not event_handlers:
            self.set_status(status=self.S_ERROR, msg='Empty list of event handlers. Cannot fix abi_critical errors')
            return 0

        count, done = 0, len(event_handlers) * [0]

        report = self.get_event_report()
        if report is None:
            self.set_status(status=self.S_ERROR, msg='get_event_report returned None')
            return 0

        # Note we have loop over all possible events (slow, I know)
        # because we can have handlers for Error, Bug or Warning
        # (ideally only for CriticalWarnings but this is not done yet)
        for event in report:
            for i, handler in enumerate(self.event_handlers):

                if handler.can_handle(event) and not done[i]:
                    logger.info("handler %s will try to fix event %s" % (handler, event))
                    try:
                        d = handler.handle_task_event(self, event)
                        if d:
                            done[i] += 1
                            count += 1

                    except Exception as exc:
                        logger.critical(str(exc))

        if count:
            self.reset_from_scratch()
            return 1

        self.set_status(status=self.S_ERROR, msg='We encountered AbiCritical events that could not be fixed')
        return 0

    #@check_spectator
    def fix_queue_critical(self):
        """
        This function tries to fix critical events originating from the queue submission system.

        General strategy, first try to increase resources in order to fix the problem,
        if this is not possible, call a task specific method to attempt to decrease the demands.

        Returns:
            1 if task has been fixed else 0.
        """
        from pymatgen.io.abinit.scheduler_error_parsers import NodeFailureError, MemoryCancelError, TimeCancelError
        #assert isinstance(self.manager, TaskManager)

        self.history.info('fixing queue critical')
        ret = "task.fix_queue_critical: "

        if not self.queue_errors:
            # TODO
            # paral_kgb = 1 leads to nasty sigegv that are seen as Qcritical errors!
            # Try to fallback to the conjugate gradient.
            #if self.uses_paral_kgb(1):
            #    logger.critical("QCRITICAL with PARAL_KGB==1. Will try CG!")
            #    self.set_vars(paral_kgb=0)
            #    self.reset_from_scratch()
            #    return
            # queue error but no errors detected, try to solve by increasing ncpus if the task scales
            # if resources are at maximum the task is definitively turned to errored
            if self.mem_scales or self.load_scales:
                try:
                    self.manager.increase_resources()  # acts either on the policy or on the qadapter
                    self.reset_from_scratch()
                    ret += "increased resources"
                    return ret
                except ManagerIncreaseError:
                    self.set_status(self.S_ERROR, msg='unknown queue error, could not increase resources any further')
                    raise FixQueueCriticalError
            else:
                self.set_status(self.S_ERROR, msg='unknown queue error, no options left')
                raise FixQueueCriticalError

        else:
            print("Fix_qcritical: received %d queue_errors" % len(self.queue_errors))
            print("type_list: %s" % list(type(qe) for qe in self.queue_errors))

            for error in self.queue_errors:
                self.history.info('fixing: %s' % str(error))
                ret += str(error)
                if isinstance(error, NodeFailureError):
                    # if the problematic node is known, exclude it
                    if error.nodes is not None:
                        try:
                            self.manager.exclude_nodes(error.nodes)
                            self.reset_from_scratch()
                            self.set_status(self.S_READY, msg='excluding nodes')
                        except:
                            raise FixQueueCriticalError
                    else:
                        self.set_status(self.S_ERROR, msg='Node error but no node identified.')
                        raise FixQueueCriticalError

                elif isinstance(error, MemoryCancelError):
                    # ask the qadapter to provide more resources, i.e. more cpu's so more total memory if the code
                    # scales this should fix the memeory problem
                    # increase both max and min ncpu of the autoparalel and rerun autoparalel
                    if self.mem_scales:
                        try:
                            self.manager.increase_ncpus()
                            self.reset_from_scratch()
                            self.set_status(self.S_READY, msg='increased ncps to solve memory problem')
                            return
                        except ManagerIncreaseError:
                            self.history.warning('increasing ncpus failed')

                    # if the max is reached, try to increase the memory per cpu:
                    try:
                        self.manager.increase_mem()
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='increased mem')
                        return
                    except ManagerIncreaseError:
                        self.history.warning('increasing mem failed')

                    # if this failed ask the task to provide a method to reduce the memory demand
                    try:
                        self.reduce_memory_demand()
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='decreased mem demand')
                        return
                    except DecreaseDemandsError:
                        self.history.warning('decreasing demands failed')

                    msg = ('Memory error detected but the memory could not be increased neither could the\n'
                           'memory demand be decreased. Unrecoverable error.')
                    self.set_status(self.S_ERROR, msg)
                    raise FixQueueCriticalError

                elif isinstance(error, TimeCancelError):
                    # ask the qadapter to provide more time
                    print('trying to increase time')
                    try:
                        self.manager.increase_time()
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='increased wall time')
                        return
                    except ManagerIncreaseError:
                        self.history.warning('increasing the waltime failed')

                    # if this fails ask the qadapter to increase the number of cpus
                    if self.load_scales:
                        try:
                            self.manager.increase_ncpus()
                            self.reset_from_scratch()
                            self.set_status(self.S_READY, msg='increased number of cpus')
                            return
                        except ManagerIncreaseError:
                            self.history.warning('increase ncpus to speed up the calculation to stay in the walltime failed')

                    # if this failed ask the task to provide a method to speed up the task
                    try:
                        self.speed_up()
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='task speedup')
                        return
                    except DecreaseDemandsError:
                        self.history.warning('decreasing demands failed')

                    msg = ('Time cancel error detected but the time could not be increased neither could\n'
                           'the time demand be decreased by speedup of increasing the number of cpus.\n'
                           'Unrecoverable error.')
                    self.set_status(self.S_ERROR, msg)

                else:
                    msg = 'No solution provided for error %s. Unrecoverable error.' % error.name
                    self.set_status(self.S_ERROR, msg)

        return 0

    def parse_timing(self):
        """
        Parse the timer data in the main output file of Abinit.
        Requires timopt /= 0 in the input file (usually timopt = -1)

        Return: :class:`AbinitTimerParser` instance, None if error.
        """
        from .abitimer import AbinitTimerParser
        parser = AbinitTimerParser()
        read_ok = parser.parse(self.output_file.path)
        if read_ok:
            return parser
        return None


class ProduceHist:
    """
    Mixin class for an :class:`AbinitTask` producing a HIST file.
    Provide the method `open_hist` that reads and return a HIST file.
    """
    @property
    def hist_path(self):
        """Absolute path of the HIST file. Empty string if file is not present."""
        # Lazy property to avoid multiple calls to has_abiext.
        try:
            return self._hist_path
        except AttributeError:
            path = self.outdir.has_abiext("HIST")
            if path: self._hist_path = path
            return path

    def open_hist(self):
        """
        Open the HIST file located in the in self.outdir.
        Returns :class:`HistFile` object, None if file could not be found or file is not readable.
        """
        if not self.hist_path:
            if self.status == self.S_OK:
                logger.critical("%s reached S_OK but didn't produce a HIST file in %s" % (self, self.outdir))
            return None

        # Open the HIST file
        from abipy.dynamics.hist import HistFile
        try:
            return HistFile(self.hist_path)
        except Exception as exc:
            logger.critical("Exception while reading HIST file at %s:\n%s" % (self.hist_path, str(exc)))
            return None


class GsTask(AbinitTask):
    """
    Base class for ground-state tasks. A ground state task produces a GSR file
    Provides the method `open_gsr` that reads and returns a GSR file.
    """
    @property
    def gsr_path(self):
        """Absolute path of the GSR file. Empty string if file is not present."""
        # Lazy property to avoid multiple calls to has_abiext.
        try:
            return self._gsr_path
        except AttributeError:
            path = self.outdir.has_abiext("GSR")
            if path: self._gsr_path = path
            return path

    def open_gsr(self):
        """
        Open the GSR file located in the in self.outdir.
        Returns :class:`GsrFile` object, None if file could not be found or file is not readable.
        """
        gsr_path = self.gsr_path
        if not gsr_path:
            if self.status == self.S_OK:
                logger.critical("%s reached S_OK but didn't produce a GSR file in %s" % (self, self.outdir))
            return None

        # Open the GSR file.
        from abipy.electrons.gsr import GsrFile
        try:
            return GsrFile(gsr_path)
        except Exception as exc:
            logger.critical("Exception while reading GSR file at %s:\n%s" % (gsr_path, str(exc)))
            return None


class ScfTask(GsTask):
    """
    Self-consistent ground-state calculations.
    Provide support for in-place restart via (WFK|DEN) files
    """
    CRITICAL_EVENTS = [
        events.ScfConvergenceWarning,
    ]

    color_rgb = np.array((255, 0, 0)) / 255

    def restart(self):
        """SCF calculations can be restarted if we have either the WFK file or the DEN file."""
        # Prefer WFK over DEN files since we can reuse the wavefunctions.
        for ext in ("WFK", "DEN"):
            restart_file = self.outdir.has_abiext(ext)
            irdvars = irdvars_for_ext(ext)
            if restart_file: break
        else:
            raise self.RestartError("%s: Cannot find WFK or DEN file to restart from." % self)

        # Move out --> in.
        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        self.set_vars(irdvars)

        # Now we can resubmit the job.
        self.history.info("Will restart from %s", restart_file)
        return self._restart()

    def inspect(self, **kwargs):
        """
        Plot the SCF cycle results with matplotlib.

        Returns
            `matplotlib` figure, None if some error occurred.
        """
        try:
            scf_cycle = abiinspect.GroundStateScfCycle.from_file(self.output_file.path)
        except IOError:
            return None

        if scf_cycle is not None:
            if "title" not in kwargs: kwargs["title"] = str(self)
            return scf_cycle.plot(**kwargs)

        return None

    def get_results(self, **kwargs):
        results = super(ScfTask, self).get_results(**kwargs)

        # Open the GSR file and add its data to results.out
        with self.open_gsr() as gsr:
            results["out"].update(gsr.as_dict())
            # Add files to GridFS
            results.register_gridfs_files(GSR=gsr.filepath)

        return results


class CollinearThenNonCollinearScfTask(ScfTask):
    """
    A specialized ScfTaks that performs an initial SCF run with nsppol = 2.
    The spin polarized WFK file is then used to start a non-collinear SCF run (nspinor == 2)
    initialized from the previous WFK file.
    """
    def __init__(self, input, workdir=None, manager=None, deps=None):
        super(CollinearThenNonCollinearScfTask, self).__init__(input, workdir=workdir, manager=manager, deps=deps)
        # Enforce nspinor = 1, nsppol = 2 and prtwf = 1.
        self._input = self.input.deepcopy()
        self.input.set_spin_mode("polarized")
        self.input.set_vars(prtwf=1)
        self.collinear_done = False

    def _on_ok(self):
        results = super(CollinearThenNonCollinearScfTask, self)._on_ok()
        if not self.collinear_done:
            self.input.set_spin_mode("spinor")
            self.collinear_done = True
            self.finalized = False
            self.restart()

        return results


class NscfTask(GsTask):
    """
    Non-Self-consistent GS calculation. Provide in-place restart via WFK files
    """
    CRITICAL_EVENTS = [
        events.NscfConvergenceWarning,
    ]

    color_rgb = np.array((200, 80, 100)) / 255

    def setup(self):
        """
        NSCF calculations should use the same FFT mesh as the one employed in the GS task
        (in principle, it's possible to interpolate inside Abinit but tests revealed some numerical noise
        Here we change the input file of the NSCF task to have the same FFT mesh.
        """
        for dep in self.deps:
            if "DEN" in dep.exts:
                parent_task = dep.node
                break
        else:
            raise RuntimeError("Cannot find parent node producing DEN file")

        with parent_task.open_gsr() as gsr:
            den_mesh = 3 * [None]
            den_mesh[0] = gsr.reader.read_dimvalue("number_of_grid_points_vector1")
            den_mesh[1] = gsr.reader.read_dimvalue("number_of_grid_points_vector2")
            den_mesh[2] = gsr.reader.read_dimvalue("number_of_grid_points_vector3")
            if self.ispaw:
                self.set_vars(ngfftdg=den_mesh)
            else:
                self.set_vars(ngfft=den_mesh)

        super(NscfTask, self).setup()

    def restart(self):
        """NSCF calculations can be restarted only if we have the WFK file."""
        ext = "WFK"
        restart_file = self.outdir.has_abiext(ext)
        if not restart_file:
            raise self.RestartError("%s: Cannot find the WFK file to restart from." % self)

        # Move out --> in.
        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        irdvars = irdvars_for_ext(ext)
        self.set_vars(irdvars)

        # Now we can resubmit the job.
        self.history.info("Will restart from %s", restart_file)
        return self._restart()

    def get_results(self, **kwargs):
        results = super(NscfTask, self).get_results(**kwargs)

        # Read the GSR file.
        with self.open_gsr() as gsr:
            results["out"].update(gsr.as_dict())
            # Add files to GridFS
            results.register_gridfs_files(GSR=gsr.filepath)

        return results


class RelaxTask(GsTask, ProduceHist):
    """
    Task for structural optimizations.
    """
    # TODO possible ScfConvergenceWarning?
    CRITICAL_EVENTS = [
        events.RelaxConvergenceWarning,
    ]

    color_rgb = np.array((255, 61, 255)) / 255

    def get_final_structure(self):
        """Read the final structure from the GSR file."""
        try:
            with self.open_gsr() as gsr:
                return gsr.structure
        except AttributeError:
            raise RuntimeError("Cannot find the GSR file with the final structure to restart from.")

    def restart(self):
        """
        Restart the structural relaxation.

        Structure relaxations can be restarted only if we have the WFK file or the DEN or the GSR file
        from which we can read the last structure (mandatory) and the wavefunctions (not mandatory but useful).
        Prefer WFK over other files since we can reuse the wavefunctions.

        .. note::

            The problem in the present approach is that some parameters in the input
            are computed from the initial structure and may not be consistent with
            the modification of the structure done during the structure relaxation.
        """
        restart_file = None

        # Try to restart from the WFK file if possible.
        # FIXME: This part has been disabled because WFK=IO is a mess if paral_kgb == 1
        # This is also the reason why I wrote my own MPI-IO code for the GW part!
        wfk_file = self.outdir.has_abiext("WFK")
        if False and wfk_file:
            irdvars = irdvars_for_ext("WFK")
            restart_file = self.out_to_in(wfk_file)

        # Fallback to DEN file. Note that here we look for out_DEN instead of out_TIM?_DEN
        # This happens when the previous run completed and task.on_done has been performed.
        # ********************************************************************************
        # Note that it's possible to have an undetected error if we have multiple restarts
        # and the last relax died badly. In this case indeed out_DEN is the file produced
        # by the last run that has executed on_done.
        # ********************************************************************************
        if restart_file is None:
            for ext in ("", ".nc"):
                out_den = self.outdir.path_in("out_DEN" + ext)
                if os.path.exists(out_den):
                    irdvars = irdvars_for_ext("DEN")
                    restart_file = self.out_to_in(out_den)
                    break

        if restart_file is None:
            # Try to restart from the last TIM?_DEN file.
            # This should happen if the previous run didn't complete in clean way.
            # Find the last TIM?_DEN file.
            last_timden = self.outdir.find_last_timden_file()
            if last_timden is not None:
                if last_timden.path.endswith(".nc"):
                    ofile = self.outdir.path_in("out_DEN.nc")
                else:
                    ofile = self.outdir.path_in("out_DEN")
                os.rename(last_timden.path, ofile)
                restart_file = self.out_to_in(ofile)
                irdvars = irdvars_for_ext("DEN")

        if restart_file is None:
            # Don't raise RestartError as we can still change the structure.
            self.history.warning("Cannot find the WFK|DEN|TIM?_DEN file to restart from.")
        else:
            # Add the appropriate variable for restarting.
            self.set_vars(irdvars)
            self.history.info("Will restart from %s", restart_file)

        # FIXME Here we should read the HIST file but restartxf if broken!
        #self.set_vars({"restartxf": -1})

        # Read the relaxed structure from the GSR file and change the input.
        self._change_structure(self.get_final_structure())

        # Now we can resubmit the job.
        return self._restart()

    def inspect(self, **kwargs):
        """
        Plot the evolution of the structural relaxation with matplotlib.

        Args:
            what: Either "hist" or "scf". The first option (default) extracts data
                from the HIST file and plot the evolution of the structural
                parameters, forces, pressures and energies.
                The second option, extracts data from the main output file and
                plot the evolution of the SCF cycles (etotal, residuals, etc).

        Returns:
            `matplotlib` figure, None if some error occurred.
        """
        what = kwargs.pop("what", "hist")

        if what == "hist":
            # Read the hist file to get access to the structure.
            with self.open_hist() as hist:
                return hist.plot(**kwargs) if hist else None

        elif what == "scf":
            # Get info on the different SCF cycles
            relaxation = abiinspect.Relaxation.from_file(self.output_file.path)
            if "title" not in kwargs: kwargs["title"] = str(self)
            return relaxation.plot(**kwargs) if relaxation is not None else None

        else:
            raise ValueError("Wrong value for what %s" % what)

    def get_results(self, **kwargs):
        results = super(RelaxTask, self).get_results(**kwargs)

        # Open the GSR file and add its data to results.out
        with self.open_gsr() as gsr:
            results["out"].update(gsr.as_dict())
            # Add files to GridFS
            results.register_gridfs_files(GSR=gsr.filepath)

        return results

    def reduce_dilatmx(self, target=1.01):
        actual_dilatmx = self.get_inpvar('dilatmx', 1.)
        new_dilatmx = actual_dilatmx - min((actual_dilatmx-target), actual_dilatmx*0.05)
        self.set_vars(dilatmx=new_dilatmx)

    def fix_ofiles(self):
        """
        Note that ABINIT produces lots of out_TIM1_DEN files for each step.
        Here we list all TIM*_DEN files, we select the last one and we rename it in out_DEN

        This change is needed so that we can specify dependencies with the syntax {node: "DEN"}
        without having to know the number of iterations needed to converge the run in node!
        """
        super(RelaxTask, self).fix_ofiles()

        # Find the last TIM?_DEN file.
        last_timden = self.outdir.find_last_timden_file()
        if last_timden is None:
            logger.warning("Cannot find TIM?_DEN files")
            return

        # Rename last TIMDEN with out_DEN.
        ofile = self.outdir.path_in("out_DEN")
        if last_timden.path.endswith(".nc"): ofile += ".nc"
        self.history.info("Renaming last_denfile %s --> %s" % (last_timden.path, ofile))
        os.rename(last_timden.path, ofile)


class DfptTask(AbinitTask):
    """
    Base class for DFPT tasks (Phonons, DdeTask, DdkTask, ElasticTask ...)
    Mainly used to implement methods that are common to DFPT calculations with Abinit.
    Provide the method `open_ddb` that reads and return a Ddb file.

    .. warning::

        This class should not be instantiated directly.
    """
    # TODO:
    # for the time being we don't discern between GS and PhononCalculations.
    CRITICAL_EVENTS = [
        events.ScfConvergenceWarning,
    ]

    def __repr__(self):
        # Get info about DFT perturbation from input file.
        qpt = self.input.get("qpt", [0, 0, 0])
        rfphon = self.input.get("rfphon", 0)
        rfatpol = self.input.get("rfatpol", [1, 1])
        rfelfd = self.input.get("rfelfd", 0)
        rfstrs = self.input.get("rfstrs", 0)
        rfdir = self.input.get("rfdir", [0, 0, 0])
        irdddk = self.input.get("irdddk", 0)

        dfpt_info = ""
        if rfphon != 0:
            dfpt_info = "qpt: {}, rfphon: {}, rfatpol: {}, rfdir: {}, irdddk: {}".format(
                    qpt, rfphon, rfatpol, rfdir, irdddk)

        elif rfelfd != 0:
            dfpt_info = "qpt: {}, rfelfd: {} rfdir: {}, irdddk: {}".format(
                    qpt, rfelfd, rfdir, irdddk)

        elif rfstrs != 0:
            dfpt_info = "qpt: {}, rfstrs: {}, rfdir: {}, irdddk: {}".format(
                    qpt, rfstrs, rfdir, irdddk)

        try:
            return "<%s, node_id=%s, workdir=%s, %s>" % (
                self.__class__.__name__, self.node_id, self.relworkdir, dfpt_info)
        except AttributeError:
            # this usually happens when workdir has not been initialized
            return "<%s, node_id=%s, workdir=None, %s>" % (
                self.__class__.__name__, self.node_id, dfpt_info)

    @property
    def ddb_path(self):
        """Absolute path of the DDB file. Empty string if file is not present."""
        # Lazy property to avoid multiple calls to has_abiext.
        try:
            return self._ddb_path
        except AttributeError:
            path = self.outdir.has_abiext("DDB")
            if path: self._ddb_path = path
            return path

    def open_ddb(self):
        """
        Open the DDB file located in the in self.outdir.
        Returns :class:`DdbFile` object, None if file could not be found or file is not readable.
        """
        ddb_path = self.ddb_path
        if not ddb_path:
            if self.status == self.S_OK:
                logger.critical("%s reached S_OK but didn't produce a DDB file in %s" % (self, self.outdir))
            return None

        # Open the DDB file.
        from abipy.dfpt.ddb import DdbFile
        try:
            return DdbFile(ddb_path)
        except Exception as exc:
            logger.critical("Exception while reading DDB file at %s:\n%s" % (ddb_path, str(exc)))
            return None

    def make_links(self):
        """
        Replace the default behaviour of make_links. More specifically, this method
        implements the logic required to connect DFPT calculation to `DDK` files.
        Remember that DDK is an extension introduced in AbiPy to deal with the
        irdddk input variable and the fact that the 3 files with du/dk produced by Abinit
        have a file extension constructed from the number of atom (e.g. 1WF[3natom +1]).

        AbiPy uses the user-friendly syntax deps={node: "DDK"} to specify that
        the children will read the DDK from `node` but this also means that
        we have to implement extract logic to handle this case at runtime.
        """
        for dep in self.deps:
            for d in dep.exts:

                if d == "DDK":
                    ddk_task = dep.node
                    out_ddk = ddk_task.outdir.has_abiext("DDK")
                    if not out_ddk:
                        raise RuntimeError("%s didn't produce the DDK file" % ddk_task)

                    # Get (fortran) idir and costruct the name of the 1WF expected by Abinit
                    rfdir = list(ddk_task.input["rfdir"])
                    if rfdir.count(1) != 1:
                        raise RuntimeError("Only one direction should be specifned in rfdir but rfdir = %s" % rfdir)

                    idir = rfdir.index(1) + 1
                    ddk_case = idir + 3 * len(ddk_task.input.structure)

                    infile = self.indir.path_in("in_1WF%d" % ddk_case)
                    os.symlink(out_ddk, infile)

                elif d in ("WFK", "WFQ"):
                    gs_task = dep.node
                    out_wfk = gs_task.outdir.has_abiext(d)
                    if not out_wfk:
                        raise RuntimeError("%s didn't produce the %s file" % (gs_task, d))

                    if d == "WFK":
                        bname = "in_WFK"
                    elif d == "WFQ":
                        bname = "in_WFQ"
                    else:
                        raise ValueError("Don't know how to handle `%s`" % d)

                    if not os.path.exists(self.indir.path_in(bname)):
                            os.symlink(out_wfk, self.indir.path_in(bname))

                elif d == "DEN":
                    gs_task = dep.node
                    out_wfk = gs_task.outdir.has_abiext("DEN")
                    if not out_wfk:
                        raise RuntimeError("%s didn't produce the DEN file" % gs_task)
                    if not os.path.exists(self.indir.path_in("in_DEN")):
                        os.symlink(out_wfk, self.indir.path_in("in_DEN"))

                elif d == "1WF":
                    gs_task = dep.node
                    out_wfk = gs_task.outdir.has_abiext("1WF")
                    if not out_wfk:
                        raise RuntimeError("%s didn't produce the 1WF file" % gs_task)
                    dest = self.indir.path_in("in_" + out_wfk.split("_")[-1])
                    if not os.path.exists(dest):
                        os.symlink(out_wfk, dest)

                elif d == "1DEN":
                    gs_task = dep.node
                    out_wfk = gs_task.outdir.has_abiext("DEN")
                    if not out_wfk:
                        raise RuntimeError("%s didn't produce the 1DEN file" % gs_task)
                    dest = self.indir.path_in("in_" + out_wfk.split("_")[-1])
                    if not os.path.exists(dest):
                        os.symlink(out_wfk, dest)

                else:
                    raise ValueError("Don't know how to handle extension: %s" % str(dep.exts))

    def restart(self):
        """
        DFPT calculations can be restarted only if we have the 1WF file or the 1DEN file.
        from which we can read the first-order wavefunctions or the first order density.
        Prefer 1WF over 1DEN since we can reuse the wavefunctions.
        """
        # Abinit adds the idir-ipert index at the end of the file and this breaks the extension
        # e.g. out_1WF4, out_DEN4. find_1wf_files and find_1den_files returns the list of files found
        restart_file, irdvars = None, None

        # Highest priority to the 1WF file because restart is more efficient.
        wf_files = self.outdir.find_1wf_files()
        if wf_files is not None:
            restart_file = wf_files[0].path
            irdvars = irdvars_for_ext("1WF")
            if len(wf_files) != 1:
                restart_file = None
                self.history.critical("Found more than one 1WF file in outdir. Restart is ambiguous!")

        if restart_file is None:
            den_files = self.outdir.find_1den_files()
            if den_files is not None:
                restart_file = den_files[0].path
                irdvars = {"ird1den": 1}
                if len(den_files) != 1:
                    restart_file = None
                    self.history.critical("Found more than one 1DEN file in outdir. Restart is ambiguous!")

        if restart_file is None:
            # Raise because otherwise restart is equivalent to a run from scratch --> infinite loop!
            raise self.RestartError("%s: Cannot find the 1WF|1DEN file to restart from." % self)

        # Move file.
        self.history.info("Will restart from %s", restart_file)
        restart_file = self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        self.set_vars(irdvars)

        # Now we can resubmit the job.
        return self._restart()


class DdeTask(DfptTask):
    """Task for DDE calculations (perturbation wrt electric field)."""

    color_rgb = np.array((61, 158, 255)) / 255

    def get_results(self, **kwargs):
        results = super(DdeTask, self).get_results(**kwargs)
        return results.register_gridfs_file(DDB=(self.outdir.has_abiext("DDE"), "t"))


class DteTask(DfptTask):
    """Task for DTE calculations."""
    color_rgb = np.array((204, 0, 204)) / 255

    # @check_spectator
    def start(self, **kwargs):
        kwargs['autoparal'] = False
        return super(DteTask, self).start(**kwargs)

    def get_results(self, **kwargs):
        results = super(DteTask, self).get_results(**kwargs)
        return results.register_gridfs_file(DDB=(self.outdir.has_abiext("DDE"), "t"))


class DdkTask(DfptTask):
    """Task for DDK calculations."""
    color_rgb = np.array((0, 204, 204)) / 255

    #@check_spectator
    def _on_ok(self):
        super(DdkTask, self)._on_ok()
        # Client code expects to find du/dk in DDK file.
        # Here I create a symbolic link out_1WF13 --> out_DDK
        # so that we can use deps={ddk_task: "DDK"} in the high-level API.
        # The price to pay is that we have to handle the DDK extension in make_links.
        # See DfptTask.make_links
        self.outdir.symlink_abiext('1WF', 'DDK')

    def get_results(self, **kwargs):
        results = super(DdkTask, self).get_results(**kwargs)
        return results.register_gridfs_file(DDK=(self.outdir.has_abiext("DDK"), "t"))


class BecTask(DfptTask):
    """
    Task for the calculation of Born effective charges.

    bec_deps = {ddk_task: "DDK" for ddk_task in ddk_tasks}
    bec_deps.update({scf_task: "WFK"})
    """
    color_rgb = np.array((122, 122, 255)) / 255


class PhononTask(DfptTask):
    """
    DFPT calculations for a single atomic perturbation.
    Provide support for in-place restart via (1WF|1DEN) files
    """
    color_rgb = np.array((0, 150, 250)) / 255

    def inspect(self, **kwargs):
        """
        Plot the Phonon SCF cycle results with matplotlib.

        Returns:
            `matplotlib` figure, None if some error occurred.
        """
        scf_cycle = abiinspect.PhononScfCycle.from_file(self.output_file.path)
        if scf_cycle is not None:
            if "title" not in kwargs: kwargs["title"] = str(self)
            return scf_cycle.plot(**kwargs)

    def get_results(self, **kwargs):
        results = super(PhononTask, self).get_results(**kwargs)
        return results.register_gridfs_files(DDB=(self.outdir.has_abiext("DDB"), "t"))


class ElasticTask(DfptTask):
    """
    DFPT calculations for a single strain perturbation (uniaxial or shear strain).
    Provide support for in-place restart via (1WF|1DEN) files
    """
    color_rgb = np.array((255, 204, 255)) / 255


class EphTask(AbinitTask):
    """
    Class for electron-phonon calculations.
    """
    color_rgb = np.array((255, 128, 0)) / 255


class ManyBodyTask(AbinitTask):
    """
    Base class for Many-body tasks (Screening, Sigma, Bethe-Salpeter)
    Mainly used to implement methods that are common to MBPT calculations with Abinit.

    .. warning::

        This class should not be instantiated directly.
    """
    def reduce_memory_demand(self):
        """
        Method that can be called by the scheduler to decrease the memory demand of a specific task.
        Returns True in case of success, False in case of Failure.
        """
        # The first digit governs the storage of W(q), the second digit the storage of u(r)
        # Try to avoid the storage of u(r) first since reading W(q) from file will lead to a drammatic slowdown.
        prev_gwmem = int(self.get_inpvar("gwmem", default=11))
        first_dig, second_dig = prev_gwmem // 10, prev_gwmem % 10

        if second_dig == 1:
            self.set_vars(gwmem="%.2d" % (10 * first_dig))
            return True

        if first_dig == 1:
            self.set_vars(gwmem="%.2d" % 00)
            return True

        # gwmem 00 d'oh!
        return False


class ScrTask(ManyBodyTask):
    """Tasks for SCREENING calculations """

    color_rgb = np.array((255, 128, 0)) / 255

    #def inspect(self, **kwargs):
    #    """Plot graph showing the number of q-points computed and the wall-time used"""

    @property
    def scr_path(self):
        """Absolute path of the SCR file. Empty string if file is not present."""
        # Lazy property to avoid multiple calls to has_abiext.
        try:
            return self._scr_path
        except AttributeError:
            path = self.outdir.has_abiext("SCR.nc")
            if path: self._scr_path = path
            return path

    def open_scr(self):
        """
        Open the SIGRES file located in the in self.outdir.
        Returns :class:`ScrFile` object, None if file could not be found or file is not readable.
        """
        scr_path = self.scr_path

        if not scr_path:
            logger.critical("%s didn't produce a SCR.nc file in %s" % (self, self.outdir))
            return None

        # Open the GSR file and add its data to results.out
        from abipy.electrons.scr import ScrFile
        try:
            return ScrFile(scr_path)
        except Exception as exc:
            logger.critical("Exception while reading SCR file at %s:\n%s" % (scr_path, str(exc)))
            return None


class SigmaTask(ManyBodyTask):
    """
    Tasks for SIGMA calculations. Provides support for in-place restart via QPS files
    """
    CRITICAL_EVENTS = [
        events.QPSConvergenceWarning,
    ]

    color_rgb = np.array((0, 255, 0)) / 255

    def restart(self):
        # G calculations can be restarted only if we have the QPS file
        # from which we can read the results of the previous step.
        ext = "QPS"
        restart_file = self.outdir.has_abiext(ext)
        if not restart_file:
            raise self.RestartError("%s: Cannot find the QPS file to restart from." % self)

        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        irdvars = irdvars_for_ext(ext)
        self.set_vars(irdvars)

        # Now we can resubmit the job.
        self.history.info("Will restart from %s", restart_file)
        return self._restart()

    #def inspect(self, **kwargs):
    #    """Plot graph showing the number of k-points computed and the wall-time used"""

    @property
    def sigres_path(self):
        """Absolute path of the SIGRES file. Empty string if file is not present."""
        # Lazy property to avoid multiple calls to has_abiext.
        try:
            return self._sigres_path
        except AttributeError:
            path = self.outdir.has_abiext("SIGRES")
            if path: self._sigres_path = path
            return path

    def open_sigres(self):
        """
        Open the SIGRES file located in the in self.outdir.
        Returns :class:`SigresFile` object, None if file could not be found or file is not readable.
        """
        sigres_path = self.sigres_path

        if not sigres_path:
            logger.critical("%s didn't produce a SIGRES file in %s" % (self, self.outdir))
            return None

        # Open the SIGRES file and add its data to results.out
        from abipy.electrons.gw import SigresFile
        try:
            return SigresFile(sigres_path)
        except Exception as exc:
            logger.critical("Exception while reading SIGRES file at %s:\n%s" % (sigres_path, str(exc)))
            return None

    def get_scissors_builder(self):
        """
        Returns an instance of :class:`ScissorsBuilder` from the SIGRES file.

        Raise:
            `RuntimeError` if SIGRES file is not found.
        """
        from abipy.electrons.scissors import ScissorsBuilder
        if self.sigres_path:
            return ScissorsBuilder.from_file(self.sigres_path)
        else:
            raise RuntimeError("Cannot find SIGRES file!")

    def get_results(self, **kwargs):
        results = super(SigmaTask, self).get_results(**kwargs)

        # Open the SIGRES file and add its data to results.out
        with self.open_sigres() as sigres:
            #results["out"].update(sigres.as_dict())
            results.register_gridfs_files(SIGRES=sigres.filepath)

        return results


class BseTask(ManyBodyTask):
    """
    Task for Bethe-Salpeter calculations.

    .. note::

        The BSE codes provides both iterative and direct schemes for the computation of the dielectric function.
        The direct diagonalization cannot be restarted whereas Haydock and CG support restarting.
    """
    CRITICAL_EVENTS = [
        events.HaydockConvergenceWarning,
        #events.BseIterativeDiagoConvergenceWarning,
    ]

    color_rgb = np.array((128, 0, 255)) / 255

    def restart(self):
        """
        BSE calculations with Haydock can be restarted only if we have the
        excitonic Hamiltonian and the HAYDR_SAVE file.
        """
        # TODO: This version seems to work but the main output file is truncated
        # TODO: Handle restart if CG method is used
        # TODO: restart should receive a list of critical events
        # the log file is complete though.
        irdvars = {}

        # Move the BSE blocks to indata.
        # This is done only once at the end of the first run.
        # Successive restarts will use the BSR|BSC files in the indir directory
        # to initialize the excitonic Hamiltonian
        count = 0
        for ext in ("BSR", "BSC"):
            ofile = self.outdir.has_abiext(ext)
            if ofile:
                count += 1
                irdvars.update(irdvars_for_ext(ext))
                self.out_to_in(ofile)

        if not count:
            # outdir does not contain the BSR|BSC file.
            # This means that num_restart > 1 and the files should be in task.indir
            count = 0
            for ext in ("BSR", "BSC"):
                ifile = self.indir.has_abiext(ext)
                if ifile:
                    count += 1

            if not count:
                raise self.RestartError("%s: Cannot find BSR|BSC files in %s" % (self, self.indir))

        # Rename HAYDR_SAVE files
        count = 0
        for ext in ("HAYDR_SAVE", "HAYDC_SAVE"):
            ofile = self.outdir.has_abiext(ext)
            if ofile:
                count += 1
                irdvars.update(irdvars_for_ext(ext))
                self.out_to_in(ofile)

        if not count:
            raise self.RestartError("%s: Cannot find the HAYDR_SAVE file to restart from." % self)

        # Add the appropriate variable for restarting.
        self.set_vars(irdvars)

        # Now we can resubmit the job.
        #self.history.info("Will restart from %s", restart_file)
        return self._restart()

    #def inspect(self, **kwargs):
    #    """
    #    Plot the Haydock iterations with matplotlib.
    #
    #    Returns
    #        `matplotlib` figure, None if some error occurred.
    #    """
    #    haydock_cycle = abiinspect.HaydockIterations.from_file(self.output_file.path)
    #    if haydock_cycle is not None:
    #        if "title" not in kwargs: kwargs["title"] = str(self)
    #        return haydock_cycle.plot(**kwargs)

    @property
    def mdf_path(self):
        """Absolute path of the MDF file. Empty string if file is not present."""
        # Lazy property to avoid multiple calls to has_abiext.
        try:
            return self._mdf_path
        except AttributeError:
            path = self.outdir.has_abiext("MDF.nc")
            if path: self._mdf_path = path
            return path

    def open_mdf(self):
        """
        Open the MDF file located in the in self.outdir.
        Returns :class:`MdfFile` object, None if file could not be found or file is not readable.
        """
        mdf_path = self.mdf_path
        if not mdf_path:
            logger.critical("%s didn't produce a MDF file in %s" % (self, self.outdir))
            return None

        # Open the DFF file and add its data to results.out
        from abipy.electrons.bse import MdfFile
        try:
            return MdfFile(mdf_path)
        except Exception as exc:
            logger.critical("Exception while reading MDF file at %s:\n%s" % (mdf_path, str(exc)))
            return None

    def get_results(self, **kwargs):
        results = super(BseTask, self).get_results(**kwargs)

        with self.open_mdf() as mdf:
            #results["out"].update(mdf.as_dict())
            #epsilon_infinity optical_gap
            results.register_gridfs_files(MDF=mdf.filepath)

        return results


class OpticTask(Task):
    """
    Task for the computation of optical spectra with optic i.e.
    RPA without local-field effects and velocity operator computed from DDK files.
    """
    color_rgb = np.array((255, 204, 102)) / 255

    def __init__(self, optic_input, nscf_node, ddk_nodes, use_ddknc=False, workdir=None, manager=None):
        """
        Create an instance of :class:`OpticTask` from an string containing the input.

        Args:
            optic_input: :class:`OpticInput` object with optic variables.
            nscf_node: The task that will produce the WFK file with the KS energies or path to the WFK file.
            ddk_nodes: List of :class:`DdkTask` nodes that will produce the DDK files or list of DDK filepaths.
                Order (x, y, z)
            workdir: Path to the working directory.
            manager: :class:`TaskManager` object.
        """
        # Convert paths to FileNodes
        self.nscf_node = Node.as_node(nscf_node)
        self.ddk_nodes = [Node.as_node(n) for n in ddk_nodes]
        assert len(ddk_nodes) == 3
        #print(self.nscf_node, self.ddk_nodes)

        # Use DDK extension instead of 1WF
        if use_ddknc:
            deps = {n: "DDK.nc" for n in self.ddk_nodes}
        else:
            deps = {n: "1WF" for n in self.ddk_nodes}

        deps.update({self.nscf_node: "WFK"})

        super(OpticTask, self).__init__(optic_input, workdir=workdir, manager=manager, deps=deps)

    def set_workdir(self, workdir, chroot=False):
        """Set the working directory of the task."""
        super(OpticTask, self).set_workdir(workdir, chroot=chroot)
        # Small hack: the log file of optics is actually the main output file.
        self.output_file = self.log_file

    def set_vars(self, *args, **kwargs):
        """
        Optic does not use `get` or `ird` variables hence we should never try
        to change the input when we connect this task
        """
        kwargs.update(dict(*args))
        self.history.info("OpticTask intercepted set_vars with args %s" % kwargs)

        if "autoparal" in kwargs: self.input.set_vars(autoparal=kwargs["autoparal"])
        if "max_ncpus" in kwargs: self.input.set_vars(max_ncpus=kwargs["max_ncpus"])

    @property
    def executable(self):
        """Path to the executable required for running the :class:`OpticTask`."""
        try:
            return self._executable
        except AttributeError:
            return "optic"

    @property
    def filesfile_string(self):
        """String with the list of files and prefixes needed to execute ABINIT."""
        lines = []
        app = lines.append

        #optic.in     ! Name of input file
        #optic.out    ! Unused
        #optic        ! Root name for all files that will be produced
        app(self.input_file.path)                           # Path to the input file
        app(os.path.join(self.workdir, "unused"))           # Path to the output file
        app(os.path.join(self.workdir, self.prefix.odata))  # Prefix for output data

        return "\n".join(lines)

    @property
    def wfk_filepath(self):
        """Returns (at runtime) the absolute path of the WFK file produced by the NSCF run."""
        return self.nscf_node.outdir.has_abiext("WFK")

    @property
    def ddk_filepaths(self):
        """Returns (at runtime) the absolute path of the DDK files produced by the DDK runs."""
        # This to support new version of optic that used DDK.nc
        paths = [ddk_task.outdir.has_abiext("DDK.nc") for ddk_task in self.ddk_nodes]
        if all(p for p in paths):
            return paths

        # This is deprecated and can be removed when new version of Abinit is released.
        return [ddk_task.outdir.has_abiext("1WF") for ddk_task in self.ddk_nodes]

    def make_input(self):
        """Construct and write the input file of the calculation."""
        # Set the file paths.
        all_files ={"ddkfile_" + str(n + 1): ddk for n, ddk in enumerate(self.ddk_filepaths)}
        all_files.update({"wfkfile": self.wfk_filepath})
        files_nml = {"FILES": all_files}
        files= nmltostring(files_nml)

        # Get the input specified by the user
        user_file = nmltostring(self.input.as_dict())

        # Join them.
        return files + user_file

    def setup(self):
        """Public method called before submitting the task."""

    def make_links(self):
        """
        Optic allows the user to specify the paths of the input file.
        hence we don't need to create symbolic links.
        """

    def get_results(self, **kwargs):
        return super(OpticTask, self).get_results(**kwargs)

    def fix_abicritical(self):
        """
        Cannot fix abicritical errors for optic
        """
        return 0

    #@check_spectator
    def reset_from_scratch(self):
        """
        restart from scratch, this is to be used if a job is restarted with more resources after a crash
        """
        # Move output files produced in workdir to _reset otherwise check_status continues
        # to see the task as crashed even if the job did not run
        # Create reset directory if not already done.
        reset_dir = os.path.join(self.workdir, "_reset")
        reset_file = os.path.join(reset_dir, "_counter")
        if not os.path.exists(reset_dir):
            os.mkdir(reset_dir)
            num_reset = 1
        else:
            with open(reset_file, "rt") as fh:
                num_reset = 1 + int(fh.read())

        # Move files to reset and append digit with reset index.
        def move_file(f):
            if not f.exists: return
            try:
                f.move(os.path.join(reset_dir, f.basename + "_" + str(num_reset)))
            except OSError as exc:
                logger.warning("Couldn't move file {}. exc: {}".format(f, str(exc)))

        for fname in ("output_file", "log_file", "stderr_file", "qout_file", "qerr_file", "mpiabort_file"):
            move_file(getattr(self, fname))

        with open(reset_file, "wt") as fh:
            fh.write(str(num_reset))

        self.start_lockfile.remove()

        # Reset datetimes
        self.datetimes.reset()

        return self._restart(submit=False)

    def fix_queue_critical(self):
        """
        This function tries to fix critical events originating from the queue submission system.

        General strategy, first try to increase resources in order to fix the problem,
        if this is not possible, call a task specific method to attempt to decrease the demands.

        Returns:
            1 if task has been fixed else 0.
        """
        from pymatgen.io.abinit.scheduler_error_parsers import NodeFailureError, MemoryCancelError, TimeCancelError
        #assert isinstance(self.manager, TaskManager)

        if not self.queue_errors:
            if self.mem_scales or self.load_scales:
                try:
                    self.manager.increase_resources()  # acts either on the policy or on the qadapter
                    self.reset_from_scratch()
                    return
                except ManagerIncreaseError:
                    self.set_status(self.S_ERROR, msg='unknown queue error, could not increase resources any further')
                    raise FixQueueCriticalError
            else:
                self.set_status(self.S_ERROR, msg='unknown queue error, no options left')
                raise FixQueueCriticalError

        else:
            for error in self.queue_errors:
                logger.info('fixing: %s' % str(error))

                if isinstance(error, NodeFailureError):
                    # if the problematic node is known, exclude it
                    if error.nodes is not None:
                        try:
                            self.manager.exclude_nodes(error.nodes)
                            self.reset_from_scratch()
                            self.set_status(self.S_READY, msg='excluding nodes')
                        except:
                            raise FixQueueCriticalError
                    else:
                        self.set_status(self.S_ERROR, msg='Node error but no node identified.')
                        raise FixQueueCriticalError

                elif isinstance(error, MemoryCancelError):
                    # ask the qadapter to provide more resources, i.e. more cpu's so more total memory if the code
                    # scales this should fix the memeory problem
                    # increase both max and min ncpu of the autoparalel and rerun autoparalel
                    if self.mem_scales:
                        try:
                            self.manager.increase_ncpus()
                            self.reset_from_scratch()
                            self.set_status(self.S_READY, msg='increased ncps to solve memory problem')
                            return
                        except ManagerIncreaseError:
                            logger.warning('increasing ncpus failed')

                    # if the max is reached, try to increase the memory per cpu:
                    try:
                        self.manager.increase_mem()
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='increased mem')
                        return
                    except ManagerIncreaseError:
                        logger.warning('increasing mem failed')

                    # if this failed ask the task to provide a method to reduce the memory demand
                    try:
                        self.reduce_memory_demand()
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='decreased mem demand')
                        return
                    except DecreaseDemandsError:
                        logger.warning('decreasing demands failed')

                    msg = ('Memory error detected but the memory could not be increased neither could the\n'
                           'memory demand be decreased. Unrecoverable error.')
                    self.set_status(self.S_ERROR, msg)
                    raise FixQueueCriticalError

                elif isinstance(error, TimeCancelError):
                    # ask the qadapter to provide more time
                    try:
                        self.manager.increase_time()
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='increased wall time')
                        return
                    except ManagerIncreaseError:
                        logger.warning('increasing the walltime failed')

                    # if this fails ask the qadapter to increase the number of cpus
                    if self.load_scales:
                        try:
                            self.manager.increase_ncpus()
                            self.reset_from_scratch()
                            self.set_status(self.S_READY, msg='increased number of cpus')
                            return
                        except ManagerIncreaseError:
                            logger.warning('increase ncpus to speed up the calculation to stay in the walltime failed')

                    # if this failed ask the task to provide a method to speed up the task
                    try:
                        self.speed_up()
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='task speedup')
                        return
                    except DecreaseDemandsError:
                        logger.warning('decreasing demands failed')

                    msg = ('Time cancel error detected but the time could not be increased neither could\n'
                           'the time demand be decreased by speedup of increasing the number of cpus.\n'
                           'Unrecoverable error.')
                    self.set_status(self.S_ERROR, msg)

                else:
                    msg = 'No solution provided for error %s. Unrecoverable error.' % error.name
                    self.set_status(self.S_ERROR, msg)

        return 0

    def autoparal_run(self):
        """
        Find an optimal set of parameters for the execution of the Optic task
        This method can change the submission parameters e.g. the number of CPUs for MPI and OpenMp.

        Returns 0 if success
        """
        policy = self.manager.policy

        if policy.autoparal == 0: # or policy.max_ncpus in [None, 1]:
            logger.info("Nothing to do in autoparal, returning (None, None)")
            return 0

        if policy.autoparal != 1:
            raise NotImplementedError("autoparal != 1")

        ############################################################################
        # Run ABINIT in sequential to get the possible configurations with max_ncpus
        ############################################################################

        # Set the variables for automatic parallelization
        # Will get all the possible configurations up to max_ncpus
        # Return immediately if max_ncpus == 1
        max_ncpus = self.manager.max_cores
        if max_ncpus == 1: return 0

        autoparal_vars = dict(autoparal=policy.autoparal, max_ncpus=max_ncpus)
        self.set_vars(autoparal_vars)

        # Run the job in a shell subprocess with mpi_procs = 1
        # we don't want to make a request to the queue manager for this simple job!
        # Return code is always != 0
        process = self.manager.to_shell_manager(mpi_procs=1).launch(self)
        self.history.pop()
        retcode = process.wait()
        # To avoid: ResourceWarning: unclosed file <_io.BufferedReader name=87> in py3k
        process.stderr.close()
        #process.stdout.close()

        # Remove the variables added for the automatic parallelization
        self.input.remove_vars(list(autoparal_vars.keys()))

        ##############################################################
        # Parse the autoparal configurations from the main output file
        ##############################################################
        parser = ParalHintsParser()
        try:
            pconfs = parser.parse(self.output_file.path)
        except parser.Error:
            # In principle Abinit should have written a complete log file
            # because we called .wait() but sometimes the Yaml doc is incomplete and
            # the parser raises. Let's wait 5 secs and then try again.
            time.sleep(5)
            try:
                pconfs = parser.parse(self.output_file.path)
            except parser.Error:
                logger.critical("Error while parsing Autoparal section:\n%s" % straceback())
                return 2

        ######################################################
        # Select the optimal configuration according to policy
        ######################################################
        #optconf = self.find_optconf(pconfs)
        # Select the partition on which we'll be running and set MPI/OMP cores.
        optconf = self.manager.select_qadapter(pconfs)

        ####################################################
        # Change the input file and/or the submission script
        ####################################################
        self.set_vars(optconf.vars)

        # Write autoparal configurations to JSON file.
        d = pconfs.as_dict()
        d["optimal_conf"] = optconf
        json_pretty_dump(d, os.path.join(self.workdir, "autoparal.json"))

        ##############
        # Finalization
        ##############
        # Reset the status, remove garbage files ...
        self.set_status(self.S_INIT, msg='finished auto paralell')

        # Remove the output file since Abinit likes to create new files
        # with extension .outA, .outB if the file already exists.
        os.remove(self.output_file.path)
        #os.remove(self.log_file.path)
        os.remove(self.stderr_file.path)

        return 0


class AnaddbTask(Task):
    """Task for Anaddb runs (post-processing of DFPT calculations)."""

    color_rgb = np.array((204, 102, 255)) / 255

    def __init__(self, anaddb_input, ddb_node,
                 gkk_node=None, md_node=None, ddk_node=None, workdir=None, manager=None):
        """
        Create an instance of :class:`AnaddbTask` from a string containing the input.

        Args:
            anaddb_input: string with the anaddb variables.
            ddb_node: The node that will produce the DDB file. Accept :class:`Task`, :class:`Work` or filepath.
            gkk_node: The node that will produce the GKK file (optional). Accept :class:`Task`, :class:`Work` or filepath.
            md_node: The node that will produce the MD file (optional). Accept `Task`, `Work` or filepath.
            gkk_node: The node that will produce the GKK file (optional). Accept `Task`, `Work` or filepath.
            workdir: Path to the working directory (optional).
            manager: :class:`TaskManager` object (optional).
        """
        # Keep a reference to the nodes.
        self.ddb_node = Node.as_node(ddb_node)
        deps = {self.ddb_node: "DDB"}

        self.gkk_node = Node.as_node(gkk_node)
        if self.gkk_node is not None:
            deps.update({self.gkk_node: "GKK"})

        # I never used it!
        self.md_node = Node.as_node(md_node)
        if self.md_node is not None:
            deps.update({self.md_node: "MD"})

        self.ddk_node = Node.as_node(ddk_node)
        if self.ddk_node is not None:
            deps.update({self.ddk_node: "DDK"})

        super(AnaddbTask, self).__init__(input=anaddb_input, workdir=workdir, manager=manager, deps=deps)

    @classmethod
    def temp_shell_task(cls, inp, ddb_node, mpi_procs=1,
                        gkk_node=None, md_node=None, ddk_node=None, workdir=None, manager=None):
        """
        Build a :class:`AnaddbTask` with a temporary workdir. The task is executed via
        the shell with 1 MPI proc. Mainly used for post-processing the DDB files.

        Args:
            mpi_procs: Number of MPI processes to use.
            anaddb_input: string with the anaddb variables.
            ddb_node: The node that will produce the DDB file. Accept :class:`Task`, :class:`Work` or filepath.

        See `AnaddbInit` for the meaning of the other arguments.
        """
        # Build a simple manager to run the job in a shell subprocess
        import tempfile
        workdir = tempfile.mkdtemp() if workdir is None else workdir
        if manager is None: manager = TaskManager.from_user_config()

        # Construct the task and run it
        return cls(inp, ddb_node,
                   gkk_node=gkk_node, md_node=md_node, ddk_node=ddk_node,
                   workdir=workdir, manager=manager.to_shell_manager(mpi_procs=mpi_procs))

    @property
    def executable(self):
        """Path to the executable required for running the :class:`AnaddbTask`."""
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
        app(self.outdir.path_join("out"))  # 6) Base name for elphon output files e.g. t13
        app(self.ddk_filepath)             # 7) File containing ddk filenames for elphon/transport.

        return "\n".join(lines)

    @property
    def ddb_filepath(self):
        """Returns (at runtime) the absolute path of the input DDB file."""
        # This is not very elegant! A possible approach could to be path self.ddb_node.outdir!
        if isinstance(self.ddb_node, FileNode): return self.ddb_node.filepath
        path = self.ddb_node.outdir.has_abiext("DDB")
        return path if path else "DDB_FILE_DOES_NOT_EXIST"

    @property
    def md_filepath(self):
        """Returns (at runtime) the absolute path of the input MD file."""
        if self.md_node is None: return "MD_FILE_DOES_NOT_EXIST"
        if isinstance(self.md_node, FileNode): return self.md_node.filepath

        path = self.md_node.outdir.has_abiext("MD")
        return path if path else "MD_FILE_DOES_NOT_EXIST"

    @property
    def gkk_filepath(self):
        """Returns (at runtime) the absolute path of the input GKK file."""
        if self.gkk_node is None: return "GKK_FILE_DOES_NOT_EXIST"
        if isinstance(self.gkk_node, FileNode): return self.gkk_node.filepath

        path = self.gkk_node.outdir.has_abiext("GKK")
        return path if path else "GKK_FILE_DOES_NOT_EXIST"

    @property
    def ddk_filepath(self):
        """Returns (at runtime) the absolute path of the input DKK file."""
        if self.ddk_node is None: return "DDK_FILE_DOES_NOT_EXIST"
        if isinstance(self.ddk_node, FileNode): return self.ddk_node.filepath

        path = self.ddk_node.outdir.has_abiext("DDK")
        return path if path else "DDK_FILE_DOES_NOT_EXIST"

    def setup(self):
        """Public method called before submitting the task."""

    def make_links(self):
        """
        Anaddb allows the user to specify the paths of the input file.
        hence we don't need to create symbolic links.
        """

    def open_phbst(self):
        """Open PHBST file produced by Anaddb and returns :class:`PhbstFile` object."""
        from abipy.dfpt.phonons import PhbstFile
        phbst_path = os.path.join(self.workdir, "run.abo_PHBST.nc")
        if not phbst_path:
            if self.status == self.S_OK:
                logger.critical("%s reached S_OK but didn't produce a PHBST file in %s" % (self, self.outdir))
            return None

        try:
            return PhbstFile(phbst_path)
        except Exception as exc:
            logger.critical("Exception while reading GSR file at %s:\n%s" % (phbst_path, str(exc)))
            return None

    def open_phdos(self):
        """Open PHDOS file produced by Anaddb and returns :class:`PhdosFile` object."""
        from abipy.dfpt.phonons import PhdosFile
        phdos_path = os.path.join(self.workdir, "run.abo_PHDOS.nc")
        if not phdos_path:
            if self.status == self.S_OK:
                logger.critical("%s reached S_OK but didn't produce a PHBST file in %s" % (self, self.outdir))
            return None

        try:
            return PhdosFile(phdos_path)
        except Exception as exc:
            logger.critical("Exception while reading GSR file at %s:\n%s" % (phdos_path, str(exc)))
            return None

    def get_results(self, **kwargs):
        results = super(AnaddbTask, self).get_results(**kwargs)
        return results