# coding: utf-8
"""This module provides functions and classes related to Task objects."""
from __future__ import division, print_function, unicode_literals

import os
import time
import datetime
import shutil
import collections
import abc
import copy
import yaml
import six
import numpy as np

from pprint import pprint
from six.moves import map, zip, StringIO
from pydispatch import dispatcher
from monty.string import is_string, list_strings
from monty.collections import AttrDict
from monty.functools import lazy_property, return_none_if_raise
from monty.json import MontyDecoder
from monty.fnmatch import WildCard
from pymatgen.core.units import Memory
from pymatgen.serializers.json_coders import json_pretty_dump, pmg_serialize
from .utils import File, Directory, irdvars_for_ext, abi_splitext, abi_extensions, FilepathFixer, Condition, SparseHistogram
from .strategies import StrategyWithInput, OpticInput
from .qadapters import make_qadapter, QueueAdapter, slurm_parse_timestr
from .db import DBConnector
from .nodes import Status, Node, NodeResults
from . import abiinspect
from . import events


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
    "DdkTask",
    "PhononTask",
    "SigmaTask",
    "OpticTask",
    "AnaddbTask",
]

import logging
logger = logging.getLogger(__name__)


# Tools and helper functions.

def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


class TaskResults(NodeResults):

    JSON_SCHEMA = NodeResults.JSON_SCHEMA.copy() 
    JSON_SCHEMA["properties"] = {
        "executable": {"type": "string", "required": True},
    }

    @classmethod
    def from_node(cls, task):
        """Initialize an instance from an AbinitTask instance."""
        new = super(TaskResults, cls).from_node(task)

        new.update(
            executable=task.executable,
            #executable_version:
            #task_events=
            pseudos=task.strategy.pseudos.as_dict()
            #input=task.strategy
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

    # TODO: Change name in abinit
    # Remove tot_ncpus from Abinit
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


class ParalHintsParser(object):

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
                d = yaml.load(doc.text_notag)
                return ParalHints(info=d["info"], confs=d["configurations"])
            except:
                import traceback
                sexc = traceback.format_exc()
                err_msg = "Wrong YAML doc:\n%s\n\nException:\n%s" % (doc.text, sexc)
                self.add_error(err_msg)
                logger.critical(err_msg)
                raise self.Error(err_msg)


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

    def sort_by_speedup(self, reverse=True):
        """Sort the configurations in place. items with highest speedup come first"""
        self._confs.sort(key=lambda c: c.speedup, reverse=reverse)

    def sort_by_mem_per_proc(self, reverse=False):
        """Sort the configurations in place. items with lowest memory per proc come first."""
        # Avoid sorting if mem_per_cpu is not available.
        if any(c.mem_per_proc > 0.0 for c in self):
            self._confs.sort(key=lambda c: c.mem_per_proc, reverse=reverse)

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
            getattr(hints, "sort_by_" + policy.autoparal_priorities[0])()
        else:
            hints = hints.multidimensional_optimization(priorities=policy.autoparal_priorities)
            if len(hints) == 0: raise ValueError("len(hints) == 0")

        #TODO: make sure that num_cores == 1 is never selected when we have more than one configuration
        #if len(hints) > 1:
        #    hints.select_with_condition(dict(num_cores={"$eq": 1)))

        # Return final (orderded ) list of configurations (best first).
        return hints


class TaskPolicy(object):
    """
    This object stores the parameters used by the :class:`TaskManager` to 
    create the submission script and/or to modify the ABINIT variables 
    governing the parallel execution. A `TaskPolicy` object contains 
    a set of variables that specify the launcher, as well as the options
    and the condition used to select the optimal configuration for the parallel run 
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
            elif isinstance(obj, collections.Mapping):
                return cls(**obj) 
            else:
                raise TypeError("Don't know how to convert type %s to %s" % (type(obj), cls))

    @classmethod
    def autodoc(cls):
        return """
    autoparal: 0 to disable the autoparal feature (default 1 i.e. autoparal is on)
    condition: condition used to filter the autoparal configurations (Mongodb-like syntax). 
               Default: empty 
    vars_condition: condition used to filter the list of ABINIT variables reported autoparal 
                    (Mongodb-like syntax). Default: empty
    frozen_timeout: A job is considered to be frozen and its status is set to Error if no change to 
                    the output file has been done for frozen_timeout seconds. Accepts int with seconds or 
                    string in slurm form i.e. days-hours:minutes:seconds. Default: 1 hour.
    precedence:
    autoparal_priorities:
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
        self.frozen_timeout = slurm_parse_timestr(kwargs.pop("frozen_timeout", "0-1"))

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


class TaskManager(object):
    """
    A `TaskManager` is responsible for the generation of the job script and the submission 
    of the task, as well as for the specification of the parameters passed to the resource manager
    (e.g. Slurm, PBS ...) and/or the run-time specification of the ABINIT variables governing the parallel execution. 
    A `TaskManager` delegates the generation of the submission script and the submission of the task to the :class:`QueueAdapter`. 
    A `TaskManager` has a :class:`TaskPolicy` that governs the specification of the parameters for the parallel executions.
    Ideally, the TaskManager should be the **main entry point** used by the task to deal with job submission/optimization
    """
    YAML_FILE = "manager.yml"
    USER_CONFIG_DIR = os.path.join(os.getenv("HOME"), ".abinit", "abipy")

    ENTRIES = {"policy", "qadapters", "db_connector"}

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

##########################################
# Individual entries are documented below:
##########################################

"""
        s += "policy: " + TaskPolicy.autodoc() + "\n"
        s += "qadapter: " + QueueAdapter.autodoc() + "\n"
        s += "db_connector: " + DBConnector.autodoc()
        return s

    @classmethod
    def from_user_config(cls):
        """
        Initialize the :class:`TaskManager` from the YAML file 'taskmanager.yaml'.
        Search first in the working directory and then in the abipy configuration directory.

        Raises:
            RuntimeError if file is not found.
        """
        # Try in the current directory.
        path = os.path.join(os.getcwd(), cls.YAML_FILE)
        if os.path.exists(path):
            return cls.from_file(path)

        # Try in the configuration directory.
        path = os.path.join(cls.USER_CONFIG_DIR, cls.YAML_FILE)
        if os.path.exists(path):
            return cls.from_file(path)
    
        raise RuntimeError("Cannot locate %s neither in current directory nor in %s" % (cls.YAML_FILE, path))

    @classmethod
    def from_file(cls, filename):
        """Read the configuration parameters from the Yaml file filename."""
        try:
            with open(filename, "r") as fh:
                return cls.from_dict(yaml.load(fh))
        except Exception as exc:
            print("Error while reading TaskManager parameters from file %s\n" % filename)
            raise 

    @classmethod
    def from_string(cls, s):
        """Create an instance from string s containing a YAML dictionary."""
        return cls.from_dict(yaml.load(s))

    @classmethod
    def as_manager(cls, obj):
        """
        Convert obj into TaskManager instance. Accepts string, filepath, dictionary, TaskManager object.
        """
        if isinstance(obj, cls): return obj

        if is_string(obj):
            if os.path.exists(obj):
                return cls.from_file(obj)
            else:
                return cls.from_string(obj)

        elif isinstance(obj, collections.Mapping):
            return cls.from_dict(obj)
        else:
            raise TypeError("Don't know how to convert type %s to TaskManager" % type(obj))

    @classmethod
    def from_dict(cls, d):
        """Create an instance from a dictionary."""
        return cls(**{k: v for k, v in d.items() if k in cls.ENTRIES})

    def as_dict(self):
        return self._kwargs

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

        if kwargs:
            raise ValueError("Found invalid keywords in the taskmanager file:\n %s" % str(list(kwargs.keys())))

    def to_shell_manager(self, mpi_procs=1):
        """
        Returns a new `TaskManager` with the same parameters as self but replace the :class:`QueueAdapter`
        with a :class:`ShellAdapter` with mpi_procs so that we can submit the job without passing through the queue.
        """
        my_kwargs = copy.deepcopy(self._kwargs)
        my_kwargs["policy"] = TaskPolicy(autoparal=0)

        for d in my_kwargs["qadapters"]:
            #print("before", d["queue"]["qtype"])
            d["queue"]["qtype"] = "shell"
            d["limits"]["min_cores"] = mpi_procs
            d["limits"]["max_cores"] = mpi_procs

        #print(my_kwargs)
        new = self.__class__(**my_kwargs)
        new.set_mpi_procs(mpi_procs)

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
                for pconf in pconfs:
                    if qad.can_run_pconf(pconf):
                        self._use_qadpos_pconf(qadpos, pconf)
                        return pconf

        elif policy.precedence == "autoparal_conf":
            # Try to run on the first pconf irrespectively of the priority of the qadapter.
            for pconf in pconfs:
                for qadpos, qad in enumerate(self.qads):
                    if qad.can_run_pconf(pconf):
                        self._use_qadpos_pconf(qadpos, pconf)
                        return pconf

        else:
            raise ValueError("Wrong value of policy.precedence = %s" % policy.precedence)

        # No qadapter could be found
        raise RuntimeError("Cannot find qadapter for this run!")

    def _use_qadpos_pconf(self, qadpos, pconf):
        """This function is called when we have accepted the :class:`ParalConf` pconf""" 
        self._qadpos = qadpos 

        # Change the number of MPI/OMP cores.
        self.set_mpi_procs(pconf.mpi_procs)
        if self.has_omp: self.set_omp_threads(pconf.omp_threads)
                                                                      
        # Set memory per proc.
        self.set_mem_per_proc(pconf.mem_per_proc)

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
        """Set the number of MPI nodes to use."""
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
        return max(q.max_cores for q in self.qads)

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

    def write_jobfile(self, task):
        """Write the submission script. return the path of the script"""
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
        with open(task.job_file.path, "w") as fh:
            fh.write(script)
            os.chmod(0o740)
            return task.job_file.path

    def launch(self, task):
        """
        Build the input files and submit the task via the :class:`Qadapter` 

        Args:
            task: :class:`TaskObject`
        
        Returns:
            Process object.
        """
        if task.status == task.S_LOCKED:
            raise ValueError("You shall not try to submit a locked task!")

        # Build the task 
        task.build()

        # Write the submission script
        script_file = self.write_jobfile(task)
        task.set_status(task.S_SUB)

        # Submit the task and save the queue id.
        try:
            qjob, process = self.qadapter.submit_to_queue(script_file)
            task.set_qjob(qjob)
            return process
        except self.qadapter.MaxNumLaunchesError:
            # TODO: Here we should try to switch to another qadapter
            # 1) Find a new parallel configuration in those stores in task.pconfs
            # 2) Change the input file.
            # 3) Regenerate the submission script
            # 4) Relaunch
            raise

    def get_collection(self, **kwargs):
        """Return the MongoDB collection used to store the results."""
        return self.db_connector.get_collection(**kwargs)

    def increase_resources(self):
        # with GW calculations in mind with GW mem = 10, 
        # the response fuction is in memory and not distributed
        # we need to increas memory if jobs fail ...
        return self.qadapter.more_mem_per_proc()


class FakeProcess(object):
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


class TaskDateTimes(object):
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


class TaskError(Exception):
    """Base Exception for :class:`Task` methods"""
        

class TaskRestartError(TaskError):
    """Exception raised while trying to restart the :class:`Task`."""


class Task(six.with_metaclass(abc.ABCMeta, Node)):
    """A Task is a node that performs some kind of calculation."""
    # Use class attributes for TaskErrors so that we don't have to import them.
    Error = TaskError
    RestartError = TaskRestartError

    # List of `AbinitEvent` subclasses that are tested in the check_status method. 
    # Subclasses should provide their own list if they need to check the converge status.
    CRITICAL_EVENTS = [
    ]

    # List of task specific default event handlers
    # Subclasses should provide their own list
    #EVENT_HANDLERS = []

    # Prefixes for Abinit (input, output, temporary) files.
    Prefix = collections.namedtuple("Prefix", "idata odata tdata")
    pj = os.path.join

    prefix = Prefix(pj("indata", "in"), pj("outdata", "out"), pj("tmpdata", "tmp"))
    del Prefix, pj

    def __init__(self, strategy, workdir=None, manager=None, deps=None):
        """
        Args:
            strategy: Input file or :class:`Strategy` instance defining the calculation.
            workdir: Path to the working directory.
            manager: :class:`TaskManager` object.
            deps: Dictionary specifying the dependency of this node.
                  None means that this Task has no dependency.
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
            self.add_deps(deps)

        # Date-time associated to submission, start and end.
        self.datetimes = TaskDateTimes()

        # Count the number of restarts.
        self.num_restarts = 0

        self._qjob = None
        self.queue_errors = []
        self.abi_errors = []

    def __getstate__(self):
        """
        Return state is pickled as the contents for the instance.
                                                                                      
        In this case we just remove the process since Subprocess objects cannot be pickled.
        This is the reason why we have to store the returncode in self._returncode instead
        of using self.process.returncode.
        """
        return {k: v for k, v in self.__dict__.items() if k not in ["_process"]}

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

        # stderr and output file of the queue manager. Note extensions.
        self.qerr_file = File(os.path.join(self.workdir, "queue.qerr"))
        self.qout_file = File(os.path.join(self.workdir, "queue.qout"))

    def set_manager(self, manager):
        """Set the :class:`TaskManager` to use to launch the Task."""
        self.manager = manager.deepcopy()

        # TODO
        # Select adapters associated to the Task class
        #keep = []
        #for i, qad in enumerate(self.manager.qads):
        #    if self.__class__.__name__ in qad.task_classes:
        #        keep.append(i)
        #if keep:
        #    self._qads = [self._qads[i] for i in keep]
        #    self._qadpos = 0

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

    def get_inpvar(self, varname):
        """Return the value of the ABINIT variable varname, None if not present.""" 
        if hasattr(self.strategy, "abinit_input"):
            try:
                return self.strategy.abinit_input[0][varname]
            except KeyError:
                return None
        else:
            raise NotImplementedError("get_var for HTC interface!")

    def _set_inpvar(self, varname, value):
        """Set the value of the ABINIT variable varname. Return old value."""
        old_value = self.get_inpvar(varname)
        self.strategy.abinit_input[0][varname] = value
        return old_value

    @property
    def initial_structure(self):
        """Initial structure of the task."""
        if hasattr(self.strategy, "abinit_input"):
            return self.strategy.abinit_input.structure
            #return self.strategy.abinit_input[0].structure
        else:
            return self.strategy.structure

    def make_input(self, with_header=False):
        """Construct the input file of the calculation."""
        s = self.strategy.make_input()
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

    @property
    def is_completed(self):
        """True if the task has been executed."""
        return self.status >= self.S_DONE

    @property
    def can_run(self):
        """The task can run if its status is < S_SUB and all the other dependencies (if any) are done!"""
        all_ok = all([stat == self.S_OK for stat in self.deps_status])
        return self.status < self.S_SUB and self.status != self.S_LOCKED and all_ok

    def cancel(self):
        """Cancel the job. Returns 1 if job was cancelled."""
        if self.queue_id is None: return 0 
        if self.status >= self.S_DONE: return 0 

        exit_status = self.manager.cancel(self.queue_id)
        if exit_status != 0: return 0

        # Remove output files and reset the status.
        self.reset()
        return 1

    def _on_done(self):
        self.fix_ofiles()

    def _on_ok(self):
        # Fix output file names.
        self.fix_ofiles()

        # Get results
        results = self.on_ok()

        # Set internal flag.
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
        return dict(returncode=0, message="Calling on_all_ok of the base class!")

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

        if submit:
            # Remove the lock file
            self.start_lockfile.remove()
            # Relaunch the task.
            fired = self.start()
            if not fired: self.history.warning("Restart failed")
        else:
            fired = False

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
        self.set_status(self.S_ERROR)
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

        # Remove output files otherwise the EventParser will think the job is still running
        self.output_file.remove()
        self.log_file.remove()
        self.stderr_file.remove()
        self.start_lockfile.remove()
        self.qerr_file.remove()
        self.qout_file.remove()

        self.set_status(self.S_INIT, msg="Reset on %s" % time.asctime())
        self.set_qjob(None)

        # TODO send a signal to the flow 
        #self.work.check_status()
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
        """
        Lock the file, source is the :class:`Node` that is locking the file.
        """
        if self.status != self.S_INIT:
            raise ValueError("Trying to lock a task with status %s" % self.status)

        self._status = self.S_LOCKED
        self.history.info("Locked by node %s", source_node)

    def unlock(self, check_status=True):
        """
        Unlock the task, set its status to `S_READY` so that the scheduler can submit it.
        Call task.check_status if check_status is True.
        """
        if self.status != self.S_LOCKED:
            raise RuntimeError("Trying to unlock a task with status %s" % self.status)

        self._status = self.S_READY
        if check_status: self.check_status()

    def set_status(self, status, msg=None):
        """
        Set and return the status of the task.

        Args:
            status: Status object or string representation of the status
            msg: string with human-readable message used in the case of errors (optional)
        """
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
                self.history.info("Submitted with MPI=%s, Omp=%s, Memproc=%.1f [Gb]" % (
                    self.mpi_procs, self.omp_threads, self.mem_per_proc.to("Gb")))

            if status == self.S_OK:
                self.history.info("Task Completed")

            if status == self.S_ABICRITICAL:
                self.history.info("Status set to S_ABI_CRITICAL due to: %s", msg)

        if status == self.S_DONE:
            # Execute the callback
            self._on_done()
                                                                                
        if status == self.S_OK:
            # Finalize the task.
            if not self.finalized:
                self._on_ok()
                # here we remove the output files of the task and of its parents.
                if self.cleanup_exts: self.clean_output_files()
                                                                                
            logger.debug("Task %s broadcasts signal S_OK" % self)
            dispatcher.send(signal=self.S_OK, sender=self)

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
        if self.status in black_list: return self.status

        # 2) Check the returncode of the process (the process of submitting the job) first.
        # this point type of problem should also be handled by the scheduler error parser
        if self.returncode != 0:
            # The job was not submitted properly
            return self.set_status(self.S_QCRITICAL, msg="return code %s" % self.returncode)

        # Analyze the stderr file for Fortran runtime errors.
        err_msg = None
        if self.stderr_file.exists:
            err_msg = self.stderr_file.read()

        # Analyze the stderr file of the resource manager runtime errors.
        err_info = None
        if self.qerr_file.exists:
            err_info = self.qerr_file.read()

        # Start to check ABINIT status if the output file has been created.
        if self.output_file.exists:
            try:
                report = self.get_event_report()
            except Exception as exc:
                msg = "%s exception while parsing event_report:\n%s" % (self, exc)
                logger.critical(msg)
                return self.set_status(self.S_ABICRITICAL, msg=msg)

            if report.run_completed:
                # Here we  set the correct timing data reported by Abinit
                self.datetimes.start = report.start_datetime
                self.datetimes.end = report.end_datetime

                # Check if the calculation converged.
                not_ok = report.filter_types(self.CRITICAL_EVENTS)
                if not_ok:
                    return self.set_status(self.S_UNCONVERGED)
                else:
                    return self.set_status(self.S_OK)

            # Calculation still running or errors?
            if report.errors or report.bugs:
                # Abinit reported problems
                if report.errors:
                    logger.debug('Found errors in report')
                    for error in report.errors:
                        logger.debug(str(error))
                        try:
                            self.abi_errors.append(error)
                        except AttributeError:
                            self.abi_errors = [error]
                if report.bugs:
                    logger.debug('Found bugs in report:')
                    for bug in report.bugs:
                        logger.debug(str(bug))

                # The job is unfixable due to ABINIT errors
                logger.debug("%s: Found Errors or Bugs in ABINIT main output!" % self)
                msg = "\n".join(map(repr, report.errors + report.bugs))
                return self.set_status(self.S_ABICRITICAL, msg=msg)

            # 5)
            if self.stderr_file.exists and not err_info:
                if self.qerr_file.exists and not err_msg:
                    # there is output and no errors
                    # The job still seems to be running
                    return self.set_status(self.S_RUN)

        # 6)
        if not self.output_file.exists:
            logger.debug("output_file does not exists")
            if not self.stderr_file.exists and not self.qerr_file.exists:     
                # No output at allThe job is still in the queue.
                return self.status
                
        # 7) Analyze the files of the resource manager and abinit and execution err (mvs)
        if err_info:
            from pymatgen.io.abinitio.scheduler_error_parsers import get_parser
            scheduler_parser = get_parser(self.manager.qadapter.QTYPE, err_file=self.qerr_file.path,
                                          out_file=self.qout_file.path, run_err_file=self.stderr_file.path)

            if scheduler_parser is None:
                return self.set_status(self.S_QCRITICAL, 
                                       msg="Cannot find scheduler_parser for qtype %s" % self.manager.qadapter.QTYPE)
                
            scheduler_parser.parse()

            if scheduler_parser.errors:
                self.queue_errors = scheduler_parser.errors
                # the queue errors in the task
                msg = "scheduler errors found:\n%s" % str(scheduler_parser.errors)
                logger.critical(msg)
                return self.set_status(self.S_QCRITICAL, msg=msg)
                # The job is killed or crashed and we know what happened
            else:
                if len(err_info) > 0:
                    logger.debug('found unknown queue error: %s' % str(err_info))
                    return self.set_status(self.S_QCRITICAL, msg=err_info)
                    # The job is killed or crashed but we don't know what happened
                    # it is set to QCritical, we will attempt to fix it by running on more resources

        # 8) analizing the err files and abinit output did not identify a problem
        # but if the files are not empty we do have a problem but no way of solving it:
        if err_msg is not None and len(err_msg) > 0:
            logger.debug('found error message:\n %s' % str(err_msg))
            return self.set_status(self.S_QCRITICAL, msg=err_info)
            # The job is killed or crashed but we don't know what happend
            # it is set to QCritical, we will attempt to fix it by running on more resources

        # 9) if we still haven't returned there is no indication of any error and the job can only still be running
        # but we should actually never land here, or we have delays in the file system ....
        # print('the job still seems to be running maybe it is hanging without producing output... ')

        # Check time of last modification.
        if self.output_file.exists and \
           (time.time() - self.output_file.get_stat().st_mtime > self.manager.policy.frozen_timeout):
            msg = "Task seems to be frozen, last change more than %s [s] ago" % self.manager.policy.frozen_timeout
            return self.set_status(self.S_ERROR, msg)

        return self.set_status(self.S_RUN)

    def reduce_memory_demand(self):
        """
        Method that can be called by the flow to decrease the memory demand of a specific task.
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
        logger.debug("Linking path %s --> %s" % (filepath, infile))

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
                    # Try netcdf file. TODO: this case should be treated in a cleaner way.
                    path += "-etsf.nc"
                    if os.path.exists(path): dest += "-etsf.nc"

                if not os.path.exists(path):
                    raise self.Error("%s: %s is needed by this task but it does not exist" % (self, path))

                # Link path to dest if dest link does not exist.
                # else check that it points to the expected file.
                logger.debug("Linking path %s --> %s" % (path, dest))

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

    def get_event_report(self, source="log"):
        """
        Analyzes the main output file for possible Errors or Warnings.

        Args:
            source: "output" for the main output file,"log" for the log file.

        Returns:
            :class:`EventReport` instance or None if the main output file does not exist.
        """
        # TODO: For the time being, we inspect the log file,
        # We will start to use the output file when the migration to YAML is completed
        ofile = {
            "output": self.output_file,
            "log": self.log_file}[source]

        if not ofile.exists:
            return None

        # Don't parse source file if we already have its report and the source didn't change.
        #if not hasattr(self, "_prev_reports"): self._prev_reports = {}
        #old_report = self._prev_reports.get(source, None)
        #if False and old_report is not None and old_report.stat.st_mtime == ofile.get_stat().st_mtime:
        #    print("Returning old_report")
        #    return old_report

        parser = events.EventsParser()
        try:
            report = parser.parse(ofile.path)
            #self._prev_reports[source] = report
            return report

        except parser.Error as exc:
            # Return a report with an error entry with info on the exception.
            logger.critical("%s: Exception while parsing ABINIT events:\n %s" % (ofile, str(exc)))
            self.set_status(self.S_ABICRITICAL, msg=str(exc))
            return parser.report_exception(ofile.path, exc)

    def get_results(self, **kwargs):
        """
        Returns :class:`NodeResults` instance.
        Subclasses should extend this method (if needed) by adding 
        specialized code that performs some kind of post-processing.
        """
        # Check whether the process completed.
        if self.returncode is None:
            raise self.Error("return code is None, you should call wait, communitate or poll")

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
        exts = self.cleanup_exts.difference(except_exts)
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
            exts = self.cleanup_exts.difference(except_exts)
            #print("%s removes extensions %s from parent node %s" % (self, exts, parent))
            paths += parent.outdir.remove_exts(exts)

        #print("%s:\n Files removed: %s" % (self, paths))
        return paths

    def setup(self):
        """Base class does not provide any hook."""

    def start(self, **kwargs):
        """
        Starts the calculation by performing the following steps:

            - build dirs and files
            - call the _setup method
            - execute the job file by executing/submitting the job script.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        autoparal       False to skip the autoparal step (default True)
        ==============  ==============================================================

        Returns:
            1 if task was started, 0 otherwise.
            
        """
        if self.status >= self.S_SUB:
            raise self.Error("Task status: %s" % str(self.status))

        if self.start_lockfile.exists:
            logger.warning("Found lock file: %s" % self.start_lockfile.path)
            return 0

        self.start_lockfile.write("Started on %s" % time.asctime())

        self.build()
        self._setup()

        # Add the variables needed to connect the node.
        for d in self.deps:
            cvars = d.connecting_vars()
            logger.debug("Adding connecting vars %s " % cvars)
            self.strategy.add_extra_abivars(cvars)

            # Get (python) data from other nodes
            d.apply_getters(self)

        # Automatic parallelization
        if kwargs.pop("autoparal", True) and hasattr(self, "autoparal_run"):
            try:
                self.autoparal_run()
            except Exception as exc:
                # Sometimes autparal_run fails because Abinit aborts
                # at the level of the parser e.g. cannot find the spacegroup
                # due to some numerical noise in the structure.
                # In this case we call fix_abi_critical and then we try to run autoparal again.
                self.history.critical("First call to autoparal failed with `%s`. Will try fix_abi_critical" % exc)
                msg = "autoparal_fake_run raised:\n%s" % straceback()
                logger.critical(msg)

                fixed = self.fix_abi_critical()
                if not fixed:
                    self.set_status(self.S_ABICRITICAL, msg="fix_abi_critical could not solve the problem")
                    return 0

                try:
                    self.autoparal_run()
                    self.history.info("Second call to autoparal succeeded!")

                except Exception as exc:
                    self.history.critical("Second call to autoparal failed with %s. Cannot recover!", exc)
                    msg = "Tried autoparal again but got:\n%s" % straceback()
                    logger.critical(msg)

                    self.set_status(self.S_ABICRITICAL, msg=msg)
                    return 0

        # Start the calculation in a subprocess and return.
        self._process = self.manager.launch(self)
        return 1

    def start_and_wait(self, *args, **kwargs):
        """
        Helper method to start the task and wait for completetion.

        Mainly used when we are submitting the task via the shell without passing through a queue manager.
        """
        self.start(*args, **kwargs)
        retcode = self.wait()
        return retcode

    @pmg_serialize
    def as_dict(self):
        d = {}
        d['strategy'] = self.strategy.as_dict()
        d['workdir'] = getattr(self, 'workdir', None)
        if hasattr(self, 'manager'):
            d['manager'] = self.manager.as_dict()
        else:
            d['manager'] = None
        return d

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        strategy = dec.process_decoded(d['strategy'])
        manager = TaskManager.from_dict(d['manager'])
        return cls(strategy=strategy, workdir=d['workdir'], manager=manager)


class AbinitTask(Task):
    """
    Base class defining an ABINIT calculation
    """
    Results = TaskResults

    @classmethod
    def from_input(cls, ainput, workdir=None, manager=None):
        """
        Create an instance of `AbinitTask` from an ABINIT input.
    
        Args:
            ainput: `AbinitInput` object.
            workdir: Path to the working directory.
            manager: :class:`TaskManager` object.
        """
        # TODO: Find a better way to do this. I will likely need to refactor the Strategy object
        strategy = StrategyWithInput(ainput, deepcopy=True)

        return cls(strategy, workdir=workdir, manager=manager)

    def setup(self):
        """
        Abinit has the very *bad* habit of changing the file extension by appending the characters in [A,B ..., Z] 
        to the output file, and this breaks a lot of code that relies of the use of a unique file extension.
        Here we fix this issue by renaming run.abo to run.abo_[number] if the output file "run.abo" already
        exists. A few lines of code in python, a lot of problems if you try to implement this trick in Fortran90. 
        """
        if self.output_file.exists:
            # Find the index of the last file (if any) and push.
            # TODO: Maybe it's better to use run.abo --> run(1).abo
            fnames = [f for f in os.listdir(self.workdir) if f.startswith(self.output_file.basename)]
            nums = [int(f) for f in [f.split("_")[-1] for f in fnames] if f.isdigit()]
            last = max(nums) if nums else 0
            new_path = self.output_file.path + "_" + str(last+1)

            logger.info("Will rename %s to %s" % (self.output_file.path, new_path))
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
        znucl = self.strategy.structure.to_abivars()["znucl"]

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
        logger.info("in autoparal_run")
        policy = self.manager.policy

        if policy.autoparal == 0: # or policy.max_ncpus in [None, 1]:
            logger.info("Nothing to do in autoparal, returning (None, None)")
            return 1

        if policy.autoparal != 1:
            raise NotImplementedError("autoparal != 1")

        ############################################################################
        # Run ABINIT in sequential to get the possible configurations with max_ncpus
        ############################################################################

        # Set the variables for automatic parallelization
        # Will get all the possible configurations up to max_ncpus
        max_ncpus = self.manager.max_cores
        autoparal_vars = dict(autoparal=policy.autoparal, max_ncpus=max_ncpus)
        self.strategy.add_extra_abivars(autoparal_vars)

        # Run the job in a shell subprocess with mpi_procs = 1
        # we don't want to make a request to the queue manager for this simple job!
        # Return code is always != 0 
        process = self.manager.to_shell_manager(mpi_procs=1).launch(self)
        logger.info("fake run launched")
        self.history.pop()
        retcode = process.wait()

        # Remove the variables added for the automatic parallelization
        self.strategy.remove_extra_abivars(autoparal_vars.keys())

        ##############################################################
        # Parse the autoparal configurations from the main output file
        ##############################################################
        parser = ParalHintsParser()
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
        self.strategy.add_extra_abivars(optconf.vars)

        # Write autoparal configurations to JSON file.
        d = pconfs.as_dict()
        d["optimal_conf"] = optconf
        json_pretty_dump(d, os.path.join(self.workdir, "autoparal.json"))

        ##############
        # Finalization
        ##############
        # Reset the status, remove garbage files ...
        self.set_status(self.S_INIT)

        # Remove the output file since Abinit likes to create new files 
        # with extension .outA, .outB if the file already exists.
        os.remove(self.output_file.path)
        os.remove(self.log_file.path)
        os.remove(self.stderr_file.path)

        return 0

    def find_optconf(self, pconfs):
        # Save pconfs for future reference.
        self.set_pconfs(pconfs)
                                                                                
        # Select the partition on which we'll be running and set MPI/OMP cores.
        optconf = self.manager.select_qadapter(pconfs)
        return optconf

    def get_ibz(self, ngkpt=None, shiftk=None):
        """
        This function calls ABINIT to compute the list of k-points in the IBZ.
        By default, we use the value of ngkpt and shiftk in the strategy but it's also
        possible to change temporary these values.

        Returns:
            kpoints: numpy array with the reduced coordinates of the k-points in the IBZ.
            weights: numpy array with the weights normalized to one.
        """
        logger.info("in get_ibz")

        #########################################
        # Run ABINIT in sequential to get the IBZ
        #########################################

        # Set the variables for automatic parallelization
        ibz_vars = dict(prtkpt=-2)
        if ngkpt is not None: ibz_vars["ngkpt"] = ngkpt
        if shiftk is not None:
            shiftk = np.resphape(shiftk, (-1,3))
            ibz_vars["shiftk"] = shiftk
            ibz_vars["nshiftk"] = len(shiftk)

        self.strategy.add_extra_abivars(ibz_vars)
        # Build a simple manager to run the job in a shell subprocess
        # we don't want to make a request to the queue manager for this simple job!
        seq_manager = self.manager.to_shell_manager(mpi_procs=1)

        # Return code is always != 0
        process = seq_manager.launch(self)
        retcode = process.wait()

        # Remove the variables added for the automatic parallelization
        self.strategy.remove_extra_abivars(ibz_vars.keys())

        ################################################
        # Read the list of k-points from the netcdf file
        ################################################
        from pymatgen.io.abinitio import NetcdfReader
        with NetcdfReader(self.outdir.path_in("kpts.nc")) as r:
            kpoints = r.read_value("reduced_coordinates_of_kpoints")
            weights = r.read_value("kpoint_weights")

        self.set_status(self.S_INIT)

        # Remove the output file since Abinit likes to create new files
        # with extension .outA, .outB if the file already exists.
        os.remove(self.output_file.path)
        os.remove(self.log_file.path)
        os.remove(self.stderr_file.path)

        return kpoints, weights

    def restart(self):
        """
        general restart used when scheduler problems have been taken care of
        """
        return self._restart()

    def reset_from_scratch(self):
        """
        restart from scratch, reuse of output
        this is to be used if a job is restarted with more resources after a crash
        """
        # remove all 'error', else the job will be seen as crashed in the next check status
        # even if the job did not run
        self.output_file.remove()
        self.log_file.remove()
        self.stderr_file.remove()
        self.start_lockfile.remove()

        return self._restart(submit=False)

    def fix_abi_critical(self):
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

        # Note we have loop over all possible events (slow, I know)
        # because we can have handlers for Error, Bug or Warning 
        # (ideally only for CriticalWarnings but this is not done yet) 
        for event in report:
            for i, handler in enumerate(self.event_handlers):

                if handler.can_handle(event) and not done[i]:
                    logger.info("handler", handler, "will try to fix", event)
                    try:
                        d = handler.handle(self, event)
                        if d: 
                            done[i] += 1
                            count += 1

                    except Exception as exc:
                        logger.critical(str(exc))

        if count:
            self.reset_from_scratch()
            return 1

        msg = 'We encountered AbiCritical events that could not be fixed'
        logger.critical(msg)
        self.set_status(status=self.S_ERROR, msg=msg)
        return 0

    def fix_queue_critical(self):
        """
        This function tries to fix critical events originating from the queue submission system.

        General strategy, first try to increase resources in order to fix the problem,
        if this is not possible, call a task specific method to attempt to decrease the demands.

        Returns:
            1 if task has been fixed else 0.
        """
        from pymatgen.io.abinitio.scheduler_error_parsers import NodeFailureError, MemoryCancelError, TimeCancelError

        if not self.queue_errors:
            # queue error but no errors detected, try to solve by increasing resources
            # if resources are at maximum the task is definitively turned to errored
            if self.manager.increase_resources():  # acts either on the policy or on the qadapter
                self.reset_from_scratch()
                return 1 
            else:
                self.set_status(self.S_ERROR, msg='unknown queue error, could not increase resources any further')

        else:
            for error in self.queue_errors:
                logger.info('fixing: %s' % str(error))

                if isinstance(error, NodeFailureError):
                    # if the problematic node is known, exclude it
                    if error.nodes is not None:
                        self.manager.qadapter.exclude_nodes(error.nodes)
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='increased resources')
                        return 1
                    else:
                        self.set_status(self.S_ERROR, msg='Node error but no node identified.')

                elif isinstance(error, MemoryCancelError):
                    # ask the qadapter to provide more resources, i.e. more cpu's so more total memory
                    if self.manager.increase_resources():
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='increased mem')
                        return 1

                    # if the max is reached, try to increase the memory per cpu:
                    elif self.manager.qadapter.increase_mem():
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='increased mem')
                        return 1

                    # if this failed ask the task to provide a method to reduce the memory demand
                    elif self.reduce_memory_demand():
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='decreased mem demand')
                        return 1

                    else:
                        msg = ('Memory error detected but the memory could not be increased neigther could the\n'
                              'memory demand be decreased. Unrecoverable error.')
                        return self.set_status(self.S_ERROR, msg)

                elif isinstance(error, TimeCancelError):
                    # ask the qadapter to provide more memory
                    if self.manager.qadapter.increase_time():
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='increased wall time')
                        return 1

                    # if this fails ask the qadapter to increase the number of cpus
                    elif self.manager.increase_resources():
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='increased number of cpus')
                        return 1

                    # if this failed ask the task to provide a method to speed up the task
                    # MG TODO: Remove this
                    elif self.speed_up():
                        self.reset_from_scratch()
                        self.set_status(self.S_READY, msg='task speedup')
                        return 1

                    else:
                        msg = ('Time cancel error detected but the time could not be increased neither could\n'
                               'the time demand be decreased by speedup of increasing the number of cpus.\n'
                               'Unrecoverable error.')
                        self.set_status(self.S_ERROR, msg)
                else:
                    msg = 'No solution provided for error %s. Unrecoverable error.' % error.name
                    self.set_status(self.S_ERROR, msg)

        return 0

# TODO
# Enable restarting capabilites:
# Before doing so I need:
#   1) Preliminary standardization of the ABINT events and critical WARNINGS (YAML)
#   2) Change the parser so that we can use strings in the input file.
#      We need this change for restarting structural relaxations so that we can read 
#      the initial structure from file.


class ProduceGsr(object):
    """
    Mixin class for an :class:`AbinitTask` producing a GSR file.
    Provide the method `open_gsr` that reads and return a GSR file.
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


class ProduceHist(object):
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


class ProduceDdb(object):
    """
    Mixin class for :an class:`AbinitTask` producing a DDB file.
    Provide the method `open_ddb` that reads and return a Ddb file.
    """
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

        # Open the GSR file.
        from abipy.dfpt.ddb import DdbFile
        try:
            return DdbFile(ddb_path)
        except Exception as exc:
            logger.critical("Exception while reading DDB file at %s:\n%s" % (ddb_path, str(exc)))
            return None


class ScfTask(AbinitTask, ProduceGsr):
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
        for ext in ("WFK", "DEN"):
            restart_file = self.outdir.has_abiext(ext)
            irdvars = irdvars_for_ext(ext)
            if restart_file:
                break
        else:
            raise self.RestartError("Cannot find WFK or DEN file to restart from.")

        # Move out --> in.
        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

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


class NscfTask(AbinitTask, ProduceGsr):
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
        if not restart_file:
            raise self.RestartError("Cannot find the WFK file to restart from.")

        # Move out --> in.
        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        irdvars = irdvars_for_ext(ext)
        self.strategy.add_extra_abivars(irdvars)

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


class RelaxTask(AbinitTask, ProduceGsr, ProduceHist):
    """
    Task for structural optimizations.
    """
    # TODO possible ScfConvergenceWarning?
    CRITICAL_EVENTS = [
        events.RelaxConvergenceWarning,
    ]

    #EVENT_HANDLERS = [events.DilatmxErrorHandler(), events.TolSymErrorHandler()]

    def _change_structure(self, new_structure):
        """Change the input structure."""
        # Compare new and old structure for logging purpose.
        # TODO: Write method of structure to compare self and other and return a dictionary
        old_structure = self.strategy.abinit_input.structure
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

        self.strategy.abinit_input.set_structure(new_structure)
        #assert self.strategy.abinit_input.structure == new_structure

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

        Structure relaxations can be restarted only if we have the WFK file or the DEN or the GSR file.
        from which we can read the last structure (mandatory) and the wavefunctions (not mandatory but useful).
        Prefer WFK over other files since we can reuse the wavefunctions.

        .. note::

            The problem in the present approach is that some parameters in the input
            are computed from the initial structure and may not be consisten with
            the modification of the structure done during the structure relaxation.
        """
        for ext in ("WFK", "DEN"):
            ofile = self.outdir.has_abiext(ext)
            if ofile:
                irdvars = irdvars_for_ext(ext)
                restart_file = self.out_to_in(ofile)
                break
        else:
            raise self.RestartError("Cannot find the WFK|DEN file to restart from.")

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

        # Read the relaxed structure from the GSR file.
        structure = self.get_final_structure()
                                                           
        # Change the structure.
        self._change_structure(structure)

        # Now we can resubmit the job.
        self.history.info("Will restart from %s", restart_file)
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


class DdeTask(AbinitTask, ProduceDdb):
    """Task for DDE calculations."""

    def get_results(self, **kwargs):
        results = super(DdeTask, self).get_results(**kwargs)
        return results.register_gridfs_file(DDB=(self.outdir.has_abiext("DDE"), "t"))


class DdkTask(AbinitTask, ProduceDdb):
    """Task for DDK calculations."""

    def _on_ok(self):
        super(DdkTask, self)._on_ok()
        # Copy instead of removing, otherwise optic tests fail
        # Fixing this proble requires a rationalization of file extensions.
        #if self.outdir.rename_abiext('1WF', 'DDK') > 0:
        #if self.outdir.copy_abiext('1WF', 'DDK') > 0:
        if self.outdir.symlink_abiext('1WF', 'DDK') > 0:
            raise RuntimeError

    def get_results(self, **kwargs):
        results = super(DdkTask, self).get_results(**kwargs)
        return results.register_gridfs_file(DDK=(self.outdir.has_abiext("DDK"), "t"))


class PhononTask(AbinitTask, ProduceDdb):
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
        """
        Phonon calculations can be restarted only if we have the 1WF file or the 1DEN file.
        from which we can read the first-order wavefunctions or the first order density.
        Prefer 1WF over 1DEN since we can reuse the wavefunctions.
        """
        #self.fix_ofiles()
        for ext in ("1WF", "1DEN"):
            restart_file = self.outdir.has_abiext(ext)
            irdvars = irdvars_for_ext(ext)
            if restart_file:
                break
        else:
            raise self.RestartError("Cannot find the 1WF|1DEN|file to restart from.")

        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

        # Now we can resubmit the job.
        self.history.info("Will restart from %s", restart_file)
        return self._restart()

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

    def make_links(self):
        super(PhononTask, self).make_links()
        # fix the problem that abinit uses hte 1WF extension for the DDK output file but reads it with the irdddk flag
        #if self.indir.has_abiext('DDK'):
        #    self.indir.rename_abiext('DDK', '1WF')


class ScrTask(AbinitTask):
    """Tasks for SCREENING calculations """
    #def inspect(self, **kwargs):
    #    """Plot graph showing the number of q-points computed and the wall-time used"""


class SigmaTask(AbinitTask):
    """
    Tasks for SIGMA calculations. Provides support for in-place restart via QPS files
    """
    CRITICAL_EVENTS = [
        events.QPSConvergenceWarning,
    ]

    def restart(self):
        # G calculations can be restarted only if we have the QPS file 
        # from which we can read the results of the previous step.
        restart_file = self.outdir.has_abiext("QPS")
        if not restart_file:
            raise self.RestartError("Cannot find the QPS file to restart from.")

        self.out_to_in(restart_file)

        # Add the appropriate variable for restarting.
        irdvars = irdvars_for_ext(ext)
        self.strategy.add_extra_abivars(irdvars)

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
        Returns SigresFile object, None if file could not be found or file is not readable.
        """
        sigres_path = self.sigres_path

        if not sigres_path:
            logger.critical("%s didn't produce a SIGRES file in %s" % (self, self.outdir))
            return None

        # Open the GSR file and add its data to results.out
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
            RuntimeError if SIGRES file is not found
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


class BseTask(AbinitTask):
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
                raise self.RestartError("Cannot find BSR|BSC files in %s" % self.indir)

        # Rename HAYDR_SAVE files
        count = 0
        for ext in ("HAYDR_SAVE", "HAYDC_SAVE"):
            ofile = self.outdir.has_abiext(ext)
            if ofile:
                count += 1
                irdvars.update(irdvars_for_ext(ext))
                self.out_to_in(ofile)

        if not count:
            raise self.RestartError("Cannot find the HAYDR_SAVE file to restart from.")

        # Add the appropriate variable for restarting.
        self.strategy.add_extra_abivars(irdvars)

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
    def __init__(self, optic_input, nscf_node, ddk_nodes, workdir=None, manager=None):
        """
        Create an instance of :class:`OpticTask` from an string containing the input.
    
        Args:
            optic_input: string with the optic variables (filepaths will be added at run time).
            nscf_node: The NSCF task that will produce thw WFK file or string with the path of the WFK file.
            ddk_nodes: List of :class:`DdkTask` nodes that will produce the DDK files or list of DDF paths.
            workdir: Path to the working directory.
            manager: :class:`TaskManager` object.
        """
        # Convert paths to FileNodes
        self.nscf_node = Node.as_node(nscf_node)
        self.ddk_nodes = [Node.as_node(n) for n in ddk_nodes]
        assert len(ddk_nodes) == 3
        #print(self.nscf_node, self.ddk_nodes)

        # Use DDK extension instead of 1WF
        deps = {n: "1WF" for n in self.ddk_nodes}
        #deps = {n: "DDK" for n in self.ddk_nodes}
        deps.update({self.nscf_node: "WFK"})

        strategy = OpticInput(optic_input)
        super(OpticTask, self).__init__(strategy=strategy, workdir=workdir, manager=manager, deps=deps)

    def set_workdir(self, workdir, chroot=False):
        """Set the working directory of the task."""
        super(OpticTask, self).set_workdir(workdir, chroot=chroot)
        # Small hack: the log file of optics is actually the main output file. 
        self.output_file = self.log_file

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

    def get_results(self, **kwargs):
        results = super(OpticTask, self).get_results(**kwargs)
        #results.update(
        #"epsilon_infinity":
        #))
        return results


class AnaddbTask(Task):
    """Task for Anaddb runs (post-processing of DFPT calculations)."""

    def __init__(self, anaddb_input, ddb_node,
                 gkk_node=None, md_node=None, ddk_node=None, workdir=None, manager=None):
        """
        Create an instance of :class:`AnaddbTask` from an string containing the input.

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

        super(AnaddbTask, self).__init__(strategy=anaddb_input, workdir=workdir, manager=manager, deps=deps)

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
