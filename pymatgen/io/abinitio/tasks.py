# coding: utf-8
"""Classes defining Abinit calculations."""
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

from pprint import pprint
from atomicfile import AtomicFile
from six.moves import map, zip, StringIO
from monty.termcolor import colored
from monty.serialization import loadfn
from monty.string import is_string, list_strings
from monty.io import FileLock
from monty.collections import AttrDict, Namespace
from monty.functools import lazy_property, return_none_if_raise
from pymatgen.core.units import  Memory #Time,
from pymatgen.util.string_utils import WildCard
from pymatgen.util.num_utils import maxloc
from pymatgen.serializers.json_coders import PMGSONable, json_pretty_dump, pmg_serialize
from .utils import File, Directory, irdvars_for_ext, abi_splitext, abi_extensions, FilepathFixer, Condition
from .strategies import StrategyWithInput, OpticInput
from .qadapters import make_qadapter, QueueAdapter
from .db import DBConnector
from . import abiinspect
from . import events

try:
    from pydispatch import dispatcher
except ImportError:
    pass


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


class GridFsFile(AttrDict):
    def __init__(self, path, fs_id=None, mode="b"):
        super(GridFsFile, self).__init__(path=path, fs_id=fs_id, mode=mode)


class NodeResults(dict, PMGSONable):
    """
    Dictionary used to store the most important results produced by a Node.
    """
    JSON_SCHEMA = {
        "type": "object",
        "properties": {
            "node_id": {"type": "integer", "required": True},
            "node_finalized": {"type": "boolean", "required": True},
            "node_history": {"type": "array", "required": True},
            "node_class": {"type": "string", "required": True},
            "node_name": {"type": "string", "required": True},
            "node_status": {"type": "string", "required": True},
            "in": {"type": "object", "required": True, "description": "dictionary with input parameters"},
            "out": {"type": "object", "required": True, "description": "dictionary with the output results"},
            "exceptions": {"type": "array", "required": True},
            "files": {"type": "object", "required": True},
        },
    }

    @classmethod
    def from_node(cls, node):
        """Initialize an instance of `NodeResults` from a `Node` subclass."""
        kwargs = dict(
            node_id=node.node_id,
            node_finalized=node.finalized,
            node_history=list(node.history),
            node_name=node.name, 
            node_class=node.__class__.__name__,
            node_status=str(node.status),
        )

        return node.Results(node, **kwargs)

    def __init__(self, node, **kwargs):
        super(NodeResults, self).__init__(**kwargs)
        self.node = node

        if "in" not in self: self["in"] = Namespace()
        if "out" not in self: self["out"] = Namespace()
        if "exceptions" not in self: self["exceptions"] = []
        if "files" not in self: self["files"] = Namespace()

    @property
    def exceptions(self):
        return self["exceptions"]

    @property
    def gridfs_files(self):
        """List with the absolute paths of the files to be put in GridFs."""
        return self["files"]

    def add_gridfs_files(self, **kwargs):
        """
        This function registers the files that will be saved in GridFS.
        kwargs is a dictionary mapping the key associated to the file (usually the extension)
        to the absolute path. By default, files are assumed to be in binary form, for formatted files
        one should pass a tuple ("filepath", "t").

        Example::

            add_gridfs(GSR="path/to/GSR.nc", text_file=("/path/to/txt_file", "t"))

        The GSR file is a binary file, whereas text_file is a text file.
        """
        d = {}
        for k, v in kwargs.items():
            mode = "b" 
            if isinstance(v, (list, tuple)): v, mode = v
            d[k] = GridFsFile(path=v, mode=mode)
            
        self["files"].update(d)
        return self

    def push_exceptions(self, *exceptions):
        for exc in exceptions:
            newstr = str(exc)
            if newstr not in self.exceptions:
                self["exceptions"] += [newstr,]

    @pmg_serialize
    def as_dict(self):
        return self.copy()
                                                                                
    @classmethod
    def from_dict(cls, d):
        return cls({k: v for k, v in d.items() if k not in ("@module", "@class")})

    def json_dump(self, filename):
        json_pretty_dump(self.as_dict(), filename)

    @classmethod
    def json_load(cls, filename):
        return cls.from_dict(loadfn(filename))

    def validate_json_schema(self):
        import validictory
        d = self.as_dict()
        try:
            validictory.validate(d, self.JSON_SCHEMA)
            return True
        except ValueError as exc:
            pprint(d)
            print(exc)
            return False

    def update_collection(self, collection):
        """
        Update a mongodb collection.
        """
        node = self.node 
        flow = node if node.is_flow else node.flow

        # Build the key used to store the entry in the document.
        key = node.name
        if node.is_task:
            key = "w" + str(node.pos[0]) + "_t" + str(node.pos[1])
        elif node.is_work:
            key = "w" + str(node.pos)

        db = collection.database

        # Save files with GridFs first in order to get the ID.
        if self.gridfs_files:
            import gridfs
            fs = gridfs.GridFS(db)
            for ext, gridfile in self.gridfs_files.items():
                logger.info("gridfs: about to put file:", str(gridfile))
                # Here we set gridfile.fs_id that will be stored in the mondodb document
                try:
                    with open(gridfile.path, "r" + gridfile.mode) as f:
                        gridfile.fs_id = fs.put(f, filename=gridfile.path)
                except IOError as exc:
                    logger.critical(str(exc))

        if flow.mongo_id is None:
            # Flow does not have a mongo_id, allocate doc for the flow and save its id.
            flow.mongo_id = collection.insert({})
            print("Creating flow.mongo_id", flow.mongo_id, type(flow.mongo_id))

        # Get the document from flow.mongo_id and update it.
        doc = collection.find_one({"_id": flow.mongo_id})
        if key in doc:
            raise ValueError("%s is already in doc!" % key)
        doc[key] = self.as_dict()

        collection.save(doc)
        #collection.update({'_id':mongo_id}, {"$set": doc}, upsert=False)


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

        new.add_gridfs_files(
            run_abi=(task.input_file.path, "t"),
            run_abo=(task.output_file.path, "t"),
        )

        return new


class SparseHistogram(object):

    def __init__(self, items, key=None, num=None, step=None):
        if num is None and step is None:
            raise ValueError("Either num or step must be specified")
        import numpy as np
        from collections import defaultdict, OrderedDict

        values = [key(item) for item in items] if key is not None else items
        start, stop = min(values), max(values)
        if num is None:
            num = int((stop - start) / step)
            if num == 0: num = 1
        mesh = np.linspace(start, stop, num, endpoint=False)

        from monty.bisect import find_le

        hist = defaultdict(list)
        for item, value in zip(items, values):
            # Find rightmost value less than or equal to x.
            # hence each bin contains all items whose value is >= value
            pos = find_le(mesh, value)
            hist[mesh[pos]].append(item)

        #new = OrderedDict([(pos, hist[pos]) for pos in sorted(hist.keys(), reverse=reverse)])
        self.binvals = sorted(hist.keys())
        self.values = [hist[pos] for pos in self.binvals]
        self.start, self.stop, self.num = start, stop, num

    from pymatgen.util.plotting_utils import add_fig_kwargs
    @add_fig_kwargs
    def plot(self, **kwargs):
        """
        Plot the histogram with matplotlib, returns `matplotlib figure
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()

        ax = fig.add_subplot(1,1,1)

        yy = [len(v) for v in self.values]
        ax.plot(self.binvals, yy, **kwargs)
        return fig


import unittest
class SparseHistogramTest(unittest.TestCase):
    def test_sparse(self):
        items = [1, 2, 2.9, 4]
        hist = SparseHistogram(items, step=1)
        assert hist.binvals == [1.0, 2.0, 3.0] 
        assert hist.values == [[1], [2, 2.9], [4]]
        #hist.plot()

        hist = SparseHistogram([iv for iv in enumerate(items)], key=lambda t: t[1], step=1)
        assert hist.binvals == [1.0, 2.0, 3.0] 
        assert hist.values == [[(0, 1)], [(1, 2), (2, 2.9)], [(3, 4)]]


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
                err_msg = "Wrong YAML doc:\n%s\n\nException" % (doc.text, sexc)
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
                condition: `Condition` object with operators expressed with a Mongodb-like syntax
                key: Selects the sub-dictionary on which condition is applied, e.g. key="vars"
                    if we have to filter the configurations depending on the values in vars
        """
        new_confs = []

        for conf in self:
            # Select the object on which condition is applied
            obj = conf if key is None else AttrDict(conf[key])
            add_it = condition(obj=obj)
            #if key is "vars": print("conf", conf, "added:", add_it)
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

    def sort_by_mem_per_proc(self, reverse=False):
        """
        Sort the configurations in place so that conf with lowest memory per core
        appears in the first positions.
        """
        # Avoid sorting if mem_per_cpu is not available.
        if any(c.mem_per_proc > 0.0 for c in self):
            self._confs.sort(key=lambda c: c.mem_per_proc, reverse=reverse)

    def filter(self, qadapter):
        """Return a new list of configurations that can be executed on the `QueueAdapter` qadapter."""
        new_confs = [pconf for pconf in self if qadapter.can_run_pconf(pconf)]
        return self.__class__(info=self.info, confs=new_confs)

    #def hist_efficiency(self, step=0.1):
    #    return sparse_histogram(self._confs, key=lambda c: c.efficiency, step=step)
    #def hist_speedup(self, step=1.0):
    #    return sparse_histogram(self._confs, key=lambda c: c.speedup, step=step)
    #def hist_memory(self, step=102.4):
    #    return sparse_histogram(self._confs, key=lambda c: c.speedup, step=step)

    def select_optimal_conf(self, policy, max_ncpus):
        """
        Find the optimal configuration according to the :class:`TaskPolicy` policy.
        """
        # Make a copy since we are gonna change the object in place.
        hints = ParalHints(self.info, confs=[c for c in self if c.num_cores <= max_ncpus])

        if False:
            # Mapping property --> options passed to sparse_histogram
            opts = dict(speedup=dict(step=1.0),
                        efficiency=dict(step=0.1), 
                        memory=dict(memory=102.4))

            confs = self._confs
            for priority in ["speedup", "efficiency", "memory"]:
                histogram = SparseHistogram(confs, key=lambda c: getattr(c, priority), **opts[priority])
                pos = 0 if priority == "memory" else -1
                confs = histogram.values[pos]

        # First select the configurations satisfying the condition specified by the user (if any)
        if policy.condition:
            #logger.info("condition %s" % str(policy.condition))
            hints.select_with_condition(policy.condition)
            #logger.info("after condition %s" % str(hints))

            # If no configuration fullfills the requirements, 
            # we return the one with the highest speedup.
            if not hints:
                logger.warning("empty list of configurations after policy.condition")
                hints = self.copy()
                hints.sort_by_speedup()
                return hints[-1].copy()

        # Now filter the configurations depending on the values in vars
        if policy.vars_condition:
            logger.info("vars_condition %s" % str(policy.vars_condition))
            hints.select_with_condition(policy.vars_condition, key="vars")
            logger.info("After vars_condition %s" % str(hints))

            # If no configuration fullfills the requirements,
            # we return the one with the highest speedup.
            if not hints:
                logger.warning("empty list of configurations after policy.vars_condition")
                hints = self.copy()
                hints.sort_by_speedup()
                return hints[-1].copy()

        hints.sort_by_speedup()

        logger.info('speedup hints: \n' + str(hints) + '\n')

        #hints.sort_by_efficiency()

        #logger.info('efficiency hints: \n' + str(hints) + '\n')

        # Find the optimal configuration according to policy.mode.
        #if policy.mode in ["default", "aggressive"]:
        #    hints.sort_by_spedup()
        #elif policy.mode == "conservative":
        #    hints.sort_by_efficiency()
        #    # Remove num_cores == 1
        #    hints.pop(num_cores==1)
        #else:
        #    raise ValueError("Wrong value for policy.mode: %s" % str(policy.mode))
        #if not hints:

        # Return a copy of the configuration.
        optimal = hints[-1].copy()
        logger.info("Will relaunch the job with optimized parameters:\n %s" % optimal)

        return optimal


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

    def __init__(self, **kwargs):
        """
        Args:
            autoparal: Value of ABINIT autoparal input variable. 0 to disable the autoparal feature (default)
            mode: Select the algorith to select the optimal configuration for the parallel execution.
                Possible values: ["default", "aggressive", "conservative"]
            condition: condition used to filter the autoparal configuration (Mongodb-like syntax)
            vars_condition: condition used to filter the list of Abinit variables suggested by autoparal (Mongodb-like syntax)
        """
        self.autoparal = kwargs.pop("autoparal", 1)
        #self.autoparal = kwargs.pop("autoparal", 0)
        #self.mode = kwargs.pop("mode", "default")
        #self.max_ncpus = kwargs.pop("max_ncpus", None)
        self.condition = Condition(kwargs.pop("condition", {}))
        self.vars_condition = Condition(kwargs.pop("vars_condition", {}))

        #if self.autoparal and self.max_ncpus is None:
        #    raise ValueError("When autoparal is not zero, max_ncpus must be specified.")

        if kwargs:
            raise ValueError("Found invalid keywords in policy section:\n %s" % str(kwargs.keys()))

    def __str__(self):
        lines = []
        app = lines.append
        for k, v in self.__dict__.items():
            if k.startswith("_"):
                continue
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
# Main options:

policy: 

qadapters:  
    # List of qadapters objects (mandatory)
    -  # qadapter_1
    -  # qadapter_2

db_connector: # Connection to MongoDB database (optional)

##########################################
# Individual entries are documented below:
##########################################

"""
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
            #raise RuntimeError(msg + str(exc))

    @classmethod
    def from_string(cls, s):
        """Create an instance from string s containing a YAML dictionary."""
        return cls.from_dict(yaml.load(s))

    @classmethod
    def from_dict(cls, d):
        """Create an instance from a dictionary."""
        return cls(**{k: v for k, v in d.items() if k in cls.ENTRIES}) 

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
        if len(qads) != 1:
            raise NotImplementedError("For the time being multiple qadapters are not supported! Please use one adapter")

        # Order qdapters according to priority.
        qads = sorted(qads, key=lambda q: q.priority)
        priorities = [q.priority for q in qads]
        if len(priorities) != len(set(priorities)):
            raise ValueError("Two or more qadapters have same priority. This is not allowed. Check taskmanager.yml")

        self._qads, self._qid = qads, 0

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
        """List of :class:`QueueAdapter` objects."""
        return self._qads

    @property
    def qadapter(self):
        """The qadapter used to submit jobs."""
        return self._qads[self._qid]

    def select_qadapter(self, pconf):
        """
        Select the qadatper to submit a run with the parallel configuration pconf
        Return False if no qadapter could be found.
        """
        scores = [q.get_score(pconf) for q in self._qads]
        #print("scores", scores)
        if all(sc < 0 for sc in scores): return False
        self._qid = maxloc(scores)

        # Change the number of MPI/OMP cores.
        self.set_mpi_procs(pconf.mpi_procs)
        if self.has_omp: self.set_omp_threads(pconf.omp_threads)
                                                                      
        self.set_mem_per_proc(pconf.mem_per_proc)

        return True

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        #app("[Task policy]\n%s" % str(self.policy))

        for i, qad in enumerate(self.qads):
            app("[Qadapter %d]\n%s" % (i, str(qad)))
        app("Qadapter selected: %d" % self._qid)

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
            return task.job_file.path

    def launch(self, task):
        """
        Build the input files and submit the task via the `Qadapter` 

        Args:
            task:
                `TaskObject`
        
        Returns:
            Process object.
        """
        # Build the task 
        task.build()

        # Submit the task and save the queue id.
        task.set_status(task.S_SUB)
        script_file = self.write_jobfile(task)

        qjob, process = self.qadapter.submit_to_queue(script_file)
        task.set_qjob(qjob)

        return process

    def get_collection(self, **kwargs):
        """Return the MongoDB collection used to store the results."""
        return self.db_connector.get_collection(**kwargs)

    def increase_resources(self):
        # with GW calculations in mind with GW mem = 10, 
        # the response fuction is in memory and not distributed
        # we need to increas memory if jobs fail ...
        return self.qadapter.more_mem_per_proc()


# The code below initializes a counter from a file when the module is imported 
# and save the counter's updated value automatically when the program terminates 
# without relying on the application making an explicit call into this module at termination.
conf_dir = os.path.join(os.getenv("HOME"), ".abinit", "abipy")

if not os.path.exists(conf_dir):
    os.makedirs(conf_dir)

_COUNTER_FILE = os.path.join(conf_dir, "nodecounter")
del conf_dir

try:
    with open(_COUNTER_FILE, "r") as _fh:
        _COUNTER = int(_fh.read())

except IOError:
    _COUNTER = -1


def get_newnode_id():
    """
    Returns a new node identifier used both for `Task` and `Work` objects.

    .. warnings:
        The id is unique inside the same python process so be careful when 
        Works and Tasks are constructed at run-time or when threads are used.
    """
    global _COUNTER
    _COUNTER += 1
    return _COUNTER


def save_lastnode_id():
    """Save the id of the last node created."""
    with FileLock(_COUNTER_FILE):
        with AtomicFile(_COUNTER_FILE, mode="w") as fh:
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
    This file is needed to start another `Task` or another `Work`.
    """
    def __init__(self, ext, path):
        """
        Args:
            ext: ABINIT file extension
            path: (asbolute) filepath
        """
        if ext not in abi_extensions():
            raise ValueError("Extension %s has not been registered in the internal database" % str(ext))

        self.ext = ext
        self.file = File(path)

    @classmethod
    def from_file(cls, filepath):
        """Build a :class:`Product` instance from a filepath."""
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
    that are used by the other nodes (`Task` or `Work`) to start the calculation.
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
            node: The task or the worfklow associated to the dependency or string with a filepath.
            exts: Extensions of the output files that are needed for running the other tasks.
        """
        self._node = Node.as_node(node)

        if exts and is_string(exts):
            exts = exts.split()

        self.exts = exts or []

    def __hash__(self):
        return hash(self._node)

    def __repr__(self):
        return "node %s will produce: %s " % (repr(self.node), repr(self.exts))

    def __str__(self):
        return "node %s will produce: %s " % (str(self.node), str(self.exts))

    @property
    def info(self):
        return str(self.node)

    @property
    def node(self):
        """The :class:`Node` associated to the dependency."""
        return self._node

    @property
    def status(self):
        """The status of the dependency, i.e. the status of the :class:`Node`."""
        return self.node.status

    @lazy_property
    def products(self):
        """List of output files produces by self."""
        _products = []
        for ext in self.exts:
            prod = Product(ext, self.node.opath_from_ext(ext))
            _products.append(prod)

        return _products

    def connecting_vars(self):
        """
        Returns a dictionary with the variables that must be added to the 
        input file in order to connect this :class:`Node` to its dependencies.
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


def _2attrs(item):
        return item if item is None or isinstance(list, tuple) else (item,)


class Status(int):
    """This object is an integer representing the status of the `Node`."""

    # Possible status of the node. See monty.termocolor for the meaning of color, on_color and attrs.
    _STATUS_INFO = [
        #(value, name, color, on_color, attrs)
        (1,  "Initialized",   None     , None, None),         # Node has been initialized
        (2,  "Locked",        None     , None, None),         # Task is locked an must be explicitly unlocked by an external subject (Work).
        (3,  "Ready",         None     , None, None),         # Node is ready i.e. all the depencies of the node have status S_OK
        (4,  "Submitted",     "blue"   , None, None),         # Node has been submitted (The `Task` is running or we have started to finalize the Work)
        (5,  "Running",       "magenta", None, None),         # Node is running.
        (6,  "Done",          None     , None, None),         # Node done, This does not imply that results are ok or that the calculation completed successfully
        (7,  "AbiCritical",   "red"    , None, None),         # Node raised an Error by ABINIT.
        (8,  "QCritical",     "red"    , "on_white", None),   # Node raised an Error by submitting submission script, or by executing it
        (9,  "Unconverged",   "red"    , "on_yellow", None),  # This usually means that an iterative algorithm didn't converge.
        (10, "Error",         "red"    , None, None),         # Node raised an unrecoverable error, usually raised when an attempt to fix one of other types failed.
        (11, "Completed",     "green"  , None, None),         # Execution completed successfully.
        #(11, "Completed",     "green"  , None, "underline"),   
    ]
    _STATUS2STR = collections.OrderedDict([(t[0], t[1]) for t in _STATUS_INFO])
    _STATUS2COLOR_OPTS = collections.OrderedDict([(t[0], {"color": t[2], "on_color": t[3], "attrs": _2attrs(t[4])}) for t in _STATUS_INFO])

    def __repr__(self):
        return "<%s: %s, at %s>" % (self.__class__.__name__, str(self), id(self))

    def __str__(self):
        """String representation."""
        return self._STATUS2STR[self]

    @classmethod
    def as_status(cls, obj):
        """Convert obj into Status."""
        return obj if isinstance(obj, cls) else cls.from_string(obj)

    @classmethod
    def from_string(cls, s):
        """Return a `Status` instance from its string representation."""
        for num, text in cls._STATUS2STR.items():
            if text == s:
                return cls(num)
        else:
            raise ValueError("Wrong string %s" % s)

    @property
    def is_critical(self):
        """True if status is critical."""
        return str(self) in ("AbiCritical", "QCritical", "Uncoverged", "Error") 

    @property
    def color_opts(self):
        return self._STATUS2COLOR_OPTS[self]

    @property
    def colored(self):
        """Return colorized text used to print the status if the stream supports it."""
        return colored(str(self), **self.color_opts) 


class Node(six.with_metaclass(abc.ABCMeta, object)):
    """
    Abstract base class defining the interface that must be 
    implemented by the nodes of the calculation.

    Nodes are hashable and can be tested for equality
    (hash uses the node identifier, whereas eq uses workdir).
    """
    Results = NodeResults

    # Possible status of the node.
    S_INIT = Status.from_string("Initialized")
    S_LOCKED = Status.from_string("Locked")
    S_READY = Status.from_string("Ready")
    S_SUB = Status.from_string("Submitted")
    S_RUN = Status.from_string("Running")
    S_DONE = Status.from_string("Done")
    S_ABICRITICAL = Status.from_string("AbiCritical")
    S_QCRITICAL = Status.from_string("QCritical")
    S_UNCONVERGED = Status.from_string("Unconverged")
    S_ERROR = Status.from_string("Error")
    S_OK = Status.from_string("Completed")

    ALL_STATUS = [
        S_INIT,
        S_LOCKED,
        S_READY,
        S_SUB,
        S_RUN,
        S_DONE,
        S_ABICRITICAL,
        S_QCRITICAL,
        S_UNCONVERGED,
        S_ERROR,
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

        self._status = self.S_INIT

    def __eq__(self, other):
        if not isinstance(other, Node): return False

        #return self.node_id == other.node_id and 
        return (self.__class__ == other.__class__ and 
                self.workdir == other.workdir)

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

    @classmethod
    def as_node(cls, obj):
        """
        Convert obj into a Node instance.

        Return:
            obj if obj is a Node instance,
            cast obj to :class:`FileNode` instance of obj is a string.
            None if obj is None
        """
        if isinstance(obj, cls):
            return obj
        elif is_string(obj):
            # Assume filepath.
            return FileNode(obj)
        elif obj is None:
            return obj
        else:
            raise TypeError("Don't know how to convert %s to Node instance." % obj)

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
        """True if the `Work` has been finalized."""
        return self._finalized

    @finalized.setter
    def finalized(self, boolean):
        self._finalized = boolean
        self.history.append("Finalized on %s" % time.asctime())

    @property
    def str_history(self):
        """String representation of history."""
        return "\n".join(self.history)

    @property
    def is_file(self):
        """True if this node is a file"""
        return isinstance(self, FileNode)

    @property
    def is_task(self):
        """True if this node is a Task"""
        return isinstance(self, Task)

    @property
    def is_work(self):
        """True if this node is a Work"""
        from .works import Work
        return isinstance(self, Work)

    @property
    def is_flow(self):
        """True if this node is a Flow"""
        from .flows import Flow
        return isinstance(self, Flow)

    @property
    def has_subnodes(self):
        """True if self contains sub-nodes e.g. `Work` object."""
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
            deps: List of :class:`Dependency` objects specifying the dependencies of the node.
                  or dictionary mapping nodes to file extensions e.g. {task: "DEN"}
        """
        if isinstance(deps, collections.Mapping):
            # Convert dictionary into list of dependencies.
            deps = [Dependency(node, exts) for node, exts in deps.items()]

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
            deps: List of :class:`Dependency` objects specifying the  dependencies of the node.
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

    def get_parents(self):
        """Return the list of nodes in the :class:`Flow` required by this :class:`Node`"""
        parents = []
        for work in self.flow:
            if self.depends_on(work): parents.append(work)
            for task in work:
                if self.depends_on(task): parents.append(task)
        return parents

    def get_children(self):
        """Return the list of nodes in the :class:`Flow` that depends on this :class:`Node`"""
        children = []
        for work in self.flow:
            if work.depends_on(self): children.append(work)
            for task in work:
                if task.depends_on(self): children.append(task)
        return children

    def str_deps(self):
        """Return the string representation of the dependencies of the node."""
        lines = []
        app = lines.append

        app("Dependencies of node %s:" % str(self))
        for i, dep in enumerate(self.deps):
            app("%d) %s, status=%s" % (i, dep.info, str(dep.status)))

        return "\n".join(lines)

    def set_user_info(self, *args, **kwargs):
        """
        Store additional info provided by the user in self.user_info

        .. warning::

            The objects stored in the dict must support pickle.
        """
        if not hasattr(self, "_user_info"): self._user_info = {}
        self._user_info.update(*args, **kwargs)

    @property
    def user_info(self):
        """Returns an :class:`AttrDict` with the variables stored in self._user_info."""
        try:
            return AttrDict(**self._user_info)
        except AttributeError:
            return {}

    #@abc.abstractmethod
    #def set_status(self, status,  info_msg=None):
    #    """
    #    Set and return the status of the None
    #                                                                                     
    #    Args:
    #        status: Status object or string representation of the status
    #        info_msg: string with human-readable message used in the case of errors (optional)
    #    """

    @abc.abstractproperty
    def status(self):
        """The status of the `Node`."""

    @abc.abstractmethod
    def check_status(self):
        """Check the status of the `Node`."""


class FileNode(Node):
    """
    A Node that consists of a file. May be not yet existing

    Mainly used to connect `Tasks` to external files produced in previous runs
    """
    def __init__(self, filename):
        super(FileNode, self).__init__()
        self.filepath = os.path.abspath(filename)

        # Directories with input|output|temporary data.
        self.workdir = os.path.dirname(self.filepath)

        self.indir = Directory(self.workdir)
        self.outdir = Directory(self.workdir)
        self.tmpdir = Directory(self.workdir)

    @property
    def products(self):
        return [Product.from_file(self.filepath)]

    def opath_from_ext(self, ext):
        return self.filepath

    @property
    def status(self):
        return self.S_OK if os.path.exists(self.filepath) else self.S_ERROR

    def check_status(self):
        return self.status

    def get_results(self, **kwargs):
        results = super(FileNode, self).get_results(**kwargs)
        #results.add_gridfs_files(self.filepath=self.filepath)
        return results


class TaskError(Exception):
    """Base Exception for `Task` methods"""


class TaskRestartError(TaskError):
    """Exception raised while trying to restart the `Task`."""


class Task(six.with_metaclass(abc.ABCMeta, Node)):
    """A Task is a node that performs some kind of calculation."""
    # Use class attributes for TaskErrors so that we don't have to import them.
    Error = TaskError
    RestartError = TaskRestartError

    # List of `AbinitEvent` subclasses that are tested in the not_converged method. 
    # Subclasses should provide their own list if they need to check the converge status.
    CRITICAL_EVENTS = [
    ]

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
                  None means that this obj has no dependency.
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

        # Use to compute the wall-time
        self.start_datetime, self.stop_datetime = None, None

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
        """Set the `TaskManager` to use to launch the Task."""
        self.manager = manager.deepcopy()

        # TODO
        # Select adapters associated to the Task class
        #keep = []
        #for i, qad in enumerate(self.manager.qads):
        #    if self.__class__.__name__ in qad.task_classes:
        #        keep.append(i)
        #if keep:
        #    self._qads = [self._qads[i] for i in keep]
        #    self._qid = 0

    @property
    def work(self):
        """The Work containing this `Task`."""
        return self._work

    def set_work(self, work):
        """Set the Work associated to this `Task`."""
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
                return (self.work.pos, i)
        raise ValueError("Cannot find the position of %s in flow %s" % (self, self.flow))

    @property
    def pos_str(self):
        """String representation of self.pos"""
        return "w" + str(self.pos[0]) + "_t" + str(self.pos[1])

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

    @property
    def is_completed(self):
        """True if the task has been executed."""
        return self.status >= self.S_DONE

    @property
    def can_run(self):
        """The task can run if its status is < S_SUB and all the other dependencies (if any) are done!"""
        all_ok = all([stat == self.S_OK for stat in self.deps_status])
        #print("can_run: all_ok ==  ",all_ok)
        return self.status < self.S_SUB and all_ok

    def not_converged(self):
        """Return True if the calculation is not converged."""
        report = self.get_event_report()
        return report.filter_types(self.CRITICAL_EVENTS)

    def run_etime(self):
        """
        String with the wall-time

        ...note::

            The clock starts when self.status becomes S_RUN.
            thus run_etime does not correspond to the effective wall-time.
        """
        # FIXME: This does not work as expected!
        s = "None"
        if self.start_datetime is not None:
            stop = self.stop_datetime
            if stop is None:
                stop = datetime.datetime.now()

            # Compute time-delta, convert to string and remove microseconds (in any)
            s = str(stop - self.start_datetime)
            microsec = s.find(".")
            if microsec != -1: s = s[:microsec]

        return s

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
        # Read timing data.
        #self.read_timing()
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

    def _restart(self, no_submit=False):
        """
        Called by restart once we have finished preparing the task for restarting.

        Return True if task has been restarted
        """
        self.set_status(self.S_READY, info_msg="Restarted on %s" % time.asctime())

        # Increase the counter.
        self.num_restarts += 1
        self.history.append("Restarted on %s, num_restarts %d" % (time.asctime(), self.num_restarts))

        if not no_submit:
            # Remove the lock file
            self.start_lockfile.remove()
            # Relaunch the task.
            fired = self.start()
            if not fired:
                self.history.append("[%s], restart failed" % time.asctime())
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

        self.set_status(self.S_INIT, info_msg="Reset on %s" % time.asctime())
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

    def set_status(self, status, info_msg=None):
        """
        Set and return the status of the task.

        Args:
            status: Status object or string representation of the status
            info_msg: string with human-readable message used in the case of errors (optional)
        """
        status = Status.as_status(status)

        changed = True
        if hasattr(self, "_status"):
            changed = (status != self._status)

        self._status = status

        if status == self.S_RUN:
            # Set start_datetime when the task enters S_RUN
            if self.start_datetime is None:
                self.start_datetime = datetime.datetime.now()

        # Add new entry to history only if the status has changed.
        if changed:
            if status == self.S_SUB: 
                self._submission_time = time.time()
                self.history.append("Submitted on %s" % time.asctime())

            if status == self.S_OK:
                self.history.append("Completed on %s" % time.asctime())

            if status == self.S_ABICRITICAL:
                self.history.append("Error info:\n %s" % str(info_msg))

        if status == self.S_DONE:
            self.stop_datetime = datetime.datetime.now()

            # Execute the callback
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
        This function checks the status of the task by inspecting the output and the
        error files produced by the application and by the queue manager.

        The process
        1) see it the job is blocked
        2) see if an error occured at submitting the job the job was submitted, TODO these problems can be solved
        3) see if there is output
        4) see if abinit reports problems
        5) see if both err files exist and are empty
        6) no output and no err files, the job must still be running
        7) try to find out what caused the problems
        8) there is a problem but we did not figure out what ...
        9) the only way of landing here is if there is a output file but no err files...
        """

        # 1) A locked task can only be unlocked by calling set_status explicitly.
        # an errored task, should not end up here but just to be sure
        black_list = [self.S_LOCKED, self.S_ERROR]
        if self.status in black_list:
            return

        # 2) Check the returncode of the process (the process of submitting the job) first.
        # this point type of problem should also be handled by the scheduler error parser
        if self.returncode != 0:
            # The job was not submitter properly
            info_msg = "return code %s" % self.returncode
            return self.set_status(self.S_QCRITICAL, info_msg=info_msg)           

#        err_msg = None
#=======
#            if not self.stderr_file.exists and not self.qerr_file.exists:
#                # The job is still in the queue.
#                return self.status
#
#            else:
#                # Analyze the standard error of the executable:
#                if self.stderr_file.exists:
#                    err_msg = self.stderr_file.read()
#                    if err_msg:
#                        logger.critical("%s: executable stderr:\n %s" % (self, err_msg))
#                        return self.set_status(self.S_ERROR, info_msg=err_msg)
#
#                # Analyze the error file of the resource manager.
#                if self.qerr_file.exists:
#                    err_msg = self.qerr_file.read()
#                    if err_msg:
#                        logger.critical("%s: queue stderr:\n %s" % (self, err_msg))
#                        return self.set_status(self.S_ERROR, info_msg=err_msg)
#
#                return self.status
#
#        # Check if the run completed successfully.
#        report = self.get_event_report()
#
#        if report.run_completed:
#            # Check if the calculation converged.
#            not_ok = self.not_converged()

#            if not_ok:
#                return self.set_status(self.S_UNCONVERGED)
#            else:
#                return self.set_status(self.S_OK)

#       # This is the delicate part since we have to discern among different possibilities:
        #
        # 1) Calculation stopped due to an Abinit Error or Bug.
        #
        # 2) Segmentation fault that (by definition) was not handled by ABINIT.
        #    In this case we check if the ABINIT standard error is not empty.
        #    hoping that nobody has written to stderr (e.g. libraries in debug mode)
        #
        # 3) Problem with the resource manager and/or the OS (walltime error, resource error, phase of the moon ...)
        #    In this case we check if the error file of the queue manager is not empty.
        #    Also in this case we *assume* that there's something wrong if the stderr of the queue manager is not empty
        # 
        # 4) Calculation is still running!
        #
        # Point 2) and 3) are the most complicated since there's no standard!

        # 1) Search for possible errors or bugs in the ABINIT **output** file.
#        if report.errors or report.bugs:
#            logger.critical("%s: Found Errors or Bugs in ABINIT main output!" % self)
#            return self.set_status(self.S_ERROR, info_msg=str(report.errors) + str(report.bugs))

        # 2) Analyze the stderr file for Fortran runtime errors.
#       >>>>>>> pymatgen-matteo/master

        err_msg = None
        if self.stderr_file.exists:
            err_msg = self.stderr_file.read()

        err_info = None
        if self.qerr_file.exists:
            err_info = self.qerr_file.read()

        # 3) Start to check if the output file has been created.
        if self.output_file.exists:
            try:
                report = self.get_event_report()
            except Exception as exc:
                info_msg = "%s exception while parsing event_report:\n%s" % (self, exc)
                logger.critical(info_msg)
                return self.set_status(self.S_ABICRITICAL, info_msg=info_msg)

            if report.run_completed:
                # Check if the calculation converged.
                not_ok = self.not_converged()
                if not_ok:
                    return self.set_status(self.S_UNCONVERGED)
                else:
                    return self.set_status(self.S_OK)

            # 4)
            if report.errors or report.bugs:
                if report.errors:
                    logger.debug('"Found errors in report')
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
                # Abinit reports problems
                logger.critical("%s: Found Errors or Bugs in ABINIT main output!" % self)
                info_msg = str(report.errors) + str(report.bugs)
                return self.set_status(self.S_ABICRITICAL, info_msg=info_msg)
                # The job is unfixable due to ABINIT errors

            # 5)
            if self.stderr_file.exists and not err_info:
                if self.qerr_file.exists and not err_msg:
                    # there is output and no errors
                    # Check if the run completed successfully.
#                    if report.run_completed:
#                        # Check if the calculation converged.
#                        not_ok = self.not_converged()
#                        if not_ok:
#                            return self.set_status(self.S_UNCONVERGED)
#                            # The job finished but did not converge
#                        else:
#                            return self.set_status(self.S_OK)
#                            # The job finished properly

                    return self.set_status(self.S_RUN)
                    # The job still seems to be running

        # 6)
        if not self.output_file.exists:
            logger.debug("output_file does not exists")
            if not self.stderr_file.exists and not self.qerr_file.exists:     # No output at all
                return self.status
                # The job is still in the queue.

        # 7) Analyze the files of the resource manager and abinit and execution err (mvs)
        if self.qerr_file.exists:
            from pymatgen.io.abinitio.scheduler_error_parsers import get_parser
            scheduler_parser = get_parser(self.manager.qadapter.QTYPE, err_file=self.qerr_file.path,
                                          out_file=self.qout_file.path, run_err_file=self.stderr_file.path)
            scheduler_parser.parse()

            if scheduler_parser.errors:
                # the queue errors in the task
                logger.debug('scheduler errors found:')
                logger.debug(str(scheduler_parser.errors))
                self.queue_errors = scheduler_parser.errors
                return self.set_status(self.S_QCRITICAL)
                # The job is killed or crashed and we know what happened
            else:
                if len(err_info) > 0:
                    logger.debug('found unknown queue error: %s' % str(err_info))
                    return self.set_status(self.S_QCRITICAL, info_msg=err_info)
                    # The job is killed or crashed but we don't know what happened
                    # it is set to QCritical, we will attempt to fix it by running on more resources

        # 8) analizing the err files and abinit output did not identify a problem
        # but if the files are not empty we do have a problem but no way of solving it:
        if err_msg is not None and len(err_msg) > 0:
            logger.debug('found error message:\n %s' % str(err_msg))
            return self.set_status(self.S_QCRITICAL, info_msg=err_info)
            # The job is killed or crashed but we don't know what happend
            # it is set to QCritical, we will attempt to fix it by running on more resources

        # 9) if we still haven't returned there is no indication of any error and the job can only still be running
        # but we should actually never land here, or we have delays in the file system ....
        # print('the job still seems to be running maybe it is hanging without producing output... ')

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

        ..warning:
            
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
            source: "output" for the main output file,"log" for the log file.

        Returns:
            :class:`EventReport` instance or None if the main output file does not exist.
        """
        ofile = {
            "output": self.output_file,
            "log": self.log_file}[source]

        if not ofile.exists:
            return None

        parser = events.EventsParser()
        try:
            return parser.parse(ofile.path)
        except parser.Error as exc:
            # Return a report with an error entry with info on the exception.
            logger.critical("%s: Exception while parsing ABINIT events:\n %s" % (ofile, str(exc)))
            self.set_status(self.S_ABICRITICAL, info_msg=str(exc))
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

    def clean_outfiles(self):
        # Remove all files in tmpdir.
        self.tmpdir.clean()

        # Remove files in outdir that are not needed by other nodes
        # Find the file extensions that should be preserved
        exts_toclean = ("WFK", "SUS", "SCR")

        keep_exts = set()
        for node in self.get_children():
            if node.status != self.S_OK:
                i = [d.node for d in node.deps].index(self)
                keep_exts.update(deps[i].exts)
        print("keep_exts: ", keep_exts)
        #self.outdir.remove_abiexts()

        # Remove files in the outdir of the other tasks if the dependency has been fulfilled.
        #for node in self.get_parents():
        #   for child in node.get_children():

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
            logger.warning("Found lock file: %s" % self.start_lockfile.relpath)
            return 0

        self.start_lockfile.write("Started on %s" % time.asctime())

        self.build()
        self._setup()

        # Add the variables needed to connect the node.
        for d in self.deps:
            cvars = d.connecting_vars()
            logger.debug("Adding connecting vars %s " % cvars)
            self.strategy.add_extra_abivars(cvars)

        # Automatic parallelization
        if kwargs.pop("autoparal", True) and hasattr(self, "autoparal_run"):
            try:
                self.autoparal_run()
            except:
                # Log the exception and continue with the parameters specified by the user.
                logger.critical("autoparal_fake_run raised:\n%s" % straceback())
                self.set_status(self.S_ABICRITICAL)
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

    @property
    def pconfs(self):
        """List of autoparal configurations."""
        try:
            return self._pconfs
        except AttributeError:
            return None

    def set_pconfs(self, pconfs):
        """Set the list of autoparal configurations."""
        self._pconfs = pconfs

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
            #self.all_autoparal_confs = pconfs
            logger.info('speedup hints: \n' + str(pconfs) + '\n')
            # print("pconfs", pconfs)
        except parser.Error:
            logger.critical("Error while parsing Autoparal section:\n%s" % straceback())
            return 2

        ######################################################
        # Select the optimal configuration according to policy
        ######################################################
        self.set_pconfs(pconfs)
        optconf = pconfs.select_optimal_conf(policy, max_ncpus)
        #print("optimal autoparal conf:\n %s" % optconf)

        # Select the partition on which we'll be running and set MPI/OMP cores.
        optconfs = [optconf]
        for i, c in enumerate(optconfs):
            if self.manager.select_qadapter(c):
                optconf = optconfs[i]
                break
        else:
            raise RuntimeError("Cannot find qadapter for this run!")

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

    def get_ibz(self, ngkpt=None, shiftk=None):
        """
        Returns:
        """
        logger.info("in get_ibz")

        #########################################
        # Run ABINIT in sequential to get the IBZ
        #########################################

        # Set the variables for automatic parallelization
        ibz_vars = dict(prtkpt=-2)
        if ngkpt is not None: ibz_vars["ngkpt"] = ngkpt
        if shiftk is not None:
            import numpy as np
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
        print('reset_from_scatch', self)
        self.output_file.remove()
        self.log_file.remove()
        self.stderr_file.remove()
        self.start_lockfile.remove()

        return self._restart(no_submit=True)

    def fix_abicritical(self):
        """
        method to fix crashes/error caused by abinit
        currently:
            try to rerun with more resources, last resort if all else fails
        ideas:
            upon repetative no converging iscf > 2 / 12
        """
        # the crude, no idea what to do but this may work, solution.
        if self.manager.increase_resources():
            self.reset_from_scratch()
            return True
        else:
            self.set_status(self.S_ERROR, info_msg='could not increase resources any further')
            return False


# TODO
# Enable restarting capabilites:
# Before doing so I need:
#   1) Preliminary standardization of the ABINT events and critical WARNINGS (YAML)
#   2) Change the parser so that we can use strings in the input file.
#      We need this change for restarting structural relaxations so that we can read 
#      the initial structure from file.

class ProduceGsr(object):
    """
    Mixin class for AbinitTasks producing a GSR file.
    Provice the method `open_gsr` that reads and return a GSR file.
    """
    def open_gsr(self):
        """
        Open the GSR file located in the in self.outdir.
        Returns :class:`GsrFile` object, None if file could not be found or file is not readable.
        """
        gsr_path = self.outdir.has_abiext("GSR")
        if not gsr_path:
            if self.status == self.S_OK:
                logger.critical("%s reached S_OK but didn't produce a GSR file in %s" % (self, self.outdir))
            return None

        # Open the GSR file and add its data to results.out
        from abipy.electrons.gsr import GsrFile
        try:
            return GsrFile(gsr_path)
        except Exception as exc:
            logger.critical("Exception while reading GSR file at %s:\n%s" % (gsr_path, str(exc)))
            return None


class ProduceDdb(object):
    """
    Mixin class for AbinitTasks producing a DDB file.
    Provice the method `open_ddb` that reads and return a Ddb file.
    """
    def open_ddb(self):
        """
        Open the DDB file located in the in self.outdir.
        Returns :class:`DdbFile` object, None if file could not be found or file is not readable.
        """
        ddb_path = self.outdir.has_abiext("DDB")
        if not ddb_path:
            if self.status == self.S_OK:
                logger.critical("%s reached S_OK but didn't produce a DDB file in %s" % (self, self.outdir))
            return None

        # Open the GSR file and add its data to results.out
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
        return self._restart()

    def inspect(self, **kwargs):
        """
        Plot the SCF cycle results with matplotlib.

        Returns
            `matplotlib` figure, None if some error occurred.
        """
        scf_cycle = abiinspect.GroundStateScfCycle.from_file(self.output_file.path)
        if scf_cycle is not None:
            if "title" not in kwargs: kwargs["title"] = str(self)
            return scf_cycle.plot(**kwargs)

    def get_results(self, **kwargs):
        results = super(ScfTask, self).get_results(**kwargs)

        # Open the GSR file and add its data to results.out
        with self.open_gsr() as gsr:
            results["out"].update(gsr.as_dict())
            # Add files to GridFS
            results.add_gridfs_files(GSR=gsr.filepath)

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
        return self._restart()

    def get_results(self, **kwargs):
        results = super(NscfTask, self).get_results(**kwargs)

        # Read the GSR file.
        with  self.open_gsr() as gsr:
            results["out"].update(gsr.as_dict())
            # Add files to GridFS
            results.add_gridfs_files(GSR=gsr.filepath)

        return results


class RelaxTask(AbinitTask, ProduceGsr):
    """
    Task for structural optimizations.
    """
    # What about a possible ScfConvergenceWarning?
    CRITICAL_EVENTS = [
        events.RelaxConvergenceWarning,
    ]

    def change_structure(self, structure):
        """Change the input structure."""
        print("changing structure")
        print("old:\n" + str(self.strategy.abinit_input.structure) + "\n")
        print("new:\n" + str(structure) + "\n")
        self.strategy.abinit_input.set_structure(structure)

    def read_final_structure(self):
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
                infile = self.out_to_in(ofile)
                break
        else:
            raise self.RestartError("Cannot find the WFK|DEN file to restart from.")

        # Read the relaxed structure from the GSR file.
        structure = self.read_final_structure()
                                                           
        # Change the structure.
        self.change_structure(structure)

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
            if "title" not in kwargs: kwargs["title"] = str(self)
            return relaxation.plot(**kwargs)

    def get_results(self, **kwargs):
        results = super(RelaxTask, self).get_results(**kwargs)

        # Open the GSR file and add its data to results.out
        with self.open_gsr() as gsr:
            results["out"].update(gsr.as_dict())
            # Add files to GridFS
            results.add_gridfs_files(GSR=gsr.filepath)

        return results


class DdeTask(AbinitTask, ProduceDdb):
    """Task for DDE calculations."""

    def get_results(self, **kwargs):
        results = super(DdeTask, self).get_results(**kwargs)
        return results.add_gridfs_file(DDB=(self.outdir.has_abiext("DDE"), "t"))


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
        return results.add_gridfs_file(DDK=(self.outdir.has_abiext("DDK"), "t"))


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
        return results.add_gridfs_file(DDB=(self.outdir.has_abiext("DDB"), "t"))

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
    Tasks for SIGMA calculations employing the self-consistent G approximation 
    Provide support for in-place restart via QPS files
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
        return self._restart()

    #def inspect(self, **kwargs):
    #    """Plot graph showing the number of k-points computed and the wall-time used"""

    def open_sigres(self):
        """
        Open the SIGRES file located in the in self.outdir. 
        Returns SigresFile object, None if file could not be found or file is not readable.
        """
        sigres_path = self.outdir.has_abiext("SIGRES")

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

    def get_results(self, **kwargs):
        results = super(SigmaTask, self).get_results(**kwargs)

        # Open the SIGRES file and add its data to results.out
        with self.open_sigres() as sigres:
            #results["out"].update(sigres.as_dict())
            results.add_gridfs_files(SIGRES=sigres.filepath)

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
        return self._restart()

    #def inspect(self, **kwargs):
    #    """
    #    Plot the Haydock iterations with matplotlib.
    #    Returns
    #        `matplotlib` figure, None if some error occurred.
    #    """
    #    haydock_cycle = abiinspect.HaydockIterations.from_file(self.output_file.path)
    #    if haydock_cycle is not None:
    #        if "title" not in kwargs: kwargs["title"] = str(self)
    #        return haydock_cycle.plot(**kwargs)

    def open_mdf(self):
        """
        Open the MDF file located in the in self.outdir.
        Returns `MdfFile` object, None if file could not be found or file is not readable.
        """
        mdf_path = self.outdir.has_abiext("MDF.nc")
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
            results.add_gridfs_files(MDF=mdf.filepath)

        return results


class OpticTask(Task):
    """
    Task for the computation of optical spectra with optic i.e.
    RPA without local-field effects and velocity operator computed from DDK files.
    """
    def __init__(self, optic_input, nscf_node, ddk_nodes, workdir=None, manager=None):
        """
        Create an instance of `OpticTask` from an string containing the input.
    
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
        app(self.input_file.path)                 # Path to the input file
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
