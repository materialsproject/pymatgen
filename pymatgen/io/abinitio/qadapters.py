# coding: utf-8
"""
Part of this code is based on a similar implementation present in FireWorks (https://pypi.python.org/pypi/FireWorks).
Work done by D. Waroquiers, A. Jain, and M. Kocher.

The main difference wrt the Fireworks implementation is that the QueueAdapter
objects provide a programmatic interface for setting important attributes 
such as the number of MPI nodes, the number of OMP threads and the memory requirements.
This programmatic interface is used by the `TaskManager` for optimizing the parameters
of the run before submitting the job (Abinit provides the autoparal option that 
allows one to get a list of parallel configuration and their expected efficiency).
"""
from __future__ import print_function, division, unicode_literals

import os
import abc
import string
import copy
import getpass
import warnings
import six
import json

from collections import namedtuple
from subprocess import Popen, PIPE, check_output
from atomicfile import AtomicFile
from monty.string import is_string
from monty.collections import AttrDict, MongoDict
from monty.functools import lazy_property
from monty.io import FileLock
from monty.subprocess import Command
from pymatgen.core.units import Time, Memory
from .utils import Condition
from .launcher import ScriptEditor

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "parse_slurm_timestr",
    "MpiRunner",
    "make_qadapter",
]


def parse_slurm_timestr(s):
    """
    A slurm time parser. Accepts a string in one the following forms:

        # "days-hours",
        # "days-hours:minutes",
        # "days-hours:minutes:seconds".
        # "minutes",
        # "minutes:seconds",
        # "hours:minutes:seconds",

    Returns:
        Time in seconds.
    Raises:
        ValueError if string is not valid.
    """
    days, hours, minutes, seconds = 0, 0, 0, 0
    if '-' in s:
        # "days-hours",
        # "days-hours:minutes",                                        
        # "days-hours:minutes:seconds".                                
        days, s = s.split("-")
        days = int(days)

        if ':' not in s:
            hours = int(float(s))
        elif s.count(':') == 1:
            hours, minutes = map(int, s.split(':'))
        elif s.count(':') == 2:
            hours, minutes, seconds = map(int, s.split(':'))
        else:
            raise ValueError("More that 2 ':' in string!")

    else:
        # "minutes",
        # "minutes:seconds",
        # "hours:minutes:seconds",
        if ':' not in s:
            minutes = int(float(s))
        elif s.count(':') == 1:
            minutes, seconds = map(int, s.split(':'))
        elif s.count(':') == 2:
            hours, minutes, seconds = map(int, s.split(':'))
        else:
            raise ValueError("More than 2 ':' in string!")

    return Time((days*24 + hours)*3600 + minutes*60 + seconds, "s")


def time2slurm(timeval, unit="s"):
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the slurm convention: "days-hours:minutes:seconds".

    >>> assert time2slurm(61) == '0-0:1:1' and time2slurm(60*60+1) == '0-1:0:1'
    >>> assert time2slurm(0.5, unit="h") == '0-0:30:0'
    """
    d, h, m, s = 24*3600, 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    days, hours = divmod(timeval, d)
    hours, minutes = divmod(hours, h)
    minutes, secs = divmod(minutes, m)

    return "%d-%d:%d:%d" % (days, hours, minutes, secs)


def time2pbspro(timeval, unit="s"):
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the PbsPro convention: "hours:minutes:seconds".

    >>> assert time2pbspro(2, unit="d") == '48:0:0' 
    """
    h, m, s = 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    hours, minutes = divmod(timeval, h)
    minutes, secs = divmod(minutes, m)

    return "%d:%d:%d" % (hours, minutes, secs)


def timelimit_parser(s):
    """Convert a float or a string into time in seconds."""
    try:
        return Time(float(s), "s")
    except ValueError:
        return parse_slurm_timestr(s)


class MpiRunner(object):
    """
    This object provides an abstraction for the mpirunner provided 
    by the different MPI libraries. It's main task is handling the
    different syntax and options supported by the different mpirunners.
    """
    def __init__(self, name, type=None, options=""):
        self.name = name
        self.type = None
        self.options = options

    def string_to_run(self, executable, mpi_procs, stdin=None, stdout=None, stderr=None):
        stdin = "< " + stdin if stdin is not None else ""
        stdout = "> " + stdout if stdout is not None else ""
        stderr = "2> " + stderr if stderr is not None else ""

        if self.has_mpirun:
            if self.type is None:
                # TODO: better treatment of mpirun syntax.
                #se.add_line('$MPIRUN -n $MPI_PROCS $EXECUTABLE < $STDIN > $STDOUT 2> $STDERR')
                num_opt = "-n " + str(mpi_procs)
                cmd = " ".join([self.name, num_opt, executable, stdin, stdout, stderr])
            else:
                raise NotImplementedError("type %s is not supported!")
        else:
            #assert mpi_procs == 1
            cmd = " ".join([executable, stdin, stdout, stderr])

        return cmd

    @property
    def has_mpirun(self):
        """True if we are running via mpirun, mpiexec ..."""
        return self.name is not None


class QHardwareError(Exception):
    """Exceptions raised by QHardware."""


class DistribError(QHardwareError):
    """Exception raised by distribute."""


class QHardware(object):
    """
    This object collects information on a partition (a la slurm)
    Partitions can be thought of as a set of resources and parameters around their use.
    A partition has a ``QueueAdapter``.

    Basic definition::

        * A node refers to the physical box, i.e. cpu sockets with north/south switches connecting memory systems 
          and extension cards, e.g. disks, nics, and accelerators

        * A cpu socket is the connector to these systems and the cpu cores

        * A cpu core is an independent computing with its own computing pipeline, logical units, and memory controller. 
          Each cpu core will be able to service a number of cpu threads, each having an independent instruction stream 
          but sharing the cores memory controller and other logical units.
    """
    class Entry(object):
        def __init__(self, type, default=None, mandatory=False, parser=None, help="No help available"):
            self.type, self.default, self.parser, self.mandatory = type, default, parser, mandatory
            if callable(default): self.default = default()

        def eval(self, value):
            if self.type is not object: value = self.type(value)
            if self.parser is not None: value = self.parser(value)
            return value
                
    ENTRIES = dict(
        # mandatory
        #qname=Entry(type=str, mandatory=True, help="Name of the partition"),
        num_nodes=Entry(type=int, mandatory=True, help="Number of nodes"),
        sockets_per_node=Entry(type=int, mandatory=True, help="Number of sockets per node"),
        cores_per_socket=Entry(type=int, mandatory=True, help="Number of cores per node"),
        mem_per_node=Entry(type=str, mandatory=True, help="Memory per node", parser=Memory.from_string),
        #min_cores=Entry(type=int, mandatory=True, help="Minimum number of cores that can be used"),
        #max_cores=Entry(type=int, mandatory=True, help="Maximum number of nodes that can be used"),
        #timelimit=Entry(type=str, mandatory=True, help="Time limit", parser=timelimit_parser),
        #priority=Entry(type=int, mandatory=True, help="Priority level, integer number > 0"),
        # optional
        #allocate_nodes=Entry(type=bool, default=False, help="True if we must allocate entire nodes"),
        #condition=Entry(type=object, default=Condition, help="Condition object (dictionary)", parser=Condition),
    )
    del Entry

    def __init__(self, **kwargs):
        """The possible arguments are documented in Partition.ENTRIES."""
        for key, entry in self.ENTRIES.items():
            try:
                value = entry.eval(kwargs.pop(key)) #; print(key, value)
                setattr(self, key, value)
            except KeyError:
                if entry.mandatory: raise ValueError("key %s must be specified" % key)
                setattr(self, key, entry.default)

        if kwargs:
            raise ValueError("Found invalid keywords in the partition section:\n %s" % str(list(kwargs.keys())))

        # Convert memory to megabytes.
        self.mem_per_node = float(self.mem_per_node.to("Mb"))
        if self.mem_per_node <= 0:
            raise ValueError("mem_per_node %s <= 0" % self.mem_per_node)

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        app("   num_nodes: %d, sockets_per_node: %d, cores_per_socket: %d, mem_per_node %s," % 
            (self.num_nodes, self.sockets_per_node, self.cores_per_socket, self.mem_per_node))
        return "\n".join(lines)

    @property
    def num_cores(self):
        """Total number of cores available in the partition."""
        return self.cores_per_socket * self.sockets_per_node * self.num_nodes

    @property
    def cores_per_node(self):
        """Number of cores per node."""
        return self.cores_per_socket * self.sockets_per_node

    @property
    def mem_per_core(self):
        """Memory available on a single node."""
        return self.mem_per_node / self.cores_per_node

    def can_use_omp_threads(self, omp_threads):
        """True if omp_threads fit in a node."""
        return self.cores_per_node >= omp_threads

    def divmod_node(self, mpi_procs, omp_threads):
        """
        Use divmod to compute (num_nodes, rest_cores)
        """
        return divmod(mpi_procs * omp_threads, self.cores_per_node)


class _ExcludeNodesFile(object):
    """
    This file contains the list of nodes to be excluded. 
    Nodes are indexed by queue name.
    """ 
    DIRPATH = os.path.join(os.getenv("HOME"), ".abinit", "abipy")
    FILEPATH = os.path.join(DIRPATH, "exclude_nodes.json")

    def __init__(self):
        if not os.path.exists(self.FILEPATH):
            if not os.path.exists(self.DIRPATH): os.makedirs(self.DIRPATH)
            with FileLock(self.FILEPATH):
                with open(self.FILEPATH, "w") as fh:
                    json.dump({}, fh)

    def read_nodes(self, qname):
        with open(self.FILEPATH, "w") as fh:
            return json.load(fh).get(qname, [])

    def add_nodes(self, qname, nodes):
        nodes = (nodes,) if not isinstance(nodes, (tuple, list)) else nodes
        with FileLock(self.FILEPATH):
            with AtomicFile(filepath, mode="w+") as fh:
                d = json.load(fh)
                if qname in d:
                    d["qname"].extend(nodes)
                    d["qname"] = list(set(d["qname"]))
                else:
                    d["qname"] = nodes
                json.dump(d, fh)

_EXCL_NODES_FILE = _ExcludeNodesFile()


def make_qadapter(**kwargs):
    """
    Return the concrete `Adapter` class from a string.
    """
    qtype = kwargs["queue"].pop("qtype")
    return {"shell": ShellAdapter,
            "slurm": SlurmAdapter,
            "pbspro": PbsProAdapter,
            "torque": TorqueAdapter,
            "sge": SGEAdapter,
            "moab": MOABAdapter,
            }[qtype.lower()](**kwargs)


class QueueAdapterError(Exception):
    """Error class for exceptions raise by QueueAdapter."""

class QueueAdapterDistribError(Exception):
    """Raised by distribute."""


class QueueAdapter(six.with_metaclass(abc.ABCMeta, object)):
    """
    The QueueAdapter is responsible for all interactions with a specific
    queue management system. This includes handling all details of queue
    script format as well as queue submission and management.

    This is the Abstract base class defining the methods that 
    must be implemented by the concrete classes.
    A user should extend this class with implementations that work on
    specific queue systems.
    """
    Error = QueueAdapterError
    DistribError = QueueAdapterDistribError

    def __init__(self, **kwargs):
                 #qparams=None, setup=None, modules=None, shell_env=None, omp_env=None, 
                 #pre_run=None, post_run=None, mpi_runner=None,):
        """
        Args:
            qname:
            qparams:
                Dictionary with the paramenters used in the template.
            setup:
                String or list of commands to execute during the initial setup.
            modules:
                String or list of modules to load before running the application.
            shell_env:
                Dictionary with the environment variables to export
                before running the application.
            omp_env:
                Dictionary with the OpenMP variables.
            pre_run:
                String or list of commands to execute before launching the calculation.
            post_run:
                String or list of commands to execute once the calculation is completed.
            mpi_runner:
                Path to the MPI runner or `MpiRunner` instance. None if not used
            max_num_attempts:
                Default to 2
            qverbatim
        """
        # Make defensive copies so that we can change the values at runtime.
        kwargs = copy.deepcopy(kwargs)
        self.priority = int(kwargs.pop("priority"))

        self.hw = QHardware(**kwargs.pop("hardware"))
        self._parse_queue(kwargs.pop("queue"))
        self._parse_limits(kwargs.pop("limits"))
        self._parse_job(kwargs.pop("job"))

        if kwargs:
            raise ValueError("Found unknown keywords:\n%s" % kwargs.keys())

        self.validate()

        # List of dictionaries with the parameters used to submit jobs
        # The launcher will use this information to increase the resources
        self.attempts, self.max_num_attempts = [], kwargs.pop("max_num_attempts", 2)

        # Initialize some values from the info reported in the partition.
        self.set_mpi_procs(self.min_cores)
        self.set_mem_per_proc(self.min_mem_per_core)

        # Final consistency check.
        self.validate()

    def validate(self):
        # Parse the template so that we know the list of supported options.
        err_msg = ""
        for param in self.qparams:
            if param not in self.supported_qparams:
                err_msg += "Unsupported QUEUE parameter name %s\n" % param
                err_msg += "Supported are: \n"
                for param_sup in self.supported_qparams:
                    err_msg += "    %s \n" % param_sup
        if err_msg:
            raise ValueError(err_msg)

        errors = []
        app = errors.append
        if self.priority <= 0: app("priority must be > 0")
        if not (1 <= self.min_cores <= self.hw.num_cores >= self.max_cores):
            app("1 <= min_cores <= hardware num_cores >= max_cores not satisfied")

        if errors:
            raise ValueError("\n".join(errors))

    def _parse_limits(self, d):
        self.set_timelimit(timelimit_parser(d.pop("timelimit")))
        self.min_cores = int(d.pop("min_cores"))
        self.max_cores = int(d.pop("max_cores"))
        self.min_mem_per_core = d.pop("min_mem_per_core", self.hw.mem_per_core)
        self.max_mem_per_core = d.pop("max_mem_per_core", self.hw.mem_per_node)
        self.allocate_nodes = bool(d.pop("allocate_nodes", False))
        self.condition = Condition(d.pop("condition", {}))

        if d:
            raise ValueError("Found unknown keyword(s) in limits section:\n %s" % d.keys())

    def _parse_job(self, d):
        setup = d.pop("setup", None)
        if is_string(setup): setup = [setup]
        self.setup = setup[:] if setup is not None else []

        omp_env = d.pop("omp_env", None)
        self.omp_env = omp_env.copy() if omp_env is not None else {}

        modules = d.pop("modules", None)
        if is_string(modules): modules = [modules]
        self.modules = modules[:] if modules is not None else []

        shell_env = d.pop("shell_env", None)
        self.shell_env = shell_env.copy() if shell_env is not None else {}

        self.mpi_runner = d.pop("mpi_runner", None)
        if not isinstance(self.mpi_runner, MpiRunner):
            self.mpi_runner = MpiRunner(self.mpi_runner)

        pre_run = d.pop("pre_run", None)
        if is_string(pre_run): pre_run = [pre_run]
        self.pre_run = pre_run[:] if pre_run is not None else []

        post_run = d.pop("post_run", None)
        if is_string(post_run): post_run = [post_run]
        self.post_run = post_run[:] if post_run is not None else []

        if d:
            raise ValueError("Found unknown keyword(s) in job section:\n %s" % d.keys())

    def _parse_queue(self, d):
        # Init params
        qparams = d.pop("qparams", None)
        self._qparams = qparams.copy() if qparams is not None else {}

        self.set_qname(d.pop("qname"))
        if d:
            raise ValueError("Found unknown keyword(s) in queue section:\n %s" % d.keys())

    def __str__(self):
        lines = ["%s:%s" % (self.__class__.__name__, self.qname)]
        app = lines.append
        app("Hardware:\n" + str(self.hw))
        #lines.extend(["qparams:\n", str(self.qparams)])
        if self.has_omp: app(str(self.omp_env))

        return "\n".join(lines)

    @property
    def qparams(self):
        return self._qparams

    @lazy_property
    def supported_qparams(self):
        """
        Dictionary with the supported parameters that can be passed to the 
        queue manager (obtained by parsing QTEMPLATE).
        """ 
        import re
        return re.findall("\$\$\{(\w+)\}", self.QTEMPLATE)

    @property
    def has_mpi(self):
        return self.has_mpirun
    
    @property
    #@deprecated(has_mpi)
    def has_mpirun(self):
        """True if we are using a mpirunner"""
        return bool(self.mpi_runner)

    @property
    def has_omp(self):
        """True if we are using OpenMP threads"""
        return hasattr(self, "omp_env") and bool(getattr(self, "omp_env"))

    @property
    def num_cores(self):
        """Total number of cores employed"""
        return self.mpi_procs * self.omp_threads 

    @property
    def omp_threads(self):
        """Number of OpenMP threads."""
        if self.has_omp:
            return self.omp_env["OMP_NUM_THREADS"]
        else:
            return 1

    @property
    def pure_mpi(self):
        """True if only MPI is used."""
        return self.has_mpi and not self.has_omp

    @property
    def pure_omp(self):
        """True if only OpenMP is used."""
        return self.has_omp and not self.has_mpi

    @property
    def hybrid_mpi_omp(self):
        """True if we are running in MPI+Openmp mode."""
        return self.has_omp and self.has_mpi

    @property
    def run_info(self):
        """String with info on the run."""
        return "MPI: %d, OMP: %d" % (self.mpi_procs, self.omp_threads)

    def deepcopy(self):
        return copy.deepcopy(self)

    def record_attempt(self, queue_id): # retcode):
        self.attempts.append(AttrDict(queue_id=queue_id, mpi_procs=self.mpi_procs, omp_threads=self.omp_threads,
                                      mem_per_proc=self.mem_per_proc, timelimit=self.timelimit)) # retcode=retcode)) 
        return len(self.attempts)

    def remove_attempt(self, index):
        self.attempts.pop(index)

    @property
    def num_attempts(self):
        return len(self.attempts)

    @property
    def last_attempt(self):
        if len(self.attempts) > 0:
            return self.attempts[-1]
        else:
            return None

    def set_omp_threads(self, omp_threads):
        """Set the number of OpenMP threads."""
        #if not self.max_cores >= self.mpi_procs * omp_threads >= self.min_cores:
        #    print(self.max_cores, self.mpi_procs, omp_threads, self.min_cores)
        #    raise self.Error("self.max_cores >= mpi_procs * omp_threads >= self.min_cores not satisfied")
        self.omp_env["OMP_NUM_THREADS"] = omp_threads

    @property
    def mpi_procs(self):
        """Number of CPUs used for MPI."""
        return self._mpi_procs

    def set_mpi_procs(self, mpi_procs):
        """Set the number of MPI processes."""
        #if not self.max_cores >= mpi_procs * self.omp_threads >= self.min_cores:
        #    print(self.max_cores, mpi_procs, self.omp_threads, self.min_cores)
        #    raise self.Error("self.max_cores >= mpi_procs * omp_threads >= self.min_cores not satisfied")
        self._mpi_procs = mpi_procs

    @property
    def qname(self):
        return self._qname

    def set_qname(self, qname):
        self._qname = qname

    @property
    def timelimit(self):
        """Returns the walltime in seconds."""
        return self._timelimit

    def set_timelimit(self, timelimit):
        """Set the walltime in seconds."""
        self._timelimit = timelimit

    @property
    def mem_per_proc(self):
        """The memory per process in Megabytes."""
        return self._mem_per_proc
                                                
    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in Megabytes"""
        self._mem_per_proc = mem_mb

    @property
    def total_mem(self):
        """Total memory required by the job in Megabytes."""
        return Memory(self.mem_per_proc * self.mpi_procs, "Mb")

    @abc.abstractmethod
    def cancel(self, job_id):
        """
        Cancel the job. 

        Args:
            job_id:
                (in) Job identifier.

        Returns:
            Exit status.
        """

    def can_run_pconf(self, pconf):
        """True if the qadapter in principle is able to run the ``ParalConf`` pconf"""
        if not self.max_cores >= pconf.num_cores >= self.min_cores: return False
        if not self.hw.can_use_omp_threads(self.omp_threads): return False
        if pconf.mem_per_proc > self.hw.mem_per_node: return False

        try:
            self.distribute(pconf.mpi_procs, pconf.omp_threads, pconf.mem_per_proc)
        except self.DistribError:
            return False

        return self.condition(pconf)

    def distribute(self, mpi_procs, omp_threads, mem_per_proc):
        """
        Returns (num_nodes, mpi_per_node)

        Aggressive: When Open MPI thinks that it is in an exactly- or under-subscribed mode
        (i.e., the number of running processes is equal to or less than the number of available processors),
        MPI processes will automatically run in aggressive mode, meaning that they will never voluntarily give
        up the processor to other processes. With some network transports, this means that Open MPI will spin
        in tight loops attempting to make message passing progress, effectively causing other processes to not get
        any CPU cycles (and therefore never make any progress)
        """
        Distrib = namedtuple("Distrib", "num_nodes mpi_per_node exact") # mem_per_node
        hw = self.hw

        # TODO: Add check on user-memory
        if mem_per_proc < hw.mem_per_node:
            # Try to use all then cores in the node.
            num_nodes, rest_cores = hw.divmod_node(mpi_procs, omp_threads)
            if num_nodes == 0 and mpi_procs * mem_per_proc <= hw.mem_per_node:
                return Distrib(num_nodes=1, mpi_per_node=mpi_procs, exact=True)

            mpi_per_node = mpi_procs // num_nodes
            if rest_cores == 0 and mpi_per_node * mem_per_proc <= hw.mem_per_node:
                return Distrib(num_nodes=num_nodes, mpi_per_node=mpi_per_node, exact=True)

        if mem_per_proc <= 0:
            logger.warning("mem_per_proc <= 0")
            mem_per_proc =  hw.mem_per_core

        # Try first to pack MPI processors in a node as much as possible
        mpi_per_node = int(hw.mem_per_node / mem_per_proc)
        if mpi_per_node == 0:
            raise self.DistribError(
                "mem_pre_proc > mem_per_node.\n Cannot distribute mpi_procs %d, omp_threads %d, mem_per_proc %s" %
                             (mpi_procs, omp_threads, mem_per_proc))

        num_nodes = (mpi_procs * omp_threads) // mpi_per_node
        #print(mpi_per_node, num_nodes)

        if mpi_per_node * omp_threads <= hw.cores_per_node and mem_per_proc <= hw.mem_per_node:
            return Distrib(num_nodes=num_nodes, mpi_per_node=mpi_per_node, exact=False)

        if (mpi_procs * omp_threads) % mpi_per_node != 0:
            # Have to reduce the number of MPI procs per node
            for mpi_per_node in reversed(range(1, mpi_per_node)):
                if mpi_per_node > hw.cores_per_node: continue
                num_nodes = (mpi_procs * omp_threads) // mpi_per_node
                if (mpi_procs * omp_threads) % mpi_per_node == 0 and mpi_per_node * mem_per_proc <= hw.mem_per_node:
                    return Distrib(num_nodes=num_nodes, mpi_per_node=mpi_per_node, exact=False)
        else:
            raise self.DistribError("Cannot distribute mpi_procs %d, omp_threads %d, mem_per_proc %s" %
                            (mpi_procs, omp_threads, mem_per_proc))

    def optimize_params(self):
        logger.debug("optimize_params of baseclass --> no optimization available!!!")
        return {}

    def get_subs_dict(self):
        """
        Return substitution dict for replacements into the template
        Subclasses may want to customize this method.
        """ 
        d = self.qparams.copy()
        d.update(self.optimize_params())
        # clean null values
        return {k: v for k, v in d.items() if v is not None}

    def _make_qheader(self, job_name, qout_path, qerr_path):
        """Return a string with the options that are passed to the resource manager."""
        # get substitution dict for replacements into the template 
        subs_dict = self.get_subs_dict()

        # Set job_name and the names for the stderr and stdout of the 
        # queue manager (note the use of the extensions .qout and .qerr
        # so that we can easily locate this file.
        subs_dict['job_name'] = job_name.replace('/', '_')
        subs_dict['_qout_path'] = qout_path
        subs_dict['_qerr_path'] = qerr_path

        qtemplate = QScriptTemplate(self.QTEMPLATE)
        # might contain unused parameters as leftover $$.
        unclean_template = qtemplate.safe_substitute(subs_dict)  

        # Remove lines with leftover $$.
        clean_template = []
        for line in unclean_template.split('\n'):
            if '$$' not in line:
                clean_template.append(line)

        return '\n'.join(clean_template)

    def get_script_str(self, job_name, launch_dir, executable, qout_path, qerr_path,
                       stdin=None, stdout=None, stderr=None):
        """
        Returns a (multi-line) String representing the queue script, e.g. PBS script.
        Uses the template_file along with internal parameters to create the script.

        Args:
            job_name:
                Name of the job.
            launch_dir: 
                (str) The directory the job will be launched in.
            executable:
                String with the name of the executable to be executed.
            qout_path
                Path of the Queue manager output file.
            qerr_path:
                Path of the Queue manager error file.
        """
        # PbsPro does not accept job_names longer than 15 chars.
        if len(job_name) > 14 and isinstance(self, PbsProAdapter):
            job_name = job_name[:14]

        # Construct the header for the Queue Manager.
        qheader = self._make_qheader(job_name, qout_path, qerr_path)

        # Add the bash section.
        se = ScriptEditor()

        if self.setup:
            se.add_comment("Setup section")
            se.add_lines(self.setup)
            se.add_emptyline()

        if self.modules:
            se.add_comment("Load Modules")
            se.add_line("module purge")
            se.load_modules(self.modules)
            se.add_emptyline()

        if self.has_omp:
            se.add_comment("OpenMp Environment")
            se.declare_vars(self.omp_env)
            se.add_emptyline()

        if self.shell_env:
            se.add_comment("Shell Environment")
            se.declare_vars(self.shell_env)
            se.add_emptyline()

        # Cd to launch_dir
        se.add_line("cd " + os.path.abspath(launch_dir))

        if self.pre_run:
            se.add_comment("Commands before execution")
            se.add_lines(self.pre_run)
            se.add_emptyline()

        # Construct the string to run the executable with MPI and mpi_procs.
        line = self.mpi_runner.string_to_run(executable, self.mpi_procs, 
                                             stdin=stdin, stdout=stdout, stderr=stderr)
        se.add_line(line)

        if self.post_run:
            se.add_emptyline()
            se.add_comment("Commands after execution")
            se.add_lines(self.post_run)

        return qheader + se.get_script_str() + "\n"

    def submit_to_queue(self, script_file):
        """
        Public API: wraps the concrete implementation _submit_to_queue

        Raises:
            QueueAdapterError if we have already tried to submit the job max_num_attempts
        """
        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        if self.num_attempts == self.max_num_attempts:
            raise self.Error("num_attempts %s == max_num_attempts" % (self.num_attempts, self.max_num_attempts))

        process, queue_id = self._submit_to_queue(script_file)
        self.record_attempt(queue_id)
        return process, queue_id

    @abc.abstractmethod
    def _submit_to_queue(self, script_file):
        """
        Submits the job to the queue, probably using subprocess or shutil
        This method must be provided by the concrete classes and will be called by submit_to_queue

        Args:
            script_file: 
                (str) name of the script file to use (String)
        Returns:
            process, queue_id
        """

    @abc.abstractmethod
    def get_njobs_in_queue(self, username=None):
        """
        returns the number of jobs in the queue, probably using subprocess or shutil to
        call a command like 'qstat'. returns None when the number of jobs cannot be determined.

        Args:
            username: (str) the username of the jobs to count (default is to autodetect)
        """

    # Methods to fix problems
    def add_exclude_nodes(self, nodes):
        return _EXCL_NODES_FILE.add_nodes(self.qname, nodes)

    def get_exclude_nodes(self):
        return _EXCL_NODES_FILE.read_nodes(self.qname)

    @abc.abstractmethod
    def exclude_nodes(self, nodes):
        """
        Method to exclude nodes in the calculation
        """

    def more_mem_per_proc(self, factor):
        """
        Method to increase the amount of memory asked for, by factor.
        Return True if success.
        """
        base_increase = 2000
        old_mem = self.mem_per_proc
        new_mem = old_mem + factor*base_increase

        if new_mem < self.max_mem_per_node:
            self.set_mem_per_proc(new_mem)
            return True

        logger.warning('could not increase mem_per_proc further')
        return False

    def more_mpi_procs(self, factor):
        """
        Method to increase the number of MPI procs. Return True if success.
        """
        base_increase = 12
        new_cpus = self.mpi_procs + factor * base_increase

        if new_cpus * self.omp_threads < self.max_cores:
            self.set_mpi_procs(new_cpus)
            return True

        logger.warning('more_mpi_procs reached the limit')
        return False

    def get_score(self, pconf):
        """
        Receives a ``ParalConf`` object, pconf, and returns a number that will be used
        to select the partion on the cluster on which the task will be submitted.
        Returns -inf if paral_conf cannot be executed on this partition.
        """
        minf = float("-inf")
        if not self.can_run_pconf(pconf): return minf
        return self.priority

####################
# Concrete classes #
####################


class ShellAdapter(QueueAdapter):
    QTYPE = "shell"

    QTEMPLATE = """\
#!/bin/bash
$${qverbatim}
"""

    def cancel(self, job_id):
        return os.system("kill -9 %d" % job_id)

    def _submit_to_queue(self, script_file):
        try:
            # submit the job, return process and pid.
            process = Popen(("/bin/bash", script_file), stderr=PIPE)
            return process, process.pid
        except Exception as exc:
            raise self.Error(str(exc))

    def get_njobs_in_queue(self, username=None):
        return None

    def exclude_nodes(self, nodes):
        return False


class SlurmAdapter(QueueAdapter):
    QTYPE = "slurm"

    QTEMPLATE = """\
#!/bin/bash

#SBATCH --partition=$${partition}
#SBATCH --job-name=$${job_name}
#SBATCH --ntasks=$${ntasks}
#SBATCH --ntasks-per-node=$${ntasks_per_node}
#SBATCH --cpus-per-task=$${cpus_per_task}
#SBATCH --mem=$${mem}
#SBATCH --mem-per-cpu=$${mem_per_cpu}
#SBATCH --time=$${time}
#SBATCH	--exclude=$${exclude_nodes}
#SBATCH	--nodes=$${nodes}
#SBATCH --account=$${account}
#SBATCH --mail-user=$${mail_user}
#SBATCH --mail-type=$${mail_type}
#SBATCH --constraint=$${constraint}
#SBATCH --gres=$${gres}
#SBATCH --requeue=$${requeue}
#SBATCH --nodelist=$${nodelist}
#SBATCH --propagate=$${propagate}
#SBATCH --licenses=$${licenses}
#SBATCH --output=$${_qout_path}
#SBATCH --error=$${_qerr_path}
$${qverbatim}
"""

    def set_qname(self, qname):
        super(self.__class__, self).set_qname(qname)
        self.qparams["partition"] = qname

    def set_mpi_procs(self, mpi_procs):
        """Set the number of CPUs used for MPI."""
        super(self.__class__, self).set_mpi_procs(mpi_procs)
        self.qparams["ntasks"] = mpi_procs

    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in Megabytes"""
        super(self.__class__, self).set_mem_per_proc(mem_mb)
        self.qparams["mem_per_cpu"] = int(mem_mb)
        # Remove mem if it's defined.
        self.qparams.pop("mem", None)

    def set_timelimit(self, timelimit):
        super(self.__class__, self).set_timelimit(timelimit)
        self.qparams["time"] = time2slurm(timelimit)

    def cancel(self, job_id):
        return os.system("scancel %d" % job_id)

    def _submit_to_queue(self, script_file, submit_err_file="sbatch.err"):
        submit_err_file = os.path.join(os.path.dirname(script_file), submit_err_file)

        # submit the job
        try:
            cmd = ['sbatch', script_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            # write the err output to file, a error parser may read it and a fixer may know what to do ...

            with open(submit_err_file, mode='w') as f:
                f.write('sbatch submit process stderr:')
                f.write(str(process.stderr.read()))
                f.write('qparams:')
                f.write(str(self.qparams))

            process.wait()

            # grab the returncode. SLURM returns 0 if the job was successful
            if process.returncode == 0:
                try:
                    # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
                    queue_id = int(process.stdout.read().split()[3])
                    logger.info('Job submission was successful and queue_id is {}'.format(queue_id))
                except:
                    # probably error parsing job code
                    queue_id = None
                    logger.warning('Could not parse job id following slurm...')

                finally:
                    return process, queue_id

            else:
                # some qsub error, e.g. maybe wrong queue specified, don't have permission to submit, etc...
                err_msg = ("Error in job submission with SLURM file {f} and cmd {c}\n".format(f=script_file, c=cmd) + 
                           "The error response reads:\n {c}".format(c=process.stderr.read()))
                raise self.Error(err_msg)

        except Exception as details:
            msg = 'Error while submitting job:\n' + str(details)
            logger.critical(msg)
            with open(submit_err_file, mode='a') as f:
                f.write(msg)

            try:
                print('sometimes we land here, no idea what is happening ... Michiel')
                print("details:\n", details, "cmd\n", cmd, "\nprocess.returcode:", process.returncode)
            except:
                pass

            # random error, e.g. no qsub on machine!
            raise self.Error('Running sbatch caused an error...')

    def exclude_nodes(self, nodes):
        try:
            if 'exclude_nodes' not in self.qparams:
                self.qparams.update({'exclude_nodes': 'node' + nodes[0]})
                print('excluded node %s' % nodes[0])

            for node in nodes[1:]:
                self.qparams['exclude_nodes'] += ',node' + node
                print('excluded node %s' % node)

            return True

        except (KeyError, IndexError):
            return False

    def get_njobs_in_queue(self, username=None):
        if username is None:
            username = getpass.getuser()

        cmd = ['squeue', '-o "%u"', '-u', username]
        process = Popen(cmd, shell=False, stdout=PIPE)
        process.wait()

        # parse the result
        if process.returncode == 0:
            # lines should have this form
            # username
            # count lines that include the username in it

            outs = process.stdout.readlines()
            njobs = len([line.split() for line in outs if username in line])
            logger.info('The number of jobs currently in the queue is: {}'.format(njobs))
            return njobs

        # there's a problem talking to squeue server?
        err_msg = ('Error trying to get the number of jobs in the queue using squeue service' + 
                   'The error response reads:\n {}'.format(process.stderr.read()))
        logger.critical(err_msg)

        return None

    def get_job_info(self, job_id):
        # See https://computing.llnl.gov/linux/slurm/sacct.html
        #If SLURM job ids are reset, some job numbers will        
	    #probably appear more than once refering to different jobs.
	    #Without this option only the most recent jobs will be displayed.          

        #state
        #Displays the job status, or state.
        #Output can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, 
        #PREEMPTED or NODE_FAIL. If more information is available on the job state than will fit 
        #into the current field width (for example, the uid that CANCELLED a job) the state will be followed by a "+". 
        # You can increase the size of the displayed state using the "%NUMBER" format modifier described earlier. 

        #gmatteo@master2:~
        #sacct --job 112367 --format=jobid,exitcode,state --allocations --parsable2
        #JobID|ExitCode|State
        #112367|0:0|RUNNING

        #output = check_output(["ls", "-l", "/dev/null"]
        cmd = "sacct --job %d --format=jobid,exitcode,state --allocations --parsable2" % job_id
        output = str(check_output([cmd]))

        # Parse output.
        qid, exitcode, state = output.split("|")
        qid = int(qid)
        assert qid == job_id
        if ":" in exitcode:
            exitcode, signal = map(int, exitcode.split(":"))
        else:
            exitcode, signal = int(exitcode), None

        i = state.find("+")
        if i != -1: state = state[:i]

        #class JobInfo(namedtuple("JobInfo", "queue_id exitcode signal state")):
        #    def __bool__(self):
        #        return self.state != "CannotDected"
        #    __notzero__ = __bool_
        #    def completed(self):
        #    def cancelled(self):
        #    def failed(self):
        #    def timeout(self):
        #    def node_fail(self):
        #return jobinfo()

    #def get_start_time(self, job_id)
    #    squeue  --start -j  116791           
    #      JOBID PARTITION     NAME     USER  ST           START_TIME  NODES NODELIST(REASON)
    #     116791      defq gs6q2wop cyildiri  PD  2014-11-04T09:27:15     16 (QOSResourceLimit)
    # For more info
    #login1$ scontrol show job 1676354


class PbsProAdapter(QueueAdapter):
    QTYPE = "pbspro"

#PBS -l select=$${select}:ncpus=$${ncpus}:vmem=$${vmem}mb:mpiprocs=$${mpiprocs}:ompthreads=$${ompthreads}
#PBS -l select=$${select}:ncpus=1:vmem=$${vmem}mb:mpiprocs=1:ompthreads=$${ompthreads}
####PBS -l select=$${select}:ncpus=$${ncpus}:vmem=$${vmem}mb:mpiprocs=$${mpiprocs}:ompthreads=$${ompthreads}
####PBS -l pvmem=$${pvmem}mb

    QTEMPLATE = """\
#!/bin/bash

#PBS -q $${queue}
#PBS -N $${job_name}
#PBS -A $${account}
#PBS -l select=$${select}
#PBS -l walltime=$${walltime}
#PBS -l model=$${model}
#PBS -l place=$${place}
#PBS -W group_list=$${group_list}
#PBS -M $${mail_user}
#PBS -m $${mail_type}
# Submission environment
####PBS -V
#PBS -o $${_qout_path}
#PBS -e $${_qerr_path}
$${qverbatim}
"""

    def set_qname(self, qname):
        super(self.__class__, self).set_qname(qname)
        self.qparams["queue"] = qname

    def set_timelimit(self, timelimit):
        super(self.__class__, self).set_timelimit(timelimit)
        self.qparams["walltime"] = time2pbspro(timelimit)

    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in Megabytes"""
        super(self.__class__, self).set_mem_per_proc(mem_mb)
        #self.qparams["pvmem"] = int(mem_mb)
        #self.qparams["vmem"] = int(mem_mb)

    def cancel(self, job_id):
        return os.system("qdel %d" % job_id)

    def optimize_params(self):
        return {"select": self.get_select()}

    def get_select(self, ret_dict=False):
        """
        Select is not the most intuitive command. For more info see:

            * http://www.cardiff.ac.uk/arcca/services/equipment/User-Guide/pbs.html
            * https://portal.ivec.org/docs/Supercomputers/PBS_Pro
        """
        mem_per_proc = int(self.mem_per_proc)
        hw = self.hw
        #dist = part.distribute(self.mpi_procs, self.omp_threads, mem_per_proc)

        if self.pure_mpi:
            num_nodes, rest_cores = hw.divmod_node(self.mpi_procs, self.omp_threads)

            if num_nodes == 0:
                logger.info("IN_CORE PURE MPI: %s" % self.run_info)
                chunks = 1
                ncpus = rest_cores
                mpiprocs = rest_cores
                vmem = mem_per_proc * ncpus
                ompthreads = 1

            elif rest_cores == 0:
                # Can allocate entire nodes because self.mpi_procs is divisible by cores_per_node.
                logger.info("PURE MPI run commensurate with cores_per_node %s" % self.run_info)
                chunks = num_nodes
                ncpus = hw.cores_per_node
                mpiprocs = hw.cores_per_node
                vmem = ncpus * mem_per_proc
                ompthreads = 1

            else:
                logger.info("OUT-OF-CORE PURE MPI (not commensurate with cores_per_node): %s" % self.run_info)
                chunks = self.mpi_procs
                ncpus = 1
                mpiprocs = 1
                vmem = mem_per_proc
                ompthreads = 1

        elif self.pure_omp:
            # Pure OMP run.
            logger.info("PURE OPENMP run: %s" % self.run_info)
            assert hw.can_use_omp_threads(self.omp_threads)
            chunks = 1
            ncpus = self.omp_threads
            mpiprocs = 1
            vmem = mem_per_proc
            ompthreads = self.omp_threads

        elif self.hybrid_mpi_omp:
            assert hw.can_use_omp_threads(self.omp_threads)
            num_nodes, rest_cores = hw.divmod_node(self.mpi_procs, self.omp_threads)
            #print(num_nodes, rest_cores)
            # TODO: test this

            if rest_cores == 0 or num_nodes == 0:  
                logger.info("HYBRID MPI-OPENMP run, perfectly divisible among nodes: %s" % self.run_info)
                chunks = max(num_nodes, 1)
                mpiprocs = self.mpi_procs // chunks

                chunks = chunks
                ncpus = mpiprocs * self.omp_threads
                mpiprocs = mpiprocs
                vmem = mpiprocs * mem_per_proc 
                ompthreads = self.omp_threads

            else:
                logger.info("HYBRID MPI-OPENMP, NOT commensurate with nodes: %s" % self.run_info)
                chunks=self.mpi_procs
                ncpus=self.omp_threads
                mpiprocs=1
                vmem= mem_per_proc
                ompthreads=self.omp_threads

        else:
            raise RuntimeError("You should not be here")

        select_params = AttrDict(chunks=chunks, ncpus=ncpus, mpiprocs=mpiprocs, vmem=int(vmem), ompthreads=ompthreads)

        if not self.has_omp:
            s = "{chunks}:ncpus={ncpus}:vmem={vmem}mb:mpiprocs={mpiprocs}".format(**select_params)
        else:
            s = "{chunks}:ncpus={ncpus}:vmem={vmem}mb:mpiprocs={mpiprocs}:ompthreads={ompthreads}".format(**select_params)
                                                                                                            
        if ret_dict:
            return s, select_params
        return s

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        # submit the job
        try:
            cmd = ['qsub', script_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            process.wait()

            # grab the return code. PBS returns 0 if the job was successful
            if process.returncode == 0:
                try:
                    # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
                    queue_id = int(process.stdout.read().split('.')[0])
                    logger.info('Job submission was successful and queue_id is {}'.format(queue_id))

                except:
                    # probably error parsing job code
                    logger.warning("Could not parse job id following qsub...")
                    queue_id = None

                finally:
                    return process, queue_id

            else:
                # some qsub error, e.g. maybe wrong queue specified, don't have permission to submit, etc...
                msg = ('Error in job submission with PBS file {f} and cmd {c}\n'.format(f=script_file, c=cmd) + 
                       'The error response reads:\n {}'.format(process.stderr.read()))
                raise self.Error(msg)

        except Exception as exc:
            # random error, e.g. no qsub on machine!
            raise self.Error("Running qsub caused an error...\n%s" % str(exc))

    def get_njobs_in_queue(self, username=None):
        # Initialize username
        if username is None:
            username = getpass.getuser()

        # run qstat
        try:
            qstat = Command(['qstat', '-a', '-u', username]).run(timeout=5)

            # parse the result
            if qstat.retcode == 0:
                # lines should have this form
                # '1339044.sdb          username  queuename    2012-02-29-16-43  20460   --   --    --  00:20 C 00:09'
                # count lines that include the username in it

                # TODO: only count running or queued jobs. or rather, *don't* count jobs that are 'C'.
                outs = qstat.output.split('\n')
                njobs = len([line.split() for line in outs if username in line])
                logger.info('The number of jobs currently in the queue is: {}'.format(njobs))
                return njobs
        except:
            # there's a problem talking to qstat server?
            print(qstat.output.split('\n'))
            err_msg = ('Error trying to get the number of jobs in the queue using qstat service\n' +
                       'The error response reads:\n {}'.format(qstat.error))
            logger.critical(err_msg)
            return None

    def exclude_nodes(self, nodes):
        logger.warning('exluding nodes, not implemented yet in pbs')
        return False


class TorqueAdapter(PbsProAdapter):
    """Adapter for Torque."""

    QTYPE = "torque"

    QTEMPLATE = """\
#!/bin/bash

#PBS -q $${queue}
#PBS -N $${job_name}
#PBS -A $${account}
#PBS -l pmem=$${pmem}mb
####PBS -l mppwidth=$${mppwidth}
#PBS -l nodes=$${nodes}:ppn=$${ppn} 
#PBS -l walltime=$${walltime}
#PBS -l model=$${model}
#PBS -l place=$${place}
#PBS -W group_list=$${group_list}
#PBS -M $${mail_user}
#PBS -m $${mail_type}
# Submission environment
#PBS -V
#PBS -o $${_qout_path}
#PBS -e $${_qerr_path}
$${qverbatim}
"""
    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in Megabytes"""
        QueueAdapter.set_mem_per_proc(self, mem_mb)
        self.qparams["pmem"] = mem_mb
        self.qparams["mem"] = mem_mb

    @property
    def mpi_procs(self):
        """Number of MPI processes."""
        return self.qparams.get("nodes", 1)*self.qparams.get("ppn", 1)

    def set_mpi_procs(self, mpi_procs):
        """Set the number of CPUs used for MPI."""
        QueueAdapter.set_mpi_procs(mpi_procs)
        self.qparams["nodes"] = 1
        self.qparams["ppn"] = mpi_procs


class SGEAdapter(QueueAdapter):
    """
    Adapter for Sun Grid Engine (SGE) task submission software.
    """
    QTYPE = "sge"

    QTEMPLATE = """\
#!/bin/bash

#$ -A $${account}
#$ -N $${job_name}
#$ -l h rt=$${walltime}
#$ -pe $${queue} $${ncpus}
#$ -cwd
#$ -j y
#$ -m n
#$ -e $${_qerr_path}
#$ -o $${_qout_path}
#$ -S /bin/bash
$${qverbatim}
"""
    def set_mpi_procs(self, mpi_procs):
        """Set the number of CPUs used for MPI."""
        super(self.__class__, self).set_mpi_procs(mpi_procs)
        self.qparams["ncpus"] = mpi_procs

    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in Megabytes"""
        super(self.__class__, self).set_mem_per_proc(mem_mb)
        # TODO
        #raise NotImplementedError("")
        #self.qparams["mem_per_cpu"] = mem_mb
        ## Remove mem if it's defined.
        #self.qparams.pop("mem", None)

    def cancel(self, job_id):
        return os.system("qdel %d" % job_id)

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        # submit the job
        try:
            cmd = ['qsub', script_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            process.wait()

            # grab the returncode. SGE returns 0 if the job was successful
            if process.returncode == 0:
                try:
                    # output should of the form 
                    # Your job 1659048 ("NAME_OF_JOB") has been submitted 
                    queue_id = int(process.stdout.read().split(' ')[2])
                    logger.info('Job submission was successful and queue_id is {}'.format(queue_id))

                except:
                    # probably error parsing job code
                    logger.warning("Could not parse job id following qsub...")
                    queue_id = None

                finally:
                    return process, queue_id

            else:
                # some qsub error, e.g. maybe wrong queue specified, don't have permission to submit, etc...
                msg = ('Error in job submission with PBS file {f} and cmd {c}\n'.format(f=script_file, c=cmd) + 
                       'The error response reads:\n {}'.format(process.stderr.read()))
                raise self.Error(msg)

        except:
            # random error, e.g. no qsub on machine!
            raise self.Error("Running qsub caused an error...")

    def get_njobs_in_queue(self, username=None):
        # Initialize username
        if username is None:
            username = getpass.getuser()

        # run qstat
        qstat = Command(['qstat', '-u', username]).run(timeout=5)

        # parse the result
        if qstat.retcode == 0:
            # lines should contain username
            # count lines that include the username in it

            # TODO: only count running or queued jobs. or rather, *don't* count jobs that are 'C'.
            outs = qstat.output.split('\n')
            njobs = len([line.split() for line in outs if username in line])
            logger.info('The number of jobs currently in the queue is: {}'.format(njobs))

            return njobs

        # there's a problem talking to qstat server?
        err_msg = ('Error trying to get the number of jobs in the queue using qstat service\n' + 
                   'The error response reads:\n {}'.format(qstat.error))
        logger.critical(err_msg)

        return None

    def exclude_nodes(self, nodes):
        """Method to exclude nodes in the calculation"""
        logger.warning('exluding nodes, not implemented yet in SGE')
        return False


class MOABAdapter(QueueAdapter):
    """https://computing.llnl.gov/tutorials/moab/"""
    QTYPE = "moab"

    QTEMPLATE = """\
#!/bin/bash

#MSUB -a $${eligible_date}
#MSUB -A $${account}
#MSUB -c $${checkpoint_interval}
#MSUB -l feature=$${feature}
#MSUB -l gres=$${gres}
#MSUB -l nodes=$${nodes}
#MSUB -l partition=$${partition}
#MSUB -l procs=$${procs}
#MSUB -l ttc=$${ttc}
#MSUB -l walltime=$${walltime}
#MSUB -l $${resources}
#MSUB -p $${priority}
#MSUB -q $${queue}
#MSUB -S $${shell}
#MSUB -N $${job_name}
#MSUB -v $${variable_list}

#MSUB -o $${_qout_path}
#MSUB -e $${_qerr_path}
$${qverbatim}
"""
    def set_mpi_procs(self, mpi_procs):
        """Set the number of CPUs used for MPI."""
        super(self.__class__, self).set_mpi_procs(mpi_procs)
        self.qparams["procs"] = mpi_procs

    def set_timelimit(self, timelimit):
        super(self.__class__, self).set_timelimit(timelimit)
        self.qparams["walltime"] = time2pbspro(timelimit)

    def set_mem_per_proc(self, mem_mb):
        super(self.__class__, self).set_mem_per_proc(mem_mb)
        #TODO
        #raise NotImplementedError("set_mem_per_cpu")

    def exclude_nodes(self, nodes):
        logger.warning('exluding nodes, not implemented yet in MOAB')
        return False

    def cancel(self, job_id):
        return os.system("canceljob %d" % job_id)

    def _submit_to_queue(self, script_file, submit_err_file="sbatch.err"):
        """Submit a job script to the queue."""
        submit_err_file = os.path.join(os.path.dirname(script_file), submit_err_file)

        # submit the job
        try:
            cmd = ['msub', script_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            # write the err output to file, a error parser may read it and a fixer may know what to do ...

            with open(submit_err_file, mode='w') as f:
                f.write('msub submit process stderr:')
                f.write(str(process.stderr.read()))
                f.write('qparams:')
                f.write(str(self.qparams))

            process.wait()

            # grab the returncode. MOAB returns 0 if the job was successful
            if process.returncode == 0:
                try:
                    # output should be the queue_id
                    queue_id = int(process.stdout.read().split()[0])
                    logger.info('Job submission was successful and queue_id is {}'.format(queue_id))
                except:
                    # probably error parsing job code
                    queue_id = None
                    logger.warning('Could not parse job id following msub...')

                finally:
                    return process, queue_id

            else:
                # some qsub error, e.g. maybe wrong queue specified, don't have permission to submit, etc...
                err_msg = ("Error in job submission with MOAB file {f} and cmd {c}\n".format(f=script_file, c=cmd) + 
                           "The error response reads:\n {c}".format(c=process.stderr.read()))
                raise self.Error(err_msg)

        except Exception as details:
            msg = 'Error while submitting job:\n' + str(details)
            logger.critical(msg)
            with open(submit_err_file, mode='a') as f:
                f.write(msg)

            try:
                print('sometimes we land here, no idea what is happening ... Michiel')
                print("details:\n", details, "cmd\n", cmd, "\nprocess.returcode:", process.returncode)
            except:
                pass

            # random error, e.g. no qsub on machine!
            raise self.Error('Running msub caused an error...')

    def get_njobs_in_queue(self, username=None):
        if username is None:
            username = getpass.getuser()

        cmd = ['showq', '-s -u', username]
        process = Popen(cmd, shell=False, stdout=PIPE)
        process.wait()

        # parse the result
        if process.returncode == 0:
            # lines should have this form:
            ## 
            ## active jobs: N  eligible jobs: M  blocked jobs: P
            ##
            ## Total job:  1
            ##
            # Split the output string and return the last element.

            outs = process.stdout.readlines()
            njobs = int(outs.split()[-1])
            logger.info('The number of jobs currently in the queue is: {}'.format(njobs))
            return njobs

        # there's a problem talking to squeue server?
        err_msg = ('Error trying to get the number of jobs in the queue using showq service' + 
                   'The error response reads:\n {}'.format(process.stderr.read()))
        logger.critical(err_msg)

        return None


class QScriptTemplate(string.Template):
    delimiter = '$$'


