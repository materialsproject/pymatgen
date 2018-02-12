# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
The initial version of this module was based on a similar implementation
present in FireWorks (https://pypi.python.org/pypi/FireWorks).
Work done by D. Waroquiers, A. Jain, and M. Kocher.

The main difference wrt the Fireworks implementation is that the QueueAdapter
objects provide a programmatic interface for setting important attributes
such as the number of MPI nodes, the number of OMP threads and the memory requirements.
This programmatic interface is used by the `TaskManager` for optimizing the parameters
of the run before submitting the job (Abinit provides the autoparal option that
allows one to get a list of parallel configuration and their expected efficiency).
"""
from __future__ import print_function, division, unicode_literals

import sys
import os
import abc
import string
import copy
import getpass
import six
import json
import math
from . import qutils as qu

from collections import namedtuple
from subprocess import Popen, PIPE
from pymatgen.util.io_utils import AtomicFile
from monty.string import is_string, list_strings
from monty.collections import AttrDict
from monty.functools import lazy_property
from monty.inspect import all_subclasses
from monty.io import FileLock
from monty.json import MSONable
from pymatgen.core.units import Memory
from .utils import Condition
from .launcher import ScriptEditor
from .qjobs import QueueJob

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "make_qadapter",
]


__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


class SubmitResults(namedtuple("SubmitResult", "qid, out, err, process")):
    """
    named tuple createc by the concrete implementation of _submit_to_que to pass the results of the process of
    submitting the jobfile to the que.
    qid: queue id of the submission
    out: stdout of the submission
    err: stdrr of the submisison
    process: process object of the submission
    """


class MpiRunner(object):
    """
    This object provides an abstraction for the mpirunner provided
    by the different MPI libraries. It's main task is handling the
    different syntax and options supported by the different mpirunners.
    """
    def __init__(self, name, type=None, options=""):
        """
        Args:
            name (str): Name of the mpirunner e.g. mpirun, mpiexec, srun ...
            type: Type of the mpirunner (not used at present)
            options (str): String with options passed to the mpi runner e.g. "--bind-to None"
        """
        self.name = name if name else ""
        self.type = None
        self.options = str(options)

    def string_to_run(self, qad, executable, stdin=None, stdout=None, stderr=None, exec_args=None):
        """
        Build and return a string with the command required to launch `executable` with the qadapter `qad`.

        Args
            qad: Qadapter instance.
            executable (str): Executable name or path
            stdin (str): Name of the file to be used as standard input. None means no redirection.
            stdout (str): Name of the file to be used as standard output. None means no redirection.
            stderr (str): Name of the file to be used as standard error. None means no redirection.
            exec_args: Optional list of strings with options passed to `executable`.

        Return:
            String with command to execute.
        """
        stdin = "< " + stdin if stdin is not None else ""
        stdout = "> " + stdout if stdout is not None else ""
        stderr = "2> " + stderr if stderr is not None else ""

        if exec_args:
            executable = executable + " " + " ".join(list_strings(exec_args))

        basename = os.path.basename(self.name)
        if basename in ["mpirun", "mpiexec", "srun"]:
            if self.type is None:
                # $MPIRUN -n $MPI_PROCS $EXECUTABLE < $STDIN > $STDOUT 2> $STDERR
                num_opt = "-n " + str(qad.mpi_procs)
                cmd = " ".join([self.name, self.options, num_opt, executable, stdin, stdout, stderr])
            else:
                raise NotImplementedError("type %s is not supported!" % self.type)

        elif basename == "runjob":
            #runjob --ranks-per-node 2 --exp-env OMP_NUM_THREADS --exe $ABINIT < $STDIN > $STDOUT 2> $STDERR
            #runjob -n 2 --exp-env=OMP_NUM_THREADS --exe $ABINIT < $STDIN > $STDOUT 2> $STDERR
            # exe must be absolute path or relative to cwd.
            bg_size, rpn = qad.bgsize_rankspernode()
            #num_opt = "-n " + str(qad.mpi_procs)
            num_opt = "--ranks-per-node " + str(rpn)
            cmd = " ".join([self.name, self.options, num_opt, "--exp-env OMP_NUM_THREADS",
                           "--exe `which " + executable + "` ", stdin, stdout, stderr])
        else:
            if qad.mpi_procs != 1:
                raise ValueError("Cannot use mpi_procs > when mpi_runner basename=%s" % basename)
            cmd = " ".join([executable, stdin, stdout, stderr])

        return cmd

    #@property
    #def has_mpirun(self):
    #    """True if we are running via mpirun, mpiexec ..."""
    #    return self.name in ("mpirun", "mpiexec", "srun", "runjob")


class OmpEnv(AttrDict):
    """
    Dictionary with the OpenMP environment variables
    see https://computing.llnl.gov/tutorials/openMP/#EnvironmentVariables
    """
    _KEYS = [
        "OMP_SCHEDULE",
        "OMP_NUM_THREADS",
        "OMP_DYNAMIC",
        "OMP_PROC_BIND",
        "OMP_NESTED",
        "OMP_STACKSIZE",
        "OMP_WAIT_POLICY",
        "OMP_MAX_ACTIVE_LEVELS",
        "OMP_THREAD_LIMIT",
        "OMP_STACKSIZE",
        "OMP_PROC_BIND",
    ]

    @classmethod
    def as_ompenv(cls, obj):
        """Convert an object into a OmpEnv"""
        if isinstance(obj, cls): return obj
        if obj is None: return cls()
        return cls(**obj)

    def __init__(self, *args, **kwargs):
        """
        Constructor method inherited from dictionary:

        >>> assert OmpEnv(OMP_NUM_THREADS=1).OMP_NUM_THREADS == 1

        To create an instance from an INI file, use:
           OmpEnv.from_file(filename)
        """
        super(OmpEnv, self).__init__(*args, **kwargs)

        err_msg = ""
        for key, value in self.items():
            self[key] = str(value)
            if key not in self._KEYS:
                err_msg += "unknown option %s\n" % key

        if err_msg:
            raise ValueError(err_msg)

    def export_str(self):
        """Return a string with the bash statements needed to setup the OMP env."""
        return "\n".join("export %s=%s" % (k, v) for k, v in self.items())


class Hardware(object):
    """
    This object collects information on the hardware available in a given queue.

    Basic definitions:

        - A node refers to the physical box, i.e. cpu sockets with north/south switches connecting memory systems
          and extension cards, e.g. disks, nics, and accelerators

        - A cpu socket is the connector to these systems and the cpu cores

        - A cpu core is an independent computing with its own computing pipeline, logical units, and memory controller.
          Each cpu core will be able to service a number of cpu threads, each having an independent instruction stream
          but sharing the cores memory controller and other logical units.
    """
    def __init__(self, **kwargs):
        self.num_nodes = int(kwargs.pop("num_nodes"))
        self.sockets_per_node = int(kwargs.pop("sockets_per_node"))
        self.cores_per_socket = int(kwargs.pop("cores_per_socket"))

        # Convert memory to megabytes.
        m = str(kwargs.pop("mem_per_node"))
        self.mem_per_node = int(Memory.from_string(m).to("Mb"))

        if self.mem_per_node <= 0 or self.sockets_per_node <= 0 or self.cores_per_socket <= 0:
            raise ValueError("invalid parameters: %s" % kwargs)

        if kwargs:
            raise ValueError("Found invalid keywords in the partition section:\n %s" % list(kwargs.keys()))

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        app("   num_nodes: %d, sockets_per_node: %d, cores_per_socket: %d, mem_per_node %s," %
            (self.num_nodes, self.sockets_per_node, self.cores_per_socket, self.mem_per_node))
        return "\n".join(lines)

    @property
    def num_cores(self):
        """Total number of cores available"""
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
        """Use divmod to compute (num_nodes, rest_cores)"""
        return divmod(mpi_procs * omp_threads, self.cores_per_node)

    def as_dict(self):
        return {'num_nodes': self.num_nodes,
                'sockets_per_node': self.sockets_per_node,
                'cores_per_socket': self.cores_per_socket,
                'mem_per_node': str(Memory(val=self.mem_per_node, unit='Mb'))}

    @classmethod
    def from_dict(cls, dd):
        return cls(num_nodes=dd['num_nodes'],
                   sockets_per_node=dd['sockets_per_node'],
                   cores_per_socket=dd['cores_per_socket'],
                   mem_per_node=dd['mem_per_node'])


class _ExcludeNodesFile(object):
    """
    This file contains the list of nodes to be excluded.
    Nodes are indexed by queue name.
    """
    DIRPATH = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")
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
            with AtomicFile(self.FILEPATH, mode="w+") as fh:
                d = json.load(fh)
                if qname in d:
                    d["qname"].extend(nodes)
                    d["qname"] = list(set(d["qname"]))
                else:
                    d["qname"] = nodes
                json.dump(d, fh)

_EXCL_NODES_FILE = _ExcludeNodesFile()


def show_qparams(qtype, stream=sys.stdout):
    """Print to the given stream the template of the :class:`QueueAdapter` of type `qtype`."""
    for cls in all_subclasses(QueueAdapter):
        if cls.QTYPE == qtype: return stream.write(cls.QTEMPLATE)

    raise ValueError("Cannot find class associated to qtype %s" % qtype)


def all_qtypes():
    """Return sorted list with all qtypes supported."""
    return sorted([cls.QTYPE for cls in all_subclasses(QueueAdapter)])


def make_qadapter(**kwargs):
    """
    Return the concrete :class:`QueueAdapter` class from a string.
    Note that one can register a customized version with:

    .. example::

        from qadapters import SlurmAdapter

        class MyAdapter(SlurmAdapter):
            QTYPE = "myslurm"
            # Add your customized code here

        # Register your class.
        SlurmAdapter.register(MyAdapter)

        make_qadapter(qtype="myslurm", **kwargs)

    .. warning::

        MyAdapter should be pickleable, hence one should declare it
        at the module level so that pickle can import it at run-time.
    """
    # Get all known subclasses of QueueAdapter.
    d = {c.QTYPE: c for c in all_subclasses(QueueAdapter)}

    # Preventive copy before pop
    kwargs = copy.deepcopy(kwargs)
    qtype = kwargs["queue"].pop("qtype")

    return d[qtype](**kwargs)


class QScriptTemplate(string.Template):
    delimiter = '$$'


class QueueAdapterError(Exception):
    """Base Error class for exceptions raise by QueueAdapter."""


class MaxNumLaunchesError(QueueAdapterError):
    """Raised by `submit_to_queue` if we try to submit more than `max_num_launches` times."""


class QueueAdapter(six.with_metaclass(abc.ABCMeta, MSONable)):
    """
    The `QueueAdapter` is responsible for all interactions with a specific queue management system.
    This includes handling all details of queue script format as well as queue submission and management.

    This is the **abstract** base class defining the methods that must be implemented by the concrete classes.
    Concrete classes should extend this class with implementations that work on specific queue systems.

    .. note::

        A `QueueAdapter` has a handler (:class:`QueueJob`) defined in qjobs.py that allows one
        to contact the resource manager to get info about the status of the job.
        Each concrete implementation of `QueueAdapter` should have a corresponding `QueueJob`.
    """
    Error = QueueAdapterError

    MaxNumLaunchesError = MaxNumLaunchesError

    @classmethod
    def all_qtypes(cls):
        """Return sorted list with all qtypes supported."""
        return sorted([subcls.QTYPE for subcls in all_subclasses(cls)])

    @classmethod
    def autodoc(cls):
        return """
# Dictionary with info on the hardware available on this queue.
hardware:
    num_nodes:           # Number of nodes available on this queue (integer, MANDATORY).
    sockets_per_node:    # Number of sockets per node (integer, MANDATORY).
    cores_per_socket:    # Number of cores per socket (integer, MANDATORY).
                         # The total number of cores available on this queue is
                         # `num_nodes * sockets_per_node * cores_per_socket`.

# Dictionary with the options used to prepare the enviroment before submitting the job
job:
    setup:                # List of commands (strings) executed before running (DEFAULT: empty)
    omp_env:              # Dictionary with OpenMP environment variables (DEFAULT: empty i.e. no OpenMP)
    modules:              # List of modules to be imported before running the code (DEFAULT: empty).
                          # NB: Error messages produced by module load are redirected to mods.err
    shell_env:            # Dictionary with shell environment variables.
    mpi_runner:           # MPI runner. Possible values in ["mpirun", "mpiexec", "srun", None]
                          # DEFAULT: None i.e. no mpirunner is used.
    mpi_runner_options    # String with optional options passed to the `mpi_runner` e.g. "--bind-to None"
    shell_runner:         # Used for running small sequential jobs on the front-end. Set it to None
                          # if mpirun or mpiexec are not available on the fron-end. If not
                          # given, small sequential jobs are executed with `mpi_runner`.
    shell_runner_options  # Similar to mpi_runner_options but for the runner used on the front-end.
    pre_run:              # List of commands (strings) executed before the run (DEFAULT: empty)
    post_run:             # List of commands (strings) executed after the run (DEFAULT: empty)

# dictionary with the name of the queue and optional parameters
# used to build/customize the header of the submission script.
queue:
    qtype:                # String defining the qapapter type e.g. slurm, shell ...
    qname:                # Name of the submission queue (string, MANDATORY)
    qparams:              # Dictionary with values used to generate the header of the job script
                          # We use the *normalized* version of the options i.e dashes in the official name
                          # are replaced by underscores e.g. ``--mail-type`` becomes ``mail_type``
                          # See pymatgen.io.abinit.qadapters.py for the list of supported values.
                          # Use ``qverbatim`` to pass additional options that are not included in the template.

# dictionary with the constraints that must be fulfilled in order to run on this queue.
limits:
    min_cores:             # Minimum number of cores (integer, DEFAULT: 1)
    max_cores:             # Maximum number of cores (integer, MANDATORY). Hard limit to hint_cores:
                           # it's the limit beyond which the scheduler will not accept the job (MANDATORY).
    hint_cores:            # The limit used in the initial setup of jobs.
                           # Fix_Critical method may increase this number until max_cores is reached
    min_mem_per_proc:      # Minimum memory per MPI process in Mb, units can be specified e.g. 1.4 Gb
                           # (DEFAULT: hardware.mem_per_core)
    max_mem_per_proc:      # Maximum memory per MPI process in Mb, units can be specified e.g. `1.4Gb`
                           # (DEFAULT: hardware.mem_per_node)
    timelimit:             # Initial time-limit. Accepts time according to slurm-syntax i.e:
                           # "days-hours" or "days-hours:minutes" or "days-hours:minutes:seconds" or
                           # "minutes" or "minutes:seconds" or "hours:minutes:seconds",
    timelimit_hard:        # The hard time-limit for this queue. Same format as timelimit.
                           # Error handlers could try to submit jobs with increased timelimit
                           # up to timelimit_hard. If not specified, timelimit_hard == timelimit
    condition:             # MongoDB-like condition (DEFAULT: empty, i.e. not used)
    allocation:            # String defining the policy used to select the optimal number of CPUs.
                           # possible values are in ["nodes", "force_nodes", "shared"]
                           # "nodes" means that we should try to allocate entire nodes if possible.
                           # This is a soft limit, in the sense that the qadapter may use a configuration
                           # that does not fulfill this requirement. In case of failure, it will try to use the
                           # smallest number of nodes compatible with the optimal configuration.
                           # Use `force_nodes` to enfore entire nodes allocation.
                           # `shared` mode does not enforce any constraint (DEFAULT: shared).
    max_num_launches:      # Limit to the number of times a specific task can be restarted (integer, DEFAULT: 5)
"""

    def __init__(self, **kwargs):
        """
        Args:
            qname: Name of the queue.
            qparams: Dictionary with the parameters used in the template.
            setup: String or list of commands to execute during the initial setup.
            modules: String or list of modules to load before running the application.
            shell_env: Dictionary with the environment variables to export before running the application.
            omp_env: Dictionary with the OpenMP variables.
            pre_run: String or list of commands to execute before launching the calculation.
            post_run: String or list of commands to execute once the calculation is completed.
            mpi_runner: Path to the MPI runner or :class:`MpiRunner` instance. None if not used
            mpi_runner_options: Optional string with options passed to the mpi_runner.
            max_num_launches: Maximum number of submissions that can be done for a specific task. Defaults to 5
            qverbatim:
            min_cores, max_cores, hint_cores: Minimum, maximum, and hint limits of number of cores that can be used
            min_mem_per_proc=Minimum memory per process in megabytes.
            max_mem_per_proc=Maximum memory per process in megabytes.
            timelimit: initial time limit in seconds
            timelimit_hard: hard limelimit for this queue
            priority: Priority level, integer number > 0
            condition: Condition object (dictionary)

        .. note::

            priority is a non-negative integer used to order the qadapters. The :class:`TaskManager` will
                try to run jobs on the qadapter with the highest priority if possible
        """
        # TODO
        #task_classes

        # Make defensive copies so that we can change the values at runtime.
        kwargs = copy.deepcopy(kwargs)
        self.priority = int(kwargs.pop("priority"))

        self.hw = Hardware(**kwargs.pop("hardware"))
        self._parse_queue(kwargs.pop("queue"))
        self._parse_limits(kwargs.pop("limits"))
        self._parse_job(kwargs.pop("job"))

        self.set_master_mem_overhead(kwargs.pop("master_mem_overhead", 0))

        # List of dictionaries with the parameters used to submit jobs
        # The launcher will use this information to increase the resources
        self.launches = []

        if kwargs:
            raise ValueError("Found unknown keywords:\n%s" % list(kwargs.keys()))

        self.validate_qparams()

        # Initialize some values from the info reported in the partition.
        self.set_mpi_procs(self.min_cores)
        self.set_mem_per_proc(self.min_mem_per_proc)

        # Final consistency check.
        self.validate_qparams()

    def as_dict(self):
        """
        Provides a simple though not complete dict serialization of the object (OMP missing, not all limits are
        kept in the dictionary, ... other things to be checked)

        Raise:
            `ValueError` if errors.
        """
        if self.has_omp:
            raise NotImplementedError('as_dict method of QueueAdapter not yet implemented when OpenMP is activated')
        return {'@module': self.__class__.__module__,
                '@class': self.__class__.__name__,
                'priority': self.priority,
                'hardware': self.hw.as_dict(),
                'queue': {'qtype': self.QTYPE,
                          'qname': self._qname,
                          'qnodes': self.qnodes,
                          'qparams': self._qparams},
                'limits': {'timelimit_hard': self._timelimit_hard,
                           'timelimit': self._timelimit,
                           'min_cores': self.min_cores,
                           'max_cores': self.max_cores,
                           'min_mem_per_proc': self.min_mem_per_proc,
                           'max_mem_per_proc': self.max_mem_per_proc,
                           'memory_policy': self.memory_policy
                           },
                'job': {},
                'mpi_procs': self._mpi_procs,
                'mem_per_proc': self._mem_per_proc,
                'master_mem_overhead': self._master_mem_overhead
                }

    @classmethod
    def from_dict(cls, dd):
        priority = dd.pop('priority')
        hardware = dd.pop('hardware')
        queue = dd.pop('queue')
        limits = dd.pop('limits')
        job = dd.pop('job')
        qa = make_qadapter(priority=priority, hardware=hardware, queue=queue, limits=limits, job=job)
        qa.set_mpi_procs(dd.pop('mpi_procs'))
        qa.set_mem_per_proc(dd.pop('mem_per_proc'))
        qa.set_master_mem_overhead(dd.pop('master_mem_overhead', 0))
        timelimit = dd.pop('timelimit', None)
        if timelimit is not None:
            qa.set_timelimit(timelimit=timelimit)
        dd.pop('@module', None)
        dd.pop('@class', None)
        if dd:
            raise ValueError("Found unknown keywords:\n%s" % list(dd.keys()))
        return qa

    def validate_qparams(self):
        """
        Check if the keys specified by the user in qparams are supported.

        Raise:
            `ValueError` if errors.
        """
        # No validation for ShellAdapter.
        if isinstance(self, ShellAdapter): return

        # Parse the template so that we know the list of supported options.
        err_msg = ""
        for param in self.qparams:
            if param not in self.supported_qparams:
                err_msg += "Unsupported QUEUE parameter name %s\n" % param
                err_msg += "Supported parameters:\n"
                for param_sup in self.supported_qparams:
                    err_msg += "    %s \n" % param_sup

        if err_msg:
            raise ValueError(err_msg)

    def _parse_limits(self, d):
        # Time limits.
        self.set_timelimit(qu.timelimit_parser(d.pop("timelimit")))
        tl_hard = d.pop("timelimit_hard",None)
        tl_hard = qu.timelimit_parser(tl_hard) if tl_hard is not None else self.timelimit
        self.set_timelimit_hard(tl_hard)

        # Cores
        self.min_cores = int(d.pop("min_cores", 1))
        self.max_cores = int(d.pop("max_cores"))
        self.hint_cores = int(d.pop("hint_cores", self.max_cores))
        self.memory_policy = d.pop("memory_policy", "mem")
        if self.min_cores > self.max_cores:
            raise ValueError("min_cores %s cannot be greater than max_cores %s" % (self.min_cores, self.max_cores))

        # Memory
        # FIXME: Neeed because autoparal 1 with paral_kgb 1 is not able to estimate memory
        self.min_mem_per_proc = qu.any2mb(d.pop("min_mem_per_proc", self.hw.mem_per_core))
        self.max_mem_per_proc = qu.any2mb(d.pop("max_mem_per_proc", self.hw.mem_per_node))

        # Misc
        self.max_num_launches = int(d.pop("max_num_launches", 5))
        self.condition = Condition(d.pop("condition", {}))
        self.allocation = d.pop("allocation", "shared")
        if self.allocation not in ("nodes", "force_nodes", "shared"):
            raise ValueError("Wrong value for `allocation` option")

        if d:
            raise ValueError("Found unknown keyword(s) in limits section:\n %s" % list(d.keys()))

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

        mpi_options = d.pop("mpi_runner_options", "")
        self.mpi_runner = d.pop("mpi_runner", None)
        if not isinstance(self.mpi_runner, MpiRunner):
            self.mpi_runner = MpiRunner(self.mpi_runner, options=mpi_options)

        self.shell_runner = d.pop("shell_runner", None)
        shell_runner_options = d.pop("shell_runner_options", "")
        if self.shell_runner is not None:
            self.shell_runner = MpiRunner(self.shell_runner, options=shell_runner_options)

        pre_run = d.pop("pre_run", None)
        if is_string(pre_run): pre_run = [pre_run]
        self.pre_run = pre_run[:] if pre_run is not None else []

        post_run = d.pop("post_run", None)
        if is_string(post_run): post_run = [post_run]
        self.post_run = post_run[:] if post_run is not None else []

        if d:
            raise ValueError("Found unknown keyword(s) in job section:\n %s" % list(d.keys()))

    def _parse_queue(self, d):
        # Init params
        qparams = d.pop("qparams", None)
        self._qparams = copy.deepcopy(qparams) if qparams is not None else {}

        self.set_qname(d.pop("qname", ""))
        self.qnodes = d.pop("qnodes", "standard")
        if self.qnodes not in ["standard", "shared", "exclusive"]:
            raise ValueError("Nodes must be either in standard, shared or exclusive mode "
                             "while qnodes parameter was {}".format(self.qnodes))
        if d:
            raise ValueError("Found unknown keyword(s) in queue section:\n %s" % list(d.keys()))

    def __str__(self):
        lines = ["%s:%s" % (self.__class__.__name__, self.qname)]
        app = lines.append
        app("Hardware:\n" + str(self.hw))
        #lines.extend(["qparams:\n", str(self.qparams)])
        if self.has_omp: app(str(self.omp_env))

        return "\n".join(lines)

    @property
    def qparams(self):
        """Dictionary with the parameters used to construct the header."""
        return self._qparams

    @lazy_property
    def supported_qparams(self):
        """
        Dictionary with the supported parameters that can be passed to the
        queue manager (obtained by parsing QTEMPLATE).
        """
        import re
        return re.findall(r"\$\$\{(\w+)\}", self.QTEMPLATE)

    @property
    def has_mpi(self):
        """True if we are using MPI"""
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
        """Deep copy of the object."""
        return copy.deepcopy(self)

    def record_launch(self, queue_id): # retcode):
        """Save submission"""
        self.launches.append(
            AttrDict(queue_id=queue_id, mpi_procs=self.mpi_procs, omp_threads=self.omp_threads,
                     mem_per_proc=self.mem_per_proc, timelimit=self.timelimit))
        return len(self.launches)

    def remove_launch(self, index):
        """Remove launch with the given index."""
        self.launches.pop(index)

    @property
    def num_launches(self):
        """Number of submission tried with this adapter so far."""
        return len(self.launches)

    @property
    def last_launch(self):
        """Return the last launch."""
        if len(self.launches) > 0:
            return self.launches[-1]
        else:
            return None

    def validate(self):
        """Validate the parameters of the run. Raises self.Error if invalid parameters."""
        errors = []
        app = errors.append

        if not self.hint_cores >= self.mpi_procs * self.omp_threads >= self.min_cores:
            app("self.hint_cores >= mpi_procs * omp_threads >= self.min_cores not satisfied")

        if self.omp_threads > self.hw.cores_per_node:
            app("omp_threads > hw.cores_per_node")

        if self.mem_per_proc > self.hw.mem_per_node:
            app("mem_mb >= self.hw.mem_per_node")

        if not self.max_mem_per_proc >= self.mem_per_proc >= self.min_mem_per_proc:
            app("self.max_mem_per_proc >= mem_mb >= self.min_mem_per_proc not satisfied")

        if self.priority <= 0:
            app("priority must be > 0")

        if not (1 <= self.min_cores <= self.hw.num_cores >= self.hint_cores):
            app("1 <= min_cores <= hardware num_cores >= hint_cores not satisfied")

        if errors:
            raise self.Error(str(self) + "\n".join(errors))

    def set_omp_threads(self, omp_threads):
        """Set the number of OpenMP threads."""
        self.omp_env["OMP_NUM_THREADS"] = omp_threads

    @property
    def mpi_procs(self):
        """Number of CPUs used for MPI."""
        return self._mpi_procs

    def set_mpi_procs(self, mpi_procs):
        """Set the number of MPI processes to mpi_procs"""
        self._mpi_procs = mpi_procs

    @property
    def qname(self):
        """The name of the queue."""
        return self._qname

    def set_qname(self, qname):
        """Set the name of the queue."""
        self._qname = qname

    # todo this assumes only one wall time. i.e. the one in the mananager file is the one always used.
    # we should use the standard walltime to start with but also allow to increase the walltime

    @property
    def timelimit(self):
        """Returns the walltime in seconds."""
        return self._timelimit

    @property
    def timelimit_hard(self):
        """Returns the walltime in seconds."""
        return self._timelimit_hard

    def set_timelimit(self, timelimit):
        """Set the start walltime in seconds, fix method may increase this one until timelimit_hard is reached."""
        self._timelimit = timelimit

    def set_timelimit_hard(self, timelimit_hard):
        """Set the maximal possible walltime in seconds."""
        self._timelimit_hard = timelimit_hard

    @property
    def mem_per_proc(self):
        """The memory per process in megabytes."""
        return self._mem_per_proc

    @property
    def master_mem_overhead(self):
        """The memory overhead for the master process in megabytes."""
        return self._master_mem_overhead

    def set_mem_per_proc(self, mem_mb):
        """
        Set the memory per process in megabytes. If mem_mb <=0, min_mem_per_proc is used.
        """
        # Hack needed because abinit is still not able to estimate memory.
        # COMMENTED by David.
        # This is not needed anymore here because the "hack" is performed directly in select_qadapter/_use_qadpos_pconf
        # methods of TaskManager. Moreover, this hack should be performed somewhere else (this part should be
        # independent of abinit ... and if we want to have less memory than the average memory available per node, we
        # have to allow it!)
        #if mem_mb <= self.min_mem_per_proc: mem_mb = self.min_mem_per_proc
        self._mem_per_proc = int(mem_mb)

    def set_master_mem_overhead(self, mem_mb):
        """
        Set the memory overhead for the master process in megabytes.
        """
        if mem_mb < 0:
            raise ValueError("Memory overhead for the master process should be >= 0")
        self._master_mem_overhead = int(mem_mb)

    @property
    def total_mem(self):
        """Total memory required by the job in megabytes."""
        return Memory(self.mem_per_proc * self.mpi_procs + self.master_mem_overhead, "Mb")

    @abc.abstractmethod
    def cancel(self, job_id):
        """
        Cancel the job.

        Args:
            job_id: Job identifier.

        Returns:
            Exit status.
        """

    def can_run_pconf(self, pconf):
        """True if the qadapter in principle is able to run the :class:`ParalConf` pconf"""
        if not self.hint_cores >= pconf.num_cores >= self.min_cores: return False
        if not self.hw.can_use_omp_threads(self.omp_threads): return False
        if pconf.mem_per_proc > self.hw.mem_per_node: return False
        if self.allocation == "force_nodes" and pconf.num_cores % self.hw.cores_per_node != 0:
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
        class Distrib(namedtuple("Distrib", "num_nodes mpi_per_node exact")):
            pass
            #@property
            #def mem_per_node
            #    return self.mpi_per_node * mem_per_proc
            #def set_nodes(self, nodes):

        hw = self.hw

        # TODO: Add check on user-memory
        if mem_per_proc <= 0:
            logger.warning("mem_per_proc <= 0")
            mem_per_proc = hw.mem_per_core

        if mem_per_proc > hw.mem_per_node:
            raise self.Error(
                "mem_per_proc > mem_per_node.\n Cannot distribute mpi_procs %d, omp_threads %d, mem_per_proc %s" %
                             (mpi_procs, omp_threads, mem_per_proc))

        # Try to use all then cores in the node.
        num_nodes, rest_cores = hw.divmod_node(mpi_procs, omp_threads)

        if num_nodes == 0 and mpi_procs * mem_per_proc <= hw.mem_per_node:
            # One node is enough
            return Distrib(num_nodes=1, mpi_per_node=mpi_procs, exact=True)

        if num_nodes == 0: num_nodes = 2
        mpi_per_node = mpi_procs // num_nodes
        if mpi_per_node * mem_per_proc <= hw.mem_per_node and rest_cores == 0:
            # Commensurate with nodes.
            return Distrib(num_nodes=num_nodes, mpi_per_node=mpi_per_node, exact=True)

        #if mode == "block", "cyclic"

        # Try first to pack MPI processors in a node as much as possible
        mpi_per_node = int(hw.mem_per_node / mem_per_proc)
        assert mpi_per_node != 0
        num_nodes = (mpi_procs * omp_threads) // mpi_per_node
        print("exact --> false", num_nodes, mpi_per_node)

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
            raise self.Error("Cannot distribute mpi_procs %d, omp_threads %d, mem_per_proc %s" %
                            (mpi_procs, omp_threads, mem_per_proc))

    def optimize_params(self, qnodes=None):
        """
        This method is called in get_subs_dict. Return a dict with parameters to be added to qparams
        Subclasses may provide a specialized version.
        """
        #logger.debug("optimize_params of baseclass --> no optimization available!!!")
        return {}

    def get_subs_dict(self, qnodes=None):
        """
        Return substitution dict for replacements into the template
        Subclasses may want to customize this method.
        """
        #d = self.qparams.copy()
        d = self.qparams
        d.update(self.optimize_params(qnodes=qnodes))
        # clean null values
        subs_dict = {k: v for k, v in d.items() if v is not None}
        #print("subs_dict:", subs_dict)
        return subs_dict

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
                       stdin=None, stdout=None, stderr=None, exec_args=None):
        """
        Returns a (multi-line) String representing the queue script, e.g. PBS script.
        Uses the template_file along with internal parameters to create the script.

        Args:
            job_name: Name of the job.
            launch_dir: (str) The directory the job will be launched in.
            executable: String with the name of the executable to be executed or list of commands
            qout_path Path of the Queue manager output file.
            qerr_path: Path of the Queue manager error file.
            exec_args: List of arguments passed to executable (used only if executable is a string, default: empty)
        """
        # PbsPro does not accept job_names longer than 15 chars.
        if len(job_name) > 14 and isinstance(self, PbsProAdapter):
            job_name = job_name[:14]

        # Construct the header for the Queue Manager.
        qheader = self._make_qheader(job_name, qout_path, qerr_path)

        # Add the bash section.
        se = ScriptEditor()

        # Cd to launch_dir immediately.
        se.add_line("cd " + os.path.abspath(launch_dir))

        if self.setup:
            se.add_comment("Setup section")
            se.add_lines(self.setup)
            se.add_emptyline()

        if self.modules:
            # stderr is redirected to mods.err file.
            # module load 2>> mods.err
            se.add_comment("Load Modules")
            se.add_line("module purge")
            se.load_modules(self.modules)
            se.add_emptyline()

        se.add_comment("OpenMp Environment")
        if self.has_omp:
            se.declare_vars(self.omp_env)
            se.add_emptyline()
        else:
            se.declare_vars({"OMP_NUM_THREADS": 1})

        if self.shell_env:
            se.add_comment("Shell Environment")
            se.declare_vars(self.shell_env)
            se.add_emptyline()

        if self.pre_run:
            se.add_comment("Commands before execution")
            se.add_lines(self.pre_run)
            se.add_emptyline()

        # Construct the string to run the executable with MPI and mpi_procs.
        if is_string(executable):
            line = self.mpi_runner.string_to_run(self, executable,
                                                 stdin=stdin, stdout=stdout, stderr=stderr, exec_args=exec_args)
            se.add_line(line)
        else:
            assert isinstance(executable, (list, tuple))
            se.add_lines(executable)

        if self.post_run:
            se.add_emptyline()
            se.add_comment("Commands after execution")
            se.add_lines(self.post_run)

        return qheader + se.get_script_str() + "\n"

    def submit_to_queue(self, script_file):
        """
        Public API: wraps the concrete implementation _submit_to_queue

        Raises:
            `self.MaxNumLaunchesError` if we have already tried to submit the job max_num_launches
            `self.Error` if generic error
        """
        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        if self.num_launches == self.max_num_launches:
            raise self.MaxNumLaunchesError("num_launches %s == max_num_launches %s" % (self.num_launches, self.max_num_launches))

        # Call the concrete implementation.
        s = self._submit_to_queue(script_file)
        self.record_launch(s.qid)

        if s.qid is None:
            raise self.Error("Error in job submission with %s. file %s \n" %
                            (self.__class__.__name__, script_file) +
                             "The error response reads:\n %s \n " % s.err +
                             "The out response reads:\n %s \n" % s.out)

        # Here we create a concrete instance of QueueJob
        return QueueJob.from_qtype_and_id(self.QTYPE, s.qid, self.qname), s.process

    @abc.abstractmethod
    def _submit_to_queue(self, script_file):
        """
        Submits the job to the queue, probably using subprocess or shutil
        This method must be provided by the concrete classes and will be called by submit_to_queue

        Args:
            script_file:  (str) name of the script file to use (String)

        Returns:
            queue_id, process
        """

    def get_njobs_in_queue(self, username=None):
        """
        returns the number of jobs in the queue, probably using subprocess or shutil to
        call a command like 'qstat'. returns None when the number of jobs cannot be determined.

        Args:
            username: (str) the username of the jobs to count (default is to autodetect)
        """
        if username is None: username = getpass.getuser()
        njobs, process = self._get_njobs_in_queue(username=username)

        if process is not None and process.returncode != 0:
            # there's a problem talking to squeue server?
            err_msg = ('Error trying to get the number of jobs in the queue' +
                       'The error response reads:\n {}'.format(process.stderr.read()))
            logger.critical(err_msg)

        if not isinstance(self, ShellAdapter):
            logger.info('The number of jobs currently in the queue is: {}'.format(njobs))

        return njobs

    @abc.abstractmethod
    def _get_njobs_in_queue(self, username):
        """
        Concrete Subclasses must implement this method. Return (njobs, process)
        """

    # Methods to fix problems
    def add_exclude_nodes(self, nodes):
        return _EXCL_NODES_FILE.add_nodes(self.qname, nodes)

    def get_exclude_nodes(self):
        return _EXCL_NODES_FILE.read_nodes(self.qname)

    @abc.abstractmethod
    def exclude_nodes(self, nodes):
        """Method to exclude nodes in the calculation. Return True if nodes have been excluded"""

    def more_mem_per_proc(self, factor=1):
        """
        Method to increase the amount of memory asked for, by factor.
        Return: new memory if success, 0 if memory cannot be increased.
        """
        base_increase = 2000
        old_mem = self.mem_per_proc
        new_mem = old_mem + factor*base_increase

        if new_mem < self.hw.mem_per_node:
            self.set_mem_per_proc(new_mem)
            return new_mem

        raise self.Error('could not increase mem_per_proc further')

    def more_master_mem_overhead(self, mem_increase_mb=1000):
        """
        Method to increase the amount of memory overheaded asked for the master node.
        Return: new master memory overhead if success, 0 if it cannot be increased.
        """
        old_master_mem_overhead = self.master_mem_overhead
        new_master_mem_overhead = old_master_mem_overhead + mem_increase_mb
        if new_master_mem_overhead + self.mem_per_proc < self.hw.mem_per_node:
            self.set_master_mem_overhead(new_master_mem_overhead)
            return new_master_mem_overhead

        raise self.Error('could not increase master_mem_overhead further')

    def more_cores(self, factor=1):
        """
        Method to increase the number of MPI procs.
        Return: new number of processors if success, 0 if processors cannot be increased.
        """
        # TODO : find a formula that works for all max_cores
        if self.max_cores > 40:
          base_increase = 4 * int(self.max_cores / 40)
        else:
          base_increase = 4

        new_cores = self.hint_cores + factor * base_increase

        if new_cores < self.max_cores:
            self.hint_cores = new_cores
            return new_cores

        raise self.Error('%s hint_cores reached limit on max_core %s' % (new_cores, self.max_cores))

    def more_time(self, factor=1):
        """
        Method to increase the wall time
        """
        base_increase = int(self.timelimit_hard / 10)

        new_time = self.timelimit + base_increase*factor
        print('qadapter: trying to increase time')
        if new_time < self.timelimit_hard:
            self.set_timelimit(new_time)
            print('new time set: ', new_time)
            return new_time

        self.priority = -1

        raise self.Error("increasing time is not possible, the hard limit has been reached")

####################
# Concrete classes #
####################


class ShellAdapter(QueueAdapter):
    """Simple Adapter used to submit runs through the shell."""
    QTYPE = "shell"

    QTEMPLATE = """\
#!/bin/bash
$${qverbatim}
"""

    def cancel(self, job_id):
        return os.system("kill -9 %d" % job_id)

    def _submit_to_queue(self, script_file):
        # submit the job, return process and pid.
        process = Popen(("/bin/bash", script_file), stderr=PIPE)
        return SubmitResults(qid=process.pid, out='no out in shell submission', err='no err in shell submission', process=process)

    def _get_njobs_in_queue(self, username):
        return None, None

    def exclude_nodes(self, nodes):
        return False


class SlurmAdapter(QueueAdapter):
    """Adapter for SLURM."""
    QTYPE = "slurm"

    QTEMPLATE = """\
#!/bin/bash

#SBATCH --partition=$${partition}
#SBATCH --job-name=$${job_name}
#SBATCH --nodes=$${nodes}
#SBATCH --total_tasks=$${total_tasks}
#SBATCH --ntasks=$${ntasks}
#SBATCH --ntasks-per-node=$${ntasks_per_node}
#SBATCH --cpus-per-task=$${cpus_per_task}
#####SBATCH --mem=$${mem}
#SBATCH --mem-per-cpu=$${mem_per_cpu}
#SBATCH --hint=$${hint}
#SBATCH --time=$${time}
#SBATCH	--exclude=$${exclude_nodes}
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
        super(SlurmAdapter, self).set_qname(qname)
        if qname:
            self.qparams["partition"] = qname

    def set_mpi_procs(self, mpi_procs):
        """Set the number of CPUs used for MPI."""
        super(SlurmAdapter, self).set_mpi_procs(mpi_procs)
        self.qparams["ntasks"] = mpi_procs

    def set_omp_threads(self, omp_threads):
        super(SlurmAdapter, self).set_omp_threads(omp_threads)
        self.qparams["cpus_per_task"] = omp_threads

    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in megabytes"""
        super(SlurmAdapter, self).set_mem_per_proc(mem_mb)
        self.qparams["mem_per_cpu"] = self.mem_per_proc
        # Remove mem if it's defined.
        #self.qparams.pop("mem", None)

    def set_timelimit(self, timelimit):
        super(SlurmAdapter, self).set_timelimit(timelimit)
        self.qparams["time"] = qu.time2slurm(timelimit)

    def cancel(self, job_id):
        return os.system("scancel %d" % job_id)

    def optimize_params(self, qnodes=None):
        params = {}
        if self.allocation == "nodes":
            # run on the smallest number of nodes compatible with the configuration
            params["nodes"] = max(int(math.ceil(self.mpi_procs / self.hw.cores_per_node)),
                                  int(math.ceil(self.total_mem / self.hw.mem_per_node)))
        return params

        #dist = self.distribute(self.mpi_procs, self.omp_threads, self.mem_per_proc)
        ##print(dist)

        #if False and dist.exact:
        #    # Can optimize parameters
        #    self.qparams["nodes"] = dist.num_nodes
        #    self.qparams.pop("ntasks", None)
        #    self.qparams["ntasks_per_node"] = dist.mpi_per_node
        #    self.qparams["cpus_per_task"] = self.omp_threads
        #    self.qparams["mem"] = dist.mpi_per_node * self.mem_per_proc
        #    self.qparams.pop("mem_per_cpu", None)
        #else:
        #    # Delegate to slurm.
        #    self.qparams["ntasks"] = self.mpi_procs
        #    self.qparams.pop("nodes", None)
        #    self.qparams.pop("ntasks_per_node", None)
        #    self.qparams["cpus_per_task"] = self.omp_threads
        #    self.qparams["mem_per_cpu"] = self.mem_per_proc
        #    self.qparams.pop("mem", None)
        #return {}

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        if sys.version_info[0] < 3:
            process = Popen(['sbatch', script_file], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['sbatch', script_file], stdout=PIPE, stderr=PIPE, universal_newlines=True)

        out, err = process.communicate()

        # grab the returncode. SLURM returns 0 if the job was successful
        queue_id = None
        if process.returncode == 0:
            try:
                # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
                queue_id = int(out.split()[3])
                logger.info('Job submission was successful and queue_id is {}'.format(queue_id))
            except:
                # probably error parsing job code
                logger.critical('Could not parse job id following slurm...')
        return SubmitResults(qid=queue_id, out=out, err=err, process=process)

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
            raise self.Error('qadapter failed to exclude nodes')

    def _get_njobs_in_queue(self, username):
        if sys.version_info[0] < 3:
            process = Popen(['squeue', '-o "%u"', '-u', username], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['squeue', '-o "%u"', '-u', username], stdout=PIPE, stderr=PIPE,
                            universal_newlines=True)

        out, err = process.communicate()
        njobs = None
        if process.returncode == 0:
            # parse the result. lines should have this form:
            # username
            # count lines that include the username in it
            outs = out.splitlines()
            njobs = len([line.split() for line in outs if username in line])

        return njobs, process


class PbsProAdapter(QueueAdapter):
    """Adapter for PbsPro"""
    QTYPE = "pbspro"

#PBS -l select=$${select}:ncpus=$${ncpus}:mem=$${mem}mb:mpiprocs=$${mpiprocs}:ompthreads=$${ompthreads}
#PBS -l select=$${select}:ncpus=1:mem=$${mem}mb:mpiprocs=1:ompthreads=$${ompthreads}

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
#PBS -o $${_qout_path}
#PBS -e $${_qerr_path}
$${qverbatim}
"""

    def set_qname(self, qname):
        super(PbsProAdapter, self).set_qname(qname)
        if qname:
            self.qparams["queue"] = qname

    def set_timelimit(self, timelimit):
        super(PbsProAdapter, self).set_timelimit(timelimit)
        self.qparams["walltime"] = qu.time2pbspro(timelimit)

    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in megabytes"""
        super(PbsProAdapter, self).set_mem_per_proc(mem_mb)
        #self.qparams["mem"] = self.mem_per_proc

    def cancel(self, job_id):
        return os.system("qdel %d" % job_id)

    def optimize_params(self, qnodes=None):
        return {"select": self.get_select(qnodes=qnodes)}

    def get_select(self, ret_dict=False, qnodes=None, memory_policy=None):
        """
        Select is not the most intuitive command. For more info see:

            * http://www.cardiff.ac.uk/arcca/services/equipment/User-Guide/pbs.html
            * https://portal.ivec.org/docs/Supercomputers/PBS_Pro
        """
        hw, mem_per_proc = self.hw, int(self.mem_per_proc)
        #dist = self.distribute(self.mpi_procs, self.omp_threads, mem_per_proc)
        """
        if self.pure_mpi:
            num_nodes, rest_cores = hw.divmod_node(self.mpi_procs, self.omp_threads)

            if num_nodes == 0:
                logger.info("IN_CORE PURE MPI: %s" % self.run_info)
                chunks = 1
                ncpus = rest_cores
                mpiprocs = rest_cores
                mem = mem_per_proc * ncpus
                ompthreads = 1

            elif rest_cores == 0:
                # Can allocate entire nodes because self.mpi_procs is divisible by cores_per_node.
                logger.info("PURE MPI run commensurate with cores_per_node %s" % self.run_info)
                chunks = num_nodes
                ncpus = hw.cores_per_node
                mpiprocs = hw.cores_per_node
                mem = ncpus * mem_per_proc
                ompthreads = 1

            else:
                logger.info("OUT-OF-CORE PURE MPI (not commensurate with cores_per_node): %s" % self.run_info)
                chunks = self.mpi_procs
                ncpus = 1
                mpiprocs = 1
                mem = mem_per_proc
                ompthreads = 1

        elif self.pure_omp:
            # Pure OMP run.
            logger.info("PURE OPENMP run: %s" % self.run_info)
            assert hw.can_use_omp_threads(self.omp_threads)
            chunks = 1
            ncpus = self.omp_threads
            mpiprocs = 1
            mem = mem_per_proc
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
                mem = mpiprocs * mem_per_proc
                ompthreads = self.omp_threads

            else:
                logger.info("HYBRID MPI-OPENMP, NOT commensurate with nodes: %s" % self.run_info)
                chunks=self.mpi_procs
                ncpus=self.omp_threads
                mpiprocs=1
                mem= mem_per_proc
                ompthreads=self.omp_threads

        else:
            raise RuntimeError("You should not be here")
        """
        if memory_policy is None:
            memory_policy = self.memory_policy
        if qnodes is None:
            qnodes = self.qnodes
        else:
            if qnodes not in ["standard", "shared", "exclusive"]:
                raise ValueError("Nodes must be either in standard, shared or exclusive mode "
                                 "while qnodes parameter was {}".format(self.qnodes))
        if qnodes == "standard":
            return self._get_select_standard(ret_dict=ret_dict, memory_policy=memory_policy)
        else:
            return self._get_select_with_master_mem_overhead(ret_dict=ret_dict, qnodes=qnodes,
                                                             memory_policy=memory_policy)

    def _get_select_with_master_mem_overhead(self, ret_dict=False, qnodes=None, memory_policy='mem'):
        if self.has_omp:
            raise NotImplementedError("select with master mem overhead not yet implemented with has_omp")
        if qnodes is None:
            qnodes = self.qnodes
        else:
            if qnodes not in ["standard", "shared", "exclusive"]:
                raise ValueError("Nodes must be either in standard, shared or exclusive mode "
                                 "while qnodes parameter was {}".format(self.qnodes))
        if qnodes == "exclusive":
            return self._get_select_with_master_mem_overhead_exclusive(ret_dict=ret_dict, memory_policy=memory_policy)
        elif qnodes == "shared":
            return self._get_select_with_master_mem_overhead_shared(ret_dict=ret_dict, memory_policy=memory_policy)
        else:
            raise ValueError("Wrong value of qnodes parameter : {}".format(self.qnodes))

    def _get_select_with_master_mem_overhead_shared(self, ret_dict=False, memory_policy='mem'):
        chunk_master, ncpus_master, vmem_master, mpiprocs_master = 1, 1, self.mem_per_proc+self.master_mem_overhead, 1
        if self.mpi_procs > 1:
            chunks_slaves, ncpus_slaves, vmem_slaves, mpiprocs_slaves = self.mpi_procs - 1, 1, self.mem_per_proc, 1
            select_params = AttrDict(chunk_master=chunk_master, ncpus_master=ncpus_master,
                                     mpiprocs_master=mpiprocs_master, vmem_master=int(vmem_master),
                                     chunks_slaves=chunks_slaves, ncpus_slaves=ncpus_slaves,
                                     mpiprocs_slaves=mpiprocs_slaves, vmem_slaves=int(vmem_slaves))
            if memory_policy == 'vmem':
                s = "{chunk_master}:ncpus={ncpus_master}:vmem={vmem_master}mb:mpiprocs={mpiprocs_master}+" \
                    "{chunks_slaves}:ncpus={ncpus_slaves}:vmem={vmem_slaves}mb:" \
                    "mpiprocs={mpiprocs_slaves}".format(**select_params)
            elif memory_policy == 'mem':
                s = "{chunk_master}:ncpus={ncpus_master}:mem={vmem_master}mb:mpiprocs={mpiprocs_master}+" \
                    "{chunks_slaves}:ncpus={ncpus_slaves}:mem={vmem_slaves}mb:" \
                    "mpiprocs={mpiprocs_slaves}".format(**select_params)
            tot_ncpus = chunk_master*ncpus_master + chunks_slaves*ncpus_slaves
            if tot_ncpus != self.mpi_procs:
                raise ValueError('Total number of cpus is different from mpi_procs ...')
        else:
            select_params = AttrDict(chunk_master=chunk_master, ncpus_master=ncpus_master,
                                     mpiprocs_master=mpiprocs_master, vmem_master=int(vmem_master))
            if memory_policy == 'vmem':
                s = "{chunk_master}:ncpus={ncpus_master}:vmem={vmem_master}mb:" \
                    "mpiprocs={mpiprocs_master}".format(**select_params)
            elif memory_policy == 'mem':
                s = "{chunk_master}:ncpus={ncpus_master}:mem={vmem_master}mb:" \
                    "mpiprocs={mpiprocs_master}".format(**select_params)
        if ret_dict:
            return s, select_params
        return s

    def _get_select_with_master_mem_overhead_exclusive(self, ret_dict=False, memory_policy='mem'):
        max_ncpus_master = min(self.hw.cores_per_node,
                               int((self.hw.mem_per_node-self.mem_per_proc-self.master_mem_overhead)
                                   / self.mem_per_proc) + 1)
        if max_ncpus_master >= self.mpi_procs:
            chunk, ncpus, mem, mpiprocs = 1, self.mpi_procs, self.hw.mem_per_node, self.mpi_procs
            if memory_policy == 'vmem':
                select_params = AttrDict(chunks=chunk, ncpus=ncpus, mpiprocs=mpiprocs, vmem=int(mem))
                s = "{chunks}:ncpus={ncpus}:vmem={vmem}mb:mpiprocs={mpiprocs}".format(**select_params)
            elif memory_policy == 'mem':
                select_params = AttrDict(chunks=chunk, ncpus=ncpus, mpiprocs=mpiprocs, mem=int(mem))
                s = "{chunks}:ncpus={ncpus}:mem={mem}mb:mpiprocs={mpiprocs}".format(**select_params)
            tot_ncpus = chunk*ncpus
        else:
            ncpus_left = self.mpi_procs-max_ncpus_master
            max_ncpus_per_slave_node = min(self.hw.cores_per_node, int(self.hw.mem_per_node/self.mem_per_proc))
            nslaves_float = float(ncpus_left)/float(max_ncpus_per_slave_node)
            ncpus_per_slave = max_ncpus_per_slave_node
            mpiprocs_slaves = max_ncpus_per_slave_node
            chunk_master = 1
            mem_slaves = self.hw.mem_per_node
            explicit_last_slave = False
            chunk_last_slave, ncpus_last_slave, mem_last_slave, mpiprocs_last_slave = None, None, None, None
            if nslaves_float > int(nslaves_float):
                chunks_slaves = int(nslaves_float) + 1
                pot_ncpus_all_slaves = chunks_slaves*ncpus_per_slave
                if pot_ncpus_all_slaves >= self.mpi_procs-1:
                    explicit_last_slave = True
                    chunks_slaves = chunks_slaves-1
                    chunk_last_slave = 1
                    ncpus_master = 1
                    ncpus_last_slave = self.mpi_procs - 1 - chunks_slaves*ncpus_per_slave
                    mem_last_slave = self.hw.mem_per_node
                    mpiprocs_last_slave = ncpus_last_slave
                else:
                    ncpus_master = self.mpi_procs-pot_ncpus_all_slaves
                if ncpus_master > max_ncpus_master:
                    raise ValueError('ncpus for the master node exceeds the maximum ncpus for the master ... this'
                                     'should not happen ...')
                if ncpus_master < 1:
                    raise ValueError('ncpus for the master node is 0 ... this should not happen ...')
            elif nslaves_float == int(nslaves_float):
                chunks_slaves = int(nslaves_float)
                ncpus_master = max_ncpus_master
            else:
                raise ValueError('nslaves_float < int(nslaves_float) ...')
            mem_master, mpiprocs_master = self.hw.mem_per_node, ncpus_master
            if explicit_last_slave:
                if memory_policy == 'vmem':
                    select_params = AttrDict(chunk_master=chunk_master, ncpus_master=ncpus_master,
                                             mpiprocs_master=mpiprocs_master, vmem_master=int(mem_master),
                                             chunks_slaves=chunks_slaves, ncpus_per_slave=ncpus_per_slave,
                                             mpiprocs_slaves=mpiprocs_slaves, vmem_slaves=int(mem_slaves),
                                             chunk_last_slave=chunk_last_slave, ncpus_last_slave=ncpus_last_slave,
                                             vmem_last_slave=int(mem_last_slave),
                                             mpiprocs_last_slave=mpiprocs_last_slave)
                    s = "{chunk_master}:ncpus={ncpus_master}:vmem={vmem_master}mb:mpiprocs={mpiprocs_master}+" \
                        "{chunks_slaves}:ncpus={ncpus_per_slave}:vmem={vmem_slaves}mb:mpiprocs={mpiprocs_slaves}+" \
                        "{chunk_last_slave}:ncpus={ncpus_last_slave}:vmem={vmem_last_slave}mb:" \
                        "mpiprocs={mpiprocs_last_slave}".format(**select_params)
                elif memory_policy == 'mem':
                    select_params = AttrDict(chunk_master=chunk_master, ncpus_master=ncpus_master,
                                             mpiprocs_master=mpiprocs_master, mem_master=int(mem_master),
                                             chunks_slaves=chunks_slaves, ncpus_per_slave=ncpus_per_slave,
                                             mpiprocs_slaves=mpiprocs_slaves, mem_slaves=int(mem_slaves),
                                             chunk_last_slave=chunk_last_slave, ncpus_last_slave=ncpus_last_slave,
                                             mem_last_slave=int(mem_last_slave),
                                             mpiprocs_last_slave=mpiprocs_last_slave)
                    s = "{chunk_master}:ncpus={ncpus_master}:mem={mem_master}mb:mpiprocs={mpiprocs_master}+" \
                        "{chunks_slaves}:ncpus={ncpus_per_slave}:mem={mem_slaves}mb:mpiprocs={mpiprocs_slaves}+" \
                        "{chunk_last_slave}:ncpus={ncpus_last_slave}:mem={mem_last_slave}mb:" \
                        "mpiprocs={mpiprocs_last_slave}".format(**select_params)
                tot_ncpus = chunk_master*ncpus_master+chunks_slaves*ncpus_per_slave+chunk_last_slave*ncpus_last_slave
            else:
                if memory_policy == 'vmem':
                    select_params = AttrDict(chunk_master=chunk_master, ncpus_master=ncpus_master,
                                             mpiprocs_master=mpiprocs_master, vmem_master=int(mem_master),
                                             chunks_slaves=chunks_slaves, ncpus_per_slave=ncpus_per_slave,
                                             mpiprocs_slaves=mpiprocs_slaves, vmem_slaves=int(mem_slaves))
                    s = "{chunk_master}:ncpus={ncpus_master}:vmem={vmem_master}mb:mpiprocs={mpiprocs_master}+" \
                        "{chunks_slaves}:ncpus={ncpus_per_slave}:vmem={vmem_slaves}mb:" \
                        "mpiprocs={mpiprocs_slaves}".format(**select_params)
                elif memory_policy == 'mem':
                    select_params = AttrDict(chunk_master=chunk_master, ncpus_master=ncpus_master,
                                             mpiprocs_master=mpiprocs_master, mem_master=int(mem_master),
                                             chunks_slaves=chunks_slaves, ncpus_per_slave=ncpus_per_slave,
                                             mpiprocs_slaves=mpiprocs_slaves, mem_slaves=int(mem_slaves))
                    s = "{chunk_master}:ncpus={ncpus_master}:mem={mem_master}mb:mpiprocs={mpiprocs_master}+" \
                        "{chunks_slaves}:ncpus={ncpus_per_slave}:mem={mem_slaves}mb:" \
                        "mpiprocs={mpiprocs_slaves}".format(**select_params)
                tot_ncpus = chunk_master*ncpus_master + chunks_slaves*ncpus_per_slave

        if tot_ncpus != self.mpi_procs:
            raise ValueError('Total number of cpus is different from mpi_procs ...')
        if ret_dict:
            return s, select_params
        return s

    def _get_select_standard(self, ret_dict=False, memory_policy='mem'):
        if not self.has_omp:
            chunks, ncpus, mem, mpiprocs = self.mpi_procs, 1, self.mem_per_proc, 1
            if memory_policy == 'vmem':
                select_params = AttrDict(chunks=chunks, ncpus=ncpus, mpiprocs=mpiprocs, vmem=int(mem))
                s = "{chunks}:ncpus={ncpus}:vmem={vmem}mb:mpiprocs={mpiprocs}".format(**select_params)
            elif memory_policy == 'mem':
                select_params = AttrDict(chunks=chunks, ncpus=ncpus, mpiprocs=mpiprocs, mem=int(mem))
                s = "{chunks}:ncpus={ncpus}:mem={mem}mb:mpiprocs={mpiprocs}".format(**select_params)
        else:
            chunks, ncpus, mem, mpiprocs, ompthreads = self.mpi_procs, self.omp_threads, self.mem_per_proc, 1, self.omp_threads
            if memory_policy == 'vmem':
                select_params = AttrDict(chunks=chunks, ncpus=ncpus, mpiprocs=mpiprocs, vmem=int(mem),
                                         ompthreads=ompthreads)
                s = "{chunks}:ncpus={ncpus}:vmem={vmem}mb:mpiprocs={mpiprocs}:ompthreads={ompthreads}".format(**select_params)
            elif memory_policy == 'mem':
                select_params = AttrDict(chunks=chunks, ncpus=ncpus, mpiprocs=mpiprocs, mem=int(mem),
                                         ompthreads=ompthreads)
                s = "{chunks}:ncpus={ncpus}:mem={mem}mb:mpiprocs={mpiprocs}:ompthreads={ompthreads}".format(
                    **select_params)

        if ret_dict:
            return s, select_params
        return s

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        if sys.version_info[0] < 3:
            process = Popen(['qsub', script_file], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['qsub', script_file], stdout=PIPE, stderr=PIPE, universal_newlines=True)

        out, err = process.communicate()
        # grab the return code. PBS returns 0 if the job was successful
        queue_id = None
        if process.returncode == 0:
            try:
                # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
                queue_id = int(out.split('.')[0])
            except:
                # probably error parsing job code
                logger.critical("Could not parse job id following qsub...")
        return SubmitResults(qid=queue_id, out=out, err=err, process=process)

    def _get_njobs_in_queue(self, username):
        if sys.version_info[0] < 3:
            process = Popen(['qstat', '-a', '-u', username], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['qstat', '-a', '-u', username], stdout=PIPE, stderr=PIPE, universal_newlines=True)

        out, err = process.communicate()
        njobs = None
        if process.returncode == 0:
            # parse the result
            # lines should have this form
            # '1339044.sdb          username  queuename    2012-02-29-16-43  20460   --   --    --  00:20 C 00:09'
            # count lines that include the username in it

            # TODO: only count running or queued jobs. or rather, *don't* count jobs that are 'C'.
            outs = out.split('\n')
            njobs = len([line.split() for line in outs if username in line])

        return njobs, process

    def exclude_nodes(self, nodes):
        return False


class TorqueAdapter(PbsProAdapter):
    """Adapter for Torque."""
    QTYPE = "torque"

    QTEMPLATE = """\
#!/bin/bash

#PBS -q $${queue}
#PBS -N $${job_name}
#PBS -A $${account}
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
        """Set the memory per process in megabytes"""
        QueueAdapter.set_mem_per_proc(self, mem_mb)
        #self.qparams["mem"] = self.mem_per_proc

    #@property
    #def mpi_procs(self):
    #    """Number of MPI processes."""
    #    return self.qparams.get("nodes", 1) * self.qparams.get("ppn", 1)

    def set_mpi_procs(self, mpi_procs):
        """Set the number of CPUs used for MPI."""
        QueueAdapter.set_mpi_procs(self, mpi_procs)
        self.qparams["nodes"] = 1
        self.qparams["ppn"] = mpi_procs

    def exclude_nodes(self, nodes):
        raise self.Error('qadapter failed to exclude nodes, not implemented yet in torque')


class SGEAdapter(QueueAdapter):
    """
    Adapter for Sun Grid Engine (SGE) task submission software.

    See also:

        * https://www.wiki.ed.ac.uk/display/EaStCHEMresearchwiki/How+to+write+a+SGE+job+submission+script
        * http://www.uibk.ac.at/zid/systeme/hpc-systeme/common/tutorials/sge-howto.html
    """
    QTYPE = "sge"

    QTEMPLATE = """\
#!/bin/bash

#$ -account_name $${account_name}
#$ -N $${job_name}
#$ -q $${queue_name}
#$ -pe $${parallel_environment} $${ncpus}
#$ -l h_rt=$${walltime}
# request a per slot memory limit of size bytes.
##$ -l h_vmem=$${mem_per_slot}
##$ -l mf=$${mem_per_slot}
###$ -j no
#$ -M $${mail_user}
#$ -m $${mail_type}
# Submission environment
##$ -S /bin/bash
###$ -cwd                       # Change to current working directory
###$ -V                         # Export environment variables into script
#$ -e $${_qerr_path}
#$ -o $${_qout_path}
$${qverbatim}
"""
    def set_qname(self, qname):
        super(SGEAdapter, self).set_qname(qname)
        if qname:
            self.qparams["queue_name"] = qname

    def set_mpi_procs(self, mpi_procs):
        """Set the number of CPUs used for MPI."""
        super(SGEAdapter, self).set_mpi_procs(mpi_procs)
        self.qparams["ncpus"] = mpi_procs

    def set_omp_threads(self, omp_threads):
        super(SGEAdapter, self).set_omp_threads(omp_threads)
        logger.warning("Cannot use omp_threads with SGE")

    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in megabytes"""
        super(SGEAdapter, self).set_mem_per_proc(mem_mb)
        self.qparams["mem_per_slot"] = str(int(self.mem_per_proc)) + "M"

    def set_timelimit(self, timelimit):
        super(SGEAdapter, self).set_timelimit(timelimit)
        # Same convention as pbspro e.g. [hours:minutes:]seconds
        self.qparams["walltime"] = qu.time2pbspro(timelimit)

    def cancel(self, job_id):
        return os.system("qdel %d" % job_id)

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        if sys.version_info[0] < 3:
            process = Popen(['qsub', script_file], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['qsub', script_file], stdout=PIPE, stderr=PIPE, universal_newlines=True)

        out, err = process.communicate()
        # grab the returncode. SGE returns 0 if the job was successful
        queue_id = None
        if process.returncode == 0:
            try:
                # output should of the form
                # Your job 1659048 ("NAME_OF_JOB") has been submitted
                queue_id = int(out.split(' ')[2])
            except:
                # probably error parsing job code
                logger.critical("Could not parse job id following qsub...")
        return SubmitResults(qid=queue_id, out=out, err=err, process=process)

    def exclude_nodes(self, nodes):
        """Method to exclude nodes in the calculation"""
        raise self.Error('qadapter failed to exclude nodes, not implemented yet in sge')

    def _get_njobs_in_queue(self, username):
        if sys.version_info[0] < 3:
            process = Popen(['qstat', '-u', username], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['qstat', '-u', username], stdout=PIPE, stderr=PIPE, universal_newlines=True)

        out, err = process.communicate()
        njobs = None
        if process.returncode == 0:
            # parse the result
            # lines should contain username
            # count lines that include the username in it

            # TODO: only count running or queued jobs. or rather, *don't* count jobs that are 'C'.
            outs = out.splitlines()
            njobs = len([line.split() for line in outs if username in line])

        return njobs, process


class MOABAdapter(QueueAdapter):
    """Adapter for MOAB. See https://computing.llnl.gov/tutorials/moab/"""
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
        super(MOABAdapter, self).set_mpi_procs(mpi_procs)
        self.qparams["procs"] = mpi_procs

    def set_timelimit(self, timelimit):
        super(MOABAdapter, self).set_timelimit(timelimit)
        self.qparams["walltime"] = qu.time2slurm(timelimit)

    def set_mem_per_proc(self, mem_mb):
        super(MOABAdapter, self).set_mem_per_proc(mem_mb)
        #TODO
        #raise NotImplementedError("set_mem_per_cpu")

    def exclude_nodes(self, nodes):
        raise self.Error('qadapter failed to exclude nodes, not implemented yet in moad')

    def cancel(self, job_id):
        return os.system("canceljob %d" % job_id)

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        if sys.version_info[0] < 3:
            process = Popen(['msub', script_file], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['msub', script_file], stdout=PIPE, stderr=PIPE, universal_newlines=True)

        out, err = process.communicate()
        queue_id = None
        if process.returncode == 0:
            # grab the returncode. MOAB returns 0 if the job was successful
            try:
                # output should be the queue_id
                queue_id = int(out.split()[0])
            except:
                # probably error parsing job code
                logger.critical('Could not parse job id following msub...')

        return SubmitResults(qid=queue_id, out=out, err=err, process=process)

    def _get_njobs_in_queue(self, username):
        if sys.version_info[0] < 3:
            process = Popen(['showq', '-s -u', username], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['showq', '-s -u', username], stdout=PIPE, stderr=PIPE, universal_newlines=True)

        out, err = process.communicate()
        njobs = None
        if process.returncode == 0:
            # parse the result
            # lines should have this form:
            ##
            ## active jobs: N  eligible jobs: M  blocked jobs: P
            ##
            ## Total job:  1
            ##
            # Split the output string and return the last element.
            out = out.splitlines()[-1]
            njobs = int(out.split()[-1])

        return njobs, process


class BlueGeneAdapter(QueueAdapter):
    """
    Adapter for LoadLever on BlueGene architectures.

    See:
        http://www.prace-ri.eu/best-practice-guide-blue-gene-q-html/#id-1.5.4.8
        https://www.lrz.de/services/compute/supermuc/loadleveler/
    """
    QTYPE = "bluegene"

    QTEMPLATE = """\
#!/bin/bash
# @ job_name = $${job_name}
# @ class = $${class}
# @ error = $${_qout_path}
# @ output = $${_qerr_path}
# @ wall_clock_limit = $${wall_clock_limit}
# @ notification = $${notification}
# @ notify_user = $${mail_user}
# @ environment = $${environment}
# @ account_no = $${account_no}
# @ job_type = bluegene
# @ bg_connectivity = $${bg_connectivity}
# @ bg_size = $${bg_size}
$${qverbatim}
# @ queue
"""

    def set_qname(self, qname):
        super(BlueGeneAdapter, self).set_qname(qname)
        if qname:
            self.qparams["class"] = qname

    #def set_mpi_procs(self, mpi_procs):
    #    """Set the number of CPUs used for MPI."""
    #    super(BlueGeneAdapter, self).set_mpi_procs(mpi_procs)
    #    #self.qparams["ntasks"] = mpi_procs

    #def set_omp_threads(self, omp_threads):
    #    super(BlueGeneAdapter, self).set_omp_threads(omp_threads)
    #    #self.qparams["cpus_per_task"] = omp_threads

    #def set_mem_per_proc(self, mem_mb):
    #    """Set the memory per process in megabytes"""
    #    super(BlueGeneAdapter, self).set_mem_per_proc(mem_mb)
    #    #self.qparams["mem_per_cpu"] = self.mem_per_proc

    def set_timelimit(self, timelimit):
        """Limits are specified with the format hh:mm:ss (hours:minutes:seconds)"""
        super(BlueGeneAdapter, self).set_timelimit(timelimit)
        self.qparams["wall_clock_limit"] = qu.time2loadlever(timelimit)

    def cancel(self, job_id):
        return os.system("llcancel %d" % job_id)

    def bgsize_rankspernode(self):
	    """Return (bg_size, ranks_per_node) from mpi_procs and omp_threads."""
	    bg_size = int(math.ceil((self.mpi_procs * self.omp_threads)/ self.hw.cores_per_node))
	    bg_size = max(bg_size, 32) # TODO hardcoded
	    ranks_per_node = int(math.ceil(self.mpi_procs / bg_size))

	    return bg_size, ranks_per_node

    def optimize_params(self, qnodes=None):
        params = {}
        bg_size, rpn = self.bgsize_rankspernode()
        print("in optimize params")
        print("mpi_procs:", self.mpi_procs, "omp_threads:",self.omp_threads)
        print("bg_size:",bg_size,"ranks_per_node",rpn)

        return {"bg_size": bg_size}

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        if sys.version_info[0] < 3:
            process = Popen(['llsubmit', script_file], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['llsubmit', script_file], stdout=PIPE, stderr=PIPE, universal_newlines=True)

        out, err = process.communicate()
        # grab the return code. llsubmit returns 0 if the job was successful
        queue_id = None
        if process.returncode == 0:
            try:
                # on JUQUEEN, output should of the form
                #llsubmit: The job "juqueen1c1.zam.kfa-juelich.de.281506" has been submitted.
                token = out.split()[3]
                s = token.split(".")[-1].replace('"', "")
                queue_id = int(s)
            except:
                # probably error parsing job code
                logger.critical("Could not parse job id following llsubmit...")
                raise

        return SubmitResults(qid=queue_id, out=out, err=err, process=process)

    def _get_njobs_in_queue(self, username):
        if sys.version_info[0] < 3:
            process = Popen(['llq', '-u', username], stdout=PIPE, stderr=PIPE)
        else:
            # need string not bytes so must use universal_newlines
            process = Popen(['llq', '-u', username], stdout=PIPE, stderr=PIPE, universal_newlines=True)

        out, err = process.communicate()
        njobs = None
        if process.returncode == 0:
            # parse the result. lines should have this form:
            #
            # Id                       Owner      Submitted   ST PRI Class        Running On
            # ------------------------ ---------- ----------- -- --- ------------ -----------
            # juqueen1c1.281508.0      paj15530    1/23 13:20 I  50  n001
            # 1 job step(s) in query, 1 waiting, 0 pending, 0 running, 0 held, 0 preempted
            #
            # count lines that include the username in it
            outs = out.split('\n')
            njobs = len([line.split() for line in outs if username in line])

        return njobs, process

    def exclude_nodes(self, nodes):
        return False
