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


import sys
import os
import abc
import string
import shlex
import copy
import getpass
import six
import json

from collections import namedtuple, OrderedDict, defaultdict
from subprocess import Popen, PIPE
from atomicfile import AtomicFile
from monty.string import is_string
from monty.collections import AttrDict
from monty.functools import lazy_property
from monty.io import FileLock
from pymatgen.core.units import Time, Memory
from .utils import Condition
from .launcher import ScriptEditor

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "MpiRunner",
    "make_qadapter",
]


def slurm_parse_timestr(s):
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
        `ValueError` if string is not valid.
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
        return slurm_parse_timestr(s)


def any2mb(s):
    """Convert string or number to memory in megabytes."""
    if is_string(s):
        return int(Memory.from_string(s).to("Mb"))
    else:
        return int(s)


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

    Basic definition::

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
            raise ValueError("Found invalid keywords in the partition section:\n %s" % kwargs.keys())

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
            with AtomicFile(self.FILEPATH, mode="w+") as fh:
                d = json.load(fh)
                if qname in d:
                    d["qname"].extend(nodes)
                    d["qname"] = list(set(d["qname"]))
                else:
                    d["qname"] = nodes
                json.dump(d, fh)

_EXCL_NODES_FILE = _ExcludeNodesFile()


class JobStatus(int):
    """
    This object is an integer representing the status of a :class:`QueueJob`.

    Slurm API, see `man squeue`.

    JOB STATE CODES
       Jobs typically pass through several states in the course of their execution.  The typical  states  are
       PENDING, RUNNING, SUSPENDED, COMPLETING, and COMPLETED.  An explanation of each state follows.

    BF  BOOT_FAIL       Job  terminated  due  to launch failure, typically due to a hardware failure (e.g.
                        unable to boot the node or block and the job can not be requeued).
    CA  CANCELLED       Job was explicitly cancelled by the user or system administrator.
                        The job may or may not have been initiated.
    CD  COMPLETED       Job has terminated all processes on all nodes.
    CF  CONFIGURING     Job has been allocated resources, but are waiting for them to become ready for use (e.g. booting).
    CG  COMPLETING      Job is in the process of completing. Some processes on some  nodes may still be active.
    F   FAILED          Job terminated with non-zero exit code or other failure condition.
    NF  NODE_FAIL       Job terminated due to failure of one or more allocated nodes.
    PD  PENDING         Job is awaiting resource allocation.
    PR  PREEMPTED       Job terminated due to preemption.
    R   RUNNING         Job currently has an allocation.
    S   SUSPENDED       Job has an allocation, but execution has been suspended.
    TO  TIMEOUT         Job terminated upon reaching its time limit.
    SE SPECIAL_EXIT     The job was requeued in a special state. This state can be set by users, typically
                        in EpilogSlurmctld, if the job has terminated with a particular exit value.
    """

    _STATUS_TABLE = OrderedDict([
        (-1, "UNKNOWN"),
        (0, "PENDING"),
        (1, "RUNNING"),
        (2, "RESIZING"),
        (3, "SUSPENDED"),
        (4, "COMPLETED"),
        (5, "CANCELLED"),
        (6, "FAILED"),
        (7, "TIMEOUT"),
        (8, "PREEMPTED"),
        (9, "NODEFAIL"),
    ])

    def __repr__(self):
        return "<%s: %s, at %s>" % (self.__class__.__name__, str(self), id(self))

    def __str__(self):
        """String representation."""
        return self._STATUS_TABLE[self]

    @classmethod
    def from_string(cls, s):
        """Return a :class`JobStatus` instance from its string representation."""
        for num, text in cls._STATUS_TABLE.items():
            if text == s: return cls(num)
        else:
            #raise ValueError("Wrong string %s" % s)
            logger.warning("Got unknown status: %s" % s)
            return cls.from_string("UNKNOWN")


class QueueJob(object):
    """
    This object provides methods to contact the resource manager to get info on the status
    of the job and useful statistics. This is an abstract class.
    """
    # Used to handle other resource managers.
    S_UNKNOWN   = JobStatus.from_string("UNKNOWN")
    # Slurm status
    S_PENDING   = JobStatus.from_string("PENDING")
    S_RUNNING   = JobStatus.from_string("RUNNING")
    S_RESIZING  = JobStatus.from_string("RESIZING")
    S_SUSPENDED = JobStatus.from_string("SUSPENDED")
    S_COMPLETED = JobStatus.from_string("COMPLETED")
    S_CANCELLED = JobStatus.from_string("CANCELLED")
    S_FAILED    = JobStatus.from_string("FAILED")
    S_TIMEOUT   = JobStatus.from_string("TIMEOUT")
    S_PREEMPTED = JobStatus.from_string("PREEMPTED")
    S_NODEFAIL  = JobStatus.from_string("NODEFAIL")

    def __init__(self, queue_id, qname, qout_path=None, qerr_path=None):
        self.qid, self.qname = queue_id, qname
        self.qout_path, self.qerr_path = qout_path, qerr_path

        # Initialize properties.
        self.status, self.exitcode, self.signal = None, None, None

    #def __str__(self):

    def __bool__(self):
        return self.qid is not None

    __nonzero__ = __bool__

    @property
    def is_completed(self):
        return self.status == self.S_COMPLETED

    @property
    def is_running(self):
        return self.status == self.S_RUNNING

    @property
    def is_failed(self):
        return self.status == self.S_FAILED

    @property
    def timeout(self):
        return self.status == self.S_TIMEOUT

    @property
    def has_node_failures(self):
        return self.status == self.S_NODEFAIL

    @property
    def unknown_status(self):
        return self.status == self.S_UNKNOWN

    def set_status_exitcode_signal(self, status, exitcode, signal):
        self.status, self.exitcode, self.signal = status, exitcode, signal

    def likely_code_error(self):
        """
        See http://man7.org/linux/man-pages/man7/signal.7.html

        SIGHUP        1       Term    Hangup detected on controlling terminal or death of controlling process
        SIGINT        2       Term    Interrupt from keyboard
        SIGQUIT       3       Core    Quit from keyboard
        SIGILL        4       Core    Illegal Instruction
        SIGABRT       6       Core    Abort signal from abort(3)
        SIGFPE        8       Core    Floating point exception
        SIGKILL       9       Term    Kill signal
        SIGSEGV      11       Core    Invalid memory reference
        SIGPIPE      13       Term    Broken pipe: write to pipe with no readers
        SIGALRM      14       Term    Timer signal from alarm(2)
        SIGTERM      15       Term    Termination signal
        SIGUSR1   30,10,16    Term    User-defined signal 1
        SIGUSR2   31,12,17    Term    User-defined signal 2
        SIGCHLD   20,17,18    Ign     Child stopped or terminated
        SIGCONT   19,18,25    Cont    Continue if stopped
        SIGSTOP   17,19,23    Stop    Stop process
        SIGTSTP   18,20,24    Stop    Stop typed at terminal
        SIGTTIN   21,21,26    Stop    Terminal input for background process
        SIGTTOU   22,22,27    Stop    Terminal output for background process

        The signals SIGKILL and SIGSTOP cannot be caught, blocked, or ignored.

        Next the signals not in the POSIX.1-1990 standard but described in
        SUSv2 and POSIX.1-2001.

        Signal       Value     Action   Comment
        ────────────────────────────────────────────────────────────────────
        SIGBUS      10,7,10     Core    Bus error (bad memory access)
        SIGPOLL                 Term    Pollable event (Sys V).
                                        Synonym for SIGIO
        SIGPROF     27,27,29    Term    Profiling timer expired
        SIGSYS      12,31,12    Core    Bad argument to routine (SVr4)
        SIGTRAP        5        Core    Trace/breakpoint trap
        SIGURG      16,23,21    Ign     Urgent condition on socket (4.2BSD)
        SIGVTALRM   26,26,28    Term    Virtual alarm clock (4.2BSD)
        SIGXCPU     24,24,30    Core    CPU time limit exceeded (4.2BSD)
        SIGXFSZ     25,25,31    Core    File size limit exceeded (4.2BSD)
        """
        for sig_name in ("SIGFPE",):
            if self.received_signal(sig_name): return sig_name
        return False

    def received_signal(self, sig_name):
        if self.signal is None: return False
        # Get the numeric value from signal and compare it with self.signal
        import signal
        try:
            return self.signal == getattr(signal, sig_name) 
        except AttributeError:
            # invalid sig_name or sig_name not available on this OS.
            return False

    def estimated_start_time(self):
        """Return date with estimated start time. None if it cannot be detected"""
        return None

    def get_info(self, **kwargs):
        return None

    def get_nodes(self, **kwargs):
        return None

    def get_stats(self, **kwargs):
        return None


class SlurmJob(QueueJob):
    """Handler for Slurm jobs."""

    def estimated_start_time(self):
        #squeue  --start -j  116791
        #  JOBID PARTITION     NAME     USER  ST           START_TIME  NODES NODELIST(REASON)
        # 116791      defq gs6q2wop cyildiri  PD  2014-11-04T09:27:15     16 (QOSResourceLimit)
        cmd = "squeue" "--start", "--job %d" % self.qid
        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        process.wait()

        if process.returncode != 0: return None
        lines = process.stdout.readlines()
        if len(lines) <= 2: return None

        from datetime import datetime
        for line in lines:
            tokens = line.split()
            if int(tokens[0]) == self.qid:
                date_string = tokens[5]
                if date_string == "N/A": return None
                return datetime.strptime(date_string, "%Y-%m-%dT%H:%M:%S")
        return None

    def get_info(self, **kwargs):
        # See https://computing.llnl.gov/linux/slurm/sacct.html
        #If SLURM job ids are reset, some job numbers will
        #probably appear more than once refering to different jobs.
        #Without this option only the most recent jobs will be displayed.

        #state Displays the job status, or state.
        #Output can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, 
        #PREEMPTED or NODE_FAIL. If more information is available on the job state than will fit 
        #into the current field width (for example, the uid that CANCELLED a job) the state will be followed by a "+". 

        #gmatteo@master2:~
        #sacct --job 112367 --format=jobid,exitcode,state --allocations --parsable2
        #JobID|ExitCode|State
        #112367|0:0|RUNNING
        #scontrol show job 800197 --oneliner

        # For more info
        #login1$ scontrol show job 1676354

        #cmd = "sacct --job %i --format=jobid,exitcode,state --allocations --parsable2" % self.qid
        cmd = "scontrol show job %i --oneliner" % self.qid
        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        process.wait()

        if process.returncode != 0:
            #print(process.stderr.readlines())
            return None
        line = process.stdout.read()
        #print("line", line)

        tokens = line.split()
        info = AttrDict()
        for line in tokens:
            #print(line)
            k, v = line.split("=")
            info[k] = v
            #print(info)

        qid = int(info.JobId)
        assert qid == self.qid
        exitcode = info.ExitCode
        status = info.JobState

        if ":" in exitcode:
            exitcode, signal = map(int, exitcode.split(":"))
        else:
            exitcode, signal = int(exitcode), None

        i = status.find("+")
        if i != -1: status = status[:i]

        self.set_status_exitcode_signal(JobStatus.from_string(status), exitcode, signal)
        return AttrDict(exitcode=exitcode, signal=signal, status=status)

    def get_stats(self, **kwargs):
        cmd = "sacct --long --job %s --parsable2" % self.qid
        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        process.wait()

        lines = process.stdout.readlines()
        keys = lines[0].strip().split("|")
        values = lines[1].strip().split("|")
        #print("lines0", lines[0])
        return dict(zip(keys, values))


class PbsProJob(QueueJob):
    """Handler for PbsPro Jobs"""
    # Mapping PrbPro --> Slurm. From `man qstat`
    #
    # S  The job’s state:
    #      B  Array job has at least one subjob running.
    #      E  Job is exiting after having run.
    #      F  Job is finished.
    #      H  Job is held.
    #      M  Job was moved to another server.
    #      Q  Job is queued.
    #      R  Job is running.
    #      S  Job is suspended.
    #      T  Job is being moved to new location.
    #      U  Cycle-harvesting job is suspended due to keyboard activity.
    #      W  Job is waiting for its submitter-assigned start time to be reached.
    #      X  Subjob has completed execution or has been deleted.

    PBSSTAT_TO_SLURM = defaultdict(lambda x: QueueJob.S_UNKNOWN, [
        ("E", QueueJob.S_FAILED),
        ("F", QueueJob.S_COMPLETED),
        ("Q", QueueJob.S_PENDING),
        ("R", QueueJob.S_RUNNING),
        ("S", QueueJob.S_SUSPENDED),
    ])

    def estimated_start_time(self):
        # qstat -T - Shows the estimated start time for all jobs in the queue. 
        #                                                                           Est
        #                                                            Req'd  Req'd   Start
        #Job ID          Username Queue    Jobname    SessID NDS TSK Memory Time  S Time
        #--------------- -------- -------- ---------- ------ --- --- ------ ----- - -----
        #5669001.frontal username large    gs.Pt         --   96  96    --  03:00 Q    --
        process = Popen(["qstat" "-T", str(self.qid)], stdout=PIPE, stderr=PIPE)
        process.wait()

        if process.returncode != 0: return None
        line = process.stdout.readlines()[-1]
        sdate = line.split()[-1]
        if sdate in ("--", "?"): return None
        # TODO One should convert to datetime
        return sdate

    def get_info(self, **kwargs):
        #$> qstat 5666289
        #frontal1:
        #                                                            Req'd  Req'd   Elap
        #Job ID          Username Queue    Jobname    SessID NDS TSK Memory Time  S Time
        #--------------- -------- -------- ---------- ------ --- --- ------ ----- - -----
        #5666289.frontal username main_ivy MorfeoTChk  57546   1   4    --  08:00 R 00:17
        cmd = "qstat %d" % self.qid
        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)

        if process.returncode != 0:
            #print(process.stderr.readlines())
            return None

        out = process.stdout.readlines()[-1]
        status = self.PBSSTAT_TO_SLURM[out.split()[9]]

        # Exit code and signal are not available.
        self.set_status_exitcode_signal(status, None, None)


def all_subclasses(cls):
    """
    Given a class `cls`, this recursive function returns a list with
    all subclasses, subclasses of subclasses, and so on.
    """
    subclasses = cls.__subclasses__() 
    return subclasses + [g for s in subclasses for g in all_subclasses(s)]


def show_qparams(qtype, stream=sys.stdout):
    """Print to the given stream the template of the :class:`QueueAdapter` of type `qtype`."""
    for cls in all_subclasses(QueueAdapter):
        if cls.QTYPE == qtype: return stream.write(cls.QTEMPLATE)

    raise ValueError("Cannot find class associated to qtype %s" % qtype)


def all_qtypes():
    """List of qtype supported."""
    return [cls.QTYPE for cls in all_subclasses(QueueAdapter)]


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
    qtype = kwargs["queue"].pop("qtype")

    return d[qtype](**kwargs)


class QueueAdapterError(Exception):
    """Error class for exceptions raise by QueueAdapter."""


class QueueAdapter(six.with_metaclass(abc.ABCMeta, object)):
    """
    The `QueueAdapter` is responsible for all interactions with a specific queue management system.
    This includes handling all details of queue script format as well as queue submission and management.

    This is the **abstract** base class defining the methods that must be implemented by the concrete classes.
    Concrete classes should extend this class with implementations that work on specific queue systems.
    """
    Error = QueueAdapterError

    Job = QueueJob

    @classmethod
    def autodoc(cls):
        return """
# dictionary with info on the hardware available on this particular queue.
hardware:  
    num_nodes:        # Number of nodes available on this queue. Mandatory
    sockets_per_node: # Self-explanatory. Mandatory.
    cores_per_socket: # Self-explanatory. Mandatory.

# dictionary with the options used to prepare the enviroment before submitting the job
job:
    setup:       # List of commands (str) executed before running (default empty)
    omp_env:     # Dictionary with OpenMP env variables (default empty i.e. no OpenMP)
    modules:     # List of modules to be imported (default empty)
    shell_env:   # Dictionary with shell env variables.
    mpi_runner:  # MPI runner i.e. mpirun, mpiexec, Default is None i.e. no mpirunner
    pre_run:     # List of commands executed before the run (default: empty)
    post_run:    # List of commands executed after the run (default: empty)

# dictionary with the name of the queue and optional parameters 
# used to build/customize the header of the submission script.
queue:
    qname:   # Name of the queue (mandatory)
    qparams: # Dictionary with values used to generate the header of the job script
             # See pymatgen.io.abinitio.qadapters.py for the list of supported values.

# dictionary with the constraints that must be fulfilled in order to run on this queue.
limits:
    min_cores:         # Minimum number of cores (default 1)
    max_cores:         # Maximum number of cores (mandatory)
    min_mem_per_proc:  # Minimum memory per MPI process in Mb, units can be specified e.g. 1.4 Gb
                       # (default hardware.mem_per_core)
    max_mem_per_proc:  # Maximum memory per MPI process in Mb, units can be specified e.g. `1.4Gb`
                       # (default hardware.mem_per_node)
    condition:         # MongoDB-like condition (default empty, i.e. not used)
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
            max_num_attempts: Default to 10
            qverbatim:
            min_cores, max_cores: Minimum and maximum number of cores that can be used
            min_mem_per_proc=Minimun memory per process in megabytes.
            max_mem_per_proc=Maximum memory per process in megabytes.
            timelimit: Time limit in seconds
            priority: Priority level, integer number > 0
            allocate_nodes: True if we must allocate entire nodes"
            condition: Condition object (dictionary)

        .. note::

            priority is a non-negative integer used to order the qadapters. The :class:`TaskManager` will
                try to run jobs on the qadapter with the highest priority if possible
        """
        # TODO
        #max_num_attempts:
        #task_classes

        # Make defensive copies so that we can change the values at runtime.
        kwargs = copy.deepcopy(kwargs)
        self.priority = int(kwargs.pop("priority"))

        self.hw = Hardware(**kwargs.pop("hardware"))
        self._parse_queue(kwargs.pop("queue"))
        self._parse_limits(kwargs.pop("limits"))
        self._parse_job(kwargs.pop("job"))

        if kwargs:
            raise ValueError("Found unknown keywords:\n%s" % kwargs.keys())

        self.validate_qparams()

        # List of dictionaries with the parameters used to submit jobs
        # The launcher will use this information to increase the resources
        self.attempts, self.max_num_attempts = [], kwargs.pop("max_num_attempts", 10)

        # Initialize some values from the info reported in the partition.
        self.set_mpi_procs(self.min_cores)
        self.set_mem_per_proc(self.min_mem_per_proc)

        # Final consistency check.
        self.validate_qparams()

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
                err_msg += "Supported are: \n"
                for param_sup in self.supported_qparams:
                    err_msg += "    %s \n" % param_sup

        if err_msg:
            raise ValueError(err_msg)

    def _parse_limits(self, d):
        self.set_timelimit(timelimit_parser(d.pop("timelimit")))
        self.min_cores = int(d.pop("min_cores", 1))
        self.max_cores = int(d.pop("max_cores"))
        # FIXME: Neeed because autoparal 1 with paral_kgb 1 is not able to estimate memory 
        self.min_mem_per_proc = any2mb(d.pop("min_mem_per_proc", self.hw.mem_per_core))
        self.max_mem_per_proc = any2mb(d.pop("max_mem_per_proc", self.hw.mem_per_node))
        #self.allocate_nodes = bool(d.pop("allocate_nodes", False))
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
        self._qparams = copy.deepcopy(qparams) if qparams is not None else {}

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
        """Dictionary with the parameters used to construct the header."""
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
        """True if we are using MPI"""
        return bool(self.mpi_runner)

    #@property
    #@deprecated(message="use has_mpi")
    #def has_mpirun(self):
    #    """True if we are using a mpirunner"""
    #    return bool(self.mpi_runner)

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

    def record_attempt(self, queue_id): # retcode):
        """Save submission"""
        self.attempts.append(
            AttrDict(queue_id=queue_id, mpi_procs=self.mpi_procs, omp_threads=self.omp_threads,
                     mem_per_proc=self.mem_per_proc, timelimit=self.timelimit))
        return len(self.attempts)

    def remove_attempt(self, index):
        """Remove attempt with the given index."""
        self.attempts.pop(index)

    @property
    def num_attempts(self):
        """Number of submission tried so far."""
        return len(self.attempts)

    @property
    def last_attempt(self):
        """Return the last attempt."""
        if len(self.attempts) > 0:
            return self.attempts[-1]
        else:
            return None

    def validate(self):
        """Validate the parameters of the run. Raises self.Error if invalid parameters."""
        errors = []
        app = errors.append

        if not self.max_cores >= self.mpi_procs * self.omp_threads >= self.min_cores:
            app("self.max_cores >= mpi_procs * omp_threads >= self.min_cores not satisfied")

        if self.omp_threads > self.hw.cores_per_node:
            app("omp_threads > hw.cores_per_node")

        if self.mem_per_proc > self.hw.mem_per_node:
            app("mem_mb >= self.hw.mem_per_node")

        if not self.max_mem_per_proc >= self.mem_per_proc >= self.min_mem_per_proc:
            app("self.max_mem_per_proc >= mem_mb >= self.min_mem_per_proc not satisfied")

        if self.priority <= 0: 
            app("priority must be > 0")

        if not (1 <= self.min_cores <= self.hw.num_cores >= self.max_cores):
            app("1 <= min_cores <= hardware num_cores >= max_cores not satisfied")

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

    @property
    def timelimit(self):
        """Returns the walltime in seconds."""
        return self._timelimit

    def set_timelimit(self, timelimit):
        """Set the walltime in seconds."""
        self._timelimit = timelimit

    @property
    def mem_per_proc(self):
        """The memory per process in megabytes."""
        return self._mem_per_proc
                                                
    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in megabytes"""
        # Hack needed because abinit is still not able to estimate memory.
        if mem_mb <= 0:
            mem_mb = self.min_mem_per_proc

        self._mem_per_proc = mem_mb

    @property
    def total_mem(self):
        """Total memory required by the job in megabytes."""
        return Memory(self.mem_per_proc * self.mpi_procs, "Mb")

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
        if not self.max_cores >= pconf.num_cores >= self.min_cores: return False
        if not self.hw.can_use_omp_threads(self.omp_threads): return False
        if pconf.mem_per_proc > self.hw.mem_per_node: return False
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

    def optimize_params(self):
        """
        This method is called in get_subs_dict. Return a dict with parameters to be added to qparams
        Subclasses may provide a specialized version.
        """
        logger.debug("optimize_params of baseclass --> no optimization available!!!")
        return {}

    def get_subs_dict(self):
        """
        Return substitution dict for replacements into the template
        Subclasses may want to customize this method.
        """ 
        #d = self.qparams.copy()
        d = self.qparams
        d.update(self.optimize_params())
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
                       stdin=None, stdout=None, stderr=None):
        """
        Returns a (multi-line) String representing the queue script, e.g. PBS script.
        Uses the template_file along with internal parameters to create the script.

        Args:
            job_name: Name of the job.
            launch_dir: (str) The directory the job will be launched in.
            executable: String with the name of the executable to be executed.
            qout_path Path of the Queue manager output file.
            qerr_path: Path of the Queue manager error file.
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
            `QueueAdapterError` if we have already tried to submit the job max_num_attempts
        """
        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        if self.num_attempts == self.max_num_attempts:
            raise self.Error("num_attempts %s == max_num_attempts %s" % (self.num_attempts, self.max_num_attempts))

        # Call the concrete implementation.
        queue_id, process = self._submit_to_queue(script_file)
        self.record_attempt(queue_id)

        if queue_id is None:
            submit_err_file = script_file + ".err"
            err = str(process.stderr.read())
            # Dump the error and raise
            with open(submit_err_file, mode='w') as f:
                f.write("sbatch submit process stderr:\n" + err)
                f.write("qparams:\n" + str(self.qparams))

            try:
                args = process.args
            except AttributeError:
                args = ["Unknown",]

            raise self.Error("Error in job submission with %s. file %s and args %s\n" % 
                             (self.__class__.__name__, script_file, args) + 
                             "The error response reads:\n %s" % err)

        # Here we create a concrete instance of QueueJob
        return self.Job(queue_id, self.qname), process

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
        """
        Method to exclude nodes in the calculation
        """

    def more_mem_per_proc(self, factor=1):
        """
        Method to increase the amount of memory asked for, by factor.
        Return: True if success.
        """
        base_increase = 2000
        old_mem = self.mem_per_proc
        new_mem = old_mem + factor*base_increase

        if new_mem < self.hw.mem_per_node:
            self.set_mem_per_proc(new_mem)
            return True

        logger.warning('could not increase mem_per_proc further')
        return False

    def more_mpi_procs(self, factor=1):
        """
        Method to increase the number of MPI procs.
        Return: True if success.
        """
        base_increase = 12
        new_cpus = self.mpi_procs + factor * base_increase

        if new_cpus * self.omp_threads < self.max_cores:
            self.set_mpi_procs(new_cpus)
            return True

        logger.warning('more_mpi_procs reached the limit')
        return False


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
        return process.pid, process

    def _get_njobs_in_queue(self, username):
        return None, None

    def exclude_nodes(self, nodes):
        return False


class SlurmAdapter(QueueAdapter):
    """Adapter for SLURM."""
    QTYPE = "slurm"

    Job = SlurmJob

    QTEMPLATE = """\
#!/bin/bash

#SBATCH --partition=$${partition}
#SBATCH --job-name=$${job_name}
#SBATCH --nodes=$${nodes}
#SBATCH --total_tasks=$${total_tasks}
#SBATCH --ntasks=$${ntasks}
#SBATCH --ntasks-per-node=$${ntasks_per_node}
#SBATCH --cpus-per-task=$${cpus_per_task}
#SBATCH --mem=$${mem}
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
        self.qparams["mem_per_cpu"] = int(mem_mb)
        # Remove mem if it's defined.
        self.qparams.pop("mem", None)

    def set_timelimit(self, timelimit):
        super(SlurmAdapter, self).set_timelimit(timelimit)
        self.qparams["time"] = time2slurm(timelimit)

    def cancel(self, job_id):
        return os.system("scancel %d" % job_id)

    def optimize_params(self):
        return {}
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
        process = Popen(['sbatch', script_file], stdout=PIPE, stderr=PIPE)
        process.wait()

        # grab the returncode. SLURM returns 0 if the job was successful
        queue_id = None
        if process.returncode == 0:
            try:
                # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
                queue_id = int(process.stdout.read().split()[3])
                logger.info('Job submission was successful and queue_id is {}'.format(queue_id))
            except:
                # probably error parsing job code
                logger.critical('Could not parse job id following slurm...')

        return queue_id, process

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

    def _get_njobs_in_queue(self, username):
        process = Popen(['squeue', '-o "%u"', '-u', username], stdout=PIPE, stderr=PIPE)
        process.wait()

        njobs = None
        if process.returncode == 0:
            # parse the result. lines should have this form:
            # username
            # count lines that include the username in it
            outs = process.stdout.readlines()
            njobs = len([line.split() for line in outs if username in line])

        return njobs, process


class PbsProAdapter(QueueAdapter):
    """Adapter for PbsPro"""
    QTYPE = "pbspro"

#PBS -l select=$${select}:ncpus=$${ncpus}:vmem=$${vmem}mb:mpiprocs=$${mpiprocs}:ompthreads=$${ompthreads}
#PBS -l select=$${select}:ncpus=1:vmem=$${vmem}mb:mpiprocs=1:ompthreads=$${ompthreads}
####PBS -l select=$${select}:ncpus=$${ncpus}:vmem=$${vmem}mb:mpiprocs=$${mpiprocs}:ompthreads=$${ompthreads}
####PBS -l pvmem=$${pvmem}mb

    Job = PbsProJob

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
        self.qparams["queue"] = qname

    def set_timelimit(self, timelimit):
        super(PbsProAdapter, self).set_timelimit(timelimit)
        self.qparams["walltime"] = time2pbspro(timelimit)

    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in megabytes"""
        super(PbsProAdapter, self).set_mem_per_proc(mem_mb)
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
        """

        if not self.has_omp:
            chunks, ncpus, vmem, mpiprocs = self.mpi_procs, 1, self.mem_per_proc, 1
            select_params = AttrDict(chunks=chunks, ncpus=ncpus, mpiprocs=mpiprocs, vmem=int(vmem))
            s = "{chunks}:ncpus={ncpus}:vmem={vmem}mb:mpiprocs={mpiprocs}".format(**select_params)
        else:
            chunks, ncpus, vmem, mpiprocs, ompthreads = self.mpi_procs, self.omp_threads, self.mem_per_proc, 1, self.omp_threads
            select_params = AttrDict(chunks=chunks, ncpus=ncpus, mpiprocs=mpiprocs, vmem=int(vmem), ompthreads=ompthreads)
            s = "{chunks}:ncpus={ncpus}:vmem={vmem}mb:mpiprocs={mpiprocs}:ompthreads={ompthreads}".format(**select_params)

        if ret_dict:
            return s, select_params
        return s

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        process = Popen(['qsub', script_file], stdout=PIPE, stderr=PIPE)
        process.wait()

        # grab the return code. PBS returns 0 if the job was successful
        queue_id = None
        if process.returncode == 0:
            try:
                # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
                queue_id = int(process.stdout.read().split('.')[0])
            except:
                # probably error parsing job code
                logger.critical("Could not parse job id following qsub...")

        return queue_id, process 

    def _get_njobs_in_queue(self, username):
        process = Popen(['qstat', '-a', '-u', username], stdout=PIPE, stderr=PIPE)
        process.wait()

        njobs = None
        if process.returncode == 0:
            # parse the result
            # lines should have this form
            # '1339044.sdb          username  queuename    2012-02-29-16-43  20460   --   --    --  00:20 C 00:09'
            # count lines that include the username in it

            # TODO: only count running or queued jobs. or rather, *don't* count jobs that are 'C'.
            outs = process.stdout.read().split('\n')
            njobs = len([line.split() for line in outs if username in line])

        return njobs, process

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
        """Set the memory per process in megabytes"""
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
    """Adapter for Sun Grid Engine (SGE) task submission software."""
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
        super(SGEAdapter, self).set_mpi_procs(mpi_procs)
        self.qparams["ncpus"] = mpi_procs

    def set_mem_per_proc(self, mem_mb):
        """Set the memory per process in megabytes"""
        super(SGEAdapter, self).set_mem_per_proc(mem_mb)
        # TODO
        #raise NotImplementedError("")
        #self.qparams["mem_per_cpu"] = mem_mb
        ## Remove mem if it's defined.
        #self.qparams.pop("mem", None)

    def cancel(self, job_id):
        return os.system("qdel %d" % job_id)

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        process = Popen(['qsub', script_file], stdout=PIPE, stderr=PIPE)
        process.wait()

        # grab the returncode. SGE returns 0 if the job was successful
        queue_id = None
        if process.returncode == 0:
            try:
                # output should of the form 
                # Your job 1659048 ("NAME_OF_JOB") has been submitted 
                queue_id = int(process.stdout.read().split(' ')[2])
            except:
                # probably error parsing job code
                logger.critical("Could not parse job id following qsub...")

        return queue_id, process

    def _get_njobs_in_queue(self, username):
        process = Popen(['qstat', '-u', username], stdout=PIPE, stderr=PIPE)
        process.wait()

        njobs = None
        if process.returncode == 0:
            # parse the result
            # lines should contain username
            # count lines that include the username in it

            # TODO: only count running or queued jobs. or rather, *don't* count jobs that are 'C'.
            outs = process.stdout.readlines()
            njobs = len([line.split() for line in outs if username in line])

        return njobs, process

    def exclude_nodes(self, nodes):
        """Method to exclude nodes in the calculation"""
        logger.warning('exluding nodes, not implemented yet in SGE')
        return False


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
        self.qparams["walltime"] = time2slurm(timelimit)

    def set_mem_per_proc(self, mem_mb):
        super(MOABAdapter, self).set_mem_per_proc(mem_mb)
        #TODO
        #raise NotImplementedError("set_mem_per_cpu")

    def exclude_nodes(self, nodes):
        logger.warning('exluding nodes, not implemented yet in MOAB')
        return False

    def cancel(self, job_id):
        return os.system("canceljob %d" % job_id)

    def _submit_to_queue(self, script_file):
        """Submit a job script to the queue."""
        process = Popen(['msub', script_file], stdout=PIPE, stderr=PIPE)
        process.wait()

        queue_id = None
        if process.returncode == 0:
            # grab the returncode. MOAB returns 0 if the job was successful
            try:
                # output should be the queue_id
                queue_id = int(process.stdout.read().split()[0])
            except:
                # probably error parsing job code
                logger.critical('Could not parse job id following msub...')

        return queue_id, process

    def _get_njobs_in_queue(self, username):
        process = Popen(['showq', '-s -u', username], stdout=PIPE, stderr=PIPE)
        process.wait()

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
            out = process.stdout.readlines()[-1]
            njobs = int(out.split()[-1])

        return njobs, process


class QScriptTemplate(string.Template):
    delimiter = '$$'
