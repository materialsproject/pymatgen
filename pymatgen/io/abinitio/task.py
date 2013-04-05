"""
Classes defining Abinit calculations and workflows
"""
from __future__ import division, print_function

import sys
import os
import os.path
import shutil
import collections
import abc
import numpy as np

from pymatgen.core.design_patterns import Enum, AttrDict
from pymatgen.util.string_utils import stream_has_colours
from pymatgen.serializers.json_coders import MSONable, json_load, json_pretty_dump #, PMGJSONDecoder

from pymatgen.io.abinitio.utils import abinit_output_iscomplete, File
from pymatgen.io.abinitio.launcher import TaskLauncher
from pymatgen.io.abinitio.events import EventParser

#import logging
#logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

__all__ = [
"AbinitTask",
]

##########################################################################################
class FakeProcess(object):
    def poll(self):
        return None

    def wait(self):
        raise RuntimeError("Cannot wait a Fake Process")

    def communicate(self, input=None):
        raise RuntimeError("Cannot communicate with Fake Process")

    def kill(self):
        raise RuntimeError("Cannot kill a Fake Process")

    @property
    def returncode(self):
        return None

##########################################################################################

class TaskError(Exception):
    pass

class Task(object):
    __metaclass__ = abc.ABCMeta

    Error = TaskError

    # Possible status of the task.
    # Rank associated to the different status (in order of increasing priority)
    S_DONE  = 1  # Task completed, note that this does not imply that results are ok or that the calculation completed succesfully
    S_READY = 2 # Task is ready for submission.
    S_SUB   = 4  # Task has been submitted.
    S_RUN   = 8  # Task is running.

    _status2str = {
      S_DONE : "done",
      S_READY: "waiting",
      S_SUB  : "submitted", 
      S_RUN  : "running", 
    }

    def __init__(self):
        self.set_status(Task.S_READY)

    # Interface modeled after subprocess.Popen
    @abc.abstractproperty
    def process(self):
        "Return an object that supports the subprocess.Popen protocol."

    def poll(self):
        "Check if child process has terminated. Set and return returncode attribute."
        returncode = self.process.poll()
        print("In task poll returncode %s" % returncode)
        if returncode is not None:
            print("setting status for task %s" % self)
            self.set_status(Task.S_DONE)
        return returncode

    def wait(self):
        "Wait for child process to terminate. Set and return returncode attribute."
        returncode = self.process.wait()
        self.set_status(Task.S_DONE)
        return returncode

    def communicate(self, input=None):
        """
        Interact with process: Send data to stdin. Read data from stdout and stderr, until end-of-file is reached. 
        Wait for process to terminate. The optional input argument should be a string to be sent to the 
        child process, or None, if no data should be sent to the child.

        communicate() returns a tuple (stdoutdata, stderrdata).
        """
        stdoutdata, stderrdata = self.process.communicate(input=input)
        self.set_status(Task.S_DONE)
        return stdoutdata, stderrdata 

    def kill(self):
        "Kill the child"
        self.process.kill()

    @property
    def returncode(self):
        """
        The child return code, set by poll() and wait() (and indirectly by communicate()). 
        A None value indicates that the process hasn't terminated yet.
        A negative value -N indicates that the child was terminated by signal N (Unix only).
        """
        return self.process.returncode

    @property
    def id(self):
        "Task identifier"
        return self._id

    def set_id(self, id):
        "Set the task identifier"
        self._id = id

    @property
    def tot_ncpus(self):
        "Total number of CPUs used to run the task."
        return self.launcher.tot_ncpus

    @property
    def status(self):
        return self._status

    def set_status(self, status):
        if status not in Task._status2str:
            raise RunTimeError("Unknown status: %s" % status)
        self._status = status

        # Notify the observers
        #self.subject.notify_observers(status)

    @property
    def links_status(self):
        "Returns a list wit the status of the links"
        if not self._links:
            return [Task.S_DONE,]
        return [l.status for l in self._links]

    @abc.abstractmethod
    def setup(self, *args, **kwargs):
        """Method called before submitting the task."""

    def _setup(self, *args, **kwargs):
        #
        # Create symbolic links to the output files produced by the other tasks
        for link in self._links:
            filepaths, exts = link.get_filepaths_and_exts()
            for (path, ext) in zip(filepaths, exts):
                dst = self.idata_path_from_ext(ext)
                print("Will make %s point to %s" % (dst, path))

                if not os.path.exists(path):
                    raise self.Error("%s is needed by this task by it does not exist" % path)
                else:
                    # Link path to dst
                    os.symlink(path, dst)

        self.setup(*args, **kwargs)

    @property
    def results(self):
        "The results produced by the task. Set by get_results"
        try:
            return self._results
        except AttributeError:
            return None

    def get_results(self, *args, **kwargs):
        """
        Method called once the calculation is completed, 
        Updates self._results and returns TaskResults instance.
        Subclasses should extend this method (if needed) by adding 
        specialized code that perform some kind of post-processing.
        """
        # Check whether the process completed.
        if self.returncode is None:
            raise self.Error("return code is None, you should call wait, communitate or poll")

        if self.status != Task.S_DONE:
            raise self.Error("Task is not completed")

        # Analyze the main output and the log file for possible errors or warnings.
        main_events = EventParser().parse(self.output_file.path)

        results = TaskResults({
            "task_name"      : self.name,
            "task_returncode": self.returncode,
            "task_status"    : self.status,
            "task_events"    : main_events.to_dict
            })

        if not hasattr(self, "_results"):
            self._results = TaskResults()

        self._results.update(results)
        return results

##########################################################################################

class AbinitTask(Task):
    """
    Base class defining an abinit calculation
    """
    # Prefixes for Abinit (input, output, temporary) files.
    Prefix = collections.namedtuple("Prefix", "idata odata tdata")
    pj = os.path.join

    prefix = Prefix(pj("indata", "in"), pj("outdata","out"), pj("tmpdata","tmp"))
    del Prefix, pj

    # Basenames for Abinit input and output files.
    Basename = collections.namedtuple("Basename", "input output files_file log_file stderr_file jobfile lockfile")
    basename = Basename("run.in", "run.out", "run.files", "log", "stderr", "job.sh", "__lock__")
    del Basename

    def __init__(self, input, workdir, runmode, 
                 task_id  = 1,
                 links    = None, 
                 **kwargs
                ):
        """
        Args:
            input: 
                AbinitInput instance.
            workdir:
                Path to the working directory.
            runmode:
                RunMode instance.
            links:
                List of WorkLink objects specifying the dependencies of self.
        """
        super(AbinitTask, self).__init__()

        self.workdir = os.path.abspath(workdir)

        self.runmode = runmode

        self.set_id(task_id)

        self.input = input.copy()

        # Connect this task to the other tasks 
        # (needed if are creating a Work instance with dependencies
        self._links = []
        if links is not None: 
            self._links = links
            for link in links:
                abivars = link.get_abivars()
                print("abivars", abivars)
                self.input = self.input.add_vars_to_control(abivars)

        # Files required for the execution.
        self.input_file  = File(os.path.join(self.workdir, AbinitTask.basename.input))
        self.output_file = File(os.path.join(self.workdir, AbinitTask.basename.output))
        self.files_file  = File(os.path.join(self.workdir, AbinitTask.basename.files_file))
        self.log_file    = File(os.path.join(self.workdir, AbinitTask.basename.log_file))
        self.stderr_file = File(os.path.join(self.workdir, AbinitTask.basename.stderr_file))

        # Find number of processors ....
        #runhints = self.get_runhints()
        self.set_launcher(TaskLauncher.from_task(self))

    #def __str__(self):

    def __repr__(self):
        return "<%s at %s, workdir = %s>" % (
            self.__class__.__name__, id(self), os.path.basename(self.workdir))

    @property
    def name(self):
        return self.workdir

    @property
    def process(self):
        try:
            return self._process
        except AttributeError:
            return FakeProcess()

    @property
    def launcher(self):
        return self._launcher

    def set_launcher(self, launcher):
        self._launcher = launcher

    @property
    def jobfile_path(self):
        "Absolute path of the job file (shell script)."
        return os.path.join(self.workdir, self.basename.jobfile)    

    @property
    def indata_dir(self):
        head, tail = os.path.split(self.prefix.idata)
        return os.path.join(self.workdir, head)

    @property
    def outdata_dir(self):
        head, tail = os.path.split(self.prefix.odata)
        return os.path.join(self.workdir, head)

    @property
    def tmpdata_dir(self):
        head, tail = os.path.split(self.prefix.tdata)
        return os.path.join(self.workdir, head)

    def idata_path_from_ext(self, ext):
        "Return the path of the input file with extension ext"
        return os.path.join(self.workdir, (self.prefix.idata + ext))

    def odata_path_from_ext(self, ext):
        "Return the path of the output file with extension ext"
        return os.path.join(self.workdir, (self.prefix.odata + ext))

    @property
    def pseudos(self):
        "List of pseudos used in the calculation."
        return self.input.pseudos

    @property
    def filesfile_string(self):
        "String with the list of files and prefixes needed to execute abinit."
        lines = []
        app = lines.append
        pj = os.path.join

        app(self.input_file.path)                 # Path to the input file
        app(self.output_file.path)                # Path to the output file
        app(pj(self.workdir, self.prefix.idata))  # Prefix for in data
        app(pj(self.workdir, self.prefix.odata))  # Prefix for out data
        app(pj(self.workdir, self.prefix.tdata))  # Prefix for tmp data

        # Paths to the pseudopotential files.
        for pseudo in self.pseudos:
            app(pseudo.path)

        return "\n".join(lines)

    @property
    def to_dict(self):
        raise NotimplementedError("")
        d = {k: v.to_dict for k, v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["input"] = self.input 
        d["workdir"] = workdir 
        return d
                                                                    
    @staticmethod
    def from_dict(d):
        raise NotimplementedError("")

    def isdone(self):
        "True if the output file is complete."
        return abinit_output_iscomplete(self.output_file.path)

    @property
    def isnc(self):
        "True if norm-conserving calculation"
        return all(p.isnc for p in self.pseudos)

    @property
    def ispaw(self):
        "True if PAW calculation"
        return all(p.ispaw for p in self.pseudos)

    def outfiles(self):
        "Return all the output data files produced."
        files = list()
        for file in os.listdir(self.outdata_dir):
            if file.startswith(os.path.basename(self.prefix.odata)):
                files.append(os.path.join(self.outdata_dir, file))
        return files
                                                                  
    def tmpfiles(self):
        "Return all the input data files produced."
        files = list()
        for file in os.listdir(self.tmpdata_dir):
            if file.startswith(os.path.basename(self.prefix.tdata)):
                files.append(os.path.join(self.tmpdata_dir, file))
        return files

    def path_in_workdir(self, filename):
        "Create the absolute path of filename in the top-level working directory."
        return os.path.join(self.workdir, filename)

    def path_in_indatadir(self, filename):
        "Create the absolute path of filename in the indata directory."
        return os.path.join(self.indata_dir, filename)

    def path_in_outdatadir(self, filename):
        "Create the absolute path of filename in the outdata directory."
        return os.path.join(self.outdata_dir, filename)

    def build(self, *args, **kwargs):
        """
        Writes Abinit input files and directory.
        Do not overwrite files if they already exist.
        """
        # Top-level dir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # Create dirs for input, output and tmp data.
        if not os.path.exists(self.indata_dir):
            os.makedirs(self.indata_dir)

        if not os.path.exists(self.outdata_dir):
            os.makedirs(self.outdata_dir)

        if not os.path.exists(self.tmpdata_dir):
            os.makedirs(self.tmpdata_dir)

        # Write input file and files file.
        if not self.input_file.exists:
            self.input_file.write(str(self.input))

        if not self.files_file.exists:
            self.files_file.write(self.filesfile_string)

        #if not os.path.exists(self.jobfile_path):
        #    with open(self.jobfile_path, "w") as f:
        #            f.write(str(self.jobfile))

    def rmtree(self, *args, **kwargs):
        """
        Remove all files and directories.
                                                                                   
        Keyword arguments:
            force: (False)
                Do not ask confirmation.
            verbose: (0)
                Print message if verbose is not zero.
        """
        if kwargs.pop('verbose', 0):
            print('Removing directory tree: %s' % self.workdir)

        shutil.rmtree(self.workdir)

    def get_runhints(self):
        """
        Run abinit in sequential to obtain a set of possible
        configuration for the number of processors.
        RunMode is used to select the paramenters of the run.
        """
        raise NotImplementedError("")

        # number of MPI processors, number of OpenMP processors, memory estimate in Gb.
        return runhints

    def setup(self, *args, **kwargs):
        pass

    def start(self, *args, **kwargs):
        """
        Start the calculation by performing the following steps: 
            - build dirs and files
            - call the _setup method
            - execute the job file via the shell or other means.
            - return Results object

        Keyword arguments:
            verbose: (0)
                Print message if verbose is not zero.

        .. warning::
             This method must be thread safe since we may want to run several indipendent
             calculations with different python threads. 
        """
        self.build(*args, **kwargs)

        self._setup(*args, **kwargs)

        # Start the calculation in a subprocess and return
        self._process = self.launcher.launch(self)

        self.set_status(Task.S_RUN)

##########################################################################################

class TaskResults(AttrDict, MSONable):
    """
    Dictionary used to store some of the results produce by a Task object
    """
    _mandatory_keys = [
        "task_name",
        "task_returncode",
        "task_status",
        "task_events",
    ]

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

    def json_dump(self, filename):
        json_pretty_dump(self.to_dict, filename) 

    @classmethod
    def json_load(cls, filename):
        return cls.from_dict(json_load(filename))    

##########################################################################################

class RunMode(AttrDict, MSONable):
    """
    User-specified options controlling the execution of the run (submission, 
    coarse- and fine-grained parallelism ...)
    """
    _defaults = {
        "launcher"      : "shell",    # ["shell", "slurm", "slurm-fw", "pbs", "pbs-fw",]
        "policy"        : "default",  # Policy used to select the number of MPI processes when there are ambiguities.
        "mpirun"        : "mpirun",   # path to the mpirunnner.
        "max_npcus"     : 0,          # Max number of cpus that can be used. If 0, no maximum limit is enforced
        "omp_numthreads": 0,          # Number of OpenMP threads, 0 if OMP is not used 
        "chunk_size"    : 1,          # Used when we don't have a queue_manager and several jobs to run
                                      # In this case we run at most chunk_size tasks when work.start is called.
        "queue_params"  : {},         # Parameters passed to the firework QueueAdapter.
    }

    @staticmethod
    def load_user_configuration():
        fname = "runmode.json"
        cwd = os.getcwd()

        path = os.path.join(cwd, fname)
        if os.path.exists(path):
            return RunMode.from_file(path)
        else:
            raise NotImplementedError("Try in the pymatgen dir")
        #else:
        #    return RunMode.sequential()

    @classmethod
    def from_defaults(cls):
        return cls(cls._defaults.copy())

    @classmethod
    def sequential(cls, chunk_size=1, launcher=None):
        d = cls.from_defaults()
        d["chunk_size"] = chunk_size
        if launcher is not None:
            d["launcher"] = launcher
        return cls(d)

    @classmethod
    def from_file(cls, filename):
        """
        Initialize an instance of RunMode from the configuration file filename (JSON format)
        """
        defaults = RunMode._defaults.copy() 
        d = json_load(filename) 

        # Put default values if they are not in d.
        for (k,v) in defaults.items():
            if k not in d:
                d[k] = v
        return cls(d)

    @property
    def has_queue_manager(self):
        "True if job submission is handled by a resource manager."
        return self.launcher not in ["shell",]

    def get_chunk_size(self):
        if self.has_queue_manager:
            return None
        return self.chunk_size

    @property
    def to_dict(self):
        d = self.copy()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d
                                                 
    @classmethod
    def from_dict(cls, d):
        return RunMode({k:v for (k, v) in d.items() if k not in ("@module", "@class")})

    #def with_max_ncpus(self, max_ncpus):
    #    d = RunMode._defaults.copy()
    #    d["max_ncpus"] = 1
    #    return RunMode(**d)

##########################################################################################

class RunHints(collections.OrderedDict):
    """
    Dictionary with the hints for the parallel execution reported by abinit.

    Example:
        <RUN_HINTS, max_ncpus = "108", autoparal="3">
            "1" : {"tot_ncpus": 2,           # Total number of CPUs
                   "mpi_ncpus": 2,           # Number of MPI processes.
                   "omp_ncpus": 1,           # Number of OMP threads (0 if not used)
                   "memory_gb": 10,          # Estimated memory requirement in Gigabytes (total or per proc?)
                   "weight"   : 0.4,         # 1.0 corresponds to an "expected" optimal efficiency.
                   "variables": {            # Dictionary with the variables that should be added to the input.
                        "varname1": varvalue1,
                        "varname2": varvalue2,
                        }
                }
        </RUN_HINTS>

    For paral_kgb we have:
    nproc     npkpt  npspinor    npband     npfft    bandpp    weight   
       108       1         1        12         9         2        0.25
       108       1         1       108         1         2       27.00
        96       1         1        24         4         1        1.50
        84       1         1        12         7         2        0.25
    """
    @classmethod
    def from_file(cls, filename):
        """
        Read the <RUN_HINTS> section from file filename
        Assumes the file contains only one section.
        """
        with open(filaname, "r") as fh:
            lines = fh.readlines()

        try:
            start = lines.index("<RUN_HINTS>\n")
            stop  = lines.index("</RUN_HINTS>\n")
        except ValueError:
            raise ValueError("%s does not contain any valid RUN_HINTS section" % filename)

        runhints = json_loads("".join(l for l in lines[start+1:stop]))

        # Sort the dict so that configuration with minimum (weight-1.0) appears in the first positions
        return cls(sorted(runhints, key=lambda t : abs(t[1]["weight"] - 1.0), reverse=False))

    def __init__(self, *args, **kwargs):
        super(RunHints, self).__init__(*args, **kwargs)

    def get_hint(self, policy):
        # Find the optimal configuration depending on policy
        raise NotImplementedError("")

        # Copy the dictionary and return AttrDict instance.
        odict = self[okey].copy()
        odict["_policy"] = policy

        return AttrDict(odict)

##########################################################################################

class RunParams(AttrDict):
    """
    Dictionary with the parameters governing the execution of the task.
    In contains the parameters passed by the user via RunMode and the 
    parameteres suggested by abinit.
    """
    def __init__(self, *args, **kwargs):
        super(TaskRunParams, self).__init__(*args, **kwargs)

##########################################################################################
