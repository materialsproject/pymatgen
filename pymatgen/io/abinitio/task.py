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
import json

from pymatgen.core.design_patterns import Enum
from pymatgen.util.string_utils import stream_has_colours
#from pymatgen.util.filelock import FileLock

from .utils import parse_ewc, abinit_output_iscomplete, File
from .jobfile import JobFile

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
class BaseTask(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def process(self):
        "Return an object that support the subprocess.Popen protocol."

    # Interface modeled after subprocess.Popen
    def poll(self):
        "Check if child process has terminated. Set and return returncode attribute."
        return self.process.poll()

    def wait(self):
        "Wait for child process to terminate. Set and return returncode attribute."
        return self.process.wait()

    def communicate(self, input=None):
        """
        Interact with process: Send data to stdin. Read data from stdout and stderr, until end-of-file is reached. 
        Wait for process to terminate. The optional input argument should be a string to be sent to the 
        child process, or None, if no data should be sent to the child.

        communicate() returns a tuple (stdoutdata, stderrdata).
        """
        return self.process.communicate(input=input)

    @property
    def returncode(self):
        """
        The child return code, set by poll() and wait() (and indirectly by communicate()). 
        A None value indicates that the process hasn't terminated yet.
        A negative value -N indicates that the child was terminated by signal N (Unix only).
        """
        return self.process.returncode

    def kill(self):
        "Kill the child"
        self.process.kill()

    @abc.abstractmethod
    def setup(self, *args, **kwargs):
        """Method called before submitting the task."""

    def _setup(self, *args, **kwargs):
        self.setup(*args, **kwargs)

    def get_results(self, *args, **kwargs):
        """
        Method called once the calculation is completed, returns TaskResults instance.
        Subclasses should extend this method (if needed) by adding 
        specialized code that perform some kind of post-processing.
        """
        # Check whether the process completed.
        if self.returncode is None:
            raise ValueError("return code is None, you should call wait, communitate or poll")

        return TaskResults({
            "task_name": self.name,
            "task_returncode" : self.returncode,
            #"task_status": self.get_status()
            })

##########################################################################################

class AbinitTask(BaseTask):
    """
    Base class defining a generic abinit calculation
    """
    # Prefixes for Abinit (input, output, temporary) files.
    Prefix = collections.namedtuple("Prefix", "idata odata tdata")
    pj = os.path.join

    prefix = Prefix("in", pj("output","out"), pj("temporary","tmp"))
    del Prefix, pj

    # Basenames for Abinit input and output files.
    Basename = collections.namedtuple("Basename", "input output files_file log_file stderr_file jobfile lockfile")
    basename = Basename("run.input", "run.output", "run.files", "log", "stderr", "job.sh", "__lock__")
    del Basename

    def __init__(self, input, workdir, 
                 varpaths = None, 
                 **kwargs
                ):
        """
        Args:
            input: 
                AbinitInput instance.
            workdir:
                Path to the working directory.
            varpaths:
        """
        self.workdir = os.path.abspath(workdir)

        self.input = input.copy()

        self.need_files = []
        if varpaths is not None: 
            self.need_files.extend(varpaths.values())
            self.input = self.input.add_vars_to_control(varpaths)

        # Files required for the execution.
        self.input_file  = File(self.basename.input      , dirname=self.workdir)
        self.output_file = File(self.basename.output     , dirname=self.workdir)
        self.files_file  = File(self.basename.files_file , dirname=self.workdir)
        self.log_file    = File(self.basename.log_file   , dirname=self.workdir)
        self.stderr_file = File(self.basename.stderr_file, dirname=self.workdir)

        # Find number of processors ....
        #self.paral_hint = self.get_runhints(max_numproc)

        # Set jobfile variables.
        # TODO Rename JobFile to ShellScript(File)
        #self.jobfile = File(self.workdir, self.basename.jobfile)
        #excecutable = "abinit"
        self.jobfile = JobFile(name   = self.jobfile_path, 
                               input  = self.files_file.path, 
                               log    = self.log_file.path,
                               stderr = self.stderr_file.path,
                              )
    def __str__(self):
        lines = []
        app = lines.append

        #basenames = ["input", "files_file", "jobfile"]
        #for s in basenames:
        #    apath = getattr(self, s + "_path") 
        #    app(s + ": " + apath)

        app("outfiles_dir: " + self.outfiles_dir)
        app("tmpfiles_dir: " + self.tmpfiles_dir)

        return "\n".join(lines)

    def __repr__(self):
        return "<%s at %s, task_workdir = %s>" % (
            self.__class__.__name__, id(self), os.path.basename(self.workdir))

    @property
    def name(self):
        return self.workdir

    @property
    def process(self):
        return self._process

    @property
    def jobfile_path(self):
        "Absolute path of the job file (shell script)."
        return os.path.join(self.workdir, self.basename.jobfile)    

    @property
    def outfiles_dir(self):
        head, tail = os.path.split(self.prefix.odata)
        return os.path.join(self.workdir, head)

    @property
    def tmpfiles_dir(self):
        head, tail = os.path.split(self.prefix.tdata)
        return os.path.join(self.workdir, head)

    def odata_path_from_ext(self, ext):
        "Return the path of the output file with extension ext"
        if ext[0] != "_": ext = "_" + ext
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
        return AbinitTask(**d)

    def iscomplete(self):
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
        for file in os.listdir(self.outfiles_dir):
            if file.startswith(os.path.basename(self.prefix.odata)):
                files.append(os.path.join(self.outfiles_dir, file))
        return files
                                                                  
    def tmpfiles(self):
        "Return all the input data files produced."
        files = list()
        for file in os.listdir(self.tmpfiles_dir):
            if file.startswith(os.path.basename(self.prefix.tdata)):
                files.append(os.path.join(self.tmpfiles_dir, file))
        return files

    def path_in_workdir(self, filename):
        "Create the absolute path of filename in the workind directory."
        return os.path.join(self.workdir, filename)

    def build(self, *args, **kwargs):
        """
        Writes Abinit input files and directory.
        Do not overwrite files if they already exist.
        """
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        if not os.path.exists(self.outfiles_dir):
            os.makedirs(self.outfiles_dir)

        if not os.path.exists(self.tmpfiles_dir):
            os.makedirs(self.tmpfiles_dir)

        if not self.input_file.exists:
            self.input_file.write(str(self.input))

        if not self.files_file.exists:
            self.files_file.write(self.filesfile_string)

        if not os.path.exists(self.jobfile_path):
            with open(self.jobfile_path, "w") as f:
                    f.write(str(self.jobfile))

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

    def read_mainlog_ewc(self, nafter=5):
       """
       Read errors, warnings and comments from the main output and the log file.
                                                                                        
       :return: Two namedtuple instances: main, log. 
                The lists of strings with the corresponding messages are
                available in main.errors, main.warnings, main.comments, log.errors etc.
       """
       main = parse_ewc(self.output_file.path, nafter=nafter)
       log  = parse_ewc(self.log_file.path, nafter=nafter)

       return main, log

    def get_status(self):
        return Status(self)

    def show_status(self, stream=sys.stdout, *args, **kwargs):
        self.get_status().show(stream=stream, *args, **kwargs)

    def get_runhints(self, max_numproc):
        """
        Run abinit in sequential to obtain a set of possible
        configuration for the number of processors
        """
        raise NotImplementedError("")

        # number of MPI processors, number of OpenMP processors, memory estimate in Gb.
        RunHints = collections.namedtuple('RunHints', "mpi_nproc omp_nproc memory_gb")

        return hints

    def setup(self, *args, **kwargs):
        pass

    def _submit(self, *args, **kwargs):
        # Start the calculation in a subprocess and return
        from subprocess import Popen, PIPE
        self._process = Popen((self.jobfile.shell, self.jobfile.path), cwd=self.workdir, stderr=PIPE)

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
        if kwargs.pop('verbose', 0):
            print('Running ' + self.input_path)

        self.build(*args, **kwargs)

        self._setup(*args, **kwargs)

        self._submit(*args, **kwargs)

##########################################################################################

import cPickle as pickle
class Pickleable(object):

    def write_to_pickle_file(self, filename, protocol=-1):
        """
        Writes the pickle representation to a file.

        Args:
            filename:
                filename to write to. It is recommended that the file extension be ".pickle".
        """
        with open(filename, "w") as f:
            pickle.dump(self, filename)

    def from_pickle_file(filename):
        """
        Return the object from the pickle representation stored in a file.
                                                                                              
        Args:
            filename:
                filename to read. 
        """
        with open(filename, "r") as f:
            return pickle.load(filename)
#TODO
from pymatgen.serializers.json_coders import MSONable #, PMGJSONDecoder

class TaskResults(dict, MSONable):
    """
    Dictionary used to store some of the results produce by a Task object
    """
    _mandatory_keys = [
        "task_name",
        "task_returncode",
        #"status",
    ]

    def __init__(self, *args, **kwargs):
        super(TaskResults, self).__init__(*args, **kwargs)

        for k in self._mandatory_keys:
            if k not in self:
                raise ValueError("%s must be specified" % k)

    def __getattribute__(self, name):
        try:
            # Default behaviour
            return super(TaskResults, self).__getattribute__(name)
        except AttributeError:
            try:
                # Try in the dictionary.
                return self[name]
            except KeyError as exc:
                raise AttributeError(str(exc))

    def __setitem__(self, key, val):
        "Add key-value pair to self."
        if key in self:
            try:
                if val == self[key]:
                    return
            except:
                raise ValueError("key %s already exists and cannot check for equality" % key)

        super(TaskResult, self).__setitem__(key, val)

    def update(self, *args, **kwargs):
        for (k, v) in dict(*args, **kwargs).iteritems():
            self[k] = v

    @property
    def to_dict(self):
        d = self
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        my_dict = {k: v for k in d if k not in ["@module", "@class",]}
        return cls(**my_dict)

    #def json_load(self, filename):
    #    with open(filename, "r") as fh
    #       return Results.from_dict(json.load(fh))

    #def json_dump(self, filename):
    #    with open(filename, "w") as fh
    #        json.dump(self.to_dict, fh, indent=4, sort_keys=4)
    
##########################################################################################

class Status(object):
    """
    Object used to inquire the status of the task and to have access 
    to the error and warning messages
    """
    S_OK   = 1 
    S_WAIT = 2 
    S_RUN  = 4
    S_ERR  = 8

    #: Rank associated to the different status (in order of increasing priority)
    _level2rank = {
      "completed" : S_OK, 
      "unstarted" : S_WAIT, 
      "running"   : S_RUN, 
      "error"     : S_ERR,
    }

    levels = Enum(s for s in _level2rank)

    def __init__(self, task):

        self._task = task

        if task.iscomplete():
            # The main output file seems completed.
            level_str = 'completed'

            #TODO
            # Read the number of errors, warnings and comments 
            # for the (last) main output and the log file.
            #main, log =  task.read_mainlog_ewc()
                                                                
            #main_info = main.tostream(stream)
            #log_info  = log.tostream(stream)

        # TODO err file!
        elif task.output_file.exists and task.log_file.exists:
            # Inspect the main output and the log file for ERROR messages.
            main, log = task.read_mainlog_ewc()

            level_str = 'running'
            if main.errors or log.errors:
                level_str = 'error'

            with open(task.stderr_file.path, "r") as f:
                lines = f.readlines()
                if lines: level_str = 'error'
        else:
            level_str = 'unstarted'

        self._level_str = level_str

        assert self._level_str in Status.levels

    @property
    def rank(self):
        return self._level2rank[self._level_str]

    @property
    def is_completed(self):
        return self._level_str == "completed"

    @property
    def is_unstarted(self):
        return self._level_str == "unstarted"

    @property
    def is_running(self):
        return self._level_str == "running"

    @property
    def is_error(self):
        return self._level_str == "error"

    def __repr__(self):
        return "<%s at %s, task = %s, level = %s>" % (
            self.__class__.__name__, id(self), repr(self._task), self._level_str)

    def __str__(self):
        return str(self._task) + self._level_str

    # Rich comparison support. Mainly used for selecting the 
    # most critical status when we have several tasks executed inside a workflow.
    def __lt__(self, other):
        return self.rank < other.rank

    def __le__(self, other):
        return self.rank <= other.rank

    def __eq__(self, other):
        return self.rank == other.rank

    def __ne__(self, other):
        return not self == other

    def __gt__(self, other):
        return self.rank > other.rank

    def __ge__(self, other):
        return self.rank >= other.rank

    def show(self, stream=sys.stdout, *args, **kwargs):

        from ..tools import StringColorizer
        str_colorizer = StringColorizer(stream)
                                                                       
        _status2txtcolor = {
          "completed" : lambda string : str_colorizer(string, "green"),
          "error"     : lambda string : str_colorizer(string, "red"),
          "running"   : lambda string : str_colorizer(string, "blue"),
          "unstarted" : lambda string : str_colorizer(string, "cyan"),
        }

        color = lambda stat : str
        if stream_has_colours(stream):
            color = lambda stat : _status2txtcolor[stat](stat)

        stream.write(self._task.name + ': ' + color(self._level_str) + "\n")

    def print_ewc_messages(self, stream=sys.stdout, **kwargs):

        verbose = kwargs.pop("verbose", 0)

        # Read the number of errors, warnings and comments 
        # for the (last) main output and the log file.
        main, log =  self._task.read_mainlog_ewc()
                                                            
        main_info = main.tostream(stream)
        log_info  = log.tostream(stream)
                                                            
        stream.write("\n".join([main_info, log_info]) + "\n")

        if verbose:
            for w in main.warnings: 
                stream.write(w + "\n")
            if verbose > 1:
                for w in log.warnings: 
                    stream.write(w + "\n")

        if self.is_error:
            for e in log.errors: 
                    stream.write(e + "\n")

            if not log.errors: # Read the standard error.
                with open(self._task.stderr_file.path, "r") as f:
                    lines = f.readlines()
                    stream.write("\n".join(lines))

##########################################################################################

class TaskDependencies(object):
    """
    This object describes the dependencies among Task instances.
    """

    _ext2getvars = {
        "DEN" : "getden_path",
        "WFK" : "getwfk_path",
        "SCR" : "getscr_path",
        "QPS" : "getqps_path",
    }

    def __init__(self, task, task_id, *odata_required):
        """
        Args:
            task: 
                AbinitTask instance
            task_id: 
                Task identifier 
            odata_required:
                Output data required for running the task.
        """
        self.task    = task
        self.task_id = task_id

        self.odata_required = []
        if odata_required:
            self.odata_required = odata_required

    def with_odata(self, *odata_required):
        return TaskDependencies(self.task, self.task_id, *odata_required)

    def get_varpaths(self):
        vars = {}
        for ext in self.odata_required:
            varname = self._ext2getvars[ext]
            path = self.task.odata_path_from_ext(ext)
            vars.update({varname : path})
        return vars


##########################################################################################
class TaskLauncher(object): 
    __metaclass__ = abc.ABCMeta

    def __init__(self, paral_info=None, env=None, modules=None, pre_cmds=None, post_cmds=None):

        self.paral_info = paral_info if paral_info is not None else {}
        self.env = env if env is not None else {}
        self.modules = modules if modules is not None else []
        self.pre_cmds = pre_cmds if pre_cmds is not None else []
        self.post_cmds = post_cmds if post_cmds is not None else []

    @abc.abstractmethod
    def write_jobfile(self, task):
        "Write the submission script"

    @abc.abstractmethod
    def submit_task(self, task):
        "Submit a task, set task._process (process-like object)"

    def launch_task(self, task):
        self.write_jobfile(task)
        self.submit_task(task)

##########################################################################################

class ShellLauncher(TaskLauncher): 

    def submit_task(self, task):
        # Start the calculation in a subprocess and return
        from subprocess import Popen, PIPE
        task._process = Popen((task.jobfile.shell, task.jobfile.path), cwd=task.workdir, stderr=PIPE)

##########################################################################################

class SlurmLauncher(TaskLauncher): 
    pass

##########################################################################################
