"""
Classes defining Abinit calculations and workflows
"""
from __future__ import division, print_function

import sys
import os
import os.path
import collections
import subprocess
import numpy as np

from pymatgen.core.design_patterns import Enum
from pymatgen.util.string_utils import stream_has_colours
from pymatgen.util.filelock import FileLock

from .utils import parse_ewc, abinit_output_iscomplete
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
# Helper functions.

def _remove_tree(*paths, **kwargs):
    import shutil
    ignore_errors = kwargs.pop("ignore_errors", True)
    onerror       = kwargs.pop("onerror", None)
                                                                           
    for path in paths:
        shutil.rmtree(path, ignore_errors=ignore_errors, onerror=onerror)

class File(object):
    """
    Very simple class used to store file basenames, absolute paths and directory names.

    Provides wrappers for the most commonly used os.path functions.
    """
    def __init__(self, basename, dirname="."):
        self.basename = basename
        self.dirname = os.path.abspath(dirname)
        self.path = os.path.join(self.dirname, self.basename)

    def __str__(self):
       return self.read()

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    @property
    def exists(self):
        "True if file exists."
        return os.path.exists(self.path)

    @property
    def isncfile(self):
        "True if self is a NetCDF file"
        return self.basename.endswith(".nc")

    def read(self):
        with open(self.path, "r") as f: 
            return f.read()

    def readlines(self):
        with open(self.path, "r") as f: 
            return f.readlines()

    def write(self, string):
        with open(self.path, "w") as f: 
            return f.write(string)
                                        
    def writelines(self, lines):
        with open(self.path, "w") as f: 
            return f.writelines()

##########################################################################################

class AbinitTask(object):

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
                 varpaths    = None, 
                 user_options= None,
                ):
        """
        Args:
            input: 
                AbinitInput instance.
            workdir:
                Name of the working directory.
            varpaths:
            user_options:
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
        return self.input_file.basename

    @property
    def return_code(self):
        try:
            return self._return_code
        except AttributeError:
            return 666

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
        return all(p.isnc for p in self.pseudos)

    @property
    def ispaw(self):
        return all(p.ispaw for p in self.pseudos)

    #def in_files(self):
    #    "Return all the input data files."
    #    files = list()
    #    for file in os.listdir(dirname(self.idat_root)):
    #        if file.startswith(basename(self.idat_root)):
    #            files.append(join(dirname(self.idat_root), file))
    #    return files
                                                                  
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

    def destroy(self, *args, **kwargs):
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

        _remove_tree(self.workdir, **kwargs)

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

    def setup(self, *args, **kwargs):
        """
        Method called before running the calculations. 
        The default implementation creates the workspace directory and the input files.
        """

    def teardown(self, *args, **kwargs):
        """
        Method called once the calculations completed.
        The default implementation does nothing.
        """

    def get_runhints(self, max_numproc):
        """
        Run abinit in sequential to obtain a set of possible
        configuration for the number of processors
        """
        raise NotImplementedError("")

        # number of MPI processors, number of OpenMP processors, memory estimate in Gb.
        RunHints = collections.namedtuple('RunHints', "mpi_nproc omp_nproc memory_gb")

        return hints

    def run(self, *args, **kwargs):
        """
        Run the calculation by executing the job file from the shell.

        Keyword arguments:
            verbose: (0)
                Print message if verbose is not zero.

        .. warning::
             This method must be thread safe since we may want to run several indipendent
             calculations with different python threads. 
        """
        if kwargs.get('verbose'):
            print('Running ' + self.input_path)

        lock = FileLock(self.path_in_workdir(AbinitTask.basename.lockfile))

        try:
            lock.acquire()
        except FileLock.Exception:
            raise

        self.setup(*args, **kwargs)

        self._return_code = subprocess.call((self.jobfile.shell, self.jobfile.path), cwd=self.workdir)

        if self._return_code == 0:
            self.teardown(*args, **kwargs)

        lock.release()

        return self._return_code

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

