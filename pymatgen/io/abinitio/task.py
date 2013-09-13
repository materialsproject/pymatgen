"""
Classes defining Abinit calculations and workflows
"""
from __future__ import division, print_function

import os
import shutil
import collections
import abc
import numpy as np

from pymatgen.core.design_patterns import Enum, AttrDict
from pymatgen.util.string_utils import stream_has_colours
from pymatgen.serializers.json_coders import MSONable, json_load, json_pretty_dump

from pymatgen.io.abinitio.utils import abinit_output_iscomplete, File
from pymatgen.io.abinitio.launcher import TaskLauncher
from pymatgen.io.abinitio.events import EventParser

#import logging
#logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

__all__ = [
    "task_factory",
]


class FakeProcess(object):
    """
    This object is attached to a Task instance if the task has not been submitted
    This trick allows us to simulate a process that is still running so that we can safely poll task.process.
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


def task_factory(strategy, workdir, runmode, task_id=1, links=None, **kwargs):
    """Factory function for Task instances."""
    # TODO
    # Instanciate subclasses depending on the runlevel.
    classes = {
       "scf": AbinitTask,
       "nscf": AbinitTask,
       "relax": AbinitTask,
       "dfpt": AbinitTask,
       "screening": AbinitTask,
       "sigma": AbinitTask,
       "bse": AbinitTask,
    }

    return classes[strategy.runlevel](strategy, workdir, runmode, task_id=task_id, links=links, **kwargs)


class TaskError(Exception):
    """Base Exception for Task methods"""


class Task(object):
    __metaclass__ = abc.ABCMeta

    Error = TaskError

    # Possible status of the task.
    S_READY = 1    # Task is ready for submission.
    S_SUB = 2      # Task has been submitted.
    S_RUN = 4      # Task is running.
    S_DONE = 8     # Task done, This does not imply that results are ok or that the calculation completed successfully
    S_ERROR = 16   # Task raised some kind of Error (either the process or ABINIT).
    S_OK = 32      # Task completed successfully.

    STATUS2STR = {
        S_READY: "ready",
        S_SUB: "submitted",
        S_RUN: "running",
        S_DONE: "done",
        S_ERROR: "error",
        S_OK: "completed",
    }

    def __init__(self):
        # Set the initial status.
        self.set_status(Task.S_READY)

        # Used to push additional info on the task during its execution. 
        #self._history = collections.OrderedDict()

    def __getstate__(self):
        """
        Return state is pickled as the contents for the instance.
                                                                                      
        In this case we just remove the process since Subprocess objects cannot be pickled.
        This is the reason why we have to store the returncode in self_returncode instead
        of using self.process.returncode.
        """
        return {k:v for k,v in self.__dict__.items() if k not in ["_process",]}

    @abc.abstractproperty
    def executable(self):
        """
        Path to the executable associated to the task (internally stored in self._executable).
        """

    def set_executable(self, executable):
        self._executable = executable

    # Interface modeled after subprocess.Popen
    @abc.abstractproperty
    def process(self):
        """Return an object that supports the `subprocess.Popen` protocol."""

    def poll(self):
        """Check if child process has terminated. Set and return returncode attribute."""
        self._returncode = self.process.poll()

        if self._returncode is not None:
            self.set_status(Task.S_DONE)

        return self._returncode

    def wait(self):
        """Wait for child process to terminate. Set and return returncode attribute."""
        self._returncode = self.process.wait()
        self.set_status(Task.S_DONE)

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
        self.set_status(Task.S_DONE)

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

    @property
    def id(self):
        """Task identifier."""
        return self._id

    def set_id(self, id):
        """Set the task identifier."""
        self._id = id

    @property
    def tot_ncpus(self):
        """Total number of CPUs used to run the task."""
        return self.launcher.tot_ncpus

    @property
    def status(self):
        """Gives the status of the task."""
        return self._status

    @property
    def str_status(self):
        """String representation of the status."""
        return self.STATUS2STR[self.status]

    def set_status(self, status):
        """Set the status of the task."""
        if status not in self.STATUS2STR:
            raise RuntimeError("Unknown status: %s" % status)

        self._status = status

        # Notify the observers
        #self.subject.notify_observers(self)

    def recheck_status(self):

        if self.status not in [self.S_SUB, self.S_RUN, self.S_DONE]: 
            return

        # Check the returncode of the process first.
        if self.returncode != 0:
            self.set_status(self.S_ERROR)

        else:
            # Check if the run completed successfully.
            report = EventParser().parse(self.output_file.path)
                                                                                                     
            if report.run_completed:
                self.set_status(self.S_OK)

            else:
                # TODO
                # This is the delicate part since we have to discern among different possibilities:
                #
                # 1) Calculation stopped due to an Abinit Error or Bug.
                # 2) Segmentation fault that (by definition) was not handled by ABINIT.
                # 3) Problem with the resource manager and/or the OS (walltime error, resource error, phase of the moon ...)
                # 4) Calculation is still running!
                #
                # Point 2) and 3) are the most complicated since there's no standard!

                # 1) Search for possible errors or bugs in the ABINIT **output** file.
                if report.errors or report.bugs:
                    self._status = self.S_ERROR
                    return 

                # 2) Analyze the stderr file for Fortran runtime errors.
                #   err_lines = self.stderr_file.readlines()
                #   if err_lines:
                #       self._status = self.S_ERROR
                #       return 

                # 3) Analyze the error file of the resource manager.
                #   self._status = self.S_ERROR
                #   err_lines = self.qerr_file.readlines()
                #   if err_lines:
                #       self._status = self.S_ERROR
                #       return 

                # 4)
                self._status = self.S_OK
                return 

    @property
    def links(self):
        return self._links

    @property
    def links_status(self):
        """Returns a list with the status of the links."""
        if not self.links:
            return [Task.S_DONE,]

        return [l.status for l in self.links]

    @property
    def is_allocated(self):
        """True if the task has been allocated, i.e. if it has been submitted or if it's running."""
        return self.status in [Task.S_SUB, Task.S_RUN]

    @property
    def is_completed(self):
        """True if the task has been executed."""
        return self.status >= Task.S_DONE

    @property
    def can_run(self):
        """The task can run if its status is S_READY and all the other links (if any) are done!"""
        #return (self.status == Task.S_READY) and all([link_stat==Task.S_DONE for link_stat in task.links_status])
        return (self.status == Task.S_READY) and all([stat >= Task.S_DONE for stat in self.links_status])

    def connect(self):
        """Create symbolic links to the output files produced by the other tasks."""
        for link in self.links:
            filepaths, exts = link.get_filepaths_and_exts()
            for (path, ext) in zip(filepaths, exts):
                dst = self.idata_path_from_ext(ext)

                if not os.path.exists(path): 
                    # Try netcdf file.
                    # TODO: this case should be treated in a cleaner way.
                    path += "-etsf.nc"
                    if os.path.exists(path):
                        dst += "-etsf.nc"

                if not os.path.exists(path):
                    raise self.Error("%s is needed by this task but it does not exist" % path)

                # Link path to dst
                print("Linking ", path," --> ",dst)
                os.symlink(path, dst)

    @abc.abstractmethod
    def setup(self, *args, **kwargs):
        """Public method called before submitting the task."""

    def _setup(self, *args, **kwargs):
        """
        This method calls self.setup after having performed additional operations
        such as the creation of the symbolic links needed to connect different tasks.
        """
        self.connect()
        self.setup(*args, **kwargs)

    def _parse_events(self):
        """Analyzes the main output for possible errors or warnings."""
        parser = EventParser()
        try:
            return parser.parse(self.output_file.path)

        except parser.Error as exc:
            raise
            # TODO: Handle possible errors in the parser by generating a custom EventList object
            #return EventList(self.output_file.path, events=[Error(str(exc))])

    @property
    def events(self):
        """List of errors or warnings reported by ABINIT."""
        if self.status is None or self.status < self.S_DONE:
            raise self.Error(
                "Task %s is not completed.\nYou cannot access its events now, use _parse_events" % repr(self))

        try:
            return self._events

        except AttributeError:
            self._events = self._parse_events()
            return self._events

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
            "task_events"    : self.events.to_dict
        })

    def move(self, dst, isabspath=False):
        """
        Recursively move self.workdir to another location. This is similar to the Unix "mv" command.
        The destination path must not already exist. If the destination already exists
        but is not a directory, it may be overwritten depending on os.rename() semantics.

        Be default, dst is located in the parent directory of self.workdir, use isabspath=True
        to specify an absolute path.
        """
        if not isabspath:
            dst = os.path.join(os.path.dirname(self.workdir), dst)

        shutil.move(self.workdir, dst)


class AbinitTask(Task):
    """
    Base class defining an ABINIT calculation
    """
    # Prefixes for Abinit (input, output, temporary) files.
    Prefix = collections.namedtuple("Prefix", "idata odata tdata")
    pj = os.path.join

    prefix = Prefix(pj("indata", "in"), pj("outdata", "out"), pj("tmpdata", "tmp"))
    del Prefix, pj

    # Basenames for Abinit input and output files.
    Basename = collections.namedtuple("Basename", "input output files_file log_file stderr_file jobfile lockfile")
    basename = Basename("run.abi", "run.abo", "run.files", "run.log", "stderr", "job.sh", "__abilock__")
    del Basename

    def __init__(self, strategy, workdir, runmode, task_id=1, links=None, **kwargs):
        """
        Args:
            strategy: 
                `Strategy` instance describing the calculation.
            workdir:
                Path to the working directory.
            runmode:
                `RunMode` instance or string "sequential"
            task_id:
                Task identifier (must be unique if self belongs to a `Workflow`).
            links:
                List of `WorkLink` objects specifying the dependencies of the task.
                Used for tasks belonging to a `Workflow`.
            kwargs:
                keyword arguments (not used here)
        """
        super(AbinitTask, self).__init__()

        self.strategy = strategy

        self.workdir = os.path.abspath(workdir)

        self.runmode = RunMode.asrunmode(runmode).copy()

        self.set_id(task_id)

        # Connect this task to the other tasks 
        # (needed if are creating a Work instance with dependencies
        self._links = []
        if links is not None: 
            if not isinstance(links, collections.Iterable):
                links = [links,]

            self._links = links
            for link in links:
                print("Adding ", link.get_abivars())
                self.strategy.add_extra_abivars(link.get_abivars())

        # Files required for the execution.
        self.input_file = File(os.path.join(self.workdir, AbinitTask.basename.input))
        self.output_file = File(os.path.join(self.workdir, AbinitTask.basename.output))
        self.files_file = File(os.path.join(self.workdir, AbinitTask.basename.files_file))
        self.log_file = File(os.path.join(self.workdir, AbinitTask.basename.log_file))
        self.stderr_file = File(os.path.join(self.workdir, AbinitTask.basename.stderr_file))

        # Build the launcher and store it.
        self.set_launcher(self.runmode.make_launcher(self))

    def __repr__(self):
        return "<%s at %s, workdir = %s>" % (
            self.__class__.__name__, id(self), os.path.basename(self.workdir))

    def __str__(self):
        return self.__repr__()

    @classmethod
    def from_input(cls, abinit_input, workdir, runmode, task_id=1, links=None, **kwargs):
        """
        Create an instance of `AbinitTask` from an ABINIT input.

        Args:
            abinit_input:
                `AbinitInput` object.
            workdir:
                Path to the working directory.
            runmode:
                `RunMode` instance or string "sequential"
            task_id:
                Task identifier (must be unique if self belongs to a `Workflow`).
            links:
                List of `WorkLink` objects specifying the dependencies of the task.
                Used for tasks belonging to a `Workflow`.
            kwargs:
                keyword arguments (not used here)
        """
        # TODO: Find a better way to do this. I will likely need to refactor the Strategy object
        from pymatgen.io.abinitio.strategies import StrategyWithInput
        strategy = StrategyWithInput(abinit_input)
        return cls(strategy, workdir, runmode, task_id=task_id, links=links, **kwargs)

    @property
    def executable(self):
        try:
            return self._executable
        except AttributeError:
            return "abinit"

    @property
    def name(self):
        return self.workdir

    @property
    def process(self):
        try:
            return self._process
        except AttributeError:
            # Attach a fake process so that we can poll it.
            return FakeProcess()

    @property
    def launcher(self):
        return self._launcher

    def set_launcher(self, launcher):
        self._launcher = launcher

    @property
    def jobfile_path(self):
        """Absolute path of the job file (shell script)."""
        return os.path.join(self.workdir, self.basename.jobfile)    

    @property
    def indata_dir(self):
        """Directory with the input data."""
        head, tail = os.path.split(self.prefix.idata)
        return os.path.join(self.workdir, head)

    @property
    def outdata_dir(self):
        """Directory with the output data."""
        head, tail = os.path.split(self.prefix.odata)
        return os.path.join(self.workdir, head)

    @property
    def tmpdata_dir(self):
        """Directory with the temporary data."""
        head, tail = os.path.split(self.prefix.tdata)
        return os.path.join(self.workdir, head)

    def idata_path_from_ext(self, ext):
        """Returns the path of the input file with extension ext."""
        return os.path.join(self.workdir, self.prefix.idata + ext)

    def odata_path_from_ext(self, ext):
        """Returns the path of the output file with extension ext"""
        return os.path.join(self.workdir, self.prefix.odata + ext)

    @property
    def pseudos(self):
        """List of pseudos used in the calculation."""
        return self.strategy.pseudos

    @property
    def filesfile_string(self):
        """String with the list of files and prefixex needed to execute ABINIT."""
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
        raise NotImplementedError("")
        d = {k: v.to_dict for k, v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["strategy"] = self.strategy 
        d["workdir"] = workdir 
        return d
                                                                    
    @staticmethod
    def from_dict(d):
        raise NotImplementedError("")

    def isdone(self):
        """True if the output file is complete."""
        return abinit_output_iscomplete(self.output_file.path)

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return all(p.isnc for p in self.pseudos)

    @property
    def ispaw(self):
        """True if PAW calculation"""
        return all(p.ispaw for p in self.pseudos)

    def make_input(self):
        """Construct and write the input file of the calculation."""
        return self.strategy.make_input()

    def outfiles(self):
        """Return all the output data files produced."""
        files = []
        for fname in os.listdir(self.outdata_dir):
            if fname.startswith(os.path.basename(self.prefix.odata)):
                files.append(os.path.join(self.outdata_dir, fname))

        return files
                                                                  
    def tmpfiles(self):
        """Return all the input data files produced."""
        files = []
        for fname in os.listdir(self.tmpdata_dir):
            if file.startswith(os.path.basename(self.prefix.tdata)):
                files.append(os.path.join(self.tmpdata_dir, fname))

        return files

    def path_in_workdir(self, filename):
        """Create the absolute path of filename in the top-level working directory."""
        return os.path.join(self.workdir, filename)

    def path_in_indatadir(self, filename):
        """Create the absolute path of filename in the indata directory."""
        return os.path.join(self.indata_dir, filename)

    def path_in_outdatadir(self, filename):
        """Create the absolute path of filename in the outdata directory."""
        return os.path.join(self.outdata_dir, filename)

    def path_in_tmpdatadir(self, filename):
        """Create the absolute path of filename in the tmp directory."""
        return os.path.join(self.tmpdata_dir, filename)

    def rename(self, src_basename, dest_basename, datadir="outdir"):
        """
        Rename a file located in datadir.

        src_basename and dest_basename are the basename of the source file
        and of the destination file, respectively.
        """
        if datadir == "outdir":
            src = self.path_in_outdatadir(src_basename)
            dest = self.path_in_outdatadir(dest_basename)

        elif datadir == "tmpdir":
            src = self.path_in_tmpdatadir(src_basename)
            dest = self.path_in_tmpdatadir(dest_basename)

        else:
            raise ValueError("Wrong datadir %s" % datadir)

        os.rename(src, dest)

    def build(self, *args, **kwargs):
        """
        Creates the working directory and the input files of the `Task`.
        It does not overwrite files if they already exist.
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
            self.input_file.write(self.make_input())

        if not self.files_file.exists:
            self.files_file.write(self.filesfile_string)

        #if not os.path.exists(self.jobfile_path):
        #    with open(self.jobfile_path, "w") as f:
        #            f.write(str(self.jobfile))

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
            pats = exclude_wildcard.split("|")

            import fnmatch

            def keep(fname):
                """True if fname must be preserved."""
                for pat in pats:
                    if fnmatch.fnmatch(fname, pat):
                        return True
                return False

            for dirpath, dirnames, filenames in os.walk(self.workdir):
                for fname in filenames:
                    filepath = os.path.join(dirpath, fname)
                    if not keep(fname):
                        os.remove(filepath)

    def remove_files(self, *filenames):
        """Remove all the files listed in filenames."""
        if isinstance(filenames, str):
            filenames = [filenames]

        for dirpath, dirnames, fnames in os.walk(self.workdir):
            for fname in fnames:
                if fname in filenames:
                    filepath = os.path.join(dirpath, fname)
                    os.remove(filepath)

    def rm_tmpdatadir(self):
        """Remove the directory with the temporary files."""
        shutil.rmtree(self.tmpdata_dir)

    def rm_indatadir(self):
        """Remove the directory with the input files (indata dir)."""
        shutil.rmtree(self.indata_dir)

    def setup(self, *args, **kwargs):
        pass

    def start(self, *args, **kwargs):
        """
        Starts the calculation by performing the following steps:

            - build dirs and files
            - call the _setup method
            - execute the job file via the shell or other means.
            - return Results object
        """
        self.build(*args, **kwargs)

        self._setup(*args, **kwargs)

        # Start the calculation in a subprocess and return.
        self._process = self.launcher.launch(self)

    def start_and_wait(self, *args, **kwargs):
        """Helper method to start the task and wait."""
        self.start(*args, **kwargs)
        return self.wait()


class TaskResults(dict, MSONable):
    """
    Dictionary used to store the most important results produced by a Task.
    """
    _MANDATORY_KEYS = [
        "task_name",
        "task_returncode",
        "task_status",
        "task_events",
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
            assert (self["task_returncode"] == 0 and self["task_status"] == self.S_DONE)
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
        return cls.from_dict(json_load(filename))    


class RunHintsError(Exception):
    """Base error class for `RunHints`."""


class RunHints(collections.OrderedDict):
    """
    Dictionary with the hints for the parallel execution reported by ABINIT.

    Example:
        <RUN_HINTS, max_ncpus = "108", autoparal="3">
            "1" : {"tot_ncpus": 2,         # Total number of CPUs
                   "mpi_ncpus": 2,         # Number of MPI processes.
                   "omp_ncpus": 1,         # Number of OMP threads (0 if not used)
                   "memory_gb": 10,        # Estimated memory requirement in Gigabytes (total or per proc?)
                   "weight"   : 0.4,       # 1.0 corresponds to an "expected" optimal efficiency.
                   "abivars": {            # Dictionary with the variables that should be added to the input.
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
    Error = RunHintsError

    def __init__(self, *args, **kwargs):
        super(RunHints, self).__init__(*args, **kwargs)

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
            raise self.Error("%s does not contain any valid RUN_HINTS section" % filename)

        runhints = json_loads("".join(l for l in lines[start+1:stop]))

        # Sort the dict so that configuration with minimum (weight-1.0) appears in the first positions
        return cls(sorted(runhints, key=lambda t : abs(t[1]["weight"] - 1.0), reverse=False))

    def select(self, policy):
        """Find the optimal configuration according on policy."""
        raise NotImplementedError("")
        #okey = ??
        #self["mpi_ncpus"] = min(mpi_ncpus, self["max_ncpus"])
        #self["omp_ncpus"] = min(omp_ncpus, self["max_ncpus"])

        # Copy the dictionary and return AttrDict instance.
        return AttrDict(self[okey].copy())


class RunMode(dict, MSONable):
    """
    This object contains the user-specified settings controlling the execution 
    of the run (submission, coarse- and fine-grained parallelism ...)
    It acts as a factory for `Launcher` objects.

    A `TaskManager` is responsible for the generation and the submission 
    script, as well as of the specification of the parameters passed to the Resource Manager
    (e.g. Slurm, PBS ...) and/or the dynamic specification of the ABINIT variables governing the 
    parallel execution. A `TaskManager` delegates the generation of the submission
    script and the submission of the task to the `QueueAdapter`. 
    A `TaskManager` has a `TaskPolicy` object that governs the specification of the 
    parameters for the parallel executions.
    """
    # Default values (correspond to the sequential mode).
    _DEFAULTS = {
        "launcher_type": "shell",    # ["shell", "slurm", "pbs",]
        "use_fw"       : False,      # True if we are using fireworks.
        "policy"       : "default",  # Policy used to select the number of MPI processes when there are ambiguities.
        "max_ncpus"    : np.inf,     # Max number of CPUs that can be used (DEFAULT: no limit is enforced)
    }

    def __init__(self, *args, **kwargs):
        super(RunMode, self).__init__(*args, **kwargs)

    def copy(self):
        return RunMode(**super(RunMode, self).copy())
    
    @classmethod
    def asrunmode(cls, obj):
        """
        Returns a `RunMode` instance from obj. Accepts:

            - `RunMode` instances
            - A string with "sequential" to indicate that Abinit will be executed in sequential mode.
        """
        if isinstance(obj, cls):
            return obj

        if isinstance(obj, str) and obj == "sequential":
            return cls.sequential()

        raise ValueError("Don't know how to convert obj %s into a %s instance" % (obj, obj.__class__.__name__))

    @classmethod
    def sequential(cls, launcher_type=None):
        d = cls._DEFAULTS.copy()
        if launcher_type is not None:
            d["launcher_type"] = launcher_type

        return cls(d)

    @classmethod
    def mpi_parallel(cls, mpi_ncpus, launcher_type=None):
        d = cls._DEFAULTS.copy()
        if launcher_type is not None:
            d["launcher_type"] = launcher_type
        new = cls(d)
        new.set_mpi_ncpus(mpi_ncpus)
        return new

    def set_mpi_ncpus(self, mpi_ncpus):
        self["mpi_ncpus"] = min(mpi_ncpus, self["max_ncpus"])

    #def set_omp_ncpus(self, omp_ncpus):
    #     self["omp_ncpus"] = min(omp_ncpus, self["max_ncpus"])

    @property
    def to_dict(self):
        d = self.copy()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d
                                                 
    @classmethod
    def from_dict(cls, d):
        return cls({k: v for (k, v) in d.items() if k not in ("@module", "@class")})

    def make_launcher(self, task):
        """Factory function returning a launcher instance."""

        #if self.use_fw:
        #   raise NotImplementedError("Firework not yet supported")
                                                       
        # Find the subclass to instanciate.
        cls = TaskLauncher.from_launcher_type(task.runmode["launcher_type"])

        return cls(task.jobfile_path, 
                   stdin  = task.files_file.path, 
                   stdout = task.log_file.path,
                   stderr = task.stderr_file.path,
                   **task.runmode
                  )


class TaskPolicy(dict):
    """
    This object stores the parameters used by the `TaskManager` to 
    create the submission script and/or the modification of the 
    ABINIT variables governing the parallel execution.
    """

    # Default values.
    _DEFAULTS = dict(
        use_fw=False,      # True if we are using fireworks.
        max_ncpus=np.inf,  # Max number of CPUs that can be used (DEFAULT: no limit is enforced)
    )

    def __init__(self, *args, **kwargs):
        super(TaskPolicy, self).__init__(*args, **kwargs)

        err_msg = ""
        for k, v in self:
            if k not in self._DEFAULTS:
                err_msg += "Unknow key %s\n" % k

        if err_msg:
            raise ValueError(err_msg)

        for k, v in self._DEFAULTS:
            if k not in self:
                self[k] = v

    @classmethod
    def default_policy(cls):
        return cls(**cls._DEFAULTS.copy())


class TaskManager(object):
    """
    A `TaskManager` is responsible for the generation and the submission 
    script, as well as of the specification of the parameters passed to the Resource Manager
    (e.g. Slurm, PBS ...) and/or the dynamic specification of the ABINIT variables governing the 
    parallel execution. A `TaskManager` delegates the generation of the submission
    script and the submission of the task to the `QueueAdapter`. 
    A `TaskManager` has a `TaskPolicy` object that governs the specification of the 
    parameters for the parallel executions.
    """

    def __init__(self, qadapter_class, policy=None):
        self.qadapter = qadapter_class()
        self.policy = policy if policy is not None else TaskPolicy.default_policy()

    #@classmethod 
    #def simple_manager(cls, mpi_ncpus=1):
    #    """
    #    Simplest task manager that submits jobs with a simple shell script.
    #    Assume the shell environment is already properly initialized.
    #    """
    #    return cls(qadapter_class=ShellAdapter, policy=None)

    def optimize_params(self, task):
        """
        Find optimal parameters for the execution of the task 
        using the options specified in self.policy.
     
        This method can change the ABINIT input variables and/or the 
        parameters passed to the Resource Manager e.g. the number of CPUs.
        """
        raise NotImplementedError("")

        # 1) Run ABINIT in sequential to get the possible configurations.
        retcode = self._dummy_run(task)

        if retcode != 0:
            return 

        # 2) Parse the configurations
        try:
            hints = RunHints.from_file(task.log_file.path)

        except RunHints.Error:
            raise

        # 3) Select the optimal configuration according to policy
        optimal = hints.select(self.policy)

        # 4) Change the input file and/or the submission script
        task.strategy.add_abivars(optimal.abivars)

        #self.qadapter.set_mpi_nprocs(optimal.mpi_nprocs)
        #self.qadapter.set_omp_nprocs(optimal.omp_nprocs)

    def _dummy_run(self, task):
        retcode = 0
        return retcode

    def submit_task(self, task):
        """Submit the task."""

        #if self.policy.use_fw:
        #   raise NotImplementedError("Firework not yet supported")
                                                       
        # Find the subclass to instanciate.
        #cls = TaskLauncher.from_launcher_type(task.runmode["launcher_type"])
        #return cls(task.jobfile_path, 
        #           stdin  = task.files_file.path, 
        #           stdout = task.log_file.path,
        #           stderr = task.stderr_file.path,
        #           **task.runmode
        #          )

        # NEW VERSION based of queue_adapter.
        #qad = self.qadapter
        #script = qad.get_script_str(
        #    job_name=task.name, 
        #    launch_dir=task.workdir, 
        #    executable=task.executable,
        #    stdin =task.files_file.path, 
        #    stdout=task.log_file.path,
        #    stderr=task.stderr_file.path,
        #)

        #script_file = task.job_file.path
        #with open(script_file, "w") as fh:
        #    fh.write(script)

        #process = qad.submit_to_queue(script_file)
        #return process
