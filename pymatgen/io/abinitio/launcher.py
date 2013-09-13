#!/usr/bin/env python
from __future__ import division, print_function

import os 
import abc

from subprocess import Popen, PIPE

from pymatgen.io.abinitio.utils import File


__all__ = [
    "ScriptEditor",
    "ShellLauncher",
]


class ScriptEditor(object):
    """Simple editor that simplifies the writing of shell scripts"""
    _shell = '/bin/bash'

    def __init__(self):
        self._lines = []

    @property
    def shell(self):
        return self._shell

    def _add(self, text, pre=""):
        if isinstance(text, str):
            self._lines.append(pre+text)
        else:
            self._lines.extend([pre + t for t in text])

    def reset(self):
        """Reset the editor."""
        try:
            del self._lines
        except AttributeError:
            pass

    def shebang(self):
        """Adds the shebang line."""
        self._lines.append('#!' + self.shell)

    def declare_var(self, key, val):
        """Declare a env variable. If val is None the variable is unset."""
        if val is not None:
            line = "export " + key + '=' + str(val)
        else:
            line = "unset " + key 

        self._add(line)

    def declare_vars(self, d):
        """Declare the variables defined in the dictionary d."""
        for k,v in d.items():
            self.declare_var(k, v)

    def export_envar(self, key, val):
        """Export an environment variable."""
        line = "export " + key + "=" + str(val)
        self._add(line)

    def export_envars(self, env):
        """Export the environment variables contained in the dict env."""
        for k,v in env.items():
            self.export_envar(k, v)

    def add_emptyline(self):
        self._add("", pre="")

    def add_comment(self, comment):
        self._add(comment, pre="# ")

    def load_modules(self, modules):
        for module in modules:
            self.load_module(module)

    def load_module(self, module):
        self._add('module load ' + module)

    def add_line(self, line):
        self._add(line)

    def add_lines(self, lines):
        self._add(lines)

    def get_script_str(self, reset=True):
        """Returns a string with the script and reset the editor if reset is True"""
        s = "\n".join(l for l in self._lines)
        if reset:
            self.reset()
        return s


class OmpEnv(dict):
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

    def __init__(self, *args, **kwargs):
        """ 
        Constructor method inherited from dictionary:
                                                   
        >>> OmpEnv(OMP_NUM_THREADS=1)
        {'OMP_NUM_THREADS': '1'}
                                                   
        To create an instance from an INI file, use:
           OmpEnv.from_file(filename)
        """
        self.update(*args, **kwargs)

        err_msg = ""
        for key, value in self.items():
            self[key] = str(value)
            if key not in self._KEYS:
                err_msg += "unknown option %s\n" % key

        if err_msg: 
            raise ValueError(err_msg)

    @staticmethod
    def from_file(filename, allow_empty=False):
        from ConfigParser import SafeConfigParser, NoOptionError
        parser = SafeConfigParser()
        parser.read(filename)

        obj = OMPEnv()

        # Consistency check. Note that we only check if the option name is correct, 
        # we do not check whether the value is correct or not.
        if "openmp" not in parser.sections():
            if not allow_empty:
                raise ValueError("%s does not contain any [openmp] section" % filename) 
            return obj

        err_msg = ""
        for key in parser.options("openmp"):
            if key.upper() not in self._KEYS:
                err_msg += "unknown option %s, maybe a typo" % key

        if err_msg: 
            raise ValueError(err_msg)

        for key in self._KEYS:
            try:
                obj[key] = str(parser.get("openmp", key))
            except NoOptionError:
                try:
                    obj[key] = str(parser.get("openmp", key.lower()))
                except NoOptionError:
                    pass

        if not allow_empty and not obj:
            raise ValueError("Refusing to return with an empty dict") 

        return obj


class TaskLauncher(object):
    """Abstract class for a task launcher."""
    __metaclass__ = abc.ABCMeta

    @staticmethod
    def from_launcher_type(launcher_type):
        """Returns the subclass from a string giving its type."""
        classes = []
        for cls in TaskLauncher.__subclasses__():
           if cls._type == launcher_type:
               classes.append(cls)
        if len(classes) != 1:
            raise ValueError("Cannot find class from launcher_type %s" % launcher_type)

        return classes[0]

    def __init__(self, path, **kwargs):
        """
        Args:
            path:
                Path to the script file.
            stdin:
                The input file to feed in the executable as the standard input. Mandatory.
            stdout:
                The file into which the standard output is redirected. Default is 'log'.
            stderr:
                The file into which the standard error is redirected. Default is 'stderr'.
            bindir:
                The directory in which to look for binaries. Default is "".
            modules:
                The modules which will be loaded with 'module load'. Default is None
            pre_lines:
                Lines before the main execution. Default is []
            post_lines:
                Lines after the main execution. Default is [].
            exe:
                The name of the binary to be executed 
                (located in bindir or in $PATH if not bindir). Default is abinit.
            mpirun:
                The mpi runner.
                E.g. 'mpirun -n 6'. Default is None.
            mpi_ncpus
                Number of MPI processes (default 1)
            omp_env
                Dictionary with the OMP environment variables. Default is {}
            vars:
                Dictionary {varname:varvalue} with variable declaration.
        """
        self.jobfile    = File(path)
        self.stdin      = kwargs.get("stdin",  "abi.files")
        self.stdout     = kwargs.get("stdout", "abi.log")
        self.stderr     = kwargs.get("stderr", "abi.err")
        self.bindir     = kwargs.get("bindir", "")
        self.exe        = os.path.join(self.bindir, kwargs.get("exe", "abinit"))
        self.mpirun     = kwargs.get("mpirun", "")
        self.mpi_ncpus  = kwargs.get("mpi_ncpus", 1)
        self.omp_env    = kwargs.get("omp_env", {})
        if self.omp_env:
            self.omp_env = OMPEnv(self.omp_env)
        self.vars         = kwargs.get("vars", {})
        self.modules      = kwargs.get("modules", [])
        self.envars       = kwargs.get("envars", {})
        self.pre_lines    = kwargs.get("pre_lines", [])
        self.post_lines   = kwargs.get("post_lines", [])
        self.queue_params = kwargs.get("queue_params", {})

    @property
    def tot_ncpus(self):
        """Total number of CPUs employed"""
        return self.mpi_ncpus * self.omp_ncpus 

    @property
    def has_mpirun(self):
        """True if we are using a mpirunner"""
        return bool(self.mpirun)

    @property
    def has_omp(self):
        """True if we are using OMP threads"""
        return hasattr(self,"omp_env") and bool(getattr(self, "omp_env"))
                                                      
    @property
    def omp_ncpus(self):
        """Number of OMP threads."""
        if self.has_omp:
            return self.omp_env["OMP_NUM_THREADS"]
        else:
            return 1

    def get_script_str(self):
        """Returns a string wth the script."""
        se = ScriptEditor()

        se.shebang()

        # Submission instructions for the queue manager.
        se.add_lines(self.make_rsmheader())

        se.add_comment("Variable declarations")
        vars = self.vars.copy()

        vars.update({
                "EXECUTABLE": self.exe,
                "STDIN"     : self.stdin,
                "STDOUT"    : self.stdout,
                "STDERR"    : self.stderr,
                "MPIRUN"    : self.mpirun,
                "MPI_NCPUS" : self.mpi_ncpus,
               })

        se.declare_vars(vars)

        se.add_comment("Modules")
        se.load_modules(self.modules)

        se.add_comment("Environment")
        se.export_envars(self.envars)

        if self.has_omp:
            se.export_envars(self.omp_env)

        se.add_comment("Commands before execution")
        se.add_lines(self.pre_lines)

        # Execution line
        # TODO: better treatment of mpirun syntax.
        if self.mpirun:
            se.add_line('$MPIRUN -n $MPI_NCPUS $EXECUTABLE < $STDIN > $STDOUT 2> $STDERR')
        else:
            se.add_line('$EXECUTABLE < $STDIN > $STDOUT 2> $STDERR')

        se.add_comment("Commands after execution")
        se.add_lines(self.post_lines)

        return se.get_script_str()

    def write_script(self, overwrite=True):
        """Writes the script file."""
        #if not self.jobfile.exists:
        #    os.makedirs(self.jobfile.dirname)

        #if self.jobfile.exists and not overwrite:
        #    raise ValueError("%s already exists, cannot overwrite" % self.jobfile.path)

        with open(self.jobfile.path, 'w') as f:
            f.write(self.get_script_str())

    # TODO
    #@abc.abstractproperty
    #def has_resource_manager(self):
    #    "True if job submission is handled by a resource manager."
    #    #return self.launcher not in ["shell",]

    @abc.abstractmethod
    def make_rsmheader(self):
        """Return a list of string with the options passed to the resource manager"""

    @abc.abstractmethod
    def launch(self, task, *args, **kwargs):
        """
        Run the script in a subprocess, returns a Popen-like object.

        Args:
            task:
                Task instance. launch should call task.set_status to modify the status 
                of the object. 
                task.set_status(Task.S_SUB) is called before submitting the task
                task.set_status(Task.S_RUN) is called once the task has started to run.
        """

##########################################################################################
# Concrete classes

class ShellLauncher(TaskLauncher):
    _type = "shell"

    def make_rsmheader(self):
        return []

    def launch(self, task, *args, **kwargs):
        self.write_script()

        task.set_status(task.S_SUB)
        process = Popen(("/bin/bash", self.jobfile.path), cwd=self.jobfile.dirname, stderr=PIPE)
        task.set_status(task.S_RUN)

        return process


class ResourceManagerError(Exception):
    """Base error class for `SimpleResourceManager."""


class SimpleResourceManager(object):

    Error = ResourceManagerError

    def __init__(self, work, max_ncpus, sleep_time=20, verbose=0):
        """
            Args:
                work:
                    Work instance.
                max_ncpus: 
                    The maximum number of CPUs that can be used.
                sleep_time:
                    Time delay (seconds) before trying to start a new task.
        """
        #from pymatgen.io.abinitio.workflow import Workflow
        #if not isinstance(work, Workflow):
        #    work = Workflow.from_task(work)

        self.work = work

        self.max_ncpus = max_ncpus 
        self.sleep_time = sleep_time
        self.verbose = verbose

        for task in self.work:
            if task.tot_ncpus > self.max_ncpus:
                err_msg = "Task %s requires %s CPUs, but max_ncpus is %d" % (
                    repr(task), task.tot_ncpus, max_ncpus)
                raise self.Error(err_msg)

    def run(self, *args, **kwargs):
        """Call the start method of the object contained in work."""
        while True:
            polls = self.work.poll()
            # Fetch the first task that is ready to run
            try:
                task = self.work.fetch_task_to_run()
            except StopIteration:
                break

            if task is None:
                import time
                time.sleep(self.sleep_time)
            else:
                # Check that we don't exceed the number of cpus employed, before starting.
                if self.verbose:
                    print("work polls %s" % polls)
                    print("work status %s" % self.work.get_status())

                if (task.tot_ncpus + self.work.ncpus_reserved <= self.max_ncpus): 
                    if self.verbose: print("Starting task %s" % task)
                    task.start()

        # Wait until all tasks are completed.
        self.work.wait()

        #if any([t.status != t.S_DONE for t in self]):
        #   for task in self:
        #       print([link for link in task._links])
        #       #print([link_stat==task.S_DONE for link_stat in task.links_status])
        #       msg = "Deadlock, likely due to task dependencies: status %s" % str([t.status for t in self]))
        #       self.Error(msg)

        return self.work.returncodes

