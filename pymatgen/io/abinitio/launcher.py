"""Tools for the submission of Tasks."""
from __future__ import division, print_function

import os 
import time

from subprocess import Popen, PIPE
from pymatgen.util.string_utils import is_string
from pymatgen.io.abinitio.utils import File

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "ScriptEditor",
    "PyLauncher",
]

class ScriptEditor(object):
    """
    Simple editor that simplifies the writing of shell scripts
    """
    _shell = '/bin/bash'

    def __init__(self):
        self._lines = []

    @property
    def shell(self):
        return self._shell

    def _add(self, text, pre=""):
        if is_string(text):
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


class ResourceManagerError(Exception):
    """Base error class for `SimpleResourceManager."""


class PyResourceManager(object):

    Error = ResourceManagerError

    def __init__(self, work, max_ncpus, sleep_time=20):
        """
        This object submits the tasks contained in a `Workflow`
        inside an infinite loop. Therefore it is mainly used 
        when we are running in interative mode and we have to 
        execute several small jobs without having to pass throuh
        the queue resource manager. Its main goal is organizing
        the execution of the tasks so that we don't exceed the 
        computing resources at hand.

        Args:
            work:
                `Workflow` instance.
            max_ncpus: 
                The maximum number of CPUs that can be used at any given time.
            sleep_time:
                Time delay (seconds) before trying to start a new `Task`.
        """
        self.work = work

        self.max_ncpus = max_ncpus 
        self.sleep_time = sleep_time

        for task in self.work:
            if task.tot_ncpus > self.max_ncpus:
                err_msg = "Task %s requires %s CPUs, but max_ncpus is %d" % (
                    repr(task), task.tot_ncpus, max_ncpus)
                raise self.Error(err_msg)

    def run(self, *args, **kwargs):
        """Call the start method of the object contained in work."""

        if self.work.all_done:
            return self.work.returncodes

        while True:
            polls = self.work.poll()
            # Fetch the first task that is ready to run
            try:
                task = self.work.fetch_task_to_run()
            except StopIteration:
                break

            if task is None:

                time.sleep(self.sleep_time)
                self.work.check_status()

            else:
                # Check that we don't exceed the number of cpus employed, before starting.
                logger.info("work polls %s" % polls)
                logger.info("work status %s" % self.work.get_all_status())

                if (task.tot_ncpus + self.work.ncpus_allocated <= self.max_ncpus): 
                    logger.info("Starting task %s" % task)
                    task.start()

        # Wait until all tasks are completed.
        self.work.wait()

        return self.work.returncodes


class PyLauncherError(Exception):
    """Error class for PyLauncher."""

class PyLauncher(object):
    """
    This object handle the sumbmission of the `Task` in a `Workflow`
    """

    Error = PyLauncherError

    def __init__(self, work):
        """
        Initialize the object

        Args:
            work:
                `Workflow` object
        """
        self.work = work

    def single_shot(self):
        """
        Run the fist `Task than is ready for execution.
        
        Returns:
            Number of jobs launched.
        """
        work = self.work
        nlaunch = 0

        try:
            task = work.fetch_task_to_run()

            if task is None:
                logger.debug("No task to run! Possible deadlock")
                                                            
            else:
                logger.debug("Starting task %s" % task)
                task.start()
                nlaunch += 1
                                                            
        except StopIteration:
            logger.debug("Out of tasks.")

        finally:
            work.pickle_dump()

        return nlaunch 

    def rapidfire(self):
        """
        Keep on submitting `Tasks` until we are out of jobs or 
        no job is ready to run.
        
        Returns:
            The number of tasks launched.
        """
        nlaunch, launched = 0, []

        work = self.work
        while True:
            try:
                task = work.fetch_task_to_run()
                #work.check_status()

                if task is None:
                    logger.debug("fetch_task_to_run returned None.")
                    break

                logger.debug("Starting task %s" % task)
                assert task not in launched
                launched.append(task)
                task.start()

                nlaunch += 1

            except StopIteration:
                logger.debug("Out of tasks.")
                break

        work.pickle_dump()

        return nlaunch 
