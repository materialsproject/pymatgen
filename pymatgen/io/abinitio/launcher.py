"""Tools for the submission of Tasks."""
from __future__ import division, print_function

import os 
import time
import yaml

from subprocess import Popen, PIPE
from pymatgen.core.design_patterns import AttrDict
from pymatgen.util.string_utils import is_string
from pymatgen.io.abinitio.utils import File

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "ScriptEditor",
    "PyLauncher",
    "PyFlowsScheduler",
]

def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


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
        """Add an empty line."""
        self._add("", pre="")

    def add_comment(self, comment):
        """Add a comment"""
        self._add(comment, pre="# ")

    def load_modules(self, modules):
        """Load the list of specified modules."""
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
        """Reads the OpenMP variables from a INI file."""
        if filename.endswith(".ini"):
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

        else:
            raise NotImplementedError("Don't how how to read data from %s" % filename)


class PyLauncherError(Exception):
    """Error class for PyLauncher."""


class PyLauncher(object):
    """
    This object handle the submission of the tasks contained in a `AbinitFlow`
    """
    Error = PyLauncherError

    def __init__(self, flow):
        """
        Initialize the object

        Args:
            flow:
                `AbinitFlow` object
        """
        self.flow = flow

    def single_shot(self):
        """
        Run the first `Task` than is ready for execution.
        
        Returns:
            Number of jobs launched.
        """
        num_launched = 0

        # Get the tasks that can be executed in each workflow.
        tasks = []
        for work in self.flow:
            try:
                task = work.fetch_task_to_run()

                if task is not None:
                    tasks.append(task)
                else:
                    logger.debug("No task to run! Possible deadlock")

            except StopIteration:
                logger.info("Out of tasks.")

        # Submit the tasks and update the database.
        if tasks:
            tasks[0].start()
            num_launched += 1
            
            self.flow.pickle_dump()

        return num_launched 

    def rapidfire(self, max_nlaunch=-1):  # nlaunches=0, max_loops=-1, sleep_time=None
        """
        Keeps submitting `Tasks` until we are out of jobs or no job is ready to run.

        Args:
            max_nlaunch:
                Maximum number of launches. default: no limit.

        nlaunches: 
            0 means 'until completion', -1 or "infinite" means to loop forever
        max_loops: 
            maximum number of loops
        sleep_time: 
            secs to sleep between rapidfire loop iterations
        
        Returns:
            The number of tasks launched.
        """
        #sleep_time = sleep_time if sleep_time else FWConfig().RAPIDFIRE_SLEEP_SECS
        #nlaunches = -1 if nlaunches == 'infinite' else int(nlaunches)
        #num_loops = 0
        num_launched, launched = 0, []

        for work in self.flow:
            if num_launched == max_nlaunch: 
                break

            while True:
                try:
                    task = work.fetch_task_to_run()

                    if task is None:
                        logger.debug("fetch_task_to_run returned None.")
                        break

                    logger.debug("Starting task %s" % task)
                    if task in launched:
                        err_msg = "task %s in already in launched %s" % (str(task), str(launched))
                        raise RuntimeError(err_msg)

                    task.start()
                    launched.append(task)

                    num_launched += 1

                    if num_launched == max_nlaunch:
                        break

                except StopIteration:
                    logger.debug("Out of tasks.")
                    break

        # Update the database.
        #if num_launched:
        self.flow.check_status()
        self.flow.pickle_dump()

        return num_launched 


class PyFlowsScheduler(object):
    """
    This object schedules the submission of the tasks in an `AbinitFlow`.
    """
    # There are two types of errors that migh happen during the execution of the jobs: Python exceptions and Abinit Errors.
    # Python exceptions are easy to detect and are usually due to a bug in abinitio or random errors such as IOError.
    # The set of Abinit Error is much much broader. It includes wrong input data, segmentation
    # faults, problem with the resource manager, etc. Abinitio tries to handle the most common cases
    # but there's still a lot of room for improvement.
    # The scheduler will shutdown automatically if the number of python exceptions is > MAX_NUM_PYEXC
    # or if the number of Abinit Errors (i.e. the number of tasks whose status is S_ERROR) is > MAX_NUM_ERRORS
    MAX_NUM_PYEXCS = 2
    MAX_NUM_ABIERRS = 0

    YAML_FILE = "scheduler.yml"

    def __init__(self, weeks=0, days=0, hours=0, minutes=0, seconds=0, start_date=None):
        """
        Args:
            weeks:
                number of weeks to wait
            days:
                number of days to wait
            hours:
                number of hours to wait
            minutes:
                number of minutes to wait
            seconds:
                number of seconds to wait
            start_date:
                when to first execute the job and start the counter (default is after the given interval)
        """
        # Options passed to the scheduler.
        self.sched_options = AttrDict(
             weeks=weeks,
             days=days,
             hours=hours,
             minutes=minutes,
             seconds=seconds,
             start_date=start_date,
        )

        from apscheduler.scheduler import Scheduler
        self.sched = Scheduler(standalone=True)

        self._pidfiles2flows = {}

        self.exceptions = []

    @classmethod
    def from_file(cls, filepath):
        """Read the configuration parameters from a Yaml file."""
        with open(filepath, "r") as fh:
            d = yaml.load(fh)
            return cls(**d)

    @classmethod
    def from_user_config(cls):
        """
        Initialize the `PyFlowsScheduler` from the YAML file 'scheduler.yml'.
        Search first in the working directory and then in the configuration
        directory of abipy.

        Raises:
            RuntimeError if file is not found.
        """
        # Try in the current directory.
        path = os.path.join(os.getcwd(), cls.YAML_FILE)

        if os.path.exists(path):
            return cls.from_file(path)

        # Try in the configuration directory.
        home = os.getenv("HOME")
        dirpath = os.path.join(home, ".abinit", "abipy")
        path = os.path.join(dirpath, cls.YAML_FILE)

        if os.path.exists(path):
            return cls.from_file(path)
    
        raise RuntimeError("Cannot locate %s neither in current directory nor in %s" % (cls.YAML_FILE, dirpath))

    def __str__(self):
        """String representation."""
        lines = [self.__class__.__name__ + ", Pid: %d" % self.pid]
        app = lines.append

        app("Scheduler options: %s" % str(self.sched_options)) 

        for flow in self.flows:
            app(80 * "=")
            app(str(flow))

        return "\n".join(lines)

    @property
    def pid(self):
        """Pid of the process"""
        try:
            return self._pid

        except AttributeError:
            self._pid = os.getpid()
            return self._pid

    @property
    def pid_files(self):
        """List of files with the pid. The files are located in the workdir of the flows"""
        return self._pidfiles2flows.keys()
                                                                                            
    @property
    def flows(self):
        """List of `AbinitFlows`."""
        return self._pidfiles2flows.values()

    def add_flow(self, flow):
        """Add an `AbinitFlow` to the scheduler."""
        pid_file = os.path.join(flow.workdir, "_PyFlowsScheduler.pid")

        if pid_file in self._pidfiles2flows:
            raise ValueError("Cannot add the same flow twice!")

        if os.path.isfile(pid_file):
            err_msg = (
                "pid_file %s already exists\n"
                "There are two possibilities:\n\n"
                "       1) There's an another instance of PyFlowsScheduler running.\n"
                "       2) The previous scheduler didn't exit in a clean way.\n\n"
                "To solve case 1:\n"
                "       Kill the previous scheduler (use `kill pid` where pid is the number reported in the file)\n"
                "       Then you run can start the new scheduler.\n\n"
                "To solve case 2:\n"
                "   Remove the pid_file and rerun the scheduler.\n\n"
                "Exiting\n" % pid_file
                )

            raise RuntimeError(err_msg)

        else:
            with open(pid_file, "w") as fh:
                fh.write(str(self.pid))

        self._pidfiles2flows[pid_file] = flow

    def start(self):
        """
        Starts the scheduler in a new thread.
        In standalone mode, this method will block until there are no more scheduled jobs.
        """
        self.sched.add_interval_job(self._runem_all, **self.sched_options)
        self.sched.start()

    def _runem_all(self):
        """The function that will be executed by the scheduler."""
        nlaunch, max_nlaunch, exceptions = 0, -1, []

        for flow in self.flows:
            flow.check_status()

            try:
                nlaunch += PyLauncher(flow).rapidfire(max_nlaunch=max_nlaunch)

            except Exception as exc:
                exceptions.append(straceback())

        print("Number of launches: ", nlaunch)
        flow.show_status()
        
        if exceptions:
            logger.critical("Scheduler exceptions: %s" % str(exceptions))

        self.exceptions.extend(exceptions)

        # Too many exceptions. Shutdown the scheduler.
        if len(self.exceptions) > self.MAX_NUM_PYEXCS:
            msg = "Number of exceptions %s > %s.\nWill shutdown the scheduler and exit" % (
                len(self.exceptions), self.MAX_NUM_PYEXCS)
            msg = 5*"=" + msg + 5*"="
                                                                                      
            print(len(msg)* "=")
            print(msg)
            print(len(msg)* "=")

            self.shutdown()

        # Count the number of tasks with status == S_ERROR:
        num_errors = sum([flow.num_tasks_with_error for flow in self.flows])
        if num_errors > self.MAX_NUM_ABIERRS:
            msg = "Number of task with ERROR status %s > %s. Will shutdown the scheduler and exit" % (
                num_errors, self.MAX_NUM_ABIERRS)
            msg = 5*"=" + msg + 5*"="
                                                                                      
            print(len(msg)* "=")
            print(msg)
            print(len(msg)* "=")
                                                                                            
            self.shutdown()

        # Mission accomplished. Shutdown the scheduler.
        if all(flow.all_ok for flow in self.flows):
            msg = "All tasks have reached S_OK. Will shutdown the scheduler and exit"
            msg = 5*"=" + msg + 5*"="

            print(len(msg)* "=")
            print(msg)
            print(len(msg)* "=")

            self.shutdown()

        return len(exceptions) 

    def shutdown(self):
        """Shutdown the scheduler."""
        for pid_file in self.pid_files:
            try:
                os.unlink(pid_file)
            except:
                pass
                                                                    
        # Save the final status of the flow.
        for flow in self.flows:
            flow.pickle_dump()
                                                                    
        # Shutdown the scheduler thus allowing the process to exit.
        self.sched.shutdown(wait=False)
