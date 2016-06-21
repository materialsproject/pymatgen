# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""Tools for the submission of Tasks."""
from __future__ import unicode_literals, division, print_function

import os
import time
import yaml
import pickle

from collections import deque
from datetime import timedelta
from six.moves import cStringIO
from monty.io import get_open_fds
from monty.string import boxed, is_string
from monty.os.path import which
from monty.collections import AttrDict, dict2namedtuple
from monty.termcolor import cprint
from .utils import as_bool, File, Directory
from . import qutils as qu
from pymatgen.util.io_utils import ask_yesno

try:
    import apscheduler
    has_apscheduler = True
    has_sched_v3 = apscheduler.version >= "3.0.0"
except ImportError:
    has_apscheduler = False

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "ScriptEditor",
    "PyLauncher",
    "PyFlowScheduler",
]


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


class ScriptEditor(object):
    """Simple editor that simplifies the writing of shell scripts"""
    _shell = '/bin/bash'

    def __init__(self):
        self._lines = []

    @property
    def shell(self):
        return self._shell

    def _add(self, text, pre=""):
        if is_string(text):
            self._lines.append(pre + text)
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
        for k, v in d.items():
            self.declare_var(k, v)

    def export_envar(self, key, val):
        """Export an environment variable."""
        line = "export " + key + "=" + str(val)
        self._add(line)

    def export_envars(self, env):
        """Export the environment variables contained in the dict env."""
        for k, v in env.items():
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
        self._add('module load ' + module + " 2>> mods.err")

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


class PyLauncherError(Exception):
    """Error class for PyLauncher."""


class PyLauncher(object):
    """This object handle the submission of the tasks contained in a :class:`Flow`"""
    Error = PyLauncherError

    def __init__(self, flow, **kwargs):
        """
        Initialize the object

        Args:
            flow: :class:`Flow` object
            max_njobs_inqueue: The launcher will stop submitting jobs when the
                number of jobs in the queue is >= Max number of jobs
        """
        self.flow = flow
        self.max_njobs_inqueue = kwargs.get("max_njobs_inqueue", 200)

        #self.flow.check_pid_file()

    def single_shot(self):
        """
        Run the first :class:`Task` than is ready for execution.

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
                    # No task found, this usually happens when we have dependencies.
                    # Beware of possible deadlocks here!
                    logger.debug("No task to run! Possible deadlock")

            except StopIteration:
                logger.info("All tasks completed.")

        # Submit the tasks and update the database.
        if tasks:
            tasks[0].start()
            num_launched += 1

            self.flow.pickle_dump()

        return num_launched

    def rapidfire(self, max_nlaunch=-1, max_loops=1, sleep_time=5):
        """
        Keeps submitting `Tasks` until we are out of jobs or no job is ready to run.

        Args:
            max_nlaunch: Maximum number of launches. default: no limit.
            max_loops: Maximum number of loops
            sleep_time: seconds to sleep between rapidfire loop iterations

        Returns:
            The number of tasks launched.
        """
        num_launched, do_exit, launched = 0, False, []

        for count in range(max_loops):
            if do_exit:
                break
            if count > 0:
                time.sleep(sleep_time)

            tasks = self.fetch_tasks_to_run()

            # I don't know why but we receive duplicated tasks.
            if any(task in launched for task in tasks):
                logger.critical("numtasks %d already in launched list:\n%s" % (len(tasks), launched))

            # Preventive test.
            tasks = [t for t in tasks if t not in launched]

            if not tasks:
                continue

            for task in tasks:
                fired = task.start()
                if fired:
                    launched.append(task)
                    num_launched += 1

                if num_launched >= max_nlaunch > 0:
                    logger.info('num_launched >= max_nlaunch, going back to sleep')
                    do_exit = True
                    break

        # Update the database.
        self.flow.pickle_dump()

        return num_launched

    def fetch_tasks_to_run(self):
        """
        Return the list of tasks that can be submitted.
        Empty list if no task has been found.
        """
        tasks_to_run = []

        for work in self.flow:
            tasks_to_run.extend(work.fetch_alltasks_to_run())

        return tasks_to_run


class PyFlowSchedulerError(Exception):
    """Exceptions raised by `PyFlowScheduler`."""


class PyFlowScheduler(object):
    """
    This object schedules the submission of the tasks in a :class:`Flow`.
    There are two types of errors that might occur during the execution of the jobs:

        #. Python exceptions
        #. Errors in the ab-initio code

    Python exceptions are easy to detect and are usually due to a bug in the python code or random errors such as IOError.
    The set of errors in the ab-initio is much much broader. It includes wrong input data, segmentation
    faults, problems with the resource manager, etc. The flow tries to handle the most common cases
    but there's still a lot of room for improvement.
    Note, in particular, that `PyFlowScheduler` will shutdown automatically in the following cases:

        #. The number of python exceptions is > max_num_pyexcs

        #. The number of task errors (i.e. the number of tasks whose status is S_ERROR) is > max_num_abierrs

        #. The number of jobs launched becomes greater than (`safety_ratio` * total_number_of_tasks).

        #. The scheduler will send an email to the user (specified by `mailto`) every `remindme_s` seconds.
           If the mail cannot be sent, the scheduler will shutdown automatically.
           This check prevents the scheduler from being trapped in an infinite loop.
    """
    # Configuration file.
    YAML_FILE = "scheduler.yml"
    USER_CONFIG_DIR = os.path.join(os.getenv("HOME"), ".abinit", "abipy")

    Error = PyFlowSchedulerError

    @classmethod
    def autodoc(cls):
        i = cls.__init__.__doc__.index("Args:")
        return cls.__init__.__doc__[i+5:]

    def __init__(self, **kwargs):
        """
        Args:
            weeks: number of weeks to wait (DEFAULT: 0).
            days: number of days to wait (DEFAULT: 0).
            hours: number of hours to wait (DEFAULT: 0).
            minutes: number of minutes to wait (DEFAULT: 0).
            seconds: number of seconds to wait (DEFAULT: 0).
            mailto: The scheduler will send an email to `mailto` every `remindme_s` seconds.
                (DEFAULT: None i.e. not used).
            verbose: (int) verbosity level. (DEFAULT: 0)
            use_dynamic_manager: "yes" if the :class:`TaskManager` must be re-initialized from
                file before launching the jobs. (DEFAULT: "no")
            max_njobs_inqueue: Limit on the number of jobs that can be present in the queue. (DEFAULT: 200)
            remindme_s: The scheduler will send an email to the user specified by `mailto` every `remindme_s` seconds.
                (int, DEFAULT: 1 day).
            max_num_pyexcs: The scheduler will exit if the number of python exceptions is > max_num_pyexcs
                (int, DEFAULT: 0)
            max_num_abierrs: The scheduler will exit if the number of errored tasks is > max_num_abierrs
                (int, DEFAULT: 0)
            safety_ratio: The scheduler will exits if the number of jobs launched becomes greater than
               `safety_ratio` * total_number_of_tasks_in_flow. (int, DEFAULT: 5)
            max_nlaunches: Maximum number of tasks launched in a single iteration of the scheduler.
                (DEFAULT: -1 i.e. no limit)
            debug: Debug level. Use 0 for production (int, DEFAULT: 0)
            fix_qcritical: "yes" if the launcher should try to fix QCritical Errors (DEFAULT: "yes")
            rmflow: If "yes", the scheduler will remove the flow directory if the calculation
                completed successfully. (DEFAULT: "no")
            killjobs_if_errors: "yes" if the scheduler should try to kill all the runnnig jobs
                before exiting due to an error. (DEFAULT: "yes")
        """
        # Options passed to the scheduler.
        self.sched_options = AttrDict(
            weeks=kwargs.pop("weeks", 0),
            days=kwargs.pop("days", 0),
            hours=kwargs.pop("hours", 0),
            minutes=kwargs.pop("minutes", 0),
            seconds=kwargs.pop("seconds", 0),
            #start_date=kwargs.pop("start_date", None),
        )
        if all(not v for v in self.sched_options.values()):
            raise self.Error("Wrong set of options passed to the scheduler.")

        self.mailto = kwargs.pop("mailto", None)
        self.verbose = int(kwargs.pop("verbose", 0))
        self.use_dynamic_manager = as_bool(kwargs.pop("use_dynamic_manager", False))
        self.max_njobs_inqueue = kwargs.pop("max_njobs_inqueue", 200)
        self.max_ncores_used = kwargs.pop("max_ncores_used", None)
        self.contact_resource_manager = as_bool(kwargs.pop("contact_resource_manager", False))

        self.remindme_s = float(kwargs.pop("remindme_s", 1 * 24 * 3600))
        self.max_num_pyexcs = int(kwargs.pop("max_num_pyexcs", 0))
        self.max_num_abierrs = int(kwargs.pop("max_num_abierrs", 0))
        self.safety_ratio = int(kwargs.pop("safety_ratio", 5))
        #self.max_etime_s = kwargs.pop("max_etime_s", )
        self.max_nlaunches = kwargs.pop("max_nlaunches", -1)
        self.debug = kwargs.pop("debug", 0)
        self.fix_qcritical = as_bool(kwargs.pop("fix_qcritical", True))
        self.rmflow = as_bool(kwargs.pop("rmflow", False))
        self.killjobs_if_errors = as_bool(kwargs.pop("killjobs_if_errors", True))

        self.customer_service_dir = kwargs.pop("customer_service_dir", None)
        if self.customer_service_dir is not None:
            self.customer_service_dir = Directory(self.customer_service_dir)
            self._validate_customer_service()

        if kwargs:
            raise self.Error("Unknown arguments %s" % kwargs)

        if not has_apscheduler:
            raise RuntimeError("Install apscheduler with pip")

        if has_sched_v3:
            logger.warning("Using scheduler v>=3.0.0")
            from apscheduler.schedulers.blocking import BlockingScheduler
            self.sched = BlockingScheduler()
        else:
            from apscheduler.scheduler import Scheduler
            self.sched = Scheduler(standalone=True)

        self.nlaunch = 0
        self.num_reminders = 1

        # Used to keep track of the exceptions raised while the scheduler is running
        self.exceptions = deque(maxlen=self.max_num_pyexcs + 10)

        # Used to push additional info during the execution.
        self.history = deque(maxlen=100)

    @classmethod
    def from_file(cls, filepath):
        """Read the configuration parameters from a Yaml file."""
        with open(filepath, "r") as fh:
            return cls(**yaml.load(fh))

    @classmethod
    def from_string(cls, s):
        """Create an istance from string s containing a YAML dictionary."""
        stream = cStringIO(s)
        stream.seek(0)
        return cls(**yaml.load(stream))

    @classmethod
    def from_user_config(cls):
        """
        Initialize the :class:`PyFlowScheduler` from the YAML file 'scheduler.yml'.
        Search first in the working directory and then in the configuration directory of abipy.

        Raises:
            `RuntimeError` if file is not found.
        """
        # Try in the current directory.
        path = os.path.join(os.getcwd(), cls.YAML_FILE)

        if os.path.exists(path):
            return cls.from_file(path)

        # Try in the configuration directory.
        path = os.path.join(cls.USER_CONFIG_DIR, cls.YAML_FILE)

        if os.path.exists(path):
            return cls.from_file(path)

        raise cls.Error("Cannot locate %s neither in current directory nor in %s" % (cls.YAML_FILE, path))

    def __str__(self):
        """String representation."""
        lines = [self.__class__.__name__ + ", Pid: %d" % self.pid]
        app = lines.append
        app("Scheduler options: %s" % str(self.sched_options))

        if self.flow is not None:
            app(80 * "=")
            app(str(self.flow))

        return "\n".join(lines)

    @property
    def pid(self):
        """The pid of the process associated to the scheduler."""
        try:
            return self._pid
        except AttributeError:
            self._pid = os.getpid()
            return self._pid

    @property
    def pid_file(self):
        """
        Absolute path of the file with the pid.
        The file is located in the workdir of the flow
        """
        return self._pid_file

    @property
    def flow(self):
        """`Flow`."""
        try:
            return self._flow
        except AttributeError:
            return None

    @property
    def num_excs(self):
        """Number of exceptions raised so far."""
        return len(self.exceptions)

    def get_delta_etime(self):
        """Returns a `timedelta` object representing with the elapsed time."""
        return timedelta(seconds=(time.time() - self.start_time))

    def add_flow(self, flow):
        """
        Add an :class:`Flow` flow to the scheduler.
        """
        if hasattr(self, "_flow"):
            raise self.Error("Only one flow can be added to the scheduler.")

        # Check if we are already using a scheduler to run this flow
        flow.check_pid_file()
        flow.set_spectator_mode(False)

        # Build dirs and files (if not yet done)
        flow.build()

        with open(flow.pid_file, "w") as fh:
            fh.write(str(self.pid))

        self._pid_file = flow.pid_file
        self._flow = flow

    def _validate_customer_service(self):
        """
        Validate input parameters if customer service is on then
        create directory for tarball files with correct premissions for user and group.
        """
        direc = self.customer_service_dir
        if not direc.exists:
            mode = 0o750
            print("Creating customer_service_dir %s with mode %s" % (direc, mode))
            direc.makedirs()
            os.chmod(direc.path, mode)

        if self.mailto is None:
            raise RuntimeError("customer_service_dir requires mailto option in scheduler.yml")

    def _do_customer_service(self):
        """
        This method is called before the shutdown of the scheduler.
        If customer_service is on and the flow didn't completed successfully,
        a lightweight tarball file with inputs and the most important output files
        is created in customer_servide_dir.
        """
        if self.customer_service_dir is None: return
        doit = self.exceptions or not self.flow.all_ok
        doit = True
        if not doit: return

        prefix = os.path.basename(self.flow.workdir) + "_"

        import tempfile, datetime
        suffix = str(datetime.datetime.now()).replace(" ", "-")
        # Remove milliseconds
        i = suffix.index(".")
        if i != -1: suffix = suffix[:i]
        suffix += ".tar.gz"

        #back = os.getcwd()
        #os.chdir(self.customer_service_dir.path)

        _, tmpname = tempfile.mkstemp(suffix="_" + suffix, prefix=prefix,
                                      dir=self.customer_service_dir.path, text=False)

        print("Dear customer,\n We are about to generate a tarball in\n  %s" % tmpname)
        self.flow.make_light_tarfile(name=tmpname)
        #os.chdir(back)

    def start(self):
        """
        Starts the scheduler in a new thread. Returns 0 if success.
        In standalone mode, this method will block until there are no more scheduled jobs.
        """
        self.history.append("Started on %s" % time.asctime())
        self.start_time = time.time()

        if not has_apscheduler:
            raise RuntimeError("Install apscheduler with pip")

        if has_sched_v3:
            self.sched.add_job(self.callback, "interval", **self.sched_options)
        else:
            self.sched.add_interval_job(self.callback, **self.sched_options)

        errors = self.flow.look_before_you_leap()
        if errors:
            self.exceptions.append(errors)
            return 1

        # Try to run the job immediately. If something goes wrong return without initializing the scheduler.
        self._runem_all()

        if self.exceptions:
            self.cleanup()
            self.send_email(msg="Error while trying to run the flow for the first time!\n %s" % self.exceptions)
            return 1

        try:
            self.sched.start()
            return 0

        except KeyboardInterrupt:
            self.shutdown(msg="KeyboardInterrupt from user")
            if ask_yesno("Do you want to cancel all the jobs in the queue? [Y/n]"):
                print("Number of jobs cancelled:", self.flow.cancel())

            self.flow.pickle_dump()
            return -1

    def _runem_all(self):
        """
        This function checks the status of all tasks,
        tries to fix tasks that went unconverged, abicritical, or queuecritical
        and tries to run all the tasks that can be submitted.+
        """
        excs = []
        flow = self.flow

        # Allow to change the manager at run-time
        if self.use_dynamic_manager:
            from pymatgen.io.abinit.tasks import TaskManager
            new_manager = TaskManager.from_user_config()
            for work in flow:
                work.set_manager(new_manager)

        nqjobs = 0
        if self.contact_resource_manager:
            # This call is expensive and therefore it's optional
            nqjobs = flow.get_njobs_in_queue()
            if nqjobs is None:
                nqjobs = 0
                if flow.manager.has_queue: logger.warning('Cannot get njobs_inqueue')

            if nqjobs >= self.max_njobs_inqueue:
                print("Too many jobs in the queue: %s, returning" % nqjobs)
                return

        if self.max_nlaunches == -1:
            max_nlaunch = self.max_njobs_inqueue - nqjobs
        else:
            max_nlaunch = min(self.max_njobs_inqueue - nqjobs, self.max_nlaunches)

        # check status.
        flow.check_status(show=False)

        # This check is not perfect, we should make a list of tasks to sumbit
        # and select only the subset so that we don't exceeed mac_ncores_used
        # Many sections of this code should be rewritten.
        #if self.max_ncores_used is not None and flow.ncores_used > self.max_ncores_used:
        if self.max_ncores_used is not None and flow.ncores_allocated > self.max_ncores_used:
            print("Cannot exceed max_ncores_use:d %s" % self.max_ncores_used)
            return

        # Try to restart the unconverged tasks
        # TODO: do not fire here but prepare for fireing in rapidfire
        for task in self.flow.unconverged_tasks:
            try:
                logger.info("Flow will try restart task %s" % task)
                fired = task.restart()
                if fired:
                    self.nlaunch += 1
                    max_nlaunch -= 1
                    if max_nlaunch == 0:
                        logger.info("Restart: too many jobs in the queue, returning")
                        flow.pickle_dump()
                        return

            except task.RestartError:
                excs.append(straceback())

        # Temporarily disable by MG because I don't know if fix_critical works after the
        # introduction of the new qadapters
        # reenabled by MsS disable things that do not work at low level
        # fix only prepares for restarting, and sets to ready
        if self.fix_qcritical:
            nfixed = flow.fix_queue_critical()
            if nfixed: print("Fixed %d QCritical error(s)" % nfixed)

        nfixed = flow.fix_abicritical()
        if nfixed: print("Fixed %d AbiCritical error(s)" % nfixed)

        # update database
        flow.pickle_dump()

        # Submit the tasks that are ready.
        try:
            nlaunch = PyLauncher(flow).rapidfire(max_nlaunch=max_nlaunch, sleep_time=10)
            self.nlaunch += nlaunch

            if nlaunch:
                print("[%s] Number of launches: %d" % (time.asctime(), nlaunch))

        except Exception:
            excs.append(straceback())

        # check status.
        flow.show_status()

        if excs:
            logger.critical("*** Scheduler exceptions:\n *** %s" % "\n".join(excs))
            self.exceptions.extend(excs)

    def callback(self):
        """The function that will be executed by the scheduler."""
        try:
            return self._callback()
        except:
            # All exceptions raised here will trigger the shutdown!
            s = straceback()
            self.exceptions.append(s)

            # This is useful when debugging
            #try:
            #    print("Exception in callback, will cancel all tasks")
            #    for task in self.flow.iflat_tasks():
            #        task.cancel()
            #except Exception:
            #    pass

            self.shutdown(msg="Exception raised in callback!\n" + s)

    def _callback(self):
        """The actual callback."""
        if self.debug:
            # Show the number of open file descriptors
            print(">>>>> _callback: Number of open file descriptors: %s" % get_open_fds())

        self._runem_all()

        # Mission accomplished. Shutdown the scheduler.
        all_ok = self.flow.all_ok
        if all_ok:
            return self.shutdown(msg="All tasks have reached S_OK. Will shutdown the scheduler and exit")

        # Handle failures.
        err_lines = []

        # Shall we send a reminder to the user?
        delta_etime = self.get_delta_etime()

        if delta_etime.total_seconds() > self.num_reminders * self.remindme_s:
            self.num_reminders += 1
            msg = ("Just to remind you that the scheduler with pid %s, flow %s\n has been running for %s " %
                  (self.pid, self.flow, delta_etime))
            retcode = self.send_email(msg, tag="[REMINDER]")

            if retcode:
                # Cannot send mail, shutdown now!
                msg += ("\nThe scheduler tried to send an e-mail to remind the user\n" +
                        " but send_email returned %d. Aborting now" % retcode)
                err_lines.append(msg)

        #if delta_etime.total_seconds() > self.max_etime_s:
        #    err_lines.append("\nExceeded max_etime_s %s. Will shutdown the scheduler and exit" % self.max_etime_s)

        # Too many exceptions. Shutdown the scheduler.
        if self.num_excs > self.max_num_pyexcs:
            msg = "Number of exceptions %s > %s. Will shutdown the scheduler and exit" % (
                self.num_excs, self.max_num_pyexcs)
            err_lines.append(boxed(msg))

        # Paranoid check: disable the scheduler if we have submitted
        # too many jobs (it might be due to some bug or other external reasons
        # such as race conditions between difference callbacks!)
        if self.nlaunch > self.safety_ratio * self.flow.num_tasks:
            msg = "Too many jobs launched %d. Total number of tasks = %s, Will shutdown the scheduler and exit" % (
                self.nlaunch, self.flow.num_tasks)
            err_lines.append(boxed(msg))

        # Count the number of tasks with status == S_ERROR.
        if self.flow.num_errored_tasks > self.max_num_abierrs:
            msg = "Number of tasks with ERROR status %s > %s. Will shutdown the scheduler and exit" % (
                self.flow.num_errored_tasks, self.max_num_abierrs)
            err_lines.append(boxed(msg))

        # Test on the presence of deadlocks.
        g = self.flow.find_deadlocks()
        if g.deadlocked:
            # Check the flow again so that status are updated.
            self.flow.check_status()

            g = self.flow.find_deadlocks()
            print("deadlocked:\n", g.deadlocked, "\nrunnables:\n", g.runnables, "\nrunning\n", g.running)
            if g.deadlocked and not g.runnables and not g.running:
                err_lines.append("No runnable job with deadlocked tasks:\n%s." % str(g.deadlocked))

        if not g.runnables and not g.running:
            # Check the flow again so that status are updated.
            self.flow.check_status()
            g = self.flow.find_deadlocks()
            if not g.runnables and not g.running:
                err_lines.append("No task is running and cannot find other tasks to submit.")

        # Something wrong. Quit
        if err_lines:
            # Cancel all jobs.
            if self.killjobs_if_errors:
                cprint("killjobs_if_errors set to 'yes' in scheduler file. Will kill jobs before exiting.", "yellow")
                try:
                    num_cancelled = 0
                    for task in self.flow.iflat_tasks():
                        num_cancelled += task.cancel()
                    cprint("Killed %d tasks" % num_cancelled, "yellow")
                except Exception as exc:
                    cprint("Exception while trying to kill jobs:\n%s" % str(exc), "red")

            self.shutdown("\n".join(err_lines))

        return len(self.exceptions)

    def cleanup(self):
        """Cleanup routine: remove the pid file and save the pickle database"""
        try:
            os.remove(self.pid_file)
        except OSError as exc:
            logger.critical("Could not remove pid_file: %s", exc)

        # Save the final status of the flow.
        self.flow.pickle_dump()

    def shutdown(self, msg):
        """Shutdown the scheduler."""
        try:
            self.cleanup()

            self.history.append("Completed on: %s" % time.asctime())
            self.history.append("Elapsed time: %s" % self.get_delta_etime())

            if self.debug:
                print(">>>>> shutdown: Number of open file descriptors: %s" % get_open_fds())

            retcode = self.send_email(msg)
            if self.debug:
                print("send_mail retcode", retcode)

            # Write file with the list of exceptions:
            if self.exceptions:
                dump_file = os.path.join(self.flow.workdir, "_exceptions")
                with open(dump_file, "w") as fh:
                    fh.writelines(self.exceptions)
                    fh.write("Shutdown message:\n%s" % msg)

            lines = []
            app = lines.append
            app("Submitted on: %s" % time.ctime(self.start_time))
            app("Completed on: %s" % time.asctime())
            app("Elapsed time: %s" % str(self.get_delta_etime()))

            if self.flow.all_ok:
                app("Flow completed successfully")
            else:
                app("Flow %s didn't complete successfully" % repr(self.flow.workdir))
                app("use `abirun.py FLOWDIR debug` to analyze the problem.")
                app("Shutdown message:\n%s" % msg)

            print("")
            print("\n".join(lines))
            print("")

            self._do_customer_service()

            if self.flow.all_ok:
                print("Calling flow.finalize()...")
                self.flow.finalize()
                if self.rmflow:
                    app("Flow directory will be removed...")
                    try:
                        self.flow.rmtree()
                    except Exception:
                        logger.warning("Ignoring exception while trying to remove flow dir.")

        finally:
            # Shutdown the scheduler thus allowing the process to exit.
            logger.debug('This should be the shutdown of the scheduler')

            # Unschedule all the jobs before calling shutdown
            #self.sched.print_jobs()
            if not has_sched_v3:
                for job in self.sched.get_jobs():
                    self.sched.unschedule_job(job)
            #self.sched.print_jobs()

            self.sched.shutdown()
            # Uncomment the line below if shutdown does not work!
            #os.system("kill -9 %d" % os.getpid())

    def send_email(self, msg, tag=None):
        """
        Send an e-mail before completing the shutdown.
        Returns 0 if success.
        """
        try:
            return self._send_email(msg, tag)
        except:
            self.exceptions.append(straceback())
            return -2

    def _send_email(self, msg, tag):
        if self.mailto is None:
            return -1

        header = msg.splitlines()
        app = header.append

        app("Submitted on: %s" % time.ctime(self.start_time))
        app("Completed on: %s" % time.asctime())
        app("Elapsed time: %s" % str(self.get_delta_etime()))
        app("Number of errored tasks: %d" % self.flow.num_errored_tasks)
        app("Number of unconverged tasks: %d" % self.flow.num_unconverged_tasks)

        strio = cStringIO()
        strio.writelines("\n".join(header) + 4 * "\n")

        # Add the status of the flow.
        self.flow.show_status(stream=strio)

        if self.exceptions:
            # Report the list of exceptions.
            strio.writelines(self.exceptions)

        if tag is None:
            tag = " [ALL OK]" if self.flow.all_ok else " [WARNING]"

        return sendmail(subject=self.flow.name + tag, text=strio.getvalue(), mailto=self.mailto)


def sendmail(subject, text, mailto, sender=None):
    """
    Sends an e-mail with unix sendmail.

    Args:
        subject: String with the subject of the mail.
        text: String with the body of the mail.
        mailto: String or list of string with the recipients.
        sender: string with the sender address.
            If sender is None, username@hostname is used.

    Returns:
        Exit status
    """
    def user_at_host():
        from socket import gethostname
        return os.getlogin() + "@" + gethostname()

    # Body of the message.
    try:
        sender = user_at_host() if sender is None else sender
    except OSError:
        sender = 'abipyscheduler@youknowwhere'

    if is_string(mailto): mailto = [mailto]

    from email.mime.text import MIMEText
    mail = MIMEText(text)
    mail["Subject"] = subject
    mail["From"] = sender
    mail["To"] = ", ".join(mailto)

    msg = mail.as_string()

    # sendmail works much better than the python interface.
    # Note that sendmail is available only on Unix-like OS.
    from subprocess import Popen, PIPE

    sendmail = which("sendmail")
    if sendmail is None: return -1
    p = Popen([sendmail, "-t"], stdin=PIPE, stderr=PIPE)

    outdata, errdata = p.communicate(msg)
    return len(errdata)


def __test_sendmail():
    retcode = sendmail("sendmail_test", text="hello\nworld", mailto="nobody@nowhere.com")
    print("Retcode", retcode)
    assert retcode == 0


class BatchLauncherError(Exception):
    """Exceptions raised by :class:`BatchLauncher`."""


class BatchLauncher(object):
    """
    This object automates the execution of multiple flow. It generates a job script
    that uses abirun.py to run each flow stored in self with a scheduler.
    The execution of the flows is done in sequential but each scheduler will start
    to submit the tasks of the flow in autoparal mode.

    The `BatchLauncher` is pickleable, hence one can reload it, check if all flows are completed
    and rerun only those that are not completed due to the timelimit.
    """
    PICKLE_FNAME = "__BatchLauncher__.pickle"

    Error = BatchLauncherError

    @classmethod
    def from_dir(cls, top, workdir=None, name=None, manager=None, max_depth=2):
        """
        Find all flows located withing the directory `top` and build the `BatchLauncher`.

        Args:
            top: Top level directory or list of directories.
            workdir: Batch workdir.
            name:
            manager: :class:`TaskManager` object. If None, the manager is read from `manager.yml`
                In this case the YAML file must provide the entry `batch_manager` that defined
                the queue adapter used to submit the batch script.
            max_depth: Search in directory only if it is N or fewer levels below top
        """
        from .flows import Flow
        def find_pickles(dirtop):
            # Walk through each directory inside path and find the pickle database.
            paths = []
            for dirpath, dirnames, filenames in os.walk(dirtop):
                fnames = [f for f in filenames if f == Flow.PICKLE_FNAME]
                paths.extend([os.path.join(dirpath, f) for f in fnames])
            return paths

        if is_string(top):
            pickle_paths = find_pickles(top)
        else:
            # List of directories.
            pickle_paths = []
            for p in top:
                pickle_paths.extend(find_pickles(p))

        #workdir = os.path.join(top, "batch") if workdir is None else workdir
        workdir = "batch" if workdir is None else workdir
        new = cls(workdir, name=name, manager=manager)

        for path in pickle_paths:
            new.add_flow(path)

        return new

    @classmethod
    def pickle_load(cls, filepath):
        """
        Loads the object from a pickle file.

        Args:
            filepath: Filename or directory name. It filepath is a directory, we
                scan the directory tree starting from filepath and we
                read the first pickle database. Raise RuntimeError if multiple
                databases are found.
        """
        if os.path.isdir(filepath):
            # Walk through each directory inside path and find the pickle database.
            for dirpath, dirnames, filenames in os.walk(filepath):
                fnames = [f for f in filenames if f == cls.PICKLE_FNAME]
                if fnames:
                    if len(fnames) == 1:
                        filepath = os.path.join(dirpath, fnames[0])
                        break  # Exit os.walk
                    else:
                        err_msg = "Found multiple databases:\n %s" % str(fnames)
                        raise RuntimeError(err_msg)
            else:
                err_msg = "Cannot find %s inside directory %s" % (cls.PICKLE_FNAME, filepath)
                raise ValueError(err_msg)

        with open(filepath, "rb") as fh:
            new = pickle.load(fh)

        # new.flows is a list of strings with the workdir of the flows (see __getstate__).
        # Here we read the Flow from the pickle file so that we have
        # and up-to-date version and we set the flow in visitor_mode
        from .flows import Flow
        flow_workdirs, new.flows = new.flows, []
        for flow in map(Flow.pickle_load, flow_workdirs):
            new.add_flow(flow)

        return new

    def pickle_dump(self):
        """Save the status of the object in pickle format."""
        with open(os.path.join(self.workdir, self.PICKLE_FNAME), mode="wb") as fh:
            pickle.dump(self, fh)

    def __getstate__(self):
        """
        Return state is pickled as the contents for the instance.

        Here we replace the flow objects with their workdir because we are observing
        the flows and we want to have the updated version when we reload the `BatchLauncher` from pickle.
        """
        d = {k: v for k, v in self.__dict__.items() if k not in ["flows"]}
        d["flows"] = [flow.workdir for flow in self.flows]
        return d

    def __init__(self, workdir, name=None, flows=None, manager=None, timelimit=None):
        """
        Args:
            workdir: Working directory
            name: Name assigned to the `BatchLauncher`.
            flows:  List of `Flow` objects.
            manager: :class:`TaskManager` object responsible for the submission of the jobs.
                     If manager is None, the object is initialized from the yaml file
                     located either in the working directory or in the user configuration dir.
            timelimit: Time limit (int with seconds or string with time given with
                the slurm convention: "days-hours:minutes:seconds".
                If timelimit is None, the default value specified in the `batch_adapter` is taken.
        """
        self.workdir = os.path.abspath(workdir)

        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        else:
            pass
            #raise RuntimeError("Directory %s already exists. Use BatchLauncher.pickle_load()" % self.workdir)

        self.name = os.path.basename(self.workdir) if name is None else name
        self.script_file = File(os.path.join(self.workdir, "run.sh"))
        self.qerr_file = File(os.path.join(self.workdir, "queue.qerr"))
        self.qout_file = File(os.path.join(self.workdir, "queue.qout"))
        self.log_file = File(os.path.join(self.workdir, "run.log"))
        self.batch_pidfile = File(os.path.join(self.workdir, "batch.pid"))

        from .tasks import TaskManager
        manager = TaskManager.as_manager(manager)

        # Extract the qadapater to be used for the batch script.
        try:
            self.qadapter = qad = manager.batch_adapter
        except AttributeError:
            raise RuntimeError("Your manager.yml file does not define an entry for the batch_adapter")

        if qad is None:
            raise RuntimeError("Your manager.yml file does not define an entry for the batch_adapter")

        # Set mpi_procs to 1 just to be on the safe side
        # Then allow the user to change the timelimit via __init__
        qad.set_mpi_procs(1)
        if timelimit is not None:
            self.set_timelimit(timelimit)
            # FIXME: Remove me!
            self.set_timelimit(36000)

        # Initialize list of flows.
        if flows is None: flows = []
        if not isinstance(flows, (list, tuple)): flows = [flows]
        self.flows = flows

    def set_timelimit(self, timelimit):
        """
        Set the timelimit of the batch launcher.

        Args:
            timelimit: Time limit (int with seconds or string with time given
                with the slurm convention: "days-hours:minutes:seconds".
        """
        self.qad.set_timelimit(qu.timelimit_parser(timelimit))

    def to_string(self, **kwargs):
        lines = []
        lines.extend(str(self.qadapter).splitlines())

        for i, flow in enumerate(self.flows):
            lines.append("Flow [%d] " % i + str(flow))

        return "\n".join(lines)

    def __str__(self):
        return self.to_string()

    def add_flow(self, flow):
        """
        Add a flow. Accept filepath or :class:`Flow` object. Return 1 if flow was added else 0.
        """
        from .flows import Flow
        flow = Flow.as_flow(flow)

        if flow in self.flows:
            raise self.Error("Cannot add same flow twice!")

        if not flow.allocated:
            # Set the workdir of the flow here. Create a dir in self.workdir with name flow.name
            flow_workdir = os.path.join(self.workdir, os.path.basename(flow.name))
            if flow_workdir in (flow.workdir for flow in self.flows):
                raise self.Error("Two flows have the same name and hence the same workdir!")
            flow.allocate(workdir=flow_workdir)

        # Check if we are already using a scheduler to run this flow
        flow.check_pid_file()
        flow.set_spectator_mode(False)

        flow.check_status(show=False)

        #if flow.all_ok:
        #    print("flow.all_ok: Ignoring %s" % flow)
        #    return 0

        self.flows.append(flow)
        #print("Flow %s added to the BatchLauncher" % flow)

        return 1

    def submit(self, **kwargs):
        """
        Submit a job script that will run the schedulers with `abirun.py`.

        Args:
            verbose: Verbosity level
            dry_run: Don't submit the script if dry_run. Default: False

        Returns:
            namedtuple with attributes:
                retcode: Return code as returned by the submission script.
                qjob: :class:`QueueJob` object.
                num_flows_inbatch: Number of flows executed by the batch script

            Return code of the job script submission.
        """
        verbose, dry_run = kwargs.pop("verbose", 0), kwargs.pop("dry_run", False)

        if not self.flows:
            print("Cannot submit an empty list of flows!")
            return 0

        if hasattr(self, "qjob"):
            # This usually happens when we have loaded the object from pickle
            # and we have already submitted to batch script to the queue.
            # At this point we need to understand if the previous batch job
            # is still running before trying to submit it again. There are three cases:
            #
            # 1) The batch script has completed withing timelimit and therefore
            #    the pid_file has been removed by the script. In this case, we
            #    should not try to submit it again.

            # 2) The batch script has been killed due to timelimit (other reasons are possible
            #    but we neglect them). In this case the pid_file exists but there's no job with
            #    this pid runnig and we can resubmit it again.

            # 3) The batch script is still running.
            print("BatchLauncher has qjob %s" % self.qjob)

            if not self.batch_pid_file.exists:
                print("It seems that the batch script reached the end. Wont' try to submit it again")
                return 0

            msg = ("Here I have to understand if qjob is in the queue."
                   " but I need an abstract API that can retrieve info from the queue id")
            raise RuntimeError(msg)

            # TODO: Temptative API
            if self.qjob.in_status("Running|Queued"):
                print("Job is still running. Cannot submit")
            else:
                del self.qjob

        script, num_flows_inbatch = self._get_script_nflows()

        if num_flows_inbatch == 0:
            print("All flows have reached all_ok! Batch script won't be submitted")
            return 0

        if verbose:
            print("*** submission script ***")
            print(script)

        # Write the script.
        self.script_file.write(script)
        self.script_file.chmod(0o740)

        # Builf the flow.
        for flow in self.flows:
            flow.build_and_pickle_dump()

        # Submit the task and save the queue id.
        if dry_run: return -1

        print("Will submit %s flows in batch script" % len(self.flows))
        self.qjob, process = self.qadapter.submit_to_queue(self.script_file.path)

        # Save the queue id in the pid file
        # The file will be removed by the job script if execution is completed.
        self.batch_pidfile.write(str(self.qjob.qid))

        self.pickle_dump()
        process.wait()

        return dict2namedtuple(retcode=process.returncode, qjob=self.qjob,
                               num_flows_inbatch=num_flows_inbatch)

    def _get_script_nflows(self):
        """
        Write the submission script. Return (script, num_flows_in_batch)
        """
        flows_torun = [f for f in self.flows if not f.all_ok]
        if not flows_torun:
            return "", 0

        executable = [
            'export _LOG=%s' % self.log_file.path,
            'date1=$(date +"%s")',
            'echo Running abirun.py in batch mode > ${_LOG}',
            " ",
        ]
        app = executable.append

        # Build list of abirun commands and save the name of the log files.
        self.sched_logs, num_flows = [], len(flows_torun)
        for i, flow in enumerate(flows_torun):

            logfile = os.path.join(self.workdir, "log_" + os.path.basename(flow.workdir))

            app("echo Starting flow %d/%d on: `date` >> ${LOG}" % (i+1, num_flows))
            app("\nabirun.py %s scheduler > %s" % (flow.workdir, logfile))
            app("echo Returning from abirun on `date` with retcode $? >> ${_LOG}")

            assert logfile not in self.sched_logs
            self.sched_logs.append(logfile)

        # Remove the batch pid_file and compute elapsed time.
        executable.extend([
            " ",
            "# Remove batch pid file",
            'rm %s' % self.batch_pidfile.path,
            " ",
            "# Compute elapsed time",
            'date2=$(date +"%s")',
            'diff=$(($date2-$date1))',
            'echo $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed. >> ${_LOG}'
        ])

        return self.qadapter.get_script_str(
            job_name=self.name,
            launch_dir=self.workdir,
            executable=executable,
            qout_path=self.qout_file.path,
            qerr_path=self.qerr_file.path,
        ), num_flows

    def show_summary(self, **kwargs):
        """
        Show a summary with the status of the flows.
        """
        for flow in self.flows:
            flow.show_summary()

    def show_status(self, **kwargs):
        """
        Report the status of the flows.

        Args:
            stream: File-like object, Default: sys.stdout
            verbose: Verbosity level (default 0). > 0 to show only the works that are not finalized.
        """
        for flow in self.flows:
            flow.show_status(**kwargs)
