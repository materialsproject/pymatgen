# coding: utf-8
"""Tools for the submission of Tasks."""
from __future__ import unicode_literals, division, print_function

import os
import time
import collections
import yaml

from six.moves import cStringIO
from datetime import timedelta
from monty.io import get_open_fds
from monty.string import boxed, is_string
from monty.os.path import which
from monty.collections import AttrDict
from .utils import as_bool, File, Directory

try:
    import apscheduler
    has_sched_v3 = apscheduler.version >= "3.0.0"
except ImportError:
    pass

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


def ask_yesno(prompt, default=True):
    import six
    # Fix python 2.x.
    if six.PY2:
        my_input = raw_input
    else:
        my_input = input

    try:
        answer = my_input(prompt)
        return answer.lower().strip() in ["y", "yes"]
    except EOFError:
        return default


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
            kwargs:
                max_njobs_inqueue:
                    The launcher will stop submitting jobs when the
                    number of jobs in the queue is >= Max number of jobs
        """
        self.flow = flow
        self.max_njobs_inqueue = kwargs.get("max_njobs_inqueue", 200)

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
                logger.critical("numtasks %d already in launched list:\n%s" % (len(task), launched))

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
    This object schedules the submission of the tasks in an :class:`Flow`.
    There are two types of errors that might occur during the execution of the jobs:

        #. Python exceptions
        #. Abinit Errors.

    Python exceptions are easy to detect and are usually due to a bug in abinitio or random errors such as IOError.
    The set of Abinit Errors is much much broader. It includes wrong input data, segmentation
    faults, problems with the resource manager, etc. Abinitio tries to handle the most common cases
    but there's still a lot of room for improvement.
    Note, in particular, that `PyFlowScheduler` will shutdown automatically if

        #. The number of python exceptions is > MAX_NUM_PYEXC

        #. The number of Abinit Errors (i.e. the number of tasks whose status is S_ERROR) is > MAX_NUM_ERRORS

        #. The number of jobs launched becomes greater than (`safety_ratio` * total_number_of_tasks).

        #. The scheduler will send an email to the user (specified by `mailto`) every `remindme_s` seconds.
           If the mail cannot be sent, it will shutdown automatically.
           This check prevents the scheduler from being trapped in an infinite loop.
    """
    # Configuration file.
    YAML_FILE = "scheduler.yml"
    USER_CONFIG_DIR = os.path.join(os.getenv("HOME"), ".abinit", "abipy")

    Error = PyFlowSchedulerError

    def __init__(self, **kwargs):
        """
        Args:
            weeks: number of weeks to wait
            days: number of days to wait
            hours: number of hours to wait
            minutes: number of minutes to wait
            seconds: number of seconds to wait
            verbose: (int) verbosity level
            max_njobs_inque: Limit on the number of jobs that can be present in the queue
            use_dynamic_manager: True if the :class:`TaskManager` must be re-initialized from
                file before launching the jobs. Default: False
            max_nlaunches: Maximum number of tasks launched by radpifire (default -1 i.e. no limit)
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
        self.use_dynamic_manager = kwargs.pop("use_dynamic_manager", False)
        self.max_njobs_inqueue = kwargs.pop("max_njobs_inqueue", 200)
        self.max_ncores_used = kwargs.pop("max_ncores_used", None)
        self.contact_resource_manager = as_bool(kwargs.pop("contact_resource_manager", False))

        self.remindme_s = float(kwargs.pop("remindme_s", 4 * 24 * 3600))
        self.max_num_pyexcs = int(kwargs.pop("max_num_pyexcs", 0))
        self.max_num_abierrs = int(kwargs.pop("max_num_abierrs", 0))
        self.safety_ratio = int(kwargs.pop("safety_ratio", 5))
        #self.max_etime_s = kwargs.pop("max_etime_s", )
        self.max_nlaunches = kwargs.pop("max_nlaunches", -1)
        self.debug = kwargs.pop("debug", 0)

        self.customer_service_dir = kwargs.pop("customer_service_dir", None)
        if self.customer_service_dir is not None:
            self.customer_service_dir = Directory(self.customer_service_dir)
            self._validate_customer_service()

        if kwargs:
            raise self.Error("Unknown arguments %s" % kwargs)

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
        self.exceptions = collections.deque(maxlen=self.max_num_pyexcs + 10)

        # Used to push additional info during the execution.
        self.history = collections.deque(maxlen=100)

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
        return self._flow

    @property
    def num_excs(self):
        """Number of exceptions raised so far."""
        return len(self.exceptions)

    def get_delta_etime(self):
        """Returns a `timedelta` object representing with the elapsed time."""
        return timedelta(seconds=(time.time() - self.start_time))

    def add_flow(self, flow):
        """Add an :class:`Flow` flow to the scheduler."""
        if hasattr(self, "_flow"):
            raise self.Error("Only one flow can be added to the scheduler.")

        if os.path.isfile(flow.pid_file):
            flow.show_status()

            raise self.Error("""\
                pid_file %s already exists
                There are two possibilities:

                   1) There's an another instance of PyFlowScheduler running
                   2) The previous scheduler didn't exit in a clean way

                To solve case 1:
                   Kill the previous scheduler (use 'kill pid' where pid is the number reported in the file)
                   Then you can restart the new scheduler.

                To solve case 2:
                   Remove the pid_file and restart the scheduler.

                Exiting""" % flow.pid_file)

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
        Starts the scheduler in a new thread. Returns True if success.
        In standalone mode, this method will block until there are no more scheduled jobs.
        """
        self.history.append("Started on %s" % time.asctime())
        self.start_time = time.time()

        if has_sched_v3:
            self.sched.add_job(self.callback, "interval", **self.sched_options)
        else:
            self.sched.add_interval_job(self.callback, **self.sched_options)

        errors = self.flow.look_before_you_leap()
        if errors:
            self.exceptions.append(errors)
            return False

        # Try to run the job immediately. If something goes wrong return without initializing the scheduler.
        self._runem_all()

        if self.exceptions:
            self.cleanup()
            self.send_email(msg="Error while trying to run the flow for the first time!\n %s" % self.exceptions)
            return False

        try:
            self.sched.start()
            return True

        except KeyboardInterrupt:
            self.shutdown(msg="KeyboardInterrupt from user")
            if ask_yesno("Do you want to cancel all the jobs in the queue? [Y/n]"): 
                print("Number of jobs cancelled %s", self.flow.cancel())

            self.flow.pickle_dump()
            return False

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
            from pymatgen.io.abinitio.tasks import TaskManager
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
                logger.info("Too many jobs in the queue, returning")
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
            print("Cannot exceed max_ncores_used %s" % self.max_ncores_used)
            logger.info("Cannot exceed max_ncores_used %s" % self.max_ncores_used)
            return

        # fix problems
        # Try to restart the unconverged tasks
        # todo donot fire here but prepare for fireing in rapidfire
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

        # move here from withing rapid fire ...
        # fix only prepares for restarting, and sets to ready
        nfixed = flow.fix_abi_critical()
        if nfixed: print("Fixed %d AbiCritical errors" % nfixed)

        # Temporarily disable by MG because I don't know if fix_critical works after the
        # introduction of the new qadapters
        if False:
            nfixed = flow.fix_queue_critical()
            if nfixed: print("Fixed %d QueueCritical errors" % nfixed)

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
            self.exceptions.append(straceback())
            self.shutdown(msg="Exception raised in callback!")

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
        err_msg = ""

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
                err_msg += msg

        #if delta_etime.total_seconds() > self.max_etime_s:
        #    err_msg += "\nExceeded max_etime_s %s. Will shutdown the scheduler and exit" % self.max_etime_s

        # Too many exceptions. Shutdown the scheduler.
        if self.num_excs > self.max_num_pyexcs:
            msg = "Number of exceptions %s > %s. Will shutdown the scheduler and exit" % (
                self.num_excs, self.max_num_pyexcs)
            err_msg += boxed(msg)

        # Paranoid check: disable the scheduler if we have submitted
        # too many jobs (it might be due to some bug or other external reasons 
        # such as race conditions between difference callbacks!)
        if self.nlaunch > self.safety_ratio * self.flow.num_tasks:
            msg = "Too many jobs launched %d. Total number of tasks = %s, Will shutdown the scheduler and exit" % (
                self.nlaunch, self.flow.num_tasks)
            err_msg += boxed(msg)

        # Count the number of tasks with status == S_ERROR.
        if self.flow.num_errored_tasks > self.max_num_abierrs:
            msg = "Number of tasks with ERROR status %s > %s. Will shutdown the scheduler and exit" % (
                self.flow.num_errored_tasks, self.max_num_abierrs)
            err_msg += boxed(msg)

        # TODO: This is not correct because I should also check tasks that are submitted
        # Test on the presence of deadlocks.
        """
        deadlocked, runnables, running = self.flow.deadlocked_runnables_running()
        if deadlocked: 
            # Check the flow agains to that status are updated. 
            self.flow.check_status()

            deadlocked, runnables, running = self.flow.deadlocked_runnables_running()
            print("deadlocked:\n", deadlocked, "\nrunnables:\n", runnables, "\nrunning\n", running)
            if deadlocked and not runnables and not running:
                err_msg += "No runnable job with deadlocked tasks:\n%s." % str(deadlocked)

        if not runnables and not running:
            # Check the flow again to that status are updated. 
            self.flow.check_status()
            deadlocked, runnables, running = self.flow.deadlocked_runnables_running()
            if not runnables and not running:
                err_msg += "No task is running and cannot find other tasks to sumbmit."
        """

        if err_msg:
            # Something wrong. Quit
            self.shutdown(err_msg)

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

            self.history.append("Completed on %s" % time.asctime())
            self.history.append("Elapsed time %s" % self.get_delta_etime())

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
            app("Submitted on %s" % time.ctime(self.start_time))
            app("Completed on %s" % time.asctime())
            app("Elapsed time %s" % str(self.get_delta_etime()))
            if self.flow.all_ok:
                app("Flow completed successfully")
            else:
                app("Flow didn't complete successfully")
                app("Shutdown message:\n%s" % msg)

            print("")
            print("\n".join(lines))

            self._do_customer_service()

        finally:
            # Shutdown the scheduler thus allowing the process to exit.
            logger.debug('this should be the shutdown of the scheduler')

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

        app("Submitted on %s" % time.ctime(self.start_time))
        app("Completed on %s" % time.asctime())
        app("Elapsed time %s" % str(self.get_delta_etime()))
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
    sender = user_at_host() if sender is None else sender
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


class BatchLauncher(object):

    PICKLE_FNAME = "__BatchLauncher__.pickle"

    @classmethod
    def from_dir(cls, top, max_depth=2):
        """
        Find all flows located withing the directory `top` and build the `BatchLauncher`.

        Args:
            max_depth: Search in directory only if it is N or fewer levels below top
        """
        # Walk through each directory inside path and find the pickle database.
        from .flows import Flow
        for dirpath, dirnames, filenames in os.walk(top):
            fnames = [f for f in filenames if f == Flow.PICKLE_FNAME]

        #new = cls(workdir, name=None, flows=None, manager=None, timelimit=1200)
        #return new

    @classmethod
    def pickle_load(self, filepath):
        """
        Loads the object from a pickle file and performs initial setup.

        Args:
            filepath: Filename or directory name. It filepath is a directory, we
                scan the directory tree starting from filepath and we
                read the first pickle database. Raise RuntimeError if multiple
                databases are found.
        """
        with open(filepath, "rb") as fh:
            new = pickle.load(fh)
        return new

    def pickle_dump(self):
        """
        Save the status of the object in pickle format.
        """
        import pickle
        with open(os.path.join(self.workdir, self.PICKLE_FNAME), mode="wb") as fh:
            pickle.dump(self, fh)

    def __init__(self, workdir, name=None, flows=None, manager=None):
        """
        Args:
            workdir: Working directory
            name: Name assigned to the `BatchLauncher`.
            manager: :class:`TaskManager` object responsible for the submission of the jobs.
                     If manager is None, the object is initialized from the yaml file
                     located either in the working directory or in the user configuration dir.
            flows:  List of `Flow` objects.
        """
        self.workdir = os.path.abspath(workdir)

        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        else:
            #raise RuntimeError("Directory %s already esists. Use BatchLauncher.pickle_load()" % self.workdir)
            pass

        self.name = os.path.basename(self.workdir) if name is None else name
        self.qerr_file = File(os.path.join(self.workdir, "queue.qerr"))
        self.qout_file = File(os.path.join(self.workdir, "queue.qout"))

        from .tasks import TaskManager
        manager = TaskManager.as_manager(manager)

        # Extract the qadapater to be used for the batch script.
        self.qadapter = qad = manager.qads[0]
        qad.set_mpi_procs(1)
        qad.set_timelimit(36000)

        # Initialize list of flows.
        if flows is None: flows = []
        if not isinstance(flows, (list, tuple)): flows = [flows]
        self.flows = flows

    def __str__(self):
        lines = []
        lines.extend(str(self.qadapter).splitlines())

        for i, flow in enumerate(self.flows):
            lines.append("Flow [%d] " % i + str(flow))

        return "\n".join(lines)

    def add_flow(self, flow):
        """Add a flow. Accept filepath or :class:`Flow` object."""
        from .flows import Flow
        flow = Flow.as_flow(flow)

        if flow in self.flows:
            raise ValueError("Cannot add same flow twice!")

        # Check if we are already using a scheduler to run this flow
        if os.path.exists(flow.pid_file):
            msg = ("Pid file %s already exists.\n" % flow.pid_file +
                   "Either the flow is already in execution or the scheduler didn't exit in a clean way.\n" +
                   "Make sure no scheduler is running and then remove the pid file")
            raise RuntimeError(flow.pid_file)

        self.flows.append(flow)

    def submit(self, verbose=0):
        if not self.flows:
            raise RuntimeError("Cannot submit an empty list of flows!")

        if hasattr(self, "qjob"):
            msg = "Got qjob %s" % self.qjob
            raise RuntimeError(msg)

        script = self._get_job_script_str()
        if verbose:
            print("*** submission script ***")
            print(script)

        # Write the script.
        script_file = os.path.join(self.workdir, "run.sh")
        with open(script_file, "wt") as fh:
            fh.write(script)
            os.chmod(script_file, 0o740)

        # Submit the task and save the queue id.
        self.qjob, process = self.qadapter.submit_to_queue(script_file)
        #process.wait()
        self.pickle_dump()

        return process.returncode

    def _get_job_script_str(self):
        """Write the submission script. return the path of the script"""
        executable = []
        app = executable.append

        # Build list of abirun commands and save the name of the log files.
        self.sched_logs = []
        for i, flow in enumerate(self.flows):
            logfile = "log_flow_" + os.path.basename(flow.workdir)
            app("abirun.py %s scheduler > %s" % (flow.workdir, logfile))
            assert logfile not in self.sched_logs
            self.sched_logs.append(logfile)

        script = self.qadapter.get_script_str(
            job_name=self.name, 
            launch_dir=self.workdir,
            executable=executable,
            qout_path=self.qout_file.path,
            qerr_path=self.qerr_file.path,
        )

        return script
