"""Tools for the submission of Tasks."""
from __future__ import division, print_function

import os
import time
import collections
import yaml
import cStringIO as StringIO
from monty.os.path import which
from pymatgen.io.abinitio import myaml

from datetime import timedelta
from pymatgen.core.design_patterns import AttrDict
from pymatgen.util.string_utils import is_string
from pymatgen.io.abinitio.utils import File

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

    @classmethod
    def from_file(cls, filename, allow_empty=False):
        """Reads the OpenMP variables from a INI file."""
        if filename.endswith(".ini"):
            from ConfigParser import SafeConfigParser, NoOptionError

            parser = SafeConfigParser()
            parser.read(filename)

            obj = OmpEnv()

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

            for key in cls._KEYS:
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

    def rapidfire(self, max_nlaunch=-1, max_loops=1, sleep_time=None): # nlaunches=0,
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
        sleep_time = sleep_time if sleep_time else 5
        #nlaunches = -1 if nlaunches == 'infinite' else int(nlaunches)
        num_loops, num_launched, launched = 0, 0, []
        got_empty_list = 0

        while num_loops != max_loops:
            tasks = self.fetch_tasks_to_run()

            # I don't know why but we receive duplicated tasks.
            for task in tasks:
                if task in launched:
                    err_msg = "task %s already in launched list:\n%s" % (task, launched)
                    logger.critical(err_msg)

            # Preventive test.
            tasks = [t for t in tasks if t not in launched]

            if not tasks:
                got_empty_list += 1

            if got_empty_list == max_loops:
                break

            for task in tasks:
                fired = task.start()

                if fired:
                    launched.append(task)
                    num_launched += 1

                if num_launched == max_nlaunch:
                    # Exit the outermst loop.
                    num_loops = max_loops
                    break

            #time.sleep(sleep_time)
            num_loops += 1

        # Update the database.
        self.flow.check_status()
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
            #try:
            #    task = work.fetch_task_to_run()

            #    if task is not None:
            #        tasks_to_run.append(task)
            #    else:
            #        # No task found, this usually happens when we have dependencies.
            #        # Beware of possible deadlocks here!
            #        logger.debug("fetch_task_to_run returned None.")

            #except StopIteration:
            #    # All the tasks in work are done.
            #    logger.debug("Out of tasks in work %s" % work)

        return tasks_to_run


def boxed(msg, ch="=", pad=5):
    if pad > 0:
        msg = pad * ch + msg + pad * ch

    return "\n".join([len(msg) * ch,
                      msg,
                      len(msg) * ch,
                      "", ])


class PyFlowScheduler(object):
    """
    This object schedules the submission of the tasks in an `AbinitFlow`.
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

        #. The number of jobs launched becomes greater than (SAFETY_RATIO * total_number_of_tasks).

        #. The scheduler will send an email to the user (specified by mailto) every REMINDME_S seconds.
           If the mail cannot be sent, it will shutdown automatically.
           This check prevents the scheduler from being trapped in an infinite loop.
    """
    # Configuration file.
    YAML_FILE = "scheduler.yml"

    def __init__(self, **kwargs):
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
            verbose:
                (int) verbosity level
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
            raise ValueError("Wrong set of options passed to the scheduler.")

        self.mailto = kwargs.pop("mailto", None)
        self.verbose = int(kwargs.pop("verbose", 0))

        self.REMINDME_S = float(kwargs.pop("REMINDME_S", 24 * 3600))
        self.MAX_NUM_PYEXCS = int(kwargs.pop("MAX_NUM_PYEXCS", 0))
        self.MAX_NUM_ABIERRS = int(kwargs.pop("MAX_NUM_ABIERRS", 0))
        self.SAFETY_RATIO = int(kwargs.pop("SAFETY_RATIO", 3))
        #self.MAX_ETIME_S = kwargs.pop("MAX_ETIME_S", )

        if kwargs:
            raise ValueError("Unknown arguments %s" % kwargs)

        from apscheduler.scheduler import Scheduler

        self.sched = Scheduler(standalone=True)

        self.nlaunch = 0
        self.num_reminders = 1

        # Used to keep track of the exceptions raised while the scheduler is running
        self.exceptions = collections.deque(maxlen=self.MAX_NUM_PYEXCS + 10)

        # Used to push additional info during the execution.
        self.history = collections.deque(maxlen=100)

    @classmethod
    def from_file(cls, filepath):
        """Read the configuration parameters from a Yaml file."""
        with open(filepath, "r") as fh:
            d = myaml.load(fh)
            return cls(**d)

    @classmethod
    def from_user_config(cls):
        """
        Initialize the `PyFlowScheduler` from the YAML file 'scheduler.yml'.
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

        err_msg = "Cannot locate %s neither in current directory nor in %s" % (cls.YAML_FILE, dirpath)
        raise RuntimeError(err_msg)

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
        """`AbinitFlow`."""
        return self._flow

    @property
    def num_excs(self):
        """Number of exceptions raised so far."""
        return len(self.exceptions)

    def get_delta_etime(self):
        """Returns a `timedelta` object representing with the elapsed time."""
        return timedelta(seconds=(time.time() - self.start_time))

    def add_flow(self, flow):
        """Add an `AbinitFlow` to the scheduler."""
        if hasattr(self, "_flow"):
            raise ValueError("Only one flow can be added to the scheduler.")

        pid_file = os.path.join(flow.workdir, "_PyFlowScheduler.pid")

        if os.path.isfile(pid_file):
            flow.show_status()

            err_msg = (
                "pid_file %s already exists\n"
                "There are two possibilities:\n\n"
                "   1) There's an another instance of PyFlowScheduler running.\n"
                "   2) The previous scheduler didn't exit in a clean way.\n\n"
                "To solve case 1:\n"
                "   Kill the previous scheduler (use 'kill pid' where pid is the number reported in the file)\n"
                "   Then you run can start the new scheduler.\n\n"
                "To solve case 2:\n"
                "   Remove the pid_file and restart the scheduler.\n\n"
                "Exiting\n" % pid_file
            )

            raise RuntimeError(err_msg)

        with open(pid_file, "w") as fh:
            fh.write(str(self.pid))

        self._pid_file = pid_file
        self._flow = flow

    def start(self):
        """
        Starts the scheduler in a new thread.
        In standalone mode, this method will block until there are no more scheduled jobs.
        """
        self.history.append("Started on %s" % time.asctime())
        self.start_time = time.time()

        self.sched.add_interval_job(self.callback, **self.sched_options)

        # Try to run the job immediately. If something goes wrong
        # return without initializing the scheduler.
        retcode = self._runem_all()

        if retcode:
            self.cleanup()
            self.send_email(msg="Error while trying to run the flow for the first time!")

        self.sched.start()

    def _runem_all(self):
        """
        This function tries to run all the tasks that can be submitted.
        """
        max_nlaunch, excs = -1, []

        flow = self.flow
        flow.check_status()

        # Try to restart the unconverged tasks
        for task in self.flow.unconverged_tasks:
            try:
                msg = "AbinitFlow will try restart task %s" % task
                print(msg)
                logger.info(msg)

                fired = task.restart()
                if fired: self.nlaunch += 1

            except Exception:
                excs.append(straceback())

        flow.pickle_dump()

        #if self.num_restarts == self.max_num_restarts:
        #    info_msg = "Reached maximum number of restarts. Cannot restart anymore Returning"
        #    logger.info(info_msg)
        #    self.history.append(info_msg)
        #    return 1

        # Submit the tasks that are ready.
        try:
            nlaunch = PyLauncher(flow).rapidfire(max_nlaunch=max_nlaunch)
            self.nlaunch += nlaunch

            if nlaunch:
                print("[%s] Number of launches: %d" % (time.asctime(), nlaunch))

        except Exception:
            excs.append(straceback())

        show_status = (self.verbose or flow.num_errored_tasks or flow.num_unconverged_tasks)
        if show_status: flow.show_status()

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
        self._runem_all()

        # Mission accomplished. Shutdown the scheduler.
        all_ok = self.flow.all_ok
        if self.verbose: print("all_ok", all_ok)
        if all_ok:
            self.shutdown(msg="All tasks have reached S_OK. Will shutdown the scheduler and exit")

        # Handle failures.
        err_msg = ""

        # Shall we send a reminder to the user?
        delta_etime = self.get_delta_etime()

        if delta_etime.total_seconds() > self.num_reminders * self.REMINDME_S:
            self.num_reminders += 1
            msg = ("Just to remind you that the scheduler with pid %s, flow %s\n has been running for %s " %
                   (self.pid, self.flow, delta_etime))
            retcode = self.send_email(msg, tag="[REMINDER]")

            if retcode:
                # Cannot send mail, shutdown now!
                msg += "\nThe scheduler tried to send an e-mail to remind user but send_email returned %d. Aborting now" % retcode
                err_msg += msg

        #if delta_etime.total_seconds() > self.MAX_ETIME_S:
        #    err_msg += "\nExceeded MAX_ETIME_S %s. Will shutdown the scheduler and exit" % self.MAX_ETIME_S

        # Too many exceptions. Shutdown the scheduler.
        if self.num_excs > self.MAX_NUM_PYEXCS:
            msg = "Number of exceptions %s > %s. Will shutdown the scheduler and exit" % (
            self.num_excs, self.MAX_NUM_PYEXCS)
            err_msg += boxed(msg)

        # Paranoid check: disable the scheduler if we have submitted
        # too many jobs (it might be due to some bug or other external reasons such as race conditions betwee difference callbacks.!)
        if self.nlaunch > self.SAFETY_RATIO * self.flow.num_tasks:
            msg = "Too many jobs launched %d. Total number of tasks = %s, Will shutdown the scheduler and exit" % (
            self.nlaunch, self.flow.num_tasks)
            err_msg += boxed(msg)

        # Count the number of tasks with status == S_ERROR.
        if self.flow.num_errored_tasks > self.MAX_NUM_ABIERRS:
            msg = "Number of tasks with ERROR status %s > %s. Will shutdown the scheduler and exit" % (
                self.flow.num_errored_tasks, self.MAX_NUM_ABIERRS)
            err_msg += boxed(msg)

        # Count the number of tasks with status == S_UNCONVERGED.
        #if self.flow.num_unconverged_tasks:
        #    # TODO: this is needed to avoid deadlocks, automatic restarting is not available yet
        #    msg = ("Found %d unconverged tasks."
        #           "Automatic restarting is not available yet. Will shutdown the scheduler and exit"
        #           % self.flow.num_unconverged_tasks)
        #    err_msg += boxed(msg)

        #deadlocks = self.detect_deadlocks()
        #if deadlocks:
        #    msg = ("Detected deadlocks in flow. Will shutdown the scheduler and exit"
        #           % self.flow.num_unconverged_tasks)
        #    err_msg += boxed(msg)

        if err_msg:
            # Something wrong. Quit
            self.shutdown(err_msg)

        return len(self.exceptions)

    def cleanup(self):
        """
        Cleanup routine: remove the pid file and save the pickle database
        """
        try:
            os.remove(self.pid_file)
        except OSError:
            logger.critical("Could not remove pid_file")
            pass

        # Save the final status of the flow.
        self.flow.pickle_dump()

    def shutdown(self, msg):
        """Shutdown the scheduler."""
        try:
            self.cleanup()

            self.history.append("Completed on %s" % time.asctime())
            self.history.append("Elapsed time %s" % self.get_delta_etime())

            retcode = self.send_email(msg)
            print("send_mail retcode", retcode)

        finally:
            # Write file with the list of exceptions:
            if self.exceptions:
                dump_file = os.path.join(self.flow.workdir, "_exceptions")
                with open(dump_file, "w") as fh:
                    fh.writelines(self.exceptions)

            # Shutdown the scheduler thus allowing the process to exit.
            self.sched.shutdown(wait=False)

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

        strio = StringIO.StringIO()
        strio.writelines("\n".join(header) + 4 * "\n")

        # Add the status of the flow.
        self.flow.show_status(stream=strio)

        if self.exceptions:
            # Report the list of exceptions.
            strio.writelines(self.exceptions)

        text = strio.getvalue()
        #print("text", text)

        if tag is None:
            tag = " [ALL OK]" if self.flow.all_ok else " [WARNING]"

        return sendmail(subject=self.flow.name + tag, text=text, mailto=self.mailto)


def sendmail(subject, text, mailto, sender=None):
    """
    Sends an e-mail either with unix sendmail.

    Args:
        subject:
            String with the subject of the mail.
        text:
            String with the body of the mail.
        mailto:
            String or list of string with the recipients.
        sender:
            string with the sender address.
            If sender is None, username@hostname is used.

    Returns:
        exit status
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


# Test for sendmail
#def sendmail_test():
#    text = "hello\nworld"""
#    mailto = "matteo.giantomassi@uclouvain.be"
#    retcode = sendmail("sendmail_test", text, mailto)
#    print("Retcode", retcode)
