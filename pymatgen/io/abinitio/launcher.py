"""Tools for the submission of Tasks."""
from __future__ import division, print_function

import os 
import time
import collections
import yaml
import cStringIO as StringIO

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
                    # No task found, this usually happens when we have dependencies. 
                    # Beware of possible deadlocks here!
                    logger.debug("No task to run! Possible deadlock")

            except StopIteration:
                logger.info("All tasks completeeed.")

        # Submit the tasks and update the database.
        if tasks:
            tasks[0].start()
            num_launched += 1
            
            self.flow.pickle_dump()

        return num_launched 

    def rapidfire(self, max_nlaunch=-1, max_loops=10, sleep_time=None): # nlaunches=0, 
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
                    err_msg = "task %s in already in launched:\n%s" % (task, launched)
                    logger.warning(err_msg)
                    #raise RuntimeError(err_msg)

            # Preventive test.
            tasks = [t for t in tasks if t not in launched]

            if not tasks:
                got_empty_list += 1

            if got_empty_list == max_loops:
                break

            for task in tasks:
                task.start()
                launched.append(task)
                num_launched += 1
                                                
                if num_launched == max_nlaunch:
                    # Exit the outermst loop.
                    num_loops = max_loops
                    break

            #time.sleep(sleep_time)
            num_loops += 1

        #for work in self.flow:
        #    if num_launched == max_nlaunch: 
        #        break

        #    while True:
        #        try:
        #            task = work.fetch_task_to_run()

        #            if task is None:
        #                logger.debug("fetch_task_to_run returned None.")
        #                break

        #            logger.debug("Starting task %s" % task)
        #            if task in launched:
        #                err_msg = "task %s in already in launched %s" % (str(task), str(launched))
        #                raise RuntimeError(err_msg)

        #            task.start()
        #            launched.append(task)
        #            num_launched += 1

        #            if num_launched == max_nlaunch:
        #                break

        #        except StopIteration:
        #            logger.debug("Out of tasks.")
        #            break

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
            try:
                task = work.fetch_task_to_run()
                                                                                                   
                if task is not None:
                    tasks_to_run.append(task)
                else:
                    # No task found, this usually happens when we have dependencies. 
                    # Beware of possible deadlocks here!
                    logger.debug("fetch_task_to_run returned None.")

            except StopIteration:
                # All the tasks in work are done.
                logger.debug("Out of tasks in work %s" % work)

        return tasks_to_run



class PyFlowsScheduler(object):
    """
    This object schedules the submission of the tasks in an `AbinitFlow`.
     There are two types of errors that migh happen during the execution of the jobs: 
    
        #. Python exceptions 
        #. Abinit Errors.
    
    Python exceptions are easy to detect and are usually due to a bug in abinitio or random errors such as IOError.
    The set of Abinit Error is much much broader. It includes wrong input data, segmentation
    faults, problem with the resource manager, etc. Abinitio tries to handle the most common cases
    but there's still a lot of room for improvement.
    The scheduler will shutdown automatically if the number of python exceptions is > MAX_NUM_PYEXC
    or if the number of Abinit Errors (i.e. the number of tasks whose status is S_ERROR) is > MAX_NUM_ERRORS

    There's another check that prevents the scheduler from running infinitely, 
    If the number of launches becomes greater than SAFETY_RATIO * number_of_tasks, the scheduler
    will automatically shutdow.

    TODO: detection of Deadlocks NLAUNCH_FLOW
    """
    MAX_NUM_PYEXCS = 2
    MAX_NUM_ABIERRS = 0

    SAFETY_RATIO = 3
    #MAX_ETIME = 

    # Configuration file.
    YAML_FILE = "scheduler.yml"

    def __init__(self, weeks=0, days=0, hours=0, minutes=0, seconds=0, start_date=None, mailto=None):
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
             start_date=start_date)

        self.mailto = mailto
        self.nlaunch = 0

        from apscheduler.scheduler import Scheduler
        self.sched = Scheduler(standalone=True)

        # Used to keep track of the exceptions raised 
        # while the scheduler is running
        self.exceptions = collections.deque(maxlen=self.MAX_NUM_PYEXCS + 10)

        # Used to push additional info during the execution. 
        self.history = collections.deque(maxlen=100)

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

        app(80 * "=")
        app(str(self.flow))

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
    def pid_file(self):
        """Absolute path of the file with the pid. The file is located in the workdir of the flow"""
        return self._pidfile
                                                                                            
    @property
    def flow(self):
        """`AbinitFlows`."""
        return self._flow

    @property
    def num_excs(self):
        return len(self.exceptions)

    def add_flow(self, flow):
        """Add an `AbinitFlow` to the scheduler."""
        pid_file = os.path.join(flow.workdir, "_PyFlowsScheduler.pid")

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

        with open(pid_file, "w") as fh:
            fh.write(str(self.pid))

        self._pid_file = pid_file
        self._flow = flow

    def start(self):
        """
        Starts the scheduler in a new thread.
        In standalone mode, this method will block until there are no more scheduled jobs.
        """
        self.sched.add_interval_job(self.callback, **self.sched_options)

        # Try to run the job immediately. If something goes wrong 
        # return without initializing the scheduler.
        retcode = self._runem_all()
        if retcode:
            self.cleanup()
            return retcode

        self.history.append("Started on %s" % time.asctime())
        self.start_time = time.time()

        self.sched.start()

    def _runem_all(self):
        max_nlaunch, excs = -1, []

        self.flow.check_status()

        try:
            nlaunch = PyLauncher(self.flow).rapidfire(max_nlaunch=max_nlaunch)
            print("Number of launches: ", nlaunch)
            self.nlaunch += nlaunch

        except Exception:
            excs.append(straceback())

        self.flow.show_status()
        
        if excs:
            logger.critical("*** Scheduler exceptions:\n *** %s" % "\n".join(excs))
            self.exceptions.extend(excs)

    def callback(self):
        """The function that will be executed by the scheduler."""
        try:
            return self._callback()
        except:
            self.exceptions.append(straceback())
            return self.num_excs

    def _callback(self):
        """The actual callback."""
        self._runem_all()

        # Too many exceptions. Shutdown the scheduler.
        if self.num_excs > self.MAX_NUM_PYEXCS:
            msg = "Number of exceptions %s > %s.\nWill shutdown the scheduler and exit" % (self.num_excs, self.MAX_NUM_PYEXCS)
            msg = 5*"=" + msg + 5*"="
                                                                                      
            print(len(msg)* "=")
            print(msg)
            print(len(msg)* "=")

            self.shutdown()

        # Paranoid check: disable the scheduler if we have submitted 
        # too many jobs (it might be due to some bug or other external reasons!) 
        if self.nlaunch > self.SAFETY_RATIO * self.flow.num_tasks:
            msg = "Too many jobs launched %d.\nWill shutdown the scheduler and exit" % self.nlaunch
            msg = 5*"=" + msg + 5*"="
                                                                                      
            print(len(msg)* "=")
            print(msg)
            print(len(msg)* "=")
                                                                                                                               
            self.shutdown()

        # Too long elapsed time
        #if time.time() - self.start_time) > MAX_ETIME:
        #    self.shutdown()

        # Count the number of tasks with status == S_ERROR.
        num_errors = self.flow.num_tasks_with_error

        # Count the number of tasks with status == S_UNCONVERGED.
        num_unconverged = self.flow.num_tasks_unconverged

        message = ""
        if num_errors > self.MAX_NUM_ABIERRS:
            err_msg = "Number of task with ERROR status %s > %s. Will shutdown the scheduler and exit" % (
                num_errors, self.MAX_NUM_ABIERRS)
            err_msg = 5*"=" + msg + 5*"="

            message += err_msg

        if num_unconverged:
            # TODO: this is needed to avoid deadlocks, automatic restarting is not available yet
            err_msg = ("Found %d unconverged tasks." 
                "Automatic restarting is not available yet. Will shutdown the scheduler and exit" 
                % num_unconverged)
            err_msg = 5*"=" + err_msg + 5*"="

            message += "\n" + err_msg

        if message:
            # Something wrong. Quit
            print(len(message)* "=")
            print(message)
            print(len(message)* "=")
                                                                                            
            self.shutdown()

        # Mission accomplished. Shutdown the scheduler.
        if self.flow.all_ok:
            message = "All tasks have reached S_OK. Will shutdown the scheduler and exit"
            message = 5*"=" + message + 5*"="

            print(len(message)* "=")
            print(message)
            print(len(message)* "=")

            self.shutdown()

        return len(self.exceptions) 

    def cleanup(self):
        """
        Cleanup routine: remove the pid files and save the pickle database
        """
        try:
            os.unlink(self.pid_file)
        except:
            pass
                                                                    
        # Save the final status of the flow.
        self.flow.pickle_dump()

    def shutdown(self):
        """Shutdown the scheduler."""
        self.cleanup()

        from datetime import timedelta
        self.etime = timedelta(seconds=(time.time() - self.start_time))

        self.history.append("Completed on %s" % time.asctime())
        self.history.append("Elapsed time %s" % str(self.etime))
        print("history", self.history)

        try:
            self.send_email()
        except:
            exc_msg = straceback()
            print("exc_msg", exc_msg)
            self.exceptions.append(exc_msg)

        # Write text file with the list of exceptions:
        if self.exceptions:
            dump_file = os.path.join(self.flow.workdir, "launcher.log")
            with open(dump_file,"w") as fh:
                fh.writelines(self.exceptions)

        # Shutdown the scheduler thus allowing the process to exit.
        self.sched.shutdown(wait=False)

    def send_email(self):
        """Send an e-mail before terminating the shutdown."""
        if self.mailto is None: return 0
        print("mailto", self.mailto)

        # The status of the flow.
        strio = StringIO.StringIO()
        self.flow.show_status(stream=strio)

        #self.history.append("Elapsed time %s" % str(self.etime))
        #self.history.append("Submitted on %s" % time.asctime())
        #self.start_time = time.time()
        #self.history.append("Completed on %s" % time.asctime())
        #self.history.append("Elapsed time %s" % str(self.etime)))

        strio.write("Number of exceptions: %d" % self.num_excs)
        if self.exceptions:
            strio.writelines(self.exceptions)

        strio.seek(0)
        text = strio.read()
        print("text", text)

        return sendmail(subject=self.flow.name,
                        text=text,
                        mailto=self.mailto,
                        sender=None)


def sendmail(subject, text, mailto, sender=None, test_sendmail=False):
    """
    Sends an e-mail either with smtplib.SMTP or with the unix sendmail.

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
    """
    print("in sendmail")
    def user_at_host():
        from socket import gethostname
        return os.getlogin() + "@" + gethostname()

    import smtplib
    from email.mime.text import MIMEText

    # Body of the message
    sender = user_at_host() if sender is None else sender 
    if is_string(mailto): mailto = [mailto]

    mail = MIMEText(text)
    mail["Subject"] = subject
    mail["From"] = sender
    mail["To"] = ", ".join(mailto)

    msg = mail.as_string()
    print("msg", msg)

    try:
        if test_sendmail:
            a += 1 
        # Send the message via our own SMTP server.
        server = smtplib.SMTP("localhost")
        server.sendmail(sender, mailto, msg)
        return server.quit()

    except:
        # Fallback to sendmail (Unix only).
        print("Using sendmail")
        from subprocess import Popen, PIPE
        SENDMAIL = "/usr/sbin/sendmail" 
        p = Popen([SENDMAIL, "-t"], stdin=PIPE, stderr=PIPE)

        outdata, errdata = p.communicate(msg)
        print(outdata, errdata)
        return len(errdata)


# Test for sendmail
if __name__ == "__main__":
    text = \
"""hello
world
"""
    mailto = "matteo.giantomassi@uclouvain.be"
    retcode = sendmail("sendmail_test", text, mailto, sender=None, test_sendmail=False)
    print("Retcode", retcode)

