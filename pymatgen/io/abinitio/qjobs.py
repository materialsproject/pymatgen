# coding: utf-8
"""
Objects and methods to contact the resource manager to get info on the status of the job and useful statistics. 
Note that this is not a wrapper for the C API but a collection of simple wrappers around the shell commands 
provided by the resource manager  (qsub, qdel and qstat for PBS, sinfo, squeue... for Slurm).
The main goal indeed is providing a simplified common interface for different resource managers without
having to rely on external libraries.
"""
from __future__ import print_function, division, unicode_literals

import shlex

from collections import OrderedDict, defaultdict
from subprocess import Popen, PIPE
from monty.collections import AttrDict
from monty.inspect import all_subclasses

import logging
logger = logging.getLogger(__name__)


class JobStatus(int):
    """
    This object is an integer representing the status of a :class:`QueueJob`.

    Slurm API, see `man squeue`.

    JOB STATE CODES
       Jobs typically pass through several states in the course of their execution.  The typical  states  are
       PENDING, RUNNING, SUSPENDED, COMPLETING, and COMPLETED.  An explanation of each state follows.

    BF  BOOT_FAIL       Job  terminated  due  to launch failure, typically due to a hardware failure (e.g.
                        unable to boot the node or block and the job can not be requeued).
    CA  CANCELLED       Job was explicitly cancelled by the user or system administrator.
                        The job may or may not have been initiated.
    CD  COMPLETED       Job has terminated all processes on all nodes.
    CF  CONFIGURING     Job has been allocated resources, but are waiting for them to become ready for use (e.g. booting).
    CG  COMPLETING      Job is in the process of completing. Some processes on some  nodes may still be active.
    F   FAILED          Job terminated with non-zero exit code or other failure condition.
    NF  NODE_FAIL       Job terminated due to failure of one or more allocated nodes.
    PD  PENDING         Job is awaiting resource allocation.
    PR  PREEMPTED       Job terminated due to preemption.
    R   RUNNING         Job currently has an allocation.
    S   SUSPENDED       Job has an allocation, but execution has been suspended.
    TO  TIMEOUT         Job terminated upon reaching its time limit.
    SE SPECIAL_EXIT     The job was requeued in a special state. This state can be set by users, typically
                        in EpilogSlurmctld, if the job has terminated with a particular exit value.
    """

    _STATUS_TABLE = OrderedDict([
        (-1, "UNKNOWN"),
        (0, "PENDING"),
        (1, "RUNNING"),
        (2, "RESIZING"),
        (3, "SUSPENDED"),
        (4, "COMPLETED"),
        (5, "CANCELLED"),
        (6, "FAILED"),
        (7, "TIMEOUT"),
        (8, "PREEMPTED"),
        (9, "NODEFAIL"),
    ])

    def __repr__(self):
        return "<%s: %s, at %s>" % (self.__class__.__name__, str(self), id(self))

    def __str__(self):
        """String representation."""
        return self._STATUS_TABLE[self]

    @classmethod
    def from_string(cls, s):
        """Return a :class:`JobStatus` instance from its string representation."""
        for num, text in cls._STATUS_TABLE.items():
            if text == s: return cls(num)
        else:
            #raise ValueError("Wrong string %s" % s)
            logger.warning("Got unknown status: %s" % s)
            return cls.from_string("UNKNOWN")


class QueueJob(object):
    """
    This object provides methods to contact the resource manager to get info on the status
    of the job and useful statistics. This is an abstract class.
    """
    QTYPE = None

    # Used to handle other resource managers.
    S_UNKNOWN   = JobStatus.from_string("UNKNOWN")
    # Slurm status
    S_PENDING   = JobStatus.from_string("PENDING")
    S_RUNNING   = JobStatus.from_string("RUNNING")
    S_RESIZING  = JobStatus.from_string("RESIZING")
    S_SUSPENDED = JobStatus.from_string("SUSPENDED")
    S_COMPLETED = JobStatus.from_string("COMPLETED")
    S_CANCELLED = JobStatus.from_string("CANCELLED")
    S_FAILED    = JobStatus.from_string("FAILED")
    S_TIMEOUT   = JobStatus.from_string("TIMEOUT")
    S_PREEMPTED = JobStatus.from_string("PREEMPTED")
    S_NODEFAIL  = JobStatus.from_string("NODEFAIL")

    @staticmethod
    def from_qtype_and_id(qtype, queue_id, qname=None):
        """
        Return a new istance of the appropriate subclass.

        Args:
            qtype: String specifying the Resource manager type.
            queue_id: Job identifier.
            qname: Name of the queue (optional).
        """
        for cls in all_subclasses(QueueJob):
            if cls.QTYPE == qtype: break
        else:
            logger.critical("Cannot find QueueJob subclass registered for qtype %s" % qtype)
            cls = QueueJob

        return cls(queue_id, qname=qname)

    def __init__(self, queue_id, qname="UnknownQueue"):
        """
        Args:
            queue_id: Job identifier.
            qname: Name of the queue (optional).
        """
        self.qid, self.qname = queue_id, qname

        # Initialize properties.
        self.status, self.exitcode, self.signal = None, None, None

    def __repr__(self):
        return "<%s, qid=%s, status=%s, exit_code=%s>" % (
            self.__class__.__name__, self.qid, self.status, self.exitcode)

    def __bool__(self):
        return self.qid is not None

    __nonzero__ = __bool__

    
    #In many cases, we only need to know if job is terminated or not
    #def is_terminated()

    @property
    def is_completed(self):
        return self.status == self.S_COMPLETED

    @property
    def is_running(self):
        return self.status == self.S_RUNNING

    @property
    def is_failed(self):
        return self.status == self.S_FAILED

    @property
    def timeout(self):
        return self.status == self.S_TIMEOUT

    @property
    def has_node_failures(self):
        return self.status == self.S_NODEFAIL

    @property
    def unknown_status(self):
        return self.status == self.S_UNKNOWN

    def set_status_exitcode_signal(self, status, exitcode, signal):
        self.status, self.exitcode, self.signal = status, exitcode, signal

    def likely_code_error(self):
        """
        See http://man7.org/linux/man-pages/man7/signal.7.html

        SIGHUP        1       Term    Hangup detected on controlling terminal or death of controlling process
        SIGINT        2       Term    Interrupt from keyboard
        SIGQUIT       3       Core    Quit from keyboard
        SIGILL        4       Core    Illegal Instruction
        SIGABRT       6       Core    Abort signal from abort(3)
        SIGFPE        8       Core    Floating point exception
        SIGKILL       9       Term    Kill signal
        SIGSEGV      11       Core    Invalid memory reference
        SIGPIPE      13       Term    Broken pipe: write to pipe with no readers
        SIGALRM      14       Term    Timer signal from alarm(2)
        SIGTERM      15       Term    Termination signal
        SIGUSR1   30,10,16    Term    User-defined signal 1
        SIGUSR2   31,12,17    Term    User-defined signal 2
        SIGCHLD   20,17,18    Ign     Child stopped or terminated
        SIGCONT   19,18,25    Cont    Continue if stopped
        SIGSTOP   17,19,23    Stop    Stop process
        SIGTSTP   18,20,24    Stop    Stop typed at terminal
        SIGTTIN   21,21,26    Stop    Terminal input for background process
        SIGTTOU   22,22,27    Stop    Terminal output for background process

        The signals SIGKILL and SIGSTOP cannot be caught, blocked, or ignored.

        Next the signals not in the POSIX.1-1990 standard but described in
        SUSv2 and POSIX.1-2001.

        Signal       Value     Action   Comment
        ────────────────────────────────────────────────────────────────────
        SIGBUS      10,7,10     Core    Bus error (bad memory access)
        SIGPOLL                 Term    Pollable event (Sys V).
                                        Synonym for SIGIO
        SIGPROF     27,27,29    Term    Profiling timer expired
        SIGSYS      12,31,12    Core    Bad argument to routine (SVr4)
        SIGTRAP        5        Core    Trace/breakpoint trap
        SIGURG      16,23,21    Ign     Urgent condition on socket (4.2BSD)
        SIGVTALRM   26,26,28    Term    Virtual alarm clock (4.2BSD)
        SIGXCPU     24,24,30    Core    CPU time limit exceeded (4.2BSD)
        SIGXFSZ     25,25,31    Core    File size limit exceeded (4.2BSD)
        """
        for sig_name in ("SIGFPE",):
            if self.received_signal(sig_name): return sig_name
        return False

    def received_signal(self, sig_name):
        if self.signal is None: return False
        # Get the numeric value from signal and compare it with self.signal
        import signal
        try:
            return self.signal == getattr(signal, sig_name) 
        except AttributeError:
            # invalid sig_name or sig_name not available on this OS.
            return False

    def estimated_start_time(self):
        """Return date with estimated start time. None if it cannot be detected"""
        return None

    def get_info(self, **kwargs):
        return None

    def get_nodes(self, **kwargs):
        return None

    def get_stats(self, **kwargs):
        return None


class ShellJob(QueueJob):
    """Handler for Shell jobs."""
    QTYPE = "shell"


class SlurmJob(QueueJob):
    """Handler for Slurm jobs."""
    QTYPE = "slurm"

    def estimated_start_time(self):
        #squeue  --start -j  116791
        #  JOBID PARTITION     NAME     USER  ST           START_TIME  NODES NODELIST(REASON)
        # 116791      defq gs6q2wop username  PD  2014-11-04T09:27:15     16 (QOSResourceLimit)
        cmd = "squeue" "--start", "--job %d" % self.qid
        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        out, err = process.communicate()

        if process.returncode != 0: 
            logger.critical(err)
            return None

        lines = out.splitlines()
        if len(lines) <= 2: return None

        from datetime import datetime
        for line in lines:
            tokens = line.split()
            if int(tokens[0]) == self.qid:
                date_string = tokens[5]
                if date_string == "N/A": return None
                return datetime.strptime(date_string, "%Y-%m-%dT%H:%M:%S")

        return None

    def get_info(self, **kwargs):
        # See https://computing.llnl.gov/linux/slurm/sacct.html
        #If SLURM job ids are reset, some job numbers will
        #probably appear more than once refering to different jobs.
        #Without this option only the most recent jobs will be displayed.

        #state Displays the job status, or state.
        #Output can be RUNNING, RESIZING, SUSPENDED, COMPLETED, CANCELLED, FAILED, TIMEOUT, 
        #PREEMPTED or NODE_FAIL. If more information is available on the job state than will fit 
        #into the current field width (for example, the uid that CANCELLED a job) the state will be followed by a "+". 

        #gmatteo@master2:~
        #sacct --job 112367 --format=jobid,exitcode,state --allocations --parsable2
        #JobID|ExitCode|State
        #112367|0:0|RUNNING
        #scontrol show job 800197 --oneliner

        # For more info
        #login1$ scontrol show job 1676354

        #cmd = "sacct --job %i --format=jobid,exitcode,state --allocations --parsable2" % self.qid
        cmd = "scontrol show job %i --oneliner" % self.qid
        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        out, err = process.communicate()

        if process.returncode != 0:
            logger.critical(err)
            return None

        tokens = out.splitlines()
        info = AttrDict()
        for line in tokens:
            #print(line)
            k, v = line.split("=")
            info[k] = v
            #print(info)

        qid = int(info.JobId)
        assert qid == self.qid
        exitcode = info.ExitCode
        status = info.JobState

        if ":" in exitcode:
            exitcode, signal = map(int, exitcode.split(":"))
        else:
            exitcode, signal = int(exitcode), None

        i = status.find("+")
        if i != -1: status = status[:i]

        self.set_status_exitcode_signal(JobStatus.from_string(status), exitcode, signal)
        return AttrDict(exitcode=exitcode, signal=signal, status=status)

    def get_stats(self, **kwargs):
        cmd = "sacct --long --job %s --parsable2" % self.qid
        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        out, err = process.communicate()

        if process.returncode != 0: 
            logger.critical(err)
            return {}

        lines = out.splitlines()
        keys = lines[0].strip().split("|")
        values = lines[1].strip().split("|")
        #print("lines0", lines[0])
        return dict(zip(keys, values))


class PbsProJob(QueueJob):
    """
    Handler for PbsPro Jobs.

    See also https://github.com/plediii/pbs_util for a similar project.
    """
    QTYPE = "pbspro"
    # Mapping PrbPro --> Slurm. From `man qstat`
    #
    # S  The job’s state:
    #      B  Array job has at least one subjob running.
    #      E  Job is exiting after having run.
    #      F  Job is finished.
    #      H  Job is held.
    #      M  Job was moved to another server.
    #      Q  Job is queued.
    #      R  Job is running.
    #      S  Job is suspended.
    #      T  Job is being moved to new location.
    #      U  Cycle-harvesting job is suspended due to keyboard activity.
    #      W  Job is waiting for its submitter-assigned start time to be reached.
    #      X  Subjob has completed execution or has been deleted.

    PBSSTAT_TO_SLURM = defaultdict(lambda x: QueueJob.S_UNKNOWN, [
        ("E", QueueJob.S_FAILED),
        ("F", QueueJob.S_COMPLETED),
        ("Q", QueueJob.S_PENDING),
        ("R", QueueJob.S_RUNNING),
        ("S", QueueJob.S_SUSPENDED),
    ])

    def estimated_start_time(self):
        # qstat -T - Shows the estimated start time for all jobs in the queue. 
        #                                                                           Est
        #                                                            Req'd  Req'd   Start
        #Job ID          Username Queue    Jobname    SessID NDS TSK Memory Time  S Time
        #--------------- -------- -------- ---------- ------ --- --- ------ ----- - -----
        #5669001.frontal username large    gs.Pt         --   96  96    --  03:00 Q    --
        cmd = "qstat %s -T" % self.qid
        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        out, err = process.communicate()

        if process.returncode != 0: 
            logger.critical(err)
            return None

        line = out.splitlines()[-1]
        sdate = line.split()[-1]
        if sdate in ("--", "?"): 
            return None

        # TODO One should convert to datetime
        return sdate

    def get_info(self, **kwargs):

        # See also qstat -f 
        #http://sc.tamu.edu/help/origins/batch.shtml#qstat

        #$> qstat 5666289
        #frontal1:
        #                                                            Req'd  Req'd   Elap
        #Job ID          Username Queue    Jobname    SessID NDS TSK Memory Time  S Time
        #--------------- -------- -------- ---------- ------ --- --- ------ ----- - -----
        #5666289.frontal username main_ivy MorfeoTChk  57546   1   4    --  08:00 R 00:17

        cmd = "qstat %d" % self.qid
        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        out, err = process.communicate()

        if process.returncode != 0:
            # qstat: 5904257.frontal1 Job has finished, use -x or -H to obtain historical job information\n
            cmd = "qstat %d -x" % self.qid
            process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
            out, err = process.communicate()

            if process.returncode != 0:
                logger.critical(out)
                logger.critical(err)
                return None

        # Here I don't know what's happeing but I get an output that differs from the one obtained in the terminal.
        # Job id            Name             User              Time Use S Queue
        # ----------------  ---------------- ----------------  -------- - -----
        # 5905011.frontal1  t0               gmatteo           01:37:08 F main_wes  
        #print(out)

        line = out.splitlines()[-1]
        #print(line.split())
        status = self.PBSSTAT_TO_SLURM[line.split()[4]]

        # Exit code and signal are not available.
        # Once could use tracejob....
        # See also http://docs.adaptivecomputing.com/torque/3-0-5/a.gprologueepilogue.php
        self.set_status_exitcode_signal(status, None, None)


#################################
# Unsupported resource managers #
#################################

class TorqueJob(QueueJob):
    """Not supported"""
    QTYPE = "torque"


class SgeJob(QueueJob):
    """Not supported"""
    QTYPE = "sge"


class MoabJob(QueueJob):
    """Not supported"""
    QTYPE = "moab"
