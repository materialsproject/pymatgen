"""
Part of this code is based on a similar implementation present in FireWorks (https://pypi.python.org/pypi/FireWorks).
Work done by D. Waroquiers, A. Jain, and M. Kocher.

The main difference wrt the Fireworks implementation is that the QueueAdapter
objects provide a programmatic interface for setting important attributes 
such as the number of MPI nodes, the number of OMP threads and the memory requirements.
This programmatic interface is used by the `TaskManager` for optimizing the parameters
of the run before submitting the job (Abinit provides the autoparal option that 
allows one to get a list of parallel configuration and their expected efficiency).
"""
from __future__ import print_function, division

import os
import abc
import string
import copy
import functools
import getpass

from subprocess import Popen, PIPE
from pymatgen.io.abinitio.launcher import ScriptEditor
from pymatgen.util.string_utils import is_string

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "MpiRunner",
    "qadapter_class",
]

class Command(object):
    """
    From https://gist.github.com/kirpit/1306188

    Enables to run subprocess commands in a different thread with TIMEOUT option.

    Based on jcollado's solution:
    http://stackoverflow.com/questions/1191374/subprocess-with-timeout/4825933#4825933
    """
    command = None
    process = None
    status = None
    output, error = '', ''

    def __init__(self, command):
        if is_string(command):
            import shlex
            command = shlex.split(command)

        self.command = command

    def run(self, timeout=None, **kwargs):
        """ Run a command then return: (status, output, error). """

        def target(**kwargs):
            try:
                self.process = Popen(self.command, **kwargs)
                self.output, self.error = self.process.communicate()
                self.status = self.process.returncode
            except:
                import traceback
                self.error = traceback.format_exc()
                self.status = -1

        # default stdout and stderr
        if 'stdout' not in kwargs:
            kwargs['stdout'] = PIPE
        if 'stderr' not in kwargs:
            kwargs['stderr'] = PIPE
        # thread
        import threading
        thread = threading.Thread(target=target, kwargs=kwargs)
        thread.start()
        thread.join(timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()

        return self.status, self.output, self.error


class MpiRunner(object):
    """
    This object provides an abstraction for the mpirunner provided 
    by the different MPI libraries. It's main task is handling the
    different syntax and options supported by the different mpirunners.
    """
    def __init__(self, name, type=None, options=""):
        self.name = name
        self.type = None
        self.options = options

    def string_to_run(self, executable, mpi_ncpus, stdin=None, stdout=None, stderr=None):
        stdin = "< " + stdin if stdin is not None else ""
        stdout = "> " + stdout if stdout is not None else ""
        stderr = "2> " + stderr if stderr is not None else ""

        if self.has_mpirun:

            if self.type is None:
                # TODO: better treatment of mpirun syntax.
                #se.add_line('$MPIRUN -n $MPI_NCPUS $EXECUTABLE < $STDIN > $STDOUT 2> $STDERR')
                num_opt = "-n " + str(mpi_ncpus)
                cmd = " ".join([self.name, num_opt, executable, stdin, stdout, stderr])

            else:
                raise NotImplementedError("type %s is not supported!")

        else:
            #assert mpi_ncpus == 1
            cmd = " ".join([executable, stdin, stdout, stderr])

        return cmd

    @property
    def has_mpirun(self):
        """True if we are running via mpirun, mpiexec ..."""
        return self.name is not None


def qadapter_class(qtype):
    """Return the concrete `Adapter` class from a string."""
    return {"shell": ShellAdapter,
            "slurm": SlurmAdapter,
            "pbs": PbsAdapter,
            }[qtype.lower()]


class QueueAdapterError(Exception):
    """Error class for exceptions raise by QueueAdapter."""


class AbstractQueueAdapter(object):
    """
    The QueueAdapter is responsible for all interactions with a specific
    queue management system. This includes handling all details of queue
    script format as well as queue submission and management.

    This is the Abstract base class defining the methods that 
    must be implemented by the concrete classes.
    A user should extend this class with implementations that work on
    specific queue systems.
    """
    __metaclass__ = abc.ABCMeta

    Error = QueueAdapterError

    def __init__(self, qparams=None, setup=None, modules=None, shell_env=None, omp_env=None, 
                 pre_run=None, post_run=None, mpi_runner=None):
        """
        Args:
            setup:
                String or list of commands to execute during the initial setup.
            modules:
                String or list of modules to load before running the application.
            shell_env:
                Dictionary with the environment variables to export
                before running the application.
            omp_env:
                Dictionary with the OpenMP variables.
            pre_run:
                String or list of commands to execute before launching the calculation.
            post_run:
                String or list of commands to execute once the calculation is completed.
            mpi_runner:
                Path to the MPI runner or `MpiRunner` instance. None if not used
        """
        # Make defensive copies so that we can change the values at runtime.
        self.qparams = qparams.copy() if qparams is not None else {}

        if is_string(setup):
            setup = [setup]
        self.setup = setup[:] if setup is not None else []

        self.omp_env = omp_env.copy() if omp_env is not None else {}

        if is_string(modules):
            modules = [modules]
        self.modules = modules[:] if modules is not None else []

        self.shell_env = shell_env.copy() if shell_env is not None else {}

        self.mpi_runner = mpi_runner
        if not isinstance(mpi_runner, MpiRunner):
            self.mpi_runner = MpiRunner(mpi_runner)

        if is_string(pre_run):
            pre_run = [pre_run]
        self.pre_run = pre_run[:] if pre_run is not None else []

        if is_string(post_run):
            post_run = [post_run]
        self.post_run = post_run[:] if post_run is not None else []

        # Parse the template so that we know the list of supported options.
        cls = self.__class__
        if hasattr(cls, "QTEMPLATE"): 
            # Consistency check.
            err_msg = ""
            for param in self.qparams:
                if param not in self.supported_qparams:
                    err_msg += "Unsupported QUEUE parameter name %s\n"  % param

            if err_msg:
                raise ValueError(err_msg)

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    @property
    def supported_qparams(self):
        """
        Dictionary with the supported parameters that can be passed to the 
        queue manager (obtained by parsing QTEMPLATE).
        """ 
        try:
            return self._supported_qparams

        except AttributeError:
            import re
            self._supported_qparams = re.findall("\$\$\{(\w+)\}", self.QTEMPLATE)
            return self._supported_qparams
    
    @property
    def has_mpirun(self):
        """True if we are using a mpirunner"""
        return bool(self.mpi_runner)

    @property
    def has_omp(self):
        """True if we are using OpenMP threads"""
        return hasattr(self,"omp_env") and bool(getattr(self, "omp_env"))

    @property
    def tot_ncpus(self):
        """Total number of CPUs employed"""
        return self.mpi_ncpus * self.omp_ncpus 

    @property
    def omp_ncpus(self):
        """Number of OpenMP threads."""
        if self.has_omp:
            return self.omp_env["OMP_NUM_THREADS"]
        else:
            return 1

    @abc.abstractproperty
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""

    @abc.abstractmethod
    def set_mpi_ncpus(self, mpi_ncpus):
        """Set the number of CPUs used for MPI."""

    #@abc.abstractproperty
    #def queue_walltime(self):
    #    """Returns the walltime in seconds."""

    #@abc.abstractmethod
    #def set_queue_walltime(self):
    #    """Set the walltime in seconds."""

    #@abc.abstractproperty
    #def mem_per_cpu(self):
    #    """The memory per CPU in Megabytes."""
                                                
    @abc.abstractmethod
    def set_mem_per_cpu(self, mem_mb):
        """Set the memory per CPU in Megabytes"""

    #@property
    #def tot_mem(self):
    #    """Total memory required by the job n Megabytes."""
    #    return self.mem_per_cpu * self.mpi_ncpus

    def _make_qheader(self, job_name, qout_path, qerr_path):
        """Return a string with the options that are passed to the resource manager."""
        qtemplate = QScriptTemplate(self.QTEMPLATE)

        # set substitution dict for replacements into the template and clean null values
        subs_dict = {k: v for k,v in self.qparams.items() if v is not None}  

        # Set job_name and the names for the stderr and stdout of the 
        # queue manager (note the use of the extensions .qout and .qerr
        # so that we can easily locate this file.
        subs_dict['job_name'] = job_name 
        subs_dict['_qout_path'] = qout_path
        subs_dict['_qerr_path'] = qerr_path

        # might contain unused parameters as leftover $$.
        unclean_template = qtemplate.safe_substitute(subs_dict)  

        # Remove lines with leftover $$.
        clean_template = []
        for line in unclean_template.split('\n'):
            if '$$' not in line:
                clean_template.append(line)

        return '\n'.join(clean_template)

    def get_script_str(self, job_name, launch_dir, executable, qout_path, qerr_path, stdin=None, stdout=None, stderr=None):
        """
        Returns a (multi-line) String representing the queue script, e.g. PBS script.
        Uses the template_file along with internal parameters to create the script.

        Args:
            launch_dir: 
                (str) The directory the job will be launched in.
            qout_path
                Path of the Queue manager output file.
            qerr_path:
                Path of the Queue manager error file.
        """
        # Construct the header for the Queue Manager.
        qheader = self._make_qheader(job_name, qout_path, qerr_path)

        # Add the bash section.
        se = ScriptEditor()

        if self.setup:
            se.add_comment("Setup section")
            se.add_lines(self.setup)

        if self.modules:
            se.add_comment("Load Modules")
            se.add_line("module purge")
            se.load_modules(self.modules)

        if self.has_omp:
            se.add_comment("OpenMp Environment")
            se.declare_vars(self.omp_env)

        if self.shell_env:
            se.add_comment("Shell Environment")
            se.declare_vars(self.shell_env)

        # Cd to launch_dir
        se.add_line("cd " + os.path.abspath(launch_dir))

        if self.pre_run:
            se.add_comment("Commands before execution")
            se.add_lines(self.pre_run)

        # Construct the string to run the executable with MPI and mpi_ncpus.
        mpi_ncpus = self.mpi_ncpus

        line = self.mpi_runner.string_to_run(executable, mpi_ncpus, stdin=stdin, stdout=stdout, stderr=stderr)
        se.add_line(line)

        if self.post_run:
            se.add_comment("Commands after execution")
            se.add_lines(self.post_run)

        shell_text = se.get_script_str()

        return qheader + shell_text

    @abc.abstractmethod
    def submit_to_queue(self, script_file):
        """
        Submits the job to the queue, probably using subprocess or shutil

        Args:
            script_file: 
                (str) name of the script file to use (String)
        Returns:
            process, queue_id
        """

    @abc.abstractmethod
    def get_njobs_in_queue(self, username=None):
        """
        returns the number of jobs in the queue, probably using subprocess or shutil to
        call a command like 'qstat'. returns None when the number of jobs cannot be determined.

        Args:
            username: (str) the username of the jobs to count (default is to autodetect)
        """

####################
# Concrete classes #
####################

class ShellAdapter(AbstractQueueAdapter):
    QTYPE = "shell"

    QTEMPLATE = """\
#!/bin/bash

export MPI_NCPUS=$${MPI_NCPUS}

"""

    @property
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""
        return self.qparams.get("MPI_NCPUS", 1)
                                                    
    def set_mpi_ncpus(self, mpi_ncpus):
        """Set the number of CPUs used for MPI."""
        self.qparams["MPI_NCPUS"] = mpi_ncpus

    def set_mem_per_cpu(self, mem_mb):
        """mem_per_cpu is not available in ShellAdapter."""

    def submit_to_queue(self, script_file):

        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        # submit the job
        try:
            process = Popen(("/bin/bash", script_file), stderr=PIPE)
            queue_id = process.pid
            return process, queue_id

        except:
            # random error
            raise self.Error("Random Error ...!")

    def get_njobs_in_queue(self, username=None):
        return None


class SlurmAdapter(AbstractQueueAdapter):
    QTYPE = "slurm"

    QTEMPLATE = """\
#!/bin/bash

#SBATCH --ntasks=$${ntasks}
#SBATCH --ntasks-per-node=$${ntasks_per_node}
#SBATCH --cpus-per-task=$${cpus_per_task}
#SBATCH --time=$${time}
#SBATCH --partition=$${partition}
#SBATCH --account=$${account}
#SBATCH --job-name=$${job_name}
#SBATCH	--nodes=$${nodes}
#SBATCH --mem=$${mem}
#SBATCH --mem-per-cpu=$${mem_per_cpu}
#SBATCH --mail-user=$${mail_user}
#SBATCH --mail-type=$${mail_type}
#SBATCH --constraint=$${constraint}
#SBATCH --gres=$${gres}
#SBATCH --requeue=$${requeue}
#SBATCH --nodelist=$${nodelist}
#SBATCH --propagate=$${propagate}

#SBATCH --output=$${_qerr_path}
#SBATCH --error=$${_qout_path}

"""

    @property
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""
        return self.qparams.get("ntasks", 1)

    def set_mpi_ncpus(self, mpi_ncpus):
        """Set the number of CPUs used for MPI."""
        self.qparams["ntasks"] = mpi_ncpus

    def set_mem_per_cpu(self, mem_mb):
        """Set the memory per CPU in Megabytes"""
        self.qparams["mem_per_cpu"] = mem_mb
        # Remove mem if it's defined.
        self.qparams.pop("mem", None)

    def submit_to_queue(self, script_file):

        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        # submit the job
        try:
            cmd = ['sbatch', script_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            process.wait()

            # grab the returncode. SLURM returns 0 if the job was successful
            if process.returncode == 0:
                try:
                    # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
                    queue_id = int(process.stdout.read().split()[3])
                    logger.info('Job submission was successful and queue_id is {}'.format(queue_id))

                except:
                    # probably error parsing job code
                    queue_id = None
                    logger.warning('Could not parse job id following slurm...')

                finally:
                    return process, queue_id

            else:
                # some qsub error, e.g. maybe wrong queue specified, don't have permission to submit, etc...
                err_msg = ("Error in job submission with SLURM file {f} and cmd {c}\n".format(f=script_file, c=cmd) + 
                           "The error response reads: {}".format(process.stderr.read()))
                raise self.Error(err_msg)

        except:
            # random error, e.g. no qsub on machine!
            raise self.Error('Running sbatch caused an error...')

    def get_njobs_in_queue(self, username=None):
        if username is None:
            username = getpass.getuser()

        cmd = ['squeue', '-o "%u"', '-u', username]
        process = Popen(cmd, shell=False, stdout=PIPE)
        process.wait()

        # parse the result
        if process.returncode == 0:
            # lines should have this form
            # username
            # count lines that include the username in it

            outs = process.stdout.readlines()
            njobs = len([line.split() for line in outs if username in line])
            logger.info('The number of jobs currently in the queue is: {}'.format(njobs))
            return njobs

        # there's a problem talking to squeue server?
        err_msg = ('Error trying to get the number of jobs in the queue using squeue service' + 
                   'The error response reads: {}'.format(process.stderr.read()))
        logger.critical(err_msg)

        return None


class PbsAdapter(AbstractQueueAdapter):
    QTYPE = "pbs"

    QTEMPLATE = """\
#!/bin/bash

#PBS -A $${account}
#PBS -l walltime=$${walltime}
#PBS -q $${queue}
#PBS -l mppwidth=$${mppwidth}
#PBS -l nodes=$${nodes}:ppn=$${ppn}
#PBS -N $${job_name}
#PBS -o $${_qerr_path}
#PBS -e $${_qout_path}

"""
    @property
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""
        return self.qparams.get("nodes", 1) * self.qparams.get("ppn", 1)
                                                    
    def set_mpi_ncpus(self, mpi_ncpus):
        """Set the number of CPUs used for MPI."""
        if "ppn" not in self.qparams: self.qparams["ppn"] = 1

        ppnode = self.qparams.get("ppn")
        self.qparams["nodes"] = mpi_ncpus // ppnode 

    def set_mem_per_cpu(self, mem_mb):
        """Set the memory per CPU in Megabytes"""
        raise NotImplementedError("")
        #self.qparams["mem_per_cpu"] = mem_mb
        ## Remove mem if it's defined.
        #self.qparams.pop("mem", None)

    def submit_to_queue(self, script_file):

        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        # submit the job
        try:
            cmd = ['qsub', scriprocesst_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            process.wait()

            # grab the returncode. PBS returns 0 if the job was successful
            if process.returncode == 0:
                try:
                    # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
                    queue_id = int(process.stdout.read().split('.')[0])
                    logger.info('Job submission was successful and queue_id is {}'.format(queue_id))

                except:
                    # probably error parsing job code
                    logger.warning("Could not parse job id following qsub...")
                    queue_id = None

                finally:
                    return process, queue_id

            else:
                # some qsub error, e.g. maybe wrong queue specified, don't have permission to submit, etc...
                msg = ('Error in job submission with PBS file {f} and cmd {c}\n'.format(f=script_file, c=cmd) + 
                      'The error response reads: {}'.format(process.stderr.read()))

        except:
            # random error, e.g. no qsub on machine!
            raise self.Error("Running qsub caused an error...")


    def get_njobs_in_queue(self, username=None):
        # Initialize username
        if username is None:
            username = getpass.getuser()

        # run qstat
        qstat = Command(['qstat', '-a', '-u', username])
        process = qstat.run(timeout=5)

        # parse the result
        if process[0] == 0:
            # lines should have this form
            # '1339044.sdb          username  queuename    2012-02-29-16-43  20460   --   --    --  00:20 C 00:09'
            # count lines that include the username in it

            # TODO: only count running or queued jobs. or rather, *don't* count jobs that are 'C'.
            outs = process[1].split('\n')
            njobs = len([line.split() for line in outs if username in line])
            logger.info('The number of jobs currently in the queue is: {}'.format(njobs))

            return njobs

        # there's a problem talking to qstat server?
        err_msg = ('Error trying to get the number of jobs in the queue using qstat service\n' + 
                   'The error response reads: {}'.format(process[2]))
        logger.critical(err_msg)

        return None


class QScriptTemplate(string.Template):
    delimiter = '$$'
