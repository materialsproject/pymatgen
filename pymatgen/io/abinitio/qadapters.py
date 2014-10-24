# coding: utf-8
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
from __future__ import print_function, division, unicode_literals

import os
import abc
import string
import copy
import getpass
import warnings
import six

from subprocess import Popen, PIPE
from monty.string import is_string, boxed
from .launcher import ScriptEditor

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



#class ClusterPartition(object):
#    """
#    This object collects information of a partition
#    # Based on https://computing.llnl.gov/linux/slurm/sinfo.html
#PartitionName=debug TotalNodes=5 TotalCPUs=40 RootOnly=NO
#   Default=YES Shared=FORCE:4 Priority=1 State=UP
#   MaxTime=00:30:00 Hidden=NO
#   MinNodes=1 MaxNodes=26 DisableRootJobs=NO AllowGroups=ALL
#   Nodes=adev[1-5] NodeIndices=0-4
#    """
#    #def from_dict(cls, d):
#
#    def __init__(self, name, num_nodes, timelimit, min_nodes, max_nodes)
#
#        self.name = name
#        #datetime.datetime.strptime("1:00:00", "%H:%M:%S")
#        #self.timelimit = timelimit
#        self.num_sockets = 
#        self.num_nodes = int(num_nodes)
#        self.min_nodes = int(min_nodes)
#        self.max_nodes = int(max_nodes)



def qadapter_class(qtype):
    """Return the concrete `Adapter` class from a string."""
    return {"shell": ShellAdapter,
            "slurm": SlurmAdapter,
            "pbs": PbsAdapter,
            "sge": SGEAdapter,
            "moab": MOABAdapter,
            }[qtype.lower()]


class QueueAdapterError(Exception):
    """Error class for exceptions raise by QueueAdapter."""


class AbstractQueueAdapter(six.with_metaclass(abc.ABCMeta, object)):
    """
    The QueueAdapter is responsible for all interactions with a specific
    queue management system. This includes handling all details of queue
    script format as well as queue submission and management.

    This is the Abstract base class defining the methods that 
    must be implemented by the concrete classes.
    A user should extend this class with implementations that work on
    specific queue systems.
    """
    Error = QueueAdapterError

    # the limits for certain parameters set on the cluster. currently hard coded, should be read at init
    # the increase functions will not increase beyond thise limits
    LIMITS = []

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
        self._verbatim = []

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
                    err_msg += "Unsupported QUEUE parameter name %s\n" % param
                    err_msg += "Supported are: \n"
                    for param_sup in self.supported_qparams:
                        err_msg += "    %s \n" % param_sup
            if err_msg:
                raise ValueError(err_msg)

    def __str__(self):
        lines = [self.__class__.__name__]
        app = lines.append

        if self.has_omp: app(str(self.omp_env))

        return "\n".join(lines)

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
        return hasattr(self, "omp_env") and bool(getattr(self, "omp_env"))

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

    @abc.abstractmethod
    def set_omp_ncpus(self, omp_ncpus):
        """Set the number of OpenMP threads."""

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

    @abc.abstractmethod
    def cancel(self, job_id):
        """
        Cancel the job. 

        Args:
            job_id:
                (in) Job identifier.

        Returns:
            Exit status.
        """

    def add_verbatim(self, lines):
        """
        Add a list of lines or just a string to the header.
        No programmatic interface to change these options is provided
        """
        if is_string(lines): lines = [lines]
        self._verbatim.extend(lines)

    def _make_qheader(self, job_name, qout_path, qerr_path):
        """Return a string with the options that are passed to the resource manager."""
        qtemplate = QScriptTemplate(self.QTEMPLATE)

        # set substitution dict for replacements into the template and clean null values
        subs_dict = {k: v for k, v in self.qparams.items() if v is not None}

        # Set job_name and the names for the stderr and stdout of the 
        # queue manager (note the use of the extensions .qout and .qerr
        # so that we can easily locate this file.
        subs_dict['job_name'] = job_name.replace('/', '_')
        subs_dict['_qout_path'] = qout_path
        subs_dict['_qerr_path'] = qerr_path

        # might contain unused parameters as leftover $$.
        unclean_template = qtemplate.safe_substitute(subs_dict)  

        # Remove lines with leftover $$.
        clean_template = []
        for line in unclean_template.split('\n'):
            if '$$' not in line:
                clean_template.append(line)

        # Add verbatim lines
        if self._verbatim:
            clean_template.extend(self._verbatim)

        return '\n'.join(clean_template)

    def get_script_str(self, job_name, launch_dir, executable, qout_path, qerr_path, stdin=None, stdout=None, stderr=None):
        """
        Returns a (multi-line) String representing the queue script, e.g. PBS script.
        Uses the template_file along with internal parameters to create the script.

        Args:
            job_name:
                Name of the job.
            launch_dir: 
                (str) The directory the job will be launched in.
            qout_path
                Path of the Queue manager output file.
            qerr_path:
                Path of the Queue manager error file.
        """
        # PBS does not accept job_names longer than 15 chars.
        if len(job_name) > 14 and isinstance(self, PbsAdapter):
            job_name = job_name[:14]

        # Construct the header for the Queue Manager.
        qheader = self._make_qheader(job_name, qout_path, qerr_path)

        # Add the bash section.
        se = ScriptEditor()

        if self.setup:
            se.add_comment("Setup section")
            se.add_lines(self.setup)
            se.add_emptyline()

        if self.modules:
            se.add_comment("Load Modules")
            se.add_line("module purge")
            se.load_modules(self.modules)
            se.add_emptyline()

        if self.has_omp:
            se.add_comment("OpenMp Environment")
            se.declare_vars(self.omp_env)
            se.add_emptyline()

        if self.shell_env:
            se.add_comment("Shell Environment")
            se.declare_vars(self.shell_env)
            se.add_emptyline()

        # Cd to launch_dir
        se.add_line("cd " + os.path.abspath(launch_dir))

        if self.pre_run:
            se.add_comment("Commands before execution")
            se.add_lines(self.pre_run)
            se.add_emptyline()

        # Construct the string to run the executable with MPI and mpi_ncpus.
        mpi_ncpus = self.mpi_ncpus

        line = self.mpi_runner.string_to_run(executable, mpi_ncpus, stdin=stdin, stdout=stdout, stderr=stderr)
        se.add_line(line)

        if self.post_run:
            se.add_emptyline()
            se.add_comment("Commands after execution")
            se.add_lines(self.post_run)

        shell_text = se.get_script_str()

        return qheader + shell_text + "\n"

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

    #some method to fix problems

    @abc.abstractmethod
    def exclude_nodes(self, nodes):
        """
        Method to exclude nodes in the calculation
        """

    @abc.abstractmethod
    def increase_mem(self, factor):
        """
        Method to increase the amount of memory asked for, by factor.
        """

    @abc.abstractmethod
    def increase_time(self, factor):
        """
        Method to increase the available wall time asked for, by factor.
        """

    @abc.abstractmethod
    def increase_cpus(self, factor):
        """
        Method to increase the number of cpus asked for.
        """

    # @abc.abstractmethod
    def increase_resources(self):
        """
        Method to generally increase resources.
        """
        return False


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

    def set_omp_ncpus(self, omp_ncpus):
        """Set the number of OpenMP threads."""
        self.omp_env["OMP_NUM_THREADS"] = omp_ncpus

    def set_mem_per_cpu(self, mem_mb):
        """mem_per_cpu is not available in ShellAdapter."""

    def cancel(self, job_id):
        return os.system("kill -9 %d" % job_id)

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
        return 1

    def exclude_nodes(self, nodes):
        return False

    def increase_mem(self, factor):
        return False

    def increase_time(self, factor):
        return False

    def increase_cpus(self, factor):
        return False


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
#SBATCH	--exclude=$${exclude_nodes}
#SBATCH --mem=$${mem}
#SBATCH --mem-per-cpu=$${mem_per_cpu}
#SBATCH --mail-user=$${mail_user}
#SBATCH --mail-type=$${mail_type}
#SBATCH --constraint=$${constraint}
#SBATCH --gres=$${gres}
#SBATCH --requeue=$${requeue}
#SBATCH --nodelist=$${nodelist}
#SBATCH --propagate=$${propagate}

#SBATCH --output=$${_qout_path}
#SBATCH --error=$${_qerr_path}
"""

    LIMITS = {'max_total_tasks': 544, 'max_cpus_per_node': 16, 'mem': 6400000, 'mem_per_cpu': 64000, 'time': 2880}

    @property
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""
        return self.qparams.get("ntasks", 1)

    def set_mpi_ncpus(self, mpi_ncpus):
        """Set the number of CPUs used for MPI."""
        self.qparams["ntasks"] = mpi_ncpus

    def set_omp_ncpus(self, omp_ncpus):
        """Set the number of OpenMP threads."""
        self.omp_env["OMP_NUM_THREADS"] = omp_ncpus
        warnings.warn("set_omp_ncpus not availabe for %s" % self.__class__.__name__)

    def set_mem_per_cpu(self, mem_mb):
        """Set the memory per CPU in Megabytes"""
        self.qparams["mem_per_cpu"] = int(mem_mb)
        # Remove mem if it's defined.
        self.qparams.pop("mem", None)

    def cancel(self, job_id):
        return os.system("scancel %d" % job_id)

    def submit_to_queue(self, script_file, submit_err_file="sbatch.err"):

        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        submit_err_file = os.path.join(os.path.dirname(script_file), submit_err_file)

        # submit the job
        try:
            cmd = ['sbatch', script_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            # write the err output to file, a error parser may read it and a fixer may know what to do ...

            with open(submit_err_file, mode='w') as f:
                f.write('sbatch submit process stderr:')
                f.write(str(process.stderr.read()))
                f.write('qparams:')
                f.write(str(self.qparams))

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
                           "The error response reads: {c}".format(c=process.stderr.read()))
                raise self.Error(err_msg)

        except Exception as details:
            msg = 'Error while submitting job:\n' + str(details)
            logger.critical(msg)
            with open(submit_err_file, mode='a') as f:
                f.write(msg)

            try:
                print('sometimes we land here, no idea what is happening ... Michiel')
                print("details:\n", details, "cmd\n", cmd, "\nprocess.returcode:", process.returncode)
            except:
                pass

            # random error, e.g. no qsub on machine!
            raise self.Error('Running sbatch caused an error...')

    def exclude_nodes(self, nodes):
        try:
            if 'exclude_nodes' not in self.qparams.keys():
                self.qparams.update({'exclude_nodes': 'node'+nodes[0]})
                print('excluded node %s' % nodes[0])
            for node in nodes[1:]:
                self.qparams['exclude_nodes'] += ',node'+node
                print('excluded node %s' % node)
            #print(self.qparams)
            return True
        except (KeyError, IndexError):
            return False

    def increase_cpus(self, factor=1.5):
        logger.info('increasing cpus')
        try:
            if self.qparams['ntasks'] > 1:
                # mpi parallel
                n = int(self.qparams['ntasks'] * factor)
                if n < self.limits['max_total_tasks']:
                    self.qparams['ntasks'] = n
                    logger.info('increased ntasks to %s' % n)
                    return True
                else:
                    raise QueueAdapterError
            elif self.qparams['ntasks'] == 1 and self.qparams['cpus_per_task'] > 1:
                # open mp parallel
                n = int(self.qparams['cpus_per_task'] * factor)
                if n < self.limits['max_cpus_per_node']:
                    self.qparams['cpus_per_task'] = n
                    return True
                else:
                    raise QueueAdapterError
            else:
                raise QueueAdapterError
        except (KeyError, QueueAdapterError):
            return False

    def increase_mem(self, factor=1.5):
        logger.info('increasing memory')
        try:
            if 'mem' in self.qparams.keys():
                n = int(self.qparams['mem'] * factor)
                if n < self.limits['mem']:
                    self.qparams['mem'] = n
                    logger.info('increased mem to %s' % n)
                    return True
                else:
                    raise QueueAdapterError
            elif 'mem_per_cpu' in self.qparams.keys():
                n = int(self.qparams['mem_per_cpu'] * factor)
                if n < self.limits['mem_per_cpu']:
                    self.qparams['mem'] = n
                    logger.info('increased mem_per_cpu to %s' % n)
                    return True
                else:
                    raise QueueAdapterError
            else:
                raise QueueAdapterError
        except (KeyError, IndexError, QueueAdapterError):
            return False

    def increase_time(self, factor=1.5):
        logger.info('increasing time')
        days, hours, minutes = 0, 0, 0
        try:
            # a slurm time parser ;-) forgetting about seconds
            # feel free to pull this out and mak time in minutes always
            if '-' not in self.qparams['time']:
                # "minutes",
                # "minutes:seconds",
                # "hours:minutes:seconds",
                if ':' not in self.qparams['time']:
                    minutes = int(float(self.qparams['time']))
                elif self.qparams['time'].count(':') == 1:
                    minutes = int(float(self.qparams['time'].split(':')[0]))
                else:
                    minutes = int(float(self.qparams['time'].split(':')[1]))
                    hours = int(float(self.qparams['time'].split(':')[0]))
            else:
                # "days-hours",
                # "days-hours:minutes",
                # "days-hours:minutes:seconds".
                days = int(float(self.qparams['time'].split('-')[0]))
                hours = int(float(self.qparams['time'].split('-')[1].split(':')[0]))
                try:
                    minutes = int(float(self.qparams['time'].split('-')[1].split(':')[1]))
                except IndexError:
                    pass
            time = (days * 24 + hours) * 60 + minutes
            time *= factor
            if time < self.limits['time']:
                self.qparams['time'] = time
                logger.info('increased time to %s' % time)
                return True
            else:
                raise QueueAdapterError

        except (KeyError, QueueAdapterError):
            return False

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
#PBS -N $${job_name}
#PBS -l walltime=$${walltime}
#PBS -q $${queue}
#PBS -l model=$${model}
#PBS -l place=$${place}
#PBS -W group_list=$${group_list}
#PBS -l select=$${select}:ncpus=1:vmem=$${vmem}mb:mpiprocs=1:ompthreads=$${ompthreads}
#PBS -l pvmem=$${pvmem}mb
#PBS -r y
####PBS -l mppwidth=$${mppwidth}
####PBS -l nodes=$${nodes}:ppn=$${ppn}  # OLD SYNTAX
#PBS -o $${_qout_path}
#PBS -e $${_qerr_path}
"""
    LIMITS = {'max_total_tasks': 3888, 'time': 48, 'max_select': 300, 'mem': 16000}

    @property
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""
        return self.qparams.get("select", 1)
                                                    
    def set_mpi_ncpus(self, mpi_ncpus):
        """Set the number of CPUs used for MPI."""
        self.qparams["select"] = mpi_ncpus

    def set_omp_ncpus(self, omp_ncpus):
        """Set the number of OpenMP threads."""
        self.omp_env["OMP_NUM_THREADS"] = omp_ncpus
        self.qparams["ompthreads"] = omp_ncpus

    def set_mem_per_cpu(self, mem_mb):
        """Set the memory per CPU in Megabytes"""
        self.qparams["pvmem"] = mem_mb
        self.qparams["vmem"] = mem_mb

    def cancel(self, job_id):
        return os.system("qdel %d" % job_id)

    def submit_to_queue(self, script_file):

        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        # submit the job
        try:
            cmd = ['qsub', script_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            process.wait()

            # grab the return code. PBS returns 0 if the job was successful
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
                raise self.Error(msg)

        except Exception as exc:
            # random error, e.g. no qsub on machine!
            raise self.Error("Running qsub caused an error...\n%s" % str(exc))

    def get_njobs_in_queue(self, username=None):
        # Initialize username
        return 0

        if username is None:
            username = getpass.getuser()

        # run qstat
        try:
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
        except:
            # there's a problem talking to qstat server?
            print(process[1].split('\n'))
            err_msg = ('Error trying to get the number of jobs in the queue using qstat service\n' +
                       'The error response reads: {}'.format(process[2]))
            logger.critical(boxed(err_msg))
            return None

    # no need to raise an error, if False is returned the fixer may try something else, we don't need to kill the
    # scheduler just yet

    def do(self):
        return 'this is not FORTRAN'

    def exclude_nodes(self, nodes):
        logger.warning('exluding nodes, not implemented yet pbs')
        return False

    def increase_mem(self, factor=1):
        base_increase = 2000
        new_mem = self.qparams["pvmem"] + factor*base_increase
        if new_mem < self.limits['mem']:
            self.set_mem_per_cpu(new_mem)
            return True
        else:
            logger.warning('could not increase mem further')
            return False

    def increase_time(self, factor=1.5):
        days, hours, minutes = 0, 0, 0
        try:
            # a pbe time parser [HH:MM]:SS
            # feel free to pull this out and mak time in minutes always
            n = str(self.qparams['time']).count(':')
            if n == 0:
                hours = int(float(self.qparams['time']))
            elif n > 1:
                hours = int(float(self.qparams['time'].split(':')[0]))
                minutes = int(float(self.qparams['time'].split(':')[1]))
            time = hours * 60 + minutes
            time *= factor
            if time < self.limits['time']:
                self.qparams['time'] = str(int(time / 60)) + ':' + str(int(time - 60 * int(time / 60))) + ':00'
                logger.info('increased time to %s minutes' % time)
                return True
            else:
                raise QueueAdapterError
        except (KeyError, QueueAdapterError):
            return False

    def increase_cpus(self, factor):
        base_increase = 12
        new_cpus = self.qparams['select'] + factor * base_increase
        if new_cpus < self.limits['max_select']:
            self.qparams['select'] = new_cpus
            return True
        else:
            logger.warning('increasing cpus reached the limit')
            return False

    # moved to the level of the manager:
    #def increase_resources(self):
    #    """
    #    Method to generally increase resources. On typical large machines we only increas cpu's since we use all
    #    mem per cpu per core
    #    """
    #    if self.increase_cpus(1):
    #        return True
    #    else:
    #        return False


class SGEAdapter(AbstractQueueAdapter):
    """
    Adapter for Sun Grid Engine (SGE) task submission software.
    """
    QTYPE = "sge"

    QTEMPLATE = """\
#!/bin/bash

#$ -A $${account}
#$ -N $${job_name}
#$ -l h rt=$${walltime}
#$ -pe $${queue} $${ncpus}
#$ -cwd
#$ -j y
#$ -m n
#$ -e $${_qerr_path}
#$ -o $${_qout_path}
#$ -S /bin/bash
"""
    @property
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""
        return self.qparams.get("ncpus", 1) 
                                                    
    def set_mpi_ncpus(self, mpi_ncpus):
        """Set the number of CPUs used for MPI."""
        self.qparams["ncpus"] = mpi_ncpus

    def set_omp_ncpus(self, omp_ncpus):
        """Set the number of OpenMP threads."""
        self.omp_env["OMP_NUM_THREADS"] = omp_ncpus
        warnings.warn("set_omp_ncpus not availabe for %s" % self.__class__.__name__)

    def set_mem_per_cpu(self, mem_mb):
        """Set the memory per CPU in Megabytes"""
        raise NotImplementedError("")
        #self.qparams["mem_per_cpu"] = mem_mb
        ## Remove mem if it's defined.
        #self.qparams.pop("mem", None)

    def cancel(self, job_id):
        return os.system("qdel %d" % job_id)

    def submit_to_queue(self, script_file):

        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        # submit the job
        try:
            cmd = ['qsub', script_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            process.wait()

            # grab the returncode. SGE returns 0 if the job was successful
            if process.returncode == 0:
                try:
                    # output should of the form 
                    # Your job 1659048 ("NAME_OF_JOB") has been submitted 
                    queue_id = int(process.stdout.read().split(' ')[2])
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
                raise self.Error(msg)

        except:
            # random error, e.g. no qsub on machine!
            raise self.Error("Running qsub caused an error...")

    def get_njobs_in_queue(self, username=None):
        # Initialize username
        if username is None:
            username = getpass.getuser()

        # run qstat
        qstat = Command(['qstat', '-u', username])
        process = qstat.run(timeout=5)

        # parse the result
        if process[0] == 0:
            # lines should contain username
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

    def exclude_nodes(self, nodes):
        """
        Method to exclude nodes in the calculation
        """
        raise NotImplementedError("exclude_nodes")
                                                                         
    def increase_mem(self, factor):
        """
        Method to increase the amount of memory asked for, by factor.
        """
        raise NotImplementedError("increase_mem")
                                                                         
    def increase_time(self, factor):
        """
        Method to increase the available wall time asked for, by factor.
        """
        raise NotImplementedError("increase_time")

    def increase_cpus(self, factor):
        raise NotImplementedError("increase_cpus")


class MOABAdapter(AbstractQueueAdapter):
# https://computing.llnl.gov/tutorials/moab/
    QTYPE = "moab"

    QTEMPLATE = """\
#!/bin/bash

#MSUB -a $${eligible_date}
#MSUB -A $${account}
#MSUB -c $${checkpoint_interval}
#MSUB -l feature=$${feature}
#MSUB -l gres=$${gres}
#MSUB -l nodes=$${nodes}
#MSUB -l partition=$${partition}
#MSUB -l procs=$${procs}
#MSUB -l ttc=$${ttc}
#MSUB -l walltime=$${walltime}
#MSUB -l $${resources}
#MSUB -p $${priority}
#MSUB -q $${queue}
#MSUB -S $${shell}
#MSUB -N $${job_name}
#MSUB -v $${variable_list}

#MSUB -o $${_qout_path}
#MSUB -e $${_qerr_path}

"""

    @property
    def mpi_ncpus(self):
        """Number of CPUs used for MPI."""
        return self.qparams.get("procs", 1)

    def set_mpi_ncpus(self, mpi_ncpus):
        """Set the number of CPUs used for MPI."""
        self.qparams["procs"] = mpi_ncpus

    def set_omp_ncpus(self, omp_ncpus):
        """Set the number of OpenMP threads."""
        self.omp_env["OMP_NUM_THREADS"] = omp_ncpus

    def cancel(self, job_id):
        return os.system("canceljob %d" % job_id)

    def submit_to_queue(self, script_file, submit_err_file="sbatch.err"):

        if not os.path.exists(script_file):
            raise self.Error('Cannot find script file located at: {}'.format(script_file))

        submit_err_file = os.path.join(os.path.dirname(script_file), submit_err_file)

        # submit the job
        try:
            cmd = ['msub', script_file]
            process = Popen(cmd, stdout=PIPE, stderr=PIPE)
            # write the err output to file, a error parser may read it and a fixer may know what to do ...

            with open(submit_err_file, mode='w') as f:
                f.write('msub submit process stderr:')
                f.write(str(process.stderr.read()))
                f.write('qparams:')
                f.write(str(self.qparams))

            process.wait()

            # grab the returncode. MOAB returns 0 if the job was successful
            if process.returncode == 0:
                try:
                    # output should be the queue_id
                    queue_id = int(process.stdout.read().split()[0])
                    logger.info('Job submission was successful and queue_id is {}'.format(queue_id))
                except:
                    # probably error parsing job code
                    queue_id = None
                    logger.warning('Could not parse job id following msub...')

                finally:
                    return process, queue_id

            else:
                # some qsub error, e.g. maybe wrong queue specified, don't have permission to submit, etc...
                err_msg = ("Error in job submission with MOAB file {f} and cmd {c}\n".format(f=script_file, c=cmd) + 
                           "The error response reads: {c}".format(c=process.stderr.read()))
                raise self.Error(err_msg)

        except Exception as details:
            msg = 'Error while submitting job:\n' + str(details)
            logger.critical(msg)
            with open(submit_err_file, mode='a') as f:
                f.write(msg)

            try:
                print('sometimes we land here, no idea what is happening ... Michiel')
                print("details:\n", details, "cmd\n", cmd, "\nprocess.returcode:", process.returncode)
            except:
                pass

            # random error, e.g. no qsub on machine!
            raise self.Error('Running msub caused an error...')

    def get_njobs_in_queue(self, username=None):
        if username is None:
            username = getpass.getuser()

        cmd = ['showq', '-s -u', username]
        process = Popen(cmd, shell=False, stdout=PIPE)
        process.wait()

        # parse the result
        if process.returncode == 0:
            # lines should have this form:
            ## 
            ## active jobs: N  eligible jobs: M  blocked jobs: P
            ##
            ## Total job:  1
            ##
            # Split the output string and return the last element.

            outs = process.stdout.readlines()
            njobs = int(outs.split()[-1])
            logger.info('The number of jobs currently in the queue is: {}'.format(njobs))
            return njobs

        # there's a problem talking to squeue server?
        err_msg = ('Error trying to get the number of jobs in the queue using showq service' + 
                   'The error response reads: {}'.format(process.stderr.read()))
        logger.critical(err_msg)

        return None
    
    def exclude_nodes(self, nodes):
        raise NotImplementedError("exclude_nodes")
                                                                         
    def increase_mem(self, factor):
        raise NotImplementedError("increase_mem")
                                                                         
    def increase_time(self, factor):
        raise NotImplementedError("increase_time")

    def increase_cpus(self, factor):
        raise NotImplementedError("increase_cpus")
        
    def set_mem_per_cpu(self, factor):
        raise NotImplementedError("set_mem_per_cpu")

class QScriptTemplate(string.Template):
    delimiter = '$$'
