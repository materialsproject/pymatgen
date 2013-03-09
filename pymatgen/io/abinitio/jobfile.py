from os import makedirs
from os.path import basename, dirname, join, abspath, exists, realpath
from copy import deepcopy

__all__ = ['JobFile', 'PBSJobFile', 'SGEJobFile', 'SlurmJobFile']

# =========================================================================== #

__author__ = "Gabriel Antonius"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"

class JobFile(object):
    """
    The job file is organized as follow::

        >>> 
        #!/bin/csh                           # 1) Shell specification.
                                             #
        #PBS -N jobname                      # 2) Submission commands.
        #PBS -l walltime=48:00:00            #    Depends on the subclass used.
        #PBS -l nodes=4:ppn=12               #
                                             #
        set MPIRUN="mpiexec"                 # 3) Declarations.
        set EXECUTABLE=/path/to/executable   #    These are also properties
        set INPUT=calculation.in             #    of the jobfile object.
        set LOG=calculation.log              #
                                             #
        module load intel-compilers          # 4) Modules.
        module load MPI/Intel/mvapich2       #
                                             #
        cd ${PBS_O_WORKDIR}                  # 5) Lines before execution.
        limit stacksize unlimited            #
                                             #
        $MPIRUN $EXECUTABLE < $INPUT > $LOG  # 6) Execution line.
                                             #
        echo "Job done!"                     # 7) Lines after execution.
        date                                 #

    Base class properties:

        shell:
            The shell binary to be used.
            E.g. '/bin/csh'.
            Default is '/bin/bash'.
        mpirun:
            The mpi runner.
            E.g. 'mpiexec -npernode 6'.
            Default is none.
        executable:
            The binary to be executed.
            Default is abinit.
        bindir:
            The directory in which to look for binaries.
            Default is none.
        input:
            The input file to feed in the executable as the standard input.
            Mandatory.
        log:
            The file into which the standard output is redirected.
            Default is 'log'.
        stderr:
            The file into which the standard error is redirected.
            Default is 'stderr'.
        modules:
            The modules which will be loaded with 'module load'.
            Default is none.
        lines_before:
            Lines before the main execution.
            Default is none.
        lines_after:
            Lines after the main execution.
            Default is none.
        other_lines:
            Other lines your job submission script would like to have.
            Must be preceded by the approbriate tag (#!, #PBS).
            Default is none.
        submission_command:
            The command which should be used to launch the job.
            E.g. 'qsub', 'bqsub', 'sbatch'.
            Default depends on the job type.
    """

    _executable = 'abinit'

    def __init__(self, name='job.sh', **kwargs):

        # Name
        self.name = name
        self.absdir = realpath(dirname(self.name))
        self.absname = realpath(self.name)

        # Shell
        self.shell = '/bin/bash'

        # Execution lines
        self.input = ''
        self.log = 'log'
        self.stderr = 'stderr'
        self.executable = 'abinit'
        self.bindir = ''
        self.mpirun = ''

        # Modules
        self.modules = list()

        # Other lines
        self.other_lines = list()
        self.lines_before = list()
        self.lines_after = list()

        # Command used to submit the job
        self.submission_command = 'qsub'

        # Set attributes
        for (arg, val) in kwargs.iteritems():
            try:
                getattr(self, 'set_' + arg)(val)
            except:
                pass

    def __str__(self):
        lines = []
        app = lines.append

        # Shell line
        app('#!' + self.shell)
        app('')

        # Submission instructions
        lines.extend(self._get_command_lines())

        # Other submission inscrutions
        lines.extend(self.other_lines)
        app('')

        # Declarations
        for (key, value) in [ ('MPIRUN', self.mpirun),
                              ('EXECUTABLE', self.executable),
                              ('INPUT', self.input),
                              ('LOG', self.log),
                              ('STDERR', self.stderr),
                            ]:
            app(self._declaration_line(key, value))
        app('')

        # Modules
        for module in self.modules:
            app('module load ' + module)
        app('')

        # Lines before execution
        lines.extend(self.lines_before)
        app('')

        # Execution lines
        app('$MPIRUN $EXECUTABLE < $INPUT > $LOG 2> $STDERR')
        app('')

        # Lines after execution
        lines.extend(self.lines_after)
        app('')

        return "\n".join(lines)

    def _get_command_lines(self):
        """Return the lines specifying instructions for job submission."""
        return list()

    def _declaration_line(self, key, val):
        """Return a lines setting a variable."""
        if 'csh' in self.shell:
            declare = 'set '
        elif 'tcsh' in self.shell:
            declare = 'set '
        else:
            declare = ''
        return declare + key + '=' + val + '\n'

    def _set_property(self, name, *args, **kwargs):
        """Set a property through the corresponding 'set_' function."""
        return getattr(self, 'set_' + key)(*args, **kwargs)

    @property
    def dirname(self):
        """The directory containing the file."""
        return dirname(self.name)

    @property
    def path(self):
        """The path of the file."""
        return abspath(self.name)

    @property
    def basename(self):
        """The base name of the file."""
        return basename(self.name)

    @classmethod
    def properties(cls):
        """Return the list of properties with a 'set_' function."""
        funcs = filter(lambda s: s.startswith('set_'), dir(cls))
        return [ f.split('set_', 1)[-1] for f in funcs ]

    def write(self, name=None):
        """Write the file."""
        if name is None:
            name = self.name

        if not exists(self.dirname):
            makedirs(self.dirname)

        with open(name, 'w') as f:
            f.write(str(self))

    @property
    def executable(self):
        return join(self.bindir, self._executable)

    @executable.setter
    def executable(self, executable):
        self._executable = basename(executable)
        if basename(executable) != executable:
            self.set_bindir(dirname(executable))

    def set_shell(self, shell):
        """
        Sets the shell type.  The argument can either be an absolute path,
        or just the shell type e.g. bash, csh, tcsh, in which case
        the executable is assumed to be located in /bin/.
        The shell also determine how a variable is declared.
        """
        if shell == basename(shell):
            self.shell = join('/bin', shell)
        else:
            self.shell = abspath(shell)

    def set_mpirun(self, mpirun):
        """
        Set the mpi runner to execute the program.
        E.g. 'mpiexec -npernode 6', 'mpirun -np 12', ''.
        """
        if not (mpirun.startswith("'") or mpirun.startswith('"')):
            mpirun = '"' + mpirun + '"'
        self.mpirun = mpirun

    def set_bindir(self, bindir):
        """Set the directory for binaries (abinit, mrgscr...)."""
        self.bindir = realpath(bindir)

    def set_executable(self, executable):
        """Set the executable to use."""
        self.executable = executable

    def set_input(self, input):
        """Set the input file for the main executable."""
        self.input = input

    def set_log(self, log):
        """Set the log file to collect standard output of the executable."""
        self.log = log

    def set_stderr(self, stderr):
        """Set the log file to collect standard output of the executable."""
        self.stderr = stderr

    def set_modules(self, *modules):
        """Set one or many modules to be loaded."""
        if not modules:
            modules = list()
        elif '__iter__' in dir(modules[0]):
            modules = tuple(modules[0])

        self.modules = deepcopy(modules)

    def set_lines_before(self, *lines):
        """Set one or many lines to be executed before the main execution."""
        if not lines:
            lines = list()
        elif '__iter__' in dir(lines[0]):
            lines = tuple(lines[0])

        self.lines_before = list()
        for line in lines:
            self.lines_before.append(line.rstrip('\n'))

    def set_lines_after(self, *lines):
        """Set one or many lines to be executed before the main execution."""
        if not lines:
            lines = list()
        elif '__iter__' in dir(lines[0]):
            lines = tuple(lines[0])

        self.lines_after = list()
        for line in lines:
            self.lines_after.append(line.rstrip('\n'))

    def set_submission_command(self, command):
        """
        Sets the command used for job submission,
        e.g. qsub, bqsub, sbatch, ...
        """
        self.submission_command = command

# =========================================================================== #


class PBSJobFile(JobFile):
    """
    Portable Batch System.

    Properties:

        jobname:
            Name of the job.
        runtime:
            Maximum time for the job.
        nodes:
            Number of nodes on which to run the job.
        ppn:
            Number of processors per node.
        memory:
            Memory per node. E.g. '48G'.
        queue:
            The queue to which the job is submitted.
        mail:
            The mail to which a notification will be sent.
        mail_options:
            The conditions under which a mail will be sent.
            E.G. 'abe'.
        submission_command:
            default is 'qsub'.

    See man qsub for more info.
    """
    __doc__ += "\n" + JobFile.__doc__

    _command = "#PBS "

    def __inti__(self, **kwargs):

        kwargs.setdefault(submission_command='qsub')
        JobFile.__init__(self, **kwargs)

    nodes = None
    def set_nodes(self, val):
        self.nodes = val

    ppn = None
    def set_ppn(self, val):
        self.ppn = val

    memory = None
    def set_memory(self, val):
        self.memory = val

    runtime = None
    def set_runtime(self, val):
        """Either set the numer of hours, or a triplet for (hours,min,sec)."""
        if isinstance(val, int):
            val = [val, 0, 0]
        self.runtime = val

    jobname = None
    def set_jobname(self, val):
        self.jobname = val

    queue = None
    def set_queue(self, val):
        self.queue = val

    mail = None
    def set_mail(self, val):
        self.mail = val

    mail_options = None
    def set_mail_options(self, val):
        self.mail_options = val

    def _get_command_lines(self):
        """Return the lines specifying instructions for job submission."""
        lines = list()
        def add(line):
            lines.append(self._command + line + '\n')

        if self.jobname:
            add('-N ' + str(self.jobname))

        if self.runtime:
            add('-l walltime={}:{}:{}'.format(*self.runtime))

        if self.nodes and self.ppn:
            add('-l nodes=' + str(self.nodes) + ':ppn=' + str(self.ppn))

        if self.memory:
            add('-l mem=' + str(self.memory))

        if self.queue:
            add('-q ' + self.queue)

        if self.mail:
            add('-M ' + self.mail)

        if self.mail_options:
            add('-m ' + self.mail_options)

        return lines

# =========================================================================== #


class SGEJobFile(JobFile):
    """
    Sun Grid Engine.

    Properties:

        jobname:
            Name of the job.
        runtime:
            Maximum time for the job.
        nproc:
            Number of processors.
        queue:
            The queue to which the job is submitted.
        environment:
            The parallel environment under which the job is ran.
        memory:
            The requested memory, in M.
        mail:
            The mail to which a notification will be sent.
        mail_options:
            The conditions under which a mail will be sent.
            E.G. 'abe'.
        submission_command:
            default is 'qsub'.

    See man qsub for more info.
    """
    __doc__ += "\n" + JobFile.__doc__

    _command = "#$ "

    def __inti__(self, **kwargs):

        kwargs.setdefault(submission_command='qsub')
        JobFile.__init__(self, **kwargs)

    jobname = None
    def set_jobname(self, val):
        self.jobname = val

    runtime = None
    def set_runtime(self, val):
        """Either set the numer of hours, or a triplet for (hours,min,sec)."""
        if isinstance(val, int):
            val = [val, 0, 0]
        self.runtime = val

    nproc = None
    def set_nproc(self, val):
        self.nproc = val

    queue = None
    def set_queue(self, val):
        self.queue = val

    environment = None
    def set_environment(self, val):
        self.environment = val

    memory = None
    def set_memory(self, val):
        self.memory = val

    mail = None
    def set_mail(self, val):
        self.mail = val

    mail_options = None
    def set_mail_options(self, val):
        self.mail_options = val

    def _get_command_lines(self):
        """Return the lines specifying instructions for job submission."""
        lines = list()
        def add(line):
            lines.append(self._command + line + '\n')

        if self.jobname:
            add('-N ' + str(self.jobname))

        if self.runtime:
            add('-l h_rt={}:{}:{}'.format(*self.runtime))

        if self.environment and self.nproc:
            line = '-pe ' + self.environment + ' ' + str(self.nproc)
            if self.memory:
                line += ' -l mem=' + str(self.memory)
            add(line)

        if self.queue:
            add('-q ' + self.queue)

        if self.mail:
            add('-M ' + self.mail)

        if self.mail_options:
            add('-m ' + self.mail_options)

        return lines

# =========================================================================== #


class SlurmJobFile(JobFile):
    """
    Simple Linux Utility for Resource Management.

    Properties:

        jobname:
            Name of the job.
        time:
            Maximum time for the job.
        ntasks:
            The number of processes.
        cpus_per_task:
            The number of cpus per process.
        mem_per_cpu:
            The memory per cpu.
        partition:
            The partition...
        mail_user:
            The mail to which a notification is sent.
        mail_type:
            The conditions unde which to send a mail.
        submission_command:
            default is 'sbatch'.
    """
    __doc__ += "\n" + JobFile.__doc__

    _command = "#SBATCH "

    def __inti__(self, **kwargs):

        kwargs.setdefault(submission_command='sbatch')
        JobFile.__init__(self, **kwargs)

    jobname = None
    def set_jobname(self, val):
        self.jobname = val

    time = None
    def set_time(self, val):
        """Either set the numer of hours, or a triplet for (hours,min,sec)."""
        if isinstance(val, int):
            val = [val, 0, 0]
        self.time = val

    def set_runtime(self, val): self.set_time(val)

    ntasks = None
    def set_ntasks(self, val):
        self.ntasks = val

    ntasks_per_node = None
    def set_ntasks_per_node(self, val):
        self.ntasks_per_node = val

    cpus_per_task = None
    def set_cpus_per_task(self, val):
        self.cpus_per_task = val

    mem_per_cpu = None
    def set_mem_per_cpu(self, val):
        self.mem_per_cpu = val

    partition = None
    def set_partition(self, val):
        self.partition = val

    mail_user = None
    def set_mail_user(self, val):
        self.mail_user = val

    mail_type = None
    def set_mail_type(self, val):
        self.mail_type = val

    def _get_command_lines(self):
        """Return the lines specifying instructions for job submission."""
        lines = list()
        def add(line):
            lines.append(self._command + line + '\n')

        if self.jobname:
            add('--job-name=' + str(self.jobname))

        if self.time:
            add('--time={}:{}:{}\n'.format(*self.time))

        if self.ntasks:
            add('--ntasks=' + str(self.ntasks))

        if self.partition:
            add('--partition=' + self.partition)

        if self.ntasks_per_node:
            add('--ntasks-per-node=' + str(self.ntasks_per_node))

        if self.cpus_per_task:
            add('--cpus-per-task=' + str(self.cpus_per_task))

        if self.mem_per_cpu:
            add('--mem-per-cpu=' + str(self.mem_per_cpu))

        if self.mail_user:
            add('--mail-user=' + self.mail_user)

        if self.mail_type:
            add('--mail-type=' + self.mail_type)

        return lines
