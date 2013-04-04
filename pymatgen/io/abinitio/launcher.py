#!/usr/bin/env python
from __future__ import division, print_function

import os 
import os.path
import copy 
import collections
import abc

from subprocess import Popen, PIPE
from pymatgen.io.abinitio.utils import File

__all__ = ['AbinitJobScript',]

# =========================================================================== #

class ShellEditor(object):
    "Simple editor that simplifies the writing of shell scripts"
    _shell = '/bin/bash'

    @property
    def shell(self):
        return self._shell

    def _add(self, text, pre=""):
        if isinstance(text, str):
            self._lines.append(pre+text)
        else:
            self._lines.extend([pre + t for t in text])

    def shabang(self):
        self._lines = ['#!' + self.shell,]

    def declare_var(self, key, val):
        "Return a lines setting a variable."
        line = key + '=' + str(val)
        self._add(line)

    def declare_vars(self, d):
        for k,v in d.items():
            self.declare_var(k, v)

    def export_envar(self, key, val):
        line = "export " + key + "=" + str(val)
        self._add(line)

    def export_envars(self, env):
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

    def get_script_str(self):
        "Returns a string with the script and reset the editor"
        s = "\n".join(l for l in self._lines)
        self.clear()
        return s

    def clear(self):
        del self._lines

##########################################################################################

class TaskLauncher(object):
    __metaclass__ = abc.ABCMeta

    @staticmethod
    def from_task(task):
         # shell, slurm, pbs, shell-fw, slurm-fw, pbs-fw
         queue_type, have_fw = task.runmode["launcher"], False
         if "-" in queue_type:
            queue_type, fw = queue_type.split("-")
            have_fw = (fw == "fw")

         if have_fw:
            raise ValueError("Firework not yet supported")
                                                        
         # Find the subclass to instanciate.
         for c in TaskLauncher.__subclasses__():
            if c._queue_type == queue_type:
                cls = c

        # Find number of processors ....
        #runhints = self.get_runhints()

         return cls(task.jobfile_path, 
                    abi_stdin  = task.files_file.path, 
                    abi_stdout = task.log_file.path,
                    abi_stderr = task.stderr_file.path,
                    **task.runmode
                   )
                                                        
    def __init__(self, path, **kwargs):
        """
        Args:
            path:
                Path to the script file.
            abi_stdin:
                The input file to feed in the executable as the standard input. Mandatory.
            abi_stdout:
                The file into which the standard output is redirected. Default is 'log'.
            abi_stderr:
                The file into which the standard error is redirected. Default is 'stderr'.
            bindir:
                The directory in which to look for binaries. Default is "".
            modules:
                The modules which will be loaded with 'module load'. Default is None
            pre_lines:
                Lines before the main execution. Default is []
            post_lines:
                Lines after the main execution. Default is [].
            exe:
                The name of the binary to be executed (located in bindir or in $PATH if not bindir). Default is abinit.
            mpirun:
                The mpi runner.
                E.g. 'mpiexec -npernode 6'. Default is None.
            mpi_ncpus
                Number of MPI processes (default 1)
            omp_ncpus
                Number of OpenMP threads (default 0, i.e not used)
            vars:
                Dictionary {varname:varvalue} containing variable declaration.
        """
        self.jobfile      = File(path)
        self.abi_stdin    = kwargs.pop("abi_stdin",  "abi.files")
        self.abi_stdout   = kwargs.pop("abi_stdout", "abi.log")
        self.abi_stderr   = kwargs.pop("abi_stderr", "abi.err")
        self.bindir       = kwargs.pop("bindir", "")
        self.exe          = os.path.join(self.bindir, kwargs.pop("exe", "abinit"))
        self.mpirun       = kwargs.pop("mpirun", "")
        self.mpi_ncpus    = kwargs.pop("mpi_ncpus", 1)
        self.omp_ncpus    = kwargs.pop("omp_ncpus", 0)
        self.vars         = kwargs.pop("vars", {})
        self.modules      = kwargs.pop("modules", [])
        self.envars       = kwargs.pop("envars", {})
        self.pre_lines    = kwargs.pop("pre_lines", [])
        self.post_lines   = kwargs.pop("post_lines", [])
        self.queue_params = kwargs.pop("queue_params", {})

    def get_script_str(self):
        "Returns a string wth the script."
        se = ShellEditor()

        se.shabang()

        # Submission instructions for the queue manager.
        se.add_lines(self.get_queue_header())

        se.add_comment("Variable declarations")
        vars = self.vars.copy()
        vars.update({
                "EXECUTABLE": self.exe,
                "ABI_STDIN" : self.abi_stdin,
                "ABI_STDOUT": self.abi_stdout,
                "ABI_STDERR": self.abi_stderr,
                "MPIRUN"    : self.mpirun,
                "MPI_NCPUS" : self.mpi_ncpus,
               })

        se.declare_vars(vars)

        se.add_comment("Modules")
        se.load_modules(self.modules)

        se.add_comment("Environment")
        se.export_envars(self.envars)

        if self.omp_ncpus:
            se.export_envars("OMP_NUM_THREADS", self.omp_ncpus)

        se.add_comment("Commands before execution")
        se.add_lines(self.pre_lines)

        # Execution line
        # TODO: better treatment of mpirun syntax.
        if self.mpirun:
            se.add_line('$MPIRUN -n $MPI_NCPUS $EXECUTABLE < $ABI_STDIN > $ABI_STDOUT 2> $ABI_STDERR')
        else:
            se.add_line('$EXECUTABLE < $ABI_STDIN > $ABI_STDOUT 2> $ABI_STDERR')

        se.add_comment("Commands after execution")
        se.add_lines(self.post_lines)

        return se.get_script_str()

    def write_script(self, overwrite=True):
        "Writes the script file."
        #if not self.jobfile.exists:
        #    os.makedirs(self.jobfile.dirname)

        #if self.jobfile.exists and not overwrite:
        #    raise ValueError("%s already exists, cannot overwrite" % self.jobfile.path)

        with open(self.jobfile.path, 'w') as f:
            f.write(self.get_script_str())

    @abc.abstractmethod
    def get_queue_header(self):
        "Return a list of string with the options passes to the queue manager"

    @abc.abstractmethod
    def launch(self, *args, **kwargs):
        "Run the script in a subprocess, returns a Popen-like object."

##########################################################################################
# Concrete classes

class ShellLauncher(TaskLauncher):
    _queue_type = "shell"

    def get_queue_header(self):
        return []

    def launch(self, *args, **kwargs):
        self.write_script()
        return Popen(("/bin/bash", self.jobfile.path), cwd=self.jobfile.dirname, stderr=PIPE)

##########################################################################################

class SlurmLauncher(TaskLauncher):
    _queue_type = "slurm"

    def get_queue_header(self):
        raise NotImplementedError()

    def launch(self, *args, **kwargs):
        """
        Run the script in a subprocess, returns a Popen-like object.
        This is a very delicate part. Hopefully fireworks will help solve the problem!
        """
        raise NotImplementedError()

##########################################################################################

if __name__ == "__main__":
    se = ShellEditor()
    se.shabang()
    se.declare_var("FOO", "BAR")
    se.add_emptyline()
    se.add_comment("This is a comment")
    se.declare_vars({"FOO1": "BAR1"})
    se.load_modules(["module1", "module2"])
    print(se.get_script_str())

    #launcher = TaskLauncher.from_task(self)

    launcher = ShellLauncher("job.sh")
    print(launcher.get_script_str())
    #process = job.launch()
