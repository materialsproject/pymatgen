# coding: utf-8
"""Wrappers for ABINIT main executables"""
from __future__ import unicode_literals, division, print_function

import os

from subprocess import Popen, PIPE
from monty.os.path import which
from pymatgen.util.string_utils import list_strings
from six.moves import map, cStringIO

import logging
logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

__all__ = [
    "Mrgscr",
    "Mrggkk",
    "Mrgddb",
]


class ExecError(Exception):
    """Error class raised by :class`ExecWrapper`"""


class ExecWrapper(object):
    """This class runs an executable in a subprocess."""
    Error = ExecError

    def __init__(self, executable=None, verbose=0):
        """
        Args:
            executable: path to the executable.
            verbose: Verbosity level.
        """
        if executable is None:
            executable = self.name

        self.executable = which(executable)

        self.verbose = int(verbose)

        if self.executable is None:
            raise self.Error("Cannot find %s in $PATH\n Use export PATH=/dir_with_exec:$PATH" % executable)

        assert os.path.basename(self.executable) == self.name

    def __str__(self):
        return "%s" % self.executable

    def set_mpi_runner(self, mpi_runner="mpirun"):
        # TODO better treatment of mpirunner syntax.
        self._mpi_runner = mpi_runner

    @property
    def mpi_runner(self):
        try:
            return self._mpi_runner
        except AttributeError:
            return ""

    @property
    def name(self):
        return self._name

    def execute(self, cwd=None):
        # Try to execute binary without and with mpirun.
        try:
            self._execute(cwd=cwd, with_mpirun=True)
        except self.Error:
            self._execute(cwd=cwd, with_mpirun=False)

    def _execute(self, cwd=None, with_mpirun=False):
        """
        Execute the executable in a subprocess.
        """
        args = [self.executable, "<", self.stdin_fname, ">", self.stdout_fname, "2>", self.stderr_fname]

        if self.mpi_runner and with_mpirun:
            args.insert(0, self.mpi_runner)

        self.cmd_str = " ".join(args)

        p = Popen(self.cmd_str, shell=True, stdout=PIPE, stderr=PIPE, cwd=cwd)

        self.stdout_data, self.stderr_data = p.communicate()

        self.returncode = p.returncode

        if self.returncode != 0:
            with open(self.stdout_fname, "r") as out, open(self.stderr_fname, "r") as err:
                self.stdout_data = out.read()
                self.stderr_data = err.read()

            if self.verbose > 3:
                print("*** stdout: ***\n", self.stdout_data)
                print("*** stderr  ***\n", self.stderr_data)

            raise self.Error("%s returned %s\n cmd_str: %s" % (self, self.returncode, self.cmd_str))


class Mrgscr(ExecWrapper):
    _name = "mrgscr"

    def merge_qpoints(self, files_to_merge, out_prefix, cwd=None):
        """
        Execute mrgscr in a subprocess to merge files_to_merge. Produce new file with prefix out_prefix
        If cwd is not None, the child's current directory will be changed to cwd before it is executed.
        """
        # We work with absolute paths.
        files_to_merge = [os.path.abspath(s) for s in list_strings(files_to_merge)]
        nfiles = len(files_to_merge)

        if self.verbose:
            print("Will merge %d files with output_prefix %s" % (nfiles, out_prefix))
            for (i, f) in enumerate(files_to_merge):
                print(" [%d] %s" % (i, f))

        if nfiles == 1:
            raise self.Error("merge_qpoints does not support nfiles == 1")

        self.stdin_fname, self.stdout_fname, self.stderr_fname = (
            "mrgscr.stdin", "mrgscr.stdout", "mrgscr.stderr")

        if cwd is not None:
            self.stdin_fname, self.stdout_fname, self.stderr_fname = \
                map(os.path.join, 3 * [cwd], [self.stdin_fname, self.stdout_fname, self.stderr_fname])

        inp = cStringIO()

        inp.write(str(nfiles) + "\n")     # Number of files to merge.
        inp.write(out_prefix + "\n")      # Prefix for the final output file:

        for filename in files_to_merge:
            inp.write(filename + "\n")   # List with the files to merge.

        inp.write("1\n")                 # Option for merging q-points.

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)

        self.execute(cwd=cwd)


class Mrggkk(ExecWrapper):
    _name = "mrggkk"

    def merge(self, gswfk_file, dfpt_files, gkk_files, out_gkk, binascii=0, cwd=None):
        """
        Merge GGK files, return the absolute path of the new database.

        Args:
            gswfk_file: Ground-state WFK filename
            dfpt_files: List of 1WFK files to merge.
            gkk_files: List of GKK files to merge.
            out_gkk: Name of the output GKK file
            binascii: Integer flat. 0 --> binary output, 1 --> ascii formatted output
            cwd: Directory where the subprocess will be executed.
        """
        raise NotImplementedError("This method should be tested")

        out_gkk = out_gkk if cwd is None else os.path.join(os.path.abspath(cwd), out_gkk)

        # We work with absolute paths.
        gswfk_file = absath(gswfk_file)
        dfpt_files = [os.path.abspath(s) for s in list_strings(dfpt_files)]
        gkk_files = [os.path.abspath(s) for s in list_strings(gkk_files)]

        print("Will merge %d 1WF files, %d GKK file in output %s" %
              (len(dfpt_nfiles), len_gkk_files, out_gkk))

        if self.verbose:
            for i, f in enumerate(dfpt_files): print(" [%d] 1WF %s" % (i, f))
            for i, f in enumerate(gkk_files): print(" [%d] GKK %s" % (i, f))

        self.stdin_fname, self.stdout_fname, self.stderr_fname = (
            "mrggkk.stdin", "mrggkk.stdout", "mrggkk.stderr")

        if cwd is not None:
            self.stdin_fname, self.stdout_fname, self.stderr_fname = \
                map(os.path.join, 3 * [cwd], [self.stdin_fname, self.stdout_fname, self.stderr_fname])

        inp = cStringIO()

        inp.write(out_gkk + "\n")        # Name of the output file
        inp.write(str(binascii) + "\n")  # Integer flag: 0 --> binary output, 1 --> ascii formatted output
        inp.write(gswfk_file + "\n")     # Name of the groud state wavefunction file WF

        #dims = len(dfpt_files, gkk_files, ?)
        dims = " ".join([str(d) for d in dims])
        inp.write(dims + "\n")             # Number of 1WF, of GKK files, and number of 1WF files in all the GKK files

        # Names of the 1WF files...
        for fname in dfpt_files:
            inp.write(fname + "\n")

        # Names of the GKK files...
        for fname in gkk_files:
            inp.write(fname + "\n")

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)

        self.execute(cwd=cwd)

        return out_gkk


class Mrgddb(ExecWrapper):
    _name = "mrgddb"

    def merge(self, ddb_files, out_ddb, description, cwd=None):
        """Merge DDB file, return the absolute path of the new database."""
        # We work with absolute paths.
        ddb_files = [os.path.abspath(s) for s in list_strings(ddb_files)]

        out_ddb = out_ddb if cwd is None else os.path.join(os.path.abspath(cwd), out_ddb)

        print("Will merge %d files into output DDB %s" % (len(ddb_files), out_ddb))
        if self.verbose:
            for i, f in enumerate(ddb_files):
                print(" [%d] %s" % (i, f))

        # Handle the case of a single file since mrgddb uses 1 to denote GS files!
        if len(ddb_files) == 1:
            with open(ddb_files[0], "r") as inh, open(out_ddb, "w") as out:
                for line in inh:
                    out.write(line)
            return out_ddb

        self.stdin_fname, self.stdout_fname, self.stderr_fname = "mrgddb.stdin", "mrgddb.stdout", "mrgddb.stderr"

        if cwd is not None:
            self.stdin_fname, self.stdout_fname, self.stderr_fname = \
                map(os.path.join, 3 * [cwd], [self.stdin_fname, self.stdout_fname, self.stderr_fname])

        inp = cStringIO()

        inp.write(out_ddb + "\n")              # Name of the output file.
        inp.write(str(description) + "\n")     # Description.
        inp.write(str(len(ddb_files)) + "\n")  # Number of input DDBs.

        # Names of the DDB files.
        for fname in ddb_files:
            inp.write(fname + "\n")

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)

        self.execute(cwd=cwd)

        return out_ddb
