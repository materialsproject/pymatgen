# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""Wrappers for ABINIT main executables"""
from __future__ import unicode_literals, division, print_function

import os

from monty.string import list_strings
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
    "Mrgdvdb",
]


class ExecError(Exception):
    """Error class raised by :class:`ExecWrapper`"""


class ExecWrapper(object):
    """Base class that runs an executable in a subprocess."""
    Error = ExecError

    def __init__(self, manager=None, executable=None, verbose=0):
        """
        Args:
            manager: :class:`TaskManager` object responsible for the submission of the jobs.
                if manager is None, the default manager is used.
            executable: path to the executable.
            verbose: Verbosity level.
        """
        from .tasks import TaskManager
        self.manager = manager if manager is not None else TaskManager.from_user_config()
        self.manager = self.manager.to_shell_manager(mpi_procs=1)

        self.executable = executable if executable is not None else self.name
        assert os.path.basename(self.executable) == self.name
        self.verbose = int(verbose)

    def __str__(self):
        return "%s" % self.executable

    @property
    def name(self):
        return self._name

    def execute(self, workdir, exec_args=None):
        # Try to execute binary without and with mpirun.
        try:
            return self._execute(workdir, with_mpirun=True, exec_args=exec_args)
        except self.Error:
            return self._execute(workdir, with_mpirun=False, exec_args=exec_args)

    def _execute(self, workdir, with_mpirun=False, exec_args=None):
        """
        Execute the executable in a subprocess inside workdir.

        Some executables fail if we try to launch them with mpirun.
        Use with_mpirun=False to run the binary without it.
        """
        qadapter = self.manager.qadapter
        if not with_mpirun: qadapter.name = None

        script = qadapter.get_script_str(
            job_name=self.name,
            launch_dir=workdir,
            executable=self.executable,
            qout_path="qout_file.path",
            qerr_path="qerr_file.path",
            stdin=self.stdin_fname,
            stdout=self.stdout_fname,
            stderr=self.stderr_fname,
            exec_args=exec_args
        )

        # Write the script.
        script_file = os.path.join(workdir, "run" + self.name + ".sh")
        with open(script_file, "w") as fh:
            fh.write(script)
            os.chmod(script_file, 0o740)

        qjob, process = qadapter.submit_to_queue(script_file)
        self.stdout_data, self.stderr_data = process.communicate()
        self.returncode = process.returncode
        #raise self.Error("%s returned %s\n cmd_str: %s" % (self, self.returncode, self.cmd_str))

        return self.returncode


class Mrgscr(ExecWrapper):
    _name = "mrgscr"

    def merge_qpoints(self, workdir, files_to_merge, out_prefix):
        """
        Execute mrgscr inside directory `workdir` to merge `files_to_merge`.
        Produce new file with prefix `out_prefix`
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

        self.stdin_fname, self.stdout_fname, self.stderr_fname = \
            map(os.path.join, 3 * [workdir], ["mrgscr.stdin", "mrgscr.stdout", "mrgscr.stderr"])

        inp = cStringIO()
        inp.write(str(nfiles) + "\n")     # Number of files to merge.
        inp.write(out_prefix + "\n")      # Prefix for the final output file:

        for filename in files_to_merge:
            inp.write(filename + "\n")   # List with the files to merge.

        inp.write("1\n")                 # Option for merging q-points.

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)

        self.execute(workdir)


class Mrggkk(ExecWrapper):
    _name = "mrggkk"

    def merge(self, workdir, gswfk_file, dfpt_files, gkk_files, out_gkk, binascii=0):
        """
        Merge GGK files, return the absolute path of the new database.

        Args:
            gswfk_file: Ground-state WFK filename
            dfpt_files: List of 1WFK files to merge.
            gkk_files: List of GKK files to merge.
            out_gkk: Name of the output GKK file
            binascii: Integer flat. 0 --> binary output, 1 --> ascii formatted output
        """
        raise NotImplementedError("This method should be tested")
        #out_gkk = out_gkk if cwd is None else os.path.join(os.path.abspath(cwd), out_gkk)

        # We work with absolute paths.
        gswfk_file = os.path.absath(gswfk_file)
        dfpt_files = [os.path.abspath(s) for s in list_strings(dfpt_files)]
        gkk_files = [os.path.abspath(s) for s in list_strings(gkk_files)]

        print("Will merge %d 1WF files, %d GKK file in output %s" %
              (len(dfpt_files), len(gkk_files), out_gkk))

        if self.verbose:
            for i, f in enumerate(dfpt_files): print(" [%d] 1WF %s" % (i, f))
            for i, f in enumerate(gkk_files): print(" [%d] GKK %s" % (i, f))

        self.stdin_fname, self.stdout_fname, self.stderr_fname = \
            map(os.path.join, 3 * [workdir], ["mrggkk.stdin", "mrggkk.stdout", "mrggkk.stderr"])

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

        self.execute(workdir)

        return out_gkk


class Mrgddb(ExecWrapper):
    _name = "mrgddb"

    def merge(self, workdir, ddb_files, out_ddb, description, delete_source_ddbs=True):
        """Merge DDB file, return the absolute path of the new database in workdir."""
        # We work with absolute paths.
        ddb_files = [os.path.abspath(s) for s in list_strings(ddb_files)]
        if not os.path.isabs(out_ddb):
            out_ddb = os.path.join(os.path.abspath(workdir), os.path.basename(out_ddb))

        if self.verbose:
            print("Will merge %d files into output DDB %s" % (len(ddb_files), out_ddb))
            for i, f in enumerate(ddb_files):
                print(" [%d] %s" % (i, f))

        # Handle the case of a single file since mrgddb uses 1 to denote GS files!
        if len(ddb_files) == 1:
            with open(ddb_files[0], "r") as inh, open(out_ddb, "w") as out:
                for line in inh:
                    out.write(line)
            return out_ddb

        self.stdin_fname, self.stdout_fname, self.stderr_fname = \
            map(os.path.join, 3 * [os.path.abspath(workdir)], ["mrgddb.stdin", "mrgddb.stdout", "mrgddb.stderr"])

        inp = cStringIO()
        inp.write(out_ddb + "\n")              # Name of the output file.
        inp.write(str(description) + "\n")     # Description.
        inp.write(str(len(ddb_files)) + "\n")  # Number of input DDBs.

        # Names of the DDB files.
        for fname in ddb_files:
            inp.write(fname + "\n")

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "wt") as fh:
            fh.writelines(self.stdin_data)

        retcode = self.execute(workdir, exec_args=['--nostrict'])
        if retcode == 0 and delete_source_ddbs:
            # Remove ddb files.
            for f in ddb_files:
                try:
                    os.remove(f)
                except IOError:
                    pass

        return out_ddb


class Mrgdvdb(ExecWrapper):
    _name = "mrgdv"

    def merge(self, workdir, pot_files, out_dvdb, delete_source=True):
        """
        Merge POT files containing 1st order DFPT potential
        return the absolute path of the new database in workdir.
        """
        # We work with absolute paths.
        pot_files = [os.path.abspath(s) for s in list_strings(pot_files)]
        if not os.path.isabs(out_dvdb):
            out_dvdb = os.path.join(os.path.abspath(workdir), os.path.basename(out_dvdb))

        if self.verbose:
            print("Will merge %d files into output DVDB %s" % (len(pot_files), out_dvdb))
            for i, f in enumerate(pot_files):
                print(" [%d] %s" % (i, f))

        # Handle the case of a single file since mrgddb uses 1 to denote GS files!
        if len(pot_files) == 1:
            with open(pot_files[0], "r") as inh, open(out_dvdb, "w") as out:
                for line in inh:
                    out.write(line)
            return out_dvdb

        self.stdin_fname, self.stdout_fname, self.stderr_fname = \
            map(os.path.join, 3 * [workdir], ["mrgdvdb.stdin", "mrgdvdb.stdout", "mrgdvdb.stderr"])

        inp = cStringIO()
        inp.write(out_dvdb + "\n")              # Name of the output file.
        inp.write(str(len(pot_files)) + "\n")  # Number of input POT files.

        # Names of the POT files.
        for fname in pot_files:
            inp.write(fname + "\n")

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "wt") as fh:
            fh.writelines(self.stdin_data)

        retcode = self.execute(workdir)
        if retcode == 0 and delete_source:
            # Remove pot files.
            for f in pot_files:
                try:
                    os.remove(f)
                except IOError:
                    pass

        return out_dvdb


class Cut3D(ExecWrapper):
    _name = "cut3d"

    def cut3d(self, cut3d_input, workdir):
        """
        Runs cut3d with a Cut3DInput

        Args:
            cut3d_input: a Cut3DInput object.
            workdir: directory where cut3d is executed.

        Returns:
            (string) absolute path to the standard output of the cut3d execution.
            (string) absolute path to the output filepath. None if output is required.
        """
        self.stdin_fname, self.stdout_fname, self.stderr_fname = \
            map(os.path.join, 3 * [os.path.abspath(workdir)], ["cut3d.stdin", "cut3d.stdout", "cut3d.stderr"])

        cut3d_input.write(self.stdin_fname)

        retcode = self._execute(workdir, with_mpirun=False)

        if retcode != 0:
            raise RuntimeError("Error while running cut3d.")

        output_filepath = cut3d_input.output_filepath

        if output_filepath is not None:
            if not os.path.isabs(output_filepath):
                output_filepath = os.path.abspath(os.path.join(workdir, output_filepath))

            if not os.path.isfile(output_filepath):
                raise RuntimeError("The file was not converted correctly.")

        return self.stdout_fname, output_filepath
