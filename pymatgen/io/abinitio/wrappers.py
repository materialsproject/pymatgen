from __future__ import division, print_function

import os
import cStringIO as StringIO

from subprocess import Popen, PIPE
from pymatgen.util.io_utils import which
from pymatgen.util.string_utils import list_strings

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
    "Anaddb",
]


class ExecWrapper(object):
    """This class runs an executable in a subprocess."""
    def __init__(self, executable=None, verbose=0):
        """
        Args:
            executable:
                path to the executable.
            verbose:
                Verbosity level.
        """
        if executable is None: 
            executable = self.name

        self.executable = which(executable)

        self.verbose = int(verbose)

        if self.executable is None:
            msg = "Cannot find executable %s is PATH\n Use export PATH=/dir_with_exec:$PATH" % executable
            raise self.Error(msg)

        assert os.path.basename(self.executable) == self._name

    def __str__(self):
        return "%s" % self.executable

    @property
    def name(self):
        return self._name

    def execute(self, cwd=None, **kwargs):
        """Execute the executable in a subprocess."""
        args = [self.executable, "<", self.stdin_fname, ">", self.stdout_fname, "2>", self.stderr_fname]

        self.cmd_str = " ".join(args)

        p = Popen(self.cmd_str, shell=True, stdout=PIPE, stderr=PIPE, cwd=cwd)

        (self.stdout_data, self.stderr_data) = p.communicate()

        self.returncode = p.returncode

        if self.returncode != 0:
            with open(self.stdout_fname, "r") as out, open(self.stderr_fname, "r") as err:
                self.stdout_data = out.read()
                self.stderr_data = err.read()

            if self.verbose:
                print("*** stdout: ***\n", self.stdout_data)
                print("*** stderr  ***\n", self.stderr_data)

            raise self.Error("%s returned %s\n cmd_str: %s" % (self, self.returncode, self.cmd_str))


class MrgscrError(Exception):
    """Error class for Mrgscr"""


class Mrgscr(ExecWrapper):
    _name = "mrgscr"

    Error = MrgscrError

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

        inp = StringIO.StringIO()

        inp.write(str(nfiles) + "\n")      # Number of files to merge.
        inp.write(out_prefix + "\n")      # Prefix for the final output file:

        for filename in files_to_merge:
            inp.write(filename + "\n")   # List with the files to merge.

        inp.write("1\n")                 # Option for merging q-points.

        inp.seek(0)
        self.stdin_data = [s for s in inp]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)

        try:
            self.execute(cwd=cwd)
        except self.Error:
            raise


class MrggkkError(Exception):
    """Error class for Mrggkk."""


class Mrggkk(ExecWrapper):
    _name = "mrggkk"

    Error = MrggkkError

    def merge(self, gswfk_file, dfpt_files, gkk_files, out_gkk, binascii=0, cwd=None):
        """
        Merge GGK files, return the absolute path of the new database.

        Args:
            gswfk_file: 
                Ground-state WFK filename
            dfpt_files:
                List of 1WFK files to merge.
            gkk_files:
                List of GKK files to merge.
            out_gkk:
                Name of the output GKK file
            binascii:
                Integer flat. 0 --> binary output, 1 --> ascii formatted output
            cwd:
                Directory where the subprocess will be executed.
        """
        raise NotImplementedError("This method should be tested")

        out_gkk = out_gkk if cwd is None else os.path.join(os.path.abspath(cwd), out_gkk)

        # We work with absolute paths.
        gswfk_file = absath(gswfk_file)
        dfpt_files = [os.path.abspath(s) for s in list_strings(dfpt_files)]
        gkk_files = [os.path.abspath(s) for s in list_strings(gkk_files)]

        if self.verbose:
            print("Will merge %d 1WF files, %d GKK file in output %s" %
                  (len(dfpt_nfiles), (len_gkk_files), out_gkk))

            for (i, f) in enumerate(dfpt_files): 
                print(" [%d] 1WF %s" % (i, f))

            for (i, f) in enumerate(gkk_files):  
                print(" [%d] GKK %s" % (i, f))

        self.stdin_fname, self.stdout_fname, self.stderr_fname = (
            "mrggkk.stdin", "mrggkk.stdout", "mrggkk.stderr")

        if cwd is not None:
            self.stdin_fname, self.stdout_fname, self.stderr_fname = \
                map(os.path.join, 3 * [cwd], [self.stdin_fname, self.stdout_fname, self.stderr_fname])

        inp = StringIO.StringIO()

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

        inp.seek(0)
        self.stdin_data = [s for s in inp]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)

        try:
            self.execute(cwd=cwd)
        except self.Error:
            raise

        return out_gkk


class MrgddbError(Exception):
    """Error class for Mrgddb."""


class Mrgddb(ExecWrapper):
    _name = "mrgddb"

    Error = MrgddbError

    def merge(self, ddb_files, out_ddb, description, cwd=None):
        """Merge DDB file, return the absolute path of the new database."""

        # We work with absolute paths.
        ddb_files = [os.path.abspath(s) for s in list_strings(ddb_files)]

        out_ddb = out_ddb if cwd is None else os.path.join(os.path.abspath(cwd), out_ddb)

        if self.verbose:
            print("Will merge %d files into output DDB %s" % (len(ddb_files), out_ddb))
            for (i, f) in enumerate(ddb_files):
                print(" [%d] %s" % (i, f))

        # Handle the case of a single file since mrgddb uses 1 to denote GS files!
        if len(ddb_files) == 1:
            with open(ddb_files[0], "r") as inh, open(out_ddb, "w") as out:
                for line in inh:
                    out.write(line)
            return out_ddb

        self.stdin_fname, self.stdout_fname, self.stderr_fname = (
            "mrgddb.stdin", "mrgddb.stdout", "mrgddb.stderr")

        if cwd is not None:
            self.stdin_fname, self.stdout_fname, self.stderr_fname = \
                map(os.path.join, 3 * [cwd], [self.stdin_fname, self.stdout_fname, self.stderr_fname])

        inp = StringIO.StringIO()

        inp.write(out_ddb + "\n")            # Name of the output file.
        inp.write(str(description) + "\n")     # Description.
        inp.write(str(len(ddb_files)) + "\n")  # Number of input DDBs.

        # Names of the DDB files.
        for fname in ddb_files:
            inp.write(fname + "\n")

        inp.seek(0)
        self.stdin_data = [s for s in inp]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)

        try:
            self.execute(cwd=cwd)
        except self.Error:
            raise

        return out_ddb


class AnaddbError(Exception):
    """Error class for Anaddb."""


class Anaddb(ExecWrapper):
    _name = "anaddb"

    Error = AnaddbError

    #def make_stdin(self):
    #    # Files file
    #    inp = StringIO.StringIO()
    #    inp.write(self.input_fname + "\n")     # Input file.
    #    inp.write(self.stdout_fname + "\n")    # Output file.
    #    inp.write(ddb_file + "\n")             # DDB file
    #    inp.write("dummy_band2eps" + "\n")
    #    inp.write("dummy1" + "\n")
    #    inp.write("dummy2" + "\n")
    #    inp.write("dummy3" + "\n")
    #    inp.seek(0)
    #    self.stdin_data = [s for s in inp]

    def diagonalize_1q(self, ddb_file, cwd=None):

        # We work with absolute paths.
        ddb_file = os.path.abspath(ddb_file)

        self.stdin_fname, self.input_fname, self.stdout_fname, self.stderr_fname = (
            "anaddb.stdin", "anaddb.input", "anaddb.stdout", "anaddb.stderr")

        if cwd is not None:
            self.stdin_fname, self.input_fname, self.stdout_fname, self.stderr_fname = \
                map(os.path.join, 3 * [cwd], [self.stdin_fname, self.inp_fname, self.stdout_fname, self.stderr_fname])

        # Files file
        inp = StringIO.StringIO()

        inp.write(self.input_fname + "\n")     # Input file.
        inp.write(self.stdout_fname + "\n")    # Output file.
        inp.write(ddb_file + "\n")             # DDB file
        inp.write("dummy_band2eps" + "\n")
        inp.write("dummy1" + "\n")
        inp.write("dummy2" + "\n")
        inp.write("dummy3" + "\n")

        inp.seek(0)
        self.stdin_data = [s for s in inp]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)

        # Get the q-point from the DDB file
        with open(ddb_file, "r") as fh:
            nfound = 0
            tag = " qpt "
            for line in fh:
                print(line)
                if line.startswith(tag):
                    nfound += 1
                    # Coordinates of the q-points.
                    qcoords_str = line.split()[1:4]
                    #qcoords_str = [ s.replace("D", "E") for s in qcoords_str]
                    qpoint = map(float, qcoords_str)

        if nfound != 1:
            raise self.Error("Found %s occurrences of tag %s in file %s" % (nfound, tag, ddb_file))

        # Write simple input file for the anaddb code.
        with open(self.input_fname, "w") as inp:
            inp.write('# Flags\n')
            inp.write(' ifcflag   1 # Interatomic force constant flag\n\n')
            inp.write('# Wavevector grid number 1 (coarse grid, from DDB)\n\n')
            inp.write(' brav    1         # Bravais Lattice : 1-S.C., 2-F.C., 3-B.C., 4-Hex.\n')
            inp.write(' ngqpt   1  1  1   # Q-mesh\n')
            inp.write(' nqshft  1         # number of q-shifts\n')
            inp.write(' q1shft  %f %f %f' % tuple(qpoint))

            #inp.write('# Effective charges
            #inp.write('     asr   1     ! Acoustic Sum Rule. 1 => imposed asymetrically
            #inp.write('  chneut   1     ! Charge neutrality requirement for effective charges.
            #inp.write('# Interatomic force constant info
            #inp.write('  dipdip  1      ! Dipole-dipole interaction treatment
            #inp.write('  ifcana  1      ! Analysis of the IFCs
            #inp.write('  ifcout 20      ! Number of IFC's written in the output, per atom
            #inp.write('  natifc  1      ! Number of atoms in the cell for which ifc's are analysed
            #inp.write('   atifc  1      ! List of atoms
            #inp.write('
            #inp.write('# This line added when defaults were changed (v5.3) to keep the previous, old behaviour
            #inp.write('#  symdynmat 0

        if self.verbose:
            print("Will diagonalize DDB file : %s" % ddb_file)

        try:
            self.execute(cwd=cwd)
        except self.Error:
            raise

        # Get frequencies from the output file
        # TODO
        #with open(self.stdout_fname, "r") as out:
            #print(out.readlines())
            #for line in out:
            #    if line: raise
            #return frequencies
