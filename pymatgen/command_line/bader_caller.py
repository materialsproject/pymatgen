#!/usr/bin/env python

"""
This module implements an interface to the Henkelmann et al.'s excellent
Fortran code for calculating a Bader charge analysis.

This module depends on a compiled bader executable available in the path.
Please download the library at http://theory.cm.utexas.edu/vasp/bader/ and
follow the instructions to compile the executable.

If you use this module, please cite the following:

G. Henkelman, A. Arnaldsson, and H. Jonsson, "A fast and robust algorithm for
Bader decomposition of charge density", Comput. Mater. Sci. 36, 254-360 (2006).
"""

from __future__ import division

__author__ = "shyuepingong"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Beta"
__date__ = "4/5/13"

import os
import subprocess
import tempfile
import shutil


class BaderAnalysis(object):
    """
    Bader analysis for a CHGCAR.

    .. attribute: data

        Atomic data parsed from bader analysis. Essentially a list of dicts
        of the form::

        [
            {
                "dist": 8.769,
                "min": 0.8753,
                "charge": 7.4168,
                "y": 1.1598,
                "x": 0.0079,
                "z": 0.8348
            },
            ...
        ]

    .. attribute: vacuum_volume

        Vacuum volume of the Bader analysis.

    .. attribute: vacuum_charge

        Vacuum charge of the Bader analysis.

    .. attribute: nelectrons

        Number of electrons of the Bader analysis.
    """

    def __init__(self, filename):
        """
        Args:
            filename:
                The filename of the CHGCAR.
        """
        temp_dir = tempfile.mkdtemp()
        try:
            shutil.copy(filename, os.path.join(temp_dir, "CHGCAR"))
            current_dir = os.getcwd()
            os.chdir(temp_dir)
            rs = subprocess.Popen(["bader", "CHGCAR"],
                                  stdout=subprocess.PIPE,
                                  stdin=subprocess.PIPE, close_fds=True)
            rs.communicate()
            data = []
            with open("ACF.dat") as f:
                raw = f.readlines()
                headers = [s.lower() for s in raw.pop(0).split()]
                raw.pop(0)
                while True:
                    l = raw.pop(0).strip()
                    if l.startswith("-"):
                        break
                    vals = map(float, l.split()[1:])
                    data.append(dict(zip(headers[1:], vals)))
                for l in raw:
                    toks = l.strip().split(":")
                    if toks[0] == "VACUUM CHARGE":
                        self.vacuum_charge = float(toks[1])
                    elif toks[0] == "VACUUM VOLUME":
                        self.vacuum_volume = float(toks[1])
                    elif toks[0] == "NUMBER OF ELECTRONS":
                        self.nelectrons = float(toks[1])
            self.data = data
            os.chdir(current_dir)
        except Exception as ex:
            print str(ex)
        finally:
            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    a = BaderAnalysis("../../test_files/CHGCAR.noncubic")