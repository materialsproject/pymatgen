# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

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

from six.moves import map
from six.moves import zip

__author__ = "shyuepingong"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Beta"
__date__ = "4/5/13"

import os
import subprocess
import shutil

from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.io.vasp.inputs import Potcar
from monty.os.path import which
from monty.dev import requires
from monty.tempfile import ScratchDir


@requires(which("bader"),
          "BaderAnalysis requires the executable bader to be in the path."
          " Please download the library at http://theory.cm.utexas"
          ".edu/vasp/bader/ and compile the executable.")
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

    .. attribute: chgcar

        Chgcar object associated with input CHGCAR file.

    .. attribute: potcar

        Potcar object associated with POTCAR used for calculation (used for
        calculating charge transferred).
    """


    def __init__(self, chgcar_filename, potcar_filename=None):
        """
        Initializes the Bader caller.

        Args:
            chgcar_filename: The filename of the CHGCAR.
            potcar_filename: Optional: the filename of the corresponding
                POTCAR file. Used for calculating the charge transfer. If
                None, the get_charge_transfer method will raise a ValueError.
        """
        self.chgcar = Chgcar.from_file(chgcar_filename)
        self.potcar = Potcar.from_file(potcar_filename) \
            if potcar_filename is not None else None
        self.natoms = self.chgcar.poscar.natoms
        chgcarpath = os.path.abspath(chgcar_filename)

        with ScratchDir(".") as temp_dir:
            shutil.copy(chgcarpath, os.path.join(temp_dir, "CHGCAR"))

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

    def get_charge(self, atom_index):
        """
        Convenience method to get the charge on a particular atom.

        Args:
            atom_index:
                Index of atom.

        Returns:
            Charge associated with atom from the Bader analysis.
        """
        return self.data[atom_index]["charge"]

    def get_charge_transfer(self, atom_index):
        """
        Returns the charge transferred for a particular atom. Requires POTCAR
        to be supplied.

        Args:
            atom_index:
                Index of atom.

        Returns:
            Charge transfer associated with atom from the Bader analysis.
            Given by final charge on atom - nelectrons in POTCAR for
            associated atom.
        """
        if self.potcar is None:
            raise ValueError("POTCAR must be supplied in order to calculate "
                             "charge transfer!")
        potcar_indices = []
        for i, v in enumerate(self.natoms):
            potcar_indices += [i] * v
        nelect = self.potcar[potcar_indices[atom_index]].nelectrons
        return self.data[atom_index]["charge"] - nelect

    def get_oxidation_state_decorated_structure(self):
        """
        Returns an oxidation state decorated structure.

        Returns:
            Returns an oxidation state decorated structure. Requires POTCAR
            to be supplied.
        """
        structure = self.chgcar.structure
        charges = [self.get_charge_transfer(i) for i in range(len(structure))]
        structure.add_oxidation_state_by_site(charges)
        return structure
