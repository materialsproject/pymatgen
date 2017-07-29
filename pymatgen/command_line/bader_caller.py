# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import os
import subprocess
import shutil
import warnings

from six.moves import map, zip
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.io.vasp.inputs import Potcar
from monty.os.path import which
from monty.dev import requires
from monty.tempfile import ScratchDir
from monty.json import MSONable

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

__author__ = "shyuepingong"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Beta"
__date__ = "4/5/13"


@requires(which("bader") or which("bader.exe"),
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
                "atomic_vol": 8.769,
                "min_dist": 0.8753,
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

    .. attribute: chgcar_ref
    
        Chgcar reference which calculated by AECCAR0 + AECCAR2.
        (See http://theory.cm.utexas.edu/henkelman/code/bader/ for details.)
    """

    def __init__(self, chgcar_filename, potcar_filename=None,
                 chgref_filename=None):
        """
        Initializes the Bader caller.

        Args:
            chgcar_filename (str): The filename of the CHGCAR.
            potcar_filename (str): Optional: the filename of the corresponding
                POTCAR file. Used for calculating the charge transfer. If
                None, the get_charge_transfer method will raise a ValueError.
            chgref_filename (str): Optional. The filename of the reference
                CHGCAR, which calculated by AECCAR0 + AECCAR2. (See
                http://theory.cm.utexas.edu/henkelman/code/bader/ for details.)
        """
        self.chgcar = Chgcar.from_file(chgcar_filename)
        self.potcar = Potcar.from_file(potcar_filename) \
            if potcar_filename is not None else None
        self.natoms = self.chgcar.poscar.natoms
        chgcarpath = os.path.abspath(chgcar_filename)
        chgrefpath = os.path.abspath(chgref_filename) if chgref_filename else None
        self.reference_used = True if chgref_filename else False
        with ScratchDir(".") as temp_dir:
            shutil.copy(chgcarpath, os.path.join(temp_dir, "CHGCAR"))
            args = ["bader", "CHGCAR"]
            if chgref_filename:
                shutil.copy(chgrefpath, os.path.join(temp_dir, "CHGCAR_ref"))
                args += ['-ref', 'CHGCAR_ref']
            rs = subprocess.Popen(args,
                                  stdout=subprocess.PIPE,
                                  stdin=subprocess.PIPE, close_fds=True)
            stdout, stderr = rs.communicate()
            if rs.returncode != 0:
                raise RuntimeError("bader exited with return code %d. "
                                   "Please check your bader installation."
                                   % rs.returncode)

            try:
                self.version = float(stdout.split()[5])
            except:
                self.version = -1  # Unknown
            if self.version < 1.0:
                warnings.warn('Your installed version of Bader is outdated, '
                              'calculation of vacuum charge may be incorrect.')

            data = []
            with open("ACF.dat") as f:
                raw = f.readlines()
                headers = ('x', 'y', 'z', 'charge', 'min_dist', 'atomic_vol')
                raw.pop(0)
                raw.pop(0)
                while True:
                    l = raw.pop(0).strip()
                    if l.startswith("-"):
                        break
                    vals = map(float, l.split()[1:])
                    data.append(dict(zip(headers, vals)))
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
        charges = [-self.get_charge_transfer(i)
                   for i in range(len(structure))]
        structure.add_oxidation_state_by_site(charges)
        return structure

    @property
    def result(self):

        summary = {
            "min_dist": [d['min_dist'] for d in self.data],
            "charge": [d['charge'] for d in self.data],
            "volume": [d['atomic_vol'] for d in self.data],
            "vacuum_charge": self.vacuum_charge,
            "vacuum_volume": self.vacuum_volume,
            "reference_used": self.reference_used,
            "bader_version": self.version,
        }

        if self.potcar:
            charge_transfer = [self.get_charge_transfer(i) for i in range(len(self.data))]
            summary['charge_transfer'] = charge_transfer

        return BaderResult(summary)


class BaderResult(MSONable):

    def __init__(self, summary):
        """
        A class to hold results of a BaderAnalysis. Simply a dict
        containing lists that can be used as site properties for
        the corresponding input structure.

        Use `from_path` or `from_objects` to conveniently obtain
        a BaderResult. To manually specify arguments and call the
        bader binary directly, use BaderAnalysis class instead.

        :param summary: BaderResult
        """

        for key in summary:
            setattr(self, key, summary[key])

    @classmethod
    def from_dict(cls, d):
        del d["@module"]
        del d["@class"]
        return cls(d)

    @classmethod
    def as_dict(self):
        d = self.__dict__
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_path(cls, path):
        """
        Convenience method to run Bader analysis on a folder containing
        typical VASP output files.

        This method will:

        1. Look for files CHGCAR, AECAR0, AECAR2, POTCAR or their gzipped
        counterparts.
        2. If AECCAR* files are present, constructs a temporary reference
        file as AECCAR0 + AECCAR2
        3. Runs Bader analysis twice: once for charge, and a second time
        for the charge difference (magnetization density).

        :param path: path to folder to search in
        :return: BaderResult
        """

        if os.path.isfile(os.path.join(path, 'CHGCAR')):
            chgcar = Chgcar.from_file(os.path.join(path, 'CHGCAR'))
        elif os.path.isfile(os.path.join(path, 'CHGCAR.gz')):
            chgcar = Chgcar.from_file(os.path.join(path, 'CHGCAR.gz'))
        else:
            raise IOError('Could not find CHGCAR.')

        if os.path.isfile(os.path.join(path, 'AECCAR0')):
            aeccar0 = Chgcar.from_file(os.path.join(path, 'AECCAR0'))
        elif os.path.isfile(os.path.join(path, 'AECCAR0.gz')):
            aeccar0 = Chgcar.from_file(os.path.join(path, 'AECCAR0.gz'))
        else:
            warnings.warn('Could not find AECCAR0, interpret Bader results with caution.')
            aeccar0 = None

        if os.path.isfile(os.path.join(path, 'AECCAR2')):
            aeccar2 = Chgcar.from_file(os.path.join(path, 'AECCAR2'))
        elif os.path.isfile(os.path.join(path, 'AECCAR2.gz')):
            aeccar2 = Chgcar.from_file(os.path.join(path, 'AECCAR2.gz'))
        else:
            warnings.warn('Could not find AECCAR2, interpret Bader results with caution.')
            aeccar2 = None

        if os.path.isfile(os.path.join(path, 'POTCAR')):
            potcar = Potcar.from_file(os.path.join(path, 'POTCAR'))
        elif os.path.isfile(os.path.join(path, 'POTCAR.gz')):
            potcar = Potcar.from_file(os.path.join(path, 'POTCAR.gz'))
        else:
            warnings.warn('Could not find POTCAR, cannot calculate charge transfer.')
            potcar = None

        return cls.from_objects(chgcar, potcar, aeccar0, aeccar2)


    @classmethod
    def from_objects(cls, chgcar, potcar=None, aeccar0=None, aeccar2=None):
        """
        Convenience method to run Bader analysis from a set
        of pymatgen Chgcar and Potcar objects.

        This method will:

        1. If aeccar objects are present, constructs a temporary reference
        file as AECCAR0 + AECCAR2
        2. Runs Bader analysis twice: once for charge, and a second time
        for the charge difference (magnetization density).

        :param chgcar: Chgcar object
        :param potcar: (optional) Potcar object
        :param aeccar0: (optional) Chgcar object from aeccar0 file
        :param aeccar2: (optional) Chgcar object from aeccar2 file
        :return: BaderResult
        """

        with ScratchDir(".") as temp_dir:

            if aeccar0 and aeccar2:
                # construct reference file
                chgref = aeccar0.linear_add(aeccar2)
                chgref_path = os.path.join(temp_dir, 'CHGCAR_ref')
                chgref.write_file(chgref_path)
            else:
                chgref_path = None

            chgcar.write_file('CHGCAR')
            chgcar_path = os.path.join(temp_dir, 'CHGCAR')

            if potcar:
                potcar.write_file('POTCAR')
                potcar_path = os.path.join(temp_dir, 'POTCAR')
            else:
                potcar_path = None

            ba = BaderAnalysis(chgcar_path, potcar_filename=potcar_path, chgref_filename=chgref_path)

            summary = {
                "min_dist": [d['min_dist'] for d in ba.data],
                "charge": [d['charge'] for d in ba.data],
                "volume": [d['atomic_vol'] for d in ba.data],
                "vacuum_charge": ba.vacuum_charge,
                "vacuum_volume": ba.vacuum_volume,
                "reference_used": True if chgref_path else False,
                "bader_version": ba.version,
            }

            if potcar:
                charge_transfer = [ba.get_charge_transfer(i) for i in range(len(ba.data))]
                summary['charge_transfer'] = charge_transfer

            if chgcar.is_spin_polarized:

                # write a CHGCAR containing magnetization density only
                chgcar.data['total'] = chgcar.data['diff']
                chgcar.is_spin_polarized = False
                chgcar.write_file('CHGCAR_mag')

                chgcar_mag_path = os.path.join(temp_dir, 'CHGCAR_mag')
                ba = BaderAnalysis(chgcar_mag_path, potcar_filename=potcar_path, chgref_filename=chgref_path)
                summary["magmom"] = [d['charge'] for d in ba.data]

            return cls(summary)