# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

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

import os
import subprocess
import shutil
import warnings
import glob
import numpy as np

from pymatgen.io.cube import Cube
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.io.vasp.inputs import Potcar
from monty.dev import requires
from monty.os.path import which
from monty.tempfile import ScratchDir
from monty.io import zopen

__author__ = "shyuepingong"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Beta"
__date__ = "4/5/13"

BADEREXE = which("bader") or which("bader.exe")


class BaderAnalysis:
    """
    Bader analysis for Cube files and VASP outputs.

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

    .. attribute: atomic_densities

        list of charge densities for each atom centered on the atom
        excess 0's are removed from the array to reduce the size of the array
        the charge densities are dicts with the charge density map,
        the shift vector applied to move the data to the center, and the original dimension of the charge density map
        charge:
            {
            "data": charge density array
            "shift": shift used to center the atomic charge density
            "dim": dimension of the original charge density map
            }
    """

    @requires(which("bader") or which("bader.exe"),
              "BaderAnalysis requires the executable bader to be in the path."
              " Please download the library at http://theory.cm.utexas"
              ".edu/vasp/bader/ and compile the executable.")
    def __init__(self,
                 chgcar_filename=None,
                 potcar_filename=None,
                 chgref_filename=None,
                 parse_atomic_densities=False,
                 cube_filename=None):
        """
        Initializes the Bader caller.

        Args:
            chgcar_filename (str): The filename of the CHGCAR.

            parse_atomic_densities (bool): Optional. turns on atomic partition of the charge density
                charge densities are atom centered

        """
        if not BADEREXE:
            raise RuntimeError(
                "BaderAnalysis requires the executable bader to be in the path."
                " Please download the library at http://theory.cm.utexas"
                ".edu/vasp/bader/ and compile the executable.")

        if not (cube_filename or chgcar_filename):
            raise ValueError("You must provide a file! Either a cube file or a CHGCAR")
        if cube_filename and chgcar_filename:
            raise ValueError("You cannot parse a cube and a CHGCAR at the same time!")

        self.parse_atomic_densities = parse_atomic_densities

        if chgcar_filename:
            fpath = os.path.abspath(chgcar_filename)
            self.is_vasp = True
            self.chgcar = Chgcar.from_file(chgcar_filename)
            self.structure = self.chgcar.structure
            self.potcar = Potcar.from_file(potcar_filename) \
                if potcar_filename is not None else None
            self.natoms = self.chgcar.poscar.natoms
            chgrefpath = os.path.abspath(chgref_filename) if chgref_filename else None
            self.reference_used = True if chgref_filename else False

            # List of nelects for each atom from potcar
            potcar_indices = []
            for i, v in enumerate(self.natoms):
                potcar_indices += [i] * v
            self.nelects = [self.potcar[potcar_indices[i]].nelectrons for i in range(len(self.structure))] \
                if self.potcar else []

        else:
            fpath = os.path.abspath(cube_filename)
            self.is_vasp = False
            self.cube = Cube(fpath)
            self.structure = self.cube.structure
            self.nelects = self.structure.site_properties.get('nelect', [])  # For cube, see if struc has nelects

        tmpfile = 'CHGCAR' if chgcar_filename else 'CUBE'
        with ScratchDir("."):
            with zopen(fpath, 'rt') as f_in:
                with open(tmpfile, "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            args = [BADEREXE, tmpfile]
            if chgref_filename:
                with zopen(chgrefpath, 'rt') as f_in:
                    with open("CHGCAR_ref", "wt") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                args += ['-ref', 'CHGCAR_ref']
            if parse_atomic_densities:
                args += ['-p', 'all_atom']
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
            except ValueError:
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

            if self.parse_atomic_densities:
                # convert the charge denisty for each atom spit out by Bader into Chgcar objects for easy parsing
                atom_chgcars = [Chgcar.from_file("BvAt{}.dat".format(str(i).zfill(4))) for i in
                                range(1, len(self.chgcar.structure) + 1)]

                atomic_densities = []
                # For each atom in the structure
                for atom, loc, chg in zip(self.chgcar.structure,
                                          self.chgcar.structure.frac_coords,
                                          atom_chgcars):
                    # Find the index of the atom in the charge density atom
                    index = np.round(np.multiply(loc, chg.dim))

                    data = chg.data['total']
                    # Find the shift vector in the array
                    shift = (np.divide(chg.dim, 2) - index).astype(int)

                    # Shift the data so that the atomic charge density to the center for easier manipulation
                    shifted_data = np.roll(data, shift, axis=(0, 1, 2))

                    # Slices a central window from the data array
                    def slice_from_center(data, xwidth, ywidth, zwidth):
                        x, y, z = data.shape
                        startx = x // 2 - (xwidth // 2)
                        starty = y // 2 - (ywidth // 2)
                        startz = z // 2 - (zwidth // 2)
                        return data[startx:startx + xwidth, starty:starty + ywidth, startz:startz + zwidth]

                    # Finds the central encompassing volume which holds all the data within a precision
                    def find_encompassing_vol(data, prec=1e-3):
                        total = np.sum(data)
                        for i in range(np.max(data.shape)):
                            sliced_data = slice_from_center(data, i, i, i)
                            if total - np.sum(sliced_data) < 0.1:
                                return sliced_data
                        return None

                    d = {
                        "data": find_encompassing_vol(shifted_data),
                        "shift": shift,
                        "dim": self.chgcar.dim
                    }
                    atomic_densities.append(d)
                self.atomic_densities = atomic_densities

    def get_charge(self, atom_index):
        """
        Convenience method to get the charge on a particular atom. If the cube file
        is a spin-density file, then this will return the spin density per atom with
        positive being spin up and negative being spin down.

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
        if not self.nelects:
            raise ValueError("No NELECT info! Need POTCAR for VASP, or a structure object"
                             "with nelect as a site property.")
        return self.data[atom_index]["charge"] - self.nelects[atom_index]

    def get_charge_decorated_structure(self):
        """
        Returns an charge decorated structure

        Note, this assumes that the Bader analysis was correctly performed on a file
        with electron densities

        """
        charges = [-self.get_charge(i) for i in range(len(self.structure))]
        struc = self.structure.copy()
        struc.add_site_property('charge', charges)
        return struc

    def get_oxidation_state_decorated_structure(self):
        """
        Returns an oxidation state decorated structure based on bader analysis results.

        Note, this assumes that the Bader analysis was correctly performed on a file
        with electron densities
        """
        charges = [-self.get_charge_transfer(i) for i in range(len(self.structure))]
        struc = self.structure.copy()
        struc.add_oxidation_state_by_site(charges)
        return struc

    def get_spin_state_decorated_structure(self):
        """
        Returns a structure decorated with spins.

        Note, this assumes that the Bader analysis was correctly performed on a file
        with spin densities
        """
        spins = [self.get_charge(i) for i in range(len(self.structure))]
        struc = self.structure.copy()
        struc.add_spin_by_site(spins)
        return struc

    @property
    def summary(self):
        """
        :return: Dict summary of key analysis, e.g., atomic volume, charge, etc.
        """

        summary = {
            "min_dist": [d['min_dist'] for d in self.data],
            "charge": [d['charge'] for d in self.data],
            "atomic_volume": [d['atomic_vol'] for d in self.data],
            "vacuum_charge": self.vacuum_charge,
            "vacuum_volume": self.vacuum_volume,
            "reference_used": self.reference_used,
            "bader_version": self.version,
        }

        if self.parse_atomic_densities:
            summary["charge_densities"] = self.atomic_densities

        if self.potcar:
            charge_transfer = [self.get_charge_transfer(i) for i in range(len(self.data))]
            summary['charge_transfer'] = charge_transfer

        return summary

    @classmethod
    def from_path(cls, path, suffix=""):
        """
        Convenient constructor that takes in the path name of VASP run
        to perform Bader analysis.

        Args:
            path (str): Name of directory where VASP output files are
                stored.
            suffix (str): specific suffix to look for (e.g. '.relax1'
                for 'CHGCAR.relax1.gz').

        """

        def _get_filepath(filename):
            name_pattern = filename + suffix + '*' if filename != 'POTCAR' \
                else filename + '*'
            paths = glob.glob(os.path.join(path, name_pattern))
            fpath = None
            if len(paths) >= 1:
                # using reverse=True because, if multiple files are present,
                # they likely have suffixes 'static', 'relax', 'relax2', etc.
                # and this would give 'static' over 'relax2' over 'relax'
                # however, better to use 'suffix' kwarg to avoid this!
                paths.sort(reverse=True)
                warning_msg = "Multiple files detected, using %s" \
                              % os.path.basename(paths[0]) if len(paths) > 1 \
                    else None
                fpath = paths[0]
            else:
                warning_msg = "Could not find %s" % filename
                if filename in ['AECCAR0', 'AECCAR2']:
                    warning_msg += ", cannot calculate charge transfer."
                elif filename == "POTCAR":
                    warning_msg += ", interpret Bader results with caution."
            if warning_msg:
                warnings.warn(warning_msg)
            return fpath

        chgcar_filename = _get_filepath("CHGCAR")
        if chgcar_filename is None:
            raise IOError("Could not find CHGCAR!")
        potcar_filename = _get_filepath("POTCAR")
        aeccar0 = _get_filepath("AECCAR0")
        aeccar2 = _get_filepath("AECCAR2")
        if (aeccar0 and aeccar2):
            # `chgsum.pl AECCAR0 AECCAR2` equivalent to obtain chgref_file
            chgref = Chgcar.from_file(aeccar0) + Chgcar.from_file(aeccar2)
            chgref_filename = "CHGREF"
            chgref.write_file(chgref_filename)
        else:
            chgref_filename = None
        return cls(chgcar_filename=chgcar_filename, potcar_filename=potcar_filename,
                   chgref_filename=chgref_filename)


def get_filepath(filename, warning, path, suffix):
    """
    Args:
        filename: Filename
        warning: Warning message
        path: Path to search
        suffix: Suffixes to search.
    """
    paths = glob.glob(os.path.join(path, filename + suffix + '*'))
    if not paths:
        warnings.warn(warning)
        return None
    if len(paths) > 1:
        # using reverse=True because, if multiple files are present,
        # they likely have suffixes 'static', 'relax', 'relax2', etc.
        # and this would give 'static' over 'relax2' over 'relax'
        # however, better to use 'suffix' kwarg to avoid this!
        paths.sort(reverse=True)
        warnings.warn('Multiple files detected, using {}'.format(os.path.basename(path)))
    path = paths[0]
    return path


def bader_analysis_from_path(path, suffix=''):
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
    :param suffix: specific suffix to look for (e.g. '.relax1' for 'CHGCAR.relax1.gz'
    :return: summary dict
    """

    def _get_filepath(filename, warning, path=path, suffix=suffix):
        paths = glob.glob(os.path.join(path, filename + suffix + '*'))
        if not paths:
            warnings.warn(warning)
            return None
        if len(paths) > 1:
            # using reverse=True because, if multiple files are present,
            # they likely have suffixes 'static', 'relax', 'relax2', etc.
            # and this would give 'static' over 'relax2' over 'relax'
            # however, better to use 'suffix' kwarg to avoid this!
            paths.sort(reverse=True)
            warnings.warn('Multiple files detected, using {}'.format(os.path.basename(path)))
        path = paths[0]
        return path

    chgcar_path = _get_filepath('CHGCAR', 'Could not find CHGCAR!')
    chgcar = Chgcar.from_file(chgcar_path)

    aeccar0_path = _get_filepath('AECCAR0', 'Could not find AECCAR0, interpret Bader results with caution.')
    aeccar0 = Chgcar.from_file(aeccar0_path) if aeccar0_path else None

    aeccar2_path = _get_filepath('AECCAR2', 'Could not find AECCAR2, interpret Bader results with caution.')
    aeccar2 = Chgcar.from_file(aeccar2_path) if aeccar2_path else None

    potcar_path = _get_filepath('POTCAR', 'Could not find POTCAR, cannot calculate charge transfer.')
    potcar = Potcar.from_file(potcar_path) if potcar_path else None

    return bader_analysis_from_objects(chgcar, potcar, aeccar0, aeccar2)


def bader_analysis_from_objects(chgcar, potcar=None, aeccar0=None, aeccar2=None):
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
    :return: summary dict
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

        ba = BaderAnalysis(chgcar_filename=chgcar_path, potcar_filename=potcar_path, chgref_filename=chgref_path)

        summary = {
            "min_dist": [d['min_dist'] for d in ba.data],
            "charge": [d['charge'] for d in ba.data],
            "atomic_volume": [d['atomic_vol'] for d in ba.data],
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
            ba = BaderAnalysis(chgcar_filename=chgcar_mag_path,
                               potcar_filename=potcar_path,
                               chgref_filename=chgref_path)
            summary["magmom"] = [d['charge'] for d in ba.data]

        return summary
