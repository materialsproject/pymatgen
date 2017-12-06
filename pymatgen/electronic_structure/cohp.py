# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import re
import warnings
import numpy as np

from monty.io import zopen
from monty.json import MSONable
from pymatgen.electronic_structure.core import Spin
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.units import Ry_to_eV
from pymatgen.io.lmto import LMTOCopl
from pymatgen.io.lobster import Cohpcar
from pymatgen.util.coord_utils import get_linear_interpolated_value

"""
This module defines classes to represent crystal orbital Hamilton
populations (COHP) and integrated COHP (ICOHP), but can also be used
for crystal orbital overlap populations (COOP).

TODO: Write code for orbital-resolved COHPs/COOPs.
"""

__author__ = "Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Marco Esters"
__email__ = "esters@uoregon.edu"
__date__ = "Nov 30, 2017"


class Cohp(MSONable):
    """
    Basic COHP object.

    Args/attributes:
        are_coops: Indicates whether this object describes COHPs or COOPs.

        efermi: Fermi energy.

        energies: A sequence of energies.

        cohp ({Spin: np.array}): representing the COHP for each spin.

        icohp ({Spin: np.array}): representing the ICOHP for each spin.
    """
    def __init__(self, efermi, energies, cohp, are_coops=False, icohp=None):
        self.are_coops = are_coops
        self.efermi = efermi
        self.energies = np.array(energies)
        self.cohp = cohp
        self.icohp = icohp

    def __add__(self, other):
        """
        Adds two COHP together. Checks that energy scales are the same.
        Otherwise, it raises a ValueError. It also adds ICOHP if present.
        If ICOHP is only present in one object, it displays a warning and
        will not add ICOHP.

        Args:
            other: Another COHP object.

        Returns:
            Sum of the two COHPs as a COHP object.
        """
        if not all(np.equal(self.energies, other.energies)):
            raise ValueError("Energies of both COHP are not compatible.")
        populations = {spin: self.populations[spin] + other.populations[spin]
                       for spin in self.cohp}
        if self.icohp is not None and other.icohp is not None:
            int_pop = {spin: self.icohp[spin] + other.icohp[spin]
                       for spin in self.icohp}
        else:
            if self.icohp is not None or other.icohp is not None:
                warnings.warn("One of the COHP objects does not contain "
                              "ICOHPs. Setting ICOHP to None.")
            int_pop = None
        return Cohp(self.efermi, self.energies, populations, icohp=int_pop)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """
        Returns a string that can be easily plotted (e.g. using gnuplot).
        """
        cohpstring = "COOP" if self.are_coops else "COHP"
        spins = {Spin.up: "Up", Spin.down: "Down"}
        header = ["Energy", cohpstring + "Up"]
        data = [self.energies, self.cohp[Spin.up]]
        if Spin.down in self.cohp:
            header.append(cohpstring + "Down")
            data.append(self.cohp[Spin.down])
        if self.icohp:
            header.append("I" + cohpstring + "Up")
            data.append(self.icohp[Spin.up])
            if Spin.down in self.cohp:
                header.append("I" + cohpstring + "Down")
                data.append(self.icohp[Spin.down])
        formatheader = "#" + " ".join(["{:15s}" for __ in header])
        formatdata = " ".join(["{:.5f}" for __ in header])
        stringarray = [formatheader.format(*header)]
        for i, __ in enumerate(self.energies):
            stringarray.append(formatdata.format(*[d[i] for d in data]))
        return "\n".join(stringarray)

    def as_dict(self):
        """
        Json-serializable dict representation of COHP.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "are_coops": self.are_coops,
             "efermi": self.efermi,
             "energies": self.energies.tolist(),
             "COHP": {str(spin): pops.tolist()
                      for spin, pops in self.cohp.items()}}
        if self.icohp:
            d["ICOHP"] = {str(spin): pops.tolist()
                          for spin, pops in self.icohp.items()}
        return d

    def get_cohp(self, spin=None, integrated=False):
        """
        Returns the COHP or ICOHP for a particular spin.

        Args:
            spin: Spin. Can be parsed as spin object, integer (-1/1)
                or str ("up"/"down")
            integrated: Return COHP (False) or ICOHP (True)

        Returns:
            Returns the CHOP or ICOHP for the input spin. If Spin is
            None and both spins are present, both spins will be returned
            as a dictionary.
        """
        if not integrated:
            populations = self.cohp
        else:
            populations = self.icohp

        if populations is None:
            return None
        elif spin is None:
            return populations
        else:
            if type(spin) == int:
                spin = Spin(spin)
            elif type(spin) == str:
                s = {"up": 1, "down": -1}[spin.lower()]
                spin = Spin(s)
            return {spin: populations[spin]}

    def get_icohp(self, spin=None):
        """
        Convenient alternative to get the ICOHP for a particular spin.
        """
        return self.get_cohp(spin=spin, integrated=True)

    def get_interpolated_value(self, energy, integrated=False):
        """
        Returns the COHP for a particular energy.

        Args:
            energy: Energy to return the COHP value for.
        """
        inter = {}
        for spin in self.cohp:
            if not integrated:
                inter[spin] = get_linear_interpolated_value(self.energies,
                                                            self.cohp[spin],
                                                            energy)
            elif self.icohp is not None:
                inter[spin] = get_linear_interpolated_value(self.energies,
                                                            self.icohp[spin],
                                                            energy)
            else:
                raise ValueError("ICOHP is empty.")
        return inter

    @classmethod
    def from_dict(cls, d):
        """
        Returns a COHP object from a dict representation of the COHP.

        """
        if "ICOHP" in d:
            icohp = {Spin(int(key)): np.array(val)
                     for key, val in d["ICOHP"].items()}
        else:
            icohp = None
        return Cohp(d["efermi"], d["energies"],
                    {Spin(int(key)): np.array(val)
                     for key, val in d["COHP"].items()},
                    icohp=icohp, are_coops=d["are_coops"])


class CompleteCohp(Cohp):
    """
    A wrapper class that defines an average COHP, and individual COHPs.

    Args:
        structure: Structure assosciated with this COHP.

        avg_cohp: The average cohp as a COHP object.

        cohps: A dict of COHP objects for individual bonds of the form
            {label: COHP}

        bonds: A dict containing information on the bonds of the form
            {label: {key: val}}. The key-val pair can be any information
            the user wants to put in, but typically contains the sites,
            the bond length, and the number of bonds. If nothing is
            supplied, it will default to an empty dict.

        are_coops: indicates whether the Cohp objects are COHPs or COOPs.
            Defauls to False for COHPs.


    .. attribute: are_coops

         Indicates whether the object is of COOPs or COHPs.

    .. attribute: efermi

         Fermi energy

    .. attribute: energies

         Sequence of energies

    .. attribute: structure

         Structure associated with the COHPs.

    .. attribute: cohp, icohp

         The average COHP/ICOHP.

    .. attribute: all_cohps

         A dict of COHPs for individual bonds of the form {label: COHP}
    """
    def __init__(self, structure, avg_cohp,
                 cohp_dict, bonds=None, are_coops=False):
        super(CompleteCohp, self).__init__(avg_cohp.efermi,
                                           avg_cohp.energies,
                                           avg_cohp.cohp,
                                           are_coops=are_coops,
                                           icohp=avg_cohp.icohp)
        self.structure = structure
        self.are_coops = are_coops
        self.all_cohps = cohp_dict
        if bonds is None:
            self.bonds = {label: {} for label in self.all_cohps.keys()}
        else:
            self.bonds = bonds

    def __str__(self):
        if self.are_coops:
            return "Complete COOPs for " + str(self.structure)
        else:
            return "Complete COHPs for " + str(self.structure)

    def as_dict(self):
        """
        Json-serializable dict representation of CompleteCohp.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "are_coops": self.are_coops,
             "efermi": self.efermi,
             "structure": self.structure.as_dict(),
             "energies": self.energies.tolist(),
             "COHP": {"average": {str(spin): pops.tolist()
                                  for spin, pops in
                                  self.cohp.items()}}}

        if self.icohp is not None:
            d["ICOHP"] = {"average": {str(spin): pops.tolist()
                                      for spin, pops in
                                      self.icohp.items()}}

        for label in self.all_cohps.keys():
            d["COHP"].update({label: {str(spin): pops.tolist()
                                      for spin, pops in
                                      self.all_cohps[label].cohp.items()}})
            if self.all_cohps[label].icohp is not None:
                if "ICOHP" not in d.keys():
                    d["ICOHP"] = {label: {str(spin): pops.tolist()
                                          for spin, pops in
                                          self.all_cohps[label].icohp.items()}}
                else:
                    d["ICOHP"].update({label: {str(spin): pops.tolist()
                                               for spin, pops in
                                               self.all_cohps[label].icohp.items()}})
        if False in [bond_dict == {} for bond_dict in self.bonds.values()]:
            d["bonds"] = {bond: {"length": self.bonds[bond]["length"],
                                 "sites": [site.as_dict() for site
                                           in self.bonds[bond]["sites"]]}
                                 for bond in self.bonds}
        return d

    def get_cohp(self, label, spin=None, integrated=False):
        """
        Get specific COHP or ICOHP.
        """
        if label.lower() == "average":
            return self.cohp.get_cohp(spin=spin, integrated=integrated)
        else:
            return self.all_cohps[label].get_cohp(spin=spin,
                                                  integrated=integrated)

    def get_icohp(self, label, spin=None):
        """
        Convenient alternative to get a specific ICOHP.
        """
        return self.get_cohp(label, spin=spin, integrated=True)

    @classmethod
    def from_dict(cls, d):
        """
        Returns CompleteCohp object from dict representation.
        """
        cohp_dict = {}
        efermi = d["efermi"]
        energies = d["energies"]
        structure = Structure.from_dict(d["structure"])
        if "bonds" in d.keys():
            bonds = {bond: {"length": d["bonds"][bond]["length"],
                            "sites": tuple(PeriodicSite.from_dict(site)
                                           for site in d["bonds"][bond]["sites"])}
                     for bond in d["bonds"]}
        else:
            bonds = None
        for label in d["COHP"]:
            cohp = {Spin(int(spin)): np.array(d["COHP"][label][spin])
                    for spin in d["COHP"][label]}
            try:
                icohp = {Spin(int(spin)): np.array(d["ICOHP"][label][spin])
                         for spin in d["ICOHP"][label]}
            except KeyError:
                icohp = None
            if label == "average":
                avg_cohp = Cohp(efermi, energies, cohp, icohp=icohp)
            else:
                cohp_dict[label] = Cohp(efermi, energies, cohp, icohp=icohp)
        if "average" not in d["COHP"].keys():
            # calculate average
            cohp = np.array([np.array(c)
                             for c in d["COHP"].values()]).mean(axis=0)
            try:
                icohp = np.array([np.array(c)
                                  for c in d["ICOHP"].values()]).mean(axis=0)
            except KeyError:
                icohp = None
            avg_cohp = Cohp(efermi, energies, cohp, icohp=icohp)

        return CompleteCohp(structure, avg_cohp, cohp_dict, bonds=bonds,
                            are_coops=d["are_coops"])

    @classmethod
    def from_file(cls, fmt, filename=None,
                  structure_file=None, are_coops=False):
        """
        Creates a CompleteCohp object from an output file of a COHP
        calculation. Valid formats are either LMTO (for the Stuttgart
        LMTO-ASA code) or LOBSTER (for the LOBSTER code).

        Args:
            cohp_file: Name of the COHP output file. Defaults to COPL
                for LMTO and COHPCAR.lobster/COOPCAR.lobster for LOBSTER.

            are_coops: Indicates whether the populations are COOPs or
                COHPs. Defaults to False for COHPs.

            fmt: A string for the code that was used to calculate
                the COHPs so that the output file can be handled
                correctly. Can take the values "LMTO" or "LOBSTER".

            structure_file: Name of the file containing the structure.
                If no file name is given, use CTRL for LMTO and POSCAR
                for LOBSTER.

        Returns:
            A CompleteCohp object.
        """
        fmt = fmt.upper()
        if fmt == "LMTO":
            # LMTO COOPs cannot be handled yet.
            are_coops = False
            if structure_file is None:
                structure_file = "CTRL"
            if filename is None:
                filename = "COPL"
            cohp_file = LMTOCopl(filename=filename, to_eV=True)
        elif fmt == "LOBSTER":
            if structure_file is None:
                structure_file = "POSCAR"
            if filename is None:
                filename = "COOPCAR.lobster" if are_coops \
                            else "COHPCAR.lobster"
            cohp_file = Cohpcar(filename=filename, are_coops=are_coops)
        else:
            raise ValueError("Unknown format %s. Valid formats are LMTO "
                             "and LOBSTER." % fmt)

        structure = Structure.from_file(structure_file)
        efermi = cohp_file.efermi
        cohp_data = cohp_file.cohp_data
        energies = cohp_file.energies

        # Lobster shifts the energies so that the Fermi energy is at zero.
        # Shifting should be done by the plotter object though.

        if fmt == "LOBSTER":
            energies += efermi

        if fmt == "LMTO":
            # Calculate the average COHP for the LMTO file to be
            # consistent with LOBSTER output.
            spins = [Spin.up, Spin.down] if cohp_file.is_spin_polarized \
                      else [Spin.up]
            avg_data = {"COHP": {}, "ICOHP": {}}
            for i in avg_data:
                for spin in spins:
                        rows = np.array([cohp_data[label][i][spin]
                                         for label in cohp_data])
                        avg = np.average(rows, axis=0)
                        # LMTO COHPs have 5 significant figures
                        avg_data[i].update({spin: np.array(["{:.5g}".format(a)
                                            for a in avg], dtype=float)})
            avg_cohp = Cohp(efermi, energies,
                            avg_data["COHP"],
                            icohp=avg_data["ICOHP"])
        else:
            avg_cohp = Cohp(efermi, energies,
                            cohp_data["average"]["COHP"],
                            icohp=cohp_data["average"]["COHP"],
                            are_coops=are_coops)
            del cohp_data["average"]

        cohp_dict = {label: Cohp(efermi, energies,
                                 cohp_data[label]["COHP"],
                                 icohp=cohp_data[label]["ICOHP"],
                                 are_coops=are_coops)
                     for label in cohp_data}

        bond_dict = {label: {"length": cohp_data[label]["length"],
                             "sites": [structure.sites[site]
                                       for site in cohp_data[label]["sites"]]}
                     for label in cohp_data}

        return CompleteCohp(structure, avg_cohp, cohp_dict, bonds=bond_dict,
                            are_coops=are_coops)
