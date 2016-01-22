# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import numpy as np

import six

from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.structure import Structure
from pymatgen.util.coord_utils import get_linear_interpolated_value
from monty.json import MSONable

"""
This module defines classes to represent the density of states, etc.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 20, 2012"


class Dos(MSONable):
    """
    Basic DOS object. All other DOS objects are extended versions of this
    object.

    Args:
        efermi: Fermi level energy
        energies: A sequences of energies
        densities ({Spin: np.array}): representing the density of states
            for each Spin.

    .. attribute: energies

        The sequence of energies

    .. attribute: densities

        A dict of spin densities, e.g., {Spin.up: [...], Spin.down: [...]}

    .. attribute: efermi

        Fermi level
    """

    def __init__(self, efermi, energies, densities):
        self.efermi = efermi
        self.energies = np.array(energies)
        self.densities = {k: np.array(d) for k, d in densities.items()}

    def get_densities(self, spin=None):
        """
        Returns the density of states for a particular spin.

        Args:
            spin: Spin

        Returns:
            Returns the density of states for a particular spin. If Spin is
            None, the sum of all spins is returned.
        """
        if self.densities is None:
            result = None
        elif spin is None:
            if Spin.down in self.densities:
                result = self.densities[Spin.up] + self.densities[Spin.down]
            else:
                result = self.densities[Spin.up]
        else:
            result = self.densities[spin]
        return result

    def get_smeared_densities(self, sigma):
        """
        Returns the Dict representation of the densities, {Spin: densities},
        but with a Gaussian smearing of std dev sigma applied about the fermi
        level.

        Args:
            sigma: Std dev of Gaussian smearing function.

        Returns:
            Dict of Gaussian-smeared densities.
        """
        from scipy.ndimage.filters import gaussian_filter1d
        smeared_dens = {}
        diff = [self.energies[i + 1] - self.energies[i]
                for i in range(len(self.energies) - 1)]
        avgdiff = sum(diff) / len(diff)
        for spin, dens in self.densities.items():
            smeared_dens[spin] = gaussian_filter1d(dens, sigma / avgdiff)
        return smeared_dens

    def __add__(self, other):
        """
        Adds two DOS together. Checks that energy scales are the same.
        Otherwise, a ValueError is thrown.

        Args:
            other: Another DOS object.

        Returns:
            Sum of the two DOSs.
        """
        if not all(np.equal(self.energies, other.energies)):
            raise ValueError("Energies of both DOS are not compatible!")
        densities = {spin: self.densities[spin] + other.densities[spin]
                     for spin in self.densities.keys()}
        return Dos(self.efermi, self.energies, densities)

    def get_interpolated_value(self, energy):
        """
        Returns interpolated density for a particular energy.

        Args:
            energy: Energy to return the density for.
        """
        f = {}
        for spin in self.densities.keys():
            f[spin] = get_linear_interpolated_value(self.energies,
                                                    self.densities[spin],
                                                    energy)
        return f

    def get_interpolated_gap(self, tol=0.001, abs_tol=False, spin=None):
        """
        Expects a DOS object and finds the gap

        Args:
            tol: tolerance in occupations for determining the gap
            abs_tol: Set to True for an absolute tolerance and False for a
                relative one.
            spin: Possible values are None - finds the gap in the summed
                densities, Up - finds the gap in the up spin channel,
                Down - finds the gap in the down spin channel.

        Returns:
            (gap, cbm, vbm):
                Tuple of floats in eV corresponding to the gap, cbm and vbm.
        """

        tdos = self.get_densities(spin)
        if not abs_tol:
            tol = tol * tdos.sum() / tdos.shape[0]
        energies = self.energies
        below_fermi = [i for i in range(len(energies))
                       if energies[i] < self.efermi and tdos[i] > tol]
        above_fermi = [i for i in range(len(energies))
                       if energies[i] > self.efermi and tdos[i] > tol]
        vbm_start = max(below_fermi)
        cbm_start = min(above_fermi)
        if vbm_start == cbm_start:
            return 0.0, self.efermi, self.efermi
        else:
            # Interpolate between adjacent values
            terminal_dens = tdos[vbm_start:vbm_start + 2][::-1]
            terminal_energies = energies[vbm_start:vbm_start + 2][::-1]
            start = get_linear_interpolated_value(terminal_dens,
                                                  terminal_energies, tol)
            terminal_dens = tdos[cbm_start - 1:cbm_start + 1]
            terminal_energies = energies[cbm_start - 1:cbm_start + 1]
            end = get_linear_interpolated_value(terminal_dens,
                                                terminal_energies, tol)
            return end - start, end, start

    def get_cbm_vbm(self, tol=0.001, abs_tol=False, spin=None):
        """
        Expects a DOS object and finds the cbm and vbm.

        Args:
            tol: tolerance in occupations for determining the gap
            abs_tol: An absolute tolerance (True) and a relative one (False)
            spin: Possible values are None - finds the gap in the summed
                densities, Up - finds the gap in the up spin channel,
                Down - finds the gap in the down spin channel.

        Returns:
            (cbm, vbm): float in eV corresponding to the gap
        """
        #determine tolerance
        tdos = self.get_densities(spin)
        if not abs_tol:
            tol = tol * tdos.sum() / tdos.shape[0]

        # find index of fermi energy
        i_fermi = 0
        while self.energies[i_fermi] <= self.efermi:
            i_fermi += 1

        # work backwards until tolerance is reached
        i_gap_start = i_fermi
        while i_gap_start - 1 >= 0 and tdos[i_gap_start - 1] <= tol:
            i_gap_start -= 1

        # work forwards until tolerance is reached
        i_gap_end = i_gap_start
        while i_gap_end < tdos.shape[0] and tdos[i_gap_end] <= tol:
            i_gap_end += 1
        i_gap_end -= 1
        return self.energies[i_gap_end], self.energies[i_gap_start]

    def get_gap(self, tol=0.001, abs_tol=False, spin=None):
        """
        Expects a DOS object and finds the gap.

        Args:
            tol: tolerance in occupations for determining the gap
            abs_tol: An absolute tolerance (True) and a relative one (False)
            spin: Possible values are None - finds the gap in the summed
                densities, Up - finds the gap in the up spin channel,
                Down - finds the gap in the down spin channel.

        Returns:
            gap in eV
        """
        (cbm, vbm) = self.get_cbm_vbm(tol, abs_tol, spin)
        return max(cbm - vbm, 0.0)

    def __str__(self):
        """
        Returns a string which can be easily plotted (using gnuplot).
        """
        if Spin.down in self.densities:
            stringarray = ["#{:30s} {:30s} {:30s}".format("Energy",
                                                          "DensityUp",
                                                          "DensityDown")]
            for i, energy in enumerate(self.energies):
                stringarray.append("{:.5f} {:.5f} {:.5f}"
                                   .format(energy, self.densities[Spin.up][i],
                                           self.densities[Spin.down][i]))
        else:
            stringarray = ["#{:30s} {:30s}".format("Energy", "DensityUp")]
            for i, energy in enumerate(self.energies):
                stringarray.append("{:.5f} {:.5f}"
                                   .format(energy, self.densities[Spin.up][i]))
        return "\n".join(stringarray)

    @classmethod
    def from_dict(cls, d):
        """
        Returns Dos object from dict representation of Dos.
        """
        return Dos(d["efermi"], d["energies"],
                   {Spin(int(k)): v
                    for k, v in d["densities"].items()})

    def as_dict(self):
        """
        Json-serializable dict representation of Dos.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__, "efermi": self.efermi,
                "energies": list(self.energies),
                "densities": {str(spin): list(dens)
                              for spin, dens in self.densities.items()}}


class CompleteDos(Dos):
    """
    This wrapper class defines a total dos, and also provides a list of PDos.
    Mainly used by pymatgen.io.vasp.Vasprun to create a complete Dos from
    a vasprun.xml file. You are unlikely to try to generate this object
    manually.

    Args:
        structure: Structure associated with this particular DOS.
        total_dos: total Dos for structure
        pdoss: The pdoss are supplied as an {Site:{Orbital:{
            Spin:Densities}}}

    .. attribute:: structure

        Structure associated with the CompleteDos.

    .. attribute:: pdos

        Dict of partial densities of the form {Site:{Orbital:{Spin:Densities}}}
    """
    def __init__(self, structure, total_dos, pdoss):
        super(CompleteDos, self).__init__(
            total_dos.efermi, energies=total_dos.energies,
            densities={k: np.array(d) for k, d in total_dos.densities.items()})
        self.pdos = pdoss
        self.structure = structure

    def get_site_orbital_dos(self, site, orbital):
        """
        Get the Dos for a particular orbital of a particular site.

        Args:
            site: Site in Structure associated with CompleteDos.
            orbital: Orbital in the site.

        Returns:
            Dos containing densities for orbital of site.
        """
        return Dos(self.efermi, self.energies, self.pdos[site][orbital])

    def get_site_dos(self, site):
        """
        Get the total Dos for a site (all orbitals).

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            Dos containing summed orbital densities for site.
        """
        site_dos = six.moves.reduce(add_densities, self.pdos[site].values())
        return Dos(self.efermi, self.energies, site_dos)

    def get_site_spd_dos(self, site):
        """
        Get orbital projected Dos of a particular site

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            dict of {orbital: Dos}, e.g. {"s": Dos object, ...}
        """
        spd_dos = dict()
        for orb, pdos in self.pdos[site].items():
            orbital_type = _get_orb_type(orb)
            if orbital_type in spd_dos:
                spd_dos[orbital_type] = add_densities(spd_dos[orbital_type], pdos)
            else:
                spd_dos[orbital_type] = pdos
        return {orb: Dos(self.efermi, self.energies, densities)
                for orb, densities in spd_dos.items()}

    def get_site_t2g_eg_resolved_dos(self, site):
        """
        Get the t2g, eg projected DOS for a particular site.

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            A dict {"e_g": Dos, "t2g": Dos} containing summed e_g and t2g DOS
            for the site.
        """
        t2g_dos = []
        eg_dos = []
        for s, atom_dos in self.pdos.items():
            if s == site:
                for orb, pdos in atom_dos.items():
                    if orb in (Orbital.dxy, Orbital.dxz, Orbital.dyz):
                        t2g_dos.append(pdos)
                    elif orb in (Orbital.dx2, Orbital.dz2):
                        eg_dos.append(pdos)
        return {"t2g": Dos(self.efermi, self.energies,
                           six.moves.reduce(add_densities, t2g_dos)),
                "e_g": Dos(self.efermi, self.energies,
                           six.moves.reduce(add_densities, eg_dos))}

    def get_spd_dos(self):
        """
        Get orbital projected Dos.

        Returns:
            dict of {orbital: Dos}, e.g. {"s": Dos object, ...}
        """
        spd_dos = {}
        for atom_dos in self.pdos.values():
            for orb, pdos in atom_dos.items():
                orbital_type = _get_orb_type(orb)
                if orbital_type not in spd_dos:
                    spd_dos[orbital_type] = pdos
                else:
                    spd_dos[orbital_type] = \
                        add_densities(spd_dos[orbital_type], pdos)
        return {orb: Dos(self.efermi, self.energies, densities)
                for orb, densities in spd_dos.items()}

    def get_element_dos(self):
        """
        Get element projected Dos.

        Returns:
            dict of {Element: Dos}
        """

        el_dos = {}
        for site, atom_dos in self.pdos.items():
            el = site.specie
            for pdos in atom_dos.values():
                if el not in el_dos:
                    el_dos[el] = pdos
                else:
                    el_dos[el] = add_densities(el_dos[el], pdos)
        return {el: Dos(self.efermi, self.energies, densities)
                for el, densities in el_dos.items()}

    def get_element_spd_dos(self, el):
        """
        Get element and spd projected Dos

        Args:
            el: Element in Structure.composition associated with CompleteDos

        Returns:
            dict of {Element: {"S": densities, "P": densities, "D": densities}}
        """
        el = get_el_sp(el)
        el_dos = {}
        for site, atom_dos in self.pdos.items():
            if site.specie == el:
                for orb, pdos in atom_dos.items():
                    orbital_type = _get_orb_type(orb)
                    if orbital_type not in el_dos:
                        el_dos[orbital_type] = pdos
                    else:
                        el_dos[orbital_type] = \
                            add_densities(el_dos[orbital_type], pdos)

        return {orb: Dos(self.efermi, self.energies, densities)
                for orb, densities in el_dos.items()}

    @classmethod
    def from_dict(cls, d):
        """
        Returns CompleteDos object from dict representation.
        """
        tdos = Dos.from_dict(d)
        struct = Structure.from_dict(d["structure"])
        pdoss = {}
        for i in range(len(d["pdos"])):
            at = struct[i]
            orb_dos = {}
            for orb_str, odos in d["pdos"][i].items():
                orb = Orbital[orb_str]
                orb_dos[orb] = {Spin(int(k)): v
                                for k, v in odos["densities"].items()}
            pdoss[at] = orb_dos
        return CompleteDos(struct, tdos, pdoss)

    def as_dict(self):
        """
        Json-serializable dict representation of CompleteDos.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__, "efermi": self.efermi,
             "structure": self.structure.as_dict(),
             "energies": list(self.energies),
             "densities": {str(spin): list(dens)
                           for spin, dens in self.densities.items()},
             "pdos": []}
        if len(self.pdos) > 0:
            for at in self.structure:
                dd = {}
                for orb, pdos in self.pdos[at].items():
                    dd[str(orb)] = {"densities": {str(int(spin)): list(dens)
                                                  for spin,
                                                  dens in pdos.items()}}
                d["pdos"].append(dd)
            d["atom_dos"] = {str(at): dos.as_dict() for at,
                             dos in self.get_element_dos().items()}
            d["spd_dos"] = {str(orb): dos.as_dict() for orb,
                            dos in self.get_spd_dos().items()}
        return d

    def __str__(self):
        return "Complete DOS for " + str(self.structure)


def add_densities(density1, density2):
    """
    Method to sum two densities.

    Args:
        density1: First density.
        density2: Second density.

    Returns:
        Dict of {spin: density}.
    """
    return {spin: np.array(density1[spin]) + np.array(density2[spin])
            for spin in density1.keys()}


def _get_orb_type(orb):
    try:
        return orb.orbital_type
    except AttributeError:
        return orb