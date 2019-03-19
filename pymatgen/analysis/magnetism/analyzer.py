# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import warnings
import numpy as np
import os

from enum import Enum, unique
from collections import namedtuple

from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema

from pymatgen.core.structure import Specie, Structure
from pymatgen.electronic_structure.core import Magmom
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import (
    AutoOxiStateDecorationTransformation,
)
from pymatgen.analysis.bond_valence import BVAnalyzer
from monty.serialization import loadfn

from typing import Union, List

"""
This module provides some useful functions for dealing with magnetic Structures
(e.g. Structures with associated magmom tags).
"""

__author__ = "Matthew Horton"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Development"
__date__ = "Feb 2017"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

try:
    DEFAULT_MAGMOMS = loadfn(os.path.join(MODULE_DIR, "default_magmoms.yaml"))
except:
    warnings.warn(
        "Could not load default_magmoms.yaml, " "falling back to VASPIncarBase.yaml"
    )
    DEFAULT_MAGMOMS = loadfn(
        os.path.join(MODULE_DIR, "../../io/vasp/VASPIncarBase.yaml")
    )
    DEFAULT_MAGMOMS = DEFAULT_MAGMOMS["MAGMOM"]


@unique
class Ordering(Enum):
    FM = "FM"  # Ferromagnetic
    AFM = "AFM"  # Antiferromagnetic
    FiM = "FiM"  # Ferrimagnetic
    NM = "NM"  # Non-magnetic
    Unknown = "Unknown"


@unique
class MagneticOrderingEnumerationStrategy(Enum):
    pass


@unique
class OverwriteMagmomMode(Enum):
    none = "none"
    respect_sign = "respect_sign"
    respect_zero = "respect_zeros"
    replace_all = "replace_all"
    normalize = "normalize"


class CollinearMagneticStructureAnalyzer:
    def __init__(
        self,
        structure: Structure,
        overwrite_magmom_mode: Union[OverwriteMagmomMode, str] = "none",
        round_magmoms: bool = False,
        detect_valences: bool = False,
        make_primitive: bool = True,
        default_magmoms: bool = None,
        set_net_positive: bool = True,
        threshold: float = 0.1,
    ):
        """
        A class which provides a few helpful methods to analyze
        collinear magnetic structures.

        If magnetic moments are not defined, moments will be
        taken either from default_magmoms.yaml (similar to the
        default magmoms in MPRelaxSet, with a few extra definitions)
        or from a specie:magmom dict provided by the default_magmoms
        kwarg.

        Input magmoms can be replaced using the 'overwrite_magmom_mode'
        kwarg. This can be:
        * "none" to do nothing,
        * "respect_sign" which will overwrite existing magmoms with
          those from default_magmoms but will keep sites with positive magmoms
          positive, negative magmoms negative and zero magmoms zero,
        * "respect_zeros", which will give a ferromagnetic structure
          (all positive magmoms from default_magmoms) but still keep sites with
          zero magmoms as zero,
        * "replace_all" which will try to guess initial magmoms for
          all sites in the structure irrespective of input structure
          (this is most suitable for an initial DFT calculation),
        * "replace_all_if_undefined" is the same as "replace_all" but only if
          no magmoms are defined in input structure, otherwise it will respect
          existing magmoms.
        * "normalize" will normalize magmoms to unity, but will respect sign
          (used for comparing orderings), magmoms < theshold will be set to zero

        :param structure: Structure object
        :param overwrite_magmom_mode (str): default "none"
        :param round_magmoms (int or bool): will round input magmoms to
        specified number of decimal places if integer is supplied, if set
        to a float will try and group magmoms together using a kernel density
        estimator of provided width, and extracting peaks of the estimator
        :param detect_valences (bool): if True, will attempt to assign valences
        to input structure
        :param make_primitive (bool): if True, will transform to primitive
        magnetic cell
        :param default_magmoms (dict): (optional) dict specifying default magmoms
        :param set_net_positive (bool): if True, will change sign of magnetic
        moments such that the net magnetization is positive. Argument will be
        ignored if mode "respect_sign" is used.
        :param threshold (float): number (in Bohr magnetons) below which magmoms
        will be rounded to zero, default of 0.1 can probably be increased for many
        magnetic systems, depending on your application
        """

        if default_magmoms:
            self.default_magmoms = default_magmoms
        else:
            self.default_magmoms = DEFAULT_MAGMOMS

        structure = structure.copy()

        # check for disorder
        if not structure.is_ordered:
            raise NotImplementedError(
                "Not implemented for disordered structures, "
                "make ordered approximation first."
            )

        if detect_valences:
            trans = AutoOxiStateDecorationTransformation()
            bva = BVAnalyzer()
            try:
                structure = trans.apply_transformation(structure)
            except ValueError:
                warnings.warn(
                    "Could not assign valences "
                    "for {}".format(structure.composition.reduced_formula)
                )

        # check to see if structure has magnetic moments
        # on site properties or species spin properties,
        # prioritize site properties

        has_magmoms = bool(structure.site_properties.get("magmom", False))

        has_spin = False
        for comp in structure.species_and_occu:
            for sp, occu in comp.items():
                if getattr(sp, "spin", False):
                    has_spin = True

        # perform input sanitation ...
        # rest of class will assume magnetic moments
        # are stored on site properties:
        # this is somewhat arbitrary, arguments can
        # be made for both approaches

        if has_magmoms and has_spin:
            raise ValueError(
                "Structure contains magnetic moments on both "
                "magmom site properties and spin species "
                "properties. This is ambiguous. Remove one or "
                "the other."
            )
        elif has_magmoms:
            if None in structure.site_properties["magmom"]:
                warnings.warn(
                    "Be careful with mixing types in your magmom "
                    "site properties. Any 'None' magmoms have been "
                    "replaced with zero."
                )
            magmoms = [m if m else 0 for m in structure.site_properties["magmom"]]
        elif has_spin:
            magmoms = [getattr(sp, "spin", 0) for sp in structure.species]
            structure.remove_spin()
        else:
            # no magmoms present, add zero magmoms for now
            magmoms = [0] * len(structure)
            # and overwrite magmoms with default magmoms later unless otherwise stated
            if overwrite_magmom_mode == "replace_all_if_undefined":
                overwrite_magmom_mode = "replace_all"

        # test to see if input structure has collinear magmoms
        self.is_collinear = Magmom.are_collinear(magmoms)

        if not self.is_collinear:
            warnings.warn(
                "This class is not designed to be used with "
                "non-collinear structures. If your structure is "
                "only slightly non-collinear (e.g. canted) may still "
                "give useful results, but use with caution."
            )

        # this is for collinear structures only, make sure magmoms
        # are all floats
        magmoms = list(map(float, magmoms))

        # set properties that should be done /before/ we process input magmoms
        self.total_magmoms = sum(magmoms)
        self.magnetization = sum(magmoms) / structure.volume

        # round magmoms below threshold to zero
        magmoms = [m if abs(m) > threshold else 0 for m in magmoms]

        # overwrite existing magmoms with default_magmoms
        if overwrite_magmom_mode not in (
            "none",
            "respect_sign",
            "respect_zeros",
            "replace_all",
            "replace_all_if_undefined",
            "normalize",
        ):
            raise ValueError("Unsupported mode.")

        for idx, site in enumerate(structure):

            if site.species_string in self.default_magmoms:
                # look for species first, e.g. Fe2+
                default_magmom = self.default_magmoms[site.species_string]
            elif (
                isinstance(site.specie, Specie)
                and str(site.specie.element) in self.default_magmoms
            ):
                # look for element, e.g. Fe
                default_magmom = self.default_magmoms[str(site.specie.element)]
            else:
                default_magmom = 0

            # overwrite_magmom_mode = "respect_sign" will change magnitude of
            # existing moments only, and keep zero magmoms as
            # zero: it will keep the magnetic ordering intact

            if overwrite_magmom_mode == "respect_sign":
                set_net_positive = False
                if magmoms[idx] > 0:
                    magmoms[idx] = default_magmom
                elif magmoms[idx] < 0:
                    magmoms[idx] = -default_magmom

            # overwrite_magmom_mode = "respect_zeros" will give a ferromagnetic
            # structure but will keep zero magmoms as zero

            elif overwrite_magmom_mode == "respect_zeros":
                if magmoms[idx] != 0:
                    magmoms[idx] = default_magmom

            # overwrite_magmom_mode = "replace_all" will ignore input magmoms
            # and give a ferromagnetic structure with magnetic
            # moments on *all* atoms it thinks could be magnetic

            elif overwrite_magmom_mode == "replace_all":
                magmoms[idx] = default_magmom

            # overwrite_magmom_mode = "normalize" set magmoms magnitude to 1

            elif overwrite_magmom_mode == "normalize":
                if magmoms[idx] != 0:
                    magmoms[idx] = int(magmoms[idx] / abs(magmoms[idx]))

        # round magmoms, used to smooth out computational data
        magmoms = (
            self._round_magmoms(magmoms, round_magmoms) if round_magmoms else magmoms
        )

        if set_net_positive:
            sign = np.sum(magmoms)
            if sign < 0:
                magmoms = -np.array(magmoms)

        structure.add_site_property("magmom", magmoms)

        if make_primitive:
            structure = structure.get_primitive_structure(use_site_props=True)

        self.structure = structure

    @staticmethod
    def _round_magmoms(magmoms, round_magmoms_mode: Union[int, float]):
        """
        If round_magmoms_mode is an integer, simply round to that number
        of decimal places, else if set to a float will try and round
        intelligently by grouping magmoms.
        """

        if isinstance(round_magmoms_mode, int):

            # simple rounding to number of decimal places
            magmoms = np.around(magmoms, decimals=round_magmoms_mode)

        elif isinstance(round_magmoms_mode, float):

            try:

                # get range of possible magmoms, pad by 50% just to be safe
                range_m = max([max(magmoms), abs(min(magmoms))]) * 1.5

                # construct kde, here "round_magmoms_mode" is the width of the kde
                kernel = gaussian_kde(magmoms, bw_method=round_magmoms_mode)

                # with a linearly spaced grid 1000x finer than width
                xgrid = np.linspace(
                    -range_m, range_m, 1000 * range_m / round_magmoms_mode
                )

                # and evaluate the kde on this grid, extracting the maxima of the kde peaks
                kernel_m = kernel.evaluate(xgrid)
                extrema = xgrid[argrelextrema(kernel_m, comparator=np.greater)]

                # round magmoms to these extrema
                magmoms = [extrema[(np.abs(extrema - m)).argmin()] for m in magmoms]

            except Exception as e:

                # TODO: typically a singular matrix warning, investigate this
                warnings.warn(
                    "Failed to round magmoms intelligently, "
                    "falling back to simple rounding."
                )
                warnings.warn(e)

            # and finally round roughly to the number of significant figures in our kde width
            num_decimals = len(str(round_magmoms_mode).split(".")[1]) + 1
            magmoms = np.around(magmoms, decimals=num_decimals)

        return magmoms

    def get_structure_with_spin(self):
        """
        Returns a Structure with species decorated with spin values instead
        of using magmom site properties.
        :return: Structure
        """

        structure = self.structure.copy()
        structure.add_spin_by_site(structure.site_properties["magmom"])
        structure.remove_site_property("magmom")

        return structure

    def get_structure_with_only_magnetic_atoms(self, make_primitive=True):
        """
        Returns a Structure with only magnetic atoms present.
        :return: Structure
        """

        sites = [site for site in self.structure if abs(site.properties["magmom"]) > 0]

        structure = Structure.from_sites(sites)

        if make_primitive:
            structure = structure.get_primitive_structure(use_site_props=True)

        return structure

    def get_nonmagnetic_structure(self, make_primitive=True):
        """
        Returns a Structure without magnetic moments defined.
        :param make_primitive (bool): Return a primitive
        structure, defaults to True.
        :return: Structure
        """

        structure = self.structure.copy()
        structure.remove_site_property("magmom")

        if make_primitive:
            structure = structure.get_primitive_structure()

        return structure

    def get_ferromagnetic_structure(self, make_primitive=True):
        """
        Returns a Structure with all magnetic moments positive
        or zero.
        :param make_primitive (bool): Return a primitive
        structure, defaults to True.
        :return: Structure
        """

        structure = self.structure.copy()

        structure.add_site_property("magmom", [abs(m) for m in self.magmoms])

        if make_primitive:
            structure = structure.get_primitive_structure(use_site_props=True)

        return structure

    @property
    def is_magnetic(self):
        """
        Convenience property, returns True if any non-zero magmoms present.
        :return:
        """
        return any(map(abs, self.structure.site_properties["magmom"]))

    @property
    def magmoms(self):
        """
        Convenience property, returns magmoms as a numpy array.
        :return: np.array
        """

        return np.array(self.structure.site_properties["magmom"])

    @property
    def types_of_magnetic_specie(self):
        """
        Equivalent to Structure.types_of_specie but only returns
        magnetic species.
        :return: types of Specie
        """
        if self.number_of_magnetic_sites > 0:
            structure = self.get_structure_with_only_magnetic_atoms()
            return structure.types_of_specie
        else:
            return []

    @property
    def magnetic_species_and_magmoms(self):
        """
        Returns a dict of magnetic species and the magnitude of
        their associated magmoms. Will return a set if there are
        multiple magmoms per species.

        :return: dict of magnetic species and magmoms
        """

        structure = self.get_ferromagnetic_structure()

        magtypes = {
            str(site.specie): set()
            for site in structure
            if site.properties["magmom"] != 0
        }

        for site in structure:
            if site.properties["magmom"] != 0:
                magtypes[str(site.specie)].add(site.properties["magmom"])

        for sp, magmoms in magtypes.items():
            if len(magmoms) == 1:
                magtypes[sp] = magmoms.pop()
            else:
                magtypes[sp] = sorted(list(magmoms))

        return magtypes

    @property
    def number_of_magnetic_sites(self):
        """
        :return (int): Number of magnetic sites present in structure.
        """
        return np.sum([abs(m) > 0 for m in self.magmoms])

    def number_of_unique_magnetic_sites(self, symprec=1e-3, angle_tolerance=5):
        """
        :param symprec (float): same as in SpacegroupAnalyzer
        :param angle_tolerance (float): same as in SpacegroupAnalyzer
        :return (int): Number of symmetrically-distinct magnetic sites present
        in structure.
        """

        structure = self.get_nonmagnetic_structure()

        sga = SpacegroupAnalyzer(
            structure, symprec=symprec, angle_tolerance=angle_tolerance
        )

        symm_structure = sga.get_symmetrized_structure()

        num_unique_mag_sites = 0

        for group_of_sites in symm_structure.equivalent_sites:
            if group_of_sites[0].specie in self.types_of_magnetic_specie:
                num_unique_mag_sites += 1

        return num_unique_mag_sites

    @property
    def ordering(self):
        """
        Applies heuristics to return a magnetic ordering for a collinear
        magnetic structure. Result is not guaranteed for correctness.
        :return: Ordering Enum ('FiM' is used as the abbreviation for
        ferrimagnetic)
        """

        if not self.is_collinear:
            warnings.warn(
                "Detecting ordering in non-collinear structures not yet implemented."
            )
            return Ordering.Unknown

        if "magmom" not in self.structure.site_properties:
            # maybe this was a non-spin-polarized calculation, or we've
            # lost the magnetic moment information
            return Ordering.Unknown

        magmoms = self.magmoms

        max_magmom = max(magmoms)

        total_magnetization = abs(sum(magmoms))

        is_potentially_ferromagnetic = np.all(magmoms >= 0) or np.all(magmoms <= 0)

        if total_magnetization > 0 and is_potentially_ferromagnetic:
            return Ordering.FM
        elif total_magnetization > 0:
            return Ordering.FiM
        elif max_magmom > 0:
            return Ordering.AFM
        else:
            return Ordering.NM

    def get_exchange_group_info(self, symprec=1e-2, angle_tolerance=5.0):
        """
        Returns the information on the symmetry of the Hamiltonian
        describing the exchange energy of the system, taking into
        account relative direction of magnetic moments but not their
        absolute direction.

        This is not strictly accurate (e.g. some/many atoms will
        have zero magnetic moments), but defining symmetry this
        way is a useful way of keeping track of distinct magnetic
        orderings within pymatgen.

        :param symprec: same as SpacegroupAnalyzer
        :param angle_tolerance: same as SpacegroupAnalyzer
        :return: spacegroup_symbol, international_number
        """

        structure = self.get_structure_with_spin()

        return structure.get_space_group_info(
            symprec=symprec, angle_tolerance=angle_tolerance
        )

    def matches_ordering(self, other):
        """
        Compares the magnetic orderings of one structure with another.
        :param other: Structure
        :return (bool):
        """

        a = CollinearMagneticStructureAnalyzer(
            self.structure, overwrite_magmom_mode="normalize"
        ).get_structure_with_spin()

        # sign of spins doesn't matter, so we're comparing both
        # positive and negative versions of the structure
        # this code is possibly redundant, but is included out of
        # an abundance of caution
        b_positive = CollinearMagneticStructureAnalyzer(
            other, overwrite_magmom_mode="normalize", make_primitive=False
        )

        b_negative = b_positive.structure.copy()
        b_negative.add_site_property(
            "magmom", np.multiply(-1, b_negative.site_properties["magmom"])
        )

        b_negative = CollinearMagneticStructureAnalyzer(
            b_negative, overwrite_magmom_mode="normalize", make_primitive=False
        )

        b_positive = b_positive.get_structure_with_spin()
        b_negative = b_negative.get_structure_with_spin()

        if a.matches(b_positive) or a.matches(
            b_negative
        ):  # sometimes returns None (bug?)
            return True
        else:
            return False

    @property
    def propagation_vector(self):
        return NotImplementedError

    def __str__(self):
        """
        Sorts a Structure (by fractional co-ordinate), and
        prints sites with magnetic information. This is
        useful over Structure.__str__ because sites are in
        a consistent order, which makes visual comparison between
        two identical Structures with different magnetic orderings
        easier.
        :return:
        """

        frac_coords = self.structure.frac_coords
        sorted_indices = np.lexsort(
            (frac_coords[:, 2], frac_coords[:, 1], frac_coords[:, 0])
        )
        s = Structure.from_sites([self.structure[idx] for idx in sorted_indices])

        # adapted from Structure.__repr__
        outs = ["Structure Summary", repr(s.lattice)]
        outs.append("Magmoms Sites")
        for site in s:
            if site.properties["magmom"] != 0:
                prefix = "{:+.2f}   ".format(site.properties["magmom"])
            else:
                prefix = "        "
            outs.append(prefix + repr(site))
        return "\n".join(outs)


def magnetic_deformation(structure_A, structure_B):
    """
    Calculates 'magnetic deformation proxy',
    a measure of deformation (norm of finite strain)
    between 'non-magnetic' (non-spin-polarized) and
    ferromagnetic structures.

    Adapted from Bocarsly et al. 2017,
    doi: 10.1021/acs.chemmater.6b04729

    :param structure_A: Structure
    :param structure_B: Structure
    :return:
    """

    # retrieve orderings of both input structures
    ordering_a = CollinearMagneticStructureAnalyzer(
        structure_A, overwrite_magmom_mode="none"
    ).ordering
    ordering_b = CollinearMagneticStructureAnalyzer(
        structure_B, overwrite_magmom_mode="none"
    ).ordering

    # get a type string, this is either 'NM-FM' for between non-magnetic
    # and ferromagnetic, as in Bocarsly paper, or e.g. 'FM-AFM'
    type_str = "{}-{}".format(ordering_a.value, ordering_b.value)

    lattice_a = structure_A.lattice.matrix.T
    lattice_b = structure_B.lattice.matrix.T
    lattice_a_inv = np.linalg.inv(lattice_a)
    p = np.dot(lattice_a_inv, lattice_b)
    eta = 0.5 * (np.dot(p.T, p) - np.identity(3))
    w, v = np.linalg.eig(eta)
    deformation = 100 * (1.0 / 3.0) * np.sqrt(w[0] ** 2 + w[1] ** 2 + w[2] ** 2)

    MagneticDeformation = namedtuple("MagneticDeformation", "type deformation")

    return MagneticDeformation(deformation=deformation, type=type_str)
