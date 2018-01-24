# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import warnings
import numpy as np
from pymatgen.core.structure import Specie, Structure
from pymatgen.electronic_structure.core import Magmom
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import AutoOxiStateDecorationTransformation
from pymatgen.analysis.bond_valence import BVAnalyzer
from monty.serialization import loadfn
from enum import Enum, unique
import itertools
import os

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
    warnings.warn("Could not load default_magmoms.yaml, "
                  "falling back to VASPIncarBase.yaml")
    DEFAULT_MAGMOMS = loadfn(os.path.join(MODULE_DIR, "../../io/vasp/VASPIncarBase.yaml"))
    DEFAULT_MAGMOMS = DEFAULT_MAGMOMS['MAGMOM']

@unique
class Ordering(Enum):
    FM = 'FM'  # Ferromagnetic
    AFM = 'AFM'  # Antiferromagnetic
    FiM = 'FiM'  # Ferrimagnetic
    NM = 'NM'  # Non-magnetic
    Unknown = 'Unknown'

class CollinearMagneticStructureAnalyzer:
    def __init__(self, structure,
                 overwrite_magmom_mode="none",
                 round_magmoms=False,
                 detect_valences=False,
                 make_primitive=True,
                 default_magmoms=None,
                 threshold=0.1):
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

        :param structure: Structure object
        :param overwrite_magmom_mode (str): default "none"
        :param round_magmoms (int): will round input magmoms to
        specified number of decimal places, suggest value of 1 or False
        for typical DFT calculations depending on application
        :param detect_valences (bool): if True, will attempt to assign valences
        to input structure
        :param make_primitive (bool): if True, will transform to primitive
        magnetic cell
        :param default_magmoms (dict): (optional) dict specifying default magmoms
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
            raise NotImplementedError("Not implemented for disordered structures, "
                                      "make ordered approximation first.")

        if detect_valences:
            trans = AutoOxiStateDecorationTransformation()
            bva = BVAnalyzer()
            try:
                structure = trans.apply_transformation(structure)
            except ValueError:
                warnings.warn("Could not assign valences "
                              "for {}".format(structure.composition.reduced_formula))

        # check to see if structure has magnetic moments
        # on site properties or species spin properties,
        # prioritize site properties

        has_magmoms = bool(structure.site_properties.get('magmom', False))

        has_spin = False
        for comp in structure.species_and_occu:
            for sp, occu in comp.items():
                if getattr(sp, 'spin', False):
                    has_spin = True

        # perform input sanitation ...
        # rest of class will assume magnetic moments
        # are stored on site properties:
        # this is somewhat arbitrary, arguments can
        # be made for both approaches

        if has_magmoms and has_spin:
            raise ValueError("Structure contains magnetic moments on both "
                             "magmom site properties and spin species "
                             "properties. This is ambiguous. Remove one or "
                             "the other.")
        elif has_magmoms:
            if None in structure.site_properties['magmom']:
                warnings.warn("Be careful with mixing types in your magmom "
                              "site properties. Any 'None' magmoms have been "
                              "replaced with zero.")
            magmoms = [m if m else 0 for m in structure.site_properties['magmom']]
        elif has_spin:
            magmoms = [getattr(sp, 'spin', 0) for sp
                       in structure.species]
            structure.remove_spin()
        else:
            # no magmoms present, add zero magmoms for now
            magmoms = [0]*len(structure)
            # and overwrite magmoms with default magmoms later unless otherwise stated
            if overwrite_magmom_mode == "replace_all_if_undefined":
                overwrite_magmom_mode = "replace_all"

        # test to see if input structure has collinear magmoms
        self.is_collinear = Magmom.are_collinear(magmoms)

        if not self.is_collinear:
            warnings.warn("This class is not designed to be used with "
                          "non-collinear structures. If your structure is "
                          "only slightly non-collinear (e.g. canted) may still "
                          "give useful results, but use with caution.")

        # this is for collinear structures only, make sure magmoms
        # are all floats
        magmoms = list(map(float, magmoms))

        # set properties that should be done /before/ we process input magmoms
        self.total_magmoms = sum(magmoms)
        self.magnetization = sum(magmoms)/structure.volume

        # round magmoms below threshold to zero
        magmoms = [m if abs(m) > threshold else 0 for m in magmoms]

        # overwrite existing magmoms with default_magmoms
        if overwrite_magmom_mode not in ("none", "respect_sign",
                                         "respect_zeros", "replace_all",
                                         "replace_all_if_undefined"):
            raise ValueError("Unsupported mode.")

        for idx, site in enumerate(structure):

            if site.species_string in self.default_magmoms:
                # look for species first, e.g. Fe2+
                default_magmom = self.default_magmoms[site.species_string]
            elif isinstance(site.specie, Specie) and \
                    str(site.specie.element) in self.default_magmoms:
                # look for element, e.g. Fe
                default_magmom = self.default_magmoms[str(site.specie.element)]
            else:
                default_magmom = 0

            # overwrite_magmom_mode = "respect_sign" will change magnitude of
            # existing moments only, and keep zero magmoms as
            # zero: it will keep the magnetic ordering intact

            if overwrite_magmom_mode == "respect_sign":
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

        # round magmoms to specified number of
        # decimal places, used to smooth out
        # computational data
        # TODO: be a bit smarter about rounding magmoms!
        if round_magmoms:
            magmoms = np.around(structure.site_properties['magmom'],
                                decimals=round_magmoms)
            structure.add_site_property(magmoms)

        structure.add_site_property('magmom', magmoms)

        if make_primitive:
            structure = structure.get_primitive_structure(use_site_props=True)

        self.structure = structure

    def get_structure_with_spin(self):
        """
        Returns a Structure with species decorated with spin values instead
        of using magmom site properties.
        :return: Structure
        """

        structure = self.structure.copy()
        structure.add_spin_by_site(structure.site_properties['magmom'])
        structure.remove_site_property('magmom')

        return structure

    def get_structure_with_only_magnetic_atoms(self, make_primitive=True):
        """
        Returns a Structure with only magnetic atoms present.
        :return: Structure
        """

        sites = [site for site in self.structure
                 if abs(site.properties['magmom']) > 0]

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
        structure.remove_site_property('magmom')

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

        structure.add_site_property('magmom',
                                    [abs(m) for m in self.magmoms])

        if make_primitive:
            structure = structure.get_primitive_structure(use_site_props=True)

        return structure

    @property
    def is_magnetic(self):
        """
        Convenience property, returns True if any non-zero magmoms present.
        :return:
        """
        return any(map(abs, self.structure.site_properties['magmom']))

    @property
    def magmoms(self):
        """
        Convenience property, returns magmoms as a numpy array.
        :return: np.array
        """

        return np.array(self.structure.site_properties['magmom'])

    @property
    def types_of_magnetic_specie(self):
        """
        Equivalent to Structure.types_of_specie but only returns
        magnetic species.
        :return: types of Specie
        """
        structure = self.get_structure_with_only_magnetic_atoms()
        return structure.types_of_specie

    @property
    def magnetic_species_and_magmoms(self):
        """
        Returns a dict of magnetic species and the magnitude of
        their associated magmoms. Implicitly assumes the magnetic
        moment is the same magnitude for a given species.
        :return: dict of magnetic species and magmoms
        """

        # TODO: improve detection when magnitude of magmoms varies

        structure = self.get_ferromagnetic_structure()

        magtypes = {str(site.specie): site.properties['magmom'] for site in structure
                    if site.properties['magmom'] > 0}

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

        sga = SpacegroupAnalyzer(structure, symprec=symprec,
                                 angle_tolerance=angle_tolerance)

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
            warnings.warn('Detecting ordering in non-collinear structures not yet implemented.')
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

        return structure.get_space_group_info(symprec=symprec,
                                              angle_tolerance=angle_tolerance)

    def matches_ordering(self, other):
        """
        Compares the magnetic orderings of one structure with another.
        :param other: Structure
        :return (bool):
        """

        a = CollinearMagneticStructureAnalyzer(self.structure,
                                               overwrite_magmom_mode="respect_sign")\
            .get_structure_with_spin()

        # sign of spins doesn't matter, so we're comparing both
        # positive and negative versions of the structure
        # this code is possibly redundant, but is included out of
        # an abundance of caution
        b_positive = CollinearMagneticStructureAnalyzer(other,
                                                        overwrite_magmom_mode="respect_sign")

        b_negative = b_positive.structure.copy()
        b_negative.add_site_property('magmom',
                                     np.multiply(-1, b_negative.site_properties['magmom']))

        b_negative = CollinearMagneticStructureAnalyzer(b_negative,
                                                        overwrite_magmom_mode="respect_sign")

        b_positive = b_positive.get_structure_with_spin()
        b_negative = b_negative.get_structure_with_spin()

        if a.matches(b_positive) or a.matches(b_negative):  # sometimes returns None (bug?)
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
        sorted_indices = np.lexsort((frac_coords[:, 2],
                                     frac_coords[:, 1],
                                     frac_coords[:, 0]))
        s = Structure.from_sites([self.structure[idx] for idx in sorted_indices])

        # adapted from Structure.__repr__
        outs = ["Structure Summary", repr(s.lattice)]
        outs.append("Magmoms Sites")
        for site in s:
            if site.properties['magmom'] != 0:
                prefix = "{:+.2f}   ".format(site.properties['magmom'])
            else:
                prefix = "        "
            outs.append(prefix+repr(site))
        return "\n".join(outs)
