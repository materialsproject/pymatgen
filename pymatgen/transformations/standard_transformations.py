#!/usr/bin/env python

"""
This module defines standard transformations which transforms a structure into
another structure. Standard transformations operate in a structure-wide manner,
rather than site-specific manner.
All transformations should inherit the AbstractTransformation ABC.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Will Richards"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.2"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"

import itertools
import numpy as np
from operator import itemgetter
import logging

from pymatgen.core.periodic_table import smart_element_or_specie
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure_modifier import StructureEditor, SupercellMaker
from pymatgen.analysis.ewald import EwaldSummation, EwaldMinimizer
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.transformations.site_transformations import \
    PartialRemoveSitesTransformation

logger = logging.getLogger(__name__)


class IdentityTransformation(AbstractTransformation):
    """
    This is a demo transformation which does nothing, i.e. just returns a copy
    of the same structure.
    """

    def __init__(self):
        pass

    def apply_transformation(self, structure):
        return Structure(structure.lattice, structure.species_and_occu,
                         structure.frac_coords)

    def __str__(self):
        return "Identity Transformation"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return self

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "init_args": {},
                "version": __version__, "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class RotationTransformation(AbstractTransformation):
    """
    The RotationTransformation applies a rotation to a structure.
    """

    def __init__(self, axis, angle, angle_in_radians=False):
        """
        Args:
            axis:
                Axis of rotation, e.g., [1, 0, 0]
            angle:
                Angle to rotate
            angle_in_radians:
                Set to True if angle is supplied in radians. Else degrees are
                assumed.
        """
        self._axis = axis
        self._angle = angle
        self._angle_in_radians = angle_in_radians
        self._symmop = SymmOp.from_axis_angle_and_translation(self._axis,
                                        self._angle, self._angle_in_radians)

    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        editor.apply_operation(self._symmop)
        return editor.modified_structure

    def __str__(self):
        return "Rotation Transformation about axis " + \
               "{} with angle = {:.4f} {}".format(self._axis, self._angle,
                            "radians" if self._angle_in_radians else "degrees")

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return RotationTransformation(self._axis, -self._angle,
                                      self._angle_in_radians)

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"axis": self._axis, "angle": self._angle,
                              "angle_in_radians": self._angle_in_radians},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class OxidationStateDecorationTransformation(AbstractTransformation):
    """
    This transformation decorates a structure with oxidation states.
    """

    def __init__(self, oxidation_states):
        """
        Args:
            oxidation_states
                Oxidation states supplied as a dict, e.g., {"Li":1, "O":-2}
        """
        self.oxi_states = oxidation_states

    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        editor.add_oxidation_state_by_element(self.oxi_states)
        return editor.modified_structure

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"oxidation_states": self.oxi_states},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class AutoOxiStateDecorationTransformation(AbstractTransformation):
    """
    This transformation automatically decorates a structure with oxidation
    states using a bond valence approach.
    """

    def __init__(self, symm_tol=0.1, max_radius=4, max_permutations=100000,
                 distance_scale_factor=1.015):
        """
        Args:
            symm_tol:
                Symmetry tolerance used to determine which sites are
                symmetrically equivalent. Set to 0 to turn off symmetry.
            max_radius:
                Maximum radius in Angstrom used to find nearest neighbors.
            max_permutations:
                The maximum number of permutations of oxidation states to test.
            distance_scale_factor:
                A scale factor to be applied. This is useful for scaling
                distances, esp in the case of calculation-relaxed structures
                which may tend to under (GGA) or over bind (LDA). The default
                of 1.015 works for GGA. For experimental structure, set this to
                1.
        """
        self.analyzer = BVAnalyzer(symm_tol, max_radius, max_permutations,
                                   distance_scale_factor)

    def apply_transformation(self, structure):
        return self.analyzer.get_oxi_state_decorated_structure(structure)

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"symm_tol": self.analyzer.symm_tol,
                              "max_radius": self.analyzer.max_radius,
                              "max_permutations": self.analyzer
                              .max_permutations,
                              "distance_scale_factor":
                                  self.analyzer.dist_scale_factor},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class OxidationStateRemovalTransformation(AbstractTransformation):
    """
    This transformation removes oxidation states from a structure
    """
    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        editor.remove_oxidation_states()
        return editor.modified_structure

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {}, "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class SupercellTransformation(AbstractTransformation):
    """
    The RotationTransformation applies a rotation to a structure.
    """

    def __init__(self, scaling_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1))):
        """
        Args:
            scaling_matrix:
                a matrix of transforming the lattice vectors. Defaults to the
                identity matrix. Has to be all integers. e.g.,
                [[2,1,0],[0,3,0],[0,0,1]] generates a new structure with
                lattice vectors a" = 2a + b, b" = 3b, c" = c where a, b, and c
                are the lattice vectors of the original structure.
        """
        self._matrix = scaling_matrix

    @staticmethod
    def from_scaling_factors(scale_a=1, scale_b=1, scale_c=1):
        """
        Convenience method to get a SupercellTransformation from a simple
        series of three numbers for scaling each lattice vector. Equivalent to
        calling the normal with [[scale_a, 0, 0], [0, scale_b, 0],
        [0, 0, scale_c]]

        Args:
            scale_a:
                Scaling factor for lattice direction a. Defaults to 1.
            scale_b:
                Scaling factor for lattice direction b. Defaults to 1.
            scale_c:
                Scaling factor for lattice direction c. Defaults to 1.
        """
        return SupercellTransformation([[scale_a, 0, 0], [0, scale_b, 0],
                                        [0, 0, scale_c]])

    def apply_transformation(self, structure):
        maker = SupercellMaker(structure, self._matrix)
        return maker.modified_structure

    def __str__(self):
        return "Supercell Transformation with scaling matrix " + \
            "{}".format(self._matrix)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        raise NotImplementedError()

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"scaling_matrix": self._matrix},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class SubstitutionTransformation(AbstractTransformation):
    """
    This transformation substitutes species for one another.
    """
    def __init__(self, species_map):
        """
        Args:
            species_map:
                A dict containing the species mapping in string-string pairs.
                E.g., { "Li":"Na"} or {"Fe2+","Mn2+"}. Multiple substitutions
                can be done. Overloaded to accept sp_and_occu dictionary
                E.g. {"Si: {"Ge":0.75, "C":0.25} }, which substitutes a single
                species with multiple species to generate a disordered
                structure.
        """
        self._species_map = species_map

    def apply_transformation(self, structure):
        species_map = {}
        for k, v in self._species_map.items():
            if isinstance(v, dict):
                value = {smart_element_or_specie(x): y for x, y in v.items()}
            else:
                value = smart_element_or_specie(v)
            species_map[smart_element_or_specie(k)] = value
        editor = StructureEditor(structure)
        editor.replace_species(species_map)
        return editor.modified_structure

    def __str__(self):
        return "Substitution Transformation :" + \
            ", ".join([k + "->" + str(v)
                       for k, v in self._species_map.items()])

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        inverse_map = {v: k for k, v in self._species_map.items()}
        return SubstitutionTransformation(inverse_map)

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"species_map": self._species_map},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class RemoveSpeciesTransformation(AbstractTransformation):
    """
    Remove all occurrences of some species from a structure.
    """
    def __init__(self, species_to_remove):
        """
        Args:
            species_to_remove:
                List of species to remove. E.g., ["Li", "Mn"]
        """
        self._species = species_to_remove

    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        map(editor.remove_species, [[smart_element_or_specie(sp)]
                                    for sp in self._species])
        return editor.modified_structure

    def __str__(self):
        return "Remove Species Transformation :" + ", ".join(self._species)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"species_to_remove": self._species},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class PartialRemoveSpecieTransformation(AbstractTransformation):
    """
    Remove fraction of specie from a structure.

    Requires an oxidation state decorated structure for ewald sum to be
    computed.

    Given that the solution to selecting the right removals is NP-hard, there
    are several algorithms provided with varying degrees of accuracy and speed.
    Please see
    :class:`pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation`.
    """

    ALGO_FAST = 0
    ALGO_COMPLETE = 1
    ALGO_BEST_FIRST = 2
    ALGO_ENUMERATE = 3

    def __init__(self, specie_to_remove, fraction_to_remove, algo=ALGO_FAST):
        """
        Args:
            specie_to_remove:
                Specie to remove. Must have oxidation state E.g., "Li1+"
            fraction_to_remove:
                Fraction of specie to remove. E.g., 0.5
            algo:
                This parameter allows you to choose the algorithm to perform
                ordering. Use one of PartialRemoveSpecieTransformation.ALGO_*
                variables to set the algo.
        """
        self._specie = specie_to_remove
        self._frac = fraction_to_remove
        self._algo = algo

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Apply the transformation.

        Args:
            structure:
                input structure
            return_ranked_list:
                Boolean stating whether or not multiple structures are
                returned. If return_ranked_list is an int, that number of
                structures is returned.

        Returns:
            Depending on returned_ranked list, either a transformed structure
            or
            a list of dictionaries, where each dictionary is of the form
            {"structure" = .... , "other_arguments"}
            the key "transformation" is reserved for the transformation that
            was actually applied to the structure.
            This transformation is parsed by the alchemy classes for generating
            a more specific transformation history. Any other information will
            be stored in the transformation_parameters dictionary in the
            transmuted structure class.
        """
        sp = smart_element_or_specie(self._specie)
        specie_indices = [i for i in xrange(len(structure)) \
                          if structure[i].specie == sp]
        trans = PartialRemoveSitesTransformation([specie_indices],
                                                 [self._frac], algo=self._algo)
        return trans.apply_transformation(structure, return_ranked_list)

    @property
    def is_one_to_many(self):
        return True

    def __str__(self):
        spec_str = ["Species = {}".format(self._specie),
                    "Fraction to remove = {}".format(self._frac),
                    "ALGO = {}".format(self._algo)]
        return "PartialRemoveSpecieTransformation : " + ", ".join(spec_str)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"specie_to_remove": self._specie,
                              "fraction_to_remove": self._frac,
                              "algo": self._algo},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class OrderDisorderedStructureTransformation(AbstractTransformation):
    """
    Order a disordered structure. The disordered structure must be oxidation
    state decorated for ewald sum to be computed. No attempt is made to perform
    symmetry determination to reduce the number of combinations.

    Hence, attempting to performing ordering on a large number of disordered
    sites may be extremely expensive. The time scales approximately with the
    number of possible combinations. The algorithm can currently compute
    approximately 5,000,000 permutations per minute.

    Also, simple rounding of the occupancies are performed, with no attempt
    made to achieve a target composition.  This is usually not a problem for
    most ordering problems, but there can be times where rounding errors may
    result in structures that do not have the desired composition.
    This second step will be implemented in the next iteration of the code.

    If multiple fractions for a single species are found for different sites,
    these will be treated separately if the difference is above a threshold
    tolerance. currently this is .1

    For example, if a fraction of .25 Li is on sites 0,1,2,3  and .5 on sites
    4, 5, 6, 7 1 site from [0,1,2,3] will be filled and 2 sites from [4,5,6,7]
    will be filled, even though a lower energy combination might be found by
    putting all lithium in sites [4,5,6,7].

    USE WITH CARE.
    """

    ALGO_FAST = 0
    ALGO_COMPLETE = 1
    ALGO_BEST_FIRST = 2

    def __init__(self, algo=ALGO_FAST):
        """
        Args:
            num_structures:
                maximum number of structures to return
            mev_cutoff:
                maximum mev per atom above the minimum energy ordering for a
                structure to be returned
        """
        self._algo = algo
        self._all_structures = []

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        For this transformation, the apply_transformation method will return
        only the ordered structure with the lowest Ewald energy, to be
        consistent with the method signature of the other transformations.
        However, all structures are stored in the  all_structures attribute in
        the transformation object for easy access.

        Args:
            structure:
                Oxidation state decorated disordered structure to order
            return_ranked_list:
                Boolean stating whether or not multiple structures are
                returned. If return_ranked_list is a number, that number of
                structures is returned.

        Returns:
            Depending on returned_ranked list, either a transformed structure
            or
            a list of dictionaries, where each dictionary is of the form
            {"structure" = .... , "other_arguments"}
            the key "transformation" is reserved for the transformation that
            was actually applied to the structure.
            This transformation is parsed by the alchemy classes for generating
            a more specific transformation history. Any other information will
            be stored in the transformation_parameters dictionary in the
            transmuted structure class.
        """
        ordered_sites = []
        sites_to_order = {}

        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        num_to_return = max(1, num_to_return)

        sites = list(structure.sites)
        for i in range(len(structure)):
            site = sites[i]
            if sum(site.species_and_occu.values()) == 1 and \
                    len(site.species_and_occu) == 1:
                ordered_sites.append(site)
            else:
                species = tuple([sp for sp, occu
                                 in site.species_and_occu.items()])
                #group the sites by the list of species on that site
                for sp, occu in site.species_and_occu.items():
                    if species not in sites_to_order:
                        sites_to_order[species] = {}
                    if sp not in sites_to_order[species]:
                        sites_to_order[species][sp] = [[occu, i]]
                    else:
                        sites_to_order[species][sp].append([occu, i])

                total_occu = sum(site.species_and_occu.values())
                #if the total occupancy on a site is less than one, add
                #a list with None as the species (for removal)
                if total_occu < 1:
                    if None not in sites_to_order[species]:
                        sites_to_order[species][None] = [[1 - total_occu, i]]
                    else:
                        sites_to_order[species][None].append([1 - total_occu,
                                                              i])

        """
        Create a list of [multiplication fraction, number of replacements,
        [indices], replacement species]
        """

        m_list = []
        se = StructureEditor(structure)

        for species in sites_to_order.values():
            initial_sp = None
            sorted_keys = sorted(species.keys(),
                    key=lambda x: x is not None and -abs(x.oxi_state) or 1000)
            for sp in sorted_keys:
                if initial_sp is None:
                    initial_sp = sp
                    for site in species[sp]:
                        se.replace_site(site[1], initial_sp)
                else:
                    if sp is None:
                        oxi = 0
                    else:
                        oxi = float(sp.oxi_state)

                    manipulation = [oxi / initial_sp.oxi_state, 0, [], sp]
                    site_list = species[sp]
                    site_list.sort(key=itemgetter(0))

                    prev_fraction = site_list[0][0]
                    for site in site_list:
                        if site[0] - prev_fraction > .1:
                            """
                            tolerance for creating a new group of sites.
                            if site occupancies are similar, they will be put
                            in a group where the fraction has to be consistent
                            over the whole.
                            """
                            manipulation[1] = int(round(manipulation[1]))
                            m_list.append(manipulation)
                            manipulation = [oxi / initial_sp.oxi_state, 0, [],
                                            sp]
                        prev_fraction = site[0]
                        manipulation[1] += site[0]
                        manipulation[2].append(site[1])
                    #if the # of atoms to remove isn"t within .25 of an integer
                    if abs(manipulation[1] - round(manipulation[1])) > .25:
                        raise ValueError("Occupancy fractions not consistent "
                                        "with size of unit cell")

                    manipulation[1] = int(round(manipulation[1]))
                    m_list.append(manipulation)

        structure = se.modified_structure
        matrix = EwaldSummation(structure).total_energy_matrix
        ewald_m = EwaldMinimizer(matrix, m_list, num_to_return, self._algo)

        self._all_structures = []

        lowest_energy = ewald_m.output_lists[0][0]
        num_atoms = sum(structure.composition.values())

        for output in ewald_m.output_lists:
            se = StructureEditor(structure)
            # do deletions afterwards because they screw up the indices of the
            # structure
            del_indices = []

            for manipulation in output[1]:
                if manipulation[1] is None:
                    del_indices.append(manipulation[0])
                else:
                    se.replace_site(manipulation[0], manipulation[1])
            se.delete_sites(del_indices)
            self._all_structures.append({"energy": output[0],
                "energy_above_minimum": (output[0]
                                         - lowest_energy) / num_atoms,
                "structure": se.modified_structure.get_sorted_structure()})

        if return_ranked_list:
            return self._all_structures
        else:
            return self._all_structures[0]["structure"]

    def __str__(self):
        return "Order disordered structure transformation"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return True

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"algo": self._algo},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @property
    def lowest_energy_structure(self):
        return self._all_structures[0]["structure"]


class PrimitiveCellTransformation(AbstractTransformation):
    """
    This class finds the primitive cell of the input structure.
    It returns a structure that is not necessarily orthogonalized
    Author: Will Richards
    """
    def __init__(self, tolerance=0.5):
        """
        Args:
            tolerance:
                Tolerance for each coordinate of a particular site. For
                example, [0.5, 0, 0.5] in cartesian coordinates will be
                considered to be on the same coordinates as [0, 0, 0] for a
                tolerance of 0.5. Defaults to 0.5.
        """
        self._tolerance = tolerance

    def apply_transformation(self, structure):
        """
        Returns most primitive cell for structure.

        Args:
            structure:
                A structure

        Returns:
            The most primitive structure found. The returned structure is
            guaranteed to have len(new structure) <= len(structure).
        """
        return structure.get_primitive_structure(tolerance=self._tolerance)

    def __str__(self):
        return "Primitive cell transformation"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"tolerance": self._tolerance},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class PerturbStructureTransformation(AbstractTransformation):
    """
    This transformation perturbs a structure by a specified distance in random
    directions. Used for breaking symmetries.
    """

    def __init__(self, amplitude=0.01):
        """
        Args:
            amplitude:
                Amplitude of perturbation in angstroms. All sites will be
                perturbed by exactly that amplitude in a random direction.
        """
        self._amp = amplitude

    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        editor.perturb_structure(self._amp)
        return editor.modified_structure

    def __str__(self):
        return "PerturbStructureTransformation : " + \
            "Amplitude = {}".format(self._amp)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"amplitude": self._amp},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}
