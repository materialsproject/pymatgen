# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines standard transformations which transforms a structure into
another structure. Standard transformations operate in a structure-wide manner,
rather than site-specific manner.
All transformations should inherit the AbstractTransformation ABC.
"""

from __future__ import annotations

import logging
from fractions import Fraction

from numpy import around

from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.analysis.ewald import EwaldMinimizer, EwaldSummation
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.composition import Composition
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.site_transformations import (
    PartialRemoveSitesTransformation,
)
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.util.typing import SpeciesLike

logger = logging.getLogger(__name__)


class RotationTransformation(AbstractTransformation):
    """
    The RotationTransformation applies a rotation to a structure.
    """

    def __init__(self, axis, angle, angle_in_radians=False):
        """
        Args:
            axis (3x1 array): Axis of rotation, e.g., [1, 0, 0]
            angle (float): Angle to rotate
            angle_in_radians (bool): Set to True if angle is supplied in radians.
                Else degrees are assumed.
        """
        self.axis = axis
        self.angle = angle
        self.angle_in_radians = angle_in_radians
        self._symmop = SymmOp.from_axis_angle_and_translation(self.axis, self.angle, self.angle_in_radians)

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Args:
            structure (Structure): Input Structure

        Returns:
            Rotated Structure.
        """
        s = structure.copy()
        s.apply_operation(self._symmop)
        return s

    def __str__(self):
        return (
            f"Rotation Transformation about axis {self.axis} with angle = "
            f"{self.angle:.4f} {'radians' if self.angle_in_radians else 'degrees'}"
        )

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns:
            Inverse Transformation.
        """
        return RotationTransformation(self.axis, -self.angle, self.angle_in_radians)

    @property
    def is_one_to_many(self):
        """Returns: False"""
        return False


class OxidationStateDecorationTransformation(AbstractTransformation):
    """
    This transformation decorates a structure with oxidation states.
    """

    def __init__(self, oxidation_states):
        """
        Args:
            oxidation_states (dict): Oxidation states supplied as a dict,
            e.g., {"Li":1, "O":-2}
        """
        self.oxidation_states = oxidation_states

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Args:
            structure (Structure): Input Structure

        Returns:
            Oxidation state decorated Structure.
        """
        s = structure.copy()
        s.add_oxidation_state_by_element(self.oxidation_states)
        return s

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class AutoOxiStateDecorationTransformation(AbstractTransformation):
    """
    This transformation automatically decorates a structure with oxidation
    states using a bond valence approach.
    """

    def __init__(
        self,
        symm_tol=0.1,
        max_radius=4,
        max_permutations=100000,
        distance_scale_factor=1.015,
    ):
        """
        Args:
            symm_tol (float): Symmetry tolerance used to determine which sites are
                symmetrically equivalent. Set to 0 to turn off symmetry.
            max_radius (float): Maximum radius in Angstrom used to find nearest
                neighbors.
            max_permutations (int): Maximum number of permutations of oxidation
                states to test.
            distance_scale_factor (float): A scale factor to be applied. This is
                useful for scaling distances, esp in the case of
                calculation-relaxed structures, which may tend to under (GGA) or
                over bind (LDA). The default of 1.015 works for GGA. For
                experimental structure, set this to 1.
        """
        self.symm_tol = symm_tol
        self.max_radius = max_radius
        self.max_permutations = max_permutations
        self.distance_scale_factor = distance_scale_factor
        self.analyzer = BVAnalyzer(symm_tol, max_radius, max_permutations, distance_scale_factor)

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Args:
            structure (Structure): Input Structure

        Returns:
            Oxidation state decorated Structure.
        """
        return self.analyzer.get_oxi_state_decorated_structure(structure)

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class OxidationStateRemovalTransformation(AbstractTransformation):
    """
    This transformation removes oxidation states from a structure.
    """

    def __init__(self):
        """
        No arg needed.
        """

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Args:
            structure (Structure): Input Structure

        Returns:
            Non-oxidation state decorated Structure.
        """
        s = structure.copy()
        s.remove_oxidation_states()
        return s

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class SupercellTransformation(AbstractTransformation):
    """
    The RotationTransformation applies a rotation to a structure.
    """

    def __init__(self, scaling_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1))):
        """
        Args:
            scaling_matrix: A matrix of transforming the lattice vectors.
                Defaults to the identity matrix. Has to be all integers. e.g.,
                [[2,1,0],[0,3,0],[0,0,1]] generates a new structure with
                lattice vectors a" = 2a + b, b" = 3b, c" = c where a, b, and c
                are the lattice vectors of the original structure.
        """
        self.scaling_matrix = scaling_matrix

    @staticmethod
    def from_scaling_factors(scale_a=1, scale_b=1, scale_c=1):
        """
        Convenience method to get a SupercellTransformation from a simple
        series of three numbers for scaling each lattice vector. Equivalent to
        calling the normal with [[scale_a, 0, 0], [0, scale_b, 0],
        [0, 0, scale_c]]

        Args:
            scale_a: Scaling factor for lattice direction a. Defaults to 1.
            scale_b: Scaling factor for lattice direction b. Defaults to 1.
            scale_c: Scaling factor for lattice direction c. Defaults to 1.

        Returns:
            SupercellTransformation.
        """
        return SupercellTransformation([[scale_a, 0, 0], [0, scale_b, 0], [0, 0, scale_c]])

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Args:
            structure (Structure): Input Structure

        Returns:
            Supercell Structure.
        """
        return structure * self.scaling_matrix

    def __str__(self):
        return f"Supercell Transformation with scaling matrix {self.scaling_matrix}"

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Raises: NotImplementedError
        """
        raise NotImplementedError()

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class SubstitutionTransformation(AbstractTransformation):
    """
    This transformation substitutes species for one another.
    """

    def __init__(
        self,
        species_map: dict[SpeciesLike, SpeciesLike | dict[SpeciesLike, float]] | list[tuple[SpeciesLike, SpeciesLike]],
    ) -> None:
        """
        Args:
            species_map: A dict or list of tuples containing the species mapping in
                string-string pairs. E.g., {"Li": "Na"} or [("Fe2+","Mn2+")].
                Multiple substitutions can be done. Overloaded to accept
                sp_and_occu dictionary E.g. {"Si: {"Ge":0.75, "C":0.25}},
                which substitutes a single species with multiple species to
                generate a disordered structure.
        """
        self.species_map = species_map
        self._species_map = dict(species_map)
        for k, v in self._species_map.items():
            if isinstance(v, (tuple, list)):
                self._species_map[k] = dict(v)  # type: ignore

    def apply_transformation(self, structure: Structure) -> Structure:
        """
        Apply the transformation.

        Args:
            structure (Structure): Input Structure

        Returns:
            Substituted Structure.
        """
        species_map = {}
        for k, v in self._species_map.items():
            if isinstance(v, dict):
                value = {get_el_sp(x): y for x, y in v.items()}
            else:
                value = get_el_sp(v)  # type: ignore
            species_map[get_el_sp(k)] = value
        s = structure.copy()
        s.replace_species(species_map)
        return s

    def __str__(self):
        return "Substitution Transformation :" + ", ".join(
            [str(k) + "->" + str(v) for k, v in self._species_map.items()]
        )

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns:
            Inverse Transformation.
        """
        inverse_map = {v: k for k, v in self._species_map.items()}
        return SubstitutionTransformation(inverse_map)

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class RemoveSpeciesTransformation(AbstractTransformation):
    """
    Remove all occurrences of some species from a structure.
    """

    def __init__(self, species_to_remove):
        """
        Args:
            species_to_remove: List of species to remove. E.g., ["Li", "Mn"]
        """
        self.species_to_remove = species_to_remove

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Args:
            structure (Structure): Input Structure

        Returns:
            Structure with species removed.
        """
        s = structure.copy()
        for sp in self.species_to_remove:
            s.remove_species([get_el_sp(sp)])
        return s

    def __str__(self):
        return "Remove Species Transformation :" + ", ".join(self.species_to_remove)

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class PartialRemoveSpecieTransformation(AbstractTransformation):
    """
    Remove fraction of specie from a structure.

    Requires an oxidation state decorated structure for Ewald sum to be
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
            specie_to_remove: Species to remove. Must have oxidation state E.g.,
                "Li+"
            fraction_to_remove: Fraction of specie to remove. E.g., 0.5
            algo: This parameter allows you to choose the algorithm to perform
                ordering. Use one of PartialRemoveSpecieTransformation.ALGO_*
                variables to set the algo.
        """
        self.specie_to_remove = specie_to_remove
        self.fraction_to_remove = fraction_to_remove
        self.algo = algo

    def apply_transformation(self, structure: Structure, return_ranked_list=False):
        """
        Apply the transformation.

        Args:
            structure: input structure
            return_ranked_list (bool/int): Boolean stating whether or not
                multiple structures are returned. If return_ranked_list is
                an int, that number of structures is returned.

        Returns:
            Depending on returned_ranked list, either a transformed structure
            or a list of dictionaries, where each dictionary is of the form
            {"structure" = .... , "other_arguments"}
            the key "transformation" is reserved for the transformation that
            was actually applied to the structure.
            This transformation is parsed by the alchemy classes for generating
            a more specific transformation history. Any other information will
            be stored in the transformation_parameters dictionary in the
            transmuted structure class.
        """
        sp = get_el_sp(self.specie_to_remove)
        specie_indices = [i for i in range(len(structure)) if structure[i].species == Composition({sp: 1})]
        trans = PartialRemoveSitesTransformation([specie_indices], [self.fraction_to_remove], algo=self.algo)
        return trans.apply_transformation(structure, return_ranked_list)

    @property
    def is_one_to_many(self):
        """
        Returns: True
        """
        return True

    def __str__(self):
        spec_str = [
            f"Species = {self.specie_to_remove}",
            f"Fraction to remove = {self.fraction_to_remove}",
            f"ALGO = {self.algo}",
        ]
        return "PartialRemoveSpecieTransformation : " + ", ".join(spec_str)

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None


class OrderDisorderedStructureTransformation(AbstractTransformation):
    """
    Order a disordered structure. The disordered structure must be oxidation
    state decorated for Ewald sum to be computed. No attempt is made to perform
    symmetry determination to reduce the number of combinations.

    Hence, attempting to performing ordering on a large number of disordered
    sites may be extremely expensive. The time scales approximately with the
    number of possible combinations. The algorithm can currently compute
    approximately 5,000,000 permutations per minute.

    Also, simple rounding of the occupancies are performed, with no attempt
    made to achieve a target composition. This is usually not a problem for
    most ordering problems, but there can be times where rounding errors may
    result in structures that do not have the desired composition.
    This second step will be implemented in the next iteration of the code.

    If multiple fractions for a single species are found for different sites,
    these will be treated separately if the difference is above a threshold
    tolerance. currently this is .1

    For example, if a fraction of .25 Li is on sites 0,1,2,3  and .5 on sites
    4, 5, 6, 7 then 1 site from [0,1,2,3] will be filled and 2 sites from [4,5,6,7]
    will be filled, even though a lower energy combination might be found by
    putting all lithium in sites [4,5,6,7].

    USE WITH CARE.
    """

    ALGO_FAST = 0
    ALGO_COMPLETE = 1
    ALGO_BEST_FIRST = 2

    def __init__(self, algo=ALGO_FAST, symmetrized_structures=False, no_oxi_states=False):
        """
        Args:
            algo (int): Algorithm to use.
            symmetrized_structures (bool): Whether the input structures are
                instances of SymmetrizedStructure, and that their symmetry
                should be used for the grouping of sites.
            no_oxi_states (bool): Whether to remove oxidation states prior to
                ordering.
        """
        self.algo = algo
        self._all_structures = []
        self.no_oxi_states = no_oxi_states
        self.symmetrized_structures = symmetrized_structures

    def apply_transformation(self, structure: Structure, return_ranked_list=False):
        """
        For this transformation, the apply_transformation method will return
        only the ordered structure with the lowest Ewald energy, to be
        consistent with the method signature of the other transformations.
        However, all structures are stored in the  all_structures attribute in
        the transformation object for easy access.

        Args:
            structure: Oxidation state decorated disordered structure to order
            return_ranked_list (bool): Whether or not multiple structures are
                returned. If return_ranked_list is a number, that number of
                structures is returned.

        Returns:
            Depending on returned_ranked list, either a transformed structure
            or a list of dictionaries, where each dictionary is of the form
            {"structure" = .... , "other_arguments"}
            the key "transformation" is reserved for the transformation that
            was actually applied to the structure.
            This transformation is parsed by the alchemy classes for generating
            a more specific transformation history. Any other information will
            be stored in the transformation_parameters dictionary in the
            transmuted structure class.
        """

        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        num_to_return = max(1, num_to_return)

        if self.no_oxi_states:
            structure = Structure.from_sites(structure)
            for i, site in enumerate(structure):
                structure[i] = {f"{k.symbol}0+": v for k, v in site.species.items()}

        equivalent_sites: list[list[int]] = []
        exemplars: list[PeriodicSite] = []
        # generate list of equivalent sites to order
        # equivalency is determined by sp_and_occu and symmetry
        # if symmetrized structure is true
        for i, site in enumerate(structure):
            if site.is_ordered:
                continue
            for j, ex in enumerate(exemplars):
                sp = ex.species
                if not site.species.almost_equals(sp):
                    continue
                if self.symmetrized_structures:
                    sym_equiv = structure.find_equivalent_sites(ex)
                    sym_test = site in sym_equiv
                else:
                    sym_test = True
                if sym_test:
                    equivalent_sites[j].append(i)
                    break
            else:
                equivalent_sites.append([i])
                exemplars.append(site)

        # generate the list of manipulations and input structure
        s = Structure.from_sites(structure)

        m_list = []
        for g in equivalent_sites:
            total_occupancy = sum((structure[i].species for i in g), Composition())
            total_occupancy = dict(total_occupancy.items())
            # round total occupancy to possible values
            for k, v in total_occupancy.items():
                if abs(v - round(v)) > 0.25:
                    raise ValueError("Occupancy fractions not consistent with size of unit cell")
                total_occupancy[k] = int(round(v))
            # start with an ordered structure
            initial_sp = max(total_occupancy, key=lambda x: abs(x.oxi_state))
            for i in g:
                s[i] = initial_sp
            # determine the manipulations
            for k, v in total_occupancy.items():
                if k == initial_sp:
                    continue
                m = [
                    k.oxi_state / initial_sp.oxi_state if initial_sp.oxi_state else 0,
                    v,
                    list(g),
                    k,
                ]
                m_list.append(m)
            # determine the number of empty sites
            empty = len(g) - sum(total_occupancy.values())
            if empty > 0.5:
                m_list.append([0, empty, list(g), None])

        matrix = EwaldSummation(s).total_energy_matrix
        ewald_m = EwaldMinimizer(matrix, m_list, num_to_return, self.algo)

        self._all_structures = []

        lowest_energy = ewald_m.output_lists[0][0]
        num_atoms = sum(structure.composition.values())

        for output in ewald_m.output_lists:
            s_copy = s.copy()
            # do deletions afterwards because they screw up the indices of the
            # structure
            del_indices = []
            for manipulation in output[1]:
                if manipulation[1] is None:
                    del_indices.append(manipulation[0])
                else:
                    s_copy[manipulation[0]] = manipulation[1]
            s_copy.remove_sites(del_indices)

            if self.no_oxi_states:
                s_copy.remove_oxidation_states()

            self._all_structures.append(
                {
                    "energy": output[0],
                    "energy_above_minimum": (output[0] - lowest_energy) / num_atoms,
                    "structure": s_copy.get_sorted_structure(),
                }
            )

        if return_ranked_list:
            return self._all_structures[:num_to_return]
        return self._all_structures[0]["structure"]

    def __str__(self):
        return "Order disordered structure transformation"

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: True
        """
        return True

    @property
    def lowest_energy_structure(self):
        """
        :return: Lowest energy structure found.
        """
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
            tolerance (float): Tolerance for each coordinate of a particular
                site. For example, [0.5, 0, 0.5] in Cartesian coordinates will be
                considered to be on the same coordinates as [0, 0, 0] for a
                tolerance of 0.5. Defaults to 0.5.
        """
        self.tolerance = tolerance

    def apply_transformation(self, structure):
        """
        Returns most primitive cell for structure.

        Args:
            structure: A structure

        Returns:
            The most primitive structure found. The returned structure is
            guaranteed to have len(new structure) <= len(structure).
        """
        return structure.get_primitive_structure(tolerance=self.tolerance)

    def __str__(self):
        return "Primitive cell transformation"

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class ConventionalCellTransformation(AbstractTransformation):
    """
    This class finds the conventional cell of the input structure.
    """

    def __init__(self, symprec=0.01, angle_tolerance=5, international_monoclinic=True):
        """
        Args:
            symprec (float): tolerance as in SpacegroupAnalyzer
            angle_tolerance (float): angle tolerance as in SpacegroupAnalyzer
            international_monoclinic (bool): whether to use beta (True) or alpha (False)
        as the non-right-angle in the unit cell
        """
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance
        self.international_monoclinic = international_monoclinic

    def apply_transformation(self, structure):
        """
        Returns most primitive cell for structure.

        Args:
            structure: A structure

        Returns:
            The same structure in a conventional standard setting
        """
        sga = SpacegroupAnalyzer(structure, symprec=self.symprec, angle_tolerance=self.angle_tolerance)
        return sga.get_conventional_standard_structure(international_monoclinic=self.international_monoclinic)

    def __str__(self):
        return "Conventional cell transformation"

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class PerturbStructureTransformation(AbstractTransformation):
    """
    This transformation perturbs a structure by a specified distance in random
    directions. Used for breaking symmetries.
    """

    def __init__(
        self,
        distance: float = 0.01,
        min_distance: int | float | None = None,
    ):
        """
        Args:
            distance: Distance of perturbation in angstroms. All sites
                will be perturbed by exactly that distance in a random
                direction.
            min_distance: if None, all displacements will be equidistant. If int
                or float, perturb each site a distance drawn from the uniform
                distribution between 'min_distance' and 'distance'.

        """
        self.distance = distance
        self.min_distance = min_distance

    def apply_transformation(self, structure: Structure) -> Structure:
        """
        Apply the transformation.

        Args:
            structure: Input Structure

        Returns:
            Structure with sites perturbed.
        """
        s = structure.copy()
        s.perturb(self.distance, min_distance=self.min_distance)
        return s

    def __str__(self):
        return f"PerturbStructureTransformation : Min_distance = {self.min_distance}"

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class DeformStructureTransformation(AbstractTransformation):
    """
    This transformation deforms a structure by a deformation gradient matrix
    """

    def __init__(self, deformation=((1, 0, 0), (0, 1, 0), (0, 0, 1))):
        """
        Args:
            deformation (array): deformation gradient for the transformation
        """
        self._deform = Deformation(deformation)
        self.deformation = self._deform.tolist()

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Args:
            structure (Structure): Input Structure

        Returns:
            Deformed Structure.
        """
        return self._deform.apply_to_structure(structure)

    def __str__(self):
        return f"DeformStructureTransformation : Deformation = {self.deformation}"

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns:
            Inverse Transformation.
        """
        return DeformStructureTransformation(self._deform.inv)

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class DiscretizeOccupanciesTransformation(AbstractTransformation):
    """
    Discretizes the site occupancies in a disordered structure; useful for
    grouping similar structures or as a pre-processing step for order-disorder
    transformations.
    """

    def __init__(self, max_denominator=5, tol=None, fix_denominator=False):
        """
        Args:
            max_denominator:
                An integer maximum denominator for discretization. A higher
                denominator allows for finer resolution in the site occupancies.
            tol:
                A float that sets the maximum difference between the original and
                discretized occupancies before throwing an error. If None, it is
                set to 1 / (4 * max_denominator).
            fix_denominator(bool):
                If True, will enforce a common denominator for all species.
                This prevents a mix of denominators (for example, 1/3, 1/4)
                that might require large cell sizes to perform an enumeration.
                'tol' needs to be > 1.0 in some cases.
        """
        self.max_denominator = max_denominator
        self.tol = tol if tol is not None else 1 / (4 * max_denominator)
        self.fix_denominator = fix_denominator

    def apply_transformation(self, structure):
        """
        Discretizes the site occupancies in the structure.

        Args:
            structure: disordered Structure to discretize occupancies

        Returns:
            A new disordered Structure with occupancies discretized
        """
        if structure.is_ordered:
            return structure

        species = [dict(sp) for sp in structure.species_and_occu]

        for sp in species:
            for k in sp:
                old_occ = sp[k]
                new_occ = float(Fraction(old_occ).limit_denominator(self.max_denominator))
                if self.fix_denominator:
                    new_occ = around(old_occ * self.max_denominator) / self.max_denominator
                if round(abs(old_occ - new_occ), 6) > self.tol:
                    raise RuntimeError("Cannot discretize structure within tolerance!")
                sp[k] = new_occ

        return Structure(structure.lattice, species, structure.frac_coords)

    def __str__(self):
        return "DiscretizeOccupanciesTransformation"

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class ChargedCellTransformation(AbstractTransformation):
    """
    The ChargedCellTransformation applies a charge to a structure (or defect
    object).
    """

    def __init__(self, charge=0):
        """
        Args:
            charge: A integer charge to apply to the structure.
                Defaults to zero. Has to be a single integer. e.g. 2
        """
        self.charge = charge

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Args:
            structure (Structure): Input Structure

        Returns:
            Charged Structure.
        """
        s = structure.copy()
        s.set_charge(self.charge)
        return s

    def __str__(self):
        return f"Structure with charge {self.charge}"

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Raises: NotImplementedError
        """
        raise NotImplementedError()

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False


class ScaleToRelaxedTransformation(AbstractTransformation):
    """
    Takes the unrelaxed and relaxed structure and applies its site and volume
    relaxation to a structurally similar structures (e.g. bulk: NaCl and PbTe
    (rock-salt), slab: Sc(10-10) and Mg(10-10) (hcp), GB: Mo(001) sigma 5 GB,
    Fe(001) sigma 5). Useful for finding an initial guess of a set of similar
    structures closer to its most relaxed state.
    """

    def __init__(self, unrelaxed_structure, relaxed_structure, species_map=None):
        """
        Args:
            unrelaxed_structure (Structure): Initial, unrelaxed structure
            relaxed_structure (Structure): Relaxed structure
            species_map (dict): A dict or list of tuples containing the species mapping in
                string-string pairs. The first species corresponds to the relaxed
                structure while the second corresponds to the species in the
                structure to be scaled. E.g., {"Li":"Na"} or [("Fe2+","Mn2+")].
                Multiple substitutions can be done. Overloaded to accept
                sp_and_occu dictionary E.g. {"Si: {"Ge":0.75, "C":0.25}},
                which substitutes a single species with multiple species to
                generate a disordered structure.
        """
        # Get the ratio matrix for lattice relaxation which can be
        # applied to any similar structure to simulate volumetric relaxation
        relax_params = list(relaxed_structure.lattice.abc)
        relax_params.extend(relaxed_structure.lattice.angles)
        unrelax_params = list(unrelaxed_structure.lattice.abc)
        unrelax_params.extend(unrelaxed_structure.lattice.angles)

        self.params_percent_change = []
        for i, p in enumerate(relax_params):
            self.params_percent_change.append(p / unrelax_params[i])

        self.unrelaxed_structure = unrelaxed_structure
        self.relaxed_structure = relaxed_structure
        self.species_map = species_map

    def apply_transformation(self, structure):
        """
        Returns a copy of structure with lattice parameters
        and sites scaled to the same degree as the relaxed_structure.

        Arg:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.
        """

        if self.species_map is None:
            match = StructureMatcher()
            s_map = match.get_best_electronegativity_anonymous_mapping(self.unrelaxed_structure, structure)
        else:
            s_map = self.species_map

        params = list(structure.lattice.abc)
        params.extend(structure.lattice.angles)
        new_lattice = Lattice.from_parameters(*(p * self.params_percent_change[i] for i, p in enumerate(params)))
        species, frac_coords = [], []
        for site in self.relaxed_structure:
            species.append(s_map[site.specie])
            frac_coords.append(site.frac_coords)

        return Structure(new_lattice, species, frac_coords)

    def __str__(self):
        return "ScaleToRelaxedTransformation"

    def __repr__(self):
        return str(self)

    @property
    def inverse(self):
        """
        Returns: None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False
