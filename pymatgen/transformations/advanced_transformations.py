# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements more advanced transformations.
"""
import logging
import math
import warnings
from fractions import Fraction
from itertools import groupby, product
from math import gcd
from string import ascii_lowercase
from typing import Dict, Optional

import numpy as np
from monty.dev import requires
from monty.fractions import lcm
from monty.json import MSONable

from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.energy_models import SymmetryModel
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.analysis.gb.grain import GrainBoundaryGenerator
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.analysis.structure_matcher import SpinComparator, StructureMatcher
from pymatgen.analysis.structure_prediction.substitution_probability import (
    SubstitutionPredictor,
)
from pymatgen.command_line.enumlib_caller import EnumError, EnumlibAdaptor
from pymatgen.command_line.mcsqs_caller import run_mcsqs
from pymatgen.core.periodic_table import DummySpecies, Element, Species, get_el_sp
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import (
    OrderDisorderedStructureTransformation,
    SubstitutionTransformation,
    SupercellTransformation,
)
from pymatgen.transformations.transformation_abc import AbstractTransformation

try:
    import hiphive  # type: ignore
except ImportError:
    hiphive = None

__author__ = "Shyue Ping Ong, Stephen Dacek, Anubhav Jain, Matthew Horton, Alex Ganose"

logger = logging.getLogger(__name__)


class ChargeBalanceTransformation(AbstractTransformation):
    """
    This is a transformation that disorders a structure to make it charge
    balanced, given an oxidation state-decorated structure.
    """

    def __init__(self, charge_balance_sp):
        """
        Args:
            charge_balance_sp: specie to add or remove. Currently only removal
                is supported
        """
        self.charge_balance_sp = str(charge_balance_sp)

    def apply_transformation(self, structure):
        """
        Applies the transformation.

        Args:
            structure: Input Structure

        Returns:
            Charge balanced structure.
        """
        charge = structure.charge
        specie = get_el_sp(self.charge_balance_sp)
        num_to_remove = charge / specie.oxi_state
        num_in_structure = structure.composition[specie]
        removal_fraction = num_to_remove / num_in_structure
        if removal_fraction < 0:
            raise ValueError("addition of specie not yet supported by ChargeBalanceTransformation")
        trans = SubstitutionTransformation({self.charge_balance_sp: {self.charge_balance_sp: 1 - removal_fraction}})
        return trans.apply_transformation(structure)

    def __str__(self):
        return f"Charge Balance Transformation : Species to remove = {self.charge_balance_sp}"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: False"""
        return False


class SuperTransformation(AbstractTransformation):
    """
    This is a transformation that is inherently one-to-many. It is constructed
    from a list of transformations and returns one structure for each
    transformation. The primary use for this class is extending a transmuter
    object.
    """

    def __init__(self, transformations, nstructures_per_trans=1):
        """
        Args:
            transformations ([transformations]): List of transformations to apply
                to a structure. One transformation is applied to each output
                structure.
            nstructures_per_trans (int): If the transformations are one-to-many and,
                nstructures_per_trans structures from each transformation are
                added to the full list. Defaults to 1, i.e., only best structure.
        """

        self._transformations = transformations
        self.nstructures_per_trans = nstructures_per_trans

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Applies the transformation.

        Args:
            structure: Input Structure
            return_ranked_list: Number of structures to return.

        Returns:
            Structures with all transformations applied.
        """
        if not return_ranked_list:
            raise ValueError("SuperTransformation has no single best structure output. Must use return_ranked_list")
        structures = []
        for t in self._transformations:
            if t.is_one_to_many:
                for d in t.apply_transformation(structure, return_ranked_list=self.nstructures_per_trans):
                    d["transformation"] = t
                    structures.append(d)
            else:
                structures.append(
                    {
                        "transformation": t,
                        "structure": t.apply_transformation(structure),
                    }
                )
        return structures

    def __str__(self):
        return "Super Transformation : Transformations = " + "{}".format(
            " ".join([str(t) for t in self._transformations])
        )

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True


class MultipleSubstitutionTransformation:
    """
    Performs multiple substitutions on a structure. For example, can do a
    fractional replacement of Ge in LiGePS with a list of species, creating one
    structure for each substitution. Ordering is done using a dummy element so
    only one ordering must be done per substitution oxidation state. Charge
    balancing of the structure is optionally performed.

    .. note::
        There are no checks to make sure that removal fractions are possible
        and rounding may occur. Currently charge balancing only works for
        removal of species.
    """

    def __init__(
        self,
        sp_to_replace,
        r_fraction,
        substitution_dict,
        charge_balance_species=None,
        order=True,
    ):
        """
        Performs multiple fractional substitutions on a transmuter.

        Args:
            sp_to_replace: species to be replaced
            r_fraction: fraction of that specie to replace
            substitution_dict: dictionary of the format
                {2: ["Mg", "Ti", "V", "As", "Cr", "Ta", "N", "Nb"],
                3: ["Ru", "Fe", "Co", "Ce", "As", "Cr", "Ta", "N", "Nb"],
                4: ["Ru", "V", "Cr", "Ta", "N", "Nb"],
                5: ["Ru", "W", "Mn"]
                }
                The number is the charge used for each of the list of elements
                (an element can be present in multiple lists)
            charge_balance_species: If specified, will balance the charge on
                the structure using that specie.
        """
        self.sp_to_replace = sp_to_replace
        self.r_fraction = r_fraction
        self.substitution_dict = substitution_dict
        self.charge_balance_species = charge_balance_species
        self.order = order

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Applies the transformation.

        Args:
            structure: Input Structure
            return_ranked_list: Number of structures to return.

        Returns:
            Structures with all substitutions applied.
        """
        if not return_ranked_list:
            raise ValueError(
                "MultipleSubstitutionTransformation has no single"
                " best structure output. Must use"
                " return_ranked_list."
            )
        outputs = []
        for charge, el_list in self.substitution_dict.items():
            mapping = {}
            if charge > 0:
                sign = "+"
            else:
                sign = "-"
            dummy_sp = f"X{charge}{sign}"
            mapping[self.sp_to_replace] = {
                self.sp_to_replace: 1 - self.r_fraction,
                dummy_sp: self.r_fraction,
            }
            trans = SubstitutionTransformation(mapping)
            dummy_structure = trans.apply_transformation(structure)
            if self.charge_balance_species is not None:
                cbt = ChargeBalanceTransformation(self.charge_balance_species)
                dummy_structure = cbt.apply_transformation(dummy_structure)
            if self.order:
                trans = OrderDisorderedStructureTransformation()
                dummy_structure = trans.apply_transformation(dummy_structure)

            for el in el_list:
                if charge > 0:
                    sign = "+"
                else:
                    sign = "-"
                st = SubstitutionTransformation({f"X{charge}+": f"{el}{charge}{sign}"})
                new_structure = st.apply_transformation(dummy_structure)
                outputs.append({"structure": new_structure})
        return outputs

    def __str__(self):
        return "Multiple Substitution Transformation : Substitution on " + f"{self.sp_to_replace}"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True


class EnumerateStructureTransformation(AbstractTransformation):
    """
    Order a disordered structure using enumlib. For complete orderings, this
    generally produces fewer structures that the OrderDisorderedStructure
    transformation, and at a much faster speed.
    """

    def __init__(
        self,
        min_cell_size=1,
        max_cell_size=1,
        symm_prec=0.1,
        refine_structure=False,
        enum_precision_parameter=0.001,
        check_ordered_symmetry=True,
        max_disordered_sites=None,
        sort_criteria="ewald",
        timeout=None,
    ):
        """
        Args:
            min_cell_size:
                The minimum cell size wanted. Must be an int. Defaults to 1.
            max_cell_size:
                The maximum cell size wanted. Must be an int. Defaults to 1.
            symm_prec:
                Tolerance to use for symmetry.
            refine_structure:
                This parameter has the same meaning as in enumlib_caller.
                If you are starting from a structure that has been relaxed via
                some electronic structure code, it is usually much better to
                start with symmetry determination and then obtain a refined
                structure. The refined structure have cell parameters and
                atomic positions shifted to the expected symmetry positions,
                which makes it much less sensitive precision issues in enumlib.
                If you are already starting from an experimental cif, refinement
                should have already been done and it is not necessary. Defaults
                to False.
            enum_precision_parameter (float): Finite precision parameter for
                enumlib. Default of 0.001 is usually ok, but you might need to
                tweak it for certain cells.
            check_ordered_symmetry (bool): Whether to check the symmetry of
                the ordered sites. If the symmetry of the ordered sites is
                lower, the lowest symmetry ordered sites is included in the
                enumeration. This is important if the ordered sites break
                symmetry in a way that is important getting possible
                structures. But sometimes including ordered sites
                slows down enumeration to the point that it cannot be
                completed. Switch to False in those cases. Defaults to True.
            max_disordered_sites (int):
                An alternate parameter to max_cell size. Will sequentially try
                larger and larger cell sizes until (i) getting a result or (ii)
                the number of disordered sites in the cell exceeds
                max_disordered_sites. Must set max_cell_size to None when using
                this parameter.
            sort_criteria (str): Sort by Ewald energy ("ewald", must have oxidation
                states and slow) or by number of sites ("nsites", much faster).
            timeout (float): timeout in minutes to pass to EnumlibAdaptor
        """
        self.symm_prec = symm_prec
        self.min_cell_size = min_cell_size
        self.max_cell_size = max_cell_size
        self.refine_structure = refine_structure
        self.enum_precision_parameter = enum_precision_parameter
        self.check_ordered_symmetry = check_ordered_symmetry
        self.max_disordered_sites = max_disordered_sites
        self.sort_criteria = sort_criteria
        self.timeout = timeout

        if max_cell_size and max_disordered_sites:
            raise ValueError("Cannot set both max_cell_size and max_disordered_sites!")

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Returns either a single ordered structure or a sequence of all ordered
        structures.

        Args:
            structure: Structure to order.
            return_ranked_list (bool): Whether or not multiple structures are
                returned. If return_ranked_list is a number, that number of
                structures is returned.

        Returns:
            Depending on returned_ranked list, either a transformed structure
            or a list of dictionaries, where each dictionary is of the form
            {"structure" = .... , "other_arguments"}

            The list of ordered structures is ranked by ewald energy / atom, if
            the input structure is an oxidation state decorated structure.
            Otherwise, it is ranked by number of sites, with smallest number of
            sites first.
        """
        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        if self.refine_structure:
            finder = SpacegroupAnalyzer(structure, self.symm_prec)
            structure = finder.get_refined_structure()

        contains_oxidation_state = all(
            hasattr(sp, "oxi_state") and sp.oxi_state != 0 for sp in structure.composition.elements
        )

        structures = None

        if structure.is_ordered:
            warnings.warn(
                "Enumeration skipped for structure with composition {} "
                "because it is ordered".format(structure.composition)
            )
            structures = [structure.copy()]

        if self.max_disordered_sites:
            ndisordered = sum(1 for site in structure if not site.is_ordered)
            if ndisordered > self.max_disordered_sites:
                raise ValueError(f"Too many disordered sites! ({ndisordered} > {self.max_disordered_sites})")
            max_cell_sizes = range(
                self.min_cell_size,
                int(math.floor(self.max_disordered_sites / ndisordered)) + 1,
            )
        else:
            max_cell_sizes = [self.max_cell_size]

        for max_cell_size in max_cell_sizes:
            adaptor = EnumlibAdaptor(
                structure,
                min_cell_size=self.min_cell_size,
                max_cell_size=max_cell_size,
                symm_prec=self.symm_prec,
                refine_structure=False,
                enum_precision_parameter=self.enum_precision_parameter,
                check_ordered_symmetry=self.check_ordered_symmetry,
                timeout=self.timeout,
            )
            try:
                adaptor.run()
                structures = adaptor.structures
                if structures:
                    break
            except EnumError:
                warnings.warn(f"Unable to enumerate for max_cell_size = {max_cell_size}")

        if structures is None:
            raise ValueError("Unable to enumerate")

        original_latt = structure.lattice
        inv_latt = np.linalg.inv(original_latt.matrix)
        ewald_matrices = {}
        all_structures = []
        for s in structures:
            new_latt = s.lattice
            transformation = np.dot(new_latt.matrix, inv_latt)
            transformation = tuple(tuple(int(round(cell)) for cell in row) for row in transformation)
            if contains_oxidation_state and self.sort_criteria == "ewald":
                if transformation not in ewald_matrices:
                    s_supercell = structure * transformation
                    ewald = EwaldSummation(s_supercell)
                    ewald_matrices[transformation] = ewald
                else:
                    ewald = ewald_matrices[transformation]
                energy = ewald.compute_sub_structure(s)
                all_structures.append({"num_sites": len(s), "energy": energy, "structure": s})
            else:
                all_structures.append({"num_sites": len(s), "structure": s})

        def sort_func(s):
            return (
                s["energy"] / s["num_sites"]
                if contains_oxidation_state and self.sort_criteria == "ewald"
                else s["num_sites"]
            )

        self._all_structures = sorted(all_structures, key=sort_func)

        if return_ranked_list:
            return self._all_structures[0:num_to_return]
        return self._all_structures[0]["structure"]

    def __str__(self):
        return "EnumerateStructureTransformation"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True


class SubstitutionPredictorTransformation(AbstractTransformation):
    """
    This transformation takes a structure and uses the structure
    prediction module to find likely site substitutions.
    """

    def __init__(self, threshold=1e-2, scale_volumes=True, **kwargs):
        r"""
        Args:
            threshold: Threshold for substitution.
            scale_volumes: Whether to scale volumes after substitution.
            **kwargs: Args for SubstitutionProbability class lambda_table, alpha
        """
        self.kwargs = kwargs
        self.threshold = threshold
        self.scale_volumes = scale_volumes
        self._substitutor = SubstitutionPredictor(threshold=threshold, **kwargs)

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Applies the transformation.

        Args:
            structure: Input Structure
            return_ranked_list: Number of structures to return.

        Returns:
            Predicted Structures.
        """
        if not return_ranked_list:
            raise ValueError("SubstitutionPredictorTransformation doesn't support returning 1 structure")

        preds = self._substitutor.composition_prediction(structure.composition, to_this_composition=False)
        preds.sort(key=lambda x: x["probability"], reverse=True)

        outputs = []
        for pred in preds:
            st = SubstitutionTransformation(pred["substitutions"])
            output = {
                "structure": st.apply_transformation(structure),
                "probability": pred["probability"],
                "threshold": self.threshold,
                "substitutions": {},
            }

            # dictionary keys have to be converted to strings for JSON
            for key, value in pred["substitutions"].items():
                output["substitutions"][str(key)] = str(value)
            outputs.append(output)
        return outputs

    def __str__(self):
        return "SubstitutionPredictorTransformation"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True


class MagOrderParameterConstraint(MSONable):
    """
    This class can be used to supply MagOrderingTransformation
    to just a specific subset of species or sites that satisfy the
    provided constraints. This can be useful for setting an order
    parameters for, for example, ferrimagnetic structures which
    might order on certain motifs, with the global order parameter
    dependent on how many sites satisfy that motif.
    """

    def __init__(
        self,
        order_parameter,
        species_constraints=None,
        site_constraint_name=None,
        site_constraints=None,
    ):
        """
        :param order_parameter (float): any number from 0.0 to 1.0,
            typically 0.5 (antiferromagnetic) or 1.0 (ferromagnetic)
        :param species_constraint (list): str or list of strings
            of Species symbols that the constraint should apply to
        :param site_constraint_name (str): name of the site property
            that the constraint should apply to, e.g. "coordination_no"
        :param site_constraints (list): list of values of the site
            property that the constraints should apply to
        """

        # validation
        if site_constraints and site_constraints != [None] and not site_constraint_name:
            raise ValueError("Specify the name of the site constraint.")
        if not site_constraints and site_constraint_name:
            raise ValueError("Please specify some site constraints.")
        if not isinstance(species_constraints, list):
            species_constraints = [species_constraints]
        if not isinstance(site_constraints, list):
            site_constraints = [site_constraints]

        if order_parameter > 1 or order_parameter < 0:
            raise ValueError("Order parameter must lie between 0 and 1")
        if order_parameter != 0.5:
            warnings.warn(
                "Use care when using a non-standard order parameter, "
                "though it can be useful in some cases it can also "
                "lead to unintended behavior. Consult documentation."
            )

        self.order_parameter = order_parameter
        self.species_constraints = species_constraints
        self.site_constraint_name = site_constraint_name
        self.site_constraints = site_constraints

    def satisfies_constraint(self, site):
        """
        Checks if a periodic site satisfies the constraint.
        """
        if not site.is_ordered:
            return False

        satisfies_constraints = self.species_constraints and str(site.specie) in self.species_constraints

        if self.site_constraint_name and self.site_constraint_name in site.properties:
            prop = site.properties[self.site_constraint_name]
            satisfies_constraints = prop in self.site_constraints

        return satisfies_constraints


class MagOrderingTransformation(AbstractTransformation):
    """
    This transformation takes a structure and returns a list of collinear
    magnetic orderings. For disordered structures, make an ordered
    approximation first.
    """

    def __init__(self, mag_species_spin, order_parameter=0.5, energy_model=SymmetryModel(), **kwargs):
        """
        :param mag_species_spin: A mapping of elements/species to their
            spin magnitudes, e.g. {"Fe3+": 5, "Mn3+": 4}
        :param order_parameter (float or list): if float, a specifies a
            global order parameter and can take values from 0.0 to 1.0
            (e.g. 0.5 for antiferromagnetic or 1.0 for ferromagnetic), if
            list has to be a list of
            :class: `pymatgen.transformations.advanced_transformations.MagOrderParameterConstraint`
            to specify more complicated orderings, see documentation for
            MagOrderParameterConstraint more details on usage
        :param energy_model: Energy model to rank the returned structures,
            see :mod: `pymatgen.analysis.energy_models` for more information (note
            that this is not necessarily a physical energy). By default, returned
            structures use SymmetryModel() which ranks structures from most
            symmetric to least.
        :param kwargs: Additional kwargs that are passed to
        :class:`EnumerateStructureTransformation` such as min_cell_size etc.
        """

        # checking for sensible order_parameter values
        if isinstance(order_parameter, float):
            # convert to constraint format
            order_parameter = [
                MagOrderParameterConstraint(
                    order_parameter=order_parameter,
                    species_constraints=list(mag_species_spin.keys()),
                )
            ]
        elif isinstance(order_parameter, list):
            ops = [isinstance(item, MagOrderParameterConstraint) for item in order_parameter]
            if not any(ops):
                raise ValueError("Order parameter not correctly defined.")
        else:
            raise ValueError("Order parameter not correctly defined.")

        self.mag_species_spin = mag_species_spin
        # store order parameter constraints as dicts to save implementing
        # to/from dict methods for MSONable compatibility
        self.order_parameter = [op.as_dict() for op in order_parameter]
        self.energy_model = energy_model
        self.enum_kwargs = kwargs

    @staticmethod
    def determine_min_cell(disordered_structure):
        """
        Determine the smallest supercell that is able to enumerate
        the provided structure with the given order parameter
        """

        def lcm(n1, n2):
            """
            Find least common multiple of two numbers
            """
            return n1 * n2 / gcd(n1, n2)

        # assumes all order parameters for a given species are the same
        mag_species_order_parameter = {}
        mag_species_occurrences = {}
        for idx, site in enumerate(disordered_structure):
            if not site.is_ordered:
                op = max(site.species.values())
                # this very hacky bit of code only works because we know
                # that on disordered sites in this class, all species are the same
                # but have different spins, and this is comma-delimited
                sp = str(list(site.species.keys())[0]).split(",", maxsplit=1)[0]
                if sp in mag_species_order_parameter:
                    mag_species_occurrences[sp] += 1
                else:
                    mag_species_order_parameter[sp] = op
                    mag_species_occurrences[sp] = 1

        smallest_n = []

        for sp, order_parameter in mag_species_order_parameter.items():
            denom = Fraction(order_parameter).limit_denominator(100).denominator
            num_atom_per_specie = mag_species_occurrences[sp]
            n_gcd = gcd(denom, num_atom_per_specie)
            smallest_n.append(lcm(int(n_gcd), denom) / n_gcd)

        return max(smallest_n)

    @staticmethod
    def _add_dummy_species(structure, order_parameters):
        """
        :param structure: ordered Structure
        :param order_parameters: list of MagOrderParameterConstraints
        :return: A structure decorated with disordered
            DummySpecies on which to perform the enumeration.
            Note that the DummySpecies are super-imposed on
            to the original sites, to make it easier to
            retrieve the original site after enumeration is
            performed (this approach is preferred over a simple
            mapping since multiple species may have the same
            DummySpecies, depending on the constraints specified).
            This approach can also preserve site properties even after
            enumeration.
        """

        dummy_struct = structure.copy()

        def generate_dummy_specie():
            """
            Generator which returns DummySpecies symbols Mma, Mmb, etc.
            """
            subscript_length = 1
            while True:
                for subscript in product(ascii_lowercase, repeat=subscript_length):
                    yield "Mm" + "".join(subscript)
                subscript_length += 1

        dummy_species_gen = generate_dummy_specie()

        # one dummy species for each order parameter constraint
        dummy_species_symbols = [next(dummy_species_gen) for i in range(len(order_parameters))]
        dummy_species = [
            {
                DummySpecies(symbol, properties={"spin": Spin.up}): constraint.order_parameter,
                DummySpecies(symbol, properties={"spin": Spin.down}): 1 - constraint.order_parameter,
            }
            for symbol, constraint in zip(dummy_species_symbols, order_parameters)
        ]

        for idx, site in enumerate(dummy_struct):
            satisfies_constraints = [c.satisfies_constraint(site) for c in order_parameters]
            if satisfies_constraints.count(True) > 1:
                # site should either not satisfy any constraints, or satisfy
                # one constraint
                raise ValueError(f"Order parameter constraints conflict for site: {site.specie}, {site.properties}")
            if any(satisfies_constraints):
                dummy_specie_idx = satisfies_constraints.index(True)
                dummy_struct.append(dummy_species[dummy_specie_idx], site.coords, site.lattice)

        return dummy_struct

    @staticmethod
    def _remove_dummy_species(structure):
        """
        :return: Structure with dummy species removed, but
        their corresponding spin properties merged with the
        original sites. Used after performing enumeration.
        """
        if not structure.is_ordered:
            raise Exception("Something went wrong with enumeration.")

        sites_to_remove = []
        logger.debug(f"Dummy species structure:\n{structure}")
        for idx, site in enumerate(structure):
            if isinstance(site.specie, DummySpecies):
                sites_to_remove.append(idx)
                spin = site.specie._properties.get("spin", None)
                neighbors = structure.get_neighbors(
                    site,
                    0.05,  # arbitrary threshold, needs to be << any bond length
                    # but >> floating point precision issues
                    include_index=True,
                )
                if len(neighbors) != 1:
                    raise Exception(f"This shouldn't happen, found neighbors: {neighbors}")
                orig_site_idx = neighbors[0][2]
                orig_specie = structure[orig_site_idx].specie
                new_specie = Species(
                    orig_specie.symbol,
                    getattr(orig_specie, "oxi_state", None),
                    properties={"spin": spin},
                )
                structure.replace(
                    orig_site_idx,
                    new_specie,
                    properties=structure[orig_site_idx].properties,
                )
        structure.remove_sites(sites_to_remove)
        logger.debug(f"Structure with dummy species removed:\n{structure}")
        return structure

    def _add_spin_magnitudes(self, structure):
        """
        Replaces Spin.up/Spin.down with spin magnitudes specified
        by mag_species_spin.
        :param structure:
        :return:
        """
        for idx, site in enumerate(structure):
            if getattr(site.specie, "_properties", None):
                spin = site.specie._properties.get("spin", None)
                sign = int(spin) if spin else 0
                if spin:
                    new_properties = site.specie._properties.copy()
                    # this very hacky bit of code only works because we know
                    # that on disordered sites in this class, all species are the same
                    # but have different spins, and this is comma-delimited
                    sp = str(site.specie).split(",", maxsplit=1)[0]
                    new_properties.update({"spin": sign * self.mag_species_spin.get(sp, 0)})
                    new_specie = Species(
                        site.specie.symbol,
                        getattr(site.specie, "oxi_state", None),
                        new_properties,
                    )
                    structure.replace(idx, new_specie, properties=site.properties)
        logger.debug(f"Structure with spin magnitudes:\n{structure}")
        return structure

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Apply MagOrderTransformation to an input structure.
        :param structure: Any ordered structure.
        :param return_ranked_list: As in other Transformations.
        :return:
        """

        if not structure.is_ordered:
            raise ValueError("Create an ordered approximation of your  input structure first.")

        # retrieve order parameters
        order_parameters = [MagOrderParameterConstraint.from_dict(op_dict) for op_dict in self.order_parameter]
        # add dummy species on which to perform enumeration
        structure = self._add_dummy_species(structure, order_parameters)

        # trivial case
        if structure.is_ordered:
            structure = self._remove_dummy_species(structure)
            return [structure] if return_ranked_list > 1 else structure

        enum_kwargs = self.enum_kwargs.copy()

        enum_kwargs["min_cell_size"] = max(int(self.determine_min_cell(structure)), enum_kwargs.get("min_cell_size", 1))

        if enum_kwargs.get("max_cell_size", None):
            if enum_kwargs["min_cell_size"] > enum_kwargs["max_cell_size"]:
                warnings.warn(
                    "Specified max cell size ({}) is smaller "
                    "than the minimum enumerable cell size ({}), "
                    "changing max cell size to {}".format(
                        enum_kwargs["max_cell_size"],
                        enum_kwargs["min_cell_size"],
                        enum_kwargs["min_cell_size"],
                    )
                )
                enum_kwargs["max_cell_size"] = enum_kwargs["min_cell_size"]
        else:
            enum_kwargs["max_cell_size"] = enum_kwargs["min_cell_size"]

        t = EnumerateStructureTransformation(**enum_kwargs)

        alls = t.apply_transformation(structure, return_ranked_list=return_ranked_list)

        # handle the fact that EnumerateStructureTransformation can either
        # return a single Structure or a list
        if isinstance(alls, Structure):
            # remove dummy species and replace Spin.up or Spin.down
            # with spin magnitudes given in mag_species_spin arg
            alls = self._remove_dummy_species(alls)
            alls = self._add_spin_magnitudes(alls)
        else:
            for idx, _ in enumerate(alls):
                alls[idx]["structure"] = self._remove_dummy_species(alls[idx]["structure"])
                alls[idx]["structure"] = self._add_spin_magnitudes(alls[idx]["structure"])

        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        if num_to_return == 1 or not return_ranked_list:
            return alls[0]["structure"] if num_to_return else alls

        # remove duplicate structures and group according to energy model
        m = StructureMatcher(comparator=SpinComparator())

        def key(x):
            return SpacegroupAnalyzer(x, 0.1).get_space_group_number()

        out = []
        for _, g in groupby(sorted((d["structure"] for d in alls), key=key), key):
            g = list(g)
            grouped = m.group_structures(g)
            out.extend([{"structure": g[0], "energy": self.energy_model.get_energy(g[0])} for g in grouped])

        self._all_structures = sorted(out, key=lambda d: d["energy"])

        return self._all_structures[0:num_to_return]

    def __str__(self):
        return "MagOrderingTransformation"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True


def _find_codopant(target, oxidation_state, allowed_elements=None):
    """
    Finds the element from "allowed elements" that (i) possesses the desired
    "oxidation state" and (ii) is closest in ionic radius to the target specie

    Args:
        target: (Species) provides target ionic radius.
        oxidation_state: (float) codopant oxidation state.
        allowed_elements: ([str]) List of allowed elements. If None,
            all elements are tried.

    Returns:
        (Species) with oxidation_state that has ionic radius closest to
        target.
    """
    ref_radius = target.ionic_radius
    candidates = []
    symbols = allowed_elements or [el.symbol for el in Element]
    for sym in symbols:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sp = Species(sym, oxidation_state)
                r = sp.ionic_radius
                if r is not None:
                    candidates.append((r, sp))
        except Exception:
            pass
    return min(candidates, key=lambda l: abs(l[0] / ref_radius - 1))[1]


class DopingTransformation(AbstractTransformation):
    """
    A transformation that performs doping of a structure.
    """

    def __init__(
        self,
        dopant,
        ionic_radius_tol=float("inf"),
        min_length=10,
        alio_tol=0,
        codopant=False,
        max_structures_per_enum=100,
        allowed_doping_species=None,
        **kwargs,
    ):
        r"""
        Args:
            dopant (Species-like): E.g., Al3+. Must have oxidation state.
            ionic_radius_tol (float): E.g., Fractional allowable ionic radii
                mismatch for dopant to fit into a site. Default of inf means
                that any dopant with the right oxidation state is allowed.
            min_Length (float): Min. lattice parameter between periodic
                images of dopant. Defaults to 10A for now.
            alio_tol (int): If this is not 0, attempt will be made to dope
                sites with oxidation_states +- alio_tol of the dopant. E.g.,
                1 means that the ions like Ca2+ and Ti4+ are considered as
                potential doping sites for Al3+.
            codopant (bool): If True, doping will be carried out with a
                codopant to maintain charge neutrality. Otherwise, vacancies
                will be used.
            max_structures_per_enum (float): Maximum number of structures to
                return per enumeration. Note that there can be more than one
                candidate doping site, and each site enumeration will return at
                max max_structures_per_enum structures. Defaults to 100.
            allowed_doping_species (list): Species that are allowed to be
                doping sites. This is an inclusionary list. If specified,
                any sites which are not
            **kwargs:
                Same keyword args as :class:`EnumerateStructureTransformation`,
                i.e., min_cell_size, etc.
        """
        self.dopant = get_el_sp(dopant)
        self.ionic_radius_tol = ionic_radius_tol
        self.min_length = min_length
        self.alio_tol = alio_tol
        self.codopant = codopant
        self.max_structures_per_enum = max_structures_per_enum
        self.allowed_doping_species = allowed_doping_species
        self.kwargs = kwargs

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Args:
            structure (Structure): Input structure to dope

        Returns:
            [{"structure": Structure, "energy": float}]
        """
        comp = structure.composition
        logger.info(f"Composition: {comp}")

        for sp in comp:
            try:
                sp.oxi_state
            except AttributeError:
                analyzer = BVAnalyzer()
                structure = analyzer.get_oxi_state_decorated_structure(structure)
                comp = structure.composition
                break

        ox = self.dopant.oxi_state
        radius = self.dopant.ionic_radius

        compatible_species = [
            sp for sp in comp if sp.oxi_state == ox and abs(sp.ionic_radius / radius - 1) < self.ionic_radius_tol
        ]

        if (not compatible_species) and self.alio_tol:
            # We only consider aliovalent doping if there are no compatible
            # isovalent species.
            compatible_species = [
                sp
                for sp in comp
                if abs(sp.oxi_state - ox) <= self.alio_tol
                and abs(sp.ionic_radius / radius - 1) < self.ionic_radius_tol
                and sp.oxi_state * ox >= 0
            ]

        if self.allowed_doping_species is not None:
            # Only keep allowed doping species.
            compatible_species = [
                sp for sp in compatible_species if sp in [get_el_sp(s) for s in self.allowed_doping_species]
            ]

        logger.info(f"Compatible species: {compatible_species}")

        lengths = structure.lattice.abc
        scaling = [max(1, int(round(math.ceil(self.min_length / x)))) for x in lengths]
        logger.info(f"Lengths are {str(lengths)}")
        logger.info(f"Scaling = {str(scaling)}")

        all_structures = []
        t = EnumerateStructureTransformation(**self.kwargs)

        for sp in compatible_species:
            supercell = structure * scaling
            nsp = supercell.composition[sp]
            if sp.oxi_state == ox:
                supercell.replace_species({sp: {sp: (nsp - 1) / nsp, self.dopant: 1 / nsp}})
                logger.info(f"Doping {sp} for {self.dopant} at level {1 / nsp:.3f}")
            elif self.codopant:
                codopant = _find_codopant(sp, 2 * sp.oxi_state - ox)
                supercell.replace_species({sp: {sp: (nsp - 2) / nsp, self.dopant: 1 / nsp, codopant: 1 / nsp}})
                logger.info(f"Doping {sp} for {self.dopant} + {codopant} at level {1 / nsp:.3f}")
            elif abs(sp.oxi_state) < abs(ox):
                # Strategy: replace the target species with a
                # combination of dopant and vacancy.
                # We will choose the lowest oxidation state species as a
                # vacancy compensation species as it is likely to be lower in
                # energy
                sp_to_remove = min(
                    (s for s in comp if s.oxi_state * ox > 0),
                    key=lambda ss: abs(ss.oxi_state),
                )

                if sp_to_remove == sp:
                    common_charge = lcm(int(abs(sp.oxi_state)), int(abs(ox)))
                    ndopant = common_charge / abs(ox)
                    nsp_to_remove = common_charge / abs(sp.oxi_state)
                    logger.info("Doping %d %s with %d %s." % (nsp_to_remove, sp, ndopant, self.dopant))
                    supercell.replace_species(
                        {
                            sp: {
                                sp: (nsp - nsp_to_remove) / nsp,
                                self.dopant: ndopant / nsp,
                            }
                        }
                    )
                else:
                    ox_diff = int(abs(round(sp.oxi_state - ox)))
                    vac_ox = int(abs(sp_to_remove.oxi_state))
                    common_charge = lcm(vac_ox, ox_diff)
                    ndopant = common_charge / ox_diff
                    nx_to_remove = common_charge / vac_ox
                    nx = supercell.composition[sp_to_remove]
                    logger.info(
                        "Doping %d %s with %s and removing %d %s."
                        % (ndopant, sp, self.dopant, nx_to_remove, sp_to_remove)
                    )
                    supercell.replace_species(
                        {
                            sp: {sp: (nsp - ndopant) / nsp, self.dopant: ndopant / nsp},
                            sp_to_remove: {sp_to_remove: (nx - nx_to_remove) / nx},
                        }
                    )
            elif abs(sp.oxi_state) > abs(ox):
                # Strategy: replace the target species with dopant and also
                # remove some opposite charged species for charge neutrality
                if ox > 0:
                    sp_to_remove = max(supercell.composition.keys(), key=lambda el: el.X)
                else:
                    sp_to_remove = min(supercell.composition.keys(), key=lambda el: el.X)
                # Confirm species are of opposite oxidation states.
                assert sp_to_remove.oxi_state * sp.oxi_state < 0

                ox_diff = int(abs(round(sp.oxi_state - ox)))
                anion_ox = int(abs(sp_to_remove.oxi_state))
                nx = supercell.composition[sp_to_remove]
                common_charge = lcm(anion_ox, ox_diff)
                ndopant = common_charge / ox_diff
                nx_to_remove = common_charge / anion_ox
                logger.info(
                    "Doping %d %s with %s and removing %d %s." % (ndopant, sp, self.dopant, nx_to_remove, sp_to_remove)
                )
                supercell.replace_species(
                    {
                        sp: {sp: (nsp - ndopant) / nsp, self.dopant: ndopant / nsp},
                        sp_to_remove: {sp_to_remove: (nx - nx_to_remove) / nx},
                    }
                )

            ss = t.apply_transformation(supercell, return_ranked_list=self.max_structures_per_enum)
            logger.info(f"{len(ss)} distinct structures")
            all_structures.extend(ss)

        logger.info(f"Total {len(all_structures)} doped structures")
        if return_ranked_list:
            return all_structures[:return_ranked_list]

        return all_structures[0]["structure"]

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True


class SlabTransformation(AbstractTransformation):
    """
    A transformation that creates a slab from a structure.
    """

    def __init__(
        self,
        miller_index,
        min_slab_size,
        min_vacuum_size,
        lll_reduce=False,
        center_slab=False,
        in_unit_planes=False,
        primitive=True,
        max_normal_search=None,
        shift=0,
        tol=0.1,
    ):
        """
        Args:
            miller_index (3-tuple or list): miller index of slab
            min_slab_size (float): minimum slab size in angstroms
            min_vacuum_size (float): minimum size of vacuum
            lll_reduce (bool): whether to apply LLL reduction
            center_slab (bool): whether to center the slab
            primitive (bool): whether to reduce slabs to most primitive cell
            max_normal_search (int): maximum index to include in linear
                combinations of indices to find c lattice vector orthogonal
                to slab surface
            shift (float): shift to get termination
            tol (float): tolerance for primitive cell finding
        """
        self.miller_index = miller_index
        self.min_slab_size = min_slab_size
        self.min_vacuum_size = min_vacuum_size
        self.lll_reduce = lll_reduce
        self.center_slab = center_slab
        self.in_unit_planes = in_unit_planes
        self.primitive = primitive
        self.max_normal_search = max_normal_search
        self.shift = shift
        self.tol = tol

    def apply_transformation(self, structure):
        """
        Applies the transformation.

        Args:
            structure: Input Structure

        Returns:
            Slab Structures.
        """
        sg = SlabGenerator(
            structure,
            self.miller_index,
            self.min_slab_size,
            self.min_vacuum_size,
            self.lll_reduce,
            self.center_slab,
            self.in_unit_planes,
            self.primitive,
            self.max_normal_search,
        )
        slab = sg.get_slab(self.shift, self.tol)
        return slab

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: False"""
        return False


class DisorderOrderedTransformation(AbstractTransformation):
    """
    Not to be confused with OrderDisorderedTransformation,
    this transformation attempts to obtain a
    *disordered* structure from an input ordered structure.
    This may or may not be physically plausible, further
    inspection of the returned structures is advised.
    The main purpose for this transformation is for structure
    matching to crystal prototypes for structures that have
    been derived from a parent prototype structure by
    substitutions or alloying additions.
    """

    def __init__(self, max_sites_to_merge=2):
        """
        Args:
            max_sites_to_merge: only merge this number of sites together
        """
        self.max_sites_to_merge = max_sites_to_merge

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Args:
            structure: ordered structure
            return_ranked_list: as in other pymatgen Transformations

        Returns:
            Transformed disordered structure(s)
        """

        if not structure.is_ordered:
            raise ValueError("This transformation is for disordered structures only.")

        partitions = self._partition_species(structure.composition, max_components=self.max_sites_to_merge)
        disorder_mappings = self._get_disorder_mappings(structure.composition, partitions)

        disordered_structures = []
        for mapping in disorder_mappings:
            disordered_structure = structure.copy()
            disordered_structure.replace_species(mapping)
            disordered_structures.append({"structure": disordered_structure, "mapping": mapping})

        if len(disordered_structures) == 0:
            return None
        if not return_ranked_list:
            return disordered_structures[0]["structure"]
        if len(disordered_structures) > return_ranked_list:
            disordered_structures = disordered_structures[0:return_ranked_list]
        return disordered_structures

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True

    @staticmethod
    def _partition_species(composition, max_components=2):
        """
        Private method to split a list of species into
        various partitions.
        """

        def _partition(collection):
            # thanks https://stackoverflow.com/a/30134039

            if len(collection) == 1:
                yield [collection]
                return

            first = collection[0]
            for smaller in _partition(collection[1:]):
                # insert `first` in each of the subpartition's subsets
                for n, subset in enumerate(smaller):
                    yield smaller[:n] + [[first] + subset] + smaller[n + 1 :]
                # put `first` in its own subset
                yield [[first]] + smaller

        def _sort_partitions(partitions_to_sort):
            """
            Sort partitions by those we want to check first
            (typically, merging two sites into one is the
            one to try first).
            """

            partition_indices = [(idx, [len(p) for p in partition]) for idx, partition in enumerate(partitions_to_sort)]

            # sort by maximum length of partition first (try smallest maximums first)
            # and secondarily by number of partitions (most partitions first, i.e.
            # create the 'least disordered' structures first)
            partition_indices = sorted(partition_indices, key=lambda x: (max(x[1]), -len(x[1])))

            # merge at most max_component sites,
            # e.g. merge at most 2 species into 1 disordered site
            partition_indices = [x for x in partition_indices if max(x[1]) <= max_components]

            partition_indices.pop(0)  # this is just the input structure

            sorted_partitions = [partitions_to_sort[x[0]] for x in partition_indices]

            return sorted_partitions

        collection = list(composition.keys())
        partitions = list(_partition(collection))
        partitions = _sort_partitions(partitions)

        return partitions

    @staticmethod
    def _get_disorder_mappings(composition, partitions):
        """
        Private method to obtain the mapping to create
        a disordered structure from a given partition.
        """

        def _get_replacement_dict_from_partition(partition):
            d = {}  # to be passed to Structure.replace_species()
            for sp_list in partition:
                if len(sp_list) > 1:
                    total_occ = sum(composition[sp] for sp in sp_list)
                    merged_comp = {sp: composition[sp] / total_occ for sp in sp_list}
                    for sp in sp_list:
                        d[sp] = merged_comp
            return d

        disorder_mapping = [_get_replacement_dict_from_partition(p) for p in partitions]

        return disorder_mapping


class GrainBoundaryTransformation(AbstractTransformation):
    """
    A transformation that creates a gb from a bulk structure.
    """

    def __init__(
        self,
        rotation_axis,
        rotation_angle,
        expand_times=4,
        vacuum_thickness=0.0,
        ab_shift=None,
        normal=False,
        ratio=True,
        plane=None,
        max_search=20,
        tol_coi=1.0e-8,
        rm_ratio=0.7,
        quick_gen=False,
    ):
        """
        Args:
            rotation_axis (list): Rotation axis of GB in the form of a list of integer
                e.g.: [1, 1, 0]
            rotation_angle (float, in unit of degree): rotation angle used to generate GB.
                Make sure the angle is accurate enough. You can use the enum* functions
                in this class to extract the accurate angle.
                e.g.: The rotation angle of sigma 3 twist GB with the rotation axis
                [1, 1, 1] and GB plane (1, 1, 1) can be 60.000000000 degree.
                If you do not know the rotation angle, but know the sigma value, we have
                provide the function get_rotation_angle_from_sigma which is able to return
                all the rotation angles of sigma value you provided.
            expand_times (int): The multiple times used to expand one unit grain to larger grain.
                This is used to tune the grain length of GB to warrant that the two GBs in one
                cell do not interact with each other. Default set to 4.
            vacuum_thickness (float): The thickness of vacuum that you want to insert between
                two grains of the GB. Default to 0.
            ab_shift (list of float, in unit of a, b vectors of Gb): in plane shift of two grains
            normal (logic):
                determine if need to require the c axis of top grain (first transformation matrix)
                perperdicular to the surface or not.
                default to false.
            ratio (list of integers): lattice axial ratio.
                If True, will try to determine automatically from structure.
                For cubic system, ratio is not needed and can be set to None.
                For tetragonal system, ratio = [mu, mv], list of two integers,
                that is, mu/mv = c2/a2. If it is irrational, set it to None.
                For orthorhombic system, ratio = [mu, lam, mv], list of three integers,
                    that is, mu:lam:mv = c2:b2:a2. If irrational for one axis, set it to None.
                e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
                For rhombohedral system, ratio = [mu, mv], list of two integers,
                that is, mu/mv is the ratio of (1+2*cos(alpha))/cos(alpha).
                If irrational, set it to None.
                For hexagonal system, ratio = [mu, mv], list of two integers,
                that is, mu/mv = c2/a2. If it is irrational, set it to none.
            plane (list): Grain boundary plane in the form of a list of integers
                e.g.: [1, 2, 3]. If none, we set it as twist GB. The plane will be perpendicular
                to the rotation axis.
            max_search (int): max search for the GB lattice vectors that give the smallest GB
                lattice. If normal is true, also max search the GB c vector that perpendicular
                to the plane. For complex GB, if you want to speed up, you can reduce this value.
                But too small of this value may lead to error.
            tol_coi (float): tolerance to find the coincidence sites. When making approximations to
                the ratio needed to generate the GB, you probably need to increase this tolerance to
                obtain the correct number of coincidence sites. To check the number of coincidence
                sites are correct or not, you can compare the generated Gb object's sigma with enum*
                sigma values (what user expected by input).
            rm_ratio (float): the criteria to remove the atoms which are too close with each other.
                rm_ratio * bond_length of bulk system is the criteria of bond length, below which the atom
                will be removed. Default to 0.7.
            quick_gen (bool): whether to quickly generate a supercell, if set to true, no need to
                find the smallest cell.
        Returns:
           Grain boundary structure (gb (Structure) object).
        """
        self.rotation_axis = rotation_axis
        self.rotation_angle = rotation_angle
        self.expand_times = expand_times
        self.vacuum_thickness = vacuum_thickness
        self.ab_shift = ab_shift or [0, 0]
        self.normal = normal
        self.ratio = ratio
        self.plane = plane
        self.max_search = max_search
        self.tol_coi = tol_coi
        self.rm_ratio = rm_ratio
        self.quick_gen = quick_gen

    def apply_transformation(self, structure):
        """
        Applies the transformation.

        Args:
            structure: Input Structure
            return_ranked_list: Number of structures to return.

        Returns:
            Grain boundary Structures.
        """
        gbg = GrainBoundaryGenerator(structure)
        gb_struct = gbg.gb_from_parameters(
            self.rotation_axis,
            self.rotation_angle,
            expand_times=self.expand_times,
            vacuum_thickness=self.vacuum_thickness,
            ab_shift=self.ab_shift,
            normal=self.normal,
            ratio=gbg.get_ratio() if self.ratio is True else self.ratio,
            plane=self.plane,
            max_search=self.max_search,
            tol_coi=self.tol_coi,
            rm_ratio=self.rm_ratio,
            quick_gen=self.quick_gen,
        )
        return gb_struct

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: False"""
        return False


class CubicSupercellTransformation(AbstractTransformation):
    """
    A transformation that aims to generate a nearly cubic supercell structure
    from a structure.

    The algorithm solves for a transformation matrix that makes the supercell
    cubic. The matrix must have integer entries, so entries are rounded (in such
    a way that forces the matrix to be nonsingular). From the supercell
    resulting from this transformation matrix, vector projections are used to
    determine the side length of the largest cube that can fit inside the
    supercell. The algorithm will iteratively increase the size of the supercell
    until the largest inscribed cube's side length is at least 'min_length'
    and the number of atoms in the supercell falls in the range
    ``min_atoms < n < max_atoms``.
    """

    def __init__(
        self,
        min_atoms: Optional[int] = None,
        max_atoms: Optional[int] = None,
        min_length: float = 15.0,
        force_diagonal: bool = False,
    ):
        """
        Args:
            max_atoms: Maximum number of atoms allowed in the supercell.
            min_atoms: Minimum number of atoms allowed in the supercell.
            min_length: Minimum length of the smallest supercell lattice vector.
            force_diagonal: If True, return a transformation with a diagonal
                transformation matrix.
        """
        self.min_atoms = min_atoms if min_atoms else -np.Inf
        self.max_atoms = max_atoms if max_atoms else np.Inf
        self.min_length = min_length
        self.force_diagonal = force_diagonal
        self.transformation_matrix = None

    def apply_transformation(self, structure: Structure) -> Structure:
        """
        The algorithm solves for a transformation matrix that makes the
        supercell cubic. The matrix must have integer entries, so entries are
        rounded (in such a way that forces the matrix to be nonsingular). From
        the supercell resulting from this transformation matrix, vector
        projections are used to determine the side length of the largest cube
        that can fit inside the supercell. The algorithm will iteratively
        increase the size of the supercell until the largest inscribed cube's
        side length is at least 'num_nn_dists' times the nearest neighbor
        distance and the number of atoms in the supercell falls in the range
        defined by min_atoms and max_atoms.

        Returns:
            supercell: Transformed supercell.
        """

        lat_vecs = structure.lattice.matrix

        # boolean for if a sufficiently large supercell has been created
        sc_not_found = True

        if self.force_diagonal:
            scale = self.min_length / np.array(structure.lattice.abc)
            self.transformation_matrix = np.diag(np.ceil(scale).astype(int))
            st = SupercellTransformation(self.transformation_matrix)
            return st.apply_transformation(structure)

        # target_threshold is used as the desired cubic side lengths
        target_sc_size = self.min_length
        while sc_not_found:
            target_sc_lat_vecs = np.eye(3, 3) * target_sc_size
            self.transformation_matrix = target_sc_lat_vecs @ np.linalg.inv(lat_vecs)

            # round the entries of T and force T to be nonsingular
            self.transformation_matrix = _round_and_make_arr_singular(self.transformation_matrix)  # type: ignore

            proposed_sc_lat_vecs = self.transformation_matrix @ lat_vecs  # type: ignore

            # Find the shortest dimension length and direction
            a = proposed_sc_lat_vecs[0]
            b = proposed_sc_lat_vecs[1]
            c = proposed_sc_lat_vecs[2]

            length1_vec = c - _proj(c, a)  # a-c plane
            length2_vec = a - _proj(a, c)
            length3_vec = b - _proj(b, a)  # b-a plane
            length4_vec = a - _proj(a, b)
            length5_vec = b - _proj(b, c)  # b-c plane
            length6_vec = c - _proj(c, b)
            length_vecs = np.array(
                [
                    length1_vec,
                    length2_vec,
                    length3_vec,
                    length4_vec,
                    length5_vec,
                    length6_vec,
                ]
            )

            # Get number of atoms
            st = SupercellTransformation(self.transformation_matrix)
            superstructure = st.apply_transformation(structure)
            num_at = superstructure.num_sites

            # Check if constraints are satisfied
            if (
                np.min(np.linalg.norm(length_vecs, axis=1)) >= self.min_length
                and self.min_atoms <= num_at <= self.max_atoms
            ):
                return superstructure
            # Increase threshold until proposed supercell meets requirements
            target_sc_size += 0.1
            if num_at > self.max_atoms:
                raise AttributeError(
                    "While trying to solve for the supercell, the max "
                    "number of atoms was exceeded. Try lowering the number"
                    "of nearest neighbor distances."
                )
        raise AttributeError("Unable to find cubic supercell")

    @property
    def inverse(self):
        """
        Returns:
            None
        """
        return None

    @property
    def is_one_to_many(self):
        """
        Returns:
            False
        """
        return False


class AddAdsorbateTransformation(AbstractTransformation):
    """
    Create absorbate structures.
    """

    def __init__(
        self,
        adsorbate,
        selective_dynamics=False,
        height=0.9,
        mi_vec=None,
        repeat=None,
        min_lw=5.0,
        translate=True,
        reorient=True,
        find_args=None,
    ):
        """
        Use AdsorbateSiteFinder to add an absorbate to a slab.

        Args:
            adsorbate (Molecule): molecule to add as adsorbate
            selective_dynamics (bool): flag for whether to assign
                non-surface sites as fixed for selective dynamics
            height (float): height criteria for selection of surface sites
            mi_vec : vector corresponding to the vector
                concurrent with the miller index, this enables use with
                slabs that have been reoriented, but the miller vector
                must be supplied manually
            repeat (3-tuple or list): repeat argument for supercell generation
            min_lw (float): minimum length and width of the slab, only used
                if repeat is None
            translate (bool): flag on whether to translate the molecule so
                that its CoM is at the origin prior to adding it to the surface
            reorient (bool): flag on whether or not to reorient adsorbate
                along the miller index
            find_args (dict): dictionary of arguments to be passed to the
                call to self.find_adsorption_sites, e.g. {"distance":2.0}
        """
        self.adsorbate = adsorbate
        self.selective_dynamics = selective_dynamics
        self.height = height
        self.mi_vec = mi_vec
        self.repeat = repeat
        self.min_lw = min_lw
        self.translate = translate
        self.reorient = reorient
        self.find_args = find_args

    def apply_transformation(self, structure, return_ranked_list=False):
        """

        Args:
            structure: Must be a Slab structure
            return_ranked_list:  Whether or not multiple structures are
                returned. If return_ranked_list is a number, up to that number of
                structures is returned.

        Returns: Slab with adsorbate

        """

        sitefinder = AdsorbateSiteFinder(
            structure,
            selective_dynamics=self.selective_dynamics,
            height=self.height,
            mi_vec=self.mi_vec,
        )

        structures = sitefinder.generate_adsorption_structures(
            self.adsorbate,
            repeat=self.repeat,
            min_lw=self.min_lw,
            translate=self.translate,
            reorient=self.reorient,
            find_args=self.find_args,
        )

        if not return_ranked_list:
            return structures[0]
        return [{"structure": structure} for structure in structures[:return_ranked_list]]

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True


def _round_and_make_arr_singular(arr: np.ndarray) -> np.ndarray:
    """
    This function rounds all elements of a matrix to the nearest integer,
    unless the rounding scheme causes the matrix to be singular, in which
    case elements of zero rows or columns in the rounded matrix with the
    largest absolute valued magnitude in the unrounded matrix will be
    rounded to the next integer away from zero rather than to the
    nearest integer.

    The transformation is as follows. First, all entries in 'arr' will be
    rounded to the nearest integer to yield 'arr_rounded'. If 'arr_rounded'
    has any zero rows, then one element in each zero row of 'arr_rounded'
    corresponding to the element in 'arr' of that row with the largest
    absolute valued magnitude will be rounded to the next integer away from
    zero (see the '_round_away_from_zero(x)' function) rather than the
    nearest integer. This process is then repeated for zero columns. Also
    note that if 'arr' already has zero rows or columns, then this function
    will not change those rows/columns.

    Args:
        arr: Input matrix

    Returns:
        Transformed matrix.
    """

    def round_away_from_zero(x):
        """
        Returns 'x' rounded to the next integer away from 0.
        If 'x' is zero, then returns zero.
        E.g. -1.2 rounds to -2.0. 1.2 rounds to 2.0.
        """
        abs_x = abs(x)
        return math.ceil(abs_x) * (abs_x / x) if x != 0 else 0

    arr_rounded = np.around(arr)

    # Zero rows in 'arr_rounded' make the array singular, so force zero rows to
    # be nonzero
    if (~arr_rounded.any(axis=1)).any():
        # Check for zero rows in T_rounded

        # indices of zero rows
        zero_row_idxs = np.where(~arr_rounded.any(axis=1))[0]

        for zero_row_idx in zero_row_idxs:  # loop over zero rows
            zero_row = arr[zero_row_idx, :]

            # Find the element of the zero row with the largest absolute
            # magnitude in the original (non-rounded) array (i.e. 'arr')
            matches = np.absolute(zero_row) == np.amax(np.absolute(zero_row))
            col_idx_to_fix = np.where(matches)[0]

            # Break ties for the largest absolute magnitude
            r_idx = np.random.randint(len(col_idx_to_fix))
            col_idx_to_fix = col_idx_to_fix[r_idx]

            # Round the chosen element away from zero
            arr_rounded[zero_row_idx, col_idx_to_fix] = round_away_from_zero(arr[zero_row_idx, col_idx_to_fix])

    # Repeat process for zero columns
    if (~arr_rounded.any(axis=0)).any():

        # Check for zero columns in T_rounded
        zero_col_idxs = np.where(~arr_rounded.any(axis=0))[0]
        for zero_col_idx in zero_col_idxs:
            zero_col = arr[:, zero_col_idx]
            matches = np.absolute(zero_col) == np.amax(np.absolute(zero_col))
            row_idx_to_fix = np.where(matches)[0]

            for i in row_idx_to_fix:
                arr_rounded[i, zero_col_idx] = round_away_from_zero(arr[i, zero_col_idx])
    return arr_rounded.astype(int)


class SubstituteSurfaceSiteTransformation(AbstractTransformation):
    """
    Use AdsorptionSiteFinder to perform substitution-type doping on the surface
    and returns all possible configurations where one dopant is substituted
    per surface. Can substitute one surface or both.
    """

    def __init__(
        self,
        atom,
        selective_dynamics=False,
        height=0.9,
        mi_vec=None,
        target_species=None,
        sub_both_sides=False,
        range_tol=1e-2,
        dist_from_surf=0,
    ):
        """
        Args:
            atom (str): atom corresponding to substitutional dopant
            selective_dynamics (bool): flag for whether to assign
                non-surface sites as fixed for selective dynamics
            height (float): height criteria for selection of surface sites
            mi_vec : vector corresponding to the vector
                concurrent with the miller index, this enables use with
                slabs that have been reoriented, but the miller vector
                must be supplied manually
            target_species:  List of specific species to substitute
            sub_both_sides (bool): If true, substitute an equivalent
                site on the other surface
            range_tol (float): Find viable substitution sites at a specific
                distance from the surface +- this tolerance
            dist_from_surf (float): Distance from the surface to find viable
                substitution sites, defaults to 0 to substitute at the surface
        """
        self.atom = atom
        self.selective_dynamics = selective_dynamics
        self.height = height
        self.mi_vec = mi_vec
        self.target_species = target_species
        self.sub_both_sides = sub_both_sides
        self.range_tol = range_tol
        self.dist_from_surf = dist_from_surf

    def apply_transformation(self, structure, return_ranked_list=False):
        """

        Args:
            structure: Must be a Slab structure
            return_ranked_list:  Whether or not multiple structures are
                returned. If return_ranked_list is a number, up to that number of
                structures is returned.

        Returns: Slab with sites substituted

        """

        sitefinder = AdsorbateSiteFinder(
            structure,
            selective_dynamics=self.selective_dynamics,
            height=self.height,
            mi_vec=self.mi_vec,
        )

        structures = sitefinder.generate_substitution_structures(
            self.atom,
            target_species=self.target_species,
            sub_both_sides=self.sub_both_sides,
            range_tol=self.range_tol,
            dist_from_surf=self.dist_from_surf,
        )

        if not return_ranked_list:
            return structures[0]
        return [{"structure": structure} for structure in structures[:return_ranked_list]]

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True


def _proj(b, a):
    """
    Returns vector projection (np.ndarray) of vector b (np.ndarray)
    onto vector a (np.ndarray)
    """
    return (b.T @ (a / np.linalg.norm(a))) * (a / np.linalg.norm(a))


class SQSTransformation(AbstractTransformation):
    """
    A transformation that creates a special quasirandom structure (SQS) from a structure with partial occupancies.
    """

    def __init__(
        self,
        scaling,
        cluster_size_and_shell=None,
        search_time=60,
        directory=None,
        instances=None,
        temperature=1,
        wr=1,
        wn=1,
        wd=0.5,
        tol=1e-3,
        best_only=True,
        remove_duplicate_structures=True,
        reduction_algo="LLL",
    ):
        """
        Args:
            structure (Structure): Disordered pymatgen Structure object
            scaling (int or list): Scaling factor to determine supercell. Two options are possible:
                    a. (preferred) Scales number of atoms, e.g., for a structure with 8 atoms,
                       scaling=4 would lead to a 32 atom supercell
                    b. A sequence of three scaling factors, e.g., [2, 1, 1], which
                       specifies that the supercell should have dimensions 2a x b x c
            cluster_size_and_shell (Optional[Dict[int, int]]): Dictionary of cluster interactions with entries in
                the form number of atoms: nearest neighbor shell
        Keyword Args:
            search_time (float): Time spent looking for the ideal SQS in minutes (default: 60)
            directory (str): Directory to run mcsqs calculation and store files (default: None
                runs calculations in a temp directory)
            instances (int): Specifies the number of parallel instances of mcsqs to run
                (default: number of cpu cores detected by Python)
            temperature (int or float): Monte Carlo temperature (default: 1), "T" in atat code
            wr (int or float): Weight assigned to range of perfect correlation match in objective
                function (default = 1)
            wn (int or float): Multiplicative decrease in weight per additional point in cluster (default: 1)
            wd (int or float): Exponent of decay in weight as function of cluster diameter (default: 0)
            tol (int or float): Tolerance for matching correlations (default: 1e-3)
            best_only (bool): only return structures with lowest objective function
            remove_duplicate_structures (bool): only return unique structures
            reduction_algo (str): The lattice reduction algorithm to use.
                Currently supported options are "niggli" or "LLL".
                "False" does not reduce structure.
        """
        self.scaling = scaling
        self.search_time = search_time
        self.cluster_size_and_shell = cluster_size_and_shell
        self.directory = directory
        self.instances = instances
        self.temperature = temperature
        self.wr = wr
        self.wn = wn
        self.wd = wd
        self.tol = tol
        self.best_only = best_only
        self.remove_duplicate_structures = remove_duplicate_structures
        self.reduction_algo = reduction_algo

    @staticmethod
    def _get_max_neighbor_distance(struc, shell):
        """
        Calculate maximum nearest neighbor distance
        Args:
            struc: pymatgen Structure object
            shell: nearest neighbor shell, such that shell=1 is the first nearest
                neighbor, etc.

        Returns:
            maximum nearest neighbor distance, in angstroms
        """

        mdnn = MinimumDistanceNN()
        distances = []

        for site_num, site in enumerate(struc):
            shell_info = mdnn.get_nn_shell_info(struc, site_num, shell)
            for entry in shell_info:
                image = entry["image"]
                distance = site.distance(struc[entry["site_index"]], jimage=image)
                distances.append(distance)

        return max(distances)

    @staticmethod
    def _get_disordered_substructure(struc_disordered):
        """
        Converts disordered structure into a substructure consisting of only disordered sites
        Args:
            struc_disordered: pymatgen disordered Structure object
        Returns:
            pymatgen Structure object representing a substructure of disordered sites
        """

        disordered_substructure = struc_disordered.copy()

        idx_to_remove = []
        for idx, site in enumerate(disordered_substructure.sites):
            if site.is_ordered:
                idx_to_remove.append(idx)
        disordered_substructure.remove_sites(idx_to_remove)

        return disordered_substructure

    @staticmethod
    def _sqs_cluster_estimate(struc_disordered, cluster_size_and_shell: Optional[Dict[int, int]] = None):
        """
        Set up an ATAT cluster.out file for a given structure and set of constraints
        Args:
            struc_disordered: disordered pymatgen Structure object
            cluster_size_and_shell: dict of integers {cluster: shell}
        Returns:
            dict of {cluster size: distance in angstroms} for mcsqs calculation
        """

        cluster_size_and_shell = cluster_size_and_shell or {2: 3, 3: 2, 4: 1}

        disordered_substructure = SQSTransformation._get_disordered_substructure(struc_disordered)

        clusters = {}
        for cluster_size, shell in cluster_size_and_shell.items():
            max_distance = SQSTransformation._get_max_neighbor_distance(disordered_substructure, shell)
            clusters[cluster_size] = max_distance + 0.01  # add small tolerance

        return clusters

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Applies SQS transformation
        Args:
            structure (pymatgen Structure): pymatgen Structure with partial occupancies
            return_ranked_list (bool): number of structures to return
        Returns:
            pymatgen Structure which is an SQS of the input structure
        """

        if return_ranked_list and self.instances is None:
            raise ValueError("mcsqs has no instances, so cannot return a ranked list")
        if (
            isinstance(return_ranked_list, int)
            and isinstance(self.instances, int)
            and return_ranked_list > self.instances
        ):
            raise ValueError("return_ranked_list cannot be less that number of instances")

        clusters = self._sqs_cluster_estimate(structure, self.cluster_size_and_shell)

        # useful for debugging and understanding
        self._last_used_clusters = clusters

        sqs = run_mcsqs(
            structure=structure,
            clusters=clusters,
            scaling=self.scaling,
            search_time=self.search_time,
            directory=self.directory,
            instances=self.instances,
            temperature=self.temperature,
            wr=self.wr,
            wn=self.wn,
            wd=self.wd,
            tol=self.tol,
        )

        return self._get_unique_bestsqs_strucs(
            sqs,
            best_only=self.best_only,
            return_ranked_list=return_ranked_list,
            remove_duplicate_structures=self.remove_duplicate_structures,
            reduction_algo=self.reduction_algo,
        )

    @staticmethod
    def _get_unique_bestsqs_strucs(sqs, best_only, return_ranked_list, remove_duplicate_structures, reduction_algo):
        """
        Gets unique sqs structures with lowest objective function. Requires an mcsqs output that has been run
            in parallel, otherwise returns Sqs.bestsqs
        Args:
            sqs (Sqs): Sqs class object.
            best_only (bool): only return structures with lowest objective function.
            return_ranked_list (bool): Number of structures to return.
            remove_duplicate_structures (bool): only return unique structures.
            reduction_algo (str): The lattice reduction algorithm to use.
                Currently supported options are "niggli" or "LLL".
                "False" does not reduce structure.
        Returns:
            list of dicts of the form {'structure': Structure, 'objective_function': ...}, unless run in serial
                (returns a single structure Sqs.bestsqs)
        """

        if not return_ranked_list:
            return_struc = sqs.bestsqs

            # reduce structure
            if reduction_algo:
                return_struc = return_struc.get_reduced_structure(reduction_algo=reduction_algo)

            # return just the structure
            return return_struc

        strucs = []
        for d in sqs.allsqs:
            # filter for best structures only if enabled, else use full sqs.all_sqs list
            if (not best_only) or (best_only and d["objective_function"] == sqs.objective_function):
                struc = d["structure"]
                # add temporary objective_function attribute to access objective_function after grouping
                struc.objective_function = d["objective_function"]
                strucs.append(struc)

        if remove_duplicate_structures:
            matcher = StructureMatcher()
            # sort by unique structures ... can take a while for a long list of strucs
            unique_strucs_grouped = matcher.group_structures(strucs)
            # get unique structures only
            strucs = [group[0] for group in unique_strucs_grouped]

        # sort structures by objective function
        strucs.sort(key=lambda x: x.objective_function if isinstance(x.objective_function, float) else -np.inf)

        to_return = [{"structure": struc, "objective_function": struc.objective_function} for struc in strucs]

        for d in to_return:

            # delete temporary objective_function attribute
            del d["structure"].objective_function

            # reduce structure
            if reduction_algo:
                d["structure"] = d["structure"].get_reduced_structure(reduction_algo=reduction_algo)

        if isinstance(return_ranked_list, int):
            return to_return[:return_ranked_list]
        return to_return

    @property
    def inverse(self):
        """Returns: None"""
        return None

    @property
    def is_one_to_many(self):
        """Returns: True"""
        return True


class MonteCarloRattleTransformation(AbstractTransformation):
    r"""
    Uses a Monte Carlo rattle procedure to randomly perturb the sites in a
    structure.

    This class requires the hiPhive package to be installed.

    Rattling atom `i` is carried out as a Monte Carlo move that is accepted with
    a probability determined from the minimum interatomic distance
    :math:`d_{ij}`.  If :math:`\\min(d_{ij})` is smaller than :math:`d_{min}`
    the move is only accepted with a low probability.

    This process is repeated for each atom a number of times meaning
    the magnitude of the final displacements is not *directly*
    connected to `rattle_std`.
    """

    @requires(hiphive, "hiphive is required for MonteCarloRattleTransformation")
    def __init__(self, rattle_std: float, min_distance: float, seed: Optional[int] = None, **kwargs):
        """
        Args:
            rattle_std: Rattle amplitude (standard deviation in normal
                distribution). Note: this value is not *directly* connected to the
                final average displacement for the structures
            min_distance: Interatomic distance used for computing the probability
                for each rattle move.
            seed: Seed for setting up NumPy random state from which random numbers
                are generated. If ``None``, a random seed will be generated
                (default). This option allows the output of this transformation
                to be deterministic.
            **kwargs: Additional keyword arguments to be passed to the hiPhive
                mc_rattle function.
        """
        self.rattle_std = rattle_std
        self.min_distance = min_distance
        self.seed = seed

        if not seed:
            # if seed is None, use a random RandomState seed but make sure
            # we store that the original seed was None
            seed = np.random.randint(1, 1000000000)

        self.random_state = np.random.RandomState(seed)  # pylint: disable=E1101
        self.kwargs = kwargs

    def apply_transformation(self, structure: Structure) -> Structure:
        """
        Apply the transformation.

        Args:
            structure: Input Structure

        Returns:
            Structure with sites perturbed.
        """
        from hiphive.structure_generation.rattle import mc_rattle  # type: ignore

        atoms = AseAtomsAdaptor.get_atoms(structure)
        seed = self.random_state.randint(1, 1000000000)
        displacements = mc_rattle(atoms, self.rattle_std, self.min_distance, seed=seed, **self.kwargs)

        transformed_structure = Structure(
            structure.lattice,
            structure.species,
            structure.cart_coords + displacements,
            coords_are_cartesian=True,
        )

        return transformed_structure

    def __str__(self):
        return f"{__name__} : rattle_std = {self.rattle_std}"

    def __repr__(self):
        return self.__str__()

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
