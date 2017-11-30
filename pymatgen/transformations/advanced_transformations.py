# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import numpy as np
from fractions import Fraction
try:
    from math import gcd
except ImportError:
    from fractions import gcd
from itertools import groupby, product
from string import ascii_lowercase
from warnings import warn
import logging
import math

import warnings
from monty.fractions import lcm
from monty.json import MSONable

from pymatgen.core.periodic_table import Element, Specie, get_el_sp, DummySpecie
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.transformations.standard_transformations import \
    SubstitutionTransformation, OrderDisorderedStructureTransformation
from pymatgen.command_line.enumlib_caller import EnumlibAdaptor, EnumError
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_prediction.substitution_probability import \
    SubstitutionPredictor
from pymatgen.analysis.structure_matcher import StructureMatcher, \
    SpinComparator
from pymatgen.analysis.energy_models import SymmetryModel
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.core.surface import SlabGenerator
from pymatgen.electronic_structure.core import Spin

"""
This module implements more advanced transformations.
"""

__author__ = "Shyue Ping Ong, Stephen Dacek, Anubhav Jain, Matthew Horton"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 24, 2012"


logger = logging.getLogger(__name__)


class ChargeBalanceTransformation(AbstractTransformation):
    """
    This is a transformation that disorders a structure to make it charge
    balanced, given an oxidation state-decorated structure.

    Args:
        charge_balance_sp: specie to add or remove. Currently only removal
            is supported
    """
    def __init__(self, charge_balance_sp):
        self.charge_balance_sp = str(charge_balance_sp)

    def apply_transformation(self, structure):
        charge = structure.charge
        specie = get_el_sp(self.charge_balance_sp)
        num_to_remove = charge / specie.oxi_state
        num_in_structure = structure.composition[specie]
        removal_fraction = num_to_remove / num_in_structure
        if removal_fraction < 0:
            raise ValueError("addition of specie not yet supported by "
                             "ChargeBalanceTransformation")
        trans = SubstitutionTransformation(
            {self.charge_balance_sp: {
                self.charge_balance_sp: 1 - removal_fraction}})
        return trans.apply_transformation(structure)

    def __str__(self):
        return "Charge Balance Transformation : " + \
               "Species to remove = {}".format(str(self.charge_balance_sp))

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False


class SuperTransformation(AbstractTransformation):
    """
    This is a transformation that is inherently one-to-many. It is constructed
    from a list of transformations and returns one structure for each
    transformation. The primary use for this class is extending a transmuter
    object.

    Args:
        transformations ([transformations]): List of transformations to apply
            to a structure. One transformation is applied to each output
            structure.
        nstructures_per_trans (int): If the transformations are one-to-many and,
            nstructures_per_trans structures from each transformation are
            added to the full list. Defaults to 1, i.e., only best structure.
    """

    def __init__(self, transformations, nstructures_per_trans=1):
        self._transformations = transformations
        self.nstructures_per_trans = nstructures_per_trans

    def apply_transformation(self, structure, return_ranked_list=False):
        if not return_ranked_list:
            raise ValueError("SuperTransformation has no single best structure"
                             " output. Must use return_ranked_list")
        structures = []
        for t in self._transformations:
            if t.is_one_to_many:
                for d in t.apply_transformation(
                        structure,
                        return_ranked_list=self.nstructures_per_trans):
                    d["transformation"] = t
                    structures.append(d)
            else:
                structures.append(
                    {"transformation": t,
                     "structure": t.apply_transformation(structure)})
        return structures

    def __str__(self):
        return "Super Transformation : Transformations = " + \
               "{}".format(" ".join([str(t) for t in self._transformations]))

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return True


class MultipleSubstitutionTransformation(object):
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

    def __init__(self, sp_to_replace, r_fraction, substitution_dict,
                 charge_balance_species=None, order=True):
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
        if not return_ranked_list:
            raise ValueError("MultipleSubstitutionTransformation has no single"
                             " best structure output. Must use"
                             " return_ranked_list.")
        outputs = []
        for charge, el_list in self.substitution_dict.items():
            mapping = {}
            if charge > 0:
                sign = "+"
            else:
                sign = "-"
            dummy_sp = "X{}{}".format(str(charge), sign)
            mapping[self.sp_to_replace] = {
                self.sp_to_replace: 1 - self.r_fraction,
                dummy_sp: self.r_fraction}
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
                st = SubstitutionTransformation(
                    {"X{}+".format(str(charge)): "{}{}{}".format(el, charge,
                                                                 sign)})
                new_structure = st.apply_transformation(dummy_structure)
                outputs.append({"structure": new_structure})
        return outputs

    def __str__(self):
        return "Multiple Substitution Transformation : Substitution on " + \
               "{}".format(self.sp_to_replace)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return True


class EnumerateStructureTransformation(AbstractTransformation):
    """
    Order a disordered structure using enumlib. For complete orderings, this
    generally produces fewer structures that the OrderDisorderedStructure
    transformation, and at a much faster speed.

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
            If you are already starting from an experimental cif, refinment
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
    """

    def __init__(self, min_cell_size=1, max_cell_size=1, symm_prec=0.1,
                 refine_structure=False, enum_precision_parameter=0.001,
                 check_ordered_symmetry=True, max_disordered_sites=None,
                 sort_criteria="ewald"):
        self.symm_prec = symm_prec
        self.min_cell_size = min_cell_size
        self.max_cell_size = max_cell_size
        self.refine_structure = refine_structure
        self.enum_precision_parameter = enum_precision_parameter
        self.check_ordered_symmetry = check_ordered_symmetry
        self.max_disordered_sites = max_disordered_sites
        self.sort_criteria = sort_criteria

        if max_cell_size and max_disordered_sites:
            raise ValueError("Cannot set both max_cell_size and "
                             "max_disordered_sites!")

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Return either a single ordered structure or a sequence of all ordered
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
            [hasattr(sp, "oxi_state") and sp.oxi_state != 0 for sp in
             structure.composition.elements]
        )

        structures = None

        if structure.is_ordered:
            warn("Enumeration skipped for structure with composition {} "
                 "because it is ordered".format(structure.composition))
            structures = [structure.copy()]

        if self.max_disordered_sites:
            ndisordered = sum([1 for site in structure if not site.is_ordered])
            if ndisordered > self.max_disordered_sites:
                raise ValueError(
                    "Too many disordered sites! ({} > {})".format(
                        ndisordered, self.max_disordered_sites))
            max_cell_sizes = range(self.min_cell_size, int(
                    math.floor(self.max_disordered_sites / ndisordered)) + 1)

        else:
            max_cell_sizes = [self.max_cell_size]

        for max_cell_size in max_cell_sizes:
            adaptor = EnumlibAdaptor(
                structure, min_cell_size=self.min_cell_size,
                max_cell_size=max_cell_size,
                symm_prec=self.symm_prec, refine_structure=False,
                enum_precision_parameter=self.enum_precision_parameter,
                check_ordered_symmetry=self.check_ordered_symmetry)
            try:
                adaptor.run()
            except EnumError:
                warn("Unable to enumerate for max_cell_size = %d".format(
                    max_cell_size))
            structures = adaptor.structures
            if structures:
                break

        if structures is None:
            raise ValueError("Unable to enumerate")

        original_latt = structure.lattice
        inv_latt = np.linalg.inv(original_latt.matrix)
        ewald_matrices = {}
        all_structures = []
        for s in structures:
            new_latt = s.lattice
            transformation = np.dot(new_latt.matrix, inv_latt)
            transformation = tuple([tuple([int(round(cell)) for cell in row])
                                    for row in transformation])
            if contains_oxidation_state and self.sort_criteria == "ewald":
                if transformation not in ewald_matrices:
                    s_supercell = structure * transformation
                    ewald = EwaldSummation(s_supercell)
                    ewald_matrices[transformation] = ewald
                else:
                    ewald = ewald_matrices[transformation]
                energy = ewald.compute_sub_structure(s)
                all_structures.append({"num_sites": len(s), "energy": energy,
                                       "structure": s})
            else:
                all_structures.append({"num_sites": len(s), "structure": s})

        def sort_func(s):
            return s["energy"] / s["num_sites"] \
                if contains_oxidation_state and self.sort_criteria == "ewald" \
                else s["num_sites"]

        self._all_structures = sorted(all_structures, key=sort_func)

        if return_ranked_list:
            return self._all_structures[0:num_to_return]
        else:
            return self._all_structures[0]["structure"]

    def __str__(self):
        return "EnumerateStructureTransformation"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return True


class SubstitutionPredictorTransformation(AbstractTransformation):
    """
    This transformation takes a structure and uses the structure
    prediction module to find likely site substitutions.

    Args:
        threshold: Threshold for substitution.
        **kwargs: Args for SubstitutionProbability class lambda_table, alpha
    """

    def __init__(self, threshold=1e-2, **kwargs):
        self.kwargs = kwargs
        self.threshold = threshold
        self._substitutor = SubstitutionPredictor(threshold=threshold,
                                                  **kwargs)

    def apply_transformation(self, structure, return_ranked_list=False):
        if not return_ranked_list:
            raise ValueError("SubstitutionPredictorTransformation doesn't"
                             " support returning 1 structure")

        preds = self._substitutor.composition_prediction(
            structure.composition, to_this_composition=False)
        preds.sort(key=lambda x: x['probability'], reverse=True)

        outputs = []
        for pred in preds:
            st = SubstitutionTransformation(pred['substitutions'])
            output = {'structure': st.apply_transformation(structure),
                      'probability': pred['probability'],
                      'threshold': self.threshold, 'substitutions': {}}
            # dictionary keys have to be converted to strings for JSON
            for key, value in pred['substitutions'].items():
                output['substitutions'][str(key)] = str(value)
            outputs.append(output)
        return outputs

    def __str__(self):
        return "SubstitutionPredictorTransformation"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return True


class MagOrderParameterConstraint(MSONable):

    def __init__(self, order_parameter,
                 species_constraints=None,
                 site_constraint_name=None,
                 site_constraints=None):
        """
        This class can be used to supply MagOrderingTransformation
        to just a specific subset of species or sites that satisfy the
        provided constraints. This can be useful for setting an order
        parameters for, for example, ferrimagnetic structures which
        might order on certain motifs, with the global order parameter
        dependent on how many sites satisfy that motif.

        :param order_parameter (float): any number from 0.0 to 1.0,
        typically 0.5 (antiferromagnetic) or 1.0 (ferromagnetic)
        :param species_constraint (list): str or list of strings
        of Specie symbols that the constraint should apply to
        :param site_constraint_name (str): name of the site property
        that the constraint should apply to, e.g. "coordination_no"
        :param site_constraints (list): list of values of the site
        property that the constraints should apply to
        """

        # validation
        if site_constraints and site_constraints != [None] \
                and not site_constraint_name:
                raise ValueError("Specify the name of the site constraint.")
        elif not site_constraints and site_constraint_name:
            raise ValueError("Please specify some site constraints.")
        if not isinstance(species_constraints, list):
            species_constraints = [species_constraints]
        if not isinstance(site_constraints, list):
            site_constraints = [site_constraints]

        if order_parameter > 1 or order_parameter < 0:
            raise ValueError('Order parameter must lie between 0 and 1')
        elif order_parameter != 0.5:
            warnings.warn("Use care when using a non-standard order parameter, "
                          "though it can be useful in some cases it can also "
                          "lead to unintended behavior. Consult documentation.")

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

        if self.species_constraints \
                and str(site.specie) in self.species_constraints:
            satisfies_constraints = True
        else:
            satisfies_constraints = False

        if self.site_constraint_name \
                and self.site_constraint_name in site.properties:
            prop = site.properties[self.site_constraint_name]
            if prop in self.site_constraints:
                satisfies_constraints = True
            else:
                satisfies_constraints = False

        return satisfies_constraints


class MagOrderingTransformation(AbstractTransformation):

    def __init__(self, mag_species_spin, order_parameter=0.5,
                 energy_model=SymmetryModel(), **kwargs):
        """
        This transformation takes a structure and returns a list of collinear
        magnetic orderings. For disordered structures, make an ordered
        approximation first.

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
            order_parameter = [MagOrderParameterConstraint(order_parameter=order_parameter,
                                                           species_constraints=
                                                           list(mag_species_spin.keys()))]
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
                op = max(site.species_and_occu.values())
                # this very hacky bit of code only works because we know
                # that on disordered sites in this class, all species are the same
                # but have different spins, and this is comma-delimited
                sp = str(list(site.species_and_occu.keys())[0]).split(",")[0]
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
        DummySpecie, depending on the constraints specified).
        This approach can also preserve site properties even after
        enumeration.
        """

        dummy_struct = structure.copy()

        def generate_dummy_specie():
            """
            Generator which returns DummySpecie symbols Mma, Mmb, etc.
            """
            subscript_length = 1
            while True:
                for subscript in product(ascii_lowercase, repeat=subscript_length):
                    yield "Mm"+"".join(subscript)
                subscript_length += 1
        dummy_species_gen = generate_dummy_specie()

        # one dummy species for each order parameter constraint
        dummy_species_symbols = [next(dummy_species_gen) for i in range(len(order_parameters))]
        dummy_species = [{
            DummySpecie(symbol, properties={'spin': Spin.up}): constraint.order_parameter,
            DummySpecie(symbol, properties={'spin': Spin.down}): 1-constraint.order_parameter
        } for symbol, constraint in zip(dummy_species_symbols, order_parameters)]

        sites_to_add = []

        for idx, site in enumerate(dummy_struct):
            satisfies_constraints = [c.satisfies_constraint(site) for c in order_parameters]
            if satisfies_constraints.count(True) > 1:
                # site should either not satisfy any constraints, or satisfy
                # one constraint
                raise ValueError("Order parameter constraints conflict for site: {}, {}"
                                 .format(str(site.specie), site.properties))
            elif any(satisfies_constraints):
                dummy_specie_idx = satisfies_constraints.index(True)
                dummy_struct.append(
                    dummy_species[dummy_specie_idx],
                    site.coords,
                    site.lattice
                )

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
        logger.debug('Dummy species structure:\n{}'.format(str(structure)))
        for idx, site in enumerate(structure):
            if isinstance(site.specie, DummySpecie):
                sites_to_remove.append(idx)
                spin = site.specie._properties.get('spin', None)
                neighbors = structure.get_neighbors(
                    site,
                    0.05, # arbitrary threshold, needs to be << any bond length
                    # but >> floating point precision issues
                    include_index=True
                )
                if len(neighbors) != 1:
                    raise Exception("This shouldn't happen, found neighbors: {}"
                                    .format(neighbors))
                orig_site_idx = neighbors[0][2]
                orig_specie = structure[orig_site_idx].specie
                new_specie = Specie(orig_specie.symbol,
                                    getattr(orig_specie, 'oxi_state', None),
                                    properties={'spin': spin})
                structure.replace(orig_site_idx,
                                  new_specie,
                                  properties=structure[orig_site_idx].properties)
        structure.remove_sites(sites_to_remove)
        logger.debug('Structure with dummy species removed:\n{}'.format(str(structure)))
        return structure

    def _add_spin_magnitudes(self, structure):
        """
        Replaces Spin.up/Spin.down with spin magnitudes specified
        by mag_species_spin.
        :param structure:
        :return:
        """
        for idx, site in enumerate(structure):
            if getattr(site.specie, '_properties', None):
                spin = site.specie._properties.get('spin', None)
                sign = int(spin) if spin else 0
                if spin:
                    new_properties = site.specie._properties.copy()
                    # this very hacky bit of code only works because we know
                    # that on disordered sites in this class, all species are the same
                    # but have different spins, and this is comma-delimited
                    sp = str(site.specie).split(",")[0]
                    new_properties.update({
                        'spin': sign*self.mag_species_spin.get(sp, 0)
                    })
                    new_specie = Specie(site.specie.symbol,
                                        getattr(site.specie, 'oxi_state', None),
                                        new_properties)
                    structure.replace(idx, new_specie,
                                      properties=site.properties)
        logger.debug('Structure with spin magnitudes:\n{}'.format(str(structure)))
        return structure

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Apply MagOrderTransformation to an input structure.
        :param structure: Any ordered structure.
        :param return_ranked_list: As in other Transformations.
        :return:
        """

        if not structure.is_ordered:
            raise ValueError("Create an ordered approximation of "
                             "your  input structure first.")

        # retrieve order parameters
        order_parameters = [MagOrderParameterConstraint.from_dict(op_dict)
                            for op_dict in self.order_parameter]
        # add dummy species on which to perform enumeration
        structure = self._add_dummy_species(structure, order_parameters)

        # trivial case
        if structure.is_ordered:
            structure = self._remove_dummy_species(structure)
            return [structure] if return_ranked_list > 1 else structure

        enum_kwargs = self.enum_kwargs.copy()

        enum_kwargs["min_cell_size"] = max(
            int(self.determine_min_cell(structure)),
            enum_kwargs.get("min_cell_size", 1)
        )

        max_cell = enum_kwargs.get('max_cell_size')
        if max_cell:
            if enum_kwargs["min_cell_size"] > max_cell:
                raise ValueError('Specified max cell size is smaller'
                                 ' than the minimum enumerable cell size')
        else:
            enum_kwargs["max_cell_size"] = enum_kwargs["min_cell_size"]

        t = EnumerateStructureTransformation(**enum_kwargs)

        alls = t.apply_transformation(structure,
                                      return_ranked_list=return_ranked_list)

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
        key = lambda x: SpacegroupAnalyzer(x, 0.1).get_space_group_number()
        out = []
        for _, g in groupby(sorted([d["structure"] for d in alls],
                                   key=key), key):
            g = list(g)
            grouped = m.group_structures(g)
            out.extend([{"structure": g[0],
                         "energy": self.energy_model.get_energy(g[0])}
                        for g in grouped])

        self._all_structures = sorted(out, key=lambda d: d["energy"])

        return self._all_structures[0:num_to_return]

    def __str__(self):
        return "MagOrderingTransformation"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return True


def _find_codopant(target, oxidation_state, allowed_elements=None):
    """
    Finds the element from "allowed elements" that (i) possesses the desired
    "oxidation state" and (ii) is closest in ionic radius to the target specie

    Args:
        target: (Specie) provides target ionic radius.
        oxidation_state: (float) codopant oxidation state.
        allowed_elements: ([str]) List of allowed elements. If None,
            all elements are tried.

    Returns:
        (Specie) with oxidation_state that has ionic radius closest to
        target.
    """
    ref_radius = target.ionic_radius
    candidates = []
    symbols = allowed_elements or [el.symbol for el in Element]
    for sym in symbols:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sp = Specie(sym, oxidation_state)
                r = sp.ionic_radius
                if r is not None:
                    candidates.append((r, sp))
        except:
            pass
    return min(candidates, key=lambda l: abs(l[0]/ref_radius - 1))[1]


class DopingTransformation(AbstractTransformation):
    """
    A transformation that performs doping of a structure.
    """

    def __init__(self, dopant, ionic_radius_tol=float("inf"), min_length=10,
                 alio_tol=0, codopant=False, max_structures_per_enum=100,
                 allowed_doping_species=None, **kwargs):
        """
        Args:
            dopant (Specie-like): E.g., Al3+. Must have oxidation state.
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
            \\*\\*kwargs:
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
        logger.info("Composition: %s" % comp)

        for sp in comp:
            try:
                sp.oxi_state
            except AttributeError:
                analyzer = BVAnalyzer()
                structure = analyzer.get_oxi_state_decorated_structure(
                    structure)
                comp = structure.composition
                break

        ox = self.dopant.oxi_state
        radius = self.dopant.ionic_radius

        compatible_species = [
            sp for sp in comp if sp.oxi_state == ox and
            abs(sp.ionic_radius / radius - 1) < self.ionic_radius_tol]

        if (not compatible_species) and self.alio_tol:
            # We only consider aliovalent doping if there are no compatible
            # isovalent species.
            compatible_species = [
                sp for sp in comp
                if abs(sp.oxi_state - ox) <= self.alio_tol and
                abs(sp.ionic_radius / radius - 1) < self.ionic_radius_tol and
                sp.oxi_state * ox >= 0]

        if self.allowed_doping_species is not None:
            # Only keep allowed doping species.
            compatible_species = [
                sp for sp in compatible_species
                if sp in [get_el_sp(s) for s in self.allowed_doping_species]]

        logger.info("Compatible species: %s" % compatible_species)

        lengths = structure.lattice.abc
        scaling = [max(1, int(round(math.ceil(self.min_length/x))))
                   for x in lengths]
        logger.info("Lengths are %s" % str(lengths))
        logger.info("Scaling = %s" % str(scaling))

        all_structures = []
        t = EnumerateStructureTransformation(**self.kwargs)

        for sp in compatible_species:
            supercell = structure * scaling
            nsp = supercell.composition[sp]
            if sp.oxi_state == ox:
                supercell.replace_species({sp: {sp: (nsp - 1)/nsp,
                                                self.dopant: 1/nsp}})
                logger.info("Doping %s for %s at level %.3f" % (
                    sp, self.dopant, 1 / nsp))
            elif self.codopant:
                codopant = _find_codopant(sp, 2 * sp.oxi_state - ox)
                supercell.replace_species({sp: {sp: (nsp - 2) / nsp,
                                                self.dopant: 1 / nsp,
                                                codopant: 1 / nsp}})
                logger.info("Doping %s for %s + %s at level %.3f" % (
                    sp, self.dopant, codopant, 1 / nsp))
            elif abs(sp.oxi_state) < abs(ox):
                # Strategy: replace the target species with a
                # combination of dopant and vacancy.
                # We will choose the lowest oxidation state species as a
                # vacancy compensation species as it is likely to be lower in
                # energy
                sp_to_remove = min([s for s in comp if s.oxi_state * ox > 0],
                                    key=lambda ss: abs(ss.oxi_state))

                if sp_to_remove == sp:
                    common_charge = lcm(int(abs(sp.oxi_state)), int(abs(ox)))
                    ndopant = common_charge / abs(ox)
                    nsp_to_remove = common_charge / abs(sp.oxi_state)
                    logger.info("Doping %d %s with %d %s." %
                                (nsp_to_remove, sp, ndopant, self.dopant))
                    supercell.replace_species(
                        {sp: {sp: (nsp - nsp_to_remove) / nsp,
                              self.dopant: ndopant / nsp}})
                else:
                    ox_diff = int(abs(round(sp.oxi_state - ox)))
                    vac_ox = int(abs(sp_to_remove.oxi_state))
                    common_charge = lcm(vac_ox, ox_diff)
                    ndopant = common_charge / ox_diff
                    nx_to_remove = common_charge / vac_ox
                    nx = supercell.composition[sp_to_remove]
                    logger.info("Doping %d %s with %s and removing %d %s." %
                                (ndopant, sp, self.dopant,
                                 nx_to_remove, sp_to_remove))
                    supercell.replace_species(
                        {sp: {sp: (nsp - ndopant) / nsp,
                              self.dopant: ndopant / nsp},
                         sp_to_remove: {
                             sp_to_remove: (nx - nx_to_remove) / nx}})
            elif abs(sp.oxi_state) > abs(ox):
                # Strategy: replace the target species with dopant and also
                # remove some opposite charged species for charge neutrality
                if ox > 0:
                    sp_to_remove = max(supercell.composition.keys(),
                                       key=lambda el: el.X)
                else:
                    sp_to_remove = min(supercell.composition.keys(),
                                       key=lambda el: el.X)
                # Confirm species are of opposite oxidation states.
                assert sp_to_remove.oxi_state * sp.oxi_state < 0

                ox_diff = int(abs(round(sp.oxi_state - ox)))
                anion_ox = int(abs(sp_to_remove.oxi_state))
                nx = supercell.composition[sp_to_remove]
                common_charge = lcm(anion_ox, ox_diff)
                ndopant = common_charge / ox_diff
                nx_to_remove = common_charge / anion_ox
                logger.info("Doping %d %s with %s and removing %d %s." %
                            (ndopant, sp, self.dopant,
                             nx_to_remove, sp_to_remove))
                supercell.replace_species(
                    {sp: {sp: (nsp - ndopant) / nsp,
                          self.dopant: ndopant / nsp},
                     sp_to_remove: {sp_to_remove: (nx - nx_to_remove)/nx}})

            ss = t.apply_transformation(
                supercell, return_ranked_list=self.max_structures_per_enum)
            logger.info("%s distinct structures" % len(ss))
            all_structures.extend(ss)

        logger.info("Total %s doped structures" % len(all_structures))
        if return_ranked_list:
            return all_structures[:return_ranked_list]

        return all_structures[0]["structure"]

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return True


class SlabTransformation(AbstractTransformation):
    """
    A transformation that creates a slab from a structure.

    """
    def __init__(self, miller_index, min_slab_size, min_vacuum_size,
                 lll_reduce=False, center_slab=False, primitive=True,
                 max_normal_search=None, shift=0, tol=0.1):
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
        self.primitive = primitive
        self.max_normal_search = max_normal_search
        self.shift = shift
        self.tol = 0.1

    def apply_transformation(self, structure):
        sg = SlabGenerator(structure, self.miller_index, self.min_slab_size,
                           self.min_vacuum_size, self.lll_reduce,
                           self.center_slab, self.primitive,
                           self.max_normal_search)
        slab = sg.get_slab(self.shift, self.tol)
        return slab

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return None
