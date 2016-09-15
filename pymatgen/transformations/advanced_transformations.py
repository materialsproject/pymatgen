# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import numpy as np
from fractions import gcd, Fraction
from itertools import groupby
from warnings import warn
import logging
import math

import six
from monty.json import MontyDecoder
from monty.fractions import lcm

from pymatgen.core.structure import Composition
from pymatgen.core.periodic_table import Element, Specie, get_el_sp
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.transformations.standard_transformations import \
    SubstitutionTransformation, OrderDisorderedStructureTransformation
from pymatgen.command_line.enumlib_caller import EnumlibAdaptor
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.structure_prediction.substitution_probability import \
    SubstitutionPredictor
from pymatgen.analysis.structure_matcher import StructureMatcher, \
    SpinComparator
from pymatgen.analysis.energy_models import SymmetryModel
from pymatgen.analysis.bond_valence import BVAnalyzer


"""
This module implements more advanced transformations.
"""

__author__ = "Shyue Ping Ong, Stephen Dacek, Anubhav Jain"
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
    """

    def __init__(self, min_cell_size=1, max_cell_size=1, symm_prec=0.1,
                 refine_structure=False, enum_precision_parameter=0.001,
                 check_ordered_symmetry=True):
        self.symm_prec = symm_prec
        self.min_cell_size = min_cell_size
        self.max_cell_size = max_cell_size
        self.refine_structure = refine_structure
        self.enum_precision_parameter = enum_precision_parameter
        self.check_ordered_symmetry = check_ordered_symmetry

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

        if structure.is_ordered:
            warn("Enumeration skipped for structure with composition {} "
                 "because it is ordered".format(structure.composition))
            structures = [structure.copy()]
        else:
            adaptor = EnumlibAdaptor(
                structure, min_cell_size=self.min_cell_size,
                max_cell_size=self.max_cell_size,
                symm_prec=self.symm_prec, refine_structure=False,
                enum_precision_parameter=self.enum_precision_parameter,
                check_ordered_symmetry=self.check_ordered_symmetry)
            adaptor.run()
            structures = adaptor.structures

        original_latt = structure.lattice
        inv_latt = np.linalg.inv(original_latt.matrix)
        ewald_matrices = {}
        all_structures = []
        for s in structures:
            new_latt = s.lattice
            transformation = np.dot(new_latt.matrix, inv_latt)
            transformation = tuple([tuple([int(round(cell)) for cell in row])
                                    for row in transformation])
            if contains_oxidation_state:
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
            return s["energy"] / s["num_sites"] if contains_oxidation_state \
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


class MagOrderingTransformation(AbstractTransformation):
    """
    This transformation takes a structure and returns a list of magnetic
    orderings. Currently only works for ordered structures.

    Args:
        mag_elements_spin:
            A mapping of elements/species to magnetically order to spin
            magnitudes. E.g., {"Fe3+": 5, "Mn3+": 4}
        order_parameter:
            degree of magnetization. 0.5 corresponds to
            antiferromagnetic order
        energy_model:
            Energy model used to rank the structures. Some models are
            provided in :mod:`pymatgen.analysis.energy_models`.
        **kwargs:
            Same keyword args as :class:`EnumerateStructureTransformation`,
            i.e., min_cell_size, etc.
    """

    def __init__(self, mag_species_spin, order_parameter=0.5,
                 energy_model=SymmetryModel(), **kwargs):
        self.mag_species_spin = mag_species_spin
        if order_parameter > 1 or order_parameter < 0:
            raise ValueError('Order Parameter must lie between 0 and 1')
        else:
            self.order_parameter = order_parameter
        self.energy_model = energy_model
        self.kwargs = kwargs

    @classmethod
    def determine_min_cell(cls, structure, mag_species_spin, order_parameter):
        """
        Determine the smallest supercell that is able to enumerate
        the provided structure with the given order parameter
        """

        def lcm(n1, n2):
            """
            Find least common multiple of two numbers
            """
            return n1 * n2 / gcd(n1, n2)

        denom = Fraction(order_parameter).limit_denominator(100).denominator
        atom_per_specie = [structure.composition[m]
                           for m in mag_species_spin.keys()]
        n_gcd = six.moves.reduce(gcd, atom_per_specie)

        if not n_gcd:
            raise ValueError(
                'The specified species do not exist in the structure'
                ' to be enumerated')

        return lcm(n_gcd, denom) / n_gcd

    def apply_transformation(self, structure, return_ranked_list=False):
        # Make a mutable structure first
        mods = Structure.from_sites(structure)
        for sp, spin in self.mag_species_spin.items():
            sp = get_el_sp(sp)
            oxi_state = getattr(sp, "oxi_state", 0)
            if spin:
                up = Specie(sp.symbol, oxi_state, {"spin": abs(spin)})
                down = Specie(sp.symbol, oxi_state, {"spin": -abs(spin)})
                mods.replace_species(
                    {sp: Composition({up: self.order_parameter,
                                      down: 1 - self.order_parameter})})
            else:
                mods.replace_species(
                    {sp: Specie(sp.symbol, oxi_state, {"spin": spin})})

        if mods.is_ordered:
            return [mods] if return_ranked_list > 1 else mods

        enum_args = self.kwargs

        enum_args["min_cell_size"] = max(int(
            MagOrderingTransformation.determine_min_cell(
                structure, self.mag_species_spin,
                self.order_parameter)),
            enum_args.get("min_cell_size", 1))

        max_cell = enum_args.get('max_cell_size')
        if max_cell:
            if enum_args["min_cell_size"] > max_cell:
                raise ValueError('Specified max cell size is smaller'
                                 ' than the minimum enumerable cell size')
        else:
            enum_args["max_cell_size"] = enum_args["min_cell_size"]

        t = EnumerateStructureTransformation(**enum_args)

        alls = t.apply_transformation(mods,
                                      return_ranked_list=return_ranked_list)

        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        if num_to_return == 1 or not return_ranked_list:
            return alls[0]["structure"] if num_to_return else alls

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
            \*\*kwargs:
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
