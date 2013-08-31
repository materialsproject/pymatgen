#!/usr/bin/env python

"""
This module implements more advanced transformations.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 24, 2012"

import numpy as np

from pymatgen.core.periodic_table import smart_element_or_specie
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.transformations.standard_transformations import \
    SubstitutionTransformation, OrderDisorderedStructureTransformation
from pymatgen.command_line.enumlib_caller import EnumlibAdaptor
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.core.structure import Structure
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.structure_prediction.substitution_probability import \
    SubstitutionPredictor


class ChargeBalanceTransformation(AbstractTransformation):
    """
    This is a transformation that disorders a structure to make it charge
    balanced, given an oxidation state-decorated structure.
    """

    def __init__(self, charge_balance_sp):
        """
        Args:
            charge_balance_sp
                specie to add or remove. Currently only removal is supported
        """
        self._charge_balance_sp = str(charge_balance_sp)

    def apply_transformation(self, structure):
        charge = structure.charge
        specie = smart_element_or_specie(self._charge_balance_sp)
        num_to_remove = charge / specie.oxi_state
        num_in_structure = structure.composition[specie]
        removal_fraction = num_to_remove / num_in_structure
        if removal_fraction < 0:
            raise ValueError("addition of specie not yet supported by "
                             "ChargeBalanceTransformation")
        trans = SubstitutionTransformation({self._charge_balance_sp:
                                            {self._charge_balance_sp:
                                             1 - removal_fraction}})
        return trans.apply_transformation(structure)

    def __str__(self):
        return "Charge Balance Transformation : " + \
            "Species to remove = {}".format(str(self._charge_balance_sp))

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
                "init_args": {"charge_balance_sp": self._charge_balance_sp},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class SuperTransformation(AbstractTransformation):
    """
    This is a transformation that is inherently one-to-many. It is constructed
    from a list of transformations and returns one structure for each
    transformation. The primary use for this class is extending a transmuter
    object.
    """

    def __init__(self, transformations):
        """
        Args:
            transformations:
                list of transformations to apply to a structure. One
                transformation is applied to each output structure.
        """
        self._transformations = transformations

    def apply_transformation(self, structure, return_ranked_list=False):
        if not return_ranked_list:
            raise ValueError("SuperTransformation has no single best structure"
                             " output. Must use return_ranked_list")
        structures = []
        for t in self._transformations:
            structures.append({"transformation": t,
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

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"transformations": self._transformations},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


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
            sp_to_replace
                species to be replaced
            r_fraction
                fraction of that specie to replace
            substitution_dict
                dictionary of the format
                {2: ["Mg", "Ti", "V", "As", "Cr", "Ta", "N", "Nb"],
                3: ["Ru", "Fe", "Co", "Ce", "As", "Cr", "Ta", "N", "Nb"],
                4: ["Ru", "V", "Cr", "Ta", "N", "Nb"],
                5: ["Ru", "W", "Mn"]
                }
                The number is the charge used for each of the list of elements
                (an element can be present in multiple lists)
            charge_balance_species:
                If specified, will balance the charge on the structure using
                that specie.
        """
        self._sp_to_replace = sp_to_replace
        self._r_fraction = r_fraction
        self._substitution_dict = substitution_dict
        self._charge_balance_species = charge_balance_species
        self._order = order

    def apply_transformation(self, structure, return_ranked_list=False):
        if not return_ranked_list:
            raise ValueError("MultipleSubstitutionTransformation has no single"
                             " best structure output. Must use"
                             " return_ranked_list.")
        outputs = []
        for charge, el_list in self._substitution_dict.items():
            mapping = {}
            if charge > 0:
                sign = "+"
            else:
                sign = "-"
            dummy_sp = "X{}{}".format(str(charge), sign)
            mapping[self._sp_to_replace] = {self._sp_to_replace:
                                            1 - self._r_fraction,
                                            dummy_sp: self._r_fraction}
            trans = SubstitutionTransformation(mapping)
            dummy_structure = trans.apply_transformation(structure)
            if self._charge_balance_species is not None:
                cbt = ChargeBalanceTransformation(self._charge_balance_species)
                dummy_structure = cbt.apply_transformation(dummy_structure)
            if self._order:
                trans = OrderDisorderedStructureTransformation()
                dummy_structure = trans.apply_transformation(dummy_structure)

            for el in el_list:
                if charge > 0:
                    sign = "+"
                else:
                    sign = "-"
                st = SubstitutionTransformation({"X{}+".format(str(charge)):
                                                 "{}{}{}".format(el, charge,
                                                                 sign)})
                new_structure = st.apply_transformation(dummy_structure)
                outputs.append({"structure": new_structure})
        return outputs

    def __str__(self):
        return "Multiple Substitution Transformation : Substitution on " + \
               "{}".format(self._sp_to_replace)

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
                "init_args": {"sp_to_replace": self._sp_to_replace,
                              "r_fraction": self._r_fraction,
                              "substitution_dict": self._substitution_dict,
                              "charge_balance_species":
                              self._charge_balance_species},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class EnumerateStructureTransformation(AbstractTransformation):
    """
    Order a disordered structure using enumlib. For complete orderings, this
    generally produces fewer structures that the OrderDisorderedStructure
    transformation, and at a much faster speed.
    """

    def __init__(self, min_cell_size=1, max_cell_size=1, symm_prec=0.1,
                 refine_structure=False):
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
                If you are already starting from an experimental cif, refinment
                should have already been done and it is not necessary. Defaults
                to False.
        """
        self.symm_prec = symm_prec
        self.min_cell_size = min_cell_size
        self.max_cell_size = max_cell_size
        self.refine_structure = refine_structure

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Return either a single ordered structure or a sequence of all ordered
        structures.

        Args:
            structure:
                Structure to order.
            return_ranked_list:
                Boolean stating whether or not multiple structures are
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

        if structure.is_ordered:
            raise ValueError("Enumeration can be carried out only on "
                             "disordered structures!")

        if self.refine_structure:
            finder = SymmetryFinder(structure, self.symm_prec)
            structure = finder.get_refined_structure()

        contains_oxidation_state = True
        for sp in structure.composition.elements:
            if not hasattr(sp, "oxi_state"):
                contains_oxidation_state = False
                break

        adaptor = EnumlibAdaptor(structure, min_cell_size=self.min_cell_size,
                                 max_cell_size=self.max_cell_size,
                                 symm_prec=self.symm_prec,
                                 refine_structure=False)
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
                    s_supercell = Structure.from_sites(structure.sites)
                    s_supercell.make_supercell(transformation)
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

    @property
    def to_dict(self):
        return {"name": self.__class__.__name__, "version": __version__,
                "init_args": {"symm_prec": self.symm_prec,
                              "min_cell_size": self.min_cell_size,
                              "max_cell_size": self.max_cell_size,
                              "refine_structure": self.refine_structure},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class SubstitutionPredictorTransformation(AbstractTransformation):
    """
    This transformation takes a structure and uses the structure
    prediction module to find likely site substitutions.
    """

    def __init__(self, threshold=1e-2, **kwargs):
        """
        Args:
            kwargs:
                args for SubstitutionProbability class
                lambda_table, alpha
        """
        self._kwargs = kwargs
        self._threshold = threshold
        self._substitutor = SubstitutionPredictor(threshold=threshold, **kwargs)

    def apply_transformation(self, structure, return_ranked_list=False):
        if not return_ranked_list:
            raise ValueError("SubstitutionPredictorTransformation doesn't"
                             " support returning 1 structure")

        preds = self._substitutor.composition_prediction(structure.composition,
                                                         to_this_composition=False)
        preds.sort(key=lambda x: x['probability'], reverse=True)

        outputs = []
        for pred in preds:
            st = SubstitutionTransformation(pred['substitutions'])
            output = {'structure': st.apply_transformation(structure),
                      'probability': pred['probability'],
                      'threshold': self._threshold, 'substitutions': {}}
            #dictionary keys have to be converted to strings for JSON
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

    @property
    def to_dict(self):
        d = {"name": self.__class__.__name__, "version": __version__,
             "init_args": self._kwargs, "@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d["init_args"]["threshold"] = self._threshold
        return d
