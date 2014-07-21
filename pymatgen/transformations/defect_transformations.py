#!/usr/bin/env python

"""
This module defines classes for point defect transformations on structures
"""
from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "0.1"
__maintainier__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__date__ = "Jul 1 2014"


from pymatgen.core.structure import Structure
from pymatgen.analysis.defects.point_defects import Vacancy, \
    ValenceIonicRadiusEvaluator, Interstitial
from pymatgen.transformations.transformation_abc import AbstractTransformation


class VacancyTransformation(AbstractTransformation):
    """
    Generates vacancy structures
    """
    def __init__(self, supercell_dim, valences=None, radii=None):
        """
        :param supecell_dim: Supercell scaling matrix
        :return:
        """
        self.supercell_dim = supercell_dim
        self.valences = valences
        self.radii = radii

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        :param structure:
        :param return_ranked_list (Logical or integer): Use big enough
         number to return all defect structures
        :return:
            scs: Supercells with one vacancy in each structure.
        """
        if not return_ranked_list:
            raise ValueError("VacancyTransformation has no single best structure"
                             " output. Must use return_ranked_list.")
        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1
        vac = Vacancy(structure,self.valences,self.radii)
        #print vac.enumerate_defectsites()
        scs = vac.make_supercells_with_defects(self.supercell_dim)
        #if num_to_return < len(scs)-1:
        #    raise ValueError("VacancyTransformation has no ordering of best "
        #                     "structure. Must increase return_ranked_list.")
        structures = []
        num_to_return = min(num_to_return, len(scs)-1)
        for sc in scs[1:num_to_return+1]: 
            structures.append({'structure':sc})
        return structures

    def __str__(self):
        inp_args = ["Supercell scaling matrix = {}".format(self.supercell_dim),
                    "Valences of ions = {}".format(self.valences),
                    "Radii of ions = {}".format(self.radii)]
        return "Vacancy Transformation : " + ", ".join(inp_args)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        pass

    @property
    def is_one_to_many(self):
        return True

    @property
    def to_dict(self):
        return {"name":self.__class__.__name__, "version":__version__,
                "init_args":{"supercell_dim":self.supercell_dim,
                             "valences":self.valences,"radii":self.radii},
                "@module":self.__class__.__module__,
                "@class":self.__class__.__name__ }


class SubstitutionDefectTransformation(AbstractTransformation):
    """
    Generates substiutional defect structures.
    The first structure is the supercell of the original structure
    and is not a defect structure.
    """
    def __init__(self, specie_map, supercell_dim,
                valences=None, radii=None):
        """
        :param supecell_dim: Supercell scaling matrix
        :return:
        """
        #self.substitute_specie = substitute_specie
        #self.site_specie = site_specie
        self._specie_map = specie_map
        self.supercell_dim = supercell_dim
        self.valences = valences
        self.radii = radii

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        :param structure:
        :param return_ranked_list (Logical or integer): Use big enough
         number to return all defect structures
        :return:
            scs: Supercells with one substitution defect in each structure.
        """
        if not return_ranked_list:
            raise ValueError("SubstitutionDefectTransformation has no single" 
                    "best structure output. Must use return_ranked_list.")
        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        vac = Vacancy(structure,self.valences,self.radii)
        scs = vac.make_supercells_with_defects(self.supercell_dim)
        blk_sc = scs[0]
        sub_scs = []
        for i in range(1,len(scs)):
            vac_sc = scs[i]
            vac_site = list(set(blk_sc.sites) - set(vac_sc.sites))[0]
            site_specie = str(vac_site.specie)
            if site_specie in self._specie_map.keys():
                substitute_specie = self._specie_map[site_specie]
                vac_sc.append(substitute_specie, vac_site.frac_coords)
                sub_scs.append(vac_sc.get_sorted_structure())

        #if num_to_return < len(sub_scs)-1:
        #    raise ValueError("SubstitutionDefectTransformation has no ordering"
        #            " of best structure. Must increase return_ranked_list.")
        num_to_return = min(num_to_return, len(sub_scs))

        structures = []
        for sc in sub_scs[0:num_to_return]:
            structures.append({'structure':sc})
        return structures

    def __str__(self):
        inp_args = ["Specie map = {}".format(self._specie_map),
                    "Supercell scaling matrix = {}".format(self.supercell_dim),
                    "Valences of ions = {}".format(self.valences),
                    "Radii of ions = {}".format(self.radii)]
        return "Substitutional Defect Transformation : " + ", ".join(inp_args)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        pass

    @property
    def is_one_to_many(self):
        return True

    @property
    def to_dict(self):
        return {"name":self.__class__.__name__, "version":__version__,
                "init_args":{"specie_map":self._specie_map,
                             "supercell_dim":self.supercell_dim,
                             "valences":self.valences,"radii":self.radii},
                "@module":self.__class__.__module__,
                "@class":self.__class__.__name__ }


class AntisiteDefectTransformation(AbstractTransformation):
    """
    Generates antisite defect structures
    """
    def __init__(self, supercell_dim, valences=None, radii=None):
        """
        :param supecell_dim: Supercell scaling matrix
        :return:
        """
        self.supercell_dim = supercell_dim
        self.valences = valences
        self.radii = radii

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        :param structure:
        :param return_ranked_list (Logical or integer): Use big enough
         number to return all defect structures
        :return:
            scs: Supercells with one antisite defect in each structure. 
        """
        if not return_ranked_list:
            raise ValueError("AntisiteDefectTransformation has no single best"
                    "structure output. Must use return_ranked_list.")
        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        vac = Vacancy(structure,self.valences,self.radii)
        scs = vac.make_supercells_with_defects(self.supercell_dim)
        blk_sc = scs[0]
        as_scs = []
        struct_species = blk_sc.types_of_specie
        for i in range(1,len(scs)):
            vac_sc = scs[i]
            vac_site = list(set(blk_sc.sites) - set(vac_sc.sites))[0]
            for specie in set(struct_species) - set([vac_site.specie]):
                anti_struct = vac_sc.copy()
                anti_struct.append(specie, vac_site.frac_coords)
                as_scs.append(anti_struct.get_sorted_structure())

        #if num_to_return < len(as_scs)-1:
        #    raise ValueError("AntisiteDefectTransformation has no ordering "
        #            "of best structures. Must increase return_ranked_list.")
        num_to_return = min(num_to_return,len(as_scs))
        structures = []
        for sc in as_scs[0:num_to_return]:
            structures.append({'structure':sc})
        return structures

    def __str__(self):
        inp_args = ["Supercell scaling matrix = {}".format(self.supercell_dim),
                    "Valences of ions = {}".format(self.valences),
                    "Radii of ions = {}".format(self.radii)]
        return "Antisite Defect Transformation : " + ", ".join(inp_args)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        pass

    @property
    def is_one_to_many(self):
        return True

    @property
    def to_dict(self):
        return {"name":self.__class__.__name__, "version":__version__,
                "init_args":{"supercell_dim":self.supercell_dim,
                             "valences":self.valences,"radii":self.radii},
                "@module":self.__class__.__module__,
                "@class":self.__class__.__name__ }


class InterstitialTransformation(AbstractTransformation):
    """
    Generates interstitial structures from the input structure
    """
    def __init__(self, interstitial_specie, supercell_dim,
                 valences={}, radii={}):
        """
        :param supercell_dim:
        :param valences:
        :param radii:
        :return:
        """
        self.supercell_dim = supercell_dim
        self.valences = valences
        self.radii = radii
        self.inter_specie = interstitial_specie

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        :param structure:
        :param return_ranked_list (Logical or integer): Use big enough
         number to return all defect structures
        :return:
            scs: Supercells with one interstitial defect in each structure. 
        """
        if not return_ranked_list:
            raise ValueError("InterstitialTransformation has no single best "
                    "structure output. Must use return_ranked_list.")
        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        if self.radii:
            inter = Interstitial(structure, self.valences, self.radii)
        else:
            s = structure.copy()
            valrad_eval = ValenceIonicRadiusEvaluator(s)
            s = valrad_eval.structure 
            val = valrad_eval.valences
            rad = valrad_eval.radii
            inter = Interstitial(s,val,rad,oxi_state=True)

        scs = inter.make_supercells_with_defects(
            self.supercell_dim, self.inter_specie)

        #if num_to_return < len(scs)-1:
        #    raise ValueError("InterstitialTransformation has no ordering "
        #            "of best structures. Must increase return_ranked_list.")

        structures = []
        num_to_return = min(num_to_return,len(scs)-1)
        for sc in scs[1:num_to_return+1]:
            structures.append({'structure':sc})
        return structures

    def __str__(self):
        inp_args = ["Supercell scaling matrix = {}".format(self.supercell_dim),
                    "Valences of ions = {}".format(self.valences),
                    "Radii of ions = {}".format(self.radii),
                    "interstitial specie = {}".format(self.inter_specie)]
        return "Interstitial Transformation : " + ", ".join(inp_args)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        pass

    @property
    def is_one_to_many(self):
        return True

    @property
    def to_dict(self):
        return {"name":self.__class__.__name__, "version":__version__,
                "init_args":{"supercell_dim":self.supercell_dim,
                             "valences":self.valences,"radii":self.radii,
                             "interstitial_specie":self.inter_specie},
                "@module":self.__class__.__module__,
                "@class":self.__class__.__name__ }
