# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
"""
This module defines classes to define point defect objects
"""

__author__ = "Danny Broberg, Shyam Dwaraknath, Bharat Medasani, Nils E. R. Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"

import six
import copy
import logging
import numpy as np

from abc import ABCMeta, abstractmethod
from monty.json import MSONable
from monty.functools import lru_cache

from pymatgen.core.composition import Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

logger = logging.getLogger(__name__)


class Defect(six.with_metaclass(ABCMeta, MSONable)):
    """
    Abstract class for a single point defect
    """

    def __init__(self, structure, defect_site, charge=0.):
        """
        Initializes an abstract defect

        Args:
            structure: Pymatgen Structure without any defects
            charge: (int or float) defect charge
                default is zero, meaning no change to NELECT after defect is created in the structure
        """
        self._structure = structure
        self._charge = charge
        self._defect_site = defect_site

    @property
    def bulk_structure(self):
        """
        Returns the structure without any defects.
        """
        return self._structure

    @property
    def charge(self):
        """
        Returns the charge of a defect
        """
        return self._charge

    @property
    def site(self):
        """
        Returns the defect position as a site object
        """
        return self._defect_site

    @property
    @abstractmethod
    def multiplicity(self):
        """
        Returns the multiplicity of a defect site within the structure (needed for concentration analysis)
        """
        return

    @property
    @abstractmethod
    def defect_composition(self):
        """
        Returns the defect composition as a Composition object
        """
        return

    @abstractmethod
    def generate_defect_structure(self, supercell=(1, 1, 1)):
        """
        Given structure and defect_site (and type of defect) should return a defect_structure that is charged
        Args:
            supercell (int, [3x1], or [[]] (3x3)): supercell integer, vector, or scaling matrix
        """
        return

    @property
    @abstractmethod
    def name(self):
        """
        Returns a name for this defect
        """
        return

    def copy(self):
        """
        Convenience method to get a copy of the defect.

        Returns:
            A copy of the Defect.
        """
        struc = self._structure.copy()
        def_site = copy.copy(self._defect_site)
        charge = copy.copy(self.charge)
        return Defect(struc, def_site, charge=charge)

    def set_charge(self,new_charge=0.):
        """
        Sets the overall charge
        Args:
            charge (float): new charge to set
        """
        self._charge = new_charge

class Vacancy(Defect):
    """
    Subclass of Defect to capture essential information for a single Vacancy defect structure.
    """

    @property
    def defect_composition(self):
        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] -= 1
        return Composition(temp_comp)

    def generate_defect_structure(self, supercell=(1, 1, 1)):
        """
        Returns Defective Vacancy structure, decorated with charge
        Args:
            supercell (int, [3x1], or [[]] (3x3)): supercell integer, vector, or scaling matrix
        """
        defect_structure = self.bulk_structure.copy()
        defect_structure.make_supercell(supercell)
        poss_deflist = sorted(
            defect_structure.get_sites_in_sphere(self.site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]
        defect_structure.remove_sites([defindex])
        defect_structure.set_charge(self.charge)
        return defect_structure

    @property
    def multiplicity(self):
        """
        Returns the multiplicity of a defect site within the structure (needed for concentration analysis)
        """
        sga = SpacegroupAnalyzer(self.bulk_structure)
        periodic_struc = sga.get_symmetrized_structure()
        poss_deflist = sorted(
            periodic_struc.get_sites_in_sphere(self.site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]

        equivalent_sites = periodic_struc.find_equivalent_sites(self.bulk_structure[defindex])
        return len(equivalent_sites)

    @property
    def name(self):
        """
        Returns a name for this defect
        """
        return "Vac_{}_{}".format(self.site.specie,self.multiplicity)

    def copy(self):
        """
        Convenience method to get a copy of the Vacancy.

        Returns:
            A copy of the Vacancy.
        """
        struc = self._structure.copy()
        def_site = copy.copy(self._defect_site)
        charge = copy.copy(self.charge)
        return Vacancy(struc, def_site, charge=charge)


class Substitution(Defect):
    """
    Subclass of Defect to capture essential information for a single Substitution defect structure.
    """

    @property
    @lru_cache(1)
    def defect_composition(self):
        poss_deflist = sorted(
            self.bulk_structure.get_sites_in_sphere(self.site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]

        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] += 1
        temp_comp[str(self.bulk_structure[defindex].specie)] -= 1
        return Composition(temp_comp)

    def generate_defect_structure(self, supercell=(1, 1, 1)):
        """
        Returns Defective Substitution structure, decorated with charge
        Args:
            supercell (int, [3x1], or [[]] (3x3)): supercell integer, vector, or scaling matrix
        """
        defect_structure = self.bulk_structure.copy()
        defect_structure.make_supercell(supercell)
        poss_deflist = sorted(
            defect_structure.get_sites_in_sphere(self.site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]

        subsite = defect_structure.pop(defindex)
        defect_structure.append(self.site.specie.symbol, subsite.coords, coords_are_cartesian=True)
        defect_structure.set_charge(self.charge)
        return defect_structure

    @property
    def multiplicity(self):
        """
        Returns the multiplicity of a defect site within the structure (needed for concentration analysis)
        """
        sga = SpacegroupAnalyzer(self.bulk_structure)
        periodic_struc = sga.get_symmetrized_structure()
        poss_deflist = sorted(
            periodic_struc.get_sites_in_sphere(self.site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]

        equivalent_sites = periodic_struc.find_equivalent_sites(self.bulk_structure[defindex])
        return len(equivalent_sites)

    @property
    @lru_cache(1)
    def name(self):
        """
        Returns a name for this defect
        """
        poss_deflist = sorted(
            self.bulk_structure.get_sites_in_sphere(self.site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]
        return "Sub_{}_on_{}_{}".format(self.site.specie,self.bulk_structure[defindex].specie,self.multiplicity)

    def copy(self):
        """
        Convenience method to get a copy of the Substitution.

        Returns:
            A copy of the Substitution.
        """
        struc = self._structure.copy()
        def_site = copy.copy(self._defect_site)
        charge = copy.copy(self.charge)
        return Substitution(struc, def_site, charge=charge)

class Interstitial(Defect):
    """
    Subclass of Defect to capture essential information for a single Interstitial defect structure.
    """

    def __init__(self, structure, defect_site, charge=0.,multiplicity=1):
        """
        Initializes an interstial defect. 
        User must specify multiplity. Default is 1
        Args:
            structure: Pymatgen Structure without any defects
            charge: (int or float) defect charge
                default is zero, meaning no change to NELECT after defect is created in the structure
            defect_site (Site): the site for the interstial
            multiplicity (int): multiplicity
        """
        super().__init__(structure=structure,defect_site=defect_site,charge=charge)
        self._multiplicity = multiplicity

    @property
    def defect_composition(self):
        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] += 1
        return Composition(temp_comp)

    def generate_defect_structure(self, supercell=(1, 1, 1)):
        """
        Returns Defective Interstitial structure, decorated with charge
        Args:
            supercell (int, [3x1], or [[]] (3x3)): supercell integer, vector, or scaling matrix
        """
        defect_structure = self.bulk_structure.copy()
        defect_structure.make_supercell(supercell)
        defect_structure.append(self.site.specie.symbol, self.site.coords, coords_are_cartesian=True)
        defect_structure.set_charge(self.charge)
        return defect_structure

    @property
    def multiplicity(self):
        """
        Returns the multiplicity of a defect site within the structure (needed for concentration analysis)
        """
        return self._multiplicity

    @property
    def name(self):
        """
        Returns a name for this defect
        """
        return "Int_{}_{}".format(self.site.specie,self.multiplicity)

    def copy(self):
        """
        Convenience method to get a copy of the Interstitial.

        Returns:
            A copy of the Interstitial.
        """
        struc = self._structure.copy()
        def_site = copy.copy(self._defect_site)
        charge = copy.copy(self.charge)
        multiplicity = copy.copy(self._multiplicity)
        return Interstitial(struc, def_site, charge=charge, multiplicity=multiplicity)


class DefectEntry(MSONable):
    """
    An lightweight DefectEntry object containing key computed data
    for many defect analysis.
    """

    def __init__(self, defect, uncorrected_energy, corrections={}, parameters={}, entry_id=None):
        """
        Args:
            defect:
                A Defect object from pymatgen.analysis.defects.core
            uncorrected_energy (float): Energy of the defect entry. Usually the difference between
                the final calculated energy for the defect supercell - the perfect 
                supercell energy
            corrections ([Correction]):
                List of Correction classes (from pymatgen.analysis.defects.corrections)
                which correct energy due to charge (e.g. Freysoldt or Kumagai)
                or other factors (e.g. Shallow level shifts)
            parameters (dict): An optional dict of calculation parameters and data to 
                use with correction schemes
                (examples of parameter keys: supercell_size, axis_grid, bulk_planar_averages
                defect_planar_averages )
            entry_id (obj): An id to uniquely identify this defect, can be any MSONable
                type
        """
        self.defect = defect
        self.uncorrected_energy = uncorrected_energy
        self.corrections = corrections
        self.entry_id = entry_id
        self.parameters = parameters

    @property
    def bulk_structure(self):
        return self.defect.bulk_structure

    @property
    def bulk_sc_structure(self):
        bulk_sc = self.defect.bulk_structure.copy()
        if 'scaling_matrix' in self.parameters:
            bulk_sc.make_supercell( self.parameters['scaling_matrix'])
        return bulk_sc

    @property
    def defect_sc_structure(self):
        if 'scaling_matrix' in self.parameters:
            scaling_matrix = self.parameters['scaling_matrix']
        else:
            scaling_matrix = (1,1,1)
        defect_sc = self.defect.generate_defect_structure(scaling_matrix)
        return defect_sc

    @property
    def site(self):
        return self.defect.site

    @property
    def multiplicty(self):
        sc_multiplier = 1.
        if 'scaling_matrix' in self.parameters:
            for val in np.array(self.parameters['scaling_matrix']).flatten():
                if val:
                    #scale with non-zero scaling matrix entries
                    sc_multiplier *= val
        return self.defect.multiplicty * sc_multiplier

    @property
    def charge(self):
        return self.defect.charge

    @property
    def energy(self):
        """
        Returns the *corrected* energy of the entry
        """
        return self.uncorrected_energy + np.sum(list(self.corrections.values()))

    @property
    def name(self):
        """
        Returms the defect name
        """
        return self.defect.name

    def formation_energy(self, chemical_potentials = None, fermi_level=0):
        """
        Computes the formation energy for a defect taking into account a given chemical potential and fermi_level
        """
        chemical_potentials = chemical_potentials if chemical_potentials else {}
        
        chempot_correction = sum([
            chem_pot * (self.bulk_structure.composition[el] - self.defect.defect_composition[el])
            for el, chem_pot in chemical_potentials.items()
        ])

        formation_energy = self.energy + chempot_correction

        if "vbm" in self.parameters:
            formation_energy += self.charge * (self.parameters["vbm"] + fermi_level)

        return formation_energy

    def defect_concentration(self, chemical_potentials, temperature=300, fermi_level=0.0):
        """
        Get the defect concentration for a temperature and Fermi level.
        Args:
            temperature:
                the temperature in K
            fermi_level:
                the fermi level in eV (with respect to the VBM)
        Returns:
            defects concentration in cm^-3
        """
        n = self.multiplicity * 1e24 / self.defect.bulk_structure.volume
        conc = n*exp( -1.0*self.formation_energy(chemical_potentials, fermi_level=fermi_level)/(kb*temperature))

        return conc

    def __repr__(self):
        """
        Human readable string representation of this entry
        """
        output = [
            #TODO: add defect.name abilities... maybe with composition?
            # "DefectEntry {} - {}".format(self.entry_id, self.defect.name), "Energy = {:.4f}".format(self.energy),
            "DefectEntry {} - {}".format(self.entry_id, "DEFECT"), "Energy = {:.4f}".format(self.energy),
            "Correction = {:.4f}".format(np.sum(self.corrections.values())), "Parameters:"
        ]
        for k, v in self.parameters.items():
            output.append("\t{} = {}".format(k, v))
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()


class DefectCorrection(MSONable):
    """
    A Correction class modeled off the computed entry correction format
    """

    @abstractmethod
    def get_correction(self, entry):
        """
        Returns correction for a single entry.

        Args:
            entry: A DefectEntry object.

        Returns:
            A single dictionary with the format
            correction_name: energy_correction

        Raises:
            CompatibilityError if entry is not compatible.
        """
        return

    def correct_entry(self, entry):
        """
        Corrects a single entry.

        Args:
            entry: A DefectEntry object.

        Returns:
            An processed entry.

        Raises:
            CompatibilityError if entry is not compatible.
        """
        entry.correction.update(self.get_correction(entry))
        return entry
