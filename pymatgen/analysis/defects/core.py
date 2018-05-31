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
        defect_dict = self.as_dict()
        return Defect.from_dict(defect_dict)

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
        return "Vac_{}_mult{}".format(self.site.specie,self.multiplicity)

    def copy(self):
        """
        Convenience method to get a copy of the Vacancy.

        Returns:
            A copy of the Vacancy.
        """
        defect_dict = self.as_dict()
        return Vacancy.from_dict(defect_dict)


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
        return "Sub_{}_on_{}_mult{}".format(self.site.specie,self.bulk_structure[defindex].specie,self.multiplicity)

    def copy(self):
        """
        Convenience method to get a copy of the Substitution.

        Returns:
            A copy of the Substitution.
        """
        defect_dict = self.as_dict()
        return Substitution.from_dict(defect_dict)

class Interstitial(Defect):
    """
    Subclass of Defect to capture essential information for a single Interstitial defect structure.
    """

    def __init__(self, structure, defect_site, charge=0., site_name='', multiplicity=1):
        """
        Initializes an interstial defect. 
        User must specify multiplity. Default is 1
        Args:
            structure: Pymatgen Structure without any defects
            defect_site (Site): the site for the interstial
            charge: (int or float) defect charge
                default is zero, meaning no change to NELECT after defect is created in the structure
            site_name: allows user to give a unique name to defect, since Wyckoff symbol/multiplicity is
                    insufficient to categorize the defect type
                default is no name.
            multiplicity (int): multiplicity
                default is 1,
        """
        super().__init__(structure=structure,defect_site=defect_site,charge=charge)
        self._multiplicity = multiplicity
        self.site_name = site_name

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
        if self.site_name:
            return "Int_{}_{}_mult{}".format(self.site.specie,self.site_name,self.multiplicity)
        else:
            return "Int_{}_mult{}".format(self.site.specie,self.multiplicity)

    def copy(self):
        """
        Convenience method to get a copy of the Interstitial.

        Returns:
            A copy of the Interstitial.
        """
        defect_dict = self.as_dict()
        return Interstitial.from_dict(defect_dict)


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
    def site(self):
        return self.defect.site

    @property
    def multiplicty(self):
        return self.defect.multiplicty

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

    def copy(self):
        """
        Convenience method to get a copy of the DefectEntry.

        Returns:
            A copy of the DefectEntry.
        """
        defectentry_dict = self.as_dict()
        return DefectEntry.from_dict(defectentry_dict)

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
            "Correction = {:.4f}".format(np.sum(list(self.corrections.values()))), "Parameters:"
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


"""
Below is to be reviewed for accuracy (an idea for interstitial multiplicity generalization)
"""
from pymatgen.core import Molecule, PeriodicSite
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import PointGroupAnalyzer


def create_saturated_interstitial_structure(interstitial_def):
    """
    this takes a Interstitial defect object and finds multiplicity for it by decorating the defect structure with sites
    which are identical under the nearest neighbor's point group symmetry

    Process:
        1) Find nearest neighbor (nn) to defect site (just choose one) and calculate init_basis_vector from nn to defect site
        2) From NON DEFECTIVE structure create a Molecule object centered on nn,
            then generate point group operations of that molecule
            IF no centered molecule on point operations is possible -> just assume point group op is identity
        3) Create a list of basis vectors transformed which are the transformations of the init_basis_vector
            with respect to the point group operations
        4) For all equivalent sites to the nn, consider the list of basis vectors and check if an interstitial site has been
            added with this new basis vector relative to this equivalent nn site. If not - add it to the structure

    Returns:
        Structure object decorated with interstitial site equivalents
    """
    print('\nTry to find multiplicity for {}'.format(interstitial_def.name))
    #1) Find nearest neighbor to defect site (just choose one) and calculate init_basis_vector from nn to defect site
    defect_struct = interstitial_def.generate_defect_structure()
    minrad = 1000.
    nn_index = None
    init_basis_vector = None
    for pos_nn in defect_struct.get_neighbors(interstitial_def.site, 6., include_index=True):
        if pos_nn[1] < minrad:
            minrad = pos_nn[1]
            init_basis_vector = np.subtract( interstitial_def.site.coords, pos_nn[0].coords)
            nn_index = pos_nn[2]

    if nn_index == None:
        raise ValueError("Could not find nearest neighbor?")


    # 2) From NON DEFECTIVE structure create molecule centered on nn and generate point group operations of that molecule
    # 	IF no centered molecule on point operations is possible -> just assume point group op is identity
    bulk_struct = interstitial_def.bulk_structure.copy()
    center_site = defect_struct[nn_index] #this will break things if supercell is used in generate_defect_structure above...
    max_radius = 12

    range_set = np.arange(7,max_radius,.2)
    successful_centering = False
    range_ind = 0
    while not successful_centering:
        mol_coord_set = [] #zero centered site will get added
        mol_spec_set = []
        radius = range_set[range_ind]
        for test_site in bulk_struct.get_sites_in_sphere( center_site.coords, radius):
            #center the test site and append
            centered_test_site = test_site[0].coords - center_site.coords
            mol_coord_set.append(centered_test_site)
            mol_spec_set.append(test_site[0].specie)

        #check if molecule is centered
        mol = Molecule( mol_spec_set, mol_coord_set, validate_proximity=True)
        if np.linalg.norm(mol.center_of_mass) < 0.01:
            successful_centering = True

        range_ind += 1
        if range_ind == len(range_set) and not successful_centering:
            print("Unable to center molecule through {} Angstroms... just using identity operation.".format(max_radius))
            sg_ops = [SymmOp(np.eye(4))]
        elif successful_centering:
            pga = PointGroupAnalyzer(mol)
            sg_ops = pga.get_symmetry_operations()
            print('Succesfully centered molecule on nn with {} atoms; has {} symm ops'.format(len(mol), len(sg_ops)))


    # 3) Create a list of basis vectors transformed which are the transformations of the init_basis_vector
    # 	with respect to the point group operations
    basis_vector_list = []
    for sgo in sg_ops:
        new_vector = sgo.operate( init_basis_vector)
        add_vector = True
        for prev_vec in basis_vector_list:
            if np.linalg.norm( np.subtract( new_vector, prev_vec)) < 0.01:
                add_vector = False
        if add_vector:
            basis_vector_list.append( new_vector)

    print('This process ended up with {} new basis vectors'.format(len(basis_vector_list)))


    # 4) For all equivalent sites of nn, consider the list of basis vectors and check if an interstitial site has been
    # 	added with this new basis vector relative to the equivalent nn site. If not - add it to the structure
    mirrored_def_structure = defect_struct.copy()
    poss_mirr_nn_site_ind = sorted( mirrored_def_structure.get_sites_in_sphere(defect_struct[nn_index].coords, 2, include_index=True), key=lambda x: x[1])
    mirr_nn_site_ind = poss_mirr_nn_site_ind[0][2]

    sga = SpacegroupAnalyzer(mirrored_def_structure)
    mirrored_periodic_struc = sga.get_symmetrized_structure() #needed for finding equivalent images

    id_inter_list = [defect_struct[nn_index].coords]

    for equiv_site in mirrored_periodic_struc.find_equivalent_sites(mirrored_periodic_struc[mirr_nn_site_ind]):
        for basevec in basis_vector_list:
            trial_coord = equiv_site.coords + basevec
            poss_new_site = PeriodicSite(interstitial_def.site.specie, trial_coord,
                                mirrored_def_structure.lattice, to_unit_cell=True,
                                coords_are_cartesian=True)
            append_site = True
            for mdsite in mirrored_def_structure:
                if mdsite.distance( poss_new_site) < 0.5:
                    append_site = False
            if append_site:
                mirrored_def_structure.append( poss_new_site.specie, poss_new_site.coords, coords_are_cartesian=True, validate_proximity=True)
                id_inter_list.append(poss_new_site)

    print('Final multiplicity is {}'.format(len(id_inter_list)))

    return mirrored_def_structure
