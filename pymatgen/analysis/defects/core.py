# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging
import numpy as np

from abc import ABCMeta, abstractmethod
from monty.json import MSONable, MontyDecoder
from functools import lru_cache

from pymatgen.core.structure import Structure, PeriodicSite
from pymatgen.core.composition import Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.utils import kb

__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"

logger = logging.getLogger(__name__)


class Defect(MSONable, metaclass=ABCMeta):
    """
    Abstract class for a single point defect
    """

    def __init__(self, structure, defect_site, charge=0.):
        """
        Initializes an abstract defect

        Args:
            structure: Pymatgen Structure without any defects
            defect_site (Site): site for defect within structure
                must have same lattice as structure
            charge: (int or float) defect charge
                default is zero, meaning no change to NELECT after defect is created in the structure
                (assuming use_structure_charge=True in vasp input set)
        """
        self._structure = structure
        self._charge = charge
        self._defect_site = defect_site
        if structure.lattice != defect_site.lattice:
            raise ValueError("defect_site lattice must be same as structure lattice.")

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
        return self.from_dict(self.as_dict())

    def set_charge(self, new_charge=0.):
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

        #create a trivial defect structure to find where supercell transformation moves the lattice
        struct_for_defect_site = Structure( self.bulk_structure.copy().lattice,
                                             [self.site.specie],
                                             [self.site.frac_coords],
                                             to_unit_cell=True)
        struct_for_defect_site.make_supercell(supercell)
        defect_site = struct_for_defect_site[0]

        poss_deflist = sorted(
            defect_structure.get_sites_in_sphere(defect_site.coords, 2, include_index=True), key=lambda x: x[1])
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
        return "Vac_{}_mult{}".format(self.site.specie, self.multiplicity)


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

        # consider modifying velocity property to make sure defect site is decorated
        # consistently with bulk structure for final defect_structure
        defect_properties = self.site.properties.copy()
        if ('velocities' in self.bulk_structure.site_properties) and \
            'velocities' not in defect_properties:
            if all( vel == self.bulk_structure.site_properties['velocities'][0]
                    for vel in self.bulk_structure.site_properties['velocities']):
                defect_properties['velocities'] = self.bulk_structure.site_properties['velocities'][0]
            else:
                raise ValueError("No velocity property specified for defect site and "
                                 "bulk_structure velocities are not homogeneous. Please specify this "
                                 "property within the initialized defect_site object.")

        #create a trivial defect structure to find where supercell transformation moves the lattice
        site_properties_for_fake_struct = {prop: [val] for prop,val in defect_properties.items()}
        struct_for_defect_site = Structure( self.bulk_structure.copy().lattice,
                                             [self.site.specie],
                                             [self.site.frac_coords],
                                             to_unit_cell=True,
                                             site_properties = site_properties_for_fake_struct)
        struct_for_defect_site.make_supercell(supercell)
        defect_site = struct_for_defect_site[0]

        poss_deflist = sorted(
            defect_structure.get_sites_in_sphere(defect_site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]

        subsite = defect_structure.pop(defindex)
        defect_structure.append(self.site.specie.symbol, subsite.coords, coords_are_cartesian=True,
                                properties = defect_site.properties)
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
        return "Sub_{}_on_{}_mult{}".format(self.site.specie, self.bulk_structure[defindex].specie, self.multiplicity)


class Interstitial(Defect):
    """
    Subclass of Defect to capture essential information for a single Interstitial defect structure.
    """

    def __init__(self, structure, defect_site, charge=0., site_name='', multiplicity=None):
        """
        Initializes an interstial defect.
        User must specify multiplity. Default is 1
        Args:
            structure: Pymatgen Structure without any defects
            defect_site (Site): the site for the interstitial
            charge: (int or float) defect charge
                default is zero, meaning no change to NELECT after defect is created in the structure
                (assuming use_structure_charge=True in vasp input set)
            site_name: allows user to give a unique name to defect, since Wyckoff symbol/multiplicity
                is sometimes insufficient to categorize the defect type.
                 default is no name beyond multiplicity.
            multiplicity (int): multiplicity of defect within
                the supercell can be supplied by user. if not
                specified, then space group symmetry is used
                to generator interstitial sublattice

                NOTE: multiplicity generation will not work for
                interstitial complexes,
                where multiplicity may depend on additional
                factors (ex. orientation etc.)
                If defect is not a complex, then this
                process will yield the correct multiplicity,
                provided that the defect does not undergo
                significant relaxation.
        """
        super().__init__(structure=structure, defect_site=defect_site, charge=charge)
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

        # consider modifying velocity property to make sure defect site is decorated
        # consistently with bulk structure for final defect_structure
        defect_properties = self.site.properties.copy()
        if ('velocities' in self.bulk_structure.site_properties) and \
            'velocities' not in defect_properties:
            if all( vel == self.bulk_structure.site_properties['velocities'][0]
                    for vel in self.bulk_structure.site_properties['velocities']):
                defect_properties['velocities'] = self.bulk_structure.site_properties['velocities'][0]
            else:
                raise ValueError("No velocity property specified for defect site and "
                                 "bulk_structure velocities are not homogeneous. Please specify this "
                                 "property within the initialized defect_site object.")

        #create a trivial defect structure to find where supercell transformation moves the defect site
        site_properties_for_fake_struct = {prop: [val] for prop,val in defect_properties.items()}
        struct_for_defect_site = Structure( self.bulk_structure.copy().lattice,
                                             [self.site.specie],
                                             [self.site.frac_coords],
                                             to_unit_cell=True,
                                             site_properties = site_properties_for_fake_struct)
        struct_for_defect_site.make_supercell(supercell)
        defect_site = struct_for_defect_site[0]

        defect_structure.append(self.site.specie.symbol, defect_site.coords, coords_are_cartesian=True,
                                properties = defect_site.properties)
        defect_structure.set_charge(self.charge)
        return defect_structure

    @property
    def multiplicity(self):
        """
        Returns the multiplicity of a defect site within the structure (needed for concentration analysis)
        """
        if self._multiplicity is None:
            # generate multiplicity based on space group symmetry operations performed on defect coordinates
            try:
                d_structure = create_saturated_interstitial_structure(self)
            except ValueError:
                logger.debug('WARNING! Multiplicity was not able to be calculated adequately '
                             'for interstitials...setting this to 1 and skipping for now...')
                return 1

            sga = SpacegroupAnalyzer(d_structure)
            periodic_struc = sga.get_symmetrized_structure()
            poss_deflist = sorted(
                periodic_struc.get_sites_in_sphere(self.site.coords, 2, include_index=True),
                key=lambda x: x[1])
            defindex = poss_deflist[0][2]

            equivalent_sites = periodic_struc.find_equivalent_sites(periodic_struc[defindex])
            return len(equivalent_sites)

        else:
            return self._multiplicity

    @property
    def name(self):
        """
        Returns a name for this defect
        """
        if self.site_name:
            return "Int_{}_{}_mult{}".format(self.site.specie, self.site_name, self.multiplicity)
        else:
            return "Int_{}_mult{}".format(self.site.specie, self.multiplicity)


def create_saturated_interstitial_structure( interstitial_def, dist_tol=0.1):
    """
    this takes a Interstitial defect object and generates the
    sublattice for it based on the structure's space group.
    Useful for understanding multiplicity of an interstitial
    defect in thermodynamic analysis.

    NOTE: if large relaxation happens to interstitial or
        defect involves a complex then there maybe additional
        degrees of freedom that need to be considered for
        the multiplicity.

    Args:
        dist_tol: changing distance tolerance of saturated structure,
                allowing for possibly overlapping sites
                but ensuring space group is maintained

    Returns:
        Structure object decorated with interstitial site equivalents
    """
    sga = SpacegroupAnalyzer( interstitial_def.bulk_structure.copy())
    sg_ops = sga.get_symmetry_operations( cartesian=True)

    # copy bulk structure to make saturated interstitial structure out of
    # artificially lower distance_tolerance to allow for distinct interstitials
    # with lower symmetry to be replicated - This is OK because one would never
    # actually use this structure for a practical calcualtion...
    saturated_defect_struct = interstitial_def.bulk_structure.copy()
    saturated_defect_struct.DISTANCE_TOLERANCE = dist_tol

    for sgo in sg_ops:
        new_interstit_coords = sgo.operate( interstitial_def.site.coords[:])
        poss_new_site = PeriodicSite(
                interstitial_def.site.specie,
                new_interstit_coords,
                saturated_defect_struct.lattice,
                to_unit_cell=True,
                coords_are_cartesian=True)
        try:
            #will raise value error if site already exists in structure
            saturated_defect_struct.append(
                        poss_new_site.specie, poss_new_site.coords,
                        coords_are_cartesian=True, validate_proximity=True)
        except ValueError:
            pass

    # do final space group analysis to make sure symmetry not lowered by saturating defect structure
    saturated_sga = SpacegroupAnalyzer( saturated_defect_struct)
    if saturated_sga.get_space_group_number() != sga.get_space_group_number():
        raise ValueError("Warning! Interstitial sublattice generation "
                         "has changed space group symmetry. I recommend "
                         "reducing dist_tol and trying again...")

    return saturated_defect_struct


class DefectEntry(MSONable):
    """
    An lightweight DefectEntry object containing key computed data
    for many defect analysis.
    """

    def __init__(self, defect, uncorrected_energy, corrections=None, parameters=None, entry_id=None):
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

        Optional:
            note that if you intend to use this defect entry with Charge Corrections
            but the bulk_structure stored in defect is not the final supercell,
            then 'scaling_matrix' must be stored in parameters
                for example: parameters = {'scaling_matrix': [3,3,3]}
        """
        self.defect = defect
        self.uncorrected_energy = uncorrected_energy
        self.corrections = corrections if corrections else {}
        self.entry_id = entry_id
        self.parameters = parameters if parameters else {}

    @property
    def bulk_structure(self):
        return self.defect.bulk_structure

    def as_dict(self):
        """
        Json-serializable dict representation of DefectEntry
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "defect": self.defect.as_dict(),
             "uncorrected_energy": self.uncorrected_energy,
             "corrections": self.corrections,
             "parameters": self.parameters,
             "entry_id": self.entry_id}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Reconstitute a DefectEntry object from a dict representation created using
        as_dict().
         Args:
            d (dict): dict representation of DefectEntry.
         Returns:
            DefectEntry object
        """
        defect = MontyDecoder().process_decoded( d["defect"])
        uncorrected_energy = d["uncorrected_energy"]
        corrections = d.get("corrections", None)
        parameters = d.get("parameters", None)
        entry_id = d.get("entry_id", None)

        return cls(defect, uncorrected_energy, corrections=corrections,
                   parameters=parameters, entry_id=entry_id)

    @property
    def site(self):
        return self.defect.site

    @property
    def multiplicity(self):
        return self.defect.multiplicity

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

    def formation_energy(self, chemical_potentials=None, fermi_level=0):
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
        else:
            formation_energy += self.charge * fermi_level

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
        conc = n * np.exp(-1.0 * self.formation_energy(chemical_potentials, fermi_level=fermi_level) /
                          (kb * temperature))

        return conc

    def __repr__(self):
        """
        Human readable string representation of this entry
        """
        output = [
            # TODO: add defect.name abilities... maybe with composition?
            # "DefectEntry {} - {}".format(self.entry_id, self.defect.name), "Energy = {:.4f}".format(self.energy),
            "DefectEntry {} - {}".format(self.entry_id, "DEFECT"),
            "Energy = {:.4f}".format(self.energy),
            "Correction = {:.4f}".format(np.sum(list(self.corrections.values()))),
            "Parameters:"
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
