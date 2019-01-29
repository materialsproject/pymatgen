# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging
from abc import ABCMeta, abstractmethod

from monty.json import MSONable

from pymatgen.core import PeriodicSite
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.defects.core import Vacancy, Interstitial, Substitution
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.utils import StructureMotifInterstitial, TopographyAnalyzer
from pymatgen.analysis.structure_matcher import PointDefectComparator


__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"
"""
This module defines classes to generate point defect structures
"""

logger = logging.getLogger(__name__)


class DefectGenerator(MSONable, metaclass=ABCMeta):
    """
    Abstract class for point defects
    Implements generator pattern
    """

    def __iter__(self):
        """
        Return self as this should be an iterator
        """
        return self

    @abstractmethod
    def __next__(self):
        """
        Abstract method to return defects
        """
        return


class VacancyGenerator(DefectGenerator):
    """
    Simple generator for vacancies based on periodically
    equivalent sites
    """

    def __init__(self, structure, include_bv_charge=False):
        """
        Initializes a Vacancy Generator
        Args:
            structure(Structure): pymatgen structure object
        """
        self.structure = structure
        self.include_bv_charge = include_bv_charge

        # Find equivalent site list
        sga = SpacegroupAnalyzer(self.structure)
        self.symm_structure = sga.get_symmetrized_structure()
        self.equiv_site_seq = list(self.symm_structure.equivalent_sites)

        self.struct_valences = None
        if self.include_bv_charge:
            bv = BVAnalyzer()
            self.struct_valences = bv.get_valences(self.structure)

    def __next__(self):
        """
        Returns the next vacancy in the sequence or
        raises StopIteration
        """
        if len(self.equiv_site_seq) > 0:
            vac_site = self.equiv_site_seq.pop(0)
            charge = 0.0
            if self.struct_valences:
                site_index = self.structure.get_sites_in_sphere(vac_site[0].coords, 0.1, include_index=True)[0][2]
                charge = -1 * self.struct_valences[site_index]

            return Vacancy(self.structure, vac_site[0], charge=charge)
        else:
            raise StopIteration


class SubstitutionGenerator(DefectGenerator):
    """
    Simple generator for substitution based on periodically
    equivalent sites
    """

    def __init__(self, structure, element):
        """
        Initializes a Substitution Generator
        note: an Antisite is considered a type of substitution
        Args:
            structure(Structure): pymatgen structure object
            element (str or Element or Specie): element for the substitution
        """
        self.structure = structure
        self.element = element

        # Find equivalent site list
        sga = SpacegroupAnalyzer(self.structure)
        self.symm_structure = sga.get_symmetrized_structure()

        self.equiv_sub = []
        for equiv_site_set in list(self.symm_structure.equivalent_sites):
            vac_site = equiv_site_set[0]
            if isinstance(element, str):  # make sure you compare with specie symbol or Element type
                vac_specie = vac_site.specie.symbol
            else:
                vac_specie = vac_site.specie
            if element != vac_specie:
                defect_site = PeriodicSite(element, vac_site.coords, structure.lattice, coords_are_cartesian=True)
                sub = Substitution(structure, defect_site)
                self.equiv_sub.append(sub)

    def __next__(self):
        """
        Returns the next Substitution in the sequence or
        raises StopIteration
        """
        if len(self.equiv_sub) > 0:
            return self.equiv_sub.pop(0)
        else:
            raise StopIteration


class InterstitialGenerator(DefectGenerator):
    """
    Generator for interstitials at positions
    where the interstitialcy is coordinated by nearest neighbors
    in a way that resembles basic structure motifs
    (e.g., tetrahedra, octahedra).  The algorithm is called InFiT
    (Interstitialcy Finding Tool), it was introducted by
    Nils E. R. Zimmermann, Matthew K. Horton, Anubhav Jain,
    and Maciej Haranczyk (Front. Mater., 4, 34, 2017),
    and it is used by the Python Charged Defect Toolkit
    (PyCDT: D. Broberg et al., Comput. Phys. Commun., in press, 2018).
    """

    def __init__(self, structure, element):
        """
        Initializes an Interstitial generator using structure motifs
        Args:
            structure (Structure): pymatgen structure object
            element (str or Element or Specie): element for the interstitial
        """
        self.structure = structure
        self.element = element
        interstitial_finder = StructureMotifInterstitial(self.structure, self.element)

        self.unique_defect_seq = []
        # eliminate sublattice equivalent defects which may
        # have slipped through interstitial finder
        pdc = PointDefectComparator()

        for poss_site in interstitial_finder.enumerate_defectsites():
            now_defect = Interstitial( self.structure, poss_site)
            append_defect = True
            for unique_defect in self.unique_defect_seq:
                if pdc.are_equal( now_defect, unique_defect):
                    append_defect = False
            if append_defect:
                self.unique_defect_seq.append( now_defect)

        self.count_def = 0  # for counting the index of the generated defect

    def __next__(self):
        """
        Returns the next interstitial or
        raises StopIteration
        """
        if len(self.unique_defect_seq) > 0:
            inter_defect = self.unique_defect_seq.pop(0)
            inter_site = inter_defect.site
            self.count_def += 1
            site_name = 'InFiT' + str(self.count_def)
            return Interstitial(self.structure, inter_site, site_name=site_name)
        else:
            raise StopIteration


class VoronoiInterstitialGenerator(DefectGenerator):
    """
    Generator for interstitials based on a simple Voronoi analysis
    """

    def __init__(self, structure, element):
        """
        Initializes an Interstitial generator using Voronoi sites
        Args:
            structure (Structure): pymatgen structure object
            element (str or Element or Specie): element for the interstitial
        """
        self.structure = structure
        self.element = element

        framework = list(self.structure.symbol_set)
        get_voronoi = TopographyAnalyzer(self.structure, framework, [], check_volume=False)
        get_voronoi.cluster_nodes()
        get_voronoi.remove_collisions()

        # trim equivalent nodes with symmetry analysis
        struct_to_trim = self.structure.copy()
        for poss_inter in get_voronoi.vnodes:
            struct_to_trim.append(self.element, poss_inter.frac_coords, coords_are_cartesian=False)

        symmetry_finder = SpacegroupAnalyzer(struct_to_trim, symprec=1e-1)
        equiv_sites_list = symmetry_finder.get_symmetrized_structure().equivalent_sites

        # do additional screening for sublattice equivalent
        # defects which may have slipped through
        pdc = PointDefectComparator()
        self.unique_defect_seq = []
        for poss_site_list in equiv_sites_list:
            poss_site = poss_site_list[0]
            if poss_site not in self.structure:
                now_defect = Interstitial( self.structure, poss_site)
                append_defect = True
                for unique_defect in self.unique_defect_seq:
                    if pdc.are_equal( now_defect, unique_defect):
                        append_defect = False
                if append_defect:
                    self.unique_defect_seq.append( now_defect)

        self.count_def = 0  # for counting the index of the generated defect

    def __next__(self):
        """
        Returns the next interstitial or
        raises StopIteration
        """
        if len(self.unique_defect_seq) > 0:
            inter_defect = self.unique_defect_seq.pop(0)
            inter_site = inter_defect.site
            self.count_def += 1
            site_name = 'Voronoi' + str(self.count_def)
            return Interstitial( self.structure, inter_site, site_name=site_name)
        else:
            raise StopIteration


class MidpointInterstitialGenerator(DefectGenerator):
    """
    Generator for interstitials based on a simple Midpoint/sublattice analysis

    Algorithm is as follows:
        i) consider all pairs of atoms and find mid point between the pairs
        ii) trim mid point list such that there are no standard atomic overlap
            errors caused when appending the site
                < OPTIONAL > do this based on size of atom being considered
        iii) use sublattice generation approach to see if number of sublattice sites
            generated is (strictly) less than half the number of symmetry operations
            for the host's space group.
    """

    def __init__(self, structure, element, include_ionic_radius=True):
        """
        Initializes an Interstitial generator
        Args:
            structure (Structure): pymatgen structure object
            element (str or Element or Specie): element for the interstitial
            include_atomic_radius (bool) Whether to include the maximum ionic
                radius as an additional screening step for the element of interest
        """
        import itertools
        import numpy as np
        from pymatgen.core import Element
        from pymatgen.analysis.defects.core import create_saturated_interstitial_structure
        from pymatgen.analysis.structure_matcher import PointDefectComparator

        self.structure = structure
        self.element = element

        # create site_list with sites near border repeated in each periodic direction
        # doing this sloppy because it doesnt matter if we have too many sites in site_list
        site_list = [site.coords[:] for site in self.structure.sites]
        for site in self.structure.sites:
            if (site.frac_coords < 0.1).any():
                for ijk in itertools.product([-1,0,1],repeat=3):
                    if ijk != (0,0,0):
                        trans_vec = np.dot( ijk, self.structure.lattice.matrix)
                        site_list.append( site.coords[:] + trans_vec[:])

        print('Using {} sites for midpoint consideration \n'.format(len(site_list)))

        # generate initial mid point list and make sure
        # there are no failures with standard appending
        candidate_mid_point_site_list = []
        count_pairs = 0
        if include_ionic_radius:
            ir_val = Element(self.element).ionic_radii.values()
            dist_tol = 2. * max(ir_val) if ir_val else 0.5
        else:
            dist_tol = 0.5

        for s1, s2 in itertools.combinations(site_list, 2):
            count_pairs += 1
            half_way_coords= (s1 + s2) / 2.
            ps = PeriodicSite( self.element, half_way_coords,
                               self.structure.lattice,
                               to_unit_cell=True,
                               coords_are_cartesian=True)
            trial_struc = self.structure.copy()
            trial_struc.DISTANCE_TOLERANCE = dist_tol
            try:
                #will raise value error if bad overlap
                trial_struc.append( ps.specie, ps.coords,
                            coords_are_cartesian=True, validate_proximity=True)

                #now check if it already exists in our candidate list
                already_exists = False
                for already_found in candidate_mid_point_site_list:
                    if already_found == ps:
                        already_exists = True
                if not already_exists:
                    candidate_mid_point_site_list.append( ps)
            except:
                continue

        print('Found {} total pairs of elements\n'
              'With {} halfway points eligible for appending\n'
              ''.format(count_pairs, len(candidate_mid_point_site_list)))

        # use sublattice generation approach to see if number of sublattice sites
        # generated is (strictly) less than half the number of symmetry operations
        # for the host's space group.
        pdc = PointDefectComparator()
        sga = SpacegroupAnalyzer( self.structure)
        cut_off_size = len(self.structure) +  (len(sga.get_symmetry_operations())/2.)
        mid_point_defects = []
        for inter_site in candidate_mid_point_site_list:
            interstitial_def = Interstitial( self.structure, inter_site)
            sat_def_struct = create_saturated_interstitial_structure( interstitial_def)
            if len( sat_def_struct) < cut_off_size:
                already_exists = False
                for already_found_defect in mid_point_defects:
                    if pdc.are_equal(interstitial_def, already_found_defect):
                        already_exists = True
                if not already_exists:
                    # print('\tCoords = {} have multiplicity of {} in structure.'
                    #       ''.format(inter_site, len(sat_def_struct) - len(self.structure)))
                    mid_point_defects.append( interstitial_def)

        print("In end, produced {} midpoint-type interstitial defects.".format(len(mid_point_defects)))

        self.count_def = 0  # for counting the index of the generated defect
        self.unique_defect_seq = mid_point_defects

    def __next__(self):
        """
        Returns the next interstitial or
        raises StopIteration
        """
        if len(self.unique_defect_seq) > 0:
            inter_defect = self.unique_defect_seq.pop(0)
            inter_site = inter_defect.site
            self.count_def += 1
            site_name = 'MidPoint' + str(self.count_def)
            return Interstitial( self.structure, inter_site, site_name=site_name)
        else:
            raise StopIteration



class MasterInterstitialGenerator(DefectGenerator):
    """
    An interstitial defect generator which combines all existing
    interstitial defect generation methods in Materials Project.
    Input is structure and element of interest.

    Symmetry analysis is done to only produce symmetrically
        distinct interstitial sites.
    """

    def __init__(self, structure, element):
        """
        Initializes the Master Interstitial generator
        Args:
            structure (Structure): pymatgen structure object
            element (str or Element or Specie): element for the interstitial
        """
        import itertools
        from pymatgen.core import Element
        from pymatgen.analysis.structure_matcher import PointDefectComparator
        pdc = PointDefectComparator()

        self.structure = structure
        self.element = element

        print('Perform Voronoi Site finding method.')
        vor_inter = VoronoiInterstitialGenerator(self.structure, self.element)
        list_vor_inter = list(vor_inter)

        self.unique_defect_seq = [inter.copy() for inter in list_vor_inter]
        vor_unique = len(self.unique_defect_seq)
        print('\tproduced {} new interstitial defects'.format(vor_unique))


        print('Perform InFit Site finding method (time intensive).')
        infit_inter = InterstitialGenerator(self.structure, self.element)
        list_infit_inter = list(infit_inter)

        #reduce defects that were found previously
        for inter in list_infit_inter:
            already_exists = False
            for already_found_inter in self.unique_defect_seq:
                if pdc.are_equal(inter, already_found_inter):
                    already_exists = True
            if not already_exists:
                self.unique_defect_seq.append( inter.copy())

        infit_unique = len(self.unique_defect_seq) - vor_unique
        print('\tproduced {} new interstitial defects'.format(infit_unique))


        print('Perform MidPoint Site finding method.')
        midpoint_inter = MidpointInterstitialGenerator(self.structure, self.element)
        list_midpoint_inter = list(midpoint_inter)

        #reduce defects that were found previously
        for inter in list_midpoint_inter:
            already_exists = False
            for already_found_inter in self.unique_defect_seq:
                if pdc.are_equal(inter, already_found_inter):
                    already_exists = True
            if not already_exists:
                self.unique_defect_seq.append( inter.copy())

        midpoint_unique = len(self.unique_defect_seq) - vor_unique - infit_unique
        print('\tproduced {} new interstitial defects'.format(midpoint_unique))


        #now trim based on interstitial candidates which are too close to each other
        ir_val = Element(self.element).ionic_radii.values()
        dist_tol = 2. * max(ir_val) if ir_val else 0.5
        print("Removing defects that are within {} of each other."
              "".format(dist_tol))
        group_too_close = []
        for d1, d2 in itertools.combinations( self.unique_defect_seq, 2):
            if d1.site.distance_from_point(d2.site.coords) < dist_tol:
                # too close, but want to group based on sets of defects
                # which are too close and only keep the highest symmetry one
                new_defects = True

                for defect_set in group_too_close:
                    if (defect_set[0].site.distance_from_point(d1.site.coords) < dist_tol) or \
                        (defect_set[0].site.distance_from_point(d2.site.coords) < dist_tol):
                        d1_append = True
                        d2_append = True
                        for def_in_set in defect_set:
                            if pdc.are_equal(d1, def_in_set):
                                d1_append = False
                            if pdc.are_equal(d2, def_in_set):
                                d2_append = False
                        if d1_append:
                            defect_set.append( d1.copy())
                        if d2_append:
                            defect_set.append( d2.copy())
                        new_defects = False

                if new_defects:
                    group_too_close.append( [d1.copy(), d2.copy()])

        too_close = []
        for close_group in group_too_close:
            close_group = sorted( close_group, key=lambda defect: defect.multiplicity) #want to keep lowest multiplicity because this implies higher symmetry location
            for ind in range(1, len(close_group)):
                too_close.append( close_group[ind])

        print("This removed {} defects.".format(dist_tol, len(too_close)))
        self.unique_defect_seq = list(set(self.unique_defect_seq) - set(too_close))
        print('Overall found {} distinct defects.'.format(len(self.unique_defect_seq)))



    def __next__(self):
        """
        Returns the next interstitial or
        raises StopIteration
        """
        if len(self.unique_defect_seq) > 0:
            return self.unique_defect_seq.pop(0)
        else:
            raise StopIteration


class SimpleChargeGenerator(DefectGenerator):
    """
    Does an extremely simple/limited charge generation scheme (only one charge generated)

    for vacancies: use bond valence method to assign oxidation states and consider
                    negative of the vacant site's oxidation state as single charge to try
    for antisites and subs: use bond valence method to assign oxidation states and consider
                    negative of the vacant site's oxidation state as single charge to try +
                    added to likely charge of substitutional site (closest to zero)
    for interstitial: charge zero
    """

    def __init__(self, defect):
        """
        Args:
            defect(Defect): pymatgen Defect object
        """
        self.defect = defect

        try:
            bv = BVAnalyzer()
            struct_valences = bv.get_valences(self.defect.bulk_structure)
            site_index = self.defect.bulk_structure.get_sites_in_sphere(
                self.defect.site.coords, 0.1, include_index=True)[0][2]
            def_site_valence = struct_valences[site_index]
        except Exception:  # sometimes valences cant be assigned
            def_site_valence = 0

        if isinstance(defect, Vacancy):
            self.charges = [-1 * def_site_valence]
        elif isinstance(defect, Substitution):
            #(minimize difference with host site specie)
            probable_chgs = [ox - def_site_valence for ox in self.defect.site.specie.oxidation_states]
            self.charges = [min(probable_chgs, key=abs)]
        elif isinstance(defect, Interstitial):
            self.charges = [0]
        else:
            raise ValueError("Defect Type not recognized.")

    def __next__(self):
        """
        Returns the next defect type with the correct charge appended
        raises StopIteration
        """
        if len(self.charges) > 0:
            charge = self.charges.pop(0)
            defect = self.defect.copy()
            defect.set_charge(charge)
            return defect
        else:
            raise StopIteration
