# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module defines classes for point defects.
"""

__author__ = "Bharat Medasani, Nils E. R. Zimmermann"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani, Nils E. R. Zimmermann"
__email__ = "mbkumar@gmail.com, n.zimmermann@tuhh.de"
__status__ = "Production"
__date__ = "Nov 28, 2016"


import os
import abc
import json
import numpy as np
from bisect import bisect_left
import time
from math import fabs

from pymatgen.core.periodic_table import Specie, Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.analysis.structure_analyzer import OrderParameters
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, \
    SpacegroupOperations
from pymatgen.io.zeopp import get_voronoi_nodes, get_void_volume_surfarea, \
    get_high_accuracy_voronoi_nodes
from pymatgen.command_line.gulp_caller import get_energy_buckingham, \
    get_energy_relax_structure_buckingham
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder, \
    RelaxationAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher, \
    SpeciesComparator
from pymatgen.analysis.bond_valence import BVAnalyzer
import six
from six.moves import filter
from six.moves import map
from six.moves import zip

file_dir = os.path.dirname(__file__)
rad_file = os.path.join(file_dir, 'ionic_radii.json')
with open(rad_file, 'r') as fp:
    _ion_radii = json.load(fp)


class ValenceIonicRadiusEvaluator(object):
    """
    Computes site valences and ionic radii for a structure using bond valence
    analyzer

    Args:
        structure: pymatgen.core.structure.Structure
    """

    def __init__(self, structure):
        self._structure = structure.copy()
        self._valences = self._get_valences()
        self._ionic_radii = self._get_ionic_radii()

    @property
    def radii(self):
        """
        List of ionic radii of elements in the order of sites.
        """
        el = [site.species_string for site in self._structure.sites]
        radii_dict = dict(zip(el, self._ionic_radii))
        #print radii_dict
        return radii_dict

    @property
    def valences(self):
        """
        List of oxidation states of elements in the order of sites.
        """
        el = [site.species_string for site in self._structure.sites]
        valence_dict = dict(zip(el, self._valences))
        return valence_dict

    @property
    def structure(self):
        """
        Returns oxidation state decorated structure.
        """
        return self._structure.copy()


    def _get_ionic_radii(self):
        """
        Computes ionic radii of elements for all sites in the structure.
        If valence is zero, atomic radius is used.
        """
        radii = []
        coord_finder = VoronoiCoordFinder(self._structure)

        def nearest_key(sorted_vals, key):
            i = bisect_left(sorted_vals, key)
            if i == len(sorted_vals):
                return sorted_vals[-1]
            if i == 0:
                return sorted_vals[0]
            before = sorted_vals[i-1]
            after = sorted_vals[i]
            if after-key < key-before:
                return after
            else:
                return before

        for i in range(len(self._structure.sites)):
            site = self._structure.sites[i]
            if isinstance(site.specie,Element):
                radius = site.specie.atomic_radius
                # Handle elements with no atomic_radius
                # by using calculated values instead.
                if radius is None:
                    radius = site.specie.atomic_radius_calculated
                if radius is None:
                    raise ValueError(
                            "cannot assign radius to element {}".format(
                            site.specie))
                radii.append(radius)
                continue

            el = site.specie.symbol
            oxi_state = int(round(site.specie.oxi_state))
            coord_no = int(round(coord_finder.get_coordination_number(i)))
            try:
                tab_oxi_states = sorted(map(int, _ion_radii[el].keys()))
                oxi_state = nearest_key(tab_oxi_states, oxi_state)
                radius = _ion_radii[el][str(oxi_state)][str(coord_no)]
            except KeyError:
                if coord_finder.get_coordination_number(i)-coord_no > 0:
                    new_coord_no = coord_no + 1
                else:
                    new_coord_no = coord_no - 1
                try:
                    radius = _ion_radii[el][str(oxi_state)][str(new_coord_no)]
                    coord_no = new_coord_no
                except:
                    tab_coords = sorted(map(int, _ion_radii[el][str(oxi_state)].keys()))
                    new_coord_no = nearest_key(tab_coords, coord_no)
                    i = 0
                    for val in tab_coords:
                        if  val > coord_no:
                            break
                        i = i + 1
                    if i == len(tab_coords):
                        key = str(tab_coords[-1])
                        radius = _ion_radii[el][str(oxi_state)][key]
                    elif i == 0:
                        key = str(tab_coords[0])
                        radius = _ion_radii[el][str(oxi_state)][key]
                    else:
                        key = str(tab_coords[i-1])
                        radius1 = _ion_radii[el][str(oxi_state)][key]
                        key = str(tab_coords[i])
                        radius2 = _ion_radii[el][str(oxi_state)][key]
                        radius = (radius1+radius2)/2

            #implement complex checks later
            radii.append(radius)
        return radii

    def _get_valences(self):
        """
        Computes ionic valences of elements for all sites in the structure.
        """
        try:
            bv = BVAnalyzer()
            self._structure = bv.get_oxi_state_decorated_structure(self._structure)
            valences = bv.get_valences(self._structure)
        except:
            try:
                bv = BVAnalyzer(symm_tol=0.0)
                self._structure = bv.get_oxi_state_decorated_structure(self._structure)
                valences = bv.get_valences(self._structure)
            except:
                valences = []
                for site in self._structure.sites:
                    if len(site.specie.common_oxidation_states) > 0:
                        valences.append(site.specie.common_oxidation_states[0])
                    # Handle noble gas species
                    # which have no entries in common_oxidation_states.
                    else:
                        valences.append(0)
                if sum(valences):
                    valences = [0]*self._structure.num_sites
                else:
                    self._structure.add_oxidation_state_by_site(valences)
                #raise

        #el = [site.specie.symbol for site in self._structure.sites]
        #el = [site.species_string for site in self._structure.sites]
        #el = [site.specie for site in self._structure.sites]
        #valence_dict = dict(zip(el, valences))
        #print valence_dict
        return valences


class Defect(six.with_metaclass(abc.ABCMeta, object)):
    """
    Abstract class for point defects
    """

    @abc.abstractmethod
    def enumerate_defectsites(self):
        """
        Enumerates all the symmetrically distinct defects.
        """
        raise NotImplementedError()

    @property
    def structure(self):
        """
        Returns the structure without any defects
        Useful for Mott-Littleton calculations.
        """
        return self._structure

    @property
    def struct_radii(self):
        """
        Radii of elements in the structure
        """
        return self._rad_dict

    @property
    def struct_valences(self):
        """
        Valence of elements in the structure
        """
        return self._valence_dict

    def defectsite_count(self):
        """
        Returns the number of symmetrically distinct defect sites
        """
        return len(self._defect_sites)

    def get_defectsite(self, n):
        """
        Returns the defect site at the index.
        """
        return self._defect_sites[n]

    def get_defectsite_multiplicity(self, n):
        """
        Returns the symmtric multiplicity of the defect site at the index.
        """
        return self._defect_site_multiplicity[n]

    def get_defectsite_coordination_number(self, n):
        """
        Coordination number of interstitial site.

        Args:
            n: Index of interstitial list
        """
        return self._defectsite_coord_no[n]

    def get_coordinated_sites(self, n):
        """
        The sites in structure surrounding the defect site.

        Args:
            n: Index of defects list
        """
        return self._defect_coord_sites[n]

    def get_coordinated_elements(self, n):
        """
        Elements of sites in structure surrounding the defect site.

        Args:
            n: Index of defect list
        """
        coordinated_species = []
        for site in self._defect_coord_sites[n]:
            coordinated_species.append(site.specie.symbol)
        return list(set(coordinated_species))

    @abc.abstractmethod
    def make_supercells_with_defects(self, scaling_matrix):
        """
        Generate the supercell with input multipliers and create the defect.
        First supercell has no defects.
        To create unit cell with defect pass unit matrix.
        """
        raise NotImplementedError()


class Vacancy(Defect):
    """
    Subclass of Defect to generate vacancies and their analysis.

    Args:
        structure: pymatgen.core.structure.Structure
        valences: valences of elements as a dictionary
        radii: Radii of elements as a dictionary
    """

    def __init__(self, structure, valences, radii):
        self._structure = structure
        self._valence_dict = valences
        self._rad_dict = radii
        # Store symmetrically distinct sites, their coordination numbers
        # coordinated_sites, effective charge
        symm_finder = SpacegroupAnalyzer(self._structure)
        symm_structure = symm_finder.get_symmetrized_structure()
        equiv_site_seq = symm_structure.equivalent_sites

        self._defect_sites = []
        self._defect_site_multiplicity = []
        for equiv_sites in equiv_site_seq:
            self._defect_sites.append(equiv_sites[0])
            self._defect_site_multiplicity.append(len(equiv_sites))

        self._vac_site_indices = []
        for site in self._defect_sites:
            for i in range(len(self._structure.sites)):
                if site == self._structure[i]:
                    self._vac_site_indices.append(i)

        coord_finder = VoronoiCoordFinder(self._structure)
        self._defectsite_coord_no = []
        self._defect_coord_sites = []
        for i in self._vac_site_indices:
            self._defectsite_coord_no.append(
                coord_finder.get_coordination_number(i)
            )
            self._defect_coord_sites.append(
                coord_finder.get_coordinated_sites(i)
            )

        # Store the ionic radii for the elements in the structure
        # (Used to  computing the surface are and volume)
        # Computed based on valence of each element

        self._vac_eff_charges = None
        self._vol = None
        self._sa = None

    #@property
    #def valence_dict(self):
    #    return self._valence_dict

    def enumerate_defectsites(self):
        """
        Returns symmetrically distinct vacancy sites
        """
        return self._defect_sites

    def get_defectsite_structure_indices(self):
        """
        Returns indices of symmetrically distinct vacancy sites
        """
        return self._vac_site_indices

    def get_defectsite_structure_index(self, n):
        """
        index of the vacacy site in the structure.sites list

        Args:
            n:
                Index of vacancy list
        """
        return self._vac_site_indices[n]

    def get_defectsite_effective_charge(self, n):
        """
        Effective charge (In Kroger-Vink notation, cation vacancy has
        effectively -ve charge and anion vacancy has +ve charge.)

        Args:
            n: Index of vacancy list

        Returns:
            Effective charnge of defect site
        """
        # Effective charge (In Kroger-Vink notation, cation vacancy has
        # effectively -ve charge and anion vacancy has +ve charge.) Inverse
        # the BVAnalyzer.get_valences result.

        el = self.get_defectsite(n).species_string
        return -self._valence_dict[el]
        #if not self._vac_eff_charges:
        #    self._vac_eff_charges = []
        #    for site in self.enumerate_defectsites():
        #        specie = site.specie.symbol
        #        self._vac_eff_charges.append(-self._valence_dict[specie])
        #return self._vac_eff_charges[n]

    def get_coordsites_min_max_charge(self, n):
        """
        Minimum and maximum charge of sites surrounding the vacancy site.

        Args:
            n: Index of vacancy list
        """
        bv = BVAnalyzer()
        struct_valences = bv.get_valences(self._structure)
        coordinated_site_valences = []

        def _get_index(site):
            for i in range(len(self._structure.sites)):
                if site.is_periodic_image(self._structure.sites[i]):
                    return i
            raise ValueError("Site not found")

        for site in self._defect_coord_sites[n]:
            ind = _get_index(site)
            coordinated_site_valences.append(struct_valences[ind])
        coordinated_site_valences.sort()
        return coordinated_site_valences[0], coordinated_site_valences[-1]

    # deprecated
    def get_volume(self, n):
        """
        Volume of the nth vacancy

        Args:
            n: Index of symmetrically distinct vacancies list

        Returns:
            floating number representing volume of vacancy
        """
        if not self._vol:
            self._vol = []
            self._sa = []
            um = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            sc = self.make_supercells_with_defects(um)[1:]
            rad_dict = self.struct_radii
            for i in range(len(sc)):
                vol, sa = get_void_volume_surfarea(sc[i], rad_dict)
                self._vol.append(vol)
                self._sa.append(sa)

        return self._vol[n]

    # deprecated
    def get_surface_area(self, n):
        """
        Surface area of the nth vacancy

        Args:
            n: Index of symmetrically distinct vacancies list

        Returns:
            floating number representing volume of vacancy
        """
        if not self._sa:
            self._vol = []
            self._sa = []
            um = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            supercells = self.make_supercells_with_defects(um)[1:]
            rad_dict = self.struct_radii
            for sc in supercells:
                vol, sa = get_void_volume_surfarea(sc, rad_dict)
                self._vol.append(vol)
                self._sa.append(sa)

        return self._sa[n]

    def _supercell_with_defect(self, scaling_matrix, defect_site):
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        oldf_coords = defect_site.frac_coords
        coords = defect_site.lattice.get_cartesian_coords(oldf_coords)
        newf_coords = sc.lattice.get_fractional_coords(coords)
        sc_defect_site = PeriodicSite(defect_site.species_and_occu, newf_coords,
                                      sc.lattice,
                                      properties=defect_site.properties)
        for i in range(len(sc.sites)):
            #if sc_defect_site == sc.sites[i]:
            if sc_defect_site.distance(sc.sites[i]) < 1e-3:
                del sc[i]
                return sc
        raise ValueError('Something wrong if reached here')

    def make_supercells_with_defects(self, scaling_matrix, species=None,
                                     limit_return_structures=False):
        """
        Generate sequence of supercells in pymatgen.core.structure.Structure
        format, with each supercell containing one vacancy.

        Args:
            scaling_matrix: super cell scale parameters in matrix forms
            species: Species in list format only for which vacancy supercells
                are required. If not specified all the species are considered.
            limit_return_structures: Boolean or positive number
                If number, only that many structures are returned.

        Returns:
            Supercells with vacancies. First supercell has no defects.
        """
        sc_with_vac = []
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        sc_with_vac.append(sc)

        if not species:
            species = sc.symbol_set
        if not limit_return_structures:
            limit_return_structures = self.defectsite_count()
        for defect_site in self.enumerate_defectsites():
            if len(sc_with_vac) <= limit_return_structures:
                if isinstance(defect_site.specie,Specie):
                    site_specie = defect_site.specie.element.symbol
                elif isinstance(defect_site.specie,Element):
                    site_specie = defect_site.specie.symbol
                else:
                    raise TypeError("site specie is neither Specie nor Element")

                if site_specie in species:
                    sc_with_vac.append(self._supercell_with_defect(
                        scaling_matrix, defect_site))
        return sc_with_vac


class VacancyFormationEnergy(object):
    """
    Using GULP compute the vacancy formation energy.
    Works only for binary metal oxides due to the use of Buckingham Potentials
    """

    def __init__(self, vacancy):
        self._vacancy = vacancy
        self._energies = []

    def get_energy(self, n, tol=0.5):
        """
        Formation Energy for nth symmetrically distinct vacancy.
        """
        #generate defect free structure energy
        if not self._energies:
            no_vac = self._vacancy.defectsite_count()
            prev_energies = [0.0] * no_vac
            tol_flg = [False] * no_vac
            vac_gulp_kw = ('optimise', 'conp', 'qok')
            val_dict = self._vacancy.struct_valences
            for sp in range(2, 6):
                if not (False in tol_flg):
                    #print sp
                    break
                scale_mat = [[sp, 0, 0], [0, sp, 0], [0, 0, sp]]
                sc = self._vacancy.make_supercells_with_defects(scale_mat)
                blk_energy = get_energy_buckingham(sc[0])
                no = len(sc[0].sites)
                #print no
                for i in range(1, no_vac + 1):
                    if not tol_flg[i - 1]:
                        vac_energy = get_energy_buckingham(
                            sc[i], keywords=vac_gulp_kw,
                            valence_dict=val_dict
                        )
                        form_energy = vac_energy - (no - 1) / no * blk_energy
                        if abs(form_energy - prev_energies[i - 1]) < tol:
                            tol_flg[i - 1] = True
                        prev_energies[i - 1] = form_energy

            self._energies = prev_energies
            self._tol_flg = tol_flg

        if not self._tol_flg[n]:
            print("Caution: tolerance not reached for {0} vacancy".format(n))
        return self._energies[n]


class Interstitial(Defect):
    """
    Subclass of Defect to generate interstitial sites
    """

    def __init__(self, structure, valences, radii, site_type='voronoi_vertex',
                 accuracy='Normal', symmetry_flag=True, oxi_state = False):
        """
        Given a structure, generate symmetrically distinct interstitial sites.
        For a non-ionic structure, use oxi_state=True and give atomic radii.

        Args:
            structure: pymatgen.core.structure.Structure
            valences: Dictionary of oxidation states of elements in 
                {el:valence} form
            radii: Radii of elemnts in the structure
            site_type: "voronoi_vertex" uses voronoi nodes
                "voronoi_edgecenter" uses voronoi polyhedra edge centers
                "voronoi_facecenter" uses voronoi polyhedra face centers
                "all" combines vertices, edgecenters and facecenters.
                Default is "voronoi_vertex"
            accuracy: Flag denoting whether to use high accuracy version 
                of Zeo++. Options are "Normal" and "High". Default is normal.
            symmetry_flag: If True, only returns symmetrically distinct sites
            oxi_state: If False, input structure is considered devoid of 
                oxidation-state decoration. And oxi-state for each site is 
                determined. Use True, if input structure is oxi-state 
                decorated. This option is useful when the structure is 
                not electro-neutral after deleting/adding sites. In that
                case oxi-decorate the structure before deleting/adding the
                sites.
        """
        if not oxi_state:
            self._structure = ValenceIonicRadiusEvaluator(structure).structure
        else:
            self._structure = structure

        self._valence_dict = valences
        self._rad_dict = radii

        """
        Use Zeo++ to obtain the voronoi nodes. Apply symmetry reduction
        and the symmetry reduced voronoi nodes are possible candidates
        for interstitial sites.
        """
        if accuracy == "Normal":
            high_accuracy_flag = False
        elif accuracy == "High":
            high_accuracy_flag = True
        else:
            raise NotImplementedError("Accuracy setting not implemented.")

        if accuracy == "High":
            if site_type in ('voronoi_facecenter','voronoi_edgecenter','all'): 
                raise NotImplementedError(
                        "Site type not implemented for the accuracy setting")


        vor_node_sites, vor_edgecenter_sites, vor_facecenter_sites = \
            symmetry_reduced_voronoi_nodes(self._structure, self._rad_dict,
                                           high_accuracy_flag, symmetry_flag)
        
        if site_type == 'voronoi_vertex':
            possible_interstitial_sites = vor_node_sites
        elif site_type == 'voronoi_facecenter':
            possible_interstitial_sites = vor_facecenter_sites
        elif site_type == 'voronoi_edgecenter':
            possible_interstitial_sites = vor_edgecenter_sites
        elif site_type == "all":
            possible_interstitial_sites = vor_node_sites + \
                    vor_facecenter_sites + vor_edgecenter_sites
        else:
            raise ValueError("Input site type not implemented")

        #Do futher processing on possibleInterstitialSites to obtain
        #interstitial sites
        self._defect_sites = possible_interstitial_sites
        self._defectsite_coord_no = []
        self._defect_coord_sites = []
        self._defect_coord_charge = []
        self._radii = []

        for site in self._defect_sites:
            coord_no, coord_sites, chrg = self._get_coord_no_sites_chrg(site)
            self._defectsite_coord_no.append(coord_no)
            self._defect_coord_sites.append(coord_sites)
            self._defect_coord_charge.append(chrg)

        for site in self._defect_sites:
            vor_radius = site.properties.get('voronoi_radius',None)
            if vor_radius:
                vor_radius = float(vor_radius)
            self._radii.append(vor_radius)

    def _get_coord_no_sites_chrg(self, site):
        """
        Compute the coordination number and coordination charge

        Args:
            site:
                pymatgen.core.sites.Site
        """
        struct = self._structure.copy()
        struct.append(site.specie.symbol, site.frac_coords)
        coord_finder = VoronoiCoordFinder(struct)
        coord_no = coord_finder.get_coordination_number(-1)
        coord_sites = coord_finder.get_coordinated_sites(-1)

        # In some cases coordination sites to interstitials include
        # interstitials also. Filtering them.
        def no_inter(site):
            return not site.specie.symbol == 'X'
        coord_sites = filter(no_inter, coord_sites)

        coord_chrg = 0
        if self._valence_dict:
            for site, weight in coord_finder.get_voronoi_polyhedra(-1).items():
                if not site.specie.symbol == 'X':
                    coord_chrg += weight * self._valence_dict[site.species_string]

        return coord_no, coord_sites, coord_chrg

    def enumerate_defectsites(self):
        """
        Enumerate all the symmetrically distinct interstitial sites.
        The defect site has "X" as occupied specie.
        """
        return self._defect_sites

    def append_defectsite(self, site):
        """
        Append a site to list of possible interstitials

        Args:
            site: pymatgen.core.sites.Site
        """
        raise NotImplementedError()

    def delete_defectsite(self, n):
        """
        Remove a symmetrically distinct interstitial site

        Args:
            n: Index of interstitial site
        """
        del self._defect_sites[n]

    def get_coordsites_charge_sum(self, n):
        """
        Total charge of the interstitial coordinated sites.

        Args:
            n: Index of interstitial list
        """
        return self._defect_coord_charge[n]

    def get_coordsites_min_max_charge(self, n):
        """
        Minimum and maximum charge of sites surrounding the interstitial site.

        Args:
            n: Index of symmetrical distinct interstitial site
        """
        coord_site_valences = []

        for site in self._defect_coord_sites[n]:
            coord_site_valences.append(self._valence_dict[site.specie.symbol])
        coord_site_valences.sort()
        return coord_site_valences[0], coord_site_valences[-1]

    def get_radius(self, n):
        """
        Volume of the nth interstitial

        Args:
            n: Index of symmetrically distinct vacancies list

        Returns:
            floating number representing radius of interstitial sphere
        """
        return self._radii[n]

    def get_radii(self):
        return self._radii

    def reduce_defectsites(self):
        """
        If multiple defect sites have same voronoi radius, only one is kept.
        Useful if the symmetry based reduction of initial sites returned
        from Zeo++ is not working properly due to deviation in ideal lattice
        coordinates.
        """
        distinct_radii = list(set(self._radii))
        for rad in distinct_radii:
            ind = self._radii.index(rad)  # Index of first site with 'rad'
            for i in reversed(list(range(ind + 1, len(self._radii)))):
                # Backward search for remaining sites so index is not changed
                if self._radii[i] == rad:
                    self._defect_sites.pop(i)
                    self._defectsite_coord_no.pop(i)
                    self._defect_coord_sites.pop(i)
                    self._radii.pop(i)

    def radius_prune_defectsites(self, radius):
        """
        Remove all the defect sites with voronoi radius less than input radius
        """
        for i in reversed(list(range(len(self._radii)))):
            if self._radii[i] < radius:
                self._defect_sites.pop(i)
                self._defectsite_coord_no.pop(i)
                self._defect_coord_sites.pop(i)
                self._radii.pop(i)

    def prune_defectsites(self, el="C", oxi_state=4, dlta=0.1):
        """
        Prune all the defect sites which can't acoomodate the input elment
        with the input oxidation state.
        """
        rad = float(Specie(el, oxi_state).ionic_radius) - dlta
        self.radius_prune_defectsites(rad)

    def prune_close_defectsites(self, dist=0.2):
        """
        Prune the sites that are very close.
        """
        #print self.defectsite_count()
        ind = 0
        while ind < self.defectsite_count():
            #i = ind + 1
            #while i < self.defectsite_count():
            i = self.defectsite_count()-1
            #print ind, i
            while i > ind: 
                d = self._defect_sites[ind].distance(self._defect_sites[i])
                #print d, dist
                if d < dist:
                    self._defect_sites.pop(i)
                    #self._defectsite_coord_no.pop(i)
                    #self._defect_coord_sites.pop(i)
                    #self._radii.pop(i)
            #    i += 1
                i -= 1
            ind += 1
        #print self.defectsite_count()

    def _supercell_with_defect(self, scaling_matrix, defect_site, element):
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        oldf_coords = defect_site.frac_coords
        coords = defect_site.lattice.get_cartesian_coords(oldf_coords)
        #print coords
        newf_coords = sc.lattice.get_fractional_coords(coords)
        for i in range(3):
            coord = newf_coords[i]
            if coord < 0:
                while (coord < 0):
                    coord = coord+1
                newf_coords[i] = coord
            elif coord > 1:
                while (coord > 1):
                    coord = coord-1
                newf_coords[i] = coord
        #print newf_coords
        #sc_defect_site = PeriodicSite(element, newf_coords,
        #                              sc.lattice)
        try:
            sc.append(element, newf_coords, coords_are_cartesian=False,
                      validate_proximity=True)
        except ValueError:
            sc = None
        finally:
            return sc

    def make_supercells_with_defects(self, scaling_matrix, element):
        """
        Returns sequence of supercells in pymatgen.core.structure.Structure
        format, with each supercell containing an interstitial.
        First supercell has no defects.
        """
        sc_list_with_interstitial = []
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        sc_list_with_interstitial.append(sc)
        for defect_site in self.enumerate_defectsites():
            sc_with_inter = self._supercell_with_defect(
                scaling_matrix, defect_site, element
            )
            if sc_with_inter:
                sc_list_with_interstitial.append(sc_with_inter)
        return sc_list_with_interstitial


class InterstitialAnalyzer(object):
    """
    Use GULP to compute the interstitial formation energy, relaxed structures.
    Works only for metal oxides due to the use of Buckingham Potentials.

    Args:
        inter: pymatgen.defects.point_defects.Interstitial
        el: Element name in short hand notation ("El")
        oxi_state: Oxidtation state
        scd: Super cell dimension as number. The scaling is equal along xyz.
    """
    def __init__(self, inter, el, oxi_state, scd=2):
        self._inter = inter
        self._el = el
        self._oxi_state = oxi_state
        self._scd = scd
        self._relax_energies = []
        self._norelax_energies = []
        self._relax_struct = []

    def get_energy(self, n, relax=True):
        """
        Formation Energy for nth symmetrically distinct interstitial.
        """
        if relax and not self._relax_energies:
            self._relax_analysis()
        if not relax and not self._norelax_energies:
            no_inter = self._inter.defectsite_count()
            inter_gulp_kw = ('qok',)
            val_dict = self._inter.struct_valences
            val_dict[self._el] = self._oxi_state  # If element not in structure
            scd = self._scd
            scale_mat = [[scd, 0, 0], [0, scd, 0], [0, 0, scd]]
            sc = self._inter.make_supercells_with_defects(scale_mat, self._el)
            blk_energy = get_energy_buckingham(sc[0])
            for i in range(1, no_inter + 1):
                inter_energy = get_energy_buckingham(
                    sc[i], keywords=inter_gulp_kw, valence_dict=val_dict
                )
                form_energy = inter_energy - blk_energy
                self._norelax_energies.append(form_energy)

        if relax:
            return self._relax_energies[n]
        else:
            return self._norelax_energies[n]

    def _relax_analysis(self):
        """
        Optimize interstitial structures
        """

        no_inter = self._inter.defectsite_count()
        inter_gulp_kw = ('optimise', 'conp', 'qok')
        val_dict = self._inter.struct_valences
        scd = self._scd
        scale_mat = [[scd, 0, 0], [0, scd, 0], [0, 0, scd]]
        sc = self._inter.make_supercells_with_defects(scale_mat, self._el)
        blk_energy, rlx_struct = get_energy_relax_structure_buckingham(sc[0])
        self._relax_struct.append(rlx_struct)
        val_dict[self._el] = self._oxi_state  # If element not in structure
        for i in range(1, no_inter + 1):
            energy, rlx_struct = get_energy_relax_structure_buckingham(
                sc[i], keywords=inter_gulp_kw, valence_dict=val_dict
            )
            form_energy = energy - blk_energy
            self._relax_energies.append(form_energy)
            self._relax_struct.append(rlx_struct)

    def get_relaxed_structure(self, n):
        """
        Optimized interstitial structure

        Args:
            n: Symmetrically distinct interstitial index

        .. note::

            To get relaxed bulk structure pass -1.
            -ve index will not work as expected.
        """
        if not self._relax_struct:
            self._relax_analysis()
        return self._relax_struct[n + 1]

    def get_percentage_volume_change(self, n):
        """
        Volume change after the introduction of interstitial

        Args:
            n: Symmetrically distinct interstitial index
        """
        if not self._relax_struct:
            self._relax_analysis()
        blk_struct = self._relax_struct[0]
        def_struct = self._relax_struct[n + 1:n + 2][0]
        del def_struct.sites[-1]
        rv = RelaxationAnalyzer(blk_struct, def_struct)
        return rv.get_percentage_volume_change()

    def get_percentage_lattice_parameter_change(self, n):
        """
        Lattice parameter change after the introduction of interstitial

        Args:
            n: Symmetrically distinct interstitial index
        """
        if not self._relax_struct:
            self._relax_analysis()
        blk_struct = self._relax_struct[0]
        def_struct = self._relax_struct[n + 1:n + 2][0]
        del def_struct.sites[-1]
        rv = RelaxationAnalyzer(blk_struct, def_struct)
        return rv.get_percentage_lattice_parameter_changes()

    def get_percentage_bond_distance_change(self, n):
        """
        Bond distance change after the introduction of interstitial

        Args:
            n: Symmetrically distinct interstitial index
        """
        if not self._relax_struct:
            self._relax_analysis()
        blk_struct = self._relax_struct[0]
        def_struct = self._relax_struct[n + 1:n + 2][0]
        del def_struct.sites[-1]
        #print def_struct
        rv = RelaxationAnalyzer(blk_struct, def_struct)
        return rv.get_percentage_bond_dist_changes()

    def relaxed_structure_match(self, i, j):
        """
        Check if the relaxed structures of two interstitials match

        Args:
            i: Symmetrically distinct interstitial index
            j: Symmetrically distinct interstitial index

        .. note::

            To use relaxed bulk structure pass -1.
            -ve index will not work as expected
        """
        if not self._relax_struct:
            self._relax_analysis()
        sm = StructureMatcher()
        struct1 = self._relax_struct[i + 1]
        struct2 = self._relax_struct[j + 1]
        return sm.fit(struct1, struct2)


class StructureRelaxer(object):
    def __init__(self, structure):
        self._unrelax_struct = structure
        self.relax()

    def relax(self):
        energy, rlx_struct = get_energy_relax_structure_buckingham(
            self._unrelax_struct)
        self._relax_struct = rlx_struct

    def get_relaxed_structure(self):
        return self._relax_struct


class InterstitialStructureRelaxer(object):
    """
    Performs structural relaxation for each interstitial supercell.

    Args:
        interstitial: Unrelaxed interstitial
        el: Species string in short notation
        oxi_state: Oxidation state of the element
        supercell_dim: Dimensions of super cell
    """

    def __init__(self, interstitial, el, oxi_state, supercell_dim=2):
        self._inter = interstitial
        self._scd = supercell_dim
        self._el = el
        self._oxi_state = oxi_state
        self._relax_structs = []
        self._relax_energies = []

    def relax(self):
        """
        Optimize interstitial structures
        """

        no_inter = self._inter.defectsite_count()
        inter_gulp_kw = ('optimise', 'conp', 'qok')
        val_dict = self._inter.struct_valences
        scd = self._scd
        scale_mat = [[scd, 0, 0], [0, scd, 0], [0, 0, scd]]
        sc = self._inter.make_supercells_with_defects(scale_mat, self._el)
        blk_energy, rlx_struct = get_energy_relax_structure_buckingham(sc[0])
        self._relax_structs.append(rlx_struct)
        self._relax_energies.append(blk_energy)
        val_dict[self._el] = self._oxi_state  # If element not in structure
        for i in range(1, no_inter + 1):
            try:
                energy, rlx_struct = get_energy_relax_structure_buckingham(
                    sc[i], keywords=inter_gulp_kw, valence_dict=val_dict
                )
                self._relax_energies.append(energy)
                self._relax_structs.append(rlx_struct)
            except:
                self._relax_energies.append(None)
                self._relax_structs.append(None)

        def is_empty(lst):
            for value in lst:
                if value:
                    return False
            return True

        if is_empty(self._relax_energies):
            raise IOError('Relaxation failed')

    def relaxed_structure_match(self, i, j):
        """
        Check if the relaxed structures of two interstitials match

        Args:
            i: Symmetrically distinct interstitial index
            j: Symmetrically distinct interstitial index

        .. note::

            Index 0 corresponds to bulk.
        """
        if not self._relax_structs:
            self.relax()
        sm = StructureMatcher()
        struct1 = self._relax_structs[i]
        struct2 = self._relax_structs[j]
        return sm.fit(struct1, struct2)

    def relaxed_energy_match(self, i, j):
        """
        Check if the relaxed energies of two interstitials match

        Args:
            i: Symmetrically distinct interstitial index
            j: Symmetrically distinct interstitial index

        .. note::

            Index 0 corresponds to bulk.
        """
        if not self._relax_energies:
            self.relax()
        energy1 = self._relax_energies[i]
        energy2 = self._relax_energies[j]
        return energy1 == energy2

    def get_relaxed_structure(self, n):
        """
        Get the relaxed structure of nth symmetrically distinct interstitial.

        Args:
            n: Symmetrically distinct interstitial index

        .. note::

            0 corresponds to relaxed bulk structure
        """
        if not self._relax_structs:
            self.relax()
        return self._relax_structs[n]

    def get_relaxed_energy(self, n):
        """
        Get the relaxed structure of nth symmetrically distinct interstitial.

        Args:
            n: Symmetrically distinct interstitial index

        .. note::

            0 corresponds to relaxed bulk structure
        """
        if not self._relax_energies:
            self.relax()
        return self._relax_energies[n]

    def get_relaxed_interstitial(self):
        """
        Get the relaxed structure of nth symmetrically distinct interstitial.

        Args:
            n: Symmetrically distinct interstitial index
        """
        if not self._relax_energies:
            self.relax()
        energies = self._relax_energies[:]
        structs = self._relax_structs[:]
        distinct_energy_set = set(energies[1:])  # only interstitial energies
        if None in distinct_energy_set:
            distinct_energy_set.remove(None)
        distinct_structs = [structs[0]]         # bulk
        distinct_energies = [energies[0]]
        for energy in distinct_energy_set:
            ind = energies.index(energy)
            distinct_structs.append(structs[ind])
            distinct_energies.append(energies[ind])
        return RelaxedInterstitial(
            distinct_structs, distinct_energies, self._inter.struct_valences
        )


class RelaxedInterstitial(object):
    """
    Stores the relaxed supercell structures for each interstitial
    Used to compute formation energies, displacement of atoms near the
    the interstitial.

    Args:
        struct_list: List of structures(supercells). The first structure should
            represent relaxed bulk structure and the subsequent ones
            interstitial structures (with the extra interstitial site
            appended at the end).
        energy_list: List of energies for relaxed interstitial structures.
            The first energy should correspond to bulk structure
        valence_dict: Valences of elements in dictionary form
    """

    def __init__(self, struct_list, energy_list, valence_dict):
        self._blk_struct = struct_list[0]
        struct_list.pop(0)
        self._structs = struct_list
        self._blk_energy = energy_list[0]
        energy_list.pop(0)
        self._energies = energy_list
        self._valence_dict = valence_dict

        self._coord_no = []
        self._coord_sites = []
        self._coord_charge_no = []

    def formation_energy(self, n, chem_pot=0):
        """
        Compute the interstitial formation energy

        Args:
            n: Index of interstitials
            chem_pot: Chemical potential of interstitial site element.
                If not given, assumed as zero. The user is strongly
                urged to supply the chemical potential value
        """
        return self._energies[n] - self._blk_energy - chem_pot

    def get_percentage_volume_change(self, n):
        """
        Volume change after the introduction of interstitial

        Args:
            n: index of interstitials
        """
        def_struct = self._structs[n:n + 1][0]
        del def_struct.sites[-1]
        rv = RelaxationAnalyzer(self._blk_struct, def_struct)
        return rv.get_percentage_volume_change()

    def get_percentage_lattice_parameter_change(self, n):
        """
        Lattice parameter change after the introduction of interstitial

        Args:
            n: index of interstitials
        """
        def_struct = self._structs[n:n + 1][0]  # copy
        del def_struct.sites[-1]
        rv = RelaxationAnalyzer(self._blk_struct, def_struct)
        return rv.get_percentage_lattice_parameter_changes()

    def get_percentage_bond_distance_change(self, n):
        """
        Bond distance change after the introduction of interstitial.

        Args:
            n: index of interstitials
        """
        def_struct = self._structs[n:n + 1][0]  # copy
        del def_struct.sites[-1]
        rv = RelaxationAnalyzer(self._blk_struct, def_struct)
        return rv.get_percentage_bond_dist_changes()

    def get_bulk_structure(self):
        """
        Return relaxed bulk structure
        """
        return self._blk_struct

    def get_interstitial_structure(self, n):
        """
        Return relaxed bulk structure
        """
        return self._structs[n]

    def defect_count(self):
        """
        Returns the number of distinct interstitials
        """
        return len(self._structs)

    def get_defectsite(self, n):
        """
        Returns the defect site of nth interstitial.

        Args:
            n: Index of interstitial
        """
        return self._structs[n][-1]

    def get_coordination_number(self, n):
        """
        Coordination number for nth interstitial.

        Args:
            n: Index of interstitials
        """
        if not self._coord_no:
            self._coord_find()
        return self._coord_no[n]

    def get_charge_coordination_number(self, n):
        """
        Charge coordination number for nth interstitial.

        Args:
            n: Index of interstitials
        """
        if not self._coord_charge_no:
            self._coord_find()
        return self._coord_charge_no[n]

    def get_coordinated_sites(self, n):
        """
        Coordinated sites for nth interstitial.

        Args:
            n: Index of interstitials
        """
        if not self._coord_sites:
            self._coord_find()
        return self._coord_sites[n]

    def get_coordinated_bulk_sites(self, n):
        """
        Bulk sites corresponding to the coordinated sites for nth interstitial.

        Args:
            n: Index of interstitials
        """
        blk_sites = []
        for site in self.get_coordinated_sites(n):
            site_index = self._structs[n].sites.index(site)
            blk_sites.append(self._blk_struct[site_index])
        return blk_sites

    def get_coordinated_site_displacement(self, n):
        """
        Compute the total displacement of coordinated sites from the
        interstitial sites during the relaxation

        Args:
            n: Index  of defect site
        """
        coord_sites = self.get_coordinated_sites(n)
        coord_blk_sites = self.get_coordinated_bulk_sites(n)
        dist_sum = 0
        for i in range(len(coord_sites)):
            dist_sum += coord_sites[i].distance_from_point(coord_blk_sites[i])
            # How to compute the average?
        return dist_sum


    def _coord_find(self):
        """
        calls VoronoiCoordFinder to compute the coordination number,
        coordination charge
        """
        for i in range(self.defect_count()):
            struct = self._structs[i].copy()
            coord_finder = VoronoiCoordFinder(struct)
            self._coord_no.append(coord_finder.get_coordination_number(-1))
            self._coord_sites.append(coord_finder.get_coordinated_sites(-1))
            coord_chrg = 0
            for site, weight in coord_finder.get_voronoi_polyhedra(-1).items():
                coord_chrg += weight * self._valence_dict[site.species_string]
            self._coord_charge_no.append(coord_chrg)


def symmetry_reduced_voronoi_nodes(
        structure, rad_dict, high_accuracy_flag=False, symm_flag=True):
    """
    Obtain symmetry reduced voronoi nodes using Zeo++ and
    pymatgen.symmetry.finder.SpacegroupAnalyzer

    Args:
        strucutre: pymatgen Structure object
        rad_dict: Dictionary containing radii of spcies in the structure
        high_accuracy_flag: Flag denotting whether to use high accuracy version of Zeo++
        symm_flag: Flag denoting whether to return symmetrically distinct sites only

    Returns:
        Symmetrically distinct voronoi nodes as pymatgen Strucutre
    """

    def add_closest_equiv_site(dist_sites, equiv_sites):
        if not dist_sites:
            dist_sites.append(equiv_sites[0])
        else:
            avg_dists = []
            for site in equiv_sites:
                dists = [site.distance(dst_site, jimage=[0, 0, 0])
                         for dst_site in dist_sites]
                avg_dist = sum(dists) / len(dist_sites)
                avg_dists.append(avg_dist)

            min_avg_dist = min(avg_dists)
            ind = avg_dists.index(min_avg_dist)
            dist_sites.append(equiv_sites[ind])

    def cmp_memoize_last_site(f): #Compares and stores last site
        def not_duplicates(site1, site2):
            if site1.distance(site2) < 1e-5:
                return False
            else:
                return True

        cmp_memoize_last_site.cache = None
        def helper(x):
            if not cmp_memoize_last_site.cache: 
                cmp_memoize_last_site.cache = f(x)
                return True
            y = f(x)
            if not_duplicates(cmp_memoize_last_site.cache, y):
                cmp_memoize_last_site.cache = y
                return True
            else:
                return False
        return helper

    @cmp_memoize_last_site
    def check_not_duplicates(site):
        return site


    if not symm_flag:
        if not high_accuracy_flag:
            vor_node_struct, vor_edgecenter_struct, vor_facecenter_struct = \
                get_voronoi_nodes(structure, rad_dict)
            return vor_node_struct.sites, vor_edgecenter_struct.sites, \
                   vor_facecenter_struct.sites
        else:
            # Only the nodes are from high accuracy voronoi decomposition
            vor_node_struct = \
                    get_high_accuracy_voronoi_nodes(structure, rad_dict)
            # Before getting the symmetry, remove the duplicates
            vor_node_struct.sites.sort(key = lambda site: site.voronoi_radius)
            #print type(vor_node_struct.sites[0])
            dist_sites = filter(check_not_duplicates, vor_node_struct.sites)
            return dist_sites, None, None

    if not high_accuracy_flag:
        def get_dist_sites(vor_struct):
            SpgA = SpacegroupAnalyzer
            try:
                symmetry_finder = SpgA(vor_struct, symprec=1e-1)
                symm_struct = symmetry_finder.get_symmetrized_structure()
            except:
                vor_struct.merge_sites(0.1, 'delete')
                symmetry_finder = SpgA(vor_struct, symprec=1e-1)
                symm_struct = symmetry_finder.get_symmetrized_structure()
            equiv_sites_list = symm_struct.equivalent_sites

            if not equiv_sites_list:
                dist_sites = vor_struct.sites
            else:
                dist_sites = []
                for equiv_sites in equiv_sites_list:
                    add_closest_equiv_site(dist_sites, equiv_sites)
            return dist_sites

        vor_node_struct, vor_edgecenter_struct, vor_facecenter_struct = \
            get_voronoi_nodes(structure, rad_dict)

        node_dist_sites = get_dist_sites(vor_node_struct)
        edgecenter_dist_sites = get_dist_sites(vor_edgecenter_struct)
        facecenter_dist_sites = get_dist_sites(vor_facecenter_struct)

        return node_dist_sites, edgecenter_dist_sites, facecenter_dist_sites
    else:
        # Only the nodes are from high accuracy voronoi decomposition
        vor_node_struct = \
                get_high_accuracy_voronoi_nodes(structure, rad_dict)

        # Before getting the symmetry, remove the duplicates
        vor_node_struct.sites.sort(key = lambda site: site.voronoi_radius)
        #print type(vor_node_struct.sites[0])
        dist_sites = list(filter(check_not_duplicates, vor_node_struct.sites))

        # Ignore symmetry from ha voronoi nodes
        # Increase the symmetry precision to 0.25
        #spg = SpacegroupAnalyzer(structure,symprec=1e-1).get_spacegroup()
        
        # Remove symmetrically equivalent sites
        #i = 0
        #while (i < len(dist_sites)-1):
        #    sites1 = [dist_sites[i]]
        #    sites2 = [dist_sites[i+1]]
        #    if spg.are_symmetrically_equivalent(sites1,sites2):
        #        del dist_sites[i+1]
        #    else:
        #        i = i+1

        node_dist_sites = dist_sites
        return (node_dist_sites, None, None)

        #vor_edge_symmetry_finder = SpacegroupAnalyzer(
        #    vor_edgecenter_struct, symprec=1e-1)
        #vor_edge_symm_struct = vor_edge_symmetry_finder.get_symmetrized_structure()
        #edgecenter_equiv_sites_list = vor_edge_symm_struct.equivalent_sites

        #edgecenter_dist_sites = []
        #for equiv_sites in edgecenter_equiv_sites_list:
        #    add_closest_equiv_site(edgecenter_dist_sites, equiv_sites)
        #if not edgecenter_equiv_sites_list:     
        #    edgecenter_dist_sites = vor_edgecenter_struct.sites

        #vor_fc_symmetry_finder = SpacegroupAnalyzer(
        #                vor_facecenter_struct, symprec=1e-1)
        #vor_fc_symm_struct = vor_fc_symmetry_finder.get_symmetrized_structure()
        #facecenter_equiv_sites_list = vor_fc_symm_struct.equivalent_sites

        #facecenter_dist_sites = []
        #for equiv_sites in facecenter_equiv_sites_list:
        #    add_closest_equiv_site(facecenter_dist_sites, equiv_sites)
        #if not facecenter_equiv_sites_list:     
        #    facecenter_dist_sites = vor_facecenter_struct.sites

        #return node_dist_sites, edgecenter_dist_sites, facecenter_dist_sites


def get_neighbors_of_site_with_index(struct, n, p=None):
    """
    Determine the neighbors around the site that has index n in the input
    Structure object struct, given the approach defined by parameters
    p.  All supported neighbor-finding approaches and listed and
    explained in the following.  All approaches start by creating a
    tentative list of neighbors using a large cutoff radius defined in
    parameter dictionary p via key "cutoff".
    "min_dist": find nearest neighbor and its distance d_nn; consider all
            neighbors which are within a distance of d_nn * (1 + delta),
            where delta is an additional parameter provided in the
            dictionary p via key "delta".
    "scaled_VIRE": compute the radii, r_i, of all sites on the basis of
            the valence-ionic radius evaluator (VIRE); consider all
            neighbors for which the distance to the central site is less
            than the sum of the radii multiplied by an a priori chosen
            parameter, delta,
            (i.e., dist < delta * (r_central + r_neighbor)).
    "min_relative_VIRE": same approach as "min_dist", except that we
            use relative distances (i.e., distances divided by the sum of the
            atom radii from VIRE).
    "min_relative_OKeeffe": same approach as "min_relative_VIRE", except
            that we use the bond valence parameters from O'Keeffe's bond valence
            method (J. Am. Chem. Soc. 1991, 3226-3229) to calculate
            relative distances.

    Args:
        struct (Structure): input structure.
        n (int): index of site in Structure object for which
                neighbors are to be determined.
        p (dict): specification (via "approach" key; default is "min_dist")
                and parameters of neighbor-finding approach.
                Default cutoff radius is 6 Angstrom (key: "cutoff").
                Other default parameters are as follows.
                min_dist: "delta": 0.15;
                min_relative_OKeeffe: "delta": 0.05;
                min_relative_VIRE: "delta": 0.05;
                scaled_VIRE: "delta": 2.

    Returns: ([site]) list of sites that are considered to be nearest
            neighbors to site with index n in Structure object struct.
    """
    sites = []
    if p is None:
        p = {"approach": "min_dist", "delta": 0.15,
                "cutoff": 6}

    if p["approach"] not in [
            "min_relative_OKeeffe", "min_dist", "min_relative_VIRE", \
            "scaled_VIRE"]:
        raise RuntimeError("Unsupported neighbor-finding approach"
                " (\"{}\")".format(p["approach"]))

    if p["approach"] == "min_relative_OKeeffe" or p["approach"] == "min_dist":
        neighs_dists = struct.get_neighbors(struct[n], p["cutoff"])
        try:
            eln = struct[n].specie.element
        except:
            eln = struct[n].species_string
    elif p["approach"] == "scaled_VIRE" or p["approach"] == "min_relative_VIRE":
        vire = ValenceIonicRadiusEvaluator(struct)
        if np.linalg.norm(struct[n].coords-vire.structure[n].coords) > 1e-6:
            raise RuntimeError("Mismatch between input structure and VIRE structure.")
        neighs_dists = vire.structure.get_neighbors(vire.structure[n], p["cutoff"])
        rn = vire.radii[vire.structure[n].species_string]

    reldists_neighs = []
    for neigh, dist in neighs_dists:
        if p["approach"] == "scaled_VIRE":
            dscale = p["delta"] * (vire.radii[neigh.species_string] + rn)
            if dist < dscale:
                sites.append(neigh)
        elif p["approach"] == "min_relative_VIRE":
            reldists_neighs.append([dist / (
                    vire.radii[neigh.species_string] + rn), neigh])
        elif p["approach"] == "min_relative_OKeeffe":
            try:
                el2 = neigh.specie.element
            except:
                el2 = neigh.species_string
            reldists_neighs.append([dist / get_okeeffe_distance_prediction(
                    eln, el2), neigh])
        elif p["approach"] == "min_dist":
            reldists_neighs.append([dist, neigh])

    if p["approach"] == "min_relative_VIRE" or \
            p["approach"] == "min_relative_OKeeffe" or \
            p["approach"] == "min_dist":
        min_reldist = min([reldist for reldist, neigh in reldists_neighs])
        for reldist, neigh in reldists_neighs:
            if reldist / min_reldist < 1.0 + p["delta"]:
                sites.append(neigh)

    return sites



class StructureMotifInterstitial(Defect):

    """
    Subclass of Defect to generate interstitial sites at positions
    where the interstitialcy is coordinated by nearest neighbors
    in a way that resembles basic structure motifs
    (e.g., tetrahedra, octahedra).  The algorithm will be formally
    introducted in an upcoming publication
    by Nils E. R. Zimmermann, Anubhav Jain, and Maciej Haranczyk,
    and it is already used by the Python Charged Defect Toolkit
    (PyCDT, https://arxiv.org/abs/1611.07481).
    """

    __supported_types = ("tetalt", "octalt", "bcc")

    def __init__(self, struct, inter_elem,
                 motif_types=("tetalt", "octalt"),
                 op_threshs=(0.3, 0.5),
                 dl=0.2, doverlap=1.0, facmaxdl=1.01, verbose=False):
        """
        Generates symmetrically distinct interstitial sites at positions
        where the interstitial is coordinated by nearest neighbors
        in a pattern that resembles a supported structure motif
        (e.g., tetrahedra, octahedra).

        Args:
            struct (Structure): input structure for which symmetrically
                distinct interstitial sites are to be found.
            inter_elem (string): element symbol of desired interstitial.
            motif_types ([string]): list of structure motif types that are
                to be considered.  Permissible types are:
                tet (tetrahedron), oct (octahedron).
            op_threshs ([float]): threshold values for the underlying order
                parameters to still recognize a given structural motif
                (i.e., for an OP value >= threshold the coordination pattern
                match is positive, for OP < threshold the match is
                negative.
            dl (float): grid fineness in Angstrom.  The input
                structure is divided into a grid of dimension
                a/dl x b/dl x c/dl along the three crystallographic
                directions, with a, b, and c being the lengths of
                the three lattice vectors of the input unit cell.
            doverlap (float): distance that is considered
                to flag an overlap between any trial interstitial site
                and a host atom.
            facmaxdl (float): factor to be multiplied with the maximum grid
                width that is then used as a cutoff distance for the
                clustering prune step.
            verbose (bool): flag indicating whether (True) or not (False;
                default) to print additional information to screen.
        """
        # Initialize interstitial finding.
        self._structure = struct.copy()
        self._motif_types = motif_types[:]
        if len(self._motif_types) == 0:
            raise RuntimeError("no motif types provided.")
        self._op_threshs = op_threshs[:]
        for imotif, motif in enumerate(self._motif_types):
            if motif not in self.__supported_types:
                raise RuntimeError("unsupported motif type: {}.".format(
                        motif))
        self._dl = dl
        self._defect_sites = []
        self._defect_types = []
        self._defect_site_multiplicity = []
        self._defect_cns = []
        self._defect_opvals = []

        rots, trans = SpacegroupAnalyzer(
                struct)._get_symmetry()        
        nbins = [int(struct.lattice.a / dl), \
                int(struct.lattice.b / dl), \
                int(struct.lattice.c / dl)]
        dls = [struct.lattice.a / float(nbins[0]), \
                struct.lattice.b / float(nbins[1]), \
                struct.lattice.c / float(nbins[2])]
        maxdl = max(dls)
        if verbose:
            print("Grid size: {} {} {}".format(nbins[0], nbins[1], nbins[2]))
            print("dls: {} {} {}".format(dls[0], dls[1], dls[2]))
        struct_w_inter = struct.copy()
        struct_w_inter.append(inter_elem, [0, 0, 0])
        natoms = len(list(struct_w_inter.sites))
        ops = OrderParameters(motif_types, cutoff=-10.0)
        trialsites = []

        # Loop over trial positions that are based on a regular
        # grid in fractional coordinate space
        # within the unit cell.
        for ia in range(nbins[0]):
            a = (float(ia)+0.5) / float(nbins[0])
            for ib in range(nbins[1]):
                b = (float(ib)+0.5) / float(nbins[1])
                for ic in range(nbins[2]):
                    c = (float(ic)+0.5) / float(nbins[2])
                    struct_w_inter.replace(
                            natoms-1, inter_elem, coords=[a, b, c],
                            coords_are_cartesian=False)
                    if len(struct_w_inter.get_sites_in_sphere(struct_w_inter.sites[natoms-1].coords, doverlap)) == 1:
                        delta = 0.1
                        ddelta = 0.1
                        delta_end = 0.8
                        while delta < delta_end:
                            neighs = get_neighbors_of_site_with_index(
                                    struct_w_inter, natoms-1, p={
                                    "approach": "min_dist", "delta": delta,
                                    "cutoff": 6})
                            nneighs = len(neighs)
                            if nneighs > 6:
                                break
                            if nneighs not in [4, 6]:
                                delta += ddelta
                                continue

                            allsites = [s for s in neighs]
                            indeces_neighs = [i for i in range(len(allsites))]
                            allsites.append(struct_w_inter.sites[natoms-1])
                            opvals = ops.get_order_parameters(
                                        allsites, len(allsites)-1,
                                        indeces_neighs=indeces_neighs)
                            motif_type = "unrecognized"
                            if "tetalt" in motif_types:
                                if nneighs == 4 and \
                                        opvals[motif_types.index("tetalt")] > \
                                        op_threshs[motif_types.index("tetalt")]:
                                    motif_type = "tet"
                                    this_op = opvals[motif_types.index("tetalt")]
                            if "octalt" in motif_types:
                                if nneighs == 6 and \
                                        opvals[motif_types.index("octalt")] > \
                                        op_threshs[motif_types.index("octalt")]:
                                    motif_type = "oct"
                                    this_op = opvals[motif_types.index("octalt")]

                            if motif_type != "unrecognized":
                                cns = {}
                                for isite, site in enumerate(neighs):
                                    if isinstance(site.specie, Element):
                                        elem = site.specie.symbol
                                    else:
                                        elem = site.specie.element.symbol
                                    if elem in list(cns.keys()):
                                        cns[elem] = cns[elem] + 1
                                    else:
                                        cns[elem] = 1
                                trialsites.append({
                                        "mtype": motif_type,
                                        "opval": this_op,
                                        "delta": delta,
                                        "coords": struct_w_inter.sites[natoms-1].coords[:],
                                        "fracs": np.array([a, b, c]),
                                        "cns": dict(cns)})
                                break

                            delta += ddelta

        # Prune list of trial sites by clustering and find the site
        # with the largest order parameter value in each cluster.
        nintersites = len(trialsites)
        unique_motifs = []
        for ts in trialsites:
            if ts["mtype"] not in unique_motifs:
                unique_motifs.append(ts["mtype"])
        labels = {}
        connected = []
        for i in range(nintersites):
            connected.append([])
            for j in range(nintersites):
                dist, image = struct_w_inter.lattice.get_distance_and_image(
                        trialsites[i]["fracs"],
                        trialsites[j]["fracs"])
                connected[i].append(True if dist < (maxdl*facmaxdl) else False)
        include = []
        for motif in unique_motifs:
            labels[motif] = []
            for i, ts in enumerate(trialsites):
                labels[motif].append(i if ts["mtype"] == motif else -1)
            change = True
            while change:
                change = False
                for i in range(nintersites-1):
                    if change:
                        break
                    if labels[motif][i] == -1:
                        continue
                    for j in range(i+1, nintersites):
                        if labels[motif][j] == -1:
                            continue
                        if connected[i][j] and labels[motif][i] != labels[motif][j]:
                            if labels[motif][i] < labels[motif][j]:
                                labels[motif][j] = labels[motif][i]
                            else:
                                labels[motif][i] = labels[motif][j]
                            change = True
                            break
            unique_ids = []
            for l in labels[motif]:
                if l != -1 and l not in unique_ids:
                    unique_ids.append(l)
            if verbose:
                print("unique_ids {} {}".format(motif, unique_ids))
            for uid in unique_ids:
                maxq = 0.0
                imaxq = -1
                for i in range(nintersites):
                    if labels[motif][i] == uid:
                        if imaxq < 0 or trialsites[i]["opval"] > maxq:
                            imaxq = i
                            maxq = trialsites[i]["opval"]
                include.append(imaxq)

        # Prune by symmetry.
        multiplicity = {}
        discard = []
        for motif in unique_motifs:
            discard_motif = []
            for indi, i in enumerate(include):
                if trialsites[i]["mtype"] != motif or \
                        i in discard_motif:
                    continue
                multiplicity[i] = 1
                symposlist = [trialsites[i]["fracs"].dot(
                        np.array(m, dtype=float)) for m in rots]
                for t in trans:
                    symposlist.append(trialsites[i]["fracs"]+np.array(t))
                for indj in range(indi+1, len(include)):
                    j = include[indj]
                    if trialsites[j]["mtype"] != motif or \
                            j in discard_motif:
                        continue
                    for sympos in symposlist:
                        dist, image = struct.lattice.get_distance_and_image(
                                sympos, trialsites[j]["fracs"])
                        if dist < maxdl * facmaxdl:
                            discard_motif.append(j)
                            multiplicity[i] += 1
                            break
            for i in discard_motif:
                if i not in discard:
                    discard.append(i)

        if verbose:
            print("Initial trial sites: {}\nAfter clustering: {}\n"
                    "After symmetry pruning: {}".format(
                    len(trialsites), len(include),
                    len(include)-len(discard)))
        c = 0
        for i in include:
            if i not in discard:
                self._defect_sites.append(
                    PeriodicSite(
                        Element(inter_elem),
                        trialsites[i]["fracs"],
                        self._structure.lattice,
                        to_unit_cell=False,
                        coords_are_cartesian=False,
                        properties=None))
                self._defect_types.append(trialsites[i]["mtype"])
                self._defect_cns.append(trialsites[i]["cns"])
                self._defect_site_multiplicity.append(multiplicity[i])
                self._defect_opvals.append(trialsites[i]["opval"])


    def enumerate_defectsites(self):
        """
        Get all defect sites.

        Returns:
            defect_sites ([PeriodicSite]): list of periodic sites
                    representing the interstitials.
        """
        return self._defect_sites


    def get_motif_type(self, i):
        """
        Get the motif type of defect with index i (e.g., "tet").

        Returns:
            motif (string): motif type.
        """
        return self._defect_types[i]


    def get_coordinating_elements_cns(self, i):
        """
        Get element-specific coordination numbers of defect with index i.

        Returns:
            elem_cn (dict): dictionary storing the coordination numbers (int)
                    with string representation of elements as keys.
                    (i.e., {elem1 (string): cn1 (int), ...}).
        """
        return self._defect_cns[i]


    def get_op_value(self, i):
        """
        Get order-parameter value of defect with index i.

        Returns:
            opval (float): OP value.
        """
        return self._defect_opvals[i]


    def make_supercells_with_defects(self, scaling_matrix):
        """
        Generate a sequence of supercells
        in which each supercell contains a single interstitial,
        except for the first supercell in the sequence
        which is a copy of the defect-free input structure.

        Args:
            scaling_matrix (3x3 integer array): scaling matrix
                to transform the lattice vectors.
        Returns:
            scs ([Structure]): sequence of supercells.

        """
        scs = []
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        scs.append(sc)
        for ids, defect_site in enumerate(self._defect_sites):
            sc_with_inter = sc.copy()
            sc_with_inter.append(defect_site.species_string,
                    defect_site.frac_coords,
                    coords_are_cartesian=False,
                    validate_proximity=False,
                    properties=None)
            if not sc_with_inter:
                raise RuntimeError("could not generate supercell with"
                        " interstitial {}".format(ids+1))
            scs.append(sc_with_inter.copy())
        return scs


