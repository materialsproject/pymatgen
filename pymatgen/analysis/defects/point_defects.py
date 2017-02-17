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
        vor_node_struct, vor_edgecenter_struct, vor_facecenter_struct = \
            get_voronoi_nodes(structure, rad_dict)
        vor_node_symmetry_finder = SpacegroupAnalyzer(vor_node_struct, symprec=1e-1)
        vor_node_symm_struct = vor_node_symmetry_finder.get_symmetrized_structure()
        node_equiv_sites_list = vor_node_symm_struct.equivalent_sites

        node_dist_sites = []
        for equiv_sites in node_equiv_sites_list:
            add_closest_equiv_site(node_dist_sites, equiv_sites)

        vor_edge_symmetry_finder = SpacegroupAnalyzer(
            vor_edgecenter_struct, symprec=1e-1)
        vor_edge_symm_struct = vor_edge_symmetry_finder.get_symmetrized_structure()
        edgecenter_equiv_sites_list = vor_edge_symm_struct.equivalent_sites

        edgecenter_dist_sites = []
        for equiv_sites in edgecenter_equiv_sites_list:
            add_closest_equiv_site(edgecenter_dist_sites, equiv_sites)
        if not edgecenter_equiv_sites_list:     # Fix this so doesn't arise
            edgecenter_dist_sites = vor_edgecenter_struct.sites

        vor_fc_symmetry_finder = SpacegroupAnalyzer(
                        vor_facecenter_struct, symprec=1e-1)
        vor_fc_symm_struct = vor_fc_symmetry_finder.get_symmetrized_structure()
        facecenter_equiv_sites_list = vor_fc_symm_struct.equivalent_sites

        facecenter_dist_sites = []
        for equiv_sites in facecenter_equiv_sites_list:
            add_closest_equiv_site(facecenter_dist_sites, equiv_sites)
        if not facecenter_equiv_sites_list:     # Fix this so doesn't arise
            facecenter_dist_sites = vor_facecenter_struct.sites

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


class StructureMotifInterstitial(Defect):

    """
    Subclass of Defect to generate interstitial sites at positions
    where the interstitialcy is coordinated by nearest neighbors
    in a way that resembles basic structure motifs
    (tetrahedra, octahedra, bcc) or an overlay of two motifs
    (tetrahedral-octahedral).
    The algorithm will be formally introducted in an upcoming publication
    by Nils E. R. Zimmermann and Maciej Haranczyk,
    and it is already used by the Python Charged Defect Toolkit
    (PyCDT, https://arxiv.org/abs/1611.07481).
    """

    __supported_types = ("tet", "oct", "bcc", "tetoct")

    def __init__(self, structure, inter_elem,
                 motif_types=["tet", "oct", "tetoct"],
                 op_targets=[1.0, 1.0, 1.0], op_threshs=[0.5, 0.5, 0.5],
                 dl=0.2, fac_max_radius=4.5, drel_overlap=0.5,
                 write_timings=False):
        """
        Generate symmetrically distinct interstitial sites at positions
        where the interstitial is coordinated by nearest neighbors
        in a pattern that resembles a supported structure motif
        (tetrahedra, octahedra, bcc) or an overlay of some motifs
        (tetrahedral-octahedral).

        Args:
            struct (Structure): input structure for which symmetrically
                distinct interstitial sites are to be found.
            inter_elem (string): element symbol of desired interstitial.
            motif_types ([string]): list of structure motif types that are
                to be considered.  Permissible types are:
                tet (tetrahedron), oct (octahedron),
                tetoct (tetrahedron-octahedron overlay).
            op_targets ([float]): target values for the underlying order
                parameters to recognize a given structural motif.
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
            fac_max_radius (float): factor to generate a (large) neighbor
                list per grid point with a fixed large cutoff distance,
                d_cutoff.  The list is then sorted to give the
                neighbors in an ascending list in relative distance from
                the interstitial site.  The relative distance between
                an interstitial trial site and a "bulk" atom is calculated
                by the distance, d, and the (ionic) radii of the two
                species involved, (r_inter+r_bulk), yielding
                d_rel = d / (r_inter+r_bulk).
                Because d_cutoff = fac_max_radius * r_max,
                where r_max is the largest site radius encountered in the
                structure, fac_max_radius should always be larger than 2.
            drel_overlap (float): relative distance that is considered
                to flag an overlap between any two atoms.  It is used
                a) to skip a given grid point because of overlap
                between an interstitial trial site and any one atom
                from the input structure and
                b) to remove trial sites that are too close to each other.
                The latter may or may not be desirable.  A future revision
                should make this step optional.
            write_timings (boolean): flag indicating whether or not to
                write out times for sections of the interstitial search
                (default: False).
        """

        if write_timings:
            file_time = open("time_interfinding.dat", "w")
            tt0, tc0 = time.time(), time.clock()

        self._structure = structure
        sgops = SpacegroupAnalyzer(structure).get_space_group_operations() # get_spacegroup()
        self._motif_types = motif_types
        self._unique_op_types = []
        self._map_imotif_iop = []
        for imotif, motif in enumerate(self._motif_types):
            if motif not in self.__supported_types:
                raise RuntimeError("unsupported motif type: {}.".format(
                        motif))
            self._map_imotif_iop.append([])
            if motif == "tet" or motif == "oct" or motif == "bcc":
                if motif not in self._unique_op_types:
                    self._unique_op_types.append(motif)
                self._map_imotif_iop[imotif].append(
                    self._unique_op_types.index(motif))
            elif motif == "tetoct":
                if "tet" not in self._unique_op_types:
                    self._unique_op_types.append("tet")
                self._map_imotif_iop[imotif].append(
                    self._unique_op_types.index("tet"))
                if "oct" not in self._unique_op_types:
                    self._unique_op_types.append("oct")
                self._map_imotif_iop[imotif].append(
                    self._unique_op_types.index("oct"))
            else:
                raise ValueError("cannot associate motif {} to any"
                        " order parameter type.".format(motif))
        if len(motif_types) != len(op_targets) or \
                len(motif_types) != len(op_threshs):
            raise ValueError("list sizes of structure motif types,"
                    " OP target values and OP threshold values are"
                    " not equal.")
        self._op_targets = op_targets
        self._op_threshs = op_threshs
        self._dl = dl
        self._defect_sites = []
        self._defect_types = []
        self._defect_cns = []
        self._defect_site_multiplicity = []

        # Set up list of elements found in the input structure.
        elem_list = []
        natoms_uc = 0
        for site in structure:
            natoms_uc = natoms_uc + 1
            if isinstance(site.specie, Element):
                elem = site.specie.symbol
            elif isinstance(site.specie, Specie):
                elem = site.specie.element.symbol
            else:
                raise RuntimeError("unexpected instance type")
            if elem not in elem_list:
                elem_list.append(elem)
        if inter_elem not in elem_list:
            elem_list_plus_inter = elem_list + [inter_elem]
        else:
            elem_list_plus_inter = list(elem_list)
        nelems = len(elem_list)
        if nelems < 1:
            raise RuntimeError("nelems < 1")

        # Create a working supercell.
        sc_struct = self._structure.copy()
        sc_struct.make_supercell([[1,0,0], [0,1,0], [0,0,1]])
        natoms = natoms_uc
        lat = sc_struct.lattice

        if write_timings:
            tt1, tc1 = time.time(), time.clock()
            file_time.write("init {} {}\n".format(tt1-tt0, tc1-tc0))
            tt0, tc0 = time.time(), time.clock()

        # Determine sensible radii of each site in the input structure.
        vire = ValenceIonicRadiusEvaluator(sc_struct)
        radii_sites = []
        species_string_sites = []
        for isite, site in enumerate(vire.structure):
            if np.linalg.norm(
                    site.coords-sc_struct.sites[isite].coords) > 1.0e-6:
                raise RuntimeError("inconsistent structures (coords)")
            if len(sc_struct.sites[isite].species_string) == 1:
                if site.species_string[0] != \
                        sc_struct.sites[isite].species_string[0]:
                    raise RuntimeError("inconsistent structures"
                            " (species string (len 1))")
            elif len(sc_struct.sites[isite].species_string) == 2:
                if site.species_string[0] != \
                        sc_struct.sites[isite].species_string[0] or \
                        site.species_string[1] != \
                        sc_struct.sites[isite].species_string[1]:
                    raise RuntimeError("inconsistent structures"
                            " (species string (len 2))")
            else:
                raise RuntimeError("inconsistent structures"
                        " (len species string)")
            radii_sites.append(vire.radii[site.species_string])
            species_string_sites.append(site.species_string)

        if write_timings:
            tt1, tc1 = time.time(), time.clock()
            file_time.write("radiideter {} {}\n".format(tt1-tt0, tc1-tc0))
            tt0, tc0 = time.time(), time.clock()

        # Insert interstitial site into structure.
        sc_struct.append(inter_elem, [0, 0, 0])
        natoms = natoms + 1
        radii_sites.append(Element(inter_elem).average_ionic_radius)
        species_string_sites.append(Element(inter_elem))
        if len(sc_struct.sites) != natoms:
            raise RuntimeError("")
        if sc_struct.sites[natoms-1].species_string != inter_elem:
            raise RuntimeError("{} != {}".format(
                    sc_struct.sites[natoms-1].species_string, inter_elem))
        max_radius_sphere = fac_max_radius * max(radii_sites)

        # Set up index list of all possible tet- and oct-subset neighbor lists
        # for tet-oct overlay motif, given CN = 10.
        # Note that the interstitial is always placed at index 0 in the
        # site lists so that the neighbor indeces start with 1---not with 0.
        # Hence, the indeces i, j, k, and m are shifted one position "upwards"
        # before included in lists_tet_indeces and lists_oct_indeces.
        lists_tet_indeces = []
        lists_oct_indeces = []
        i_lists_tet_indeces = 0
        if "tetoct" in self._motif_types:
            for i in range(0, 6+1):
                for j in range(i+1, 7+1):
                    for k in range(j+1, 8+1):
                        for m in range(k+1, 9+1):
                            lists_tet_indeces.append([i+1, j+1, k+1, m+1])
                            lists_oct_indeces.append([])
                            for n in range(1, 10+1):
                                if n not in lists_tet_indeces[i_lists_tet_indeces]:
                                    lists_oct_indeces[i_lists_tet_indeces].append(n)
                            i_lists_tet_indeces = i_lists_tet_indeces + 1

        # Number of voxels in input structure
        # and working supercell
        # to span an equally fine grid in both structures.
        nbins = [int(self._structure.lattice.a / self._dl), \
                int(self._structure.lattice.b / self._dl), \
                int(self._structure.lattice.c / self._dl)]
        nbins_sc = [float(1 * nbins[0]), float(1 * nbins[1]), \
                float(1 * nbins[2])]

        if write_timings:
            tt1, tc1 = time.time(), time.clock()
            file_time.write("init2 {} {}\n".format(tt1-tt0, tc1-tc0))
            tt0, tc0 = time.time(), time.clock()

        # Prescreen to find the site with the largest nearest-neighbor
        # distance.  Use that location to determine a sensible
        # global r_inter.
        max_min_dist = -1.0
        coords_max_min_dist = None
        for ia in range(nbins[0]):
            a = (float(ia)+0.5) / nbins_sc[0]
            for ib in range(nbins[1]):
                b = (float(ib)+0.5) / nbins_sc[1]
                for ic in range(nbins[2]):
                    c = (float(ic)+0.5) / nbins_sc[2]
                    sc_struct.replace(
                            natoms-1, inter_elem, coords=[a, b, c],
                            coords_are_cartesian=False)
                    neighs_dists = sc_struct.get_neighbors(
                            sc_struct.sites[natoms-1],
                            max_radius_sphere, include_index=False)
                    if len(neighs_dists) > 0:
                        tmp = min([dist for (neigh, dist) in neighs_dists])
                        if tmp > max_min_dist:
                            max_min_dist = tmp
                            coords_max_min_dist = [a, b, c]

        if write_timings:
            tt1, tc1 = time.time(), time.clock()
            file_time.write("prescreen {} {}\n".format(tt1-tt0, tc1-tc0))
            tt0, tc0 = time.time(), time.clock()

        if not coords_max_min_dist:
            raise RuntimeError("could not find any neighbor during presreening"
                    " to assign hard-sphere radius to interstitial species."
                    " Suggesting to increase fac_max_radius.")
        sc_struct.replace(natoms-1, inter_elem, coords=coords_max_min_dist,
                coords_are_cartesian=False)
        vire = ValenceIonicRadiusEvaluator(sc_struct)
        radii_sites[natoms-1] = vire.radii[vire.structure[natoms-1].species_string]
        species_string_sites[natoms-1] = vire.structure[natoms-1].species_string

        # Create underlying order parameter objects.
        ops = []
        iop=0
        for optype in self._unique_op_types:
            ops.append(OrderParameters([optype], cutoff=-10.0))
            iop=iop+1

        # Create lists to store tentative interstitial positions
        # and further information.
        inter_sites = []
        inter_sites_pruned = []
        for i_elem in range(nelems+1):
            inter_sites.append({})
            inter_sites_pruned.append({})
            for t in self._motif_types:
                inter_sites[i_elem][t] = []
                inter_sites_pruned[i_elem][t] = []

        if write_timings:
            tt1, tc1 = time.time(), time.clock()
            file_time.write("init3 {} {}\n".format(tt1-tt0, tc1-tc0))
            tt0, tc0 = time.time(), time.clock()

        # Loop over trial positions that are based on a regular
        # grid in fractional coordinate space
        # within the unit cell.
        ntrials_tet = ntrials_oct = ntrials_tetoct = 0
        for ia in range(nbins[0]):
            a = (float(ia)+0.5) / nbins_sc[0]
            for ib in range(nbins[1]):
                b = (float(ib)+0.5) / nbins_sc[1]
                for ic in range(nbins[2]):
                    c = (float(ic)+0.5) / nbins_sc[2]

                    # Place interstitial at next trial position.
                    r = lat.get_cartesian_coords([a, b, c])
                    sc_struct.replace(
                            natoms-1, inter_elem, coords=[a, b, c],
                            coords_are_cartesian=False)

                    # Determine the current list of neighbors using
                    # max_radius_sphere which is the product of the
                    # largest species radius encountered in the bulk
                    # structure and the input factor fac_max_radius.
                    rinter = radii_sites[natoms-1]
                    two_rinter = 2.0 * rinter
                    neighs_dists_indeces = sc_struct.get_neighbors(
                            sc_struct.sites[natoms-1],
                            max_radius_sphere, include_index=True)

                    # Calculate relative distances between trial site
                    # and neighbors and check whether trial site
                    # overlaps with any of the neighbors using the
                    # input relative distance drel_overlap.
                    skip_this = False
                    rel_dist = []
                    for i_rel_dist, (neigh_site, dist, index) in \
                            enumerate(neighs_dists_indeces):
                        dhs_inter_neigh = radii_sites[natoms-1]+radii_sites[index]
                        # Multiplication and subsequent division with 20
                        # makes it possible that we identify neighbors that
                        # are equally far from the interstitial trial site
                        # to a reasonable numerical precision of 1/20 = 0.05
                        # (typically, Angstrom).
                        rel_dist.append(float(int(20.0 * dist / dhs_inter_neigh)) / 20.0)
                        if rel_dist[i_rel_dist] < drel_overlap:
                            skip_this = True
                            break

                    if skip_this:
                        continue

                    # Sort relative distances and determine unique values.
                    rel_dist_sorted = sorted(rel_dist)
                    change = True
                    while change:
                        change = False
                        for irds, rds in enumerate(rel_dist_sorted):
                            if irds > 0:
                                if fabs(rds - rel_dist_sorted[irds-1]) < 0.049:
                                    del rel_dist_sorted[irds]
                                    change = True
                                    break

                    # Create correspondence map between list storing
                    # unique sorted relative distances (rel_dist_sorted) and
                    # list storing neighbors, distances, and
                    # indeces of sites in Structure object (neighs_dists_indeces).
                    map_urds_ndi = []
                    ncheck = 0
                    for iurds, urds in enumerate(rel_dist_sorted):
                        map_urds_ndi.append([])
                        for ird, rd in enumerate(rel_dist):
                            if fabs(urds - rd) < 0.049:
                                ncheck = ncheck + 1
                                map_urds_ndi[iurds].append(ird)
                        if len(map_urds_ndi[iurds]) < 1: # xxx new May 24, 2016
                            raise RuntimeError("expected at least one entry"
                                    " in rel_dist per unique relative distance")
                    if ncheck != len(neighs_dists_indeces):
                        raise RuntimeError("found incorrect number ({}) of entries"
                                " (expected: {})".format(ncheck,
                                len(neighs_dists_indeces)))

                    site_lists = [[] for i in range(1+nelems)]
                    skip_elem_list = [False for i in range(1+nelems)]
                    last_iurds = -1
                    for iurds, urds in enumerate(rel_dist_sorted):
                        last_iurds = iurds

                        # Put interstitial at first place in first site list.
                        if iurds == 0:
                            # Put interstitial at first place in first site list
                            # for each (elemental sub)lattice.
                            for site_list_elem in site_lists:
                                site_list_elem.append([sc_struct.sites[natoms-1]])

                        # Add all sites associated with previous relative distance.
                        else:
                            for i_site_list_elem, site_list_elem in enumerate(site_lists):
                                site_list_elem.append([])
                                for site in site_list_elem[iurds-1]:
                                    site_list_elem[iurds].append(site)

                        # Add new sites asssociated with the current relative distance.
                        for i in map_urds_ndi[iurds]:
                            neigh_site = neighs_dists_indeces[i][0]
                            if isinstance(neigh_site.specie, Element):
                                this_elem = neigh_site.specie.symbol
                            elif isinstance(neigh_site.specie, Specie):
                                this_elem = neigh_site.specie.element.symbol
                            else:
                                raise RuntimeError("unexpected instance type")
                            if this_elem not in elem_list_plus_inter:
                                raise RuntimeError("this_elem is not in elem_list")
                            if not skip_elem_list[nelems]:
                                site_lists[nelems][iurds].append(neigh_site)
                            for i_this_elem, elem in enumerate(elem_list_plus_inter):
                                if this_elem == elem and not skip_elem_list[i_this_elem]:
                                    site_lists[i_this_elem][iurds].append(neigh_site)

                        # No structure motifs so far with > 10 neighs.
                        # Therefore, cutting off everything
                        # that has more than 10 neighbors.
                        for i_site_list_elem, site_list_elem in enumerate(site_lists):
                            if len(site_list_elem[iurds]) > 10:
                                site_list_elem[iurds] = []
                                skip_elem_list[i_site_list_elem] = True
                        if False not in skip_elem_list:
                            last_iurds = iurds - 1
                            break
                    if last_iurds <= -1:
                        continue

                    # Remove those neighbor lists of sublattices
                    # which have the same number of sites as their preceding list
                    # for that element type, thus, indicating same sites.
                    for i_site_list_elem, site_list_elem in enumerate(site_lists):
                        change = True
                        while change:
                            change = False
                            for iurds in range(1, len(site_list_elem)):
                                if site_list_elem[iurds]:
                                    if len(site_list_elem[iurds]) < len(site_list_elem[iurds-1]):
                                        raise RuntimeError("next nonempty site list"
                                                " should have more than or equal"
                                                " number of sites as previous list.")
                                    if len(site_list_elem[iurds]) == len(site_list_elem[iurds-1]):
                                        if i_site_list_elem == nelems:
                                            raise RuntimeError("there should not be any"
                                                    " consecutive neighbor lists that"
                                                    " have the same number of sites for"
                                                    " the entire lattice.")
                                        for iurds2 in range(iurds-1, len(site_list_elem)-1):
                                            site_list_elem[iurds2] = list(site_list_elem[iurds2+1])
                                        site_list_elem[len(site_list_elem)-1] = []
                                        change = True

                    # Create lists of indeces representing the sites
                    # in a given neighbor site list.
                    indeces_lists = [[] for i in range(1+nelems)]
                    for i_site_list_elem, site_list_elem in enumerate(site_lists):
                        for iurds in range(last_iurds+1):
                            if not site_list_elem[iurds]:
                                break
                            else:
                                indeces_lists[i_site_list_elem].append([])
                                for i_site, site in enumerate(site_list_elem[iurds]):
                                    if site != sc_struct.sites[natoms-1]:
                                        if i_site == 0:
                                            raise RuntimeError("unexpected index of noncentral site!")
                                        indeces_lists[i_site_list_elem][iurds].append(i_site)
                                    elif i_site != 0:
                                        raise RuntimeError("unexpected index of central site!")
                                if len(indeces_lists[i_site_list_elem][iurds])+1 != \
                                        len(site_lists[i_site_list_elem][iurds]):
                                    raise RuntimeError("inconsistency between indeces_lists"
                                                       " and site_lists!")

                    # Check whether this trial position looks like
                    # a coordination pattern-resembling interstitial site.
                    # Because we prefer to get motifs from all-element lattice
                    # we iterate backwardly (i.e., from nelems to 0 (inclusive)).
                    for i_elem in range(nelems,-1, -1):

                        for iurds in range(len(indeces_lists[i_elem])):

                            for it, t in enumerate(self._motif_types):

                                # 1st criterion (very stringent)
                                # to consider the trial position as interstitial site:
                                # does the total coordination number comply with the
                                # target structure motif (tetrahedral --> 4,
                                # octahedral --> 6, bcc --> 8, tet-oct overlay --> 10)?
                                cn_total = len(indeces_lists[i_elem][iurds])

                                if t == "tet" and cn_total != 4:
                                    continue
                                elif t == "oct" and cn_total != 6:
                                    continue
                                elif t == "bcc" and cn_total != 8:
                                    continue
                                elif t == "tetoct" and cn_total != 10:
                                    continue
                                #elif (t == "fcc" or t == "hcp") and cn_total != 12:
                                #    continue

                                # 2nd criterion: does any existing
                                # interstitial site of the current
                                # coordination pattern type
                                # that is within 2*rinter distance
                                # of the current trial site
                                # possess the same chemical environment?
                                # First, determine element-wise CNs.
                                index = -1
                                cns = {}
                                for isite, site in enumerate(site_lists[i_elem][iurds]):
                                    # Skip the interstitial itself.
                                    if isite == 0:
                                        continue
                                    if isinstance(site.specie, Element):
                                        elem = site.specie.symbol
                                    elif isinstance(site.specie, Specie):
                                        elem = site.specie.element.symbol
                                    else:
                                        raise RuntimeError("unexpected instance type")
                                    if elem in cns.keys():
                                        cns[elem] = cns[elem] + 1
                                    else:
                                        cns[elem] = 1
                                # Now, is there already another site
                                # with the same chemical environment
                                # within 2*rinter distance?
                                # If yes, register for potential replacement.
                                for i_intersite, intersite in enumerate(inter_sites[i_elem][t]):
                                    if len(intersite["cns"].keys()) == len(cns.keys()):
                                        for elem, cn in cns.items():
                                            if elem in intersite["cns"].keys():
                                                if cn == intersite["cns"][elem]:
                                                    (dist, jimage)  = \
                                                            lat.get_distance_and_image(
                                                            [a, b, c],
                                                            intersite["frac_coords"])
                                                    if dist < two_rinter:
                                                        index = i_intersite
                                                else:
                                                    break
                                            else:
                                                break

                                if t != "tetoct":
                                    if t == "tet":
                                        ntrials_tet = ntrials_tet + 1
                                    elif t == "oct":
                                        ntrials_oct = ntrials_oct + 1

                                    opvals = ops[self._map_imotif_iop[it][0]].get_order_parameters(
                                            site_lists[i_elem][iurds], 0,
                                            indeces_neighs=indeces_lists[i_elem][iurds])
                                    if not opvals:
                                        raise RuntimeError("unable to compute OPs")
                                    opval = opvals[0]
                                    # Introduce sensible threshold
                                    # to keep the number of potential sites small enough.
                                    if opval > op_threshs[it]:
                                        if not inter_sites[i_elem][t]:
                                            insert_index = 0
                                            inter_sites[i_elem][t].append({})
                                        elif index == -1:
                                            insert_index = len(inter_sites[i_elem][t])
                                            inter_sites[i_elem][t].append({})
                                        else:
                                            if fabs(opval-op_targets[it]) < \
                                                        fabs(inter_sites[i_elem][t][index]["op_value"]-op_targets[it]):
                                                insert_index = index
                                            else:
                                                insert_index = -1
                                        if insert_index > -1:
                                            if insert_index > len(inter_sites[i_elem][t])-1:
                                                raise RuntimeError("")
                                            inter_sites[i_elem][t][insert_index]["frac_coords"] = [a, b, c]
                                            inter_sites[i_elem][t][insert_index]["position"] = r
                                            inter_sites[i_elem][t][insert_index]["cns"] = cns
                                            inter_sites[i_elem][t][insert_index]["op_value"] = opval
                                            inter_sites[i_elem][t][insert_index]["site_list"] = site_lists[i_elem][iurds]

                                # Note that we only search until one combination is successful.
                                else: # if t == "tetoct"
                                    ntrials_tetoct = ntrials_tetoct + 1
                                    for i_tet_indeces, tet_indeces in enumerate(lists_tet_indeces):
                                        qtet_tetmot = ops[self._map_imotif_iop[it][0]].get_order_parameters(
                                                site_lists[i_elem][iurds], 0,
                                                indeces_neighs=tet_indeces)[0]
                                        qoct_tetmot = ops[self._map_imotif_iop[it][1]].get_order_parameters(
                                                site_lists[i_elem][iurds], 0,
                                                indeces_neighs=tet_indeces)[0]
                                        qtet_octmot = ops[self._map_imotif_iop[it][0]].get_order_parameters(
                                                site_lists[i_elem][iurds], 0,
                                                indeces_neighs=lists_oct_indeces[i_tet_indeces])[0]
                                        qoct_octmot = ops[self._map_imotif_iop[it][1]].get_order_parameters(
                                                site_lists[i_elem][iurds], 0,
                                                indeces_neighs=lists_oct_indeces[i_tet_indeces])[0]
                                        if not qtet_tetmot or not qoct_tetmot \
                                                or not qtet_octmot or \
                                                not qoct_octmot:
                                            raise RuntimeError("could not"
                                                    " compute one or more"
                                                    " order parameters.")

                                        # Introduce sensible threshold
                                        # to keep the number of potential sites small enough.
                                        if qtet_tetmot > op_threshs[0] and \
                                                qoct_tetmot < op_threshs[1] and \
                                                qoct_octmot > op_threshs[1] and \
                                                qtet_octmot < op_threshs[0]:
                                            if not inter_sites[i_elem][t]:
                                                insert_index = 0
                                                inter_sites[i_elem][t].append({})
                                            elif index == -1:
                                                insert_index = len(inter_sites[i_elem][t])
                                                inter_sites[i_elem][t].append({})
                                            else:
                                                if fabs(0.5 * (qtet_tetmot + qoct_octmot) - op_targets[it]) < \
                                                            fabs(inter_sites[i_elem][t][index]["op_value"] - op_targets[it]):
                                                    insert_index = index
                                                else:
                                                    insert_index = -1
                                            if insert_index > -1:
                                                if insert_index > len(inter_sites[i_elem][t])-1:
                                                    raise RuntimeError("")
                                                inter_sites[i_elem][t][insert_index]["frac_coords"] = [a, b, c]
                                                inter_sites[i_elem][t][insert_index]["position"] = r
                                                inter_sites[i_elem][t][insert_index]["cns"] = cns
                                                inter_sites[i_elem][t][insert_index]["op_value"] = 0.5 * (qtet_tetmot + qoct_octmot)
                                                inter_sites[i_elem][t][insert_index]["site_list"] = site_lists[i_elem][iurds]

        if write_timings:
            tt1, tc1 = time.time(), time.clock()
            file_time.write("findinters {} {}\n".format(tt1-tt0, tc1-tc0))
            tt0, tc0 = time.time(), time.clock()

        # Check whether any trial site
        # from the sublattice structure-motif search
        # overlaps with any discarded neighbor
        # using drel_overlap * (rinter+rneigh).
        dont_include = []

        for i_elem in range(nelems+1):
            for it, t in enumerate(self._motif_types):
                for i_intersite, intersite in enumerate(inter_sites[i_elem][t]):
                    sc_struct.replace(
                            natoms-1,
                            inter_elem, coords=intersite["frac_coords"],
                            coords_are_cartesian=False)
                    neighs_dists_indeces = sc_struct.get_neighbors(
                            sc_struct.sites[natoms-1],
                            max_radius_sphere, include_index=True)
                    for neighsite, dist, index in neighs_dists_indeces:
                        if dist < (radii_sites[natoms-1]+radii_sites[index]) * drel_overlap:
                            if [i_elem, it, i_intersite] not in dont_include:
                                dont_include.append([i_elem, it, i_intersite])
                                break

        if write_timings:
            tt1, tc1 = time.time(), time.clock()
            file_time.write("pruneatomoverlap {} {}\n".format(tt1-tt0, tc1-tc0))
            tt0, tc0 = time.time(), time.clock()

        # Check whether any two trial sites
        # overlap in the primitive unit cell
        # using drel_overlap * 2*rinter or
        # whether their positions yield
        # similar structure using the StructureMatcher
        # class.  If so, remove the site with the
        # OP value farthest from the respective target
        # value.

        lat = sc_struct.lattice
        overlap = 2.0*radii_sites[natoms-1] * drel_overlap
        sc_struct2 = sc_struct.copy()
        for i_elem1 in range(nelems+1):
            for it1, t1 in enumerate(self._motif_types):
                for i_intersite1, intersite1 in enumerate(inter_sites[i_elem1][t1]):
                    for i_elem2 in range(nelems+1):
                        for it2, t2 in enumerate(self._motif_types):
                            for i_intersite2, intersite2 in enumerate(inter_sites[i_elem2][t2]):
                                if i_elem1 != i_elem2 or \
                                        it1 != it2 or \
                                        i_intersite1 != i_intersite2:
                                    if [i_elem1, it1, i_intersite1] not in dont_include and \
                                            [i_elem2, it2, i_intersite2] not in dont_include:
                                        dist, jimage = lat.get_distance_and_image(
                                                intersite1["frac_coords"], intersite2["frac_coords"])
                                        if dist < overlap:
                                            if fabs(intersite1["op_value"]-op_targets[it1]) < \
                                                    fabs(intersite2["op_value"]-op_targets[it2]):
                                                dont_include.append([i_elem2, it2, i_intersite2])
                                            else:
                                                dont_include.append([i_elem1, it1, i_intersite1])

        if write_timings:
            tt1, tc1 = time.time(), time.clock()
            file_time.write("pruneinteroverlap {} {}\n".format(tt1-tt0, tc1-tc0))
            tt0, tc0 = time.time(), time.clock()

        # Repeating same loop to have the smallest number of sites
        # to be checked for symmetric equivalence
        # because that can take a lot of time.
        for i_elem1 in range(nelems+1):
            for it1, t1 in enumerate(self._motif_types):
                for i_intersite1, intersite1 in enumerate(inter_sites[i_elem1][t1]):
                    for i_elem2 in range(nelems+1):
                        for it2, t2 in enumerate(self._motif_types):
                            for i_intersite2, intersite2 in enumerate(inter_sites[i_elem2][t2]):
                                if i_elem1 != i_elem2 or \
                                        it1 != it2 or \
                                        i_intersite1 != i_intersite2:
                                    if [i_elem1, it1, i_intersite1] not in dont_include and \
                                            [i_elem2, it2, i_intersite2] not in dont_include:
                                        sc_struct.replace(
                                                natoms-1,
                                                inter_elem, coords=intersite1["frac_coords"],
                                                coords_are_cartesian=False)
                                        sc_struct2.replace(
                                                natoms-1,
                                                inter_elem, coords=intersite2["frac_coords"],
                                                coords_are_cartesian=False)
                                        if sgops.are_symmetrically_equivalent(
                                                sc_struct, sc_struct2):
                                            if fabs(intersite1["op_value"]-op_targets[it1]) < \
                                                    fabs(intersite2["op_value"]-op_targets[it2]):
                                                dont_include.append([i_elem2, it2, i_intersite2])
                                            else:
                                                dont_include.append([i_elem1, it1, i_intersite1])

        if write_timings:
            tt1, tc1 = time.time(), time.clock()
            file_time.write("prunesymequiv {} {}\n".format(tt1-tt0, tc1-tc0))
            tt0, tc0 = time.time(), time.clock()

        # Set-up final interstitial site lists, taking into account only
        # those that are not mentioned in dont_include.
        for i_elem in range(nelems+1):
            for it, t in enumerate(self._motif_types):
                for i_intersite, intersite in enumerate(inter_sites[i_elem][t]):
                    if [i_elem, it, i_intersite] not in dont_include:
                        inter_sites_pruned[i_elem][t].append(intersite)

        i = 0
        for i_elem in range(nelems+1):
            for it, t in enumerate(self._motif_types):
                for i_intersite, intersite in enumerate(inter_sites_pruned[i_elem][t]):

                    self._defect_sites.append(
                        PeriodicSite(
                            Element(inter_elem),
                            self._structure.lattice.get_fractional_coords(
                                intersite["position"]),
                            self._structure.lattice,
                            to_unit_cell=False,
                            coords_are_cartesian=False,
                            properties=None))
                    self._defect_types.append(t)
                    self._defect_cns.append(intersite["cns"])
                    # xxx still have to change this
                    self._defect_site_multiplicity.append(1)

                    i = i + 1

        if write_timings:
            file_time.close()


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


