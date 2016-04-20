# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module contains the object used to describe the possible bonded atoms based on a Voronoi analysis
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


from operator import attrgetter

import logging
import numpy as np
import time
from pyhull.voronoi import VoronoiTess
from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite
from monty.json import MSONable
from pymatgen.analysis.structure_analyzer import solid_angle

from pymatgen.analysis.chemenv.utils.chemenv_errors import ChemenvError
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import my_solid_angle
from pymatgen.analysis.chemenv.utils.defs_utils import AdditionalConditions


def from_bson_voronoi_list(bson_nb_voro_list, structure):
    """
    Returns the voronoi_list needed for the VoronoiContainer object from a bson-encoded voronoi_list (composed of
    vlist and bson_nb_voro_list).
    :param vlist: List of voronoi objects
    :param bson_nb_voro_list: List of periodic sites involved in the Voronoi
    :return: The voronoi_list needed for the VoronoiContainer (with PeriodicSites as keys of the dictionary - not
    allowed in the BSON format)
    """
    voronoi_list = [None] * len(bson_nb_voro_list)
    for isite, voro in enumerate(bson_nb_voro_list):
        if voro is None or voro == 'None':
            continue
        voronoi_list[isite] = []
        for psd, dd in voro:
            struct_site = structure[dd['index']]
            periodic_site = PeriodicSite(struct_site._species, struct_site.frac_coords + psd[1],
                                         struct_site._lattice, properties=struct_site._properties)
            voronoi_list[isite].append((periodic_site, dd))
    return voronoi_list


class DetailedVoronoiContainer(MSONable):
    """
    Class used to store the full Voronoi of a given structure.
    """
    AC = AdditionalConditions()
    default_voronoi_cutoff = 10.0

    def __init__(self, structure=None, voronoi_list=None,
                 neighbors_lists=None,
                 voronoi_cutoff=default_voronoi_cutoff, isites=None,
                 weighted_distance_tolerance=1e-5, weighted_angle_tolerance=1e-3,
                 additional_conditions=None, valences=None,
                 maximum_distance_factor=None, minimum_angle_factor=None):
        """
        Constructor for the VoronoiContainer object. Either a structure is given, in which case the Voronoi is
        computed, or the different components of the VoronoiContainer are given (used in the from_dict method)
        :param structure: Structure for which the Voronoi is computed
        :param voronoi_list: List of voronoi polyhedrons for each site
        :param neighbors_list: list of neighbors for each site
        :param voronoi_cutoff: cutoff used for the voronoi
        :param isites: indices of sites for which the Voronoi has to be computed
        :raise: RuntimeError if the Voronoi cannot be constructed
        """
        self.weighted_distance_tolerance = weighted_distance_tolerance
        self.weighted_angle_tolerance = weighted_angle_tolerance
        if additional_conditions is None:
            self.additional_conditions = [self.AC.NONE, self.AC.ONLY_ACB]
        else:
            self.additional_conditions = additional_conditions
        self.valences = valences
        self.maximum_distance_factor = maximum_distance_factor
        self.minimum_angle_factor = minimum_angle_factor
        if isites is None:
            indices = list(range(len(structure)))
        else:
            indices = isites
        self.structure = structure
        logging.info('Setting Voronoi list')
        if voronoi_list is not None:
            self.voronoi_list = voronoi_list
        else:
            self.setup_voronoi_list(indices=indices, voronoi_cutoff=voronoi_cutoff)
        logging.info('Setting neighbors distances and angles')
        t1 = time.clock()
        self.setup_neighbors_distances_and_angles(indices=indices)
        t2 = time.clock()
        logging.info('Neighbors distances and angles set up in {:.2f} seconds'.format(t2-t1))
        if neighbors_lists is None:
            self.setup_neighbors(additional_conditions=self.additional_conditions, valences=self.valences)
        else:
            self.neighbors_lists = neighbors_lists
        logging.info('Setting unique coordinations')
        t1 = time.clock()
        self.setup_unique_coordinations()
        t2 = time.clock()
        logging.info('Unique coordinations set up in {:.2f} seconds'.format(t2-t1))

    def setup_voronoi_list(self, indices, voronoi_cutoff):
        """
        Set up of the voronoi list of neighbours by calling qhull
        :param indices: indices of the sites for which the Voronoi is needed
        :param voronoi_cutoff: Voronoi cutoff for the search of neighbours
        :raise RuntimeError: If an infinite vertex is found in the voronoi construction
        """
        self.voronoi_list = [None] * len(self.structure)
        logging.info('Getting all neighbors in structure')
        struct_neighbors = self.structure.get_all_neighbors(voronoi_cutoff, include_index=True)
        t1 = time.clock()
        logging.info('Setting up Voronoi list :')

        for jj, isite in enumerate(indices):
            logging.info('  - Voronoi analysis for site #{:d} ({:d}/{:d})'.format(isite, jj+1, len(indices)))
            site = self.structure[isite]
            neighbors1 = [(site, 0.0, isite)]
            neighbors1.extend(struct_neighbors[isite])
            distances = [i[1] for i in sorted(neighbors1, key=lambda s: s[1])]
            neighbors = [i[0] for i in sorted(neighbors1, key=lambda s: s[1])]
            qvoronoi_input = [s.coords for s in neighbors]
            voro = VoronoiTess(qvoronoi_input)
            all_vertices = voro.vertices

            results = []
            maxangle = 0.0
            mindist = 10000.0
            for nn, vind in list(voro.ridges.items()):
                if 0 in nn:
                    if 0 in vind:
                        raise RuntimeError("This structure is pathological,"
                                           " infinite vertex in the voronoi "
                                           "construction")

                    facets = [all_vertices[i] for i in vind]
                    try:
                        sa = solid_angle(site.coords, facets)
                    except ValueError:
                        sa = my_solid_angle(site.coords, facets)
                    maxangle = max([sa, maxangle])
                    mindist = min([mindist, distances[nn[1]]])
                    for iii, sss in enumerate(self.structure):
                        if neighbors[nn[1]].is_periodic_image(sss):
                            myindex = iii
                            break
                    results.append((neighbors[nn[1]],
                                    {'angle': sa,
                                     'distance': distances[nn[1]],
                                     'index': myindex}))
            for (nn, dd) in results:
                dd['weighted_angle'] = dd['angle'] / maxangle
                dd['weighted_distance'] = dd['distance'] / mindist
            self.voronoi_list[isite] = results
        t2 = time.clock()
        logging.info('Voronoi list set up in {:.2f} seconds'.format(t2-t1))

    def setup_neighbors_distances_and_angles(self, indices):
        """
        Initializes the angle and distance separations
        :param indices: indices of the sites for which the Voronoi is needed
        """
        self.neighbors_distances = [None] * len(self.structure)
        self.neighbors_weighted_distances = [None] * len(self.structure)
        self.neighbors_angles = [None] * len(self.structure)
        self.neighbors_weighted_angles = [None] * len(self.structure)
        for isite in indices:
            results = self.voronoi_list[isite]
            if results is None:
                continue
            #Initializes neighbors distances and weighted distances groups
            self.neighbors_distances[isite] = []
            self.neighbors_weighted_distances[isite] = []
            weighted_distances = [dd['weighted_distance'] for (nn, dd) in results]
            isorted_distances = np.argsort(weighted_distances)
            #self.neighbors_weighted_distances[isite].append(weighted_distances[isorted_distances[0]])
            self.neighbors_weighted_distances[isite].append({'min': weighted_distances[isorted_distances[0]],
                                                             'max': weighted_distances[isorted_distances[0]]})
            self.neighbors_distances[isite].append({'min': results[isorted_distances[0]][1]['distance'],
                                                    'max': results[isorted_distances[0]][1]['distance']})
            for idist in iter(isorted_distances):
                if self.maximum_distance_factor is not None:
                    if weighted_distances[idist] > self.maximum_distance_factor:
                        self.neighbors_weighted_distances[isite][-1]['max'] = weighted_distances[idist]
                        self.neighbors_distances[isite][-1]['max'] = results[idist][1]['distance']
                        break
                if not np.isclose(weighted_distances[idist], self.neighbors_weighted_distances[isite][-1]['max'],
                                  rtol=0.0, atol=self.weighted_distance_tolerance):
                    self.neighbors_weighted_distances[isite].append({'min': weighted_distances[idist],
                                                                     'max': weighted_distances[idist]})
                    self.neighbors_distances[isite].append({'min': results[idist][1]['distance'],
                                                            'max': results[idist][1]['distance']})
                else:
                    self.neighbors_weighted_distances[isite][-1]['max'] = weighted_distances[idist]
                    self.neighbors_distances[isite][-1]['max'] = results[idist][1]['distance']
            #Initializes neighbors angles and weighted angles groups
            self.neighbors_angles[isite] = []
            self.neighbors_weighted_angles[isite] = []
            weighted_angles = [dd['weighted_angle'] for (nn, dd) in results]
            isorted_angles = np.argsort(weighted_angles)
            self.neighbors_weighted_angles[isite].append({'max': weighted_angles[isorted_angles[0]],
                                                          'min': weighted_angles[isorted_angles[0]]})
            self.neighbors_angles[isite].append({'max': results[isorted_angles[0]][1]['angle'],
                                                 'min': results[isorted_angles[0]][1]['angle']})
            for iang in iter(isorted_angles[1:]):
                if self.minimum_angle_factor is not None:
                    if weighted_angles[iang] < self.minimum_angle_factor:
                        self.neighbors_weighted_angles[isite][-1]['min'] = weighted_angles[iang]
                        self.neighbors_angles[isite][-1]['min'] = results[iang][1]['angle']
                        break
                if not np.isclose(weighted_angles[iang], self.neighbors_weighted_angles[isite][-1]['min'],
                                  rtol=0.0, atol=self.weighted_angle_tolerance):
                    self.neighbors_weighted_angles[isite].append({'max': weighted_angles[iang],
                                                                  'min': weighted_angles[iang]})
                    self.neighbors_angles[isite].append({'max': results[iang][1]['angle'],
                                                         'min': results[iang][1]['angle']})
                else:
                    self.neighbors_weighted_angles[isite][-1]['min'] = weighted_angles[iang]
                    self.neighbors_angles[isite][-1]['min'] = results[iang][1]['angle']

    def setup_neighbors(self, additional_conditions=None, valences=None):
        """
        Compute the list of neighbors for each distfactor/angfactor set of parameters from the voronoi list. The set of
        distfactor/angfactor parameters is a "square" of all combinations of distfactors with angfactors.
        :param distfactors: list of distfactors
        :param angfactors: list of angfactors
        :param only_anion_cation_bonds: Allows only neighbors that are cations (resp. anions) when the current site is
        an anion (resp. a cation).
        :param valences: Valences of all the sites in the structure, needed to check the anion-cation bond when
        only_anion_cation_bonds is set to True.
        :raise: ChemenvError if only_anion_cation_bonds is set to True and valences are not given.
        """
        if additional_conditions is None:
            additional_conditions = self.AC.ALL
        if (self.AC.ONLY_ACB in additional_conditions or self.AC.ONLY_ACB_AND_NO_E2SEB) and valences is None:
            raise ChemenvError('VoronoiContainer', 'setup_neighbors',
                               'Valences are not given while only_anion_cation_bonds are allowed. Cannot continue')
        self.neighbors_lists = [None] * len(self.voronoi_list)
        self.additional_conditions = additional_conditions

        for ivoronoi, voronoi in enumerate(self.voronoi_list):
            if voronoi is None:
                continue
            self.neighbors_lists[ivoronoi] = []
            voronoi_ac = self._precompute_additional_conditions(ivoronoi, voronoi, valences)
            distance_conditions = self._precompute_distance_conditions(ivoronoi, voronoi)
            angle_conditions = self._precompute_angle_conditions(ivoronoi, voronoi)

            for idp, dp_dict in enumerate(self.neighbors_weighted_distances[ivoronoi]):
                self.neighbors_lists[ivoronoi].append([])
                for iap, ap_dict in enumerate(self.neighbors_weighted_angles[ivoronoi]):
                    self.neighbors_lists[ivoronoi][idp].append([])
                    for iac, ac in enumerate(self.additional_conditions):
                        nlist = [(ips, vals['index'], {'weighted_distance': vals['weighted_distance'],
                                                       'weighted_angle': vals['weighted_angle'],
                                                       'distance': vals['distance'],
                                                       'angle': vals['angle']})
                                 for ips, (ps, vals) in enumerate(voronoi)
                                 if (distance_conditions[idp][ips]) and
                                 (angle_conditions[iap][ips]) and
                                 (voronoi_ac[ac][ips])]
                        self.neighbors_lists[ivoronoi][idp][iap].append(nlist)

    def _precompute_additional_conditions(self, ivoronoi, voronoi, valences):
        additional_conditions = {ac: [] for ac in self.additional_conditions}
        for ips, (ps, vals) in enumerate(voronoi):
            for ac in self.additional_conditions:
                additional_conditions[ac].append(self.AC.check_condition(condition=ac, structure=self.structure,
                                                                         parameters={'valences': valences,
                                                                                     'neighbor_index': vals['index'],
                                                                                     'site_index': ivoronoi}))
        return additional_conditions

    def _precompute_distance_conditions(self, ivoronoi, voronoi):
        distance_conditions = []
        for idp, dp_dict in enumerate(self.neighbors_weighted_distances[ivoronoi]):
            distance_conditions.append([])
            dp = dp_dict['max']
            for ips, (ps, vals) in enumerate(voronoi):
                distance_conditions[idp].append(vals['weighted_distance'] <= dp or
                                                np.isclose(vals['weighted_distance'], dp,
                                                           rtol=0.0, atol=self.weighted_distance_tolerance/2.0))
        return distance_conditions

    def _precompute_angle_conditions(self, ivoronoi, voronoi):
        angle_conditions = []
        for iap, ap_dict in enumerate(self.neighbors_weighted_angles[ivoronoi]):
            angle_conditions.append([])
            ap = ap_dict['max']
            for ips, (ps, vals) in enumerate(voronoi):
                angle_conditions[iap].append(vals['weighted_angle'] >= ap or
                                             np.isclose(vals['weighted_angle'], ap,
                                                        rtol=0.0, atol=self.weighted_angle_tolerance/2.0))
        return angle_conditions

    def setup_unique_coordinations(self):
        """
        Setup of the unique coordinations and the mapping of distfactor/angfactor parameters
        to the unique coordinations.
        """
        self._unique_coordinated_neighbors = [None] * len(self.voronoi_list)
        self._unique_coordinated_neighbors_parameters_indices = [None] * len(self.voronoi_list)
        self._parameters_to_unique_coordinated_neighbors_map = [None] * len(self.voronoi_list)
        ncond_params = len(self.additional_conditions)
        for isite, voronoi in enumerate(self.voronoi_list):

            if voronoi is None:
                continue
            ndist_params = len(self.neighbors_distances[isite])
            nang_params = len(self.neighbors_angles[isite])
            self._unique_coordinated_neighbors[isite] = {}
            self._unique_coordinated_neighbors_parameters_indices[isite] = {}
            self._parameters_to_unique_coordinated_neighbors_map[isite] = [None] * ndist_params
            for idp, dp in enumerate(self.neighbors_distances[isite]):
                self._parameters_to_unique_coordinated_neighbors_map[isite][idp] = [None] * nang_params
                for iap, ap in enumerate(self.neighbors_weighted_angles[isite]):
                    self._parameters_to_unique_coordinated_neighbors_map[isite][idp][iap] = [None] * ncond_params
                    for iac, ac in enumerate(self.additional_conditions):
                        cn = len(self.neighbors_lists[isite][idp][iap][iac])
                        if not cn in self._unique_coordinated_neighbors[isite]:
                            self._unique_coordinated_neighbors[isite][cn] = []
                            self._unique_coordinated_neighbors_parameters_indices[isite][cn] = []
                        indices_array = [nlist[1] for nlist in self.neighbors_lists[isite][idp][iap][iac]]
                        nlist_dict_array = [nlist[2] for nlist in self.neighbors_lists[isite][idp][iap][iac]]
                        ps_array = sorted([self.voronoi_list[isite][nlist[0]][0] for nlist in
                                           self.neighbors_lists[isite][idp][iap][iac]],
                                          key=attrgetter('a', 'b', 'c'))
                        found = False
                        for i_cn_neighblist, nlist_tuple in enumerate(self._unique_coordinated_neighbors[isite][cn]):
                            if indices_array == nlist_tuple[1]:
                                self._parameters_to_unique_coordinated_neighbors_map[isite][idp][iap][iac] = [cn,
                                                                                                              i_cn_neighblist]
                                (self._unique_coordinated_neighbors_parameters_indices[isite][cn][i_cn_neighblist].
                                 append((idp, iap, iac)))
                                found = True
                                break
                        if not found:
                            self._unique_coordinated_neighbors[isite][cn].append((ps_array, indices_array,
                                                                                  nlist_dict_array))
                            i_parameters_list = [(idp, iap, iac)]
                            self._unique_coordinated_neighbors_parameters_indices[isite][cn].append(i_parameters_list)
                            self._parameters_to_unique_coordinated_neighbors_map[isite][idp][iap][iac] = [cn, 0]

    def unique_coordinated_neighbors(self, isite=None, cn_map=None):
        if isite is None:
            return self._unique_coordinated_neighbors
        else:
            if cn_map is None:
                return self._unique_coordinated_neighbors[isite]
            else:
                return self._unique_coordinated_neighbors[isite][cn_map[0]][cn_map[1]]

    @property
    def parameters_to_unique_coordinated_neighbors_map(self):
        return self._parameters_to_unique_coordinated_neighbors_map

    def angles(self, isite, include_index=True):
        if not include_index:
            return [result[1]['angle'] for result in self.voronoi_list[isite]]
        else:
            return [(result[1]['angle'], result[1]['index']) for result in self.voronoi_list[isite]]

    def satisfy_condition(self, isite, cn, i_coordn_neighb, additional_condition):
        parameters_indices = self._unique_coordinated_neighbors_parameters_indices[isite][cn][i_coordn_neighb]
        additional_conditions = [pp[2] for pp in parameters_indices]
        return additional_condition in additional_conditions

    def neighbors_map(self, isite, distfactor, angfactor, additional_condition):
        if self.neighbors_weighted_distances[isite] is None:
            return None
        dist_where = np.argwhere(np.array([wd['min'] for wd in self.neighbors_weighted_distances[isite]]) <= distfactor)
        if len(dist_where) == 0:
            return None
        idist = dist_where[-1][0]
        ang_where = np.argwhere(np.array([wa['max'] for wa in self.neighbors_weighted_angles[isite]]) >= angfactor)
        if len(ang_where) == 0:
            return None
        iang = ang_where[0][0]
        if self.additional_conditions.count(additional_condition) != 1:
            return None
        i_additional_condition = self.additional_conditions.index(additional_condition)
        return {'i_distfactor': idist, 'i_angfactor': iang, 'i_additional_condition': i_additional_condition}

    def neighbors_surfaces(self, isite, surface_calculation_type=None, max_dist=2.0):
        if self.voronoi_list[isite] is None:
            return None
        #surfaces = np.zeros((len(self.neighbors_weighted_distances), len(self.neighbors_weighted_angles)), np.float)
        bounds_and_limits = self.voronoi_parameters_bounds_and_limits(isite, surface_calculation_type, max_dist)
        distance_bounds = bounds_and_limits['distance_bounds']
        angle_bounds = bounds_and_limits['angle_bounds']
        surfaces = np.zeros((len(distance_bounds), len(angle_bounds)), np.float)
        for idp in range(len(distance_bounds) - 1):
            this_dist_plateau = distance_bounds[idp + 1] - distance_bounds[idp]
            for iap in range(len(angle_bounds) - 1):
                this_ang_plateau = angle_bounds[iap + 1] - angle_bounds[iap]
                surfaces[idp][iap] = np.absolute(this_dist_plateau*this_ang_plateau)
        return surfaces

    def maps_with_condition(self, isite, additional_condition, return_parameter_indices=False):
        cn_maps = []
        ucnpi_isite = self._unique_coordinated_neighbors_parameters_indices[isite]
        if ucnpi_isite is None:
            return None
        for cn, coordnbs_list in ucnpi_isite.items():
            for i_coordnbs, coordnbs in enumerate(coordnbs_list):
                if self.satisfy_condition(isite, cn, i_coordnbs, additional_condition):
                    cn_maps.append((cn, i_coordnbs))
        result = {'cn_maps': cn_maps}
        if return_parameter_indices:
            parameter_indices = []
            for (cn, i_coordnbs) in cn_maps:
                cn_map_parameter_indices = []
                for params in self._unique_coordinated_neighbors_parameters_indices[isite][cn][i_coordnbs]:
                    if params[2] == additional_condition:
                        cn_map_parameter_indices.append(params)
                parameter_indices.append(cn_map_parameter_indices)
            result['parameter_indices'] = parameter_indices
        return result

    def maps_and_surface_vertices(self, isite, additional_condition=AC.ONLY_ACB, plot_type=None, max_dist=2.0):
        cn_maps_parameter_indices = self.maps_with_condition(isite=isite, additional_condition=additional_condition,
                                                             return_parameter_indices=True)
        if cn_maps_parameter_indices is None:
            return None
        bounds_and_limits = self.voronoi_parameters_bounds_and_limits(isite, plot_type, max_dist)
        vertices_dist_ang_indices_list = []
        vertices_dist_ang_list = []
        text_info_dist_ang_list = []
        for i_cn_map, cn_map in enumerate(cn_maps_parameter_indices['cn_maps']):
            parameter_indices_list = cn_maps_parameter_indices['parameter_indices'][i_cn_map]
            vertices_dist_ang_indices = self._get_vertices_dist_ang_indices(parameter_indices_list)
            vertices_dist_ang = []
            idist, iang = vertices_dist_ang_indices[0]
            dist = bounds_and_limits['distance_bounds'][idist]
            ang = bounds_and_limits['angle_bounds'][iang]
            vertices_dist_ang.append([dist, ang])
            idist, iang = vertices_dist_ang_indices[1]
            dist = bounds_and_limits['distance_bounds'][idist+1]
            ang = bounds_and_limits['angle_bounds'][iang]
            vertices_dist_ang.append([dist, ang])
            idist, iang = vertices_dist_ang_indices[2]
            dist = bounds_and_limits['distance_bounds'][idist+1]
            ang = bounds_and_limits['angle_bounds'][iang]
            vertices_dist_ang.append([dist, ang])
            idist, iang = vertices_dist_ang_indices[3]
            dist = bounds_and_limits['distance_bounds'][idist+1]
            ang = bounds_and_limits['angle_bounds'][iang]
            vertices_dist_ang.append([dist, ang])
            idist, iang = vertices_dist_ang_indices[4]
            dist = bounds_and_limits['distance_bounds'][idist+1]
            ang = bounds_and_limits['angle_bounds'][iang+1]
            vertices_dist_ang.append([dist, ang])
            idist, iang = vertices_dist_ang_indices[5]
            dist = bounds_and_limits['distance_bounds'][idist]
            ang = bounds_and_limits['angle_bounds'][iang+1]
            vertices_dist_ang.append([dist, ang])
            vertices_dist_ang_indices_list.append(vertices_dist_ang_indices)
            vertices_dist_ang_list.append(vertices_dist_ang)
            text_info_dist_ang = ((vertices_dist_ang[2][0] + vertices_dist_ang[5][0]) / 2.0,
                                  (vertices_dist_ang[2][1] + vertices_dist_ang[5][1]) / 2.0)
            text_info_dist_ang_list.append(text_info_dist_ang)
        result = {'cn_maps': cn_maps_parameter_indices['cn_maps'],
                  'bounds_and_limits': bounds_and_limits,
                  'vertices_dist_ang_indices': vertices_dist_ang_indices_list,
                  'vertices_dist_ang': vertices_dist_ang_list,
                  'text_info_dist_ang': text_info_dist_ang_list
                  }
        return result

    @staticmethod
    def _get_vertices_dist_ang_indices(parameter_indices_list):
        pp0 = [pp[0] for pp in parameter_indices_list]
        pp1 = [pp[1] for pp in parameter_indices_list]
        min_idist = min(pp0)
        min_iang = min(pp1)
        max_idist = max(pp0)
        max_iang = max(pp1)
        i_min_angs = np.argwhere(np.array(pp1) == min_iang)
        i_max_dists = np.argwhere(np.array(pp0) == max_idist)
        pp0_at_min_iang = [pp0[ii[0]] for ii in i_min_angs]
        pp1_at_max_idist = [pp1[ii[0]] for ii in i_max_dists]
        max_idist_at_min_iang = max(pp0_at_min_iang)
        min_iang_at_max_idist = min(pp1_at_max_idist)

        p1 = (min_idist, min_iang)
        p2 = (max_idist_at_min_iang, min_iang)
        p3 = (max_idist_at_min_iang, min_iang_at_max_idist)
        p4 = (max_idist, min_iang_at_max_idist)
        p5 = (max_idist, max_iang)
        p6 = (min_idist, max_iang)

        return [p1, p2, p3, p4, p5, p6]

    def maps_and_surfaces(self, isite, surface_calculation_type=None, max_dist=2.0, additional_conditions=None):
        if self.voronoi_list[isite] is None:
            return None
        if additional_conditions is None:
            additional_conditions = [self.AC.ONLY_ACB]
        surfaces = self.neighbors_surfaces(isite=isite, surface_calculation_type=surface_calculation_type,
                                           max_dist=max_dist)
        maps_and_surfaces = []
        for cn, value in self._unique_coordinated_neighbors_parameters_indices[isite].items():
            for imap, list_parameters_indices in enumerate(value):
                thissurf = 0.0
                for (idp, iap, iacb) in list_parameters_indices:
                    if iacb in additional_conditions:
                        thissurf += surfaces[idp, iap]
                maps_and_surfaces.append({'map': (cn, imap), 'surface': thissurf,
                                          'parameters_indices': list_parameters_indices})
        return maps_and_surfaces

    def get_neighbors(self, isite, neighbors_map):
        """
        Returns the list of neighbours of a given site given the neighbors_map indices (index of the distance and
        angle factors as well as index of the only_anion_cation_bonds). This will usually not be used outside Voronoi.
        :param isite: Index of the site for which the neighbors have to be given.
        :param neighbors_map: Mapping of the distance/angle/additional_condition as a dict with keys "i_distfactor",
                              "i_angfactor" and "i_additional_condition" and corresponding values.
        :return: Neighbors for this neighbors_map
        """
        this_site_this_map_neighbors_list = (self.neighbors_lists
                                             [isite]
                                             [neighbors_map['i_distfactor']]
                                             [neighbors_map['i_angfactor']]
                                             [neighbors_map['i_additional_condition']])
        neighbors = [self.voronoi_list[isite][nlist[0]][0] for nlist in this_site_this_map_neighbors_list]
        return neighbors

    def neighbors(self, isite, distfactor, angfactor, additional_condition):
        neighbors_map = self.neighbors_map(isite=isite, distfactor=distfactor, angfactor=angfactor,
                                           additional_condition=additional_condition)
        if neighbors_map is None:
            return []
        return self.get_neighbors(isite=isite, neighbors_map=neighbors_map)

    def get_coordination_numbers_figure(self, isite, plot_type=None, title='Coordination numbers', max_dist=2.0,
                                        figsize=None):
        """
        Plotting of the coordination numbers of a given site for all the distfactor/angfactor parameters. If the
        chemical environments are given, a color map is added to the plot, with the lowest continuous symmetry measure
        as the value for the color of that distfactor/angfactor set.
        :param isite: Index of the site for which the plot has to be done
        :param plot_type: How to plot the coordinations
        :param title: Title for the figure
        :param max_dist: Maximum distance to be plotted when the plotting of the distance is set to 'initial_normalized'
                         or 'initial_real' (Warning: this is not the same meaning in both cases! In the first case,
                         the closest atom lies at a "normalized" distance of 1.0 so that 2.0 means refers to this
                         normalized distance while in the second case, the real distance is used)
        :param figsize: Size of the figure to be plotted
        :return: The figure object to be plotted or saved to file
        """
        try:
            import matplotlib.pyplot as mpl
            from matplotlib import cm
            from matplotlib.colors import Normalize, LinearSegmentedColormap, ListedColormap
            from matplotlib.patches import Rectangle, Polygon
        except ImportError:
            print('Plotting Chemical Environments requires matplotlib ... exiting "plot" function')
            return

        #Initializes the figure
        if figsize is None:
            fig = mpl.figure()
        else:
            fig = mpl.figure(figsize=figsize)
        subplot = fig.add_subplot(111)

        #Initializes the distance and angle parameters
        bounds_and_limits = self.voronoi_parameters_bounds_and_limits(isite, plot_type, max_dist)
        if bounds_and_limits is None:
            return None
        distance_bounds = bounds_and_limits['distance_bounds']
        angle_bounds = bounds_and_limits['angle_bounds']
        dist_limits = bounds_and_limits['distance_limits']
        ang_limits = bounds_and_limits['angle_limits']

        #Plot the rectangles and coordinations
        for idp in range(len(distance_bounds) - 1):
            this_dist_plateau = distance_bounds[idp + 1] - distance_bounds[idp]
            #print('Dist Plateau : ', this_dist_plateau)
            for iap in range(len(angle_bounds) - 1):
                this_ang_plateau = angle_bounds[iap + 1] - angle_bounds[iap]
                #print('Ang Plateau : ', this_ang_plateau)
                #print('Rectangle from xy = ', distance_bounds[idp], ' ', angle_bounds[iap], 'of width ', this_dist_plateau, 'and height ', this_ang_plateau)
                r = Rectangle((distance_bounds[idp], angle_bounds[iap]),
                              this_dist_plateau,
                              this_ang_plateau, edgecolor='k', facecolor='b')
                subplot.annotate('{:d}'.format(len(self.neighbors_lists[isite][idp][iap][0])),
                                 xy=(distance_bounds[idp] + this_dist_plateau / 2.0,
                                     angle_bounds[iap] + this_ang_plateau / 2.0),
                                 ha='center', va='center', color='y', fontsize='x-small')
                subplot.add_patch(r)
        title += '\nDist: {}, Ang: {}'.format(plot_type['distance_parameter'][0], plot_type['angle_parameter'][0])
        subplot.set_title(title)
        subplot.set_xlabel('Distance factor')
        subplot.set_ylabel('Angle factor')
        subplot.set_xlim(dist_limits)
        subplot.set_ylim(ang_limits)
        return fig

    def voronoi_parameters_bounds_and_limits(self, isite, plot_type, max_dist):
        #Initializes the distance and angle parameters
        if self.voronoi_list[isite] is None:
            return None
        if plot_type is None:
            plot_type = {'distance_parameter': ('initial_inverse_opposite', None),
                         'angle_parameter': ('initial_opposite', None)}
        if plot_type['distance_parameter'][0] == 'initial_normalized':
            #dd = [dist['min'] for dist in self.neighbors_weighted_distances[isite] if dist['min'] <= max_dist]
            dd = [dist['min'] for dist in self.neighbors_weighted_distances[isite]]
        else:
            dd = [dist['min'] for dist in self.neighbors_weighted_distances[isite]]
        dd[0] = 1.0
        if plot_type['distance_parameter'][0] == 'initial_normalized':
            dd.append(max_dist)
            distance_bounds = np.array(dd)
            dist_limits = [1.0, max_dist]
        elif plot_type['distance_parameter'][0] == 'initial_inverse_opposite':
            ddinv = [1.0 / dist for dist in dd]
            ddinv.append(0.0)
            distance_bounds = np.array([1.0 - invdist for invdist in ddinv])
            dist_limits = [0.0, 1.0]
        elif plot_type['distance_parameter'][0] == 'initial_inverse3_opposite':
            ddinv = [1.0 / dist**3.0 for dist in dd]
            ddinv.append(0.0)
            distance_bounds = np.array([1.0 - invdist for invdist in ddinv])
            dist_limits = [0.0, 1.0]
        else:
            raise NotImplementedError('Plotting type "{}" '
                                      'for the distance is not implemented'.format(plot_type['distance_parameter']))
        if plot_type['angle_parameter'][0] == 'initial_normalized':
            aa = [0.0]
            aa.extend([ang['max'] for ang in self.neighbors_weighted_angles[isite]])
            angle_bounds = np.array(aa)
        elif plot_type['angle_parameter'][0] == 'initial_opposite':
            aa = [0.0]
            aa.extend([ang['max'] for ang in self.neighbors_weighted_angles[isite]])
            aa = [1.0 - ang for ang in aa]
            angle_bounds = np.array(aa)
        else:
            raise NotImplementedError('Plotting type "{}" '
                                      'for the angle is not implemented'.format(plot_type['angle_parameter']))
        ang_limits = [0.0, 1.0]
        return {'distance_bounds': distance_bounds, 'distance_limits': dist_limits,
                'angle_bounds': angle_bounds, 'angle_limits': ang_limits}

    def save_coordination_numbers_figure(self, isite, imagename='image.png', plot_type=None,
                                         title='Coordination numbers', max_dist=2.0,
                                         figsize=None):
        fig = self.get_coordination_numbers_figure(isite=isite, plot_type=plot_type, title=title, max_dist=max_dist,
                                                   figsize=figsize)
        if fig is None:
            return
        fig.savefig(imagename)

    def plot_coordination_numbers(self, isite, plot_type=None, title='Coordination numbers', max_dist=2.0,
                                  figsize=None):
        """
        Plotting of the coordination numbers of a given site for all the distfactor/angfactor parameters. If the
        chemical environments are given, a color map is added to the plot, with the lowest continuous symmetry measure
        as the value for the color of that distfactor/angfactor set.
        :param isite: Index of the site for which the plot has to be done
        :param plot_type: How to plot the coordinations
        :param title: Title for the figure
        :param max_dist: Maximum distance to be plotted when the plotting of the distance is set to 'initial_normalized'
                         or 'initial_real' (Warning: this is not the same meaning in both cases! In the first case,
                         the closest atom lies at a "normalized" distance of 1.0 so that 2.0 means refers to this
                         normalized distance while in the second case, the real distance is used)
        :param figsize: Size of the figure to be plotted
        :return: Nothing returned, just plot the figure
        """
        fig = self.get_coordination_numbers_figure(isite=isite, plot_type=plot_type, title=title, max_dist=max_dist,
                                                   figsize=figsize)
        if fig is None:
            return
        fig.show()

    def unique_coordinations(self, isite):
        """
        Returns all the unique coordinations of a given site (two different sets of distfactor/angfactor parameters
        can lead to the same coordination).
        :param isite: Site for which the unique coordinations are needed
        :return: unique coordinations for this site.
        """
        return self._unique_coordinated_neighbors[isite]

    def unique_coordinated_neighbors_parameters_indices(self, isite):
        return self._unique_coordinated_neighbors_parameters_indices[isite]

    def __eq__(self, other):
        return (self.weighted_angle_tolerance == other.weighted_angle_tolerance and
                self.weighted_distance_tolerance == other.weighted_distance_tolerance and
                self.additional_conditions == other.additional_conditions and
                self.valences == other.valences and
                self.voronoi_list == other.voronoi_list and
                self.structure == other.structure)

    def to_bson_voronoi_list(self):
        """
        Transforms the voronoi_list into a vlist + bson_nb_voro_list, that are BSON-encodable.
        :return: [vlist, bson_nb_voro_list], to be used in the as_dict method
        """
        bson_nb_voro_list = [None] * len(self.voronoi_list)
        for ivoro, voro in enumerate(self.voronoi_list):
            if voro is None or voro == 'None':
                continue
            site_voro = []
            for (ps, dd) in voro:
                #site_voro.append([ps.as_dict(), dd]) [float(c) for c in self._fcoords]
                diff = ps._fcoords - self.structure[dd['index']]._fcoords
                site_voro.append([[dd['index'], [float(c) for c in diff]],
                                  dd])
            bson_nb_voro_list[ivoro] = site_voro
        return bson_nb_voro_list

    def as_dict(self):
        """
        Bson-serializable dict representation of the VoronoiContainer.
        :return: dictionary that is BSON-encodable
        """
        bson_nb_voro_list = self.to_bson_voronoi_list()
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "bson_nb_voro_list": bson_nb_voro_list,
                "neighbors_lists": self.neighbors_lists,
                "structure": self.structure.as_dict(),
                "weighted_angle_tolerance": self.weighted_angle_tolerance,
                "weighted_distance_tolerance": self.weighted_distance_tolerance,
                "additional_conditions": self.additional_conditions,
                "valences": self.valences,
                "maximum_distance_factor": self.maximum_distance_factor,
                "minimum_angle_factor": self.minimum_angle_factor}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the VoronoiContainer object from a dict representation of the VoronoiContainer created using
        the as_dict method.
        :param d: dict representation of the VoronoiContainer object
        :return: VoronoiContainer object
        """
        structure = Structure.from_dict(d['structure'])
        voronoi_list = from_bson_voronoi_list(d['bson_nb_voro_list'], structure)
        neighbors_lists = d['neighbors_lists'] if 'neighbors_lists' in d else None
        maximum_distance_factor = d['maximum_distance_factor'] if 'maximum_distance_factor' in d else None
        minimum_angle_factor = d['minimum_angle_factor'] if 'minimum_angle_factor' in d else None
        return cls(structure=structure, voronoi_list=voronoi_list,
                   neighbors_lists=neighbors_lists,
                   weighted_angle_tolerance=d['weighted_angle_tolerance'],
                   weighted_distance_tolerance=d['weighted_distance_tolerance'],
                   additional_conditions=d['additional_conditions'],
                   valences=d['valences'],
                   maximum_distance_factor=maximum_distance_factor,
                   minimum_angle_factor=minimum_angle_factor)