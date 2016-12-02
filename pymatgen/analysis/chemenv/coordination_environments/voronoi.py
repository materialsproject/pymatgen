# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging
import numpy as np
import time
from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite
from monty.json import MSONable
from pymatgen.analysis.structure_analyzer import solid_angle
from scipy.spatial import Voronoi

from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import my_solid_angle
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import get_lower_and_upper_f
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import rectangle_surface_intersection
from pymatgen.analysis.chemenv.utils.defs_utils import AdditionalConditions

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

def from_bson_voronoi_list2(bson_nb_voro_list2, structure):
    """
    Returns the voronoi_list needed for the VoronoiContainer object from a bson-encoded voronoi_list (composed of
    vlist and bson_nb_voro_list).
    :param vlist: List of voronoi objects
    :param bson_nb_voro_list: List of periodic sites involved in the Voronoi
    :return: The voronoi_list needed for the VoronoiContainer (with PeriodicSites as keys of the dictionary - not
    allowed in the BSON format)
    """
    voronoi_list = [None] * len(bson_nb_voro_list2)
    for isite, voro in enumerate(bson_nb_voro_list2):
        if voro is None or voro == 'None':
            continue
        voronoi_list[isite] = []
        for psd, dd in voro:
            struct_site = structure[dd['index']]
            periodic_site = PeriodicSite(struct_site._species, struct_site.frac_coords + psd[1],
                                         struct_site._lattice, properties=struct_site._properties)
            dd['site'] = periodic_site
            voronoi_list[isite].append(dd)
    return voronoi_list


class DetailedVoronoiContainer(MSONable):
    """
    Class used to store the full Voronoi of a given structure.
    """
    AC = AdditionalConditions()
    default_voronoi_cutoff = 10.0

    def __init__(self, structure=None, voronoi_list=None, voronoi_list2=None,
                 # neighbors_lists=None,
                 voronoi_cutoff=default_voronoi_cutoff, isites=None,
                 normalized_distance_tolerance=1e-5, normalized_angle_tolerance=1e-3,
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
        self.normalized_distance_tolerance = normalized_distance_tolerance
        self.normalized_angle_tolerance = normalized_angle_tolerance
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
        if voronoi_list2 is not None:
            self.voronoi_list2 = voronoi_list2
        else:
            self.setup_voronoi_list(indices=indices, voronoi_cutoff=voronoi_cutoff)
        logging.info('Setting neighbors distances and angles')
        t1 = time.clock()
        self.setup_neighbors_distances_and_angles(indices=indices)
        t2 = time.clock()
        logging.info('Neighbors distances and angles set up in {:.2f} seconds'.format(t2-t1))

    def setup_voronoi_list(self, indices, voronoi_cutoff):
        """
        Set up of the voronoi list of neighbours by calling qhull
        :param indices: indices of the sites for which the Voronoi is needed
        :param voronoi_cutoff: Voronoi cutoff for the search of neighbours
        :raise RuntimeError: If an infinite vertex is found in the voronoi construction
        """
        self.voronoi_list2 = [None] * len(self.structure)
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
            voro = Voronoi(points=qvoronoi_input, qhull_options="o Fv")
            all_vertices = voro.vertices

            results2 = []
            maxangle = 0.0
            mindist = 10000.0
            for iridge, ridge_points in enumerate(voro.ridge_points):
                if 0 in ridge_points:
                    ridge_vertices_indices = voro.ridge_vertices[iridge]
                    if -1 in ridge_vertices_indices:
                        raise RuntimeError("This structure is pathological,"
                                           " infinite vertex in the voronoi "
                                           "construction")

                    ridge_point2 = max(ridge_points)
                    facets = [all_vertices[i] for i in ridge_vertices_indices]
                    try:
                        sa = solid_angle(site.coords, facets)
                    except ValueError:
                        sa = my_solid_angle(site.coords, facets)
                    maxangle = max([sa, maxangle])

                    mindist = min([mindist, distances[ridge_point2]])
                    for iii, sss in enumerate(self.structure):
                        if neighbors[ridge_point2].is_periodic_image(sss):
                            myindex = iii
                            break
                    results2.append({'site': neighbors[ridge_point2],
                                     'angle': sa,
                                     'distance': distances[ridge_point2],
                                     'index': myindex})
            for dd in results2:
                dd['normalized_angle'] = dd['angle'] / maxangle
                dd['normalized_distance'] = dd['distance'] / mindist
            self.voronoi_list2[isite] = results2
        t2 = time.clock()
        logging.info('Voronoi list set up in {:.2f} seconds'.format(t2-t1))

    def setup_neighbors_distances_and_angles(self, indices):
        """
        Initializes the angle and distance separations
        :param indices: indices of the sites for which the Voronoi is needed
        """
        self.neighbors_distances = [None] * len(self.structure)
        self.neighbors_normalized_distances = [None] * len(self.structure)
        self.neighbors_angles = [None] * len(self.structure)
        self.neighbors_normalized_angles = [None] * len(self.structure)
        for isite in indices:
            results = self.voronoi_list2[isite]
            if results is None:
                continue
            #Initializes neighbors distances and normalized distances groups
            self.neighbors_distances[isite] = []
            self.neighbors_normalized_distances[isite] = []
            normalized_distances = [nb_dict['normalized_distance'] for nb_dict in results]
            isorted_distances = np.argsort(normalized_distances)
            self.neighbors_normalized_distances[isite].append({'min': normalized_distances[isorted_distances[0]],
                                                               'max': normalized_distances[isorted_distances[0]]})
            self.neighbors_distances[isite].append({'min': results[isorted_distances[0]]['distance'],
                                                    'max': results[isorted_distances[0]]['distance']})
            icurrent = 0
            nb_indices = {int(isorted_distances[0])}
            dnb_indices = {int(isorted_distances[0])}
            for idist in iter(isorted_distances):
                wd = normalized_distances[idist]
                if self.maximum_distance_factor is not None:
                    if wd > self.maximum_distance_factor:
                        self.neighbors_normalized_distances[isite][icurrent]['nb_indices'] = list(nb_indices)
                        self.neighbors_distances[isite][icurrent]['nb_indices'] = list(nb_indices)
                        self.neighbors_normalized_distances[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                        self.neighbors_distances[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                        break
                if np.isclose(wd, self.neighbors_normalized_distances[isite][icurrent]['max'],
                              rtol=0.0, atol=self.normalized_distance_tolerance):
                    self.neighbors_normalized_distances[isite][icurrent]['max'] = wd
                    self.neighbors_distances[isite][icurrent]['max'] = results[idist]['distance']
                    dnb_indices.add(int(idist))
                else:
                    self.neighbors_normalized_distances[isite][icurrent]['nb_indices'] = list(nb_indices)
                    self.neighbors_distances[isite][icurrent]['nb_indices'] = list(nb_indices)
                    self.neighbors_normalized_distances[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                    self.neighbors_distances[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                    dnb_indices = {int(idist)}
                    self.neighbors_normalized_distances[isite].append({'min': wd,
                                                                       'max': wd})
                    self.neighbors_distances[isite].append({'min': results[idist]['distance'],
                                                            'max': results[idist]['distance']})
                    icurrent += 1
                nb_indices.add(int(idist))
            else:
                self.neighbors_normalized_distances[isite][icurrent]['nb_indices'] = list(nb_indices)
                self.neighbors_distances[isite][icurrent]['nb_indices'] = list(nb_indices)
                self.neighbors_normalized_distances[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                self.neighbors_distances[isite][icurrent]['dnb_indices'] = list(dnb_indices)
            for idist in range(len(self.neighbors_distances[isite]) - 1):
                dist_dict = self.neighbors_distances[isite][idist]
                dist_dict_next = self.neighbors_distances[isite][idist+1]
                dist_dict['next'] = dist_dict_next['min']
                ndist_dict = self.neighbors_normalized_distances[isite][idist]
                ndist_dict_next = self.neighbors_normalized_distances[isite][idist + 1]
                ndist_dict['next'] = ndist_dict_next['min']
            if self.maximum_distance_factor is not None:
                dfact = self.maximum_distance_factor
            else:
                dfact = self.default_voronoi_cutoff / self.neighbors_distances[isite][0]['min']
            self.neighbors_normalized_distances[isite][-1]['next'] = dfact
            self.neighbors_distances[isite][-1]['next'] = dfact * self.neighbors_distances[isite][0]['min']
            #Initializes neighbors angles and normalized angles groups
            self.neighbors_angles[isite] = []
            self.neighbors_normalized_angles[isite] = []
            normalized_angles = [nb_dict['normalized_angle'] for nb_dict in results]
            isorted_angles = np.argsort(normalized_angles)[::-1]
            self.neighbors_normalized_angles[isite].append({'max': normalized_angles[isorted_angles[0]],
                                                            'min': normalized_angles[isorted_angles[0]]})
            self.neighbors_angles[isite].append({'max': results[isorted_angles[0]]['angle'],
                                                 'min': results[isorted_angles[0]]['angle']})
            icurrent = 0
            nb_indices = {int(isorted_angles[0])}
            dnb_indices = {int(isorted_angles[0])}
            for iang in iter(isorted_angles):
                wa = normalized_angles[iang]
                if self.minimum_angle_factor is not None:
                    if wa < self.minimum_angle_factor:
                        self.neighbors_normalized_angles[isite][icurrent]['nb_indices'] = list(nb_indices)
                        self.neighbors_angles[isite][icurrent]['nb_indices'] = list(nb_indices)
                        self.neighbors_normalized_angles[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                        self.neighbors_angles[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                        break
                if np.isclose(wa, self.neighbors_normalized_angles[isite][icurrent]['min'],
                              rtol=0.0, atol=self.normalized_angle_tolerance):
                    self.neighbors_normalized_angles[isite][icurrent]['min'] = wa
                    self.neighbors_angles[isite][icurrent]['min'] = results[iang]['angle']
                    dnb_indices.add(int(iang))
                else:
                    self.neighbors_normalized_angles[isite][icurrent]['nb_indices'] = list(nb_indices)
                    self.neighbors_angles[isite][icurrent]['nb_indices'] = list(nb_indices)
                    self.neighbors_normalized_angles[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                    self.neighbors_angles[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                    dnb_indices = {int(iang)}
                    self.neighbors_normalized_angles[isite].append({'max': wa,
                                                                    'min': wa})
                    self.neighbors_angles[isite].append({'max': results[iang]['angle'],
                                                         'min': results[iang]['angle']})
                    icurrent += 1
                nb_indices.add(int(iang))
            else:
                self.neighbors_normalized_angles[isite][icurrent]['nb_indices'] = list(nb_indices)
                self.neighbors_angles[isite][icurrent]['nb_indices'] = list(nb_indices)
                self.neighbors_normalized_angles[isite][icurrent]['dnb_indices'] = list(dnb_indices)
                self.neighbors_angles[isite][icurrent]['dnb_indices'] = list(dnb_indices)
            for iang in range(len(self.neighbors_angles[isite]) - 1):
                ang_dict = self.neighbors_angles[isite][iang]
                ang_dict_next = self.neighbors_angles[isite][iang + 1]
                ang_dict['next'] = ang_dict_next['max']
                nang_dict = self.neighbors_normalized_angles[isite][iang]
                nang_dict_next = self.neighbors_normalized_angles[isite][iang + 1]
                nang_dict['next'] = nang_dict_next['max']
            if self.minimum_angle_factor is not None:
                afact = self.minimum_angle_factor
            else:
                afact = 0.0
            self.neighbors_normalized_angles[isite][-1]['next'] = afact
            self.neighbors_angles[isite][-1]['next'] = afact * self.neighbors_angles[isite][0]['max']

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
        for idp, dp_dict in enumerate(self.neighbors_normalized_distances[ivoronoi]):
            distance_conditions.append([])
            dp = dp_dict['max']
            for ips, (ps, vals) in enumerate(voronoi):
                distance_conditions[idp].append(vals['normalized_distance'] <= dp or
                                                np.isclose(vals['normalized_distance'], dp,
                                                           rtol=0.0, atol=self.normalized_distance_tolerance/2.0))
        return distance_conditions

    def _precompute_angle_conditions(self, ivoronoi, voronoi):
        angle_conditions = []
        for iap, ap_dict in enumerate(self.neighbors_normalized_angles[ivoronoi]):
            angle_conditions.append([])
            ap = ap_dict['max']
            for ips, (ps, vals) in enumerate(voronoi):
                angle_conditions[iap].append(vals['normalized_angle'] >= ap or
                                             np.isclose(vals['normalized_angle'], ap,
                                                        rtol=0.0, atol=self.normalized_angle_tolerance/2.0))
        return angle_conditions

    def neighbors_map(self, isite, distfactor, angfactor, additional_condition):
        if self.neighbors_normalized_distances[isite] is None:
            return None
        dist_where = np.argwhere(np.array([wd['min'] for wd in self.neighbors_normalized_distances[isite]]) <= distfactor)
        if len(dist_where) == 0:
            return None
        idist = dist_where[-1][0]
        ang_where = np.argwhere(np.array([wa['max'] for wa in self.neighbors_normalized_angles[isite]]) >= angfactor)
        if len(ang_where) == 0:
            return None
        iang = ang_where[0][0]
        if self.additional_conditions.count(additional_condition) != 1:
            return None
        i_additional_condition = self.additional_conditions.index(additional_condition)
        return {'i_distfactor': idist, 'i_angfactor': iang, 'i_additional_condition': i_additional_condition}

    def neighbors_surfaces(self, isite, surface_calculation_type=None, max_dist=2.0):
        if self.voronoi_list2[isite] is None:
            return None
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

    def neighbors_surfaces_bounded(self, isite, surface_calculation_options=None):
        if self.voronoi_list2[isite] is None:
            return None
        if surface_calculation_options is None:
            surface_calculation_options = {'type': 'standard_elliptic',
                                           'distance_bounds': {'lower': 1.2, 'upper': 1.8},
                                           'angle_bounds': {'lower': 0.1, 'upper': 0.8}}
        if surface_calculation_options['type'] in ['standard_elliptic', 'standard_diamond', 'standard_spline']:
            plot_type = {'distance_parameter': ('initial_normalized', None),
                         'angle_parameter': ('initial_normalized', None)}
        else:
            raise ValueError('Type "{}" for the surface calculation in DetailedVoronoiContainer '
                             'is invalid'.format(surface_calculation_options['type']))
        max_dist = surface_calculation_options['distance_bounds']['upper'] + 0.1
        bounds_and_limits = self.voronoi_parameters_bounds_and_limits(isite=isite,
                                                                      plot_type=plot_type,
                                                                      max_dist=max_dist)

        distance_bounds = bounds_and_limits['distance_bounds']
        angle_bounds = bounds_and_limits['angle_bounds']
        lower_and_upper_functions = get_lower_and_upper_f(surface_calculation_options=surface_calculation_options)
        mindist = surface_calculation_options['distance_bounds']['lower']
        maxdist = surface_calculation_options['distance_bounds']['upper']
        minang = surface_calculation_options['angle_bounds']['lower']
        maxang = surface_calculation_options['angle_bounds']['upper']

        f_lower = lower_and_upper_functions['lower']
        f_upper = lower_and_upper_functions['upper']
        surfaces = np.zeros((len(distance_bounds), len(angle_bounds)), np.float)
        for idp in range(len(distance_bounds) - 1):
            dp1 = distance_bounds[idp]
            dp2 = distance_bounds[idp+1]
            if dp2 < mindist or dp1 > maxdist:
                continue
            if dp1 < mindist:
                d1 = mindist
            else:
                d1 = dp1
            if dp2 > maxdist:
                d2 = maxdist
            else:
                d2 = dp2
            for iap in range(len(angle_bounds) - 1):
                ap1 = angle_bounds[iap]
                ap2 = angle_bounds[iap+1]
                if ap1 > ap2:
                    ap1 = angle_bounds[iap + 1]
                    ap2 = angle_bounds[iap]
                if ap2 < minang or ap1 > maxang:
                    continue
                intersection, interror = rectangle_surface_intersection(rectangle=((d1, d2),
                                                                                   (ap1, ap2)),
                                                                        f_lower=f_lower,
                                                                        f_upper=f_upper,
                                                                        bounds_lower=[mindist, maxdist],
                                                                        bounds_upper=[mindist, maxdist],
                                                                        check=False)
                surfaces[idp][iap] = intersection
        return surfaces

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
        if self.voronoi_list2[isite] is None:
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

    def maps_and_surfaces_bounded(self, isite, surface_calculation_options=None, additional_conditions=None):
        if self.voronoi_list2[isite] is None:
            return None
        if additional_conditions is None:
            additional_conditions = [self.AC.ONLY_ACB]
        surfaces = self.neighbors_surfaces_bounded(isite=isite, surface_calculation_options=surface_calculation_options)
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

    def neighbors(self, isite, distfactor, angfactor, additional_condition=None):
        idist = None
        dfact = None
        for iwd, wd in enumerate(self.neighbors_normalized_distances[isite]):
            if distfactor >= wd['min']:
                idist = iwd
                dfact = wd['max']
            else:
                break
        iang = None
        afact = None
        for iwa, wa in enumerate(self.neighbors_normalized_angles[isite]):
            if angfactor <= wa['max']:
                iang = iwa
                afact = wa['min']
            else:
                break
        if idist is None or iang is None:
            raise ValueError('Distance or angle parameter not found ...')

        return [nb for nb in self.voronoi_list2[isite] if
                nb['normalized_distance'] <= dfact and nb['normalized_angle'] >= afact]

    def voronoi_parameters_bounds_and_limits(self, isite, plot_type, max_dist):
        #Initializes the distance and angle parameters
        if self.voronoi_list2[isite] is None:
            return None
        if plot_type is None:
            plot_type = {'distance_parameter': ('initial_inverse_opposite', None),
                         'angle_parameter': ('initial_opposite', None)}
        dd = [dist['min'] for dist in self.neighbors_normalized_distances[isite]]
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
            aa.extend([ang['max'] for ang in self.neighbors_normalized_angles[isite]])
            angle_bounds = np.array(aa)
        elif plot_type['angle_parameter'][0] == 'initial_opposite':
            aa = [0.0]
            aa.extend([ang['max'] for ang in self.neighbors_normalized_angles[isite]])
            aa = [1.0 - ang for ang in aa]
            angle_bounds = np.array(aa)
        else:
            raise NotImplementedError('Plotting type "{}" '
                                      'for the angle is not implemented'.format(plot_type['angle_parameter']))
        ang_limits = [0.0, 1.0]
        return {'distance_bounds': distance_bounds, 'distance_limits': dist_limits,
                'angle_bounds': angle_bounds, 'angle_limits': ang_limits}

    def is_close_to(self, other, rtol=0.0, atol=1e-8):
        isclose = (np.isclose(self.normalized_angle_tolerance, other.normalized_angle_tolerance,
                              rtol=rtol, atol=atol) and
                   np.isclose(self.normalized_distance_tolerance, other.normalized_distance_tolerance,
                              rtol=rtol, atol=atol) and
                   self.additional_conditions == other.additional_conditions and
                   self.valences == other.valences)
        if not isclose:
            return isclose
        for isite, site_voronoi in enumerate(self.voronoi_list2):
            self_to_other_nbs = {}
            for inb, nb in enumerate(site_voronoi):
                if nb is None:
                    if other.voronoi_list2[isite] is None:
                        continue
                    else:
                        return False
                else:
                    if other.voronoi_list2[isite] is None:
                        return False
                nb_other = None
                for inb2, nb2 in enumerate(other.voronoi_list2[isite]):
                    if nb['site'] == nb2['site']:
                        self_to_other_nbs[inb] = inb2
                        nb_other = nb2
                        break
                if nb_other is None:
                    return False
                if not np.isclose(nb['distance'], nb_other['distance'],
                                  rtol=rtol, atol=atol):
                    return False
                if not np.isclose(nb['angle'], nb_other['angle'],
                                  rtol=rtol, atol=atol):
                    return False
                if not np.isclose(nb['normalized_distance'], nb_other['normalized_distance'],
                                  rtol=rtol, atol=atol):
                    return False
                if not np.isclose(nb['normalized_angle'], nb_other['normalized_angle'],
                                  rtol=rtol, atol=atol):
                    return False
                if nb['index'] != nb_other['index']:
                    return False
                if nb['site'] != nb_other['site']:
                    return False
        return True

    def __eq__(self, other):
        return (self.normalized_angle_tolerance == other.normalized_angle_tolerance and
                self.normalized_distance_tolerance == other.normalized_distance_tolerance and
                self.additional_conditions == other.additional_conditions and
                self.valences == other.valences and
                self.voronoi_list2 == other.voronoi_list2 and
                self.structure == other.structure)

    def __ne__(self, other):
        return not self == other

    def to_bson_voronoi_list2(self):
        """
        Transforms the voronoi_list into a vlist + bson_nb_voro_list, that are BSON-encodable.
        :return: [vlist, bson_nb_voro_list], to be used in the as_dict method
        """
        bson_nb_voro_list2 = [None] * len(self.voronoi_list2)
        for ivoro, voro in enumerate(self.voronoi_list2):
            if voro is None or voro == 'None':
                continue
            site_voro = []
            # {'site': neighbors[nn[1]],
            #  'angle': sa,
            #  'distance': distances[nn[1]],
            #  'index': myindex}
            for nb_dict in voro:
                site = nb_dict['site']
                site_dict = {key: val for key, val in nb_dict.items() if key not in ['site']}
                #site_voro.append([ps.as_dict(), dd]) [float(c) for c in self._fcoords]
                diff = site._fcoords - self.structure[nb_dict['index']]._fcoords
                site_voro.append([[nb_dict['index'], [float(c) for c in diff]],
                                  site_dict])
            bson_nb_voro_list2[ivoro] = site_voro
        return bson_nb_voro_list2

    def as_dict(self):
        """
        Bson-serializable dict representation of the VoronoiContainer.
        :return: dictionary that is BSON-encodable
        """
        bson_nb_voro_list2 = self.to_bson_voronoi_list2()
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "bson_nb_voro_list2": bson_nb_voro_list2,
                # "neighbors_lists": self.neighbors_lists,
                "structure": self.structure.as_dict(),
                "normalized_angle_tolerance": self.normalized_angle_tolerance,
                "normalized_distance_tolerance": self.normalized_distance_tolerance,
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
        voronoi_list2 = from_bson_voronoi_list2(d['bson_nb_voro_list2'], structure)
        maximum_distance_factor = d['maximum_distance_factor'] if 'maximum_distance_factor' in d else None
        minimum_angle_factor = d['minimum_angle_factor'] if 'minimum_angle_factor' in d else None
        return cls(structure=structure, voronoi_list2=voronoi_list2,
                   # neighbors_lists=neighbors_lists,
                   normalized_angle_tolerance=d['normalized_angle_tolerance'],
                   normalized_distance_tolerance=d['normalized_distance_tolerance'],
                   additional_conditions=d['additional_conditions'],
                   valences=d['valences'],
                   maximum_distance_factor=maximum_distance_factor,
                   minimum_angle_factor=minimum_angle_factor)