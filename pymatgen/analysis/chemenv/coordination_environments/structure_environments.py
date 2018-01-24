# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module contains objects that are used to describe the environments in a structure. The most detailed object
(StructureEnvironments) contains a very thorough analysis of the environments of a given atom but is difficult to
used as such. The LightStructureEnvironments object is a lighter version that is obtained by applying a "strategy"
on the StructureEnvironments object. Basically, the LightStructureEnvironments provides the coordination environment(s)
and possibly some fraction corresponding to these.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


import numpy as np
from collections import OrderedDict
from pymatgen.core.sites import PeriodicSite
from monty.json import MSONable, MontyDecoder
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Structure
from monty.json import jsanitize
from pymatgen.analysis.chemenv.coordination_environments.voronoi import DetailedVoronoiContainer
from pymatgen.analysis.chemenv.utils.chemenv_errors import ChemenvError
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.utils.defs_utils import AdditionalConditions

allcg = AllCoordinationGeometries()
symbol_cn_mapping = allcg.get_symbol_cn_mapping()

class StructureEnvironments(MSONable):
    """
    Class used to store the chemical environments of a given structure.
    """
    AC = AdditionalConditions()

    class NeighborsSet():
        """
        Class used to store a given set of neighbors of a given site (based on the detailed_voronoi).
        """
        def __init__(self, structure, isite, detailed_voronoi, site_voronoi_indices, sources=None):
            self.structure = structure
            self.isite = isite
            self.detailed_voronoi = detailed_voronoi
            self.voronoi = detailed_voronoi.voronoi_list2[isite]
            myset = set(site_voronoi_indices)
            if len(myset) != len(site_voronoi_indices):
                raise ValueError('Set of neighbors contains duplicates !')
            self.site_voronoi_indices = sorted(myset)
            if sources is None:
                self.sources = [{'origin': 'UNKNOWN'}]
            elif isinstance(sources, list):
                self.sources = sources
            else:
                self.sources = [sources]

        def get_neighb_voronoi_indices(self, permutation):
            return [self.site_voronoi_indices[ii] for ii in permutation]

        @property
        def neighb_coords(self):
            return [self.voronoi[inb]['site'].coords for inb in self.site_voronoi_indices]

        @property
        def neighb_sites(self):
            return [self.voronoi[inb]['site'] for inb in self.site_voronoi_indices]

        @property
        def neighb_sites_and_indices(self):
            return [{'site': self.voronoi[inb]['site'],
                     'index': self.voronoi[inb]['index']} for inb in self.site_voronoi_indices]

        @property
        def coords(self):
            coords = [self.structure[self.isite].coords]
            coords.extend(self.neighb_coords)
            return coords

        @property
        def normalized_distances(self):
            return [self.voronoi[inb]['normalized_distance'] for inb in self.site_voronoi_indices]

        @property
        def normalized_angles(self):
            return [self.voronoi[inb]['normalized_angle'] for inb in self.site_voronoi_indices]

        @property
        def distances(self):
            return [self.voronoi[inb]['distance'] for inb in self.site_voronoi_indices]

        @property
        def angles(self):
            return [self.voronoi[inb]['angle'] for inb in self.site_voronoi_indices]

        @property
        def sphere_fraction_angles(self):
            return [0.25*self.voronoi[inb]['angle']/np.pi for inb in self.site_voronoi_indices]

        @property
        def info(self):
            was = self.normalized_angles
            wds = self.normalized_distances
            angles = self.angles
            distances = self.distances
            return {'normalized_angles': was,
                    'normalized_distances': wds,
                    'normalized_angles_sum': np.sum(was),
                    'normalized_angles_mean': np.mean(was),
                    'normalized_angles_std': np.std(was),
                    'normalized_angles_min': np.min(was),
                    'normalized_angles_max': np.max(was),
                    'normalized_distances_mean': np.mean(wds),
                    'normalized_distances_std': np.std(wds),
                    'normalized_distances_min': np.min(wds),
                    'normalized_distances_max': np.max(wds),
                    'angles': angles,
                    'distances': distances,
                    'angles_sum': np.sum(angles),
                    'angles_mean': np.mean(angles),
                    'angles_std': np.std(angles),
                    'angles_min': np.min(angles),
                    'angles_max': np.max(angles),
                    'distances_mean': np.mean(distances),
                    'distances_std': np.std(distances),
                    'distances_min': np.min(distances),
                    'distances_max': np.max(distances)
                    }

        def distance_plateau(self):
            all_nbs_normalized_distances_sorted = sorted([nb['normalized_distance'] for nb in self.voronoi],
                                                         reverse=True)
            maxdist = np.max(self.normalized_distances)
            plateau = None
            for idist, dist in enumerate(all_nbs_normalized_distances_sorted):
                if np.isclose(dist, maxdist,
                              rtol=0.0, atol=self.detailed_voronoi.normalized_distance_tolerance):
                    if idist == 0:
                        plateau = np.inf
                    else:
                        plateau = all_nbs_normalized_distances_sorted[idist-1] - maxdist
                    break
            if plateau is None:
                raise ValueError('Plateau not found ...')
            return plateau

        def angle_plateau(self):
            all_nbs_normalized_angles_sorted = sorted([nb['normalized_angle'] for nb in self.voronoi])
            minang = np.min(self.normalized_angles)
            print('minang', minang)
            print('all_nbs_normalized_angles_sorted', all_nbs_normalized_angles_sorted)
            for nb in self.voronoi:
                print(nb)
            plateau = None
            for iang, ang in enumerate(all_nbs_normalized_angles_sorted):
                if np.isclose(ang, minang,
                              rtol=0.0, atol=self.detailed_voronoi.normalized_angle_tolerance):
                    if iang == 0:
                        plateau = minang
                    else:
                        plateau = minang - all_nbs_normalized_angles_sorted[iang-1]
                    break
            if plateau is None:
                raise ValueError('Plateau not found ...')
            return plateau

        def voronoi_grid_surface_points(self, additional_condition=1, other_origins='DO_NOTHING'):
            """
            Get the surface points in the Voronoi grid for this neighbor from the sources.
            The general shape of the points should look like a staircase such as in the following figure :

               ^
            0.0|
               |
               |      B----C
               |      |    |
               |      |    |
            a  |      k    D-------E
            n  |      |            |
            g  |      |            |
            l  |      |            |
            e  |      j            F----n---------G
               |      |                           |
               |      |                           |
               |      A----g-------h----i---------H
               |
               |
            1.0+------------------------------------------------->
              1.0              distance              2.0   ->+Inf

            :param additional_condition: Additional condition for the neighbors.
            :param other_origins: What to do with sources that do not come from the Voronoi grid (e.g. "from hints")
            """
            mysrc = []
            for src in self.sources:
                if src['origin'] == 'dist_ang_ac_voronoi':
                    if src['ac'] != additional_condition:
                        continue
                    mysrc.append(src)
                else:
                    if other_origins == 'DO_NOTHING':
                        continue
                    else:
                        raise NotImplementedError('Nothing implemented for other sources ...')
            if len(mysrc) == 0:
                return None

            dists = [src['dp_dict']['min'] for src in mysrc]
            angs = [src['ap_dict']['max'] for src in mysrc]
            next_dists = [src['dp_dict']['next'] for src in mysrc]
            next_angs = [src['ap_dict']['next'] for src in mysrc]

            points_dict = OrderedDict()

            pdists = []
            pangs = []

            for isrc in range(len(mysrc)):
                if not any(np.isclose(pdists, dists[isrc])):
                    pdists.append(dists[isrc])
                if not any(np.isclose(pdists, next_dists[isrc])):
                    pdists.append(next_dists[isrc])
                if not any(np.isclose(pangs, angs[isrc])):
                    pangs.append(angs[isrc])
                if not any(np.isclose(pangs, next_angs[isrc])):
                    pangs.append(next_angs[isrc])
                d1_indices = np.argwhere(np.isclose(pdists, dists[isrc])).flatten()
                if len(d1_indices) != 1:
                    raise ValueError('Distance parameter not found ...')
                d2_indices = np.argwhere(np.isclose(pdists, next_dists[isrc])).flatten()
                if len(d2_indices) != 1:
                    raise ValueError('Distance parameter not found ...')
                a1_indices = np.argwhere(np.isclose(pangs, angs[isrc])).flatten()
                if len(a1_indices) != 1:
                    raise ValueError('Angle parameter not found ...')
                a2_indices = np.argwhere(np.isclose(pangs, next_angs[isrc])).flatten()
                if len(a2_indices) != 1:
                    raise ValueError('Angle parameter not found ...')
                id1 = d1_indices[0]
                id2 = d2_indices[0]
                ia1 = a1_indices[0]
                ia2 = a2_indices[0]
                for id_ia in [(id1, ia1), (id1, ia2), (id2, ia1), (id2, ia2)]:
                    if id_ia not in points_dict:
                        points_dict[id_ia] = 0
                    points_dict[id_ia] += 1

            new_pts = []
            for pt, pt_nb in points_dict.items():
                if pt_nb % 2 == 1:
                    new_pts.append(pt)

            sorted_points = [(0, 0)]
            move_ap_index = True
            while True:
                last_pt = sorted_points[-1]
                if move_ap_index: # "Move" the angle parameter
                    idp = last_pt[0]
                    iap = None
                    for pt in new_pts:
                        if pt[0] == idp and pt != last_pt:
                            iap = pt[1]
                            break
                else:             # "Move" the distance parameter
                    idp = None
                    iap = last_pt[1]
                    for pt in new_pts:
                        if pt[1] == iap and pt != last_pt:
                            idp = pt[0]
                            break
                if (idp, iap) == (0, 0):
                    break
                if (idp, iap) in sorted_points:
                    raise ValueError('Error sorting points ...')
                sorted_points.append((idp, iap))
                move_ap_index = not move_ap_index

            points = [(pdists[idp], pangs[iap]) for (idp, iap) in sorted_points]
            return points

        @property
        def source(self):
            if len(self.sources) != 1:
                raise RuntimeError('Number of sources different from 1 !')
            return self.sources[0]

        def add_source(self, source):
            if source not in self.sources:
                self.sources.append(source)

        def __len__(self):
            return len(self.site_voronoi_indices)

        def __hash__(self):
            return len(self.site_voronoi_indices)

        def __eq__(self, other):
            return self.isite == other.isite and self.site_voronoi_indices == other.site_voronoi_indices

        def __ne__(self, other):
            return not self == other

        def __str__(self):
            out = 'Neighbors Set for site #{:d} :\n'.format(self.isite)
            out += ' - Coordination number : {:d}\n'.format(len(self))
            out += ' - Voronoi indices : {}\n'.format(', '.join(['{:d}'.format(site_voronoi_index)
                                                                 for site_voronoi_index in self.site_voronoi_indices]))
            return out

        def as_dict(self):
            return {'isite': self.isite,
                    'site_voronoi_indices': self.site_voronoi_indices,
                    'sources': self.sources}

        @classmethod
        def from_dict(cls, dd, structure, detailed_voronoi):
            return cls(structure=structure,
                       isite=dd['isite'],
                       detailed_voronoi=detailed_voronoi,
                       site_voronoi_indices=dd['site_voronoi_indices'],
                       sources=dd['sources'])

    def __init__(self, voronoi, valences, sites_map, equivalent_sites,
                 ce_list, structure, neighbors_sets=None, info=None):
        """
        Constructor for the StructureEnvironments object.
        :param voronoi: VoronoiContainer object for the structure
        :param valences: Valences provided
        :param sites_map: Mapping of equivalent sites to the unequivalent sites that have been computed.
        :param equivalent_sites: List of list of equivalent sites of the structure
        :param struct_sites_to_irreducible_site_list_map: Maps the index of a site to the index of the item in the
        list of equivalent sites to which the site belongs.
        :param ce_list: List of chemical environments
        :param structure: Structure object
        """
        self.voronoi = voronoi
        self.valences = valences
        self.sites_map = sites_map
        self.equivalent_sites = equivalent_sites
        #self.struct_sites_to_irreducible_site_list_map = struct_sites_to_irreducible_site_list_map
        self.ce_list = ce_list
        self.structure = structure
        if neighbors_sets is None:
            self.neighbors_sets = [None] * len(self.structure)
        else:
            self.neighbors_sets = neighbors_sets
        self.info = info

    def init_neighbors_sets(self, isite, additional_conditions=None, valences=None):
        site_voronoi = self.voronoi.voronoi_list2[isite]
        if site_voronoi is None:
            return
        if additional_conditions is None:
            additional_conditions = self.AC.ALL
        if (self.AC.ONLY_ACB in additional_conditions or
                self.AC.ONLY_ACB_AND_NO_E2SEB) and valences is None:
            raise ChemenvError('StructureEnvironments', 'init_neighbors_sets',
                               'Valences are not given while only_anion_cation_bonds are allowed. Cannot continue')
        site_distance_parameters = self.voronoi.neighbors_normalized_distances[isite]
        site_angle_parameters = self.voronoi.neighbors_normalized_angles[isite]
        # Precompute distance conditions
        distance_conditions = []
        for idp, dp_dict in enumerate(site_distance_parameters):
            distance_conditions.append([])
            for inb, voro_nb_dict in enumerate(site_voronoi):
                cond = inb in dp_dict['nb_indices']
                distance_conditions[idp].append(cond)
        # Precompute angle conditions
        angle_conditions = []
        for iap, ap_dict in enumerate(site_angle_parameters):
            angle_conditions.append([])
            for inb, voro_nb_dict in enumerate(site_voronoi):
                cond = inb in ap_dict['nb_indices']
                angle_conditions[iap].append(cond)
        # Precompute additional conditions
        precomputed_additional_conditions = {ac: [] for ac in additional_conditions}
        for inb, voro_nb_dict in enumerate(site_voronoi):
            for ac in additional_conditions:
                cond = self.AC.check_condition(condition=ac, structure=self.structure,
                                               parameters={'valences': valences,
                                                           'neighbor_index': voro_nb_dict['index'],
                                                           'site_index': isite})
                precomputed_additional_conditions[ac].append(cond)
        # Add the neighbors sets based on the distance/angle/additional parameters
        for idp, dp_dict in enumerate(site_distance_parameters):
            for iap, ap_dict in enumerate(site_angle_parameters):
                for iac, ac in enumerate(additional_conditions):
                    src = {'origin': 'dist_ang_ac_voronoi',
                           'idp': idp, 'iap': iap, 'dp_dict': dp_dict, 'ap_dict': ap_dict,
                           'iac': iac, 'ac': ac, 'ac_name': self.AC.CONDITION_DESCRIPTION[ac]}
                    site_voronoi_indices = [inb for inb, voro_nb_dict in enumerate(site_voronoi)
                                            if (distance_conditions[idp][inb] and
                                                angle_conditions[iap][inb] and
                                                precomputed_additional_conditions[ac][inb])]
                    nb_set = self.NeighborsSet(structure=self.structure,
                                               isite=isite,
                                               detailed_voronoi=self.voronoi,
                                               site_voronoi_indices=site_voronoi_indices,
                                               sources=src)
                    self.add_neighbors_set(isite=isite, nb_set=nb_set)

    def add_neighbors_set(self, isite, nb_set):
        if self.neighbors_sets[isite] is None:
            self.neighbors_sets[isite] = {}
            self.ce_list[isite] = {}
        cn = len(nb_set)
        if cn not in self.neighbors_sets[isite]:
            self.neighbors_sets[isite][cn] = []
            self.ce_list[isite][cn] = []
        try:
            nb_set_index = self.neighbors_sets[isite][cn].index(nb_set)
            self.neighbors_sets[isite][cn][nb_set_index].add_source(nb_set.source)
        except ValueError:
            self.neighbors_sets[isite][cn].append(nb_set)
            self.ce_list[isite][cn].append(None)

    def update_coordination_environments(self, isite, cn, nb_set, ce):
        if self.ce_list[isite] is None:
            self.ce_list[isite] = {}
        if cn not in self.ce_list[isite]:
            self.ce_list[isite][cn] = []
        try:
            nb_set_index = self.neighbors_sets[isite][cn].index(nb_set)
        except ValueError:
            raise ValueError('Neighbors set not found in the structure environments')
        if nb_set_index == len(self.ce_list[isite][cn]):
            self.ce_list[isite][cn].append(ce)
        elif nb_set_index < len(self.ce_list[isite][cn]):
            self.ce_list[isite][cn][nb_set_index] = ce
        else:
            raise ValueError('Neighbors set not yet in ce_list !')

    def update_site_info(self, isite, info_dict):
        if 'sites_info' not in self.info:
            self.info['sites_info'] = [{} for _ in range(len(self.structure))]
        self.info['sites_info'][isite].update(info_dict)

    def get_coordination_environments(self, isite, cn, nb_set):
        if self.ce_list[isite] is None:
            return None
        if cn not in self.ce_list[isite]:
            return None
        try:
            nb_set_index = self.neighbors_sets[isite][cn].index(nb_set)
        except ValueError:
            return None
        return self.ce_list[isite][cn][nb_set_index]

    def get_csm(self, isite, mp_symbol):
        csms = self.get_csms(isite, mp_symbol)
        if len(csms) != 1:
            raise ChemenvError('StructureEnvironments',
                               'get_csm',
                               'Number of csms for site #{} with '
                               'mp_symbol "{}" = {}'.format(str(isite),
                                                            mp_symbol,
                                                            str(len(csms))))
        return csms[0]

    def get_csms(self, isite, mp_symbol):
        """
        Returns the continuous symmetry measure(s) of site with index isite with respect to the
         perfect coordination environment with mp_symbol. For some environments, a given mp_symbol might not
         be available (if there is no voronoi parameters leading to a number of neighbours corresponding to
         the coordination number of environment mp_symbol). For some environments, a given mp_symbol might
         lead to more than one csm (when two or more different voronoi parameters lead to different neighbours
         but with same number of neighbours).
        :param isite: Index of the site
        :param mp_symbol: MP symbol of the perfect environment for which the csm has to be given
        :return: List of csms for site isite with respect to geometry mp_symbol
        """
        cn = symbol_cn_mapping[mp_symbol]
        if cn not in self.ce_list[isite]:
            return []
        else:
            return [envs[mp_symbol] for envs in self.ce_list[isite][cn]]

    def plot_csm_and_maps(self, isite, max_csm=8.0):
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
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print('Plotting Chemical Environments requires matplotlib ... exiting "plot" function')
            return
        fig = self.get_csm_and_maps(isite=isite, max_csm=max_csm)
        if fig is None:
            return
        plt.show()

    def get_csm_and_maps(self, isite, max_csm=8.0, figsize=None, symmetry_measure_type=None):
        try:
            import matplotlib.pyplot as plt
            from matplotlib.gridspec import GridSpec
        except ImportError:
            print('Plotting Chemical Environments requires matplotlib ... exiting "plot" function')
            return

        if symmetry_measure_type is None:
            symmetry_measure_type = 'csm_wcs_ctwcc'
        # Initializes the figure
        if figsize is None:
            fig = plt.figure()
        else:
            fig = plt.figure(figsize=figsize)
        gs = GridSpec(2, 1, hspace=0.0, wspace=0.0)
        subplot = fig.add_subplot(gs[:])
        subplot_distang = subplot.twinx()

        ix = 0
        cn_maps = []
        all_wds = []
        all_was = []
        max_wd = 0.0
        for cn, nb_sets in self.neighbors_sets[isite].items():
            for inb_set, nb_set in enumerate(nb_sets):
                ce = self.ce_list[isite][cn][inb_set]
                if ce is None:
                    continue
                mingeoms = ce.minimum_geometries(max_csm=max_csm)
                if len(mingeoms) == 0:
                    continue
                wds = nb_set.normalized_distances
                max_wd = max(max_wd, max(wds))
                all_wds.append(wds)
                all_was.append(nb_set.normalized_angles)
                for mp_symbol, cg_dict in mingeoms:
                    csm = cg_dict['other_symmetry_measures'][symmetry_measure_type]
                    subplot.plot(ix, csm, 'ob')
                    subplot.annotate(mp_symbol, xy = (ix, csm))
                cn_maps.append((cn, inb_set))
                ix += 1

        if max_wd < 1.225:
            ymax_wd = 1.25
            yticks_wd = np.linspace(1.0, ymax_wd, 6)
        elif max_wd < 1.36:
            ymax_wd = 1.4
            yticks_wd = np.linspace(1.0, ymax_wd, 5)
        elif max_wd < 1.45:
            ymax_wd = 1.5
            yticks_wd = np.linspace(1.0, ymax_wd, 6)
        elif max_wd < 1.55:
            ymax_wd = 1.6
            yticks_wd = np.linspace(1.0, ymax_wd, 7)
        elif max_wd < 1.75:
            ymax_wd = 1.8
            yticks_wd = np.linspace(1.0, ymax_wd, 5)
        elif max_wd < 1.95:
            ymax_wd = 2.0
            yticks_wd = np.linspace(1.0, ymax_wd, 6)
        elif max_wd < 2.35:
            ymax_wd = 2.5
            yticks_wd = np.linspace(1.0, ymax_wd, 7)
        else:
            ymax_wd = np.ceil(1.1*max_wd)
            yticks_wd = np.linspace(1.0, ymax_wd, 6)

        yticks_wa = np.linspace(0.0, 1.0, 6)

        frac_bottom = 0.05
        frac_top = 0.05
        frac_middle = 0.1
        yamin = frac_bottom
        yamax = 0.5 - frac_middle / 2
        ydmin = 0.5 + frac_middle / 2
        ydmax = 1.0 - frac_top

        def yang(wa):
            return (yamax-yamin) * np.array(wa) + yamin

        def ydist(wd):
            return (np.array(wd) - 1.0) / (ymax_wd - 1.0) * (ydmax - ydmin) + ydmin

        for ix, was in enumerate(all_was):
            subplot_distang.plot(0.2+ix*np.ones_like(was), yang(was), '<g')
            if np.mod(ix, 2) == 0:
                alpha = 0.3
            else:
                alpha = 0.1
            subplot_distang.fill_between([-0.5+ix, 0.5+ix],
                                         [1.0, 1.0], 0.0,
                                         facecolor='k', alpha=alpha, zorder=-1000)
        for ix, wds in enumerate(all_wds):
            subplot_distang.plot(0.2+ix*np.ones_like(wds), ydist(wds), 'sm')



        subplot_distang.plot([-0.5, len(cn_maps)], [0.5, 0.5], 'k--', alpha=0.5)

        yticks = yang(yticks_wa).tolist()
        yticks.extend(ydist(yticks_wd).tolist())
        yticklabels = yticks_wa.tolist()
        yticklabels.extend(yticks_wd.tolist())
        subplot_distang.set_yticks(yticks)
        subplot_distang.set_yticklabels(yticklabels)

        fake_subplot_ang = fig.add_subplot(gs[1], frame_on=False)
        fake_subplot_dist = fig.add_subplot(gs[0], frame_on=False)
        fake_subplot_ang.set_yticks([])
        fake_subplot_dist.set_yticks([])
        fake_subplot_ang.set_xticks([])
        fake_subplot_dist.set_xticks([])
        fake_subplot_ang.set_ylabel('Angle parameter', labelpad=45, rotation=-90)
        fake_subplot_dist.set_ylabel('Distance parameter', labelpad=45, rotation=-90)
        fake_subplot_ang.yaxis.set_label_position("right")
        fake_subplot_dist.yaxis.set_label_position("right")


        subplot_distang.set_ylim([0.0, 1.0])
        subplot.set_xticks(range(len(cn_maps)))
        subplot.set_ylabel('Continuous symmetry measure')
        subplot.set_xlim([-0.5, len(cn_maps)-0.5])
        subplot_distang.set_xlim([-0.5, len(cn_maps)-0.5])
        subplot.set_xticklabels([str(cn_map) for cn_map in cn_maps])

        return fig, subplot

    def get_environments_figure(self, isite, plot_type=None, title='Coordination numbers', max_dist=2.0,
                                additional_condition=AC.ONLY_ACB, colormap=None, figsize=None,
                                strategy=None):
        """
        Plotting of the coordination environments of a given site for all the distfactor/angfactor regions. The
        chemical environments with the lowest continuous symmetry measure is shown for each distfactor/angfactor
        region as the value for the color of that distfactor/angfactor region (using a colormap).
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
        if plot_type is None:
            plot_type = {'distance_parameter': ('initial_normalized', None),
                         'angle_parameter': ('initial_normalized_inverted', None)}
        if colormap is None:
            mycm = cm.jet
        else:
            mycm = colormap
        mymin = 0.0
        mymax = 10.0
        norm = Normalize(vmin=mymin, vmax=mymax)
        scalarmap = cm.ScalarMappable(norm=norm, cmap=mycm)
        dist_limits = [1.0, max_dist]
        ang_limits = [0.0, 1.0]
        if plot_type['distance_parameter'][0] == 'one_minus_inverse_alpha_power_n':
            if plot_type['distance_parameter'][1] is None:
                exponent = 3
            else:
                exponent = plot_type['distance_parameter'][1]['exponent']
            xlabel = 'Distance parameter : $1.0-\\frac{{1.0}}{{\\alpha^{{{:d}}}}}$'.format(exponent)
            def dp_func(dp):
                return 1.0-1.0/np.power(dp, exponent)
        elif plot_type['distance_parameter'][0] == 'initial_normalized':
            xlabel = 'Distance parameter : $\\alpha$'
            def dp_func(dp):
                return dp
        else:
            raise ValueError('Wrong value for distance parameter plot type "{}"'.
                             format(plot_type['distance_parameter'][0]))

        if plot_type['angle_parameter'][0] == 'one_minus_gamma':
            ylabel = 'Angle parameter : $1.0-\\gamma$'
            def ap_func(ap):
                return 1.0-ap
        elif plot_type['angle_parameter'][0] in ['initial_normalized_inverted', 'initial_normalized']:
            ylabel = 'Angle parameter : $\\gamma$'
            def ap_func(ap):
                return ap
        else:
            raise ValueError('Wrong value for angle parameter plot type "{}"'.
                             format(plot_type['angle_parameter'][0]))
        dist_limits = [dp_func(dp) for dp in dist_limits]
        ang_limits = [ap_func(ap) for ap in ang_limits]

        for cn, cn_nb_sets in self.neighbors_sets[isite].items():
            for inb_set, nb_set in enumerate(cn_nb_sets):
                nb_set_surface_pts = nb_set.voronoi_grid_surface_points()
                if nb_set_surface_pts is None:
                    continue
                ce = self.ce_list[isite][cn][inb_set]
                if ce is None:
                    mycolor = 'w'
                    myinvcolor = 'k'
                    mytext = '{:d}'.format(cn)
                else:
                    mingeom = ce.minimum_geometry()
                    if mingeom is not None:
                        mp_symbol = mingeom[0]
                        csm = mingeom[1]['symmetry_measure']
                        mycolor = scalarmap.to_rgba(csm)
                        myinvcolor = [1.0 - mycolor[0], 1.0 - mycolor[1], 1.0 - mycolor[2], 1.0]
                        mytext = '{}'.format(mp_symbol)
                    else:
                        mycolor = 'w'
                        myinvcolor = 'k'
                        mytext = '{:d}'.format(cn)
                nb_set_surface_pts = [(dp_func(pt[0]), ap_func(pt[1])) for pt in nb_set_surface_pts]
                polygon = Polygon(nb_set_surface_pts, closed=True, edgecolor='k', facecolor=mycolor, linewidth=1.2)
                subplot.add_patch(polygon)
                myipt = len(nb_set_surface_pts) / 2
                ipt = int(myipt)
                if myipt != ipt:
                    raise RuntimeError('Number of surface points not even')
                patch_center = ((nb_set_surface_pts[0][0] + min(nb_set_surface_pts[ipt][0], dist_limits[1])) / 2,
                                (nb_set_surface_pts[0][1] + nb_set_surface_pts[ipt][1]) / 2)

                if (np.abs(nb_set_surface_pts[-1][1] - nb_set_surface_pts[-2][1]) > 0.06 and
                      np.abs(min(nb_set_surface_pts[-1][0], dist_limits[1]) - nb_set_surface_pts[0][0]) > 0.125):
                    xytext = ((min(nb_set_surface_pts[-1][0], dist_limits[1]) + nb_set_surface_pts[0][0]) / 2,
                              (nb_set_surface_pts[-1][1] + nb_set_surface_pts[-2][1]) / 2)
                    subplot.annotate(mytext, xy=xytext,
                                     ha='center', va='center', color=myinvcolor, fontsize='x-small')
                elif (np.abs(nb_set_surface_pts[ipt][1] - nb_set_surface_pts[0][1]) > 0.1 and
                    np.abs(min(nb_set_surface_pts[ipt][0], dist_limits[1]) - nb_set_surface_pts[0][0]) > 0.125):
                    xytext = patch_center
                    subplot.annotate(mytext, xy=xytext,
                                     ha='center', va='center', color=myinvcolor, fontsize='x-small')

        subplot.set_title(title)
        subplot.set_xlabel(xlabel)
        subplot.set_ylabel(ylabel)

        dist_limits.sort()
        ang_limits.sort()
        subplot.set_xlim(dist_limits)
        subplot.set_ylim(ang_limits)
        if strategy is not None:
            try:
                strategy.add_strategy_visualization_to_subplot(subplot=subplot)
            except:
                pass
        if plot_type['angle_parameter'][0] == 'initial_normalized_inverted':
            subplot.axes.invert_yaxis()

        scalarmap.set_array([mymin, mymax])
        cb = fig.colorbar(scalarmap, ax=subplot, extend='max')
        cb.set_label('Continuous symmetry measure')
        return fig, subplot

    def plot_environments(self,  isite, plot_type=None, title='Coordination numbers', max_dist=2.0,
                          additional_condition=AC.ONLY_ACB, figsize=None, strategy=None):
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
        fig, subplot = self.get_environments_figure(isite=isite, plot_type=plot_type, title=title, max_dist=max_dist,
                                                    additional_condition=additional_condition, figsize=figsize,
                                                    strategy=strategy)
        if fig is None:
            return
        fig.show()

    def save_environments_figure(self,  isite, imagename='image.png', plot_type=None, title='Coordination numbers',
                                 max_dist=2.0, additional_condition=AC.ONLY_ACB, figsize=None):
        fig, subplot = self.get_environments_figure(isite=isite, plot_type=plot_type, title=title, max_dist=max_dist,
                                           additional_condition=additional_condition, figsize=figsize)
        if fig is None:
            return
        fig.savefig(imagename)

    def differences_wrt(self, other):
        differences = []
        if self.structure != other.structure:
            differences.append({'difference': 'structure',
                                'comparison': '__eq__',
                                'self': self.structure,
                                'other': other.structure})
            differences.append({'difference': 'PREVIOUS DIFFERENCE IS DISMISSIVE',
                                'comparison': 'differences_wrt'})
            return differences
        if self.valences != other.valences:
            differences.append({'difference': 'valences',
                                'comparison': '__eq__',
                                'self': self.valences,
                                'other': other.valences})
        if self.info != other.info:
            differences.append({'difference': 'info',
                                'comparison': '__eq__',
                                'self': self.info,
                                'other': other.info})
        if self.voronoi != other.voronoi:
            if self.voronoi.is_close_to(other.voronoi):
                differences.append({'difference': 'voronoi',
                                    'comparison': '__eq__',
                                    'self': self.voronoi,
                                    'other': other.voronoi})
                differences.append({'difference': 'PREVIOUS DIFFERENCE IS DISMISSIVE',
                                    'comparison': 'differences_wrt'})
                return differences
            else:
                differences.append({'difference': 'voronoi',
                                    'comparison': 'is_close_to',
                                    'self': self.voronoi,
                                    'other': other.voronoi})
                # TODO: make it possible to have "close" voronoi's
                differences.append({'difference': 'PREVIOUS DIFFERENCE IS DISMISSIVE',
                                    'comparison': 'differences_wrt'})
                return differences
        for isite, self_site_nb_sets in enumerate(self.neighbors_sets):
            other_site_nb_sets = other.neighbors_sets[isite]
            if self_site_nb_sets is None:
                if other_site_nb_sets is None:
                    continue
                else:
                    differences.append({'difference': 'neighbors_sets[isite={:d}]'.format(isite),
                                        'comparison': 'has_neighbors',
                                        'self': 'None',
                                        'other': set(other_site_nb_sets.keys())})
                    continue
            elif other_site_nb_sets is None:
                differences.append({'difference': 'neighbors_sets[isite={:d}]'.format(isite),
                                    'comparison': 'has_neighbors',
                                    'self': set(self_site_nb_sets.keys()),
                                    'other': 'None'})
                continue
            self_site_cns = set(self_site_nb_sets.keys())
            other_site_cns = set(other_site_nb_sets.keys())
            if self_site_cns != other_site_cns:
                differences.append({'difference': 'neighbors_sets[isite={:d}]'.format(isite),
                                    'comparison': 'coordination_numbers',
                                    'self': self_site_cns,
                                    'other': other_site_cns})
            common_cns = self_site_cns.intersection(other_site_cns)
            for cn in common_cns:
                other_site_cn_nb_sets = other_site_nb_sets[cn]
                self_site_cn_nb_sets = self_site_nb_sets[cn]
                set_self_site_cn_nb_sets = set(self_site_cn_nb_sets)
                set_other_site_cn_nb_sets = set(other_site_cn_nb_sets)
                if set_self_site_cn_nb_sets != set_other_site_cn_nb_sets:
                    differences.append({'difference': 'neighbors_sets[isite={:d}][cn={:d}]'.format(isite, cn),
                                        'comparison': 'neighbors_sets',
                                        'self': self_site_cn_nb_sets,
                                        'other': other_site_cn_nb_sets})
                common_nb_sets = set_self_site_cn_nb_sets.intersection(set_other_site_cn_nb_sets)
                for nb_set in common_nb_sets:
                    inb_set_self = self_site_cn_nb_sets.index(nb_set)
                    inb_set_other = other_site_cn_nb_sets.index(nb_set)
                    self_ce = self.ce_list[isite][cn][inb_set_self]
                    other_ce = other.ce_list[isite][cn][inb_set_other]
                    if self_ce != other_ce:
                        if self_ce.is_close_to(other_ce):
                            differences.append({'difference': 'ce_list[isite={:d}][cn={:d}]'
                                                              '[inb_set={:d}]'.format(isite, cn, inb_set_self),
                                                'comparison': '__eq__',
                                                'self': self_ce,
                                                'other': other_ce})
                        else:
                            differences.append({'difference': 'ce_list[isite={:d}][cn={:d}]'
                                                              '[inb_set={:d}]'.format(isite, cn, inb_set_self),
                                                'comparison': 'is_close_to',
                                                'self': self_ce,
                                                'other': other_ce})
        return differences

    def __eq__(self, other):
        if len(self.ce_list) != len(other.ce_list):
            return False
        if self.voronoi != other.voronoi:
            return False
        if len(self.valences) != len(other.valences):
            return False
        if self.sites_map != other.sites_map:
            return False
        if self.equivalent_sites != other.equivalent_sites:
            return False
        if self.structure != other.structure:
            return False
        if self.info != other.info:
            return False
        for isite, site_ces in enumerate(self.ce_list):
            site_nb_sets_self = self.neighbors_sets[isite]
            site_nb_sets_other = other.neighbors_sets[isite]
            if site_nb_sets_self != site_nb_sets_other:
                return False
            if site_ces != other.ce_list[isite]:
                return False
        return True

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        """
        Bson-serializable dict representation of the StructureEnvironments object.
        :return: Bson-serializable dict representation of the StructureEnvironments object.
        """
        ce_list_dict = [{str(cn): [ce.as_dict() if ce is not None else None for ce in ce_dict[cn]]
                         for cn in ce_dict} if ce_dict is not None else None for ce_dict in self.ce_list]
        nbs_sets_dict = [{str(cn): [nb_set.as_dict() for nb_set in nb_sets]
                          for cn, nb_sets in site_nbs_sets.items()}
                         if site_nbs_sets is not None else None
                         for site_nbs_sets in self.neighbors_sets]
        info_dict = {key: val for key, val in self.info.items() if key not in ['sites_info']}
        info_dict['sites_info'] = [{'nb_sets_info': {str(cn): {str(inb_set): nb_set_info
                                                               for inb_set, nb_set_info in cn_sets.items()}
                                                     for cn, cn_sets in site_info['nb_sets_info'].items()},
                                    'time': site_info['time']} if 'nb_sets_info' in site_info else {}
                                   for site_info in self.info['sites_info']]

        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "voronoi": self.voronoi.as_dict(),
                "valences": self.valences,
                "sites_map": self.sites_map,
                "equivalent_sites": [[ps.as_dict() for ps in psl] for psl in self.equivalent_sites],
                "ce_list": ce_list_dict,
                "structure": self.structure.as_dict(),
                "neighbors_sets": nbs_sets_dict,
                "info": info_dict}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the StructureEnvironments object from a dict representation of the StructureEnvironments created
        using the as_dict method.
        :param d: dict representation of the StructureEnvironments object
        :return: StructureEnvironments object
        """
        ce_list = [None if (ce_dict == 'None' or ce_dict is None) else {
            int(cn): [None if (ced is None or ced == 'None') else
                      ChemicalEnvironments.from_dict(ced) for ced in ce_dict[cn]]
            for cn in ce_dict} for ce_dict in d['ce_list']]
        voronoi = DetailedVoronoiContainer.from_dict(d['voronoi'])
        structure = Structure.from_dict(d['structure'])
        neighbors_sets = [{int(cn): [cls.NeighborsSet.from_dict(dd=nb_set_dict,
                                                                structure=structure,
                                                                detailed_voronoi=voronoi)
                                     for nb_set_dict in nb_sets]
                           for cn, nb_sets in site_nbs_sets_dict.items()}
                          if site_nbs_sets_dict is not None else None
                          for site_nbs_sets_dict in d['neighbors_sets']]
        info = {key: val for key, val in d['info'].items() if key not in ['sites_info']}
        if 'sites_info' in d['info']:
            info['sites_info'] = [{'nb_sets_info': {int(cn): {int(inb_set): nb_set_info
                                                              for inb_set, nb_set_info in cn_sets.items()}
                                                    for cn, cn_sets in site_info['nb_sets_info'].items()},
                                   'time': site_info['time']} if 'nb_sets_info' in site_info else {}
                                  for site_info in d['info']['sites_info']]
        return cls(voronoi=voronoi, valences=d['valences'],
                   sites_map=d['sites_map'],
                   equivalent_sites=[[PeriodicSite.from_dict(psd) for psd in psl] for psl in d['equivalent_sites']],
                   ce_list=ce_list, structure=structure,
                   neighbors_sets=neighbors_sets,
                   info=info)


class LightStructureEnvironments(MSONable):
    """
    Class used to store the chemical environments of a given structure obtained from a given ChemenvStrategy. Currently,
    only strategies leading to the determination of a unique environment for each site is allowed
    This class does not store all the information contained in the StructureEnvironments object, only the coordination
    environment found
    """
    DELTA_MAX_OXIDATION_STATE = 0.1
    DEFAULT_STATISTICS_FIELDS = ['anion_list', 'anion_atom_list', 'cation_list', 'cation_atom_list',
                                 'neutral_list', 'neutral_atom_list',
                                 'atom_coordination_environments_present',
                                 'ion_coordination_environments_present',
                                 'fraction_atom_coordination_environments_present',
                                 'fraction_ion_coordination_environments_present',
                                 'coordination_environments_atom_present',
                                 'coordination_environments_ion_present']

    class NeighborsSet():
        """
        Class used to store a given set of neighbors of a given site (based on a list of sites, the voronoi
        container is not part of the LightStructureEnvironments object).
        """
        def __init__(self, structure, isite, all_nbs_sites, all_nbs_sites_indices):
            self.structure = structure
            self.isite = isite
            self.all_nbs_sites = all_nbs_sites
            myset = set(all_nbs_sites_indices)
            if len(myset) != len(all_nbs_sites_indices):
                raise ValueError('Set of neighbors contains duplicates !')
            self.all_nbs_sites_indices = sorted(myset)
            self.all_nbs_sites_indices_unsorted = all_nbs_sites_indices
            self.all_nbs_sites_indices_and_image = []

        @property
        def neighb_coords(self):
            return [self.all_nbs_sites[inb]['site'].coords for inb in self.all_nbs_sites_indices_unsorted]

        @property
        def neighb_sites(self):
            return [self.all_nbs_sites[inb]['site'] for inb in self.all_nbs_sites_indices_unsorted]

        @property
        def neighb_sites_and_indices(self):
            return [{'site': self.all_nbs_sites[inb]['site'],
                     'index': self.all_nbs_sites[inb]['index']} for inb in self.all_nbs_sites_indices_unsorted]

        @property
        def neighb_indices_and_images(self):
            return [{'index': self.all_nbs_sites[inb]['index'],
                     'image_cell': self.all_nbs_sites[inb]['image_cell']}
                    for inb in self.all_nbs_sites_indices_unsorted]

        def __len__(self):
            return len(self.all_nbs_sites_indices)

        def __hash__(self):
            return len(self.all_nbs_sites_indices)

        def __eq__(self, other):
            return self.isite == other.isite and self.all_nbs_sites_indices == other.all_nbs_sites_indices

        def __ne__(self, other):
            return not self == other

        def __str__(self):
            out = 'Neighbors Set for site #{:d} :\n'.format(self.isite)
            out += ' - Coordination number : {:d}\n'.format(len(self))
            out += ' - Neighbors sites indices : {}' \
                   '\n'.format(', '.join(['{:d}'.format(nb_list_index)
                                          for nb_list_index in self.all_nbs_sites_indices]))
            return out

        def as_dict(self):
            return {'isite': self.isite,
                    'all_nbs_sites_indices': self.all_nbs_sites_indices_unsorted}
                    # 'all_nbs_sites_indices_unsorted': self.all_nbs_sites_indices_unsorted}

        @classmethod
        def from_dict(cls, dd, structure, all_nbs_sites):
            return cls(structure=structure,
                       isite=dd['isite'],
                       all_nbs_sites=all_nbs_sites,
                       all_nbs_sites_indices=dd['all_nbs_sites_indices'])

    def __init__(self, strategy,
                 coordination_environments=None, all_nbs_sites=None, neighbors_sets=None,
                 structure=None, valences=None, valences_origin=None):
        """
        Constructor for the LightStructureEnvironments object.
        """
        self.strategy = strategy
        self.statistics_dict = None
        self.coordination_environments = coordination_environments
        self._all_nbs_sites = all_nbs_sites
        self.neighbors_sets = neighbors_sets
        self.structure = structure
        self.valences = valences
        self.valences_origin = valences_origin

    @classmethod
    def from_structure_environments(cls, strategy, structure_environments, valences=None, valences_origin=None):
        structure = structure_environments.structure
        strategy.set_structure_environments(structure_environments=structure_environments)
        coordination_environments = [None] * len(structure)
        neighbors_sets = [None] * len(structure)
        _all_nbs_sites = []
        my_all_nbs_sites = []
        if valences is None:
            valences = structure_environments.valences
            if valences_origin is None:
                valences_origin = 'from_structure_environments'
        else:
            if valences_origin is None:
                valences_origin = 'user-specified'

        for isite, site in enumerate(structure):
            site_ces_and_nbs_list = strategy.get_site_ce_fractions_and_neighbors(site, strategy_info=True)
            if site_ces_and_nbs_list is None:
                continue
            coordination_environments[isite] = []
            neighbors_sets[isite] = []
            site_ces = []
            site_nbs_sets = []
            for ce_and_neighbors in site_ces_and_nbs_list:
                _all_nbs_sites_indices = []
                # Coordination environment
                ce_dict = {'ce_symbol': ce_and_neighbors['ce_symbol'],
                           'ce_fraction': ce_and_neighbors['ce_fraction']}
                if ce_and_neighbors['ce_dict'] is not None:
                    csm = ce_and_neighbors['ce_dict']['other_symmetry_measures'][strategy.symmetry_measure_type]
                else:
                    csm = None
                ce_dict['csm'] = csm
                ce_dict['permutation'] = ce_and_neighbors['ce_dict']['permutation']
                site_ces.append(ce_dict)
                # Neighbors
                neighbors = ce_and_neighbors['neighbors']
                for nb_site_and_index in neighbors:
                    nb_site = nb_site_and_index['site']
                    try:
                        nb_allnbs_sites_index = my_all_nbs_sites.index(nb_site)
                    except ValueError:
                        nb_index_unitcell = nb_site_and_index['index']
                        diff = nb_site.frac_coords - structure[nb_index_unitcell].frac_coords
                        rounddiff = np.round(diff)
                        if not np.allclose(diff, rounddiff):
                            raise ValueError('Weird, differences between one site in a periodic image cell is not '
                                             'integer ...')
                        nb_image_cell = np.array(rounddiff, np.int)
                        nb_allnbs_sites_index = len(_all_nbs_sites)
                        _all_nbs_sites.append({'site': nb_site,
                                               'index': nb_index_unitcell,
                                               'image_cell': nb_image_cell})
                        my_all_nbs_sites.append(nb_site)
                    _all_nbs_sites_indices.append(nb_allnbs_sites_index)

                nb_set = cls.NeighborsSet(structure=structure, isite=isite,
                                          all_nbs_sites=_all_nbs_sites,
                                          all_nbs_sites_indices=_all_nbs_sites_indices)
                site_nbs_sets.append(nb_set)
            coordination_environments[isite] = site_ces
            neighbors_sets[isite] = site_nbs_sets
        return cls(strategy=strategy,
                   coordination_environments=coordination_environments,
                   all_nbs_sites=_all_nbs_sites,
                   neighbors_sets=neighbors_sets,
                   structure=structure, valences=valences,
                   valences_origin=valences_origin)

    def setup_statistic_lists(self):
        self.statistics_dict = {'valences_origin': self.valences_origin,
                                'anion_list': {},                                                    # OK
                                'anion_number': None,                                                # OK
                                'anion_atom_list': {},                                               # OK
                                'anion_atom_number': None,                                           # OK
                                'cation_list': {},                                                   # OK
                                'cation_number': None,                                               # OK
                                'cation_atom_list': {},                                              # OK
                                'cation_atom_number': None,                                          # OK
                                'neutral_list': {},                                                  # OK
                                'neutral_number': None,                                              # OK
                                'neutral_atom_list': {},                                             # OK
                                'neutral_atom_number': None,                                         # OK
                                'atom_coordination_environments_present': {},                        # OK
                                'ion_coordination_environments_present': {},                         # OK
                                'coordination_environments_ion_present': {},                         # OK
                                'coordination_environments_atom_present': {},                        # OK
                                'fraction_ion_coordination_environments_present': {},                # OK
                                'fraction_atom_coordination_environments_present': {},               # OK
                                'fraction_coordination_environments_ion_present': {},                # OK
                                'fraction_coordination_environments_atom_present': {},               # OK
                                'count_ion_present': {},                                             # OK
                                'count_atom_present': {},                                            # OK
                                'count_coordination_environments_present': {}}
        atom_stat = self.statistics_dict['atom_coordination_environments_present']
        ce_atom_stat = self.statistics_dict['coordination_environments_atom_present']
        fraction_atom_stat = self.statistics_dict['fraction_atom_coordination_environments_present']
        fraction_ce_atom_stat = self.statistics_dict['fraction_coordination_environments_atom_present']
        count_atoms = self.statistics_dict['count_atom_present']
        count_ce = self.statistics_dict['count_coordination_environments_present']
        for isite, site in enumerate(self.structure):
            # Building anion and cation list
            site_species = []
            if self.valences != 'undefined':
                for sp, occ in site.species_and_occu.items():
                    valence = self.valences[isite]
                    strspecie = str(Specie(sp.symbol, valence))
                    if valence < 0:
                        specielist = self.statistics_dict['anion_list']
                        atomlist = self.statistics_dict['anion_atom_list']
                    elif valence > 0:
                        specielist = self.statistics_dict['cation_list']
                        atomlist = self.statistics_dict['cation_atom_list']
                    else:
                        specielist = self.statistics_dict['neutral_list']
                        atomlist = self.statistics_dict['neutral_atom_list']
                    if strspecie not in specielist:
                        specielist[strspecie] = occ
                    else:
                        specielist[strspecie] += occ
                    if sp.symbol not in atomlist:
                        atomlist[sp.symbol] = occ
                    else:
                        atomlist[sp.symbol] += occ
                    site_species.append((sp.symbol, valence, occ))
            # Building environments lists
            if self.coordination_environments[isite] is not None:
                site_envs = [(ce_piece_dict['ce_symbol'], ce_piece_dict['ce_fraction'])
                             for ce_piece_dict in self.coordination_environments[isite]]
                for ce_symbol, fraction in site_envs:
                    if fraction is None:
                        continue
                    if ce_symbol not in count_ce:
                        count_ce[ce_symbol] = 0.0
                    count_ce[ce_symbol] += fraction
                for sp, occ in site.species_and_occu.items():
                    elmt = sp.symbol
                    if elmt not in atom_stat:
                        atom_stat[elmt] = {}
                        count_atoms[elmt] = 0.0
                    count_atoms[elmt] += occ
                    for ce_symbol, fraction in site_envs:
                        if fraction is None:
                            continue
                        if ce_symbol not in atom_stat[elmt]:
                            atom_stat[elmt][ce_symbol] = 0.0

                        atom_stat[elmt][ce_symbol] += occ * fraction
                        if ce_symbol not in ce_atom_stat:
                            ce_atom_stat[ce_symbol] = {}
                        if elmt not in ce_atom_stat[ce_symbol]:
                            ce_atom_stat[ce_symbol][elmt] = 0.0
                        ce_atom_stat[ce_symbol][elmt] += occ * fraction

                if self.valences != 'undefined':
                    ion_stat = self.statistics_dict['ion_coordination_environments_present']
                    ce_ion_stat = self.statistics_dict['coordination_environments_ion_present']
                    count_ions = self.statistics_dict['count_ion_present']
                    for elmt, oxi_state, occ in site_species:
                        if elmt not in ion_stat:
                            ion_stat[elmt] = {}
                            count_ions[elmt] = {}
                        if oxi_state not in ion_stat[elmt]:
                            ion_stat[elmt][oxi_state] = {}
                            count_ions[elmt][oxi_state] = 0.0
                        count_ions[elmt][oxi_state] += occ
                        for ce_symbol, fraction in site_envs:
                            if fraction is None:
                                continue
                            if ce_symbol not in ion_stat[elmt][oxi_state]:
                                ion_stat[elmt][oxi_state][ce_symbol] = 0.0
                            ion_stat[elmt][oxi_state][ce_symbol] += occ * fraction
                            if ce_symbol not in ce_ion_stat:
                                ce_ion_stat[ce_symbol] = {}
                            if elmt not in ce_ion_stat[ce_symbol]:
                                ce_ion_stat[ce_symbol][elmt] = {}
                            if oxi_state not in ce_ion_stat[ce_symbol][elmt]:
                                ce_ion_stat[ce_symbol][elmt][oxi_state] = 0.0
                            ce_ion_stat[ce_symbol][elmt][oxi_state] += occ * fraction
        self.statistics_dict['anion_number'] = len(self.statistics_dict['anion_list'])
        self.statistics_dict['anion_atom_number'] = len(self.statistics_dict['anion_atom_list'])
        self.statistics_dict['cation_number'] = len(self.statistics_dict['cation_list'])
        self.statistics_dict['cation_atom_number'] = len(self.statistics_dict['cation_atom_list'])
        self.statistics_dict['neutral_number'] = len(self.statistics_dict['neutral_list'])
        self.statistics_dict['neutral_atom_number'] = len(self.statistics_dict['neutral_atom_list'])

        for elmt, envs in atom_stat.items():
            sumelement = count_atoms[elmt]
            fraction_atom_stat[elmt] = {env: fraction / sumelement for env, fraction in envs.items()}
        for ce_symbol, atoms in ce_atom_stat.items():
            sumsymbol = count_ce[ce_symbol]
            fraction_ce_atom_stat[ce_symbol] = {atom: fraction / sumsymbol for atom, fraction in atoms.items()}
        ion_stat = self.statistics_dict['ion_coordination_environments_present']
        fraction_ion_stat = self.statistics_dict['fraction_ion_coordination_environments_present']
        ce_ion_stat = self.statistics_dict['coordination_environments_ion_present']
        fraction_ce_ion_stat = self.statistics_dict['fraction_coordination_environments_ion_present']
        count_ions = self.statistics_dict['count_ion_present']
        for elmt, oxi_states_envs in ion_stat.items():
            fraction_ion_stat[elmt] = {}
            for oxi_state, envs in oxi_states_envs.items():
                sumspecie = count_ions[elmt][oxi_state]
                fraction_ion_stat[elmt][oxi_state] = {env: fraction / sumspecie
                                                      for env, fraction in envs.items()}
        for ce_symbol, ions in ce_ion_stat.items():
            fraction_ce_ion_stat[ce_symbol] = {}
            sum_ce = np.sum([np.sum(list(oxistates.values())) for elmt, oxistates in ions.items()])
            for elmt, oxistates in ions.items():
                fraction_ce_ion_stat[ce_symbol][elmt] = {oxistate: fraction / sum_ce
                                                         for oxistate, fraction in oxistates.items()}

    def get_site_info_for_specie_ce(self, specie, ce_symbol, min_fraction=0.0):
        element = specie.symbol
        oxi_state = specie.oxi_state
        isites = []
        csms = []
        fractions = []
        for isite, site in enumerate(self.structure):
            if element in [sp.symbol for sp in site.species_and_occu]:
                if self.valences == 'undefined' or oxi_state == self.valences[isite]:
                    for ce_dict in self.coordination_environments[isite]:
                        if ce_symbol == ce_dict['ce_symbol']:
                            isites.append(isite)
                            csms.append(ce_dict['csm'])
                            fractions.append(ce_dict['ce_fraction'])
        return {'isites': isites, 'fractions': fractions, 'csms': csms}

    def get_site_info_for_specie_allces(self, specie, min_fraction=0.0):
        allces = {}
        element = specie.symbol
        oxi_state = specie.oxi_state
        for isite, site in enumerate(self.structure):
            if element in [sp.symbol for sp in site.species_and_occu]:
                if self.valences == 'undefined' or oxi_state == self.valences[isite]:
                    if self.coordination_environments[isite] is None:
                        continue
                    for ce_dict in self.coordination_environments[isite]:
                        if ce_dict['ce_fraction'] < min_fraction:
                            continue
                        if ce_dict['ce_symbol'] not in allces:
                            allces[ce_dict['ce_symbol']] = {'isites': [], 'fractions': [], 'csms': []}
                        allces[ce_dict['ce_symbol']]['isites'].append(isite)
                        allces[ce_dict['ce_symbol']]['fractions'].append(ce_dict['ce_fraction'])
                        allces[ce_dict['ce_symbol']]['csms'].append(ce_dict['csm'])
        return allces

    def get_statistics(self, statistics_fields=DEFAULT_STATISTICS_FIELDS, bson_compatible=False):
        if self.statistics_dict is None:
            self.setup_statistic_lists()
        if statistics_fields == 'ALL':
            statistics_fields = [key for key in self.statistics_dict]
        if bson_compatible:
            dd = jsanitize({field: self.statistics_dict[field] for field in statistics_fields})
        else:
            dd = {field: self.statistics_dict[field] for field in statistics_fields}
        return dd

    def contains_only_one_anion_atom(self, anion_atom):
        return (len(self.statistics_dict['anion_atom_list']) == 1 and
                anion_atom in self.statistics_dict['anion_atom_list'])

    def contains_only_one_anion(self, anion):
        return len(self.statistics_dict['anion_list']) == 1 and anion in self.statistics_dict['anion_list']

    def site_contains_environment(self, isite, ce_symbol):
        if self.coordination_environments[isite] is None:
            return False
        return ce_symbol in [ce_dict['ce_symbol'] for ce_dict in self.coordination_environments[isite]]

    def site_has_clear_environment(self, isite, conditions=None):
        if self.coordination_environments[isite] is None:
            raise ValueError('Coordination environments have not been determined for site {:d}'.format(isite))
        if conditions is None:
            return len(self.coordination_environments[isite]) == 1
        ce = max(self.coordination_environments[isite], key=lambda x: x['ce_fraction'])
        for condition in conditions:
            target = condition['target']
            if target == 'ce_fraction':
                if ce[target] < condition['minvalue']:
                    return False
            elif target == 'csm':
                if ce[target] > condition['maxvalue']:
                    return False
            elif target == 'number_of_ces':
                if ce[target] > condition['maxnumber']:
                    return False
            else:
                raise ValueError('Target "{}" for condition of clear environment is not allowed'.format(target))
            pass
        return True

    def structure_has_clear_environments(self, conditions=None, skip_none=True, skip_empty=False):
        for isite in range(len(self.structure)):
            if self.coordination_environments[isite] is None:
                if skip_none:
                    continue
                else:
                    return False
            if len(self.coordination_environments[isite]) == 0:
                if skip_empty:
                    continue
                else:
                    return False
            if not self.site_has_clear_environment(isite=isite, conditions=conditions):
                return False
        return True

    def clear_environments(self, conditions=None):
        clear_envs_list = set()
        for isite in range(len(self.structure)):
            if self.coordination_environments[isite] is None:
                continue
            if len(self.coordination_environments[isite]) == 0:
                continue
            if self.site_has_clear_environment(isite=isite, conditions=conditions):
                ce = max(self.coordination_environments[isite], key=lambda x: x['ce_fraction'])
                clear_envs_list.add(ce['ce_symbol'])
        return list(clear_envs_list)

    def structure_contains_atom_environment(self, atom_symbol, ce_symbol):
        """
        Checks whether the structure contains a given atom in a given environment
        :param atom_symbol: Symbol of the atom
        :param ce_symbol: Symbol of the coordination environment
        :return: True if the coordination environment is found, False otherwise
        """
        for isite, site in enumerate(self.structure):
            if (Element(atom_symbol) in site.species_and_occu.
                    element_composition and self.site_contains_environment(isite, ce_symbol)):
                return True
        return False

    @property
    def uniquely_determines_coordination_environments(self):
        """
        True if the coordination environments are uniquely determined.
        """
        return self.strategy.uniquely_determines_coordination_environments

    def __eq__(self, other):
        """
        Equality method that checks if the LightStructureEnvironments object is equal to another
        LightStructureEnvironments object. Two LightStructureEnvironments objects are equal if the strategy used
        is the same, if the structure is the same, if the valences used in the strategies are the same, if the
        coordination environments and the neighbours determined by the strategy are the same
        :param other: LightStructureEnvironments object to compare with
        :return: True if both objects are equal, False otherwise
        """
        is_equal = (self.strategy == other.strategy and
                    self.structure == other.structure and
                    self.coordination_environments == other.coordination_environments and
                    self.valences == other.valences and
                    self.neighbors_sets == other.neighbors_sets)
        this_sites = [ss['site'] for ss in self._all_nbs_sites]
        other_sites = [ss['site'] for ss in other._all_nbs_sites]
        this_indices = [ss['index'] for ss in self._all_nbs_sites]
        other_indices = [ss['index'] for ss in other._all_nbs_sites]
        return (is_equal and this_sites == other_sites and this_indices == other_indices)

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        """
        Bson-serializable dict representation of the LightStructureEnvironments object.
        :return: Bson-serializable dict representation of the LightStructureEnvironments object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "strategy": self.strategy.as_dict(),
                "structure": self.structure.as_dict(),
                "coordination_environments": self.coordination_environments,
                "all_nbs_sites": [{'site': nb_site['site'].as_dict(),
                                   'index': nb_site['index'],
                                   'image_cell': [int(ii) for ii in nb_site['image_cell']]}
                                  for nb_site in self._all_nbs_sites],
                "neighbors_sets": [[nb_set.as_dict() for nb_set in site_nb_sets] if site_nb_sets is not None else None
                                   for site_nb_sets in self.neighbors_sets],
                "valences": self.valences}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the LightStructureEnvironments object from a dict representation of the
        LightStructureEnvironments created using the as_dict method.
        :param d: dict representation of the LightStructureEnvironments object
        :return: LightStructureEnvironments object
        """
        dec = MontyDecoder()
        structure = dec.process_decoded(d['structure'])
        all_nbs_sites = []
        for nb_site in d['all_nbs_sites']:
            site = dec.process_decoded(nb_site['site'])
            if 'image_cell' in nb_site:
                image_cell = np.array(nb_site['image_cell'], np.int)
            else:
                diff = site.frac_coords - structure[nb_site['index']].frac_coords
                rounddiff = np.round(diff)
                if not np.allclose(diff, rounddiff):
                    raise ValueError('Weird, differences between one site in a periodic image cell is not '
                                     'integer ...')
                image_cell = np.array(rounddiff, np.int)
            all_nbs_sites.append({'site': site,
                                  'index': nb_site['index'],
                                  'image_cell': image_cell})
        neighbors_sets = [[cls.NeighborsSet.from_dict(dd=nb_set, structure=structure,
                                                      all_nbs_sites=all_nbs_sites)
                           for nb_set in site_nb_sets] if site_nb_sets is not None else None
                          for site_nb_sets in d['neighbors_sets']]
        return cls(strategy=dec.process_decoded(d['strategy']),
                   coordination_environments=d['coordination_environments'],
                   all_nbs_sites=all_nbs_sites,
                   neighbors_sets=neighbors_sets,
                   structure=structure,
                   valences=d['valences'])

class ChemicalEnvironments(MSONable):
    """
    Class used to store all the information about the chemical environment of a given site for a given list of
    coordinated neighbours (internally called "cn_map")
    """

    def __init__(self, coord_geoms=None):
        """
        Initializes the ChemicalEnvironments object containing all the information about the chemical
        environment of a given site
        :param coord_geoms: coordination geometries to be added to the chemical environment.
        """
        if coord_geoms is None:
            self.coord_geoms = {}
        else:
            raise NotImplementedError('Constructor for ChemicalEnvironments with the coord_geoms argument is not'
                                      'yet implemented')

    def __getitem__(self, mp_symbol):
        if not mp_symbol in self.coord_geoms:
            raise IndexError()
        return self.coord_geoms[mp_symbol]

    def __len__(self):
        """
        Returns the number of coordination geometries in this ChemicalEnvironments object
        :return: Number of coordination geometries in this ChemicalEnvironments object
        """
        return len(self.coord_geoms)

    def __iter__(self):
        for cg, cg_dict in self.coord_geoms.items():
            yield (cg, cg_dict)

    def minimum_geometry(self, symmetry_measure_type=None, max_csm=None):
        """
        Returns the geometry with the minimum continuous symmetry measure of this ChemicalEnvironments
        :return: tuple (symbol, csm) with symbol being the geometry with the minimum continuous symmetry measure and
        csm being the continuous symmetry measure associted to it
        :raise: ValueError if no coordination geometry is found in this ChemicalEnvironments object
        """
        if len(self.coord_geoms) == 0:
            return None
        cglist = [cg for cg in self.coord_geoms]
        if symmetry_measure_type is None:
            csms = np.array([self.coord_geoms[cg]['other_symmetry_measures']['csm_wcs_ctwcc'] for cg in cglist])
        else:
            csms = np.array([self.coord_geoms[cg]['other_symmetry_measures'][symmetry_measure_type] for cg in cglist])
        csmlist = [self.coord_geoms[cg] for cg in cglist]
        imin = np.argmin(csms)
        if max_csm is not None:
            if csmlist[imin] > max_csm:
                return None
        return cglist[imin], csmlist[imin]

    def minimum_geometries(self, n=None, symmetry_measure_type=None, max_csm=None):
        """
        Returns a list of geometries with increasing continuous symmetry measure in this ChemicalEnvironments object
        :param n: Number of geometries to be included in the list
        :return: list of geometries with increasing continuous symmetry measure in this ChemicalEnvironments object
        :raise: ValueError if no coordination geometry is found in this ChemicalEnvironments object
        """
        cglist = [cg for cg in self.coord_geoms]
        if symmetry_measure_type is None:
            csms = np.array([self.coord_geoms[cg]['other_symmetry_measures']['csm_wcs_ctwcc'] for cg in cglist])
        else:
            csms = np.array([self.coord_geoms[cg]['other_symmetry_measures'][symmetry_measure_type] for cg in cglist])
        csmlist = [self.coord_geoms[cg] for cg in cglist]
        isorted = np.argsort(csms)
        if max_csm is not None:
            if n is None:
                return [(cglist[ii], csmlist[ii]) for ii in isorted if csms[ii] <= max_csm]
            else:
                return [(cglist[ii], csmlist[ii]) for ii in isorted[:n] if csms[ii] <= max_csm]
        else:
            if n is None:
                return [(cglist[ii], csmlist[ii]) for ii in isorted]
            else:
                return [(cglist[ii], csmlist[ii]) for ii in isorted[:n]]

    def add_coord_geom(self, mp_symbol, symmetry_measure, algo='UNKNOWN', permutation=None, override=False,
                       local2perfect_map=None, perfect2local_map=None, detailed_voronoi_index=None,
                       other_symmetry_measures=None, rotation_matrix=None, scaling_factor=None):
        """
        Adds a coordination geometry to the ChemicalEnvironments object
        :param mp_symbol: Symbol (internal) of the coordination geometry added
        :param symmetry_measure: Symmetry measure of the coordination geometry added
        :param algo: Algorithm used for the search of the coordination geometry added
        :param permutation: Permutation of the neighbors that leads to the csm stored
        :param override: If set to True, the coordination geometry will override the existent one if present
        :return: :raise: ChemenvError if the coordination geometry is already added and override is set to False
        """
        if not allcg.is_a_valid_coordination_geometry(mp_symbol=mp_symbol):
            raise ChemenvError(self.__class__,
                               'add_coord_geom',
                               'Coordination geometry with mp_symbol "{mp}" is not valid'
                               .format(mp=mp_symbol))
        if mp_symbol in list(self.coord_geoms.keys()) and not override:
            raise ChemenvError(self.__class__,
                               "add_coord_geom",
                               "This coordination geometry is already present and override is set to False")
        else:
            self.coord_geoms[mp_symbol] = {'symmetry_measure': float(symmetry_measure), 'algo': algo,
                                           'permutation': [int(i) for i in permutation],
                                           'local2perfect_map': local2perfect_map,
                                           'perfect2local_map': perfect2local_map,
                                           'detailed_voronoi_index': detailed_voronoi_index,
                                           'other_symmetry_measures': other_symmetry_measures,
                                           'rotation_matrix': rotation_matrix,
                                           'scaling_factor': scaling_factor}

    def __str__(self):
        """
        Returns a string representation of the ChemicalEnvironments object
        :return: String representation of the ChemicalEnvironments object
        """
        out = 'Chemical environments object :\n'
        if len(self.coord_geoms) == 0:
            out += ' => No coordination in it <=\n'
            return out
        for key in self.coord_geoms.keys():
            mp_symbol = key
            break
        cn = symbol_cn_mapping[mp_symbol]
        out += ' => Coordination {} <=\n'.format(cn)
        mp_symbols = list(self.coord_geoms.keys())
        csms_wcs = [self.coord_geoms[mp_symbol]['other_symmetry_measures']['csm_wcs_ctwcc'] for mp_symbol in mp_symbols]
        icsms_sorted = np.argsort(csms_wcs)
        mp_symbols = [mp_symbols[ii] for ii in icsms_sorted]
        for mp_symbol in mp_symbols:
            csm_wcs = self.coord_geoms[mp_symbol]['other_symmetry_measures']['csm_wcs_ctwcc']
            csm_wocs = self.coord_geoms[mp_symbol]['other_symmetry_measures']['csm_wocs_ctwocc']
            out += '   - {}\n'.format(mp_symbol)
            out += '      csm1 (with central site) : {}'.format(csm_wcs)
            out += '      csm2 (without central site) : {}'.format(csm_wocs)
            out += '     algo : {}'.format(self.coord_geoms[mp_symbol]['algo'])
            out += '     perm : {}\n'.format(self.coord_geoms[mp_symbol]['permutation'])
            out += '       local2perfect : {}\n'.format(str(self.coord_geoms[mp_symbol]['local2perfect_map']))
            out += '       perfect2local : {}\n'.format(str(self.coord_geoms[mp_symbol]['perfect2local_map']))
        return out

    def is_close_to(self, other, rtol=0.0, atol=1e-8):
        if set(self.coord_geoms.keys()) != set(other.coord_geoms.keys()):
            return False
        for mp_symbol, cg_dict_self in self.coord_geoms.items():
            cg_dict_other = other[mp_symbol]
            other_csms_self = cg_dict_self['other_symmetry_measures']
            other_csms_other = cg_dict_other['other_symmetry_measures']
            for csmtype in ['csm_wcs_ctwcc', 'csm_wcs_ctwocc', 'csm_wcs_csc',
                            'csm_wocs_ctwcc', 'csm_wocs_ctwocc', 'csm_wocs_csc']:
                if not np.isclose(other_csms_self[csmtype], other_csms_other[csmtype], rtol=rtol, atol=atol):
                    return False
        return True

    def __eq__(self, other):
        """
        Equality method that checks if the ChemicalEnvironments object is equal to another ChemicalEnvironments
        object.
        :param other: ChemicalEnvironments object to compare with
        :return: True if both objects are equal, False otherwise
        """
        if set(self.coord_geoms.keys()) != set(other.coord_geoms.keys()):
            return False
        for mp_symbol, cg_dict_self in self.coord_geoms.items():
            cg_dict_other = other.coord_geoms[mp_symbol]
            if cg_dict_self['symmetry_measure'] != cg_dict_other['symmetry_measure']:
                return False
            if cg_dict_self['algo'] != cg_dict_other['algo']:
                return False
            if cg_dict_self['permutation'] != cg_dict_other['permutation']:
                return False
            if cg_dict_self['detailed_voronoi_index'] != cg_dict_other['detailed_voronoi_index']:
                return False
            other_csms_self = cg_dict_self['other_symmetry_measures']
            other_csms_other = cg_dict_other['other_symmetry_measures']
            for csmtype in ['csm_wcs_ctwcc', 'csm_wcs_ctwocc', 'csm_wcs_csc',
                            'csm_wocs_ctwcc', 'csm_wocs_ctwocc', 'csm_wocs_csc']:
                if other_csms_self[csmtype] != other_csms_other[csmtype]:
                    return False
        return True

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        """
        Returns a dictionary representation of the ChemicalEnvironments object
        :return:
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "coord_geoms": jsanitize(self.coord_geoms)}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the ChemicalEnvironments object from a dict representation of the ChemicalEnvironments created
        using the as_dict method.
        :param d: dict representation of the ChemicalEnvironments object
        :return: ChemicalEnvironments object
        """
        ce = cls()
        for cg in d['coord_geoms'].keys():
            if d['coord_geoms'][cg]['local2perfect_map'] is None:
                l2p_map = None
            else:
                l2p_map = {int(key): int(val) for key, val in d['coord_geoms'][cg]['local2perfect_map'].items()}
            if d['coord_geoms'][cg]['perfect2local_map'] is None:
                p2l_map = None
            else:
                p2l_map = {int(key): int(val) for key, val in d['coord_geoms'][cg]['perfect2local_map'].items()}
            if ('other_symmetry_measures' in d['coord_geoms'][cg] and
                        d['coord_geoms'][cg]['other_symmetry_measures'] is not None):
                other_csms = d['coord_geoms'][cg]['other_symmetry_measures']
            else:
                other_csms = None
            ce.add_coord_geom(cg,
                              d['coord_geoms'][cg]['symmetry_measure'],
                              d['coord_geoms'][cg]['algo'],
                              permutation=d['coord_geoms'][cg]['permutation'],
                              local2perfect_map=l2p_map,
                              perfect2local_map=p2l_map,
                              detailed_voronoi_index=d['coord_geoms'][cg]['detailed_voronoi_index'],
                              other_symmetry_measures=other_csms,
                              rotation_matrix=d['coord_geoms'][cg]['rotation_matrix'],
                              scaling_factor=d['coord_geoms'][cg]['scaling_factor'])
        return ce