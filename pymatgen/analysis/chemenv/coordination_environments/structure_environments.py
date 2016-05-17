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

    def __init__(self, voronoi, bva_valences, sites_map, equivalent_sites,
                 ce_list, structure):
        """
        Constructor for the StructureEnvironments object.
        :param voronoi: VoronoiContainer object for the structure
        :param bva_valences: Valences obtained using the bond-valence analysis
        :param sites_map: Mapping of equivalent sites to the unequivalent sites that have been computed.
        :param equivalent_sites: List of list of equivalent sites of the structure
        :param struct_sites_to_irreducible_site_list_map: Maps the index of a site to the index of the item in the
        list of equivalent sites to which the site belongs.
        :param ce_list: List of chemical environments
        :param structure: Structure object
        """
        self.voronoi = voronoi
        self.bva_valences = bva_valences
        self.sites_map = sites_map
        self.equivalent_sites = equivalent_sites
        #self.struct_sites_to_irreducible_site_list_map = struct_sites_to_irreducible_site_list_map
        self.ce_list = ce_list
        self.structure = structure

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

    def get_environments_figure(self, isite, plot_type=None, title='Coordination numbers', max_dist=2.0,
                                additional_condition=AC.ONLY_ACB, colormap=None, figsize=None):
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
            plot_type = {'distance_parameter': ('initial_inverse_opposite', None),
                         'angle_parameter': ('initial_opposite', None)}
        bounds_and_limits = self.voronoi.voronoi_parameters_bounds_and_limits(isite, plot_type, max_dist)
        # distance_bounds = bounds_and_limits['distance_bounds']
        # angle_bounds = bounds_and_limits['angle_bounds']
        # dist_limits = bounds_and_limits['distance_limits']
        # ang_limits = bounds_and_limits['angle_limits']
        # iac = self.voronoi.additional_conditions.index(additional_condition)
        if colormap is None:
            mycm = cm.jet
        else:
            mycm = colormap
        mymin = 0.0
        mymax = 10.0
        norm = Normalize(vmin=mymin, vmax=mymax)
        scalarmap = cm.ScalarMappable(norm=norm, cmap=mycm)

        #Plot the rectangles and coordinations
        maps_and_vertices = self.voronoi.maps_and_surface_vertices(isite,
                                                                   additional_condition=additional_condition,
                                                                   plot_type=plot_type, max_dist=max_dist)
        if maps_and_vertices is None:
            return None
        cn_maps = maps_and_vertices['cn_maps']
        bounds_and_limits = maps_and_vertices['bounds_and_limits']
        dist_limits = bounds_and_limits['distance_limits']
        ang_limits = bounds_and_limits['angle_limits']
        vertices_dist_ang = maps_and_vertices['vertices_dist_ang']
        text_info_dist_ang = maps_and_vertices['text_info_dist_ang']
        for i_cn_map, cn_map in enumerate(cn_maps):
            ce = self.ce_list[isite][cn_map[0]][cn_map[1]]
            mingeom = ce.minimum_geometry()
            if mingeom is not None:
                mp_symbol = mingeom[0]
                csm = mingeom[1]['symmetry_measure']
                mycolor = scalarmap.to_rgba(csm)
                myinvcolor = [1.0 - mycolor[0], 1.0 - mycolor[1], 1.0 - mycolor[2], 1.0]
                mytext = '{}'.format(mp_symbol)
            else:
                cn = cn_map[0]
                mycolor = 'w'
                myinvcolor = 'k'
                mytext = '{:d}'.format(cn)
            polygon = Polygon(vertices_dist_ang[i_cn_map], closed=True, edgecolor='k', facecolor=mycolor)
            subplot.add_patch(polygon)
            subplot.annotate(mytext,
                             xy=text_info_dist_ang[i_cn_map],
                             ha='center', va='center', color=myinvcolor, fontsize='x-small')

        title += '\nDist: {}, Ang: {}'.format(plot_type['distance_parameter'][0], plot_type['angle_parameter'][0])
        subplot.set_title(title)
        subplot.set_xlabel('Distance factor')
        subplot.set_ylabel('Angle factor')
        subplot.set_xlim(dist_limits)
        subplot.set_ylim(ang_limits)
        # ax2 = subplot.twin()  # ax2 is responsible for "top" axis and "right" axis
        # ax2.set_xticks([0., .5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
        # ax2.set_xticklabels(["$0$", r"$\frac{1}{2}\pi$",
        #                      r"$\pi$", r"$\frac{3}{2}\pi$", r"$2\pi$"])
        #
        # ax2.axis["right"].major_ticklabels.set_visible(False)
        scalarmap.set_array([mymin, mymax])
        cb = fig.colorbar(scalarmap, ax=subplot, extend='max')
        cb.set_label('Continuous symmetry measure')
        return fig

    def plot_environments(self,  isite, plot_type=None, title='Coordination numbers', max_dist=2.0,
                          additional_condition=AC.ONLY_ACB, figsize=None):
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
        fig = self.get_environments_figure(isite=isite, plot_type=plot_type, title=title, max_dist=max_dist,
                                           additional_condition=additional_condition, figsize=figsize)
        if fig is None:
            return
        fig.show()

    def save_environments_figure(self,  isite, imagename='image.png', plot_type=None, title='Coordination numbers',
                                 max_dist=2.0, additional_condition=AC.ONLY_ACB, figsize=None):
        fig = self.get_environments_figure(isite=isite, plot_type=plot_type, title=title, max_dist=max_dist,
                                           additional_condition=additional_condition, figsize=figsize)
        if fig is None:
            return
        fig.savefig(imagename)

    def to_bv_dict(self, isite):
        out = {"neighbours_lists": [], "continuous_symmetry_measures": []}
        for cn, coordnbs_list in self.voronoi._unique_coordinated_neighbors[isite].items():
            for i_coordnbs, coordnbs in enumerate(coordnbs_list):
                neighbours = []
                for ineighb, neighb in enumerate(coordnbs[0]):
                    mydict = {'neighbour': neighb.as_dict(),
                              'distance': coordnbs[2][ineighb]['distance'],
                              'normalized_distance': coordnbs[2][ineighb]['weighted_distance'],
                              'angle': coordnbs[2][ineighb]['angle'],
                              'normalized_angle': coordnbs[2][ineighb]['weighted_angle']}
                    neighbours.append(mydict)
                ce = self.ce_list[isite][cn][i_coordnbs]
                mingeoms = ce.minimum_geometries()
                csm_dict = {ce: ce_dict['symmetry_measure'] for ce, ce_dict in mingeoms}
                out["neighbours_lists"].append(neighbours)
                out["continuous_symmetry_measures"].append(csm_dict)
        return out
      #           - "neighbour": the PeriodicSite object, as a dictionary
      # - "distance": the distance to the site X
      # - "angle": the angle from the site X
      # - "normalized_distance": the normalized distance to the site X
      # - "normalized_angle": the normalized angle from the site X

    def unique_coordinated_neighbors(self, isite=None, cn_map=None):
        return self.voronoi.unique_coordinated_neighbors(isite=isite, cn_map=cn_map)

    def __eq__(self, other):
        if len(self.ce_list) != len(other.ce_list):
            return False
        for ii in range(len(self.ce_list)):
            if self.ce_list[ii] is not None:
                if other.ce_list[ii] is None:
                    return False
                if set(self.ce_list[ii].keys()) != set(other.ce_list[ii].keys()):
                    return False
                for cn, ces in list(self.ce_list[ii].items()):
                    if ces != other.ce_list[ii][cn]:
                        return False
            else:
                if other.ce_list[ii] is not None:
                    return False
        return (self.voronoi == other.voronoi and self.bva_valences == other.bva_valences and
                self.sites_map == other.sites_map and self.equivalent_sites == other.equivalent_sites and
                self.ce_list == other.ce_list and self.structure == other.structure)

    def as_dict(self):
        """
        Bson-serializable dict representation of the StructureEnvironments object.
        :return: Bson-serializable dict representation of the StructureEnvironments object.
        """
        ce_list_dict = [{str(cn): [ce.as_dict() for ce in ce_dict[cn]]
                         for cn in ce_dict} if ce_dict is not None else None for ce_dict in self.ce_list]
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "voronoi": self.voronoi.as_dict(),
                "bva_valences": self.bva_valences,
                "sites_map": self.sites_map,
                "equivalent_sites": [[ps.as_dict() for ps in psl] for psl in self.equivalent_sites],
                "ce_list": ce_list_dict,
                "structure": self.structure.as_dict()}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the StructureEnvironments object from a dict representation of the StructureEnvironments created
        using the as_dict method.
        :param d: dict representation of the StructureEnvironments object
        :return: StructureEnvironments object
        """
        ce_list = [None if (ce_dict == 'None' or ce_dict is None) else {
            int(cn): [ChemicalEnvironments.from_dict(ced) for ced in ce_dict[cn]]
            for cn in ce_dict} for ce_dict in d['ce_list']]
        return cls(DetailedVoronoiContainer.from_dict(d['voronoi']), d['bva_valences'], d['sites_map'],
                   [[PeriodicSite.from_dict(psd) for psd in psl] for psl in d['equivalent_sites']],
                   ce_list, Structure.from_dict(d['structure']))


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
                                 'bva_ion_coordination_environments_present',
                                 'structure_ion_coordination_environments_present',
                                 'ion_coordination_environments_present',
                                 'fraction_atom_coordination_environments_present',
                                 'fraction_bva_ion_coordination_environments_present',
                                 'fraction_structure_ion_coordination_environments_present',
                                 'fraction_ion_coordination_environments_present',
                                 'coordination_environments_atom_present',
                                 'coordination_environments_ion_present']

    def __init__(self, strategy, structure_environments=None, bva_valences=None,
                 coordination_environments=None, neighbors=None, structure=None, setup_neighbors_by_indices=False,
                 neighbors_by_indices=None):
        """
        Constructor for the LightStructureEnvironments object.
        """
        self._strategy = strategy
        self._uniquely_determined_coordination_environments = strategy.uniquely_determines_coordination_environments
        self.statistics_dict = None
        self._neighbors_by_indices = None

        if coordination_environments is not None:
            self._coordination_environments = coordination_environments
            self._neighbors = neighbors
            self._bva_valences = bva_valences
            self._structure = structure
            if neighbors_by_indices is not None:
                self._neighbors_by_indices = neighbors_by_indices
            if setup_neighbors_by_indices:
                self._setup_neighbors_by_indices()
            return
        self._coordination_environments = []
        self._neighbors = []

        if structure_environments is None:
            raise RuntimeError("coordination_environments and structure_environments are both None")
        self._strategy.set_structure_environments(structure_environments)
        self._structure = structure_environments.structure
        self._bva_valences = structure_environments.bva_valences if bva_valences is None else bva_valences
        for isite, site in enumerate(structure_environments.structure):
            if structure_environments.ce_list[isite] is None:
                self._coordination_environments.append([])
                self._neighbors.append({})
                continue
            site_ce_list = self._strategy.get_site_ce_fractions_and_neighbors(site)
            if site_ce_list is None:
                self._coordination_environments.append([])
                self._neighbors.append({})
                continue
            #This contains a dictionary with keys "ce" and "neighbors".
            #  In site_ce_list["ce"], a list of dictionaries with each possible environment determined by the strategy
            #    In each dictionary, the following keys exist : mp_symbol, fraction, cn_map and csm
            #  In site_ce_list["neighbors"], a dictionary with the corresponding neighbors (as a list of pymatgen's
            #    Site objects) for each cn_map is stored
            ce_piece_dict_list = []
            neighbors_list = {}
            for ce_piece in site_ce_list['ce']:
                if ce_piece['cn_map'] not in neighbors_list:
                    neighbors_list[ce_piece['cn_map']] = site_ce_list['neighbors'][ce_piece['cn_map']]
                ce_piece_dict = {'ce_symbol': ce_piece['mp_symbol'], 'fraction': ce_piece['fraction'],
                                 'csm': ce_piece['csm'], 'cn_map': ce_piece['cn_map']}
                ce_piece_dict_list.append(ce_piece_dict)
            self._coordination_environments.append(ce_piece_dict_list)
            self._neighbors.append(neighbors_list)
        if setup_neighbors_by_indices:
            self._setup_neighbors_by_indices()

    def _setup_neighbors_by_indices(self):
        if not self._strategy.uniquely_determines_coordination_environments:
            raise ValueError('Should use a strategy that uniquely determines coordination environments !')
        self._neighbors_by_indices = []
        for isite_this_site, this_site_neighbors_all_cn_maps in enumerate(self._neighbors):
            if len(self._coordination_environments[isite_this_site]) == 0:
                self._neighbors_by_indices.append(None)
                continue
            cn_map = self._coordination_environments[isite_this_site][0]['cn_map']
            this_site_neighbors = this_site_neighbors_all_cn_maps[tuple(cn_map)]
            if this_site_neighbors is None:
                self._neighbors_by_indices.append(None)
            else:
                this_site_neighbors_indices = []
                for neighbor in this_site_neighbors:
                    ineighbor = None
                    for isite, site in enumerate(self.structure):
                        found = False
                        for tolerance in [1e-8, 1e-6, 1e-4, 1e-3]:
                            if site.is_periodic_image(neighbor, check_lattice=False, tolerance=tolerance):
                                ineighbor = isite
                                diff = neighbor.frac_coords - site.frac_coords
                                image_cell = np.array(np.round(diff), np.int)
                                #if np.allclose(image_cell, np.zeros(3, np.int)):
                                #    image_cell = 0
                                found = True
                                break
                        if found:
                            break
                    if not found:
                        raise ChemenvError('LightStructureEnvironments', '__init__', 'Site indices not found')
                    this_site_neighbors_indices.append((ineighbor, image_cell))
                self._neighbors_by_indices.append(this_site_neighbors_indices)

    def setup_statistic_lists(self):
        self.statistics_dict = {'anion_list': {},                                                    # OK
                                'anion_number': None,                                                # OK
                                'anion_atom_list': {},                                               # OK
                                'anion_atom_number': None,                                           # OK
                                'cation_list': {},                                                   # OK
                                'cation_number': None,                                               # OK
                                'cation_atom_list': {},                                              # OK
                                'cation_atom_number': None,                                          # OK
                                'structure_anion_list': {},                                          # OK
                                'structure_anion_number': None,                                      # OK
                                'structure_anion_undefined_oxidation_list': [],                      # OK
                                'structure_anion_undefined_oxidation_occupation': [],                # OK
                                'structure_anion_atom_list': {},                                     # OK
                                'structure_anion_atom_number': None,                                 # OK
                                'structure_cation_list': {},                                         # OK
                                'structure_cation_number': None,                                     # OK
                                'structure_cation_undefined_oxidation_list': [],                     # OK
                                'structure_cation_undefined_oxidation_occupation': [],               # OK
                                'structure_cation_atom_list': {},                                    # OK
                                'structure_cation_atom_number': None,                                # OK
                                'bva_anion_list': {},                                                # OK
                                'bva_anion_number': None,                                            # OK
                                'bva_anion_atom_list': {},                                           # OK
                                'bva_anion_atom_number': None,                                       # OK
                                'bva_cation_list': {},                                               # OK
                                'bva_cation_number': None,                                           # OK
                                'bva_cation_atom_list': {},                                          # OK
                                'bva_cation_atom_number': None,                                      # OK
                                'neutral_list': {},                                                  # OK
                                'neutral_number': None,                                              # OK
                                'neutral_atom_list': {},                                             # OK
                                'neutral_atom_number': None,                                         # OK
                                'structure_neutral_list': {},                                        # OK
                                'structure_neutral_number': None,                                    # OK
                                'structure_neutral_atom_list': {},                                   # OK
                                'structure_neutral_atom_number': None,                               # OK
                                'bva_neutral_list': {},                                              # OK
                                'bva_neutral_number': None,                                          # OK
                                'bva_neutral_atom_list': {},                                         # OK
                                'bva_neutral_atom_number': None,                                     # OK
                                'atom_coordination_environments_present': {},                        # OK
                                'structure_ion_coordination_environments_present': {},               # OK
                                'bva_ion_coordination_environments_present': {},                     # OK
                                'ion_coordination_environments_present': {},                         # OK
                                'coordination_environments_structure_ion_present': {},               # OK
                                'coordination_environments_ion_present': {},                         # OK
                                'coordination_environments_bva_ion_present': {},                     # OK
                                'coordination_environments_atom_present': {},                        # OK
                                'fraction_bva_ion_coordination_environments_present': {},            # OK
                                'fraction_structure_ion_coordination_environments_present': {},      # OK
                                'fraction_ion_coordination_environments_present': {},                # OK
                                'fraction_atom_coordination_environments_present': {},               # OK
                                'fraction_coordination_environments_bva_ion_present': {},            # OK
                                'fraction_coordination_environments_structure_ion_present': {},      # OK
                                'fraction_coordination_environments_ion_present': {},                # OK
                                'fraction_coordination_environments_atom_present': {},               # OK
                                'count_bva_ion_present': {},                                         # OK
                                'count_structure_ion_present': {},                                   # OK
                                'count_ion_present': {},                                             # OK
                                'count_atom_present': {},                                            # OK
                                'count_coordination_environments_present': {}}
        structure_has_species = isinstance(self._structure[0].species_and_occu.elements[0], Specie)
        atom_stat = self.statistics_dict['atom_coordination_environments_present']
        ce_atom_stat = self.statistics_dict['coordination_environments_atom_present']
        fraction_atom_stat = self.statistics_dict['fraction_atom_coordination_environments_present']
        fraction_ce_atom_stat = self.statistics_dict['fraction_coordination_environments_atom_present']
        count_atoms = self.statistics_dict['count_atom_present']
        count_ce = self.statistics_dict['count_coordination_environments_present']
        for isite, site in enumerate(self._structure):
            #Building anion and cation list
            site_species = {'bva_oxistates': [],
                            'structure_oxistates': []}
            if self._bva_valences != 'undefined':
                for sp, occ in site.species_and_occu.items():
                    valence = self._bva_valences[isite]
                    if valence < 0:
                        specielist = self.statistics_dict['bva_anion_list']
                        atomlist = self.statistics_dict['bva_anion_atom_list']
                    elif valence > 0:
                        specielist = self.statistics_dict['bva_cation_list']
                        atomlist = self.statistics_dict['bva_cation_atom_list']
                    else:
                        specielist = self.statistics_dict['bva_neutral_list']
                        atomlist = self.statistics_dict['bva_neutral_atom_list']
                    strspecie = str(Specie(sp.symbol, valence))
                    if strspecie not in specielist:
                        specielist[strspecie] = occ
                    else:
                        specielist[strspecie] += occ
                    if sp.symbol not in atomlist:
                        atomlist[sp.symbol] = occ
                    else:
                        atomlist[sp.symbol] += occ
                    site_species['bva_oxistates'].append((sp.symbol, valence, occ))
                self.statistics_dict['bva_anion_number'] = len(self.statistics_dict['bva_anion_list'])
                self.statistics_dict['bva_anion_atom_number'] = len(self.statistics_dict['bva_anion_atom_list'])
                self.statistics_dict['bva_cation_number'] = len(self.statistics_dict['bva_cation_list'])
                self.statistics_dict['bva_cation_atom_number'] = len(self.statistics_dict['bva_cation_atom_list'])
            if structure_has_species:
                for sp, occ in site.species_and_occu.items():
                    oxi_state_rounded = round(sp.oxi_state)
                    if np.isclose(sp.oxi_state, oxi_state_rounded, rtol=0.0, atol=self.DELTA_MAX_OXIDATION_STATE):
                        valence = int(oxi_state_rounded)
                        if valence < 0:
                            specielist = self.statistics_dict['structure_anion_list']
                            atomlist = self.statistics_dict['structure_anion_atom_list']
                        elif valence > 0:
                            specielist = self.statistics_dict['structure_cation_list']
                            atomlist = self.statistics_dict['structure_cation_atom_list']
                        else:
                            specielist = self.statistics_dict['structure_neutral_list']
                            atomlist = self.statistics_dict['structure_neutral_atom_list']
                        strspecie = str(Specie(sp.symbol, valence))
                        if strspecie not in specielist:
                            specielist[strspecie] = occ
                        else:
                            specielist[strspecie] += occ
                        if sp.symbol not in atomlist:
                            atomlist[sp.symbol] = occ
                        else:
                            atomlist[sp.symbol] += occ
                        site_species['structure_oxistates'].append((sp.symbol, valence, occ))
                    else:
                        valence = sp.oxi_state
                        if valence < 0:
                            specielist = self.statistics_dict['structure_anion_undefined_oxidation_list']
                            specieocc = self.statistics_dict['structure_anion_undefined_oxidation_occupation']
                            atomlist = self.statistics_dict['structure_anion_atom_list']
                        elif valence > 0:
                            specielist = self.statistics_dict['structure_cation_undefined_oxidation_list']
                            specieocc = self.statistics_dict['structure_cation_undefined_oxidation_occupation']
                            atomlist = self.statistics_dict['structure_cation_atom_list']
                        else:
                            raise RuntimeError('Should not be here ...')
                        strspecie = str(Specie(sp.symbol, valence))
                        if strspecie not in specielist:
                            specielist.append(strspecie)
                            specieocc.append(occ)
                        else:
                            specie_index = specielist.index(strspecie)
                            specieocc[specie_index] += occ
                        if sp.symbol not in atomlist:
                            atomlist[sp.symbol] = occ
                        else:
                            atomlist[sp.symbol] += occ
                self.statistics_dict['structure_anion_number'] = len(self.statistics_dict['structure_anion_list']) + len(self.statistics_dict['structure_anion_undefined_oxidation_list'])
                self.statistics_dict['structure_anion_atom_number'] = len(self.statistics_dict['structure_anion_atom_list'])
                self.statistics_dict['structure_cation_number'] = len(self.statistics_dict['structure_cation_list']) + len(self.statistics_dict['structure_cation_undefined_oxidation_list'])
                self.statistics_dict['structure_cation_atom_number'] = len(self.statistics_dict['structure_cation_atom_list'])
            # for prefix in ['structure_', 'bva_']:
            #     if ((prefix == 'bva_' and self._bva_valences != 'undefined') or
            #             (prefix == 'structure_' and structure_has_species)):
            #         for sp, occ in site.species_and_occu.items():
            #             valence = self._bva_valences[isite] if prefix == 'bva_' else sp.oxi_state
            #             if valence < 0:
            #                 specielist = self.statistics_dict['{}anion_list'.format(prefix)]
            #                 atomlist = self.statistics_dict['{}anion_atom_list'.format(prefix)]
            #             elif valence > 0:
            #                 specielist = self.statistics_dict['{}cation_list'.format(prefix)]
            #                 atomlist = self.statistics_dict['{}cation_atom_list'.format(prefix)]
            #             else:
            #                 specielist = self.statistics_dict['{}neutral_list'.format(prefix)]
            #                 atomlist = self.statistics_dict['{}neutral_atom_list'.format(prefix)]
            #             strspecie = str(Specie(sp.symbol, valence))
            #             if strspecie not in specielist:
            #                 specielist[strspecie] = occ
            #             else:
            #                 specielist[strspecie] += occ
            #             if sp.symbol not in atomlist:
            #                 atomlist[sp.symbol] = occ
            #             else:
            #                 atomlist[sp.symbol] += occ
            #             site_species['{}oxistates'.format(prefix)].append((sp.symbol, valence, occ))
            #         self.statistics_dict['{}anion_number'.format(prefix)] = len(self.statistics_dict['{}anion_list'.format(prefix)])
            #         self.statistics_dict['{}anion_atom_number'.format(prefix)] = len(self.statistics_dict['{}anion_atom_list'.format(prefix)])
            #         self.statistics_dict['{}cation_number'.format(prefix)] = len(self.statistics_dict['{}cation_list'.format(prefix)])
            #         self.statistics_dict['{}cation_atom_number'.format(prefix)] = len(self.statistics_dict['{}cation_atom_list'.format(prefix)])

            #Building environments lists
            if self._coordination_environments[isite] is not None:
                site_envs = [(ce_piece_dict['ce_symbol'], ce_piece_dict['fraction'])
                             for ce_piece_dict in self._coordination_environments[isite]]
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

                if self._bva_valences != 'undefined':
                    ion_stat = self.statistics_dict['bva_ion_coordination_environments_present']
                    ce_ion_stat = self.statistics_dict['coordination_environments_bva_ion_present']
                    count_ions = self.statistics_dict['count_bva_ion_present']
                    for elmt, oxi_state, occ in site_species['bva_oxistates']:
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
                if structure_has_species:
                    ion_stat = self.statistics_dict['structure_ion_coordination_environments_present']
                    ce_ion_stat = self.statistics_dict['coordination_environments_structure_ion_present']
                    count_ions = self.statistics_dict['count_structure_ion_present']
                    for elmt, oxi_state, occ in site_species['structure_oxistates']:
                        if elmt not in ion_stat:
                            ion_stat[elmt] = {}
                            count_ions[elmt] = {}
                        oxi_state_rounded = round(oxi_state)
                        if np.isclose(oxi_state, oxi_state_rounded, rtol=0.0, atol=self.DELTA_MAX_OXIDATION_STATE):
                            my_oxi_state = int(oxi_state_rounded)
                            if my_oxi_state not in ion_stat[elmt]:
                                ion_stat[elmt][my_oxi_state] = {}
                                count_ions[elmt][my_oxi_state] = 0.0
                            count_ions[elmt][my_oxi_state] += occ
                            for ce_symbol, fraction in site_envs:
                                if fraction is None:
                                    continue
                                if ce_symbol not in ion_stat[elmt][my_oxi_state]:
                                    ion_stat[elmt][my_oxi_state][ce_symbol] = 0.0
                                ion_stat[elmt][my_oxi_state][ce_symbol] += occ * fraction
                                if ce_symbol not in ce_ion_stat:
                                    ce_ion_stat[ce_symbol] = {}
                                if elmt not in ce_ion_stat[ce_symbol]:
                                    ce_ion_stat[ce_symbol][elmt] = {}
                                if my_oxi_state not in ce_ion_stat[ce_symbol][elmt]:
                                    ce_ion_stat[ce_symbol][elmt][my_oxi_state] = 0.0
                                ce_ion_stat[ce_symbol][elmt][my_oxi_state] += occ * fraction
                        else:
                            #TODO: make something slightly better here ... like put the fractional oxidation states
                            #somewhere ...
                            continue
                # for prefix in ['structure_', 'bva_']:
                #     if prefix == 'bva_' and self._bva_valences == 'undefined':
                #         continue
                #     elif prefix == 'structure_' and not structure_has_species:
                #         continue
                #     ion_stat = self.statistics_dict['{}ion_coordination_environments_present'.format(prefix)]
                #     ce_ion_stat = self.statistics_dict['coordination_environments_{}ion_present'.format(prefix)]
                #     count_ions = self.statistics_dict['count_{}ion_present'.format(prefix)]
                #     for elmt, oxi_state, occ in site_species['{}oxistates'.format(prefix)]:
                #         if elmt not in ion_stat:
                #             ion_stat[elmt] = {}
                #             count_ions[elmt] = {}
                #         if oxi_state not in ion_stat[elmt]:
                #             ion_stat[elmt][oxi_state] = {}
                #             count_ions[elmt][oxi_state] = 0.0
                #         count_ions[elmt][oxi_state] += occ
                #         for ce_symbol, fraction in site_envs:
                #             if fraction is None:
                #                 continue
                #             if ce_symbol not in ion_stat[elmt][oxi_state]:
                #                 ion_stat[elmt][oxi_state][ce_symbol] = 0.0
                #             ion_stat[elmt][oxi_state][ce_symbol] += occ * fraction
                #             if ce_symbol not in ce_ion_stat:
                #                 ce_ion_stat[ce_symbol] = {}
                #             if elmt not in ce_ion_stat[ce_symbol]:
                #                 ce_ion_stat[ce_symbol][elmt] = {}
                #             if oxi_state not in ce_ion_stat[ce_symbol][elmt]:
                #                 ce_ion_stat[ce_symbol][elmt][oxi_state] = 0.0
                #             ce_ion_stat[ce_symbol][elmt][oxi_state] += occ * fraction
        for elmt, envs in atom_stat.items():
            sumelement = count_atoms[elmt]
            fraction_atom_stat[elmt] = {env: fraction / sumelement for env, fraction in envs.items()}
        for ce_symbol, atoms in ce_atom_stat.items():
            sumsymbol = count_ce[ce_symbol]
            fraction_ce_atom_stat[ce_symbol] = {atom: fraction / sumsymbol for atom, fraction in atoms.items()}
        for prefix in ['structure_', 'bva_']:
            ion_stat = self.statistics_dict['{}ion_coordination_environments_present'.format(prefix)]
            fraction_ion_stat = self.statistics_dict['fraction_{}ion_coordination_environments_present'.format(prefix)]
            ce_ion_stat = self.statistics_dict['coordination_environments_{}ion_present'.format(prefix)]
            fraction_ce_ion_stat = self.statistics_dict['fraction_coordination_environments_{}ion_present'.format(prefix)]
            count_ions = self.statistics_dict['count_{}ion_present'.format(prefix)]
            if prefix == 'bva_' and self._bva_valences == 'undefined':
                continue
            elif prefix == 'structure_' and not structure_has_species:
                continue
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
        ceip = 'coordination_environments_ion_present'
        icep = 'ion_coordination_environments_present'
        fceip = 'fraction_coordination_environments_ion_present'
        ficep = 'fraction_ion_coordination_environments_present'
        if self._bva_valences != 'undefined':
            cebip = 'coordination_environments_bva_ion_present'
            bicep = 'bva_ion_coordination_environments_present'
            fcebip = 'fraction_coordination_environments_bva_ion_present'
            fbicep = 'fraction_bva_ion_coordination_environments_present'
            self.statistics_dict['anion_list'].update(self.statistics_dict['bva_anion_list'])
            self.statistics_dict['anion_atom_list'].update(self.statistics_dict['bva_anion_atom_list'])
            self.statistics_dict['cation_list'].update(self.statistics_dict['bva_cation_list'])
            self.statistics_dict['cation_atom_list'].update(self.statistics_dict['bva_cation_atom_list'])
            self.statistics_dict['anion_number'] = self.statistics_dict['bva_anion_number']
            self.statistics_dict['anion_atom_number'] = self.statistics_dict['bva_anion_atom_number']
            self.statistics_dict['cation_number'] = self.statistics_dict['bva_cation_number']
            self.statistics_dict['cation_atom_number'] = self.statistics_dict['bva_cation_atom_number']
            self.statistics_dict[ceip].update(self.statistics_dict[cebip])
            self.statistics_dict[icep].update(self.statistics_dict[bicep])
            self.statistics_dict[fceip].update(self.statistics_dict[fcebip])
            self.statistics_dict[ficep].update(self.statistics_dict[fbicep])
        else:
            cesip = 'coordination_environments_structure_ion_present'
            sicep = 'structure_ion_coordination_environments_present'
            fcesip = 'fraction_coordination_environments_structure_ion_present'
            fsicep = 'fraction_structure_ion_coordination_environments_present'
            self.statistics_dict['anion_list'].update(self.statistics_dict['structure_anion_list'])
            self.statistics_dict['anion_atom_list'].update(self.statistics_dict['structure_anion_atom_list'])
            self.statistics_dict['cation_list'].update(self.statistics_dict['structure_cation_list'])
            self.statistics_dict['cation_atom_list'].update(self.statistics_dict['structure_cation_atom_list'])
            self.statistics_dict['anion_number'] = self.statistics_dict['structure_anion_number']
            self.statistics_dict['anion_atom_number'] = self.statistics_dict['structure_anion_atom_number']
            self.statistics_dict['cation_number'] = self.statistics_dict['structure_cation_number']
            self.statistics_dict['cation_atom_number'] = self.statistics_dict['structure_cation_atom_number']
            self.statistics_dict[ceip].update(self.statistics_dict[cesip])
            self.statistics_dict[icep].update(self.statistics_dict[sicep])
            self.statistics_dict[fceip].update(self.statistics_dict[fcesip])
            self.statistics_dict[ficep].update(self.statistics_dict[fsicep])

    def site_has_clear_environment(self, isite, min_fraction=0.95):
        if len(self._coordination_environments[isite]) == 0:
            return True
        fractions = [ce_dict['fraction'] for ce_dict in self._coordination_environments[isite]]
        maxfraction = max(fractions)
        return maxfraction > min_fraction

    def structure_has_clear_environments_for_element(self, element, min_fraction=0.95):
        for isite, site_ces in enumerate(self._coordination_environments):
            if site_ces is None:
                continue
            if element in [sp.symbol for sp in self._structure[isite].species_and_occu]:
                if not self.site_has_clear_environment(isite=isite, min_fraction=min_fraction):
                    return False
        return True

    def structure_has_clear_environments(self, min_fraction=0.95):
        for isite, site_ces in enumerate(self._coordination_environments):
            if site_ces is None:
                continue
            if not self.site_has_clear_environment(isite=isite, min_fraction=min_fraction):
                return False
        return True

    def get_site_info_for_specie_ce(self, specie, ce_symbol, min_fraction=0.0):
        element = specie.symbol
        oxi_state = specie.oxi_state
        isites = []
        csms = []
        fractions = []
        for isite, site in enumerate(self._structure):
            if element in [sp.symbol for sp in site.species_and_occu]:
                if oxi_state == self._bva_valences[isite]:
                    for ce_dict in self._coordination_environments[isite]:
                        if ce_symbol == ce_dict['ce_symbol']:
                            isites.append(isite)
                            csms.append(ce_dict['ce_symbol'])
                            fractions.append(ce_dict['fraction'])

    def get_site_info_for_specie_allces(self, specie, min_fraction=0.0):
        allces = {}
        element = specie.symbol
        oxi_state = specie.oxi_state
        for isite, site in enumerate(self._structure):
            if element in [sp.symbol for sp in site.species_and_occu]:
                if oxi_state == self._bva_valences[isite]:
                    for ce_dict in self._coordination_environments[isite]:
                        if ce_dict['fraction'] < min_fraction:
                            continue
                        if ce_dict['ce_symbol'] not in allces:
                            allces[ce_dict['ce_symbol']] = {'isites': [], 'fractions': [], 'csms': []}
                        allces[ce_dict['ce_symbol']]['isites'].append(isite)
                        allces[ce_dict['ce_symbol']]['fractions'].append(ce_dict['fraction'])
                        allces[ce_dict['ce_symbol']]['csms'].append(ce_dict['csm'])
        return allces

    def get_statistics(self, statistics_fields=DEFAULT_STATISTICS_FIELDS, bson_compatible=False):
        if self.statistics_dict is None:
            self.setup_statistic_lists()
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
        return ce_symbol in [ce_dict['ce_symbol'] for ce_dict in self._coordination_environments[isite]]

    def structure_contains_atom_environment(self, atom_symbol, ce_symbol):
        """
        Checks whether the structure contains a given atom in a given environment
        :param atom_symbol: Symbol of the atom
        :param ce_symbol: Symbol of the coordination environment
        :return: True if the coordination environment is found, False otherwise
        """
        for isite, site in enumerate(self._structure):
            if (Element(atom_symbol) in site.species_and_occu.
                    element_composition and self.site_contains_environment(isite, ce_symbol)):
                return True
        return False

    @property
    def coordination_environments(self):
        """
        Coordination environments determined by the strategy
        """
        return self._coordination_environments

    @property
    def uniquely_determined_coordination_environments(self):
        """
        True if the coordination environments are uniquely determined.
        """
        return self._uniquely_determined_coordination_environments

    @property
    def neighbors(self):
        """
        Neighbors determined by the strategy
        """
        return self._neighbors

    @property
    def neighbors_by_indices(self):
        """
        Neighbors by indices and image cell determined by the strategy
        """
        if self._neighbors_by_indices is None:
            self._setup_neighbors_by_indices()
        return self._neighbors_by_indices

    @property
    def structure(self):
        """
        Structure
        :return: Structure
        """
        return self._structure

    def __eq__(self, other):
        """
        Equality method that checks if the LightStructureEnvironments object is equal to another
        LightStructureEnvironments object. Two LightStructureEnvironments objects are equal if the strategy used
        is the same, if the structure is the same, if the valences used in the strategies are the same, if the
        coordination environments and the neighbours determined by the strategy are the same
        :param other: LightStructureEnvironments object to compare with
        :return: True if both objects are equal, False otherwise
        """
        return (self._strategy == other._strategy and self._structure == other._structure and
                self._bva_valences == other._bva_valences and
                self._coordination_environments == other._coordination_environments and
                self._neighbors == other._neighbors)

    def as_dict(self):
        """
        Bson-serializable dict representation of the LightStructureEnvironments object.
        :return: Bson-serializable dict representation of the LightStructureEnvironments object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "strategy": self._strategy.as_dict(),
                "structure": self._structure.as_dict(),
                "bva_valences": self._bva_valences,
                "coordination_environments": self._coordination_environments,
                "neighbors": [{'{:d}_{:d}'.format(cn_map[0], cn_map[1]): [ps.as_dict() for ps in neighbors]
                               for cn_map, neighbors in self._neighbors[isite].items()}
                              for isite in range(len(self._structure))],
                "neighbors_by_indices": jsanitize([None if item is None else list(item)
                                                   for item in self.neighbors_by_indices])}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the LightStructureEnvironments object from a dict representation of the
        LightStructureEnvironments created using the as_dict method.
        :param d: dict representation of the LightStructureEnvironments object
        :return: LightStructureEnvironments object
        """
        dec = MontyDecoder()
        neighbors = [{(int(cnmap.split('_')[0]), int(cnmap.split('_')[1])): dec.process_decoded(cnmap_neighbours)
                      for cnmap, cnmap_neighbours in neighbs_site.items()} for neighbs_site in d['neighbors']]
        coordination_environments = [[{key: val if key != 'cn_map' else tuple(val) for key, val in item.items()}
                                      for item in ces_site] for ces_site in d['coordination_environments']]
        if 'neighbors_by_indices' in d:
            neighbors_by_indices = []
            for item in d['neighbors_by_indices']:
                if item is None:
                    neighbors_by_indices.append(None)
                else:
                    neighbors_by_indices.append([(nb[0], np.array(nb[1], np.int)) for nb in item])
            # neighbors_by_indices = [None if item is None else (item[0], np.array(item[1], np.int))
            #                         for item in d['neighbors_by_indices']]
        else:
            neighbors_by_indices = None
        return cls(dec.process_decoded(d['strategy']), structure_environments=None,
                   structure=dec.process_decoded(d['structure']),
                   bva_valences=d['bva_valences'],
                   coordination_environments=coordination_environments,
                   neighbors=neighbors,
                   neighbors_by_indices=neighbors_by_indices)


# class LightStructureEnvironmentsOld(MSONable):
#     """
#     Class used to store the chemical environments of a given structure obtained from a given ChemenvStrategy. Currently,
#     only strategies leading to the determination of a unique environment for each site is allowed
#     This class does not store all the information contained in the StructureEnvironments object, only the coordination
#     environment found
#     """
#
#     def __init__(self, strategy_used, structure_environments=None, structure=None, bva_valences=None,
#                  coordination_environments=None, neighbors=None):
#         """
#         Constructor for the LightStructureEnvironments object.
#         """
#         if not strategy_used.uniquely_determines_coordination_environments:
#             raise NotImplementedError('LightStructureEnvironments not yet implemented with complex strategies leading '
#                                       'to multiple environments of a same site ...')
#         self.strategy_used = strategy_used
#         self._uniquely_determined_coordination_environments = strategy_used.uniquely_determines_coordination_environments
#         self._neighbors = []
#         self._coordination_environments = []
#         if structure_environments is not None:
#             self.strategy_used.set_structure_environments(structure_environments)
#             self.structure = structure_environments.structure
#             self.bva_valences = structure_environments.bva_valences
#             for site in self.strategy_used.structure_environments.structure:
#                 try:
#                     this_site_neighbors = self.strategy_used.get_site_neighbors(site)
#                     self._neighbors.append(this_site_neighbors)
#                     try:
#                         (ce, ce_dict) = self.strategy_used.get_site_coordination_environment(site)
#                         self._coordination_environments.append(ce)
#                     except TypeError:
#                         cn = self.strategy_used.get_site_coordination_environment(site)
#                         self._coordination_environments.append(cn)
#                 except NeighborsNotComputedChemenvError:
#                     self._neighbors.append(None)
#                     self._coordination_environments.append(None)
#         else:
#             if structure is None or coordination_environments is None or neighbors is None or bva_valences is None:
#                 raise InitializationChemenvError(self.__class__.__name__)
#             self.structure = structure
#             self._coordination_environments = coordination_environments
#             self._neighbors = neighbors
#             self.bva_valences = bva_valences
#         self._neighbors_by_indices = []
#         for this_site_neighbors in self._neighbors:
#             if this_site_neighbors is None:
#                 self._neighbors_by_indices.append(None)
#             else:
#                 this_site_neighbors_indices = []
#                 for neighbor in this_site_neighbors:
#                     ineighbor = None
#                     for isite, site in enumerate(self.structure):
#                         found = False
#                         for tolerance in [1e-8, 1e-6, 1e-4, 1e-3]:
#                             if site.is_periodic_image(neighbor, check_lattice=False, tolerance=tolerance):
#                                 ineighbor = isite
#                                 diff = neighbor.frac_coords - site.frac_coords
#                                 image_cell = np.array(np.round(diff), np.int)
#                                 #if np.allclose(image_cell, np.zeros(3, np.int)):
#                                 #    image_cell = 0
#                                 found = True
#                                 break
#                         if found:
#                             break
#                     if not found:
#                         raise ChemenvError('LightStructureEnvironments', '__init__', 'Site indices not found')
#                     this_site_neighbors_indices.append((ineighbor, image_cell))
#                 self._neighbors_by_indices.append(this_site_neighbors_indices)
#         self.setup_statistic_lists()
#
#     def contains_only_one_anion_atom(self, anion_atom):
#         return len(self.anion_atom_list) == 1 and self.anion_atom_list[0] == anion_atom
#
#     def contains_only_one_anion(self, anion):
#         return len(self.anion_list) == 1 and self.anion_list[0] == anion
#
#     def site_contains_environment(self, isite, mp_symbol):
#         if self.uniquely_determined_coordination_environments:
#             return self.coordination_environments[isite] == mp_symbol
#         else:
#             return mp_symbol in self.coordination_environments[isite]
#
#     def structure_contains_atom_environment(self, atom_symbol, mp_symbol):
#         """
#         Checks whether the structure contains a given atom in a given environment
#         :param atom_symbol: Symbol of the atom
#         :param mp_symbol: Symbol of the coordination environment
#         :return: True if the coordination environment is found, False otherwise
#         """
#         if self.uniquely_determined_coordination_environments:
#             for isite, site in enumerate(self.structure):
#                 if Element(atom_symbol) in site.species_and_occu.element_composition:
#                     if self.coordination_environments[isite] == mp_symbol:
#                         return True
#             return False
#         else:
#             for isite, site in enumerate(self.structure):
#                 if Element(atom_symbol) in site.species_and_occu.element_composition:
#                     if mp_symbol in self.coordination_environments[isite]:
#                         return True
#             return False
#
#     def setup_statistic_lists(self):
#         self.anion_list = []
#         self.anion_atom_list = []
#         self.cation_list = []
#         self.cation_atom_list = []
#         self.coordination_environments_present = {}
#         self.site_count_with_computed_ce = 0.0
#         self.atom_coordination_environments_present = {}
#         self.atom_count_with_computed_ce = {}
#         self.structure_ion_coordination_environments_present = {}
#         self.structure_ion_count_with_computed_ce = {}
#         self.bva_ion_coordination_environments_present = {}
#         self.bva_ion_count_with_computed_ce = {}
#         self.coordination_environments_bva_ion_present = {}
#         self.fraction_bva_ion_coordination_environments_present = {}
#         self.fraction_coordination_environments_bva_ion_present = {}
#         for isite, site in enumerate(self.structure):
#             #Building anion list
#             if self.bva_valences != 'undefined':
#                 for sp, occ in site.species_and_occu.items():
#                     if self.bva_valences[isite] < 0:
#                         stranion = str(Specie(sp.symbol, self.bva_valences[isite]))
#                         if stranion not in self.anion_list:
#                             self.anion_list.append(stranion)
#                         if sp.symbol not in self.anion_atom_list:
#                             self.anion_atom_list.append(sp.symbol)
#                     else:
#                         strcation = str(Specie(sp.symbol, self.bva_valences[isite]))
#                         if strcation not in self.cation_list:
#                             self.cation_list.append(strcation)
#                         if sp.symbol not in self.cation_atom_list:
#                             self.cation_atom_list.append(sp.symbol)
#             if self.coordination_environments[isite] is not None:
#                 ce = '{}'.format(self.coordination_environments[isite])
#                 self.site_count_with_computed_ce += 1.0
#                 if not ce in self.coordination_environments_present:
#                     self.coordination_environments_present[ce] = 1
#                 else:
#                     self.coordination_environments_present[ce] += 1
#                 for specie in site.species_and_occu:
#                     #Counting atoms in specific environments
#                     if specie.symbol not in self.atom_coordination_environments_present:
#                         self.atom_coordination_environments_present[specie.symbol] = {}
#                     atom_ce = '{}'.format(self.coordination_environments[isite])
#                     if not atom_ce in self.atom_coordination_environments_present:
#                         self.atom_coordination_environments_present[specie.symbol][atom_ce] = site.species_and_occu[specie]
#                     else:
#                         self.atom_coordination_environments_present[specie.symbol][atom_ce] += site.species_and_occu[specie]
#                     #Counting atoms for which environments are determined
#                     if specie.symbol not in self.atom_count_with_computed_ce:
#                         self.atom_count_with_computed_ce[specie.symbol] = site.species_and_occu[specie]
#                     else:
#                         self.atom_count_with_computed_ce[specie.symbol] += site.species_and_occu[specie]
#                     #Counting ions in specific environments
#                     if specie.symbol not in self.structure_ion_coordination_environments_present:
#                         self.structure_ion_coordination_environments_present[specie.symbol] = {}
#                     try:
#                         str_sp_oxi_state = str(specie.oxi_state)
#                     except AttributeError:
#                         str_sp_oxi_state = 'element'
#                     str_sp_oxi_state = str_sp_oxi_state.replace('.', ',')
#                     if str_sp_oxi_state not in self.structure_ion_coordination_environments_present[specie.symbol]:
#                         self.structure_ion_coordination_environments_present[specie.symbol][str_sp_oxi_state] = {}
#                     ion_ce = '{}'.format(self.coordination_environments[isite])
#                     if ion_ce not in self.structure_ion_coordination_environments_present[specie.symbol][str_sp_oxi_state]:
#                         self.structure_ion_coordination_environments_present[specie.symbol][str_sp_oxi_state][ion_ce] = site.species_and_occu[specie]
#                     else:
#                         self.structure_ion_coordination_environments_present[specie.symbol][str_sp_oxi_state][ion_ce] += site.species_and_occu[specie]
#                     #Counting ions for which environments are determined
#                     if specie.symbol not in self.structure_ion_count_with_computed_ce:
#                         self.structure_ion_count_with_computed_ce[specie.symbol] = {}
#                     if str_sp_oxi_state not in self.structure_ion_count_with_computed_ce[specie.symbol]:
#                         self.structure_ion_count_with_computed_ce[specie.symbol][str_sp_oxi_state] = site.species_and_occu[specie]
#                     else:
#                         self.structure_ion_count_with_computed_ce[specie.symbol][str_sp_oxi_state] += site.species_and_occu[specie]
#                         #Counting ions in specific environments
#                     if specie.symbol not in self.bva_ion_coordination_environments_present:
#                         self.bva_ion_coordination_environments_present[specie.symbol] = {}
#                     if self.bva_valences != 'undefined':
#                         str_val = str(self.bva_valences[isite])
#                         if str_val not in self.bva_ion_coordination_environments_present[specie.symbol]:
#                             self.bva_ion_coordination_environments_present[specie.symbol][str_val] = {}
#                         ion_ce = '{}'.format(self.coordination_environments[isite])
#                         if ion_ce not in self.bva_ion_coordination_environments_present[specie.symbol][str_val]:
#                             self.bva_ion_coordination_environments_present[specie.symbol][str_val][ion_ce] = site.species_and_occu[specie]
#                         else:
#                             self.bva_ion_coordination_environments_present[specie.symbol][str_val][ion_ce] += site.species_and_occu[specie]
#                         if specie.symbol not in self.bva_ion_count_with_computed_ce:
#                             self.bva_ion_count_with_computed_ce[specie.symbol] = {}
#                         if str_val not in self.bva_ion_count_with_computed_ce[specie.symbol]:
#                             self.bva_ion_count_with_computed_ce[specie.symbol][str_val] = site.species_and_occu[specie]
#                         else:
#                             self.bva_ion_count_with_computed_ce[specie.symbol][str_val] += site.species_and_occu[specie]
#         for ce in self.coordination_environments_present:
#             self.coordination_environments_bva_ion_present[ce] = {}
#         for sp, spdict in self.bva_ion_coordination_environments_present.items():
#             self.fraction_bva_ion_coordination_environments_present[sp] = {}
#             for val, valdict in spdict.items():
#                 strbvasp = str(Specie(sp, int(val)))
#                 self.fraction_bva_ion_coordination_environments_present[sp][val] = {}
#                 for ion_ce, n_ion_ce in valdict.items():
#                     self.fraction_bva_ion_coordination_environments_present[sp][val][ion_ce] = float(n_ion_ce) / self.bva_ion_count_with_computed_ce[sp][val]
#                     if strbvasp in self.coordination_environments_bva_ion_present[ion_ce]:
#                         self.coordination_environments_bva_ion_present[ion_ce][strbvasp] += n_ion_ce
#                     else:
#                         self.coordination_environments_bva_ion_present[ion_ce][strbvasp] = n_ion_ce
#         for ce, cedict in self.coordination_environments_bva_ion_present.items():
#             self.fraction_coordination_environments_bva_ion_present[ce] = {}
#             for cation, ncations in cedict.items():
#                 self.fraction_coordination_environments_bva_ion_present[ce][cation] = float(ncations) / self.coordination_environments_present[ce]
#
#     @property
#     def coordination_environments(self):
#         return self._coordination_environments
#
#     @property
#     def uniquely_determined_coordination_environments(self):
#         return self._uniquely_determined_coordination_environments
#
#     @property
#     def neighbors(self):
#         return self._neighbors
#
#     @property
#     def neighbors_by_indices(self):
#         return self._neighbors_by_indices
#
#     def __eq__(self, other):
#         return (self.strategy_used == other.strategy_used and self.structure == other.structure and
#                 self.bva_valences == other.bva_valences and
#                 self._coordination_environments == other._coordination_environments and
#                 self._neighbors == other._neighbors)
#
#     def as_dict(self):
#         """
#         Bson-serializable dict representation of the LightStructureEnvironments object.
#         :return: Bson-serializable dict representation of the LightStructureEnvironments object.
#         """
#         #strategy_used, structure_environments=None, structure=None, coordination_environments=None,
#         #         neighbors=None
#         return {"@module": self.__class__.__module__,
#                 "@class": self.__class__.__name__,
#                 "strategy_used": self.strategy_used.as_dict(),
#                 "structure": self.structure.as_dict(),
#                 "bva_valences": self.bva_valences,
#                 "coordination_environments": self._coordination_environments,
#                 "neighbors": [[ps.as_dict() for ps in neighbs] if neighbs is not None else None
#                               for neighbs in self._neighbors]}
#
#     @classmethod
#     def from_dict(cls, d):
#         """
#         Reconstructs the LightStructureEnvironments object from a dict representation of the
#         LightStructureEnvironments created using the as_dict method.
#         :param d: dict representation of the LightStructureEnvironments object
#         :return: LightStructureEnvironments object
#         """
#         dec = MontyDecoder()
#         return cls(dec.process_decoded(d['strategy_used']), structure_environments=None,
#                    structure=dec.process_decoded(d['structure']),
#                    bva_valences=d['bva_valences'],
#                    coordination_environments=d['coordination_environments'],
#                    neighbors=dec.process_decoded(d['neighbors']))


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

    def minimum_csm(self):
        """
        Returns the minimum continuous symmetry measure of this ChemicalEnvironments object
        :return: The minimum CSM for this ChemicalEnvironments object
        """
        if len(self.coord_geoms) == 0:
            return None
        csms = [self.coord_geoms[cg]['symmetry_measure'] for cg in self.coord_geoms]
        return min(csms)

    def minimum_geometry(self):
        """
        Returns the geometry with the minimum continuous symmetry measure of this ChemicalEnvironments
        :return: tuple (symbol, csm) with symbol being the geometry with the minimum continuous symmetry measure and
        csm being the continuous symmetry measure associted to it
        :raise: ValueError if no coordination geometry is found in this ChemicalEnvironments object
        """
        if len(self.coord_geoms) == 0:
            return None
        cglist = [cg for cg in self.coord_geoms]
        csms = np.array([self.coord_geoms[cg]['symmetry_measure'] for cg in cglist])
        csmlist = [self.coord_geoms[cg] for cg in cglist]
        imin = np.argmin(csms)
        return cglist[imin], csmlist[imin]

    def minimum_geometries(self, n=None):
        """
        Returns a list of geometries with increasing continuous symmetry measure in this ChemicalEnvironments object
        :param n: Number of geometries to be included in the list
        :return: list of geometries with increasing continuous symmetry measure in this ChemicalEnvironments object
        :raise: ValueError if no coordination geometry is found in this ChemicalEnvironments object
        """
        cglist = [cg for cg in self.coord_geoms]
        csms = np.array([self.coord_geoms[cg]['symmetry_measure'] for cg in cglist])
        csmlist = [self.coord_geoms[cg] for cg in cglist]
        isorted = np.argsort(csms)
        if n is None:
            return [(cglist[ii], csmlist[ii]) for ii in isorted]
        else:
            return [(cglist[ii], csmlist[ii]) for ii in isorted[:n]]

    def add_coord_geom(self, mp_symbol, symmetry_measure, algo='UNKNOWN', permutation=None, override=False,
                       local2perfect_map=None, perfect2local_map=None, detailed_voronoi_index=None,
                       other_symmetry_measures=None):
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
                                           'other_symmetry_measures': other_symmetry_measures}

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
        for mp_symbol in self.coord_geoms:
            out += '   - {}\n'.format(mp_symbol)
            out += '      csm : {}'.format(self.coord_geoms[mp_symbol]['symmetry_measure'])
            out += '     algo : {}'.format(self.coord_geoms[mp_symbol]['algo'])
            out += '     perm : {}\n'.format(self.coord_geoms[mp_symbol]['permutation'])
            out += '       local2perfect : {}\n'.format(str(self.coord_geoms[mp_symbol]['local2perfect_map']))
            out += '       perfect2local : {}\n'.format(str(self.coord_geoms[mp_symbol]['perfect2local_map']))
        return out

    def __eq__(self, other):
        """
        Equality method that checks if the ChemicalEnvironments object is equal to another ChemicalEnvironments
        object.
        :param other: ChemicalEnvironments object to compare with
        :return: True if both objects are equal, False otherwise
        """
        return self.coord_geoms == other.coord_geoms

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
                              other_symmetry_measures=other_csms)
        return ce
