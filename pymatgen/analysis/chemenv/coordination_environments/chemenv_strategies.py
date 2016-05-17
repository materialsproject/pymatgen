# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module provides so-called "strategies" to determine the coordination environments of an atom in a structure.
Some strategies can favour larger or smaller environments. Some strategies uniquely identifies the environments while
some others can identify the environment as a "mix" of several environments, each of which is assigned with a given
fraction. The choice of the strategy depends on the purpose of the user.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


import abc
import os
import json
from monty.json import MSONable
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp
from pymatgen.core.sites import PeriodicSite
import numpy as np
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import UNCLEAR_ENVIRONMENT_SYMBOL
from pymatgen.analysis.chemenv.utils.func_utils import CSMFiniteRatioFunction
from pymatgen.analysis.chemenv.utils.func_utils import CSMInfiniteRatioFunction
from pymatgen.analysis.chemenv.utils.func_utils import DeltaCSMRatioFunction
from pymatgen.analysis.chemenv.utils.chemenv_errors import EquivalentSiteSearchError
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.utils.defs_utils import AdditionalConditions
from six import with_metaclass
from pymatgen.analysis.chemenv.coordination_environments.voronoi import DetailedVoronoiContainer
from collections import OrderedDict


module_dir = os.path.dirname(os.path.abspath(__file__))

MPSYMBOL_TO_CN = AllCoordinationGeometries().get_symbol_cn_mapping()
ALLCG = AllCoordinationGeometries()


class StrategyOption(with_metaclass(abc.ABCMeta, MSONable)):

    allowed_values = None

    @abc.abstractmethod
    def as_dict(self):
        """
        A JSON serializable dict representation of this strategy option.
        """
        pass


class DistanceCutoffFloat(float, StrategyOption):

    allowed_values = 'Real number between 1.0 and +infinity'

    def __new__(cls, myfloat):
        flt = float.__new__(cls, myfloat)
        if flt < 1.0:
            raise ValueError("Distance cutoff should be between 1.0 and +infinity")
        return flt

    def as_dict(self):
        return {'@module': self.__class__.__module__,
                '@class': self.__class__.__name__,
                'value': self}

    @classmethod
    def from_dict(cls, d):
        return cls(d['value'])


class AngleCutoffFloat(float, StrategyOption):

    allowed_values = 'Real number between 0.0 and 1.0'

    def __new__(cls, myfloat):
        flt = float.__new__(cls, myfloat)
        if flt < 0.0 or flt > 1.0:
            raise ValueError("Angle cutoff should be between 0.0 and 1.0")
        return flt

    def as_dict(self):
        return {'@module': self.__class__.__module__,
                '@class': self.__class__.__name__,
                'value': self}

    @classmethod
    def from_dict(cls, d):
        return cls(d['value'])


class CSMFloat(float, StrategyOption):

    allowed_values = 'Real number between 0.0 and 100.0'

    def __new__(cls, myfloat):
        flt = float.__new__(cls, myfloat)
        if flt < 0.0 or flt > 100.0:
            raise ValueError("Continuous symmetry measure limits should be between 0.0 and 100.0")
        return flt

    def as_dict(self):
        return {'@module': self.__class__.__module__,
                '@class': self.__class__.__name__,
                'value': self}

    @classmethod
    def from_dict(cls, d):
        return cls(d['value'])


class AdditionalConditionInt(int, StrategyOption):

    allowed_values = 'Integer amongst :\n'
    for integer, description in AdditionalConditions.CONDITION_DESCRIPTION.items():
        allowed_values += ' - {:d} for "{}"\n'.format(integer, description)

    def __new__(cls, integer):
        intger = int.__new__(cls, integer)
        if intger not in AdditionalConditions.ALL:
            raise ValueError("Additional condition {:d} is not allowed".format(integer))
        return intger

    def as_dict(self):
        return {'@module': self.__class__.__module__,
                '@class': self.__class__.__name__,
                'value': self}

    @classmethod
    def from_dict(cls, d):
        return cls(d['value'])


class AbstractChemenvStrategy(with_metaclass(abc.ABCMeta, MSONable)):
    """
    Class used to define a Chemenv strategy for the neighbors and coordination environment to be applied to a
    StructureEnvironments object
    """
    DETAILED_VORONOI_CONTAINER = 'DetailedVoronoiContainer'
    ALLOWED_VORONOI_CONTAINERS = []
    AC = AdditionalConditions()
    STRATEGY_OPTIONS = OrderedDict()
    STRATEGY_DESCRIPTION = None

    def __init__(self, structure_environments=None):
        """
        Abstract constructor for the all chemenv strategies.
        :param structure_environments: StructureEnvironments object containing all the information on the
        coordination of the sites in a structure
        """
        self.structure_environments = None
        if structure_environments is not None:
            self.set_structure_environments(structure_environments)

    def set_structure_environments(self, structure_environments):
        self.structure_environments = structure_environments
        if not isinstance(self.structure_environments.voronoi, DetailedVoronoiContainer):
            raise ValueError('Voronoi Container not of type "DetailedVoronoiContainer"')
        self.prepare_symmetries()

    def prepare_symmetries(self):
        self.spg_analyzer = SpacegroupAnalyzer(self.structure_environments.structure)
        self.symops = self.spg_analyzer.get_symmetry_operations()

    def equivalent_site_index(self, psite):
        return self.equivalent_site_index_and_transform(psite)[0]
        #return self.structure_environments.structure.index(psite.to_unit_cell)

    def equivalent_site_index_and_transform(self, psite):
        # Get the index of the site in the unit cell of which the PeriodicSite psite is a replica.
        try:
            isite = self.structure_environments.structure.index(psite)
        except ValueError:
            try:
                uc_psite = psite.to_unit_cell
                isite = self.structure_environments.structure.index(uc_psite)
            except ValueError:
                for isite2, site2 in enumerate(self.structure_environments.structure):
                    if psite.is_periodic_image(site2):
                        isite = isite2
                        break
        # Get the translation between psite and its corresponding site in the unit cell (Translation I)
        thissite = self.structure_environments.structure[isite]
        dthissite = psite.frac_coords - thissite.frac_coords
        # Get the translation between the equivalent site for which the neighbors have been computed and the site in
        # the unit cell that corresponds to psite (Translation II)
        equivsite = self.structure_environments.structure[self.structure_environments.sites_map[isite]].to_unit_cell
        #equivsite = self.structure_environments.structure[self.structure_environments.sites_map[isite]]
        dequivsite = (self.structure_environments.structure[self.structure_environments.sites_map[isite]].frac_coords
                      - equivsite.frac_coords)
        found = False
        # Find the symmetry that applies the site in the unit cell to the equivalent site, as well as the translation
        # that gets back the site to the unit cell (Translation III)
        #TODO: check that these tolerances are needed, now that the structures are refined before analyzing environments
        tolerances = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
        for tolerance in tolerances:
            for symop in self.symops:
                newsite = PeriodicSite(equivsite._species, symop.operate(equivsite.frac_coords), equivsite._lattice)
                if newsite.is_periodic_image(thissite, tolerance=tolerance):
                    mysym = symop
                    dthissite2 = thissite.frac_coords - newsite.frac_coords
                    found = True
                    break
            if not found:
                symops = [SymmOp.from_rotation_and_translation()]
                for symop in symops:
                    newsite = PeriodicSite(equivsite._species, symop.operate(equivsite.frac_coords), equivsite._lattice)
                    #if newsite.is_periodic_image(thissite):
                    if newsite.is_periodic_image(thissite, tolerance=tolerance):
                        mysym = symop
                        dthissite2 = thissite.frac_coords - newsite.frac_coords
                        found = True
                        break
            if found:
                break
        if not found:
            raise EquivalentSiteSearchError(psite)
        return [self.structure_environments.sites_map[isite], dequivsite, dthissite + dthissite2, mysym]

    def __str__(self):
        return self.__class__.__name__

    def apply_strategy(self):
        """
        Applies the strategy to the structure_environments object in order to define the coordination environments of
        the sites in the structure_environments object
        :param structure_environments: StructureEnvironments object containing all the information needed to define the
        coordination environment of the sites in the structure
        :return: The coordination environment(s) of the sites
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_site_neighbors(self, site):
        """
        Applies the strategy to the structure_environments object in order to get the neighbors of a given site.
        :param site: Site for which the neighbors are looked for
        :param structure_environments: StructureEnvironments object containing all the information needed to get the
        neighbors of the site
        :return: The list of neighbors of the site. For complex strategies, where one allows multiple solutions, this
        can return a list of list of neighbors
        """
        return None

    @property
    def uniquely_determines_coordination_environments(self):
        """
        Returns True if the strategy leads to a unique coordination environment, False otherwise.
        :return: True if the strategy leads to a unique coordination environment, False otherwise.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_site_coordination_environment(self, site):
        """
        Applies the strategy to the structure_environments object in order to define the coordination environment of
        a given site.
        :param site: Site for which the coordination environment is looked for
        :return: The coordination environment of the site. For complex strategies, where one allows multiple
        solutions, this can return a list of coordination environments for the site
        """
        return None

    @abc.abstractmethod
    def get_site_coordination_environments(self, site):
        """
        Applies the strategy to the structure_environments object in order to define the coordination environment of
        a given site.
        :param site: Site for which the coordination environment is looked for
        :return: The coordination environment of the site. For complex strategies, where one allows multiple
        solutions, this can return a list of coordination environments for the site
        """
        return None

    def get_site_ce_fractions_and_neighbors(self, site, full_ce_info=False):
        """
        Applies the strategy to the structure_environments object in order to get coordination environments, their
        fraction, csm, geometry_info, and neighbors
        :param site: Site for which the above information is seeked
        :return: The list of neighbors of the site. For complex strategies, where one allows multiple solutions, this
        can return a list of list of neighbors
        """
        [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        if self.uniquely_determines_coordination_environments:
            geom_and_map = self.get_site_coordination_environment(site=site, isite=isite, dequivsite=dequivsite,
                                                                  dthissite=dthissite, mysym=mysym, return_map=True)
            if geom_and_map is None:
                return None
            ce_and_neighbors = {'ce': [], 'neighbors': {}}
            mingeom = geom_and_map[0]
            cn_map = geom_and_map[1]
            if cn_map[0] > 12 or cn_map[0] == 0:
                return None
            tuple_cn_map = tuple(cn_map)

            if tuple_cn_map not in ce_and_neighbors['neighbors']:
                ce_and_neighbors['neighbors'][tuple_cn_map] = (self.structure_environments.voronoi.
                                                               unique_coordinated_neighbors(isite=isite,
                                                                                            cn_map=cn_map))[0]
            if mingeom is not None:
                geom_dict = {'mp_symbol': mingeom[0], 'fraction': 1.0,
                             'cn_map': tuple_cn_map, 'csm': mingeom[1]['symmetry_measure']}
                if full_ce_info:
                    geom_dict['coordination_geometry_info'] = mingeom[1]
            else:
                geom_dict = {'mp_symbol': None, 'fraction': None,
                             'cn_map': tuple_cn_map, 'csm': None}
            ce_and_neighbors['ce'].append(geom_dict)
        else:
            geoms_and_maps_list = self.get_site_coordination_environments_fractions(site=site, isite=isite,
                                                                                    dequivsite=dequivsite,
                                                                                    dthissite=dthissite, mysym=mysym,
                                                                                    return_maps=True)
            if geoms_and_maps_list is None:
                return None
            ce_and_neighbors = {'ce': [], 'neighbors': {}}
            for ce_symbol, ce_dict, ce_fraction, cn_map in geoms_and_maps_list:
                tuple_cn_map = tuple(cn_map)
                if tuple_cn_map not in ce_and_neighbors['neighbors']:
                    ce_and_neighbors['neighbors'][tuple_cn_map] = (self.structure_environments.
                                                                   unique_coordinated_neighbors(isite=isite,
                                                                                                cn_map=cn_map))[0]
                geom_dict = {'mp_symbol': ce_symbol, 'fraction': ce_fraction,
                             'cn_map': tuple_cn_map, 'csm': ce_dict['symmetry_measure']}
                if full_ce_info:
                    geom_dict['coordination_geometry_info'] = ce_dict
                ce_and_neighbors['ce'].append(geom_dict)
        return ce_and_neighbors

    def structure_has_environment(self, mp_symbol, unequivocal=True):
        """
        Checks whether the structure contains the environment symbolized by mp_symbol. For strategies allowing
        mixed environments, the unequivocal argument specifies whether the check should be on the most probable
        environment or on all of the possible environments.
        :param mp_symbol:
        :param unequivocal:
        :return:
        """
        symmetrized_structure = self.spg_analyzer.get_symmetrized_structure()
        for sites_group in symmetrized_structure.equivalent_sites:
            site = sites_group[0]
            if self.uniquely_determines_coordination_environments or unequivocal:
                ce, ce_dict = self.get_site_coordination_environment(site)
                if ce == mp_symbol:
                    return True
            else:
                allce = self.get_site_coordination_environments(site)
                if mp_symbol in [ce for ce, ce_dict in allce]:
                    return True
        return False

    def set_option(self, option_name, option_value):
        self.__setattr__(option_name, option_value)

    def setup_options(self, all_options_dict):
        for option_name, option_value in all_options_dict.items():
            self.set_option(self.STRATEGY_OPTIONS[option_name]['internal'], option_value)

    @abc.abstractmethod
    def __eq__(self, other):
        """
        Equality method that should be implemented for any strategy
        :param other: strategy to be compared with the current one
        :return:
        """
        return

    def __str__(self):
        out = '  Chemenv Strategy "{}"\n'.format(self.__class__.__name__)
        out += '  {}\n\n'.format('='*(19+len(self.__class__.__name__)))
        out += '  Description :\n  {}\n'.format('-'*13)
        out += self.STRATEGY_DESCRIPTION
        out += '\n\n'
        out += '  Options :\n  {}\n'.format('-'*9)
        for option_name, option_dict in self.STRATEGY_OPTIONS.items():
            out += '   - {} : {}\n'.format(option_name, str(getattr(self, option_name)))
        return out

    def as_dict(self):
        """
        Bson-serializable dict representation of the SimplestChemenvStrategy object.
        :return: Bson-serializable dict representation of the SimplestChemenvStrategy object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the SimpleAbundanceChemenvStrategy object from a dict representation of the
        SimpleAbundanceChemenvStrategy object created using the as_dict method.
        :param d: dict representation of the SimpleAbundanceChemenvStrategy object
        :return: StructureEnvironments object
        """
        return cls()


class SimplestChemenvStrategy(AbstractChemenvStrategy):
    """
    Simplest ChemenvStrategy using fixed angle and distance parameters for the definition of neighbors in the
    Voronoi approach. The coordination environment is then given as the one with the lowest continuous symmetry measure
    """

    # Default values for the distance and angle cutoffs
    DEFAULT_DISTANCE_CUTOFF = 1.4
    DEFAULT_ANGLE_CUTOFF = 0.3
    DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF = 10.0
    DEFAULT_ADDITIONAL_CONDITION = AbstractChemenvStrategy.AC.ONLY_ACB
    ALLOWED_VORONOI_CONTAINERS = [AbstractChemenvStrategy.DETAILED_VORONOI_CONTAINER]
    STRATEGY_OPTIONS = OrderedDict({'distance_cutoff': {'type': DistanceCutoffFloat, 'internal': '_distance_cutoff',
                                                        'default': DEFAULT_DISTANCE_CUTOFF},
                                    'angle_cutoff': {'type': AngleCutoffFloat, 'internal': '_angle_cutoff',
                                                     'default': DEFAULT_ANGLE_CUTOFF},
                                    'additional_condition': {'type': AdditionalConditionInt,
                                                             'internal': '_additional_condition',
                                                             'default': DEFAULT_ADDITIONAL_CONDITION},
                                    'continuous_symmetry_measure_'
                                    'cutoff': {'type': CSMFloat,
                                               'internal': '_continuous_symmetry_measure_cutoff',
                                               'default': DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF}})
    STRATEGY_DESCRIPTION = '    Simplest ChemenvStrategy using fixed angle and distance parameters \n' \
                           '    for the definition of neighbors in the Voronoi approach. \n' \
                           '    The coordination environment is then given as the one with the \n' \
                           '    lowest continuous symmetry measure.'

    def __init__(self, structure_environments=None, distance_cutoff=DEFAULT_DISTANCE_CUTOFF,
                 angle_cutoff=DEFAULT_ANGLE_CUTOFF, additional_condition=DEFAULT_ADDITIONAL_CONDITION,
                 continuous_symmetry_measure_cutoff=DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF):
        """
        Constructor for this SimplestChemenvStrategy.
        :param distance_cutoff: Distance cutoff used
        :param angle_cutoff: Angle cutoff used
        """
        AbstractChemenvStrategy.__init__(self, structure_environments)
        self._distance_cutoff = DistanceCutoffFloat(distance_cutoff)
        self._angle_cutoff = AngleCutoffFloat(angle_cutoff)
        self._additional_condition = AdditionalConditionInt(additional_condition)
        self._continuous_symmetry_measure_cutoff = CSMFloat(continuous_symmetry_measure_cutoff)

    @property
    def uniquely_determines_coordination_environments(self):
        return True

    @property
    def distance_cutoff(self):
        return self._distance_cutoff

    @property
    def angle_cutoff(self):
        return self._angle_cutoff

    @property
    def additional_condition(self):
        return self._additional_condition

    @property
    def continuous_symmetry_measure_cutoff(self):
        return self._continuous_symmetry_measure_cutoff

    def apply_strategy(self):
        return None

    def get_site_neighbors(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None):#, neighbors_map=None):
        #if neighbors_map is not None:
        #    return self.structure_environments.voronoi.get_neighbors(isite=isite, neighbors_map=neighbors_map)
        if isite is None:
            [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        eqsite_ps = self.structure_environments.voronoi.neighbors(isite=isite,
                                                                  distfactor=self.distance_cutoff,
                                                                  angfactor=self.angle_cutoff,
                                                                  additional_condition=
                                                                  self._additional_condition)
        ce = self.get_site_coordination_environment(site=site, isite=isite, dequivsite=dequivsite, dthissite=dthissite, mysym=mysym)
        uniquenbs = self.structure_environments.voronoi.unique_coordinated_neighbors(isite=isite)
        detailed_voronoi_index = ce[1]['detailed_voronoi_index']
        eqsite_ps = uniquenbs[detailed_voronoi_index['cn']][detailed_voronoi_index['index']][0]

        coordinated_neighbors = []
        for ips, ps in enumerate(eqsite_ps):
            coords = mysym.operate(ps.frac_coords + dequivsite) + dthissite
            ps_site = PeriodicSite(ps._species, coords, ps._lattice)
            coordinated_neighbors.append(ps_site)
        return coordinated_neighbors

    def get_site_coordination_environment(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None,
                                          return_map=False):
        if isite is None:
            [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        neighbors_map = self.structure_environments.voronoi.neighbors_map(isite=isite,
                                                                          distfactor=self.distance_cutoff,
                                                                          angfactor=self.angle_cutoff,
                                                                          additional_condition=self.AC.ONLY_ACB)
        if neighbors_map is None:
            return None
        cn_map = (self.structure_environments.voronoi.parameters_to_unique_coordinated_neighbors_map
                  [isite]
                  [neighbors_map['i_distfactor']][neighbors_map['i_angfactor']]
                  [neighbors_map['i_additional_condition']])
        coord_geoms = (self.structure_environments.ce_list[self.structure_environments.sites_map[isite]]
                       [cn_map[0]][cn_map[1]].coord_geoms)
        if return_map:
            if coord_geoms is None:
                return cn_map[0], cn_map
            return (self.structure_environments.ce_list[self.structure_environments.sites_map[isite]][cn_map[0]][
                cn_map[1]].minimum_geometry(), cn_map)
        else:
            if coord_geoms is None:
                return cn_map[0]
            return self.structure_environments.ce_list[self.structure_environments.sites_map[isite]][cn_map[0]][
                cn_map[1]].minimum_geometry()

    def get_site_coordination_environments(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None,
                                           return_maps=False):
        return [self.get_site_coordination_environment(site=site, isite=isite, dequivsite=dequivsite,
                                                       dthissite=dthissite, mysym=mysym, return_map=return_maps)]

    def __eq__(self, other):
        return (self.__class__.__name__ == other.__class__.__name__ and
                self._distance_cutoff == other._distance_cutoff and self._angle_cutoff == other._angle_cutoff and
                self._additional_condition == other._additional_condition and
                self._continuous_symmetry_measure_cutoff == other._continuous_symmetry_measure_cutoff)

    def as_dict(self):
        """
        Bson-serializable dict representation of the SimplestChemenvStrategy object.
        :return: Bson-serializable dict representation of the SimplestChemenvStrategy object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "distance_cutoff": float(self._distance_cutoff),
                "angle_cutoff": float(self._angle_cutoff),
                "additional_condition": int(self._additional_condition),
                "continuous_symmetry_measure_cutoff": float(self._continuous_symmetry_measure_cutoff)}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the SimplestChemenvStrategy object from a dict representation of the SimplestChemenvStrategy object
        created using the as_dict method.
        :param d: dict representation of the SimplestChemenvStrategy object
        :return: StructureEnvironments object
        """
        return cls(distance_cutoff=d["distance_cutoff"], angle_cutoff=d["angle_cutoff"],
                   additional_condition=d["additional_condition"],
                   continuous_symmetry_measure_cutoff=d["continuous_symmetry_measure_cutoff"])


class SimpleAbundanceChemenvStrategy(AbstractChemenvStrategy):
    """
    Simple ChemenvStrategy using the neighbors that are the most "abundant" in the grid of angle and distance
    parameters for the definition of neighbors in the Voronoi approach.
    The coordination environment is then given as the one with the lowest continuous symmetry measure
    """

    ALLOWED_VORONOI_CONTAINERS = [AbstractChemenvStrategy.DETAILED_VORONOI_CONTAINER]
    DEFAULT_MAX_DIST = 2.0
    DEFAULT_ADDITIONAL_CONDITION = AbstractChemenvStrategy.AC.ONLY_ACB
    STRATEGY_OPTIONS = OrderedDict({'additional_condition': {'type': AdditionalConditionInt,
                                                             'internal': '_additional_condition',
                                                             'default': DEFAULT_ADDITIONAL_CONDITION},
                                    'surface_calculation_type': {}})
    STRATEGY_DESCRIPTION = '    Simple Abundance ChemenvStrategy using the most "abundant" neighbors map \n' \
                           '    for the definition of neighbors in the Voronoi approach. \n' \
                           '    The coordination environment is then given as the one with the \n' \
                           '    lowest continuous symmetry measure.'

    def __init__(self, structure_environments=None,
                 additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB):
        """
        Constructor for the SimpleAbundanceChemenvStrategy.
        :param structure_environments: StructureEnvironments object containing all the information on the
        coordination of the sites in a structure
        """
        AbstractChemenvStrategy.__init__(self, structure_environments)
        self._additional_condition = additional_condition

    @property
    def uniquely_determines_coordination_environments(self):
        return True

    def get_site_neighbors(self, site):
        [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        cn_map = self._get_map(isite)
        eqsite_ps = (self.structure_environments.unique_coordinated_neighbors(isite, cn_map=cn_map))
        coordinated_neighbors = []
        for ips, ps in enumerate(eqsite_ps):
            coords = mysym.operate(ps.frac_coords + dequivsite) + dthissite
            ps_site = PeriodicSite(ps._species, coords, ps._lattice)
            coordinated_neighbors.append(ps_site)
        return coordinated_neighbors

    def get_site_coordination_environment(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None,
                                          return_map=False):
        if isite is None:
            [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        cn_map = self._get_map(isite)
        if cn_map is None:
            return None
        coord_geoms = (self.structure_environments.
                       ce_list[self.structure_environments.sites_map[isite]][cn_map[0]][cn_map[1]])
        if return_map:
            if coord_geoms is None:
                return cn_map[0], cn_map
            return coord_geoms.minimum_geometry(), cn_map
        else:
            if coord_geoms is None:
                return cn_map[0]
            return coord_geoms.minimum_geometry()

    def get_site_coordination_environments(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None,
                                           return_maps=False):
        return [self.get_site_coordination_environment(site=site, isite=isite, dequivsite=dequivsite,
                                                       dthissite=dthissite, mysym=mysym, return_map=return_maps)]

    def _get_map(self, isite):
        maps_and_surfaces = self._get_maps_surfaces(isite)
        if maps_and_surfaces is None:
            return None
        surface_max = 0.0
        imax = -1
        for ii, map_and_surface in enumerate(maps_and_surfaces):
            all_additional_conditions = [ac[2] for ac in map_and_surface['parameters_indices']]
            if self._additional_condition in all_additional_conditions and map_and_surface['surface'] > surface_max:
                surface_max = map_and_surface['surface']
                imax = ii
        return maps_and_surfaces[imax]['map']

    def _get_maps_surfaces(self, isite, surface_calculation_type=None):
        if surface_calculation_type is None:
            surface_calculation_type = {'distance_parameter': ('initial_normalized', None),
                                        'angle_parameter': ('initial_normalized', None)}
        return self.structure_environments.voronoi.maps_and_surfaces(isite=isite,
                                                                     surface_calculation_type=surface_calculation_type,
                                                                     max_dist=self.DEFAULT_MAX_DIST)

    def __eq__(self, other):
        return (self.__class__.__name__ == other.__class__.__name__ and
                self._additional_condition == other.additional_condition)

    def as_dict(self):
        """
        Bson-serializable dict representation of the SimpleAbundanceChemenvStrategy object.
        :return: Bson-serializable dict representation of the SimpleAbundanceChemenvStrategy object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "additional_condition": self._additional_condition}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the SimpleAbundanceChemenvStrategy object from a dict representation of the
        SimpleAbundanceChemenvStrategy object created using the as_dict method.
        :param d: dict representation of the SimpleAbundanceChemenvStrategy object
        :return: StructureEnvironments object
        """
        return cls(additional_condition=d["additional_condition"])


class TargettedPenaltiedAbundanceChemenvStrategy(SimpleAbundanceChemenvStrategy):
    """
    Simple ChemenvStrategy using the neighbors that are the most "abundant" in the grid of angle and distance
    parameters for the definition of neighbors in the Voronoi approach, with a bias for a given list of target
    environments. This can be useful in the case of, e.g. connectivity search of some given environment.
    The coordination environment is then given as the one with the lowest continuous symmetry measure
    """
    ALLOWED_VORONOI_CONTAINERS = [AbstractChemenvStrategy.DETAILED_VORONOI_CONTAINER]
    DEFAULT_TARGET_ENVIRONMENTS = ['O:6']

    def __init__(self, structure_environments=None, truncate_dist_ang=True,
                 additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB,
                 max_nabundant=5, target_environments=DEFAULT_TARGET_ENVIRONMENTS, target_penalty_type='max_csm',
                 max_csm=5.0):
        SimpleAbundanceChemenvStrategy.__init__(self, structure_environments,
                                                additional_condition=additional_condition)
        self.max_nabundant = max_nabundant
        self.target_environments = target_environments
        self.target_penalty_type = target_penalty_type
        self.max_csm = max_csm

    def get_site_coordination_environment(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None,
                                          return_map=False):
        if isite is None:
            [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        cn_map = self._get_map(isite)
        if cn_map is None:
            return None
        chemical_environments = (self.structure_environments.ce_list
                                 [self.structure_environments.sites_map[isite]][cn_map[0]][cn_map[1]])

        if return_map:
            if chemical_environments.coord_geoms is None or len(chemical_environments) == 0:
                return cn_map[0], cn_map
            return chemical_environments.minimum_geometry(), cn_map
        else:
            if chemical_environments.coord_geoms is None:
                return cn_map[0]
            return chemical_environments.minimum_geometry()

    def _get_map(self, isite):
        maps_and_surfaces = SimpleAbundanceChemenvStrategy._get_maps_surfaces(self, isite)
        if maps_and_surfaces is None:
            return SimpleAbundanceChemenvStrategy._get_map(self, isite)
        current_map = None
        current_target_env_csm = 100.0
        surfaces = [map_and_surface['surface'] for map_and_surface in maps_and_surfaces]
        order = np.argsort(surfaces)[::-1]
        target_cgs = [AllCoordinationGeometries().get_geometry_from_mp_symbol(mp_symbol)
                      for mp_symbol in self.target_environments]
        target_cns = [cg.coordination_number for cg in target_cgs]
        for ii in range(min([len(maps_and_surfaces), self.max_nabundant])):
            my_map_and_surface = maps_and_surfaces[order[ii]]
            mymap = my_map_and_surface['map']
            cn = mymap[0]
            if cn not in target_cns or cn > 12 or cn == 0:
                continue
            all_conditions = [params[2] for params in my_map_and_surface['parameters_indices']]
            if self._additional_condition not in all_conditions:
                continue
            cg, cgdict = (self.structure_environments.ce_list
                          [self.structure_environments.sites_map[isite]][mymap[0]][mymap[1]].minimum_geometry())
            if (cg in self.target_environments and cgdict['symmetry_measure'] <= self.max_csm and
                        cgdict['symmetry_measure'] < current_target_env_csm):
                current_map = mymap
                current_target_env_csm = cgdict['symmetry_measure']
        if current_map is not None:
            return current_map
        else:
            return SimpleAbundanceChemenvStrategy._get_map(self, isite)

    @property
    def uniquely_determines_coordination_environments(self):
        return True

    def as_dict(self):
        """
        Bson-serializable dict representation of the TargettedPenaltiedAbundanceChemenvStrategy object.
        :return: Bson-serializable dict representation of the TargettedPenaltiedAbundanceChemenvStrategy object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "additional_condition": self._additional_condition,
                "max_nabundant": self.max_nabundant,
                "target_environments": self.target_environments,
                "target_penalty_type": self.target_penalty_type,
                "max_csm": self.max_csm}

    def __eq__(self, other):
        return (self.__class__.__name__ == other.__class__.__name__ and
                self._additional_condition == other.additional_condition and
                self.max_nabundant == other.max_nabundant and
                self.target_environments == other.target_environments and
                self.target_penalty_type == other.target_penalty_type and
                self.max_csm == other.max_csm)

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the TargettedPenaltiedAbundanceChemenvStrategy object from a dict representation of the
        TargettedPenaltiedAbundanceChemenvStrategy object created using the as_dict method.
        :param d: dict representation of the TargettedPenaltiedAbundanceChemenvStrategy object
        :return: TargettedPenaltiedAbundanceChemenvStrategy object
        """
        return cls(additional_condition=d["additional_condition"],
                   max_nabundant=d["max_nabundant"],
                   target_environments=d["target_environments"],
                   target_penalty_type=d["target_penalty_type"],
                   max_csm=d["max_csm"])


#Read in BV parameters.
DEFAULT_CSM_CUTOFFS = {}
with open(os.path.join(module_dir, "strategy_files", "ImprovedConfidenceCutoffDefaultParameters.json"), "r") as f:
    dd = json.load(f)
    for mp_symbol, cutoff in dd['csm_cutoffs'].items():
        DEFAULT_CSM_CUTOFFS[mp_symbol] = cutoff
    DEFAULT_VORONOI_PARAMETERS_FRACTIONS = dd['voronoi_parameters_fractions']


class ImprovedConfidenceCutoffChemenvStrategy(AbstractChemenvStrategy):
    """
    ChemenvStrategy that only returns an environment when it is clear enough within a given "confidence cutoff".
    This confidence cutoff is given as a cutoff in the continuous symmetry measure, as well as something about
    the Voronoi ...
    """
    DEFAULT_ADDITIONAL_CONDITION = AbstractChemenvStrategy.AC.ONLY_ACB
    ALLOWED_VORONOI_CONTAINERS = []

    def __init__(self, structure_environments=None, additional_condition=DEFAULT_ADDITIONAL_CONDITION,
                 csm_cutoffs=None,
                 voronoi_parameters_fractions=None):
        """
        Constructor for the SimpleAbundanceChemenvStrategy.
        :param structure_environments: StructureEnvironments object containing all the information on the
        coordination of the sites in a structure
        """
        AbstractChemenvStrategy.__init__(self, structure_environments)
        self._additional_condition = additional_condition
        if csm_cutoffs is None:
            self.csm_cutoffs = DEFAULT_CSM_CUTOFFS
        else:
            self.csm_cutoffs = csm_cutoffs
        if voronoi_parameters_fractions is None:
            self.voronoi_parameters_fractions = DEFAULT_VORONOI_PARAMETERS_FRACTIONS
        else:
            self.voronoi_parameters_fractions = voronoi_parameters_fractions

    def get_site_coordination_environment(self, site):
        [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        maps, counts = self._get_maps_counts(isite)
        maps_minimum_geometries = {}
        for imap, some_map in enumerate(maps):
            if some_map is None:
                continue
            if (self.structure_environments.ce_list[
                    self.structure_environments.sites_map[isite]][some_map[0]][some_map[1]].
                        coord_geoms is None):
                continue
            maps_minimum_geometries[some_map] = (
                self.structure_environments.ce_list[self.structure_environments.sites_map[isite]][
                    some_map[0]][some_map[1]].minimum_geometry(), counts(imap))
        return UNCLEAR_ENVIRONMENT_SYMBOL, None

    def __eq__(self, other):
        return (self.__class__.__name__ == other.__class__.__name__ and
                self._additional_condition == other._additional_condition and
                self.csm_cutoffs == other.csm_cutoffs and
                self.voronoi_parameters_fractions == other.voronoi_parameters_cutoffs)

    def _get_maps(self, isite):
        return None

    def as_dict(self):
        """
        Bson-serializable dict representation of the ImprovedConfidenceCutoffChemenvStrategy object.
        :return: Bson-serializable dict representation of the ImprovedConfidenceCutoffChemenvStrategy object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "only_anion_cation_bonds": self.only_anion_cation_bonds,
                "csm_cutoff": self.csm_cutoffs,
                "voronoi_parameters_fraction": self.voronoi_parameters_fractions}

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the ImprovedConfidenceCutoffChemenvStrategy object from a dict representation of the
        ImprovedConfidenceCutoffChemenvStrategy object created using the as_dict method.
        :param d: dict representation of the ImprovedConfidenceCutoffChemenvStrategy object
        :return: StructureEnvironments object
        """
        return cls(only_anion_cation_bonds=d["only_anion_cation_bonds"],
                   csm_cutoffs=d["csm_cutoffs"],
                   voronoi_parameters_fractions=d["voronoi_parameters_fractions"])


class ComplexCSMBasedChemenvStrategy(AbstractChemenvStrategy):
    """
    ChemenvStrategy giving a percentage for each environment.
    #TODO: document how this is performed exactly
    """
    ALLOWED_VORONOI_CONTAINERS = [AbstractChemenvStrategy.DETAILED_VORONOI_CONTAINER]
    ALLOWED_CN_DELTA_MEAN_CSM_ESTIMATOR_CONCATENATORS = {'product': np.product,
                                                         'minimum': np.min,
                                                         'min': np.min}
    DEFAULT_ADDITIONAL_CONDITION = AbstractChemenvStrategy.AC.ONLY_ACB
    DEFAULT_MAX_CSM = 8.0
    DEFAULT_MEAN_CSM_ESTIMATOR = ('power2_inverse_decreasing', {'max_csm': DEFAULT_MAX_CSM})
    DEFAULT_CN_SELF_MEAN_CSM_ESTIMATOR = ('power2_decreasing_exp', {'max_csm': DEFAULT_MAX_CSM,
                                                                    'alpha': 1.0})
    DEFAULT_CN_DELTA_MEAN_CSM_ESTIMATOR = ('smootherstep', {'delta_csm_min': 1.0,
                                                            'delta_csm_max': 4.0})
    DEFAULT_CN_SELF_MEAN_CSM_ESTIMATOR_CN_SPECIFICS = {}
    DEFAULT_CN_SELF_MEAN_CSM_ESTIMATOR_CE_SPECIFICS = {}
    DEFAULT_CN_DELTA_MEAN_CSM_ESTIMATOR_CN_SPECIFICS = {}
    DEFAULT_CN_DELTA_MEAN_CSM_ESTIMATOR_CE_SPECIFICS = {}
    DEFAULT_CE_ESTIMATOR = ('power2_inverse_power2_decreasing', {'max_csm': DEFAULT_MAX_CSM})
    DEFAULT_CN_DELTA_MEAN_CSM_ESTIMATOR_CONCATENATOR = 'product'

    def __init__(self, structure_environments=None, additional_condition=DEFAULT_ADDITIONAL_CONDITION,
                 mean_csm_estimator=DEFAULT_MEAN_CSM_ESTIMATOR,
                 cn_self_mean_csm_estimator=DEFAULT_CN_SELF_MEAN_CSM_ESTIMATOR,
                 cn_self_mean_csm_estimator_cn_specifics=DEFAULT_CN_SELF_MEAN_CSM_ESTIMATOR_CN_SPECIFICS,
                 cn_self_mean_csm_estimator_ce_specifics=DEFAULT_CN_SELF_MEAN_CSM_ESTIMATOR_CE_SPECIFICS,
                 cn_delta_mean_csm_estimator=DEFAULT_CN_DELTA_MEAN_CSM_ESTIMATOR,
                 cn_delta_mean_csm_estimator_cn_specifics=DEFAULT_CN_DELTA_MEAN_CSM_ESTIMATOR_CN_SPECIFICS,
                 cn_delta_mean_csm_estimator_ce_specifics=DEFAULT_CN_DELTA_MEAN_CSM_ESTIMATOR_CE_SPECIFICS,
                 cn_delta_mean_csm_estimator_concatenator=DEFAULT_CN_DELTA_MEAN_CSM_ESTIMATOR_CONCATENATOR,
                 ce_estimator=DEFAULT_CE_ESTIMATOR):
        AbstractChemenvStrategy.__init__(self, structure_environments)
        self._additional_condition = additional_condition
        #Definition of the estimator/function used to estimate the "mean" csm used for the calculation of
        # the fractions of each cn_map
        self._mean_csm_estimator = mean_csm_estimator
        self._mean_csm_estimator_ratio_function = CSMInfiniteRatioFunction(mean_csm_estimator[0],
                                                                           options_dict=mean_csm_estimator[1])
        self._mean_csm = self._mean_csm_estimator_ratio_function.mean_estimator
        #Definition of the estimator/function used to compute the "self" contribution to the fraction of the cn_map
        self._cn_self_mean_csm_estimator = cn_self_mean_csm_estimator
        self._cn_self_mean_csm_estimator_cn_specifics = cn_self_mean_csm_estimator_cn_specifics
        self._cn_self_mean_csm_estimator_ce_specifics = cn_self_mean_csm_estimator_ce_specifics
        self._cn_self_mean_csm_estimator_ratio_function = \
            CSMFiniteRatioFunction(cn_self_mean_csm_estimator[0],
                                   options_dict=cn_self_mean_csm_estimator[1])
        self._cn_self_mean_csm_evaluate = self._cn_self_mean_csm_estimator_ratio_function.evaluate
        #Definition of the estimator/function used to compute the "delta" contribution to the fraction of the cn_map
        self._cn_delta_mean_csm_estimator = cn_delta_mean_csm_estimator
        self._cn_delta_mean_csm_estimator_ratio_function = \
            DeltaCSMRatioFunction(cn_delta_mean_csm_estimator[0],
                                  options_dict=cn_delta_mean_csm_estimator[1])
        self._cn_delta_mean_csm_evaluate = self._cn_delta_mean_csm_estimator_ratio_function.evaluate
        #Definition of the CN specific estimator/function used to compute the "delta" contribution to
        #the fraction of the cn_map
        self._cn_delta_mean_csm_estimator_cn_specifics = cn_delta_mean_csm_estimator_cn_specifics
        if len(cn_delta_mean_csm_estimator_cn_specifics) == 0:
            self._cn_delta_mean_csm_is_cn_specific = False
        else:
            self._cn_delta_mean_csm_is_cn_specific = True
        self._cn_delta_mean_csm_estimator_cn_specifics_ratio_functions = {}
        self._cn_delta_mean_csm_estimator_cn_specifics_evaluate = {}
        self._maximum_continuous_symmetry_measure = self._cn_self_mean_csm_estimator[1]['max_csm']
        for cn1, cn2_functions in self._cn_delta_mean_csm_estimator_cn_specifics.items():
            if not 'other' in cn2_functions.keys():
                raise ValueError('Should have a "other" key in cn_delta_mean_csm_estimator_cn_specifics')
            self._cn_delta_mean_csm_estimator_cn_specifics_ratio_functions[cn1] = {}
            self._cn_delta_mean_csm_estimator_cn_specifics_evaluate[cn1] = {}
            for cn2, function in cn2_functions.items():
                self._cn_delta_mean_csm_estimator_cn_specifics_ratio_functions[cn1][cn2] = (
                    DeltaCSMRatioFunction(function[0], options_dict=function[1]))
                self._cn_delta_mean_csm_estimator_cn_specifics_evaluate[cn1][cn2] = \
                    self._cn_delta_mean_csm_estimator_cn_specifics_ratio_functions[cn1][cn2].evaluate
        #Definition of the CE specific estimator/function used to compute the "delta" contribution to
        #the fraction of the cn_map
        if len(cn_delta_mean_csm_estimator_ce_specifics) > 0:
            raise NotImplementedError('CE specifics estimators are not yet implemented')
        self._cn_delta_mean_csm_estimator_ce_specifics = cn_delta_mean_csm_estimator_ce_specifics
        self._cn_delta_mean_csm_estimator_concatenator = cn_delta_mean_csm_estimator_concatenator
        self._cn_delta_mean_csm_estimator_concatenator_function = (self.
                                                                   ALLOWED_CN_DELTA_MEAN_CSM_ESTIMATOR_CONCATENATORS
                                                                   [self._cn_delta_mean_csm_estimator_concatenator])
        #Definition of the estimator/function used to compute the fractions of each coordination environment within
        # one cn_map
        self._ce_estimator = ce_estimator
        self._ce_estimator_ratio_function = CSMInfiniteRatioFunction(ce_estimator[0], options_dict=ce_estimator[1])
        self._ce_estimator_fractions = self._ce_estimator_ratio_function.fractions

    @property
    def additional_condition(self):
        return self._additional_condition

    @property
    def uniquely_determines_coordination_environments(self):
        return False

    def get_site_coordination_environments_fractions(self, site, isite=None, dequivsite=None, dthissite=None,
                                                     mysym=None, ordered=True, min_fraction=0.0, return_maps=True):
        if isite is None or dequivsite is None or dthissite is None or mysym is None:
            [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        #self._setup_mean_csms(isite)
        #self._setup_mean_csms_deltas()
        all_f_cn_maps = self.get_f_cn_maps(isite)
        f_cn_maps = all_f_cn_maps['total']
        f_cn_maps_total = np.sum(list(f_cn_maps.values()))
        cn_maps_fractions = {cn_map: f_cn_map / f_cn_maps_total for cn_map, f_cn_map in f_cn_maps.items()}
        ce_symbols = []
        ce_dicts = []
        ce_fractions = []
        ce_maps = []
        for cn_map, cn_map_fraction in cn_maps_fractions.items():
            if cn_map_fraction > 0.0:
                mingeoms = self._mingeoms_isite[cn_map]
                csms = [ce_dict['symmetry_measure'] for ce_symbol, ce_dict in mingeoms]
                fractions = self._ce_estimator_fractions(csms)
                for ifraction, fraction in enumerate(fractions):
                    if fraction > 0.0:
                        ce_symbols.append(mingeoms[ifraction][0])
                        ce_dicts.append(mingeoms[ifraction][1])
                        ce_fractions.append(cn_map_fraction*fraction)
                        ce_maps.append(cn_map)
        if ordered:
            indices = np.argsort(ce_fractions)[::-1]
        else:
            indices = list(range(len(ce_fractions)))
        if return_maps:
            return [(ce_symbols[ii], ce_dicts[ii], ce_fractions[ii], ce_maps[ii])
                    for ii in indices if ce_fractions[ii] > min_fraction]
        else:
            return [(ce_symbols[ii], ce_dicts[ii], ce_fractions[ii])
                    for ii in indices if ce_fractions[ii] > min_fraction]

    def _get_evaluate_f_cn_maps_deltas(self, cn_map, mean_csm_deltas):
        if self._cn_delta_mean_csm_is_cn_specific:
            f_cn_maps_deltas = []
            cn1 = cn_map[0]
            for cn_map2, delta_mean_csm in mean_csm_deltas.items():
                cn2 = cn_map2[0]
                if cn2 > cn1:
                    if cn1 in self._cn_delta_mean_csm_estimator_cn_specifics_evaluate:
                        if cn2 in self._cn_delta_mean_csm_estimator_cn_specifics_evaluate[cn1]:
                            fdelta = self._cn_delta_mean_csm_estimator_cn_specifics_evaluate[cn1][cn2](delta_mean_csm)
                        else:
                            fdelta = (self._cn_delta_mean_csm_estimator_cn_specifics_evaluate
                                      [cn1]['other'](delta_mean_csm))
                    else:
                        fdelta = self._cn_delta_mean_csm_evaluate(delta_mean_csm)
                    f_cn_maps_deltas.append(fdelta)
                else:
                    f_cn_maps_deltas.append(1.0)
        else:
            f_cn_maps_deltas = [self._cn_delta_mean_csm_evaluate(delta_mean_csm)
                                for cn_map2, delta_mean_csm in mean_csm_deltas.items()]
        return f_cn_maps_deltas

    def _get_f_cn_maps_self_delta(self):
        f_cn_maps_self = {}
        f_cn_maps_delta = {}
        for cn_map, mean_csm in self._mean_csms_isite.items():
            f_cn_maps_self[cn_map] = self._cn_self_mean_csm_evaluate(mean_csm)
            mean_csm_deltas = self._mean_csms_deltas_isite[cn_map]
            f_cn_maps_deltas = self._get_evaluate_f_cn_maps_deltas(cn_map, mean_csm_deltas)
            if len(f_cn_maps_deltas) == 0:
                f_cn_maps_delta[cn_map] = 1.0
            else:
                f_cn_maps_delta[cn_map] = self._cn_delta_mean_csm_estimator_concatenator_function(f_cn_maps_deltas)
        return f_cn_maps_self, f_cn_maps_delta

    def get_f_cn_maps(self, isite):
        self._setup_mean_csms(isite)
        self._setup_mean_csms_deltas()
        f_cn_maps_self, f_cn_maps_delta = self._get_f_cn_maps_self_delta()
        f_cn_maps = {}
        for cn_map, mean_csm in self._mean_csms_isite.items():
            f_cn_maps[cn_map] = f_cn_maps_self[cn_map] * f_cn_maps_delta[cn_map]
        return {'self': f_cn_maps_self, 'delta': f_cn_maps_delta, 'total': f_cn_maps}

    def _setup_mean_csms(self, isite):
        self._mean_csms_isite = {}
        self._mingeoms_isite = {}
        for cn, cn_coordnbs_list in self.structure_environments.unique_coordinated_neighbors(isite).items():
            if cn > 12 or cn == 0:
                continue
            for i_coordnbs in range(len(cn_coordnbs_list)):
                if not (self.structure_environments.voronoi.satisfy_condition(isite, cn, i_coordnbs,
                                                                              self._additional_condition)):
                    continue
                mingeoms = self.structure_environments.ce_list[isite][cn][i_coordnbs].minimum_geometries()
                csms = [ce_dict['symmetry_measure'] for ce_symbol, ce_dict in mingeoms
                        if ce_dict['symmetry_measure'] <= self._cn_self_mean_csm_estimator[1]['max_csm']]
                mean_csm = self._mean_csm(csms)
                if mean_csm is None:
                    continue
                self._mean_csms_isite[(cn, i_coordnbs)] = mean_csm
                self._mingeoms_isite[(cn, i_coordnbs)] = mingeoms

    def _setup_mean_csms_deltas(self):
        self._mean_csms_deltas_isite = {}
        for cn_map1 in self._mean_csms_isite:
            cn1 = cn_map1[0]
            self._mean_csms_deltas_isite[cn_map1] = {}
            for cn_map2 in self._mean_csms_isite:
                cn2 = cn_map2[0]
                if cn1 < cn2:
                    self._mean_csms_deltas_isite[cn_map1][cn_map2] = (self._mean_csms_isite[cn_map2] -
                                                                      self._mean_csms_isite[cn_map1])
                if cn1 == cn2 and cn_map1[1] != cn_map2[1]:
                    self._mean_csms_deltas_isite[cn_map1][cn_map2] = (self._mean_csms_isite[cn_map2] -
                                                                      self._mean_csms_isite[cn_map1])

    def get_site_neighbors(self, site):
        #TODO: do this one ...
        [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)

    def get_site_coordination_environment(self, site, return_maps=False):
        #TODO: do this one ...
        envs_ce = self.get_site_coordination_environments_fractions(site)
        return envs_ce[0]

    def get_site_coordination_environments(self, site, return_maps=False):
        #TODO: do this one ...
        envs_ce = self.get_site_coordination_environments_fractions(site)
        return envs_ce

    def __eq__(self, other):
        return (self.__class__.__name__ == other.__class__.__name__ and
                self._additional_condition == other._additional_condition and
                self._mean_csm_estimator == other._mean_csm_estimator and
                self._cn_self_mean_csm_estimator == other._cn_self_mean_csm_estimator and
                self._cn_self_mean_csm_estimator_cn_specifics == other._cn_self_mean_csm_estimator_cn_specifics and
                self._cn_self_mean_csm_estimator_ce_specifics == other._cn_self_mean_csm_estimator_ce_specifics and
                self._cn_delta_mean_csm_estimator == other._cn_delta_mean_csm_estimator and
                self._cn_delta_mean_csm_estimator_cn_specifics == other._cn_delta_mean_csm_estimator_cn_specifics and
                self._cn_delta_mean_csm_estimator_ce_specifics == other._cn_delta_mean_csm_estimator_ce_specifics and
                self._ce_estimator == other._ce_estimator)

    def as_dict(self):
        """
        Bson-serializable dict representation of the ComplexCSMBasedChemenvStrategy object.
        :return: Bson-serializable dict representation of the ComplexCSMBasedChemenvStrategy object.
        """
        dd = {'@module': self.__class__.__module__,
              '@class': self.__class__.__name__,
              'additional_condition': self._additional_condition,
              'mean_csm_estimator': self._mean_csm_estimator,
              'cn_self_mean_csm_estimator': self._cn_self_mean_csm_estimator,
              'cn_self_mean_csm_estimator_cn_specifics': self._cn_self_mean_csm_estimator_cn_specifics,
              'cn_self_mean_csm_estimator_ce_specifics': self._cn_self_mean_csm_estimator_ce_specifics,
              'cn_delta_mean_csm_estimator': self._cn_delta_mean_csm_estimator,
              'cn_delta_mean_csm_estimator_cn_specifics': self._cn_delta_mean_csm_estimator_cn_specifics,
              'cn_delta_mean_csm_estimator_ce_specifics': self._cn_delta_mean_csm_estimator_ce_specifics,
              'ce_estimator': self._ce_estimator}
        return dd

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the ComplexCSMBasedChemenvStrategy object from a dict representation of the
        ComplexCSMBasedChemenvStrategy object created using the as_dict method.
        :param d: dict representation of the ComplexCSMBasedChemenvStrategy object
        :return: ComplexCSMBasedChemenvStrategy object
        """
        if 'additional_condition' in d:
            return cls(additional_condition=d['additional_condition'],
                       mean_csm_estimator=d['mean_csm_estimator'],
                       cn_self_mean_csm_estimator=d['cn_self_mean_csm_estimator'],
                       cn_self_mean_csm_estimator_cn_specifics=d['cn_self_mean_csm_estimator_cn_specifics'],
                       cn_self_mean_csm_estimator_ce_specifics=d['cn_self_mean_csm_estimator_ce_specifics'],
                       cn_delta_mean_csm_estimator=d['cn_delta_mean_csm_estimator'],
                       cn_delta_mean_csm_estimator_cn_specifics=d['cn_delta_mean_csm_estimator_cn_specifics'],
                       cn_delta_mean_csm_estimator_ce_specifics=d['cn_delta_mean_csm_estimator_ce_specifics'],
                       ce_estimator=d['ce_estimator'])
        else:
            return cls()