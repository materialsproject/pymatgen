# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


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
from monty.json import MSONable
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp
from pymatgen.core.sites import PeriodicSite
import numpy as np
from scipy.stats import gmean
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import UNCLEAR_ENVIRONMENT_SYMBOL
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import get_lower_and_upper_f
from pymatgen.analysis.chemenv.utils.func_utils import CSMFiniteRatioFunction
from pymatgen.analysis.chemenv.utils.func_utils import CSMInfiniteRatioFunction
from pymatgen.analysis.chemenv.utils.func_utils import DeltaCSMRatioFunction
from pymatgen.analysis.chemenv.utils.func_utils import RatioFunction
from pymatgen.analysis.chemenv.utils.chemenv_errors import EquivalentSiteSearchError
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.utils.defs_utils import AdditionalConditions

from pymatgen.analysis.chemenv.coordination_environments.voronoi import DetailedVoronoiContainer
from collections import OrderedDict


module_dir = os.path.dirname(os.path.abspath(__file__))

MPSYMBOL_TO_CN = AllCoordinationGeometries().get_symbol_cn_mapping()
ALLCG = AllCoordinationGeometries()


class StrategyOption(MSONable, metaclass=abc.ABCMeta):

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
        if str(int(integer)) != str(integer):
            raise ValueError("Additional condition {} is not an integer".format(str(integer)))
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


class AbstractChemenvStrategy(MSONable, metaclass=abc.ABCMeta):
    """
    Class used to define a Chemenv strategy for the neighbors and coordination environment to be applied to a
    StructureEnvironments object
    """
    AC = AdditionalConditions()
    STRATEGY_OPTIONS = OrderedDict()
    STRATEGY_DESCRIPTION = None
    STRATEGY_INFO_FIELDS = []
    DEFAULT_SYMMETRY_MEASURE_TYPE = 'csm_wcs_ctwcc'

    def __init__(self, structure_environments=None, symmetry_measure_type=DEFAULT_SYMMETRY_MEASURE_TYPE):
        """
        Abstract constructor for the all chemenv strategies.
        :param structure_environments: StructureEnvironments object containing all the information on the
        coordination of the sites in a structure
        """
        self.structure_environments = None
        if structure_environments is not None:
            self.set_structure_environments(structure_environments)
        self._symmetry_measure_type = symmetry_measure_type

    @property
    def symmetry_measure_type(self):
        return self._symmetry_measure_type

    def set_structure_environments(self, structure_environments):
        self.structure_environments = structure_environments
        if not isinstance(self.structure_environments.voronoi, DetailedVoronoiContainer):
            raise ValueError('Voronoi Container not of type "DetailedVoronoiContainer"')
        self.prepare_symmetries()

    def prepare_symmetries(self):
        try:
            self.spg_analyzer = SpacegroupAnalyzer(self.structure_environments.structure)
            self.symops = self.spg_analyzer.get_symmetry_operations()
        except:
            self.symops = []

    def equivalent_site_index_and_transform(self, psite):
        # Get the index of the site in the unit cell of which the PeriodicSite psite is a replica.
        try:
            isite = self.structure_environments.structure.index(psite)
        except ValueError:
            try:
                uc_psite = psite.to_unit_cell()
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
        equivsite = self.structure_environments.structure[self.structure_environments.sites_map[isite]].to_unit_cell()
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
        raise NotImplementedError()

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
        raise NotImplementedError()

    @abc.abstractmethod
    def get_site_coordination_environments(self, site):
        """
        Applies the strategy to the structure_environments object in order to define the coordination environment of
        a given site.
        :param site: Site for which the coordination environment is looked for
        :return: The coordination environment of the site. For complex strategies, where one allows multiple
        solutions, this can return a list of coordination environments for the site
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_site_coordination_environments_fractions(self, site, isite=None, dequivsite=None, dthissite=None,
                                                     mysym=None, ordered=True, min_fraction=0.0, return_maps=True,
                                                     return_strategy_dict_info=False):
        """
        Applies the strategy to the structure_environments object in order to define the coordination environment of
        a given site.
        :param site: Site for which the coordination environment is looked for
        :return: The coordination environment of the site. For complex strategies, where one allows multiple
        solutions, this can return a list of coordination environments for the site
        """
        raise NotImplementedError()

    def get_site_ce_fractions_and_neighbors(self, site, full_ce_info=False, strategy_info=False):
        """
        Applies the strategy to the structure_environments object in order to get coordination environments, their
        fraction, csm, geometry_info, and neighbors
        :param site: Site for which the above information is seeked
        :return: The list of neighbors of the site. For complex strategies, where one allows multiple solutions, this
        can return a list of list of neighbors
        """
        [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        geoms_and_maps_list = self.get_site_coordination_environments_fractions(site=site, isite=isite,
                                                                                dequivsite=dequivsite,
                                                                                dthissite=dthissite, mysym=mysym,
                                                                                return_maps=True,
                                                                                return_strategy_dict_info=True)
        if geoms_and_maps_list is None:
            return None
        site_nbs_sets = self.structure_environments.neighbors_sets[isite]
        ce_and_neighbors = []
        for fractions_dict in geoms_and_maps_list:
            ce_map = fractions_dict['ce_map']
            ce_nb_set = site_nbs_sets[ce_map[0]][ce_map[1]]
            neighbors = [{'site': nb_site_and_index['site'],
                          'index': nb_site_and_index['index']}
                         for nb_site_and_index in ce_nb_set.neighb_sites_and_indices]
            fractions_dict['neighbors'] = neighbors
            ce_and_neighbors.append(fractions_dict)
        return ce_and_neighbors

    def set_option(self, option_name, option_value):
        self.__setattr__(option_name, option_value)

    def setup_options(self, all_options_dict):
        for option_name, option_value in all_options_dict.items():
            self.set_option(option_name, option_value)

    @abc.abstractmethod
    def __eq__(self, other):
        """
        Equality method that should be implemented for any strategy
        :param other: strategy to be compared with the current one
        :return:
        """
        raise NotImplementedError()

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

    @abc.abstractmethod
    def as_dict(self):
        """
        Bson-serializable dict representation of the SimplestChemenvStrategy object.
        :return: Bson-serializable dict representation of the SimplestChemenvStrategy object.
        """
        raise NotImplementedError()

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the SimpleAbundanceChemenvStrategy object from a dict representation of the
        SimpleAbundanceChemenvStrategy object created using the as_dict method.
        :param d: dict representation of the SimpleAbundanceChemenvStrategy object
        :return: StructureEnvironments object
        """
        raise NotImplementedError()


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
    STRATEGY_OPTIONS = OrderedDict()
    STRATEGY_OPTIONS['distance_cutoff'] =  {'type': DistanceCutoffFloat, 'internal': '_distance_cutoff',
                                            'default': DEFAULT_DISTANCE_CUTOFF}
    STRATEGY_OPTIONS['angle_cutoff'] = {'type': AngleCutoffFloat, 'internal': '_angle_cutoff',
                                        'default': DEFAULT_ANGLE_CUTOFF}
    STRATEGY_OPTIONS['additional_condition'] = {'type': AdditionalConditionInt,
                                                'internal': '_additional_condition',
                                                'default': DEFAULT_ADDITIONAL_CONDITION}
    STRATEGY_OPTIONS['continuous_symmetry_measure_cutoff'] = {'type': CSMFloat,
                                                              'internal': '_continuous_symmetry_measure_cutoff',
                                                              'default': DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF}

    STRATEGY_DESCRIPTION = '    Simplest ChemenvStrategy using fixed angle and distance parameters \n' \
                           '    for the definition of neighbors in the Voronoi approach. \n' \
                           '    The coordination environment is then given as the one with the \n' \
                           '    lowest continuous symmetry measure.'

    def __init__(self, structure_environments=None, distance_cutoff=DEFAULT_DISTANCE_CUTOFF,
                 angle_cutoff=DEFAULT_ANGLE_CUTOFF, additional_condition=DEFAULT_ADDITIONAL_CONDITION,
                 continuous_symmetry_measure_cutoff=DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF,
                 symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE):
        """
        Constructor for this SimplestChemenvStrategy.
        :param distance_cutoff: Distance cutoff used
        :param angle_cutoff: Angle cutoff used
        """
        AbstractChemenvStrategy.__init__(self, structure_environments, symmetry_measure_type=symmetry_measure_type)
        self.distance_cutoff = distance_cutoff
        self.angle_cutoff = angle_cutoff
        self.additional_condition = additional_condition
        self.continuous_symmetry_measure_cutoff = continuous_symmetry_measure_cutoff

    @property
    def uniquely_determines_coordination_environments(self):
        return True

    @property
    def distance_cutoff(self):
        return self._distance_cutoff

    @distance_cutoff.setter
    def distance_cutoff(self, distance_cutoff):
        self._distance_cutoff = DistanceCutoffFloat(distance_cutoff)

    @property
    def angle_cutoff(self):
        return self._angle_cutoff

    @angle_cutoff.setter
    def angle_cutoff(self, angle_cutoff):
        self._angle_cutoff = AngleCutoffFloat(angle_cutoff)

    @property
    def additional_condition(self):
        return self._additional_condition

    @additional_condition.setter
    def additional_condition(self, additional_condition):
        self._additional_condition = AdditionalConditionInt(additional_condition)

    @property
    def continuous_symmetry_measure_cutoff(self):
        return self._continuous_symmetry_measure_cutoff

    @continuous_symmetry_measure_cutoff.setter
    def continuous_symmetry_measure_cutoff(self, continuous_symmetry_measure_cutoff):
        self._continuous_symmetry_measure_cutoff = CSMFloat(continuous_symmetry_measure_cutoff)

    def get_site_neighbors(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None):#, neighbors_map=None):
        #if neighbors_map is not None:
        #    return self.structure_environments.voronoi.get_neighbors(isite=isite, neighbors_map=neighbors_map)
        if isite is None:
            [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)

        ce, cn_map = self.get_site_coordination_environment(site=site, isite=isite,
                                                            dequivsite=dequivsite, dthissite=dthissite, mysym=mysym,
                                                            return_map=True)

        nb_set = self.structure_environments.neighbors_sets[isite][cn_map[0]][cn_map[1]]
        eqsite_ps = nb_set.neighb_sites

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
        neighbors_normalized_distances = self.structure_environments.voronoi.neighbors_normalized_distances[isite]
        neighbors_normalized_angles = self.structure_environments.voronoi.neighbors_normalized_angles[isite]
        idist = None
        for iwd, wd in enumerate(neighbors_normalized_distances):
            if self.distance_cutoff >= wd['min']:
                idist = iwd
            else:
                break
        iang = None
        for iwa, wa in enumerate(neighbors_normalized_angles):
            if self.angle_cutoff <= wa['max']:
                iang = iwa
            else:
                break
        if idist is None or iang is None:
            raise ValueError('Distance or angle parameter not found ...')

        my_cn = None
        my_inb_set = None
        found = False
        for cn, nb_sets in self.structure_environments.neighbors_sets[isite].items():
            for inb_set, nb_set in enumerate(nb_sets):
                sources = [src for src in nb_set.sources
                           if src['origin'] == 'dist_ang_ac_voronoi' and src['ac'] == self.additional_condition]
                for src in sources:
                    if src['idp'] == idist and src['iap'] == iang:
                        my_cn = cn
                        my_inb_set = inb_set
                        found = True
                        break
                if found:
                    break
            if found:
                break

        if not found:
            return None

        cn_map = (my_cn, my_inb_set)
        ce = self.structure_environments.ce_list[self.structure_environments.sites_map[isite]][cn_map[0]][cn_map[1]]
        if ce is None:
            return None
        coord_geoms = ce.coord_geoms
        if return_map:
            if coord_geoms is None:
                return cn_map[0], cn_map
            return (ce.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type), cn_map)
        else:
            if coord_geoms is None:
                return cn_map[0]
            return ce.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type)

    def get_site_coordination_environments_fractions(self, site, isite=None, dequivsite=None, dthissite=None,
                                                     mysym=None, ordered=True, min_fraction=0.0, return_maps=True,
                                                     return_strategy_dict_info=False):
        if isite is None or dequivsite is None or dthissite is None or mysym is None:
            [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        site_nb_sets = self.structure_environments.neighbors_sets[isite]
        if site_nb_sets is None:
            return None
        ce_and_map = self.get_site_coordination_environment(site=site, isite=isite, dequivsite=dequivsite,
                                                            dthissite=dthissite, mysym=mysym,
                                                            return_map=True)
        if ce_and_map is None:
            return None
        ce, ce_map = ce_and_map
        if ce is None:
            ce_dict = {'ce_symbol': 'UNKNOWN:{:d}'.format(ce_map[0]), 'ce_dict': None, 'ce_fraction': 1.0}
        else:
            ce_dict = {'ce_symbol': ce[0], 'ce_dict': ce[1], 'ce_fraction': 1.0}
        if return_maps:
            ce_dict['ce_map'] = ce_map
        if return_strategy_dict_info:
            ce_dict['strategy_info'] = {}
        fractions_info_list = [ce_dict]
        return fractions_info_list

    def get_site_coordination_environments(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None,
                                           return_maps=False):
        return [self.get_site_coordination_environment(site=site, isite=isite, dequivsite=dequivsite,
                                                       dthissite=dthissite, mysym=mysym, return_map=return_maps)]

    def add_strategy_visualization_to_subplot(self, subplot, visualization_options=None, plot_type=None):
        subplot.plot(self._distance_cutoff, self._angle_cutoff, 'o', mec=None, mfc='w', markersize=12)
        subplot.plot(self._distance_cutoff, self._angle_cutoff, 'x', linewidth=2, markersize=12)

    def __eq__(self, other):
        return (self.__class__.__name__ == other.__class__.__name__ and
                self._distance_cutoff == other._distance_cutoff and self._angle_cutoff == other._angle_cutoff and
                self._additional_condition == other._additional_condition and
                self._continuous_symmetry_measure_cutoff == other._continuous_symmetry_measure_cutoff and
                self.symmetry_measure_type == other.symmetry_measure_type)

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
                "continuous_symmetry_measure_cutoff": float(self._continuous_symmetry_measure_cutoff),
                "symmetry_measure_type": self._symmetry_measure_type}

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
                   continuous_symmetry_measure_cutoff=d["continuous_symmetry_measure_cutoff"],
                   symmetry_measure_type=d["symmetry_measure_type"])


class SimpleAbundanceChemenvStrategy(AbstractChemenvStrategy):
    """
    Simple ChemenvStrategy using the neighbors that are the most "abundant" in the grid of angle and distance
    parameters for the definition of neighbors in the Voronoi approach.
    The coordination environment is then given as the one with the lowest continuous symmetry measure
    """

    DEFAULT_MAX_DIST = 2.0
    DEFAULT_ADDITIONAL_CONDITION = AbstractChemenvStrategy.AC.ONLY_ACB
    STRATEGY_OPTIONS = OrderedDict()
    STRATEGY_OPTIONS['additional_condition'] = {'type': AdditionalConditionInt,
                                                'internal': '_additional_condition',
                                                'default': DEFAULT_ADDITIONAL_CONDITION}
    STRATEGY_OPTIONS['surface_calculation_type'] = {}
    STRATEGY_DESCRIPTION = '    Simple Abundance ChemenvStrategy using the most "abundant" neighbors map \n' \
                           '    for the definition of neighbors in the Voronoi approach. \n' \
                           '    The coordination environment is then given as the one with the \n' \
                           '    lowest continuous symmetry measure.'

    def __init__(self, structure_environments=None,
                 additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB,
                 symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE):
        """
        Constructor for the SimpleAbundanceChemenvStrategy.
        :param structure_environments: StructureEnvironments object containing all the information on the
        coordination of the sites in a structure
        """
        raise NotImplementedError('SimpleAbundanceChemenvStrategy not yet implemented')
        AbstractChemenvStrategy.__init__(self, structure_environments, symmetry_measure_type=symmetry_measure_type)
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
            return coord_geoms.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type), cn_map
        else:
            if coord_geoms is None:
                return cn_map[0]
            return coord_geoms.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type)

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
    DEFAULT_TARGET_ENVIRONMENTS = ['O:6']

    def __init__(self, structure_environments=None, truncate_dist_ang=True,
                 additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB,
                 max_nabundant=5, target_environments=DEFAULT_TARGET_ENVIRONMENTS, target_penalty_type='max_csm',
                 max_csm=5.0, symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE):
        raise NotImplementedError('TargettedPenaltiedAbundanceChemenvStrategy not yet implemented')
        SimpleAbundanceChemenvStrategy.__init__(self, structure_environments,
                                                additional_condition=additional_condition,
                                                symmetry_measure_type=symmetry_measure_type)
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
            return chemical_environments.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type), cn_map
        else:
            if chemical_environments.coord_geoms is None:
                return cn_map[0]
            return chemical_environments.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type)

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
                          [self.structure_environments.sites_map[isite]]
                          [mymap[0]][mymap[1]].minimum_geometry(symmetry_measure_type=self._symmetry_measure_type))
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


class NbSetWeight(MSONable, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def as_dict(self):
        """
        A JSON serializable dict representation of this neighbors set weight.
        """
        pass

    @abc.abstractmethod
    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        pass


class AngleNbSetWeight(NbSetWeight):

    SHORT_NAME = 'AngleWeight'

    def __init__(self, aa=1.0):
        self.aa = aa
        if self.aa == 1.0:
            self.aw = self.angle_sum
        else:
            self.aw = self.angle_sumn

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        return self.aw(nb_set=nb_set)

    def angle_sum(self, nb_set):
        return np.sum(nb_set.angles) / (4.0 * np.pi)

    def angle_sumn(self, nb_set):
        return np.power(self.angle_sum(nb_set=nb_set), self.aa)

    def __eq__(self, other):
        return self.aa == other.aa

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "aa": self.aa
                }

    @classmethod
    def from_dict(cls, dd):
        return cls(aa=dd['aa'])


class NormalizedAngleDistanceNbSetWeight(NbSetWeight):

    SHORT_NAME = 'NormAngleDistWeight'

    def __init__(self, average_type, aa, bb):
        self.average_type = average_type
        if self.average_type == 'geometric':
            self.eval = self.gweight
        elif self.average_type == 'arithmetic':
            self.eval = self.aweight
        else:
            raise ValueError('Average type is "{}" while it should be '
                             '"geometric" or "arithmetic"'.format(average_type))
        self.aa = aa
        self.bb = bb
        if self.aa == 0:
            if self.bb == 1:
                self.fda = self.invdist
            elif self.bb == 0:
                raise ValueError('Both exponents are 0.')
            else:
                self.fda = self.invndist
        elif self.bb == 0:
            if self.aa == 1:
                self.fda = self.ang
            else:
                self.fda = self.angn
        else:
            if self.aa == 1:
                if self.bb == 1:
                    self.fda = self.anginvdist
                else:
                    self.fda = self.anginvndist
            else:
                if self.bb == 1:
                    self.fda = self.angninvdist
                else:
                    self.fda = self.angninvndist

    def __eq__(self, other):
        return self.average_type == other.average_type and self.aa == other.aa and self.bb == other.bb

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "average_type": self.average_type,
                "aa": self.aa,
                "bb": self.bb
                }

    @classmethod
    def from_dict(cls, dd):
        return cls(average_type=dd['average_type'], aa=dd['aa'], bb=dd['bb'])

    def invdist(self, nb_set):
        return [1.0 / dist for dist in nb_set.normalized_distances]

    def invndist(self, nb_set):
        return [1.0 / dist**self.bb for dist in nb_set.normalized_distances]

    def ang(self, nb_set):
        return nb_set.normalized_angles

    def angn(self, nb_set):
        return [ang**self.aa for ang in nb_set.normalized_angles]

    def anginvdist(self, nb_set):
        nangles = nb_set.normalized_angles
        return [nangles[ii] / dist for ii, dist in enumerate(nb_set.normalized_distances)]

    def anginvndist(self, nb_set):
        nangles = nb_set.normalized_angles
        return [nangles[ii] / dist**self.bb for ii, dist in enumerate(nb_set.normalized_distances)]

    def angninvdist(self, nb_set):
        nangles = nb_set.normalized_angles
        return [nangles[ii]**self.aa / dist for ii, dist in enumerate(nb_set.normalized_distances)]

    def angninvndist(self, nb_set):
        nangles = nb_set.normalized_angles
        return [nangles[ii]**self.aa / dist**self.bb for ii, dist in enumerate(nb_set.normalized_distances)]

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        fda_list = self.fda(nb_set=nb_set)
        return self.eval(fda_list=fda_list)

    def gweight(self, fda_list):
        return gmean(fda_list)

    def aweight(self, fda_list):
        return np.mean(fda_list)


def get_effective_csm(nb_set, cn_map, structure_environments, additional_info,
                      symmetry_measure_type, max_effective_csm, effective_csm_estimator_ratio_function):
    try:
        effective_csm = additional_info['effective_csms'][nb_set.isite][cn_map]
    except KeyError:
        site_ce_list = structure_environments.ce_list[nb_set.isite]
        site_chemenv = site_ce_list[cn_map[0]][cn_map[1]]
        if site_chemenv is None:
            effective_csm = 100.0
        else:
            mingeoms = site_chemenv.minimum_geometries(symmetry_measure_type=symmetry_measure_type,
                                                   max_csm=max_effective_csm)
            if len(mingeoms) == 0:
                effective_csm = 100.0
            else:
                csms = [ce_dict['other_symmetry_measures'][symmetry_measure_type] for mp_symbol, ce_dict in mingeoms
                        if ce_dict['other_symmetry_measures'][symmetry_measure_type] <= max_effective_csm]
                effective_csm = effective_csm_estimator_ratio_function.mean_estimator(csms)
        set_info(additional_info=additional_info, field='effective_csms',
                 isite=nb_set.isite, cn_map=cn_map, value=effective_csm)
    return effective_csm


def set_info(additional_info, field, isite, cn_map, value):
    try:
        additional_info[field][isite][cn_map] = value
    except KeyError:
        try:
            additional_info[field][isite] = {cn_map: value}
        except KeyError:
            additional_info[field] = {isite: {cn_map: value}}


class SelfCSMNbSetWeight(NbSetWeight):

    SHORT_NAME = 'SelfCSMWeight'

    DEFAULT_EFFECTIVE_CSM_ESTIMATOR = {'function': 'power2_inverse_decreasing',
                                       'options': {'max_csm': 8.0}}
    DEFAULT_WEIGHT_ESTIMATOR = {'function': 'power2_decreasing_exp',
                                'options': {'max_csm': 8.0,
                                            'alpha': 1.0}}
    DEFAULT_SYMMETRY_MEASURE_TYPE = 'csm_wcs_ctwcc'

    def __init__(self, effective_csm_estimator=DEFAULT_EFFECTIVE_CSM_ESTIMATOR,
                 weight_estimator=DEFAULT_WEIGHT_ESTIMATOR,
                 symmetry_measure_type=DEFAULT_SYMMETRY_MEASURE_TYPE):
        self.effective_csm_estimator = effective_csm_estimator
        self.effective_csm_estimator_rf = CSMInfiniteRatioFunction.from_dict(effective_csm_estimator)
        self.weight_estimator = weight_estimator
        self.weight_estimator_rf =  CSMFiniteRatioFunction.from_dict(weight_estimator)
        self.symmetry_measure_type = symmetry_measure_type
        self.max_effective_csm = self.effective_csm_estimator['options']['max_csm']

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        effective_csm = get_effective_csm(nb_set=nb_set, cn_map=cn_map,
                                          structure_environments=structure_environments,
                                          additional_info=additional_info,
                                          symmetry_measure_type=self.symmetry_measure_type,
                                          max_effective_csm=self.max_effective_csm,
                                          effective_csm_estimator_ratio_function=self.effective_csm_estimator_rf)
        weight = self.weight_estimator_rf.evaluate(effective_csm)
        set_info(additional_info=additional_info, field='self_csms_weights', isite=nb_set.isite,
                 cn_map=cn_map, value=weight)
        return weight

    def __eq__(self, other):
        return (self.effective_csm_estimator == other.effective_csm_estimator and
                self.weight_estimator == other.weight_estimator and
                self.symmetry_measure_type == other.symmetry_measure_type)

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "effective_csm_estimator": self.effective_csm_estimator,
                "weight_estimator": self.weight_estimator,
                "symmetry_measure_type": self.symmetry_measure_type
                }

    @classmethod
    def from_dict(cls, dd):
        return cls(effective_csm_estimator=dd['effective_csm_estimator'],
                   weight_estimator=dd['weight_estimator'],
                   symmetry_measure_type=dd['symmetry_measure_type'])


class DeltaCSMNbSetWeight(NbSetWeight):

    SHORT_NAME = 'DeltaCSMWeight'

    DEFAULT_EFFECTIVE_CSM_ESTIMATOR = {'function': 'power2_inverse_decreasing',
                                       'options': {'max_csm': 8.0}}
    DEFAULT_SYMMETRY_MEASURE_TYPE = 'csm_wcs_ctwcc'
    DEFAULT_WEIGHT_ESTIMATOR = {'function': 'smootherstep',
                                'options': {'delta_csm_min': 0.5,
                                            'delta_csm_max': 3.0}}

    def __init__(self, effective_csm_estimator=DEFAULT_EFFECTIVE_CSM_ESTIMATOR,
                 weight_estimator=DEFAULT_WEIGHT_ESTIMATOR,
                 delta_cn_weight_estimators=None,
                 symmetry_measure_type=DEFAULT_SYMMETRY_MEASURE_TYPE):
        self.effective_csm_estimator = effective_csm_estimator
        self.effective_csm_estimator_rf = CSMInfiniteRatioFunction.from_dict(effective_csm_estimator)
        self.weight_estimator = weight_estimator
        if self.weight_estimator is not None:
            self.weight_estimator_rf = DeltaCSMRatioFunction.from_dict(weight_estimator)
        self.delta_cn_weight_estimators = delta_cn_weight_estimators
        self.delta_cn_weight_estimators_rfs = {}
        if delta_cn_weight_estimators is not None:
            for delta_cn, dcn_w_estimator in delta_cn_weight_estimators.items():
                self.delta_cn_weight_estimators_rfs[delta_cn] = DeltaCSMRatioFunction.from_dict(dcn_w_estimator)
        self.symmetry_measure_type = symmetry_measure_type
        self.max_effective_csm = self.effective_csm_estimator['options']['max_csm']

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        effcsm = get_effective_csm(nb_set=nb_set, cn_map=cn_map,
                                   structure_environments=structure_environments,
                                   additional_info=additional_info,
                                   symmetry_measure_type=self.symmetry_measure_type,
                                   max_effective_csm=self.max_effective_csm,
                                   effective_csm_estimator_ratio_function=self.effective_csm_estimator_rf)
        cn = cn_map[0]
        inb_set = cn_map[1]
        isite = nb_set.isite
        delta_csm = None
        delta_csm_cn_map2 = None
        nb_set_weight = 1.0
        for cn2, nb_sets in structure_environments.neighbors_sets[isite].items():
            if cn2 < cn:
                continue
            for inb_set2, nb_set2 in enumerate(nb_sets):
                if cn == cn2 and inb_set == inb_set:
                    continue
                effcsm2 = get_effective_csm(nb_set=nb_set2, cn_map=(cn2, inb_set2),
                                            structure_environments=structure_environments,
                                            additional_info=additional_info,
                                            symmetry_measure_type=self.symmetry_measure_type,
                                            max_effective_csm=self.max_effective_csm,
                                            effective_csm_estimator_ratio_function=self.effective_csm_estimator_rf)
                this_delta_csm = effcsm2 - effcsm
                if cn2 == cn:
                    if this_delta_csm < 0.0:
                        set_info(additional_info=additional_info, field='delta_csms', isite=isite,
                                 cn_map=cn_map, value=this_delta_csm)
                        set_info(additional_info=additional_info, field='delta_csms_weights', isite=isite,
                                 cn_map=cn_map, value=0.0)
                        set_info(additional_info=additional_info, field='delta_csms_cn_map2', isite=isite,
                                 cn_map=cn_map, value=(cn2, inb_set2))
                        return 0.0
                else:
                    dcn = cn2 - cn
                    if dcn in self.delta_cn_weight_estimators_rfs:
                        this_delta_csm_weight = self.delta_cn_weight_estimators_rfs[dcn].evaluate(this_delta_csm)
                    else:
                        this_delta_csm_weight = self.weight_estimator_rf.evaluate(this_delta_csm)
                    if this_delta_csm_weight < nb_set_weight:
                        delta_csm = this_delta_csm
                        delta_csm_cn_map2 = (cn2, inb_set2)
                        nb_set_weight = this_delta_csm_weight
        set_info(additional_info=additional_info, field='delta_csms', isite=isite,
                 cn_map=cn_map, value=delta_csm)
        set_info(additional_info=additional_info, field='delta_csms_weights', isite=isite,
                 cn_map=cn_map, value=nb_set_weight)
        set_info(additional_info=additional_info, field='delta_csms_cn_map2', isite=isite,
                 cn_map=cn_map, value=delta_csm_cn_map2)
        return nb_set_weight

    def __eq__(self, other):
        return (self.effective_csm_estimator == other.effective_csm_estimator and
                self.weight_estimator == other.weight_estimator and
                self.delta_cn_weight_estimators == other.delta_cn_weight_estimators and
                self.symmetry_measure_type == other.symmetry_measure_type)

    def __ne__(self, other):
        return not self == other

    @classmethod
    def delta_cn_specifics(cls, delta_csm_mins=None, delta_csm_maxs=None, function='smootherstep',
                           symmetry_measure_type='csm_wcs_ctwcc',
                           effective_csm_estimator=DEFAULT_EFFECTIVE_CSM_ESTIMATOR):
        if delta_csm_mins is None or delta_csm_maxs is None:
            delta_cn_weight_estimators = {dcn: {'function': function,
                                                'options': {'delta_csm_min': 0.25+dcn*0.25,
                                                            'delta_csm_max': 5.0+dcn*0.25}} for dcn in range(1, 13)}
        else:
            delta_cn_weight_estimators = {dcn: {'function': function,
                                                'options': {'delta_csm_min': delta_csm_mins[dcn-1],
                                                            'delta_csm_max': delta_csm_maxs[dcn-1]}}
                                          for dcn in range(1, 13)}
        return cls(effective_csm_estimator=effective_csm_estimator,
                   weight_estimator={'function': function,
                                     'options': {'delta_csm_min': delta_cn_weight_estimators[12]
                                                 ['options']['delta_csm_min'],
                                                 'delta_csm_max': delta_cn_weight_estimators[12]
                                                 ['options']['delta_csm_max']}},
                   delta_cn_weight_estimators=delta_cn_weight_estimators,
                   symmetry_measure_type=symmetry_measure_type)

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "effective_csm_estimator": self.effective_csm_estimator,
                "weight_estimator": self.weight_estimator,
                "delta_cn_weight_estimators": self.delta_cn_weight_estimators,
                "symmetry_measure_type": self.symmetry_measure_type
                }

    @classmethod
    def from_dict(cls, dd):
        return cls(effective_csm_estimator=dd['effective_csm_estimator'],
                   weight_estimator=dd['weight_estimator'],
                   delta_cn_weight_estimators={int(dcn): dcn_estimator
                                               for dcn, dcn_estimator in dd['delta_cn_weight_estimators'].items()}
                   if ('delta_cn_weight_estimators' in dd and dd['delta_cn_weight_estimators'] is not None) else None,
                   symmetry_measure_type=dd['symmetry_measure_type'])


class CNBiasNbSetWeight(NbSetWeight):

    SHORT_NAME = 'CNBiasWeight'

    def __init__(self, cn_weights, initialization_options):
        self.cn_weights = cn_weights
        self.initialization_options = initialization_options

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        return self.cn_weights[len(nb_set)]

    def __eq__(self, other):
        return (self.cn_weights == other.cn_weights and
                self.initialization_options == other.initialization_options)

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "cn_weights": {str(cn): cnw for cn, cnw in self.cn_weights.items()},
                "initialization_options": self.initialization_options,
                }

    @classmethod
    def from_dict(cls, dd):
        return cls(cn_weights={int(cn): cnw for cn, cnw in dd['cn_weights'].items()},
                   initialization_options=dd['initialization_options'])

    @classmethod
    def linearly_equidistant(cls, weight_cn1, weight_cn13):
        initialization_options = {'type': 'linearly_equidistant',
                                  'weight_cn1': weight_cn1,
                                  'weight_cn13': weight_cn13
                                  }
        dw = (weight_cn13 - weight_cn1) / 12.0
        cn_weights = {cn: weight_cn1 + (cn - 1) * dw for cn in range(1, 14)}
        return cls(cn_weights=cn_weights, initialization_options=initialization_options)

    @classmethod
    def geometrically_equidistant(cls, weight_cn1, weight_cn13):
        initialization_options = {'type': 'geometrically_equidistant',
                                  'weight_cn1': weight_cn1,
                                  'weight_cn13': weight_cn13
                                  }
        factor = np.power(float(weight_cn13) / weight_cn1, 1.0 / 12.0)
        cn_weights = {cn: weight_cn1 * np.power(factor, cn - 1)  for cn in range(1, 14)}
        return cls(cn_weights=cn_weights, initialization_options=initialization_options)

    @classmethod
    def explicit(cls, cn_weights):
        initialization_options = {'type': 'explicit'}
        if set(cn_weights.keys()) != set(range(1, 14)):
            raise ValueError('Weights should be provided for CN 1 to 13')
        return cls(cn_weights=cn_weights, initialization_options=initialization_options)

    @classmethod
    def from_description(cls, dd):
        if dd['type'] == 'linearly_equidistant':
            return cls.linearly_equidistant(weight_cn1=dd['weight_cn1'], weight_cn13=dd['weight_cn13'])
        elif dd['type'] == 'geometrically_equidistant':
            return cls.geometrically_equidistant(weight_cn1=dd['weight_cn1'], weight_cn13=dd['weight_cn13'])
        elif dd['type'] == 'explicit':
            return cls.explicit(cn_weights=dd['cn_weights'])


class DistanceAngleAreaNbSetWeight(NbSetWeight):

    SHORT_NAME = 'DistAngleAreaWeight'

    AC = AdditionalConditions()
    DEFAULT_SURFACE_DEFINITION = {'type': 'standard_elliptic',
                                  'distance_bounds': {'lower': 1.2, 'upper': 1.8},
                                  'angle_bounds': {'lower': 0.1, 'upper': 0.8}}

    def __init__(self, weight_type='has_intersection', surface_definition=DEFAULT_SURFACE_DEFINITION,
                 nb_sets_from_hints='fallback_to_source', other_nb_sets='0_weight',
                 additional_condition=AC.ONLY_ACB, smoothstep_distance=None, smoothstep_angle=None):
        self.weight_type = weight_type
        if weight_type == 'has_intersection':
            self.area_weight = self.w_area_has_intersection
        elif weight_type == 'has_intersection_smoothstep':
            raise NotImplementedError()
            # self.area_weight = self.w_area_has_intersection_smoothstep
        else:
            raise ValueError('Weight type is "{}" while it should be "has_intersection"'.format(weight_type))
        self.surface_definition = surface_definition
        self.nb_sets_from_hints = nb_sets_from_hints
        self.other_nb_sets = other_nb_sets
        self.additional_condition = additional_condition
        self.smoothstep_distance = smoothstep_distance
        self.smoothstep_angle = smoothstep_angle
        if self.nb_sets_from_hints == 'fallback_to_source':
            if self.other_nb_sets == '0_weight':
                self.w_area_intersection_specific = self.w_area_intersection_nbsfh_fbs_onb0
            else:
                raise ValueError('Other nb_sets should be "0_weight"')
        else:
            raise ValueError('Nb_sets from hints should fallback to source')
        lower_and_upper_functions = get_lower_and_upper_f(surface_calculation_options=surface_definition)
        self.dmin = surface_definition['distance_bounds']['lower']
        self.dmax = surface_definition['distance_bounds']['upper']
        self.amin = surface_definition['angle_bounds']['lower']
        self.amax = surface_definition['angle_bounds']['upper']
        self.f_lower = lower_and_upper_functions['lower']
        self.f_upper = lower_and_upper_functions['upper']

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        return self.area_weight(nb_set=nb_set, structure_environments=structure_environments,
                                cn_map=cn_map, additional_info=additional_info)

    def w_area_has_intersection_smoothstep(self, nb_set, structure_environments,
                                           cn_map, additional_info):
        w_area =  self.w_area_intersection_specific(nb_set=nb_set, structure_environments=structure_environments,
                                                    cn_map=cn_map, additional_info=additional_info)
        if w_area > 0.0:
            if self.smoothstep_distance is not None:
                w_area = w_area
            if self.smoothstep_angle is not None:
                w_area = w_area
        return w_area

    def w_area_has_intersection(self, nb_set, structure_environments,
                                cn_map, additional_info):
        return self.w_area_intersection_specific(nb_set=nb_set, structure_environments=structure_environments,
                                                 cn_map=cn_map, additional_info=additional_info)

    def w_area_intersection_nbsfh_fbs_onb0(self, nb_set, structure_environments,
                                           cn_map, additional_info):
        dist_ang_sources = [src for src in nb_set.sources
                            if src['origin'] == 'dist_ang_ac_voronoi' and src['ac'] == self.additional_condition]
        if len(dist_ang_sources) > 0:
            for src in dist_ang_sources:
                d1 = src['dp_dict']['min']
                d2 = src['dp_dict']['next']
                a1 = src['ap_dict']['next']
                a2 = src['ap_dict']['max']
                if self.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2):
                    return 1.0
            return 0.0
        else:
            from_hints_sources = [src for src in nb_set.sources if src['origin'] == 'nb_set_hints']
            if len(from_hints_sources) == 0:
                return 0.0
            elif len(from_hints_sources) != 1:
                raise ValueError('Found multiple hints sources for nb_set')
            else:
                cn_map_src = from_hints_sources[0]['cn_map_source']
                nb_set_src = structure_environments.neighbors_sets[nb_set.isite][cn_map_src[0]][cn_map_src[1]]
                dist_ang_sources = [src for src in nb_set_src.sources
                                    if src['origin'] == 'dist_ang_ac_voronoi' and
                                    src['ac'] == self.additional_condition]
                if len(dist_ang_sources) == 0:
                    return 0.0
                for src in dist_ang_sources:
                    d1 = src['dp_dict']['min']
                    d2 = src['dp_dict']['next']
                    a1 = src['ap_dict']['next']
                    a2 = src['ap_dict']['max']
                    if self.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2):
                        return 1.0
                return 0.0

    def rectangle_crosses_area(self, d1, d2, a1, a2):
        # Case 1
        if d1 <= self.dmin and d2 <= self.dmin:
            return False
        # Case 6
        if d1 >= self.dmax and d2 >= self.dmax:
            return False
        # Case 2
        if d1 <= self.dmin and d2 <= self.dmax:
            ld2 = self.f_lower(d2)
            if a2 <= ld2 or a1 >= self.amax:
                return False
            return True
        # Case 3
        if d1 <= self.dmin and d2 >= self.dmax:
            if a2 <= self.amin or a1 >= self.amax:
                return False
            return True
        # Case 4
        if self.dmin <= d1 <= self.dmax and self.dmin <= d2 <= self.dmax:
            ld1 = self.f_lower(d1)
            ld2 = self.f_lower(d2)
            if a2 <= ld1 and a2 <= ld2:
                return False
            ud1 = self.f_upper(d1)
            ud2 = self.f_upper(d2)
            if a1 >= ud1 and a1 >= ud2:
                return False
            return True
        # Case 5
        if self.dmin <= d1 <= self.dmax and d2 >= self.dmax:
            ud1 = self.f_upper(d1)
            if a1 >= ud1 or a2 <= self.amin:
                return False
            return True
        raise ValueError('Should not reach this point!')

    def __eq__(self, other):
        return (self.weight_type == other.weight_type and
                self.surface_definition == other.surface_definition and
                self.nb_sets_from_hints == other.nb_sets_from_hints and
                self.other_nb_sets == other.other_nb_sets and
                self.additional_condition == other.additional_condition
                )

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "weight_type": self.weight_type,
                "surface_definition": self.surface_definition,
                "nb_sets_from_hints": self.nb_sets_from_hints,
                "other_nb_sets": self.other_nb_sets,
                "additional_condition": self.additional_condition}

    @classmethod
    def from_dict(cls, dd):
        return cls(weight_type=dd['weight_type'], surface_definition=dd['surface_definition'],
                   nb_sets_from_hints=dd['nb_sets_from_hints'], other_nb_sets=dd['other_nb_sets'],
                   additional_condition=dd['additional_condition'])


class DistancePlateauNbSetWeight(NbSetWeight):

    SHORT_NAME = 'DistancePlateauWeight'

    def __init__(self, distance_function=None, weight_function=None):
        if distance_function is None:
            self.distance_function = {'type': 'normalized_distance'}
        else:
            self.distance_function = distance_function
        if weight_function is None:
            self.weight_function = {'function': 'inverse_smootherstep', 'options': {'lower': 0.2, 'upper': 0.4}}
        else:
            self.weight_function = weight_function
        self.weight_rf = RatioFunction.from_dict(self.weight_function)

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        return self.weight_rf.eval(nb_set.distance_plateau())

    def __eq__(self, other):
        return self.__class__ == other.__class__

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "distance_function": self.distance_function,
                "weight_function": self.weight_function
                }

    @classmethod
    def from_dict(cls, dd):
        return cls(distance_function=dd['distance_function'], weight_function=dd['weight_function'])


class AnglePlateauNbSetWeight(NbSetWeight):

    SHORT_NAME = 'AnglePlateauWeight'

    def __init__(self, angle_function=None, weight_function=None):
        if angle_function is None:
            self.angle_function = {'type': 'normalized_angle'}
        else:
            self.angle_function = angle_function
        if weight_function is None:
            self.weight_function = {'function': 'inverse_smootherstep', 'options': {'lower': 0.05, 'upper': 0.15}}
        else:
            self.weight_function = weight_function
        self.weight_rf = RatioFunction.from_dict(self.weight_function)

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        return self.weight_rf.eval(nb_set.angle_plateau())

    def __eq__(self, other):
        return self.__class__ == other.__class__

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "angle_function": self.angle_function,
                "weight_function": self.weight_function
                }

    @classmethod
    def from_dict(cls, dd):
        return cls(angle_function=dd['angle_function'], weight_function=dd['weight_function'])


class DistanceNbSetWeight(NbSetWeight):

    SHORT_NAME = 'DistanceNbSetWeight'

    def __init__(self, weight_function=None, nbs_source='voronoi'):
        if weight_function is None:
            self.weight_function = {'function': 'smootherstep', 'options': {'lower': 1.2, 'upper': 1.3}}
        else:
            self.weight_function = weight_function
        self.weight_rf = RatioFunction.from_dict(self.weight_function)
        if nbs_source not in ['nb_sets', 'voronoi']:
            raise ValueError('"nbs_source" should be one of ["nb_sets", "voronoi"]')
        self.nbs_source = nbs_source

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        cn = cn_map[0]
        inb_set = cn_map[1]
        isite = nb_set.isite
        voronoi = structure_environments.voronoi.voronoi_list2[isite]
        if self.nbs_source == 'nb_sets':
            all_nbs_voro_indices = set()
            for cn2, nb_sets in structure_environments.neighbors_sets[isite].items():
                for inb_set2, nb_set2 in enumerate(nb_sets):
                    if cn == cn2 and inb_set == inb_set:
                        continue
                    all_nbs_voro_indices.update(nb_set2.site_voronoi_indices)
        elif self.nbs_source == "voronoi":
            all_nbs_voro_indices = set(range(len(voronoi)))
        else:
            raise ValueError('"nbs_source" should be one of ["nb_sets", "voronoi"]')
        all_nbs_indices_except_nb_set = all_nbs_voro_indices.difference(nb_set.site_voronoi_indices)
        normalized_distances = [voronoi[inb]['normalized_distance'] for inb in all_nbs_indices_except_nb_set]
        if len(normalized_distances) == 0:
            return 1.0
        return self.weight_rf.eval(min(normalized_distances))

    def __eq__(self, other):
        return self.__class__ == other.__class__

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "weight_function": self.weight_function,
                "nbs_source": self.nbs_source
                }

    @classmethod
    def from_dict(cls, dd):
        return cls(weight_function=dd['weight_function'], nbs_source=dd['nbs_source'])


class DeltaDistanceNbSetWeight(NbSetWeight):

    SHORT_NAME = 'DeltaDistanceNbSetWeight'

    def __init__(self, weight_function=None, nbs_source='voronoi'):
        if weight_function is None:
            self.weight_function = {'function': 'smootherstep', 'options': {'lower': 0.1, 'upper': 0.2}}
        else:
            self.weight_function = weight_function
        self.weight_rf = RatioFunction.from_dict(self.weight_function)
        if nbs_source not in ['nb_sets', 'voronoi']:
            raise ValueError('"nbs_source" should be one of ["nb_sets", "voronoi"]')
        self.nbs_source = nbs_source

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        cn = cn_map[0]
        inb_set = cn_map[1]
        isite = nb_set.isite
        voronoi = structure_environments.voronoi.voronoi_list2[isite]
        if self.nbs_source == 'nb_sets':
            all_nbs_voro_indices = set()
            for cn2, nb_sets in structure_environments.neighbors_sets[isite].items():
                for inb_set2, nb_set2 in enumerate(nb_sets):
                    if cn == cn2 and inb_set == inb_set:
                        continue
                    all_nbs_voro_indices.update(nb_set2.site_voronoi_indices)
        elif self.nbs_source == "voronoi":
            all_nbs_voro_indices = set(range(len(voronoi)))
        else:
            raise ValueError('"nbs_source" should be one of ["nb_sets", "voronoi"]')
        all_nbs_indices_except_nb_set = all_nbs_voro_indices.difference(nb_set.site_voronoi_indices)
        normalized_distances = [voronoi[inb]['normalized_distance'] for inb in all_nbs_indices_except_nb_set]
        if len(normalized_distances) == 0:
            return 1.0
        if len(nb_set) == 0:
            return 0.0
        nb_set_max_normalized_distance = max(nb_set.normalized_distances)
        return self.weight_rf.eval(min(normalized_distances)-nb_set_max_normalized_distance)

    def __eq__(self, other):
        return self.__class__ == other.__class__

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "weight_function": self.weight_function,
                "nbs_source": self.nbs_source
                }

    @classmethod
    def from_dict(cls, dd):
        return cls(weight_function=dd['weight_function'], nbs_source=dd['nbs_source'])


class WeightedNbSetChemenvStrategy(AbstractChemenvStrategy):
    """
    WeightedNbSetChemenvStrategy
    """
    STRATEGY_DESCRIPTION = '    WeightedNbSetChemenvStrategy'
    DEFAULT_CE_ESTIMATOR = {'function': 'power2_inverse_power2_decreasing',
                            'options': {'max_csm': 8.0}}

    def __init__(self, structure_environments=None,
                 additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB,
                 symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE,
                 nb_set_weights=None,
                 ce_estimator=DEFAULT_CE_ESTIMATOR):
        """
        Constructor for the WeightedNbSetChemenvStrategy.
        :param structure_environments: StructureEnvironments object containing all the information on the
        coordination of the sites in a structure
        """
        AbstractChemenvStrategy.__init__(self, structure_environments, symmetry_measure_type=symmetry_measure_type)
        self._additional_condition = additional_condition
        if nb_set_weights is None:
            raise ValueError()
        self.nb_set_weights = nb_set_weights
        self.ordered_weights = []
        for nb_set_weight in self.nb_set_weights:
            self.ordered_weights.append({'weight': nb_set_weight, 'name': nb_set_weight.SHORT_NAME})
        self.ce_estimator = ce_estimator
        self.ce_estimator_ratio_function = CSMInfiniteRatioFunction.from_dict(self.ce_estimator)
        self.ce_estimator_fractions = self.ce_estimator_ratio_function.fractions

    @property
    def uniquely_determines_coordination_environments(self):
        return False

    def get_site_coordination_environments_fractions(self, site, isite=None, dequivsite=None, dthissite=None,
                                                     mysym=None, ordered=True, min_fraction=0.0, return_maps=True,
                                                     return_strategy_dict_info=False, return_all=False):
        if isite is None or dequivsite is None or dthissite is None or mysym is None:
            [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        site_nb_sets = self.structure_environments.neighbors_sets[isite]
        if site_nb_sets is None:
            return None
        cn_maps = []
        for cn, nb_sets in site_nb_sets.items():
            for inb_set, nb_set in enumerate(nb_sets):
                #CHECK THE ADDITIONAL CONDITION HERE ?
                cn_maps.append((cn, inb_set))
        weights_additional_info = {'weights': {isite: {}}}
        for wdict in self.ordered_weights:
            cn_maps_new = []
            weight = wdict['weight']
            weight_name = wdict['name']
            for cn_map in cn_maps:
                nb_set = site_nb_sets[cn_map[0]][cn_map[1]]
                w_nb_set = weight.weight(nb_set=nb_set, structure_environments=self.structure_environments,
                                         cn_map=cn_map, additional_info=weights_additional_info)
                if cn_map not in weights_additional_info['weights'][isite]:
                    weights_additional_info['weights'][isite][cn_map] = {}
                weights_additional_info['weights'][isite][cn_map][weight_name] = w_nb_set
                if w_nb_set > 0.0:
                    cn_maps_new.append(cn_map)
            cn_maps = cn_maps_new
        for cn_map, weights in weights_additional_info['weights'][isite].items():
            weights_additional_info['weights'][isite][cn_map]['Product'] = np.product(list(weights.values()))

        w_nb_sets = {cn_map: weights['Product']
                     for cn_map, weights in weights_additional_info['weights'][isite].items()}
        w_nb_sets_total = np.sum(list(w_nb_sets.values()))
        nb_sets_fractions = {cn_map: w_nb_set / w_nb_sets_total for cn_map, w_nb_set in w_nb_sets.items()}
        for cn_map in weights_additional_info['weights'][isite]:
            weights_additional_info['weights'][isite][cn_map]['NbSetFraction'] = nb_sets_fractions[cn_map]
        ce_symbols = []
        ce_dicts = []
        ce_fractions = []
        ce_dict_fractions = []
        ce_maps = []
        site_ce_list = self.structure_environments.ce_list[isite]
        if return_all:
            for cn_map, nb_set_fraction in nb_sets_fractions.items():
                cn = cn_map[0]
                inb_set = cn_map[1]
                site_ce_nb_set = site_ce_list[cn][inb_set]
                if site_ce_nb_set is None:
                    continue
                mingeoms = site_ce_nb_set.minimum_geometries(symmetry_measure_type=self.symmetry_measure_type)
                if len(mingeoms) > 0:
                    csms = [ce_dict['other_symmetry_measures'][self.symmetry_measure_type]
                            for ce_symbol, ce_dict in mingeoms]
                    fractions = self.ce_estimator_fractions(csms)
                    if fractions is None:
                        ce_symbols.append('UNCLEAR:{:d}'.format(cn))
                        ce_dicts.append(None)
                        ce_fractions.append(nb_set_fraction)
                        all_weights = weights_additional_info['weights'][isite][cn_map]
                        dict_fractions = {wname: wvalue for wname, wvalue in all_weights.items()}
                        dict_fractions['CEFraction'] = None
                        dict_fractions['Fraction'] = nb_set_fraction
                        ce_dict_fractions.append(dict_fractions)
                        ce_maps.append(cn_map)
                    else:
                        for ifraction, fraction in enumerate(fractions):
                            ce_symbols.append(mingeoms[ifraction][0])
                            ce_dicts.append(mingeoms[ifraction][1])
                            ce_fractions.append(nb_set_fraction * fraction)
                            all_weights = weights_additional_info['weights'][isite][cn_map]
                            dict_fractions = {wname: wvalue for wname, wvalue in all_weights.items()}
                            dict_fractions['CEFraction'] = fraction
                            dict_fractions['Fraction'] = nb_set_fraction * fraction
                            ce_dict_fractions.append(dict_fractions)
                            ce_maps.append(cn_map)
                else:
                    ce_symbols.append('UNCLEAR:{:d}'.format(cn))
                    ce_dicts.append(None)
                    ce_fractions.append(nb_set_fraction)
                    all_weights = weights_additional_info['weights'][isite][cn_map]
                    dict_fractions = {wname: wvalue for wname, wvalue in all_weights.items()}
                    dict_fractions['CEFraction'] = None
                    dict_fractions['Fraction'] = nb_set_fraction
                    ce_dict_fractions.append(dict_fractions)
                    ce_maps.append(cn_map)
        else:
            for cn_map, nb_set_fraction in nb_sets_fractions.items():
                if nb_set_fraction > 0.0:
                    cn = cn_map[0]
                    inb_set = cn_map[1]
                    site_ce_nb_set = site_ce_list[cn][inb_set]
                    mingeoms = site_ce_nb_set.minimum_geometries(symmetry_measure_type=self._symmetry_measure_type)
                    csms = [ce_dict['other_symmetry_measures'][self._symmetry_measure_type]
                            for ce_symbol, ce_dict in mingeoms]
                    fractions = self.ce_estimator_fractions(csms)
                    for ifraction, fraction in enumerate(fractions):
                        if fraction > 0.0:
                            ce_symbols.append(mingeoms[ifraction][0])
                            ce_dicts.append(mingeoms[ifraction][1])
                            ce_fractions.append(nb_set_fraction * fraction)
                            all_weights = weights_additional_info['weights'][isite][cn_map]
                            dict_fractions = {wname: wvalue for wname, wvalue in all_weights.items()}
                            dict_fractions['CEFraction'] = fraction
                            dict_fractions['Fraction'] = nb_set_fraction * fraction
                            ce_dict_fractions.append(dict_fractions)
                            ce_maps.append(cn_map)
        if ordered:
            indices = np.argsort(ce_fractions)[::-1]
        else:
            indices = list(range(len(ce_fractions)))

        fractions_info_list = [
            {'ce_symbol': ce_symbols[ii], 'ce_dict': ce_dicts[ii], 'ce_fraction': ce_fractions[ii]}
            for ii in indices if ce_fractions[ii] >= min_fraction]

        if return_maps:
            for ifinfo, ii in enumerate(indices):
                if ce_fractions[ii] >= min_fraction:
                    fractions_info_list[ifinfo]['ce_map'] = ce_maps[ii]
        if return_strategy_dict_info:
            for ifinfo, ii in enumerate(indices):
                if ce_fractions[ii] >= min_fraction:
                    fractions_info_list[ifinfo]['strategy_info'] = ce_dict_fractions[ii]
        return fractions_info_list

    def get_site_coordination_environment(self, site):
        pass

    def get_site_neighbors(self, site):
        pass

    def get_site_coordination_environments(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None,
                                           return_maps=False):
        if isite is None or dequivsite is None or dthissite is None or mysym is None:
            [isite, dequivsite, dthissite, mysym] = self.equivalent_site_index_and_transform(site)
        return [self.get_site_coordination_environment(site=site, isite=isite, dequivsite=dequivsite,
                                                       dthissite=dthissite, mysym=mysym, return_map=return_maps)]

    def __eq__(self, other):
        return (self.__class__.__name__ == other.__class__.__name__ and
                self._additional_condition == other._additional_condition and
                self.symmetry_measure_type == other.symmetry_measure_type and
                self.nb_set_weights == other.nb_set_weights and
                self.ce_estimator == other.ce_estimator)

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        """
        Bson-serializable dict representation of the WeightedNbSetChemenvStrategy object.
        :return: Bson-serializable dict representation of the WeightedNbSetChemenvStrategy object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "additional_condition": self._additional_condition,
                "symmetry_measure_type": self.symmetry_measure_type,
                "nb_set_weights": [nb_set_weight.as_dict() for nb_set_weight in self.nb_set_weights],
                "ce_estimator": self.ce_estimator,
                }

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the WeightedNbSetChemenvStrategy object from a dict representation of the
        WeightedNbSetChemenvStrategy object created using the as_dict method.
        :param d: dict representation of the WeightedNbSetChemenvStrategy object
        :return: WeightedNbSetChemenvStrategy object
        """
        return cls(additional_condition=d["additional_condition"],
                   symmetry_measure_type=d["symmetry_measure_type"],
                   nb_set_weights=d["nb_set_weights"],
                   ce_estimator=d["ce_estimator"])


class MultiWeightsChemenvStrategy(WeightedNbSetChemenvStrategy):
    """
    MultiWeightsChemenvStrategy
    """
    STRATEGY_DESCRIPTION = '    Multi Weights ChemenvStrategy'
    # STRATEGY_INFO_FIELDS = ['cn_map_surface_fraction', 'cn_map_surface_weight',
    #                         'cn_map_mean_csm', 'cn_map_csm_weight',
    #                         'cn_map_delta_csm', 'cn_map_delta_csms_cn_map2', 'cn_map_delta_csm_weight',
    #                         'cn_map_cn_weight',
    #                         'cn_map_fraction', 'cn_map_ce_fraction', 'ce_fraction']
    DEFAULT_CE_ESTIMATOR = {'function': 'power2_inverse_power2_decreasing',
                            'options': {'max_csm': 8.0}}
    DEFAULT_DIST_ANG_AREA_WEIGHT = {}

    def __init__(self, structure_environments=None,
                 additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB,
                 symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE,
                 dist_ang_area_weight=None,
                 self_csm_weight=None,
                 delta_csm_weight=None,
                 cn_bias_weight=None,
                 angle_weight=None,
                 normalized_angle_distance_weight=None,
                 ce_estimator=DEFAULT_CE_ESTIMATOR
                 ):
        """
        Constructor for the MultiWeightsChemenvStrategy.
        :param structure_environments: StructureEnvironments object containing all the information on the
        coordination of the sites in a structure
        """
        self._additional_condition = additional_condition
        self.dist_ang_area_weight = dist_ang_area_weight
        self.angle_weight = angle_weight
        self.normalized_angle_distance_weight = normalized_angle_distance_weight
        self.self_csm_weight = self_csm_weight
        self.delta_csm_weight = delta_csm_weight
        self.cn_bias_weight = cn_bias_weight
        self.ordered_weights = []
        nb_sets_weights = []
        if dist_ang_area_weight is not None:
            self.ordered_weights.append({'weight': dist_ang_area_weight, 'name': 'DistAngArea'})
            nb_sets_weights.append(dist_ang_area_weight)
        if self_csm_weight is not None:
            self.ordered_weights.append({'weight': self_csm_weight, 'name': 'SelfCSM'})
            nb_sets_weights.append(self_csm_weight)
        if delta_csm_weight is not None:
            self.ordered_weights.append({'weight': delta_csm_weight, 'name': 'DeltaCSM'})
            nb_sets_weights.append(delta_csm_weight)
        if cn_bias_weight is not None:
            self.ordered_weights.append({'weight': cn_bias_weight, 'name': 'CNBias'})
            nb_sets_weights.append(cn_bias_weight)
        if angle_weight is not None:
            self.ordered_weights.append({'weight': angle_weight, 'name': 'Angle'})
            nb_sets_weights.append(angle_weight)
        if normalized_angle_distance_weight is not None:
            self.ordered_weights.append({'weight': normalized_angle_distance_weight, 'name': 'NormalizedAngDist'})
            nb_sets_weights.append(normalized_angle_distance_weight)

        self.ce_estimator = ce_estimator
        self.ce_estimator_ratio_function = CSMInfiniteRatioFunction.from_dict(self.ce_estimator)
        self.ce_estimator_fractions = self.ce_estimator_ratio_function.fractions
        WeightedNbSetChemenvStrategy.__init__(self, structure_environments,
                                              additional_condition=additional_condition,
                                              symmetry_measure_type=symmetry_measure_type,
                                              nb_set_weights=nb_sets_weights,
                                              ce_estimator=ce_estimator)

    @classmethod
    def stats_article_weights_parameters(cls):
        self_csm_weight = SelfCSMNbSetWeight(weight_estimator={'function': 'power2_decreasing_exp',
                                                               'options': {'max_csm': 8.0,
                                                                           'alpha': 1.0}})
        surface_definition = {'type': 'standard_elliptic',
                              'distance_bounds': {'lower': 1.15, 'upper': 2.0},
                              'angle_bounds': {'lower': 0.05, 'upper': 0.75}}
        da_area_weight = DistanceAngleAreaNbSetWeight(weight_type='has_intersection',
                                                      surface_definition=surface_definition,
                                                      nb_sets_from_hints='fallback_to_source',
                                                      other_nb_sets='0_weight',
                                                      additional_condition=DistanceAngleAreaNbSetWeight.AC.ONLY_ACB)
        symmetry_measure_type = 'csm_wcs_ctwcc'
        delta_weight = DeltaCSMNbSetWeight.delta_cn_specifics()
        bias_weight = None
        angle_weight = None
        nad_weight = None
        return cls(dist_ang_area_weight=da_area_weight,
                   self_csm_weight=self_csm_weight,
                   delta_csm_weight=delta_weight,
                   cn_bias_weight=bias_weight,
                   angle_weight=angle_weight,
                   normalized_angle_distance_weight=nad_weight,
                   symmetry_measure_type=symmetry_measure_type)

    @property
    def uniquely_determines_coordination_environments(self):
        return False

    def __eq__(self, other):
        return (self.__class__.__name__ == other.__class__.__name__ and
                self._additional_condition == other._additional_condition and
                self.symmetry_measure_type == other.symmetry_measure_type and
                self.dist_ang_area_weight == other.dist_ang_area_weight and
                self.self_csm_weight == other.self_csm_weight and
                self.delta_csm_weight == other.delta_csm_weight and
                self.cn_bias_weight == other.cn_bias_weight and
                self.angle_weight == other.angle_weight and
                self.normalized_angle_distance_weight == other.normalized_angle_distance_weight and
                self.ce_estimator == other.ce_estimator)

    def __ne__(self, other):
        return not self == other

    def as_dict(self):
        """
        Bson-serializable dict representation of the MultiWeightsChemenvStrategy object.
        :return: Bson-serializable dict representation of the MultiWeightsChemenvStrategy object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "additional_condition": self._additional_condition,
                "symmetry_measure_type": self.symmetry_measure_type,
                "dist_ang_area_weight": self.dist_ang_area_weight.as_dict()
                if self.dist_ang_area_weight is not None else None,
                "self_csm_weight": self.self_csm_weight.as_dict()
                if self.self_csm_weight is not None else None,
                "delta_csm_weight": self.delta_csm_weight.as_dict()
                if self.delta_csm_weight is not None else None,
                "cn_bias_weight": self.cn_bias_weight.as_dict()
                if self.cn_bias_weight is not None else None,
                "angle_weight": self.angle_weight.as_dict()
                if self.angle_weight is not None else None,
                "normalized_angle_distance_weight": self.normalized_angle_distance_weight.as_dict()
                if self.normalized_angle_distance_weight is not None else None,
                "ce_estimator": self.ce_estimator,
                }

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the MultiWeightsChemenvStrategy object from a dict representation of the
        MultipleAbundanceChemenvStrategy object created using the as_dict method.
        :param d: dict representation of the MultiWeightsChemenvStrategy object
        :return: MultiWeightsChemenvStrategy object
        """
        if d["normalized_angle_distance_weight"] is not None:
            nad_w = NormalizedAngleDistanceNbSetWeight.from_dict(d["normalized_angle_distance_weight"])
        else:
            nad_w = None
        return cls(additional_condition=d["additional_condition"],
                   symmetry_measure_type=d["symmetry_measure_type"],
                   dist_ang_area_weight=DistanceAngleAreaNbSetWeight.from_dict(d["dist_ang_area_weight"])
                   if d["dist_ang_area_weight"] is not None else None,
                   self_csm_weight=SelfCSMNbSetWeight.from_dict(d["self_csm_weight"])
                   if d["self_csm_weight"] is not None else None,
                   delta_csm_weight=DeltaCSMNbSetWeight.from_dict(d["delta_csm_weight"])
                   if d["delta_csm_weight"] is not None else None,
                   cn_bias_weight=CNBiasNbSetWeight.from_dict(d["cn_bias_weight"])
                   if d["cn_bias_weight"] is not None else None,
                   angle_weight=AngleNbSetWeight.from_dict(d["angle_weight"])
                   if d["angle_weight"] is not None else None,
                   normalized_angle_distance_weight=nad_w,
                   ce_estimator=d["ce_estimator"])
