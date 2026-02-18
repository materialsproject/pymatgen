"""
This module provides so-called "strategies" to determine the coordination environments of an atom in a structure.
Some strategies can favour larger or smaller environments. Some strategies uniquely identifies the environments while
some others can identify the environment as a "mix" of several environments, each of which is assigned with a given
fraction. The choice of the strategy depends on the purpose of the user.
"""

from __future__ import annotations

import abc
from typing import TYPE_CHECKING

import numpy as np
from monty.json import MSONable
from scipy.stats import gmean

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.coordination_environments.voronoi import DetailedVoronoiContainer
from pymatgen.analysis.chemenv.utils.chemenv_errors import EquivalentSiteSearchError
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import get_lower_and_upper_f
from pymatgen.analysis.chemenv.utils.defs_utils import AdditionalConditions
from pymatgen.analysis.chemenv.utils.func_utils import (
    CSMFiniteRatioFunction,
    CSMInfiniteRatioFunction,
    DeltaCSMRatioFunction,
    RatioFunction,
)
from pymatgen.core.operations import SymmOp
from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from typing import ClassVar, Self

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

MPSYMBOL_TO_CN = AllCoordinationGeometries().get_symbol_cn_mapping()
ALLCG = AllCoordinationGeometries()


class StrategyOption(MSONable, abc.ABC):
    """Abstract class for the options of the chemenv strategies."""

    allowed_values: str | None = None

    @abc.abstractmethod
    def as_dict(self):
        """A JSON-serializable dict representation of this strategy option."""


class DistanceCutoffFloat(float, StrategyOption):
    """Distance cutoff in a strategy."""

    allowed_values = "Real number between 1 and +infinity"

    def __new__(cls, cutoff) -> Self:
        """Special float that should be between 1 and infinity.

        Args:
            cutoff: Distance cutoff.
        """
        flt = float.__new__(cls, cutoff)
        if flt < 1:
            raise ValueError("Distance cutoff should be between 1 and +infinity")
        return flt

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "value": self,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize distance cutoff from dict.

        Args:
            dct (dict): Dict representation of the distance cutoff.
        """
        return cls(dct["value"])


class AngleCutoffFloat(float, StrategyOption):
    """Angle cutoff in a strategy."""

    allowed_values = "Real number between 0 and 1"

    def __new__(cls, cutoff) -> Self:
        """Special float that should be between 0 and 1.

        Args:
            cutoff: Angle cutoff.
        """
        flt = float.__new__(cls, cutoff)
        if not 0 <= flt <= 1:
            raise ValueError(f"Angle cutoff should be between 0 and 1, got {flt}")
        return flt

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "value": self,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize angle cutoff from dict.

        Args:
            dct (dict): Dict representation of the angle cutoff.
        """
        return cls(dct["value"])


class CSMFloat(float, StrategyOption):
    """Real number representing a Continuous Symmetry Measure."""

    allowed_values = "Real number between 0 and 100"

    def __new__(cls, cutoff) -> Self:
        """Special float that should be between 0 and 100.

        Args:
            cutoff: CSM.
        """
        flt = float.__new__(cls, cutoff)
        if not 0 <= flt <= 100:
            raise ValueError(f"Continuous symmetry measure limits should be between 0 and 100, got {flt}")
        return flt

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "value": self,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize CSM from dict.

        Args:
           dct (dict): Dict representation of the CSM.
        """
        return cls(dct["value"])


class AdditionalConditionInt(int, StrategyOption):
    """Integer representing an additional condition in a strategy."""

    allowed_values = "Integer amongst :\n"
    for integer, description in AdditionalConditions.CONDITION_DESCRIPTION.items():
        allowed_values += f" - {integer} for {description!r}\n"

    def __new__(cls, integer) -> Self:
        """Special int representing additional conditions."""
        if str(int(integer)) != str(integer):
            raise ValueError(f"Additional condition {integer} is not an integer")
        integer = int.__new__(cls, integer)
        if integer not in AdditionalConditions.ALL:
            raise ValueError(f"Additional condition {integer} is not allowed")
        return integer

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "value": self,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize additional condition from dict.

        Args:
           dct (dict): Dict representation of the additional condition.
        """
        return cls(dct["value"])


class AbstractChemenvStrategy(MSONable, abc.ABC):
    """Base class to define a Chemenv strategy for the neighbors and coordination environment
    to be applied to a StructureEnvironments object.
    """

    AC = AdditionalConditions()
    STRATEGY_OPTIONS: ClassVar[dict[str, dict]] = {}
    STRATEGY_DESCRIPTION: str | None = None
    STRATEGY_INFO_FIELDS: ClassVar[list] = []
    DEFAULT_SYMMETRY_MEASURE_TYPE = "csm_wcs_ctwcc"

    def __init__(
        self,
        structure_environments=None,
        symmetry_measure_type=DEFAULT_SYMMETRY_MEASURE_TYPE,
    ):
        """
        Abstract constructor for the all chemenv strategies.

        Args:
            structure_environments: StructureEnvironments object containing all the information on the
                coordination of the sites in a structure.
        """
        self.structure_environments = None
        if structure_environments is not None:
            self.set_structure_environments(structure_environments)
        self._symmetry_measure_type = symmetry_measure_type

    @property
    def symmetry_measure_type(self):
        """Type of symmetry measure."""
        return self._symmetry_measure_type

    def set_structure_environments(self, structure_environments):
        """Set the structure environments to this strategy.

        Args:
            structure_environments: StructureEnvironments object.
        """
        self.structure_environments = structure_environments
        if not isinstance(self.structure_environments.voronoi, DetailedVoronoiContainer):
            raise TypeError('Voronoi Container not of type "DetailedVoronoiContainer"')
        self.prepare_symmetries()

    # AL:  Creating SpacegroupAnalyzer is expensive. Cache by a cheap structure fingerprint and skip when unchanged.
    def prepare_symmetries(self):
        """Prepare the symmetries for the structure contained in the structure environments."""
        try:
            s = self.structure_environments.structure
            # cheap, deterministic fingerprint (no heavy symmetry ops)
            fp = (
                tuple(s.lattice.matrix.ravel().round(10)),
                tuple(sp.symbol for sp in s.species),
                tuple(map(tuple, np.asarray(s.frac_coords).round(10))),
            )
            if getattr(self, "_symops_cache_fp", None) == fp and getattr(self, "symops", None) is not None:
                return  # cache hit

            self.spg_analyzer = SpacegroupAnalyzer(s)
            # NOTE: cartesian=False (default) keeps operations in fractional space and is cheaper to compare later
            self.symops = self.spg_analyzer.get_symmetry_operations()
            self._symops_cache_fp = fp
        except Exception:
            self.symops = []
            self._symops_cache_fp = None

    # AL The inner loop constructs many PeriodicSite objects and calls is_periodic_image repeatedly.
    # Comparing fractional coordinates modulo 1 is enough and much cheaper.
    def equivalent_site_index_and_transform(self, psite):
        """Get the equivalent site and corresponding symmetry+translation transformations.

        Args:
            psite: Periodic site.

        Returns:
            Equivalent site in the unit cell, translations and symmetry transformation.
        """
        # Get the index of the site in the unit cell of which the PeriodicSite psite is a replica.
        site_idx = 0
        try:
            site_idx = self.structure_environments.structure.index(psite)
        except ValueError:
            try:
                uc_psite = psite.to_unit_cell()
                site_idx = self.structure_environments.structure.index(uc_psite)
            except ValueError:
                for site_idx2, site2 in enumerate(self.structure_environments.structure):
                    if psite.is_periodic_image(site2):
                        site_idx = site_idx2
                        break
        # Get the translation between psite and its corresponding site in the unit cell (Translation I)
        this_site = self.structure_environments.structure[site_idx]
        dist_this_site = psite.frac_coords - this_site.frac_coords
        # Get the translation between the equivalent site for which the neighbors have been computed and the site in
        # the unit cell that corresponds to psite (Translation II)
        equiv_site = self.structure_environments.structure[
            self.structure_environments.sites_map[site_idx]
        ].to_unit_cell()
        # equivsite = self.structure_environments.structure[self.structure_environments.sites_map[isite]]
        dist_equiv_site = (
            self.structure_environments.structure[self.structure_environments.sites_map[site_idx]].frac_coords
            - equiv_site.frac_coords
        )
        found = False
        # Find the symmetry that applies the site in the unit cell to the equivalent site, as well as the translation
        # that gets back the site to the unit cell (Translation III)
        # TODO: check that these tolerances are needed, now that the structures are refined before analyzing envs
        tolerances = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
        d_this_site2 = (0, 0, 0)
        sym_trafo = None
        eq_frac = np.asarray(equiv_site.frac_coords, dtype=float)
        this_frac = np.asarray(this_site.frac_coords, dtype=float)

        def is_same_frac(a, b, tol):
            # compare on the 3-torus: bring difference to [-0.5, 0.5) and test
            d = a - b
            d -= np.round(d)
            return np.all(np.abs(d) <= tol)

        for tol in tolerances:
            # identity fast-path
            if is_same_frac(eq_frac, this_frac, tol):
                # use identity SymmOp (downstream assumes mysym.operate exists)
                sym_trafo = SymmOp.from_rotation_and_translation(np.eye(3), [0.0, 0.0, 0.0])
                d = this_frac - eq_frac
                d -= np.round(d)
                d_this_site2 = tuple(d)
                found = True
                break

            for sym_op in self.symops:
                new_frac = sym_op.operate(eq_frac)
                if is_same_frac(new_frac, this_frac, tol):
                    sym_trafo = sym_op
                    d = this_frac - new_frac
                    d -= np.round(d)
                    d_this_site2 = tuple(d)
                    found = True
                    break
            if found:
                break

        if not found:
            raise EquivalentSiteSearchError(psite)

        equivalent_site_map = self.structure_environments.sites_map[site_idx]
        return (
            equivalent_site_map,
            dist_equiv_site,
            dist_this_site + d_this_site2,
            sym_trafo,
        )

    @abc.abstractmethod
    def get_site_neighbors(self, site):
        """
        Applies the strategy to the structure_environments object in order to get the neighbors of a given site.

        Args:
            site: Site for which the neighbors are looked for
            structure_environments: StructureEnvironments object containing all the information needed to get the
                neighbors of the site

        Returns:
            The list of neighbors of the site. For complex strategies, where one allows multiple solutions, this
            can return a list of list of neighbors.
        """
        raise NotImplementedError

    @property
    def uniquely_determines_coordination_environments(self):
        """True if the strategy leads to a unique coordination environment."""
        raise NotImplementedError

    @abc.abstractmethod
    def get_site_coordination_environment(self, site):
        """
        Applies the strategy to the structure_environments object in order to define the coordination environment of
        a given site.

        Args:
            site: Site for which the coordination environment is looked for

        Returns:
            The coordination environment of the site. For complex strategies, where one allows multiple
            solutions, this can return a list of coordination environments for the site.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def get_site_coordination_environments(self, site):
        """
        Applies the strategy to the structure_environments object in order to define the coordination environment of
        a given site.

        Args:
            site: Site for which the coordination environment is looked for

        Returns:
            The coordination environment of the site. For complex strategies, where one allows multiple
            solutions, this can return a list of coordination environments for the site.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def get_site_coordination_environments_fractions(
        self,
        site,
        isite=None,
        dequivsite=None,
        dthissite=None,
        mysym=None,
        ordered=True,
        min_fraction=0,
        return_maps=True,
        return_strategy_dict_info=False,
    ):
        """
        Applies the strategy to the structure_environments object in order to define the coordination environment of
        a given site.

        Args:
            site: Site for which the coordination environment is looked for

        Returns:
            The coordination environment of the site. For complex strategies, where one allows multiple
            solutions, this can return a list of coordination environments for the site.
        """
        raise NotImplementedError

    # AL Tight path during many-site scans. Cache locals and build neighbor list without intermediate dict churn.
    def get_site_ce_fractions_and_neighbors(self, site, full_ce_info=False, strategy_info=False):
        site_idx, dist_equiv_site, dist_this_site, mysym = self.equivalent_site_index_and_transform(site)
        geoms_and_maps_list = self.get_site_coordination_environments_fractions(
            site=site,
            isite=site_idx,
            dequivsite=dist_equiv_site,
            dthissite=dist_this_site,
            mysym=mysym,
            return_maps=True,
            return_strategy_dict_info=True,
        )
        if geoms_and_maps_list is None:
            return None
        site_nbs_sets = self.structure_environments.neighbors_sets[site_idx]
        ce_and_neighbors = []
        for fractions_dict in geoms_and_maps_list:
            cn0, inb0 = fractions_dict["ce_map"]
            nbs = site_nbs_sets[cn0][inb0].neighb_sites_and_indices
            # reuse objects, minimal dicts
            fractions_dict["neighbors"] = [{"site": x["site"], "index": x["index"]} for x in nbs]
            ce_and_neighbors.append(fractions_dict)
        return ce_and_neighbors

    def set_option(self, option_name, option_value):
        """Set up a given option for this strategy.

        Args:
            option_name: Name of the option.
            option_value: Value for this option.
        """
        setattr(self, option_name, option_value)

    def setup_options(self, all_options_dict):
        """Set up options for this strategy based on a dict.

        Args:
            all_options_dict: Dict of option_name->option_value.
        """
        for option_name, option_value in all_options_dict.items():
            self.set_option(option_name, option_value)

    @abc.abstractmethod
    def __eq__(self, other: object) -> bool:
        """
        Equality method that should be implemented for any strategy.

        Args:
            other: strategy to be compared with the current one
        """
        raise NotImplementedError

    def __str__(self):
        out = f"  Chemenv Strategy {type(self).__name__!r}\n"
        out += f"  {'=' * (19 + len(type(self).__name__))}\n\n"
        out += f"  Description :\n  {'-' * 13}\n"
        out += self.STRATEGY_DESCRIPTION
        out += "\n\n"
        out += f"  Options :\n  {'-' * 9}\n"
        for option_name in self.STRATEGY_OPTIONS:
            out += f"   - {option_name} : {getattr(self, option_name)}\n"
        return out

    @abc.abstractmethod
    def as_dict(self):
        """
        Bson-serializable dict representation of the SimplestChemenvStrategy object.

        Returns:
            Bson-serializable dict representation of the SimplestChemenvStrategy object.
        """
        raise NotImplementedError

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Reconstructs the SimpleAbundanceChemenvStrategy object from a dict representation of the
        SimpleAbundanceChemenvStrategy object created using the as_dict method.

        Args:
            dct: dict representation of the SimpleAbundanceChemenvStrategy object

        Returns:
            StructureEnvironments object.
        """
        raise NotImplementedError


class SimplestChemenvStrategy(AbstractChemenvStrategy):
    """
    Simplest ChemenvStrategy using fixed angle and distance parameters for the definition of neighbors in the
    Voronoi approach. The coordination environment is then given as the one with the lowest continuous symmetry measure.
    """

    # Default values for the distance and angle cutoffs
    DEFAULT_DISTANCE_CUTOFF = 1.4
    DEFAULT_ANGLE_CUTOFF = 0.3
    DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF = 10
    DEFAULT_ADDITIONAL_CONDITION = AbstractChemenvStrategy.AC.ONLY_ACB
    STRATEGY_OPTIONS: ClassVar[dict[str, dict]] = {
        "distance_cutoff": {
            "type": DistanceCutoffFloat,
            "internal": "_distance_cutoff",
            "default": DEFAULT_DISTANCE_CUTOFF,
        },
        "angle_cutoff": {
            "type": AngleCutoffFloat,
            "internal": "_angle_cutoff",
            "default": DEFAULT_ANGLE_CUTOFF,
        },
        "additional_condition": {
            "type": AdditionalConditionInt,
            "internal": "_additional_condition",
            "default": DEFAULT_ADDITIONAL_CONDITION,
        },
        "continuous_symmetry_measure_cutoff": {
            "type": CSMFloat,
            "internal": "_continuous_symmetry_measure_cutoff",
            "default": DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF,
        },
    }

    STRATEGY_DESCRIPTION = (
        "Simplest ChemenvStrategy using fixed angle and distance parameters \n"
        "for the definition of neighbors in the Voronoi approach. \n"
        "The coordination environment is then given as the one with the \n"
        "lowest continuous symmetry measure."
    )

    def __init__(
        self,
        structure_environments=None,
        distance_cutoff=DEFAULT_DISTANCE_CUTOFF,
        angle_cutoff=DEFAULT_ANGLE_CUTOFF,
        additional_condition=DEFAULT_ADDITIONAL_CONDITION,
        continuous_symmetry_measure_cutoff=DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF,
        symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE,
    ):
        """
        Constructor for this SimplestChemenvStrategy.

        Args:
            distance_cutoff: Distance cutoff used
            angle_cutoff: Angle cutoff used.
        """
        AbstractChemenvStrategy.__init__(self, structure_environments, symmetry_measure_type=symmetry_measure_type)
        self.distance_cutoff = distance_cutoff
        self.angle_cutoff = angle_cutoff
        self.additional_condition = additional_condition
        self.continuous_symmetry_measure_cutoff = continuous_symmetry_measure_cutoff

    @property
    def uniquely_determines_coordination_environments(self):
        """Whether this strategy uniquely determines coordination environments."""
        return True

    @property
    def distance_cutoff(self):
        """Distance cutoff used."""
        return self._distance_cutoff

    @distance_cutoff.setter
    def distance_cutoff(self, distance_cutoff):
        """Set the distance cutoff for this strategy.

        Args:
            distance_cutoff: Distance cutoff.
        """
        self._distance_cutoff = DistanceCutoffFloat(distance_cutoff)

    @property
    def angle_cutoff(self):
        """Angle cutoff used."""
        return self._angle_cutoff

    @angle_cutoff.setter
    def angle_cutoff(self, angle_cutoff):
        """Set the angle cutoff for this strategy.

        Args:
            angle_cutoff: Angle cutoff.
        """
        self._angle_cutoff = AngleCutoffFloat(angle_cutoff)

    @property
    def additional_condition(self) -> AdditionalConditionInt:
        """Additional condition for this strategy."""
        return self._additional_condition

    @additional_condition.setter
    def additional_condition(self, additional_condition):
        """Set the additional condition for this strategy.

        Args:
            additional_condition: Additional condition.
        """
        self._additional_condition = AdditionalConditionInt(additional_condition)

    @property
    def continuous_symmetry_measure_cutoff(self):
        """CSM cutoff used."""
        return self._continuous_symmetry_measure_cutoff

    @continuous_symmetry_measure_cutoff.setter
    def continuous_symmetry_measure_cutoff(self, continuous_symmetry_measure_cutoff):
        """Set the CSM cutoff for this strategy.

        Args:
            continuous_symmetry_measure_cutoff: CSM cutoff
        """
        self._continuous_symmetry_measure_cutoff = CSMFloat(continuous_symmetry_measure_cutoff)

    # AL: prebind + listcomp
    def get_site_neighbors(self, site, isite=None, dequivsite=None, dthissite=None, mysym=None):
        if isite is None:
            isite, dequivsite, dthissite, mysym = self.equivalent_site_index_and_transform(site)

        _ce, cn_map = self.get_site_coordination_environment(
            site=site,
            isite=isite,
            dequivsite=dequivsite,
            dthissite=dthissite,
            mysym=mysym,
            return_map=True,
        )

        nb_set = self.structure_environments.neighbors_sets[isite][cn_map[0]][cn_map[1]]
        eq_site_ps = nb_set.neighb_sites

        op = mysym.operate
        de = dequivsite
        dt = dthissite
        return [PeriodicSite(ps._species, op(ps.frac_coords + de) + dt, ps._lattice) for ps in eq_site_ps]

    # build (idp,iap)->(cn,inb) map once per site/AC
    def _get_cn_map_cache(self, isite):
        cache = getattr(self, "_cn_map_cache", None)
        if cache is None:
            cache = self._cn_map_cache = {}
        key = (isite, int(self._additional_condition))
        mapping = cache.get(key)
        if mapping is not None:
            return mapping

        mapping = {}
        site_nb_sets = self.structure_environments.neighbors_sets[isite]
        for cn, nb_sets in site_nb_sets.items():
            for inb_set, nb_set in enumerate(nb_sets):
                for src in nb_set.sources:
                    if src.get("origin") == "dist_ang_ac_voronoi" and src.get("ac") == self._additional_condition:
                        mapping[(src["idp"], src["iap"])] = (cn, inb_set)
        cache[key] = mapping
        return mapping

    # faster index find + cached map
    def get_site_coordination_environment(
        self,
        site,
        isite=None,
        dequivsite=None,
        dthissite=None,
        mysym=None,
        return_map=False,
    ):
        if isite is None:
            isite, *_ = self.equivalent_site_index_and_transform(site)

        # prebind locals
        se = self.structure_environments
        v = se.voronoi
        dist_cut = float(self._distance_cutoff)
        ang_cut = float(self._angle_cutoff)

        # compute i_dist (mins increasing)
        mins = np.fromiter((wd["min"] for wd in v.neighbors_normalized_distances[isite]), dtype=float)
        i_dist = int(np.searchsorted(mins, dist_cut, side="right") - 1)

        # compute i_ang (maxs decreasing) using -maxs as increasing
        maxs = np.fromiter((wa["max"] for wa in v.neighbors_normalized_angles[isite]), dtype=float)
        i_ang = int(np.searchsorted(-maxs, -ang_cut, side="right") - 1)

        if i_dist < 0 or i_ang < 0:
            raise ValueError("Distance or angle parameter not found ...")

        # map (i_dist,i_ang)->(cn,inb) from cache
        cn_map = self._get_cn_map_cache(isite).get((i_dist, i_ang))
        if cn_map is None:
            return None

        ce = se.ce_list[se.sites_map[isite]][cn_map[0]][cn_map[1]]
        if ce is None:
            return None
        coord_geoms = ce.coord_geoms

        if return_map:
            if coord_geoms is None:
                return cn_map[0], cn_map
            return (ce.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type), cn_map)

        if coord_geoms is None:
            return cn_map[0]
        return ce.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type)

    def get_site_coordination_environments_fractions(
        self,
        site,
        isite=None,
        dequivsite=None,
        dthissite=None,
        mysym=None,
        ordered=True,
        min_fraction=0,
        return_maps=True,
        return_strategy_dict_info=False,
    ):
        """Get the coordination environments of a given site and additional information.

        Args:
            site: Site for which coordination environment is needed.
            isite: Index of the site.
            dequivsite: Translation of the equivalent site.
            dthissite: Translation of this site.
            mysym: Symmetry to be applied.
            ordered: Whether to order the list by fractions.
            min_fraction: Minimum fraction to include in the list
            return_maps: Whether to return cn_maps (identifies all the NeighborsSet used).
            return_strategy_dict_info: Whether to add the info about the strategy used.

        Returns:
            List of Dict with coordination environment, fraction and additional info.
        """
        if isite is None or dequivsite is None or dthissite is None or mysym is None:
            isite, dequivsite, dthissite, mysym = self.equivalent_site_index_and_transform(site)
        site_nb_sets = self.structure_environments.neighbors_sets[isite]
        if site_nb_sets is None:
            return None
        ce_and_map = self.get_site_coordination_environment(
            site=site,
            isite=isite,
            dequivsite=dequivsite,
            dthissite=dthissite,
            mysym=mysym,
            return_map=True,
        )
        if ce_and_map is None:
            return None
        ce, ce_map = ce_and_map
        if ce is None:
            ce_dict = {
                "ce_symbol": f"UNKNOWN:{ce_map[0]}",
                "ce_dict": None,
                "ce_fraction": 1,
            }
        else:
            ce_dict = {"ce_symbol": ce[0], "ce_dict": ce[1], "ce_fraction": 1}
        if return_maps:
            ce_dict["ce_map"] = ce_map
        if return_strategy_dict_info:
            ce_dict["strategy_info"] = {}
        return [ce_dict]

    def get_site_coordination_environments(
        self,
        site,
        isite=None,
        dequivsite=None,
        dthissite=None,
        mysym=None,
        return_maps=False,
    ):
        """Get the coordination environments of a given site.

        Args:
            site: Site for which coordination environment is needed.
            isite: Index of the site.
            dequivsite: Translation of the equivalent site.
            dthissite: Translation of this site.
            mysym: Symmetry to be applied.
            return_maps: Whether to return cn_maps (identifies all the NeighborsSet used).

        Returns:
            List of coordination environment.
        """
        env = self.get_site_coordination_environment(
            site=site,
            isite=isite,
            dequivsite=dequivsite,
            dthissite=dthissite,
            mysym=mysym,
            return_map=return_maps,
        )
        return [env]

    def add_strategy_visualization_to_subplot(self, subplot, visualization_options=None, plot_type=None):
        """Add a visual of the strategy on a distance-angle plot.

        Args:
            subplot: Axes object onto the visual should be added.
            visualization_options: Options for the visual.
            plot_type: Type of distance-angle plot.
        """
        subplot.plot(
            self._distance_cutoff,
            self._angle_cutoff,
            "o",
            markeredgecolor=None,
            markerfacecolor="w",
            markersize=12,
        )
        subplot.plot(self._distance_cutoff, self._angle_cutoff, "x", linewidth=2, markersize=12)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented

        return (
            self._distance_cutoff == other._distance_cutoff
            and self._angle_cutoff == other._angle_cutoff
            and self._additional_condition == other._additional_condition
            and self._continuous_symmetry_measure_cutoff == other._continuous_symmetry_measure_cutoff
            and self.symmetry_measure_type == other.symmetry_measure_type
        )

    def as_dict(self):
        """
        Bson-serializable dict representation of the SimplestChemenvStrategy object.

        Returns:
            Bson-serializable dict representation of the SimplestChemenvStrategy object.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "distance_cutoff": float(self._distance_cutoff),
            "angle_cutoff": float(self._angle_cutoff),
            "additional_condition": int(self._additional_condition),
            "continuous_symmetry_measure_cutoff": float(self._continuous_symmetry_measure_cutoff),
            "symmetry_measure_type": self._symmetry_measure_type,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Reconstructs the SimplestChemenvStrategy object from a dict representation of the SimplestChemenvStrategy object
        created using the as_dict method.

        Args:
            dct: dict representation of the SimplestChemenvStrategy object

        Returns:
            StructureEnvironments object.
        """
        return cls(
            distance_cutoff=dct["distance_cutoff"],
            angle_cutoff=dct["angle_cutoff"],
            additional_condition=dct["additional_condition"],
            continuous_symmetry_measure_cutoff=dct["continuous_symmetry_measure_cutoff"],
            symmetry_measure_type=dct["symmetry_measure_type"],
        )


class SimpleAbundanceChemenvStrategy(AbstractChemenvStrategy):
    """
    Simple ChemenvStrategy using the neighbors that are the most "abundant" in the grid of angle and distance
    parameters for the definition of neighbors in the Voronoi approach.
    The coordination environment is then given as the one with the lowest continuous symmetry measure.
    """

    DEFAULT_MAX_DIST = 2.0
    DEFAULT_ADDITIONAL_CONDITION = AbstractChemenvStrategy.AC.ONLY_ACB
    STRATEGY_OPTIONS: ClassVar[dict[str, dict]] = {
        "surface_calculation_type": {},
        "additional_condition": {
            "type": AdditionalConditionInt,
            "internal": "_additional_condition",
            "default": DEFAULT_ADDITIONAL_CONDITION,
        },
    }
    STRATEGY_DESCRIPTION = (
        'Simple Abundance ChemenvStrategy using the most "abundant" neighbors map \n'
        "for the definition of neighbors in the Voronoi approach. \n"
        "The coordination environment is then given as the one with the \n"
        "lowest continuous symmetry measure."
    )

    # AL cache key for default surface params
    _DEFAULT_SURFACE_CALC: ClassVar = {
        "distance_parameter": ("initial_normalized", None),
        "angle_parameter": ("initial_normalized", None),
    }

    def __init__(
        self,
        structure_environments=None,
        additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB,
        symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE,
    ):
        """
        Constructor for the SimpleAbundanceChemenvStrategy.

        Args:
            structure_environments: StructureEnvironments object containing all the information on the
                coordination of the sites in a structure.
        """
        raise NotImplementedError("SimpleAbundanceChemenvStrategy not yet implemented")
        AbstractChemenvStrategy.__init__(self, structure_environments, symmetry_measure_type=symmetry_measure_type)
        self._additional_condition = additional_condition

    @property
    def uniquely_determines_coordination_environments(self):
        """Whether this strategy uniquely determines coordination environments."""
        return True

    def get_site_neighbors(self, site):
        isite, dequivsite, dthissite, mysym = self.equivalent_site_index_and_transform(site)
        cn_map = self._get_map(isite)
        eqsite_ps = self.structure_environments.unique_coordinated_neighbors(isite, cn_map=cn_map)
        op = mysym.operate
        de = dequivsite
        dt = dthissite
        return [PeriodicSite(ps._species, op(ps.frac_coords + de) + dt, ps._lattice) for ps in eqsite_ps]

    def get_site_coordination_environment(
        self,
        site,
        isite=None,
        dequivsite=None,
        dthissite=None,
        mysym=None,
        return_map=False,
    ):
        if isite is None:
            isite, *_ = self.equivalent_site_index_and_transform(site)

        cn_map = self._get_map(isite)
        if cn_map is None:
            return None

        se = self.structure_environments
        sm = se.sites_map[isite]
        cn, inb = cn_map
        ce = se.ce_list[sm][cn][inb]

        if return_map:
            if ce is None or ce.coord_geoms is None:
                return cn, cn_map
            return (ce.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type), cn_map)

        if ce is None or ce.coord_geoms is None:
            return cn
        return ce.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type)

    def get_site_coordination_environments(
        self,
        site,
        isite=None,
        dequivsite=None,
        dthissite=None,
        mysym=None,
        return_maps=False,
    ):
        return [
            self.get_site_coordination_environment(
                site=site,
                isite=isite,
                dequivsite=dequivsite,
                dthissite=dthissite,
                mysym=mysym,
                return_map=return_maps,
            )
        ]

    def _get_map(self, isite):
        # cache best map per (isite, additional_condition)
        cache = getattr(self, "_best_map_cache", None)
        if cache is None:
            cache = self._best_map_cache = {}
        key = (isite, int(self._additional_condition))
        hit = cache.get(key)
        if hit is not None:
            return hit

        maps_and_surfaces = self._get_maps_surfaces(isite)
        if not maps_and_surfaces:
            return None

        ac = self._additional_condition
        best_map = None
        best_surface = -1.0
        for entry in maps_and_surfaces:
            # fast skip if AC not present
            if not any(p[2] == ac for p in entry["parameters_indices"]):
                continue
            s = entry["surface"]
            if s > best_surface:
                best_surface = s
                best_map = entry["map"]

        cache[key] = best_map
        return best_map

    def _get_maps_surfaces(self, isite, surface_calculation_type=None):
        # cache maps/surfaces per isite when default calc is used
        if surface_calculation_type is None:
            surface_calculation_type = self._DEFAULT_SURFACE_CALC
            mcache = getattr(self, "_maps_surfaces_cache", None)
            if mcache is None:
                mcache = self._maps_surfaces_cache = {}
            hit = mcache.get(isite)
            if hit is not None:
                return hit
            res = self.structure_environments.voronoi.maps_and_surfaces(
                isite=isite,
                surface_calculation_type=surface_calculation_type,
                max_dist=self.DEFAULT_MAX_DIST,
            )
            mcache[isite] = res
            return res

        # non-default: no caching (params may vary)
        return self.structure_environments.voronoi.maps_and_surfaces(
            isite=isite,
            surface_calculation_type=surface_calculation_type,
            max_dist=self.DEFAULT_MAX_DIST,
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented

        return self._additional_condition == other.additional_condition  # type: ignore[has-type]

    def as_dict(self):
        """
        Bson-serializable dict representation of the SimpleAbundanceChemenvStrategy object.

        Returns:
            Bson-serializable dict representation of the SimpleAbundanceChemenvStrategy object.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "additional_condition": self._additional_condition,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Reconstructs the SimpleAbundanceChemenvStrategy object from a dict representation of the
        SimpleAbundanceChemenvStrategy object created using the as_dict method.

        Args:
            dct: dict representation of the SimpleAbundanceChemenvStrategy object

        Returns:
            StructureEnvironments object.
        """
        return cls(additional_condition=dct["additional_condition"])


class TargetedPenaltiedAbundanceChemenvStrategy(SimpleAbundanceChemenvStrategy):
    """
    Simple ChemenvStrategy using the neighbors that are the most "abundant" in the grid of angle and distance
    parameters for the definition of neighbors in the Voronoi approach, with a bias for a given list of target
    environments. This can be useful in the case of, e.g. connectivity search of some given environment.
    The coordination environment is then given as the one with the lowest continuous symmetry measure.
    """

    DEFAULT_TARGET_ENVIRONMENTS = ("O:6",)

    def __init__(
        self,
        structure_environments=None,
        truncate_dist_ang=True,
        additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB,
        max_nabundant=5,
        target_environments=DEFAULT_TARGET_ENVIRONMENTS,
        target_penalty_type="max_csm",
        max_csm=5.0,
        symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE,
    ):
        """Initialize strategy.

        Not yet implemented.

        Args:
            structure_environments:
            truncate_dist_ang:
            additional_condition:
            max_nabundant:
            target_environments:
            target_penalty_type:
            max_csm:
            symmetry_measure_type:
        """
        raise NotImplementedError("TargetedPenaltiedAbundanceChemenvStrategy not yet implemented")

        super().__init__(
            self,
            structure_environments,
            additional_condition=additional_condition,
            symmetry_measure_type=symmetry_measure_type,
        )
        self.max_nabundant = max_nabundant
        self.target_environments = target_environments
        self.target_penalty_type = target_penalty_type
        self.max_csm = max_csm

    def get_site_coordination_environment(
        self,
        site,
        isite=None,
        dequivsite=None,
        dthissite=None,
        mysym=None,
        return_map=False,
    ):
        if isite is None:
            isite, *_ = self.equivalent_site_index_and_transform(site)

        cn_map = self._get_map(isite)
        if cn_map is None:
            return None

        se = self.structure_environments
        sm = se.sites_map[isite]
        cn, inb = cn_map
        chemical_environments = se.ce_list[sm][cn][inb]

        if return_map:
            # keep original semantics (coord_geoms is None OR len == 0)
            if chemical_environments.coord_geoms is None or len(chemical_environments) == 0:
                return cn_map[0], cn_map
            return (chemical_environments.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type), cn_map)

        if chemical_environments.coord_geoms is None:
            return cn_map[0]
        return chemical_environments.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type)

    def _get_map(self, isite):
        # cache by inputs that affect choice
        cache = getattr(self, "_tpab_map_cache", None)
        if cache is None:
            cache = self._tpab_map_cache = {}
        cache_key = (
            isite,
            int(self._additional_condition),
            int(self.max_nabundant),
            tuple(self.target_environments),
            float(self.max_csm),
            self._symmetry_measure_type,
        )
        hit = cache.get(cache_key)
        if hit is not None:
            return hit

        maps_and_surfaces = SimpleAbundanceChemenvStrategy._get_maps_surfaces(self, isite)
        if not maps_and_surfaces:
            # fall back to parent logic if nothing computed
            res = SimpleAbundanceChemenvStrategy._get_map(self, isite)
            cache[cache_key] = res
            return res

        # precompute target CNs once per target_environments
        tcache = getattr(self, "_target_cns_cache", None)
        if tcache is None:
            tcache = self._target_cns_cache = {}
        tenv_key = tuple(self.target_environments)
        target_cns = tcache.get(tenv_key)
        if target_cns is None:
            # use global ALLCG already constructed at module load
            target_cns = [ALLCG.get_geometry_from_mp_symbol(mp).coordination_number for mp in self.target_environments]
            tcache[tenv_key] = target_cns

        # vectorized top-K by surface (desc)
        n = len(maps_and_surfaces)
        k = min(n, int(self.max_nabundant))
        surfaces = np.fromiter((e["surface"] for e in maps_and_surfaces), dtype=float, count=n)
        top_idx = np.argpartition(surfaces, -k)[-k:]
        top_idx = top_idx[np.argsort(surfaces[top_idx])[::-1]]

        ac = self._additional_condition
        best_map = None
        best_csm = 100.0

        se = self.structure_environments
        sm = se.sites_map[isite]
        smt = self._symmetry_measure_type
        max_csm = float(self.max_csm)

        for ii in top_idx:
            entry = maps_and_surfaces[ii]
            params = entry["parameters_indices"]
            # skip if AC not present
            if not any(p[2] == ac for p in params):
                continue

            my_map = entry["map"]
            cn = my_map[0]
            if cn not in target_cns or cn == 0 or cn > 12:
                continue

            ce = se.ce_list[sm][my_map[0]][my_map[1]]
            cg, cgdict = ce.minimum_geometry(symmetry_measure_type=smt)

            csm = cgdict["symmetry_measure"]
            if cg in self.target_environments and csm <= max_csm and csm < best_csm:
                best_csm = csm
                best_map = my_map

        if best_map is None:
            best_map = SimpleAbundanceChemenvStrategy._get_map(self, isite)

        cache[cache_key] = best_map
        return best_map

    # def get_site_coordination_environment(
    #     self,
    #     site,
    #     isite=None,
    #     dequivsite=None,
    #     dthissite=None,
    #     mysym=None,
    #     return_map=False,
    # ):
    #     """Get the coordination environment of a given site.

    #     Args:
    #         site: Site for which coordination environment is needed.
    #         isite: Index of the site.
    #         dequivsite: Translation of the equivalent site.
    #         dthissite: Translation of this site.
    #         mysym: Symmetry to be applied.
    #         return_map: Whether to return cn_map (identifies the NeighborsSet used).

    #     Returns:
    #         Coordination environment of site.
    #     """
    #     if isite is None:
    #         isite, *_ = self.equivalent_site_index_and_transform(site)
    #     cn_map = self._get_map(isite)
    #     if cn_map is None:
    #         return None
    #     chemical_environments = self.structure_environments.ce_list[self.structure_environments.sites_map[isite]][
    #         cn_map[0]
    #     ][cn_map[1]]
    #     if return_map:
    #         if chemical_environments.coord_geoms is None or len(chemical_environments) == 0:
    #             return cn_map[0], cn_map
    #         return (
    #             chemical_environments.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type),
    #             cn_map,
    #         )

    #     if chemical_environments.coord_geoms is None:
    #         return cn_map[0]
    #     return chemical_environments.minimum_geometry(symmetry_measure_type=self._symmetry_measure_type)

    # def _get_map(self, isite):
    #     maps_and_surfaces = SimpleAbundanceChemenvStrategy._get_maps_surfaces(self, isite)
    #     if maps_and_surfaces is None:
    #         return SimpleAbundanceChemenvStrategy._get_map(self, isite)
    #     current_map = None
    #     current_target_env_csm = 100
    #     surfaces = [map_and_surface["surface"] for map_and_surface in maps_and_surfaces]
    #     order = np.argsort(surfaces)[::-1]
    #     target_cgs = [
    #         AllCoordinationGeometries().get_geometry_from_mp_symbol(mp_symbol)
    #         for mp_symbol in self.target_environments
    #     ]
    #     target_cns = [cg.coordination_number for cg in target_cgs]
    #     for ii in range(min([len(maps_and_surfaces), self.max_nabundant])):
    #         my_map_and_surface = maps_and_surfaces[order[ii]]
    #         my_map = my_map_and_surface["map"]
    #         cn = my_map[0]
    #         if cn not in target_cns or cn > 12 or cn == 0:
    #             continue
    #         all_conditions = [params[2] for params in my_map_and_surface["parameters_indices"]]
    #         if self._additional_condition not in all_conditions:
    #             continue
    #         cg, cgdict = self.structure_environments.ce_list[self.structure_environments.sites_map[isite]][my_map[0]][
    #             my_map[1]
    #         ].minimum_geometry(symmetry_measure_type=self._symmetry_measure_type)
    #         if (
    #             cg in self.target_environments
    #             and cgdict["symmetry_measure"] <= self.max_csm
    #             and cgdict["symmetry_measure"] < current_target_env_csm
    #         ):
    #             current_map = my_map
    #             current_target_env_csm = cgdict["symmetry_measure"]
    #     if current_map is not None:
    #         return current_map
    #     return SimpleAbundanceChemenvStrategy._get_map(self, isite)

    @property
    def uniquely_determines_coordination_environments(self):
        """Whether this strategy uniquely determines coordination environments."""
        return True

    def as_dict(self):
        """
        Bson-serializable dict representation of the TargetedPenaltiedAbundanceChemenvStrategy object.

        Returns:
            Bson-serializable dict representation of the TargetedPenaltiedAbundanceChemenvStrategy object.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "additional_condition": self._additional_condition,
            "max_nabundant": self.max_nabundant,
            "target_environments": self.target_environments,
            "target_penalty_type": self.target_penalty_type,
            "max_csm": self.max_csm,
        }

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented

        return (
            self.additional_condition == other.additional_condition
            and self.max_nabundant == other.max_nabundant
            and self.target_environments == other.target_environments
            and self.target_penalty_type == other.target_penalty_type
            and self.max_csm == other.max_csm
        )

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Reconstructs the TargetedPenaltiedAbundanceChemenvStrategy object from a dict representation of the
        TargetedPenaltiedAbundanceChemenvStrategy object created using the as_dict method.

        Args:
            dct: dict representation of the TargetedPenaltiedAbundanceChemenvStrategy object

        Returns:
            TargetedPenaltiedAbundanceChemenvStrategy object.
        """
        return cls(
            additional_condition=dct["additional_condition"],
            max_nabundant=dct["max_nabundant"],
            target_environments=dct["target_environments"],
            target_penalty_type=dct["target_penalty_type"],
            max_csm=dct["max_csm"],
        )


class NbSetWeight(MSONable, abc.ABC):
    """Abstract base class for neighbor set weight estimations."""

    @abc.abstractmethod
    def as_dict(self):
        """A JSON-serializable dict representation of this neighbors set weight."""

    @abc.abstractmethod
    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        """Get the weight of a given neighbors set.

        Args:
            nb_set: Neighbors set.
            structure_environments: Structure environments used to estimate weight.
            cn_map: Mapping index for this neighbors set.
            additional_info: Additional information.

        Returns:
            float: Weight of the neighbors set.
        """


class AngleNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on the angle."""

    SHORT_NAME = "AngleWeight"

    def __init__(self, aa=1):
        """Initialize AngleNbSetWeight estimator.

        Args:
            aa: Exponent of the angle for the estimator.
        """
        self.aa = aa
        if self.aa == 1:
            self.aw = self.angle_sum
        else:
            self.aw = self.angle_sumn

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        """Get the weight of a given neighbors set.

        Args:
            nb_set: Neighbors set.
            structure_environments: Structure environments used to estimate weight.
            cn_map: Mapping index for this neighbors set.
            additional_info: Additional information.

        Returns:
            float: Weight of the neighbors set.
        """
        return self.aw(nb_set=nb_set)

    @staticmethod
    def angle_sum(nb_set):
        """Sum of all angles in a neighbors set.

        Args:
            nb_set: Neighbors set.

        Returns:
            Sum of solid angles for the neighbors set.
        """
        return np.sum(nb_set.angles) / (4.0 * np.pi)

    def angle_sumn(self, nb_set):
        """Sum of all angles to a given power in a neighbors set.

        Args:
            nb_set: Neighbors set.

        Returns:
            Sum of solid angles to the power aa for the neighbors set.
        """
        return np.power(self.angle_sum(nb_set=nb_set), self.aa)

    def __eq__(self, other: object) -> bool:
        if not hasattr(other, "aa"):
            return NotImplemented
        return self.aa == other.aa

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "aa": self.aa,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Construct AngleNbSetWeight from dict representation."""
        return cls(aa=dct["aa"])


class NormalizedAngleDistanceNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on the normalized angle/distance."""

    SHORT_NAME = "NormAngleDistWeight"

    # AL:
    ## All the branching logic is flatter and easier to follow.
    ## Dictionary lookups make average type resolution O(1) and avoid multiple comparisons.
    ## No redundant double-checking of the same conditions.
    def __init__(self, average_type, aa, bb):
        """Initialize NormalizedAngleDistanceNbSetWeight.

        Args:
            average_type: 'geometric' or 'arithmetic'.
            aa: Exponent for angle values.
            bb: Exponent for distance values.
        """
        # eval function selection
        eval_map = {"geometric": self.gweight, "arithmetic": self.aweight}
        try:
            self.eval = eval_map[average_type]
        except KeyError:
            raise ValueError(f"Average type is {average_type!r} while it should be 'geometric' or 'arithmetic'")

        self.average_type = average_type
        self.aa = aa
        self.bb = bb

        # fda function selection
        if aa == 0 and bb == 1:
            self.fda = self.invdist
        elif aa == 0 and bb != 0:
            self.fda = self.invndist
        elif bb == 0 and aa == 1:
            self.fda = self.ang
        elif bb == 0:
            self.fda = self.angn
        elif aa == 1 and bb == 1:
            self.fda = self.anginvdist
        elif aa == 1:
            self.fda = self.anginvndist
        elif bb == 1:
            self.fda = self.angninvdist
        else:
            self.fda = self.angninvndist

    def __eq__(self, other: object) -> bool:
        needed_attrs = ("average_type", "aa", "bb")
        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented

        return all(getattr(self, attr) == getattr(other, attr) for attr in needed_attrs)

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "average_type": self.average_type,
            "aa": self.aa,
            "bb": self.bb,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize from dict.

        Args:
            dct (dict): Dict representation of NormalizedAngleDistanceNbSetWeight.

        Returns:
            NormalizedAngleDistanceNbSetWeight.
        """
        return cls(average_type=dct["average_type"], aa=dct["aa"], bb=dct["bb"])

    # AL vectorized inner ops
    @staticmethod
    def invdist(nb_set):
        nd = np.asarray(nb_set.normalized_distances, dtype=float)
        return 1.0 / nd

    def invndist(self, nb_set):
        nd = np.asarray(nb_set.normalized_distances, dtype=float)
        return 1.0 / (nd**self.bb)

    @staticmethod
    def ang(nb_set):
        return np.asarray(nb_set.normalized_angles, dtype=float)

    def angn(self, nb_set):
        na = np.asarray(nb_set.normalized_angles, dtype=float)
        return na**self.aa

    @staticmethod
    def anginvdist(nb_set):
        na = np.asarray(nb_set.normalized_angles, dtype=float)
        nd = np.asarray(nb_set.normalized_distances, dtype=float)
        return na / nd

    def anginvndist(self, nb_set):
        na = np.asarray(nb_set.normalized_angles, dtype=float)
        nd = np.asarray(nb_set.normalized_distances, dtype=float)
        return na / (nd**self.bb)

    def angninvdist(self, nb_set):
        na = np.asarray(nb_set.normalized_angles, dtype=float)
        nd = np.asarray(nb_set.normalized_distances, dtype=float)
        return (na**self.aa) / nd

    def angninvndist(self, nb_set):
        na = np.asarray(nb_set.normalized_angles, dtype=float)
        nd = np.asarray(nb_set.normalized_distances, dtype=float)
        return (na**self.aa) / (nd**self.bb)

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        fda_vals = self.fda(nb_set=nb_set)
        return self.eval(fda_list=fda_vals)

    @staticmethod
    def gweight(fda_list):
        """Geometric mean of the weights.

        Args:
            fda_list: List of estimator weights for each neighbor.

        Returns:
            Geometric mean of the weights.
        """
        return gmean(fda_list)

    @staticmethod
    def aweight(fda_list):
        """Standard mean of the weights.

        Args:
            fda_list: List of estimator weights for each neighbor.

        Returns:
            Standard mean of the weights.
        """
        return np.mean(fda_list)


# AL:
# Replaced try/except cache read with a direct dict cache (setdefault + .get) -->  no exception overhead.
# Cached per-site map once (eff_map) and wrote back directly --> fewer nested dict hits.
# Prebound hot callables/attrs (mean_estimator, ce_list_site) --> fewer attribute lookups.
# Used truthiness checks (if not mingeoms) instead of len(...)==0 --> cheaper.
# Tight loop over mingeoms with a single get and threshold check --> fewer dict walks and temporaries.
# Early returns for cached values and empty/None cases --> less work on common paths.


def get_effective_csm(
    nb_set,
    cn_map,
    structure_environments,
    additional_info,
    symmetry_measure_type,
    max_effective_csm,
    effective_csm_estimator_ratio_function,
):
    """Get the effective continuous symmetry measure of a given neighbors set."""
    #  fast cache path without exceptions
    eff_map = additional_info.setdefault("effective_csms", {}).setdefault(nb_set.isite, {})
    cached = eff_map.get(cn_map)
    if cached is not None:
        return cached

    #  locals to reduce attribute lookups in hot path
    ce_list_site = structure_environments.ce_list[nb_set.isite]
    mean_estimator = effective_csm_estimator_ratio_function.mean_estimator

    site_chemenv = ce_list_site[cn_map[0]][cn_map[1]]
    if site_chemenv is None:
        eff = 100
    else:
        mingeoms = site_chemenv.minimum_geometries(
            symmetry_measure_type=symmetry_measure_type,
            max_csm=max_effective_csm,
        )
        if not mingeoms:
            eff = 100
        else:
            # build filtered list once; avoid repeated dict walks
            csms = []
            sm_key = symmetry_measure_type
            for _mp_symbol, ce_dict in mingeoms:
                v = ce_dict["other_symmetry_measures"].get(sm_key)
                if v is not None and v <= max_effective_csm:
                    csms.append(v)
            eff = 100 if not csms else mean_estimator(csms)

    # store in cache and return
    eff_map[cn_map] = eff
    return eff


def set_info(additional_info, field, isite, cn_map, value) -> None:
    """Set additional information for the weights.

    Args:
        additional_info: Additional information.
        field: Type of additional information.
        isite: Index of site to add info.
        cn_map: Mapping index of the neighbors set.
        value: Value of this additional information.
    """
    try:
        additional_info[field][isite][cn_map] = value
    except KeyError:
        try:
            additional_info[field][isite] = {cn_map: value}
        except KeyError:
            additional_info[field] = {isite: {cn_map: value}}


class SelfCSMNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on the Self CSM."""

    SHORT_NAME = "SelfCSMWeight"

    DEFAULT_EFFECTIVE_CSM_ESTIMATOR: ClassVar = {
        "function": "power2_inverse_decreasing",
        "options": {"max_csm": 8.0},
    }
    DEFAULT_WEIGHT_ESTIMATOR: ClassVar = {
        "function": "power2_decreasing_exp",
        "options": {"max_csm": 8.0, "alpha": 1},
    }
    DEFAULT_SYMMETRY_MEASURE_TYPE = "csm_wcs_ctwcc"

    def __init__(
        self,
        effective_csm_estimator=DEFAULT_EFFECTIVE_CSM_ESTIMATOR,
        weight_estimator=DEFAULT_WEIGHT_ESTIMATOR,
        symmetry_measure_type=DEFAULT_SYMMETRY_MEASURE_TYPE,
    ):
        """Initialize SelfCSMNbSetWeight.

        Args:
            effective_csm_estimator: Ratio function used for the effective CSM (comparison between neighbors sets).
            weight_estimator: Weight estimator within a given neighbors set.
            symmetry_measure_type: Type of symmetry measure to be used.
        """
        self.effective_csm_estimator = effective_csm_estimator
        self.effective_csm_estimator_rf = CSMInfiniteRatioFunction.from_dict(effective_csm_estimator)
        self.weight_estimator = weight_estimator
        self.weight_estimator_rf = CSMFiniteRatioFunction.from_dict(weight_estimator)
        self.symmetry_measure_type = symmetry_measure_type
        self.max_effective_csm = self.effective_csm_estimator["options"]["max_csm"]

    # AL: micro-optimization of Python call overhead, the big performance gains will come from
    # making get_effective_csm faster, since that's the heavy part this method calls every time.
    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        # Bind locals to cut attribute lookups in hot path
        evaluate = self.weight_estimator_rf.evaluate
        get_eff_csm = get_effective_csm

        effective_csm = get_eff_csm(
            nb_set=nb_set,
            cn_map=cn_map,
            structure_environments=structure_environments,
            additional_info=additional_info,
            symmetry_measure_type=self.symmetry_measure_type,
            max_effective_csm=self.max_effective_csm,
            effective_csm_estimator_ratio_function=self.effective_csm_estimator_rf,
        )

        weight = evaluate(effective_csm)

        # Avoid function call if no info collection
        if additional_info is not None:
            set_info(
                additional_info=additional_info,
                field="self_csms_weights",
                isite=nb_set.isite,
                cn_map=cn_map,
                value=weight,
            )
        return weight

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented

        return (
            self.effective_csm_estimator == other.effective_csm_estimator
            and self.weight_estimator == other.weight_estimator
            and self.symmetry_measure_type == other.symmetry_measure_type
        )

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "effective_csm_estimator": self.effective_csm_estimator,
            "weight_estimator": self.weight_estimator,
            "symmetry_measure_type": self.symmetry_measure_type,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize from dict.

        Args:
            dct (dict): Dict representation of SelfCSMNbSetWeight.

        Returns:
            SelfCSMNbSetWeight.
        """
        return cls(
            effective_csm_estimator=dct["effective_csm_estimator"],
            weight_estimator=dct["weight_estimator"],
            symmetry_measure_type=dct["symmetry_measure_type"],
        )


class DeltaCSMNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on the differences of CSM."""

    SHORT_NAME = "DeltaCSMWeight"

    DEFAULT_EFFECTIVE_CSM_ESTIMATOR: ClassVar = {
        "function": "power2_inverse_decreasing",
        "options": {"max_csm": 8.0},
    }
    DEFAULT_SYMMETRY_MEASURE_TYPE = "csm_wcs_ctwcc"
    DEFAULT_WEIGHT_ESTIMATOR: ClassVar = {
        "function": "smootherstep",
        "options": {"delta_csm_min": 0.5, "delta_csm_max": 3.0},
    }

    # AL
    def __init__(
        self,
        effective_csm_estimator=DEFAULT_EFFECTIVE_CSM_ESTIMATOR,
        weight_estimator=DEFAULT_WEIGHT_ESTIMATOR,
        delta_cn_weight_estimators=None,
        symmetry_measure_type=DEFAULT_SYMMETRY_MEASURE_TYPE,
    ):
        self.effective_csm_estimator = effective_csm_estimator
        self.effective_csm_estimator_rf = CSMInfiniteRatioFunction.from_dict(effective_csm_estimator)
        self.weight_estimator = weight_estimator
        if self.weight_estimator is not None:
            self.weight_estimator_rf = DeltaCSMRatioFunction.from_dict(weight_estimator)

        self.delta_cn_weight_estimators = delta_cn_weight_estimators
        # build once with a dict comprehension
        self.delta_cn_weight_estimators_rfs = (
            {dcn: DeltaCSMRatioFunction.from_dict(cfg) for dcn, cfg in delta_cn_weight_estimators.items()}
            if delta_cn_weight_estimators is not None
            else {}
        )
        self.symmetry_measure_type = symmetry_measure_type
        self.max_effective_csm = self.effective_csm_estimator["options"]["max_csm"]

    # AL
    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        # prebind hot callables/locals
        get_eff = get_effective_csm
        se = structure_environments
        isite = nb_set.isite
        cn = cn_map[0]
        smt = self.symmetry_measure_type
        max_csm = self.max_effective_csm
        eval_default = self.weight_estimator_rf.evaluate if self.weight_estimator is not None else None
        eval_per_dcn = self.delta_cn_weight_estimators_rfs

        effcsm = get_eff(
            nb_set=nb_set,
            cn_map=cn_map,
            structure_environments=se,
            additional_info=additional_info,
            symmetry_measure_type=smt,
            max_effective_csm=max_csm,
            effective_csm_estimator_ratio_function=self.effective_csm_estimator_rf,
        )

        nb_set_weight = 1.0
        delta_csm = None
        delta_csm_cn_map2 = None

        # iterate only cn2 > cn (cn2 == cn was skipped in original loop)
        for cn2, nb_sets in se.neighbors_sets[isite].items():
            if cn2 <= cn:
                continue

            dcn = cn2 - cn
            rf = eval_per_dcn.get(dcn)
            eval_fn = rf.evaluate if rf is not None else eval_default

            for inb_set2, nb_set2 in enumerate(nb_sets):
                effcsm2 = get_eff(
                    nb_set=nb_set2,
                    cn_map=(cn2, inb_set2),
                    structure_environments=se,
                    additional_info=additional_info,
                    symmetry_measure_type=smt,
                    max_effective_csm=max_csm,
                    effective_csm_estimator_ratio_function=self.effective_csm_estimator_rf,
                )
                this_delta_csm = effcsm2 - effcsm

                # evaluate weight for this competing nb-set
                this_w = eval_fn(this_delta_csm) if eval_fn is not None else 1.0
                if this_w < nb_set_weight:
                    nb_set_weight = this_w
                    delta_csm = this_delta_csm
                    delta_csm_cn_map2 = (cn2, inb_set2)

                    # early-exit if proven worst-case
                    if nb_set_weight == 0.0:
                        if additional_info is not None:
                            set_info(additional_info, "delta_csms", isite, cn_map, delta_csm)
                            set_info(additional_info, "delta_csms_weights", isite, cn_map, nb_set_weight)
                            set_info(additional_info, "delta_csms_cn_map2", isite, cn_map, delta_csm_cn_map2)
                        return 0.0

        # record info once
        if additional_info is not None:
            set_info(additional_info, "delta_csms", isite, cn_map, delta_csm)
            set_info(additional_info, "delta_csms_weights", isite, cn_map, nb_set_weight)
            set_info(additional_info, "delta_csms_cn_map2", isite, cn_map, delta_csm_cn_map2)

        return nb_set_weight

    def __eq__(self, other: object) -> bool:
        needed_attrs = [
            "effective_csm_estimator",
            "weight_estimator",
            "delta_cn_weight_estimators",
            "symmetry_measure_type",
        ]
        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented
        return all(getattr(self, attr) == getattr(other, attr) for attr in needed_attrs)

    @classmethod
    def delta_cn_specifics(
        cls,
        delta_csm_mins=None,
        delta_csm_maxs=None,
        function="smootherstep",
        symmetry_measure_type="csm_wcs_ctwcc",
        effective_csm_estimator=DEFAULT_EFFECTIVE_CSM_ESTIMATOR,
    ):
        """Initialize DeltaCSMNbSetWeight from specific coordination number differences.

        Args:
            delta_csm_mins: Minimums for each coordination number.
            delta_csm_maxs: Maximums for each coordination number.
            function: Ratio function used.
            symmetry_measure_type: Type of symmetry measure to be used.
            effective_csm_estimator: Ratio function used for the effective CSM (comparison between neighbors sets).

        Returns:
            DeltaCSMNbSetWeight.
        """
        if delta_csm_mins is None or delta_csm_maxs is None:
            delta_cn_weight_estimators = {
                dcn: {
                    "function": function,
                    "options": {
                        "delta_csm_min": 0.25 + dcn * 0.25,
                        "delta_csm_max": 5.0 + dcn * 0.25,
                    },
                }
                for dcn in range(1, 13)
            }
        else:
            delta_cn_weight_estimators = {
                dcn: {
                    "function": function,
                    "options": {
                        "delta_csm_min": delta_csm_mins[dcn - 1],
                        "delta_csm_max": delta_csm_maxs[dcn - 1],
                    },
                }
                for dcn in range(1, 13)
            }
        return cls(
            effective_csm_estimator=effective_csm_estimator,
            weight_estimator={
                "function": function,
                "options": {
                    "delta_csm_min": delta_cn_weight_estimators[12]["options"]["delta_csm_min"],
                    "delta_csm_max": delta_cn_weight_estimators[12]["options"]["delta_csm_max"],
                },
            },
            delta_cn_weight_estimators=delta_cn_weight_estimators,
            symmetry_measure_type=symmetry_measure_type,
        )

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "effective_csm_estimator": self.effective_csm_estimator,
            "weight_estimator": self.weight_estimator,
            "delta_cn_weight_estimators": self.delta_cn_weight_estimators,
            "symmetry_measure_type": self.symmetry_measure_type,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize from dict.

        Args:
            dct (dict): Dict representation of DeltaCSMNbSetWeight.

        Returns:
            DeltaCSMNbSetWeight.
        """
        return cls(
            effective_csm_estimator=dct["effective_csm_estimator"],
            weight_estimator=dct["weight_estimator"],
            delta_cn_weight_estimators=(
                {int(dcn): dcn_estimator for dcn, dcn_estimator in dct["delta_cn_weight_estimators"].items()}
                if dct.get("delta_cn_weight_estimators") is not None
                else None
            ),
            symmetry_measure_type=dct["symmetry_measure_type"],
        )


class CNBiasNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on specific biases towards specific coordination numbers."""

    SHORT_NAME = "CNBiasWeight"

    def __init__(self, cn_weights, initialization_options):
        """Initialize CNBiasNbSetWeight.

        Args:
            cn_weights: Weights for each coordination.
            initialization_options: Options for initialization.
        """
        self.cn_weights = cn_weights
        self.initialization_options = initialization_options

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        """Get the weight of a given neighbors set.

        Args:
            nb_set: Neighbors set.
            structure_environments: Structure environments used to estimate weight.
            cn_map: Mapping index for this neighbors set.
            additional_info: Additional information.

        Returns:
            float: Weight of the neighbors set.
        """
        return self.cn_weights[len(nb_set)]

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, NbSetWeight):
            return NotImplemented

        return self.cn_weights == other.cn_weights and self.initialization_options == other.initialization_options

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "cn_weights": {str(cn): cnw for cn, cnw in self.cn_weights.items()},
            "initialization_options": self.initialization_options,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize from dict.

        Args:
            dct (dict): Dict representation of CNBiasNbSetWeight.

        Returns:
            CNBiasNbSetWeight.
        """
        return cls(
            cn_weights={int(cn): cnw for cn, cnw in dct["cn_weights"].items()},
            initialization_options=dct["initialization_options"],
        )

    @classmethod
    def linearly_equidistant(cls, weight_cn1, weight_cn13):
        """Initialize linearly equidistant weights for each coordination.

        Args:
            weight_cn1: Weight of coordination 1.
            weight_cn13: Weight of coordination 13.

        Returns:
            CNBiasNbSetWeight.
        """
        initialization_options = {
            "type": "linearly_equidistant",
            "weight_cn1": weight_cn1,
            "weight_cn13": weight_cn13,
        }
        dw = (weight_cn13 - weight_cn1) / 12.0
        cn_weights = {cn: weight_cn1 + (cn - 1) * dw for cn in range(1, 14)}
        return cls(cn_weights=cn_weights, initialization_options=initialization_options)

    @classmethod
    def geometrically_equidistant(cls, weight_cn1, weight_cn13):
        """Initialize geometrically equidistant weights for each coordination.

        Arge:
            weight_cn1: Weight of coordination 1.
            weight_cn13: Weight of coordination 13.

        Returns:
            CNBiasNbSetWeight.
        """
        initialization_options = {
            "type": "geometrically_equidistant",
            "weight_cn1": weight_cn1,
            "weight_cn13": weight_cn13,
        }
        factor = np.power(float(weight_cn13) / weight_cn1, 1 / 12.0)
        cn_weights = {cn: weight_cn1 * np.power(factor, cn - 1) for cn in range(1, 14)}
        return cls(cn_weights=cn_weights, initialization_options=initialization_options)

    @classmethod
    def explicit(cls, cn_weights):
        """Initialize weights explicitly for each coordination.

        Args:
            cn_weights: Weights for each coordination.

        Returns:
            CNBiasNbSetWeight.
        """
        initialization_options = {"type": "explicit"}
        if set(cn_weights) != set(range(1, 14)):
            raise ValueError("Weights should be provided for CN 1 to 13")
        return cls(cn_weights=cn_weights, initialization_options=initialization_options)

    @classmethod
    def from_description(cls, dct: dict) -> Self:
        """Initialize weights from description.

        Args:
            dct (dict): Dictionary description.

        Returns:
            CNBiasNbSetWeight.
        """
        if dct["type"] == "linearly_equidistant":
            return cls.linearly_equidistant(weight_cn1=dct["weight_cn1"], weight_cn13=dct["weight_cn13"])
        if dct["type"] == "geometrically_equidistant":
            return cls.geometrically_equidistant(weight_cn1=dct["weight_cn1"], weight_cn13=dct["weight_cn13"])
        if dct["type"] == "explicit":
            return cls.explicit(cn_weights=dct["cn_weights"])

        raise RuntimeError("Cannot initialize Weights.")


class DistanceAngleAreaNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on the area in the distance-angle space."""

    SHORT_NAME = "DistAngleAreaWeight"

    AC = AdditionalConditions()
    DEFAULT_SURFACE_DEFINITION: ClassVar = {
        "type": "standard_elliptic",
        "distance_bounds": {"lower": 1.2, "upper": 1.8},
        "angle_bounds": {"lower": 0.1, "upper": 0.8},
    }

    def __init__(
        self,
        weight_type="has_intersection",
        surface_definition=DEFAULT_SURFACE_DEFINITION,
        nb_sets_from_hints="fallback_to_source",
        other_nb_sets="0_weight",
        additional_condition=AC.ONLY_ACB,
        smoothstep_distance=None,
        smoothstep_angle=None,
    ):
        """Initialize CNBiasNbSetWeight.

        Args:
            weight_type: Type of weight.
            surface_definition: Definition of the surface.
            nb_sets_from_hints: How to deal with neighbors sets obtained from "hints".
            other_nb_sets: What to do with other neighbors sets.
            additional_condition: Additional condition to be used.
            smoothstep_distance: Smoothstep distance.
            smoothstep_angle: Smoothstep angle.
        """
        self.weight_type = weight_type
        if weight_type == "has_intersection":
            self.area_weight = self.w_area_has_intersection
        elif weight_type == "has_intersection_smoothstep":
            raise NotImplementedError
        else:
            raise ValueError(f'Weight type is {weight_type!r} while it should be "has_intersection"')
        self.surface_definition = surface_definition
        self.nb_sets_from_hints = nb_sets_from_hints
        self.other_nb_sets = other_nb_sets
        self.additional_condition = additional_condition
        self.smoothstep_distance = smoothstep_distance
        self.smoothstep_angle = smoothstep_angle
        if self.nb_sets_from_hints == "fallback_to_source":
            if self.other_nb_sets == "0_weight":
                self.w_area_intersection_specific = self.w_area_intersection_nbsfh_fbs_onb0
            else:
                raise ValueError('Other nb_sets should be "0_weight"')
        else:
            raise ValueError("Nb_sets from hints should fallback to source")
        lower_and_upper_functions = get_lower_and_upper_f(surface_calculation_options=surface_definition)
        self.dmin = surface_definition["distance_bounds"]["lower"]
        self.dmax = surface_definition["distance_bounds"]["upper"]
        self.amin = surface_definition["angle_bounds"]["lower"]
        self.amax = surface_definition["angle_bounds"]["upper"]
        self.f_lower = lower_and_upper_functions["lower"]
        self.f_upper = lower_and_upper_functions["upper"]

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        """Get the weight of a given neighbors set.

        Args:
            nb_set: Neighbors set.
            structure_environments: Structure environments used to estimate weight.
            cn_map: Mapping index for this neighbors set.
            additional_info: Additional information.

        Returns:
            float: Weight of the neighbors set.
        """
        return self.area_weight(
            nb_set=nb_set,
            structure_environments=structure_environments,
            cn_map=cn_map,
            additional_info=additional_info,
        )

    def w_area_has_intersection(self, nb_set, structure_environments, cn_map, additional_info):
        """Get intersection of the neighbors set area with the surface.

        Args:
            nb_set: Neighbors set.
            structure_environments: Structure environments.
            cn_map: Mapping index of the neighbors set.
            additional_info: Additional information.

        Returns:
            Area intersection between neighbors set and surface.
        """
        return self.w_area_intersection_specific(
            nb_set=nb_set,
            structure_environments=structure_environments,
            cn_map=cn_map,
            additional_info=additional_info,
        )

    def w_area_intersection_nbsfh_fbs_onb0(self, nb_set, structure_environments, cn_map, additional_info):
        """Get intersection of the neighbors set area with the surface.

        Args:
            nb_set: Neighbors set.
            structure_environments: Structure environments.
            cn_map: Mapping index of the neighbors set.
            additional_info: Additional information.

        Returns:
            Area intersection between neighbors set and surface.
        """
        dist_ang_sources = [
            src
            for src in nb_set.sources
            if src["origin"] == "dist_ang_ac_voronoi" and src["ac"] == self.additional_condition
        ]
        if len(dist_ang_sources) > 0:
            for src in dist_ang_sources:
                d1 = src["dp_dict"]["min"]
                d2 = src["dp_dict"]["next"]
                a1 = src["ap_dict"]["next"]
                a2 = src["ap_dict"]["max"]
                if self.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2):
                    return 1
            return 0

        from_hints_sources = [src for src in nb_set.sources if src["origin"] == "nb_set_hints"]
        if len(from_hints_sources) == 0:
            return 0
        if len(from_hints_sources) != 1:
            raise ValueError("Found multiple hints sources for nb_set")
        cn_map_src = from_hints_sources[0]["cn_map_source"]
        nb_set_src = structure_environments.neighbors_sets[nb_set.isite][cn_map_src[0]][cn_map_src[1]]
        dist_ang_sources = [
            src
            for src in nb_set_src.sources
            if src["origin"] == "dist_ang_ac_voronoi" and src["ac"] == self.additional_condition
        ]
        if len(dist_ang_sources) == 0:
            return 0
        for src in dist_ang_sources:
            d1 = src["dp_dict"]["min"]
            d2 = src["dp_dict"]["next"]
            a1 = src["ap_dict"]["next"]
            a2 = src["ap_dict"]["max"]
            if self.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2):
                return 1
        return 0

    def rectangle_crosses_area(self, d1, d2, a1, a2):
        """Whether a given rectangle crosses the area defined by the upper and lower curves.

        Args:
            d1: lower d.
            d2: upper d.
            a1: lower a.
            a2: upper a.
        """
        # Case 1
        if d1 <= self.dmin and d2 <= self.dmin:
            return False
        # Case 6
        if d1 >= self.dmax and d2 >= self.dmax:
            return False
        # Case 2
        if d1 <= self.dmin and d2 <= self.dmax:
            ld2 = self.f_lower(d2)
            return not (a2 <= ld2 or a1 >= self.amax)
        # Case 3
        if d1 <= self.dmin and d2 >= self.dmax:
            return not (a2 <= self.amin or a1 >= self.amax)
        # Case 4
        if self.dmin <= d1 <= self.dmax and self.dmin <= d2 <= self.dmax:
            ld1 = self.f_lower(d1)
            ld2 = self.f_lower(d2)
            if a2 <= ld1 and a2 <= ld2:
                return False
            ud1 = self.f_upper(d1)
            ud2 = self.f_upper(d2)
            return not (a1 >= ud1 and a1 >= ud2)
        # Case 5
        if self.dmin <= d1 <= self.dmax and d2 >= self.dmax:
            ud1 = self.f_upper(d1)
            return not (a1 >= ud1 or a2 <= self.amin)
        raise ValueError("Should not reach this point!")

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented

        return (
            self.weight_type == other.weight_type
            and self.surface_definition == other.surface_definition
            and self.nb_sets_from_hints == other.nb_sets_from_hints
            and self.other_nb_sets == other.other_nb_sets
            and self.additional_condition == other.additional_condition
        )

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "weight_type": self.weight_type,
            "surface_definition": self.surface_definition,
            "nb_sets_from_hints": self.nb_sets_from_hints,
            "other_nb_sets": self.other_nb_sets,
            "additional_condition": self.additional_condition,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize from dict.

        Args:
            dct (dict): Dict representation of DistanceAngleAreaNbSetWeight.

        Returns:
            DistanceAngleAreaNbSetWeight.
        """
        return cls(
            weight_type=dct["weight_type"],
            surface_definition=dct["surface_definition"],
            nb_sets_from_hints=dct["nb_sets_from_hints"],
            other_nb_sets=dct["other_nb_sets"],
            additional_condition=dct["additional_condition"],
        )


class DistancePlateauNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on the distance."""

    SHORT_NAME = "DistancePlateauWeight"

    def __init__(self, distance_function=None, weight_function=None):
        """Initialize DistancePlateauNbSetWeight.

        Args:
            distance_function: Distance function to use.
            weight_function: Ratio function to use.
        """
        if distance_function is None:
            self.distance_function = {"type": "normalized_distance"}
        else:
            self.distance_function = distance_function
        if weight_function is None:
            self.weight_function = {
                "function": "inverse_smootherstep",
                "options": {"lower": 0.2, "upper": 0.4},
            }
        else:
            self.weight_function = weight_function
        self.weight_rf = RatioFunction.from_dict(self.weight_function)

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        """Get the weight of a given neighbors set.

        Args:
            nb_set: Neighbors set.
            structure_environments: Structure environments used to estimate weight.
            cn_map: Mapping index for this neighbors set.
            additional_info: Additional information.

        Returns:
            float: Weight of the neighbors set.
        """
        return self.weight_rf.eval(nb_set.distance_plateau())

    def __eq__(self, other: object) -> bool:
        return isinstance(other, type(self))

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "distance_function": self.distance_function,
            "weight_function": self.weight_function,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize from dict.

        Args:
            dct (dict): Dict representation of DistancePlateauNbSetWeight.

        Returns:
            DistancePlateauNbSetWeight.
        """
        return cls(
            distance_function=dct["distance_function"],
            weight_function=dct["weight_function"],
        )


class AnglePlateauNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on the angle."""

    SHORT_NAME = "AnglePlateauWeight"

    def __init__(self, angle_function=None, weight_function=None):
        """Initialize AnglePlateauNbSetWeight.

        Args:
            angle_function: Angle function to use.
            weight_function: Ratio function to use.
        """
        if angle_function is None:
            self.angle_function = {"type": "normalized_angle"}
        else:
            self.angle_function = angle_function
        if weight_function is None:
            self.weight_function = {
                "function": "inverse_smootherstep",
                "options": {"lower": 0.05, "upper": 0.15},
            }
        else:
            self.weight_function = weight_function
        self.weight_rf = RatioFunction.from_dict(self.weight_function)

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        """Get the weight of a given neighbors set.

        Args:
            nb_set: Neighbors set.
            structure_environments: Structure environments used to estimate weight.
            cn_map: Mapping index for this neighbors set.
            additional_info: Additional information.

        Returns:
            float: Weight of the neighbors set.
        """
        return self.weight_rf.eval(nb_set.angle_plateau())

    def __eq__(self, other: object) -> bool:
        return isinstance(other, type(self))

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "angle_function": self.angle_function,
            "weight_function": self.weight_function,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize from dict.

        Args:
            dct (dict): Dict representation of AnglePlateauNbSetWeight.

        Returns:
            AnglePlateauNbSetWeight.
        """
        return cls(angle_function=dct["angle_function"], weight_function=dct["weight_function"])


class DistanceNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on the distance."""

    SHORT_NAME = "DistanceNbSetWeight"

    def __init__(self, weight_function=None, nbs_source="voronoi"):
        """Initialize DistanceNbSetWeight.

        Args:
            weight_function: Ratio function to use.
            nbs_source: Source of the neighbors.
        """
        if weight_function is None:
            self.weight_function = {
                "function": "smootherstep",
                "options": {"lower": 1.2, "upper": 1.3},
            }
        else:
            self.weight_function = weight_function
        self.weight_rf = RatioFunction.from_dict(self.weight_function)
        if nbs_source not in ["nb_sets", "voronoi"]:
            raise ValueError('"nbs_source" should be one of ["nb_sets", "voronoi"]')
        self.nbs_source = nbs_source

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        """Get the weight of a given neighbors set."""
        cn = cn_map[0]
        isite = nb_set.isite

        se = structure_environments
        voro_site = se.voronoi.voronoi_list2[isite]
        nb_indices = nb_set.site_voronoi_indices  # assumed to be a set

        if self.nbs_source == "nb_sets":
            # collect all competing voronoi indices in one pass
            all_idx = set()
            for cn2, nb_sets in se.neighbors_sets[isite].items():
                if cn2 == cn:
                    continue
                for nb2 in nb_sets:
                    all_idx |= nb2.site_voronoi_indices
        elif self.nbs_source == "voronoi":
            all_idx = set(range(len(voro_site)))
        else:
            raise ValueError('"nbs_source" should be one of ["nb_sets", "voronoi"]')

        # exclude current nb_set indices
        candidates = all_idx - nb_indices
        if not candidates:
            return 1

        # min over generator (no list)
        try:
            mind = min(voro_site[i]["normalized_distance"] for i in candidates)
        except ValueError:
            return 1

        return self.weight_rf.eval(mind)

    def as_dict(self):
        """MSOnable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "weight_function": self.weight_function,
            "nbs_source": self.nbs_source,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize from dict.

        Args:
            dct (dict): Dict representation of DistanceNbSetWeight.

        Returns:
            DistanceNbSetWeight.
        """
        return cls(weight_function=dct["weight_function"], nbs_source=dct["nbs_source"])


class DeltaDistanceNbSetWeight(NbSetWeight):
    """Weight of neighbors set based on the difference of distances."""

    SHORT_NAME = "DeltaDistanceNbSetWeight"

    def __init__(self, weight_function=None, nbs_source="voronoi"):
        """Initialize DeltaDistanceNbSetWeight.

        Args:
            weight_function: Ratio function to use.
            nbs_source: Source of the neighbors.
        """
        if weight_function is None:
            self.weight_function = {
                "function": "smootherstep",
                "options": {"lower": 0.1, "upper": 0.2},
            }
        else:
            self.weight_function = weight_function
        self.weight_rf = RatioFunction.from_dict(self.weight_function)
        if nbs_source not in ["nb_sets", "voronoi"]:
            raise ValueError('"nbs_source" should be one of ["nb_sets", "voronoi"]')
        self.nbs_source = nbs_source

    def weight(self, nb_set, structure_environments, cn_map=None, additional_info=None):
        """Get the weight of a given neighbors set."""
        cn = cn_map[0]
        isite = nb_set.isite
        se = structure_environments
        voro_site = se.voronoi.voronoi_list2[isite]
        exclude = nb_set.site_voronoi_indices  # set
        eval_rf = self.weight_rf.eval

        # find minimum normalized_distance among candidates (excluding current nb_set)
        min_other = float("inf")

        if self.nbs_source == "nb_sets":
            for cn2, nb_sets in se.neighbors_sets[isite].items():
                if cn2 == cn:
                    continue
                for nb2 in nb_sets:
                    for idx in nb2.site_voronoi_indices:
                        if idx not in exclude:
                            d = voro_site[idx]["normalized_distance"]
                            min_other = min(min_other, d)
        elif self.nbs_source == "voronoi":
            for idx, rec in enumerate(voro_site):
                if idx not in exclude:
                    d = rec["normalized_distance"]
                    min_other = min(min_other, d)
        else:
            raise ValueError('"nbs_source" should be one of ["nb_sets", "voronoi"]')

        if min_other is float("inf"):
            return 1

        if len(nb_set) == 0:
            return 0

        nb_set_max = max(nb_set.normalized_distances)
        return eval_rf(min_other - nb_set_max)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, type(self))

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "weight_function": self.weight_function,
            "nbs_source": self.nbs_source,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Initialize from dict.

        Args:
            dct (dict): Dict representation of DeltaDistanceNbSetWeight.

        Returns:
            DeltaDistanceNbSetWeight.
        """
        return cls(weight_function=dct["weight_function"], nbs_source=dct["nbs_source"])


class WeightedNbSetChemenvStrategy(AbstractChemenvStrategy):
    """WeightedNbSetChemenvStrategy."""

    STRATEGY_DESCRIPTION = "    WeightedNbSetChemenvStrategy"
    DEFAULT_CE_ESTIMATOR: ClassVar = {
        "function": "power2_inverse_power2_decreasing",
        "options": {"max_csm": 8.0},
    }

    def __init__(
        self,
        structure_environments=None,
        additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB,
        symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE,
        nb_set_weights=None,
        ce_estimator=DEFAULT_CE_ESTIMATOR,
    ):
        if nb_set_weights is None:
            raise ValueError(f"{nb_set_weights=} must be provided")
        AbstractChemenvStrategy.__init__(self, structure_environments, symmetry_measure_type=symmetry_measure_type)
        self._additional_condition = additional_condition
        self.nb_set_weights = nb_set_weights
        # faster: single pass comprehension
        self.ordered_weights = [{"weight": w, "name": w.SHORT_NAME} for w in nb_set_weights]
        self.ce_estimator = ce_estimator
        self.ce_estimator_ratio_function = CSMInfiniteRatioFunction.from_dict(self.ce_estimator)
        self.ce_estimator_fractions = self.ce_estimator_ratio_function.fractions

    @property
    def uniquely_determines_coordination_environments(self):
        """Whether this strategy uniquely determines coordination environments."""
        return False

    def get_site_coordination_environments_fractions(
        self,
        site,
        isite=None,
        dequivsite=None,
        dthissite=None,
        mysym=None,
        ordered=True,
        min_fraction=0,
        return_maps=True,
        return_strategy_dict_info=False,
        return_all=False,
    ):
        if isite is None or dequivsite is None or dthissite is None or mysym is None:
            isite, dequivsite, dthissite, mysym = self.equivalent_site_index_and_transform(site)

        se = self.structure_environments
        site_nb_sets = se.neighbors_sets[isite]
        if site_nb_sets is None:
            return None

        # build all candidate maps once
        cn_maps = [(cn, inb) for cn, nb_sets in site_nb_sets.items() for inb in range(len(nb_sets))]

        weights_additional_info = {"weights": {isite: {}}}
        weights_isite = weights_additional_info["weights"][isite]

        # cascade weights, pruning zero-weight maps unless return_all
        for wdict in self.ordered_weights:
            weight = wdict["weight"]
            wname = wdict["name"]
            cn_maps_new = []
            for cn_map in cn_maps:
                nb_set = site_nb_sets[cn_map[0]][cn_map[1]]
                w_nb_set = weight.weight(
                    nb_set=nb_set,
                    structure_environments=se,
                    cn_map=cn_map,
                    additional_info=weights_additional_info,
                )
                # store weight
                d = weights_isite.setdefault(cn_map, {})
                d[wname] = w_nb_set
                if return_all or w_nb_set > 0:
                    cn_maps_new.append(cn_map)
            cn_maps = cn_maps_new

        # product of weights per nb_set (avoid list conversion)
        from math import prod

        for d in weights_isite.values():
            d["Product"] = prod(d.values())

        # normalize to fractions
        w_nb_sets = {cm: d["Product"] for cm, d in weights_isite.items()}
        total = float(sum(w_nb_sets.values()))
        if total == 0.0:
            return []

        nb_sets_fractions = {cm: w / total for cm, w in w_nb_sets.items()}
        for cm, d in weights_isite.items():
            d["NbSetFraction"] = nb_sets_fractions[cm]

        ce_symbols, ce_dicts, ce_fractions, ce_dict_fractions, ce_maps = [], [], [], [], []
        site_ce_list = se.ce_list[isite]
        smt = self._symmetry_measure_type

        if return_all:
            for cn_map, nb_set_fraction in nb_sets_fractions.items():
                cn, inb = cn_map
                site_ce_nb_set = site_ce_list[cn][inb]
                if site_ce_nb_set is None:
                    continue
                mingeoms = site_ce_nb_set.minimum_geometries(symmetry_measure_type=smt)
                if not mingeoms:
                    ce_symbols.append(f"UNCLEAR:{cn}")
                    ce_dicts.append(None)
                    ce_fractions.append(nb_set_fraction)
                    d = dict(weights_isite[cn_map])
                    d["CEFraction"] = None
                    d["Fraction"] = nb_set_fraction
                    ce_dict_fractions.append(d)
                    ce_maps.append(cn_map)
                    continue
                csms = [ced["other_symmetry_measures"][smt] for _sym, ced in mingeoms]
                fractions = self.ce_estimator_fractions(csms)
                if fractions is None:
                    ce_symbols.append(f"UNCLEAR:{cn}")
                    ce_dicts.append(None)
                    ce_fractions.append(nb_set_fraction)
                    d = dict(weights_isite[cn_map])
                    d["CEFraction"] = None
                    d["Fraction"] = nb_set_fraction
                    ce_dict_fractions.append(d)
                    ce_maps.append(cn_map)
                else:
                    for (sym, ced), frac in zip(mingeoms, fractions, strict=False):
                        ce_symbols.append(sym)
                        ce_dicts.append(ced)
                        frac_tot = nb_set_fraction * frac
                        ce_fractions.append(frac_tot)
                        d = dict(weights_isite[cn_map])
                        d["CEFraction"] = frac
                        d["Fraction"] = frac_tot
                        ce_dict_fractions.append(d)
                        ce_maps.append(cn_map)
        else:
            for cn_map, nb_set_fraction in nb_sets_fractions.items():
                if nb_set_fraction <= 0:
                    continue
                cn, inb = cn_map
                site_ce_nb_set = site_ce_list[cn][inb]
                mingeoms = site_ce_nb_set.minimum_geometries(symmetry_measure_type=smt)
                if not mingeoms:
                    continue
                csms = [ced["other_symmetry_measures"][smt] for _sym, ced in mingeoms]
                fractions = self.ce_estimator_fractions(csms)
                for (sym, ced), frac in zip(mingeoms, fractions, strict=False):
                    if frac > 0:
                        frac_tot = nb_set_fraction * frac
                        ce_symbols.append(sym)
                        ce_dicts.append(ced)
                        ce_fractions.append(frac_tot)
                        d = dict(weights_isite[cn_map])
                        d["CEFraction"] = frac
                        d["Fraction"] = frac_tot
                        ce_dict_fractions.append(d)
                        ce_maps.append(cn_map)

        indices = np.argsort(ce_fractions)[::-1] if ordered else list(range(len(ce_fractions)))

        fractions_info_list = [
            {"ce_symbol": ce_symbols[i], "ce_dict": ce_dicts[i], "ce_fraction": ce_fractions[i]}
            for i in indices
            if ce_fractions[i] >= min_fraction
        ]

        if return_maps:
            j = 0
            for i in indices:
                if ce_fractions[i] >= min_fraction:
                    fractions_info_list[j]["ce_map"] = ce_maps[i]
                    j += 1

        if return_strategy_dict_info:
            j = 0
            for i in indices:
                if ce_fractions[i] >= min_fraction:
                    fractions_info_list[j]["strategy_info"] = ce_dict_fractions[i]
                    j += 1

        return fractions_info_list

    def get_site_coordination_environment(self, site):
        """Get the coordination environment of a given site.

        Not implemented for this strategy
        """

    def get_site_neighbors(self, site):
        """Get the neighbors of a given site.

        Not implemented for this strategy.
        """

    def get_site_coordination_environments(
        self,
        site,
        isite=None,
        dequivsite=None,
        dthissite=None,
        mysym=None,
        return_maps=False,
    ):
        """Get the coordination environments of a given site.

        Args:
            site: Site for which coordination environment is needed.
            isite: Index of the site.
            dequivsite: Translation of the equivalent site.
            dthissite: Translation of this site.
            mysym: Symmetry to be applied.
            return_maps: Whether to return cn_maps (identifies all the NeighborsSet used).

        Returns:
            List of coordination environment.
        """
        if isite is None or dequivsite is None or dthissite is None or mysym is None:
            isite, dequivsite, dthissite, mysym = self.equivalent_site_index_and_transform(site)
        return [
            self.get_site_coordination_environment(
                site=site,
                isite=isite,
                dequivsite=dequivsite,
                dthissite=dthissite,
                mysym=mysym,
                return_map=return_maps,
            )
        ]

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, WeightedNbSetChemenvStrategy):
            return NotImplemented

        return (
            self._additional_condition == other._additional_condition
            and self.symmetry_measure_type == other.symmetry_measure_type
            and self.nb_set_weights == other.nb_set_weights
            and self.ce_estimator == other.ce_estimator
        )

    def as_dict(self):
        """
        Bson-serializable dict representation of the WeightedNbSetChemenvStrategy object.

        Returns:
            Bson-serializable dict representation of the WeightedNbSetChemenvStrategy object.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "additional_condition": self._additional_condition,
            "symmetry_measure_type": self.symmetry_measure_type,
            "nb_set_weights": [nb_set_weight.as_dict() for nb_set_weight in self.nb_set_weights],
            "ce_estimator": self.ce_estimator,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Reconstructs the WeightedNbSetChemenvStrategy object from a dict representation of the
        WeightedNbSetChemenvStrategy object created using the as_dict method.

        Args:
            dct: dict representation of the WeightedNbSetChemenvStrategy object

        Returns:
            WeightedNbSetChemenvStrategy object.
        """
        return cls(
            additional_condition=dct["additional_condition"],
            symmetry_measure_type=dct["symmetry_measure_type"],
            nb_set_weights=dct["nb_set_weights"],
            ce_estimator=dct["ce_estimator"],
        )


class MultiWeightsChemenvStrategy(WeightedNbSetChemenvStrategy):
    """MultiWeightsChemenvStrategy."""

    STRATEGY_DESCRIPTION = "Multi Weights ChemenvStrategy"
    # STRATEGY_INFO_FIELDS = ['cn_map_surface_fraction', 'cn_map_surface_weight',
    #                         'cn_map_mean_csm', 'cn_map_csm_weight',
    #                         'cn_map_delta_csm', 'cn_map_delta_csms_cn_map2', 'cn_map_delta_csm_weight',
    #                         'cn_map_cn_weight',
    #                         'cn_map_fraction', 'cn_map_ce_fraction', 'ce_fraction']
    DEFAULT_CE_ESTIMATOR: ClassVar = {
        "function": "power2_inverse_power2_decreasing",
        "options": {"max_csm": 8.0},
    }

    def __init__(
        self,
        structure_environments=None,
        additional_condition=AbstractChemenvStrategy.AC.ONLY_ACB,
        symmetry_measure_type=AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE,
        dist_ang_area_weight=None,
        self_csm_weight=None,
        delta_csm_weight=None,
        cn_bias_weight=None,
        angle_weight=None,
        normalized_angle_distance_weight=None,
        ce_estimator=DEFAULT_CE_ESTIMATOR,
    ):
        self._additional_condition = additional_condition
        self.dist_ang_area_weight = dist_ang_area_weight
        self.angle_weight = angle_weight
        self.normalized_angle_distance_weight = normalized_angle_distance_weight
        self.self_csm_weight = self_csm_weight
        self.delta_csm_weight = delta_csm_weight
        self.cn_bias_weight = cn_bias_weight

        # build ordered_weights and nb_sets_weights in one compact pass
        ordered = []
        nb_sets_weights = []
        for w, name in (
            (dist_ang_area_weight, "DistAngArea"),
            (self_csm_weight, "SelfCSM"),
            (delta_csm_weight, "DeltaCSM"),
            (cn_bias_weight, "CNBias"),
            (angle_weight, "Angle"),
            (normalized_angle_distance_weight, "NormalizedAngDist"),
        ):
            if w is not None:
                ordered.append({"weight": w, "name": name})
                nb_sets_weights.append(w)
        self.ordered_weights = ordered

        self.ce_estimator = ce_estimator
        self.ce_estimator_ratio_function = CSMInfiniteRatioFunction.from_dict(self.ce_estimator)
        self.ce_estimator_fractions = self.ce_estimator_ratio_function.fractions

        WeightedNbSetChemenvStrategy.__init__(
            self,
            structure_environments,
            additional_condition=additional_condition,
            symmetry_measure_type=symmetry_measure_type,
            nb_set_weights=nb_sets_weights,
            ce_estimator=ce_estimator,
        )

    @classmethod
    def stats_article_weights_parameters(cls):
        """Initialize strategy used in the statistics article."""
        self_csm_weight = SelfCSMNbSetWeight(
            weight_estimator={
                "function": "power2_decreasing_exp",
                "options": {"max_csm": 8.0, "alpha": 1},
            }
        )
        surface_definition = {
            "type": "standard_elliptic",
            "distance_bounds": {"lower": 1.15, "upper": 2.0},
            "angle_bounds": {"lower": 0.05, "upper": 0.75},
        }
        da_area_weight = DistanceAngleAreaNbSetWeight(
            weight_type="has_intersection",
            surface_definition=surface_definition,
            nb_sets_from_hints="fallback_to_source",
            other_nb_sets="0_weight",
            additional_condition=DistanceAngleAreaNbSetWeight.AC.ONLY_ACB,
        )
        symmetry_measure_type = "csm_wcs_ctwcc"
        delta_weight = DeltaCSMNbSetWeight.delta_cn_specifics()
        bias_weight = angle_weight = nad_weight = None
        return cls(
            dist_ang_area_weight=da_area_weight,
            self_csm_weight=self_csm_weight,
            delta_csm_weight=delta_weight,
            cn_bias_weight=bias_weight,
            angle_weight=angle_weight,
            normalized_angle_distance_weight=nad_weight,
            symmetry_measure_type=symmetry_measure_type,
        )

    @property
    def uniquely_determines_coordination_environments(self):
        """Whether this strategy uniquely determines coordination environments."""
        return False

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented

        return (
            self._additional_condition == other._additional_condition
            and self.symmetry_measure_type == other.symmetry_measure_type
            and self.dist_ang_area_weight == other.dist_ang_area_weight
            and self.self_csm_weight == other.self_csm_weight
            and self.delta_csm_weight == other.delta_csm_weight
            and self.cn_bias_weight == other.cn_bias_weight
            and self.angle_weight == other.angle_weight
            and self.normalized_angle_distance_weight == other.normalized_angle_distance_weight
            and self.ce_estimator == other.ce_estimator
        )

    def as_dict(self):
        """
        Returns:
            Bson-serializable dict representation of the MultiWeightsChemenvStrategy object.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "additional_condition": self._additional_condition,
            "symmetry_measure_type": self.symmetry_measure_type,
            "dist_ang_area_weight": (
                self.dist_ang_area_weight.as_dict() if self.dist_ang_area_weight is not None else None
            ),
            "self_csm_weight": (self.self_csm_weight.as_dict() if self.self_csm_weight is not None else None),
            "delta_csm_weight": (self.delta_csm_weight.as_dict() if self.delta_csm_weight is not None else None),
            "cn_bias_weight": (self.cn_bias_weight.as_dict() if self.cn_bias_weight is not None else None),
            "angle_weight": (self.angle_weight.as_dict() if self.angle_weight is not None else None),
            "normalized_angle_distance_weight": (
                self.normalized_angle_distance_weight.as_dict()
                if self.normalized_angle_distance_weight is not None
                else None
            ),
            "ce_estimator": self.ce_estimator,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Reconstructs the MultiWeightsChemenvStrategy object from a dict representation of the
        MultipleAbundanceChemenvStrategy object created using the as_dict method.

        Args:
            dct: dict representation of the MultiWeightsChemenvStrategy object

        Returns:
            MultiWeightsChemenvStrategy object.
        """
        if dct["normalized_angle_distance_weight"] is not None:
            nad_w = NormalizedAngleDistanceNbSetWeight.from_dict(dct["normalized_angle_distance_weight"])
        else:
            nad_w = None
        return cls(
            additional_condition=dct["additional_condition"],
            symmetry_measure_type=dct["symmetry_measure_type"],
            dist_ang_area_weight=(
                DistanceAngleAreaNbSetWeight.from_dict(dct["dist_ang_area_weight"])
                if dct["dist_ang_area_weight"] is not None
                else None
            ),
            self_csm_weight=(
                SelfCSMNbSetWeight.from_dict(dct["self_csm_weight"]) if dct["self_csm_weight"] is not None else None
            ),
            delta_csm_weight=(
                DeltaCSMNbSetWeight.from_dict(dct["delta_csm_weight"]) if dct["delta_csm_weight"] is not None else None
            ),
            cn_bias_weight=(
                CNBiasNbSetWeight.from_dict(dct["cn_bias_weight"]) if dct["cn_bias_weight"] is not None else None
            ),
            angle_weight=(AngleNbSetWeight.from_dict(dct["angle_weight"]) if dct["angle_weight"] is not None else None),
            normalized_angle_distance_weight=nad_w,
            ce_estimator=dct["ce_estimator"],
        )
