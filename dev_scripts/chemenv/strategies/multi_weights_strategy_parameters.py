# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Script to visualize the model coordination environments
"""

from __future__ import annotations

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

import copy
import json
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
    AngleNbSetWeight,
    CNBiasNbSetWeight,
    DeltaCSMNbSetWeight,
    DistanceAngleAreaNbSetWeight,
    MultiWeightsChemenvStrategy,
    NormalizedAngleDistanceNbSetWeight,
    SelfCSMNbSetWeight,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
    AllCoordinationGeometries,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import (
    AbstractGeometry,
    LocalGeometryFinder,
)
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

allcg = AllCoordinationGeometries()


class CoordinationEnvironmentMorphing:
    """
    Class to morph a coordination environment into another one.
    """

    def __init__(self, initial_environment_symbol, expected_final_environment_symbol, morphing_description):
        self.initial_environment_symbol = initial_environment_symbol
        self.expected_final_environment_symbol = expected_final_environment_symbol
        self.morphing_description = morphing_description
        self.coordination_geometry = allcg.get_geometry_from_mp_symbol(initial_environment_symbol)
        self.abstract_geometry = AbstractGeometry.from_cg(self.coordination_geometry)

    @classmethod
    def simple_expansion(cls, initial_environment_symbol, expected_final_environment_symbol, neighbors_indices):
        """
        Simple expansion of a coordination environment.

        Args:
            initial_environment_symbol (str): The initial coordination environment symbol.
            expected_final_environment_symbol (str): The expected final coordination environment symbol.
            neighbors_indices (list): The indices of the neighbors to be expanded.

        Returns:
            CoordinationEnvironmentMorphing
        """
        morphing_description = [
            {"ineighbor": i_nb, "site_type": "neighbor", "expansion_origin": "central_site"}
            for i_nb in neighbors_indices
        ]
        return cls(
            initial_environment_symbol=initial_environment_symbol,
            expected_final_environment_symbol=expected_final_environment_symbol,
            morphing_description=morphing_description,
        )

    def figure_fractions(self, weights_options: dict, morphing_factors: Sequence[float] = None) -> None:
        """
        Plot the fractions of the initial and final coordination environments as a function of the morphing factor.

        Args:
            weights_options (dict): The weights options.
            morphing_factors (list): The morphing factors.
        """
        if morphing_factors is None:
            morphing_factors = np.linspace(1.0, 2.0, 21)
        # Set up the local geometry finder
        lgf = LocalGeometryFinder()
        lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)
        # Set up the weights for the MultiWeights strategy
        weights = self.get_weights(weights_options)
        # Set up the strategy
        strategy = MultiWeightsChemenvStrategy(
            dist_ang_area_weight=weights["DistAngArea"],
            self_csm_weight=weights["SelfCSM"],
            delta_csm_weight=weights["DeltaCSM"],
            cn_bias_weight=weights["CNBias"],
            angle_weight=weights["Angle"],
            normalized_angle_distance_weight=weights["NormalizedAngDist"],
        )
        fake_valences = [-1] * (self.coordination_geometry.coordination_number + 1)
        fake_valences[0] = 1
        fractions_initial_environment = np.zeros_like(morphing_factors)
        fractions_final_environment = np.zeros_like(morphing_factors)
        for ii, morphing_factor in enumerate(morphing_factors):
            print(ii)
            struct = self.get_structure(morphing_factor=morphing_factor)
            print(struct)
            # Get the StructureEnvironments
            lgf.setup_structure(structure=struct)
            se = lgf.compute_structure_environments(only_indices=[0], valences=fake_valences)
            strategy.set_structure_environments(structure_environments=se)
            result = strategy.get_site_coordination_environments_fractions(
                site=se.structure[0], isite=0, return_strategy_dict_info=True, return_all=True
            )
            for res in result:
                if res["ce_symbol"] == self.initial_environment_symbol:
                    fractions_initial_environment[ii] = res["ce_fraction"]
                elif res["ce_symbol"] == self.expected_final_environment_symbol:
                    fractions_final_environment[ii] = res["ce_fraction"]

        fig_width_cm = 8.25
        fig_height_cm = 7.0
        fig_width = fig_width_cm / 2.54
        fig_height = fig_height_cm / 2.54

        fig = plt.figure(num=1, figsize=(fig_width, fig_height))
        subplot = fig.add_subplot(111)

        subplot.plot(
            morphing_factors,
            fractions_initial_environment,
            "b-",
            label=self.initial_environment_symbol,
            linewidth=1.5,
        )
        subplot.plot(
            morphing_factors,
            fractions_final_environment,
            "g--",
            label=self.expected_final_environment_symbol,
            linewidth=1.5,
        )

        plt.legend(fontsize=8.0, loc=7)
        plt.show()

    def get_structure(self, morphing_factor):
        lattice = Lattice.cubic(5.0)
        myspecies = ["O"] * (self.coordination_geometry.coordination_number + 1)
        myspecies[0] = "Cu"

        coords = copy.deepcopy(self.abstract_geometry.points_wcs_ctwcc())
        bare_points = self.abstract_geometry.bare_points_with_centre

        for morphing in self.morphing_description:
            if morphing["site_type"] == "neighbor":
                isite = morphing["ineighbor"] + 1
                if morphing["expansion_origin"] == "central_site":
                    origin = bare_points[0]
                vector = bare_points[isite] - origin
                coords[isite] += vector * (morphing_factor - 1.0)
            else:
                raise ValueError(f"Key \"site_type\" is {morphing['site_type']} while it can only be neighbor")

        structure = Structure(lattice=lattice, species=myspecies, coords=coords, coords_are_cartesian=True)
        return structure

    def estimate_parameters(self, dist_factor_min, dist_factor_max, symmetry_measure_type="csm_wcs_ctwcc"):
        only_symbols = [self.initial_environment_symbol, self.expected_final_environment_symbol]
        # Set up the local geometry finder
        lgf = LocalGeometryFinder()
        lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)
        # Get the StructureEnvironments
        fake_valences = [-1] * (self.coordination_geometry.coordination_number + 1)
        fake_valences[0] = 1
        # Get the StructureEnvironments for the structure with dist_factor_min
        struct = self.get_structure(morphing_factor=dist_factor_min)
        lgf.setup_structure(structure=struct)
        se = lgf.compute_structure_environments(only_indices=[0], valences=fake_valences, only_symbols=only_symbols)
        csm_info = se.get_csms(isite=0, mp_symbol=self.initial_environment_symbol)
        if len(csm_info) == 0:
            raise ValueError(f"No csm found for {self.initial_environment_symbol}")
        csm_info.sort(key=lambda x: x["other_symmetry_measures"][symmetry_measure_type])
        csm_initial_min_dist = csm_info[0]["other_symmetry_measures"][symmetry_measure_type]

        csm_info = se.get_csms(isite=0, mp_symbol=self.expected_final_environment_symbol)
        if len(csm_info) == 0:
            raise ValueError(f"No csm found for {self.initial_environment_symbol}")
        csm_info.sort(key=lambda x: x["other_symmetry_measures"][symmetry_measure_type])

        csm_final = csm_info[0]["other_symmetry_measures"][symmetry_measure_type]
        if not np.isclose(csm_final, 0.0, rtol=0.0, atol=1e-10):
            raise ValueError("Final coordination is not perfect !")
        # Get the StructureEnvironments for the structure with dist_factor_max
        struct = self.get_structure(morphing_factor=dist_factor_max)
        lgf.setup_structure(structure=struct)
        se = lgf.compute_structure_environments(only_indices=[0], valences=fake_valences, only_symbols=only_symbols)
        csm_info = se.get_csms(isite=0, mp_symbol=self.initial_environment_symbol)
        if len(csm_info) == 0:
            raise ValueError(f"No csm found for {self.initial_environment_symbol}")
        csm_info.sort(key=lambda x: x["other_symmetry_measures"][symmetry_measure_type])
        csm_initial_max_dist = csm_info[0]["other_symmetry_measures"][symmetry_measure_type]

        csm_info = se.get_csms(isite=0, mp_symbol=self.expected_final_environment_symbol)
        if len(csm_info) == 0:
            raise ValueError(f"No csm found for {self.initial_environment_symbol}")
        csm_info.sort(key=lambda x: x["other_symmetry_measures"][symmetry_measure_type])

        csm_final = csm_info[0]["other_symmetry_measures"][symmetry_measure_type]
        if not np.isclose(csm_final, 0.0, rtol=0.0, atol=1e-10):
            raise ValueError("Final coordination is not perfect !")

        return {"delta_csm_min": csm_initial_min_dist, "self_weight_max_csm": csm_initial_max_dist}

    def get_weights(self, weights_options):
        effective_csm_estimator = {"function": "power2_inverse_decreasing", "options": {"max_csm": 8.0}}

        self_weight_estimator = {
            "function": "power2_decreasing_exp",
            "options": {"max_csm": 5.4230949041608305, "alpha": 1.0},
        }

        self_csm_weight = SelfCSMNbSetWeight(
            effective_csm_estimator=effective_csm_estimator, weight_estimator=self_weight_estimator
        )

        surface_definition = {
            "type": "standard_elliptic",
            "distance_bounds": {"lower": 1.05, "upper": 2.0},
            "angle_bounds": {"lower": 0.05, "upper": 0.95},
        }

        da_area_weight = DistanceAngleAreaNbSetWeight(
            weight_type="has_intersection",
            surface_definition=surface_definition,
            nb_sets_from_hints="fallback_to_source",
            other_nb_sets="0_weight",
            additional_condition=DistanceAngleAreaNbSetWeight.AC.ONLY_ACB,
        )

        weight_estimator = {"function": "smootherstep", "options": {"delta_csm_min": 0.5, "delta_csm_max": 3.0}}

        symmetry_measure_type = "csm_wcs_ctwcc"
        delta_csm_weight = DeltaCSMNbSetWeight(
            effective_csm_estimator=effective_csm_estimator,
            weight_estimator=weight_estimator,
            symmetry_measure_type=symmetry_measure_type,
        )

        bias_weight = CNBiasNbSetWeight.linearly_equidistant(weight_cn1=1.0, weight_cn13=4.0)
        angle_weight = AngleNbSetWeight()

        nad_weight = NormalizedAngleDistanceNbSetWeight(average_type="geometric", aa=1, bb=1)

        weights = {
            "DistAngArea": da_area_weight,
            "SelfCSM": self_csm_weight,
            "DeltaCSM": delta_csm_weight,
            "CNBias": bias_weight,
            "Angle": angle_weight,
            "NormalizedAngDist": nad_weight,
        }

        return weights


if __name__ == "__main__":
    print(
        "+-------------------------------------------------------------+\n"
        "|    Development script of the ChemEnv utility of pymatgen    |\n"
        "| Definition of parameters for the MultiWeightChemenvStrategy |\n"
        "+-------------------------------------------------------------+\n"
    )

    with open("ce_pairs.json") as f:
        ce_pairs = json.load(f)
    self_weight_max_csms = {}
    self_weight_max_csms_per_cn = {}
    allselfmaxcsms = []
    delta_csm_mins = {}
    alldeltacsmmins = []
    all_cn_pairs = []
    for ii in range(1, 14):
        self_weight_max_csms_per_cn[str(ii)] = list()
        for jj in range(ii + 1, 14):
            cn_pair = f"{ii:d}_{jj:d}"
            self_weight_max_csms[cn_pair] = list()
            delta_csm_mins[cn_pair] = list()
            all_cn_pairs.append(cn_pair)
    for ce_pair_dict in ce_pairs:
        ce1 = ce_pair_dict["initial_environment_symbol"]
        ce2 = ce_pair_dict["expected_final_environment_symbol"]
        cn_pair = f"{ce2.split(':')[1]}_{ce1.split(':')[1]}"
        nb_indices = ce_pair_dict["neighbors_indices"]
        mindist = ce_pair_dict["dist_factor_min"]
        maxdist = ce_pair_dict["dist_factor_max"]
        morph = CoordinationEnvironmentMorphing.simple_expansion(
            initial_environment_symbol=ce1, expected_final_environment_symbol=ce2, neighbors_indices=nb_indices
        )
        params = morph.estimate_parameters(dist_factor_min=mindist, dist_factor_max=maxdist)
        print(f"For pair {ce1} to {ce2}, parameters are : ")
        print(params)
        self_weight_max_csms[cn_pair].append(params["self_weight_max_csm"])
        delta_csm_mins[cn_pair].append(params["delta_csm_min"])
        allselfmaxcsms.append(params["self_weight_max_csm"])
        alldeltacsmmins.append(params["delta_csm_min"])
        self_weight_max_csms_per_cn[ce1.split(":")[1]].append(params["self_weight_max_csm"])

    fig = plt.figure(1)
    subplot = fig.add_subplot(111)

    for ipair, cn_pair in enumerate(all_cn_pairs):
        if len(self_weight_max_csms[cn_pair]) == 0:
            continue
        subplot.plot(ipair * np.ones_like(self_weight_max_csms[cn_pair]), self_weight_max_csms[cn_pair], "rx")
        subplot.plot(ipair * np.ones_like(delta_csm_mins[cn_pair]), delta_csm_mins[cn_pair], "b+")

    subplot.set_xticks(range(len(all_cn_pairs)))
    subplot.set_xticklabels(all_cn_pairs, rotation="vertical")
    fig.savefig("self_delta_params.pdf")

    fig2 = plt.figure(2)
    subplot2 = fig2.add_subplot(111)

    for cn in range(1, 14):
        subplot2.plot(
            cn * np.ones_like(self_weight_max_csms_per_cn[str(cn)]), self_weight_max_csms_per_cn[str(cn)], "rx"
        )

    subplot2.set_xticks(range(1, 14))
    fig2.savefig("self_params_per_cn.pdf")
    print(np.mean(allselfmaxcsms))
    print(np.mean(alldeltacsmmins))

    fig3 = plt.figure(3, figsize=(24, 12))
    subplot3 = fig3.add_subplot(111)

    for ipair, cn_pair in enumerate(all_cn_pairs):
        if len(delta_csm_mins[cn_pair]) == 0:
            continue
        subplot3.plot(ipair * np.ones_like(delta_csm_mins[cn_pair]), delta_csm_mins[cn_pair], "b+")

    subplot3.set_xticks(range(len(all_cn_pairs)))
    subplot3.set_xticklabels(all_cn_pairs, rotation="vertical")
    fig3.savefig("delta_params_per_cn_pair.pdf")

    plt.show()
