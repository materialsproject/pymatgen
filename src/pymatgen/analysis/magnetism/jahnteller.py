"""JahnTeller distortion analysis."""

from __future__ import annotations

import math
import os
import warnings
from typing import TYPE_CHECKING, Literal, cast

import numpy as np

from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.local_env import LocalStructOrderParams, get_neighbors_of_site_with_index
from pymatgen.core import Species, get_el_sp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from typing import Any

    from pymatgen.core import Structure

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class JahnTellerAnalyzer:
    """
    Will attempt to classify if structure *may* be Jahn-Teller active.
    Class currently uses datafile of hard-coded common Jahn-Teller
    active ions.
    If structure is annotated with magnetic moments, will estimate
    if structure may be high-spin or low-spin.
    Class aims for more false-positives than false-negatives.
    """

    def __init__(self):
        """Init for JahnTellerAnalyzer."""
        self.spin_configs = {
            "oct": {  # key is number of d electrons
                0: {"high": {"e_g": 0, "t_2g": 0}, "default": "high"},
                1: {"high": {"e_g": 0, "t_2g": 1}, "default": "high"},  # weak J-T
                2: {"high": {"e_g": 0, "t_2g": 2}, "default": "high"},  # weak
                3: {"high": {"e_g": 0, "t_2g": 3}, "default": "high"},  # no J-T
                4: {
                    "high": {"e_g": 1, "t_2g": 3},
                    "low": {"e_g": 0, "t_2g": 4},
                    "default": "high",
                },  # strong high, weak low
                5: {
                    "high": {"e_g": 2, "t_2g": 3},
                    "low": {"e_g": 0, "t_2g": 5},
                    "default": "low",
                },  # no high, weak low
                6: {
                    "high": {"e_g": 2, "t_2g": 4},
                    "low": {"e_g": 0, "t_2g": 6},
                    "default": "high",
                },  # weak high, no low
                7: {
                    "high": {"e_g": 2, "t_2g": 5},
                    "low": {"e_g": 1, "t_2g": 6},
                    "default": "low",
                },  # weak high, strong low
                8: {"high": {"e_g": 2, "t_2g": 6}, "default": "high"},  # no
                9: {"high": {"e_g": 3, "t_2g": 6}, "default": "high"},  # strong
                10: {"high": {"e_g": 4, "t_2g": 6}, "default": "high"},
            },
            "tet": {  # no low spin observed experimentally in tetrahedral, all weak J-T
                0: {"high": {"e": 0, "t_2": 0}, "default": "high"},
                1: {"high": {"e": 1, "t_2": 0}, "default": "high"},
                2: {"high": {"e": 2, "t_2": 0}, "default": "high"},
                3: {"high": {"e": 2, "t_2": 1}, "default": "high"},
                4: {"high": {"e": 2, "t_2": 2}, "default": "high"},
                5: {"high": {"e": 2, "t_2": 3}, "default": "high"},
                6: {"high": {"e": 3, "t_2": 3}, "default": "high"},
                7: {"high": {"e": 4, "t_2": 3}, "default": "high"},
                8: {"high": {"e": 4, "t_2": 4}, "default": "high"},
                9: {"high": {"e": 4, "t_2": 5}, "default": "high"},
                10: {"high": {"e": 4, "t_2": 6}, "default": "high"},
            },
        }

    def get_analysis_and_structure(
        self,
        structure: Structure,
        calculate_valences: bool = True,
        guesstimate_spin: bool = False,
        op_threshold: float = 0.1,
    ) -> tuple[dict, Structure]:
        """Obtain an analysis of a given structure and if it may be Jahn-Teller
        active or not. This is a heuristic, and may give false positives and
        false negatives (false positives are preferred).

        Args:
            structure: input structure
            calculate_valences: whether to attempt to calculate valences or not, structure
                should have oxidation states to perform analysis (Default value = True)
            guesstimate_spin: whether to guesstimate spin state from magnetic moments
                or not, use with caution (Default value = False)
            op_threshold: threshold for order parameter above which to consider site
                to match an octahedral or tetrahedral motif, since Jahn-Teller structures
                can often be
                quite distorted, this threshold is smaller than one might expect

        Returns:
            analysis of structure, with key 'strength' which may be 'none', 'strong',
            'weak', or 'unknown' (Default value = 0.1) and decorated structure
        """
        structure = structure.get_primitive_structure()

        if calculate_valences:
            bva = BVAnalyzer()
            structure = bva.get_oxi_state_decorated_structure(structure)

        # no point testing multiple equivalent sites, doesn't make any difference to analysis
        # but makes returned
        symmetrized_structure = SpacegroupAnalyzer(structure).get_symmetrized_structure()

        # to detect structural motifs of a given site
        op = LocalStructOrderParams(["oct", "tet"])

        # dict of site index to the Jahn-Teller analysis of that site
        jt_sites = []
        non_jt_sites = []

        for indices in symmetrized_structure.equivalent_indices:
            idx = indices[0]
            site = symmetrized_structure[idx]

            # only interested in sites with oxidation states
            if isinstance(site.specie, Species) and site.specie.element.is_transition_metal:
                # get motif around site
                order_params = op.get_order_parameters(symmetrized_structure, idx)

                if order_params[0] > order_params[1] and order_params[0] > op_threshold:
                    motif = "oct"
                    motif_order_parameter = order_params[0]
                elif order_params[1] > op_threshold:
                    motif = "tet"
                    motif_order_parameter = order_params[1]
                else:
                    motif = "unknown"
                    motif_order_parameter = None

                if motif in ["oct", "tet"]:
                    motif = cast("Literal['oct', 'tet']", motif)  # mypy needs help

                    # guess spin of metal ion
                    if guesstimate_spin and "magmom" in site.properties:
                        # estimate if high spin or low spin
                        magmom = site.properties["magmom"]
                        spin_state = self._estimate_spin_state(site.specie, motif, magmom)
                    else:
                        spin_state = "unknown"

                    magnitude = self.get_magnitude_of_effect_from_species(site.specie, spin_state, motif)

                    if magnitude != "none":
                        ligands = get_neighbors_of_site_with_index(structure, idx, approach="min_dist", delta=0.15)
                        ligand_bond_lengths = [ligand.distance(structure[idx]) for ligand in ligands]
                        ligands_species = list({str(ligand.specie) for ligand in ligands})
                        ligand_bond_length_spread = max(ligand_bond_lengths) - min(ligand_bond_lengths)

                        def trim(f):
                            """Avoid storing to unreasonable precision, hurts readability."""
                            return float(f"{f:.4f}")

                        # to be Jahn-Teller active, all ligands have to be the same
                        if len(ligands_species) == 1:
                            jt_sites.append(
                                {
                                    "strength": magnitude,
                                    "motif": motif,
                                    "motif_order_parameter": trim(motif_order_parameter),
                                    "spin_state": spin_state,
                                    "species": str(site.specie),
                                    "ligand": ligands_species[0],
                                    "ligand_bond_lengths": [trim(length) for length in ligand_bond_lengths],
                                    "ligand_bond_length_spread": trim(ligand_bond_length_spread),
                                    "site_indices": indices,
                                }
                            )

                    # store reasons for not being J-T active
                    else:
                        non_jt_sites.append(
                            {
                                "site_indices": indices,
                                "strength": "none",
                                "reason": "Not Jahn-Teller active for this electronic configuration.",
                            }
                        )
                else:
                    non_jt_sites.append(
                        {
                            "site_indices": indices,
                            "strength": "none",
                            "reason": f"{motif=}",
                        }
                    )

        # perform aggregation of all sites
        if jt_sites:
            analysis: dict[str, Any] = {"active": True}
            # if any site could exhibit 'strong' Jahn-Teller effect
            # then mark whole structure as strong
            strong_magnitudes = [site["strength"] == "strong" for site in jt_sites]
            if any(strong_magnitudes):
                analysis["strength"] = "strong"
            else:
                analysis["strength"] = "weak"
            analysis["sites"] = jt_sites
            return analysis, structure
        return {"active": False, "sites": non_jt_sites}, structure

    def get_analysis(
        self,
        structure: Structure,
        calculate_valences: bool = True,
        guesstimate_spin: bool = False,
        op_threshold: float = 0.1,
    ) -> dict:
        """
        Convenience method, uses get_analysis_and_structure method.

        Obtain an analysis of a given structure and if it may be Jahn-Teller
        active or not. This is a heuristic, and may give false positives and
        false negatives (false positives are preferred).

        Args:
            structure: input structure
            calculate_valences: whether to attempt to calculate valences or not, structure
                should have oxidation states to perform analysis (Default value = True)
            guesstimate_spin: whether to guesstimate spin state from magnetic moments
                or not, use with caution (Default value = False)
            op_threshold: threshold for order parameter above which to consider site
                to match an octahedral or tetrahedral motif, since Jahn-Teller structures
                can often be
                quite distorted, this threshold is smaller than one might expect

        Returns:
            analysis of structure, with key 'strength' which may be 'none', 'strong',
            'weak', or 'unknown' (Default value = 0.1)
        """
        return self.get_analysis_and_structure(
            structure,
            calculate_valences=calculate_valences,
            guesstimate_spin=guesstimate_spin,
            op_threshold=op_threshold,
        )[0]

    def is_jahn_teller_active(
        self,
        structure: Structure,
        calculate_valences: bool = True,
        guesstimate_spin: bool = False,
        op_threshold: float = 0.1,
    ) -> bool:
        """
        Convenience method, uses get_analysis_and_structure method.
        Check if a given structure and if it may be Jahn-Teller
        active or not. This is a heuristic, and may give false positives and
        false negatives (false positives are preferred).

        Args:
            structure: input structure
            calculate_valences: whether to attempt to calculate valences or not, structure
                should have oxidation states to perform analysis (Default value = True)
            guesstimate_spin: whether to guesstimate spin state from magnetic moments
                or not, use with caution (Default value = False)
            op_threshold: threshold for order parameter above which to consider site
                to match an octahedral or tetrahedral motif, since Jahn-Teller structures
                can often be
                quite distorted, this threshold is smaller than one might expect

        Returns:
            bool: True if might be Jahn-Teller active, False if not
        """
        active = False

        try:
            analysis = self.get_analysis(
                structure,
                calculate_valences=calculate_valences,
                guesstimate_spin=guesstimate_spin,
                op_threshold=op_threshold,
            )
            active = analysis["active"]
        except Exception as exc:
            warnings.warn(f"Error analyzing {structure.reduced_formula}: {exc}", stacklevel=2)

        return active

    def tag_structure(
        self,
        structure: Structure,
        calculate_valences: bool = True,
        guesstimate_spin: bool = False,
        op_threshold: float = 0.1,
    ) -> Structure:
        """
        Convenience method, uses get_analysis_and_structure method.
        Add a "possible_jt_active" site property on Structure.

        Args:
            structure: input structure
            calculate_valences: whether to attempt to calculate valences or not, structure
                should have oxidation states to perform analysis (Default value = True)
            guesstimate_spin: whether to guesstimate spin state from magnetic moments
                or not, use with caution (Default value = False)
            op_threshold: threshold for order parameter above which to consider site
                to match an octahedral or tetrahedral motif, since Jahn-Teller structures
                can often be
                quite distorted, this threshold is smaller than one might expect

        Returns:
            Decorated Structure, will be in primitive setting.
        """
        try:
            analysis, structure = self.get_analysis_and_structure(
                structure,
                calculate_valences=calculate_valences,
                guesstimate_spin=guesstimate_spin,
                op_threshold=op_threshold,
            )
            jt_sites = [False] * len(structure)
            if analysis["active"]:
                for site in analysis["sites"]:
                    for index in site["site_indices"]:
                        jt_sites[index] = True
                        structure.add_site_property("possible_jt_active", jt_sites)
            return structure
        except Exception as exc:
            warnings.warn(f"Error analyzing {structure.reduced_formula}: {exc}", stacklevel=2)
            return structure

    @staticmethod
    def _get_number_of_d_electrons(species: Species) -> float:
        """Get number of d electrons of a species.

        Args:
            species: Species object

        Returns:
            int: Number of d electrons.
        """
        # TODO: replace with more generic Hund's rule algorithm?

        # taken from get_crystal_field_spin
        elec = species.element.full_electronic_structure
        if len(elec) < 4 or elec[-2][1] != "s" or elec[-1][1] != "d":
            raise AttributeError(f"Invalid element {species.symbol} for crystal field calculation.")
        n_electrons = int(elec[-1][2] + elec[-2][2] - species.oxi_state)  # type: ignore[operator]
        if n_electrons < 0 or n_electrons > 10:
            raise AttributeError(f"Invalid oxidation state {species.oxi_state} for element {species.symbol}")

        return n_electrons

    def get_magnitude_of_effect_from_species(self, species: str | Species, spin_state: str, motif: str) -> str:
        """Get magnitude of Jahn-Teller effect from provided species, spin state and motif.

        Args:
            species: e.g. Fe2+
            spin_state: "high" or "low"
            motif: "oct" or "tet"

        Returns:
            str: "none", "weak" or "strong"
        """
        magnitude = "none"

        sp = get_el_sp(species)

        # has to be Species; we need to know the oxidation state
        if isinstance(sp, Species) and sp.element.is_transition_metal:
            d_electrons = self._get_number_of_d_electrons(sp)

            if motif in self.spin_configs:
                if spin_state not in self.spin_configs[motif][d_electrons]:
                    spin_state = self.spin_configs[motif][d_electrons]["default"]
                spin_config = self.spin_configs[motif][d_electrons][spin_state]
                magnitude = JahnTellerAnalyzer.get_magnitude_of_effect_from_spin_config(motif, spin_config)
        else:
            warnings.warn("No data for this species.", stacklevel=2)

        return magnitude

    @staticmethod
    def get_magnitude_of_effect_from_spin_config(motif: str, spin_config: dict[str, float]) -> str:
        """
        Roughly, the magnitude of Jahn-Teller distortion will be:
        * in octahedral environments, strong if e_g orbitals
        unevenly occupied but weak if t_2g orbitals unevenly
        occupied
        * in tetrahedral environments always weaker.

        Args:
            motif: "oct" or "tet"
            spin_config: dict of 'e' (e_g) and 't' (t2_g) with number of electrons in each state

        Returns:
            str: "none", "weak" or "strong"
        """
        magnitude = "none"
        if motif == "oct":
            e_g = spin_config["e_g"]
            t_2g = spin_config["t_2g"]
            if (e_g % 2 != 0) or (t_2g % 3 != 0):
                magnitude = "weak"
                if e_g % 2 == 1:
                    magnitude = "strong"
        elif motif == "tet":
            e = spin_config["e"]
            t_2 = spin_config["t_2"]
            if (e % 3 != 0) or (t_2 % 2 != 0):
                magnitude = "weak"
        return magnitude

    @staticmethod
    def _estimate_spin_state(
        species: str | Species, motif: Literal["oct", "tet"], known_magmom: float
    ) -> Literal["undefined", "low", "high", "unknown"]:
        """Simple heuristic to estimate spin state. If magnetic moment
        is sufficiently close to that predicted for a given spin state,
        we assign it that state. If we only have data for one spin
        state then that's the one we use (e.g. we assume all tetrahedral
        complexes are high-spin, since this is typically the case).

        Args:
            species: str or Species
            motif ("oct" | "tet"): Tetrahedron or octahedron crystal site coordination
            known_magmom: magnetic moment in Bohr magnetons

        Returns:
            "undefined" (if only one spin state possible), "low", "high" or "unknown"
        """
        mu_so_high = JahnTellerAnalyzer.mu_so(species, motif=motif, spin_state="high")
        mu_so_low = JahnTellerAnalyzer.mu_so(species, motif=motif, spin_state="low")
        if mu_so_high == mu_so_low:
            return "undefined"  # undefined or only one spin state possible
        if mu_so_high is None:
            return "low"
        if mu_so_low is None:
            return "high"
        diff = mu_so_high - mu_so_low
        # WARNING! this heuristic has not been robustly tested or benchmarked
        # using 'diff*0.25' as arbitrary measure, if known magmom is
        # too far away from expected value, we don't try to classify it
        if known_magmom > mu_so_high or math.isclose(mu_so_high, known_magmom, abs_tol=diff * 0.25, rel_tol=0):
            return "high"
        if known_magmom < mu_so_low or math.isclose(mu_so_low, known_magmom, abs_tol=diff * 0.25, rel_tol=0):
            return "low"
        return "unknown"

    @staticmethod
    def mu_so(
        species: str | Species,
        motif: Literal["oct", "tet"],
        spin_state: Literal["high", "low"],
    ) -> float | None:
        """Calculate the spin-only magnetic moment for a given species. Only supports transition metals.

        Args:
            species: Species
            motif ("oct" | "tet"): Tetrahedron or octahedron crystal site coordination
            spin_state ("low" | "high"): Whether the species is in a high or low spin state

        Returns:
            float: Spin-only magnetic moment in Bohr magnetons or None if
                species crystal field not defined
        """
        try:
            sp = get_el_sp(species)
            unpaired_spins = sp.get_crystal_field_spin(coordination=motif, spin_config=spin_state)
            # calculation spin-only magnetic moment for this number of unpaired spins
            return np.sqrt(unpaired_spins * (unpaired_spins + 2))
        except AttributeError:
            return None
