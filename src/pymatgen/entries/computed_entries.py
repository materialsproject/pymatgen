"""This module implements equivalents of the basic ComputedEntry objects, which
is the basic entity that can be used to perform many analyses. ComputedEntries
contain calculated information, typically from VASP or other electronic
structure codes. For example, ComputedEntries can be used as inputs for phase
diagram analysis.
"""

from __future__ import annotations

import os
from itertools import combinations
from typing import TYPE_CHECKING

import numpy as np
import orjson
from monty.json import MontyDecoder
from scipy.interpolate import interp1d

from pymatgen.core.composition import Composition
from pymatgen.core.entries import *  # noqa: F403
from pymatgen.core.entries import ComputedStructureEntry
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing import Literal, Self

    from pymatgen.analysis.phase_diagram import PhaseDiagram
    from pymatgen.core import Structure

__author__ = "Ryan Kingsbury, Matt McDermott, Shyue Ping Ong, Anubhav Jain"
__copyright__ = "Copyright 2011-2020, The Materials Project"
__version__ = "1.1"
__date__ = "April 2020"

with open(os.path.join(os.path.dirname(__file__), "data/g_els.json"), "rb") as file:
    G_ELEMS = orjson.loads(file.read())
with open(os.path.join(os.path.dirname(__file__), "data/nist_gas_gf.json"), "rb") as file:
    G_GASES = orjson.loads(file.read())


class GibbsComputedStructureEntry(ComputedStructureEntry):  # type: ignore[name-defined]
    """An extension to ComputedStructureEntry which includes the estimated Gibbs
    free energy of formation via a machine-learned model.
    """

    def __init__(
        self,
        structure: Structure,
        formation_enthalpy_per_atom: float,
        temp: float = 300,
        gibbs_model: Literal["SISSO"] = "SISSO",
        composition: Composition | None = None,
        correction: float = 0.0,
        energy_adjustments: list | None = None,
        parameters: dict | None = None,
        data: dict | None = None,
        entry_id: str | None = None,
    ) -> None:
        """
        Args:
            structure (Structure): The pymatgen Structure object of an entry.
            formation_enthalpy_per_atom (float): Formation enthalpy of the entry;
            must be
                calculated using phase diagram construction (eV)
            temp (float): Temperature in Kelvin. If temperature is not selected from
                one of [300, 400, 500, ... 2000 K], then free energies will
                be interpolated. Defaults to 300 K.
            gibbs_model ('SISSO'): Model for Gibbs Free energy. "SISSO", the descriptor
                created by Bartel et al. (2018) -- see reference in documentation, is
                currently the only supported option.
            composition (Composition): The composition of the entry. Defaults to None.
            correction (float): A correction to be applied to the energy. Defaults to 0.
            energy_adjustments (list): A list of energy adjustments to be applied to
                the energy. Defaults to None.
            parameters (dict): An optional dict of parameters associated with
                the entry. Defaults to None.
            data (dict): An optional dict of any additional data associated
                with the entry. Defaults to None.
            entry_id: An optional id to uniquely identify the entry.
        """
        if temp < 300 or temp > 2000:
            raise ValueError("Temperature must be selected from range: [300, 2000] K.")

        integer_formula, _ = structure.composition.get_integer_formula_and_factor()

        self.experimental = False
        if integer_formula in G_GASES:
            self.experimental = True
            if "Experimental" not in str(entry_id):
                entry_id = f"{entry_id} (Experimental)"

        super().__init__(
            structure,
            energy=0,  # placeholder, energy reassigned at end of __init__
            composition=composition,
            correction=correction,
            energy_adjustments=energy_adjustments,
            parameters=parameters,
            data=data,
            entry_id=entry_id,
        )

        self.temp = temp
        self.gibbs_model = gibbs_model
        self.formation_enthalpy_per_atom = formation_enthalpy_per_atom

        self.interpolated = False
        if self.temp % 100:
            self.interpolated = True

        if gibbs_model.lower() == "sisso":
            self.gibbs_fn = self.gf_sisso
        else:
            raise ValueError(f"{gibbs_model} not a valid model. The only currently available model is 'SISSO'.")

        self._energy = self.gibbs_fn()

    @due.dcite(Doi("10.1038/s41467-018-06682-4", "Gibbs free energy SISSO descriptor"))
    def gf_sisso(self) -> float:
        """Gibbs Free Energy of formation as calculated by SISSO descriptor from Bartel
        et al. (2018). Units: eV (not normalized).

        WARNING: This descriptor only applies to solids. The implementation here
        attempts to detect and use downloaded NIST-JANAF data for common
        experimental gases (e.g. CO2) where possible. Note that experimental data is
        only for Gibbs Free Energy of formation, so expt. entries will register as
        having a formation enthalpy of 0.

        Reference: Bartel, C. J., Millican, S. L., Deml, A. M., Rumptz, J. R.,
        Tumas, W., Weimer, A. W., … Holder, A. M. (2018). Physical descriptor for
        the Gibbs energy of inorganic crystalline solids and
        temperature-dependent materials chemistry. Nature Communications, 9(1),
        4168. https://doi.org/10.1038/s41467-018-06682-4

        Returns:
            float: the difference between formation enthalpy (T=0 K, Materials
            Project) and the predicted Gibbs free energy of formation (eV)
        """
        comp = self.composition

        if comp.is_element:
            return 0

        integer_formula, factor = comp.get_integer_formula_and_factor()
        if self.experimental:
            data = G_GASES[integer_formula]

            if self.interpolated:
                g_interp = interp1d([int(t) for t in data], list(data.values()))
                energy = g_interp(self.temp)
            else:
                energy = data[str(self.temp)]

            gibbs_energy = energy * factor
        else:
            n_atoms = len(self.structure)
            vol_per_atom = self.structure.volume / n_atoms
            reduced_mass = self._reduced_mass(self.structure)

            gibbs_energy = (
                comp.num_atoms
                * (self.formation_enthalpy_per_atom + self._g_delta_sisso(vol_per_atom, reduced_mass, self.temp))
                - self._sum_g_i()
            )

        return gibbs_energy

    def _sum_g_i(self) -> float:
        """Sum of the stoichiometrically weighted chemical potentials of the elements
        at specified temperature, as acquired from "g_els.json".

        Returns:
            float: sum of weighted chemical potentials [eV]
        """
        elems = self.composition.get_el_amt_dict()

        if self.interpolated:
            sum_g_i = 0
            for elem, amt in elems.items():
                g_interp = interp1d(
                    [float(t) for t in G_ELEMS],
                    [g_dict[elem] for g_dict in G_ELEMS.values()],
                )
                sum_g_i += amt * g_interp(self.temp)
        else:
            sum_g_i = sum(amt * G_ELEMS[str(int(self.temp))][elem] for elem, amt in elems.items())

        return sum_g_i

    @staticmethod
    def _reduced_mass(structure) -> float:
        """Reduced mass as calculated via Eq. 6 in Bartel et al. (2018).

        Args:
            structure (Structure): The pymatgen Structure object of the entry.

        Returns:
            float: reduced mass (amu)
        """
        reduced_comp = structure.composition.reduced_composition
        n_elems = len(reduced_comp.elements)
        elem_dict = reduced_comp.get_el_amt_dict()

        denominator = (n_elems - 1) * reduced_comp.num_atoms

        all_pairs = combinations(elem_dict.items(), 2)
        mass_sum = 0

        for pair in all_pairs:
            m_i = Composition(pair[0][0]).weight
            m_j = Composition(pair[1][0]).weight
            alpha_i = pair[0][1]
            alpha_j = pair[1][1]

            mass_sum += (alpha_i + alpha_j) * (m_i * m_j) / (m_i + m_j)

        return (1 / denominator) * mass_sum

    @staticmethod
    def _g_delta_sisso(vol_per_atom, reduced_mass, temp) -> float:
        """G^delta as predicted by SISSO-learned descriptor from Eq. (4) in
        Bartel et al. (2018).

        Args:
            vol_per_atom (float): volume per atom [Å^3/atom]
            reduced_mass (float): as calculated with pair-wise sum formula
                [amu]
            temp (float): Temperature [K]

        Returns:
            float: G^delta [eV/atom]
        """
        return (
            (-2.48e-4 * np.log(vol_per_atom) - 8.94e-5 * reduced_mass / vol_per_atom) * temp
            + 0.181 * np.log(temp)
            - 0.882
        )

    @classmethod
    def from_pd(
        cls,
        pd: PhaseDiagram,
        temp: float = 300,
        gibbs_model: Literal["SISSO"] = "SISSO",
    ) -> list[Self]:
        """Constructor method for initializing a list of GibbsComputedStructureEntry
        objects from an existing T = 0 K phase diagram composed of
        ComputedStructureEntry objects, as acquired from a thermochemical database;
        (e.g.. The Materials Project).

        Args:
            pd (PhaseDiagram): T = 0 K phase diagram as created in pymatgen. Must
                contain ComputedStructureEntry objects.
            temp (float): Temperature [K] for estimating Gibbs free energy of formation.
            gibbs_model (str): Gibbs model to use; currently the only option is "SISSO".

        Returns:
            [GibbsComputedStructureEntry]: list of new entries which replace the orig.
                entries with inclusion of Gibbs free energy of formation at the
                specified temperature.
        """
        return [
            cls(
                entry.structure,
                formation_enthalpy_per_atom=pd.get_form_energy_per_atom(entry),
                temp=temp,
                correction=0,
                gibbs_model=gibbs_model,
                data=entry.data,
                entry_id=entry.entry_id,
            )
            for entry in pd.all_entries
            if entry in pd.el_refs.values() or not entry.composition.is_element
        ]

    @classmethod
    def from_entries(
        cls,
        entries: list,
        temp: float = 300,
        gibbs_model: Literal["SISSO"] = "SISSO",
    ) -> list[Self]:
        """Constructor method for initializing GibbsComputedStructureEntry objects from
        T = 0 K ComputedStructureEntry objects, as acquired from a thermochemical
        database e.g. The Materials Project.

        Args:
            entries ([ComputedStructureEntry]): List of ComputedStructureEntry objects,
                as downloaded from The Materials Project API.
            temp (float): Temperature [K] for estimating Gibbs free energy of formation.
            gibbs_model (str): Gibbs model to use; currently the only option is "SISSO".

        Returns:
            list[GibbsComputedStructureEntry]: new entries which replace the orig.
                entries with inclusion of Gibbs free energy of formation at the
                specified temperature.
        """
        from pymatgen.analysis.phase_diagram import PhaseDiagram

        pd = PhaseDiagram(entries)
        return cls.from_pd(pd, temp, gibbs_model)

    def as_dict(self) -> dict:
        """MSONable dict."""
        dct = super().as_dict()
        dct["formation_enthalpy_per_atom"] = self.formation_enthalpy_per_atom
        dct["temp"] = self.temp
        dct["gibbs_model"] = self.gibbs_model
        dct["interpolated"] = self.interpolated
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            GibbsComputedStructureEntry
        """
        dec = MontyDecoder()
        return cls(
            dec.process_decoded(dct["structure"]),
            dct["formation_enthalpy_per_atom"],
            dct["temp"],
            dct["gibbs_model"],
            composition=dct.get("composition"),
            correction=dct["correction"],
            energy_adjustments=[dec.process_decoded(e) for e in dct.get("energy_adjustments", {})],
            parameters={k: dec.process_decoded(v) for k, v in dct.get("parameters", {}).items()},
            data={k: dec.process_decoded(v) for k, v in dct.get("data", {}).items()},
            entry_id=dct.get("entry_id"),
        )

    def __repr__(self) -> str:
        return (
            f"GibbsComputedStructureEntry {self.entry_id} - {self.formula}\n"
            f"Gibbs Free Energy (Formation) = {self.energy:.4f}"
        )
