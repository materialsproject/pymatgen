"""Module containing class to create an ion."""

from __future__ import annotations

import re
from copy import deepcopy
from typing import TYPE_CHECKING

from monty.json import MSONable

from pymatgen.core.composition import Composition, reduce_formula
from pymatgen.util.string import Stringify, charge_string, formula_double_format

if TYPE_CHECKING:
    from typing_extensions import Self


class Ion(Composition, MSONable, Stringify):
    """Just a Composition object with an additional variable to store charge.

    The net charge can either be represented as Mn++, Mn+2, Mn[2+], Mn[++], or
    Mn[+2]. Note the order of the sign and magnitude in each representation.
    """

    def __init__(self, composition: Composition, charge: float = 0.0, **kwargs) -> None:
        """Flexible Ion construction, similar to Composition.
        For more information, please see pymatgen.core.Composition.
        """
        super().__init__(composition, **kwargs)
        self._charge = charge

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        if self.composition != other.composition:
            return False
        return self.charge == other.charge

    def __add__(self, other) -> Self:
        """Addition of two ions."""
        new_composition = self.composition + other.composition
        new_charge = self.charge + other.charge
        return type(self)(new_composition, new_charge)

    def __sub__(self, other) -> Self:
        """Subtraction of two ions."""
        new_composition = self.composition - other.composition
        new_charge = self.charge - other.charge
        return type(self)(new_composition, new_charge)

    def __mul__(self, other) -> Self:
        """Multiplication of an Ion with a factor."""
        new_composition = self.composition * other
        new_charge = self.charge * other
        return type(self)(new_composition, new_charge)

    def __hash__(self) -> int:
        return hash((self.composition, self.charge))

    def __str__(self) -> str:
        return self.formula

    def __repr__(self) -> str:
        return f"Ion: {self.formula}"

    @classmethod
    def from_formula(cls, formula: str) -> Self:
        """Create Ion from formula. The net charge can either be represented as
        Mn++, Mn+2, Mn[2+], Mn[++], or Mn[+2]. Note the order of the sign and
        magnitude in each representation.

        Also note that (aq) can be included in the formula, e.g. "NaOH (aq)".

        Args:
            formula (str): The formula to create ion from.

        Returns:
            Ion
        """
        # Strip (aq), if present
        if match := re.search(r"\(aq\)", formula):
            formula = formula.replace(match.group(), "", 1)

        # Check for charge in brackets
        charge = 0.0
        if match := re.search(r"\[([^\[\]]+)\]", formula):
            if m_chg := re.search(r"([\.\d]*)([+-]*)([\.\d]*)", match[1]):
                if m_chg[1] != "":
                    if m_chg[3] != "":
                        raise ValueError("Invalid formula")
                    charge += float(m_chg[1]) * (float(f"{m_chg[2]}1"))
                elif m_chg[3] != "":
                    charge += float(m_chg[3]) * (float(f"{m_chg[2]}1"))
                else:
                    for val in re.findall("[+-]", m_chg[2]):
                        charge += float(f"{val}1")

            formula = formula.replace(match.group(), "", 1)

        # If no brackets, parse trailing +/-
        for m_chg in re.finditer(r"([+-])([\.\d]*)", formula):
            sgn = float(f"{m_chg[1]}1")
            charge += float(m_chg[2]) * sgn if m_chg[2].strip() != "" else sgn
            formula = formula.replace(m_chg.group(), "", 1)

        return cls(Composition(formula), charge)

    @property
    def formula(self) -> str:
        """A formula string with appended charge. The
        charge is written with the sign preceding the magnitude, e.g.
        'Ca1 +2'. Uncharged species have "(aq)" appended, e.g. "O2 (aq)".
        """
        formula = super().formula
        return f"{formula} {charge_string(self.charge, brackets=False)}"

    @property
    def anonymized_formula(self) -> str:
        """An anonymized formula. Appends charge to the end
        of anonymized composition.
        """
        anon_formula = super().anonymized_formula
        chg_str = charge_string(self._charge, brackets=False)
        return anon_formula + chg_str

    def get_reduced_formula_and_factor(
        self,
        iupac_ordering: bool = False,
        hydrates: bool = False,
    ) -> tuple[str, float]:
        """Calculate a reduced formula and factor.

        Similar to Composition.get_reduced_formula_and_factor except that O-H formulas
        receive special handling to differentiate between hydrogen peroxide and OH-.
        Formulas containing HO are written with oxygen first (e.g. 'Fe(OH)2' rather than
        'Fe(HO)2'), and special formulas that apply to solids (e.g. Li2O2 instead of LiO)
        are not used.

        Note that the formula returned by this method does not contain a charge.
        To include charge, use formula or reduced_formula instead.

        Args:
            iupac_ordering (bool, optional): Whether to order the
                formula by the iupac "electronegativity" series, defined in
                Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
                Recommendations 2005)". This ordering effectively follows
                the groups and rows of the periodic table, except the
                Lanthanides, Actinides and hydrogen. Note that polyanions
                will still be determined based on the true electronegativity of
                the elements.
            hydrates: If True (default), attempt to recognize hydrated metal
                complexes and separate out the H2O in the reduced formula.
                For example, Zr(OH)4 becomes ZrO2.2H2O. Applies only to
                Ions containing metals.

        Returns:
            tuple[str, float]: A pretty normalized formula and a multiplicative factor, i.e.,
                H4O4 returns ('H2O2', 2.0).
        """
        all_int = all(abs(val - round(val)) < Composition.amount_tolerance for val in self.values())
        if not all_int:
            return self.formula.replace(" ", ""), 1

        comp = self.composition
        nH2O = 0
        if hydrates:
            # Detect hydrated metal complexes
            nH = comp.get("H", 0)
            nO = comp.get("O", 0)
            if nO > 0 and any(e.is_metal for e in comp):
                nH2O = int(nO) if nH >= 2 * nO else int(nH) // 2
                comp = self.composition - nH2O * Composition("H2O")

        el_amt_dict = {k: int(round(v)) for k, v in comp.get_el_amt_dict().items()}
        formula, factor = reduce_formula(el_amt_dict, iupac_ordering=iupac_ordering)

        # This line checks specifically that the contains an equal amount of O and H. When that is the case,
        # they should be displayed as "OH" rather than "HO".
        if self.composition.get("H") == self.composition.get("O"):
            formula = formula.replace("HO", "OH")
        if nH2O > 0:
            formula += f".{nH2O}H2O"

        # Special handling of peroxide, acetic acid / acetate, and small alcohols
        if formula == "OH" and self.charge == 0:
            formula = "H2O2"
            factor /= 2
        # acetic acid
        elif formula == "H2CO":
            formula = "CH3COOH"
            factor /= 2
        # phosphoric acid system
        elif formula == "PH3O4":
            formula = "H3PO4"
        elif formula == "PHO4":
            formula = "HPO4"
        elif formula == "P(HO2)2":
            formula = "H2PO4"
        # acetate
        elif formula == "H3(CO)2":
            formula = "CH3COO"
        # methanol
        elif formula == "H4CO":
            formula = "CH3OH"
        # ethanol
        elif formula == "H6C2O":
            formula = "C2H5OH"
        # propanol
        elif formula == "H8C3O":
            formula = "C3H7OH"
        # butanol
        elif formula == "H10C4O":
            formula = "C4H9OH"
        elif formula == "O" and factor % 3 == 0:
            formula = "O3"
            factor /= 3
        # ammonia
        elif formula == "H4N":
            formula = "NH4"
        elif formula == "H3N":
            formula = "NH3"
        # methane
        elif formula == "H4C":
            formula = "CH4"
        # thiocyanate
        elif formula == "CSN":
            formula = "SCN"
        # triiodide, nitride, an phosphide
        elif (formula in ["N", "P"] and self.charge == -1) or (formula == "I" and self.charge == 1 / 3):
            formula += "3"
            factor /= 3
        # formate # codespell:ignore
        elif formula == "HCOO":
            formula = "HCO2"
        # oxalate
        elif formula == "CO2":
            formula = "C2O4"
            factor /= 2
        # diatomic gases
        elif formula in {"O", "N", "F", "Cl", "H"} and factor % 2 == 0:
            formula += "2"
            factor /= 2

        return formula, factor

    @property
    def reduced_formula(self) -> str:
        """A reduced formula string with appended charge. The
        charge is placed in brackets with the sign preceding the magnitude, e.g.
        'Ca[+2]'. Uncharged species have "(aq)" appended, e.g. "O2(aq)".
        """
        formula, factor = self.get_reduced_formula_and_factor()
        charge = self._charge / factor
        chg_str = charge_string(charge)
        return formula + chg_str

    @property
    def alphabetical_formula(self) -> str:
        """A formula string, with elements sorted by alphabetically and
        appended charge.
        """
        alph_formula = self.composition.alphabetical_formula
        return f"{alph_formula} {charge_string(self.charge, brackets=False)}"

    @property
    def charge(self) -> float:
        """Charge of the ion."""
        return self._charge

    def as_dict(self) -> dict[str, float]:
        """
        Returns:
            dict with composition, as well as charge.
        """
        dct = super().as_dict()
        dct["charge"] = self.charge
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Generate an ion object from a dict created by as_dict().

        Args:
            dct: {symbol: amount} dict.
        """
        dct_copy = deepcopy(dct)
        charge = dct_copy.pop("charge")
        composition = Composition(dct_copy)
        return cls(composition, charge)

    @property
    def to_reduced_dict(self) -> dict:
        """
        Returns:
            dict with element symbol and reduced amount e.g.
                {"Fe": 2.0, "O":3.0}.
        """
        dct = self.composition.to_reduced_dict
        dct["charge"] = self.charge
        return dct

    @property
    def composition(self) -> Composition:
        """Composition of ion."""
        return Composition(self._data)

    def oxi_state_guesses(
        self,
        oxi_states_override: dict | None = None,
        all_oxi_states: bool = False,
        max_sites: int | None = None,
    ) -> list[dict[str, float]]:
        """Check if the composition is charge-balanced and returns all
        charge-balanced oxidation state combinations. Composition must have
        integer values. Note that more num_atoms in the composition gives
        more degrees of freedom. e.g. if possible oxidation states of
        element X are [2,4] and Y are [-3], then XY is not charge balanced
        but X2Y2 is. Results are returned from most to least probable based
        on ICSD statistics. Use max_sites to improve performance if needed.

        Args:
            oxi_states_override (dict): dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
            all_oxi_states (bool): if True, an element defaults to
                all oxidation states in pymatgen Element.icsd_oxidation_states.
                Otherwise, default is Element.common_oxidation_states. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results.
            max_sites (int): if possible, will reduce Compositions to at most
                this many sites to speed up oxidation state guesses. If the
                composition cannot be reduced to this many sites a ValueError
                will be raised. Set to -1 to just reduce fully. If set to a
                number less than -1, the formula will be fully reduced but a
                ValueError will be thrown if the number of atoms in the reduced
                formula is greater than abs(max_sites).

        Returns:
            A list of dicts - each dict reports an element symbol and average
                oxidation state across all sites in that composition. If the
                composition is not charge balanced, an empty list is returned.
        """
        return self._get_oxi_state_guesses(all_oxi_states, max_sites, oxi_states_override, self.charge)[0]  # type: ignore[return-value]

    def to_pretty_string(self) -> str:
        """Pretty string with proper superscripts."""
        _str = super().reduced_formula
        if val := formula_double_format(self.charge, ignore_ones=False):
            _str += f"^{val:+}"
        return _str
