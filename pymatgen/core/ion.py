# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Module containing class to create an ion
"""

import re
from copy import deepcopy
from typing import Dict, Tuple

from monty.json import MSONable

from pymatgen.core.composition import Composition, reduce_formula
from pymatgen.util.string import Stringify, charge_string, formula_double_format


class Ion(Composition, MSONable, Stringify):
    """
    Ion object. Just a Composition object with an additional variable to store
    charge.

    The net charge can either be represented as Mn++, Mn+2, Mn[2+], Mn[++], or
    Mn[+2]. Note the order of the sign and magnitude in each representation.
    """

    def __init__(self, composition, charge=0.0, properties=None):
        """
        Flexible Ion construction, similar to Composition.
        For more information, please see pymatgen.core.Composition
        """
        super().__init__(composition)
        self._charge = charge

    @classmethod
    def from_formula(cls, formula: str) -> "Ion":
        """
        Creates Ion from formula. The net charge can either be represented as
        Mn++, Mn+2, Mn[2+], Mn[++], or Mn[+2]. Note the order of the sign and
        magnitude in each representation.

        Also note that (aq) can be included in the formula, e.g. "NaOH (aq)".

        :param formula:
        :return: Ion
        """
        charge = 0.0
        f = formula
        # strip (aq), if present
        m = re.search(r"\(aq\)", f)
        if m:
            f = f.replace(m.group(), "", 1)
        # check for charge in brackets
        m = re.search(r"\[([^\[\]]+)\]", f)
        if m:
            m_chg = re.search(r"([\.\d]*)([+-]*)([\.\d]*)", m.group(1))
            if m_chg:
                if m_chg.group(1) != "":
                    if m_chg.group(3) != "":
                        raise ValueError("Invalid formula")
                    charge += float(m_chg.group(1)) * (float(m_chg.group(2) + "1"))
                elif m_chg.group(3) != "":
                    charge += float(m_chg.group(3)) * (float(m_chg.group(2) + "1"))
                else:
                    for i in re.findall("[+-]", m_chg.group(2)):
                        charge += float(i + "1")

            f = f.replace(m.group(), "", 1)

        # if no brackets, parse trailing +/-
        for m_chg in re.finditer(r"([+-])([\.\d]*)", f):
            sign = m_chg.group(1)
            sgn = float(str(sign + "1"))
            if m_chg.group(2).strip() != "":
                charge += float(m_chg.group(2)) * sgn
            else:
                charge += sgn
            f = f.replace(m_chg.group(), "", 1)
        composition = Composition(f)
        return cls(composition, charge)

    @property
    def formula(self) -> str:
        """
        Returns a formula string, with elements sorted by electronegativity,
        e.g., Li4 Fe4 P4 O16.
        """
        formula = super().formula
        return formula + " " + charge_string(self.charge, brackets=False)

    @property
    def anonymized_formula(self) -> str:
        """
        An anonymized formula. Appends charge to the end
        of anonymized composition
        """
        anon_formula = super().anonymized_formula
        chg_str = charge_string(self._charge, brackets=False)
        return anon_formula + chg_str

    def get_reduced_formula_and_factor(self, iupac_ordering: bool = False, hydrates: bool = True) -> Tuple[str, float]:
        """
        Calculates a reduced formula and factor.

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
                Lanthanides, Actanides and hydrogen. Note that polyanions
                will still be determined based on the true electronegativity of
                the elements.
            hydrates: If True (default), attempt to recognize hydrated metal
                complexes and separate out the H2O in the reduced formula.
                For example, Zr(OH)4 becomes ZrO2.2H2O. Applies only to
                Ions containing metals.

        Returns:
            A pretty normalized formula and a multiplicative factor, i.e.,
            H4O4 returns ('H2O2', 2.0).
        """
        all_int = all(abs(x - round(x)) < Composition.amount_tolerance for x in self.values())
        if not all_int:
            return self.formula.replace(" ", ""), 1

        comp = self.composition
        nH2O = 0
        if hydrates:
            # detect hydrated metal complexes
            nH = comp.get("H", 0)
            nO = comp.get("O", 0)
            if nO > 0 and any(e.is_metal for e in comp):
                nH2O = int(nO) if nH >= 2 * nO else int(nH) // 2
                comp = self.composition - nH2O * Composition("H2O")

        d = {k: int(round(v)) for k, v in comp.get_el_amt_dict().items()}
        (formula, factor) = reduce_formula(d, iupac_ordering=iupac_ordering)

        if "HO" in formula:
            formula = formula.replace("HO", "OH")

        if nH2O > 0:
            formula += f".{nH2O}H2O"

        # special handling of peroxide, acetic acid / acetate, and small alcohols
        if formula == "OH" and self.charge == 0:
            formula = "H2O2"
            factor /= 2
        # acetic acid
        elif formula == "H2CO":
            formula = "CH3COOH"
            factor /= 2
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

        return formula, factor

    @property
    def reduced_formula(self) -> str:
        """
        Returns a reduced formula string with appended charge. The
        charge is placed in brackets with the sign preceding the magnitude, e.g.,
        'Ca[+2]'.
        """
        reduced_formula = super().reduced_formula
        charge = self._charge / self.get_reduced_composition_and_factor()[1]
        chg_str = charge_string(charge)
        return reduced_formula + chg_str

    @property
    def alphabetical_formula(self) -> str:
        """
        Returns a formula string, with elements sorted by alphabetically and
        appended charge
        """
        alph_formula = self.composition.alphabetical_formula
        return alph_formula + " " + charge_string(self.charge, brackets=False)

    @property
    def charge(self) -> float:
        """
        Charge of the ion
        """
        return self._charge

    def as_dict(self) -> Dict[str, float]:
        """
        Returns:
            dict with composition, as well as charge
        """
        d = super().as_dict()
        d["charge"] = self.charge
        return d

    @classmethod
    def from_dict(cls, d) -> "Ion":
        """
        Generates an ion object from a dict created by as_dict().

        Args:
            d:
                {symbol: amount} dict.
        """
        input = deepcopy(d)
        charge = input.pop("charge")
        composition = Composition(input)
        return Ion(composition, charge)

    @property
    def to_reduced_dict(self) -> dict:
        """
        Returns:
            dict with element symbol and reduced amount e.g.,
            {"Fe": 2.0, "O":3.0}.
        """
        d = self.composition.to_reduced_dict
        d["charge"] = self.charge
        return d

    @property
    def composition(self) -> Composition:
        """Composition of ion."""
        return Composition(self._data)

    def __eq__(self, other):
        if self.composition != other.composition:
            return False
        if self.charge != other.charge:
            return False
        return True

    def __add__(self, other):
        """
        Addition of two ions.
        """
        new_composition = self.composition + other.composition
        new_charge = self.charge + other.charge
        return Ion(new_composition, new_charge)

    def __sub__(self, other):
        """
        Subtraction of two ions
        """
        new_composition = self.composition - other.composition
        new_charge = self.charge - other.charge
        return Ion(new_composition, new_charge)

    def __mul__(self, other):
        """
        Multiplication of an Ion with a factor
        """
        new_composition = self.composition * other
        new_charge = self.charge * other
        return Ion(new_composition, new_charge)

    def __hash__(self):
        return hash((self.composition, self.charge))

    def __str__(self):
        return self.formula

    def __repr__(self):
        return "Ion: " + self.formula

    def to_pretty_string(self) -> str:
        """
        :return: Pretty string with proper superscripts.
        """
        str_ = super().reduced_formula
        if self.charge > 0:
            str_ += "^+" + formula_double_format(self.charge, False)
        elif self._charge < 0:
            str_ += "^" + formula_double_format(self.charge, False)
        return str_
