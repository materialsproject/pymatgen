# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Module containing class to create an ion
"""

import re
from copy import deepcopy

import numpy as np
from monty.json import MSONable

from pymatgen.core.composition import Composition
from pymatgen.util.string import formula_double_format, Stringify


class Ion(Composition, MSONable, Stringify):
    """
    Basic ion object. It is just a Composition object with an additional
    variable to store charge.
    The net charge can either be represented as Mn++, or Mn+2, or Mn[2+].
    Note the order of the sign and magnitude in each representation.
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
        Creates Ion from formula.

        :param formula:
        :return: Ion
        """
        charge = 0.0
        f = formula
        m = re.search(r"\[([^\[\]]+)\]", f)
        if m:
            m_chg = re.search(r"([\.\d]*)([+-])", m.group(1))
            if m_chg:
                if m_chg.group(1) != "":
                    charge += float(m_chg.group(1)) * (float(m_chg.group(2) + "1"))
                else:
                    charge += float(m_chg.group(2) + "1")
            f = f.replace(m.group(), "", 1)
        m = re.search(r"\(aq\)", f)
        if m:
            f = f.replace(m.group(), "", 1)
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
    def formula(self):
        """
        Returns a formula string, with elements sorted by electronegativity,
        e.g., Li4 Fe4 P4 O16.
        """
        formula = super().formula
        chg_str = ""
        if self.charge > 0:
            chg_str = " +" + formula_double_format(self.charge, False)
        elif self._charge < 0:
            chg_str = " " + formula_double_format(self.charge, False)
        return formula + chg_str

    @property
    def anonymized_formula(self):
        """
        An anonymized formula. Appends charge to the end
        of anonymized composition
        """
        anon_formula = super().anonymized_formula
        chg = self._charge
        chg_str = ""
        if chg > 0:
            chg_str += "{}{}".format("+", str(int(chg)))
        elif chg < 0:
            chg_str += "{}{}".format("-", str(int(np.abs(chg))))
        return anon_formula + chg_str

    @property
    def reduced_formula(self):
        """
        Returns a reduced formula string with appended charge.
        """
        reduced_formula = super().reduced_formula
        charge = self._charge / self.get_reduced_composition_and_factor()[1]
        if charge > 0:
            if abs(charge) == 1:
                chg_str = "[+]"
            else:
                chg_str = "[" + formula_double_format(charge, False) + "+]"
        elif charge < 0:
            if abs(charge) == 1:
                chg_str = "[-]"
            else:
                chg_str = "[{}-]".format(formula_double_format(abs(charge), False))
        else:
            chg_str = "(aq)"
        return reduced_formula + chg_str

    @property
    def alphabetical_formula(self):
        """
        Returns a reduced formula string with appended charge
        """
        alph_formula = super().alphabetical_formula
        chg_str = ""
        if self.charge > 0:
            chg_str = " +" + formula_double_format(self.charge, False)
        elif self.charge < 0:
            chg_str = " " + formula_double_format(self.charge, False)
        return alph_formula + chg_str

    @property
    def charge(self):
        """
        Charge of the ion
        """
        return self._charge

    def as_dict(self):
        """
        Returns:
            dict with composition, as well as charge
        """
        d = super().as_dict()
        d["charge"] = self.charge
        return d

    @classmethod
    def from_dict(cls, d):
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
    def to_reduced_dict(self):
        """
        Returns:
            dict with element symbol and reduced amount e.g.,
            {"Fe": 2.0, "O":3.0}.
        """
        d = self.composition.to_reduced_dict
        d["charge"] = self.charge
        return d

    @property
    def composition(self):
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
        str_ = super().formula
        if self.charge > 0:
            str_ += "^+" + formula_double_format(self.charge, False)
        elif self._charge < 0:
            str_ += "^" + formula_double_format(self.charge, False)
        return str_
