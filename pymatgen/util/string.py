# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides utility classes for string operations.
"""
from __future__ import annotations

import re
from fractions import Fraction

SUBSCRIPT_UNICODE = {
    "0": "₀",
    "1": "₁",
    "2": "₂",
    "3": "₃",
    "4": "₄",
    "5": "₅",
    "6": "₆",
    "7": "₇",
    "8": "₈",
    "9": "₉",
}

SUPERSCRIPT_UNICODE = {
    "0": "⁰",
    "1": "¹",
    "2": "²",
    "3": "³",
    "4": "⁴",
    "5": "⁵",
    "6": "⁶",
    "7": "⁷",
    "8": "⁸",
    "9": "⁹",
    "+": "⁺",
    "-": "⁻",
}

# TODO: make standalone functions in this module use the same implementation as Stringify
# Note: previous deprecations of standalone functions in this module were removed due to
# a community need.


class Stringify:
    """
    Mix-in class for string formatting, e.g. superscripting numbers and symbols or superscripting.
    """

    STRING_MODE = "SUBSCRIPT"

    def to_pretty_string(self) -> str:
        """
        :return: A pretty string representation. By default, the __str__ output is used, but this method can be
            overridden if a different representation from default is desired.
        """
        return str(self)

    def to_latex_string(self) -> str:
        """
        Generates a LaTeX formatted string. The mode is set by the class variable STRING_MODE, which defaults to
        "SUBSCRIPT". E.g., Fe2O3 is transformed to Fe$_{2}$O$_{3}$. Setting STRING_MODE to "SUPERSCRIPT" creates
        superscript, e.g., Fe2+ becomes Fe^{2+}. The initial string is obtained from the class's __str__ method.

        :return: String for display as in LaTeX with proper superscripts and subscripts.
        """
        str_ = self.to_pretty_string()
        # First we process strings that already have _ and ^ by escaping the relevant parts.
        str_ = re.sub(r"_(\d+)", r"$_{\1}$", str_)
        str_ = re.sub(r"\^([\d\+\-]+)", r"$^{\1}$", str_)
        if self.STRING_MODE == "SUBSCRIPT":
            return re.sub(r"([A-Za-z\(\)])([\d\+\-\.]+)", r"\1$_{\2}$", str_)
        if self.STRING_MODE == "SUPERSCRIPT":
            return re.sub(r"([A-Za-z\(\)])([\d\+\-\.]+)", r"\1$^{\2}$", str_)
        return str_

    def to_html_string(self) -> str:
        """
        Generates a HTML formatted string. This uses the output from to_latex_string to generate a HTML output.
        :return: HTML formatted string.
        """
        str_ = re.sub(r"\$_\{([^}]+)\}\$", r"<sub>\1</sub>", self.to_latex_string())
        str_ = re.sub(r"\$\^\{([^}]+)\}\$", r"<sup>\1</sup>", str_)
        return re.sub(r"\$\\overline\{([^}]+)\}\$", r'<span style="text-decoration:overline">\1</span>', str_)

    def to_unicode_string(self):
        """
        :return: Unicode string with proper sub and superscripts. Note that this works only with systems where the sub
            and superscripts are pure integers.
        """
        str_ = self.to_latex_string()
        for m in re.finditer(r"\$_\{(\d+)\}\$", str_):
            s1 = m.group()
            s2 = [SUBSCRIPT_UNICODE[s] for s in m.group(1)]
            str_ = str_.replace(s1, "".join(s2))
        for m in re.finditer(r"\$\^\{([\d\+\-]+)\}\$", str_):
            s1 = m.group()
            s2 = [SUPERSCRIPT_UNICODE[s] for s in m.group(1)]
            str_ = str_.replace(s1, "".join(s2))
        return str_


def str_delimited(results, header=None, delimiter="\t"):
    """
    Given a tuple of tuples, generate a delimited string form.
    >>> results = [["a","b","c"],["d","e","f"],[1,2,3]]
    >>> print(str_delimited(results,delimiter=","))
    a,b,c
    d,e,f
    1,2,3

    Args:
        result: 2d sequence of arbitrary types.
        header: optional header

    Returns:
        Aligned string output in a table-like format.
    """
    returnstr = ""
    if header is not None:
        returnstr += delimiter.join(header) + "\n"
    return returnstr + "\n".join(delimiter.join([str(m) for m in result]) for result in results)


def formula_double_format(afloat, ignore_ones=True, tol: float = 1e-8):
    """
    This function is used to make pretty formulas by formatting the amounts.
    Instead of Li1.0 Fe1.0 P1.0 O4.0, you get LiFePO4.

    Args:
        afloat (float): a float
        ignore_ones (bool): if true, floats of 1 are ignored.
        tol (float): Tolerance to round to nearest int. i.e. 2.0000000001 -> 2

    Returns:
        A string representation of the float for formulas.
    """
    if ignore_ones and afloat == 1:
        return ""
    if abs(afloat - int(afloat)) < tol:
        return str(int(afloat))
    return str(round(afloat, 8))


def charge_string(charge, brackets=True, explicit_one=True):
    """
    Returns a string representing the charge of an Ion. By default, the
    charge is placed in brackets with the sign preceding the magnitude, e.g.,
    '[+2]'. For uncharged species, the string returned is '(aq)'

    Args:
        brackets: whether to enclose the charge in brackets, e.g. [+2]. Default: True
        explicit_one: whether to include the number one for monovalent ions, e.g.
            +1 rather than +. Default: True
    """
    if charge > 0:
        chg_str = f"+{formula_double_format(charge, False)}"
    elif charge < 0:
        chg_str = f"-{formula_double_format(abs(charge), False)}"
    else:
        chg_str = "(aq)"

    if chg_str in ["+1", "-1"] and not explicit_one:
        chg_str = chg_str.replace("1", "")

    if chg_str != "(aq)" and brackets:
        chg_str = "[" + chg_str + "]"

    return chg_str


def latexify(formula):
    """
    Generates a LaTeX formatted formula. E.g., Fe2O3 is transformed to
    Fe$_{2}$O$_{3}$.

    Note that Composition now has a to_latex_string() method that may
    be used instead.

    Args:
        formula (str): Input formula.

    Returns:
        Formula suitable for display as in LaTeX with proper subscripts.
    """
    return re.sub(r"([A-Za-z\(\)])([\d\.]+)", r"\1$_{\2}$", formula)


def htmlify(formula):
    """
    Generates a HTML formatted formula, e.g. Fe2O3 is transformed to
    Fe<sub>2</sub>O</sub>3</sub>

    Note that Composition now has a to_html_string() method that may
    be used instead.

    :param formula:
    :return:
    """
    return re.sub(r"([A-Za-z\(\)])([\d\.]+)", r"\1<sub>\2</sub>", formula)


def unicodeify(formula):
    """
    Generates a formula with unicode subscripts, e.g. Fe2O3 is transformed
    to Fe₂O₃. Does not support formulae with decimal points.

    Note that Composition now has a to_unicode_string() method that may
    be used instead.

    :param formula:
    :return:
    """

    if "." in formula:
        raise ValueError("No unicode character exists for subscript period.")

    for original_subscript, subscript_unicode in SUBSCRIPT_UNICODE.items():
        formula = formula.replace(str(original_subscript), subscript_unicode)

    return formula


def latexify_spacegroup(spacegroup_symbol):
    r"""
    Generates a latex formatted spacegroup. E.g., P2_1/c is converted to
    P2$_{1}$/c and P-1 is converted to P$\\overline{1}$.

    Note that SymmetryGroup now has a to_latex_string() method that may
    be called instead.

    Args:
        spacegroup_symbol (str): A spacegroup symbol

    Returns:
        A latex formatted spacegroup with proper subscripts and overlines.
    """
    sym = re.sub(r"_(\d+)", r"$_{\1}$", spacegroup_symbol)
    return re.sub(r"-(\d)", r"$\\overline{\1}$", sym)


def unicodeify_spacegroup(spacegroup_symbol):
    r"""
    Generates a unicode formatted spacegroup. E.g., P2$_{1}$/c is converted to
    P2₁/c and P$\\overline{1}$ is converted to P̅1.

    Note that SymmetryGroup now has a to_unicode_string() method that
    may be called instead.

    Args:
        spacegroup_symbol (str): A spacegroup symbol as LaTeX

    Returns:
        A unicode spacegroup with proper subscripts and overlines.
    """

    if not spacegroup_symbol:
        return ""

    symbol = latexify_spacegroup(spacegroup_symbol)

    for number, unicode_number in SUBSCRIPT_UNICODE.items():
        symbol = symbol.replace("$_{" + str(number) + "}$", unicode_number)
        symbol = symbol.replace("_" + str(number), unicode_number)

    overline = "\u0305"  # u"\u0304" (macron) is also an option

    symbol = symbol.replace("$\\overline{", "")
    symbol = symbol.replace("$", "")
    symbol = symbol.replace("{", "")
    # overline unicode symbol comes after the character with the overline
    symbol = symbol.replace("}", overline)

    return symbol


def unicodeify_species(specie_string):
    """
    Generates a unicode formatted species string, with appropriate
    superscripts for oxidation states.

    Note that Species now has a to_unicode_string() method that
    may be used instead.

    Args:
        specie_string (str): Species string, e.g. O2-

    Returns:
        Species string, e.g. O²⁻
    """

    if not specie_string:
        return ""

    for character, unicode_character in SUPERSCRIPT_UNICODE.items():
        specie_string = specie_string.replace(character, unicode_character)

    return specie_string


def stream_has_colours(stream):
    """
    True if stream supports colours. Python cookbook, #475186
    """
    if not hasattr(stream, "isatty"):
        return False

    if not stream.isatty():
        return False  # auto color only on TTYs
    try:
        import curses

        curses.setupterm()
        return curses.tigetnum("colors") > 2
    except Exception:
        return False  # guess false in case of error


def transformation_to_string(matrix, translation_vec=(0, 0, 0), components=("x", "y", "z"), c="", delim=","):
    """
    Convenience method. Given matrix returns string, e.g. x+2y+1/4
    :param matrix
    :param translation_vec
    :param components: either ('x', 'y', 'z') or ('a', 'b', 'c')
    :param c: optional additional character to print (used for magmoms)
    :param delim: delimiter
    :return: xyz string
    """
    parts = []
    for i in range(3):
        s = ""
        m = matrix[i]
        t = translation_vec[i]
        for j, dim in enumerate(components):
            if m[j] != 0:
                f = Fraction(m[j]).limit_denominator()
                if s != "" and f >= 0:
                    s += "+"
                if abs(f.numerator) != 1:
                    s += str(f.numerator)
                elif f < 0:
                    s += "-"
                s += c + dim
                if f.denominator != 1:
                    s += "/" + str(f.denominator)
        if t != 0:
            s += ("+" if (t > 0 and s != "") else "") + str(Fraction(t).limit_denominator())
        if s == "":
            s += "0"
        parts.append(s)
    return delim.join(parts)


def disordered_formula(disordered_struct, symbols=("x", "y", "z"), fmt="plain"):
    """
    Returns a formula of a form like AxB1-x (x=0.5)
    for disordered structures. Will only return a
    formula for disordered structures with one
    kind of disordered site at present.

    Args:
        disordered_struct: a disordered structure
        symbols: a tuple of characters to use for
        subscripts, by default this is ('x', 'y', 'z')
        but if you have more than three disordered
        species more symbols will need to be added
        fmt (str): 'plain', 'HTML' or 'LaTeX'

    Returns (str): a disordered formula string
    """

    # this is in string utils and not in
    # Composition because we need to have access
    # to site occupancies to calculate this, so
    # have to pass the full structure as an argument
    # (alternatively this could be made a method on
    # Structure)
    from pymatgen.core.composition import Composition
    from pymatgen.core.periodic_table import get_el_sp

    if disordered_struct.is_ordered:
        raise ValueError("Structure is not disordered, so disordered formula not defined.")

    disordered_site_compositions = {site.species for site in disordered_struct if not site.is_ordered}

    if len(disordered_site_compositions) > 1:
        # this probably won't happen too often
        raise ValueError(
            "Ambiguous how to define disordered formula when more than one type of disordered site is present."
        )
    disordered_site_composition = disordered_site_compositions.pop()

    disordered_species = {str(sp) for sp, occu in disordered_site_composition.items()}

    if len(disordered_species) > len(symbols):
        # this probably won't happen too often either
        raise ValueError(f"Not enough symbols to describe disordered composition: {symbols}")
    symbols = list(symbols)[0 : len(disordered_species) - 1]

    comp = disordered_struct.composition.get_el_amt_dict().items()
    # sort by electronegativity, as per composition
    comp = sorted(comp, key=lambda x: get_el_sp(x[0]).X)

    disordered_comp = []
    variable_map = {}

    total_disordered_occu = sum(occu for sp, occu in comp if str(sp) in disordered_species)

    # composition to get common factor
    factor_comp = disordered_struct.composition.as_dict()
    factor_comp["X"] = total_disordered_occu
    for sp in disordered_species:
        del factor_comp[str(sp)]
    factor_comp = Composition.from_dict(factor_comp)
    factor = factor_comp.get_reduced_formula_and_factor()[1]

    total_disordered_occu /= factor
    remainder = formula_double_format(total_disordered_occu, ignore_ones=False) + "-" + "-".join(symbols)

    for sp, occu in comp:
        sp = str(sp)
        if sp not in disordered_species:
            disordered_comp.append((sp, formula_double_format(occu / factor)))
        else:
            if len(symbols) > 0:
                symbol = symbols.pop(0)
                disordered_comp.append((sp, symbol))
                variable_map[symbol] = occu / total_disordered_occu / factor
            else:
                disordered_comp.append((sp, remainder))

    if fmt == "LaTeX":
        sub_start = "_{"
        sub_end = "}"
    elif fmt == "HTML":
        sub_start = "<sub>"
        sub_end = "</sub>"
    elif fmt != "plain":
        raise ValueError("Unsupported output format, choose from: LaTeX, HTML, plain")

    disordered_formula = []
    for sp, occu in disordered_comp:
        disordered_formula.append(sp)
        if occu:  # can be empty string if 1
            if fmt != "plain":
                disordered_formula.append(sub_start)
            disordered_formula.append(occu)
            if fmt != "plain":
                disordered_formula.append(sub_end)
    disordered_formula.append(" ")
    disordered_formula += [f"{k}={formula_double_format(v)} " for k, v in variable_map.items()]

    return "".join(map(str, disordered_formula))[0:-1]
