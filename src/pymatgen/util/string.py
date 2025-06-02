"""This module provides utility classes for string operations.

TODO: make standalone functions in this module use the same implementation as Stringify
Note: previous deprecations of standalone functions in this module were removed due to
a community need.
"""

from __future__ import annotations

import math
import re
from fractions import Fraction
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, Literal, TextIO

    from numpy.typing import ArrayLike

    from pymatgen.core import Structure

SUBSCRIPT_UNICODE: dict[str, str] = {
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

SUPERSCRIPT_UNICODE: dict[str, str] = {
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


class Stringify:
    """Mix-in class for string formatting, e.g. superscripting numbers and symbols or superscripting."""

    STRING_MODE = "SUBSCRIPT"

    def to_pretty_string(self) -> str:
        """A pretty string representation. By default, the __str__ output is used, but this method can be
        overridden if a different representation from default is desired.
        """
        return str(self)

    def to_latex_string(self) -> str:
        """Generate a LaTeX formatted string. The mode is set by the class variable STRING_MODE, which defaults to
        "SUBSCRIPT". e.g. Fe2O3 is transformed to Fe$_{2}$O$_{3}$. Setting STRING_MODE to "SUPERSCRIPT" creates
        superscript, e.g. Fe2+ becomes Fe^{2+}. The initial string is obtained from the class's __str__ method.

        Returns:
            str: for LaTeX display with proper sub-/superscripts.
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
        """Generate a HTML formatted string. This uses the output from to_latex_string to generate a HTML output.

        Returns:
            HTML formatted string.
        """
        str_ = re.sub(r"\$_\{([^}]+)\}\$", r"<sub>\1</sub>", self.to_latex_string())
        str_ = re.sub(r"\$\^\{([^}]+)\}\$", r"<sup>\1</sup>", str_)
        return re.sub(
            r"\$\\overline\{([^}]+)\}\$",
            r'<span style="text-decoration:overline">\1</span>',
            str_,
        )

    def to_unicode_string(self) -> str:
        """Unicode string with proper sub and superscripts. Note that this works only
        with systems where the sub and superscripts are pure integers.
        """
        str_ = self.to_latex_string()
        for match in re.finditer(r"\$_\{(\d+)\}\$", str_):
            s1 = match.group()
            s2 = [SUBSCRIPT_UNICODE[s] for s in match[1]]
            str_ = str_.replace(s1, "".join(s2))
        for match in re.finditer(r"\$\^\{([\d\+\-]+)\}\$", str_):
            s1 = match.group()
            s2 = [SUPERSCRIPT_UNICODE[s] for s in match[1]]
            str_ = str_.replace(s1, "".join(s2))
        return str_


def str_delimited(
    results: Sequence[Sequence[Any]],
    header: Sequence[str] | None = None,
    delimiter: str = "\t",
) -> str:
    r"""Given a tuple of tuples, generate a delimited string form.
    >>> results = (("a", "b", "c"), ("d", "e", "f"), (1, 2, 3))
    >>> print(str_delimited(results, delimiter=","))
        a,b,c
        d,e,f
        1,2,3.

    Args:
        results (Sequence[Sequence[Any]]): 2D sequence of arbitrary types.
        header (Sequence[str]): optional headers.
        delimiter (str): Defaults to "\t" for tab-delimited output.

    Returns:
        str: Aligned string output in a table-like format.
    """
    out = ""
    if header is not None:
        out += f"{delimiter.join(header)}\n"
    return out + "\n".join(delimiter.join([str(m) for m in result]) for result in results)


def formula_double_format(
    afloat: float,
    ignore_ones: bool = True,
    tol: float = 1e-8,
) -> float | Literal[""]:
    """Format a float for pretty formulas.
    E.g. "Li1.0 Fe1.0 P1.0 O4.0" -> "LiFePO4".

    Args:
        afloat (float): The float to be formatted.
        ignore_ones (bool): if true, floats of 1.0 are ignored.
        tol (float): Absolute tolerance to round to nearest int. i.e. (2 + 1E-9) -> 2.

    Returns:
        float | "": Formatted float for formulas.
    """
    if ignore_ones and math.isclose(afloat, 1, abs_tol=tol):
        return ""

    if math.isclose(afloat, round(afloat), abs_tol=tol, rel_tol=0):
        return round(afloat)
    return round(afloat, 8)


def charge_string(charge: float, brackets: bool = True, explicit_one: bool = True) -> str:
    """Get a string representing the charge of an Ion. By default, the
    charge is placed in brackets with the sign preceding the magnitude, e.g.
    '[+2]'. For uncharged species, the string returned is '(aq)'.

    Args:
        charge (float): The charge of the ion.
        brackets (bool): Whether to enclose the charge in brackets, e.g. [+2]. Default is True.
        explicit_one (bool): whether to include the number one for monovalent ions,
            e.g. "+1" rather than "+". Default is True.
    """
    chg_str = "(aq)" if charge == 0 else f"{formula_double_format(charge, ignore_ones=False):+}"

    if chg_str in ["+1", "-1"] and not explicit_one:
        chg_str = chg_str.replace("1", "")

    if chg_str != "(aq)" and brackets:
        chg_str = f"[{chg_str}]"

    return chg_str


def latexify(formula: str, bold: bool = False) -> str:
    """Generate a LaTeX formatted formula. e.g. Fe2O3 is transformed to
    Fe$_{2}$O$_{3}$.

    Note that Composition now has `to_latex_string` method that may
    be used instead.

    Args:
        formula (str): Input formula.
        bold (bool): Whether to make the subscripts bold. Defaults to False.

    Returns:
        str: Formula suitable for display as in LaTeX with proper subscripts.
    """
    return re.sub(
        r"([A-Za-z\(\)])([\d\.]+)",
        r"\1$_{\\mathbf{\2}}$" if bold else r"\1$_{\2}$",
        formula,
    )


def htmlify(formula: str) -> str:
    """Generate a HTML formatted formula, e.g. Fe2O3 is transformed to
    Fe<sub>2</sub>O</sub>3</sub>.

    Note that Composition now has a to_html_string() method that may
    be used instead.

    Args:
        formula (str): The string to format.
    """
    return re.sub(r"([A-Za-z\(\)])([\d\.]+)", r"\1<sub>\2</sub>", formula)


def unicodeify(formula: str) -> str:
    """Generate a formula with unicode subscripts, e.g. Fe2O3 is transformed
    to Fe₂O₃. Does not support formulae with decimal points.

    Note that Composition now has a to_unicode_string() method that may
    be used instead.

    Args:
        formula (str): The string to format.
    """
    if "." in formula:
        raise ValueError("No unicode character exists for subscript period.")

    for original_subscript, subscript_unicode in SUBSCRIPT_UNICODE.items():
        formula = formula.replace(str(original_subscript), subscript_unicode)

    return formula


def latexify_spacegroup(spacegroup_symbol: str) -> str:
    r"""Generate a latex formatted spacegroup. e.g. P2_1/c is converted to
    P2$_{1}$/c and P-1 is converted to P$\\overline{1}$.

    Note that SymmetryGroup now has a to_latex_string() method that may
    be called instead.

    Args:
        spacegroup_symbol (str): A spacegroup symbol

    Returns:
        str: A latex formatted spacegroup with proper subscripts and overlines.
    """
    sym = re.sub(r"_(\d+)", r"$_{\1}$", spacegroup_symbol)
    return re.sub(r"-(\d)", r"$\\overline{\1}$", sym)


def unicodeify_spacegroup(spacegroup_symbol: str) -> str:
    r"""Generate a unicode formatted spacegroup. e.g. P2$_{1}$/c is converted to
    P2₁/c and P$\\overline{1}$ is converted to P̅1.

    Note that SymmetryGroup now has a to_unicode_string() method that
    may be called instead.

    Args:
        spacegroup_symbol (str): A spacegroup symbol as LaTeX.

    Returns:
        str: A unicode spacegroup with proper subscripts and overlines.
    """
    if not spacegroup_symbol:
        return ""

    symbol = latexify_spacegroup(spacegroup_symbol)

    for num, unicode_number in SUBSCRIPT_UNICODE.items():
        symbol = symbol.replace(f"$_{{{num}}}$", unicode_number)
        symbol = symbol.replace(f"_{num}", unicode_number)

    overline = "\u0305"  # u"\u0304" (macron) is also an option

    symbol = symbol.replace("$\\overline{", "")
    symbol = symbol.replace("$", "")
    symbol = symbol.replace("{", "")
    # overline unicode symbol comes after the character with the overline
    return symbol.replace("}", overline)


def unicodeify_species(specie_string: str) -> str:
    """Generate a unicode formatted species string, with appropriate
    superscripts for oxidation states.

    Note that Species now has `to_unicode_string` method that
    may be used instead.

    Args:
        specie_string (str): Species string, e.g. "O2-"

    Returns:
        str: Species string, e.g. "O²⁻"
    """
    if not specie_string:
        return ""

    for char, unicode_char in SUPERSCRIPT_UNICODE.items():
        specie_string = specie_string.replace(char, unicode_char)

    return specie_string


def stream_has_colors(stream: TextIO) -> bool:
    """True if stream supports colors. Python cookbook, #475186."""
    if not hasattr(stream, "isatty"):
        return False

    if not stream.isatty():
        return False  # auto color only on TTYs
    try:
        import curses

        curses.setupterm()
    except Exception:
        return False  # guess false in case of error
    else:
        return curses.tigetnum("colors") > 2


def transformation_to_string(
    matrix: ArrayLike,
    translation_vec: ArrayLike = (0, 0, 0),
    components: tuple[str, str, str] = ("x", "y", "z"),
    c: str = "",
    delim: str = ",",
) -> str:
    """Convenience method. Given matrix returns string, e.g. x+2y+1/4.

    Args:
        matrix (ArrayLike): A 3x3 matrix.
        translation_vec (tuple[float, float, float]): The translation vector. Defaults to (0, 0, 0).
        components(tuple[str, str, str]): The components. Either ('x', 'y', 'z') or ('a', 'b', 'c').
            Defaults to ('x', 'y', 'z').
        c (str): An optional additional character to print (used for magmoms). Defaults to "".
        delim (str): A delimiter. Defaults to ",".

    Returns:
        str: xyz string.
    """
    parts = []
    matrix = np.asarray(matrix)
    translation_vec = np.asarray(translation_vec)
    for idx in range(3):
        string = ""
        mat = matrix[idx]
        offset = translation_vec[idx]
        for j, dim in enumerate(components):
            if mat[j] != 0:
                f = Fraction(mat[j]).limit_denominator()
                if string != "" and f >= 0:
                    string += "+"
                if abs(f.numerator) != 1:
                    string += str(f.numerator)
                elif f < 0:
                    string += "-"
                string += c + dim
                if f.denominator != 1:
                    string += f"/{f.denominator}"
        if offset != 0:
            string += ("+" if (offset > 0 and string != "") else "") + str(Fraction(offset).limit_denominator())
        if string == "":
            string += "0"
        parts.append(string)
    return delim.join(parts)


def disordered_formula(
    disordered_struct: Structure,
    symbols: Sequence[str] = ("x", "y", "z"),
    fmt: Literal["plain", "HTML", "LaTex"] = "plain",
) -> str:
    """Get a formula of a form like AxB1-x (x=0.5)
    for disordered structures. Will only return a
    formula for disordered structures with one
    kind of disordered site at present.

    Args:
        disordered_struct (Structure): a disordered structure.
        symbols (Sequence[str]): Characters to use for subscripts,
            by default this is ('x', 'y', 'z') but if you have more than three
            disordered species more symbols will need to be added.
        fmt (str): 'plain', 'HTML' or 'LaTeX',

    Returns:
        str: a disordered formula string
    """
    # this is in string utils and not in Composition because we need to have access to
    # site occupancies to calculate this, so have to pass the full structure as an
    # argument (alternatively this could be made a method on Structure)
    from pymatgen.core import Composition, get_el_sp

    if disordered_struct.is_ordered:
        raise ValueError("Structure is not disordered, so disordered formula not defined.")

    disordered_site_compositions = {site.species for site in disordered_struct if not site.is_ordered}

    if len(disordered_site_compositions) > 1:
        # this probably won't happen too often
        raise ValueError(
            "Ambiguous how to define disordered formula when more than one type of disordered site is present."
        )
    disordered_site_composition = disordered_site_compositions.pop()

    disordered_species = {str(sp) for sp, _occu in disordered_site_composition.items()}

    if len(disordered_species) > len(symbols):
        # this probably won't happen too often either
        raise ValueError(f"Not enough symbols to describe disordered composition: {symbols}")
    symbols = list(symbols)[: len(disordered_species) - 1]

    comp = disordered_struct.composition.get_el_amt_dict().items()
    # sort by electronegativity, as per composition
    comp = sorted(comp, key=lambda x: get_el_sp(x[0]).X)  # type:ignore[assignment]

    disordered_comp: list[tuple[str, str]] = []
    variable_map = {}

    total_disordered_occu = sum(occu for sp, occu in comp if str(sp) in disordered_species)

    # composition to get common factor
    factor_comp = disordered_struct.composition.as_dict()
    factor_comp["X"] = total_disordered_occu
    for sp in disordered_species:
        del factor_comp[str(sp)]
    factor_comp = Composition.from_dict(factor_comp)  # type:ignore[assignment]
    factor = factor_comp.get_reduced_formula_and_factor()[1]

    total_disordered_occu /= factor
    remainder = f"{formula_double_format(total_disordered_occu, ignore_ones=False)}-{'-'.join(symbols)}"

    for sp, occu in comp:
        species = str(sp)
        if species not in disordered_species:
            disordered_comp.append((species, str(formula_double_format(occu / factor))))
        elif len(symbols) > 0:
            symbol = symbols.pop(0)
            disordered_comp.append((species, symbol))
            variable_map[symbol] = occu / total_disordered_occu / factor
        else:
            disordered_comp.append((species, remainder))

    if fmt == "LaTeX":
        sub_start = "_{"
        sub_end = "}"
    elif fmt == "HTML":
        sub_start = "<sub>"
        sub_end = "</sub>"
    elif fmt == "plain":
        sub_start = ""
        sub_end = ""
    else:
        raise ValueError("Unsupported output format, choose from: LaTeX, HTML, plain")

    disordered_formulas = []
    for sp, occu in disordered_comp:  # type:ignore[assignment]
        disordered_formulas.append(sp)
        if occu:  # can be empty string if 1
            if fmt != "plain":
                disordered_formulas.append(sub_start)
            disordered_formulas.append(occu)  # type:ignore[arg-type]
            if fmt != "plain":
                disordered_formulas.append(sub_end)
    disordered_formulas.append(" ")
    disordered_formulas += [f"{key}={formula_double_format(val)} " for key, val in variable_map.items()]

    return "".join(map(str, disordered_formulas))[:-1]
