# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides utility classes for string operations.
"""
from __future__ import unicode_literals
import re
from fractions import Fraction

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "$Sep 23, 2011M$"


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
    return returnstr + "\n".join([delimiter.join([str(m) for m in result])
                                  for result in results])


def formula_double_format(afloat, ignore_ones=True, tol=1e-8):
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
    elif abs(afloat - int(afloat)) < tol:
        return str(int(afloat))
    else:
        return str(round(afloat, 8))


def latexify(formula):
    """
    Generates a latex formatted formula. E.g., Fe2O3 is transformed to
    Fe$_{2}$O$_{3}$.

    Args:
        formula (str): Input formula.

    Returns:
        Formula suitable for display as in LaTeX with proper subscripts.
    """
    return re.sub(r"([A-Za-z\(\)])([\d\.]+)", r"\1$_{\2}$", formula)


def latexify_spacegroup(spacegroup_symbol):
    """
    Generates a latex formatted spacegroup. E.g., P2_1/c is converted to
    P2$_{1}$/c and P-1 is converted to P$\\overline{1}$.

    Args:
        spacegroup_symbol (str): A spacegroup symbol

    Returns:
        A latex formatted spacegroup with proper subscripts and overlines.
    """
    sym = re.sub(r"_(\d+)", r"$_{\1}$", spacegroup_symbol)
    return re.sub(r"-(\d)", r"$\\overline{\1}$", sym)


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
    except:
        return False  # guess false in case of error


def transformation_to_string(matrix, translation_vec=(0, 0, 0), components=('x', 'y', 'z'), c='', delim=','):
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
        s = ''
        m = matrix[i]
        t = translation_vec[i]
        for j, dim in enumerate(components):
            if m[j] != 0:
                f = Fraction(m[j]).limit_denominator()
                if s != '' and f >= 0:
                    s += '+'
                if abs(f.numerator) != 1:
                    s += str(f.numerator)
                elif f < 0:
                    s += '-'
                s += c + dim
                if f.denominator != 1:
                    s += '/' + str(f.denominator)
        if t != 0:
            s += ('+' if (t > 0 and s != '') else '') + str(Fraction(t).limit_denominator())
        if s == '':
            s += '0'
        parts.append(s)
    return delim.join(parts)


class StringColorizer(object):
    colours = {"default": "",
               "blue": "\x1b[01;34m",
               "cyan": "\x1b[01;36m",
               "green": "\x1b[01;32m",
               "red": "\x1b[01;31m",
               # lighting colours.
               #"lred":    "\x1b[01;05;37;41m"
               }

    def __init__(self, stream):
        self.has_colours = stream_has_colours(stream)

    def __call__(self, string, colour):
        if self.has_colours:
            code = self.colours.get(colour.lower(), "")
            if code:
                return code + string + "\x1b[00m"
            else:
                return string
        else:
            return string


if __name__ == "__main__":
    import doctest
    doctest.testmod()
