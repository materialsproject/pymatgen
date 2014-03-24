#!/usr/bin/env python

"""
This module provides utility classes for string operations.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "$Sep 23, 2011M$"

import re
import sys
import fnmatch


def generate_latex_table(results, header=None, caption=None, label=None):
    """
    Generates a string latex table from a sequence of sequence.

    Args:
        result: 2d sequence of arbitrary types.
        header: optional header

    Returns:
        String representation of Latex table with data.
    """
    body = []
    if header is not None:
        body.append(" & ".join(header) + "\\\\")
        body.append("\\hline")
    maxlength = 0
    for result in results:
        maxlength = max(maxlength, len(result))
        body.append(" & ".join([str(m) for m in result]) + "\\\\")
    colstr = "c" * maxlength
    output = ["\\begin{table}[H]",
              "\\caption{{{}}}".format(caption if caption else
              "Caption"), "\\label{{{}}}".format(label if label else "Label"),
              "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}" +
              colstr + "}", "\\hline", "\n".join(body), "\\hline",
              "\\end{tabular*}", "\\end{table}"]
    return "\n".join(output)


def str_delimited(results, header=None, delimiter="\t"):
    """
    Given a tuple of tuples, generate a delimited string form.
    >>> results = [["a","b","c"],["d","e","f"],[1,2,3]]
    >>> print str_delimited(results,delimiter=",")
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


def str_aligned(results, header=None):
    """
    Given a tuple, generate a nicely aligned string form.
    >>> results = [["a","b","cz"],["d","ez","f"],[1,2,3]]
    >>> print str_aligned(results)
    a    b   cz
    d   ez    f
    1    2    3

    Args:
        result: 2d sequence of arbitrary types.
        header: optional header

    Returns:
        Aligned string output in a table-like format.
    """
    k = list(zip(*results))
    stringlengths = list()
    count = 0
    for i in k:
        col_max_len = max([len(str(m)) for m in i])
        if header is not None:
            col_max_len = max([len(str(header[count])), col_max_len])
        stringlengths.append(col_max_len)
        count += 1
    format_string = "   ".join(["%" + str(d) + "s" for d in stringlengths])
    returnstr = ""
    if header is not None:
        header_str = format_string % tuple(header)
        returnstr += header_str + "\n"
        returnstr += "-" * len(header_str) + "\n"
    return returnstr + "\n".join([format_string % tuple(result)
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
        return str(afloat)


def latexify(formula):
    """
    Generates a latex formatted formula. E.g., Fe2O3 is transformed to
    Fe$_{2}$O$_{3}$.

    Args:
        formula (str): Input formula.

    Returns:
        Formula suitable for display as in LaTeX with proper subscripts.
    """
    return re.sub(r"([A-Za-z\(\)])(\d+)", r"\1$_{\2}$", formula)


def latexify_spacegroup(spacegroup_symbol):
    """
    Generates a latex formatted spacegroup. E.g., P2_1/c is converted to
    P2$_{1}$/c and P-1 is converted to P$\overline{1}$.

    Args:
        spacegroup_symbol (str): A spacegroup symbol

    Returns:
        A latex formatted spacegroup with proper subscripts and overlines.
    """
    sym = re.sub(r"_(\d+)", r"$_{\1}$", spacegroup_symbol)
    return re.sub(r"-(\d)", r"$\overline{\1}$", sym)


def pprint_table(table, out=sys.stdout, rstrip=False):
    """
    Prints out a table of data, padded for alignment
    Each row must have the same number of columns.

    Args:
        table: The table to print. A list of lists.
        out: Output stream (file-like object)
        rstrip: if True, trailing withespaces are removed from the entries.
    """
    def max_width_col(table, col_idx):
        """
        Get the maximum width of the given column index
        """
        return max([len(row[col_idx]) for row in table])

    if rstrip:
        for row_idx, row in enumerate(table):
            table[row_idx] = [c.rstrip() for c in row]

    col_paddings = []
    ncols = len(table[0])
    for i in range(ncols):
        col_paddings.append(max_width_col(table, i))

    for row in table:
        # left col
        out.write(row[0].ljust(col_paddings[0] + 1))
        # rest of the cols
        for i in range(1, len(row)):
            col = row[i].rjust(col_paddings[i] + 2)
            out.write(col)
        out.write("\n")


def is_string(s):
    """True if s behaves like a string (duck typing test)."""
    try:
        s + " "
        return True

    except TypeError:
        return False


def list_strings(arg):
    """
    Always return a list of strings, given a string or list of strings as
    input.

    :Examples:

    >>> list_strings('A single string')
    ['A single string']

    >>> list_strings(['A single string in a list'])
    ['A single string in a list']

    >>> list_strings(['A','list','of','strings'])
    ['A', 'list', 'of', 'strings']
    """
    if is_string(arg):
        return [arg]
    else:
        return arg


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


class WildCard(object):
    """
    This object provides an easy-to-use interface for
    filename matching with shell patterns (fnmatch).

    .. example:

    >>> w = WildCard("*.nc|*.pdf")
    >>> w.filter(["foo.nc", "bar.pdf", "hello.txt"])
    ['foo.nc', 'bar.pdf']

    >>> w.filter("foo.nc")
    ['foo.nc']
    """
    def __init__(self, wildcard, sep="|"):
        """
        Initializes a WildCard.

        Args:
            wildcard (str): String of tokens separated by sep. Each token
                represents a pattern.
            sep (str): Separator for shell patterns.
        """
        self.pats = ["*"]
        if wildcard:
            self.pats = wildcard.split(sep)

    def __str__(self):
        return "<%s, patterns = %s>" % (self.__class__.__name__, self.pats)

    def filter(self, names):
        """
        Returns a list with the names matching the pattern.
        """
        names = list_strings(names)

        fnames = []
        for f in names:
            for pat in self.pats:
                if fnmatch.fnmatch(f, pat):
                    fnames.append(f)

        return fnames

    def match(self, name):
        """
        Returns True if name matches one of the patterns.
        """
        for pat in self.pats:
            if fnmatch.fnmatch(name, pat):
                return True

        return False


if __name__ == "__main__":
    import doctest
    doctest.testmod()
