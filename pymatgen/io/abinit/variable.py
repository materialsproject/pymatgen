"""Support for Abinit input variables."""
from __future__ import annotations

import collections
import collections.abc
import string

import numpy as np

__all__ = [
    "InputVariable",
]

_SPECIAL_DATASET_INDICES = (":", "+", "?")

_DATASET_INDICES = "".join(list(string.digits) + list(_SPECIAL_DATASET_INDICES))

_INTERNAL_DATASET_INDICES = ("__s", "__i", "__a")

_SPECIAL_CONVERSION = zip(_INTERNAL_DATASET_INDICES, _SPECIAL_DATASET_INDICES)

_UNITS = {
    "bohr": 1.0,
    "angstrom": 1.8897261328856432,
    "hartree": 1.0,
    "Ha": 1.0,
    "eV": 0.03674932539796232,
}


class InputVariable:
    """
    An Abinit input variable.
    """

    def __init__(self, name, value, units="", valperline=3):
        """
        Args:
            name: Name of the variable.
            value: Value of the variable.
            units: String specifying one of the units supported by Abinit. Default: atomic units.
            valperline: Number of items printed per line.
        """
        self._name = name
        self.value = value
        self._units = units

        # Maximum number of values per line.
        self.valperline = valperline
        if name in ["bdgw"]:
            self.valperline = 2

        if is_iter(self.value) and isinstance(self.value[-1], str) and self.value[-1] in _UNITS:
            self.value = list(self.value)
            self._units = self.value.pop(-1)

    def get_value(self):
        """Return the value."""
        if self.units:
            return list(self.value) + [self.units]
        return self.value

    @property
    def name(self):
        """Name of the variable."""
        return self._name

    @property
    def basename(self):
        """Return the name trimmed of any dataset index."""
        basename = self.name
        return basename.rstrip(_DATASET_INDICES)

    @property
    def dataset(self):
        """Return the dataset index in string form."""
        return self.name.split(self.basename)[-1]

    @property
    def units(self):
        """Return the units."""
        return self._units

    def __str__(self):
        """Declaration of the variable in the input file."""
        value = self.value
        if value is None or not str(value):
            return ""

        var = self.name
        line = " " + var

        # By default, do not impose a number of decimal points
        floatdecimal = 0

        # For some inputs, enforce number of decimal points...
        if any(inp in var for inp in ("xred", "xcart", "rprim", "qpt", "kpt")):
            floatdecimal = 16

        # ...but not for those
        if any(inp in var for inp in ("ngkpt", "kptrlatt", "ngqpt", "ng2qpt")):
            floatdecimal = 0

        if isinstance(value, np.ndarray):
            n = 1
            for i in np.shape(value):
                n *= i
            value = np.reshape(value, n)
            value = list(value)

        # values in lists
        if isinstance(value, (list, tuple)):

            # Reshape a list of lists into a single list
            if all(isinstance(v, (list, tuple)) for v in value):
                line += self.format_list2d(value, floatdecimal)

            else:
                line += self.format_list(value, floatdecimal)

        # scalar values
        else:
            line += " " + str(value)

        # Add units
        if self.units:
            line += " " + self.units

        return line

    @staticmethod
    def format_scalar(val, floatdecimal=0):
        """
        Format a single numerical value into a string
        with the appropriate number of decimal.
        """
        sval = str(val)
        if sval.lstrip("-").lstrip("+").isdigit() and floatdecimal == 0:
            return sval

        try:
            fval = float(val)
        except Exception:
            return sval

        if fval == 0 or (1e-3 < abs(fval) < 1e4):
            form = "f"
            addlen = 5
        else:
            form = "e"
            addlen = 8

        ndec = max(len(str(fval - int(fval))) - 2, floatdecimal)
        ndec = min(ndec, 10)

        sval = f"{fval:>{ndec + addlen}.{ndec}{form}}"

        sval = sval.replace("e", "d")

        return sval

    @staticmethod
    def format_list2d(values, floatdecimal=0):
        """Format a list of lists."""
        lvals = flatten(values)

        # Determine the representation
        if all(isinstance(v, int) for v in lvals):
            type_all = int
        else:
            try:
                for v in lvals:
                    float(v)
                type_all = float
            except Exception:
                type_all = str

        # Determine the format
        width = max(len(str(s)) for s in lvals)
        if type_all == int:
            formatspec = f">{width}d"
        elif type_all == str:
            formatspec = f">{width}"
        else:

            # Number of decimal
            maxdec = max(len(str(f - int(f))) - 2 for f in lvals)
            ndec = min(max(maxdec, floatdecimal), 10)

            if all(f == 0 or (abs(f) > 1e-3 and abs(f) < 1e4) for f in lvals):
                formatspec = f">{ndec + 5}.{ndec}f"
            else:
                formatspec = f">{ndec + 8}.{ndec}e"

        line = "\n"
        for L in values:
            for val in L:
                line += f" {val:{{formatspec}}}"
            line += "\n"

        return line.rstrip("\n")

    def format_list(self, values, floatdecimal=0):
        """
        Format a list of values into a string.
        The result might be spread among several lines.
        """
        line = ""

        # Format the line declaring the value
        for i, val in enumerate(values):
            line += " " + self.format_scalar(val, floatdecimal)
            if self.valperline is not None and (i + 1) % self.valperline == 0:
                line += "\n"

        # Add a carriage return in case of several lines
        if "\n" in line.rstrip("\n"):
            line = "\n" + line

        return line.rstrip("\n")


def is_iter(obj) -> bool:
    """Return True if the argument is list-like."""
    return hasattr(obj, "__iter__")


def flatten(iterable):
    """Make an iterable flat, i.e. a 1d iterable object."""
    iterator = iter(iterable)
    array, stack = collections.deque(), collections.deque()
    while True:
        try:
            value = next(iterator)
        except StopIteration:
            if not stack:
                return tuple(array)
            iterator = stack.pop()
        else:
            if not isinstance(value, str) and isinstance(value, collections.abc.Iterable):
                stack.append(iterator)
                iterator = iter(value)
            else:
                array.append(value)
