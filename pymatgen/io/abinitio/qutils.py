# coding: utf-8
"""
Collection of low-level tools that faciliate the interface with resource managers.

The preferred way of importing this module is:

    import qutils as qu
"""
from __future__ import print_function, division, unicode_literals

from monty.string import is_string
from pymatgen.core.units import Time, Memory

import logging
logger = logging.getLogger(__name__)


def slurm_parse_timestr(s):
    """
    A slurm time parser. Accepts a string in one the following forms:

        # "days-hours",
        # "days-hours:minutes",
        # "days-hours:minutes:seconds".
        # "minutes",
        # "minutes:seconds",
        # "hours:minutes:seconds",

    Returns:
        Time in seconds.

    Raises:
        `ValueError` if string is not valid.
    """
    days, hours, minutes, seconds = 0, 0, 0, 0
    if '-' in s:
        # "days-hours",
        # "days-hours:minutes",                                        
        # "days-hours:minutes:seconds".                                
        days, s = s.split("-")
        days = int(days)

        if ':' not in s:
            hours = int(float(s))
        elif s.count(':') == 1:
            hours, minutes = map(int, s.split(':'))
        elif s.count(':') == 2:
            hours, minutes, seconds = map(int, s.split(':'))
        else:
            raise ValueError("More that 2 ':' in string!")

    else:
        # "minutes",
        # "minutes:seconds",
        # "hours:minutes:seconds",
        if ':' not in s:
            minutes = int(float(s))
        elif s.count(':') == 1:
            minutes, seconds = map(int, s.split(':'))
        elif s.count(':') == 2:
            hours, minutes, seconds = map(int, s.split(':'))
        else:
            raise ValueError("More than 2 ':' in string!")

    return Time((days*24 + hours)*3600 + minutes*60 + seconds, "s")


def time2slurm(timeval, unit="s"):
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the slurm convention: "days-hours:minutes:seconds".

    >>> assert time2slurm(61) == '0-0:1:1' and time2slurm(60*60+1) == '0-1:0:1'
    >>> assert time2slurm(0.5, unit="h") == '0-0:30:0'
    """
    d, h, m, s = 24*3600, 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    days, hours = divmod(timeval, d)
    hours, minutes = divmod(hours, h)
    minutes, secs = divmod(minutes, m)

    return "%d-%d:%d:%d" % (days, hours, minutes, secs)


def time2pbspro(timeval, unit="s"):
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the PbsPro convention: "hours:minutes:seconds".

    >>> assert time2pbspro(2, unit="d") == '48:0:0' 
    """
    h, m, s = 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    hours, minutes = divmod(timeval, h)
    minutes, secs = divmod(minutes, m)

    return "%d:%d:%d" % (hours, minutes, secs)


def timelimit_parser(s):
    """Convert a float or a string into time in seconds."""
    try:
        return Time(float(s), "s")
    except ValueError:
        return slurm_parse_timestr(s)


def any2mb(s):
    """Convert string or number to memory in megabytes."""
    if is_string(s):
        return int(Memory.from_string(s).to("Mb"))
    else:
        return int(s)
