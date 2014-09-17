# coding: utf-8

from __future__ import division, unicode_literals

"""
This module defines some useful design patterns.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"


class Enum(set):
    """
    Creates an enum out of a set.
    """
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError


class NullFile(object):
    """A file object that is associated to /dev/null."""
    def __new__(cls):
        import os
        return open(os.devnull, 'w')

    def __init__(self):
        """no-op"""


class NullStream(object):
    """A fake stream with a no-op write.."""
    def write(*args):
        """no-op"""

