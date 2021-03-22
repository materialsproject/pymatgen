# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module contains the error classes for the chemenv package.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


class AbstractChemenvError(Exception):
    """
    Abstract class for Chemenv errors.
    """

    def __init__(self, cls, method, msg):
        """
        :param cls:
        :param method:
        :param msg:
        """
        self.cls = cls
        self.method = method
        self.msg = msg

    def __str__(self):
        return str(self.cls) + ": " + self.method + "\n" + repr(self.msg)


class NeighborsNotComputedChemenvError(AbstractChemenvError):
    """
    Neighbors not computed error.
    """

    def __init__(self, site):
        """
        :param site:
        """
        self.site = site

    def __str__(self):
        return "The neighbors were not computed for the following site : \n" + str(self.site)


class EquivalentSiteSearchError(AbstractChemenvError):
    """
    Equivalent site search error.
    """

    def __init__(self, site):
        """
        :param site:
        """
        self.site = site

    def __str__(self):
        return "Equivalent site could not be found for the following site : {}".format(str(self.site))


class SolidAngleError(AbstractChemenvError):
    """
    Solid angle error.
    """

    def __init__(self, cosinus):
        """
        :param cosinus:
        """
        self.cosinus = cosinus

    def __str__(self):
        return "Value of cosinus ({}) from which an angle should be retrieved" "is not between -1.0 and 1.0".format(
            self.cosinus
        )


class ChemenvError(Exception):
    """
    Chemenv error.
    """

    def __init__(self, cls, method, msg):
        """
        :param cls:
        :param method:
        :param msg:
        """
        self.cls = cls
        self.method = method
        self.msg = msg

    def __str__(self):
        return str(self.cls) + ": " + self.method + "\n" + repr(self.msg)
