# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

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
    def __init__(self, cls, method, msg):
        self.cls = cls
        self.method = method
        self.msg = msg

    def __str__(self):
        return str(self.cls) + ': ' + self.method + '\n' + repr(self.msg)


class ClassTypeChemenvError(AbstractChemenvError):
    def __init__(self, object, classtype):
        self.object = object
        self.classtype = classtype

    def __str__(self):
        return '"{}" expected, {} found\n'.format(self.classtype.__name__, self.object.__class__.__name__)


class NeighborsNotComputedChemenvError(AbstractChemenvError):
    def __init__(self, site):
        self.site = site

    def __str__(self):
        return 'The neighbors were not computed for the following site : \n' + str(self.site)


class BVAValencesNotFoundChemenvError(AbstractChemenvError):
    def __init__(self, structure):
        self.structure = structure

    def __str__(self):
        return 'The valences were not found for the following structure : \n' +\
               self.structure.composition.reduced_formula


class ChemenvStrategyError(AbstractChemenvError):
    def __init__(self, cls, method, msg):
        self.cls = cls
        self.method = method
        self.msg = msg

    def __str__(self):
        return str(self.cls) + ': ' + self.method + '\n' + repr(self.msg)


class InitializationChemenvError(AbstractChemenvError):
    def __init__(self, cls):
        self.cls = cls

    def __str__(self):
        return 'There is some missing arguments for the initialization of a {} object'.format(self.cls)


class EquivalentSiteSearchError(AbstractChemenvError):
    def __init__(self, site):
        self.site = site

    def __str__(self):
        return 'Equivalent site could not be found for the following site : {}'.format(str(self.site))


class VoronoiParametersError(AbstractChemenvError):
    def __init__(self, vp):
        self.vp = vp

    def __str__(self):
        return 'The list of Voronoi parameters does not contain the following set of parameters :\n' \
               ' - distfactor : {},\n' \
               ' - angfactor : {},\n' \
               ' - only_anion_cation_bond : {}'.format(self.vp.distance_parameter, self.vp.angle_parameter,
                                                       self.vp.only_anion_cation_bond)


class SolidAngleError(AbstractChemenvError):
    def __init__(self, cosinus):
        self.cosinus = cosinus

    def __str__(self):
        return 'Value of cosinus ({}) from which an angle should be retrieved' \
               'is not between -1.0 and 1.0'.format(self.cosinus)


class RatioFunctionError(AbstractChemenvError):
    def __init__(self, function):
        self.function = function

    def __str__(self):
        return 'Function "{}" is not allowed as a ratio function.'.format(self.function)


class PenaltyFunctionError(AbstractChemenvError):
    def __init__(self, function):
        self.function = function

    def __str__(self):
        return 'Function "{}" is not possible.'.format(self.function)


class ChemenvError(Exception):
    def __init__(self, cls, method, msg):
        self.cls = cls
        self.method = method
        self.msg = msg

    def __str__(self):
        return str(self.cls) + ': ' + self.method + '\n' + repr(self.msg)