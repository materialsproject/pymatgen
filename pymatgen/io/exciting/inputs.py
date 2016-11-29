# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import os
import re
import itertools
import warnings
import logging

import six
import numpy as np
from numpy.linalg import det
from collections import OrderedDict, namedtuple
from hashlib import md5

from monty.io import zopen
from monty.os.path import zpath
from monty.json import MontyDecoder

import xml.etree.cElementTree as ET
from enum import Enum
from tabulate import tabulate

import scipy.constants as const

from pymatgen import SETTINGS
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element, get_el_sp
from monty.design_patterns import cached_class
from pymatgen.util.string_utils import str_delimited
from pymatgen.util.io_utils import clean_lines
from monty.json import MSONable

"""
Classes for reading/manipulating/writing exciting input files.
"""

__author__ = "Christian Vorwerk"
__copyright__ = "Copyright 2016"
__version__ = "1.0"
__maintainer__ = "Christian Vorwerk"
__email__ = "vorwerk@physik.hu-berlin.de"
__status__ = "Development"
__date__ = "Nov 28, 2016"

class input(MSONable):
    """
        Object for representing the data stored in the structure part of the
        exciting input.

        Args:
            structure (Structure):  Structure object.
            title (str): Optional title for exciting input. Defaults to unit
                cell formula of structure. Defaults to None.
            lockxyz (Nx3 array): bool values for selective dynamics,
                where N is number of sites. Defaults to None.

        .. attribute:: structure

            Associated Structure.

        .. attribute:: title

            Optional title string.

        .. attribute:: lockxyz

            Lockxyz attribute for each site if available. A Nx3 array of
            booleans.
    """
    def __init__(self, structure, title=None, lockxyz=None):
        if structure.is_ordered:
            site_properties = {}
            if lockxyz:
                site_properties["selective_dynamics"] = lockxyz
            self.structure = structure.copy(site_properties=site_properties)
            self.title = structure.formula if title is None else title
        else:
            raise ValueError("Structure with partial occupancies cannot be "
                             "converted into exciting input!")
    @property
    def lockxyz(self):
        return self.structure.site_properties.get("selective_dynamics")
    @lockxyz.setter
    def lockxyz(self, lockxyz):
        self.structure.add_site_property("selective_dynamics",
                                         lockxyz)
    @staticmethod
    def from_string(data):
        """
        Reads the exciting input from a string
        """
        root=ET.fromstring(data)
        speciesnode=root.find('structure').iter('species')
        elements = []
        positions = []
        vectors=[]
        lockxyz=[]
        # get title
        title_in=str(root.find('title').text)
        # Read elements and coordinates
        for nodes in speciesnode:
            symbol = nodes.get('speciesfile').split('.')[0]
            if len(symbol.split('_'))==2:
              symbol=symbol.split('_')[0]
            if Element.is_valid_symbol(symbol):
                # Try to recognize the element symbol
                element = symbol
            else:
                raise NLValueError("Unknown element!")
            natoms = nodes.getiterator('atom')
            for atom in natoms:
                x, y, z = atom.get('coord').split()
                positions.append([float(x), float(y), float(z)])
                elements.append(element)
                # Obtain lockxyz for each atom
                if atom.get('lockxyz') is not None:
                    lxy=[]
                    for l in atom.get('lockxyz').split():
                        if l=='True' or l=='true':
                            lxyz.append(True)
                        else:
                            lxyz.append(False)
                    lockxyz.append(lxyz)
                else:
                    lockxyz.append([False, False, False])
        #check the atomic positions type
        if 'cartesian' in root.find('structure').attrib.keys():
          if root.find('structure').attrib['cartesian']:
            cartesian=True
        else:
          cartesian=False
        # get the scale attribute
        scale_in=root.find('structure').find('crystal').get('scale')
        if scale_in:
            scale=float(scale_in)
        else:
            scale=1.0
        # define conversion factor between Bohr radius and Angstrom
        bohr2ang=const.value('Bohr radius')/const.value('Angstrom star')
        # get the stretch attribute
        stretch_in=root.find('structure').find('crystal').get('stretch')
        if stretch_in:
          stretch=np.array([float(a) for a in stretch_in])
        else:
          stretch=np.array([1.0,1.0,1.0])
        # get basis vectors and scale them accordingly
        basisnode=root.find('structure').find('crystal').iter('basevect')
        for vect in basisnode:
          x, y, z=vect.text.split()
          vectors.append([float(x)*stretch[0],
                          float(y)*stretch[1],
                          float(z)*stretch[2]])
        # create lattice and structure object
        lattice_in=Lattice(vectors)
        structure_in=Structure(lattice_in,elements,positions,
                               False,cartesian,False)

        return input(structure_in, title_in, lockxyz)
    @staticmethod
    def from_file(filename):
        with zopen(filename, 'rt') as f:
            data=f.read().replace('\n','')
        return input.from_string(data)
