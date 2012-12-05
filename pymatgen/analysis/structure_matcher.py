#!/usr/bin/env python

"""
This module provides classes to perform fitting of two structures.
"""

from __future__ import division

__author__ = "Stephen Dacek, William Davidson Richards, Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Stephen Dacek"
__email__ = "sdacek@mit.edu"
__status__ = "Beta"
__date__ = "Dec 3, 2012"

from pymatgen.core.structure_modifier import BasisChange
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.core.lattice import Lattice
from pymatgen.transformations.standard_transformations import \
    PrimitiveCellTransformation
from pymatgen.util.coord_utils import find_in_coord_list_pbc

import numpy as np
import itertools


def basis_change(coords, new_lattice):
    #Convert cartesian coordinates to new basis
    frac_coords = []
    for i in range(len(coords)):
        frac_coords.append(new_lattice.get_fractional_coords(coords[i]))

    return frac_coords


def cmp_struct(s1, s2, frac_tol):
    #compares the fractional coordinates
    for s1_coords, s2_coords in zip(s1, s2):
        #Available vectors
        avail = [1] * len(s1_coords)
        for coord in s1_coords:
            ind = find_in_coord_list_pbc(s2_coords, coord, frac_tol)
            #if more than one match found, take closest
            if len(ind) > 1:
                fcoords = np.tile(coord, (len(s2_coords[ind]), 1))
                fdist = np.mod(s2_coords[ind], 1) - np.mod(fcoords, 1)
                fdist -= np.round(fdist)
                avail[np.where(
                    np.all(np.abs(fdist) == np.min(np.abs(fdist)), axis=1))[
                      0]] = 0
            elif len(ind) and avail[ind]:
                avail[ind] = 0
            else:
                return False
    return True


def compare_structures(struct1, struct2, ltol=0.2, atol=5,
                       primitive_cell=False,
                       equivalent_species=[]):
    # Same number of sites
    if struct1.num_sites != struct2.num_sites:
        return False

    #primitive cell transformation - Needs work
    if primitive_cell:
        prim = PrimitiveCellTransformation()
        struct1 = prim.apply_transformation(struct1)
        struct2 = prim.apply_transformation(struct2)

    #compute niggli lattices, 
    nl1 = struct1.lattice.get_niggli_reduced_lattice()
    nl2 = struct2.lattice.get_niggli_reduced_lattice()
    struct1 = BasisChange(struct1, nl1).modified_structure
    struct2 = BasisChange(struct2, nl2).modified_structure

    #rescale lattice to same volume
    se = StructureEditor(struct2)
    nl2 = Lattice(nl2.matrix * (nl1.volume / nl2.volume) ** (1.0 / 3))
    se.modify_lattice(nl2)
    struct2 = se.modified_structure

    #fractional tolerance
    frac_tol = np.array([ltol / i for i in struct1.lattice.abc])

    #get possible new lattice vectors
    ds = Structure(struct2.lattice, ['X'], [[0, 0, 0]])
    nv = []
    for i in range(3):
        l = struct1.lattice.abc[i]
        nvi = []
        vs = ds.get_neighbors_in_shell([0, 0, 0], l, ltol * l)
        for site, dist in vs:
            nvi.append(site.coords)
        nv.append(nvi)

    #generate structure coordinate lists
    species_list = []
    s1 = []
    for site in struct1:
        ind = None
        for i, species in enumerate(species_list):
            if site.species_and_occu == species:
                ind = i
                break

        if ind is not None:
            s1[ind] = np.append(s1[ind], [site.frac_coords], axis=0)
        else:
            s1.append(np.array([site.frac_coords]))
            species_list.append(site.species_and_occu)

    zipped = sorted(zip(s1, species_list), key=lambda x: len(x[0]))
    s1 = [x[0] for x in zipped]
    species_list = [x[1] for x in zipped]
    s2_cart = [[] for i in s1]

    for site in struct2:
        ind = None
        for i, species in enumerate(species_list):
            if site.species_and_occu == species:
                ind = i
                break
            #get cart coords for s2
        if s2_cart[ind] != []:
            s2_cart[ind] = np.append(s2_cart[ind], [site.coords], axis=0)
        else:
            s2_cart[ind] = np.array([site.coords])

    #translate s1
    s1_translation = s1[0][0]

    for i in range(len(species_list)):
        s1[i] = (s1[i] - s1_translation) % 1

    #do permutations of vectors, check for equality
    for a, b, c in itertools.product(nv[0], nv[1], nv[2]):
        if np.linalg.det([a, b, c]) == 0: #invalid lattice
            continue

        nl = Lattice([a, b, c])
        if np.allclose(nl.angles, struct1.lattice.angles, rtol=0, atol=atol):
            #Basis Change into new lattice
            s2 = basis_change(s2_cart, nl)

            for coord in s2[0]:
                t_s2 = []
                for coords in s2:
                    t_s2.append((coords - coord) % 1)
                if cmp_struct(s1, t_s2, frac_tol):
                    return True
    return False
