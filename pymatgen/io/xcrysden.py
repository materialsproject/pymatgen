# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from __future__ import division, unicode_literals, print_function

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


class XSF(object):
    """
    Class for parsing XCrysden files.
    """

    def __init__(self, structure):
        self.structure = structure

    def to_string(self):
        """
        Returns a string with the structure in XSF format
        See http://www.xcrysden.org/doc/XSF.html
        """
        lines = []
        app = lines.append

        app("CRYSTAL")
        app("# Primitive lattice vectors in Angstrom")
        app("PRIMVEC")
        cell = self.structure.lattice.matrix
        for i in range(3):
            app(' %.14f %.14f %.14f' % tuple(cell[i]))

        cart_coords = self.structure.cart_coords
        app("# Cartesian coordinates in Angstrom.")
        app("PRIMCOORD")
        app(" %d 1" % len(cart_coords))

        for a in range(len(cart_coords)):
            sp = "%d" % self.structure.atomic_numbers[a]
            app(sp + ' %20.14f %20.14f %20.14f' % tuple(cart_coords[a]))

        return "\n".join(lines)

    @classmethod
    def from_string(cls, input_string, cls_=None):
        """
        Initialize a `Structure` object from a string with data in XSF format.

        Args:
            input_string: String with the structure in XSF format.
                See http://www.xcrysden.org/doc/XSF.html
            cls_: Structure class to be created. default: pymatgen structure

        """
        # CRYSTAL                                        see (1)
        # these are primitive lattice vectors (in Angstroms)
        # PRIMVEC
        #    0.0000000    2.7100000    2.7100000         see (2)
        #    2.7100000    0.0000000    2.7100000
        #    2.7100000    2.7100000    0.0000000

        # these are conventional lattice vectors (in Angstroms)
        # CONVVEC
        #    5.4200000    0.0000000    0.0000000         see (3)
        #    0.0000000    5.4200000    0.0000000
        #    0.0000000    0.0000000    5.4200000

        # these are atomic coordinates in a primitive unit cell  (in Angstroms)
        # PRIMCOORD
        # 2 1                                            see (4)
        # 16      0.0000000     0.0000000     0.0000000  see (5)
        # 30      1.3550000    -1.3550000    -1.3550000

        lattice, coords, species = [], [], []
        lines = input_string.splitlines()

        for i in range(len(lines)):
            if "PRIMVEC" in lines[i]:
                for j in range(i+1, i+4):
                    lattice.append([float(c) for c in lines[j].split()])

            if "PRIMCOORD" in lines[i]:
                num_sites = int(lines[i+1].split()[0])

                for j in range(i+2, i+2+num_sites):
                    tokens = lines[j].split()
                    species.append(int(tokens[0]))
                    coords.append([float(j) for j in tokens[1:4]])
                break
        else:
            raise ValueError("Invalid XSF data")

        if cls_ is None:
            from pymatgen.core.structure import Structure
            cls_ = Structure

        s = cls_(lattice, species, coords, coords_are_cartesian=True)
        return XSF(s)
