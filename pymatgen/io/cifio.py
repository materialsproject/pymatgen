#!/usr/bin/env python

"""
Wrapper classes for Cif input and output from Structures.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"


import re
import StringIO
import math
import warnings
from collections import OrderedDict

import CifFile
import numpy as np

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.util.io_utils import zopen
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure, Composition
from pymatgen.core.operations import SymmOp


class CifParser(object):
    """
    A wrapper class around PyCifRW to read Cif and convert into a pymatgen
    Structure object.
    """

    def __init__(self, filename, occupancy_tolerance=1.):
        """
        Args:
            filename:
                Cif filename. bzipped or gzipped cifs are fine too.
            occupancy_tolerance:
                If total occupancy of a site is between 1 and
                occupancy_tolerance, the occupancies will be scaled down to 1.
        """
        self._occupancy_tolerance = occupancy_tolerance
        if isinstance(filename, basestring):
            with zopen(filename, "r") as f:
                self._cif = CifFile.ReadCif(f)
        else:
            self._cif = CifFile.ReadCif(filename)

    @staticmethod
    def from_string(cif_string, occupancy_tolerance=1.):
        output = StringIO.StringIO()
        output.write(cif_string)
        output.seek(0)
        return CifParser(output, occupancy_tolerance)

    def _unique_coords(self, coord_in, primitive, lattice, primlattice):
        """
        Generate unique coordinates using coord and symmetry positions.
        """
        coords = list()
        for op in self.symmetry_operations:
            coord = op.operate(coord_in)
            if primitive:
                cart = lattice.get_cartesian_coords(np.array(coord))
                coord = primlattice.get_fractional_coords(cart)
            coord = np.array([i - math.floor(i) for i in coord])
            if not coord_in_list(coord, coords, 1e-3):
                coords.append(coord)
        return coords

    def _get_structure(self, data, primitive):
        """
        Generate structure from part of the cif.
        """
        spacegroup = data.get("_symmetry_space_group_name_H-M", "P1")

        if len(spacegroup) == 0:
            latt_type = "P"
        else:
            latt_type = spacegroup[0]
        lengths = [float_from_str(data["_cell_length_" + i])
                   for i in ["a", "b", "c"]]
        angles = [float_from_str(data["_cell_angle_" + i])
                  for i in ["alpha", "beta", "gamma"]]
        lattice = Lattice.from_lengths_and_angles(lengths, angles)
        primlattice = lattice.get_primitive_lattice(latt_type)
        try:
            sympos = data["_symmetry_equiv_pos_as_xyz"]
        except:
            try:
                sympos = data["_symmetry_equiv_pos_as_xyz_"]
            except:
                warnings.warn("No _symmetry_equiv_pos_as_xyz type key found. "
                              "Defaulting to P1.")
                sympos
        self.symmetry_operations = parse_symmetry_operations(sympos)

        def parse_symbol(sym):
            m = re.search("([A-Z][a-z]*)", sym)
            if m:
                return m.group(1)
            return ""

        #oxi_states = None
        try:
            oxi_states = dict()
            for i in xrange(len(data["_atom_type_symbol"])):
                oxi_states[data["_atom_type_symbol"][i]] = \
                    float_from_str(data["_atom_type_oxidation_number"][i])
        except:
            oxi_states = None

        coord_to_species = OrderedDict()

        for i in xrange(len(data["_atom_site_type_symbol"])):
            symbol = parse_symbol(data["_atom_site_type_symbol"][i])
            if oxi_states is not None:
                el = Specie(symbol,
                            oxi_states[data["_atom_site_type_symbol"][i]])
            else:
                el = Element(symbol)
            x = float_from_str(data["_atom_site_fract_x"][i])
            y = float_from_str(data["_atom_site_fract_y"][i])
            z = float_from_str(data["_atom_site_fract_z"][i])
            try:
                occu = float_from_str(data["_atom_site_occupancy"][i])
            except:
                occu = 1
            if occu > 0:
                coord = (x, y, z)
                if coord not in coord_to_species:
                    coord_to_species[coord] = {el: occu}
                else:
                    coord_to_species[coord][el] = occu

        allspecies = list()
        allcoords = list()

        for coord, species in coord_to_species.items():
            coords = self._unique_coords(coord, primitive, lattice,
                                         primlattice)
            allcoords.extend(coords)
            allspecies.extend(len(coords) * [species])

        #rescale occupancies if necessary
        for species in allspecies:
            totaloccu = sum(species.values())
            if  1 < totaloccu <= self._occupancy_tolerance:
                for key, value in species.iteritems():
                    species[key] = value / totaloccu

        if primitive:
            struct = Structure(primlattice, allspecies, allcoords)
        else:
            struct = Structure(lattice, allspecies, allcoords)
        return struct.get_sorted_structure()

    def get_structures(self, primitive=True):
        """
        Return list of structures in CIF file. primitive boolean sets whether a
        conventional cell structure or primitive cell structure is returned.

        Args:
            primitive:
                Set to False to return conventional unit cells. Defaults to
                True.

        Returns:
            List of Structures.
        """
        return [self._get_structure(v, primitive)
                for k, v in self._cif.items()]

    @property
    def to_dict(self):
        d = OrderedDict()
        for k, v in self._cif.items():
            d[k] = {}
            for k2, v2 in v.items():
                d[k][k2] = v2
        return d


class CifWriter:
    """
    A wrapper around PyCifRW to write CIF files from pymatgen structures.
    """

    def __init__(self, struct):
        """
        Args:
            struct:
                A pymatgen.core.structure.Structure object.
        """
        block = CifFile.CifBlock()
        latt = struct.lattice
        comp = struct.composition
        no_oxi_comp = Composition.from_formula(comp.formula)
        block["_symmetry_space_group_name_H-M"] = "P 1"
        block["_cell_length_a"] = str(latt.a)
        block["_cell_length_b"] = str(latt.b)
        block["_cell_length_c"] = str(latt.c)
        block["_cell_angle_alpha"] = str(latt.alpha)
        block["_cell_angle_beta"] = str(latt.beta)
        block["_cell_angle_gamma"] = str(latt.gamma)
        block["_chemical_name_systematic"] = "Generated by pymatgen"
        block["_symmetry_Int_Tables_number"] = 1
        block["_chemical_formula_structural"] = str(no_oxi_comp
                                                    .reduced_formula)
        block["_chemical_formula_sum"] = str(no_oxi_comp.formula)
        block["_cell_volume"] = str(latt.volume)

        reduced_comp = Composition.from_formula(no_oxi_comp.reduced_formula)
        el = no_oxi_comp.elements[0]
        amt = comp[el]
        fu = int(amt / reduced_comp[Element(el.symbol)])

        block["_cell_formula_units_Z"] = str(fu)
        block.AddCifItem(([["_symmetry_equiv_pos_site_id",
                            "_symmetry_equiv_pos_as_xyz"]],
                          [[["1"], ["x, y, z"]]]))

        contains_oxidation = True
        symbol_to_oxinum = dict()
        try:
            symbol_to_oxinum = {str(el): el.oxi_state for el in comp.elements}
        except:
            symbol_to_oxinum = {el.symbol: 0 for el in comp.elements}
            contains_oxidation = False
        if contains_oxidation:
            block.AddCifItem(([["_atom_type_symbol",
                                "_atom_type_oxidation_number"]],
                              [[symbol_to_oxinum.keys(),
                                symbol_to_oxinum.values()]]))

        atom_site_type_symbol = []
        atom_site_symmetry_multiplicity = []
        atom_site_fract_x = []
        atom_site_fract_y = []
        atom_site_fract_z = []
        atom_site_attached_hydrogens = []
        atom_site_B_iso_or_equiv = []
        atom_site_label = []
        atom_site_occupancy = []
        count = 1
        for site in struct:
            for sp, occu in site.species_and_occu.items():
                atom_site_type_symbol.append(str(sp))
                atom_site_symmetry_multiplicity.append("1")
                atom_site_fract_x.append(str("{0:f}".format(site.a)))
                atom_site_fract_y.append(str("{0:f}".format(site.b)))
                atom_site_fract_z.append(str("{0:f}".format(site.c)))
                atom_site_attached_hydrogens.append("0")
                atom_site_B_iso_or_equiv.append(".")
                atom_site_label.append(str(sp.symbol) + str(count))
                atom_site_occupancy.append(str(occu))
                count += 1

        block["_atom_site_type_symbol"] = atom_site_type_symbol
        block.AddToLoop("_atom_site_type_symbol",
                        {"_atom_site_label": atom_site_label})
        block.AddToLoop("_atom_site_type_symbol",
                        {"_atom_site_symmetry_multiplicity":
                         atom_site_symmetry_multiplicity})
        block.AddToLoop("_atom_site_type_symbol",
                        {"_atom_site_fract_x": atom_site_fract_x})
        block.AddToLoop("_atom_site_type_symbol",
                        {"_atom_site_fract_y": atom_site_fract_y})
        block.AddToLoop("_atom_site_type_symbol",
                        {"_atom_site_fract_z": atom_site_fract_z})
        block.AddToLoop("_atom_site_type_symbol",
                        {"_atom_site_attached_hydrogens":
                         atom_site_attached_hydrogens})
        block.AddToLoop("_atom_site_type_symbol",
                        {"_atom_site_B_iso_or_equiv":
                         atom_site_B_iso_or_equiv})
        block.AddToLoop("_atom_site_type_symbol",
                        {"_atom_site_occupancy": atom_site_occupancy})

        self._cf = CifFile.CifFile()
        # AJ says: CIF Block names cannot be more than 75 characters or you
        # get an Exception
        self._cf[comp.reduced_formula[0:74]] = block

    def __str__(self):
        """
        Returns the cif as a string.
        """
        return str(self._cf)

    def write_file(self, filename):
        """
        Write the cif file.
        """
        with open(filename, "w") as f:
            f.write(self.__str__())


def around_diff_num(a, b):
    """
    Used to compare differences in fractional coordinates, taking into account
    PBC.
    """
    diff_num = abs(a - b)
    return diff_num if diff_num < 0.5 else abs(1 - diff_num)


def coord_in_list(coord, coord_list, tol):
    """
    Helper method to check if coord is already in a list of coords, subject to
    a tolerance.
    """
    for c in coord_list:
        diff = np.array([around_diff_num(c[i], coord[i]) for i in xrange(3)])
        if (diff < tol).all():
            return True
    return False


def float_from_str(text):
    """
    Remove uncertainty brackets from strings and return the float.
    """
    return float(re.sub("\(\d+\)", "", text))


def parse_symmetry_operations(symmops_str):
    ops = []
    for op_str in symmops_str:
        rot_matrix = np.zeros((3, 3))
        trans = np.zeros(3)
        toks = op_str.strip().split(",")
        for i, tok in enumerate(toks):
            for m in re.finditer("([\+\-]*)\s*([x-z\d]+)\/*(\d*)", tok):
                factor = -1 if m.group(1) == "-" else 1
                if m.group(2) in ("x", "y", "z"):
                    j = ord(m.group(2)) - 120
                    rot_matrix[i, j] = factor
                else:
                    num = float(m.group(2))
                    if m.group(3) != "":
                        num = num / float(m.group(3))
                    trans[i] = factor * num
        op = SymmOp.from_rotation_and_translation(rot_matrix, trans)
        ops.append(op)
    return ops
