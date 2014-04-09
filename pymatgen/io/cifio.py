#!/usr/bin/env python

"""
Wrapper classes for Cif input and output from Structures.
"""

from __future__ import division
from pymatgen.symmetry.finder import SymmetryFinder

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"


import re
import cStringIO
import math
import warnings
from collections import OrderedDict

import CifFile
import numpy as np

from pymatgen.core.periodic_table import Element, Specie
from monty.io import zopen
from pymatgen.util.coord_utils import in_coord_list_pbc
from monty.string import remove_non_ascii
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.core.operations import SymmOp


class CifParser(object):
    """
    A wrapper class around PyCifRW to read Cif and convert into a pymatgen
    Structure object.

    Args:
        filename (str): Cif filename. bzipped or gzipped cifs are fine too.
        occupancy_tolerance (float): If total occupancy of a site is between 1
            and occupancy_tolerance, the occupancies will be scaled down to 1.
    """

    def __init__(self, filename, occupancy_tolerance=1.):
        self._occupancy_tolerance = occupancy_tolerance
        if isinstance(filename, basestring):
            with zopen(filename, "r") as f:
                # We use this round-about way to clean up the CIF first.
                stream = cStringIO.StringIO(_clean_cif(f.read()))
                self._cif = CifFile.ReadCif(stream)
        else:
            self._cif = CifFile.ReadCif(filename)

    @staticmethod
    def from_string(cif_string, occupancy_tolerance=1.):
        """
        Creates a CifParser from a string.

        Args:
            cif_string (str): String representation of a CIF.
            occupancy_tolerance (float): If total occupancy of a site is
                between 1 and occupancy_tolerance, the occupancies will be
                scaled down to 1.

        Returns:
            CifParser
        """
        stream = cStringIO.StringIO(_clean_cif(cif_string))
        return CifParser(stream, occupancy_tolerance)

    def _unique_coords(self, coord_in):
        """
        Generate unique coordinates using coord and symmetry positions.
        """
        coords = []
        for op in self.symmetry_operations:
            coord = op.operate(coord_in)
            coord = np.array([i - math.floor(i) for i in coord])
            if not in_coord_list_pbc(coords, coord, atol=1e-3):
                coords.append(coord)
        return coords

    def _get_structure(self, data, primitive):
        """
        Generate structure from part of the cif.
        """
        lengths = [str2float(data["_cell_length_" + i])
                   for i in ["a", "b", "c"]]
        angles = [str2float(data["_cell_angle_" + i])
                  for i in ["alpha", "beta", "gamma"]]
        lattice = Lattice.from_lengths_and_angles(lengths, angles)
        try:
            sympos = data["_symmetry_equiv_pos_as_xyz"]
        except KeyError:
            try:
                sympos = data["_symmetry_equiv_pos_as_xyz_"]
            except KeyError:
                warnings.warn("No _symmetry_equiv_pos_as_xyz type key found. "
                              "Defaulting to P1.")
                sympos = ['x, y, z']
        self.symmetry_operations = parse_symmetry_operations(sympos)

        def parse_symbol(sym):
            m = re.search("([A-Z][a-z]*)", sym)
            if m:
                return m.group(1)
            return ""

        try:
            oxi_states = {
                data["_atom_type_symbol"][i]:
                str2float(data["_atom_type_oxidation_number"][i])
                for i in xrange(len(data["_atom_type_symbol"]))}
        except (ValueError, KeyError):
            oxi_states = None

        coord_to_species = OrderedDict()

        for i in xrange(len(data["_atom_site_type_symbol"])):
            symbol = parse_symbol(data["_atom_site_type_symbol"][i])
            if oxi_states is not None:
                el = Specie(symbol,
                            oxi_states[data["_atom_site_type_symbol"][i]])
            else:
                el = Element(symbol)
            x = str2float(data["_atom_site_fract_x"][i])
            y = str2float(data["_atom_site_fract_y"][i])
            z = str2float(data["_atom_site_fract_z"][i])
            try:
                occu = str2float(data["_atom_site_occupancy"][i])
            except (KeyError, ValueError):
                occu = 1
            if occu > 0:
                coord = (x, y, z)
                if coord not in coord_to_species:
                    coord_to_species[coord] = {el: occu}
                else:
                    coord_to_species[coord][el] = occu

        allspecies = []
        allcoords = []

        for coord, species in coord_to_species.items():
            coords = self._unique_coords(coord)
            allcoords.extend(coords)
            allspecies.extend(len(coords) * [species])

        #rescale occupancies if necessary
        for species in allspecies:
            totaloccu = sum(species.values())
            if 1 < totaloccu <= self._occupancy_tolerance:
                for key, value in species.iteritems():
                    species[key] = value / totaloccu

        struct = Structure(lattice, allspecies, allcoords)
        if primitive:
            struct = struct.get_primitive_structure().get_reduced_structure()
        return struct.get_sorted_structure()

    def get_structures(self, primitive=True):
        """
        Return list of structures in CIF file. primitive boolean sets whether a
        conventional cell structure or primitive cell structure is returned.

        Args:
            primitive (bool): Set to False to return conventional unit cells.
                Defaults to True.

        Returns:
            List of Structures.
        """
        structures = []
        for k, v in self._cif.items():
            try:
                structures.append(self._get_structure(v, primitive))
            except KeyError as exc:
                # Warn the user (Errors should never pass silently)
                # A user reported a problem with cif files produced by Avogadro
                # in which the atomic coordinates are in Cartesian coords.
                warnings.warn(str(exc))

        return structures

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

    Args:
        struct (Structure): A pymatgen.core.structure.Structure object.
        find_spacegroup (bool): Whether to find spacegroup.
            If so, spacegroup information is written.
    """

    def __init__(self, struct, find_spacegroup=False):
        block = CifFile.CifBlock()
        latt = struct.lattice
        comp = struct.composition
        no_oxi_comp = Composition(comp.formula)
        spacegroup = ("P 1", 1)
        if find_spacegroup:
            sf = SymmetryFinder(struct, 0.001)
            spacegroup = (sf.get_spacegroup_symbol(), sf.get_spacegroup_number())
        block["_symmetry_space_group_name_H-M"] = spacegroup[0]
        for cell_attr in ['a', 'b', 'c']:
            block["_cell_length_" + cell_attr] = str(getattr(latt, cell_attr))
        for cell_attr in ['alpha', 'beta', 'gamma']:
            block["_cell_angle_" + cell_attr] = float(getattr(latt, cell_attr))
        block["_chemical_name_systematic"] = "Generated by pymatgen"
        block["_symmetry_Int_Tables_number"] = spacegroup[1]
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
        try:
            symbol_to_oxinum = {str(el): float(el.oxi_state)
                                for el in comp.elements}
        except AttributeError:
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
                atom_site_fract_x.append("{0:f}".format(site.a))
                atom_site_fract_y.append("{0:f}".format(site.b))
                atom_site_fract_z.append("{0:f}".format(site.c))
                atom_site_attached_hydrogens.append("0")
                atom_site_B_iso_or_equiv.append(".")
                atom_site_label.append("{}{}".format(sp.symbol, count))
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


def _clean_cif(s):
    """
    Removes non-ASCII and some unsupported _cgraph fields from the cif
    string
    """
    clean = []
    lines = s.split("\n")
    skip = False
    while len(lines) > 0:
        l = lines.pop(0)
        if skip:
            if l.strip().startswith("_") or l.strip() == "loop_":
                skip = False
            else:
                continue

        if l.strip().startswith("_cgraph"):
            skip = True
        elif not l.strip().startswith("_eof"):
            clean.append(remove_non_ascii(l))

    return "\n".join(clean)


def str2float(text):
    """
    Remove uncertainty brackets from strings and return the float.
    """
    return float(re.sub("\(.+\)", "", text))


def parse_symmetry_operations(symmops_str_list):
    """
    Helper method to parse the symmetry operations.

    Args:
        symmops_str_list ([str]): List of symmops strings of the form
            ['x, y, z', '-x, -y, z', '-y+1/2, x+1/2, z+1/2', ...]

    Returns:
        List of SymmOps
    """
    ops = []
    for op_str in symmops_str_list:
        rot_matrix = np.zeros((3, 3))
        trans = np.zeros(3)
        toks = op_str.strip().split(",")
        for i, tok in enumerate(toks):
            for m in re.finditer("([\+\-]*)\s*([x-z\d]+)/*(\d*)", tok):
                factor = -1 if m.group(1) == "-" else 1
                if m.group(2) in ("x", "y", "z"):
                    j = ord(m.group(2)) - 120
                    rot_matrix[i, j] = factor
                else:
                    num = float(m.group(2))
                    if m.group(3) != "":
                        num /= float(m.group(3))
                    trans[i] = factor * num
        op = SymmOp.from_rotation_and_translation(rot_matrix, trans)
        ops.append(op)
    return ops
