# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import math
import re
import os
import textwrap
import warnings
from collections import OrderedDict, deque

import six
from six.moves import zip, cStringIO

import numpy as np
from functools import partial
from inspect import getargspec
from itertools import groupby
from pymatgen.core.periodic_table import Element, Specie, get_el_sp
from monty.io import zopen
from pymatgen.util.coord_utils import in_coord_list_pbc, pbc_diff
from monty.string import remove_non_ascii
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.groups import SpaceGroup, SYMM_DATA
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

"""
Wrapper classes for Cif input and output from Structures.
"""

__author__ = "Shyue Ping Ong, Will Richards"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "3.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"

sub_spgrp = partial(re.sub, r"[\s_]", "")

space_groups = {sub_spgrp(k): k
                for k in SYMM_DATA['space_group_encoding'].keys()}

space_groups.update({sub_spgrp(k): k
                     for k in SYMM_DATA['space_group_encoding'].keys()})

_COD_DATA = None


def _get_cod_data():
    global _COD_DATA
    if _COD_DATA is None:
        import pymatgen
        with open(os.path.join(pymatgen.symmetry.__path__[0],
                               "symm_ops.json")) \
                as f:
            import json
            _COD_DATA = json.load(f)

    return _COD_DATA


class CifBlock(object):
    maxlen = 70  # not quite 80 so we can deal with semicolons and things

    def __init__(self, data, loops, header):
        """
        Object for storing cif data. All data is stored in a single dictionary.
        Data inside loops are stored in lists in the data dictionary, and
        information on which keys are grouped together are stored in the loops
        attribute.

        Args:
            data: dict or OrderedDict of data to go into the cif. Values should
                    be convertible to string, or lists of these if the key is
                    in a loop
            loops: list of lists of keys, grouped by which loop they should
                    appear in
            header: name of the block (appears after the data_ on the first
                line)
        """
        self.loops = loops
        self.data = data
        # AJ says: CIF Block names cannot be more than 75 characters or you
        # get an Exception
        self.header = header[:74]

    def __eq__(self, other):
        return self.loops == other.loops \
               and self.data == other.data \
               and self.header == other.header

    def __getitem__(self, key):
        return self.data[key]

    def __str__(self):
        """
        Returns the cif string for the data block
        """
        s = ["data_{}".format(self.header)]
        keys = self.data.keys()
        written = []
        for k in keys:
            if k in written:
                continue
            for l in self.loops:
                # search for a corresponding loop
                if k in l:
                    s.append(self._loop_to_string(l))
                    written.extend(l)
                    break
            if k not in written:
                # k didn't belong to a loop
                v = self._format_field(self.data[k])
                if len(k) + len(v) + 3 < self.maxlen:
                    s.append("{}   {}".format(k, v))
                else:
                    s.extend([k, v])
        return "\n".join(s)

    def _loop_to_string(self, loop):
        s = "loop_"
        for l in loop:
            s += '\n ' + l
        for fields in zip(*[self.data[k] for k in loop]):
            line = "\n"
            for val in map(self._format_field, fields):
                if val[0] == ";":
                    s += line + "\n" + val
                    line = "\n"
                elif len(line) + len(val) + 2 < self.maxlen:
                    line += "  " + val
                else:
                    s += line
                    line = '\n  ' + val
            s += line
        return s

    def _format_field(self, v):
        v = v.__str__().strip()
        if len(v) > self.maxlen:
            return ';\n' + textwrap.fill(v, self.maxlen) + '\n;'
        # add quotes if necessary
        if v == '':
            return '""'
        if (" " in v or v[0] == "_") \
                and not (v[0] == "'" and v[-1] == "'") \
                and not (v[0] == '"' and v[-1] == '"'):
            if "'" in v:
                q = '"'
            else:
                q = "'"
            v = q + v + q
        return v

    @classmethod
    def _process_string(cls, string):
        # remove comments
        string = re.sub("(\s|^)#.*$", "", string, flags=re.MULTILINE)
        # remove empty lines
        string = re.sub("^\s*\n", "", string, flags=re.MULTILINE)
        # remove non_ascii
        string = remove_non_ascii(string)

        # since line breaks in .cif files are mostly meaningless,
        # break up into a stream of tokens to parse, rejoining multiline
        # strings (between semicolons)
        q = deque()
        multiline = False
        ml = []
        # this regex splits on spaces, except when in quotes.
        # starting quotes must not be preceded by non-whitespace
        # (these get eaten by the first expression)
        # ending quotes must not be followed by non-whitespace
        p = re.compile(r'''([^'"\s][\S]*)|'(.*?)'(?!\S)|"(.*?)"(?!\S)''')
        for l in string.splitlines():
            if multiline:
                if l.startswith(";"):
                    multiline = False
                    q.append(('', '', '', ' '.join(ml)))
                    ml = []
                    l = l[1:].strip()
                else:
                    ml.append(l)
                    continue
            if l.startswith(";"):
                multiline = True
                ml.append(l[1:].strip())
            else:
                for s in p.findall(l):
                    # s is tuple. location of the data in the tuple
                    # depends on whether it was quoted in the input
                    q.append(s)
        return q

    @classmethod
    def from_string(cls, string):
        q = cls._process_string(string)
        header = q.popleft()[0][5:]
        data = OrderedDict()
        loops = []
        while q:
            s = q.popleft()
            # cif keys aren't in quotes, so show up in s[0]
            if s[0] == "_eof":
                break
            if s[0].startswith("_"):
                data[s[0]] = "".join(q.popleft())
            elif s[0].startswith("loop_"):
                columns = []
                items = []
                while q:
                    s = q[0]
                    if s[0].startswith("loop_") or not s[0].startswith("_"):
                        break
                    columns.append("".join(q.popleft()))
                    data[columns[-1]] = []
                while q:
                    s = q[0]
                    if s[0].startswith("loop_") or s[0].startswith("_"):
                        break
                    items.append("".join(q.popleft()))
                n = len(items) // len(columns)
                assert len(items) % n == 0
                loops.append(columns)
                for k, v in zip(columns * n, items):
                    data[k].append(v.strip())
            elif "".join(s).strip() != "":
                warnings.warn("Possible error in cif format"
                              " error at {}".format("".join(s).strip()))
        return cls(data, loops, header)


class CifFile(object):
    """
    Reads and parses CifBlocks from a .cif file
    """

    def __init__(self, data, orig_string=None, comment=None):
        """
        Args:
            data (OrderedDict): Of CifBlock objects.å
            orig_string (str): The original cif string.
            comment (str): Comment string.
        """
        self.data = data
        self.orig_string = orig_string
        self.comment = comment or "# generated using pymatgen"

    def __str__(self):
        s = ["%s" % v for v in self.data.values()]
        return self.comment + "\n" + "\n".join(s) + "\n"

    @classmethod
    def from_string(cls, string):
        d = OrderedDict()
        for x in re.split("^\s*data_", "x\n" + string,
                          flags=re.MULTILINE | re.DOTALL)[1:]:

            # Skip over Cif block that contains powder diffraction data.
            # Some elements in this block were missing from CIF files in Springer materials/Pauling file DBs.
            # This block anyway does not contain any structure information, and CifParser was also not parsing it.
            if 'powder_pattern' in re.split("\n", x, 1)[0]:
                continue
            c = CifBlock.from_string("data_" + x)
            d[c.header] = c
        return cls(d, string)

    @classmethod
    def from_file(cls, filename):
        with zopen(filename, "rt") as f:
            return cls.from_string(f.read())


class CifParser(object):
    """
    Parses a cif file

    Args:
        filename (str): Cif filename. bzipped or gzipped cifs are fine too.
        occupancy_tolerance (float): If total occupancy of a site is between 1
            and occupancy_tolerance, the occupancies will be scaled down to 1.
        site_tolerance (float): This tolerance is used to determine if two
            sites are sitting in the same position, in which case they will be
            combined to a single disordered site. Defaults to 1e-4.
    """

    def __init__(self, filename, occupancy_tolerance=1., site_tolerance=1e-4):
        self._occupancy_tolerance = occupancy_tolerance
        self._site_tolerance = site_tolerance
        if isinstance(filename, six.string_types):
            self._cif = CifFile.from_file(filename)
        else:
            self._cif = CifFile.from_string(filename.read())

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
        stream = cStringIO(cif_string)
        return CifParser(stream, occupancy_tolerance)

    def _unique_coords(self, coords_in):
        """
        Generate unique coordinates using coord and symmetry positions.
        """
        coords = []
        for tmp_coord in coords_in:
            for op in self.symmetry_operations:
                coord = op.operate(tmp_coord)
                coord = np.array([i - math.floor(i) for i in coord])
                if not in_coord_list_pbc(coords, coord,
                                         atol=self._site_tolerance):
                    coords.append(coord)
        return coords

    def get_lattice(self, data, length_strings=("a", "b", "c"),
                    angle_strings=("alpha", "beta", "gamma"),
                    lattice_type=None):
        """
        Generate the lattice from the provided lattice parameters. In
        the absence of all six lattice parameters, the crystal system
        and necessary parameters are parsed
        """
        try:

            lengths = [str2float(data["_cell_length_" + i])
                       for i in length_strings]
            angles = [str2float(data["_cell_angle_" + i])
                      for i in angle_strings]
            if not lattice_type:
                return Lattice.from_lengths_and_angles(lengths, angles)

            else:
                return getattr(Lattice, lattice_type)(*(lengths + angles))

        except KeyError:
            # Missing Key search for cell setting
            for lattice_lable in ["_symmetry_cell_setting",
                                  "_space_group_crystal_system"]:
                if data.data.get(lattice_lable):
                    lattice_type = data.data.get(lattice_lable).lower()
                    try:

                        required_args = getargspec(
                            getattr(Lattice, lattice_type)).args

                        lengths = (l for l in length_strings
                                   if l in required_args)
                        angles = (a for a in angle_strings
                                  if a in required_args)
                        return self.get_lattice(data, lengths, angles,
                                                lattice_type=lattice_type)
                    except AttributeError as exc:
                        warnings.warn(exc)

                else:
                    return None

    def get_symops(self, data):
        """
        In order to generate symmetry equivalent positions, the symmetry
        operations are parsed. If the symops are not present, the space
        group symbol is parsed, and symops are generated.
        """
        symops = []
        for symmetry_label in ["_symmetry_equiv_pos_as_xyz",
                               "_symmetry_equiv_pos_as_xyz_",
                               "_space_group_symop_operation_xyz",
                               "_space_group_symop_operation_xyz_"]:
            if data.data.get(symmetry_label):
                xyz = data.data.get(symmetry_label)
                if isinstance(xyz, six.string_types):
                    warnings.warn("A 1-line symmetry op P1 CIF is detected!")
                    xyz = [xyz]
                try:
                    symops = [SymmOp.from_xyz_string(s)
                              for s in xyz]
                    break
                except ValueError:
                    continue
        if not symops:
            # Try to parse symbol
            for symmetry_label in ["_symmetry_space_group_name_H-M",
                                   "_symmetry_space_group_name_H_M",
                                   "_symmetry_space_group_name_H-M_",
                                   "_symmetry_space_group_name_H_M_",
                                   "_space_group_name_Hall",
                                   "_space_group_name_Hall_",
                                   "_space_group_name_H-M_alt",
                                   "_space_group_name_H-M_alt_",
                                   "_symmetry_space_group_name_hall",
                                   "_symmetry_space_group_name_hall_",
                                   "_symmetry_space_group_name_h-m",
                                   "_symmetry_space_group_name_h-m_"]:
                sg = data.data.get(symmetry_label)

                if sg:
                    sg = sub_spgrp(sg)
                    try:
                        spg = space_groups.get(sg)
                        if spg:
                            symops = SpaceGroup(spg).symmetry_ops
                            warnings.warn(
                                "No _symmetry_equiv_pos_as_xyz type key found. "
                                "Spacegroup from %s used." % symmetry_label)
                            break
                    except ValueError:
                        # Ignore any errors
                        pass

                    try:
                        for d in _get_cod_data():
                            if sg == re.sub("\s+", "",
                                            d["hermann_mauguin"]) :
                                xyz = d["symops"]
                                symops = [SymmOp.from_xyz_string(s)
                                          for s in xyz]
                                warnings.warn(
                                    "No _symmetry_equiv_pos_as_xyz type key found. "
                                    "Spacegroup from %s used." % symmetry_label)
                                break
                    except Exception as ex:
                        continue

                    if symops:
                        break
        if not symops:
            # Try to parse International number
            for symmetry_label in ["_space_group_IT_number",
                                   "_space_group_IT_number_",
                                   "_symmetry_Int_Tables_number",
                                   "_symmetry_Int_Tables_number_"]:
                if data.data.get(symmetry_label):
                    try:
                        i = int(str2float(data.data.get(symmetry_label)))
                        symops = SpaceGroup.from_int_number(i).symmetry_ops
                        break
                    except ValueError:
                        continue

        if not symops:
            warnings.warn("No _symmetry_equiv_pos_as_xyz type key found. "
                          "Defaulting to P1.")
            symops = [SymmOp.from_xyz_string(s) for s in ['x', 'y', 'z']]

        return symops

    def parse_oxi_states(self, data):
        """
        Parse oxidation states from data dictionary
        """
        try:
            oxi_states = {
                data["_atom_type_symbol"][i]:
                    str2float(data["_atom_type_oxidation_number"][i])
                for i in range(len(data["_atom_type_symbol"]))}
            # attempt to strip oxidation state from _atom_type_symbol
            # in case the label does not contain an oxidation state
            for i, symbol in enumerate(data["_atom_type_symbol"]):
                oxi_states[re.sub(r"\d?[\+,\-]?$", "", symbol)] = \
                    str2float(data["_atom_type_oxidation_number"][i])

        except (ValueError, KeyError):
            oxi_states = None

        return oxi_states

    def _get_structure(self, data, primitive, substitution_dictionary=None):
        """
        Generate structure from part of the cif.
        """
        # Symbols often representing
        # common representations for elements/water in cif files
        special_symbols = {"D": "D", "Hw": "H", "Ow": "O", "Wat": "O",
                           "wat": "O"}
        elements = [el.symbol for el in Element]

        lattice = self.get_lattice(data)
        self.symmetry_operations = self.get_symops(data)
        oxi_states = self.parse_oxi_states(data)

        coord_to_species = OrderedDict()

        def parse_symbol(sym):

            if substitution_dictionary:
                return substitution_dictionary.get(sym)
            elif sym in ['OH', 'OH2']:
                warnings.warn("Symbol '{}' not recognized".format(sym))
                return ""
            else:
                m = re.findall(r"w?[A-Z][a-z]*", sym)
                if m and m != "?":
                    return m[0]
                return ""

        def get_matching_coord(coord):
            for op in self.symmetry_operations:
                c = op.operate(coord)
                for k in coord_to_species.keys():
                    if np.allclose(pbc_diff(c, k), (0, 0, 0),
                                   atol=self._site_tolerance):
                        return tuple(k)
            return False

        ############################################################
        """
        This part of the code deals with handling formats of data as found in
        CIF files extracted from the Springer Materials/Pauling File
        databases, and that are different from standard ICSD formats.
        """

        # Check to see if "_atom_site_type_symbol" exists, as some test CIFs do
        # not contain this key.
        if "_atom_site_type_symbol" in data.data.keys():

            # Keep a track of which data row needs to be removed.
            # Example of a row: Nb,Zr '0.8Nb + 0.2Zr' .2a .m-3m 0 0 0 1 14
            # 'rhombic dodecahedron, Nb<sub>14</sub>'
            # Without this code, the above row in a structure would be parsed
            # as an ordered site with only Nb (since
            # CifParser would try to parse the first two characters of the
            # label "Nb,Zr") and occupancy=1.
            # However, this site is meant to be a disordered site with 0.8 of
            # Nb and 0.2 of Zr.
            idxs_to_remove = []

            for idx, el_row in enumerate(data["_atom_site_label"]):

                # CIF files from the Springer Materials/Pauling File have
                # switched the label and symbol. Thus, in the
                # above shown example row, '0.8Nb + 0.2Zr' is the symbol.
                # Below, we split the strings on ' + ' to
                # check if the length (or number of elements) in the label and
                # symbol are equal.
                if len(data["_atom_site_type_symbol"][idx].split(' + ')) > \
                        len(data["_atom_site_label"][idx].split(' + ')):

                    # Dictionary to hold extracted elements and occupancies
                    els_occu = {}

                    # parse symbol to get element names and occupancy and store
                    # in "els_occu"
                    symbol_str = data["_atom_site_type_symbol"][idx]
                    symbol_str_lst = symbol_str.split(' + ')
                    for elocc_idx in range(len(symbol_str_lst)):
                        # Remove any bracketed items in the string
                        symbol_str_lst[elocc_idx] = re.sub(
                            '\([0-9]*\)', '', symbol_str_lst[elocc_idx].strip())

                        # Extract element name and its occupancy from the
                        # string, and store it as a
                        # key-value pair in "els_occ".
                        els_occu[str(re.findall('\D+', symbol_str_lst[
                            elocc_idx].strip())[1]).replace('<sup>', '')] = \
                            float('0' + re.findall('\.?\d+', symbol_str_lst[
                                elocc_idx].strip())[1])

                    x = str2float(data["_atom_site_fract_x"][idx])
                    y = str2float(data["_atom_site_fract_y"][idx])
                    z = str2float(data["_atom_site_fract_z"][idx])

                    coord = (x, y, z)
                    # Add each partially occupied element on the site coordinate
                    for et in els_occu:
                        match = get_matching_coord(coord)
                        if not match:
                            coord_to_species[coord] = Composition(
                                {parse_symbol(et): els_occu[parse_symbol(et)]})
                        else:
                            coord_to_species[match] += {
                                parse_symbol(et): els_occu[parse_symbol(et)]}
                    idxs_to_remove.append(idx)

            # Remove the original row by iterating over all keys in the CIF data looking for lists, which indicates
            # multiple data items, one for each row, and remove items from the list that corresponds to the removed row,
            # so that it's not processed by the rest of this function (which would result in an error).
            for cif_key in data.data:
                if type(data.data[cif_key]) == list:
                    for id in sorted(idxs_to_remove, reverse=True):
                        del data.data[cif_key][id]

        ############################################################

        for i in range(len(data["_atom_site_label"])):
            symbol = parse_symbol(data["_atom_site_label"][i])

            if symbol:
                if symbol not in elements and symbol not in special_symbols:
                    symbol = symbol[:2]
            else:
                continue
            # make sure symbol was properly parsed from _atom_site_label
            # otherwise get it from _atom_site_type_symbol
            try:
                if symbol in special_symbols:
                    get_el_sp(special_symbols.get(symbol))
                else:
                    Element(symbol)
            except (KeyError, ValueError):
                # sometimes the site doesn't have the type_symbol.
                # we then hope the type_symbol can be parsed from the label
                if "_atom_site_type_symbol" in data.data.keys():
                    symbol = data["_atom_site_type_symbol"][i]

            if oxi_states is not None:
                if symbol in special_symbols:
                    el = get_el_sp(special_symbols.get(symbol) +
                                   str(oxi_states[symbol]))
                else:
                    el = Specie(symbol, oxi_states.get(symbol, 0))
            else:

                el = get_el_sp(special_symbols.get(symbol, symbol))

            x = str2float(data["_atom_site_fract_x"][i])
            y = str2float(data["_atom_site_fract_y"][i])
            z = str2float(data["_atom_site_fract_z"][i])
            try:
                occu = str2float(data["_atom_site_occupancy"][i])
            except (KeyError, ValueError):
                occu = 1

            if occu > 0:
                coord = (x, y, z)
                match = get_matching_coord(coord)
                if not match:
                    coord_to_species[coord] = Composition({el: occu})
                else:
                    coord_to_species[match] += {el: occu}

        if any([sum(c.values()) > 1 for c in coord_to_species.values()]):
            warnings.warn("Some occupancies sum to > 1! If they are within "
                          "the tolerance, they will be rescaled.")

        allspecies = []
        allcoords = []

        if coord_to_species.items():
            for species, group in groupby(
                    sorted(list(coord_to_species.items()), key=lambda x: x[1]),
                    key=lambda x: x[1]):
                tmp_coords = [site[0] for site in group]

                coords = self._unique_coords(tmp_coords)

                allcoords.extend(coords)
                allspecies.extend(len(coords) * [species])

            # rescale occupancies if necessary
            for i, species in enumerate(allspecies):
                totaloccu = sum(species.values())
                if 1 < totaloccu <= self._occupancy_tolerance:
                    allspecies[i] = species / totaloccu

        if allspecies and len(allspecies) == len(allcoords):
            struct = Structure(lattice, allspecies, allcoords)
            struct = struct.get_sorted_structure()

            if primitive:
                struct = struct.get_primitive_structure()
                struct = struct.get_reduced_structure()
            return struct

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
        for d in self._cif.data.values():
            try:
                s = self._get_structure(d, primitive)
                if s:
                    structures.append(s)
            except (KeyError, ValueError) as exc:
                # Warn the user (Errors should never pass silently)
                # A user reported a problem with cif files produced by Avogadro
                # in which the atomic coordinates are in Cartesian coords.
                warnings.warn(str(exc))
        if len(structures) == 0:
            raise ValueError("Invalid cif file with no structures!")
        return structures

    def as_dict(self):
        d = OrderedDict()
        for k, v in self._cif.data.items():
            d[k] = {}
            for k2, v2 in v.data.items():
                d[k][k2] = v2
        return d


class CifWriter(object):
    """
    A wrapper around CifFile to write CIF files from pymatgen structures.

    Args:
        struct (Structure): structure to write
        symprec (float): If not none, finds the symmetry of the structure
            and writes the cif with symmetry information. Passes symprec
            to the SpacegroupAnalyzer
    """

    def __init__(self, struct, symprec=None):
        format_str = "{:.8f}"

        block = OrderedDict()
        loops = []
        spacegroup = ("P 1", 1)
        if symprec is not None:
            sf = SpacegroupAnalyzer(struct, symprec)
            spacegroup = (sf.get_space_group_symbol(),
                          sf.get_space_group_number())
            # Needs the refined struture when using symprec. This converts
            # primitive to conventional structures, the standard for CIF.
            struct = sf.get_refined_structure()

        latt = struct.lattice
        comp = struct.composition
        no_oxi_comp = comp.element_composition
        block["_symmetry_space_group_name_H-M"] = spacegroup[0]
        for cell_attr in ['a', 'b', 'c']:
            block["_cell_length_" + cell_attr] = format_str.format(
                getattr(latt, cell_attr))
        for cell_attr in ['alpha', 'beta', 'gamma']:
            block["_cell_angle_" + cell_attr] = format_str.format(
                getattr(latt, cell_attr))
        block["_symmetry_Int_Tables_number"] = spacegroup[1]
        block["_chemical_formula_structural"] = no_oxi_comp.reduced_formula
        block["_chemical_formula_sum"] = no_oxi_comp.formula
        block["_cell_volume"] = latt.volume.__str__()

        reduced_comp, fu = no_oxi_comp.get_reduced_composition_and_factor()
        block["_cell_formula_units_Z"] = str(int(fu))

        if symprec is None:
            block["_symmetry_equiv_pos_site_id"] = ["1"]
            block["_symmetry_equiv_pos_as_xyz"] = ["x, y, z"]
        else:
            sf = SpacegroupAnalyzer(struct, symprec)

            symmops = []
            for op in sf.get_symmetry_operations():
                v = op.translation_vector
                symmops.append(SymmOp.from_rotation_and_translation(
                    op.rotation_matrix, v))

            ops = [op.as_xyz_string() for op in symmops]
            block["_symmetry_equiv_pos_site_id"] = \
                ["%d" % i for i in range(1, len(ops) + 1)]
            block["_symmetry_equiv_pos_as_xyz"] = ops

        loops.append(["_symmetry_equiv_pos_site_id",
                      "_symmetry_equiv_pos_as_xyz"])

        contains_oxidation = True
        try:
            symbol_to_oxinum = OrderedDict([
                                               (el.__str__(),
                                                float(el.oxi_state))
                                               for el in sorted(comp.elements)])
        except AttributeError:
            symbol_to_oxinum = OrderedDict([(el.symbol, 0) for el in
                                            sorted(comp.elements)])
            contains_oxidation = False
        if contains_oxidation:
            block["_atom_type_symbol"] = symbol_to_oxinum.keys()
            block["_atom_type_oxidation_number"] = symbol_to_oxinum.values()
            loops.append(["_atom_type_symbol", "_atom_type_oxidation_number"])

        atom_site_type_symbol = []
        atom_site_symmetry_multiplicity = []
        atom_site_fract_x = []
        atom_site_fract_y = []
        atom_site_fract_z = []
        atom_site_label = []
        atom_site_occupancy = []
        count = 1
        if symprec is None:
            for site in struct:
                for sp, occu in site.species_and_occu.items():
                    atom_site_type_symbol.append(sp.__str__())
                    atom_site_symmetry_multiplicity.append("1")
                    atom_site_fract_x.append("{0:f}".format(site.a))
                    atom_site_fract_y.append("{0:f}".format(site.b))
                    atom_site_fract_z.append("{0:f}".format(site.c))
                    atom_site_label.append("{}{}".format(sp.symbol, count))
                    atom_site_occupancy.append(occu.__str__())
                    count += 1
        else:
            # The following just presents a deterministic ordering.
            unique_sites = [
                (sorted(sites, key=lambda s: tuple([abs(x) for x in
                                                    s.frac_coords]))[0],
                 len(sites))
                for sites in sf.get_symmetrized_structure().equivalent_sites
                ]
            for site, mult in sorted(
                    unique_sites,
                    key=lambda t: (t[0].species_and_occu.average_electroneg,
                                   -t[1], t[0].a, t[0].b, t[0].c)):
                for sp, occu in site.species_and_occu.items():
                    atom_site_type_symbol.append(sp.__str__())
                    atom_site_symmetry_multiplicity.append("%d" % mult)
                    atom_site_fract_x.append("{0:f}".format(site.a))
                    atom_site_fract_y.append("{0:f}".format(site.b))
                    atom_site_fract_z.append("{0:f}".format(site.c))
                    atom_site_label.append("{}{}".format(sp.symbol, count))
                    atom_site_occupancy.append(occu.__str__())
                    count += 1

        block["_atom_site_type_symbol"] = atom_site_type_symbol
        block["_atom_site_label"] = atom_site_label
        block["_atom_site_symmetry_multiplicity"] = \
            atom_site_symmetry_multiplicity
        block["_atom_site_fract_x"] = atom_site_fract_x
        block["_atom_site_fract_y"] = atom_site_fract_y
        block["_atom_site_fract_z"] = atom_site_fract_z
        block["_atom_site_occupancy"] = atom_site_occupancy
        loops.append(["_atom_site_type_symbol",
                      "_atom_site_label",
                      "_atom_site_symmetry_multiplicity",
                      "_atom_site_fract_x",
                      "_atom_site_fract_y",
                      "_atom_site_fract_z",
                      "_atom_site_occupancy"])
        d = OrderedDict()
        d[comp.reduced_formula] = CifBlock(block, loops, comp.reduced_formula)
        self._cf = CifFile(d)

    def __str__(self):
        """
        Returns the cif as a string.
        """
        return self._cf.__str__()

    def write_file(self, filename):
        """
        Write the cif file.
        """
        with zopen(filename, "wt") as f:
            f.write(self.__str__())


def str2float(text):
    """
    Remove uncertainty brackets from strings and return the float.
    """
    try:
        return float(re.sub("\(.+\)", "", text))
    except TypeError:
        if isinstance(text, list) and len(text) == 1:
            return float(re.sub("\(.+\)", "", text[0]))
