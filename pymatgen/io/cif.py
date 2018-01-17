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

try:
    from inspect import getfullargspec as getargspec
except ImportError:
    from inspect import getargspec
from itertools import groupby
from pymatgen.core.periodic_table import Element, Specie, get_el_sp, DummySpecie
from monty.io import zopen
from pymatgen.util.coord import in_coord_list_pbc, pbc_diff, \
    find_in_coord_list_pbc
from monty.string import remove_non_ascii
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.groups import SpaceGroup, SYMM_DATA
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.electronic_structure.core import Magmom
from pymatgen.core.operations import MagSymmOp
from pymatgen.symmetry.maggroups import MagneticSpaceGroup

try:
    from pybtex.database import BibliographyData, Entry
except ImportError:
    warnings.warn("Please install optional dependency pybtex if you"
                  "want to extract references from CIF files.")

"""
Wrapper classes for Cif input and output from Structures.
"""

__author__ = "Shyue Ping Ong, Will Richards, Matthew Horton"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "4.0"
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
        string = re.sub(r"(\s|^)#.*$", "", string, flags=re.MULTILINE)
        # remove empty lines
        string = re.sub(r"^\s*\n", "", string, flags=re.MULTILINE)
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
    Reads and parses CifBlocks from a .cif file or string
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
        for x in re.split(r"^\s*data_", "x\n" + string,
                          flags=re.MULTILINE | re.DOTALL)[1:]:

            # Skip over Cif block that contains powder diffraction data.
            # Some elements in this block were missing from CIF files in
            # Springer materials/Pauling file DBs.
            # This block anyway does not contain any structure information, and
            # CifParser was also not parsing it.
            if 'powder_pattern' in re.split(r"\n", x, 1)[0]:
                continue
            c = CifBlock.from_string("data_" + x)
            d[c.header] = c
        return cls(d, string)

    @classmethod
    def from_file(cls, filename):
        with zopen(filename, "rt", errors="replace") as f:
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

        # store if CIF contains features from non-core CIF dictionaries
        # e.g. magCIF
        self.feature_flags = {}

        def is_magcif():
            """
            Checks to see if file appears to be a magCIF file (heuristic).
            """
            # Doesn't seem to be a canonical way to test if file is magCIF or
            # not, so instead check for magnetic symmetry datanames
            prefixes = ['_space_group_magn', '_atom_site_moment',
                        '_space_group_symop_magn']
            for d in self._cif.data.values():
                for k in d.data.keys():
                    for prefix in prefixes:
                        if prefix in k:
                            return True
            return False

        self.feature_flags['magcif'] = is_magcif()

        def is_magcif_incommensurate():
            """
            Checks to see if file contains an incommensurate magnetic
            structure (heuristic).
            """
            # Doesn't seem to be a canonical way to test if magCIF file
            # describes incommensurate strucure or not, so instead check
            # for common datanames
            if not self.feature_flags["magcif"]:
                return False
            prefixes = ['_cell_modulation_dimension', '_cell_wave_vector']
            for d in self._cif.data.values():
                for k in d.data.keys():
                    for prefix in prefixes:
                        if prefix in k:
                            return True
            return False

        self.feature_flags['magcif_incommensurate'] = is_magcif_incommensurate()

        for k in self._cif.data.keys():
            # pass individual CifBlocks to _sanitize_data
            self._cif.data[k] = self._sanitize_data(self._cif.data[k])

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

    def _sanitize_data(self, data):
        """
        Some CIF files do not conform to spec. This function corrects
        known issues, particular in regards to Springer materials/
        Pauling files.

        This function is here so that CifParser can assume its
        input conforms to spec, simplifying its implementation.
        :param data: CifBlock
        :return: data CifBlock
        """

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

            new_atom_site_label = []
            new_atom_site_type_symbol = []
            new_atom_site_occupancy = []
            new_fract_x = []
            new_fract_y = []
            new_fract_z = []

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
                            r'\([0-9]*\)', '',
                            symbol_str_lst[elocc_idx].strip())

                        # Extract element name and its occupancy from the
                        # string, and store it as a
                        # key-value pair in "els_occ".
                        els_occu[str(re.findall(r'\D+', symbol_str_lst[
                            elocc_idx].strip())[1]).replace('<sup>', '')] = \
                            float('0' + re.findall(r'\.?\d+', symbol_str_lst[
                                elocc_idx].strip())[1])

                    x = str2float(data["_atom_site_fract_x"][idx])
                    y = str2float(data["_atom_site_fract_y"][idx])
                    z = str2float(data["_atom_site_fract_z"][idx])

                    for et, occu in els_occu.items():
                        # new atom site labels have 'fix' appended
                        new_atom_site_label.append(
                            et + '_fix' + str(len(new_atom_site_label)))
                        new_atom_site_type_symbol.append(et)
                        new_atom_site_occupancy.append(str(occu))
                        new_fract_x.append(str(x))
                        new_fract_y.append(str(y))
                        new_fract_z.append(str(z))

                    idxs_to_remove.append(idx)

            # Remove the original row by iterating over all keys in the CIF
            # data looking for lists, which indicates
            # multiple data items, one for each row, and remove items from the
            # list that corresponds to the removed row,
            # so that it's not processed by the rest of this function (which
            # would result in an error).
            for original_key in data.data:
                if isinstance(data.data[original_key], list):
                    for id in sorted(idxs_to_remove, reverse=True):
                        del data.data[original_key][id]

            if len(idxs_to_remove) > 0:
                data.data["_atom_site_label"] += new_atom_site_label
                data.data["_atom_site_type_symbol"] += new_atom_site_type_symbol
                data.data["_atom_site_occupancy"] += new_atom_site_occupancy
                data.data["_atom_site_fract_x"] += new_fract_x
                data.data["_atom_site_fract_y"] += new_fract_y
                data.data["_atom_site_fract_z"] += new_fract_z

        """
        This fixes inconsistencies in naming of several magCIF tags
        as a result of magCIF being in widespread use prior to
        specification being finalized (on advice of Branton Campbell).
        """

        if self.feature_flags["magcif"]:

            # CIF-1 style has all underscores, interim standard
            # had period before magn instead of before the final
            # component (e.g. xyz)
            # we want to standardize on a specific key, to simplify
            # parsing code
            correct_keys = ["_space_group_symop_magn_operation.xyz",
                            "_space_group_symop_magn_centering.xyz",
                            "_space_group_magn.name_BNS",
                            "_space_group_magn.number_BNS",
                            "_atom_site_moment_crystalaxis_x",
                            "_atom_site_moment_crystalaxis_y",
                            "_atom_site_moment_crystalaxis_z",
                            "_atom_site_moment_label"]

            # cannot mutate OrderedDict during enumeration,
            # so store changes we want to make
            changes_to_make = {}

            for original_key in data.data:
                for correct_key in correct_keys:
                    # convert to all underscore
                    trial_key = "_".join(correct_key.split("."))
                    test_key = "_".join(original_key.split("."))
                    if trial_key == test_key:
                        changes_to_make[correct_key] = original_key

            # make changes
            for correct_key, original_key in changes_to_make.items():
                data.data[correct_key] = data.data[original_key]

            # renamed_keys maps interim_keys to final_keys
            renamed_keys = {
                "_magnetic_space_group.transform_to_standard_Pp_abc":
                    "_space_group_magn.transform_BNS_Pp_abc"}
            changes_to_make = {}

            for interim_key, final_key in renamed_keys.items():
                if data.data.get(interim_key):
                    changes_to_make[final_key] = interim_key
            for final_key, interim_key in changes_to_make.items():
                data.data[final_key] = data.data[interim_key]

        return data

    def _unique_coords(self, coords_in, magmoms_in=None, lattice=None):
        """
        Generate unique coordinates using coord and symmetry positions
        and also their corresponding magnetic moments, if supplied.
        """
        coords = []
        if magmoms_in:
            magmoms = []
            if len(magmoms_in) != len(coords_in):
                raise ValueError
            for tmp_coord, tmp_magmom in zip(coords_in, magmoms_in):
                for op in self.symmetry_operations:
                    coord = op.operate(tmp_coord)
                    coord = np.array([i - math.floor(i) for i in coord])
                    if isinstance(op, MagSymmOp):
                        # Up to this point, magmoms have been defined relative
                        # to crystal axis. Now convert to Cartesian and into
                        # a Magmom object.
                        magmom = Magmom.from_moment_relative_to_crystal_axes(
                            op.operate_magmom(tmp_magmom),
                            lattice=lattice
                        )
                    else:
                        magmom = Magmom(tmp_magmom)
                    if not in_coord_list_pbc(coords, coord,
                                             atol=self._site_tolerance):
                        coords.append(coord)
                        magmoms.append(magmom)
            return coords, magmoms
        else:
            for tmp_coord in coords_in:
                for op in self.symmetry_operations:
                    coord = op.operate(tmp_coord)
                    coord = np.array([i - math.floor(i) for i in coord])
                    if not in_coord_list_pbc(coords, coord,
                                             atol=self._site_tolerance):
                        coords.append(coord)
            return coords, [Magmom(0)] * len(coords)  # return dummy magmoms

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
                            if sg == re.sub(r"\s+", "",
                                            d["hermann_mauguin"]):
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

    def get_magsymops(self, data):
        """
        Equivalent to get_symops except for magnetic symmetry groups.
        Separate function since additional operation for time reversal symmetry
        (which changes magnetic moments on sites) needs to be returned.
        """
        magsymmops = []

        # check to see if magCIF file explicitly contains magnetic symmetry operations
        if data.data.get("_space_group_symop_magn_operation.xyz"):

            xyzt = data.data.get("_space_group_symop_magn_operation.xyz")
            if isinstance(xyzt, six.string_types):
                xyzt = [xyzt]
            magsymmops = [MagSymmOp.from_xyzt_string(s) for s in xyzt]

            if data.data.get("_space_group_symop_magn_centering.xyz"):

                xyzt = data.data.get("_space_group_symop_magn_centering.xyz")
                if isinstance(xyzt, six.string_types):
                    xyzt = [xyzt]
                centering_symops = [MagSymmOp.from_xyzt_string(s) for s in xyzt]

                all_ops = []
                for op in magsymmops:
                    for centering_op in centering_symops:
                        new_translation = [i - np.floor(i) for i
                                           in
                                           op.translation_vector + centering_op.translation_vector]
                        new_time_reversal = op.time_reversal * centering_op.time_reversal
                        all_ops.append(
                            MagSymmOp.from_rotation_and_translation_and_time_reversal(
                                rotation_matrix=op.rotation_matrix,
                                translation_vec=new_translation,
                                time_reversal=new_time_reversal))
                magsymmops = all_ops

        # else check to see if it specifies a magnetic space group
        elif data.data.get("_space_group_magn.name_BNS") or data.data.get(
                "_space_group_magn.number_BNS"):

            if data.data.get("_space_group_magn.name_BNS"):
                # get BNS label for MagneticSpaceGroup()
                id = data.data.get("_space_group_magn.name_BNS")
            else:
                # get BNS number for MagneticSpaceGroup()
                # by converting string to list of ints
                id = list(map(int, (
                data.data.get("_space_group_magn.number_BNS").split("."))))

            msg = MagneticSpaceGroup(id)

            if data.data.get("_space_group_magn.transform_BNS_Pp_abc"):
                if data.data.get(
                        "_space_group_magn.transform_BNS_Pp_abc") != "a,b,c;0,0,0":
                    return NotImplementedError(
                        "Non-standard settings not currently supported.")
            elif data.data.get("_space_group_magn.transform_BNS_Pp"):
                return NotImplementedError(
                    "Incomplete specification to implement.")

            magsymmops = msg.symmetry_ops

        if not magsymmops:
            warnings.warn(
                "No magnetic symmetry detected, using primitive symmetry.")
            magsymmops = [MagSymmOp.from_xyzt_string("x, y, z, 1")]

        return magsymmops

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

    def parse_magmoms(self, data, lattice=None):
        """
        Parse atomic magnetic moments from data dictionary
        """
        if lattice is None:
            raise Exception(
                'Magmoms given in terms of crystal axes in magCIF spec.')
        try:
            magmoms = {
                data["_atom_site_moment_label"][i]:
                    np.array(
                        [str2float(data["_atom_site_moment_crystalaxis_x"][i]),
                         str2float(data["_atom_site_moment_crystalaxis_y"][i]),
                         str2float(data["_atom_site_moment_crystalaxis_z"][i])]
                    )
                for i in range(len(data["_atom_site_moment_label"]))
            }
        except (ValueError, KeyError):
            return None
        return magmoms

    def _get_structure(self, data, primitive):
        """
        Generate structure from part of the cif.
        """

        def parse_symbol(sym):
            # Common representations for elements/water in cif files
            # TODO: fix inconsistent handling of water
            special = {"D": "D", "Hw": "H", "Ow": "O", "Wat": "O",
                       "wat": "O", "OH": "", "OH2": ""}
            m = re.findall(r"w?[A-Z][a-z]*", sym)
            if m and m != "?":
                if sym in special:
                    v = special[sym]
                else:
                    v = special.get(m[0], m[0])
                if len(m) > 1 or (m[0] in special):
                    warnings.warn("{} parsed as {}".format(sym, v))
                return v

        def get_num_implicit_hydrogens(sym):
            num_h = {"Wat": 2, "wat": 2, "O-H": 1}
            return num_h.get(sym[:3], 0)

        lattice = self.get_lattice(data)

        # if magCIF, get magnetic symmetry moments and magmoms
        # else standard CIF, and use empty magmom dict
        if self.feature_flags["magcif_incommensurate"]:
            raise NotImplementedError(
                "Incommensurate structures not currently supported.")
        elif self.feature_flags["magcif"]:
            self.symmetry_operations = self.get_magsymops(data)
            magmoms = self.parse_magmoms(data, lattice=lattice)
        else:
            self.symmetry_operations = self.get_symops(data)
            magmoms = {}

        oxi_states = self.parse_oxi_states(data)

        coord_to_species = OrderedDict()
        coord_to_magmoms = OrderedDict()

        def get_matching_coord(coord):
            keys = list(coord_to_species.keys())
            coords = np.array(keys)
            for op in self.symmetry_operations:
                c = op.operate(coord)
                inds = find_in_coord_list_pbc(coords, c,
                                              atol=self._site_tolerance)
                # cant use if inds, because python is dumb and np.array([0]) evaluates
                # to False
                if len(inds):
                    return keys[inds[0]]
            return False

        for i in range(len(data["_atom_site_label"])):
            try:
                # If site type symbol exists, use it. Otherwise, we use the
                # label.
                symbol = parse_symbol(data["_atom_site_type_symbol"][i])
                num_h = get_num_implicit_hydrogens(
                    data["_atom_site_type_symbol"][i])
            except KeyError:
                symbol = parse_symbol(data["_atom_site_label"][i])
                num_h = get_num_implicit_hydrogens(data["_atom_site_label"][i])
            if not symbol:
                continue

            if oxi_states is not None:
                o_s = oxi_states.get(symbol, 0)
                # use _atom_site_type_symbol if possible for oxidation state
                if "_atom_site_type_symbol" in data.data.keys():
                    oxi_symbol = data["_atom_site_type_symbol"][i]
                    o_s = oxi_states.get(oxi_symbol, o_s)
                try:
                    el = Specie(symbol, o_s)
                except:
                    el = DummySpecie(symbol, o_s)
            else:
                el = get_el_sp(symbol)

            x = str2float(data["_atom_site_fract_x"][i])
            y = str2float(data["_atom_site_fract_y"][i])
            z = str2float(data["_atom_site_fract_z"][i])
            magmom = magmoms.get(data["_atom_site_label"][i],
                                 np.array([0, 0, 0]))

            try:
                occu = str2float(data["_atom_site_occupancy"][i])
            except (KeyError, ValueError):
                occu = 1

            if occu > 0:
                coord = (x, y, z)
                match = get_matching_coord(coord)
                comp_d = {el: occu}
                if num_h > 0:
                    comp_d["H"] = num_h
                comp = Composition(comp_d)
                if not match:
                    coord_to_species[coord] = comp
                    coord_to_magmoms[coord] = magmom
                else:
                    coord_to_species[match] += comp
                    # disordered magnetic not currently supported
                    coord_to_magmoms[match] = None

        sum_occu = [sum(c.values()) for c in coord_to_species.values()
                    if not set(c.elements) == {Element("O"), Element("H")}]
        if any([o > 1 for o in sum_occu]):
            warnings.warn(
                "Some occupancies (%s) sum to > 1! If they are within "
                "the tolerance, they will be rescaled." % str(sum_occu))

        allspecies = []
        allcoords = []
        allmagmoms = []
        allhydrogens = []

        # check to see if magCIF file is disordered
        if self.feature_flags["magcif"]:
            for k, v in coord_to_magmoms.items():
                if v is None:
                    # Proposed solution to this is to instead store magnetic
                    # moments as Specie 'spin' property, instead of site
                    # property, but this introduces ambiguities for end user
                    # (such as unintended use of `spin` and Specie will have
                    # fictious oxidation state).
                    raise NotImplementedError(
                        'Disordered magnetic structures not currently supported.')

        if coord_to_species.items():
            for comp, group in groupby(
                    sorted(list(coord_to_species.items()), key=lambda x: x[1]),
                    key=lambda x: x[1]):
                tmp_coords = [site[0] for site in group]
                tmp_magmom = [coord_to_magmoms[tmp_coord] for tmp_coord in
                              tmp_coords]

                if self.feature_flags["magcif"]:
                    coords, magmoms = self._unique_coords(tmp_coords,
                                                          magmoms_in=tmp_magmom,
                                                          lattice=lattice)
                else:
                    coords, magmoms = self._unique_coords(tmp_coords)

                if set(comp.elements) == {Element("O"), Element("H")}:
                    # O with implicit hydrogens
                    im_h = comp["H"]
                    species = Composition({"O": comp["O"]})
                else:
                    im_h = 0
                    species = comp

                allhydrogens.extend(len(coords) * [im_h])
                allcoords.extend(coords)
                allspecies.extend(len(coords) * [species])
                allmagmoms.extend(magmoms)

            # rescale occupancies if necessary
            for i, species in enumerate(allspecies):
                totaloccu = sum(species.values())
                if 1 < totaloccu <= self._occupancy_tolerance:
                    allspecies[i] = species / totaloccu

        if allspecies and len(allspecies) == len(allcoords) \
                and len(allspecies) == len(allmagmoms):
            site_properties = dict()
            if any(allhydrogens):
                assert len(allhydrogens) == len(allcoords)
                site_properties["implicit_hydrogens"] = allhydrogens

            if self.feature_flags["magcif"]:
                site_properties["magmom"] = allmagmoms

            if len(site_properties) == 0:
                site_properties = None

            struct = Structure(lattice, allspecies, allcoords,
                               site_properties=site_properties)

            struct = struct.get_sorted_structure()

            if primitive and self.feature_flags['magcif']:
                struct = struct.get_primitive_structure(use_site_props=True)
            elif primitive:
                struct = struct.get_primitive_structure()
                struct = struct.get_reduced_structure()

            return struct

    def get_structures(self, primitive=True):
        """
        Return list of structures in CIF file. primitive boolean sets whether a
        conventional cell structure or primitive cell structure is returned.

        Args:
            primitive (bool): Set to False to return conventional unit cells.
                Defaults to True. With magnetic CIF files, will return primitive
                magnetic cell which may be larger than nuclear primitive cell.

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

    def get_bibtex_string(self):
        """
        Get BibTeX reference from CIF file.
        :param data:
        :return: BibTeX string
        """

        bibtex_keys = {'author': ('_publ_author_name', '_citation_author_name'),
                       'title': ('_publ_section_title', '_citation_title'),
                       'journal': ('_journal_name_full', '_journal_name_abbrev',
                                   '_citation_journal_full', '_citation_journal_abbrev'),
                       'volume': ('_journal_volume', '_citation_journal_volume'),
                       'year': ('_journal_year', '_citation_year'),
                       'number': ('_journal_number', '_citation_number'),
                       'page_first': ('_journal_page_first', '_citation_page_first'),
                       'page_last': ('_journal_page_last', '_citation_page_last'),
                       'doi': ('_journal_DOI', '_citation_DOI')}

        entries = {}

        # TODO: parse '_publ_section_references' when it exists?
        # TODO: CIF specification supports multiple citations.

        for idx, data in enumerate(self._cif.data.values()):

            # convert to lower-case keys, some cif files inconsistent
            data = {k.lower(): v for k, v in data.data.items()}

            bibtex_entry = {}

            for field, tags in bibtex_keys.items():
                for tag in tags:
                    if tag in data:
                        bibtex_entry[field] = data[tag]

            # convert to bibtex author format ('and' delimited)
            if 'author' in bibtex_entry:
                bibtex_entry['author'] = ' and '.join(bibtex_entry['author'])

            # convert to bibtex page range format, use empty string if not specified
            if ('page_first' in bibtex_entry) or ('page_last' in bibtex_entry):
                bibtex_entry['pages'] = '{0}--{1}'.format(bibtex_entry.get('page_first', ''),
                                                          bibtex_entry.get('page_last', ''))
                bibtex_entry.pop('page_first', None)  # and remove page_first, page_list if present
                bibtex_entry.pop('page_last', None)


            # cite keys are given as cif-reference-idx in order they are found
            entries['cif-reference-{}'.format(idx)] = Entry('article', list(bibtex_entry.items()))

        return BibliographyData(entries).to_string(bib_format='bibtex')

    def as_dict(self):
        d = OrderedDict()
        for k, v in self._cif.data.items():
            d[k] = {}
            for k2, v2 in v.data.items():
                d[k][k2] = v2
        return d


class CifWriter(object):
    def __init__(self, struct, symprec=None, write_magmoms=False):
        """
        A wrapper around CifFile to write CIF files from pymatgen structures.

        Args:
            struct (Structure): structure to write
            symprec (float): If not none, finds the symmetry of the structure
                and writes the cif with symmetry information. Passes symprec
                to the SpacegroupAnalyzer
            write_magmoms (bool): If True, will write magCIF file. Incompatible
                with symprec
        """

        if write_magmoms and symprec:
            warnings.warn(
                "Magnetic symmetry cannot currently be detected by pymatgen,"
                "disabling symmetry detection.")
            symprec = None

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
        block["_cell_volume"] = "%.8f" % latt.volume

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
        atom_site_moment_label = []
        atom_site_moment_crystalaxis_x = []
        atom_site_moment_crystalaxis_y = []
        atom_site_moment_crystalaxis_z = []
        count = 1
        if symprec is None:
            for site in struct:
                for sp, occu in sorted(site.species_and_occu.items()):
                    atom_site_type_symbol.append(sp.__str__())
                    atom_site_symmetry_multiplicity.append("1")
                    atom_site_fract_x.append("{0:f}".format(site.a))
                    atom_site_fract_y.append("{0:f}".format(site.b))
                    atom_site_fract_z.append("{0:f}".format(site.c))
                    atom_site_label.append("{}{}".format(sp.symbol, count))
                    atom_site_occupancy.append(occu.__str__())

                    magmom = Magmom(
                        site.properties.get('magmom', getattr(sp, 'spin', 0)))
                    if write_magmoms and abs(magmom) > 0:
                        moment = Magmom.get_moment_relative_to_crystal_axes(
                            magmom, latt)
                        atom_site_moment_label.append(
                            "{}{}".format(sp.symbol, count))
                        atom_site_moment_crystalaxis_x.append("%.5f" % moment[0])
                        atom_site_moment_crystalaxis_y.append("%.5f" % moment[1])
                        atom_site_moment_crystalaxis_z.append("%.5f" % moment[2])

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
        if write_magmoms:
            block["_atom_site_moment_label"] = atom_site_moment_label
            block[
                "_atom_site_moment_crystalaxis_x"] = atom_site_moment_crystalaxis_x
            block[
                "_atom_site_moment_crystalaxis_y"] = atom_site_moment_crystalaxis_y
            block[
                "_atom_site_moment_crystalaxis_z"] = atom_site_moment_crystalaxis_z
            loops.append(["_atom_site_moment_label",
                          "_atom_site_moment_crystalaxis_x",
                          "_atom_site_moment_crystalaxis_y",
                          "_atom_site_moment_crystalaxis_z"])
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
        return float(re.sub(r"\(.+\)", "", text))
    except TypeError:
        if isinstance(text, list) and len(text) == 1:
            return float(re.sub(r"\(.+\)", "", text[0]))
    except ValueError as ex:
        if text.strip() == ".":
            return 0
        raise ex
