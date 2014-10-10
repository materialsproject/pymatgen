# coding: utf-8

from __future__ import division, unicode_literals

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


import math
import re
import textwrap
import warnings
from collections import OrderedDict, deque

import six
from six.moves import zip, cStringIO

import numpy as np

from pymatgen.core.periodic_table import Element, Specie
from monty.io import zopen
from pymatgen.util.coord_utils import in_coord_list_pbc
from monty.string import remove_non_ascii
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


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
                #search for a corresponding loop
                if k in l:
                    s.append(self._loop_to_string(l))
                    written.extend(l)
                    break
            if k not in written:
                #k didn't belong to a loop
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
        v = str(v).strip()
        if len(v) > self.maxlen:
            return ';\n' + textwrap.fill(v, self.maxlen) + '\n;'
        #add quotes if necessary
        if " " in v and not (v[0] == "'" and v[-1] == "'") \
                and not (v[0] == '"' and v[-1] == '"'):
            if "'" in v:
                q = '"'
            else:
                q = "'"
            v = q + v + q
        return v

    @classmethod
    def _process_string(cls, string):
        #remove comments
        string = re.sub("#.*", "", string)
        #remove empty lines
        string = re.sub("^\s*\n", "", string, flags=re.MULTILINE)
        #remove whitespaces at beginning of lines
        string = re.sub("^\s*", "", string, flags=re.MULTILINE)
        #remove non_ascii
        string = remove_non_ascii(string)
        
        #since line breaks in .cif files are mostly meaningless,
        #break up into a stream of tokens to parse, rejoining multiline
        #strings (between semicolons)
        q = deque()
        multiline = False
        ml = []
        #this regex splits on spaces, except when in quotes.
        #it also ignores single quotes when surrounded by non-whitespace
        #since they are sometimes used in author names
        p = re.compile(r'''([^'"\s]+)|'((?:\S'\S|[^'])*)'|"([^"]*)"''')
        for l in string.splitlines():
            if multiline:
                if l.startswith(";"):
                    multiline = False
                    q.append(" ".join(ml))
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
                    q.append(''.join(s))
        return q

    @classmethod
    def from_string(cls, string):
        q = cls._process_string(string)
        header = q.popleft()[5:]
        data = OrderedDict()
        loops = []
        while q:
            s = q.popleft()
            if s == "_eof":
                break
            if s.startswith("_"):
                data[s] = q.popleft()
            elif s.startswith("loop_"):
                columns = []
                items = []
                while q:
                    s = q[0]
                    if s.startswith("loop_") or not s.startswith("_"):
                        break
                    columns.append(q.popleft())
                    data[columns[-1]] = []
                while q:
                    s = q[0]
                    if s.startswith("loop_") or s.startswith("_"):
                        break
                    items.append(q.popleft())
                n = len(items) // len(columns)
                assert len(items) % n == 0
                loops.append(columns)
                for k, v in zip(columns * n, items):
                    data[k].append(v.strip())
            elif s.strip() != "":
                warnings.warn("Possible error in cif format"
                              " error at {}".format(s.strip()))
        return cls(data, loops, header)


class CifFile(object):
    """
    Reads and parses CifBlocks from a .cif file
    """

    def __init__(self, data, orig_string=None):
        """
        Args:
            data: OrderedDict of CifBlock objects
            string: The original cif string
        """
        self.data = data
        self.orig_string = orig_string

    def __str__(self):
        s = ["%s" % v for v in self.data.values()]
        comment = "#generated using pymatgen\n"
        return comment + "\n".join(s)+"\n"

    @classmethod
    def from_string(cls, string):
        d = OrderedDict()
        for x in re.split("^data_", "x\n"+string, 
                          flags=re.MULTILINE | re.DOTALL)[1:]:
            c = CifBlock.from_string("data_"+x)
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
    """

    def __init__(self, filename, occupancy_tolerance=1.):
        self._occupancy_tolerance = occupancy_tolerance
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
        self.symmetry_operations = [SymmOp.from_xyz_string(s) for s in sympos]

        def parse_symbol(sym):
            m = re.search("([A-Z][a-z]*)", sym)
            if m:
                return m.group(1)
            return ""

        try:
            oxi_states = {
                data["_atom_type_symbol"][i]:
                str2float(data["_atom_type_oxidation_number"][i])
                for i in range(len(data["_atom_type_symbol"]))}
        except (ValueError, KeyError):
            oxi_states = None

        coord_to_species = OrderedDict()

        for i in range(len(data["_atom_site_type_symbol"])):
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
                for key, value in six.iteritems(species):
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
        for d in self._cif.data.values():
            try:
                structures.append(self._get_structure(d, primitive))
            except KeyError as exc:
                # Warn the user (Errors should never pass silently)
                # A user reported a problem with cif files produced by Avogadro
                # in which the atomic coordinates are in Cartesian coords.
                warnings.warn(str(exc))

        return structures

    def as_dict(self):
        d = OrderedDict()
        for k, v in self._cif.data.items():
            d[k] = {}
            for k2, v2 in v.data.items():
                d[k][k2] = v2
        return d


class CifWriter:
    """
    A wrapper around CifFile to write CIF files from pymatgen structures.

    Args:
        struct (Structure): A pymatgen.core.structure.Structure object.
        find_spacegroup (bool): Whether to find spacegroup.
            If so, spacegroup information is written.
    """

    def __init__(self, struct, find_spacegroup=False, symprec=None):
        """
        Args:
            struct (Structure): structure to write
            find_spacegroup (bool): whether to try to determine the spacegroup
            symprec (float): If not none, finds the symmetry of the structure and
                writes the cif with symmetry information. Passes symprec to the
                SpacegroupAnalyzer
        """
        format_str = "{:.8f}"
        
        block = OrderedDict()
        loops = []
        latt = struct.lattice
        comp = struct.composition
        no_oxi_comp = comp.element_composition
        spacegroup = ("P 1", 1)
        if find_spacegroup:
            sf = SpacegroupAnalyzer(struct, 0.001)
            spacegroup = (sf.get_spacegroup_symbol(),
                          sf.get_spacegroup_number())
        block["_symmetry_space_group_name_H-M"] = spacegroup[0]
        for cell_attr in ['a', 'b', 'c']:
            block["_cell_length_" + cell_attr] = format_str.format(
                getattr(latt, cell_attr))
        for cell_attr in ['alpha', 'beta', 'gamma']:
            block["_cell_angle_" + cell_attr] = format_str.format(
                getattr(latt, cell_attr))
        block["_symmetry_Int_Tables_number"] = spacegroup[1]
        block["_chemical_formula_structural"] = str(no_oxi_comp
                                                    .reduced_formula)
        block["_chemical_formula_sum"] = str(no_oxi_comp.formula)
        block["_cell_volume"] = str(latt.volume)

        reduced_comp = no_oxi_comp.reduced_composition
        el = no_oxi_comp.elements[0]
        amt = comp[el]
        fu = int(amt / reduced_comp[Element(el.symbol)])

        block["_cell_formula_units_Z"] = str(fu)

        if symprec is None:
            block["_symmetry_equiv_pos_site_id"] = ["1"]
            block["_symmetry_equiv_pos_as_xyz"] = ["x, y, z"]
        else:
            sf = SpacegroupAnalyzer(struct, symprec)
            ops = [op.as_xyz_string() for op in sf.get_symmetry_operations()]
            block["_symmetry_equiv_pos_site_id"] = \
                ["%d" % i for i in range(1, len(ops) + 1)]
            block["_symmetry_equiv_pos_as_xyz"] = ops

        loops.append(["_symmetry_equiv_pos_site_id",
                      "_symmetry_equiv_pos_as_xyz"])

        contains_oxidation = True
        try:
            symbol_to_oxinum = OrderedDict([
                (str(el), float(el.oxi_state))
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
                    atom_site_type_symbol.append(str(sp))
                    atom_site_symmetry_multiplicity.append("1")
                    atom_site_fract_x.append("{0:f}".format(site.a))
                    atom_site_fract_y.append("{0:f}".format(site.b))
                    atom_site_fract_z.append("{0:f}".format(site.c))
                    atom_site_label.append("{}{}".format(sp.symbol, count))
                    atom_site_occupancy.append(str(occu))
                    count += 1
        else:
            for group in sf.get_symmetrized_structure().equivalent_sites:
                site = group[0]
                for sp, occu in site.species_and_occu.items():
                    atom_site_type_symbol.append(str(sp))
                    atom_site_symmetry_multiplicity.append(str(len(group)))
                    atom_site_fract_x.append("{0:f}".format(site.a))
                    atom_site_fract_y.append("{0:f}".format(site.b))
                    atom_site_fract_z.append("{0:f}".format(site.c))
                    atom_site_label.append("{}{}".format(sp.symbol, count))
                    atom_site_occupancy.append(str(occu))
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
    return float(re.sub("\(.+\)", "", text))
