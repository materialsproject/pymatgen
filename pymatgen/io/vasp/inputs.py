# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import os
import re
import itertools
import warnings
import logging
import math
import glob

import six
import numpy as np
from numpy.linalg import det
from collections import OrderedDict, namedtuple
from hashlib import md5

from monty.io import zopen
from monty.os.path import zpath
from monty.json import MontyDecoder

from enum import Enum
from tabulate import tabulate

import scipy.constants as const

from pymatgen import SETTINGS
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element, get_el_sp
from pymatgen.electronic_structure.core import Magmom
from monty.design_patterns import cached_class
from pymatgen.util.string import str_delimited
from pymatgen.util.io_utils import clean_lines
from monty.json import MSONable

"""
Classes for reading/manipulating/writing VASP input files. All major VASP input
files.
"""

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Rickard Armiento, " + \
             "Vincent L Chevrier, Stephen Dacek"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Jul 16, 2012"


logger = logging.getLogger(__name__)


class Poscar(MSONable):
    """
    Object for representing the data in a POSCAR or CONTCAR file.
    Please note that this current implementation. Most attributes can be set
    directly.

    Args:
        structure (Structure):  Structure object.
        comment (str): Optional comment line for POSCAR. Defaults to unit
            cell formula of structure. Defaults to None.
        selective_dynamics (Nx3 array): bool values for selective dynamics,
            where N is number of sites. Defaults to None.
        true_names (bool): Set to False is the names in the POSCAR are not
            well-defined and ambiguous. This situation arises commonly in
            vasp < 5 where the POSCAR sometimes does not contain element
            symbols. Defaults to True.
        velocities (Nx3 array): Velocities for the POSCAR. Typically parsed
            in MD runs or can be used to initialize velocities.
        predictor_corrector (Nx3 array): Predictor corrector for the POSCAR.
            Typically parsed in MD runs.

    .. attribute:: structure

        Associated Structure.

    .. attribute:: comment

        Optional comment string.

    .. attribute:: true_names

        Boolean indication whether Poscar contains actual real names parsed
        from either a POTCAR or the POSCAR itself.

    .. attribute:: selective_dynamics

        Selective dynamics attribute for each site if available. A Nx3 array of
        booleans.

    .. attribute:: velocities

        Velocities for each site (typically read in from a CONTCAR). A Nx3
        array of floats.

    .. attribute:: predictor_corrector

        Predictor corrector coordinates and derivatives for each site; i.e.
        a list of three 1x3 arrays for each site (typically read in from a MD 
        CONTCAR).

    .. attribute:: predictor_corrector_preamble

        Predictor corrector preamble contains the predictor-corrector key,
        POTIM, and thermostat parameters that precede the site-specic predictor 
        corrector data in MD CONTCAR

    .. attribute:: temperature

        Temperature of velocity Maxwell-Boltzmann initialization. Initialized
        to -1 (MB hasn"t been performed).
    """

    def __init__(self, structure, comment=None, selective_dynamics=None,
                 true_names=True, velocities=None, predictor_corrector=None,
                 predictor_corrector_preamble=None):
        if structure.is_ordered:
            site_properties = {}
            if selective_dynamics:
                site_properties["selective_dynamics"] = selective_dynamics
            if velocities:
                site_properties["velocities"] = velocities
            if predictor_corrector:
                site_properties["predictor_corrector"] = predictor_corrector
            structure = Structure.from_sites(structure)
            self.structure = structure.copy(site_properties=site_properties)
            self.true_names = true_names
            self.comment = structure.formula if comment is None else comment
            self.predictor_corrector_preamble = predictor_corrector_preamble
        else:
            raise ValueError("Structure with partial occupancies cannot be "
                             "converted into POSCAR!")

        self.temperature = -1

    @property
    def velocities(self):
        return self.structure.site_properties.get("velocities")

    @property
    def selective_dynamics(self):
        return self.structure.site_properties.get("selective_dynamics")

    @property
    def predictor_corrector(self):
        return self.structure.site_properties.get("predictor_corrector")

    @velocities.setter
    def velocities(self, velocities):
        self.structure.add_site_property("velocities", velocities)

    @selective_dynamics.setter
    def selective_dynamics(self, selective_dynamics):
        self.structure.add_site_property("selective_dynamics",
                                         selective_dynamics)

    @predictor_corrector.setter
    def predictor_corrector(self, predictor_corrector):
        self.structure.add_site_property("predictor_corrector",
                                         predictor_corrector)

    @property
    def site_symbols(self):
        """
        Sequence of symbols associated with the Poscar. Similar to 6th line in
        vasp 5+ POSCAR.
        """
        syms = [site.specie.symbol for site in self.structure]
        return [a[0] for a in itertools.groupby(syms)]

    @property
    def natoms(self):
        """
        Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.
        """
        syms = [site.specie.symbol for site in self.structure]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]

    def __setattr__(self, name, value):
        if name in ("selective_dynamics", "velocities"):
            if value is not None and len(value) > 0:
                value = np.array(value)
                dim = value.shape
                if dim[1] != 3 or dim[0] != len(self.structure):
                    raise ValueError(name + " array must be same length as" +
                                     " the structure.")
                value = value.tolist()
        super(Poscar, self).__setattr__(name, value)

    @staticmethod
    def from_file(filename, check_for_POTCAR=True, read_velocities=True):
        """
        Reads a Poscar from a file.

        The code will try its best to determine the elements in the POSCAR in
        the following order:
        1. If check_for_POTCAR is True, the code will try to check if a POTCAR
        is in the same directory as the POSCAR and use elements from that by
        default. (This is the VASP default sequence of priority).
        2. If the input file is Vasp5-like and contains element symbols in the
        6th line, the code will use that if check_for_POTCAR is False or there
        is no POTCAR found.
        3. Failing (2), the code will check if a symbol is provided at the end
        of each coordinate.

        If all else fails, the code will just assign the first n elements in
        increasing atomic number, where n is the number of species, to the
        Poscar. For example, H, He, Li, ....  This will ensure at least a
        unique element is assigned to each site and any analysis that does not
        require specific elemental properties should work fine.

        Args:
            filename (str): File name containing Poscar data.
            check_for_POTCAR (bool): Whether to check if a POTCAR is present
                in the same directory as the POSCAR. Defaults to True.
            read_velocities (bool): Whether to read or not velocities if they
                are present in the POSCAR. Default is True.

        Returns:
            Poscar object.
        """
        dirname = os.path.dirname(os.path.abspath(filename))
        names = None
        if check_for_POTCAR:
            potcars = glob.glob(os.path.join(dirname, "*POTCAR*"))
            if potcars:
                try:
                    potcar = Potcar.from_file(sorted(potcars)[0])
                    names = [sym.split("_")[0] for sym in potcar.symbols]
                    [get_el_sp(n) for n in names]  # ensure valid names
                except:
                    names = None
        with zopen(filename, "rt") as f:
            return Poscar.from_string(f.read(), names,
                                      read_velocities=read_velocities)

    @staticmethod
    def from_string(data, default_names=None, read_velocities=True):
        """
        Reads a Poscar from a string.

        The code will try its best to determine the elements in the POSCAR in
        the following order:
        1. If default_names are supplied and valid, it will use those. Usually,
        default names comes from an external source, such as a POTCAR in the
        same directory.
        2. If there are no valid default names but the input file is Vasp5-like
        and contains element symbols in the 6th line, the code will use that.
        3. Failing (2), the code will check if a symbol is provided at the end
        of each coordinate.

        If all else fails, the code will just assign the first n elements in
        increasing atomic number, where n is the number of species, to the
        Poscar. For example, H, He, Li, ....  This will ensure at least a
        unique element is assigned to each site and any analysis that does not
        require specific elemental properties should work fine.

        Args:
            data (str): String containing Poscar data.
            default_names ([str]): Default symbols for the POSCAR file,
                usually coming from a POTCAR in the same directory.
            read_velocities (bool): Whether to read or not velocities if they
                are present in the POSCAR. Default is True.

        Returns:
            Poscar object.
        """
        # "^\s*$" doesn't match lines with no whitespace
        chunks = re.split(r"\n\s*\n", data.rstrip(), flags=re.MULTILINE)
        try:
            if chunks[0] == "":
                chunks.pop(0)
                chunks[0] = "\n" + chunks[0]
        except IndexError:
            raise ValueError("Empty POSCAR")

        # Parse positions
        lines = tuple(clean_lines(chunks[0].split("\n"), False))
        comment = lines[0]
        scale = float(lines[1])
        lattice = np.array([[float(i) for i in line.split()]
                            for line in lines[2:5]])
        if scale < 0:
            # In vasp, a negative scale factor is treated as a volume. We need
            # to translate this to a proper lattice vector scaling.
            vol = abs(det(lattice))
            lattice *= (-scale / vol) ** (1 / 3)
        else:
            lattice *= scale

        vasp5_symbols = False
        try:
            natoms = [int(i) for i in lines[5].split()]
            ipos = 6
        except ValueError:
            vasp5_symbols = True
            symbols = lines[5].split()

            """
            Atoms and number of atoms in POSCAR written with vasp appear on 
            multiple lines when atoms of the same type are not grouped together 
            and more than 20 groups are then defined ...
            
            Example :
            
            Cr16 Fe35 Ni2
               1.00000000000000
                 8.5415010000000002   -0.0077670000000000   -0.0007960000000000
                -0.0077730000000000    8.5224019999999996    0.0105580000000000
                -0.0007970000000000    0.0105720000000000    8.5356889999999996
               Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Ni   Fe   Cr   Fe   Cr
               Fe   Ni   Fe   Cr   Fe
                 1     1     2     4     2     1     1     1     2     1     1     1     4     1     1     1     5     3     6     1
                 2     1     3     2     5
            Direct
              ...
            """
            nlines_symbols = 1
            for nlines_symbols in range(1, 11):
                try:
                    int(lines[5+nlines_symbols].split()[0])
                    break
                except ValueError:
                    pass
            for iline_symbols in range(6, 5+nlines_symbols):
                symbols.extend(lines[iline_symbols].split())
            natoms = []
            iline_natoms_start = 5+nlines_symbols
            for iline_natoms in range(iline_natoms_start,
                                      iline_natoms_start+nlines_symbols):
                natoms.extend([int(i) for i in lines[iline_natoms].split()])
            atomic_symbols = list()
            for i in range(len(natoms)):
                atomic_symbols.extend([symbols[i]] * natoms[i])
            ipos = 5+2*nlines_symbols

        postype = lines[ipos].split()[0]

        sdynamics = False
        # Selective dynamics
        if postype[0] in "sS":
            sdynamics = True
            ipos += 1
            postype = lines[ipos].split()[0]

        cart = postype[0] in "cCkK"
        nsites = sum(natoms)

        # If default_names is specified (usually coming from a POTCAR), use
        # them. This is in line with Vasp"s parsing order that the POTCAR
        # specified is the default used.
        if default_names:
            try:
                atomic_symbols = []
                for i in range(len(natoms)):
                    atomic_symbols.extend([default_names[i]] * natoms[i])
                vasp5_symbols = True
            except IndexError:
                pass

        if not vasp5_symbols:
            ind = 3 if not sdynamics else 6
            try:
                # Check if names are appended at the end of the coordinates.
                atomic_symbols = [l.split()[ind]
                                  for l in lines[ipos + 1:ipos + 1 + nsites]]
                # Ensure symbols are valid elements
                if not all([Element.is_valid_symbol(sym)
                            for sym in atomic_symbols]):
                    raise ValueError("Non-valid symbols detected.")
                vasp5_symbols = True
            except (ValueError, IndexError):
                # Defaulting to false names.
                atomic_symbols = []
                for i in range(len(natoms)):
                    sym = Element.from_Z(i + 1).symbol
                    atomic_symbols.extend([sym] * natoms[i])
                warnings.warn("Elements in POSCAR cannot be determined. "
                              "Defaulting to false names %s." %
                              " ".join(atomic_symbols))

        # read the atomic coordinates
        coords = []
        selective_dynamics = list() if sdynamics else None
        for i in range(nsites):
            toks = lines[ipos + 1 + i].split()
            crd_scale = scale if cart else 1
            coords.append([float(j) * crd_scale for j in toks[:3]])
            if sdynamics:
                selective_dynamics.append([tok.upper()[0] == "T"
                                           for tok in toks[3:6]])

        struct = Structure(lattice, atomic_symbols, coords,
                           to_unit_cell=False, validate_proximity=False,
                           coords_are_cartesian=cart)

        if read_velocities:
            # Parse velocities if any
            velocities = []
            if len(chunks) > 1:
                for line in chunks[1].strip().split("\n"):
                    velocities.append([float(tok) for tok in line.split()])

            # Parse the predictor-corrector data
            predictor_corrector = []
            predictor_corrector_preamble = None

            if len(chunks) > 2:
                lines = chunks[2].strip().split("\n")
                # There are 3 sets of 3xN Predictor corrector parameters
                # So can't be stored as a single set of "site_property"

                # First line in chunk is a key in CONTCAR
                # Second line is POTIM
                # Third line is the thermostat parameters
                predictor_corrector_preamble = (lines[0] + "\n" + lines[1]
                                                + "\n" + lines[2])
                # Rest is three sets of parameters, each set contains
                # x, y, z predictor-corrector parameters for every atom in orde
                lines = lines[3:]
                for st in range(nsites):
                    d1 = [float(tok) for tok in lines[st].split()]
                    d2 = [float(tok) for tok in lines[st+nsites].split()]
                    d3 = [float(tok) for tok in lines[st+2*nsites].split()]
                    predictor_corrector.append([d1,d2,d3])
        else:
            velocities = None
            predictor_corrector = None
            predictor_corrector_preamble = None

        return Poscar(struct, comment, selective_dynamics, vasp5_symbols,
                      velocities=velocities,
                      predictor_corrector=predictor_corrector,
                      predictor_corrector_preamble=predictor_corrector_preamble)

    def get_string(self, direct=True, vasp4_compatible=False,
                   significant_figures=6):
        """
        Returns a string to be written as a POSCAR file. By default, site
        symbols are written, which means compatibility is for vasp >= 5.

        Args:
            direct (bool): Whether coordinates are output in direct or
                cartesian. Defaults to True.
            vasp4_compatible (bool): Set to True to omit site symbols on 6th
                line to maintain backward vasp 4.x compatibility. Defaults
                to False.
            significant_figures (int): No. of significant figures to
                output all quantities. Defaults to 6. Note that positions are
                output in fixed point, while velocities are output in
                scientific format.

        Returns:
            String representation of POSCAR.
        """

        # This corrects for VASP really annoying bug of crashing on lattices
        # which have triple product < 0. We will just invert the lattice
        # vectors.
        latt = self.structure.lattice
        if np.linalg.det(latt.matrix) < 0:
            latt = Lattice(-latt.matrix)

        format_str = "{{:.{0}f}}".format(significant_figures)
        lines = [self.comment, "1.0"]
        for v in latt.matrix:
            lines.append(" ".join([format_str.format(c) for c in v]))

        if self.true_names and not vasp4_compatible:
            lines.append(" ".join(self.site_symbols))
        lines.append(" ".join([str(x) for x in self.natoms]))
        if self.selective_dynamics:
            lines.append("Selective dynamics")
        lines.append("direct" if direct else "cartesian")

        for (i, site) in enumerate(self.structure):
            coords = site.frac_coords if direct else site.coords
            line = " ".join([format_str.format(c) for c in coords])
            if self.selective_dynamics is not None:
                sd = ["T" if j else "F" for j in self.selective_dynamics[i]]
                line += " %s %s %s" % (sd[0], sd[1], sd[2])
            line += " " + site.species_string
            lines.append(line)

        if self.velocities:
            try:
                lines.append("")
                for v in self.velocities:
                    lines.append(" ".join([format_str.format(i) for i in v]))
            except:
                warnings.warn("Velocities are missing or corrupted.")

        if self.predictor_corrector:
            lines.append("")
            if self.predictor_corrector_preamble:
                lines.append(self.predictor_corrector_preamble)
                pred = np.array(self.predictor_corrector)
                for col in range(3):
                    for z in pred[:,col]:
                        lines.append(" ".join([format_str.format(i) for i in z]))
            else:
                warnings.warn(
                    "Preamble information missing or corrupt. " 
                    "Writing Poscar with no predictor corrector data.")

        return "\n".join(lines) + "\n"

    def __repr__(self):
        return self.get_string()

    def __str__(self):
        """
        String representation of Poscar file.
        """
        return self.get_string()

    def write_file(self, filename, **kwargs):
        """
        Writes POSCAR to a file. The supported kwargs are the same as those for
        the Poscar.get_string method and are passed through directly.
        """
        with zopen(filename, "wt") as f:
            f.write(self.get_string(**kwargs))

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "structure": self.structure.as_dict(),
                "true_names": self.true_names,
                "selective_dynamics": np.array(
                    self.selective_dynamics).tolist(),
                "velocities": self.velocities,
                "predictor_corrector": self.predictor_corrector,
                "comment": self.comment}

    @classmethod
    def from_dict(cls, d):
        return Poscar(Structure.from_dict(d["structure"]),
                      comment=d["comment"],
                      selective_dynamics=d["selective_dynamics"],
                      true_names=d["true_names"],
                      velocities=d.get("velocities", None),
                      predictor_corrector=d.get("predictor_corrector", None))

    def set_temperature(self, temperature):
        """
        Initializes the velocities based on Maxwell-Boltzmann distribution.
        Removes linear, but not angular drift (same as VASP)

        Scales the energies to the exact temperature (microcanonical ensemble)
        Velocities are given in A/fs. This is the vasp default when
        direct/cartesian is not specified (even when positions are given in
        direct coordinates)

        Overwrites imported velocities, if any.

        Args:
            temperature (float): Temperature in Kelvin.
        """
        # mean 0 variance 1
        velocities = np.random.randn(len(self.structure), 3)

        # in AMU, (N,1) array
        atomic_masses = np.array([site.specie.atomic_mass.to("kg")
                                  for site in self.structure])
        dof = 3 * len(self.structure) - 3

        # scale velocities due to atomic masses
        # mean 0 std proportional to sqrt(1/m)
        velocities /= atomic_masses[:, np.newaxis] ** (1 / 2)

        # remove linear drift (net momentum)
        velocities -= np.average(atomic_masses[:, np.newaxis] * velocities,
                                 axis=0) / np.average(atomic_masses)

        # scale velocities to get correct temperature
        energy = np.sum(1 / 2 * atomic_masses *
                        np.sum(velocities ** 2, axis=1))
        scale = (temperature * dof / (2 * energy / const.k)) ** (1 / 2)

        velocities *= scale * 1e-5  # these are in A/fs

        self.temperature = temperature
        try:
            del self.structure.site_properties["selective_dynamics"]
        except KeyError:
            pass

        try:
            del self.structure.site_properties["predictor_corrector"]
        except KeyError:
            pass
        # returns as a list of lists to be consistent with the other
        # initializations

        self.structure.add_site_property("velocities", velocities.tolist())


class Incar(dict, MSONable):
    """
    INCAR object for reading and writing INCAR files. Essentially consists of
    a dictionary with some helper functions
    """

    def __init__(self, params=None):
        """
        Creates an Incar object.

        Args:
            params (dict): A set of input parameters as a dictionary.
        """
        super(Incar, self).__init__()
        if params:

            # if Incar contains vector-like magmoms given as a list
            # of floats, convert to a list of lists
            if (params.get("MAGMOM") and isinstance(params["MAGMOM"][0], (int, float))) \
                    and (params.get("LSORBIT") or params.get("LNONCOLLINEAR")):
                val = []
                for i in range(len(params["MAGMOM"])//3):
                    val.append(params["MAGMOM"][i*3:(i+1)*3])
                params["MAGMOM"] = val

            self.update(params)

    def __setitem__(self, key, val):
        """
        Add parameter-val pair to Incar.  Warns if parameter is not in list of
        valid INCAR tags. Also cleans the parameter and val by stripping
        leading and trailing white spaces.
        """
        super(Incar, self).__setitem__(
            key.strip(), Incar.proc_val(key.strip(), val.strip())
            if isinstance(val, six.string_types) else val)

    def as_dict(self):
        d = dict(self)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        if d.get("MAGMOM") and isinstance(d["MAGMOM"][0], dict):
            d["MAGMOM"] = [Magmom.from_dict(m) for m in d["MAGMOM"]]
        return Incar({k: v for k, v in d.items() if k not in ("@module",
                                                              "@class")})

    def get_string(self, sort_keys=False, pretty=False):
        """
        Returns a string representation of the INCAR.  The reason why this
        method is different from the __str__ method is to provide options for
        pretty printing.

        Args:
            sort_keys (bool): Set to True to sort the INCAR parameters
                alphabetically. Defaults to False.
            pretty (bool): Set to True for pretty aligned output. Defaults
                to False.
        """
        keys = self.keys()
        if sort_keys:
            keys = sorted(keys)
        lines = []
        for k in keys:
            if k == "MAGMOM" and isinstance(self[k], list):
                value = []

                if (isinstance(self[k][0], list) or isinstance(self[k][0], Magmom)) and \
                        (self.get("LSORBIT") or self.get("LNONCOLLINEAR")):
                    value.append(" ".join(str(i) for j in self[k] for i in j))
                elif self.get("LSORBIT") or self.get("LNONCOLLINEAR"):
                    for m, g in itertools.groupby(self[k]):
                        value.append("3*{}*{}".format(len(tuple(g)), m))
                else:
                    # float() to ensure backwards compatibility between
                    # float magmoms and Magmom objects
                    for m, g in itertools.groupby(self[k], lambda x: float(x)):
                        value.append("{}*{}".format(len(tuple(g)), m))

                lines.append([k, " ".join(value)])
            elif isinstance(self[k], list):
                lines.append([k, " ".join([str(i) for i in self[k]])])
            else:
                lines.append([k, self[k]])

        if pretty:
            return str(tabulate([[l[0], "=", l[1]] for l in lines],
                                tablefmt="plain"))
        else:
            return str_delimited(lines, None, " = ") + "\n"

    def __str__(self):
        return self.get_string(sort_keys=True, pretty=False)

    def write_file(self, filename):
        """
        Write Incar to a file.

        Args:
            filename (str): filename to write to.
        """
        with zopen(filename, "wt") as f:
            f.write(self.__str__())

    @staticmethod
    def from_file(filename):
        """
        Reads an Incar object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            Incar object
        """
        with zopen(filename, "rt") as f:
            return Incar.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads an Incar object from a string.

        Args:
            string (str): Incar string

        Returns:
            Incar object
        """
        lines = list(clean_lines(string.splitlines()))
        params = {}
        for line in lines:
            for sline in line.split(';'):
                m = re.match(r'(\w+)\s*=\s*(.*)', sline.strip())
                if m:
                    key = m.group(1).strip()
                    val = m.group(2).strip()
                    val = Incar.proc_val(key, val)
                    params[key] = val
        return Incar(params)

    @staticmethod
    def proc_val(key, val):
        """
        Static helper method to convert INCAR parameters to proper types, e.g.,
        integers, floats, lists, etc.

        Args:
            key: INCAR parameter key
            val: Actual value of INCAR parameter.
        """
        list_keys = ("LDAUU", "LDAUL", "LDAUJ", "MAGMOM", "DIPOL",
                     "LANGEVIN_GAMMA", "QUAD_EFG", "EINT")
        bool_keys = ("LDAU", "LWAVE", "LSCALU", "LCHARG", "LPLANE", "LUSE_VDW",
                     "LHFCALC", "ADDGRID", "LSORBIT", "LNONCOLLINEAR")
        float_keys = ("EDIFF", "SIGMA", "TIME", "ENCUTFOCK", "HFSCREEN",
                      "POTIM", "EDIFFG", "AGGAC", "PARAM1", "PARAM2")
        int_keys = ("NSW", "NBANDS", "NELMIN", "ISIF", "IBRION", "ISPIN",
                    "ICHARG", "NELM", "ISMEAR", "NPAR", "LDAUPRINT", "LMAXMIX",
                    "ENCUT", "NSIM", "NKRED", "NUPDOWN", "ISPIND", "LDAUTYPE",
                    "IVDW")

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)

        try:
            if key in list_keys:
                output = []
                toks = re.findall(
                    r"(-?\d+\.?\d*)\*?(-?\d+\.?\d*)?\*?(-?\d+\.?\d*)?", val)
                for tok in toks:
                    if tok[2] and "3" in tok[0]:
                        output.extend(
                            [smart_int_or_float(tok[2])] * int(tok[0])
                            * int(tok[1]))
                    elif tok[1]:
                        output.extend([smart_int_or_float(tok[1])] *
                                      int(tok[0]))
                    else:
                        output.append(smart_int_or_float(tok[0]))
                return output
            if key in bool_keys:
                m = re.match(r"^\.?([T|F|t|f])[A-Za-z]*\.?", val)
                if m:
                    if m.group(1) == "T" or m.group(1) == "t":
                        return True
                    else:
                        return False
                raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*[e|E]?-?\d*", val).group(0))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))

        except ValueError:
            pass

        # Not in standard keys. We will try a hierarchy of conversions.
        try:
            val = int(val)
            return val
        except ValueError:
            pass

        try:
            val = float(val)
            return val
        except ValueError:
            pass

        if "true" in val.lower():
            return True

        if "false" in val.lower():
            return False

        return val.strip().capitalize()

    def diff(self, other):
        """
        Diff function for Incar.  Compares two Incars and indicates which
        parameters are the same and which are not. Useful for checking whether
        two runs were done using the same parameters.

        Args:
            other (Incar): The other Incar object to compare to.

        Returns:
            Dict of the following format:
            {"Same" : parameters_that_are_the_same,
            "Different": parameters_that_are_different}
            Note that the parameters are return as full dictionaries of values.
            E.g. {"ISIF":3}
        """
        similar_param = {}
        different_param = {}
        for k1, v1 in self.items():
            if k1 not in other:
                different_param[k1] = {"INCAR1": v1, "INCAR2": None}
            elif v1 != other[k1]:
                different_param[k1] = {"INCAR1": v1, "INCAR2": other[k1]}
            else:
                similar_param[k1] = v1
        for k2, v2 in other.items():
            if k2 not in similar_param and k2 not in different_param:
                if k2 not in self:
                    different_param[k2] = {"INCAR1": None, "INCAR2": v2}
        return {"Same": similar_param, "Different": different_param}

    def __add__(self, other):
        """
        Add all the values of another INCAR object to this object.
        Facilitates the use of "standard" INCARs.
        """
        params = {k: v for k, v in self.items()}
        for k, v in other.items():
            if k in self and v != self[k]:
                raise ValueError("Incars have conflicting values!")
            else:
                params[k] = v
        return Incar(params)


class Kpoints_supported_modes(Enum):
    Automatic = 0
    Gamma = 1
    Monkhorst = 2
    Line_mode = 3
    Cartesian = 4
    Reciprocal = 5

    def __str__(self):
        return self.name

    @staticmethod
    def from_string(s):
        c = s.lower()[0]
        for m in Kpoints_supported_modes:
            if m.name.lower()[0] == c:
                return m
        raise ValueError("Can't interprete Kpoint mode %s" % s)


class Kpoints(MSONable):
    """
    KPOINT reader/writer.
    """
    supported_modes = Kpoints_supported_modes

    def __init__(self, comment="Default gamma", num_kpts=0,
                 style=supported_modes.Gamma,
                 kpts=((1, 1, 1),), kpts_shift=(0, 0, 0),
                 kpts_weights=None, coord_type=None, labels=None,
                 tet_number=0, tet_weight=0, tet_connections=None):
        """
        Highly flexible constructor for Kpoints object.  The flexibility comes
        at the cost of usability and in general, it is recommended that you use
        the default constructor only if you know exactly what you are doing and
        requires the flexibility.  For most usage cases, the three automatic
        schemes can be constructed far more easily using the convenience static
        constructors (automatic, gamma_automatic, monkhorst_automatic) and it
        is recommended that you use those.

        Args:
            comment (str): String comment for Kpoints
            num_kpts: Following VASP method of defining the KPOINTS file, this
                parameter is the number of kpoints specified. If set to 0
                (or negative), VASP automatically generates the KPOINTS.
            style: Style for generating KPOINTS.  Use one of the
                Kpoints.supported_modes enum types.
            kpts (2D array): 2D array of kpoints.  Even when only a single
                specification is required, e.g. in the automatic scheme,
                the kpts should still be specified as a 2D array. e.g.,
                [[20]] or [[2,2,2]].
            kpts_shift (3x1 array): Shift for Kpoints.
            kpts_weights: Optional weights for kpoints.  Weights should be
                integers. For explicit kpoints.
            coord_type: In line-mode, this variable specifies whether the
                Kpoints were given in Cartesian or Reciprocal coordinates.
            labels: In line-mode, this should provide a list of labels for
                each kpt. It is optional in explicit kpoint mode as comments for
                k-points.
            tet_number: For explicit kpoints, specifies the number of
                tetrahedrons for the tetrahedron method.
            tet_weight: For explicit kpoints, specifies the weight for each
                tetrahedron for the tetrahedron method.
            tet_connections: For explicit kpoints, specifies the connections
                of the tetrahedrons for the tetrahedron method.
                Format is a list of tuples, [ (sym_weight, [tet_vertices]),
                ...]

        The default behavior of the constructor is for a Gamma centered,
        1x1x1 KPOINTS with no shift.
        """
        if num_kpts > 0 and (not labels) and (not kpts_weights):
            raise ValueError("For explicit or line-mode kpoints, either the "
                             "labels or kpts_weights must be specified.")

        self.comment = comment
        self.num_kpts = num_kpts
        self.kpts = kpts
        self.style = style
        self.coord_type = coord_type
        self.kpts_weights = kpts_weights
        self.kpts_shift = kpts_shift
        self.labels = labels
        self.tet_number = tet_number
        self.tet_weight = tet_weight
        self.tet_connections = tet_connections

    @property
    def style(self):
        return self._style

    @style.setter
    def style(self, style):
        if isinstance(style, six.string_types):
            style = Kpoints.supported_modes.from_string(style)

        if style in (Kpoints.supported_modes.Automatic,
                     Kpoints.supported_modes.Gamma,
                     Kpoints.supported_modes.Monkhorst) and len(self.kpts) > 1:
            raise ValueError("For fully automatic or automatic gamma or monk "
                             "kpoints, only a single line for the number of "
                             "divisions is allowed.")

        self._style = style

    @staticmethod
    def automatic(subdivisions):
        """
        Convenient static constructor for a fully automatic Kpoint grid, with
        gamma centered Monkhorst-Pack grids and the number of subdivisions
        along each reciprocal lattice vector determined by the scheme in the
        VASP manual.

        Args:
            subdivisions: Parameter determining number of subdivisions along
                each reciprocal lattice vector.

        Returns:
            Kpoints object
        """
        return Kpoints("Fully automatic kpoint scheme", 0,
                       style=Kpoints.supported_modes.Automatic,
                       kpts=[[subdivisions]])

    @staticmethod
    def gamma_automatic(kpts=(1, 1, 1), shift=(0, 0, 0)):
        """
        Convenient static constructor for an automatic Gamma centered Kpoint
        grid.

        Args:
            kpts: Subdivisions N_1, N_2 and N_3 along reciprocal lattice
                vectors. Defaults to (1,1,1)
            shift: Shift to be applied to the kpoints. Defaults to (0,0,0).

        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0,
                       Kpoints.supported_modes.Gamma, kpts=[kpts],
                       kpts_shift=shift)

    @staticmethod
    def monkhorst_automatic(kpts=(2, 2, 2), shift=(0, 0, 0)):
        """
        Convenient static constructor for an automatic Monkhorst pack Kpoint
        grid.

        Args:
            kpts: Subdivisions N_1, N_2 and N_3 along reciprocal lattice
                vectors. Defaults to (2,2,2)
            shift: Shift to be applied to the kpoints. Defaults to (0,0,0).

        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0,
                       Kpoints.supported_modes.Monkhorst, kpts=[kpts],
                       kpts_shift=shift)

    @staticmethod
    def automatic_density(structure, kppa, force_gamma=False):
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
        density. Uses Gamma centered meshes for hexagonal cells and
        Monkhorst-Pack grids otherwise.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.

        Args:
            structure (Structure): Input structure
            kppa (int): Grid density
            force_gamma (bool): Force a gamma centered mesh (default is to
                use gamma only for hexagonal cells or odd meshes)

        Returns:
            Kpoints
        """
        comment = "pymatgen 4.7.6+ generated KPOINTS with grid density = " + \
            "%.0f / atom" % kppa
        if math.fabs((math.floor(kppa ** (1 / 3) + 0.5)) ** 3 - kppa) < 1:
            kppa += kppa * 0.01
        latt = structure.lattice
        lengths = latt.abc
        ngrid = kppa / structure.num_sites
        mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3)

        num_div = [int(math.floor(max(mult / l, 1))) for l in lengths]

        is_hexagonal = latt.is_hexagonal()

        has_odd = any([i % 2 == 1 for i in num_div])
        if has_odd or is_hexagonal or force_gamma:
            style = Kpoints.supported_modes.Gamma
        else:
            style = Kpoints.supported_modes.Monkhorst

        return Kpoints(comment, 0, style, [num_div], [0, 0, 0])

    @staticmethod
    def automatic_gamma_density(structure, kppa):
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
        density. Uses Gamma centered meshes always. For GW.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.

        Args:
            structure:
                Input structure
            kppa:
                Grid density
        """

        latt = structure.lattice
        lengths = latt.abc
        ngrid = kppa / structure.num_sites

        mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3)
        num_div = [int(round(mult / l)) for l in lengths]

        # ensure that numDiv[i] > 0
        num_div = [i if i > 0 else 1 for i in num_div]

        # VASP documentation recommends to use even grids for n <= 8 and odd
        # grids for n > 8.
        num_div = [i + i % 2 if i <= 8 else i - i % 2 + 1 for i in num_div]

        style = Kpoints.supported_modes.Gamma

        comment = "pymatgen 4.7.6+ generated KPOINTS with grid density = " + \
                  "{} / atom".format(kppa)
        num_kpts = 0
        return Kpoints(comment, num_kpts, style, [num_div], [0, 0, 0])

    @staticmethod
    def automatic_density_by_vol(structure, kppvol, force_gamma=False):
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
        density per inverse Angstrom of reciprocal cell.

        Algorithm:
            Same as automatic_density()

        Args:
            structure (Structure): Input structure
            kppvol (int): Grid density per Angstrom^(-3) of reciprocal cell
            force_gamma (bool): Force a gamma centered mesh

        Returns:
            Kpoints
        """
        vol = structure.lattice.reciprocal_lattice.volume
        kppa = kppvol * vol * structure.num_sites
        return Kpoints.automatic_density(structure, kppa,
                                         force_gamma=force_gamma)

    @staticmethod
    def automatic_linemode(divisions, ibz):
        """
        Convenient static constructor for a KPOINTS in mode line_mode.
        gamma centered Monkhorst-Pack grids and the number of subdivisions
        along each reciprocal lattice vector determined by the scheme in the
        VASP manual.

        Args:
            divisions: Parameter determining the number of k-points along each
                hight symetry lines.
            ibz: HighSymmKpath object (pymatgen.symmetry.bandstructure)

        Returns:
            Kpoints object
        """
        kpoints = list()
        labels = list()
        for path in ibz.kpath["path"]:
            kpoints.append(ibz.kpath["kpoints"][path[0]])
            labels.append(path[0])
            for i in range(1, len(path) - 1):
                kpoints.append(ibz.kpath["kpoints"][path[i]])
                labels.append(path[i])
                kpoints.append(ibz.kpath["kpoints"][path[i]])
                labels.append(path[i])

            kpoints.append(ibz.kpath["kpoints"][path[-1]])
            labels.append(path[-1])

        return Kpoints("Line_mode KPOINTS file",
                       style=Kpoints.supported_modes.Line_mode,
                       coord_type="Reciprocal",
                       kpts=kpoints,
                       labels=labels,
                       num_kpts=int(divisions))

    @staticmethod
    def from_file(filename):
        """
        Reads a Kpoints object from a KPOINTS file.

        Args:
            filename (str): filename to read from.

        Returns:
            Kpoints object
        """
        with zopen(filename, "rt") as f:
            return Kpoints.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads a Kpoints object from a KPOINTS string.

        Args:
            string (str): KPOINTS string.

        Returns:
            Kpoints object
        """
        lines = [line.strip() for line in string.splitlines()]

        comment = lines[0]
        num_kpts = int(lines[1].split()[0].strip())
        style = lines[2].lower()[0]

        # Fully automatic KPOINTS
        if style == "a":
            return Kpoints.automatic(int(lines[3]))

        coord_pattern = re.compile(r'^\s*([\d+.\-Ee]+)\s+([\d+.\-Ee]+)\s+'
                                   r'([\d+.\-Ee]+)')

        # Automatic gamma and Monk KPOINTS, with optional shift
        if style == "g" or style == "m":
            kpts = [int(i) for i in lines[3].split()]
            kpts_shift = (0, 0, 0)
            if len(lines) > 4 and coord_pattern.match(lines[4]):
                try:
                    kpts_shift = [float(i) for i in lines[4].split()]
                except ValueError:
                    pass
            return Kpoints.gamma_automatic(kpts, kpts_shift) if style == "g" \
                else Kpoints.monkhorst_automatic(kpts, kpts_shift)

        # Automatic kpoints with basis
        if num_kpts <= 0:
            style = Kpoints.supported_modes.Cartesian if style in "ck" \
                else Kpoints.supported_modes.Reciprocal
            kpts = [[float(j) for j in lines[i].split()] for i in range(3, 6)]
            kpts_shift = [float(i) for i in lines[6].split()]
            return Kpoints(comment=comment, num_kpts=num_kpts, style=style,
                           kpts=kpts, kpts_shift=kpts_shift)

        # Line-mode KPOINTS, usually used with band structures
        if style == "l":
            coord_type = "Cartesian" if lines[3].lower()[0] in "ck" \
                else "Reciprocal"
            style = Kpoints.supported_modes.Line_mode
            kpts = []
            labels = []
            patt = re.compile(r'([e0-9.\-]+)\s+([e0-9.\-]+)\s+([e0-9.\-]+)'
                              r'\s*!*\s*(.*)')
            for i in range(4, len(lines)):
                line = lines[i]
                m = patt.match(line)
                if m:
                    kpts.append([float(m.group(1)), float(m.group(2)),
                                 float(m.group(3))])
                    labels.append(m.group(4).strip())
            return Kpoints(comment=comment, num_kpts=num_kpts, style=style,
                           kpts=kpts, coord_type=coord_type, labels=labels)

        # Assume explicit KPOINTS if all else fails.
        style = Kpoints.supported_modes.Cartesian if style in "ck" \
            else Kpoints.supported_modes.Reciprocal
        kpts = []
        kpts_weights = []
        labels = []
        tet_number = 0
        tet_weight = 0
        tet_connections = None

        for i in range(3, 3 + num_kpts):
            toks = lines[i].split()
            kpts.append([float(j) for j in toks[0:3]])
            kpts_weights.append(float(toks[3]))
            if len(toks) > 4:
                labels.append(toks[4])
            else:
                labels.append(None)
        try:
            # Deal with tetrahedron method
            if lines[3 + num_kpts].strip().lower()[0] == "t":
                toks = lines[4 + num_kpts].split()
                tet_number = int(toks[0])
                tet_weight = float(toks[1])
                tet_connections = []
                for i in range(5 + num_kpts, 5 + num_kpts + tet_number):
                    toks = lines[i].split()
                    tet_connections.append((int(toks[0]),
                                            [int(toks[j])
                                             for j in range(1, 5)]))
        except IndexError:
            pass

        return Kpoints(comment=comment, num_kpts=num_kpts,
                       style=Kpoints.supported_modes[str(style)],
                       kpts=kpts, kpts_weights=kpts_weights,
                       tet_number=tet_number, tet_weight=tet_weight,
                       tet_connections=tet_connections, labels=labels)

    def write_file(self, filename):
        """
        Write Kpoints to a file.

        Args:
            filename (str): Filename to write to.
        """
        with zopen(filename, "wt") as f:
            f.write(self.__str__())

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        lines = [self.comment, str(self.num_kpts), self.style.name]
        style = self.style.name.lower()[0]
        if style == "l":
            lines.append(self.coord_type)
        for i in range(len(self.kpts)):
            lines.append(" ".join([str(x) for x in self.kpts[i]]))
            if style == "l":
                lines[-1] += " ! " + self.labels[i]
                if i % 2 == 1:
                    lines[-1] += "\n"
            elif self.num_kpts > 0:
                if self.labels is not None:
                    lines[-1] += " %i %s" % (self.kpts_weights[i],
                                             self.labels[i])
                else:
                    lines[-1] += " %i" % (self.kpts_weights[i])

        # Print tetrahedron parameters if the number of tetrahedrons > 0
        if style not in "lagm" and self.tet_number > 0:
            lines.append("Tetrahedron")
            lines.append("%d %f" % (self.tet_number, self.tet_weight))
            for sym_weight, vertices in self.tet_connections:
                lines.append("%d %d %d %d %d" % (sym_weight, vertices[0],
                                                 vertices[1], vertices[2],
                                                 vertices[3]))

        # Print shifts for automatic kpoints types if not zero.
        if self.num_kpts <= 0 and tuple(self.kpts_shift) != (0, 0, 0):
            lines.append(" ".join([str(x) for x in self.kpts_shift]))
        return "\n".join(lines) + "\n"

    def as_dict(self):
        """json friendly dict representation of Kpoints"""
        d = {"comment": self.comment, "nkpoints": self.num_kpts,
             "generation_style": self.style.name, "kpoints": self.kpts,
             "usershift": self.kpts_shift,
             "kpts_weights": self.kpts_weights, "coord_type": self.coord_type,
             "labels": self.labels, "tet_number": self.tet_number,
             "tet_weight": self.tet_weight,
             "tet_connections": self.tet_connections}

        optional_paras = ["genvec1", "genvec2", "genvec3", "shift"]
        for para in optional_paras:
            if para in self.__dict__:
                d[para] = self.__dict__[para]
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        comment = d.get("comment", "")
        generation_style = d.get("generation_style")
        kpts = d.get("kpoints", [[1, 1, 1]])
        kpts_shift = d.get("usershift", [0, 0, 0])
        num_kpts = d.get("nkpoints", 0)
        return cls(comment=comment, kpts=kpts, style=generation_style,
                   kpts_shift=kpts_shift, num_kpts=num_kpts,
                   kpts_weights=d.get("kpts_weights"),
                   coord_type=d.get("coord_type"),
                   labels=d.get("labels"), tet_number=d.get("tet_number", 0),
                   tet_weight=d.get("tet_weight", 0),
                   tet_connections=d.get("tet_connections"))


def parse_string(s):
    return "{}".format(s.strip())


def parse_bool(s):
    m = re.match(r"^\.?([TFtf])[A-Za-z]*\.?", s)
    if m:
        if m.group(1) == "T" or m.group(1) == "t":
            return True
        else:
            return False
    raise ValueError(s + " should be a boolean type!")


def parse_float(s):
    return float(re.search(r"^-?\d*\.?\d*[eE]?-?\d*", s).group(0))


def parse_int(s):
    return int(re.match(r"^-?[0-9]+", s).group(0))


def parse_list(s):
    return [float(y) for y in re.split(r"\s+", s.strip()) if not y.isalpha()]


Orbital = namedtuple('Orbital', ['n', 'l', 'j', 'E', 'occ'])
OrbitalDescription = namedtuple('OrbitalDescription',
                                ['l', 'E', 'Type', "Rcut", "Type2", "Rcut2"])


class PotcarSingle(object):
    """
    Object for a **single** POTCAR. The builder assumes the complete string is
    the POTCAR contains the complete untouched data in "data" as a string and
    a dict of keywords.

    Args:
        data:
            Complete and single potcar file as a string.

    .. attribute:: data

        POTCAR data as a string.

    .. attribute:: keywords

        Keywords parsed from the POTCAR as a dict. All keywords are also
        accessible as attributes in themselves. E.g., potcar.enmax,
        potcar.encut, etc.
    """
    functional_dir = {"PBE": "POT_GGA_PAW_PBE",
                      "PBE_52": "POT_GGA_PAW_PBE_52",
                      "PBE_54": "POT_GGA_PAW_PBE_54",
                      "LDA": "POT_LDA_PAW",
                      "LDA_52": "POT_LDA_PAW_52",
                      "LDA_54": "POT_LDA_PAW_54",
                      "PW91": "POT_GGA_PAW_PW91",
                      "LDA_US": "POT_LDA_US",
                      "PW91_US": "POT_GGA_US_PW91"}

    functional_tags = {"pe": {"name": "PBE", "class": "GGA"},
                       "91": {"name": "PW91", "class": "GGA"},
                       "rp": {"name": "revPBE", "class": "GGA"},
                       "am": {"name": "AM05", "class": "GGA"},
                       "ps": {"name": "PBEsol", "class": "GGA"},
                       "pw": {"name": "PW86", "class": "GGA"},
                       "lm": {"name": "Langreth-Mehl-Hu", "class": "GGA"},
                       "pb": {"name": "Perdew-Becke", "class": "GGA"},
                       "ca": {"name": "Perdew-Zunger81", "class": "LDA"},
                       "hl": {"name": "Hedin-Lundquist", "class": "LDA"},
                       "wi": {"name": "Wigner Interpoloation", "class": "LDA"}}

    parse_functions = {"LULTRA": parse_bool,
                       "LUNSCR": parse_bool,
                       "LCOR": parse_bool,
                       "LPAW": parse_bool,
                       "EATOM": parse_float,
                       "RPACOR": parse_float,
                       "POMASS": parse_float,
                       "ZVAL": parse_float,
                       "RCORE": parse_float,
                       "RWIGS": parse_float,
                       "ENMAX": parse_float,
                       "ENMIN": parse_float,
                       "EMMIN": parse_float,
                       "EAUG": parse_float,
                       "DEXC": parse_float,
                       "RMAX": parse_float,
                       "RAUG": parse_float,
                       "RDEP": parse_float,
                       "RDEPT": parse_float,
                       "QCUT": parse_float,
                       "QGAM": parse_float,
                       "RCLOC": parse_float,
                       "IUNSCR": parse_int,
                       "ICORE": parse_int,
                       "NDATA": parse_int,
                       "VRHFIN": parse_string,
                       "LEXCH": parse_string,
                       "TITEL": parse_string,
                       "STEP": parse_list,
                       "RRKJ": parse_list,
                       "GGA": parse_list}

    def __init__(self, data):
        self.data = data  # raw POTCAR as a string

        # Vasp parses header in vasprun.xml and this differs from the titel
        self.header = data.split("\n")[0].strip()

        search_lines = re.search(r"(?s)(parameters from PSCTR are:"
                                 r".*?END of PSCTR-controll parameters)",
                                 data).group(1)

        self.keywords = {}
        for key, val in re.findall(r"(\S+)\s*=\s*(.*?)(?=;|$)",
                                   search_lines, flags=re.MULTILINE):
            try:
                self.keywords[key] = self.parse_functions[key](val)
            except KeyError:
                warnings.warn("Ignoring unknown variable type %s" % key)

        PSCTR = OrderedDict()

        array_search = re.compile(r"(-*[0-9.]+)")
        orbitals = []
        descriptions = []
        atomic_configuration = re.search(r"Atomic configuration\s*\n?"
                                         r"(.*?)Description", search_lines)
        if atomic_configuration:
            lines = atomic_configuration.group(1).splitlines()
            num_entries = re.search(r"([0-9]+)", lines[0]).group(1)
            num_entries = int(num_entries)
            PSCTR['nentries'] = num_entries
            for line in lines[1:]:
                orbit = array_search.findall(line)
                if orbit:
                    orbitals.append(self.Orbital(int(orbit[0]),
                                                 int(orbit[1]),
                                                 float(orbit[2]),
                                                 float(orbit[3]),
                                                 float(orbit[4])))
            PSCTR['Orbitals'] = tuple(orbitals)

        description_string = re.search(r"(?s)Description\s*\n"
                                       r"(.*?)Error from kinetic"
                                       r" energy argument \(eV\)",
                                       search_lines)
        if description_string:
            for line in description_string.group(1).splitlines():
                description = array_search.findall(line)
                if description:
                    descriptions.append(
                        OrbitalDescription(
                            int(description[0]), float(description[1]),
                            int(description[2]), float(description[3]),
                            int(description[4]) if len(description) > 4 else None,
                            float(description[5]) if len(description) > 4 else None))

        if descriptions:
            PSCTR['OrbitalDescriptions'] = tuple(descriptions)

        rrkj_kinetic_energy_string = re.search(
            r"(?s)Error from kinetic energy argument \(eV\)\s*\n"
            r"(.*?)END of PSCTR-controll parameters",
            search_lines)
        rrkj_array = []
        if rrkj_kinetic_energy_string:
            for line in rrkj_kinetic_energy_string.group(1).splitlines():
                if "=" not in line:
                    rrkj_array += parse_list(line.strip('\n'))
            if rrkj_array:
                PSCTR['RRKJ'] = tuple(rrkj_array)

        PSCTR.update(self.keywords)
        self.PSCTR = OrderedDict(sorted(PSCTR.items(), key=lambda x: x[0]))
        self.hash = self.get_potcar_hash()

    def __str__(self):
        return self.data + "\n"

    @property
    def electron_configuration(self):
        el = Element.from_Z(self.atomic_no)
        full_config = el.full_electronic_structure
        nelect = self.nelectrons
        config = []
        while nelect > 0:
            e = full_config.pop(-1)
            config.append(e)
            nelect -= e[-1]
        return config

    def write_file(self, filename):
        with zopen(filename, "wt") as f:
            f.write(self.__str__())

    @staticmethod
    def from_file(filename):
        try:
            with zopen(filename, "rt") as f:
                return PotcarSingle(f.read())
        except UnicodeDecodeError:
            warnings.warn("POTCAR contains invalid unicode errors. "
                          "We will attempt to read it by ignoring errors.")
            import codecs
            with codecs.open(filename, "r", encoding="utf-8",
                             errors="ignore") as f:
                return PotcarSingle(f.read())

    @staticmethod
    def from_symbol_and_functional(symbol, functional=None):
        if functional is None:
            functional = SETTINGS.get("PMG_DEFAULT_FUNCTIONAL", "PBE")
        funcdir = PotcarSingle.functional_dir[functional]
        d = SETTINGS.get("PMG_VASP_PSP_DIR")
        if d is None:
            raise ValueError(
                "No POTCAR for %s with functional %s found. "
                "Please set the PMG_VASP_PSP_DIR environment in "
                ".pmgrc.yaml, or you may need to set "
                "PMG_DEFAULT_FUNCTIONAL to PBE_52 or PBE_54 if you "
                "are using newer psps from VASP." % (symbol, functional))
        paths_to_try = [os.path.join(d, funcdir, "POTCAR.{}".format(symbol)),
                        os.path.join(d, funcdir, symbol, "POTCAR")]
        for p in paths_to_try:
            p = os.path.expanduser(p)
            p = zpath(p)
            if os.path.exists(p):
                return PotcarSingle.from_file(p)
        raise IOError("You do not have the right POTCAR with functional " +
                      "{} and label {} in your VASP_PSP_DIR".format(functional,
                                                                    symbol))

    @property
    def symbol(self):
        """
        Symbol of POTCAR, e.g., Fe_pv
        """
        return self.keywords["TITEL"].split(" ")[1].strip()

    @property
    def element(self):
        """
        Attempt to return the atomic symbol based on the VRHFIN keyword.
        """
        element = self.keywords["VRHFIN"].split(":")[0].strip()
        try:
            return Element(element).symbol
        except ValueError:
            # VASP incorrectly gives the element symbol for Xe as "X"
            # Some potentials, e.g., Zr_sv, gives the symbol as r.
            if element == "X":
                return "Xe"
            return Element(self.symbol.split("_")[0]).symbol

    @property
    def atomic_no(self):
        """
        Attempt to return the atomic number based on the VRHFIN keyword.
        """
        return Element(self.element).Z

    @property
    def nelectrons(self):
        return self.zval

    @property
    def potential_type(self):
        if self.lultra:
            return "US"
        elif self.lpaw:
            return "PAW"
        else:
            return "NC"

    @property
    def functional(self):
        return self.functional_tags.get(self.LEXCH.lower(), {}).get('name')

    @property
    def functional_class(self):
        return self.functional_tags.get(self.LEXCH.lower(), {}).get('class')

    def get_potcar_hash(self):
        hash_str = ""
        for k, v in self.PSCTR.items():
            hash_str += "{}".format(k)
            if isinstance(v, int):
                hash_str += "{}".format(v)
            elif isinstance(v, float):
                hash_str += "{:.3f}".format(v)
            elif isinstance(v, bool):
                hash_str += "{}".format(bool)
            elif isinstance(v, (tuple, list)):
                for item in v:
                    if isinstance(item, float):
                        hash_str += "{:.3f}".format(item)
                    elif isinstance(item, (Orbital, OrbitalDescription)):
                        for item_v in item:
                            if isinstance(item_v, (int, str)):
                                hash_str += "{}".format(item_v)
                            elif isinstance(item_v, float):
                                hash_str += "{:.3f}".format(item_v)
                            else:
                                hash_str += "{}".format(item_v) if item_v else ""
            else:
                hash_str += v.replace(" ", "")

        self.hash_str = hash_str
        return md5(hash_str.lower().encode('utf-8')).hexdigest()

    def __getattr__(self, a):
        """
        Delegates attributes to keywords. For example, you can use
        potcarsingle.enmax to get the ENMAX of the POTCAR.

        For float type properties, they are converted to the correct float. By
        default, all energies in eV and all length scales are in Angstroms.
        """

        try:
            return self.keywords[a.upper()]
        except:
            raise AttributeError(a)


class Potcar(list, MSONable):
    """
    Object for reading and writing POTCAR files for calculations. Consists of a
    list of PotcarSingle.

    Args:
        symbols ([str]): Element symbols for POTCAR. This should correspond
            to the symbols used by VASP. E.g., "Mg", "Fe_pv", etc.
        functional (str): Functional used. To know what functional options
            there are, use Potcar.FUNCTIONAL_CHOICES. Note that VASP has
            different versions of the same functional. By default, the old
            PBE functional is used. If you want the newer ones, use PBE_52 or
            PBE_54. Note that if you intend to compare your results with the
            Materials Project, you should use the default setting. You can also
            override the default by setting PMG_DEFAULT_FUNCTIONAL in your
            .pmgrc.yaml.
        sym_potcar_map (dict): Allows a user to specify a specific element
            symbol to raw POTCAR mapping.
    """

    FUNCTIONAL_CHOICES = list(PotcarSingle.functional_dir.keys())

    def __init__(self, symbols=None, functional=None, sym_potcar_map=None):
        if functional is None:
            functional = SETTINGS.get("PMG_DEFAULT_FUNCTIONAL", "PBE")
        super(Potcar, self).__init__()
        self.functional = functional
        if symbols is not None:
            self.set_symbols(symbols, functional, sym_potcar_map)

    def as_dict(self):
        return {"functional": self.functional, "symbols": self.symbols,
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @classmethod
    def from_dict(cls, d):
        return Potcar(symbols=d["symbols"], functional=d["functional"])

    @staticmethod
    def from_file(filename):
        try:
            with zopen(filename, "rt") as f:
                fdata = f.read()
        except UnicodeDecodeError:
            warnings.warn("POTCAR contains invalid unicode errors. "
                          "We will attempt to read it by ignoring errors.")
            import codecs
            with codecs.open(filename, "r", encoding="utf-8",
                             errors="ignore") as f:
                fdata = f.read()

        potcar = Potcar()
        potcar_strings = re.compile(r"\n?(\s*.*?End of Dataset)",
                                    re.S).findall(fdata)
        functionals = []
        for p in potcar_strings:
            single = PotcarSingle(p)
            potcar.append(single)
            functionals.append(single.functional)
        if len(set(functionals)) != 1:
            raise ValueError("File contains incompatible functionals!")
        else:
            potcar.functional = functionals[0]
        return potcar

    def __str__(self):
        return "\n".join([str(potcar).strip("\n") for potcar in self]) + "\n"

    def write_file(self, filename):
        """
        Write Potcar to a file.

        Args:
            filename (str): filename to write to.
        """
        with zopen(filename, "wt") as f:
            f.write(self.__str__())

    @property
    def symbols(self):
        """
        Get the atomic symbols of all the atoms in the POTCAR file.
        """
        return [p.symbol for p in self]

    @symbols.setter
    def symbols(self, symbols):
        self.set_symbols(symbols, functional=self.functional)

    @property
    def spec(self):
        """
        Get the atomic symbols and hash of all the atoms in the POTCAR file.
        """
        return [{"symbol": p.symbol, "hash": p.get_potcar_hash()} for p in self]

    def set_symbols(self, symbols, functional=None,
                    sym_potcar_map=None):
        """
        Initialize the POTCAR from a set of symbols. Currently, the POTCARs can
        be fetched from a location specified in .pmgrc.yaml. Use pmg config
        to add this setting.

        Args:
            symbols ([str]): A list of element symbols
            functional (str): The functional to use. If None, the setting
                PMG_DEFAULT_FUNCTIONAL in .pmgrc.yaml is used, or if this is
                not set, it will default to PBE.
            sym_potcar_map (dict): A map of symbol:raw POTCAR string. If
                sym_potcar_map is specified, POTCARs will be generated from
                the given map data rather than the config file location.
        """
        del self[:]
        if sym_potcar_map:
            for el in symbols:
                self.append(PotcarSingle(sym_potcar_map[el]))
        else:
            for el in symbols:
                p = PotcarSingle.from_symbol_and_functional(el, functional)
                self.append(p)


class VaspInput(dict, MSONable):
    """
    Class to contain a set of vasp input objects corresponding to a run.

    Args:
        incar: Incar object.
        kpoints: Kpoints object.
        poscar: Poscar object.
        potcar: Potcar object.
        optional_files: Other input files supplied as a dict of {
            filename: object}. The object should follow standard pymatgen
            conventions in implementing a as_dict() and from_dict method.
    """

    def __init__(self, incar, kpoints, poscar, potcar, optional_files=None,
                 **kwargs):
        super(VaspInput, self).__init__(**kwargs)
        self.update({'INCAR': incar,
                     'KPOINTS': kpoints,
                     'POSCAR': poscar,
                     'POTCAR': potcar})
        if optional_files is not None:
            self.update(optional_files)

    def __str__(self):
        output = []
        for k, v in self.items():
            output.append(k)
            output.append(str(v))
            output.append("")
        return "\n".join(output)

    def as_dict(self):
        d = {k: v.as_dict() for k, v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        sub_d = {"optional_files": {}}
        for k, v in d.items():
            if k in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
                sub_d[k.lower()] = dec.process_decoded(v)
            elif k not in ["@module", "@class"]:
                sub_d["optional_files"][k] = dec.process_decoded(v)
        return cls(**sub_d)

    def write_input(self, output_dir=".", make_dir_if_not_present=True):
        """
        Write VASP input to a directory.

        Args:
            output_dir (str): Directory to write to. Defaults to current
                directory (".").
            make_dir_if_not_present (bool): Create the directory if not
                present. Defaults to True.
        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.items():
            with zopen(os.path.join(output_dir, k), "wt") as f:
                f.write(v.__str__())

    @staticmethod
    def from_directory(input_dir, optional_files=None):
        """
        Read in a set of VASP input from a directory. Note that only the
        standard INCAR, POSCAR, POTCAR and KPOINTS files are read unless
        optional_filenames is specified.

        Args:
            input_dir (str): Directory to read VASP input from.
            optional_files (dict): Optional files to read in as well as a
                dict of {filename: Object type}. Object type must have a
                static method from_file.
        """
        sub_d = {}
        for fname, ftype in [("INCAR", Incar), ("KPOINTS", Kpoints),
                             ("POSCAR", Poscar), ("POTCAR", Potcar)]:
            fullzpath = zpath(os.path.join(input_dir, fname))
            sub_d[fname.lower()] = ftype.from_file(fullzpath)
        sub_d["optional_files"] = {}
        if optional_files is not None:
            for fname, ftype in optional_files.items():
                sub_d["optional_files"][fname] = \
                    ftype.from_file(os.path.join(input_dir, fname))
        return VaspInput(**sub_d)
