#!/usr/bin/env python

"""
Classes for reading/manipulating/writing VASP files. All major VASP input 
files, plus the more common Vasp output files are available.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Rickard Armiento, Vincent L Chevrier"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import os
import glob
import re
import math
import itertools
import warnings
import xml.sax.handler
import StringIO
from collections import defaultdict
import ConfigParser
import logging

import numpy as np
from numpy.linalg import det

from pymatgen.core.physical_constants import AMU_TO_KG, BOLTZMANN_CONST
from pymatgen.core.design_patterns import Enum
from pymatgen.io.io_abc import VaspInput
from pymatgen.util.string_utils import str_aligned, str_delimited
from pymatgen.util.io_utils import file_open_zip_aware, clean_lines, micro_pyawk, clean_json
from pymatgen.core.structure import Structure, Composition
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine, get_reconstructed_band_structure
from pymatgen.core.lattice import Lattice
import pymatgen


logger = logging.getLogger(__name__)


class Poscar(VaspInput):
    """
    Object for representing the data in a POSCAR or CONTCAR file.
    Please note that this current implementation. Most attributes can be set
    directly.
    
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
        
        Predictor corrector coordinates for each site (typically read in from a 
        MD CONTCAR).
        
    .. attribute:: temperature
    
        Temperature of velocity Maxwell-Boltzmann initialization. Initialized
        to -1 (MB hasn't been performed).
    """

    def __init__(self, structure, comment=None, selective_dynamics=None,
                 true_names=True, velocities=None, predictor_corrector=None):
        """        
        Args:
            structure:
                Structure object. See pymatgen.core.structure.Structure.
            comment:
                Optional comment line for POSCAR. Defaults to unit cell
                formula of structure. Defaults to None.
            selective_dynamics:
                Nx3 2D array of boolean values for selective dynamics, where N
                is number of sites. Defaults to None.
            true_names:
                Set to False is the names in the POSCAR are not well-defined
                and ambiguous. This situation arises commonly in vasp < 5 where
                the POSCAR sometimes does not contain element symbols. Defaults
                to True.
            velocities:
                Velocities for the POSCAR. Typically parsed in MD runs or can be
                used to initialize velocities.
            predictor_corrector:
                Velocities for the POSCAR. Typically parsed in MD runs or can be
                used to initialize velocities.
        """

        if structure.is_ordered:
            self.structure = structure
            self.true_names = true_names
            self.selective_dynamics = selective_dynamics
            self.comment = structure.formula if comment == None else comment
            self.velocities = velocities
            self.predictor_corrector = predictor_corrector
        else:
            raise ValueError("Structure with partial occupancies cannot be converted into POSCAR!")

        self.temperature = -1

    @property
    def struct(self):
        """
        .. deprecated:: 1.9.1
    
            For backwards compatibility.
        """
        return self.structure

    @property
    def site_symbols(self):
        """
        Sequence of symbols associated with the Poscar. Similar to 6th line in
        vasp 5+ POSCAR.
        """
        syms = [site.specie.symbol for site in self.structure]
        return [s for (s, data) in itertools.groupby(syms)]

    @property
    def natoms(self):
        """
        Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.
        """
        syms = [site.specie.symbol for site in self.structure]
        return [len(tuple(data)) for (s, data) in itertools.groupby(syms)]

    def __setattr__(self, name, value):
        if name in ("selective_dynamics", "velocities"):
            if value:
                dim = np.array(value).shape
                if dim[1] != 3 or dim[0] != len(self.structure):
                    raise ValueError("{} array must be same length as the structure.".format(name))
        elif name == "structure":
            #If we set a new structure, we should discard the velocities and
            #predictor_corrector and selective dynamics.
            self.velocities = None
            self.predictor_corrector = None
            self.selective_dynamics = None
        super(Poscar, self).__setattr__(name, value)

    @staticmethod
    def from_file(filename, check_for_POTCAR=True):
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
            filename:
                File name containing Poscar data.
            check_for_POTCAR:
                Whether to check if a POTCAR is present in the same directory
                as the POSCAR. Defaults to True.
                
        Returns:
            Poscar object.
        """

        dirname = os.path.dirname(os.path.abspath(filename))
        names = None
        if check_for_POTCAR:
            for f in os.listdir(dirname):
                if f == "POTCAR":
                    try:
                        potcar = Potcar.from_file(os.path.join(dirname, f))
                        names = [sym.split("_")[0] for sym in potcar.symbols]
                    except:
                        names = None
        with file_open_zip_aware(filename, "r") as f:
            return Poscar.from_string(f.read(), names)

    @staticmethod
    def from_string(data, default_names=None):
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
            data:
                string containing Poscar data.
            default_names:
                default symbols for the POSCAR file, usually coming from a
                POTCAR in the same directory.
                
        Returns:
            Poscar object.
        """

        chunks = re.split("^\s*$", data.strip(), flags=re.MULTILINE)

        #Parse positions
        lines = tuple(clean_lines(chunks[0].split("\n"), False))
        comment = lines[0]
        scale = float(lines[1])
        lattice = np.array([[float(s) for s in line.split()] for line in lines[2:5]])
        if scale < 0:
            # In vasp, a negative scale factor is treated as a volume. We need
            # to translate this to a proper lattice vector scaling.
            vol = abs(det(lattice))
            lattice = (-scale / vol) ** (1 / 3) * lattice
        else:
            lattice = scale * lattice
        lattice = Lattice(lattice)

        vasp5_symbols = False
        try:
            natoms = [int(s) for s in lines[5].split()]
            ipos = 6
        except:
            vasp5_symbols = True
            symbols = lines[5].split()
            natoms = [int(s) for s in lines[6].split()]
            atomic_symbols = list()
            for i in xrange(len(natoms)):
                atomic_symbols.extend([symbols[i]] * natoms[i])
            ipos = 7

        postype = lines[ipos].split()[0]

        sdynamics = False
        # Selective dynamics
        if postype[0] in 'sS':
            sdynamics = True
            ipos += 1
            postype = lines[ipos].split()[0]

        cart = postype[0] in 'cCkK'
        nsites = sum(natoms)

        #If default_names is specified (usually coming from a POTCAR), use them.
        #This is in line with Vasp's parsing order that the POTCAR specified 
        #is the default used.
        if default_names:
            try:
                atomic_symbols = list()
                for i in xrange(len(natoms)):
                    atomic_symbols.extend([default_names[i]] * natoms[i])
                vasp5_symbols = True
            except:
                pass
        if not vasp5_symbols:
            ind = 3 if not sdynamics else 6
            try:
                #check if names are appended at the end of the coordinates.
                atomic_symbols = [l.split()[ind] for l in lines[ipos + 1:ipos + 1 + nsites]]
                #Ensure symbols are valid elements
                if not all([Element.is_valid_symbol(sym) for sym in atomic_symbols]):
                    raise ValueError("Non-valid symbols detected.")
                vasp5_symbols = True
            except:
                #Defaulting to false names.
                atomic_symbols = list()
                for i in xrange(len(natoms)):
                    sym = Element.from_Z(i + 1).symbol
                    atomic_symbols.extend([sym] * natoms[i])
                warnings.warn("Elements in POSCAR cannot be determined. Defaulting to false names, " + " ".join(atomic_symbols) + ".")

        # read the atomic coordinates
        coords = []
        selective_dynamics = list() if sdynamics else None
        for i in xrange(nsites):
            toks = lines[ipos + 1 + i].split()
            coords.append([float(s) for s in toks[:3]])
            if sdynamics:
                selective_dynamics.append([True if tok.upper()[0] == 'T' else False for tok in toks[3:6]])

        struct = Structure(lattice, atomic_symbols, coords, False, False, cart)

        #parse velocities if any
        velocities = []
        if len(chunks) > 1:
            for line in chunks[1].strip().split("\n"):
                velocities.append([float(tok) for tok in re.split("\s+", line.strip())])

        predictor_corrector = []
        if len(chunks) > 2:
            lines = chunks[2].strip().split("\n")
            predictor_corrector.append([int(lines[0])])
            for line in lines[1:]:
                predictor_corrector.append([float(tok) for tok in re.split("\s+", line.strip())])

        return Poscar(struct, comment, selective_dynamics, vasp5_symbols,
                      velocities=velocities,
                      predictor_corrector=predictor_corrector)

    def get_string(self, direct=True, vasp4_compatible=False,
                   significant_figures=6):
        """
        Returns a string to be written as a POSCAR file. By default, site
        symbols are written, which means compatibility is for vasp >= 5.
        
        Args:
            direct:
                Whether coordinates are output in direct or cartesian. Defaults
                to True.
            vasp4_compatible:
                Set to True to omit site symbols on 6th line to maintain
                backward vasp 4.x compatibility. Defaults to False.
            significant_figures:
                Number of significant figures to output all quantities. Defaults
                to 6. Note that positions are output in fixed point, while
                velocities are output in scientific format.
                
        Returns:
            String representation of POSCAR.
        """
        lines = [self.comment]
        lines.append("1.0")
        lines.append(str(self.structure.lattice))
        if self.true_names and not vasp4_compatible:
            lines.append(" ".join(self.site_symbols))
        lines.append(" ".join([str(x) for x in self.natoms]))
        if self.selective_dynamics:
            lines.append("Selective dynamics")
        lines.append('direct' if direct else 'cartesian')

        format_str = "{{:.{0}f}}".format(significant_figures)

        for (i, site) in enumerate(self.structure):
            coords = site.frac_coords if direct else site.coords
            line = " ".join([format_str.format(c) for c in coords])
            if self.selective_dynamics:
                sd = ['T' if j else 'F' for j in self.selective_dynamics[i]]
                line += " %s %s %s" % (sd[0], sd[1], sd[2])
            line += " " + site.species_string
            lines.append(line)

        if self.velocities:
            lines.append("")
            for v in self.velocities:
                lines.append(" ".join([format_str.format(i) for i in v]))

        if self.predictor_corrector:
            lines.append("")
            lines.append(str(self.predictor_corrector[0][0]))
            lines.append(str(self.predictor_corrector[1][0]))
            for v in self.predictor_corrector[2:]:
                lines.append(" ".join([format_str.format(i) for i in v]))

        return "\n".join(lines)

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
        with open(filename, 'w') as f:
            f.write(self.get_string(**kwargs) + "\n")

    @property
    def to_dict(self):
        d = {}
        d['module'] = self.__class__.__module__
        d['class'] = self.__class__.__name__
        d['structure'] = self.structure.to_dict
        d['true_names'] = self.true_names
        d['selective_dynamics'] = self.selective_dynamics
        d['velocities'] = self.velocities
        d['predictor_corrector'] = self.predictor_corrector
        d['comment'] = self.comment
        return d

    @staticmethod
    def from_dict(d):
        return Poscar(Structure.from_dict(d['structure']), comment=d['comment'],
                      selective_dynamics=d['selective_dynamics'],
                      true_names=d['true_names'],
                      velocities=d.get('velocities', None),
                      predictor_corrector=d.get('predictor_corrector', None))

    def set_temperature(self, temperature):
        '''
        Initializes the velocities based on Maxwell-Boltzmann distribution.
        Removes linear, but not angular drift (same as VASP)
        
        Scales the energies to the exact temperature (microcanonical ensemble)
        Velocities are given in A/fs. This is the vasp default when
        direct/cartesian is not specified (even when positions are given in
        direct coordinates)
        
        Overwrites imported velocities, if any.
        
        Args:
            temperature:
                Temperature in Kelvin.
        '''
        # mean 0 variance 1
        velocities = np.random.randn(len(self.structure), 3)

        #in AMU, (N,1) array
        atomic_masses = np.array([site.specie.atomic_mass for site in self.structure])
        atomic_masses *= AMU_TO_KG #in Kg
        dof = 3 * len(self.structure) - 3

        #scale velocities due to atomic masses
        #mean 0 std proportional to sqrt(1/m)
        velocities = velocities / atomic_masses[:, np.newaxis] ** (1 / 2)

        #remove linear drift (net momentum)
        velocities -= np.average(atomic_masses[:, np.newaxis] * velocities, axis=0) / np.average(atomic_masses)

        #scale velocities to get correct temperature
        energy = np.sum(1 / 2 * atomic_masses * np.sum(velocities ** 2, axis=1))
        scale = (temperature * dof / (2 * energy / BOLTZMANN_CONST)) ** (1 / 2) #boltzmann

        velocities *= scale * 1e-5  #these are in A/fs

        self.temperature = temperature
        self.selective_dynamics = None
        self.predictor_corrector = None
        self.velocities = velocities.tolist() #returns as a list of lists to be consistent with the other initializations


"""**Non-exhaustive** list of valid INCAR tags"""
VALID_INCAR_TAGS = ('NGX', 'NGY', 'NGZ', 'NGXF', 'NGYF', 'NGZF', 'NBANDS',
                    'NBLK', 'SYSTEM', 'NWRITE', 'ENCUT', 'ENAUG', 'PREC',
                    'ISPIN', 'MAGMOM', 'ISTART', 'ICHARG', 'INIWAV', 'NELM',
                    'NELMIN', 'NELMDL', 'EDIFF', 'EDIFFG', 'NSW', 'NBLOCK',
                    'KBLOCK', 'IBRION', 'NFREE', 'POTIM', 'ISIF', 'PSTRESS',
                    'IWAVPR', 'ISYM', 'SYMPREC', 'LCORR', 'TEBEG', 'TEEND',
                    'SMASS', 'NPACO', 'APACO', 'POMASS', 'ZVAL', 'RWIGS',
                    'LORBIT', 'NELECT', 'NUPDOWN', 'EMIN', 'EMAX', 'NEDOS',
                    'ISMEAR', 'SIGMA', 'FERWE', 'FERDO', 'SMEARINGS', 'LREAL',
                    'ROPT', 'GGA', 'VOSKOWN', 'LASPH', 'ALGO', 'IALGO',
                    'LDIAG', 'NSIM', 'IMIX', 'INIMIX', 'MAXMIX', 'AMIX',
                    'BMIX', 'AMIX_MAG', 'BMIX_MAG', 'AMIN', 'MIXPRE', 'WC',
                    'WEIMIN', 'EBREAK', 'DEPER', 'TIME', 'LWAVE', 'LCHARG',
                    'LVTOT', 'LELF', 'NPAR', 'LPLANE', 'LASYNC', 'LSCALAPACK',
                    'LSCALU', 'ISPIND', 'HFSCREEN', 'LHFCALC', 'ENCUTFOCK',
                    'NKRED', 'LMAXMIX', 'PRECFOCK', 'AEXX', 'AGGAX', 'AGGAC',
                    'ALDAC', 'LMAXFOCK', 'LMAXFOCKAE', 'LTHOMAS', 'NKREDX',
                    'NKREDY', 'NKREDZ', 'EVENONLY', 'ODDONLY', 'LDAU', 'LDAUJ',
                    'LDAUL', 'LDAUPRINT', 'LDAUTYPE', 'LDAUU', 'LPEAD',
                    'LCALCPOL', 'LCALCEPS', 'LEFG', 'EFIELD_PEAD',
                    'LNONCOLLINEAR', 'LSORBIT', 'IDIPOL', 'DIPOL', 'LMONO',
                    'LDIPOL', 'EPSILON', 'EFIELD', 'LBERRY', 'IGPAR', 'NPPSTR',
                    'IPEAD', 'I_CONSTRAINED_M', 'LAMBDA', 'M_CONSTR', 'IMAGES',
                    'SPRING', 'LOPTICS', 'CSHIFT', 'LNABLA', 'LEPSILON', 'LRPA',
                    'NOMEGA', 'NOMEGAR', 'LSPECTRAL', 'OMEGAMAX', 'OMEGATL',
                    'ENCUTGW', 'ENCUTGWSOFT', 'ODDONLYGW', 'EVENONLYGW',
                    'LSELFENERGY', 'LRHFATM', 'METAGGA', 'LMAXTAU', 'LCOMPAT',
                    'ENMAX', 'LMAXPAW', 'LSPIRAL', 'LZEROZ', 'LMETAGGA',
                    'ENINI', 'NRMM', 'MREMOVE', 'ADDGRID', 'EFERMI', 'LPARD',
                    'LSCAAWARE', 'IDIOT', 'LMUSIC', 'LREAL_COMPAT',
                    'GGA_COMPAT', 'ICORELEVEL', 'LHFONE', 'LRHFCALC',
                    'LMODELHF', 'ENCUT4O', 'EXXOEP', 'FOURORBIT', 'HFALPHA',
                    'ALDAX', 'SHIFTRED', 'NMAXFOCKAE', 'HFSCREENC', 'MODEL_GW',
                    'MODEL_EPS0', 'MODEL_ALPHA', 'LVEL', 'SAXIS', 'QSPIRAL',
                    'STM', 'KINTER', 'ORBITALMAG', 'LMAGBLOCH', 'LCHIMAG',
                    'LGAUGE', 'MAGATOM', 'MAGDIPOL', 'AVECCONST', 'LTCTE',
                    'LTETE', 'L2ORDER', 'LGWLF', 'ENCUTLF', 'LMAXMP2',
                    'SCISSOR', 'NBANDSGW', 'NBANDSLF', 'DIM', 'ANTIRES',
                    'LUSEW', 'OMEGAGRID', 'SELFENERGY', 'NKREDLFX', 'NKREDLFY',
                    'NKREDLFZ', 'MAXMEM', 'TELESCOPE', 'LCRITICAL_MEM', 'GGA2',
                    'TURBO', 'QUAD_EFG', 'IRESTART', 'NREBOOT', 'NMIN', 'EREF',
                    'KSPACING', 'KGAMMA', 'LSUBROT', 'SCALEE', 'LVHAR',
                    'LORBITALREAL', 'DARWINR', 'DARWINV', 'LFOCKAEDFT',
                    'NUCIND', 'MAGPOS', 'LNICSALL', 'LADDER', 'LHARTREE',
                    'IBSE', 'NBANDSO', 'NBANDSV', 'OPTEMAX', 'LIP')


class Incar(dict, VaspInput):
    """
    INCAR object for reading and writing INCAR files. Essentially consists of
    a dictionary with some helper functions
    """

    def __init__(self, params=dict()):
        """
        Creates an Incar object.
        
        Args:
            params:
                A set of input parameters as a dictionary.
        """
        super(Incar, self).__init__()
        self.update(params)

    def __setitem__(self, key, val):
        """
        Add parameter-val pair to Incar.  Warns if parameter is not in list of
        valid INCAR tags. Also cleans the parameter and val by stripping
        leading and trailing white spaces.
        """
        if key.strip().upper() not in VALID_INCAR_TAGS:
            warnings.warn(key.strip() + " not in VALID_INCAR_TAGS")
        super(Incar, self).__setitem__(key.strip(), Incar.proc_val(key.strip(), val.strip()) if isinstance(val, basestring) else val)

    @property
    def to_dict(self):
        d = {k: v for k, v in self.items()}
        d['module'] = self.__class__.__module__
        d['class'] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        i = Incar()
        for k, v in d.items():
            if k not in ("module", "class"):
                i[k] = v
        return i

    def get_string(self, sort_keys=False, pretty=False):
        """
        Returns a string representation of the INCAR.  The reason why this
        method is different from the __str__ method is to provide options for
        pretty printing.
        
        Args:
            sort_keys:
                Set to True to sort the INCAR parameters alphabetically.
                Defaults to False.
            pretty:
                Set to True for pretty aligned output. Defaults to False.
        """
        keys = self.keys()
        if sort_keys:
            keys = sorted(keys)
        lines = []
        for k in keys:
            if k == "MAGMOM" and isinstance(self[k], list):
                value = []
                for m, g in itertools.groupby(self[k]):
                    value.append("{}*{}".format(len(tuple(g)), m))
                lines.append([k, " ".join(value)])
            elif isinstance(self[k], list):
                lines.append([k, " ".join([str(i) for i in self[k]])])
            else:
                lines.append([k, self[k]])

        if pretty:
            return str_aligned(lines)
        else:
            return str_delimited(lines, None, " = ")

    def __str__(self):
        return self.get_string(sort_keys=True, pretty=False)

    def write_file(self, filename):
        """
        Write Incar to a file.
        
        Args:
            filename:
                filename to write to.
        """
        with open(filename, 'w') as f:
            f.write(self.__str__() + "\n")

    @staticmethod
    def from_file(filename):
        """
        Reads an Incar object from a file.
        
        Args:
            filename - Filename for file
            
        Returns:
            Incar object
        """
        with file_open_zip_aware(filename, "r") as f:
            lines = list(clean_lines(f.readlines()))
        params = {}
        for line in lines:
            m = re.match("(\w+)\s*=\s*(.*)", line)
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
            key:
                INCAR parameter key
            val:
                Actual value of INCAR parameter.
        """
        list_keys = ('LDAUU', 'LDAUL', 'LDAUJ', 'LDAUTYPE', 'MAGMOM')
        bool_keys = ('LDAU', 'LWAVE', 'LSCALU', 'LCHARG', 'LPLANE', 'LHFCALC')
        float_keys = ("EDIFF", "SIGMA", 'TIME', 'ENCUTFOCK', 'HFSCREEN')
        int_keys = ('NSW', 'NELMIN', 'ISIF', 'IBRION', "ISPIN", "ICHARG",
                         "NELM", "ISMEAR", "NPAR", "LDAUPRINT", 'LMAXMIX',
                         'ENCUT', 'NSIM', 'NKRED', 'NUPDOWN', 'ISPIND')

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)
        try:
            if key in list_keys:
                output = list()
                toks = re.split("\s+", val)

                for tok in toks:
                    m = re.match("(\d+)\*([\d\.\-\+]+)", tok)
                    if m:
                        output.extend([smart_int_or_float(m.group(2))] * int(m.group(1)))
                    else:
                        output.append(smart_int_or_float(tok))
                return output
            if key in bool_keys:
                m = re.search("^\W+([TtFf])", val)
                if m:
                    if m.group(1) == "T" or m.group(1) == "t":
                        return True
                    else:
                        return False
                raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(val)

            if key in int_keys:
                return int(val)

        except:
            return val.capitalize()

        return val.capitalize()

    def diff(self, other):
        """
        Diff function for Incar.  Compares two Incars and indicates which
        parameters are the same and which are not. Useful for checking whether
        two runs were done using the same parameters.
        
        Args:
            other:
                The other Incar object to compare to.
        
        Returns:
            Dict of the following format:
            {'Same' : parameters_that_are_the_same,
            'Different': parameters_that_are_different}
            Note that the parameters are return as full dictionaries of values.
            E.g. {'ISIF':3}
        """
        similar_param = {}
        different_param = {}
        for k1, v1 in self.items():
            if k1 not in other:
                different_param[k1] = {"INCAR1": v1, "INCAR2": 'Default'}
            elif v1 != other[k1]:
                different_param[k1] = {"INCAR1": v1, "INCAR2": other[k1]}
            else:
                similar_param[k1] = v1
        for k2, v2 in other.items():
            if k2 not in similar_param and k2 not in different_param:
                if k2 not in self:
                    different_param[k2] = {"INCAR1": 'Default', "INCAR2": v2}
        return {'Same' : similar_param, 'Different': different_param}

    def __add__(self, other):
        """
        Add all the values of another INCAR object to this object.
        Facilitates the use of "standard" INCARs.
        """
        params = {k:v for k, v in self.items()}
        for k, v in other.items():
            if k in self and v != self[k]:
                raise ValueError("Incars have conflicting values!")
            else:
                params[k] = v
        return Incar(params)


class Kpoints(VaspInput):
    """
    KPOINT reader/writer.
    """
    supported_modes = Enum(("Gamma", "Monkhorst", "Automatic", "Line_mode",
                            "Cartesian", "Reciprocal"))

    def __init__(self, comment="Default gamma", num_kpts=0,
                 style=supported_modes.Gamma,
                 kpts=[[1, 1, 1]], kpts_shift=(0, 0, 0),
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
            comment:
                String comment for Kpoints
            num_kpts:
                Following VASP method of defining the KPOINTS file, this
                parameter is the number of kpoints specified. If set to 0
                (or negative), VASP automatically generates the KPOINTS.
            style:
                Style for generating KPOINTS.  Use one of the
                Kpoints.supported_modes enum types.
            kpts:
                2D array of kpoints.  Even when only a single specification is
                required, e.g. in the automatic scheme, the kpts should still
                be specified as a 2D array. e.g., [[20]] or [[2,2,2]].
            kpts_shift:
                Shift for Kpoints.
            kpts_weights:
                Optional weights for kpoints.  For explicit kpoints.
            coord_type:
                In line-mode, this variable specifies whether the Kpoints were
                given in Cartesian or Reciprocal coordinates.
            labels:
                In line-mode, this should provide a list of labels for each kpt.
            tet_number:
                For explicit kpoints, specifies the number of tetrahedrons for
                the tetrahedron method.
            tet_weight:
                For explicit kpoints, specifies the weight for each tetrahedron
                for the tetrahedron method.
            tet_connections:
                For explicit kpoints, specifies the connections of the
                tetrahedrons for the tetrahedron method.
                Format is a list of tuples, [ (sym_weight, [tet_vertices]), ...]           
        
        The default behavior of the constructor is for a Gamma centered,
        1x1x1 KPOINTS with no shift.
        """
        if num_kpts > 0 and (not labels) and (not kpts_weights):
            raise ValueError("For explicit or line-mode kpoints, either the labels or kpts_weights must be specified.")
        if style in (Kpoints.supported_modes.Automatic, Kpoints.supported_modes.Gamma, Kpoints.supported_modes.Monkhorst) and len(kpts) > 1:
            raise ValueError("For fully automatic or automatic gamma or monk kpoints, only a single line for the number of divisions is allowed.")

        self.comment = comment
        self.num_kpts = num_kpts
        self.style = style
        self.coord_type = coord_type
        self.kpts = kpts
        self.kpts_weights = kpts_weights
        self.kpts_shift = kpts_shift
        self.labels = labels
        self.tet_number = tet_number
        self.tet_weight = tet_weight
        self.tet_connections = tet_connections

    @staticmethod
    def automatic(subdivisions):
        """
        Convenient static constructor for a fully automatic Kpoint grid, with 
        gamma centered Monkhorst-Pack grids and the number of subdivisions
        along each reciprocal lattice vector determined by the scheme in the
        VASP manual.
        
        Args:
            subdivisions:
                 Parameter determining number of subdivisions along each
                 reciprocal lattice vector.
                 
        Returns:
            Kpoints object
        """
        return Kpoints("Fully automatic kpoint scheme", 0, style=Kpoints.supported_modes.Automatic, kpts=[[subdivisions]])

    @staticmethod
    def gamma_automatic(kpts=(1, 1, 1), shift=(0, 0, 0)):
        """
        Convenient static constructor for an automatic Gamma centered Kpoint
        grid.
        
        Args:
            kpts:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
                Defaults to (1,1,1)
            shift:
                Shift to be applied to the kpoints. Defaults to (0,0,0).
                
        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0, Kpoints.supported_modes.Gamma, kpts=[kpts], kpts_shift=shift)

    @staticmethod
    def monkhorst_automatic(kpts=(2, 2, 2), shift=(0, 0, 0)):
        """
        Convenient static constructor for an automatic Monkhorst pack Kpoint
        grid.
        
        Args:
            kpts:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
                Defaults to (2,2,2)
            shift:
                Shift to be applied to the kpoints. Defaults to (0,0,0).
                
        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0, Kpoints.supported_modes.Monkhorst, kpts=[kpts], kpts_shift=shift)

    @staticmethod
    def automatic_density(structure, kppa):
        '''
        Returns an automatic Kpoint object based on a structure and a kpoint 
        density. Uses Gamma centered meshes for hexagonal cells and 
        Monkhorst-Pack grids otherwise.
        
        Algorithm: 
            Uses a simple approach scaling the number of divisions along each 
            reciprocal lattice vector proportional to its length. 
            
        Args:
            structure:
                Input structure
            kppa:
                Grid density
        '''

        latt = structure.lattice
        lengths = latt.abc
        ngrid = kppa / structure.num_sites

        mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3)

        num_div = [int(round(1 / lengths[i] * mult)) for i in xrange(3)]
        #ensure that numDiv[i] > 0
        num_div = [i if i > 0 else 1 for i in num_div]

        angles = latt.angles
        hex_angle_tol = 5 #in degrees
        hex_length_tol = 0.01 #in angstroms
        right_angles = [i for i in xrange(3) if abs(angles[i] - 90) < hex_angle_tol]
        hex_angles = [i for i in xrange(3) if abs(angles[i] - 60) < hex_angle_tol or abs(angles[i] - 120) < hex_angle_tol]

        is_hexagonal = (len(right_angles) == 2 and len(hex_angles) == 1 and abs(lengths[right_angles[0]] == lengths[right_angles[1]]) < hex_length_tol)

        style = Kpoints.supported_modes.Gamma
        if not is_hexagonal:
            num_div = [i + i % 2 for i in num_div]
            style = Kpoints.supported_modes.Monkhorst
        comment = "pymatgen generated Materials Project kpoints with grid density = " + str(kppa) + ' per atom.'
        num_kpts = 0
        return Kpoints(comment, num_kpts, style, [num_div], [0, 0, 0])

    @staticmethod
    def from_file(filename):
        """
        Reads a Kpoints object from a KPOINTS file.
        
        Args:
            filename:
                filename to read from.
                
        Returns:
            Kpoints object
        """
        with file_open_zip_aware(filename) as f:
            lines = [line.strip() for line in f.readlines()]
        comment = lines[0]
        num_kpts = int(lines[1].split()[0].strip())
        style = lines[2].lower()[0]

        #Fully automatic KPOINTS
        if style == "a":
            return Kpoints.automatic(int(lines[3]))

        coord_pattern = re.compile("^\s*([\d+\.\-Ee]+)\s+([\d+\.\-Ee]+)\s+([\d+\.\-Ee]+)")

        #Automatic gamma and Monk KPOINTS, with optional shift
        if style == "g" or style == "m":
            kpts = [int(x) for x in lines[3].split()]
            kpts_shift = (0, 0, 0)
            if len(lines) > 4 and coord_pattern.match(lines[4]):
                try:
                    kpts_shift = [int(x) for x in lines[4].split()]
                except:
                    pass
            return Kpoints.gamma_automatic(kpts, kpts_shift) if style == "g" else Kpoints.monkhorst_automatic(kpts, kpts_shift)

        #Automatic kpoints with basis
        if num_kpts <= 0:
            style = Kpoints.supported_modes.Cartesian if style in "ck" else Kpoints.supported_modes.Reciprocal
            kpts = [[float(x) for x in lines[i].split()] for i in xrange(3, 6)]
            kpts_shift = [float(x) for x in lines[6].split()]
            return Kpoints(comment=comment, num_kpts=num_kpts, style=style, kpts=kpts, kpts_shift=kpts_shift)

        #Line-mode KPOINTS, usually used with band structures
        if style == 'l':
            coord_type = 'Cartesian' if lines[3].lower()[0] in 'ck' else 'Reciprocal'
            style = Kpoints.supported_modes.Line_mode
            kpts = []
            labels = []
            patt = re.compile('([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s*!\s*(.*)')
            for i in range(4, len(lines)):
                line = lines[i]
                m = patt.match(line)
                if m:
                    kpts.append([float(m.group(1)), float(m.group(2)), float(m.group(3))])
                    labels.append(m.group(4).strip())
            return Kpoints(comment=comment, num_kpts=num_kpts, style=style,
                 kpts=kpts, coord_type=coord_type, labels=labels)

        #Assume explicit KPOINTS if all else fails.
        style = Kpoints.supported_modes.Cartesian if style == "ck" else Kpoints.supported_modes.Reciprocal
        kpts = []
        kpts_weights = []
        tet_number = 0
        tet_weight = 0
        tet_connections = None

        for i in xrange(3, 3 + num_kpts):
            toks = re.split("\s+", lines[i])
            kpts.append([float(toks[0]), float(toks[1]), float(toks[2])])
            kpts_weights.append(float(toks[3]))
        try:
            #Deal with tetrahedron method
            if lines[3 + num_kpts].strip().lower()[0] == 't':
                toks = lines[4 + num_kpts].split()
                tet_number = int(toks[0])
                tet_weight = float(toks[1])
                tet_connections = []
                for i in xrange(5 + num_kpts, 5 + num_kpts + tet_number):
                    toks = lines[i].split()
                    tet_connections.append((int(toks[0]), [int(toks[j]) for j in xrange(1, 5)]))
        except:
            pass

        return Kpoints(comment=comment, num_kpts=num_kpts, style=style,
                 kpts=kpts, kpts_weights=kpts_weights,
                 tet_number=tet_number, tet_weight=tet_weight, tet_connections=tet_connections)

    def write_file(self, filename):
        """
        Write Kpoints to a file.
        
        Args:
            filename:
                filename to write to.
        """
        with open(filename, 'w') as f:
            f.write(self.__str__() + "\n")

    def __str__(self):
        lines = []
        lines.append(self.comment)
        lines.append(str(self.num_kpts))
        lines.append(self.style)
        style = self.style.lower()[0]
        if style == "Line-mode":
            lines.append(self.coord_type)
        for i in xrange(len(self.kpts)):
            lines.append(" ".join([str(x) for x in self.kpts[i]]))
            if style == "l":
                lines[-1] += " ! " + self.labels[i]
            elif self.num_kpts > 0:
                lines[-1] += " %f" % (self.kpts_weights[i])

        #Print tetrahedorn parameters if the number of tetrahedrons > 0
        if style not in "lagm" and self.tet_number > 0:
            lines.append("Tetrahedron")
            lines.append("%d %f" % (self.tet_number, self.tet_weight))
            for sym_weight, vertices in self.tet_connections:
                lines.append("%d %d %d %d %d" % (sym_weight, vertices[0], vertices[1], vertices[2], vertices[3]))

        #Print shifts for automatic kpoints types if not zero.
        if self.num_kpts <= 0 and tuple(self.kpts_shift) != (0, 0, 0):
            lines.append(" ".join([str(x) for x in self.kpts_shift]))
        return "\n".join(lines)

    @property
    def to_dict(self):
        """json friendly dict representation of Kpoints"""
        d = {'comment': self.comment, 'nkpoints' : self.num_kpts, 'generation_style' : self.style, 'kpoints': self.kpts, 'usershift': self.kpts_shift}
        optional_paras = ['genvec1', 'genvec2', 'genvec3', 'shift']
        for para in optional_paras:
            if para in self.__dict__:
                d[para] = self.__dict__[para]
        d['module'] = self.__class__.__module__
        d['class'] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        comment = d.get('comment', '')
        generation_style = d.get('generation_style')
        kpts = d.get('kpoints', [[1, 1, 1]])
        kpts_shift = d.get('usershift', [0, 0, 0])
        num_kpts = d.get('nkpoints', 0)
        #coord_type = d.get('coord_type', None)
        return Kpoints(comment=comment, kpts=kpts, style=generation_style, kpts_shift=kpts_shift, num_kpts=num_kpts)



#set the global VASP_PSP_DIR
VASP_PSP_DIR = None
if 'VASP_PSP_DIR' in os.environ:
    VASP_PSP_DIR = os.environ['VASP_PSP_DIR']
elif os.path.exists(os.path.join(os.path.dirname(pymatgen.__file__), "pymatgen.cfg")):
    module_dir = os.path.dirname(pymatgen.__file__)
    config = ConfigParser.SafeConfigParser()
    config.readfp(open(os.path.join(module_dir, "pymatgen.cfg")))
    VASP_PSP_DIR = config.get('VASP', 'pspdir')


class PotcarSingle(object):
    """
    Object for a **single** POTCAR. The builder assumes the complete string is
    the POTCAR contains the complete untouched data in "data" as a string and
    a dict of keywords.
    
    .. attribute:: data
    
        POTCAR data as a string.
        
    .. attribute:: keywords
    
        Keywords parsed from the POTCAR as a dict. All keywords are also
        accessible as attributes in themselves. E.g., potcar.enmax,
        potcar.encut, etc.
    """
    functional_dir = {'PBE':'POT_GGA_PAW_PBE', 'LDA':'POT_LDA_PAW',
                      'PW91':'POT_GGA_PAW_PW91'}


    def __init__(self, data):
        """
        Args:
            data:
                Complete and single potcar file as a string.
        """
        self.data = data  # raw POTCAR as a string

        # AJ (5/18/2012) - only search on relevant portion of POTCAR, should
        # fail gracefully if string not found
        search_string = self.data[0:self.data.find("END of PSCTR-controll parameters")]
        keypairs = re.compile(r";*\s*(.+?)\s*=\s*([^;\n]+)\s*", re.M).findall(search_string)
        self.keywords = dict(keypairs)

    def __str__(self):
        return self.data

    def write_file(self, filename):
        writer = open(filename, 'w')
        writer.write(self.__str__() + "\n")
        writer.close()

    @staticmethod
    def from_file(filename):
        with file_open_zip_aware(filename, "rb") as f:
            return PotcarSingle(f.read())

    @staticmethod
    def from_symbol_and_functional(symbol, functional="PBE"):
        funcdir = PotcarSingle.functional_dir[functional]
        paths_to_try = [os.path.join(VASP_PSP_DIR, funcdir, "POTCAR.{}.gz".format(symbol)),
                        os.path.join(VASP_PSP_DIR, funcdir, symbol, "POTCAR")]
        for p in paths_to_try:
            if os.path.exists(p):
                return PotcarSingle.from_file(p)
        raise IOError("You do not have the right POTCAR with functional {} and label {} in your VASP_PSP_DIR".format(functional, symbol))

    @property
    def symbol(self):
        """
        Symbol of POTCAR, e.g., Fe_pv
        """
        return self.keywords['TITEL'].split(" ")[1].strip()

    @property
    def element(self):
        """
        Attempt to return the atomic symbol based on the VRHFIN keyword.
        """
        element = self.keywords['VRHFIN'].split(":")[0].strip()
        #VASP incorrectly gives the element symbol for Xe as 'X'
        return 'Xe' if element == 'X' else element

    @property
    def atomic_no(self):
        """
        Attempt to return the atomic number based on the VRHFIN keyword.
        """
        return Element(self.element).Z

    @property
    def nelectrons(self):
        return self.zval

    def __getattr__(self, a):
        """
        Delegates attributes to keywords. For example, you can use
        potcarsingle.enmax to get the ENMAX of the POTCAR.
        
        For float type properties, they are converted to the correct float. By
        default, all energies in eV and all length scales are in Angstroms.
        """
        floatkeywords = ['DEXC', 'RPACOR', 'ENMAX', 'QCUT', 'EAUG', 'RMAX',
                         'ZVAL', 'EATOM', 'NDATA', 'QGAM', 'ENMIN', 'RCLOC',
                         'RCORE', 'RDEP', 'RAUG', 'POMASS', 'RWIGS']
        a_caps = a.upper()
        if a_caps in self.keywords:
            return self.keywords[a_caps] if a_caps not in floatkeywords else float(self.keywords[a_caps].split()[0])
        raise AttributeError(a)


class Potcar(list, VaspInput):
    """
    Object for reading and writing POTCAR files for calculations. Consists of a
    list of PotcarSingle.
    """

    DEFAULT_FUNCTIONAL = "PBE"

    """
    Cache for PotcarSingles. Results in orders of magnitude faster generation
    of output when doing high-throughput run generation.
    """
    _cache = {}

    def __init__(self, symbols=None, functional=DEFAULT_FUNCTIONAL,
                 sym_potcar_map=None):
        """
        Args:
            symbols:
                Element symbols for POTCAR
            functional:
                Functional used.
            sym_potcar_map:
                Allows a user to specify a specific element symbol to POTCAR
                symbol mapping. For example, you can have {'Fe':'Fe_pv'} to
                specify that the Fe_pv psuedopotential should be used for Fe.
                Default is None, which uses a pre-determined mapping used in
                the Materials Project.
        """
        if symbols is not None:
            self.functional = functional
            self.set_symbols(symbols, functional, sym_potcar_map)

    @property
    def to_dict(self):
        d = {'functional': self.functional, 'symbols': self.symbols}
        d['module'] = self.__class__.__module__
        d['class'] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        functional = d['functional']
        symbols = d['symbols']
        return Potcar(symbols=symbols, functional=functional)

    @staticmethod
    def from_file(filename):
        with file_open_zip_aware(filename, "r") as reader:
            fdata = reader.read()
        potcar = Potcar()
        potcar_strings = re.compile(r"\n{0,1}\s*(.*?End of Dataset)", re.S).findall(fdata)
        for p in potcar_strings:
            potcar.append(PotcarSingle(p))
        return potcar

    def __str__(self):
        return "".join([str(potcar) for potcar in self])    #line break not used because there is already one at the end of str(potcar) and it causes VASP issues

    def write_file(self, filename):
        """
        Write Potcar to a file.
        
        Args:
            filename:
                filename to write to.
        """
        with open(filename, 'w') as f:
            f.write(self.__str__() + "\n")

    @property
    def symbols(self):
        """
        Get the atomic symbols of all the atoms in the POTCAR file.
        """
        return [p.symbol for p in self]

    def set_symbols(self, symbols, functional=DEFAULT_FUNCTIONAL, sym_potcar_map=None):
        '''
        Initialize the POTCAR from a set of symbols. Currently, the POTCARs can
        be fetched from a location specified in the environment variable 
        VASP_PSP_DIR or in a pymatgen.cfg or specified explicitly in a map.
        
        Args:
            symbols:
                A list of element symbols
            functional:
                (optional) the functional to use from the config file
            sym_potcar_map:
                (optional) a map of symbol:raw POTCAR string. If sym_potcar_map
                is specified, POTCARs will be generated from the given map data
                rather than the config file location.
        '''
        del self[:]
        if sym_potcar_map:
            for el in symbols:
                self.append(PotcarSingle(sym_potcar_map[el]))
        else:
            for el in symbols:
                if (el, functional) in Potcar._cache:
                    self.append(Potcar._cache[(el, functional)])
                else:
                    p = PotcarSingle.from_symbol_and_functional(el, functional)
                    self.append(p)
                    Potcar._cache[(el, functional)] = p


class Vasprun(object):
    """
    Vastly improved sax-based parser for vasprun.xml files.
    Speedup over Dom is at least 2x for smallish files (~1Mb) to orders of
    magnitude for larger files (~10Mb). All data is stored as attributes, which
    are delegated to the VasprunHandler object. Note that the results would
    differ depending on whether the read_electronic_structure option is set to
    True.
    
    **Vasp results**
        
    .. attribute:: ionic_steps
    
        All ionic steps in the run as a list of
        {'structure': structure at end of run,
        'electronic_steps': {All electronic step data in vasprun file},
        'stresses': stress matrix}
    
    .. attribute:: structures
    
        List of Structure objects for the structure at each ionic step.
    
    .. attribute:: tdos
     
        Total dos calculated at the end of run.
    
    .. attribute:: idos 
        
        Integrated dos calculated at the end of run.
    
    .. attribute:: pdos 
        
        List of list of PDos objects. Access as pdos[atomindex][orbitalindex]
    
    .. attribute:: efermi
        
        Fermi energy
    
    .. attribute:: eigenvalues 
        
        Available only if parse_eigen=True. Final eigenvalues as a dict of
        {(spin, kpoint index):[[eigenvalue, occu]]}.
        This representation is based on actual ordering in VASP and is meant as
        an intermediate representation to be converted into proper objects. The
        kpoint index is 0-based (unlike the 1-based indexing in VASP).
  
    .. attribute:: projected_eigenvalues 
        
        Final projected eigenvalues as a dict of
        {(atom index, band index, kpoint index, Orbital, Spin):float}
        This representation is based on actual ordering in VASP and is meant as
        an intermediate representation to be converted into proper objects. The
        kpoint, band and atom indices are 0-based (unlike the 1-based indexing
        in VASP).
         
    **Vasp inputs**
        
    .. attribute:: incar
        
        Incar object for parameters specified in INCAR file.
    
    .. attribute:: parameters
        
        Incar object with parameters that vasp actually used, including all
        defaults.
    
    .. attribute:: kpoints
        
        Kpoints object for KPOINTS specified in run.
    
    .. attribute:: actual_kpoints
        
        List of actual kpoints, e.g.,
        [[0.25, 0.125, 0.08333333], [-0.25, 0.125, 0.08333333],
        [0.25, 0.375, 0.08333333], ....]
    
    .. attribute:: actual_kpoints_weights
        
        List of kpoint weights, E.g.,
        [0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, ....]
    
    .. attribute:: atomic_symbols
        
        List of atomic symbols, e.g., ['Li', 'Fe', 'Fe', 'P', 'P', 'P']
    
    .. attribute:: potcar_symbols
        
        List of POTCAR symbols. e.g., 
        ['PAW_PBE Li 17Jan2003', 'PAW_PBE Fe 06Sep2000', ..]
            
    Author: Shyue Ping Ong
    """
    supported_properties = ['lattice_rec', 'vasp_version', 'incar',
                            'parameters', 'potcar_symbols', 'atomic_symbols',
                            'kpoints', 'actual_kpoints', 'structures',
                            'actual_kpoints_weights', 'dos_energies',
                            'eigenvalues', 'tdos', 'idos', 'pdos', 'efermi',
                            'ionic_steps', 'dos_has_errors',
                            'projected_eigenvalues']

    def __init__(self, filename, ionic_step_skip=None,
                 parse_dos=True, parse_eigen=True, parse_projected_eigen=False):
        """
        Args:
            filename:
                Filename to parse
            ionic_step_skip:
                If ionic_step_skip is a number > 1, only every ionic_step_skip
                ionic steps will be read for structure and energies. This is
                very useful if you are parsing very large vasprun.xml files and
                you are not interested in every single ionic step. Note that the
                initial and final structure of all runs will always be read,
                regardless of the ionic_step_skip.
            parse_dos:
                Whether to parse the dos. Defaults to True. Set
                to False to shave off significant time from the parsing if you
                are not interested in getting those data.
            parse_eigen:
                Whether to parse the eigenvalues. Defaults to True. Set
                to False to shave off significant time from the parsing if you
                are not interested in getting those data.
            parse_projected_eigen:
                Whether to parse the projected eigenvalues. Defaults to False.
                Set to True to obtain projected eigenvalues. **Note that this
                can take an extreme amount of time and memory.** So use this
                wisely.
        """
        self.filename = filename

        with file_open_zip_aware(filename) as f:
            self._handler = VasprunHandler(filename, parse_dos=parse_dos,
                                    parse_eigen=parse_eigen,
                                    parse_projected_eigen=parse_projected_eigen)
            if ionic_step_skip == None:
                self._parser = xml.sax.parse(f, self._handler)
            else:
                #remove parts of the xml file and parse the string
                run = f.read()
                steps = run.split('<calculation>')
                new_steps = steps[::int(ionic_step_skip)]
                #add the last step from the run
                if steps[-1] != new_steps[-1]:
                    new_steps.append(steps[-1])
                self._parser = xml.sax.parseString('<calculation>'.join(new_steps), self._handler)
            for k in Vasprun.supported_properties:
                setattr(self, k, getattr(self._handler, k))

    @property
    def converged(self):
        """
        True if a relaxation run is converged.  Always True for a static run.
        """
        return len(self.structures) - 2 < self.parameters['NSW'] or self.parameters['NSW'] == 0

    @property
    def final_energy(self):
        """
        Final energy from the vasp run.
        """
        return self.ionic_steps[-1]['electronic_steps'][-1]['e_wo_entrp']

    @property
    def final_structure(self):
        """
        Final structure from vasprun.
        """
        return self.structures[-1]

    @property
    def initial_structure(self):
        """
        Initial structure from vasprun.
        """
        return self.structures[0]

    @property
    def complete_dos(self):
        """
        A complete dos object which incorporates the total dos and all projected dos.
        """
        final_struct = self.final_structure
        pdoss = {final_struct[i]:pdos for i, pdos in enumerate(self.pdos)}
        return CompleteDos(self.final_structure, self.tdos, pdoss)

    @property
    def hubbards(self):
        """
        Hubbard U values used if a vasprun is a GGA+U run. Empty dict otherwise.
        """
        symbols = [re.split("\s+", s)[1] for s in self.potcar_symbols]
        symbols = [re.split("_", s)[0] for s in symbols]
        if not self.incar.get('LDAU', False):
            return {}
        us = self.incar.get('LDAUU', self.parameters.get('LDAUU'))
        js = self.incar.get('LDAUJ', self.parameters.get('LDAUJ'))
        if len(us) == len(symbols):
            return { symbols[i] : us[i] - js[i] for i in xrange(len(symbols))}
        elif sum(us) == 0 and sum(js) == 0:
            return {}
        else:
            raise VaspParserError("Length of U value parameters and atomic symbols are mismatched")

    @property
    def run_type(self):
        """
        Returns the run type. Currently supports only GGA and HF calcs. 
        
        TODO: Fix for other functional types like LDA, PW91, etc.
        """
        if self.is_hubbard:
            return "GGA+U"
        elif self.parameters.get('LHFCALC', False):
            return "HF"
        else:
            return "GGA"

    @property
    def is_hubbard(self):
        """
        True if run is a DFT+U run.
        """
        if len(self.hubbards) == 0:
            return False
        return sum(self.hubbards.values()) > 0

    @property
    def is_spin(self):
        """
        True if run is spin-polarized.
        """
        return True if self.incar.get('ISPIN', 1) == 2 else False

    def get_band_structure(self, kpoints_filename=None, efermi=None):
        """
        Returns the band structure as a BandStructureSymmLine object
        
        Args:
            kpoints_filename:
                Full path of the KPOINTS file from which the band structure is
                generated.
                If none is provided, the code will try to intelligently
                determine the appropriate KPOINTS file by substituting the
                filename of the vasprun.xml with KPOINTS.
                The latter is the default behavior.
            efermi:
                If you want to specify manually the fermi energy this is where
                you should do it. By default, the None value means the code
                will get it from the vasprun.
                
        Returns:
            a BandStructure object (or more specifically a BandStructureSymmLine object if 
            the run is detected to be a run along symmetry lines)
        
        TODO:
            - make a bit more general for non Symm Line band structures
            - make a decision on the convention with 2*pi or not 
        """
        if not kpoints_filename:
            kpoints_filename = self.filename.replace('vasprun.xml', 'KPOINTS')
        if not os.path.exists(kpoints_filename):
            raise VaspParserError('KPOINTS file needed to obtain band structure.')
        if not self.incar['ICHARG'] == 11:
            raise VaspParserError('band structure runs have to be non-self consistent (ICHARG=11)')

        if efermi == None:
            efermi = self.efermi

        kpoint_file = Kpoints.from_file(kpoints_filename)
        lattice_new = Lattice(self.lattice_rec.matrix * 2 * math.pi)
        #lattice_rec=[self.lattice_rec.matrix[i][j] for i,j in range(3)]

        kpoints = [np.array(self.actual_kpoints[i]) for i in range(len(self.actual_kpoints))]
        dict_eigen = self.to_dict['output']['eigenvalues']

        eigenvals = {}
        if dict_eigen['1'].has_key('up') and dict_eigen['1'].has_key('down') and self.incar['ISPIN'] == 2:
            eigenvals = {Spin.up:[], Spin.down:[]}
        else:
            eigenvals = {Spin.up:[]}

        neigenvalues = [len(v['up']) for v in dict_eigen.values()]
        min_eigenvalues = min(neigenvalues)

        for i in range(min_eigenvalues):
            eigenvals[Spin.up].append([dict_eigen[str(j)]['up'][i][0] for j in range(len(kpoints))]);
        if eigenvals.has_key(Spin.down):
            for i in range(min_eigenvalues):
                eigenvals[Spin.down].append([dict_eigen[str(j)]['down'][i][0] for j in range(len(kpoints))]);
        if kpoint_file.style == "Line_mode":
            labels_dict = dict(zip(kpoint_file.labels, kpoint_file.kpts))
            return BandStructureSymmLine(kpoints, eigenvals, lattice_new, efermi, labels_dict)
        else:
            return BandStructure(kpoints, eigenvals, lattice_new, efermi)

    @property
    def eigenvalue_band_properties(self):
        """
        Band properties from the eigenvalues as a tuple,
        (band gap, cbm, vbm, is_band_gap_direct).
        """
        vbm = -float('inf')
        vbm_kpoint = None
        cbm = float('inf')
        cbm_kpoint = None
        for k, val in self.eigenvalues.items():
            for (eigenval, occu) in val:
                if occu > 1e-8 and eigenval > vbm:
                    vbm = eigenval
                    vbm_kpoint = k[0]
                elif occu <= 1e-8 and eigenval < cbm:
                    cbm = eigenval
                    cbm_kpoint = k[0]
        return (cbm - vbm, cbm, vbm, vbm_kpoint == cbm_kpoint)

    @property
    def to_dict(self):
        """
        Json-serializable dict representation.
        """
        d = {}
        d['vasp_version'] = self.vasp_version
        d['has_vasp_completed'] = self.converged
        d['nsites'] = len(self.final_structure)
        comp = self.final_structure.composition
        d['unit_cell_formula'] = comp.to_dict
        d['reduced_cell_formula'] = Composition(comp.reduced_formula).to_dict
        d['pretty_formula'] = comp.reduced_formula
        symbols = [re.split("\s+", s)[1] for s in self.potcar_symbols]
        symbols = [re.split("_", s)[0] for s in symbols]
        d['elements'] = symbols
        d['nelements'] = len(symbols)
        d['is_hubbard'] = self.incar.get('LDAU', False)
        if d['is_hubbard']:
            us = self.incar.get('LDAUU', self.parameters.get('LDAUU'))
            js = self.incar.get('LDAUJ', self.parameters.get('LDAUJ'))
            if len(us) == len(symbols):
                d['hubbards'] = {symbols[i]:us[i] - js[i] for i in xrange(len(symbols))}
            elif sum(us) == 0 and sum(js) == 0:
                d['is_hubbard'] = False
                d['hubbards'] = {}
            else:
                raise VaspParserError("Length of U value parameters and atomic symbols are mismatched")
        else:
            d['hubbards'] = {}

        d['run_type'] = self.run_type

        vasp_input = {}
        vasp_input['incar'] = {k:v for k, v in self.incar.items()}
        vasp_input['crystal'] = self.initial_structure.to_dict
        vasp_input['kpoints'] = self.kpoints.to_dict
        vasp_input['kpoints']['actual_points'] = [{'abc':list(self.actual_kpoints[i]), 'weight':self.actual_kpoints_weights[i]} for i in xrange(len(self.actual_kpoints))]
        vasp_input['potcar'] = [s.split(" ")[1] for s in self.potcar_symbols]
        vasp_input['parameters'] = {k:v for k, v in self.parameters.items()}
        vasp_input['lattice_rec'] = self.lattice_rec.to_dict
        d['input'] = vasp_input

        vasp_output = {}
        vasp_output['ionic_steps'] = self.ionic_steps
        vasp_output['final_energy'] = self.final_energy
        vasp_output['final_energy_per_atom'] = self.final_energy / len(self.final_structure)
        vasp_output['crystal'] = self.final_structure.to_dict
        vasp_output['efermi'] = self.efermi
        vasp_output['eigenvalues'] = {}
        for (spin, index), values in self.eigenvalues.items():
            if str(index) not in vasp_output['eigenvalues']:
                vasp_output['eigenvalues'][str(index)] = {str(spin):values}
            else:
                vasp_output['eigenvalues'][str(index)][str(spin)] = values

        (gap, cbm, vbm, is_direct) = self.eigenvalue_band_properties
        vasp_output.update(dict(bandgap=gap, cbm=cbm, vbm=vbm, is_gap_direct=is_direct))
        d['output'] = vasp_output

        return clean_json(d, strict=True)


class VasprunHandler(xml.sax.handler.ContentHandler):
    """
    Sax handler for vasprun.xml. Attributes are mirrored into Vasprun object.
    Generally should not be initiatized on its own.
    """

    def __init__(self, filename, parse_dos=True, parse_eigen=True,
                 parse_projected_eigen=False):
        self.filename = filename
        self.parse_dos = parse_dos
        self.parse_eigen = parse_eigen
        self.parse_projected_eigen = parse_projected_eigen

        self.step_count = 0
        # variables to be filled
        self.vasp_version = None
        self.incar = Incar()
        self.parameters = Incar()
        self.potcar_symbols = []
        self.atomic_symbols = []
        self.kpoints = Kpoints()
        self.actual_kpoints = []
        self.actual_kpoints_weights = []
        self.dos_energies = None

        #  will  be  {(spin, kpoint index): [[energy, occu]]}
        self.eigenvalues = {}

        #{(spin, kpoint_index, band_index, atom_ind, orb):float}
        self.projected_eigenvalues = {}

        self.tdos = {}
        self.idos = {}
        self.pdos = {}
        self.efermi = None
        self.ionic_steps = [] # should be a list of dict
        self.structures = []
        self.lattice_rec = []
        self.stress = []

        self.input_read = False
        self.read_structure = False
        self.read_rec_lattice = False
        self.read_calculation = False
        self.read_eigen = False
        self.read_projected_eigen = False
        self.read_dos = False
        self.in_efermi = False
        self.read_atoms = False
        self.read_lattice = False
        self.read_positions = False
        self.incar_param = None

        #Intermediate variables
        self.dos_energies_val = []
        self.dos_val = []
        self.idos_val = []
        self.raw_data = []

        #will be set to true if there is an error parsing the Dos.
        self.dos_has_errors = False
        self.state = defaultdict(bool)

    def startElement(self, name, attributes):
        self.state[name] = attributes.get('name', True)
        self.read_val = False

        #Nested if loops makes reading much faster.
        if not self.input_read: #reading input parameters
            self._init_input(name, attributes)
        else: #reading calculations and structures and eigenvalues.
            self._init_calc(name, attributes)
        if self.read_val:
            self.val = StringIO.StringIO()

    def _init_input(self, name, attributes):
        if (name == "i" or name == "v") and (self.state['incar'] or self.state['parameters']):
            self.incar_param = attributes['name']
            self.param_type = 'float' if 'type' not in attributes else attributes['type']
            self.read_val = True
        elif name == "v" and self.state['kpoints']:
            self.read_val = True
        elif name == "generation" and self.state['kpoints']:
            self.kpoints.comment = "Kpoints from vasprun.xml"
            self.kpoints.num_kpts = 0
            self.kpoints.style = attributes['param']
            self.kpoints.kpts = []
            self.kpoints.kpts_shift = [0, 0, 0]
        elif name == "c" and (self.state['array'] == "atoms" or self.state['array'] == "atomtypes"):
            self.read_val = True
        elif name == "i" and self.state['i'] == "version" and self.state['generator']:
            self.read_val = True

    def _init_calc(self, name, attributes):
        if self.read_structure and name == "v":
            if self.state['varray'] == 'basis':
                self.read_lattice = True
            elif self.state['varray'] == 'positions':
                self.read_positions = True
            elif self.state['varray'] == 'rec_basis':
                self.read_rec_lattice = True
        elif self.read_calculation:
            if name == "i" and self.state['scstep']:
                logger.debug('Reading scstep...')
                self.read_val = True
            elif name == "v" and (self.state['varray'] == "forces" or self.state['varray'] == "stress"):
                self.read_positions = True
            elif name == "dos" and self.parse_dos:
                logger.debug('Reading dos...')
                self.dos_energies = None
                self.tdos = {}
                self.idos = {}
                self.pdos = {}
                self.efermi = None
                self.read_dos = True
            elif name == "eigenvalues" and self.parse_eigen and (not self.state['projected']):
                logger.debug('Reading eigenvalues. Projected = {}'.format(self.state['projected']))
                self.eigenvalues = {}
                self.read_eigen = True
            elif name == "eigenvalues" and self.parse_projected_eigen and self.state['projected']:
                logger.debug('Reading projected eigenvalues...')
                self.projected_eigen = {}
                self.read_projected_eigen = True
            elif self.read_eigen or self.read_projected_eigen:
                if name == "r" and self.state["set"]:
                    self.read_val = True
                elif name == "set" and "comment" in attributes:
                    comment = attributes["comment"]
                    self.state["set"] = comment
                    if comment.startswith("spin"):
                        self.eigen_spin = Spin.up if self.state["set"] in ["spin 1", "spin1"] else Spin.down
                        logger.debug("Reading spin {}".format(self.eigen_spin))
                    elif comment.startswith("kpoint"):
                        self.eigen_kpoint = int(comment.split(" ")[1])
                        logger.debug("Reading kpoint {}".format(self.eigen_kpoint))
                    elif comment.startswith("band"):
                        self.eigen_band = int(comment.split(" ")[1])
                        logger.debug("Reading band {}".format(self.eigen_band))
            elif self.read_dos:
                if (name == "i" and self.state["i"] == "efermi") or (name == "r" and self.state["set"]):
                    self.read_val = True
                elif name == "set" and "comment" in attributes:
                    comment = attributes["comment"]
                    self.state["set"] = comment
                    if self.state['partial']:
                        if comment.startswith("ion"):
                            self.pdos_ion = int(comment.split(" ")[1])
                        elif comment.startswith("spin"):
                            self.pdos_spin = Spin.up if self.state["set"] == "spin 1" else Spin.down

        if name == "calculation":
            self.step_count += 1
            self.scdata = []
            self.read_calculation = True
        elif name == "scstep":
            self.scstep = {}
        elif name == 'structure':
            self.latticestr = StringIO.StringIO()
            self.latticerec = StringIO.StringIO()
            self.posstr = StringIO.StringIO()
            self.read_structure = True
        elif name == 'varray' and (self.state['varray'] == "forces" or self.state['varray'] == "stress"):
            self.posstr = StringIO.StringIO()

    def characters(self, data):
        if self.read_val:
            self.val.write(data)
        if self.read_lattice:
            self.latticestr.write(data)
        elif self.read_positions:
            self.posstr.write(data)
        elif self.read_rec_lattice:
            self.latticerec.write(data)

    def _read_input(self, name):
        state = self.state
        if name == "i":
            if self.state['incar']:
                self.incar[self.incar_param] = parse_parameters(self.param_type, self.val.getvalue().strip())
            elif state['parameters']:
                self.parameters[self.incar_param] = parse_parameters(self.param_type, self.val.getvalue().strip())
            elif state['generator'] and state["i"] == "version":
                self.vasp_version = self.val.getvalue().strip()
            self.incar_param = None
        elif name == "set":
            if state['array'] == "atoms":
                self.atomic_symbols = self.atomic_symbols[::2]
                self.atomic_symbols = [sym if sym != "X" else "Xe" for sym in self.atomic_symbols]
            elif state['array'] == "atomtypes":
                self.potcar_symbols = self.potcar_symbols[4::5]
                self.input_read = True
        elif name == "c":
            if state['array'] == "atoms":
                self.atomic_symbols.append(self.val.getvalue().strip())
            elif state['array'] == "atomtypes":
                self.potcar_symbols.append(self.val.getvalue().strip())
        elif name == "v":
            if state['incar']:
                self.incar[self.incar_param] = parse_v_parameters(self.param_type, self.val.getvalue().strip(), self.filename, self.incar_param)
                self.incar_param = None
            elif state['parameters']:
                self.parameters[self.incar_param] = parse_v_parameters(self.param_type, self.val.getvalue().strip(), self.filename, self.incar_param)
            elif state['kpoints']:
                if state['varray'] == 'kpointlist':
                    self.actual_kpoints.append([float(x) for x in re.split("\s+", self.val.getvalue().strip())])
                if state['varray'] == 'weights':
                    self.actual_kpoints_weights.append(float(self.val.getvalue()))
                if state['v'] == "divisions":
                    self.kpoints.kpts = [[int(x) for x in re.split("\s+", self.val.getvalue().strip())]]
                elif state['v'] == "usershift":
                    self.kpoints.kpts_shift = [float(x) for x in re.split("\s+", self.val.getvalue().strip())]
                elif state['v'] == "genvec1" or state['v'] == "genvec2" or state['v'] == "genvec3" or state['v'] == "shift":
                    setattr(self.kpoints, state['v'], [float(x) for x in re.split("\s+", self.val.getvalue().strip())])

    def _read_calc(self, name):
        state = self.state
        if name == "i" and state['scstep']:
            self.scstep[state['i']] = float(self.val.getvalue())
        elif name == 'scstep':
            self.scdata.append(self.scstep)
            logger.debug('Finished reading scstep...')
        elif name == 'varray' and state['varray'] == "forces":
            self.forces = np.array([float(x) for x in re.split("\s+", self.posstr.getvalue().strip())])
            self.forces.shape = (len(self.atomic_symbols), 3)
            self.read_positions = False
        elif name == 'varray' and state['varray'] == "stress":
            self.stress = np.array([float(x) for x in re.split("\s+", self.posstr.getvalue().strip())])
            self.stress.shape = (3, 3)
            self.read_positions = False
        elif name == "calculation":
            self.ionic_steps.append({'electronic_steps':self.scdata, 'structure':self.structures[-1], 'forces': self.forces, 'stress':self.stress})
            self.read_calculation = False

    def _read_structure(self, name):
        if name == "v":
            self.read_positions = False
            self.read_lattice = False
            self.read_rec_lattice = False
        elif name == "structure":
            self.lattice = np.array([float(x) for x in re.split("\s+", self.latticestr.getvalue().strip())])
            self.lattice.shape = (3, 3)
            self.pos = np.array([float(x) for x in re.split("\s+", self.posstr.getvalue().strip())])
            self.pos.shape = (len(self.atomic_symbols), 3)
            self.structures.append(Structure(self.lattice, self.atomic_symbols, self.pos))
            self.lattice_rec = Lattice([float(x) for x in re.split("\s+", self.latticerec.getvalue().strip())])
            self.read_structure = False
            self.read_positions = False
            self.read_lattice = False
            self.read_rec_lattice = False

    def _read_dos(self, name):
        state = self.state
        try:
            if name == "i" and state["i"] == "efermi":
                self.efermi = float(self.val.getvalue().strip())
            elif name == "r" and state["total"]  and str(state["set"]).startswith("spin"):
                tok = re.split("\s+", self.val.getvalue().strip())
                self.dos_energies_val.append(float(tok[0]))
                self.dos_val.append(float(tok[1]))
                self.idos_val.append(float(tok[2]))
            elif name == "r" and state["partial"]  and str(state["set"]).startswith("spin"):
                tok = re.split("\s+", self.val.getvalue().strip())
                self.raw_data.append([float(i) for i in tok[1:]])
            elif name == "set" and state["total"] and str(state["set"]).startswith("spin"):
                spin = Spin.up if state["set"] == "spin 1" else Spin.down
                self.tdos[spin] = self.dos_val
                self.idos[spin] = self.dos_val
                self.dos_energies = self.dos_energies_val
                self.dos_energies_val = []
                self.dos_val = []
                self.idos_val = []
            elif name == "set" and state["partial"] and str(state["set"]).startswith("spin"):
                spin = Spin.up if state["set"] == "spin 1" else Spin.down
                self.norbitals = len(self.raw_data[0])
                for i in xrange(self.norbitals):
                    self.pdos[(self.pdos_ion, i, spin)] = [row[i] for row in self.raw_data]
                self.raw_data = []
            elif name == "partial":
                all_pdos = []
                natom = len(self.atomic_symbols)
                for iatom in xrange(1, natom + 1):
                    all_pdos.append(defaultdict())
                    for iorbital in xrange(self.norbitals):
                        updos = self.pdos[(iatom, iorbital, Spin.up)]
                        downdos = None if (iatom, iorbital, Spin.down) not in self.pdos else self.pdos[(iatom, iorbital, Spin.down)]
                        if downdos:
                            all_pdos[-1][Orbital.from_vasp_index(iorbital)] = {Spin.up:updos, Spin.down:downdos}
                        else:
                            all_pdos[-1][Orbital.from_vasp_index(iorbital)] = {Spin.up:updos}
                self.pdos = all_pdos
            elif name == "total":
                self.tdos = Dos(self.efermi, self.dos_energies, self.tdos)
                self.idos = Dos(self.efermi, self.dos_energies, self.idos)
            elif name == "dos":
                self.read_dos = False
        except:
            self.dos_has_errors = True

    def _read_eigen(self, name):
        state = self.state
        if name == "r" and str(state["set"]).startswith("kpoint"):
            tok = re.split("\s+", self.val.getvalue().strip())
            self.raw_data.append([float(i) for i in tok])
        elif name == "set" and str(state["set"]).startswith("kpoint"):
            self.eigenvalues[(self.eigen_spin, self.eigen_kpoint - 1)] = self.raw_data
            self.raw_data = []
        elif name == "eigenvalues":
            logger.debug("Finished reading eigenvalues. No. eigen = {}".format(len(self.eigenvalues)))
            self.read_eigen = False

    def _read_projected_eigen(self, name):
        state = self.state
        if name == "r" and str(state["set"]).startswith("band"):
            tok = re.split("\s+", self.val.getvalue().strip())
            self.raw_data.append({Orbital.from_vasp_index(i): float(val) for i, val in enumerate(tok)})
        elif name == "set" and str(state["set"]).startswith("band"):
            logger.debug("Processing projected eigenvalues for band {}, kpoint {}, spin {}.".format(self.eigen_band - 1, self.eigen_kpoint - 1, self.eigen_spin))
            for atom_ind, data in enumerate(self.raw_data):
                for orb, val in data.items():
                    self.projected_eigenvalues[(self.eigen_spin, self.eigen_kpoint - 1, self.eigen_band - 1, atom_ind, orb)] = val
            self.raw_data = []
        elif name == "projected":
            logger.debug("Finished reading projected eigenvalues. No. eigen = {}".format(len(self.eigenvalues)))
            self.read_projected_eigen = False

    def endElement(self, name):
        if not self.input_read:
            self._read_input(name)
        else:
            if self.read_structure:
                self._read_structure(name)
            elif self.read_dos:
                self._read_dos(name)
            elif self.read_eigen:
                self._read_eigen(name)
            elif self.read_projected_eigen:
                self._read_projected_eigen(name)
            elif self.read_calculation:
                self._read_calc(name)
        self.state[name] = False


def parse_parameters(val_type, val):
    """
    Helper function to convert a Vasprun parameter into the proper type.
    Boolean, int and float types are converted.
    
    Args:
        val_type : Value type parsed from vasprun.xml.
        val : Actual string value parsed for vasprun.xml.
    """
    if val_type == "logical":
        return (val == "T")
    elif val_type == "int":
        return int(val)
    elif val_type == "string":
        return val.strip()
    else:
        return float(val)


def parse_v_parameters(val_type, val, filename, param_name):
    """
    Helper function to convert a Vasprun array-type parameter into the proper type.
    Boolean, int and float types are converted.
    
    Args:
        val_type: 
            Value type parsed from vasprun.xml.
        val: 
            Actual string value parsed for vasprun.xml.
        filename: 
            Fullpath of vasprun.xml. Used for robust error handling.  E.g.,
            if vasprun.xml contains \*\*\* for some Incar parameters, the code
            will try to read from an INCAR file present in the same directory.
        param_name: 
            Name of parameter.
            
    Returns:
        Parsed value.
    """
    if val_type == "logical":
        val = [True if i == "T" else False for i in re.split("\s+", val)]
    elif val_type == "int":
        try:
            val = [int(i) for i in re.split("\s+", val)]
        except ValueError:
            # Fix for stupid error in vasprun sometimes which displays
            # LDAUL/J as 2****
            val = parse_from_incar(filename, param_name)
            if val == None:
                raise IOError("Error in parsing vasprun.xml")
    elif val_type == "string":
        val = [i for i in re.split("\s+", val)]
    else:
        try:
            val = [float(i) for i in re.split("\s+", val)]
        except ValueError:
            # Fix for stupid error in vasprun sometimes which displays
            # MAGMOM as 2****
            val = parse_from_incar(filename, param_name)
            if val == None:
                raise IOError("Error in parsing vasprun.xml")
    return val


def parse_from_incar(filename, key):
    """
    Helper function to parse a parameter from the INCAR.
    """
    dirname = os.path.dirname(filename)
    for f in os.listdir(dirname):
        if re.search("INCAR", f):
            warnings.warn("INCAR found. Using " + key + " from INCAR.")
            incar = Incar.from_file(os.path.join(dirname, f))
            if key in incar:
                return incar[key]
            else:
                return None
    return None


class Outcar(object):
    """
    Parser for data in OUTCAR that is not available in Vasprun.xml

    Note, this class works a bit differently than most of the other VaspObjects,
    since the OUTCAR can be very different depending on which "type of run"
    performed.

    Creating the OUTCAR class with a filename reads "regular parameters" that
    are always present.
    
    .. attribute:: magnetization
    
        Magnetization on each ion as a tuple of dict, e.g.,
        ({'d': 0.0, 'p': 0.003, 's': 0.002, 'tot': 0.005}, ... )
        Note that this data is not always present.  LORBIT must be set to some
        other value than the default.
    
    .. attribute:: charge
        
        Charge on each ion as a tuple of dict, e.g.,
        ({'p': 0.154, 's': 0.078, 'd': 0.0, 'tot': 0.232}, ...)
        Note that this data is not always present.  LORBIT must be set to some
        other value than the default.
    
    .. attribute:: is_stopped
    
        True if OUTCAR is from a stopped run (using STOPCAR, see Vasp Manual).
    
    .. attribute:: run_stats
        
        Various useful run stats as a dict including 'System time (sec)', 
        'Total CPU time used (sec)', 'Elapsed time (sec)',
        'Maximum memory used (kb)', 'Average memory used (kb)',
        'User time (sec)'.
    
    One can then call a specific reader depending on the type of run being
    performed. These are currently: read_igpar(), read_lepsilon() and
    read_lcalcpol().

    See the documentation of those methods for more documentation.
    
    Authors: Rickard Armiento, Shyue Ping Ong
    """
    def __init__(self, filename):
        self.filename = filename
        self.is_stopped = False
        with file_open_zip_aware(filename, "r") as f:
            lines = f.readlines()

        read_charge = False
        read_mag = False
        charge = []
        mag = []
        header = []
        run_stats = {}
        total_mag = None
        nelect = None
        efermi = None
        for line in lines:
            clean = line.strip()
            if clean == "total charge":
                read_charge = True
                charge = []
            elif clean == "magnetization (x)":
                read_mag = True
                mag = []
            elif read_charge or read_mag:
                if clean.startswith("# of ion"):
                    header = re.split("\s{2,}", line.strip())
                elif clean.startswith("tot"):
                    read_charge = False
                    read_mag = False
                else:
                    m = re.match("\s*(\d+)\s+(([\d\.\-]+)\s+)+", clean)
                    if m:
                        to_append = charge if read_charge else mag
                        data = re.findall("[\d\.\-]+", clean)
                        to_append.append({header[i]:float(data[i]) for i in xrange(1, len(header))})
            elif line.find('soft stop encountered!  aborting job') != -1:
                self.is_stopped = True
            elif re.search("\((sec|kb)\):", line):
                tok = line.strip().split(":")
                run_stats[tok[0].strip()] = float(tok[1].strip())
            elif re.search("E-fermi\s+:", clean):
                m = re.search("E-fermi\s*:\s*(\S+)", clean)
                efermi = float(m.group(1))
            elif re.search("number of electron\s+\S+", clean):
                m = re.search("number of electron\s+(\S+)\s+magnetization\s+(\S+)", clean)
                nelect = float(m.group(1))
                total_mag = float(m.group(2))

        self.run_stats = run_stats
        self.magnetization = tuple(mag)
        self.charge = tuple(charge)
        self.efermi = efermi
        self.nelect = nelect
        self.total_mag = total_mag

    def read_igpar(self):
        """ 
        Renders accessible:
            er_ev = e<r>_ev (dictionary with Spin.up/Spin.down as keys)
            er_bp = e<r>_bp (dictionary with Spin.up/Spin.down as keys)
            er_ev_tot = spin up + spin down summed
            er_bp_tot = spin up + spin down summed
            p_elc = spin up + spin down summed
            p_ion = spin up + spin down summed
        
        (See VASP section 'LBERRY,  IGPAR,  NPPSTR,  DIPOL tags' for info on
        what these are).
        """

        # variables to be filled
        self.er_ev = {}  #  will  be  dict (Spin.up/down) of array(3*float)
        self.er_bp = {}  #  will  be  dics (Spin.up/down) of array(3*float)
        self.er_ev_tot = None # will be array(3*float)
        self.er_bp_tot = None # will be array(3*float)
        self.p_elec = None
        self.p_ion = None
        try:
            search = []
            # Nonspin cases
            def er_ev(results, match): results.er_ev[Spin.up] = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))]) / 2.0; results.er_ev[Spin.down] = results.er_ev[Spin.up]; results.context = 2
            search.append(['^ *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', None, er_ev])

            def er_bp(results, match): results.er_bp[Spin.up] = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))]) / 2.0; results.er_bp[Spin.down] = results.er_bp[Spin.up]
            search.append(['^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', lambda results, line: results.context == 2, er_bp])

            # Spin cases
            def er_ev_up(results, match): results.er_ev[Spin.up] = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))]); results.context = Spin.up
            search.append(['^.*Spin component 1 *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', None, er_ev_up])

            def er_bp_up(results, match): results.er_bp[Spin.up] = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))])
            search.append(['^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', lambda results, line: results.context == Spin.up, er_bp_up])

            def er_ev_dn(results, match): results.er_ev[Spin.down] = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))]); results.context = Spin.down
            search.append(['^.*Spin component 2 *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', None, er_ev_dn])

            def er_bp_dn(results, match): results.er_bp[Spin.down] = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))])
            search.append(['^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', lambda results, line: results.context == Spin.down, er_bp_dn])

            # Always present spin/non-spin
            def p_elc(results, match): results.p_elc = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))])
            search.append(['^.*Total electronic dipole moment: *p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', None, p_elc])

            def p_ion(results, match): results.p_ion = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))])
            search.append(['^.*ionic dipole moment: *p\[ion\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', None, p_ion])

            self.context = None
            self.er_ev = {Spin.up: None, Spin.down: None}
            self.er_bp = {Spin.up: None, Spin.down: None}

            micro_pyawk(self.filename, search, self)

            if self.er_ev[Spin.up] != None and self.er_ev[Spin.down] != None:
                self.er_ev_tot = self.er_ev[Spin.up] + self.er_ev[Spin.down]

            if self.er_bp[Spin.up] != None and self.er_bp[Spin.down] != None:
                self.er_bp_tot = self.er_bp[Spin.up] + self.er_bp[Spin.down]

        except:
            self.er_ev_tot = None
            self.er_bp_tot = None
            raise Exception("IGPAR OUTCAR could not be parsed.")

    def read_lepsilon(self):
        # variables to be filled
        try:
            search = []

            def dielectric_section_start(results, match): results.dielectric_index = -1;
            search.append(['MACROSCOPIC STATIC DIELECTRIC TENSOR', None, dielectric_section_start])

            def dielectric_section_start2(results, match): results.dielectric_index = 0
            search.append(['-------------------------------------', lambda results, line: results.dielectric_index == -1, dielectric_section_start2])

            def dielectric_data(results, match): results.dielectric_tensor[results.dielectric_index, :] = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))]); results.dielectric_index += 1
            search.append(['^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$', lambda results, line: results.dielectric_index >= 0, dielectric_data])

            def dielectric_section_stop(results, match): results.dielectric_index = None
            search.append(['-------------------------------------', lambda results, line: results.dielectric_index >= 1, dielectric_section_stop])

            self.dielectric_index = None
            self.dielectric_tensor = np.zeros((3, 3))

            def piezo_section_start(results, match): results.piezo_index = 0;
            search.append(['PIEZOELECTRIC TENSOR  for field in x, y, z        \(e  Angst\)', None, piezo_section_start])

            def piezo_data(results, match): results.piezo_tensor[results.piezo_index, :] = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3)), float(match.group(4)), float(match.group(5)), float(match.group(6))]); results.piezo_index += 1
            search.append(['^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)*$', lambda results, line: results.piezo_index >= 0, piezo_data])

            def piezo_section_stop(results, match): results.piezo_index = None
            search.append(['-------------------------------------', lambda results, line: results.piezo_index >= 1, piezo_section_stop])

            self.piezo_index = None
            self.piezo_tensor = np.zeros((3, 6))

            def born_section_start(results, match): results.born_ion = -1;
            search.append(['BORN EFFECTIVE CHARGES \(in e, cummulative output\)', None, born_section_start])

            def born_ion(results, match): results.born_ion = int(match.group(1)) - 1; results.born[results.born_ion] = np.zeros((3, 3));
            search.append(['ion +([0-9]+)', lambda results, line: results.born_ion != None, born_ion])

            def born_data(results, match): results.born[results.born_ion][int(match.group(1)) - 1, :] = np.array([float(match.group(2)), float(match.group(3)), float(match.group(4))]);
            search.append(['^ *([1-3]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)$', lambda results, line: results.born_ion >= 0, born_data])

            def born_section_stop(results, match): results.born_index = None
            search.append(['-------------------------------------', lambda results, line: results.born_ion >= 1, born_section_stop])

            self.born_ion = None
            self.born = {}

            micro_pyawk(self.filename, search, self)

        except:
            raise Exception("LEPSILON OUTCAR could not be parsed.")

    def read_lcalcpol(self):
        # variables to be filled
        self.p_elec = None
        self.p_ion = None
        try:
            search = []

            # Always present spin/non-spin
            def p_elc(results, match): results.p_elc = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))])
            search.append(['^.*Total electronic dipole moment: *p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', None, p_elc])

            def p_ion(results, match): results.p_ion = np.array([float(match.group(1)), float(match.group(2)), float(match.group(3))])
            search.append(['^.*Ionic dipole moment: *p\[ion\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)', None, p_ion])

            micro_pyawk(self.filename, search, self)

        except:
            raise Exception("CLACLCPOL OUTCAR could not be parsed.")

    @property
    def to_dict(self):
        d = {}
        d['module'] = self.__class__.__module__
        d['class'] = self.__class__.__name__
        d['efermi'] = self.efermi
        d['run_stats'] = self.run_stats
        d['magnetization'] = self.magnetization
        d['charge'] = self.charge
        d['total_magnetization'] = self.total_mag
        d['nelect'] = self.nelect
        d['is_stopped'] = self.is_stopped
        return d

class VolumetricData(object):
    """
    Simple volumetric object for reading LOCPOT and CHGCAR type files.
    
    .. attribute:: structure
    
        Structure associated with the Volumetric Data object
    
    ..attribute:: is_spin_polarized
    
        True if run is spin polarized
    
    ..attribute:: dim
        
        Tuple of dimensions of volumetric grid in each direction (nx, ny, nz).
    
    ..attribute:: data
        
        Actual data as a dict of {string: np.array}. The string are "total"
        and "diff", in accordance to the output format of vasp LOCPOT and
        CHGCAR files where the total spin density is written first, followed
        by the difference spin density.
    
    .. attribute:: ngridpts
        
        Total number of grid points in volumetric data.
    """
    def __init__(self, structure, data, distance_matrix=None):
        """
        Typically, this constructor is not used directly and the static
        from_file constructor is used. This constructor is designed to allow
        summation and other operations between VolumetricData objects.
        
        Args:
            structure:
                Structure associated with the volumetric data
            data:
                Actual volumetric data.
            distance_matrix:
                A pre-computed distance matrix if available. Useful so pass
                distance_matrices between sums, shortcircuiting an otherwise
                expensive operation.
        """
        self.structure = structure
        self.is_spin_polarized = len(data) == 2
        self.dim = data['total'].shape
        self.data = data
        self.ngridpts = self.dim[0] * self.dim[1] * self.dim[2]
        #lazy init the spin data since this is not always needed.
        self._spin_data = {}
        self._distance_matrix = {} if not distance_matrix else distance_matrix

    @property
    def spin_data(self):
        """
        The data decomposed into actual spin data as {spin: data}.
        Essentially, this provides the actual Spin.up and Spin.down data
        instead of the total and diff.  Note that by definition, a
        non-spin-polarized run would have Spin.up data == Spin.down data.
        """
        if not self._spin_data:
            spin_data = dict()
            spin_data[Spin.up] = 0.5 * (self.data['total'] + self.data.get('diff', 0))
            spin_data[Spin.down] = 0.5 * (self.data['total'] - self.data.get('diff', 0))
            self._spin_data = spin_data
        return self._spin_data

    def get_axis_grid(self, ind):
        """
        Returns the grid for a particular axis.
        
        Args:
            ind:
                Axis index.
        """
        ng = self.dim
        num_pts = ng[ind]
        lengths = self.structure.lattice.abc
        return [i / num_pts * lengths[ind] for i in xrange(num_pts)]

    def __add__(self, other):
        return self.linear_add(other, 1.0)

    def __sub__(self, other):
        return self.linear_add(other, -1.0)

    def linear_add(self, other, scale_factor=1.0):
        '''
        Method to do a linear sum of volumetric objects. Used by + and -
        operators as well. Returns a VolumetricData object containing the
        linear sum.
        
        Args:
            other:
                Another VolumetricData object
            scale_factor:
                Factor to scale the other data by.
                
        Returns:
            VolumetricData corresponding to self + scale_factor * other.
        '''
        if self.structure != other.structure:
            raise ValueError("Adding or subtraction operations can only be performed for volumetric data with the exact same structure.")
        #To add checks
        data = {}
        for k in self.data.keys():
            data[k] = self.data[k] + scale_factor * other.data[k]
        return VolumetricData(self.structure, data, self._distance_matrix)

    @staticmethod
    def parse_file(filename):
        """
        Convenience method to parse a generic volumetric data file in the vasp
        like format. Used by subclasses for parsing file.
        
        Args:
            filename:
                Path of file to parse
        
        Returns:
            (poscar, data)
        """

        with file_open_zip_aware(filename) as f:
            contents = f.read()
            (poscar_string, grid_data) = re.split("^\s*$", contents, flags=re.MULTILINE)
            poscar = Poscar.from_string(poscar_string)

            grid_data = grid_data.strip()
            ind = grid_data.find("\n")
            dim_line = grid_data[0:ind].strip()
            a = [int(i) for i in dim_line.split()]

            pieces = grid_data[ind:].strip().split(dim_line)

            spinpolarized = (len(pieces) == 2)

            def parse_data(piece):
                pot = np.zeros(a)
                count = 0
                data = re.sub("\n", " ", piece).split()
                for (z, y, x) in itertools.product(xrange(a[2]), xrange(a[1]), xrange(a[0])):
                    pot[x, y, z] = float(data[count])
                    count += 1
                return pot

            totalpot = parse_data(pieces[0].strip())
            if spinpolarized:
                diffpot = parse_data(pieces[1].strip())
                data = {"total":totalpot, "diff":diffpot}
            else:
                data = {"total":totalpot}

            return (poscar, data)

    def write_file(self, file_name, vasp4_compatible=False):
        """
        Write the VolumetricData object to a vasp compatible file.
            
        Args:
            file_name:
                the path to a file
            vasp4_compatible:
                True if the format is vasp4 compatible
        """

        f = open(file_name, 'w')
        p = Poscar(self.structure)
        f.write(p.get_string(vasp4_compatible=vasp4_compatible) + "\n")
        a = self.dim
        f.write("\n")
        def write_spin(data_type):
            lines = []
            count = 0
            f.write("{} {} {}\n".format(a[0], a[1], a[2]))
            for (k, j, i) in itertools.product(xrange(a[2]), xrange(a[1]), xrange(a[0])):
                lines.append('%0.11e' % self.data[data_type][i, j, k])
                count += 1
                if count % 5 == 0:
                    f.write(''.join(lines) + "\n")
                    lines = []
                else:
                    lines.append(" ")
            f.write(''.join(lines) + "\n")
        write_spin("total")
        if self.is_spin_polarized:
            f.write("\n")
            write_spin("diff")
        f.close()

    def _calculate_distance_matrix(self, ind):
        site = self.structure[ind]
        a = self.dim
        distances = dict()
        for (x, y, z) in itertools.product(xrange(a[0]), xrange(a[1]), xrange(a[2])):
            pt = np.array([x / a[0], y / a[1] , z / a[2]])
            distances[(x, y, z)] = site.distance_and_image_from_frac_coords(pt)[0]
        self._distance_matrix[ind] = distances

    def get_integrated_diff(self, ind, radius):
        """
        Get integrated difference of atom index ind up to radius.
        
        Args:
            ind:
                Index of atom.
            radius:
                Radius of integration.
            
        Returns:
            Differential integrated charge.
        """
        #Shortcircuit the data since by definition, this will be zero for non
        #spin-polarized runs.
        if not self.is_spin_polarized:
            return 0
        if ind not in self._distance_matrix:
            self._calculate_distance_matrix(ind)
        a = self.dim
        intchg = 0
        for (x, y, z) in itertools.product(xrange(a[0]), xrange(a[1]), xrange(a[2])):
            if self._distance_matrix[ind][(x, y, z)] < radius:
                intchg += self.data["diff"][x, y, z]
        return intchg / self.ngridpts

    def get_average_along_axis(self, ind):
        """
        Get the averaged total of the volumetric data a certain axis direction.
        For example, useful for visualizing Hartree Potentials from a LOCPOT
        fike.
        
        Args:
            ind : Index of axis.
            
        Returns:
            Average total along axis
        """
        m = self.data["total"]

        ng = self.dim
        avg = []
        for i in xrange(ng[ind]):
            subtotal = 0
            for j in xrange(ng[(ind + 1) % 3]):
                for k in xrange(ng[(ind + 2) % 3]):
                    if ind == 0:
                        subtotal += m[i, j, k]
                    if ind == 1:
                        subtotal += m[k, i, j]
                    if ind == 2:
                        subtotal += m[j, k, i]
            avg.append(subtotal)
        avg = np.array(avg) / ng[(ind + 1) % 3] / ng[(ind + 2) % 3]
        return avg


class Locpot(VolumetricData):
    """
    Simple object for reading a LOCPOT file.
    """

    def __init__(self, poscar, data):
        """
        Args:
            poscar:
                Poscar object containing structure.
            data:
                Actual data.
        """
        VolumetricData.__init__(self, poscar.structure, data)
        self.name = poscar.comment

    @staticmethod
    def from_file(filename):
        (poscar, data) = VolumetricData.parse_file(filename)
        return Locpot(poscar, data)



class Chgcar(VolumetricData):
    """
    Simple object for reading a CHGCAR file.
    """

    def __init__(self, poscar, data):
        """
        Args:
            poscar:
                Poscar object containing structure.
            data:
                Actual data.
        """
        VolumetricData.__init__(self, poscar.structure, data)
        self.poscar = poscar
        self.name = poscar.comment
        self._distance_matrix = {}

    @staticmethod
    def from_file(filename):
        (poscar, data) = VolumetricData.parse_file(filename)
        return Chgcar(poscar, data)


class Procar(object):
    """
    Object for reading a PROCAR file
    """
    def __init__(self, filename):
        """
        Args:
            filename:
                Name of file containing PROCAR.
        """
        #create and return data object containing the information of a PROCAR type file
        self.name = ""
        self.data = dict()
        self._read_file(filename)

    def get_d_occupation(self, atomNo):
        row = self.data[atomNo]
        return sum(row[4:9])

    def _read_file(self, filename):
        reader = file_open_zip_aware(filename, "r")
        lines = clean_lines(reader.readlines())
        reader.close()
        self.name = lines[0]
        kpointexpr = re.compile("^\s*k-point\s+(\d+).*weight = ([0-9\.]+)")
        expr = re.compile('^\s*([0-9]+)\s+')
        dataexpr = re.compile('[\.0-9]+')
        currentKpoint = 0
        weight = 0
        for l in lines:
            if kpointexpr.match(l):
                m = kpointexpr.match(l)
                currentKpoint = int(m.group(1))
                weight = float(m.group(2))
                if currentKpoint == 1:
                    self.data = dict()
            if expr.match(l):
                linedata = dataexpr.findall(l)
                linefloatdata = map(float, linedata)
                index = int(linefloatdata.pop(0))
                if index in self.data:
                    self.data[index] = self.data[index] + np.array(linefloatdata) * weight
                else:
                    self.data[index] = np.array(linefloatdata) * weight


class Oszicar(object):
    """
    A basic parser for an OSZICAR output from VASP.  In general, while the
    OSZICAR is useful for a quick look at the output from a VASP run, we
    recommend that you use the Vasprun parser instead, which gives far richer
    information about a run.
    
    .. attribute:: electronic_steps
    
            All electronic steps as a list of list of dict. e.g., 
            [[{'rms': 160.0, 'E': 4507.24605593, 'dE': 4507.2, 'N': 1,
            'deps': -17777.0, 'ncg': 16576}, ...], [....]
            where electronic_steps[index] refers the list of electronic steps
            in one ionic_step, electronic_steps[index][subindex] refers to a
            particular electronic step at subindex in ionic step at index. The
            dict of properties depends on the type of VASP run, but in general,
            "E", "dE" and "rms" should be present in almost all runs.
    
    .. attribute:: ionic_steps:
    
            All ionic_steps as a list of dict, e.g.,
            [{'dE': -526.36, 'E0': -526.36024, 'mag': 0.0, 'F': -526.36024},
            ...]
            This is the typical output from VASP at the end of each ionic step.  
    """

    def __init__(self, filename):
        """
        Args:
            filename:
                Filename of file to parse
        """
        electronic_steps = []
        ionic_steps = []
        ionic_pattern = re.compile("(\d+)\s+F=\s*([\d\-\.E\+]+)\s+E0=\s*([\d\-\.E\+]+)\s+d\s*E\s*=\s*([\d\-\.E\+]+)\s+mag=\s*([\d\-\.E\+]+)")
        electronic_pattern = re.compile("\s*\w+\s*:(.*)")
        def smart_convert(header, num):
            if header == "N" or header == "ncg":
                return int(num)
            return float(num)
        header = []
        with open(filename, 'r') as fid:
            for line in fid.readlines():
                m = electronic_pattern.match(line)
                if m:
                    toks = re.split("\s+", m.group(1).strip())
                    data = {header[i]:smart_convert(header[i], toks[i]) for i in xrange(len(toks))}
                    if toks[0] == '1':
                        electronic_steps.append([data])
                    else:
                        electronic_steps[-1].append(data)
                elif ionic_pattern.match(line.strip()):
                    m = ionic_pattern.match(line.strip())
                    ionic_steps.append({'F':float(m.group(2)), 'E0':float(m.group(3)), 'dE':float(m.group(4)), 'mag':float(m.group(5))})
                elif re.match("^\s*N\s+E\s*", line):
                    header = re.split("\s+", line.strip().replace("d eps", "deps"))
        self.electronic_steps = electronic_steps
        self.ionic_steps = ionic_steps

    @property
    def all_energies(self):
        """
        Compilation of all energies from all electronic steps and ionic steps
        as a tuple of list of energies, e.g.,
        ((4507.24605593, 143.824705755, -512.073149912, ...), ...)
        """
        all_energies = []
        for i in xrange(len(self.electronic_steps)):
            energies = [step['E'] for step in self.electronic_steps[i]]
            energies.append(self.ionic_steps[i]['F'])
            all_energies.append(tuple(energies))
        return tuple(all_energies)

    @property
    def final_energy(self):
        """
        Final energy from run.
        """
        return self.ionic_steps[-1]['F']

    @property
    def to_dict(self):
        return {'electronic_steps':self.electronic_steps,
                'ionic_steps':self.ionic_steps}


class VaspParserError(Exception):
    '''
    Exception class for Structure.
    Raised when the structure has problems, e.g., atoms that are too close.
    '''

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return "VaspParserError : " + self.msg


def get_band_structure_from_vasp_multiple_branches(dir_name, efermi=None):
    """
    this method is used to get band structure info from a VASP directory. It
    takes into account that the run can be divided in several branches named
    "branch_x". If the run has not been divided in branches the method will
    turn to parsing vasprun.xml directly.

    The method returns None is there's a parsing error

    Args:
        dir_name:
            Directory containing all bandstructure runs.
        efermi:
            Efermi for bandstructure.
    
    Returns:
        (bandstructure_up, bandstructure_down)
    """
    #ToDo: Add better error handling!!!
    if os.path.exists(os.path.join(dir_name, "branch_0")):
        #get all branch dir names
        branch_dir_names = [os.path.abspath(d) for d in glob.glob("{i}/branch_*".format(i=dir_name)) if os.path.isdir(d)]

        #sort by the directory name (e.g, branch_10)
        sort_by = lambda x: int(x.split('_')[-1])
        sorted_branch_dir_names = sorted(branch_dir_names, key=sort_by)

        # populate branches with Bandstructure instances
        branches = []

        for dir_name in sorted_branch_dir_names:
            xml_file = os.path.join(dir_name, 'vasprun.xml')
            if os.path.exists(xml_file):
                run = Vasprun(xml_file)
                branches.append(run.get_band_structure(efermi=efermi))
            else:
                # It might be better to throw an exception
                warnings.warn("Skipping {d}. Unable to find {f}".format(d=dir_name, f=xml_file))
        return get_reconstructed_band_structure(branches, efermi)
    else:
        xml_file = os.path.join(dir_name, 'vasprun.xml')
        #Better handling of Errors
        if os.path.exists(xml_file):
            return Vasprun(xml_file).get_band_structure(kpoints_filename=None, efermi=efermi)
        else:
            return None


def get_adjusted_fermi_level(run_static, band_structure):
    """
    When running a band structure computations the fermi level needs to be taken
    from the static run that gave the charge density used for the non-self
    consistent band structure run. Sometimes this fermi level is however a
    little too low because of the mismatch between the uniform grid used in
    the static run and the band structure k-points (e.g., the VBM is on Gamma
    and the Gamma point is not in the uniform mesh). Here we use a procedure
    consisting in looking for energy levels higher than the static fermi level
    (but lower than the LUMO) if any of these levels make the band structure
    appears insulating and not metallic anymore, we keep this adjusted fermi
    level. This procedure has shown to detect correctly most insulators.
    
    Args:
        run_static: 
            a Vasprun object for the static run
        run_bandstructure: 
            a band_structure object
            
    Returns:
        a new adjusted fermi level
    """
    #make a working copy of band_structure
    bs_working = BandStructureSymmLine.from_dict(band_structure.to_dict)
    if bs_working.is_metal():
        e = run_static.efermi
        while e < run_static.eigenvalue_band_properties[1]:
            e += 0.01
            bs_working._efermi = e
            if not bs_working.is_metal():
                return e
    return run_static.efermi

