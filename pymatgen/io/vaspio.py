#!/usr/bin/env python

"""
Classes for reading/manipulating/writing VASP files.
All major VASP input files, plus the more common Vasp output files are available.
"""

from __future__ import division

__author__="Shyue Ping Ong, Rickard Armiento, Vincent L Chevrier"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="Sep 23, 2011"

import os
import re
import math
import itertools
import warnings
import xml.sax.handler
import StringIO
from collections import defaultdict
import ConfigParser

import numpy as np
from numpy.linalg import det

import pymatgen.command_line.aconvasp_caller
import pymatgen.core.electronic_structure
from pymatgen.core.design_patterns import Enum
from pymatgen.io.io_abc import VaspInput
from pymatgen.util.string_utils import str_aligned, str_delimited
from pymatgen.util.io_utils import file_open_zip_aware, clean_lines, micro_pyawk, clean_json
from pymatgen.core.structure import Structure, Composition
from pymatgen.core.periodic_table import Element
from pymatgen.core.electronic_structure import CompleteDos, Dos, PDos, Spin, Orbital, Bandstructure
from pymatgen.core.lattice import Lattice


coord_pattern = re.compile("^\s*([\d+\.\-Ee]+)\s+([\d+\.\-Ee]+)\s+([\d+\.\-Ee]+)")    

class Poscar(VaspInput):
    """
    Object for representing the data in a POSCAR or CONTCAR file.
    Please note that this current implementation does not parse selective
    dynamcics POSCAR files.
    """
    
    def __init__(self, struct, comment = None, selective_dynamics = None):
        """        
        Arguments:
            struct:
                Structure object. See pymatgen.core.structure.Structure.
            comment:
                Optional comment line for POSCAR. Defaults to unit cell formula of structure. Defaults to None.
            selective_dynamics:
                Nx3 2D array of boolean values for selective dynamics, where N is number of sites. Defaults to None.
        """

        if struct.is_ordered:
            self._struct = struct
            self._site_symbols = []
            self._natoms = []
            syms = [site.specie.symbol for site in struct]
            for (s, data) in itertools.groupby(syms):
                self._site_symbols.append(s)
                self._natoms.append(len(tuple(data)))
            self._true_names = True
            if selective_dynamics:
                self.set_selective_dynamics(selective_dynamics)
            else:
                self._selective_dynamics = None
            self.comment = struct.formula if comment == None else comment    
        else:
            raise ValueError("Structure with partial occupancies cannot be converted into POSCAR!")

    def set_selective_dynamics(self, selective_dynamics):
        if selective_dynamics and len(selective_dynamics) != len(self._struct):
            raise ValueError("Selective dynamics array must be same length as the structure.")      
        self._selective_dynamics = selective_dynamics
    
    @property
    def selective_dynamics(self):
        return self._selective_dynamics

    @property
    def site_symbols(self):
        """
        Symbols for each site in POSCAR.
        """
        return self._site_symbols
    
    @property
    def struct(self):
        """
        Structure associated with the Poscar file.
        """
        return self._struct

    @staticmethod
    def from_file(filename):
        """
        Reads a Poscar from a file.
        The code will try its best to determine the elements in the POSCAR in the following order:
        1. Ideally, if the input file is Vasp5-like and contains element symbols in the 6th line, the code will use that.
        2. Failing (1), the code will check if a symbol is provided at the end of each coordinate.
        3. Failing (i) and (ii), the code will try to check if a POTCAR is in the same directory as the POSCAR and use elements from that.
        If all else fails, the code will just assign the first n elements in increasing atomic number, where n is the number of species, 
        to the Poscar.  For example, H, He, Li, ....  This will ensure at least a unique element is assigned to each site and any analysis 
        that does not require specific elemental properties should work fine.
        
        Arguments:
            filename - file name containing Poscar data.
        
        Returns:
            Poscar object.
        """
        
        dirname = os.path.dirname(os.path.abspath(filename))
        names = None
        for f in os.listdir(dirname):
            if re.search("POTCAR.*",f):
                try:
                    potcar = Potcar.from_file(os.path.join(dirname, f))
                    names = [sym.split("_")[0] for sym in potcar.symbols]
                except:
                    names = None    
        with file_open_zip_aware(filename, "r") as f:
            return Poscar.from_string(f.read(), names)
    
    @staticmethod
    def from_string(data, default_names = None):
        """
        Reads a Poscar from a string.
        The code will try its best to determine the elements in the POSCAR in the following order:
        1) Ideally, if the input file is Vasp5-like and contains element symbols in the 6th line, the code will use that.
        2) Failing (i), the code will check if a symbol is provided at the end of each coordinate.
        If all else fails, the code will just assign the first n elements in increasing atomic number, where n is the number of species, 
        to the Poscar.  For example, H, He, Li, ....  This will ensure at least a unique element is assigned to each site and any analysis 
        that does not require specific elemental properties should work fine.
        
        Arguments:
            data:
                string containing Poscar data.
                
        Returns:
            Poscar object.
        """
        
        lines = tuple(clean_lines(data.split("\n"), False))

        comment = lines[0]
        scale = float(lines[1])
        lattice = np.array([[float(s) for s in line.split()] for line in lines[2:5]])
        if scale < 0:
            vol = abs(det(lattice))
            lattice = (-scale / vol) ** (1/3) * lattice
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
        #This is in line with Vasp's parsing order that the POTCAR specified is the default used.
        if default_names:
            atomic_symbols = list()
            for i in xrange(len(natoms)):
                atomic_symbols.extend([default_names[i]] * natoms[i])
        elif not vasp5_symbols:
            ind = 3 if not sdynamics else 6
            try: #check if names are appended at the end of the POSCAR coordinates
                atomic_symbols = [l.split()[ind] for l in lines[ipos + 1:ipos + 1 + nsites]]
                [Element(sym) for sym in atomic_symbols] #Ensure symbols are valid elements
            except:
                #Defaulting to false names.
                atomic_symbols = list()
                for i in xrange(len(natoms)):
                    sym = Element.from_Z(i+1).symbol
                    atomic_symbols.extend([sym] * natoms[i])
                warnings.warn("Elements in POSCAR cannot be determined. Defaulting to false names, " + " ".join(atomic_symbols)+".")

        # read the atomic coordinates
        coords = []
        selective_dynamics = list() if sdynamics else None
        for i in xrange(nsites):
            toks = lines[ipos + 1 + i].split()
            coords.append([float(s) for s in toks[:3]])
            if sdynamics:
                selective_dynamics.append([True if tok.upper()[0] == 'T' else False for tok in toks[3:6]])
            
        struct = Structure(lattice, atomic_symbols, coords, False, False, cart)

        return Poscar(struct,comment, selective_dynamics)
    
    def get_string(self, direct = True, vasp4_compatible = False):
        """
        Returns a string to be written as a POSCAR file. By default, site symbols are written
        which means compatibilty is for vasp >= 5.
        
        Arguments:
            direct:
                Whether coordinates are output in direct or cartesian. Defaults to True.
            vasp4_compatible:
                Set to True to omit site symbols on 6th line to maintain backward vasp 4.x compatibility. Defaults to False.
                
        Returns:
            String representation of POSCAR.
        """
        lines = [self.comment]
        lines.append("1.0")
        lines.append(str(self._struct.lattice))
        if self._true_names == True and not vasp4_compatible:
            lines.append(" ".join(self._site_symbols))
        lines.append(" ".join([str(x) for x in self._natoms]))
        if self._selective_dynamics:
            lines.append("Selective dynamics")
        lines.append('direct' if direct else 'cartesian')         
        
        for i in xrange(len(self._struct)):
            site = self._struct[i]
            coords = site.frac_coords if direct else site.coords
            line = "%.6f %.6f %.6f" % (coords[0], coords[1], coords[2])
            if self._selective_dynamics:
                sd = ['T' if j else 'F' for j in self._selective_dynamics[i]]
                line += " %s %s %s" % (sd[0], sd[1], sd[2])
            line += " " + site.species_string
            lines.append(line)

        return "\n".join(lines)    

    def __str__(self):
        """
        String representation of Poscar file.
        """
        return self.get_string()

    def set_site_symbols(self, symbols):
        self._site_symbols = symbols
        self._true_names = True
        #update the Structure as well
        elements = list()
        for i in range(len(symbols)):
            elements.extend([symbols[i]]*self._natoms[i])
        self._struct = Structure(self._struct.lattice, elements, self._struct.frac_coords)

    def write_file(self, filename):
        with open(filename, 'w') as f:
            f.write(str(self) + "\n")

"""**Non-exhaustive** list of valid INCAR tags"""
VALID_INCAR_TAGS = ("NGX", "NGY", "NGZ", "NGXF", "NGYF", "NGZF", "NBANDS", "NBLK", "SYSTEM", "NWRITE", "ENCUT", "ENAUG",
"PREC", "ISPIN", "MAGMOM", "ISTART", "ICHARG", "INIWAV", "NELM", "NELMIN", "NELMDL", "EDIFF", "EDIFFG", "NSW", "NBLOCK",
"KBLOCK", "IBRION", "NFREE", "POTIM", "ISIF", "PSTRESS", "IWAVPR", "ISYM", "SYMPREC", "LCORR", "TEBEG", "TEEND", "SMASS",
"NPACO", "APACO", "POMASS", "ZVAL", "RWIGS", "LORBIT", "NELECT", "NUPDOWN", "EMIN", "EMAX", "NEDOS", "ISMEAR", "SIGMA", 
"FERWE", "FERDO", "SMEARINGS", "LREAL", "ROPT", "GGA", "VOSKOWN", "LASPH", "ALGO", "IALGO", "LDIAG", "NSIM", "IMIX", "INIMIX", 
"MAXMIX", "AMIX", "BMIX", "AMIX_MAG", "BMIX_MAG", "AMIN", "MIXPRE", "WC", "WEIMIN", "EBREAK", "DEPER", "TIME", "LWAVE", "LCHARG",
"LVTOT", "LELF", "NPAR", "LPLANE","LASYNC", "LSCALAPACK", "LSCALU", "ISPIND", "HFSCREEN", "LHFCALC", "ENCUTFOCK", "NKRED", "LMAXMIX",
"PRECFOCK", "AEXX", "AGGAX", "AGGAC", "ALDAC", "LMAXFOCK", "LMAXFOCKAE", "LTHOMAS", "NKREDX", "NKREDY", "NKREDZ", "EVENONLY", "ODDONLY", "LDAU", "LDAUJ", "LDAUL", "LDAUPRINT", "LDAUTYPE", "LDAUU", "LPEAD", "LCALCPOL", "LCALCEPS", "LEFG", "EFIELD_PEAD", "LNONCOLLINEAR",
"LSORBIT", "IDIPOL", "DIPOL", "LMONO", "LDIPOL", "EPSILON", "EFIELD", "LBERRY", "IGPAR", "NPPSTR", "IPEAD", "I_CONSTRAINED_M", "LAMBDA", "M_CONSTR",
"IMAGES", "SPRING", "LOPTICS", "CSHIFT", "LNABLA", "LEPSILON", "LRPA", "NOMEGA", "NOMEGAR", "LSPECTRAL", "OMEGAMAX", "OMEGATL", "ENCUTGW",
"ENCUTGWSOFT", "ODDONLYGW", "EVENONLYGW", "LSELFENERGY", 'LRHFATM', 'METAGGA', 'LMAXTAU', 'LCOMPAT','ENMAX', 'LMAXPAW', 'LSPIRAL', 'LZEROZ',
'LMETAGGA','ENINI', 'NRMM', 'MREMOVE', 'ADDGRID', 'EFERMI', 'LPARD', 'LSCAAWARE', 'IDIOT', 'LMUSIC', 'LREAL_COMPAT', 'GGA_COMPAT', 'ICORELEVEL', 'LHFONE',
'LRHFCALC', 'LMODELHF', 'ENCUT4O', 'EXXOEP', 'FOURORBIT', 'HFALPHA', 'ALDAX', 'SHIFTRED', 'NMAXFOCKAE', 'HFSCREENC', 'MODEL_GW', 'MODEL_EPS0', 'MODEL_ALPHA',
'LVEL', 'SAXIS', 'QSPIRAL', 'STM', 'KINTER', 'ORBITALMAG', 'LMAGBLOCH', 'LCHIMAG', 'LGAUGE', 'MAGATOM', 'MAGDIPOL', 'AVECCONST', 'LTCTE', 'LTETE',
'L2ORDER', 'LGWLF', 'ENCUTLF', 'LMAXMP2', 'SCISSOR', 'NBANDSGW', 'NBANDSLF', 'DIM', 'ANTIRES', 'LUSEW', 'OMEGAGRID', 'SELFENERGY', 'NKREDLFX', 'NKREDLFY',
'NKREDLFZ', 'MAXMEM', 'TELESCOPE', 'LCRITICAL_MEM', 'GGA2',
'TURBO', 'QUAD_EFG','IRESTART','NREBOOT','NMIN','EREF','KSPACING','KGAMMA','LSUBROT','SCALEE','LVHAR','LORBITALREAL','DARWINR','DARWINV','LFOCKAEDFT','NUCIND','MAGPOS','LNICSALL','LADDER','LHARTREE','IBSE','NBANDSO','NBANDSV','OPTEMAX')

class Incar(dict, VaspInput):
    """
    INCAR object for reading and writing INCAR files
    essentially consists of a dictionary with some helper functions
    """
    
    def __init__(self,params = dict()):
        """
        Creates an Incar object.
        
        Arguments:
            params - A set of input parameters as a dictionary.
        """
        super(Incar,self).__init__()
        self.update(params)

    def __setitem__(self, key, val):
        """
        Add parameter-val pair to Incar.  Warns if parameter is not in list of valid INCAR tags.
        Also cleans the parameter and val by stripping leading and trailing white spaces.
        """
        if key.strip().upper() not in VALID_INCAR_TAGS:
            warnings.warn(key.strip() + " not in VALID_INCAR_TAGS")
        super(Incar,self).__setitem__(key.strip(), Incar.proc_val(key.strip(), val.strip()) if isinstance(val,basestring) else val)

    def get_string(self, sort_keys = False, pretty = False):
        """
        Returns a string representation of the INCAR.  The reason why this method is
        different from the __str__ method is to provide options for pretty printing.
        
        Arguments:
            sort_keys:
                Set to True to sort the INCAR parameters alphabetically. Defaults to False.
            pretty:
                Set to True for pretty aligned output. Defaults to False.
        """
        keys = self.keys()
        if sort_keys:
            keys = sorted(keys)
        lines = []
        for k in keys:
            if isinstance(self[k], list):
                lines.append([k," ".join([str(i) for i in self[k]])])
            else:
                lines.append([k,self[k]])
        
        if pretty:
            return str_aligned(lines)
        else:
            return str_delimited(lines, None," = ")
    
    def __str__(self):
        return self.get_string(sort_keys = True, pretty = False)

    def write_file(self, filename):
        """
        Write Incar to a file.
        
        Arguments:
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
        Static helper method to convert INCAR parameters to proper types, e.g. integers, floats, lists, etc.
        
        Arguments:
            key:
                INCAR parameter key
            val:
                Actual value of INCAR parameter.
        """
        list_type_keys = ('LDAUU', 'LDAUL', 'LDAUJ', 'LDAUTYPE','MAGMOM')
        boolean_type_keys = ('LDAU', 'LWAVE','LSCALU', 'LCHARG', 'LPLANE', 'LHFCALC')
        float_type_keys = ("EDIFF", "SIGMA", 'TIME','ENCUTFOCK', 'HFSCREEN')
        int_type_keys = ('NSW', 'NELMIN', 'ISIF', 'IBRION', "ISPIN", "ICHARG", "NELM", 
                         "ISMEAR", "NPAR", "LDAUPRINT", 'LMAXMIX', 'ENCUT', 'NSIM',
                         'NKRED', 'NUPDOWN', 'ISPIND')

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)
        try:
            if key in list_type_keys:
                output = list()
                toks = re.split("\s+", val)
                
                for tok in toks:
                    m = re.match("(\d+)\*([\d\.\-\+]+)", tok)
                    if m:
                        output.extend([smart_int_or_float(m.group(2))] * int(m.group(1)))
                    else:
                        output.append(smart_int_or_float(tok))
                return output
            if key in boolean_type_keys:
                m = re.search("^\W+([TtFf])", val)
                if m:
                    if m.group(1) == "T" or m.group(1) == "t":
                        return True
                    else:
                        return False
                raise ValueError(key + " should be a boolean type!")
            
            if key in float_type_keys:
                return float(val)
            
            if key in int_type_keys:
                return int(val)
            
        except:
            return val.capitalize()
              
        return val.capitalize()
            
    def diff(self, other):
        """
        Diff function for Incar.  Compares two Incars and indicates which parameters are the same and which are not.
        Useful for checking whether two runs were done using the same parameters.
        
        Args:
            other : The other Incar object to compare to.
        
        Returns:
            Dict of the following format - {'Same' : parameters_that_are_the_same, 'Different': parameters_that_are_different}
            Note that the parameters are return as full dictionaries of values. E.g. {'ISIF':3}
        """
        similar_param = {}
        different_param = {}
        for k1,v1 in self.items():
            if k1 not in other:
                different_param[k1] = {"INCAR1": v1, "INCAR2": 'Default'}
            elif v1 != other[k1]:
                different_param[k1] = {"INCAR1": v1, "INCAR2": other[k1]}
            else:
                similar_param[k1] = v1
        for k2,v2 in other.items():
            if k2 not in similar_param and k2 not in different_param:
                if k2 not in self:
                    different_param[k2] = {"INCAR1": 'Default', "INCAR2": v2}
        return {'Same' : similar_param, 'Different': different_param}
                
    def __add__(self, other):
        """
        Add all the values of another INCAR object to this object
        Facilitates the use of "standard" INCARs
        """
        params = {k:v for k,v in self.items()}
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
    supported_modes = Enum(("Gamma", "Monkhorst", "Automatic", "Line_mode", "Cartesian", "Reciprocal"))
    
    def __init__(self, comment = "Default gamma", num_kpts = 0, style = supported_modes.Gamma, 
                 kpts = [[1,1,1]], kpts_shift = (0,0,0), 
                 kpts_weights=None, coord_type = None, labels = None, 
                 tet_number = 0, tet_weight = 0, tet_connections = None):
        """
        Highly flexible constructor for Kpoints object.  The flexibility comes at the cost of usability and in 
        general, it is recommended that you use the default constructor only if you know exactly what you are doing
        and requires the flexibility.  For most usage cases, the three automatic schemes can be constructed far 
        more easily using the convenience static constructors (automatic, gamma_automatic, monkhorst_automatic) and it 
        is recommended that you use those.
        
        Arguments:
            comment:
                String comment for Kpoints
            num_kpts:
                Following VASP method of defining the KPOINTS file, this parameter is the number of kpoints specified.
                If set to 0 (or negative), VASP automatically generates the KPOINTS.
            style:
                Style for generating KPOINTS.  Use one of the Kpoints.supported_modes enum types.
            kpts:
                2D array of kpoints.  Even when only a single specification is required, e.g. in the automatic scheme, 
                the kpts should still be specified as a 2D array. e.g., [[20]] or [[2,2,2]].
            kpts_shift:
                Shift for Kpoints.
            kpts_weights:
                Optional weights for kpoints.  For explicit kpoints.
            coord_type:
                In line-mode, this variable specifies whether the Kpoints were given in Cartesian or Reciprocal coordinates
            labels:
                In line-mode, this should provide a list of labels for each kpts.
            tet_number:
                For explicit kpoints, specifies the number of tetrahedrons for the tetrahedron method.
            tet_weight:
                For explicit kpoints, specifies the weight for each tetrahedron for the tetrahedron method.
            tet_connections:
                For explicit kpoints, specifies the connections of the tetrahedrons for the tetrahedron method.
                Format is a list of tuples, [ (sym_weight, [tet_vertices]), ...]           
        
        The default behavior of the constructor is for a Gamma centered, 1x1x1 KPOINTS with no shift.
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
        gamma centered Monkhorst-Pack grids and the number of subdivisions along
        each reciprocal lattice vector determined by the scheme in the VASP manual.
        
        Arguments:
            subdivisions:
                 Parameter determining number of subdivisions along each reciprocal lattice vector
                 
        Returns:
            Kpoints object
        """
        return Kpoints("Fully automatic kpoint scheme", 0, style = Kpoints.supported_modes.Automatic, kpts=[[subdivisions]])
    
    @staticmethod
    def gamma_automatic(kpts = (1,1,1), shift = (0,0,0)):
        """
        Convenient static constructor for an automatic Gamma centered Kpoint grid.
        
        Arguments:
            kpts:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors. Defaults to (1,1,1)
            shift:
                Shift to be applied to the kpoints. Defaults to (0,0,0)
                
        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0, Kpoints.supported_modes.Gamma, kpts = [kpts], kpts_shift = shift)

    @staticmethod
    def monkhorst_automatic(kpts = (2,2,2), shift = (0,0,0)):
        """
        Convenient static constructor for an automatic Monkhorst pack Kpoint grid.
        
        Arguments:
            kpts:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors. Defaults to (2,2,2)
            shift:
                Shift to be applied to the kpoints. Defaults to (0,0,0)
                
        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0, Kpoints.supported_modes.Monkhorst, kpts = [kpts], kpts_shift = shift)
    
    @staticmethod
    def monkhorst_automatic_density(structure, kppa):
        """
            The approach we took here is to use aconvasp to get the kpoint divisions
        """
        div=pymatgen.command_line.aconvasp_caller.get_num_division_kpoints(structure, kppa)
        #check if hexagonal!!!!
        if(Kpoints._is_hexagonal(structure)==True):
            return Kpoints.gamma_automatic(div, shift = (0,0,0))
    
    @staticmethod
    def _is_hexagonal(structure):
        proto1=[math.pi/2.0,math.pi/2.0,2.0*math.pi/3.0]
        proto2=[math.pi/3.0,math.pi/2.0,math.pi/2.0]
        angCopy=[c for c in structure.lattice.angles]
        angCopy.sort()
        anglesCheck=False
        for i in range(len(angCopy)):
            if(math.fabs(angCopy[i]-proto1[i])<0.1):
                anglesCheck=True
            if(math.fabs(angCopy[i]-proto2[i])<0.1):
                anglesCheck=True

        num_equals=0
        if(math.fabs(structure.lattice.abc[0]-structure.lattice.abc[1])<0.001):
            num_equals=num_equals+1
        if(math.fabs(structure.lattice.abc[0]-structure.lattice.abc[2])<0.001):
            num_equals=num_equals+1
        if(math.fabs(structure.lattice.abc[1]-structure.lattice.abc[2])<0.001):
            num_equals=num_equals+1

        if(num_equals > 0 and anglesCheck):
            return True

        return False


    @staticmethod
    def from_file(filename):
        """
        Reads a Kpoints object from a KPOINTS file.
        
        Arguments:
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
        
        #Automatic gamma and Monk KPOINTS, with optional shift
        if style == "g" or style == "m":
            kpts = [int(x) for x in lines[3].split()]
            kpts_shift = (0,0,0)
            if len(lines)>4 and coord_pattern.match(lines[4]):
                try:
                    kpts_shift = [int(x) for x in lines[4].split()]
                except:
                    pass
            return Kpoints.gamma_automatic(kpts, kpts_shift) if style == "g" else Kpoints.monkhorst_automatic(kpts, kpts_shift) 

        #Automatic kpoints with basis
        if num_kpts <= 0:
            style = Kpoints.supported_modes.Cartesian if style in "ck" else Kpoints.supported_modes.Reciprocal
            kpts = [[float(x) for x in lines[i].split()] for i in xrange(3,6)]
            kpts_shift = [float(x) for x in lines[6].split()]
            return Kpoints(comment = comment, num_kpts = num_kpts, style = style, kpts = kpts, kpts_shift = kpts_shift)

        #Line-mode KPOINTS, usually used with band structures
        if style == 'l':
            coord_type = 'Cartesian' if lines[3].lower()[0] in 'ck' else 'Reciprocal'
            style = Kpoints.supported_modes.Line_mode
            kpts = []
            labels = []
            patt = re.compile('([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s*!\s*(.*)')
            for i in range(4,len(lines)):
                line = lines[i]
                m = patt.match(line)
                if m:
                    kpts.append([float(m.group(1)),float(m.group(2)),float(m.group(3))])
                    labels.append(m.group(4).strip())
            return Kpoints(comment = comment, num_kpts = num_kpts, style = style, 
                 kpts = kpts, coord_type = coord_type, labels = labels)
        
        #Assume explicit KPOINTS if all else fails.
        style = Kpoints.supported_modes.Cartesian if style == "ck" else Kpoints.supported_modes.Reciprocal
        kpts = []
        kpts_weights = []
        tet_number = 0
        tet_weight = 0
        tet_connections = None

        for i in xrange(3, 3 + num_kpts):
            toks = re.split("\s+", lines[i])
            kpts.append([float(toks[0]),float(toks[1]),float(toks[2])])
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
                    tet_connections.append((int(toks[0]), [int(toks[j]) for j in xrange(1,5)]))
        except:
            pass
        
        return Kpoints(comment = comment, num_kpts = num_kpts, style = style, 
                 kpts = kpts, kpts_weights = kpts_weights, 
                 tet_number = tet_number, tet_weight = tet_weight, tet_connections = tet_connections)
        
    def write_file(self, filename):
        """
        Write Kpoints to a file.
        
        Arguments:
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
        if self.num_kpts <= 0 and tuple(self.kpts_shift) != (0,0,0):
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
        return d 
    
class PotcarSingle(VaspInput):
    """
    Object for a **single** POTCAR.
    The builder assumes the complete string is the POTCAR
    contains the complete untouched data in "data" as a string and
    a dictionary of keywords
    """
    def __init__(self, data):
        """
        Expects a complete and single potcar file as a string in "data"
        """
        self.data = data # raw POTCAR as a string
        keyValPairs = re.compile(r";*\s*(.+?)\s*=\s*([^;\n]+)\s*", re.M).findall(data)
        self.keywords = dict(keyValPairs) # all key = val found in the POTCAR as a dictionary all keys and vals are strings

    def __str__(self):
        return self.data

    def write_file(self, filename):
        writer = open(filename, 'w')
        writer.write(self.__str__() + "\n")
        writer.close()

    @property
    def symbol(self):
        """
        Full name of POTCAR, e.g., Fe_pv
        """
        return self.keywords['TITEL'].split(" ")[1].strip()        
    
    @property
    def element(self):
        """
        Attempt to return the atomic symbol based on the VRHFIN keyword
        """
        return self.keywords['VRHFIN'].split(":")[0].strip()

    @property
    def atomic_no(self):
        """
        Attempt to return the atomic number based on the VRHFIN keyword
        """
        return Element(self.element).Z

class Potcar(list,VaspInput):
    """
    Object for reading and writing POTCAR files for
    calculations.
    """
    functional_dir = {'PBE':'POT_GGA_PAW_PBE', 'LDA':'POT_LDA_PAW', 'PW91':'POT_GGA_PAW_PW91'}
    
    def __init__(self, symbols = None):
        if symbols != None:
            self.set_symbols(symbols)

    @staticmethod
    def from_file(filename):
        with file_open_zip_aware(filename, "r") as reader:
            fData = reader.read()
        potcar = Potcar()
        potcar_strings = re.compile(r"\n{0,1}\s*(.*?End of Dataset)", re.S).findall(fData)
        for p in potcar_strings:
            potcar.append(PotcarSingle(p))
        return potcar

    def __str__(self):
        return "".join([str(potcar) for potcar in self])    #line break not used because there is already one at the end of str(potcar) and it causes VASP issues

    def write_file(self, filename):
        """
        Write Potcar to a file.
        
        Arguments:
            filename:
                filename to write to.
        """
        with open(filename, 'w') as f:
            f.write(self.__str__() + "\n")
        
    @property
    def symbols(self):
        """
        Get the atomic symbols of all the atoms in the POTCAR file
        """
        return [p.symbol for p in self]

    def set_symbols(self, elements, functional = 'PBE'):
        module_dir = os.path.dirname(pymatgen.__file__)
        config = ConfigParser.SafeConfigParser()
        config.readfp(open(os.path.join(module_dir, "pymatgen.cfg")))
        VASP_PSP_DIR = os.path.join(config.get('VASP', 'pspdir'), Potcar.functional_dir[functional])
        del self[:]
        for el in elements:
            with file_open_zip_aware(os.path.join(VASP_PSP_DIR, "POTCAR." + el + ".gz"), 'rb') as f:
                self.append(PotcarSingle(f.read()))

class Vasprun(object):
    """
    Vastly improved sax-based parser for vasprun.xml files.
    Speedup over Dom is at least 2x for smallish files (~1Mb) to orders of magnitude for larger files (~10Mb).
    All data is stored as attributes, which are delegated to the VasprunHandler object.
    
    Attributes:
    
        **Vasp results**
        
        ionic_steps: 
            All ionic steps in the run as a list of {'structure': structure at end of run, 'electronic_steps' : {All electronic step data in vasprun file}, 'stresses' : stress matrix}
        structures: 
            List of Structure objects for the structure at each ionic step.
        tdos: 
            Total dos calculated at the end of run.
        idos: 
            Integrated dos calculated at the end of run.
        pdos: 
            List of list of PDos objects. Access as pdos[atomindex][orbitalindex]
        efermi: 
            Fermi energy
        eigenvalues: 
            Final eigenvalues as a dict of {(kpoint index, Spin.up):[[eigenvalue, occu]]}. 
            This representation is probably not ideal, but since this is not used anywhere else for now, I leave it as such.
            Future developers who need to work with this should refactored the object into a sensible structure.
        
        **Vasp inputs**
        
        incar:
            Incar object for parameters specified in INCAR file.
        parameters:
            Incar object with parameters that vasp actually used, including all defaults.
        kpoints:
            Kpoints object for KPOINTS specified in run.
        actual_kpoints:
            List of actual kpoints, e.g., [[0.25, 0.125, 0.08333333], [-0.25, 0.125, 0.08333333], [0.25, 0.375, 0.08333333], ....]
        actual_kpoints_weights:
            List of kpoint weights, E.g., [0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, ....]
        atomic_symbols:
            List of atomic symbols, e.g., [u'Li', u'Fe', u'Fe', u'P', u'P', u'P']
        potcar_symbols:
            List of POTCAR symbols. E.g., [u'PAW_PBE Li 17Jan2003', u'PAW_PBE Fe 06Sep2000', ..]
    
        **Convenience attributes**
        
        final_energy: 
            Final energy from the run.
        final_structure: 
            Final relaxed structure.
        initial_structure: 
            Initial input structure.
        complete_dos: 
            CompleteDos object containing both the total and projected Dos from the run.
    
    Author: Shyue Ping Ong
    """
    supported_properties = ['lattice_rec','vasp_version', 'incar', 'parameters', 'potcar_symbols', 'atomic_symbols', 'kpoints', 'actual_kpoints', 'structures',
                            'actual_kpoints_weights', 'dos_energies', 'eigenvalues', 'tdos', 'idos', 'pdos', 'efermi', 'ionic_steps', 'dos_has_errors']
    
    def __init__(self, filename):
        self._filename = filename    
        with file_open_zip_aware(filename) as f:
            self._handler = VasprunHandler(filename)
            self._parser = xml.sax.parse(f, self._handler)
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
        return CompleteDos(self.final_structure, self.tdos, self.pdos)
    
    def get_band_structure(self, kpoints_filename = None):
        """
        Returns the band structure.
        
        Arguments:
            kpoints_filename:
                Full path of the KPOINTS file from which the band structure is generated.
                If none is provided, the code will try to intelligently determine the appropriate
                KPOINTS file by substituting the filename of the vasprun.xml with KPOINTS.
                The latter is the default behavior.
                
        Returns:
            Band structure object
        """
        if not kpoints_filename:
            kpoints_filename = self._filename.replace('vasprun.xml', 'KPOINTS')
        if not os.path.exists(kpoints_filename):
            raise VaspParserError('KPOINTS file needed to obtain band structure.')
        k = Kpoints.from_file(kpoints_filename)
        labels_dict = dict(zip(k.labels, k.kpts))
        lattice_rec = Lattice(self.lattice_rec)
        #make the labels_dict to work with cartesian
        if k.coord_type == "Reciprocal":
            for c in labels_dict:
                labels_dict[c]=lattice_rec.get_cartesian_coords(2 * math.pi * np.array(labels_dict[c]))
        kpoints = [lattice_rec.get_cartesian_coords(2 * math.pi * np.array(k)) for k in self.actual_kpoints]
        eig = self.eigenvalues
        eigenvals=[]
        for i in range(len(eig[(1,Spin.up)])):
            eigenvals.append({'energy':[eig[(j+1, Spin.up)][i][0] for j in range(len(kpoints))], 
                              'occup': [eig[(j+1, Spin.up)][i][1] for j in range(len(kpoints))]})
        return Bandstructure(kpoints,eigenvals,labels_dict, lattice_rec, self.final_structure)
    
    @property
    def eigenvalue_band_properties(self):
        """
        Returns band properties from the eigenvalues.
        Returns:
            (band gap, cbm, vbm, is_band_gap_direct)
        """
        vbm = - float('inf')
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
        json-friendly dict representation for Vasprun for transferring between different applications.
        """
        d = {}
        d['vasp_version'] = self.vasp_version
        d['has_vasp_completed'] = self.converged
        d['nsites'] = len(self.final_structure)
        d['unit_cell_formula'] = self.final_structure.composition.to_dict
        comp = self.final_structure.composition
        d['reduced_cell_formula'] = Composition.from_formula(comp.reduced_formula).to_dict
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
                d['hubbards'] = { symbols[i] : us[i] - js[i] for i in xrange(len(symbols))}
            elif sum(us) == 0 and sum(js) == 0:
                d['is_hubbard'] = False
                d['hubbards'] = {}
            else:
                raise VaspParserError("Length of U value parameters and atomic symbols are mismatched")
        else:
            d['hubbards'] = {}
        if d['is_hubbard']:
            d['run_type'] = "GGA+U"
        elif self.parameters.get('LHFCALC', False):
            d['run_type'] = "HF"
        else:
            d['run_type'] = "GGA"
                   
        d['input'] = {}
        d['input']['incar'] = {k:v for k,v in self.incar.items()}
        d['input']['crystal'] = self.initial_structure.to_dict
        d['input']['kpoints'] = self.kpoints.to_dict
        d['input']['kpoints']['actual_points'] = [{'abc':list(self.actual_kpoints[i]), 'weight':self.actual_kpoints_weights[i]} for i in xrange(len(self.actual_kpoints))] 
        d['input']['potcar'] = [s.split(" ")[1] for s in self.potcar_symbols]
        d['input']['parameters'] = {k:v for k,v in self.parameters.items()}
        d['input']['lattice_rec'] = self.lattice_rec.to_dict
        
        d['output'] = {}
        d['output']['ionic_steps'] = clean_json(self.ionic_steps)
        d['output']['final_energy'] = self.final_energy
        d['output']['final_energy_per_atom'] = self.final_energy / len(self.final_structure)
        d['output']['crystal'] = self.final_structure.to_dict
        d['output']['efermi'] = self.efermi
        #{(kpoint index, Spin.up):array(float)}
        
        d['output']['eigenvalues'] = {}
        for (index, spin), values in self.eigenvalues.items():
            if str(index) not in d['output']['eigenvalues']:
                d['output']['eigenvalues'][str(index)] = {str(spin):values}
            else:
                d['output']['eigenvalues'][str(index)][str(spin)] = values
        
        (gap, cbm, vbm, is_direct) = self.eigenvalue_band_properties
        d['output'].update(dict(bandgap = gap, cbm = cbm, vbm = vbm, is_gap_direct = is_direct))
        
        return d

 
class VasprunHandler(xml.sax.handler.ContentHandler):
    """
    Sax handler for vasprun.xml.
    Attributes are mirrored into Vasprun object.
    
    Author: Shyue Ping Ong
    """
        
    def __init__(self, filename):
        self.filename = filename
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
        self.eigenvalues = {}#  will  be  {(kpoint index, Spin.up):array(float)}
        self.tdos = {}
        self.idos = {}
        self.pdos = {}
        self.efermi = None 
        self.ionic_steps = [] # should be a list of dict
        self.structures = []
        self.lattice_rec = []
        
        self.input_read = False
        self.all_calculations_read = False
        self.read_structure = False
        self.read_rec_lattice = False
        self.read_calculation = False
        self.read_eigen = False
        self.read_dos= False
        self.in_efermi = False
        self.read_atoms = False
        self.read_lattice = False
        self.read_positions = False
        self.incar_param = None
        self.dos_energies_val = []
        self.dos_val = []
        self.idos_val = []
        self.raw_data = []
        self.dos_has_errors = False #will be set to true if there is an error parsing the Dos.
        self.state = defaultdict(bool)
        
    
    def in_all(self, xml_tags):
        return all([getattr(self, 'in_'+tag, None) for tag in xml_tags])
    
    def startElement(self, name, attributes):
        
        self.state[name] = True if 'name' not in attributes else attributes['name']
        self.read_val = False
        
        #Nested if loops makes reading much faster.
        if not self.input_read: #reading input parameters
            if (name == "i" or name == "v") and (self.state['incar'] or self.state['parameters']):
                self.incar_param = attributes['name']
                self.param_type = 'float' if 'type' not in attributes else attributes['type']
                self.read_val = True
            elif name == "v" and self.state['kpoints']:
                self.read_val = True
            elif name == "generation" and self.state['kpoints']:
                self.kpoints.comment   = "Kpoints from vasprun.xml"
                self.kpoints.num_kpts  = 0
                self.kpoints.style = attributes['param']
                self.kpoints.kpts = []
                self.kpoints.kpts_shift = [0,0,0]
            elif name == "c" and (self.state['array'] == "atoms" or self.state['array'] == "atomtypes"):
                self.read_val = True
            elif name == "i" and self.state['i'] == "version" and self.state['generator']:
                self.read_val = True
                
        else: #reading calculations and structures.
            if self.read_calculation:
                if name == "i" and self.state['scstep']:
                    self.read_val = True
                elif name == "v" and (self.state['varray'] == "forces" or self.state['varray'] == "stress"):
                    self.read_positions = True
            if self.read_structure:
                if name == "v" and self.state['varray'] == 'basis':
                    self.read_lattice = True
                elif name == "v" and self.state['varray'] == 'positions':
                    self.read_positions = True
                if name == 'v' and self.state['varray'] == 'rec_basis':
                    self.read_rec_lattice = True
            if name == "calculation":
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
            elif name == "eigenvalues":
                self.all_calculations_read = True
            if self.read_eigen:
                if name == "r" and self.state["set"]:
                    self.read_val = True
                elif name == "set" and "comment" in attributes:
                    comment = attributes["comment"]
                    self.state["set"] = comment
                    if comment.startswith("spin"):
                        self.eigen_spin = Spin.up if self.state["set"] == "spin 1" else Spin.down
                    if comment.startswith("kpoint"):
                        self.eigen_kpoint = int(comment.split(" ")[1])
            elif self.read_dos:
                if (name == "i" and self.state["i"] == "efermi") or (name == "r" and self.state["set"] ):
                    self.read_val = True
                elif name == "set" and "comment" in attributes:
                    comment = attributes["comment"]
                    self.state["set"] = comment
                    if self.state['partial']:
                        if comment.startswith("ion"):
                            self.pdos_ion = int(comment.split(" ")[1])
                        elif comment.startswith("spin"):
                            self.pdos_spin = Spin.up if self.state["set"] == "spin 1" else Spin.down
            elif name == "dos":
                self.dos_energies = None
                self.tdos = {}
                self.idos = {}
                self.pdos = {}
                self.efermi = None 
                self.read_dos = True
            elif name == "eigenvalues":
                self.eigenvalues = {}#  will  be  {(kpoint index, Spin.up):array(float)}
                self.read_eigen = True
            
        if self.read_val:
            self.val = StringIO.StringIO()

    def characters(self, data):
        if self.read_val:
            self.val.write(data)
        if self.read_lattice:
            self.latticestr.write(data)
        elif self.read_positions:
            self.posstr.write(data)
        elif self.read_rec_lattice:
            self.latticerec.write(data)

    
    #To correct for stupid vasp bug which names Xenon as X!!
    EL_MAPPINGS = {'X':'Xe'}
    
    def endElement(self, name):
        
        if not self.input_read:
            if name == "i":
                if self.state['incar']:
                    self.incar[self.incar_param] = parse_parameters(self.param_type, self.val.getvalue().strip())
                elif self.state['parameters']:
                    self.parameters[self.incar_param] = parse_parameters(self.param_type, self.val.getvalue().strip())
                elif self.state['generator'] and self.state["i"] == "version":
                    self.vasp_version = self.val.getvalue().strip()
                self.incar_param = None
            elif name == "set":
                if self.state['array'] == "atoms":
                    self.atomic_symbols = self.atomic_symbols[::2]
                    self.atomic_symbols = [sym if sym not in VasprunHandler.EL_MAPPINGS else VasprunHandler.EL_MAPPINGS[sym] for sym in self.atomic_symbols]
                elif self.state['array'] == "atomtypes":
                    self.potcar_symbols = self.potcar_symbols[4::5]
                    self.input_read = True
            elif name == "c":
                if self.state['array'] == "atoms":
                    self.atomic_symbols.append(self.val.getvalue().strip())
                elif self.state['array'] == "atomtypes":
                    self.potcar_symbols.append(self.val.getvalue().strip())
            elif name == "v":
                if self.state['incar']:
                    self.incar[self.incar_param] = _parse_v_parameters(self.param_type, self.val.getvalue().strip(), self.filename, self.incar_param)
                    self.incar_param = None
                elif self.state['parameters']:
                    self.parameters[self.incar_param] = _parse_v_parameters(self.param_type, self.val.getvalue().strip(), self.filename, self.incar_param)
                elif self.state['kpoints']:
                    if self.state['varray'] == 'kpointlist':
                        self.actual_kpoints.append([float(x) for x in re.split("\s+",self.val.getvalue().strip())])
                    if self.state['varray'] == 'weights':
                        self.actual_kpoints_weights.append(float(self.val.getvalue()))
                    if self.state['v'] == "divisions":
                        self.kpoints.kpts   = [[int(x) for x in re.split("\s+",self.val.getvalue().strip())]]
                    elif self.state['v'] == "usershift":
                        self.kpoints.kpts_shift   = [float(x) for x in re.split("\s+",self.val.getvalue().strip())]
                    elif self.state['v'] == "genvec1" or self.state['v'] == "genvec2" or self.state['v'] == "genvec3" or self.state['v'] == "shift":
                        setattr(self.kpoints, self.state['v'], [float(x) for x in re.split("\s+",self.val.getvalue().strip())])
                        
        else:
            if self.read_calculation:
                if name == "i" and self.state['scstep']:
                    self.scstep[self.state['i']] = float(self.val.getvalue())
                elif name == 'scstep':
                    self.scdata.append(self.scstep)
                elif name == 'varray' and self.state['varray'] == "forces":
                    self.forces = np.array([float(x) for x in re.split("\s+",self.posstr.getvalue().strip())])
                    self.forces.shape = (len(self.atomic_symbols), 3)
                elif name == 'varray' and self.state['varray'] == "stress":
                    self.stress = np.array([float(x) for x in re.split("\s+",self.posstr.getvalue().strip())])
                    self.stress.shape = (3, 3)
                elif name == "calculation":
                    self.ionic_steps.append({'electronic_steps':self.scdata, 'structure':self.structures[-1], 'forces': self.forces, 'stress':self.stress})
                    self.read_calculation = False
            if self.read_structure:
                if name == "v":
                    self.read_positions = False
                    self.read_lattice = False
                    self.read_lattice_rec = False
                    self.read_rec_lattice = False
                elif name == "structure":
                    self.lattice = np.array([float(x) for x in re.split("\s+",self.latticestr.getvalue().strip())])
                    self.lattice.shape = (3,3)
                    self.pos = np.array([float(x) for x in re.split("\s+",self.posstr.getvalue().strip())])
                    self.pos.shape = (len(self.atomic_symbols), 3)
                    self.structures.append(Structure(self.lattice, self.atomic_symbols, self.pos))
                    self.lattice_rec = Lattice([float(x) for x in re.split("\s+",self.latticerec.getvalue().strip())])
                    self.read_structure = False
            elif self.read_dos:
                try:
                    if name == "i" and self.state["i"] == "efermi":
                        self.efermi = float(self.val.getvalue().strip())
                    elif name == "r" and self.state["total"]  and str(self.state["set"]).startswith("spin"):
                        tok = re.split("\s+", self.val.getvalue().strip())
                        self.dos_energies_val.append(float(tok[0]))
                        self.dos_val.append(float(tok[1]))
                        self.idos_val.append(float(tok[2]))
                    elif name == "r" and self.state["partial"]  and str(self.state["set"]).startswith("spin"):
                        tok = re.split("\s+", self.val.getvalue().strip())
                        self.raw_data.append([float(i) for i in tok[1:]])
                    elif name == "set" and self.state["total"] and str(self.state["set"]).startswith("spin"):
                        spin = Spin.up if self.state["set"] == "spin 1" else Spin.down
                        self.tdos[spin] = self.dos_val
                        self.idos[spin] = self.dos_val
                        self.dos_energies = self.dos_energies_val
                        self.dos_energies_val = []
                        self.dos_val = []
                        self.idos_val = []
                    elif name == "set" and self.state["partial"] and str(self.state["set"]).startswith("spin"):
                        spin = Spin.up if self.state["set"] == "spin 1" else Spin.down
                        self.norbitals = len(self.raw_data[0])
                        for i in xrange(self.norbitals):
                            self.pdos[(self.pdos_ion, i, spin)] = [row[i] for row in self.raw_data]
                        self.raw_data = []
                    elif name == "partial":
                        all_pdos = []
                        natom = len(self.atomic_symbols)
                        for iatom in xrange(1,natom+1):
                            all_pdos.append(list())
                            for iorbital in xrange(self.norbitals):
                                updos = self.pdos[(iatom, iorbital, Spin.up)]
                                downdos = None if (iatom, iorbital, Spin.down) not in self.pdos else self.pdos[(iatom, iorbital, Spin.down)]
                                if downdos:
                                    all_pdos[-1].append(PDos(self.efermi, self.dos_energies, {Spin.up:updos, Spin.down:downdos}, Orbital.from_vasp_index(iorbital)))
                                else:
                                    all_pdos[-1].append(PDos(self.efermi, self.dos_energies, {Spin.up:updos}, Orbital.from_vasp_index(iorbital)))
                        self.pdos = all_pdos
                    elif name == "total":
                        self.tdos = Dos(self.efermi, self.dos_energies, self.tdos)
                        self.idos = Dos(self.efermi, self.dos_energies, self.idos)
                    elif name == "dos":
                        self.read_dos = False
                except:
                    self.dos_has_errors = True
            elif self.read_eigen:
                if name == "r" and str(self.state["set"]).startswith("kpoint"):
                    tok = re.split("\s+", self.val.getvalue().strip())
                    self.raw_data.append([float(i) for i in tok])
                elif name == "set" and str(self.state["set"]).startswith("kpoint"):
                    self.eigenvalues[(self.eigen_kpoint, self.eigen_spin)] = self.raw_data
                    self.raw_data = []
                elif name == "eigenvalues":
                    self.read_eigen = False
                
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

def _parse_v_parameters(val_type, val, filename, param_name):
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
            if vasprun.xml contains *** for some Incar parameters, the code will try
            to read from an INCAR file present in the same directory.
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
            val = _parse_from_incar(filename, param_name)
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
            val = _parse_from_incar(filename, param_name)
            if val == None:
                raise IOError("Error in parsing vasprun.xml")
    return val

def _parse_from_incar(filename, key):
    """
    Helper function to parse a parameter from the INCAR.
    """
    dirname = os.path.dirname(filename)
    for f in os.listdir(dirname):
        if re.search("INCAR",f):
            warnings.warn("INCAR found. Using "+key+" from INCAR.")
            incar = Incar.from_file(os.path.join(dirname, f))
            if key in incar:
                return incar[key]
            else:
                return None
    return None


class Outcar(object):
    """
    Parser for data in OUTCAR that is not available in Vasprun.xml

    Note, this class works a bit differently than most of the other VaspObjects, since the OUTCAR can
    be very different depending on which "type of run" performed.

    Creating the OUTCAR class with a filename reads "regular parameters" that are always present.
    
    Default attributes:
        magnetization:
            Magnetization on each ion as a tuple of dict, e.g., ({'d': 0.0, 'p': 0.003, 's': 0.002, 'tot': 0.005}, ... )
            Note that this data is not always present.  LORBIT must be set to some other value than the default.
        charge:
            Charge on each ion as a tuple of dict, e.g., ({'p': 0.154, 's': 0.078, 'd': 0.0, 'tot': 0.232}, ...)
            Note that this data is not always present.  LORBIT must be set to some other value than the default.
        is_stopped:
            True if OUTCAR is from a stopped run (using STOPCAR, see Vasp Manual).
        run_stats:
            Various useful run stats as a dict including 'System time (sec)', 'Total CPU time used (sec)'
            'Elapsed time (sec)', 'Maximum memory used (kb)', 'Average memory used (kb)', 'User time (sec)'.
                
            
    One can then call a specific reader depending on the type of run being perfromed. These are currently:
       read_igpar()
       read_lepsilon()
       read_lcalcpol()

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
                        to_append.append({header[i]:float(data[i]) for i in xrange(1,len(header))})
            elif line.find('soft stop encountered!  aborting job') != -1:
                self.is_stopped = True
            elif re.search("\((sec|kb)\):", line):
                tok = line.strip().split(":")
                run_stats[tok[0].strip()] = float(tok[1].strip())
            self.run_stats = run_stats
            self.magnetization = tuple(mag)
            self.charge = tuple(charge)

    def read_igpar(self):
        """ 
        Renders accessible:
            er_ev = e<r>_ev (dictionary with Spin.up/Spin.down as keys)
            er_bp = e<r>_bp (dictionary with Spin.up/Spin.down as keys)
            er_ev_tot = spin up + spin down summed
            er_bp_tot = spin up + spin down summed
            p_elc = spin up + spin down summed
            p_ion = spin up + spin down summed
        
        (See VASP section 'LBERRY,  IGPAR,  NPPSTR,  DIPOL tags' for info on what these are)"""

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

            def born_ion(results, match): results.born_ion = int(match.group(1))-1; results.born[results.born_ion] = np.zeros((3, 3));
            search.append(['ion +([0-9]+)', lambda results, line: results.born_ion != None, born_ion])

            def born_data(results, match): results.born[results.born_ion][int(match.group(1))-1, :] = np.array([float(match.group(2)), float(match.group(3)), float(match.group(4))]);
            search.append(['^ *([1-3]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)$', lambda results, line: results.born_ion >= 0, born_data])

            def born_section_stop(results, match): results.born_index = None
            search.append(['-------------------------------------', lambda results, line: results.born_ion >= 1, born_section_stop])

            self.born_ion = None
            self.born = {}

            #def debug_print(results,match): print "MATCH:",match.group(0),':',results.born_ion
            #micro_pyawk(filename,search,self,debug=debug_print, postdebug=debug_print)

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

class VolumetricData(object):
    """
    Simple volumetric object for reading LOCPOT and CHGCAR type files.
    
    Attributes:
        name : The name from the comment line.
        poscar : Poscar object
        spinpolarized : True if run is spin polarized
        dim: Tuple of dimensions of volumetric grid in each direction (nx, ny, nz)
        data : Actual data as a dict of {grid coordinate: value}.  Grid coordinate is a (x,y,z) tuple.
        ngridpts: Total number of grid points in volumetric data.
    """
    def __init__(self, filename):
        self.name   = str()
        self.poscar = None
        self.spinpolarized = False
        self.dim = None
        self.data = dict()
        self.ngridpts = 0
        self._read_file(filename)            

    def __add__(self, other):
        return self.linear_add(other, 1.0)

    def __sub__(self, other):
        return self.linear_add(other, -1.0)

    def linear_add(self, other, scalefactor=1.0):
        '''
        Method to do a linear sum of volumetric objects.  Use by + and - operators as well.
        '''
        #To add checks
        summed = VolumetricData()
        summed.name = self.name
        summed.poscar = self.poscar
        summed.spinpolarized = self.spinpolarized
        summed.dim = self.dim
        summed.numpts = self.ngridpts
        for spin in self.data.keys():
            summed.data[spin] = self.data[spin] + scalefactor* other.data[spin]
        return summed

    def _read_file(self, filename):

        reader = file_open_zip_aware(filename)
        lines = reader.readlines()
        reader.close()

        self.poscar = Poscar.from_file(filename)
        
        # Skip whitespace between POSCAR and LOCPOT data
        i = 0
        while lines[i].strip() != "":
            i += 1
        while lines[i].strip() == "":
            i += 1

        dimensionline = lines[i].strip()
        i += 1
        spinpolarized = False
        # Search for the second dimension line, where the next spin starts
        for j in xrange(i, len(lines)):
            if(dimensionline == lines[j].strip()):
                spinpolarized = True
                break

        if not spinpolarized:
            j = j + 2

        self.spinpolarized = spinpolarized

        # Read three numbers that is the dimension
        dimensionexpr = re.compile('([0-9]+) +([0-9]+) +([0-9]+)')
        m = dimensionexpr.match(dimensionline.strip())
        a = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
        data = (" ".join(lines[i:j-1])).split()
        data = data[:(a[0] * a[1] * a[2])]

        self.dim = a
        self.ngridpts = self.dim[0] * self.dim[1] * self.dim[2]
        uppot = np.zeros((a[0], a[1], a[2]))
        count = 0
        for z in xrange(a[2]):
            for y in xrange(a[1]):
                for x in xrange(a[0]):
                    uppot[x, y, z] = float(data[count])
                    count += 1
        if spinpolarized:
            data = (" ".join(lines[j + 1:])).split()
            data = data[:(a[0] * a[1] * a[2])]
            downpot = np.zeros((a[0], a[1], a[2]))
            count = 0
            for (z,y,x) in itertools.product(xrange(a[2]), xrange(a[1]), xrange(a[0])):
                downpot[x, y, z] = float(data[count]) 
                count += 1
            self.data = {Spin.up:uppot, Spin.down:downpot}
        else:
            self.data = {Spin.up:uppot}


class Locpot(VolumetricData):
    """
    Simple object for reading a LOCPOT file
    """
    def __init__(self, filename):
        super(Locpot,self).__init__(filename)

    def get_avg_potential_along_axis(self, ind):
        """
        Get the averaged LOCPOT along a certain axis direction. Useful for visualizing Hartree Potentials.
        
        Args:
            ind : Index of axis.
            
        Returns:
            Average Hatree potential along axis
        """
        m = self.data[Spin.up]
        
        ng = self.dim
        avg = np.zeros((ng[ind], 1))
        for i in xrange(ng[ind]):
            mysum = 0
            for j in xrange(ng[(ind + 1) % 3]):
                for k in xrange(ng[(ind + 2) % 3]):
                    if ind == 0:
                        mysum += m[i,j,k]
                    if ind == 1:
                        mysum += m[k,i,j]
                    if ind == 2:
                        mysum += m[j,k,i]

            avg[i] = mysum / (ng[(ind + 1) % 3] * 1.0) / (ng[(ind + 2) % 3] * 1.0)
        return avg

class Chgcar(VolumetricData):
    """
    Simple object for reading a CHGCAR file
    """
    def __init__(self, filename):
        super(Chgcar,self).__init__(filename)
        #Chgcar format is total density in first set, and moment density in second set.
        # need to split them into up and down.
        updowndata = dict()
        updowndata[Spin.up] = 0.5 * (self.data[Spin.up] + self.data[Spin.down])
        updowndata[Spin.down] = 0.5 * (self.data[Spin.up] - self.data[Spin.down])
        self.data = updowndata
        self._distance_matrix = dict()

    def _calculate_distance_matrix(self, ind):
        structure = self.poscar.struct
        a = self.dim
        distances = dict()
        for (x,y,z) in itertools.product(xrange(a[0]), xrange(a[1]), xrange(a[2])):
            pt = np.array([x / a[0], y / a[1] ,z / a[2]])
            distances[(x,y,z)] = structure[ind].distance_and_image_from_frac_coords(pt)[0]
        self._distance_matrix[ind] = distances

    def get_diff_int_charge(self, ind, radius):
        """
        Get differential integrated charge of atom index ind up to radius.
        
        Args:
            ind : Index of atom.
            radius : Radius of integration.
            
        Returns:
            Differential integrated charge.
        """
        if ind not in self._distance_matrix:
            self._calculate_distance_matrix(ind)
        a = self.dim
        intchg = 0
        for (x,y,z) in itertools.product(xrange(a[0]), xrange(a[1]), xrange(a[2])):
            if self._distance_matrix[ind][(x,y,z)] < radius:
                intchg += self.data[Spin.up][x, y, z] - self.data[Spin.down][x, y, z]
        return intchg / self.ngridpts
    
    def get_diff_int_charge_slow(self, ind, radius):
        """
        Deprecated.  **Much** slower algorithm for finding differential integrated charge.  Used mainly for testing purposes.
        
        Args:
            ind : Index of atom.
            radius : Radius of integration.
            
        Returns:
            Differential integrated charge.       
        """
        st = self.poscar.struct
        a = self.dim
        intchg = 0
        ioncoord = st[ind].frac_coords
        iongridpt = [int(round(ioncoord[i]*a[i])) for i in xrange(3)]
        max_grid_pts = [min(int(round(radius/st.lattice.abc[i] * a[i]))+1,int(round(a[i]/2))) for i in xrange(3)]
        
        for x in xrange(iongridpt[0]-max_grid_pts[0],iongridpt[0]+max_grid_pts[0]):
            for y in xrange(iongridpt[1]-max_grid_pts[1],iongridpt[1]+max_grid_pts[1]):
                for z in xrange(iongridpt[2]-max_grid_pts[2],iongridpt[2]+max_grid_pts[2]):
                    modx = x % a[0]
                    mody = y % a[1]
                    modz = z % a[2]
                    pt = np.array([modx *1.0 / a[0], 1.0 * mody / a[1] , 1.0 * modz / a[2]])
                    dist = st[ind].distance_and_image_from_frac_coords(pt)[0]
                    if dist < radius:
                        intchg += self.data[Spin.up][modx, mody, modz] - self.data[Spin.down][modx, mody, modz]
        return intchg / self.ngridpts

class Procar(object):

    """
    Object for reading a PROCAR file
    """
    def __init__(self, filename):
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
    A basic parser for an OSZICAR output from VASP.  In general, while the OSZICAR is useful for a quick look
    at the output from a VASP run, we recommend that you use the Vasprun parser instead, which gives far richer information
    about a run.
    
    Attributes:
        electronic_steps:
            All electronic steps as a list of list of dict. e.g., 
            [[{'rms': 160.0, 'E': 4507.24605593, 'dE': 4507.2, 'N': 1, 'deps': -17777.0, 'ncg': 16576}, ...], [....]
            where electronic_steps[index] refers the list of electronic steps in one ionic_step, electronic_steps[index][subindex]
            refers to a particular electronic step at subindex in ionic step at index.  The dict of properties depends on the type
            of VASP run, but in general, "E", "dE" and "rms" should be present in almost all runs.
        ionic_steps:
            All ionic_steps as a list of dict, e.g.,
            [{'dE': -526.36, 'E0': -526.36024, 'mag': 0.0, 'F': -526.36024}, ...]
            This is the typical output from VASP at the end of each ionic step.

    Please refer to the vasp manual for the definition for each of the terms.
            
    In addition, two convenience properties, all_energies and final_energy are provided for quick access to the commonly used 
    energetic output from a run.  Please refer to the doc for those two methods for details.    
    """
    
    def __init__(self, filename):
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
                    ionic_steps.append({'F':float(m.group(2)),'E0':float(m.group(3)),'dE':float(m.group(4)),'mag':float(m.group(5))})
                elif re.match("^\s*N\s+E\s*", line):
                    header = re.split("\s+", line.strip().replace("d eps", "deps"))
        self.electronic_steps = electronic_steps
        self.ionic_steps = ionic_steps

    @property
    def all_energies(self):
        """
        Compilation of all energies from all electronic steps and ionic steps as a list of list of energies, e.g.,
        [[4507.24605593, 143.824705755, -512.073149912, -547.713139455, ...], ...]
        """
        all_energies = []
        for i in xrange(len(self.electronic_steps)):
            energies = [step['E'] for step in self.electronic_steps[i]]
            energies.append(self.ionic_steps[i]['F'])
            all_energies.append(energies)
        return all_energies
    
    @property
    def final_energy(self):
        """
        Final energy from run.
        """
        return self.ionic_steps[-1]['F']
    
class VaspParserError(Exception):
    '''
    Exception class for Structure.
    Raised when the structure has problems, e.g., atoms that are too close.
    '''

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return "VaspParserError : " + self.msg
    
def get_band_structure_from_vasp(path):
    """
    method taking a directory with a band structure vasp run
    and returning the corresponding Bandstructure Object
    also takes into account runs that have been separated in several branches
    """
    if(os.path.exists(path+"/branch_0")):
        #get all branches in a list of BandStructurs
        list_branches=[]
        listdir_clean=[]
        for f in os.listdir(path):
            if "branch" in f:
                listdir_clean.append(f)
        for i in range(len(listdir_clean)):
            for f in listdir_clean:
                if(int(f.split("_")[1])==i):
                    list_branches.append(get_band_structure_from_vasp_individual(path+"/"+f))
        return pymatgen.core.electronic_structure.get_reconstructed_band_structure(list_branches)
    else:
        return get_band_structure_from_vasp_individual(path)
    
def get_band_structure_from_vasp_individual(path):
    run=Vasprun(path+"/vasprun.xml")
    labels_dict=parse_kpoint_labels(path+"/KPOINTS")
    lattice_rec=run.final_structure.lattice.reciprocal_lattice
    for c in labels_dict:
        labels_dict[c]=lattice_rec.get_cartesian_coords(labels_dict[c])
    kpoints=[lattice_rec.get_cartesian_coords(np.array(run.actual_kpoints[i])) for i in range(len(run.actual_kpoints))]
    dict_eigen=run.to_dict['output']['eigenvalues']
    eigenvals=[]
    max_band=int(math.floor(len(dict_eigen['1']['up'])*0.9))
    for i in range(max_band):
        eigenvals.append({'energy':[dict_eigen[str(j+1)]['up'][i][0] for j in range(len(kpoints))]})
        eigenvals[i]['occup']=[dict_eigen[str(j+1)]['up'][i][1] for j in range(len(kpoints))]
    bands=Bandstructure(kpoints,eigenvals,labels_dict, run.final_structure, run.efermi)
    return bands

def parse_kpoint_labels(file_kpoints):
    """
    helper method to parse a kpoint file and get the kpoint labels
    """
    with file_open_zip_aware(file_kpoints, "r") as f:
        lines = f.readlines()
    count=-1
    dict_label_kpoints={}
    for line in lines:
        count=count+1
        if(count<4):
            continue
        
        tokens=re.split(" +",line)
        if len(tokens)<6:
            continue
        array=np.array([float(tokens[1]),float(tokens[2]),float(tokens[3])])
        dict_label_kpoints[tokens[5]]=array
        
    return dict_label_kpoints
