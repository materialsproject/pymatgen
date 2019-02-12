# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import re

import numpy as np
from monty.io import zopen

# from monty.re import regrep
from collections import defaultdict

from pymatgen.core.periodic_table import Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.io_utils import clean_lines

"""
This module implements input and output processing from PWSCF.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "3/27/15"


class PWInput:
    """
    Base input file class. Right now, only supports no symmetry and is
    very basic.
    """

    def __init__(self, structure, pseudo=None, control=None, system=None,
                 electrons=None, ions=None, cell=None, kpoints_mode="automatic",
                 kpoints_grid=(1, 1, 1),kpoints_shift=(0, 0, 0)):
        """
        Initializes a PWSCF input file.

        Args:
            structure (Structure): Input structure. For spin-polarized calculation,
                properties (e.g. {"starting_magnetization": -0.5, 
                "pseudo": "Mn.pbe-sp-van.UPF"}) on each site is needed instead of 
                pseudo (dict).
            pseudo (dict): A dict of the pseudopotentials to use. Default to None.
            control (dict): Control parameters. Refer to official PWSCF doc
                on supported parameters. Default to {"calculation": "scf"}
            system (dict): System parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            electrons (dict): Electron parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            ions (dict): Ions parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            cell (dict): Cell parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            kpoints_mode (str): Kpoints generation mode. Default to automatic.
            kpoints_grid (sequence): The kpoint grid. Default to (1, 1, 1).
            kpoints_shift (sequence): The shift for the kpoints. Defaults to
                (0, 0, 0).
        """
        self.structure = structure
        sections = {}
        sections["control"] = control or {"calculation": "scf"}
        sections["system"] = system or {}
        sections["electrons"] = electrons or {}
        sections["ions"] = ions or {}
        sections["cell"] = cell or {}
        if pseudo == None:
            for site in structure:
                try:
                    site.properties['pseudo']
                except KeyError:
                    raise PWInputError("Missing %s in pseudo specification!" 
                                       % site)
        else:
            for species in self.structure.composition.keys():
                if species.symbol not in pseudo:
                    raise PWInputError("Missing %s in pseudo specification!" 
                                       % species.symbol)
        self.pseudo = pseudo

        self.sections = sections
        self.kpoints_mode = kpoints_mode
        self.kpoints_grid = kpoints_grid
        self.kpoints_shift = kpoints_shift

    def __str__(self):
        out = []
        site_descriptions = {}

        if self.pseudo != None:
            site_descriptions = self.pseudo
        else:
            c = 1
            for site in self.structure:
                name = None
                for k, v in site_descriptions.items():
                    if site.properties == v:
                        name = k

                if name == None:
                    name = site.specie.symbol+str(c)
                    site_descriptions[name] = site.properties
                    c += 1

        def to_str(v):
            if isinstance(v, str):
                return "'%s'" % v
            elif isinstance(v, float):
                return "%s" % str(v).replace("e", "d")
            elif isinstance(v, bool):
                if v:
                    return ".TRUE."
                else:
                    return ".FALSE."
            return v

        for k1 in ["control", "system", "electrons", "ions", "cell"]:
            v1 = self.sections[k1]
            out.append("&%s" % k1.upper())
            sub = []
            for k2 in sorted(v1.keys()):
                if isinstance(v1[k2], list):
                    n = 1
                    for l in v1[k2][:len(site_descriptions)]:
                        sub.append("  %s(%d) = %s" % (k2, n, to_str(v1[k2][n-1])))
                        n += 1
                else:
                    sub.append("  %s = %s" % (k2, to_str(v1[k2])))
            if k1 == "system":
                if 'ibrav' not in self.sections[k1]:
                    sub.append("  ibrav = 0")
                if 'nat' not in self.sections[k1]:
                    sub.append("  nat = %d" % len(self.structure))
                if 'ntyp' not in self.sections[k1]:
                    sub.append("  ntyp = %d" % len(site_descriptions))
            sub.append("/")
            out.append(",\n".join(sub))

        out.append("ATOMIC_SPECIES")
        for k, v in sorted(site_descriptions.items(), key=lambda i: i[0]):
            e = re.match(r"[A-Z][a-z]?", k).group(0)
            if self.pseudo is not None:
                p = v
            else:
                p = v['pseudo']
            out.append("  %s  %.4f %s" % (k, Element(e).atomic_mass, p)) 

        out.append("ATOMIC_POSITIONS crystal")
        if self.pseudo is not None:
            for site in self.structure:
                out.append("  %s %.6f %.6f %.6f" % (site.specie.symbol, site.a,
                                                    site.b, site.c))
        else:
            for site in self.structure:
                name = None
                for k, v in sorted(site_descriptions.items(),
                                   key=lambda i: i[0]):
                    if v == site.properties:
                        name = k
                out.append("  %s %.6f %.6f %.6f" % (name, site.a, site.b, site.c))

        out.append("K_POINTS %s" % self.kpoints_mode)
        kpt_str = ["%s" % i for i in self.kpoints_grid]
        kpt_str.extend(["%s" % i for i in self.kpoints_shift])
        out.append("  %s" % " ".join(kpt_str))
        out.append("CELL_PARAMETERS angstrom")
        for vec in self.structure.lattice.matrix:
            out.append("  %f %f %f" % (vec[0], vec[1], vec[2]))
        return "\n".join(out)
    
    def as_dict(self):
        """
        Create a dictionary representation of a PWInput object
        
        Returns:
            dict
        """
        pwinput_dict = {'structure': self.structure.as_dict(),
                        'pseudo': self.pseudo,
                        'sections': self.sections,
                        'kpoints_mode': self.kpoints_mode,
                        'kpoints_grid': self.kpoints_grid,
                        'kpoints_shift': self.kpoints_shift}
        return pwinput_dict
    
    @classmethod
    def from_dict(cls, pwinput_dict):
        """
        Load a PWInput object from a dictionary.
        
        Args:
            pwinput_dict (dict): dictionary with PWInput data
            
        Returns:
            PWInput object
        """
        pwinput = cls(structure=Structure.from_dict(pwinput_dict['structure']),
                          pseudo=pwinput_dict['pseudo'],
                          control=pwinput_dict['sections']['control'],
                          system=pwinput_dict['sections']['system'],
                          electrons=pwinput_dict['sections']['electrons'],
                          ions=pwinput_dict['sections']['ions'],
                          cell=pwinput_dict['sections']['cell'],
                          kpoints_mode=pwinput_dict['kpoints_mode'],
                          kpoints_grid=pwinput_dict['kpoints_grid'],
                          kpoints_shift=pwinput_dict['kpoints_shift'])
        return pwinput

    def write_file(self, filename):
        """
        Write the PWSCF input file.

        Args:
            filename (str): The string filename to output to.
        """
        with open(filename, "w") as f:
            f.write(self.__str__())

    @staticmethod
    def from_file(filename):
        """
        Reads an PWInput object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            PWInput object
        """
        with zopen(filename, "rt") as f:
            return PWInput.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads an PWInput object from a string.

        Args:
            string (str): PWInput string

        Returns:
            PWInput object
        """
        lines = list(clean_lines(string.splitlines()))

        def input_mode(line):
            if line[0] == "&":
                return ("sections", line[1:].lower())
            elif "ATOMIC_SPECIES" in line:
                return ("pseudo", )
            elif "K_POINTS" in line:
                return ("kpoints", line.split("{")[1][:-1])
            elif "CELL_PARAMETERS" in line or "ATOMIC_POSITIONS" in line:
                return ("structure", line.split("{")[1][:-1])
            elif line == "/":
                return None
            else:
                return mode

        sections = {"control": {}, "system": {}, "electrons": {}, 
                    "ions": {}, "cell":{}}
        pseudo = {}
        pseudo_index = 0
        lattice = []
        species = []
        coords = []
        structure = None
        site_properties = {"pseudo":[]}
        mode = None
        for line in lines:
            mode = input_mode(line)
            if mode == None:
                pass
            elif mode[0] == "sections":
                section = mode[1]
                m = re.match(r'(\w+)\(?(\d*?)\)?\s*=\s*(.*)', line)
                if m:
                    key = m.group(1).strip()
                    key_ = m.group(2).strip()
                    val = m.group(3).strip()
                    if key_ != "":
                        if sections[section].get(key, None) == None:
                            val_ = [0.0]*20 # MAX NTYP DEFINITION
                            val_[int(key_)-1] = PWInput.proc_val(key, val)
                            sections[section][key] = val_

                            site_properties[key] = []
                        else:
                            sections[section][key][int(key_)-1] = PWInput.proc_val(key, val) 
                    else:
                        sections[section][key] = PWInput.proc_val(key, val)

            elif mode[0] == "pseudo":
                m = re.match(r'(\w+)\s+(\d*.\d*)\s+(.*)', line)
                if m:
                    pseudo[m.group(1).strip()] = {}
                    pseudo[m.group(1).strip()]["index"] = pseudo_index
                    pseudo[m.group(1).strip()]["pseudopot"] = m.group(3).strip()
                    pseudo_index += 1
            elif mode[0] == "kpoints":
                m = re.match(r'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
                if m:
                    kpoints_grid = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
                    kpoints_shift = (int(m.group(4)), int(m.group(5)), int(m.group(6)))
                else:
                    kpoints_mode = mode[1]
            elif mode[0] == "structure":
                m_l = re.match(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', line)
                m_p = re.match(r'(\w+)\s+(-?\d+\.\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', line)
                if m_l:
                    lattice += [ float(m_l.group(1)), float(m_l.group(2)), float(m_l.group(3)) ]
                elif m_p:
                    site_properties["pseudo"].append(pseudo[m_p.group(1)]["pseudopot"])
                    species += [pseudo[m_p.group(1)]["pseudopot"].split(".")[0]]
                    coords += [[float(m_p.group(2)), float(m_p.group(3)), float(m_p.group(4))]]

                    for k, v in site_properties.items():
                        if k != "pseudo":
                            site_properties[k].append(sections['system'][k][pseudo[m_p.group(1)]["index"]])
                if mode[1] == "angstrom":
                    coords_are_cartesian = True
                elif mode[1] == "crystal":
                    coords_are_cartesian = False

        structure = Structure(Lattice(lattice), species, coords, 
                              coords_are_cartesian=coords_are_cartesian,
                              site_properties=site_properties)
        return PWInput(structure=structure, control=sections["control"],
                       system=sections["system"], electrons=sections["electrons"], 
                       ions=sections["ions"], cell=sections["cell"], kpoints_mode=kpoints_mode,
                       kpoints_grid=kpoints_grid, kpoints_shift=kpoints_shift)

    def proc_val(key, val):
        """
        Static helper method to convert PWINPUT parameters to proper type, e.g.,
        integers, floats, etc.

        Args:
            key: PWINPUT parameter key
            val: Actual value of PWINPUT parameter.
        """
        float_keys = ('etot_conv_thr','forc_conv_thr','conv_thr','Hubbard_U','Hubbard_J0','defauss',
                      'starting_magnetization',)

        int_keys = ('nstep','iprint','nberrycyc','gdir','nppstr','ibrav','nat','ntyp','nbnd','nr1',
                    'nr2','nr3','nr1s','nr2s','nr3s','nspin','nqx1','nqx2','nqx3','lda_plus_u_kind',
                    'edir','report','esm_nfit','space_group','origin_choice','electron_maxstep',
                    'mixing_ndim','mixing_fixed_ns','ortho_para','diago_cg_maxiter','diago_david_ndim',
                    'nraise','bfgs_ndim','if_pos','nks','nk1','nk2','nk3','sk1','sk2','sk3','nconstr')

        bool_keys = ('wf_collect','tstress','tprnfor','lkpoint_dir','tefield','dipfield','lelfield',
                     'lorbm','lberry','lfcpopt','monopole','nosym','nosym_evc','noinv','no_t_rev',
                     'force_symmorphic','use_all_frac','one_atom_occupations','starting_spin_angle',
                     'noncolin','x_gamma_extrapolation','lda_plus_u','lspinorb','london',
                     'ts_vdw_isolated','xdm','uniqueb','rhombohedral','realxz','block',
                     'scf_must_converge','adaptive_thr','diago_full_acc','tqr','remove_rigid_rot',
                     'refold_pos')

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)

        try:
            if key in bool_keys:
                if val.lower() == ".true.":
                    return True
                elif val.lower() == ".false.":
                    return False
                else:
                    raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*d?-?\d*", val.lower()).group(0).replace("d", "e"))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))

        except ValueError:
            pass

        try:
            val = val.replace("d","e")
            return smart_int_or_float(val)
        except ValueError:
            pass

        if "true" in val.lower():
            return True
        if "false" in val.lower():
            return False

        m = re.match(r"^[\"|'](.+)[\"|']$", val)
        if m:
            return m.group(1)



class PWInputError(BaseException):
    pass



def _stress_post(match):
    lines = match.strip().split('\n')
    matrix = np.array([[float(stress)
                        for stress in line.strip().split()]
                       for line in lines])
    stress_kbar = matrix[:,[3,4,5]].tolist()
    return stress_kbar

def _atomic_positions_post(match):
    lines = match.strip().split('\n')
    positions = []
    for line in lines:
        line = line.strip()
        if len(line.split()) == 4:
            specie, x, y, z = line.strip().split()
            position = (specie, [float(pos)
                                 for pos in [x, y, z]])
            positions.append(position)
    return positions

def _cell_parameters_post(match):
    lines = match.strip().split('\n')
    matrix = [[float(param)
               for param in line.strip().split()]
              for line in lines]
    return matrix

def _force_post(match):
    atom = int(match[0])
    type_ = int(match[1])
    force = [float(f) for f in match[2].strip().split()]
    return {'atom': atom, 'type': type_, 'force': force}

def _bands_post(match):
    k_re = re.compile('k\s+=\s+([\s\d\.\-]+)')
    bands_re = re.compile('bands\s+\(ev\)\:\n+([\s\d\-\.]+)', re.MULTILINE)
    occupations_re = re.compile('occupation numbers\s+\n([\s\d\.]+)', re.MULTILINE)

    k = k_re.findall(match)
    k = [[float(j) for j in i.replace('-', ' -').strip().split()] for i in k]
    bands = bands_re.findall(match)
    bands = [[float(j) for j in i.strip().split()] for i in bands]
    occupations = occupations_re.findall(match)
    occupations = [[float(j) for j in i.strip().split()] for i in occupations]

    return {'kpoints': k, 'bands': bands, 'occupations': occupations}

def _kpoints_post(match):
    lines = match.strip().split('\n')

    kpoints = []
    for line in lines:
        line = line.strip().split()
        index = int(line[1].strip(')'))
        x, y, z = [float(i.strip('),')) for i in line[4:7]]
        weight = float(line[-1])
        kpoints.append({'index': index,
                        'coords': [x, y, z],
                        'weight': weight})
    return kpoints

def _initial_atomic_positions_post(match):
    lines = match.strip().split('\n')
    positions = []
    for line in lines:
        line = line.split()
        site = int(line[0])
        specie = line[1]
        coords = [float(i) for i in [6, 7, 8]]
        positions.append((specie, coords))

    return positions



class PWOutput:
    
    patterns = {
        'energy': {
            'pattern': r'total energy\s+=\s+([\d\.\-]+)\s+Ry',
            'flags': [],
            'postprocess': float,
        },
        'final_energy': {
            'pattern': r'!\s+total energy\s+=\s+([\d\.\-]+)\s+Ry',
            'flags': [],
            'postprocess': float,
        },
        'enthalpy': {
            'pattern': r'enthalpy new\s+=\s+([\d\.\-]+)\s+Ry',
            'flags': [],
            'postprocess': float,
        },
        'final_enthalpy': {
            'pattern': r'Final enthalpy\s+=\s+([\d\.\-]+)\s+Ry',
            'flags': [],
            'postprocess': float,
        },
        'density': {
            'pattern': r'density\s+=\s+([\d\.]+)\s+g\/cm\^3',
            'flags': [],
            'postprocess': float,
        },
        'warning': {
            'pattern': r'Warning:\s+([\w\s\/\&]+)$',
            'flags': [re.MULTILINE],
            'postprocess': lambda x: str(x).strip(),
        },
        'total_stress': {
            'pattern': r'total\s+stress\s+\(Ry\/bohr\*\*3\)\s+\(kbar\)\s+P=\s+([\d\.\-]+)',
            'flags': [],
            'postprocess': float,
        },
        'stress': {
            'pattern': r'total\s+stress\s+\(Ry\/bohr\*\*3\)\s+\(kbar\)\s+P=\s+[\d\.\-]+\n([\s\d\.\-]+)\n',
            'flags': [re.MULTILINE],
            'postprocess': _stress_post,
        },
        'total_force': {
            'pattern': r'Total force\s+=\s+([\d\.\-]+)',
            'flags': [],
            'postprocess': float,
        },
        'force': {
            'pattern': r'Forces acting on atoms \(cartesian axes, Ry\/au\)\:\n\n\s+atom\s+([\d]+)\s+type\s+([\d]+)\s+force\s+=\s+([\s\d\.\-]+)',
            'flags': [re.MULTILINE],
            'postprocess': _force_post,
        },
        'bands_data': {
            'pattern': r'End of self\-consistent calculation\n\n(.*?)the Fermi energy is',
            'flags': [re.DOTALL],
            'postprocess': _bands_post,
        },
        'fermi_energy': {
            'pattern': r'the Fermi energy is\s+([\d\.]+) ev',
            'flags': [],
            'postprocess': float,
        },
        'conv_iters': {
            'pattern': r'convergence has been achieved in\s+([\d+]) iterations',
            'flags': [],
            'postprocess': int,
        },
        'cell_parameters': {
            'pattern': r'CELL\_PARAMETERS\s+\(angstrom\)\s+([\s\d\.\-]+)^$',
            'flags': [re.MULTILINE],
            'postprocess': _cell_parameters_post,
        },
        'atomic_positions': {
            'pattern': r'ATOMIC_POSITIONS\s+\(crystal\)\s+([\w\s\d\.\-\n]+)^$',
            'flags': [re.MULTILINE],
            'postprocess': _atomic_positions_post,
        },
        'version': {
            'pattern': r'Program PWSCF v.([\d\.]+)',
            'flags': [],
            'postprocess': float,
        },
        'date': {
            'pattern': r'Program PWSCF v.[\d\.]+ starts on\s+([\S]+)',
            'flags': [],
            'postprocess': lambda x: str(x).strip(),
        },
        'time': {
            'pattern': r'Program PWSCF v.[\d\.]+ starts on\s+\S+\s+at\s+([\d\:\s]+)$',
            'flags': [re.MULTILINE],
            'postprocess': lambda x: str(x).strip(),
        },
        'lattice_type': {
            'pattern': r'bravais\-lattice index\s+=\s+(\d+)',
            'flags': [],
            'postprocess': int,
        },
        'lattice_parameter': {
            'pattern': r'lattice parameter \(alat\)\s+=\s+([\d\.]+)',
            'flags': [],
            'postprocess': float,
        },
        'unit_cell_volume': {
            'pattern': r'unit-cell volume\s+=\s+([\d\.]+)',
            'flags': [],
            'postprocess': float,
        },
        'nat': {
            'pattern': r'number of atoms\/cell\s+=\s+(\d+)',
            'flags': [],
            'postprocess': int,
        },
        'ntype': {
            'pattern': r'number of atomic types\s+=\s+(\d+)',
            'flags': [],
            'postprocess': int,
        },
        'nelectrons': {
            'pattern': r'number of electrons\s+=\s+([\d\.]+)',
            'flags': [],
            'postprocess': float,
        },
        'nks_states': {
            'pattern': r'number of Kohn\-Sham states=\s+(\d+)',
            'flags': [],
            'postprocess': int,
        },
        'ecutwfc': {
            'pattern': r'kinetic\-energy cutoff\s+=\s+([\d\.\-]+)\s+Ry',
            'flags': [],
            'postprocess': float,
        },
        'echutrho': {
            'pattern': r'charge density cutoff\s+=\s+([\d\.\-]+)\s+Ry',
            'flags': [],
            'postprocess': float,
        },
        'conv_thr': {
            'pattern': r'convergence threshold\s+=\s+([\d\.\-E]+)',
            'flags': [],
            'postprocess': float,
        },
        'mixing_beta': {
            'pattern': r'mixing beta\s+=\s+([\d\.]+)',
            'flags': [],
            'postprocess': float,
        },
        'niter': {
            'pattern': r'number of iterations used\s+=\s+([\d\w\s]+)$',
            'flags': [re.MULTILINE],
            'postprocess': lambda x: str(x).strip(),
        },
        'exc': {
            'pattern': r'Exchange\-correlation\s+=\s+([\d\w\s\(\)]+)$',
            'flags': [re.MULTILINE],
            'postprocess': lambda x: str(x).strip(),
        },
        'celldm1': {
            'pattern': r'celldm\(1\)=\s+([\d\.]+)\s',
            'flags': [],
            'postprocess': float,
        },
        'celldm2': {
            'pattern': r'celldm\(2\)=\s+([\d\.]+)\s',
            'flags': [],
            'postprocess': float,
        },
        'celldm3': {
            'pattern': r'celldm\(3\)=\s+([\d\.]+)\s',
            'flags': [],
            'postprocess': float,
        },
        'celldm4': {
            'pattern': r'celldm\(4\)=\s+([\d\.]+)\s',
            'flags': [],
            'postprocess': float,
        },
        'celldm5': {
            'pattern': r'celldm\(5\)=\s+([\d\.]+)\s',
            'flags': [],
            'postprocess': float,
        },
        'celldm6': {
            'pattern': r'celldm\(6\)=\s+([\d\.]+)\s',
            'flags': [],
            'postprocess': float,
        },
        'a1': {
            'pattern': r'a\(1\)\s+=\s+\(\s+([\d\s\.\-]+)\s+\)',
            'flags': [],
            'postprocess': lambda x: [float(i) for i in str(x).split()],
        },
        'a2': {
            'pattern': r'a\(2\)\s+=\s+\(\s+([\d\s\.\-]+)\s+\)',
            'flags': [],
            'postprocess': lambda x: [float(i) for i in str(x).split()],
        },
        'a3': {
            'pattern': r'a\(3\)\s+=\s+\(\s+([\d\s\.\-]+)\s+\)',
            'flags': [],
            'postprocess': lambda x: [float(i) for i in str(x).split()],
        },
        'b1': {
            'pattern': r'b\(1\)\s+=\s+\(\s+([\d\s\.\-]+)\s+\)',
            'flags': [],
            'postprocess': lambda x: [float(i) for i in str(x).split()],
        },
        'b2': {
            'pattern': r'b\(2\)\s+=\s+\(\s+([\d\s\.\-]+)\s+\)',
            'flags': [],
            'postprocess': lambda x: [float(i) for i in str(x).split()],
        },
        'b3': {
            'pattern': r'b\(3\)\s+=\s+\(\s+([\d\s\.\-]+)\s+\)',
            'flags': [],
            'postprocess': lambda x: [float(i) for i in str(x).split()],
        },
        'nsymop': {
            'pattern': r'^\s+([\d]+)\s+Sym\. Ops\.',
            'flags': [re.MULTILINE],
            'postprocess': int,
        },
        # TODO: symmetry operations (frac)
        # TODO: symmetry operations (cart)
        'initial_atomic_positions_cart': {
            'pattern': r'Cartesian axes[\s\n]+site n\.\s+atom\s+positions \(alat units\)\n(.*?)Crystallographic',
            'flags': [re.DOTALL],
            'postprocess': _initial_atomic_positions_post,
        },
        'initial_atomic_positions_frac': {
            'pattern': r'Crystallographic axes[\s\n]+site n\.\s+atom\s+positions \(cryst\. coord\.\)\n(.*?)number',
            'flags': [re.DOTALL],
            'postprocess': _initial_atomic_positions_post,
        },
        'nkpts': {
            'pattern': r'number of k points=\s+([\d]+)',
            'flags': [],
            'postprocess': int,
        },
        'smearing': {
            'pattern': r'number of k points=\s+[\d]+\s+([\S]+) smearing',
            'flags': [],
            'postprocess': lambda x:x,
        },
        'degauss': {
            'pattern': r'number of k points=\s+\d+\s+\S+ smearing, width\s+\(Ry\)=\s+([\d\.]+)',
            'flags': [],
            'postprocess': float,
        },
        'kpoints_cart': {
            'pattern': r'cart\. coord\. in units 2pi\/alat(.*?)cryst\. coord\.',
            'flags': [re.DOTALL],
            'postprocess': _kpoints_post,
        },
        'kpoints_frac': {
            'pattern': r'k\( .*?cryst\.\s+coord\.(.*?)Dense\s+grid',
            'flags': [re.DOTALL],
            'postprocess': _kpoints_post,
        },
        'job_done': {
            'pattern': r'JOB DONE\.',
            'flags': [],
            'postprocess': lambda x:x,
        },
    }
    
    def __init__(self, filename, data=defaultdict(list)):
        self.filename = filename
        self.data = data
        self.read_pattern(PWOutput.patterns)
        
    def read_pattern(self, patterns):
        '''
        General pattern reading uses stdlib's re module.
        
        Args:
            patterns (dict): A dict of patterns, e.g.,
                {"energy": {
                    "pattern": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)",
                    "flags": [],
                    "postprocess": str}
                }
                
        Renders accessible:
            Any attribute in patterns. For example, the energy example above
            will set the value of self.data["energy"] = [-1234, -3453, ...],
            to the results from regex and postprocess
        '''
        with zopen(self.filename, 'rt') as file_:
            out = file_.read()
        all_matches = {}
        for key, value in patterns.items():
            pattern = re.compile(value['pattern'], *value['flags'])
            matches = pattern.findall(out)
            matches = [value['postprocess'](match) for match in matches]
            all_matches[key] = matches
        self.data.update(all_matches)
    
    def get_celldm(self, i):
        return self.data['celldm%d' % i]
    
    # FIXME: add initial structures which are
    #    printed in a different format at the beginning
    #    of the calculation and the final scf step
    @property
    def structures(self):
        cells = self.data['cell_parameters']
        atpos = self.data['atomic_positions']

        structures = []
        for cell, at in zip(cells, atpos):
            species = [a[0] for a in at]
            coords = [a[1] for a in at]
            structure = Structure(cell, species, coords)
            structures.append(structure)

        return structures
    
    # TODO: band structure
    
    # TODO: symmetry operations
    
    @property
    def final_energy(self):
        return self.data['final_energy'][-1]
    
    @property
    def lattice_type(self):
        return self.data['lattice_type']
    
    def as_dict(self):
        dict_ = {
            'filename': self.filename,
            'data': self.data,
        }
        return dict_
    
    @classmethod
    def from_dict(cls, dict_):
        return cls(**dict_)
