#!/usr/bin/env python

"""
This module defines the VaspInputSet abstract base class and a concrete
implementation for the parameters used by the Materials Project and the MIT
high throughput project.  The basic concept behind an input set is to specify
a scheme to generate a consistent set of Vasp inputs from a structure
without further user intervention. This ensures comparability across
runs.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Wei Chen, Will Richards"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 16, 2011"

import os
import abc
import ConfigParser
import json
import re

from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vaspio.vasp_output import Vasprun, Outcar
from pymatgen.serializers.json_coders import MSONable
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.util.decorators import deprecated
import traceback
import numpy as np
import shutil


class AbstractVaspInputSet(MSONable):
    """
    Abstract base class representing a set of Vasp input parameters.
    The idea is that using a VaspInputSet, a complete set of input files
    (INPUT, KPOINTS, POSCAR and POTCAR) can be generated in an automated
    fashion for any structure.
    """
    __metaclass__ = abc.ABCMeta

    def get_poscar(self, structure):
        """
        Returns Poscar from a structure.
        """
        return Poscar(structure)

    @abc.abstractmethod
    def get_kpoints(self, structure):
        """
        Returns Kpoints from a structure.

        Args:
            structure:
                Structure object

        Returns:
            Kpoints object
        """
        return

    @abc.abstractmethod
    def get_incar(self, structure):
        """
        Returns Incar from a structure.

        Args:
            structure:
                Structure object

        Returns:
            Incar object
        """
        return

    @abc.abstractmethod
    def get_potcar(self, structure):
        """
        Returns Potcar from a structure.

        Args:
            structure:
                Structure object

        Returns:
            Potcar object
        """
        return

    @abc.abstractmethod
    def get_potcar_symbols(self, structure):
        """
        Returns list of POTCAR symbols from a structure.

        Args:
            structure:
                Structure object

        Returns:
            List of POTCAR symbols
        """
        return

    def get_all_vasp_input(self, structure, generate_potcar=True):
        """
        Returns all input files as a dict of {filename: vaspio object}

        Args:
            structure:
                Structure object
            generate_potcar:
                Set to False to generate a POTCAR.spec file instead of a
                POTCAR, which contains the POTCAR labels but not the actual
                POTCAR. Defaults to True.

        Returns:
            dict of {filename: file_as_string}, e.g., {'INCAR':'EDIFF=1e-4...'}
        """
        d = {'INCAR': self.get_incar(structure),
             'KPOINTS': self.get_kpoints(structure),
             'POSCAR': self.get_poscar(structure)}
        if generate_potcar:
            d['POTCAR'] = self.get_potcar(structure)
        else:
            d['POTCAR.spec'] = "\n".join(self.get_potcar_symbols(structure))
        return d

    def write_input(self, structure, output_dir, make_dir_if_not_present=True):
        """
        Writes a set of VASP input to a directory.

        Args:
            structure:
                Structure object
            output_dir:
                Directory to output the VASP input files
            make_dir_if_not_present:
                Set to True if you want the directory (and the whole path) to
                be created if it is not present.
        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.get_all_vasp_input(structure).items():
            v.write_file(os.path.join(output_dir, k))


class DictVaspInputSet(AbstractVaspInputSet):
    """
    Concrete implementation of VaspInputSet that is initialized from a dict
    settings. This allows arbitrary settings to be input. In general,
    this is rarely used directly unless there is a source of settings in JSON
    format (e.g., from a REST interface). It is typically used by other
    VaspInputSets for initialization.

    Special consideration should be paid to the way the MAGMOM initialization
    for the INCAR is done. The initialization differs depending on the type of
    structure and the configuration settings. The order in which the magmom is
    determined is as follows:

    1. If the site itself has a magmom setting, that is used.
    2. If the species on the site has a spin setting, that is used.
    3. If the species itself has a particular setting in the config file, that
       is used, e.g., Mn3+ may have a different magmom than Mn4+.
    4. Lastly, the element symbol itself is checked in the config file. If
       there are no settings, VASP's default of 0.6 is used.


    """
    def __init__(self, name, config_dict, constrain_total_magmom=False):
        """
        Args:
            name:
                The name in the config file.
            config_dict:
                The config dictionary to use.
            user_incar_settings:
                User INCAR settings. This allows a user to override INCAR
                settings, e.g., setting a different MAGMOM for various elements
                or species.
            constrain_total_magmom:
                Whether to constrain the total magmom (NUPDOWN in INCAR) to be
                the sum of the expected MAGMOM for all species. Defaults to
                False.
        """
        self.name = name
        self.potcar_settings = config_dict["POTCAR"]
        self.kpoints_settings = config_dict['KPOINTS']
        self.incar_settings = config_dict['INCAR']
        self.set_nupdown = constrain_total_magmom

    def get_incar(self, structure):
        incar = Incar()
        comp = structure.composition
        elements = sorted([el for el in comp.elements if comp[el] > 0],
                          key=lambda el: el.X)
        most_electroneg = elements[-1].symbol
        poscar = Poscar(structure)
        for key, setting in self.incar_settings.items():
            if key == "MAGMOM":
                mag = []
                for site in structure:
                    if hasattr(site, 'magmom'):
                        mag.append(site.magmom)
                    elif hasattr(site.specie, 'spin'):
                        mag.append(site.specie.spin)
                    elif str(site.specie) in setting:
                        mag.append(setting.get(str(site.specie)))
                    else:
                        mag.append(setting.get(site.specie.symbol, 0.6))
                incar[key] = mag
            elif key in ('LDAUU', 'LDAUJ', 'LDAUL'):
                if most_electroneg in setting.keys():
                    incar[key] = [setting[most_electroneg].get(sym, 0)
                                  for sym in poscar.site_symbols]
                else:
                    incar[key] = [0] * len(poscar.site_symbols)
            elif key == "EDIFF":
                incar[key] = float(setting) * structure.num_sites
            else:
                incar[key] = setting

        has_u = ("LDAUU" in incar and sum(incar['LDAUU']) > 0)
        if has_u:
            # modify LMAXMIX if LSDA+U and you have d or f electrons
            # note that if the user explicitly sets LMAXMIX in settings it will
            # override this logic.
            if 'LMAXMIX' not in self.incar_settings.keys():
                # contains f-electrons
                if any([el.Z > 56 for el in structure.composition]):
                    incar['LMAXMIX'] = 6
                # contains d-electrons
                elif any([el.Z > 20 for el in structure.composition]):
                    incar['LMAXMIX'] = 4
        else:
            for key in incar.keys():
                if key.startswith('LDAU'):
                    del incar[key]

        if self.set_nupdown:
            nupdown = sum([mag if mag > 0.6 else 0 for mag in incar['MAGMOM']])
            incar['NUPDOWN'] = nupdown

        return incar

    def get_potcar(self, structure):
        return Potcar(self.get_potcar_symbols(structure))

    def get_potcar_symbols(self, structure):
        p = self.get_poscar(structure)
        elements = p.site_symbols
        potcar_symbols = []
        for el in elements:
            potcar_symbols.append(self.potcar_settings[el]
                                  if el in self.potcar_settings else el)
        return potcar_symbols

    def get_kpoints(self, structure):
        """
        Writes out a KPOINTS file using the fully automated grid method. Uses
        Gamma centered meshes  for hexagonal cells and Monk grids otherwise.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.
        """
        dens = int(self.kpoints_settings['grid_density'])
        return Kpoints.automatic_density(structure, dens)

    def __str__(self):
        return self.name

    def __repr__(self):
        output = [self.name, ""]
        section_names = ['INCAR settings', 'KPOINTS settings',
                         'POTCAR settings']
        count = 0
        for d in [self.incar_settings, self.kpoints_settings,
                  self.potcar_settings]:
            output.append(section_names[count])
            for k, v in d.items():
                output.append("%s = %s" % (k, str(v)))
            output.append("")
            count += 1
        return "\n".join(output)

    @property
    def to_dict(self):
        config_dict = {
            "INCAR": self.incar_settings,
            "KPOINTS": self.kpoints_settings,
            "POTCAR": self.potcar_settings
        }
        return {
            "name": self.name,
            "config_dict": config_dict,
            "constrain_total_magmom": self.set_nupdown,
            "@class": self.__class__.__name__,
            "@module": self.__module__.__name__,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(d["name"], d["config_dict"],
                   d["constrain_total_magmom"])


class VaspInputSet(DictVaspInputSet):
    """
    Standard implementation of VaspInputSet that uses a config file to
    initialize settings. See DictVaspInputSet for specific details regarding
    how MAGMOM, LDAU settings are set.
    """
    def __init__(self, name, config_file, user_incar_settings=None,
                 constrain_total_magmom=False):
        """
        Args:
            name:
                The name in the config file.
            config_file:
                The config file to use. If None (the default), a default config
                file containing Materials Project and MIT parameters is used.
            user_incar_settings:
                User INCAR settings. This allows a user to override INCAR
                settings, e.g., setting a different MAGMOM for various elements
                or species.
            constrain_total_magmom:
                Whether to constrain the total magmom (NUPDOWN in INCAR) to be
                the sum of the expected MAGMOM for all species. Defaults to
                False.
        """
        self.name = name
        self.config_file = config_file
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        with open(config_file, "r") as f:
            self._config.readfp(f)

        self.user_incar_settings = user_incar_settings
        potcar_settings = dict(self._config.items(self.name + 'POTCAR'))
        kpoints_settings = dict(self._config.items(self.name + 'KPOINTS'))
        incar_settings = dict(self._config.items(self.name + 'INCAR'))
        for key in ['MAGMOM', 'LDAUU', 'LDAUJ', 'LDAUL']:
            if key in incar_settings:
                incar_settings[key] = json.loads(incar_settings[key])
        if user_incar_settings:
            incar_settings.update(user_incar_settings)
        DictVaspInputSet.__init__(
            self, name, {"INCAR": incar_settings, "KPOINTS": kpoints_settings,
                         "POTCAR": potcar_settings},
            constrain_total_magmom=constrain_total_magmom)

    @property
    def to_dict(self):
        return {
            "name": self.name,
            "config_file": self.config_file,
            "constrain_total_magmom": self.set_nupdown,
            "user_incar_settings": self.user_incar_settings,
            "@class": self.__class__.__name__,
            "@module": self.__module__.__name__,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(d["name"], d["config_file"],
                   user_incar_settings=d["user_incar_settings"],
                   constrain_total_magmom=d["constrain_total_magmom"])


class MITVaspInputSet(VaspInputSet):
    """
    Standard implementation of VaspInputSet utilizing parameters in the MIT
    High-throughput project.
    The parameters are chosen specifically for a high-throughput project,
    which means in general pseudopotentials with fewer electrons were chosen.

    Please refer to A Jain, G. Hautier, C. Moore, S. P. Ong, C. Fischer,
    T. Mueller, K. A. Persson, G. Ceder. A high-throughput infrastructure for
    density functional theory calculations. Computational Materials Science,
    2011, 50(8), 2295-2310.
    doi:10.1016/j.commatsci.2011.02.023 for more information.
    """
    def __init__(self, user_incar_settings=None, constrain_total_magmom=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(module_dir, "MITVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MIT", json.load(f),
                constrain_total_magmom=constrain_total_magmom)
        self.user_incar_settings = user_incar_settings
        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

    @property
    def to_dict(self):
        return {
            "name": self.name,
            "constrain_total_magmom": self.set_nupdown,
            "user_incar_settings": self.user_incar_settings,
            "@class": self.__class__.__name__,
            "@module": self.__class__.__module__,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(user_incar_settings=d["user_incar_settings"],
                   constrain_total_magmom=d["constrain_total_magmom"])


class MITGGAVaspInputSet(VaspInputSet):
    """
    Typical implementation of input set for a GGA run based on MIT parameters.
    """
    def __init__(self, user_incar_settings=None, constrain_total_magmom=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(module_dir, "MITVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MIT GGA", json.load(f),
                constrain_total_magmom=constrain_total_magmom)
        self.user_incar_settings = user_incar_settings
        # INCAR settings to override for GGA runs, since we are basing off a
        # GGA+U inputset
        self.incar_settings['LDAU'] = False
        if 'LDAUU' in self.incar_settings:
            # technically not needed, but clarifies INCAR
            del self.incar_settings['LDAUU']
            del self.incar_settings['LDAUJ']
            del self.incar_settings['LDAUL']
        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

    @property
    def to_dict(self):
        return {
            "constrain_total_magmom": self.set_nupdown,
            "user_incar_settings": self.user_incar_settings,
            "@class": self.__class__.__name__,
            "@module": self.__class__.__module__,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(user_incar_settings=d["user_incar_settings"],
                   constrain_total_magmom=d["constrain_total_magmom"])


class MITHSEVaspInputSet(VaspInputSet):
    """
    Typical implementation of input set for a HSE run.
    """
    def __init__(self, user_incar_settings=None, constrain_total_magmom=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(module_dir, "MITHSEVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MIT HSE", json.load(f),
                constrain_total_magmom=constrain_total_magmom)
        self.user_incar_settings = user_incar_settings
        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

    @property
    def to_dict(self):
        return {
            "constrain_total_magmom": self.set_nupdown,
            "user_incar_settings": self.user_incar_settings,
            "@class": self.__class__.__name__,
            "@module": self.__class__.__module__,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(user_incar_settings=d["user_incar_settings"],
                   constrain_total_magmom=d["constrain_total_magmom"])


class MITMDVaspInputSet(VaspInputSet):
    """
    Class for writing a vasp md run. This DOES NOT do multiple stage
    runs.
    """

    def __init__(self, start_temp, end_temp, nsteps, time_step=2,
                 prec="Low", ggau=False, user_incar_settings=None):
        """
        Args:
            start_temp:
                Starting temperature.
            end_temp:
                Final temperature.
            nsteps:
                Number of time steps for simulations. The NSW parameter.
            time_step:
                The time step for the simulation. The POTIM parameter.
                Defaults to 2fs.
            prec:
                precision - Normal or LOW. Defaults to Low.
            ggau:
                whether to use +U or not. Defaults to False.
            user_incar_settings:
                dictionary of incar settings to override
        """
        module_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(module_dir, "MITVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MITMD", json.load(f))
        self.start_temp = start_temp
        self.end_temp = end_temp
        self.nsteps = nsteps
        self.time_step = time_step
        self.prec = prec
        self.ggau = ggau
        self.user_incar_settings = user_incar_settings

        #Optimized parameters for MD simulations.
        incar_settings = {'TEBEG': start_temp, 'TEEND': end_temp,
                          'NSW': nsteps, 'PREC': prec,
                          'EDIFF': 0.000001, 'LSCALU': False,
                          'LCHARG': False, 'LPLANE': False,
                          'LWAVE': False, "ICHARG": 1, "ISMEAR": 0,
                          "SIGMA": 0.05, "NELMIN": 4, "LREAL": True,
                          "BMIX": 1, "MAXMIX": 20, "NELM": 500, "NSIM": 4,
                          "ISYM": 0, "ISIF": 0, "IBRION": 0, "NBLOCK": 1,
                          "KBLOCK": 100, "SMASS": 0, "POTIM": time_step}

        self.incar_settings.update(incar_settings)

        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

        if not ggau:
            self.incar_settings['LDAU'] = False
            if 'LDAUU' in self.incar_settings:
                # technically not needed, but clarifies INCAR
                del self.incar_settings['LDAUU']
                del self.incar_settings['LDAUJ']
                del self.incar_settings['LDAUL']

    def get_kpoints(self, structure):
        return Kpoints.gamma_automatic()

    @property
    def to_dict(self):
        return {
            "start_temp": self.start_temp,
            "end_temp": self.end_temp,
            "nsteps": self.nsteps,
            "time_step": self.time_step,
            "ggau": self.ggau,
            "prec": self.prec,
            "user_incar_settings": self.user_incar_settings,
            "@class": self.__class__.__name__,
            "@module": self.__class__.__module__,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(start_temp=d["start_temp"], end_temp=d["end_temp"],
                   nsteps=d["nsteps"], time_step=d["time_step"],
                   prec=d["prec"], ggau=d["ggau"],
                   user_incar_settings=d["user_incar_settings"])


class MPVaspInputSet(DictVaspInputSet):
    """
    Implementation of VaspInputSet utilizing parameters in the public
    Materials Project. Typically, the pseudopotentials chosen contain more
    electrons than the MIT parameters, and the k-point grid is ~50% more dense.
    The LDAUU parameters are also different due to the different psps used,
    which result in different fitted values (even though the methodology of
    fitting is exactly the same as the MIT scheme).
    """
    def __init__(self, user_incar_settings=None, constrain_total_magmom=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(module_dir, "MPVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MaterialsProject", json.load(f),
                constrain_total_magmom=constrain_total_magmom)
        self.user_incar_settings = user_incar_settings
        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

    @property
    def to_dict(self):
        return {
            "constrain_total_magmom": self.set_nupdown,
            "user_incar_settings": self.user_incar_settings,
            "@class": self.__class__.__name__,
            "@module": self.__class__.__module__,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            user_incar_settings=d["user_incar_settings"],
            constrain_total_magmom=d["constrain_total_magmom"])


class MPGGAVaspInputSet(DictVaspInputSet):
    """
    Same as the MaterialsProjectVaspInput set, but the +U is enforced to be
    turned off.
    """
    def __init__(self, user_incar_settings=None, constrain_total_magmom=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(module_dir, "MPVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MaterialsProject GGA", json.load(f),
                constrain_total_magmom=constrain_total_magmom)

        self.user_incar_settings = user_incar_settings
        # INCAR settings to override for GGA runs, since we are basing off a
        # GGA+U inputset
        self.incar_settings['LDAU'] = False
        if 'LDAUU' in self.incar_settings:
            # technically not needed, but clarifies INCAR
            del self.incar_settings['LDAUU']
            del self.incar_settings['LDAUJ']
            del self.incar_settings['LDAUL']
        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

    @property
    def to_dict(self):
        return {
            "constrain_total_magmom": self.set_nupdown,
            "user_incar_settings": self.user_incar_settings,
            "@class": self.__class__.__name__,
            "@module": self.__class__.__module__,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(user_incar_settings=d["user_incar_settings"],
                   constrain_total_magmom=d["constrain_total_magmom"])


class MPStaticVaspInputSet(MPVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for static calculations that typically follow relaxation runs.
    It is recommended to use the static from_previous_run method to construct
    the input set to inherit most of the functions.
    """
    def __init__(self, user_incar_settings=None, constrain_total_magmom=False):
        """
        Args:
            user_incar_settings:
                A dict specify customized settings for INCAR
        """
        module_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(module_dir, "MPVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MaterialsProject Static", json.load(f),
                constrain_total_magmom=constrain_total_magmom)

        self.user_incar_settings = user_incar_settings
        self.incar_settings.update(
            {"IBRION": -1, "ISMEAR": -5, "LAECHG": True, "LCHARG": True,
             "LORBIT": 11, "LVHAR": True, "LWAVE": False, "NSW": 0})
        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

    def get_kpoints(self, structure, kpoints_density=90):
        """
        Get a KPOINTS file using the fully automated grid method. Uses
        Gamma centered meshes for hexagonal cells and Monk grids otherwise.

        Args:
            kpoints_density:
                kpoints density for the reciprocal cell of structure.
                Might need to increase the default value when calculating
                metallic materials.
        """
        kpoints_density = kpoints_density
        self.kpoints_settings['grid_density'] = kpoints_density * \
            structure.lattice.reciprocal_lattice.volume * structure.num_sites
        return super(MPStaticVaspInputSet, self).get_kpoints(structure)

    @staticmethod
    def get_structure(vasp_run, outcar=None, initial_structure=False,
                      additional_info=False):
        """
        Process structure for static calculations from previous run.

        Args:
            vasp_run:
                Vasprun object that contains the final structure from previous
                run.
            outcar:
                Outcar object that contains the magnetization info from
                previous run.
            initial_structure:
                Whether to return the structure from previous run. Default is
                False.
            additional_info:
                Whether to return additional symmetry info related to the
                structure. If True, return a list of the refined structure (
                conventional cell), the conventional standard structure,
                the symmetry dataset and symmetry operations of the structure
                (see SymmetryFinder doc for details)

        Returns:
            Returns the magmom-decorated structure that can be passed to get
            Vasp input files, e.g. get_kpoints.
        """
        #TODO: fix magmom for get_*_structures
        if vasp_run.is_spin:
            if outcar and outcar.magnetization:
                magmom = {"magmom": [i['tot'] for i in outcar.magnetization]}
            else:
                magmom = {
                    "magmom": vasp_run.to_dict['input']['parameters']
                    ['MAGMOM']}
        else:
            magmom = None
        structure = vasp_run.final_structure
        if magmom:
            structure = structure.copy(site_properties=magmom)
        sym_finder = SymmetryFinder(structure, symprec=0.01)
        if initial_structure:
            return structure
        elif additional_info:
            info = [sym_finder.get_refined_structure(),
                    sym_finder.get_conventional_standard_structure(),
                    sym_finder.get_symmetry_dataset(),
                    sym_finder.get_symmetry_operations()]
            return [sym_finder.get_primitive_standard_structure(),
                    info]
        else:
            return sym_finder.get_primitive_standard_structure()

    @staticmethod
    def from_previous_vasp_run(previous_vasp_dir, output_dir='.',
                               user_incar_settings=None,
                               make_dir_if_not_present=True):
        """
        Generate a set of Vasp input files for static calculations from a
        directory of previous Vasp run.

        Args:
            previous_vasp_dir:
                The directory contains the outputs(vasprun.xml and OUTCAR) of
                previous vasp run.
            output_dir:
                The directory to write the VASP input files for the static
                calculations. Default to write in the current directory.
            make_dir_if_not_present:
                Set to True if you want the directory (and the whole path) to
                be created if it is not present.
        """

        try:
            vasp_run = Vasprun(os.path.join(previous_vasp_dir, "vasprun.xml"),
                               parse_dos=False, parse_eigen=None)
            outcar = Outcar(os.path.join(previous_vasp_dir, "OUTCAR"))
        except:
            traceback.format_exc()
            raise RuntimeError("Can't get valid results from previous run")

        structure = MPStaticVaspInputSet.get_structure(
            vasp_run, outcar)
        mpsvip = MPStaticVaspInputSet(
            user_incar_settings=user_incar_settings)
        mpsvip.write_input(structure, output_dir, make_dir_if_not_present)


class MPNonSCFVaspInputSet(MPStaticVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for non self-consistent field (NonSCF) calculation that follows
    a static run to calculate bandstructure, density of states(DOS) and etc.
    It is recommended to use the NonSCF from_previous_run method to construct
    the input set to inherit most of the functions.
    """
    def __init__(self, user_incar_settings, mode="Line",
                 constrain_total_magmom=False):
        """
        Args:
            user_incar_settings:
                A dict specify customized settings for INCAR.
                Must contain a NBANDS value, suggest to use
                1.2*(NBANDS from static run).
            mode:
                Line: Generate k-points along symmetry lines for bandstructure
                Uniform: Generate uniform k-points grids for DOS
        """
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self.mode = mode
        if mode not in ["Line", "Uniform"]:
            raise ValueError("Supported modes for NonSCF runs are 'Line' and "
                             "'Uniform'!")
        with open(os.path.join(module_dir, "MPVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MaterialsProject Static", json.load(f),
                constrain_total_magmom=constrain_total_magmom)
        self.user_incar_settings = user_incar_settings
        self.incar_settings.update(
            {"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001, "LCHARG": False,
             "LORBIT": 11, "LWAVE": False, "NSW": 0, "ISYM": 0, "ICHARG": 11})
        if mode == "Uniform":
            # Set smaller steps for DOS output
            self.incar_settings.update({"NEDOS": 601})
        if "NBANDS" not in user_incar_settings:
            raise KeyError("For NonSCF runs, NBANDS value from SC runs is "
                           "required!")
        else:
            self.incar_settings.update(user_incar_settings)

    def get_kpoints(self, structure, kpoints_density=1000):
        """
        Get a KPOINTS file for NonSCF calculation. In "Line" mode, kpoints are
        generated along high symmetry lines. In "Uniform" mode, kpoints are
        Gamma-centered mesh grid. Kpoints are written explicitly in both cases.

        Args:
            kpoints_density:
                kpoints density for the reciprocal cell of structure.
                Suggest to use a large kpoints_density.
                Might need to increase the default value when calculating
                metallic materials.
        """
        if self.mode == "Line":
            kpath = HighSymmKpath(structure)
            cart_k_points, k_points_labels = kpath.get_kpoints()
            frac_k_points = [kpath._prim_rec.get_fractional_coords(k)
                             for k in cart_k_points]
            return Kpoints(comment="Non SCF run along symmetry lines",
                           style="Reciprocal", num_kpts=len(frac_k_points),
                           kpts=frac_k_points, labels=k_points_labels,
                           kpts_weights=[1]*len(cart_k_points))
        else:
            num_kpoints = kpoints_density * \
                structure.lattice.reciprocal_lattice.volume
            kpoints = Kpoints.automatic_density(
                structure, num_kpoints * structure.num_sites)
            mesh = kpoints.kpts[0]
            ir_kpts = SymmetryFinder(structure, symprec=0.01)\
                .get_ir_reciprocal_mesh(mesh)
            kpts = []
            weights = []
            for k in ir_kpts:
                kpts.append(k[0])
                weights.append(int(k[1]))
            return Kpoints(comment="Non SCF run on uniform grid",
                           style="Reciprocal", num_kpts=len(ir_kpts),
                           kpts=kpts, kpts_weights=weights)

    @staticmethod
    def get_incar_settings(vasp_run, outcar=None):
        """
        Helper method to get necessary user_incar_settings from previous run.
        """
        # Turn off spin when magmom for every site is smaller than 0.02.
        if outcar and outcar.magnetization:
            site_magmom = np.array([i['tot'] for i in outcar.magnetization])
            ispin = 2 if np.any(site_magmom[np.abs(site_magmom) > 0.02]) else 1
        elif vasp_run.is_spin:
            ispin = 2
        else:
            ispin = 1
        nbands = int(np.ceil(vasp_run.to_dict["input"]["parameters"]["NBANDS"]
                             * 1.2))
        return {"ISPIN": ispin, "NBANDS": nbands}

    def get_incar(self, structure):
        incar = super(MPNonSCFVaspInputSet, self).get_incar(structure)
        incar.pop("MAGMOM", None)
        return incar

    @staticmethod
    def from_previous_vasp_run(previous_vasp_dir, output_dir='.',
                               mode="Uniform", user_incar_settings=None,
                               copy_chgcar=True, make_dir_if_not_present=True):
        """
        Generate a set of Vasp input files for NonSCF calculations from a
        directory of previous static Vasp run.

        Args:
            previous_vasp_dir:
                The directory contains the outputs(vasprun.xml and OUTCAR) of
                previous vasp run.
            output_dir:
                The directory to write the VASP input files for the NonSCF
                calculations. Default to write in the current directory.
            copy_chgcar:
                Default to copy CHGCAR from SC run
            make_dir_if_not_present:
                Set to True if you want the directory (and the whole path) to
                be created if it is not present.
        """
        try:
            vasp_run = Vasprun(os.path.join(previous_vasp_dir, "vasprun.xml"),
                               parse_dos=False, parse_eigen=None)
            outcar = Outcar(os.path.join(previous_vasp_dir, "OUTCAR"))
        except:
            traceback.format_exc()
            raise RuntimeError("Can't get valid results from previous run")

        #Get a Magmom-decorated structure
        structure = MPNonSCFVaspInputSet.get_structure(vasp_run, outcar)
        user_incar_settings = MPNonSCFVaspInputSet.get_incar_settings(vasp_run,
                                                                      outcar)
        mpnscfvip = MPNonSCFVaspInputSet(user_incar_settings, mode)
        mpnscfvip.write_input(structure, output_dir, make_dir_if_not_present)
        if copy_chgcar:
            try:
                shutil.copyfile(os.path.join(previous_vasp_dir, "CHGCAR"),
                                os.path.join(output_dir, "CHGCAR"))
            except Exception as e:
                traceback.format_exc()
                raise RuntimeError("Can't copy CHGCAR from SC run"+'\n'
                                   + str(e))


@deprecated(MPVaspInputSet)
class MaterialsProjectVaspInputSet(MPVaspInputSet):
    """
    A direct subclass of MPVaspInputSet (for backwards compatibility).

    .. deprecated:: v2.6.7

        Use MPVaspInputSet instead.
    """
    pass


@deprecated(MPGGAVaspInputSet)
class MaterialsProjectGGAVaspInputSet(MPGGAVaspInputSet):
    """
    A direct subclass of MPGGAVaspInputSet (for backwards compatibility).

    .. deprecated:: v2.6.7

        Use MPGGAVaspInputSet instead.
    """
    pass


def batch_write_vasp_input(structures, vasp_input_set, output_dir,
                           make_dir_if_not_present=True, subfolder=None,
                           sanitize=False, include_cif=False):
    """
    Batch write vasp input for a sequence of structures to
    output_dir, following the format output_dir/{group}/{formula}_{number}.

    Args:
        structures:
            Sequence of Structures.
        vasp_input_set:
            pymatgen.io.vaspio_set.VaspInputSet like object that creates
            vasp input files from structures
        output_dir:
            Directory to output files
        make_dir_if_not_present:
            Create the directory if not present. Defaults to True.
        subfolder:
            function to create subdirectory name from structure.
            Defaults to simply "formula_count".
        sanitize:
            Boolean indicating whether to sanitize the structure before
            writing the VASP input files. Sanitized output are generally easier
            for viewing and certain forms of analysis. Defaults to False.
        include_cif:
            Boolean indication whether to output a CIF as well. CIF files are
            generally better supported in visualization programs.
    """
    for i, s in enumerate(structures):
        formula = re.sub("\s+", "", s.formula)
        if subfolder is not None:
            subdir = subfolder(s)
            dirname = os.path.join(output_dir, subdir)
        else:
            dirname = os.path.join(output_dir, '{}_{}'.format(formula, i))
        if sanitize:
            s = s.copy(sanitize=True)
        vasp_input_set.write_input(
            s, dirname, make_dir_if_not_present=make_dir_if_not_present
        )
        if include_cif:
            from pymatgen.io.cifio import CifWriter
            writer = CifWriter(s)
            writer.write_file(os.path.join(dirname, "{}_{}.cif"
                                           .format(formula, i)))
