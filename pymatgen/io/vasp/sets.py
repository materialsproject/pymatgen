# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module defines the VaspInputSet abstract base class and a concrete
implementation for the parameters used by the Materials Project and the MIT
high throughput project.  The basic concept behind an input set is to specify
a scheme to generate a consistent set of VASP inputs from a structure
without further user intervention. This ensures comparability across
runs.
"""

__author__ = "Shyue Ping Ong, Wei Chen, Will Richards, Geoffroy Hautier"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Nov 16, 2011"

import os
import abc

import re
import traceback
import shutil
from functools import partial

import six
import numpy as np

from monty.serialization import loadfn

from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.serializers.json_coders import PMGSONable
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class AbstractVaspInputSet(six.with_metaclass(abc.ABCMeta, PMGSONable)):
    """
    Abstract base class representing a set of Vasp input parameters.
    The idea is that using a VaspInputSet, a complete set of input files
    (INPUT, KPOINTS, POSCAR and POTCAR) can be generated in an automated
    fashion for any structure.
    """

    @abc.abstractmethod
    def get_poscar(self, structure):
        """
        Returns Poscar from a structure.
        """
        return

    @abc.abstractmethod
    def get_kpoints(self, structure):
        """
        Returns Kpoints from a structure.

        Args:
            structure (Structure/IStructure): Structure to generate kpoints
                for.

        Returns:
            Kpoints object
        """
        return

    @abc.abstractmethod
    def get_incar(self, structure):
        """
        Returns Incar from a structure.

        Args:
            structure (Structure/IStructure): Structure to generate Incar for.

        Returns:
            Incar object
        """
        return

    @abc.abstractmethod
    def get_potcar(self, structure):
        """
        Returns Potcar from a structure.

        Args:
            structure (Structure/IStructure): Structure to generate potcar
                for.

        Returns:
            Potcar object
        """
        return

    @abc.abstractmethod
    def get_potcar_symbols(self, structure):
        """
        Returns list of POTCAR symbols from a structure.

        Args:
            structure (Structure/IStructure): Structure to generate potcar
                symbols for.

        Returns:
            List of POTCAR symbols
        """
        return

    def get_all_vasp_input(self, structure, generate_potcar=True):
        """
        Returns all input files as a dict of {filename: vasp object}

        Args:
            structure (Structure/IStructure): Structure to generate vasp
                input for.
            generate_potcar (bool): Set to False to generate a POTCAR.spec
                file instead of a POTCAR, which contains the POTCAR labels
                but not the actual POTCAR. Defaults to True.

        Returns:
            dict of {filename: file_as_string}, e.g., {'INCAR':'EDIFF=1e-4...'}
        """
        kpoints = self.get_kpoints(structure)
        incar = self.get_incar(structure)
        if np.product(kpoints.kpts) < 4 and incar.get("ISMEAR", 0) == -5:
            incar["ISMEAR"] = 0

        d = {'INCAR': incar,
             'KPOINTS': kpoints,
             'POSCAR': self.get_poscar(structure)}

        if generate_potcar:
            d['POTCAR'] = self.get_potcar(structure)
        else:
            d['POTCAR.spec'] = "\n".join(self.get_potcar_symbols(structure))
        return d

    def write_input(self, structure, output_dir,
                    make_dir_if_not_present=True, include_cif=False):
        """
        Writes a set of VASP input to a directory.

        Args:
            structure (Structure/IStructure): Structure to write VASP input
                files for.
            output_dir (str): Directory to output the VASP input files
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            include_cif (bool): Whether to write a CIF file in the output
                directory for easier opening by VESTA.
        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.get_all_vasp_input(structure).items():
            v.write_file(os.path.join(output_dir, k))
            if k == "POSCAR" and include_cif:
                v.structure.to(
                    filename=os.path.join(output_dir,
                                          "%s.cif" % v.structure.formula))


class DictVaspInputSet(AbstractVaspInputSet):
    """
    Concrete implementation of VaspInputSet that is initialized from a dict
    settings. This allows arbitrary settings to be input. In general,
    this is rarely used directly unless there is a source of settings in yaml
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

    Args:
        name (str): A name fo the input set.
        config_dict (dict): The config dictionary to use.
        hubbard_off (bool): Whether to turn off Hubbard U if it is specified in
            config_dict. Defaults to False, i.e., follow settings in
            config_dict.
        user_incar_settings (dict): User INCAR settings. This allows a user
            to override INCAR settings, e.g., setting a different MAGMOM for
            various elements or species.
        constrain_total_magmom (bool): Whether to constrain the total magmom
            (NUPDOWN in INCAR) to be the sum of the expected MAGMOM for all
            species. Defaults to False.
        sort_structure (bool): Whether to sort the structure (using the
            default sort order of electronegativity) before generating input
            files. Defaults to True, the behavior you would want most of the
            time. This ensures that similar atomic species are grouped
            together.
        ediff_per_atom (bool): Whether the EDIFF is specified on a per atom
            basis. This is generally desired, though for some calculations (
            e.g. NEB) this should be turned off (and an appropriate EDIFF
            supplied in user_incar_settings)
        potcar_functional (str): Functional to use. Default (None) is to use
            the functional in Potcar.DEFAULT_FUNCTIONAL. Valid values:
            "PBE", "LDA", "PW91", "LDA_US"
        force_gamma (bool): Force gamma centered kpoint generation. Default
            (False) is to use the Automatic Density kpoint scheme, which
            will use the Gamma centered generation scheme for hexagonal
            cells, and Monkhorst-Pack otherwise.
        reduce_structure (None/str): Before generating the input files,
            generate the reduced structure. Default (None), does not
            alter the structure. Valid values: None, "niggli", "LLL"
    """

    def __init__(self, name, config_dict, hubbard_off=False,
                 user_incar_settings=None,
                 constrain_total_magmom=False, sort_structure=True,
                 ediff_per_atom=True, potcar_functional=None,
                 force_gamma=False, reduce_structure=None):
        self.name = name
        self.potcar_settings = config_dict["POTCAR"]
        self.kpoints_settings = config_dict['KPOINTS']
        self.incar_settings = config_dict['INCAR']
        self.set_nupdown = constrain_total_magmom
        self.sort_structure = sort_structure
        self.ediff_per_atom = ediff_per_atom
        self.hubbard_off = hubbard_off
        self.potcar_functional = potcar_functional
        self.force_gamma = force_gamma
        self.reduce_structure = reduce_structure
        if hubbard_off:
            for k in list(self.incar_settings.keys()):
                if k.startswith("LDAU"):
                    del self.incar_settings[k]
        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

    def get_incar(self, structure):
        incar = Incar()
        if self.reduce_structure:
            structure = structure.get_reduced_structure(self.reduce_structure)
        if self.sort_structure:
            structure = structure.get_sorted_structure()
        comp = structure.composition
        elements = sorted([el for el in comp.elements if comp[el] > 0],
                          key=lambda e: e.X)
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
                if self.ediff_per_atom:
                    incar[key] = float(setting) * structure.num_sites
                else:
                    incar[key] = float(setting)
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
            for key in list(incar.keys()):
                if key.startswith('LDAU'):
                    del incar[key]

        if self.set_nupdown:
            nupdown = sum([mag if abs(mag) > 0.6 else 0
                           for mag in incar['MAGMOM']])
            incar['NUPDOWN'] = nupdown

        return incar

    def get_poscar(self, structure):
        if self.reduce_structure:
            structure = structure.get_reduced_structure(self.reduce_structure)
        if self.sort_structure:
            structure = structure.get_sorted_structure()
        return Poscar(structure)

    def get_potcar(self, structure, check_hash=False):
        if self.reduce_structure:
            structure = structure.get_reduced_structure(self.reduce_structure)
        if self.sort_structure:
            structure = structure.get_sorted_structure()
        if self.potcar_functional:
            p = Potcar(self.get_potcar_symbols(structure),
                          functional=self.potcar_functional)
        else:
            p = Potcar(self.get_potcar_symbols(structure))

        if check_hash:
            hash_check = [ps.hash == self.potcar_settings[ps.element][
                                                        'hash'] for ps in p]
            if all(hash_check):
                return p
            else:
                wrong_hashes = [p.symbols[i] for i, tf in enumerate(
                                                    hash_check) if not tf]
                raise ValueError("Potcars {} have different hashes "
                                 "than those specified in the config "
                                 "dictionary".format(wrong_hashes))
        else:
            return p

    def get_nelect(self, structure):
        """
        Gets the default number of electrons for a given structure.
        """
        n = 0
        for ps in self.get_potcar(structure):
            n += structure.composition[ps.element] * ps.ZVAL
        return n

    def get_potcar_symbols(self, structure):
        if self.reduce_structure:
            structure = structure.get_reduced_structure(self.reduce_structure)
        if self.sort_structure:
            structure = structure.get_sorted_structure()
        p = self.get_poscar(structure)
        elements = p.site_symbols
        potcar_symbols = []

        if isinstance(self.potcar_settings[elements[-1]], dict):
            for el in elements:
                potcar_symbols.append(self.potcar_settings[el]['symbol']
                                      if el in self.potcar_settings else el)
        else:
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
        if self.reduce_structure:
            structure = structure.get_reduced_structure(self.reduce_structure)
        if self.sort_structure:
            structure = structure.get_sorted_structure()

        # If grid_density is in the kpoints_settings use Kpoints.automatic_density
        if self.kpoints_settings.get('grid_density'):
            return Kpoints.automatic_density(
                structure, int(self.kpoints_settings['grid_density']),
                self.force_gamma)

        # If length is in the kpoints_settings use Kpoints.automatic
        elif self.kpoints_settings.get('length'):
            return Kpoints.automatic(self.kpoints_settings['length'])

        # Raise error. Unsure of which kpoint generation to use
        else:
            raise ValueError(
                "Invalid KPoint Generation algo : Supported Keys are "
                "grid_density: for Kpoints.automatic_density generation "
                "and length  : for Kpoints.automatic generation")

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

    def as_dict(self):
        config_dict = {
            "INCAR": self.incar_settings,
            "KPOINTS": self.kpoints_settings,
            "POTCAR": self.potcar_settings
        }
        return {
            "name": self.name,
            "config_dict": config_dict,
            "hubbard_off": self.hubbard_off,
            "constrain_total_magmom": self.set_nupdown,
            "sort_structure": self.sort_structure,
            "potcar_functional": self.potcar_functional,
            "@class": self.__class__.__name__,
            "@module": self.__class__.__module__,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(d["name"], d["config_dict"],
                   hubbard_off=d.get("hubbard_off", False),
                   constrain_total_magmom=d["constrain_total_magmom"],
                   sort_structure=d.get("sort_structure", True),
                   potcar_functional=d.get("potcar_functional", None))

    @staticmethod
    def from_file(name, filename, **kwargs):
        """
        Creates a DictVaspInputSet from a yaml/json file.

        Args:
            name (str): A name for the input set.
            filename (str): Path to a yaml/json file containing the settings.
            \*\*kwargs: Same kwargs as in the constructor.

        Returns:
            DictVaspInputSet
        """
        return DictVaspInputSet(name, loadfn(filename), **kwargs)


MITVaspInputSet = partial(DictVaspInputSet.from_file, "MIT",
                          os.path.join(MODULE_DIR, "MITVaspInputSet.yaml"))
"""
Standard implementation of VaspInputSet utilizing parameters in the MIT
High-throughput project.
The parameters are chosen specifically for a high-throughput project,
which means in general pseudopotentials with fewer electrons were chosen.

Please refer::

    A Jain, G. Hautier, C. Moore, S. P. Ong, C. Fischer, T. Mueller,
    K. A. Persson, G. Ceder. A high-throughput infrastructure for density
    functional theory calculations. Computational Materials Science,
    2011, 50(8), 2295-2310. doi:10.1016/j.commatsci.2011.02.023
"""

MITGGAVaspInputSet = partial(DictVaspInputSet.from_file, "MIT GGA",
                             os.path.join(MODULE_DIR, "MITVaspInputSet.yaml"),
                             hubbard_off=True)
"""
GGA (no U) version of MITVaspInputSet.
"""

MITHSEVaspInputSet = partial(
    DictVaspInputSet.from_file, "MIT HSE",
    os.path.join(MODULE_DIR, "MITHSEVaspInputSet.yaml"))
"""
Typical implementation of input set for a HSE run using MIT parameters.
"""


class MITNEBVaspInputSet(DictVaspInputSet):
    """
    Class for writing NEB inputs. Note that EDIFF is not on a per atom
    basis for this input set.

    Args:
        nimages (int): Number of NEB images (excluding start and ending
            structures).
        \*\*kwargs: Other kwargs supported by :class:`DictVaspInputSet`.
    """

    def __init__(self, nimages=8, user_incar_settings=None, **kwargs):
        #NEB specific defaults
        defaults = {'IMAGES': nimages, 'IBRION': 1, 'NFREE': 2, 'ISYM': 0,
                    'LORBIT': 0, 'LCHARG': False}
        if user_incar_settings:
            defaults.update(user_incar_settings)

        super(MITNEBVaspInputSet, self).__init__(
            "MIT NEB",
            loadfn(os.path.join(MODULE_DIR, "MITVaspInputSet.yaml")),
            user_incar_settings=defaults, ediff_per_atom=False, **kwargs)
        self.nimages = nimages

    def _process_structures(self, structures):
        """
        Remove any atom jumps across the cell
        """
        input_structures = structures
        structures = [input_structures[0]]
        for s in input_structures[1:]:
            prev = structures[-1]
            for i in range(len(s)):
                t = np.round(prev[i].frac_coords - s[i].frac_coords)
                if np.sum(t) > 0.5:
                    s.translate_sites([i], t, to_unit_cell=False)
            structures.append(s)
        return structures

    def write_input(self, structures, output_dir, make_dir_if_not_present=True,
                    write_cif=False):
        """
        NEB inputs has a special directory structure where inputs are in 00,
        01, 02, ....

        Args:
            structures ([Structure]): nimages + 2 structures (including
                start and end structures).
            output_dir (str): Directory to output the VASP input files
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            write_cif (bool): If true, writes a cif along with each POSCAR.
        """
        if len(structures) != self.incar_settings['IMAGES'] + 2:
            raise ValueError('incorrect number of structures')

        structures = self._process_structures(structures)

        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        s0 = structures[0]
        self.get_incar(s0).write_file(os.path.join(output_dir, 'INCAR'))
        self.get_kpoints(s0).write_file(os.path.join(output_dir, 'KPOINTS'))
        self.get_potcar(s0).write_file(os.path.join(output_dir, 'POTCAR'))
        for i, s in enumerate(structures):
            d = os.path.join(output_dir, str(i).zfill(2))
            if make_dir_if_not_present and not os.path.exists(d):
                os.makedirs(d)
            self.get_poscar(s).write_file(os.path.join(d, 'POSCAR'))
            if write_cif:
                s.to(filename=os.path.join(d, '{}.cif'.format(i)))

    def as_dict(self):
        d = super(MITNEBVaspInputSet, self).as_dict()
        d["nimages"] = self.nimages
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(user_incar_settings=d.get("user_incar_settings", None),
                   constrain_total_magmom=d["constrain_total_magmom"],
                   sort_structure=d.get("sort_structure", True),
                   hubbard_off=d.get("hubbard_off", False),
                   nimages=d["nimages"])


class MITMDVaspInputSet(DictVaspInputSet):
    """
    Class for writing a vasp md run. This DOES NOT do multiple stage
    runs.

    Args:
        start_temp (int): Starting temperature.
        end_temp (int): Final temperature.
        nsteps (int): Number of time steps for simulations. The NSW parameter.
        time_step (int): The time step for the simulation. The POTIM
            parameter. Defaults to 2fs.
        hubbard_off (bool): Whether to turn off Hubbard U. Defaults to
            *True* (different behavior from standard input sets) for MD runs.
        spin_polarized (bool): Whether to do spin polarized calculations.
            The ISPIN parameter. Defaults to False.
        sort_structure (bool): Whether to sort structure. Defaults to False
            (different behavior from standard input sets).
        **kwargs:
            Other kwargs supported by :class:`DictVaspInputSet`.
    """

    def __init__(self, start_temp, end_temp, nsteps, time_step=2,
                 hubbard_off=True, spin_polarized=False,
                 sort_structure=False, user_incar_settings=None,
                 **kwargs):

        #MD default settings
        defaults = {'TEBEG': start_temp, 'TEEND': end_temp, 'NSW': nsteps,
                    'EDIFF': 0.000001, 'LSCALU': False, 'LCHARG': False,
                    'LPLANE': False, 'LWAVE': True, 'ICHARG': 0, 'ISMEAR': 0,
                    'SIGMA': 0.05, 'NELMIN': 4, 'LREAL': True, 'BMIX': 1,
                    'MAXMIX': 20, 'NELM': 500, 'NSIM': 4, 'ISYM': 0,
                    'ISIF': 0, 'IBRION': 0, 'NBLOCK': 1, 'KBLOCK': 100,
                    'SMASS': 0, 'POTIM': time_step, 'PREC': 'Normal',
                    'ISPIN': 2 if spin_polarized else 1}

        #override default settings with user supplied settings
        if user_incar_settings:
            defaults.update(user_incar_settings)
        super(MITMDVaspInputSet, self).__init__(
            "MIT MD",
            loadfn(os.path.join(MODULE_DIR, "MITVaspInputSet.yaml")),
            hubbard_off=hubbard_off, sort_structure=sort_structure,
            user_incar_settings=defaults, **kwargs)

        self.start_temp = start_temp
        self.end_temp = end_temp
        self.nsteps = nsteps
        self.time_step = time_step
        self.spin_polarized = spin_polarized
        self.user_incar_settings = user_incar_settings or {}

        #use VASP default ENCUT
        if 'ENCUT' not in self.user_incar_settings:
            del self.incar_settings['ENCUT']

        if defaults['ISPIN'] == 1:
            del self.incar_settings['MAGMOM']

    def get_kpoints(self, structure):
        return Kpoints.gamma_automatic()

    def as_dict(self):
        d = super(MITMDVaspInputSet, self).as_dict()
        d.update({
            "start_temp": self.start_temp,
            "end_temp": self.end_temp,
            "nsteps": self.nsteps,
            "time_step": self.time_step,
            "spin_polarized": self.spin_polarized,
            "user_incar_settings": self.user_incar_settings
        })
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(start_temp=d["start_temp"], end_temp=d["end_temp"],
                   nsteps=d["nsteps"], time_step=d["time_step"],
                   hubbard_off=d.get("hubbard_off", False),
                   user_incar_settings=d["user_incar_settings"],
                   spin_polarized=d.get("spin_polarized", False),
                   constrain_total_magmom=d["constrain_total_magmom"],
                   sort_structure=d.get("sort_structure", True))


MPVaspInputSet = partial(DictVaspInputSet.from_file, "MP",
                         os.path.join(MODULE_DIR, "MPVaspInputSet.yaml"))
"""
Implementation of VaspInputSet utilizing parameters in the public
Materials Project. Typically, the pseudopotentials chosen contain more
electrons than the MIT parameters, and the k-point grid is ~50% more dense.
The LDAUU parameters are also different due to the different psps used,
which result in different fitted values.
"""

MPGGAVaspInputSet = partial(DictVaspInputSet.from_file, "MP GGA",
                            os.path.join(MODULE_DIR, "MPVaspInputSet.yaml"),
                            hubbard_off=True)
"""
Same as the MPVaspInput set, but the +U is enforced to be turned off.
"""


MPHSEVaspInputSet = partial(DictVaspInputSet.from_file, "MP HSE",
                            os.path.join(MODULE_DIR, "MPHSEVaspInputSet.yaml"))
"""
Same as the MPVaspInput set, but with HSE parameters.
"""


class MPStaticVaspInputSet(DictVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for static calculations that typically follow relaxation runs.
    It is recommended to use the static from_previous_run method to construct
    the input set to inherit most of the functions.

    Args:
        kpoints_density (int): kpoints density for the reciprocal cell of
            structure. Might need to increase the default value when
            calculating metallic materials.
        sym_prec (float): Tolerance for symmetry finding

    kwargs:
        hubbard_off (bool): Whether to turn off Hubbard U if it is specified in
            config_dict ("MP Static"). Defaults to False, i.e., follow settings
            in config_dict.
        user_incar_settings (dict): User INCAR settings. This allows a user
            to override INCAR settings, e.g., setting a different MAGMOM for
            various elements or species.
        constrain_total_magmom (bool): Whether to constrain the total magmom
            (NUPDOWN in INCAR) to be the sum of the expected MAGMOM for all
            species. Defaults to False.
        sort_structure (bool): Whether to sort the structure (using the
            default sort order of electronegativity) before generating input
            files. Defaults to True, the behavior you would want most of the
            time. This ensures that similar atomic species are grouped
            together.
        ediff_per_atom (bool): Whether the EDIFF is specified on a per atom
            basis.
    """

    def __init__(self, kpoints_density=90, sym_prec=0.1, **kwargs):
        super(MPStaticVaspInputSet, self).__init__(
            "MP Static",
            loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml")),
            **kwargs)
        self.incar_settings.update(
            {"IBRION": -1, "ISMEAR": -5, "LAECHG": True, "LCHARG": True,
             "LORBIT": 11, "LVHAR": True, "LWAVE": False, "NSW": 0,
             "ICHARG": 0, "EDIFF": 0.000001, "ALGO": "Normal"})
        self.kpoints_settings.update({"kpoints_density": kpoints_density})
        self.sym_prec = sym_prec

    def get_kpoints(self, structure, primitive_standard=False):
        """
        Get a KPOINTS file using the fully automated grid method. Uses
        Gamma centered meshes for hexagonal cells and Monk grids otherwise.

        Args:
            structure (Structure/IStructure): structure to get kpoints
            primitive_standard (Bool): whether the input structure is
            a primitive standardized cell
        """
        if not primitive_standard:
            structure = self.get_poscar(structure).structure
        self.kpoints_settings['grid_density'] = \
            self.kpoints_settings["kpoints_density"] * \
            structure.lattice.reciprocal_lattice.volume * \
            structure.num_sites
        return super(MPStaticVaspInputSet, self).get_kpoints(structure)

    def get_poscar(self, structure):
        """
        Get a POSCAR file with a primitive standardized cell of
        the giving structure.

        Args:
            structure (Structure/IStructure): structure to get POSCAR
        """
        sym_finder = SpacegroupAnalyzer(structure, symprec=self.sym_prec)
        return Poscar(sym_finder.get_primitive_standard_structure(False))

    @staticmethod
    def get_structure(vasp_run, outcar=None, initial_structure=False,
                      additional_info=False, sym_prec=0.1):
        """
        Process structure for static calculations from previous run.

        Args:
            vasp_run (Vasprun): Vasprun that contains the final structure
                from previous run.
            outcar (Outcar): Outcar that contains the magnetization info from
                previous run.
            initial_structure (bool): Whether to return the structure from
                previous run. Default is False.
            additional_info (bool):
                Whether to return additional symmetry info related to the
                structure. If True, return a list of the refined structure (
                conventional cell), the conventional standard structure,
                the symmetry dataset and symmetry operations of the
                structure (see SpacegroupAnalyzer doc for details).
            sym_prec (float): Tolerance for symmetry finding

        Returns:
            Returns the magmom-decorated structure that can be passed to get
            Vasp input files, e.g. get_kpoints.
        """
        if vasp_run.is_spin:
            if outcar and outcar.magnetization:
                magmom = {"magmom": [i['tot'] for i in outcar.magnetization]}
            else:
                magmom = {
                    "magmom": vasp_run.as_dict()['input']['parameters']
                    ['MAGMOM']}
        else:
            magmom = None
        structure = vasp_run.final_structure
        if magmom:
            structure = structure.copy(site_properties=magmom)
        sym_finder = SpacegroupAnalyzer(structure, symprec=sym_prec)
        if initial_structure:
            return structure
        elif additional_info:
            info = [sym_finder.get_refined_structure(),
                    sym_finder.get_conventional_standard_structure(False),
                    sym_finder.get_symmetry_dataset(),
                    sym_finder.get_symmetry_operations()]
            return [sym_finder.get_primitive_standard_structure(False),
                    info]
        else:
            return sym_finder.get_primitive_standard_structure(False)

    @staticmethod
    def from_previous_vasp_run(previous_vasp_dir, output_dir='.',
                               user_incar_settings=None,
                               make_dir_if_not_present=True,
                               kpoints_density=90, sym_prec=0.1):
        """
        Generate a set of Vasp input files for static calculations from a
        directory of previous Vasp run.

        Args:
            previous_vasp_dir (str): Directory containing the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            output_dir (str): Directory to write the VASP input files for
                the static calculations. Defaults to current directory.
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            kpoints_density (int): kpoints density for the reciprocal cell
                of structure. Might need to increase the default value when
                calculating metallic materials.
            sym_prec (float): Tolerance for symmetry finding
        """
        # Read input and output from previous run
        try:
            vasp_run = Vasprun(os.path.join(previous_vasp_dir, "vasprun.xml"),
                               parse_dos=False, parse_eigen=None)
            outcar = Outcar(os.path.join(previous_vasp_dir, "OUTCAR"))
            previous_incar = vasp_run.incar
            previous_kpoints = vasp_run.kpoints
        except:
            traceback.print_exc()
            raise RuntimeError("Can't get valid results from previous run. prev dir: {}".format(previous_vasp_dir))

        mpsvip = MPStaticVaspInputSet(kpoints_density=kpoints_density,
                                      sym_prec=sym_prec)
        structure = mpsvip.get_structure(vasp_run, outcar)

        mpsvip.write_input(structure, output_dir, make_dir_if_not_present)
        new_incar = mpsvip.get_incar(structure)

        # Use previous run INCAR and override necessary parameters
        previous_incar.update({"IBRION": -1, "ISMEAR": -5, "LAECHG": True,
                               "LCHARG": True, "LORBIT": 11, "LVHAR": True,
                               "LWAVE": False, "NSW": 0, "ICHARG": 0,
                               "ALGO": "Normal"})

        for incar_key in ["MAGMOM", "NUPDOWN"]:
            if new_incar.get(incar_key, None):
                previous_incar.update({incar_key: new_incar[incar_key]})
            else:
                previous_incar.pop(incar_key, None)

        # use new LDAUU when possible b/c the Poscar might have changed
        # representation
        if previous_incar.get('LDAU'):
            u = previous_incar.get('LDAUU', [])
            j = previous_incar.get('LDAUJ', [])
            if sum([u[x] - j[x] for x, y in enumerate(u)]) > 0:
                for tag in ('LDAUU', 'LDAUL', 'LDAUJ'):
                    previous_incar.update({tag: new_incar[tag]})
            # ensure to have LMAXMIX for GGA+U static run
            if "LMAXMIX" not in previous_incar:
                previous_incar.update({"LMAXMIX": new_incar["LMAXMIX"]})

        # Compare ediff between previous and staticinputset values,
        # choose the tighter ediff
        previous_incar.update({"EDIFF": min(previous_incar.get("EDIFF", 1),
                                            new_incar["EDIFF"])})

        # add user settings
        if user_incar_settings:
            previous_incar.update(user_incar_settings)
        previous_incar.write_file(os.path.join(output_dir, "INCAR"))

        # Perform checking on INCAR parameters
        if any([previous_incar.get("NSW", 0) != 0,
                previous_incar["IBRION"] != -1,
                previous_incar["LCHARG"] is not True,
               any([sum(previous_incar["LDAUU"]) <= 0,
                    previous_incar["LMAXMIX"] < 4])
               if previous_incar.get("LDAU") else False]):
            raise ValueError("Incompatible INCAR parameters!")

        # Prefer to use k-point scheme from previous run
        new_kpoints = mpsvip.get_kpoints(structure)
        if previous_kpoints.style[0] != new_kpoints.style[0]:
            if previous_kpoints.style[0] == "M" and \
                    SpacegroupAnalyzer(structure, 0.1).get_lattice_type() != \
                    "hexagonal":
                k_div = (kp + 1 if kp % 2 == 1 else kp
                         for kp in new_kpoints.kpts[0])
                Kpoints.monkhorst_automatic(k_div). \
                    write_file(os.path.join(output_dir, "KPOINTS"))
            else:
                Kpoints.gamma_automatic(new_kpoints.kpts[0]). \
                    write_file(os.path.join(output_dir, "KPOINTS"))
        else:
            new_kpoints.write_file(os.path.join(output_dir, "KPOINTS"))


class MPStaticDielectricDFPTVaspInputSet(DictVaspInputSet):
    """
    Using MP parameters to compute a static dielectric constant
    with DFPT. This includes the electronic and ionic contributions
    to the static dielectric constant.

    Args:
        user_incar_settings (dict): A dict specifying additional incar
            settings
        ionic: a boolean telling if we clamp the ions (False) or we
        add the ionic part to the dielectric constant (True default)
    """

    def __init__(self, user_incar_settings=None, ionic=True):
        super(MPStaticDielectricDFPTVaspInputSet, self).__init__(
            "Materials Project Static Dielectric DFPT",
            loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml")))
        self.user_incar_settings = user_incar_settings if \
            user_incar_settings is not None else {}
        self.incar_settings.update(self.user_incar_settings)
        if ionic:
            self.incar_settings.update(
                {"IBRION": 8, "LEPSILON": True, 'LREAL':False})
        else:
            self.incar_settings.update(
                {"LEPSILON": True, 'LREAL': False})
        if 'NPAR' in self.incar_settings:
            del self.incar_settings['NPAR']
        if 'NSW' in self.incar_settings:
            del self.incar_settings['NSW']


class MPBSHSEVaspInputSet(DictVaspInputSet):
    """
    Implementation of a VaspInputSet for HSE band structure computations
    Remember that HSE band structures cannot be non-self consistent. So, a band
    structure along syymetry lines for instance needs a uniform grid with
    appropriate weights + weight 0 path in reciprocal space

    Args:
        user_incar_settings(dict): A dict specifying additional incar
            settings
        added_kpoints: a list of kpoints (list of 3 number list) with weight 0
            added to the run. The k-points are in fractional coordinates
        kpoints_density: the kpoint density of the uniform grid used for the
            band structure run (By default: this is the same as in
            MPHSEVaspInputSet, i.e. 1000 / atom). Note that the uniform grid is
            always Gamma centered for now (this might be changed ?).
        mode: Line: Generate k-points along symmetry lines for
            bandstructure. Uniform: Generate uniform k-points grid

    """

    def __init__(self, user_incar_settings=None, added_kpoints=None, mode="Line",
                 kpoints_density=None, kpoints_line_density=20):
        super(MPBSHSEVaspInputSet, self).__init__(
            "Materials Project HSE Band Structure",
            loadfn(os.path.join(MODULE_DIR, "MPHSEVaspInputSet.yaml")))
        self.user_incar_settings = user_incar_settings if \
            user_incar_settings is not None else {}
        self.incar_settings.update(
            {"NSW": 0, "ISMEAR": 0, "SIGMA": 0.05, "ISYM": 0, "LCHARG": False})
        self.incar_settings.update(self.user_incar_settings)
        self.added_kpoints = added_kpoints if added_kpoints is not None else []
        self.mode = mode
        self.kpoints_density = (kpoints_density if kpoints_density is not None
                                else self.kpoints_settings['grid_density'])
        self.kpoints_line_density = kpoints_line_density

    def get_kpoints(self, structure):
        self.kpoints_settings['grid_density'] = self.kpoints_density
        grid = super(MPBSHSEVaspInputSet, self).get_kpoints(structure).kpts
        if self.mode == "Line":
            ir_kpts = SpacegroupAnalyzer(structure, symprec=0.1)\
                .get_ir_reciprocal_mesh(grid[0])
            kpath = HighSymmKpath(structure)
            frac_k_points, labels = kpath.get_kpoints(line_density=self.kpoints_line_density,
                                                      coords_are_cartesian=False)
            kpts = []
            weights = []
            all_labels = []
            for k in ir_kpts:
                kpts.append(k[0])
                weights.append(int(k[1]))
                all_labels.append(None)
            for k in range(len(frac_k_points)):
                kpts.append(frac_k_points[k])
                weights.append(0.0)
                all_labels.append(labels[k])
            return Kpoints(comment="HSE run along symmetry lines",
                           style="Reciprocal", num_kpts=len(kpts),
                           kpts=kpts, kpts_weights=weights, labels=all_labels)

        elif self.mode == "Uniform":
            ir_kpts = SpacegroupAnalyzer(structure, symprec=0.1)\
                .get_ir_reciprocal_mesh(grid[0])
            kpts = []
            weights = []
            for k in ir_kpts:
                kpts.append(k[0])
                weights.append(int(k[1]))
            for k in self.added_kpoints:
                kpts.append(k)
                weights.append(0.0)
            return Kpoints(comment="HSE run on uniform grid",
                           style="Reciprocal", num_kpts=len(kpts),
                           kpts=kpts, kpts_weights=weights)

    def as_dict(self):
        d = super(MPBSHSEVaspInputSet, self).as_dict()
        d['added_kpoints'] = self.added_kpoints
        d['mode'] = self.mode
        d['kpoints_density'] = self.kpoints_density
        d['kpoints_line_density'] = self.kpoints_line_density
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(user_incar_settings=d.get("user_incar_settings", None),
                   added_kpoints=d.get("added_kpoints", []),
                   mode=d.get("mode", "Line"),
                   kpoints_density=d.get("kpoints_density", None),
                   kpoints_line_density=d.get("kpoints_line_density", 20))


class MPNonSCFVaspInputSet(MPStaticVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for non self-consistent field (NonSCF) calculation that follows
    a static run to calculate bandstructure, density of states(DOS) and etc.
    It is recommended to use the NonSCF from_previous_run method to construct
    the input set to inherit most of the functions.

    Args:
        user_incar_settings (dict): A dict specify customized settings
            for INCAR. Must contain a NBANDS value, suggest to use
            1.2*(NBANDS from static run).
        mode: Line: Generate k-points along symmetry lines for
            bandstructure. Uniform: Generate uniform k-points
            grids for DOS.
        constrain_total_magmom (bool): Whether to constrain the total
            magmom (NUPDOWN in INCAR) to be the sum of the expected
            MAGMOM for all species. Defaults to False.
        kpoints_density (int): kpoints density for the reciprocal cell
            of structure. Might need to increase the default value when
            calculating metallic materials.
        kpoints_line_density (int): kpoints density to use in line-mode.
            Might need to increase the default value when calculating
            metallic materials.
        sort_structure (bool): Whether to sort structure. Defaults to
            False.
        sym_prec (float): Tolerance for symmetry finding
    """

    def __init__(self, user_incar_settings, mode="Line",
                 constrain_total_magmom=False, sort_structure=False,
                 kpoints_density=1000, sym_prec=0.1, kpoints_line_density=20):
        self.mode = mode
        self.sym_prec = sym_prec
        self.kpoints_line_density = kpoints_line_density
        if mode not in ["Line", "Uniform"]:
            raise ValueError("Supported modes for NonSCF runs are 'Line' and "
                             "'Uniform'!")
        DictVaspInputSet.__init__(self,
            "Materials Project Static",
            loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml")),
            constrain_total_magmom=constrain_total_magmom,
            sort_structure=sort_structure)
        self.user_incar_settings = user_incar_settings
        self.incar_settings.update(
            {"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001, "LCHARG": False,
             "LORBIT": 11, "LWAVE": False, "NSW": 0, "ISYM": 0, "ICHARG": 11})
        self.kpoints_settings.update({"kpoints_density": kpoints_density})
        if mode == "Uniform":
            # Set smaller steps for DOS output
            self.incar_settings.update({"NEDOS": 601})
        if "NBANDS" not in user_incar_settings:
            raise KeyError("For NonSCF runs, NBANDS value from SC runs is "
                           "required!")
        else:
            self.incar_settings.update(user_incar_settings)

    def get_kpoints(self, structure):
        """
        Get a KPOINTS file for NonSCF calculation. In "Line" mode, kpoints are
        generated along high symmetry lines. In "Uniform" mode, kpoints are
        Gamma-centered mesh grid. Kpoints are written explicitly in both cases.

        Args:
            structure (Structure/IStructure): structure to get Kpoints
        """
        if self.mode == "Line":
            kpath = HighSymmKpath(structure)
            frac_k_points, k_points_labels = kpath.get_kpoints(line_density=self.kpoints_line_density,
                                                               coords_are_cartesian=False)
            return Kpoints(comment="Non SCF run along symmetry lines",
                           style="Reciprocal", num_kpts=len(frac_k_points),
                           kpts=frac_k_points, labels=k_points_labels,
                           kpts_weights=[1] * len(frac_k_points))
        else:
            num_kpoints = self.kpoints_settings["kpoints_density"] * \
                structure.lattice.reciprocal_lattice.volume
            kpoints = Kpoints.automatic_density(
                structure, num_kpoints * structure.num_sites)
            mesh = kpoints.kpts[0]
            ir_kpts = SpacegroupAnalyzer(structure, symprec=self.sym_prec) \
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

        Args:
            vasp_run (Vasprun): Vasprun that contains the final
                structure from previous run.
            outcar (Outcar): Outcar that contains the magnetization info
                from previous run.

        """
        # Turn off spin when magmom for every site is smaller than 0.02.
        if outcar and outcar.magnetization:
            site_magmom = np.array([i['tot'] for i in outcar.magnetization])
            ispin = 2 if np.any(site_magmom[np.abs(site_magmom) > 0.02]) else 1
        elif vasp_run.is_spin:
            ispin = 2
        else:
            ispin = 1
        nbands = int(np.ceil(vasp_run.as_dict()["input"]["parameters"]["NBANDS"]
                             * 1.2))
        incar_settings = {"ISPIN": ispin, "NBANDS": nbands}
        for grid in ["NGX", "NGY", "NGZ"]:
            if vasp_run.incar.get(grid):
                incar_settings.update({grid: vasp_run.incar.get(grid)})
        return incar_settings

    def get_incar(self, structure):
        incar = super(MPNonSCFVaspInputSet, self).get_incar(structure)
        incar.pop("MAGMOM", None)
        return incar

    def get_poscar(self, structure, get_primitive_standard=False):
        """
        Get a POSCAR file of the giving structure.

        Args:
            structure (Structure/IStructure): structure to get POSCAR
            get_primitive_standard (bool): if convert the input structure to a
            primitive standard structure
        """
        if get_primitive_standard:
            sym_finder = SpacegroupAnalyzer(structure, symprec=self.sym_prec)
            return Poscar(sym_finder.get_primitive_standard_structure(False))
        else:
            return Poscar(structure)

    @staticmethod
    def from_previous_vasp_run(previous_vasp_dir, output_dir='.',
                               mode="Uniform", user_incar_settings=None,
                               copy_chgcar=True, make_dir_if_not_present=True,
                               kpoints_density=1000, kpoints_line_density=20):
        """
        Generate a set of Vasp input files for NonSCF calculations from a
        directory of previous static Vasp run.

        Args:
            previous_vasp_dir (str): The directory contains the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            output_dir (str): The directory to write the VASP input files
                for the NonSCF calculations. Default to write in the current
                directory.
            mode (str): Line: Generate k-points along symmetry lines for
                bandstructure. Uniform: Generate uniform k-points
                grids for DOS.
            user_incar_settings (dict): A dict specify customized settings
                for INCAR. Must contain a NBANDS value, suggest to use
                1.2*(NBANDS from static run).
            copy_chgcar (bool): Default to copy CHGCAR from SC run
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            kpoints_density (int): kpoints density for the reciprocal cell
                of structure. Might need to increase the default value when
                calculating metallic materials.
            kpoints_line_density (int): kpoints density to use in line-mode.
                Might need to increase the default value when calculating
                metallic materials.
        """
        user_incar_settings = user_incar_settings or {}

        try:
            vasp_run = Vasprun(os.path.join(previous_vasp_dir, "vasprun.xml"),
                               parse_dos=False, parse_eigen=None)
            outcar = Outcar(os.path.join(previous_vasp_dir, "OUTCAR"))
            previous_incar = vasp_run.incar
        except:
            traceback.print_exc()
            raise RuntimeError("Can't get valid results from previous run: {}"
                               .format(previous_vasp_dir))

        #Get a Magmom-decorated structure
        structure = MPNonSCFVaspInputSet.get_structure(vasp_run, outcar,
                                                       initial_structure=True)
        nscf_incar_settings = MPNonSCFVaspInputSet.get_incar_settings(vasp_run,
                                                                      outcar)
        mpnscfvip = MPNonSCFVaspInputSet(nscf_incar_settings, mode,
                                         kpoints_density=kpoints_density,
                                         kpoints_line_density=kpoints_line_density)
        mpnscfvip.write_input(structure, output_dir, make_dir_if_not_present)
        if copy_chgcar:
            try:
                shutil.copyfile(os.path.join(previous_vasp_dir, "CHGCAR"),
                                os.path.join(output_dir, "CHGCAR"))
            except Exception as e:
                traceback.print_exc()
                raise RuntimeError("Can't copy CHGCAR from SC run" + '\n'
                                   + str(e))

        #Overwrite necessary INCAR parameters from previous runs
        previous_incar.update({"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001,
                               "LCHARG": False, "LORBIT": 11, "LWAVE": False,
                               "NSW": 0, "ISYM": 0, "ICHARG": 11})
        previous_incar.update(nscf_incar_settings)
        previous_incar.update(user_incar_settings)
        previous_incar.pop("MAGMOM", None)
        previous_incar.write_file(os.path.join(output_dir, "INCAR"))

        # Perform checking on INCAR parameters
        if any([previous_incar.get("NSW", 0) != 0,
                previous_incar["IBRION"] != -1,
                previous_incar["ICHARG"] != 11,
               any([sum(previous_incar["LDAUU"]) <= 0,
                    previous_incar["LMAXMIX"] < 4])
               if previous_incar.get("LDAU") else False]):
            raise ValueError("Incompatible INCAR parameters!")


class MPOpticsNonSCFVaspInputSet(MPNonSCFVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for non self-consistent field (NonSCF) calculation with the computation
    of the dielectric function that follows a static run
    It is recommended to use the NonSCF from_previous_run method to construct
    the input set to inherit most of the functions.

    Args:
        user_incar_settings (dict): A dict specify customized settings
            for INCAR. Must contain a NBANDS value, suggest to use
            factor*(NBANDS from static run) with factor between 5 and 10.
        mode: Line: Generate k-points along symmetry lines for
            bandstructure. Uniform: Generate uniform k-points
            grids for DOS.
        constrain_total_magmom (bool): Whether to constrain the total
            magmom (NUPDOWN in INCAR) to be the sum of the expected
            MAGMOM for all species. Defaults to False.
        kpoints_density (int): kpoints density for the reciprocal cell
            of structure. Might need to increase the default value when
            calculating metallic materials.
        sort_structure (bool): Whether to sort structure. Defaults to
            False.
        sym_prec (float): Tolerance for symmetry finding
    """

    def __init__(self, user_incar_settings,
                 constrain_total_magmom=False, sort_structure=False,
                 kpoints_density=1000, sym_prec=0.1, nedos=2001):
        self.sym_prec = sym_prec
        self.nedos = nedos
        DictVaspInputSet.__init__(
            user_incar_settings, mode="Uniform",
            constrain_total_magmom=constrain_total_magmom,
            sort_structure=sort_structure,
            kpoints_density=kpoints_density, sym_prec=sym_prec)
        self.incar_settings.update({"NEDOS": nedos})
        self.incar_settings.update({"LOPTICS": True})

    @staticmethod
    def from_previous_vasp_run(previous_vasp_dir, output_dir='.',
                               user_incar_settings=None,
                               copy_chgcar=True, make_dir_if_not_present=True,
                               nbands_factor=5.0, nedos=2001):
        """
        Generate a set of Vasp input files for NonSCF calculations from a
        directory of previous static Vasp run.

        Args:
            previous_vasp_dir (str): The directory contains the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            output_dir (str): The directory to write the VASP input files
                for the NonSCF calculations. Default to write in the current
                directory.
            user_incar_settings (dict): A dict specify customized settings
                for INCAR. Must contain a NBANDS value, suggest to use
                1.2*(NBANDS from static run).
            copy_chgcar (bool): Default to copy CHGCAR from SC run
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            nbands_factor (float): Factor by which the number of bands is to be
                multiplied. Typical calculations of dielectric functions need a
                total number of bands of 5 to 10 times the number of valence
                bands.
        """
        user_incar_settings = user_incar_settings or {}

        try:
            vasp_run = Vasprun(os.path.join(previous_vasp_dir, "vasprun.xml"),
                               parse_dos=False, parse_eigen=False)
            outcar = Outcar(os.path.join(previous_vasp_dir, "OUTCAR"))
            previous_incar = vasp_run.incar
        except:
            traceback.print_exc()
            raise RuntimeError("Can't get valid results from previous run. prev dir: {}".format(previous_vasp_dir))

        #Get a Magmom-decorated structure
        structure = MPNonSCFVaspInputSet.get_structure(vasp_run, outcar,
                                                       initial_structure=True)
        nscf_incar_settings = MPNonSCFVaspInputSet.get_incar_settings(vasp_run,
                                                                      outcar)
        spin_band_settings = MPOpticsNonSCFVaspInputSet.get_ispin_nbands(
            vasp_run, outcar, nbands_factor=nbands_factor)
        mpnscfvip = MPNonSCFVaspInputSet(nscf_incar_settings, "Uniform")
        mpnscfvip.incar_settings.update(spin_band_settings)
        mpnscfvip.write_input(structure, output_dir, make_dir_if_not_present)
        if copy_chgcar:
            try:
                shutil.copyfile(os.path.join(previous_vasp_dir, "CHGCAR"),
                                os.path.join(output_dir, "CHGCAR"))
            except Exception as e:
                traceback.print_exc()
                raise RuntimeError("Can't copy CHGCAR from SC run" + '\n'
                                   + str(e))

        #Overwrite necessary INCAR parameters from previous runs
        previous_incar.update({"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001,
                               "LCHARG": False, "LORBIT": 11, "LWAVE": False,
                               "NSW": 0, "ISYM": 0, "ICHARG": 11,
                               "LOPTICS": True, "NEDOS": nedos})
        previous_incar.update(nscf_incar_settings)
        previous_incar.update(user_incar_settings)
        previous_incar.pop("MAGMOM", None)
        previous_incar.write_file(os.path.join(output_dir, "INCAR"))

        # Perform checking on INCAR parameters
        if any([previous_incar.get("NSW", 0) != 0,
                previous_incar["IBRION"] != -1,
                previous_incar["ICHARG"] != 11,
               any([sum(previous_incar["LDAUU"]) <= 0,
                    previous_incar["LMAXMIX"] < 4])
               if previous_incar.get("LDAU") else False]):
            raise ValueError("Incompatible INCAR parameters!")

    @staticmethod
    def get_ispin_nbands(vasp_run, outcar=None, nbands_factor=5.0):
        """
        Helper method to get necessary user_incar_settings from previous run.

        Args:
            vasp_run (Vasprun): Vasprun that contains the final
                structure from previous run.
            outcar (Outcar): Outcar that contains the magnetization info
                from previous run.

        """
        # Turn off spin when magmom for every site is smaller than 0.02.
        if outcar and outcar.magnetization:
            site_magmom = np.array([i['tot'] for i in outcar.magnetization])
            ispin = 2 if np.any(site_magmom[np.abs(site_magmom) > 0.02]) else 1
        elif vasp_run.is_spin:
            ispin = 2
        else:
            ispin = 1
        nbands = int(np.ceil(vasp_run.as_dict()["input"]["parameters"]["NBANDS"]
                             * nbands_factor))
        incar_settings = {"ISPIN": ispin, "NBANDS": nbands}
        for grid in ["NGX", "NGY", "NGZ"]:
            if vasp_run.incar.get(grid):
                incar_settings.update({grid: vasp_run.incar.get(grid)})
        return incar_settings


def batch_write_vasp_input(structures, vasp_input_set, output_dir,
                           make_dir_if_not_present=True, subfolder=None,
                           sanitize=False, include_cif=False):
    """
    Batch write vasp input for a sequence of structures to
    output_dir, following the format output_dir/{group}/{formula}_{number}.

    Args:
        structures ([Structure]): Sequence of Structures.
        vasp_input_set (VaspInputSet): VaspInputSet that creates
            vasp input files from structures
        output_dir (str): Directory to output files
        make_dir_if_not_present (bool): Create the directory if not present.
            Defaults to True.
        subfolder (callable): Function to create subdirectory name from
            structure. Defaults to simply "formula_count".
        sanitize (bool): Boolean indicating whether to sanitize the
            structure before writing the VASP input files. Sanitized output
            are generally easier for viewing and certain forms of analysis.
            Defaults to False.
        include_cif (bool): Whether to output a CIF as well. CIF files are
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
            s, dirname, make_dir_if_not_present=make_dir_if_not_present,
            include_cif=include_cif
        )
