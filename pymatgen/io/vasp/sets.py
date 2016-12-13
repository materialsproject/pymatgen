# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import abc
import re

import shutil
import warnings
from itertools import chain
from copy import deepcopy
try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path
import six
import numpy as np

from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from monty.json import MSONable
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.sites import PeriodicSite

"""
This module defines the VaspInputSet abstract base class and a concrete
implementation for the parameters developed and tested by the core team
of pymatgen, including the Materials Virtual Lab, Materials Project and the MIT
high throughput project.  The basic concept behind an input set is to specify
a scheme to generate a consistent set of VASP inputs from a structure
without further user intervention. This ensures comparability across
runs.

Read the following carefully before implementing new input sets:

1. 99% of what needs to be done can be done by specifying user_incar_settings
   to override some of the defaults of various input sets. Unless there is an
   extremely good reason to add a new set, DO NOT add one. E.g., if you want
   to turn the hubbard U off, just set "LDAU": False as a user_incar_setting.
2. All derivative input sets should inherit from one of the usual MPRelaxSet or
   MITRelaxSet, and proper superclass delegation should be used where possible.
   In particular, you are not supposed to implement your own as_dict or
   from_dict for derivative sets unless you know what you are doing.
   Improper overriding the as_dict and from_dict protocols is the major
   cause of implementation headaches. If you need an example, look at how the
   MPStaticSet or MPNonSCFSets are constructed.

The above are recommendations. The following are UNBREAKABLE rules:
1. All input sets must take in a structure or list of structures as the first
   argument.
2. user_incar_settings and user_kpoints_settings are absolute. Any new sets you
   implement must obey this. If a user wants to override your settings,
   you assume he knows what he is doing. Do not magically override user
   supplied settings. You can issue a warning if you think the user is wrong.
3. All input sets must save all supplied args and kwargs as instance variables.
   E.g., self.my_arg = my_arg and self.kwargs = kwargs in the __init__. This
   ensures the as_dict and from_dict work correctly.
"""

__author__ = "Shyue Ping Ong, Wei Chen, Will Richards, Geoffroy Hautier, Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "May 28 2016"


MODULE_DIR = Path(__file__).parent


class VaspInputSet(six.with_metaclass(abc.ABCMeta, MSONable)):
    """
    Base class representing a set of Vasp input parameters with a structure
    supplied as init parameters. Typically, you should not inherit from this
    class. Start from DictSet or MPRelaxSet or MITRelaxSet.
    """

    @abc.abstractproperty
    def incar(self):
        """Incar object"""
        pass

    @abc.abstractproperty
    def kpoints(self):
        """Kpoints object"""
        pass

    @abc.abstractproperty
    def poscar(self):
        """Poscar object"""
        pass

    @property
    def potcar_symbols(self):
        """
        List of POTCAR symbols.
        """
        elements = self.poscar.site_symbols
        potcar_symbols = []
        settings = self.config_dict["POTCAR"]

        if isinstance(settings[elements[-1]], dict):
            for el in elements:
                potcar_symbols.append(settings[el]['symbol']
                                      if el in settings else el)
        else:
            for el in elements:
                potcar_symbols.append(settings.get(el, el))

        return potcar_symbols

    @property
    def potcar(self):
        """
        Potcar object.
        """
        return Potcar(self.potcar_symbols, functional=self.potcar_functional)

    @property
    def all_input(self):
        """
        Returns all input files as a dict of {filename: vasp object}

        Returns:
            dict of {filename: object}, e.g., {'INCAR': Incar object, ...}
        """
        kpoints = self.kpoints
        incar = self.incar
        if np.product(kpoints.kpts) < 4 and incar.get("ISMEAR", 0) == -5:
            incar["ISMEAR"] = 0

        return {'INCAR': incar,
                'KPOINTS': kpoints,
                'POSCAR': self.poscar,
                'POTCAR': self.potcar}

    def write_input(self, output_dir,
                    make_dir_if_not_present=True, include_cif=False):
        """
        Writes a set of VASP input to a directory.

        Args:
            output_dir (str): Directory to output the VASP input files
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            include_cif (bool): Whether to write a CIF file in the output
                directory for easier opening by VESTA.
        """
        p = Path(output_dir)
        if make_dir_if_not_present and not p.exists():
            p.mkdir(parents=True)
        for k, v in self.all_input.items():
            v.write_file(str(p / k))
        if include_cif:
            s = self.all_input["POSCAR"].structure
            fname = str(p / ("%s.cif" % re.sub("\s", "", s.formula)))
            s.to(filename=fname)

    def as_dict(self, verbosity=2):
        d = MSONable.as_dict(self)
        if verbosity == 1:
            d.pop("structure", None)
        return d


class DictSet(VaspInputSet):
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
        structure (Structure): The Structure to create inputs for.
        config_dict (dict): The config dictionary to use.
        files_to_transfer (dict): A dictionary of {filename: filepath}. This
            allows the transfer of files from a previous calculation.
        user_incar_settings (dict): User INCAR settings. This allows a user
            to override INCAR settings, e.g., setting a different MAGMOM for
            various elements or species. Note that in the new scheme,
            ediff_per_atom and hubbard_u are no longer args. Instead, the
            config_dict supports EDIFF_PER_ATOM and EDIFF keys. The former
            scales with # of atoms, the latter does not. If both are
            present, EDIFF is preferred. To force such settings, just supply
            user_incar_settings={"EDIFF": 1e-5, "LDAU": False} for example.
        user_kpoints_settings (dict): Allow user to override kpoints setting by
            supplying a dict. E.g., {"reciprocal_density": 1000}. Default is
            None.
        constrain_total_magmom (bool): Whether to constrain the total magmom
            (NUPDOWN in INCAR) to be the sum of the expected MAGMOM for all
            species. Defaults to False.
        sort_structure (bool): Whether to sort the structure (using the
            default sort order of electronegativity) before generating input
            files. Defaults to True, the behavior you would want most of the
            time. This ensures that similar atomic species are grouped
            together.
        potcar_functional (str): Functional to use. Default (None) is to use
            the functional in Potcar.DEFAULT_FUNCTIONAL. Valid values:
            "PBE", "LDA", "PW91", "LDA_US"
        force_gamma (bool): Force gamma centered kpoint generation. Default
            (False) is to use the Automatic Density kpoint scheme, which
            will use the Gamma centered generation scheme for hexagonal
            cells, and Monkhorst-Pack otherwise.
        reduce_structure (None/str): Before generating the input files,
            generate the reduced structure. Default (None), does not
            alter the structure. Valid values: None, "niggli", "LLL".
    """

    def __init__(self, structure, config_dict,
                 files_to_transfer=None, user_incar_settings=None,
                 user_kpoints_settings=None,
                 constrain_total_magmom=False, sort_structure=True,
                 potcar_functional="PBE", force_gamma=False,
                 reduce_structure=None):
        if reduce_structure:
            structure = structure.get_reduced_structure(reduce_structure)
        if sort_structure:
            structure = structure.get_sorted_structure()
        self.structure = structure
        self.config_dict = deepcopy(config_dict)
        self.files_to_transfer = files_to_transfer or {}
        self.constrain_total_magmom = constrain_total_magmom
        self.sort_structure = sort_structure
        self.potcar_functional = potcar_functional
        self.force_gamma = force_gamma
        self.reduce_structure = reduce_structure
        self.user_incar_settings = user_incar_settings or {}
        self.user_kpoints_settings = user_kpoints_settings

    @property
    def incar(self):
        settings = dict(self.config_dict["INCAR"])
        settings.update(self.user_incar_settings)
        structure = self.structure
        incar = Incar()
        comp = structure.composition
        elements = sorted([el for el in comp.elements if comp[el] > 0],
                          key=lambda e: e.X)
        most_electroneg = elements[-1].symbol
        poscar = Poscar(structure)
        hubbard_u = settings.get("LDAU", False)
        for k, v in settings.items():
            if k == "MAGMOM":
                mag = []
                for site in structure:
                    if hasattr(site, 'magmom'):
                        mag.append(site.magmom)
                    elif hasattr(site.specie, 'spin'):
                        mag.append(site.specie.spin)
                    elif str(site.specie) in v:
                        mag.append(v.get(str(site.specie)))
                    else:
                        mag.append(v.get(site.specie.symbol, 0.6))
                incar[k] = mag
            elif k in ('LDAUU', 'LDAUJ', 'LDAUL'):
                if hubbard_u:
                    if hasattr(structure[0], k.lower()):
                        m = dict([(site.specie.symbol, getattr(site, k.lower()))
                                  for site in structure])
                        incar[k] = [m[sym] for sym in poscar.site_symbols]
                    elif most_electroneg in v.keys():
                        incar[k] = [v[most_electroneg].get(sym, 0)
                                    for sym in poscar.site_symbols]
                    else:
                        incar[k] = [0] * len(poscar.site_symbols)
            elif k.startswith("EDIFF") and k != "EDIFFG":
                if "EDIFF" not in settings and k == "EDIFF_PER_ATOM":
                    incar["EDIFF"] = float(v) * structure.num_sites
                else:
                    incar["EDIFF"] = float(settings["EDIFF"])
            else:
                incar[k] = v

        has_u = hubbard_u and sum(incar['LDAUU']) > 0
        if has_u:
            # modify LMAXMIX if LSDA+U and you have d or f electrons
            # note that if the user explicitly sets LMAXMIX in settings it will
            # override this logic.
            if 'LMAXMIX' not in settings.keys():
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

        if self.constrain_total_magmom:
            nupdown = sum([mag if abs(mag) > 0.6 else 0
                           for mag in incar['MAGMOM']])
            incar['NUPDOWN'] = nupdown

        return incar

    @property
    def poscar(self):
        return Poscar(self.structure)

    @property
    def nelect(self):
        """
        Gets the default number of electrons for a given structure.
        """
        n = 0
        for ps in self.potcar:
            n += self.structure.composition[ps.element] * ps.ZVAL
        return n

    @property
    def kpoints(self):
        """
        Writes out a KPOINTS file using the fully automated grid method. Uses
        Gamma centered meshes  for hexagonal cells and Monk grids otherwise.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.
        """
        settings = self.user_kpoints_settings or self.config_dict["KPOINTS"]

        # If grid_density is in the kpoints_settings use
        # Kpoints.automatic_density
        if settings.get('grid_density'):
            return Kpoints.automatic_density(
                self.structure, int(settings['grid_density']),
                self.force_gamma)

        # If reciprocal_density is in the kpoints_settings use
        # Kpoints.automatic_density_by_vol
        elif settings.get('reciprocal_density'):
            return Kpoints.automatic_density_by_vol(
                self.structure, int(settings['reciprocal_density']),
                self.force_gamma)

        # If length is in the kpoints_settings use Kpoints.automatic
        elif settings.get('length'):
            return Kpoints.automatic(settings['length'])

        # Raise error. Unsure of which kpoint generation to use
        else:
            raise ValueError(
                "Invalid KPoint Generation algo : Supported Keys are "
                "grid_density: for Kpoints.automatic_density generation, "
                "reciprocal_density: for KPoints.automatic_density_by_vol "
                "generation, and length  : for Kpoints.automatic generation")

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__class__.__name__

    def write_input(self, output_dir,
                    make_dir_if_not_present=True, include_cif=False):
        super(DictSet, self).write_input(
            output_dir=output_dir,
            make_dir_if_not_present=make_dir_if_not_present,
            include_cif=include_cif)
        p = Path(output_dir)
        for k, v in self.files_to_transfer.items():
            shutil.copy(v, str(p / k))


class MITRelaxSet(DictSet):
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
    CONFIG = loadfn(str(MODULE_DIR / "MITRelaxSet.yaml"))

    def __init__(self, structure, **kwargs):
        super(MITRelaxSet, self).__init__(
            structure, MITRelaxSet.CONFIG, **kwargs)
        self.kwargs = kwargs


class MPRelaxSet(DictSet):
    """
    Implementation of VaspInputSet utilizing parameters in the public
    Materials Project. Typically, the pseudopotentials chosen contain more
    electrons than the MIT parameters, and the k-point grid is ~50% more dense.
    The LDAUU parameters are also different due to the different psps used,
    which result in different fitted values.
    """
    CONFIG = loadfn(str(MODULE_DIR / "MPRelaxSet.yaml"))

    def __init__(self, structure, **kwargs):
        super(MPRelaxSet, self).__init__(
            structure, MPRelaxSet.CONFIG, **kwargs)
        self.kwargs = kwargs


class MPHSERelaxSet(DictSet):
    """
    Same as the MPRelaxSet, but with HSE parameters.
    """
    CONFIG = loadfn(str(MODULE_DIR / "MPHSERelaxSet.yaml"))

    def __init__(self, structure, **kwargs):
        super(MPHSERelaxSet, self).__init__(
            structure, MPHSERelaxSet.CONFIG, **kwargs)
        self.kwargs = kwargs


class MPStaticSet(MPRelaxSet):

    def __init__(self, structure, prev_incar=None, prev_kpoints=None,
                 lepsilon=False, reciprocal_density=100, **kwargs):
        """
        Run a static calculation.

        Args:
            structure (Structure): Structure from previous run.
            prev_incar (Incar): Incar file from previous run.
            prev_kpoints (Kpoints): Kpoints from previous run.
            lepsilon (bool): Whether to add static dielectric calculation
            reciprocal_density (int): For static calculations,
                we usually set the reciprocal density by voluyme. This is a
                convenience arg to change that, rather than using
                user_kpoints_settings. Defaults to 100.
            \*\*kwargs: kwargs supported by MPRelaxSet.
        """
        super(MPStaticSet, self).__init__(structure, **kwargs)
        self.prev_incar = prev_incar
        self.prev_kpoints = prev_kpoints
        self.reciprocal_density = reciprocal_density
        self.structure = structure
        self.kwargs = kwargs
        self.lepsilon = lepsilon

    @property
    def incar(self):
        parent_incar = super(MPStaticSet, self).incar
        incar = Incar(self.prev_incar) if self.prev_incar is not None else \
            Incar(parent_incar)

        incar.update(
            {"IBRION": -1, "ISMEAR": -5, "LAECHG": True, "LCHARG": True,
             "LORBIT": 11, "LVHAR": True, "LWAVE": False, "NSW": 0,
             "ICHARG": 0, "ALGO": "Normal"})

        if self.lepsilon:
            incar["IBRION"] = 8
            incar["LEPSILON"] = True
            # Note that DFPT calculations MUST unset NSW. NSW = 0 will fail
            # to output ionic.
            incar.pop("NSW", None)
            incar.pop("NPAR", None)

        for k in ["MAGMOM", "NUPDOWN"] + list(self.kwargs.get(
                "user_incar_settings", {}).keys()):
            # For these parameters as well as user specified settings, override
            # the incar settings.
            if parent_incar.get(k, None) is not None:
                incar[k] = parent_incar[k]
            else:
                incar.pop(k, None)

        # use new LDAUU when possible b/c the Poscar might have changed
        # representation
        if incar.get('LDAU'):
            u = incar.get('LDAUU', [])
            j = incar.get('LDAUJ', [])
            if sum([u[x] - j[x] for x, y in enumerate(u)]) > 0:
                for tag in ('LDAUU', 'LDAUL', 'LDAUJ'):
                    incar.update({tag: parent_incar[tag]})
            # ensure to have LMAXMIX for GGA+U static run
            if "LMAXMIX" not in incar:
                incar.update({"LMAXMIX": parent_incar["LMAXMIX"]})

        # Compare ediff between previous and staticinputset values,
        # choose the tighter ediff
        incar["EDIFF"] = min(incar.get("EDIFF", 1), parent_incar["EDIFF"])
        return incar

    @property
    def kpoints(self):
        self.config_dict["KPOINTS"]["reciprocal_density"] = \
            self.reciprocal_density
        kpoints = super(MPStaticSet, self).kpoints
        # Prefer to use k-point scheme from previous run
        if self.prev_kpoints and self.prev_kpoints.style != kpoints.style:
            if self.prev_kpoints.style == Kpoints.supported_modes.Monkhorst:
                k_div = [kp + 1 if kp % 2 == 1 else kp
                         for kp in kpoints.kpts[0]]
                kpoints = Kpoints.monkhorst_automatic(k_div)
            else:
                kpoints = Kpoints.gamma_automatic(kpoints.kpts[0])
        return kpoints

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, standardize=False, sym_prec=0.1,
                       international_monoclinic=True, reciprocal_density=100,
                       small_gap_multiply=None, **kwargs):
        """
        Generate a set of Vasp input files for static calculations from a
        directory of previous Vasp run.

        Args:
            prev_calc_dir (str): Directory containing the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            standardize (float): Whether to standardize to a primitive
                standard cell. Defaults to False.
            sym_prec (float): Tolerance for symmetry finding. If not 0,
                the final structure from the previous run will be symmetrized
                to get a primitive standard cell. Set to 0 if you don't want
                that.
            international_monoclinic (bool): Whether to use international
                    convention (vs Curtarolo) for monoclinic. Defaults True.
            reciprocal_density (int): density of k-mesh by reciprocal
                                    volume (defaults to 100)
            small_gap_multiply ([float, float]): If the gap is less than
                1st index, multiply the default reciprocal_density by the 2nd
                index.
            \*\*kwargs: All kwargs supported by MPStaticSet,
                other than prev_incar and prev_structure and prev_kpoints which
                are determined from the prev_calc_dir.
        """
        vasprun, outcar = get_vasprun_outcar(prev_calc_dir)

        prev_incar = vasprun.incar
        prev_kpoints = vasprun.kpoints

        # We will make a standard structure for the given symprec.
        prev_structure = get_structure_from_prev_run(
            vasprun, outcar, sym_prec=standardize and sym_prec,
            international_monoclinic=international_monoclinic)

        # multiply the reciprocal density if needed:
        if small_gap_multiply:
            gap = vasprun.eigenvalue_band_properties[0]
            if gap <= small_gap_multiply[0]:
                reciprocal_density = reciprocal_density * small_gap_multiply[1]

        return MPStaticSet(
            structure=prev_structure, prev_incar=prev_incar,
            prev_kpoints=prev_kpoints,
            reciprocal_density=reciprocal_density, **kwargs)


class MPHSEBSSet(MPHSERelaxSet):

    def __init__(self, structure, user_incar_settings=None, added_kpoints=None,
                 mode="Uniform", reciprocal_density=None,
                 kpoints_line_density=20, **kwargs):
        """
        Implementation of a VaspInputSet for HSE band structure computations.
        Remember that HSE band structures must be self-consistent in VASP. A
        band structure along symmetry lines for instance needs BOTH a uniform
        grid with appropriate weights AND a path along the lines with weight 0.

        Thus, the "Uniform" mode is just like regular static SCF but allows
        adding custom kpoints (e.g., corresponding to known VBM/CBM) to the
        uniform grid that have zero weight (e.g., for better gap estimate).

        The "Line" mode is just like Uniform mode, but additionally adds
        k-points along symmetry lines with zero weight.

        Args:
            structure (Structure): Structure to compute
            user_incar_settings (dict): A dict specifying additional incar
                settings
            added_kpoints (list): a list of kpoints (list of 3 number list)
                added to the run. The k-points are in fractional coordinates
            mode (str): "Line" - generate k-points along symmetry lines for
                bandstructure. "Uniform" - generate uniform k-points grid
            reciprocal_density (int): k-point density to use for uniform mesh
            kpoints_line_density (int): k-point density for high symmetry lines
            **kwargs (dict): Any other parameters to pass into DictVaspInputSet

        """
        super(MPHSEBSSet, self).__init__(structure, **kwargs)
        self.structure = structure
        self.user_incar_settings = user_incar_settings or {}
        self.config_dict["INCAR"].update(
            {"NSW": 0, "ISMEAR": 0, "SIGMA": 0.05, "ISYM": 3, "LCHARG": False, "NELMIN": 5})
        self.added_kpoints = added_kpoints if added_kpoints is not None else []
        self.mode = mode
        self.reciprocal_density = reciprocal_density or \
            self.kpoints_settings['reciprocal_density']
        self.kpoints_line_density = kpoints_line_density

    @property
    def kpoints(self):
        kpts = []
        weights = []
        all_labels = []

        # for both modes, include the Uniform mesh w/standard weights
        grid = Kpoints.automatic_density_by_vol(self.structure,
                                                self.reciprocal_density).kpts
        ir_kpts = SpacegroupAnalyzer(self.structure, symprec=0.1)\
            .get_ir_reciprocal_mesh(grid[0])
        for k in ir_kpts:
            kpts.append(k[0])
            weights.append(int(k[1]))
            all_labels.append(None)

        # for both modes, include any user-added kpoints w/zero weight
        for k in self.added_kpoints:
            kpts.append(k)
            weights.append(0.0)
            all_labels.append("user-defined")

        # for line mode only, add the symmetry lines w/zero weight
        if self.mode.lower() == "line":
            kpath = HighSymmKpath(self.structure)
            frac_k_points, labels = kpath.get_kpoints(
                line_density=self.kpoints_line_density,
                coords_are_cartesian=False)

            for k in range(len(frac_k_points)):
                kpts.append(frac_k_points[k])
                weights.append(0.0)
                all_labels.append(labels[k])

        comment = "HSE run along symmetry lines" if self.mode.lower() == "line" \
            else "HSE run on uniform grid"

        return Kpoints(comment=comment,
                       style=Kpoints.supported_modes.Reciprocal,
                       num_kpts=len(kpts), kpts=kpts, kpts_weights=weights,
                       labels=all_labels)

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, mode="gap",
                       reciprocal_density=50, copy_chgcar=True, **kwargs):
        """
        Generate a set of Vasp input files for HSE calculations from a
        directory of previous Vasp run. if mode=="gap", it explicitly adds VBM and CBM
        of the prev. run to the k-point list of this run.

        Args:
            prev_calc_dir (str): Directory containing the outputs
                (vasprun.xml and OUTCAR) of previous vasp run.
            mode (str): Either "uniform", "gap" or "line"
            reciprocal_density (int): density of k-mesh
            copy_chgcar (bool): whether to copy CHGCAR of previous run
            \*\*kwargs: All kwargs supported by MPHSEBSStaticSet,
                other than prev_structure which is determined from the previous
                calc dir.
        """

        vasprun, outcar = get_vasprun_outcar(prev_calc_dir)

        # note: don't standardize the cell because we want to retain k-points
        prev_structure = get_structure_from_prev_run(vasprun, outcar,
                                                     sym_prec=0)

        added_kpoints = []
        bs = vasprun.get_band_structure()
        if mode.lower() == "gap":
            vbm, cbm = bs.get_vbm()["kpoint"], bs.get_cbm()["kpoint"]
            if vbm:
                added_kpoints.append(vbm.frac_coords)
            if cbm:
                added_kpoints.append(cbm.frac_coords)

        files_to_transfer = {}
        if copy_chgcar:
            chgcars = sorted(Path(prev_calc_dir).glob("CHGCAR*"))
            if chgcars:
                files_to_transfer["CHGCAR"] = str(chgcars[-1])

        return MPHSEBSSet(
            structure=prev_structure,
            added_kpoints=added_kpoints, reciprocal_density=reciprocal_density,
            mode=mode, files_to_transfer=files_to_transfer, **kwargs)


class MPNonSCFSet(MPRelaxSet):

    def __init__(self, structure, prev_incar=None,
                 mode="line", nedos=601, reciprocal_density=100, sym_prec=0.1,
                 kpoints_line_density=20, optics=False, **kwargs):
        """
        Init a MPNonSCFSet. Typically, you would use the classmethod
        from_prev_calc to initialize from a previous SCF run.

        Args:
            structure (Structure): Structure to compute
            prev_incar (Incar): Incar file from previous run.
            mode (str): Line or Uniform mode supported.
            nedos (int): nedos parameter. Default to 601.
            reciprocal_density (int): density of k-mesh by reciprocal
                                    volume (defaults to 100)
            sym_prec (float): Symmetry precision (for Uniform mode).
            kpoints_line_density (int): Line density for Line mode.
            optics (bool): whether to add dielectric function
            \*\*kwargs: kwargs supported by MPVaspInputSet.
        """
        super(MPNonSCFSet, self).__init__(structure, **kwargs)
        self.prev_incar = prev_incar
        self.kwargs = kwargs
        self.nedos = nedos
        self.reciprocal_density = reciprocal_density
        self.sym_prec = sym_prec
        self.kpoints_line_density = kpoints_line_density
        self.optics = optics
        self.mode = mode.lower()

        if self.mode.lower() not in ["line", "uniform"]:
            raise ValueError("Supported modes for NonSCF runs are 'Line' and "
                             "'Uniform'!")
        if (self.mode.lower() != "uniform" or nedos < 2000) and optics:
            warnings.warn("It is recommended to use Uniform mode with a high "
                          "NEDOS for optics calculations.")

    @property
    def incar(self):
        incar = super(MPNonSCFSet, self).incar
        if self.prev_incar is not None:
            incar.update({k: v for k, v in self.prev_incar.items()
                         if k not in self.kwargs.get("user_incar_settings",
                                                     {})})

        # Overwrite necessary INCAR parameters from previous runs
        incar.update({"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001,
                      "LCHARG": False, "LORBIT": 11, "LWAVE": False,
                      "NSW": 0, "ISYM": 0, "ICHARG": 11})

        if self.mode.lower() == "uniform":
            # Set smaller steps for DOS output
            incar["NEDOS"] = self.nedos

        if self.optics:
            incar["LOPTICS"] = True

        incar.pop("MAGMOM", None)

        return incar

    @property
    def kpoints(self):
        if self.mode == "line":
            kpath = HighSymmKpath(self.structure)
            frac_k_points, k_points_labels = kpath.get_kpoints(
                line_density=self.kpoints_line_density,
                coords_are_cartesian=False)
            kpoints = Kpoints(
                comment="Non SCF run along symmetry lines",
                style=Kpoints.supported_modes.Reciprocal,
                num_kpts=len(frac_k_points),
                kpts=frac_k_points, labels=k_points_labels,
                kpts_weights=[1] * len(frac_k_points))
        else:
            kpoints = Kpoints.automatic_density_by_vol(self.structure,
                                                       self.reciprocal_density)
            mesh = kpoints.kpts[0]
            ir_kpts = SpacegroupAnalyzer(
                self.structure,
                symprec=self.sym_prec).get_ir_reciprocal_mesh(mesh)
            kpts = []
            weights = []
            for k in ir_kpts:
                kpts.append(k[0])
                weights.append(int(k[1]))
            kpoints = Kpoints(comment="Non SCF run on uniform grid",
                              style=Kpoints.supported_modes.Reciprocal,
                              num_kpts=len(ir_kpts),
                              kpts=kpts, kpts_weights=weights)
        return kpoints

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, copy_chgcar=True,
                       nbands_factor=1.2, standardize=False, sym_prec=0.1,
                       international_monoclinic=True, reciprocal_density=100,
                       kpoints_line_density=20, small_gap_multiply=None,
                       **kwargs):
        """
        Generate a set of Vasp input files for NonSCF calculations from a
        directory of previous static Vasp run.

        Args:
            prev_calc_dir (str): The directory contains the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            copy_chgcar: Whether to copy the old CHGCAR. Defaults to True.
            nbands_factor (float): Multiplicative factor for NBANDS. Choose a
                higher number if you are doing an LOPTICS calculation.
            standardize (float): Whether to standardize to a primitive
                standard cell. Defaults to False.
            sym_prec (float): Tolerance for symmetry finding. If not 0,
                the final structure from the previous run will be symmetrized
                to get a primitive standard cell. Set to 0 if you don't want
                that.
            international_monoclinic (bool): Whether to use international
                convention (vs Curtarolo) for monoclinic. Defaults True.
            reciprocal_density (int): density of k-mesh by reciprocal
                volume in uniform mode (defaults to 100)
            kpoints_line_density (int): density of k-mesh in line mode
                (defaults to 20)
            small_gap_multiply ([float, float]): If the gap is less than
                1st index, multiply the default reciprocal_density by the 2nd
                index.
            \*\*kwargs: All kwargs supported by MPNonSCFSet,
                other than structure, prev_incar and prev_chgcar which
                are determined from the prev_calc_dir.
        """
        vasprun, outcar = get_vasprun_outcar(prev_calc_dir)

        incar = vasprun.incar
        # Get a Magmom-decorated structure
        structure = get_structure_from_prev_run(
            vasprun, outcar, sym_prec=standardize and sym_prec,
            international_monoclinic=international_monoclinic)
        # Turn off spin when magmom for every site is smaller than 0.02.
        if outcar and outcar.magnetization:
            site_magmom = np.array([i['tot'] for i in outcar.magnetization])
            ispin = 2 if np.any(site_magmom[np.abs(site_magmom) > 0.02]) else 1
        elif vasprun.is_spin:
            ispin = 2
        else:
            ispin = 1
        nbands = int(np.ceil(vasprun.parameters["NBANDS"] * nbands_factor))
        incar.update({"ISPIN": ispin, "NBANDS": nbands})

        files_to_transfer = {}
        if copy_chgcar:
            chgcars = sorted(Path(prev_calc_dir).glob("CHGCAR*"))
            if chgcars:
                files_to_transfer["CHGCAR"] = str(chgcars[-1])

        # multiply the reciprocal density if needed:
        if small_gap_multiply:
            gap = vasprun.eigenvalue_band_properties[0]
            if gap <= small_gap_multiply[0]:
                reciprocal_density = reciprocal_density * small_gap_multiply[1]
                kpoints_line_density = kpoints_line_density * \
                    small_gap_multiply[1]

        return MPNonSCFSet(structure=structure, prev_incar=incar,
                           reciprocal_density=reciprocal_density,
                           kpoints_line_density=kpoints_line_density,
                           files_to_transfer=files_to_transfer, **kwargs)


class MPSOCSet(MPStaticSet):

    def __init__(self, structure, saxis=(0, 0, 1), prev_incar=None,
                 reciprocal_density=100, **kwargs):
        """
        Init a MPSOCSet.

        Args:
            structure (Structure): the structure must have the 'magmom' site
                property and each magnetic moment value must have 3
                components. eg:- magmom = [[0,0,2], ...]
            saxis (tuple): magnetic moment orientation
            prev_incar (Incar): Incar file from previous run.
            reciprocal_density (int): density of k-mesh by reciprocal
                                    volume (defaults to 100)
            \*\*kwargs: kwargs supported by MPVaspInputSet.
        """
        if not hasattr(structure[0], "magmom") and \
                not isinstance(structure[0].magmom, list):
            raise ValueError("The structure must have the 'magmom' site "
                             "property and each magnetic moment value must have 3 "
                             "components. eg:- magmom = [0,0,2]")
        self.saxis = saxis
        super(MPSOCSet, self).__init__(
            structure, prev_incar=prev_incar,
            reciprocal_density=reciprocal_density, **kwargs)

    @property
    def incar(self):
        incar = super(MPSOCSet, self).incar
        if self.prev_incar is not None:
            incar.update({k: v for k, v in self.prev_incar.items()
                         if k not in self.kwargs.get("user_incar_settings",
                                                     {})})

        # Overwrite necessary INCAR parameters from previous runs
        incar.update({"ISYM": -1, "LSORBIT": "T", "ICHARG": 11,
                      "SAXIS": list(self.saxis)})

        return incar

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, copy_chgcar=True,
                       nbands_factor=1.2, standardize=False, sym_prec=0.1,
                       international_monoclinic=True, reciprocal_density=100,
                       small_gap_multiply=None, **kwargs):
        """
        Generate a set of Vasp input files for SOC calculations from a
        directory of previous static Vasp run. SOC calc requires all 3
        components for MAGMOM for each atom in the structure.

        Args:
            prev_calc_dir (str): The directory contains the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            copy_chgcar: Whether to copy the old CHGCAR. Defaults to True.
            nbands_factor (float): Multiplicative factor for NBANDS. Choose a
                higher number if you are doing an LOPTICS calculation.
            standardize (float): Whether to standardize to a primitive
                standard cell. Defaults to False.
            sym_prec (float): Tolerance for symmetry finding. If not 0,
                the final structure from the previous run will be symmetrized
                to get a primitive standard cell. Set to 0 if you don't want
                that.
            international_monoclinic (bool): Whether to use international
                convention (vs Curtarolo) for monoclinic. Defaults True.
            reciprocal_density (int): density of k-mesh by reciprocal
                volume (defaults to 100)
            small_gap_multiply ([float, float]): If the gap is less than
                1st index, multiply the default reciprocal_density by the 2nd
                index.
            \*\*kwargs: All kwargs supported by MPSOCSet,
                other than structure, prev_incar and prev_chgcar which
                are determined from the prev_calc_dir.
        """
        vasprun, outcar = get_vasprun_outcar(prev_calc_dir)

        incar = vasprun.incar
        # Get a magmom-decorated structure
        structure = get_structure_from_prev_run(
            vasprun, outcar, sym_prec=standardize and sym_prec,
            international_monoclinic=international_monoclinic)
        # override magmom if provided
        if kwargs.get("magmom", None):
            structure = structure.copy(
                site_properties={"magmom": kwargs["magmom"]})
            kwargs.pop("magmom", None)
        # magmom has to be 3D for SOC calculation.
        if hasattr(structure[0], "magmom"):
            if not isinstance(structure[0].magmom, list):
                structure = structure.copy(site_properties={
                    "magmom": [[0, 0, site.magmom] for site in structure]})
        else:
            raise ValueError("Neither the previous structure has mamgom "
                             "property nor magmom provided")

        nbands = int(np.ceil(vasprun.parameters["NBANDS"] * nbands_factor))
        incar.update({"NBANDS": nbands})

        files_to_transfer = {}
        if copy_chgcar:
            chgcars = sorted(Path(prev_calc_dir).glob("CHGCAR*"))
            if chgcars:
                files_to_transfer["CHGCAR"] = str(chgcars[-1])

        # multiply the reciprocal density if needed:
        if small_gap_multiply:
            gap = vasprun.eigenvalue_band_properties[0]
            if gap <= small_gap_multiply[0]:
                reciprocal_density = reciprocal_density * small_gap_multiply[1]

        return MPSOCSet(structure, prev_incar=incar,
                        files_to_transfer=files_to_transfer,
                        reciprocal_density=reciprocal_density, **kwargs)


class MVLElasticSet(MPRelaxSet):
    """
    MVL denotes VASP input sets that are implemented by the Materials Virtual
    Lab (http://www.materialsvirtuallab.org) for various research.

    This input set is used to calculate elastic constants in VASP. It is used
    in the following work::

        Z. Deng, Z. Wang, I.-H. Chu, J. Luo, S. P. Ong.
        “Elastic Properties of Alkali Superionic Conductor Electrolytes
        from First Principles Calculations”, J. Electrochem. Soc.
        2016, 163(2), A67-A74. doi: 10.1149/2.0061602jes

    To read the elastic constants, you may use the Outcar class which parses the
    elastic constants.

    Args:
        scale (float): POTIM parameter. The default of 0.015 is usually fine,
            but some structures may require a smaller step.
        user_incar_settings (dict): A dict specifying additional incar
            settings.
    """

    def __init__(self, structure, potim=0.015, **kwargs):
        super(MVLElasticSet, self).__init__(structure, **kwargs)
        self.config_dict["INCAR"].update({"IBRION": 6, "NFREE": 2,
                                          "POTIM": potim})
        self.config_dict["INCAR"].pop("NPAR", None)


class MVLSlabSet(MPRelaxSet):
    """
    Class for writing a set of slab vasp runs,
    including both slabs (along the c direction) and orient unit cells (bulk),
    to ensure the same KPOINTS, POTCAR and INCAR criterion.

    Args:
        k_product: default to 50, kpoint number * length for a & b directions,
            also for c direction in bulk calculations
        bulk (bool): Set to True for bulk calculation. Defaults to False.
        **kwargs:
            Other kwargs supported by :class:`DictSet`.
    """
    def __init__(self, structure, k_product=50, bulk=False, **kwargs):
        super(MVLSlabSet, self).__init__(structure, **kwargs)
        self.structure = structure
        self.k_product = k_product
        self.bulk = bulk
        self.kwargs = kwargs

        slab_incar = {"EDIFF": 1e-6, "EDIFFG": -0.01, "ENCUT": 400,
                      "ISMEAR": 0, "SIGMA": 0.05, "ISIF": 3}
        if not self.bulk:
            slab_incar["ISIF"] = 2
            slab_incar["AMIN"] = 0.01
            slab_incar["AMIX"] = 0.2
            slab_incar["BMIX"] = 0.001
            slab_incar["NELMIN"] = 8

        self.config_dict["INCAR"].update(slab_incar)

    @property
    def kpoints(self):
        """
        k_product, default to 50, is kpoint number * length for a & b
            directions, also for c direction in bulk calculations
        Automatic mesh & Gamma is the default setting.
        """

        # To get input sets, the input structure has to has the same number
        # of required parameters as a Structure object (ie. 4). Slab
        # attributes aren't going to affect the VASP inputs anyways so
        # converting the slab into a structure should not matter

        kpt = super(MVLSlabSet, self).kpoints
        kpt.comment = "Automatic mesh"
        kpt.style = 'Gamma'

        # use k_product to calculate kpoints, k_product = kpts[0][0] * a
        abc = self.structure.lattice.abc
        kpt_calc = [int(self.k_product/abc[0]+0.5),
                    int(self.k_product/abc[1]+0.5), 1]
        self.kpt_calc = kpt_calc
        # calculate kpts (c direction) for bulk. (for slab, set to 1)
        if self.bulk:
            kpt_calc[2] = int(self.k_product/abc[2]+0.5)

        kpt.kpts[0] = kpt_calc

        return kpt


class MITNEBSet(MITRelaxSet):
    """
    Class for writing NEB inputs. Note that EDIFF is not on a per atom
    basis for this input set.

    Args:
        unset_encut (bool): Whether to unset ENCUT.
        \*\*kwargs: Other kwargs supported by :class:`DictSet`.
    """

    def __init__(self, structures, unset_encut=False, **kwargs):
        if len(structures) < 3:
            raise ValueError("You need at least 3 structures for an NEB.")
        kwargs["sort_structure"] = False
        super(MITNEBSet, self).__init__(structures[0], **kwargs)
        self.structures = self._process_structures(structures)
        self.unset_encut = False
        if unset_encut:
            self.config_dict["INCAR"].pop("ENCUT", None)

        if "EDIFF" not in self.config_dict["INCAR"]:
            self.config_dict["INCAR"]["EDIFF"] = self.config_dict[
                "INCAR"].pop("EDIFF_PER_ATOM")

        # NEB specific defaults
        defaults = {'IMAGES': len(structures) - 2, 'IBRION': 1, 'ISYM': 0,
                    'LCHARG': False, "LDAU": False}
        self.config_dict["INCAR"].update(defaults)

    @property
    def poscar(self):
        return Poscar(self.structures[0])

    @property
    def poscars(self):
        return [Poscar(s) for s in self.structures]

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
                if np.any(np.abs(t)>0.5):
                    s.translate_sites([i], t, to_unit_cell=False)
            structures.append(s)
        return structures

    def write_input(self, output_dir, make_dir_if_not_present=True,
                    write_cif=False, write_path_cif=False,
                    write_endpoint_inputs=False):
        """
        NEB inputs has a special directory structure where inputs are in 00,
        01, 02, ....

        Args:
            output_dir (str): Directory to output the VASP input files
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            write_cif (bool): If true, writes a cif along with each POSCAR.
            write_path_cif (bool): If true, writes a cif for each image.
            write_endpoint_inputs (bool): If true, writes input files for
                running endpoint calculations.
        """
        path = Path(output_dir)
        if make_dir_if_not_present and not path.exists():
            path.mkdir(parents=True)
        self.incar.write_file(str(path / 'INCAR'))
        self.kpoints.write_file(str(path / 'KPOINTS'))
        self.potcar.write_file(str(path / 'POTCAR'))

        for i, p in enumerate(self.poscars):
            d = path / str(i).zfill(2)
            if not d.exists():
                d.mkdir()
            p.write_file(str(d / 'POSCAR'))
            if write_cif:
                p.structure.to(filename=str(d / '{}.cif'.format(i)))
        if write_endpoint_inputs:
            end_point_param = MITRelaxSet(
                self.structures[0],
                user_incar_settings=self.user_incar_settings)

            for image in ['00', str(len(self.structures) - 1).zfill(2)]:
                end_point_param.incar.write_file(str(path / image / 'INCAR'))
                end_point_param.kpoints.write_file(str(path / image / 'KPOINTS'))
                end_point_param.potcar.write_file(str(path / image / 'POTCAR'))
        if write_path_cif:
            sites = set()
            l = self.structures[0].lattice
            for site in chain(*(s.sites for s in self.structures)):
                sites.add(PeriodicSite(site.species_and_occu, site.frac_coords, l))
            nebpath = Structure.from_sites(sorted(sites))
            nebpath.to(filename=str(path / 'path.cif'))


class MITMDSet(MITRelaxSet):
    """
    Class for writing a vasp md run. This DOES NOT do multiple stage
    runs.

    Args:
        start_temp (int): Starting temperature.
        end_temp (int): Final temperature.
        nsteps (int): Number of time steps for simulations. The NSW parameter.
        time_step (int): The time step for the simulation. The POTIM
            parameter. Defaults to 2fs.
        spin_polarized (bool): Whether to do spin polarized calculations.
            The ISPIN parameter. Defaults to False.
        sort_structure (bool): Whether to sort structure. Defaults to False
            (different behavior from standard input sets).
        \*\*kwargs: Other kwargs supported by :class:`DictSet`.
    """

    def __init__(self, structure, start_temp, end_temp, nsteps, time_step=2,
                 spin_polarized=False, **kwargs):

        # MD default settings
        defaults = {'TEBEG': start_temp, 'TEEND': end_temp, 'NSW': nsteps,
                    'EDIFF_PER_ATOM': 0.000001, 'LSCALU': False,
                    'LCHARG': False,
                    'LPLANE': False, 'LWAVE': True, 'ISMEAR': 0,
                    'NELMIN': 4, 'LREAL': True, 'BMIX': 1,
                    'MAXMIX': 20, 'NELM': 500, 'NSIM': 4, 'ISYM': 0,
                    'ISIF': 0, 'IBRION': 0, 'NBLOCK': 1, 'KBLOCK': 100,
                    'SMASS': 0, 'POTIM': time_step, 'PREC': 'Normal',
                    'ISPIN': 2 if spin_polarized else 1,
                    "LDAU": False}

        super(MITMDSet, self).__init__(structure, **kwargs)

        self.start_temp = start_temp
        self.end_temp = end_temp
        self.nsteps = nsteps
        self.time_step = time_step
        self.spin_polarized = spin_polarized
        self.kwargs = kwargs

        # use VASP default ENCUT
        self.config_dict["INCAR"].pop('ENCUT', None)

        if defaults['ISPIN'] == 1:
            self.config_dict["INCAR"].pop('MAGMOM', None)
        self.config_dict["INCAR"].update(defaults)

    @property
    def kpoints(self):
        return Kpoints.gamma_automatic()


def get_vasprun_outcar(path, parse_dos=True, parse_eigen=True):
    p = Path(path)
    vruns = list(p.glob("vasprun.xml*"))
    outcars = list(p.glob("OUTCAR*"))

    if len(vruns) == 0 or len(outcars) == 0:
        raise ValueError(
            "Unable to get vasprun.xml/OUTCAR from prev calculation in %s" %
            path)
    vsfile_fullpath = p / "vasprun.xml"
    outcarfile_fullpath = p / "OUTCAR"
    vsfile = vsfile_fullpath if vsfile_fullpath in vruns else sorted(vruns)[-1]
    outcarfile = outcarfile_fullpath if outcarfile_fullpath in outcars else sorted(outcars)[-1]
    return Vasprun(str(vsfile), parse_dos=parse_dos, parse_eigen=parse_eigen), \
           Outcar(str(outcarfile))


def get_structure_from_prev_run(vasprun, outcar=None, sym_prec=0.1,
                                international_monoclinic=True):
    """
    Process structure from previous run.

    Args:
        vasprun (Vasprun): Vasprun that contains the final structure
            from previous run.
        outcar (Outcar): Outcar that contains the magnetization info from
            previous run.
        sym_prec (float): Tolerance for symmetry finding for standardization. If
            no standardization is desired, set to 0 or a False.
        international_monoclinic (bool): Whether to use international
            convention (vs Curtarolo) for monoclinic. Defaults True.

    Returns:
        Returns the magmom-decorated structure that can be passed to get
        Vasp input files, e.g. get_kpoints.
    """
    structure = vasprun.final_structure

    site_properties = {}
    # magmom
    if vasprun.is_spin:
        if outcar and outcar.magnetization:
            site_properties.update({"magmom": [i['tot']
                                               for i in outcar.magnetization]})
        else:
            site_properties.update({"magmom": vasprun.parameters['MAGMOM']})
    # ldau
    if vasprun.parameters.get("LDAU", False):
        for k in ("LDAUU", "LDAUJ", "LDAUL"):
            vals = vasprun.incar[k]
            m = {}
            l = []
            s = 0
            for site in structure:
                if site.specie.symbol not in m:
                    m[site.specie.symbol] = vals[s]
                    s += 1
                l.append(m[site.specie.symbol])
            if len(l) == len(structure):
                site_properties.update({k.lower(): l})
            else:
                raise ValueError("length of list {} not the same as"
                                 "structure".format(l))

    structure = structure.copy(site_properties=site_properties)

    if sym_prec:
        sym_finder = SpacegroupAnalyzer(structure, symprec=sym_prec)
        new_structure = sym_finder.get_primitive_standard_structure(
            international_monoclinic=international_monoclinic)
        # the primitive structure finding has had several bugs in the past
        # defend through validation
        vpa_old = structure.volume / structure.num_sites
        vpa_new = new_structure.volume / new_structure.num_sites
        if abs(vpa_old - vpa_new) / vpa_old > 0.02:
            raise ValueError(
                "Standardizing cell failed! VPA old: {}, VPA new: {}".format(
                    vpa_old, vpa_new))
        sm = StructureMatcher()
        if not sm.fit(structure, new_structure):
            raise ValueError(
                "Standardizing cell failed! Old structure doesn't match new.")
        structure = new_structure

    return structure


def batch_write_input(structures, vasp_input_set=MPRelaxSet, output_dir=".",
                      make_dir_if_not_present=True, subfolder=None,
                      sanitize=False, include_cif=False, **kwargs):
    """
    Batch write vasp input for a sequence of structures to
    output_dir, following the format output_dir/{group}/{formula}_{number}.

    Args:
        structures ([Structure]): Sequence of Structures.
        vasp_input_set (VaspInputSet): VaspInputSet class that creates
            vasp input files from structures. Note that a class should be
            supplied. Defaults to MPRelaxSet.
        output_dir (str): Directory to output files. Defaults to current
            directory ".".
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
        \*\*kwargs: Additional kwargs are passed to the vasp_input_set class in
            addition to structure.
    """
    for i, s in enumerate(structures):
        formula = re.sub("\s+", "", s.formula)
        d = Path(output_dir)
        if subfolder is not None:
            subdir = subfolder(s)
            d /= subdir
        else:
            d /= '{}_{}'.format(formula, i)
        if sanitize:
            s = s.copy(sanitize=True)
        v = vasp_input_set(s, **kwargs)
        v.write_input(str(d), make_dir_if_not_present=make_dir_if_not_present,
                      include_cif=include_cif)