# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import os
import abc

import re
import traceback
import shutil
from functools import partial
from glob import glob
import warnings

import six
import numpy as np

from monty.serialization import loadfn

from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from monty.json import MSONable, MontyDecoder
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.analysis.structure_matcher import StructureMatcher

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


MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

from pymatgen.io.vasp.sets_deprecated import *


class VaspInputSet(six.with_metaclass(abc.ABCMeta, MSONable)):
    """
    Base class representing a set of DERIVED Vasp input parameters,
    which means that a previous run, including a structure is supplied as
    init parameters. The different with AbstractVaspInputSet is that these
    are simpler because a structure is not supplied to each of the abstract
    methods but is instead already provided. See examples below.
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

    @abc.abstractproperty
    def potcar(self):
        """Potcar object"""
        pass

    def get_all_vasp_input(self):
        """
        Returns all input files as a dict of {filename: vasp object}

        Returns:
            dict of {filename: file_as_string}, e.g., {'INCAR':'EDIFF=1e-4...'}
        """
        return {'INCAR': self.incar,
                'KPOINTS': self.kpoints,
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
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.get_all_vasp_input().items():
            v.write_file(os.path.join(output_dir, k))
            if k == "POSCAR" and include_cif:
                v.structure.to(
                    filename=os.path.join(output_dir,
                                          "%s.cif" % v.structure.formula))

    def as_dict(self):
        d = MSONable.as_dict(self)
        if hasattr(self, "kwargs"):
            d.update(**self.kwargs)
        return d

    @classmethod
    def from_dict(cls, d):
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in d.items()
                   if not k.startswith("@")}
        return cls(**decoded)


class MPStaticSet(VaspInputSet):

    def __init__(self, structure, prev_incar=None, prev_kpoints=None,
                 lepsilon=False, reciprocal_density=100, **kwargs):
        """
        Init a MPStaticSet. Typically, you would use the classmethod
        from_prev_calc instead.

        Args:
            structure (Structure): Structure from previous run.
            prev_incar (Incar): Incar file from previous run.
            prev_kpoints (Kpoints): Kpoints from previous run.
            reciprocal_density (int): density of k-mesh by reciprocal
                volume (defaults to 100)
            \*\*kwargs: kwargs supported by MPVaspInputSet.
        """

        self.prev_incar = prev_incar
        self.prev_kpoints = prev_kpoints
        self.reciprocal_density = reciprocal_density
        self.structure = structure
        self.kwargs = kwargs
        self.lepsilon = lepsilon
        self.parent_vis = MPVaspInputSet(**self.kwargs)

    @property
    def incar(self):

        parent_incar = self.parent_vis.get_incar(self.structure)
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
        if self.parent_vis.kpoints_settings.get("grid_density"):
            del self.parent_vis.kpoints_settings["grid_density"]
        self.parent_vis.kpoints_settings["reciprocal_density"] = \
            self.reciprocal_density
        kpoints = self.parent_vis.get_kpoints(self.structure)

        # Prefer to use k-point scheme from previous run
        if self.prev_kpoints and self.prev_kpoints.style != kpoints.style:
            if self.prev_kpoints.style == Kpoints.supported_modes.Monkhorst:
                k_div = [kp + 1 if kp % 2 == 1 else kp
                         for kp in kpoints.kpts[0]]
                kpoints = Kpoints.monkhorst_automatic(k_div)
            else:
                kpoints = Kpoints.gamma_automatic(kpoints.kpts[0])
        return kpoints

    @property
    def poscar(self):
        return Poscar(self.structure)

    @property
    def potcar(self):
        return self.parent_vis.get_potcar(self.structure)

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
            small_gap_multiply ([float, float]) - if the gap is less than 1st index,
                                multiply the default reciprocal_density by the 2nd index
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


class MPNonSCFSet(VaspInputSet):

    def __init__(self, structure, prev_incar=None, prev_chgcar_path=None,
                 mode="line", nedos=601, reciprocal_density=100, sym_prec=0.1,
                 kpoints_line_density=20, optics=False, **kwargs):
        """
        Init a MPNonSCFSet. Typically, you would use the classmethod
        from_prev_calc instead.

        Args:
            prev_incar (Incar): Incar file from previous run.
            prev_structure (Structure): Structure from previous run.
            prev_chgcar_path (str): Path to CHGCAR from previous run.
            mode (str): Line or Uniform mode supported.
            nedos (int): nedos parameter. Default to 601.
            reciprocal_density (int): density of k-mesh by reciprocal
                                    volume (defaults to 100)
            sym_prec (float): Symmetry precision (for Uniform mode).
            kpoints_line_density (int): Line density for Line mode.
            \*\*kwargs: kwargs supported by MPVaspInputSet.
        """
        self.structure = structure
        self.prev_incar = prev_incar
        self.prev_chgcar_path = prev_chgcar_path
        self.kwargs = kwargs

        self.nedos = nedos
        self.reciprocal_density = reciprocal_density
        self.sym_prec = sym_prec
        self.kpoints_line_density = kpoints_line_density
        self.optics = optics
        self.mode = mode.lower()

        if self.mode not in ["line", "uniform"]:
            raise ValueError("Supported modes for NonSCF runs are 'Line' and "
                             "'Uniform'!")
        if (self.mode != "uniform" or nedos < 2000) and optics:
            warnings.warn("It is recommended to use Uniform mode with a high "
                          "NEDOS for optics calculations.")

        self.parent_vis = MPVaspInputSet(**kwargs)


    @property
    def incar(self):
        incar = self.parent_vis.get_incar(self.structure)
        if self.prev_incar is not None:
            incar.update({k: v for k, v in self.prev_incar.items()
                         if k not in self.kwargs.get("user_incar_settings",
                                                     {})})

        # Overwrite necessary INCAR parameters from previous runs
        incar.update({"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001,
                      "LCHARG": False, "LORBIT": 11, "LWAVE": False,
                      "NSW": 0, "ISYM": 0, "ICHARG": 11})

        if self.mode == "uniform":
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

    @property
    def poscar(self):
        return Poscar(self.structure)

    @property
    def potcar(self):
        return self.parent_vis.get_potcar(self.structure)

    def write_input(self, output_dir,
                    make_dir_if_not_present=True, include_cif=False):
        super(MPNonSCFSet, self).write_input(output_dir,
            make_dir_if_not_present=make_dir_if_not_present,
            include_cif=include_cif)
        if self.prev_chgcar_path:
            shutil.copy(self.prev_chgcar_path,
                        os.path.join(output_dir, "CHGCAR"))

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, copy_chgcar=True,
                       nbands_factor=1.2, standardize=False, sym_prec=0.1,
                       international_monoclinic=True, reciprocal_density=100,
                       small_gap_multiply=None, **kwargs):
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
                volume (defaults to 100)
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

        chgcar_path = None
        if copy_chgcar:
            chgcars = glob(os.path.join(prev_calc_dir, "CHGCAR*"))
            if chgcars:
                chgcar_path = sorted(chgcars)[-1]

        # multiply the reciprocal density if needed:
        if small_gap_multiply:
            gap = vasprun.eigenvalue_band_properties[0]
            if gap <= small_gap_multiply[0]:
                reciprocal_density = reciprocal_density * small_gap_multiply[1]

        return MPNonSCFSet(structure=structure, prev_incar=incar,
                           prev_chgcar_path=chgcar_path,
                           reciprocal_density=reciprocal_density, **kwargs)


class MPSOCSet(MPStaticSet):

    def __init__(self, structure, saxis=(0, 0, 1), prev_incar=None,
                 prev_chgcar_path=None, reciprocal_density=100, **kwargs):
        """
        Init a MPSOCSet. Typically, you would use the classmethod
        from_prev_calc instead.

        Args:
            structure (Structure): the structure must have the 'magmom' site
                property and each magnetic moment value must have 3
                components. eg:- magmom = [[0,0,2], ...]
            saxis (tuple): magnetic moment orientation
            prev_incar (Incar): Incar file from previous run.
            prev_chgcar_path (str): Path to CHGCAR from previous run.
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
        self.prev_chgcar_path = prev_chgcar_path
        super(MPSOCSet, self).__init__(structure, prev_incar=prev_incar,
                                      reciprocal_density=reciprocal_density,
                                      **kwargs)

    @property
    def incar(self):
        incar = self.parent_vis.get_incar(self.structure)
        if self.prev_incar is not None:
            incar.update({k: v for k, v in self.prev_incar.items()
                         if k not in self.kwargs.get("user_incar_settings",
                                                     {})})

        # Overwrite necessary INCAR parameters from previous runs
        incar.update({"ISYM": -1, "LSORBIT": "T", "ICHARG": 11,
                      "SAXIS": list(self.saxis)})

        return incar

    def write_input(self, output_dir,
                    make_dir_if_not_present=True, include_cif=False):
        super(MPSOCSet, self).write_input(output_dir,
            make_dir_if_not_present=make_dir_if_not_present,
            include_cif=include_cif)
        if self.prev_chgcar_path:
            shutil.copy(self.prev_chgcar_path,
                        os.path.join(output_dir, "CHGCAR"))

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

        chgcar_path = None
        if copy_chgcar:
            chgcars = glob(os.path.join(prev_calc_dir, "CHGCAR*"))
            if chgcars:
                chgcar_path = sorted(chgcars)[-1]

        # multiply the reciprocal density if needed:
        if small_gap_multiply:
            gap = vasprun.eigenvalue_band_properties[0]
            if gap <= small_gap_multiply[0]:
                reciprocal_density = reciprocal_density * small_gap_multiply[1]

        return MPSOCSet(structure, prev_incar=incar,
                        prev_chgcar_path=chgcar_path,
                        reciprocal_density=reciprocal_density, **kwargs)


class MVLSlabSet(DictVaspInputSet):
    """
    Class for writing a set of slab vasp runs,
    including both slabs (along the c direction) and orient unit cells (bulk),
    to ensure the same KPOINTS, POTCAR and INCAR criterion.

    Args:
        user_incar_settings(dict): A dict specifying additional incar
            settings, default to None, ediff_per_atom=False
        k_product: default to 50, kpoint number * length for a & b directions,
            also for c direction in bulk calculations
        potcar_functional: default to PBE
        bulk (bool): Set to True for bulk calculation. Defaults to False.
        **kwargs:
            Other kwargs supported by :class:`DictVaspInputSet`.
    """
    def __init__(self, user_incar_settings=None, gpu=False, k_product=50,
                 potcar_functional='PBE', bulk=False, **kwargs):

        user_incar_settings = user_incar_settings or {}
        vis = MPVaspInputSet(ediff_per_atom=False).as_dict()
        DictVaspInputSet.__init__(self, "MVLSlabSet",
                                  vis["config_dict"],
                                  **kwargs)
        incar_settings_basic = {
            "EDIFF": 1e-6, "EDIFFG": -0.01, "ENCUT": 400,
            "ISMEAR": 0, "SIGMA": 0.05, "ISIF": 3}

        if bulk:
             self.incar_settings.update(incar_settings_basic)
        else:
            incar_settings_basic["ISIF"] = 2
            incar_settings_basic["AMIN"] = 0.01
            incar_settings_basic["AMIX"] = 0.2
            incar_settings_basic["BMIX"] = 0.001
            incar_settings_basic["NELMIN"] = 8
            self.incar_settings.update(incar_settings_basic)

        if gpu:
            # Sets KPAR to allow for the use of VASP-gpu
            if "KPAR" not in user_incar_settings.keys():
                # KPAR setting of 1 is the safest setting,
                # but for optimal performance, we normally
                # use the square root of the number of KPOINTS
                self.incar_settings["KPAR"] = 1
            if "NPAR" in self.incar_settings.keys():
                del self.incar_settings["NPAR"]


        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

        self.user_incar_settings = user_incar_settings
        self.k_product = k_product
        self.potcar_functional = potcar_functional
        self.bulk = bulk

    def get_kpoints(self, structure):
        """
        k_product, default to 50, is kpoint number * length for a & b directions,
            also for c direction in bulk calculations
        Automatic mesh & Gamma is the default setting.
        """

        # To get input sets, the input structure has to has the same number
        # of required parameters as a Structure object (ie. 4). Slab
        # attributes aren't going to affect the VASP inputs anyways so
        # converting the slab into a structure should not matter

        kpt = super(MVLSlabSet, self).get_kpoints(structure)
        kpt.comment = "Automatic mesh"
        kpt.style = 'Gamma'

        # use k_product to calculate kpoints, k_product = kpts[0][0] * a
        abc = structure.lattice.abc
        kpt_calc = [int(self.k_product/abc[0]+0.5),
                    int(self.k_product/abc[1]+0.5), 1]
        self.kpt_calc = kpt_calc
        # calculate kpts (c direction) for bulk. (for slab, set to 1)
        if self.bulk:
            kpt_calc[2] = int(self.k_product/abc[2]+0.5)

        kpt.kpts[0] = kpt_calc

        return kpt

    def get_incar(self, structure):

        # To get input sets, the input structure has to has the same number
        # of required parameters as a Structure object (ie. 4). Slab
        # attributes aren't going to affect the VASP inputs anyways so
        # converting the slab into a structure should not matter

        abc = structure.lattice.abc
        kpt_calc = [int(self.k_product/abc[0]+0.5),
                    int(self.k_product/abc[1]+0.5),
                    int(self.k_product/abc[1]+0.5)]


        kpts = kpt_calc

        if kpts[0]<5 and kpts[1]<5:
            if not self.bulk:
                self.incar_settings.update(
                    {"ISMEAR": 0})
            else:
                if kpts[2]<5:
                    self.incar_settings.update(
                        {"ISMEAR": 0})
        if self.user_incar_settings:
                self.incar_settings.update(self.user_incar_settings)

        incr = super(MVLSlabSet, self).get_incar(structure)

        return incr

    def as_dict(self):
        d = super(MVLSlabSet, self).as_dict()
        d.update({
            "potcar_functional": self.potcar_functional,
            "user_incar_settings": self.user_incar_settings
        })
        return d

    def from_dict(cls, d):
        return cls(
                   user_incar_settings=d.get("user_incar_settings", None),
                   potcar_functional=d.get("potcar_functional", None))

    def get_all_vasp_input(self, structure):
        """
        Returns all input files as a dict of {filename: vaspio object}

        Args:
            structure (Structure/IStructure): Structure to generate vasp
                input for.
        Returns:
            dict of {filename: file_as_string}, e.g., {'INCAR':'EDIFF=1e-4...'}
        """

        # To get input sets, the input structure has to has the same number
        # of required parameters as a Structure object (ie. 4). Slab
        # attributes aren't going to affect the VASP inputs anyways so
        # converting the slab into a structure should not matter

        data = {'INCAR': self.get_incar(structure),
                'KPOINTS': self.get_kpoints(structure),
                'POSCAR': self.get_poscar(structure),
                'POTCAR': self.get_potcar(structure)}
        return data


def get_vasprun_outcar(path):
    vruns = glob(os.path.join(path, "vasprun.xml*"))
    outcars = glob(os.path.join(path, "OUTCAR*"))

    if len(vruns) == 0 or len(outcars) == 0:
        raise ValueError(
            "Unable to get vasprun.xml/OUTCAR from prev calculation in %s" %
            path)
    vsfile_fullpath = os.path.join(path, "vasprun.xml")
    outcarfile_fullpath = os.path.join(path, "OUTCAR")
    vsfile = vsfile_fullpath if vsfile_fullpath in vruns else sorted(vruns)[-1]
    outcarfile = outcarfile_fullpath if outcarfile_fullpath in outcars else sorted(outcars)[-1]
    return Vasprun(vsfile, parse_dos=False, parse_eigen=None), Outcar(
        outcarfile)


def get_structure_from_prev_run(vasprun, outcar=None, sym_prec=0.1,
                                international_monoclinic=True):
    """
    Process structure from previous run.

    Args:
        vasp_run (Vasprun): Vasprun that contains the final structure
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
            m[structure[0].specie.symbol] = vals.pop(0)
            for site in structure:
                if site.specie.symbol not in m:
                    m[site.specie.symbol] = vals.pop(0)
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
