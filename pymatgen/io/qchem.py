# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals
import textwrap

"""
This module implements input and output processing from QChem.
"""

import copy
import re
import numpy as np
from string import Template

import six

from monty.io import zopen
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Molecule
from pymatgen.core.units import Energy, FloatWithUnit
from monty.json import MSONable
from pymatgen.util.coord_utils import get_angle

__author__ = "Xiaohui Qu"
__copyright__ = "Copyright 2013, The Electrolyte Genome Project"
__version__ = "0.1"
__maintainer__ = "Xiaohui Qu"
__email__ = "xhqu1981@gmail.com"
__date__ = "11/4/13"


class QcTask(MSONable):
    """
    An object representing a QChem input file.

    Args:
        molecule: The input molecule. If it is None of string "read",
            QChem will read geometry from checkpoint file. If it is a
            Molecule object, QcInput will convert it into Cartesian
            coordinates. Valid values: pymatgen Molecule object, "read", None
            Defaults to None.
        charge (int): Charge of the molecule. If None, charge on molecule is
            used. Defaults to None.
        spin_multiplicity (int): Spin multiplicity of molecule. Defaults to
            None, which means that the spin multiplicity is set to 1 if the
            molecule has no unpaired electrons and to 2 if there are
            unpaired electrons.
        jobtype (str): The type the QChem job. "SP" for Single Point Energy,
            "opt" for geometry optimization, "freq" for
            vibrational frequency.
        title (str): Comments for the job. Defaults to None. Which means the
            $comment section will be discarded.
        exchange (str): The exchange methods of the theory. Examples including:
            "B" (in pure BLYP), "PW91", "PBE", "TPSS".
            Defaults to "HF".
            This parameter can also be common names of hybrid
            functionals, such as B3LYP, TPSSh, XYGJOS. In such cases,
            the correlation parameter should be left as None.
        correlation (str): The correlation level of the theory. Example
            including: "MP2", "RI-MP2", "CCSD(T)", "LYP", "PBE", "TPSS"
            Defaults to None.
        basis_set (str/dict): The basis set.
            If it is a dict, each element can use different basis set.
        aux_basis_set (str/dict): Auxiliary basis set. For methods,
            like RI-MP2, XYG3, OXYJ-OS, auxiliary basis set is required.
            If it is a dict, each element can use different auxiliary
            basis set.
        ecp: Effective core potential (ECP) to be used.
            If it is a dict, each element can use different ECP.
        rem_params (dict): The parameters supposed to write in the $rem
            section. Dict of key/value pairs.
            Example: {"scf_algorithm": "diis_gdm", "scf_max_cycles": 100}
        optional_params (dict): The parameter for keywords other than $rem
            section. Dict of key/value pairs.
            Example: {"basis": {"Li": "cc-PVTZ", "B": "aug-cc-PVTZ",
            "F": "aug-cc-PVTZ"} "ecp": {"Cd": "srsc", "Br": "srlc"}}
        ghost_atom (list): List of ghost atoms indices. Indices start from 0.
            The ghost atom will be represented in of the form of @element_symmbol
    """

    optional_keywords_list = {"basis", "basis2", "ecp", "empirical_dispersion",
                              "external_charges", "force_field_params",
                              "intracule", "isotopes", "aux_basis",
                              "localized_diabatization", "multipole_field",
                              "nbo", "occupied", "swap_occupied_virtual", "opt",
                              "pcm", "pcm_solvent", "solvent", "plots", "qm_atoms", "svp",
                              "svpirf", "van_der_waals", "xc_functional",
                              "cdft", "efp_fragments", "efp_params", "alist"}
    alternative_keys = {"job_type": "jobtype",
                        "symmetry_ignore": "sym_ignore",
                        "scf_max_cycles": "max_scf_cycles"}
    alternative_values = {"optimization": "opt",
                          "frequency": "freq"}
    zmat_patt = re.compile("^(\w+)*([\s,]+(\w+)[\s,]+(\w+))*[\-\.\s,\w]*$")
    xyz_patt = re.compile("^(\w+)[\s,]+([\d\.eE\-]+)[\s,]+([\d\.eE\-]+)[\s,]+"
                          "([\d\.eE\-]+)[\-\.\s,\w.]*$")

    def __init__(self, molecule=None, charge=None, spin_multiplicity=None,
                 jobtype='SP', title=None, exchange="HF", correlation=None,
                 basis_set="6-31+G*", aux_basis_set=None, ecp=None,
                 rem_params=None, optional_params=None, ghost_atoms=None):
        self.mol = copy.deepcopy(molecule) if molecule else "read"
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        if isinstance(self.mol, six.string_types):
            self.mol = self.mol.lower()
            if self.mol != "read":
                raise ValueError('The only accept text value for mol is "read"')
        elif isinstance(self.mol, list):
            for m in self.mol:
                if not isinstance(m, Molecule):
                    raise ValueError("In case of type list, every element of mol must be a pymatgen Molecule")
            if self.charge is None or self.spin_multiplicity is None:
                raise ValueError("For fragments molecule section input, charge and spin_multiplicity "
                                 "must be specificed")
            total_charge = sum([m.charge for m in self.mol])
            total_unpaired_electron = sum([m.spin_multiplicity-1 for m in self.mol])
            if total_charge != self.charge:
                raise ValueError("The charge of the molecule doesn't equal to the sum of the fragment charges")
            if total_unpaired_electron % 2 != (self.spin_multiplicity - 1) % 2:
                raise ValueError("Spin multiplicity of molecule and fragments doesn't match")
        elif isinstance(self.mol, Molecule):
            self.charge = charge if charge is not None else self.mol.charge
            ghost_nelectrons = 0
            if ghost_atoms:
                for i in ghost_atoms:
                    site = self.mol.sites[i]
                    for sp, amt in site.species_and_occu.items():
                        ghost_nelectrons += sp.Z * amt
            nelectrons = self.mol.charge + self.mol.nelectrons - ghost_nelectrons - self.charge
            if spin_multiplicity is not None:
                self.spin_multiplicity = spin_multiplicity
                if (nelectrons + spin_multiplicity) % 2 != 1:
                    raise ValueError("Charge of {} and spin multiplicity of {} "
                                     "is not possible for this molecule"
                                     .format(self.charge, spin_multiplicity))
            else:
                self.spin_multiplicity = 1 if nelectrons % 2 == 0 else 2
        else:
            raise ValueError("The molecule must be a pymatgen Molecule "
                             "object or read/None or list of pymatgen Molecule")
        if (self.charge is None) != (self.spin_multiplicity is None):
            raise ValueError("spin multiplicity must be set together")
        if self.charge is not None and isinstance(self.mol, Molecule) and not ghost_atoms:
            self.mol.set_charge_and_spin(self.charge, self.spin_multiplicity)
        self.params = dict()
        if title is not None:
            self.params["comment"] = self._wrap_comment(title)
        if "rem" not in self.params:
            self.params["rem"] = dict()
        self.params["rem"]["exchange"] = exchange.lower()
        available_jobtypes = {"sp", "opt", "ts", "freq", "force", "rpath",
                              "nmr", "bsse", "eda", "pes_scan", "fsm", "aimd",
                              "pimc", "makeefp"}
        jt = jobtype.lower()
        if jt in self.alternative_values:
            jt = self.alternative_values[jt]
        if jt not in available_jobtypes:
            raise ValueError("Job type " + jobtype + " is not supported yet")
        self.params["rem"]["jobtype"] = jobtype.lower()
        if correlation is not None:
            self.params["rem"]["correlation"] = correlation.lower()
        if rem_params is not None:
            for k, v in rem_params.items():
                k = k.lower()
                if k in self.alternative_keys:
                    k = self.alternative_keys[k]
                if isinstance(v, six.string_types):
                    v = str(v).lower()
                    if v in self.alternative_values:
                        # noinspection PyTypeChecker
                        v = self.alternative_values[v]
                    self.params["rem"][k] = v
                elif isinstance(v, int) or isinstance(v, float):
                    self.params["rem"][k] = v
                else:
                    raise ValueError("The value in $rem can only be Integer "
                                     "or string")
        if optional_params:
            op_key = set([k.lower() for k in optional_params.keys()])
            if len(op_key - self.optional_keywords_list) > 0:
                invalid_keys = op_key - self.optional_keywords_list
                raise ValueError(','.join(['$' + k for k in invalid_keys]) +
                                 'is not a valid optional section')
            self.params.update(optional_params)

        self.set_basis_set(basis_set)

        if aux_basis_set is None:
            if self._aux_basis_required():
                if isinstance(self.params["rem"]["basis"], six.string_types):
                    if self.params["rem"]["basis"].startswith("6-31+g"):
                        self.set_auxiliary_basis_set("rimp2-aug-cc-pvdz")
                    elif self.params["rem"]["basis"].startswith("6-311+g"):
                        self.set_auxiliary_basis_set("rimp2-aug-cc-pvtz")
                if "aux_basis" not in self.params["rem"]:
                    raise ValueError("Auxiliary basis set is missing")
        else:
            self.set_auxiliary_basis_set(aux_basis_set)

        if ecp:
            self.set_ecp(ecp)
        self.ghost_atoms = ghost_atoms
        if self.ghost_atoms:
            if not isinstance(self.ghost_atoms, list):
                raise ValueError("ghost_atoms must be a list of integers")
            for atom in self.ghost_atoms:
                if not isinstance(atom, int):
                    raise ValueError("Each element of ghost atom list must an integer")

    def _aux_basis_required(self):
        if self.params["rem"]["exchange"] in ['xygjos', 'xyg3', 'lxygjos']:
            return True
        if 'correlation' in self.params["rem"]:
            if self.params["rem"]["correlation"].startswith("ri"):
                return True

    def set_basis_set(self, basis_set):
        if isinstance(basis_set, six.string_types):
            self.params["rem"]["basis"] = str(basis_set).lower()
            if basis_set.lower() not in ["gen", "mixed"]:
                self.params.pop("basis", None)
        elif isinstance(basis_set, dict):
            self.params["rem"]["basis"] = "gen"
            bs = dict()
            for element, basis in basis_set.items():
                bs[element.strip().capitalize()] = basis.lower()
            self.params["basis"] = bs
            if self.mol:
                mol_elements = set([site.species_string for site
                                    in self.mol.sites])
                basis_elements = set(self.params["basis"].keys())
                if len(mol_elements - basis_elements) > 0:
                    raise ValueError("The basis set for elements " +
                                     ", ".join(
                                         list(mol_elements - basis_elements)) +
                                     " is missing")
                if len(basis_elements - mol_elements) > 0:
                    raise ValueError("Basis set error: the molecule "
                                     "doesn't contain element " +
                                     ", ".join(basis_elements - mol_elements))
        elif isinstance(basis_set, list):
            self.params["rem"]["basis"] = "mixed"
            bs = [(a[0].capitalize(), a[1].lower()) for a in basis_set]
            self.params["basis"] = bs
            if len(self.mol) != len(basis_set):
                raise ValueError("Must specific a basis set for every atom")
            mol_elements = [site.species_string for site in self.mol.sites]
            basis_elements = [a[0] for a in bs]
            if mol_elements != basis_elements:
                raise ValueError("Elements in molecule and mixed basis set don't match")
        else:
            raise Exception('Can\'t handle type "{}"'.format(type(basis_set)))

    def set_partial_hessian_atoms(self, alist, phess=1):
        for a in alist:
            if not isinstance(a, int):
                raise ValueError("the parament alist must a list of atom indices")
        self.params["rem"]["n_sol"] = len(alist)
        if phess == 1:
            self.params["rem"]["phess"] = True
        else:
            self.params["rem"]["phess"] = phess
        self.params["rem"]["jobtype"] = "freq"
        self.params["alist"] = alist

    def set_basis2(self, basis2_basis_set):
        if isinstance(basis2_basis_set, six.string_types):
            self.params["rem"]["basis2"] = basis2_basis_set.lower()
            if basis2_basis_set.lower() not in ["basis2_gen", "basis2_mixed"]:
                self.params.pop("basis2", None)
        elif isinstance(basis2_basis_set, dict):
            self.params["rem"]["basis2"] = "basis2_gen"
            bs = dict()
            for element, basis in basis2_basis_set.items():
                bs[element.strip().capitalize()] = basis.lower()
            self.params["basis2"] = bs
            if self.mol:
                mol_elements = set([site.species_string for site
                                    in self.mol.sites])
                basis_elements = set(self.params["basis2"].keys())
                if len(mol_elements - basis_elements) > 0:
                    raise ValueError("The BASIS2 basis set for "
                                     "elements " +
                                     ", ".join(
                                         list(mol_elements - basis_elements)) +
                                     " is missing")
                if len(basis_elements - mol_elements) > 0:
                    raise ValueError("BASIS2 basis set error: the "
                                     "molecule doesn't contain element " +
                                     ", ".join(basis_elements - mol_elements))
        elif isinstance(basis2_basis_set, list):
            self.params["rem"]["basis2"] = "basis2_mixed"
            bs = [(a[0].capitalize(), a[1].lower()) for a in basis2_basis_set]
            self.params["basis2"] = bs
            if len(self.mol) != len(basis2_basis_set):
                raise ValueError("Must specific a BASIS2 basis set for every atom")
            mol_elements = [site.species_string for site in self.mol.sites]
            basis_elements = [a[0] for a in bs]
            if mol_elements != basis_elements:
                raise ValueError("Elements in molecule and mixed basis set don't match")
        else:
            raise Exception('Can\'t handle type "{}"'.format(type(basis2_basis_set)))

    def set_auxiliary_basis_set(self, aux_basis_set):
        if isinstance(aux_basis_set, six.string_types):
            self.params["rem"]["aux_basis"] = aux_basis_set.lower()
            if aux_basis_set.lower() not in ["gen", "mixed"]:
                self.params.pop("aux_basis", None)
        elif isinstance(aux_basis_set, dict):
            self.params["rem"]["aux_basis"] = "gen"
            bs = dict()
            for element, basis in aux_basis_set.items():
                bs[element.strip().capitalize()] = basis.lower()
            self.params["aux_basis"] = bs
            if self.mol:
                mol_elements = set([site.species_string for site
                                    in self.mol.sites])
                basis_elements = set(self.params["aux_basis"].keys())
                if len(mol_elements - basis_elements) > 0:
                    raise ValueError("The auxiliary basis set for "
                                     "elements " +
                                     ", ".join(
                                         list(mol_elements - basis_elements)) +
                                     " is missing")
                if len(basis_elements - mol_elements) > 0:
                    raise ValueError("Auxiliary asis set error: the "
                                     "molecule doesn't contain element " +
                                     ", ".join(basis_elements - mol_elements))
        elif isinstance(aux_basis_set, list):
            self.params["rem"]["aux_basis"] = "mixed"
            bs = [(a[0].capitalize(), a[1].lower()) for a in aux_basis_set]
            self.params["aux_basis"] = bs
            if len(self.mol) != len(aux_basis_set):
                raise ValueError("Must specific a auxiliary basis set for every atom")
            mol_elements = [site.species_string for site in self.mol.sites]
            basis_elements = [a[0] for a in bs]
            if mol_elements != basis_elements:
                raise ValueError("Elements in molecule and mixed basis set don't match")
        else:
            raise Exception('Can\'t handle type "{}"'.format(type(aux_basis_set)))

    def set_ecp(self, ecp):
        if isinstance(ecp, six.string_types):
            self.params["rem"]["ecp"] = ecp.lower()
        elif isinstance(ecp, dict):
            self.params["rem"]["ecp"] = "gen"
            potentials = dict()
            for element, p in ecp.items():
                potentials[element.strip().capitalize()] = p.lower()
            self.params["ecp"] = potentials
            if self.mol:
                mol_elements = set([site.species_string for site
                                    in self.mol.sites])
                ecp_elements = set(self.params["ecp"].keys())
                if len(ecp_elements - mol_elements) > 0:
                    raise ValueError("ECP error: the molecule "
                                     "doesn't contain element " +
                                     ", ".join(ecp_elements - mol_elements))

    @property
    def molecule(self):
        return self.mol

    def set_memory(self, total=None, static=None):
        """
        Set the maxium allowed memory.

        Args:
            total: The total memory. Integer. Unit: MBytes. If set to None,
                this parameter will be neglected.
            static: The static memory. Integer. Unit MBytes. If set to None,
                this parameterwill be neglected.
        """
        if total:
            self.params["rem"]["mem_total"] = total
        if static:
            self.params["rem"]["mem_static"] = static

    def set_max_num_of_scratch_files(self, num=16):
        """
        In QChem, the size of a single scratch is limited 2GB. By default,
        the max number of scratich is 16, which is cooresponding to 32GB
        scratch space. If you want to use more scratch disk space, you need
        to increase the number of scratch files:

        Args:
            num: The max number of the scratch files. (Integer)
        """
        self.params["rem"]["max_sub_file_num"] = num

    def set_scf_algorithm_and_iterations(self, algorithm="diis",
                                         iterations=50):
        """
        Set algorithm used for converging SCF and max number of SCF iterations.

        Args:
            algorithm: The algorithm used for converging SCF. (str)
            iterations: The max number of SCF iterations. (Integer)
        """
        available_algorithms = {"diis", "dm", "diis_dm", "diis_gdm", "gdm",
                                "rca", "rca_diis", "roothaan"}
        if algorithm.lower() not in available_algorithms:
            raise ValueError("Algorithm " + algorithm +
                             " is not available in QChem")
        self.params["rem"]["scf_algorithm"] = algorithm.lower()
        self.params["rem"]["max_scf_cycles"] = iterations

    def set_scf_convergence_threshold(self, exponent=8):
        """
        SCF is considered converged when the wavefunction error is less than
        10**(-exponent).
        In QChem, the default values are:
        5	For single point energy calculations.
        7	For geometry optimizations and vibrational analysis.
        8	For SSG calculations

        Args:
            exponent: The exponent of the threshold. (Integer)
        """
        self.params["rem"]["scf_convergence"] = exponent

    def set_integral_threshold(self, thresh=12):
        """
        Cutoff for neglect of two electron integrals. 10−THRESH (THRESH <= 14).
        In QChem, the default values are:
        8	For single point energies.
        10	For optimizations and frequency calculations.
        14	For coupled-cluster calculations.

        Args:
            thresh: The exponent of the threshold. (Integer)
        """
        self.params["rem"]["thresh"] = thresh

    def set_dft_grid(self, radical_points=128, angular_points=302,
                     grid_type="Lebedev"):
        """
        Set the grid for DFT numerical integrations.

        Args:
            radical_points: Radical points. (Integer)
            angular_points: Angular points. (Integer)
            grid_type: The type of of the grid. There are two standard grids:
                SG-1 and SG-0. The other two supported grids are "Lebedev" and
                "Gauss-Legendre"
        """
        available_lebedev_angular_points = {6, 18, 26, 38, 50, 74, 86, 110, 146,
                                            170, 194, 230, 266, 302, 350, 434,
                                            590, 770, 974, 1202, 1454, 1730,
                                            2030, 2354, 2702, 3074, 3470, 3890,
                                            4334, 4802, 5294}
        if grid_type.lower() == "sg-0":
            self.params["rem"]["xc_grid"] = 0
        elif grid_type.lower() == "sg-1":
            self.params["rem"]["xc_grid"] = 1
        elif grid_type.lower() == "lebedev":
            if angular_points not in available_lebedev_angular_points:
                raise ValueError(str(angular_points) + " is not a valid "
                                 "Lebedev angular points number")
            self.params["rem"]["xc_grid"] = "{rp:06d}{ap:06d}".format(
                rp=radical_points, ap=angular_points)
        elif grid_type.lower() == "gauss-legendre":
            self.params["rem"]["xc_grid"] = "-{rp:06d}{ap:06d}".format(
                rp=radical_points, ap=angular_points)
        else:
            raise ValueError("Grid type " + grid_type + " is not supported "
                                                        "currently")

    def set_scf_initial_guess(self, guess="SAD"):
        """
        Set initial guess method to be used for SCF

        Args:
            guess: The initial guess method. (str)
        """
        availabel_guesses = {"core", "sad", "gwh", "read", "fragmo"}
        if guess.lower() not in availabel_guesses:
            raise ValueError("The guess method " + guess + " is not supported "
                                                           "yet")
        self.params["rem"]["scf_guess"] = guess.lower()

    def set_geom_max_iterations(self, iterations):
        """
        Set the max iterations of geometry optimization.

        Args:
            iterations: the maximum iterations of geometry optimization.
            (Integer)
        """
        self.params["rem"]["geom_opt_max_cycles"] = iterations

    def set_geom_opt_coords_type(self, coords_type="internal_switch"):
        """
        Set the coordinates system used in geometry optimization.
        "cartesian"       --- always cartesian coordinates.
        "internal"        --- always internal coordinates.
        "internal-switch" --- try internal coordinates first, if fails, switch
        to cartesian coordinates.
        "z-matrix"        --- always z-matrix coordinates.
        "z-matrix-switch" --- try z-matrix first, if fails, switch to
        cartesian coordinates.

        Args:
            coords_type: The type of the coordinates. (str)
        """
        coords_map = {"cartesian": 0, "internal": 1, "internal-switch": -1,
                      "z-matrix": 2, "z-matrix-switch": -2}
        if coords_type.lower() not in set(coords_map.keys()):
            raise ValueError("Coodinate system " + coords_type + " is not "
                             "supported yet")
        else:
            self.params["rem"]["geom_opt_coords"] = \
                coords_map[coords_type.lower()]

    def scale_geom_opt_threshold(self, gradient=0.1, displacement=0.1,
                                 energy=0.1):
        """
        Adjust the convergence criteria of geometry optimization.

        Args:
            gradient: the scale factor for gradient criteria. If less than
                1.0, you are tightening the threshold. The base value is
                300 × 10E−6
            displacement: the scale factor for atomic displacement. If less
                then 1.0, you are tightening the threshold. The base value is
                1200 × 10E−6
            energy: the scale factor for energy change between successive
                iterations. If less than 1.0, you are tightening the
                threshold. The base value is 100 × 10E−8.
        """
        if gradient < 1.0/(300-1) or displacement < 1.0/(1200-1) or \
                energy < 1.0/(100-1):
            raise ValueError("The geometry optimization convergence criteria "
                             "is too tight")
        self.params["rem"]["geom_opt_tol_gradient"] = int(gradient * 300)
        self.params["rem"]["geom_opt_tol_displacement"] = int(displacement *
                                                              1200)
        self.params["rem"]["geom_opt_tol_energy"] = int(energy * 100)

    def set_geom_opt_use_gdiis(self, subspace_size=None):
        """
        Use GDIIS algorithm in geometry optimization.

        Args:
            subspace_size: The size of the DIIS subsapce. None for default
                value. The default value is min(NDEG, NATOMS, 4) NDEG = number
                of moleculardegrees of freedom.
        """
        subspace_size = subspace_size if subspace_size is not None else -1
        self.params["rem"]["geom_opt_max_diis"] = subspace_size

    def disable_symmetry(self):
        """
        Turn the symmetry off.
        """
        self.params["rem"]["sym_ignore"] = True
        self.params["rem"]["symmetry"] = False

    def use_cosmo(self, dielectric_constant=78.4):
        """
        Set the solvent model to COSMO.

        Args:
            dielectric_constant: the dielectric constant for the solvent.
        """
        self.params["rem"]["solvent_method"] = "cosmo"
        self.params["rem"]["solvent_dielectric"] = dielectric_constant

    def use_pcm(self, pcm_params=None, solvent_key="solvent", solvent_params=None,
                radii_force_field=None):
        """
        Set the solvent model to PCM. Default parameters are trying to comply to
        gaussian default value

        Args:
            pcm_params (dict): The parameters of "$pcm" section.
            solvent_key (str): for versions < 4.2 the section name is "pcm_solvent"
            solvent_params (dict): The parameters of solvent_key section
            radii_force_field (str): The force fied used to set the solute
                radii. Default to UFF.
        """
        self.params["pcm"] = dict()
        self.params[solvent_key] = dict()
        default_pcm_params = {"Theory": "SSVPE",
                              "vdwScale": 1.1,
                              "Radii": "UFF"}
        if not solvent_params:
            solvent_params = {"Dielectric": 78.3553}
        if pcm_params:
            for k, v in pcm_params.items():
                self.params["pcm"][k.lower()] = v.lower() \
                    if isinstance(v, six.string_types) else v

        for k, v in default_pcm_params.items():
            if k.lower() not in self.params["pcm"].keys():
                self.params["pcm"][k.lower()] = v.lower() \
                    if isinstance(v, six.string_types) else v
        for k, v in solvent_params.items():
            self.params[solvent_key][k.lower()] = v.lower() \
                if isinstance(v, six.string_types) else copy.deepcopy(v)
        self.params["rem"]["solvent_method"] = "pcm"
        if radii_force_field:
            self.params["pcm"]["radii"] = "bondi"
            self.params["rem"]["force_fied"] = radii_force_field.lower()

    def __str__(self):
        sections = ["comment", "molecule", "rem"] + \
            sorted(list(self.optional_keywords_list))
        lines = []
        for sec in sections:
            if sec in self.params or sec == "molecule":
                foramt_sec = self.__getattribute__("_format_" + sec)
                lines.append("$" + sec)
                lines.extend(foramt_sec())
                lines.append("$end")
                lines.append('\n')
        return '\n'.join(lines)

    @classmethod
    def _wrap_comment(cls, comment):
        ml_section_start = comment.find('<')
        if ml_section_start >= 0:
            title_section = comment[0:ml_section_start]
            ml_section = comment[ml_section_start:]
        else:
            title_section = comment
            ml_section = ''
        wrapped_title_lines = textwrap.wrap(title_section.strip(), width=70, initial_indent=' ')
        wrapped_ml_lines = []
        for l in ml_section.splitlines():
            if len(l) > 70:
                wrapped_ml_lines.extend(textwrap.wrap(l.strip(), width=70, initial_indent=' '))
            else:
                wrapped_ml_lines.append(l)
        return '\n'.join(wrapped_title_lines + wrapped_ml_lines)

    def _format_comment(self):
        return self._wrap_comment(self.params["comment"]).splitlines()

    def _format_alist(self):
        return [" {}".format(x) for x in self.params["alist"]]

    def _format_molecule(self):
        lines = []

        def inner_format_mol(m2, index_base):
            mol_lines = []
            for i, site in enumerate(m2.sites):
                ghost = "@" if self.ghost_atoms \
                    and i + index_base in self.ghost_atoms else ""
                atom = "{ghost:s}{element:s}".format(ghost=ghost,
                                                     element=site.species_string)
                mol_lines.append(" {atom:<4} {x:>17.8f} {y:>17.8f} "
                                 "{z:>17.8f}".format(atom=atom, x=site.x,
                                                     y=site.y, z=site.z))
            return mol_lines

        if self.charge is not None:
            lines.append(" {charge:d}  {multi:d}".format(charge=self
                         .charge, multi=self.spin_multiplicity))
        if isinstance(self.mol, six.string_types) and self.mol == "read":
            lines.append(" read")
        elif isinstance(self.mol, list):
            starting_index = 0
            for m in self.mol:
                lines.append("--")
                lines.append(" {charge:d}  {multi:d}".format(
                    charge=m.charge, multi=m.spin_multiplicity))
                lines.extend(inner_format_mol(m, starting_index))
                starting_index += len(m)
        else:
            lines.extend(inner_format_mol(self.mol, 0))
        return lines

    def _format_rem(self):
        rem_format_template = Template("  {name:>$name_width} = "
                                       "{value}")
        name_width = 0
        for name, value in self.params["rem"].items():
            if len(name) > name_width:
                name_width = len(name)
        rem = rem_format_template.substitute(name_width=name_width)
        lines = []
        all_keys = set(self.params["rem"].keys())
        priority_keys = ["jobtype", "exchange", "basis"]
        additional_keys = all_keys - set(priority_keys)
        ordered_keys = priority_keys + sorted(list(additional_keys))
        for name in ordered_keys:
            value = self.params["rem"][name]
            lines.append(rem.format(name=name, value=value))
        return lines

    def _format_basis(self):
        lines = []
        if isinstance(self.params["basis"], dict):
            for element in sorted(self.params["basis"].keys()):
                basis = self.params["basis"][element]
                lines.append(" " + element)
                lines.append(" " + basis)
                lines.append(" ****")
        elif isinstance(self.params["basis"], list):
            for i, (element, bs) in enumerate(self.params["basis"]):
                lines.append(" {element:2s} {number:3d}".format(element=element, number=i+1))
                lines.append(" {}".format(bs))
                lines.append(" ****")
        return lines

    def _format_aux_basis(self):
        lines = []
        if isinstance(self.params["aux_basis"], dict):
            for element in sorted(self.params["aux_basis"].keys()):
                basis = self.params["aux_basis"][element]
                lines.append(" " + element)
                lines.append(" " + basis)
                lines.append(" ****")
        else:
            for i, (element, bs) in enumerate(self.params["aux_basis"]):
                lines.append(" {element:2s} {number:3d}".format(element=element, number=i+1))
                lines.append(" {}".format(bs))
                lines.append(" ****")
        return lines

    def _format_basis2(self):
        lines = []
        if isinstance(self.params["basis2"], dict):
            for element in sorted(self.params["basis2"].keys()):
                basis = self.params["basis2"][element]
                lines.append(" " + element)
                lines.append(" " + basis)
                lines.append(" ****")
        else:
            for i, (element, bs) in enumerate(self.params["basis2"]):
                lines.append(" {element:2s} {number:3d}".format(element=element, number=i+1))
                lines.append(" {}".format(bs))
                lines.append(" ****")
        return lines

    def _format_ecp(self):
        lines = []
        for element in sorted(self.params["ecp"].keys()):
            ecp = self.params["ecp"][element]
            lines.append(" " + element)
            lines.append(" " + ecp)
            lines.append(" ****")
        return lines

    def _format_pcm(self):
        pcm_format_template = Template("  {name:>$name_width}   "
                                       "{value}")
        name_width = 0
        for name in self.params["pcm"].keys():
            if len(name) > name_width:
                name_width = len(name)
        rem = pcm_format_template.substitute(name_width=name_width)
        lines = []
        for name in sorted(self.params["pcm"].keys()):
            value = self.params["pcm"][name]
            lines.append(rem.format(name=name, value=value))
        return lines

    def _format_pcm_solvent(self, key="pcm_solvent"):
        pp_format_template = Template("  {name:>$name_width}   "
                                      "{value}")
        name_width = 0
        for name in self.params[key].keys():
            if len(name) > name_width:
                name_width = len(name)
        rem = pp_format_template.substitute(name_width=name_width)
        lines = []
        all_keys = set(self.params[key].keys())
        priority_keys = []
        for k in ["dielectric", "nonels", "nsolventatoms", "solventatom"]:
            if k in all_keys:
                priority_keys.append(k)
        additional_keys = all_keys - set(priority_keys)
        ordered_keys = priority_keys + sorted(list(additional_keys))
        for name in ordered_keys:
            value = self.params[key][name]
            if name == "solventatom":
                for v in copy.deepcopy(value):
                    value = "{:<4d} {:<4d} {:<4d} {:4.2f}".format(*v)
                    lines.append(rem.format(name=name, value=value))
                continue
            lines.append(rem.format(name=name, value=value))
        return lines

    def _format_solvent(self):
        return self._format_pcm_solvent(key="solvent")

    def _format_opt(self):
        # lines is a list of all opt keywords
        lines = []
        # only constraints added at this point
        constraint_lines = ['CONSTRAINT']
        for index in range(len(self.params['opt'])):
            vals = self.params['opt'][index]
            if vals[0] in ['outp', 'tors', 'linc', 'linp']:
                constraint_lines.append("{vals[0]} {vals[1]} {vals[2]} {vals[3]} {vals[4]} {vals[5]}".format(vals=vals))
            elif vals[0] == 'stre':
                constraint_lines.append("{vals[0]} {vals[1]} {vals[2]}".format(vals=vals))
            elif vals[0] == 'bend':
                constraint_lines.append("{vals[0]} {vals[1]} {vals[2]} {vals[3]} {vals[4]}".format(vals=vals))
        constraint_lines.append('ENDCONSTRAINT')
        lines.extend(constraint_lines)
        # another opt keyword can be added by extending lines
        return lines

    def as_dict(self):
        if isinstance(self.mol, six.string_types):
            mol_dict = self.mol
        elif isinstance(self.mol, Molecule):
            mol_dict = self.mol.as_dict()
        elif isinstance(self.mol, list):
            mol_dict = [m.as_dict() for m in self.mol]
        else:
            raise ValueError('Unknow molecule type "{}"'.format(type(self.mol)))
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "molecule": mol_dict,
             "charge": self.charge,
             "spin_multiplicity": self.spin_multiplicity,
             "params": self.params}
        if self.ghost_atoms:
            d["ghost_atoms"] = self.ghost_atoms
        return d

    @classmethod
    def from_dict(cls, d):
        if d["molecule"] == "read":
            mol = "read"
        elif isinstance(d["molecule"], dict):
            mol = Molecule.from_dict(d["molecule"])
        elif isinstance(d["molecule"], list):
            mol = [Molecule.from_dict(m) for m in d["molecule"]]
        else:
            raise ValueError('Unknow molecule type "{}"'.format(type(d["molecule"])))
        jobtype = d["params"]["rem"]["jobtype"]
        title = d["params"].get("comment", None)
        exchange = d["params"]["rem"]["exchange"]
        correlation = d["params"]["rem"].get("correlation", None)
        basis_set = d["params"]["rem"]["basis"]
        aux_basis_set = d["params"]["rem"].get("aux_basis", None)
        ecp = d["params"]["rem"].get("ecp", None)
        ghost_atoms = d.get("ghost_atoms", None)
        optional_params = None
        op_keys = set(d["params"].keys()) - {"comment", "rem"}
        if len(op_keys) > 0:
            optional_params = dict()
            for k in op_keys:
                optional_params[k] = d["params"][k]
        return QcTask(molecule=mol, charge=d["charge"],
                      spin_multiplicity=d["spin_multiplicity"],
                      jobtype=jobtype, title=title,
                      exchange=exchange, correlation=correlation,
                      basis_set=basis_set, aux_basis_set=aux_basis_set,
                      ecp=ecp, rem_params=d["params"]["rem"],
                      optional_params=optional_params,
                      ghost_atoms=ghost_atoms)

    def write_file(self, filename):
        with zopen(filename, "wt") as f:
            f.write(self.__str__())

    @classmethod
    def from_file(cls, filename):
        with zopen(filename, "rt") as f:
            return cls.from_string(f.read())

    @classmethod
    def from_string(cls, contents):
        """
        Creates QcInput from a string.

        Args:
            contents: String representing a QChem input file.

        Returns:
            QcInput object
        """
        mol = None
        charge = None
        spin_multiplicity = None
        params = dict()
        lines = contents.split('\n')
        parse_section = False
        section_name = None
        section_text = []
        ghost_atoms = None
        for line_num, line in enumerate(lines):
            l = line.strip().lower()

            if len(l) == 0:
                continue
            if (not parse_section) and (l == "$end" or not l.startswith("$")):
                raise ValueError("Format error, parsing failed")
            if parse_section and l != "$end":
                section_text.append(line)
            if l.startswith("$") and not parse_section:
                parse_section = True
                section_name = l[1:]
                available_sections = ["comment", "molecule", "rem"] + \
                    sorted(list(cls.optional_keywords_list))
                if section_name not in available_sections:
                    raise ValueError("Unrecognized keyword " + line.strip() +
                                     " at line " + str(line_num))
                if section_name in params:
                    raise ValueError("duplicated keyword " + line.strip() +
                                     "at line " + str(line_num))
            if parse_section and l == "$end":
                func_name = "_parse_" + section_name
                if func_name not in QcTask.__dict__:
                    raise Exception(func_name + " is not implemented yet, "
                                    "please implement it")
                parse_func = QcTask.__dict__[func_name].__get__(None, QcTask)
                if section_name == "molecule":
                    mol, charge, spin_multiplicity, ghost_atoms = parse_func(section_text)
                else:
                    d = parse_func(section_text)
                    params[section_name] = d
                parse_section = False
                section_name = None
                section_text = []
        if parse_section:
            raise ValueError("Format error. " + section_name + " is not "
                             "terminated")
        jobtype = params["rem"]["jobtype"]
        title = params.get("comment", None)
        exchange = params["rem"].get("exchange", "hf")
        correlation = params["rem"].get("correlation", None)
        basis_set = params["rem"]["basis"]
        aux_basis_set = params["rem"].get("aux_basis", None)
        ecp = params["rem"].get("ecp", None)
        optional_params = None
        op_keys = set(params.keys()) - {"comment", "rem"}
        if len(op_keys) > 0:
            optional_params = dict()
            for k in op_keys:
                optional_params[k] = params[k]
        return QcTask(molecule=mol, charge=charge,
                      spin_multiplicity=spin_multiplicity,
                      jobtype=jobtype, title=title,
                      exchange=exchange, correlation=correlation,
                      basis_set=basis_set, aux_basis_set=aux_basis_set,
                      ecp=ecp, rem_params=params["rem"],
                      optional_params=optional_params,
                      ghost_atoms=ghost_atoms)

    @classmethod
    def _parse_comment(cls, contents):
        return '\n'.join(contents).strip()

    @classmethod
    def _parse_coords(cls, coord_lines):
        """
        Helper method to parse coordinates. Copied from GaussianInput class.
        """
        paras = {}
        var_pattern = re.compile("^([A-Za-z]+\S*)[\s=,]+([\d\-\.]+)$")
        for l in coord_lines:
            m = var_pattern.match(l.strip())
            if m:
                paras[m.group(1)] = float(m.group(2))

        species = []
        coords = []
        # Stores whether a Zmatrix format is detected. Once a zmatrix format
        # is detected, it is assumed for the remaining of the parsing.
        zmode = False
        for l in coord_lines:
            l = l.strip()
            if not l:
                break
            if (not zmode) and cls.xyz_patt.match(l):
                m = cls.xyz_patt.match(l)
                species.append(m.group(1))
                toks = re.split("[,\s]+", l.strip())
                if len(toks) > 4:
                    coords.append(list(map(float, toks[2:5])))
                else:
                    coords.append(list(map(float, toks[1:4])))
            elif cls.zmat_patt.match(l):
                zmode = True
                toks = re.split("[,\s]+", l.strip())
                species.append(toks[0])
                toks.pop(0)
                if len(toks) == 0:
                    coords.append(np.array([0.0, 0.0, 0.0]))
                else:
                    nn = []
                    parameters = []
                    while len(toks) > 1:
                        ind = toks.pop(0)
                        data = toks.pop(0)
                        try:
                            nn.append(int(ind))
                        except ValueError:
                            nn.append(species.index(ind) + 1)
                        try:
                            val = float(data)
                            parameters.append(val)
                        except ValueError:
                            if data.startswith("-"):
                                parameters.append(-paras[data[1:]])
                            else:
                                parameters.append(paras[data])
                    if len(nn) == 1:
                        coords.append(np.array(
                            [0.0, 0.0, float(parameters[0])]))
                    elif len(nn) == 2:
                        coords1 = coords[nn[0] - 1]
                        coords2 = coords[nn[1] - 1]
                        bl = parameters[0]
                        angle = parameters[1]
                        axis = [0, 1, 0]
                        op = SymmOp.from_origin_axis_angle(coords1, axis,
                                                           angle, False)
                        coord = op.operate(coords2)
                        vec = coord - coords1
                        coord = vec * bl / np.linalg.norm(vec) + coords1
                        coords.append(coord)
                    elif len(nn) == 3:
                        coords1 = coords[nn[0] - 1]
                        coords2 = coords[nn[1] - 1]
                        coords3 = coords[nn[2] - 1]
                        bl = parameters[0]
                        angle = parameters[1]
                        dih = parameters[2]
                        v1 = coords3 - coords2
                        v2 = coords1 - coords2
                        axis = np.cross(v1, v2)
                        op = SymmOp.from_origin_axis_angle(
                            coords1, axis, angle, False)
                        coord = op.operate(coords2)
                        v1 = coord - coords1
                        v2 = coords1 - coords2
                        v3 = np.cross(v1, v2)
                        adj = get_angle(v3, axis)
                        axis = coords1 - coords2
                        op = SymmOp.from_origin_axis_angle(
                            coords1, axis, dih - adj, False)
                        coord = op.operate(coord)
                        vec = coord - coords1
                        coord = vec * bl / np.linalg.norm(vec) + coords1
                        coords.append(coord)

        def parse_species(sp_str):
            """
            The species specification can take many forms. E.g.,
            simple integers representing atomic numbers ("8"),
            actual species string ("C") or a labelled species ("C1").
            Sometimes, the species string is also not properly capitalized,
            e.g, ("c1"). This method should take care of these known formats.
            """
            try:
                return int(sp_str)
            except ValueError:
                sp = re.sub("\d", "", sp_str)
                return sp.capitalize()

        species = list(map(parse_species, species))

        return Molecule(species, coords)

    @classmethod
    def _parse_molecule(cls, contents):
        def parse_ghost_indices(coord_text_lines):
            no_ghost_text = [l.replace("@", "") for l in coord_text_lines]
            ghosts = []
            for index, l in enumerate(coord_text_lines):
                l = l.strip()
                if not l:
                    break
                if "@" in l:
                    ghosts.append(index)
            return ghosts, no_ghost_text

        text = copy.deepcopy(contents[:2])
        charge_multi_pattern = re.compile('\s*(?P<charge>'
                                          '[-+]?\d+)\s+(?P<multi>\d+)')
        line = text.pop(0)
        m = charge_multi_pattern.match(line)
        if m:
            charge = int(m.group("charge"))
            spin_multiplicity = int(m.group("multi"))
            line = text.pop(0)
        else:
            charge = None
            spin_multiplicity = None
        if line.strip().lower() == "read":
            return "read", charge, spin_multiplicity, None
        elif charge is None or spin_multiplicity is None:
            raise ValueError("Charge or spin multiplicity is not found")
        else:
            if contents[1].strip()[0:2] == "--":
                chunks = "\n".join(contents[2:]).split("--\n")
                mol = []
                ghost_atoms = []
                starting_index = 0
                for chunk in chunks:
                    frag_contents = chunk.split("\n")
                    m = charge_multi_pattern.match(frag_contents[0])
                    if m:
                        fragment_charge = int(m.group("charge"))
                        fragment_spin_multiplicity = int(m.group("multi"))
                    else:
                        raise Exception("charge and spin multiplicity must be specified for each fragment")
                    gh, coord_lines = parse_ghost_indices(frag_contents[1:])
                    fragment = cls._parse_coords(coord_lines)
                    fragment.set_charge_and_spin(fragment_charge, fragment_spin_multiplicity)
                    mol.append(fragment)
                    ghost_atoms.extend([i+starting_index for i in gh])
                    starting_index += len(fragment)
            else:
                ghost_atoms, coord_lines = parse_ghost_indices(contents[1:])
                mol = cls._parse_coords(coord_lines)
                if len(ghost_atoms) == 0:
                    mol.set_charge_and_spin(charge, spin_multiplicity)
            ghost_atoms = ghost_atoms if len(ghost_atoms) > 0 else None
            return mol, charge, spin_multiplicity, ghost_atoms

    @classmethod
    def _parse_rem(cls, contents):
        d = dict()
        int_pattern = re.compile('^[-+]?\d+$')
        float_pattern = re.compile('^[-+]?\d+\.\d+([eE][-+]?\d+)?$')

        for line in contents:
            tokens = line.strip().replace("=", ' ').split()
            if len(tokens) < 2:
                raise ValueError("Can't parse $rem section, there should be "
                                 "at least two field: key and value!")
            k1, v = tokens[:2]
            k2 = k1.lower()
            if k2 in cls.alternative_keys:
                k2 = cls.alternative_keys[k2]
            if v in cls.alternative_values:
                v = cls.alternative_values
            if k2 == "xc_grid":
                d[k2] = v
            elif v == "True":
                d[k2] = True
            elif v == "False":
                d[k2] = False
            elif int_pattern.match(v):
                d[k2] = int(v)
            elif float_pattern.match(v):
                d[k2] = float(v)
            else:
                d[k2] = v.lower()
        return d

    @classmethod
    def _parse_aux_basis(cls, contents):
        if len(contents) % 3 != 0:
            raise ValueError("Auxiliary basis set section format error")
        chunks = zip(*[iter(contents)]*3)
        t = contents[0].split()
        if len(t) == 2 and int(t[1]) > 0:
            bs = []
            for i, ch in enumerate(chunks):
                element, number = ch[0].split()
                basis = ch[1]
                if int(number) != i+1:
                    raise ValueError("Atom order number doesn't match in $aux_basis section")
                bs.append((element.strip().capitalize(), basis.strip().lower()))
        else:
            bs = dict()
            for ch in chunks:
                element, basis = ch[:2]
                bs[element.strip().capitalize()] = basis.strip().lower()
        return bs

    @classmethod
    def _parse_basis2(cls, contents):
        if len(contents) % 3 != 0:
            raise ValueError("Auxiliary basis set section format error")
        chunks = zip(*[iter(contents)]*3)
        t = contents[0].split()
        if len(t) == 2 and int(t[1]) > 0:
            bs = []
            for i, ch in enumerate(chunks):
                element, number = ch[0].split()
                basis = ch[1]
                if int(number) != i+1:
                    raise ValueError("Atom order number doesn't match in $aux_basis section")
                bs.append((element.strip().capitalize(), basis.strip().lower()))
        else:
            bs = dict()
            for ch in chunks:
                element, basis = ch[:2]
                bs[element.strip().capitalize()] = basis.strip().lower()
        return bs

    @classmethod
    def _parse_basis(cls, contents):
        if len(contents) % 3 != 0:
            raise ValueError("Basis set section format error")
        chunks = zip(*[iter(contents)]*3)
        t = contents[0].split()
        if len(t) == 2 and int(t[1]) > 0:
            bs = []
            for i, ch in enumerate(chunks):
                element, number = ch[0].split()
                basis = ch[1]
                if int(number) != i+1:
                    raise ValueError("Atom order number doesn't match in $basis section")
                bs.append((element.strip().capitalize(), basis.strip().lower()))
        else:
            bs = dict()
            for ch in chunks:
                element, basis = ch[:2]
                bs[element.strip().capitalize()] = basis.strip().lower()
        return bs

    @classmethod
    def _parse_ecp(cls, contents):
        if len(contents) % 3 != 0:
            raise ValueError("ECP section format error")
        chunks = zip(*[iter(contents)]*3)
        d = dict()
        for ch in chunks:
            element, ecp = ch[:2]
            d[element.strip().capitalize()] = ecp.strip().lower()
        return d

    @classmethod
    def _parse_alist(cls, contents):
        atom_list = []
        for line in contents:
            atom_list.extend([int(x) for x in line.split()])
        return atom_list

    @classmethod
    def _parse_pcm(cls, contents):
        d = dict()
        int_pattern = re.compile('^[-+]?\d+$')
        float_pattern = re.compile('^[-+]?\d+\.\d+([eE][-+]?\d+)?$')

        for line in contents:
            tokens = line.strip().replace("=", ' ').split()
            if len(tokens) < 2:
                raise ValueError("Can't parse $pcm section, there should be "
                                 "at least two field: key and value!")
            k1, v = tokens[:2]
            k2 = k1.lower()
            if k2 in cls.alternative_keys:
                k2 = cls.alternative_keys[k2]
            if v in cls.alternative_values:
                v = cls.alternative_values
            if v == "True":
                d[k2] = True
            elif v == "False":
                d[k2] = False
            elif int_pattern.match(v):
                d[k2] = int(v)
            elif float_pattern.match(v):
                d[k2] = float(v)
            else:
                d[k2] = v.lower()
        return d

    @classmethod
    def _parse_pcm_solvent(cls, contents):
        d = dict()
        int_pattern = re.compile('^[-+]?\d+$')
        float_pattern = re.compile('^[-+]?\d+\.\d+([eE][-+]?\d+)?$')

        for line in contents:
            tokens = line.strip().replace("=", ' ').split()
            if len(tokens) < 2:
                raise ValueError("Can't parse $pcm_solvent section, "
                                 "there should be at least two field: "
                                 "key and value!")
            k1, v = tokens[:2]
            k2 = k1.lower()
            if k2 in cls.alternative_keys:
                k2 = cls.alternative_keys[k2]
            if v in cls.alternative_values:
                v = cls.alternative_values
            if k2 == "solventatom":
                v = [int(i) for i in tokens[1:4]]
                # noinspection PyTypeChecker
                v.append(float(tokens[4]))
                if k2 not in d:
                    d[k2] = [v]
                else:
                    d[k2].append(v)
            elif v == "True":
                d[k2] = True
            elif v == "False":
                d[k2] = False
            elif int_pattern.match(v):
                d[k2] = int(v)
            elif float_pattern.match(v):
                d[k2] = float(v)
            else:
                d[k2] = v.lower()
        return d

    @classmethod
    def _parse_solvent(cls, contents):
        return cls._parse_pcm_solvent(contents)

    @classmethod
    def _parse_opt(cls, contents):
        #only parses opt constraints
        opt_list = []
        constraints = False
        int_pattern = re.compile('^[-+]?\d+$')
        float_pattern = re.compile('^[-+]?\d+\.\d+([eE][-+]?\d+)?$')
        for line in contents:
            tokens = line.strip().split()
            if re.match('ENDCONSTRAINT', line):
                constraints = False
            elif constraints:
                vals = []
                for val in tokens:
                    if int_pattern.match(val):
                        vals.append(int(val))
                    elif float_pattern.match(val):
                        vals.append(float(val))
                    else:
                        vals.append(val)
                opt_list.append(vals)
            elif re.match('CONSTRAINT', line):
                constraints = True
        return opt_list


class QcInput(MSONable):
    """
    An object representing a multiple step QChem input file.

    Args:
        jobs: The QChem jobs (List of QcInput object)
    """

    def __init__(self, jobs):
        jobs = jobs if isinstance(jobs, list) else [jobs]
        for j in jobs:
            if not isinstance(j, QcTask):
                raise ValueError("jobs must be a list QcInput object")
            self.jobs = jobs

    def __str__(self):
        return "\n@@@\n\n\n".join([str(j) for j in self.jobs])

    def write_file(self, filename):
        with zopen(filename, "wt") as f:
            f.write(self.__str__())

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "jobs": [j.as_dict() for j in self.jobs]}

    @classmethod
    def from_dict(cls, d):
        jobs = [QcTask.from_dict(j) for j in d["jobs"]]
        return QcInput(jobs)

    @classmethod
    def from_string(cls, contents):
        qc_contents = contents.split("@@@")
        jobs = [QcTask.from_string(cont) for cont in qc_contents]
        return QcInput(jobs)

    @classmethod
    def from_file(cls, filename):
        with zopen(filename, "rt") as f:
            return cls.from_string(f.read())


class QcOutput(object):

    kcal_per_mol_2_eV = 4.3363E-2

    def __init__(self, filename):
        self.filename = filename
        split_pattern = "\n\nRunning Job \d+ of \d+ \S+|" \
                        "[*]{61}\nJob \d+ of \d+ \n[*]{61}|" \
                        "\n.*time.*\nRunning Job \d+ of \d+ \S+"
        try:
            with zopen(filename, "rt") as f:
                data = f.read()
        except UnicodeDecodeError:
            with zopen(filename, "rb") as f:
                data = f.read().decode("latin-1")
        try:
            chunks = re.split(split_pattern, data)
            # noinspection PyTypeChecker
            self.data = list(map(self._parse_job, chunks))
        except UnicodeDecodeError:
            data = data.decode("latin-1")
            chunks = re.split(split_pattern, data)
            # noinspection PyTypeChecker
            self.data = list(map(self._parse_job, chunks))

    @property
    def final_energy(self):
        return self.data[-1]["energies"][-1][-1]

    @property
    def final_structure(self):
        return self.data[-1]["molecules"][-1]

    @classmethod
    def _expected_successful_pattern(cls, qctask):
        text = ["Convergence criterion met"]
        if "correlation" in qctask.params["rem"]:
            if "ccsd" in qctask.params["rem"]["correlation"]\
                    or "qcisd" in qctask.params["rem"]["correlation"]:
                text.append('CC.*converged')
        if qctask.params["rem"]["jobtype"] == "opt"\
                or qctask.params["rem"]["jobtype"] == "ts":
            text.append("OPTIMIZATION CONVERGED")
        if qctask.params["rem"]["jobtype"] == "freq":
            text.append("VIBRATIONAL ANALYSIS")
        if qctask.params["rem"]["jobtype"] == "gradient":
            text.append("Gradient of SCF Energy")
        return text

    @classmethod
    def _parse_job(cls, output):
        scf_energy_pattern = re.compile("Total energy in the final basis set ="
                                        "\s+(?P<energy>-\d+\.\d+)")
        corr_energy_pattern = re.compile("(?P<name>[A-Z\-\(\)0-9]+)\s+"
                                         "([tT]otal\s+)?[eE]nergy\s+=\s+"
                                         "(?P<energy>-\d+\.\d+)")
        coord_pattern = re.compile("\s*\d+\s+(?P<element>[A-Z][a-zH]*)\s+"
                                   "(?P<x>\-?\d+\.\d+)\s+"
                                   "(?P<y>\-?\d+\.\d+)\s+"
                                   "(?P<z>\-?\d+\.\d+)")
        num_ele_pattern = re.compile("There are\s+(?P<alpha>\d+)\s+alpha "
                                     "and\s+(?P<beta>\d+)\s+beta electrons")
        total_charge_pattern = re.compile("Sum of atomic charges ="
                                          "\s+(?P<charge>\-?\d+\.\d+)")
        scf_iter_pattern = re.compile("\d+\s*(?P<energy>\-\d+\.\d+)\s+"
                                      "(?P<diis_error>\d+\.\d+E[-+]\d+)")
        zpe_pattern = re.compile("Zero point vibrational energy:"
                                 "\s+(?P<zpe>\d+\.\d+)\s+kcal/mol")
        thermal_corr_pattern = re.compile("(?P<name>\S.*\S):\s+"
                                          "(?P<correction>\d+\.\d+)\s+"
                                          "k?cal/mol")
        detailed_charge_pattern = re.compile("(Ground-State )?(?P<method>\w+)( Net)?"
                                             " Atomic Charges")
        nbo_charge_pattern = re.compile("(?P<element>[A-Z][a-z]{0,2})\s*(?P<no>\d+)\s+(?P<charge>\-?\d\.\d+)"
                                        "\s+(?P<core>\-?\d+\.\d+)\s+(?P<valence>\-?\d+\.\d+)"
                                        "\s+(?P<rydberg>\-?\d+\.\d+)\s+(?P<total>\-?\d+\.\d+)"
                                        "(\s+(?P<spin>\-?\d\.\d+))?")
        nbo_wavefunction_type_pattern = re.compile("This is an? (?P<type>\w+\-\w+) NBO calculation")
        bsse_pattern = re.compile("DE, kJ/mol\s+(?P<raw_be>\-?\d+\.?\d+([eE]\d+)?)\s+"
                                  "(?P<corrected_be>\-?\d+\.?\d+([eE]\d+)?)")
        float_pattern = re.compile("\-?\d+\.?\d+([eE]\d+)?$")

        error_defs = (
            (re.compile("Convergence failure"), "Bad SCF convergence"),
            (re.compile("Coordinates do not transform within specified "
                        "threshold"), "autoz error"),
            (re.compile("MAXIMUM OPTIMIZATION CYCLES REACHED"),
                "Geometry optimization failed"),
            (re.compile("\s+[Nn][Aa][Nn]\s+"), "NAN values"),
            (re.compile("energy\s+=\s*(\*)+"), "Numerical disaster"),
            (re.compile("NewFileMan::OpenFile\(\):\s+nopenfiles=\d+\s+"
                        "maxopenfiles=\d+s+errno=\d+"), "Open file error"),
            (re.compile("Application \d+ exit codes: 1[34]\d+"), "Exit Code 134"),
            (re.compile("Negative overlap matrix eigenvalue. Tighten integral "
                        "threshold \(REM_THRESH\)!"), "Negative Eigen"),
            (re.compile("Unable to allocate requested memory in mega_alloc"),
                "Insufficient static memory"),
            (re.compile("Application \d+ exit signals: Killed"),
                "Killed"),
            (re.compile("UNABLE TO DETERMINE Lamda IN FormD"),
                "Lamda Determination Failed"),
            (re.compile("Job too small. Please specify .*CPSCF_NSEG"),
                "Freq Job Too Small"),
            (re.compile("Not enough total memory"),
                "Not Enough Total Memory"),
            (re.compile("Use of \$pcm_solvent section has been deprecated starting in Q-Chem"),
                "pcm_solvent deprecated")
        )

        energies = []
        scf_iters = []
        coords = []
        species = []
        molecules = []
        gradients = []
        freqs = []
        vib_freqs = []
        vib_modes = []
        grad_comp = None
        errors = []
        parse_input = False
        parse_coords = False
        parse_scf_iter = False
        parse_gradient = False
        parse_freq = False
        parse_modes = False
        qctask_lines = []
        qctask = None
        jobtype = None
        charge = None
        spin_multiplicity = None
        thermal_corr = dict()
        properly_terminated = False
        pop_method = None
        parse_charge = False
        nbo_available = False
        nbo_charge_header = None
        parse_nbo_charge = False
        charges = dict()
        scf_successful = False
        opt_successful = False
        parse_alpha_homo = False
        parse_alpha_lumo = False
        parse_beta_homo = False
        parse_beta_lumo = False
        current_alpha_homo = None
        current_alpha_lumo = None
        current_beta_homo = None
        homo_lumo = []
        bsse = None
        hiershfiled_pop = False
        for line in output.split("\n"):
            for ep, message in error_defs:
                if ep.search(line):
                    if message == "NAN values":
                        if "time" in line:
                            continue
                    errors.append(message)
            if parse_input:
                if "-" * 50 in line:
                    if len(qctask_lines) == 0:
                        continue
                    else:
                        qctask = QcTask.from_string('\n'.join(qctask_lines))
                        jobtype = qctask.params["rem"]["jobtype"]
                        parse_input = False
                        continue
                qctask_lines.append(line)
            elif parse_coords:
                if "-" * 50 in line:
                    if len(coords) == 0:
                        continue
                    else:
                        if qctask and qctask.ghost_atoms:
                            if isinstance(qctask.mol, Molecule):
                                for i in qctask.ghost_atoms:
                                    species[i] = qctask.mol.sites[i].specie.symbol
                        molecules.append(Molecule(species, coords))
                        coords = []
                        species = []
                        parse_coords = False
                        continue
                if "Atom" in line:
                    continue
                m = coord_pattern.match(line)
                coords.append([float(m.group("x")), float(m.group("y")),
                              float(m.group("z"))])
                species.append(m.group("element"))
            elif parse_scf_iter:
                if "SCF time:  CPU" in line:
                    parse_scf_iter = False
                    continue
                if 'Convergence criterion met' in line:
                    scf_successful = True
                m = scf_iter_pattern.search(line)
                if m:
                    scf_iters[-1].append((float(m.group("energy")),
                                          float(m.group("diis_error"))))
            elif parse_gradient:
                if "Max gradient component" in line:
                    gradients[-1]["max_gradient"] = \
                        float(line.split("=")[1])
                    if grad_comp:
                        if len(grad_comp) == 3:
                            gradients[-1]["gradients"].extend(zip(*grad_comp))
                        else:
                            raise Exception("Gradient section parsing failed")
                    continue
                elif "RMS gradient" in line:
                    gradients[-1]["rms_gradient"] = \
                        float(line.split("=")[1])
                    parse_gradient = False
                    grad_comp = None
                    continue
                elif "." not in line:
                    if grad_comp:
                        if len(grad_comp) == 3:
                            gradients[-1]["gradients"].extend(zip(*grad_comp))
                        else:
                            raise Exception("Gradient section parsing failed")
                    grad_comp = []
                else:
                    grad_line_token = list(line)
                    grad_crowd = False
                    grad_line_final = line
                    for i in range(5, len(line), 12):
                        c = grad_line_token[i]
                        if not c.isspace():
                            grad_crowd = True
                            if ' ' in grad_line_token[i+1: i+6+1] or \
                                    len(grad_line_token[i+1: i+6+1]) < 6:
                                continue
                            grad_line_token[i-1] = ' '
                    if grad_crowd:
                        grad_line_final = ''.join(grad_line_token)
                    grad_comp.append([float(x) for x
                                      in grad_line_final.strip().split()[1:]])
            elif parse_freq:
                if parse_modes:
                    if "TransDip" in line:
                        parse_modes = False
                        for freq, mode in zip(vib_freqs, zip(*vib_modes)):
                            freqs.append({"frequency": freq,
                                          "vib_mode": mode})
                        vib_modes = []
                        continue
                    dis_flat = [float(x) for x in line.strip().split()[1:]]
                    dis_atom = zip(*([iter(dis_flat)]*3))
                    vib_modes.append(dis_atom)
                if "STANDARD THERMODYNAMIC QUANTITIES" in line\
                        or "Imaginary Frequencies" in line:
                    parse_freq = False
                    continue
                if "Frequency:" in line:
                    vib_freqs = [float(vib) for vib
                                 in line.strip().strip().split()[1:]]
                elif "X      Y      Z" in line:
                    parse_modes = True
                    continue
            elif parse_charge:
                if '-'*20 in line:
                    if len(charges[pop_method]) == 0:
                        continue
                    else:
                        pop_method = None
                        parse_charge = False
                else:
                    if len(line.strip()) == 0 or\
                            'Atom' in line:
                        continue
                    else:
                        charges[pop_method].append(float(line.split()[2]))
            elif parse_nbo_charge:
                if '-'*20 in line:
                    if len(charges[pop_method]) == 0:
                        continue
                elif "="*20 in line:
                    pop_method = None
                    parse_nbo_charge = False
                else:
                    m = nbo_charge_pattern.search(line)
                    if m:
                        charges[pop_method].append(float(m.group("charge")))
                    else:
                        raise Exception("Can't find NBO charges")
            elif parse_alpha_homo:
                if "-- Occupied --" in line:
                    continue
                elif "-- Virtual --" in line:
                    parse_alpha_homo = False
                    parse_alpha_lumo = True
                    continue
                else:
                    tokens = line.split()
                    m = float_pattern.search(tokens[-1])
                    if m:
                        current_alpha_homo = float(m.group(0))
                    continue
            elif parse_alpha_lumo:
                current_alpha_lumo = float(line.split()[0])
                parse_alpha_lumo = False
                continue
            elif parse_beta_homo:
                if "-- Occupied --" in line:
                    continue
                elif "-- Virtual --" in line:
                    parse_beta_homo = False
                    parse_beta_lumo = True
                    continue
                else:
                    tokens = line.split()
                    m = float_pattern.search(tokens[-1])
                    if m:
                        current_beta_homo = float(m.group(0))
                    continue
            elif parse_beta_lumo:
                current_beta_lumo = float(line.split()[0])
                parse_beta_lumo = False
                if isinstance(current_alpha_homo, float) and isinstance(current_beta_homo, float):
                    current_homo = max([current_alpha_homo, current_beta_homo])
                else:
                    current_homo = 0.0
                if isinstance(current_alpha_lumo, float) and isinstance(current_beta_lumo, float):
                    current_lumo = min([current_alpha_lumo, current_beta_lumo])
                else:
                    current_lumo = 0.0
                homo_lumo.append([Energy(current_homo, "Ha").to("eV"),
                                  Energy(current_lumo, "Ha").to("eV")])
                current_alpha_homo = None
                current_alpha_lumo = None
                current_beta_homo = None
                continue
            elif "-" * 50 in line and not (current_alpha_lumo is None):
                homo_lumo.append([Energy(current_alpha_homo, "Ha").to("eV"),
                                  Energy(current_alpha_lumo, "Ha").to("eV")])
                current_alpha_homo = None
                current_alpha_lumo = None
                current_beta_homo = None
                continue
            else:
                if spin_multiplicity is None:
                    m = num_ele_pattern.search(line)
                    if m:
                        spin_multiplicity = int(m.group("alpha")) - \
                            int(m.group("beta")) + 1
                if charge is None:
                    m = total_charge_pattern.search(line)
                    if m:
                        charge = int(float(m.group("charge")))
                if jobtype and jobtype == "freq":
                    m = zpe_pattern.search(line)
                    if m:
                        zpe = float(m.group("zpe"))
                        thermal_corr["ZPE"] = zpe
                    m = thermal_corr_pattern.search(line)
                    if m:
                        thermal_corr[m.group("name")] = \
                            float(m.group("correction"))
                m = bsse_pattern.search(line)
                if m:
                    raw_be = float(m.group("raw_be"))
                    corrected_be = float(m.group("corrected_be"))
                    bsse_fwu = FloatWithUnit(raw_be - corrected_be, "kJ mol^-1")
                    bsse = bsse_fwu.to('eV atom^-1').real
                name = None
                energy = None
                m = scf_energy_pattern.search(line)
                if m:
                    name = "SCF"
                    energy = Energy(m.group("energy"), "Ha").to("eV")
                m = corr_energy_pattern.search(line)
                if m and m.group("name") != "SCF":
                    name = m.group("name")
                    energy = Energy(m.group("energy"), "Ha").to("eV")
                m = detailed_charge_pattern.search(line)
                if m:
                    pop_method = m.group("method").lower()
                    parse_charge = True
                    charges[pop_method] = []
                if nbo_available:
                    if nbo_charge_header is None:
                        m = nbo_wavefunction_type_pattern.search(line)
                        if m:
                            nbo_wavefunction_type = m.group("type")
                            nbo_charge_header_dict = {
                                "closed-shell": "Atom No    Charge        Core      "
                                                "Valence    Rydberg      Total",
                                "open-shell": "Atom No    Charge        Core      "
                                              "Valence    Rydberg      Total      Density"}
                            nbo_charge_header = nbo_charge_header_dict[nbo_wavefunction_type]
                        continue
                    if nbo_charge_header in line:
                        pop_method = "nbo"
                        parse_nbo_charge = True
                        charges[pop_method] = []
                if "N A T U R A L   B O N D   O R B I T A L   A N A L Y S I S" in line:
                    nbo_available = True
                if name and energy:
                    energies.append(tuple([name, energy]))
                if "User input:" in line:
                    parse_input = True
                elif "Standard Nuclear Orientation (Angstroms)" in line:
                    parse_coords = True
                elif "Performing Hirshfeld population analysis" in line:
                    hiershfiled_pop = True
                elif "Hirshfeld: atomic densities completed" in line:
                    hiershfiled_pop = False
                elif ("Cycle       Energy         DIIS Error" in line
                        or "Cycle       Energy        RMS Gradient" in line)\
                        and not hiershfiled_pop:
                    parse_scf_iter = True
                    scf_iters.append([])
                    scf_successful = False
                elif "Gradient of SCF Energy" in line:
                    parse_gradient = True
                    gradients.append({"gradients": []})
                elif "VIBRATIONAL ANALYSIS" in line:
                    parse_freq = True
                elif "Alpha MOs" in line:
                    parse_alpha_homo = True
                    parse_alpha_lumo = False
                elif "Beta MOs" in line:
                    parse_beta_homo = True
                    parse_beta_lumo = False
                elif "Thank you very much for using Q-Chem." in line:
                    properly_terminated = True
                elif "OPTIMIZATION CONVERGED" in line:
                    opt_successful = True
        if charge is None:
            errors.append("Molecular charge is not found")
        elif spin_multiplicity is None:
            errors.append("Molecular spin multipilicity is not found")
        else:
            for mol in molecules:
                if qctask is None or qctask.ghost_atoms is None:
                    mol.set_charge_and_spin(charge, spin_multiplicity)
        for k in thermal_corr.keys():
            v = thermal_corr[k]
            if "Entropy" in k:
                v *= cls.kcal_per_mol_2_eV * 1.0E-3
            else:
                v *= cls.kcal_per_mol_2_eV
            thermal_corr[k] = v

        solvent_method = "NA"
        if qctask:
            if "solvent_method" in qctask.params["rem"]:
                solvent_method = qctask.params["rem"]["solvent_method"]
        else:
            errors.append("No input text")

        if not scf_successful:
            if 'Bad SCF convergence' not in errors:
                errors.append('Bad SCF convergence')

        if jobtype == 'opt':
            if not opt_successful:
                if 'Geometry optimization failed' not in errors:
                    errors.append('Geometry optimization failed')

        if len(errors) == 0:
            for text in cls._expected_successful_pattern(qctask):
                success_pattern = re.compile(text)
                if not success_pattern.search(output):
                    errors.append("Can't find text to indicate success")

        data = {
            "jobtype": jobtype,
            "energies": energies,
            "HOMO/LUMOs": homo_lumo,
            "bsse": bsse,
            'charges': charges,
            "corrections": thermal_corr,
            "molecules": molecules,
            "errors": errors,
            "has_error": len(errors) > 0,
            "frequencies": freqs,
            "gradients": gradients,
            "input": qctask,
            "gracefully_terminated": properly_terminated,
            "scf_iteration_energies": scf_iters,
            "solvent_method": solvent_method
        }
        return data
