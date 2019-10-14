# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements input and output processing from Gaussian.
"""

import re

import numpy as np
import warnings

from pymatgen.core.operations import SymmOp
from pymatgen import Element, Molecule, Composition
from monty.io import zopen
from pymatgen.core.units import Ha_to_eV
from pymatgen.util.coord import get_angle
import scipy.constants as cst

from pymatgen.electronic_structure.core import Spin

__author__ = 'Shyue Ping Ong, Germain  Salvato-Vallverdu, Xin Chen'
__copyright__ = 'Copyright 2013, The Materials Virtual Lab'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = '8/1/15'


float_patt = re.compile(r"\s*([+-]?\d+\.\d+)")


def read_route_line(route):
    """
    read route line in gaussian input/output and return functional basis_set
    and a dictionary of other route parameters

    Args:
        route (str) : the route line

    return
        functional (str) : the method (HF, PBE ...)
        basis_set (str) : the basis set
        route (dict) : dictionary of parameters
    """
    scrf_patt = re.compile(r"^([sS][cC][rR][fF])\s*=\s*(.+)")
    multi_params_patt = re.compile(r"^([A-z]+[0-9]*)[\s=]+\((.*)\)$")
    functional = None
    basis_set = None
    route_params = {}
    dieze_tag = None
    if route:
        if "/" in route:
            tok = route.split("/")
            functional = tok[0].split()[-1]
            basis_set = tok[1].split()[0]
            for tok in [functional, basis_set, "/"]:
                route = route.replace(tok, "")

        for tok in route.split():
            if scrf_patt.match(tok):
                m = scrf_patt.match(tok)
                route_params[m.group(1)] = m.group(2)
            elif tok.upper() in ["#", "#N", "#P", "#T"]:
                # does not store # in route to avoid error in input
                if tok == "#":
                    dieze_tag = "#N"
                else:
                    dieze_tag = tok
                continue
            else:
                m = re.match(multi_params_patt, tok.strip("#"))
                if m:
                    pars = {}
                    for par in m.group(2).split(","):
                        p = par.split("=")
                        pars[p[0]] = None if len(p) == 1 else p[1]
                    route_params[m.group(1)] = pars
                else:
                    d = tok.strip("#").split("=")
                    route_params[d[0]] = None if len(d) == 1 else d[1]

    return functional, basis_set, route_params, dieze_tag


class GaussianInput:
    """
    An object representing a Gaussian input file.
    """

    # Commonly used regex patterns
    _zmat_patt = re.compile(r"^(\w+)*([\s,]+(\w+)[\s,]+(\w+))*[\-\.\s,\w]*$")
    _xyz_patt = re.compile(r"^(\w+)[\s,]+([\d\.eE\-]+)[\s,]+([\d\.eE\-]+)[\s,]+"
                           r"([\d\.eE\-]+)[\-\.\s,\w.]*$")

    def __init__(self, mol, charge=None, spin_multiplicity=None, title=None,
                 functional="HF", basis_set="6-31G(d)", route_parameters=None,
                 input_parameters=None, link0_parameters=None, dieze_tag="#P",
                 gen_basis=None):
        """
        Args:
            mol: Input molecule. If molecule is a single string, it is used as a
                direct input to the geometry section of the Gaussian input
                file.
            charge: Charge of the molecule. If None, charge on molecule is used.
                Defaults to None. This allows the input file to be set a
                charge independently from the molecule itself.
            spin_multiplicity: Spin multiplicity of molecule. Defaults to None,
                which means that the spin multiplicity is set to 1 if the
                molecule has no unpaired electrons and to 2 if there are
                unpaired electrons.
            title: Title for run. Defaults to formula of molecule if None.
            functional: Functional for run.
            basis_set: Basis set for run.
            route_parameters: Additional route parameters as a dict. For example,
                {'SP':"", "SCF":"Tight"}
            input_parameters: Additional input parameters for run as a dict. Used
                for example, in PCM calculations.  E.g., {"EPS":12}
            link0_parameters: Link0 parameters as a dict. E.g., {"%mem": "1000MW"}
            dieze_tag: # preceding the route line. E.g. "#p"
            gen_basis: allows a user-specified basis set to be used in a Gaussian
                calculation. If this is not None, the attribute ``basis_set`` will
                be set to "Gen".
        """
        self._mol = mol
        self.charge = charge if charge is not None else mol.charge
        nelectrons = - self.charge + mol.charge + mol.nelectrons
        if spin_multiplicity is not None:
            self.spin_multiplicity = spin_multiplicity
            if (nelectrons + spin_multiplicity) % 2 != 1:
                raise ValueError(
                    "Charge of {} and spin multiplicity of {} is"
                    " not possible for this molecule".format(
                        self.charge, spin_multiplicity))
        else:
            self.spin_multiplicity = 1 if nelectrons % 2 == 0 else 2
        self.functional = functional
        self.basis_set = basis_set
        self.link0_parameters = link0_parameters if link0_parameters else {}
        self.route_parameters = route_parameters if route_parameters else {}
        self.input_parameters = input_parameters if input_parameters else {}
        self.title = title if title else self._mol.composition.formula
        self.dieze_tag = dieze_tag if dieze_tag[0] == "#" else "#" + dieze_tag
        self.gen_basis = gen_basis
        if gen_basis is not None:
            self.basis_set = "Gen"

    @property
    def molecule(self):
        """
        Returns molecule associated with this GaussianInput.
        """
        return self._mol

    @staticmethod
    def _parse_coords(coord_lines):
        """
        Helper method to parse coordinates.
        """
        paras = {}
        var_pattern = re.compile(r"^([A-Za-z]+\S*)[\s=,]+([\d\-\.]+)$")
        for l in coord_lines:
            m = var_pattern.match(l.strip())
            if m:
                paras[m.group(1).strip("=")] = float(m.group(2))

        species = []
        coords = []
        # Stores whether a Zmatrix format is detected. Once a zmatrix format
        # is detected, it is assumed for the remaining of the parsing.
        zmode = False
        for l in coord_lines:
            l = l.strip()
            if not l:
                break
            if (not zmode) and GaussianInput._xyz_patt.match(l):
                m = GaussianInput._xyz_patt.match(l)
                species.append(m.group(1))
                toks = re.split(r"[,\s]+", l.strip())
                if len(toks) > 4:
                    coords.append([float(i) for i in toks[2:5]])
                else:
                    coords.append([float(i) for i in toks[1:4]])
            elif GaussianInput._zmat_patt.match(l):
                zmode = True
                toks = re.split(r"[,\s]+", l.strip())
                species.append(toks[0])
                toks.pop(0)
                if len(toks) == 0:
                    coords.append(np.array([0, 0, 0]))
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
                        coords.append(np.array([0, 0, parameters[0]]))
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

        def _parse_species(sp_str):
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
                sp = re.sub(r"\d", "", sp_str)
                return sp.capitalize()

        species = [_parse_species(sp) for sp in species]

        return Molecule(species, coords)

    @staticmethod
    def from_string(contents):
        """
        Creates GaussianInput from a string.

        Args:
            contents: String representing an Gaussian input file.

        Returns:
            GaussianInput object
        """
        lines = [l.strip() for l in contents.split("\n")]

        link0_patt = re.compile(r"^(%.+)\s*=\s*(.+)")
        link0_dict = {}
        for i, l in enumerate(lines):
            if link0_patt.match(l):
                m = link0_patt.match(l)
                link0_dict[m.group(1).strip("=")] = m.group(2)

        route_patt = re.compile(r"^#[sSpPnN]*.*")
        route = ""
        route_index = None
        for i, l in enumerate(lines):
            if route_patt.match(l):
                route += " " + l
                route_index = i
            # This condition allows for route cards spanning multiple lines
            elif (l == "" or l.isspace()) and route_index:
                break
        functional, basis_set, route_paras, dieze_tag = read_route_line(route)
        ind = 2
        title = []
        while lines[route_index + ind].strip():
            title.append(lines[route_index + ind].strip())
            ind += 1
        title = ' '.join(title)
        ind += 1
        toks = re.split(r"[,\s]+", lines[route_index + ind])
        charge = int(float(toks[0]))
        spin_mult = int(toks[1])
        coord_lines = []
        spaces = 0
        input_paras = {}
        ind += 1
        for i in range(route_index + ind, len(lines)):
            if lines[i].strip() == "":
                spaces += 1
            if spaces >= 2:
                d = lines[i].split("=")
                if len(d) == 2:
                    input_paras[d[0]] = d[1]
            else:
                coord_lines.append(lines[i].strip())
        mol = GaussianInput._parse_coords(coord_lines)
        mol.set_charge_and_spin(charge, spin_mult)

        return GaussianInput(mol, charge=charge, spin_multiplicity=spin_mult,
                             title=title, functional=functional,
                             basis_set=basis_set,
                             route_parameters=route_paras,
                             input_parameters=input_paras,
                             link0_parameters=link0_dict,
                             dieze_tag=dieze_tag)

    @staticmethod
    def from_file(filename):
        """
        Creates GaussianInput from a file.

        Args:
            filename: Gaussian input filename

        Returns:
            GaussianInput object
        """
        with zopen(filename, "r") as f:
            return GaussianInput.from_string(f.read())

    def _find_nn_pos_before_site(self, siteindex):
        """
        Returns index of nearest neighbor atoms.
        """
        alldist = [(self._mol.get_distance(siteindex, i), i)
                   for i in range(siteindex)]
        alldist = sorted(alldist, key=lambda x: x[0])
        return [d[1] for d in alldist]

    def get_zmatrix(self):
        """
        Returns a z-matrix representation of the molecule.
        """
        output = []
        outputvar = []
        for i, site in enumerate(self._mol):
            if i == 0:
                output.append("{}".format(site.specie))
            elif i == 1:
                nn = self._find_nn_pos_before_site(i)
                bondlength = self._mol.get_distance(i, nn[0])
                output.append("{} {} B{}".format(self._mol[i].specie,
                                                 nn[0] + 1, i))
                outputvar.append("B{}={:.6f}".format(i, bondlength))
            elif i == 2:
                nn = self._find_nn_pos_before_site(i)
                bondlength = self._mol.get_distance(i, nn[0])
                angle = self._mol.get_angle(i, nn[0], nn[1])
                output.append("{} {} B{} {} A{}".format(self._mol[i].specie,
                                                        nn[0] + 1, i,
                                                        nn[1] + 1, i))
                outputvar.append("B{}={:.6f}".format(i, bondlength))
                outputvar.append("A{}={:.6f}".format(i, angle))
            else:
                nn = self._find_nn_pos_before_site(i)
                bondlength = self._mol.get_distance(i, nn[0])
                angle = self._mol.get_angle(i, nn[0], nn[1])
                dih = self._mol.get_dihedral(i, nn[0], nn[1], nn[2])
                output.append("{} {} B{} {} A{} {} D{}"
                              .format(self._mol[i].specie, nn[0] + 1, i,
                                      nn[1] + 1, i, nn[2] + 1, i))
                outputvar.append("B{}={:.6f}".format(i, bondlength))
                outputvar.append("A{}={:.6f}".format(i, angle))
                outputvar.append("D{}={:.6f}".format(i, dih))
        return "\n".join(output) + "\n\n" + "\n".join(outputvar)

    def get_cart_coords(self):
        """
        Return the cartesian coordinates of the molecule
        """
        def to_s(x):
            return "%0.6f" % x

        outs = []
        for i, site in enumerate(self._mol):
            outs.append(" ".join([site.species_string,
                                  " ".join([to_s(j) for j in site.coords])]))
        return "\n".join(outs)

    def __str__(self):
        return self.to_string()

    def to_string(self, cart_coords=False):
        """
        Return GaussianInput string

        Option: whe cart_coords sets to True return the cartesian coordinates
                instead of the z-matrix

        """
        def para_dict_to_string(para, joiner=" "):
            para_str = []
            # sorted is only done to make unittests work reliably
            for par, val in sorted(para.items()):
                if val is None or val == "":
                    para_str.append(par)
                elif isinstance(val, dict):
                    val_str = para_dict_to_string(val, joiner=",")
                    para_str.append("{}=({})".format(par, val_str))
                else:
                    para_str.append("{}={}".format(par, val))
            return joiner.join(para_str)

        output = []
        if self.link0_parameters:
            output.append(para_dict_to_string(self.link0_parameters, "\n"))
        output.append("{diez} {func}/{bset} {route}"
                      .format(diez=self.dieze_tag, func=self.functional,
                              bset=self.basis_set,
                              route=para_dict_to_string(self.route_parameters))
                      )
        output.append("")
        output.append(self.title)
        output.append("")
        output.append("%d %d" % (self.charge, self.spin_multiplicity))
        if isinstance(self._mol, Molecule):
            if cart_coords is True:
                output.append(self.get_cart_coords())
            else:
                output.append(self.get_zmatrix())
        else:
            output.append(str(self._mol))
        output.append("")
        if self.gen_basis is not None:
            output.append("{:s}\n".format(self.gen_basis))
        output.append(para_dict_to_string(self.input_parameters, "\n"))
        output.append("\n")
        return "\n".join(output)

    def write_file(self, filename, cart_coords=False):
        """
        Write the input string into a file

        Option: see __str__ method
        """
        with zopen(filename, "w") as f:
            f.write(self.to_string(cart_coords))

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "molecule": self.molecule.as_dict(),
                "functional": self.functional,
                "basis_set": self.basis_set,
                "route_parameters": self.route_parameters,
                "title": self.title,
                "charge": self.charge,
                "spin_multiplicity": self.spin_multiplicity,
                "input_parameters": self.input_parameters,
                "link0_parameters": self.link0_parameters,
                "dieze_tag": self.dieze_tag}

    @classmethod
    def from_dict(cls, d):
        """
        :param d: dict
        :return: GaussianInput
        """
        return GaussianInput(mol=Molecule.from_dict(d["molecule"]),
                             functional=d["functional"],
                             basis_set=d["basis_set"],
                             route_parameters=d["route_parameters"],
                             title=d["title"],
                             charge=d["charge"],
                             spin_multiplicity=d["spin_multiplicity"],
                             input_parameters=d["input_parameters"],
                             link0_parameters=d["link0_parameters"])


class GaussianOutput:
    """
    Parser for Gaussian output files.

    .. note::

        Still in early beta.

    Attributes:

    .. attribute:: structures

        All structures from the calculation in the standard orientation. If the
        symmetry is not considered, the standard orientation is not printed out
        and the input orientation is used instead. Check the `standard_orientation`
        attribute.

    .. attribute:: structures_input_orientation

        All structures from the calculation in the input orientation.

    .. attribute:: energies

        All energies from the calculation.

    .. attribute:: eigenvalues

        List of eigenvalues for the last geometry

    .. attribute:: MO_coefficients

        Matrix of MO coefficients for the last geometry

    .. attribute:: cart_forces

        All cartesian forces from the calculation.

    .. attribute:: frequencies

        A list for each freq calculation and for each mode of a dict with
        {
            "frequency": freq in cm-1,
            "symmetry": symmetry tag
            "r_mass": Reduce mass,
            "f_constant": force constant,
            "IR_intensity": IR Intensity,
            "mode": normal mode
         }

        The normal mode is a 1D vector of dx, dy dz of each atom.

    .. attribute:: hessian

        Matrix of second derivatives of the energy with respect to cartesian
        coordinates in the **input orientation** frame. Need #P in the
        route section in order to be in the output.

    .. attribute:: properly_terminated

        True if run has properly terminated

    .. attribute:: is_pcm

        True if run is a PCM run.

    .. attribute:: is_spin

        True if it is an unrestricted run

    .. attribute:: stationary_type

        If it is a relaxation run, indicates whether it is a minimum (Minimum)
        or a saddle point ("Saddle").

    .. attribute:: corrections

        Thermochemical corrections if this run is a Freq run as a dict. Keys
        are "Zero-point", "Thermal", "Enthalpy" and "Gibbs Free Energy"

    .. attribute:: functional

        Functional used in the run.

    .. attribute:: basis_set

        Basis set used in the run

    .. attribute:: route

        Additional route parameters as a dict. For example,
            {'SP':"", "SCF":"Tight"}

    .. attribute:: dieze_tag

        # preceding the route line, e.g. "#P"

    .. attribute:: link0

        Link0 parameters as a dict. E.g., {"%mem": "1000MW"}

    .. attribute:: charge

        Charge for structure

    .. attribute:: spin_multiplicity

        Spin multiplicity for structure

    .. attribute:: num_basis_func

        Number of basis functions in the run.

    .. attribute:: electrons

        number of alpha and beta electrons as (N alpha, N beta)

    .. attribute:: pcm

        PCM parameters and output if available.

    .. attribute:: errors

        error if not properly terminated (list to be completed in error_defs)

    .. attribute:: Mulliken_charges

        Mulliken atomic charges

    .. attribute:: eigenvectors

        Matrix of shape (num_basis_func, num_basis_func). Each column is an
        eigenvectors and contains AO coefficients of an MO.

        eigenvectors[Spin] = mat(num_basis_func, num_basis_func)

    .. attribute:: molecular_orbital

        MO development coefficients on AO in a more convenient array dict
        for each atom and basis set label.

        mo[Spin][OM j][atom i] = {AO_k: coeff, AO_k: coeff ... }

    .. attribute:: atom_basis_labels

        Labels of AO for each atoms. These labels are those used in the output
        of molecular orbital coefficients (POP=Full) and in the
        molecular_orbital array dict.

        atom_basis_labels[iatom] = [AO_k, AO_k, ...]

    .. attribute:: resumes

        List of gaussian data resume given at the end of the output file before
        the quotation. The resumes are given as string.

    .. attribute:: title

        Title of the gaussian run.

    .. attribute:: standard_orientation

        If True, the geometries stored in the structures are in the standard
        orientation. Else, the geometries are in the input orientation.

    .. attribute:: bond_orders

        Dict of bond order values read in the output file such as:
        {(0, 1): 0.8709, (1, 6): 1.234, ...}

        The keys are the atom indexes and the values are the Wiberg bond indexes
        that are printed using `pop=NBOREAD` and `$nbo bndidx $end`.

    Methods:

    .. method:: to_input()

        Return a GaussianInput object using the last geometry and the same
        calculation parameters.

    .. method:: read_scan()

        Read a potential energy surface from a gaussian scan calculation.

    .. method:: get_scan_plot()

        Get a matplotlib plot of the potential energy surface

    .. method:: save_scan_plot()

        Save a matplotlib plot of the potential energy surface to a file

    """

    def __init__(self, filename):
        """
        Args:
            filename: Filename of Gaussian output file.
        """
        self.filename = filename
        self._parse(filename)

    @property
    def final_energy(self):
        """
        :return: Final energy in Gaussian output.
        """
        return self.energies[-1]

    @property
    def final_structure(self):
        """
        :return: Final structure in Gaussian output.
        """
        return self.structures[-1]

    def _parse(self, filename):
        start_patt = re.compile(r" \(Enter \S+l101\.exe\)")
        route_patt = re.compile(r" #[pPnNtT]*.*")
        link0_patt = re.compile(r"^\s(%.+)\s*=\s*(.+)")
        charge_mul_patt = re.compile(r"Charge\s+=\s*([-\d]+)\s+"
                                     r"Multiplicity\s+=\s*(\d+)")
        num_basis_func_patt = re.compile(r"([0-9]+)\s+basis functions")
        num_elec_patt = re.compile(r"(\d+)\s+alpha electrons\s+(\d+)\s+beta electrons")
        pcm_patt = re.compile(r"Polarizable Continuum Model")
        stat_type_patt = re.compile(r"imaginary frequencies")
        scf_patt = re.compile(r"E\(.*\)\s*=\s*([-\.\d]+)\s+")
        mp2_patt = re.compile(r"EUMP2\s*=\s*(.*)")
        oniom_patt = re.compile(r"ONIOM:\s+extrapolated energy\s*=\s*(.*)")
        termination_patt = re.compile(r"(Normal|Error) termination")
        error_patt = re.compile(
            r"(! Non-Optimized Parameters !|Convergence failure)")
        mulliken_patt = re.compile(
            r"^\s*(Mulliken charges|Mulliken atomic charges)")
        mulliken_charge_patt = re.compile(
            r'^\s+(\d+)\s+([A-Z][a-z]?)\s*(\S*)')
        end_mulliken_patt = re.compile(
            r'(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)')
        std_orientation_patt = re.compile(r"Standard orientation")
        input_orientation_patt = re.compile(r"Input orientation")
        orbital_patt = re.compile(r"(Alpha|Beta)\s*\S+\s*eigenvalues --(.*)")
        thermo_patt = re.compile(r"(Zero-point|Thermal) correction(.*)="
                                 r"\s+([\d\.-]+)")
        forces_on_patt = re.compile(
            r"Center\s+Atomic\s+Forces\s+\(Hartrees/Bohr\)")
        forces_off_patt = re.compile(r"Cartesian\s+Forces:\s+Max.*RMS.*")
        forces_patt = re.compile(
            r"\s+(\d+)\s+(\d+)\s+([0-9\.-]+)\s+([0-9\.-]+)\s+([0-9\.-]+)")

        freq_on_patt = re.compile(
            r"Harmonic\sfrequencies\s+\(cm\*\*-1\),\sIR\sintensities.*Raman.*")

        normal_mode_patt = re.compile(
            r"\s+(\d+)\s+(\d+)\s+([0-9\.-]{4,5})\s+([0-9\.-]{4,5}).*")

        mo_coeff_patt = re.compile(r"Molecular Orbital Coefficients:")
        mo_coeff_name_patt = re.compile(r"\d+\s((\d+|\s+)\s+([a-zA-Z]{1,2}|\s+))\s+(\d+\S+)")

        hessian_patt = re.compile(r"Force constants in Cartesian coordinates:")
        resume_patt = re.compile(r"^\s1\\1\\GINC-\S*")
        resume_end_patt = re.compile(r"^\s.*\\\\@")

        bond_order_patt = re.compile(r"Wiberg bond index matrix in the NAO basis:")

        self.properly_terminated = False
        self.is_pcm = False
        self.stationary_type = "Minimum"
        self.corrections = {}
        self.energies = []
        self.pcm = None
        self.errors = []
        self.Mulliken_charges = {}
        self.link0 = {}
        self.cart_forces = []
        self.frequencies = []
        self.eigenvalues = []
        self.is_spin = False
        self.hessian = None
        self.resumes = []
        self.title = None
        self.bond_orders = {}

        read_coord = 0
        read_mulliken = False
        read_eigen = False
        eigen_txt = []
        parse_stage = 0
        num_basis_found = False
        terminated = False
        parse_forces = False
        forces = []
        parse_freq = False
        frequencies = []
        read_mo = False
        parse_hessian = False
        routeline = ""
        standard_orientation = False
        parse_bond_order = False
        input_structures = list()
        std_structures = list()
        geom_orientation = None

        with zopen(filename) as f:
            for line in f:
                if parse_stage == 0:
                    if start_patt.search(line):
                        parse_stage = 1
                    elif link0_patt.match(line):
                        m = link0_patt.match(line)
                        self.link0[m.group(1)] = m.group(2)
                    elif route_patt.search(line) or routeline != "":
                        if set(line.strip()) == {"-"}:
                            params = read_route_line(routeline)
                            self.functional = params[0]
                            self.basis_set = params[1]
                            self.route_parameters = params[2]
                            route_lower = {k.lower(): v
                                           for k, v in
                                           self.route_parameters.items()}
                            self.dieze_tag = params[3]
                            parse_stage = 1
                        else:
                            routeline += line.strip()
                elif parse_stage == 1:
                    if set(line.strip()) == {"-"} and self.title is None:
                        self.title = ""
                    elif self.title == "":
                        self.title = line.strip()
                    elif charge_mul_patt.search(line):
                        m = charge_mul_patt.search(line)
                        self.charge = int(m.group(1))
                        self.spin_multiplicity = int(m.group(2))
                        parse_stage = 2
                elif parse_stage == 2:

                    if self.is_pcm:
                        self._check_pcm(line)

                    if "freq" in route_lower and thermo_patt.search(line):
                        m = thermo_patt.search(line)
                        if m.group(1) == "Zero-point":
                            self.corrections["Zero-point"] = float(m.group(3))
                        else:
                            key = m.group(2).strip(" to ")
                            self.corrections[key] = float(m.group(3))

                    if read_coord:
                        [f.readline() for i in range(3)]
                        line = f.readline()
                        sp = []
                        coords = []
                        while set(line.strip()) != {"-"}:
                            toks = line.split()
                            sp.append(Element.from_Z(int(toks[1])))
                            coords.append([float(x) for x in toks[3:6]])
                            line = f.readline()

                        read_coord = False
                        if geom_orientation == "input":
                            input_structures.append(Molecule(sp, coords))
                        elif geom_orientation == "standard":
                            std_structures.append(Molecule(sp, coords))

                    if parse_forces:
                        m = forces_patt.search(line)
                        if m:
                            forces.extend([float(_v)
                                           for _v in m.groups()[2:5]])
                        elif forces_off_patt.search(line):
                            self.cart_forces.append(forces)
                            forces = []
                            parse_forces = False

                    # read molecular orbital eigenvalues
                    if read_eigen:
                        m = orbital_patt.search(line)
                        if m:
                            eigen_txt.append(line)
                        else:
                            read_eigen = False
                            self.eigenvalues = {Spin.up: []}
                            for eigenline in eigen_txt:
                                if "Alpha" in eigenline:
                                    self.eigenvalues[Spin.up] += [float(e) for e in float_patt.findall(eigenline)]
                                elif "Beta" in eigenline:
                                    if Spin.down not in self.eigenvalues:
                                        self.eigenvalues[Spin.down] = []
                                    self.eigenvalues[Spin.down] += [float(e) for e in float_patt.findall(eigenline)]
                            eigen_txt = []

                    # read molecular orbital coefficients
                    if (not num_basis_found) and num_basis_func_patt.search(line):
                        m = num_basis_func_patt.search(line)
                        self.num_basis_func = int(m.group(1))
                        num_basis_found = True
                    elif read_mo:
                        # build a matrix with all coefficients
                        all_spin = [Spin.up]
                        if self.is_spin:
                            all_spin.append(Spin.down)

                        mat_mo = {}
                        for spin in all_spin:
                            mat_mo[spin] = np.zeros((self.num_basis_func,
                                                     self.num_basis_func))
                            nMO = 0
                            end_mo = False
                            while nMO < self.num_basis_func and not end_mo:
                                f.readline()
                                f.readline()
                                self.atom_basis_labels = []
                                for i in range(self.num_basis_func):
                                    line = f.readline()

                                    # identify atom and OA labels
                                    m = mo_coeff_name_patt.search(line)
                                    if m.group(1).strip() != "":
                                        iat = int(m.group(2)) - 1
                                        # atname = m.group(3)
                                        self.atom_basis_labels.append([m.group(4)])
                                    else:
                                        self.atom_basis_labels[iat].append(m.group(4))

                                    # MO coefficients
                                    coeffs = [float(c) for c in
                                              float_patt.findall(line)]
                                    for j in range(len(coeffs)):
                                        mat_mo[spin][i, nMO + j] = coeffs[j]

                                nMO += len(coeffs)
                                line = f.readline()
                                # manage pop=regular case (not all MO)
                                if nMO < self.num_basis_func and \
                                    ("Density Matrix:" in line or
                                     mo_coeff_patt.search(line)):
                                    end_mo = True
                                    warnings.warn("POP=regular case, matrix "
                                                  "coefficients not complete")
                            f.readline()

                        self.eigenvectors = mat_mo
                        read_mo = False

                        # build a more convenient array dict with MO
                        # coefficient of each atom in each MO.
                        # mo[Spin][OM j][atom i] =
                        # {AO_k: coeff, AO_k: coeff ... }
                        mo = {}
                        for spin in all_spin:
                            mo[spin] = [[{} for iat in
                                         range(len(self.atom_basis_labels))]
                                        for j in range(self.num_basis_func)]
                            for j in range(self.num_basis_func):
                                i = 0
                                for iat in range(len(self.atom_basis_labels)):
                                    for label in self.atom_basis_labels[iat]:
                                        mo[spin][j][iat][label] = self.eigenvectors[spin][i, j]
                                        i += 1

                        self.molecular_orbital = mo

                    elif parse_freq:
                        while line.strip() != "":  #  blank line
                            ifreqs = [int(val) - 1 for val in line.split()]
                            for ifreq in ifreqs:
                                frequencies.append({"frequency": None,
                                                    "r_mass": None,
                                                    "f_constant": None,
                                                    "IR_intensity": None,
                                                    "symmetry": None,
                                                    "mode": []})
                            # read freq, intensity, masses, symmetry ...
                            while "Atom  AN" not in line:
                                if "Frequencies --" in line:
                                    freqs = map(float,
                                                float_patt.findall(line))
                                    for ifreq, freq in zip(ifreqs, freqs):
                                        frequencies[ifreq]["frequency"] = freq
                                elif "Red. masses --" in line:
                                    r_masses = map(float,
                                                   float_patt.findall(line))
                                    for ifreq, r_mass in zip(ifreqs, r_masses):
                                        frequencies[ifreq]["r_mass"] = r_mass
                                elif "Frc consts  --" in line:
                                    f_consts = map(float,
                                                   float_patt.findall(line))
                                    for ifreq, f_const in zip(ifreqs, f_consts):
                                        frequencies[ifreq]["f_constant"] = f_const
                                elif "IR Inten    --" in line:
                                    IR_intens = map(float,
                                                    float_patt.findall(line))
                                    for ifreq, intens in zip(ifreqs, IR_intens):
                                        frequencies[ifreq]["IR_intensity"] = intens
                                else:
                                    syms = line.split()[:3]
                                    for ifreq, sym in zip(ifreqs, syms):
                                        frequencies[ifreq]["symmetry"] = sym
                                line = f.readline()

                            #  read normal modes
                            line = f.readline()
                            while normal_mode_patt.search(line):
                                values = list(map(float,
                                                  float_patt.findall(line)))
                                for i, ifreq in zip(range(0, len(values), 3),
                                                    ifreqs):
                                    frequencies[ifreq]["mode"].extend(values[i:i+3])
                                line = f.readline()

                        parse_freq = False
                        self.frequencies.append(frequencies)
                        frequencies = []

                    elif parse_hessian:
                        #  read Hessian matrix under "Force constants in Cartesian coordinates"
                        #  Hessian matrix is in the input  orientation framework
                        # WARNING : need #P in the route line
                        parse_hessian = False
                        ndf = 3 * len(input_structures[0])
                        self.hessian = np.zeros((ndf, ndf))
                        j_indices = range(5)
                        jndf = 0
                        while jndf < ndf:
                            for i in range(jndf, ndf):
                                line = f.readline()
                                vals = re.findall(r"\s*([+-]?\d+\.\d+[eEdD]?[+-]\d+)", line)
                                vals = [float(val.replace("D", "E"))
                                        for val in vals]
                                for jval, val in enumerate(vals):
                                    j = j_indices[jval]
                                    self.hessian[i, j] = val
                                    self.hessian[j, i] = val
                            jndf += len(vals)
                            line = f.readline()
                            j_indices = [j + 5 for j in j_indices]

                    elif parse_bond_order:
                        # parse Wiberg bond order
                        line = f.readline()
                        line = f.readline()
                        nat = len(input_structures[0])
                        matrix = list()
                        for iat in range(nat):
                            line = f.readline()
                            matrix.append([float(v) for v in line.split()[2:]])

                        self.bond_orders = dict()
                        for iat in range(nat):
                            for jat in range(iat + 1, nat):
                                self.bond_orders[(iat, jat)] = matrix[iat][jat]
                        parse_bond_order = False

                    elif termination_patt.search(line):
                        m = termination_patt.search(line)
                        if m.group(1) == "Normal":
                            self.properly_terminated = True
                            terminated = True
                    elif error_patt.search(line):
                        error_defs = {
                            "! Non-Optimized Parameters !": "Optimization "
                                                            "error",
                            "Convergence failure": "SCF convergence error"
                        }
                        m = error_patt.search(line)
                        self.errors.append(error_defs[m.group(1)])
                    elif num_elec_patt.search(line):
                        m = num_elec_patt.search(line)
                        self.electrons = (int(m.group(1)), int(m.group(2)))
                    elif (not self.is_pcm) and pcm_patt.search(line):
                        self.is_pcm = True
                        self.pcm = {}
                    elif "freq" in route_lower and "opt" in route_lower and \
                            stat_type_patt.search(line):
                        self.stationary_type = "Saddle"
                    elif mp2_patt.search(line):
                        m = mp2_patt.search(line)
                        self.energies.append(float(m.group(1).replace("D",
                                                                      "E")))
                    elif oniom_patt.search(line):
                        m = oniom_patt.matcher(line)
                        self.energies.append(float(m.group(1)))
                    elif scf_patt.search(line):
                        m = scf_patt.search(line)
                        self.energies.append(float(m.group(1)))
                    elif std_orientation_patt.search(line):
                        standard_orientation = True
                        geom_orientation = "standard"
                        read_coord = True
                    elif input_orientation_patt.search(line):
                        geom_orientation = "input"
                        read_coord = True
                    elif not read_eigen and orbital_patt.search(line):
                        eigen_txt.append(line)
                        read_eigen = True
                    elif mulliken_patt.search(line):
                        mulliken_txt = []
                        read_mulliken = True
                    elif not parse_forces and forces_on_patt.search(line):
                        parse_forces = True
                    elif freq_on_patt.search(line):
                        parse_freq = True
                        [f.readline() for i in range(3)]
                    elif mo_coeff_patt.search(line):
                        if "Alpha" in line:
                            self.is_spin = True
                        read_mo = True
                    elif hessian_patt.search(line):
                        parse_hessian = True
                    elif resume_patt.search(line):
                        resume = []
                        while not resume_end_patt.search(line):
                            resume.append(line)
                            line = f.readline()
                            #  security if \\@ not in one line !
                            if line == "\n":
                                break
                        resume.append(line)
                        resume = "".join([r.strip() for r in resume])
                        self.resumes.append(resume)
                    elif bond_order_patt.search(line):
                        parse_bond_order = True

                    if read_mulliken:
                        if not end_mulliken_patt.search(line):
                            mulliken_txt.append(line)
                        else:
                            m = end_mulliken_patt.search(line)
                            mulliken_charges = {}
                            for line in mulliken_txt:
                                if mulliken_charge_patt.search(line):
                                    m = mulliken_charge_patt.search(line)
                                    dic = {int(m.group(1)):
                                           [m.group(2), float(m.group(3))]}
                                    mulliken_charges.update(dic)
                            read_mulliken = False
                            self.Mulliken_charges = mulliken_charges

        # store the structures. If symmetry is considered, the standard orientation
        # is used. Else the input orientation is used.
        if standard_orientation:
            self.structures = std_structures
            self.structures_input_orientation = input_structures
        else:
            self.structures = input_structures
            self.structures_input_orientation = input_structures

        if not terminated:
            warnings.warn("\n" + self.filename +
                          ": Termination error or bad Gaussian output file !")

    def _check_pcm(self, line):
        energy_patt = re.compile(r"(Dispersion|Cavitation|Repulsion) energy"
                                 r"\s+\S+\s+=\s+(\S*)")
        total_patt = re.compile(r"with all non electrostatic terms\s+\S+\s+"
                                r"=\s+(\S*)")
        parameter_patt = re.compile(r"(Eps|Numeral density|RSolv|Eps"
                                    r"\(inf[inity]*\))\s+=\s*(\S*)")

        if energy_patt.search(line):
            m = energy_patt.search(line)
            self.pcm['{} energy'.format(m.group(1))] = float(m.group(2))
        elif total_patt.search(line):
            m = total_patt.search(line)
            self.pcm['Total energy'] = float(m.group(1))
        elif parameter_patt.search(line):
            m = parameter_patt.search(line)
            self.pcm[m.group(1)] = float(m.group(2))

    def as_dict(self):
        """
        Json-serializable dict representation.
        """
        structure = self.final_structure
        d = {"has_gaussian_completed": self.properly_terminated,
             "nsites": len(structure)}
        comp = structure.composition
        d["unit_cell_formula"] = comp.as_dict()
        d["reduced_cell_formula"] = Composition(comp.reduced_formula).as_dict()
        d["pretty_formula"] = comp.reduced_formula
        d["is_pcm"] = self.is_pcm
        d["errors"] = self.errors
        d["Mulliken_charges"] = self.Mulliken_charges

        unique_symbols = sorted(list(d["unit_cell_formula"].keys()))
        d["elements"] = unique_symbols
        d["nelements"] = len(unique_symbols)
        d["charge"] = self.charge
        d["spin_multiplicity"] = self.spin_multiplicity

        vin = {"route": self.route_parameters, "functional": self.functional,
               "basis_set": self.basis_set,
               "nbasisfunctions": self.num_basis_func,
               "pcm_parameters": self.pcm}

        d["input"] = vin

        nsites = len(self.final_structure)

        vout = {
            "energies": self.energies,
            "final_energy": self.final_energy,
            "final_energy_per_atom": self.final_energy / nsites,
            "molecule": structure.as_dict(),
            "stationary_type": self.stationary_type,
            "corrections": self.corrections
        }

        d['output'] = vout
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__

        return d

    def read_scan(self):
        """
        Read a potential energy surface from a gaussian scan calculation.

        Returns:

            A dict: {"energies": [ values ],
                     "coords": {"d1": [ values ], "A2", [ values ], ... }}

            "energies" are the energies of all points of the potential energy
            surface. "coords" are the internal coordinates used to compute the
            potential energy surface and the internal coordinates optimized,
            labelled by their name as defined in the calculation.
        """

        def floatList(l):
            """ return a list of float from a list of string """
            return [float(v) for v in l]

        scan_patt = re.compile(r"^\sSummary of the potential surface scan:")
        optscan_patt = re.compile(r"^\sSummary of Optimized Potential Surface Scan")

        # data dict return
        data = {"energies": list(), "coords": dict()}

        # read in file
        with zopen(self.filename, "r") as f:
            line = f.readline()

            while line != "":
                if optscan_patt.match(line):
                    f.readline()
                    line = f.readline()
                    endScan = False
                    while not endScan:
                        data["energies"] += floatList(float_patt.findall(line))
                        line = f.readline()
                        while not re.search(r"(^\s+(\d+)|^\s-+)", line):
                            icname = line.split()[0].strip()
                            if icname in data["coords"]:
                                data["coords"][icname] += floatList(float_patt.findall(line))
                            else:
                                data["coords"][icname] = floatList(float_patt.findall(line))
                            line = f.readline()
                        if re.search(r"^\s-+", line):
                            endScan = True
                        else:
                            line = f.readline()

                elif scan_patt.match(line):
                    line = f.readline()
                    data["coords"] = {icname: list()
                                      for icname in line.split()[1:-1]}
                    f.readline()
                    line = f.readline()
                    while not re.search(r"^\s-+", line):
                        values = floatList(line.split())
                        data["energies"].append(values[-1])
                        for i, icname in enumerate(data["coords"]):
                            data["coords"][icname].append(values[i+1])
                        line = f.readline()
                else:
                    line = f.readline()

        return data

    def get_scan_plot(self, coords=None):
        """
        Get a matplotlib plot of the potential energy surface.

        Args:
            coords: internal coordinate name to use as abcissa.
        """
        from pymatgen.util.plotting import pretty_plot

        plt = pretty_plot(12, 8)

        d = self.read_scan()

        if coords and coords in d["coords"]:
            x = d["coords"][coords]
            plt.xlabel(coords)
        else:
            x = range(len(d["energies"]))
            plt.xlabel("points")

        plt.ylabel("Energy (eV)")

        e_min = min(d["energies"])
        y = [(e - e_min) * Ha_to_eV for e in d["energies"]]

        plt.plot(x, y, "ro--")
        return plt

    def save_scan_plot(self, filename="scan.pdf",
                       img_format="pdf", coords=None):
        """
        Save matplotlib plot of the potential energy surface to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            coords: internal coordinate name to use as abcissa.
        """
        plt = self.get_scan_plot(coords)
        plt.savefig(filename, format=img_format)

    def read_excitation_energies(self):
        """
        Read a excitation energies after a TD-DFT calculation.

        Returns:

            A list: A list of tuple for each transition such as
                    [(energie (eV), lambda (nm), oscillatory strength), ... ]
        """

        transitions = list()

        # read in file
        with zopen(self.filename, "r") as f:
            line = f.readline()
            td = False
            while line != "":
                if re.search(r"^\sExcitation energies and oscillator strengths:", line):
                    td = True

                if td:
                    if re.search(r"^\sExcited State\s*\d", line):
                        val = [float(v) for v in float_patt.findall(line)]
                        transitions.append(tuple(val[0:3]))
                line = f.readline()
        return transitions

    def get_spectre_plot(self, sigma=0.05, step=0.01):
        """
        Get a matplotlib plot of the UV-visible xas. Transition are plotted
        as vertical lines and as a sum of normal functions with sigma with. The
        broadening is applied in energy and the xas is plotted as a function
        of the wavelength.

        Args:
            sigma: Full width at half maximum in eV for normal functions.
            step: bin interval in eV

        Returns:
            A dict: {"energies": values, "lambda": values, "xas": values}
                    where values are lists of abscissa (energies, lamba) and
                    the sum of gaussian functions (xas).
            A matplotlib plot.
        """
        from pymatgen.util.plotting import pretty_plot
        from matplotlib.mlab import normpdf
        plt = pretty_plot(12, 8)

        transitions = self.read_excitation_energies()

        minval = min([val[0] for val in transitions]) - 5.0 * sigma
        maxval = max([val[0] for val in transitions]) + 5.0 * sigma
        npts = int((maxval - minval) / step) + 1

        eneval = np.linspace(minval, maxval, npts)  # in eV
        lambdaval = [cst.h * cst.c / (val * cst.e) * 1.e9
                     for val in eneval]  # in nm

        # sum of gaussian functions
        spectre = np.zeros(npts)
        for trans in transitions:
            spectre += trans[2] * normpdf(eneval, trans[0], sigma)
        spectre /= spectre.max()
        plt.plot(lambdaval, spectre, "r-", label="spectre")

        data = {"energies": eneval, "lambda": lambdaval, "xas": spectre}

        # plot transitions as vlines
        plt.vlines([val[1] for val in transitions],
                   0.,
                   [val[2] for val in transitions],
                   color="blue",
                   label="transitions",
                   linewidth=2)

        plt.xlabel("$\\lambda$ (nm)")
        plt.ylabel("Arbitrary unit")
        plt.legend()

        return data, plt

    def save_spectre_plot(self, filename="spectre.pdf", img_format="pdf",
                          sigma=0.05, step=0.01):
        """
        Save matplotlib plot of the spectre to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            sigma: Full width at half maximum in eV for normal functions.
            step: bin interval in eV
        """
        d, plt = self.get_spectre_plot(sigma, step)
        plt.savefig(filename, format=img_format)

    def to_input(self, mol=None,  charge=None,
                 spin_multiplicity=None, title=None, functional=None,
                 basis_set=None, route_parameters=None, input_parameters=None,
                 link0_parameters=None, dieze_tag=None, cart_coords=False):
        """
        Create a new input object using by default the last geometry read in
        the output file and with the same calculation parameters. Arguments
        are the same as GaussianInput class.

        Returns
            gaunip (GaussianInput) : the gaussian input object
        """
        if not mol:
            mol = self.final_structure

        if charge is None:
            charge = self.charge

        if spin_multiplicity is None:
            spin_multiplicity = self.spin_multiplicity

        if not title:
            title = self.title

        if not functional:
            functional = self.functional

        if not basis_set:
            basis_set = self.basis_set

        if not route_parameters:
            route_parameters = self.route_parameters

        if not link0_parameters:
            link0_parameters = self.link0

        if not dieze_tag:
            dieze_tag = self.dieze_tag

        return GaussianInput(mol=mol,
                             charge=charge,
                             spin_multiplicity=spin_multiplicity,
                             title=title,
                             functional=functional,
                             basis_set=basis_set,
                             route_parameters=route_parameters,
                             input_parameters=input_parameters,
                             link0_parameters=link0_parameters,
                             dieze_tag=dieze_tag)
