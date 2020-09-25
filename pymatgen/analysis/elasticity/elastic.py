# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module provides a class used to describe the elastic tensor,
including methods used to fit the elastic tensor from linear response
stress-strain data
"""

import itertools
import warnings
from collections import OrderedDict

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import root
from scipy.special import factorial

from monty.dev import deprecated

from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.core.tensors import Tensor, TensorCollection, get_uvec, SquareTensor, DEFAULT_QUAD
from pymatgen.core.units import Unit

__author__ = "Joseph Montoya"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = ("Maarten de Jong, Ian Winter, Shyam Dwaraknath, "
               "Mark Asta, Anubhav Jain")
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Production"
__date__ = "July 24, 2018"


class NthOrderElasticTensor(Tensor):
    """
    An object representing an nth-order tensor expansion
    of the stress-strain constitutive equations
    """
    GPa_to_eV_A3 = Unit("GPa").get_conversion_factor(Unit("eV ang^-3"))
    symbol = "C"

    def __new__(cls, input_array, check_rank=None, tol=1e-4):
        """
        Args:
            input_array ():
            check_rank ():
            tol ():
        """
        obj = super().__new__(
            cls, input_array, check_rank=check_rank)
        if obj.rank % 2 != 0:
            raise ValueError("ElasticTensor must have even rank")
        if not obj.is_voigt_symmetric(tol):
            warnings.warn("Input elastic tensor does not satisfy "
                          "standard voigt symmetries")
        return obj.view(cls)

    @property
    def order(self):
        """
        Order of the elastic tensor
        """
        return self.rank // 2

    def calculate_stress(self, strain):
        """
        Calculate's a given elastic tensor's contribution to the
        stress using Einstein summation

        Args:
            strain (3x3 array-like): matrix corresponding to strain
        """
        strain = np.array(strain)
        if strain.shape == (6,):
            strain = Strain.from_voigt(strain)
        assert strain.shape == (3, 3), "Strain must be 3x3 or voigt-notation"
        stress_matrix = self.einsum_sequence([strain] * (self.order - 1)) / factorial(self.order - 1)
        return Stress(stress_matrix)

    def energy_density(self, strain, convert_GPa_to_eV=True):
        """
        Calculates the elastic energy density due to a strain
        """
        e_density = np.sum(self.calculate_stress(strain) * strain) / self.order
        if convert_GPa_to_eV:
            e_density *= self.GPa_to_eV_A3  # Conversion factor for GPa to eV/A^3
        return e_density

    @classmethod
    def from_diff_fit(cls, strains, stresses, eq_stress=None,
                      order=2, tol=1e-10):
        """

        Args:
            strains ():
            stresses ():
            eq_stress ():
            order ():
            tol ():

        Returns:

        """
        return cls(diff_fit(strains, stresses, eq_stress, order, tol)[order - 2])


def raise_error_if_unphysical(f):
    """
    Wrapper for functions or properties that should raise an error
    if tensor is unphysical.
    """

    def wrapper(self, *args, **kwargs):
        """
        Args:
            self ():
            *args ():
            **kwargs ():

        Returns:

        """
        if self.k_vrh < 0 or self.g_vrh < 0:
            raise ValueError("Bulk or shear modulus is negative, property "
                             "cannot be determined")
        return f(self, *args, **kwargs)

    return wrapper


class ElasticTensor(NthOrderElasticTensor):
    """
    This class extends Tensor to describe the 3x3x3x3
    second-order elastic tensor, C_{ijkl}, with various
    methods for estimating other properties derived from
    the second order elastic tensor
    """

    def __new__(cls, input_array, tol=1e-4):
        """
        Create an ElasticTensor object.  The constructor throws an error if
        the shape of the input_matrix argument is not 3x3x3x3, i. e. in true
        tensor notation.  Issues a warning if the input_matrix argument does
        not satisfy standard symmetries.  Note that the constructor uses
        __new__ rather than __init__ according to the standard method of
        subclassing numpy ndarrays.

        Args:
            input_array (3x3x3x3 array-like): the 3x3x3x3 array-like
                representing the elastic tensor

            tol (float): tolerance for initial symmetry test of tensor
        """

        obj = super().__new__(cls, input_array,
                              check_rank=4, tol=tol)
        return obj.view(cls)

    @property
    def compliance_tensor(self):
        """
        returns the Voigt-notation compliance tensor,
        which is the matrix inverse of the
        Voigt-notation elastic tensor
        """
        s_voigt = np.linalg.inv(self.voigt)
        return ComplianceTensor.from_voigt(s_voigt)

    @property
    def k_voigt(self):
        """
        returns the K_v bulk modulus
        """
        return self.voigt[:3, :3].mean()

    @property
    def g_voigt(self):
        """
        returns the G_v shear modulus
        """
        return (2. * self.voigt[:3, :3].trace() -
                np.triu(self.voigt[:3, :3]).sum() +
                3 * self.voigt[3:, 3:].trace()) / 15.

    @property
    def k_reuss(self):
        """
        returns the K_r bulk modulus
        """
        return 1. / self.compliance_tensor.voigt[:3, :3].sum()

    @property
    def g_reuss(self):
        """
        returns the G_r shear modulus
        """
        return 15. / (8. * self.compliance_tensor.voigt[:3, :3].trace() -
                      4. * np.triu(self.compliance_tensor.voigt[:3, :3]).sum() +
                      3. * self.compliance_tensor.voigt[3:, 3:].trace())

    @property
    def k_vrh(self):
        """
        returns the K_vrh (Voigt-Reuss-Hill) average bulk modulus
        """
        return 0.5 * (self.k_voigt + self.k_reuss)

    @property
    def g_vrh(self):
        """
        returns the G_vrh (Voigt-Reuss-Hill) average shear modulus
        """
        return 0.5 * (self.g_voigt + self.g_reuss)

    @property
    def y_mod(self):
        """
        Calculates Young's modulus (in SI units) using the
        Voigt-Reuss-Hill averages of bulk and shear moduli
        """
        return 9.e9 * self.k_vrh * self.g_vrh / (3. * self.k_vrh + self.g_vrh)

    def directional_poisson_ratio(self, n, m, tol=1e-8):
        """
        Calculates the poisson ratio for a specific direction
        relative to a second, orthogonal direction

        Args:
            n (3-d vector): principal direction
            m (3-d vector): secondary direction orthogonal to n
            tol (float): tolerance for testing of orthogonality
        """
        n, m = get_uvec(n), get_uvec(m)
        if not np.abs(np.dot(n, m)) < tol:
            raise ValueError("n and m must be orthogonal")
        v = self.compliance_tensor.einsum_sequence([n] * 2 + [m] * 2)
        v *= -1 / self.compliance_tensor.einsum_sequence([n] * 4)
        return v

    def directional_elastic_mod(self, n):
        """
        Calculates directional elastic modulus for a specific vector
        """
        n = get_uvec(n)
        return self.einsum_sequence([n] * 4)

    @raise_error_if_unphysical
    def trans_v(self, structure):
        """
        Calculates transverse sound velocity (in SI units) using the
        Voigt-Reuss-Hill average bulk modulus

        Args:
            structure: pymatgen structure object

        Returns: transverse sound velocity (in SI units)

        """
        nsites = structure.num_sites
        volume = structure.volume
        natoms = structure.composition.num_atoms
        weight = float(structure.composition.weight)
        mass_density = 1.6605e3 * nsites * weight / (natoms * volume)
        if self.g_vrh < 0:
            raise ValueError("k_vrh or g_vrh is negative, "
                             "sound velocity is undefined")
        return (1e9 * self.g_vrh / mass_density) ** 0.5

    @raise_error_if_unphysical
    def long_v(self, structure):
        """
        Calculates longitudinal sound velocity (in SI units)
        using the Voigt-Reuss-Hill average bulk modulus

        Args:
            structure: pymatgen structure object

        Returns: longitudinal sound velocity (in SI units)

        """
        nsites = structure.num_sites
        volume = structure.volume
        natoms = structure.composition.num_atoms
        weight = float(structure.composition.weight)
        mass_density = 1.6605e3 * nsites * weight / (natoms * volume)
        if self.g_vrh < 0:
            raise ValueError("k_vrh or g_vrh is negative, "
                             "sound velocity is undefined")
        return (1e9 * (self.k_vrh + 4. / 3. * self.g_vrh) / mass_density) ** 0.5

    @raise_error_if_unphysical
    def snyder_ac(self, structure):
        """
        Calculates Snyder's acoustic sound velocity (in SI units)

        Args:
            structure: pymatgen structure object

        Returns: Snyder's acoustic sound velocity (in SI units)

        """
        nsites = structure.num_sites
        volume = structure.volume
        natoms = structure.composition.num_atoms
        num_density = 1e30 * nsites / volume
        tot_mass = sum([e.atomic_mass for e in structure.species])
        avg_mass = 1.6605e-27 * tot_mass / natoms
        return 0.38483 * avg_mass * ((self.long_v(structure) + 2. * self.trans_v(structure)) / 3.) ** 3. \
            / (300. * num_density ** (-2. / 3.) * nsites ** (1. / 3.))

    @raise_error_if_unphysical
    def snyder_opt(self, structure):
        """
        Calculates Snyder's optical sound velocity (in SI units)

        Args:
            structure: pymatgen structure object

        Returns: Snyder's optical sound velocity (in SI units)

        """
        nsites = structure.num_sites
        volume = structure.volume
        num_density = 1e30 * nsites / volume
        return 1.66914e-23 * (self.long_v(structure) + 2. * self.trans_v(structure)) / 3. \
            / num_density ** (-2. / 3.) * (1 - nsites ** (-1. / 3.))

    @raise_error_if_unphysical
    def snyder_total(self, structure):
        """
        Calculates Snyder's total sound velocity (in SI units)

        Args:
            structure: pymatgen structure object

        Returns: Snyder's total sound velocity (in SI units)

        """
        return self.snyder_ac(structure) + self.snyder_opt(structure)

    @raise_error_if_unphysical
    def clarke_thermalcond(self, structure):
        """
        Calculates Clarke's thermal conductivity (in SI units)

        Args:
            structure: pymatgen structure object

        Returns: Clarke's thermal conductivity (in SI units)

        """
        nsites = structure.num_sites
        volume = structure.volume
        tot_mass = sum([e.atomic_mass for e in structure.species])
        natoms = structure.composition.num_atoms
        weight = float(structure.composition.weight)
        avg_mass = 1.6605e-27 * tot_mass / natoms
        mass_density = 1.6605e3 * nsites * weight / (natoms * volume)
        return 0.87 * 1.3806e-23 * avg_mass ** (-2. / 3.) * mass_density ** (1. / 6.) * self.y_mod ** 0.5

    @raise_error_if_unphysical
    def cahill_thermalcond(self, structure):
        """
        Calculates Cahill's thermal conductivity (in SI units)

        Args:
            structure: pymatgen structure object

        Returns: Cahill's thermal conductivity (in SI units)

        """
        nsites = structure.num_sites
        volume = structure.volume
        num_density = 1e30 * nsites / volume
        return 1.3806e-23 / 2.48 * num_density ** (2. / 3.) * (self.long_v(structure) + 2 * self.trans_v(structure))

    @raise_error_if_unphysical
    def debye_temperature(self, structure):
        """
        Estimates the debye temperature from longitudinal and
        transverse sound velocities

        Args:
            structure: pymatgen structure object

        Returns: debye temperature (in SI units)

        """
        v0 = (structure.volume * 1e-30 / structure.num_sites)
        vl, vt = self.long_v(structure), self.trans_v(structure)
        vm = 3 ** (1. / 3.) * (1 / vl ** 3 + 2 / vt ** 3) ** (-1. / 3.)
        td = 1.05457e-34 / 1.38065e-23 * vm * (6 * np.pi ** 2 / v0) ** (1. / 3.)
        return td

    @deprecated("debye_temperature_from_sound_velocities is now the default"
                "debye_temperature function, this one will be removed.")
    @raise_error_if_unphysical
    def debye_temperature_from_sound_velocities(self, structure):
        """
        Estimates Debye temperature from sound velocities
        """
        return self.debye_temperature(structure)

    @property
    def universal_anisotropy(self):
        """
        returns the universal anisotropy value
        """
        return 5. * self.g_voigt / self.g_reuss + self.k_voigt / self.k_reuss - 6.

    @property
    def homogeneous_poisson(self):
        """
        returns the homogeneous poisson ratio
        """
        return (1. - 2. / 3. * self.g_vrh / self.k_vrh) / (2. + 2. / 3. * self.g_vrh / self.k_vrh)

    def green_kristoffel(self, u):
        """
        Returns the Green-Kristoffel tensor for a second-order tensor
        """
        return self.einsum_sequence([u, u], "ijkl,i,l")

    @property
    def property_dict(self):
        """
        returns a dictionary of properties derived from the elastic tensor
        """
        props = ["k_voigt", "k_reuss", "k_vrh", "g_voigt", "g_reuss", "g_vrh",
                 "universal_anisotropy", "homogeneous_poisson", "y_mod"]
        return {prop: getattr(self, prop) for prop in props}

    def get_structure_property_dict(self, structure, include_base_props=True,
                                    ignore_errors=False):
        """
        returns a dictionary of properties derived from the elastic tensor
        and an associated structure

        Args:
            structure (Structure): structure object for which to calculate
                associated properties
            include_base_props (bool): whether to include base properties,
                like k_vrh, etc.
            ignore_errors (bool): if set to true, will set problem properties
                that depend on a physical tensor to None, defaults to False
        """
        s_props = ["trans_v", "long_v", "snyder_ac", "snyder_opt",
                   "snyder_total", "clarke_thermalcond", "cahill_thermalcond",
                   "debye_temperature"]
        if ignore_errors and (self.k_vrh < 0 or self.g_vrh < 0):
            sp_dict = {prop: None for prop in s_props}
        else:
            sp_dict = {prop: getattr(self, prop)(structure) for prop in s_props}
        sp_dict["structure"] = structure
        if include_base_props:
            sp_dict.update(self.property_dict)
        return sp_dict

    @classmethod
    def from_pseudoinverse(cls, strains, stresses):
        """
        Class method to fit an elastic tensor from stress/strain
        data.  Method uses Moore-Penrose pseudoinverse to invert
        the s = C*e equation with elastic tensor, stress, and
        strain in voigt notation

        Args:
            stresses (Nx3x3 array-like): list or array of stresses
            strains (Nx3x3 array-like): list or array of strains
        """
        # convert the stress/strain to Nx6 arrays of voigt-notation
        warnings.warn("Pseudoinverse fitting of Strain/Stress lists may yield "
                      "questionable results from vasp data, use with caution.")
        stresses = np.array([Stress(stress).voigt for stress in stresses])
        with warnings.catch_warnings(record=True):
            strains = np.array([Strain(strain).voigt for strain in strains])

        voigt_fit = np.transpose(np.dot(np.linalg.pinv(strains), stresses))
        return cls.from_voigt(voigt_fit)

    @classmethod
    def from_independent_strains(cls, strains, stresses, eq_stress=None,
                                 vasp=False, tol=1e-10):
        """
        Constructs the elastic tensor least-squares fit of independent strains
        Args:
            strains (list of Strains): list of strain objects to fit
            stresses (list of Stresses): list of stress objects to use in fit
                corresponding to the list of strains
            eq_stress (Stress): equilibrium stress to use in fitting
            vasp (boolean): flag for whether the stress tensor should be
                converted based on vasp units/convention for stress
            tol (float): tolerance for removing near-zero elements of the
                resulting tensor
        """
        strain_states = [tuple(ss) for ss in np.eye(6)]
        ss_dict = get_strain_state_dict(strains, stresses, eq_stress=eq_stress)
        if not set(strain_states) <= set(ss_dict.keys()):
            raise ValueError("Missing independent strain states: "
                             "{}".format(set(strain_states) - set(ss_dict)))
        if len(set(ss_dict.keys()) - set(strain_states)) > 0:
            warnings.warn("Extra strain states in strain-stress pairs "
                          "are neglected in independent strain fitting")
        c_ij = np.zeros((6, 6))
        for i in range(6):
            istrains = ss_dict[strain_states[i]]["strains"]
            istresses = ss_dict[strain_states[i]]["stresses"]
            for j in range(6):
                c_ij[i, j] = np.polyfit(istrains[:, i], istresses[:, j], 1)[0]
        if vasp:
            c_ij *= -0.1  # Convert units/sign convention of vasp stress tensor
        c = cls.from_voigt(c_ij)
        c = c.zeroed(tol)
        return c


class ComplianceTensor(Tensor):
    """
    This class represents the compliance tensor, and exists
    primarily to keep the voigt-conversion scheme consistent
    since the compliance tensor has a unique vscale
    """

    def __new__(cls, s_array):
        """
        Args:
            s_array ():
        """
        vscale = np.ones((6, 6))
        vscale[3:] *= 2
        vscale[:, 3:] *= 2
        obj = super().__new__(cls, s_array, vscale=vscale)
        return obj.view(cls)


class ElasticTensorExpansion(TensorCollection):
    """
    This class is a sequence of elastic tensors corresponding
    to an elastic tensor expansion, which can be used to
    calculate stress and energy density and inherits all
    of the list-based properties of TensorCollection
    (e. g. symmetrization, voigt conversion, etc.)
    """

    def __init__(self, c_list):
        """
        Initialization method for ElasticTensorExpansion

        Args:
            c_list (list or tuple): sequence of Tensor inputs
                or tensors from which the elastic tensor
                expansion is constructed.
        """
        c_list = [NthOrderElasticTensor(c, check_rank=4 + i * 2)
                  for i, c in enumerate(c_list)]
        super().__init__(c_list)

    @classmethod
    def from_diff_fit(cls, strains, stresses, eq_stress=None,
                      tol=1e-10, order=3):
        """
        Generates an elastic tensor expansion via the fitting function
        defined below in diff_fit
        """
        c_list = diff_fit(strains, stresses, eq_stress, order, tol)
        return cls(c_list)

    @property
    def order(self):
        """
        Order of the elastic tensor expansion, i. e. the order of the
        highest included set of elastic constants
        """
        return self[-1].order

    def calculate_stress(self, strain):
        """
        Calculate's a given elastic tensor's contribution to the
        stress using Einstein summation
        """
        return sum([c.calculate_stress(strain) for c in self])

    def energy_density(self, strain, convert_GPa_to_eV=True):
        """
        Calculates the elastic energy density due to a strain
        """
        return sum([c.energy_density(strain, convert_GPa_to_eV)
                    for c in self])

    def get_ggt(self, n, u):
        """
        Gets the Generalized Gruneisen tensor for a given
        third-order elastic tensor expansion.

        Args:
            n (3x1 array-like): normal mode direction
            u (3x1 array-like): polarization direction
        """
        gk = self[0].einsum_sequence([n, u, n, u])
        result = -(2 * gk * np.outer(u, u) + self[0].einsum_sequence([n, n])
                   + self[1].einsum_sequence([n, u, n, u])) / (2 * gk)
        return result

    def get_tgt(self, temperature=None, structure=None, quad=None):
        """
        Gets the thermodynamic Gruneisen tensor (TGT) by via an
        integration of the GGT weighted by the directional heat
        capacity.

        See refs:
            R. N. Thurston and K. Brugger, Phys. Rev. 113, A1604 (1964).
            K. Brugger Phys. Rev. 137, A1826 (1965).

        Args:
            temperature (float): Temperature in kelvin, if not specified
                will return non-cv-normalized value
            structure (float): Structure to be used in directional heat
                capacity determination, only necessary if temperature
                is specified
            quad (dict): quadrature for integration, should be
                dictionary with "points" and "weights" keys defaults
                to quadpy.sphere.Lebedev(19) as read from file
        """
        if temperature and not structure:
            raise ValueError("If using temperature input, you must also "
                             "include structure")

        quad = quad if quad else DEFAULT_QUAD
        points = quad['points']
        weights = quad['weights']
        num, denom, c = np.zeros((3, 3)), 0, 1
        for p, w in zip(points, weights):
            gk = ElasticTensor(self[0]).green_kristoffel(p)
            rho_wsquareds, us = np.linalg.eigh(gk)
            us = [u / np.linalg.norm(u) for u in np.transpose(us)]
            for u in us:
                # TODO: this should be benchmarked
                if temperature:
                    c = self.get_heat_capacity(temperature, structure, p, u)
                num += c * self.get_ggt(p, u) * w
                denom += c * w
        return SquareTensor(num / denom)

    def get_gruneisen_parameter(self, temperature=None, structure=None,
                                quad=None):
        """
        Gets the single average gruneisen parameter from the TGT.

        Args:
            temperature (float): Temperature in kelvin, if not specified
                will return non-cv-normalized value
            structure (float): Structure to be used in directional heat
                capacity determination, only necessary if temperature
                is specified
            quad (dict): quadrature for integration, should be
                dictionary with "points" and "weights" keys defaults
                to quadpy.sphere.Lebedev(19) as read from file
        """
        return np.trace(self.get_tgt(temperature, structure, quad)) / 3.

    def get_heat_capacity(self, temperature, structure, n, u, cutoff=1e2):
        """
        Gets the directional heat capacity for a higher order tensor
        expansion as a function of direction and polarization.

        Args:
            temperature (float): Temperature in kelvin
            structure (float): Structure to be used in directional heat
                capacity determination
            n (3x1 array-like): direction for Cv determination
            u (3x1 array-like): polarization direction, note that
                no attempt for verification of eigenvectors is made
            cutoff (float): cutoff for scale of kt / (hbar * omega)
                if lower than this value, returns 0
        """
        k = 1.38065e-23
        kt = k * temperature
        hbar_w = 1.05457e-34 * self.omega(structure, n, u)
        if hbar_w > kt * cutoff:
            return 0.0
        c = k * (hbar_w / kt) ** 2
        c *= np.exp(hbar_w / kt) / (np.exp(hbar_w / kt) - 1) ** 2
        return c * 6.022e23

    def omega(self, structure, n, u):
        """
        Finds directional frequency contribution to the heat
        capacity from direction and polarization

        Args:
            structure (Structure): Structure to be used in directional heat
                capacity determination
            n (3x1 array-like): direction for Cv determination
            u (3x1 array-like): polarization direction, note that
                no attempt for verification of eigenvectors is made
        """
        l0 = np.dot(np.sum(structure.lattice.matrix, axis=0), n)
        l0 *= 1e-10  # in A
        weight = float(structure.composition.weight) * 1.66054e-27  # in kg
        vol = structure.volume * 1e-30  # in m^3
        vel = (1e9 * self[0].einsum_sequence([n, u, n, u])
               / (weight / vol)) ** 0.5
        return vel / l0

    def thermal_expansion_coeff(self, structure, temperature, mode="debye"):
        """
        Gets thermal expansion coefficient from third-order constants.

        Args:
            temperature (float): Temperature in kelvin, if not specified
                will return non-cv-normalized value
            structure (Structure): Structure to be used in directional heat
                capacity determination, only necessary if temperature
                is specified
            mode (string): mode for finding average heat-capacity,
                current supported modes are 'debye' and 'dulong-petit'
        """
        soec = ElasticTensor(self[0])
        v0 = (structure.volume * 1e-30 / structure.num_sites)
        if mode == "debye":
            td = soec.debye_temperature(structure)
            t_ratio = temperature / td

            def integrand(x):
                return (x ** 4 * np.exp(x)) / (np.exp(x) - 1) ** 2
            cv = 9 * 8.314 * t_ratio ** 3 * quad(integrand, 0, t_ratio ** -1)[0]
        elif mode == "dulong-petit":
            cv = 3 * 8.314
        else:
            raise ValueError("Mode must be debye or dulong-petit")
        tgt = self.get_tgt(temperature, structure)
        alpha = np.einsum('ijkl,ij', soec.compliance_tensor, tgt)
        alpha *= cv / (1e9 * v0 * 6.022e23)
        return SquareTensor(alpha)

    def get_compliance_expansion(self):
        """
        Gets a compliance tensor expansion from the elastic
        tensor expansion.
        """
        # TODO: this might have a general form
        if not self.order <= 4:
            raise ValueError("Compliance tensor expansion only "
                             "supported for fourth-order and lower")
        ce_exp = [ElasticTensor(self[0]).compliance_tensor]
        einstring = "ijpq,pqrsuv,rskl,uvmn->ijklmn"
        ce_exp.append(np.einsum(einstring, -ce_exp[-1], self[1],
                                ce_exp[-1], ce_exp[-1]))
        if self.order == 4:
            # Four terms in the Fourth-Order compliance tensor
            # pylint: disable=E1130
            einstring_1 = "pqab,cdij,efkl,ghmn,abcdefgh"
            tensors_1 = [ce_exp[0]] * 4 + [self[-1]]
            temp = -np.einsum(einstring_1, *tensors_1)
            einstring_2 = "pqab,abcdef,cdijmn,efkl"
            einstring_3 = "pqab,abcdef,efklmn,cdij"
            einstring_4 = "pqab,abcdef,cdijkl,efmn"
            for es in [einstring_2, einstring_3, einstring_4]:
                temp -= np.einsum(es, ce_exp[0], self[-2], ce_exp[1], ce_exp[0])
            ce_exp.append(temp)
        return TensorCollection(ce_exp)

    def get_strain_from_stress(self, stress):
        """
        Gets the strain from a stress state according
        to the compliance expansion corresponding to the
        tensor expansion.
        """
        compl_exp = self.get_compliance_expansion()
        strain = 0
        for n, compl in enumerate(compl_exp):
            strain += compl.einsum_sequence([stress] * (n + 1)) / factorial(n + 1)
        return strain

    def get_effective_ecs(self, strain, order=2):
        """
        Returns the effective elastic constants
        from the elastic tensor expansion.

        Args:
            strain (Strain or 3x3 array-like): strain condition
                under which to calculate the effective constants
            order (int): order of the ecs to be returned
        """
        ec_sum = 0
        for n, ecs in enumerate(self[order - 2:]):
            ec_sum += ecs.einsum_sequence([strain] * n) / factorial(n)
        return ec_sum

    def get_wallace_tensor(self, tau):
        """
        Gets the Wallace Tensor for determining yield strength
        criteria.

        Args:
            tau (3x3 array-like): stress at which to evaluate
                the wallace tensor
        """
        b = 0.5 * (np.einsum("ml,kn->klmn", tau, np.eye(3)) +
                   np.einsum("km,ln->klmn", tau, np.eye(3)) +
                   np.einsum("nl,km->klmn", tau, np.eye(3)) +
                   np.einsum("kn,lm->klmn", tau, np.eye(3)) +
                   -2 * np.einsum("kl,mn->klmn", tau, np.eye(3)))
        strain = self.get_strain_from_stress(tau)
        b += self.get_effective_ecs(strain)
        return b

    def get_symmetric_wallace_tensor(self, tau):
        """
        Gets the symmetrized wallace tensor for determining
        yield strength criteria.

        Args:
            tau (3x3 array-like): stress at which to evaluate
                the wallace tensor.
        """
        wallace = self.get_wallace_tensor(tau)
        return Tensor(0.5 * (wallace + np.transpose(wallace, [2, 3, 0, 1])))

    def get_stability_criteria(self, s, n):
        """
        Gets the stability criteria from the symmetric
        Wallace tensor from an input vector and stress
        value.

        Args:
            s (float): Stress value at which to evaluate
                the stability criteria
            n (3x1 array-like): direction of the applied
                stress
        """
        n = get_uvec(n)
        stress = s * np.outer(n, n)
        sym_wallace = self.get_symmetric_wallace_tensor(stress)
        return np.linalg.det(sym_wallace.voigt)

    def get_yield_stress(self, n):
        """
        Gets the yield stress for a given direction

        Args:
            n (3x1 array-like): direction for which to find the
                yield stress
        """
        # TODO: root finding could be more robust
        comp = root(self.get_stability_criteria, -1, args=n)
        tens = root(self.get_stability_criteria, 1, args=n)
        return (comp.x, tens.x)


# TODO: abstract this for other tensor fitting procedures
def diff_fit(strains, stresses, eq_stress=None, order=2, tol=1e-10):
    """
    nth order elastic constant fitting function based on
    central-difference derivatives with respect to distinct
    strain states.  The algorithm is summarized as follows:

    1. Identify distinct strain states as sets of indices
       for which nonzero strain values exist, typically
       [(0), (1), (2), (3), (4), (5), (0, 1) etc.]
    2. For each strain state, find and sort strains and
       stresses by strain value.
    3. Find first, second .. nth derivatives of each stress
       with respect to scalar variable corresponding to
       the smallest perturbation in the strain.
    4. Use the pseudoinverse of a matrix-vector expression
       corresponding to the parameterized stress-strain
       relationship and multiply that matrix by the respective
       calculated first or second derivatives from the
       previous step.
    5. Place the calculated nth-order elastic
       constants appropriately.

    Args:
        order (int): order of the elastic tensor set to return
        strains (nx3x3 array-like): Array of 3x3 strains
            to use in fitting of ECs
        stresses (nx3x3 array-like): Array of 3x3 stresses
            to use in fitting ECs.  These should be PK2 stresses.
        eq_stress (3x3 array-like): stress corresponding to
            equilibrium strain (i. e. "0" strain state).
            If not specified, function will try to find
            the state in the list of provided stresses
            and strains.  If not found, defaults to 0.
        tol (float): value for which strains below
            are ignored in identifying strain states.

    Returns:
        Set of tensors corresponding to nth order expansion of
        the stress/strain relation
    """
    strain_state_dict = get_strain_state_dict(
        strains, stresses, eq_stress=eq_stress, tol=tol,
        add_eq=True, sort=True)

    # Collect derivative data
    c_list = []
    dei_dsi = np.zeros((order - 1, 6, len(strain_state_dict)))
    for n, (strain_state, data) in enumerate(strain_state_dict.items()):
        hvec = data["strains"][:, strain_state.index(1)]
        for i in range(1, order):
            coef = get_diff_coeff(hvec, i)
            dei_dsi[i - 1, :, n] = np.dot(coef, data["stresses"])

    m, absent = generate_pseudo(list(strain_state_dict.keys()), order)
    for i in range(1, order):
        cvec, carr = get_symbol_list(i + 1)
        svec = np.ravel(dei_dsi[i - 1].T)
        cmap = dict(zip(cvec, np.dot(m[i - 1], svec)))
        c_list.append(v_subs(carr, cmap))
    return [Tensor.from_voigt(c) for c in c_list]


def find_eq_stress(strains, stresses, tol=1e-10):
    """
    Finds stress corresponding to zero strain state in stress-strain list

    Args:
        strains (Nx3x3 array-like): array corresponding to strains
        stresses (Nx3x3 array-like): array corresponding to stresses
        tol (float): tolerance to find zero strain state
    """
    stress_array = np.array(stresses)
    strain_array = np.array(strains)
    eq_stress = stress_array[np.all(abs(strain_array) < tol, axis=(1, 2))]

    if eq_stress.size != 0:
        all_same = (abs(eq_stress - eq_stress[0]) < 1e-8).all()
        if len(eq_stress) > 1 and not all_same:
            raise ValueError("Multiple stresses found for equilibrium strain"
                             " state, please specify equilibrium stress or  "
                             " remove extraneous stresses.")
        eq_stress = eq_stress[0]
    else:
        warnings.warn("No eq state found, returning zero voigt stress")
        eq_stress = Stress(np.zeros((3, 3)))
    return eq_stress


def get_strain_state_dict(strains, stresses, eq_stress=None,
                          tol=1e-10, add_eq=True, sort=True):
    """
    Creates a dictionary of voigt-notation stress-strain sets
    keyed by "strain state", i. e. a tuple corresponding to
    the non-zero entries in ratios to the lowest nonzero value,
    e.g. [0, 0.1, 0, 0.2, 0, 0] -> (0,1,0,2,0,0)
    This allows strains to be collected in stencils as to
    evaluate parameterized finite difference derivatives

    Args:
        strains (Nx3x3 array-like): strain matrices
        stresses (Nx3x3 array-like): stress matrices
        eq_stress (Nx3x3 array-like): equilibrium stress
        tol (float): tolerance for sorting strain states
        add_eq (bool): flag for whether to add eq_strain
            to stress-strain sets for each strain state
        sort (bool): flag for whether to sort strain states

    Returns:
        OrderedDict with strain state keys and dictionaries
        with stress-strain data corresponding to strain state
    """
    # Recast stress/strains
    vstrains = np.array([Strain(s).zeroed(tol).voigt for s in strains])
    vstresses = np.array([Stress(s).zeroed(tol).voigt for s in stresses])
    # Collect independent strain states:
    independent = {tuple(np.nonzero(vstrain)[0].tolist()) for vstrain in vstrains}
    strain_state_dict = OrderedDict()
    if add_eq:
        if eq_stress is not None:
            veq_stress = Stress(eq_stress).voigt
        else:
            veq_stress = find_eq_stress(strains, stresses).voigt

    for n, ind in enumerate(independent):
        # match strains with templates
        template = np.zeros(6, dtype=bool)
        np.put(template, ind, True)
        template = np.tile(template, [vstresses.shape[0], 1])
        mode = (template == (np.abs(vstrains) > 1e-10)).all(axis=1)
        mstresses = vstresses[mode]
        mstrains = vstrains[mode]
        # Get "strain state", i.e. ratio of each value to minimum strain
        min_nonzero_ind = np.argmin(np.abs(np.take(mstrains[-1], ind)))
        min_nonzero_val = np.take(mstrains[-1], ind)[min_nonzero_ind]
        strain_state = mstrains[-1] / min_nonzero_val
        strain_state = tuple(strain_state)

        if add_eq:
            # add zero strain state
            mstrains = np.vstack([mstrains, np.zeros(6)])
            mstresses = np.vstack([mstresses, veq_stress])
        # sort strains/stresses by strain values
        if sort:
            mstresses = mstresses[mstrains[:, ind[0]].argsort()]
            mstrains = mstrains[mstrains[:, ind[0]].argsort()]
        strain_state_dict[strain_state] = {"strains": mstrains,
                                           "stresses": mstresses}
    return strain_state_dict


def generate_pseudo(strain_states, order=3):
    """
    Generates the pseudoinverse for a given set of strains.

    Args:
        strain_states (6xN array like): a list of voigt-notation
            "strain-states", i. e. perturbed indices of the strain
            as a function of the smallest strain e. g. (0, 1, 0, 0, 1, 0)
        order (int): order of pseudoinverse to calculate

    Returns:
        mis: pseudo inverses for each order tensor, these can
            be multiplied by the central difference derivative
            of the stress with respect to the strain state
        absent_syms: symbols of the tensor absent from the PI
            expression
    """
    s = sp.Symbol('s')
    nstates = len(strain_states)
    ni = np.array(strain_states) * s
    mis, absent_syms = [], []
    for degree in range(2, order + 1):
        cvec, carr = get_symbol_list(degree)
        sarr = np.zeros((nstates, 6), dtype=object)
        for n, strain_v in enumerate(ni):
            # Get expressions
            exps = carr.copy()
            for i in range(degree - 1):
                exps = np.dot(exps, strain_v)
            exps /= np.math.factorial(degree - 1)
            sarr[n] = [sp.diff(exp, s, degree - 1) for exp in exps]
        svec = sarr.ravel()
        present_syms = set.union(*[exp.atoms(sp.Symbol) for exp in svec])
        absent_syms += [set(cvec) - present_syms]
        m = np.zeros((6 * nstates, len(cvec)))
        for n, c in enumerate(cvec):
            m[:, n] = v_diff(svec, c)
        mis.append(np.linalg.pinv(m))
    return mis, absent_syms


def get_symbol_list(rank, dim=6):
    """
    Returns a symbolic representation of the voigt-notation
    tensor that places identical symbols for entries related
    by index transposition, i. e. C_1121 = C_1211 etc.

    Args:
        dim (int): dimension of matrix/tensor, e. g. 6 for
            voigt notation and 3 for standard
        rank (int): rank of tensor, e. g. 3 for third-order ECs

    Returns:
        c_vec (array): array representing distinct indices
        c_arr (array): array representing tensor with equivalent
            indices assigned as above
    """
    indices = list(
        itertools.combinations_with_replacement(range(dim), r=rank))
    c_vec = np.zeros(len(indices), dtype=object)
    c_arr = np.zeros([dim] * rank, dtype=object)
    for n, idx in enumerate(indices):
        c_vec[n] = sp.Symbol('c_' + ''.join([str(i) for i in idx]))
        for perm in itertools.permutations(idx):
            c_arr[perm] = c_vec[n]
    return c_vec, c_arr


def subs(entry, cmap):
    """
    Sympy substitution function, primarily for the purposes
    of numpy vectorization

    Args:
        entry (symbol or exp): sympy expr to undergo subs
        cmap (dict): map for symbols to values to use in subs

    Returns:
        Evaluated expression with substitution
    """
    return entry.subs(cmap)


# Vectorized functions
v_subs = np.vectorize(subs)
v_diff = np.vectorize(sp.diff)


def get_diff_coeff(hvec, n=1):
    """
    Helper function to find difference coefficients of an
    derivative on an arbitrary mesh.

    Args:
        hvec (1D array-like): sampling stencil
        n (int): degree of derivative to find
    """
    hvec = np.array(hvec, dtype=np.float)
    acc = len(hvec)
    exp = np.column_stack([np.arange(acc)] * acc)
    a = np.vstack([hvec] * acc) ** exp
    b = np.zeros(acc)
    b[n] = factorial(n)
    return np.linalg.solve(a, b)
