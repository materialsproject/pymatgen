# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals
from __future__ import absolute_import

from pymatgen.analysis.elasticity.tensors import Tensor, \
        voigt_map as vmap, TensorCollection
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity.strain import Strain
from scipy.misc import factorial
from collections import OrderedDict
import numpy as np
import warnings
import itertools
import string

import sympy as sp

"""
This module provides a class used to describe the elastic tensor,
including methods used to fit the elastic tensor from linear response
stress-strain data
"""


__author__ = "Maarten de Jong, Joseph Montoya"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = ("Ian Winter, Shyam Dwaraknath, "
               "Mark Asta, Anubhav Jain")
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Development"
__date__ = "March 22, 2012"


class NthOrderElasticTensor(Tensor):
    """
    An object representing an nth-order tensor expansion 
    of the stress-strain constitutive equations
    """
    def __new__(cls, input_array, check_rank=None, tol=1e-4):
        obj = super(NthOrderElasticTensor, cls).__new__(
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
        lc = string.ascii_lowercase[:self.rank-2]
        lc_pairs = map(''.join, zip(*[iter(lc)]*2))
        einsum_string = "ij" + lc + ',' + ','.join(lc_pairs) + "->ij"
        einsum_args = [self] + [strain] * (self.order - 1)
        stress_matrix = np.einsum(einsum_string, *einsum_args) \
            / factorial(self.order - 1)
        return Stress(stress_matrix)

    def energy_density(self, strain, convert_GPa_to_eV=True):
        """
        Calculates the elastic energy density due to a strain
        """
        e_density = np.sum(self.calculate_stress(strain)*strain) / self.order
        if convert_GPa_to_eV:
            e_density *= 0.000624151  # Conversion factor for GPa to eV/A^3
        return e_density

    @classmethod
    def from_diff_fit(cls, strains, stresses, eq_stress=None,
                      order=2, tol=1e-10):
        return cls(diff_fit(strains, stresses, eq_stress, order, tol)[order-2])


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

        obj = super(ElasticTensor, cls).__new__(cls, input_array,
                                                check_rank=4, tol=tol)
        return obj.view(cls)

    @property
    def compliance_tensor(self):
        """
        returns the Voigt-notation compliance tensor, 
        which is the matrix inverse of the
        Voigt-notation elastic tensor
        """
        return np.linalg.inv(self.voigt)

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
        return 1. / self.compliance_tensor[:3, :3].sum()

    @property
    def g_reuss(self):
        """
        returns the G_r shear modulus
        """
        return 15. / (8. * self.compliance_tensor[:3, :3].trace() -
                      4. * np.triu(self.compliance_tensor[:3, :3]).sum() +
                      3. * self.compliance_tensor[3:, 3:].trace())

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
        Calculates Young's modulus (in SI units) using the Voigt-Reuss-Hill
        averages of bulk and shear moduli
        """
        return 9.e9 * self.k_vrh * self.g_vrh / (3. * self.k_vrh + self.g_vrh)

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
        weight = structure.composition.weight
        mass_density = 1.6605e3 * nsites * weight / (natoms * volume)
        return (1e9 * self.g_vrh / mass_density) ** 0.5

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
        weight = structure.composition.weight
        mass_density = 1.6605e3 * nsites * weight / (natoms * volume)
        return (1e9 * (self.k_vrh + 4./3. * self.g_vrh) / mass_density) ** 0.5

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
        return 0.38483*avg_mass * \
            ((self.long_v(structure) + 2.*self.trans_v(structure))/3.) ** 3.\
            / (300.*num_density ** (-2./3.) * nsites ** (1./3.))

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
        return 1.66914e-23 * \
            (self.long_v(structure) + 2.*self.trans_v(structure))/3. \
            / num_density ** (-2./3.) * (1 - nsites ** (-1./3.))

    def snyder_total(self, structure):
        """
        Calculates Snyder's total sound velocity (in SI units)

        Args:
            structure: pymatgen structure object

        Returns: Snyder's total sound velocity (in SI units)

        """
        return self.snyder_ac(structure) + self.snyder_opt(structure)

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
        weight = structure.composition.weight
        avg_mass = 1.6605e-27 * tot_mass / natoms
        mass_density = 1.6605e3 * nsites * weight / (natoms * volume)
        return 0.87 * 1.3806e-23 * avg_mass**(-2./3.) \
            * mass_density**(1./6.) * self.y_mod**0.5

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
        return 1.3806e-23 / 2.48 * num_density**(2./3.) \
            * (self.long_v(structure) + 2 * self.trans_v(structure))

    def debye_temperature(self, structure):
        """
        Calculates the debye temperature (in SI units)

        Args:
            structure: pymatgen structure object

        Returns: debye temperature (in SI units)

        """
        nsites = structure.num_sites
        volume = structure.volume
        tot_mass = sum([e.atomic_mass for e in structure.species])
        natoms = structure.composition.num_atoms
        weight = structure.composition.weight
        avg_mass = 1.6605e-27 * tot_mass / natoms
        mass_density = 1.6605e3 * nsites * weight / (natoms * volume)
        return 2.589e-11 * avg_mass**(-1./3.) * mass_density**(-1./6.) \
            * self.y_mod**0.5

    def debye_temperature_gibbs(self, structure):
        """
        Calculates the debye temperature accordings to the GIBBS
        formulation (in SI units)

        Args:
            structure: pymatgen structure object

        Returns: debye temperature (in SI units)

        """
        volume = structure.volume
        tot_mass = sum([e.atomic_mass for e in structure.species])
        natoms = structure.composition.num_atoms
        avg_mass = 1.6605e-27 * tot_mass / natoms
        t = self.homogeneous_poisson
        f = (3.*(2.*(2./3.*(1. + t)/(1. - 2.*t))**1.5 +
                 (1./3.*(1. + t)/(1. - t))**1.5)**-1) ** (1./3.)
        return 2.9772e-11 * avg_mass**(-1./2.) * (volume / natoms) ** (-1./6.) \
            * f * self.k_vrh ** 0.5

    @property
    def universal_anisotropy(self):
        """
        returns the universal anisotropy value
        """
        return 5. * self.g_voigt / self.g_reuss + \
            self.k_voigt / self.k_reuss - 6.

    @property
    def homogeneous_poisson(self):
        """
        returns the homogeneous poisson ratio
        """
        return (1. - 2. / 3. * self.g_vrh / self.k_vrh) / \
               (2. + 2. / 3. * self.g_vrh / self.k_vrh)

    @property
    def property_dict(self):
        """
        returns a dictionary of properties derived from the elastic tensor
        """
        props = ["k_voigt", "k_reuss", "k_vrh", "g_voigt", "g_reuss", "g_vrh",
                 "universal_anisotropy", "homogeneous_poisson", "y_mod"]
        return {prop: getattr(self, prop) for prop in props}

    def get_structure_property_dict(self, structure, include_base_props=True):
        """
        returns a dictionary of properties derived from the elastic tensor
        and an associated structure
        """
        s_props = ["trans_v", "long_v", "snyder_ac", "snyder_opt",
                   "snyder_total", "clarke_thermalcond", "cahill_thermalcond",
                   "debye_temperature", "debye_temperature_gibbs"]
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
    def from_stress_dict(cls, stress_dict, vasp=True):
        """
        Constructs the elastic tensor from IndependentStrain-Stress dictionary
        corresponding to legacy behavior of elasticity package.

        Args:
            stress_dict (dict): dictionary of stresses indexed by corresponding
                IndependentStrain objects.
            vasp (boolean): flag for whether the stress tensor should be
                converted based on vasp units/convention for stress
        """
        c_ij = np.zeros((6, 6))
        for i, j in itertools.product(range(6), repeat=2):
            strains = [s for s in stress_dict.keys()
                       if s.ij == vmap[i]]
            xy = [(s[vmap[i]], stress_dict[s][vmap[j]]) for s in strains]
            if len(xy) == 0:
                raise ValueError("No ind. strains for vgt index {}".format(i))
            elif len(xy) == 1:
                xy += [(0, 0)]  # Fit through 0
            c_ij[i, j] = np.polyfit(*zip(*xy), deg=1)[0]
        if vasp:
            c_ij *= -0.1  # Convert units/sign convention of vasp stress tensor
        c_ij[0:, 3:] = 0.5 * c_ij[0:, 3:]  # for vgt doubling of e4,e5,e6
        c = cls.from_voigt(c_ij)
        c = c.zeroed()
        return c


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
        c_list = [NthOrderElasticTensor(c, check_rank=4+i*2)
                  for i, c in enumerate(c_list)]
        super(ElasticTensorExpansion, self).__init__(c_list)

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
            dei_dsi[i-1, :, n] = np.dot(coef, data["stresses"])

    m, absent = generate_pseudo(list(strain_state_dict.keys()), order)
    for i in range(1, order):
        cvec, carr = get_symbol_list(i+1)
        svec = np.ravel(dei_dsi[i-1].T)
        cmap = dict(zip(cvec, np.dot(m[i-1], svec)))
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
    eq_stress = stress_array[np.all(abs(strain_array)<tol, axis=(1,2))]

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
    independent = set([tuple(np.nonzero(vstrain)[0].tolist())
                       for vstrain in vstrains])
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

        if add_eq:
            # add zero strain state
            mstrains = np.vstack([mstrains, np.zeros(6)])
            mstresses = np.vstack([mstresses, veq_stress])
        # sort strains/stresses by strain values
        if sort:
            mstresses = mstresses[mstrains[:, ind[0]].argsort()]
            mstrains = mstrains[mstrains[:, ind[0]].argsort()]
        # Get "strain state", i.e. ratio of each value to minimum strain
        strain_state = mstrains[-1] / np.min(np.take(mstrains[-1], ind))
        strain_state = tuple(strain_state)
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
    ni = np.array(strain_states)*s
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
        m = np.zeros((6*nstates, len(cvec)))
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
    c_arr = np.zeros([dim]*rank, dtype=object)
    for n, idx in enumerate(indices):
        c_vec[n] = sp.Symbol('c_'+''.join([str(i) for i in idx]))
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
    exp = np.column_stack([np.arange(acc)]*acc)
    a = np.vstack([hvec] * acc) ** exp
    b = np.zeros(acc)
    b[n] = factorial(n)
    return np.linalg.solve(a, b)
