# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals
from __future__ import absolute_import

"""
This module provides a class used to describe the elastic tensor,
including methods used to fit the elastic tensor from linear response
stress-strain data
"""

from pymatgen.analysis.elasticity.tensors import TensorBase, \
        voigt_map as vmap
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity.strain import Strain
import numpy as np
import warnings
import itertools
from six.moves import range
from monty.dev import requires

__author__ = "Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = ("Joseph Montoya, Ian Winter, Shyam Dwaraknath, "
               "Mark Asta, Anubhav Jain")
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Development"
__date__ = "March 22, 2012"

try:
    import sympy as sp
    sympy_found = True
except ImportError:
    sympy_found = False

class ElasticTensor(TensorBase):
    """
    This class extends TensorBase to describe the 3x3x3x3
    elastic tensor, C_{ij}, in Voigt-notation
    """

    def __new__(cls, input_array, tol=1e-3):
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

        obj = TensorBase(input_array).view(cls)
        if obj.shape != (3, 3, 3, 3):
            raise ValueError("Default elastic tensor constructor requires "
                             "input to be the true 3x3x3x3 representation. "
                             "To construct from an elastic tensor from "
                             "6x6 Voigt array, use ElasticTensor.from_voigt")

        if not ((obj - np.transpose(obj, (1, 0, 2, 3)) < tol).all() and
                    (obj - np.transpose(obj, (0, 1, 3, 2)) < tol).all() and
                    (obj - np.transpose(obj, (1, 0, 3, 2)) < tol).all() and
                    (obj - np.transpose(obj, (3, 2, 0, 1)) < tol).all()):
            warnings.warn("Input elasticity tensor does "
                          "not satisfy standard symmetries")

        return obj


    @property
    def compliance_tensor(self):
        """
        returns the compliance tensor, which is the matrix inverse of the
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
    def kg_average(self):
        """
        returns a list of Voigt, Reuss, and Voigt-Reuss-Hill averages of bulk
        and shear moduli similar to legacy behavior
        """
        return [self.k_voigt, self.g_voigt, self.k_reuss, self.g_reuss,
                self.k_vrh, self.g_vrh]

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
        nsites = structure.num_sites
        volume = structure.volume
        tot_mass = sum([e.atomic_mass for e in structure.species])
        natoms = structure.composition.num_atoms
        avg_mass = 1.6605e-27 * tot_mass / natoms
        t = self.homogeneous_poisson
        f = (3.*(2.*(2./3.*(1. + t)/(1. - 2.*t))**(1.5) + \
                 (1./3.*(1. + t)/(1. - t))**(1.5))**-1) ** (1./3.)
        return 2.9772e-11 * avg_mass**(-1./2.) * (volume / natoms) ** (-1./6.) \
                * f * self.k_vrh**(0.5)

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

    def energy_density(self, strain):
        """
        Calculates the elastic energy density due to a strain
        """
        # Conversion factor for GPa to eV/Angstrom^3
        GPA_EV = 0.000624151

        with warnings.catch_warnings(record=True):
            e_density = np.dot(np.transpose(Strain(strain).voigt),
                               np.dot(self.voigt, Strain(strain).voigt)) / 2 * GPA_EV

        return e_density

    @classmethod
    def from_strain_stress_list(cls, strains, stresses):
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
        warnings.warn("Linear fitting of Strain/Stress lists may yield "
                      "questionable results from vasp data, use with caution.")
        stresses = np.array([Stress(stress).voigt for stress in stresses])
        with warnings.catch_warnings(record=True):
            strains = np.array([Strain(strain).voigt for strain in strains])

        voigt_fit = np.transpose(np.dot(np.linalg.pinv(strains), stresses))
        return cls.from_voigt(voigt_fit)

    @classmethod
    def from_stress_dict(cls, stress_dict, tol=0.1, vasp=True, symmetry=False):
        """
        Constructs the elastic tensor from IndependentStrain-Stress dictionary
        corresponding to legacy behavior of elasticity package.

        Args:
            stress_dict (dict): dictionary of stresses indexed by corresponding
                IndependentStrain objects.
            tol (float): tolerance for zeroing small values of the tensor
            vasp (boolean): flag for whether the stress tensor should be
                converted based on vasp units/convention for stress
            symmetry (boolean): flag for whether or not the elastic tensor
                should fit from data based on symmetry
        """
        c_ij = np.zeros((6, 6))
        for i, j in itertools.product(range(6), repeat=2):
            strains = [s for s in stress_dict.keys() 
                       if (s.i, s.j) == vmap[i]]
            xy = [(s[vmap[i]], stress_dict[s][vmap[j]]) for s in strains]
            if len(xy) == 0:
                raise ValueError("No ind. strains for vgt index {}".format(i))
            elif len(xy) == 1:
                xy += [(0, 0)]
            c_ij[i, j] = np.polyfit(*zip(*xy), deg=1)[0]
        if vasp:
            c_ij *= -0.1  # Convert units/sign convention of vasp stress tensor
        c_ij[0:, 3:] = 0.5 * c_ij[0:, 3:]  # account for voigt doubling of e4,e5,e6
        c = cls.from_voigt(c_ij)
        c = c.zeroed()
        return c

    @property
    def voigt_symmetrized(self):
        """
        Reconstructs the elastic tensor by symmetrizing the voigt
        notation tensor, to allow for legacy behavior
        """

        v = self.voigt
        new_v = 0.5 * (np.transpose(v) + v)
        return ElasticTensor.from_voigt(new_v)

@requires(sympy_found, "TOEC fitter requires sympy")
def toec_fit(strains, stresses, eq_stress = None, zero_crit=1e-10):
    """
    A third-order elastic constant fitting function based on 
    central-difference derivatives with respect to distinct
    strain states.  The algorithm is summarized as follows:

    1. Identify distinct strain states as sets of indices 
       for which nonzero strain values exist, typically
       [(0), (1), (2), (3), (4), (5), (0, 1) etc.]
    2. For each strain state, find and sort strains and
       stresses by strain value.
    3. Find first and second derivatives of each stress
       with respect to scalar variable corresponding to
       the smallest perturbation in the strain.
    4. Use the pseudoinverse of a matrix-vector expression 
       corresponding to the parameterized stress-strain
       relationship and multiply that matrix by the respective 
       calculated first or second derivatives from the
       previous step.
    5. Place the calculated second and third-order elastic 
       constants appropriately.

    Args:
        strains (nx3x3 array-like): Array of 3x3 strains
            to use in fitting of TOEC and SOEC
        stresses (nx3x3 array-like): Array of 3x3 stresses
            to use in fitting of TOEC and SOEC.  These
            should be PK2 stresses.
        eq_stress (3x3 array-like): stress corresponding to
            equilibrium strain (i. e. "0" strain state).
            If not specified, function will try to find
            the state in the list of provided stresses
            and strains.  If not found, defaults to 0.
        zero_crit (float): value for which strains below
            are ignored in identifying strain states.
    """

    if len(stresses) != len(strains):
        raise ValueError("Length of strains and stresses are not equivalent")
    vstresses = np.array([Stress(stress).voigt for stress in stresses])
    vstrains = np.array([Strain(strain).voigt for strain in strains])
    vstrains[np.abs(vstrains) < zero_crit] = 0

    # Try to find eq_stress if not specified
    if eq_stress is not None:
        veq_stress = Stress(eq_stress).voigt
    else:
        veq_stress = vstresses[np.all(vstrains==0, axis=1)]
        if veq_stress:
            if np.shape(veq_stress) > 1 and not \
               (abs(veq_stress - veq_stress[0]) < 1e-8).all():
                raise ValueError("Multiple stresses found for equilibrium strain"
                                 " state, please specify equilibrium stress or  "
                                 " remove extraneous stresses.")
            veq_stress = veq_stress[0]
        else:
            veq_stress = np.zeros(6)

    # Collect independent strain states:
    independent = set([tuple(np.nonzero(vstrain)[0].tolist())
                       for vstrain in vstrains])
    
    strain_states = []
    dsde = np.zeros((6, len(independent)))
    d2sde2 = np.zeros((6, len(independent)))
    for n, ind in enumerate(independent):
        # match strains with templates
        template = np.zeros(6, dtype=bool)
        np.put(template, ind, True)
        template = np.tile(template, [vstresses.shape[0], 1])
        mode = (template == (np.abs(vstrains) > 1e-10)).all(axis=1)
        mstresses = vstresses[mode]
        mstrains = vstrains[mode]
        # add zero strain state
        mstrains = np.vstack([mstrains, np.zeros(6)])
        mstresses = np.vstack([mstresses, np.zeros(6)])
        # sort strains/stresses by strain values
        mstresses = mstresses[mstrains[:, ind[0]].argsort()]
        mstrains = mstrains[mstrains[:, ind[0]].argsort()]
        strain_states.append(mstrains[-1] / \
                             np.min(mstrains[-1][np.nonzero(mstrains[0])]))
        diff = np.diff(mstrains, axis=0)
        if not (abs(diff - diff[0]) < 1e-8).all():
            raise ValueError("Stencil for strain state {} must be odd-sampling"
                             " centered at 0.".format(ind))
        h = np.min(diff[np.nonzero(diff)])
        coef1 = central_diff(1, len(mstresses))
        coef2 = central_diff(2, len(mstresses))
        if eq_stress is not None:
            mstresses[3] = veq_stress
        dsde[:, n] = np.dot(np.transpose(mstresses), coef1) / h
        d2sde2[:, n] = np.dot(np.transpose(mstresses), coef2) / h**2

    m2i, m3i = generate_pseudo(strain_states)
    s2vec = np.ravel(dsde.T)
    c2vec = np.dot(m2i, s2vec)
    c2 = np.zeros((6, 6))
    c2[np.triu_indices(6)] = c2vec
    c2 = c2 + c2.T - np.diag(np.diag(c2))
    c3 = np.zeros((6, 6, 6))
    s3vec = np.ravel(d2sde2.T)
    c3vec = np.dot(m3i, s3vec)
    list_indices = list(itertools.combinations_with_replacement(range(6), r=3))
    indices_ij = itertools.combinations_with_replacement(range(6), r=3)

    indices = list(itertools.combinations_with_replacement(range(6), r=3))
    for n, (i, j, k) in enumerate(indices):
        c3[i,j,k] = c3[i,k,j] = c3[j,i,k] = c3[j,k,i] = \
                c3[k,i,j] = c3[k,j,i] = c3vec[n]
    return TensorBase.from_voigt(c2), TensorBase.from_voigt(c3)


def central_diff(k, n):
    """
    Generates central difference operator
    """
    A = np.array([(np.linspace(-1, 1, n) * (n-1) / 2)**i \
                  / np.math.factorial(i) for i in range(n)])
    b = np.zeros(n)
    b[k] = 1
    return np.linalg.solve(A, b)

def generate_pseudo(strain_states):
    """
    Generates the pseudoinverse for a given set of strains
    """
    s = sp.Symbol('s')
    nstates = len(strain_states)
    ni = np.array(strain_states)*s
    c2vec, c2arr = get_symbol_list(6, 2)
    c3vec, c3arr = get_symbol_list(6, 3)
    s2arr = np.zeros((nstates, 6), dtype=object)
    s3arr = np.zeros((nstates, 6), dtype=object)
    v_diff = np.vectorize(sp.diff)
    for n, strain_v in enumerate(ni):
        s2arr[n] = [sp.diff(exp, s) for exp in np.dot(c2arr, strain_v)]
        s3arr[n] = [sp.diff(exp, s, 2) 
                    for exp in np.dot(np.dot(c3arr, strain_v), strain_v) / 2]
    s2vec, s3vec = s2arr.ravel(), s3arr.ravel()
    m2 = np.zeros((6*nstates, len(c2vec)))
    m3 = np.zeros((6*nstates, len(c3vec)))
    for n, c in enumerate(c2vec):
        m2[:, n] = v_diff(s2vec, c)
    for n, c in enumerate(c3vec):
        m3[:, n] = v_diff(s3vec, c)
    m2i = np.linalg.pinv(m2)
    m3i = np.linalg.pinv(m3)
    return m2i, m3i

def get_symbol_list(dim, rank):
    indices = list(
        itertools.combinations_with_replacement(range(dim), r=rank))
    c_vec = np.zeros(len(indices), dtype=object)
    c_arr = np.zeros([dim]*rank, dtype=object)
    for n, idx in enumerate(indices):
        c_vec[n] = sp.Symbol('c_'+''.join([str(i) for i in idx]))
        for perm in itertools.permutations(idx):
            c_arr[perm] = c_vec[n]
    return c_vec, c_arr
