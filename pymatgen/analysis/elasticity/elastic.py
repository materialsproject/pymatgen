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
        voigt_map, reverse_voigt_map
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity.strain import Strain
import numpy as np
import sympy as sp
import warnings
import itertools
from six.moves import range

__author__ = "Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Joseph Montoya, Shyam Dwaraknath, Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Development"
__date__ = "March 22, 2012"


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
        inds = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
        c_ij = np.array([[np.polyfit([strain[ind1] for strain in list(stress_dict.keys())
                                      if (strain.i, strain.j) == ind1],
                                     [stress_dict[strain][ind2] for strain
                                      in list(stress_dict.keys())
                                      if (strain.i, strain.j) == ind1], 1)[0]
                          for ind1 in inds] for ind2 in inds])
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


nf_eval = 1
import scipy.optimize as opt
class ToecFitter(object):
    """
    Optimizer class for third-order elastic constants
    """
    
    def __init__(self, strains, stresses, structure=None):
        self.strains = strains
        self.stresses = stresses
        self.structure = None
        self.c3_indices = list(itertools.combinations_with_replacement(range(6), r=3))

    # How much does enforcing symmetry save?  What are the appropriate symmetries?
    def opt_func(self, super_vec, strains, stresses):
        c3_vec = super_vec[:56]
        c2_vec = super_vec[56:]
        assert len(c2_vec) == 21
        C3 = np.zeros((6, 6, 6))
        C2 = np.zeros((6, 6))
        
        # Construct c_ijkl
        C2[np.triu_indices(6)] = c2_vec
        C2 = C2 + C2.T - np.diag(np.diag(C2))
        c_ijkl = TensorBase.from_voigt(C2)

        # Construct c_ijlkmn
        for n, (i, j, k) in enumerate(self.c3_indices):
            C3[i,j,k] = C3[i,k,j] = C3[j,i,k] = C3[j,k,i] = \
                    C3[k,i,j] = C3[k,j,i] = c3_vec[n]
        c_ijklmn = TensorBase.from_voigt(C3)
        #total_resid = np.zeros((3, 3))
        total_norm = 0
        for stress, strain in zip(stresses, strains):
            resid =  np.einsum("ijkl,kl->ij", c_ijkl, strain) \
                    + 0.5*np.einsum("ijklmn,kl,mn->ij",c_ijklmn, strain, strain) \
                    - stress
            total_norm += np.linalg.norm(resid)
        return total_norm

    def get_coeff(self, strains, stresses, symm_init=False):
        guess = self.gen_init()#symm=symm_init)
        result = opt.minimize(self.opt_func, guess, args = (strains, stresses),
                              callback=self.callback_f)
        if result.success:
            return result.x
        else:
            raise ValueError("Optimizer failed with message: {}".format(result.message))

    def gen_init(self, symm=True):
        t1 = 50*np.ones(56)
        t2 = 50*np.ones(21)
        """
        if symm:
            t1 = t1.fit_to_structure(self.structure)
            t2 = t2.fit_to_structure(self.structure)
        """
        return np.concatenate((t1.ravel(), t2.ravel()))

    def callback_f(self, resid):
        global nf_eval
        print("{}: {}".format(nf_eval, resid))
        nf_eval += 1


def new_fit(strains, stresses, structure = None, 
            output=None, eq_stress = None, zero_crit=1e-10):
    """
    Temporary home for I. Winter's TOEC fitting function.

    Args:
        strains nx3x3 array-like: Array of 3x3 strains
            to use in fitting of TOEC and SOEC
        stresses nx3x3 array-like: Array of 3x3 stresses
            to use in fitting of TOEC and SOEC.  These
            should be PK2 stresses.
        structure: Structure
        output: name for file output, if none doesn't output file
    """

    # perturbed strain states
    #M2I = np.genfromtxt("SM2I.csv", delimiter=",")
    #M3I = np.genfromtxt("SM3I.csv", delimiter=",")
    if len(stresses) != len(strains):
        raise ValueError("Length of strains and stresses are not equivalent")
    vstresses = np.array([Stress(stress).voigt for stress in stresses])
    vstrains = np.array([Strain(strain).voigt for strain in strains])
    vstrains[np.abs(vstrains) < zero_crit] = 0

    # Try to find eq_stress if not specified
    if eq_stress:
        veq_stress = Stress(eq_stress).voigt
    else:
        veq_stress = vstresses[np.all(vstrains==0, axis=0)]
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
    #jj = [0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 3, 3, 4, 1]
    #kk = [0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 4, 5, 5, 2]
    #gamma = np.linspace(-0.05, 0.05, 7)
    #gdiag = gamma.copy()

    # h = np.diff(gamma)[0]
    sig = np.zeros([6, 14, 7])
    #coef1 = central_diff(1, 7)
    #coef2 = central_diff(2, 7)
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
        if not (abs(diff - diff[0]) < 1e-10).all():
            raise ValueError("Stencil for strain state {} must be odd-sampling"
                             " centered at 0.".format(ind))
        h = np.min(diff[np.nonzero(diff)])
        coef1 = central_diff(1, len(mstresses))
        coef2 = central_diff(2, len(mstresses))

        # This is a bit nebulous, I've checked it using an assertion statement
        # but should maybe refactor to be a bit clearer
        #import pdb; pdb.set_trace()
        # HACK
        if eq_stress is not None:
            # eq_stress in standard notation
            mstresses[3] = veq_stress
        sig[:, n, :] = np.transpose(mstresses)
        #import pdb; pdb.set_trace()
        #if n==5:
            #import pdb; pdb.set_trace()
        #assert (sig2[:, n, :] == sig[:, n, :]).all()
        #import pdb; pdb.set_trace()
        dsde[:, n] = np.dot(np.transpose(mstresses), coef1) / h
        d2sde2[:, n] = np.dot(np.transpose(mstresses), coef2) / h**2
        #import pdb; pdb.set_trace()

    M2I, M3I = generate_pseudo(strain_states)
    s2vec = np.ravel(dsde.T)
    c2vec = np.dot(M2I, s2vec)
    C2 = np.zeros((6, 6))
    C2[np.triu_indices(6)] = c2vec
    C2 = C2 + C2.T - np.diag(np.diag(C2))
    C3 = np.zeros((6, 6, 6))
    s3vec = np.ravel(d2sde2.T)
    c3vec = np.dot(M3I, s3vec)
    #import pdb; pdb.set_trace()
    list_indices = list(itertools.combinations_with_replacement(range(6), r=3))
    indices_ij = itertools.combinations_with_replacement(range(6), r=3)

    indices = list(itertools.combinations_with_replacement(range(6), r=3))
    txt = ''
    for n, (i, j, k) in enumerate(indices):
        C3[i,j,k] = C3[i,k,j] = C3[j,i,k] = C3[j,k,i] = \
                C3[k,i,j] = C3[k,j,i] = c3vec[n]
        txt += '\nc_'+''.join(str(m+1) for m in [i, j, k])\
                + ' = {}'.format(c3vec[n])
    if output:
        with open("C_ijk_{}".format(output), 'w') as f:
            f.write(txt)

    c3_tens = TensorBase.from_voigt(C3)
    if structure:
        C3_sym = c3_tens.fit_to_structure(structure).voigt
        sym_txt = ''
        for (i, j, k) in indices:
            sym_txt += '\nc_'+''.join(str(m+1) for m in [i, j, k])\
                + ' = {}'.format(c3vec[n])
        with open("C_ijk_sym_{}".format(sys.argv[1].split('_')[0]), "w") as f:
            f.write(sym_txt)

    #import pdb; pdb.set_trace()
    return C2, C3

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
        s2arr[n] = v_diff(np.dot(c2arr, strain_v), s)
        s3arr[n] = v_diff(np.dot(np.dot(c3arr, strain_v), strain_v) / 2, s, 2)
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


if __name__=="__main__":
    import json, sys, pdb, traceback
    with open(sys.argv[1]) as f:
        sdict = json.load(f)
    strains = sdict["strains"]
    stresses = sdict["stresses"]
    pk_stresses = sdict["pk_stresses"]
    try:
        C2, C3 = new_fit(strains, pk_stresses)
        cijkl = TensorBase.from_voigt(C2)
        cijklmn = TensorBase.from_voigt(C3)
        model_strain = Strain(strains[0])
        model_stress = 0.5 * np.einsum("ijklmn,kl,mn->ij", cijklmn, model_strain, model_strain) \
                + np.einsum("ijkl,kl->ij", cijkl, model_strain)
        resid = np.array(pk_stresses[0]) - model_stress
        #generate_pseudo()
    except:
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
    """
    try:
        toec_fitter = ToecFitter(strains, pk_stresses)
        vec = toec_fitter.get_coeff(strains, pk_stresses)
    except:
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
        """


