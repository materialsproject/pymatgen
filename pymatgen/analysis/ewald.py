"""This module provides classes for calculating the Ewald sum of a structure."""

from __future__ import annotations

import bisect
from copy import copy, deepcopy
from datetime import datetime
from math import log, pi, sqrt
from typing import Any
from warnings import warn

import numpy as np
from monty.json import MSONable
from scipy import constants
from scipy.special import comb, erfc

from pymatgen.core.structure import Structure
from pymatgen.util.due import Doi, due

__author__ = "Shyue Ping Ong, William Davidson Richard"
__copyright__ = "Copyright 2011, The Materials Project"
__credits__ = "Christopher Fischer"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Aug 1 2012"


@due.dcite(
    Doi("10.1016/0010-4655(96)00016-1"),
    description="Ewald summation techniques in perspective: a survey",
    path="pymatgen.analysis.ewald.EwaldSummation",
)
class EwaldSummation(MSONable):
    """
    Calculates the electrostatic energy of a periodic array of charges using
    the Ewald technique.

    Ref:
        Ewald summation techniques in perspective: a survey
        Abdulnour Y. Toukmaji and John A. Board Jr.
        DOI: 10.1016/0010-4655(96)00016-1
        URL: http://www.ee.duke.edu/~ayt/ewaldpaper/ewaldpaper.html

    This matrix can be used to do fast calculations of Ewald sums after species
    removal.

    E = E_recip + E_real + E_point

    Atomic units used in the code, then converted to eV.
    """

    # Converts unit of q*q/r into eV
    CONV_FACT = 1e10 * constants.e / (4 * pi * constants.epsilon_0)

    def __init__(
        self,
        structure,
        real_space_cut=None,
        recip_space_cut=None,
        eta=None,
        acc_factor=12.0,
        w=1 / 2**0.5,
        compute_forces=False,
    ):
        """
        Initializes and calculates the Ewald sum. Default convergence
        parameters have been specified, but you can override them if you wish.

        Args:
            structure (Structure): Input structure that must have proper
                Species on all sites, i.e. Element with oxidation state. Use
                Structure.add_oxidation_state... for example.
            real_space_cut (float): Real space cutoff radius dictating how
                many terms are used in the real space sum. Defaults to None,
                which means determine automatically using the formula given
                in gulp 3.1 documentation.
            recip_space_cut (float): Reciprocal space cutoff radius.
                Defaults to None, which means determine automatically using
                the formula given in gulp 3.1 documentation.
            eta (float): The screening parameter. Defaults to None, which means
                determine automatically.
            acc_factor (float): No. of significant figures each sum is
                converged to.
            w (float): Weight parameter, w, has been included that represents
                the relative computational expense of calculating a term in
                real and reciprocal space. Default of 0.7 reproduces result
                similar to GULP 4.2. This has little effect on the total
                energy, but may influence speed of computation in large
                systems. Note that this parameter is used only when the
                cutoffs are set to None.
            compute_forces (bool): Whether to compute forces. False by
                default since it is usually not needed.
        """
        self._struct = structure
        self._charged = abs(structure.charge) > 1e-8
        self._vol = structure.volume
        self._compute_forces = compute_forces

        self._acc_factor = acc_factor
        # set screening length
        self._eta = eta or (len(structure) * w / (self._vol**2)) ** (1 / 3) * pi
        self._sqrt_eta = sqrt(self._eta)

        # acc factor used to automatically determine the optimal real and
        # reciprocal space cutoff radii
        self._accf = sqrt(log(10**acc_factor))

        self._rmax = real_space_cut or self._accf / self._sqrt_eta
        self._gmax = recip_space_cut or 2 * self._sqrt_eta * self._accf

        # The next few lines pre-compute certain quantities and store them.
        # Ewald summation is rather expensive, and these shortcuts are
        # necessary to obtain several factors of improvement in speedup.
        self._oxi_states = [compute_average_oxidation_state(site) for site in structure]

        self._coords = np.array(self._struct.cart_coords)

        # Define the private attributes to lazy compute reciprocal and real
        # space terms.
        self._initialized = False
        self._recip = self._real = self._point = self._forces = None

        # Compute the correction for a charged cell
        self._charged_cell_energy = (
            -EwaldSummation.CONV_FACT / 2 * np.pi / structure.volume / self._eta * structure.charge**2
        )

    def compute_partial_energy(self, removed_indices):
        """
        Gives total Ewald energy for certain sites being removed, i.e. zeroed
        out.
        """
        total_energy_matrix = self.total_energy_matrix.copy()
        for idx in removed_indices:
            total_energy_matrix[idx, :] = 0
            total_energy_matrix[:, idx] = 0
        return sum(sum(total_energy_matrix))

    def compute_sub_structure(self, sub_structure, tol: float = 1e-3):
        """
        Gives total Ewald energy for an sub structure in the same
        lattice. The sub_structure must be a subset of the original
        structure, with possible different charges.

        Args:
            substructure (Structure): Substructure to compute Ewald sum for.
            tol (float): Tolerance for site matching in fractional coordinates.

        Returns:
            Ewald sum of substructure.
        """
        total_energy_matrix = self.total_energy_matrix.copy()

        def find_match(site):
            for test_site in sub_structure:
                frac_diff = abs(np.array(site.frac_coords) - np.array(test_site.frac_coords)) % 1
                frac_diff = [abs(a) < tol or abs(a) > 1 - tol for a in frac_diff]
                if all(frac_diff):
                    return test_site
            return None

        matches = []
        for i, site in enumerate(self._struct):
            matching_site = find_match(site)
            if matching_site:
                new_charge = compute_average_oxidation_state(matching_site)
                old_charge = self._oxi_states[i]
                scaling_factor = new_charge / old_charge
                matches.append(matching_site)
            else:
                scaling_factor = 0
            total_energy_matrix[i, :] *= scaling_factor
            total_energy_matrix[:, i] *= scaling_factor

        if len(matches) != len(sub_structure):
            output = ["Missing sites."]
            for site in sub_structure:
                if site not in matches:
                    output.append(f"unmatched = {site}")
            raise ValueError("\n".join(output))

        return sum(sum(total_energy_matrix))

    @property
    def reciprocal_space_energy(self):
        """The reciprocal space energy."""
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True
        return sum(sum(self._recip))

    @property
    def reciprocal_space_energy_matrix(self):
        """
        The reciprocal space energy matrix. Each matrix element (i, j)
        corresponds to the interaction energy between site i and site j in
        reciprocal space.
        """
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True
        return self._recip

    @property
    def real_space_energy(self):
        """The real space energy."""
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True
        return sum(sum(self._real))

    @property
    def real_space_energy_matrix(self):
        """
        The real space energy matrix. Each matrix element (i, j) corresponds to
        the interaction energy between site i and site j in real space.
        """
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True
        return self._real

    @property
    def point_energy(self):
        """The point energy."""
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True
        return sum(self._point)

    @property
    def point_energy_matrix(self):
        """
        The point space matrix. A diagonal matrix with the point terms for each
        site in the diagonal elements.
        """
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True
        return self._point

    @property
    def total_energy(self):
        """The total energy."""
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True
        return sum(sum(self._recip)) + sum(sum(self._real)) + sum(self._point) + self._charged_cell_energy

    @property
    def total_energy_matrix(self):
        """
        The total energy matrix. Each matrix element (i, j) corresponds to the
        total interaction energy between site i and site j.

        Note that this does not include the charged-cell energy, which is only important
        when the simulation cell is not charge balanced.
        """
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True

        total_energy = self._recip + self._real
        for idx, energy in enumerate(self._point):
            total_energy[idx, idx] += energy
        return total_energy

    @property
    def forces(self):
        """
        The forces on each site as a Nx3 matrix. Each row corresponds to a
        site.
        """
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True

        if not self._compute_forces:
            raise AttributeError("Forces are available only if compute_forces is True!")
        return self._forces

    def get_site_energy(self, site_index):
        """Compute the energy for a single site in the structure.

        Args:
            site_index (int): Index of site
        ReturnS:
        (float) - Energy of that site
        """
        if not self._initialized:
            self._calc_ewald_terms()
            self._initialized = True

        if self._charged:
            warn("Per atom energies for charged structures not supported in EwaldSummation")
        return np.sum(self._recip[:, site_index]) + np.sum(self._real[:, site_index]) + self._point[site_index]

    def _calc_ewald_terms(self):
        """Calculates and sets all Ewald terms (point, real and reciprocal)."""
        self._recip, recip_forces = self._calc_recip()
        self._real, self._point, real_point_forces = self._calc_real_and_point()
        if self._compute_forces:
            self._forces = recip_forces + real_point_forces

    def _calc_recip(self):
        """
        Perform the reciprocal space summation. Calculates the quantity
        E_recip = 1/(2PiV) sum_{G < Gmax} exp(-(G.G/4/eta))/(G.G) S(G)S(-G)
        where
        S(G) = sum_{k=1,N} q_k exp(-i G.r_k)
        S(G)S(-G) = |S(G)|**2.

        This method is heavily vectorized to utilize numpy's C backend for
        speed.
        """
        n_sites = len(self._struct)
        prefactor = 2 * pi / self._vol
        e_recip = np.zeros((n_sites, n_sites), dtype=np.float_)
        forces = np.zeros((n_sites, 3), dtype=np.float_)
        coords = self._coords
        rcp_latt = self._struct.lattice.reciprocal_lattice
        recip_nn = rcp_latt.get_points_in_sphere([[0, 0, 0]], [0, 0, 0], self._gmax)

        frac_coords = [fcoords for (fcoords, dist, i, img) in recip_nn if dist != 0]

        gs = rcp_latt.get_cartesian_coords(frac_coords)
        g2s = np.sum(gs**2, 1)
        exp_vals = np.exp(-g2s / (4 * self._eta))
        grs = np.sum(gs[:, None] * coords[None, :], 2)

        oxi_states = np.array(self._oxi_states)

        # create array where q_2[i,j] is qi * qj
        qiqj = oxi_states[None, :] * oxi_states[:, None]

        # calculate the structure factor
        s_reals = np.sum(oxi_states[None, :] * np.cos(grs), 1)
        s_imags = np.sum(oxi_states[None, :] * np.sin(grs), 1)

        for g, g2, gr, expval, sreal, simag in zip(gs, g2s, grs, exp_vals, s_reals, s_imags):
            # Uses the identity sin(x)+cos(x) = 2**0.5 sin(x + pi/4)
            m = (gr[None, :] + pi / 4) - gr[:, None]
            np.sin(m, m)
            m *= expval / g2

            e_recip += m

            if self._compute_forces:
                pref = 2 * expval / g2 * oxi_states
                factor = prefactor * pref * (sreal * np.sin(gr) - simag * np.cos(gr))

                forces += factor[:, None] * g[None, :]

        forces *= EwaldSummation.CONV_FACT
        e_recip *= prefactor * EwaldSummation.CONV_FACT * qiqj * 2**0.5
        return e_recip, forces

    def _calc_real_and_point(self):
        """Determines the self energy -(eta/pi)**(1/2) * sum_{i=1}^{N} q_i**2."""
        frac_coords = self._struct.frac_coords
        force_pf = 2 * self._sqrt_eta / sqrt(pi)
        coords = self._coords
        n_sites = len(self._struct)
        e_real = np.empty((n_sites, n_sites), dtype=np.float_)

        forces = np.zeros((n_sites, 3), dtype=np.float_)

        qs = np.array(self._oxi_states)

        epoint = -(qs**2) * sqrt(self._eta / pi)

        for i in range(n_sites):
            nf_coords, rij, js, _ = self._struct.lattice.get_points_in_sphere(
                frac_coords, coords[i], self._rmax, zip_results=False
            )

            # remove the rii term
            inds = rij > 1e-8
            js = js[inds]
            rij = rij[inds]
            nf_coords = nf_coords[inds]

            qi = qs[i]
            qj = qs[js]

            erfc_val = erfc(self._sqrt_eta * rij)
            new_ereals = erfc_val * qi * qj / rij

            # insert new_ereals
            for key in range(n_sites):
                e_real[key, i] = np.sum(new_ereals[js == key])

            if self._compute_forces:
                nc_coords = self._struct.lattice.get_cartesian_coords(nf_coords)

                fijpf = qj / rij**3 * (erfc_val + force_pf * rij * np.exp(-self._eta * rij**2))
                forces[i] += np.sum(
                    np.expand_dims(fijpf, 1) * (np.array([coords[i]]) - nc_coords) * qi * EwaldSummation.CONV_FACT,
                    axis=0,
                )

        e_real *= 0.5 * EwaldSummation.CONV_FACT
        epoint *= EwaldSummation.CONV_FACT
        return e_real, epoint, forces

    @property
    def eta(self):
        """Returns: eta value used in Ewald summation."""
        return self._eta

    def __str__(self):
        output = [
            f"Real = {self.real_space_energy}",
            f"Reciprocal = {self.reciprocal_space_energy}",
            f"Point = {self.point_energy}",
            f"Total = {self.total_energy}",
            f"Forces:\n{self.forces}" if self._compute_forces else "Forces were not computed",
        ]
        return "\n".join(output)

    def as_dict(self, verbosity: int = 0) -> dict:
        """
        Json-serialization dict representation of EwaldSummation.

        Args:
            verbosity (int): Verbosity level. Default of 0 only includes the
                matrix representation. Set to 1 for more details.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "structure": self._struct.as_dict(),
            "compute_forces": self._compute_forces,
            "eta": self._eta,
            "acc_factor": self._acc_factor,
            "real_space_cut": self._rmax,
            "recip_space_cut": self._gmax,
            "_recip": None if self._recip is None else self._recip.tolist(),
            "_real": None if self._real is None else self._real.tolist(),
            "_point": None if self._point is None else self._point.tolist(),
            "_forces": None if self._forces is None else self._forces.tolist(),
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any], fmt: str | None = None, **kwargs) -> EwaldSummation:
        """Create an EwaldSummation instance from JSON-serialized dictionary.

        Args:
            d (dict): Dictionary representation
            fmt (str, optional): Unused. Defaults to None.

        Returns:
            EwaldSummation: class instance
        """
        summation = cls(
            structure=Structure.from_dict(d["structure"]),
            real_space_cut=d["real_space_cut"],
            recip_space_cut=d["recip_space_cut"],
            eta=d["eta"],
            acc_factor=d["acc_factor"],
            compute_forces=d["compute_forces"],
        )

        # set previously computed private attributes
        if d["_recip"] is not None:
            summation._recip = np.array(d["_recip"])
            summation._real = np.array(d["_real"])
            summation._point = np.array(d["_point"])
            summation._forces = np.array(d["_forces"])
            summation._initialized = True

        return summation


class EwaldMinimizer:
    """
    This class determines the manipulations that will minimize an Ewald matrix,
    given a list of possible manipulations. This class does not perform the
    manipulations on a structure, but will return the list of manipulations
    that should be done on one to produce the minimal structure. It returns the
    manipulations for the n lowest energy orderings. This class should be used
    to perform fractional species substitution or fractional species removal to
    produce a new structure. These manipulations create large numbers of
    candidate structures, and this class can be used to pick out those with the
    lowest Ewald sum.

    An alternative (possibly more intuitive) interface to this class is the
    order disordered structure transformation.

    Author - Will Richards
    """

    ALGO_FAST = 0
    ALGO_COMPLETE = 1
    ALGO_BEST_FIRST = 2

    """
    ALGO_TIME_LIMIT: Slowly increases the speed (with the cost of decreasing
    accuracy) as the minimizer runs. Attempts to limit the run time to
    approximately 30 minutes.
    """
    ALGO_TIME_LIMIT = 3

    def __init__(self, matrix, m_list, num_to_return=1, algo=ALGO_FAST):
        """
        Args:
            matrix: A matrix of the Ewald sum interaction energies. This is stored
                in the class as a diagonally symmetric array and so
                self._matrix will not be the same as the input matrix.
            m_list: list of manipulations. each item is of the form
                (multiplication fraction, number_of_indices, indices, species)
                These are sorted such that the first manipulation contains the
                most permutations. this is actually evaluated last in the
                recursion since I'm using pop.
            num_to_return: The minimizer will find the number_returned lowest
                energy structures. This is likely to return a number of duplicate
                structures so it may be necessary to overestimate and then
                remove the duplicates later. (duplicate checking in this
                process is extremely expensive).
        """
        # Setup and checking of inputs
        self._matrix = copy(matrix)
        # Make the matrix diagonally symmetric (so matrix[i,:] == matrix[:,j])
        for ii in range(len(self._matrix)):
            for jj in range(ii, len(self._matrix)):
                value = (self._matrix[ii, jj] + self._matrix[jj, ii]) / 2
                self._matrix[ii, jj] = value
                self._matrix[jj, ii] = value

        # sort the m_list based on number of permutations
        self._m_list = sorted(m_list, key=lambda x: comb(len(x[2]), x[1]), reverse=True)

        for mlist in self._m_list:
            if mlist[0] > 1:
                raise ValueError("multiplication fractions must be <= 1")
        self._current_minimum = float("inf")
        self._num_to_return = num_to_return
        self._algo = algo
        if algo == EwaldMinimizer.ALGO_COMPLETE:
            raise NotImplementedError("Complete algo not yet implemented for EwaldMinimizer")

        self._output_lists = []
        # Tag that the recurse function looks at each level. If a method
        # sets this to true it breaks the recursion and stops the search.
        self._finished = False

        self._start_time = datetime.utcnow()

        self.minimize_matrix()

        self._best_m_list = self._output_lists[0][1]
        self._minimized_sum = self._output_lists[0][0]

    def minimize_matrix(self):
        """
        This method finds and returns the permutations that produce the lowest
        Ewald sum calls recursive function to iterate through permutations.
        """
        if self._algo in (EwaldMinimizer.ALGO_FAST, EwaldMinimizer.ALGO_BEST_FIRST):
            return self._recurse(self._matrix, self._m_list, set(range(len(self._matrix))))
        return None

    def add_m_list(self, matrix_sum, m_list):
        """
        This adds an m_list to the output_lists and updates the current
        minimum if the list is full.
        """
        if self._output_lists is None:
            self._output_lists = [[matrix_sum, m_list]]
        else:
            bisect.insort(self._output_lists, [matrix_sum, m_list])
        if self._algo == EwaldMinimizer.ALGO_BEST_FIRST and len(self._output_lists) == self._num_to_return:
            self._finished = True
        if len(self._output_lists) > self._num_to_return:
            self._output_lists.pop()
        if len(self._output_lists) == self._num_to_return:
            self._current_minimum = self._output_lists[-1][0]

    def best_case(self, matrix, m_list, indices_left):
        """
        Computes a best case given a matrix and manipulation list.

        Args:
            matrix: the current matrix (with some permutations already
                performed)
            m_list: [(multiplication fraction, number_of_indices, indices,
                species)] describing the manipulation
            indices: Set of indices which haven't had a permutation
                performed on them.
        """
        m_indices = []
        fraction_list = []
        for m in m_list:
            m_indices.extend(m[2])
            fraction_list.extend([m[0]] * m[1])

        indices = list(indices_left.intersection(m_indices))

        interaction_matrix = matrix[indices, :][:, indices]

        fractions = np.zeros(len(interaction_matrix)) + 1
        fractions[: len(fraction_list)] = fraction_list
        fractions = np.sort(fractions)

        # Sum associated with each index (disregarding interactions between
        # indices)
        sums = 2 * np.sum(matrix[indices], axis=1)
        sums = np.sort(sums)

        # Interaction corrections. Can be reduced to (1-x)(1-y) for x,y in
        # fractions each element in a column gets multiplied by (1-x), and then
        # the sum of the columns gets multiplied by (1-y) since fractions are
        # less than 1, there is no effect of one choice on the other
        step1 = np.sort(interaction_matrix) * (1 - fractions)
        step2 = np.sort(np.sum(step1, axis=1))
        step3 = step2 * (1 - fractions)
        interaction_correction = np.sum(step3)

        if self._algo == self.ALGO_TIME_LIMIT:
            elapsed_time = datetime.utcnow() - self._start_time
            speedup_parameter = elapsed_time.total_seconds() / 1800
            avg_int = np.sum(interaction_matrix, axis=None)
            avg_frac = np.average(np.outer(1 - fractions, 1 - fractions))
            average_correction = avg_int * avg_frac

            interaction_correction = average_correction * speedup_parameter + interaction_correction * (
                1 - speedup_parameter
            )

        return np.sum(matrix) + np.inner(sums[::-1], fractions - 1) + interaction_correction

    @classmethod
    def get_next_index(cls, matrix, manipulation, indices_left):
        """
        Returns an index that should have the most negative effect on the
        matrix sum.
        """
        # pylint: disable=E1126
        f = manipulation[0]
        indices = list(indices_left.intersection(manipulation[2]))
        sums = np.sum(matrix[indices], axis=1)
        return indices[sums.argmax(axis=0)] if f < 1 else indices[sums.argmin(axis=0)]

    def _recurse(self, matrix, m_list, indices, output_m_list=None):
        """
        This method recursively finds the minimal permutations using a binary
        tree search strategy.

        Args:
            matrix: The current matrix (with some permutations already
                performed).
            m_list: The list of permutations still to be performed
            indices: Set of indices which haven't had a permutation
                performed on them.
        """
        # check to see if we've found all the solutions that we need
        if self._finished:
            return

        if output_m_list is None:
            output_m_list = []

        # if we're done with the current manipulation, pop it off.
        while m_list[-1][1] == 0:
            m_list = copy(m_list)
            m_list.pop()
            # if there are no more manipulations left to do check the value
            if not m_list:
                matrix_sum = np.sum(matrix)
                if matrix_sum < self._current_minimum:
                    self.add_m_list(matrix_sum, output_m_list)
                return

        # if we won't have enough indices left, return
        if m_list[-1][1] > len(indices.intersection(m_list[-1][2])):
            return

        if (len(m_list) == 1 or m_list[-1][1] > 1) and self.best_case(matrix, m_list, indices) > self._current_minimum:
            return

        index = self.get_next_index(matrix, m_list[-1], indices)

        m_list[-1][2].remove(index)

        # Make the matrix and new m_list where we do the manipulation to the
        # index that we just got
        matrix2 = np.copy(matrix)
        m_list2 = deepcopy(m_list)
        output_m_list2 = copy(output_m_list)

        matrix2[index, :] *= m_list[-1][0]
        matrix2[:, index] *= m_list[-1][0]
        output_m_list2.append([index, m_list[-1][3]])
        indices2 = copy(indices)
        indices2.remove(index)
        m_list2[-1][1] -= 1

        # recurse through both the modified and unmodified matrices

        self._recurse(matrix2, m_list2, indices2, output_m_list2)
        self._recurse(matrix, m_list, indices, output_m_list)

    @property
    def best_m_list(self):
        """Returns: Best m_list found."""
        return self._best_m_list

    @property
    def minimized_sum(self):
        """Returns: Minimized sum."""
        return self._minimized_sum

    @property
    def output_lists(self):
        """Returns: output lists."""
        return self._output_lists


def compute_average_oxidation_state(site):
    """
    Calculates the average oxidation state of a site.

    Args:
        site: Site to compute average oxidation state

    Returns:
        Average oxidation state of site.
    """
    try:
        return sum(sp.oxi_state * occu for sp, occu in site.species.items() if sp is not None)
    except AttributeError:
        pass
    try:
        return site.charge
    except AttributeError:
        raise ValueError(
            "Ewald summation can only be performed on structures "
            "that are either oxidation state decorated or have "
            "site charges."
        )
