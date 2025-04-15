"""
This module implements the Quasi-harmonic Debye approximation that can
be used to compute thermal properties.

See the following papers for more info:

    https://doi.org/10.1016/j.comphy.2003.12.001 (2004)
    https://doi.org/10.1103/PhysRevB.90.174107 (2014)
"""

from __future__ import annotations

import logging
from collections import defaultdict

import numpy as np
from monty.dev import deprecated
from scipy import linalg
from scipy.constants import physical_constants
from scipy.integrate import quad
from scipy.optimize import minimize

from pymatgen.analysis.eos import EOS, PolynomialEOS
from pymatgen.core.units import FloatWithUnit
from pymatgen.util.due import Doi, due

__author__ = "Kiran Mathew, Brandon Bocklund"
__credits__ = "Cormac Toher"


logger = logging.getLogger(__name__)


cite_gibbs = due.dcite(
    Doi("10.1016/j.comphy.2003.12.001"),
    description="GIBBS: isothermal-isobaric thermodynamics of solids from energy curves using a "
    "quasi-harmonic Debye model",
    path="pymatgen.analysis.quasiharmonic",
)


@cite_gibbs
@due.dcite(
    Doi("10.1103/PhysRevB.90.174107"),
    description="High-throughput computational screening of thermal conductivity, Debye "
    "temperature, and Grüneisen parameter using a quasiharmonic Debye model",
    path="pymatgen.analysis.quasiharmonic",
)
class QuasiHarmonicDebyeApprox:
    """Quasi-harmonic approximation."""

    def __init__(
        self,
        energies,
        volumes,
        structure,
        t_min=300.0,
        t_step=100,
        t_max=300.0,
        eos="vinet",
        pressure=0.0,
        poisson=0.25,
        use_mie_gruneisen=False,
        anharmonic_contribution=False,
    ) -> None:
        """
        Args:
            energies (list): list of DFT energies in eV
            volumes (list): list of volumes in Ang^3
            structure (Structure): pymatgen structure object
            t_min (float): min temperature
            t_step (float): temperature step
            t_max (float): max temperature
            eos (str): equation of state used for fitting the energies and the volumes.
                options supported by pymatgen: "quadratic", "murnaghan", "birch",
                    "birch_murnaghan", "pourier_tarantola", "vinet", "deltafactor", "numerical_eos"
            pressure (float): in GPa, optional.
            poisson (float): poisson ratio.
            use_mie_gruneisen (bool): whether or not to use the mie-gruneisen
                formulation to compute the gruneisen parameter.
                The default is the slater-gamma formulation.
            anharmonic_contribution (bool): whether or not to consider the anharmonic
                contribution to the Debye temperature. Cannot be used with
                use_mie_gruneisen. Defaults to False.
        """
        self.energies = energies
        self.volumes = volumes
        self.structure = structure
        self.temperature_min = t_min
        self.temperature_max = t_max
        self.temperature_step = t_step
        self.eos_name = eos
        self.pressure = pressure
        self.poisson = poisson
        self.use_mie_gruneisen = use_mie_gruneisen
        self.anharmonic_contribution = anharmonic_contribution
        if self.use_mie_gruneisen and self.anharmonic_contribution:
            raise ValueError(
                "The Mie-Gruneisen formulation and anharmonic contribution are circular referenced and "
                "cannot be used together."
            )
        self.mass = sum(spec.atomic_mass for spec in self.structure.species)
        self.natoms = self.structure.composition.num_atoms
        self.avg_mass = physical_constants["atomic mass constant"][0] * self.mass / self.natoms  # kg
        self.kb = physical_constants["Boltzmann constant in eV/K"][0]
        self.hbar = physical_constants["Planck constant over 2 pi in eV s"][0]
        self.gpa_to_ev_ang = 1 / 160.21766208  # 1 GPa in eV/Ang^3
        self.gibbs_free_energy: list[float] = []  # optimized values, eV
        # list of temperatures for which the optimized values are available, K
        self.temperatures: list[float] = []
        self.optimum_volumes: list[float] = []  # in Ang^3
        # fit E and V and get the bulk modulus(used to compute the Debye
        # temperature)
        logger.info("Fitting E and V")
        self.eos = EOS(eos)
        self.ev_eos_fit = self.eos.fit(volumes, energies)
        self.bulk_modulus = self.ev_eos_fit.b0_GPa  # in GPa
        self.optimize_gibbs_free_energy()

    def optimize_gibbs_free_energy(self):
        """Evaluate the Gibbs free energy as a function of V, T and P i.e
        G(V, T, P), minimize G(V, T, P) w.r.t. V for each T and store the
        optimum values.

        Note: The data points for which the equation of state fitting fails
            are skipped.
        """
        temperatures = np.linspace(
            self.temperature_min,
            self.temperature_max,
            int(np.ceil((self.temperature_max - self.temperature_min) / self.temperature_step) + 1),
        )

        for temp in temperatures:
            G_opt = V_opt = None
            try:
                G_opt, V_opt = self.optimizer(temp)
            except Exception:
                if len(temperatures) <= 1:
                    raise
                logger.info(f"EOS fitting failed, skipping {temp=}")
            self.gibbs_free_energy.append(G_opt)
            self.temperatures.append(temp)
            self.optimum_volumes.append(V_opt)

    def optimizer(self, temperature):
        """Evaluate G(V, T, P) at the given temperature(and pressure) and minimize it w.r.t. V.

        1. Compute the vibrational Helmholtz free energy, A_vib.
        2. Compute the Gibbs free energy as a function of volume, temperature
            and pressure, G(V,T,P).
        3. Perform an equation of state fit to get the functional form of
            Gibbs free energy:G(V, T, P).
        4. Finally G(V, P, T) is minimized with respect to V.

        Args:
            temperature (float): temperature in K

        Returns:
            float, float: G_opt(V_opt, T, P) in eV and V_opt in Ang^3.
        """
        G_V = []  # G for each volume
        # G = E(V) + PV + A_vib(V, T)
        for idx, v in enumerate(self.volumes):
            G_V.append(
                self.energies[idx]
                + self.pressure * v * self.gpa_to_ev_ang
                + self.vibrational_free_energy(temperature, v)
            )

        # fit equation of state, G(V, T, P)
        eos_fit = self.eos.fit(self.volumes, G_V)
        # minimize the fit EoS w.r.t. volume
        # Note: the ref energy and the ref volume(E0 and V0) not necessarily
        # the same as minimum energy and min volume.
        volume_guess = eos_fit.volumes[np.argmin(eos_fit.energies)]
        min_wrt_vol = minimize(eos_fit.func, volume_guess)
        # G_opt=G(V_opt, T, P), V_opt
        return min_wrt_vol.fun, min_wrt_vol.x[0]

    @cite_gibbs
    def vibrational_free_energy(self, temperature, volume):
        """
        Vibrational Helmholtz free energy, A_vib(V, T).
        Eq(4) in doi.org/10.1016/j.comphy.2003.12.001.

        Args:
            temperature (float): temperature in K
            volume (float)

        Returns:
            float: vibrational free energy in eV
        """
        y = self.debye_temperature(volume) / temperature
        return (
            self.kb * self.natoms * temperature * (9.0 / 8.0 * y + 3 * np.log(1 - np.exp(-y)) - self.debye_integral(y))
        )

    @cite_gibbs
    def vibrational_internal_energy(self, temperature, volume):
        """
        Vibrational internal energy, U_vib(V, T).
        Eq(4) in doi.org/10.1016/j.comphy.2003.12.001.

        Args:
            temperature (float): temperature in K
            volume (float): in Ang^3

        Returns:
            float: vibrational internal energy in eV
        """
        y = self.debye_temperature(volume) / temperature
        return self.kb * self.natoms * temperature * (9.0 / 8.0 * y + 3 * self.debye_integral(y))

    @cite_gibbs
    def debye_temperature(self, volume: float) -> float:
        """
        Calculates the Debye temperature.
        Eq(6) in doi.org/10.1016/j.comphy.2003.12.001. Thanks to Joey.

        Eq(6) above is equivalent to Eq(3) in doi.org/10.1103/PhysRevB.37.790
        which does not consider anharmonic effects. Eq(20) in the same paper
        and Eq(18) in doi.org/10.1016/j.commatsci.2009.12.006 both consider
        anharmonic contributions to the Debye temperature through the Gruneisen
        parameter at 0K (Gruneisen constant).

        The anharmonic contribution is toggled by setting the anharmonic_contribution
        to True or False in the QuasiHarmonicDebyeApprox constructor.

        Args:
            volume (float): in Ang^3

        Returns:
            float: Debye temperature in K
        """
        term1 = (2.0 / 3.0 * (1.0 + self.poisson) / (1.0 - 2.0 * self.poisson)) ** 1.5
        term2 = (1 / 3 * (1.0 + self.poisson) / (1.0 - self.poisson)) ** 1.5
        f_sigma = (3.0 / (2.0 * term1 + term2)) ** (1 / 3)
        debye = (
            2.9772e-11 * (volume / self.natoms) ** (-1.0 / 6.0) * f_sigma * np.sqrt(self.bulk_modulus / self.avg_mass)
        )
        if self.anharmonic_contribution:
            gamma = self.gruneisen_parameter(0, self.ev_eos_fit.v0)  # 0K equilibrium Gruneisen parameter
            return debye * (self.ev_eos_fit.v0 / volume) ** (gamma)
        return debye

    @cite_gibbs
    @staticmethod
    def debye_integral(y):
        """
        Debye integral. Eq(5) in doi.org/10.1016/j.comphy.2003.12.001.

        Args:
            y (float): Debye temperature / T, upper limit

        Returns:
            float: unitless
        """
        # floating point limit is reached around y=155, so values beyond that
        # are set to the limiting value(T-->0, y --> \infty) of
        # 6.4939394 (from wolfram alpha).
        factor = 3.0 / y**3
        if y < 155:
            integral = quad(lambda x: x**3 / (np.exp(x) - 1.0), 0, y)
            return next(iter(integral)) * factor
        return 6.493939 * factor

    @cite_gibbs
    def gruneisen_parameter(self, temperature, volume):
        """Slater-gamma formulation(the default):
            gruneisen parameter = - d log(theta)/ d log(V) = - (1/6 + 0.5 d log(B)/ d log(V))
                                = - (1/6 + 0.5 V/B dB/dV), where dB/dV = d^2E/dV^2 + V * d^3E/dV^3.

        Mie-gruneisen formulation:
            Eq(31) in doi.org/10.1016/j.comphy.2003.12.001
            Eq(7) in Blanco et. al. Journal of Molecular Structure (Theochem)
                368 (1996) 245-255
            Also see J.P. Poirier, Introduction to the Physics of the Earth's
                Interior, 2nd ed. (Cambridge University Press, Cambridge,
                2000) Eq(3.53)

        Args:
            temperature (float): temperature in K
            volume (float): in Ang^3

        Returns:
            float: unitless
        """
        if isinstance(self.eos, PolynomialEOS):
            p = np.poly1d(self.eos.eos_params)
            # first derivative of energy at 0K w.r.t. volume evaluated at the
            # given volume, in eV/Ang^3
            dEdV = np.polyder(p, 1)(volume)
            # second derivative of energy at 0K w.r.t. volume evaluated at the
            # given volume, in eV/Ang^6
            d2EdV2 = np.polyder(p, 2)(volume)
            # third derivative of energy at 0K w.r.t. volume evaluated at the
            # given volume, in eV/Ang^9
            d3EdV3 = np.polyder(p, 3)(volume)
        else:

            def central_diff_weights(Np, ndiv):
                """
                Notes:
                    This function is taken from SciPy's deprecated implementation of
                    central difference weights. The original code can be found at:
                    https://github.com/scipy/scipy/blob/ea916c6f7f487bd53e98de08264
                    9d542cc6106ed/scipy/_lib/_finite_differences.py

                    License: This code is distributed under the BSD 3-Clause license.

                    Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
                    All rights reserved.

                Return weights for an Np-point central derivative.

                Assumes equally-spaced function points.

                If weights are in the vector w, then
                derivative is w[0] * f(x-ho*dx) + ... + w[-1] * f(x+h0*dx)

                Args:
                    Np (int): Number of points for the central derivative.
                    ndiv (int, optional): Number of divisions. Default is 1.

                Returns:
                    w (ndarray): Weights for an Np-point central derivative.
                        Its size is `Np`.
                """

                ho = Np >> 1
                x = np.arange(-ho, ho + 1.0)
                x = x[:, np.newaxis]
                X = x**0.0
                for k in range(1, Np):
                    X = np.hstack([X, x**k])
                return np.prod(np.arange(1, ndiv + 1), axis=0) * linalg.inv(X)[ndiv]

            def derivative(func, x0, dx=1.0, n=1, args=(), order=3):
                """
                Notes:
                    This function is taken from SciPy's deprecated implementation of
                    central difference weights. The original code can be found at:
                    https://github.com/scipy/scipy/blob/ea916c6f7f487bd53e98de08264
                    9d542cc6106ed/scipy/_lib/_finite_differences.py

                    License: This code is distributed under the BSD 3-Clause license.

                    Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
                    All rights reserved.

                Find the nth derivative of a function at a point.

                Given a function, use a central difference formula with spacing `dx` to
                compute the nth derivative at `x0`.

                Args:
                    func (function): Input function.
                    x0 (float): The point at which the nth derivative is found.
                    dx (float, optional): Spacing.
                    n (int, optional): Order of the derivative. Default is 1.
                    args (tuple, optional): Arguments
                    order (int, optional): Number of points to use, must be odd.
                """
                if order < n + 1:
                    raise ValueError(
                        "'order' (the number of points used to compute the derivative), "
                        "must be at least the derivative order 'n' + 1."
                    )
                if order % 2 == 0:
                    raise ValueError("'order' (the number of points used to compute the derivative) must be odd.")
                # pre-computed for n=1 and 2 and low-order for speed.
                if n == 1:
                    if order == 3:
                        weights = np.array([-1, 0, 1]) / 2.0
                    elif order == 5:
                        weights = np.array([1, -8, 0, 8, -1]) / 12.0
                    elif order == 7:
                        weights = np.array([-1, 9, -45, 0, 45, -9, 1]) / 60.0
                    elif order == 9:
                        weights = np.array([3, -32, 168, -672, 0, 672, -168, 32, -3]) / 840.0
                    else:
                        weights = central_diff_weights(order, 1)
                elif n == 2:
                    if order == 3:
                        weights = np.array([1, -2.0, 1])
                    elif order == 5:
                        weights = np.array([-1, 16, -30, 16, -1]) / 12.0
                    elif order == 7:
                        weights = np.array([2, -27, 270, -490, 270, -27, 2]) / 180.0
                    elif order == 9:
                        weights = np.array([-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9]) / 5040.0
                    else:
                        weights = central_diff_weights(order, 2)
                else:
                    weights = central_diff_weights(order, n)
                val = 0.0
                ho = order >> 1
                for k in range(order):
                    val += weights[k] * func(x0 + (k - ho) * dx, *args)
                return val / np.prod((dx,) * n, axis=0)

            func = self.ev_eos_fit.func
            dEdV = derivative(func, volume, dx=1e-3, n=1)
            d2EdV2 = derivative(func, volume, dx=1e-3, n=2, order=5)
            d3EdV3 = derivative(func, volume, dx=1e-3, n=3, order=7)

        # Mie-gruneisen formulation
        if self.use_mie_gruneisen:
            p0 = dEdV
            return (
                self.gpa_to_ev_ang
                * volume
                * (self.pressure + p0 / self.gpa_to_ev_ang)
                / self.vibrational_internal_energy(temperature, volume)
            )

        # Slater-gamma formulation
        # first derivative of bulk modulus w.r.t. volume, eV/Ang^6
        dBdV = d2EdV2 + d3EdV3 * volume
        return -(1.0 / 6.0 + 0.5 * volume * dBdV / FloatWithUnit(self.ev_eos_fit.b0_GPa, "GPa").to("eV ang^-3"))

    def thermal_conductivity(self, temperature: float, volume: float) -> float:
        """Eq(17) in 10.1103/PhysRevB.90.174107.

        Args:
            temperature (float): temperature in K
            volume (float): in Ang^3

        Returns:
            float: thermal conductivity in W/K/m
        """
        gamma = self.gruneisen_parameter(temperature, volume)
        theta_d = self.debye_temperature(volume)  # K
        theta_a = theta_d * self.natoms ** (-1 / 3)  # K
        prefactor = (0.849 * 3 * 4 ** (1 / 3)) / (20.0 * np.pi**3)
        # kg/K^3/s^3
        prefactor = prefactor * (self.kb / self.hbar) ** 3 * self.avg_mass
        kappa = prefactor / (gamma**2 - 0.514 * gamma + 0.228)
        # kg/K/s^3 * Ang = (kg m/s^2)/(Ks)*1e-10
        # = N/(Ks)*1e-10 = Nm/(Kms)*1e-10 = W/K/m*1e-10
        return kappa * theta_a**2 * volume ** (1 / 3) * 1e-10

    def get_summary_dict(self):
        """Get a dict with a summary of the computed properties."""
        dct = defaultdict(list)
        dct["pressure"] = self.pressure
        dct["poisson"] = self.poisson
        dct["mass"] = self.mass
        dct["natoms"] = int(self.natoms)
        dct["bulk_modulus"] = self.bulk_modulus
        dct["gibbs_free_energy"] = self.gibbs_free_energy
        dct["temperatures"] = self.temperatures
        dct["optimum_volumes"] = self.optimum_volumes
        for v, t in zip(self.optimum_volumes, self.temperatures, strict=True):
            dct["debye_temperature"].append(self.debye_temperature(v))
            dct["gruneisen_parameter"].append(self.gruneisen_parameter(t, v))
            dct["thermal_conductivity"].append(self.thermal_conductivity(t, v))
        return dct


@deprecated(
    QuasiHarmonicDebyeApprox,
    message="Deprecated on 2024-03-27.",
    deadline=(2025, 3, 27),
)
class QuasiharmonicDebyeApprox(QuasiHarmonicDebyeApprox):
    pass
