#!/usr/bin/env python

"""
A module to perform diffusion analyses (e.g. calculating diffusivity from
mean square displacements etc.). If you use this module, please consider
citing the following papers::

    Ong, S. P., Mo, Y., Richards, W. D., Miara, L., Lee, H. S., & Ceder, G.
    (2013). Phase stability, electrochemical stability and ionic conductivity
    of the Li10+-1MP2X12 (M = Ge, Si, Sn, Al or P, and X = O, S or Se) family
    of superionic conductors. Energy & Environmental Science, 6(1), 148.
    doi:10.1039/c2ee23355j

    Mo, Y., Ong, S. P., & Ceder, G. (2012). First Principles Study of the
    Li10GeP2S12 Lithium Super Ionic Conductor Material. Chemistry of Materials,
    24(1), 15-17. doi:10.1021/cm203303y
"""

from __future__ import division

__author__ = "Will Richards"
__version__ = "0.1"
__maintainer__ = "Will Richards"
__email__ = "wrichard@mit.edu"
__status__ = "Beta"
__date__ = "5/2/13"

import math

import numpy as np

from pymatgen.core import Structure, smart_element_or_specie
import pymatgen.core.physical_constants as phyc
from pymatgen.serializers.json_coders import MSONable
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.util.coord_utils import pbc_diff


class DiffusionAnalyzer(MSONable):
    """
    Class for performing diffusion analysis.

    .. attribute: diffusivity

        Diffusivity in cm^2 / cm

    .. attribute: conductivity

        Conductivity in mS / cm

    .. attribute: diffusivity_components

        A vector with diffusivity in the a, b and c directions in cm^2 / cm

    .. attribute: conductivity components

        A vector with conductivity in the a, b and c directions in mS / cm

    .. attribute: max_framework_displacement

        The maximum (drift adjusted) distance of any framework atom from its
        starting location in A

    """

    def __init__(self, structure, displacements, specie, temperature,
                 time_step, step_skip=10, min_obs=30, weighted=True):
        """
        This constructor is meant to be used with pre-processed data.
        Other convenient constructors are provided as class methods (see
        from_vaspruns and from_files).

        Given a matrix of displacements (see arguments below for expected
        format), the diffusivity is given by:

        D = 1 / 2dt * <mean square displacement>

        where d is the dimensionality, t is the time. To obtain a reliable
        diffusion estimate, a least squares regression of the MSD against
        time to obtain the slope, which is then related to the diffusivity.

        Args:
            structure:
                Initial structure.
            displacements:
                Numpy array of with shape [site, time step, axis]
            specie:
                Specie to calculate diffusivity for as a String. E.g., "Li".
            temperature:
                Temperature of the diffusion run in Kelvin.
            time_step:
                Time step between measurements.
            step_skip:
                Sampling frequency of the displacements (time_step is
                multiplied by this number to get the real time between
                measurements)
            min_obs:
                Minimum number of observations to have before including in
                the MSD vs dt calculation. E.g. If a structure has 10
                diffusing atoms, and min_obs = 30, the MSD vs dt will be
                calculated up to dt = total_run_time / 3, so that each
                diffusing atom is measured at least 3 uncorrelated times.
            weighted:
                Uses a weighted least squares to fit the MSD vs dt. Weights are
                proportional to 1/dt, since the number of observations are
                also proportional to 1/dt (and hence the variance is
                proportional to dt)
        """
        self.s = structure
        self.disp = displacements
        self.sp = specie
        self.temperature = temperature
        self.time_step = time_step
        self.step_skip = step_skip
        self.min_obs = min_obs
        self.weighted = weighted

        self.indices = []
        self.framework_indices = []

        for i, site in enumerate(structure):
            if site.specie.symbol == specie:
                self.indices.append(i)
            else:
                self.framework_indices.append(i)
        if self.disp.shape[1] < 2:
            self.diffusivity = 0.
            self.conductivity = 0.
            self.diffusivity_components = np.array([0., 0., 0.])
            self.conductivity_components = np.array([0., 0., 0.])
            self.max_framework_displacement = 0
        else:
            framework_disp = self.disp[self.framework_indices]
            drift = np.average(framework_disp, axis=0)[None, :, :]

            #drift corrected position
            dc_x = self.disp[self.indices] - drift
            dc_framework = self.disp[self.framework_indices] - drift
            self.max_framework_displacement = \
                np.max(np.sum(dc_framework ** 2, axis=-1) ** 0.5)
            df_x = self.s.lattice.get_fractional_coords(dc_x)

            #limit the number of sampled timesteps to 200
            min_dt = int(1000 / (self.step_skip * self.time_step))
            max_dt = min(dc_x.shape[0] * dc_x.shape[1] // self.min_obs,
                         dc_x.shape[1])
            if min_dt >= max_dt:
                raise ValueError('Not enough data to calculate diffusivity')
            timesteps = np.arange(min_dt, max_dt,
                                  max(int((max_dt - min_dt) / 200), 1))
            self.dt = timesteps * self.time_step * self.step_skip

            #calculate the smoothed msd values
            self.s_msd = np.zeros_like(self.dt, dtype=np.double)
            self.s_msd_components = np.zeros(self.dt.shape + (3,))
            lengths = np.array(self.s.lattice.abc)[None, None, :]
            for i, n in enumerate(timesteps):
                dx = dc_x[:, n:, :] - dc_x[:, :-n, :]
                self.s_msd[i] = 3 * np.average(dx ** 2)
                dcomponents = (df_x[:, n:, :] - df_x[:, :-n, :]) * lengths
                self.s_msd_components[i] = np.average(
                    np.average(dcomponents ** 2, axis=1), axis=0)

            #run the regression on the msd components
            if weighted:
                w = 1 / self.dt
            else:
                w = np.ones_like(self.dt)

            #weighted least squares
            def weighted_lstsq(a, b, w):
                w_root = w ** 0.5
                x, res, rank, s = np.linalg.lstsq(a * w_root[:, None],
                                                  b * w_root)
                return x

            m_components = np.zeros(3)
            a = np.ones((len(self.dt), 2))
            a[:, 0] = self.dt
            for i in range(3):
                (m, c) = weighted_lstsq(a, self.s_msd_components[:, i], w)
                m_components[i] = max(m, 1e-15)

            (m, c) = weighted_lstsq(a, self.s_msd, w)
            #m shouldn't be negative
            m = max(m, 1e-15)

            #factor of 10 is to convert from A^2/fs to cm^2/s
            #factor of 6 is for dimensionality
            conv_factor = get_conversion_factor(self.s, self.sp,
                                                self.temperature)
            self.diffusivity = m / 60
            self.conductivity = self.diffusivity * conv_factor

            self.diffusivity_components = m_components / 20
            self.conductivity_components = self.diffusivity_components * \
                conv_factor

    def plot_smoothed_msd(self):
        """
        Plot the smoothed msd vs time graph. Useful for checking convergence.
        """
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)
        plt.plot(self.dt, self.s_msd, 'k')
        plt.plot(self.dt, self.s_msd_components[:, 0], 'r')
        plt.plot(self.dt, self.s_msd_components[:, 1], 'g')
        plt.plot(self.dt, self.s_msd_components[:, 2], 'b')
        plt.legend(["Overall", "a", "b", "c"], loc=2, prop={"size": 20})
        plt.xlabel("Timestep")
        plt.ylabel("MSD")
        plt.tight_layout()
        plt.show()

    @classmethod
    def from_vaspruns(cls, vaspruns, specie, min_obs=30, weighted=True):
        """
        Convenient constructor that takes in a list of Vasprun objects to
        perform diffusion analysis.

        Args:
            vaspruns:
                list of Vasprun objects (must be ordered in sequence of run).
                 E.g., you may have performed sequential VASP runs to obtain
                 sufficient statistics.
            specie:
                Specie to calculate diffusivity for as a String. E.g., "Li".
            min_obs:
                Minimum number of observations to have before including in
                the MSD vs dt calculation. E.g. If a structure has 10
                diffusing atoms, and min_obs = 30, the MSD vs dt will be
                calculated up to dt = total_run_time / 3, so that each
                diffusing atom is measured at least 3 uncorrelated times.
            weighted:
                Uses a weighted least squares to fit the MSD vs dt. Weights are
                proportional to 1/dt, since the number of observations are
                also proportional to 1/dt (and hence the variance is
                proportional to dt)
        """
        structure = vaspruns[0].initial_structure
        step_skip = vaspruns[0].ionic_step_skip or 1

        p = []
        final_structure = vaspruns[0].initial_structure
        for vr in vaspruns:
            #check that the runs are continuous
            fdist = pbc_diff(vr.initial_structure.frac_coords, 
                             final_structure.frac_coords)
            if np.any(fdist > 0.001):
                raise ValueError('initial and final structures do not '
                                 'match.')
            final_structure = vr.final_structure
            
            assert (vr.ionic_step_skip or 1) == step_skip
            p.extend([np.array(s['structure'].frac_coords)[:, None]
                      for s in vr.ionic_steps])
        p = np.concatenate(p, axis=1)
        dp = p[:, 1:] - p[:, :-1]
        dp = np.concatenate([np.zeros_like(dp[:, (0,)]), dp], axis=1)
        dp = dp - np.round(dp)
        f_disp = np.cumsum(dp, axis=1)
        disp = structure.lattice.get_cartesian_coords(f_disp)

        temperature = vaspruns[0].parameters['TEEND']
        time_step = vaspruns[0].parameters['POTIM']

        return cls(structure, disp, specie, temperature,
                   time_step, step_skip=step_skip, min_obs=min_obs,
                   weighted=weighted)

    @classmethod
    def from_files(cls, filepaths, specie, step_skip=10, min_obs=30,
                   weighted=True, ncores=None):
        """
        Convenient constructor that takes in a list of vasprun.xml paths to
        perform diffusion analysis.

        Args:
            filepaths:
                List of paths to vasprun.xml files of runs. (must be
                ordered in sequence of run). For example,
                you may have done sequential VASP runs and they are in run1,
                run2, run3, etc. You should then pass in
                ["run1/vasprun.xml", "run2/vasprun.xml", ...].
            specie:
                Specie to calculate diffusivity for as a String. E.g., "Li".
            step_skip:
                Sampling frequency of the displacements (time_step is
                multiplied by this number to get the real time between
                measurements). E.g., you may not want to sample every single
                time step.
            min_obs:
                Minimum number of observations to have before including in
                the MSD vs dt calculation. E.g. If a structure has 10
                diffusing atoms, and min_obs = 30, the MSD vs dt will be
                calculated up to dt = total_run_time / 3, so that each
                diffusing atom is measured at least 3 uncorrelated times.
            weighted:
                Uses a weighted least squares to fit the MSD vs dt. Weights are
                proportional to 1/dt, since the number of observations are
                also proportional to 1/dt (and hence the variance is
                proportional to dt)
            ncores:
                Numbers of cores to use for multiprocessing. Can speed up
                vasprun parsing considerable. Defaults to None,
                which means serial.
        """
        func = map
        if ncores is not None:
            import multiprocessing
            p = multiprocessing.Pool(ncores)
            func = p.map
        vaspruns = func(_get_vasprun, [(p, step_skip) for p in filepaths])
        return cls.from_vaspruns(vaspruns, min_obs=min_obs,
                                 weighted=weighted, specie=specie)

    @property
    def to_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "structure": self.s.to_dict,
            "displacements": self.disp.tolist(),
            "specie": self.sp,
            "temperature": self.temperature,
            "time_step": self.time_step,
            "step_skip": self.step_skip,
            "min_obs": self.min_obs,
            "weighted": self.weighted
        }

    @classmethod
    def from_dict(cls, d):
        structure = Structure.from_dict(d["structure"])
        return cls(structure, np.array(d["displacements"]), specie=d["specie"],
                   temperature=d["temperature"], time_step=d["time_step"],
                   step_skip=d["step_skip"], min_obs=d["min_obs"],
                   weighted=d["weighted"])


def get_conversion_factor(structure, species, temperature):
    """
    Conversion factor to convert between cm^2/s diffusivity measurements and
    mS/cm conductivity measurements based on number of atoms of diffusing
    species. Note that the charge is based on the oxidation state of the
    species (where available), or else the number of valence electrons
    (usually a good guess, esp for main group ions).

    Args:
        structure:
            Input structure.
        species:
            Diffusing species.
        temperature:
            Temperature of the diffusion run in Kelvin.

    Returns:
        Conversion factor.
        Conductivity (in mS/cm) = Conversion Factor * Diffusivity (in cm^2/s)
    """
    df_sp = smart_element_or_specie(species)
    if hasattr(df_sp, "oxi_state"):
        z = df_sp.oxi_state
    else:
        z = df_sp.full_electronic_structure[-1][2]

    n = structure.composition[species]

    V = structure.volume * 1e-24  # units cm^3
    return 1000 * n / (V * phyc.N_a) * z ** 2 * phyc.F ** 2\
        / (phyc.R * temperature)


def _get_vasprun(args):
    """
    Internal method to support multiprocessing.
    """
    return Vasprun(args[0], ionic_step_skip=args[1])


def get_arrhenius_plot(temps, diffusivites, **kwargs):
    """
    Returns an Arrhenius plot.

    Args:
        temps:
            A sequence of temperatures.
        diffusivities:
            A sequence of diffusivities (e.g., from DiffusionAnalyzer
            .diffusivity).
        \*\*kwargs:
            Any keyword args supported by matplotlib.pyplot.plot.

    Returns:
        A matplotlib.pyplot object. Do plt.show() to show the plot.
    """
    t_1 = 1000 / np.array(temps)
    logd = np.log10(diffusivites)
    #Do a least squares regression of log(D) vs 1000/T
    A = np.array([t_1, np.ones(len(temps))]).T
    w = np.array(np.linalg.lstsq(A, logd)[0])
    from pymatgen.util.plotting_utils import get_publication_quality_plot
    plt = get_publication_quality_plot(12, 8)
    plt.plot(t_1, logd, 'ko', t_1, np.dot(A, w), 'k--', markersize=10,
             **kwargs)
    # Calculate the activation energy in meV = negative of the slope,
    # * kB (/ electron charge to convert to eV), * 1000 (inv. temperature
    # scale), * 1000 (eV -> meV), * math.log(10) (the regression is carried
    # out in base 10 for easier reading of the diffusivity scale,
    # but the Arrhenius relationship is in base e).
    actv_energy = - w[0] * phyc.k_b / phyc.e * 1e6 * math.log(10)
    plt.annotate("E$_a$ = {:.0f} meV".format(actv_energy),
                 (t_1[-1], w[0] * t_1[-1] + w[1]), xytext=(100, 0),
                 xycoords='data', textcoords='offset points', fontsize=30)
    plt.ylabel("log(D (cm$^2$/s))")
    plt.xlabel("1000/T (K$^{-1}$)")
    plt.tight_layout()
    return plt
