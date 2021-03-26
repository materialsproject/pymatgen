# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

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


import multiprocessing
import warnings

import numpy as np
import scipy.constants as const
from monty.json import MSONable

from pymatgen.analysis.structure_matcher import (
    OrderDisorderElementComparator,
    StructureMatcher,
)
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.util.coord import pbc_diff

__author__ = "Will Richards, Shyue Ping Ong"
__version__ = "0.2"
__maintainer__ = "Will Richards"
__email__ = "wrichard@mit.edu"
__status__ = "Beta"
__date__ = "5/2/13"

warnings.warn(
    "All code in pymatgen.analysis.diffusion_analyzer has been moved to the separate add-on package"
    "pymatgen-diffusion, which also includes a lot more functionality for analyzing diffusion"
    "calculations. This module here is retained for backwards compatibility. It will be removed from"
    "2022.1.1.",
    FutureWarning,
)


class DiffusionAnalyzer(MSONable):
    """
    Class for performing diffusion analysis.

    .. attribute: diffusivity

        Diffusivity in cm^2 / s

    .. attribute: chg_diffusivity

        Charge diffusivity in cm^2 / s

    .. attribute: conductivity

        Conductivity in mS / cm

    .. attribute: chg_conductivity

        Conductivity derived from Nernst-Einstein equation using charge
        diffusivity, in mS / cm

    .. attribute: diffusivity_components

        A vector with diffusivity in the a, b and c directions in cm^2 / s

    .. attribute: conductivity_components

        A vector with conductivity in the a, b and c directions in mS / cm

    .. attribute: diffusivity_std_dev

        Std dev in diffusivity in cm^2 / s. Note that this makes sense only
        for non-smoothed analyses.

    .. attribute: chg_diffusivity_std_dev

        Std dev in charge diffusivity in cm^2 / s. Note that this makes sense only
        for non-smoothed analyses.

    .. attribute: conductivity_std_dev

        Std dev in conductivity in mS / cm. Note that this makes sense only
        for non-smoothed analyses.

    .. attribute: diffusivity_components_std_dev

        A vector with std dev. in diffusivity in the a, b and c directions in
        cm^2 / cm. Note that this makes sense only for non-smoothed analyses.

    .. attribute: conductivity_components_std_dev

        A vector with std dev. in conductivity in the a, b and c directions
        in mS / cm. Note that this makes sense only for non-smoothed analyses.

    .. attribute: max_framework_displacement

        The maximum (drift adjusted) distance of any framework atom from its
        starting location in A.

    .. attribute: max_ion_displacements

        nions x 1 array of the maximum displacement of each individual ion.

    .. attribute: msd

        nsteps x 1 array of the mean square displacement of specie.

    .. attribute: mscd

        nsteps x 1 array of the mean square charge displacement of specie.

    .. attribute: msd_components

        nsteps x 3 array of the MSD in each lattice direction of specie.

    .. attribute: sq_disp_ions

        The square displacement of all ion (both specie and other ions) as a
        nions x nsteps array.

    .. attribute: dt

        Time coordinate array.

    .. attribute: haven_ratio
        Haven ratio defined as diffusivity / chg_diffusivity.
    """

    def __init__(
        self,
        structure,
        displacements,
        specie,
        temperature,
        time_step,
        step_skip,
        smoothed="max",
        min_obs=30,
        avg_nsteps=1000,
        lattices=None,
    ):
        """
        This constructor is meant to be used with pre-processed data.
        Other convenient constructors are provided as class methods (see
        from_vaspruns and from_files).

        Given a matrix of displacements (see arguments below for expected
        format), the diffusivity is given by::

            D = 1 / 2dt * <mean square displacement>

        where d is the dimensionality, t is the time. To obtain a reliable
        diffusion estimate, a least squares regression of the MSD against
        time to obtain the slope, which is then related to the diffusivity.

        For traditional analysis, use smoothed=False and weighted=False.

        Args:
            structure (Structure): Initial structure.
            displacements (array): Numpy array of with shape [site,
                time step, axis]
            specie (Element/Species): Species to calculate diffusivity for as a
                String. E.g., "Li".
            temperature (float): Temperature of the diffusion run in Kelvin.
            time_step (int): Time step between measurements.
            step_skip (int): Sampling frequency of the displacements (
                time_step is multiplied by this number to get the real time
                between measurements)
            smoothed (str): Whether to smooth the MSD, and what mode to smooth.
                Supported modes are:

                i. "max", which tries to use the maximum #
                   of data points for each time origin, subject to a
                   minimum # of observations given by min_obs, and then
                   weights the observations based on the variance
                   accordingly. This is the default.
                ii. "constant", in which each timestep is averaged over
                    the number of time_steps given by min_steps.
                iii. None / False / any other false-like quantity. No
                   smoothing.

            min_obs (int): Used with smoothed="max". Minimum number of
                observations to have before including in the MSD vs dt
                calculation. E.g. If a structure has 10 diffusing atoms,
                and min_obs = 30, the MSD vs dt will be
                calculated up to dt = total_run_time / 3, so that each
                diffusing atom is measured at least 3 uncorrelated times.
                Only applies in smoothed="max".
            avg_nsteps (int): Used with smoothed="constant". Determines the
                number of time steps to average over to get the msd for each
                timestep. Default of 1000 is usually pretty good.
            lattices (array): Numpy array of lattice matrix of every step. Used
                for NPT-AIMD. For NVT-AIMD, the lattice at each time step is
                set to the lattice in the "structure" argument.
        """
        self.structure = structure
        self.disp = displacements
        self.specie = specie
        self.temperature = temperature
        self.time_step = time_step
        self.step_skip = step_skip
        self.min_obs = min_obs
        self.smoothed = smoothed
        self.avg_nsteps = avg_nsteps
        self.lattices = lattices

        if lattices is None:
            self.lattices = np.array([structure.lattice.matrix.tolist()])

        indices = []
        framework_indices = []
        for i, site in enumerate(structure):
            if site.specie.symbol == specie:
                indices.append(i)
            else:
                framework_indices.append(i)
        if self.disp.shape[1] < 2:
            self.diffusivity = 0.0
            self.conductivity = 0.0
            self.diffusivity_components = np.array([0.0, 0.0, 0.0])
            self.conductivity_components = np.array([0.0, 0.0, 0.0])
            self.max_framework_displacement = 0
        else:
            framework_disp = self.disp[framework_indices]
            drift = np.average(framework_disp, axis=0)[None, :, :]

            # drift corrected position
            dc = self.disp - drift

            nions, nsteps, dim = dc.shape

            if not smoothed:
                timesteps = np.arange(0, nsteps)
            elif smoothed == "constant":
                if nsteps <= avg_nsteps:
                    raise ValueError("Not enough data to calculate diffusivity")
                timesteps = np.arange(0, nsteps - avg_nsteps)
            else:
                # limit the number of sampled timesteps to 200
                min_dt = int(1000 / (self.step_skip * self.time_step))
                max_dt = min(len(indices) * nsteps // self.min_obs, nsteps)
                if min_dt >= max_dt:
                    raise ValueError("Not enough data to calculate diffusivity")
                timesteps = np.arange(min_dt, max_dt, max(int((max_dt - min_dt) / 200), 1))

            dt = timesteps * self.time_step * self.step_skip

            # calculate the smoothed msd values
            msd = np.zeros_like(dt, dtype=np.double)
            sq_disp_ions = np.zeros((len(dc), len(dt)), dtype=np.double)
            msd_components = np.zeros(dt.shape + (3,))

            # calculate mean square charge displacement
            mscd = np.zeros_like(msd, dtype=np.double)

            for i, n in enumerate(timesteps):
                if not smoothed:
                    dx = dc[:, i : i + 1, :]
                    dcomponents = dc[:, i : i + 1, :]
                elif smoothed == "constant":
                    dx = dc[:, i : i + avg_nsteps, :] - dc[:, 0:avg_nsteps, :]
                    dcomponents = dc[:, i : i + avg_nsteps, :] - dc[:, 0:avg_nsteps, :]
                else:
                    dx = dc[:, n:, :] - dc[:, :-n, :]
                    dcomponents = dc[:, n:, :] - dc[:, :-n, :]

                # Get msd
                sq_disp = dx ** 2
                sq_disp_ions[:, i] = np.average(np.sum(sq_disp, axis=2), axis=1)
                msd[i] = np.average(sq_disp_ions[:, i][indices])

                msd_components[i] = np.average(dcomponents[indices] ** 2, axis=(0, 1))

                # Get mscd
                sq_chg_disp = np.sum(dx[indices, :, :], axis=0) ** 2
                mscd[i] = np.average(np.sum(sq_chg_disp, axis=1), axis=0) / len(indices)

            def weighted_lstsq(a, b):
                if smoothed == "max":
                    # For max smoothing, we need to weight by variance.
                    w_root = (1 / dt) ** 0.5
                    return np.linalg.lstsq(a * w_root[:, None], b * w_root, rcond=None)
                return np.linalg.lstsq(a, b, rcond=None)

            # Get self diffusivity
            m_components = np.zeros(3)
            m_components_res = np.zeros(3)
            a = np.ones((len(dt), 2))
            a[:, 0] = dt
            for i in range(3):
                (m, c), res, rank, s = weighted_lstsq(a, msd_components[:, i])
                m_components[i] = max(m, 1e-15)
                m_components_res[i] = res[0]

            (m, c), res, rank, s = weighted_lstsq(a, msd)
            # m shouldn't be negative
            m = max(m, 1e-15)

            # Get also the charge diffusivity
            (m_chg, c_chg), res_chg, _, _ = weighted_lstsq(a, mscd)
            # m shouldn't be negative
            m_chg = max(m_chg, 1e-15)

            # factor of 10 is to convert from A^2/fs to cm^2/s
            # factor of 6 is for dimensionality
            conv_factor = get_conversion_factor(self.structure, self.specie, self.temperature)
            self.diffusivity = m / 60
            self.chg_diffusivity = m_chg / 60

            # Calculate the error in the diffusivity using the error in the
            # slope from the lst sq.
            # Variance in slope = n * Sum Squared Residuals / (n * Sxx - Sx
            # ** 2) / (n-2).
            n = len(dt)

            # Pre-compute the denominator since we will use it later.
            # We divide dt by 1000 to avoid overflow errors in some systems (
            # e.g., win). This is subsequently corrected where denom is used.
            denom = (n * np.sum((dt / 1000) ** 2) - np.sum(dt / 1000) ** 2) * (n - 2)
            self.diffusivity_std_dev = np.sqrt(n * res[0] / denom) / 60 / 1000
            self.chg_diffusivity_std_dev = np.sqrt(n * res_chg[0] / denom) / 60 / 1000
            self.conductivity = self.diffusivity * conv_factor
            self.chg_conductivity = self.chg_diffusivity * conv_factor
            self.conductivity_std_dev = self.diffusivity_std_dev * conv_factor

            self.diffusivity_components = m_components / 20
            self.diffusivity_components_std_dev = np.sqrt(n * m_components_res / denom) / 20 / 1000
            self.conductivity_components = self.diffusivity_components * conv_factor
            self.conductivity_components_std_dev = self.diffusivity_components_std_dev * conv_factor

            # Drift and displacement information.
            self.drift = drift
            self.corrected_displacements = dc
            self.max_ion_displacements = np.max(np.sum(dc ** 2, axis=-1) ** 0.5, axis=1)
            self.max_framework_displacement = np.max(self.max_ion_displacements[framework_indices])
            self.msd = msd
            self.mscd = mscd
            self.haven_ratio = self.diffusivity / self.chg_diffusivity
            self.sq_disp_ions = sq_disp_ions
            self.msd_components = msd_components
            self.dt = dt
            self.indices = indices
            self.framework_indices = framework_indices

    def get_drift_corrected_structures(self, start=None, stop=None, step=None):
        """
        Returns an iterator for the drift-corrected structures. Use of
        iterator is to reduce memory usage as # of structures in MD can be
        huge. You don't often need all the structures all at once.

        Args:
            start, stop, step (int): applies a start/stop/step to the iterator.
                Faster than applying it after generation, as it reduces the
                number of structures created.
        """
        coords = np.array(self.structure.cart_coords)
        species = self.structure.species_and_occu
        lattices = self.lattices
        nsites, nsteps, dim = self.corrected_displacements.shape

        for i in range(start or 0, stop or nsteps, step or 1):
            latt = lattices[0] if len(lattices) == 1 else lattices[i]
            yield Structure(
                latt,
                species,
                coords + self.corrected_displacements[:, i, :],
                coords_are_cartesian=True,
            )

    def get_summary_dict(self, include_msd_t=False, include_mscd_t=False):
        """
        Provides a summary of diffusion information.

        Args:
            include_msd_t (bool): Whether to include mean square displace and
                time data with the data.
            include_msd_t (bool): Whether to include mean square charge displace and
                time data with the data.

        Returns:
            (dict) of diffusion and conductivity data.
        """
        d = {
            "D": self.diffusivity,
            "D_sigma": self.diffusivity_std_dev,
            "D_charge": self.chg_diffusivity,
            "D_charge_sigma": self.chg_diffusivity_std_dev,
            "S": self.conductivity,
            "S_sigma": self.conductivity_std_dev,
            "S_charge": self.chg_conductivity,
            "D_components": self.diffusivity_components.tolist(),
            "S_components": self.conductivity_components.tolist(),
            "D_components_sigma": self.diffusivity_components_std_dev.tolist(),
            "S_components_sigma": self.conductivity_components_std_dev.tolist(),
            "specie": str(self.specie),
            "step_skip": self.step_skip,
            "time_step": self.time_step,
            "temperature": self.temperature,
            "max_framework_displacement": self.max_framework_displacement,
            "Haven_ratio": self.haven_ratio,
        }
        if include_msd_t:
            d["msd"] = self.msd.tolist()
            d["msd_components"] = self.msd_components.tolist()
            d["dt"] = self.dt.tolist()
        if include_mscd_t:
            d["mscd"] = self.mscd.tolist()
        return d

    def get_framework_rms_plot(self, plt=None, granularity=200, matching_s=None):
        """
        Get the plot of rms framework displacement vs time. Useful for checking
        for melting, especially if framework atoms can move via paddle-wheel
        or similar mechanism (which would show up in max framework displacement
        but doesn't constitute melting).

        Args:
            plt (matplotlib.pyplot): If plt is supplied, changes will be made
                to an existing plot. Otherwise, a new plot will be created.
            granularity (int): Number of structures to match
            matching_s (Structure): Optionally match to a disordered structure
                instead of the first structure in the analyzer. Required when
                a secondary mobile ion is present.
        Notes:
            The method doesn't apply to NPT-AIMD simulation analysis.
        """
        from pymatgen.util.plotting import pretty_plot

        if self.lattices is not None and len(self.lattices) > 1:
            warnings.warn("Note the method doesn't apply to NPT-AIMD " "simulation analysis!")

        plt = pretty_plot(12, 8, plt=plt)
        step = (self.corrected_displacements.shape[1] - 1) // (granularity - 1)
        f = (matching_s or self.structure).copy()
        f.remove_species([self.specie])
        sm = StructureMatcher(
            primitive_cell=False,
            stol=0.6,
            comparator=OrderDisorderElementComparator(),
            allow_subset=True,
        )
        rms = []
        for s in self.get_drift_corrected_structures(step=step):
            s.remove_species([self.specie])
            d = sm.get_rms_dist(f, s)
            if d:
                rms.append(d)
            else:
                rms.append((1, 1))
        max_dt = (len(rms) - 1) * step * self.step_skip * self.time_step
        if max_dt > 100000:
            plot_dt = np.linspace(0, max_dt / 1000, len(rms))
            unit = "ps"
        else:
            plot_dt = np.linspace(0, max_dt, len(rms))
            unit = "fs"
        rms = np.array(rms)
        plt.plot(plot_dt, rms[:, 0], label="RMS")
        plt.plot(plot_dt, rms[:, 1], label="max")
        plt.legend(loc="best")
        plt.xlabel("Timestep ({})".format(unit))
        plt.ylabel("normalized distance")
        plt.tight_layout()
        return plt

    def get_msd_plot(self, plt=None, mode="specie"):
        """
        Get the plot of the smoothed msd vs time graph. Useful for
        checking convergence. This can be written to an image file.

        Args:
            plt: A plot object. Defaults to None, which means one will be
                generated.
            mode (str): Determines type of msd plot. By "species", "sites",
                or direction (default). If mode = "mscd", the smoothed mscd vs.
                time will be plotted.
        """
        from pymatgen.util.plotting import pretty_plot

        plt = pretty_plot(12, 8, plt=plt)
        if np.max(self.dt) > 100000:
            plot_dt = self.dt / 1000
            unit = "ps"
        else:
            plot_dt = self.dt
            unit = "fs"

        if mode == "species":
            for sp in sorted(self.structure.composition.keys()):
                indices = [i for i, site in enumerate(self.structure) if site.specie == sp]
                sd = np.average(self.sq_disp_ions[indices, :], axis=0)
                plt.plot(plot_dt, sd, label=sp.__str__())
            plt.legend(loc=2, prop={"size": 20})
        elif mode == "sites":
            for i, site in enumerate(self.structure):
                sd = self.sq_disp_ions[i, :]
                plt.plot(plot_dt, sd, label="%s - %d" % (site.specie.__str__(), i))
            plt.legend(loc=2, prop={"size": 20})
        elif mode == "mscd":
            plt.plot(plot_dt, self.mscd, "r")
            plt.legend(["Overall"], loc=2, prop={"size": 20})
        else:
            # Handle default / invalid mode case
            plt.plot(plot_dt, self.msd, "k")
            plt.plot(plot_dt, self.msd_components[:, 0], "r")
            plt.plot(plot_dt, self.msd_components[:, 1], "g")
            plt.plot(plot_dt, self.msd_components[:, 2], "b")
            plt.legend(["Overall", "a", "b", "c"], loc=2, prop={"size": 20})

        plt.xlabel("Timestep ({})".format(unit))
        if mode == "mscd":
            plt.ylabel("MSCD ($\\AA^2$)")
        else:
            plt.ylabel("MSD ($\\AA^2$)")
        plt.tight_layout()
        return plt

    def plot_msd(self, mode="default"):
        """
        Plot the smoothed msd vs time graph. Useful for checking convergence.

        Args:
            mode (str): Can be "default" (the default, shows only the MSD for
                the diffusing specie, and its components), "ions" (individual
                square displacements of all ions), "species" (mean square
                displacement by specie), or "mscd" (overall mean square charge
                displacement for diffusing specie).
        """
        self.get_msd_plot(mode=mode).show()

    def export_msdt(self, filename):
        """
        Writes MSD data to a csv file that can be easily plotted in other
        software.

        Args:
            filename (str): Filename. Supported formats are csv and dat. If
                the extension is csv, a csv file is written. Otherwise,
                a dat format is assumed.
        """
        fmt = "csv" if filename.lower().endswith(".csv") else "dat"
        delimiter = ", " if fmt == "csv" else " "
        with open(filename, "wt") as f:
            if fmt == "dat":
                f.write("# ")
            f.write(delimiter.join(["t", "MSD", "MSD_a", "MSD_b", "MSD_c", "MSCD"]))
            f.write("\n")
            for dt, msd, msdc, mscd in zip(self.dt, self.msd, self.msd_components, self.mscd):
                f.write(delimiter.join(["%s" % v for v in [dt, msd] + list(msdc) + [mscd]]))
                f.write("\n")

    @classmethod
    def from_structures(
        cls, structures, specie, temperature, time_step, step_skip, initial_disp=None, initial_structure=None, **kwargs
    ):
        r"""
        Convenient constructor that takes in a list of Structure objects to
        perform diffusion analysis.

        Args:
            structures ([Structure]): list of Structure objects (must be
                ordered in sequence of run). E.g., you may have performed
                sequential VASP runs to obtain sufficient statistics.
            specie (Element/Species): Species to calculate diffusivity for as a
                String. E.g., "Li".
            temperature (float): Temperature of the diffusion run in Kelvin.
            time_step (int): Time step between measurements.
            step_skip (int): Sampling frequency of the displacements (
                time_step is multiplied by this number to get the real time
                between measurements)
            initial_disp (np.ndarray): Sometimes, you need to iteratively
                compute estimates of the diffusivity. This supplies an
                initial displacement that will be added on to the initial
                displacements. Note that this makes sense only when
                smoothed=False.
            initial_structure (Structure): Like initial_disp, this is used
                for iterative computations of estimates of the diffusivity. You
                typically need to supply both variables. This stipulates the
                initial structure from which the current set of displacements
                are computed.
            \\*\\*kwargs: kwargs supported by the :class:`DiffusionAnalyzer`_.
                Examples include smoothed, min_obs, avg_nsteps.
        """
        p, l = [], []
        for i, s in enumerate(structures):
            if i == 0:
                structure = s
            p.append(np.array(s.frac_coords)[:, None])
            l.append(s.lattice.matrix)
        if initial_structure is not None:
            p.insert(0, np.array(initial_structure.frac_coords)[:, None])
            l.insert(0, initial_structure.lattice.matrix)
        else:
            p.insert(0, p[0])
            l.insert(0, l[0])

        p = np.concatenate(p, axis=1)
        dp = p[:, 1:] - p[:, :-1]
        dp = dp - np.round(dp)
        f_disp = np.cumsum(dp, axis=1)
        c_disp = []
        for i in f_disp:
            c_disp.append([np.dot(d, m) for d, m in zip(i, l[1:])])
        disp = np.array(c_disp)

        # If is NVT-AIMD, clear lattice data.
        if np.array_equal(l[0], l[-1]):
            l = np.array([l[0]])
        else:
            l = np.array(l)
        if initial_disp is not None:
            disp += initial_disp[:, None, :]

        return cls(structure, disp, specie, temperature, time_step, step_skip=step_skip, lattices=l, **kwargs)

    @classmethod
    def from_vaspruns(cls, vaspruns, specie, initial_disp=None, initial_structure=None, **kwargs):
        r"""
        Convenient constructor that takes in a list of Vasprun objects to
        perform diffusion analysis.

        Args:
            vaspruns ([Vasprun]): List of Vaspruns (must be ordered  in
                sequence of MD simulation). E.g., you may have performed
                sequential VASP runs to obtain sufficient statistics.
            specie (Element/Species): Species to calculate diffusivity for as a
                String. E.g., "Li".
            initial_disp (np.ndarray): Sometimes, you need to iteratively
                compute estimates of the diffusivity. This supplies an
                initial displacement that will be added on to the initial
                displacements. Note that this makes sense only when
                smoothed=False.
            initial_structure (Structure): Like initial_disp, this is used
                for iterative computations of estimates of the diffusivity. You
                typically need to supply both variables. This stipulates the
                initial stricture from which the current set of displacements
                are computed.
            \\*\\*kwargs: kwargs supported by the :class:`DiffusionAnalyzer`_.
                Examples include smoothed, min_obs, avg_nsteps.
        """

        def get_structures(vaspruns):
            for i, vr in enumerate(vaspruns):
                if i == 0:
                    step_skip = vr.ionic_step_skip or 1
                    final_structure = vr.initial_structure
                    temperature = vr.parameters["TEEND"]
                    time_step = vr.parameters["POTIM"]
                    yield step_skip, temperature, time_step
                # check that the runs are continuous
                fdist = pbc_diff(vr.initial_structure.frac_coords, final_structure.frac_coords)
                if np.any(fdist > 0.001):
                    raise ValueError("initial and final structures do not " "match.")
                final_structure = vr.final_structure

                assert (vr.ionic_step_skip or 1) == step_skip
                for s in vr.ionic_steps:
                    yield s["structure"]

        s = get_structures(vaspruns)
        step_skip, temperature, time_step = next(s)

        return cls.from_structures(
            structures=list(s),
            specie=specie,
            temperature=temperature,
            time_step=time_step,
            step_skip=step_skip,
            initial_disp=initial_disp,
            initial_structure=initial_structure,
            **kwargs,
        )

    @classmethod
    def from_files(
        cls, filepaths, specie, step_skip=10, ncores=None, initial_disp=None, initial_structure=None, **kwargs
    ):
        r"""
        Convenient constructor that takes in a list of vasprun.xml paths to
        perform diffusion analysis.

        Args:
            filepaths ([str]): List of paths to vasprun.xml files of runs. (
                must be ordered in sequence of MD simulation). For example,
                you may have done sequential VASP runs and they are in run1,
                run2, run3, etc. You should then pass in
                ["run1/vasprun.xml", "run2/vasprun.xml", ...].
            specie (Element/Species): Species to calculate diffusivity for as a
                String. E.g., "Li".
            step_skip (int): Sampling frequency of the displacements (
                time_step is multiplied by this number to get the real time
                between measurements)
            ncores (int): Numbers of cores to use for multiprocessing. Can
                speed up vasprun parsing considerably. Defaults to None,
                which means serial. It should be noted that if you want to
                use multiprocessing, the number of ionic steps in all vasprun
                .xml files should be a multiple of the ionic_step_skip.
                Otherwise, inconsistent results may arise. Serial mode has no
                such restrictions.
            initial_disp (np.ndarray): Sometimes, you need to iteratively
                compute estimates of the diffusivity. This supplies an
                initial displacement that will be added on to the initial
                displacements. Note that this makes sense only when
                smoothed=False.
            initial_structure (Structure): Like initial_disp, this is used
                for iterative computations of estimates of the diffusivity. You
                typically need to supply both variables. This stipulates the
                initial structure from which the current set of displacements
                are computed.
            \\*\\*kwargs: kwargs supported by the :class:`DiffusionAnalyzer`_.
                Examples include smoothed, min_obs, avg_nsteps.
        """
        if ncores is not None and len(filepaths) > 1:
            p = multiprocessing.Pool(ncores)
            vaspruns = p.imap(_get_vasprun, [(fp, step_skip) for fp in filepaths])
            analyzer = cls.from_vaspruns(
                vaspruns, specie=specie, initial_disp=initial_disp, initial_structure=initial_structure, **kwargs
            )
            p.close()
            p.join()
            return analyzer

        def vr(filepaths):
            offset = 0
            for p in filepaths:
                v = Vasprun(p, ionic_step_offset=offset, ionic_step_skip=step_skip)
                yield v
                # Recompute offset.
                offset = (-(v.nionic_steps - offset)) % step_skip

        return cls.from_vaspruns(
            vr(filepaths), specie=specie, initial_disp=initial_disp, initial_structure=initial_structure, **kwargs
        )

    def as_dict(self):
        """
        Returns: MSONable dict
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "structure": self.structure.as_dict(),
            "displacements": self.disp.tolist(),
            "specie": self.specie,
            "temperature": self.temperature,
            "time_step": self.time_step,
            "step_skip": self.step_skip,
            "min_obs": self.min_obs,
            "smoothed": self.smoothed,
            "avg_nsteps": self.avg_nsteps,
            "lattices": self.lattices.tolist(),
        }

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict representation

        Returns: DiffusionAnalyzer
        """
        structure = Structure.from_dict(d["structure"])
        return cls(
            structure,
            np.array(d["displacements"]),
            specie=d["specie"],
            temperature=d["temperature"],
            time_step=d["time_step"],
            step_skip=d["step_skip"],
            min_obs=d["min_obs"],
            smoothed=d.get("smoothed", "max"),
            avg_nsteps=d.get("avg_nsteps", 1000),
            lattices=np.array(d.get("lattices", [d["structure"]["lattice"]["matrix"]])),
        )


def get_conversion_factor(structure, species, temperature):
    """
    Conversion factor to convert between cm^2/s diffusivity measurements and
    mS/cm conductivity measurements based on number of atoms of diffusing
    species. Note that the charge is based on the oxidation state of the
    species (where available), or else the number of valence electrons
    (usually a good guess, esp for main group ions).

    Args:
        structure (Structure): Input structure.
        species (Element/Species): Diffusing species.
        temperature (float): Temperature of the diffusion run in Kelvin.

    Returns:
        Conversion factor.
        Conductivity (in mS/cm) = Conversion Factor * Diffusivity (in cm^2/s)
    """
    df_sp = get_el_sp(species)
    if hasattr(df_sp, "oxi_state"):
        z = df_sp.oxi_state
    else:
        z = df_sp.full_electronic_structure[-1][2]

    n = structure.composition[species]

    vol = structure.volume * 1e-24  # units cm^3
    return 1000 * n / (vol * const.N_A) * z ** 2 * (const.N_A * const.e) ** 2 / (const.R * temperature)


def _get_vasprun(args):
    """
    Internal method to support multiprocessing.
    """
    return Vasprun(args[0], ionic_step_skip=args[1], parse_dos=False, parse_eigen=False)


def fit_arrhenius(temps, diffusivities):
    """
    Returns Ea, c, standard error of Ea from the Arrhenius fit:
        D = c * exp(-Ea/kT)

    Args:
        temps ([float]): A sequence of temperatures. units: K
        diffusivities ([float]): A sequence of diffusivities (e.g.,
            from DiffusionAnalyzer.diffusivity). units: cm^2/s
    """
    t_1 = 1 / np.array(temps)
    logd = np.log(diffusivities)
    # Do a least squares regression of log(D) vs 1/T
    a = np.array([t_1, np.ones(len(temps))]).T
    w, res, _, _ = np.linalg.lstsq(a, logd, rcond=None)
    w = np.array(w)
    n = len(temps)
    if n > 2:
        std_Ea = (res[0] / (n - 2) / (n * np.var(t_1))) ** 0.5 * const.k / const.e
    else:
        std_Ea = None
    return -w[0] * const.k / const.e, np.exp(w[1]), std_Ea


def get_extrapolated_diffusivity(temps, diffusivities, new_temp):
    """
    Returns (Arrhenius) extrapolated diffusivity at new_temp

    Args:
        temps ([float]): A sequence of temperatures. units: K
        diffusivities ([float]): A sequence of diffusivities (e.g.,
            from DiffusionAnalyzer.diffusivity). units: cm^2/s
        new_temp (float): desired temperature. units: K

    Returns:
        (float) Diffusivity at extrapolated temp in mS/cm.
    """
    Ea, c, _ = fit_arrhenius(temps, diffusivities)
    return c * np.exp(-Ea / (const.k / const.e * new_temp))


def get_extrapolated_conductivity(temps, diffusivities, new_temp, structure, species):
    """
    Returns extrapolated mS/cm conductivity.

    Args:
        temps ([float]): A sequence of temperatures. units: K
        diffusivities ([float]): A sequence of diffusivities (e.g.,
            from DiffusionAnalyzer.diffusivity). units: cm^2/s
        new_temp (float): desired temperature. units: K
        structure (structure): Structure used for the diffusivity calculation
        species (string/Species): conducting species

    Returns:
        (float) Conductivity at extrapolated temp in mS/cm.
    """
    return get_extrapolated_diffusivity(temps, diffusivities, new_temp) * get_conversion_factor(
        structure, species, new_temp
    )


def get_arrhenius_plot(temps, diffusivities, diffusivity_errors=None, **kwargs):
    r"""
    Returns an Arrhenius plot.

    Args:
        temps ([float]): A sequence of temperatures.
        diffusivities ([float]): A sequence of diffusivities (e.g.,
            from DiffusionAnalyzer.diffusivity).
        diffusivity_errors ([float]): A sequence of errors for the
            diffusivities. If None, no error bar is plotted.
        \\*\\*kwargs:
            Any keyword args supported by matplotlib.pyplot.plot.

    Returns:
        A matplotlib.pyplot object. Do plt.show() to show the plot.
    """
    Ea, c, _ = fit_arrhenius(temps, diffusivities)

    from pymatgen.util.plotting import pretty_plot

    plt = pretty_plot(12, 8)

    # log10 of the arrhenius fit
    arr = c * np.exp(-Ea / (const.k / const.e * np.array(temps)))

    t_1 = 1000 / np.array(temps)

    plt.plot(t_1, diffusivities, "ko", t_1, arr, "k--", markersize=10, **kwargs)
    if diffusivity_errors is not None:
        n = len(diffusivity_errors)
        plt.errorbar(
            t_1[0:n],
            diffusivities[0:n],
            yerr=diffusivity_errors,
            fmt="ko",
            ecolor="k",
            capthick=2,
            linewidth=2,
        )
    ax = plt.axes()
    ax.set_yscale("log")
    plt.text(
        0.6,
        0.85,
        "E$_a$ = {:.0f} meV".format(Ea * 1000),
        fontsize=30,
        transform=plt.axes().transAxes,
    )
    plt.ylabel("D (cm$^2$/s)")
    plt.xlabel("1000/T (K$^{-1}$)")
    plt.tight_layout()
    return plt
