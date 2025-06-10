"""This module provides classes to handle the calculation of the IR spectra
This implementation is adapted from Abipy
https://github.com/abinit/abipy
where it was originally done by Guido Petretto and Matteo Giantomassi.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import orjson
from monty.json import MSONable

from pymatgen.core.spectrum import Spectrum
from pymatgen.core.structure import Structure
from pymatgen.util.plotting import add_fig_kwargs
from pymatgen.vis.plotters import SpectrumPlotter

if TYPE_CHECKING:
    from collections.abc import Sequence

    from matplotlib.axes import Axes
    from numpy.typing import ArrayLike
    from typing_extensions import Self

    from pymatgen.util.typing import PathLike

__author__ = "Henrique Miranda, Guido Petretto, Matteo Giantomassi"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Henrique Miranda"
__email__ = "miranda.henrique@gmail.com"
__date__ = "Oct 31, 2018"


class IRDielectricTensor(MSONable):
    """Handle the Ionic Dielectric Tensor
    The implementation is adapted from Abipy
    See the definitions Eq.(53-54) in :cite:`Gonze1997` PRB55, 10355 (1997).
    """

    def __init__(
        self,
        oscillator_strength: ArrayLike,
        ph_freqs_gamma: ArrayLike,
        epsilon_infinity: ArrayLike,
        structure: Structure,
    ) -> None:
        """
        Args:
            oscillator_strength: IR oscillator strengths as defined in Eq. 54 in
                :cite:`Gonze1997` PRB55, 10355 (1997).
            ph_freqs_gamma: Phonon frequencies at the Gamma point
            epsilon_infinity: electronic susceptibility as defined in Eq. 29.
            structure: A Structure object corresponding to the structure used for the calculation.
        """
        self.structure = structure
        self.oscillator_strength = np.array(oscillator_strength).real
        self.ph_freqs_gamma = np.array(ph_freqs_gamma)
        self.epsilon_infinity = np.array(epsilon_infinity)

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Get IRDielectricTensor from dict representation."""
        structure = Structure.from_dict(dct["structure"])
        oscillator_strength = dct["oscillator_strength"]
        ph_freqs_gamma = dct["ph_freqs_gamma"]
        epsilon_infinity = dct["epsilon_infinity"]
        return cls(oscillator_strength, ph_freqs_gamma, epsilon_infinity, structure)

    @property
    def max_phfreq(self) -> float:
        """Maximum phonon frequency."""
        return max(self.ph_freqs_gamma)

    @property
    def nph_freqs(self) -> int:
        """Number of phonon frequencies."""
        return len(self.ph_freqs_gamma)

    def as_dict(self) -> dict:
        """JSON-serializable dict representation of IRDielectricTensor."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "oscillator_strength": self.oscillator_strength.tolist(),
            "ph_freqs_gamma": self.ph_freqs_gamma.tolist(),
            "structure": self.structure.as_dict(),
            "epsilon_infinity": self.epsilon_infinity.tolist(),
        }

    def write_json(self, filename: PathLike) -> None:
        """Save a JSON file with this data."""
        with open(filename, "wb") as file:
            file.write(orjson.dumps(self.as_dict()))

    def get_ir_spectra(
        self,
        broad: list | float = 0.00005,
        emin: float = 0,
        emax: float | None = None,
        divs: int = 500,
    ) -> tuple:
        """The IR spectra is obtained for the different directions.

        Args:
            broad: a list of broadenings or a single broadening for the phonon peaks
            emin (float): minimum energy in which to obtain the spectra. Defaults to 0.
            emax (float): maximum energy in which to obtain the spectra. Defaults to None.
            divs: number of frequency samples between emin and emax

        Returns:
            frequencies: divs array with the frequencies at which the
                         dielectric tensor is calculated
            dielectric_tensor: divsx3x3 numpy array with the dielectric tensor
                         for the range of frequencies
        """
        if isinstance(broad, float):
            broad = [broad] * self.nph_freqs
        if isinstance(broad, list) and len(broad) != self.nph_freqs:
            raise ValueError("The number of elements in the broad_list is not the same as the number of frequencies")

        if emax is None:
            emax = self.max_phfreq + max(broad) * 20
        frequencies = np.linspace(emin, emax, divs)

        dielectric_tensor = np.zeros((divs, 3, 3), dtype=complex)
        for i in range(3, len(self.ph_freqs_gamma)):
            g = broad[i] * self.ph_freqs_gamma[i]
            num = self.oscillator_strength[i, :, :]
            den = self.ph_freqs_gamma[i] ** 2 - frequencies[:, None, None] ** 2 - 1j * g
            dielectric_tensor += num / den
        dielectric_tensor += self.epsilon_infinity[None, :, :]

        return frequencies, dielectric_tensor

    @add_fig_kwargs
    def plot(
        self,
        components: Sequence = ("xx",),
        reim: str = "reim",
        show_phonon_frequencies: bool = True,
        xlim: float | None = None,
        ylim: float | None = None,
        **kwargs,
    ) -> Axes:
        """Helper function to generate the Spectrum plotter and directly plot the results.

        Arguments:
            components: A list with the components of the dielectric tensor to plot.
                Can be either two indexes or a string like 'xx' to plot the (0,0) component
            reim: If 're' (im) is present in the string plots the real (imaginary) part of the dielectric tensor
            show_phonon_frequencies: plot a dot where the phonon frequencies are to help identify IR inactive modes
            xlim: x-limits of the plot. Defaults to None for automatic determination.
            ylim: y-limits of the plot. Defaults to None for automatic determination.
            kwargs: keyword arguments passed to the plotter
        """
        plotter = self.get_plotter(components=components, reim=reim, **kwargs)
        ax = plotter.get_plot(xlim=xlim, ylim=ylim)

        if show_phonon_frequencies:
            ph_freqs_gamma = self.ph_freqs_gamma[3:]
            ax.scatter(ph_freqs_gamma * 1000, np.zeros_like(ph_freqs_gamma))
        ax.set_xlabel(r"$\epsilon(\omega)$")
        ax.set_xlabel(r"Frequency (meV)")
        return ax

    def get_spectrum(
        self,
        component: Sequence | str,
        reim: str,
        broad: list | float = 0.00005,
        emin: float = 0,
        emax: float | None = None,
        divs: int = 500,
        label=None,
    ) -> Spectrum:
        """component: either two indexes or a string like 'xx' to plot the (0,0) component
        reim: only "re" or "im"
        broad: a list of broadenings or a single broadening for the phonon peaks.
        """
        # some check on component and reim value? but not really necessary maybe

        directions_map = {"x": 0, "y": 1, "z": 2, 0: 0, 1: 1, 2: 2}
        functions_map = {"re": lambda x: x.real, "im": lambda x: x.imag}
        reim_label = {"re": "Re", "im": "Im"}
        i, j = (directions_map[direction] for direction in component)
        label = rf"{reim_label[reim]}{{$\epsilon_{{{'xyz'[i]}{'xyz'[j]}}}$}}"

        frequencies, dielectric_tensor = self.get_ir_spectra(broad=broad, emin=emin, emax=emax, divs=divs)
        y = functions_map[reim](dielectric_tensor[:, i, j])

        return Spectrum(frequencies * 1000, y, label=label)

    def get_plotter(
        self,
        components: Sequence = ("xx",),
        reim: str = "reim",
        broad: list | float = 0.00005,
        emin: float = 0,
        emax: float | None = None,
        divs: int = 500,
        **kwargs,
    ) -> SpectrumPlotter:
        """Return an instance of the Spectrum plotter containing the different requested components.

        Arguments:
            components: A list with the components of the dielectric tensor to plot.
                        Can be either two indexes or a string like 'xx' to plot the (0,0) component
            reim: If 're' (im) is present in the string plots the real (imaginary) part of the dielectric tensor
            broad (float): a list of broadenings or a single broadening for the phonon peaks. Defaults to 0.00005.
            emin (float): minimum energy in which to obtain the spectra. Defaults to 0.
            emax (float): maximum energy in which to obtain the spectra. Defaults to None.
            divs: number of frequency samples between emin and emax
            **kwargs: Passed to IRDielectricTensor.get_spectrum()
        """
        directions_map = {"x": 0, "y": 1, "z": 2, 0: 0, 1: 1, 2: 2}
        reim_label = {"re": "Re", "im": "Im"}

        plotter = SpectrumPlotter()
        for component in components:
            i, j = (directions_map[direction] for direction in component)
            for re_im in ("re", "im"):
                if re_im in reim:
                    label = rf"{reim_label[re_im]}{{$\epsilon_{{{'xyz'[i]}{'xyz'[j]}}}$}}"
                    spectrum = self.get_spectrum(
                        component,
                        re_im,
                        broad=broad,
                        emin=emin,
                        emax=emax,
                        divs=divs,
                        **kwargs,
                    )
                    spectrum.XLABEL = r"Frequency (meV)"
                    spectrum.YLABEL = r"$\epsilon(\omega)$"
                    plotter.add_spectrum(label, spectrum)

        return plotter
