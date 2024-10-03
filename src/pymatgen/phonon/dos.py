"""This module defines classes to represent the phonon density of states."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, NamedTuple

import numpy as np
import scipy.constants as const
from monty.functools import lazy_property
from monty.json import MSONable
from packaging import version
from scipy.ndimage import gaussian_filter1d
from scipy.stats import wasserstein_distance

from pymatgen.core.structure import Structure
from pymatgen.util.coord import get_linear_interpolated_value

if version.parse(np.__version__) < version.parse("2.0.0"):
    np.trapezoid = np.trapz  # noqa: NPY201

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray
    from typing_extensions import Self

BOLTZ_THZ_PER_K = const.value("Boltzmann constant in Hz/K") / const.tera  # Boltzmann constant in THz/K
THZ_TO_J = const.value("hertz-joule relationship") * const.tera


class PhononDos(MSONable):
    """Basic DOS object. All other DOS objects are extended versions of this object."""

    def __init__(self, frequencies: Sequence, densities: Sequence) -> None:
        """
        Args:
            frequencies: A sequence of frequencies in THz
            densities: A sequence representing the density of states.
        """
        self.frequencies = np.array(frequencies)
        self.densities = np.array(densities)

    def get_smeared_densities(self, sigma: float) -> np.ndarray:
        """Get the densities, but with a Gaussian smearing of
        std dev sigma applied.

        Args:
            sigma: Std dev of Gaussian smearing function. In units of
                THz. Common values are 0.01 - 0.1 THz.

        Returns:
            np.array: Gaussian-smeared DOS densities.
        """
        if sigma == 0:
            return self.densities
        diff = [self.frequencies[idx + 1] - self.frequencies[idx] for idx in range(len(self.frequencies) - 1)]
        avg_diff = sum(diff) / len(diff)

        return gaussian_filter1d(self.densities, sigma / avg_diff)

    def __add__(self, other: PhononDos) -> PhononDos:
        """Add two DOS together. Pads densities with zeros to make frequencies matching.

        Args:
            other: Another DOS object.

        Returns:
            Sum of the two DOSs.
        """
        if isinstance(other, int | float):
            return PhononDos(self.frequencies, self.densities + other)
        if not all(np.equal(self.frequencies, other.frequencies)):
            raise ValueError("Frequencies of both DOS are not compatible!")
        densities = self.densities + other.densities
        return PhononDos(self.frequencies, densities)

    def __sub__(self, other: PhononDos) -> PhononDos:
        """Subtracts two DOS together. Pads densities with zeros to make frequencies matching.

        Args:
            other: Another DOS object.

        Returns:
            Difference of the two DOSs.
        """
        return self + (-other)

    def __mul__(self, scalar: float) -> PhononDos:
        """Multiplies the DOS by a scalar.

        Args:
            scalar: A scalar to multiply by.

        Returns:
            A new DOS multiplied by a scalar.
        """
        return PhononDos(self.frequencies, self.densities * scalar)

    def __neg__(self) -> PhononDos:
        """Inverts the DOS.

        Returns:
            A new DOS with densities inverted. Useful for subtracting from a total DOS.
        """
        return PhononDos(self.frequencies, -self.densities)

    __radd__ = __add__
    __rmul__ = __mul__

    def __eq__(self, other: object) -> bool:
        """Two DOS are equal if their densities are equal.

        Args:
            other: Another DOS object.

        Returns:
            True if densities are equal.
        """
        if not isinstance(other, PhononDos):
            return NotImplemented
        return np.allclose(self.densities, other.densities)

    def __repr__(self) -> str:
        frequencies, densities = self.frequencies.shape, self.densities.shape
        n_positive_freqs = len(self._positive_frequencies)
        return f"{type(self).__name__}({frequencies=}, {densities=}, {n_positive_freqs=})"

    def get_interpolated_value(self, frequency) -> float:
        """Get interpolated density for a particular frequency.

        Args:
            frequency: frequency to return the density for.
        """
        return get_linear_interpolated_value(self.frequencies, self.densities, frequency)

    def __str__(self) -> str:
        """Get a string which can be easily plotted (using gnuplot)."""
        str_arr = [f"#{'Frequency':30s} {'Density':30s}"]
        for idx, freq in enumerate(self.frequencies):
            str_arr.append(f"{freq:.5f} {self.densities[idx]:.5f}")
        return "\n".join(str_arr)

    @classmethod
    def from_dict(cls, dct: dict[str, Sequence]) -> Self:
        """Get PhononDos object from dict representation of PhononDos."""
        return cls(dct["frequencies"], dct["densities"])

    def as_dict(self) -> dict:
        """JSON-serializable dict representation of PhononDos."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "frequencies": list(self.frequencies),
            "densities": list(self.densities),
        }

    @lazy_property
    def ind_zero_freq(self) -> int:
        """Index of the first point for which the frequencies are >= 0."""
        ind = np.searchsorted(self.frequencies, 0)
        if ind >= len(self.frequencies):
            raise ValueError("No positive frequencies found")
        return ind

    @lazy_property
    def _positive_frequencies(self) -> np.ndarray:
        """Numpy array containing the list of positive frequencies."""
        return self.frequencies[self.ind_zero_freq :]

    @lazy_property
    def _positive_densities(self) -> np.ndarray:
        """Numpy array containing the list of densities corresponding to positive frequencies."""
        return self.densities[self.ind_zero_freq :]

    def cv(self, temp: float | None = None, structure: Structure | None = None, **kwargs) -> float:
        """Constant volume specific heat C_v at temperature T obtained from the integration of the DOS.
        Only positive frequencies will be used.
        Result in J/(K*mol-c). A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/(K*mol).

        Args:
            temp: a temperature in K
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units
            **kwargs: allows passing in deprecated t parameter for temp

        Returns:
            float: Constant volume specific heat C_v
        """
        temp = kwargs.get("t", temp)
        if temp == 0:
            return 0

        freqs = self._positive_frequencies
        dens = self._positive_densities

        def csch2(x):
            return 1.0 / (np.sinh(x) ** 2)

        wd2kt = freqs / (2 * BOLTZ_THZ_PER_K * temp)
        cv = np.trapezoid(wd2kt**2 * csch2(wd2kt) * dens, x=freqs)
        cv *= const.Boltzmann * const.Avogadro

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            cv /= formula_units

        return cv

    def entropy(self, temp: float | None = None, structure: Structure | None = None, **kwargs) -> float:
        """Vibrational entropy at temperature T obtained from the integration of the DOS.
        Only positive frequencies will be used.
        Result in J/(K*mol-c). A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/(K*mol).

        Args:
            temp: a temperature in K
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units
            **kwargs: allows passing in deprecated t parameter for temp

        Returns:
            float: Vibrational entropy
        """
        temp = kwargs.get("t", temp)
        if temp == 0:
            return 0

        freqs = self._positive_frequencies
        dens = self._positive_densities

        wd2kt = freqs / (2 * BOLTZ_THZ_PER_K * temp)
        entropy = np.trapezoid((wd2kt * 1 / np.tanh(wd2kt) - np.log(2 * np.sinh(wd2kt))) * dens, x=freqs)

        entropy *= const.Boltzmann * const.Avogadro

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            entropy /= formula_units

        return entropy

    def internal_energy(self, temp: float | None = None, structure: Structure | None = None, **kwargs) -> float:
        """Phonon contribution to the internal energy at temperature T obtained from the integration of the DOS.
        Only positive frequencies will be used.
        Result in J/mol-c. A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/mol.

        Args:
            temp: a temperature in K
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units
            **kwargs: allows passing in deprecated t parameter for temp

        Returns:
            float: Phonon contribution to the internal energy
        """
        temp = kwargs.get("t", temp)
        if temp == 0:
            return self.zero_point_energy(structure=structure)

        freqs = self._positive_frequencies
        dens = self._positive_densities

        wd2kt = freqs / (2 * BOLTZ_THZ_PER_K * temp)
        e_phonon = np.trapezoid(freqs * 1 / np.tanh(wd2kt) * dens, x=freqs) / 2

        e_phonon *= THZ_TO_J * const.Avogadro

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            e_phonon /= formula_units

        return e_phonon

    def helmholtz_free_energy(self, temp: float | None = None, structure: Structure | None = None, **kwargs) -> float:
        """Phonon contribution to the Helmholtz free energy at temperature T obtained from the integration of the DOS.
        Only positive frequencies will be used.
        Result in J/mol-c. A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/mol.

        Args:
            temp: a temperature in K
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units
            **kwargs: allows passing in deprecated t parameter for temp

        Returns:
            float: Phonon contribution to the Helmholtz free energy
        """
        temp = kwargs.get("t", temp)
        if temp == 0:
            return self.zero_point_energy(structure=structure)

        freqs = self._positive_frequencies
        dens = self._positive_densities

        wd2kt = freqs / (2 * BOLTZ_THZ_PER_K * temp)
        e_free = np.trapezoid(np.log(2 * np.sinh(wd2kt)) * dens, x=freqs)

        e_free *= const.Boltzmann * const.Avogadro * temp

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            e_free /= formula_units

        return e_free

    def zero_point_energy(self, structure: Structure | None = None) -> float:
        """Zero point energy of the system. Only positive frequencies will be used.
        Result in J/mol-c. A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/mol.

        Args:
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units

        Returns:
            Phonon contribution to the internal energy
        """
        freqs = self._positive_frequencies
        dens = self._positive_densities

        zpe = 0.5 * np.trapezoid(freqs * dens, x=freqs)
        zpe *= THZ_TO_J * const.Avogadro

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            zpe /= formula_units

        return zpe

    def mae(self, other: PhononDos, two_sided: bool = True) -> float:
        """Mean absolute error between two DOSs.

        Args:
            other (PhononDos): Another phonon DOS
            two_sided (bool): Whether to calculate the two-sided MAE meaning interpolate each DOS to the
                other's frequencies and averaging the two MAEs. Defaults to True.

        Returns:
            float: Mean absolute error.
        """
        # Interpolate other.densities to align with self.frequencies
        self_interpolated = np.interp(self.frequencies, other.frequencies, other.densities)
        self_mae = np.abs(self.densities - self_interpolated).mean()

        if two_sided:
            other_interpolated = np.interp(other.frequencies, self.frequencies, self.densities)
            other_mae = np.abs(other.densities - other_interpolated).mean()
            return (self_mae + other_mae) / 2

        return self_mae

    def r2_score(self, other: PhononDos) -> float:
        """R^2 score between two DOSs.

        Args:
            other (PhononDos): Another phonon DOS

        Returns:
            float: R^2 score
        """
        var = self.densities.var()
        if var == 0:
            return 0
        mse = ((self.densities - other.densities) ** 2).mean()
        return 1 - mse / var

    def get_last_peak(self, threshold: float = 0.05) -> float:
        """Find the last peak in the phonon DOS defined as the highest frequency with a DOS
        value at least threshold * height of the overall highest DOS peak.
        A peak is any local maximum of the DOS as a function of frequency.
        Use dos.get_interpolated_value(peak_freq) to get density at peak_freq.

        TODO method added by @janosh on 2023-12-18. seems to work in most cases but
        was not extensively tested. PRs with improvements welcome!

        Args:
            threshold (float, optional): Minimum ratio of the height of the last peak
                to the height of the highest peak. Defaults to 0.05 = 5%. In case no peaks
                are high enough to match, the threshold is reset to half the height of the
                second-highest peak.

        Returns:
            float: last DOS peak frequency (in THz)
        """
        first_deriv = np.gradient(self.densities, self.frequencies)
        second_deriv = np.gradient(first_deriv, self.frequencies)

        maxima = (  # maxima indices of the first DOS derivative w.r.t. frequency
            (first_deriv[:-1] > 0) & (first_deriv[1:] < 0) & (second_deriv[:-1] < 0)
        )
        # get mean of the two nearest frequencies around the maximum as better estimate
        maxima_freqs = (self.frequencies[:-1][maxima] + self.frequencies[1:][maxima]) / 2

        # filter maxima based on the threshold
        max_dos = max(self.densities)
        threshold *= max_dos
        filtered_maxima_freqs = maxima_freqs[self.densities[:-1][maxima] >= threshold]

        if len(filtered_maxima_freqs) == 0:
            # if no maxima reach the threshold (i.e. 1 super high peak and all other peaks
            # tiny), use half the height of second highest peak as threshold
            second_highest_peak = sorted(self.densities)[-2]
            threshold = second_highest_peak / 2
            filtered_maxima_freqs = maxima_freqs[self.densities[:-1][maxima] >= threshold]

        return max(filtered_maxima_freqs)

    def get_dos_fp(
        self,
        binning: bool = True,
        min_f: float | None = None,
        max_f: float | None = None,
        n_bins: int = 256,
        normalize: bool = True,
    ) -> PhononDosFingerprint:
        """Generate the DOS fingerprint.

        Args:
            binning (bool): If true, the DOS fingerprint is binned using np.linspace and n_bins.
                Default is True.
            min_f (float): The minimum mode frequency to include in the fingerprint (default is None)
            max_f (float): The maximum mode frequency to include in the fingerprint (default is None)
            n_bins (int): Number of bins to be used in the fingerprint (default is 256)
            normalize (bool): If true, normalizes the area under fp to equal to 1. Default is True.

        Returns:
            PhononDosFingerprint: The phonon density of states fingerprint
                of format (frequencies, densities, n_bins, bin_width)
        """
        frequencies = self.frequencies

        if max_f is None:
            max_f = np.max(frequencies)

        if min_f is None:
            min_f = np.min(frequencies)

        densities = self.densities

        if len(frequencies) < n_bins:
            inds = np.where((frequencies >= min_f) & (frequencies <= max_f))
            return PhononDosFingerprint(frequencies[inds], densities[inds], len(frequencies), np.diff(frequencies)[0])

        if binning:
            freq_bounds = np.linspace(min_f, max_f, n_bins + 1)
            freq = freq_bounds[:-1] + (freq_bounds[1] - freq_bounds[0]) / 2.0
            bin_width = np.diff(freq)[0]
        else:
            freq_bounds = np.array(frequencies)
            freq = np.append(frequencies, [frequencies[-1] + np.abs(frequencies[-1]) / 10])
            n_bins = len(frequencies)
            bin_width = np.diff(frequencies)[0]

        dos_rebin = np.zeros(freq.shape)

        for ii, e1, e2 in zip(range(len(freq)), freq_bounds[:-1], freq_bounds[1:], strict=False):
            inds = np.where((frequencies >= e1) & (frequencies < e2))
            dos_rebin[ii] = np.sum(densities[inds])
        if normalize:  # scale DOS bins to make area under histogram equal 1
            area = np.sum(dos_rebin * bin_width)
            dos_rebin_sc = dos_rebin / area
        else:
            dos_rebin_sc = dos_rebin

        return PhononDosFingerprint(np.array([freq]), dos_rebin_sc, n_bins, bin_width)

    @staticmethod
    def fp_to_dict(fp: PhononDosFingerprint) -> dict:
        """Convert a fingerprint into a dictionary.

        Args:
            fp: The DOS fingerprint to be converted into a dictionary

        Returns:
            dict: A dict of the fingerprint Keys=type, Values=np.ndarray(frequencies, densities)
        """
        fp_dict = {}
        fp_dict[fp[2]] = np.array([fp[0], fp[1]], dtype="object").T

        return fp_dict

    @staticmethod
    def get_dos_fp_similarity(
        fp1: PhononDosFingerprint,
        fp2: PhononDosFingerprint,
        col: int = 1,
        pt: int | str = "All",
        normalize: bool = False,
        metric: Literal["tanimoto", "wasserstein", "cosine-sim"] = "tanimoto",
    ) -> float:
        """Calculate the similarity index between two fingerprints.

        Args:
            fp1 (PhononDosFingerprint): The 1st dos fingerprint object
            fp2 (PhononDosFingerprint): The 2nd dos fingerprint object
            col (int): The item in the fingerprints (0:frequencies,1: densities) to compute
                the similarity index of  (default is 1)
            pt (int or str) : The index of the point that the dot product is to be taken (default is All)
            normalize (bool): If True normalize the scalar product to 1 (default is False)
            metric (Literal): Metric used to compute similarity default is "tanimoto".

        Raises:
            ValueError: If metric other than "tanimoto", "wasserstein" and "cosine-sim" is requested.
            ValueError: If normalize is set to True along with the metric.

        Returns:
            float: Similarity index given by the dot product
        """
        valid_metrics = ("tanimoto", "wasserstein", "cosine-sim")
        if metric not in valid_metrics:
            raise ValueError(f"Invalid {metric=}, choose from {valid_metrics}.")

        fp1_dict = CompletePhononDos.fp_to_dict(fp1) if not isinstance(fp1, dict) else fp1

        fp2_dict = CompletePhononDos.fp_to_dict(fp2) if not isinstance(fp2, dict) else fp2

        if pt == "All":
            vec1 = np.array([pt[col] for pt in fp1_dict.values()]).flatten()
            vec2 = np.array([pt[col] for pt in fp2_dict.values()]).flatten()
        else:
            vec1 = fp1_dict[fp1[2][pt]][col]  # type: ignore[index]
            vec2 = fp2_dict[fp2[2][pt]][col]  # type: ignore[index]

        if not normalize and metric == "tanimoto":
            rescale = np.linalg.norm(vec1) ** 2 + np.linalg.norm(vec2) ** 2 - np.dot(vec1, vec2)
            return np.dot(vec1, vec2) / rescale

        if not normalize and metric == "wasserstein":
            return wasserstein_distance(
                u_values=np.cumsum(vec1 * fp1.bin_width), v_values=np.cumsum(vec2 * fp2.bin_width)
            )

        if normalize and metric == "cosine-sim":
            rescale = np.linalg.norm(vec1) * np.linalg.norm(vec2)
            return np.dot(vec1, vec2) / rescale

        if not normalize and metric == "cosine-sim":
            rescale = 1.0
            return np.dot(vec1, vec2) / rescale

        raise ValueError("Cannot compute similarity index. When normalize=True, then please set metric=cosine-sim")


class PhononDosFingerprint(NamedTuple):
    """
    Represents a Phonon Density of States (DOS) fingerprint.

    This named tuple is used to store information related to the Density of States (DOS)
    in a material. It includes the frequencies, densities, number of bins, and bin width.

    Args:
        frequencies: The frequency values associated with the DOS.
        densities: The corresponding density values for each energy.
        n_bins: The number of bins used in the fingerprint.
        bin_width: The width of each bin in the DOS fingerprint.
    """

    frequencies: NDArray
    densities: NDArray
    n_bins: int
    bin_width: float


class CompletePhononDos(PhononDos):
    """This wrapper class defines a total dos, and also provides a list of PDos.

    Attributes:
        pdos (dict): Dict of partial densities of the form {Site:Densities}.
            Densities are a dict of {Orbital:Values} where Values are a list of floats.
            Site is a pymatgen.core.sites.Site object.
    """

    def __init__(self, structure: Structure, total_dos, ph_doses: dict) -> None:
        """
        Args:
            structure: Structure associated with this particular DOS.
            total_dos: total Dos for structure
            ph_doses: The phonon DOSes are supplied as a dict of {Site: Densities}.
        """
        super().__init__(frequencies=total_dos.frequencies, densities=total_dos.densities)
        self.pdos = {site: np.array(dens) for site, dens in ph_doses.items()}
        self.structure = structure

    def get_site_dos(self, site) -> PhononDos:
        """Get the Dos for a site.

        Args:
            site: Site in Structure associated with CompletePhononDos.

        Returns:
            PhononDos: containing summed orbital densities for site.
        """
        return PhononDos(self.frequencies, self.pdos[site])

    def get_element_dos(self) -> dict:
        """Get element projected Dos.

        Returns:
            dict of {Element: Dos}
        """
        el_dos = {}
        for site, atom_dos in self.pdos.items():
            el = site.specie
            if el not in el_dos:
                el_dos[el] = np.array(atom_dos)
            else:
                el_dos[el] += np.array(atom_dos)
        return {el: PhononDos(self.frequencies, densities) for el, densities in el_dos.items()}

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Get CompleteDos object from dict representation."""
        total_dos = PhononDos.from_dict(dct)
        struct = Structure.from_dict(dct["structure"])
        ph_doses = dict(zip(struct, dct["pdos"], strict=True))

        return cls(struct, total_dos, ph_doses)

    def as_dict(self):
        """JSON-serializable dict representation of CompletePhononDos."""
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "structure": self.structure.as_dict(),
            "frequencies": list(self.frequencies),
            "densities": list(self.densities),
            "pdos": [],
        }
        if len(self.pdos) > 0:
            for site in self.structure:
                dct["pdos"].append(list(self.pdos[site]))
        return dct

    def __str__(self) -> str:
        return f"Complete phonon DOS for {self.structure}"
