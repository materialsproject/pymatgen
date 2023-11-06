"""This module defines an excitation spectrum class."""

from __future__ import annotations

from pymatgen.core.spectrum import Spectrum


class ExcitationSpectrum(Spectrum):
    """
    Basic excitation spectrum object.

    Attributes:
        x (Sequence[float]): The sequence of energies.
        y (Sequence[float]): The sequence of mu(E).
    """

    XLABEL = "Energy (eV)"
    YLABEL = "Intensity"

    def __init__(self, x, y):
        """
        Args:
            x: A sequence of x-ray energies in eV
            y: A sequence of intensity values.
        """
        super().__init__(x, y)
