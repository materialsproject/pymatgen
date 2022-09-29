"""
This is a module for XPS analysis. It is modelled after the Galore package (https://github.com/SMTG-UCL/galore), but
with some modifications for easier analysis from pymatgen itself. Please cite the following original work if you use
this::

    Adam J. Jackson, Alex M. Ganose, Anna Regoutz, Russell G. Egdell, David O. Scanlon (2018). Galore: Broadening and
    weighting for simulation of photoelectron spectroscopy. Journal of Open Source Software, 3(26), 773,
    doi: 10.21105/joss.007733

You may wish to look at the optional dependency galore for more functionality such as plotting and other cross-sections.
Note that the atomic_subshell_photoionization_cross_sections.csv has been reparsed from the original compilation::

    Yeh, J. J.; Lindau, I. Atomic Subshell Photoionization Cross Sections and Asymmetry Parameters: 1 ⩽ Z ⩽ 103.
    Atomic Data and Nuclear Data Tables 1985, 32 (1), 1-155. https://doi.org/10.1016/0092-640X(85)90016-6.

This version contains all detailed information for all orbitals.
"""

import collections
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from pymatgen.core.periodic_table import Element
from pymatgen.core.spectrum import Spectrum
from pymatgen.electronic_structure.dos import CompleteDos


def _load_cross_sections(fname):
    data = pd.read_csv(fname)

    d = collections.defaultdict(dict)
    for row in data.itertuples():
        sym = row.element
        el = Element(sym)
        if el.Z > 92:
            continue
        orb = row.orbital
        shell = int(orb[0])
        orbtype = orb[1]
        nelect = None
        for l in el.full_electronic_structure:
            if l[0] == shell and l[1] == orbtype:
                nelect = l[2]
                break
        if nelect is not None:
            d[sym][orbtype] = row.weight / nelect
    return d


CROSS_SECTIONS = _load_cross_sections(Path(__file__).parent / "atomic_subshell_photoionization_cross_sections.csv")


class XPS(Spectrum):
    """
    Class representing an X-ray photoelectron spectra.
    """

    XLABEL = "Binding Energy (eV)"
    YLABEL = "Intensity"

    @classmethod
    def from_dos(cls, dos: CompleteDos):
        """
        :param dos: CompleteDos object with project element-orbital DOS. Can be obtained from Vasprun.get_complete_dos.
        :param sigma: Smearing for Gaussian.
        :return: XPS
        """
        total = np.zeros(dos.energies.shape)
        for el in dos.structure.composition:
            spd_dos = dos.get_element_spd_dos(el)
            for orb, pdos in spd_dos.items():
                weight = CROSS_SECTIONS[el.symbol].get(str(orb), None)
                if weight is not None:
                    total += pdos.get_densities() * weight
                else:
                    warnings.warn(f"No cross-section for {el}{orb}")
        return XPS(-dos.energies, total / np.max(total))
