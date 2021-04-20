"""
This is a module for XPS analysis. It is modelled after the Galore package (https://github.com/SMTG-UCL/galore), but
with some modifications for easier analysis from pymatgen itself. Please cite the following original work if you use
this::

    Adam J. Jackson, Alex M. Ganose, Anna Regoutz, Russell G. Egdell, David O. Scanlon (2018). Galore: Broadening and
    weighting for simulation of photoelectron spectroscopy. Journal of Open Source Software, 3(26), 773,
    doi: 10.21105/joss.007733
"""

from pathlib import Path

import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

from monty.serialization import loadfn

from pymatgen.core.spectrum import Spectrum
from pymatgen.electronic_structure.dos import CompleteDos

CROSS_SECTIONS = loadfn(Path(__file__).parent / "cross_sections_AlKalpha.json")


class XPS(Spectrum):
    """
    Class representing an X-ray photoelectron spectra.
    """

    XLABEL = "Binding Energy (eV)"
    YLABEL = "Intensity"

    @classmethod
    def from_dos(cls, dos: CompleteDos, sigma: float = 0):
        """
        :param dos: CompleteDos object with project element-orbital DOS. Can be obtained from Vasprun.get_complete_dos.
        :param sigma: Smearing for Gaussian.
        :return:
        """
        total = np.zeros(dos.energies.shape)
        for el in dos.structure.composition.keys():
            spd_dos = dos.get_element_spd_dos(el)
            print(el)
            for orb, pdos in spd_dos.items():
                weight = CROSS_SECTIONS[el.symbol].get(str(orb), None)
                if weight is not None:
                    total += pdos.get_densities() * weight

        if sigma != 0:
            diff = [dos.energies[i + 1] - dos.energies[i] for i in range(len(dos.energies) - 1)]
            avgdiff = sum(diff) / len(diff)
            total = gaussian_filter1d(total, sigma / avgdiff)
        return XPS(-dos.energies, total / np.max(total))


if __name__ == "__main__":
    from pymatgen.ext.matproj import MPRester

    mpr = MPRester()
    dos_SnO2 = mpr.get_dos_by_material_id("mp-856")
    xps = XPS.from_dos(dos_SnO2, 0.3)
    import matplotlib.pyplot as plt

    plt.plot(xps.x, xps.y, "k-")
    plt.xlim((-5, 10))
    plt.ylim((0, 0.04))
    plt.show()
