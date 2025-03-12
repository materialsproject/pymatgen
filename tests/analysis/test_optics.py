from __future__ import annotations

from pymatgen.analysis.optics import DielectricAnalysis
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/io/vasp/outputs"


def test_dielectric_analysis():
    analysis = DielectricAnalysis.from_vasprun(Vasprun(f"{TEST_DIR}/vasprun.dielectric_6.0.8.xml.gz"))
    print(analysis.data["n_1"])
    from pymatgen.util.plotting import pretty_plot

    ax = pretty_plot(12, 8)
    ax.plot(analysis.data["wavelength"], analysis.data["n_1"], label="n_1")
    import matplotlib.pyplot as plt

    plt.show()
