import json
import os
import unittest
from io import open

from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos
from pymatgen.phonon.plotter import PhononBSPlotter, PhononDosPlotter, ThermoPlotter
from pymatgen.util.testing import PymatgenTest


class PhononDosPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "NaCl_complete_ph_dos.json"), "r") as f:
            self.dos = CompletePhononDos.from_dict(json.load(f))
            self.plotter = PhononDosPlotter(sigma=0.2, stack=True)
            self.plotter_nostack = PhononDosPlotter(sigma=0.2, stack=False)

    def test_add_dos_dict(self):
        d = self.plotter.get_dos_dict()
        self.assertEqual(len(d), 0)
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        self.assertEqual(len(d), 2)

    def test_get_dos_dict(self):
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        for el in ["Na", "Cl"]:
            self.assertIn(el, d)

    def test_plot(self):
        # Disabling latex for testing.
        from matplotlib import rc

        rc("text", usetex=False)
        self.plotter.add_dos("Total", self.dos)
        self.plotter.get_plot(units="mev")
        self.plotter_nostack.add_dos("Total", self.dos)
        self.plotter_nostack.get_plot(units="mev")


class PhononBSPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "NaCl_phonon_bandstructure.json"), "r") as f:
            d = json.loads(f.read())
            self.bs = PhononBandStructureSymmLine.from_dict(d)
            self.plotter = PhononBSPlotter(self.bs)

    def test_bs_plot_data(self):
        self.assertEqual(
            len(self.plotter.bs_plot_data()["distances"][0]),
            51,
            "wrong number of distances in the first branch",
        )
        self.assertEqual(len(self.plotter.bs_plot_data()["distances"]), 4, "wrong number of branches")
        self.assertEqual(
            sum([len(e) for e in self.plotter.bs_plot_data()["distances"]]),
            204,
            "wrong number of distances",
        )
        self.assertEqual(self.plotter.bs_plot_data()["ticks"]["label"][4], "Y", "wrong tick label")
        self.assertEqual(
            len(self.plotter.bs_plot_data()["ticks"]["label"]),
            8,
            "wrong number of tick labels",
        )

    def test_plot(self):
        # Disabling latex for testing.
        from matplotlib import rc

        rc("text", usetex=False)
        self.plotter.get_plot(units="mev")

    def test_plot_compare(self):
        # Disabling latex for testing.
        from matplotlib import rc

        rc("text", usetex=False)
        self.plotter.plot_compare(self.plotter, units="mev")


class ThermoPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "NaCl_complete_ph_dos.json"), "r") as f:
            self.dos = CompletePhononDos.from_dict(json.load(f))
            self.plotter = ThermoPlotter(self.dos, self.dos.structure)

    def test_plot_functions(self):
        # Disabling latex for testing.
        from matplotlib import rc

        rc("text", usetex=False)
        self.plotter.plot_cv(5, 100, 5, show=False)
        self.plotter.plot_entropy(5, 100, 5, show=False)
        self.plotter.plot_internal_energy(5, 100, 5, show=False)
        self.plotter.plot_helmholtz_free_energy(5, 100, 5, show=False)
        self.plotter.plot_thermodynamic_properties(5, 100, 5, show=False, fig_close=True)


if __name__ == "__main__":
    unittest.main()
