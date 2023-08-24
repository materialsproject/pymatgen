from __future__ import annotations

import unittest

import matplotlib.pyplot as plt
from pytest import approx

from pymatgen.io.phonopy import get_gruneisen_ph_bs_symm_line, get_gruneisenparameter
from pymatgen.phonon.gruneisen import GruneisenParameter
from pymatgen.phonon.plotter import GruneisenPhononBandStructureSymmLine, GruneisenPhononBSPlotter, GruneisenPlotter
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

try:
    import phonopy
    from phonopy.phonon.dos import TotalDos
except ImportError as exc:
    print(exc)
    phonopy = TotalDos = None


class TestGruneisenPhononBandStructureSymmLine(PymatgenTest):
    def setUp(self) -> None:
        self.bs_symm_line = get_gruneisen_ph_bs_symm_line(
            gruneisen_path=f"{TEST_FILES_DIR}/gruneisen/gruneisen_eq_plus_minus_InP.yaml",
            structure_path=f"{TEST_FILES_DIR}/gruneisen/eq/POSCAR_InP",
            fit=True,
        )

    def test_plot(self):
        plotter = GruneisenPhononBSPlotter(bs=self.bs_symm_line)
        ax = plotter.get_plot_gs()
        assert isinstance(ax, plt.Axes)

    def test_as_dict_from_dict(self):
        new_dict = self.bs_symm_line.as_dict()
        self.new_bs_symm_line = GruneisenPhononBandStructureSymmLine.from_dict(new_dict)
        plotter = GruneisenPhononBSPlotter(bs=self.new_bs_symm_line)
        ax = plotter.get_plot_gs()
        assert isinstance(ax, plt.Axes)


@unittest.skipIf(TotalDos is None, "Phonopy not present")
class TestGruneisenParameter(PymatgenTest):
    def setUp(self) -> None:
        self.gruneisen_obj = get_gruneisenparameter(
            f"{TEST_FILES_DIR}/gruneisen/gruneisen_mesh_InP.yaml",
            structure_path=f"{TEST_FILES_DIR}/gruneisen/eq/POSCAR_InP",
        )
        self.gruneisen_obj_small = get_gruneisenparameter(
            f"{TEST_FILES_DIR}/gruneisen/gruneisen_mesh_only_one_q_InP.yaml",
            structure_path=f"{TEST_FILES_DIR}/gruneisen/eq/POSCAR_InP",
        )
        self.gruneisen_obj_Si = get_gruneisenparameter(
            f"{TEST_FILES_DIR}/gruneisen/gruneisen_mesh_Si.yaml",
            structure_path=f"{TEST_FILES_DIR}/gruneisen/eq/POSCAR_Si",
        )

    def test_plot(self):
        plotter = GruneisenPlotter(self.gruneisen_obj)
        ax = plotter.get_plot(units="mev")
        assert isinstance(ax, plt.Axes)

    def test_fromdict_asdict(self):
        new_dict = self.gruneisen_obj.as_dict()
        self.gruneisen_obj2 = GruneisenParameter.from_dict(new_dict)

    def test_frequencies(self):
        assert self.gruneisen_obj_small.frequencies == approx(
            [0.12642146, 0.12642146, 0.25272004, 8.85202452, 8.85202452, 9.66016595]
        )

    def test_multi(self):
        assert self.gruneisen_obj_small.multiplicities[0] == 1
        assert self.gruneisen_obj.multiplicities[0] == 2

    def test_gruneisen(self):
        assert self.gruneisen_obj_small.gruneisen[0] == approx(-0.6176464482)
        assert self.gruneisen_obj_small.gruneisen[5] == approx(1.7574050911)

    def test_tdos(self):
        tdos = self.gruneisen_obj.tdos
        assert isinstance(tdos, phonopy.phonon.dos.TotalDos)

    def test_phdos(self):
        assert self.gruneisen_obj.phdos.cv(298.15) == approx(45.17772584681599)

    def test_average_gruneisen(self):
        assert self.gruneisen_obj.average_gruneisen() == approx(1.164231026696211)
        assert self.gruneisen_obj.average_gruneisen(squared=False) == approx(0.849759667411049)
        assert self.gruneisen_obj.average_gruneisen(limit_frequencies="debye") == approx(0.848865124114612)
        assert self.gruneisen_obj.average_gruneisen(limit_frequencies="acoustic") == approx(1.283180896570312)
        assert self.gruneisen_obj_Si.average_gruneisen() == approx(1.1090815951892143)

    def test_thermal_conductivity_slack(self):
        assert self.gruneisen_obj.thermal_conductivity_slack() == approx(77.97582174520458)
        assert self.gruneisen_obj.thermal_conductivity_slack(t=300) == approx(88.94562145031158)
        assert self.gruneisen_obj_Si.thermal_conductivity_slack(t=300) == approx(127.69008331982265)

    def test_debye_temp_phonopy(self):
        # This is the correct conversion when starting from THz in the debye_freq
        assert self.gruneisen_obj_small.debye_temp_phonopy() == approx(473.31932718764284)

    def test_acoustic_debye_temp(self):
        assert self.gruneisen_obj_small.acoustic_debye_temp == approx(317.54811309631845)
        assert self.gruneisen_obj.acoustic_debye_temp == approx(342.2046198151735)
        assert self.gruneisen_obj_Si.acoustic_debye_temp == approx(526.0725636300882)
