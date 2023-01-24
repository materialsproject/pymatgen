from __future__ import annotations

import os
import unittest

from pytest import approx

try:
    import phonopy
    from phonopy.phonon.dos import TotalDos
except ImportError as ex:
    print(ex)
    phonopy = None
    TotalDos = None

from pymatgen.io.phonopy import get_gruneisen_ph_bs_symm_line, get_gruneisenparameter
from pymatgen.phonon.gruneisen import GruneisenParameter
from pymatgen.phonon.plotter import (
    GruneisenPhononBandStructureSymmLine,
    GruneisenPhononBSPlotter,
    GruneisenPlotter,
)
from pymatgen.util.testing import PymatgenTest


class GruneisenPhononBandStructureSymmLineTest(PymatgenTest):
    def setUp(self) -> None:
        self.bs_symm_line = get_gruneisen_ph_bs_symm_line(
            gruneisen_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_eq_plus_minus_InP.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR_InP"),
            fit=True,
        )

    def test_plot(self):
        plotter = GruneisenPhononBSPlotter(bs=self.bs_symm_line)
        plt = plotter.get_plot_gs()
        assert str(type(plt)) == "<class 'module'>"

    def test_as_dict_from_dict(self):
        new_dict = self.bs_symm_line.as_dict()
        self.new_bs_symm_line = GruneisenPhononBandStructureSymmLine.from_dict(new_dict)
        plotter = GruneisenPhononBSPlotter(bs=self.new_bs_symm_line)
        plt = plotter.get_plot_gs()
        assert str(type(plt)) == "<class 'module'>"


@unittest.skipIf(TotalDos is None, "Phonopy not present")
class GruneisenParameterTest(PymatgenTest):
    def setUp(self) -> None:
        self.gruneisen_obj = get_gruneisenparameter(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_mesh_InP.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR_InP"),
        )
        self.gruneisen_obj_small = get_gruneisenparameter(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_mesh_only_one_q_InP.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR_InP"),
        )
        self.gruneisen_obj_Si = get_gruneisenparameter(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_mesh_Si.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR_Si"),
        )

    def test_plot(self):
        plotter = GruneisenPlotter(self.gruneisen_obj)
        plt = plotter.get_plot(units="mev")
        assert str(type(plt)) == "<class 'module'>"

    def test_fromdict_asdict(self):
        new_dict = self.gruneisen_obj.as_dict()
        self.gruneisen_obj2 = GruneisenParameter.from_dict(new_dict)

    def test_frequencies(self):
        assert self.gruneisen_obj_small.frequencies[0] == approx(0.1264214687)
        assert self.gruneisen_obj_small.frequencies[1] == approx(0.1264214687)
        assert self.gruneisen_obj_small.frequencies[2] == approx(0.2527200484)
        assert self.gruneisen_obj_small.frequencies[3] == approx(8.8520245263)
        assert self.gruneisen_obj_small.frequencies[4] == approx(8.8520245263)
        assert self.gruneisen_obj_small.frequencies[5] == approx(9.6601659578)

    def test_multi(self):
        assert self.gruneisen_obj_small.multiplicities[0] == approx(1)
        assert self.gruneisen_obj.multiplicities[0] == approx(2)

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


if __name__ == "__main__":
    unittest.main()
