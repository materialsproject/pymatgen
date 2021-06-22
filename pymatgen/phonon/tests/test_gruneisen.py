import copy
import os
import unittest

import phonopy

from pymatgen.core import Structure
from pymatgen.io.phonopy import get_gruneisen_ph_bs_symm_line
from pymatgen.io.phonopy import get_gruneisenparamter
from pymatgen.phonon.gruneisen import GruneisenParameter
from pymatgen.phonon.plotter import GruneisenPhononBSPlotter, GruneisenPhononBandStructureSymmLine, GruneisenPlotter
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.util.testing import PymatgenTest


class GruneisenPhononBandStructureSymmLineTest(PymatgenTest):
    def setUp(self) -> None:
        kpath = HighSymmKpath(Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR")),
                              symprec=0.01)
        structure = kpath.prim
        kpath_dict, kpath_concrete = self.get_kpath(structure)

        self.bs_symm_line = get_gruneisen_ph_bs_symm_line(
            gruneisen_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_eq_plus_minus.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR"),
            labels_dict=kpath_dict["kpoints"],
            fit=True)

    @staticmethod
    def get_kpath(structure: Structure):
        """
        get high-symmetry points in k-space
        Args:
            structure: Structure Object
        Returns:
        """
        kpath = HighSymmKpath(structure, symprec=0.01)
        kpath_save = kpath.kpath
        labels = copy.deepcopy(kpath_save["path"])
        path = copy.deepcopy(kpath_save["path"])

        for ilabelset, labelset in enumerate(labels):
            for ilabel, label in enumerate(labelset):
                path[ilabelset][ilabel] = kpath_save["kpoints"][label]
        return kpath_save, path

    def test_plot(self):
        plotter = GruneisenPhononBSPlotter(bs=self.bs_symm_line)
        plotter.get_plot_gs().show()

    def test_as_dict_from_dict(self):
        new_dict = self.bs_symm_line.as_dict()
        self.new_bs_symm_line = GruneisenPhononBandStructureSymmLine.from_dict(new_dict)
        plotter = GruneisenPhononBSPlotter(bs=self.new_bs_symm_line)
        # TODO: exclude show somehow?
        plotter.get_plot_gs().show()

    def test_fitting_procedure(self):
        pass


class GruneisenParameterTest(PymatgenTest):
    def setUp(self) -> None:
        self.gruneisenobject = get_gruneisenparamter(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_mesh.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR"))
        self.gruneisenobject_small = get_gruneisenparamter(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_mesh_only_one_q.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR"))

    def test_plot(self):
        plotter = GruneisenPlotter(self.gruneisenobject)
        plt = plotter.get_plot(units="mev")
        self.assertEqual(str(type(plt)), "<class 'module'>")

    def test_fromdict_asdict(self):
        new_dict = self.gruneisenobject.as_dict()
        self.gruneisenobject2 = GruneisenParameter.from_dict(new_dict)

    def test_frequencies(self):
        self.assertAlmostEqual(self.gruneisenobject_small.frequencies[0], 0.1264214687)
        self.assertAlmostEqual(self.gruneisenobject_small.frequencies[1], 0.1264214687)
        self.assertAlmostEqual(self.gruneisenobject_small.frequencies[2], 0.2527200484)
        self.assertAlmostEqual(self.gruneisenobject_small.frequencies[3], 8.8520245263)
        self.assertAlmostEqual(self.gruneisenobject_small.frequencies[4], 8.8520245263)
        self.assertAlmostEqual(self.gruneisenobject_small.frequencies[5], 9.6601659578)

    def test_multi(self):
        self.assertAlmostEqual(self.gruneisenobject_small.multiplicities[0], 1)
        self.assertAlmostEqual(self.gruneisenobject.multiplicities[0], 2)

    def test_gruneisen(self):
        self.assertAlmostEqual(self.gruneisenobject_small.gruneisenparamters[0], -0.6176464482)
        self.assertAlmostEqual(self.gruneisenobject_small.gruneisenparamters[5], 1.7574050911)

    def test_tdos(self):
        tdos = self.gruneisenobject.tdos
        self.assertEqual(type(tdos), phonopy.phonon.dos.TotalDos)

    def test_phdos(self):
        self.assertAlmostEqual(self.gruneisenobject.phdos.cv(298.15), 45.17772584681599)

    def test_average_gruneisen(self):
        self.assertAlmostEqual(self.gruneisenobject.average_gruneisen(), 1.164231026696211)
        self.assertAlmostEqual(self.gruneisenobject.average_gruneisen(squared=False), 0.8497596674110489)

    def test_thermal_conductivity_slack(self):
        self.assertAlmostEqual(self.gruneisenobject.thermal_conductivity_slack(), 77.97582174520458)
        self.assertAlmostEqual(self.gruneisenobject.thermal_conductivity_slack(t=300), 88.94562145031158)

    def test_debye_temp_phonopy(self):
        # This is the correct conversion when starting from THz in the debye_freq
        self.assertAlmostEqual(self.gruneisenobject_small.debye_temp_phonopy, 473.31932718764284)

    def test_acoustic_debye_temp(self):
        self.assertAlmostEqual(self.gruneisenobject_small.acoustic_debye_temp, 317.54811309631845)
        self.assertAlmostEqual(self.gruneisenobject.acoustic_debye_temp, 342.2046198151735)


if __name__ == "__main__":
    unittest.main()
