# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os
import json

from pymatgen.core.structure import Molecule
from pymatgen.io.nwchem import NwTask, NwInput, NwInputError, NwOutput

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "nwchem")

coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]
mol = Molecule(["C", "H", "H", "H", "H"], coords)


class NwTaskTest(unittest.TestCase):

    def setUp(self):
        self.task = NwTask(0, 1, basis_set={"H": "6-31g"}, theory="dft",
                           theory_directives={"xc": "b3lyp"})
        self.task_cosmo = NwTask(0, 1, basis_set={"H": "6-31g"}, theory="dft",
                                 theory_directives={"xc": "b3lyp"},
                                 alternate_directives={'cosmo': "cosmo"})
        self.task_esp = NwTask(0, 1, basis_set={"H": "6-31g"}, theory="esp")

    def test_multi_bset(self):
        t = NwTask.from_molecule(
            mol, theory="dft", basis_set={"C": "6-311++G**",
                                          "H": "6-31++G**"},
            theory_directives={"xc": "b3lyp"})
        ans = """title "H4C1 dft optimize"
charge 0
basis cartesian
 C library "6-311++G**"
 H library "6-31++G**"
end
dft
 xc b3lyp
end
task dft optimize"""
        self.assertEqual(str(t), ans)

    def test_str_and_from_string(self):
        ans = """title "dft optimize"
charge 0
basis cartesian
 H library "6-31g"
end
dft
 xc b3lyp
end
task dft optimize"""
        self.assertEqual(str(self.task), ans)

    def test_to_from_dict(self):
        d = self.task.as_dict()
        t = NwTask.from_dict(d)
        self.assertIsInstance(t, NwTask)

    def test_init(self):
        self.assertRaises(NwInputError, NwTask, 0, 1, {"H": "6-31g"},
                          theory="bad")
        self.assertRaises(NwInputError, NwTask, 0, 1, {"H": "6-31g"},
                          operation="bad")

    def test_dft_task(self):
        task = NwTask.dft_task(mol, charge=1, operation="energy")
        ans = """title "H4C1 dft energy"
charge 1
basis cartesian
 C library "6-31g"
 H library "6-31g"
end
dft
 mult 2
 xc b3lyp
end
task dft energy"""
        self.assertEqual(str(task), ans)

    def test_dft_cosmo_task(self):
        task = NwTask.dft_task(
            mol, charge=mol.charge, operation="energy",
            xc="b3lyp", basis_set="6-311++G**",
            alternate_directives={'cosmo': {"dielec": 78.0}})
        ans = """title "H4C1 dft energy"
charge 0
basis cartesian
 C library "6-311++G**"
 H library "6-311++G**"
end
dft
 mult 1
 xc b3lyp
end
cosmo
 dielec 78.0
end
task dft energy"""
        self.assertEqual(str(task), ans)

    def test_esp_task(self):
        task = NwTask.esp_task(mol, charge=mol.charge, operation="",
                               basis_set="6-311++G**")
        ans = """title "H4C1 esp "
charge 0
basis cartesian
 C library "6-311++G**"
 H library "6-311++G**"
end

task esp """
        self.assertEqual(str(task), ans)


class NwInputTest(unittest.TestCase):
    def setUp(self):
        tasks = [
            NwTask.dft_task(mol, operation="optimize", xc="b3lyp",
                            basis_set="6-31++G*"),
            NwTask.dft_task(mol, operation="freq", xc="b3lyp",
                            basis_set="6-31++G*"),
            NwTask.dft_task(mol, operation="energy", xc="b3lyp",
                            basis_set="6-311++G**"),
            NwTask.dft_task(mol, charge=mol.charge + 1, operation="energy",
                            xc="b3lyp", basis_set="6-311++G**"),
            NwTask.dft_task(mol, charge=mol.charge - 1, operation="energy",
                            xc="b3lyp", basis_set="6-311++G**")
        ]

        self.nwi = NwInput(mol, tasks,
                           geometry_options=["units", "angstroms", "noautoz"],
                           memory_options="total 1000 mb")
        self.nwi_symm = NwInput(mol, tasks,
                                geometry_options=["units", "angstroms",
                                                  "noautoz"],
                                symmetry_options=["c1"])

    def test_str(self):
        ans = """memory total 1000 mb
geometry units angstroms noautoz
 C 0.0 0.0 0.0
 H 0.0 0.0 1.089
 H 1.026719 0.0 -0.363
 H -0.51336 -0.889165 -0.363
 H -0.51336 0.889165 -0.363
end

title "H4C1 dft optimize"
charge 0
basis cartesian
 C library "6-31++G*"
 H library "6-31++G*"
end
dft
 mult 1
 xc b3lyp
end
task dft optimize

title "H4C1 dft freq"
charge 0
basis cartesian
 C library "6-31++G*"
 H library "6-31++G*"
end
dft
 mult 1
 xc b3lyp
end
task dft freq

title "H4C1 dft energy"
charge 0
basis cartesian
 C library "6-311++G**"
 H library "6-311++G**"
end
dft
 mult 1
 xc b3lyp
end
task dft energy

title "H4C1 dft energy"
charge 1
basis cartesian
 C library "6-311++G**"
 H library "6-311++G**"
end
dft
 mult 2
 xc b3lyp
end
task dft energy

title "H4C1 dft energy"
charge -1
basis cartesian
 C library "6-311++G**"
 H library "6-311++G**"
end
dft
 mult 2
 xc b3lyp
end
task dft energy
"""
        self.assertEqual(str(self.nwi), ans)

        ans_symm = """geometry units angstroms noautoz
 symmetry c1
 C 0.0 0.0 0.0
 H 0.0 0.0 1.089
 H 1.026719 0.0 -0.363
 H -0.51336 -0.889165 -0.363
 H -0.51336 0.889165 -0.363
end

title "H4C1 dft optimize"
charge 0
basis cartesian
 C library "6-31++G*"
 H library "6-31++G*"
end
dft
 mult 1
 xc b3lyp
end
task dft optimize

title "H4C1 dft freq"
charge 0
basis cartesian
 C library "6-31++G*"
 H library "6-31++G*"
end
dft
 mult 1
 xc b3lyp
end
task dft freq

title "H4C1 dft energy"
charge 0
basis cartesian
 C library "6-311++G**"
 H library "6-311++G**"
end
dft
 mult 1
 xc b3lyp
end
task dft energy

title "H4C1 dft energy"
charge 1
basis cartesian
 C library "6-311++G**"
 H library "6-311++G**"
end
dft
 mult 2
 xc b3lyp
end
task dft energy

title "H4C1 dft energy"
charge -1
basis cartesian
 C library "6-311++G**"
 H library "6-311++G**"
end
dft
 mult 2
 xc b3lyp
end
task dft energy
"""

        self.assertEqual(str(self.nwi_symm), ans_symm)

    def test_to_from_dict(self):
        d = self.nwi.as_dict()
        nwi = NwInput.from_dict(d)
        self.assertIsInstance(nwi, NwInput)
        # Ensure it is json-serializable.
        json.dumps(d)
        d = self.nwi_symm.as_dict()
        nwi_symm = NwInput.from_dict(d)
        self.assertIsInstance(nwi_symm, NwInput)
        json.dumps(d)

    def test_from_string_and_file(self):
        nwi = NwInput.from_file(os.path.join(test_dir, "ch4.nw"))
        self.assertEqual(nwi.tasks[0].theory, "dft")
        self.assertEqual(nwi.memory_options, "total 1000 mb stack 400 mb")
        self.assertEqual(nwi.tasks[0].basis_set["C"], "6-31++G*")
        self.assertEqual(nwi.tasks[-1].basis_set["C"], "6-311++G**")
        # Try a simplified input.
        str_inp = """start H4C1
geometry units angstroms
 C 0.0 0.0 0.0
 H 0.0 0.0 1.089
 H 1.026719 0.0 -0.363
 H -0.51336 -0.889165 -0.363
 H -0.51336 0.889165 -0.363
end

title "H4C1 dft optimize"
charge 0
basis cartesian
 H library "6-31++G*"
 C library "6-31++G*"
end
dft
 xc b3lyp
 mult 1
end
task scf optimize

title "H4C1 dft freq"
charge 0
task scf freq

title "H4C1 dft energy"
charge 0
basis cartesian
 H library "6-311++G**"
 C library "6-311++G**"
end
task dft energy

title "H4C1 dft energy"
charge 1
dft
 xc b3lyp
 mult 2
end
task dft energy

title "H4C1 dft energy"
charge -1
task dft energy
"""
        nwi = NwInput.from_string(str_inp)
        self.assertEqual(nwi.geometry_options, ['units', 'angstroms'])
        self.assertEqual(nwi.tasks[0].theory, "scf")
        self.assertEqual(nwi.tasks[0].basis_set["C"], "6-31++G*")
        self.assertEqual(nwi.tasks[-1].theory, "dft")
        self.assertEqual(nwi.tasks[-1].basis_set["C"], "6-311++G**")

        str_inp_symm = str_inp.replace("geometry units angstroms",
                                       "geometry units angstroms\n symmetry "
                                       "c1")

        nwi_symm = NwInput.from_string(str_inp_symm)
        self.assertEqual(nwi_symm.geometry_options, ['units', 'angstroms'])
        self.assertEqual(nwi_symm.symmetry_options, ['c1'])
        self.assertEqual(nwi_symm.tasks[0].theory, "scf")
        self.assertEqual(nwi_symm.tasks[0].basis_set["C"], "6-31++G*")
        self.assertEqual(nwi_symm.tasks[-1].theory, "dft")
        self.assertEqual(nwi_symm.tasks[-1].basis_set["C"], "6-311++G**")


class NwOutputTest(unittest.TestCase):
    def test_read(self):
        nwo = NwOutput(os.path.join(test_dir, "CH4.nwout"))
        nwo_cosmo = NwOutput(os.path.join(test_dir, "N2O4.nwout"))

        self.assertEqual(0, nwo[0]["charge"])
        self.assertEqual(-1, nwo[-1]["charge"])
        self.assertEqual(len(nwo), 5)
        self.assertAlmostEqual(-1102.6224491715582, nwo[0]["energies"][-1], 2)
        self.assertAlmostEqual(-1102.9986291578023, nwo[2]["energies"][-1], 3)
        self.assertAlmostEqual(-11156.354030653656,
                               nwo_cosmo[5]["energies"][0]["cosmo scf"], 3)
        self.assertAlmostEqual(-11153.374133394364,
                               nwo_cosmo[5]["energies"][0]["gas phase"], 3)
        self.assertAlmostEqual(-11156.353632962995,
                               nwo_cosmo[5]["energies"][0]["sol phase"], 2)
        self.assertAlmostEqual(-11168.818934311605,
                               nwo_cosmo[6]["energies"][0]["cosmo scf"], 2)
        self.assertAlmostEqual(-11166.3624424611462,
                               nwo_cosmo[6]["energies"][0]['gas phase'], 2)
        self.assertAlmostEqual(-11168.818934311605,
                               nwo_cosmo[6]["energies"][0]['sol phase'], 2)
        self.assertAlmostEqual(-11165.227959110889,
                               nwo_cosmo[7]["energies"][0]['cosmo scf'], 2)
        self.assertAlmostEqual(-11165.025443612385,
                               nwo_cosmo[7]["energies"][0]['gas phase'], 2)
        self.assertAlmostEqual(-11165.227959110154,
                               nwo_cosmo[7]["energies"][0]['sol phase'], 2)

        self.assertAlmostEqual(nwo[1]["hessian"][0][0], 4.60187e+01)
        self.assertAlmostEqual(nwo[1]["hessian"][1][2], -1.14030e-08)
        self.assertAlmostEqual(nwo[1]["hessian"][2][3], 2.60819e+01)
        self.assertAlmostEqual(nwo[1]["hessian"][6][6], 1.45055e+02)
        self.assertAlmostEqual(nwo[1]["hessian"][11][14], 1.35078e+01)

        # CH4.nwout, line 722
        self.assertAlmostEqual(nwo[0]["forces"][0][3], -0.001991)

        # N2O4.nwout, line 1071
        self.assertAlmostEqual(nwo_cosmo[0]["forces"][0][4], 0.011948)

        # There should be four DFT gradients.
        self.assertEqual(len(nwo_cosmo[0]["forces"]), 4)

        ie = (nwo[4]["energies"][-1] - nwo[2]["energies"][-1])
        ea = (nwo[2]["energies"][-1] - nwo[3]["energies"][-1])
        self.assertAlmostEqual(0.7575358648355177, ie)
        self.assertAlmostEqual(-14.997877958701338, ea, 3)
        self.assertEqual(nwo[4]["basis_set"]["C"]["description"],
                         "6-311++G**")

        nwo = NwOutput(os.path.join(test_dir, "H4C3O3_1.nwout"))
        self.assertTrue(nwo[-1]["has_error"])
        self.assertEqual(nwo[-1]["errors"][0], "Bad convergence")

        nwo = NwOutput(os.path.join(test_dir, "CH3CH2O.nwout"))
        self.assertTrue(nwo[-1]["has_error"])
        self.assertEqual(nwo[-1]["errors"][0], "Bad convergence")

        nwo = NwOutput(os.path.join(test_dir, "C1N1Cl1_1.nwout"))
        self.assertTrue(nwo[-1]["has_error"])
        self.assertEqual(nwo[-1]["errors"][0], "autoz error")

        nwo = NwOutput(os.path.join(test_dir,
                                    "anthrachinon_wfs_16_ethyl.nwout"))
        self.assertTrue(nwo[-1]["has_error"])
        self.assertEqual(nwo[-1]["errors"][0],
                         "Geometry optimization failed")
        nwo = NwOutput(os.path.join(test_dir,
                                    "anthrachinon_wfs_15_carboxyl.nwout"))
        self.assertEqual(nwo[1]['frequencies'][0][0], -70.47)
        self.assertEqual(len(nwo[1]['frequencies'][0][1]), 27)
        self.assertEqual(nwo[1]['frequencies'][-1][0], 3696.74)
        self.assertEqual(nwo[1]['frequencies'][-1][1][-1],
                         (0.20498, -0.94542, -0.00073))
        self.assertEqual(nwo[1]["normal_frequencies"][1][0], -70.72)
        self.assertEqual(nwo[1]["normal_frequencies"][3][0], -61.92)
        self.assertEqual(nwo[1]["normal_frequencies"][1][1][-1],
                         (0.00056, 0.00042, 0.06781))

    def test_parse_tddft(self):
        nwo = NwOutput(os.path.join(test_dir, "phen_tddft.log"))
        roots = nwo.parse_tddft()
        self.assertEqual(len(roots["singlet"]), 20)
        self.assertAlmostEqual(roots["singlet"][0]["energy"], 3.9291)
        self.assertAlmostEqual(roots["singlet"][0]["osc_strength"], 0.0)
        self.assertAlmostEqual(roots["singlet"][1]["osc_strength"], 0.00177)

    def test_get_excitation_spectrum(self):
        nwo = NwOutput(os.path.join(test_dir, "phen_tddft.log"))
        spectrum = nwo.get_excitation_spectrum()
        self.assertEqual(len(spectrum.x), 2000)
        self.assertAlmostEqual(spectrum.x[0], 1.9291)
        self.assertAlmostEqual(spectrum.y[0], 0.0)
        self.assertAlmostEqual(spectrum.y[1000], 0.0007423569947114812)


if __name__ == "__main__":
    unittest.main()
