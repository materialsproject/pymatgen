#!/usr/bin/env python

"""
TODO: Modify module doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "6/6/13"


import unittest
import os
import json

from pymatgen.core.structure import Molecule
from pymatgen.io.nwchemio import NwTask, NwInput, NwInputError, NwOutput


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules")


coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]
mol = Molecule(["C", "H", "H", "H", "H"], coords)


class NwTaskTest(unittest.TestCase):

    def setUp(self):
        self.task = NwTask(mol, theory_directives={"xc": "b3lyp"})

    def test_str_and_from_string(self):
        ans = """title "H4C1 dft optimize"
charge 0
basis
 H library 6-31++G**
 C library 6-31++G**
end
dft
 xc b3lyp
end
task dft optimize"""
        self.assertEqual(str(self.task), ans)

    def test_to_from_dict(self):
        d = self.task.to_dict
        t = NwTask.from_dict(d)
        self.assertIsInstance(t, NwTask)

    def test_init(self):
        self.assertRaises(NwInputError, NwTask, mol, theory="bad")
        self.assertRaises(NwInputError, NwTask, mol, operation="bad")

    def test_dft_task(self):
        task = NwTask.dft_task(mol, charge=1, operation="energy")
        ans = """title "H4C1 dft energy"
charge 1
basis
 H library 6-31++G**
 C library 6-31++G**
end
dft
 xc b3lyp
 mult 2
end
task dft energy"""

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
        self.nwi = NwInput(mol, tasks)

    def test_str(self):
        ans = """start H4C1
geometry units angstroms
 C 0.0 0.0 0.0
 H 0.0 0.0 1.089
 H 1.026719 0.0 -0.363
 H -0.51336 -0.889165 -0.363
 H -0.51336 0.889165 -0.363
end

title "H4C1 dft optimize"
charge 0
basis
 H library 6-31++G*
 C library 6-31++G*
end
dft
 xc b3lyp
 mult 1
end
task dft optimize

title "H4C1 dft freq"
charge 0
basis
 H library 6-31++G*
 C library 6-31++G*
end
dft
 xc b3lyp
 mult 1
end
task dft freq

title "H4C1 dft energy"
charge 0
basis
 H library 6-311++G**
 C library 6-311++G**
end
dft
 xc b3lyp
 mult 1
end
task dft energy

title "H4C1 dft energy"
charge 1
basis
 H library 6-311++G**
 C library 6-311++G**
end
dft
 xc b3lyp
 mult 2
end
task dft energy

title "H4C1 dft energy"
charge -1
basis
 H library 6-311++G**
 C library 6-311++G**
end
dft
 xc b3lyp
 mult 2
end
task dft energy
"""
        self.assertEqual(str(self.nwi), ans)

    def test_to_from_dict(self):
        d = self.nwi.to_dict
        nwi = NwInput.from_dict(d)
        self.assertIsInstance(nwi, NwInput)
        #Ensure it is json-serializable.
        json.dumps(d)


class NwOutputTest(unittest.TestCase):

    def test_read(self):
        nwo = NwOutput(os.path.join(test_dir, "CH4.nwout"))
        print nwo.data

