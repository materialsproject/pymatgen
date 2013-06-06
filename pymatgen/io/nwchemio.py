#!/usr/bin/env python

"""
This module implements input and output processing from Nwchem.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "6/5/13"


import re

from pymatgen.core import Molecule
from pymatgen.util.io_utils import zopen
from pymatgen.serializers.json_coders import MSONable


class NwTask(MSONable):
    """
    Base task for Nwchem. Very flexible arguments to support many types of
    potential setups. Users should use more friendly subclasses unless they
    need the flexibility.
    """

    def __init__(self, mol, charge=None, spin_multiplicity=None,
                 title=None, theory="dft", operation="optimize",
                 basis_set="6-31++G**", theory_directives=None):
        self.mol = mol
        self.title = title if title is not None else "{} {} {}".format(
            re.sub("\s", "", mol.formula), theory, operation)

        self.charge = charge if charge is not None else mol.charge
        nelectrons = - self.charge + mol.charge + mol.nelectrons
        if spin_multiplicity is not None:
            self.spin_multiplicity = spin_multiplicity
            if (nelectrons + spin_multiplicity) % 2 != 1:
                raise ValueError(
                    "Charge of {} and spin multiplicity of {} is"
                    " not possible for this molecule".format(
                        charge, spin_multiplicity))
        else:
            self.spin_multiplicity = 1 if nelectrons % 2 == 0 else 2

        elements = set(mol.composition.get_el_amt_dict().keys())

        self.theory = theory
        self.basis_set = basis_set
        self.elements = elements
        self.operation = operation
        self.theory_directives = theory_directives \
            if theory_directives is not None else {}

    def __str__(self):
        o = ["title \"{}\"".format(self.title),
             "charge {}".format(self.charge),
             "basis"]
        for el in self.elements:
            o.append(" {} library {}".format(el, self.basis_set))
        o.append("end")
        if self.theory_directives:
            o.append("{}".format(self.theory))
            for k, v in self.theory_directives.items():
                o.append(" {} {}".format(k, v))
            o.append("end")
        o.append("task {} {}".format(self.theory, self.operation))
        return "\n".join(o)

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__, "mol": self.mol.to_dict,
                "charge": self.charge,
                "spin_multiplicity": self.spin_multiplicity,
                "title": self.title, "theory": self.theory,
                "operation": self.operation, "basis_set": self.basis_set,
                "theory_directives": self.theory_directives}

    @classmethod
    def from_dict(cls, d):
        mol = Molecule.from_dict(d["mol"])
        return NwTask(mol, charge=d["charge"],
                      spin_multiplicity=d["spin_multiplicity"],
                      title=d["title"], theory=d["theory"],
                      operation=d["operation"], basis_set=d["basis_set"],
                      theory_directives=d["theory_directives"])


class DFTTask(NwTask):
    """
    Creates a DFT task from a molecule.
    """

    def __init__(self, mol, xc="b3lyp", charge=None, spin_multiplicity=None,
                 title=None, operation="optimize", basis_set="6-31++G*",
                 theory_directives=None):
        super(DFTTask, self).__init__(
            mol, charge=charge, spin_multiplicity=spin_multiplicity,
            title=title, theory="dft", operation=operation,
            basis_set=basis_set, theory_directives=theory_directives)

        self.theory_directives.update({"xc": xc,
                                       "mult": self.spin_multiplicity})


class NwInput(MSONable):
    """
    An object representing a Nwchem input file.
    """

    def __init__(self, mol, tasks, directives=None):
        """
        Args:
            mol:
                Input molecule. If molecule is a single string, it is used as a
                direct input to the geometry section of the Gaussian input
                file.
            tasks:
                List of NwTasks.
            directives:
                List of root level directives as tuple. E.g.,
                [("start", "water"), ("print", "high")]
        """
        self._mol = mol
        if directives is None:
            self.directives = [("start", re.sub("\s", "", self._mol.formula))]
        else:
            self.directives = directives
        self.tasks = tasks

    @property
    def molecule(self):
        """
        Returns molecule associated with this GaussianInput.
        """
        return self._mol

    def __str__(self):
        o = []
        for d in self.directives:
            o.append("{} {}".format(d[0], d[1]))
        o.append("geometry units angstroms")
        for site in self._mol:
            o.append(" {} {} {} {}".format(site.specie.symbol, site.x, site.y,
                                           site.z))
        o.append("end\n")
        for t in self.tasks:
            o.append(str(t))
            o.append("")
        return "\n".join(o)

    def write_file(self, filename):
        with zopen(filename, "w") as f:
            f.write(self.__str__())

    @property
    def to_dict(self):
        return {
            "mol": self._mol.to_dict,
            "tasks": [t.to_dict for t in self.tasks],
            "directives": self.directives
        }

    @classmethod
    def from_dict(cls, d):
        return NwInput(Molecule.from_dict(d["mol"]),
                       [NwTask.from_dict(dt) for dt in d["tasks"]],
                       d["directives"])


import unittest
import os
import json


test_dir = os.path.join(os.path.dirname(__file__), "..", "..",
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


class DFTTaskTest(unittest.TestCase):

    def setUp(self):
        self.task = DFTTask(mol, charge=1, operation="energy")

    def test_str_and_from_string(self):
        ans = """title "H4C1 dft energy"
charge 1
basis
 H library 6-31++G*
 C library 6-31++G*
end
dft
 xc b3lyp
 mult 2
end
task dft energy"""
        self.assertEqual(str(self.task), ans)


class NwInputTest(unittest.TestCase):

    def setUp(self):
        tasks = [
            DFTTask(mol, operation="optimize", xc="b3lyp",
                    basis_set="6-31++G*"),
            DFTTask(mol, operation="freq", xc="b3lyp", basis_set="6-31++G*"),
            DFTTask(mol, operation="energy", xc="b3lyp",
                    basis_set="6-311++G**"),
            DFTTask(mol, charge=mol.charge + 1, operation="energy", xc="b3lyp",
                    basis_set="6-311++G**"),
            DFTTask(mol, charge=mol.charge - 1, operation="energy", xc="b3lyp",
                    basis_set="6-311++G**")
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
        json.dumps(d)


if __name__ == "__main__":
    unittest.main()
