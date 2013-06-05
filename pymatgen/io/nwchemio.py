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


class NwTask(object):

    def __init__(self, elements, title=None, theory="dft",
                 operation="optimize", basis_set="6-31++G**",
                 theory_directives=None):
        self.title = title if title is not None else "{} {} {}".format(
            "".join(map(str, elements)), theory, operation)
        self.theory = theory
        self.basis_set = basis_set
        self.elements = elements
        self.operation = operation
        self.directives = theory_directives if theory_directives is not None \
            else {}

    def __str__(self):
        o = ["title \"{}\"".format(self.title),
             "basis"]
        for el in self.elements:
            o.append(" {} library {}".format(el, self.basis_set))
        o.append("end")
        if self.directives:
            o.append("{}".format(self.theory))
            for k, v in self.directives.items():
                o.append(" {} {}".format(k, v))
            o.append("end")
        o.append("task {} {}".format(self.theory, self.operation))
        return "\n".join(o)



class NwInput(object):
    """
    An object representing a Nwchem input file.
    """

    def __init__(self, mol, tasks, directives=None, charge=None,
                 spin_multiplicity=None):
        """
        Args:
            mol:
                Input molecule. If molecule is a single string, it is used as a
                direct input to the geometry section of the Gaussian input
                file.
            tasks:
                List of NwTasks.
            job_title:
                Title of job. Typically
            charge:
                Charge of the molecule. If None, charge on molecule is used.
                Defaults to None. This allows the input file to be set a
                charge independently from the molecule itself.
            spin_multiplicity:
                Spin multiplicity of molecule. Defaults to None,
                which means that the spin multiplicity is set to 1 if the
                molecule has no unpaired electrons and to 2 if there are
                unpaired electrons.
        """
        self._mol = mol
        self.charge = charge if charge is not None else mol.charge
        nelectrons = - self.charge + mol.charge + mol.nelectrons
        if spin_multiplicity is not None:
            self.spin_multiplicity = spin_multiplicity
            if (nelectrons + spin_multiplicity) % 2 != 1:
                raise ValueError(
                    "Charge of {} and spin multiplicity of {} is"
                    " not possible for this molecule".format(
                        self.charge, spin_multiplicity))
        else:
            self.spin_multiplicity = 1 if nelectrons % 2 == 0 else 2
        self.tasks = tasks
        if self.spin_multiplicity != 1:
            for t in self.tasks:
                if t.theory.lower() == "dft":
                    t.directives["mult"] = self.spin_multiplicity

        if directives is None:
            self.directives = [("start", re.sub("\s", "", self._mol.formula))]
        else:
            self.directives = self.directives

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
        o.append("charge {}".format(self.charge))
        o.append("geometry units angstroms")
        for site in self._mol:
            o.append(" {} {} {} {}".format(site.specie.symbol, site.x, site.y,
                                           site.z))
        o.append("end")
        for t in self.tasks:
            o.append(str(t))
        return "\n".join(o)

    def write_file(self, filename):
        with zopen(filename, "w") as f:
            f.write(self.__str__())


import unittest
import os


test_dir = os.path.join(os.path.dirname(__file__), "..", "..",
                        'test_files', "molecules")


class NwTaskTest(unittest.TestCase):

    def setUp(self):
        self.task = NwTask(["H", "O"], theory_directives={"xc": "b3lyp"})

    def test_str_and_from_string(self):
        ans = """title "HO dft optimize"
basis
 H library 6-31++G**
 O library 6-31++G**
end
dft
 xc b3lyp
end
task dft optimize"""
        self.assertEqual(str(self.task), ans)


class NwInputTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        mol = Molecule(["C", "H", "H", "H", "H"], coords)
        tasks = [
            NwTask(["H", "C"], operation="optimize",
                   theory_directives={"xc": "b3lyp"}),
            NwTask(["H", "C"], operation="freq",
                   theory_directives={"xc": "b3lyp"})
        ]
        self.nwi = NwInput(
            mol, tasks, charge=1)

    def test_str(self):
        print self.nwi
        self.nwi.write_file("CH4.nw")


if __name__ == "__main__":
    unittest.main()
