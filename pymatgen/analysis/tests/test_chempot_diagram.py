# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import warnings
from pathlib import Path
from pymatgen.core.composition import Element
from pymatgen.util.testing import PymatgenTest
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.analysis.chempot_diagram import ChemicalPotentialDiagram


module_dir = Path(__file__).absolute().parent


class ChemicalPotentialDiagramTest(PymatgenTest):
    def setUp(self):
        self.entries = EntrySet.from_csv(str(module_dir / "pdentries_test.csv"))
        self.cpd = ChemicalPotentialDiagram(entries=self.entries,
                                            default_min_limit=-25)
        print(self.cpd)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_domains(self):
        pass

    def test_dim(self):
        pass

    def test_2d_plot(self):
        pass

    def test_3d_plot(self):
        pass

    def test_border_hyperplanes(self):
        pass

    def test_hyperplanes_and_entries(self):
        pass

    def test_get_min_entries_and_el_refs(self):
        pass

    def test_domain_simplices(self):
        pass

    def test_lims(self):
        pass

    def test_el_refs(self):
        el_refs = {elem: entry.energy for elem, entry in
                   self.cpd.el_refs.items()}
        print(el_refs)

        elems = [Element("Li"), Element("Fe"), Element("O")]
        energies = [-1.89932649, -6.44218495, -8.49740766]
        correct_el_refs = dict(zip(elems, energies))

        self.assertDictsAlmostEqual(el_refs, correct_el_refs)

    def test_pca(self):
        pass

    def test_centroid(self):
        pass


if __name__ == "__main__":
    unittest.main()