# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.util.testing import PymatgenTest
from pymatgen.ext.optimade import OptimadeRester


class OptimadeTest(PymatgenTest):
    def test_get_structures_mp(self):

        with OptimadeRester("mp") as optimade:

            structs = optimade.get_structures(elements=["Ga", "N"], nelements=2)

        test_struct = next(iter(structs.values()))

        self.assertEqual([str(el) for el in test_struct.types_of_species], ["Ga", "N"])

    def test_get_structures_mcloud_2dstructures(self):

        with OptimadeRester("mcloud.2dstructures") as optimade:

            structs = optimade.get_structures(elements=["B", "N"], nelements=2)

        test_struct = next(iter(structs.values()))

        self.assertEqual([str(el) for el in test_struct.types_of_species], ["B", "N"])

    def test_update_aliases(self):

        with OptimadeRester() as optimade:
            optimade.refresh_aliases()

        self.assertIn("mp", optimade.aliases)

        from pprint import pprint

        pprint(optimade.aliases)
