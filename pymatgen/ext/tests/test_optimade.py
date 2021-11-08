# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.util.testing import PymatgenTest
from pymatgen.ext.optimade import OptimadeRester


class OptimadeTest(PymatgenTest):
    def test_get_structures_mp(self):

        with OptimadeRester("mp") as optimade:

            structs = optimade.get_structures(elements=["Ga", "N"], nelements=2)

        with OptimadeRester("mp") as optimade:

            _filter = 'elements HAS ALL "Ga", "N" AND nelements=2'
            raw_filter_structs = optimade.get_structures_with_filter(_filter)

            # skip if query fails to return any results (e.g. server down or inaccessible)
            if ("mp" in structs) and ("mp" in raw_filter_structs):

                test_struct = next(iter(structs["mp"].values()))
                self.assertEqual([str(el) for el in test_struct.types_of_species], ["Ga", "N"])

                self.assertEqual(
                    len(structs["mp"]),
                    len(raw_filter_structs["mp"]),
                    msg="Raw filter {_filter} did not return the same number of results as the query builder.",
                )

    # Tests fail in CI for unknown reason, use for development only.
    # def test_get_structures_mcloud_2dstructures(self):
    #
    #     with OptimadeRester("mcloud.2dstructures") as optimade:
    #
    #         structs = optimade.get_structures(elements=["B", "N"], nelements=2)
    #
    #     test_struct = next(iter(structs["mcloud.2dstructures"].values()))
    #
    #     self.assertEqual([str(el) for el in test_struct.types_of_species], ["B", "N"])

    # def test_update_aliases(self):
    #
    #     with OptimadeRester() as optimade:
    #         optimade.refresh_aliases()
    #
    #     self.assertIn("mp", optimade.aliases)
