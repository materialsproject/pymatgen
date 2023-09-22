from __future__ import annotations

import unittest

import requests

from pymatgen.core import SETTINGS
from pymatgen.ext.optimade import OptimadeRester
from pymatgen.util.testing import PymatgenTest

try:
    website_down = requests.get("https://materialsproject.org").status_code != 200
except requests.exceptions.ConnectionError:
    website_down = True


class TestOptimade(PymatgenTest):
    @unittest.skipIf(
        not SETTINGS.get("PMG_MAPI_KEY") or website_down,
        "PMG_MAPI_KEY environment variable not set or MP is down.",
    )
    def test_get_structures_mp(self):
        with OptimadeRester("mp") as optimade:
            structs = optimade.get_structures(elements=["Ga", "N"], nelements=2)

        with OptimadeRester("mp") as optimade:
            _filter = 'elements HAS ALL "Ga", "N" AND nelements=2'
            raw_filter_structs = optimade.get_structures_with_filter(_filter)

            # skip if query fails to return any results (e.g. server down or inaccessible)
            if ("mp" in structs) and ("mp" in raw_filter_structs):
                test_struct = next(iter(structs["mp"].values()))
                assert [str(el) for el in test_struct.types_of_species] == ["Ga", "N"]

                assert len(structs["mp"]) == len(
                    raw_filter_structs["mp"]
                ), f"Raw filter {_filter} did not return the same number of results as the query builder."

    @unittest.skipIf(
        not SETTINGS.get("PMG_MAPI_KEY") or website_down,
        "PMG_MAPI_KEY environment variable not set or MP is down.",
    )
    def test_get_snls_mp(self):
        base_query = dict(elements=["Ga", "N"], nelements=2, nsites=[2, 6])
        with OptimadeRester("mp") as optimade:
            structs = optimade.get_snls(**base_query)

        field_set = {"nsites", "nelements"}
        with OptimadeRester("mp") as optimade:
            extra_fields_single = optimade.get_snls(**base_query, additional_response_fields="nsites")
            extra_fields_set = optimade.get_snls(**base_query, additional_response_fields=field_set)

        if "mp" in extra_fields_single and "mp" in structs:
            assert len(structs["mp"]) == len(extra_fields_single["mp"])
            struct_nl = next(iter(extra_fields_single["mp"].values()))
            assert "nsites" in struct_nl.data["_optimade"]

        if "mp" in extra_fields_set and "mp" in structs:
            assert len(structs["mp"]) == len(extra_fields_set["mp"])

            # Check that the requested response fields appear in the SNL metadata
            struct_nl_set = next(iter(extra_fields_set["mp"].values()))
            assert field_set <= {*struct_nl_set.data["_optimade"]}

    # Tests fail in CI for unknown reason, use for development only.
    # def test_get_structures_mcloud_2dstructures(self):
    #     with OptimadeRester("mcloud.2dstructures") as optimade:
    #         structs = optimade.get_structures(elements=["B", "N"], nelements=2)

    #     test_struct = next(iter(structs["mcloud.2dstructures"].values()))

    #     assert [str(el) for el in test_struct.types_of_species] == ["B", "N"]

    # def test_update_aliases(self):
    #
    #     with OptimadeRester() as optimade:
    #         optimade.refresh_aliases()
    #
    #     self.assertIn("mp", optimade.aliases)

    def test_build_filter(self):
        with OptimadeRester("mp") as optimade:
            assert optimade._build_filter(
                elements=["Ga", "N"],
                nelements=2,
                nsites=(1, 100),
                chemical_formula_anonymous="A2B",
                chemical_formula_hill="GaN",
            ) == (
                '(elements HAS ALL "Ga", "N")'
                " AND (nsites>=1 AND nsites<=100)"
                " AND (nelements=2)"
                " AND (chemical_formula_anonymous='A2B')"
                " AND (chemical_formula_hill='GaN')"
            )

            assert optimade._build_filter(
                elements=["C", "H", "O"],
                nelements=(3, 4),
                nsites=(1, 100),
                chemical_formula_anonymous="A4B3C",
                chemical_formula_hill="C4H3O",
            ) == (
                '(elements HAS ALL "C", "H", "O")'
                " AND (nsites>=1 AND nsites<=100)"
                " AND (nelements>=3 AND nelements<=4)"
                " AND (chemical_formula_anonymous='A4B3C')"
                " AND (chemical_formula_hill='C4H3O')"
            )
