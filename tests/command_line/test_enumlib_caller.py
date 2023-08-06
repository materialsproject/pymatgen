from __future__ import annotations

import unittest
from shutil import which

import pytest
from pytest import approx

from pymatgen.command_line.enumlib_caller import EnumError, EnumlibAdaptor
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.site_transformations import RemoveSitesTransformation
from pymatgen.transformations.standard_transformations import SubstitutionTransformation
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

enum_cmd = which("enum.x") or which("multienum.x")
makestr_cmd = which("makestr.x") or which("makeStr.x") or which("makeStr.py")
enumlib_present = enum_cmd and makestr_cmd


@unittest.skipIf(not enumlib_present, "enum_lib not present.")
class TestEnumlibAdaptor(PymatgenTest):
    _multiprocess_shared_ = True

    def test_init(self):
        struct = self.get_structure("LiFePO4")
        subtrans = SubstitutionTransformation({"Li": {"Li": 0.5}})
        adaptor = EnumlibAdaptor(subtrans.apply_transformation(struct), 1, 2)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 86
        for s in structures:
            assert s.composition.get_atomic_fraction(Element("Li")) == approx(0.5 / 6.5)
        adaptor = EnumlibAdaptor(subtrans.apply_transformation(struct), 1, 2, refine_structure=True)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 52

        subtrans = SubstitutionTransformation({"Li": {"Li": 0.25}})
        adaptor = EnumlibAdaptor(subtrans.apply_transformation(struct), 1, 1, refine_structure=True)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 1
        for s in structures:
            assert s.composition.get_atomic_fraction(Element("Li")) == approx(0.25 / 6.25)

        # Make sure it works for completely disordered structures.
        struct = Structure([[10, 0, 0], [0, 10, 0], [0, 0, 10]], [{"Fe": 0.5}], [[0, 0, 0]])
        adaptor = EnumlibAdaptor(struct, 1, 2)
        adaptor.run()
        assert len(adaptor.structures) == 3

        # Make sure it works properly when symmetry is broken by ordered sites.
        struct = self.get_structure("LiFePO4")
        subtrans = SubstitutionTransformation({"Li": {"Li": 0.25}})
        s = subtrans.apply_transformation(struct)
        # REmove some ordered sites to break symmetry.
        removetrans = RemoveSitesTransformation([4, 7])
        s = removetrans.apply_transformation(s)
        adaptor = EnumlibAdaptor(s, 1, 1, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 4

        struct = Structure(
            [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
            [{"Si": 0.5}] * 2,
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        adaptor = EnumlibAdaptor(struct, 1, 3, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 10

        struct = Structure.from_file(f"{TEST_FILES_DIR}/EnumerateTest.json")
        adaptor = EnumlibAdaptor(struct, 1, 1)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 2

    def test_rounding_errors(self):
        # It used to be that a rounding issue would result in this structure
        # showing that Cu3Te2 satisfies an ordering of this structure.
        # This has been fixed by multiplying the base by 100.
        struct = Structure.from_file(f"{TEST_FILES_DIR}/Cu7Te5.cif")
        adaptor = EnumlibAdaptor(struct, 1, 2)
        with pytest.raises(EnumError, match="Unable to enumerate structure"):
            adaptor.run()
        adaptor = EnumlibAdaptor(struct, 1, 5)
        adaptor.run()
        assert len(adaptor.structures) == 197

    def test_partial_disorder(self):
        struct = Structure.from_file(filename=f"{TEST_FILES_DIR}/garnet.cif")
        spga = SpacegroupAnalyzer(struct, 0.1)
        prim = spga.find_primitive()
        struct = prim.copy()
        struct["Al3+"] = {"Al3+": 0.5, "Ga3+": 0.5}
        adaptor = EnumlibAdaptor(struct, 1, 1, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 7
        for struct in structures:
            assert struct.formula == "Ca12 Al4 Ga4 Si12 O48"
        struct = prim.copy()
        struct["Ca2+"] = {"Ca2+": 1 / 3, "Mg2+": 2 / 3}
        adaptor = EnumlibAdaptor(struct, 1, 1, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 20
        for struct in structures:
            assert struct.formula == "Ca4 Mg8 Al8 Si12 O48"

        struct = prim.copy()
        struct["Si4+"] = {"Si4+": 1 / 3, "Ge4+": 2 / 3}
        adaptor = EnumlibAdaptor(struct, 1, 1, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 18
        for struct in structures:
            assert struct.formula == "Ca12 Al8 Si4 Ge8 O48"

    @unittest.skip("Fails seemingly at random.")
    def test_timeout(self):
        struct = Structure.from_file(filename=f"{TEST_FILES_DIR}/garnet.cif")
        SpacegroupAnalyzer(struct, 0.1)
        struct["Al3+"] = {"Al3+": 0.5, "Ga3+": 0.5}
        adaptor = EnumlibAdaptor(struct, 1, 1, enum_precision_parameter=0.01, timeout=0.0000000000001)
        with pytest.raises(TimeoutError, match="Enumeration took too long"):
            adaptor._run_multienum()
