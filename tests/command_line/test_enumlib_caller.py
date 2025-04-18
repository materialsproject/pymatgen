from __future__ import annotations

from shutil import which

import numpy as np
import pytest
from pytest import approx

from pymatgen.command_line.enumlib_caller import EnumError, EnumlibAdaptor
from pymatgen.core import Element, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.site_transformations import RemoveSitesTransformation
from pymatgen.transformations.standard_transformations import SubstitutionTransformation
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

ENUM_CMD = which("enum.x") or which("multienum.x")
MAKESTR_CMD = which("makestr.x") or which("makeStr.x") or which("makeStr.py")

ENUMLIB_TEST_FILES_DIR: str = f"{TEST_FILES_DIR}/command_line/enumlib"


@pytest.mark.skipif(not (ENUM_CMD and MAKESTR_CMD), reason="enumlib not present.")
class TestEnumlibAdaptor(MatSciTest):
    def test_init(self):
        struct = self.get_structure("LiFePO4")
        sub_trans = SubstitutionTransformation({"Li": {"Li": 0.5}})
        adaptor = EnumlibAdaptor(sub_trans.apply_transformation(struct), 1, 2)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 86
        for struct_trafo in structures:
            assert struct_trafo.composition.get_atomic_fraction(Element("Li")) == approx(0.5 / 6.5)
        adaptor = EnumlibAdaptor(sub_trans.apply_transformation(struct), 1, 2, refine_structure=True)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 52

        sub_trans = SubstitutionTransformation({"Li": {"Li": 0.25}})
        adaptor = EnumlibAdaptor(sub_trans.apply_transformation(struct), refine_structure=True)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 1
        for struct_trafo in structures:
            assert struct_trafo.composition.get_atomic_fraction(Element("Li")) == approx(0.25 / 6.25)

        # Make sure it works for completely disordered structures.
        struct = Structure(np.eye(3) * 10, [{"Fe": 0.5}], [[0, 0, 0]])
        adaptor = EnumlibAdaptor(struct, 1, 2)
        adaptor.run()
        assert len(adaptor.structures) == 3

        # Make sure it works properly when symmetry is broken by ordered sites.
        struct = self.get_structure("LiFePO4")
        sub_trans = SubstitutionTransformation({"Li": {"Li": 0.25}})
        struct_trafo = sub_trans.apply_transformation(struct)
        # REmove some ordered sites to break symmetry.
        remove_trans = RemoveSitesTransformation([4, 7])
        struct_trafo = remove_trans.apply_transformation(struct_trafo)
        adaptor = EnumlibAdaptor(struct_trafo, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 4

        struct = Structure(
            np.eye(3) * 3,
            [{"Si": 0.5}] * 2,
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        adaptor = EnumlibAdaptor(struct, 1, 3, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 10

        struct = Structure.from_file(f"{ENUMLIB_TEST_FILES_DIR}/EnumerateTest.json.gz")
        adaptor = EnumlibAdaptor(struct)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 2

    def test_rounding_errors(self):
        # It used to be that a rounding issue would result in this structure
        # showing that Cu3Te2 satisfies an ordering of this structure.
        # This has been fixed by multiplying the base by 100.
        struct = Structure.from_file(f"{TEST_FILES_DIR}/cif/Cu7Te5.cif")
        adaptor = EnumlibAdaptor(struct, 1, 2)
        with pytest.raises(EnumError, match="Unable to enumerate structure"):
            adaptor.run()
        adaptor = EnumlibAdaptor(struct, 1, 5)
        adaptor.run()
        assert len(adaptor.structures) == 197

    def test_partial_disorder(self):
        struct = Structure.from_file(filename=f"{TEST_FILES_DIR}/cif/garnet.cif")
        spga = SpacegroupAnalyzer(struct, 0.1)
        prim = spga.find_primitive()
        struct = prim.copy()
        struct["Al3+"] = {"Al3+": 0.5, "Ga3+": 0.5}
        adaptor = EnumlibAdaptor(struct, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 7
        for struct in structures:
            assert struct.formula == "Ca12 Al4 Ga4 Si12 O48"
        struct = prim.copy()
        struct["Ca2+"] = {"Ca2+": 1 / 3, "Mg2+": 2 / 3}
        adaptor = EnumlibAdaptor(struct, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 20
        for struct in structures:
            assert struct.formula == "Ca4 Mg8 Al8 Si12 O48"

        struct = prim.copy()
        struct["Si4+"] = {"Si4+": 1 / 3, "Ge4+": 2 / 3}
        adaptor = EnumlibAdaptor(struct, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        assert len(structures) == 18
        for struct in structures:
            assert struct.formula == "Ca12 Al8 Si4 Ge8 O48"

    def test_timeout(self):
        struct = Structure.from_file(f"{ENUMLIB_TEST_FILES_DIR}/test_timeout.json.gz")

        adaptor = EnumlibAdaptor(struct, max_cell_size=10, timeout=0.05)  # timeout in minute

        with pytest.raises(TimeoutError, match="Enumeration took more than timeout 0.05 minutes"):
            adaptor.run()
