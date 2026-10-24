"""Tests for borg."""

from __future__ import annotations

import gzip
import os
import shutil
from typing import TYPE_CHECKING

import orjson
import pytest
from monty.json import MontyDecoder
from pytest import approx

from pymatgen.apps.borg.hive import (
    GaussianToComputedEntryDrone,
    SimpleVaspToComputedEntryDrone,
    VaspToComputedEntryDrone,
    _get_transformation_history,
)
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.vasp.inputs import Incar, Poscar
from tests.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from pathlib import Path

TEST_DIR = f"{TEST_FILES_DIR}/apps/borg"
VASP_IN_DIR = f"{TEST_FILES_DIR}/io/vasp/inputs"
VASP_OUT_DIR = f"{TEST_FILES_DIR}/io/vasp/outputs"
VASP_FIXTURES = f"{TEST_FILES_DIR}/io/vasp/fixtures"
GAUSSIAN_LOG = f"{TEST_FILES_DIR}/analysis/quasirrho/quasirrho_gaufreq.log"

XE_VASPRUN = f"{TEST_DIR}/vasprun.xml.xe.gz"
LI_VASPRUN = f"{VASP_FIXTURES}/relaxation/vasprun.xml.gz"  # Li, GGA
LIH_VASPRUN = f"{VASP_FIXTURES}/scan_relaxation/vasprun.xml.gz"  # LiH, R2SCAN

OSZICAR_ENERGY = -526.63928


def _walk_tuple(path: Path) -> tuple[str, list[str], list[str]]:
    """Return the first os.walk tuple for path."""
    return next(os.walk(path))


@pytest.fixture
def relax_dir(tmp_path: Path) -> Path:
    """An aflow-style double relaxation: relax1 (Li) and relax2 (LiH)."""
    for sub, src in (("relax1", LI_VASPRUN), ("relax2", LIH_VASPRUN)):
        (tmp_path / sub).mkdir()
        shutil.copy(src, tmp_path / sub / "vasprun.xml.gz")
    return tmp_path


@pytest.fixture
def simple_vasp_dir(tmp_path: Path) -> Path:
    """A LiFePO4 run with INCAR, POTCAR, POSCAR, CONTCAR (+2% volume) and OSZICAR."""
    structure = Poscar.from_file(f"{VASP_IN_DIR}/POSCAR_LiFePO4", check_for_potcar=False).structure
    # Order species to match POTCAR_LiFePO4 (Li Fe P O)
    order = {"Li": 0, "Fe": 1, "P": 2, "O": 3}
    structure = structure.get_sorted_structure(key=lambda site: order[site.specie.symbol])
    Poscar(structure).write_file(tmp_path / "POSCAR")
    structure.scale_lattice(structure.volume * 1.02)
    Poscar(structure).write_file(tmp_path / "CONTCAR")
    Incar({"LDAU": True, "LDAUU": [0, 5.3, 0, 0], "ENCUT": 520}).write_file(tmp_path / "INCAR")
    shutil.copy(f"{VASP_IN_DIR}/POTCAR_LiFePO4.gz", tmp_path / "POTCAR.gz")
    shutil.copy(f"{VASP_OUT_DIR}/OSZICAR", tmp_path / "OSZICAR")
    return tmp_path


class TestVaspToComputedEntryDrone:
    def setup_method(self):
        self.drone = VaspToComputedEntryDrone(data=["efermi"])
        self.structure_drone = VaspToComputedEntryDrone(inc_structure=True)

    def test_str(self):
        assert str(self.drone) == "VaspToComputedEntryDrone"

    def test_get_valid_paths(self, tmp_path):
        assert self.drone.get_valid_paths(_walk_tuple(TEST_DIR)) == [TEST_DIR]

        # POSCAR + OSZICAR without vasprun.xml is also considered valid
        (tmp_path / "POSCAR").touch()
        assert self.drone.get_valid_paths(_walk_tuple(tmp_path)) == []
        (tmp_path / "OSZICAR").touch()
        assert self.drone.get_valid_paths(_walk_tuple(tmp_path)) == [str(tmp_path)]

    def test_get_valid_paths_relax_dirs(self, relax_dir):
        walked = {parent: (parent, subdirs, files) for parent, subdirs, files in os.walk(relax_dir)}
        # The parent of relax1/relax2 is the only valid path
        assert self.drone.get_valid_paths(walked[str(relax_dir)]) == [str(relax_dir)]
        for sub in ("relax1", "relax2"):
            assert self.drone.get_valid_paths(walked[str(relax_dir / sub)]) == []

    def test_get_valid_paths_empty_dir(self, tmp_path):
        assert self.drone.get_valid_paths(_walk_tuple(tmp_path)) == []

    def test_assimilate(self):
        """Test assimilate data from "vasprun.xml.xe.gz" file."""

        entry = self.drone.assimilate(f"{TEST_DIR}")

        for param in ("hubbards", "is_hubbard", "potcar_spec", "run_type"):
            assert param in entry.parameters
        assert entry.data["efermi"] == approx(-6.62148548)
        assert entry.reduced_formula == "Xe"
        assert entry.energy == approx(0.5559329)

        entry = self.structure_drone.assimilate(f"{TEST_DIR}")
        assert entry.reduced_formula == "Xe"
        assert entry.energy == approx(0.5559329)
        assert isinstance(entry, ComputedStructureEntry)
        assert entry.structure is not None

    def test_assimilate_parameters(self):
        drone = VaspToComputedEntryDrone(parameters=["incar"], data=["efermi", "final_energy"])
        entry = drone.assimilate(TEST_DIR)
        assert {"incar", "is_hubbard", "potcar_symbols", "run_type"} <= set(entry.parameters)
        assert entry.data["final_energy"] == approx(0.5559329)
        assert entry.data["efermi"] == approx(-6.62148548)

    def test_assimilate_relax_dirs(self, relax_dir):
        """Only relax2 is parsed for relax1/relax2 runs."""
        entry = self.structure_drone.assimilate(relax_dir)
        assert isinstance(entry, ComputedStructureEntry)
        assert entry.reduced_formula == "LiH"
        assert entry.energy == approx(-6.82078391)
        assert entry.parameters["run_type"] == "R2SCAN"

    def test_assimilate_multiple_vaspruns(self, tmp_path):
        shutil.copy(LI_VASPRUN, tmp_path / "vasprun.xml.relax1.gz")
        shutil.copy(LIH_VASPRUN, tmp_path / "vasprun.xml.relax2.gz")
        with pytest.warns(UserWarning, match="2 vasprun.xml.* found.*vasprun.xml.relax2.gz is being parsed"):
            entry = self.drone.assimilate(tmp_path)
        assert isinstance(entry, ComputedEntry)
        assert not isinstance(entry, ComputedStructureEntry)
        assert entry.reduced_formula == "LiH"
        assert entry.energy == approx(-6.82078391)

    def test_assimilate_bad_vasprun(self, tmp_path):
        with gzip.open(tmp_path / "vasprun.xml.gz", "wt", encoding="utf-8") as file:
            file.write("<modeling><not valid xml")
        assert self.drone.assimilate(tmp_path) is None

    def test_assimilate_no_vasprun(self, tmp_path):
        """A POSCAR+OSZICAR dir passes get_valid_paths but has no vasprun.xml."""
        (tmp_path / "POSCAR").touch()
        (tmp_path / "OSZICAR").touch()
        assert self.drone.get_valid_paths(_walk_tuple(tmp_path)) == [str(tmp_path)]
        assert self.drone.assimilate(tmp_path) is None

    def test_as_from_dict(self):
        drone = VaspToComputedEntryDrone(inc_structure=True, parameters=["incar"], data=["efermi"])
        dct = drone.as_dict()
        assert dct["@class"] == "VaspToComputedEntryDrone"
        assert dct["init_args"]["inc_structure"] is True
        assert "incar" in dct["init_args"]["parameters"]
        assert dct["init_args"]["data"] == ["efermi"]

        # Must survive a JSON round trip
        round_tripped = MontyDecoder().process_decoded(orjson.loads(drone.to_json()))
        assert isinstance(round_tripped, VaspToComputedEntryDrone)
        assert round_tripped.as_dict() == dct
        assert isinstance(VaspToComputedEntryDrone.from_dict(dct), VaspToComputedEntryDrone)


class TestSimpleVaspToComputedEntryDrone:
    def setup_method(self):
        self.drone = SimpleVaspToComputedEntryDrone()
        self.structure_drone = SimpleVaspToComputedEntryDrone(inc_structure=True)

    def test_str(self):
        assert str(self.drone) == "SimpleVaspToComputedEntryDrone"

    def test_get_valid_paths(self, simple_vasp_dir):
        assert self.drone.get_valid_paths(_walk_tuple(simple_vasp_dir)) == [str(simple_vasp_dir)]

    def test_assimilate_miss_core_file(self):
        with pytest.raises(ValueError, match="not all necessary files are present"):
            self.drone.assimilate(".")

    @pytest.mark.parametrize("inc_structure", [False, True])
    def test_assimilate(self, simple_vasp_dir, inc_structure):
        drone = SimpleVaspToComputedEntryDrone(inc_structure=inc_structure)
        entry = drone.assimilate(simple_vasp_dir)
        assert isinstance(entry, ComputedStructureEntry) is inc_structure
        assert entry.composition == Composition("Li4Fe4P4O16")
        assert entry.energy == approx(OSZICAR_ENERGY)
        assert entry.parameters["hubbards"] == {"Li": 0, "Fe": 5.3, "P": 0, "O": 0}
        assert entry.parameters["is_hubbard"] is True
        assert entry.parameters["run_type"] is None
        assert [spec["titel"] for spec in entry.parameters["potcar_spec"]] == [
            "PAW_PBE Li 17Jan2003",
            "PAW_PBE Fe 06Sep2000",
            "PAW_PBE P 17Jan2003",
            "PAW_PBE O 08Apr2002",
        ]
        assert entry.data["delta_volume"] == approx(0.02)
        assert entry.data["filename"] == simple_vasp_dir
        assert "phonon_frequencies" not in entry.data

    def test_assimilate_no_hubbard(self, simple_vasp_dir):
        Incar({"ENCUT": 520}).write_file(simple_vasp_dir / "INCAR")
        entry = self.drone.assimilate(simple_vasp_dir)
        assert entry.parameters["hubbards"] == {}
        assert entry.parameters["is_hubbard"] is False

    def test_assimilate_dynmat(self, simple_vasp_dir):
        shutil.copy(f"{VASP_OUT_DIR}/DYNMAT", simple_vasp_dir / "DYNMAT")
        entry = self.drone.assimilate(simple_vasp_dir)
        assert entry.data["phonon_frequencies"][:3] == approx([58.40051742, 58.37077149, 58.01005265])

    def test_assimilate_multiple_files(self, simple_vasp_dir):
        """With several CONTCAR/OSZICAR/POSCARs, the first POSCAR and last CONTCAR/OSZICAR are used."""
        shutil.copy(simple_vasp_dir / "POSCAR", simple_vasp_dir / "POSCAR.relax2")
        shutil.copy(simple_vasp_dir / "POSCAR", simple_vasp_dir / "CONTCAR.relax2")  # zero volume change
        shutil.copy(simple_vasp_dir / "OSZICAR", simple_vasp_dir / "OSZICAR.relax2")
        with pytest.warns(UserWarning, match="2 files found"):
            entry = self.drone.assimilate(simple_vasp_dir)
        assert entry.energy == approx(OSZICAR_ENERGY)
        assert entry.data["delta_volume"] == approx(0)

    def test_assimilate_relax_dirs(self, simple_vasp_dir, tmp_path):
        """INCAR/POTCAR/POSCAR come from relax1, CONTCAR/OSZICAR from relax2."""
        run_dir = tmp_path / "run"
        relax1, relax2 = run_dir / "relax1", run_dir / "relax2"
        relax1.mkdir(parents=True)
        relax2.mkdir()
        for filename in ("INCAR", "POTCAR.gz", "POSCAR"):
            shutil.copy(simple_vasp_dir / filename, relax1 / filename)
        for filename in ("CONTCAR", "OSZICAR"):
            shutil.copy(simple_vasp_dir / filename, relax2 / filename)
        assert self.drone.get_valid_paths(_walk_tuple(run_dir)) == [str(run_dir)]

        entry = self.structure_drone.assimilate(run_dir)
        assert isinstance(entry, ComputedStructureEntry)
        assert entry.composition.reduced_formula == "LiFePO4"
        assert entry.energy == approx(OSZICAR_ENERGY)
        assert entry.data["delta_volume"] == approx(0.02)

    def test_assimilate_bad_file(self, simple_vasp_dir):
        """Parse errors other than missing files are swallowed and None is returned."""
        (simple_vasp_dir / "OSZICAR").write_text("garbage\n", encoding="utf-8")
        assert self.drone.assimilate(simple_vasp_dir) is None

    def test_as_from_dict(self):
        dct = self.structure_drone.as_dict()
        assert dct["init_args"] == {"inc_structure": True}
        drone = SimpleVaspToComputedEntryDrone.from_dict(dct)
        assert isinstance(drone, SimpleVaspToComputedEntryDrone)
        assert drone.as_dict() == dct


class TestGaussianToComputedEntryDrone:
    def setup_method(self):
        self.drone = GaussianToComputedEntryDrone(data=["num_basis_func"])
        self.structure_drone = GaussianToComputedEntryDrone(inc_structure=True, parameters=["title"])

    def test_str(self):
        assert str(self.drone) == "GaussianToComputedEntryDrone"

    def test_get_valid_paths(self, tmp_path):
        for filename in ("mol.log", "mol.out", "notes.txt"):
            (tmp_path / filename).touch()
        assert self.drone.get_valid_paths(_walk_tuple(tmp_path)) == [str(tmp_path / "mol.log")]

        drone = GaussianToComputedEntryDrone(file_extensions=(".log", ".out"))
        assert sorted(drone.get_valid_paths(_walk_tuple(tmp_path))) == [
            str(tmp_path / "mol.log"),
            str(tmp_path / "mol.out"),
        ]

    def test_assimilate(self):
        entry = self.drone.assimilate(GAUSSIAN_LOG)
        assert type(entry) is ComputedEntry
        assert entry.composition == Composition("C14H28N2O3")
        assert entry.energy == approx(-885.155226974)
        gau_out = GaussianOutput(GAUSSIAN_LOG)
        for param in ("functional", "basis_set", "charge", "spin_multiplicity", "route_parameters"):
            assert entry.parameters[param] == getattr(gau_out, param)
        assert entry.parameters["basis_set"] == "6-311++g(2df,2p)"
        assert entry.parameters["charge"] == 2
        assert entry.parameters["spin_multiplicity"] == 1
        assert entry.data["stationary_type"] == "Minimum"
        assert entry.data["properly_terminated"] is True
        assert entry.data["num_basis_func"] > 0

        entry = self.structure_drone.assimilate(GAUSSIAN_LOG)
        assert isinstance(entry, ComputedStructureEntry)
        assert len(entry.structure) == 47
        assert entry.energy == approx(-885.155226974)
        assert "title" in entry.parameters

    def test_assimilate_bad_file(self, tmp_path):
        assert self.drone.assimilate(tmp_path / "does_not_exist.log") is None

    def test_as_from_dict(self):
        drone = GaussianToComputedEntryDrone(
            inc_structure=True, parameters=["title"], data=["num_basis_func"], file_extensions=(".out",)
        )
        dct = drone.as_dict()
        assert dct["init_args"]["inc_structure"] is True
        assert "title" in dct["init_args"]["parameters"]
        assert "num_basis_func" in dct["init_args"]["data"]
        assert dct["init_args"]["file_extensions"] == [".out"]

        round_tripped = MontyDecoder().process_decoded(orjson.loads(drone.to_json()))
        assert isinstance(round_tripped, GaussianToComputedEntryDrone)
        assert round_tripped.as_dict() == dct


class TestGetTransformationHistory:
    def test_no_file(self, tmp_path):
        assert _get_transformation_history(tmp_path) is None

    @pytest.mark.parametrize("filename", ["transformations.json", "transformations.json.gz"])
    def test_history(self, tmp_path, filename):
        history = [{"@class": "SubstitutionTransformation", "species_map": {"Li": "Na"}}]
        payload = orjson.dumps({"history": history})
        if filename.endswith(".gz"):
            payload = gzip.compress(payload)
        (tmp_path / filename).write_bytes(payload)
        assert _get_transformation_history(tmp_path) == history

    @pytest.mark.parametrize("content", [b"not json", b'{"no_history": []}'])
    def test_bad_file(self, tmp_path, content):
        (tmp_path / "transformations.json").write_bytes(content)
        assert _get_transformation_history(tmp_path) is None
