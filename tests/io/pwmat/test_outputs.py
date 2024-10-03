from __future__ import annotations

from pymatgen.io.pwmat.outputs import DosSpin, Movement, OutFermi, Report
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/io/pwmat"


class TestMovement(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_DIR}/MOVEMENT.lzma"
        movement = Movement(filepath)
        assert movement.n_ionic_steps == len(movement.chunk_sizes) == 1
        assert movement.n_ionic_steps == len(movement.chunk_starts)
        assert movement.n_ionic_steps == len(movement.atom_configs)
        assert "atom_config" in movement.ionic_steps[0]
        assert "e_tot" in movement.ionic_steps[0]
        assert "atom_forces" in movement.ionic_steps[0]
        assert "virial" in movement.ionic_steps[0]
        assert "atom_energies" in movement.ionic_steps[0]
        assert movement.e_tots == -357677.2281
        assert movement.atom_configs[0] == movement.ionic_steps[0]["atom_config"]
        assert list(movement.e_atoms) == []
        assert movement.atom_forces.shape == (1, 72, 3)
        assert movement.virials.shape == (1, 3, 3)
        assert movement.ionic_steps[0]["e_tot"] == -357677.2281


class TestOutFermi(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_DIR}/OUT.FERMI.lzma"
        out_fermi = OutFermi(filepath)
        assert out_fermi.e_fermi == -2.359


class TestReport(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_DIR}/REPORT"
        report = Report(filepath)
        assert report.spin == 1
        assert report.n_kpoints == 1
        assert report.n_bands == 10
        assert report.eigenvalues.shape == (1, 1, 10)
        assert report.kpoints.shape == (1, 3)
        assert report.kpoints_weight.shape == (1,)
        assert report.hsps == {}


class TestDosSpin(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_DIR}/DOS.spinup_projected"
        dos_spin = DosSpin(filepath)
        assert dos_spin.dos.shape == (20, 21)
        assert dos_spin.labels == [
            "Energy",
            "Total",
            "Cr-3S",
            "Cr-3Pz",
            "Cr-3Px",
            "Cr-3Py",
            "Cr-4S",
            "Cr-3Dz2",
            "Cr-3Dxz",
            "Cr-3Dyz",
            "Cr-3D(x^2-y^2)",
            "Cr-3Dxy",
            "I-4Dz2",
            "I-4Dxz",
            "I-4Dyz",
            "I-4D(x^2-y^2)",
            "I-4Dxy",
            "I-5S",
            "I-5Pz",
            "I-5Px",
            "I-5Py",
        ]
