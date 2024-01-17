from __future__ import annotations

from pymatgen.io.pwmat.outputs import DosSpin, Movement, OutFermi, Report
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestMovement(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/MOVEMENT.lzma"
        movement = Movement(filepath)
        assert movement.n_ionic_steps == len(movement.chunk_sizes)
        assert movement.n_ionic_steps == len(movement.chunk_starts)
        assert movement.n_ionic_steps == len(movement.atom_configs)
        assert movement.atom_configs[0].structure.num_sites == movement.atom_configs[0].structure.num_sites
        assert "atom_config" in movement.ionic_steps[0]
        assert "etot" in movement.ionic_steps[0]
        assert "fatoms" in movement.ionic_steps[0]
        # assert "virial" in movement.ionic_steps
        # assert "eatoms" in movement.ionic_steps
        assert hasattr(movement, "atom_configs")
        assert hasattr(movement, "etots")
        assert hasattr(movement, "fatoms")
        # assert hasattr(movement, "eatoms")
        # assert hasattr(movement, "virials")


class TestOutFermi(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/OUT.FERMI.lzma"
        out_fermi = OutFermi(filepath)
        assert isinstance(out_fermi.efermi, float)


class TestReport(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/REPORT"
        report = Report(filepath)
        assert hasattr(report, "spin")
        assert hasattr(report, "n_kpoints")
        assert hasattr(report, "n_bands")
        assert hasattr(report, "eigenvalues")
        assert hasattr(report, "kpoints")
        assert hasattr(report, "kpoints_weight")
        assert hasattr(report, "hsps")


class TestDosSpin(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/DOS.spinup_projected.lzma"
        dos_spin = DosSpin(filepath)
        assert dos_spin.dos.shape == (200, 21)
        assert dos_spin.labels == []
