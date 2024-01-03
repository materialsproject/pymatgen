from __future__ import annotations

from pymatgen.io.pwmat.outputs import DosSpin, Movement, OutFermi, Report
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestMovement(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/MOVEMENT"
        movement = Movement(filepath)
        assert movement.n_ionic_steps == len(movement.chunk_sizes)
        assert movement.n_ionic_steps == len(movement.chunk_starts)
        assert movement.n_ionic_steps == len(movement.atom_configs)
        assert movement.atom_configs[0].structure.num_sites == movement.atom_configs[0].structure.num_sites
        assert hasattr(movement.n_ionic_steps[0], "atom_configs")
        assert hasattr(movement.n_ionic_steps[0], "etot")
        assert hasattr(movement.n_ionic_steps[0], "fatoms")
        # assert hasattr(movement.nionic_steps[0], "virial")
        # assert hasattr(movement.nionic_steps[0], "eatoms")
        assert hasattr(movement, "atom_configs")
        assert hasattr(movement, "etots")
        assert hasattr(movement, "fatoms")
        # assert hasattr(movement, "eatoms")
        # assert hasattr(movement, "virials")


class TestOutFermi(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/OUT.FERMI"
        out_fermi = OutFermi(filepath)
        assert isinstance(out_fermi.efermi, float)


class TestReport(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/REPORT"
        report = Report(filepath)
        assert hasattr(report, "spin")
        assert hasattr(report, "nkpoints")
        assert hasattr(report, "nbands")
        assert hasattr(report, "eigenvalues")
        assert hasattr(report, "kpoints")
        assert hasattr(report, "kpoints_weight")
        assert hasattr(report, "hsps")


class TestDosspin(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/DOS.spinup_projected"
        dosspin = DosSpin(filepath)
        assert hasattr(dosspin, "dos")
        # assert type(dosspin.dos) == np.ndarray
        assert hasattr(dosspin, "label")
