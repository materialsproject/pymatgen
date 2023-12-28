from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.io.pwmat.outputs import(
    Movement,
    OutFermi,
    Report,
    Dosspin
)
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

if TYPE_CHECKING:
    from pathlib import Path
    
    
class TestMovement(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/MOVEMENT"
        movement = Movement(filepath)
        assert movement.nionic_steps == len(movement.chunksizes)
        assert movement.nionic_steps == len(movement.chunkstarts)
        assert movement.nionic_steps == len(movement.atom_configs)
        assert movement.atom_configs[0].structure.num_sites == movement.atom_configs[0].structure.num_sites
        assert hasattr(movement.nionic_steps[0], "atom_configs")
        assert hasattr(movement.nionic_steps[0], "etot")
        assert hasattr(movement.nionic_steps[0], "fatoms")
        #assert hasattr(movement.nionic_steps[0], "virial")
        #assert hasattr(movement.nionic_steps[0], "eatoms")
        assert hasattr(movement, "atom_configs")
        assert hasattr(movement, "etots")
        assert hasattr(movement, "fatoms")
        #assert hasattr(movement, "eatoms")
        #assert hasattr(movement, "virials")
    

class TestOutFermi(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/OUT.FERMI"
        out_fermi = OutFermi(filepath)
        assert type(out_fermi.efermi) == float


class TestReport(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/REPORT"
        report = Report(filepath)
        assert hasattr(report, "spin")
        assert hasattr(report, "num_kpts")
        assert hasattr(report, "num_bands")
        assert hasattr(report, "eigenvalues")
        assert hasattr(report, "kpts")
        assert hasattr(report, "kpts_weight")
        assert hasattr(report, "hsps")


class TestDosspin(PymatgenTest):
    def test_init_and_properties(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/DOS.spinup_projected"
        dosspin = Dosspin(filepath)
        assert hasattr(dosspin, "dos")
        #assert type(dosspin.dos) == np.ndarray
        assert hasattr(dosspin, "label")
        