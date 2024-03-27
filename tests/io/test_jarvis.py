from __future__ import annotations

import pytest

from pymatgen.core import Structure
from pymatgen.io.jarvis import Atoms, JarvisAtomsAdaptor
from pymatgen.util.testing import VASP_IN_DIR


@pytest.mark.skipif(Atoms is None, reason="JARVIS-tools not loaded.")
class TestJarvisAtomsAdaptor:
    def test_get_atoms_from_structure(self):
        struct = Structure.from_file(f"{VASP_IN_DIR}/POSCAR")
        atoms = JarvisAtomsAdaptor.get_atoms(struct)
        jarvis_composition = atoms.composition.reduced_formula
        assert jarvis_composition == struct.reduced_formula
        assert atoms.lattice_mat is not None
        assert atoms.lattice_mat.any()

    def test_get_structure(self):
        struct = Structure.from_file(f"{VASP_IN_DIR}/POSCAR")
        atoms = JarvisAtomsAdaptor.get_atoms(struct)
        assert len(atoms.frac_coords) == len(struct) == 24
        round_trip = JarvisAtomsAdaptor.get_structure(atoms)
        assert round_trip.formula == struct.formula == "Fe4 P4 O16"
