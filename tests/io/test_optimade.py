from __future__ import annotations

import numpy as np

from pymatgen.core import Structure
from pymatgen.io.optimade import OptimadeStructureAdapter
from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR

STRUCTURE = Structure.from_file(f"{VASP_IN_DIR}/POSCAR")
XYZ_STRUCTURE = f"{TEST_FILES_DIR}/io/xyz/acetylene.xyz"


def test_get_optimade_structure_roundtrip():
    optimade_structure = OptimadeStructureAdapter.get_optimade_structure(STRUCTURE)

    assert optimade_structure["attributes"]["nsites"] == len(STRUCTURE)
    assert optimade_structure["attributes"]["elements"] == ["Fe", "O", "P"]
    assert optimade_structure["attributes"]["nelements"] == 3
    assert optimade_structure["attributes"]["chemical_formula_reduced"] == "FeO4P"
    assert optimade_structure["attributes"]["species_at_sites"] == 4 * ["Fe"] + 4 * ["P"] + 16 * ["O"]
    np.testing.assert_array_almost_equal(
        np.abs(optimade_structure["attributes"]["lattice_vectors"]), np.abs(STRUCTURE.lattice.matrix)
    )

    # Set an OPTIMADE ID and some custom properties and ensure they are preserved in the properties
    test_id = "test_id"
    optimade_structure["id"] = test_id
    custom_properties = {"_custom_field": "test_custom_field", "_custom_band_gap": 2.2}
    optimade_structure["attributes"].update(custom_properties)

    roundtrip_structure = OptimadeStructureAdapter.get_structure(optimade_structure)
    assert roundtrip_structure.properties["optimade_id"] == test_id
    assert roundtrip_structure.properties["optimade_attributes"] == custom_properties

    # Delete the properties before the check for equality
    roundtrip_structure.properties = {}
    assert roundtrip_structure == STRUCTURE
