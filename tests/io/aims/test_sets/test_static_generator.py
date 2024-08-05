from __future__ import annotations

from pathlib import Path

import numpy as np

from pymatgen.io.aims.sets.core import StaticSetGenerator
from pymatgen.util.testing.aims import O2, Si, comp_system

module_dir = Path(__file__).resolve().parents[1]
ref_path = (module_dir / "aims_input_generator_ref").resolve()


def test_static_si(tmp_path):
    parameters = {
        "species_dir": "light",
        "k_grid": [2, 2, 2],
    }
    comp_system(Si, parameters, "static-si", tmp_path, ref_path, StaticSetGenerator)


def test_static_si_no_kgrid(tmp_path):
    parameters = {"species_dir": "light"}
    Si_supercell = Si.make_supercell(np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]]), in_place=False)
    comp_system(Si_supercell, parameters, "static-no-kgrid-si", tmp_path, ref_path, StaticSetGenerator)


def test_static_o2(tmp_path):
    parameters = {"species_dir": "light"}
    comp_system(O2, parameters, "static-o2", tmp_path, ref_path, StaticSetGenerator)
