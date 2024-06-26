"""Tests the MD input set generator"""

from __future__ import annotations

from pathlib import Path

import pytest

from pymatgen.io.aims.sets.core import MDSetGenerator
from pymatgen.util.testing.aims import Si, compare_files

module_dir = Path(__file__).resolve().parents[1]
species_dir = module_dir / "species_directory"
ref_path = (module_dir / "aims_input_generator_ref").resolve()


def test_si_md(tmp_path):
    # default behaviour
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [2, 2, 2],
    }
    test_name = "md-si"
    # First do the exceptions
    with pytest.raises(ValueError, match="Ensemble something not valid"):
        MDSetGenerator(ensemble="something").get_input_set(Si)
    with pytest.raises(ValueError, match="Type parrinello is not valid for nve ensemble"):
        MDSetGenerator(ensemble="nve", ensemble_specs={"type": "parrinello"}).get_input_set(Si)
    with pytest.raises(ValueError, match="Velocities must be initialized"):
        MDSetGenerator(ensemble="nve", init_velocities=False).get_input_set(Si)
    with pytest.raises(ValueError, match="Temperature must be set"):
        MDSetGenerator(ensemble="nve").get_input_set(Si)
    with pytest.raises(ValueError, match="Temperature must be set"):
        MDSetGenerator(ensemble="nvt").get_input_set(Si)
    with pytest.raises(ValueError, match="parameter is not defined"):
        MDSetGenerator(ensemble="nve", ensemble_specs={"type": "damped"}, temp=300).get_input_set(Si)
    # then do the actual input set
    generator = MDSetGenerator(
        ensemble="nvt",
        ensemble_specs={"type": "parrinello", "parameter": 0.4},
        temp=300,
        time=10.0,
        time_step=0.002,
        user_params=parameters,
    )
    input_set = generator.get_input_set(Si)
    input_set.write_input(tmp_path / test_name)

    return compare_files(test_name, tmp_path, ref_path)
