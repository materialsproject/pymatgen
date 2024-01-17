from __future__ import annotations

import gzip
import json
from glob import glob
from pathlib import Path

from pymatgen.core import Lattice, Molecule, Structure
from pymatgen.io.aims.sets.core import StaticSetGenerator


def compare_files(test_name, work_dir, ref_dir):
    for file in glob(f"{work_dir / test_name}/*in"):
        with open(file) as test_file:
            test_lines = [line.strip() for line in test_file.readlines()[4:] if len(line.strip()) > 0]

        with gzip.open(f"{ref_dir / test_name / Path(file).name}.gz", "rt") as ref_file:
            ref_lines = [line.strip() for line in ref_file.readlines()[4:] if len(line.strip()) > 0]

        assert test_lines == ref_lines

    with open(f"{ref_dir / test_name}/parameters.json") as ref_file:
        ref = json.load(ref_file)
    ref.pop("species_dir", None)

    with open(f"{work_dir / test_name}/parameters.json") as check_file:
        check = json.load(check_file)
    check.pop("species_dir", None)

    assert ref == check


def comp_system(atoms, user_params, test_name, work_path, ref_path):
    generator = StaticSetGenerator(user_parameters=user_params)
    input_set = generator.get_input_set(atoms)
    input_set.write_input(work_path / test_name)
    compare_files(test_name, work_path, ref_path)


O2 = Molecule(species=["O", "O"], coords=[[0, 0, 0.622978], [0, 0, -0.622978]])

Si = Structure(
    lattice=Lattice([[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]),
    species=["Si", "Si"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)
species_dir = Path(__file__).resolve().parents[1] / "species_directory"

module_dir = Path(__file__).resolve().parents[1]
ref_path = (module_dir / "aims_input_generator_ref").resolve()


def test_static_si(tmp_path):
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [2, 2, 2],
    }
    comp_system(Si, parameters, "static-si", tmp_path, ref_path)


def test_static_si_no_kgrid(tmp_path):
    parameters = {"species_dir": str(species_dir / "light")}
    comp_system(Si, parameters, "static-no-kgrid-si", tmp_path, ref_path)


def test_static_default_species_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("AIMS_SPECIES_DIR", str(species_dir / "light"))
    parameters = {"k_grid": [2, 2, 2]}

    comp_system(Si, parameters, "static-default-sd-si", tmp_path, ref_path)


def test_static_o2(tmp_path):
    parameters = {"species_dir": str(species_dir / "light")}
    comp_system(O2, parameters, "static-o2", tmp_path, ref_path)


def test_static_default_species_dir_o2(tmp_path, monkeypatch):
    monkeypatch.setenv("AIMS_SPECIES_DIR", str(species_dir / "light"))
    parameters = {"k_grid": [2, 2, 2]}

    comp_system(O2, parameters, "static-default-sd-o2", tmp_path, ref_path)
