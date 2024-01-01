from __future__ import annotations

import os
import unittest
import zipfile
from typing import TYPE_CHECKING

import pytest

from pymatgen.command_line.chargemol_caller import ChargemolAnalysis
from pymatgen.core import Element
from pymatgen.util.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from pathlib import Path


class TestChargemolAnalysis(unittest.TestCase):
    def test_parse_chargemol(self):
        test_dir = f"{TEST_FILES_DIR}/chargemol/spin_unpolarized"
        ca = ChargemolAnalysis(path=test_dir, run_chargemol=False)
        assert ca.ddec_charges == [0.8432, -0.8432]
        assert ca.get_partial_charge(0) == 0.8432
        assert ca.get_partial_charge(0, charge_type="cm5") == 0.420172
        assert ca.get_charge_transfer(0) == -0.8432
        assert ca.get_charge_transfer(0, charge_type="cm5") == -0.420172
        assert ca.get_charge(0, nelect=1) == 1 - 0.8432
        assert ca.get_charge(0, nelect=1, charge_type="cm5") == 1 - 0.420172
        assert ca.dipoles == [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        assert ca.ddec_spin_moments is None
        assert ca.bond_order_sums == [0.53992, 0.901058]
        assert ca.ddec_rsquared_moments == [8.261378, 34.237274]
        assert ca.ddec_rcubed_moments == [14.496002, 88.169236]
        assert ca.ddec_rfourth_moments == [37.648248, 277.371929]
        assert ca.cm5_charges == [0.420172, -0.420172]
        summary = ca.summary["ddec"]
        assert summary["partial_charges"] == ca.ddec_charges
        assert summary["dipoles"] == ca.dipoles
        assert summary["bond_order_sums"] == ca.bond_order_sums
        assert summary["rsquared_moments"] == ca.ddec_rsquared_moments
        assert summary["rcubed_moments"] == ca.ddec_rcubed_moments
        assert summary["rfourth_moments"] == ca.ddec_rfourth_moments
        assert ca.summary["cm5"]["partial_charges"] == ca.cm5_charges
        assert summary["bond_order_dict"] == ca.bond_order_dict
        assert summary.get("spin_moments") is None
        assert ca.natoms == [1, 1]
        assert ca.structure is not None
        assert len(ca.bond_order_dict) == 2
        assert ca.bond_order_dict[0]["bonded_to"][0] == {
            "spin_polarization": 0.0,
            "index": 1,
            "direction": (-1, -1, 0),
            "bond_order": 0.0882,
            "element": Element("Cl"),
        }
        assert ca.bond_order_dict[1]["bonded_to"][-1]["direction"] == (-1, 0, 0)
        assert ca.get_property_decorated_structure().site_properties["partial_charge_ddec6"] == ca.ddec_charges

    def test_parse_chargemol2(self):
        test_dir = f"{TEST_FILES_DIR}/chargemol/spin_polarized"
        ca = ChargemolAnalysis(path=test_dir, run_chargemol=False)
        assert ca.ddec_spin_moments == [0.201595, 0.399203, 0.399203]
        assert ca.summary["ddec"]["bond_order_dict"][0]["bonded_to"][0]["spin_polarization"] == 0.0490
        assert ca.summary["ddec"]["spin_moments"] == ca.ddec_spin_moments
        assert ca.natoms is None
        assert ca.structure is None


@pytest.fixture()
def mock_xyz(tmp_path: Path):
    test_dir = tmp_path / "chargemol" / "spin_unpolarized"
    test_dir.mkdir(parents=True)
    xyz_file_content = """2

C       0.00000       0.00000       0.00000       0.8432
O       1.20000       0.00000       0.00000      -0.8432"""

    (test_dir / "DDEC6_even_tempered_net_atomic_charges.xyz").write_text(xyz_file_content)
    return test_dir


def test_parse_chargemol(monkeypatch, mock_xyz):
    monkeypatch.setattr(
        ChargemolAnalysis,
        "_download_and_unzip_atomic_densities",
        lambda self, version, verbose: fake_download(self),
    )

    ca = ChargemolAnalysis(path=str(mock_xyz), run_chargemol=False)
    ca._atomic_densities_path = os.path.expanduser("~/.cache/pymatgen/ddec/fake_file")

    assert ca.ddec_charges == [0.8432, -0.8432]


def fake_download(self, version: str = "latest", verbose: bool = True) -> None:
    extraction_path = "~/.cache/pymatgen/ddec"
    os.makedirs(os.path.expanduser(extraction_path), exist_ok=True)

    # Create a fake zipfile with a folder named "atomic_densities"
    fake_zip_path = os.path.expanduser(f"{extraction_path}/fake_file.zip")
    fake_folder_path = os.path.expanduser(f"{extraction_path}/atomic_densities")

    # Create a fake file inside the "atomic_densities" folder
    fake_file_path = os.path.join(fake_folder_path, "fake_file.txt")
    os.makedirs(fake_folder_path, exist_ok=True)
    with open(fake_file_path, "w") as f:
        f.write("fake content")

    # Create a zipfile containing the "atomic_densities" folder
    with zipfile.ZipFile(fake_zip_path, "w") as fake_zip:
        fake_zip.write(fake_file_path, arcname="atomic_densities/fake_file.txt")

    self._atomic_densities_path = os.path.expanduser(extraction_path)


def test_fake_download_and_modify_path(monkeypatch):
    # mock os.path.exists to simulate atomic densities are not found
    # monkeypatch.setattr(os.path, "exists", lambda *args, **kwargs: False)

    # mock subprocess.Popen to simulate a successful Chargemol execution
    monkeypatch.setattr("subprocess.Popen", lambda *args, **kwargs: type("", (), {"returncode": 0})())

    # monkeypatch urllib.request.urlretrieve to fake download
    import urllib.request

    monkeypatch.setattr(urllib.request, "urlretrieve", lambda url, path: None)

    # monkeypatch the _download_and_unzip_atomic_densities method
    monkeypatch.setattr(
        ChargemolAnalysis,
        "_download_and_unzip_atomic_densities",
        lambda self, version, verbose: fake_download(self),
    )

    # Create an instance of ChargemolAnalysis
    test_dir = f"{TEST_FILES_DIR}/chargemol/spin_unpolarized"
    ca = ChargemolAnalysis(
        path=test_dir,
        # atomic_densities_path=os.path.expanduser("~/.cache/pymatgen/ddec"),
        run_chargemol=False,
    )

    # Your assertions
    assert ca._atomic_densities_path == os.getcwd()
    # assert ca._atomic_densities_path == os.path.expanduser("~/.cache/pymatgen/ddec")
    # TODO: Add more tests to ensure Chargemol was run correctly
