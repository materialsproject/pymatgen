from __future__ import annotations

import os
import subprocess
import zipfile
from unittest.mock import patch
from urllib.request import urlretrieve

import pytest

from pymatgen.command_line.chargemol_caller import ChargemolAnalysis
from pymatgen.core import Element, Structure
from pymatgen.util.testing import TEST_FILES_DIR

ddec_url = "https://sourceforge.net/projects/ddec/files/latest/download"
download_path = TEST_FILES_DIR / "command_line" / "chargemol" / "ddec-latest.gz"
extract_path = TEST_FILES_DIR / "command_line" / "chargemol"

if not os.path.exists(download_path):
    urlretrieve(ddec_url, download_path)

if not os.path.exists(extract_path / "chargemol_09_26_2017"):
    with zipfile.ZipFile(download_path, "r") as zip_ref:
        zip_ref.extractall(extract_path)

exe_path = extract_path / "chargemol_09_26_2017/chargemol_FORTRAN_09_26_2017/compiled_binaries/linux"

current_path = os.environ.get("PATH", "")
if str(exe_path) not in current_path:
    os.environ["PATH"] = str(exe_path) + os.pathsep + current_path

subprocess.call(["chmod", "+x", str(exe_path / "Chargemol_09_26_2017_linux_serial")])
subprocess.call(["chmod", "+x", str(exe_path / "Chargemol_09_26_2017_linux_parallel")])


class TestChargemolAnalysis:
    def test_parse_chargemol(self):
        test_dir = f"{TEST_FILES_DIR}/command_line/chargemol/spin_unpolarized"
        chg_mol = ChargemolAnalysis(path=test_dir, run_chargemol=False)
        assert chg_mol.ddec_charges == [0.8432, -0.8432]
        assert chg_mol.get_partial_charge(0) == 0.8432
        assert chg_mol.get_partial_charge(0, charge_type="cm5") == 0.420172
        assert chg_mol.get_charge_transfer(0) == -0.8432
        assert chg_mol.get_charge_transfer(0, charge_type="cm5") == -0.420172
        assert chg_mol.get_charge(0, nelect=1) == 1 - 0.8432
        assert chg_mol.get_charge(0, nelect=1, charge_type="cm5") == 1 - 0.420172
        assert chg_mol.dipoles == [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        assert chg_mol.ddec_spin_moments is None
        assert chg_mol.bond_order_sums == [0.53992, 0.901058]
        assert chg_mol.ddec_rsquared_moments == [8.261378, 34.237274]
        assert chg_mol.ddec_rcubed_moments == [14.496002, 88.169236]
        assert chg_mol.ddec_rfourth_moments == [37.648248, 277.371929]
        assert chg_mol.cm5_charges == [0.420172, -0.420172]
        assert chg_mol.summary["ddec"]["partial_charges"] == chg_mol.ddec_charges
        assert chg_mol.summary["ddec"]["dipoles"] == chg_mol.dipoles
        assert chg_mol.summary["ddec"]["bond_order_sums"] == chg_mol.bond_order_sums
        assert chg_mol.summary["ddec"]["rsquared_moments"] == chg_mol.ddec_rsquared_moments
        assert chg_mol.summary["ddec"]["rcubed_moments"] == chg_mol.ddec_rcubed_moments
        assert chg_mol.summary["ddec"]["rfourth_moments"] == chg_mol.ddec_rfourth_moments
        assert chg_mol.summary["cm5"]["partial_charges"] == chg_mol.cm5_charges
        assert chg_mol.summary["ddec"]["bond_order_dict"] == chg_mol.bond_order_dict
        assert chg_mol.summary["ddec"].get("spin_moments") is None
        assert chg_mol.natoms == [1, 1]
        assert isinstance(chg_mol.structure, Structure)
        assert len(chg_mol.structure) == 2
        assert chg_mol.structure.formula == "Na1 Cl1"
        assert len(chg_mol.bond_order_dict) == 2
        assert chg_mol.bond_order_dict[0]["bonded_to"][0] == {
            "spin_polarization": 0.0,
            "index": 1,
            "direction": (-1, -1, 0),
            "bond_order": 0.0882,
            "element": Element("Cl"),
        }
        assert chg_mol.bond_order_dict[1]["bonded_to"][-1]["direction"] == (-1, 0, 0)
        # check that partial charges are written to structure site properties
        charge_decorated_struct = chg_mol.get_property_decorated_structure()
        struct_partial_charge = charge_decorated_struct.site_properties["partial_charge_ddec6"]
        assert struct_partial_charge == chg_mol.ddec_charges

    def test_parse_chargemol2(self):
        test_dir = f"{TEST_FILES_DIR}/command_line/chargemol/spin_polarized"
        chg_mol = ChargemolAnalysis(path=test_dir, run_chargemol=False)
        assert chg_mol.ddec_spin_moments == [0.201595, 0.399203, 0.399203]
        assert chg_mol.summary["ddec"]["bond_order_dict"][0]["bonded_to"][0]["spin_polarization"] == 0.0490
        assert chg_mol.summary["ddec"]["spin_moments"] == chg_mol.ddec_spin_moments
        assert chg_mol.natoms is None
        assert chg_mol.structure is None

    def test_missing_exe_error(self):
        # monkeypatch CHARGEMOL_EXE to raise in ChargemolAnalysis.__init__
        patch.dict("os.environ", {"CHARGEMOL_EXE": "non_existent"})

        test_dir = f"{TEST_FILES_DIR}/command_line/chargemol/spin_unpolarized"
        with pytest.raises(
            OSError, match="ChargemolAnalysis requires the Chargemol executable to be in PATH. Please download"
        ):
            ChargemolAnalysis(path=test_dir, run_chargemol=True)

    def test_chargemol_custom_path(self):
        chg_mol = ChargemolAnalysis(
            path=TEST_FILES_DIR / "command_line" / "bader",
            atomic_densities_path=extract_path / "chargemol_09_26_2017" / "atomic_densities",
            run_chargemol=True,
        )
        print(chg_mol.ddec_charges)
