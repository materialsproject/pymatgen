"""TODOs:
- also check type (floatwithunit vs float, str vs float)
- factor seems to cause some rounding issues
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from compare_json import compare_jsons

ROOT = Path(__file__).parent.resolve()

VENV_PYPI = ROOT / ".venv_pypi"
VENV_PR = ROOT / ".venv_pr"

OUTPUT_PYPI = ROOT / "element_properties_pypi.json"
OUTPUT_PR = ROOT / "element_properties_pr.json"

PR_URL = "https://github.com/DanielYang59/pymatgen.git"
PR_BRANCH = "directly-overwrite-json-in-core"

PYPI_VERSION = "2025.4.24"


def create_venv(path):
    subprocess.run([sys.executable, "-m", "venv", path], check=True)


def install_in_venv(venv_path, install_target):
    pip = venv_path / "bin" / "pip"
    subprocess.run([pip, "install", install_target], check=True)


def run_dump_in_venv(venv_path, output_path):
    python = venv_path / "bin" / "python"
    subprocess.run([python, "dump_all_properties.py"], check=True)
    shutil.move("element_properties.json", output_path)


if __name__ == "__main__":
    # print("🔧 Creating virtualenvs...")
    # create_venv(VENV_PYPI)
    # create_venv(VENV_PR)

    # print("⬇️ Installing PyPI pymatgen...")
    # install_in_venv(VENV_PYPI, f"pymatgen=={PYPI_VERSION}")

    # print("⬇️ Installing PR pymatgen...")
    # install_in_venv(VENV_PR, f"git+{PR_URL}@{PR_BRANCH}")

    print("📦 Dumping properties from PyPI...")
    run_dump_in_venv(VENV_PYPI, OUTPUT_PYPI)

    print("📦 Dumping properties from PR...")
    run_dump_in_venv(VENV_PR, OUTPUT_PR)

    print("🔍 Comparing results...")
    compare_jsons(OUTPUT_PYPI, OUTPUT_PR, log_path="diff.log")
