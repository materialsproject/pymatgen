"""
Test the import time of several important modules.
"""
# ruff: noqa: T201 (check for print statement)

from __future__ import annotations

import json
import os
import subprocess
import time
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from typing import Literal

# Toggle this to generate reference import times
GEN_REF_TIME = True

MODULES_TO_TEST = (
    "from pymatgen.core.bonds import CovalentBond",
    "from pymatgen.core.composition import Composition",
    "from pymatgen.core.interface import Interface",
    "from pymatgen.core.ion import Ion",
    "from pymatgen.core.lattice import Lattice",
    "from pymatgen.core.libxcfunc import LibxcFunc",
    "from pymatgen.core.molecular_orbitals import MolecularOrbitals",
    "from pymatgen.core.operations import SymmOp",
    "from pymatgen.core.periodic_table import Element",
    "from pymatgen.core.sites import Site",
    "from pymatgen.core.spectrum import Spectrum",
    "from pymatgen.core.structure import Structure",
    "from pymatgen.core.surface import Slab",
    "from pymatgen.core.tensors import Tensor",
    "from pymatgen.core.trajectory import Trajectory",
    "from pymatgen.core.units import Unit",
    "from pymatgen.core.xcfunc import XcFunc",
)

# Get runner OS and reference file
RUNNER_OS: Literal["linux", "windows", "macos"] = os.getenv("RUNNER_OS", "").lower()  # type: ignore[assignment]
assert RUNNER_OS in {"linux", "windows", "macos"}
REF_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"import_time_{RUNNER_OS}.json")


@pytest.mark.skipif(not GEN_REF_TIME, reason="Set GEN_REF_TIME to generate reference import time.")
def test_get_ref_import_time() -> None:
    """A dummy test that would always fail, used to generate copyable reference time."""
    import_times = {
        module_import_cmd: _measure_import_time_in_ms(module_import_cmd) for module_import_cmd in MODULES_TO_TEST
    }

    # Print a copyable JSON format for easy reference updating
    print("\nCopyable import time JSON:")
    print(json.dumps(import_times, indent=4))

    pytest.fail("Reference import times generated. Copy from output to update JSON file.")


@pytest.mark.skipif(GEN_REF_TIME, reason="Generating reference import time.")
def test_import_time(grace_percent: float = 0.20, hard_percent: float = 0.50) -> None:
    """Test the import time of core modules to avoid performance regression.

    Args:
        grace_percent (float): Maximum allowed percentage increase in import time
            before a warning is raised.
        hard_percent (float): Maximum allowed percentage increase in import time
            before the test fails.
    """
    try:
        with open(REF_FILE, encoding="utf-8") as file:
            ref_import_times = json.load(file)
    except FileNotFoundError:
        pytest.fail(f"Reference file {REF_FILE} not found. Please generate it.")

    for module_import_cmd, ref_time in ref_import_times.items():
        current_time = _measure_import_time_in_ms(module_import_cmd)

        # Calculate thresholds for grace and hard limits
        grace_threshold = ref_time * (1 + grace_percent)
        hard_threshold = ref_time * (1 + hard_percent)

        if current_time > grace_threshold:
            if current_time > hard_threshold:
                pytest.fail(f"{module_import_cmd} import too slow! {hard_threshold=:.2f} ms")
            else:
                pytest.warns(
                    UserWarning,
                    f"{module_import_cmd} import slightly slower than reference: {grace_threshold=:.2f} ms",
                )


def _measure_import_time_in_ms(module_import_cmd: str, count: int = 3) -> float:
    """Measure import time of a module in milliseconds across several runs.

    Args:
        module_import_cmd (str): The module import command.
        count (int): Number of runs to average.

    Returns:
        float: Import time in milliseconds.
    """
    total_time = 0.0

    for _ in range(count):
        start_time = time.perf_counter()
        subprocess.run(["python", "-c", f"{module_import_cmd}"], check=True)
        total_time += time.perf_counter() - start_time

    return (total_time / count) * 1000
