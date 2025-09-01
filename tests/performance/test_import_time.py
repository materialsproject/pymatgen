"""
Test the import time of several important modules.

NOTE:
    - Toggle the "GEN_REF_TIME" to generate reference import time.

Last update: 2024-10-26 (PR 4128)

Runner specs:
Linux:   4 CPU
Windows: 4 CPU
macOS:   3 CPU (M1)
"""

from __future__ import annotations

import os
import subprocess
import time
import warnings
from typing import TYPE_CHECKING

import orjson
import pytest

from pymatgen.util.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from typing import Literal

if not os.getenv("CI"):
    pytest.skip("ref time only comparable in CI runner", allow_module_level=True)

# NOTE: Toggle this to generate reference import time
GEN_REF_TIME: bool = False

MODULES_TO_TEST: tuple[str, ...] = (
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

# Skip test on macOS due to inconsistent import time results
if RUNNER_OS == "macos":
    pytest.skip("Import time measurements are unstable on macOS", allow_module_level=True)

REF_FILE: str = f"{TEST_FILES_DIR}/performance/import_time_{RUNNER_OS}.json"


@pytest.mark.skipif(not GEN_REF_TIME, reason="Set GEN_REF_TIME to generate reference import time.")
def test_get_ref_import_time() -> None:
    """A dummy test that would always fail, used to generate copyable reference time."""
    import_times: dict[str, float] = {
        module_import_cmd: _measure_import_time_in_ms(module_import_cmd) for module_import_cmd in MODULES_TO_TEST
    }

    # Print a copyable JSON format for easy reference updating
    print("\nCopyable import time JSON:")
    print(orjson.dumps(import_times, option=orjson.OPT_INDENT_2).decode())

    pytest.fail("Reference import times generated. Copy from output to update JSON file.")


@pytest.mark.skipif(GEN_REF_TIME, reason="Generating reference import time.")
@pytest.mark.parametrize(
    ("grace_percent", "hard_percent"),
    [(0.5, 1.0)],
)
def test_import_time(grace_percent: float, hard_percent: float) -> None:
    """Test the import time of core modules to avoid performance regression.

    Args:
        grace_percent (float): Maximum allowed percentage increase in import time
            before a warning is raised.
        hard_percent (float): Maximum allowed percentage increase in import time
            before the test fails.
    """

    with open(REF_FILE, "rb") as file:
        ref_import_times: dict[str, float] = orjson.loads(file.read())

    for module_import_cmd, ref_time in ref_import_times.items():
        current_time: float = _measure_import_time_in_ms(module_import_cmd)

        # Calculate thresholds for grace and hard limits
        grace_threshold = ref_time * (1 + grace_percent)
        hard_threshold = ref_time * (1 + hard_percent)

        if current_time > grace_threshold:
            if current_time > hard_threshold:
                pytest.fail(f"{module_import_cmd} import too slow at {current_time:.2f} ms! {hard_threshold=:.2f} ms")
            else:
                warnings.warn(
                    f"{module_import_cmd} import slightly slower than reference: {grace_threshold=:.2f} ms",
                    stacklevel=2,
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
        start_time = time.perf_counter_ns()
        subprocess.run(["python", "-c", f"{module_import_cmd}"], check=True)
        total_time += time.perf_counter_ns() - start_time

    return (total_time / count) / 1e6
