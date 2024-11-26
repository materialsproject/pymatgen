"""Check the content of source code and test build source distribution.

Binary distribution (wheel) is checked in test workflow.
"""

from __future__ import annotations

import subprocess
import tarfile
from glob import glob
from pathlib import Path
from typing import TYPE_CHECKING

from monty.tempfile import ScratchDir

if TYPE_CHECKING:
    from pymatgen.util.typing import PathLike

PROJECT_ROOT = Path(__file__).parent.parent
NAMESPACE_PKGS = {"analysis", "ext", "io"}  # TODO: double check


def test_source_code():
    """Check the source code in the working directory."""
    src_txt_path = PROJECT_ROOT / "src/pymatgen.egg-info/SOURCES.txt"

    _check_src_txt_is_complete(PROJECT_ROOT, src_txt_path)

    # Check existence of `py.typed` file in each sub-package
    _check_py_typed_files(PROJECT_ROOT / "src/pymatgen")


def test_build_source_distribution():
    """Test build the source distribution (sdist), also check `py.typed` files."""
    with ScratchDir("."):
        # Build the source distribution
        subprocess.run(["python", "-m", "pip", "install", "--upgrade", "build"], check=True)
        subprocess.run(["python", "-m", "build", "--sdist", PROJECT_ROOT, "--outdir", ".", "-C--quiet"], check=True)

        # Decompress sdist
        sdist_file = next(Path(".").glob("*.tar.gz"))
        sdist_dir = Path(sdist_file.name.removesuffix(".tar.gz"))
        with tarfile.open(sdist_file, "r:gz") as tar:
            # TODO: remove attr check after only 3.12+
            if hasattr(tarfile, "data_filter"):
                tar.extractall("", filter="data")
            else:
                tar.extractall("")  # noqa: S202

        # Check existence of `py.typed` file in each sub-package
        _check_py_typed_files(sdist_dir / "src/pymatgen")


def _check_src_txt_is_complete(project_root: PathLike, src_txt_path: PathLike) -> None:
    """Check that all source code and data files are listed in given SOURCES.txt.

    Args:
        project_root (PathLike): Path to the directory containing "src/pymatgen/".
        src_txt_path (PathLike): Path to the "SOURCES.txt" file.
    """
    project_root = Path(project_root)
    src_txt_path = Path(src_txt_path)

    assert project_root.is_dir(), f"{project_root} is not a directory"
    assert src_txt_path.is_file(), f"{src_txt_path} doesn't exist"

    sources = src_txt_path.read_text(encoding="utf-8")

    # Check that all files listed in "SOURCES.txt" exist
    for src_file in sources.splitlines():
        assert (project_root / src_file).is_file(), f"{src_file!r} does not exist!"

    # Check that all files in src/pymatgen/ are listed in SOURCES.txt
    for ext in ("py", "json", "json.*", "yaml", "csv"):
        for filepath in glob(f"{project_root}/src/pymatgen/**/*.{ext}", recursive=True):
            unix_path = Path(filepath).relative_to(project_root).as_posix()

            if unix_path not in sources:
                raise ValueError(f"{unix_path} not found in {src_txt_path}, check package data config")


def _check_py_typed_files(pkg_root: PathLike) -> None:
    """
    Ensure all sub-packages contain a `py.typed` file.

    Args:
        pkg_root (PathLike): Path to the namespace package's root directory.

    Raises:
        FileNotFoundError: If any sub-package is missing the `py.typed` file.
    """
    if not (pkg_root := Path(pkg_root)).is_dir():
        raise NotADirectoryError(f"Provided path {pkg_root} is not a valid directory.")

    # Iterate through all directories under the namespace package
    for sub_pkg in pkg_root.glob("*/"):
        if sub_pkg.name.startswith(".") or sub_pkg.name.startswith("__"):
            continue

        # Check for __init__.py to ensure it's not a namespace package
        sub_pkg_path = Path(sub_pkg)
        if (sub_pkg_path / "__init__.py").exists():
            if not (sub_pkg_path / "py.typed").exists():
                raise FileNotFoundError(f"Missing py.typed in sub-package: {sub_pkg_path}")

        elif sub_pkg_path.name not in NAMESPACE_PKGS:
            raise ValueError(f"Unexpected namespace package {sub_pkg_path.name}")
