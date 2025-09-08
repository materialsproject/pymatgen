from __future__ import annotations

from pymatgen.cli import pmg
from pymatgen.util.testing import VASP_IN_DIR


def test_pmg_diff():
    pmg.main(["diff", "--incar", f"{VASP_IN_DIR}/INCAR", f"{VASP_IN_DIR}/INCAR_2"])


def test_pmg_query():
    # TODO: add test
    pass


def test_pmg_view():
    # TODO: add test
    pass
