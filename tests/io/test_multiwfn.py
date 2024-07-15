from __future__ import annotations

import numpy as np
import pytest
from monty.json import MontyDecoder, jsanitize

from pymatgen.core.structure import Molecule
from pymatgen.io.multiwfn import(
    QTAIM_CONDITIONALS,
    extract_info_from_cp_text,
    parse_cp,
    get_qtaim_descs,
    separate_cps_by_type,
    match_atom_cp,
    map_atoms_cps,
    add_atoms,
    process_multiwfn_qtaim,
)
from pymatgen.util.testing import TEST_FILES_DIR


base_dir = TEST_FILES_DIR / "io" / "multiwfn"


def test_parse_single_cp():
    pass


def test_parse_cps():
    pass


def test_separate_by_type():
    pass


def test_atom_matching():
    pass


def test_add_atoms():
    pass


def test_process_multiwfn_qtaim():
    pass