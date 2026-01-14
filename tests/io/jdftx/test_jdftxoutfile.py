from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from pymatgen.io.jdftx._output_utils import read_outfile_slices
from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice
from pymatgen.io.jdftx.outputs import JDFTXOutfile, _jof_atr_from_last_slice

from .outputs_test_utils import (
    etot_etype_outfile_known_simple,
    etot_etype_outfile_path,
    example_ionmin_outfile_known,
    example_ionmin_outfile_known_simple,
    example_ionmin_outfile_path,
    example_latmin_outfile_known,
    example_latmin_outfile_known_simple,
    example_latmin_outfile_path,
    example_sp_outfile_known,
    example_sp_outfile_known_simple,
    example_sp_outfile_path,
    example_vib_modes_known,
    example_vib_nrg_components,
    example_vib_outfile_path,
    jdftxoutfile_fromfile_matches_known,
    jdftxoutfile_fromfile_matches_known_simple,
    noeigstats_outfile_known_simple,
    noeigstats_outfile_path,
    partial_lattice_init_outfile_known_lattice,
    partial_lattice_init_outfile_path,
    problem2_outfile_known_simple,
    problem2_outfile_path,
)
from .shared_test_utils import assert_same_value

if TYPE_CHECKING:
    from pathlib import Path


@pytest.mark.parametrize(
    ("filename", "known"),
    [
        (example_sp_outfile_path, example_sp_outfile_known),
        (example_latmin_outfile_path, example_latmin_outfile_known),
        (example_ionmin_outfile_path, example_ionmin_outfile_known),
    ],
)
def test_JDFTXOutfile_fromfile(filename: Path, known: dict):
    jdftxoutfile_fromfile_matches_known(filename, known)


@pytest.mark.parametrize(
    ("filename", "known"),
    [
        (example_sp_outfile_path, example_sp_outfile_known_simple),
        (example_latmin_outfile_path, example_latmin_outfile_known_simple),
        (example_ionmin_outfile_path, example_ionmin_outfile_known_simple),
        (noeigstats_outfile_path, noeigstats_outfile_known_simple),
        (problem2_outfile_path, problem2_outfile_known_simple),
        (etot_etype_outfile_path, etot_etype_outfile_known_simple),
    ],
)
def test_JDFTXOutfile_fromfile_simple(filename: Path, known: dict):
    jdftxoutfile_fromfile_matches_known_simple(filename, known)


@pytest.mark.parametrize(
    ("filename", "latknown"),
    [
        (partial_lattice_init_outfile_path, partial_lattice_init_outfile_known_lattice),
    ],
)
def test_JDFTXOutfile_default_struc_inheritance(filename: Path, latknown: dict):
    jof = JDFTXOutfile.from_file(filename)
    lattice = jof.lattice
    for i in range(3):
        for j in range(3):
            assert pytest.approx(lattice[i][j]) == latknown[f"{i}{j}"]


# Make sure all possible exceptions are caught when none_on_error is True
@pytest.mark.parametrize(("ex_outfile_path"), [(partial_lattice_init_outfile_path)])
def test_none_on_partial(ex_outfile_path: Path):
    texts = read_outfile_slices(str(ex_outfile_path))
    texts0 = texts[:-1]

    slices = [
        JDFTXOutfileSlice._from_out_slice(text, is_bgw=False, none_on_error=False) for i, text in enumerate(texts0)
    ]
    outfile1 = JDFTXOutfile(slices=slices)
    slices.append(None)
    outfile2 = JDFTXOutfile(slices=slices)
    assert isinstance(outfile2, JDFTXOutfile)
    for var in _jof_atr_from_last_slice:
        assert_same_value(
            getattr(outfile1, var),
            getattr(outfile2, var),
        )


@pytest.mark.parametrize(
    ("example_vib_outfile_path", "example_vib_modes_known", "example_vib_nrg_components"),
    [
        (example_vib_outfile_path, example_vib_modes_known, example_vib_nrg_components),
    ],
)
def test_vib_parse(
    example_vib_outfile_path: Path, example_vib_modes_known: list[dict], example_vib_nrg_components: dict
):
    """
    Test that the vibration modes are parsed correctly from the outfile.
    """
    jdftxoutfile = JDFTXOutfile.from_file(example_vib_outfile_path)
    assert_same_value(jdftxoutfile.slices[-1].vibrational_modes, example_vib_modes_known)
    assert_same_value(jdftxoutfile.slices[-1].vibrational_energy_components, example_vib_nrg_components)
