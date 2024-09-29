from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from pymatgen.io.jdftx.jdftxoutfile import JDFTXOutfile

from .conftest import (
    ex_files_dir,
    example_ionmin_outfile_known,
    example_ionmin_outfile_known_simple,
    example_ionmin_outfile_path,
    example_latmin_outfile_known,
    example_latmin_outfile_known_simple,
    example_latmin_outfile_path,
    example_sp_outfile_known,
    example_sp_outfile_known_simple,
    example_sp_outfile_path,
    jdftxoutfile_fromfile_matches_known,
    jdftxoutfile_fromfile_matches_known_simple,
    noeigstats_outfile_known_simple,
    noeigstats_outfile_path,
)

if TYPE_CHECKING:
    from pymatgen.util.typing import PathLike


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
    ],
)
def test_JDFTXOutfile_fromfile_simple(filename: Path, known: dict):
    jdftxoutfile_fromfile_matches_known_simple(filename, known)


empty_slice_exception_varnames = [
    "prefix",
    "jstrucs",
    "jsettings_fluid",
    "jsettings_electronic",
    "jsettings_lattice",
    "jsettings_ionic",
    "xc_func",
    "lattice_initial",
    "lattice_final",
    "lattice",
    "a",
    "b",
    "c",
    "fftgrid",
    "geom_opt",
    "geom_opt_type",
    "efermi",
    "egap",
    "emin",
    "emax",
    "homo",
    "lumo",
    "homo_filling",
    "lumo_filling",
    "is_metal",
    "etype",
    "broadening_type",
    "broadening",
    "kgrid",
    "truncation_type",
    "truncation_radius",
    "pwcut",
    "rhocut",
    "pp_type",
    "total_electrons",
    "semicore_electrons",
    "valence_electrons",
    "total_electrons_uncharged",
    "semicore_electrons_uncharged",
    "valence_electrons_uncharged",
    "nbands",
    "atom_elements",
    "atom_elements_int",
    "atom_types",
    "spintype",
    "nspin",
    "nat",
    "atom_coords_initial",
    "atom_coords_final",
    "atom_coords",
    "has_solvation",
    "fluid",
    "is_gc",
    "eiter_type",
    "elecmindata",
    "stress",
    "strain",
    "niter",
    "e",
    "grad_k",
    "alpha",
    "linmin",
    "nelectrons",
    "abs_magneticmoment",
    "tot_magneticmoment",
    "mu",
    "elec_niter",
    "elec_e",
    "elec_grad_k",
    "elec_alpha",
    "elec_linmin",
]


@pytest.mark.parametrize(
    ("dummy_filename", "none_exceptions"),
    [
        (
            ex_files_dir / Path("example_sp.out"),
            [
                "niter",
                "stress",
                "strain",
                "linmin",
                "alpha",
                "truncation_radius",
                "grad_k",
                "abs_magneticmoment",
                "tot_magneticmoment",
                "elec_grad_k",
                "elec_alpha",
                "elec_linmin",
            ],
        )
    ],
)
def test_JDFTXOutfile_expected_exceptions_empty_slices(
    dummy_filename: PathLike,
    none_exceptions: list[str],
):
    for var in empty_slice_exception_varnames:
        expected_exceptions_empty_slices_test_varname(dummy_filename, var, none_exceptions)


def expected_exceptions_empty_slices_test_varname(jout_filename: str, varname: str, none_exceptions: list[str]):
    jout = JDFTXOutfile.from_file(jout_filename)
    # First test that the attribute can be called
    val = getattr(jout, varname)
    # Next test it was properly called, or that its None if None is expected (but not both)
    v1 = bool(val is not None)
    v2 = bool(varname in none_exceptions)
    assert v1 != v2
    # assert bool(val is not None) != bool(val in none_exceptions)
    # Next test the error when slices is empty
    jout.slices = []
    with pytest.raises(AttributeError):
        getattr(jout, varname)
