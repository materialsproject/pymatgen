from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from pytest import approx

from pymatgen.core.units import Ha_to_eV
from pymatgen.io.jdftx.jdftxoutfile import JDFTXOutfile
from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice
from pymatgen.util.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from pymatgen.util.typing import PathLike

ex_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "example_files"

test_read = JDFTXOutfile.from_file(ex_files_dir / Path("problem1.out"))

example_sp_known = {
    "Nspin": 1,
    "spintype": "no-spin",
    "broadening_type": "MP1",
    "broadening": 0.00367493,
    "truncation_type": "slab",
    "pwcut": 30 * Ha_to_eV,
    "fftgrid": (54, 54, 224),
    "kgrid": (6, 6, 1),
    "Emin": -3.836283 * Ha_to_eV,
    "HOMO": -0.212435 * Ha_to_eV,
    "EFermi": -0.209509 * Ha_to_eV,
    "LUMO": -0.209424 * Ha_to_eV,
    "Emax": 0.113409 * Ha_to_eV,
    "Egap": 0.003011 * Ha_to_eV,
    "is_metal": True,
    "fluid": "None",
    "total_electrons": 288.0,
    "nbands": 174,
    "nat": 16,
    "F": -1940.762261217305650 * Ha_to_eV,
    "TS": -0.0001776512106456 * Ha_to_eV,
    "Etot": -1940.7624388685162558 * Ha_to_eV,
    "KE": 593.1822417205943339 * Ha_to_eV,
    "Exc": -185.5577583222759870 * Ha_to_eV,
    "Epulay": 0.0000125227478554 * Ha_to_eV,
    "Enl": 174.1667582919756114 * Ha_to_eV,
    "Eloc": 29663.3545152997867262 * Ha_to_eV,
    "EH": -15284.4385436602351547 * Ha_to_eV,
    "Eewald": -16901.4696647211094387 * Ha_to_eV,
    "nSlices": 1,
    "t_s": 165.87,
    "iter_type": None,
    "prefix": "jdft",
    "etype": "F",
    "converged": True,
}

example_latmin_known = {
    "Nspin": 2,
    "spintype": "z-spin",
    "broadening_type": "Fermi",
    "broadening": 0.001,
    "truncation_type": "periodic",
    "pwcut": 20 * Ha_to_eV,
    "fftgrid": (28, 80, 28),
    "kgrid": (6, 2, 7),
    "Emin": -1.780949 * Ha_to_eV,
    "HOMO": 0.704289 * Ha_to_eV,
    "EFermi": 0.704399 * Ha_to_eV,
    "LUMO": 0.704651 * Ha_to_eV,
    "Emax": 0.949497 * Ha_to_eV,
    "Egap": 0.000362 * Ha_to_eV,
    "is_metal": True,
    "fluid": "None",
    "total_electrons": 64.0,
    "nbands": 42,
    "nat": 8,
    "F": -246.5310423967243025 * Ha_to_eV,
    "TS": 0.0003221374940495 * Ha_to_eV,
    "Etot": -246.5307202592302644 * Ha_to_eV,
    "KE": 89.2073662863590755 * Ha_to_eV,
    "Exc": -90.7880124097588208 * Ha_to_eV,
    "Enl": -69.0117974720974559 * Ha_to_eV,
    "Eloc": -40.0429414587348518 * Ha_to_eV,
    "EH": 28.5721759138337354 * Ha_to_eV,
    "Eewald": -214.7213057123609019 * Ha_to_eV,
    "nSlices": 7,
    "t_s": 314.16,
    "iter_type": "LatticeMinimize",
    "prefix": "$VAR",
    "etype": "F",
    "converged": True,
}

example_ionmin_known = {
    "Nspin": 2,
    "spintype": "z-spin",
    "broadening_type": "Fermi",
    "broadening": 0.001,
    "truncation_type": "slab",
    "pwcut": 25 * Ha_to_eV,
    "fftgrid": (56, 56, 320),
    "kgrid": (4, 4, 1),
    "Emin": -2.488051 * Ha_to_eV,
    "HOMO": -0.190949 * Ha_to_eV,
    "EFermi": -0.190000 * Ha_to_eV,
    "LUMO": -0.189724 * Ha_to_eV,
    "Emax": -0.042437 * Ha_to_eV,
    "Egap": 0.001225 * Ha_to_eV,
    "is_metal": False,  # Oh god oh god oh god
    "fluid": "LinearPCM",
    "total_electrons": 325.541406,
    "nbands": 195,
    "nat": 41,
    "G": -1059.062593502930213 * Ha_to_eV,
    "F": -1120.9154606162035179 * Ha_to_eV,
    "TS": 0.0014609776617570 * Ha_to_eV,
    "Etot": -1120.9139996385417817 * Ha_to_eV,
    "KE": 421.4844651353773770 * Ha_to_eV,
    "Exc": -796.7101488293942566 * Ha_to_eV,
    "Enl": -270.1618154209642739 * Ha_to_eV,
    "Eloc": -79647.5920994735934073 * Ha_to_eV,
    "EH": 39775.3166089357473538 * Ha_to_eV,
    "Eewald": 38803.1912795634780196 * Ha_to_eV,
    "nSlices": 1,
    "t_s": 2028.57,
    "iter_type": "IonicMinimize",
    "prefix": "$VAR",
    "etype": "G",
    "converged": True,
}


@pytest.mark.parametrize(
    ("filename", "known"),
    [
        (ex_files_dir / Path("example_sp.out"), example_sp_known),
        (ex_files_dir / Path("example_latmin.out"), example_latmin_known),
        (ex_files_dir / Path("example_ionmin.out"), example_ionmin_known),
    ],
)
def test_JDFTXOutfile_fromfile(filename: PathLike, known: dict):
    # filename = ex_files_dir / Path("jdftx.out")
    jout = JDFTXOutfile.from_file(filename)
    assert jout.nspin == known["Nspin"]
    assert jout["nspin"] == known["Nspin"]  # Once for the coverage
    assert isinstance(jout[-1], JDFTXOutfileSlice)
    with pytest.raises(TypeError):
        jout[{}]
    assert jout.spintype == known["spintype"]
    assert jout.broadening_type == known["broadening_type"]
    assert jout.broadening == approx(known["broadening"])
    assert jout.truncation_type == known["truncation_type"]
    assert jout.pwcut == approx(known["pwcut"])
    # Don't bully me, I'm testing this way in case we flip-flop between lists and tuples
    for i in range(3):
        assert jout.fftgrid[i] == known["fftgrid"][i]
    for i in range(3):
        assert jout.kgrid[i] == known["kgrid"][i]
    assert jout.emin == approx(known["Emin"])
    assert approx(known["HOMO"]) == jout.homo
    assert jout.efermi == approx(known["EFermi"])
    assert approx(known["LUMO"]) == jout.lumo
    assert jout.emax == approx(known["Emax"])
    assert jout.egap == approx(known["Egap"])
    # TODO: filling tests
    # assert jout.HOMO_filling == approx(None)
    # assert jout.LUMO_filling == approx(None)
    assert jout.is_metal == known["is_metal"]
    assert jout.fluid == known["fluid"]
    assert jout.total_electrons == approx(known["total_electrons"])
    assert jout.nbands == known["nbands"]
    assert jout.nat == known["nat"]
    for listlike in (
        jout.atom_coords,
        jout.atom_coords_final,
        jout.atom_coords_initial,
        jout.atom_elements,
        jout.atom_elements_int,
    ):
        assert len(listlike) == known["nat"]
    assert jout.ecomponents["F"] == approx(known["F"])
    assert jout.ecomponents["TS"] == approx(known["TS"])
    assert jout.ecomponents["Etot"] == approx(known["Etot"])
    assert jout.ecomponents["KE"] == approx(known["KE"])
    assert jout.ecomponents["Exc"] == approx(known["Exc"])
    assert jout.ecomponents["Enl"] == approx(known["Enl"])
    assert jout.ecomponents["Eloc"] == approx(known["Eloc"])
    assert jout.ecomponents["EH"] == approx(known["EH"])
    assert jout.ecomponents["Eewald"] == approx(known["Eewald"])
    assert len(jout.slices) == known["nSlices"]
    assert jout.t_s == approx(known["t_s"])
    assert jout.jstrucs.iter_type == known["iter_type"]
    assert jout.prefix == known["prefix"]
    assert jout.etype == known["etype"]
    assert jout.e == approx(known[known["etype"]])
    assert jout.converged == known["converged"]
    # Not testing values yet, just testing they dont raise errors
    assert jout.trajectory is not None
    assert jout.electronic_output is not None
    assert jout.structure is not None
    jout[-1].jstrucs = None
    assert jout.is_converged is None


noeigstats_known = {
    "mu": -0.050095169 * Ha_to_eV,
    "efermi": -0.050095169 * Ha_to_eV,
}


@pytest.mark.parametrize(
    ("outfilefname", "knowndict"),
    [
        (ex_files_dir / Path("noeigstats.out"), noeigstats_known),
    ],
)
def test_jdftxoutfile_knowns_simple(outfilefname, knowndict):
    jout = JDFTXOutfile.from_file(outfilefname)
    for k, v in knowndict.items():
        assert hasattr(jout, k)
        val = getattr(jout, k)
        if isinstance(v, float):
            assert val == pytest.approx(v)
        elif v is None:
            assert val is None
        else:
            assert val == v


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
    "iter",
    "e",
    "grad_k",
    "alpha",
    "linmin",
    "nelectrons",
    "abs_magneticmoment",
    "tot_magneticmoment",
    "mu",
    "elec_iter",
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
                "iter",
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


# Exception - trauncation_radius only set for some calculations and is thus None when unset
# Linmin also None sometimes

# Some only should be set for geometric optimizations and excluded from single point out file
# stress, strain, iter
