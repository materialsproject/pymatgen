from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from pymatgen.core.units import Ha_to_eV
from pymatgen.io.jdftx.jdftxoutfile import JDFTXOutfile
from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice
from pymatgen.util.testing import TEST_FILES_DIR

################################################################################
# General methods and variables
################################################################################

ex_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "example_files"
dump_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "new_files"


def object_hasall_known_simple(obj: Any, knowndict: dict):
    for k in knowndict:
        assert hasattr(obj, k)


def object_matchall_known_simple(obj: Any, knowndict: dict):
    for k, v in knowndict.items():
        val = getattr(obj, k)
        _is_same_value(val, v)


def _is_same_value(testval, knownval):
    if type(testval) not in [tuple, list]:
        assert isinstance(testval, type(knownval))
        if isinstance(testval, float):
            assert testval == pytest.approx(knownval)
        elif isinstance(testval, dict):
            for k in knownval:
                assert k in testval
                _is_same_value(testval[k], knownval[k])
        elif testval is None:
            assert knownval is None
        else:
            assert testval == knownval
    else:
        assert len(testval) == len(knownval)
        for i in range(len(testval)):
            _is_same_value(testval[i], knownval[i])


################################################################################
# JDFTXOutfile test methods and ref/known pairs
################################################################################


def jdftxoutfile_fromfile_matches_known_simple(outfilefname: Path, knowndict: dict):
    joutfile = JDFTXOutfile.from_file(outfilefname)
    jdftxoutfile_matches_known_simple(joutfile, knowndict)
    del joutfile


def jdftxoutfile_matches_known_simple(joutfile: JDFTXOutfile, knowndict: dict):
    object_hasall_known_simple(joutfile, knowndict)
    object_matchall_known_simple(joutfile, knowndict)


# @pytest.fixture(autouse=True)
def jdftxoutfile_fromfile_matches_known(filename: Path, known: dict):
    joutfile = JDFTXOutfile.from_file(filename)
    jdftxoutfile_matches_known(joutfile, known)
    del joutfile


def jdftxoutfile_matches_known(joutfile: JDFTXOutfile, known: dict):
    assert isinstance(joutfile[-1], JDFTXOutfileSlice)
    with pytest.raises(TypeError):
        joutfile[{}]
    for listlike in (
        joutfile.atom_coords,
        joutfile.atom_coords_final,
        joutfile.atom_coords_initial,
        joutfile.atom_elements,
        joutfile.atom_elements_int,
    ):
        assert len(listlike) == known["nat"]
    assert len(joutfile.slices) == known["nSlices"]
    # Not testing values yet, just testing they dont raise errors
    assert joutfile.trajectory is not None
    assert joutfile.electronic_output is not None
    assert joutfile.structure is not None
    joutfile[-1].jstrucs = None
    assert joutfile.is_converged is None


example_sp_outfile_path = ex_files_dir / Path("example_sp.out")
example_sp_outfile_known = {
    "nat": 16,
    "nSlices": 1,
}
example_sp_outfile_known_ecomp = {
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
}
example_sp_outfile_known_simple = {
    "nspin": 1,
    "spintype": "no-spin",
    "broadening_type": "MP1",
    "broadening": 0.00367493,
    "truncation_type": "slab",
    "pwcut": 30 * Ha_to_eV,
    "fftgrid": (54, 54, 224),
    "kgrid": (6, 6, 1),
    "emin": -3.836283 * Ha_to_eV,
    "homo": -0.212435 * Ha_to_eV,
    "efermi": -0.209509 * Ha_to_eV,
    "lumo": -0.209424 * Ha_to_eV,
    "emax": 0.113409 * Ha_to_eV,
    "egap": 0.003011 * Ha_to_eV,
    "is_metal": True,
    "fluid": "None",
    "total_electrons": 288.0,
    "nbands": 174,
    "nat": 16,
    "t_s": 165.87,
    "iter_type": None,
    "prefix": "jdft",
    "etype": "F",
    "converged": True,
    "ecomponents": example_sp_outfile_known_ecomp,
}

example_latmin_outfile_path = ex_files_dir / Path("example_latmin.out")
example_latmin_outfile_known = {
    "nat": 8,
    "nSlices": 7,
}
example_latmin_outfile_known_ecomp = {
    "F": -246.5310423967243025 * Ha_to_eV,
    "TS": 0.0003221374940495 * Ha_to_eV,
    "Etot": -246.5307202592302644 * Ha_to_eV,
    "KE": 89.2073662863590755 * Ha_to_eV,
    "Exc": -90.7880124097588208 * Ha_to_eV,
    "Enl": -69.0117974720974559 * Ha_to_eV,
    "Eloc": -40.0429414587348518 * Ha_to_eV,
    "EH": 28.5721759138337354 * Ha_to_eV,
    "Eewald": -214.7213057123609019 * Ha_to_eV,
}
example_latmin_outfile_known_simple = {
    "nspin": 2,
    "spintype": "z-spin",
    "broadening_type": "Fermi",
    "broadening": 0.001,
    "truncation_type": "periodic",
    "pwcut": 20 * Ha_to_eV,
    "fftgrid": (28, 80, 28),
    "kgrid": (6, 2, 7),
    "emin": -1.780949 * Ha_to_eV,
    "homo": 0.704289 * Ha_to_eV,
    "efermi": 0.704399 * Ha_to_eV,
    "lumo": 0.704651 * Ha_to_eV,
    "emax": 0.949497 * Ha_to_eV,
    "egap": 0.000362 * Ha_to_eV,
    "is_metal": True,
    "fluid": "None",
    "total_electrons": 64.0,
    "nbands": 42,
    "nat": 8,
    "t_s": 314.16,
    "iter_type": "LatticeMinimize",
    "prefix": "$VAR",
    "etype": "F",
    "converged": True,
    "ecomponents": example_latmin_outfile_known_ecomp,
}

example_ionmin_outfile_path = ex_files_dir / Path("example_ionmin.out")
example_ionmin_outfile_known = {
    "nat": 41,
    "nSlices": 1,
}
example_ionmin_outfile_known_ecomp = {
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
}
example_ionmin_outfile_known_simple = {
    "nspin": 2,
    "spintype": "z-spin",
    "broadening_type": "Fermi",
    "broadening": 0.001,
    "truncation_type": "slab",
    "pwcut": 25 * Ha_to_eV,
    "fftgrid": (56, 56, 320),
    "kgrid": (4, 4, 1),
    "emin": -2.488051 * Ha_to_eV,
    "homo": -0.190949 * Ha_to_eV,
    "efermi": -0.190000 * Ha_to_eV,
    "lumo": -0.189724 * Ha_to_eV,
    "emax": -0.042437 * Ha_to_eV,
    "egap": 0.001225 * Ha_to_eV,
    "is_metal": False,  # Oh god oh god oh god
    "fluid": "LinearPCM",
    "total_electrons": 325.541406,
    "nbands": 195,
    "nat": 41,
    "t_s": 2028.57,
    "iter_type": "IonicMinimize",
    "prefix": "$VAR",
    "etype": "G",
    "converged": True,
    "ecomponents": example_ionmin_outfile_known_ecomp,
}

noeigstats_outfile_path = ex_files_dir / Path("noeigstats.out")
noeigstats_outfile_known_simple = {
    "mu": -0.050095169 * Ha_to_eV,
    "efermi": -0.050095169 * Ha_to_eV,
}
