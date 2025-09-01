"""Shared test utilities for JDFTx output files.

This module contains shared testing functions and example data + known values for JDFTx output files.
This module will be combined with shared_test_utils.py upon final implementation.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest

from pymatgen.core.units import Ha_to_eV, bohr_to_ang
from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice
from pymatgen.io.jdftx.outputs import JDFTXOutfile
from pymatgen.util.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from collections.abc import Callable


from .shared_test_utils import assert_same_value

# def write_mt_file(fname: str, write_dir: Path = dump_files_dir):
#     filepath = write_dir / fname
#     with open(filepath, "w", encoding="utf-8") as f:
#         f.write("if you're reading this yell at ben")
#     f.close()


def object_hasall_known_simple(obj: Any, knowndict: dict):
    for k in knowndict:
        assert hasattr(obj, k)


def object_matchall_known_simple(obj: Any, knowndict: dict):
    for k, v in knowndict.items():
        val = getattr(obj, k)
        assert_same_value(val, v)


def assert_slices_attribute_error(
    init_meth: Callable, init_var: Any, varname: str, slicename: str, assert_2layer_error: bool = False
):
    """Assert raises AttributeError upon certain conditions for slices attribute.

    Assert that varname property inherited from object's 'slices' type class variable
    will always raise AttributeError if the 'slices' attribute is empty. If assert_2layer_error,
    also assert that that final slice does not return None for that varname.

    Parameters:
    ----------
    init_meth: callable
        The method to initialize the object with.
    init_var: Any
        The variable to initialize the object with.
    varname: str
        The name of the attribute to test.
    slicename: str
        The name of the attribute that is a list of slices.
    assert_2layer_error: bool
        If True, assert that varname will raise AttributeError if the last slice's varname is None.
        (The attribute will not accept None from the last slice)
        (If returning None is okay as long as the slices are set, set this to False)
    """
    obj = init_meth(init_var)
    getattr(obj, varname)  # No freakout here
    setattr(getattr(obj, slicename)[-1], varname, None)
    if assert_2layer_error:
        with pytest.raises(AttributeError):
            getattr(obj, varname)
    else:
        getattr(obj, varname)  # No freakout here
    setattr(obj, slicename, [])
    with pytest.raises(AttributeError):
        getattr(obj, varname)


def assert_slices_1layer_attribute_error(init_meth: Callable, init_var: Any, varname: str, slicename: str):
    assert_slices_attribute_error(init_meth, init_var, varname, slicename, assert_2layer_error=False)


def assert_slices_2layer_attribute_error(init_meth: Callable, init_var: Any, varname: str, slicename: str):
    assert_slices_attribute_error(init_meth, init_var, varname, slicename, assert_2layer_error=True)


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


ex_out_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "test_jdftx_out_files"
ex_out_file_sections_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "test_jdftx_out_file_sections"
ex_calc_dirs_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "test_jdftx_calc_dirs"

n2_ex_calc_dir = ex_calc_dirs_dir / Path("N2")
n2_ex_calc_dir_known_paths = {
    "bandProjections": n2_ex_calc_dir / Path("bandProjections"),
    "eigenvals": n2_ex_calc_dir / Path("eigenvals"),
}
n2_ex_calc_dir_bandprojections_metadata = {
    "orb_label_list": ["N#1(s)", "N#1(py)", "N#1(pz)", "N#1(px)", "N#2(s)", "N#2(py)", "N#2(pz)", "N#2(px)"],
    "shape": (54, 15, 8),
    "first val": -0.1331527 + 0.5655596j,
}


nh3_ex_calc_dir = ex_calc_dirs_dir / Path("NH3")
nh3_ex_calc_dir_known_paths = {
    "bandProjections": nh3_ex_calc_dir / Path("bandProjections"),
    "eigenvals": nh3_ex_calc_dir / Path("eigenvals"),
}
nh3_ex_calc_dir_bandprojections_metadata = {
    "orb_label_list": ["N#1(s)", "N#1(py)", "N#1(pz)", "N#1(px)", "H#1(s)", "H#2(s)", "H#3(s)"],
    "shape": (16, 14, 7),
    "first val": -0.0688767 + 0.9503786j,
}

example_sp_outfile_path = ex_out_files_dir / Path("example_sp.out")
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
    "fluid": None,
    "total_electrons": 288.0,
    "nbands": 174,
    "nat": 16,
    "t_s": 165.87,
    "geom_opt_type": "single point",
    "prefix": "jdft",
    "etype": "F",
    "converged": True,
    "ecomponents": example_sp_outfile_known_ecomp,
}

example_latmin_outfile_path = ex_out_files_dir / Path("example_latmin.out")
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
    "fluid": None,
    "total_electrons": 64.0,
    "nbands": 42,
    "nat": 8,
    "t_s": 314.16,
    "geom_opt_type": "lattice",
    "prefix": "$VAR",
    "etype": "F",
    "converged": True,
    "ecomponents": example_latmin_outfile_known_ecomp,
}

example_ionmin_outfile_path = ex_out_files_dir / Path("example_ionmin.out")
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
    "geom_opt_type": "ionic",
    "prefix": "$VAR",
    "etype": "G",
    "converged": True,
    "ecomponents": example_ionmin_outfile_known_ecomp,
}

noeigstats_outfile_path = ex_out_files_dir / Path("noeigstats.out")
noeigstats_outfile_known_simple = {
    "mu": -0.050095169 * Ha_to_eV,
    "efermi": -0.050095169 * Ha_to_eV,
}

problem2_outfile_path = ex_out_files_dir / Path("problem2.out")
problem2_outfile_known_simple = {
    "mu": 0.464180124 * Ha_to_eV,
}

etot_etype_outfile_path = ex_out_files_dir / Path("etot_etype.out")
etot_etype_outfile_known_simple = {
    "e": -17.265553748795949 * Ha_to_eV,
    "elec_grad_k": 2.991e-07,
}

partial_lattice_init_outfile_path = ex_out_files_dir / Path("partial_lattice_init.out")
partial_lattice_init_outfile_known_lattice = {
    "00": 13.850216000000000 * bohr_to_ang,
    "01": 0.000000000000000 * bohr_to_ang,
    "02": -0.297459000000000 * bohr_to_ang,
    "10": -4.625257000000000 * bohr_to_ang,
    "11": 13.055094000000000 * bohr_to_ang,
    "12": -0.297459000000000 * bohr_to_ang,
    "20": 0.000000000000000 * bohr_to_ang,
    "21": 0.000000000000000 * bohr_to_ang,
    "22": 54.648857000000000 * bohr_to_ang,
}


ex_outfileslice1_fname = ex_out_file_sections_dir / "ex_out_slice_latmin"
ex_outfileslice2_fname = ex_out_file_sections_dir / "ex_out_slice_ionmin"
with open(ex_outfileslice1_fname, encoding="utf-8") as f:
    ex_outfileslice1 = list.copy(list(f))
with open(ex_outfileslice2_fname, encoding="utf-8") as f:
    ex_outfileslice2 = list.copy(list(f))
ex_outfileslice1_known = {
    "mu0_0": 0.713855355 * Ha_to_eV,
    "mu0_-1": 0.703866408 * Ha_to_eV,
    "nEminSteps0": 18,
    "etype0": "F",
    "E0": -246.531007900240667 * Ha_to_eV,
    "conv0": True,
    "mu-1_0": 0.704400512 * Ha_to_eV,
    "mu-1_-1": 0.704399109 * Ha_to_eV,
    "nEminSteps-1": 4,
    "etype-1": "F",
    "E-1": -246.531042396724303 * Ha_to_eV,
    "nGeomSteps": 7,
    "conv-1": True,
    "nelec0_0": 64.0,
    "nelec0_-1": 64.0,
    "nelec-1_0": 64.0,
    "nelec-1_-1": 64.0,
}
ex_outfileslice2_known = {
    "mu0_0": -0.190000000 * Ha_to_eV,
    "mu0_-1": -0.190000000 * Ha_to_eV,
    "nEminSteps0": 101,
    "etype0": "G",
    "E0": -1058.990493255521415 * Ha_to_eV,
    "conv0": False,
    "mu-1_0": -0.190000000 * Ha_to_eV,
    "mu-1_-1": -0.190000000 * Ha_to_eV,
    "nEminSteps-1": 13,
    "etype-1": "G",
    "E-1": -1059.062593502930213 * Ha_to_eV,
    "nGeomSteps": 7,
    "conv-1": True,
    "nelec0_0": 325.000000,
    "nelec0_-1": 325.610434,
    "nelec-1_0": 325.541001,
    "nelec-1_-1": 325.541406,
}

ex_jstruc_slice_fname1 = ex_out_file_sections_dir / "ex_text_slice_forJAtoms_latmin"
ex_jstruc_slice1 = []
with open(ex_jstruc_slice_fname1, encoding="utf-8") as f:
    ex_jstruc_slice1 = list.copy(list(f))

ex_jstruc_slice1_known = {
    "opt_type": "lattice",
    "nstep": 0,
    "etype": "F",
    "E": -246.5310079002406667 * Ha_to_eV,
    "Eewald": -214.6559882144248945 * Ha_to_eV,
    "EH": 28.5857387723713110 * Ha_to_eV,
    "Eloc": -40.1186842665999635 * Ha_to_eV,
    "Enl": -69.0084493129606642 * Ha_to_eV,
    "EvdW": -0.1192533377321287 * Ha_to_eV,
    "Exc": -90.7845534796796727 * Ha_to_eV,
    "Exc_core": 50.3731883713289008 * Ha_to_eV,
    "KE": 89.1972709081141488 * Ha_to_eV,
    "Etot": -246.5307305595829348 * Ha_to_eV,
    "TS": 0.0002773406577414 * Ha_to_eV,
    "F": -246.5310079002406667 * Ha_to_eV,
    "mu0": 0.713855355 * Ha_to_eV,
    "mu-1": 0.703866408 * Ha_to_eV,
    "E0": -246.455370884127575 * Ha_to_eV,
    "E-1": -246.531007900240667 * Ha_to_eV,
    "nEminSteps": 18,
    "EconvReason": "|Delta F|<1.000000e-07 for 2 iters",
    "conv": True,
    "cell_00": 6.16844 * bohr_to_ang,
    "strain_00": 10.0,
    "stress_00": -1.69853e-06 * Ha_to_eV / bohr_to_ang**3,
    "nAtoms": 8,
    "posn0": np.array([0.000011000000000, 2.394209000000000, 1.474913000000000]) * bohr_to_ang,
    "force0": np.array([0.000003219385226, 0.000024941936105, -0.000004667309539]) * (Ha_to_eV / bohr_to_ang),
    "posn-1": np.array([0.000007000000000, 9.175312000000002, 4.423851000000000]) * bohr_to_ang,
    "force-1": np.array([0.000000021330734, -0.000015026361853, -0.000010315177459]) * (Ha_to_eV / bohr_to_ang),
    "ox0": 0.048,
    "mag0": 0.000,
    "ox-1": -0.034,
    "mag-1": 0.000,
}


ex_jstruc_slice_fname2 = ex_out_file_sections_dir / "ex_text_slice_forJAtoms_latmin2"
ex_jstruc_slice2 = []
with open(ex_jstruc_slice_fname2, encoding="utf-8") as f:
    ex_jstruc_slice2 = list.copy(list(f))


ex_jstruc_slice2_known = {
    "opt_type": "lattice",
    "nstep": 9,
    "etype": "F",
    "E": -246.5310079002406667 * Ha_to_eV,
    "Eewald": -214.6559882144248945 * Ha_to_eV,
    "EH": 28.5857387723713110 * Ha_to_eV,
    "Eloc": -40.1186842665999635 * Ha_to_eV,
    "Enl": -69.0084493129606642 * Ha_to_eV,
    "EvdW": -0.1192533377321287 * Ha_to_eV,
    "Exc": -90.7845534796796727 * Ha_to_eV,
    "Exc_core": 50.3731883713289008 * Ha_to_eV,
    "KE": 89.1972709081141488 * Ha_to_eV,
    "Etot": -246.5307305595829348 * Ha_to_eV,
    "TS": 0.0002773406577414 * Ha_to_eV,
    "F": -246.5310079002406667 * Ha_to_eV,
    "mu0": 1.713855355 * Ha_to_eV,
    "mu-1": 0.703866408 * Ha_to_eV,
    "E0": -246.455370884127575 * Ha_to_eV,
    "E-1": -246.531007900240667 * Ha_to_eV,
    "nEminSteps": 18,
    "EconvReason": "|Delta F|<1.000000e-07 for 2 iters",
    "conv": True,
    "cell_00": 6.16844 * bohr_to_ang,
    "strain_00": 10.0,
    "stress_00": -1.69853e-06 * Ha_to_eV / bohr_to_ang**3,
    "nAtoms": 8,
    "posn0": np.array([0.000011000000000, 2.394209000000000, 1.474913000000000]) * bohr_to_ang,
    "force0": np.array([0.000003219385226, 0.000024941936105, -0.000004667309539]) * (Ha_to_eV / bohr_to_ang),
    "posn-1": np.array([0.000007000000000, 9.175312000000002, 4.423851000000000]) * bohr_to_ang,
    "force-1": np.array([0.000000021330734, -0.000015026361853, -0.000010315177459]) * (Ha_to_eV / bohr_to_ang),
    "ox0": 0.048,
    "mag0": 0.000,
    "ox-1": -0.034,
    "mag-1": 0.100,
}

# JESteps test knowns


ex_fillings_line1 = "FillingsUpdate:  mu: +0.714406772  \
    nElectrons: 64.000000  magneticMoment: [ Abs: 0.00578  Tot: -0.00141 ]"
ex_fillings_line1_known = {
    "mu": 0.714406772 * Ha_to_eV,
    "nelectrons": 64.0,
    "abs_magneticmoment": 0.00578,
    "tot_magneticmoment": -0.00141,
}

ex_fillings_line2 = "FillingsUpdate:  mu: +0.814406772  \
    nElectrons: 60.000000  magneticMoment: [ Abs: 0.0578  Tot: -0.0141 ]"
ex_fillings_line2_known = {
    "mu": 0.814406772 * Ha_to_eV,
    "nelectrons": 60.0,
    "abs_magneticmoment": 0.0578,
    "tot_magneticmoment": -0.0141,
}

ex_subspace_line1 = "SubspaceRotationAdjust: set factor to 0.229"
ex_subspace_line1_known = {"subspacerotationadjust": 0.229}

ex_subspace_line2 = "SubspaceRotationAdjust: set factor to 0.329"
ex_subspace_line2_known = {"subspacerotationadjust": 0.329}

ex_iter_line1 = "ElecMinimize: Iter:   6  F: -246.531038317370076\
        |grad|_K:  6.157e-08  alpha:  5.534e-01  linmin: -4.478e-06\
              t[s]:    248.68"
ex_iter_line1_known = {
    "nstep": 6,
    "e": -246.531038317370076 * Ha_to_eV,
    "grad_k": 6.157e-08,
    "alpha": 5.534e-01,
    "linmin": -4.478e-06,
    "t_s": 248.68,
}

ex_iter_line2 = "ElecMinimize: Iter:   7  F: -240.531038317370076\
        |grad|_K:  6.157e-07  alpha:  5.534e-02  linmin: -5.478e-06\
                t[s]:    48.68"
ex_iter_line2_known = {
    "nstep": 7,
    "e": -240.531038317370076 * Ha_to_eV,
    "grad_k": 6.157e-07,
    "alpha": 5.534e-02,
    "linmin": -5.478e-06,
    "t_s": 48.68,
}


ex_jstep_lines1 = [ex_fillings_line1, ex_subspace_line1, ex_iter_line1]
ex_jstep_lines2 = [ex_fillings_line2, ex_subspace_line2, ex_iter_line2]
ex_jstep_known1 = {}
for known1 in [ex_fillings_line1_known, ex_iter_line1_known, ex_subspace_line1_known]:
    ex_jstep_known1.update(known1)
ex_jstep_known2 = {}
for known2 in [ex_fillings_line2_known, ex_iter_line2_known, ex_subspace_line2_known]:
    ex_jstep_known2.update(known2)


example_vib_outfile_path = ex_out_files_dir / Path("vib.out")
example_vib_nrg_components = {
    "T": 298.0,
    "ZPE": 0.012162 * Ha_to_eV,
    "Evib": 0.012641 * Ha_to_eV,
    "TSvib": 0.000686 * Ha_to_eV,
    "Avib": 0.011954 * Ha_to_eV,
}
example_vib_modes_known = [
    {
        "Type": "Imaginary",
        "Type index": 1,
        "Frequency": 0.003211 * 1j * Ha_to_eV,
        "Degeneracy": 1,
        "IR intensity": 0.1015,
        "Displacements": np.array(
            [
                np.array([-0.001502301500021, 0.016139807700462, 0.004794115920733]),
                np.array([-0.003227666066075, -0.003444853384899, -0.005623654997270]),
                np.array([0.003563062362181, -0.009568755195267, -0.010017876395276]),
            ]
        )
        * bohr_to_ang,
    },
    {
        "Type": "Imaginary",
        "Type index": 2,
        "Frequency": 0.001816 * 1j * Ha_to_eV,
        "Degeneracy": 1,
        "IR intensity": 0.0014,
        "Displacements": np.array(
            [
                np.array([0.011897355539145, 0.001393779704454, -0.010870152434402]),
                np.array([-0.014950015017923, 0.003952780094546, 0.004794677356827]),
                np.array([0.002882750977633, -0.003471024246565, 0.000366070625853]),
            ]
        )
        * bohr_to_ang,
    },
    {
        "Type": "Imaginary",
        "Type index": 3,
        "Frequency": 0.000986 * 1j * Ha_to_eV,
        "Degeneracy": 1,
        "IR intensity": 0.0003,
        "Displacements": np.array(
            [
                np.array([0.001831896775059, -0.000519493708404, -0.000826114063229]),
                np.array([0.008215382766313, 0.013135969275973, 0.003002047268263]),
                np.array([-0.008612872298526, -0.014707651990282, 0.000628692658427]),
            ]
        )
        * bohr_to_ang,
    },
    {
        "Type": "Imaginary",
        "Type index": 4,
        "Frequency": 0.000584 * 1j * Ha_to_eV,
        "Degeneracy": 1,
        "IR intensity": 0.0009,
        "Displacements": np.array(
            [
                np.array([-0.002873587752642, -0.000357270159722, 0.009507035658568]),
                np.array([-0.002773621141597, -0.002278652939269, 0.008444977841195]),
                np.array([0.010719232053596, -0.008451829956361, 0.013226919526156]),
            ]
        )
        * bohr_to_ang,
    },
    {
        "Type": "Zero",
        "Type index": 1,
        "Frequency": 0.000010 * 1j * Ha_to_eV,
        "Degeneracy": 2,
        "IR intensity": 0.0000,
        "Displacements": np.array(
            [
                np.array([0.007597327520746, 0.010712391019933, 0.002704566793247]),
                np.array([0.007881150132850, 0.010845838899770, 0.002620013897312]),
                np.array([0.007765190577341, 0.010591653551291, 0.002319182104746]),
            ]
        )
        * bohr_to_ang,
    },
    {
        "Type": "Zero",
        "Type index": 2,
        "Frequency": 0.000050 * Ha_to_eV,
        "Degeneracy": 2,
        "IR intensity": 0.0000,
        "Displacements": np.array(
            [
                np.array([0.010928992861829, -0.007092818205113, -0.003874657684275]),
                np.array([0.010923114633010, -0.007525781466854, -0.004091667059576]),
                np.array([0.010694926666458, -0.006058360112810, -0.003964322783956]),
            ]
        )
        * bohr_to_ang,
    },
    {
        "Type": "Real",
        "Type index": 1,
        "Frequency": 0.001795 * Ha_to_eV,
        "Degeneracy": 1,
        "IR intensity": 0.0050,
        "Displacements": np.array(
            [
                np.array([0.004689955512831, -0.004561918723060, 0.010223393583369]),
                np.array([-0.001077029836943, -0.001525888489046, 0.014555134211414]),
                np.array([-0.002282652762669, 0.001952724669842, -0.013136592620777]),
            ]
        )
        * bohr_to_ang,
    },
    {
        "Type": "Real",
        "Type index": 2,
        "Frequency": 0.002701 * Ha_to_eV,
        "Degeneracy": 1,
        "IR intensity": 0.0896,
        "Displacements": np.array(
            [
                np.array([0.011841920013071, 0.006471479676558, 0.004038350981094]),
                np.array([0.002080302463227, -0.010601817998144, 0.000214161706513]),
                np.array([-0.012420751778051, 0.000399317651251, 0.008639012103301]),
            ]
        )
        * bohr_to_ang,
    },
    {
        "Type": "Real",
        "Type index": 3,
        "Frequency": 0.019827 * Ha_to_eV,
        "Degeneracy": 1,
        "IR intensity": 0.0376,
        "Displacements": np.array(
            [
                np.array([0.007032223801816, -0.007326012017446, 0.012995955668716]),
                np.array([-0.006945332786266, 0.007072807313742, -0.013183707786736]),
                np.array([-0.000064297833090, 0.000203820782562, 0.000350655054206]),
            ]
        )
        * bohr_to_ang,
    },
]
