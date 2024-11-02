from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest

from pymatgen.core.trajectory import Trajectory
from pymatgen.core.units import Ha_to_eV, ang_to_bohr
from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice
from pymatgen.util.testing import TEST_FILES_DIR

ex_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "example_files"

ex_slice1_fname = ex_files_dir / "ex_out_slice_latmin"
ex_slice2_fname = ex_files_dir / "ex_out_slice_ionmin"
with open(ex_slice1_fname) as f:
    ex_slice1 = list.copy(list(f))
with open(ex_slice2_fname) as f:
    ex_slice2 = list.copy(list(f))


# @pytest.mark.parametrize(
#     ("out_slice", "varname"),
#     [
#         (ex_slice1, "structure"),
#         (ex_slice1, "eiter_type"),
#         (ex_slice1, "elecmindata"),
#         (ex_slice1, "stress"),
#         (ex_slice1, "strain"),
#         (ex_slice1, "nstep"),
#         (ex_slice1, "e"),
#         (ex_slice1, "grad_k"),
#         (ex_slice1, "alpha"),
#         (ex_slice1, "linmin"),
#         (ex_slice1, "abs_magneticmoment"),
#         (ex_slice1, "tot_magneticmoment"),
#         (ex_slice1, "mu"),
#         (ex_slice1, "elec_nstep"),
#         (ex_slice1, "elec_e"),
#         (ex_slice1, "elec_grad_k"),
#         (ex_slice1, "elec_alpha"),
#         (ex_slice1, "elec_linmin"),
#         (ex_slice1, "nelectrons"),
#     ],
# )
# def test_jdftxoutfileslice_has_1layer_jstrucs_freakout(out_slice: list[str], varname: str):
#     joutslice = JDFTXOutfileSlice.from_out_slice(out_slice)
#     getattr(joutslice, varname)  # No freakout here
#     joutslice.jstrucs = None
#     with pytest.raises(AttributeError):
#         getattr(joutslice, varname)  # Freakout here


def test_jdftxoutfileslice_stringify():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    out_str = str(joutslice)
    assert isinstance(out_str, str)
    assert len(out_str)


def test_jdftxoutfileslice_converge():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    assert joutslice.is_converged


def test_jdftxoutfileslice_trajectory():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    traj = joutslice.trajectory
    assert isinstance(traj, Trajectory)
    del joutslice.jsettings_lattice.params["niterations"]
    # joutslice.jsettings_lattice.params["niterations"] = None
    with pytest.raises(ValueError, match=re.escape("Unknown issue due to partial initialization of settings objects.")):
        traj = joutslice.trajectory


def test_get_broadeningvars():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    btype = "btype"
    bval = 1.0
    text = [f"elec-smearing {btype} {bval}"]
    broadening_type, broadening = joutslice.get_broadeningvars(text)
    assert broadening_type == btype
    assert broadening == pytest.approx(bval)
    broadening_type, broadening = joutslice.get_broadeningvars([])
    assert broadening_type is None
    assert broadening == pytest.approx(0.0)


def test_get_truncationvars():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    # with pytest.raises(ValueError, match="No truncation type found in out file."):
    #     joutslice.get_truncationvars([])
    joutslice.is_bgw = True
    with pytest.raises(ValueError, match="BGW slab Coulomb truncation must be along z!"):
        joutslice.get_truncationvars(["coulomb-interaction Slab 010"])
    with pytest.raises(ValueError, match="BGW wire Coulomb truncation must be periodic in z!"):
        joutslice.get_truncationvars(["coulomb-interaction Cylindrical 010"])
    with pytest.raises(ValueError, match="Problem with this truncation!"):
        joutslice.get_truncationvars(["coulomb-interaction barbie 010"])
    truncation_type, truncation_radius = joutslice.get_truncationvars(
        ["coulomb-interaction Spherical", "Initialized spherical truncation of radius 1.0"]
    )
    assert truncation_type == "spherical"
    assert truncation_radius == pytest.approx(1.0 / ang_to_bohr)


def test_get_rho_cutoff():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    text = ["elec-cutoff 1.0"]
    joutslice.pwcut = None
    rhocut = joutslice.get_rho_cutoff(text)
    assert joutslice.pwcut == pytest.approx(1.0 * Ha_to_eV)
    assert rhocut == pytest.approx(joutslice.pwcut * 4)


def test_get_eigstats_varsdict():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    evardict = joutslice.get_eigstats_varsdict([], "$VAR")
    for key in evardict:
        assert evardict[key] is None
    joutslice.set_eigvars([])
    for key in evardict:
        if key != "efermi":
            assert getattr(joutslice, key) is None
    for key in ["efermi", "mu"]:
        assert getattr(joutslice, key) is not None


def test_get_pp_type():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    # with pytest.raises(ValueError, match="Could not determine pseudopotential type from file name root/PAW"):
    #     joutslice.get_pp_type(["Reading pseudopotential file root/PAW:"])
    assert joutslice.get_pp_type(["Reading pseudopotential file root/PAW:"]) is None
    assert joutslice.get_pp_type(["Reading pseudopotential file not_SG15/GBRV"]) == "GBRV"
    assert joutslice.get_pp_type(["Reading pseudopotential file not_GBRV/SG15"]) == "SG15"


def test_set_pseudo_vars_t1():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    # Just need more bound sets than there are atom types
    text = [
        "Reading pseudopotential file not_SG15/GBRV",
        "10 valence electrons ",
        "",
        "Reading pseudopotential file not_SG15/GBRV",
        "10 valence electrons ",
        "",
        "Reading pseudopotential file not_SG15/GBRV",
        "10 valence electrons ",
        "",
        "Reading pseudopotential file not_SG15/GBRV",
        "10 valence electrons ",
        "",
        "Reading pseudopotential file not_SG15/GBRV",
        "10 valence electrons ",
        "",
        "Reading pseudopotential file not_SG15/GBRV",
        "10 valence electrons ",
        "",
    ]
    joutslice.total_electrons_backup = None
    joutslice.jstrucs = None
    with pytest.raises(ValueError, match="Total electrons and semicore electrons must be set."):
        joutslice.set_pseudo_vars_t1(text)
    joutslice.atom_elements = None
    with pytest.raises(ValueError, match="Atom elements not set yet."):
        joutslice.set_pseudo_vars_t1(text)
    joutslice.atom_types = None
    with pytest.raises(ValueError, match="Pseuopotential data cannot be allocated without atom types."):
        joutslice.set_pseudo_vars_t1(text)


def test_set_geomopt_vars():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    joutslice.jsettings_ionic = None
    with pytest.raises(ValueError, match="Unknown issue in setting settings objects"):
        joutslice.set_geomopt_vars([])


def test_set_orb_fillings_nobroad():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    joutslice.set_orb_fillings_nobroad(1)
    assert joutslice.homo_filling == pytest.approx(2)
    assert joutslice.lumo_filling == pytest.approx(0)


def test_set_orb_fillings_broad():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    joutslice.lumo = None
    with pytest.raises(ValueError, match="Cannot set orbital fillings with broadening with self.lumo as None"):
        joutslice.set_orb_fillings()
    joutslice.homo = None
    with pytest.raises(ValueError, match="Cannot set orbital fillings with broadening with self.homo as None"):
        joutslice.set_orb_fillings()
    joutslice.efermi = None
    with pytest.raises(ValueError, match="Cannot set orbital fillings with broadening with self.efermi as None"):
        joutslice.set_orb_fillings()
    joutslice.broadening = None
    with pytest.raises(ValueError, match="Cannot set orbital fillings with broadening with self.broadening as None"):
        joutslice.set_orb_fillings()
    joutslice.broadening_type = None
    joutslice.nspin = None
    with pytest.raises(ValueError, match="Cannot set homo/lumo filling with self.nspin as None"):
        joutslice.set_orb_fillings()


def test_set_lattice_vars():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    joutslice.jstrucs = None
    with pytest.raises(ValueError, match="No structures found in out file."):
        joutslice.set_lattice_vars([])


def test_set_ecomponents():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    joutslice.jstrucs = None
    with pytest.raises(ValueError, match="No structures found in out file."):
        joutslice.set_ecomponents([])


def test_calculate_filling():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    broadening = 1.0
    eig = 0.5
    efermi = 0.6
    x = (eig - efermi) / (2.0 * broadening)
    assert 0.5 * (1 - np.tanh(x)) == pytest.approx(joutslice.calculate_filling("Fermi", broadening, eig, efermi))
    assert 0.5 * (1 - math.erf(x)) - x * np.exp(-1 * x**2) / (2 * np.pi**0.5) == pytest.approx(
        joutslice.calculate_filling("MP1", broadening, eig, efermi)
    )
    assert 0.5 * (1 - math.erf(x)) == pytest.approx(joutslice.calculate_filling("Gauss", broadening, eig, efermi))
    assert 0.5 * (1 - math.erf(x + 0.5**0.5)) + np.exp(-1 * (x + 0.5**0.5) ** 2) / (2 * np.pi) ** 0.5 == pytest.approx(
        joutslice.calculate_filling("Cold", broadening, eig, efermi)
    )
    with pytest.raises(NotImplementedError, match="Have not added other broadening types"):
        joutslice.calculate_filling("Unknown", broadening, eig, efermi)


def test_determine_is_metal():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    for varname in ["lumo_filling", "homo_filling", "nspin"]:
        setattr(joutslice, varname, None)
        with pytest.raises(ValueError, match=f"Cannot determine if system is metal - self.{varname} undefined"):
            joutslice.determine_is_metal()


def test_write():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    with pytest.raises(NotImplementedError):
        joutslice.write()


def test_to_dict():
    joutslice = JDFTXOutfileSlice.from_out_slice(ex_slice1)
    out_dict = joutslice.to_dict()
    assert isinstance(out_dict, dict)
