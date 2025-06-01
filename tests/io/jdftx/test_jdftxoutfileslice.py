from __future__ import annotations

import math

import numpy as np
import pytest

from pymatgen.core.trajectory import Trajectory
from pymatgen.core.units import Ha_to_eV, ang_to_bohr
from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice

from .outputs_test_utils import ex_outfileslice1 as ex_slice1


def test_jdftxoutfileslice_stringify():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    out_str = str(joutslice)
    assert isinstance(out_str, str)
    assert out_str


def test_jdftxoutfileslice_converge():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    assert joutslice.converged


def test_jdftxoutfileslice_trajectory():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    traj = joutslice.trajectory
    assert isinstance(traj, Trajectory)
    assert len(traj) == len(joutslice.jstrucs)


def test_get_broadeningvars():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    btype = "btype"
    bval = 1.0
    text = [f"elec-smearing {btype} {bval}"]
    broadening_type, broadening = joutslice._get_broadeningvars(text)
    assert broadening_type == btype
    assert broadening == pytest.approx(bval)
    broadening_type, broadening = joutslice._get_broadeningvars([])
    assert broadening_type is None
    assert broadening == pytest.approx(0.0)


def test_get_truncationvars():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    joutslice.is_bgw = True
    with pytest.raises(ValueError, match="BGW slab Coulomb truncation must be along z!"):
        joutslice._get_truncationvars(["coulomb-interaction Slab 010"])
    with pytest.raises(ValueError, match="BGW wire Coulomb truncation must be periodic in z!"):
        joutslice._get_truncationvars(["coulomb-interaction Cylindrical 010"])
    with pytest.raises(ValueError, match="Problem with this truncation!"):
        joutslice._get_truncationvars(["coulomb-interaction barbie 010"])
    truncation_type, truncation_radius = joutslice._get_truncationvars(
        ["coulomb-interaction Spherical", "Initialized spherical truncation of radius 1.0"]
    )
    assert truncation_type == "spherical"
    assert truncation_radius == pytest.approx(1.0 / ang_to_bohr)


def test_get_rho_cutoff():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    text = ["elec-cutoff 1.0"]
    joutslice.pwcut = None
    rhocut = joutslice._get_rho_cutoff(text)
    assert joutslice.pwcut == pytest.approx(1.0 * Ha_to_eV)
    assert rhocut == pytest.approx(joutslice.pwcut * 4)


def test_get_eigstats_varsdict():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    evardict = joutslice._get_eigstats_varsdict([], "$VAR")
    for key in evardict:
        assert evardict[key] is None
    # Initializing eigvars with no data will set all to None EXCEPT efermi which has "mu" as a backup reference
    joutslice._set_eigvars([])
    for key in evardict:
        if key != "efermi":
            assert getattr(joutslice, key) is None
    # Setting mu, _mu_backup, and jstrucs to None before _set_eigvars([]) will set efermi to None
    joutslice.mu = None
    joutslice._mu_backup = None
    joutslice.jstrucs = None
    joutslice._set_eigvars([])
    for key in evardict:
        assert getattr(joutslice, key) is None


def test_get_pp_type():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    assert joutslice._get_pp_type(["Reading pseudopotential file root/PAW:"]) is None
    assert joutslice._get_pp_type(["Reading pseudopotential file not_SG15/GBRV"]) == "GBRV"
    assert joutslice._get_pp_type(["Reading pseudopotential file not_GBRV/SG15"]) == "SG15"


def test_set_pseudo_vars_t1():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    # 6 instances of "reading pseudopotential file" since all possible test files have less than 6 atom types
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
    joutslice.total_electrons = None
    with pytest.raises(ValueError, match="Total electrons and semicore electrons must be set."):
        joutslice._set_pseudo_vars_t1(text)
    joutslice.atom_elements = None
    with pytest.raises(ValueError, match="Atom elements not set yet."):
        joutslice._set_pseudo_vars_t1(text)
    joutslice.atom_types = None
    with pytest.raises(ValueError, match="Pseudopotential data cannot be allocated without atom types."):
        joutslice._set_pseudo_vars_t1(text)


def test_set_orb_fillings_nobroad():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    joutslice._set_orb_fillings_nobroad(1)
    assert joutslice.homo_filling == pytest.approx(2)
    assert joutslice.lumo_filling == pytest.approx(0)


def test_set_orb_fillings_broad():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    joutslice.lumo = None
    with pytest.raises(ValueError, match="Cannot set orbital fillings with broadening with self.lumo as None"):
        joutslice._set_orb_fillings()
    joutslice.homo = None
    with pytest.raises(ValueError, match="Cannot set orbital fillings with broadening with self.homo as None"):
        joutslice._set_orb_fillings()
    joutslice.efermi = None
    with pytest.raises(ValueError, match="Cannot set orbital fillings with broadening with self.efermi as None"):
        joutslice._set_orb_fillings()
    joutslice.broadening = None
    with pytest.raises(ValueError, match="Cannot set orbital fillings with broadening with self.broadening as None"):
        joutslice._set_orb_fillings()
    joutslice.broadening_type = None
    joutslice.nspin = None
    with pytest.raises(ValueError, match="Cannot set homo/lumo filling with self.nspin as None"):
        joutslice._set_orb_fillings()


def test_set_lattice_vars():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    joutslice.jstrucs = None
    with pytest.raises(ValueError, match="No structures found in out file."):
        joutslice._set_lattice_vars([])


def test_set_ecomponents():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    joutslice.jstrucs = None
    with pytest.raises(ValueError, match="No structures found in out file."):
        joutslice._set_ecomponents([])


def test_calculate_filling():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    broadening = 1.0
    eig = 0.5
    efermi = 0.6
    x = (eig - efermi) / (2.0 * broadening)
    assert 0.5 * (1 - np.tanh(x)) == pytest.approx(joutslice._calculate_filling("Fermi", broadening, eig, efermi))
    assert 0.5 * (1 - math.erf(x)) - x * np.exp(-1 * x**2) / (2 * np.pi**0.5) == pytest.approx(
        joutslice._calculate_filling("MP1", broadening, eig, efermi)
    )
    assert 0.5 * (1 - math.erf(x)) == pytest.approx(joutslice._calculate_filling("Gauss", broadening, eig, efermi))
    assert 0.5 * (1 - math.erf(x + 0.5**0.5)) + np.exp(-1 * (x + 0.5**0.5) ** 2) / (2 * np.pi) ** 0.5 == pytest.approx(
        joutslice._calculate_filling("Cold", broadening, eig, efermi)
    )
    with pytest.raises(NotImplementedError, match="Have not added other broadening types"):
        joutslice._calculate_filling("Unknown", broadening, eig, efermi)


def test_determine_is_metal():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    for varname in ["lumo_filling", "homo_filling", "nspin"]:
        setattr(joutslice, varname, None)
        with pytest.raises(ValueError, match=f"Cannot determine if system is metal - self.{varname} undefined"):
            joutslice.determine_is_metal()


def test_write():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    with pytest.raises(NotImplementedError):
        joutslice.write()


def test_as_dict():
    joutslice = JDFTXOutfileSlice._from_out_slice(ex_slice1)
    out_dict = joutslice.as_dict()
    assert isinstance(out_dict, dict)


def should_be_parsable_out_slice(out_slice: list[str]):
    return any("ElecMinimize: Iter:" in line for line in out_slice[::-1])


# Make sure all possible exceptions are caught when none_on_error is True
@pytest.mark.parametrize(("ex_slice"), [(ex_slice1)])
def test_none_on_partial(ex_slice: list[str]):
    # freq = 1 takes about 5 seconds so cutting down the number of tests is needed
    freq = 5
    for i in range(int(len(ex_slice) / freq)):
        test_slice = ex_slice[: -(i * freq)]
        joutslice = JDFTXOutfileSlice._from_out_slice(test_slice, none_on_error=True)
        if should_be_parsable_out_slice(test_slice):
            assert isinstance(joutslice, JDFTXOutfileSlice | None)
        else:
            assert isinstance(joutslice, JDFTXOutfileSlice | None)
