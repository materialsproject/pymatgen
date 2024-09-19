from pathlib import Path

import pytest
from pymatgen.core.units import Ha_to_eV
from pytest import approx

from atomate2.jdftx.io.joutstructures import JOutStructure, JOutStructures

ex_files_dir = Path(__file__).parents[0] / "example_files"
ex_outslice_fname1 = ex_files_dir / "ex_out_slice_latmin"
ex_outslice1 = []
with open(ex_outslice_fname1) as f:
    ex_outslice1 = list.copy(list(f))
ex_outslice1_known = {
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
ex_outslice_fname2 = ex_files_dir / "ex_out_slice_ionmin"
with open(ex_outslice_fname2) as f:
    ex_outslice2 = list.copy(list(f))
ex_outslice2_known = {
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


@pytest.mark.parametrize(
    "ex_slice, ex_slice_known,iter_type",
    [(ex_outslice1, ex_outslice1_known, "lattice")],
)
def test_jstructures(
    ex_slice: list[str], ex_slice_known: dict[str, float], iter_type: str
):
    jstruct = JOutStructures.from_out_slice(ex_slice, iter_type=iter_type)
    assert isinstance(jstruct, JOutStructures)
    assert isinstance(jstruct[0], JOutStructure)
    assert jstruct[0].elecmindata[0].mu == approx(ex_slice_known["mu0_0"])
    assert jstruct[0].elecmindata[-1].mu == approx(ex_slice_known["mu0_-1"])
    assert jstruct[-1].elecmindata[0].mu == approx(ex_slice_known["mu-1_0"])
    assert jstruct[-1].elecmindata[-1].mu == approx(ex_slice_known["mu-1_-1"])
    assert jstruct.elecmindata[-1].mu == approx(ex_slice_known["mu-1_-1"])
    assert jstruct[0].elecmindata[0].nelectrons == approx(ex_slice_known["nelec0_0"])
    assert jstruct[0].elecmindata[-1].nelectrons == approx(ex_slice_known["nelec0_-1"])
    assert jstruct[-1].elecmindata[0].nelectrons == approx(ex_slice_known["nelec-1_0"])
    assert jstruct[-1].elecmindata[-1].nelectrons == approx(
        ex_slice_known["nelec-1_-1"]
    )
    assert jstruct.elecmindata[-1].nelectrons == approx(ex_slice_known["nelec-1_-1"])
    assert len(jstruct[0].elecmindata) == ex_slice_known["nEminSteps0"]
    assert len(jstruct[-1].elecmindata) == ex_slice_known["nEminSteps-1"]
    assert len(jstruct.elecmindata) == ex_slice_known["nEminSteps-1"]
    assert jstruct[0].etype == ex_slice_known["etype0"]
    assert approx(ex_slice_known["E0"]) == jstruct[0].E
    assert jstruct[-1].etype == ex_slice_known["etype-1"]
    assert jstruct.etype == ex_slice_known["etype-1"]
    assert approx(ex_slice_known["E-1"]) == jstruct[-1].E
    assert jstruct[0].elecmindata.converged == ex_slice_known["conv0"]
    assert jstruct[-1].elecmindata.converged == ex_slice_known["conv-1"]
    assert jstruct.elecmindata.converged == ex_slice_known["conv-1"]
    assert len(jstruct) == ex_slice_known["nGeomSteps"]
