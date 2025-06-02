"""Tests for the JOutStructures class."""

from __future__ import annotations

import pytest
from pytest import approx

from pymatgen.io.jdftx.joutstructures import JOutStructure, JOutStructures

from .outputs_test_utils import ex_outfileslice1 as ex_outslice1
from .outputs_test_utils import ex_outfileslice1_known as ex_outslice1_known


@pytest.mark.parametrize(
    ("ex_slice", "ex_slice_known", "opt_type"),
    [(ex_outslice1, ex_outslice1_known, "lattice")],
)
def test_jstructures(ex_slice: list[str], ex_slice_known: dict[str, float], opt_type: str):
    """General testing of the JOutStructures class.

    Args:
        ex_slice: The example slice to test.
        ex_slice_known: The known values of the example slice.
        opt_type: The optimization type
    """
    jstruct = JOutStructures._from_out_slice(ex_slice, opt_type=opt_type)
    assert isinstance(jstruct, JOutStructures)
    # Testing indexing
    assert isinstance(jstruct[0], JOutStructure)
    # Indexing should be equivalent to indexing the slices attribute
    assert jstruct[0] is jstruct.slices[0]
    # Testing initialization matches known values stored in ex_slice_known
    assert jstruct[0].elecmindata[0].mu == approx(ex_slice_known["mu0_0"])
    assert jstruct[0].elecmindata[-1].mu == approx(ex_slice_known["mu0_-1"])
    assert jstruct[-1].elecmindata[0].mu == approx(ex_slice_known["mu-1_0"])
    assert jstruct[-1].elecmindata[-1].mu == approx(ex_slice_known["mu-1_-1"])
    assert jstruct.elecmindata[-1].mu == approx(ex_slice_known["mu-1_-1"])
    assert jstruct.elecmindata.mu == approx(ex_slice_known["mu-1_-1"])
    assert jstruct.mu == approx(ex_slice_known["mu-1_-1"])
    assert jstruct["mu"] == approx(ex_slice_known["mu-1_-1"])
    assert jstruct[0].elecmindata[0].nelectrons == approx(ex_slice_known["nelec0_0"])
    assert jstruct[0].elecmindata[-1].nelectrons == approx(ex_slice_known["nelec0_-1"])
    assert jstruct[-1].elecmindata[0].nelectrons == approx(ex_slice_known["nelec-1_0"])
    assert jstruct[-1].elecmindata[-1].nelectrons == approx(ex_slice_known["nelec-1_-1"])
    assert jstruct.elecmindata[-1].nelectrons == approx(ex_slice_known["nelec-1_-1"])
    assert len(jstruct[0].elecmindata) == ex_slice_known["nEminSteps0"]
    assert len(jstruct[-1].elecmindata) == ex_slice_known["nEminSteps-1"]
    assert len(jstruct.elecmindata) == ex_slice_known["nEminSteps-1"]
    assert jstruct[0].etype == ex_slice_known["etype0"]
    assert approx(ex_slice_known["E0"]) == jstruct[0].e
    assert jstruct[-1].etype == ex_slice_known["etype-1"]
    assert jstruct.etype == ex_slice_known["etype-1"]
    assert approx(ex_slice_known["E-1"]) == jstruct[-1].e
    assert jstruct[0].elecmindata.converged == ex_slice_known["conv0"]
    assert jstruct[-1].elecmindata.converged == ex_slice_known["conv-1"]
    assert jstruct.elecmindata.converged == ex_slice_known["conv-1"]
    assert len(jstruct) == ex_slice_known["nGeomSteps"]
    # Lazy testing inherited attribute
    assert isinstance(jstruct.selective_dynamics, list)
