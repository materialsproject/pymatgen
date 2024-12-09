"""Tests for the JOutStructure class."""

from __future__ import annotations

from pathlib import Path

import pytest
from pytest import approx

from pymatgen.core import Structure
from pymatgen.io.jdftx.joutstructure import JOutStructure
from pymatgen.util.testing import TEST_FILES_DIR

from .outputs_test_utils import ex_jstruc_slice1 as ex_slice1
from .outputs_test_utils import ex_jstruc_slice1_known, ex_jstruc_slice2_known
from .outputs_test_utils import ex_jstruc_slice2 as ex_slice2

ex_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "example_files"


@pytest.mark.parametrize(
    ("eslice", "eknowns"), [(ex_slice1, ex_jstruc_slice1_known), (ex_slice2, ex_jstruc_slice2_known)]
)
def test_jstructure(eslice: list[str], eknowns: dict):
    """Test the JOutStructure class.

    Args:
        eslice: The text slice to test the JOutStructure class.
        eknowns: The known values to test against.
    """
    jst = JOutStructure._from_text_slice(eslice, opt_type=eknowns["opt_type"])
    assert jst.nstep == eknowns["nstep"]
    assert jst.etype == eknowns["etype"]
    assert approx(eknowns["E"]) == jst.e
    assert jst.ecomponents["Eewald"] == approx(eknowns["Eewald"])
    assert jst.ecomponents["EH"] == approx(eknowns["EH"])
    assert jst.ecomponents["Eloc"] == approx(eknowns["Eloc"])
    assert jst.ecomponents["Enl"] == approx(eknowns["Enl"])
    assert jst.ecomponents["EvdW"] == approx(eknowns["EvdW"])
    assert jst.ecomponents["Exc"] == approx(eknowns["Exc"])
    assert jst.ecomponents["Exc_core"] == approx(eknowns["Exc_core"])
    assert jst.ecomponents["KE"] == approx(eknowns["KE"])
    assert jst.ecomponents["Etot"] == approx(eknowns["Etot"])
    assert jst.ecomponents["TS"] == approx(eknowns["TS"])
    assert jst.ecomponents["F"] == approx(eknowns["F"])
    assert jst.elecmindata[0].mu == approx(eknowns["mu0"])
    assert jst.elecmindata[-1].mu == approx(eknowns["mu-1"])
    assert approx(eknowns["E0"]) == jst.elecmindata[0].e
    assert approx(eknowns["E-1"]) == jst.elecmindata[-1].e
    assert len(jst.elecmindata) == eknowns["nEminSteps"]
    assert len(jst.forces) == eknowns["nAtoms"]
    assert len(jst.cart_coords) == eknowns["nAtoms"]
    assert jst.elecmindata.converged_reason == eknowns["EconvReason"]
    assert jst.elecmindata.converged == eknowns["conv"]
    assert jst.lattice.matrix[0, 0] == approx(eknowns["cell_00"])
    assert jst.strain[0, 0] == approx(eknowns["strain_00"])
    assert jst.stress[0, 0] == approx(eknowns["stress_00"])
    for i in range(3):
        assert jst.cart_coords[0][i] == approx(eknowns["posn0"][i])
        assert jst.forces[0][i] == approx(eknowns["force0"][i])
        assert jst.cart_coords[-1][i] == approx(eknowns["posn-1"][i])
        assert jst.forces[-1][i] == approx(eknowns["force-1"][i])
    assert jst.charges[0] == approx(eknowns["ox0"])
    assert jst.magnetic_moments[0] == approx(eknowns["mag0"])
    assert jst.charges[-1] == approx(eknowns["ox-1"])
    assert jst.magnetic_moments[-1] == approx(eknowns["mag-1"])


@pytest.mark.parametrize(
    ("eslice", "eknowns"), [(ex_slice1, ex_jstruc_slice1_known), (ex_slice2, ex_jstruc_slice2_known)]
)
def test_jstructure_structure(eslice: list[str], eknowns: dict):
    jst = JOutStructure._from_text_slice(eslice, opt_type=eknowns["opt_type"])
    assert isinstance(jst.structure, Structure)
    assert not isinstance(jst.structure, JOutStructure)
