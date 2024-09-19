from pathlib import Path

import numpy as np
import pytest
from pymatgen.core.units import Ha_to_eV, bohr_to_ang
from pytest import approx

from pymatgen.io.jdftx.joutstructure import JOutStructure

ex_files_dir = Path(__file__).parents[0] / "example_files"
ex_slice_fname1 = ex_files_dir / "ex_text_slice_forJAtoms_latmin"
ex_slice1 = []
with open(ex_slice_fname1) as f:
    ex_slice1 = list.copy(list(f))
ex_slice1_known = {
    "iter": 0,
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
    "stress_00": -1.69853e-06,
    "nAtoms": 8,
    "posn0": np.array([0.000011000000000, 2.394209000000000, 1.474913000000000])
    * bohr_to_ang,
    "force0": np.array([0.000003219385226, 0.000024941936105, -0.000004667309539])
    * (Ha_to_eV / bohr_to_ang),
    "posn-1": np.array([0.000007000000000, 9.175312000000002, 4.423851000000000])
    * bohr_to_ang,
    "force-1": np.array([0.000000021330734, -0.000015026361853, -0.000010315177459])
    * (Ha_to_eV / bohr_to_ang),
    "ox0": 0.048,
    "mag0": 0.000,
    "ox-1": -0.034,
    "mag-1": 0.000,
}
ex_slice_fname2 = ex_files_dir / "ex_text_slice_forJAtoms_latmin2"
ex_slice2 = []
with open(ex_slice_fname2) as f:
    ex_slice2 = list.copy(list(f))
ex_slice2_known = {
    "iter": 9,
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
    "stress_00": -1.69853e-06,
    "nAtoms": 8,
    "posn0": np.array([0.000011000000000, 2.394209000000000, 1.474913000000000])
    * bohr_to_ang,
    "force0": np.array([0.000003219385226, 0.000024941936105, -0.000004667309539])
    * (Ha_to_eV / bohr_to_ang),
    "posn-1": np.array([0.000007000000000, 9.175312000000002, 4.423851000000000])
    * bohr_to_ang,
    "force-1": np.array([0.000000021330734, -0.000015026361853, -0.000010315177459])
    * (Ha_to_eV / bohr_to_ang),
    "ox0": 0.048,
    "mag0": 0.000,
    "ox-1": -0.034,
    "mag-1": 0.100,
}


@pytest.mark.parametrize(
    "eslice,eknowns", [(ex_slice1, ex_slice1_known), (ex_slice2, ex_slice2_known)]
)
def test_jstructure(eslice: list[str], eknowns: dict):
    jst = JOutStructure.from_text_slice(eslice, iter_type="lattice")
    assert jst.iter == eknowns["iter"]
    assert jst.etype == eknowns["etype"]
    assert approx(eknowns["E"]) == jst.E
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
    assert approx(eknowns["E0"]) == jst.elecmindata[0].E
    assert approx(eknowns["E-1"]) == jst.elecmindata[-1].E
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
