# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import unittest
from pathlib import Path

from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cp2k.sets import (
    SETTINGS,
    Cp2kValidationError,
    DftSet,
    GaussianTypeOrbitalBasisSet,
    GthPotential,
)
from pymatgen.util.testing import PymatgenTest

Si_structure = Structure(
    lattice=[[0, 2.734364, 2.734364], [2.734364, 0, 2.734364], [2.734364, 2.734364, 0]],
    species=["Si", "Si"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)
molecule = Molecule(species=["Si"], coords=[[0, 0, 0]])


class SetTest(PymatgenTest):
    def setUp(self) -> None:
        self.TEST_FILES_DIR = Path.joinpath(self.TEST_FILES_DIR, "cp2k")
        SETTINGS["PMG_CP2K_DATA_DIR"] = self.TEST_FILES_DIR
        self.setkwargs = {
            "print_pdos": False,
            "print_dos": False,
            "print_v_hartree": False,
            "print_e_density": False,
        }
        return super().setUp()

    def test_dft_set(self):

        # Basis sets / potentials searching
        basis_and_potential = {"basis_type": "SZV", "potential_type": "Pseudopotential", "functional": None}
        ss = DftSet(Si_structure, basis_and_potential=basis_and_potential, xc_functionals="PBE")

        # Basis sets / potentials by name
        basis_and_potential = {"Si": {"basis": "SZV-GTH-q4", "potential": "GTH-PBE-q4"}}
        ss = DftSet(Si_structure, basis_and_potential=basis_and_potential, xc_functionals="PBE")

        # Basis sets / potentials by name with ADMM
        basis_and_potential = {"Si": {"basis": "SZV-GTH-q4", "potential": "GTH-PBE-q4", "aux_basis": "cFIT3"}}
        ss = DftSet(Si_structure, basis_and_potential=basis_and_potential, xc_functionals="PBE")
        basis_sets = ss["force_eval"]["subsys"]["Si_1"].get("basis_set")
        self.assertTrue(any("AUX_FIT" in b.values for b in basis_sets))
        self.assertTrue(any("cFIT3" in b.values for b in basis_sets))

        # Basis sets / potentials by hash value
        basis_and_potential = {
            "Si": {"basis": "30767c18f6e7e46c1b56c1d34ff6007d", "potential": "21e2f468a18404ff6119fe801da81e43"}
        }
        ss = DftSet(Si_structure, basis_and_potential=basis_and_potential, xc_functionals="PBE")

        # Basis set / potential with objects
        gto = """
         Si SZV-MOLOPT-GTH SZV-MOLOPT-GTH-q4
            1
            2 0 1 6 1 1
                2.693604434572  0.015333179500 -0.005800105400
                1.359613855428 -0.283798205000 -0.059172026000
                0.513245176029 -0.228939692700  0.121487149900
                0.326563011394  0.728834000900  0.423382421100
                0.139986977410  0.446205299300  0.474592116300
                0.068212286977  0.122025292800  0.250129397700
        """
        pot = """Si GTH-BLYP-q4 GTH-BLYP
                2    2
                0.44000000    1    -6.25958674
                2
                0.44465247    2     8.31460936    -2.33277947
                                                    3.01160535
                0.50279207    1     2.33241791"""
        basis_and_potential = {
            "Si": {"basis": GaussianTypeOrbitalBasisSet.from_string(gto), "potential": GthPotential.from_string(pot)}
        }
        ss = DftSet(Si_structure, basis_and_potential=basis_and_potential, xc_functionals="PBE", **self.setkwargs)
        self.assertAlmostEqual(ss.cutoff, 150)

        # Test that printing will activate sections
        self.assertFalse(ss.check("motion"))
        ss.activate_motion()
        self.assertTrue(ss.check("motion"))
        self.assertFalse(ss.check("force_eval/dft/print/pdos"))
        ss.print_pdos()
        self.assertTrue(ss.check("force_eval/dft/print/pdos"))
        self.assertFalse(ss.check("force_eval/dft/print/v_hartree_cube"))
        ss.print_v_hartree()
        self.assertTrue(ss.check("force_eval/dft/print/v_hartree_cube"))

        # Test property activators
        self.assertFalse(ss.check("force_eval/properties"))
        ss.activate_nmr()
        ss.activate_epr()
        ss.activate_hyperfine()
        ss.activate_polar()
        ss.activate_tddfpt()
        self.assertTrue(ss.check("force_eval/properties/linres/localize"))
        self.assertTrue(ss.check("force_eval/properties/linres/nmr/print/chi_tensor"))
        self.assertTrue(ss.check("force_eval/properties/linres/epr/print/g_tensor"))
        self.assertTrue(ss.check("force_eval/properties/tddfpt"))
        self.assertTrue(ss.check("force_eval/dft/print/hyperfine_coupling_tensor"))

        # For at least up to v2022.1, DOS doesn't work without kpoints
        self.assertFalse(ss.check("force_eval/dft/print/dos"))
        ss.print_dos()
        self.assertFalse(ss.check("force_eval/dft/print/dos"))

        self.assertFalse(ss.check("force_eval/dft/xc/hf"))
        ss.activate_hybrid()
        self.assertTrue(ss.check("force_eval/dft/xc/hf"))
        self.assertTrue(ss.check("force_eval/dft/auxiliary_density_matrix_method"))

        # Validator will trip for kpoints + hfx
        ss.update({"force_eval": {"dft": {"kpoints": {}}}})
        with self.assertRaises(Cp2kValidationError):
            ss.validate()

        ss = DftSet(molecule, basis_and_potential=basis_and_potential, xc_functionals="PBE")
        self.assertTrue(ss.check("force_eval/dft/poisson"))
        self.assertEqual(ss["force_eval"]["dft"]["poisson"].get("periodic").values[0].upper(), "NONE")


if __name__ == "__main__":
    unittest.main()
