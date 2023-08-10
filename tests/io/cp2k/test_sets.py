from __future__ import annotations

import pytest
from pytest import approx

from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cp2k.sets import SETTINGS, Cp2kValidationError, DftSet, GaussianTypeOrbitalBasisSet, GthPotential
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

Si_structure = Structure(
    lattice=[[0, 2.734364, 2.734364], [2.734364, 0, 2.734364], [2.734364, 2.734364, 0]],
    species=["Si", "Si"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)
molecule = Molecule(species=["Si"], coords=[[0, 0, 0]])


class TestSet(PymatgenTest):
    def setUp(self) -> None:
        SETTINGS["PMG_CP2K_DATA_DIR"] = f"{TEST_FILES_DIR}/cp2k"
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
        assert any("AUX_FIT" in b.values for b in basis_sets)  # noqa: PD011
        assert any("cFIT3" in b.values for b in basis_sets)  # noqa: PD011

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
            "Si": {"basis": GaussianTypeOrbitalBasisSet.from_str(gto), "potential": GthPotential.from_str(pot)}
        }
        ss = DftSet(Si_structure, basis_and_potential=basis_and_potential, xc_functionals="PBE", **self.setkwargs)
        assert ss.cutoff == approx(150)

        # Test that printing will activate sections
        assert not ss.check("motion")
        ss.activate_motion()
        assert ss.check("motion")
        assert not ss.check("force_eval/dft/print/pdos")
        ss.print_pdos()
        assert ss.check("force_eval/dft/print/pdos")
        assert not ss.check("force_eval/dft/print/v_hartree_cube")
        ss.print_v_hartree()
        assert ss.check("force_eval/dft/print/v_hartree_cube")

        # Test property activators
        assert not ss.check("force_eval/properties")
        ss.activate_nmr()
        ss.activate_epr()
        ss.activate_hyperfine()
        ss.activate_polar()
        ss.activate_tddfpt()
        assert ss.check("force_eval/properties/linres/localize")
        assert ss.check("force_eval/properties/linres/nmr/print/chi_tensor")
        assert ss.check("force_eval/properties/linres/epr/print/g_tensor")
        assert ss.check("force_eval/properties/tddfpt")
        assert ss.check("force_eval/dft/print/hyperfine_coupling_tensor")

        # For at least up to v2022.1, DOS doesn't work without kpoints
        assert not ss.check("force_eval/dft/print/dos")
        ss.print_dos()
        assert not ss.check("force_eval/dft/print/dos")

        assert not ss.check("force_eval/dft/xc/hf")
        ss.activate_hybrid()
        assert ss.check("force_eval/dft/xc/hf")
        assert ss.check("force_eval/dft/auxiliary_density_matrix_method")

        # Validator will trip for kpoints + hfx
        ss.update({"force_eval": {"dft": {"kpoints": {}}}})
        with pytest.raises(Cp2kValidationError, match="CP2K v2022.1: Does not support hartree fock with kpoints"):
            ss.validate()

        ss = DftSet(molecule, basis_and_potential=basis_and_potential, xc_functionals="PBE")
        assert ss.check("force_eval/dft/poisson")
        assert ss["force_eval"]["dft"]["poisson"].get("periodic").values[0].upper() == "NONE"  # noqa: PD011
