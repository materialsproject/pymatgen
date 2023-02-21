# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import os
import unittest

from monty.serialization import dumpfn, loadfn

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.outputs import QCOutput, check_for_structure_changes
from pymatgen.util.testing import PymatgenTest

try:
    from openbabel import openbabel
except ImportError:
    openbabel = None

__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath, Evan Spotte-Smith, Ryan Kingsbury"
__copyright__ = "Copyright 2018-2022, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

single_job_dict = loadfn(os.path.join(os.path.dirname(__file__), "single_job.json"))
multi_job_dict = loadfn(os.path.join(os.path.dirname(__file__), "multi_job.json"))

property_list = {
    "errors",
    "multiple_outputs",
    "completion",
    "unrestricted",
    "using_GEN_SCFMAN",
    "final_energy",
    "S2",
    "optimization",
    "energy_trajectory",
    "opt_constraint",
    "frequency_job",
    "charge",
    "multiplicity",
    "species",
    "initial_geometry",
    "initial_molecule",
    "SCF",
    "Mulliken",
    "optimized_geometry",
    "optimized_zmat",
    "molecule_from_optimized_geometry",
    "last_geometry",
    "molecule_from_last_geometry",
    "geometries",
    "gradients",
    "frequency_mode_vectors",
    "walltime",
    "cputime",
    "point_group",
    "frequencies",
    "IR_intens",
    "IR_active",
    "g_electrostatic",
    "g_cavitation",
    "g_dispersion",
    "g_repulsion",
    "total_contribution_pcm",
    "ZPE",
    "trans_enthalpy",
    "vib_enthalpy",
    "rot_enthalpy",
    "gas_constant",
    "trans_entropy",
    "vib_entropy",
    "rot_entropy",
    "total_entropy",
    "total_enthalpy",
    "warnings",
    "SCF_energy_in_the_final_basis_set",
    "Total_energy_in_the_final_basis_set",
    "solvent_method",
    "solvent_data",
    "using_dft_d3",
    "single_point_job",
    "force_job",
    "pcm_gradients",
    "CDS_gradients",
    "RESP",
    "trans_dip",
    "transition_state",
    "scan_job",
    "optimized_geometries",
    "molecules_from_optimized_geometries",
    "scan_energies",
    "scan_constraint_sets",
    "hf_scf_energy",
    "mp2_energy",
    "ccsd_correlation_energy",
    "ccsd_total_energy",
    "ccsd(t)_correlation_energy",
    "ccsd(t)_total_energy",
    "cdft_becke_excess_electrons",
    "cdft_becke_population",
    "cdft_becke_net_spin",
    "direct_coupling_Hif_Hartree",
    "direct_coupling_Sif_Hartree",
    "direct_coupling_Hii_Hartree",
    "direct_coupling_Sii_Hartree",
    "direct_coupling_Hff_Hartree",
    "direct_coupling_Sff_Hartree",
    "direct_coupling_eV",
    "almo_coupling_states",
    "almo_diabat_energies_Hartree",
    "almo_adiabat_energies_Hartree",
    "almo_hamiltonian",
    "almo_overlap_matrix",
    "almo_s2_matrix",
    "almo_diabat_basis_coeff",
    "almo_h_coupling_matrix",
    "almo_coupling_eV",
    "pod_coupling_eV",
    "fodft_had_eV",
    "fodft_hda_eV",
    "fodft_coupling_eV",
    "alpha_fock_matrix",
    "beta_fock_matrix",
    "alpha_eigenvalues",
    "beta_eigenvalues",
    "alpha_coeff_matrix",
    "beta_coeff_matrix",
    "final_soln_phase_e",
    "solute_internal_e",
    "total_solvation_free_e",
    "change_solute_internal_e",
    "reaction_field_free_e",
    "isosvp_dielectric",
    "dispersion_e",
    "exchange_e",
    "min_neg_field_e",
    "max_pos_field_e",
    "norm_of_stepsize",
    "version",
    "dipoles",
    "gap_info",
}

if openbabel is not None:
    property_list.add("structure_change")

single_job_out_names = {
    "unable_to_determine_lambda_in_geom_opt.qcout",
    "thiophene_wfs_5_carboxyl.qcout",
    "hf.qcout",
    "hf_opt_failed.qcout",
    "no_reading.qcout",
    "exit_code_134.qcout",
    "negative_eigen.qcout",
    "insufficient_memory.qcout",
    "freq_seg_too_small.qcout",
    "crowd_gradient_number.qcout",
    "quinoxaline_anion.qcout",
    "tfsi_nbo.qcout",
    "crowd_nbo_charges.qcout",
    "h2o_aimd.qcout",
    "quinoxaline_anion.qcout",
    "crowd_gradient_number.qcout",
    "bsse.qcout",
    "thiophene_wfs_5_carboxyl.qcout",
    "time_nan_values.qcout",
    "pt_dft_180.0.qcout",
    "qchem_energies/hf-rimp2.qcout",
    "qchem_energies/hf_b3lyp.qcout",
    "qchem_energies/hf_ccsd(t).qcout",
    "qchem_energies/hf_cosmo.qcout",
    "qchem_energies/hf_hf.qcout",
    "qchem_energies/hf_lxygjos.qcout",
    "qchem_energies/hf_mosmp2.qcout",
    "qchem_energies/hf_mp2.qcout",
    "qchem_energies/hf_qcisd(t).qcout",
    "qchem_energies/hf_riccsd(t).qcout",
    "qchem_energies/hf_tpssh.qcout",
    "qchem_energies/hf_xyg3.qcout",
    "qchem_energies/hf_xygjos.qcout",
    "qchem_energies/hf_wb97xd_gen_scfman.qcout",
    "new_qchem_files/pt_n2_n_wb_180.0.qcout",
    "new_qchem_files/pt_n2_trip_wb_90.0.qcout",
    "new_qchem_files/pt_n2_gs_rimp2_pvqz_90.0.qcout",
    "new_qchem_files/VC_solv_eps10.2.qcout",
    "crazy_scf_values.qcout",
    "new_qchem_files/N2.qcout",
    "new_qchem_files/julian.qcout.gz",
    "new_qchem_files/Frequency_no_equal.qout",
    "new_qchem_files/gdm.qout",
    "new_qchem_files/DinfH.qout",
    "new_qchem_files/mpi_error.qout",
    "new_qchem_files/molecule_read_error.qout",
    "new_qchem_files/basis_not_supported.qout",
    "new_qchem_files/lebdevpts.qout",
    "new_qchem_files/Optimization_no_equal.qout",
    "new_qchem_files/2068.qout",
    "new_qchem_files/2620.qout",
    "new_qchem_files/1746.qout",
    "new_qchem_files/1570.qout",
    "new_qchem_files/1570_2.qout",
    "new_qchem_files/single_point.qout",
    "new_qchem_files/roothaan_diis_gdm.qout",
    "new_qchem_files/pes_scan_single_variable.qout",
    "new_qchem_files/pes_scan_double_variable.qout",
    "new_qchem_files/ts.out",
    "new_qchem_files/ccsd.qout",
    "new_qchem_files/ccsdt.qout",
    "new_qchem_files/almo.out",
    "new_qchem_files/cdft_simple.qout",
    "new_qchem_files/fodft.out",
    "new_qchem_files/fodft_2.out",
    "new_qchem_files/fodft_3.out",
    "new_qchem_files/pod1.out",
    "new_qchem_files/pod2_gs.out",
    "extra_scf_print.qcout",
    "new_qchem_files/cmirs_benzene_single.qcout",
    "new_qchem_files/cmirs_dielst10_single.qcout",
    "new_qchem_files/cmirs_water_single.qcout",
    "new_qchem_files/isosvp_water_single.qcout",
    "new_qchem_files/isosvp_dielst10_single.qcout",
    "new_qchem_files/custom_gdm_gdmqls_opt.qout",
    "new_qchem_files/unable.qout",
    "new_qchem_files/unexpected_ts.out",
    "new_qchem_files/svd_failed.qout",
    "new_qchem_files/v6_old_driver.out",
    "new_qchem_files/gap.qout",
    "new_qchem_files/3C.qout",
    "new_qchem_files/hyper.qout",
    "new_qchem_files/os_gap.qout",
}

multi_job_out_names = {
    "not_enough_total_memory.qcout",
    "new_qchem_files/VC_solv_eps10.qcout",
    "new_qchem_files/MECLi_solv_eps10.qcout",
    "pcm_solvent_deprecated.qcout",
    "qchem43_batch_job.qcout",
    "ferrocenium_1pos.qcout",
    "CdBr2.qcout",
    "killed.qcout",
    "aux_mpi_time_mol.qcout",
    "new_qchem_files/VCLi_solv_eps10.qcout",
    "new_qchem_files/cdft_dc.qout",
    "new_qchem_files/cmirs_benzene.qcout",
    "new_qchem_files/cmirs_dielst10.qcout",
    "new_qchem_files/isosvp_water.qcout",
    "new_qchem_files/isosvp_dielst10.qcout",
}


class TestQCOutput(PymatgenTest):
    @staticmethod
    def generate_single_job_dict():
        """
        Used to generate test dictionary for single jobs.
        """
        single_job_dict = {}
        for file in single_job_out_names:
            single_job_dict[file] = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", file)).data
        dumpfn(single_job_dict, "single_job.json")

    @staticmethod
    def generate_multi_job_dict():
        """
        Used to generate test dictionary for multiple jobs.
        """
        multi_job_dict = {}
        for file in multi_job_out_names:
            outputs = QCOutput.multiple_outputs_from_file(
                os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", file), keep_sub_files=False
            )
            data = []
            for sub_output in outputs:
                data.append(sub_output.data)
            multi_job_dict[file] = data
        dumpfn(multi_job_dict, "multi_job.json")

    def _test_property(self, key, single_outs, multi_outs):
        for name, outdata in single_outs.items():
            try:
                assert outdata.get(key) == single_job_dict[name].get(key)
            except ValueError:
                try:
                    self.assertArrayEqual(outdata.get(key), single_job_dict[name].get(key))
                except AssertionError:
                    raise RuntimeError("Issue with file: " + name + " Exiting...")
            except AssertionError:
                raise RuntimeError("Issue with file: " + name + " Exiting...")
        for name, outputs in multi_outs.items():
            for ii, sub_output in enumerate(outputs):
                try:
                    assert sub_output.data.get(key) == multi_job_dict[name][ii].get(key)
                except ValueError:
                    self.assertArrayEqual(sub_output.data.get(key), multi_job_dict[name][ii].get(key))

    @unittest.skipIf(openbabel is None, "OpenBabel not installed.")
    def test_all(self):
        self.maxDiff = None
        single_outs = {}
        for file in single_job_out_names:
            single_outs[file] = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", file)).data

        multi_outs = {}
        for file in multi_job_out_names:
            multi_outs[file] = QCOutput.multiple_outputs_from_file(
                os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", file), keep_sub_files=False
            )

        for key in property_list:
            print("Testing ", key)
            self._test_property(key, single_outs, multi_outs)

    @unittest.skipIf((openbabel is None), "OpenBabel not installed.")
    def test_structural_change(self):
        t1 = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "structural_change", "t1.xyz"))
        t2 = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "structural_change", "t2.xyz"))
        t3 = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "structural_change", "t3.xyz"))

        thio_1 = Molecule.from_file(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "structural_change", "thiophene1.xyz")
        )
        thio_2 = Molecule.from_file(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "structural_change", "thiophene2.xyz")
        )

        frag_1 = Molecule.from_file(
            os.path.join(
                PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "test_structure_change", "frag_1.xyz"
            )
        )
        frag_2 = Molecule.from_file(
            os.path.join(
                PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "test_structure_change", "frag_2.xyz"
            )
        )

        assert check_for_structure_changes(t1, t1) == "no_change"
        assert check_for_structure_changes(t2, t3) == "no_change"
        assert check_for_structure_changes(t1, t2) == "fewer_bonds"
        assert check_for_structure_changes(t2, t1) == "more_bonds"

        assert check_for_structure_changes(thio_1, thio_2) == "unconnected_fragments"

        assert check_for_structure_changes(frag_1, frag_2) == "bond_change"

    def test_NBO_parsing(self):
        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "nbo.qout")).data
        assert len(data["nbo_data"]["natural_populations"]) == 3
        assert len(data["nbo_data"]["hybridization_character"]) == 6
        assert len(data["nbo_data"]["perturbation_energy"]) == 2
        assert data["nbo_data"]["natural_populations"][0]["Density"][5] == -0.08624
        assert data["nbo_data"]["hybridization_character"][4]["atom 2 pol coeff"][35] == "-0.7059"
        next_to_last = list(data["nbo_data"]["perturbation_energy"][-1]["fock matrix element"])[-2]
        assert data["nbo_data"]["perturbation_energy"][-1]["fock matrix element"][next_to_last] == 0.071
        assert data["nbo_data"]["perturbation_energy"][0]["acceptor type"][0] == "RY*"

    def test_NBO7_parsing(self):
        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "nbo7_1.qout")).data
        assert data["nbo_data"]["perturbation_energy"][0]["perturbation energy"][9] == 15.73
        assert len(data["nbo_data"]["perturbation_energy"][0]["donor bond index"]) == 84
        assert len(data["nbo_data"]["perturbation_energy"][1]["donor bond index"]) == 29

        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "nbo7_2.qout")).data
        assert data["nbo_data"]["perturbation_energy"][0]["perturbation energy"][13] == 32.93
        assert data["nbo_data"]["perturbation_energy"][0]["acceptor type"][13] == "LV"
        assert data["nbo_data"]["perturbation_energy"][0]["acceptor type"][12] == "RY"
        assert data["nbo_data"]["perturbation_energy"][0]["acceptor atom 1 symbol"][12] == "Mg"

        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "nbo7_3.qout")).data
        assert data["nbo_data"]["perturbation_energy"][0]["perturbation energy"][13] == 34.54
        assert data["nbo_data"]["perturbation_energy"][0]["acceptor type"][13] == "BD*"
        assert data["nbo_data"]["perturbation_energy"][0]["acceptor atom 1 symbol"][13] == "B"
        assert data["nbo_data"]["perturbation_energy"][0]["acceptor atom 2 symbol"][13] == "Mg"
        assert data["nbo_data"]["perturbation_energy"][0]["acceptor atom 2 number"][13] == 3

    def test_NBO5_vs_NBO7_hybridization_character(self):
        data5 = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "nbo5_1.qout")).data
        data7 = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "nbo7_1.qout")).data
        assert len(data5["nbo_data"]["hybridization_character"]) == len(data7["nbo_data"]["hybridization_character"])
        assert (
            data5["nbo_data"]["hybridization_character"][4]["atom 2 pol coeff"][9]
            == data7["nbo_data"]["hybridization_character"][4]["atom 2 pol coeff"][9]
        )
        assert (
            data5["nbo_data"]["hybridization_character"][0]["s"][0]
            == data7["nbo_data"]["hybridization_character"][0]["s"][0]
        )
        assert data5["nbo_data"]["hybridization_character"][1]["bond index"][7] == "149"
        assert data7["nbo_data"]["hybridization_character"][1]["bond index"][7] == "21"

    def test_NBO7_infinite_e2pert(self):
        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "nbo7_inf.qout")).data
        assert data["nbo_data"]["perturbation_energy"][0]["perturbation energy"][0] == float("inf")

    def test_cdft_parsing(self):
        data = QCOutput(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "cdft_simple.qout")
        ).data
        assert data["cdft_becke_excess_electrons"][0][0] == 0.432641
        assert len(data["cdft_becke_population"][0]) == 12
        assert data["cdft_becke_net_spin"][0][6] == -0.000316

    def test_cdft_dc_parsing(self):
        data = QCOutput.multiple_outputs_from_file(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "cdft_dc.qout"),
            keep_sub_files=False,
        )[-1].data
        assert data["direct_coupling_eV"] == 0.0103038246

    def test_almo_msdft2_parsing(self):
        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "almo.out")).data
        self.assertListEqual(data["almo_coupling_states"], [[[1, 2], [0, 1]], [[0, 1], [1, 2]]])
        assert data["almo_hamiltonian"][0][0] == -156.62929
        self.assertAlmostEqual(data["almo_coupling_eV"], 0.26895)

    def test_pod_parsing(self):
        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "pod2_gs.out")).data
        assert data["pod_coupling_eV"] == 0.247818

    def test_fodft_parsing(self):
        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "fodft.out")).data
        assert data["fodft_coupling_eV"] == 0.268383

    def test_isosvp_water(self):
        data = QCOutput(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "isosvp_water_single.qcout")
        ).data
        assert data["solvent_method"] == "ISOSVP"
        # ISOSVP parameters
        assert data["solvent_data"]["isosvp"]["isosvp_dielectric"] == 78.39
        assert data["solvent_data"]["isosvp"]["final_soln_phase_e"] == -40.4850599393
        assert data["solvent_data"]["isosvp"]["solute_internal_e"] == -40.4846329762
        assert data["solvent_data"]["isosvp"]["change_solute_internal_e"] == 0.0000121967
        assert data["solvent_data"]["isosvp"]["reaction_field_free_e"] == -0.0004269631
        assert data["solvent_data"]["isosvp"]["total_solvation_free_e"] == -0.0004147664

        # CMIRS parameters
        assert data["solvent_data"]["cmirs"]["CMIRS_enabled"] is False

    def test_isosvp_dielst10(self):
        data = QCOutput(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "isosvp_dielst10_single.qcout")
        ).data
        assert data["solvent_method"] == "ISOSVP"

        # ISOSVP parameters
        assert data["solvent_data"]["isosvp"]["isosvp_dielectric"] == 10
        assert data["solvent_data"]["isosvp"]["final_soln_phase_e"] == -40.4850012952
        assert data["solvent_data"]["isosvp"]["solute_internal_e"] == -40.4846362547
        assert data["solvent_data"]["isosvp"]["change_solute_internal_e"] == 0.0000089182
        assert data["solvent_data"]["isosvp"]["reaction_field_free_e"] == -0.0003650405
        assert data["solvent_data"]["isosvp"]["total_solvation_free_e"] == -0.0003561223

        # CMIRS parameters
        assert data["solvent_data"]["cmirs"]["CMIRS_enabled"] is False

    def test_cmirs_benzene(self):
        data = QCOutput(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "cmirs_benzene_single.qcout")
        ).data
        assert data["solvent_method"] == "ISOSVP"
        assert data["solvent_data"]["isosvp"]["isosvp_dielectric"] == 2.28
        assert data["solvent_data"]["cmirs"]["CMIRS_enabled"] is True
        assert data["solvent_data"]["cmirs"]["dispersion_e"] == 0.6955542829
        assert data["solvent_data"]["cmirs"]["exchange_e"] == 0.2654553686
        assert data["solvent_data"]["cmirs"]["min_neg_field_e"] == 0.0006019665
        assert data["solvent_data"]["cmirs"]["max_pos_field_e"] == 0.0178177740

    def test_cmirs_dielst10(self):
        data = QCOutput(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "cmirs_dielst10_single.qcout")
        ).data
        assert data["solvent_method"] == "ISOSVP"
        assert data["solvent_data"]["isosvp"]["isosvp_dielectric"] == 10
        assert data["solvent_data"]["cmirs"]["CMIRS_enabled"] is True
        assert data["solvent_data"]["cmirs"]["dispersion_e"] == 0.6955550107
        assert data["solvent_data"]["cmirs"]["exchange_e"] == 0.2652679507
        assert data["solvent_data"]["cmirs"]["min_neg_field_e"] == 0.0005235850
        assert data["solvent_data"]["cmirs"]["max_pos_field_e"] == 0.0179866718

    def test_cmirs_water(self):
        data = QCOutput(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "cmirs_water_single.qcout")
        ).data
        assert data["solvent_method"] == "ISOSVP"

        # ISOSVP parameters
        assert data["solvent_data"]["isosvp"]["isosvp_dielectric"] == 78.39
        assert data["solvent_data"]["isosvp"]["final_soln_phase_e"] == -40.4752415075
        assert data["solvent_data"]["isosvp"]["solute_internal_e"] == -40.4748535587
        assert data["solvent_data"]["isosvp"]["change_solute_internal_e"] == 0.0000122982
        assert data["solvent_data"]["isosvp"]["reaction_field_free_e"] == -0.0003879488
        assert data["solvent_data"]["isosvp"]["total_solvation_free_e"] == 0.0037602703

        # CMIRS parameters
        assert data["solvent_data"]["cmirs"]["CMIRS_enabled"] is True
        assert data["solvent_data"]["cmirs"]["dispersion_e"] == 0.6722278965
        assert data["solvent_data"]["cmirs"]["exchange_e"] == 0.2652032616
        assert data["solvent_data"]["cmirs"]["min_neg_field_e"] == 0.0004967767
        assert data["solvent_data"]["cmirs"]["max_pos_field_e"] == 0.0180445935

    def test_NBO_hyperbonds(self):
        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "hyper.qout")).data
        assert len(data["nbo_data"]["hyperbonds"][0]["hyperbond index"].keys()) == 2
        assert data["nbo_data"]["hyperbonds"][0]["BD(A-B)"][1] == 106
        assert data["nbo_data"]["hyperbonds"][0]["bond atom 2 symbol"][0] == "C"
        assert data["nbo_data"]["hyperbonds"][0]["occ"][1] == 3.0802

    def test_NBO_3C(self):
        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "3C.qout")).data
        assert len(data["nbo_data"]["hybridization_character"]) == 3
        assert data["nbo_data"]["hybridization_character"][2]["type"][0] == "3C"
        assert data["nbo_data"]["hybridization_character"][2]["type"][10] == "3Cn"
        assert data["nbo_data"]["hybridization_character"][2]["type"][20] == "3C*"
        assert data["nbo_data"]["hybridization_character"][2]["atom 3 pol coeff"][15] == "0.3643"
        assert data["nbo_data"]["hybridization_character"][2]["atom 3 polarization"][8] == "56.72"
        assert data["nbo_data"]["hybridization_character"][2]["atom 3 symbol"][3] == "B"
        assert data["nbo_data"]["perturbation_energy"][0]["donor atom 2 number"][2592] == 36
        assert data["nbo_data"]["perturbation_energy"][0]["donor atom 2 symbol"][2125] == "B12"
        assert data["nbo_data"]["perturbation_energy"][0]["donor atom 2 number"][2593] == "info_is_from_3C"
        assert data["nbo_data"]["perturbation_energy"][0]["acceptor type"][723] == "3C*"
        assert data["nbo_data"]["perturbation_energy"][0]["perturbation energy"][3209] == 3.94


if __name__ == "__main__":
    # TestQCOutput.generate_single_job_dict()
    # TestQCOutput.generate_multi_job_dict()
    unittest.main()
