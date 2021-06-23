# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest

from monty.os.path import which
from monty.serialization import dumpfn, loadfn

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.outputs import QCOutput, check_for_structure_changes
from pymatgen.util.testing import PymatgenTest

try:
    from openbabel import openbabel

    have_babel = True
except ImportError:
    have_babel = False

__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath, Evan Spotte-Smith"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"

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
}

if have_babel:
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
                QCOutput, os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", file), keep_sub_files=False
            )
            data = []
            for sub_output in outputs:
                data.append(sub_output.data)
            multi_job_dict[file] = data
        dumpfn(multi_job_dict, "multi_job.json")

    def _test_property(self, key, single_outs, multi_outs):
        for name, outdata in single_outs.items():
            try:
                self.assertEqual(outdata.get(key), single_job_dict[name].get(key))
            except ValueError:
                self.assertArrayEqual(outdata.get(key), single_job_dict[name].get(key))
        for name, outputs in multi_outs.items():
            for ii, sub_output in enumerate(outputs):
                try:
                    self.assertEqual(sub_output.data.get(key), multi_job_dict[name][ii].get(key))
                except ValueError:
                    self.assertArrayEqual(sub_output.data.get(key), multi_job_dict[name][ii].get(key))

    def test_all(self):
        self.maxDiff = None
        single_outs = dict()
        for file in single_job_out_names:
            single_outs[file] = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", file)).data

        multi_outs = dict()
        for file in multi_job_out_names:
            multi_outs[file] = QCOutput.multiple_outputs_from_file(
                QCOutput, os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", file), keep_sub_files=False
            )

        for key in property_list:
            print("Testing ", key)
            self._test_property(key, single_outs, multi_outs)

    @unittest.skipIf((not (have_babel)), "OpenBabel not installed.")
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

        self.assertEqual(check_for_structure_changes(t1, t1), "no_change")
        self.assertEqual(check_for_structure_changes(t2, t3), "no_change")
        self.assertEqual(check_for_structure_changes(t1, t2), "fewer_bonds")
        self.assertEqual(check_for_structure_changes(t2, t1), "more_bonds")

        self.assertEqual(check_for_structure_changes(thio_1, thio_2), "unconnected_fragments")

        self.assertEqual(check_for_structure_changes(frag_1, frag_2), "bond_change")

    def test_NBO_parsing(self):
        data = QCOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files", "nbo.qout")).data
        self.assertEqual(len(data["nbo_data"]["natural_populations"]), 3)
        self.assertEqual(len(data["nbo_data"]["hybridization_character"]), 4)
        self.assertEqual(len(data["nbo_data"]["perturbation_energy"]), 2)
        self.assertEqual(data["nbo_data"]["natural_populations"][0]["Density"][5], -0.08624)
        self.assertEqual(data["nbo_data"]["hybridization_character"][-1]["atom 2 pol coeff"][35], "-0.7059")
        self.assertEqual(data["nbo_data"]["perturbation_energy"][-1]["fock matrix element"][104], 0.071)


if __name__ == "__main__":
    unittest.main()
