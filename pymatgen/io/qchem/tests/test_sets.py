# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os

from pymatgen.io.qchem.sets import *
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.qchem.inputs import QCInput

__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"

test_dir = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "..", 'test_files',
    "molecules")


class QChemDictSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_DictSet = QChemDictSet(
            molecule=test_molecule,
            job_type='opt',
            basis_set='6-31G*',
            scf_algorithm='diis')
        self.assertEqual(
            test_DictSet.rem, {
                'job_type': 'opt',
                'gen_scfman': 'true',
                'basis': '6-31g*',
                'max_scf_cycles': 200,
                'method': 'wb97xd',
                'scf_algorithm': 'diis',
                'geom_opt_max_cycles': 200
            })
        self.assertEqual(test_DictSet.pcm, {})
        self.assertEqual(test_DictSet.solvent, {})
        self.assertEqual(test_DictSet.molecule, test_molecule)

    def test_full_init(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_DictSet = QChemDictSet(
            molecule=test_molecule,
            job_type='opt',
            basis_set='6-31g*',
            scf_algorithm='diis',
            dft_rung=1,
            pcm_dielectric=10.0,
            max_scf_cycles=35)
        self.assertEqual(
            test_DictSet.rem, {
                'job_type': 'opt',
                'gen_scfman': 'true',
                'basis': '6-31g*',
                'max_scf_cycles': 35,
                'exchange': 'b3lyp',
                'geom_opt_max_cycles': 200,
                'scf_algorithm': 'diis',
                'solvent_method': 'pcm'
            })
        self.assertEqual(
            test_DictSet.pcm, {
                'heavypoints': '194',
                'hpoints': '194',
                'radii': 'uff',
                'theory': 'cpcm',
                'vdwscale': '1.1'
            })
        self.assertEqual(test_DictSet.solvent, {'dielectric': 10.0})
        self.assertEqual(test_DictSet.molecule, test_molecule)

    def test_overwrite_input_addition(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        overwrite_inputs = {"rem": {'thresh': 14}}
        test_OptSet = OptSet(
            molecule=test_molecule, overwrite_inputs=overwrite_inputs)
        act_rem = {
            'job_type': 'opt',
            'gen_scfman': 'true',
            'basis': '6-311++g*',
            'max_scf_cycles': 200,
            'method': 'wb97xd',
            'scf_algorithm': 'gdm',
            'geom_opt_max_cycles': 200,
            'thresh': 14
        }
        self.assertDictEqual(act_rem, test_OptSet.rem)

    def test_overwrite_input(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        overwrite_inputs = {
            "rem": {
                'method': 'b3lyp',
                'basis': '6-31g*',
                'thresh': 14
            }
        }
        test_OptSet = OptSet(
            molecule=test_molecule, overwrite_inputs=overwrite_inputs)
        act_rem = {
            'job_type': 'opt',
            'gen_scfman': 'true',
            'basis': '6-31g*',
            'max_scf_cycles': 200,
            'method': 'b3lyp',
            'scf_algorithm': 'gdm',
            'geom_opt_max_cycles': 200,
            'thresh': 14
        }
        self.assertDictEqual(act_rem, test_OptSet.rem)


class OptSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_OptSet = OptSet(molecule=test_molecule)
        self.assertEqual(
            test_OptSet.rem, {
                'job_type': 'opt',
                'gen_scfman': 'true',
                'basis': '6-311++g*',
                'max_scf_cycles': 200,
                'method': 'wb97xd',
                'scf_algorithm': 'gdm',
                'geom_opt_max_cycles': 200
            })
        self.assertEqual(test_OptSet.pcm, {})
        self.assertEqual(test_OptSet.solvent, {})
        self.assertEqual(test_OptSet.molecule, test_molecule)

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_OptSet = OptSet(molecule=test_molecule, pcm_dielectric=10.0)
        self.assertEqual(
            test_OptSet.rem, {
                'job_type': 'opt',
                'gen_scfman': 'true',
                'basis': '6-311++g*',
                'max_scf_cycles': 200,
                'method': 'wb97xd',
                'geom_opt_max_cycles': 200,
                'scf_algorithm': 'gdm',
                'solvent_method': 'pcm'
            })
        self.assertEqual(
            test_OptSet.pcm, {
                'heavypoints': '194',
                'hpoints': '194',
                'radii': 'uff',
                'theory': 'cpcm',
                'vdwscale': '1.1'
            })
        self.assertEqual(test_OptSet.solvent, {'dielectric': 10.0})
        self.assertEqual(test_OptSet.molecule, test_molecule)


class SinglePointSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_SPSet = SinglePointSet(molecule=test_molecule)
        self.assertEqual(
            test_SPSet.rem, {
                'job_type': 'sp',
                'gen_scfman': 'true',
                'basis': '6-311++g*',
                'max_scf_cycles': 200,
                'method': 'wb97xd',
                'scf_algorithm': 'gdm'
            })
        self.assertEqual(test_SPSet.pcm, {})
        self.assertEqual(test_SPSet.solvent, {})
        self.assertEqual(test_SPSet.molecule, test_molecule)

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_SPSet = SinglePointSet(
            molecule=test_molecule, pcm_dielectric=10.0)
        self.assertEqual(
            test_SPSet.rem, {
                'job_type': 'sp',
                'gen_scfman': 'true',
                'basis': '6-311++g*',
                'max_scf_cycles': 200,
                'method': 'wb97xd',
                'scf_algorithm': 'gdm',
                'solvent_method': 'pcm'
            })
        self.assertEqual(
            test_SPSet.pcm, {
                'heavypoints': '194',
                'hpoints': '194',
                'radii': 'uff',
                'theory': 'cpcm',
                'vdwscale': '1.1'
            })
        self.assertEqual(test_SPSet.solvent, {'dielectric': 10.0})
        self.assertEqual(test_SPSet.molecule, test_molecule)


class FreqSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_FreqSet = FreqSet(molecule=test_molecule)
        self.assertEqual(
            test_FreqSet.rem, {
                'job_type': 'freq',
                'gen_scfman': 'true',
                'basis': '6-311++g*',
                'max_scf_cycles': 200,
                'method': 'wb97xd',
                'scf_algorithm': 'gdm'
            })
        self.assertEqual(test_FreqSet.pcm, {})
        self.assertEqual(test_FreqSet.solvent, {})
        self.assertEqual(test_FreqSet.molecule, test_molecule)

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(
            os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_FreqSet = FreqSet(molecule=test_molecule, pcm_dielectric=10.0)
        self.assertEqual(
            test_FreqSet.rem, {
                'job_type': 'freq',
                'gen_scfman': 'true',
                'basis': '6-311++g*',
                'max_scf_cycles': 200,
                'method': 'wb97xd',
                'scf_algorithm': 'gdm',
                'solvent_method': 'pcm'
            })
        self.assertEqual(
            test_FreqSet.pcm, {
                'heavypoints': '194',
                'hpoints': '194',
                'radii': 'uff',
                'theory': 'cpcm',
                'vdwscale': '1.1'
            })
        self.assertEqual(test_FreqSet.solvent, {'dielectric': 10.0})
        self.assertEqual(test_FreqSet.molecule, test_molecule)


if __name__ == '__main__':
    unittest.main()
