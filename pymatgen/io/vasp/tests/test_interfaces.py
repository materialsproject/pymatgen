"""
Created on Aug 7, 2015
"""


__author__ = "David Cossey"
__copyright__ = "Approved for public release: AFRL DSRC Case Number 88ABW-2016-0260"
__version__ = "0.1"
__maintainer__ = "David Cossey"
__email__ = "dcossey014@gmail.com"
__date__ = "2016-01-29"

import unittest
import shutil
import pkgutil
import os

if pkgutil.find_loader('fireworks'):
    from pymatgen.io.vasp.interfaces import VaspInputInterface, VaspFirework, VaspWorkflow, VaspTransferTask
    from fireworks import ScriptTask

from pymatgen import Structure

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
test_dir = os.path.join(MODULE_DIR, '..', 'InterfaceExamples', 'test_files')

@unittest.skipUnless(1 if pkgutil.find_loader("fireworks") else 0, "Requires Fireworks")
class VaspInputInterfaceTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff=None
        self.mol_dict = {u'@class': 'Structure',
            u'@module': 'pymatgen.core.structure',
            u'lattice': {u'a': 3.09681658,
                u'alpha': 60.00000352,
                u'b': 3.09681633,
                u'beta': 60.000006230000004,
                u'c': 3.096817,
                u'gamma': 60.00000143999999,
                u'matrix': [[2.6819219975054747, 0.0, 1.5484079983838697],
                            [0.8939740004171242, 2.5285400003050134, 1.5484080002345615],
                            [0.0, 0.0, 3.096817]],
                u'volume': 21.00059082235557},
            u'sites': [{u'abc': [0.25, 0.25, 0.25],
                        u'species': [{u'element': u'Si', u'occu': 1.0}]},
                        {u'abc': [0.0, 0.0, 0.0], 
                        u'species': [{u'element': u'C', u'occu': 1.0}]}]}

        self.mol_file = os.path.join(test_dir, 'SiC_0.cif')
        self.config_file = os.path.join(test_dir, 'vasp_interface_defaults.yaml')

        self.isp = {'ALGO': 'Fast', 'LOPTICS':False}
        self.custom_params = {'kpts': [8,8,8], 'kpts_shift':[1,0,0]}

        self.gdic = {u'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
                'custom_params': {'user_kpts_settings': {'kpts': [6, 6, 6],
                    'kpts_shift': [0, 0, 0],
                    'kpts_style': 'Gamma'}},
                'input_set_params': {'user_incar_settings': {'ALGO': 'Normal',
                    'CSHIFT': 0.1,
                    'LOPTICS': '.TRUE.',
                    'NEDOS': 4096,
                    'SIGMA': 0.01}},
                'structure': {'@class': 'Structure',
                    '@module': 'pymatgen.core.structure',
                    'lattice': {'a': 3.09681658,
                        'alpha': 60.00000352,
                        'b': 3.09681633,
                        'beta': 60.000006230000004,
                        'c': 3.096817,
                        'gamma': 60.00000143999999,
                    'matrix': [[2.6819219975054747, 0.0, 1.5484079983838697],
                        [0.8939740004171242, 2.5285400003050134, 1.5484080002345615],
                        [0.0, 0.0, 3.096817]],
                    'volume': 21.00059082235557},
                    'sites': [{'abc': [0.25, 0.25, 0.25],
                        'species': [{'element': 'Si', 'occu': 1.0}]},
                        {'abc': [0.0, 0.0, 0.0], 'species': [{'element': 'C', 'occu': 1.0}]}]},
                'vasp_input_set': 'MinimalVaspInputSet'}

    def tearDown(self):
        pass

    def test_input_from_dict(self):
        input = VaspInputInterface(s=self.mol_dict, config_file=self.config_file)
        inp_dic = input.input.as_dict()
        self.assertDictEqual(inp_dic, self.gdic, msg=None)

    def test_input_from_file(self):
        input = VaspInputInterface(s=self.mol_file, config_file=self.config_file)
        inp_dic = input.input.as_dict()
        cdic = {u'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
                'custom_params': {'user_kpts_settings': {'kpts': [6, 6, 6],
                    'kpts_shift': [0, 0, 0],
                    'kpts_style': 'Gamma'}},
                'input_set_params': {'user_incar_settings': {'ALGO': 'Normal',
                    'CSHIFT': 0.1,
                    'LOPTICS': '.TRUE.',
                    'NEDOS': 4096,
                    'SIGMA': 0.01}},
                'structure': self.mol_file,
                'vasp_input_set': 'MinimalVaspInputSet'}
        self.assertDictEqual(inp_dic, cdic, msg=None)

    def test_input_from_structure(self):
        s = Structure.from_file(self.mol_file)
        input = VaspInputInterface(s=s, config_file=self.config_file)
        inp_dic = input.input.as_dict()
        self.assertDictEqual(inp_dic, self.gdic, msg=None)

    def test_custom_incar_and_kpoints(self):
        input = VaspInputInterface(s=self.mol_dict, isp=self.isp,
                            custom_params=self.custom_params,
                            config_file=self.config_file)
        inp_dic = input.input.as_dict()
        cdic = {u'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
            'custom_params': {'user_kpts_settings': {'kpts': [8, 8, 8],
                'kpts_shift': [1, 0, 0],
                'kpts_style': 'Gamma'}},
            'input_set_params': {'user_incar_settings': {'ALGO': 'Fast',
                'CSHIFT': 0.1,
                'LOPTICS': False,
                'NEDOS': 4096,
                'SIGMA': 0.01}},
            'structure': {'@class': 'Structure',
                '@module': 'pymatgen.core.structure',
                'lattice': {'a': 3.09681658,
                'alpha': 60.00000352,
                'b': 3.09681633,
                'beta': 60.000006230000004,
                'c': 3.096817,
                'gamma': 60.00000143999999,
            'matrix': [[2.6819219975054747, 0.0, 1.5484079983838697],
                [0.8939740004171242, 2.5285400003050134, 1.5484080002345615],
                [0.0, 0.0, 3.096817]],
            'volume': 21.00059082235557},
            'sites': [{'abc': [0.25, 0.25, 0.25],
                'species': [{'element': 'Si', 'occu': 1.0}]},
                {'abc': [0.0, 0.0, 0.0], 'species': [{'element': 'C', 'occu': 1.0}]}]},
            'vasp_input_set': 'MinimalVaspInputSet'}
        self.assertDictEqual(inp_dic, cdic, msg=None)

    def test_user_friendly_custom_incar(self):
        input = VaspInputInterface(s=self.mol_file, config_file=self.config_file)
        input.LOPTICS=False
        input.ALGO='Exact'
        input.EDIFF=1e-05
        input.kpts=[8,8,8]
        input.kpts_shift=[1,0,0]
        inp_dic = input.input.as_dict()

        cdic = {u'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
            'custom_params': {'user_kpts_settings': {'kpts': [8, 8, 8],
                'kpts_shift': [1, 0, 0],
                'kpts_style': 'Gamma'}},
            'input_set_params': {'user_incar_settings': {'ALGO': 'Exact',
                'CSHIFT': 0.1,
                'EDIFF': 1e-05,
                'LOPTICS': False,
                'NEDOS': 4096,
                'SIGMA': 0.01}},
            'structure': self.mol_file,
            'vasp_input_set': 'MinimalVaspInputSet'}
        self.assertDictEqual(inp_dic, cdic, msg=None)

    def test_get_nelect_from_interface(self):
        input = VaspInputInterface(s=self.mol_file, config_file=self.config_file)
        nelect = input.get_nelect()
        self.assertEqual(nelect, 8.0, msg=None)


@unittest.skipUnless(1 if pkgutil.find_loader("fireworks") else 0, "Requires Fireworks")
class VaspFireworkTest(unittest.TestCase):
    def setUp(self):
        self.mol_file = 'SiC_0.cif'
        self.config_file = os.path.join(test_dir, 'vasp_interface_defaults.yaml')
        self.input = VaspInputInterface(s=self.mol_file, config_file=self.config_file)
        self.misc_task = ScriptTask.from_str("echo 'Hello World!'")

    def tearDown(self):
        pass


    def test_add_task_VaspInputInterface(self):
        fw = VaspFirework(self.input, config_file=self.config_file)
        fw.add_task(self.input)
        fw_dic = fw.Firework.as_dict()
        for i in ['created_on', 'updated_on', 'fw_id']:
            del fw_dic[i]

        gdic = {'name': 'vaspfw',
            'spec': {'_tasks': [{'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
                'custom_params': {'user_kpts_settings': {'kpts': [6, 6, 6],
                    'kpts_shift': [0, 0, 0],
                    'kpts_style': 'Gamma'}},
                'input_set_params': {'user_incar_settings': {'ALGO': 'Normal',
                    'CSHIFT': 0.1,
                    'LOPTICS': '.TRUE.',
                    'NEDOS': 4096,
                    'SIGMA': 0.01}},
                'structure': 'SiC_0.cif',
                'vasp_input_set': 'MinimalVaspInputSet'},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.RunCustodianTask}}',
                'handler_params': {},
                'handlers': []},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
                'custom_params': {'user_kpts_settings': {'kpts': [6, 6, 6],
                    'kpts_shift': [0, 0, 0],
                    'kpts_style': 'Gamma'}},
                'input_set_params': {'user_incar_settings': {'ALGO': 'Normal',
                    'CSHIFT': 0.1,
                    'LOPTICS': '.TRUE.',
                    'NEDOS': 4096,
                    'SIGMA': 0.01}},
                'structure': 'SiC_0.cif',
                'vasp_input_set': 'MinimalVaspInputSet'},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.RunCustodianTask}}',
                'handler_params': {},
                'handlers': []}]}}

        self.assertDictEqual(fw_dic, gdic, msg=None)

    def test_add_task_misc_task(self):
        fw = VaspFirework(self.input, config_file=self.config_file)
        fw.add_task(self.misc_task)
        fw_dic = fw.Firework.as_dict()
        for i in ['created_on', 'updated_on', 'fw_id']:
            del fw_dic[i]

        gdic = {'name': 'vaspfw',
            'spec': {'_tasks': [{'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
                'custom_params': {'user_kpts_settings': {'kpts': [6, 6, 6],
                    'kpts_shift': [0, 0, 0],
                    'kpts_style': 'Gamma'}},
                'input_set_params': {'user_incar_settings': {'ALGO': 'Normal',
                    'CSHIFT': 0.1,
                    'LOPTICS': '.TRUE.',
                    'NEDOS': 4096,
                    'SIGMA': 0.01}},
                'structure': 'SiC_0.cif',
                'vasp_input_set': 'MinimalVaspInputSet'},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.RunCustodianTask}}',
                'handler_params': {},
                'handlers': []},
            {'_fw_name': 'ScriptTask',
                'script': ["echo 'Hello World!'"],
                'use_shell': True}]}}

        self.assertDictEqual(fw_dic, gdic, msg=None)

    def test_copy_from_previous(self):
        fw = VaspFirework(self.input, config_file=self.config_file)
        fw.copy_files_from_previous('WAVECAR', 'WAVEDER', mode='copy', ignore_errors=True)
        fw_dic = fw.Firework.as_dict()
        for i in ['created_on', 'updated_on', 'fw_id']:
            del fw_dic[i]
        
        gdic = {'name': 'vaspfw',
            'spec': {'_tasks': [{'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
                'custom_params': {'user_kpts_settings': {'kpts': [6, 6, 6],
                    'kpts_shift': [0, 0, 0],
                    'kpts_style': 'Gamma'}},
                'input_set_params': {'user_incar_settings': {'ALGO': 'Normal',
                    'CSHIFT': 0.1,
                    'LOPTICS': '.TRUE.',
                    'NEDOS': 4096,
                    'SIGMA': 0.01}},
                'structure': 'SiC_0.cif',
                'vasp_input_set': 'MinimalVaspInputSet'},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.VaspTransferTask}}',
                'files': ['WAVECAR', 'WAVEDER'],
                'ignore_errors': True,
                'mode': 'copy'},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.RunCustodianTask}}',
                'handler_params': {},
                'handlers': []}]}}

        self.assertDictEqual(fw_dic, gdic, msg=None)

    def test_copy_from_dir(self):
        out_name = 'FW_copy_from_dir.yaml'
        fw = VaspFirework(self.input, config_file=self.config_file)
        fw.copy_files_from_previous('WAVECAR', 'WAVEDER', mode='copy', 
                                    dir='../../test_files/VASP_INTERFACES_TEST_FILES', 
                                    ignore_errors=True)
        fw_dic = fw.Firework.as_dict()
        for i in ['created_on', 'updated_on', 'fw_id']:
            del fw_dic[i]
        
        gdic = {'name': 'vaspfw',
            'spec': {'_tasks': [{'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
                'custom_params': {'user_kpts_settings': {'kpts': [6, 6, 6],
                    'kpts_shift': [0, 0, 0],
                    'kpts_style': 'Gamma'}},
                'input_set_params': {'user_incar_settings': {'ALGO': 'Normal',
                    'CSHIFT': 0.1,
                    'LOPTICS': '.TRUE.',
                    'NEDOS': 4096,
                    'SIGMA': 0.01}},
                'structure': 'SiC_0.cif',
                'vasp_input_set': 'MinimalVaspInputSet'},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.VaspTransferTask}}',
                'dir': '../../test_files/VASP_INTERFACES_TEST_FILES',
                'files': ['WAVECAR', 'WAVEDER'],
                'ignore_errors': True,
                'mode': 'copy'},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.RunCustodianTask}}',
                'handler_params': {},
                'handlers': []}]}}

        self.assertDictEqual(fw_dic, gdic, msg=None)

    def test_combo_misc_add_task(self):
        out_name = 'FW_combo_add_tasks.yaml'
        fw = VaspFirework(self.input, config_file=self.config_file)
        fw.add_task(self.misc_task)
        fw.add_task(self.input)
        fw_dic = fw.Firework.as_dict()
        for i in ['created_on', 'updated_on', 'fw_id']:
            del fw_dic[i]

        gdic = {'name': 'vaspfw',
            'spec': {'_tasks': [{'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
                'custom_params': {'user_kpts_settings': {'kpts': [6, 6, 6],
                    'kpts_shift': [0, 0, 0],
                    'kpts_style': 'Gamma'}},
                'input_set_params': {'user_incar_settings': {'ALGO': 'Normal',
                    'CSHIFT': 0.1,
                    'LOPTICS': '.TRUE.',
                    'NEDOS': 4096,
                    'SIGMA': 0.01}},
                'structure': 'SiC_0.cif',
                'vasp_input_set': 'MinimalVaspInputSet'},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.RunCustodianTask}}',
                'handler_params': {},
                'handlers': []},
            {'_fw_name': 'ScriptTask',
                'script': ["echo 'Hello World!'"],
                'use_shell': True},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.WriteVaspInputTask}}',
                'custom_params': {'user_kpts_settings': {'kpts': [6, 6, 6],
                    'kpts_shift': [0, 0, 0],
                    'kpts_style': 'Gamma'}},
                'input_set_params': {'user_incar_settings': {'ALGO': 'Normal',
                    'CSHIFT': 0.1,
                    'LOPTICS': '.TRUE.',
                    'NEDOS': 4096,
                    'SIGMA': 0.01}},
                'structure': 'SiC_0.cif',
                'vasp_input_set': 'MinimalVaspInputSet'},
            {'_fw_name': '{{pymatgen.io.vasp.interfaces.RunCustodianTask}}',
                'handler_params': {},
                'handlers': []}]}}

        self.assertDictEqual(fw_dic, gdic, msg=None)


@unittest.skipUnless(1 if pkgutil.find_loader("fireworks") else 0, "Requires Fireworks")
class VaspWorkflowTest(unittest.TestCase):
    def setUp(self):
        self.mol_file = 'SiC_0.cif'
        self.config_file = os.path.join(test_dir, 'vasp_interface_defaults.yaml')
        self.input = VaspInputInterface(s=self.mol_file, config_file=self.config_file)
        self.fw1 = VaspFirework(self.input, config_file=self.config_file)
        self.fw2 = VaspFirework(self.input, config_file=self.config_file)
        self.fw3 = VaspFirework(self.input, config_file=self.config_file)

        self.links1 = {'-1': [-2], '-2': [-3], '-3': []}
        self.links2 = {'-1': [-2, -3], '-2': [], '-3': []}

    def tearDown(self):
        pass


    def test_make_seq_workflow_from_multi_fws(self):
        wf1 = VaspWorkflow(self.fw1, self.fw2, self.fw3)
        wf_dic = wf1.wf.as_dict()
        ldic = wf_dic['links']

        # Rebuild WF1 Dictionary to reflect FWs '-1', '-2', and '-3' instead
        # of higher FW_ID numbers
        min_fw_id = int(min(ldic.keys()))
        modifier = -( min_fw_id + 1 )
        link_dic = {str(min_fw_id + modifier): [ldic[str(min_fw_id)][0] + modifier],
                str(min_fw_id + modifier - 1): [ldic[str(min_fw_id-1)][0] + modifier],
                str(min_fw_id + modifier - 2): []}

        self.assertDictEqual(link_dic, self.links1, msg=None)

    def test_add_fw(self):
        wf1 = VaspWorkflow(self.fw1)
        wf1.add_fw(self.fw2)
        wf1.add_fw(self.fw3)

        wf_dic = wf1.wf.as_dict()
        ldic = wf_dic['links']

        # Rebuild WF1 Dictionary to reflect FWs '-1', '-2', and '-3' instead
        # of higher FW_ID numbers
        min_fw_id = int(min(ldic.keys()))
        modifier = -( min_fw_id + 1 )
        link_dic = {str(min_fw_id + modifier): [ldic[str(min_fw_id)][0] + modifier],
                str(min_fw_id + modifier - 1): [ldic[str(min_fw_id-1)][0] + modifier],
                str(min_fw_id + modifier - 2): []}
        
        self.assertDictEqual(link_dic, self.links1, msg=None)

    def test_multifw_wf_custom_deps(self):
        wf1 = VaspWorkflow(self.fw1, self.fw2, self.fw3, 
                        deps_dict={self.fw1: [self.fw2, self.fw3]})
         
        wf_dic = wf1.wf.as_dict()
        ldic = wf_dic['links']

        # Rebuild WF1 Dictionary to reflect FWs '-1', '-2', and '-3' instead
        # of higher FW_ID numbers
        min_fw_id = int(min(ldic.keys()))
        modifier = -( min_fw_id + 1 )
        link_dic = {str(min_fw_id + modifier): [
                    ldic[str(min_fw_id)][0] + modifier, ldic[str(min_fw_id)][1] + modifier],
                str(min_fw_id + modifier - 1): [],
                str(min_fw_id + modifier - 2): []}

        self.assertDictEqual(link_dic, self.links2, msg=None)


@unittest.skipUnless(1 if pkgutil.find_loader("fireworks") else 0, "Requires Fireworks")
class VaspTransferTaskTest(unittest.TestCase):
    def setUp(self):
        self.src_dir  = os.path.join(test_dir, 'TRNFR_TEST_FILES')
        self.dest_dir = os.path.join(os.environ['HOME'], 'TRNFR_TEST')
        if not os.path.exists(self.dest_dir):
            os.mkdir(self.dest_dir)

    def tearDown(self):
        os.chdir(MODULE_DIR)
        if os.path.exists(self.dest_dir):
            shutil.rmtree(self.dest_dir)


    def test_file_transfer_fail1(self):
        os.chdir(self.dest_dir)
        fw_fail1 = VaspTransferTask(mode='copy', dir=self.src_dir,
                        files=['test_f1.txt', 'test_f2.txt', 'test_f3.txt'], 
                        ignore_errors=False)
        self.assertRaises(ValueError, fw_fail1.run_task, fw_spec={})

    def test_file_transfer_fail2(self):
        os.chdir(self.dest_dir)
        fw_fail2 = VaspTransferTask(mode='copy', 
                        files=['test_f1.txt'], ignore_errors=False)
        self.assertRaises(RuntimeError, fw_fail2.run_task, fw_spec={})

    def test_file_transfer_fail3(self):
        os.chdir(self.dest_dir)
        fw_fail3  = VaspTransferTask(mode='copy',
                        files=['test_f1.txt', 'test_f2.txt'],
                        ignore_errors=False)
        self.assertRaises(TypeError, fw_fail3.run_task)

    def test_file_transfer_succeed1(self):
        fw_succ1 = VaspTransferTask(mode='copy', dir=self.src_dir,
                        files=['test_f1.txt', 'test_f2.txt', 'test_f3.txt'],
                        ignore_errors=True)
        os.chdir(self.dest_dir)
        fw_succ1.run_task(fw_spec={})

        self.assertListEqual(os.listdir(self.dest_dir), os.listdir(self.src_dir))

    def test_file_transfer_succeed2(self):
        fw_succ2 = VaspTransferTask(mode='copy', dir=self.src_dir,
                        files=['test_f1.txt', 'test_f2.txt'],
                        ignore_errors=False)
        os.chdir(self.dest_dir)
        fw_succ2.run_task(fw_spec={})

        self.assertListEqual(os.listdir(self.dest_dir), os.listdir(self.src_dir))

    def test_file_transfer_from_previous(self):
        fw1 = VaspTransferTask(mode='copy', 
                        files=['test_f1.txt', 'test_f2.txt'], 
                        ignore_errors=False)
        os.chdir(self.dest_dir)
        fw1.run_task(fw_spec= {'PREV_DIR': self.src_dir})

        self.assertListEqual(os.listdir(self.dest_dir), os.listdir(self.src_dir))
        

if __name__ == '__main__':
    unittest.main()
