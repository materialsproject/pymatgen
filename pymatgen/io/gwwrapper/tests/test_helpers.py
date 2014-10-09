from __future__ import division, print_function, unicode_literals

__author__ = 'setten'

import os

import unittest
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
from pymatgen.matproj.rest import MPRester, MPRestError
from pymatgen.io.gwwrapper.datastructures import GWSpecs, GWConvergenceData, get_spec
from pymatgen.io.gwwrapper.codeinterfaces import AbinitInterface, AbstractCodeInterface, VaspInterface, get_code_interface
from pymatgen.io.gwwrapper.GWworkflows import GWG0W0VaspInputSet, SingleAbinitGWWorkFlow
from pymatgen.io.gwwrapper.helpers import refine_structure, clean, load_ps, read_extra_abivars, read_grid_from_file
from pymatgen.io.gwwrapper.helpers import expand

#test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",'test_files')

have_abinit_ps_ext = True
try:
    os.environ['ABINIT_PS_EXT']
except KeyError:
    have_abinit_ps_ext = False


structure_dict = {'lattice': {'a': 3.866974622849504,
                              'gamma': 60.0000000241681,
                              'c': 3.86697462,
                              'b': 3.866974623775052,
                              'matrix': [[3.34889826, 0.0, 1.93348731], [1.11629942, 3.15737156, 1.93348731], [0.0, 0.0, 3.86697462]],
                              'volume': 40.888291888494884,
                              'alpha': 60.000000032293386,
                              'beta': 60.00000002437586},
                  'sites': [{'properties': {u'coordination_no': 5, u'forces': [0.0, 0.0, 0.0]},
                             'abc': [0.875, 0.875, 0.875],
                             'xyz': [3.9070479700000003, 2.762700115, 6.767205585],
                             'species': [{'occu': 1, 'element': 'Si'}], 'label': 'Si'},
                            {'properties': {u'coordination_no': 5, u'forces': [0.0, 0.0, 0.0]},
                             'abc': [0.125, 0.125, 0.125], 'xyz': [0.55814971, 0.394671445, 0.966743655],
                             'species': [{'occu': 1, 'element': 'Si'}], 'label': 'Si'}],
                  '@class': 'Structure', '@module': 'pymatgen.core.structure'}
structure = Structure.from_dict(structure_dict)


class GWTestHelpers(PymatgenTest):
    #def test_refine_structure(self):
    #    test = refine_structure(structure)
    #    self.assertIsInstance(test, Structure)

    def test_clean(self):
        string = 'MddmmdDD  '
        self.assertEqual(clean(string), 'mddmmddd')

    def test_read_extra_abivars(self):
        vars_out = {'ecut': 40}
        f = open('extra_abivars', 'w')
        f.write(str(vars_out))
        f.close()
        vars_in = read_extra_abivars()
        self.assertEqual(vars_out, vars_in)
        os.remove('extra_abivars')


    @unittest.skipIf(not have_abinit_ps_ext, "Requires ABINIT_PS_EXT env variable")
    def test_expand(self):
        spec = get_spec('GW')
        tests = SingleAbinitGWWorkFlow(structure, spec).convs
        tests_out = {'nscf_nbands': {'test_range': (50,),
                                     'control': 'gap', 'method': 'set_bands', 'level': 'nscf'},
                     'ecut': {'test_range': (28, 32, 36, 40, 44),
                              'control': 'e_ks_max', 'method': 'direct', 'level': 'scf'},
                     'ecuteps': {'test_range': (4, 8, 12, 16, 20),
                                 'control': 'gap', 'method': 'direct', 'level': 'sigma'}}
        self.assertEqual(expand(tests, 1), tests_out)
        spec.data['code'] = 'VASP'
        spec.update_code_interface()
        tests = GWG0W0VaspInputSet(structure, spec).convs
        tests_out = {'ENCUTGW': {'test_range': (200, 400, 600, 800), 'control': 'gap', 'method': 'incar_settings'}}
        self.assertEqual(expand(tests, 1), tests_out)
