
import os
import unittest

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.Fstar.fstar import FStarDiagram
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.cif import CifParser

test_dir = os.path.join(os.path.dirname(__file__))



class test_fstardiagram(PymatgenTest):

    def setUp(self):

        self.cif_list = [i for i in os.listdir(test_dir) if i.endswith('.cif')]
        self.struct_list = [CifParser(os.path.join(test_dir, file)).get_structures(primitive=False)[0]
                            for file in self.cif_list]
        self.fstar = FStarDiagram(structure_objects=self.struct_list)

    def test_edit_fstar_diagram(self):

        self.assertEqual(self.fstar.site_labels, ['[0. 0. 0.]Li', '[0.  0.  0.5]Co', '[0.   0.   0.25]O'])
        new = FStarDiagram(structure_objects=self.struct_list)
        self.assertEqual(self.fstar.plot, new.plot)
        new.edit_fstar_diagram(combine_list=[['[0.  0.  0.5]Co', '[0. 0. 0.]Li']])
        self.assertEqual(new.site_labels, ['[0. 0. 0.]Li', '[0.  0.  0.5]Co', '[0.   0.   0.25]O',
                                           "['[0.  0.  0.5]Co', '[0. 0. 0.]Li']"])
        self.assertEqual(list(new.coords["['[0.  0.  0.5]Co', '[0. 0. 0.]Li']"].values),
                         list(self.fstar.coords['[0. 0. 0.]Li'].values+self.fstar.coords['[0.  0.  0.5]Co'].values))
        self.assertEqual(self.fstar.plot, new.plot)
        new.edit_fstar_diagram(plot_list=['[0.  0.  0.5]Co', '[0.   0.   0.25]O', '[0. 0. 0.]Li'])
        self.assertTrue(self.fstar.plot != new.plot)


if __name__ == '__main__':
    unittest.main()
