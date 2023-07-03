
import os

from pymatgen.analysis.fstar.fstar import FStarDiagram
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.cif import CifParser


class Test_FStarDiagram(PymatgenTest):

    def setUp(self):

        self.cif_list = [file for file in os.listdir(self.TEST_FILES_DIR / "fstar_data")]
        self.struct_list = [CifParser(self.TEST_FILES_DIR / "fstar_data" / file).get_structures(primitive=False)[0]
                            for file in self.cif_list]
        self.fstar = FStarDiagram(structures=self.struct_list)

    def test_edit_fstar_diagram(self):

        self.assertEqual(self.fstar.site_labels, ['[0. 0. 0.]Li', '[0.  0.  0.5]Co', '[0.   0.   0.25]O'])
        new = FStarDiagram(structures=self.struct_list)
        self.assertEqual(self.fstar.plot, new.plot)
        new.edit_fstar_diagram(combine_list=[['[0.  0.  0.5]Co', '[0. 0. 0.]Li']])
        self.assertEqual(new.site_labels, ['[0. 0. 0.]Li', '[0.  0.  0.5]Co', '[0.   0.   0.25]O',
                                           "['[0.  0.  0.5]Co', '[0. 0. 0.]Li']"])
        self.assertEqual(list(new.coords["['[0.  0.  0.5]Co', '[0. 0. 0.]Li']"].values),
                         list(self.fstar.coords['[0. 0. 0.]Li'].values+self.fstar.coords['[0.  0.  0.5]Co'].values))
        self.assertEqual(self.fstar.plot, new.plot)
        new.edit_fstar_diagram(plot_list=['[0.  0.  0.5]Co', '[0.   0.   0.25]O', '[0. 0. 0.]Li'])
        self.assertTrue(self.fstar.plot != new.plot)