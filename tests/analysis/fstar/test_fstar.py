
import os

from pymatgen.analysis.fstar.fstar import FStarDiagram
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest
from pymatgen.io.cif import CifParser


class Test_FStarDiagram(PymatgenTest):

    def setUp(self):

        self.cif_list = [file for file in os.listdir(f"{TEST_FILES_DIR}/analysis/fstar") if file.endswith(".cif")]
        self.struct_list = [CifParser(f"{TEST_FILES_DIR}/analysis/fstar/"+file).get_structures(primitive=False, symmetrized=True, check_occu=False)[0]
                            for file in self.cif_list]
        self.fstar = FStarDiagram(structures=self.struct_list)

    def test_edit_fstar_diagram(self):

        self.assertEqual(self.fstar.site_labels, ['[0. 0. 0.]Li', '[0.  0.  0.5]Co', '[0.   0.   0.25]O'])
        new = FStarDiagram(structures=self.struct_list)
        self.assertEqual(self.fstar.plot, new.plot)
        new.edit_fstar_diagram(combine_list=[['[0.  0.  0.5]Co', '[0. 0. 0.]Li']])
        self.assertEqual(new.site_labels, ['[0. 0. 0.]Li', '[0.  0.  0.5]Co', '[0.   0.   0.25]O',
                                           "['[0.  0.  0.5]Co', '[0. 0. 0.]Li']"])
        self.assertEqual(list(new.coords["['[0.  0.  0.5]Co', '[0. 0. 0.]Li']"].to_numpy()),
                         list(self.fstar.coords['[0. 0. 0.]Li'].values+self.fstar.coords['[0.  0.  0.5]Co'].to_numpy()))
        self.assertEqual(self.fstar.plot, new.plot)
        new.edit_fstar_diagram(plot_list=['[0.  0.  0.5]Co', '[0.   0.   0.25]O', '[0. 0. 0.]Li'])
        self.assertTrue(self.fstar.plot != new.plot)