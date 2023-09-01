from __future__ import annotations

import tarfile

from pymatgen.analysis.fstar.fstar import FStarDiagram
from pymatgen.io.cif import CifParser
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class Test_FStarDiagram(PymatgenTest):
    def setUp(self):
        self.cif_gz = tarfile.open(f"{TEST_FILES_DIR}/fstar/fstar.tar.gz", "r")
        self.struct_list = [
            CifParser.from_str(self.cif_gz.extractfile(file).read().decode("utf-8")).get_structures(
                primitive=False, symmetrized=True, check_occu=False
            )[0]
            for file in self.cif_gz.getnames()
        ]
        self.fstar = FStarDiagram(structures=self.struct_list)

    def test_edit_fstar_diagram(self):
        assert self.fstar.site_labels == ["[0.   0.   0.25]O", "[0.  0.  0.5]Co", "[0. 0. 0.]Li"]
        new = FStarDiagram(structures=self.struct_list)
        assert self.fstar.plot == new.plot
        new.combine_sites(site_lists=[["[0.  0.  0.5]Co", "[0. 0. 0.]Li"]])
        assert new.site_labels == [
            "[0.   0.   0.25]O",
            "[0.  0.  0.5]Co",
            "[0. 0. 0.]Li",
            "['[0.  0.  0.5]Co', '[0. 0. 0.]Li']",
        ]
        assert list(new.fstar_coords["['[0.  0.  0.5]Co', '[0. 0. 0.]Li']"].to_numpy()) == list(
            self.fstar.fstar_coords["[0. 0. 0.]Li"].to_numpy() + self.fstar.fstar_coords["[0.  0.  0.5]Co"].to_numpy()
        )
        new.set_plot_list(site_list=["[0.  0.  0.5]Co", "[0.   0.   0.25]O", "[0. 0. 0.]Li"])
        assert self.fstar.plot_list != new.plot_list
        new.make_plot()
        assert self.fstar.plot != new.plot
