from __future__ import annotations

from pymatgen.analysis.fstar.fstar import FStarDiagram
from pymatgen.io.cif import CifParser
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class Test_FStarDiagram(PymatgenTest):
    def setUp(self):
        self.struct_list = [
            CifParser(f"{TEST_FILES_DIR}/rhomb_3478.cif").get_structures(primitive=False, symmeterized=True)[0]
        ]
        self.fstar = FStarDiagram(structures=self.struct_list)

    def test_edit_fstar_diagram(self):
        assert self.fstar.site_labels == ["[0.890001 0.890001 0.890001]O", "[0. 0. 0.]Cu", "[0.5 0.5 0.5]Al"]
        new = FStarDiagram(structures=self.struct_list)
        assert self.fstar.plot == new.plot
        new.combine_sites(site_lists=[["[0. 0. 0.]Cu", "[0.5 0.5 0.5]Al"]])
        assert new.site_labels == [
            "[0.890001 0.890001 0.890001]O",
            "[0. 0. 0.]Cu",
            "[0.5 0.5 0.5]Al",
            "['[0. 0. 0.]Cu', '[0.5 0.5 0.5]Al']",
        ]
        assert list(new.fstar_coords["['[0. 0. 0.]Cu', '[0.5 0.5 0.5]Al']"].to_numpy()) == list(
            self.fstar.fstar_coords["[0. 0. 0.]Cu"].to_numpy() + self.fstar.fstar_coords["[0.5 0.5 0.5]Al"].to_numpy()
        )
        new.set_plot_list(site_list=["[0. 0. 0.]Cu", "[0.890001 0.890001 0.890001]O", "[0.5 0.5 0.5]Al"])
        assert self.fstar.plot_list != new.plot_list
        new.make_plot()
        assert self.fstar.plot != new.plot
