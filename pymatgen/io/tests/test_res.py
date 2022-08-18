import os

from pymatgen.io.res import ResParser
from pymatgen.util.testing import PymatgenTest


class TestResRead:
    # results from AIRSS runs
    res_coc = "coc-115925-9326-14.res"

    def test_airss_result_as_structure(self):
        file = os.path.join(PymatgenTest.TEST_FILES_DIR, "res", self.res_coc)
        res = ResParser._parse_filename(file)
        assert res.TITL is not None
        assert res.TITL.seed == self.res_coc[:-4]
        assert res.TITL.energy - -3.90427411e003 < 0.000001
        assert res.TITL.spacegroup_label == "R3"
        assert res.TITL.pressure - 15 < 1
        assert res.CELL.alpha - 49.32125 < 0.000001
        assert any(["AIRSS Version 0.9.1" in rem for rem in res.REMS])
