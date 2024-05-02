from __future__ import annotations

from pymatgen.util.io_utils import micro_pyawk
from pymatgen.util.testing import VASP_OUT_DIR, PymatgenTest


class TestFunc(PymatgenTest):
    def test_micro_pyawk(self):
        data = []

        def f(_x, y):
            data.append(y.group(1).strip())

        def f2(_x, y):
            return y

        micro_pyawk(f"{VASP_OUT_DIR}/OUTCAR.gz", [["POTCAR:(.*)", f2, f]])
        assert len(data) == 6
