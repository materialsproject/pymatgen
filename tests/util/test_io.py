from __future__ import annotations

from pymatgen.util.io_utils import micro_pyawk
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestFunc(PymatgenTest):
    def test_micro_pyawk(self):
        filename = f"{TEST_FILES_DIR}/OUTCAR"
        data = []

        def f(x, y):
            data.append(y.group(1).strip())

        def f2(x, y):
            return y

        micro_pyawk(filename, [["POTCAR:(.*)", f2, f]])
        assert len(data) == 6
