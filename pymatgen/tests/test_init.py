import os
import unittest
import warnings

from pymatgen import SETTINGS, get_structure_from_mp, loadfn
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import PymatgenTest


class SettingsTestCase(unittest.TestCase):

    # def test_something(self):
    #     SETTINGS = _load_pmg_settings()
    #     if os.path.exists(SETTINGS_FILE):
    #         with open(SETTINGS_FILE, "rt") as f:
    #             d = yaml.safe_load(f)
    #             for k, v in d.items():
    #                 self.assertEqual(v, SETTINGS[k])
    #     else:
    #         for k, v in SETTINGS.items():
    #             self.assertEqual(v, os.environ.get(k))

    @unittest.skipIf(not SETTINGS.get("PMG_MAPI_KEY"), "PMG_MAPI_KEY environment variable not set.")
    def test_get_structure_from_mp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertEqual(get_structure_from_mp("Li2O").formula, "Li2 O1")
            self.assertRaises(ValueError, get_structure_from_mp, "LiNaKCs")

    def test_loadfn(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            obj = loadfn(os.path.join(PymatgenTest.TEST_FILES_DIR, "Li2O.cif"))
            self.assertIsInstance(obj, Structure)
            obj = loadfn(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
            self.assertIsInstance(obj, Structure)
            obj = loadfn(os.path.join(PymatgenTest.TEST_FILES_DIR, "LiFePO4.vasp"))
            self.assertIsInstance(obj, Structure)
            obj = loadfn(os.path.join(PymatgenTest.TEST_FILES_DIR, "vasprun.xml"))
            self.assertIsInstance(obj, Vasprun)


if __name__ == "__main__":
    unittest.main()
