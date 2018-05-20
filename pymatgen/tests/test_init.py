import unittest

import os
import ruamel.yaml as yaml
from pymatgen import SETTINGS_FILE, _load_pmg_settings, get_structure_from_mp, \
    SETTINGS
import warnings


class SettingsTestCase(unittest.TestCase):

    def test_something(self):
        SETTINGS = _load_pmg_settings()
        if os.path.exists(SETTINGS_FILE):
            with open(SETTINGS_FILE, "rt") as f:
                d = yaml.safe_load(f)
                for k, v in d.items():
                    self.assertEqual(v, SETTINGS[k])
        else:
            for k, v in SETTINGS.items():
                self.assertEqual(v, os.environ.get(k))

    @unittest.skipIf(not SETTINGS.get("PMG_MAPI_KEY"),
                     "PMG_MAPI_KEY environment variable not set.")
    def test_get_structure_from_mp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertEqual(get_structure_from_mp("Li2O").formula, "Li2 O1")
            self.assertRaises(ValueError, get_structure_from_mp, "LiNaKCs")


if __name__ == '__main__':
    unittest.main()
