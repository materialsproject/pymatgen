import unittest

import os
import yaml
from pymatgen import SETTINGS_FILE, _load_pmg_settings


class SettingsTestCase(unittest.TestCase):

    def test_something(self):
        SETTINGS = _load_pmg_settings()
        if os.path.exists(SETTINGS_FILE):
            with open(SETTINGS_FILE, "rt") as f:
                d = yaml.load(f)
                for k, v in d.items():
                    self.assertEqual(v, SETTINGS[k])
        else:
            for k, v in SETTINGS.items():
                self.assertEqual(v, os.environ.get(k))


if __name__ == '__main__':
    unittest.main()
