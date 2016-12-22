import unittest

import os
import yaml
from pymatgen import SETTINGS_FILE, SETTINGS


class SettingsTestCase(unittest.TestCase):

    def test_something(self):
        if SETTINGS_FILE.exists():
            with SETTINGS_FILE.open("rt") as f:
                d = yaml.load(f)
                for k, v in d.items():
                    self.assertEqual(v, SETTINGS[k])
        else:
            for k, v in SETTINGS.items():
                self.assertEqual(v, os.environ.get(k))


if __name__ == '__main__':
    unittest.main()
