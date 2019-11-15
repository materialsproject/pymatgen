#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest
import os
import shutil
from pymatgen.analysis.chemenv.utils.chemenv_config import ChemEnvConfig
from pymatgen import SETTINGS
from monty.tempfile import ScratchDir

config_file_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                               'test_files', "chemenv", "config")


class ChemenvConfigTest(unittest.TestCase):

    def test_chemenv_config(self):
        with ScratchDir("."):
            config = ChemEnvConfig()

            if SETTINGS.get("PMG_MAPI_KEY", "") != "":
                self.assertTrue(config.has_materials_project_access)
            else:
                self.assertFalse(config.has_materials_project_access)

            package_options = ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS
            package_options['default_max_distance_factor'] = 1.8

            config = ChemEnvConfig(package_options=package_options)

            self.assertEqual(config.package_options_description(),
                             'Package options :\n'
                             ' - Maximum distance factor : 1.8000\n'
                             ' - Default strategy is "SimplestChemenvStrategy" :\n'
                             '    Simplest ChemenvStrategy using fixed angle and distance parameters \n'
                             '    for the definition of neighbors in the Voronoi approach. \n'
                             '    The coordination environment is then given as the one with the \n'
                             '    lowest continuous symmetry measure.\n'
                             '   with options :\n'
                             '     - distance_cutoff : 1.4\n'
                             '     - angle_cutoff : 0.3\n'
                             '     - additional_condition : 1\n'
                             '     - continuous_symmetry_measure_cutoff : 10.0\n')

            config.save(root_dir='tmp_dir')

            config = config.auto_load(root_dir='tmp_dir')

            self.assertEqual(config.package_options, package_options)


if __name__ == "__main__":
    unittest.main()
