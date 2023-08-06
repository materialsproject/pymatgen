from __future__ import annotations

from pymatgen.analysis.chemenv.utils.chemenv_config import ChemEnvConfig
from pymatgen.core import SETTINGS
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "waroquiers"

config_file_dir = f"{TEST_FILES_DIR}/chemenv/config"


class TestChemenvConfig(PymatgenTest):
    def test_chemenv_config(self):
        config = ChemEnvConfig()

        assert config.has_materials_project_access == bool(SETTINGS.get("PMG_MAPI_KEY"))

        package_options = {**ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS, "default_max_distance_factor": 1.8}

        config = ChemEnvConfig(package_options=package_options)

        assert (
            config.package_options_description() == "Package options :\n"
            " - Maximum distance factor : 1.8000\n"
            ' - Default strategy is "SimplestChemenvStrategy" :\n'
            "    Simplest ChemenvStrategy using fixed angle and distance parameters \n"
            "    for the definition of neighbors in the Voronoi approach. \n"
            "    The coordination environment is then given as the one with the \n"
            "    lowest continuous symmetry measure.\n"
            "   with options :\n"
            "     - distance_cutoff : 1.4\n"
            "     - angle_cutoff : 0.3\n"
            "     - additional_condition : 1\n"
            "     - continuous_symmetry_measure_cutoff : 10\n"
        )

        config.save(root_dir=self.tmp_path)

        config = config.auto_load(root_dir=self.tmp_path)

        assert config.package_options == package_options
