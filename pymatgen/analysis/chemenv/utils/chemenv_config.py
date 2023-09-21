"""This module contains the classes for configuration of the chemenv package."""

from __future__ import annotations

import json
from os import makedirs
from os.path import exists, expanduser

from pymatgen.analysis.chemenv.utils.scripts_utils import strategies_class_lookup
from pymatgen.core import SETTINGS

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


class ChemEnvConfig:
    """
    Class used to store the configuration of the chemenv package :
     - Materials project access
     - ICSD database access
     - Default options (strategies, ...).
    """

    DEFAULT_PACKAGE_OPTIONS = dict(
        default_strategy={
            "strategy": "SimplestChemenvStrategy",
            "strategy_options": {
                "distance_cutoff": strategies_class_lookup["SimplestChemenvStrategy"].DEFAULT_DISTANCE_CUTOFF,
                "angle_cutoff": strategies_class_lookup["SimplestChemenvStrategy"].DEFAULT_ANGLE_CUTOFF,
                "additional_condition": strategies_class_lookup["SimplestChemenvStrategy"].DEFAULT_ADDITIONAL_CONDITION,
                "continuous_symmetry_measure_cutoff": strategies_class_lookup[
                    "SimplestChemenvStrategy"
                ].DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF,
            },
        },
        default_max_distance_factor=1.5,
    )

    def __init__(self, package_options=None):
        """:param package_options:"""
        if SETTINGS.get("PMG_MAPI_KEY"):
            self.materials_project_configuration = SETTINGS.get("PMG_MAPI_KEY")
        else:
            self.materials_project_configuration = None

        if package_options is None:
            self.package_options = self.DEFAULT_PACKAGE_OPTIONS
        else:
            self.package_options = package_options

    def setup(self):
        """Setup the class."""
        while True:
            print("\n=> Configuration of the ChemEnv package <=")
            print("Current configuration :")
            if self.has_materials_project_access:
                print(" - Access to materials project is configured (add test ?)")
            else:
                print(" - No access to materials project")
            print(" - Package options :")
            for key, val in self.package_options.items():
                print(f"     {key}   :   {val}")
            print("\nChoose in the following :")
            print(" <1> + <ENTER> : configuration of the package options (strategy, ...)")
            print(" <q> + <ENTER> : quit without saving configuration")
            test = input(" <S> + <ENTER> : save configuration and quit\n ... ")
            if test == "1":
                self.setup_package_options()
            elif test == "q":
                break
            elif test == "S":
                config_file = self.save()
                break
            else:
                print(" ... wrong key, try again ...")
            print()
        if test == "S":
            print(f"Configuration has been saved to file {config_file!r}")

    @property
    def has_materials_project_access(self):
        """
        Whether MP access is enabled.
        """
        return self.materials_project_configuration is not None

    def setup_package_options(self):
        """Setup the package options."""
        self.package_options = self.DEFAULT_PACKAGE_OPTIONS
        print("Choose between the following strategies : ")
        strategies = list(strategies_class_lookup)
        for idx, strategy in enumerate(strategies, 1):
            print(f" <{idx}> : {strategy}")
        test = input(" ... ")
        self.package_options["default_strategy"] = {
            "strategy": strategies[int(test) - 1],
            "strategy_options": {},
        }
        strategy_class = strategies_class_lookup[strategies[int(test) - 1]]
        if len(strategy_class.STRATEGY_OPTIONS) > 0:
            for option, option_dict in strategy_class.STRATEGY_OPTIONS.items():
                while True:
                    print(f"  => Enter value for option {option!r} (<ENTER> for default = {option_dict['default']})\n")
                    print("     Valid options are :\n")
                    print(f"       {option_dict['type'].allowed_values}")
                    test = input("     Your choice : ")
                    if test == "":
                        self.package_options["default_strategy"]["strategy_options"][option] = option_dict["type"](
                            strategy_class.STRATEGY_OPTIONS[option]["default"]
                        )
                        break
                    try:
                        self.package_options["default_strategy"]["strategy_options"][option] = option_dict["type"](test)
                        break
                    except ValueError:
                        print(f"Wrong input for {option=}")

    def package_options_description(self):
        """Describe package options."""
        out = "Package options :\n"
        out += f" - Maximum distance factor : {self.package_options['default_max_distance_factor']:.4f}\n"
        out += f" - Default strategy is \"{self.package_options['default_strategy']['strategy']}\" :\n"
        strategy_class = strategies_class_lookup[self.package_options["default_strategy"]["strategy"]]
        out += f"{strategy_class.STRATEGY_DESCRIPTION}\n"
        out += "   with options :\n"
        for option in strategy_class.STRATEGY_OPTIONS:
            out += f"     - {option} : {self.package_options['default_strategy']['strategy_options'][option]}\n"
        return out

    def save(self, root_dir=None):
        """
        Save the options.
        :param root_dir:
        """
        if root_dir is None:
            home = expanduser("~")
            root_dir = f"{home}/.chemenv"
        if not exists(root_dir):
            makedirs(root_dir)
        config_dict = {"package_options": self.package_options}
        config_file = f"{root_dir}/config.json"
        if exists(config_file):
            test = input("Overwrite existing configuration ? (<Y> + <ENTER> to confirm)")
            if test != "Y":
                print("Configuration not saved")
                return config_file
        with open(config_file, "w") as f:
            json.dump(config_dict, f)
        print("Configuration saved")
        return config_file

    @classmethod
    def auto_load(cls, root_dir=None):
        """
        Autoload options.
        :param root_dir:
        """
        if root_dir is None:
            home = expanduser("~")
            root_dir = f"{home}/.chemenv"
        config_file = f"{root_dir}/config.json"
        try:
            with open(config_file) as f:
                config_dict = json.load(f)
            return ChemEnvConfig(package_options=config_dict["package_options"])

        except OSError:
            print(f"Unable to load configuration from file {config_file!r} ...")
            print(" ... loading default configuration")
            return ChemEnvConfig()
