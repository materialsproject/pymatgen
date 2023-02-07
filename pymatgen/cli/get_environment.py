from __future__ import annotations

import logging
from argparse import ArgumentParser

from pymatgen.analysis.chemenv.utils.chemenv_config import ChemEnvConfig
from pymatgen.analysis.chemenv.utils.defs_utils import chemenv_citations
from pymatgen.analysis.chemenv.utils.scripts_utils import compute_environments

__author__ = "waroquiers"


def main():
    """
    Main function for get_environment CLI.
    """
    parser = ArgumentParser(description="Welcome to the Chemical Environment Package.")
    setup_help = "Used to setup the configuration of the package "
    setup_help += "(MaterialsProject access, ICSD database access, package options, ...)"
    parser.add_argument("-s", "--setup", help=setup_help, action="store_true")
    parser.add_argument(
        "-m",
        "--message-level",
        help="Message level (DEBUG, INFO, WARNING, ERROR or CRITICAL - default : WARNING)",
        default="WARNING",
    )
    args = parser.parse_args()
    if args.setup:
        chemenv_config = ChemEnvConfig.auto_load()
        chemenv_config.setup()
        print("\n Setup completed")
    else:
        chemenv_config = ChemEnvConfig.auto_load()

    # Show welcome message
    print("Chemical Environment package (ChemEnv)")
    print(chemenv_citations())
    print(chemenv_config.package_options_description())

    logging.basicConfig(
        format="%(levelname)s:%(module)s:%(funcName)s:%(message)s",
        level=args.message_level,
    )
    compute_environments(chemenv_config)

    print("Thank you for using the ChemEnv package")
    print(chemenv_citations())


if __name__ == "__main__":
    main()
