#!/usr/bin/env python

"""
Implementation for get_environment CLI.
"""

import logging
from argparse import ArgumentParser

from pymatgen.analysis.chemenv.utils.chemenv_config import ChemEnvConfig
from pymatgen.analysis.chemenv.utils.scripts_utils import (
    compute_environments,
    thankyou,
    welcome,
)

__author__ = "waroquiers"


def main():
    """
    Main function.
    """
    m_description = "Welcome to the Chemical Environment Package."
    parser = ArgumentParser(description=m_description)
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
    welcome(chemenv_config)
    logging.basicConfig(
        format="%(levelname)s:%(module)s:%(funcName)s:%(message)s",
        level=args.message_level,
    )
    compute_environments(chemenv_config)
    thankyou()


if __name__ == "__main__":
    main()
