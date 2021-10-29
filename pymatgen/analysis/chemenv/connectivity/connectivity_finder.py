"""
Module implementing connectivity finding.
"""

import logging

import numpy as np

from pymatgen.analysis.chemenv.connectivity.structure_connectivity import (
    StructureConnectivity,
)

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "1.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "June 25, 2019"


class ConnectivityFinder:
    """
    Main class used to find the structure connectivity of a structure
    """

    def __init__(self, multiple_environments_choice=None):
        """
        Constructor for the ConnectivityFinder.

        :param multiple_environments_choice: defines the procedure to apply when
        the environment of a given site is described as a "mix" of more than one
        coordination environments.
        """
        self.setup_parameters(multiple_environments_choice=multiple_environments_choice)

    def get_structure_connectivity(self, light_structure_environments):
        """
        Get the structure connectivity from the coordination environments provided
        as an input.

        :param light_structure_environments: LightStructureEnvironments with the
        relevant coordination environments in the structure
        :return: a StructureConnectivity object describing the connectivity of
        the environments in the structure
        """
        logging.info("Setup of structure connectivity graph")
        structure_connectivity = StructureConnectivity(light_structure_environments)
        structure_connectivity.add_sites()
        for isite, site in enumerate(light_structure_environments.structure):
            site_neighbors_sets = light_structure_environments.neighbors_sets[isite]
            if site_neighbors_sets is None:
                continue
            if len(site_neighbors_sets) > 1:
                if self.multiple_environments_choice is None:
                    raise ValueError(
                        "Local environment of site {:d} is a mix and nothing is asked about it".format(isite)
                    )
                if self.multiple_environments_choice == "TAKE_HIGHEST_FRACTION":
                    imax = np.argmax(
                        [ee["ce_fraction"] for ee in light_structure_environments.coordination_environments[isite]]
                    )
                    print("IMAX {:d}".format(imax))
                    site_neighbors_set = site_neighbors_sets[imax]
                else:
                    raise RuntimeError("Should not be here")
            else:
                site_neighbors_set = site_neighbors_sets[0]
            structure_connectivity.add_bonds(isite, site_neighbors_set)
        return structure_connectivity

    def setup_parameters(self, multiple_environments_choice):
        """
        Setup of the parameters for the connectivity finder.
        """
        if multiple_environments_choice is not None:
            if multiple_environments_choice not in ["TAKE_HIGHEST_FRACTION"]:
                raise ValueError(
                    'Option "{}" for multiple_environments_choice is '
                    "not allowed".format(multiple_environments_choice)
                )
        self.multiple_environments_choice = multiple_environments_choice
