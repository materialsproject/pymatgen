# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import re
import logging
import os

from monty.io import zopen
from monty.json import MSONable
from pymatgen.core import Molecule
from pymatgen.io.qchem_io.inputs import QCInput

from .utils import read_table_pattern, read_pattern
"""
Classes for reading/manipulating/writing QChem ouput files.
"""

__author__ = "Samuel Blau, Brandon Woods, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)


class QChemDictSet(QCInput):
    """
    Args:

    """

    def __init__(self, molecule, ):
        


class OptSet(QChemDictSet):
    """
    QChemDictSet for a geometry optimization
    """

    CONFIG = loadfn(os.path.join(MODULE_DIR, "OptSet.yaml"))

    def __init__(self, ):
        """
        Args:

        """
        super(OptSet, self).__init__()



class SinglePointSet(QChemDictSet):
    """
    QChemDictSet for a single point calculation
    """

    CONFIG = loadfn(os.path.join(MODULE_DIR, "SinglePointSet.yaml"))

    def __init__(self, ):
        """
        Args:

        """
        super(OptSet, self).__init__()