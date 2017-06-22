# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import requests
from pymatgen import Structure

"""
This module provides classes to interface with the Crystallography Open 
Database.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"


def get_structure(cod_id, **kwargs):
    cif = requests.get("http://www.crystallography.net/cod/%s.cif" % cod_id)
    return Structure.from_str(cif.text, fmt="cif", **kwargs)