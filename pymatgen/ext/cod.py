# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import requests
from pymatgen import Structure
import subprocess
import tabulate
import re

"""
This module provides classes to interface with the Crystallography Open 
Database.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"


class COD(object):
    """
    An interface to the Crystallography Open Database.
    """

    def __init__(self):
        pass

    def get_cod_ids(self, formula):
        """

        :param formula:
        :return:
        """
        query = 'select file from data where formula="Li"'
        r = subprocess.check_output(["mysql", "-u", "cod_reader", "-h",
                                     "www.crystallography.net", "-e",
                                     query, "cod"])
        text = r.decode("utf-8").split("\n")
        cod_ids = []
        for l in text:
            m = re.search("(\d+)", l)
            if m:
                cod_ids.append(int(m.group(1)))
        return cod_ids

    def get_structure_by_id(self, cod_id, **kwargs):
        r = requests.get("http://www.crystallography.net/cod/%s.cif" % cod_id)
        return Structure.from_str(r.text, fmt="cif", **kwargs)

