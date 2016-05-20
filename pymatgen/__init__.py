from __future__ import unicode_literals

__author__ = "Pymatgen Development Team"
__email__ ="pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ ="shyuep@gmail.com"
__date__ = "May 15 2016"
__version__ = "3.6.0"


# Useful aliases for commonly used objects and modules.
# Allows from pymatgen import <class> for quick usage.

from .core import *
from .serializers.json_coders import pmg_dump, pmg_load
from .electronic_structure.core import Spin, Orbital
from .io.smart import read_structure, write_structure, read_mol, write_mol
from .matproj.rest import MPRester
from monty.json import MontyEncoder, MontyDecoder, MSONable
