from __future__ import unicode_literals

__author__ = "Pymatgen Development Team"
__email__ ="pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ ="shyuep@gmail.com"
__date__ = "Jul 8 2016"
__version__ = "4.0.2"


# Useful aliases for commonly used objects and modules.
# Allows from pymatgen import <class> for quick usage.

from .core import *
from .electronic_structure.core import Spin, Orbital
from .matproj.rest import MPRester
from monty.json import MontyEncoder, MontyDecoder, MSONable
