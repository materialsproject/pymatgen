from __future__ import unicode_literals

__author__ = "Pymatgen Development Team"
__email__ ="pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ ="shyuep@gmail.com"
__date__ = "Jul 14 2016"
__version__ = "4.1.0"


# Order of imports is important on some systems to avoid 
# failures when loading shared libraries.
import spglib
from . import optimization, util
del(spglib, optimization, util)

# Useful aliases for commonly used objects and modules.
# Allows from pymatgen import <class> for quick usage.

from .core import *
from .electronic_structure.core import Spin, Orbital
from .matproj.rest import MPRester
from monty.json import MontyEncoder, MontyDecoder, MSONable
