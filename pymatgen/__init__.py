from __future__ import unicode_literals

__author__ = ", ".join(["Shyue Ping Ong", "Anubhav Jain", "Geoffroy Hautier",
                        "William Davidson Richard", "Stephen Dacek",
                        "Sai Jayaraman", "Michael Kocher", "Dan Gunter",
                        "Shreyas Cholia", "Vincent L Chevrier",
                        "Rickard Armiento"])
__date__ = "Sep 23 2014"
__version__ = "3.0.5"

#Useful aliases for commonly used objects and modules.

from .core import *
from .serializers.json_coders import pmg_dump, pmg_load
from .electronic_structure.core import Spin, Orbital
from .io.smartio import read_structure, write_structure, read_mol, write_mol
from .matproj.rest import MPRester
from monty.json import MontyEncoder, MontyDecoder, MSONable
