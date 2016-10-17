from __future__ import unicode_literals

import os
import warnings

__author__ = "Pymatgen Development Team"
__email__ ="pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ ="shyuep@gmail.com"
__date__ = "Oct 17 2016"
__version__ = "4.4.5"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings():
    if os.path.exists(SETTINGS_FILE):
        try:
            from monty.serialization import loadfn
            return loadfn(SETTINGS_FILE)
        except:
            # If there are any errors, default to using environment variables
            # if present.
            pass
    d = {}
    for k in ["VASP_PSP_DIR", "MAPI_KEY"]:
        if k in os.environ:
            warnings.warn("You have %s set in the env. From pmg 5, all "
                          "settings should be in the .pmgrc.yaml file." % k)
        d[k] = os.environ.get(k)
    return d

SETTINGS = _load_pmg_settings()

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
