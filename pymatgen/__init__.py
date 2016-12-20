from __future__ import unicode_literals

import os
try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

__author__ = "Pymatgen Development Team"
__email__ ="pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ ="shyuep@gmail.com"
__date__ = "Dec 20 2016"
__version__ = "4.5.5"


SETTINGS_FILE = Path("~/.pmgrc.yaml").expanduser()


def _load_pmg_settings():
    print(SETTINGS_FILE.exists())
    if SETTINGS_FILE.exists():
        try:
            import yaml
            with SETTINGS_FILE.open("rt") as f:
                return yaml.load(f)
        except Exception as ex:
            # If there are any errors, default to using environment variables
            # if present.
            pass
    d = {}
    for k in ["VASP_PSP_DIR", "MAPI_KEY"]:
        d[k] = os.environ.get(k)
    return d

SETTINGS = _load_pmg_settings()

# Order of imports is important on some systems to avoid
# failures when loading shared libraries.
# import spglib
# from . import optimization, util
# del(spglib, optimization, util)

# Useful aliases for commonly used objects and modules.
# Allows from pymatgen import <class> for quick usage.

from pymatgen.core import *
from .electronic_structure.core import Spin, Orbital
from .matproj.rest import MPRester
from monty.json import MontyEncoder, MontyDecoder, MSONable
