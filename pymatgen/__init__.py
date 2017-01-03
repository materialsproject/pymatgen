from __future__ import unicode_literals

import os

__author__ = "Pymatgen Development Team"
__email__ ="pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ ="shyuep@gmail.com"
__date__ = "Jan 2 2017"
__version__ = "4.5.7"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings():
    try:
        import yaml
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.load(f)
    except IOError:
        # If there are any errors, default to using environment variables
        # if present.
        d = {}
        for k, v in os.environ.items():
            if k.startswith("PMG_"):
                d[k] = v
            elif k in ["VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"]:
                d["PMG_" + k] = v
    clean_d = {}
    for k, v in d.items():
        if not k.startswith("PMG_"):
            import warnings
            warnings.warn('With effect from pmg 5.0, all pymatgen settings are'
                          ' prefixed with a "PMG_". E.g., "PMG_VASP_PSP_DIR" '
                          'instead of "VASP_PSP_DIR".')
            clean_d["PMG_" + k] = v
        else:
            clean_d[k] = v
    return clean_d

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
