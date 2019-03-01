
import sys
import os
import warnings
import ruamel.yaml as yaml
from fnmatch import fnmatch

__author__ = "Pymatgen Development Team"
__email__ ="pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ ="shyuep@gmail.com"
__version__ = "2019.2.28"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings():
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
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
from .ext.matproj import MPRester
from monty.json import MontyEncoder, MontyDecoder, MSONable


def get_structure_from_mp(formula):
    """
    Convenience method to get a crystal from the Materials Project database via
    the API. Requires PMG_MAPI_KEY to be set.

    Args:
        formula (str): A formula

    Returns:
        (Structure) The lowest energy structure in Materials Project with that
            formula.
    """
    m = MPRester()
    entries = m.get_entries(formula, inc_structure="final")
    if len(entries) == 0:
        raise ValueError("No structure with formula %s in Materials Project!" %
                         formula)
    elif len(entries) > 1:
        warnings.warn("%d structures with formula %s found in Materials "
                      "Project. The lowest energy structure will be returned." %
                      (len(entries), formula))
    return min(entries, key=lambda e: e.energy_per_atom).structure


if sys.version_info < (3, 5):
    warnings.warn("""
Pymatgen will drop Py2k support from v2019.1.1. Pls consult the documentation
at https://www.pymatgen.org for more details.""")


def loadfn(fname):
    """
    Convenience method to perform quick loading of data from a filename. The
    type of object returned depends the file type.

    Args:
        fname (string): A filename.

    Returns:
        Note that fname is matched using unix-style, i.e., fnmatch.
        (Structure) if *POSCAR*/*CONTCAR*/*.cif
        (Vasprun) *vasprun*
        (obj) if *json* (passthrough to monty.serialization.loadfn)
    """
    if (fnmatch(fname, "*POSCAR*") or fnmatch(fname, "*CONTCAR*") or
            ".cif" in fname.lower()) or fnmatch(fname, "*.vasp"):
        return Structure.from_file(fname)
    elif fnmatch(fname, "*vasprun*"):
        from pymatgen.io.vasp import Vasprun
        return Vasprun(fname)
    elif fnmatch(fname, "*.json*"):
        from monty.serialization import loadfn
        return loadfn(fname)