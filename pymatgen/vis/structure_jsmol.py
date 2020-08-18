"""
Visualization for structures using jupyter-jsmol.
"""

import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from monty.dev import requires

try:
    from jupyter_jsmol import JsmolView
    jupyter_jsmol_loaded = True
except ImportError:
    jupyter_jsmol_loaded = False


@requires(jupyter_jsmol_loaded, "To use quick_view, you need to have jupyter_jsmol installed.")
def quick_view(structure, conventional=False, transform=None, show_box=True):
    """
    A function to visualize pymatgen Structure objects in jupyter notebook using chemview package.

    Args:
        structure: pymatgen Structure
        conventional: (bool) use conventional cell. Defaults to False.
        transform: (list) can be used to make supercells with pymatgen.Structure.make_supercell method
    Returns:
        A JsmolView object.
    """

    s = structure.copy()
    if conventional:
        s = SpacegroupAnalyzer(s).get_conventional_standard_structure()

    if transform:
        s.make_supercell(transform)

    return JsmolView.from_str(s.to('cif'))


def close_all():
    JsmolView.close_all()
