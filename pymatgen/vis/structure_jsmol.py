"""Visualization for structures using jupyter-jsmol.
"""

try:
    from jupyter_jsmol import JsmolView
    jupyter_jsmol_loaded = True
except ImportError:
    jupyter_jsmol_loaded = False

from monty.dev import requires

from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


@requires(jupyter_jsmol_loaded, "To use quick_view, you need to have jupyter_jsmol installed.")
def quick_view(structure: Structure, conventional: bool = False, transform: list = None) -> JsmolView:
    """A function to visualize pymatgen Structure objects in jupyter notebook using jupyter_jsmol package.

    Args:
        structure: pymatgen Structure object.
        conventional: use conventional cell. Defaults to False.
        transform: can be used to make supercells with pymatgen.Structure.make_supercell method.

    Returns:
        A jupyter widget object.
    """

    s = structure.copy()
    if conventional:
        s = SpacegroupAnalyzer(s).get_conventional_standard_structure()

    if transform:
        s.make_supercell(transform)

    return JsmolView.from_str(s.to('cif'))
