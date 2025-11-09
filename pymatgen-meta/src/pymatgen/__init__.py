"""
Pymatgen meta-package.

This package automatically imports and re-exports from pymatgen-core, pymatgen-analysis,
and pymatgen-apps to maintain backward compatibility.
"""

from __future__ import annotations

# Import everything from pymatgen-core
from pymatgen.core import *  # noqa: F403, F401
from pymatgen.core import __version__  # noqa: F401

# Import analysis modules
try:
    from pymatgen.analysis import *  # noqa: F403, F401
except ImportError:
    # pymatgen-analysis not installed
    pass

# Import apps and cli modules
try:
    from pymatgen.apps import *  # noqa: F403, F401
except ImportError:
    # pymatgen-apps not installed
    pass

try:
    from pymatgen.cli import *  # noqa: F403, F401
except ImportError:
    # pymatgen-apps not installed
    pass

__all__ = [
    "__version__",
]

