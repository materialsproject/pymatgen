"""
Backward compatibility stub for pymatgen.analysis.structure_matcher.

This module has been moved to pymatgen.core.structure_matcher.
This stub imports everything from the new location to maintain backward compatibility.

Note: This module requires pymatgen-core to be installed. Since pymatgen-analysis
depends on pymatgen-core, this should always be available.
"""

from __future__ import annotations

try:
    # Import from the new location in pymatgen.core
    from pymatgen.core.structure_matcher import (  # type: ignore[import-untyped]
        AbstractComparator,
        ElementComparator,
        FrameworkComparator,
        OccupancyComparator,
        OrderDisorderElementComparator,
        SiteOrderedIStructure,
        SpeciesComparator,
        SpinComparator,
        StructureMatcher,
        get_linear_assignment_solution,
    )
except ImportError as e:
    raise ImportError(
        "pymatgen.analysis.structure_matcher has been moved to pymatgen.core.structure_matcher. "
        "Please install pymatgen-core or use: from pymatgen.core.structure_matcher import ..."
    ) from e

__all__ = [
    "AbstractComparator",
    "ElementComparator",
    "FrameworkComparator",
    "OccupancyComparator",
    "OrderDisorderElementComparator",
    "SiteOrderedIStructure",
    "SpeciesComparator",
    "SpinComparator",
    "StructureMatcher",
    "get_linear_assignment_solution",
]
