#!/usr/bin/env python3
"""
Script to create a backward-compatibility stub for structure_matcher in analysis.

This creates a stub file in pymatgen-analysis that imports from pymatgen.core.structure_matcher
to maintain backward compatibility.
"""

from __future__ import annotations

from pathlib import Path

STUB_CONTENT = '''"""
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
'''


def main():
    base = Path(__file__).parent
    stub_path = base / "pymatgen-analysis" / "src" / "pymatgen" / "analysis" / "structure_matcher.py"

    # Create directory if it doesn't exist
    stub_path.parent.mkdir(parents=True, exist_ok=True)

    # Write the stub file
    stub_path.write_text(STUB_CONTENT)
    print(f"Created backward compatibility stub at {stub_path}")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
