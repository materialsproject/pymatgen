#!/usr/bin/env python3
"""
Script to make imports from pymatgen.analysis optional in pymatgen-core.

This script updates files in pymatgen-core to use try/except blocks for
imports from pymatgen.analysis, allowing pymatgen-core to work without
pymatgen-analysis installed.
"""

import re
from pathlib import Path

# Files that need to be updated
FILES_TO_UPDATE = {
    "pymatgen-core/src/pymatgen/core/tensors.py": [
        (r"from pymatgen\.analysis\.structure_matcher import StructureMatcher", 
         "try:\n    from pymatgen.analysis.structure_matcher import StructureMatcher\nexcept ImportError:\n    StructureMatcher = None  # type: ignore"),
    ],
    "pymatgen-core/src/pymatgen/core/interface.py": [
        (r"from pymatgen\.analysis\.adsorption import AdsorbateSiteFinder",
         "try:\n    from pymatgen.analysis.adsorption import AdsorbateSiteFinder\nexcept ImportError:\n    AdsorbateSiteFinder = None  # type: ignore"),
    ],
    "pymatgen-core/src/pymatgen/core/surface.py": [
        (r"from pymatgen\.analysis\.structure_matcher import StructureMatcher",
         "try:\n    from pymatgen.analysis.structure_matcher import StructureMatcher\nexcept ImportError:\n    StructureMatcher = None  # type: ignore"),
    ],
    "pymatgen-core/src/pymatgen/transformations/standard_transformations.py": [
        (r"from pymatgen\.analysis\.bond_valence import BVAnalyzer",
         "try:\n    from pymatgen.analysis.bond_valence import BVAnalyzer\nexcept ImportError:\n    BVAnalyzer = None  # type: ignore"),
        (r"from pymatgen\.analysis\.elasticity\.strain import Deformation",
         "try:\n    from pymatgen.analysis.elasticity.strain import Deformation\nexcept ImportError:\n    Deformation = None  # type: ignore"),
        (r"from pymatgen\.analysis\.ewald import EwaldMinimizer, EwaldSummation",
         "try:\n    from pymatgen.analysis.ewald import EwaldMinimizer, EwaldSummation\nexcept ImportError:\n    EwaldMinimizer = None  # type: ignore\n    EwaldSummation = None  # type: ignore"),
        (r"from pymatgen\.analysis\.structure_matcher import StructureMatcher",
         "try:\n    from pymatgen.analysis.structure_matcher import StructureMatcher\nexcept ImportError:\n    StructureMatcher = None  # type: ignore"),
    ],
    "pymatgen-core/src/pymatgen/entries/compatibility.py": [
        (r"from pymatgen\.analysis\.structure_analyzer import oxide_type, sulfide_type",
         "try:\n    from pymatgen.analysis.structure_analyzer import oxide_type, sulfide_type\nexcept ImportError:\n    oxide_type = None  # type: ignore\n    sulfide_type = None  # type: ignore"),
    ],
}

def update_file(file_path: Path, replacements: list[tuple[str, str]]) -> None:
    """Update a file with optional imports."""
    if not file_path.exists():
        print(f"Warning: {file_path} does not exist, skipping")
        return
    
    content = file_path.read_text()
    modified = False
    
    for pattern, replacement in replacements:
        if re.search(pattern, content):
            content = re.sub(pattern, replacement, content)
            modified = True
            print(f"  Updated import in {file_path.name}")
    
    if modified:
        file_path.write_text(content)
        print(f"  ✓ Updated {file_path}")

def main():
    base = Path(__file__).parent
    
    print("Making imports from pymatgen.analysis optional...")
    for rel_path, replacements in FILES_TO_UPDATE.items():
        file_path = base / rel_path
        update_file(file_path, replacements)
    
    print("\nDone! Note: You may need to add runtime checks for None values where these imports are used.")

if __name__ == "__main__":
    main()

