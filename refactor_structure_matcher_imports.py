#!/usr/bin/env python3
"""
Script to refactor StructureMatcher imports in pymatgen-core files.

This script updates imports from pymatgen.analysis.structure_matcher to
pymatgen.core.structure_matcher in files that will be part of pymatgen-core.
"""

from __future__ import annotations

import re
from pathlib import Path

# Files in pymatgen-core that need to be updated
# These are relative to src/pymatgen/
FILES_TO_UPDATE = [
    "core/tensors.py",
    "core/surface.py",
    "core/structure.py",
    "io/vasp/sets.py",
    "entries/mixing_scheme.py",
    "entries/entry_tools.py",
    "transformations/standard_transformations.py",
    "phonon/thermal_displacements.py",
]

# Import patterns to replace
REPLACEMENTS = [
    # Simple StructureMatcher import
    (
        r"from pymatgen\.analysis\.structure_matcher import StructureMatcher",
        "from pymatgen.core.structure_matcher import StructureMatcher",
    ),
    # Multiple imports (SpeciesComparator, StructureMatcher, etc.)
    (
        r"from pymatgen\.analysis\.structure_matcher import ([A-Za-z, ]+)",
        lambda m: f"from pymatgen.core.structure_matcher import {m.group(1)}",
    ),
]


def update_file(file_path: Path) -> bool:
    """Update a file with new imports. Returns True if file was modified."""
    if not file_path.exists():
        print(f"Warning: {file_path} does not exist, skipping")
        return False

    content = file_path.read_text()
    modified = False

    # Apply replacements
    for pattern, replacement in REPLACEMENTS:
        if callable(replacement):
            # For regex with groups, use a function
            def replacer(match):
                return replacement(match)

            new_content = re.sub(pattern, replacer, content)
        else:
            new_content = re.sub(pattern, replacement, content)

        if new_content != content:
            content = new_content
            modified = True

    if modified:
        file_path.write_text(content)
        print(f"  ✓ Updated {file_path}")
        return True
    print(f"  - No changes needed in {file_path}")
    return False


def main():
    base = Path(__file__).parent
    src_pymatgen = base / "src" / "pymatgen"

    if not src_pymatgen.exists():
        print(f"Error: {src_pymatgen} does not exist")
        print("Note: This script should be run before reorganization, or update paths to pymatgen-core")
        return 1

    print("Refactoring StructureMatcher imports in pymatgen-core files...")
    print(f"Base directory: {src_pymatgen}\n")

    updated_count = 0
    for rel_path in FILES_TO_UPDATE:
        file_path = src_pymatgen / rel_path
        if update_file(file_path):
            updated_count += 1

    print(f"\n✓ Refactoring complete! Updated {updated_count} files.")
    print("\nNote: After reorganization, these files will be in pymatgen-core")
    print("      and will use pymatgen.core.structure_matcher directly.")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
