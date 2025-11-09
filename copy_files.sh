#!/bin/bash
# Script to copy files for reorganization

set -e

BASE=$(pwd)
SRC_PYMATGEN="$BASE/src/pymatgen"

# Copy analysis to pymatgen-analysis
echo "Copying analysis directory..."
cp -r "$SRC_PYMATGEN/analysis" "$BASE/pymatgen-analysis/src/pymatgen/"

# Copy everything else to pymatgen-core
echo "Copying non-analysis modules..."
for item in "$SRC_PYMATGEN"/*; do
    name=$(basename "$item")
    if [ "$name" != "analysis" ] && [ "$name" != "__pycache__" ] && [[ ! "$name" =~ ^\. ]]; then
        echo "  Copying $name"
        cp -r "$item" "$BASE/pymatgen-core/src/pymatgen/"
    fi
done

echo "Done!"
