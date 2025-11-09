#!/usr/bin/env python3
"""Script to reorganize pymatgen into mono-repo structure."""

import os
import shutil
import sys
from pathlib import Path

def main():
    # Base paths
    base = Path(__file__).parent
    src_pymatgen = base / "src" / "pymatgen"

    if not src_pymatgen.exists():
        print(f"Error: {src_pymatgen} does not exist", file=sys.stderr)
        sys.exit(1)

    # Create new structure
    pymatgen_core = base / "pymatgen-core" / "src" / "pymatgen"
    pymatgen_analysis = base / "pymatgen-analysis" / "src" / "pymatgen" / "analysis"
    pymatgen_apps = base / "pymatgen-apps" / "src" / "pymatgen"
    pymatgen_meta = base / "pymatgen-meta" / "src" / "pymatgen"

    # Create directories
    pymatgen_core.mkdir(parents=True, exist_ok=True)
    pymatgen_analysis.parent.mkdir(parents=True, exist_ok=True)
    pymatgen_apps.mkdir(parents=True, exist_ok=True)
    pymatgen_meta.mkdir(parents=True, exist_ok=True)

    print("Created directory structure")

    # Copy analysis to pymatgen-analysis (skip if already exists)
    # But exclude structure_matcher.py as it will be moved to core
    analysis_src = src_pymatgen / "analysis"
    if analysis_src.exists() and not pymatgen_analysis.exists():
        print(f"Copying {analysis_src} to {pymatgen_analysis} (excluding structure_matcher.py)")
        # Create the directory first
        pymatgen_analysis.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(analysis_src, pymatgen_analysis, ignore=shutil.ignore_patterns("structure_matcher.py"))
        print("Copied analysis directory")
    elif pymatgen_analysis.exists():
        print("Analysis directory already exists, skipping")
    else:
        print(f"Warning: {analysis_src} does not exist")
    
    # Move structure_matcher.py from analysis to core/core directory
    structure_matcher_src = src_pymatgen / "analysis" / "structure_matcher.py"
    structure_matcher_dest = pymatgen_core / "core" / "structure_matcher.py"
    if structure_matcher_src.exists():
        print(f"Moving {structure_matcher_src} to {structure_matcher_dest}")
        structure_matcher_dest.parent.mkdir(parents=True, exist_ok=True)
        if structure_matcher_dest.exists():
            structure_matcher_dest.unlink()
        shutil.copy2(structure_matcher_src, structure_matcher_dest)
        print("Moved structure_matcher.py to pymatgen.core")

    # Copy apps and cli to pymatgen-apps
    apps_src = src_pymatgen / "apps"
    cli_src = src_pymatgen / "cli"
    
    if apps_src.exists():
        apps_dest = pymatgen_apps / "apps"
        if apps_dest.exists():
            shutil.rmtree(apps_dest)
        print(f"Copying {apps_src} to {apps_dest}")
        shutil.copytree(apps_src, apps_dest)
        print("Copied apps directory")
    
    if cli_src.exists():
        cli_dest = pymatgen_apps / "cli"
        if cli_dest.exists():
            shutil.rmtree(cli_dest)
        print(f"Copying {cli_src} to {cli_dest}")
        shutil.copytree(cli_src, cli_dest)
        print("Copied cli directory")

    # Copy everything else to pymatgen-core (excluding analysis, apps, cli)
    print("Copying non-analysis/apps/cli modules to pymatgen-core...")
    copied = []
    excluded = {"analysis", "apps", "cli", "__pycache__"}
    for item in sorted(src_pymatgen.iterdir()):
        if item.name in excluded:
            continue
        if item.name.startswith("."):
            continue
        dest = pymatgen_core / item.name
        if item.is_dir():
            print(f"  Copying {item.name}/")
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(item, dest)
            copied.append(item.name)
        else:
            print(f"  Copying {item.name}")
            shutil.copy2(item, dest)
            copied.append(item.name)

    print(f"\nReorganization complete! Copied {len(copied)} items to pymatgen-core")
    return 0

if __name__ == "__main__":
    sys.exit(main())

