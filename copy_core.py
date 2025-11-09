#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys

base = "/Users/shyue/repos/pymatgen"
src = f"{base}/src/pymatgen"
dest = f"{base}/pymatgen-core/src/pymatgen"

# Use rsync to copy everything except analysis
result = subprocess.run(
    ["rsync", "-av", "--exclude=analysis", "--exclude=__pycache__", f"{src}/", f"{dest}/"],
    check=False,
    capture_output=True,
    text=True,
)

print(result.stdout)
if result.stderr:
    print(result.stderr, file=sys.stderr)
sys.exit(result.returncode)
