"""Test fixes for pymatgen."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent))
print(sys.path)
TEST_DIR = Path(__file__).absolute().parent
