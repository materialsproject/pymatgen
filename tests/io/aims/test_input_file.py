from __future__ import annotations

from pymatgen.io.aims.sets import AimsInputFile


def test_input_file():
    input_file = AimsInputFile.from_str("testing")
    assert input_file.get_str() == "testing"

    input_file = AimsInputFile.from_string("testing")
    assert input_file.get_string() == "testing"
