from __future__ import annotations

from pathlib import Path

import pytest

from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice
from pymatgen.util.testing import TEST_FILES_DIR

ex_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "example_files"

ex_slice1_fname = ex_files_dir / "ex_out_slice_latmin"
ex_slice2_fname = ex_files_dir / "ex_out_slice_ionmin"
with open(ex_slice1_fname) as f:
    ex_slice1 = list.copy(list(f))
with open(ex_slice2_fname) as f:
    ex_slice2 = list.copy(list(f))


@pytest.mark.parametrize(
    ("out_slice", "varname"),
    [
        (ex_slice1, "structure"),
        (ex_slice1, "eiter_type"),
        (ex_slice1, "elecmindata"),
        (ex_slice1, "stress"),
        (ex_slice1, "strain"),
        (ex_slice1, "iter"),
        (ex_slice1, "e"),
        (ex_slice1, "grad_k"),
        (ex_slice1, "alpha"),
        (ex_slice1, "linmin"),
        (ex_slice1, "abs_magneticmoment"),
        (ex_slice1, "tot_magneticmoment"),
        (ex_slice1, "mu"),
        (ex_slice1, "elec_iter"),
        (ex_slice1, "elec_e"),
        (ex_slice1, "elec_grad_k"),
        (ex_slice1, "elec_alpha"),
        (ex_slice1, "elec_linmin"),
    ],
)
def test_jdftxoutfileslice_has_1layer_jstrucs_freakout(out_slice: list[str], varname: str):
    joutslice = JDFTXOutfileSlice.from_out_slice(out_slice)
    getattr(joutslice, varname)  # No freakout here
    joutslice.jstrucs = None
    with pytest.raises(AttributeError):
        getattr(joutslice, varname)  # Freakout here


# @pytest.mark.parametrize(
#     ("out_slice", "varname"),
#     [
#         (ex_text_slice, "e"),
#     ]
# )
# def jeiters_has_2layer_slice_freakout(out_slice: list[str], varname: str):
#     jstrucs = JOutStructures.from_out_slice(out_slice)
#     getattr(jstrucs, varname) # No freakout here
#     setattr(jstrucs.slices[-1], varname, None)
#     with pytest.raises(ValueError):
#         jstrucs.iter # Freakout here
#     # Reset
#     jstrucs = JOutStructures.from_out_slice(out_slice)
#     getattr(jstrucs, varname) # No freakout here
#     jstrucs.slices = []
#     with pytest.raises(ValueError):
#         jstrucs.iter # Freakout here
