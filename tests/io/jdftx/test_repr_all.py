from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pytest

from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice
from pymatgen.io.jdftx.jelstep import JElStep, JElSteps
from pymatgen.io.jdftx.joutstructures import JOutStructure, JOutStructures
from pymatgen.io.jdftx.outputs import JDFTXOutfile

from .conftest import ex_jstep_lines1, ex_jstep_lines2, ex_jstruc_slice1, ex_outfileslice1, example_sp_outfile_path

if TYPE_CHECKING:
    from collections.abc import Callable


@pytest.mark.parametrize(
    ("init_meth", "init_var", "add_checks"),
    [
        (JDFTXOutfile, example_sp_outfile_path, lambda dir_repr: None),
        (JDFTXOutfileSlice._from_out_slice, ex_outfileslice1, lambda dir_repr: None),
        (lambda x: JOutStructures._from_out_slice(x, opt_type="lattice"), ex_outfileslice1, lambda dir_repr: None),
        (lambda x: JOutStructure._from_text_slice(x, opt_type="lattice"), ex_jstruc_slice1, lambda dir_repr: None),
        (lambda x: JElStep._from_lines_collect(x, "ElecMinimize", "F"), ex_jstep_lines1, lambda dir_repr: None),
        (
            lambda x: JElSteps._from_text_slice(x, opt_type="ElecMinimize", etype="F"),
            [line for exl in [ex_jstep_lines1, ex_jstep_lines2] for line in exl],
            lambda dir_repr: None,
        ),
    ],
)
def test_dir_repr(init_meth: Callable, init_var: Any, add_checks: Callable) -> None:
    dirrable = init_meth(init_var)
    dir_repr = dir(dirrable)
    add_checks(dir_repr)


@pytest.mark.parametrize(
    ("init_meth", "init_var", "add_checks"),
    [
        (JDFTXOutfile, example_sp_outfile_path, lambda dir_repr: None),
        (JDFTXOutfileSlice._from_out_slice, ex_outfileslice1, lambda dir_repr: None),
        (lambda x: JOutStructures._from_out_slice(x, opt_type="lattice"), ex_outfileslice1, lambda dir_repr: None),
        (lambda x: JOutStructure._from_text_slice(x, opt_type="lattice"), ex_jstruc_slice1, lambda dir_repr: None),
        (lambda x: JElStep._from_lines_collect(x, "ElecMinimize", "F"), ex_jstep_lines1, lambda dir_repr: None),
        (
            lambda x: JElSteps._from_text_slice(x, opt_type="ElecMinimize", etype="F"),
            [line for exl in [ex_jstep_lines1, ex_jstep_lines2] for line in exl],
            lambda dir_repr: None,
        ),
    ],
)
def test_repr_repr(init_meth: Callable, init_var: Any, add_checks: Callable) -> None:
    reprrable = init_meth(init_var)
    repr_repr = repr(reprrable)
    add_checks(repr_repr)


@pytest.mark.parametrize(
    ("init_meth", "init_var", "add_checks"),
    [
        (JDFTXOutfile, example_sp_outfile_path, lambda dir_repr: None),
        (JDFTXOutfileSlice._from_out_slice, ex_outfileslice1, lambda dir_repr: None),
        (lambda x: JOutStructures._from_out_slice(x, opt_type="lattice"), ex_outfileslice1, lambda dir_repr: None),
        (lambda x: JOutStructure._from_text_slice(x, opt_type="lattice"), ex_jstruc_slice1, lambda dir_repr: None),
        (lambda x: JElStep.from_lines_collect(x, "ElecMinimize", "F"), ex_jstep_lines1, lambda dir_repr: None),
        (
            lambda x: JElSteps._from_text_slice(x, opt_type="ElecMinimize", etype="F"),
            [line for exl in [ex_jstep_lines1, ex_jstep_lines2] for line in exl],
            lambda dir_repr: None,
        ),
    ],
)
def test_str_repr(init_meth: Callable, init_var: Any, add_checks: Callable) -> None:
    strrable = init_meth(init_var)
    str_repr = str(strrable)
    add_checks(str_repr)
