"""This module tests the __dir__, __repr__, __str__, and to_dict methods of the objects in the outputs module.

The tests test that the desired objects do not raise errors upon calling the methods, and any additional checks
that are needed for the objects to pass the tests.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pytest

from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice
from pymatgen.io.jdftx.jelstep import JElStep, JElSteps
from pymatgen.io.jdftx.joutstructures import JOutStructure, JOutStructures
from pymatgen.io.jdftx.outputs import JDFTXOutfile

from .outputs_test_utils import (
    ex_jstep_lines1,
    ex_jstep_lines2,
    ex_jstruc_slice1,
    ex_outfileslice1,
    example_sp_outfile_path,
)

if TYPE_CHECKING:
    from collections.abc import Callable


@pytest.mark.parametrize(
    ("init_meth", "init_var", "add_checks"),
    [
        (lambda x: JDFTXOutfile.from_file(x), example_sp_outfile_path, lambda dir_repr: None),
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
    """Test the __dir__ method of the object.

    Args:
        init_meth: The method to initialize the object.
        init_var: The variable to initialize the object.
        add_checks: The checks to add to the test in case the output must pass additional checks.
    """
    dirrable = init_meth(init_var)
    dir_repr = dir(dirrable)
    add_checks(dir_repr)


@pytest.mark.parametrize(
    ("init_meth", "init_var", "add_checks"),
    [
        (lambda x: JDFTXOutfile.from_file(x), example_sp_outfile_path, lambda dir_repr: None),
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
    """Test the __repr__ method of the object.

    Args:
        init_meth: The method to initialize the object.
        init_var: The variable to initialize the object.
        add_checks: The checks to add to the test in case the output must pass additional checks.
    """
    reprrable = init_meth(init_var)
    repr_repr = repr(reprrable)
    add_checks(repr_repr)


@pytest.mark.parametrize(
    ("init_meth", "init_var", "add_checks"),
    [
        (lambda x: JDFTXOutfile.from_file(x), example_sp_outfile_path, lambda dir_repr: None),
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
def test_str_repr(init_meth: Callable, init_var: Any, add_checks: Callable) -> None:
    strrable = init_meth(init_var)
    str_repr = str(strrable)
    add_checks(str_repr)


@pytest.mark.parametrize(
    ("init_meth", "init_var", "add_checks"),
    [
        (lambda x: JDFTXOutfile.from_file(x), example_sp_outfile_path, lambda dir_repr: None),
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
def test_to_dict_repr(init_meth: Callable, init_var: Any, add_checks: Callable) -> None:
    """Test the to_dict method of the object.

    Args:
        init_meth: The method to initialize the object.
        init_var: The variable to initialize the object.
        add_checks: The checks to add to the test in case the output must pass additional checks.
    """
    dictable = init_meth(init_var)
    dict_repr = dictable.to_dict()
    add_checks(dict_repr)
