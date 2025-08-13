from __future__ import annotations

from typing import Any

import pytest
from pytest import approx

from pymatgen.io.jdftx.jelstep import JElStep, JElSteps

from .outputs_test_utils import (
    ex_fillings_line1,
    ex_fillings_line1_known,
    ex_iter_line1,
    ex_iter_line1_known,
    ex_subspace_line1,
    ex_subspace_line1_known,
)
from .outputs_test_utils import ex_jstep_known1 as ex_known1
from .outputs_test_utils import ex_jstep_known2 as ex_known2
from .outputs_test_utils import ex_jstep_lines1 as ex_lines1
from .outputs_test_utils import ex_jstep_lines2 as ex_lines2


def is_right_known(val: Any, ex_known_val: Any):
    if isinstance(val, type(ex_known_val)):
        result = ex_known_val == approx(val) if isinstance(val, float) else ex_known_val == val
    else:
        result = False
    return result


@pytest.mark.parametrize(
    (
        "exfill_line",
        "exfill_known",
        "exiter_line",
        "exiter_known",
        "exsubspace_line",
        "exsubspace_known",
        "etype",
        "eitertype",
    ),
    [
        (
            ex_fillings_line1,
            ex_fillings_line1_known,
            ex_iter_line1,
            ex_iter_line1_known,
            ex_subspace_line1,
            ex_subspace_line1_known,
            "F",
            "ElecMinimize",
        )
    ],
)
def test_JElStep_known(
    exfill_line: str,
    exfill_known: dict[str, float],
    exiter_line: str,
    exiter_known: dict[str, float],
    exsubspace_line: str,
    exsubspace_known: dict[str, float],
    etype: str,
    eitertype,
):
    ex_lines_collect = [exiter_line, exfill_line, exsubspace_line, ""]  # Empty line added for coverage
    jei = JElStep._from_lines_collect(ex_lines_collect, eitertype, etype)
    str(jei)
    ex_known = {}
    for dictlike in [exfill_known, exiter_known, exsubspace_known]:
        ex_known.update(dictlike)
    for var in [
        "mu",
        "nelectrons",
        "abs_magneticmoment",
        "tot_magneticmoment",
        "nstep",
        "e",
        "grad_k",
        "alpha",
        "linmin",
        "t_s",
        "subspacerotationadjust",
    ]:
        val = getattr(jei, var)
        assert is_right_known(val, ex_known[var])


@pytest.mark.parametrize(
    ("ex_lines", "ex_knowns", "etype", "eitertype"),
    [
        (
            [ex_lines1, ex_lines2],
            [ex_known1, ex_known2],
            "F",
            "ElecMinimize",
        )
    ],
)
def test_JElSteps_known(
    ex_lines: list[list[str]],
    ex_knowns: list[dict],
    etype: str,
    eitertype,
):
    text_slice = [line for exl in ex_lines for line in exl]
    jeis = JElSteps._from_text_slice(text_slice, opt_type=eitertype, etype=etype)
    for var in [
        "mu",
        "nelectrons",
        "abs_magneticmoment",
        "tot_magneticmoment",
        "nstep",
        "e",
        "grad_k",
        "alpha",
        "linmin",
        "t_s",
        "subspacerotationadjust",
    ]:
        val = getattr(jeis, var)
        assert is_right_known(val, ex_knowns[-1][var])
        for i in range(len(ex_lines)):
            val2 = getattr(jeis[i], var)
            assert is_right_known(val2, ex_knowns[i][var])
