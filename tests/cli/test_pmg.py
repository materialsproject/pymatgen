from __future__ import annotations

import textwrap

import pytest

from pymatgen.cli import pmg
from pymatgen.util.testing import VASP_IN_DIR


def test_main(capsys):
    with pytest.raises(SystemExit) as exc:
        pmg.main([])

    assert str(exc.value) == "Please specify a command."

    out, _err = capsys.readouterr()
    assert "`pmg` is a convenient script that uses `pymatgen`" in out


def test_pmg_diff(capsys):
    pmg.main(["diff", "--incar", f"{VASP_IN_DIR}/INCAR", f"{VASP_IN_DIR}/INCAR_2"])
    captured = capsys.readouterr()

    # The header (INCAR path) is machine dependent
    body = "SAME PARAMS" + captured.out.split("SAME PARAMS", 1)[1]

    expected = textwrap.dedent("""\
        SAME PARAMS
        ---------------
        EDIFF             0.0001                         0.0001
        IBRION            2                              2
        ISIF              3                              3
        ISPIN             2                              2
        LMAXMIX           4                              4
        LORBIT            11                             11
        LREAL             Auto                           Auto
        PREC              Accurate                       Accurate
        SIGMA             0.05                           0.05

        DIFFERENT PARAMS
        ----------------
        ALGO              Damped                         Fast
        ENCUT             500.0
        ENCUTFOCK         0.0
        HFSCREEN          0.207
        ICHARG                                           1
        ISMEAR            0                              -5
        ISPIND            2
        LCHARG            True
        LDAUPRINT                                        1
        LHFCALC           True
        LPLANE            True
        LSCALU            False
        LWAVE             True                           False
        MAGMOM            1*6.00 2*-6.00 1*6.00 20*0.60
        NELM                                             100
        NELMIN                                           3
        NKRED             2
        NPAR              8                              1
        NSIM              1
        NSW               99                             51
        NUPDOWN           0.0
        TIME              0.4
    """)

    def normalize(line: str) -> str:
        # collapse multiple spaces into one
        return " ".join(line.split())

    expected_lines = [normalize(line) for line in expected.strip().splitlines()]
    body_lines = [normalize(line) for line in body.strip().splitlines()]

    for exp in expected_lines:
        assert exp in body_lines


def test_pmg_query():
    # TODO: add test
    pass


def test_pmg_view(monkeypatch):
    pytest.importorskip("vtk", reason="vtk is not available")

    called = {}

    def fake_show(self):
        called["show_called"] = True

    monkeypatch.setattr("pymatgen.vis.structure_vtk.StructureVis.show", fake_show)

    filename = f"{VASP_IN_DIR}/POSCAR"

    pmg.main(["view", str(filename)])
    assert called["show_called"]
