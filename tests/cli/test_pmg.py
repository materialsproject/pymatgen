from __future__ import annotations

import runpy
import sys
import textwrap
import types

import pytest

from pymatgen.cli import pmg
from tests.testing import TEST_FILES_DIR

VASP_IN_DIR = TEST_FILES_DIR / "io" / "vasp" / "inputs"


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


@pytest.mark.parametrize(
    "argv",
    [["--help"], ["analyze", "--help"], ["structure", "--help"], ["config", "--help"]],
)
def test_help(capsys, argv):
    with pytest.raises(SystemExit) as exc:
        pmg.main(argv)
    assert exc.value.code == 0
    assert "usage:" in capsys.readouterr().out


def test_unknown_subcommand(capsys):
    with pytest.raises(SystemExit) as exc:
        pmg.main(["not-a-command"])
    assert exc.value.code == 2
    assert "invalid choice: 'not-a-command'" in capsys.readouterr().err


def test_argcomplete_hook(monkeypatch):
    """argcomplete is optional; when importable, the parser is registered with it."""
    calls = []
    monkeypatch.setitem(sys.modules, "argcomplete", types.SimpleNamespace(autocomplete=calls.append))
    pmg.main(["diff", "--incar", str(VASP_IN_DIR / "INCAR"), str(VASP_IN_DIR / "INCAR_2")])
    assert len(calls) == 1
    assert calls[0].prog


def test_dunder_main(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["pmg", "diff", "-i", str(VASP_IN_DIR / "INCAR"), str(VASP_IN_DIR / "INCAR")])
    with pytest.raises(SystemExit) as exc:
        runpy.run_path(pmg.__file__, run_name="__main__")
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "SAME PARAMS" in out
    # Diffing a file against itself yields no different params.
    assert out.rstrip().endswith("----------------")


@pytest.mark.parametrize(
    ("extra_args", "excluded"),
    [([], []), (["-e", "Li,O"], ["Li", "O"])],
)
def test_pmg_view_fake_vis(monkeypatch, extra_args, excluded):
    """Exercise `pmg view` without VTK by substituting a fake StructureVis."""
    calls = {}

    class FakeStructureVis:
        def __init__(self, excluded_bonding_elements):
            calls["excluded"] = excluded_bonding_elements

        def set_structure(self, structure):
            calls["structure"] = structure

        def show(self):
            calls["shown"] = True

    fake_module = types.ModuleType("pymatgen.vis.structure_vtk")
    fake_module.StructureVis = FakeStructureVis
    monkeypatch.setitem(sys.modules, "pymatgen.vis.structure_vtk", fake_module)

    assert pmg.main(["view", str(VASP_IN_DIR / "POSCAR_Li2O"), *extra_args]) == 0
    assert calls["excluded"] == excluded
    assert calls["structure"].reduced_formula == "Li2O"
    assert calls["shown"]


def test_pmg_view(monkeypatch):
    pytest.importorskip("vtk", reason="vtk is not available")

    called = {}

    def fake_show(self):
        called["show_called"] = True

    monkeypatch.setattr("pymatgen.vis.structure_vtk.StructureVis.show", fake_show)

    filename = f"{VASP_IN_DIR}/POSCAR"

    pmg.main(["view", str(filename)])
    assert called["show_called"]
