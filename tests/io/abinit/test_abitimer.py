from __future__ import annotations

import math

import pytest

from pymatgen.io.abinit.abitimer import AbinitTimerParser
from pymatgen.util.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/io/abinit"


def abiref_file(filename):
    """Return absolute path to filename in ~pymatgen/tests/files/abinit."""
    return f"{TEST_DIR}/{filename}"


class TestAbinitTimer:
    def test_timer_parser(self):
        tparser = AbinitTimerParser()
        parsed = tparser.parse([abiref_file("mgb2_scf.abo")])
        assert len(parsed) == 1

        tparser = AbinitTimerParser()
        parsed = tparser.parse([abiref_file("si_scf_v10.2.7.abo")])
        assert len(parsed) == 1

    def test_timer(self):
        tparser = AbinitTimerParser()
        tparser.parse([abiref_file("mgb2_scf.abo")])

        assert len(tparser.timers()) == 1
        timer = tparser.timers()[0]
        assert math.isclose(timer.cpu_time, 5.1)
        assert timer.get_dataframe() is not None

        assert timer.to_table()

        assert len(timer.get_values("cpu_time")) == 30

        assert timer.sum_sections("cpu_time") == pytest.approx(3.298999999)

        assert timer.cpuwall_histogram(show=False)
        assert timer.pie(show=False)
        assert timer.scatter_hist(show=False)

    def test_timer_10_2_7(self):
        tparser = AbinitTimerParser()
        tparser.parse([abiref_file("si_scf_v10.2.7.abo")])

        assert len(tparser.timers()) == 1
        timer = tparser.timers()[0]
        assert math.isclose(timer.cpu_time, 0.3)
        assert timer.get_dataframe() is not None

        assert timer.to_table()

        assert len(timer.get_values("cpu_time")) == 30

        assert timer.sum_sections("cpu_time") == pytest.approx(0.307)

        assert timer.cpuwall_histogram(show=False)
        assert timer.pie(show=False)
        assert timer.scatter_hist(show=False)
