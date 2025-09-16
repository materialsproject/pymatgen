from __future__ import annotations

import copy
import os
import pickle
import re
import warnings
from shutil import copyfile
from unittest.mock import patch

import numpy as np
import orjson
import pytest
import scipy.constants as const
from monty.io import zopen
from monty.serialization import loadfn
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.core import SETTINGS
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.vasp.inputs import (
    POTCAR_STATS_PATH,
    BadIncarWarning,
    BadPoscarWarning,
    Incar,
    Kpoints,
    KpointsSupportedModes,
    PmgVaspPspDirError,
    Poscar,
    Potcar,
    PotcarSingle,
    UnknownPotcarWarning,
    VaspInput,
    _gen_potcar_summary_stats,
)
from pymatgen.util.testing import FAKE_POTCAR_DIR, TEST_FILES_DIR, VASP_IN_DIR, VASP_OUT_DIR, MatSciTest

# Filter some expected warnings
warnings.filterwarnings(
    "ignore", message=r"POTCAR data with symbol .* is not known to pymatgen", category=UnknownPotcarWarning
)
warnings.filterwarnings("ignore", message=r"missing .* POTCAR directory", category=UserWarning)


# make sure _gen_potcar_summary_stats runs and works with all tests in this file
_SUMM_STATS = _gen_potcar_summary_stats(append=False, vasp_psp_dir=str(FAKE_POTCAR_DIR), summary_stats_filename=None)


@pytest.fixture(autouse=True)
def _mock_complete_potcar_summary_stats(monkeypatch: pytest.MonkeyPatch) -> None:
    # Override POTCAR library to use fake scrambled POTCARs
    monkeypatch.setitem(SETTINGS, "PMG_VASP_PSP_DIR", str(FAKE_POTCAR_DIR))
    monkeypatch.setattr(PotcarSingle, "_potcar_summary_stats", _SUMM_STATS)

    # The fake POTCAR library is pretty big even with just a few sub-libraries
    # just copying over entries to work with PotcarSingle.is_valid
    for func in PotcarSingle.functional_dir:
        if func in _SUMM_STATS:
            continue
        if "pbe" in func.lower() or "pw91" in func.lower():
            # Generate POTCAR hashes on the fly
            _SUMM_STATS[func] = _SUMM_STATS["PBE_54_W_HASH"].copy()
        elif "lda" in func.lower() or "perdew_zunger81" in func.lower():
            _SUMM_STATS[func] = _SUMM_STATS["LDA_64"].copy()


@pytest.mark.filterwarnings(
    "ignore:POTCAR data with symbol .* is not known to pymatgen:pymatgen.io.vasp.inputs.UnknownPotcarWarning"
)
class TestPoscar(MatSciTest):
    def test_init(self):
        comp = Structure.from_file(f"{VASP_IN_DIR}/POSCAR").composition
        assert comp == Composition("Fe4P4O16")

        # VASP 4 type with symbols at the end.
        poscar_str = """Test1
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
direct
0.000000 0.000000 0.000000 Si
0.750000 0.500000 0.750000 F
"""
        poscar = Poscar.from_str(poscar_str)
        assert poscar.structure.composition == Composition("SiF")

        poscar_str = ""
        with pytest.raises(ValueError, match="Empty POSCAR"):
            Poscar.from_str(poscar_str)

        # VASP 4 style file with default names, i.e. no element symbol found.
        poscar_str = """Test2
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
direct
0.000000 0.000000 0.000000
0.750000 0.500000 0.750000
"""
        with pytest.warns(BadPoscarWarning, match="Elements in POSCAR cannot be determined"):
            poscar = Poscar.from_str(poscar_str)
        assert poscar.structure.composition == Composition("HHe")

        # VASP 4 style file with default names, i.e. no element symbol found.
        poscar_str = """Test3
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
Selective dynamics
direct
0.000000 0.000000 0.000000 T T T Si
0.750000 0.500000 0.750000 F F F O
"""
        poscar = Poscar.from_str(poscar_str)
        selective_dynamics = [list(x) for x in poscar.selective_dynamics]

        assert selective_dynamics == [[True, True, True], [False, False, False]]
        self.selective_poscar = poscar

    def test_from_file(self):
        filepath = f"{VASP_IN_DIR}/POSCAR_symbols_natoms_multilines"
        poscar = Poscar.from_file(filepath, check_for_potcar=False, read_velocities=False)
        ordered_expected_elements = [
            "Fe",
            "Cr",
            "Fe",
            "Fe",
            "Cr",
            "Cr",
            "Cr",
            "Cr",
            "Fe",
            "Fe",
            "Cr",
            "Fe",
            "Cr",
            "Fe",
            "Fe",
            "Cr",
            "Fe",
            "Cr",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Cr",
            "Fe",
            "Ni",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Cr",
            "Cr",
            "Cr",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Cr",
            "Fe",
            "Fe",
            "Ni",
            "Fe",
            "Fe",
            "Fe",
            "Cr",
            "Cr",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
        ]
        assert [site.specie.symbol for site in poscar.structure] == ordered_expected_elements

    def test_from_file_potcar_overwrite_elements(self):
        # Make sure overwriting elements triggers warning
        # The following POTCAR has elements as ["Fe", "P", "O"]
        copyfile(f"{VASP_IN_DIR}/POSCAR_Fe3O4", tmp_poscar_path := f"{self.tmp_path}/POSCAR")
        copyfile(f"{VASP_IN_DIR}/fake_potcars/POTCAR.gz", f"{self.tmp_path}/POTCAR.gz")

        with pytest.warns(UserWarning, match="Elements in POSCAR would be overwritten"):
            _poscar = Poscar.from_file(tmp_poscar_path)

        # Make sure ValueError is raised when POTCAR has fewer elements
        copyfile(f"{VASP_IN_DIR}/POSCAR_Fe3O4", tmp_poscar_path := f"{self.tmp_path}/POSCAR")
        copyfile(f"{VASP_IN_DIR}/POTCAR_C2.gz", f"{self.tmp_path}/POTCAR.gz")

        with pytest.raises(ValueError, match="has fewer elements than POSCAR"):
            _poscar = Poscar.from_file(tmp_poscar_path)

        # Make sure no warning/error is triggered for compatible elements
        copyfile(f"{VASP_IN_DIR}/POSCAR_C2", tmp_poscar_path := f"{self.tmp_path}/POSCAR")
        copyfile(f"{VASP_IN_DIR}/POTCAR_C2.gz", f"{self.tmp_path}/POTCAR.gz")

        with warnings.catch_warnings(record=True) as record:
            _poscar = Poscar.from_file(tmp_poscar_path)
        assert not any("Elements in POSCAR would be overwritten" in str(warning.message) for warning in record)

    def test_from_str_default_names(self):
        """Similar to test_from_file_bad_potcar, ensure "default_names"
        (likely read from POTCAR) triggers correct warning/error.
        """
        poscar_str = """bad default names
1.1
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
Si F
1 1
cart
0.000000   0.00000000   0.00000000
3.840198   1.50000000   2.35163175
"""
        # ValueError if default_names shorter than POSCAR
        with pytest.raises(ValueError, match=re.escape("default_names=['Si'] (likely from POTCAR) has fewer elements")):
            _poscar = Poscar.from_str(poscar_str, default_names=["Si"])

        # Warning if overwrite triggered
        with pytest.warns(UserWarning, match="Elements in POSCAR would be overwritten"):
            poscar = Poscar.from_str(poscar_str, default_names=["Si", "O"])
        assert poscar.site_symbols == ["Si", "O"]

        # Assert no warning if using the same elements (or superset)
        with warnings.catch_warnings(record=True) as record:
            poscar = Poscar.from_str(poscar_str, default_names=["Si", "F"])
            assert poscar.site_symbols == ["Si", "F"]

            poscar = Poscar.from_str(poscar_str, default_names=["Si", "F", "O"])
            assert poscar.site_symbols == ["Si", "F"]
        assert not any("Elements in POSCAR would be overwritten" in str(warning.message) for warning in record)

        # Make sure it could be bypassed (by using None, when not check_for_potcar)
        with warnings.catch_warnings(record=True) as record:
            _poscar = Poscar.from_str(poscar_str, default_names=None)
        assert not any("Elements in POSCAR would be overwritten" in str(warning.message) for warning in record)

    def test_from_str_default_names_vasp4(self):
        """Poscar.from_str with default_names given could also be used to
        convert a VASP 4 formatted POSCAR to VASP 5/6, by inserting
        elements to the 5th (0-indexed) line.
        """
        poscar_str_vasp4 = """vasp 4
1.1
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
cart
0.000000   0.00000000   0.00000000
3.840198   1.50000000   2.35163175
"""
        # Test overwrite warning
        with pytest.warns(UserWarning, match="Elements in POSCAR would be overwritten"):
            poscar_vasp4 = Poscar.from_str(poscar_str_vasp4, default_names=["H", "He"])
        assert poscar_vasp4.site_symbols == ["H", "He"]

        # Test default names longer than need but should work the same
        poscar_vasp4 = Poscar.from_str(poscar_str_vasp4, default_names=["H", "He", "Li"])
        assert poscar_vasp4.site_symbols == ["H", "He"]

        with pytest.raises(ValueError, match=re.escape("default_names=['H'] (likely from POTCAR) has fewer elements")):
            _poscar_vasp4 = Poscar.from_str(poscar_str_vasp4, default_names=["H"])

    def test_as_from_dict(self):
        poscar_str = """Test3
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
Selective dynamics
direct
0.000000 0.000000 0.000000 T T T Si
0.750000 0.500000 0.750000 F F F O
"""
        poscar = Poscar.from_str(poscar_str)
        dct = poscar.as_dict()
        poscar2 = Poscar.from_dict(dct)
        assert poscar2.comment == "Test3"
        assert all(poscar2.selective_dynamics[0])
        assert not all(poscar2.selective_dynamics[1])

    def test_cart_scale(self):
        poscar_str = """Test1
1.1
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
Si F
1 1
cart
0.000000   0.00000000   0.00000000
3.840198   1.50000000   2.35163175
"""
        poscar = Poscar.from_str(poscar_str)
        site = poscar.structure[1]
        assert_allclose(site.coords, np.array([3.840198, 1.5, 2.35163175]) * 1.1)

    def test_significant_figures(self):
        si = 14
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]

        # Silicon structure for testing.
        lattice = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, [si, si], coords)
        poscar = Poscar(struct)
        expected_str = """Si2
1.0
   3.84    0.00    0.00
   1.92    3.33    0.00
   0.00   -2.22    3.14
Si
2
direct
   0.00    0.00    0.00 Si
   0.75    0.50    0.75 Si
"""

        actual_str = poscar.get_str(significant_figures=2)
        assert actual_str == expected_str, "Wrong POSCAR output!"

    def test_str(self):
        si = 14
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]

        # Silicon structure for testing.
        lattice = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, [si, si], coords)
        poscar = Poscar(struct)
        expected_str = """Si2
1.0
   3.8401979336999998    0.0000000000000000    0.0000000000000000
   1.9200989667999999    3.3257101909000002    0.0000000000000000
   0.0000000000000000   -2.2171384942999999    3.1355090603000000
Si
2
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Si
   0.7500000000000000    0.5000000000000000    0.7500000000000000 Si
"""

        assert str(poscar) == expected_str, "Wrong POSCAR output!"

        # VASP 4 type with symbols at the end.
        poscar_str = """Test1
1.0
-3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
direct
0.000000 0.000000 0.000000 Si
0.750000 0.500000 0.750000 F
"""

        expected = """Test1
1.0
   3.8401980000000000   -0.0000000000000000   -0.0000000000000000
  -1.9200990000000000   -3.3257099999999999   -0.0000000000000000
  -0.0000000000000000    2.2171379999999998   -3.1355089999999999
Si F
1 1
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Si
   0.7500000000000000    0.5000000000000000    0.7500000000000000 F
"""
        poscar = Poscar.from_str(poscar_str)
        assert str(poscar) == expected

    def test_from_md_run(self):
        # Parsing from an MD type run with velocities and predictor corrector data
        poscar = Poscar.from_file(f"{VASP_OUT_DIR}/CONTCAR.MD", check_for_potcar=False)
        assert np.sum(poscar.velocities) == approx(0.0065417961324)
        assert poscar.predictor_corrector[0][0][0] == approx(0.33387820e00)
        assert poscar.predictor_corrector[0][1][1] == approx(-0.10583589e-02)
        assert poscar.lattice_velocities is None

        # Parsing from an MD type run with velocities, predictor corrector data and lattice velocities
        poscar = Poscar.from_file(f"{VASP_OUT_DIR}/CONTCAR.MD.npt", check_for_potcar=False)
        assert np.sum(poscar.velocities) == approx(-0.06193299494)
        assert poscar.predictor_corrector[0][0][0] == approx(0.63981833)
        assert poscar.lattice_velocities.sum() == approx(16.49411358474)

    def test_write_md_poscar(self):
        # Parsing from an MD type run with velocities and predictor corrector data
        # And writing a new POSCAR from the new structure
        poscar = Poscar.from_file(f"{VASP_OUT_DIR}/CONTCAR.MD", check_for_potcar=False)

        path = f"{self.tmp_path}/POSCAR.testing.md"
        poscar.write_file(path)
        p3 = Poscar.from_file(path)

        assert_allclose(poscar.structure.lattice.abc, p3.structure.lattice.abc, 5)
        assert_allclose(poscar.velocities, p3.velocities, 5)
        assert_allclose(poscar.predictor_corrector, p3.predictor_corrector, 5)
        assert poscar.predictor_corrector_preamble == p3.predictor_corrector_preamble

        # Same as above except also has lattice velocities
        poscar = Poscar.from_file(f"{VASP_OUT_DIR}/CONTCAR.MD.npt", check_for_potcar=False)

        poscar.write_file(path)

        # check output produced for lattice velocities has required format and spaces
        # added in https://github.com/materialsproject/pymatgen/pull/3433
        with open(path, encoding="utf-8") as file:
            lines = file.readlines()
        pattern = (r"  [-| ]?\d\.\d{7}E[+-]\d{2}" * 3)[1:]
        for line in lines[18:24]:
            assert re.match(pattern, line.rstrip())

        p3 = Poscar.from_file(path)

        assert_allclose(poscar.structure.lattice.abc, p3.structure.lattice.abc, 5)
        assert_allclose(poscar.velocities, p3.velocities, 5)
        assert_allclose(poscar.predictor_corrector, p3.predictor_corrector, 5)
        assert poscar.predictor_corrector_preamble == p3.predictor_corrector_preamble
        assert_allclose(poscar.lattice_velocities, p3.lattice_velocities, 5)

    def test_setattr(self):
        filepath = f"{VASP_IN_DIR}/POSCAR"
        poscar = Poscar.from_file(filepath, check_for_potcar=False)
        with pytest.raises(ValueError, match="velocities array must be same length as the structure"):
            poscar.velocities = [[0, 0, 0]]
        poscar.selective_dynamics = np.array([[True, False, False]] * 24)
        expected = """Fe4P4O16
1.0
  10.4117668699494264    0.0000000000000000    0.0000000000000000
   0.0000000000000000    6.0671718799705294    0.0000000000000000
   0.0000000000000000    0.0000000000000000    4.7594895399768813
Fe P O
4 4 16
Selective dynamics
direct
   0.2187282200000000    0.7500000000000000    0.4748671100000000 T F F Fe
   0.2812717800000000    0.2500000000000000    0.9748671100000000 T F F Fe
   0.7187282200000000    0.7500000000000000    0.0251328900000000 T F F Fe
   0.7812717800000000    0.2500000000000000    0.5251328900000000 T F F Fe
   0.0946130900000000    0.2500000000000000    0.4182432700000000 T F F P
   0.4053869100000000    0.7500000000000000    0.9182432699999999 T F F P
   0.5946130900000000    0.2500000000000000    0.0817567300000000 T F F P
   0.9053869100000000    0.7500000000000000    0.5817567300000001 T F F P
   0.0433723100000000    0.7500000000000000    0.7071376700000001 T F F O
   0.0966424400000000    0.2500000000000000    0.7413203500000000 T F F O
   0.1657097400000000    0.0460723300000000    0.2853839400000000 T F F O
   0.1657097400000000    0.4539276700000000    0.2853839400000000 T F F O
   0.3342902600000000    0.5460723300000000    0.7853839400000000 T F F O
   0.3342902600000000    0.9539276700000000    0.7853839400000000 T F F O
   0.4033575600000000    0.7500000000000000    0.2413203500000000 T F F O
   0.4566276900000000    0.2500000000000000    0.2071376700000000 T F F O
   0.5433723100000000    0.7500000000000000    0.7928623299999999 T F F O
   0.5966424400000000    0.2500000000000000    0.7586796500000000 T F F O
   0.6657097400000000    0.0460723300000000    0.2146160600000000 T F F O
   0.6657097400000000    0.4539276700000000    0.2146160600000000 T F F O
   0.8342902600000000    0.5460723300000000    0.7146160600000000 T F F O
   0.8342902600000000    0.9539276700000000    0.7146160600000000 T F F O
   0.9033575600000000    0.7500000000000000    0.2586796500000000 T F F O
   0.9566276900000000    0.2500000000000000    0.2928623300000000 T F F O"""
        assert str(poscar).strip() == expected.strip()
        poscar.velocities = np.ones((24, 3))
        assert "velocities" in poscar.structure.site_properties

    def test_velocities(self):
        si = 14
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]

        # Silicon structure for testing.
        lattice = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, [si, si], coords)
        poscar = Poscar(struct)
        poscar.set_temperature(900)

        v = np.array(poscar.velocities)

        for x in np.sum(v, axis=0):
            assert x == approx(0, abs=1e-7)

        temperature = struct[0].specie.atomic_mass.to("kg") * np.sum(v**2) / (3 * const.k) * 1e10
        assert temperature == approx(900, abs=1e-4), "Temperature instantiated incorrectly"

        poscar.set_temperature(700)
        v = np.array(poscar.velocities)
        for x in np.sum(v, axis=0):
            assert x == approx(0, abs=1e-7), "Velocities initialized with a net momentum"

        temperature = struct[0].specie.atomic_mass.to("kg") * np.sum(v**2) / (3 * const.k) * 1e10
        assert temperature == approx(700, abs=1e-4), "Temperature instantiated incorrectly"

    def test_write(self):
        filepath = f"{VASP_IN_DIR}/POSCAR"
        poscar = Poscar.from_file(filepath)
        tmp_file = f"{self.tmp_path}/POSCAR.testing"
        poscar.write_file(tmp_file)
        poscar = Poscar.from_file(tmp_file)
        assert_allclose(poscar.structure.lattice.abc, poscar.structure.lattice.abc, 5)

    @pytest.mark.filterwarnings("ignore:Elements in POSCAR would be overwritten")
    def test_selective_dynamics(self):
        # Previously, this test relied on the existence of a file named POTCAR
        # that was sorted to the top of a list of POTCARs for the test to work.
        # That's far too brittle - isolating requisite files here
        copyfile(f"{VASP_IN_DIR}/POSCAR_Fe3O4", tmp_poscar_path := f"{self.tmp_path}/POSCAR")
        copyfile(f"{VASP_IN_DIR}/fake_potcars/POTCAR.gz", f"{self.tmp_path}/POTCAR.gz")

        poscar = Poscar.from_file(tmp_poscar_path)
        structure = poscar.structure

        # Fix bottom half
        fixed_indices = structure.frac_coords[:, 2] >= 0.5

        poscar = Poscar(structure, selective_dynamics=np.tile(fixed_indices.reshape(-1, 1), [1, 3]))
        selective_dynamics = [list(x) for x in poscar.selective_dynamics]

        assert selective_dynamics == [
            [True, True, True],
            [False, False, False],
            [False, False, False],
            [True, True, True],
            [False, False, False],
            [True, True, True],
            [True, True, True],
            [False, False, False],
            [False, False, False],
            [True, True, True],
            [True, True, True],
            [False, False, False],
            [True, True, True],
            [False, False, False],
        ]

    def test_invalid_selective_dynamics(self):
        """
        Check invalid selective dynamics info. The POSCAR string
        'invalid_poscar_str' represents a case with incorrect
        placement of selective dynamics information (Comment like 'Si' should
        be followed by selective dynamics values 'T' or 'F').
        """
        invalid_poscar_str = """POSCAR with invalid selective dynamics info
1.1
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
Si F
1 1
Selective dynamics
Cartesian
0.000000   0.00000000   0.00000000 Si T T F
3.840198   1.50000000   2.35163175 F T T F
"""
        with (
            pytest.warns(
                BadPoscarWarning,
                match="Selective dynamics values must be either 'T' or 'F'.",
            ),
            pytest.warns(BadPoscarWarning, match="Selective dynamics toggled with Fluorine"),
        ):
            Poscar.from_str(invalid_poscar_str)

    def test_selective_dynamics_with_fluorine(self):
        """
        Check ambiguous selective dynamics info when Fluorine(F) is
        included and position lines include comments.
        """
        poscar_str_with_fluorine = """Selective dynamics toggled with Fluorine
1.1
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
Si F
1 1
Selective dynamics
Cartesian
0.000000   0.00000000   0.00000000 Si T T F
3.840198   1.50000000   2.35163175 F T T F
"""
        with (
            pytest.warns(
                BadPoscarWarning,
                match=(
                    "Selective dynamics toggled with Fluorine element detected. "
                    "Make sure the 4th-6th entry each position line is selective dynamics info."
                ),
            ),
            pytest.warns(BadPoscarWarning, match="Selective dynamics values must be either"),
        ):
            Poscar.from_str(poscar_str_with_fluorine)

    def test_all_DOFs_relaxed(self):
        """
        A warning should be issued when selective dynamics is toggled
        while ALL degrees of freedom are relaxed.
        """
        poscar_str_all_dof_relaxed = """All degrees of freedom relaxed
1.1
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
Si O
1 1
Selective dynamics
Cartesian
0.000000   0.00000000   0.00000000 T T T
3.840198   1.50000000   2.35163175 T T T
"""
        with pytest.warns(
            BadPoscarWarning,
            match="Ignoring selective dynamics tag, as no ionic degrees of freedom were fixed.",
        ):
            Poscar.from_str(poscar_str_all_dof_relaxed)

    def test_vasp_6_4_2_format(self):
        # As of vasp 6.4.2, when using POTCARs with SHAs, there can
        # be a slash in the element names
        # Test that Poscar works for these too
        poscar_str = ""
        with open(f"{VASP_IN_DIR}/POSCAR_LiFePO4", encoding="utf-8") as file:
            for idx, line in enumerate(file):
                if idx == 5:
                    line = " ".join(f"{x}/" for x in line.split()) + "\n"
                poscar_str += line
        poscar = Poscar.from_str(poscar_str)
        assert poscar.structure.formula == "Li4 Fe4 P4 O16"


class TestIncar(MatSciTest):
    def setup_method(self):
        self.incar = Incar.from_file(f"{VASP_IN_DIR}/INCAR")

    def test_init(self):
        incar = self.incar
        incar["LDAU"] = "T"
        assert incar["ALGO"] == "Damped", "Wrong Algo"
        assert float(incar["EDIFF"]) == approx(1e-4), "Wrong EDIFF"
        assert isinstance(incar["LORBIT"], int)

    def test_lattice_constaints(self):
        incar_str = """
        ALGO = Fast
EDIFF = 0.00045000000000000004
ENCUT = 520.0
IBRION = 2
ISIF = 3
ISMEAR = -5
ISPIN = 2
LASPH = True
LATTICE_CONSTRAINTS = False False True
LORBIT = 11
LREAL = Auto
LWAVE = False
MAGMOM = 9*0.6
NELM = 100
NSW = 99
PREC = Accurate
SIGMA = 0.05"""
        incar = Incar.from_str(incar_str)
        assert incar["LATTICE_CONSTRAINTS"] == [False, False, True]

    def test_check_for_duplicate(self):
        incar_str: str = """encut = 400
        ENCUT = 500
        """
        with pytest.warns(
            BadIncarWarning,
            match=re.escape("Duplicate keys found (case-insensitive): ['ENCUT']"),
        ):
            Incar.from_str(incar_str)

        incar_dict = {"ALGO": "Fast", "algo": "fast"}
        with pytest.warns(
            BadIncarWarning,
            match=re.escape("Duplicate keys found (case-insensitive): ['ALGO']"),
        ):
            Incar.from_dict(incar_dict)

    def test_key_case_insensitive(self):
        """Verify that keys are case-insensitive by internally converting
        all keys to upper case. This includes operations such as:
        - set/get: Keys can be set and retrieved with any case.
        - update: Keys in updates are case-insensitive.
        - setdefault: Defaults are set and retrieved case-insensitively.
        """
        test_tag: str = "ENCUT"

        incar_str: str = f"""ALGO = Fast
        {test_tag} = 480
        EDIFF = 1e-07
        """

        # Test setter and getter
        incar: Incar = Incar.from_str(incar_str)
        incar[test_tag.lower()] = 490
        assert incar[test_tag.lower()] == 490
        assert incar[test_tag.upper()] == 490
        assert incar.get(test_tag.lower()) == 490
        assert incar.get(test_tag.upper()) == 490

        incar[test_tag.upper()] = 500
        assert incar[test_tag.lower()] == 500

        # Test delete
        del incar["algo"]
        assert "ALGO" not in incar

        # Test membership check
        assert test_tag.upper() in incar
        assert test_tag.lower() in incar

        # Test update
        incar.update({test_tag.lower(): 510})
        assert incar[test_tag] == 510

        incar.update({test_tag.upper(): 520})
        assert incar[test_tag] == 520

        # Test setdefault
        incar.setdefault("ismear", 0)
        assert incar["ISMEAR"] == 0

        incar.setdefault("NPAR", 4)
        assert incar["npar"] == 4

        # Test pop
        assert incar.pop(test_tag.lower()) == 520
        assert test_tag.upper() not in incar
        assert test_tag.lower() not in incar

        # Test pop with default value
        assert incar.pop("missing_key", 42) == 42

    def test_copy(self):
        incar2 = self.incar.copy()
        assert isinstance(incar2, Incar), f"Expected Incar, got {type(incar2).__name__}"
        assert incar2 == self.incar
        # modify incar2 and check that incar1 is not modified
        incar2["LDAU"] = "F"
        assert incar2["LDAU"] is False
        assert self.incar.get("LDAU") is None

    def test_diff(self):
        incar1 = Incar.from_file(f"{VASP_IN_DIR}/INCAR")
        incar2 = Incar.from_file(f"{VASP_IN_DIR}/INCAR_2")

        assert incar1.diff(incar2) == {
            "Different": {
                "NELM": {"INCAR1": None, "INCAR2": 100},
                "ISPIND": {"INCAR1": 2, "INCAR2": None},
                "LWAVE": {"INCAR1": True, "INCAR2": False},
                "LDAUPRINT": {"INCAR1": None, "INCAR2": 1},
                "MAGMOM": {
                    "INCAR1": [
                        6,
                        -6,
                        -6,
                        6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                    ],
                    "INCAR2": None,
                },
                "NELMIN": {"INCAR1": None, "INCAR2": 3},
                "ENCUTFOCK": {"INCAR1": 0.0, "INCAR2": None},
                "HFSCREEN": {"INCAR1": 0.207, "INCAR2": None},
                "LSCALU": {"INCAR1": False, "INCAR2": None},
                "ENCUT": {"INCAR1": 500, "INCAR2": None},
                "NSIM": {"INCAR1": 1, "INCAR2": None},
                "ICHARG": {"INCAR1": None, "INCAR2": 1},
                "NSW": {"INCAR1": 99, "INCAR2": 51},
                "NKRED": {"INCAR1": 2, "INCAR2": None},
                "NUPDOWN": {"INCAR1": 0, "INCAR2": None},
                "LCHARG": {"INCAR1": True, "INCAR2": None},
                "LPLANE": {"INCAR1": True, "INCAR2": None},
                "ISMEAR": {"INCAR1": 0, "INCAR2": -5},
                "NPAR": {"INCAR1": 8, "INCAR2": 1},
                "SYSTEM": {
                    "INCAR1": "id=[0] dblock_code=[97763-ICSD] formula=[Li Mn (P O4)] sg_name=[P n m a]",
                    "INCAR2": (
                        "id=[91090] dblock_code=[20070929235612LiNiO-59.53134651-VASP] "
                        "formula=[Li3 Ni3 O6] sg_name=[R-3m]"
                    ),
                },
                "ALGO": {"INCAR1": "Damped", "INCAR2": "Fast"},
                "LHFCALC": {"INCAR1": True, "INCAR2": None},
                "TIME": {"INCAR1": 0.4, "INCAR2": None},
            },
            "Same": {
                "IBRION": 2,
                "PREC": "Accurate",
                "ISIF": 3,
                "LMAXMIX": 4,
                "LREAL": "Auto",
                "ISPIN": 2,
                "EDIFF": 0.0001,
                "LORBIT": 11,
                "SIGMA": 0.05,
            },
        }

    def test_as_dict_and_from_dict(self):
        dct = self.incar.as_dict()
        incar2 = Incar.from_dict(dct)
        assert self.incar == incar2
        dct["MAGMOM"] = [Magmom([1, 2, 3]).as_dict()]
        incar3 = Incar.from_dict(dct)
        assert incar3["MAGMOM"] == [Magmom([1, 2, 3])]

    def test_from_file_and_from_dict(self):
        """
        Init from file (from str) should yield the same results as from dict.

        Previously init Incar from dict would bypass the proc_val method for
        float/int, and might yield values in wrong type.
        """
        # Init from dict
        incar_dict = {"ENCUT": 500, "GGA": "PS", "NELM": 60.0, "SYSTEM": "This should NOT BE capitalized"}
        incar_from_dict = Incar(incar_dict)

        # Init from file (from string)
        incar_str = """\
        ENCUT = 500
        GGA = PS
        NELM = 60.0
        SYSTEM = This should NOT BE capitalized
        """

        with open("INCAR", "w", encoding="utf-8") as f:
            f.write(incar_str)

        incar_from_file = Incar.from_file("INCAR")

        # Make sure int/float is cast to correct type when init from dict
        assert incar_from_dict["GGA"] == "Ps"
        assert isinstance(incar_from_dict["ENCUT"], float)
        assert isinstance(incar_from_dict["NELM"], int)

        assert incar_from_dict == incar_from_file
        for key in incar_from_dict:
            assert type(incar_from_dict[key]) is type(incar_from_file[key])
        assert incar_from_file["SYSTEM"] == "This should NOT BE capitalized"

    def test_write(self):
        tmp_file = f"{self.tmp_path}/INCAR.testing"
        self.incar.write_file(tmp_file)
        incar = Incar.from_file(tmp_file)
        assert incar == self.incar

    def test_get_str(self):
        incar_str = self.incar.get_str(pretty=True, sort_keys=True)
        expected = """ALGO       =  Damped
EDIFF      =  0.0001
ENCUT      =  500.0
ENCUTFOCK  =  0.0
HFSCREEN   =  0.207
IBRION     =  2
ISIF       =  3
ISMEAR     =  0
ISPIN      =  2
ISPIND     =  2
LCHARG     =  True
LHFCALC    =  True
LMAXMIX    =  4
LORBIT     =  11
LPLANE     =  True
LREAL      =  Auto
LSCALU     =  False
LWAVE      =  True
MAGMOM     =  1*6.0 2*-6.0 1*6.0 20*0.6
NKRED      =  2
NPAR       =  8
NSIM       =  1
NSW        =  99
NUPDOWN    =  0.0
PREC       =  Accurate
SIGMA      =  0.05
SYSTEM     =  id=[0] dblock_code=[97763-ICSD] formula=[Li Mn (P O4)] sg_name=[P n m a]
TIME       =  0.4"""
        assert incar_str == expected

    def test_lsorbit_magmom(self):
        magmom1 = [[0.0, 0.0, 3.0], [0, 1, 0], [2, 1, 2]]
        magmom2 = [-1, -1, -1, 0, 0, 0, 0, 0]
        magmom4 = [Magmom([1.0, 2.0, 2.0])]

        ans_string1 = "LANGEVIN_GAMMA = 10 10 10\nLSORBIT = True\nMAGMOM = 0.0 0.0 3.0 0 1 0 2 1 2\n"
        ans_string2 = "LANGEVIN_GAMMA = 10\nLSORBIT = True\nMAGMOM = 3*3*-1 3*5*0\n"
        ans_string3 = "LSORBIT = False\nMAGMOM = 2*-1 2*9\n"
        ans_string4_nolsorbit = "LANGEVIN_GAMMA = 10\nLSORBIT = False\nMAGMOM = 1*3.0\n"
        ans_string4_lsorbit = "LANGEVIN_GAMMA = 10\nLSORBIT = True\nMAGMOM = 1.0 2.0 2.0\n"

        incar = Incar({})
        incar["MAGMOM"] = magmom1
        incar["LSORBIT"] = "T"
        incar["LANGEVIN_GAMMA"] = [10, 10, 10]
        assert ans_string1 == str(incar)

        incar["MAGMOM"] = magmom2
        incar["LSORBIT"] = "T"
        incar["LANGEVIN_GAMMA"] = 10
        assert ans_string2 == str(incar)

        incar["MAGMOM"] = magmom4
        incar["LSORBIT"] = "F"
        assert ans_string4_nolsorbit == str(incar)
        incar["LSORBIT"] = "T"
        assert ans_string4_lsorbit == str(incar)

        incar = Incar.from_str(ans_string1)
        assert_allclose(incar["MAGMOM"], [[0.0, 0.0, 3.0], [0, 1, 0], [2, 1, 2]])
        assert incar["LANGEVIN_GAMMA"] == [10, 10, 10]

        incar = Incar.from_str(ans_string2)
        assert incar["MAGMOM"] == [
            [-1, -1, -1],
            [-1, -1, -1],
            [-1, -1, -1],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
        ]
        assert incar["LANGEVIN_GAMMA"] == [10]

        incar = Incar.from_str(ans_string3)
        assert not incar["LSORBIT"]
        assert incar["MAGMOM"] == [-1, -1, 9, 9]

    def test_quad_efg(self):
        incar1 = Incar({})
        incar1["LEFG"] = True
        incar1["QUAD_EFG"] = [0.0, 146.6, -25.58]
        ans_string1 = "LEFG = True\nQUAD_EFG = 0.0 146.6 -25.58\n"
        assert ans_string1 == str(incar1)
        incar2 = Incar.from_str(ans_string1)
        assert ans_string1 == str(incar2)

    def test_types(self):
        incar_str = """ALGO = Fast
ENCUT = 510
EDIFF = 1e-07
EINT = -0.85 0.85
IBRION = -1
ICHARG = 11
ISIF = 3
ISMEAR = 1
ISPIN = 1
LPARD = True
NBMOD = -3
PREC = Accurate
SIGMA = 0.1"""
        incar = Incar.from_str(incar_str)
        assert isinstance(incar["EINT"], list)
        assert incar["EINT"][0] == approx(-0.85)

        incar_str += "\nLHFCALC = .TRUE. ; HFSCREEN = 0.2"
        incar_str += "\nALGO = All;"
        incar = Incar.from_str(incar_str)
        assert incar["LHFCALC"]
        assert incar["HFSCREEN"] == approx(0.2)
        assert incar["ALGO"] == "All"

    def test_proc_types(self):
        assert Incar.proc_val("HELLO", "-0.85 0.85") == "-0.85 0.85"
        assert Incar.proc_val("ML_MODE", "train") == "train"
        assert Incar.proc_val("ML_MODE", "RUN") == "run"
        assert Incar.proc_val("ALGO", "fast") == "Fast"

    def test_check_params(self):
        # Triggers warnings when running into invalid parameters
        incar = Incar(
            {
                "ADDGRID": True,
                "ALGO": "Normal",
                "AMIN": 0.01,
                "ICHARG": 1,
                "MAGMOM": [1, 2, 4, 5],
                "ML_MODE": "RUN",  # lower case string
                "SYSTEM": "Hello world",  # as is string
                "ENCUT": 500,  # make sure float key is casted
                "GGA": "PS",  # test string case insensitivity
                "LREAL": True,  # special case: Union type
                "NBAND": 250,  # typo in tag
                "METAGGA": "SCAM",  # typo in value
                "EDIFF": 5 + 1j,  # value should be float
                "ISIF": 9,  # value not unknown
                "LASPH": 5,  # value should be bool
                "PHON_TLIST": "is_a_str",  # value should be a list
            }
        )
        with pytest.warns(BadIncarWarning) as record:
            incar.check_params()

        assert record[0].message.args[0] == "Cannot find NBAND in the list of INCAR tags"
        assert record[1].message.args[0] == "METAGGA: Cannot find Scam in the list of values"
        assert record[2].message.args[0] == "EDIFF: (5+1j) is not a float"
        assert record[3].message.args[0] == "ISIF: Cannot find 9 in the list of values"
        assert record[4].message.args[0] == "LASPH: 5 is not a bool"
        assert record[5].message.args[0] == "PHON_TLIST: Is_a_str is not a list"


class TestKpointsSupportedModes:
    @pytest.mark.parametrize("mode", ["Automatic", "Gamma", "Monkhorst", "Line_mode", "Cartesian", "Reciprocal"])
    def test_from_str(self, mode):
        expected = getattr(KpointsSupportedModes, mode)
        assert KpointsSupportedModes.from_str(mode) == expected
        assert KpointsSupportedModes.from_str(mode.lower()) == expected  # case insensitive
        assert KpointsSupportedModes.from_str(mode[0]) == expected  # only first letter matters

    def test_invalid_mode(self):
        mode = "InvalidMode"
        with pytest.raises(ValueError, match=f"Invalid Kpoint {mode=}"):
            KpointsSupportedModes.from_str(mode)


class TestKpoints:
    @pytest.mark.filterwarnings("ignore:Please use INCAR KSPACING tag")
    def test_init(self):
        # Automatic KPOINT grid
        filepath = f"{VASP_IN_DIR}/KPOINTS_auto"
        with pytest.warns(DeprecationWarning, match="Please use INCAR KSPACING tag"):
            kpoints = Kpoints.from_file(filepath)
        assert kpoints.comment == "Automatic mesh"
        assert kpoints.kpts == [(10,)], "Wrong kpoint lattice read"

        filepath = f"{VASP_IN_DIR}/KPOINTS_cartesian"
        kpoints = Kpoints.from_file(filepath)
        (
            assert_allclose(
                kpoints.kpts,
                [
                    (0.25, 0, 0),
                    (0, 0.25, 0),
                    (0, 0, 0.25),
                ],
            ),
            "Wrong kpoint lattice read",
        )
        assert_allclose(kpoints.kpts_shift, (0.5, 0.5, 0.5))

        # Gamma-centered Kpoint grid
        filepath = f"{VASP_IN_DIR}/KPOINTS_gamma"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.comment == "Gamma centered mesh"
        assert kpoints.kpts == [(4, 4, 4)]

        # Monkhorst-Pack Kpoint grid
        filepath = f"{VASP_IN_DIR}/KPOINTS"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.comment == "Auto-generated kpoints file: VaspIO.writeKPOINTS(XX,XX,500)"
        self.kpoints = kpoints
        assert kpoints.kpts == [(2, 4, 6)]

        # Line mode
        filepath = f"{VASP_IN_DIR}/KPOINTS_band"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.labels is not None
        assert kpoints.style == Kpoints.supported_modes.Line_mode
        assert str(kpoints).split("\n")[3] == "Reciprocal"

        # Explicit K-point mode
        filepath = f"{VASP_IN_DIR}/KPOINTS_explicit"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.kpts_weights is not None
        expected_kpt_str = """Example file
4
Cartesian
0.0 0.0 0.0 1 None
0.0 0.0 0.5 1 None
0.0 0.5 0.5 2 None
0.5 0.5 0.5 4 None"""
        assert str(kpoints).strip() == expected_kpt_str

        # Explicit tetrahedra method
        filepath = f"{VASP_IN_DIR}/KPOINTS_explicit_tet"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.comment == "Example file"
        assert kpoints.tet_connections == [(6, [1, 2, 3, 4])]

    def test_property_kpts(self):
        kpoints_0 = Kpoints(kpts=[[1, 1, 1]])
        assert kpoints_0.kpts == [(1, 1, 1)]

        kpoints_1 = Kpoints(kpts=[(1, 1, 1)])
        assert kpoints_1.kpts == [(1, 1, 1)]

        kpoints_2 = Kpoints(kpts=[np.array((1, 1, 1))])
        assert kpoints_2.kpts == [(1, 1, 1)]

        kpoints_3 = Kpoints(
            style=Kpoints.supported_modes.Line_mode,
            kpts=[[1, 1, 1], (2, 2, 2), np.array([3, 3, 3])],
        )
        assert kpoints_3.kpts == [(1, 1, 1), (2, 2, 2), (3, 3, 3)]

        kpoints_4 = Kpoints(kpts=[[1]])
        assert kpoints_4.kpts == [(1,)]

        kpoints_5 = Kpoints(kpts=[1, 1, 1])
        assert kpoints_5.kpts == [(1, 1, 1)]

    @pytest.mark.parametrize(
        "invalid_kpts",
        [
            (("1", "1", "1")),  # invalid data type
            ((1, 1)),  # length not 1 or 3
        ],
    )
    def test_property_kpts_invalid(self, invalid_kpts):
        with pytest.raises(ValueError, match="Invalid Kpoint"):
            Kpoints(kpts=invalid_kpts)

    def test_property_style(self):
        filepath = f"{VASP_IN_DIR}/KPOINTS"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.style == Kpoints.supported_modes.Monkhorst
        kpoints.style = "G"
        assert kpoints.style == Kpoints.supported_modes.Gamma

    def test_static_constructors(self):
        kpoints = Kpoints.gamma_automatic((3, 3, 3), [0, 0, 0])
        assert kpoints.style == Kpoints.supported_modes.Gamma
        assert kpoints.kpts == [(3, 3, 3)]
        assert all(isinstance(kpt, int) for kpt in kpoints.kpts[0])

        kpoints = Kpoints.monkhorst_automatic((2, 2, 2), [0, 0, 0])
        assert kpoints.style == Kpoints.supported_modes.Monkhorst
        assert kpoints.kpts == [(2, 2, 2)]
        assert all(isinstance(kpt, int) for kpt in kpoints.kpts[0])

        with pytest.deprecated_call(match="Please use INCAR KSPACING tag"):
            kpoints = Kpoints.automatic(100)
        assert kpoints.style == Kpoints.supported_modes.Automatic
        assert kpoints.kpts == [(100,)]
        assert all(isinstance(kpt, int) for kpt in kpoints.kpts[0])

        filepath = f"{VASP_IN_DIR}/POSCAR"
        struct = Structure.from_file(filepath)
        kpoints = Kpoints.automatic_density(struct, 500)
        assert kpoints.kpts == [(1, 3, 3)]
        assert all(isinstance(kpt, int) for kpt in kpoints.kpts[0])
        assert kpoints.style == Kpoints.supported_modes.Gamma

        kpoints = Kpoints.automatic_density(struct, 500, force_gamma=True)
        assert kpoints.style == Kpoints.supported_modes.Gamma
        assert all(isinstance(kpt, int) for kpt in kpoints.kpts[0])

        kpoints = Kpoints.automatic_density_by_vol(struct, 1000)
        assert kpoints.kpts == [(6, 10, 13)]
        assert all(isinstance(kpt, int) for kpt in kpoints.kpts[0])
        assert kpoints.style == Kpoints.supported_modes.Gamma

        kpoints = Kpoints.automatic_density_by_lengths(struct, [50, 50, 1], force_gamma=True)
        assert kpoints.kpts == [(5, 9, 1)]
        assert all(isinstance(kpt, int) for kpt in kpoints.kpts[0]), kpoints.kpts
        assert kpoints.style == Kpoints.supported_modes.Gamma

        struct.make_supercell(3)
        kpoints = Kpoints.automatic_density(struct, 500)
        assert kpoints.kpts == [(1, 1, 1)]
        assert all(isinstance(kpt, int) for kpt in kpoints.kpts[0])
        assert kpoints.style == Kpoints.supported_modes.Gamma
        kpoints = Kpoints.from_str(
            """k-point mesh
            0
            G
            10 10 10
            0.5 0.5 0.5
            """
        )
        assert_allclose(kpoints.kpts_shift, [0.5, 0.5, 0.5])

    def test_as_dict_from_dict(self):
        kpts = Kpoints.monkhorst_automatic((2, 2, 2), [0, 0, 0])
        dct = kpts.as_dict()
        kpts_from_dict = Kpoints.from_dict(dct)
        assert kpts.kpts == kpts_from_dict.kpts
        assert kpts.style == kpts_from_dict.style
        assert kpts.kpts_shift == kpts_from_dict.kpts_shift

    def test_kpt_bands_as_dict_from_dict(self):
        file_name = f"{VASP_IN_DIR}/KPOINTS_band"
        kpts = Kpoints.from_file(file_name)
        dct = kpts.as_dict()

        assert orjson.dumps(dct).decode()
        # This doesn't work
        k2 = Kpoints.from_dict(dct)
        assert kpts.kpts == k2.kpts
        assert kpts.style == k2.style
        assert kpts.kpts_shift == k2.kpts_shift
        assert kpts.num_kpts == k2.num_kpts

    def test_pickle(self):
        kpts = Kpoints.gamma_automatic()
        pickle.dumps(kpts)

    def test_eq(self):
        auto_g_kpts = Kpoints.gamma_automatic()
        assert auto_g_kpts == auto_g_kpts  # noqa: PLR0124
        assert auto_g_kpts == Kpoints.gamma_automatic()
        file_kpts = Kpoints.from_file(f"{VASP_IN_DIR}/KPOINTS")
        assert file_kpts == Kpoints.from_file(f"{VASP_IN_DIR}/KPOINTS")
        assert auto_g_kpts != file_kpts
        auto_m_kpts = Kpoints.monkhorst_automatic((2, 2, 2), [0, 0, 0])
        assert auto_m_kpts == Kpoints.monkhorst_automatic((2, 2, 2), [0, 0, 0])
        assert auto_g_kpts != auto_m_kpts

    def test_copy(self):
        kpts = Kpoints.gamma_automatic()
        kpt_copy = kpts.copy()
        assert kpts == kpt_copy
        kpt_copy.style = Kpoints.supported_modes.Monkhorst
        assert kpts != kpt_copy

    def test_automatic_kpoint(self):
        poscar = Poscar.from_str(
            """Al1
1.0
2.473329 0.000000 1.427977
0.824443 2.331877 1.427977
0.000000 0.000000 2.855955
Al
1
direct
0.000000 0.000000 0.000000 Al"""
        )
        kpoints = Kpoints.automatic_density(poscar.structure, 1000)
        assert_allclose(kpoints.kpts[0], [10, 10, 10])

    def test_automatic_density_by_lengths(self):
        # Load a structure from a POSCAR file
        filepath = f"{VASP_IN_DIR}/POSCAR"
        structure = Structure.from_file(filepath)

        # test different combos of length densities and expected kpoints
        # TODO should test Monkhorst style case and force_gamma=True case
        for length_densities, expected_kpts, expected_style in [
            ([50, 50, 1], [(5, 9, 1)], Kpoints.supported_modes.Gamma),
            ([25, 50, 3], [(3, 9, 1)], Kpoints.supported_modes.Gamma),
            ([24, 48, 2], [(3, 8, 1)], Kpoints.supported_modes.Gamma),
        ]:
            kpoints = Kpoints.automatic_density_by_lengths(structure, length_densities)

            assert kpoints.kpts == expected_kpts

            assert kpoints.style == expected_style

        with pytest.raises(ValueError, match="The dimensions of length_densities must be 3, not 2"):
            Kpoints.automatic_density_by_lengths(structure, [50, 50])

    def test_automatic_monkhorst_vs_gamma_style_selection(self):
        structs = {key: Structure.from_file(f"{VASP_IN_DIR}/POSCAR_{key}") for key in ("bcc", "fcc", "hcp")}

        # bcc structures should allow both Monkhorst and Gamma
        for struct_type, struct in structs.items():
            for density in (500, 600, 700):
                kpoints = Kpoints.automatic_density(struct, density)
                if struct_type == "bcc" and density in (500, 600):
                    assert kpoints.style == Kpoints.supported_modes.Monkhorst
                else:
                    assert kpoints.style == Kpoints.supported_modes.Gamma

        # Kpoints.automatic_density_by_lengths
        for struct_type, struct in structs.items():
            for lengths in [50, 50, 50], [53, 53, 53], [56, 56, 56]:
                kpoints = Kpoints.automatic_density_by_lengths(struct, lengths)
                if struct_type == "bcc" and all(length % 2 == 0 for length in lengths):
                    assert kpoints.style == Kpoints.supported_modes.Monkhorst
                else:
                    assert kpoints.style == Kpoints.supported_modes.Gamma

        # Overkill test to make sure these methods always set the style to Gamma
        for len_density in range(1, 50):
            for struct_type, struct in structs.items():
                if struct_type != "bcc":
                    kpoints = Kpoints.automatic_density_by_lengths(struct, [len_density] * 3)
                    assert kpoints.style == Kpoints.supported_modes.Gamma

    def test_non_ascii_comment(self):
        """Non-ASCII comment like 'Γ' might not be encoded correctly in Windows."""
        kpoints = Kpoints.from_file(f"{VASP_IN_DIR}/KPOINTS_band")
        assert kpoints.labels[0] == "Γ", f"Γ is not encoded correctly, got {kpoints.labels[0]}"


@pytest.mark.filterwarnings(
    "ignore:POTCAR data with symbol .* is not known to pymatgen:pymatgen.io.vasp.inputs.UnknownPotcarWarning"
)
class TestPotcarSingle:
    def setup_method(self):
        self.psingle_Mn_pv = PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE/POTCAR.Mn_pv.gz")
        self.psingle_Fe = PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE/POTCAR.Fe.gz")
        self.psingle_Fe_54 = PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe.gz")

        self.Mn_pv_attrs = {
            "DEXC": -0.003,
            "EATOM": 2024.8347,
            "EAUG": 569.085,
            "ENMAX": 269.865,
            "ENMIN": 202.399,
            "IUNSCR": 1,
            "LCOR": True,
            "LEXCH": "PE",
            "LPAW": True,
            "LULTRA": False,
            "NDATA": 70,
            "POMASS": 54.938,
            "QCUT": -4.454,
            "QGAM": 8.907,
            "RAUG": 1.3,
            "RCLOC": 1.725,
            "RCORE": 2.3,
            "RDEP": 2.338,
            "RMAX": 2.807,
            "RPACOR": 2.08,
            "RWIGS": 1.323,
            "STEP": [25.286, 0.183],
            "TITEL": "PAW_PBE Mn_pv 07Sep2000",
            "VRHFIN": "Mn: 3p4s3d",
            "ZVAL": 13.0,
        }

    def test_keywords(self):
        for key, val in self.Mn_pv_attrs.items():
            assert self.psingle_Mn_pv.keywords[key] == val

        psingle = self.psingle_Fe_54
        data = {
            "nentries": 9,
            "Orbitals": (
                (1, 0, 0.5, -6993.844, 2.0),
                (2, 0, 0.5, -814.6047, 2.0),
                (2, 1, 1.5, -693.3689, 6.0),
                (3, 0, 0.5, -89.4732, 2.0),
                (3, 1, 1.5, -55.6373, 6.0),
                (3, 2, 2.5, -3.8151, 7.0),
                (4, 0, 0.5, -4.2551, 1.0),
                (4, 1, 1.5, -3.4015, 0.0),
                (4, 3, 2.5, -1.3606, 0.0),
            ),
            "OrbitalDescriptions": (
                (2, -3.8151135, 23, 2.3, None, None),
                (2, -5.1756961, 23, 2.3, None, None),
                (0, -4.2550963, 23, 2.3, None, None),
                (0, 7.2035603, 23, 2.3, None, None),
                (1, -2.7211652, 23, 2.3, None, None),
                (1, 18.4316424, 23, 2.3, None, None),
            ),
        }

        for key, val in data.items():
            assert psingle.keywords[key] == val

    def test_nelectrons(self):
        assert self.psingle_Mn_pv.nelectrons == 13
        assert self.psingle_Fe.nelectrons == 8

    def test_electron_configuration(self):
        def assert_config_equal(actual_config, expected_config) -> None:
            """
            Helper function to assert that the electron configuration is equal.
            Each configuration contains: (n: int, l: str, occ: float).
            """
            assert len(actual_config) == len(expected_config), "Configurations have different lengths"

            for expected, actual in zip(expected_config, actual_config, strict=False):
                assert expected[0] == actual[0], f"Principal quantum number mismatch: {expected[0]} != {actual[0]}"
                assert expected[1] == actual[1], f"Subshell mismatch: {expected[1]} != {actual[1]}"

                assert expected[2] == approx(actual[2]), f"Occupation number mismatch: {expected[2]} != {actual[2]}"

        # Test s-block  (Li: 2s1)
        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Li.gz").electron_configuration,
            [
                (2.0, "s", 1.0),
            ],
        )

        # Test p-block  (O: 2s2 sp4)
        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.O.gz").electron_configuration,
            [
                (2, "s", 2.0),
                (2, "p", 4.0),
            ],
        )

        # Test d-block (Fe: 3d7 4s1)
        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe.gz").electron_configuration,
            [(3, "d", 7.0), (4, "s", 1.0)],
        )

        # Test f-block (Ce: 5s2 6s2 5p6 5d1 4f1)
        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Ce.gz").electron_configuration,
            [
                (5, "s", 2),
                (6, "s", 2),
                (5, "p", 6),
                (5, "d", 1),
                (4, "f", 1),
            ],
        )

        # Test "sv" POTCARs (Ca_sv: 3s2 4s2 3p6)
        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Ca_sv.gz").electron_configuration,
            [
                (3, "s", 2),
                (4, "s", 2),
                (3, "p", 6),
            ],
        )

        # Test "pv" POTCARs (Fe_pv: 3p6 3d7 4s1)
        assert PotcarSingle.from_file(
            f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe_pv.gz"
        ).electron_configuration == [
            (3, "p", 6),
            (3, "d", 7),
            (4, "s", 1),
        ]

        # Test non-integer occupancy (Be: 2s1.99 2p0.01)
        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Be.gz").electron_configuration,
            [
                (2, "s", 1.99),
                (2, "p", 0.01),
            ],
        )

        # Test another non-integer occupancy (H.25: 1s0.25)
        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.H.25.gz").electron_configuration,
            [
                (1, "s", 0.25),
            ],
        )

        # Test occupancy tolerance (Be: 2s1.99 2p0.01)
        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Be.gz").get_electron_configuration(
                tol=0.1
            ),
            [
                (2, "s", 1.99),
            ],
        )

        # Test POT_PAW_PBE_64 PSPs
        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_PAW_PBE_64/POTCAR.O.gz").electron_configuration,
            [
                (2, "s", 2.0),
                (2, "p", 4.0),
            ],
        )

        assert_config_equal(
            PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_PAW_PBE_64/POTCAR.Fe_pv.gz").electron_configuration,
            [
                (3, "p", 6.0),
                (3, "d", 7.0),
                (4, "s", 1.0),
            ],
        )

    def test_attributes(self):
        for key, val in self.Mn_pv_attrs.items():
            assert getattr(self.psingle_Mn_pv, key) == val
            assert isinstance(getattr(self.psingle_Fe, key), type(val))

    def test_found_unknown_key(self):
        with pytest.raises(KeyError, match="BAD_KEY"):
            PotcarSingle.parse_functions["BAD_KEY"]

    def test_bad_value(self):
        with pytest.raises(ValueError, match="could not convert string to float"):
            PotcarSingle.parse_functions["ENMAX"]("this should be a float")

    def test_functional_types(self):
        assert self.psingle_Mn_pv.functional == "PBE"
        assert self.psingle_Mn_pv.functional_class == "GGA"
        assert self.psingle_Mn_pv.potential_type == "PAW"

        psingle = PotcarSingle.from_file(f"{FAKE_POTCAR_DIR}/POT_LDA_PAW/POTCAR.Fe.gz")
        assert psingle.functional == "Perdew-Zunger81"
        assert psingle.functional_class == "LDA"
        assert psingle.potential_type == "PAW"
        assert self.psingle_Mn_pv.symbol == "Mn_pv"

    def test_is_valid(self):
        assert self.psingle_Fe.is_valid
        assert self.psingle_Fe_54.is_valid
        assert self.psingle_Mn_pv.is_valid

        # corrupt the file
        psingle = copy.deepcopy(self.psingle_Fe_54)
        assert psingle.keywords["RCORE"] == approx(2.3)
        psingle.keywords["RCORE"] = 2.2
        assert not psingle.is_valid

        psingle = copy.deepcopy(self.psingle_Fe_54)
        psingle.keywords.pop("RCORE")
        assert not psingle.is_valid

        psingle = copy.deepcopy(self.psingle_Fe_54)
        old_data = psingle.data
        psingle.data = psingle.data.replace("RCORE  =    2.3", "RCORE = 2.2")
        assert old_data != psingle.data
        # TODO: should arguably be False but since header is parsed at instantiation time and not reparsed
        # in is_valid, changing the data string in the header section does not currently invalidate POTCAR
        assert psingle.is_valid

        # this POTCAR is valid because the header is only modified in a way that is
        # irrelevant to how FORTRAN reads files, i.e. treated by Fortran as a comment
        filename = f"{FAKE_POTCAR_DIR}/modified_potcars_header/POT_GGA_PAW_PBE/POTCAR.Fe_pv.gz"
        psingle = PotcarSingle.from_file(filename)
        assert psingle.is_valid

    def test_unknown_potcar_warning(self):
        filename = f"{FAKE_POTCAR_DIR}/modified_potcars_data/POT_GGA_PAW_PBE/POTCAR.Fe_pv.gz"
        with pytest.warns(
            UnknownPotcarWarning,
            match="POTCAR data with symbol Fe_pv is not known to pymatgen. ",
        ):
            PotcarSingle.from_file(filename)

    def test_faulty_potcar_has_wrong_hash(self):
        filename = f"{FAKE_POTCAR_DIR}/modified_potcars_data/POT_GGA_PAW_PBE_54/POTCAR.Fe_pv_with_hash.gz"
        psingle = PotcarSingle.from_file(filename)
        assert not psingle.is_valid
        assert psingle.sha256_computed_file_hash != psingle.hash_sha256_from_file

    def test_verify_correct_potcar_with_sha256(self):
        filename = f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe_pv_with_hash.gz"
        psingle = PotcarSingle.from_file(filename)
        assert psingle.sha256_computed_file_hash == psingle.hash_sha256_from_file

    def test_multi_potcar_with_and_without_sha256(self):
        filename = f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe_O.gz"
        potcars = Potcar.from_file(filename)
        # Still need to test the if POTCAR can be read.
        # No longer testing for hashes
        for psingle in potcars:
            if psingle.hash_sha256_from_file:
                assert psingle.sha256_computed_file_hash == psingle.hash_sha256_from_file
            else:
                assert psingle.is_valid

    def test_default_functional(self):
        with patch.dict(SETTINGS, PMG_DEFAULT_FUNCTIONAL="PBE"):
            potcar = PotcarSingle.from_symbol_and_functional("Fe")
            assert potcar.functional_class == "GGA"
        with patch.dict(SETTINGS, PMG_DEFAULT_FUNCTIONAL="LDA"):
            SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "LDA"
            potcar = PotcarSingle.from_symbol_and_functional("Fe")
            assert potcar.functional_class == "LDA"

    def test_from_symbol_and_functional_raises(self):
        # test FileNotFoundError on non-existent PMG_VASP_PSP_DIR in SETTINGS
        PMG_VASP_PSP_DIR = "missing-dir"
        symbol, functional = "Fe", "PBE_64"
        with (  # test PMG_VASP_PSP_DIR not set in SETTINGS
            patch.dict(SETTINGS, PMG_VASP_PSP_DIR=None),
            pytest.raises(
                PmgVaspPspDirError,
                match=re.escape("Set PMG_VASP_PSP_DIR=<directory-path> in .pmgrc.yaml (needed to find POTCARs)"),
            ),
        ):
            PotcarSingle.from_symbol_and_functional(symbol, functional)

        with (  # test FileNotFoundError on non-existent PMG_VASP_PSP_DIR in SETTINGS
            patch.dict(SETTINGS, PMG_VASP_PSP_DIR=PMG_VASP_PSP_DIR),
            pytest.raises(FileNotFoundError, match=f"{PMG_VASP_PSP_DIR=} does not exist."),
        ):
            PotcarSingle.from_symbol_and_functional(symbol, functional)

        # test different FileNotFoundError on non-existent POTCAR sub-directory
        PMG_VASP_PSP_DIR = SETTINGS["PMG_VASP_PSP_DIR"]
        err_msg = f"You do not have the right POTCAR with {functional=} and {symbol=}\nin your {PMG_VASP_PSP_DIR=}"

        with (
            patch.dict(SETTINGS, PMG_VASP_PSP_SUB_DIRS={"PBE_64": "PBE_64_FOO"}),
            pytest.raises(FileNotFoundError, match=re.escape(err_msg)),
        ):
            PotcarSingle.from_symbol_and_functional(symbol, functional)

    def test_repr(self):
        expected_repr = (
            "PotcarSingle(symbol='Mn_pv', functional='PBE', TITEL='PAW_PBE Mn_pv 07Sep2000', "
            "VRHFIN='Mn: 3p4s3d', n_valence_elec=13)"
        )
        assert repr(self.psingle_Mn_pv) == expected_repr

    def test_hash(self):
        assert self.psingle_Mn_pv.md5_header_hash == "b45747d8ceeee91c3b27e8484db32f5a"
        assert self.psingle_Fe.md5_header_hash == "adcc7d2abffa088eccc74948a68235d6"

    def test_potcar_file_hash(self):
        assert self.psingle_Mn_pv.md5_computed_file_hash == "e66e5662ec6e46d6f10ce0bb07b3b742"
        assert self.psingle_Fe.md5_computed_file_hash == "ae761615a0734cc5a2a1db0d5919f12d"

    def test_sha256_file_hash(self):
        assert (
            self.psingle_Mn_pv.sha256_computed_file_hash
            == "3890fe92124e18500817b565a6048a317968613e226ab7b7c2a2d4ca62451e3a"
        )
        assert (
            self.psingle_Fe.sha256_computed_file_hash
            == "7bcf5ad80200e5d74ba63b45d87825b31e6cae2bcd03cebda2f1cbec9870c1cf"
        )

    def test_eq(self):
        assert self.psingle_Mn_pv == self.psingle_Mn_pv
        assert self.psingle_Fe == self.psingle_Fe
        assert self.psingle_Mn_pv != self.psingle_Fe
        assert self.psingle_Mn_pv != self.psingle_Fe_54

    def test_copy(self):
        psingle = self.psingle_Mn_pv.copy()
        assert psingle == self.psingle_Mn_pv
        assert psingle is not self.psingle_Mn_pv

    def test_spec(self):
        for psingle in [self.psingle_Fe, self.psingle_Fe_54, self.psingle_Mn_pv]:
            expected_spec = {
                "titel": psingle.TITEL,
                "hash": psingle.md5_header_hash,
                "summary_stats": psingle._summary_stats,
                "symbol": psingle.symbol,
            }
            assert expected_spec == psingle.spec(extra_spec=["symbol"])


@pytest.mark.filterwarnings(
    "ignore:POTCAR data with symbol .* is not known to pymatgen:pymatgen.io.vasp.inputs.UnknownPotcarWarning"
)
class TestPotcar(MatSciTest):
    def setup_method(self):
        SETTINGS.setdefault("PMG_VASP_PSP_DIR", str(TEST_FILES_DIR))
        self.filepath = f"{FAKE_POTCAR_DIR}/POTCAR.gz"
        self.potcar = Potcar.from_file(self.filepath)

    def test_init(self):
        assert self.potcar.symbols == [
            "Fe",
            "P",
            "O",
        ], "Wrong symbols read in for POTCAR"
        potcar = Potcar(["Fe_pv", "O"])
        assert potcar[0].enmax == approx(293.238)

    def test_from_file(self):
        assert {d.header for d in self.potcar} == {
            "PAW_PBE O 08Apr2002",
            "PAW_PBE P 17Jan2003",
            "PAW_PBE Fe 06Sep2000",
        }

    def test_potcar_map(self):
        fe_potcar = zopen(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE/POTCAR.Fe_pv.gz", mode="rt", encoding="utf-8").read()
        # specify V instead of Fe - this makes sure the test won't pass if the
        # code just grabs the POTCAR from the config file (the config file would
        # grab the V POTCAR)
        potcar = Potcar(["V"], sym_potcar_map={"V": fe_potcar})
        assert potcar.symbols == ["Fe_pv"], "Wrong symbols read in for POTCAR"

    def test_as_from_dict(self):
        dct = self.potcar.as_dict()
        potcar = Potcar.from_dict(dct)
        assert potcar.symbols == ["Fe", "P", "O"]

    def test_write(self):
        tmp_file = f"{self.tmp_path}/POTCAR.testing"
        self.potcar.write_file(tmp_file)
        potcar = Potcar.from_file(tmp_file)
        assert potcar.symbols == self.potcar.symbols

        with (
            zopen(self.filepath, mode="rt", encoding="utf-8") as f_ref,
            open(tmp_file, encoding="utf-8") as f_new,
        ):
            ref_potcar = f_ref.readlines()
            new_potcar = f_new.readlines()

        assert len(ref_potcar) == len(new_potcar), f"wrong POTCAR line count: {len(ref_potcar)} != {len(new_potcar)}"

        # check line by line
        for line1, line2 in zip(ref_potcar, new_potcar, strict=True):
            assert line1.strip() == line2.strip(), f"wrong POTCAR line: {line1} != {line2}"

    def test_set_symbol(self):
        assert self.potcar.symbols == ["Fe", "P", "O"]
        assert self.potcar[0].nelectrons == 8
        self.potcar.symbols = ["Fe_pv", "O"]
        assert self.potcar.symbols == ["Fe_pv", "O"]
        assert self.potcar[0].nelectrons == 14

    @pytest.mark.xfail(reason="TODO: need someone to fix this")
    def test_default_functional(self):
        potcar = Potcar(["Fe", "P"])
        assert potcar[0].functional_class == "GGA"
        assert potcar[1].functional_class == "GGA"
        SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "LDA"
        potcar = Potcar(["Fe", "P"])
        assert potcar[0].functional_class == "LDA"
        assert potcar[1].functional_class == "LDA"

    def test_pickle(self):
        pickle.dumps(self.potcar)

    def test_from_spec(self):
        orig_potcar = Potcar(symbols=["Fe", "P", "O"], functional="PBE")
        new_potcar = Potcar.from_spec(orig_potcar.spec)
        assert str(new_potcar) == str(orig_potcar)
        assert all(
            [PotcarSingle.compare_potcar_stats(p._summary_stats, orig_potcar[ip]._summary_stats)]
            for ip, p in enumerate(new_potcar)
            if p.TITEL == orig_potcar[ip].TITEL
        )

    # def tearDown(self):
    #     SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "PBE"


@pytest.mark.filterwarnings("ignore:Please use INCAR KSPACING tag:DeprecationWarning")
@pytest.mark.filterwarnings(
    "ignore:POTCAR data with symbol .* is not known to pymatgen:pymatgen.io.vasp.inputs.UnknownPotcarWarning"
)
class TestVaspInput(MatSciTest):
    def setup_method(self):
        filepath = f"{VASP_IN_DIR}/INCAR"
        incar = Incar.from_file(filepath)
        filepath = f"{VASP_IN_DIR}/POSCAR"
        poscar = Poscar.from_file(filepath, check_for_potcar=False)
        os.environ.setdefault("PMG_VASP_PSP_DIR", str(TEST_FILES_DIR))
        filepath = f"{FAKE_POTCAR_DIR}/POTCAR.gz"
        potcar = Potcar.from_file(filepath)
        filepath = f"{VASP_IN_DIR}/KPOINTS_auto"
        kpoints = Kpoints.from_file(filepath)
        self.vasp_input = VaspInput(incar, kpoints, poscar, potcar)

    def test_as_from_dict(self):
        dct = self.vasp_input.as_dict()
        vasp_input = VaspInput.from_dict(dct)
        comp = vasp_input["POSCAR"].structure.composition
        assert comp == Composition("Fe4P4O16")

    def test_write(self):
        tmp_dir = f"{self.tmp_path}/VaspInput.testing"
        self.vasp_input.write_input(tmp_dir)

        incar = Incar.from_file(f"{tmp_dir}/INCAR")
        assert incar["NSW"] == 99

        assert {*os.listdir(tmp_dir)} == {"INCAR", "KPOINTS", "POSCAR", "POTCAR"}

    def test_copy(self):
        vasp_input2 = self.vasp_input.copy(deep=True)
        assert isinstance(vasp_input2, VaspInput)
        # make copy and original serialize to the same dict
        assert vasp_input2.as_dict() == self.vasp_input.as_dict()
        # modify the copy and make sure the original is not modified
        vasp_input2["INCAR"]["NSW"] = 100
        assert vasp_input2["INCAR"]["NSW"] == 100, f"{vasp_input2['INCAR']['NSW']=}"
        orig_nsw_val = self.vasp_input["INCAR"]["NSW"]
        assert orig_nsw_val == 99, f"{orig_nsw_val=}"

        # make a shallow copy and make sure the original is modified
        vasp_input3 = self.vasp_input.copy(deep=False)
        vasp_input3["INCAR"]["NSW"] = 100
        assert vasp_input3["INCAR"]["NSW"] == 100, f"{vasp_input3['INCAR']['NSW']=}"
        orig_nsw_val = self.vasp_input["INCAR"]["NSW"]
        assert orig_nsw_val == 99, f"{orig_nsw_val=}"

    def test_run_vasp(self):
        self.vasp_input.run_vasp(".", vasp_cmd=["cat", "INCAR"])
        with open("vasp.out", encoding="utf-8") as file:
            output = file.read()
            assert output.split("\n")[0] == "ALGO = Damped"

    def test_from_directory(self):
        # Previously, this test relied on the existence of a file named POTCAR
        # that was sorted to the top of a list of POTCARs for the test to work.
        # That's far too brittle - isolating requisite files here
        for file in ("INCAR", "KPOINTS", "POSCAR_Li2O"):
            copyfile(f"{VASP_IN_DIR}/{file}", f"{self.tmp_path}/{file.split('_')[0]}")

        Potcar(symbols=["Li_sv", "O"], functional="PBE").write_file(f"{self.tmp_path}/POTCAR")

        copyfile(f"{VASP_OUT_DIR}/CONTCAR_Li2O", f"{self.tmp_path}/CONTCAR_Li2O")

        vi = VaspInput.from_directory(self.tmp_path, optional_files={"CONTCAR_Li2O": Poscar})

        assert vi["INCAR"]["ALGO"] == "Damped"
        assert "CONTCAR_Li2O" in vi

        vi.as_dict()

        vasp_input = VaspInput.from_dict(vi.as_dict())
        assert "CONTCAR_Li2O" in vasp_input

    def test_input_attr(self):
        # test attributes match dict keys
        assert all(v == getattr(self.vasp_input, k.lower()) for k, v in self.vasp_input.items())

        vis_potcar_spec = VaspInput(
            self.vasp_input.incar,
            self.vasp_input.kpoints,
            self.vasp_input.poscar,
            "\n".join(self.vasp_input.potcar.symbols),
            potcar_spec=True,
        )
        # test has expected keys
        assert {*vis_potcar_spec} == {"INCAR", "KPOINTS", "POSCAR", "POTCAR.spec"}
        # test values match
        assert all(self.vasp_input[k] == getattr(vis_potcar_spec, k.lower()) for k in ("INCAR", "KPOINTS", "POSCAR"))
        assert isinstance(vis_potcar_spec.potcar, str)

        # test incar can be updated in place, see https://github.com/materialsproject/pymatgen/issues/4051
        assert vis_potcar_spec.incar["NSW"] == 99
        vis_potcar_spec.incar["NSW"] = 100
        assert vis_potcar_spec.incar["NSW"] == 100


def test_potcar_summary_stats() -> None:
    potcar_summary_stats = loadfn(POTCAR_STATS_PATH)

    assert len(potcar_summary_stats) == 16
    n_potcars_per_functional = {
        "PBE": 251,
        "PBE_52": 303,
        "PBE_54": 326,
        "PBE_64": 343,
        "LDA": 292,
        "LDA_52": 274,
        "LDA_54": 295,
        "PW91": 169,
        "LDA_US": 74,
        "PW91_US": 75,
        "Perdew_Zunger81": 292,
        "PBE_52_W_HASH": 304,
        "PBE_54_W_HASH": 327,
        "LDA_52_W_HASH": 275,
        "LDA_54_W_HASH": 295,
        "LDA_64": 297,
    }
    assert {*potcar_summary_stats} == {*n_potcars_per_functional}

    for key, expected in n_potcars_per_functional.items():
        actual = len(potcar_summary_stats[key])
        assert actual == expected, f"{key=}, {expected=}, {actual=}"


def test_gen_potcar_summary_stats() -> None:
    assert set(_SUMM_STATS) == set(PotcarSingle.functional_dir)

    expected_funcs = [x for x in os.listdir(str(FAKE_POTCAR_DIR)) if x in PotcarSingle.functional_dir]

    for func in expected_funcs:
        bdir = f"{FAKE_POTCAR_DIR}/{PotcarSingle.functional_dir[func]}"
        valid_elements = [x for x in os.listdir(f"{bdir}") if x[0] != "." and os.path.isdir(f"{bdir}/{x}")]
        for element in valid_elements:
            assert PotcarSingle.from_file(f"{bdir}/POTCAR.{element}.gz").is_valid
