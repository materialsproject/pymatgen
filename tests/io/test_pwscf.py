from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.io.pwscf import PWInput, PWInputError, PWOutput
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestPWInput(PymatgenTest):
    def test_init(self):
        struct = self.get_structure("Li2O")
        with pytest.raises(PWInputError, match="Missing O2- in pseudo specification"):
            PWInput(
                struct,
                control={"calculation": "scf", "pseudo_dir": "./"},
                pseudo={"Li": "Li.pbe-n-kjpaw_psl.0.1.UPF"},
            )

    def test_str_mixed_oxidation(self):
        struct = self.get_structure("Li2O")
        struct.remove_oxidation_states()
        struct[1] = "Li1"
        pw = PWInput(
            struct,
            control={"calculation": "scf", "pseudo_dir": "./"},
            pseudo={
                "Li": "Li.pbe-n-kjpaw_psl.0.1.UPF",
                "Li+": "Li.pbe-n-kjpaw_psl.0.1.UPF",
                "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
            },
            system={"ecutwfc": 50},
        )
        expected = """&CONTROL
  calculation = 'scf',
  pseudo_dir = './',
/
&SYSTEM
  ecutwfc = 50,
  ibrav = 0,
  nat = 3,
  ntyp = 3,
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Li  6.9410 Li.pbe-n-kjpaw_psl.0.1.UPF
  Li+  6.9410 Li.pbe-n-kjpaw_psl.0.1.UPF
  O  15.9994 O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
  O 0.000000 0.000000 0.000000
  Li+ 0.750178 0.750178 0.750178
  Li 0.249822 0.249822 0.249822
K_POINTS automatic
  1 1 1 0 0 0
CELL_PARAMETERS angstrom
  2.917389 0.097894 1.520005
  0.964634 2.755036 1.520005
  0.133206 0.097894 3.286918
"""
        assert str(pw).strip() == expected.strip()

    def test_str_without_oxidation(self):
        struct = self.get_structure("Li2O")
        struct.remove_oxidation_states()
        pw = PWInput(
            struct,
            control={"calculation": "scf", "pseudo_dir": "./"},
            pseudo={
                "Li": "Li.pbe-n-kjpaw_psl.0.1.UPF",
                "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
            },
            system={"ecutwfc": 50},
        )
        expected = """&CONTROL
  calculation = 'scf',
  pseudo_dir = './',
/
&SYSTEM
  ecutwfc = 50,
  ibrav = 0,
  nat = 3,
  ntyp = 2,
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Li  6.9410 Li.pbe-n-kjpaw_psl.0.1.UPF
  O  15.9994 O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
  O 0.000000 0.000000 0.000000
  Li 0.750178 0.750178 0.750178
  Li 0.249822 0.249822 0.249822
K_POINTS automatic
  1 1 1 0 0 0
CELL_PARAMETERS angstrom
  2.917389 0.097894 1.520005
  0.964634 2.755036 1.520005
  0.133206 0.097894 3.286918
"""
        assert str(pw).strip() == expected.strip()

    def test_str_with_oxidation(self):
        struct = self.get_structure("Li2O")

        pw = PWInput(
            struct,
            control={"calculation": "scf", "pseudo_dir": "./"},
            pseudo={
                "Li+": "Li.pbe-n-kjpaw_psl.0.1.UPF",
                "O2-": "O.pbe-n-kjpaw_psl.0.1.UPF",
            },
            system={"ecutwfc": 50},
        )
        expected = """&CONTROL
  calculation = 'scf',
  pseudo_dir = './',
/
&SYSTEM
  ecutwfc = 50,
  ibrav = 0,
  nat = 3,
  ntyp = 2,
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Li+  6.9410 Li.pbe-n-kjpaw_psl.0.1.UPF
  O2-  15.9994 O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
  O2- 0.000000 0.000000 0.000000
  Li+ 0.750178 0.750178 0.750178
  Li+ 0.249822 0.249822 0.249822
K_POINTS automatic
  1 1 1 0 0 0
CELL_PARAMETERS angstrom
  2.917389 0.097894 1.520005
  0.964634 2.755036 1.520005
  0.133206 0.097894 3.286918
"""
        assert str(pw).strip() == expected.strip()

    def test_write_str_with_kpoints(self):
        struct = self.get_structure("Li2O")
        struct.remove_oxidation_states()
        kpoints = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.0], [0.0, 0.0, 0.5], [0.5, 0.5, 0.5]]
        pw = PWInput(
            struct,
            control={"calculation": "scf", "pseudo_dir": "./"},
            pseudo={
                "Li": "Li.pbe-n-kjpaw_psl.0.1.UPF",
                "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
            },
            system={"ecutwfc": 50},
            kpoints_mode="crystal_b",
            kpoints_grid=kpoints,
        )
        expected = """
&CONTROL
  calculation = 'scf',
  pseudo_dir = './',
/
&SYSTEM
  ecutwfc = 50,
  ibrav = 0,
  nat = 3,
  ntyp = 2,
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Li  6.9410 Li.pbe-n-kjpaw_psl.0.1.UPF
  O  15.9994 O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
  O 0.000000 0.000000 0.000000
  Li 0.750178 0.750178 0.750178
  Li 0.249822 0.249822 0.249822
K_POINTS crystal_b
 5
 0.0000 0.0000 0.0000
 0.0000 0.5000 0.5000
 0.5000 0.0000 0.0000
 0.0000 0.0000 0.5000
 0.5000 0.5000 0.5000
CELL_PARAMETERS angstrom
  2.917389 0.097894 1.520005
  0.964634 2.755036 1.520005
  0.133206 0.097894 3.286918
"""
        assert str(pw).strip() == expected.strip()

    def test_proc_val(self):
        inputs = {
            "degauss": ("7.3498618000d-03", 7.3498618000e-03),
            "nat": ("2", 2),
            "nosym": (".TRUE.", True),
            "smearing": ("'cold'", "cold"),
        }
        for key, (input_str, expected) in inputs.items():
            value = PWInput.proc_val(key, input_str)
            assert value == expected

    def test_read_str(self):
        string = """
&CONTROL
  calculation = 'scf'
  pseudo_dir = './'
  wf_collect = .TRUE.
/
&SYSTEM
  ibrav = 0,
  nat = 53
  ntyp = 2
  input_dft = 'PBE'
  ecutwfc = 80
  nspin = 1
  nbnd = 280
  smearing = 'cold'
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Mg  24.3050 Mg_ONCV_PBE-1.2.upf
  O  15.9994 O_ONCV_PBE-1.2.upf
ATOMIC_POSITIONS crystal
Mg      -0.000000000   0.000000000  -0.000000000
Mg       0.000000000   1.000000000   0.333134366
Mg       0.000000000   1.000000000   0.666865634
Mg      -0.000000000   0.333134366   1.000000000
Mg       0.000037606   0.333320465   0.333320465
Mg       0.000000000   0.333134366   0.666865634
Mg       0.000000000   0.666865634   1.000000000
Mg       0.000000000   0.666865634   0.333134366
Mg      -0.000037606   0.666679535   0.666679535
Mg       0.333134366   0.000000000   0.000000000
Mg       0.333320465   0.000037606   0.333320465
Mg       0.333134366   1.000000000   0.666865634
Mg       0.333320465   0.333320465   0.000037606
Mg       0.333320465   0.333320465   0.333320465
Mg       0.331436170   0.331436170   0.668563830
Mg       0.333134366   0.666865634  -0.000000000
Mg       0.331436170   0.668563830   0.331436170
Mg       0.331436170   0.668563830   0.668563830
Mg       0.666865634   0.000000000   0.000000000
Mg       0.666865634   0.000000000   0.333134366
Mg       0.666679535  -0.000037606   0.666679535
Mg       0.666865634   0.333134366  -0.000000000
Mg       0.668563830   0.331436170   0.331436170
Mg       0.668563830   0.331436170   0.668563830
Mg       0.666679535   0.666679535  -0.000037606
Mg       0.668563830   0.668563830   0.331436170
Mg       0.666679535   0.666679535   0.666679535
O        0.166588534   0.166588534   0.166588534
O        0.166588534   0.166588534   0.500235399
O        0.166465543   0.166465543   0.833534457
O        0.166588534   0.500235399   0.166588534
O        0.166169242   0.500000000   0.500000000
O        0.166169242   0.500000000   0.833830758
O        0.166465543   0.833534457   0.166465543
O        0.166169242   0.833830758   0.500000000
O        0.166465543   0.833534457   0.833534457
O        0.500235399   0.166588534   0.166588534
O        0.500000000   0.166169242   0.500000000
O        0.500000000   0.166169242   0.833830758
O        0.500000000   0.500000000   0.166169242
O        0.500000000   0.500000000   0.833830758
O        0.500000000   0.833830758   0.166169242
O        0.500000000   0.833830758   0.500000000
O        0.499764601   0.833411466   0.833411466
O        0.833534457   0.166465543   0.166465543
O        0.833830758   0.166169242   0.500000000
O        0.833534457   0.166465543   0.833534457
O        0.833830758   0.500000000   0.166169242
O        0.833830758   0.500000000   0.500000000
O        0.833411466   0.499764601   0.833411466
O        0.833534457   0.833534457   0.166465543
O        0.833411466   0.833411466   0.499764601
O        0.833411466   0.833411466   0.833411466
K_POINTS gamma
CELL_PARAMETERS angstrom
  0.000000 6.373854 6.373854
  6.373854 0.000000 6.373854
  6.373854 6.373854 0.000000
        """
        lattice = np.array([[0.0, 6.373854, 6.373854], [6.373854, 0.0, 6.373854], [6.373854, 6.373854, 0.0]])

        sites = np.array(
            [
                [0.0, 0.0, 0.0],
                [8.49720381, 2.12334981, 6.373854],
                [10.62435819, 4.25050419, 6.373854],
                [8.49720381, 6.373854, 2.12334981],
                [4.24907196, 2.12477567, 2.12477567],
                [6.373854, 4.25050419, 2.12334981],
                [10.62435819, 6.373854, 4.25050419],
                [6.373854, 2.12334981, 4.25050419],
                [8.49863604, 4.24907833, 4.24907833],
                [0.0, 2.12334981, 2.12334981],
                [2.12477567, 4.24907196, 2.12477567],
                [10.62435819, 6.373854, 8.49720381],
                [2.12477567, 2.12477567, 4.24907196],
                [4.24907196, 4.24907196, 4.24907196],
                [6.373854, 6.373854, 4.22505152],
                [4.25050419, 2.12334981, 6.373854],
                [6.373854, 4.22505152, 6.373854],
                [8.52265648, 6.373854, 6.373854],
                [0.0, 4.25050419, 4.25050419],
                [2.12334981, 6.373854, 4.25050419],
                [4.24907833, 8.49863604, 4.24907833],
                [2.12334981, 4.25050419, 6.373854],
                [4.22505152, 6.373854, 6.373854],
                [6.373854, 8.52265648, 6.373854],
                [4.24907833, 4.24907833, 8.49863604],
                [6.373854, 6.373854, 8.52265648],
                [8.49863604, 8.49863604, 8.49863604],
                [2.12362199, 2.12362199, 2.12362199],
                [4.25023839, 4.25023839, 2.12362199],
                [6.373854, 6.373854, 2.12205413],
                [4.25023839, 2.12362199, 4.25023839],
                [6.373854, 4.24606549, 4.24606549],
                [8.50164251, 6.373854, 4.24606549],
                [6.373854, 2.12205413, 6.373854],
                [8.50164251, 4.24606549, 6.373854],
                [10.62565387, 6.373854, 6.373854],
                [2.12362199, 4.25023839, 4.25023839],
                [4.24606549, 6.373854, 4.24606549],
                [6.373854, 8.50164251, 4.24606549],
                [4.24606549, 4.24606549, 6.373854],
                [8.50164251, 8.50164251, 6.373854],
                [6.373854, 4.24606549, 8.50164251],
                [8.50164251, 6.373854, 8.50164251],
                [10.62408601, 8.49746961, 8.49746961],
                [2.12205413, 6.373854, 6.373854],
                [4.24606549, 8.50164251, 6.373854],
                [6.373854, 10.62565387, 6.373854],
                [4.24606549, 6.373854, 8.50164251],
                [6.373854, 8.50164251, 8.50164251],
                [8.49746961, 10.62408601, 8.49746961],
                [6.373854, 6.373854, 10.62565387],
                [8.49746961, 8.49746961, 10.62408601],
                [10.62408601, 10.62408601, 10.62408601],
            ]
        )

        pwin = PWInput.from_str(string)

        # generate list of coords
        pw_sites = np.array([list(site.coords) for site in pwin.structure])

        assert_allclose(sites, pw_sites)

        assert_allclose(lattice, pwin.structure.lattice.matrix)
        assert pwin.sections["system"]["smearing"] == "cold"


class TestPWOuput(PymatgenTest):
    def setUp(self):
        self.pwout = PWOutput(f"{TEST_FILES_DIR}/Si.pwscf.out")

    def test_properties(self):
        assert self.pwout.final_energy == approx(-93.45259708)

    def test_get_celldm(self):
        assert self.pwout.get_celldm(1) == approx(10.323)
        for i in range(2, 7):
            assert self.pwout.get_celldm(i) == approx(0)
