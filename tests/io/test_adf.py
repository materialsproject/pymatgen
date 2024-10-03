from __future__ import annotations

from pytest import approx

from pymatgen.core.structure import Molecule
from pymatgen.io.adf import AdfInput, AdfKey, AdfOutput, AdfTask
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "Xin Chen, chenxin13@mails.tsinghua.edu.cn"


TEST_DIR = f"{TEST_FILES_DIR}/io/adf"

geometry_string = """GEOMETRY
smooth conservepoints
optim all cartesian
iterations 250
step rad=0.15 angle=10.0
hessupd BFGS
converge e=0.001 grad=0.0003 rad=0.01 angle=0.5
END
"""

zlm_fit_string = """ZLMFIT
AtomDepQuality
10 good
12 normal
subend
END
"""

atoms_string = """ATOMS
O      -0.90293455        0.66591421        0.00000000
H       0.05706545        0.66591421        0.00000000
H      -1.22338913        1.57085004        0.00000000
END
"""

h2o_xyz = """3
0.0
O -0.90293455 0.66591421 0.0
H  0.05706545 0.66591421 0.0
H -1.22338913 1.57085004 0.0
"""

rhb18_xyz = """19
0.0
Rh       -0.453396   -0.375115    0.000000
B         0.168139    3.232791    0.000000
B        -0.270938    1.639058    0.000000
B         0.206283    2.604044    1.459430
B         0.404410    1.880136    2.866764
B        -0.103309    0.887485    1.655272
B         0.436856    0.371367    3.299887
B         0.016593   -0.854959    1.930982
B         0.563233   -1.229713    3.453066
B         0.445855   -2.382027    2.415013
B         0.206283    2.604044   -1.459430
B         0.404410    1.880136   -2.866764
B        -0.103309    0.887485   -1.655272
B         0.436856    0.371367   -3.299887
B         0.563233   -1.229713   -3.453066
B         0.016593   -0.854959   -1.930982
B         0.200456   -2.309538   -0.836316
B         0.200456   -2.309538    0.836316
B         0.445855   -2.382027   -2.415013
"""


def readfile(file_object):
    """`
    Return the content of the file as a string.

    Args:
        file_object (file or str): The file to read. This can be either a File object or a file path.

    Returns:
        content (str): The content of the file.
    """
    if hasattr(file_object, "read"):
        return file_object.read()
    with open(file_object) as file:
        return file.read()


class TestAdfKey:
    def test_simple(self):
        unrestricted = AdfKey("unrestricted")
        assert str(unrestricted).strip() == "UNRESTRICTED"

    def test_options(self):
        charge = AdfKey("charge", [-1, 0])
        charge_string = "CHARGE -1 0\n"
        assert str(charge) == "CHARGE -1 0\n"
        assert str(AdfKey.from_dict(charge.as_dict())) == charge_string

    def test_subkeys(self):
        smooth = AdfKey("smooth", ["conservepoints"])
        optim = AdfKey("optim", ["all", "cartesian"])
        iterations = AdfKey("iterations", [250])
        step = AdfKey("step", [("rad", 0.15), ("angle", 10.0)])
        hessupd = AdfKey("hessupd", ["BFGS"])
        converge = AdfKey(
            "converge",
            [("e", 1.0e-3), ("grad", 3.0e-4), ("rad", 1.0e-2), ("angle", 0.5)],
        )
        geo = AdfKey("geometry", subkeys=[smooth, optim, iterations, step, hessupd, converge])
        assert str(geo) == geometry_string
        assert str(AdfKey.from_dict(geo.as_dict())) == geometry_string
        assert geo.has_subkey("optim")

    def test_end(self):
        geo = AdfKey("Geometry")
        assert str(geo) == "GEOMETRY\nEND\n"

    def test_subkeys_subkeys(self):
        atom_dep_quality = AdfKey("AtomDepQuality", subkeys=[AdfKey("10", ["good"]), AdfKey("12", ["normal"])])
        zlmfit = AdfKey("zlmfit", subkeys=[atom_dep_quality])
        assert str(zlmfit) == zlm_fit_string
        assert str(AdfKey.from_dict(zlmfit.as_dict())) == zlm_fit_string

    def test_from_str(self):
        k1 = AdfKey.from_str("CHARGE -1 0")
        assert k1.key == "CHARGE"
        assert k1.options == [-1, 0]

        k2 = AdfKey.from_str("step rad=0.15 angle=10.0")
        assert k2.key == "step"
        assert k2.options[0] == ["rad", 0.15]
        assert k2.options[1] == ["angle", 10.0]

        k3 = AdfKey.from_str("GEOMETRY\noptim all\niterations 100\nEND\n")
        assert k3.key == "GEOMETRY"
        assert k3.subkeys[0].options[0] == "all"
        assert k3.subkeys[1].options[0] == 100

        k4 = AdfKey.from_str(
            """SCF
            iterations 300
            converge 1.0e-7 1.0e-7
            mixing 0.2
            diis n=100 ok=0.0001 cyc=100 cx=5.0 cxx=10.0
            END"""
        )
        assert k4.key == "SCF"
        assert k4.subkeys[0].key == "iterations"
        assert k4.subkeys[1].key == "converge"
        assert k4.subkeys[1].options[0] == 1e-7
        assert k4.subkeys[2].options[0] == 0.2

    def test_option_operations(self):
        k1 = AdfKey("Charge", [-1, 0])
        k1.add_option(2)
        assert k1.options == [-1, 0, 2]
        k1.remove_option(0)
        assert k1.options == [0, 2]

        k2 = AdfKey.from_str("step rad=0.15 angle=10.0")
        k2.add_option(["length", 0.1])
        assert k2.options[2] == ["length", 0.1]
        k2.remove_option("rad")
        assert k2.options[0] == ["angle", 10.0]

    def test_atom_block_key(self):
        block = AdfKey("atoms")
        mol = Molecule.from_str(h2o_xyz, "xyz")
        for site in mol:
            block.add_subkey(AdfKey(str(site.specie), list(site.coords)))
        assert str(block) == atoms_string


energy_task = """TITLE ADF_RUN

UNITS
length angstrom
angle degree
END

XC
GGA PBE
END

BASIS
type DZ
core small
END

SCF
iterations 300
END

GEOMETRY SinglePoint
END

"""


class TestAdfTask:
    def test_energy(self):
        task = AdfTask()
        assert str(task) == energy_task

    def test_serialization(self):
        task = AdfTask()
        adf_task = AdfTask.from_dict(task.as_dict())
        assert task.title == adf_task.title
        assert task.basis_set == adf_task.basis_set
        assert task.scf == adf_task.scf
        assert task.geo == adf_task.geo
        assert task.operation == adf_task.operation
        assert task.units == adf_task.units
        assert str(task) == str(adf_task)


rhb18 = {
    "title": "RhB18",
    "basis_set": AdfKey.from_str("BASIS\ntype TZP\ncore small\nEND"),
    "xc": AdfKey.from_str("XC\nHybrid PBE0\nEND"),
    "units": AdfKey.from_str("UNITS\nlength angstrom\nEND"),
    "other_directives": [
        AdfKey.from_str("SYMMETRY"),
        AdfKey.from_str("RELATIVISTIC scalar zora"),
        AdfKey.from_str("INTEGRATION 6.0 6.0 6.0"),
        AdfKey.from_str("SAVE TAPE21"),
        AdfKey.from_str("A1FIT 10.0"),
    ],
    "geo_subkeys": [
        AdfKey.from_str("optim all"),
        AdfKey.from_str("iterations 300"),
        AdfKey.from_str("step rad=0.15 angle=10.0"),
        AdfKey.from_str("hessupd BFGS"),
    ],
    "scf": AdfKey.from_str(
        """SCF
             iterations 300
             converge 1.0e-7 1.0e-7
             mixing 0.2
             lshift 0.0
             diis n=100 ok=0.0001 cyc=100 cx=5.0 cxx=10.0
             END"""
    ),
}


class TestAdfInput(PymatgenTest):
    def test_main(self):
        tmp_file = f"{self.tmp_path}/adf.temp"
        mol = Molecule.from_str(rhb18_xyz, "xyz")
        mol.set_charge_and_spin(-1, 3)
        task = AdfTask("optimize", **rhb18)
        inp = AdfInput(task)
        inp.write_file(mol, tmp_file)
        expected = readfile(f"{TEST_DIR}/RhB18_adf.inp")
        assert readfile(tmp_file) == expected


class TestAdfOutput:
    def test_analytical_freq(self):
        filename = f"{TEST_DIR}/analytical_freq/adf.out"
        adf_out = AdfOutput(filename)
        assert adf_out.final_energy == approx(-0.54340325)
        assert len(adf_out.energies) == 4
        assert len(adf_out.structures) == 4
        assert adf_out.frequencies[0] == approx(1553.931)
        assert adf_out.frequencies[2] == approx(3793.086)
        assert adf_out.normal_modes[0][2] == approx(0.071)
        assert adf_out.normal_modes[0][6] == approx(0.000)
        assert adf_out.normal_modes[0][7] == approx(-0.426)
        assert adf_out.normal_modes[0][8] == approx(-0.562)

    def test_numerical_freq(self):
        filename = f"{TEST_DIR}/numerical_freq/adf.out"
        adf_out = AdfOutput(filename)
        assert adf_out.freq_type == "Numerical"
        assert len(adf_out.final_structure) == 4
        assert len(adf_out.frequencies) == 6
        assert len(adf_out.normal_modes) == 6
        assert adf_out.frequencies[0] == approx(938.21)
        assert adf_out.frequencies[3] == approx(3426.64)
        assert adf_out.frequencies[4] == approx(3559.35)
        assert adf_out.frequencies[5] == approx(3559.35)
        assert adf_out.normal_modes[1][0] == approx(0.067)
        assert adf_out.normal_modes[1][3] == approx(-0.536)
        assert adf_out.normal_modes[1][7] == approx(0.000)
        assert adf_out.normal_modes[1][9] == approx(-0.536)

    def test_single_point(self):
        filename = f"{TEST_DIR}/sp/adf.out"
        adf_out = AdfOutput(filename)
        assert adf_out.final_energy == approx(-0.74399276)
        assert len(adf_out.final_structure) == 4
