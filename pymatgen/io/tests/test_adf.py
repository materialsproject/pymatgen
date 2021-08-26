import os
import unittest
from os.path import join

from pymatgen.core.structure import Molecule
from pymatgen.io.adf import AdfInput, AdfKey, AdfOutput, AdfTask
from pymatgen.util.testing import PymatgenTest

__author__ = "Xin Chen, chenxin13@mails.tsinghua.edu.cn"

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules")

geometry_string = """GEOMETRY
smooth conservepoints
optim all cartesian
iterations 250
step rad=0.15 angle=10.0
hessupd BFGS
converge e=0.001 grad=0.0003 rad=0.01 angle=0.5
END
"""

zlmfit_string = """ZLMFIT
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

h2oxyz = """3
0.0
O -0.90293455 0.66591421 0.0
H  0.05706545 0.66591421 0.0
H -1.22338913 1.57085004 0.0
"""

rhb18xyz = """19
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
    """
    Return the content of the file as a string.

    Parameters
    ----------
    file_object : file or str
        The file to read. This can be either a File object or a file path.

    Returns
    -------
    content : str
        The content of the file.

    """
    if hasattr(file_object, "read"):
        return file_object.read()
    elif isinstance(file_object, str):
        f = open(file_object, "r")
        content = f.read()
        f.close()
        return content
    else:
        raise ValueError("``file_object`` must be a string or a file object!")


class AdfKeyTest(unittest.TestCase):
    def test_simple(self):
        unrestricted = AdfKey("unrestricted")
        self.assertEqual(str(unrestricted).strip(), "UNRESTRICTED")

    def test_options(self):
        charge = AdfKey("charge", [-1, 0])
        charge_string = "CHARGE -1 0\n"
        self.assertEqual(str(charge), "CHARGE -1 0\n")
        self.assertEqual(str(AdfKey.from_dict(charge.as_dict())), charge_string)

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
        self.assertEqual(str(geo), geometry_string)
        self.assertEqual(str(AdfKey.from_dict(geo.as_dict())), geometry_string)
        self.assertTrue(geo.has_subkey("optim"))

    def test_end(self):
        geo = AdfKey("Geometry")
        self.assertEqual(str(geo), "GEOMETRY\nEND\n")

    def test_subkeys_subkeys(self):
        atom_dep_quality = AdfKey("AtomDepQuality", subkeys=[AdfKey("10", ["good"]), AdfKey("12", ["normal"])])
        zlmfit = AdfKey("zlmfit", subkeys=[atom_dep_quality])
        self.assertEqual(str(zlmfit), zlmfit_string)
        self.assertEqual(str(AdfKey.from_dict(zlmfit.as_dict())), zlmfit_string)

    def test_from_string(self):
        k1 = AdfKey.from_string("CHARGE -1 0")
        self.assertEqual(k1.key, "CHARGE")
        self.assertListEqual(k1.options, [-1, 0])

        k2 = AdfKey.from_string("step rad=0.15 angle=10.0")
        self.assertEqual(k2.key, "step")
        self.assertListEqual(k2.options[0], ["rad", 0.15])
        self.assertListEqual(k2.options[1], ["angle", 10.0])

        k3 = AdfKey.from_string("GEOMETRY\noptim all\niterations 100\nEND\n")
        self.assertEqual(k3.key, "GEOMETRY")
        self.assertEqual(k3.subkeys[0].options[0], "all")
        self.assertEqual(k3.subkeys[1].options[0], 100)

        k4 = AdfKey.from_string(
            """SCF
            iterations 300
            converge 1.0e-7 1.0e-7
            mixing 0.2
            diis n=100 ok=0.0001 cyc=100 cx=5.0 cxx=10.0
            END"""
        )
        self.assertEqual(k4.key, "SCF")
        self.assertEqual(k4.subkeys[0].key, "iterations")
        self.assertEqual(k4.subkeys[1].key, "converge")
        self.assertEqual(k4.subkeys[1].options[0], 1e-7)
        self.assertEqual(k4.subkeys[2].options[0], 0.2)

    def test_option_operations(self):
        k1 = AdfKey("Charge", [-1, 0])
        k1.add_option(2)
        self.assertListEqual(k1.options, [-1, 0, 2])
        k1.remove_option(0)
        self.assertListEqual(k1.options, [0, 2])

        k2 = AdfKey.from_string("step rad=0.15 angle=10.0")
        k2.add_option(["length", 0.1])
        self.assertListEqual(k2.options[2], ["length", 0.1])
        k2.remove_option("rad")
        self.assertListEqual(k2.options[0], ["angle", 10.0])

    def test_atom_block_key(self):
        block = AdfKey("atoms")
        o = Molecule.from_str(h2oxyz, "xyz")
        for site in o:
            block.add_subkey(AdfKey(str(site.specie), list(site.coords)))
        self.assertEqual(str(block), atoms_string)


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


class AdfTaskTest(unittest.TestCase):
    def test_energy(self):
        task = AdfTask()
        self.assertEqual(str(task), energy_task)

    def test_serialization(self):
        task = AdfTask()
        o = AdfTask.from_dict(task.as_dict())
        self.assertEqual(task.title, o.title)
        self.assertEqual(task.basis_set, o.basis_set)
        self.assertEqual(task.scf, o.scf)
        self.assertEqual(task.geo, o.geo)
        self.assertEqual(task.operation, o.operation)
        self.assertEqual(task.units, o.units)
        self.assertEqual(str(task), str(o))


rhb18 = {
    "title": "RhB18",
    "basis_set": AdfKey.from_string("BASIS\ntype TZP\ncore small\nEND"),
    "xc": AdfKey.from_string("XC\nHybrid PBE0\nEND"),
    "units": AdfKey.from_string("UNITS\nlength angstrom\nEND"),
    "other_directives": [
        AdfKey.from_string("SYMMETRY"),
        AdfKey.from_string("RELATIVISTIC scalar zora"),
        AdfKey.from_string("INTEGRATION 6.0 6.0 6.0"),
        AdfKey.from_string("SAVE TAPE21"),
        AdfKey.from_string("A1FIT 10.0"),
    ],
    "geo_subkeys": [
        AdfKey.from_string("optim all"),
        AdfKey.from_string("iterations 300"),
        AdfKey.from_string("step rad=0.15 angle=10.0"),
        AdfKey.from_string("hessupd BFGS"),
    ],
    "scf": AdfKey.from_string(
        """SCF
             iterations 300
             converge 1.0e-7 1.0e-7
             mixing 0.2
             lshift 0.0
             diis n=100 ok=0.0001 cyc=100 cx=5.0 cxx=10.0
             END"""
    ),
}


class AdfInputTest(unittest.TestCase):
    def setUp(self):
        self.tempfile = "./adf.temp"

    def test_main(self):
        o = Molecule.from_str(rhb18xyz, "xyz")
        o.set_charge_and_spin(-1, 3)
        task = AdfTask("optimize", **rhb18)
        inp = AdfInput(task)
        inp.write_file(o, self.tempfile)
        s = readfile(join(test_dir, "adf", "RhB18_adf.inp"))
        self.assertEqual(readfile(self.tempfile), s)

    def tearDown(self):
        if os.path.isfile(self.tempfile):
            os.remove(self.tempfile)


class AdfOutputTest(unittest.TestCase):
    def test_analytical_freq(self):
        filename = join(test_dir, "adf", "analytical_freq", "adf.out")
        o = AdfOutput(filename)
        self.assertAlmostEqual(o.final_energy, -0.54340325)
        self.assertEqual(len(o.energies), 4)
        self.assertEqual(len(o.structures), 4)
        self.assertAlmostEqual(o.frequencies[0], 1553.931)
        self.assertAlmostEqual(o.frequencies[2], 3793.086)
        self.assertAlmostEqual(o.normal_modes[0][2], 0.071)
        self.assertAlmostEqual(o.normal_modes[0][6], 0.000)
        self.assertAlmostEqual(o.normal_modes[0][7], -0.426)
        self.assertAlmostEqual(o.normal_modes[0][8], -0.562)

    def test_numerical_freq(self):
        filename = join(test_dir, "adf", "numerical_freq", "adf.out")
        o = AdfOutput(filename)
        self.assertEqual(o.freq_type, "Numerical")
        self.assertEqual(o.final_structure.num_sites, 4)
        self.assertEqual(len(o.frequencies), 6)
        self.assertEqual(len(o.normal_modes), 6)
        self.assertAlmostEqual(o.frequencies[0], 938.21)
        self.assertAlmostEqual(o.frequencies[3], 3426.64)
        self.assertAlmostEqual(o.frequencies[4], 3559.35)
        self.assertAlmostEqual(o.frequencies[5], 3559.35)
        self.assertAlmostEqual(o.normal_modes[1][0], 0.067)
        self.assertAlmostEqual(o.normal_modes[1][3], -0.536)
        self.assertAlmostEqual(o.normal_modes[1][7], 0.000)
        self.assertAlmostEqual(o.normal_modes[1][9], -0.536)

    def test_single_point(self):
        filename = join(test_dir, "adf", "sp", "adf.out")
        o = AdfOutput(filename)
        self.assertAlmostEqual(o.final_energy, -0.74399276)
        self.assertEqual(len(o.final_structure), 4)


if __name__ == "__main__":
    unittest.main()
