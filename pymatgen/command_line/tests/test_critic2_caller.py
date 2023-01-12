# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import os
import unittest
from shutil import which

from pymatgen.command_line.critic2_caller import Critic2Analysis, Critic2Caller
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

__author__ = "Matthew Horton"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Production"
__date__ = "July 2017"


@unittest.skipIf(not which("critic2"), "critic2 executable not present")
class Critic2CallerTest(unittest.TestCase):
    def test_from_path(self):
        # uses chgcars
        test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "bader")

        c2c = Critic2Caller.from_path(test_dir)

        # check we have some results!
        self.assertGreaterEqual(len(c2c._stdout), 500)

        c2o = c2c.output

        # check we get our structure graph
        sg = c2o.structure_graph(include_critical_points=None)
        self.assertEqual(sg.get_coordination_of_site(0), 4)

        # check yt integration
        self.assertAlmostEqual(c2o.structure.site_properties["bader_volume"][0], 66.0148355)
        self.assertAlmostEqual(c2o.structure.site_properties["bader_charge"][0], 12.2229131)

        # test zpsp functionality
        # this is normally picked up from POTCARs, but since POTCARs not checked in with the
        # test suite, setting manually here
        c2o_dict = c2o.as_dict()
        c2o_dict["zpsp"] = {"Fe": 8.0, "O": 6.0}
        c2o = Critic2Analysis.from_dict(c2o_dict)
        # note: these values don't seem sensible physically, but seem to be correct with
        # respect to the input files (possibly bad/underconverged source data)
        self.assertAlmostEqual(c2o.structure.site_properties["bader_charge_transfer"][0], 4.2229131)

        # alternatively, can also set when we do the analysis, but note that this will change
        # the analysis performed since augmentation charges are added in core regions
        c2c = Critic2Caller.from_path(test_dir, zpsp={"Fe": 8.0, "O": 6.0})

        # check yt integration
        self.assertAlmostEqual(c2o.structure.site_properties["bader_volume"][0], 66.0148355)
        self.assertAlmostEqual(c2o.structure.site_properties["bader_charge"][0], 12.2229131)
        self.assertAlmostEqual(c2o.structure.site_properties["bader_charge_transfer"][0], 4.2229131)

    def test_from_structure(self):
        # uses promolecular density
        structure = Structure.from_file(
            os.path.join(
                PymatgenTest.TEST_FILES_DIR,
                "critic2/MoS2.cif",
            )
        )

        c2c = Critic2Caller.from_chgcar(structure)

        # check we have some results!
        self.assertGreaterEqual(len(c2c._stdout), 500)


class Critic2AnalysisTest(unittest.TestCase):
    _multiprocess_shared_ = True

    def setUp(self):
        stdout_file = os.path.join(PymatgenTest.TEST_FILES_DIR, "critic2/MoS2_critic2_stdout.txt")
        stdout_file_new_format = os.path.join(PymatgenTest.TEST_FILES_DIR, "critic2/MoS2_critic2_stdout_new_format.txt")
        with open(stdout_file) as f:
            reference_stdout = f.read()
        with open(stdout_file_new_format) as f:
            reference_stdout_new_format = f.read()

        structure = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "critic2/MoS2.cif"))

        self.c2o = Critic2Analysis(structure, reference_stdout)
        self.c2o_new_format = Critic2Analysis(structure, reference_stdout_new_format)

    def test_properties_to_from_dict(self):

        self.assertEqual(len(self.c2o.critical_points), 6)
        self.assertEqual(len(self.c2o.nodes), 14)
        self.assertEqual(len(self.c2o.edges), 10)

        self.assertEqual(len(self.c2o_new_format.critical_points), 6)
        self.assertEqual(len(self.c2o_new_format.nodes), 14)
        self.assertEqual(len(self.c2o_new_format.edges), 10)

        # reference dictionary for c2o.critical_points[0].as_dict()
        # {'@class': 'CriticalPoint',
        #  '@module': 'pymatgen.command_line.critic2_caller',
        #  'coords': None,
        #  'field': 93848.0413,
        #  'field_gradient': 0.0,
        #  'field_hessian': [[-2593274446000.0, -3.873587547e-19, -1.704530713e-08],
        #                    [-3.873587547e-19, -2593274446000.0, 1.386877485e-18],
        #                    [-1.704530713e-08, 1.386877485e-18, -2593274446000.0]],
        #  'frac_coords': [0.333333, 0.666667, 0.213295],
        #  'index': 0,
        #  'multiplicity': 1.0,
        #  'point_group': 'D3h',
        #  'type': < CriticalPointType.nucleus: 'nucleus' >}

        self.assertEqual(str(self.c2o.critical_points[0].type), "CriticalPointType.nucleus")

        # test connectivity
        self.assertDictEqual(
            self.c2o.edges[3],
            {"from_idx": 1, "from_lvec": (0, 0, 0), "to_idx": 0, "to_lvec": (1, 0, 0)},
        )
        # test as/from dict
        d = self.c2o.as_dict()
        self.assertEqual(
            set(d),
            {
                "@module",
                "@class",
                "@version",
                "structure",
                "stdout",
                "stderr",
                "cpreport",
                "yt",
                "zpsp",
            },
        )
        self.c2o.from_dict(d)

    def test_graph_output(self):

        sg = self.c2o.structure_graph()
        self.assertEqual(str(sg.structure[3].specie), "Xbcp")
        self.assertSetEqual(
            set(list(sg.graph.edges(data=True))[0][2]),
            {"to_jimage", "weight", "field", "laplacian", "ellipticity", "frac_coords"},
        )


if __name__ == "__main__":

    unittest.main()
