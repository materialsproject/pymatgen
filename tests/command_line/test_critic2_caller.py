from __future__ import annotations

from shutil import which
from unittest import TestCase

import pytest
from pytest import approx

from pymatgen.command_line.critic2_caller import Critic2Analysis, Critic2Caller
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR

__author__ = "Matthew Horton"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Production"
__date__ = "July 2017"

TEST_DIR = f"{TEST_FILES_DIR}/command_line/critic2"
BADER_TEST_DIR = f"{TEST_FILES_DIR}/command_line/bader"


@pytest.mark.skipif(not which("critic2"), reason="critic2 executable not present")
class TestCritic2Caller:
    def test_from_path(self):
        # uses CHGCARs
        c2c = Critic2Caller.from_path(BADER_TEST_DIR)

        # check we have some results!
        assert len(c2c._stdout) >= 500

        c2o = c2c.output

        # check we get our structure graph
        sg = c2o.structure_graph(include_critical_points=None)
        assert sg.get_coordination_of_site(0) == 4

        # check yt integration
        assert c2o.structure.site_properties["bader_volume"][0] == approx(66.0148355)
        assert c2o.structure.site_properties["bader_charge"][0] == approx(12.2229131)

        # test zpsp functionality
        # this is normally picked up from POTCARs, but since POTCARs not checked in with the
        # test suite, setting manually here
        c2o_dict = c2o.as_dict()
        c2o_dict["zpsp"] = {"Fe": 8.0, "O": 6.0}
        c2o = Critic2Analysis.from_dict(c2o_dict)
        # note: these values don't seem sensible physically, but seem to be correct with
        # respect to the input files (possibly bad/underconverged source data)
        assert c2o.structure.site_properties["bader_charge_transfer"][0] == approx(4.2229131)

        # alternatively, can also set when we do the analysis, but note that this will change
        # the analysis performed since augmentation charges are added in core regions
        c2c = Critic2Caller.from_path(BADER_TEST_DIR, zpsp={"Fe": 8.0, "O": 6.0})

        # check yt integration
        assert c2o.structure.site_properties["bader_volume"][0] == approx(66.0148355)
        assert c2o.structure.site_properties["bader_charge"][0] == approx(12.2229131)
        assert c2o.structure.site_properties["bader_charge_transfer"][0] == approx(4.2229131)

    def test_from_structure(self):
        # uses pro-molecular density
        structure = Structure.from_file(f"{TEST_FILES_DIR}/cif/MoS2.cif")

        c2c = Critic2Caller.from_chgcar(structure)

        # check we have some results!
        assert len(c2c._stdout) >= 500

        # test with chgcar and zpsp to ensure zval is formatted as int
        # https://github.com/materialsproject/pymatgen/issues/3501
        c2c = Critic2Caller.from_chgcar(structure, zpsp={"Mo": 6.0, "S": 6.0}, chgcar=f"{BADER_TEST_DIR}/CHGCAR.gz")

        assert "ERROR : load int.CHGCAR id chg_int zpsp Mo 6 S 6" in c2c._input_script


class TestCritic2Analysis(TestCase):
    def setUp(self):
        stdout_file = f"{TEST_DIR}/MoS2_critic2_stdout.txt"
        stdout_file_new_format = f"{TEST_DIR}/MoS2_critic2_stdout_new_format.txt"
        with open(stdout_file) as file:
            reference_stdout = file.read()
        with open(stdout_file_new_format) as file:
            reference_stdout_new_format = file.read()

        structure = Structure.from_file(f"{TEST_DIR}/MoS2.cif")

        self.c2o = Critic2Analysis(structure, reference_stdout)
        self.c2o_new_format = Critic2Analysis(structure, reference_stdout_new_format)

    def test_to_from_dict(self):
        """
        reference dictionary for c2o.critical_points[0].as_dict()
        {'@class': 'CriticalPoint',
         '@module': 'pymatgen.command_line.critic2_caller',
         'coords': None,
         'field': 93848.0413,
         'field_gradient': 0.0,
         'field_hessian': [[-2593274446000.0, -3.873587547e-19, -1.704530713e-08],
                           [-3.873587547e-19, -2593274446000.0, 1.386877485e-18],
                           [-1.704530713e-08, 1.386877485e-18, -2593274446000.0]],
         'frac_coords': [0.333333, 0.666667, 0.213295],
         'index': 0,
         'multiplicity': 1.0,
         'point_group': 'D3h',
         'type': < CriticalPointType.nucleus: 'nucleus' >}
        """
        assert len(self.c2o.critical_points) == 6
        assert len(self.c2o.nodes) == 14
        assert len(self.c2o.edges) == 10

        assert len(self.c2o_new_format.critical_points) == 6
        assert len(self.c2o_new_format.nodes) == 14
        assert len(self.c2o_new_format.edges) == 10

        assert str(self.c2o.critical_points[0].type) == "CriticalPointType.nucleus"

        # test connectivity
        assert self.c2o.edges[3] == {"from_idx": 1, "from_lvec": (0, 0, 0), "to_idx": 0, "to_lvec": (1, 0, 0)}
        # test as/from dict
        dct = self.c2o.as_dict()
        assert set(dct) == {
            "@module",
            "@class",
            "@version",
            "structure",
            "stdout",
            "stderr",
            "cpreport",
            "yt",
            "zpsp",
        }
        self.c2o.from_dict(dct)

    def test_graph_output(self):
        sg = self.c2o.structure_graph()
        assert str(sg.structure[3].specie) == "Xbcp"
        assert set(next(iter(sg.graph.edges(data=True)))[2]) == {
            "to_jimage",
            "weight",
            "field",
            "laplacian",
            "ellipticity",
            "frac_coords",
        }
