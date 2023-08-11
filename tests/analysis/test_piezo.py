"""Test for the piezo tensor class."""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from pymatgen.analysis.piezo import PiezoTensor
from pymatgen.util.testing import PymatgenTest

__author__ = "Shyam Dwaraknath"
__version__ = "0.1"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "4/1/16"


class TestPiezo(PymatgenTest):
    def setUp(self):
        self.piezo_struc = self.get_structure("BaNiO3")
        self.voigt_matrix = np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.03839, 0.0],
                [0.0, 0.0, 0.0, 0.03839, 0.0, 0.0],
                [6.89822, 6.89822, 27.46280, 0.0, 0.0, 0.0],
            ]
        )
        self.vasp_matrix = np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.03839],
                [0.0, 0.0, 0.0, 0.0, 0.03839, 0.0],
                [6.89822, 6.89822, 27.46280, 0.0, 0.0, 0.0],
            ]
        )
        self.full_tensor_array = [
            [[0.0, 0.0, 0.03839], [0.0, 0.0, 0.0], [0.03839, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.03839], [0.0, 0.03839, 0.0]],
            [[6.89822, 0.0, 0.0], [0.0, 6.89822, 0.0], [0.0, 0.0, 27.4628]],
        ]

    def test_new(self):
        pt = PiezoTensor(self.full_tensor_array)
        assert np.allclose(pt, self.full_tensor_array)
        bad_dim_array = np.zeros((3, 3))
        with pytest.raises(ValueError, match="PiezoTensor input must be rank 3"):
            PiezoTensor(bad_dim_array)

    def test_from_voigt(self):
        bad_voigt = np.zeros((3, 7))
        pt = PiezoTensor.from_voigt(self.voigt_matrix)
        assert_array_equal(pt, self.full_tensor_array)
        with pytest.raises(ValueError, match="Invalid shape for Voigt matrix"):
            PiezoTensor.from_voigt(bad_voigt)
        assert_array_equal(self.voigt_matrix, pt.voigt)

    def test_from_vasp_voigt(self):
        bad_voigt = np.zeros((3, 7))
        pt = PiezoTensor.from_vasp_voigt(self.vasp_matrix)
        assert_array_equal(pt, self.full_tensor_array)
        with pytest.raises(ValueError, match="Invalid shape for Voigt matrix"):
            PiezoTensor.from_voigt(bad_voigt)
        assert_array_equal(self.voigt_matrix, pt.voigt)
