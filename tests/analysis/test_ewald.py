from __future__ import annotations

from unittest import TestCase

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.ewald import EwaldMinimizer, EwaldSummation
from pymatgen.core.structure import Structure
from pymatgen.util.testing import VASP_IN_DIR


class TestEwaldSummation(TestCase):
    def setUp(self):
        filepath = f"{VASP_IN_DIR}/POSCAR"
        self.original_struct = Structure.from_file(filepath)
        self.struct = self.original_struct.copy()
        self.struct.add_oxidation_state_by_element({"Li": 1, "Fe": 2, "P": 5, "O": -2})

    def test_init(self):
        ham = EwaldSummation(self.struct, compute_forces=True)
        assert ham.real_space_energy == approx(-502.23549897772602, abs=1e-4)
        assert ham.reciprocal_space_energy == approx(6.1541071599534654, abs=1e-4)
        assert ham.point_energy == approx(-620.22598358035918, abs=1e-4)
        assert ham.total_energy == approx(-1123.00766, abs=1e-1)
        assert ham.forces[0, 0] == approx(-1.98818620e-01, abs=1e-4)
        assert sum(sum(abs(ham.forces))) == approx(915.925354346, abs=1e-4), "Forces incorrect"
        assert sum(sum(ham.real_space_energy_matrix)) == approx(ham.real_space_energy, abs=1e-4)
        assert sum(sum(ham.reciprocal_space_energy_matrix)) == approx(ham.reciprocal_space_energy, abs=1e-4)
        assert sum(ham.point_energy_matrix) == approx(ham.point_energy, abs=1e-4)
        assert sum(sum(ham.total_energy_matrix)) + ham._charged_cell_energy == approx(ham.total_energy, abs=1e-2)

        with pytest.raises(
            ValueError,
            match="Ewald summation can only be performed on structures that are either oxidation state decorated",
        ):
            EwaldSummation(self.original_struct)
        # try sites with charge.
        charges = []
        for site in self.original_struct:
            if site.specie.symbol == "Li":
                charges.append(1)
            elif site.specie.symbol == "Fe":
                charges.append(2)
            elif site.specie.symbol == "P":
                charges.append(5)
            else:
                charges.append(-2)

        self.original_struct.add_site_property("charge", charges)
        ham2 = EwaldSummation(self.original_struct)
        assert ham2.real_space_energy == approx(-502.23549897772602, abs=1e-4)

    def test_from_dict(self):
        ham = EwaldSummation(self.struct, compute_forces=True)
        ham2 = EwaldSummation.from_dict(ham.as_dict())
        assert ham._real is None
        assert not ham._initialized
        assert ham2._real is None
        assert not ham2._initialized
        assert np.array_equal(ham.total_energy_matrix, ham2.total_energy_matrix)
        # check lazy eval
        assert ham.total_energy == approx(-1123.00766, abs=1e-1)
        assert ham._real is not None
        assert ham._initialized
        ham2 = EwaldSummation.from_dict(ham.as_dict())
        assert ham2._real is not None
        assert ham2._initialized
        assert np.array_equal(ham.total_energy_matrix, ham2.total_energy_matrix)

    def test_as_dict(self):
        ham = EwaldSummation(self.struct, compute_forces=True)
        dct = ham.as_dict()
        assert dct["compute_forces"]
        assert dct["eta"] == ham._eta
        assert dct["acc_factor"] == ham._acc_factor
        assert dct["real_space_cut"] == ham._rmax
        assert dct["recip_space_cut"] == ham._gmax
        assert ham.as_dict() == EwaldSummation.from_dict(dct).as_dict()


class TestEwaldMinimizer(TestCase):
    def test_init(self):
        matrix = np.array(
            [
                [-3.0, 3.0, 4.0, -0.0, 3.0, 3.0, 1.0, 14.0, 9.0, -4.0],
                [1.0, -3.0, -3.0, 12.0, -4.0, -1.0, 5.0, 11.0, 1.0, 12.0],
                [14.0, 7.0, 13.0, 15.0, 13.0, 5.0, -5.0, 10.0, 14.0, -2.0],
                [9.0, 13.0, 4.0, 1.0, 3.0, -4.0, 7.0, 0.0, 6.0, -4.0],
                [4.0, -4.0, 6.0, 1.0, 12.0, -4.0, -2.0, 13.0, 0.0, 6.0],
                [13.0, 7.0, -4.0, 12.0, -2.0, 9.0, 8.0, -5.0, 3.0, 1.0],
                [8.0, 1.0, 10.0, -4.0, -2.0, 4.0, 13.0, 12.0, -3.0, 13.0],
                [2.0, 11.0, 8.0, 1.0, -1.0, 5.0, -3.0, 4.0, 5.0, 0.0],
                [-0.0, 14.0, 4.0, 3.0, -1.0, -5.0, 7.0, -1.0, -1.0, 3.0],
                [2.0, -2.0, 10.0, 1.0, 6.0, -5.0, -3.0, 12.0, 0.0, 13.0],
            ]
        )

        m_list = [[0.9, 4, [1, 2, 3, 4, 8], "a"], [-1, 2, [5, 6, 7], "b"]]

        e_min = EwaldMinimizer(matrix, m_list, 50)

        assert len(e_min.output_lists) == 15, "Wrong number of permutations returned"
        assert e_min.minimized_sum == approx(111.63, abs=1e-3), "Returned wrong minimum value"
        assert len(e_min.best_m_list) == 6, "Returned wrong number of permutations"

    def test_site(self):
        """Test that uses an uncharged structure."""
        filepath = f"{VASP_IN_DIR}/POSCAR"
        struct = Structure.from_file(filepath)
        struct = struct.copy()
        struct.add_oxidation_state_by_element({"Li": 1, "Fe": 3, "P": 5, "O": -2})

        # Comparison to LAMMPS result
        ham = EwaldSummation(struct, compute_forces=True)
        assert approx(ham.total_energy, abs=1e-3) == -1226.3335
        assert approx(ham.get_site_energy(0), abs=1e-3) == -45.8338
        assert approx(ham.get_site_energy(8), abs=1e-3) == -27.2978
