# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import warnings
from pathlib import Path

import numpy as np
from plotly.graph_objects import Figure

from pymatgen.analysis.chempot_diagram import (
    ChemicalPotentialDiagram,
    simple_pca,
    get_centroid_2d,
    get_2d_orthonormal_vector,
)
from pymatgen.core.composition import Element
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.util.testing import PymatgenTest

module_dir = Path(__file__).absolute().parent


class ChemicalPotentialDiagramTest(PymatgenTest):
    def setUp(self):
        self.entries = EntrySet.from_csv(str(module_dir / "pdentries_test.csv"))
        self.cpd_ternary = ChemicalPotentialDiagram(entries=self.entries, default_min_limit=-25)
        elements = [Element("Fe"), Element("O")]
        binary_entries = list(
            filter(
                lambda e: set(e.composition.elements).issubset(elements),
                self.entries,
            )
        )
        self.cpd_binary = ChemicalPotentialDiagram(entries=binary_entries, default_min_limit=-25)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_dim(self):
        self.assertEqual(self.cpd_binary.dim, 2)
        self.assertEqual(self.cpd_ternary.dim, 3)

    def test_el_refs(self):
        el_refs = {elem: entry.energy for elem, entry in self.cpd_ternary.el_refs.items()}

        elems = [Element("Li"), Element("Fe"), Element("O")]
        energies = [-1.91301487, -6.5961471, -25.54966885]
        correct_el_refs = dict(zip(elems, energies))

        self.assertDictsAlmostEqual(el_refs, correct_el_refs)

    def test_border_hyperplanes(self):
        desired = np.array(
            [[-1, 0, 0, -25], [1, 0, 0, 0], [0, -1, 0, -25], [0, 1, 0, 0], [0, 0, -1, -25], [0, 0, 1, 0]]
        )
        self.assertArrayAlmostEqual(self.cpd_ternary.border_hyperplanes, desired)

    def test_lims(self):
        desired_lims = np.array([[-25, 0], [-25, 0], [-25, 0]])
        self.assertArrayAlmostEqual(self.cpd_ternary.lims, desired_lims)

    def test_pca(self):
        points_3d = np.array(
            [
                [-25.0, -6.5961471, -7.11535414],
                [-25.0, -6.74159386, -6.96990738],
                [-4.07706195, -6.74159386, -6.96990738],
                [-3.93161519, -6.5961471, -7.11535414],
            ]
        )

        points_2d_desired = np.array(
            [
                [10.49782722, 0.10320265],
                [10.4978342, -0.10249014],
                [-10.42510384, -0.10320018],
                [-10.57055758, 0.10248767],
            ]
        )

        points_2d, _, _ = simple_pca(points_3d, k=2)

        self.assertArrayAlmostEqual(points_2d, points_2d_desired)

    def test_centroid(self):
        vertices = np.array(
            [
                [3.30046498, 0.1707431],
                [-0.63480672, 0.21578376],
                [-1.37717499, 0.13347465],
                [-1.61148314, 0.0409329],
                [-1.77880975, -0.25825963],
                [2.10180962, -0.30267477],
            ]
        )

        centroid = get_centroid_2d(vertices)
        centroid_desired = np.array([-0.00069433, -0.00886174])

        self.assertArrayAlmostEqual(centroid, centroid_desired)

    def test_get_2d_orthonormal_vector(self):
        pts_1 = np.array([[1, 1], [2, 2]])
        pts_2 = np.array([[-2, -5], [-4, 6]])

        vec_1 = get_2d_orthonormal_vector(pts_1)
        vec_2 = get_2d_orthonormal_vector(pts_2)

        vec_1_desired = np.array([0.70710678, 0.70710678])
        vec_2_desired = np.array([0.98386991, 0.17888544])

        self.assertArrayAlmostEqual(vec_1, vec_1_desired)
        self.assertArrayAlmostEqual(vec_2, vec_2_desired)

    def test_get_plot(self):
        fig_2d = self.cpd_binary.get_plot()
        fig_3d = self.cpd_ternary.get_plot()

        self.assertEqual(type(fig_2d), Figure)
        self.assertEqual(fig_2d["data"][0]["type"], "scatter")

        self.assertEqual(type(fig_3d), Figure)
        self.assertEqual(fig_3d["data"][0]["type"], "scatter3d")

    def test_domains(self):
        correct_domains = {
            "Fe3O4": np.array(
                [
                    [-25.0, -7.29639011, -6.55381019],
                    [-25.0, -6.74159386, -6.96990738],
                    [-4.35446008, -7.29639011, -6.55381019],
                    [-4.07706195, -6.74159386, -6.96990738],
                ]
            ),
            "LiFeO2": np.array(
                [
                    [-4.35446008, -7.29639011, -6.55381019],
                    [-4.07706195, -6.74159386, -6.96990738],
                    [-3.93161519, -6.5961471, -7.11535414],
                    [-3.62500226, -6.5961471, -7.2686606],
                    [-4.66220934, -9.70776834, -5.19424644],
                    [-5.40627456, -10.45183356, -4.45018123],
                ]
            ),
            "Fe2O3": np.array(
                [
                    [-25.0, -10.73968818, -4.25827814],
                    [-25.0, -7.29639011, -6.55381019],
                    [-4.35446008, -7.29639011, -6.55381019],
                    [-5.55020187, -10.73968818, -4.25827814],
                    [-5.40627456, -10.45183356, -4.45018123],
                ]
            ),
            "Li2FeO3": np.array(
                [
                    [-4.73988663, -10.25150936, -4.96121458],
                    [-4.66220934, -9.70776834, -5.19424644],
                    [-5.44282307, -10.95444579, -4.25827814],
                    [-5.55020187, -10.73968818, -4.25827814],
                    [-5.40627456, -10.45183356, -4.45018123],
                ]
            ),
            "Li5FeO4": np.array(
                [
                    [-4.61251054, -10.37888545, -5.08859067],
                    [-3.35159776, -6.5961471, -7.61041623],
                    [-3.62500226, -6.5961471, -7.2686606],
                    [-4.73988663, -10.25150936, -4.96121458],
                    [-4.66220934, -9.70776834, -5.19424644],
                ]
            ),
            "O2": np.array(
                [
                    [-25.0, -10.73968818, -4.25827814],
                    [-25.0, -25.0, -4.25827814],
                    [-5.44282307, -25.0, -4.25827814],
                    [-5.44282307, -10.95444579, -4.25827814],
                    [-5.55020187, -10.73968818, -4.25827814],
                ]
            ),
            "FeO": np.array(
                [
                    [-25.0, -6.5961471, -7.11535414],
                    [-25.0, -6.74159386, -6.96990738],
                    [-4.07706195, -6.74159386, -6.96990738],
                    [-3.93161519, -6.5961471, -7.11535414],
                ]
            ),
            "Li2O2": np.array(
                [
                    [-5.44282307, -25.0, -4.25827814],
                    [-4.61251054, -25.0, -5.08859067],
                    [-4.61251054, -10.37888545, -5.08859067],
                    [-4.73988663, -10.25150936, -4.96121458],
                    [-5.44282307, -10.95444579, -4.25827814],
                ]
            ),
            "Li2O": np.array(
                [
                    [-1.91301487, -25.0, -10.48758201],
                    [-1.91301487, -6.5961471, -10.48758201],
                    [-4.61251054, -25.0, -5.08859067],
                    [-4.61251054, -10.37888545, -5.08859067],
                    [-3.35159776, -6.5961471, -7.61041623],
                ]
            ),
            "Li": np.array(
                [
                    [-1.91301487, -6.5961471, -25.0],
                    [-1.91301487, -25.0, -25.0],
                    [-1.91301487, -25.0, -10.48758201],
                    [-1.91301487, -6.5961471, -10.48758201],
                ]
            ),
            "Fe": np.array(
                [
                    [-25.0, -6.5961471, -25.0],
                    [-1.91301487, -6.5961471, -25.0],
                    [-25.0, -6.5961471, -7.11535414],
                    [-1.91301487, -6.5961471, -10.48758201],
                    [-3.93161519, -6.5961471, -7.11535414],
                    [-3.35159776, -6.5961471, -7.61041623],
                    [-3.62500226, -6.5961471, -7.2686606],
                ]
            ),
        }
        for formula, domain in correct_domains.items():
            self.assertArrayAlmostEqual(domain, self.cpd_ternary.domains[formula])


if __name__ == "__main__":
    unittest.main()
