from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
from plotly.graph_objects import Figure
from pytest import approx

from pymatgen.analysis.chempot_diagram import (
    ChemicalPotentialDiagram,
    get_2d_orthonormal_vector,
    get_centroid_2d,
    simple_pca,
)
from pymatgen.core.composition import Element
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.util.testing import PymatgenTest

module_dir = Path(__file__).absolute().parent


class TestChemicalPotentialDiagram(PymatgenTest):
    def setUp(self):
        self.entries = EntrySet.from_csv(str(module_dir / "pd_entries_test.csv"))
        self.cpd_ternary = ChemicalPotentialDiagram(entries=self.entries, default_min_limit=-25, formal_chempots=False)
        self.cpd_ternary_formal = ChemicalPotentialDiagram(
            entries=self.entries, default_min_limit=-25, formal_chempots=True
        )
        elements = [Element("Fe"), Element("O")]
        binary_entries = list(
            filter(
                lambda e: set(e.elements).issubset(elements),
                self.entries,
            )
        )
        self.cpd_binary = ChemicalPotentialDiagram(entries=binary_entries, default_min_limit=-25, formal_chempots=False)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_dim(self):
        assert self.cpd_binary.dim == 2
        assert self.cpd_ternary.dim == 3
        assert self.cpd_ternary_formal.dim == 3

    def test_el_refs(self):
        el_refs = {elem: entry.energy for elem, entry in self.cpd_ternary.el_refs.items()}

        elems = [Element("Li"), Element("Fe"), Element("O")]
        energies = [-1.91301487, -6.5961471, -25.54966885]
        correct_el_refs = dict(zip(elems, energies))

        assert el_refs == approx(correct_el_refs)

    def test_el_refs_formal(self):
        el_refs = {elem: entry.energy for elem, entry in self.cpd_ternary_formal.el_refs.items()}
        elems = [Element("Li"), Element("Fe"), Element("O")]
        energies = [0, 0, 0]
        correct_el_refs = dict(zip(elems, energies))
        assert el_refs == approx(correct_el_refs)

    def test_border_hyperplanes(self):
        desired = np.array(
            [[-1, 0, 0, -25], [1, 0, 0, 0], [0, -1, 0, -25], [0, 1, 0, 0], [0, 0, -1, -25], [0, 0, 1, 0]]
        )
        assert self.cpd_ternary.border_hyperplanes == approx(desired)
        assert self.cpd_ternary_formal.border_hyperplanes == approx(desired)

    def test_lims(self):
        desired_lims = np.array([[-25, 0], [-25, 0], [-25, 0]])
        assert self.cpd_ternary.lims == approx(desired_lims)
        assert self.cpd_ternary_formal.lims == approx(desired_lims)

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

        assert points_2d == approx(points_2d_desired)

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

        assert centroid == approx(centroid_desired, abs=1e-6)

    def test_get_2d_orthonormal_vector(self):
        pts_1 = np.array([[1, 1], [2, 2]])
        pts_2 = np.array([[-2, -5], [-4, 6]])

        vec_1 = get_2d_orthonormal_vector(pts_1)
        vec_2 = get_2d_orthonormal_vector(pts_2)

        vec_1_desired = np.array([0.70710678, 0.70710678])
        vec_2_desired = np.array([0.98386991, 0.17888544])

        assert vec_1 == approx(vec_1_desired)
        assert vec_2 == approx(vec_2_desired)

    def test_get_plot(self):
        fig_2d = self.cpd_binary.get_plot()
        fig_3d = self.cpd_ternary.get_plot()
        fig_3d_formal = self.cpd_ternary_formal.get_plot()

        assert isinstance(fig_2d, Figure)
        assert fig_2d.data[0].type == "scatter"

        assert isinstance(fig_3d, Figure)
        assert fig_3d.data[0].type == "scatter3d"

        assert isinstance(fig_3d_formal, Figure)
        assert fig_3d_formal.data[0].type == "scatter3d"
        assert fig_3d_formal.data[0].mode == "lines"
        assert fig_3d_formal.layout.plot_bgcolor == "rgba(0,0,0,0)"
        assert fig_3d_formal.layout.scene.annotations[0].text == "FeO"

    def test_domains(self):
        correct_domains = {
            "Fe": np.array(
                [
                    [-25.0, -6.596147, -25.0],
                    [-25.0, -6.596147, -7.115354],
                    [-3.931615, -6.596147, -7.115354],
                    [-3.625002, -6.596147, -7.268661],
                    [-3.351598, -6.596147, -7.610416],
                    [-1.913015, -6.596147, -25.0],
                    [-1.913015, -6.596147, -10.487582],
                ]
            ),
            "Fe2O3": np.array(
                [
                    [-25.0, -10.739688, -4.258278],
                    [-25.0, -7.29639, -6.55381],
                    [-5.550202, -10.739688, -4.258278],
                    [-5.406275, -10.451834, -4.450181],
                    [-4.35446, -7.29639, -6.55381],
                ]
            ),
            "Fe3O4": np.array(
                [
                    [-25.0, -7.29639, -6.55381],
                    [-25.0, -6.741594, -6.969907],
                    [-4.35446, -7.29639, -6.55381],
                    [-4.077062, -6.741594, -6.969907],
                ]
            ),
            "FeO": np.array(
                [
                    [-25.0, -6.741594, -6.969907],
                    [-25.0, -6.596147, -7.115354],
                    [-4.077062, -6.741594, -6.969907],
                    [-3.931615, -6.596147, -7.115354],
                ]
            ),
            "Li": np.array(
                [
                    [-1.913015, -25.0, -25.0],
                    [-1.913015, -25.0, -10.487582],
                    [-1.913015, -6.596147, -25.0],
                    [-1.913015, -6.596147, -10.487582],
                ]
            ),
            "Li2FeO3": np.array(
                [
                    [-5.550202, -10.739688, -4.258278],
                    [-5.442823, -10.954446, -4.258278],
                    [-5.406275, -10.451834, -4.450181],
                    [-4.739887, -10.251509, -4.961215],
                    [-4.662209, -9.707768, -5.194246],
                ]
            ),
            "Li2O": np.array(
                [
                    [-4.612511, -25.0, -5.088591],
                    [-4.612511, -10.378885, -5.088591],
                    [-3.351598, -6.596147, -7.610416],
                    [-1.913015, -25.0, -10.487582],
                    [-1.913015, -6.596147, -10.487582],
                ]
            ),
            "Li2O2": np.array(
                [
                    [-5.442823, -25.0, -4.258278],
                    [-5.442823, -10.954446, -4.258278],
                    [-4.739887, -10.251509, -4.961215],
                    [-4.612511, -25.0, -5.088591],
                    [-4.612511, -10.378885, -5.088591],
                ]
            ),
            "Li5FeO4": np.array(
                [
                    [-4.739887, -10.251509, -4.961215],
                    [-4.662209, -9.707768, -5.194246],
                    [-4.612511, -10.378885, -5.088591],
                    [-3.625002, -6.596147, -7.268661],
                    [-3.351598, -6.596147, -7.610416],
                ]
            ),
            "LiFeO2": np.array(
                [
                    [-5.406275, -10.451834, -4.450181],
                    [-4.662209, -9.707768, -5.194246],
                    [-4.35446, -7.29639, -6.55381],
                    [-4.077062, -6.741594, -6.969907],
                    [-3.931615, -6.596147, -7.115354],
                    [-3.625002, -6.596147, -7.268661],
                ]
            ),
            "O2": np.array(
                [
                    [-25.0, -25.0, -4.258278],
                    [-25.0, -10.739688, -4.258278],
                    [-5.550202, -10.739688, -4.258278],
                    [-5.442823, -25.0, -4.258278],
                    [-5.442823, -10.954446, -4.258278],
                ]
            ),
        }

        for formula, domain in correct_domains.items():
            d = self.cpd_ternary.domains[formula]
            d = d.round(6)  # to get rid of numerical errors from qhull
            actual_domain_sorted = d[np.lexsort((d[:, 2], d[:, 1], d[:, 0]))]
            assert actual_domain_sorted == approx(domain)

        formal_domains = {
            "FeO": np.array(
                [
                    [-2.50000000e01, 3.55271368e-15, -2.85707600e00],
                    [-2.01860032e00, 3.55271368e-15, -2.85707600e00],
                    [-2.50000000e01, -1.45446765e-01, -2.71162923e00],
                    [-2.16404709e00, -1.45446765e-01, -2.71162923e00],
                ]
            ),
            "Fe2O3": np.array(
                [
                    [-25.0, -4.14354109, 0.0],
                    [-3.637187, -4.14354108, 0.0],
                    [-3.49325969, -3.85568646, -0.19190308],
                    [-25.0, -0.70024301, -2.29553205],
                    [-2.44144521, -0.70024301, -2.29553205],
                ]
            ),
            "Fe3O4": np.array(
                [
                    [-25.0, -0.70024301, -2.29553205],
                    [-25.0, -0.14544676, -2.71162923],
                    [-2.44144521, -0.70024301, -2.29553205],
                    [-2.16404709, -0.14544676, -2.71162923],
                ]
            ),
            "LiFeO2": np.array(
                [
                    [-3.49325969e00, -3.85568646e00, -1.91903083e-01],
                    [-2.01860032e00, 3.55271368e-15, -2.85707600e00],
                    [-2.44144521e00, -7.00243005e-01, -2.29553205e00],
                    [-2.16404709e00, -1.45446765e-01, -2.71162923e00],
                    [-1.71198739e00, 3.55271368e-15, -3.01038246e00],
                    [-2.74919447e00, -3.11162124e00, -9.35968300e-01],
                ]
            ),
            "Li2O": np.array(
                [
                    [0.00000000e00, -2.50000000e01, -6.22930387e00],
                    [-2.69949567e00, -2.50000000e01, -8.30312528e-01],
                    [3.55271368e-15, 3.55271368e-15, -6.22930387e00],
                    [-1.43858289e00, 3.55271368e-15, -3.35213809e00],
                    [-2.69949567e00, -3.78273835e00, -8.30312528e-01],
                ]
            ),
            "Li2O2": np.array(
                [
                    [-3.52980820e00, -2.50000000e01, 0.00000000e00],
                    [-2.69949567e00, -2.50000000e01, -8.30312528e-01],
                    [-3.52980820e00, -4.35829869e00, 3.55271368e-15],
                    [-2.69949567e00, -3.78273835e00, -8.30312528e-01],
                    [-2.82687176e00, -3.65536226e00, -7.02936437e-01],
                ]
            ),
            "Li2FeO3": np.array(
                [
                    [-3.52980820e00, -4.35829869e00, 3.55271368e-15],
                    [-3.63718700e00, -4.14354108e00, 0.00000000e00],
                    [-3.49325969e00, -3.85568646e00, -1.91903083e-01],
                    [-2.74919447e00, -3.11162124e00, -9.35968300e-01],
                    [-2.82687176e00, -3.65536226e00, -7.02936437e-01],
                ]
            ),
            "Li5FeO4": np.array(
                [
                    [-1.43858289e00, 3.55271368e-15, -3.35213809e00],
                    [-1.71198739e00, 3.55271368e-15, -3.01038246e00],
                    [-2.74919447e00, -3.11162124e00, -9.35968300e-01],
                    [-2.69949567e00, -3.78273835e00, -8.30312528e-01],
                    [-2.82687176e00, -3.65536226e00, -7.02936437e-01],
                ]
            ),
            "O2": np.array(
                [
                    [-2.50000000e01, -2.50000000e01, 3.55271368e-15],
                    [-3.52980820e00, -2.50000000e01, 0.00000000e00],
                    [-2.50000000e01, -4.14354109e00, 0.00000000e00],
                    [-3.52980820e00, -4.35829869e00, 3.55271368e-15],
                    [-3.63718700e00, -4.14354108e00, 0.00000000e00],
                ]
            ),
            "Fe": np.array(
                [
                    [0.00000000e00, 0.00000000e00, -2.50000000e01],
                    [-2.50000000e01, 0.00000000e00, -2.50000000e01],
                    [3.55271368e-15, 3.55271368e-15, -6.22930387e00],
                    [-2.50000000e01, 3.55271368e-15, -2.85707600e00],
                    [-2.01860032e00, 3.55271368e-15, -2.85707600e00],
                    [-1.43858289e00, 3.55271368e-15, -3.35213809e00],
                    [-1.71198739e00, 3.55271368e-15, -3.01038246e00],
                ]
            ),
            "Li": np.array(
                [
                    [3.55271368e-15, -2.50000000e01, -2.50000000e01],
                    [0.00000000e00, -2.50000000e01, -6.22930387e00],
                    [0.00000000e00, 0.00000000e00, -2.50000000e01],
                    [3.55271368e-15, 3.55271368e-15, -6.22930387e00],
                ]
            ),
        }

        for formula, domain in formal_domains.items():
            d = self.cpd_ternary_formal.domains[formula]
            d = d.round(6)  # to get rid of numerical errors from qhull
            assert d == approx(domain, abs=1e-5)
