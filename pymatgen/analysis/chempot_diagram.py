# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module implements the construction and plotting of chemical potential diagrams
from a list of entries within a chemical system containing 3 or more elements. The
chemical potential diagram is the mathematical dual to the traditional compositional
phase diagram.

For more information, please reference the paper below:

Todd, Paul K., McDermott, M.J., et al. “Selectivity in yttrium manganese oxide
synthesis via local chemical potentials in hyperdimensional phase space.”
ArXiv:2104.05986 [Cond-Mat], Apr. 2021. arXiv.org, http://arxiv.org/abs/2104.05986.
"""

import json
import os
from itertools import groupby
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import plotly.express as px
from monty.json import MSONable
from plotly.graph_objects import Scatter3d, Mesh3d, Figure
from scipy.spatial import ConvexHull, HalfspaceIntersection
from scipy.spatial.qhull import QhullError

from pymatgen.analysis.phase_diagram import PDPlotter, PhaseDiagram, PDEntry
from pymatgen.core.composition import Composition, Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.coord import Simplex

with open(os.path.join(os.path.dirname(__file__), "..", "util", "plotly_chempot_layouts.json")) as f:
    plotly_layouts = json.load(f)


class ChemicalPotentialDiagram(MSONable):
    """
    The chemical potential diagram is the mathematical dual to the
    traditional compositional phase diagram. To create the diagram, convex
    minimization is performed in E vs. μ space by taking the lower convex envelope of
    hyperplanes. Accordingly, "points" on the compositional phase diagram
    become N-dimensional convex polytopes (domains) in chemical potential space.

    For more information, please reference the paper below:

    Todd, Paul K., McDermott, M.J., et al. “Selectivity in yttrium manganese oxide
    synthesis via local chemical potentials in hyperdimensional Phase Space.”
    ArXiv:2104.05986 [Cond-Mat], Apr. 2021. arXiv.org, http://arxiv.org/abs/2104.05986.
    """

    def __init__(
        self,
        entries: List[PDEntry],
        limits: Optional[Dict[Element, float]] = None,
        default_limit: Optional[float] = -20.0,
    ):
        """
        Args:
            entries: List of PDEntry-like objects containing a composition and energy
            limits: Chemical potential value chosen for bordering elemental
                hyperplanes; these constrain the space over which the domains are
                calculated
            default_limit (float): Default minimum chemical potential limit for
                unspecified elemental limits
        """
        self.entries = sorted(entries, key=lambda e: e.composition.reduced_composition)
        self.limits = limits
        self.default_limit = default_limit

        self.elements = list(sorted({els for e in self.entries for els in
                                     e.composition.elements}))
        self.dim = len(self.elements)
        self._min_entries, self.el_refs = self._get_min_entries_and_el_refs(entries)
        self._entry_dict = {e.composition.reduced_formula: e for e in self._min_entries}
        self._border_hyperplanes = self._get_border_hyperplanes()
        self._hyperplanes, self._hyperplane_entries = self._get_hyperplanes_and_entries()

        if self.dim < 3:
            raise ValueError("ChemicalPotentialDiagram currently requires phase diagrams "
                             "with 3 or more elements!")

        if len(self.el_refs) != self.dim:
            missing = set(self.elements).difference(self.el_refs.keys())
            raise ValueError(f"There are no entries for the terminal elements: {missing}")

    def get_plot(
        self,
        elements: Optional[List[Union[Element, str]]] = None,
        label_stable: Optional[bool] = True,
        formulas_to_draw: Optional[List[str]] = None,
        draw_formula_meshes: Optional[bool] = True,
        draw_formula_lines: Optional[bool] = True,
        formula_colors: Optional[List[str]] = px.colors.qualitative.Dark2,
    ) -> Scatter3d:
        """

        Args:
            elements:
            label_stable:
            formulas_to_draw:
            draw_formula_meshes:
            draw_formula_lines:
            formula_colors:

        Returns:

        """
        if not elements:
            elements = self.elements[:3]  # type: ignore
        if not formulas_to_draw:
            formulas_to_draw = []

        elems = [Element(e) for e in elements]  # type: ignore
        elem_indices = [self.elements.index(e) for e in elems]

        domains = self.domains.copy()
        draw_domains = {}
        draw_comps = [Composition(formula).reduced_composition for formula in
                      formulas_to_draw]
        annotations = []

        for formula, points in domains.items():
            entry = self.entry_dict[formula]
            points_3d = np.array(points[:, elem_indices])
            contains_target_elems = set(entry.composition.elements).issubset(elems)

            if formulas_to_draw:
                if entry.composition.reduced_composition in draw_comps:
                    domains[formula] = None
                    draw_domains[formula] = points_3d

                    if contains_target_elems:
                        domains[formula] = points_3d
                    else:
                        continue

            if not contains_target_elems:
                domains[formula] = None
                continue

            simplices, ann_loc = self._get_domain_simplices_and_ann_loc(points_3d)

            ann_formula = formula
            if hasattr(entry, "original_entry"):
                ann_formula = entry.original_entry.composition.reduced_formula
            annotation = self._get_annotation(ann_loc, ann_formula)
            annotations.append(annotation)

            domains[formula] = simplices

        layout = plotly_layouts["default_layout_3d"].copy()
        layout["scene"].update(self._get_axis_layout_dict(elements))
        layout["scene"]["annotations"] = None

        if label_stable:
            layout["scene"].update({"annotations": annotations})
        layout["scene_camera"] = dict(eye=dict(x=0, y=0, z=2.0), projection=dict(type="orthographic"))

        data = self._get_domain_lines(domains)

        if draw_formula_lines:
            data.extend(self._get_formula_lines(draw_domains, formula_colors))
        if draw_formula_meshes:
            data.extend(self._get_formula_meshes(draw_domains, formula_colors))

        figure = Figure(data, layout)

        return figure

    def _get_domains(self):
        hyperplanes = self._hyperplanes
        border_hyperplanes = self._border_hyperplanes
        entries = self._hyperplane_entries

        hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])
        interior_point = np.average(self.lims, axis=1).tolist()
        hs_int = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

        domains = {entry.composition.reduced_formula: [] for entry in entries}

        for intersection, facet in zip(hs_int.intersections, hs_int.dual_facets):
            for v in facet:
                if v < len(entries):
                    this_entry = entries[v]
                    formula = this_entry.composition.reduced_formula
                    domains[formula].append(intersection)

        return {k: np.array(v) for k, v in domains.items() if v}

    def _get_border_hyperplanes(self):
        border_hyperplanes = np.array(([[0] * (self.dim + 1)] * (2 * self.dim)))

        for idx, limit in enumerate(self.lims):
            border_hyperplanes[2 * idx, idx] = -1
            border_hyperplanes[2 * idx, -1] = limit[0]
            border_hyperplanes[(2 * idx) + 1, idx] = 1
            border_hyperplanes[(2 * idx) + 1, -1] = limit[1]

        return border_hyperplanes

    def _get_hyperplanes_and_entries(self):
        data = np.array(
            [[e.composition.get_atomic_fraction(el) for el in self.elements] + [
                e.energy_per_atom] for e in self._min_entries]
        )
        vec = [self.el_refs[el].energy_per_atom for el in self.elements] + [-1]
        form_e = -np.dot(data, vec)
        inds = np.where(form_e < -PhaseDiagram.formation_energy_tol)[0].tolist()
        inds.extend([self._min_entries.index(el) for el in self.el_refs.values()])

        hyperplanes = data[inds]
        hyperplanes[:, -1] = hyperplanes[:, -1] * -1
        hyperplane_entries = [self._min_entries[i] for i in inds]

        return hyperplanes, hyperplane_entries

    @staticmethod
    def _get_min_entries_and_el_refs(entries):
        el_refs = {}
        min_entries = []

        for c, g in groupby(entries, key=lambda e: e.composition.reduced_composition):
            g = list(g)
            min_entry = min(g, key=lambda e: e.energy_per_atom)
            if c.is_element:
                el_refs[c.elements[0]] = min_entry
            min_entries.append(min_entry)

        return min_entries, el_refs

    @staticmethod
    def _get_domain_simplices_and_ann_loc(points_3d):
        try:
            domain = ConvexHull(points_3d)
            ann_loc = np.mean(points_3d.T, axis=1)
        except QhullError:
            points_2d, v, w = simple_pca(points_3d, k=2)
            domain = ConvexHull(points_2d)
            centroid_2d = get_centroid_2d(points_2d[domain.vertices])
            ann_loc = centroid_2d @ w.T + np.mean(points_3d.T,
                                                  axis=1)  # recover orig 3D coords from eigenvectors

        simplices = [Simplex(points_3d[indices]) for indices in domain.simplices]

        return simplices, ann_loc

    def _get_formula_meshes(self, draw_domains, formula_colors) -> List[Mesh3d]:
        meshes = []
        for idx, (formula, coords) in enumerate(draw_domains.items()):
            points_3d = coords[:, :3]
            mesh = Mesh3d(
                    x=points_3d[:, 0],
                    y=points_3d[:, 1],
                    z=points_3d[:, 2],
                    alphahull=0,
                    showlegend=True,
                    lighting=dict(fresnel=1.0),
                    color=formula_colors[idx],
                    name=f"{formula} (mesh)",
                    opacity=0.13,
                )
            meshes.append(mesh)
        return meshes

    def _get_formula_lines(self, draw_domains, formula_colors) -> List[Scatter3d]:
        lines = []
        for idx, (formula, coords) in enumerate(draw_domains.items()):
            points_3d = coords[:, :3]
            domain = ConvexHull(points_3d[:, :-1])

            simplexes = [Simplex(points_3d[indices]) for indices in domain.simplices]
            x, y, z = [], [], []
            for s in simplexes:
                x.extend(s.coords[:, 0].tolist() + [None])
                y.extend(s.coords[:, 1].tolist() + [None])
                z.extend(s.coords[:, 2].tolist() + [None])

            line = Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode="lines",
                    line={"width": 8, "color": formula_colors[idx]},
                    opacity=1.0,
                    name=f"{formula} (lines)",
                )
            lines.append(line)
        return lines

    @staticmethod
    def _get_domain_lines(domains) -> List[Scatter3d]:
        x, y, z = [], [], []
        for phase, simplexes in domains.items():
            if simplexes:
                for s in simplexes:
                    x.extend(s.coords[:, 0].tolist() + [None])
                    y.extend(s.coords[:, 1].tolist() + [None])
                    z.extend(s.coords[:, 2].tolist() + [None])

        lines = [
            Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                line=dict(color="black", width=4.5),
                showlegend=False,
            )
        ]
        return lines

    @property
    def domains(self) -> Dict[str, np.array]:
        """
        Mapping of formulas to array of domain boundary points
        """
        return self._get_domains()

    @property
    def lims(self):
        """ Returns array of limits used in constructing hyperplanes"""
        lims = np.array([[self.default_limit, 0]] * self.dim)
        for idx, elem in enumerate(self.elements):
            if self.limits and elem in self.limits:
                lims[idx, :] = self.limits[elem]
        return lims

    @property
    def entry_dict(self) -> Dict[str, ComputedEntry]:
        """Mapping between reduced formula and ComputedEntry"""
        return self._entry_dict

    @property
    def hyperplanes(self) -> np.array:
        """Returns array of hyperplane data"""
        return self._hyperplanes

    @property
    def hyperplane_entries(self) -> List[PDEntry]:
        """Returns list of entries corresponding to hyperplanes"""
        return self._hyperplane_entries

    @property
    def border_hyperplanes(self) -> np.array:
        """Returns bordering hyperplanes"""
        return self._border_hyperplanes

    @property
    def chemical_system(self) -> str:
        """Returns the chemical system (A-B-C-...) of diagram object"""
        return "-".join(sorted([e.symbol for e in self.elements]))

    @staticmethod
    def _get_chempot_axis_title(element) -> str:
        return f"μ<sub>{str(element)}</sub> - μ<sub>" \
               f"{str(element)}</sub><sup>o</sup> (eV)"

    @staticmethod
    def _get_annotation(ann_loc, formula) -> Dict[str, Union[str, float]]:
        formula = PDPlotter._htmlize_formula(formula)
        annotation = plotly_layouts["default_annotation_layout"].copy()
        annotation.update(
            {
                "x": ann_loc[0],
                "y": ann_loc[1],
                "z": ann_loc[2],
                "text": formula,
            }
        )
        return annotation

    def _get_axis_layout_dict(self, elements):
        axes = ["xaxis", "yaxis", "zaxis"]
        axes_layout = {}
        for ax, el in zip(axes, elements):
            layout = plotly_layouts["default_axis_layout"].copy()
            layout["title"] = self._get_chempot_axis_title(el)
            axes_layout[ax] = layout

        return axes_layout

    def __repr__(self):
        return f"ChemicalPotentialDiagram for {self.chemical_system} with {len(self.entries)} " \
               f"entries"


def simple_pca(data: np.array, k: int = 2) -> Tuple[np.array, np.array, np.array]:
    """
    A barebones implementation of principal component analysis (PCA) used in the
    ChemicalPotentialDiagram class for plotting.

    Args:
        data: array of observations
        k: Number of principal components returned

    Returns:
        tuple: Projected data, eigenvalues, eigenvectors
    """
    data = data - np.mean(data.T, axis=1)  # centering the data
    cov = np.cov(data.T)  # calculating covariance matrix
    v, w = np.linalg.eig(cov)  # performing eigendecomposition
    idx = v.argsort()[::-1]  # sorting the components
    v = v[idx]
    w = w[:, idx]
    scores = data.dot(w[:, :k])

    return scores, v[:k], w[:, :k]


def get_centroid_2d(vertices: np.array) -> np.array:
    """
    A barebones implementation of the formula for calculating the centroid of a 2D
    polygon.

    **NOTE**: vertices must be ordered circumfrentially!

    Args:
        vertices:

    Returns:

    """
    n = len(vertices)
    cx = 0
    cy = 0
    a = 0

    for i in range(0, n - 1):
        xi = vertices[i, 0]
        yi = vertices[i, 1]
        xi_p = vertices[i + 1, 0]
        yi_p = vertices[i + 1, 1]
        common_term = xi * yi_p - xi_p * yi

        cx += (xi + xi_p) * common_term
        cy += (yi + yi_p) * common_term
        a += common_term

    prefactor = 0.5 / (6 * a)
    return np.array([prefactor * cx, prefactor * cy])
