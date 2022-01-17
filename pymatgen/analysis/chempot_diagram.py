# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module implements the construction and plotting of chemical potential diagrams
from a list of entries within a chemical system containing 2 or more elements. The
chemical potential diagram is the mathematical dual to the traditional compositional
phase diagram.

For more information, please cite/reference the paper below:

Todd, Paul K., McDermott, M.J., et al. “Selectivity in yttrium manganese oxide
synthesis via local chemical potentials in hyperdimensional phase space.”
ArXiv:2104.05986 [Cond-Mat], Aug. 2021. arXiv.org, http://arxiv.org/abs/2104.05986.

Please also consider referencing the original 1999 paper by H. Yokokawa,
who outlined many of its possible uses:

Yokokawa, H. "Generalized chemical potential diagram and its applications to
chemical reactions at interfaces between dissimilar materials." JPE 20,
258 (1999). https://doi.org/10.1361/105497199770335794
"""

import json
import os
from itertools import groupby
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import plotly.express as px
from monty.json import MSONable
from plotly.graph_objects import Figure, Mesh3d, Scatter, Scatter3d
from scipy.spatial import ConvexHull, HalfspaceIntersection

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core.composition import Composition, Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.coord import Simplex
from pymatgen.util.string import htmlify

with open(os.path.join(os.path.dirname(__file__), "..", "util", "plotly_chempot_layouts.json")) as f:
    plotly_layouts = json.load(f)


class ChemicalPotentialDiagram(MSONable):
    """
    The chemical potential diagram is the mathematical dual to the  compositional
    phase diagram. To create the diagram, convex minimization is
    performed in energy (E) vs. chemical potential (μ) space by taking the lower convex
    envelope of hyperplanes. Accordingly, "points" on the compositional phase diagram
    become N-dimensional convex polytopes (domains) in chemical potential space.

    For more information on this specific implementation of the algorithm,
    please cite/reference the paper below:

    Todd, Paul K., McDermott, M.J., et al. “Selectivity in yttrium manganese oxide
    synthesis via local chemical potentials in hyperdimensional phase space.”
    ArXiv:2104.05986 [Cond-Mat], Apr. 2021. arXiv.org, http://arxiv.org/abs/2104.05986
    """

    def __init__(
        self,
        entries: List[PDEntry],
        limits: Optional[Dict[Element, float]] = None,
        default_min_limit: Optional[float] = -20.0,
    ):
        """
        Args:
            entries: List of PDEntry-like objects containing a composition and
                energy. Must contain elemental references and be suitable for typical
                phase diagram construction. Entries must be within a chemical system
                of with 2+ elements
            limits: Bounds of elemental chemical potentials (min, max), which are
                used to construct the border hyperplanes used in the
                HalfSpaceIntersection algorithm; these constrain the space over which the
                domains are calculated and also determine the size of the plotted
                diagram. Any elemental limits not specified are covered in the
                default_min_limit argument
            default_min_limit (float): Default minimum chemical potential limit for
                unspecified elements within the "limits" argument. This results in
                default limits of (default_min_limit, 0)
        """
        self.entries = list(sorted(entries, key=lambda e: e.composition.reduced_composition))
        self.limits = limits
        self.default_min_limit = default_min_limit
        self.elements = list(sorted({els for e in self.entries for els in e.composition.elements}))
        self.dim = len(self.elements)
        self._min_entries, self._el_refs = self._get_min_entries_and_el_refs(self.entries)
        self._entry_dict = {e.composition.reduced_formula: e for e in self._min_entries}
        self._border_hyperplanes = self._get_border_hyperplanes()
        (
            self._hyperplanes,
            self._hyperplane_entries,
        ) = self._get_hyperplanes_and_entries()

        if self.dim < 2:
            raise ValueError("ChemicalPotentialDiagram currently requires phase diagrams with 2 or more elements!")

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
        formula_colors: List[str] = px.colors.qualitative.Dark2,
    ) -> Figure:
        """
        Plot the 2-dimensional or 3-dimensional chemical potential diagram using an
        interactive Plotly interface.

        Elemental axes can be specified; if none provided, will automatically default
        to first 2-3 elements within the "elements" attribute.

        In 3D, this method also allows for plotting of lower-dimensional "slices" of
        hyperdimensional polytopes (e.g., the LiMnO2 domain within a Y-Mn-O diagram).
        This allows for visualization of some of the phase boundaries that can only
        be seen fully in high dimensional space.

        Args:
            elements: list of elements to use as axes in the diagram. If None,
                automatically defaults to the first 2 or elements within the
                object's "elements" attribute.
            label_stable: whether or not to label stable phases by their reduced
                formulas. Defaults to True.
            formulas_to_draw: for 3-dimensional diagrams, an optional list of
                formulas to plot on the diagram; if these are from a different
                chemical system a 3-d polyhedron "slice" will be plotted.
            draw_formula_meshes: whether or not to draw a colored mesh for the
                optionally specified formulas_to_draw. Defaults to True.
            draw_formula_lines: whether or not to draw bounding lines for the
                optionally specified formulas_to_draw. Defaults to True.
            formula_colors: a list of colors to use in the plotting of the optionally
                specified formulas_to-draw. Defaults to a Plotly color scheme.

        Returns:
            A Plotly Figure object
        """
        if elements:
            elems = [Element(str(e)) for e in elements]
        else:
            elems = self.elements.copy()
            if self.dim > 3:
                elems = elems[:3]  # default to first three elements

        if len(elems) == 2 and self.dim == 2:
            fig = self._get_2d_plot(elements=elems, label_stable=label_stable)
        elif len(elems) == 2 and self.dim > 2:
            entries = [e for e in self.entries if set(e.composition.elements).issubset(elems)]
            cpd = ChemicalPotentialDiagram(
                entries=entries,
                limits=self.limits,
                default_min_limit=self.default_min_limit,
            )
            fig = cpd.get_plot(elements=elems, label_stable=label_stable)  # type: ignore
        else:
            fig = self._get_3d_plot(
                elements=elems,
                label_stable=label_stable,
                formulas_to_draw=formulas_to_draw,
                draw_formula_meshes=draw_formula_meshes,
                draw_formula_lines=draw_formula_lines,
                formula_colors=formula_colors,
            )

        return fig

    def _get_domains(self) -> Dict[str, np.ndarray]:
        """Returns a dictionary of domains as {formula: np.ndarray}"""
        hyperplanes = self._hyperplanes
        border_hyperplanes = self._border_hyperplanes
        entries = self._hyperplane_entries

        hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])
        interior_point = np.average(self.lims, axis=1).tolist()
        hs_int = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

        domains = {entry.composition.reduced_formula: [] for entry in entries}  # type: ignore

        for intersection, facet in zip(hs_int.intersections, hs_int.dual_facets):
            for v in facet:
                if v < len(entries):
                    this_entry = entries[v]
                    formula = this_entry.composition.reduced_formula
                    domains[formula].append(intersection)

        return {k: np.array(v) for k, v in domains.items() if v}

    def _get_border_hyperplanes(self) -> np.ndarray:
        """Returns an array of the bounding hyperplanes given by elemental limits"""
        border_hyperplanes = np.array([[0] * (self.dim + 1)] * (2 * self.dim))

        for idx, limit in enumerate(self.lims):
            border_hyperplanes[2 * idx, idx] = -1
            border_hyperplanes[2 * idx, -1] = limit[0]
            border_hyperplanes[(2 * idx) + 1, idx] = 1
            border_hyperplanes[(2 * idx) + 1, -1] = limit[1]

        return border_hyperplanes

    def _get_hyperplanes_and_entries(self) -> Tuple[np.ndarray, List[PDEntry]]:
        """Returns both the array of hyperplanes, as well as a list of the minimum
        entries"""
        data = np.array(
            [
                [e.composition.get_atomic_fraction(el) for el in self.elements] + [e.energy_per_atom]
                for e in self._min_entries
            ]
        )
        vec = [self.el_refs[el].energy_per_atom for el in self.elements] + [-1]
        form_e = -np.dot(data, vec)

        inds = np.where(form_e < -PhaseDiagram.formation_energy_tol)[0].tolist()

        inds.extend([self._min_entries.index(el) for el in self.el_refs.values()])

        hyperplanes = data[inds]
        hyperplanes[:, -1] = hyperplanes[:, -1] * -1
        hyperplane_entries = [self._min_entries[i] for i in inds]

        return hyperplanes, hyperplane_entries

    def _get_2d_plot(self, elements: List[Element], label_stable: Optional[bool]) -> Figure:
        """Returns a Plotly figure for a 2-dimensional chemical potential diagram"""
        domains = self.domains.copy()
        elem_indices = [self.elements.index(e) for e in elements]

        annotations = []
        for formula, pts in domains.items():
            formula_elems = set(Composition(formula).elements)
            if not formula_elems.issubset(elements):
                continue
            pts_2d = np.array(pts[:, elem_indices])
            entry = self.entry_dict[formula]
            ann_formula = formula
            if hasattr(entry, "original_entry"):
                ann_formula = entry.original_entry.composition.reduced_formula

            center = pts_2d.mean(axis=0)
            normal = get_2d_orthonormal_vector(pts_2d)
            ann_loc = center + 0.8 * normal
            annotation = self._get_annotation(ann_loc, ann_formula)
            annotations.append(annotation)

        layout = plotly_layouts["default_layout_2d"].copy()
        layout.update(self._get_axis_layout_dict(elements))
        if label_stable:
            layout.update({"annotations": annotations})

        data = self._get_2d_domain_lines(domains, elements)

        fig = Figure(data, layout)

        return fig

    def _get_3d_plot(
        self,
        elements: List[Element],
        label_stable: Optional[bool],
        formulas_to_draw: Optional[List[str]],
        draw_formula_meshes: Optional[bool],
        draw_formula_lines: Optional[bool],
        formula_colors: Optional[List[str]] = px.colors.qualitative.Dark2,
    ) -> Figure:
        """Returns a Plotly figure for a 3-dimensional chemical potential diagram"""
        if not formulas_to_draw:
            formulas_to_draw = []

        elem_indices = [self.elements.index(e) for e in elements]

        domains = self.domains.copy()
        draw_domains = {}
        draw_comps = [Composition(formula).reduced_composition for formula in formulas_to_draw]
        annotations = []

        for formula, pts in domains.items():
            entry = self.entry_dict[formula]
            pts_3d = np.array(pts[:, elem_indices])
            contains_target_elems = set(entry.composition.elements).issubset(elements)

            if formulas_to_draw:
                if entry.composition.reduced_composition in draw_comps:
                    domains[formula] = None
                    draw_domains[formula] = pts_3d

                    if contains_target_elems:
                        domains[formula] = pts_3d
                    else:
                        continue

            if not contains_target_elems:
                domains[formula] = None
                continue

            simplexes, ann_loc = self._get_3d_domain_simplexes_and_ann_loc(pts_3d)

            ann_formula = formula
            if hasattr(entry, "original_entry"):
                ann_formula = entry.original_entry.composition.reduced_formula

            annotation = self._get_annotation(ann_loc, ann_formula)
            annotations.append(annotation)

            domains[formula] = simplexes

        layout = plotly_layouts["default_layout_3d"].copy()
        layout["scene"].update(self._get_axis_layout_dict(elements))
        layout["scene"]["annotations"] = None

        if label_stable:
            layout["scene"].update({"annotations": annotations})
        layout["scene_camera"] = dict(eye=dict(x=0, y=0, z=2.0), projection=dict(type="orthographic"))

        data = self._get_3d_domain_lines(domains)

        if draw_formula_lines:
            data.extend(self._get_3d_formula_lines(draw_domains, formula_colors))
        if draw_formula_meshes:
            data.extend(self._get_3d_formula_meshes(draw_domains, formula_colors))

        fig = Figure(data, layout)

        return fig

    @staticmethod
    def _get_2d_domain_lines(domains: Dict[str, np.ndarray], elements: List[Element]) -> List[Scatter]:
        """Returns a list of Scatter objects tracing the domain lines on a
        2-dimensional chemical potential diagram"""
        x, y = [], []
        elems = set(elements)
        for formula, pts in domains.items():
            formula_elems = set(Composition(formula).elements)
            if not formula_elems.issubset(elems):
                continue

            x.extend(pts[:, 0].tolist() + [None])
            y.extend(pts[:, 1].tolist() + [None])

        lines = [
            Scatter(
                x=x,
                y=y,
                mode="lines+markers",
                line=dict(color="black", width=3),
                showlegend=False,
            )
        ]
        return lines

    @staticmethod
    def _get_3d_domain_lines(domains: Dict[str, List[Simplex]]) -> List[Scatter3d]:
        """Returns a list of Scatter3d objects tracing the domain lines on a
        3-dimensional chemical potential diagram"""
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

    @staticmethod
    def _get_3d_domain_simplexes_and_ann_loc(
        points_3d: np.ndarray,
    ) -> Tuple[List[Simplex], np.ndarray]:
        """Returns a list of Simplex objects and coordinates of annotation for one
        domain in a 3-d chemical potential diagram. Uses PCA to project domain
        into 2-dimensional space so that ConvexHull can be used to identify the
        bounding polygon"""
        points_2d, v, w = simple_pca(points_3d, k=2)
        domain = ConvexHull(points_2d)
        centroid_2d = get_centroid_2d(points_2d[domain.vertices])
        ann_loc = centroid_2d @ w.T + np.mean(points_3d.T, axis=1)

        simplexes = [Simplex(points_3d[indices]) for indices in domain.simplices]

        return simplexes, ann_loc

    @staticmethod
    def _get_3d_formula_meshes(
        draw_domains: Dict[str, np.ndarray],
        formula_colors: Optional[List[str]],
    ) -> List[Mesh3d]:
        """Returns a list of Mesh3d objects for the domains specified by the
        user (i.e., draw_domains)"""
        meshes = []
        if formula_colors is None:
            formula_colors = px.colors.qualitative.Dark2

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

    @staticmethod
    def _get_3d_formula_lines(
        draw_domains: Dict[str, np.ndarray],
        formula_colors: Optional[List[str]],
    ) -> List[Scatter3d]:
        """Returns a list of Scatter3d objects defining the bounding polyhedra"""
        if formula_colors is None:
            formula_colors = px.colors.qualitative.Dark2

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
    def _get_min_entries_and_el_refs(
        entries: List[PDEntry],
    ) -> Tuple[List[PDEntry], Dict[Element, PDEntry]]:
        """Returns a list of the minimum-energy entries at each composition and the
        entries corresponding to the elemental references"""
        el_refs = {}
        min_entries = []

        for c, g in groupby(entries, key=lambda e: e.composition.reduced_formula):
            c = Composition(c)
            group = list(g)
            min_entry = min(group, key=lambda e: e.energy_per_atom)
            if c.is_element:
                el_refs[c.elements[0]] = min_entry
            min_entries.append(min_entry)

        return min_entries, el_refs

    @staticmethod
    def _get_annotation(ann_loc: np.ndarray, formula: str) -> Dict[str, Union[str, float]]:
        """Returns a Plotly annotation dict given a formula and location"""
        formula = htmlify(formula)
        annotation = plotly_layouts["default_annotation_layout"].copy()
        annotation.update({"x": ann_loc[0], "y": ann_loc[1], "text": formula})
        if len(ann_loc) == 3:
            annotation.update({"z": ann_loc[2]})
        return annotation

    @staticmethod
    def _get_axis_layout_dict(elements: List[Element]) -> Dict[str, str]:
        """Returns a Plotly layout dict for either 2-d or 3-d axes"""
        axes = ["xaxis", "yaxis"]
        layout_name = "default_2d_axis_layout"

        if len(elements) == 3:
            axes.append("zaxis")
            layout_name = "default_3d_axis_layout"

        def get_chempot_axis_title(element) -> str:
            return f"μ<sub>{element}</sub> - μ<sub>{element}</sub><sup>o</sup> (eV)"

        axes_layout = {}
        for ax, el in zip(axes, elements):
            layout = plotly_layouts[layout_name].copy()
            layout["title"] = get_chempot_axis_title(el)
            axes_layout[ax] = layout

        return axes_layout

    @property
    def domains(self) -> Dict[str, np.ndarray]:
        """
        Mapping of formulas to array of domain boundary points
        """
        return self._get_domains()

    @property
    def lims(self) -> np.ndarray:
        """Returns array of limits used in constructing hyperplanes"""
        lims = np.array([[self.default_min_limit, 0]] * self.dim)
        for idx, elem in enumerate(self.elements):
            if self.limits and elem in self.limits:
                lims[idx, :] = self.limits[elem]
        return lims

    @property
    def entry_dict(self) -> Dict[str, ComputedEntry]:
        """Mapping between reduced formula and ComputedEntry"""
        return self._entry_dict

    @property
    def hyperplanes(self) -> np.ndarray:
        """Returns array of hyperplane data"""
        return self._hyperplanes

    @property
    def hyperplane_entries(self) -> List[PDEntry]:
        """Returns list of entries corresponding to hyperplanes"""
        return self._hyperplane_entries

    @property
    def border_hyperplanes(self) -> np.ndarray:
        """Returns bordering hyperplanes"""
        return self._border_hyperplanes

    @property
    def el_refs(self) -> Dict[Element, PDEntry]:
        """Returns a dictionary of elements and reference entries"""
        return self._el_refs

    @property
    def chemical_system(self) -> str:
        """Returns the chemical system (A-B-C-...) of diagram object"""
        return "-".join(sorted(e.symbol for e in self.elements))

    def __repr__(self):
        return f"ChemicalPotentialDiagram for {self.chemical_system} with {len(self.entries)} entries"


def simple_pca(data: np.ndarray, k: int = 2) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    A barebones implementation of principal component analysis (PCA) used in the
    ChemicalPotentialDiagram class for plotting.

    Args:
        data: array of observations
        k: Number of principal components returned

    Returns:
        Tuple of projected data, eigenvalues, eigenvectors
    """
    data = data - np.mean(data.T, axis=1)  # centering the data
    cov = np.cov(data.T)  # calculating covariance matrix
    v, w = np.linalg.eig(cov)  # performing eigendecomposition
    idx = v.argsort()[::-1]  # sorting the components
    v = v[idx]
    w = w[:, idx]
    scores = data.dot(w[:, :k])

    return scores, v[:k], w[:, :k]


def get_centroid_2d(vertices: np.ndarray) -> np.ndarray:
    """
    A barebones implementation of the formula for calculating the centroid of a 2D
    polygon. Useful for calculating the location of an annotation on a chemical
    potential domain within a 3D chemical potential diagram.

    **NOTE**: vertices must be ordered circumfrentially!

    Args:
        vertices: array of 2-d coordinates corresponding to a polygon, ordered
            circumfrentially

    Returns:
        Array giving 2-d centroid coordinates
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
    centroid = np.array([prefactor * cx, prefactor * cy])

    return centroid


def get_2d_orthonormal_vector(line_pts: np.ndarray) -> np.ndarray:
    """
    Calculates a vector that is orthonormal to a line given by a set of points. Used
    for determining the location of an annotation on a 2-d chemical potential diagram.

    Args:
        line_pts: a 2x2 array in the form of [[x0, y0], [x1, y1]] giving the
            coordinates of a line

    Returns:

    """
    x = line_pts[:, 0]
    y = line_pts[:, 1]

    x_diff = abs(x[1] - x[0])
    y_diff = abs(y[1] - y[0])

    if np.isclose(x_diff, 0):
        theta = np.pi / 2
    else:
        theta = np.arctan(y_diff / x_diff)

    vec = np.array([np.sin(theta), np.cos(theta)])

    return vec
