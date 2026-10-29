"""This module contains classes to wrap Python VTK to make nice molecular plots.

The default styling follows VESTA's ball-and-stick model: spheres drawn at 40% of
the ionic radius of each species, bonds split into two halves colored by their
end atoms, and translucent coordination polyhedra colored by the central atom
with their edges outlined.
"""

from __future__ import annotations

import itertools
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from monty.dev import requires
from monty.serialization import loadfn
from scipy.spatial import ConvexHull, QhullError

from pymatgen.core import Element, PeriodicSite, Species, Structure
from pymatgen.util.coord import in_coord_list

try:
    # Importing vtkRenderingOpenGL2, vtkRenderingFreeType and vtkInteractionStyle
    # registers the OpenGL, text and interactor-style factory overrides.
    import vtkmodules.vtkRenderingOpenGL2  # noqa: F401
    from vtkmodules.vtkCommonCore import vtkPoints, vtkUnsignedCharArray
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.vtkFiltersCore import vtkFeatureEdges, vtkTubeFilter
    from vtkmodules.vtkFiltersHybrid import vtkRenderLargeImage
    from vtkmodules.vtkFiltersSources import vtkLineSource, vtkSphereSource
    from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera as TrackballCamera
    from vtkmodules.vtkIOImage import vtkBMPWriter, vtkJPEGWriter, vtkPNGWriter, vtkTIFFWriter
    from vtkmodules.vtkRenderingCore import (
        vtkActor,
        vtkActor2D,
        vtkCellPicker,
        vtkFollower,
        vtkPolyDataMapper,
        vtkRenderer,
        vtkRenderWindow,
        vtkRenderWindowInteractor,
        vtkTextMapper,
    )
    from vtkmodules.vtkRenderingFreeType import vtkVectorText

    vtk_available = True
except ImportError:
    # VTK not present. TrackballCamera is set to object so the interactor classes can still be defined.
    vtk_available = False
    TrackballCamera = object
    vtkPoints = vtkUnsignedCharArray = vtkCellArray = vtkPolyData = None
    vtkFeatureEdges = vtkTubeFilter = None
    vtkRenderLargeImage = vtkLineSource = vtkSphereSource = vtkVectorText = None
    vtkBMPWriter = vtkJPEGWriter = vtkPNGWriter = vtkTIFFWriter = None
    vtkActor = vtkActor2D = vtkCellPicker = vtkFollower = vtkPolyDataMapper = None
    vtkRenderer = vtkRenderWindow = vtkRenderWindowInteractor = vtkTextMapper = None

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, ClassVar

    from numpy.typing import ArrayLike

    from pymatgen.core import Molecule, Site

MODULE_DIR = Path(__file__).resolve().parent
EL_COLORS = loadfn(MODULE_DIR / "ElementColorSchemes.yaml")

# VESTA renders atoms in its ball-and-stick model at 40% of their radii.
ATOM_RADIUS_SCALE = 0.4
# Radius (Å) for species with no tabulated ionic or atomic radius.
DEFAULT_IONIC_RADIUS = 1.0
# Fractional-coordinate tolerance for treating a site as lying on a cell boundary.
BOUNDARY_TOL = 1e-3


def get_ionic_radius(specie: Species | Element, anion_symbols: set[str] | None = None) -> float:
    """Ionic radius (Å) of a species, used to size its sphere.

    Species with an oxidation state use their Shannon ionic radius. Bare elements
    use their average anionic radius if they are listed in anion_symbols and their
    average cationic radius otherwise. Both fall back to the average ionic radius,
    then the atomic radius, then DEFAULT_IONIC_RADIUS.

    Args:
        specie: Element or Species.
        anion_symbols: Symbols of elements to treat as anions when no oxidation
            state is set. Defaults to none.

    Returns:
        float: Ionic radius in Å.
    """
    try:
        if isinstance(specie, Species) and specie.oxi_state and specie.ionic_radius:
            return float(specie.ionic_radius)
        element = Element(specie.symbol)
    except (AttributeError, ValueError):
        # DummySpecies and other non-element species have no tabulated radii.
        return DEFAULT_IONIC_RADIUS
    is_anion = specie.symbol in (anion_symbols or set())
    primary = element.average_anionic_radius if is_anion else element.average_cationic_radius
    for radius in (primary, element.average_ionic_radius, element.atomic_radius):
        if radius:
            return float(radius)
    return DEFAULT_IONIC_RADIUS


def _get_anion_symbols(structure: Structure | Molecule) -> set[str]:
    """Symbol of the most electronegative element, treated as the anion for sphere sizes."""
    elements = {Element(el.symbol) for el in structure.elements}
    if len(elements) < 2:
        return set()
    return {max(elements, key=lambda el: el.X).symbol}


def _make_polydata(polygons: Sequence[Sequence[ArrayLike]]) -> vtkPolyData:
    """Polydata with one polygon cell per sequence of Cartesian vertex coordinates.

    Coincident vertices are merged into one point so that adjacent polygons share edges.
    """
    points = vtkPoints()
    polys = vtkCellArray()
    point_ids: dict[tuple[float, ...], int] = {}
    for polygon in polygons:
        polys.InsertNextCell(len(polygon))
        for coord in polygon:
            coord = np.asarray(coord, dtype=float)
            key = tuple(np.round(coord, 6))
            if key not in point_ids:
                point_ids[key] = points.InsertNextPoint(coord)
            polys.InsertCellPoint(point_ids[key])
    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(polys)
    return polydata


def _get_boundary_images(site: PeriodicSite) -> list[PeriodicSite]:
    """Periodic images of a site lying on a cell face, edge or corner that also lie on the cell boundary.

    For example, a site at the origin has seven further images at the other corners.
    """
    shifts = []
    for frac_coord in site.frac_coords:
        if abs(frac_coord) < BOUNDARY_TOL:
            shifts.append((0, 1))
        elif abs(frac_coord - 1) < BOUNDARY_TOL:
            shifts.append((0, -1))
        else:
            shifts.append((0,))
    return [
        PeriodicSite(
            site.species,
            site.frac_coords + np.array(shift),
            site.lattice,
            to_unit_cell=False,
            properties=site.properties,
        )
        for shift in itertools.product(*shifts)
        if any(shift)
    ]


class StructureVis:
    """Structure visualization using VTK."""

    @requires(vtk_available, "Visualization requires the installation of VTK with Python bindings.")
    def __init__(
        self,
        element_color_mapping: dict[str, Sequence[int]] | None = None,
        show_unit_cell: bool = True,
        show_bonds: bool = False,
        show_polyhedron: bool = True,
        poly_radii_tol_factor: float = 0.5,
        excluded_bonding_elements: Sequence[str] | None = None,
        offscreen: bool = False,
        atom_radius_scale: float = ATOM_RADIUS_SCALE,
        bond_radius: float = 0.1,
        polyhedron_opacity: float = 0.5,
    ) -> None:
        """
        Constructs a Structure Visualization.

        Args:
            element_color_mapping: Optional color mapping for the elements,
                as a dict of {symbol: rgb tuple}. For example, {"Fe": (255,
                123,0), ....} If None is specified, VESTA's color scheme is used.
            show_unit_cell: Set to False to not show the unit cell
                boundaries. Defaults to True.
            show_bonds: Set to True to show bonds. Defaults to False.
            show_polyhedron: Set to True to show polyhedrons. Defaults to
                True.
            poly_radii_tol_factor: The polyhedron and bonding code uses the
                ionic radii of the elements or species to determine if two
                atoms are bonded. This specifies a tolerance scaling factor
                such that atoms which are (1 + poly_radii_tol_factor) * sum
                of ionic radii apart are still considered as bonded.
            excluded_bonding_elements: List of atom types to exclude from
                bonding determination. Defaults to an empty list. Useful
                when trying to visualize a certain atom type in the
                framework (e.g., Li in a Li-ion battery cathode material).
            offscreen: Render without opening a window, e.g., to write images
                in scripts or on headless machines. show() is unavailable when
                True. Defaults to False.
            atom_radius_scale: Spheres are drawn with radius atom_radius_scale *
                ionic radius. Defaults to 0.4, as in VESTA's ball-and-stick model.
            bond_radius: Radius of bond cylinders in Å. Defaults to 0.1.
            polyhedron_opacity: Opacity of coordination polyhedra. Defaults to 0.5.

        Useful keyboard shortcuts implemented.
            h : Show help
            A/a : Increase/decrease cell by one unit vector in a-direction
            B/b : Increase/decrease cell by one unit vector in b-direction
            C/c : Increase/decrease cell by one unit vector in c-direction
            # : Toggle showing of polyhedrons
            - : Toggle showing of bonds
            [ : Decrease poly_radii_tol_factor by 0.05
            ] : Increase poly_radii_tol_factor by 0.05
            r : Reset camera direction
            o : Orthogonalize structure
            Up/Down : Rotate view along Up direction by 90 clock/anticlockwise
            Left/right : Rotate view along camera direction by 90
            clock/anticlockwise
        """
        # create a rendering window and renderer
        self.offscreen = offscreen
        self.ren = vtkRenderer()
        self.ren.SetBackground(1, 1, 1)
        self.ren.UseFXAAOn()
        self.ren_win = vtkRenderWindow()
        self.ren_win.SetOffScreenRendering(offscreen)
        self.ren_win.SetSize(800, 800)
        self.ren_win.AddRenderer(self.ren)
        self.title = "Structure Visualizer"
        self.iren = vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.ren_win)
        self.mapper_map: dict = {}
        self.structure: Structure | Molecule | None = None
        self._anion_symbols: set[str] = set()

        self.el_color_mapping = element_color_mapping or EL_COLORS["VESTA"]
        self.show_unit_cell = show_unit_cell
        self.show_bonds = show_bonds
        self.show_polyhedron = show_polyhedron
        self.poly_radii_tol_factor = poly_radii_tol_factor
        self.excluded_bonding_elements = list(excluded_bonding_elements or [])
        self.atom_radius_scale = atom_radius_scale
        self.bond_radius = bond_radius
        self.polyhedron_opacity = polyhedron_opacity
        self.show_help = True
        self.supercell = np.eye(3)
        self.redraw()

        # Keep a reference: re-wrapping via iren.GetInteractorStyle() loses the Python attributes.
        self.interactor_style = StructureInteractorStyle(self)
        self.iren.SetInteractorStyle(self.interactor_style)
        self.ren.parent = self

    def rotate_view(self, axis_ind: int = 0, angle: float = 0) -> None:
        """
        Rotate the camera view.

        Args:
            axis_ind: 0 rolls the camera about the view direction, 1 rotates it
                about the view-up vector (azimuth) and 2 about the horizontal
                screen axis (pitch). Defaults to 0.
            angle: Angle to rotate by in degrees. Defaults to 0.
        """
        camera = self.ren.GetActiveCamera()
        if axis_ind == 0:
            camera.Roll(angle)
        elif axis_ind == 1:
            camera.Azimuth(angle)
        else:
            camera.Pitch(angle)
        self.ren_win.Render()

    def write_image(
        self,
        filename: str | Path = "image.png",
        magnification: int = 1,
        image_format: str | None = None,
    ) -> None:
        """Save render window to an image.

        Args:
            filename: file to save to. Defaults to image.png.
            magnification: Use it to render high res images.
            image_format: One of png, jpeg (jpg), tiff (tif) or bmp. Defaults to
                None, which infers the format from the extension of filename
                and falls back to png.

        Raises:
            ValueError: If the image format is not supported.
        """
        writers = {
            "png": vtkPNGWriter,
            "jpg": vtkJPEGWriter,
            "jpeg": vtkJPEGWriter,
            "tif": vtkTIFFWriter,
            "tiff": vtkTIFFWriter,
            "bmp": vtkBMPWriter,
        }
        fmt = (image_format or Path(filename).suffix.lstrip(".") or "png").lower()
        if fmt not in writers:
            raise ValueError(f"Unsupported image format {fmt!r}. Choose from {sorted(writers)}.")
        writer = writers[fmt]()
        if fmt in {"jpg", "jpeg"}:
            writer.SetQuality(95)

        render_large = vtkRenderLargeImage()
        render_large.SetInput(self.ren)
        render_large.SetMagnification(magnification)
        writer.SetFileName(str(filename))
        writer.SetInputConnection(render_large.GetOutputPort())
        self.ren_win.Render()
        writer.Write()

    def redraw(self, reset_camera: bool = False) -> None:
        """
        Redraw the render window.

        Args:
            reset_camera: Set to True to reset the camera to a
                pre-determined default for each structure. Defaults to False.
        """
        self.ren.RemoveAllViewProps()
        self.picker = None
        self.add_picker_fixed()
        self.helptxt_mapper = vtkTextMapper()
        tprops = self.helptxt_mapper.GetTextProperty()
        tprops.SetFontSize(14)
        tprops.SetFontFamilyToArial()
        tprops.SetColor(0, 0, 0)
        self.helptxt_actor = vtkActor2D()
        self.helptxt_actor.SetMapper(self.helptxt_mapper)
        self.helptxt_actor.SetPosition(10, 10)
        self.helptxt_actor.VisibilityOff()

        if self.structure is not None:
            self.set_structure(self.structure, reset_camera)

        self.ren_win.Render()

    def orthogonalize_structure(self) -> None:
        """Orthogonalize the structure."""
        if self.structure is not None:
            self.set_structure(self.structure.copy(sanitize=True))
        self.ren_win.Render()

    def display_help(self) -> None:
        """Display the help for various keyboard shortcuts."""
        help_text = [
            "h : Toggle help",
            "A/a, B/b or C/c : Increase/decrease cell by one a, b or c unit vector",
            "# : Toggle showing of polyhedrons",
            "-: Toggle showing of bonds",
            "r : Reset camera direction",
            f"[/]: Decrease or increase poly_radii_tol_factor by 0.05. Value = {self.poly_radii_tol_factor}",
            "Up/Down: Rotate view along Up direction by 90 clockwise/anticlockwise",
            "Left/right: Rotate view along camera direction by 90 clockwise/anticlockwise",
            "s: Save view to image.png",
            "o: Orthogonalize structure",
        ]
        self.helptxt_mapper.SetInput("\n".join(help_text))
        self.helptxt_actor.SetPosition(10, 10)
        self.helptxt_actor.VisibilityOn()

    def set_structure(
        self,
        structure: Structure | Molecule,
        reset_camera: bool = True,
        to_unit_cell: bool = True,
    ) -> None:
        """
        Add a structure to the visualizer.

        Args:
            structure: structure to visualize
            reset_camera: Set to True to reset the camera to a default
                determined based on the structure.
            to_unit_cell: Whether or not to fall back sites into the unit cell.
        """
        self.ren.RemoveAllViewProps()

        has_lattice = hasattr(structure, "lattice")

        if has_lattice:
            struct = Structure.from_sites(structure, to_unit_cell=to_unit_cell)
            struct.make_supercell(self.supercell, to_unit_cell=to_unit_cell)
        else:
            struct = structure
        self._anion_symbols = _get_anion_symbols(struct)

        inc_coords = []
        for site in struct:
            # As in VESTA, sites on cell faces, edges and corners are drawn at all equivalent boundary positions.
            images = _get_boundary_images(site) if has_lattice else []
            for image in [site, *images]:
                self.add_site(image)
                inc_coords.append(image.coords)

        labels = ["a", "b", "c"]
        colors = [(1, 0, 0), (0, 0.6, 0), (0, 0, 1)]
        cell_color = (0.1, 0.1, 0.1)

        matrix = struct.lattice.matrix if has_lattice else None

        if self.show_unit_cell and has_lattice and matrix is not None:
            self.add_text([0, 0, 0], "o")
            for vec, label, color in zip(matrix, labels, colors, strict=True):
                self.add_line((0, 0, 0), vec, cell_color)
                self.add_text(vec, label, color)
            for vec1, vec2 in itertools.permutations(matrix, 2):
                self.add_line(vec1, vec1 + vec2, cell_color)
            for vec1, vec2, vec3 in itertools.permutations(matrix, 3):
                self.add_line(vec1 + vec2, vec1 + vec2 + vec3, cell_color)

        if self.show_bonds or self.show_polyhedron:
            elements = sorted(struct.elements, key=lambda a: a.X)
            anion = elements[-1]

            def contains_anion(site):
                return any(sp.symbol == anion.symbol for sp in site.species)

            anion_radius = anion.average_ionic_radius
            for site in struct:
                exclude = False
                max_radius = 0
                color = np.zeros(3)
                for sp, occu in site.species.items():
                    # Compare symbols so that oxidation-state-decorated anions are excluded too.
                    if sp.symbol in self.excluded_bonding_elements or sp.symbol == anion.symbol:
                        exclude = True
                        break
                    max_radius = max(max_radius, sp.average_ionic_radius)  # type:ignore[type-var,assignment]
                    color += occu * np.array(self.el_color_mapping.get(sp.symbol, [0, 0, 0]))

                if not exclude:
                    max_radius = (1 + self.poly_radii_tol_factor) * (max_radius + anion_radius)
                    nn = structure.get_neighbors(site, float(max_radius))
                    nn_sites = []
                    for neighbor in nn:
                        if contains_anion(neighbor):
                            nn_sites.append(neighbor)
                            if not in_coord_list(inc_coords, neighbor.coords):
                                self.add_site(neighbor)
                    if self.show_bonds:
                        self.add_bonds(nn_sites, site, radius=self.bond_radius)
                    if self.show_polyhedron:
                        self.add_polyhedron(
                            nn_sites,
                            site,
                            color / 255,
                            opacity=self.polyhedron_opacity,
                            draw_edges=True,
                            edges_color=(0.25, 0.25, 0.25),
                            edges_linewidth=1.5,
                        )

        self.ren.AddActor(self.helptxt_actor)
        if self.show_help:
            self.display_help()
        else:
            self.helptxt_actor.VisibilityOff()

        camera = self.ren.GetActiveCamera()
        if reset_camera:
            if has_lattice and matrix is not None:
                # Adjust the camera for best viewing
                lengths = struct.lattice.abc
                pos = (matrix[1] + matrix[2]) * 0.5 + matrix[0] * max(lengths) / lengths[0] * 3.5
                camera.SetPosition(pos)
                camera.SetViewUp(matrix[2])
                camera.SetFocalPoint((matrix[0] + matrix[1] + matrix[2]) * 0.5)
            else:
                origin = struct.center_of_mass
                max_site = max(struct, key=lambda site: site.distance_from_point(origin))
                camera.SetPosition(origin + 5 * (max_site.coords - origin))
                camera.SetFocalPoint(struct.center_of_mass)
            # Keep the viewing direction set above, but move the camera so everything fits in view.
            self.ren.ResetCamera()

        self.structure = structure
        self.title = struct.formula

    def zoom(self, factor: float) -> None:
        """Zoom the camera view by a factor."""
        camera = self.ren.GetActiveCamera()
        camera.Zoom(factor)
        self.ren_win.Render()

    def show(self) -> None:
        """Display the visualizer.

        Raises:
            RuntimeError: If the visualizer was created with offscreen=True.
        """
        if self.offscreen:
            raise RuntimeError("show() is unavailable for an offscreen visualizer. Use write_image() instead.")
        self.iren.Initialize()
        self.ren_win.SetWindowName(self.title)
        self.ren_win.Render()
        self.iren.Start()

    def get_site_radius(self, site: Site, anion_symbols: set[str] | None = None) -> float:
        """Radius (Å) of the sphere drawn for a site.

        The occupancy-weighted mean ionic radius of the species on the site,
        scaled by atom_radius_scale.

        Args:
            site: Site to size.
            anion_symbols: Symbols of elements treated as anions when they carry
                no oxidation state. Defaults to the most electronegative element
                of the last structure passed to set_structure.

        Returns:
            float: Sphere radius in Å.
        """
        anion_symbols = self._anion_symbols if anion_symbols is None else anion_symbols
        total_occu = sum(site.species.values())
        if total_occu == 0:
            return self.atom_radius_scale * DEFAULT_IONIC_RADIUS
        radius = sum(occu * get_ionic_radius(sp, anion_symbols) for sp, occu in site.species.items()) / total_occu
        return self.atom_radius_scale * radius

    def get_site_color(self, site: Site) -> list[float]:
        """RGB color (0-1) of the species with the highest occupancy on a site."""
        specie = max(site.species.items(), key=lambda item: item[1])[0]
        return [c / 255 for c in self.el_color_mapping.get(specie.symbol, (255, 255, 255))]

    def add_site(self, site: Site) -> None:
        """
        Add a site to the render window. The site is displayed as a sphere
        sized by the ionic radii of its species (see get_site_radius). Partially
        occupied sites are drawn as wedges of each species' color in proportion
        to its occupancy, with any vacancy fraction in white.

        Args:
            site: Site to add.
        """
        vis_radius = self.get_site_radius(site)
        start_angle = total_occu = 0

        for specie, occu in site.species.items():
            color = [i / 255 for i in self.el_color_mapping.get(specie.symbol, (255, 255, 255))]

            mapper = self.add_partial_sphere(site.coords, vis_radius, color, start_angle, start_angle + 360 * occu)
            self.mapper_map[mapper] = [site]
            start_angle += 360 * occu
            total_occu += occu

        if total_occu < 1:
            mapper = self.add_partial_sphere(
                site.coords,
                vis_radius,
                (1, 1, 1),
                start_angle,
                start_angle + 360 * (1 - total_occu),
            )
            self.mapper_map[mapper] = [site]

    def _add_actor(
        self,
        mapper: vtkPolyDataMapper,
        color: Sequence[float] | None = None,
        opacity: float | None = None,
        shiny: bool = True,
    ) -> vtkActor:
        """Add an actor for mapper with VESTA-like glossy shading."""
        actor = vtkActor()
        actor.SetMapper(mapper)
        prop = actor.GetProperty()
        if color is not None:
            prop.SetColor(color)
        if opacity is not None:
            prop.SetOpacity(opacity)
        if shiny:
            prop.SetInterpolationToPhong()
            prop.SetAmbient(0.15)
            prop.SetDiffuse(0.8)
            prop.SetSpecular(0.4)
            prop.SetSpecularPower(40)
        self.ren.AddActor(actor)
        return actor

    def add_partial_sphere(
        self,
        coords: ArrayLike,
        radius: float,
        color: Sequence[float],
        start: float = 0,
        end: float = 360,
        opacity: float = 1.0,
    ) -> vtkPolyDataMapper:
        """
        Adding a partial sphere (to display partial occupancies.

        Args:
            coords (NDArray): Coordinates
            radius (float): Radius of sphere
            color (tuple): RGB color of sphere
            start (float): Starting angle.
            end (float): Ending angle.
            opacity (float): Opacity.

        Returns:
            vtkPolyDataMapper: The mapper of the sphere.
        """
        sphere = vtkSphereSource()
        sphere.SetCenter(coords)
        sphere.SetRadius(radius)
        sphere.SetThetaResolution(24)
        sphere.SetPhiResolution(24)
        sphere.SetStartTheta(start)
        sphere.SetEndTheta(end)

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
        self._add_actor(mapper, color, opacity)
        return mapper

    def add_text(self, coords: ArrayLike, text: str, color: Sequence[float] = (0, 0, 0)) -> None:
        """
        Add text at a coordinate.

        Args:
            coords: Coordinates to add text at.
            text: Text to place.
            color: Color for text as RGB. Defaults to black.
        """
        source = vtkVectorText()
        source.SetText(text)
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        follower = vtkFollower()
        follower.SetMapper(mapper)
        follower.GetProperty().SetColor(color)
        follower.SetPosition(coords)
        follower.SetScale(0.5)
        self.ren.AddActor(follower)
        follower.SetCamera(self.ren.GetActiveCamera())

    def add_line(
        self,
        start: ArrayLike,
        end: ArrayLike,
        color: Sequence[float] = (0.5, 0.5, 0.5),
        width: float = 1,
    ) -> None:
        """
        Adds a line.

        Args:
            start: Starting coordinates for line.
            end: Ending coordinates for line.
            color: Color for text as RGB. Defaults to grey.
            width: Width of line. Defaults to 1.
        """
        source = vtkLineSource()
        source.SetPoint1(start)
        source.SetPoint2(end)

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        actor = self._add_actor(mapper, color, shiny=False)
        actor.GetProperty().SetLineWidth(width)

    def add_polyhedron(
        self,
        neighbors: Sequence[Site],
        center: Site,
        color: Sequence[float] | str,
        opacity: float = 1.0,
        draw_edges: bool = False,
        edges_color: Sequence[float] = (0.0, 0.0, 0.0),
        edges_linewidth: float = 2,
    ) -> None:
        """
        Adds a polyhedron, drawn as the convex hull of its vertices. Coplanar
        vertices are drawn as a flat polygon. Fewer than three non-collinear
        vertices draw nothing.

        Args:
            neighbors: Neighbors of the polyhedron (the vertices).
            center: The atom in the center of the polyhedron.
            color: Color for the polyhedron as RGB, or "element" to use the color
                of the species with the highest occupancy on the central atom.
            opacity: Opacity of the polyhedron
            draw_edges: If set to True, a line will be drawn at each edge
            edges_color: Color of the line for the edges
            edges_linewidth: Width of the line drawn for the edges
        """
        coords = np.array([neighbor.coords for neighbor in neighbors], dtype=float).reshape(-1, 3)
        if len(coords) < 3:
            return
        centroid = coords.mean(axis=0)
        _, singular_values, vt = np.linalg.svd(coords - centroid)
        rank = int(np.sum(singular_values > 1e-3))
        if rank < 2:
            return
        try:
            if rank == 3:
                polygons = []
                for simplex in ConvexHull(coords).simplices:
                    triangle = coords[simplex]
                    # Orient every face outwards so that coplanar neighbors share a normal.
                    normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
                    if np.dot(normal, triangle.mean(axis=0) - centroid) < 0:
                        triangle = triangle[::-1]
                    polygons.append(triangle)
            else:
                # Coplanar vertices: order them around their 2D convex hull in the best-fit plane.
                polygons = [coords[ConvexHull((coords - centroid) @ vt[:2].T).vertices]]
        except QhullError:
            return

        polydata = _make_polydata(polygons)
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(polydata)
        self.mapper_map[mapper] = [center, *neighbors]
        if isinstance(color, str) and color == "element":
            color = self.get_site_color(center)
        self._add_actor(mapper, color, opacity)

        if draw_edges:
            # Feature edges skip the diagonals that split coplanar faces into triangles.
            edges = vtkFeatureEdges()
            edges.SetInputData(polydata)
            edges.BoundaryEdgesOn()
            edges.FeatureEdgesOn()
            edges.SetFeatureAngle(10)
            edges.ManifoldEdgesOff()
            edges.NonManifoldEdgesOff()
            edges.ColoringOff()
            edges_mapper = vtkPolyDataMapper()
            edges_mapper.SetInputConnection(edges.GetOutputPort())
            actor = self._add_actor(edges_mapper, edges_color, shiny=False)
            actor.GetProperty().SetLineWidth(edges_linewidth)

    def _add_triangles(
        self,
        triangles: Sequence[Sequence[ArrayLike]],
        color: Sequence[float],
        opacity: float,
    ) -> vtkActor:
        """Add a set of triangles, each given as three Cartesian coordinates, as one actor."""
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(_make_polydata(triangles))
        return self._add_actor(mapper, color, opacity)

    def add_triangle(
        self,
        neighbors: Sequence[Site],
        color: Sequence[float] | str,
        center: Site | None = None,
        opacity: float = 0.4,
        draw_edges: bool = False,
        edges_color: Sequence[float] = (0.0, 0.0, 0.0),
        edges_linewidth: float = 2,
    ) -> None:
        """
        Adds a triangular surface between three atoms.

        Args:
            neighbors: Atoms between which a triangle will be drawn.
            color: Color for triangle as RGB, or "element" to use the color of
                the central atom.
            center: The "central atom" of the triangle
            opacity: opacity of the triangle
            draw_edges: If set to True, the a line will be drawn at each edge
            edges_color: Color of the line for the edges
            edges_linewidth: Width of the line drawn for the edges

        Raises:
            ValueError: If color is "element" and center is not provided.
        """
        if isinstance(color, str) and color == "element":
            if center is None:
                raise ValueError(
                    "Color should be chosen according to the central atom, and central atom is not provided"
                )
            color = self.get_site_color(center)
        actor = self._add_triangles([[neighbor.coords for neighbor in neighbors[:3]]], color, opacity)
        if draw_edges:
            actor.GetProperty().SetEdgeColor(edges_color)
            actor.GetProperty().SetLineWidth(edges_linewidth)
            actor.GetProperty().EdgeVisibilityOn()

    def add_faces(self, faces: Sequence[Sequence[ArrayLike]], color: Sequence[float], opacity: float = 0.35) -> None:
        """
        Adding face of polygon. Faces with more than three vertices are drawn as
        a fan of triangles around their centroid.

        Args:
            faces (list): Coordinates of the faces.
            color (tuple): RGB color.
            opacity (float): Opacity

        Raises:
            ValueError: If a face has fewer than three vertices.
        """
        triangles = []
        for face in faces:
            if len(face) < 3:
                raise ValueError("Number of points for a face should be >= 3")
            if len(face) == 3:
                triangles.append(face)
                continue
            center = np.mean(np.asarray(face, dtype=float), axis=0)
            triangles.extend([face[idx], face[(idx + 1) % len(face)], center] for idx in range(len(face)))
        if triangles:
            self._add_triangles(triangles, color, opacity)

    def add_edges(
        self,
        edges: Sequence[Sequence[Sequence[float]]],
        type: str = "line",  # noqa: A002
        linewidth: float = 2,
        color: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> None:
        """
        Args:
            edges (Sequence): List of edges. Each edge is a list of two points.
            type (str): Type of the edge. Defaults to "line". Unused.
            linewidth (float): Width of the line.
            color (tuple[float, float, float]): RGB color.
        """
        points = vtkPoints()
        lines = vtkCellArray()
        for iedge, edge in enumerate(edges):
            points.InsertPoint(2 * iedge, edge[0])
            points.InsertPoint(2 * iedge + 1, edge[1])
            lines.InsertNextCell(2)
            lines.InsertCellPoint(2 * iedge)
            lines.InsertCellPoint(2 * iedge + 1)
        polydata = vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(lines)
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(polydata)
        actor = self._add_actor(mapper, color, shiny=False)
        actor.GetProperty().SetLineWidth(linewidth)

    def add_bonds(
        self,
        neighbors: Sequence[Site],
        center: Site,
        color: Sequence[float] | None = None,
        opacity: float | None = None,
        radius: float = 0.1,
    ) -> None:
        """
        Adds bonds for a site.

        Args:
            neighbors: Neighbors of the site.
            center: The site in the center for all bonds.
            color: Color of the tubes representing the bonds. Defaults to None,
                which colors each half of a bond by the atom at that end, as in
                VESTA.
            opacity: Opacity of the tubes representing the bonds
            radius: Radius of tube s representing the bonds
        """
        center_coords = np.asarray(center.coords, dtype=float)
        center_rgb = [round(255 * c) for c in self.get_site_color(center)]
        points = vtkPoints()
        lines = vtkCellArray()
        colors = vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("Colors")
        for neighbor in neighbors:
            neighbor_coords = np.asarray(neighbor.coords, dtype=float)
            midpoint = (center_coords + neighbor_coords) / 2
            neighbor_rgb = [round(255 * c) for c in self.get_site_color(neighbor)]
            for start, end, rgb in ((center_coords, midpoint, center_rgb), (midpoint, neighbor_coords, neighbor_rgb)):
                lines.InsertNextCell(2)
                lines.InsertCellPoint(points.InsertNextPoint(start))
                lines.InsertCellPoint(points.InsertNextPoint(end))
                colors.InsertNextTuple3(*rgb)
        polydata = vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(lines)
        polydata.GetCellData().SetScalars(colors)

        tube = vtkTubeFilter()
        tube.SetInputData(polydata)
        tube.SetRadius(radius)
        tube.SetNumberOfSides(16)
        tube.CappingOn()

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(tube.GetOutputPort())
        if color is None:
            mapper.SetScalarModeToUseCellData()
        else:
            mapper.ScalarVisibilityOff()
        self._add_actor(mapper, color, opacity)

    def _format_site_info(self, mapper: vtkPolyDataMapper) -> str:
        """Species and fractional/Cartesian coordinates of the sites drawn by mapper."""
        output = []
        for site in self.mapper_map[mapper]:
            frac_coords = getattr(site, "frac_coords", None)
            row = [f"{site.species_string} - "]
            if frac_coords is not None:
                row.append(", ".join(f"{c:.3f}" for c in frac_coords))
            row.append("[" + ", ".join(f"{c:.3f}" for c in site.coords) + "]")
            output.append("".join(row))
        return "\n".join(output)

    def add_picker_fixed(self) -> None:
        """Create a cell picker that shows information on picked sites in the corner of the window."""
        picker = vtkCellPicker()

        # Create a Python function to create the text for the text mapper used
        # to display the results of picking.

        def annotate_pick(obj, event):
            if picker.GetCellId() < 0 and not self.show_help:
                self.helptxt_actor.VisibilityOff()
            else:
                mapper = picker.GetMapper()
                if mapper in self.mapper_map:
                    self.helptxt_mapper.SetInput(self._format_site_info(mapper))
                    self.helptxt_actor.SetPosition(10, 10)
                    self.helptxt_actor.VisibilityOn()
                    self.show_help = False

        self.picker = picker
        picker.AddObserver("EndPickEvent", annotate_pick)
        self.iren.SetPicker(picker)

    def add_picker(self) -> None:
        """Create a cell picker that shows information on picked sites next to them."""
        picker = vtkCellPicker()
        # Create a Python function to create the text for the text mapper used
        # to display the results of picking.
        source = vtkVectorText()
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        follower = vtkFollower()
        follower.SetMapper(mapper)
        follower.GetProperty().SetColor((0, 0, 0))
        follower.SetScale(0.2)
        self.ren.AddActor(follower)
        follower.SetCamera(self.ren.GetActiveCamera())
        follower.VisibilityOff()

        def annotate_pick(obj, event):
            if picker.GetCellId() < 0:
                follower.VisibilityOff()
            else:
                pick_pos = picker.GetPickPosition()
                mapper = picker.GetMapper()
                if mapper in self.mapper_map:
                    source.SetText(self._format_site_info(mapper))
                    follower.SetPosition(pick_pos)
                    follower.VisibilityOn()

        picker.AddObserver("EndPickEvent", annotate_pick)
        self.picker = picker
        self.iren.SetPicker(picker)


class StructureInteractorStyle(TrackballCamera):
    """A custom interactor style for visualizing structures."""

    def __init__(self, parent: StructureVis) -> None:
        """Initialize StructureInteractorStyle."""
        self.parent = parent
        self.mouse_motion = 0
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
        self.AddObserver("MouseMoveEvent", self.mouseMoveEvent)
        self.AddObserver("LeftButtonReleaseEvent", self.leftButtonReleaseEvent)
        self.AddObserver("KeyPressEvent", self.keyPressEvent)

    def leftButtonPressEvent(self, obj, event):
        """Left mouse button press event."""
        self.mouse_motion = 0
        self.OnLeftButtonDown()

    def mouseMoveEvent(self, obj, event):
        """Mouse move event."""
        self.mouse_motion = 1
        self.OnMouseMove()

    def leftButtonReleaseEvent(self, obj, event):
        """Left mouse button release event."""
        if self.mouse_motion == 0:
            iren = self.parent.iren
            pos = iren.GetEventPosition()
            iren.GetPicker().Pick(pos[0], pos[1], 0, self.parent.ren)
        self.OnLeftButtonUp()

    def keyPressEvent(self, obj, _event):
        """Key press event."""
        parent = self.parent
        sym = parent.iren.GetKeySym()

        if sym in {"A", "B", "C"}:
            axis = "ABC".index(sym)
            parent.supercell[axis][axis] += 1
            parent.redraw()
        elif sym in {"a", "b", "c"}:
            axis = "abc".index(sym)
            parent.supercell[axis][axis] = max(parent.supercell[axis][axis] - 1, 1)
            parent.redraw()
        elif sym == "numbersign":
            parent.show_polyhedron = not parent.show_polyhedron
            parent.redraw()
        elif sym == "minus":
            parent.show_bonds = not parent.show_bonds
            parent.redraw()
        elif sym == "bracketleft":
            parent.poly_radii_tol_factor -= 0.05 if parent.poly_radii_tol_factor > 0 else 0
            parent.redraw()
        elif sym == "bracketright":
            parent.poly_radii_tol_factor += 0.05
            parent.redraw()
        elif sym == "h":
            parent.show_help = not parent.show_help
            parent.redraw()
        elif sym == "r":
            parent.redraw(reset_camera=True)
        elif sym == "s":
            parent.write_image("image.png")
        elif sym == "Up":
            parent.rotate_view(1, 90)
        elif sym == "Down":
            parent.rotate_view(1, -90)
        elif sym == "Left":
            parent.rotate_view(0, -90)
        elif sym == "Right":
            parent.rotate_view(0, 90)
        elif sym == "o":
            parent.orthogonalize_structure()
            parent.redraw()

        self.OnKeyPress()


def make_movie(
    structures: Sequence[Structure],
    output_filename: str | Path = "movie.mp4",
    zoom: float = 1.0,
    fps: int = 20,
    bitrate: str = "10000k",
    quality: int = 1,
    **kwargs: Any,
) -> None:
    """
    Generate a movie from a sequence of structures using vtk and ffmpeg.

    Frames are rendered offscreen into a temporary directory. The camera is set
    from the first structure and kept fixed for the rest.

    Args:
        structures ([Structure]): sequence of structures
        output_filename (str): filename for structure output. defaults to
            movie.mp4
        zoom (float): A zoom to be applied to the visualizer. Defaults to 1.0.
        fps (int): Frames per second for the movie. Defaults to 20.
        bitrate (str): Video bitrate. Defaults to "10000k" (fairly high
            quality).
        quality (int): A quality scale. Defaults to 1.
        kwargs: Any kwargs supported by StructureVis to modify the images
            generated.

    Raises:
        RuntimeError: If ffmpeg is not on the PATH.
    """
    if shutil.which("ffmpeg") is None:
        raise RuntimeError("make_movie requires ffmpeg on the PATH.")
    kwargs.setdefault("offscreen", True)
    vis = StructureVis(**kwargs)
    vis.show_help = False
    vis.redraw()
    n_digits = len(str(len(structures)))
    with tempfile.TemporaryDirectory() as tmp_dir:
        for idx, struct in enumerate(structures):
            vis.set_structure(struct, reset_camera=idx == 0)
            if idx == 0:
                vis.zoom(zoom)
            vis.write_image(Path(tmp_dir) / f"image{idx:0{n_digits}d}.png", 3)
        args = [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            str(Path(tmp_dir) / f"image%0{n_digits}d.png"),
            "-q:v",
            str(quality),
            "-b:v",
            str(bitrate),
            str(output_filename),
        ]
        subprocess.run(args, check=True)


class MultiStructuresVis(StructureVis):
    """Visualization for multiple structures."""

    DEFAULT_ANIMATED_MOVIE_OPTIONS: ClassVar[dict[str, str | float]] = {
        "time_between_frames": 0.1,
        "looping_type": "restart",
        "number_of_loops": 1,
        "time_between_loops": 1.0,
    }

    def __init__(
        self,
        element_color_mapping: dict[str, Sequence[int]] | None = None,
        show_unit_cell: bool = True,
        show_bonds: bool = False,
        show_polyhedron: bool = False,
        poly_radii_tol_factor: float = 0.5,
        excluded_bonding_elements: Sequence[str] | None = None,
        animated_movie_options: dict[str, str | float] | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            element_color_mapping: Optional color mapping for the elements,
                as a dict of {symbol: rgb tuple}. For example, {"Fe": (255,
                123,0), ....} If None is specified, VESTA's color scheme is used.
            show_unit_cell: Set to False to not show the unit cell
                boundaries. Defaults to True.
            show_bonds: Set to True to show bonds. Defaults to False.
            show_polyhedron: Set to True to show polyhedrons. Defaults to
                False.
            poly_radii_tol_factor: The polyhedron and bonding code uses the
                ionic radii of the elements or species to determine if two
                atoms are bonded. This specifies a tolerance scaling factor
                such that atoms which are (1 + poly_radii_tol_factor) * sum
                of ionic radii apart are still considered as bonded.
            excluded_bonding_elements: List of atom types to exclude from
                bonding determination. Defaults to an empty list. Useful
                when trying to visualize a certain atom type in the
                framework (e.g., Li in a Li-ion battery cathode material).
            animated_movie_options (dict): Options for animated movie. Defaults
                to DEFAULT_ANIMATED_MOVIE_OPTIONS.
            kwargs: Other kwargs supported by StructureVis, e.g., offscreen.
        """
        self.structures: Sequence[Structure] = []
        self.tags: list[dict] = []
        self.all_radii: list[list[float]] = []
        self.all_vis_radii: list[list[float]] = []
        self.istruct = 0
        self.current_structure: Structure | None = None
        super().__init__(
            element_color_mapping=element_color_mapping,
            show_unit_cell=show_unit_cell,
            show_bonds=show_bonds,
            show_polyhedron=show_polyhedron,
            poly_radii_tol_factor=poly_radii_tol_factor,
            excluded_bonding_elements=excluded_bonding_elements,
            **kwargs,
        )
        self.warning_txt_actor = vtkActor2D()
        self.info_txt_actor = vtkActor2D()
        self.interactor_style = MultiStructuresInteractorStyle(self)
        self.iren.SetInteractorStyle(self.interactor_style)
        self.set_animated_movie_options(animated_movie_options=animated_movie_options)

    def set_structures(self, structures: Sequence[Structure], tags: list[dict] | None = None) -> None:
        """Add list of structures to the visualizer.

        Args:
            structures (list[Structures]): structures to be visualized.
            tags (list[dict]): List of tags to be applied to the structures.
        """
        self.structures = structures
        self.istruct = 0
        self.current_structure = self.structures[self.istruct]
        self.tags = tags if tags is not None else []
        self.all_radii = []
        self.all_vis_radii = []
        for struct in self.structures:
            anion_symbols = _get_anion_symbols(struct)
            struct_vis_radii = [self.get_site_radius(site, anion_symbols) for site in struct]
            self.all_vis_radii.append(struct_vis_radii)
            self.all_radii.append([radius / self.atom_radius_scale for radius in struct_vis_radii])

        self.set_structure(self.current_structure, reset_camera=True, to_unit_cell=False)

    def set_structure(
        self,
        structure: Structure | Molecule,
        reset_camera: bool = True,
        to_unit_cell: bool = False,
    ) -> None:
        """
        Add a structure to the visualizer.

        Args:
            structure: structure to visualize
            reset_camera: Set to True to reset the camera to a default
                determined based on the structure.
            to_unit_cell: Whether or not to fall back sites into the unit cell.
        """
        super().set_structure(structure=structure, reset_camera=reset_camera, to_unit_cell=to_unit_cell)
        self.apply_tags()

    def apply_tags(self) -> None:
        """Apply tags."""
        tags = {}
        for tag in self.tags:
            istruct = tag.get("istruct", "all")
            if istruct not in ("all", self.istruct):
                continue
            site_index = tag["site_index"]
            color = tag.get("color", [0.5, 0.5, 0.5])
            opacity = tag.get("opacity", 0.5)
            if site_index == "unit_cell_all":
                struct_radii = self.all_vis_radii[self.istruct]
                for site_idx in range(len(self.current_structure)):
                    vis_radius = 1.5 * tag.get("radius", struct_radii[site_idx])
                    tags[site_idx, (0, 0, 0)] = {
                        "radius": vis_radius,
                        "color": color,
                        "opacity": opacity,
                    }
                continue
            cell_index = tag["cell_index"]
            if "radius" in tag:
                vis_radius = tag["radius"]
            elif "radius_factor" in tag:
                vis_radius = tag["radius_factor"] * self.all_vis_radii[self.istruct][site_index]
            else:
                vis_radius = 1.5 * self.all_vis_radii[self.istruct][site_index]
            tags[site_index, cell_index] = {
                "radius": vis_radius,
                "color": color,
                "opacity": opacity,
            }
        for site_and_cell_index, tag_style in tags.items():
            site_idx, cell_index = site_and_cell_index
            site = self.current_structure[site_idx]
            if cell_index == (0, 0, 0):
                coords = site.coords
            else:
                frac_coords = site.frac_coords + np.array(cell_index)
                site_image = PeriodicSite(
                    site.species,
                    frac_coords,
                    self.current_structure.lattice,
                    to_unit_cell=False,
                    coords_are_cartesian=False,
                    properties=site.properties,
                )
                self.add_site(site_image)
                coords = site_image.coords
            vis_radius = tag_style["radius"]
            color = tag_style["color"]
            opacity = tag_style["opacity"]
            self.add_partial_sphere(
                coords=coords,
                radius=vis_radius,
                color=color,
                start=0,
                end=360,
                opacity=opacity,
            )

    def set_animated_movie_options(self, animated_movie_options: dict[str, str | float] | None = None) -> None:
        """
        Args:
            animated_movie_options (dict): Options for animated movie.

        Raises:
            ValueError: If an option is not one of DEFAULT_ANIMATED_MOVIE_OPTIONS.
        """
        self.animated_movie_options = self.DEFAULT_ANIMATED_MOVIE_OPTIONS.copy()
        if animated_movie_options is None:
            return
        for key in animated_movie_options:
            if key not in self.DEFAULT_ANIMATED_MOVIE_OPTIONS:
                raise ValueError("Wrong option for animated movie")
        self.animated_movie_options.update(animated_movie_options)

    def display_help(self) -> None:
        """Display the help for various keyboard shortcuts."""
        helptxt = [
            "h : Toggle help",
            "A/a, B/b or C/c : Increase/decrease cell by one a, b or c unit vector",
            "# : Toggle showing of polyhedrons",
            "-: Toggle showing of bonds",
            "r : Reset camera direction",
            f"[/]: Decrease or increase poly_radii_tol_factor by 0.05. Value = {self.poly_radii_tol_factor}",
            "Up/Down: Rotate view along Up direction by 90 clockwise/anticlockwise",
            "Left/right: Rotate view along camera direction by 90 clockwise/anticlockwise",
            "s: Save view to image.png",
            "o: Orthogonalize structure",
            "n: Move to next structure",
            "p: Move to previous structure",
            "m: Animated movie of the structures",
        ]
        self.helptxt_mapper.SetInput("\n".join(helptxt))
        self.helptxt_actor.SetPosition(10, 10)
        self.helptxt_actor.VisibilityOn()

    def display_warning(self, warning: str) -> None:
        """
        Args:
            warning (str): Warning.
        """
        self.warning_txt_mapper = vtkTextMapper()
        text_props = self.warning_txt_mapper.GetTextProperty()
        text_props.SetFontSize(14)
        text_props.SetFontFamilyToArial()
        text_props.SetColor(1, 0, 0)
        text_props.BoldOn()
        text_props.SetJustificationToRight()
        self.warning_txt = f"WARNING : {warning}"
        self.warning_txt_actor = vtkActor2D()
        self.warning_txt_actor.VisibilityOn()
        self.warning_txt_actor.SetMapper(self.warning_txt_mapper)
        self.ren.AddActor(self.warning_txt_actor)
        self.warning_txt_mapper.SetInput(self.warning_txt)
        win_size = self.ren_win.GetSize()
        self.warning_txt_actor.SetPosition(win_size[0] - 10, 10)
        self.warning_txt_actor.VisibilityOn()

    def erase_warning(self) -> None:
        """Remove warnings."""
        self.warning_txt_actor.VisibilityOff()

    def display_info(self, info: str) -> None:
        """
        Args:
            info (str): Information.
        """
        self.info_txt_mapper = vtkTextMapper()
        t_prop = self.info_txt_mapper.GetTextProperty()
        t_prop.SetFontSize(14)
        t_prop.SetFontFamilyToArial()
        t_prop.SetColor(0, 0, 1)
        t_prop.BoldOn()
        t_prop.SetVerticalJustificationToTop()
        self.info_txt = f"INFO : {info}"
        self.info_txt_actor = vtkActor2D()
        self.info_txt_actor.VisibilityOn()
        self.info_txt_actor.SetMapper(self.info_txt_mapper)
        self.ren.AddActor(self.info_txt_actor)
        self.info_txt_mapper.SetInput(self.info_txt)
        win_size = self.ren_win.GetSize()
        self.info_txt_actor.SetPosition(10, win_size[1] - 10)
        self.info_txt_actor.VisibilityOn()

    def erase_info(self) -> None:
        """Erase all info."""
        self.info_txt_actor.VisibilityOff()


class MultiStructuresInteractorStyle(StructureInteractorStyle):
    """Interactor for MultiStructureVis."""

    def __init__(self, parent: MultiStructuresVis) -> None:
        """Initialize MultiStructuresInteractorStyle."""
        StructureInteractorStyle.__init__(self, parent=parent)

    def _show_structure(self, istruct: int) -> None:
        parent = self.parent
        parent.istruct = istruct
        parent.current_structure = parent.structures[istruct]
        parent.set_structure(parent.current_structure, reset_camera=False, to_unit_cell=False)

    def keyPressEvent(self, obj, event):
        """Key press event."""
        parent = self.parent
        sym = parent.iren.GetKeySym()

        if sym == "n":
            if parent.istruct == len(parent.structures) - 1:
                parent.display_warning("LAST STRUCTURE")
            else:
                self._show_structure(parent.istruct + 1)
                parent.erase_warning()
            parent.ren_win.Render()
        elif sym == "p":
            if parent.istruct == 0:
                parent.display_warning("FIRST STRUCTURE")
            else:
                self._show_structure(parent.istruct - 1)
                parent.erase_warning()
            parent.ren_win.Render()
        elif sym == "m":
            n_structs = len(parent.structures)
            looping_type = parent.animated_movie_options["looping_type"]
            if looping_type == "restart":
                loop_istructs = list(range(n_structs))
            elif looping_type == "palindrome":
                loop_istructs = list(range(n_structs)) + list(range(n_structs - 2, -1, -1))
            else:
                raise ValueError('"looping_type" should be "restart" or "palindrome"')
            self._show_structure(0)
            parent.erase_warning()
            parent.ren_win.Render()
            nloops = parent.animated_movie_options["number_of_loops"]
            tstep = parent.animated_movie_options["time_between_frames"]
            tloops = parent.animated_movie_options["time_between_loops"]
            for iloop in range(nloops):
                for istruct in loop_istructs:
                    time.sleep(tstep)
                    self._show_structure(istruct)
                    parent.display_info(
                        f"Animated movie : structure {istruct + 1}/{n_structs} (loop {iloop + 1}/{nloops})"
                    )
                    parent.ren_win.Render()
                time.sleep(tloops)
            parent.erase_info()
            parent.display_info("Ended animated movie ...")
            parent.ren_win.Render()

        StructureInteractorStyle.keyPressEvent(self, obj, event)
