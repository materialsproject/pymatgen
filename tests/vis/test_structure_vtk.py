from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pymatgen.core import DummySpecies, Element, Lattice, Molecule, PeriodicSite, Species, Structure
from pymatgen.vis import structure_vtk
from pymatgen.vis.structure_vtk import (
    ATOM_RADIUS_SCALE,
    DEFAULT_IONIC_RADIUS,
    MultiStructuresVis,
    StructureVis,
    _get_boundary_images,
    get_ionic_radius,
    make_movie,
    vtk_available,
)
from tests.testing import TEST_FILES_DIR

requires_vtk = pytest.mark.skipif(not vtk_available, reason="vtk is not installed")

IMAGE_MAGIC = {".png": b"\x89PNG", ".jpg": b"\xff\xd8", ".bmp": b"BM", ".tiff": (b"II*\x00", b"MM\x00*")}


@pytest.fixture
def nacl() -> Structure:
    """Conventional rock-salt NaCl with oxidation states."""
    struct = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.64), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    struct.add_oxidation_state_by_element({"Na": 1, "Cl": -1})
    return struct


@pytest.fixture
def lfp() -> Structure:
    return Structure.from_file(TEST_FILES_DIR / "io/vasp/inputs/POSCAR_LiFePO4")


@pytest.fixture
def vis() -> StructureVis:
    return StructureVis(offscreen=True)


def _polyhedra(vis: StructureVis) -> list[list]:
    """Sites (center first) of every polyhedron drawn, keyed off the picker map."""
    return [sites for sites in vis.mapper_map.values() if len(sites) > 1]


def _sites(xyz: list[list[float]]) -> list[PeriodicSite]:
    lattice = Lattice.cubic(10)
    return [PeriodicSite("O", coords, lattice, coords_are_cartesian=True) for coords in xyz]


def _last_actor(vis: StructureVis):
    actors = vis.ren.GetActors()
    return actors.GetItemAsObject(actors.GetNumberOfItems() - 1)


def _n_actors(vis: StructureVis) -> int:
    return vis.ren.GetActors().GetNumberOfItems()


class TestGetIonicRadius:
    def test_species_uses_shannon_radius(self):
        assert get_ionic_radius(Species("Fe", 3)) == pytest.approx(float(Species("Fe", 3).ionic_radius))

    @pytest.mark.parametrize(
        ("anion_symbols", "attr"),
        [({"Cl"}, "average_anionic_radius"), (set(), "average_cationic_radius"), (None, "average_cationic_radius")],
    )
    def test_element_anion_or_cation(self, anion_symbols, attr):
        assert get_ionic_radius(Element("Cl"), anion_symbols) == pytest.approx(float(getattr(Element("Cl"), attr)))

    def test_falls_back_to_atomic_radius(self):
        # Na has no anionic radius, so an "anion" Na falls back to its average ionic radius.
        assert get_ionic_radius(Element("Na"), {"Na"}) == pytest.approx(float(Element("Na").average_ionic_radius))

    @pytest.mark.parametrize("specie", [Element("He"), DummySpecies("X", 1)])
    def test_default(self, specie):
        assert get_ionic_radius(specie) == DEFAULT_IONIC_RADIUS


class TestBoundaryImages:
    @pytest.mark.parametrize(
        ("frac_coords", "n_images"),
        [([0, 0, 0], 7), ([0.5, 0.5, 0], 1), ([0, 0.5, 0.5], 1), ([0.5, 0, 0], 3), ([0.25, 0.5, 0.75], 0)],
    )
    def test_count(self, frac_coords, n_images):
        site = PeriodicSite("Na", frac_coords, Lattice.cubic(4))
        images = _get_boundary_images(site)
        assert len(images) == n_images
        for image in images:
            assert np.all(np.abs(image.frac_coords - 0.5) <= 0.5 + 1e-8)
            assert image.species == site.species

    def test_upper_boundary(self):
        (image,) = _get_boundary_images(PeriodicSite("Na", [0.5, 0.5, 1], Lattice.cubic(4), to_unit_cell=False))
        assert image.frac_coords == pytest.approx([0.5, 0.5, 0])


@requires_vtk
class TestStructureVis:
    def test_site_radius(self, vis, nacl):
        vis.set_structure(nacl)
        for site in nacl:
            expected = ATOM_RADIUS_SCALE * float(site.specie.ionic_radius)
            assert vis.get_site_radius(site) == pytest.approx(expected)

    def test_site_radius_partial_occupancy(self, vis):
        lattice = Lattice.cubic(4)
        full = PeriodicSite(Species("Li", 1), [0, 0, 0], lattice)
        half = PeriodicSite({Species("Li", 1): 0.5}, [0, 0, 0], lattice)
        # A vacancy fraction does not shrink the sphere.
        assert vis.get_site_radius(half) == pytest.approx(vis.get_site_radius(full))

    def test_atom_radius_scale(self, nacl):
        vis = StructureVis(offscreen=True, atom_radius_scale=1.0)
        assert vis.get_site_radius(nacl[0]) == pytest.approx(float(nacl[0].specie.ionic_radius))

    def test_boundary_sites_drawn(self, nacl):
        vis = StructureVis(offscreen=True, show_polyhedron=False)
        vis.set_structure(nacl)
        # 8 corner + 6 face-center Na, and 12 edge-center + 1 body-center Cl.
        assert len(vis.mapper_map) == 27

    def test_polyhedra_exclude_decorated_anion(self, vis, nacl):
        vis.set_structure(nacl)
        polyhedra = _polyhedra(vis)
        assert len(polyhedra) == 4
        for center, *neighbors in polyhedra:
            assert center.specie.symbol == "Na"
            assert len(neighbors) == 6

    def test_polyhedra_excluded_bonding_elements(self, lfp):
        vis = StructureVis(offscreen=True, excluded_bonding_elements=["Li"])
        vis.set_structure(lfp)
        centers = [sites[0].specie.symbol for sites in _polyhedra(vis)]
        assert sorted(set(centers)) == ["Fe", "P"]

    def test_polyhedra_partial_occupancy(self, vis):
        struct = Structure.from_spacegroup(
            "Fm-3m", Lattice.cubic(4.2), [{"Li": 0.5, "Ni": 0.3}, "O"], [[0, 0, 0], [0.5, 0.5, 0.5]]
        )
        vis.set_structure(struct)
        assert len(_polyhedra(vis)) == 4

    def test_bonds(self, nacl):
        vis = StructureVis(offscreen=True, show_polyhedron=False, show_bonds=True, bond_radius=0.2)
        vis.set_structure(nacl)
        mapper = _last_actor(vis).GetMapper()
        tube = mapper.GetInputAlgorithm()
        assert tube.GetRadius() == pytest.approx(0.2)
        assert mapper.GetScalarVisibility()

    @pytest.mark.parametrize("color", [None, (0, 1, 0)])
    def test_add_bonds(self, vis, nacl, color):
        vis.set_structure(nacl)
        neighbors = nacl.get_neighbors(nacl[0], 3)
        vis.add_bonds(neighbors, nacl[0], color=color, opacity=0.5)
        actor = _last_actor(vis)
        assert bool(actor.GetMapper().GetScalarVisibility()) is (color is None)
        assert actor.GetProperty().GetOpacity() == pytest.approx(0.5)
        if color is None:
            # Two half-bond segments, each colored by the atom at its end.
            polydata = actor.GetMapper().GetInputAlgorithm().GetInput()
            colors = polydata.GetCellData().GetScalars()
            assert polydata.GetNumberOfLines() == colors.GetNumberOfTuples() == 2 * len(neighbors)
            na_rgb = [round(255 * c) for c in vis.get_site_color(nacl[0])]
            cl_rgb = [round(255 * c) for c in vis.get_site_color(neighbors[0])]
            assert list(colors.GetTuple3(0)) == na_rgb
            assert list(colors.GetTuple3(1)) == cl_rgb

    @pytest.mark.parametrize(
        ("vertices", "n_faces", "n_edges"),
        [
            ([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], 8, 12),
            ([[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)], 12, 12),
            ([[1, 1, 0], [-1, 1, 0], [-1, -1, 0], [1, -1, 0]], 1, 4),
        ],
        ids=["octahedron", "cube", "square"],
    )
    def test_add_polyhedron(self, vis, vertices, n_faces, n_edges):
        center = _sites([[0, 0, 0]])[0]
        vis.add_polyhedron(_sites(vertices), center, "element", opacity=0.3, draw_edges=True)
        edges_actor = _last_actor(vis)
        edges = edges_actor.GetMapper().GetInputAlgorithm()
        edges.Update()
        # Feature edges skip the diagonals that triangulate the square faces of the cube.
        assert edges.GetOutput().GetNumberOfLines() == n_edges
        (sites,) = _polyhedra(vis)
        assert sites[0] is center
        surface_mapper = next(mapper for mapper, mapped in vis.mapper_map.items() if len(mapped) > 1)
        assert surface_mapper.GetInput().GetNumberOfPolys() == n_faces

    @pytest.mark.parametrize("vertices", [[[1, 0, 0], [2, 0, 0], [3, 0, 0]], [[1, 0, 0], [2, 0, 0]]])
    def test_add_polyhedron_degenerate(self, vis, vertices):
        n_actors = _n_actors(vis)
        vis.add_polyhedron(_sites(vertices), _sites([[0, 0, 0]])[0], (1, 0, 0))
        assert _n_actors(vis) == n_actors

    def test_add_triangle(self, vis):
        sites = _sites([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        n_actors = _n_actors(vis)
        vis.add_triangle(sites, (1, 0, 0), draw_edges=True)
        vis.add_triangle(sites, "element", center=sites[0])
        assert _n_actors(vis) == n_actors + 2
        with pytest.raises(ValueError, match="central atom is not provided"):
            vis.add_triangle(sites, "element")

    def test_add_faces(self, vis):
        n_actors = _n_actors(vis)
        triangle = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        square = [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
        vis.add_faces([triangle, square], (0, 0, 1))
        assert _n_actors(vis) == n_actors + 1
        # The square is drawn as a fan of 4 triangles around its centroid.
        assert _last_actor(vis).GetMapper().GetInput().GetNumberOfPolys() == 5
        with pytest.raises(ValueError, match="should be >= 3"):
            vis.add_faces([[[0, 0, 0], [1, 0, 0]]], (0, 0, 1))

    def test_add_edges_line_text(self, vis):
        n_actors = _n_actors(vis)
        vis.add_edges([[[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [0, 1, 0]]], color=(1, 0, 0))
        assert _last_actor(vis).GetMapper().GetInput().GetNumberOfLines() == 2
        vis.add_line([0, 0, 0], [0, 0, 1], width=3)
        assert _last_actor(vis).GetProperty().GetLineWidth() == pytest.approx(3)
        vis.add_text([0, 0, 0], "X")
        assert _n_actors(vis) == n_actors + 3

    @pytest.mark.parametrize(
        ("suffix", "image_format"), [(".png", None), (".jpg", None), (".bmp", None), (".tiff", None), (".dat", "png")]
    )
    def test_write_image(self, vis, nacl, suffix, image_format):
        vis.set_structure(nacl)
        path = Path(f"image{suffix}")
        vis.write_image(path, image_format=image_format)
        magic = IMAGE_MAGIC.get(suffix, IMAGE_MAGIC[".png"])
        assert path.read_bytes().startswith(magic)

    def test_write_image_magnification(self, vis, nacl):
        vis.set_structure(nacl)
        vis.write_image("image.png", magnification=2)
        header = Path("image.png").read_bytes()[16:24]
        width, height = int.from_bytes(header[:4], "big"), int.from_bytes(header[4:], "big")
        assert (width, height) == (1600, 1600)

    def test_write_image_bad_format(self, vis):
        with pytest.raises(ValueError, match="Unsupported image format 'gif'"):
            vis.write_image("image.gif")

    def test_molecule(self, vis):
        mol = Molecule(["C", "O", "O"], [[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]])
        vis.set_structure(mol)
        assert vis.title == "C1 O2"
        assert len(vis.mapper_map) == 3

    def test_camera_and_view(self, vis, nacl):
        vis.set_structure(nacl)
        camera = vis.ren.GetActiveCamera()
        view_angle = camera.GetViewAngle()
        vis.zoom(2)
        assert camera.GetViewAngle() < view_angle
        view_up = camera.GetViewUp()
        vis.rotate_view(0, 90)
        assert camera.GetViewUp() != pytest.approx(view_up)
        position = camera.GetPosition()
        vis.rotate_view(1, 90)
        assert camera.GetPosition() != pytest.approx(position)
        focal_point = camera.GetFocalPoint()
        vis.rotate_view(2, 30)
        assert camera.GetFocalPoint() != pytest.approx(focal_point)

    def test_orthogonalize_structure(self, vis, nacl):
        vis.orthogonalize_structure()  # no structure is a no-op
        vis.set_structure(nacl)
        vis.orthogonalize_structure()
        assert vis.structure.lattice.abc == pytest.approx(nacl.lattice.abc)

    def test_show_offscreen(self, vis):
        with pytest.raises(RuntimeError, match="offscreen"):
            vis.show()

    def test_help(self, vis, nacl):
        vis.set_structure(nacl)
        assert vis.helptxt_actor.GetVisibility()
        assert "Toggle help" in vis.helptxt_mapper.GetInput()
        vis.show_help = False
        vis.redraw()
        assert not vis.helptxt_actor.GetVisibility()

    def test_pick_with_help_off(self, vis, nacl):
        vis.show_help = False
        vis.set_structure(nacl)
        vis.picker.InvokeEvent("EndPickEvent")
        assert not vis.helptxt_actor.GetVisibility()

    def test_format_site_info(self, vis, nacl):
        vis.set_structure(nacl)
        mapper = next(iter(vis.mapper_map))
        assert vis._format_site_info(mapper) == "Na+ - 0.000, 0.000, 0.000[0.000, 0.000, 0.000]"
        vis.mapper_map[mapper] = [Molecule(["C"], [[0, 0, 0]])[0]]
        assert vis._format_site_info(mapper) == "C - [0.000, 0.000, 0.000]"

    def test_pick_site(self, nacl):
        vis = StructureVis(offscreen=True, show_polyhedron=False, show_unit_cell=False)
        vis.set_structure(nacl)
        vis.ren_win.Render()
        site_infos = {vis._format_site_info(mapper) for mapper in vis.mapper_map}
        # An atom on the a-axis line of sight through the cell center fills the middle of the default view.
        vis.picker.Pick(400, 400, 0, vis.ren)
        assert vis.helptxt_mapper.GetInput() in site_infos
        assert not vis.show_help

    def test_add_picker(self, nacl):
        vis = StructureVis(offscreen=True, show_polyhedron=False, show_unit_cell=False)
        vis.set_structure(nacl)
        vis.add_picker()
        vis.ren_win.Render()
        vis.picker.Pick(400, 400, 0, vis.ren)
        follower = _last_actor(vis)
        assert follower.GetVisibility()
        site_infos = {vis._format_site_info(mapper) for mapper in vis.mapper_map}
        assert follower.GetMapper().GetInputAlgorithm().GetText() in site_infos
        vis.picker.Pick(0, 0, 0, vis.ren)
        assert not follower.GetVisibility()

    def test_elemental_structure(self, vis):
        fe = Structure(Lattice.cubic(2.87), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        vis.set_structure(fe)
        expected = ATOM_RADIUS_SCALE * float(Element("Fe").average_cationic_radius)
        assert vis.get_site_radius(fe[0]) == pytest.approx(expected)
        assert not _polyhedra(vis)

    def test_mouse_pick(self, vis, nacl):
        vis.set_structure(nacl)
        style = vis.interactor_style
        style.leftButtonPressEvent(style, "LeftButtonPressEvent")
        assert style.mouse_motion == 0
        style.leftButtonReleaseEvent(style, "LeftButtonReleaseEvent")
        style.mouseMoveEvent(style, "MouseMoveEvent")
        assert style.mouse_motion == 1

    def test_key_presses(self, vis, nacl):
        vis.set_structure(nacl)
        style = vis.interactor_style

        def press(sym):
            vis.iren.SetKeySym(sym)
            style.keyPressEvent(style, "KeyPressEvent")

        for sym, axis in zip("ABC", range(3), strict=True):
            press(sym)
            assert vis.supercell[axis][axis] == 2
            press(sym.lower())
            press(sym.lower())
            assert vis.supercell[axis][axis] == 1
        press("numbersign")
        assert not vis.show_polyhedron
        press("minus")
        assert vis.show_bonds
        press("bracketright")
        assert vis.poly_radii_tol_factor == pytest.approx(0.55)
        press("bracketleft")
        assert vis.poly_radii_tol_factor == pytest.approx(0.5)
        press("h")
        assert not vis.show_help
        press("s")
        assert Path("image.png").is_file()
        for sym in ("r", "Up", "Down", "Left", "Right", "o", ""):
            press(sym)

    def test_offscreen(self):
        # An offscreen visualizer renders without a window, so it needs no display.
        assert StructureVis(offscreen=True).ren_win.GetOffScreenRendering()


@requires_vtk
class TestMultiStructuresVis:
    @pytest.fixture
    def structures(self, nacl) -> list[Structure]:
        strained = nacl.copy()
        strained.apply_strain(0.02)
        return [nacl, strained, nacl.copy()]

    @pytest.fixture
    def multi(self, monkeypatch) -> MultiStructuresVis:
        monkeypatch.setattr(structure_vtk.time, "sleep", lambda _seconds: None)
        return MultiStructuresVis(offscreen=True, animated_movie_options={"number_of_loops": 2})

    def test_set_structure_before_set_structures(self, multi, nacl):
        multi.set_structure(nacl)
        assert multi.structure is nacl

    def test_set_structures_with_tags(self, multi, structures):
        tags = [
            {"site_index": 0, "cell_index": (1, 0, 0), "color": [1, 0, 0]},
            {"site_index": 1, "cell_index": (0, 0, 0), "radius_factor": 2.0},
            {"site_index": 2, "cell_index": (0, 0, 0), "radius": 0.3, "istruct": 0},
            {"site_index": "unit_cell_all", "istruct": 0},
            {"site_index": 0, "cell_index": (0, 0, 0), "istruct": 1},
        ]
        multi.set_structures(structures, tags=tags)
        assert len(multi.all_vis_radii) == len(structures)
        assert multi.all_vis_radii[0][0] == pytest.approx(multi.get_site_radius(structures[0][0]))
        assert multi.all_radii[0][0] == pytest.approx(float(structures[0][0].specie.ionic_radius))

    def test_animated_movie_options(self, multi):
        assert multi.animated_movie_options["number_of_loops"] == 2
        assert multi.animated_movie_options["looping_type"] == "restart"
        multi.set_animated_movie_options()
        assert multi.animated_movie_options == MultiStructuresVis.DEFAULT_ANIMATED_MOVIE_OPTIONS
        with pytest.raises(ValueError, match="Wrong option for animated movie"):
            multi.set_animated_movie_options({"bad_option": 1})

    @pytest.mark.parametrize("looping_type", ["restart", "palindrome"])
    def test_key_presses(self, multi, structures, looping_type):
        multi.set_structures(structures)
        multi.set_animated_movie_options({"looping_type": looping_type})
        style = multi.interactor_style

        def press(sym):
            multi.iren.SetKeySym(sym)
            style.keyPressEvent(style, "KeyPressEvent")

        press("p")
        assert multi.istruct == 0
        assert multi.warning_txt == "WARNING : FIRST STRUCTURE"
        press("n")
        press("n")
        assert multi.istruct == 2
        press("n")
        assert multi.warning_txt == "WARNING : LAST STRUCTURE"
        press("p")
        assert multi.istruct == 1
        assert not multi.warning_txt_actor.GetVisibility()
        press("m")
        assert multi.istruct == (2 if looping_type == "restart" else 0)
        assert multi.info_txt == "INFO : Ended animated movie ..."

    def test_bad_looping_type(self, multi, structures):
        multi.set_structures(structures)
        multi.animated_movie_options["looping_type"] = "bounce"
        multi.iren.SetKeySym("m")
        style = multi.interactor_style
        with pytest.raises(ValueError, match="looping_type"):
            style.keyPressEvent(style, "KeyPressEvent")

    def test_info_and_help(self, multi, structures):
        multi.set_structures(structures)
        assert "Move to next structure" in multi.helptxt_mapper.GetInput()
        multi.display_info("hello")
        assert multi.info_txt_actor.GetVisibility()
        multi.erase_info()
        assert not multi.info_txt_actor.GetVisibility()


@requires_vtk
class TestMakeMovie:
    def test_make_movie(self, monkeypatch, nacl):
        calls = {}

        def fake_run(args, check):
            pattern = Path(args[args.index("-i") + 1])
            calls["frames"] = sorted(path.name for path in pattern.parent.glob("image*.png"))
            calls["args"] = args
            calls["check"] = check

        monkeypatch.setattr(structure_vtk.shutil, "which", lambda _name: "/usr/bin/ffmpeg")
        monkeypatch.setattr(structure_vtk.subprocess, "run", fake_run)
        structures = [nacl] * 10
        make_movie(structures, output_filename="out.mp4", fps=10, show_polyhedron=False)
        assert calls["frames"] == [f"image{idx:02d}.png" for idx in range(10)]
        assert calls["args"][-1] == "out.mp4"
        assert calls["args"][calls["args"].index("-framerate") + 1] == "10"
        assert calls["check"] is True
        # Frames are written to a temporary directory that is cleaned up.
        assert not list(Path().glob("image*.png"))

    def test_make_movie_requires_ffmpeg(self, monkeypatch, nacl):
        monkeypatch.setattr(structure_vtk.shutil, "which", lambda _name: None)
        with pytest.raises(RuntimeError, match="ffmpeg"):
            make_movie([nacl])
