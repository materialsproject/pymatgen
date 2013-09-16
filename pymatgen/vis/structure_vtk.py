#!/usr/bin/env python

"""
This module contains classes to wrap Python VTK to make nice molecular plots.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Nov 27, 2011"

import os
import ConfigParser
import itertools
import math
import subprocess

import numpy as np
import vtk

from pymatgen.util.coord_utils import in_coord_list
from pymatgen.core.periodic_table import Specie
from pymatgen.core.structure import Structure


class StructureVis(object):
    """
    Provides Structure object visualization using VTK.
    """

    def __init__(self, element_color_mapping=None, show_unit_cell=True,
                 show_bonds=False, show_polyhedron=True,
                 poly_radii_tol_factor=0.5, excluded_bonding_elements=None):
        """
        Args:
            element_color_mapping:
                Optional color mapping for the elements, as a dict of
                {symbol: rgb tuple}. For example, {"Fe": (255,123,0), ....}
                If None is specified, a default based on Jmol"s color scheme
                is used.
            show_unit_cell:
                Set to False to not show the unit cell boundaries. Defaults to
                True.
            show_bonds:
                Set to True to show bonds. Defaults to True.
            show_polyhedron:
                Set to True to show polyhedrons. Defaults to False.
            poly_radii_tol_factor:
                The polyhedron and bonding code uses the ionic radii of the
                elements or species to determine if two atoms are bonded. This
                specifies a tolerance scaling factor such that atoms which are
                (1 + poly_radii_tol_factor) * sum of ionic radii apart are
                still considered as bonded.
            excluded_bonding_elements:
                List of atom types to exclude from bonding determination.
                Defaults to an empty list. Useful when trying to visualize a
                certain atom type in the framework (e.g., Li in a Li-ion
                battery cathode material).

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
        self.ren = vtk.vtkRenderer()
        self.ren_win = vtk.vtkRenderWindow()
        self.ren_win.AddRenderer(self.ren)
        self.ren.SetBackground(1, 1, 1)
        self.title = "Structure Visualizer"
        # create a renderwindowinteractor
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.ren_win)
        self.mapper_map = {}
        self.structure = None

        if element_color_mapping:
            self.el_color_mapping = element_color_mapping
        else:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            config = ConfigParser.SafeConfigParser()
            config.optionxform = str
            config.readfp(open(os.path.join(module_dir,
                                            "ElementColorSchemes.cfg")))

            self.el_color_mapping = {}
            for (el, color) in config.items("VESTA"):
                self.el_color_mapping[el] = [int(i) for i in color.split(",")]
        self.show_unit_cell = show_unit_cell
        self.show_bonds = show_bonds
        self.show_polyhedron = show_polyhedron
        self.poly_radii_tol_factor = poly_radii_tol_factor
        self.excluded_bonding_elements = excluded_bonding_elements if \
            excluded_bonding_elements else []
        self.show_help = True
        self.supercell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.redraw()

        style = StructureInteractorStyle(self)
        self.iren.SetInteractorStyle(style)
        self.ren.parent = self

    def rotate_view(self, axis_ind=0, angle=0):
        """
        Rotate the camera view.

        Args:
            axis_ind:
                Index of axis to rotate. Defaults to 0, i.e., a-axis.
            angle:
                Angle to rotate by. Defaults to 0.
        """
        camera = self.ren.GetActiveCamera()
        if axis_ind == 0:
            camera.Roll(angle)
        elif axis_ind == 1:
            camera.Azimuth(angle)
        else:
            camera.Pitch(angle)
        self.ren_win.Render()

    def write_image(self, filename="image.png", magnification=1,
                    image_format="png"):
        """
        Save render window to an image.

        Arguments:
            filename:
                filename to save to. Defaults to image.png.
            magnification:
                magnification. Use it to render high res images.
            image_format:
                choose between jpeg, png.  Png is the default.
        """
        render_large = vtk.vtkRenderLargeImage()
        render_large.SetInput(self.ren)
        if image_format == "jpeg":
            writer = vtk.vtkJPEGWriter()
            writer.SetQuality(80)
        else:
            writer = vtk.vtkPNGWriter()

        render_large.SetMagnification(magnification)
        writer.SetFileName(filename)

        writer.SetInputConnection(render_large.GetOutputPort())
        self.ren_win.Render()
        writer.Write()
        del render_large

    def redraw(self, reset_camera=False):
        """
        Redraw the render window.

        Args:
            reset_camera:
                Set to True to reset the camera to a pre-determined default for
                each structure.  Defaults to False.
        """
        self.ren.RemoveAllViewProps()
        self.picker = None
        self.add_picker_fixed()
        self.helptxt_mapper = vtk.vtkTextMapper()
        tprops = self.helptxt_mapper.GetTextProperty()
        tprops.SetFontSize(14)
        tprops.SetFontFamilyToTimes()
        tprops.SetColor(0, 0, 0)

        if self.structure is not None:
            self.set_structure(self.structure, reset_camera)

        self.ren_win.Render()

    def orthongonalize_structure(self):
        if self.structure is not None:
            self.set_structure(self.structure.copy(sanitize=True))
        self.ren_win.Render()

    def display_help(self):
        """
        Display the help for various keyboard shortcuts.
        """
        helptxt = ["h : Toggle help",
                   "A/a, B/b or C/c : Increase/decrease cell by one a,"
                   " b or c unit vector", "# : Toggle showing of polyhedrons",
                   "-: Toggle showing of bonds", "r : Reset camera direction",
                   "[/]: Decrease or increase poly_radii_tol_factor "
                   "by 0.05. Value = " + str(self.poly_radii_tol_factor),
                   "Up/Down: Rotate view along Up direction by 90 "
                   "clockwise/anticlockwise",
                   "Left/right: Rotate view along camera direction by "
                   "90 clockwise/anticlockwise", "s: Save view to image.png",
                   "o: Orthogonalize structure"]
        self.helptxt_mapper.SetInput("\n".join(helptxt))
        self.helptxt_actor.SetPosition(10, 10)
        self.helptxt_actor.VisibilityOn()

    def set_structure(self, structure, reset_camera=True):
        """
        Add a structure to the visualizer.

        Args:
            structure:
                structure to visualize
            reset_camera:
                Set to True to reset the camera to a default determined based
                on the structure.
        """
        self.ren.RemoveAllViewProps()

        has_lattice = hasattr(structure, "lattice")

        if has_lattice:
            s = Structure.from_sites(structure, to_unit_cell=True)
            s.make_supercell(self.supercell)
        else:
            s = structure

        inc_coords = []
        for site in s:
            self.add_site(site)
            inc_coords.append(site.coords)

        count = 0
        labels = ["a", "b", "c"]
        colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

        if has_lattice:
            matrix = s.lattice.matrix

        if self.show_unit_cell and has_lattice:
            #matrix = s.lattice.matrix
            self.add_text([0, 0, 0], "o")
            for vec in matrix:
                self.add_line((0, 0, 0), vec, colors[count])
                self.add_text(vec, labels[count], colors[count])
                count += 1
            for (vec1, vec2) in itertools.permutations(matrix, 2):
                self.add_line(vec1, vec1 + vec2)
            for (vec1, vec2, vec3) in itertools.permutations(matrix, 3):
                self.add_line(vec1 + vec2, vec1 + vec2 + vec3)

        if self.show_bonds or self.show_polyhedron:
            elements = sorted(s.composition.elements, key=lambda a: a.X)
            anion = elements[-1]
            def contains_anion(site):
                for sp in site.species_and_occu.keys():
                    if sp.symbol == anion.symbol:
                        return True
                return False

            anion_radius = anion.average_ionic_radius
            for site in s:
                exclude = False
                max_radius = 0
                color = np.array([0, 0, 0])
                for sp, occu in site.species_and_occu.items():
                    if sp.symbol in self.excluded_bonding_elements \
                            or sp == anion:
                        exclude = True
                        break
                    max_radius = max(max_radius, sp.average_ionic_radius)
                    color = color + \
                        occu * np.array(self.el_color_mapping.get(sp.symbol,
                                                                  [0, 0, 0]))

                if not exclude:
                    max_radius = (1 + self.poly_radii_tol_factor) * \
                        (max_radius + anion_radius)
                    nn = structure.get_neighbors(site, max_radius)
                    nn_sites = []
                    for nnsite, dist in nn:
                        if contains_anion(nnsite):
                            nn_sites.append(nnsite)
                            if not in_coord_list(inc_coords, nnsite.coords):
                                self.add_site(nnsite)
                    if self.show_bonds:
                        self.add_bonds(nn_sites, site)
                    if self.show_polyhedron:
                        color = [i / 255 for i in color]
                        self.add_polyhedron(nn_sites, site, color)

        if self.show_help:
            self.helptxt_actor = vtk.vtkActor2D()
            self.helptxt_actor.VisibilityOn()
            self.helptxt_actor.SetMapper(self.helptxt_mapper)
            self.ren.AddActor(self.helptxt_actor)
            self.display_help()

        camera = self.ren.GetActiveCamera()
        if reset_camera:
            if has_lattice:
                #Adjust the camera for best viewing
                lengths = s.lattice.abc
                pos = (matrix[1] + matrix[2]) * 0.5 + matrix[0] * max(lengths) / \
                    lengths[0] * 3.5
                camera.SetPosition(pos)
                camera.SetViewUp(matrix[2])
                camera.SetFocalPoint((matrix[0] + matrix[1] + matrix[2]) * 0.5)
            else:
                origin = s.center_of_mass
                max_site = max(s, key=lambda site: site.distance_from_point(origin))
                camera.SetPosition(origin + 5 * (max_site.coords - origin))
                camera.SetFocalPoint(s.center_of_mass)

        self.structure = structure
        self.title = s.composition.formula

    def zoom(self, factor):
        """
        Zoom the camera view by a factor.
        """
        camera = self.ren.GetActiveCamera()
        camera.Zoom(factor)
        self.ren_win.Render()

    def show(self):
        """
        Display the visualizer.
        """
        self.iren.Initialize()
        self.ren_win.SetSize(800, 800)
        self.ren_win.SetWindowName(self.title)
        self.ren_win.Render()
        self.iren.Start()

    def add_site(self, site):
        """
        Add a site to the render window. The site is displayed as a sphere, the
        color of which is determined based on the element. Partially occupied
        sites are displayed as a single element color, though the site info
        still shows the partial occupancy.

        Args:
            site:
                Site to add.
        """
        start_angle = 0
        radius = 0
        total_occu = 0

        for specie, occu in site.species_and_occu.items():
            radius += occu * (specie.ionic_radius
                              if isinstance(specie, Specie)
                              and specie.ionic_radius
                              else specie.average_ionic_radius)
            total_occu += occu

        def add_partial_sphere(start, end, specie):
            sphere = vtk.vtkSphereSource()
            sphere.SetCenter(site.coords)
            sphere.SetRadius(0.2 + 0.002 * radius)
            sphere.SetThetaResolution(18)
            sphere.SetPhiResolution(18)
            sphere.SetStartTheta(start)
            sphere.SetEndTheta(end)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphere.GetOutputPort())
            self.mapper_map[mapper] = [site]
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            color = (0, 0, 0)
            if not specie:
                color = (1, 1, 1)
            elif specie.symbol in self.el_color_mapping:
                color = [i / 255 for i in self.el_color_mapping[specie.symbol]]
            actor.GetProperty().SetColor(color)
            self.ren.AddActor(actor)

        for specie, occu in site.species_and_occu.items():
            add_partial_sphere(start_angle, start_angle + 360 * occu, specie)
            start_angle += 360 * occu

        if total_occu < 1:
            add_partial_sphere(start_angle, start_angle
                               + 360 * (1 - total_occu), None)

    def add_text(self, coords, text, color=(0, 0, 0)):
        """
        Add text at a coordinate.

        Args:
            coords:
                Coordinates to add text at.
            text:
                Text to place.
            color:
                Color for text as RGB. Defaults to black.
        """
        source = vtk.vtkVectorText()
        source.SetText(text)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        follower = vtk.vtkFollower()
        follower.SetMapper(mapper)
        follower.GetProperty().SetColor(color)
        follower.SetPosition(coords)
        follower.SetScale(0.5)
        self.ren.AddActor(follower)
        follower.SetCamera(self.ren.GetActiveCamera())

    def add_line(self, start, end, color=(0.5, 0.5, 0.5), width=1):
        """
        Adds a line.

        Args:
            start:
                Starting coordinates for line.
            end:
                Ending coordinates for line.
            color:
                Color for text as RGB. Defaults to grey.
            width:
                Width of line. Defaults to 1.
        """
        source = vtk.vtkLineSource()
        source.SetPoint1(start)
        source.SetPoint2(end)

        vertexIDs = vtk.vtkStringArray()
        vertexIDs.SetNumberOfComponents(1)
        vertexIDs.SetName("VertexIDs")
        # Set the vertex labels
        vertexIDs.InsertNextValue("a")
        vertexIDs.InsertNextValue("b")
        source.GetOutput().GetPointData().AddArray(vertexIDs)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetLineWidth(width)
        self.ren.AddActor(actor)

    def add_polyhedron(self, neighbors, center, color, opacity=0.4,
                       draw_edges=False, edges_color=[0.0, 0.0, 0.0],
                       edges_linewidth=2):
        """
        Adds a polyhedron.

        Args:
            neighbors:
                Neighbors of the polyhedron (the vertices).
            center:
                The atom in the center of the polyhedron.
            color:
                Color for text as RGB.
            opacity:
                Opacity of the polyhedron
            draw_edges:
                If set to True, the a line will be drawn at each edge
            edges_color:
                Color of the line for the edges
            edges_linewidth:
                Width of the line drawn for the edges
        """
        points = vtk.vtkPoints()
        conv = vtk.vtkConvexPointSet()
        for i in range(len(neighbors)):
            x, y, z = neighbors[i].coords
            points.InsertPoint(i, x, y, z)
            conv.GetPointIds().InsertId(i, i)
        grid = vtk.vtkUnstructuredGrid()
        grid.Allocate(1, 1)
        grid.InsertNextCell(conv.GetCellType(), conv.GetPointIds())
        grid.SetPoints(points)

        dsm = vtk.vtkDataSetMapper()
        polysites = [center]
        polysites.extend(neighbors)
        self.mapper_map[dsm] = polysites
        if vtk.VTK_MAJOR_VERSION <= 5:
            dsm.SetInputConnection(grid.GetProducerPort())
        else:
            dsm.SetInputData(grid)
        ac = vtk.vtkActor()
        #ac.SetMapper(mapHull)
        ac.SetMapper(dsm)
        ac.GetProperty().SetOpacity(opacity)
        if color == 'element':
            # If partial occupations are involved, the color of the specie with
            # the highest occupation is used
            myoccu = 0.0
            for specie, occu in center.species_and_occu.items():
                if occu > myoccu:
                    myspecie = specie
                    myoccu = occu
            color = [i / 255 for i in self.el_color_mapping[myspecie.symbol]]
            ac.GetProperty().SetColor(color)
        else:
            ac.GetProperty().SetColor(color)
        if draw_edges:
            ac.GetProperty().SetEdgeColor(edges_color)
            ac.GetProperty().SetLineWidth(edges_linewidth)
            ac.GetProperty().EdgeVisibilityOn()
        self.ren.AddActor(ac)

    def add_triangle(self, neighbors, color, center=None, opacity=0.4,
                     draw_edges=False, edges_color=[0.0, 0.0, 0.0],
                     edges_linewidth=2):
        """
        Adds a triangular surface between three atoms.

        Args:
            atoms:
                Atoms between which a triangle will be drawn.
            color:
                Color for triangle as RGB.
            center:
                The "central atom" of the triangle
            opacity:
                opacity of the triangle
            draw_edges:
                If set to True, the a line will be drawn at each edge
            edges_color:
                Color of the line for the edges
            edges_linewidth:
                Width of the line drawn for the edges
        """
        points = vtk.vtkPoints()
        triangle = vtk.vtkTriangle()
        for ii in range(3):
            points.InsertNextPoint(neighbors[ii].x, neighbors[ii].y,
                                   neighbors[ii].z)
            triangle.GetPointIds().SetId(ii,ii)
        triangles = vtk.vtkCellArray()
        triangles.InsertNextCell(triangle)

        # polydata object
        trianglePolyData = vtk.vtkPolyData()
        trianglePolyData.SetPoints( points )
        trianglePolyData.SetPolys( triangles )

        # mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(trianglePolyData)

        ac = vtk.vtkActor()
        ac.SetMapper(mapper)
        ac.GetProperty().SetOpacity(opacity)
        if color == 'element':
            if center is None:
                raise ValueError(
                    'Color should be chosen according to the central atom, '
                    'and central atom is not provided')
            # If partial occupations are involved, the color of the specie with
            # the highest occupation is used
            myoccu = 0.0
            for specie, occu in center.species_and_occu.items():
                if occu > myoccu:
                    myspecie = specie
                    myoccu = occu
            color = [i / 255 for i in self.el_color_mapping[myspecie.symbol]]
            ac.GetProperty().SetColor(color)
        else:
            ac.GetProperty().SetColor(color)
        if draw_edges:
            ac.GetProperty().SetEdgeColor(edges_color)
            ac.GetProperty().SetLineWidth(edges_linewidth)
            ac.GetProperty().EdgeVisibilityOn()
        self.ren.AddActor(ac)

    def add_faces(self, faces, color, opacity=0.35):
        for face in faces:
            if len(face) == 3:
                points = vtk.vtkPoints()
                triangle = vtk.vtkTriangle()
                for ii in range(3):
                    points.InsertNextPoint(face[ii][0], face[ii][1], face[ii][2])
                    triangle.GetPointIds().SetId(ii, ii)
                triangles = vtk.vtkCellArray()
                triangles.InsertNextCell(triangle)
                trianglePolyData = vtk.vtkPolyData()
                trianglePolyData.SetPoints(points)
                trianglePolyData.SetPolys(triangles)
                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInput(trianglePolyData)
                ac = vtk.vtkActor()
                ac.SetMapper(mapper)
                ac.GetProperty().SetOpacity(opacity)
                ac.GetProperty().SetColor(color)
                self.ren.AddActor(ac)
            elif False and len(face) == 4:
                points = vtk.vtkPoints()
                for ii in range(4):
                    points.InsertNextPoint(face[ii][0], face[ii][1], face[ii][2])
                line1 = vtk.vtkLine()
                line1.GetPointIds().SetId(0, 0)
                line1.GetPointIds().SetId(1, 2)
                line2 = vtk.vtkLine()
                line2.GetPointIds().SetId(0, 3)
                line2.GetPointIds().SetId(1, 1)
                lines = vtk.vtkCellArray()
                lines.InsertNextCell(line1)
                lines.InsertNextCell(line2)
                polydata = vtk.vtkPolyData()
                polydata.SetPoints(points)
                polydata.SetLines(lines)
                ruledSurfaceFilter = vtk.vtkRuledSurfaceFilter()
                ruledSurfaceFilter.SetInput(polydata)
                ruledSurfaceFilter.SetResolution(15, 15)
                ruledSurfaceFilter.SetRuledModeToResample()
                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInput(ruledSurfaceFilter.GetOutput())
                ac = vtk.vtkActor()
                ac.SetMapper(mapper)
                ac.GetProperty().SetOpacity(opacity)
                ac.GetProperty().SetColor(color)
                self.ren.AddActor(ac)
            elif len(face) > 3:
                center = np.zeros(3, np.float)
                for site in face:
                    center += site
                center /= np.float(len(face))
                for ii in range(len(face)):
                    points = vtk.vtkPoints()
                    triangle = vtk.vtkTriangle()
                    points.InsertNextPoint(face[ii][0], face[ii][1], face[ii][2])
                    ii2 = np.mod(ii+1, len(face))
                    points.InsertNextPoint(face[ii2][0], face[ii2][1], face[ii2][2])
                    points.InsertNextPoint(center[0], center[1], center[2])
                    for ii in range(3):
                        triangle.GetPointIds().SetId(ii, ii)
                    triangles = vtk.vtkCellArray()
                    triangles.InsertNextCell(triangle)
                    trianglePolyData = vtk.vtkPolyData()
                    trianglePolyData.SetPoints(points)
                    trianglePolyData.SetPolys(triangles)
                    mapper = vtk.vtkPolyDataMapper()
                    mapper.SetInput(trianglePolyData)
                    ac = vtk.vtkActor()
                    ac.SetMapper(mapper)
                    ac.GetProperty().SetOpacity(opacity)
                    ac.GetProperty().SetColor(color)
                    self.ren.AddActor(ac)
            else:
                raise ValueError("Number of points for a face should be >= 3")

    def add_edges(self, edges, type='line', linewidth=2, color=[0.0, 0.0, 0.0]):
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        for iedge, edge in enumerate(edges):
            points.InsertPoint(2*iedge, edge[0])
            points.InsertPoint(2*iedge + 1, edge[1])
            lines.InsertNextCell(2)
            lines.InsertCellPoint(2*iedge)
            lines.InsertCellPoint(2*iedge + 1)
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(lines)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(polydata)
        ac = vtk.vtkActor()
        ac.SetMapper(mapper)
        ac.GetProperty().SetColor(color)
        ac.GetProperty().SetLineWidth(linewidth)
        self.ren.AddActor(ac)

    def add_bonds(self, neighbors, center, color=None, opacity=None,
                  radius=0.1):
        """
        Adds bonds for a site.

        Args:
            neighbors:
                Neighbors of the site.
            center:
                The site in the center for all bonds.
            color:
                Color of the tubes representing the bonds
            opacity:
                Opacity of the tubes representing the bonds
            radius:
                radius of tube s representing the bonds
        """
        points = vtk.vtkPoints()
        points.InsertPoint(0, center.x, center.y, center.z)
        n = len(neighbors)
        lines = vtk.vtkCellArray()
        for i in range(n):
            points.InsertPoint(i + 1, neighbors[i].coords)
            lines.InsertNextCell(2)
            lines.InsertCellPoint(0)
            lines.InsertCellPoint(i + 1)
        pd = vtk.vtkPolyData()
        pd.SetPoints(points)
        pd.SetLines(lines)

        tube = vtk.vtkTubeFilter()
        if vtk.VTK_MAJOR_VERSION <= 5:
            tube.SetInputConnection(pd.GetProducerPort())
        else:
            tube.SetInputData(pd)
        tube.SetRadius(radius)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(tube.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        if opacity is not None:
            actor.GetProperty().SetOpacity(opacity)
        if color is not None:
            actor.GetProperty().SetColor(color)
        self.ren.AddActor(actor)

    def add_picker_fixed(self):
        # Create a cell picker.
        picker = vtk.vtkCellPicker()
        # Create a Python function to create the text for the text mapper used
        # to display the results of picking.

        def annotate_pick(obj, event):
            if picker.GetCellId() < 0 and not self.show_help:
                self.helptxt_actor.VisibilityOff()
            else:
                mapper = picker.GetMapper()
                if mapper in self.mapper_map:
                    output = []
                    for site in self.mapper_map[mapper]:
                        row = ["{} - ".format(site.species_string),
                               ", ".join(["{:.3f}".format(c)
                                          for c in site.frac_coords]),
                               "[" + ", ".join(["{:.3f}".format(c)
                                                for c in site.coords]) +
                               "]"]
                        output.append("".join(row))
                    self.helptxt_mapper.SetInput("\n".join(output))
                    self.helptxt_actor.SetPosition(10, 10)
                    self.helptxt_actor.VisibilityOn()
                    self.show_help = False
        self.picker = picker
        picker.AddObserver("EndPickEvent", annotate_pick)
        self.iren.SetPicker(picker)

    def add_picker(self):
        # Create a cell picker.
        picker = vtk.vtkCellPicker()
        # Create a Python function to create the text for the text mapper used
        # to display the results of picking.
        source = vtk.vtkVectorText()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        follower = vtk.vtkFollower()
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
                    site = self.mapper_map[mapper]
                    output = [site.species_string, "Frac. coords: " +
                                                   " ".join(["{:.4f}".format(c)
                                                             for c in
                                                             site.frac_coords])]
                    source.SetText("\n".join(output))
                    follower.SetPosition(pick_pos)
                    follower.VisibilityOn()
        picker.AddObserver("EndPickEvent", annotate_pick)
        self.picker = picker
        self.iren.SetPicker(picker)


class StructureInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
    """
    A custom interactor style for visualizing structures.
    """

    def __init__(self, parent):
        self.parent = parent
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
        self.AddObserver("MouseMoveEvent", self.mouseMoveEvent)
        self.AddObserver("LeftButtonReleaseEvent", self.leftButtonReleaseEvent)
        self.AddObserver("KeyPressEvent", self.keyPressEvent)

    def leftButtonPressEvent(self, obj, event):
        self.mouse_motion = 0
        self.OnLeftButtonDown()
        return

    def mouseMoveEvent(self, obj, event):
        self.mouse_motion = 1
        self.OnMouseMove()
        return

    def leftButtonReleaseEvent(self, obj, event):
        ren = obj.GetCurrentRenderer()
        iren = ren.GetRenderWindow().GetInteractor()
        if self.mouse_motion == 0:
            pos = iren.GetEventPosition()
            iren.GetPicker().Pick(pos[0], pos[1], 0, ren)
        self.OnLeftButtonUp()
        return

    def keyPressEvent(self, obj, event):
        parent = obj.GetCurrentRenderer().parent
        sym = parent.iren.GetKeySym()

        if sym in "ABCabc":
            if sym == "A":
                parent.supercell[0][0] += 1
            elif sym == "B":
                parent.supercell[1][1] += 1
            elif sym == "C":
                parent.supercell[2][2] += 1
            elif sym == "a":
                parent.supercell[0][0] = max(parent.supercell[0][0] - 1, 1)
            elif sym == "b":
                parent.supercell[1][1] = max(parent.supercell[1][1] - 1, 1)
            elif sym == "c":
                parent.supercell[2][2] = max(parent.supercell[2][2] - 1, 1)
            parent.redraw()
        elif sym == "numbersign":
            parent.show_polyhedron = not parent.show_polyhedron
            parent.redraw()
        elif sym == "minus":
            parent.show_bonds = not parent.show_bonds
            parent.redraw()
        elif sym == "bracketleft":
            parent.poly_radii_tol_factor -= 0.05 \
                if parent.poly_radii_tol_factor > 0 else 0
            parent.redraw()
        elif sym == "bracketright":
            parent.poly_radii_tol_factor += 0.05
            parent.redraw()
        elif sym == "h":
            parent.show_help = not parent.show_help
            parent.redraw()
        elif sym == "r":
            parent.redraw(True)
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
            parent.orthongonalize_structure()
            parent.redraw()

        self.OnKeyPress()


def make_movie(structures, output_filename="movie.mp4", zoom=1.0, fps=20,
               bitrate=10000, quality=5):
    """
    Generate a movie from a sequence of structures using vtk and ffmpeg.

    Args:
        structures:
            sequence of structures
        output_filename:
            filename for structure output. defaults to movie.mp4
        zoom:
            A zoom to be applied to the visulizer. Defaults to 1.0
        fps:
            Frames per second for the movie. Defaults to 20.
        bitrate:
            Video bitate.  Defaults to 10000 (fairly high quality).
        quality:
            A quality scale. Defaults to 5.

    """
    vis = StructureVis()
    vis.show_help = False
    vis.redraw()
    vis.zoom(zoom)
    sigfig = int(math.floor(math.log10(len(structures))) + 1)
    filename = "image{0:0" + str(sigfig) + "d}.png"
    for i, s in enumerate(structures):
        vis.set_structure(s)
        vis.write_image(filename.format(i), 3)
    filename = "image%0" + str(sigfig) + "d.png"
    args = ["ffmpeg", "-y", "-qscale", str(quality), "-r", str(fps), "-b",
            str(bitrate), "-i", filename, output_filename]
    subprocess.Popen(args)
