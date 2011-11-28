#!/usr/bin/env python

'''
This module contains classes to wrap Python VTK to make nice molecular plots.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 27, 2011"

import os
import ConfigParser
import itertools

import vtk

from pymatgen.util.coord_utils import in_coord_list
from pymatgen.core.periodic_table import Specie

class StructureVis(object):
    """
    Provides Structure object visualization using VTK.
    """
    
    def __init__(self, element_color_mapping = None):
        """
        Arguments:
            element_color_mapping:
                Optional color mapping for the elements, as a dict of {symbol: rgb tuple}
                For example, {"Fe": (255,123,0), ....}
                If None is specified, a default based on Jmol's color scheme is used.
        """
        # create a rendering window and renderer
        self.ren = vtk.vtkRenderer()
        self.ren_win = vtk.vtkRenderWindow()
        self.ren_win.AddRenderer(self.ren)
        self.ren.SetBackground(1,1,1)
        self.title = "Structure Visualizer"
        # create a renderwindowinteractor
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.ren_win)
        
        if element_color_mapping:
            self.el_color_mapping = element_color_mapping
        else:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            self._config = ConfigParser.SafeConfigParser()
            self._config.optionxform = str
            self._config.readfp(open(os.path.join(module_dir, "ElementColorSchemes.cfg")))
            
            self.el_color_mapping = {}
            for (el, color) in self._config.items('VESTA'):
                self.el_color_mapping[el] = [int(i) for i in color.split(",")]

    def set_structure(self, structure, show_unit_cell = True, show_polyhedron = True, poly_radii_tol_factor = 0.2, excluded_bonding_elements = []):
        """
        Add a structure to the visualizer.
        
        Arguments:
            structure:
                structure to visualize
            show_unit_cell:
                Set to False to not show the unit cell boundaries. Defaults to True.
            show_polyhedron:
                Set to False to not show polyhedrons. Defaults to True.
            poly_radii_tol_factor:
                The polyhedron and bonding code uses the ionic radii of the elements or species to determine if two atoms are bonded.
                This specifies a tolerance scaling factor such that atoms which are 
                (1 + poly_radii_tol_factor) * sum of ionic radii apart are still considered as bonded.
            excluded_bonding_elements:
                List of atom types to exclude from bonding determination. Defaults to an empty list. Useful when trying to visualize 
                a certain atom type in the framework (e.g., Li in a Li-ion battery cathode material).
        """
        
        inc_coords = []
        for site in structure:
            self.add_atom(site.coords, site.species_and_occu.keys()[0])
            inc_coords.append(site.coords)
        
        count = 0
        labels = ['a','b','c']
        colors = [(1,0,0), (0,1,0), (0,0,1)]
        if show_unit_cell:
            self.add_text([0,0,0],"o")
            for vec in structure.lattice.matrix:
                self.add_line((0,0,0), vec, colors[count])
                self.add_text(vec,labels[count], colors[count])
                count += 1
            for (vec1, vec2) in itertools.permutations(structure.lattice.matrix, 2):
                self.add_line(vec1, vec1+vec2)
            for (vec1, vec2, vec3) in itertools.permutations(structure.lattice.matrix, 3):
                self.add_line(vec1 + vec2, vec1 + vec2 + vec3)
        
        if show_polyhedron:
            elements = sorted(structure.composition.elements)
            anion = elements[-1]
            anion_radius = anion.average_ionic_radius
            for site in structure:
                sp = site.species_and_occu.keys()[0]
                if sp != anion and sp.symbol not in excluded_bonding_elements:
                    max_radius = (1 + poly_radii_tol_factor) * (sp.average_ionic_radius + anion_radius) / 100
                    nn = structure.get_neighbors(site, max_radius)
                    nn_coords = []
                    for nnsite, dist in nn:
                        nn_coords.append(nnsite.coords)
                        if not in_coord_list(inc_coords, nnsite.coords):
                            self.add_atom(nnsite.coords, nnsite.species_and_occu.keys()[0])
                    self.add_bonds(nn_coords, site.x, site.y, site.z)
                    color = [i/255 for i in self.el_color_mapping.get(sp.symbol, [0,0,0])]
                    self.add_polyhedron(nn_coords, color)
                    
        #Adjust the camera for best viewing
        lattice = structure.lattice
        matrix = lattice.matrix
        lengths = lattice.abc
        pos = (matrix[0] + matrix[1]) * 0.5 + matrix[2] * max(lengths[0], lengths[1]) / lengths[2] * 3
        camera = self.ren.GetActiveCamera()
        camera.SetPosition(pos)           
        camera.SetFocalPoint((matrix[0] + matrix[1] + matrix[2]) * 0.5)      
        
        self.title = structure.composition.formula

    def show(self):
        # enable user interface interactor
        self.iren.Initialize()
        self.ren_win.SetSize(800,800)
        self.ren_win.SetWindowName(self.title)
        self.ren_win.Render()
        self.iren.Start()
        
    def add_atom(self, coords, specie):
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(coords)
        radius = specie.ionic_radius if isinstance(specie, Specie) else specie.average_ionic_radius
        sphere.SetRadius(0.2 + 0.002 * radius)
        sphere.SetThetaResolution(18)
        sphere.SetPhiResolution(18)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(sphere.GetOutput())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        color = (0,0,0)
        if specie.symbol in self.el_color_mapping:
            color = [i/255 for i in self.el_color_mapping[specie.symbol]]
        actor.GetProperty().SetColor(color)
        self.ren.AddActor(actor)

    def add_text(self, coords, text, color = (0,0,0)):
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

    def add_line(self, start, end, color = (0.5,0.5,0.5)):
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
        mapper.SetInput(source.GetOutput())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)
        self.ren.AddActor(actor)

    def add_polyhedron(self, neighbours, color):
        points = vtk.vtkPoints()
        conv = vtk.vtkConvexPointSet()
        for i in range(len(neighbours)):
            x,y,z = neighbours[i]
            points.InsertPoint(i,x,y,z)
            conv.GetPointIds().InsertId(i,i)
        grid=vtk.vtkUnstructuredGrid()
        grid.Allocate(1,1)
        grid.InsertNextCell(conv.GetCellType(),conv.GetPointIds())
        grid.SetPoints(points)

        dsm=vtk.vtkDataSetMapper()
        dsm.SetInput(grid)
        ac=vtk.vtkActor()
        #ac.SetMapper(mapHull)
        ac.SetMapper(dsm)
        ac.GetProperty().SetOpacity(0.4)
        ac.GetProperty().SetColor(color)
        self.ren.AddActor(ac)

    def add_bonds(self, neighbours, x, y, z):
        points=vtk.vtkPoints()
        points.InsertPoint(0, x, y, z)
        n = len(neighbours)
        lines = vtk.vtkCellArray()
        for i in range(n):
            points.InsertPoint(i+1,neighbours[i])
            lines.InsertNextCell(2)
            lines.InsertCellPoint(0)
            lines.InsertCellPoint(i+1)
        pd=vtk.vtkPolyData()
        pd.SetPoints(points)
        pd.SetLines(lines)

        tube = vtk.vtkTubeFilter()
        tube.SetInput(pd)
        tube.SetRadius(0.1)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(tube.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        self.ren.AddActor(actor)      

if __name__ == "__main__":
    #To be moved to proper unittest in future version
    from pymatgen.io.vaspio import Poscar, Vasprun
    filepath = os.path.join('..', 'io','tests', 'vasp_testfiles','POSCAR')
    poscar = Poscar.from_file(filepath)
    s = poscar.struct
    vis = StructureVis()
    vis.set_structure(s, excluded_bonding_elements=["Li"], poly_radii_tol_factor = 0.2)
    vis.show()