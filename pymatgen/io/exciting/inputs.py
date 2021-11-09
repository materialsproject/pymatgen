# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Classes for reading/manipulating/writing exciting input files.
"""

import xml.etree.cElementTree as ET

import numpy as np
import scipy.constants as const
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

__author__ = "Christian Vorwerk"
__copyright__ = "Copyright 2016"
__version__ = "1.0"
__maintainer__ = "Christian Vorwerk"
__email__ = "vorwerk@physik.hu-berlin.de"
__status__ = "Development"
__date__ = "Nov 28, 2016"


class ExcitingInput(MSONable):
    """
    Object for representing the data stored in the structure part of the
    exciting input.

    .. attribute:: structure

        Associated Structure.

    .. attribute:: title

        Optional title string.

    .. attribute:: lockxyz

        Lockxyz attribute for each site if available. A Nx3 array of
        booleans.
    """

    def __init__(self, structure, title=None, lockxyz=None):
        """
        Args:
            structure (Structure):  Structure object.
            title (str): Optional title for exciting input. Defaults to unit
                cell formula of structure. Defaults to None.
            lockxyz (Nx3 array): bool values for selective dynamics,
                where N is number of sites. Defaults to None.
        """

        if structure.is_ordered:
            site_properties = {}
            if lockxyz:
                site_properties["selective_dynamics"] = lockxyz
            self.structure = structure.copy(site_properties=site_properties)
            self.title = structure.formula if title is None else title
        else:
            raise ValueError("Structure with partial occupancies cannot be converted into exciting input!")

    # define conversion factor between Bohr radius and Angstrom
    bohr2ang = const.value("Bohr radius") / const.value("Angstrom star")

    @property
    def lockxyz(self):
        """
        :return: Selective dynamics site properties.
        """
        return self.structure.site_properties.get("selective_dynamics")

    @lockxyz.setter
    def lockxyz(self, lockxyz):
        self.structure.add_site_property("selective_dynamics", lockxyz)

    @staticmethod
    def from_string(data):
        """
        Reads the exciting input from a string
        """

        root = ET.fromstring(data)
        speciesnode = root.find("structure").iter("species")
        elements = []
        positions = []
        vectors = []
        lockxyz = []
        # get title
        title_in = str(root.find("title").text)
        # Read elements and coordinates
        for nodes in speciesnode:
            symbol = nodes.get("speciesfile").split(".")[0]
            if len(symbol.split("_")) == 2:
                symbol = symbol.split("_")[0]
            if Element.is_valid_symbol(symbol):
                # Try to recognize the element symbol
                element = symbol
            else:
                raise ValueError("Unknown element!")
            for atom in nodes.iter("atom"):
                x, y, z = atom.get("coord").split()
                positions.append([float(x), float(y), float(z)])
                elements.append(element)
                # Obtain lockxyz for each atom
                if atom.get("lockxyz") is not None:
                    lxyz = []
                    for l in atom.get("lockxyz").split():
                        if l in ("True", "true"):
                            lxyz.append(True)
                        else:
                            lxyz.append(False)
                    lockxyz.append(lxyz)
                else:
                    lockxyz.append([False, False, False])
        # check the atomic positions type
        if "cartesian" in root.find("structure").attrib.keys():
            if root.find("structure").attrib["cartesian"]:
                cartesian = True
                for i, p in enumerate(positions):
                    for j in range(3):
                        p[j] = p[j] * ExcitingInput.bohr2ang
                print(positions)
        else:
            cartesian = False
        # get the scale attribute
        scale_in = root.find("structure").find("crystal").get("scale")
        if scale_in:
            scale = float(scale_in) * ExcitingInput.bohr2ang
        else:
            scale = ExcitingInput.bohr2ang
        # get the stretch attribute
        stretch_in = root.find("structure").find("crystal").get("stretch")
        if stretch_in:
            stretch = np.array([float(a) for a in stretch_in])
        else:
            stretch = np.array([1.0, 1.0, 1.0])
        # get basis vectors and scale them accordingly
        basisnode = root.find("structure").find("crystal").iter("basevect")
        for vect in basisnode:
            x, y, z = vect.text.split()
            vectors.append(
                [
                    float(x) * stretch[0] * scale,
                    float(y) * stretch[1] * scale,
                    float(z) * stretch[2] * scale,
                ]
            )
        # create lattice and structure object
        lattice_in = Lattice(vectors)
        structure_in = Structure(lattice_in, elements, positions, coords_are_cartesian=cartesian)

        return ExcitingInput(structure_in, title_in, lockxyz)

    @staticmethod
    def from_file(filename):
        """
        :param filename: Filename
        :return: ExcitingInput
        """
        with zopen(filename, "rt") as f:
            data = f.read().replace("\n", "")
        return ExcitingInput.from_string(data)

    def write_etree(self, celltype, cartesian=False, bandstr=False, symprec=0.4, angle_tolerance=5, **kwargs):
        """
        Writes the exciting input parameters to an xml object.

        Args:
            celltype (str): Choice of unit cell. Can be either the unit cell
            from self.structure ("unchanged"), the conventional cell
            ("conventional"), or the primitive unit cell ("primitive").

            cartesian (bool): Whether the atomic positions are provided in
            Cartesian or unit-cell coordinates. Default is False.

            bandstr (bool): Whether the bandstructure path along the
            HighSymmKpath is included in the input file. Only supported if the
            celltype is set to "primitive". Default is False.

            symprec (float): Tolerance for the symmetry finding. Default is 0.4.

            angle_tolerance (float): Angle tolerance for the symmetry finding.
            Default is 5.

            **kwargs: Additional parameters for the input file.

        Returns:
            ET.Element containing the input XML structure
        """
        root = ET.Element("input")
        root.set(
            "{http://www.w3.org/2001/XMLSchema-instance}noNamespaceSchemaLocation",
            "http://xml.exciting-code.org/excitinginput.xsd",
        )
        title = ET.SubElement(root, "title")
        title.text = self.title
        if cartesian:
            structure = ET.SubElement(root, "structure", cartesian="true", speciespath="./")
        else:
            structure = ET.SubElement(root, "structure", speciespath="./")

        crystal = ET.SubElement(structure, "crystal")
        # set scale such that lattice vector can be given in Angstrom
        ang2bohr = const.value("Angstrom star") / const.value("Bohr radius")
        crystal.set("scale", str(ang2bohr))
        # determine which structure to use
        finder = SpacegroupAnalyzer(self.structure, symprec=symprec, angle_tolerance=angle_tolerance)
        if celltype == "primitive":
            new_struct = finder.get_primitive_standard_structure(international_monoclinic=False)
        elif celltype == "conventional":
            new_struct = finder.get_conventional_standard_structure(international_monoclinic=False)
        elif celltype == "unchanged":
            new_struct = self.structure
        else:
            raise ValueError("Type of unit cell not recognized!")

        # write lattice
        basis = new_struct.lattice.matrix
        for i in range(3):
            basevect = ET.SubElement(crystal, "basevect")
            basevect.text = "%16.8f %16.8f %16.8f" % (
                basis[i][0],
                basis[i][1],
                basis[i][2],
            )
        # write atomic positions for each species
        index = 0
        for i in sorted(new_struct.types_of_species, key=lambda el: el.X):
            species = ET.SubElement(structure, "species", speciesfile=i.symbol + ".xml")
            sites = new_struct.indices_from_symbol(i.symbol)

            for j in sites:
                coord = "%16.8f %16.8f %16.8f" % (
                    new_struct[j].frac_coords[0],
                    new_struct[j].frac_coords[1],
                    new_struct[j].frac_coords[2],
                )
                # obtain cartesian coords from fractional ones if needed
                if cartesian:
                    coord2 = []
                    for k in range(3):
                        inter = (
                            new_struct[j].frac_coords[k] * basis[0][k]
                            + new_struct[j].frac_coords[k] * basis[1][k]
                            + new_struct[j].frac_coords[k] * basis[2][k]
                        ) * ang2bohr
                        coord2.append(inter)
                    coord = "%16.8f %16.8f %16.8f" % (coord2[0], coord2[1], coord2[2])

                # write atomic positions
                index = index + 1
                _ = ET.SubElement(species, "atom", coord=coord)
        # write bandstructure if needed
        if bandstr and celltype == "primitive":
            kpath = HighSymmKpath(new_struct, symprec=symprec, angle_tolerance=angle_tolerance)
            prop = ET.SubElement(root, "properties")
            bandstrct = ET.SubElement(prop, "bandstructure")
            for i in range(len(kpath.kpath["path"])):
                plot = ET.SubElement(bandstrct, "plot1d")
                path = ET.SubElement(plot, "path", steps="100")
                for j in range(len(kpath.kpath["path"][i])):
                    symbol = kpath.kpath["path"][i][j]
                    coords = kpath.kpath["kpoints"][symbol]
                    coord = "%16.8f %16.8f %16.8f" % (coords[0], coords[1], coords[2])
                    if symbol == "\\Gamma":
                        symbol = "GAMMA"
                    _ = ET.SubElement(path, "point", coord=coord, label=symbol)
        elif bandstr and celltype != "primitive":
            raise ValueError(
                "Bandstructure is only implemented for the \
                              standard primitive unit cell!"
            )

        # write extra parameters from kwargs if provided
        self._dicttoxml(kwargs, root)

        return root

    def write_string(self, celltype, cartesian=False, bandstr=False, symprec=0.4, angle_tolerance=5, **kwargs):
        """
        Writes exciting input.xml as a string.

        Args:
            celltype (str): Choice of unit cell. Can be either the unit cell
            from self.structure ("unchanged"), the conventional cell
            ("conventional"), or the primitive unit cell ("primitive").

            cartesian (bool): Whether the atomic positions are provided in
            Cartesian or unit-cell coordinates. Default is False.

            bandstr (bool): Whether the bandstructure path along the
            HighSymmKpath is included in the input file. Only supported if the
            celltype is set to "primitive". Default is False.

            symprec (float): Tolerance for the symmetry finding. Default is 0.4.

            angle_tolerance (float): Angle tolerance for the symmetry finding.
            Default is 5.

            **kwargs: Additional parameters for the input file.

        Returns:
            String
        """
        try:
            root = self.write_etree(celltype, cartesian, bandstr, symprec, angle_tolerance, **kwargs)
            self._indent(root)
            # output should be a string not a bytes object
            string = ET.tostring(root).decode("UTF-8")
        except Exception:
            raise ValueError("Incorrect celltype!")
        return string

    def write_file(self, celltype, filename, cartesian=False, bandstr=False, symprec=0.4, angle_tolerance=5, **kwargs):
        """
        Writes exciting input file.

        Args:
            celltype (str): Choice of unit cell. Can be either the unit cell
            from self.structure ("unchanged"), the conventional cell
            ("conventional"), or the primitive unit cell ("primitive").

            filename (str): Filename for exciting input.

            cartesian (bool): Whether the atomic positions are provided in
            Cartesian or unit-cell coordinates. Default is False.

            bandstr (bool): Whether the bandstructure path along the
            HighSymmKpath is included in the input file. Only supported if the
            celltype is set to "primitive". Default is False.

            symprec (float): Tolerance for the symmetry finding. Default is 0.4.

            angle_tolerance (float): Angle tolerance for the symmetry finding.
            Default is 5.

            **kwargs: Additional parameters for the input file.
        """
        try:
            root = self.write_etree(celltype, cartesian, bandstr, symprec, angle_tolerance, **kwargs)
            self._indent(root)
            tree = ET.ElementTree(root)
            tree.write(filename)
        except Exception:
            raise ValueError("Incorrect celltype!")

    # Missing PrettyPrint option in the current version of xml.etree.cElementTree
    @staticmethod
    def _indent(elem, level=0):
        """
        Helper method to indent elements.

        :param elem:
        :param level:
        :return:
        """
        i = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for el in elem:
                ExcitingInput._indent(el, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i

    def _dicttoxml(self, paramdict_, element):
        for key, value in paramdict_.items():
            if isinstance(value, str) and key == "text()":
                element.text = value
            elif isinstance(value, str):
                element.attrib[key] = value
            elif isinstance(value, list):
                for item in value:
                    self._dicttoxml(item, ET.SubElement(element, key))
            elif isinstance(value, dict):
                if element.findall(key) == []:
                    self._dicttoxml(value, ET.SubElement(element, key))
                else:
                    self._dicttoxml(value, element.findall(key)[0])
            else:
                print("cannot deal with", key, "=", value)
