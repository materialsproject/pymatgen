"""
Module implementing classes and functions to use Zeo++
by Maciej Haranczyk.

If using this module, cite the following paper on Zeo++:
T.F. Willems, C.H. Rycroft, M. Kazi, J.C. Meza, and M. Haranczyk,
Algorithms and tools for high-throughput geometry-based analysis of crystalline porous materials,
Microporous and Mesoporous Materials, 149 (2012) 134-141.

Zeo++ Installation Steps:
========================
A stable version of Zeo++ can be obtained from http://zeoplusplus.org.
Instructions can be found at http://www.zeoplusplus.org/download.html

Zeo++ Post-Installation Checking:
==============================
1) Go to pymatgen/io/tests and run "python test_zeoio.py"
   If Zeo++ python bindings are properly installed, the tests should
   pass. One or two tests will be skipped.
b) Go to pymatgen/analysis/defects/tests and run
   "python test_point_defects.py". Lots of tests will be skipped if GULP
   is not installed. But there should be no errors.
"""

from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING

from monty.dev import requires
from monty.io import zopen
from monty.tempfile import ScratchDir

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cssr import Cssr
from pymatgen.io.xyz import XYZ

try:
    from zeo.cluster import prune_voronoi_network_close_node
    from zeo.netstorage import AtomNetwork

    zeo_found = True
except ImportError:
    zeo_found = False
    AtomNetwork = prune_voronoi_network_close_node = None

if TYPE_CHECKING:
    from pathlib import Path

    from typing_extensions import Self

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__data__ = "Aug 2, 2013"


class ZeoCssr(Cssr):
    """
    ZeoCssr adds extra fields to CSSR sites to conform with Zeo++
    input CSSR format. The coordinate system is rotated from xyz to zyx.
    This change aligns the pivot axis of pymatgen (z-axis) to pivot axis
    of Zeo++ (x-axis) for structural modifications.
    """

    def __init__(self, structure: Structure):
        """
        Args:
            structure: A structure to create ZeoCssr object.
        """
        super().__init__(structure)

    def __str__(self):
        """
        CSSR.__str__ method is modified to pad 0's to the CSSR site data.
        The padding is to conform with the CSSR format supported Zeo++.
        The oxidation state is stripped from site.specie
        Also coordinate system is rotated from xyz to zxy.
        """
        a, b, c = self.structure.lattice.lengths
        alpha, beta, gamma = self.structure.lattice.angles
        output = [
            f"{c:.4f} {a:.4f} {b:.4f}",
            f"{gamma:.2f} {alpha:.2f} {beta:.2f} SPGR =  1 P 1    OPT = 1",
            f"{len(self.structure)} 0",
            f"0 {self.structure.formula}",
        ]
        for idx, site in enumerate(self.structure):
            charge = getattr(site, "charge", 0)
            # specie = site.specie.symbol
            specie = site.species_string
            output.append(f"{idx + 1} {specie} {site.c:.4f} {site.a:.4f} {site.b:.4f} 0 0 0 0 0 0 0 0 {charge:.4f}")

        return "\n".join(output)

    @classmethod
    def from_str(cls, string: str) -> Self:
        """
        Reads a string representation to a ZeoCssr object.

        Args:
            string: A string representation of a ZeoCSSR.

        Returns:
            ZeoCssr object.
        """
        lines = string.split("\n")
        tokens = lines[0].split()
        lengths = [float(i) for i in tokens]
        tokens = lines[1].split()
        angles = [float(i) for i in tokens[0:3]]
        # Zeo++ takes x-axis along a and pymatgen takes z-axis along c
        a = lengths.pop(-1)
        lengths.insert(0, a)
        alpha = angles.pop(-1)
        angles.insert(0, alpha)
        lattice = Lattice.from_parameters(*lengths, *angles)

        sp = []
        coords = []
        charge = []
        for line in lines[4:]:
            match = re.match(
                r"\d+\s+(\w+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+(?:0\s+){8}([0-9\-\.]+)",
                line.strip(),
            )
            if match:
                sp.append(match.group(1))
                # coords.append([float(m.group(i)) for i in xrange(2, 5)])
                # Zeo++ takes x-axis along a and pymatgen takes z-axis along c
                coords.append([float(match.group(i)) for i in [3, 4, 2]])
                charge.append(match.group(5))
        return cls(Structure(lattice, sp, coords, site_properties={"charge": charge}))

    @classmethod
    def from_file(cls, filename: str | Path) -> Self:
        """
        Reads a CSSR file to a ZeoCssr object.

        Args:
            filename: Filename to read from.

        Returns:
            ZeoCssr object.
        """
        with zopen(filename, mode="r") as file:
            return cls.from_str(file.read())


class ZeoVoronoiXYZ(XYZ):
    """
    Class to read Voronoi Nodes from XYZ file written by Zeo++.
    The sites have an additional column representing the voronoi node radius.
    The voronoi node radius is represented by the site property voronoi_radius.
    """

    def __init__(self, mol):
        """
        Args:
            mol: Input molecule holding the voronoi node information.
        """
        super().__init__(mol)

    @classmethod
    def from_str(cls, contents: str) -> Self:
        """
        Creates Zeo++ Voronoi XYZ object from a string.
        from_string method of XYZ class is being redefined.

        Args:
            contents: String representing Zeo++ Voronoi XYZ file.

        Returns:
            ZeoVoronoiXYZ object
        """
        lines = contents.split("\n")
        num_sites = int(lines[0])
        coords = []
        sp = []
        prop = []
        coord_patt = re.compile(r"(\w+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)")
        for i in range(2, 2 + num_sites):
            m = coord_patt.search(lines[i])
            if m:
                sp.append(m.group(1))  # this is 1-indexed
                # coords.append(map(float, m.groups()[1:4]))  # this is 0-indexed
                coords.append([float(j) for j in [m.group(i) for i in [3, 4, 2]]])
                prop.append(float(m.group(5)))
        return cls(Molecule(sp, coords, site_properties={"voronoi_radius": prop}))

    @classmethod
    def from_file(cls, filename: str | Path) -> Self:
        """
        Creates XYZ object from a file.

        Args:
            filename: XYZ filename

        Returns:
            XYZ object
        """
        with zopen(filename) as file:
            return cls.from_str(file.read())

    def __str__(self) -> str:
        output = [str(len(self._mols[0])), self._mols[0].formula]
        prec = self.precision
        for site in self._mols[0]:
            x, y, z = site.coords
            symbol, voronoi_radius = site.specie.symbol, site.properties["voronoi_radius"]
            output.append(f"{symbol} {z:.{prec}f} {x:.{prec}f} {y:.{prec}f} {voronoi_radius:.{prec}f}")
        return "\n".join(output)


@requires(
    zeo_found,
    "get_voronoi_nodes requires Zeo++ cython extension to be "
    "installed. Please contact developers of Zeo++ to obtain it.",
)
def get_voronoi_nodes(structure, rad_dict=None, probe_rad=0.1):
    """
    Analyze the void space in the input structure using voronoi decomposition
    Calls Zeo++ for Voronoi decomposition.

    Args:
        structure: pymatgen Structure
        rad_dict (optional): Dictionary of radii of elements in structure.
            If not given, Zeo++ default values are used.
            Note: Zeo++ uses atomic radii of elements.
            For ionic structures, pass rad_dict with ionic radii
        probe_rad (optional): Sampling probe radius in Angstroms. Default is
            0.1 A

    Returns:
        voronoi nodes as pymatgen Structure within the unit cell defined by the lattice of
        input structure voronoi face centers as pymatgen Structure within the unit cell
        defined by the lattice of input structure
    """
    with ScratchDir("."):
        name = "temp_zeo1"
        zeo_inp_filename = f"{name}.cssr"
        ZeoCssr(structure).write_file(zeo_inp_filename)
        rad_file = None
        rad_flag = False

        if rad_dict:
            rad_file = f"{name}.rad"
            rad_flag = True
            with open(rad_file, "w+", encoding="utf-8") as file:
                for el in rad_dict:
                    file.write(f"{el} {rad_dict[el].real}\n")

        atom_net = AtomNetwork.read_from_CSSR(zeo_inp_filename, rad_flag=rad_flag, rad_file=rad_file)
        vor_net, vor_edge_centers, vor_face_centers = atom_net.perform_voronoi_decomposition()
        vor_net.analyze_writeto_XYZ(name, probe_rad, atom_net)
        voro_out_filename = f"{name}_voro.xyz"
        voro_node_mol = ZeoVoronoiXYZ.from_file(voro_out_filename).molecule

    species = ["X"] * len(voro_node_mol)
    coords = []
    prop = []
    for site in voro_node_mol:
        coords.append(list(site.coords))
        prop.append(site.properties["voronoi_radius"])

    lattice = Lattice.from_parameters(*structure.lattice.parameters)
    vor_node_struct = Structure(
        lattice,
        species,
        coords,
        coords_are_cartesian=True,
        to_unit_cell=True,
        site_properties={"voronoi_radius": prop},
    )

    # PMG-Zeo c<->a transformation for Voronoi face centers
    rot_face_centers = [(center[1], center[2], center[0]) for center in vor_face_centers]
    rot_edge_centers = [(center[1], center[2], center[0]) for center in vor_edge_centers]

    species = ["X"] * len(rot_face_centers)
    prop = [0.0] * len(rot_face_centers)  # Voronoi radius not evaluated for fc
    vor_facecenter_struct = Structure(
        lattice,
        species,
        rot_face_centers,
        coords_are_cartesian=True,
        to_unit_cell=True,
        site_properties={"voronoi_radius": prop},
    )

    species = ["X"] * len(rot_edge_centers)
    prop = [0.0] * len(rot_edge_centers)  # Voronoi radius not evaluated for fc
    vor_edgecenter_struct = Structure(
        lattice,
        species,
        rot_edge_centers,
        coords_are_cartesian=True,
        to_unit_cell=True,
        site_properties={"voronoi_radius": prop},
    )

    return vor_node_struct, vor_edgecenter_struct, vor_facecenter_struct


def get_high_accuracy_voronoi_nodes(structure, rad_dict, probe_rad=0.1):
    """
    Analyze the void space in the input structure using high accuracy
    voronoi decomposition.
    Calls Zeo++ for Voronoi decomposition.

    Args:
        structure: pymatgen Structure
        rad_dict (optional): Dictionary of radii of elements in structure.
            If not given, Zeo++ default values are used.
            Note: Zeo++ uses atomic radii of elements.
            For ionic structures, pass rad_dict with ionic radii
        probe_rad (optional): Sampling probe radius in Angstroms.
            Default is 0.1 A

    Returns:
        voronoi nodes as pymatgen Structure within the
        unit cell defined by the lattice of input structure
        voronoi face centers as pymatgen Structure within the
        unit cell defined by the lattice of input structure
    """
    with ScratchDir("."):
        name = "temp_zeo1"
        zeo_inp_filename = f"{name}.cssr"
        ZeoCssr(structure).write_file(zeo_inp_filename)
        rad_flag = True
        rad_file = f"{name}.rad"
        with open(rad_file, "w+", encoding="utf-8") as file:
            for el in rad_dict:
                print(f"{el} {rad_dict[el].real}", file=file)

        atom_net = AtomNetwork.read_from_CSSR(zeo_inp_filename, rad_flag=rad_flag, rad_file=rad_file)
        # vornet, vor_edge_centers, vor_face_centers = atom_net.perform_voronoi_decomposition()
        red_ha_vornet = prune_voronoi_network_close_node(atom_net)
        # generate_simplified_highaccuracy_voronoi_network(atom_net)
        # get_nearest_largest_diameter_highaccuracy_vornode(atom_net)
        red_ha_vornet.analyze_writeto_XYZ(name, probe_rad, atom_net)
        voro_out_filename = f"{name}_voro.xyz"
        voro_node_mol = ZeoVoronoiXYZ.from_file(voro_out_filename).molecule

    species = ["X"] * len(voro_node_mol)
    coords = []
    prop = []
    for site in voro_node_mol:
        coords.append(list(site.coords))
        prop.append(site.properties["voronoi_radius"])

    lattice = Lattice.from_parameters(*structure.lattice.parameters)
    return Structure(
        lattice,
        species,
        coords,
        coords_are_cartesian=True,
        to_unit_cell=True,
        site_properties={"voronoi_radius": prop},
    )


@requires(
    zeo_found,
    "get_voronoi_nodes requires Zeo++ cython extension to be "
    "installed. Please contact developers of Zeo++ to obtain it.",
)
def get_free_sphere_params(structure, rad_dict=None, probe_rad=0.1):
    """
    Analyze the void space in the input structure using voronoi decomposition
    Calls Zeo++ for Voronoi decomposition.

    Args:
        structure: pymatgen Structure
        rad_dict (optional): Dictionary of radii of elements in structure.
            If not given, Zeo++ default values are used.
            Note: Zeo++ uses atomic radii of elements.
            For ionic structures, pass rad_dict with ionic radii
        probe_rad (optional): Sampling probe radius in Angstroms. Default is
            0.1 A

    Returns:
        voronoi nodes as pymatgen Structure within the
        unit cell defined by the lattice of input structure
        voronoi face centers as pymatgen Structure within the
        unit cell defined by the lattice of input structure
    """
    with ScratchDir("."):
        name = "temp_zeo1"
        zeo_inp_filename = f"{name}.cssr"
        ZeoCssr(structure).write_file(zeo_inp_filename)
        rad_file = None
        rad_flag = False

        if rad_dict:
            rad_file = f"{name}.rad"
            rad_flag = True
            with open(rad_file, "w+", encoding="utf-8") as file:
                for el in rad_dict:
                    file.write(f"{el} {rad_dict[el].real}\n")

        atom_net = AtomNetwork.read_from_CSSR(zeo_inp_filename, rad_flag=rad_flag, rad_file=rad_file)
        out_file = "temp.res"
        atom_net.calculate_free_sphere_parameters(out_file)
        if os.path.isfile(out_file) and os.path.getsize(out_file) > 0:
            with open(out_file, encoding="utf-8") as file:
                output = file.readline()
        else:
            output = ""
    fields = [val.strip() for val in output.split()][1:4]
    if len(fields) == 3:
        fields = [float(field) for field in fields]
        return {
            "inc_sph_max_dia": fields[0],
            "free_sph_max_dia": fields[1],
            "inc_sph_along_free_sph_path_max_dia": fields[2],
        }
    return None
