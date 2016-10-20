# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import os
import re

from monty.io import zopen
from monty.dev import requires
from monty.tempfile import ScratchDir

from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.io.cssr import Cssr
from pymatgen.io.xyz import XYZ

try:
    from zeo.netstorage import AtomNetwork, VoronoiNetwork
    from zeo.area_volume import volume, surface_area
    from zeo.cluster import get_nearest_largest_diameter_highaccuracy_vornode, \
        generate_simplified_highaccuracy_voronoi_network, \
        prune_voronoi_network_close_node

    zeo_found = True
except ImportError:
    zeo_found = False

"""
Module implementing classes and functions to use Zeo++.
Zeo++ can be obtained from http://www.maciejharanczyk.info/Zeopp/
"""

__author__ = "Bharat Medasani"
__copyright = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Bharat Medasani"
__email__ = "bkmedasani@lbl.gov"
__data__ = "Aug 2, 2013"


class ZeoCssr(Cssr):
    """
    ZeoCssr adds extra fields to CSSR sites to conform with Zeo++
    input CSSR format. The coordinate system is rorated from xyz to zyx.
    This change aligns the pivot axis of pymatgen (z-axis) to pivot axis
    of Zeo++ (x-axis) for structurural modifications.

    Args:
        structure: A structure to create ZeoCssr object
    """

    def __init__(self, structure):
        super(ZeoCssr, self).__init__(structure)

    def __str__(self):
        """
        CSSR.__str__ method is modified to padd 0's to the CSSR site data.
        The padding is to conform with the CSSR format supported Zeo++.
        The oxidation state is stripped from site.specie
        Also coordinate system is rotated from xyz to zxy
        """
        output = [
            "{:.4f} {:.4f} {:.4f}"
                # .format(*self.structure.lattice.abc),
                .format(self.structure.lattice.c,
                        self.structure.lattice.a,
                        self.structure.lattice.b),
            "{:.2f} {:.2f} {:.2f} SPGR =  1 P 1    OPT = 1"
                # .format(*self.structure.lattice.angles),
                .format(self.structure.lattice.gamma,
                        self.structure.lattice.alpha,
                        self.structure.lattice.beta),
            "{} 0".format(len(self.structure)),
            "0 {}".format(self.structure.formula)
        ]
        for i, site in enumerate(self.structure.sites):
            # if not hasattr(site, 'charge'):
            #    charge = 0
            # else:
            #    charge = site.charge
            charge = site.charge if hasattr(site, 'charge') else 0
            # specie = site.specie.symbol
            specie = site.species_string
            output.append(
                "{} {} {:.4f} {:.4f} {:.4f} 0 0 0 0 0 0 0 0 {:.4f}"
                    .format(
                    i + 1, specie, site.c, site.a, site.b, charge
                    # i+1, site.specie, site.a, site.b, site.c, site.charge
                )
            )

        return "\n".join(output)

    @staticmethod
    def from_string(string):
        """
        Reads a string representation to a ZeoCssr object.

        Args:
            string: A string representation of a ZeoCSSR.

        Returns:
            ZeoCssr object.
        """
        lines = string.split("\n")
        toks = lines[0].split()
        lengths = [float(i) for i in toks]
        toks = lines[1].split()
        angles = [float(i) for i in toks[0:3]]
        # Zeo++ takes x-axis along a and pymatgen takes z-axis along c
        a = lengths.pop(-1)
        lengths.insert(0, a)
        alpha = angles.pop(-1)
        angles.insert(0, alpha)
        latt = Lattice.from_lengths_and_angles(lengths, angles)
        sp = []
        coords = []
        chrg = []
        for l in lines[4:]:
            m = re.match("\d+\s+(\w+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+" +
                         "([0-9\-\.]+)\s+(?:0\s+){8}([0-9\-\.]+)", l.strip())
            if m:
                sp.append(m.group(1))
                # coords.append([float(m.group(i)) for i in xrange(2, 5)])
                # Zeo++ takes x-axis along a and pymatgen takes z-axis along c
                coords.append([float(m.group(i)) for i in [3, 4, 2]])
                chrg.append(m.group(5))
        return ZeoCssr(
            Structure(latt, sp, coords, site_properties={'charge': chrg})
        )

    @staticmethod
    def from_file(filename):
        """
        Reads a CSSR file to a ZeoCssr object.

        Args:
            filename: Filename to read from.

        Returns:
            ZeoCssr object.
        """
        with zopen(filename, "r") as f:
            return ZeoCssr.from_string(f.read())


class ZeoVoronoiXYZ(XYZ):
    """
    Class to read Voronoi Nodes from XYZ file written by Zeo++.
    The sites have an additional column representing the voronoi node radius.
    The voronoi node radius is represented by the site property voronoi_radius.

    Args:
        mol: Input molecule holding the voronoi node information
    """

    def __init__(self, mol):
        super(ZeoVoronoiXYZ, self).__init__(mol)

    @staticmethod
    def from_string(contents):
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
        coord_patt = re.compile(
            "(\w+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+" +
            "([0-9\-\.]+)"
        )
        for i in range(2, 2 + num_sites):
            m = coord_patt.search(lines[i])
            if m:
                sp.append(m.group(1))  # this is 1-indexed
                # coords.append(map(float, m.groups()[1:4]))  # this is 0-indexed
                coords.append([float(j)
                               for j in [m.group(i) for i in [3, 4, 2]]])
                prop.append(float(m.group(5)))
        return ZeoVoronoiXYZ(
            Molecule(sp, coords, site_properties={'voronoi_radius': prop})
        )

    @staticmethod
    def from_file(filename):
        """
        Creates XYZ object from a file.

        Args:
            filename: XYZ filename

        Returns:
            XYZ object
        """
        with zopen(filename) as f:
            return ZeoVoronoiXYZ.from_string(f.read())

    def __str__(self):
        output = [str(len(self._mol)), self._mol.composition.formula]
        fmtstr = "{{}} {{:.{0}f}} {{:.{0}f}} {{:.{0}f}} {{:.{0}f}}".format(
            self.precision
        )
        for site in self._mol:
            output.append(fmtstr.format(
                site.specie.symbol, site.z, site.x, site.y,
                # site.specie, site.x, site.y, site.z,
                site.properties['voronoi_radius']
            ))
        return "\n".join(output)


@requires(zeo_found,
          "get_voronoi_nodes requires Zeo++ cython extension to be "
          "installed. Please contact developers of Zeo++ to obtain it.")
def get_voronoi_nodes(structure, rad_dict=None, probe_rad=0.1):
    """
    Analyze the void space in the input structure using voronoi decomposition
    Calls Zeo++ for Voronoi decomposition.

    Args:
        structure: pymatgen.core.structure.Structure
        rad_dict (optional): Dictionary of radii of elements in structure.
            If not given, Zeo++ default values are used.
            Note: Zeo++ uses atomic radii of elements.
            For ionic structures, pass rad_dict with ionic radii
        probe_rad (optional): Sampling probe radius in Angstroms. Default is
            0.1 A

    Returns:
        voronoi nodes as pymatgen.core.structure.Strucutre within the
        unit cell defined by the lattice of input structure
        voronoi face centers as pymatgen.core.structure.Strucutre within the
        unit cell defined by the lattice of input structure
    """

    with ScratchDir('.'):
        name = "temp_zeo1"
        zeo_inp_filename = name + ".cssr"
        ZeoCssr(structure).write_file(zeo_inp_filename)
        rad_file = None
        rad_flag = False

        if rad_dict:
            rad_file = name + ".rad"
            rad_flag = True
            with open(rad_file, 'w+') as fp:
                for el in rad_dict.keys():
                    fp.write("{} {}\n".format(el, rad_dict[el].real))

        atmnet = AtomNetwork.read_from_CSSR(
            zeo_inp_filename, rad_flag=rad_flag, rad_file=rad_file)
        vornet, vor_edge_centers, vor_face_centers = \
            atmnet.perform_voronoi_decomposition()
        vornet.analyze_writeto_XYZ(name, probe_rad, atmnet)
        voro_out_filename = name + '_voro.xyz'
        voro_node_mol = ZeoVoronoiXYZ.from_file(voro_out_filename).molecule

    species = ["X"] * len(voro_node_mol.sites)
    coords = []
    prop = []
    for site in voro_node_mol.sites:
        coords.append(list(site.coords))
        prop.append(site.properties['voronoi_radius'])

    lattice = Lattice.from_lengths_and_angles(
        structure.lattice.abc, structure.lattice.angles)
    vor_node_struct = Structure(
        lattice, species, coords, coords_are_cartesian=True,
        to_unit_cell=True, site_properties={"voronoi_radius": prop})

    # PMG-Zeo c<->a transformation for voronoi face centers
    rot_face_centers = [(center[1], center[2], center[0]) for center in
                        vor_face_centers]
    rot_edge_centers = [(center[1], center[2], center[0]) for center in
                        vor_edge_centers]

    species = ["X"] * len(rot_face_centers)
    prop = [0.0] * len(rot_face_centers)  # Vor radius not evaluated for fc
    vor_facecenter_struct = Structure(
        lattice, species, rot_face_centers, coords_are_cartesian=True,
        to_unit_cell=True, site_properties={"voronoi_radius": prop})

    species = ["X"] * len(rot_edge_centers)
    prop = [0.0] * len(rot_edge_centers)  # Vor radius not evaluated for fc
    vor_edgecenter_struct = Structure(
        lattice, species, rot_edge_centers, coords_are_cartesian=True,
        to_unit_cell=True, site_properties={"voronoi_radius": prop})

    return vor_node_struct, vor_edgecenter_struct, vor_facecenter_struct


def get_high_accuracy_voronoi_nodes(structure, rad_dict, probe_rad=0.1):
    """
    Analyze the void space in the input structure using high accuracy
    voronoi decomposition.
    Calls Zeo++ for Voronoi decomposition.

    Args:
        structure: pymatgen.core.structure.Structure
        rad_dict (optional): Dictionary of radii of elements in structure.
            If not given, Zeo++ default values are used.
            Note: Zeo++ uses atomic radii of elements.
            For ionic structures, pass rad_dict with ionic radii
        probe_rad (optional): Sampling probe radius in Angstroms.
            Default is 0.1 A

    Returns:
        voronoi nodes as pymatgen.core.structure.Strucutre within the
        unit cell defined by the lattice of input structure
        voronoi face centers as pymatgen.core.structure.Strucutre within the
        unit cell defined by the lattice of input structure
    """

    with ScratchDir('.'):
        name = "temp_zeo1"
        zeo_inp_filename = name + ".cssr"
        ZeoCssr(structure).write_file(zeo_inp_filename)
        rad_flag = True
        rad_file = name + ".rad"
        with open(rad_file, 'w+') as fp:
            for el in rad_dict.keys():
                print("{} {}".format(el, rad_dict[el].real), file=fp)

        atmnet = AtomNetwork.read_from_CSSR(
            zeo_inp_filename, rad_flag=rad_flag, rad_file=rad_file)
        # vornet, vor_edge_centers, vor_face_centers = \
        #        atmnet.perform_voronoi_decomposition()
        red_ha_vornet = \
            prune_voronoi_network_close_node(atmnet)
        # generate_simplified_highaccuracy_voronoi_network(atmnet)
        # get_nearest_largest_diameter_highaccuracy_vornode(atmnet)
        red_ha_vornet.analyze_writeto_XYZ(name, probe_rad, atmnet)
        voro_out_filename = name + '_voro.xyz'
        voro_node_mol = ZeoVoronoiXYZ.from_file(voro_out_filename).molecule

    species = ["X"] * len(voro_node_mol.sites)
    coords = []
    prop = []
    for site in voro_node_mol.sites:
        coords.append(list(site.coords))
        prop.append(site.properties['voronoi_radius'])

    lattice = Lattice.from_lengths_and_angles(
        structure.lattice.abc, structure.lattice.angles)
    vor_node_struct = Structure(
        lattice, species, coords, coords_are_cartesian=True,
        to_unit_cell=True, site_properties={"voronoi_radius": prop})

    return vor_node_struct


@requires(zeo_found,
          "get_voronoi_nodes requires Zeo++ cython extension to be "
          "installed. Please contact developers of Zeo++ to obtain it.")
def get_free_sphere_params(structure, rad_dict=None, probe_rad=0.1):
    """
    Analyze the void space in the input structure using voronoi decomposition
    Calls Zeo++ for Voronoi decomposition.

    Args:
        structure: pymatgen.core.structure.Structure
        rad_dict (optional): Dictionary of radii of elements in structure.
            If not given, Zeo++ default values are used.
            Note: Zeo++ uses atomic radii of elements.
            For ionic structures, pass rad_dict with ionic radii
        probe_rad (optional): Sampling probe radius in Angstroms. Default is
            0.1 A

    Returns:
        voronoi nodes as pymatgen.core.structure.Strucutre within the
        unit cell defined by the lattice of input structure
        voronoi face centers as pymatgen.core.structure.Strucutre within the
        unit cell defined by the lattice of input structure
    """

    with ScratchDir('.'):
        name = "temp_zeo1"
        zeo_inp_filename = name + ".cssr"
        ZeoCssr(structure).write_file(zeo_inp_filename)
        rad_file = None
        rad_flag = False

        if rad_dict:
            rad_file = name + ".rad"
            rad_flag = True
            with open(rad_file, 'w+') as fp:
                for el in rad_dict.keys():
                    fp.write("{} {}\n".format(el, rad_dict[el].real))

        atmnet = AtomNetwork.read_from_CSSR(
            zeo_inp_filename, rad_flag=rad_flag, rad_file=rad_file)
        out_file = "temp.res"
        atmnet.calculate_free_sphere_parameters(out_file)
        if os.path.isfile(out_file) and os.path.getsize(out_file) > 0:
            with open(out_file, "rt") as fp:
                output = fp.readline()
        else:
            output = ""
    fields = [val.strip() for val in output.split()][1:4]
    if len(fields) == 3:
        fields = [float(field) for field in fields]
        free_sphere_params = {'inc_sph_max_dia': fields[0],
                              'free_sph_max_dia': fields[1],
                              'inc_sph_along_free_sph_path_max_dia': fields[2]}
    return free_sphere_params


# Deprecated. Not needed anymore
def get_void_volume_surfarea(structure, rad_dict=None, chan_rad=0.3,
                             probe_rad=0.1):
    """
    Computes the volume and surface area of isolated void using Zeo++.
    Useful to compute the volume and surface area of vacant site.

    Args:
        structure: pymatgen Structure containing vacancy
        rad_dict(optional): Dictionary with short name of elements and their
            radii.
        chan_rad(optional): Minimum channel Radius.
        probe_rad(optional): Probe radius for Monte Carlo sampling.

    Returns:
        volume: floating number representing the volume of void
    """
    with ScratchDir('.'):
        name = "temp_zeo"
        zeo_inp_filename = name + ".cssr"
        ZeoCssr(structure).write_file(zeo_inp_filename)

        rad_file = None
        if rad_dict:
            rad_file = name + ".rad"
            with open(rad_file, 'w') as fp:
                for el in rad_dict.keys():
                    fp.write("{0}     {1}".format(el, rad_dict[el]))

        atmnet = AtomNetwork.read_from_CSSR(zeo_inp_filename, True, rad_file)
        vol_str = volume(atmnet, 0.3, probe_rad, 10000)
        sa_str = surface_area(atmnet, 0.3, probe_rad, 10000)
        vol = None
        sa = None
        for line in vol_str.split("\n"):
            if "Number_of_pockets" in line:
                fields = line.split()
                if float(fields[1]) > 1:
                    vol = -1.0
                    break
                if float(fields[1]) == 0:
                    vol = -1.0
                    break
                vol = float(fields[3])
        for line in sa_str.split("\n"):
            if "Number_of_pockets" in line:
                fields = line.split()
                if float(fields[1]) > 1:
                    # raise ValueError("Too many voids")
                    sa = -1.0
                    break
                if float(fields[1]) == 0:
                    sa = -1.0
                    break
                sa = float(fields[3])

    if not vol or not sa:
        raise ValueError("Error in zeo++ output stream")
    return vol, sa
