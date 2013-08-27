#!/usr/bin/env python

"""
Module implementing classes and functions to use Zeo++.
Zeo++ can be obtained from http://www.maciejharanczyk.info/Zeopp/
"""

from __future__ import division

__author__ = "Bharat Medasani"
__copyright = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Bharat Medasani"
__email__ = "bkmedasani@lbl.gov"
__data__ = "Aug 2, 2013"

import re
import tempfile
import os
import shutil

from pymatgen.io.cssrio import Cssr
from pymatgen.io.xyzio import XYZ
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.util.io_utils import zopen
from pymatgen.util.decorators import requires

try:
    from zeo.netstorage import AtomNetwork, VoronoiNetwork
    from zeo.area_volume import volume, surface_area
    zeo_found = True
except ImportError:
    zeo_found = False



class ZeoCssr(Cssr):
    """
    ZeoCssr adds extra fields to CSSR sites to conform with Zeo++ 
    input CSSR format. The coordinate system is rorated from xyz to zyx. 
    This change aligns the pivot axis of pymatgen (z-axis) to pivot axis 
    of Zeo++ (x-axis) for structurural modifications.
    """

    @requires(zeo_found,
              "ZeoCssr requires Zeo++ cython extension to be installed. Please "
              "contact developers of Zeo++ to obtain it.")
    def __init__(self, structure):
        """
        Args:
            structure:
                A structure to create ZeoCssr object
        """
        super(ZeoCssr, self).__init__(structure)

    def __str__(self):
        """
        CSSR.__str__ method is modified to padd 0's to the CSSR site data.
        The padding is to conform with the CSSR format supported Zeo++.
        Also coordinate system is rotated from xyz to zxy
        """
        output = [
            "{:.4f} {:.4f} {:.4f}"
            #.format(*self.structure.lattice.abc),
            .format(self.structure.lattice.c,
                    self.structure.lattice.a,
                    self.structure.lattice.b),
            "{:.2f} {:.2f} {:.2f} SPGR =  1 P 1    OPT = 1"
            #.format(*self.structure.lattice.angles),
            .format(self.structure.lattice.gamma,
                    self.structure.lattice.alpha,
                    self.structure.lattice.beta),
            "{} 0".format(len(self.structure)),
            "0 {}".format(self.structure.formula)
        ]
        for i, site in enumerate(self.structure.sites):
            if not hasattr(site, 'charge'):
                output.append(
                    "{} {} {:.4f} {:.4f} {:.4f} 0 0 0 0 0 0 0 0 {:.4f}"
                    .format(i + 1, site.specie, site.c, site.a, site.b, 0.0)
                    #.format(i+1, site.specie, site.a, site.b, site.c, 0.0)
                )
            else:
                output.append(
                    "{} {} {:.4f} {:.4f} {:.4f} 0 0 0 0 0 0 0 0 {:.4f}"
                    .format(
                        i + 1, site.specie, site.c, site.a, site.b,
                        #i+1, site.specie, site.a, site.b, site.c,
                        site.charge
                    )
                )

        return "\n".join(output)

    @staticmethod
    def from_string(string):
        """ 
        Reads a string representation to a ZeoCssr object.

        Args:
            string:
                A string representation of a ZeoCSSR.

        Returns:
            ZeoCssr object.
        """
        lines = string.split("\n")
        toks = lines[0].split()
        lengths = map(float, toks)
        toks = lines[1].split()
        angles = map(float, toks[0:3])
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
                #coords.append([float(m.group(i)) for i in xrange(2, 5)])
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
            filename:
                Filename to read from.
        
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
    """

    def __init__(self, mol):
        """
        Args:
            mol:
                Input molecule holding the voronoi node information
        """
        super(ZeoVoronoiXYZ, self).__init__(mol)

    @staticmethod
    def from_string(contents):
        """
        Creates Zeo++ Voronoi XYZ object from a string.
        from_string method of XYZ class is being redefined.

        Args:
            contents:
                String representing Zeo++ Voronoi XYZ file.

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
        for i in xrange(2, 2 + num_sites):
            m = coord_patt.search(lines[i])
            if m:
                sp.append(m.group(1))  # this is 1-indexed
                #coords.append(map(float, m.groups()[1:4]))  # this is 0-indexed
                coords.append(map(float, [m.group(i) for i in
                                          [3, 4, 2]]))  # this is 0-indexed
                prop.append(float(m.group(5)))
        return ZeoVoronoiXYZ(
            Molecule(sp, coords, site_properties={'voronoi_radius': prop})
        )

    @staticmethod
    def from_file(filename):
        """
        Creates XYZ object from a file.

        Args:
            filename:
                XYZ filename

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
                site.specie, site.z, site.x, site.y,
                #site.specie, site.x, site.y, site.z,
                site.properties['voronoi_radius']
            ))
        return "\n".join(output)


def get_voronoi_nodes(structure, rad_dict=None, probe_rad=0.1):
    """
    Analyze the void space in the input structure using voronoi decomposition
    Calls Zeo++ for Voronoi decomposition

    Args:
        structure:
            pymatgen.core.structure.Structure
        rad_dict (optional):
            Dictionary of radii of elements in structure. 
            If not given, Zeo++ default values are used.
            Note: Zeo++ uses atomic radii of elements.
            For ionic structures, pass rad_dict with ionic radii
        probe_rad (optional):
            Sampling probe radius in Angstroms. Default is 0.1 A

    Returns:
        voronoi nodes as pymatgen.core.structure.Strucutre within the 
        unit cell defined by the lattice of input structure 
    """

    temp_dir = tempfile.mkdtemp()
    current_dir = os.getcwd()
    name = "temp_zeo"
    zeo_inp_filename = name + ".cssr"
    os.chdir(temp_dir)
    ZeoCssr(structure).write_file(zeo_inp_filename)
    rad_file = None
    if rad_dict:
        rad_file = name + ".rad"
        with open(rad_file, 'w+') as fp:
            for el in rad_dict.keys():
                fp.write("{0} {1}".format(el, rad_dict[el]))

    atmnet = AtomNetwork.read_from_CSSR(zeo_inp_filename, True, rad_file)
    vornet = atmnet.perform_voronoi_decomposition()
    vornet.analyze_writeto_XYZ(name, probe_rad, atmnet)
    voronoi_out_filename = name + '_voro.xyz'
    voronoi_node_mol = ZeoVoronoiXYZ.from_file(voronoi_out_filename).molecule
    #print voronoi_node_mol
    species = ["X"] * len(voronoi_node_mol.sites)
    coords = []
    prop = []
    for site in voronoi_node_mol.sites:
        coords.append(list(site.coords))
        prop.append(site.properties['voronoi_radius'])

    lattice = Lattice.from_lengths_and_angles(
        structure.lattice.abc, structure.lattice.angles
    )
    voronoi_node_struct = Structure(
        lattice, species, coords, coords_are_cartesian=True,
        site_properties={"voronoi_radius": prop}
    )

    os.chdir(current_dir)
    shutil.rmtree(temp_dir)

    return voronoi_node_struct


def get_void_volume_surfarea(structure, rad_dict=None, chan_rad=0.2,
                             probe_rad=0.1):
    """
    Computes the volume and surface area of isolated void using Zeo++.
    Useful to compute the volume and surface area of vacant site.

    Args:
        structure:
            pymatgen Structure containing vacancy
        rad_dict(optional):
            Dictionary with short name of elements and their radii.
        chan_rad(optional):
            Minimum channel Radius. 
        probe_rad(optional)
            Probe radius for Monte Carlo sampling.

    Returns:
        volume:
            floating number representing the volume of void
    """
    temp_dir = tempfile.mkdtemp()
    current_dir = os.getcwd()
    name = "temp_zeo"
    zeo_inp_filename = name + ".cssr"
    os.chdir(temp_dir)
    ZeoCssr(structure).write_file(zeo_inp_filename)

    rad_file = None
    if rad_dict:
        rad_file = name + ".rad"
        with open(rad_file, 'w+') as fp:
            for el in rad_dict.keys():
                fp.write("{0}     {1}".format(el, rad_dict[el]))

    atmnet = AtomNetwork.read_from_CSSR(zeo_inp_filename, True, rad_file)
    vol_str = volume(atmnet, 0.3, probe_rad, 5000)
    sa_str = surface_area(atmnet, 0.3, probe_rad, 5000)
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
                #raise ValueError("Too many voids")
                sa = -1.0
                break
            if float(fields[1]) == 0:
                sa = -1.0
                break
            sa = float(fields[3])
    if not vol or not sa:
        raise ValueError("Error in zeo++ output stream")

    os.chdir(current_dir)
    shutil.rmtree(temp_dir)
    return vol, sa
