#!/usr/bin/env python

"""
Module implementing classes and functions to use Zeo++.
Zeo++ can be obtained from http://www.maciejharanczyk.info/Zeopp/
"""

from __future__ import division

import re

#from zeo.netstorage import AtomNetwork

from pymatgen.io.cssrio import Cssr
from pymatgen.io.xyzio import XYZ
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.util.io_utils import zopen

class ZeoCssr(Cssr):
    """
    ZeoCssr adds extra fields to CSSR sites to conform with Zeo++ 
    input CSSR format
    """
    def __init__(self, structure):
        """
        Args:
            structure:
                A structure to create ZeoCssr object
        """
        super(ZeoCssr, self).__init__(structure)

    def __str__(self):
        output = [
                "{:.4f} {:.4f} {:.4f}"
                .format(*self.structure.lattice.abc),
                "{:.2f} {:.2f} {:.2f} SPGR =  1 P 1    OPT = 1"
                .format(*self.structure.lattice.angles),
                "{} 0".format(len(self.structure)),
                "0 {}".format(self.structure.formula)
                ]
        for i, site in enumerate(self.structure.sites):
            if not hasattr(site, 'charge'): 
                output.append(
                    "{} {} {:.4f} {:.4f} {:.4f} 0 0 0 0 0 0 0 0 {:.4f}"
                    .format(i+1, site.specie, site.a, site.b, site.c, 0.0)
                    )
            else:
                output.append(
                    "{} {} {:.4f} {:.4f} {:.4f} 0 0 0 0 0 0 0 0 {:.4f}"
                    .format(
                            i+1, site.specie, site.a, site.b, site.c, 
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
        latt = Lattice.from_lengths_and_angles(lengths, angles)
        sp = []
        coords = []
        chrg = []
        for l in lines[4:]:
            m = re.match("\d+\s+(\w+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+" +
                         "([0-9\-\.]+)\s+(?:0\s+){8}([0-9\-\.]+)", l.strip())
            if m:
                sp.append(m.group(1))
                coords.append([float(m.group(i)) for i in xrange(2, 5)])
                chrg.append(m.group(5))
        return ZeoCssr(
            Structure(latt, sp, coords, site_properties={'charge':chrg})
            )

    @staticmethod
    def from_file(filename):
        """
        Reads a CSSR file to a Cssr object.
        
        Args:
            filename:
                Filename to read from.
        
        Returns:
            Cssr object.
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
            "(\w+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+"+
            "([0-9\-\.]+)"
        )   
        for i in xrange(2, 2 + num_sites):
            m = coord_patt.search(lines[i])
            if m:
                sp.append(m.group(1))  # this is 1-indexed
                coords.append(map(float, m.groups()[1:4]))  # this is 0-indexed
                prop.append(m.group(5))
        print prop
        return ZeoVoronoiXYZ(
            Molecule(sp, coords, site_properties={'voronoi_radius':prop})
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
        

