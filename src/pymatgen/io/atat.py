"""Classes for reading/writing mcsqs files following the rndstr.in format."""

from __future__ import annotations

import numpy as np

from pymatgen.core import Lattice, Structure, get_el_sp

__author__ = "Matthew Horton"
__copyright__ = "Copyright 2017, The Materials Project"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Production"
__date__ = "October 2017"


class Mcsqs:
    """
    Handle input/output for the crystal definition format
    used by mcsqs and other ATAT codes.
    """

    def __init__(self, structure: Structure):
        """
        Args:
            Structure: input Structure.
        """
        self.structure = structure

    def to_str(self):
        """
        Returns:
            str: a structure in mcsqs rndstr.in format.
        """
        # add lattice vectors
        mat = self.structure.lattice.matrix
        output = [f"{vec[0]:6f} {vec[1]:6f} {vec[2]:6f}" for vec in mat]

        # define coord system, use Cartesian
        output.extend(("1.0 0.0 0.0", "0.0 1.0 0.0", "0.0 0.0 1.0"))

        # add species
        for site in self.structure:
            species_str = []
            for sp, occu in sorted(site.species.items()):
                sp = str(sp)
                if ("," in sp) or ("=" in sp):
                    # Invalid species string for AT-AT input, so modify
                    sp = sp.replace(",", "__").replace("=", "___")
                species_str.append(f"{sp}={occu}")
            species_str = ",".join(species_str)
            a, b, c = site.frac_coords
            output.append(f"{a:6f} {b:6f} {c:6f} {species_str}")

        return "\n".join(output)

    @staticmethod
    def structure_from_str(data):
        """
        Parses a rndstr.in, lat.in or bestsqs.out file into pymatgen's
        Structure format.

        Args:
            data: contents of a rndstr.in, lat.in or bestsqs.out file

        Returns:
            Structure object
        """
        data = data.splitlines()
        data = [x.split() for x in data if x]  # remove empty lines

        # following specification/terminology given in manual
        if len(data[0]) == 6:  # lattice parameters
            a, b, c, alpha, beta, gamma = map(float, data[0])
            coord_system = Lattice.from_parameters(a, b, c, alpha, beta, gamma).matrix
            lattice_vecs = np.array(
                [
                    [data[1][0], data[1][1], data[1][2]],
                    [data[2][0], data[2][1], data[2][2]],
                    [data[3][0], data[3][1], data[3][2]],
                ],
                dtype=float,
            )
            first_species_line = 4
        else:
            coord_system = np.array(
                [
                    [data[0][0], data[0][1], data[0][2]],
                    [data[1][0], data[1][1], data[1][2]],
                    [data[2][0], data[2][1], data[2][2]],
                ],
                dtype=float,
            )
            lattice_vecs = np.array(
                [
                    [data[3][0], data[3][1], data[3][2]],
                    [data[4][0], data[4][1], data[4][2]],
                    [data[5][0], data[5][1], data[5][2]],
                ],
                dtype=float,
            )
            first_species_line = 6

        scaled_matrix = np.matmul(lattice_vecs, coord_system)
        lattice = Lattice(scaled_matrix)

        all_coords = []
        all_species = []
        for line in data[first_species_line:]:
            coords = np.array([line[0], line[1], line[2]], dtype=float)
            scaled_coords = np.matmul(coords, np.linalg.inv(lattice_vecs))
            all_coords.append(scaled_coords)

            species_strs = "".join(line[3:])  # join multiple strings back together
            species_strs = species_strs.replace(" ", "")  # trim any white space
            species_strs = species_strs.split(",")  # comma-delimited

            species = {}

            for species_occ in species_strs:
                # gets a species, occupancy pair
                species_occ = species_occ.split("=")

                if len(species_occ) == 1:
                    # assume occupancy is 1.0
                    species_occ = [species_occ[0], 1.0]

                if "_" in species_occ[0]:
                    # see to_str() method in this file, since , and = are not valid
                    # species names in AT-AT we replace "," with "__" and "=" with "___",
                    # for pymatgen to parse these back correctly we have to replace them back
                    species_occ[0] = species_occ[0].replace("___", "=").replace("__", ",")

                species[get_el_sp(species_occ[0])] = float(species_occ[1])

            all_species.append(species)

        return Structure(lattice, all_species, all_coords)
