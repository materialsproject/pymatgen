#!/usr/bin/env python

"""
This module implements an interface to enumlib, Gus Hart"s excellent Fortran
code for enumerating derivative structures.

This module depends on a compiled enumlib with the executables multienum.x and
makestr.x available in the path. Please download the library at
http://enum.sourceforge.net/ and follow the instructions in the README to
compile these two executables accordingly.

If you use this module, please cite the following:

Gus L. W. Hart and Rodney W. Forcade, "Algorithm for generating derivative
structures," Phys. Rev. B 77 224115 (26 June 2008)

Gus L. W. Hart and Rodney W. Forcade, "Generating derivative structures from
multilattices: Application to hcp alloys," Phys. Rev. B 80 014120 (July 2009)

Gus L. W. Hart, Lance J. Nelson, and Rodney W. Forcade, "Generating
derivative structures at a fixed concentration," Comp. Mat. Sci. 59
101-107 (March 2012)
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 16, 2012"

import os
import re
import math
import tempfile
import subprocess
import shutil
import logging

import numpy as np
import warnings

from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.util.io_utils import which
from pymatgen.util.decorators import requires

logger = logging.getLogger(__name__)


@requires(which('multienum.x') and which('makestr.x'),
          "EnumlibAdaptor requires the executables 'multienum.x' and "
          "'makestr.x' to be in the path. Please download the library at"
          "http://enum.sourceforge.net/ and follow the instructions in "
          "the README to compile these two executables accordingly.")
class EnumlibAdaptor(object):
    """
    An adaptor for enumlib.

    .. attribute:: structures

        List of all enumerated structures.
    """
    amount_tol = 1e-5

    def __init__(self, structure, min_cell_size=1, max_cell_size=1,
                 symm_prec=0.1, enum_precision_parameter=0.001,
                 refine_structure=False):
        """
        Args:
            structure:
                An input structure.
            min_cell_size:
                The minimum cell size wanted. Must be an int. Defaults to 1.
            max_cell_size:
                The maximum cell size wanted. Must be an int. Defaults to 1.
            symm_prec:
                Symmetry precision. Defaults to 0.1.
            enum_precision_parameter:
                Finite precision parameter for enumlib. Default of 0.001 is
                usually ok, but you might need to tweak it for certain cells.
            refine_structure:
                If you are starting from a structure that has been relaxed via
                some electronic structure code, it is usually much better to
                start with symmetry determination and then obtain a refined
                structure. The refined structure have cell parameters and
                atomic positions shifted to the expected symmetry positions,
                which makes it much less sensitive precision issues in enumlib.
                If you are already starting from an experimental cif, refinment
                should have already been done and it is not necessary. Defaults
                to False.
        """
        if refine_structure:
            finder = SymmetryFinder(structure, symm_prec)
            self.structure = finder.get_refined_structure()
        else:
            self.structure = structure
        self.min_cell_size = min_cell_size
        self.max_cell_size = max_cell_size
        self.symm_prec = symm_prec
        self.enum_precision_parameter = enum_precision_parameter

    def run(self):
        """
        Run the enumeration.
        """
        #Create a temporary directory for working.
        curr_dir = os.getcwd()
        temp_dir = tempfile.mkdtemp()
        logger.debug("Temp dir : {}".format(temp_dir))
        try:
            #Generate input files
            self._gen_input_file(temp_dir)
            #Perform the actual enumeration
            num_structs = self._run_multienum(temp_dir)
            #Read in the enumeration output as structures.
            if num_structs > 0:
                self.structures = self._get_structures(num_structs)
            else:
                raise ValueError("Unable to enumerate structure.")
        except Exception:
            import sys
            import traceback

            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                                      limit=2, file=sys.stdout)
        finally:
            try:
                shutil.rmtree(temp_dir)
            except:
                warnings.warn("Unable to delete temp dir "
                              "{}. Please delete manually".format(temp_dir))
            os.chdir(curr_dir)

    def _gen_input_file(self, working_dir):
        """
        Generate the necessary struct_enum.in file for enumlib. See enumlib
        documentation for details.
        """
        coord_format = "{:.6f} {:.6f} {:.6f}"

        # Using symmetry finder, get the symmetrically distinct sites.
        fitter = SymmetryFinder(self.structure, self.symm_prec)
        symmetrized_structure = fitter.get_symmetrized_structure()
        logger.debug("Spacegroup {} ({}) with {} distinct sites".format(
            fitter.get_spacegroup_symbol(),
            fitter.get_spacegroup_number(),
            len(symmetrized_structure.equivalent_sites))
        )

        """
        Enumlib doesn"t work when the number of species get too large. To
        simplify matters, we generate the input file only with disordered sites
        and exclude the ordered sites from the enumeration. The fact that
        different disordered sites with the exact same species may belong to
        different equivalent sites is dealt with by having determined the
        spacegroup earlier and labelling the species differently.
        """

        # index_species and index_amounts store mappings between the indices
        # used in the enum input file, and the actual species and amounts.
        index_species = []
        index_amounts = []

        #Let"s group and sort the sites by symmetry.
        site_symmetries = []
        for sites in symmetrized_structure.equivalent_sites:
            finder = SymmetryFinder(Structure.from_sites(sites),
                                    self.symm_prec)
            sgnum = finder.get_spacegroup_number()
            site_symmetries.append((sites, sgnum))

        site_symmetries = sorted(site_symmetries, key=lambda s: s[1])

        #Stores the ordered sites, which are not enumerated.
        min_sg_num = site_symmetries[0][1]
        ordered_sites = []
        disordered_sites = []
        coord_str = []
        min_disordered_sg = 300
        for (sites, sgnum) in site_symmetries:
            if sites[0].is_ordered:
                ordered_sites.append(sites)
            else:
                min_disordered_sg = min(min_disordered_sg, sgnum)
                sp_label = []
                species = {k: v for k, v in sites[0].species_and_occu.items()}
                if sum(species.values()) < 1 - EnumlibAdaptor.amount_tol:
                    #Let us first make add a dummy element for every single
                    #site whose total occupancies don't sum to 1.
                    species[DummySpecie("X")] = 1 - sum(species.values())
                for sp in species.keys():
                    if sp not in index_species:
                        index_species.append(sp)
                        sp_label.append(len(index_species) - 1)
                        index_amounts.append(species[sp] * len(sites))
                    else:
                        ind = index_species.index(sp)
                        sp_label.append(ind)
                        index_amounts[ind] += species[sp] * len(sites)
                sp_label = "/".join(["{}".format(i) for i in sorted(sp_label)])
                for site in sites:
                    coord_str.append("{} {}".format(
                        coord_format.format(*site.coords),
                        sp_label))
                disordered_sites.append(sites)

        #It could be that some of the ordered sites has a lower symmetry than
        #the disordered sites.  So we consider the lowest symmetry sites as
        #disordered in our enumeration.
        if min_disordered_sg > min_sg_num:
            logger.debug("Ordered sites have lower symmetry than disordered.")
            sites = ordered_sites.pop(0)
            index_species.append(sites[0].specie)
            index_amounts.append(len(sites))
            sp_label = len(index_species) - 1
            logger.debug("Lowest symmetry {} sites are included in enum."
                         .format(sites[0].specie))
            for site in sites:
                coord_str.append("{} {}".format(
                    coord_format.format(*site.coords),
                    sp_label))
            disordered_sites.append(sites)

        self.ordered_sites = []
        for sites in ordered_sites:
            self.ordered_sites.extend(sites)
        self.index_species = index_species

        lattice = self.structure.lattice

        output = [self.structure.formula, "bulk"]
        for vec in lattice.matrix:
            output.append(coord_format.format(*vec))
        output.append("{}".format(len(index_species)))
        output.append("{}".format(len(coord_str)))
        output.extend(coord_str)

        output.append("{} {}".format(self.min_cell_size, self.max_cell_size))
        output.append(str(self.enum_precision_parameter))
        output.append("partial")
        #To get a reasonable number of structures, we fix concentrations to the
        #range expected in the original structure.
        total_amounts = sum(index_amounts)
        for amt in index_amounts:
            conc = amt / total_amounts
            if abs(conc * 100 - round(conc * 100)) < 1e-5:
                output.append("{} {} {}".format(int(round(conc * 100)),
                                                int(round(conc * 100)), 100))
            else:
                min_conc = int(math.floor(conc * 100))
                output.append("{} {} {}".format(min_conc - 1, min_conc + 1,
                                                100))
        output.append("")
        logger.debug("Generated input file:\n{}".format("\n".join(output)))
        with open(os.path.join(working_dir, "struct_enum.in"), "w") as f:
            f.write("\n".join(output))

    def _run_multienum(self, working_dir):
        os.chdir(working_dir)
        p = subprocess.Popen(["multienum.x"],
                             stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE, close_fds=True)
        output = p.communicate()[0]
        count = 0
        start_count = False
        for line in output.strip().split("\n"):
            if line.strip().endswith("RunTot"):
                start_count = True
            elif start_count and re.match("\d+\s+.*", line.strip()):
                count = int(line.split()[-1])
        logger.debug("Enumeration resulted in {} structures".format(count))
        return count

    def _get_structures(self, num_structs):
        structs = []
        rs = subprocess.Popen(["makestr.x",
                               "struct_enum.out", str(0),
                               str(num_structs - 1)],
                              stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE, close_fds=True)
        rs.communicate()
        if len(self.ordered_sites) > 0:
            original_latt = self.ordered_sites[0].lattice
            # Need to strip sites of site_properties, which would otherwise
            # result in an index error. Hence Structure is reconstructed in
            # the next step.
            ordered_structure = Structure(
                original_latt,
                [site.species_and_occu for site in self.ordered_sites],
                [site.frac_coords for site in self.ordered_sites])
            inv_org_latt = np.linalg.inv(original_latt.matrix)

        for n in range(1, num_structs + 1):
            with open("vasp.{:06d}".format(n)) as f:
                data = f.read()
                data = re.sub("scale factor", "1", data)
                data = re.sub("(\d+)-(\d+)", r"\1 -\2", data)
                poscar = Poscar.from_string(data, self.index_species)
                sub_structure = poscar.structure
                #Enumeration may have resulted in a super lattice. We need to
                #find the mapping from the new lattice to the old lattice, and
                #perform supercell construction if necessary.
                new_latt = sub_structure.lattice

                sites = []

                if len(self.ordered_sites) > 0:
                    transformation = np.dot(new_latt.matrix, inv_org_latt)
                    transformation = [[int(round(cell)) for cell in row]
                                      for row in transformation]
                    logger.debug("Supercell matrix: {}".format(transformation))
                    s = Structure.from_sites(ordered_structure)
                    s.make_supercell(transformation)
                    sites.extend([site.to_unit_cell for site in s])
                    super_latt = sites[-1].lattice
                else:
                    super_latt = new_latt

                for site in sub_structure:
                    if site.specie.symbol != "X":  # We exclude vacancies.
                        sites.append(PeriodicSite(site.species_and_occu,
                                                  site.frac_coords,
                                                  super_latt).to_unit_cell)
                structs.append(Structure.from_sites(sorted(sites)))

        logger.debug("Read in a total of {} structures.".format(num_structs))
        return structs
