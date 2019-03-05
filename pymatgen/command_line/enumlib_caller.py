# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import re
import math
import subprocess
import itertools
import logging
import glob
import warnings

import numpy as np
from monty.fractions import lcm
import fractions

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import DummySpecie
from monty.os.path import which
from monty.dev import requires
from monty.tempfile import ScratchDir
from threading import Timer

"""
This module implements an interface to enumlib, Gus Hart"s excellent Fortran
code for enumerating derivative structures.

This module depends on a compiled enumlib with the executables enum.x and
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

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 16, 2012"


logger = logging.getLogger(__name__)


# Favor the use of the newer "enum.x" by Gus Hart instead of the older
# "multienum.x"
enum_cmd = which('enum.x') or which('multienum.x')
# prefer makestr.x at present
makestr_cmd = which('makestr.x') or which('makeStr.x') or which('makeStr.py')


@requires(enum_cmd and makestr_cmd,
          "EnumlibAdaptor requires the executables 'enum.x' or 'multienum.x' "
          "and 'makestr.x' or 'makeStr.py' to be in the path. Please download the "
          "library at http://enum.sourceforge.net/ and follow the instructions in "
          "the README to compile these two executables accordingly.")
class EnumlibAdaptor:
    """
    An adaptor for enumlib.

    .. attribute:: structures

        List of all enumerated structures.
    """
    amount_tol = 1e-5

    def __init__(self, structure, min_cell_size=1, max_cell_size=1,
                 symm_prec=0.1, enum_precision_parameter=0.001,
                 refine_structure=False, check_ordered_symmetry=True,
                 timeout=None):
        """
        Initializes the adapter with a structure and some parameters.

        Args:
            structure: An input structure.
            min_cell_size (int): The minimum cell size wanted. Defaults to 1.
            max_cell_size (int): The maximum cell size wanted. Defaults to 1.
            symm_prec (float): Symmetry precision. Defaults to 0.1.
            enum_precision_parameter (float): Finite precision parameter for
                enumlib. Default of 0.001 is usually ok, but you might need to
                tweak it for certain cells.
            refine_structure (bool): If you are starting from a structure that
                has been relaxed via some electronic structure code,
                it is usually much better to start with symmetry determination
                and then obtain a refined structure. The refined structure have
                cell parameters and atomic positions shifted to the expected
                symmetry positions, which makes it much less sensitive precision
                issues in enumlib. If you are already starting from an
                experimental cif, refinement should have already been done and
                it is not necessary. Defaults to False.
            check_ordered_symmetry (bool): Whether to check the symmetry of
                the ordered sites. If the symmetry of the ordered sites is
                lower, the lowest symmetry ordered sites is included in the
                enumeration. This is important if the ordered sites break
                symmetry in a way that is important getting possible
                structures. But sometimes including ordered sites
                slows down enumeration to the point that it cannot be
                completed. Switch to False in those cases. Defaults to True.
            timeout (float): If specified, will kill enumlib after specified
                time in minutes. This can be useful for gracefully handling
                enumerations in a high-throughput context, for some enumerations
                which will not terminate in a realistic length of time.
        """
        if refine_structure:
            finder = SpacegroupAnalyzer(structure, symm_prec)
            self.structure = finder.get_refined_structure()
        else:
            self.structure = structure
        self.min_cell_size = min_cell_size
        self.max_cell_size = max_cell_size
        self.symm_prec = symm_prec
        self.enum_precision_parameter = enum_precision_parameter
        self.check_ordered_symmetry = check_ordered_symmetry
        self.structures = None
        self.timeout = timeout

    def run(self):
        """
        Run the enumeration.
        """
        # Create a temporary directory for working.
        with ScratchDir(".") as d:
            logger.debug("Temp dir : {}".format(d))
            # Generate input files
            self._gen_input_file()
            # Perform the actual enumeration
            num_structs = self._run_multienum()
            # Read in the enumeration output as structures.
            if num_structs > 0:
                self.structures = self._get_structures(num_structs)
            else:
                raise EnumError("Unable to enumerate structure.")

    def _gen_input_file(self):
        """
        Generate the necessary struct_enum.in file for enumlib. See enumlib
        documentation for details.
        """
        coord_format = "{:.6f} {:.6f} {:.6f}"
        # Using symmetry finder, get the symmetrically distinct sites.
        fitter = SpacegroupAnalyzer(self.structure, self.symm_prec)
        symmetrized_structure = fitter.get_symmetrized_structure()
        logger.debug("Spacegroup {} ({}) with {} distinct sites".format(
            fitter.get_space_group_symbol(),
            fitter.get_space_group_number(),
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

        # Stores the ordered sites, which are not enumerated.
        ordered_sites = []
        disordered_sites = []
        coord_str = []
        for sites in symmetrized_structure.equivalent_sites:
            if sites[0].is_ordered:
                ordered_sites.append(sites)
            else:
                sp_label = []
                species = {k: v for k, v in sites[0].species.items()}
                if sum(species.values()) < 1 - EnumlibAdaptor.amount_tol:
                    # Let us first make add a dummy element for every single
                    # site whose total occupancies don't sum to 1.
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

        def get_sg_info(ss):
            finder = SpacegroupAnalyzer(Structure.from_sites(ss),
                                        self.symm_prec)
            return finder.get_space_group_number()

        target_sgnum = get_sg_info(symmetrized_structure.sites)
        curr_sites = list(itertools.chain.from_iterable(disordered_sites))
        sgnum = get_sg_info(curr_sites)
        ordered_sites = sorted(ordered_sites, key=lambda sites: len(sites))
        logger.debug("Disordered sites has sg # %d" % (sgnum))
        self.ordered_sites = []

        # progressively add ordered sites to our disordered sites
        # until we match the symmetry of our input structure
        if self.check_ordered_symmetry:
            while sgnum != target_sgnum and len(ordered_sites) > 0:
                sites = ordered_sites.pop(0)
                temp_sites = list(curr_sites) + sites
                new_sgnum = get_sg_info(temp_sites)
                if sgnum != new_sgnum:
                    logger.debug("Adding %s in enum. New sg # %d"
                                 % (sites[0].specie, new_sgnum))
                    index_species.append(sites[0].specie)
                    index_amounts.append(len(sites))
                    sp_label = len(index_species) - 1
                    for site in sites:
                        coord_str.append("{} {}".format(
                            coord_format.format(*site.coords),
                            sp_label))
                    disordered_sites.append(sites)
                    curr_sites = temp_sites
                    sgnum = new_sgnum
                else:
                    self.ordered_sites.extend(sites)

        for sites in ordered_sites:
            self.ordered_sites.extend(sites)

        self.index_species = index_species

        lattice = self.structure.lattice

        output = [self.structure.formula, "bulk"]
        for vec in lattice.matrix:
            output.append(coord_format.format(*vec))
        output.append("%d" % len(index_species))
        output.append("%d" % len(coord_str))
        output.extend(coord_str)

        output.append("{} {}".format(self.min_cell_size, self.max_cell_size))
        output.append(str(self.enum_precision_parameter))
        output.append("partial")

        ndisordered = sum([len(s) for s in disordered_sites])
        base = int(ndisordered*lcm(*[f.limit_denominator(ndisordered *
                                          self.max_cell_size).denominator
                                       for f in map(fractions.Fraction,
                                                    index_amounts)]))

        # This multiplicative factor of 10 is to prevent having too small bases
        # which can lead to rounding issues in the next step.
        # An old bug was that a base was set to 8, with a conc of 0.4:0.6. That
        # resulted in a range that overlaps and a conc of 0.5 satisfying this
        # enumeration. See Cu7Te5.cif test file.
        base *= 10

        # base = ndisordered #10 ** int(math.ceil(math.log10(ndisordered)))
        # To get a reasonable number of structures, we fix concentrations to the
        # range expected in the original structure.
        total_amounts = sum(index_amounts)
        for amt in index_amounts:
            conc = amt / total_amounts

            if abs(conc * base - round(conc * base)) < 1e-5:
                output.append("{} {} {}".format(int(round(conc * base)),
                                                int(round(conc * base)),
                                                base))
            else:
                min_conc = int(math.floor(conc * base))
                output.append("{} {} {}".format(min_conc - 1, min_conc + 1,
                                                base))
        output.append("")
        logger.debug("Generated input file:\n{}".format("\n".join(output)))
        with open("struct_enum.in", "w") as f:
            f.write("\n".join(output))

    def _run_multienum(self):

        p = subprocess.Popen([enum_cmd],
                             stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE, close_fds=True)

        if self.timeout:

            timed_out = False
            timer = Timer(self.timeout*60, lambda p: p.kill(), [p])

            try:
                timer.start()
                output = p.communicate()[0].decode("utf-8")
            finally:
                if not timer.is_alive():
                    timed_out = True
                timer.cancel()

            if timed_out:
                raise TimeoutError('Enumeration took too long.')

        else:

            output = p.communicate()[0].decode("utf-8")

        count = 0
        start_count = False
        for line in output.strip().split("\n"):
            if line.strip().endswith("RunTot"):
                start_count = True
            elif start_count and re.match(r"\d+\s+.*", line.strip()):
                count = int(line.split()[-1])
        logger.debug("Enumeration resulted in {} structures".format(count))
        return count

    def _get_structures(self, num_structs):
        structs = []

        if ".py" in makestr_cmd:
            options = ["-input", "struct_enum.out", str(1), str(num_structs)]
        else:
            options = ["struct_enum.out", str(0), str(num_structs - 1)]

        rs = subprocess.Popen([makestr_cmd] + options,
                              stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE, close_fds=True)
        stdout, stderr = rs.communicate()
        if stderr:
            logger.warning(stderr.decode())

        # sites retrieved from enumlib will lack site properties
        # to ensure consistency, we keep track of what site properties
        # are missing and set them to None
        # TODO: improve this by mapping ordered structure to original
        # disorded structure, and retrieving correct site properties
        disordered_site_properties = {}

        if len(self.ordered_sites) > 0:
            original_latt = self.ordered_sites[0].lattice
            # Need to strip sites of site_properties, which would otherwise
            # result in an index error. Hence Structure is reconstructed in
            # the next step.
            site_properties = {}
            for site in self.ordered_sites:
                for k, v in site.properties.items():
                    disordered_site_properties[k] = None
                    if k in site_properties:
                        site_properties[k].append(v)
                    else:
                        site_properties[k] = [v]
            ordered_structure = Structure(
                original_latt,
                [site.species for site in self.ordered_sites],
                [site.frac_coords for site in self.ordered_sites],
                site_properties=site_properties
            )
            inv_org_latt = np.linalg.inv(original_latt.matrix)

        for file in glob.glob('vasp.*'):
            with open(file) as f:
                data = f.read()
                data = re.sub(r'scale factor', "1", data)
                data = re.sub(r'(\d+)-(\d+)', r'\1 -\2', data)
                poscar = Poscar.from_string(data, self.index_species)
                sub_structure = poscar.structure
                # Enumeration may have resulted in a super lattice. We need to
                # find the mapping from the new lattice to the old lattice, and
                # perform supercell construction if necessary.
                new_latt = sub_structure.lattice

                sites = []

                if len(self.ordered_sites) > 0:
                    transformation = np.dot(new_latt.matrix, inv_org_latt)
                    transformation = [[int(round(cell)) for cell in row]
                                      for row in transformation]
                    logger.debug("Supercell matrix: {}".format(transformation))
                    s = ordered_structure * transformation
                    sites.extend([site.to_unit_cell() for site in s])
                    super_latt = sites[-1].lattice
                else:
                    super_latt = new_latt

                for site in sub_structure:
                    if site.specie.symbol != "X":  # We exclude vacancies.
                        sites.append(
                            PeriodicSite(site.species,
                                         site.frac_coords,
                                         super_latt,
                                         to_unit_cell=True,
                                         properties=disordered_site_properties)
                        )
                    else:
                        warnings.warn("Skipping sites that include species X.")
                structs.append(Structure.from_sites(sorted(sites)))

        logger.debug("Read in a total of {} structures.".format(num_structs))
        return structs


class EnumError(BaseException):
    pass
