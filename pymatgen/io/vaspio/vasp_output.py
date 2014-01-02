#!/usr/bin/env python

"""
Classes for reading/manipulating/writing VASP ouput files.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Rickard Armiento, " + \
    "Vincent L Chevrier"
__credits__ = "Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.2"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Nov 30, 2012"

import os
import glob
import re
import math
import itertools
import warnings
import xml.sax.handler
import StringIO
import logging
from collections import defaultdict

import numpy as np

from pymatgen.util.coord_utils import get_points_in_sphere_pbc
from pymatgen.util.io_utils import zopen, clean_lines, micro_pyawk, \
    clean_json, reverse_readline
from pymatgen.core.structure import Structure
from pymatgen.core.units import unitized
from pymatgen.core.composition import Composition
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.electronic_structure.bandstructure import BandStructure, \
    BandStructureSymmLine, get_reconstructed_band_structure
from pymatgen.core.lattice import Lattice
from pymatgen.io.vaspio.vasp_input import Incar, Kpoints, Poscar

logger = logging.getLogger(__name__)


class Vasprun(object):
    """
    Vastly improved sax-based parser for vasprun.xml files.
    Speedup over Dom is at least 2x for smallish files (~1Mb) to orders of
    magnitude for larger files (~10Mb). All data is stored as attributes, which
    are delegated to the VasprunHandler object. Note that the results would
    differ depending on whether the read_electronic_structure option is set to
    True.

    **Vasp results**

    .. attribute:: ionic_steps

        All ionic steps in the run as a list of
        {"structure": structure at end of run,
        "electronic_steps": {All electronic step data in vasprun file},
        "stresses": stress matrix}

    .. attribute:: structures

        List of Structure objects for the structure at each ionic step.

    .. attribute:: tdos

        Total dos calculated at the end of run.

    .. attribute:: idos

        Integrated dos calculated at the end of run.

    .. attribute:: pdos

        List of list of PDos objects. Access as pdos[atomindex][orbitalindex]

    .. attribute:: efermi

        Fermi energy

    .. attribute:: eigenvalues

        Available only if parse_eigen=True. Final eigenvalues as a dict of
        {(spin, kpoint index):[[eigenvalue, occu]]}.
        This representation is based on actual ordering in VASP and is meant as
        an intermediate representation to be converted into proper objects. The
        kpoint index is 0-based (unlike the 1-based indexing in VASP).

    .. attribute:: projected_eigenvalues

        Final projected eigenvalues as a dict of
        {(atom index, band index, kpoint index, Orbital, Spin):float}
        This representation is based on actual ordering in VASP and is meant as
        an intermediate representation to be converted into proper objects. The
        kpoint, band and atom indices are 0-based (unlike the 1-based indexing
        in VASP).

    .. attribute:: dielectric
        The real and imaginary part of the dielectric constant (e.g., computed
        by RPA) in function of the energy (frequency). Optical properties (e.g.
        absorption coefficient) can be obtained through this.
        The date is given as a tuple of 3 values containing each of them
        the energy, the real part tensor, and the imaginary part tensor
        ([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
        real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,
        imag_partxy, imag_partyz, imag_partxz]])

    **Vasp inputs**

    .. attribute:: incar

        Incar object for parameters specified in INCAR file.

    .. attribute:: parameters

        Incar object with parameters that vasp actually used, including all
        defaults.

    .. attribute:: kpoints

        Kpoints object for KPOINTS specified in run.

    .. attribute:: actual_kpoints

        List of actual kpoints, e.g.,
        [[0.25, 0.125, 0.08333333], [-0.25, 0.125, 0.08333333],
        [0.25, 0.375, 0.08333333], ....]

    .. attribute:: actual_kpoints_weights

        List of kpoint weights, E.g.,
        [0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, ....]

    .. attribute:: atomic_symbols

        List of atomic symbols, e.g., ["Li", "Fe", "Fe", "P", "P", "P"]

    .. attribute:: potcar_symbols

        List of POTCAR symbols. e.g.,
        ["PAW_PBE Li 17Jan2003", "PAW_PBE Fe 06Sep2000", ..]

    Author: Shyue Ping Ong
    """
    supported_properties = ["lattice_rec", "vasp_version", "incar",
                            "parameters", "potcar_symbols", "atomic_symbols",
                            "kpoints", "actual_kpoints", "structures",
                            "actual_kpoints_weights", "dos_energies",
                            "eigenvalues", "tdos", "idos", "pdos", "efermi",
                            "ionic_steps", "dos_has_errors",
                            "projected_eigenvalues", "dielectric"]

    def __init__(self, filename, ionic_step_skip=None, parse_dos=True,
                 parse_eigen=True, parse_projected_eigen=False):
        """
        Args:
            filename:
                Filename to parse
            ionic_step_skip:
                If ionic_step_skip is a number > 1, only every ionic_step_skip
                ionic steps will be read for structure and energies. This is
                very useful if you are parsing very large vasprun.xml files and
                you are not interested in every single ionic step. Note that
                the initial and final structure of all runs will always be
                read, regardless of the ionic_step_skip.
            parse_dos:
                Whether to parse the dos. Defaults to True. Set
                to False to shave off significant time from the parsing if you
                are not interested in getting those data.
            parse_eigen:
                Whether to parse the eigenvalues. Defaults to True. Set
                to False to shave off significant time from the parsing if you
                are not interested in getting those data.
            parse_projected_eigen:
                Whether to parse the projected eigenvalues. Defaults to False.
                Set to True to obtain projected eigenvalues. **Note that this
                can take an extreme amount of time and memory.** So use this
                wisely.
        """
        self.filename = filename
        self.ionic_step_skip = ionic_step_skip

        with zopen(filename) as f:
            handler = VasprunHandler(
                filename, parse_dos=parse_dos,
                parse_eigen=parse_eigen,
                parse_projected_eigen=parse_projected_eigen
            )
            if ionic_step_skip is None:
                xml.sax.parse(f, handler)
            else:
                #remove parts of the xml file and parse the string
                run = f.read()
                steps = run.split("<calculation>")
                new_steps = steps[::int(ionic_step_skip)]
                #add the last step from the run
                if steps[-1] != new_steps[-1]:
                    new_steps.append(steps[-1])
                xml.sax.parseString("<calculation>".join(new_steps),
                                    handler)
            for k in Vasprun.supported_properties:
                setattr(self, k, getattr(handler, k))

    @property
    def converged(self):
        """
        True if a relaxation run is converged. Checking is performed on both
        the final electronic convergence as well as whether the number of
        ionic steps is equal to the NSW setting.
        """
        if len(self.ionic_steps[-1]["electronic_steps"]) == \
                self.parameters["NELM"]:
            return False
        nsw = self.parameters.get("NSW", 0)
        return len(self.ionic_steps) < nsw or nsw < 1

    @property
    @unitized("eV")
    def final_energy(self):
        """
        Final energy from the vasp run.
        """
        return self.ionic_steps[-1]["electronic_steps"][-1]["e_wo_entrp"]

    @property
    def final_structure(self):
        """
        Final structure from vasprun.
        """
        return self.structures[-1]

    @property
    def initial_structure(self):
        """
        Initial structure from vasprun.
        """
        return self.structures[0]

    @property
    def complete_dos(self):
        """
        A complete dos object which incorporates the total dos and all
        projected dos.
        """
        final_struct = self.final_structure
        pdoss = {final_struct[i]: pdos for i, pdos in enumerate(self.pdos)}
        return CompleteDos(self.final_structure, self.tdos, pdoss)

    @property
    def hubbards(self):
        """
        Hubbard U values used if a vasprun is a GGA+U run. {} otherwise.
        """
        symbols = [s.split()[1] for s in self.potcar_symbols]
        symbols = [re.split("_", s)[0] for s in symbols]
        if not self.incar.get("LDAU", False):
            return {}
        us = self.incar.get("LDAUU", self.parameters.get("LDAUU"))
        js = self.incar.get("LDAUJ", self.parameters.get("LDAUJ"))
        if len(us) == len(symbols):
            return {symbols[i]: us[i] - js[i] for i in xrange(len(symbols))}
        elif sum(us) == 0 and sum(js) == 0:
            return {}
        else:
            raise VaspParserError("Length of U value parameters and atomic "
                                  "symbols are mismatched")

    @property
    def run_type(self):
        """
        Returns the run type. Currently supports only GGA and HF calcs.

        TODO: Fix for other functional types like LDA, PW91, etc.
        """
        if self.is_hubbard:
            return "GGA+U"
        elif self.parameters.get("LHFCALC", False):
            return "HF"
        else:
            return "GGA"

    @property
    def is_hubbard(self):
        """
        True if run is a DFT+U run.
        """
        if len(self.hubbards) == 0:
            return False
        return sum(self.hubbards.values()) > 1e-8

    @property
    def is_spin(self):
        """
        True if run is spin-polarized.
        """
        return True if self.parameters.get("ISPIN", 1) == 2 else False

    def get_band_structure(self, kpoints_filename=None, efermi=None,
                           line_mode=False):
        """
        Returns the band structure as a BandStructure object

        Args:
            kpoints_filename:
                Full path of the KPOINTS file from which the band structure is
                generated.
                If none is provided, the code will try to intelligently
                determine the appropriate KPOINTS file by substituting the
                filename of the vasprun.xml with KPOINTS.
                The latter is the default behavior.

            efermi:
                If you want to specify manually the fermi energy this is where
                you should do it. By default, the None value means the code
                will get it from the vasprun.

            line_mode:
                force the band structure to be considered as a run along
                symmetry lines

        Returns:
            a BandStructure object (or more specifically a
            BandStructureSymmLine object if the run is detected to be a run
            along symmetry lines)

            Two types of runs along symmetry lines are accepted: non-sc with
            Line-Mode in the KPOINT file or hybrid, self-consistent with a
            uniform grid+a few kpoints along symmetry lines (explicit KPOINTS
            file) (it's not possible to run a non-sc band structure with hybrid
            functionals). The explicit KPOINTS file needs to have data on the
            kpoint label as commentary.

        TODO:
            - make a bit more general for non Symm Line band structures
            - make a decision on the convention with 2*pi or not.
        """
        if not kpoints_filename:
            kpoints_filename = self.filename.replace('vasprun.xml', 'KPOINTS')
        if not os.path.exists(kpoints_filename):
            raise VaspParserError('KPOINTS needed to obtain band structure.')

        if efermi is None:
            efermi = self.efermi

        kpoint_file = Kpoints.from_file(kpoints_filename)
        lattice_new = Lattice(self.lattice_rec.matrix * 2 * math.pi)
        #lattice_rec=[self.lattice_rec.matrix[i][j] for i,j in range(3)]

        kpoints = [np.array(self.actual_kpoints[i])
                   for i in range(len(self.actual_kpoints))]
        dict_eigen = self.to_dict['output']['eigenvalues']
        dict_p_eigen = {}
        if 'projected_eigenvalues' in self.to_dict['output']:
            dict_p_eigen = self.to_dict['output']['projected_eigenvalues']

        p_eigenvals = {}
        if "1" in dict_eigen["0"] and "-1" in dict_eigen["0"]\
                and self.incar['ISPIN'] == 2:
            eigenvals = {Spin.up: [], Spin.down: []}
            if len(dict_p_eigen) != 0:
                p_eigenvals = {Spin.up: [], Spin.down: []}
        else:
            eigenvals = {Spin.up: []}
            if len(dict_p_eigen) != 0:
                p_eigenvals = {Spin.up: []}

        neigenvalues = [len(v['1']) for v in dict_eigen.values()]
        min_eigenvalues = min(neigenvalues)
        for i in range(min_eigenvalues):
            eigenvals[Spin.up].append([dict_eigen[str(j)]['1'][i][0]
                                       for j in range(len(kpoints))])
            if len(dict_p_eigen) != 0:
                p_eigenvals[Spin.up].append(
                    [{Orbital.from_string(orb):
                      dict_p_eigen[j]['1'][i][orb]
                      for orb in dict_p_eigen[j]['1'][i]}
                     for j in range(len(kpoints))])
        if Spin.down in eigenvals:
            for i in range(min_eigenvalues):
                eigenvals[Spin.down].append([dict_eigen[str(j)]['-1'][i][0]
                                             for j in range(len(kpoints))])
                if len(dict_p_eigen) != 0:
                    p_eigenvals[Spin.down].append(
                        [{Orbital.from_string(orb):
                          dict_p_eigen[j]['-1'][i][orb]
                          for orb in dict_p_eigen[j]['-1'][i]}
                         for j in range(len(kpoints))]
                    )

        #check if we have an hybrid band structure computation
        #for this we look at the presence of the LHFCALC tag and of k-points
        #that have weights=0.0
        hybrid_band = False
        if self.parameters['LHFCALC']:
            for l in kpoint_file.labels:
                if l is not None:
                    hybrid_band = True

        if kpoint_file.style == "Line_mode" or hybrid_band or line_mode:
            labels_dict = {}
            if hybrid_band:
                start_bs_index = 0
                for i in range(len(self.actual_kpoints)):
                    if self.actual_kpoints_weights[i] == 0.0:
                        start_bs_index = i
                        break
                for i in range(len(kpoint_file.kpts)):
                    if kpoint_file.labels[i] is not None:
                        labels_dict[kpoint_file.labels[i]] = \
                            kpoint_file.kpts[i]
                #remake the data only considering line band structure k-points
                #(weight = 0.0 kpoints)
                kpoints = kpoints[start_bs_index:len(kpoints)]
                up_eigen = [eigenvals[Spin.up][i][start_bs_index:
                                                  len(eigenvals[Spin.up][i])]
                            for i in range(len(eigenvals[Spin.up]))]
                if self.is_spin:
                    down_eigen = [eigenvals[Spin.down][i]
                                  [start_bs_index:
                                   len(eigenvals[Spin.down][i])]
                                  for i in range(len(eigenvals[Spin.down]))]
                    eigenvals = {Spin.up: up_eigen,
                                 Spin.down: down_eigen}
                else:
                    eigenvals = {Spin.up: up_eigen}
            else:
                labels_dict = dict(zip(kpoint_file.labels, kpoint_file.kpts))
                labels_dict.pop(None, None)
            return BandStructureSymmLine(kpoints, eigenvals, lattice_new,
                                         efermi, labels_dict,
                                         structure=self.final_structure,
                                         projections=p_eigenvals)
        else:
            return BandStructure(kpoints, eigenvals, lattice_new, efermi,
                                 structure=self.final_structure,
                                 projections=p_eigenvals)

    @property
    def eigenvalue_band_properties(self):
        """
        Band properties from the eigenvalues as a tuple,
        (band gap, cbm, vbm, is_band_gap_direct).
        """
        vbm = -float("inf")
        vbm_kpoint = None
        cbm = float("inf")
        cbm_kpoint = None
        for k, val in self.eigenvalues.items():
            for (eigenval, occu) in val:
                if occu > 1e-8 and eigenval > vbm:
                    vbm = eigenval
                    vbm_kpoint = k[0]
                elif occu <= 1e-8 and eigenval < cbm:
                    cbm = eigenval
                    cbm_kpoint = k[0]
        return max(cbm - vbm, 0), cbm, vbm, vbm_kpoint == cbm_kpoint

    @property
    def to_dict(self):
        """
        Json-serializable dict representation.
        """
        d = {"vasp_version": self.vasp_version,
             "has_vasp_completed": self.converged,
             "nsites": len(self.final_structure)}
        comp = self.final_structure.composition
        d["unit_cell_formula"] = comp.to_dict
        d["reduced_cell_formula"] = Composition(comp.reduced_formula).to_dict
        d["pretty_formula"] = comp.reduced_formula
        symbols = [s.split()[1] for s in self.potcar_symbols]
        symbols = [re.split("_", s)[0] for s in symbols]
        d["is_hubbard"] = self.incar.get("LDAU", False)
        if d["is_hubbard"]:
            us = self.incar.get("LDAUU", self.parameters.get("LDAUU"))
            js = self.incar.get("LDAUJ", self.parameters.get("LDAUJ"))
            if len(us) == len(symbols):
                d["hubbards"] = {symbols[i]: us[i] - js[i]
                                 for i in xrange(len(symbols))}
            elif sum(us) == 0 and sum(js) == 0:
                d["is_hubbard"] = False
                d["hubbards"] = {}
            else:
                raise VaspParserError("Length of U value parameters and atomic"
                                      " symbols are mismatched.")
        else:
            d["hubbards"] = {}

        unique_symbols = sorted(list(set(symbols)))
        d["elements"] = unique_symbols
        d["nelements"] = len(unique_symbols)

        d["run_type"] = self.run_type

        vin = {"incar": {k: v for k, v in self.incar.items()},
               "crystal": self.initial_structure.to_dict,
               "kpoints": self.kpoints.to_dict}
        actual_kpts = [{"abc":list(self.actual_kpoints[i]),
                        "weight":self.actual_kpoints_weights[i]}
                       for i in xrange(len(self.actual_kpoints))]
        vin["kpoints"]["actual_points"] = actual_kpts
        vin["potcar"] = [s.split(" ")[1] for s in self.potcar_symbols]
        vin["potcar_type"] = [s.split(" ")[0] for s in self.potcar_symbols]
        vin["parameters"] = {k: v for k, v in self.parameters.items()}
        vin["lattice_rec"] = self.lattice_rec.to_dict
        d["input"] = vin

        nsites = len(self.final_structure)

        vout = {"ionic_steps": self.ionic_steps,
                "final_energy": self.final_energy,
                "final_energy_per_atom": self.final_energy / nsites,
                "crystal": self.final_structure.to_dict,
                "efermi": self.efermi}

        eigen = defaultdict(dict)
        for (spin, index), values in self.eigenvalues.items():
            eigen[index][str(spin)] = values
        vout["eigenvalues"] = eigen
        vout['dielectric'] = self.dielectric

        peigen = []

        if len(self.projected_eigenvalues) != 0:
            for i in range(len(eigen)):
                peigen.append({})
                for spin in eigen[i].keys():
                    peigen[i][spin] = []
                    for j in range(len(eigen[i][spin])):
                        peigen[i][spin].append({})
            for (spin, kpoint_index, band_index, ion_index, orbital), value \
                    in self.projected_eigenvalues.items():
                beigen = peigen[kpoint_index][str(spin)][band_index]
                if orbital not in beigen:
                    beigen[orbital] = [0.0] * nsites
                beigen[orbital][ion_index] = value
        vout['projected_eigenvalues'] = peigen
        (gap, cbm, vbm, is_direct) = self.eigenvalue_band_properties
        vout.update(dict(bandgap=gap, cbm=cbm, vbm=vbm,
                         is_gap_direct=is_direct))
        d['output'] = vout

        return clean_json(d, strict=True)


class VasprunHandler(xml.sax.handler.ContentHandler):
    """
    Sax handler for vasprun.xml. Attributes are mirrored into Vasprun object.
    Generally should not be initiatized on its own.
    """

    def __init__(self, filename, parse_dos=True, parse_eigen=True,
                 parse_projected_eigen=False):
        self.filename = filename
        self.parse_dos = parse_dos
        self.parse_eigen = parse_eigen
        self.parse_projected_eigen = parse_projected_eigen

        self.step_count = 0
        # variables to be filled
        self.vasp_version = None
        self.incar = Incar()
        self.parameters = Incar()
        self.potcar_symbols = []
        self.atomic_symbols = []
        self.kpoints = Kpoints()
        self.actual_kpoints = []
        self.actual_kpoints_weights = []
        self.dos_energies = None

        #  will  be  {(spin, kpoint index): [[energy, occu]]}
        self.eigenvalues = {}

        #{(spin, kpoint_index, band_index, atom_ind, orb):float}
        self.projected_eigenvalues = {}

        self.tdos = {}
        self.idos = {}
        self.pdos = {}
        self.efermi = None
        self.ionic_steps = []  # should be a list of dict
        self.structures = []
        self.lattice_rec = []
        self.stress = []
        self.dielectric = ([], [], [])

        self.input_read = False
        self.read_structure = False
        self.read_rec_lattice = False
        self.read_calculation = False
        self.read_eigen = False
        self.read_projected_eigen = False
        self.read_diel = False
        self.read_dos = False
        self.in_efermi = False
        self.read_atoms = False
        self.read_lattice = False
        self.read_positions = False
        self.incar_param = None

        #Intermediate variables
        self.dos_energies_val = []
        self.dos_val = []
        self.idos_val = []
        self.raw_data = []
        self.diel_energies_val = []
        self.diel_val = []

        #will be set to true if there is an error parsing the Dos.
        self.dos_has_errors = False

        # The internal _state attribute stores the current state of parsing.
        # It is essentially a dict of {str: boolean}. The handler sets a
        # particular key (e.g., "calculation") to True when the start
        # element is encountered and False when the endElement is found.
        self._state = defaultdict(bool)

    def startElement(self, name, attributes):
        """
        Detect the various starting elements and call the appropriate init
        fucntions. The parsing is broken down into two basic blocks for
        clarity:

        1. The input stage, where all the run parameters are read. Note that
           the vasprun.xml contains far more information about run
           parameters than present in the INCAR.
        2. The calc stage, where all the ionic steps are read with the
           structure at each stage, with the DOS, bandstructure, etc.
        """
        self._state[name] = attributes.get("name", True)
        self.read_val = False

        #Nested if loops makes reading much faster.
        if not self.input_read:  # reading input parameters
            self._init_input(name, attributes)
        else:  # reading calculations and structures and eigenvalues.
            self._init_calc(name, attributes)
        if self.read_val:
            self.val = StringIO.StringIO()

    def characters(self, data):
        """
        Logic to read in characters when the parser is in certain states.
        """
        if self.read_val:
            self.val.write(data)
        if self.read_lattice:
            self.latticestr.write(data)
        elif self.read_positions:
            self.posstr.write(data)
        elif self.read_rec_lattice:
            self.latticerec.write(data)

    def endElement(self, name):
        """
        Calls the appropriate internal object creation methods when a
        closing element is encountered. For example, a "</structure>" will
        create a Structure object and append it to the list of structures.
        """
        if not self.input_read:
            self._read_input(name)
        else:
            if self.read_structure:
                self._read_structure(name)
            elif self.read_diel:
                self._read_diel(name)
            elif self.read_dos:
                self._read_dos(name)
            elif self.read_eigen:
                self._read_eigen(name)
            elif self.read_projected_eigen:
                self._read_projected_eigen(name)
            elif self.read_calculation:
                self._read_calc(name)
        self._state[name] = False

    def _init_input(self, name, attributes):
        state = self._state
        if (name == "i" or name == "v") and \
                (state["incar"] or state["parameters"]):
            self.incar_param = attributes["name"]
            self.param_type = "float" if "type" not in attributes \
                else attributes["type"]
            self.read_val = True
        elif name == "v" and state["kpoints"]:
            self.read_val = True
        elif name == "generation" and state["kpoints"]:
            self.kpoints.comment = "Kpoints from vasprun.xml"
            self.kpoints.num_kpts = 0
            self.kpoints.style = attributes["param"]
            self.kpoints.kpts = []
            self.kpoints.kpts_shift = [0, 0, 0]
        elif name == "c" and \
                (state["array"] == "atoms" or
                 state["array"] == "atomtypes"):
            self.read_val = True
        elif name == "i" and state["i"] == "version" and state["generator"]:
            self.read_val = True

    def _init_calc(self, name, attributes):
        state = self._state
        if self.read_structure and name == "v":
            if state["varray"] == "basis":
                self.read_lattice = True
            elif state["varray"] == "positions":
                self.read_positions = True
            elif state["varray"] == "rec_basis":
                self.read_rec_lattice = True
        elif self.read_calculation:
            if name == "i" and state["scstep"]:
                logger.debug("Reading scstep...")
                self.read_val = True
            elif name == "v" and (state["varray"] == "forces" or
                                  state["varray"] == "stress"):
                self.read_positions = True
            elif name == "dielectricfunction":
                logger.debug("Reading dielectric function...")
                self.read_diel = True
            elif name == "dos" and self.parse_dos:
                logger.debug("Reading dos...")
                self.dos_energies = None
                self.tdos = {}
                self.idos = {}
                self.pdos = {}
                self.efermi = None
                self.read_dos = True
            elif name == "eigenvalues" and self.parse_eigen and \
                    (not state["projected"]):
                logger.debug("Reading eigenvalues. Projected = {}"
                             .format(state["projected"]))
                self.eigenvalues = {}
                self.read_eigen = True
            elif name == "eigenvalues" and self.parse_projected_eigen and \
                    state["projected"]:
                logger.debug("Reading projected eigenvalues...")
                self.projected_eigen = {}
                self.read_projected_eigen = True
            elif self.read_eigen or self.read_projected_eigen:
                if name == "r" and state["set"]:
                    self.read_val = True
                elif name == "set" and "comment" in attributes:
                    comment = attributes["comment"]
                    state["set"] = comment
                    if comment.startswith("spin"):
                        self.eigen_spin = Spin.up \
                            if state["set"] in ["spin 1", "spin1"] \
                            else Spin.down
                        logger.debug("Reading spin {}".format(self.eigen_spin))
                    elif comment.startswith("kpoint"):
                        self.eigen_kpoint = int(comment.split(" ")[1])
                        logger.debug("Reading kpoint {}"
                                     .format(self.eigen_kpoint))
                    elif comment.startswith("band"):
                        self.eigen_band = int(comment.split(" ")[1])
                        logger.debug("Reading band {}"
                                     .format(self.eigen_band))
            elif self.read_diel:
                if name == "r" and (state["imag"] or state["real"]):
                    self.read_val = True
            elif self.read_dos:
                if (name == "i" and state["i"] == "efermi") or \
                   (name == "r" and state["set"]):
                    self.read_val = True
                elif name == "set" and "comment" in attributes:
                    comment = attributes["comment"]
                    state["set"] = comment
                    if state["partial"]:
                        if comment.startswith("ion"):
                            self.pdos_ion = int(comment.split(" ")[1])
                        elif comment.startswith("spin"):
                            self.pdos_spin = Spin.up \
                                if state["set"] in ["spin 1", "spin1"] \
                                else Spin.down

        if name == "calculation":
            self.step_count += 1
            self.scdata = []
            self.read_calculation = True
        elif name == "scstep":
            self.scstep = {}
        elif name == "structure":
            self.latticestr = StringIO.StringIO()
            self.latticerec = StringIO.StringIO()
            self.posstr = StringIO.StringIO()
            self.read_structure = True
        elif name == "varray" and (state["varray"] in ["forces", "stress"]):
            self.posstr = StringIO.StringIO()

    def _read_input(self, name):
        state = self._state
        if name == "i":
            if state["incar"]:
                self.incar[self.incar_param] = \
                    parse_parameters(self.param_type,
                                     self.val.getvalue().strip())
            elif state["parameters"]:
                self.parameters[self.incar_param] = \
                    parse_parameters(self.param_type,
                                     self.val.getvalue().strip())
            elif state["generator"] and state["i"] == "version":
                self.vasp_version = self.val.getvalue().strip()
            self.incar_param = None
        elif name == "set":
            if state["array"] == "atoms":
                self.atomic_symbols = self.atomic_symbols[::2]
                self.atomic_symbols = [sym if sym != "X" else "Xe"
                                       for sym in self.atomic_symbols]
            elif state["array"] == "atomtypes":
                self.potcar_symbols = self.potcar_symbols[4::5]
                self.input_read = True
        elif name == "c":
            if state["array"] == "atoms":
                self.atomic_symbols.append(self.val.getvalue().strip())
            elif state["array"] == "atomtypes":
                self.potcar_symbols.append(self.val.getvalue().strip())
        elif name == "v":
            if state["incar"]:
                self.incar[self.incar_param] = \
                    parse_v_parameters(self.param_type,
                                       self.val.getvalue().strip(),
                                       self.filename, self.incar_param)
                self.incar_param = None
            elif state["parameters"]:
                self.parameters[self.incar_param] = \
                    parse_v_parameters(self.param_type,
                                       self.val.getvalue().strip(),
                                       self.filename, self.incar_param)
            elif state["kpoints"]:
                if state["varray"] == "kpointlist":
                    self.actual_kpoints.append(map(float,
                                               self.val.getvalue().split()))
                if state["varray"] == "weights":
                    val = float(self.val.getvalue())
                    self.actual_kpoints_weights.append(val)
                if state["v"] == "divisions":
                    self.kpoints.kpts = [map(int, self.val.getvalue().split())]
                elif state["v"] == "usershift":
                    self.kpoints.kpts_shift = map(float,
                                                  self.val.getvalue().split())
                elif state["v"] == "genvec1" or state["v"] == "genvec2" or \
                        state["v"] == "genvec3" or state["v"] == "shift":
                    setattr(self.kpoints, state["v"],
                            map(float, self.val.getvalue().split()))

    def _read_calc(self, name):
        state = self._state
        if name == "i" and state["scstep"]:
            self.scstep[state["i"]] = float(self.val.getvalue())
        elif name == "scstep":
            self.scdata.append(self.scstep)
            logger.debug("Finished reading scstep...")
        elif name == "varray" and state["varray"] == "forces":
            self.forces = np.array(map(float,
                                       self.posstr.getvalue().split()))
            self.forces.shape = (len(self.atomic_symbols), 3)
            self.read_positions = False
        elif name == "varray" and state["varray"] == "stress":
            self.stress = np.array(map(float,
                                       self.posstr.getvalue().split()))
            self.stress.shape = (3, 3)
            self.read_positions = False
        elif name == "calculation":
            self.ionic_steps.append({"electronic_steps": self.scdata,
                                     "structure": self.structures[-1],
                                     "forces": self.forces,
                                     "stress": self.stress})
            self.read_calculation = False

    def _read_structure(self, name):
        if name == "v":
            self.read_positions = False
            self.read_lattice = False
            self.read_rec_lattice = False
        elif name == "structure":
            self.lattice = map(float, self.latticestr.getvalue().split())
            self.pos = np.array(map(float,
                                    self.posstr.getvalue().split()))
            self.pos.shape = (len(self.atomic_symbols), 3)
            self.structures.append(Structure(self.lattice, self.atomic_symbols,
                                             self.pos))
            self.lattice_rec = Lattice(map(float,
                                       self.latticerec.getvalue().split()))
            self.read_structure = False
            self.read_positions = False
            self.read_lattice = False
            self.read_rec_lattice = False

    def _read_diel(self, name):
        state = self._state
        if name == "r" and (state["imag"] or state["real"]):
            tok = self.val.getvalue().split()
            self.diel_energies_val.append(float(tok[0]))
            self.diel_val.append([float(tok[i])
                                  for i in range(1, 7)])
        elif name == "set":
            if state["imag"]:
                self.dielectric[0].extend(self.diel_energies_val)
                self.dielectric[2].extend(self.diel_val)
                self.diel_energies_val = []
                self.diel_val = []
            if state["real"]:
                self.dielectric[1].extend(self.diel_val)
                self.diel_energies_val = []
                self.diel_val = []
        elif name == "dielectricfunction":
            self.read_diel = False

    def _read_dos(self, name):
        state = self._state
        try:
            if name == "i" and state["i"] == "efermi":
                self.efermi = float(self.val.getvalue().strip())
            elif name == "r" and state["total"] and \
                    str(state["set"]).startswith("spin"):
                tok = self.val.getvalue().split()
                self.dos_energies_val.append(float(tok[0]))
                self.dos_val.append(float(tok[1]))
                self.idos_val.append(float(tok[2]))
            elif name == "r" and state["partial"] and \
                    str(state["set"]).startswith("spin"):
                tok = self.val.getvalue().split()
                self.raw_data.append(map(float, tok[1:]))
            elif name == "set":
                if state["total"] and str(state["set"]).startswith("spin"):
                    spin = Spin.up if state["set"] == "spin 1" else Spin.down
                    self.tdos[spin] = self.dos_val
                    self.idos[spin] = self.idos_val
                    self.dos_energies = self.dos_energies_val
                    self.dos_energies_val = []
                    self.dos_val = []
                    self.idos_val = []
                elif state["partial"] and str(state["set"]).startswith("spin"):
                    spin = Spin.up if state["set"] == "spin 1" else Spin.down
                    self.norbitals = len(self.raw_data[0])
                    for i in xrange(self.norbitals):
                        self.pdos[(self.pdos_ion, i, spin)] = \
                            [row[i] for row in self.raw_data]
                    self.raw_data = []
            elif name == "partial":
                all_pdos = []
                natom = len(self.atomic_symbols)
                for iatom in xrange(1, natom + 1):
                    all_pdos.append(defaultdict())
                    for iorbital in xrange(self.norbitals):
                        updos = self.pdos[(iatom, iorbital, Spin.up)]
                        downdos = self.pdos.get((iatom, iorbital, Spin.down),
                                                None)
                        orb = Orbital.from_vasp_index(iorbital)
                        if downdos:
                            all_pdos[-1][orb] = {Spin.up: updos,
                                                 Spin.down: downdos}
                        else:
                            all_pdos[-1][orb] = {Spin.up: updos}
                self.pdos = all_pdos
            elif name == "total":
                self.tdos = Dos(self.efermi, self.dos_energies, self.tdos)
                self.idos = Dos(self.efermi, self.dos_energies, self.idos)
            elif name == "dos":
                self.read_dos = False
        except:
            self.dos_has_errors = True

    def _read_eigen(self, name):
        state = self._state
        if name == "r" and str(state["set"]).startswith("kpoint"):
            tok = self.val.getvalue().split()
            self.raw_data.append(map(float, tok))
        elif name == "set" and str(state["set"]).startswith("kpoint"):
            self.eigenvalues[(self.eigen_spin, self.eigen_kpoint - 1)] = \
                self.raw_data
            self.raw_data = []
        elif name == "eigenvalues":
            logger.debug("Finished reading eigenvalues. "
                         "No. eigen = {}".format(len(self.eigenvalues)))
            self.read_eigen = False

    def _read_projected_eigen(self, name):
        state = self._state
        if name == "r" and str(state["set"]).startswith("band"):
            tok = self.val.getvalue().split()
            self.raw_data.append({Orbital.from_vasp_index(i): float(val)
                                  for i, val in enumerate(tok)})
        elif name == "set" and str(state["set"]).startswith("band"):
            logger.debug("Processing projected eigenvalues for " +
                         "band {}, kpoint {}, spin {}."
                         .format(self.eigen_band - 1, self.eigen_kpoint - 1,
                                 self.eigen_spin))
            for atom_ind, data in enumerate(self.raw_data):
                for orb, val in data.items():
                    self.projected_eigenvalues[(self.eigen_spin,
                                                self.eigen_kpoint - 1,
                                                self.eigen_band - 1,
                                                atom_ind, orb)] = val
            self.raw_data = []
        elif name == "projected":
            logger.debug("Finished reading projected eigenvalues. "
                         "No. eigen = {}".format(len(self.eigenvalues)))
            self.read_projected_eigen = False


def parse_parameters(val_type, val):
    """
    Helper function to convert a Vasprun parameter into the proper type.
    Boolean, int and float types are converted.

    Args:
        val_type : Value type parsed from vasprun.xml.
        val : Actual string value parsed for vasprun.xml.
    """
    if val_type == "logical":
        return val == "T"
    elif val_type == "int":
        return int(val)
    elif val_type == "string":
        return val.strip()
    else:
        return float(val)


def parse_v_parameters(val_type, val, filename, param_name):
    """
    Helper function to convert a Vasprun array-type parameter into the proper
    type. Boolean, int and float types are converted.

    Args:
        val_type:
            Value type parsed from vasprun.xml.
        val:
            Actual string value parsed for vasprun.xml.
        filename:
            Fullpath of vasprun.xml. Used for robust error handling.  E.g.,
            if vasprun.xml contains \*\*\* for some Incar parameters, the code
            will try to read from an INCAR file present in the same directory.
        param_name:
            Name of parameter.

    Returns:
        Parsed value.
    """
    if val_type == "logical":
        val = map(lambda i: i == "T", val.split())
    elif val_type == "int":
        try:
            val = map(int, val.split())
        except ValueError:
            # Fix for stupid error in vasprun sometimes which displays
            # LDAUL/J as 2****
            val = parse_from_incar(filename, param_name)
            if val is None:
                raise IOError("Error in parsing vasprun.xml")
    elif val_type == "string":
        val = val.split()
    else:
        try:
            val = map(float, val.split())
        except ValueError:
            # Fix for stupid error in vasprun sometimes which displays
            # MAGMOM as 2****
            val = parse_from_incar(filename, param_name)
            if val is None:
                raise IOError("Error in parsing vasprun.xml")
    return val


def parse_from_incar(filename, key):
    """
    Helper function to parse a parameter from the INCAR.
    """
    dirname = os.path.dirname(filename)
    for f in os.listdir(dirname):
        if re.search("INCAR", f):
            warnings.warn("INCAR found. Using " + key + " from INCAR.")
            incar = Incar.from_file(os.path.join(dirname, f))
            if key in incar:
                return incar[key]
            else:
                return None
    return None


class Outcar(object):
    """
    Parser for data in OUTCAR that is not available in Vasprun.xml

    Note, this class works a bit differently than most of the other
    VaspObjects, since the OUTCAR can be very different depending on which
    "type of run" performed.

    Creating the OUTCAR class with a filename reads "regular parameters" that
    are always present.

    .. attribute:: magnetization

        Magnetization on each ion as a tuple of dict, e.g.,
        ({"d": 0.0, "p": 0.003, "s": 0.002, "tot": 0.005}, ... )
        Note that this data is not always present.  LORBIT must be set to some
        other value than the default.

    .. attribute:: charge

        Charge on each ion as a tuple of dict, e.g.,
        ({"p": 0.154, "s": 0.078, "d": 0.0, "tot": 0.232}, ...)
        Note that this data is not always present.  LORBIT must be set to some
        other value than the default.

    .. attribute:: is_stopped

        True if OUTCAR is from a stopped run (using STOPCAR, see Vasp Manual).

    .. attribute:: run_stats

        Various useful run stats as a dict including "System time (sec)",
        "Total CPU time used (sec)", "Elapsed time (sec)",
        "Maximum memory used (kb)", "Average memory used (kb)",
        "User time (sec)".

    One can then call a specific reader depending on the type of run being
    performed. These are currently: read_igpar(), read_lepsilon() and
    read_lcalcpol().

    See the documentation of those methods for more documentation.

    Authors: Rickard Armiento, Shyue Ping Ong
    """
    def __init__(self, filename):
        self.filename = filename
        self.is_stopped = False
        with zopen(filename, "r") as f:
            read_charge_mag = False
            charge = []
            mag = []
            header = []
            run_stats = {}
            total_mag = None
            nelect = None
            efermi = None

            time_patt = re.compile("\((sec|kb)\)")
            efermi_patt = re.compile("E-fermi\s*:\s*(\S+)")
            nelect_patt = re.compile("number of electron\s+(\S+)\s+"
                                     "magnetization\s+(\S+)")
            all_lines = []
            for line in reverse_readline(f):
                clean = line.strip()
                all_lines.append(clean)
                if clean.startswith("tot ") and not (charge and mag):
                    read_charge_mag = True
                    data = []
                elif read_charge_mag:
                    if clean.startswith("# of ion"):
                        header = re.split("\s{2,}", line.strip())
                        header.pop(0)
                    elif clean == "total charge":
                        data.reverse()
                        charge = [dict(zip(header, v)) for v in data]
                        read_charge_mag = False
                    elif clean == "magnetization (x)":
                        data.reverse()
                        mag = [dict(zip(header, v)) for v in data]
                        read_charge_mag = False
                    else:
                        m = re.match("\s*(\d+)\s+(([\d\.\-]+)\s+)+", clean)
                        if m:
                            toks = [float(i) for i in re.findall("[\d\.\-]+",
                                                                 clean)]
                            toks.pop(0)
                            data.append(toks)
                elif clean.find("soft stop encountered!  aborting job") != -1:
                    self.is_stopped = True
                else:
                    if time_patt.search(line):
                        tok = line.strip().split(":")
                        run_stats[tok[0].strip()] = float(tok[1].strip())
                        continue
                    m = efermi_patt.search(clean)
                    if m:
                        try:
                            #try-catch because VASP sometimes prints
                            #'E-fermi: ********     XC(G=0):  -6.1327
                            #alpha+bet : -1.8238'
                            efermi = float(m.group(1))
                            continue
                        except ValueError:
                            efermi = None
                            continue
                    m = nelect_patt.search(clean)
                    if m:
                        nelect = float(m.group(1))
                        total_mag = float(m.group(2))
                if all([nelect, total_mag is not None, efermi is not None,
                        run_stats]):
                    break

            self.run_stats = run_stats
            self.magnetization = tuple(mag)
            self.charge = tuple(charge)
            self.efermi = efermi
            self.nelect = nelect
            self.total_mag = total_mag

    def read_igpar(self):
        """
        Renders accessible:
            er_ev = e<r>_ev (dictionary with Spin.up/Spin.down as keys)
            er_bp = e<r>_bp (dictionary with Spin.up/Spin.down as keys)
            er_ev_tot = spin up + spin down summed
            er_bp_tot = spin up + spin down summed
            p_elc = spin up + spin down summed
            p_ion = spin up + spin down summed

        (See VASP section "LBERRY,  IGPAR,  NPPSTR,  DIPOL tags" for info on
        what these are).
        """

        # variables to be filled
        self.er_ev = {}  # will  be  dict (Spin.up/down) of array(3*float)
        self.er_bp = {}  # will  be  dics (Spin.up/down) of array(3*float)
        self.er_ev_tot = None  # will be array(3*float)
        self.er_bp_tot = None  # will be array(3*float)
        self.p_elec = None
        self.p_ion = None
        try:
            search = []

            # Nonspin cases
            def er_ev(results, match):
                results.er_ev[Spin.up] = np.array(map(float,
                                                      match.groups()[1:4])) / 2
                results.er_ev[Spin.down] = results.er_ev[Spin.up]
                results.context = 2

            search.append(["^ *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           None, er_ev])

            def er_bp(results, match):
                results.er_bp[Spin.up] = np.array([float(match.group(i))
                                                   for i in xrange(1, 4)]) / 2
                results.er_bp[Spin.down] = results.er_bp[Spin.up]

            search.append(["^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           lambda results, line: results.context == 2, er_bp])

            # Spin cases
            def er_ev_up(results, match):
                results.er_ev[Spin.up] = np.array([float(match.group(i))
                                                   for i in xrange(1, 4)])
                results.context = Spin.up

            search.append(["^.*Spin component 1 *e<r>_ev=\( *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                           None, er_ev_up])

            def er_bp_up(results, match):
                results.er_bp[Spin.up] = np.array([float(match.group(1)),
                                                   float(match.group(2)),
                                                   float(match.group(3))])

            search.append(["^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           lambda results,
                           line: results.context == Spin.up, er_bp_up])

            def er_ev_dn(results, match):
                results.er_ev[Spin.down] = np.array([float(match.group(1)),
                                                     float(match.group(2)),
                                                     float(match.group(3))])
                results.context = Spin.down
            search.append(["^.*Spin component 2 *e<r>_ev=\( *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                           None, er_ev_dn])

            def er_bp_dn(results, match):
                results.er_bp[Spin.down] = np.array([float(match.group(i))
                                                     for i in xrange(1, 4)])
            search.append(["^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           lambda results,
                           line: results.context == Spin.down, er_bp_dn])

            # Always present spin/non-spin
            def p_elc(results, match):
                results.p_elc = np.array([float(match.group(i))
                                          for i in xrange(1, 4)])

            search.append(["^.*Total electronic dipole moment: "
                           "*p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)", None, p_elc])

            def p_ion(results, match):
                results.p_ion = np.array([float(match.group(i))
                                          for i in xrange(1, 4)])

            search.append(["^.*ionic dipole moment: "
                           "*p\[ion\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)", None, p_ion])

            self.context = None
            self.er_ev = {Spin.up: None, Spin.down: None}
            self.er_bp = {Spin.up: None, Spin.down: None}

            micro_pyawk(self.filename, search, self)

            if self.er_ev[Spin.up] is not None and \
                    self.er_ev[Spin.down] is not None:
                self.er_ev_tot = self.er_ev[Spin.up] + self.er_ev[Spin.down]

            if self.er_bp[Spin.up] is not None and \
                    self.er_bp[Spin.down] is not None:
                self.er_bp_tot = self.er_bp[Spin.up] + self.er_bp[Spin.down]

        except:
            self.er_ev_tot = None
            self.er_bp_tot = None
            raise Exception("IGPAR OUTCAR could not be parsed.")

    def read_lepsilon(self):
        # variables to be filled
        try:
            search = []

            def dielectric_section_start(results, match):
                results.dielectric_index = -1

            search.append(["MACROSCOPIC STATIC DIELECTRIC TENSOR", None,
                           dielectric_section_start])

            def dielectric_section_start2(results, match):
                results.dielectric_index = 0

            search.append(["-------------------------------------",
                           lambda results,
                           line: results.dielectric_index == -1,
                           dielectric_section_start2])

            def dielectric_data(results, match):
                results.dielectric_tensor[results.dielectric_index, :] = \
                    np.array([float(match.group(i)) for i in xrange(1, 4)])
                results.dielectric_index += 1

            search.append(["^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$",
                           lambda results,
                           line: results.dielectric_index >= 0,
                           dielectric_data])

            def dielectric_section_stop(results, match):
                results.dielectric_index = None

            search.append(["-------------------------------------",
                           lambda results, line: results.dielectric_index >= 1,
                           dielectric_section_stop])

            self.dielectric_index = None
            self.dielectric_tensor = np.zeros((3, 3))

            def piezo_section_start(results, match):
                results.piezo_index = 0

            search.append(["PIEZOELECTRIC TENSOR  for field in x, y, z        "
                           "\(e  Angst\)",
                           None, piezo_section_start])

            def piezo_data(results, match):
                results.piezo_tensor[results.piezo_index, :] = \
                    np.array([float(match.group(i)) for i in xrange(1, 7)])
                results.piezo_index += 1

            search.append(["^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                           " +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                           " +([-0-9.Ee+]+)*$",
                           lambda results, line: results.piezo_index >= 0,
                           piezo_data])

            def piezo_section_stop(results, match):
                results.piezo_index = None

            search.append(["-------------------------------------",
                           lambda results, line: results.piezo_index >= 1,
                           piezo_section_stop])

            self.piezo_index = None
            self.piezo_tensor = np.zeros((3, 6))

            def born_section_start(results, match):
                results.born_ion = -1

            search.append(["BORN EFFECTIVE CHARGES " +
                           "\(in e, cummulative output\)",
                           None, born_section_start])

            def born_ion(results, match):
                results.born_ion = int(match.group(1)) - 1
                results.born[results.born_ion] = np.zeros((3, 3))

            search.append(["ion +([0-9]+)", lambda results,
                           line: results.born_ion is not None, born_ion])

            def born_data(results, match):
                results.born[results.born_ion][int(match.group(1)) - 1, :] = \
                    np.array([float(match.group(i)) for i in xrange(2, 5)])

            search.append(["^ *([1-3]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) "
                           "+([-0-9.Ee+]+)$",
                           lambda results, line: results.born_ion >= 0,
                           born_data])

            def born_section_stop(results, match):
                results.born_index = None

            search.append(["-------------------------------------",
                           lambda results, line: results.born_ion >= 1,
                           born_section_stop])

            self.born_ion = None
            self.born = {}

            micro_pyawk(self.filename, search, self)

        except:
            raise Exception("LEPSILON OUTCAR could not be parsed.")

    def read_lcalcpol(self):
        # variables to be filled
        self.p_elec = None
        self.p_ion = None
        try:
            search = []

            # Always present spin/non-spin
            def p_elc(results, match):
                results.p_elc = np.array([float(match.group(1)),
                                          float(match.group(2)),
                                          float(match.group(3))])

            search.append(["^.*Total electronic dipole moment: "
                           "*p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           None, p_elc])

            def p_ion(results, match):
                results.p_ion = np.array([float(match.group(1)),
                                          float(match.group(2)),
                                          float(match.group(3))])
            search.append(["^.*Ionic dipole moment: *p\[ion\]="
                           "\( *([-0-9.Ee+]*)"
                           " *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                           None, p_ion])

            micro_pyawk(self.filename, search, self)

        except:
            raise Exception("CLACLCPOL OUTCAR could not be parsed.")

    @property
    def to_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__, "efermi": self.efermi,
             "run_stats": self.run_stats, "magnetization": self.magnetization,
             "charge": self.charge, "total_magnetization": self.total_mag,
             "nelect": self.nelect, "is_stopped": self.is_stopped}
        return d


class VolumetricData(object):
    """
    Simple volumetric object for reading LOCPOT and CHGCAR type files.

    .. attribute:: structure

        Structure associated with the Volumetric Data object

    ..attribute:: is_spin_polarized

        True if run is spin polarized

    ..attribute:: dim

        Tuple of dimensions of volumetric grid in each direction (nx, ny, nz).

    ..attribute:: data

        Actual data as a dict of {string: np.array}. The string are "total"
        and "diff", in accordance to the output format of vasp LOCPOT and
        CHGCAR files where the total spin density is written first, followed
        by the difference spin density.

    .. attribute:: ngridpts

        Total number of grid points in volumetric data.
    """
    def __init__(self, structure, data, distance_matrix=None):
        """
        Typically, this constructor is not used directly and the static
        from_file constructor is used. This constructor is designed to allow
        summation and other operations between VolumetricData objects.

        Args:
            structure:
                Structure associated with the volumetric data
            data:
                Actual volumetric data.
            distance_matrix:
                A pre-computed distance matrix if available. Useful so pass
                distance_matrices between sums, shortcircuiting an otherwise
                expensive operation.
        """
        self.structure = structure
        self.is_spin_polarized = len(data) == 2
        self.dim = data["total"].shape
        self.data = data
        self.ngridpts = self.dim[0] * self.dim[1] * self.dim[2]
        #lazy init the spin data since this is not always needed.
        self._spin_data = {}
        self._distance_matrix = {} if not distance_matrix else distance_matrix

    @property
    def spin_data(self):
        """
        The data decomposed into actual spin data as {spin: data}.
        Essentially, this provides the actual Spin.up and Spin.down data
        instead of the total and diff.  Note that by definition, a
        non-spin-polarized run would have Spin.up data == Spin.down data.
        """
        if not self._spin_data:
            spin_data = dict()
            spin_data[Spin.up] = 0.5 * (self.data["total"] +
                                        self.data.get("diff", 0))
            spin_data[Spin.down] = 0.5 * (self.data["total"] -
                                          self.data.get("diff", 0))
            self._spin_data = spin_data
        return self._spin_data

    def get_axis_grid(self, ind):
        """
        Returns the grid for a particular axis.

        Args:
            ind:
                Axis index.
        """
        ng = self.dim
        num_pts = ng[ind]
        lengths = self.structure.lattice.abc
        return [i / num_pts * lengths[ind] for i in xrange(num_pts)]

    def __add__(self, other):
        return self.linear_add(other, 1.0)

    def __sub__(self, other):
        return self.linear_add(other, -1.0)

    def linear_add(self, other, scale_factor=1.0):
        """
        Method to do a linear sum of volumetric objects. Used by + and -
        operators as well. Returns a VolumetricData object containing the
        linear sum.

        Args:
            other:
                Another VolumetricData object
            scale_factor:
                Factor to scale the other data by.

        Returns:
            VolumetricData corresponding to self + scale_factor * other.
        """
        if self.structure != other.structure:
            raise ValueError("Adding or subtraction operations can only be "
                             "performed for volumetric data with the exact "
                             "same structure.")
        #To add checks
        data = {}
        for k in self.data.keys():
            data[k] = self.data[k] + scale_factor * other.data[k]
        return VolumetricData(self.structure, data, self._distance_matrix)

    @staticmethod
    def parse_file(filename):
        """
        Convenience method to parse a generic volumetric data file in the vasp
        like format. Used by subclasses for parsing file.

        Args:
            filename:
                Path of file to parse

        Returns:
            (poscar, data)
        """
        poscar_read = False
        poscar_string = []
        dataset = []
        all_dataset = []
        dim = None
        dimline = None
        read_dataset = False
        ngrid_pts = 0
        data_count = 0
        poscar = None
        with zopen(filename) as f:
            for line in f:
                line = line.strip()
                if read_dataset:
                    toks = line.split()
                    for tok in toks:
                        if data_count < ngrid_pts:
                            #This complicated procedure is necessary because
                            #vasp outputs x as the fastest index, followed by y
                            #then z.
                            x = data_count % dim[0]
                            y = int(math.floor(data_count / dim[0])) % dim[1]
                            z = int(math.floor(data_count / dim[0] / dim[1]))
                            dataset[x, y, z] = float(tok)
                            data_count += 1
                    if data_count >= ngrid_pts:
                        read_dataset = False
                        data_count = 0
                        all_dataset.append(dataset)
                elif not poscar_read:
                    if line != "" or len(poscar_string) == 0:
                        poscar_string.append(line)
                    elif line == "":
                        poscar = Poscar.from_string("\n".join(poscar_string))
                        poscar_read = True
                elif not dim:
                    dim = map(int, line.split())
                    ngrid_pts = dim[0] * dim[1] * dim[2]
                    dimline = line
                    read_dataset = True
                    dataset = np.zeros(dim)
                elif line == dimline:
                    read_dataset = True
                    dataset = np.zeros(dim)
            if len(all_dataset) == 2:
                data = {"total": all_dataset[0], "diff": all_dataset[1]}
            else:
                data = {"total": all_dataset[0]}
            return poscar, data

    def write_file(self, file_name, vasp4_compatible=False):
        """
        Write the VolumetricData object to a vasp compatible file.

        Args:
            file_name:
                the path to a file
            vasp4_compatible:
                True if the format is vasp4 compatible
        """

        f = zopen(file_name, "w")
        p = Poscar(self.structure)
        f.write(p.get_string(vasp4_compatible=vasp4_compatible) + "\n")
        a = self.dim
        f.write("\n")

        def write_spin(data_type):
            lines = []
            count = 0
            f.write("{} {} {}\n".format(a[0], a[1], a[2]))
            for (k, j, i) in itertools.product(xrange(a[2]), xrange(a[1]),
                                               xrange(a[0])):
                lines.append("%0.11e" % self.data[data_type][i, j, k])
                count += 1
                if count % 5 == 0:
                    f.write("".join(lines) + "\n")
                    lines = []
                else:
                    lines.append(" ")
            f.write("".join(lines) + "\n")

        write_spin("total")
        if self.is_spin_polarized:
            f.write("\n")
            write_spin("diff")
        f.close()

    def get_integrated_diff(self, ind, radius, nbins=1):
        """
        Get integrated difference of atom index ind up to radius. This can be
        an extremely computationally intensive process, depending on how many
        grid points are in the VolumetricData.

        Args:
            ind:
                Index of atom.
            radius:
                Radius of integration.
            nbins:
                Number of bins. Defaults to 1. This allows one to obtain the
                charge integration up to a list of the cumulative charge
                integration values for radii for [radius/nbins,
                2 * radius/nbins, ....].

        Returns:
            Differential integrated charge as a np array of [[radius, value],
            ...]. Format is for ease of plotting. E.g., plt.plot(data[:,0],
            data[:,1])
        """
        #For non-spin-polarized runs, this is zero by definition.
        if not self.is_spin_polarized:
            radii = [radius / nbins * (i + 1) for i in xrange(nbins)]
            data = np.zeros((nbins, 2))
            data[:, 0] = radii
            return data

        struct = self.structure
        a = self.dim
        if ind not in self._distance_matrix or\
                self._distance_matrix[ind]["max_radius"] < radius:
            coords = []
            for (x, y, z) in itertools.product(*[xrange(i) for i in a]):
                coords.append([x / a[0], y / a[1], z / a[2]])
            sites_dist = get_points_in_sphere_pbc(struct.lattice, coords,
                                                  struct[ind].coords,
                                                  radius)
            self._distance_matrix[ind] = {"max_radius": radius,
                                          "data": np.array(sites_dist)}

        data = self._distance_matrix[ind]["data"]

        #Use boolean indexing to find all charges within the desired distance.
        inds = data[:, 1] <= radius
        dists = data[inds, 1]
        data_inds = np.rint(np.mod(list(data[inds, 0]), 1) *
                            np.tile(a, (len(dists), 1)))
        vals = [self.data["diff"][x, y, z] for x, y, z in data_inds]

        hist, edges = np.histogram(dists, bins=nbins,
                                   range=[0, radius],
                                   weights=vals)
        data = np.zeros((nbins, 2))
        data[:, 0] = edges[1:]
        data[:, 1] = [sum(hist[0:i + 1]) / self.ngridpts
                      for i in xrange(nbins)]
        return data

    def get_average_along_axis(self, ind):
        """
        Get the averaged total of the volumetric data a certain axis direction.
        For example, useful for visualizing Hartree Potentials from a LOCPOT
        file.

        Args:
            ind : Index of axis.

        Returns:
            Average total along axis
        """
        m = self.data["total"]
        ng = self.dim
        if ind == 0:
            total = np.sum(np.sum(m, axis=1), 1)
        elif ind == 1:
            total = np.sum(np.sum(m, axis=0), 1)
        else:
            total = np.sum(np.sum(m, axis=0), 0)
        return total / ng[(ind + 1) % 3] / ng[(ind + 2) % 3]


class Locpot(VolumetricData):
    """
    Simple object for reading a LOCPOT file.
    """

    def __init__(self, poscar, data):
        """
        Args:
            poscar:
                Poscar object containing structure.
            data:
                Actual data.
        """
        VolumetricData.__init__(self, poscar.structure, data)
        self.name = poscar.comment

    @staticmethod
    def from_file(filename):
        (poscar, data) = VolumetricData.parse_file(filename)
        return Locpot(poscar, data)


class Chgcar(VolumetricData):
    """
    Simple object for reading a CHGCAR file.
    """

    def __init__(self, poscar, data):
        """
        Args:
            poscar:
                Poscar object containing structure.
            data:
                Actual data.
        """
        VolumetricData.__init__(self, poscar.structure, data)
        self.poscar = poscar
        self.name = poscar.comment
        self._distance_matrix = {}

    @staticmethod
    def from_file(filename):
        (poscar, data) = VolumetricData.parse_file(filename)
        return Chgcar(poscar, data)


class Procar(object):
    """
    Object for reading a PROCAR file.

    .. attribute:: data

        A nested dict containing the PROCAR data of the form below. It should
        be noted that VASP uses 1-based indexing for atoms, but this is
        converted to zero-based indexing in this parser to be consistent with
        representation of structures in pymatgen::

            {
                atom_index: {
                    kpoint_index: {
                        "bands": {
                            band_index: {
                                "p": 0.002,
                                "s": 0.025,
                                "d": 0.0
                            },
                            ...
                        },
                        "weight": 0.03125
                    },
                    ...
            }
    """
    def __init__(self, filename):
        """
        Args:
            filename:
                Name of file containing PROCAR.
        """
        data = defaultdict(dict)
        headers = None
        with zopen(filename, "r") as f:
            lines = list(clean_lines(f.readlines()))
            self.name = lines[0]
            kpointexpr = re.compile("^\s*k-point\s+(\d+).*weight = ([0-9\.]+)")
            bandexpr = re.compile("^\s*band\s+(\d+)")
            ionexpr = re.compile("^ion.*")
            expr = re.compile("^\s*([0-9]+)\s+")
            dataexpr = re.compile("[\.0-9]+")
            weight = 0
            current_kpoint = 0
            current_band = 0
            for l in lines:
                if bandexpr.match(l):
                    m = bandexpr.match(l)
                    current_band = int(m.group(1))
                elif kpointexpr.match(l):
                    m = kpointexpr.match(l)
                    current_kpoint = int(m.group(1))
                    weight = float(m.group(2))
                elif headers is None and ionexpr.match(l):
                    headers = l.split()
                    headers.pop(0)
                    headers.pop(-1)
                elif expr.match(l):
                    linedata = dataexpr.findall(l)
                    num_data = map(float, linedata)
                    #Convert to zero-based indexing for atoms.
                    index = int(num_data.pop(0)) - 1
                    num_data.pop(-1)
                    if current_kpoint not in data[index]:
                        data[index][current_kpoint] = {"weight": weight,
                                                       "bands": {}}
                    data[index][current_kpoint]["bands"][current_band] = \
                        dict(zip(headers, num_data))
            self.data = data

    def get_d_occupation(self, atom_index):
        """
        .. deprecated:: v2.6.4

            Use get_occpuation instead.

        Returns the d occupation of a particular atom.

        Args:
            atom_index:
                Index of atom in PROCAR. It should be noted that VASP uses
                1-based indexing for atoms, but this is converted to
                zero-based indexing in this parser to be consistent with
                representation of structures in pymatgen.

        Returns:
            d-occupation of atom at atom_index.
        """
        warnings.warn("get_d_occupation has been deprecated. Use "
                      "get_occupation instead.", DeprecationWarning)
        return self.get_occupation(atom_index, 'd')

    def get_occupation(self, atom_index, orbital):
        """
        Returns the occupation for a particular orbital of a particular atom.

        Args:
            atom_num:
                Index of atom in the PROCAR. It should be noted that VASP uses
                1-based indexing for atoms, but this is converted to
                zero-based indexing in this parser to be consistent with
                representation of structures in pymatgen.
            orbital:
                A string representing an orbital. If it is a single
                character, e.g., s, p, d or f, the sum of all s-type,
                p-type, d-type or f-type orbitals occupations are returned
                respectively. If it is a specific orbital, e.g., px, dxy,
                etc., only the occupation of that orbital is returned.

        Returns:
            Sum occupation of orbital of atom.
        """
        total = 0
        found = False
        for kpoint, d in self.data[atom_index].items():
            wt = d["weight"]
            for band, dd in d["bands"].items():
                for orb, v in dd.items():
                    if orb.startswith(orbital):
                        found = True
                        total += v * wt
        if not found:
            raise ValueError("Invalid orbital {}".format(orbital))
        return total


class Oszicar(object):
    """
    A basic parser for an OSZICAR output from VASP.  In general, while the
    OSZICAR is useful for a quick look at the output from a VASP run, we
    recommend that you use the Vasprun parser instead, which gives far richer
    information about a run.

    .. attribute:: electronic_steps

            All electronic steps as a list of list of dict. e.g.,
            [[{"rms": 160.0, "E": 4507.24605593, "dE": 4507.2, "N": 1,
            "deps": -17777.0, "ncg": 16576}, ...], [....]
            where electronic_steps[index] refers the list of electronic steps
            in one ionic_step, electronic_steps[index][subindex] refers to a
            particular electronic step at subindex in ionic step at index. The
            dict of properties depends on the type of VASP run, but in general,
            "E", "dE" and "rms" should be present in almost all runs.

    .. attribute:: ionic_steps:

            All ionic_steps as a list of dict, e.g.,
            [{"dE": -526.36, "E0": -526.36024, "mag": 0.0, "F": -526.36024},
            ...]
            This is the typical output from VASP at the end of each ionic step.
    """

    def __init__(self, filename):
        """
        Args:
            filename:
                Filename of file to parse
        """
        electronic_steps = []
        ionic_steps = []
        ionic_pattern = re.compile("(\d+)\s+F=\s*([\d\-\.E\+]+)\s+"
                                   "E0=\s*([\d\-\.E\+]+)\s+"
                                   "d\s*E\s*=\s*([\d\-\.E\+]+)\s+"
                                   "mag=\s*([\d\-\.E\+]+)")
        electronic_pattern = re.compile("\s*\w+\s*:(.*)")

        def smart_convert(header, num):
            try:
                if header == "N" or header == "ncg":
                    v = int(num)
                    return v
                v = float(num)
                return v
            except ValueError:
                return "--"

        header = []
        with zopen(filename, "r") as fid:
            for line in fid:
                line = line.strip()
                m = electronic_pattern.match(line)
                if m:
                    toks = m.group(1).split()
                    data = {header[i]: smart_convert(header[i], toks[i])
                            for i in xrange(len(toks))}
                    if toks[0] == "1":
                        electronic_steps.append([data])
                    else:
                        electronic_steps[-1].append(data)
                elif ionic_pattern.match(line.strip()):
                    m = ionic_pattern.match(line.strip())
                    ionic_steps.append({"F": float(m.group(2)),
                                        "E0": float(m.group(3)),
                                        "dE": float(m.group(4)),
                                        "mag": float(m.group(5))})
                elif re.match("^\s*N\s+E\s*", line):
                    header = line.strip().replace("d eps", "deps").split()
        self.electronic_steps = electronic_steps
        self.ionic_steps = ionic_steps

    @property
    def all_energies(self):
        """
        Compilation of all energies from all electronic steps and ionic steps
        as a tuple of list of energies, e.g.,
        ((4507.24605593, 143.824705755, -512.073149912, ...), ...)
        """
        all_energies = []
        for i in xrange(len(self.electronic_steps)):
            energies = [step["E"] for step in self.electronic_steps[i]]
            energies.append(self.ionic_steps[i]["F"])
            all_energies.append(tuple(energies))
        return tuple(all_energies)

    @property
    @unitized("eV")
    def final_energy(self):
        """
        Final energy from run.
        """
        return self.ionic_steps[-1]["F"]

    @property
    def to_dict(self):
        return {"electronic_steps": self.electronic_steps,
                "ionic_steps": self.ionic_steps}


class VaspParserError(Exception):
    """
    Exception class for VASP parsing.
    """
    pass


def get_band_structure_from_vasp_multiple_branches(dir_name, efermi=None,
                                                   projections=False):
    """
    this method is used to get band structure info from a VASP directory. It
    takes into account that the run can be divided in several branches named
    "branch_x". If the run has not been divided in branches the method will
    turn to parsing vasprun.xml directly.

    The method returns None is there"s a parsing error

    Args:
        dir_name:
            Directory containing all bandstructure runs.
        efermi:
            Efermi for bandstructure.
        projections:
            True if you want to get the data on site projections if any
            Note that this is sometimes very large

    Returns:
        A BandStructure Object
    """
    #ToDo: Add better error handling!!!
    if os.path.exists(os.path.join(dir_name, "branch_0")):
        #get all branch dir names
        branch_dir_names = [os.path.abspath(d)
                            for d in glob.glob("{i}/branch_*"
                                               .format(i=dir_name))
                            if os.path.isdir(d)]

        #sort by the directory name (e.g, branch_10)
        sort_by = lambda x: int(x.split("_")[-1])
        sorted_branch_dir_names = sorted(branch_dir_names, key=sort_by)

        # populate branches with Bandstructure instances
        branches = []
        for dir_name in sorted_branch_dir_names:
            xml_file = os.path.join(dir_name, "vasprun.xml")
            if os.path.exists(xml_file):
                run = Vasprun(xml_file, parse_projected_eigen=projections)
                branches.append(run.get_band_structure(efermi=efermi))
            else:
                # It might be better to throw an exception
                warnings.warn("Skipping {}. Unable to find {}"
                              .format(d=dir_name, f=xml_file))

        return get_reconstructed_band_structure(branches, efermi)
    else:
        xml_file = os.path.join(dir_name, "vasprun.xml")
        #Better handling of Errors
        if os.path.exists(xml_file):
            return Vasprun(xml_file, parse_projected_eigen=projections)\
                .get_band_structure(kpoints_filename=None, efermi=efermi)
        else:
            return None


def get_adjusted_fermi_level(efermi, cbm, band_structure):
    """
    When running a band structure computations the fermi level needs to be
    take from the static run that gave the charge density used for the non-self
    consistent band structure run. Sometimes this fermi level is however a
    little too low because of the mismatch between the uniform grid used in
    the static run and the band structure k-points (e.g., the VBM is on Gamma
    and the Gamma point is not in the uniform mesh). Here we use a procedure
    consisting in looking for energy levels higher than the static fermi level
    (but lower than the LUMO) if any of these levels make the band structure
    appears insulating and not metallic anymore, we keep this adjusted fermi
    level. This procedure has shown to detect correctly most insulators.

    Args:
        cbm:
            the cbm of the static run
        efermi:
            the fermi energy of the static run
        run_bandstructure:
            a band_structure object

    Returns:
        a new adjusted fermi level
    """
    #make a working copy of band_structure
    bs_working = BandStructureSymmLine.from_dict(band_structure.to_dict)
    if bs_working.is_metal():
        e = efermi
        while e < cbm:
            e += 0.01
            bs_working._efermi = e
            if not bs_working.is_metal():
                return e
    return efermi
