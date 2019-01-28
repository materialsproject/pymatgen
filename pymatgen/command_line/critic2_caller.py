# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import subprocess
import warnings
import numpy as np
import glob

from scipy.spatial import KDTree
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.analysis.graphs import StructureGraph
from monty.os.path import which
from monty.dev import requires
from monty.json import MSONable
from monty.tempfile import ScratchDir
from enum import Enum

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
This module implements an interface to the critic2 Bader analysis code.

For most Bader analysis purposes, users are referred to
pymatgen.command_line.bader_caller instead, this module is for advanced
usage requiring identification of critical points in the charge density.

This module depends on a compiled critic2 executable available in the path.
Please follow the instructions at https://github.com/aoterodelaroza/critic2
to compile or, if using macOS and homebrew, use `brew tap homebrew/science`
followed by `brew install critic2`.

New users are *strongly* recommended to read the critic2 manual first.

In brief,
* critic2 searches for critical points in charge density
* a critical point can be one of four types: nucleus, bond, ring
or cage
* it does this by seeding locations for likely critical points
and then searching in these regions
* there are two lists of critical points in the output, a list
of non-equivalent points (with in-depth information about the
field at those points), and a full list of points generated
by the appropriate symmetry operations
* connectivity between these points is also provided when
appropriate (e.g. the two nucleus critical points linked to
 a bond critical point)
* critic2 can do many other things besides

If you use this module, please cite the following:

A. Otero-de-la-Roza, E. R. Johnson and V. Luaña,
Comput. Phys. Commun. 185, 1007-1018 (2014)
(http://dx.doi.org/10.1016/j.cpc.2013.10.026)

A. Otero-de-la-Roza, M. A. Blanco, A. Martín Pendás and
V. Luaña, Comput. Phys. Commun. 180, 157–166 (2009)
(http://dx.doi.org/10.1016/j.cpc.2008.07.018)
"""

__author__ = "Matthew Horton"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Production"
__date__ = "July 2017"


class Critic2Caller:

    @requires(which("critic2"), "Critic2Caller requires the executable critic to be in the path. "
                                "Please follow the instructions at https://github.com/aoterodelaroza/critic2.")
    def __init__(self, structure, chgcar=None, chgcar_ref=None,
                 user_input_settings=None, write_cml=False):
        """
        Run Critic2 in automatic mode on a supplied structure, charge
        density (chgcar) and reference charge density (chgcar_ref).

        The reason for a separate reference field is that in
        VASP, the CHGCAR charge density only contains valence
        electrons and may be missing substantial charge at
        nuclei leading to misleading results. Thus, a reference
        field is commonly constructed from the sum of AECCAR0
        and AECCAR2 which is the total charge density, but then
        the valence charge density is used for the final analysis.

        If chgcar_ref is not supplied, chgcar will be used as the
        reference field. If chgcar is not supplied, the promolecular
        charge density will be used as the reference field -- this can
        often still give useful results if only topological information
        is wanted.

        User settings is a dictionary that can contain:
        * GRADEPS, float (field units), gradient norm threshold
        * CPEPS, float (Bohr units in crystals), minimum distance between
          critical points for them to be equivalent
        * NUCEPS, same as CPEPS but specifically for nucleus critical
          points (critic2 default is depedent on grid dimensions)
        * NUCEPSH, same as NUCEPS but specifically for hydrogen nuclei
          since associated charge density can be significantly displaced
          from hydrogen nucleus
        * EPSDEGEN, float (field units), discard critical point if any
          element of the diagonal of the Hessian is below this value,
          useful for discarding points in vacuum regions
        * DISCARD, float (field units), discard critical points with field
          value below this value, useful for discarding points in vacuum
          regions
        * SEED, list of strings, strategies for seeding points, default
          is ['WS 1', 'PAIR 10'] which seeds critical points by
          sub-dividing the Wigner-Seitz cell and between every atom pair
          closer than 10 Bohr, see critic2 manual for more options

        :param structure: Structure to analyze
        :param chgcar: Charge density to use for analysis. If None, will
        use promolecular density.
        :param chgcar_ref: Reference charge density. If None, will use
        chgcar as reference.
        :param user_input_settings (dict):
        :param write_cml (bool): Useful for debug, if True will write all
        critical points to a file 'table.cml' in the working directory
        useful for visualization
        """

        settings = {'CPEPS': 0.1, 'SEED': ["WS", "PAIR DIST 10"]}
        if user_input_settings:
            settings.update(user_input_settings)

        # Load crystal structure
        input_script = ["crystal POSCAR"]

        # Load data to use as reference field
        if chgcar_ref:
            input_script += ["load VASPCHG CHGCAR_ref id chgcar_ref",
                             "reference chgcar_ref"]

        # Load data to use for analysis
        if chgcar:
            input_script += ["load VASPCHG CHGCAR_int id chgcar",
                             "integrable chgcar"]

        # Command to run automatic analysis
        auto = "auto "
        for k, v in settings.items():
            if isinstance(v, list):
                for item in v:
                    auto += '{} {} '.format(k, item)
            else:
                auto += '{} {} '.format(k, v)
        input_script += [auto]

        if write_cml:
            input_script += ["cpreport ../table.cml cell border graph"]

        input_script = "\n".join(input_script)

        with ScratchDir(".") as temp_dir:

            os.chdir(temp_dir)

            with open('input_script.cri', 'w') as f:
                f.write(input_script)

            structure.to(filename='POSCAR')

            if chgcar:
                chgcar.write_file('CHGCAR_int')

            if chgcar_ref:
                chgcar_ref.write_file('CHGCAR_ref')

            args = ["critic2", "input_script.cri"]
            rs = subprocess.Popen(args,
                                  stdout=subprocess.PIPE,
                                  stdin=subprocess.PIPE, close_fds=True)

            stdout, stderr = rs.communicate()
            stdout = stdout.decode()

            if stderr:
                stderr = stderr.decode()
                warnings.warn(stderr)

            if rs.returncode != 0:
                raise RuntimeError("critic2 exited with return code {}.".format(rs.returncode))

            self._stdout = stdout
            self._stderr = stderr

            self.output = Critic2Output(structure, stdout)

    @classmethod
    def from_path(cls, path, suffix=''):
        """
        Convenience method to run critic2 analysis on a folder containing
        typical VASP output files.
        This method will:
        1. Look for files CHGCAR, AECAR0, AECAR2, POTCAR or their gzipped
        counterparts.
        2. If AECCAR* files are present, constructs a temporary reference
        file as AECCAR0 + AECCAR2
        3. Runs critic2 analysis twice: once for charge, and a second time
        for the charge difference (magnetization density).
        :param path: path to folder to search in
        :param suffix: specific suffix to look for (e.g. '.relax1' for
        'CHGCAR.relax1.gz')
        :return:
        """

        def _get_filepath(filename, warning, path=path, suffix=suffix):
            paths = glob.glob(os.path.join(path, filename + suffix + '*'))
            if not paths:
                warnings.warn(warning)
                return None
            if len(paths) > 1:
                # using reverse=True because, if multiple files are present,
                # they likely have suffixes 'static', 'relax', 'relax2', etc.
                # and this would give 'static' over 'relax2' over 'relax'
                # however, better to use 'suffix' kwarg to avoid this!
                paths.sort(reverse=True)
                warnings.warn('Multiple files detected, using {}'.format(os.path.basename(path)))
            path = paths[0]
            return path

        chgcar_path = _get_filepath('CHGCAR', 'Could not find CHGCAR!')
        chgcar = Chgcar.from_file(chgcar_path)

        aeccar0_path = _get_filepath('AECCAR0', 'Could not find AECCAR0, interpret Bader results with caution.')
        aeccar0 = Chgcar.from_file(aeccar0_path) if aeccar0_path else None

        aeccar2_path = _get_filepath('AECCAR2', 'Could not find AECCAR2, interpret Bader results with caution.')
        aeccar2 = Chgcar.from_file(aeccar2_path) if aeccar2_path else None

        chgcar_ref = aeccar0.linear_add(aeccar2) if (aeccar0 and aeccar2) else None

        return cls(chgcar.structure, chgcar, chgcar_ref)


class CriticalPointType(Enum):

    nucleus = "nucleus"  # (3, -3)
    bond = "bond"  # (3, -1)
    ring = "ring"  # (3, 1)
    cage = "cage"  # (3, 3)


class CriticalPoint(MSONable):

    def __init__(self, index, type, frac_coords, point_group,
                 multiplicity, field, field_gradient,
                 coords=None, field_hessian=None):
        """
        Class to characterise a critical point from a topological
        analysis of electron charge density.

        Note this class is usually associated with a Structure, so
        has information on multiplicity/point group symmetry.

        :param index: index of point
        :param type: type of point, given as a string
        :param coords: Cartesian co-ordinates in Angstroms
        :param frac_coords: fractional co-ordinates
        :param point_group: point group associated with critical point
        :param multiplicity: number of equivalent critical points
        :param field: value of field at point (f)
        :param field_gradient: gradient of field at point (grad f)
        :param field_hessian: hessian of field at point (del^2 f)
        """
        self.index = index
        self._type = type
        self.coords = coords
        self.frac_coords = frac_coords
        self.point_group = point_group
        self.multiplicity = multiplicity
        self.field = field
        self.field_gradient = field_gradient
        self.field_hessian = field_hessian

    @property
    def type(self):
        return CriticalPointType(self._type)

    def __str__(self):
        return "Critical Point: {} ({})".format(self.type.name, self.frac_coords)

    @property
    def laplacian(self):
        return np.trace(self.field_hessian)

    @property
    def ellipticity(self):
        '''
        Most meaningful for bond critical points,
        can be physically interpreted as e.g. degree
        of pi-bonding in organic molecules. Consult
        literature for more information.
        :return:
        '''
        eig = np.linalg.eig(self.field_hessian)
        eig.sort()
        return eig[0]/eig[1] - 1


class Critic2Output(MSONable):

    def __init__(self, structure, critic2_stdout):
        """
        This class is used to store results from the Critic2Caller.

        To explore the bond graph, use the "structure_graph"
        method, which returns a user-friendly StructureGraph
        class with bonding information. By default, this returns
        a StructureGraph with edge weights as bond lengths, but
        can optionally return a graph with edge weights as any
        property supported by the `CriticalPoint` class, such as
        bond ellipticity.

        This class also provides an interface to explore just the
        non-symmetrically-equivalent critical points via the
        `critical_points` attribute, and also all critical
        points (via nodes dict) and connections between them
        (via edges dict). The user should be familiar with critic2
        before trying to understand these.

        Indexes of nucleus critical points in the nodes dict are the
        same as the corresponding sites in structure, with indices of
        other critical points arbitrarily assigned.

        :param structure: associated Structure
        :param critic2_stdout: stdout from running critic2 in automatic
        mode
        """

        self.structure = structure

        self._critic2_stdout = critic2_stdout

        self.nodes = {}
        self.edges = {}

        self._parse_stdout(critic2_stdout)

    def structure_graph(self, edge_weight="bond_length", edge_weight_units="Å"):
        """
        A StructureGraph object describing bonding information
        in the crystal. Lazily constructed.

        :param edge_weight: a value to store on the Graph edges,
        by default this is "bond_length" but other supported
        values are any of the attributes of CriticalPoint
        :return:
        """

        sg = StructureGraph.with_empty_graph(self.structure, name="bonds",
                                             edge_weight_name=edge_weight,
                                             edge_weight_units=edge_weight_units)

        edges = self.edges.copy()
        idx_to_delete = []
        # check for duplicate bonds
        for idx, edge in edges.items():
            unique_idx = self.nodes[idx]['unique_idx']
            # only check edges representing bonds, not rings
            if self.critical_points[unique_idx].type == CriticalPointType.bond:
                if idx not in idx_to_delete:
                    for idx2, edge2 in edges.items():
                        if idx != idx2 and edge == edge2:
                            idx_to_delete.append(idx2)
                            warnings.warn("Duplicate edge detected, try re-running "
                                          "critic2 with custom parameters to fix this. "
                                          "Mostly harmless unless user is also "
                                          "interested in rings/cages.")
                            logger.debug("Duplicate edge between points {} (unique point {})"
                                         "and {} ({}).".format(idx, self.nodes[idx]['unique_idx'],
                                                               idx2, self.nodes[idx2]['unique_idx']))
        # and remove any duplicate bonds present
        for idx in idx_to_delete:
            del edges[idx]

        for idx, edge in edges.items():
            unique_idx = self.nodes[idx]['unique_idx']
            # only add edges representing bonds, not rings
            if self.critical_points[unique_idx].type == CriticalPointType.bond:

                from_idx = edge['from_idx']
                to_idx = edge['to_idx']

                from_lvec = edge['from_lvec']
                to_lvec = edge['to_lvec']

                relative_lvec = np.subtract(to_lvec, from_lvec)

                if edge_weight == "bond_length":
                    weight = self.structure.get_distance(from_idx, to_idx, jimage=relative_lvec)
                else:
                    weight = getattr(self.critical_points[unique_idx],
                                     edge_weight, None)

                sg.add_edge(from_idx, to_idx,
                            from_jimage=from_lvec, to_jimage=to_lvec,
                            weight=weight)

        return sg

    def get_critical_point_for_site(self, n):
        return self.critical_points[self.nodes[n]['unique_idx']]

    def _parse_stdout(self, stdout):

        stdout = stdout.split("\n")

        # NOTE WE ARE USING 0-BASED INDEXING:
        # This is different from critic2 which
        # uses 1-based indexing, so all parsed
        # indices have 1 subtracted.

        # Parsing happens in two stages:

        # 1. We construct a list of unique critical points
        #    (i.e. non-equivalent by the symmetry of the crystal)
        #   and the properties of the field at those points

        # 2. We construct a list of nodes and edges describing
        #    all critical points in the crystal

        # Steps 1. and 2. are essentially indepdendent, except
        # that the critical points in 2. have a pointer to their
        # associated unique critical point in 1. so that more
        # information on that point can be retrieved if necessary.

        unique_critical_points = []

        # parse unique critical points
        for i, line in enumerate(stdout):
            if "* Critical point list, final report (non-equivalent cps)" in line:
                start_i = i + 4
            elif "* Analysis of system bonds" in line:
                end_i = i - 2
        # if start_i and end_i haven't been found, we
        # need to re-evaluate assumptions in this parser!

        for i, line in enumerate(stdout):
            if i >= start_i and i <= end_i:
                l = line.replace("(", "").replace(")", "").split()

                unique_idx = int(l[0]) - 1
                point_group = l[1]
                # type = l[2]  # type from definition of critical point e.g. (3, -3)
                type = l[3]  # type from name, e.g. nucleus
                frac_coords = [float(l[4]), float(l[5]), float(l[6])]
                multiplicity = float(l[7])
                # name = float(l[8])
                field = float(l[9])
                field_gradient = float(l[10])
                # laplacian = float(l[11])

                point = CriticalPoint(unique_idx, type, frac_coords, point_group,
                                      multiplicity, field, field_gradient)
                unique_critical_points.append(point)

        # TODO: may be other useful information to parse here too
        for i, line in enumerate(stdout):
            if '+ Critical point no.' in line:
                unique_idx = int(line.split()[4]) - 1
            elif "Hessian:" in line:
                l1 = list(map(float, stdout[i + 1].split()))
                l2 = list(map(float, stdout[i + 2].split()))
                l3 = list(map(float, stdout[i + 3].split()))
                hessian = [[l1[0], l1[1], l1[2]],
                           [l2[0], l2[1], l2[2]],
                           [l3[0], l3[1], l3[2]]]
                unique_critical_points[unique_idx].field_hessian = hessian

        self.critical_points = unique_critical_points

        # parse graph connecting critical points
        for i, line in enumerate(stdout):
            if "* Complete CP list, bcp and rcp connectivity table" in line:
                start_i = i + 3
            elif "* Attractor connectivity matrix" in line:
                end_i = i - 2
        # if start_i and end_i haven't been found, we
        # need to re-evaluate assumptions in this parser!

        # Order of nuclei provided by critic2 doesn't
        # necessarily match order of sites in Structure.
        # We perform a mapping from one to the other,
        # and re-index all nodes accordingly.
        node_mapping = {}  # critic2_index:structure_index
        # ensure frac coords are in [0,1] range
        frac_coords = np.array(self.structure.frac_coords) % 1
        kd = KDTree(frac_coords)
        for i, line in enumerate(stdout):
            if i >= start_i and i <= end_i:
                l = line.split()
                if l[2] == "n":
                    critic2_idx = int(l[0]) - 1
                    frac_coord = np.array([float(l[3]), float(l[4]), float(l[5])]) % 1
                    node_mapping[critic2_idx] = kd.query(frac_coord)[1]

        if len(node_mapping) != len(self.structure):
            warnings.warn("Check that all sites in input structure have "
                          "been detected by critic2.")

        def _remap(critic2_idx):
            return node_mapping.get(critic2_idx, critic2_idx)

        for i, line in enumerate(stdout):
            if i >= start_i and i <= end_i:

                l = line.replace("(", "").replace(")", "").split()

                idx = _remap(int(l[0]) - 1)
                unique_idx = int(l[1]) - 1
                frac_coords = [float(l[3]), float(l[4]), float(l[5])]

                self._add_node(idx, unique_idx, frac_coords)
                if len(l) > 6:
                    from_idx = _remap(int(l[6]) - 1)
                    to_idx = _remap(int(l[10]) - 1)
                    self._add_edge(idx, from_idx=from_idx, from_lvec=(int(l[7]), int(l[8]), int(l[9])),
                                  to_idx=to_idx, to_lvec=(int(l[11]), int(l[12]), int(l[13])))

        self._map = node_mapping

    def _add_node(self, idx, unique_idx, frac_coords):
        """
        Add information about a node describing a critical point.

        :param idx: unique index
        :param unique_idx: index of unique CriticalPoint,
        used to look up more information of point (field etc.)
        :param frac_coord: fractional co-ordinates of point
        :return:
        """
        self.nodes[idx] = {'unique_idx': unique_idx, 'frac_coords': frac_coords}

    def _add_edge(self, idx, from_idx, from_lvec, to_idx, to_lvec):
        """
        Add information about an edge linking two critical points.

        This actually describes two edges:

        from_idx ------ idx ------ to_idx

        However, in practice, from_idx and to_idx will typically be
        atom nuclei, with the center node (idx) referring to a bond
        critical point. Thus, it will be more convenient to model
        this as a single edge linking nuclei with the properties
        of the bond critical point stored as an edge attribute.

        :param idx: index of node
        :param from_idx: from index of node
        :param from_lvec: vector of lattice image the from node is in
        as tuple of ints
        :param to_idx: to index of node
        :param to_lvec:  vector of lattice image the to node is in as
        tuple of ints
        :return:
        """
        self.edges[idx] = {'from_idx': from_idx, 'from_lvec': from_lvec,
                           'to_idx': to_idx, 'to_lvec': to_lvec}
