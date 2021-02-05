# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Provides a class for interacting with KPath classes to
generate high-symmetry k-paths using different conventions.
"""

import itertools
from warnings import warn

import networkx as nx
import numpy as np


from pymatgen.symmetry.kpath import (
    KPathBase,
    KPathLatimerMunro,
    KPathSeek,
    KPathSetyawanCurtarolo,
)

from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin

__author__ = "Jason Munro"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Jason Munro"
__email__ = "jmunro@lbl.gov"
__status__ = "Development"
__date__ = "March 2020"


class HighSymmKpath(KPathBase):
    """
    This class generates path along high symmetry lines in the
    Brillouin zone according to different conventions.
    The class is designed to be used with a specific primitive
    cell setting. The definitions for the primitive cell
    used can be found in: Computational Materials Science,
    49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010.
    The space group analyzer can be used to produce the correct
    primitive structure
    (method get_primitive_standard_structure(international_monoclinic=False)).
    Ensure input structure is correct before 'get_kpoints()' method is used.
    See individual KPath classes for details on specific conventions.
    """

    def __init__(
        self,
        structure,
        has_magmoms=False,
        magmom_axis=None,
        path_type="sc",
        symprec=0.01,
        angle_tolerance=5,
        atol=1e-5,
    ):
        """
        Args:
            structure (Structure): Structure object
            has_magmoms (boolean): Whether the input structure contains
                magnetic moments as site properties with the key 'magmom.'
                Values may be in the form of 3-component vectors given in
                the basis of the input lattice vectors, in
                which case the spin axis will default to a_3, the third
                real-space lattice vector (this triggers a warning).
            magmom_axis (list or numpy array): 3-component vector specifying
                direction along which magnetic moments given as scalars
                should point. If all magnetic moments are provided as
                vectors then this argument is not used.
            path_type (string): Chooses which convention to use to generate
                the high symmetry path. Options are: 'sc', 'hin', 'lm' for the
                Setyawan & Curtarolo, Hinuma et al., and  Latimer & Munro conventions.
                Choosing 'all' will generate one path with points from all three
                conventions. Equivalent labels between each will also be generated.
                Order will always be Latimer & Munro, Setyawan & Curtarolo, and Hinuma et al.
                Lengths for each of the paths will also be generated and output
                as a list. Note for 'all' the user will have to alter the labels on
                their own for plotting.
            symprec (float): Tolerance for symmetry finding
            angle_tolerance (float): Angle tolerance for symmetry finding.
            atol (float): Absolute tolerance used to determine symmetric
                equivalence of points and lines on the BZ.
        """

        super().__init__(structure, symprec=symprec, angle_tolerance=angle_tolerance, atol=atol)

        self._path_type = path_type

        self._equiv_labels = None
        self._path_lengths = None
        self._label_index = None

        if path_type != "all":

            if path_type == "lm":
                self._kpath = self._get_lm_kpath(has_magmoms, magmom_axis, symprec, angle_tolerance, atol).kpath
            elif path_type == "sc":
                self._kpath = self._get_sc_kpath(symprec, angle_tolerance, atol).kpath
            elif path_type == "hin":
                hin_dat = self._get_hin_kpath(symprec, angle_tolerance, atol, not has_magmoms)
                self._kpath = hin_dat.kpath
                self._hin_tmat = hin_dat._tmat

        else:

            if has_magmoms:
                raise ValueError("Cannot select 'all' with non-zero magmoms.")

            lm_bs = self._get_lm_kpath(has_magmoms, magmom_axis, symprec, angle_tolerance, atol)
            rpg = lm_bs._rpg

            sc_bs = self._get_sc_kpath(symprec, angle_tolerance, atol)
            hin_bs = self._get_hin_kpath(symprec, angle_tolerance, atol, not has_magmoms)

            index = 0
            cat_points = {}
            label_index = {}
            num_path = []
            self._path_lengths = []

            for bs in [lm_bs, sc_bs, hin_bs]:
                for key, value in enumerate(bs.kpath["kpoints"]):
                    cat_points[index] = bs.kpath["kpoints"][value]
                    label_index[index] = value
                    index += 1

                total_points_path = 0
                for seg in bs.kpath["path"]:
                    total_points_path += len(seg)

                for block in bs.kpath["path"]:
                    new_block = []
                    for label in block:
                        for ind in range(
                            len(label_index) - len(bs.kpath["kpoints"]),
                            len(label_index),
                        ):
                            if label_index[ind] == label:
                                new_block.append(ind)

                    num_path.append(new_block)

                self._path_lengths.append(total_points_path)

            self._label_index = label_index

            self._kpath = {"kpoints": cat_points, "path": num_path}

            self._equiv_labels = self._get_klabels(lm_bs, sc_bs, hin_bs, rpg)

    @property
    def path_type(self):
        """
        Returns:
        The type of kpath chosen
        """
        return self._path_type

    @property
    def label_index(self):
        """
        Returns:
        The correspondance between numbers and kpoint symbols for the
        combined kpath generated when path_type = 'all'. None otherwise.
        """
        return self._label_index

    @property
    def equiv_labels(self):
        """
        Returns:
        The correspondance between the kpoint symbols in the Latimer and
        Munro convention, Setyawan and Curtarolo, and Hinuma
        conventions respectively. Only generated when path_type = 'all'.
        """
        return self._equiv_labels

    @property
    def path_lengths(self):
        """
        Returns:
        List of lengths of the Latimer and Munro, Setyawan and Curtarolo, and Hinuma
        conventions in the combined HighSymmKpath object when path_type = 'all' respectively.
        None otherwise.
        """
        return self._path_lengths

    def _get_lm_kpath(self, has_magmoms, magmom_axis, symprec, angle_tolerance, atol):
        """
        Returns:
        Latimer and Munro k-path with labels.
        """

        return KPathLatimerMunro(self._structure, has_magmoms, magmom_axis, symprec, angle_tolerance, atol)

    def _get_sc_kpath(self, symprec, angle_tolerance, atol):
        """
        Returns:
        Setyawan and Curtarolo k-path with labels.
        """
        kpath = KPathSetyawanCurtarolo(self._structure, symprec, angle_tolerance, atol)

        self.prim = kpath.prim
        self.conventional = kpath.conventional
        self.prim_rec = kpath.prim_rec
        self._rec_lattice = self.prim_rec

        return kpath

    def _get_hin_kpath(self, symprec, angle_tolerance, atol, tri):
        """
        Returns:
        Hinuma et al. k-path with labels.
        """
        bs = KPathSeek(self._structure, symprec, angle_tolerance, atol, tri)

        kpoints = bs.kpath["kpoints"]
        tmat = bs._tmat
        for key in kpoints:
            kpoints[key] = np.dot(np.transpose(np.linalg.inv(tmat)), kpoints[key])

        bs.kpath["kpoints"] = kpoints
        self._rec_lattice = self._structure.lattice.reciprocal_lattice

        warn(
            "K-path from the Hinuma et al. convention has been transformed to the basis of the reciprocal lattice \
of the input structure. Use `KPathSeek` for the path in the original author-intended basis."
        )

        return bs

    def _get_klabels(self, lm_bs, sc_bs, hin_bs, rpg):
        """
        Returns:
        labels (dict): Dictionary of equivalent labels for paths if 'all' is chosen.
            If an exact kpoint match cannot be found, symmetric equivalency will be
            searched for and indicated with an asterisk in the equivalent label.
            If an equivalent label can still not be found, or the point is not in
            the explicit kpath, its equivalent label will be set to itself in the output.
        """

        lm_path = lm_bs.kpath
        sc_path = sc_bs.kpath
        hin_path = hin_bs.kpath

        n_op = len(rpg)

        pairs = itertools.permutations([{"sc": sc_path}, {"lm": lm_path}, {"hin": hin_path}], r=2)
        labels = {"sc": {}, "lm": {}, "hin": {}}

        for (a, b) in pairs:
            [(a_type, a_path)] = list(a.items())
            [(b_type, b_path)] = list(b.items())

            sc_count = np.zeros(n_op)

            for o_num in range(0, n_op):
                a_tr_coord = []

                for (label_a, coord_a) in a_path["kpoints"].items():
                    a_tr_coord.append(np.dot(rpg[o_num], coord_a))

                for coord_a in a_tr_coord:
                    for key, value in b_path["kpoints"].items():
                        if np.allclose(value, coord_a, atol=self._atol):
                            sc_count[o_num] += 1
                            break

            a_to_b_labels = {}
            unlabeled = {}

            for (label_a, coord_a) in a_path["kpoints"].items():
                coord_a_t = np.dot(rpg[np.argmax(sc_count)], coord_a)
                assigned = False

                for (label_b, coord_b) in b_path["kpoints"].items():
                    if np.allclose(coord_b, coord_a_t, atol=self._atol):
                        a_to_b_labels[label_a] = label_b
                        assigned = True
                        break

                if not assigned:
                    unlabeled[label_a] = coord_a

            for (label_a, coord_a) in unlabeled.items():
                for op in rpg:
                    coord_a_t = np.dot(op, coord_a)
                    key = [
                        key
                        for key, value in b_path["kpoints"].items()
                        if np.allclose(value, coord_a_t, atol=self._atol)
                    ]

                    if key != []:
                        a_to_b_labels[label_a] = key[0][0] + "^{*}"
                        break

                if key == []:
                    a_to_b_labels[label_a] = label_a

            labels[a_type][b_type] = a_to_b_labels

        return labels

    @staticmethod
    def get_continuous_path(bandstructure):
        """
        Obtain a continous version of an inputted path using graph theory.
        This routine will attempt to add connections between nodes of
        odd-degree to ensure a Eulerian path can be formed. Initial
        k-path must be able to be converted to a connected graph. See
        npj Comput Mater 6, 112 (2020). 10.1038/s41524-020-00383-7
        for more details.

        Args:
        bandstructure (BandstructureSymmLine): BandstructureSymmLine object.

        Returns:
        bandstructure (BandstructureSymmLine): New BandstructureSymmLine object with continous path.
        """

        G = nx.Graph()

        labels = []
        for point in bandstructure.kpoints:
            if point.label is not None:
                labels.append(point.label)

        plot_axis = []
        for i in range(int(len(labels) / 2)):
            G.add_edges_from([(labels[2 * i], labels[(2 * i) + 1])])
            plot_axis.append((labels[2 * i], labels[(2 * i) + 1]))

        G_euler = nx.algorithms.euler.eulerize(G)

        G_euler_circuit = nx.algorithms.euler.eulerian_circuit(G_euler)

        distances_map = []
        kpath_euler = []

        for edge_euler in G_euler_circuit:
            kpath_euler.append(edge_euler)
            for edge_reg in plot_axis:
                if edge_euler == edge_reg:
                    distances_map.append((plot_axis.index(edge_reg), False))
                elif edge_euler[::-1] == edge_reg:
                    distances_map.append((plot_axis.index(edge_reg), True))

        if bandstructure.is_spin_polarized:
            spins = [Spin.up, Spin.down]
        else:
            spins = [Spin.up]

        new_kpoints = []
        new_bands = {spin: [np.array([]) for _ in range(bandstructure.nb_bands)] for spin in spins}
        new_projections = {spin: [[] for _ in range(bandstructure.nb_bands)] for spin in spins}

        for entry in distances_map:
            if not entry[1]:
                branch = bandstructure.branches[entry[0]]
                start = branch["start_index"]
                stop = branch["end_index"] + 1
                step = 1

            else:
                branch = bandstructure.branches[entry[0]]
                start = branch["end_index"]
                stop = branch["start_index"] - 1
                step = -1

            # kpoints
            new_kpoints += [point.frac_coords for point in bandstructure.kpoints[start:stop:step]]

            # eigenvals
            for spin in spins:
                for n, band in enumerate(bandstructure.bands[spin]):

                    new_bands[spin][n] = np.concatenate((new_bands[spin][n], band[start:stop:step]))

            # projections
            for spin in spins:
                for n, band in enumerate(bandstructure.projections[spin]):

                    new_projections[spin][n] += band[start:stop:step].tolist()

        for spin in spins:
            new_projections[spin] = np.array(new_projections[spin])

        new_labels_dict = {label: point.frac_coords for label, point in bandstructure.labels_dict.items()}

        new_bandstructure = BandStructureSymmLine(
            kpoints=new_kpoints,
            eigenvals=new_bands,
            lattice=bandstructure.lattice_rec,
            efermi=bandstructure.efermi,
            labels_dict=new_labels_dict,
            structure=bandstructure.structure,
            projections=new_projections,
        )

        return new_bandstructure
