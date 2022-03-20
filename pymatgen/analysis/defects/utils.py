# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Utilities for defects module.
"""

import itertools
import logging
import math
import operator
from collections import defaultdict
from copy import deepcopy

import numpy as np
import pandas as pd
from monty.dev import requires
from monty.json import MSONable
from numpy.linalg import norm
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import Voronoi
from scipy.spatial.distance import squareform

from pymatgen.analysis.local_env import (
    LocalStructOrderParams,
    MinimumDistanceNN,
    cn_opt_params,
)
from pymatgen.analysis.phase_diagram import get_facets
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.periodic_table import Element, get_el_sp
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import pbc_diff
from pymatgen.vis.structure_vtk import StructureVis

try:
    from skimage.feature import peak_local_max

    peak_local_max_found = True
except ImportError:
    peak_local_max_found = False

__author__ = "Danny Broberg, Shyam Dwaraknath, Bharat Medasani, Nils Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Danny Broberg, Shyam Dwaraknath"
__email__ = "dbroberg@berkeley.edu, shyamd@lbl.gov"
__status__ = "Development"
__date__ = "January 11, 2018"

logger = logging.getLogger(__name__)
hart_to_ev = 27.2114
ang_to_bohr = 1.8897
invang_to_ev = 3.80986
kumagai_to_V = 1.809512739e2  # = Electron charge * 1e10 / VacuumPermittivity Constant

motif_cn_op = {}
for cn, di in cn_opt_params.items():  # type: ignore
    for mot, li in di.items():
        motif_cn_op[mot] = {"cn": int(cn), "optype": li[0]}
        motif_cn_op[mot]["params"] = deepcopy(li[1]) if len(li) > 1 else None


class QModel(MSONable):
    """
    Model for the defect charge distribution.
    A combination of exponential tail and gaussian distribution is used
    (see Freysoldt (2011), DOI: 10.1002/pssb.201046289 )
    q_model(r) = q [x exp(-r/gamma) + (1-x) exp(-r^2/beta^2)]
            without normalization constants
    By default, gaussian distribution with 1 Bohr width is assumed.
    If defect charge is more delocalized, exponential tail is suggested.
    """

    def __init__(self, beta=1.0, expnorm=0.0, gamma=1.0):
        """
        Args:
            beta: Gaussian decay constant. Default value is 1 Bohr.
                  When delocalized (eg. diamond), 2 Bohr is more appropriate.
            expnorm: Weight for the exponential tail in the range of [0-1].
                     Default is 0.0 indicating no tail .
                     For delocalized charges ideal value is around 0.54-0.6.
            gamma: Exponential decay constant
        """
        self.beta = beta
        self.expnorm = expnorm
        self.gamma = gamma

        self.beta2 = beta * beta
        self.gamma2 = gamma * gamma
        if expnorm and not gamma:
            raise ValueError("Please supply exponential decay constant.")

    def rho_rec(self, g2):
        """
        Reciprocal space model charge value
        for input squared reciprocal vector.
        Args:
            g2: Square of reciprocal vector

        Returns:
            Charge density at the reciprocal vector magnitude
        """
        return self.expnorm / np.sqrt(1 + self.gamma2 * g2) + (1 - self.expnorm) * np.exp(-0.25 * self.beta2 * g2)

    @property
    def rho_rec_limit0(self):
        """
        Reciprocal space model charge value
        close to reciprocal vector 0 .
        rho_rec(g->0) -> 1 + rho_rec_limit0 * g^2
        """
        return -2 * self.gamma2 * self.expnorm - 0.25 * self.beta2 * (1 - self.expnorm)


def eV_to_k(energy):
    """
    Convert energy to reciprocal vector magnitude k via hbar*k^2/2m
    Args:
        a: Energy in eV.

    Returns:
        (double) Reciprocal vector magnitude (units of 1/Bohr).
    """
    return math.sqrt(energy / invang_to_ev) * ang_to_bohr


def genrecip(a1, a2, a3, encut):
    """
    Args:
        a1, a2, a3: lattice vectors in bohr
        encut: energy cut off in eV
    Returns:
        reciprocal lattice vectors with energy less than encut
    """
    vol = np.dot(a1, np.cross(a2, a3))  # 1/bohr^3
    b1 = (2 * np.pi / vol) * np.cross(a2, a3)  # units 1/bohr
    b2 = (2 * np.pi / vol) * np.cross(a3, a1)
    b3 = (2 * np.pi / vol) * np.cross(a1, a2)

    # create list of recip space vectors that satisfy |i*b1+j*b2+k*b3|<=encut
    G_cut = eV_to_k(encut)
    # Figure out max in all recipricol lattice directions
    i_max = int(math.ceil(G_cut / norm(b1)))
    j_max = int(math.ceil(G_cut / norm(b2)))
    k_max = int(math.ceil(G_cut / norm(b3)))

    # Build index list
    i = np.arange(-i_max, i_max)
    j = np.arange(-j_max, j_max)
    k = np.arange(-k_max, k_max)

    # Convert index to vectors using meshgrid
    indices = np.array(np.meshgrid(i, j, k)).T.reshape(-1, 3)
    # Multiply integer vectors to get recipricol space vectors
    vecs = np.dot(indices, [b1, b2, b3])
    # Calculate radii of all vectors
    radii = np.sqrt(np.einsum("ij,ij->i", vecs, vecs))

    # Yield based on radii
    for vec, r in zip(vecs, radii):
        if r < G_cut and r != 0:
            yield vec


def generate_reciprocal_vectors_squared(a1, a2, a3, encut):
    """
    Generate reciprocal vector magnitudes within the cutoff along the specified
    lattice vectors.
    Args:
        a1: Lattice vector a (in Bohrs)
        a2: Lattice vector b (in Bohrs)
        a3: Lattice vector c (in Bohrs)
        encut: Reciprocal vector energy cutoff

    Returns:
        [[g1^2], [g2^2], ...] Square of reciprocal vectors (1/Bohr)^2
        determined by a1, a2, a3 and whose magntidue is less than gcut^2.
    """
    for vec in genrecip(a1, a2, a3, encut):
        yield np.dot(vec, vec)


def closestsites(struct_blk, struct_def, pos):
    """
    Returns closest site to the input position
    for both bulk and defect structures
    Args:
        struct_blk: Bulk structure
        struct_def: Defect structure
        pos: Position
    Return: (site object, dist, index)
    """
    blk_close_sites = struct_blk.get_sites_in_sphere(pos, 5, include_index=True)
    blk_close_sites.sort(key=lambda x: x[1])
    def_close_sites = struct_def.get_sites_in_sphere(pos, 5, include_index=True)
    def_close_sites.sort(key=lambda x: x[1])

    return blk_close_sites[0], def_close_sites[0]


class StructureMotifInterstitial:
    """
    Generate interstitial sites at positions
    where the interstitialcy is coordinated by nearest neighbors
    in a way that resembles basic structure motifs
    (e.g., tetrahedra, octahedra).  The algorithm is called InFiT
    (Interstitialcy Finding Tool), it was introducted by
    Nils E. R. Zimmermann, Matthew K. Horton, Anubhav Jain,
    and Maciej Haranczyk (Front. Mater., 4, 34, 2017),
    and it is used by the Python Charged Defect Toolkit
    (PyCDT: D. Broberg et al., Comput. Phys. Commun., in press, 2018).
    """

    def __init__(
        self,
        struct,
        inter_elem,
        motif_types=("tetrahedral", "octahedral"),
        op_threshs=(0.3, 0.5),
        dl=0.2,
        doverlap=1,
        facmaxdl=1.01,
        verbose=False,
    ):
        """
        Generates symmetrically distinct interstitial sites at positions
        where the interstitial is coordinated by nearest neighbors
        in a pattern that resembles a supported structure motif
        (e.g., tetrahedra, octahedra).

        Args:
            struct (Structure): input structure for which symmetrically
                distinct interstitial sites are to be found.
            inter_elem (string): element symbol of desired interstitial.
            motif_types ([string]): list of structure motif types that are
                to be considered.  Permissible types are:
                tet (tetrahedron), oct (octahedron).
            op_threshs ([float]): threshold values for the underlying order
                parameters to still recognize a given structural motif
                (i.e., for an OP value >= threshold the coordination pattern
                match is positive, for OP < threshold the match is
                negative.
            dl (float): grid fineness in Angstrom.  The input
                structure is divided into a grid of dimension
                a/dl x b/dl x c/dl along the three crystallographic
                directions, with a, b, and c being the lengths of
                the three lattice vectors of the input unit cell.
            doverlap (float): distance that is considered
                to flag an overlap between any trial interstitial site
                and a host atom.
            facmaxdl (float): factor to be multiplied with the maximum grid
                width that is then used as a cutoff distance for the
                clustering prune step.
            verbose (bool): flag indicating whether (True) or not (False;
                default) to print additional information to screen.
        """
        # Initialize interstitial finding.
        self._structure = struct.copy()
        self._motif_types = motif_types[:]
        if len(self._motif_types) == 0:
            raise RuntimeError("no motif types provided.")
        self._op_threshs = op_threshs[:]
        self.cn_motif_lostop = {}
        self.target_cns = []
        for motif in self._motif_types:
            if motif not in list(motif_cn_op.keys()):
                raise RuntimeError(f"unsupported motif type: {motif}.")
            cn = int(motif_cn_op[motif]["cn"])
            if cn not in self.target_cns:
                self.target_cns.append(cn)
            if cn not in list(self.cn_motif_lostop.keys()):
                self.cn_motif_lostop[cn] = {}
            tmp_optype = motif_cn_op[motif]["optype"]
            if tmp_optype == "tet_max":
                tmp_optype = "tet"
            if tmp_optype == "oct_max":
                tmp_optype = "oct"
            self.cn_motif_lostop[cn][motif] = LocalStructOrderParams(
                [tmp_optype], parameters=[motif_cn_op[motif]["params"]], cutoff=-10.0
            )
        self._dl = dl
        self._defect_sites = []
        self._defect_types = []
        self._defect_site_multiplicity = []
        self._defect_cns = []
        self._defect_opvals = []

        rots, trans = SpacegroupAnalyzer(struct)._get_symmetry()
        nbins = [
            int(struct.lattice.a / dl),
            int(struct.lattice.b / dl),
            int(struct.lattice.c / dl),
        ]
        dls = [
            struct.lattice.a / float(nbins[0]),
            struct.lattice.b / float(nbins[1]),
            struct.lattice.c / float(nbins[2]),
        ]
        maxdl = max(dls)
        if verbose:
            print(f"Grid size: {nbins[0]} {nbins[1]} {nbins[2]}")
            print(f"dls: {dls[0]} {dls[1]} {dls[2]}")
        struct_w_inter = struct.copy()
        struct_w_inter.append(inter_elem, [0, 0, 0])
        natoms = len(list(struct_w_inter.sites))
        trialsites = []

        # Build index list
        i = np.arange(0, nbins[0]) + 0.5
        j = np.arange(0, nbins[1]) + 0.5
        k = np.arange(0, nbins[2]) + 0.5

        # Convert index to vectors using meshgrid
        indices = np.array(np.meshgrid(i, j, k)).T.reshape(-1, 3)
        # Multiply integer vectors to get recipricol space vectors
        vecs = np.multiply(indices, np.divide(1, nbins))

        # Loop over trial positions that are based on a regular
        # grid in fractional coordinate space
        # within the unit cell.
        for vec in vecs:
            struct_w_inter.replace(natoms - 1, inter_elem, coords=vec, coords_are_cartesian=False)
            if len(struct_w_inter.get_sites_in_sphere(struct_w_inter.sites[natoms - 1].coords, doverlap)) == 1:
                neighs_images_weigths = MinimumDistanceNN(tol=0.8, cutoff=6).get_nn_info(struct_w_inter, natoms - 1)
                neighs_images_weigths_sorted = sorted(neighs_images_weigths, key=lambda x: x["weight"], reverse=True)
                for nsite in range(1, len(neighs_images_weigths_sorted) + 1):
                    if nsite not in self.target_cns:
                        continue

                    allsites = [neighs_images_weigths_sorted[i]["site"] for i in range(nsite)]
                    indices_neighs = list(range(len(allsites)))
                    allsites.append(struct_w_inter.sites[natoms - 1])
                    for mot, ops in self.cn_motif_lostop[nsite].items():
                        opvals = ops.get_order_parameters(allsites, len(allsites) - 1, indices_neighs=indices_neighs)
                        if opvals[0] > op_threshs[motif_types.index(mot)]:
                            cns = {}
                            for isite in range(nsite):
                                site = neighs_images_weigths_sorted[isite]["site"]
                                if isinstance(site.specie, Element):
                                    elem = site.specie.symbol
                                else:
                                    elem = site.specie.element.symbol
                                if elem in list(cns.keys()):
                                    cns[elem] = cns[elem] + 1
                                else:
                                    cns[elem] = 1
                            trialsites.append(
                                {
                                    "mtype": mot,
                                    "opval": opvals[0],
                                    "coords": struct_w_inter.sites[natoms - 1].coords[:],
                                    "fracs": vec,
                                    "cns": dict(cns),
                                }
                            )
                            break

        # Prune list of trial sites by clustering and find the site
        # with the largest order parameter value in each cluster.
        nintersites = len(trialsites)
        unique_motifs = []
        for ts in trialsites:
            if ts["mtype"] not in unique_motifs:
                unique_motifs.append(ts["mtype"])
        labels = {}
        connected = []
        for i in range(nintersites):
            connected.append([])
            for j in range(nintersites):
                dist, image = struct_w_inter.lattice.get_distance_and_image(
                    trialsites[i]["fracs"], trialsites[j]["fracs"]
                )
                connected[i].append(bool(dist < (maxdl * facmaxdl)))
        include = []
        for motif in unique_motifs:
            labels[motif] = []
            for i, ts in enumerate(trialsites):
                labels[motif].append(i if ts["mtype"] == motif else -1)
            change = True
            while change:
                change = False
                for i in range(nintersites - 1):
                    if change:
                        break
                    if labels[motif][i] == -1:
                        continue
                    for j in range(i + 1, nintersites):
                        if labels[motif][j] == -1:
                            continue
                        if connected[i][j] and labels[motif][i] != labels[motif][j]:
                            if labels[motif][i] < labels[motif][j]:
                                labels[motif][j] = labels[motif][i]
                            else:
                                labels[motif][i] = labels[motif][j]
                            change = True
                            break
            unique_ids = []
            for l in labels[motif]:
                if l != -1 and l not in unique_ids:
                    unique_ids.append(l)
            if verbose:
                print(f"unique_ids {motif} {unique_ids}")
            for uid in unique_ids:
                maxq = 0.0
                imaxq = -1
                for i in range(nintersites):
                    if labels[motif][i] == uid:
                        if imaxq < 0 or trialsites[i]["opval"] > maxq:
                            imaxq = i
                            maxq = trialsites[i]["opval"]
                include.append(imaxq)

        # Prune by symmetry.
        multiplicity = {}
        discard = []
        for motif in unique_motifs:
            discard_motif = []
            for indi, i in enumerate(include):
                if trialsites[i]["mtype"] != motif or i in discard_motif:
                    continue
                multiplicity[i] = 1
                symposlist = [trialsites[i]["fracs"].dot(np.array(m, dtype=float)) for m in rots]
                for t in trans:
                    symposlist.append(trialsites[i]["fracs"] + np.array(t))
                for indj in range(indi + 1, len(include)):
                    j = include[indj]
                    if trialsites[j]["mtype"] != motif or j in discard_motif:
                        continue
                    for sympos in symposlist:
                        dist, image = struct.lattice.get_distance_and_image(sympos, trialsites[j]["fracs"])
                        if dist < maxdl * facmaxdl:
                            discard_motif.append(j)
                            multiplicity[i] += 1
                            break
            for i in discard_motif:
                if i not in discard:
                    discard.append(i)

        if verbose:
            print(
                "Initial trial sites: {}\nAfter clustering: {}\n"
                "After symmetry pruning: {}".format(len(trialsites), len(include), len(include) - len(discard))
            )
        for i in include:
            if i not in discard:
                self._defect_sites.append(
                    PeriodicSite(
                        Element(inter_elem),
                        trialsites[i]["fracs"],
                        self._structure.lattice,
                        to_unit_cell=False,
                        coords_are_cartesian=False,
                        properties=None,
                    )
                )
                self._defect_types.append(trialsites[i]["mtype"])
                self._defect_cns.append(trialsites[i]["cns"])
                self._defect_site_multiplicity.append(multiplicity[i])
                self._defect_opvals.append(trialsites[i]["opval"])

    def enumerate_defectsites(self):
        """
        Get all defect sites.

        Returns:
            defect_sites ([PeriodicSite]): list of periodic sites
                    representing the interstitials.
        """
        return self._defect_sites

    def get_motif_type(self, i):
        """
        Get the motif type of defect with index i (e.g., "tet").

        Returns:
            motif (string): motif type.
        """
        return self._defect_types[i]

    def get_defectsite_multiplicity(self, n):
        """
        Returns the symmtric multiplicity of the defect site at the index.
        """
        return self._defect_site_multiplicity[n]

    def get_coordinating_elements_cns(self, i):
        """
        Get element-specific coordination numbers of defect with index i.

        Returns:
            elem_cn (dict): dictionary storing the coordination numbers (int)
                    with string representation of elements as keys.
                    (i.e., {elem1 (string): cn1 (int), ...}).
        """
        return self._defect_cns[i]

    def get_op_value(self, i):
        """
        Get order-parameter value of defect with index i.

        Returns:
            opval (float): OP value.
        """
        return self._defect_opvals[i]

    def make_supercells_with_defects(self, scaling_matrix):
        """
        Generate a sequence of supercells
        in which each supercell contains a single interstitial,
        except for the first supercell in the sequence
        which is a copy of the defect-free input structure.

        Args:
            scaling_matrix (3x3 integer array): scaling matrix
                to transform the lattice vectors.
        Returns:
            scs ([Structure]): sequence of supercells.

        """
        scs = []
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        scs.append(sc)
        for ids, defect_site in enumerate(self._defect_sites):
            sc_with_inter = sc.copy()
            sc_with_inter.append(
                defect_site.species_string,
                defect_site.frac_coords,
                coords_are_cartesian=False,
                validate_proximity=False,
                properties=None,
            )
            if not sc_with_inter:
                raise RuntimeError(f"could not generate supercell with interstitial {ids + 1}")
            scs.append(sc_with_inter.copy())
        return scs


class TopographyAnalyzer:
    """
    This is a generalized module to perform topological analyses of a crystal
    structure using Voronoi tessellations. It can be used for finding potential
    interstitial sites. Applications including using these sites for
    inserting additional atoms or for analyzing diffusion pathways.

    Note that you typically want to do some preliminary postprocessing after
    the initial construction. The initial construction will create a lot of
    points, especially for determining potential insertion sites. Some helper
    methods are available to perform aggregation and elimination of nodes. A
    typical use is something like::

        a = TopographyAnalyzer(structure, ["O"], ["P"])
        a.cluster_nodes()
        a.remove_collisions()
    """

    def __init__(
        self,
        structure,
        framework_ions,
        cations,
        tol=0.0001,
        max_cell_range=1,
        check_volume=True,
        constrained_c_frac=0.5,
        thickness=0.5,
    ):
        """
        Init.

        Args:
            structure (Structure): An initial structure.
            framework_ions ([str]): A list of ions to be considered as a
                framework. Typically, this would be all anion species. E.g.,
                ["O", "S"].
            cations ([str]): A list of ions to be considered as non-migrating
                cations. E.g., if you are looking at Li3PS4 as a Li
                conductor, Li is a mobile species. Your cations should be [
                "P"]. The cations are used to exclude polyhedra from
                diffusion analysis since those polyhedra are already occupied.
            tol (float): A tolerance distance for the analysis, used to
                determine if something are actually periodic boundary images of
                each other. Default is usually fine.
            max_cell_range (int): This is the range of periodic images to
                construct the Voronoi tessellation. A value of 1 means that we
                include all points from (x +- 1, y +- 1, z+- 1) in the
                voronoi construction. This is because the Voronoi poly
                extends beyond the standard unit cell because of PBC.
                Typically, the default value of 1 works fine for most
                structures and is fast. But for really small unit
                cells with high symmetry, you may need to increase this to 2
                or higher.
            check_volume (bool): Set False when ValueError always happen after
                tuning tolerance.
            constrained_c_frac (float): Constraint the region where users want
                to do Topology analysis the default value is 0.5, which is the
                fractional coordinate of the cell
            thickness (float): Along with constrained_c_frac, limit the
                thickness of the regions where we want to explore. Default is
                0.5, which is mapping all the site of the unit cell.

        """
        self.structure = structure
        self.framework_ions = {get_el_sp(sp) for sp in framework_ions}
        self.cations = {get_el_sp(sp) for sp in cations}

        # Let us first map all sites to the standard unit cell, i.e.,
        # 0 ≤ coordinates < 1.
        # structure = Structure.from_sites(structure, to_unit_cell=True)
        # lattice = structure.lattice

        # We could constrain the region where we want to dope/explore by setting
        # the value of constrained_c_frac and thickness. The default mode is
        # mapping all sites to the standard unit cell
        s = structure.copy()
        constrained_sites = []
        for i, site in enumerate(s):
            if (
                site.frac_coords[2] >= constrained_c_frac - thickness
                and site.frac_coords[2] <= constrained_c_frac + thickness
            ):
                constrained_sites.append(site)
        structure = Structure.from_sites(sites=constrained_sites)
        lattice = structure.lattice

        # Divide the sites into framework and non-framework sites.
        framework = []
        non_framework = []
        for site in structure:
            if self.framework_ions.intersection(site.species.keys()):
                framework.append(site)
            else:
                non_framework.append(site)

        # We construct a supercell series of coords. This is because the
        # Voronoi polyhedra can extend beyond the standard unit cell. Using a
        # range of -2, -1, 0, 1 should be fine.
        coords = []
        cell_range = list(range(-max_cell_range, max_cell_range + 1))
        for shift in itertools.product(cell_range, cell_range, cell_range):
            for site in framework:
                shifted = site.frac_coords + shift
                coords.append(lattice.get_cartesian_coords(shifted))

        # Perform the voronoi tessellation.
        voro = Voronoi(coords)

        # Store a mapping of each voronoi node to a set of points.
        node_points_map = defaultdict(set)
        for pts, vs in voro.ridge_dict.items():
            for v in vs:
                node_points_map[v].update(pts)

        logger.debug(f"{len(voro.vertices)} total Voronoi vertices")

        # Vnodes store all the valid voronoi polyhedra. Cation vnodes store
        # the voronoi polyhedra that are already occupied by existing cations.
        vnodes = []
        cation_vnodes = []

        def get_mapping(poly):
            """
            Helper function to check if a vornoi poly is a periodic image
            of one of the existing voronoi polys.
            """
            for v in vnodes:
                if v.is_image(poly, tol):
                    return v
            return None

        # Filter all the voronoi polyhedra so that we only consider those
        # which are within the unit cell.
        for i, vertex in enumerate(voro.vertices):
            if i == 0:
                continue
            fcoord = lattice.get_fractional_coords(vertex)
            poly = VoronoiPolyhedron(lattice, fcoord, node_points_map[i], coords, i)
            if np.all([-tol <= c < 1 + tol for c in fcoord]):
                if len(vnodes) == 0:
                    vnodes.append(poly)
                else:
                    ref = get_mapping(poly)
                    if ref is None:
                        vnodes.append(poly)

        logger.debug(f"{len(vnodes)} voronoi vertices in cell.")

        # Eliminate all voronoi nodes which are closest to existing cations.
        if len(cations) > 0:
            cation_coords = [
                site.frac_coords for site in non_framework if self.cations.intersection(site.species.keys())
            ]

            vertex_fcoords = [v.frac_coords for v in vnodes]
            dist_matrix = lattice.get_all_distances(cation_coords, vertex_fcoords)
            indices = np.where(dist_matrix == np.min(dist_matrix, axis=1)[:, None])[1]
            cation_vnodes = [v for i, v in enumerate(vnodes) if i in indices]
            vnodes = [v for i, v in enumerate(vnodes) if i not in indices]

        logger.debug(f"{len(vnodes)} vertices in cell not with cation.")
        self.coords = coords
        self.vnodes = vnodes
        self.cation_vnodes = cation_vnodes
        self.framework = framework
        self.non_framework = non_framework
        if check_volume:
            self.check_volume()

    def check_volume(self):
        """
        Basic check for volume of all voronoi poly sum to unit cell volume
        Note that this does not apply after poly combination.
        """
        vol = sum(v.volume for v in self.vnodes) + sum(v.volume for v in self.cation_vnodes)
        if abs(vol - self.structure.volume) > 1e-8:
            raise ValueError(
                "Sum of voronoi volumes is not equal to original volume of "
                "structure! This may lead to inaccurate results. You need to "
                "tweak the tolerance and max_cell_range until you get a "
                "correct mapping."
            )

    def cluster_nodes(self, tol=0.2):
        """
        Cluster nodes that are too close together using a tol.

        Args:
            tol (float): A distance tolerance. PBC is taken into account.
        """
        lattice = self.structure.lattice

        vfcoords = [v.frac_coords for v in self.vnodes]

        # Manually generate the distance matrix (which needs to take into
        # account PBC.
        dist_matrix = np.array(lattice.get_all_distances(vfcoords, vfcoords))
        dist_matrix = (dist_matrix + dist_matrix.T) / 2
        for i in range(len(dist_matrix)):
            dist_matrix[i, i] = 0
        condensed_m = squareform(dist_matrix)
        z = linkage(condensed_m)
        cn = fcluster(z, tol, criterion="distance")
        merged_vnodes = []
        for n in set(cn):
            poly_indices = set()
            frac_coords = []
            for i, j in enumerate(np.where(cn == n)[0]):
                poly_indices.update(self.vnodes[j].polyhedron_indices)
                if i == 0:
                    frac_coords.append(self.vnodes[j].frac_coords)
                else:
                    fcoords = self.vnodes[j].frac_coords
                    # We need the image to combine the frac_coords properly.
                    d, image = lattice.get_distance_and_image(frac_coords[0], fcoords)
                    frac_coords.append(fcoords + image)
            merged_vnodes.append(VoronoiPolyhedron(lattice, np.average(frac_coords, axis=0), poly_indices, self.coords))
        self.vnodes = merged_vnodes
        logger.debug(f"{len(self.vnodes)} vertices after combination.")

    def remove_collisions(self, min_dist=0.5):
        """
        Remove vnodes that are too close to existing atoms in the structure

        Args:
            min_dist(float): The minimum distance that a vertex needs to be
                from existing atoms.
        """
        vfcoords = [v.frac_coords for v in self.vnodes]
        sfcoords = self.structure.frac_coords
        dist_matrix = self.structure.lattice.get_all_distances(vfcoords, sfcoords)
        all_dist = np.min(dist_matrix, axis=1)
        new_vnodes = []
        for i, v in enumerate(self.vnodes):
            if all_dist[i] > min_dist:
                new_vnodes.append(v)
        self.vnodes = new_vnodes

    def get_structure_with_nodes(self):
        """
        Get the modified structure with the voronoi nodes inserted. The
        species is set as a DummySpecies X.
        """
        new_s = Structure.from_sites(self.structure)
        for v in self.vnodes:
            new_s.append("X", v.frac_coords)
        return new_s

    def print_stats(self):
        """
        Print stats such as the MSE dist.
        """
        latt = self.structure.lattice

        def get_min_dist(fcoords):
            n = len(fcoords)
            dist = latt.get_all_distances(fcoords, fcoords)
            all_dist = [dist[i, j] for i in range(n) for j in range(i + 1, n)]
            return min(all_dist)

        voro = [s.frac_coords for s in self.vnodes]
        print(f"Min dist between voronoi vertices centers = {get_min_dist(voro):.4f}")

        def get_non_framework_dist(fcoords):
            cations = [site.frac_coords for site in self.non_framework]
            dist_matrix = latt.get_all_distances(cations, fcoords)
            min_dist = np.min(dist_matrix, axis=1)
            if len(cations) != len(min_dist):
                raise Exception("Could not calculate distance to all cations")
            return np.linalg.norm(min_dist), min(min_dist), max(min_dist)

        print(len(self.non_framework))
        print(f"MSE dist voro = {str(get_non_framework_dist(voro))}")

    def write_topology(self, fname="Topo.cif"):
        """
        Write topology to a file.

        :param fname: Filename
        """
        new_s = Structure.from_sites(self.structure)
        for v in self.vnodes:
            new_s.append("Mg", v.frac_coords)
        new_s.to(filename=fname)

    def analyze_symmetry(self, tol):
        """
        :param tol: Tolerance for SpaceGroupAnalyzer
        :return: List
        """
        s = Structure.from_sites(self.framework)
        site_to_vindex = {}
        for i, v in enumerate(self.vnodes):
            s.append("Li", v.frac_coords)
            site_to_vindex[s[-1]] = i

        print(len(s))
        finder = SpacegroupAnalyzer(s, tol)
        print(finder.get_space_group_operations())
        symm_structure = finder.get_symmetrized_structure()
        print(len(symm_structure.equivalent_sites))
        return [
            [site_to_vindex[site] for site in sites]
            for sites in symm_structure.equivalent_sites
            if sites[0].specie.symbol == "Li"
        ]

    def vtk(self):
        """
        Show VTK visualization.
        """
        if StructureVis is None:
            raise NotImplementedError("vtk must be present to view.")
        lattice = self.structure.lattice
        vis = StructureVis()
        vis.set_structure(Structure.from_sites(self.structure))
        for v in self.vnodes:
            vis.add_site(PeriodicSite("K", v.frac_coords, lattice))
            vis.add_polyhedron(
                [PeriodicSite("S", c, lattice, coords_are_cartesian=True) for c in v.polyhedron_coords],
                PeriodicSite("Na", v.frac_coords, lattice),
                color="element",
                draw_edges=True,
                edges_color=(0, 0, 0),
            )
        vis.show()


class VoronoiPolyhedron:
    """
    Convenience container for a voronoi point in PBC and its associated polyhedron.
    """

    def __init__(self, lattice, frac_coords, polyhedron_indices, all_coords, name=None):
        """
        :param lattice:
        :param frac_coords:
        :param polyhedron_indices:
        :param all_coords:
        :param name:
        """
        self.lattice = lattice
        self.frac_coords = frac_coords
        self.polyhedron_indices = polyhedron_indices
        self.polyhedron_coords = np.array(all_coords)[list(polyhedron_indices), :]
        self.name = name

    def is_image(self, poly, tol):
        """
        :param poly: VoronoiPolyhedron
        :param tol: Coordinate tolerance.
        :return: Whether a poly is an image of the current one.
        """
        frac_diff = pbc_diff(poly.frac_coords, self.frac_coords)
        if not np.allclose(frac_diff, [0, 0, 0], atol=tol):
            return False
        to_frac = self.lattice.get_fractional_coords
        for c1 in self.polyhedron_coords:
            found = False
            for c2 in poly.polyhedron_coords:
                d = pbc_diff(to_frac(c1), to_frac(c2))
                if not np.allclose(d, [0, 0, 0], atol=tol):
                    found = True
                    break
            if not found:
                return False
        return True

    @property
    def coordination(self):
        """
        :return: Coordination number
        """
        return len(self.polyhedron_indices)

    @property
    def volume(self):
        """
        :return: Volume
        """
        return calculate_vol(self.polyhedron_coords)

    def __str__(self):
        return f"Voronoi polyhedron {self.name}"


class ChargeDensityAnalyzer(MSONable):
    """
    Analyzer to find potential interstitial sites based on charge density. The
    `total` charge density is used.
    """

    def __init__(self, chgcar):
        """
        Initialization.

        Args:
            chgcar (pmg.Chgcar): input Chgcar object.
        """
        self.chgcar = chgcar
        self.structure = chgcar.structure
        self.extrema_coords = []  # list of frac_coords of local extrema
        self.extrema_type = None  # "local maxima" or "local minima"
        self._extrema_df = None  # extrema frac_coords - chg density table
        self._charge_distribution_df = None  # frac_coords - chg density table

    @classmethod
    def from_file(cls, chgcar_filename):
        """
        Init from a CHGCAR.

        :param chgcar_filename:
        :return:
        """
        chgcar = Chgcar.from_file(chgcar_filename)
        return cls(chgcar=chgcar)

    @property
    def charge_distribution_df(self):
        """
        :return: Charge distribution.
        """
        if self._charge_distribution_df is None:
            return self._get_charge_distribution_df()
        return self._charge_distribution_df

    @property
    def extrema_df(self):
        """
        :return: The extrema in charge density.
        """
        if self.extrema_type is None:
            logger.warning("Please run ChargeDensityAnalyzer.get_local_extrema first!")
        return self._extrema_df

    def _get_charge_distribution_df(self):
        """
        Return a complete table of fractional coordinates - charge density.
        """
        # Fraction coordinates and corresponding indices
        axis_grid = np.array([np.array(self.chgcar.get_axis_grid(i)) / self.structure.lattice.abc[i] for i in range(3)])
        axis_index = np.array([range(len(axis_grid[i])) for i in range(3)])

        data = {}

        for index in itertools.product(*axis_index):
            a, b, c = index
            f_coords = (axis_grid[0][a], axis_grid[1][b], axis_grid[2][c])
            data[f_coords] = self.chgcar.data["total"][a][b][c]

        # Fraction coordinates - charge density table
        df = pd.Series(data).reset_index()
        df.columns = ["a", "b", "c", "Charge Density"]
        self._charge_distribution_df = df

        return df

    def _update_extrema(self, f_coords, extrema_type, threshold_frac=None, threshold_abs=None):
        """Update _extrema_df, extrema_type and extrema_coords"""

        if threshold_frac is not None:
            if threshold_abs is not None:
                logger.warning("Filter can be either threshold_frac or threshold_abs!")  # Exit if both filter are set
                return
            if threshold_frac > 1 or threshold_frac < 0:
                raise Exception("threshold_frac range is [0, 1]!")

        # Return empty result if coords list is empty
        if len(f_coords) == 0:
            df = pd.DataFrame({}, columns=["A", "B", "C", "Chgcar"])
            self._extrema_df = df
            self.extrema_coords = []
            logger.info(f"Find {len(df)} {extrema_type}.")
            return

        data = {}
        unit = 1 / np.array(self.chgcar.dim)  # pixel along a, b, c

        for fc in f_coords:
            a, b, c = tuple(map(int, fc / unit))
            data[tuple(fc)] = self.chgcar.data["total"][a][b][c]

        df = pd.Series(data).reset_index()
        df.columns = ["a", "b", "c", "Charge Density"]
        ascending = extrema_type == "local minima"

        if threshold_abs is None:
            threshold_frac = threshold_frac if threshold_frac is not None else 1.0
            num_extrema = int(threshold_frac * len(f_coords))
            df = df.sort_values(by="Charge Density", ascending=ascending)[0:num_extrema]
            df.reset_index(drop=True, inplace=True)  # reset major index
        else:  # threshold_abs is set
            df = df.sort_values(by="Charge Density", ascending=ascending)
            df = df[df["Charge Density"] <= threshold_abs] if ascending else df[df["Charge Density"] >= threshold_abs]

        extrema_coords = []
        for row in df.iterrows():
            fc = np.array(row[1]["a":"c"])
            extrema_coords.append(fc)

        self._extrema_df = df
        self.extrema_type = extrema_type
        self.extrema_coords = extrema_coords
        logger.info(f"Find {len(df)} {extrema_type}.")

    @requires(
        peak_local_max_found,
        "get_local_extrema requires skimage.feature.peak_local_max module"
        " to be installed. Please confirm your skimage installation.",
    )
    def get_local_extrema(self, find_min=True, threshold_frac=None, threshold_abs=None):
        """
        Get all local extrema fractional coordinates in charge density,
        searching for local minimum by default. Note that sites are NOT grouped
        symmetrically.

        Args:
            find_min (bool): True to find local minimum else maximum, otherwise
                find local maximum.

            threshold_frac (float): optional fraction of extrema shown, which
                returns `threshold_frac * tot_num_extrema` extrema fractional
                coordinates based on highest/lowest intensity.

                E.g. set 0.2 to show the extrema with 20% highest or lowest
                intensity. Value range: 0 <= threshold_frac <= 1

                Note that threshold_abs and threshold_frac should not set in the
                same time.

            threshold_abs (float): optional filter. When searching for local
                minima, intensity <= threshold_abs returns; when searching for
                local maxima, intensity >= threshold_abs returns.

                Note that threshold_abs and threshold_frac should not set in the
                same time.

        Returns:
            extrema_coords (list): list of fractional coordinates corresponding
                to local extrema.
        """
        sign, extrema_type = 1, "local maxima"

        if find_min:
            sign, extrema_type = -1, "local minima"

        # Make 3x3x3 supercell
        # This is a trick to resolve the periodical boundary issue.
        total_chg = sign * self.chgcar.data["total"]
        total_chg = np.tile(total_chg, reps=(3, 3, 3))
        coordinates = peak_local_max(total_chg, min_distance=1)

        # Remove duplicated sites introduced by supercell.
        f_coords = [coord / total_chg.shape * 3 for coord in coordinates]
        f_coords = [f - 1 for f in f_coords if all(np.array(f) < 2) and all(np.array(f) >= 1)]

        # Update information
        self._update_extrema(
            f_coords,
            extrema_type,
            threshold_frac=threshold_frac,
            threshold_abs=threshold_abs,
        )

        return self.extrema_coords

    def cluster_nodes(self, tol=0.2):
        """
        Cluster nodes that are too close together using a tol.

        Args:
            tol (float): A distance tolerance. PBC is taken into account.
        """
        lattice = self.structure.lattice
        vf_coords = self.extrema_coords

        if len(vf_coords) == 0:
            if self.extrema_type is None:
                logger.warning("Please run ChargeDensityAnalyzer.get_local_extrema first!")
                return None
            new_f_coords = []
            self._update_extrema(new_f_coords, self.extrema_type)
            return new_f_coords

        # Manually generate the distance matrix (which needs to take into
        # account PBC.
        dist_matrix = np.array(lattice.get_all_distances(vf_coords, vf_coords))
        dist_matrix = (dist_matrix + dist_matrix.T) / 2

        for i in range(len(dist_matrix)):
            dist_matrix[i, i] = 0
        condensed_m = squareform(dist_matrix)
        z = linkage(condensed_m)
        cn = fcluster(z, tol, criterion="distance")
        merged_fcoords = []

        for n in set(cn):
            frac_coords = []
            for i, j in enumerate(np.where(cn == n)[0]):
                if i == 0:
                    frac_coords.append(self.extrema_coords[j])
                else:
                    f_coords = self.extrema_coords[j]
                    # We need the image to combine the frac_coords properly.
                    d, image = lattice.get_distance_and_image(frac_coords[0], f_coords)
                    frac_coords.append(f_coords + image)
            merged_fcoords.append(np.average(frac_coords, axis=0))

        merged_fcoords = [f - np.floor(f) for f in merged_fcoords]
        merged_fcoords = [f * (np.abs(f - 1) > 1e-15) for f in merged_fcoords]
        # the second line for fringe cases like
        # np.array([ 5.0000000e-01 -4.4408921e-17  5.0000000e-01])
        # where the shift to [0,1) does not work due to float precision
        self._update_extrema(merged_fcoords, extrema_type=self.extrema_type)
        logger.debug(f"{len(self.extrema_coords)} vertices after combination.")
        return None

    def remove_collisions(self, min_dist=0.5):
        """
        Remove predicted sites that are too close to existing atoms in the
        structure.

        Args:
            min_dist (float): The minimum distance (in Angstrom) that
                a predicted site needs to be from existing atoms. A min_dist
                with value <= 0 returns all sites without distance checking.
        """
        s_f_coords = self.structure.frac_coords
        f_coords = self.extrema_coords
        if len(f_coords) == 0:
            if self.extrema_type is None:
                logger.warning("Please run ChargeDensityAnalyzer.get_local_extrema first!")
                return None
            new_f_coords = []
            self._update_extrema(new_f_coords, self.extrema_type)
            return new_f_coords

        dist_matrix = self.structure.lattice.get_all_distances(f_coords, s_f_coords)
        all_dist = np.min(dist_matrix, axis=1)
        new_f_coords = []

        for i, f in enumerate(f_coords):
            if all_dist[i] > min_dist:
                new_f_coords.append(f)
        self._update_extrema(new_f_coords, self.extrema_type)

        return new_f_coords

    def get_structure_with_nodes(
        self,
        find_min=True,
        min_dist=0.5,
        tol=0.2,
        threshold_frac=None,
        threshold_abs=None,
    ):
        """
        Get the modified structure with the possible interstitial sites added.
        The species is set as a DummySpecies X.

        Args:
            find_min (bool): True to find local minimum else maximum, otherwise
                find local maximum.

            min_dist (float): The minimum distance (in Angstrom) that
                a predicted site needs to be from existing atoms. A min_dist
                with value <= 0 returns all sites without distance checking.

            tol (float): A distance tolerance of nodes clustering that sites too
                closed to other predicted sites will be merged. PBC is taken
                into account.

            threshold_frac (float): optional fraction of extrema, which returns
                `threshold_frac * tot_num_extrema` extrema fractional
                coordinates based on highest/lowest intensity.

                E.g. set 0.2 to insert DummySpecies atom at the extrema with 20%
                highest or lowest intensity.
                Value range: 0 <= threshold_frac <= 1

                Note that threshold_abs and threshold_frac should not set in the
                same time.

            threshold_abs (float): optional filter. When searching for local
                minima, intensity <= threshold_abs returns; when searching for
                local maxima, intensity >= threshold_abs returns.

                Note that threshold_abs and threshold_frac should not set in the
                same time.

        Returns:
            structure (Structure)
        """

        structure = self.structure.copy()
        self.get_local_extrema(
            find_min=find_min,
            threshold_frac=threshold_frac,
            threshold_abs=threshold_abs,
        )

        self.remove_collisions(min_dist)
        self.cluster_nodes(tol=tol)
        for fc in self.extrema_coords:
            structure.append("X", fc)

        return structure

    def sort_sites_by_integrated_chg(self, r=0.4):
        """
        Get the average charge density around each local minima in the charge density
        and store the result in _extrema_df
        Args:
            r (float): radius of sphere around each site to evaluate the average
        """

        if self.extrema_type is None:
            self.get_local_extrema()
        int_den = []
        for isite in self.extrema_coords:
            mask = self._dist_mat(isite) < r
            vol_sphere = self.chgcar.structure.volume * (mask.sum() / self.chgcar.ngridpts)
            chg_in_sphere = np.sum(self.chgcar.data["total"] * mask) / mask.size / vol_sphere
            int_den.append(chg_in_sphere)
        self._extrema_df["avg_charge_den"] = int_den
        self._extrema_df.sort_values(by=["avg_charge_den"], inplace=True)
        self._extrema_df.reset_index(drop=True, inplace=True)

    def _dist_mat(self, pos_frac):
        # return a matrix that contains the distances
        aa = np.linspace(0, 1, len(self.chgcar.get_axis_grid(0)), endpoint=False)
        bb = np.linspace(0, 1, len(self.chgcar.get_axis_grid(1)), endpoint=False)
        cc = np.linspace(0, 1, len(self.chgcar.get_axis_grid(2)), endpoint=False)
        AA, BB, CC = np.meshgrid(aa, bb, cc, indexing="ij")
        dist_from_pos = self.chgcar.structure.lattice.get_all_distances(
            fcoords1=np.vstack([AA.flatten(), BB.flatten(), CC.flatten()]).T,
            fcoords2=pos_frac,
        )
        return dist_from_pos.reshape(AA.shape)


class ChargeInsertionAnalyzer(ChargeDensityAnalyzer):
    """
    Analyze the charge density and create new candidate structures by inserting at each charge minima
    The similar inserterd structures are given the same uniqueness label.
    This works best with AECCAR data since CHGCAR data often contains spurious local minima in the core.
    However you can still use CHGCAR with an appropriate max_avg_charge value.

    Application of this for Li can be found at:
    J.-X. Shen et al.: npj Comput. Mater. 6, 1 (2020)
    https://www.nature.com/articles/s41524-020-00422-3
    """

    def __init__(
        self,
        chgcar,
        working_ion="Li",
        avg_radius=0.4,
        max_avg_charge=1.0,
        clustering_tol=0.6,
        ltol=0.2,
        stol=0.3,
        angle_tol=5,
    ):
        """
        Args:
            chgcar: The charge density object to analyze
            working_ion: The working ion to be inserted
            avg_radius: The radius used to calculate average charge density at each site
            max_avg_charge: Do no consider local minmas with avg charge above this value.
            clustering_tol: Distance tolerance for grouping sites together
            ltol: StructureMatcher ltol parameter
            stol: StructureMatcher stol parameter
            angle_tol: StructureMatcher angle_tol parameter
        """
        self.working_ion = working_ion
        self.sm = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol)
        self.max_avg_charge = max_avg_charge
        self.avg_radius = avg_radius
        self.clustering_tol = clustering_tol

        super().__init__(chgcar)

    def get_labels(self):
        """
        Populate the extrema dataframe (self._extrema_df) with the insertion structure.
        Then, group the sites by structure similarity.
        Finally store a full list of the insertion sites, with their labels as a Structure Object
        """

        self.get_local_extrema()

        if len(self._extrema_df) > 1:
            self.cluster_nodes(tol=self.clustering_tol)

        self.sort_sites_by_integrated_chg(r=self.avg_radius)

        inserted_structs = []

        self._extrema_df = self._extrema_df[self._extrema_df.avg_charge_den <= self.max_avg_charge]

        for itr, li_site in self._extrema_df.iterrows():
            if li_site["avg_charge_den"] > self.max_avg_charge:
                continue
            tmp_struct = self.chgcar.structure.copy()
            li_site = self._extrema_df.iloc[itr]
            tmp_struct.insert(
                0,
                self.working_ion,
                [li_site["a"], li_site["b"], li_site["c"]],
                properties=dict(magmom=0),
            )
            tmp_struct.sort()
            inserted_structs.append(tmp_struct)
        self._extrema_df["inserted_struct"] = inserted_structs
        site_labels = generic_groupby(self._extrema_df.inserted_struct, comp=self.sm.fit)
        self._extrema_df["site_label"] = site_labels

        # generate the structure with only Li atoms for NN analysis
        self.allsites_struct = Structure(
            self.structure.lattice,
            np.repeat(self.working_ion, len(self._extrema_df)),
            self._extrema_df[["a", "b", "c"]].values,
            site_properties={"label": self._extrema_df[["site_label"]].values.flatten()},
        )


def generic_groupby(list_in, comp=operator.eq):
    """
    Group a list of unsortable objects
    Args:
        list_in: A list of generic objects
        comp: (Default value = operator.eq) The comparator
    Returns:
        [int] list of labels for the input list
    """
    list_out = [None] * len(list_in)
    label_num = 0
    for i1, ls1 in enumerate(list_out):
        if ls1 is not None:
            continue
        list_out[i1] = label_num
        for i2, ls2 in list(enumerate(list_out))[(i1 + 1) :]:
            if comp(list_in[i1], list_in[i2]):
                if list_out[i2] is None:
                    list_out[i2] = list_out[i1]
                else:
                    list_out[i1] = list_out[i2]
                    label_num -= 1
        label_num += 1
    return list_out


def calculate_vol(coords):
    """
    Calculate volume given a set of coords.

    :param coords: List of coords.
    :return: Volume
    """
    if len(coords) == 4:
        coords_affine = np.ones((4, 4))
        coords_affine[:, 0:3] = np.array(coords)
        return abs(np.linalg.det(coords_affine)) / 6

    simplices = get_facets(coords, joggle=True)
    center = np.average(coords, axis=0)
    vol = 0
    for s in simplices:
        c = list(coords[i] for i in s)
        c.append(center)
        vol += calculate_vol(c)
    return vol


def converge(f, step, tol, max_h):
    """
    simple newton iteration based convergence function
    """
    g = f(0)
    dx = 10000
    h = step
    while dx > tol:
        g2 = f(h)
        dx = abs(g - g2)
        g = g2
        h += step

        if h > max_h:
            raise Exception(f"Did not converge before {h}")
    return g


def tune_for_gamma(lattice, epsilon):
    """
    This tunes the gamma parameter for Kumagai anisotropic
    Ewald calculation. Method is to find a gamma parameter which generates a similar
    number of reciprocal and real lattice vectors,
    given the suggested cut off radii by Kumagai and Oba
    """
    logger.debug("Converging for ewald parameter...")
    prec = 25  # a reasonable precision to tune gamma for

    gamma = (2 * np.average(lattice.abc)) ** (-1 / 2.0)
    recip_set, _, real_set, _ = generate_R_and_G_vecs(gamma, prec, lattice, epsilon)
    recip_set = recip_set[0]
    real_set = real_set[0]

    logger.debug(
        "First approach with gamma ={}\nProduced {} real vecs and {} recip "
        "vecs.".format(gamma, len(real_set), len(recip_set))
    )

    while float(len(real_set)) / len(recip_set) > 1.05 or float(len(recip_set)) / len(real_set) > 1.05:
        gamma *= (float(len(real_set)) / float(len(recip_set))) ** 0.17
        logger.debug(f"\tNot converged...Try modifying gamma to {gamma}.")
        recip_set, _, real_set, _ = generate_R_and_G_vecs(gamma, prec, lattice, epsilon)
        recip_set = recip_set[0]
        real_set = real_set[0]
        logger.debug(f"Now have {len(real_set)} real vecs and {len(recip_set)} recip vecs.")

    logger.debug(f"Converged with gamma = {gamma}")

    return gamma


def generate_R_and_G_vecs(gamma, prec_set, lattice, epsilon):
    """
    This returns a set of real and reciprocal lattice vectors
    (and real/recip summation values)
    based on a list of precision values (prec_set)

    gamma (float): Ewald parameter
    prec_set (list or number): for prec values to consider (20, 25, 30 are sensible numbers)
    lattice: Lattice object of supercell in question

    """
    if not isinstance(prec_set, list):
        prec_set = [prec_set]

    [a1, a2, a3] = lattice.matrix  # Angstrom
    volume = lattice.volume
    [b1, b2, b3] = lattice.reciprocal_lattice.matrix  # 1/ Angstrom
    invepsilon = np.linalg.inv(epsilon)
    rd_epsilon = np.sqrt(np.linalg.det(epsilon))

    # generate reciprocal vector set (for each prec_set)
    recip_set = [[] for prec in prec_set]
    recip_summation_values = [0.0 for prec in prec_set]
    recip_cut_set = [(2 * gamma * prec) for prec in prec_set]

    i_max = int(math.ceil(max(recip_cut_set) / np.linalg.norm(b1)))
    j_max = int(math.ceil(max(recip_cut_set) / np.linalg.norm(b2)))
    k_max = int(math.ceil(max(recip_cut_set) / np.linalg.norm(b3)))
    for i in np.arange(-i_max, i_max + 1):
        for j in np.arange(-j_max, j_max + 1):
            for k in np.arange(-k_max, k_max + 1):
                if not i and not j and not k:
                    continue
                gvec = i * b1 + j * b2 + k * b3
                normgvec = np.linalg.norm(gvec)
                for recip_cut_ind, recip_cut in enumerate(recip_cut_set):
                    if normgvec <= recip_cut:
                        recip_set[recip_cut_ind].append(gvec)

                        Gdotdiel = np.dot(gvec, np.dot(epsilon, gvec))
                        summand = math.exp(-Gdotdiel / (4 * (gamma**2))) / Gdotdiel
                        recip_summation_values[recip_cut_ind] += summand

    recip_summation_values = np.array(recip_summation_values)
    recip_summation_values /= volume

    # generate real vector set (for each prec_set)
    real_set = [[] for prec in prec_set]
    real_summation_values = [0.0 for prec in prec_set]
    real_cut_set = [(prec / gamma) for prec in prec_set]

    i_max = int(math.ceil(max(real_cut_set) / np.linalg.norm(a1)))
    j_max = int(math.ceil(max(real_cut_set) / np.linalg.norm(a2)))
    k_max = int(math.ceil(max(real_cut_set) / np.linalg.norm(a3)))
    for i in np.arange(-i_max, i_max + 1):
        for j in np.arange(-j_max, j_max + 1):
            for k in np.arange(-k_max, k_max + 1):
                rvec = i * a1 + j * a2 + k * a3
                normrvec = np.linalg.norm(rvec)
                for real_cut_ind, real_cut in enumerate(real_cut_set):
                    if normrvec <= real_cut:
                        real_set[real_cut_ind].append(rvec)
                        if normrvec > 1e-8:
                            sqrt_loc_res = np.sqrt(np.dot(rvec, np.dot(invepsilon, rvec)))
                            nmr = math.erfc(gamma * sqrt_loc_res)
                            real_summation_values[real_cut_ind] += nmr / sqrt_loc_res

    real_summation_values = np.array(real_summation_values)
    real_summation_values /= 4 * np.pi * rd_epsilon

    return recip_set, recip_summation_values, real_set, real_summation_values
