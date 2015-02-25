# coding: utf-8

from __future__ import division, unicode_literals
from functools import reduce

"""
This module implements representations of slabs and surfaces, as well as
algorithms for generating them.
"""

__author__ = "Richard Tran, Wenhao Sun, Zihan Xu, Shyue Ping Ong"
__copyright__ = "Copyright 2014, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "6/10/14"

from fractions import gcd
import math
import itertools
import logging

import numpy as np
from numpy.linalg import norm
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster

from monty.fractions import lcm

from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord_utils import in_coord_list, find_in_coord_list_pbc, \
    in_coord_list_pbc
from pymatgen.analysis.structure_matcher import StructureMatcher

logger = logging.getLogger(__name__)


class Slab(Structure):
    """
    Subclass of Structure representing a Slab. Implements additional
    attributes pertaining to slabs, but the init method does not
    actually implement any algorithm that creates a slab. This is a
    DUMMY class who's init method only holds information about the
    slab. Also has additional methods that returns other information
    about a slab such as the surface area, normal, and atom adsorption.

    Note that all Slabs have the surface normal oriented in the c-direction.
    This means the lattice vectors a and b are in the surface plane and the c
    vector is out of the surface plane (though not necessary perpendicular to
    the surface.)

    .. attribute:: miller_index

        Miller index of plane parallel to surface.

    .. attribute:: scale_factor

        Final computed scale factor that brings the parent cell to the
        surface cell.

    .. attribute:: shift

        The shift value in Angstrom that indicates how much this
        slab has been shifted.
    """

    def __init__(self, lattice, species, coords, miller_index,
                 oriented_unit_cell, shift, scale_factor,
                 validate_proximity=False, to_unit_cell=False,
                 coords_are_cartesian=False, site_properties=None, energy=None):
        """
        Makes a Slab structure, a structure object with additional information
        and methods pertaining to slabs.

        Args:
            lattice (Lattice/3x3 array): The lattice, either as a
                :class:`pymatgen.core.lattice.Lattice` or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species ([Specie]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
            oriented_unit_cell (Structure): The oriented_unit_cell from which
                this Slab is created (by scaling in the c-direction).
            shift (float): The shift in the c-direction applied to get the
                termination.
            scale_factor (array): scale_factor Final computed scale factor
                that brings the parent cell to the surface cell.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
            energy (float): A value for the energy.
        """
        self.oriented_unit_cell = oriented_unit_cell
        self.miller_index = tuple(miller_index)
        self.shift = shift
        self.scale_factor = scale_factor
        self.energy = energy
        super(Slab, self).__init__(
            lattice, species, coords, validate_proximity=validate_proximity,
            to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties)

    def get_sorted_structure(self, key=None, reverse=False):
        """
        Get a sorted copy of the structure. The parameters have the same
        meaning as in list.sort. By default, sites are sorted by the
        electronegativity of the species. Note that Slab has to override this
        because of the different __init__ args.

        Args:
            key: Specifies a function of one argument that is used to extract
                a comparison key from each list element: key=str.lower. The
                default value is None (compare the elements directly).
            reverse (bool): If set to True, then the list elements are sorted
                as if each comparison were reversed.
        """
        sites = sorted(self, key=key, reverse=reverse)
        s = Structure.from_sites(sites)
        return Slab(s.lattice, s.species_and_occu, s.frac_coords,
                    self.miller_index, self.oriented_unit_cell, self.shift,
                    self.scale_factor, site_properties=s.site_properties)

    @property
    def dipole(self):
        """
        Calculates the dipole of the Slab in the direction of the surface
        normal. Note that the Slab must be oxidation state-decorated for this
        to work properly. Otherwise, the Slab will always have a dipole of 0.
        """
        dipole = np.zeros(3)
        mid_pt = np.sum(self.cart_coords, axis=0) / len(self)
        normal = self.normal
        for site in self:
            charge = sum([getattr(sp, "oxi_state", 0) * amt
                          for sp, amt in site.species_and_occu.items()])
            dipole += charge * np.dot(site.coords - mid_pt, normal) * normal
        return dipole

    def is_polar(self, tol_dipole_per_unit_area=1e-3):
        """
        Checks whether the surface is polar by computing the dipole per unit
        area. Note that the Slab must be oxidation state-decorated for this
        to work properly. Otherwise, the Slab will always be non-polar.

        Args:
            tol_dipole_per_unit_area (float): A tolerance. If the dipole
                magnitude per unit area is less than this value, the Slab is
                considered non-polar. Defaults to 1e-3, which is usually
                pretty good. Normalized dipole per unit area is used as it is
                more reliable than using the total, which tends to be larger for
                slabs with larger surface areas.
        """
        dip_per_unit_area = self.dipole / self.surface_area
        return np.linalg.norm(dip_per_unit_area) > tol_dipole_per_unit_area

    @property
    def normal(self):
        """
        Calculates the surface normal vector of the slab
        """
        normal = np.cross(self.lattice.matrix[0], self.lattice.matrix[1])
        normal /= np.linalg.norm(normal)
        return normal

    @property
    def surface_area(self):
        """
        Calculates the surface area of the slab
        """
        m = self.lattice.matrix
        return np.linalg.norm(np.cross(m[0], m[1]))

    def add_adsorbate_atom(self, indices, specie, distance):
        """
        Gets the structure of single atom adsorption.
        slab structure from the Slab class(in [0, 0, 1])

        Args:
            indices ([int]): Indices of sites on which to put the absorbate.
                Absorbed atom will be displaced relative to the center of
                these sites.
            specie (Specie/Element/str): adsorbed atom species
            distance (float): between centers of the adsorbed atom and the
                given site in Angstroms.
        """
        #Let's do the work in cartesian coords
        center = np.sum([self[i].coords for i in indices], axis=0) / len(
            indices)

        coords = center + self.normal * distance / np.linalg.norm(self.normal)

        self.append(specie, coords, coords_are_cartesian=True)

    def __str__(self):
        comp = self.composition
        outs = [
            "Slab Summary (%s)" % comp.formula,
            "Reduced Formula: %s" % comp.reduced_formula,
            "Miller index: %s" % (self.miller_index, ),
            "Shift: %.4f, Scale Factor: %s" % (self.shift,
                                               self.scale_factor.__str__())]
        to_s = lambda x: "%0.6f" % x
        outs.append("abc   : " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.abc]))
        outs.append("angles: " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.angles]))
        outs.append("Sites ({i})".format(i=len(self)))
        for i, site in enumerate(self):
            outs.append(" ".join([str(i + 1), site.species_string,
                                  " ".join([to_s(j).rjust(12)
                                            for j in site.frac_coords])]))
        return "\n".join(outs)

    def as_dict(self):
        d = super(Slab, self).as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["oriented_unit_cell"] = self.oriented_unit_cell.as_dict()
        d["miller_index"] = self.miller_index
        d["shift"] = self.shift
        d["scale_factor"] = self.scale_factor
        d["energy"] = self.energy
        return d

    @classmethod
    def from_dict(cls, d):
        lattice = Lattice.from_dict(d["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in d["sites"]]
        s = Structure.from_sites(sites)

        return Slab(
            lattice=lattice,
            species=s.species_and_occu, coords=s.frac_coords,
            miller_index=d["miller_index"],
            oriented_unit_cell=Structure.from_dict(d["oriented_unit_cell"]),
            shift=d["shift"], scale_factor=d["scale_factor"],
            site_properties=s.site_properties, energy=d["energy"]
        )


class SlabGenerator(object):

    """
    This class generates different slabs using shift values determined by where a
    unique termination can be found along with other criterias such as where a
    termination doesn't break a polyhedral bond. The shift value then indicates
    where the slab layer will begin and terminate in the slab-vacuum system.

    .. attribute:: oriented_unit_cell

        A unit cell of the parent structure with the miller
        index of plane parallel to surface

    .. attribute:: parent

        Parent structure from which Slab was derived.

    .. attribute:: lll_reduce

        Whether or not the slabs will be orthogonalized

    .. attribute:: center_slab

        Whether or not the slabs will be centered between
        the vacuum layer

    .. attribute:: slab_scale_factor

        Final computed scale factor that brings the parent cell to the
        surface cell.

    .. attribute:: miller_index

        Miller index of plane parallel to surface.

    .. attribute:: min_slab_size

        Minimum size in angstroms of layers containing atoms

    .. attribute:: min_vac_size

        Minimize size in angstroms of layers containing vacuum

    """

    def __init__(self, initial_structure, miller_index, min_slab_size,
                 min_vacuum_size, lll_reduce=False, center_slab=False,
                 primitive=True):
        """
        Calculates the slab scale factor and uses it to generate a unit cell
        of the initial structure that has been oriented by its miller index.
        Also stores the initial information needed later on to generate a slab.

        Args:
            initial_structure (Structure): Initial input structure.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
            min_slab_size (float): In Angstroms
            min_vac_size (float): In Angstroms
            lll_reduce (bool): Whether to perform an LLL reduction on the
                eventual structure.
            center_slab (bool): Whether to center the slab in the cell with
                equal vacuum spacing from the top and bottom.
            primitive (bool): Whether to reduce any generated slabs to a
                primitive cell (this does **not** mean the slab is generated
                from a primitive cell, it simply means that after slab
                generation, we attempt to find shorter lattice vectors,
                which lead to less surface area and smaller cells).
        """
        latt = initial_structure.lattice
        d = abs(reduce(gcd, miller_index))
        miller_index = tuple([int(i / d) for i in miller_index])
        #Calculate the surface normal using the reciprocal lattice vector.
        recp = latt.reciprocal_lattice_crystallographic
        normal = recp.get_cartesian_coords(miller_index)
        normal /= np.linalg.norm(normal)

        slab_scale_factor = []
        non_orth_ind = []
        eye = np.eye(3, dtype=np.int)
        for i, j in enumerate(miller_index):
            if j == 0:
                # Lattice vector is perpendicular to surface normal, i.e.,
                # in plane of surface. We will simply choose this lattice
                # vector as one of the basis vectors.
                slab_scale_factor.append(eye[i])
            else:
                #Calculate projection of lattice vector onto surface normal.
                d = abs(np.dot(normal, latt.matrix[i])) / latt.abc[i]
                non_orth_ind.append((i, d))

        # We want the vector that has maximum magnitude in the
        # direction of the surface normal as the c-direction.
        # Results in a more "orthogonal" unit cell.
        c_index, dist = max(non_orth_ind, key=lambda t: t[1])

        if len(non_orth_ind) > 1:
            lcm_miller = lcm(*[miller_index[i] for i, d in non_orth_ind])
            for (i, di), (j, dj) in itertools.combinations(non_orth_ind, 2):
                l = [0, 0, 0]
                l[i] = -int(round(lcm_miller / miller_index[i]))
                l[j] = int(round(lcm_miller / miller_index[j]))
                slab_scale_factor.append(l)
                if len(slab_scale_factor) == 2:
                    break

        slab_scale_factor.append(eye[c_index])

        slab_scale_factor = np.array(slab_scale_factor)

        # Let's make sure we have a left-handed crystallographic system
        if np.linalg.det(slab_scale_factor) < 0:
            slab_scale_factor *= -1

        single = initial_structure.copy()
        single.make_supercell(slab_scale_factor)

        self.oriented_unit_cell = Structure.from_sites(single,
                                                       to_unit_cell=True)
        self.oriented_unit_cell.to(filename="LCO_111_before.cif")
        if primitive:
            self.oriented_unit_cell = _get_constrained_prim(single)
        self.parent = initial_structure
        self.lll_reduce = lll_reduce
        self.center_slab = center_slab
        self.slab_scale_factor = slab_scale_factor
        self.miller_index = miller_index
        self.min_vac_size = min_vacuum_size
        self.min_slab_size = min_slab_size
        self.primitive = primitive
        self._normal = normal
        a, b, c = self.oriented_unit_cell.lattice.matrix
        self._proj_height = abs(np.dot(normal, c))

    def get_slab(self, shift=0, tol=0.1, energy=None):
        """
        This method takes in shift value for the c lattice direction and
        generates a slab based on the given shift. You should rarely use this
        method. Instead, it is used by other generation algorithms to obtain
        all slabs.

        Arg:
            shift (float): A shift value in Angstrom that determines how much a
                slab should be shifted.
            tol (float): Tolerance to determine primitive cell.
            energy (float): An energy to assign to the slab.

        Returns:
            (Slab) A Slab object with a particular shifted oriented unit cell.
        """

        h = self._proj_height
        nlayers_slab = int(math.ceil(self.min_slab_size / h))
        nlayers_vac = int(math.ceil(self.min_vac_size / h))
        nlayers = nlayers_slab + nlayers_vac

        species = self.oriented_unit_cell.species_and_occu
        props = self.oriented_unit_cell.site_properties
        props = {k: v * nlayers_slab for k, v in props.items()}
        frac_coords = self.oriented_unit_cell.frac_coords
        frac_coords = np.array(frac_coords) +\
                      np.array([0, 0, -shift])[None, :]
        frac_coords = frac_coords - np.floor(frac_coords)
        a, b, c = self.oriented_unit_cell.lattice.matrix
        new_lattice = [a, b, nlayers * c]
        frac_coords[:, 2] = frac_coords[:, 2] / nlayers
        all_coords = []
        for i in range(nlayers_slab):
            fcoords = frac_coords.copy()
            fcoords[:, 2] += i / nlayers
            all_coords.extend(fcoords)

        slab = Structure(new_lattice, species * nlayers_slab, all_coords,
                         site_properties=props)

        scale_factor = self.slab_scale_factor
        # Whether or not to orthogonalize the structure
        if self.lll_reduce:
            lll_slab = slab.copy(sanitize=True)
            mapping = lll_slab.lattice.find_mapping(slab.lattice)
            scale_factor = np.dot(mapping[2], scale_factor)
            slab = lll_slab

        # Whether or not to center the slab layer around the vacuum
        if self.center_slab:
            avg_c = np.average([c[2] for c in slab.frac_coords])
            slab.translate_sites(list(range(len(slab))), [0, 0, 0.5 - avg_c])

        if self.primitive:
            prim = slab.get_primitive_structure(tolerance=tol)
            if energy is not None:
                energy = prim.volume / slab.volume * energy
            slab = prim

        return Slab(slab.lattice, slab.species_and_occu,
                    slab.frac_coords, self.miller_index,
                    self.oriented_unit_cell, shift,
                    scale_factor, site_properties=slab.site_properties,
                    energy=energy)

    def _calculate_possible_shifts(self, tol=0.1):
        frac_coords = self.oriented_unit_cell.frac_coords
        n = len(frac_coords)

        # We cluster the sites according to the c coordinates. But we need to
        # take into account PBC. Let's compute a fractional c-coordinate
        # distance matrix that accounts for PBC.
        dist_matrix = np.zeros((n, n))
        h = self._proj_height
        # Projection of c lattice vector in
        # direction of surface normal.
        for i, j in itertools.combinations(list(range(n)), 2):
            if i != j:
                cdist = frac_coords[i][2] - frac_coords[j][2]
                cdist = abs(cdist - round(cdist)) * h
                dist_matrix[i, j] = cdist
                dist_matrix[j, i] = cdist

        condensed_m = squareform(dist_matrix)
        z = linkage(condensed_m)
        clusters = fcluster(z, tol, criterion="distance")

        #Generate dict of cluster# to c val - doesn't matter what the c is.
        c_loc = {c: frac_coords[i][2] for i, c in enumerate(clusters)}

        #Put all c into the unit cell.
        possible_c = [c - math.floor(c) for c in sorted(c_loc.values())]

        # Calculate the shifts
        nshifts = len(possible_c)
        shifts = []
        for i in range(nshifts):
            if i == nshifts - 1:
                # There is an additional shift between the first and last c
                # coordinate. But this needs special handling because of PBC.
                shift = (possible_c[0] + 1 + possible_c[i]) * 0.5
                if shift > 1:
                    shift -= 1
            else:
                shift = (possible_c[i] + possible_c[i + 1]) * 0.5
            shifts.append(shift - math.floor(shift))
        shifts = sorted(shifts)
        return shifts

    def get_slabs(self, bonds=None, tol=0.1, max_broken_bonds=0):
        """
        This method returns a list of slabs that are generated using the list of
        shift values from the method, _calculate_possible_shifts(). Before the
        shifts are used to create the slabs however, if the user decides to take
        into account whether or not a termination will break any polyhedral
        structure (bonds != None), this method will filter out any shift values
        that do so.

        Args:
            bonds ({(specie1, specie2): max_bond_dist}: bonds are
                specified as a dict of tuples: float of specie1, specie2
                and the max bonding distance. For example, PO4 groups may be
                defined as {("P", "O"): 3}.
            tol (float): Threshold parameter in fcluster in order to check
                if two atoms are lying on the same plane. Default thresh set
                to 0.1 Angstrom in the direction of the surface normal.
            max_broken_bonds (int): Maximum number of allowable broken bonds
                for the slab. Use this to limit # of slabs (some structures
                may have a lot of slabs). Defaults to zero, which means no
                defined bonds must be broken.

        Returns:
            ([Slab]) List of all possible terminations of a particular surface.
            Slabs are sorted by the # of bonds broken.
        """
        def get_c_ranges(site1, site2, bond_dist):
            lattice = site1.lattice
            f1 = site1.frac_coords
            c_ranges = []
            for dist, image in lattice.get_all_distance_and_image(
                    f1, site2.frac_coords):
                if dist < bond_dist:
                    f2 = site2.frac_coords + image
                    # Checks if the distance between the two species
                    # is less then the user input bond distance
                    min_c = f1[2]
                    max_c = f2[2]
                    c_range = sorted([min_c, max_c])
                    if c_range[1] > 1:
                        # Takes care of PBC when c coordinate of site
                        # goes beyond the upper boundary of the cell
                        c_ranges.append((c_range[0], 1))
                        c_ranges.append((0, c_range[1] - 1))
                    elif c_range[0] < 0:
                        # Takes care of PBC when c coordinate of site
                        # is below the lower boundary of the unit cell
                        c_ranges.append((0, c_range[1]))
                        c_ranges.append((c_range[0] + 1, 1))
                    else:
                        c_ranges.append(c_range)
            return c_ranges

        bond_c_ranges = []
        if bonds is not None:
            #Convert to species first
            bonds = {(get_el_sp(s1), get_el_sp(s2)): dist for (s1, s2), dist in
                     bonds.items()}
            for s1, s2 in itertools.combinations(self.oriented_unit_cell, 2):
                # Iterates through every possible pair of species in the
                # oriented unit cell
                all_sp = set(s1.species_and_occu.keys())
                all_sp.update(s2.species_and_occu.keys())
                for species, bond_dist in bonds.items():
                    # Checks if elements in species is in all_sp
                    if all_sp.issuperset(species):
                        bond_c_ranges.extend(get_c_ranges(s1, s2, bond_dist))

        slabs = []
        for shift in self._calculate_possible_shifts(tol=tol):
            bonds_broken = 0
            for r in bond_c_ranges:
                if r[0] <= shift <= r[1]:
                    bonds_broken += 1
            if bonds_broken <= max_broken_bonds:
                # For now, set the energy to be equal to no. of broken bonds
                # per unit cell.
                slab = self.get_slab(shift, tol=tol, energy=bonds_broken)
                slabs.append(slab)

        # Further filters out any surfaces made that might be the same
        m = StructureMatcher(ltol=tol, stol=tol, primitive_cell=False,
                             scale=False)
        slabs = [g[0] for g in m.group_structures(slabs)]
        return sorted(slabs, key=lambda s: s.energy)


def get_symmetrically_distinct_miller_indices(structure, max_index):
    """
    Returns all symmetrically distinct indices below a certain max-index for
    a given structure. Analysis is based on the symmetry of the reciprocal
    lattice of the structure.

    Args:
        structure (Structure): input structure.
        max_index (int): The maximum index. For example, a max_index of 1
            means that (100), (110), and (111) are returned for the cubic
            structure. All other indices are equivalent to one of these.
    """

    recp_lattice = structure.lattice.reciprocal_lattice_crystallographic
    # Need to make sure recp lattice is big enough, otherwise symmetry
    # determination will fail. We set the overall volume to 1.
    recp_lattice = recp_lattice.scale(1)
    recp = Structure(recp_lattice, ["H"], [[0, 0, 0]])

    # Creates a function that uses the symmetry operations in the
    # structure to find Miller indices that might give repetitive slabs
    analyzer = SpacegroupAnalyzer(recp, symprec=0.001)
    symm_ops = analyzer.get_symmetry_operations()
    unique_millers = []

    def is_already_analyzed(miller_index):
        for op in symm_ops:
            if in_coord_list(unique_millers, op.operate(miller_index)):
                return True
        return False

    r = list(range(-max_index, max_index + 1))
    r.reverse()
    for miller in itertools.product(r, r, r):
        if any([i != 0 for i in miller]):
            d = abs(reduce(gcd, miller))
            miller = tuple([int(i / d) for i in miller])
            if not is_already_analyzed(miller):
                unique_millers.append(miller)
    return unique_millers


def generate_all_slabs(structure, max_index, min_slab_size, min_vacuum_size,
                       bonds=None, tol=0.1, max_broken_bonds=0,
                       lll_reduce=False,
                       center_slab=False, primitive=True):
    """
    A function that finds all different slabs up to a certain miller index.
    Slabs oriented under certain Miller indices that are equivalent to other
    slabs in other Miller indices are filtered out using symmetry operations
    to get rid of any repetitive slabs. For example, under symmetry operations,
    CsCl has equivalent slabs in the (0,0,1), (0,1,0), and (1,0,0) direction.

    Args:
        structure (Structure): Initial input structure.
        max_index (int): The maximum Miller index to go up to.
        symprec (float): Tolerance for symmetry finding. Defaults to 1e-3,
            which is fairly strict and works well for properly refined
            structures with atoms in the proper symmetry coordinates. For
            structures with slight deviations from their proper atomic
            positions (e.g., structures relaxed with electronic structure
            codes), a looser tolerance of 0.1 (the value used in Materials
            Project) is often needed.
        lll_reduce (bool): Whether to perform an LLL reduction on the
            eventual structure.
        center_slab (bool): Whether to center the slab in the cell with
            equal vacuum spacing from the top and bottom.
        primitive (bool): Whether to reduce any generated slabs to a
            primitive cell (this does **not** mean the slab is generated
            from a primitive cell, it simply means that after slab
            generation, we attempt to find shorter lattice vectors,
            which lead to less surface area and smaller cells).
    """
    all_slabs = []
    for miller in get_symmetrically_distinct_miller_indices(structure,
                                                            max_index):
        gen = SlabGenerator(structure, miller, min_slab_size,
                            min_vacuum_size, lll_reduce=lll_reduce,
                            center_slab=center_slab, primitive=primitive)
        slabs = gen.get_slabs(bonds=bonds, tol=tol,
                              max_broken_bonds=max_broken_bonds)
        if len(slabs) > 0:
            logger.debug("%s has %d slabs... " % (miller, len(slabs)))
            all_slabs.extend(slabs)

    return all_slabs


def _is_valid_vec(structure, vec, atol):
    mapping = []
    lattice = structure.lattice
    fcoords = structure.frac_coords
    for i, s in enumerate(structure):
        f = lattice.get_fractional_coords(s.coords + vec)
        indices = find_in_coord_list_pbc(fcoords, f, atol=atol)
        if len(indices) != 1 or structure[indices[0]].species_and_occu != s.species_and_occu:
            return False
        else:
            mapping.append(indices[0])

    return len(mapping) == len(set(mapping))


def _get_constrained_prim(structure, atol=1e-2):
    sorted_structure = structure.get_sorted_structure()

    sorted_structure = Structure.from_sites(sorted_structure, to_unit_cell=True)
    groups = [
        list(g) for k, g in itertools.groupby(
            sorted_structure, key=lambda site: site.species_and_occu)]
    min_sites = min(groups, key=lambda g: len(g))
    vol = structure.lattice.volume
    a, b, c = structure.lattice.matrix
    for s1, s2 in itertools.combinations(min_sites, 2):
        vec = s1.coords - s2.coords
        print s1.frac_coords
        print s2.frac_coords
        print vec
        if _is_valid_vec(structure, vec, atol):
            d = abs(np.linalg.det([a, b, vec]))

            if d < atol: #In a-b plane
                a_norm = norm(np.cross(vec, a)) / norm(vec) / norm(a)
                b_norm = norm(np.cross(vec, b)) / norm(vec) / norm(b)

                if a_norm < b_norm:
                    new_lattice = Lattice([vec, b, c])
                else:
                    new_lattice = Lattice([a, vec, c])
            else:
                new_lattice = Lattice([a, b, vec])
            print new_lattice.lengths_and_angles
            print structure.lattice.lengths_and_angles
            if new_lattice.volume < 0.1 * vol:
                continue
            print "here"
            sp = []
            fcoords = []
            for site in structure:
                f = new_lattice.get_fractional_coords(site.coords)
                if not in_coord_list_pbc(fcoords, f, atol):
                    sp.append(site.species_and_occu)
                    fcoords.append(f)
            prim = Structure(new_lattice, sp, fcoords, to_unit_cell=True)
            return _get_constrained_prim(prim)
    return structure


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    import os
    cwd = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(cwd, "..", "..", "test_files", "surface_tests",
                    "icsd_LiCoO2.cif")
    s = Structure.from_file(path, primitive=False)
    print s
    gen = SlabGenerator(s, miller_index=(1, 1, 1), min_slab_size=4,
                        min_vacuum_size=4)
    gen.oriented_unit_cell.to(filename="LCO111.cif")
    #print s
    #print _get_constrained_prim(s)
    # s = Structure.from_file("../../test_files/LiFePO4").get_structure("LiFePO4")
    # self.tei = Structure.from_file(get_path("icsd_TeI.cif"),
    #                                primitive=False)
    # self.LiCoO2 = Structure.from_file(get_path("icsd_LiCoO2.cif"),
    #                                   primitive=False)
    #
    # self.p1 = Structure(Lattice.from_parameters(3, 4, 5, 31, 43, 50),
    #                     ["H", "He"], [[0, 0, 0], [0.1, 0.2, 0.3]])