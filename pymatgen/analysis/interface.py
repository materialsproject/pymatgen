# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to store, generate, and manipulate material interfaces.
"""

import warnings
from itertools import product

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.substrate_analyzer import SubstrateAnalyzer, reduce_vectors
from pymatgen.core.operations import SymmOp
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Eric Sivonxay, Shyam Dwaraknath, and Kyle Bystrom"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Kyle Bystrom"
__email__ = "kylebystrom@gmail.com"
__date__ = "5/29/2019"
__status__ = "Prototype"


class Interface(Structure):
    """
    This class stores data for defining an interface between two structures.
    It is a subclass of pymatgen.core.structure.Structure.
    """

    def __init__(
        self,
        lattice,
        species,
        coords,
        sub_plane,
        film_plane,
        sub_init_cell,
        film_init_cell,
        modified_sub_structure,
        modified_film_structure,
        strained_sub_structure,
        strained_film_structure,
        validate_proximity=False,
        coords_are_cartesian=False,
        init_inplane_shift=None,
        charge=None,
        site_properties=None,
        to_unit_cell=False,
    ):
        """
        Makes an interface structure, a Structure object with additional
        information and methods pertaining to interfaces.

        Args:
            lattice (Lattice/3x3 array): The lattice, either as a
                :class:`pymatgen.core.lattice.Lattice` or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species ([Species]): Sequence of species on each site. Can take in
                flexible input, including:
                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Species objects.
                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            sub_plane (list): Substrate plane in the form of a list of integers
                (based on the sub_init_cell), e.g.: [1, 2, 3].
            film_plane (list): Film plane in the form of a list of integers
                (based on the film_init_cell), e.g. [1, 2, 3].
            sub_init_cell (Structure): initial bulk substrate structure
            film_init_cell (Structure): initial bulk film structure
            site_properties (dict): Properties associated with the sites as a
                dict of sequences. The sequences have to be the same length as
                the atomic species and fractional_coords. For an interface, you should
                have the 'interface_label' properties to classify the sites as
                'substrate' and 'film'.
            modified_sub_structure (Slab): substrate supercell slab.
            modified_film_structure (Slab): film supercell slab.
            strained_sub_structure (Slab): strained substrate supercell slab
            strained_film_structure (Slab): strained film supercell slab
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in cartesian coordinates. Defaults to False.
            init_inplane_shift (length-2 list of float, in Cartesian coordinates):
                The initial shift of the film relative to the substrate
                in the plane of the interface.
            charge (float, optional): overal charge of the structure
        """

        super().__init__(
            lattice,
            species,
            coords,
            validate_proximity=validate_proximity,
            to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties,
            charge=charge,
        )

        self.modified_sub_structure = modified_sub_structure
        self.modified_film_structure = modified_film_structure
        self.strained_sub_structure = strained_sub_structure
        self.strained_film_structure = strained_film_structure
        self.sub_plane = sub_plane
        self.film_plane = film_plane
        self.sub_init_cell = sub_init_cell
        self.film_init_cell = film_init_cell

        z_shift = np.min(self.film.cart_coords[:, 2]) - np.max(self.substrate.cart_coords[:, 2])

        if init_inplane_shift is None:
            init_inplane_shift = np.array([0.0, 0.0])

        self._offset_vector = np.append(init_inplane_shift, [z_shift])

    def shift_film_along_surface_lattice(self, da, db):
        """
        Given two floats da and db, adjust the shift vector
        by da * (first lattice vector) + db * (second lattice vector).
        This shift is in the plane of the interface.
        I.e. da and db are fractional coordinates.

        Args:
            da (float): shift in the first lattice vector
            db (float): shift in the second lattice vector
        """
        self.shift_film(da * self.lattice.matrix[0] + db * self.lattice.matrix[1])

    def change_z_shift(self, dz):
        """
        Adjust the spacing between the substrate and film layers by dz Angstroms

        Args:
            dz (float): shift perpendicular to the plane (in Angstroms)
        """
        self.shift_film(np.array([0.0, 0.0, dz]))

    def shift_film(self, delta):
        """
        Shift the film's position relative to the substrate.

        Args:
            delta (length-3 list of float or numpy array): Cartesian coordinate
                vector by which to shift the film. After this operation
                self.offset_vector -> self.offset_vector + delta.
        """
        if self.offset_vector[2] + delta[2] < 0 or delta[2] > self.vacuum_thickness:
            raise ValueError("The shift {} will collide the film and substrate.".format(delta))
        self._offset_vector += np.array(delta)
        self.translate_sites(self.get_film_indices(), delta, frac_coords=False, to_unit_cell=True)

    @property
    def offset_vector(self):
        """
        Displacement of the origin of the film structure relative to that
        of the substrate structure in Cartesian coordinates.
        """
        return self._offset_vector.copy()

    @offset_vector.setter
    def offset_vector(self, offset_vector):
        delta = offset_vector - self._offset_vector
        self.shift_film(delta)

    @property
    def ab_shift(self):
        """
        The 2D component of offset_vector along the interface plane
        in fractional coordinates. I.e. if ab_shift = [a, b], the
        Cartesian coordinate shift in the interface plane
        is a * (first lattice vector) + b * (second lattice vector).
        """
        return np.dot(self.offset_vector, np.linalg.inv(self.lattice.matrix))[:2]

    @ab_shift.setter
    def ab_shift(self, ab_shift):
        delta = ab_shift - self.ab_shift
        self.shift_film_along_surface_lattice(delta[0], delta[1])

    @property
    def z_shift(self):
        """
        The 1D component of offset_vector along the interface plane
        in fractional coordinates. I.e. if z_shift = z, the distance
        between the substrate and film planes is z.
        """
        return self.offset_vector[2]

    @z_shift.setter
    def z_shift(self, z_shift):
        delta = z_shift - self.z_shift
        self.change_z_shift(delta)

    @property
    def vacuum_thickness(self):
        """
        Vacuum buffer above the film.
        """
        return np.min(self.substrate.cart_coords[:, 2]) + self.lattice.c - np.max(self.film.cart_coords[:, 2])

    @property
    def substrate_sites(self):
        """
        Return the substrate sites of the interface.
        """
        sub_sites = []
        for i, tag in enumerate(self.site_properties["interface_label"]):
            if "substrate" in tag:
                sub_sites.append(self.sites[i])
        return sub_sites

    @property
    def substrate(self):
        """
        Return the substrate (Structure) of the interface.
        """
        return Structure.from_sites(self.substrate_sites)

    def get_film_indices(self):
        """
        Retrieve the indices of the film sites
        """
        film_sites = []
        for i, tag in enumerate(self.site_properties["interface_label"]):
            if "film" in tag:
                film_sites.append(i)
        return film_sites

    @property
    def film_sites(self):
        """
        Return the film sites of the interface.
        """
        film_sites = []
        for i, tag in enumerate(self.site_properties["interface_label"]):
            if "film" in tag:
                film_sites.append(self.sites[i])
        return film_sites

    @property
    def film(self):
        """
        Return the film (Structure) of the interface.
        """
        return Structure.from_sites(self.film_sites)

    def copy(self, site_properties=None):
        """
        Convenience method to get a copy of the structure, with options to add
        site properties.

        Returns:
            A copy of the Interface.
        """
        props = self.site_properties
        if site_properties:
            props.update(site_properties)
        return Interface(
            self.lattice,
            self.species_and_occu,
            self.frac_coords,
            self.sub_plane,
            self.film_plane,
            self.sub_init_cell,
            self.film_init_cell,
            self.modified_sub_structure,
            self.modified_film_structure,
            self.strained_sub_structure,
            self.strained_film_structure,
            validate_proximity=False,
            coords_are_cartesian=False,
            init_inplane_shift=self.offset_vector[:2],
            charge=self.charge,
            site_properties=self.site_properties,
        )

    def get_sorted_structure(self, key=None, reverse=False):
        """
        Get a sorted copy of the structure. The parameters have the same
        meaning as in list.sort. By default, sites are sorted by the
        electronegativity of the species.

        Args:
            key: Specifies a function of one argument that is used to extract
                a comparison key from each list element: key=str.lower. The
                default value is None (compare the elements directly).
            reverse (bool): If set to True, then the list elements are sorted
                as if each comparison were reversed.
        """
        struct_copy = self.copy()
        struct_copy.sort(key=key, reverse=reverse)
        return struct_copy

    def as_dict(self):
        """
        :return: MSONable dict
        """
        d = super().as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["sub_plane"] = self.sub_plane
        d["film_plane"] = self.film_plane
        d["sub_init_cell"] = self.sub_init_cell
        d["film_init_cell"] = self.film_init_cell
        d["modified_sub_structure"] = self.modified_sub_structure
        d["modified_film_structure"] = self.modified_film_structure
        d["strained_sub_structure"] = self.strained_sub_structure
        d["strained_film_structure"] = self.strained_film_structure
        d["init_inplane_shift"] = self.offset_vector[0:2]
        return d

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: Interface
        """
        lattice = Lattice.from_dict(d["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in d["sites"]]
        s = Structure.from_sites(sites)

        return Interface(
            lattice=lattice,
            species=s.species_and_occu,
            coords=s.frac_coords,
            sub_plane=d["sub_plane"],
            film_plane=d["film_plane"],
            sub_init_cell=d["sub_init_cell"],
            film_init_cell=d["film_init_cell"],
            modified_sub_structure=d["modified_sub_structure"],
            modified_film_structure=d["modified_film_structure"],
            strained_sub_structure=d["strained_sub_structure"],
            strained_film_structure=d["strained_film_structure"],
            site_properties=s.site_properties,
            init_inplane_shift=d["init_inplane_shift"],
        )


class InterfaceBuilder:
    """
    This class constructs the epitaxially matched interfaces between two crystalline slabs
    """

    def __init__(self, substrate_structure, film_structure):
        """
        Args:
            substrate_structure (Structure): structure of substrate
            film_structure (Structure): structure of film
        """

        # Bulk structures
        self.original_substrate_structure = substrate_structure
        self.original_film_structure = film_structure

        self.matches = []

        self.match_index = None

        # SlabGenerator objects for the substrate and film
        self.sub_sg = None
        self.substrate_layers = None
        self.film_sg = None
        self.film_layers = None

        # Structures with no vacuum
        self.substrate_structures = []
        self.film_structures = []

        # "slab" structure (with no vacuum) oriented with a direction along x-axis and ab plane normal aligned with
        # z-axis
        self.oriented_substrate = None
        self.oriented_film = None

        # Strained structures with no vacuum
        self.strained_substrate = None
        self.strained_film = None

        # Substrate with transformation/matches applied
        self.modified_substrate_structures = []
        self.modified_film_structures = []

        # Non-stoichiometric slabs with symmetric surfaces, as generated by pymatgen. Please check, this is highly
        # unreliable from tests.
        self.sym_modified_substrate_structures = []
        self.sym_modified_film_structures = []

        # Interface structures
        self.interfaces = []
        self.interface_labels = []

    def get_summary_dict(self):
        """
        Return dictionary with information about the InterfaceBuilder,
        with currently generated structures included.
        """

        d = {"match": self.matches[0]}
        d["substrate_layers"] = self.substrate_layers
        d["film_layers"] = self.film_layers

        d["bulk_substrate"] = self.original_substrate_structure
        d["bulk_film"] = self.original_film_structure

        d["strained_substrate"] = self.strained_substrate
        d["strained_film"] = self.strained_film

        d["slab_substrates"] = self.modified_substrate_structures
        d["slab_films"] = self.modified_film_structures

        d["interfaces"] = self.interfaces
        d["interface_labels"] = self.interface_labels

        return d

    def write_all_structures(self):
        """
        Write all of the structures relevant for
        the interface calculation to VASP POSCAR files.
        """

        _poscar = Poscar(self.original_substrate_structure)
        _poscar.write_file("bulk_substrate_POSCAR")

        _poscar = Poscar(self.original_film_structure)
        _poscar.write_file("bulk_film_POSCAR")

        _poscar = Poscar(self.strained_substrate)
        _poscar.write_file("strained_substrate_POSCAR")

        _poscar = Poscar(self.strained_film)
        _poscar.write_file("strained_film_POSCAR")

        for i, interface in enumerate(self.modified_substrate_structures):
            _poscar = Poscar(interface)
            _poscar.write_file("slab_substrate_%d_POSCAR" % i)

        for i, interface in enumerate(self.modified_film_structures):
            _poscar = Poscar(interface)
            _poscar.write_file("slab_film_%d_POSCAR" % i)

        for i, interface in enumerate(self.film_structures):
            _poscar = Poscar(interface)
            _poscar.write_file("slab_unit_film_%d_POSCAR" % i)

        for label, interface in zip(self.interface_labels, self.interfaces):
            _poscar = Poscar(interface)
            _poscar.write_file("interface_%s_POSCAR" % label.replace("/", "-"))

    def generate_interfaces(
        self, film_millers=None, substrate_millers=None, film_layers=3, substrate_layers=3, **kwargs
    ):
        """
        Generate a list of Interface (Structure) objects and store them to self.interfaces.

        Args:
            film_millers (list of [int]): list of film surfaces
            substrate_millers (list of [int]): list of substrate surfaces
            film_layers (int): number of layers of film to include in Interface structures.
            substrate_layers (int): number of layers of substrate to include in Interface structures.
        """
        self.get_oriented_slabs(
            lowest=True,
            film_millers=film_millers,
            substrate_millers=substrate_millers,
            film_layers=film_layers,
            substrate_layers=substrate_layers,
        )

        self.combine_slabs(**kwargs)

    def get_oriented_slabs(self, film_layers=3, substrate_layers=3, match_index=0, **kwargs):
        """
        Get a list of oriented slabs for constructing interfaces and put them
        in self.film_structures, self.substrate_structures, self.modified_film_structures,
        and self.modified_substrate_structures.
        Currently only uses first match (lowest SA) in the list of matches

        Args:
            film_layers (int): number of layers of film to include in Interface structures.
            substrate_layers (int): number of layers of substrate to include in Interface structures.
            match_index (int): ZSL match from which to construct slabs.
        """
        self.match_index = match_index
        self.substrate_layers = substrate_layers
        self.film_layers = film_layers

        if "zslgen" in kwargs.keys():
            sa = SubstrateAnalyzer(zslgen=kwargs.get("zslgen"))
            del kwargs["zslgen"]
        else:
            sa = SubstrateAnalyzer()

        # Generate all possible interface matches
        self.matches = list(sa.calculate(self.original_film_structure, self.original_substrate_structure, **kwargs))
        match = self.matches[match_index]

        # Generate substrate slab and align x axis to (100) and slab normal to (001)
        # Get no-vacuum structure for strained bulk calculation
        self.sub_sg = SlabGenerator(
            self.original_substrate_structure,
            match["sub_miller"],
            substrate_layers,
            0,
            in_unit_planes=True,
            reorient_lattice=False,
            primitive=False,
        )
        no_vac_sub_slab = self.sub_sg.get_slab()
        no_vac_sub_slab = get_shear_reduced_slab(no_vac_sub_slab)
        self.oriented_substrate = align_x(no_vac_sub_slab)
        self.oriented_substrate.sort()

        # Get slab with vacuum
        self.sub_sg = SlabGenerator(
            self.original_substrate_structure,
            match["sub_miller"],
            substrate_layers,
            1,
            in_unit_planes=True,
            reorient_lattice=False,
            primitive=False,
        )
        sub_slabs = self.sub_sg.get_slabs()
        for i, sub_slab in enumerate(sub_slabs):
            sub_slab = get_shear_reduced_slab(sub_slab)
            sub_slab = align_x(sub_slab)
            sub_slab.sort()
            sub_slabs[i] = sub_slab

        self.substrate_structures = sub_slabs

        # Generate film slab and align x axis to (100) and slab normal to (001)
        # Get no-vacuum structure for strained bulk calculation
        self.film_sg = SlabGenerator(
            self.original_film_structure,
            match["film_miller"],
            film_layers,
            0,
            in_unit_planes=True,
            reorient_lattice=False,
            primitive=False,
        )
        no_vac_film_slab = self.film_sg.get_slab()
        no_vac_film_slab = get_shear_reduced_slab(no_vac_film_slab)
        self.oriented_film = align_x(no_vac_film_slab)
        self.oriented_film.sort()

        # Get slab with vacuum
        self.film_sg = SlabGenerator(
            self.original_film_structure,
            match["film_miller"],
            film_layers,
            1,
            in_unit_planes=True,
            reorient_lattice=False,
            primitive=False,
        )
        film_slabs = self.film_sg.get_slabs()
        for i, film_slab in enumerate(film_slabs):
            film_slab = get_shear_reduced_slab(film_slab)
            film_slab = align_x(film_slab)
            film_slab.sort()
            film_slabs[i] = film_slab

        self.film_structures = film_slabs

        # Apply transformation to produce matched area and a & b vectors
        self.apply_transformations(match)

        # Get non-stoichioimetric substrate slabs
        sym_sub_slabs = []
        for sub_slab in self.modified_substrate_structures:
            sym_sub_slab = self.sub_sg.nonstoichiometric_symmetrized_slab(sub_slab)
            for slab in sym_sub_slab:
                if not slab == sub_slab:
                    sym_sub_slabs.append(slab)

        self.sym_modified_substrate_structures = sym_sub_slabs

        # Get non-stoichioimetric film slabs
        sym_film_slabs = []
        for film_slab in self.modified_film_structures:
            sym_film_slab = self.film_sg.nonstoichiometric_symmetrized_slab(film_slab)
            for slab in sym_film_slab:
                if not slab == film_slab:
                    sym_film_slabs.append(slab)

        self.sym_modified_film_structures = sym_film_slabs

        # Strained film structures (No Vacuum)
        self.strained_substrate, self.strained_film = strain_slabs(self.oriented_substrate, self.oriented_film)

    @staticmethod
    def apply_transformation(structure, matrix):
        """
        Make a supercell of structure using matrix

        Args:
            structure (Slab): Slab to make supercell of
            matrix (3x3 np.ndarray): supercell matrix

        Returns:
            (Slab) The supercell of structure
        """
        modified_substrate_structure = structure.copy()
        # Apply scaling
        modified_substrate_structure.make_supercell(matrix)

        # Reduce vectors
        new_lattice = modified_substrate_structure.lattice.matrix.copy()
        new_lattice[:2, :] = reduce_vectors(*modified_substrate_structure.lattice.matrix[:2, :])
        modified_substrate_structure = Slab(
            lattice=Lattice(new_lattice),
            species=modified_substrate_structure.species,
            coords=modified_substrate_structure.cart_coords,
            miller_index=modified_substrate_structure.miller_index,
            oriented_unit_cell=modified_substrate_structure.oriented_unit_cell,
            shift=modified_substrate_structure.shift,
            scale_factor=modified_substrate_structure.scale_factor,
            coords_are_cartesian=True,
            energy=modified_substrate_structure.energy,
            reorient_lattice=modified_substrate_structure.reorient_lattice,
            to_unit_cell=True,
        )

        return modified_substrate_structure

    def apply_transformations(self, match):
        """
        Using ZSL match, transform all of the film_structures by the ZSL
        supercell transformation.

        Args:
            match (dict): ZSL match returned by ZSLGenerator.__call__
        """
        film_transformation = match["film_transformation"]
        sub_transformation = match["substrate_transformation"]

        modified_substrate_structures = [struct.copy() for struct in self.substrate_structures]
        modified_film_structures = [struct.copy() for struct in self.film_structures]

        # Match angles in lattices with 𝛾=θ° and 𝛾=(180-θ)°
        if np.isclose(
            180 - modified_film_structures[0].lattice.gamma,
            modified_substrate_structures[0].lattice.gamma,
            atol=3,
        ):
            reflection = SymmOp.from_rotation_and_translation(((-1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 1))
            for modified_film_structure in modified_film_structures:
                modified_film_structure.apply_operation(reflection, fractional=True)
            self.oriented_film.apply_operation(reflection, fractional=True)

        sub_scaling = np.diag(np.diag(sub_transformation))

        # Turn into 3x3 Arrays
        sub_scaling = np.diag(np.append(np.diag(sub_scaling), 1))
        temp_matrix = np.diag([1, 1, 1])
        temp_matrix[:2, :2] = sub_transformation

        for modified_substrate_structure in modified_substrate_structures:
            modified_substrate_structure = self.apply_transformation(modified_substrate_structure, temp_matrix)
            self.modified_substrate_structures.append(modified_substrate_structure)

        self.oriented_substrate = self.apply_transformation(self.oriented_substrate, temp_matrix)

        film_scaling = np.diag(np.diag(film_transformation))

        # Turn into 3x3 Arrays
        film_scaling = np.diag(np.append(np.diag(film_scaling), 1))
        temp_matrix = np.diag([1, 1, 1])
        temp_matrix[:2, :2] = film_transformation

        for modified_film_structure in modified_film_structures:
            modified_film_structure = self.apply_transformation(modified_film_structure, temp_matrix)
            self.modified_film_structures.append(modified_film_structure)

        self.oriented_film = self.apply_transformation(self.oriented_film, temp_matrix)

    def combine_slabs(self):
        """
        Combine the slabs generated by get_oriented_slabs into interfaces
        """

        all_substrate_variants = []
        sub_labels = []
        for i, slab in enumerate(self.modified_substrate_structures):
            all_substrate_variants.append(slab)
            sub_labels.append(str(i))
            sg = SpacegroupAnalyzer(slab, symprec=1e-3)
            if not sg.is_laue():
                mirrored_slab = slab.copy()
                reflection_z = SymmOp.from_rotation_and_translation(((1, 0, 0), (0, 1, 0), (0, 0, -1)), (0, 0, 0))
                mirrored_slab.apply_operation(reflection_z, fractional=True)
                translation = [0, 0, -min(mirrored_slab.frac_coords[:, 2])]
                mirrored_slab.translate_sites(range(mirrored_slab.num_sites), translation)
                all_substrate_variants.append(mirrored_slab)
                sub_labels.append("%dm" % i)

        all_film_variants = []
        film_labels = []
        for i, slab in enumerate(self.modified_film_structures):
            all_film_variants.append(slab)
            film_labels.append(str(i))
            sg = SpacegroupAnalyzer(slab, symprec=1e-3)
            if not sg.is_laue():
                mirrored_slab = slab.copy()
                reflection_z = SymmOp.from_rotation_and_translation(((1, 0, 0), (0, 1, 0), (0, 0, -1)), (0, 0, 0))
                mirrored_slab.apply_operation(reflection_z, fractional=True)
                translation = [0, 0, -min(mirrored_slab.frac_coords[:, 2])]
                mirrored_slab.translate_sites(range(mirrored_slab.num_sites), translation)
                all_film_variants.append(mirrored_slab)
                film_labels.append("%dm" % i)

        # substrate first index, film second index
        self.interfaces = []
        self.interface_labels = []
        # self.interfaces = [[None for j in range(len(all_film_variants))] for i in range(len(all_substrate_variants))]
        for i, substrate in enumerate(all_substrate_variants):
            for j, film in enumerate(all_film_variants):
                self.interfaces.append(self.make_interface(substrate, film))
                self.interface_labels.append("%s/%s" % (film_labels[j], sub_labels[i]))

    def make_interface(self, slab_substrate, slab_film, offset=None):
        """
        Strain a film to fit a substrate and generate an interface.

        Args:
            slab_substrate (Slab): substrate structure supercell
            slab_film (Slab): film structure supercell
            offset ([int]): separation vector of film and substrate
        """

        # Check if lattices are equal. If not, strain them to match
        # NOTE: CHANGED THIS TO MAKE COPY OF SUBSTRATE/FILM, self.modified_film_structures NO LONGER STRAINED
        unstrained_slab_substrate = slab_substrate.copy()
        slab_substrate = slab_substrate.copy()
        unstrained_slab_film = slab_film.copy()
        slab_film = slab_film.copy()
        latt_1 = slab_substrate.lattice.matrix.copy()
        latt_1[2, :] = [0, 0, 1]
        latt_2 = slab_film.lattice.matrix.copy()
        latt_2[2, :] = [0, 0, 1]
        if Lattice(latt_1) != Lattice(latt_2):
            # Calculate lattice strained to match:
            matched_slab_substrate, matched_slab_film = strain_slabs(slab_substrate, slab_film)
        else:
            matched_slab_substrate = slab_substrate
            matched_slab_film = slab_film

        # Ensure substrate has positive c-direction:
        if matched_slab_substrate.lattice.matrix[2, 2] < 0:
            latt = matched_slab_substrate.lattice.matrix.copy()
            latt[2, 2] *= -1
            new_struct = matched_slab_substrate.copy()
            new_struct.lattice = Lattice(latt)
            matched_slab_substrate = new_struct

        # Ensure film has positive c-direction:
        if matched_slab_film.lattice.matrix[2, 2] < 0:
            latt = matched_slab_film.lattice.matrix.copy()
            latt[2, 2] *= -1
            new_struct = matched_slab_film.copy()
            new_struct.lattice = Lattice(latt)
            matched_slab_film = new_struct

        if offset is None:
            offset = (2.5, 0.0, 0.0)

        _structure = merge_slabs(matched_slab_substrate, matched_slab_film, *offset)
        orthogonal_structure = _structure.get_orthogonal_c_slab()
        orthogonal_structure.sort()

        if not orthogonal_structure.is_valid(tol=1):
            warnings.warn("Check generated structure, it may contain atoms too closely placed")

        # offset_vector = (offset[1], offset[2], offset[0])
        interface = Interface(
            orthogonal_structure.lattice.copy(),
            orthogonal_structure.species,
            orthogonal_structure.frac_coords,
            slab_substrate.miller_index,
            slab_film.miller_index,
            self.original_substrate_structure,
            self.original_film_structure,
            unstrained_slab_substrate,
            unstrained_slab_film,
            slab_substrate,
            slab_film,
            init_inplane_shift=offset[1:],
            site_properties=orthogonal_structure.site_properties,
        )

        return interface

    def visualize_interface(self, interface_index=0, show_atoms=False, n_uc=2):
        """
        Plot the film-substrate superlattice match, the film superlattice,
        and the substrate superlattice in three separate plots and show them.

        Args:
            interface_index (int, 0): Choice of interface to plot
            show_atoms (bool, False): Whether to plot atomic sites
            n_uc (int, 2): Number of 2D unit cells of the interface in each direction.
                (The unit cell of the interface is the supercell of th substrate
                that matches a supercel of the film.)
        """
        film_index = int(self.interface_labels[interface_index][0])
        sub_index = int(self.interface_labels[interface_index][2])
        visualize_interface(self.interfaces[interface_index], show_atoms, n_uc)
        visualize_superlattice(
            self.film_structures[film_index],
            self.modified_film_structures[film_index],
            film=True,
            show_atoms=show_atoms,
            n_uc=n_uc,
        )
        visualize_superlattice(
            self.substrate_structures[sub_index],
            self.modified_substrate_structures[sub_index],
            film=False,
            show_atoms=show_atoms,
            n_uc=n_uc,
        )


def visualize_interface(interface, show_atoms=False, n_uc=2):
    """
    Plot the match of the substrate and film superlattices.

    Args:
        interface (Interface): Interface object
        show_atoms (bool, False): Whether to plot atomic sites
        n_uc (int, 2): Number of 2D unit cells of the interface in each direction.
            (The unit cell of the interface is the supercell of th substrate
            that matches a supercel of the film.)
    """
    # sub_struct = interface.sub_init_cell
    # film_struct = interface.film_init_cell
    modified_sub_struct = interface.modified_sub_structure
    modified_film_struct = interface.modified_film_structure
    rotated_modified_film_structure = align_x(modified_film_struct.copy(), get_ortho_axes(modified_sub_struct))

    # Show super lattice matches
    plt.figure(dpi=150)
    legend_elements = []
    for i, j in product(range(-n_uc, n_uc), range(-n_uc, n_uc)):
        v1 = modified_sub_struct.lattice.matrix[0, :]
        v2 = modified_sub_struct.lattice.matrix[1, :]
        current_start = v1 * i + v2 * j
        plt.plot(
            [current_start[0], current_start[0] + v1[0]],
            [current_start[1], current_start[1] + v1[1]],
            "-k",
            linewidth=0.3,
        )
        plt.plot(
            [current_start[0], current_start[0] + v2[0]],
            [current_start[1], current_start[1] + v2[1]],
            "-k",
            linewidth=0.3,
        )
        if show_atoms:
            plt.plot(
                np.add(modified_sub_struct.cart_coords[:, 0], current_start[0]),
                np.add(modified_sub_struct.cart_coords[:, 1], current_start[1]),
                "or",
                markersize=0.1,
            )
    legend_elements.append(Line2D([0], [0], color="k", lw=1, label="Substrate Superlattice"))
    if show_atoms:
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                lw=1,
                label="Substrate atoms",
                markerfacecolor="r",
                markersize=3,
            )
        )

    for i, j in product(range(-n_uc, n_uc), range(-n_uc, n_uc)):
        v1 = rotated_modified_film_structure.lattice.matrix[0, :]
        v2 = rotated_modified_film_structure.lattice.matrix[1, :]
        current_start = v1 * i + v2 * j
        plt.plot(
            [current_start[0], current_start[0] + v1[0]],
            [current_start[1], current_start[1] + v1[1]],
            "-b",
            linewidth=0.3,
        )
        plt.plot(
            [current_start[0], current_start[0] + v2[0]],
            [current_start[1], current_start[1] + v2[1]],
            "-b",
            linewidth=0.3,
        )
        if show_atoms:
            plt.plot(
                np.add(rotated_modified_film_structure.cart_coords[:, 0], current_start[0]),
                np.add(rotated_modified_film_structure.cart_coords[:, 1], current_start[1]),
                "og",
                markersize=0.1,
            )
    legend_elements.append(Line2D([0], [0], color="b", lw=1, label="Film Superlattice"))
    if show_atoms:
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                lw=1,
                label="Film atoms",
                markerfacecolor="g",
                markersize=3,
            )
        )
    plt.axis("scaled")
    plt.title("Superlattice Match")
    plt.legend(handles=legend_elements)
    plt.show()


def visualize_superlattice(struct, modified_struct, film=True, show_atoms=False, n_uc=2):
    """
    Visualize the unit cell-supercell match for either the film or substrate
    (specified by film boolean tag).

    Args:
        struct (Slab): unit cell slab
        modified_struct (Slab): supercell slab
        film (bool, True): True=label plot as film, False=label plot as substrate
        show_atoms (bool, False): Whether to plot atomic sites
        n_uc (int, 2): Number of 2D unit cells of the interface in each direction.
            (The unit cell of the interface is the supercell of th substrate
            that matches a supercel of the film.)
    """
    label = "Film" if film else "Substrate"
    plt.figure(dpi=150)
    legend_elements = []
    for i, j in product(range(-n_uc, n_uc), range(-n_uc, n_uc)):
        v1 = modified_struct.lattice.matrix[0, :]
        v2 = modified_struct.lattice.matrix[1, :]
        current_start = v1 * i + v2 * j
        plt.plot(
            [current_start[0], current_start[0] + v1[0]],
            [current_start[1], current_start[1] + v1[1]],
            "-k",
            linewidth=0.3,
        )
        plt.plot(
            [current_start[0], current_start[0] + v2[0]],
            [current_start[1], current_start[1] + v2[1]],
            "-k",
            linewidth=0.3,
        )
        if show_atoms:
            plt.plot(
                np.add(modified_struct.cart_coords[:, 0], current_start[0]),
                np.add(modified_struct.cart_coords[:, 1], current_start[1]),
                "or",
                markersize=0.1,
            )
    legend_elements.append(Line2D([0], [0], color="k", lw=1, label="%s Superlattice" % label))
    if show_atoms:
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                lw=1,
                label="%s Superlattice atoms" % label,
                markerfacecolor="r",
                markersize=3,
            )
        )

    uc_v1 = struct.lattice.matrix[0, :]
    uc_v2 = struct.lattice.matrix[1, :]
    sl_v1 = modified_struct.lattice.matrix[0, :]
    sl_v2 = modified_struct.lattice.matrix[1, :]

    sl_v = (sl_v1 + sl_v2) * n_uc
    uc_v = (uc_v1 + uc_v2) * n_uc
    rx = np.abs(int(n_uc * sl_v[0] / uc_v[0]))
    ry = np.abs(int(n_uc * sl_v[1] / uc_v[1]))

    for i, j in product(range(-rx, rx), range(-ry, ry)):
        v1 = struct.lattice.matrix[0, :]
        v2 = struct.lattice.matrix[1, :]
        current_start = v1 * i + v2 * j
        plt.plot(
            [current_start[0], current_start[0] + v1[0]],
            [current_start[1], current_start[1] + v1[1]],
            "-b",
            linewidth=0.3,
        )
        plt.plot(
            [current_start[0], current_start[0] + v2[0]],
            [current_start[1], current_start[1] + v2[1]],
            "-b",
            linewidth=0.3,
        )
        if show_atoms:
            plt.plot(
                np.add(struct.cart_coords[:, 0], current_start[0]),
                np.add(struct.cart_coords[:, 1], current_start[1]),
                "og",
                markersize=0.1,
            )
    legend_elements.append(Line2D([0], [0], color="b", lw=1, label="%s Lattice" % label))
    if show_atoms:
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                lw=1,
                label="%s atoms" % label,
                markerfacecolor="g",
                markersize=3,
            )
        )
    plt.axis("scaled")
    plt.legend(handles=legend_elements)
    plt.title("%s unit cell and superlattice" % label)
    plt.show()


def merge_slabs(substrate, film, slab_offset, x_offset, y_offset, vacuum=20, **kwargs):
    """
    Given substrate and film supercells (oriented to match as closely as possible),
    strain the film to match the substrate lattice and combine the slabs.

    Args:
        slab_offset: spacing between the substrate and film
        x_offset y_offset: in-plane displacement of the film in Cartesian coordinates
        vacuum: vacuum buffer above the film

    Returns:
        combined_structure (Slab): A structure with the strained film and substrate
            combined into one structure
    """
    # strain film to match substrate
    new_latt = film.lattice.matrix.copy()
    new_latt[:2, :2] = substrate.lattice.matrix[:2, :2]
    film.lattice = Lattice(new_latt)

    combined_species = [*substrate.species, *film.species]
    if kwargs.get("cell_height"):
        height = kwargs.get("cell_height")
    else:
        added_height = vacuum + slab_offset + film.lattice.c
        height = added_height + substrate.lattice.matrix[2, 2]
    combined_lattice = substrate.lattice.matrix.copy()
    combined_lattice[2, :] *= height / substrate.lattice.matrix[2, 2]

    max_substrate = np.max(substrate.cart_coords[:, 2])
    min_substrate = np.min(film.cart_coords[:, 2])
    offset = max_substrate - min_substrate + slab_offset
    offset_film_coords = [np.add(coord, [x_offset, y_offset, offset]) for coord in film.cart_coords]
    combined_coords = [*substrate.cart_coords, *offset_film_coords]
    combined_site_properties = {}
    for key, item in substrate.site_properties.items():
        combined_site_properties[key] = [
            *substrate.site_properties[key],
            *film.site_properties[key],
        ]
    labels = ["substrate"] * len(substrate) + ["film"] * len(film)
    combined_site_properties["interface_label"] = labels

    combined_structure = Slab(
        lattice=Lattice(combined_lattice),
        species=combined_species,
        coords=combined_coords,
        miller_index=substrate.miller_index,
        oriented_unit_cell=substrate,
        shift=substrate.shift,
        scale_factor=substrate.scale_factor,
        coords_are_cartesian=True,
        energy=substrate.energy,
        reorient_lattice=False,
        to_unit_cell=True,
        site_properties=combined_site_properties,
    )
    return combined_structure


def strain_slabs(sub_slab, film_slab):
    """
    Strain the film_slab to match the sub_slab,
    orient the structures to match each other,
    and return the new matching structures.

    Args:
        sub_slab (Slab): substrate supercell slab
        film_slab (Slab): film supercell slab

    Returns:
        sub_struct (Slab): substrate structure oriented
            to match the film supercell
        film_struct (Slab): film structure strained to match
            the substrate supercell lattice.
    """
    sub_struct = sub_slab.copy()
    latt_1 = sub_struct.lattice.matrix.copy()
    film_struct = align_x(film_slab, get_ortho_axes(sub_struct)).copy()
    latt_2 = film_struct.lattice.matrix.copy()

    # Rotate film so its diagonal matches with the sub's diagonal
    diag_vec = np.add(latt_1[0, :], latt_1[1, :])
    sub_norm_diag_vec = diag_vec / np.linalg.norm(diag_vec)
    sub_b = np.cross(sub_norm_diag_vec, [0, 0, 1])
    sub_matrix = np.vstack([sub_norm_diag_vec, sub_b, [0, 0, 1]])

    diag_vec = np.add(latt_2[0, :], latt_2[1, :])
    film_norm_diag_vec = diag_vec / np.linalg.norm(diag_vec)
    film_b = np.cross(film_norm_diag_vec, [0, 0, 1])
    film_matrix = np.vstack([film_norm_diag_vec, film_b, [0, 0, 1]])

    rotation = np.dot(np.linalg.inv(film_matrix), sub_matrix)
    new_latt = Lattice(np.dot(film_struct.lattice.matrix, rotation))
    film_struct.lattice = new_latt

    # Average the two lattices (Should get equal strain?)
    mean_a = np.mean([film_struct.lattice.matrix[0, :], sub_struct.lattice.matrix[0, :]], axis=0)
    mean_b = np.mean([film_struct.lattice.matrix[1, :], sub_struct.lattice.matrix[1, :]], axis=0)
    new_latt = np.vstack([mean_a, mean_b, sub_struct.lattice.matrix[2, :]])
    sub_struct.lattice = Lattice(new_latt)
    new_latt = np.vstack([mean_a, mean_b, film_struct.lattice.matrix[2, :]])
    film_struct.lattice = Lattice(new_latt)

    return sub_struct, film_struct


def get_ortho_axes(structure):
    """
    Get an orthonormal set of axes for the structure with the first axis
    pointing along the a lattice vector.

    Args:
        structure (Structure)

    Returns:
        3x3 numpy matrix with the axes
    """
    sub_a = structure.lattice.matrix[0, :] / np.linalg.norm(structure.lattice.matrix[0, :])
    sub_c = third_vect(sub_a, structure.lattice.matrix[1, :])

    sub_b = third_vect(sub_c, sub_a)
    sub_b = sub_b / np.linalg.norm(sub_b)

    return np.vstack((sub_a, sub_b, sub_c))


def align_x(slab, orthogonal_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]):
    """
    Align the a lattice vector of slab with the x axis. Optionally specify
    an orthogonal_basis to align according to a different set of axes

    Args:
        slab (Slab): input structure
        orthogonal basis (3x3 numpy matrix): If specified, align with
            orthogonal_basis[0] rather than [1,0,0]

    Returns:
        The slab, which has been aligned with the specified axis in place.
    """
    sub_ortho_axes = get_ortho_axes(slab)
    rotation = transf_mat(sub_ortho_axes, orthogonal_basis)
    new_sub_lattice = Lattice(np.dot(slab.lattice.matrix[0:3], rotation))
    slab.lattice = new_sub_lattice
    return slab


def transf_mat(A, B):
    """
    Get the matrix to transform from the set of axes A
    to the set of axes B.

    Args:
        A (3x3 numpy array): original axis basis
        B (3x3 numpy array): new axis basis

    Returns:
        3x3 numpy array transformation between the bases
    """
    return np.dot(np.linalg.inv(A), B)


def third_vect(a, b):
    """
    Get a unit vector proportional to cross(a, b).

    Args:
        a, b (numpy arrays): 3D vectors.

    Returns:
        unit vector proportional to cross(a, b).
    """
    c = np.cross(a, b)
    return c / np.linalg.norm(c)


def get_shear_reduced_slab(slab):
    """
    Reduce the vectors of the slab plane according to the algorithm in
    substrate_analyzer, then make a new Slab with a Lattice with those
    reduced vectors.

    Args:
        slab (Slab): Slab to reduce

    Returns:
        Slab object of identical structure to the input slab
        but rduced in-plane lattice vectors
    """
    reduced_vectors = reduce_vectors(slab.lattice.matrix[0], slab.lattice.matrix[1])
    new_lattice = Lattice([reduced_vectors[0], reduced_vectors[1], slab.lattice.matrix[2]])
    return Slab(
        lattice=new_lattice,
        species=slab.species,
        coords=slab.cart_coords,
        miller_index=slab.miller_index,
        oriented_unit_cell=slab.oriented_unit_cell,
        shift=slab.shift,
        scale_factor=slab.scale_factor,
        coords_are_cartesian=True,
        energy=slab.energy,
        reorient_lattice=slab.reorient_lattice,
        to_unit_cell=True,
    )
