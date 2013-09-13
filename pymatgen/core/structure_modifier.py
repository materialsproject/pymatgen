#!/usr/bin/env python

"""
This module provides classes used to modify structures.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import abc
import itertools
import warnings
import collections

import numpy as np

from pymatgen.util.decorators import deprecated
from pymatgen.core.periodic_table import Specie, Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite, Site
from pymatgen.core.structure import Structure, Molecule
from pymatgen.util.coord_utils import get_points_in_sphere_pbc


class StructureModifier(object):
    """
    Abstract class definition for all classes that modify structures.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def modified_structure(self):
        """
        Returns the modified structure.
        """
        return

    @abc.abstractproperty
    def original_structure(self):
        """
        Returns the original structure.
        """
        return


@deprecated(replacement=Structure)
class StructureEditor(StructureModifier):
    """
    Editor for adding, removing and changing sites from a structure
    """
    DISTANCE_TOLERANCE = 0.01

    def __init__(self, structure):
        """
        Args:
            structure:
                pymatgen.core.structure Structure object.
        """
        self._original_structure = structure
        self._lattice = structure.lattice
        self._sites = list(structure.sites)

    def add_site_property(self, property_name, values):
        """
        Adds a property to a site.

        Args:
            property_name:
                The name of the property to add.
            values:
                A sequence of values. Must be same length as number of sites.
        """
        if len(values) != len(self._sites):
            raise ValueError("Values must be same length as sites.")
        for i in xrange(len(self._sites)):
            site = self._sites[i]
            props = site.properties
            if not props:
                props = {}
            props[property_name] = values[i]
            self._sites[i] = PeriodicSite(site.species_and_occu,
                                          site.frac_coords, self._lattice,
                                          properties=props)

    def replace_species(self, species_mapping):
        """
        Swap species in a structure.

        Args:
            species_mapping:
                dict of species to swap. Species can be elements too.
                e.g., {Element("Li"): Element("Na")} performs a Li for Na
                substitution. The second species can be a sp_and_occu dict.
                For example, a site with 0.5 Si that is passed the mapping
                {Element('Si): {Element('Ge'):0.75, Element('C'):0.25} } will
                have .375 Ge and .125 C.
        """
        def mod_site(site):
            new_atom_occu = collections.defaultdict(int)
            for sp, amt in site.species_and_occu.items():
                if sp in species_mapping:
                    if isinstance(species_mapping[sp], (Element, Specie)):
                        new_atom_occu[species_mapping[sp]] += amt
                    elif isinstance(species_mapping[sp], dict):
                        for new_sp, new_amt in species_mapping[sp].items():
                            new_atom_occu[new_sp] += amt * new_amt
                else:
                    new_atom_occu[sp] += amt
            return PeriodicSite(new_atom_occu, site.frac_coords, self._lattice,
                                properties=site.properties)

        self._sites = map(mod_site, self._sites)

    def replace_site(self, index, species_n_occu):
        """
        Replace a single site. Takes either a species or a dict of species and
        occupations.

        Args:
            index:
                The index of the site in the _sites list.
            species:
                A species object.
        """
        self._sites[index] = PeriodicSite(species_n_occu,
                                          self._sites[index].frac_coords,
                                          self._lattice,
                                          properties=self._sites[index].
                                          properties)

    def remove_species(self, species):
        """
        Remove all occurrences of a species from a structure.

        Args:
            species:
                species to remove.
        """
        new_sites = []
        for site in self._sites:
            new_sp_occu = {sp: amt for sp, amt in site.species_and_occu.items()
                           if sp not in species}
            if len(new_sp_occu) > 0:
                new_sites.append(PeriodicSite(new_sp_occu, site.frac_coords,
                                              self._lattice,
                                              properties=site.properties))
        self._sites = new_sites

    def append_site(self, species, coords, coords_are_cartesian=False,
                    validate_proximity=True):
        """
        Append a site to the structure at the end.

        Args:
            species:
                species of inserted site
            coords:
                coordinates of inserted site
            fractional_coord:
                Whether coordinates are cartesian. Defaults to False.
            validate_proximity:
                Whether to check if inserted site is too close to an existing
                site. Defaults to True.
        """
        self.insert_site(len(self._sites), species, coords,
                         coords_are_cartesian, validate_proximity)

    def insert_site(self, i, species, coords, coords_are_cartesian=False,
                    validate_proximity=True, properties=None):
        """
        Insert a site to the structure.

        Args:
            i:
                index to insert site
            species:
                species of inserted site
            coords:
                coordinates of inserted site
            coords_are_cartesian:
                Whether coordinates are cartesian. Defaults to False.
            validate_proximity:
                Whether to check if inserted site is too close to an existing
                site. Defaults to True.
        """
        if not coords_are_cartesian:
            new_site = PeriodicSite(species, coords, self._lattice,
                                    properties=properties)
        else:
            frac_coords = self._lattice.get_fractional_coords(coords)
            new_site = PeriodicSite(species, frac_coords, self._lattice,
                                    properties=properties)

        if validate_proximity:
            for site in self._sites:
                if site.distance(new_site) < self.DISTANCE_TOLERANCE:
                    raise ValueError("New site is too close to an existing "
                                     "site!")

        self._sites.insert(i, new_site)

    def delete_site(self, i):
        """
        Delete site at index i.

        Args:
            i:
                index of site to delete.
        """
        del(self._sites[i])

    def delete_sites(self, indices):
        """
        Delete sites with at indices.

        Args:
            indices:
                sequence of indices of sites to delete.
        """
        self._sites = [self._sites[i] for i in range(len(self._sites))
                       if i not in indices]

    def apply_operation(self, symmop):
        """
        Apply a symmetry operation to the structure and return the new
        structure. The lattice is operated by the rotation matrix only.
        Coords are operated in full and then transformed to the new lattice.

        Args:
            symmop:
                Symmetry operation to apply.
        """
        self._lattice = Lattice([symmop.apply_rotation_only(row)
                                 for row in self._lattice.matrix])

        def operate_site(site):
            new_cart = symmop.operate(site.coords)
            new_frac = self._lattice.get_fractional_coords(new_cart)
            return PeriodicSite(site.species_and_occu, new_frac, self._lattice,
                                properties=site.properties)
        self._sites = map(operate_site, self._sites)

    def modify_lattice(self, new_lattice):
        """
        Modify the lattice of the structure.  Mainly used for changing the
        basis.

        Args:
            new_lattice:
                New lattice
        """
        self._lattice = new_lattice
        new_sites = []
        for site in self._sites:
            new_sites.append(PeriodicSite(site.species_and_occu,
                                          site.frac_coords,
                                          self._lattice,
                                          properties=site.properties))
        self._sites = new_sites

    def apply_strain(self, strain):
        """
        Apply an isotropic strain to the lattice.

        Args:
            strain:
                Amount of strain to apply. E.g., 0.01 means all lattice
                vectors are increased by 1%. This is equivalent to
                calling modify_lattice with a lattice with lattice parameters
                that are 1% larger.
        """
        self.modify_lattice(Lattice(self._lattice.matrix * (1 + strain)))

    def translate_sites(self, indices, vector, frac_coords=True):
        """
        Translate specific sites by some vector, keeping the sites within the
        unit cell.

        Args:
            sites:
                List of site indices on which to perform the translation.
            vector:
                Translation vector for sites.
            frac_coords:
                Boolean stating whether the vector corresponds to fractional or
                cartesian coordinates.
        """
        for i in indices:
            site = self._sites[i]
            if frac_coords:
                fcoords = site.frac_coords + vector
            else:
                fcoords = self._lattice.get_fractional_coords(site.coords
                                                              + vector)
            new_site = PeriodicSite(site.species_and_occu, fcoords,
                                    self._lattice, to_unit_cell=True,
                                    coords_are_cartesian=False,
                                    properties=site.properties)
            self._sites[i] = new_site

    def perturb_structure(self, distance=0.1):
        """
        Performs a random perturbation of the sites in a structure to break
        symmetries.

        Args:
            distance:
                distance in angstroms by which to perturb each site.
        """
        def get_rand_vec():
            #deals with zero vectors.
            vector = np.random.randn(3)
            vnorm = np.linalg.norm(vector)
            return vector / vnorm * distance if vnorm != 0 else get_rand_vec()

        for i in range(len(self._sites)):
            self.translate_sites([i], get_rand_vec(), frac_coords=False)

    def add_oxidation_state_by_element(self, oxidation_states):
        """
        Add oxidation states to a structure.

        Args:
            structure:
                pymatgen.core.structure Structure object.
            oxidation_states:
                dict of oxidation states.
                E.g., {"Li":1, "Fe":2, "P":5, "O":-2}
        """
        try:
            for i, site in enumerate(self._sites):
                new_sp = {}
                for el, occu in site.species_and_occu.items():
                    sym = el.symbol
                    new_sp[Specie(sym, oxidation_states[sym])] = occu
                new_site = PeriodicSite(new_sp, site.frac_coords,
                                        self._lattice,
                                        coords_are_cartesian=False,
                                        properties=site.properties)
                self._sites[i] = new_site

        except KeyError:
            raise ValueError("Oxidation state of all elements must be "
                             "specified in the dictionary.")

    def add_oxidation_state_by_site(self, oxidation_states):
        """
        Add oxidation states to a structure by site.

        Args:
            oxidation_states:
                List of oxidation states.
                E.g., [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2]
        """
        try:
            for i, site in enumerate(self._sites):
                new_sp = {}
                for el, occu in site.species_and_occu.items():
                    sym = el.symbol
                    new_sp[Specie(sym, oxidation_states[i])] = occu
                new_site = PeriodicSite(new_sp, site.frac_coords,
                                        self._lattice,
                                        coords_are_cartesian=False,
                                        properties=site.properties)
                self._sites[i] = new_site

        except IndexError:
            raise ValueError("Oxidation state of all sites must be "
                             "specified in the dictionary.")

    def remove_oxidation_states(self):
        """
        Removes oxidation states from a structure.
        """
        for i, site in enumerate(self._sites):
            new_sp = collections.defaultdict(float)
            for el, occu in site.species_and_occu.items():
                sym = el.symbol
                new_sp[Element(sym)] += occu
            new_site = PeriodicSite(new_sp, site.frac_coords,
                                    self._lattice,
                                    coords_are_cartesian=False,
                                    properties=site.properties)
            self._sites[i] = new_site

    def to_unit_cell(self, tolerance=0.1):
        """
        Returns all the sites to their position inside the unit cell.
        If there is a site within the tolerance already there, the site is
        deleted instead of moved.
        """
        new_sites = []
        for site in self._sites:
            if not new_sites:
                new_sites.append(site)
                frac_coords = np.array([site.frac_coords])
                continue
            if len(get_points_in_sphere_pbc(self._lattice, frac_coords,
                                            site.coords, tolerance)):
                continue
            frac_coords = np.append(frac_coords, [site.frac_coords % 1],
                                    axis=0)
            new_sites.append(site.to_unit_cell)
        self._sites = new_sites

    @property
    def original_structure(self):
        """
        The original structure.
        """
        return self._original_structure

    @property
    def modified_structure(self):
        coords = [site.frac_coords for site in self._sites]
        species = [site.species_and_occu for site in self._sites]
        props = {}
        if self._sites[0].properties:
            for k in self._sites[0].properties.keys():
                props[k] = [site.properties[k] for site in self._sites]
        return Structure(self._lattice, species, coords, False,
                         site_properties=props)

    

@deprecated(replacement=Structure)
class SupercellMaker(StructureModifier):
    """
    Makes a supercell.
    """

    def __init__(self, structure, scaling_matrix=((1, 0, 0),
                                                  (0, 1, 0),
                                                  (0, 0, 1))):
        """
        Create a supercell.

        Args:
            structure:
                pymatgen.core.structure Structure object.
            scaling_matrix:
                a matrix of transforming the lattice vectors. Defaults to the
                identity matrix. Has to be all integers. e.g.,
                [[2,1,0],[0,3,0],[0,0,1]] generates a new structure with
                lattice vectors a' = 2a + b, b' = 3b, c' = c where a, b, and c
                are the lattice vectors of the original structure.
        """
        self._original_structure = structure
        old_lattice = structure.lattice
        scale_matrix = np.array(scaling_matrix)
        new_lattice = Lattice(np.dot(scale_matrix, old_lattice.matrix))
        new_sites = []

        def range_vec(i):
            return range(max(scale_matrix[:][:, i])
                         - min(scale_matrix[:][:, i]) + 1)

        for site in structure.sites:
            for (i, j, k) in itertools.product(range_vec(0), range_vec(1),
                                               range_vec(2)):
                fcoords = site.frac_coords + np.array([i, j, k])
                coords = old_lattice.get_cartesian_coords(fcoords)
                new_coords = new_lattice.get_fractional_coords(coords)
                new_site = PeriodicSite(site.species_and_occu, new_coords,
                                        new_lattice,
                                        properties=site.properties)
                contains_site = False
                for s in new_sites:
                    if s.is_periodic_image(new_site):
                        contains_site = True
                        break
                if not contains_site:
                    new_sites.append(new_site)
        self._modified_structure = Structure.from_sites(new_sites)

    @property
    def original_structure(self):
        return self._original_structure

    @property
    def modified_structure(self):
        return self._modified_structure


@deprecated(replacement=Structure)
class OxidationStateDecorator(StructureModifier):
    """
    .. deprecated:: v2.1.3

    Use StructureEditor's add_oxidation_state_by... instead.

    Given a dictionary of oxidation states, decorate a structure by replacing
    each Element at a site with a Specie with an oxidation state. Useful for
    higher level functions.
    """

    def __init__(self, structure, oxidation_states):
        """
        Decorates a structure with oxidation states.

        Args:
            structure:
                pymatgen.core.structure Structure object.
            oxidation_states:
                dict of oxidation states.
                E.g., {"Li":1, "Fe":2, "P":5, "O": -2}
        """
        warnings.warn("OxidationStateDecorator has been deprecated. Use "
                      "StructureEditor.remove_oxidation_states instead.")
        self._original_structure = structure
        editor = StructureEditor(structure)
        editor.add_oxidation_state_by_element(oxidation_states)
        self._modified_structure = editor.modified_structure

    @property
    def original_structure(self):
        return self._original_structure

    @property
    def modified_structure(self):
        return self._modified_structure


@deprecated(replacement=Structure)
class OxidationStateRemover(StructureModifier):
    """
    .. deprecated:: v2.1.3

    Use StructureEditor's remove_oxidation_states instead.

    Replace each Specie at a site with an element. Useful for doing structure
    comparisons after applying higher level functions.
    """

    def __init__(self, structure):
        """
        Removes oxidation states from a structure

        Args:
            structure:
                pymatgen.core.structure Structure object.
        """
        warnings.warn("OxidationStateRemover has been deprecated. Use "
                      "StructureEditor.remove_oxidation_states instead.")
        self._original_structure = structure
        new_species = [{Element(el.symbol): occu
                        for el, occu in site.species_and_occu.items()}
                       for site in structure]
        self._modified_structure = Structure(structure.lattice, new_species,
                                             structure.frac_coords, False)

    @property
    def original_structure(self):
        return self._original_structure

    @property
    def modified_structure(self):
        return self._modified_structure


@deprecated(replacement=Structure)
class BasisChange(StructureModifier):
    """
    Given a new basis, we express the structure in this new basis.
    """

    def __init__(self, structure, new_lattice):
        """
        Express a given structure in a new basis.

        Args:
            structure:
                pymatgen.core.structure Structure object.
            new_lattice:
                a pymatgen.core.Lattice object
        """
        self._original_structure = structure
        sp = [site.species_and_occu for site in structure._sites]
        coords = [site.coords for site in structure._sites]
        self._modified_structure = Structure(new_lattice, sp, coords,
                                             validate_proximity=False,
                                             to_unit_cell=True,
                                             coords_are_cartesian=True)

    @property
    def original_structure(self):
        return self._original_structure

    @property
    def modified_structure(self):
        return self._modified_structure


@deprecated(replacement=Molecule)
class MoleculeEditor(StructureModifier):
    """
    Editor for adding, removing and changing sites from a molecule.
    """
    DISTANCE_TOLERANCE = 0.01

    def __init__(self, molecule):
        """
        Args:
            molecule:
                pymatgen.core.structure Molecule object.
        """
        self._original_structure = molecule
        self._sites = list(molecule.sites)

    def add_site_property(self, property_name, values):
        """
        Adds a property to a site.

        Args:
            property_name:
                The name of the property to add.
            values:
                A sequence of values. Must be same length as number of sites.
        """
        if len(values) != len(self._sites):
            raise ValueError("Values must be same length as sites.")
        for i in xrange(len(self._sites)):
            site = self._sites[i]
            props = site.properties
            if not props:
                props = {}
            props[property_name] = values[i]
            self._sites[i] = Site(site.species_and_occu, site.coords,
                                  properties=props)

    def replace_species(self, species_mapping):
        """
        Swap species in a molecule.

        Args:
            species_mapping:
                dict of species to swap. Species can be elements too.
                e.g., {Element("Li"): Element("Na")} performs a Li for Na
                substitution. The second species can be a sp_and_occu dict.
                For example, a site with 0.5 Si that is passed the mapping
                {Element('Si): {Element('Ge'):0.75, Element('C'):0.25} } will
                have .375 Ge and .125 C.
        """

        def mod_site(site):
            new_atom_occu = dict()
            for sp, amt in site.species_and_occu.items():
                if sp in species_mapping:
                    if isinstance(species_mapping[sp], (Element, Specie)):
                        if species_mapping[sp] in new_atom_occu:
                            new_atom_occu[species_mapping[sp]] += amt
                        else:
                            new_atom_occu[species_mapping[sp]] = amt
                    elif isinstance(species_mapping[sp], dict):
                        for new_sp, new_amt in species_mapping[sp].items():
                            if new_sp in new_atom_occu:
                                new_atom_occu[new_sp] += amt * new_amt
                            else:
                                new_atom_occu[new_sp] = amt * new_amt
                else:
                    if sp in new_atom_occu:
                        new_atom_occu[sp] += amt
                    else:
                        new_atom_occu[sp] = amt
            return Site(new_atom_occu, site.coords, properties=site.properties)
        self._sites = map(mod_site, self._sites)

    def replace_site(self, index, species_n_occu):
        """
        Replace a single site. Takes either a species or a dict of occus.

        Args:
            index:
                The index of the site in the _sites list
            species:
                A species object.
        """
        self._sites[index] = Site(species_n_occu, self._sites[index].coords,
                                  properties=self._sites[index].properties)

    def remove_species(self, species):
        """
        Remove all occurrences of a species from a molecule.

        Args:
            species:
                Species to remove.
        """
        new_sites = []
        for site in self._sites:
            new_sp_occu = {sp: amt for sp, amt in site.species_and_occu.items()
                           if sp not in species}
            if len(new_sp_occu) > 0:
                new_sites.append(Site(new_sp_occu, site.coords,
                                      properties=site.properties))
        self._sites = new_sites

    def append_site(self, species, coords, validate_proximity=True):
        """
        Append a site to the structure at the end.

        Args:
            species:
                species of inserted site.
            coords:
                coordinates of inserted site.
            validate_proximity:
                Whether to check if inserted site is too close to an existing
                site. Defaults to True.
        """
        self.insert_site(len(self._sites), species, coords, validate_proximity)

    def insert_site(self, i, species, coords, validate_proximity=True,
                    properties=None):
        """
        Insert a site to the structure.

        Args:
            i:
                Index to insert site.
            species:
                Species of inserted site.
            coords:
                Coordinates of inserted site.
            validate_proximity:
                Whether to check if inserted site is too close to an existing
                site. Defaults to True.
        """
        new_site = Site(species, coords, properties=properties)

        if validate_proximity:
            for site in self._sites:
                if site.distance(new_site) < self.DISTANCE_TOLERANCE:
                    raise ValueError("New site is too close to an existing "
                                     "site!")
        self._sites.insert(i, new_site)

    def delete_site(self, i):
        """
        Delete site at index i.

        Args:
            i:
                index of site to delete.
        """
        del(self._sites[i])

    def delete_sites(self, indices):
        """
        Delete sites with at indices.

        Args:
            indices:
                sequence of indices of sites to delete.
        """
        self._sites = [self._sites[i] for i in range(len(self._sites))
                       if i not in indices]

    def translate_sites(self, indices, vector):
        """
        Translate specific sites by some vector, keeping the sites within the
        unit cell.

        Args:
            sites:
                List of site indices on which to perform the translation.
            vector:
                Translation vector for sites.
        """
        for i in indices:
            site = self._sites[i]
            new_site = Site(site.species_and_occu, site.coords + vector,
                            properties=site.properties)
            self._sites[i] = new_site

    def perturb_structure(self, distance=0.1):
        """
        Performs a random perturbation of the sites in a structure to break
        symmetries.

        Args:
            distance:
                distance in angstroms by which to perturb each site.
        """
        for i in range(len(self._sites)):
            vector = np.random.rand(3)
            vector /= np.linalg.norm(vector) / distance
            self.translate_sites([i], vector)

    @property
    def original_structure(self):
        return self._original_structure

    @property
    def modified_structure(self):
        coords = [site.coords for site in self._sites]
        species = [site.species_and_occu for site in self._sites]
        props = {}
        if self._sites[0].properties:
            for k in self._sites[0].properties.keys():
                props[k] = [site.properties[k] for site in self._sites]
        return Molecule(species, coords, False, site_properties=props)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
