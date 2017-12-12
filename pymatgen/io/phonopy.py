# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import PhononDos, CompletePhononDos
from monty.serialization import loadfn
from monty.dev import requires
try:
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from phonopy.file_IO import write_disp_yaml
except ImportError:
    Phonopy = None


@requires(Phonopy, "phonopy not installed!")
def get_pmg_structure(phonopy_structure):
    """
    Convert a PhonopyAtoms object to pymatgen Structure object.

    Args:
        phonopy_structure (PhonopyAtoms): A phonopy structure object.

    """

    lattice = phonopy_structure.get_cell()
    frac_coords = phonopy_structure.get_scaled_positions()
    symbols = phonopy_structure.get_chemical_symbols()
    masses = phonopy_structure.get_masses()
    mms = phonopy_structure.get_magnetic_moments()
    mms = mms or [0] * len(symbols)

    return Structure(lattice, symbols, frac_coords,
                     site_properties={"phonopy_masses": masses,
                                      "magnetic_moments": mms})


@requires(Phonopy, "phonopy not installed!")
def get_phonopy_structure(pmg_structure):
    """
    Convert a pymatgen Structure object to a PhonopyAtoms object.

    Args:
        pmg_structure (pymatgen Structure): A Pymatgen structure object.

    """

    symbols = [site.specie.symbol for site in pmg_structure]

    return PhonopyAtoms(symbols=symbols, cell=pmg_structure.lattice.matrix,
                        scaled_positions=pmg_structure.frac_coords)


def get_structure_from_dict(d):
    """
    Extracts a structure from the dictionary extracted from the output
    files of phonopy like phonopy.yaml or band.yaml.
    Adds "phonopy_masses" in the site_properties of the structures.
    Compatible with older phonopy versions.
    """

    species = []
    frac_coords = []
    masses = []
    if 'points' in d:
        for p in d['points']:
            species.append(p['symbol'])
            frac_coords.append(p['coordinates'])
            masses.append(p['mass'])
    elif 'atoms' in d:
        for p in d['atoms']:
            species.append(p['symbol'])
            frac_coords.append(p['position'])
            masses.append(p['mass'])
    else:
        raise ValueError('The dict does not contain structural information')

    return Structure(d['lattice'], species, frac_coords,
                     site_properties={"phonopy_masses": masses})


def eigvec_to_eigdispl(v, q, frac_coords, mass):
    """
    Converts a single eigenvector to an eigendisplacement in the primitive cell
    according to the formula::
        
        exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
    
    Compared to the modulation option in phonopy, here all the additional
    multiplicative and phase factors are set to 1.
    
    Args:
        v: the vector that should be converted. A 3D complex numpy array.
        q: the q point in fractional coordinates
        frac_coords: the fractional coordinates of the atom
        mass: the mass of the atom
    """

    c = np.exp(2j * np.pi * np.dot(frac_coords, q)) / np.sqrt(mass)

    return c*v


def get_ph_bs_symm_line_from_dict(bands_dict, has_nac=False, labels_dict=None):
    """
    Creates a pymatgen PhononBandStructure object from the dictionary
    extracted by the band.yaml file produced by phonopy. The labels
    will be extracted from the dictionary, if present. If the 'eigenvector'
    key is found the eigendisplacements will be calculated according to the 
    formula::
        
        exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
    
    and added to the object.

    Args:
        bands_dict: the dictionary extracted from the band.yaml file
        has_nac: True if the data have been obtained with the option
            --nac option. Default False.
        labels_dict: dict that links a qpoint in frac coords to a label.
            Its value will replace the data contained in the band.yaml.
    """

    structure = get_structure_from_dict(bands_dict)

    qpts = []
    frequencies = []
    eigendisplacements = []
    phonopy_labels_dict = {}
    for p in bands_dict['phonon']:
        q = p['q-position']
        qpts.append(q)
        bands = []
        eig_q = []
        for b in p['band']:
            bands.append(b['frequency'])
            if 'eigenvector' in b:
                eig_b = []
                for i, eig_a in enumerate(b['eigenvector']):
                    v = np.zeros(3, np.complex)
                    for x in range(3):
                        v[x] = eig_a[x][0] + eig_a[x][1]*1j
                    eig_b.append(eigvec_to_eigdispl(
                        v, q, structure[i].frac_coords,
                        structure.site_properties['phonopy_masses'][i]))
                eig_q.append(eig_b)
        frequencies.append(bands)
        if 'label' in p:
            phonopy_labels_dict[p['label']] = p['q-position']
        if eig_q:
            eigendisplacements.append(eig_q)

    qpts = np.array(qpts)
    # transpose to match the convention in PhononBandStructure
    frequencies = np.transpose(frequencies)
    if eigendisplacements:
        eigendisplacements = np.transpose(eigendisplacements, (1, 0, 2, 3))

    rec_latt = Lattice(bands_dict['reciprocal_lattice'])

    labels_dict = labels_dict or phonopy_labels_dict

    ph_bs = PhononBandStructureSymmLine(
        qpts, frequencies, rec_latt, has_nac=has_nac, labels_dict=labels_dict,
        structure=structure, eigendisplacements=eigendisplacements)

    return ph_bs


def get_ph_bs_symm_line(bands_path, has_nac=False, labels_dict=None):
    """
    Creates a pymatgen PhononBandStructure from a band.yaml file.
    The labels will be extracted from the dictionary, if present.
    If the 'eigenvector'  key is found the eigendisplacements will be
    calculated according to the formula:
    \\exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
     and added to the object.

    Args:
        bands_path: path to the band.yaml file
        has_nac: True if the data have been obtained with the option
            --nac option. Default False.
        labels_dict: dict that links a qpoint in frac coords to a label.
    """
    return get_ph_bs_symm_line_from_dict(loadfn(bands_path), has_nac,
                                         labels_dict)


def get_ph_dos(total_dos_path):
    """
    Creates a pymatgen PhononDos from a total_dos.dat file.

    Args:
        total_dos_path: path to the total_dos.dat file.
    """
    a = np.loadtxt(total_dos_path)
    return PhononDos(a[:, 0], a[:, 1])


def get_complete_ph_dos(partial_dos_path, phonopy_yaml_path):
    """
    Creates a pymatgen CompletePhononDos from a partial_dos.dat and
    phonopy.yaml files.
    The second is produced when generating a Dos and is needed to extract
    the structure.

    Args:
        partial_dos_path: path to the partial_dos.dat file.
        phonopy_yaml_path: path to the phonopy.yaml file.
    """
    a = np.loadtxt(partial_dos_path).transpose()
    d = loadfn(phonopy_yaml_path)

    structure = get_structure_from_dict(d['primitive_cell'])

    total_dos = PhononDos(a[0], a[1:].sum(axis=0))

    pdoss = {}
    for site, pdos in zip(structure, a[1:]):
        pdoss[site] = pdos.tolist()

    return CompletePhononDos(structure, total_dos, pdoss)


@requires(Phonopy, "phonopy not installed!")
def get_displaced_structures(pmg_structure, atom_disp=0.01,
                             supercell_matrix=None, yaml_fname=None, **kwargs):
    """
    Generate a set of symmetrically inequivalent displaced structures for
    phonon calculations.

    Args:
        pmg_structure (Structure): A pymatgen structure object.
        atom_disp (float): Atomic displacement. Default is 0.01 $\\AA$.
        supercell_matrix (3x3 array): Scaling matrix for supercell.
        yaml_fname (string): If not None, it represents the full path to
            the outputting displacement yaml file, e.g. disp.yaml.
        **kwargs: Parameters used in Phonopy.generate_displacement method.

    Return:
        A list of symmetrically inequivalent structures with displacements, in
        which the first element is the perfect supercell structure.
    """

    is_plusminus = kwargs.get("is_plusminus", "auto")
    is_diagonal = kwargs.get("is_diagonal", True)
    is_trigonal = kwargs.get("is_trigonal", False)

    ph_structure = get_phonopy_structure(pmg_structure)

    if supercell_matrix is None:
        supercell_matrix = np.eye(3) * np.array((1, 1, 1))

    phonon = Phonopy(unitcell=ph_structure, supercell_matrix=supercell_matrix)
    phonon.generate_displacements(distance=atom_disp,
                                  is_plusminus=is_plusminus,
                                  is_diagonal=is_diagonal,
                                  is_trigonal=is_trigonal)

    if yaml_fname is not None:
        displacements = phonon.get_displacements()
        directions = phonon.get_displacement_directions()
        write_disp_yaml(displacements=displacements,
                        supercell=phonon.get_supercell(),
                        directions=directions, filename=yaml_fname)

    # Supercell structures with displacement
    disp_supercells = phonon.get_supercells_with_displacements()
    # Perfect supercell structure
    init_supercell = phonon.get_supercell()
    # Structure list to be returned
    structure_list = [get_pmg_structure(init_supercell)]

    for c in disp_supercells:
        if c is not None:
            structure_list.append(get_pmg_structure(c))

    return structure_list
