# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import numpy as np
from pymatgen import Structure, Lattice
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import PhononDos, CompletePhononDos
from monty.serialization import loadfn


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
        raise ValueError('The dictionary does not contain structural information')

    return Structure(d['lattice'], species, frac_coords, site_properties={"phonopy_masses": masses})


def eigvec_to_eigdispl(v, q, frac_coords, mass):
    """
    Converts a single eigenvector to an eigendisplacement in the primitive cell
    according to the formula:
    \exp(2*pi*i*(frac_coords \dot q) / sqrt(mass) * v
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
    key is found the eigendisplacements will be calculated according to the formula:
    \exp(2*pi*i*(frac_coords \dot q) / sqrt(mass) * v
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
                    eig_b.append(eigvec_to_eigdispl(v, q, structure[i].frac_coords,
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

    ph_bs = PhononBandStructureSymmLine(qpts, frequencies, rec_latt, has_nac=has_nac, labels_dict=labels_dict,
                                        structure=structure, eigendisplacements=eigendisplacements)

    return ph_bs


def get_ph_bs_symm_line(bands_path, has_nac=False, labels_dict=None):
    """
    Creates a pymatgen PhononBandStructure from a band.yaml file.
    The labels will be extracted from the dictionary, if present.
    If the 'eigenvector'  key is found the eigendisplacements will be
    calculated according to the formula:
    \exp(2*pi*i*(frac_coords \dot q) / sqrt(mass) * v
     and added to the object.

    Args:
        bands_path: path to the band.yaml file
        has_nac: True if the data have been obtained with the option
            --nac option. Default False.
        labels_dict: dict that links a qpoint in frac coords to a label.
    """
    return get_ph_bs_symm_line_from_dict(loadfn(bands_path), has_nac, labels_dict)


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
    Creates a pymatgen CompletePhononDos from a partial_dos.dat and phonopy.yaml files.
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
