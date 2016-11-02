# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import numpy as np
from pymatgen import Structure, Lattice
from pymatgen.phonon.bandstructure import PhBandStructureSymmLine
from pymatgen.phonon.dos import PhDos, CompletePhDos
from monty.serialization import loadfn


def get_structure_from_dict(d):
    """
    Extracts a structure from the dictionary extracted from the output
    files of phonopy like phonopy.yaml or band.yaml.
    """

    species = []
    frac_coords = []
    for p in d['points']:
        species.append(p['symbol'])
        frac_coords.append(p['coordinates'])

    return Structure(d['lattice'], species, frac_coords)


def get_ph_bs_symm_line_from_dict(bands_dict, has_nac=False, labels_dict=None):
    """
    Extracts a pymatgen PhBandStructure object from the dictionary
    extracted by the band.yaml file produced by phonopy. The labels
    will be extracted from the band.yaml if present.
    Args:
        bands_dict: the dictionary extracted from the band.yaml file
        has_nac: True if the data have been obtained with the option
            --nac option. Default False.
        labels_dict: dict that links a kpoint in frac coords to a label.
            Its value will replace the data contained in the band.yaml.
    """

    qpts = []
    frequencies = []
    phonopy_labels_dict = {}
    for p in bands_dict['phonon']:
        qpts.append(p['q-position'])
        bands = []
        for b in p['band']:
            bands.append(b['frequency'])
        frequencies.append(bands)
        if 'label' in p:
            phonopy_labels_dict[p['label']] = p['q-position']

    qpts = np.array(qpts)
    # transpose to match the convention in PhBandStructure
    frequencies = np.transpose(frequencies)

    structure = get_structure_from_dict(bands_dict)

    rec_latt = Lattice(bands_dict['reciprocal_lattice'])

    labels_dict = labels_dict or phonopy_labels_dict

    ph_bs = PhBandStructureSymmLine(qpts, frequencies, rec_latt, has_nac=has_nac,
                                    labels_dict=labels_dict, structure=structure)

    return ph_bs


def get_ph_bs_symm_line(bands_path, has_nac=False, labels_dict=None):
    """
    Creates a pymatgen PhBandStructure from a band.yaml file.
    Args:
        bands_path: path to the band.yaml file
        has_nac: True if the data have been obtained with the option
            --nac option. Default False.
        labels_dict: dict that links a kpoint in frac coords to a label.
    """
    return get_ph_bs_symm_line_from_dict(loadfn(bands_path), has_nac, labels_dict)


def get_ph_dos(total_dos_path):
    """
    Creates a pymatgen PhDos from a total_dos.dat file.
    Args:
        total_dos_path: path to the total_dos.dat file.
    """
    a = np.loadtxt(total_dos_path)
    return PhDos(a[:,0], a[:,1])


def get_complete_ph_dos(partial_dos_path, phonopy_yaml_path):
    """
    Creates a pymatgen CompletePhDos from a partial_dos.dat and phonopy.yaml files.
    The second is produced when generating a Dos and is needed to extract
    the structure.
    Args:
        partial_dos_path: path to the partial_dos.dat file.
        phonopy_yaml_path: path to the phonopy.yaml file.
    """
    a = np.loadtxt(partial_dos_path).transpose()
    d = loadfn(phonopy_yaml_path)

    structure = get_structure_from_dict(d['primitive_cell'])

    total_dos = PhDos(a[0], a[1:].sum(axis=0))

    pdoss = {}
    for site, pdos in zip(structure, a[1:]):
        pdoss[site] = pdos

    return CompletePhDos(structure, total_dos, pdoss)
