"""
Usage:
    get_pourbaix.py [options]

Options:
    -m, --mp_id         material_id
    -o, --output_file   output filename
"""

import itertools
import sys

import numpy as np
from scipy.spatial import ConvexHull
from tqdm import tqdm
from docopt import docopt
import matplotlib
matplotlib.use("Agg") # Use Agg backend for macosx

from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter,\
    MultiEntry
from pymatgen.core.periodic_table import Element
from pymatgen import MPRester, Composition


ELEMENTS_HO = {Element("H"), Element("O")}
MPR = MPRester()


# TODO: a lot of duplicated code
def get_hull_in_nph_nphi_space(pourbaix_entries):
    # Get non-OH elements
    pbx_elts = set(itertools.chain.from_iterable(
        [entry.composition.elements for entry in pourbaix_entries]))
    pbx_elts = list(pbx_elts - ELEMENTS_HO)
    dim = len(pbx_elts) - 1

    # Process ions:
    solid_entries = [entry for entry in pourbaix_entries
                     if entry.phase_type == "Solid"]
    ion_entries = [entry for entry in pourbaix_entries
                   if entry.phase_type == "Ion"]

    # If a conc_dict is specified, override individual entry concentrations
    for entry in ion_entries:
        ion_elts = list(set(entry.composition.elements) - ELEMENTS_HO)
        # TODO: the logic here for ion concentration setting is in two
        #       places, in PourbaixEntry and here, should be consolidated
        if len(ion_elts) == 1:
            # Just use default concentration for now
            entry.concentration = 1e-06 * entry.normalization_factor
        elif len(ion_elts) > 1 and not entry.concentration:
            raise ValueError("Elemental concentration not compatible "
                             "with multi-element ions")

    entries = solid_entries + ion_entries

    vecs = [[entry.npH, entry.nPhi, entry.energy] +
            [entry.composition.get(elt) for elt in pbx_elts[:-1]]
            for entry in entries]
    vecs = np.array(vecs)
    norms = np.transpose([[entry.normalization_factor
                           for entry in entries]])
    vecs *= norms
    maxes = np.max(vecs[:, :3], axis=0)
    extra_point = np.concatenate([maxes, np.ones(dim) / dim], axis=0)

    # Add padding for extra point
    pad = 1000
    extra_point[2] += pad
    points = np.concatenate([vecs, np.array([extra_point])], axis=0)
    hull = ConvexHull(points, qhull_options="QJ i")
    # Create facets and remove top
    facets = [facet for facet in hull.simplices
              if not len(points) - 1 in facet]
    return facets


def get_multientries(pourbaix_entries, comp_dict=None):
    """
    Generates multi-entries for pourbaix diagram

    Args:
        pourbaix_entries ([PourbaixEntry]): list of PourbaixEntries to preprocess
            into MultiEntries
        comp_dict ({Element: float}): composition dictionary

    Returns:
        ([MultiEntry]) list of stable MultiEntry candidates

    """
    # Get non-OH elements
    pbx_elts = set(itertools.chain.from_iterable(
        [entry.composition.elements for entry in pourbaix_entries]))
    pbx_elts = list(pbx_elts - ELEMENTS_HO)
    dim = len(pbx_elts) - 1
    comp_dict = comp_dict or {k: 1 / (dim + 1) for k in pbx_elts}

    # Process ions:
    solid_entries = [entry for entry in pourbaix_entries
                   if entry.phase_type == "Solid"]
    ion_entries = [entry for entry in pourbaix_entries
                   if entry.phase_type == "Ion"]

    # If a conc_dict is specified, override individual entry concentrations
    for entry in ion_entries:
        ion_elts = list(set(entry.composition.elements) - ELEMENTS_HO)
        # TODO: the logic here for ion concentration setting is in two
        #       places, in PourbaixEntry and here, should be consolidated
        if len(ion_elts) == 1:
            # Just use default concentration for now
            entry.concentration = 1e-06 * entry.normalization_factor
        elif len(ion_elts) > 1 and not entry.concentration:
            raise ValueError("Elemental concentration not compatible "
                             "with multi-element ions")

    entries = solid_entries + ion_entries

    vecs = [[entry.npH, entry.nPhi, entry.energy] +
            [entry.composition.get(elt) for elt in pbx_elts[:-1]]
            for entry in entries]
    vecs = np.array(vecs)
    norms = np.transpose([[entry.normalization_factor
                           for entry in entries]])
    vecs *= norms
    maxes = np.max(vecs[:, :3], axis=0)
    extra_point = np.concatenate([maxes, np.ones(dim) / dim], axis=0)

    # Add padding for extra point
    pad = 1000
    extra_point[2] += pad
    points = np.concatenate([vecs, np.array([extra_point])], axis=0)
    hull = ConvexHull(points, qhull_options="QJ i")

    # Create facets and remove top
    facets = [facet for facet in hull.simplices
              if not len(points) - 1 in facet]
    facets = get_hull_in_nph_nphi_space(entries)
    if dim > 1:
        valid_facets = []
        for facet in facets:
            comps = vecs[facet][:, 3:]
            full_comps = np.concatenate([
                comps, 1 - np.sum(comps, axis=1).reshape(len(comps), 1)], axis=1)
            # Ensure an compositional interior point exists in the simplex
            if np.linalg.matrix_rank(full_comps) > dim:
                valid_facets.append(facet)

    else:
        valid_facets = facets
    combos = []
    for facet in valid_facets:
        for i in range(1, dim + 2):
            combos.append([
                frozenset(combo) for combo in itertools.combinations(facet, i)])

    all_combos = set(itertools.chain.from_iterable(combos))
    multi_entries = []
    for combo in tqdm(all_combos):
        these_entries = [entries[i] for i in combo]
        mentry = PourbaixDiagram.process_multientry(
            these_entries, Composition(comp_dict))
        if mentry:
            multi_entries.append(mentry)

    return multi_entries


def get_pourbaix_composition(comp):
    """
    Quick routine to get non-OH composition

    Args:
        comp (Composition): composition to get non-OH composition from
    """
    new_comp = comp - Composition({el: comp.get(el) for el in ['H', 'O']})
    total = sum(new_comp.values())
    new_comp = Composition({el: comp.get(el) / total for el in new_comp})
    return new_comp


def get_pourbaix_diagram_from_mp_id(mp_id, conc_dict=None):
    """

    Args:
        mp_id: material-id to retrieve
        conc_dict: dictionary of element-specific concentration values

    Returns:
        Pourbaix diagram constructed from material-id (i.e. with a
            composition from the queried material)
        List of single PourbaixEntries
    """
    # Get material composition
    struct = MPR.get_structure_by_material_id(mp_id)
    pbx_comp = get_pourbaix_composition(struct.composition)

    # Get pourbaix entries, this is still a bit slow from querying
    chemsys = [str(el) for el in pbx_comp]
    pbx_entries = MPR.get_pourbaix_entries(chemsys)

    # Modify concentrations if a conc dict is specified
    if conc_dict is not None:
        for entry in pbx_entries:
            ion_elts = list(set(entry.composition.elements) - ELEMENTS_HO)
            if entry.phase is "Ion":
                entry.concentration = conc_dict[ion_elts[0].symbol] \
                                      * entry.normalization_factor

    # Get multientries
    multi_entries = get_multientries(pbx_entries)

    # Get pourbaix diagram
    pourbaix_diagram = PourbaixDiagram(multi_entries)

    return pourbaix_diagram, pbx_entries


# TODO: more customization options
def get_material_stability_mesh(material_id, conc_dict=None):
    """
    Convenience method to get material stability mesh

    Args:
        material_id (str): material id to be analyzed
        conc_dict (dict): dictionary of concentrations to apply to ions

    Returns:
        mesh of pH/V/dG corresponding to Pourbaix stability

    """
    pourbaix_diagram, pbx_entries = get_pourbaix_diagram_from_mp_id(
        material_id, conc_dict)

    # Retrieve entry and make into a multientry
    entry = [entry for entry in pbx_entries if entry.entry_id == material_id][0]
    mentry = MultiEntry([entry])

    # Specify the pH range
    ph, v = np.mgrid[0:14:100j, -3:3:100j]
    decomp_energy = pourbaix_diagram.get_decomposition_energy(mentry, ph, v)
    return ph, v, decomp_energy, pourbaix_diagram


def get_material_stability_plot(material_id, conc_dict=None, **kwargs):
    """
    Convenience method for plotting pourbaix stability heatmaps

    Args:
        material_id (str): material id corresponding to entry to be analyzed
        conc_dict (dict): dictionary of concentrations to apply to ions
        **kwargs: kwargs for plot_entry_stability

    Returns:
        pyplot object with which plots can be made

    """
    pourbaix_diagram, pbx_entries = get_pourbaix_diagram_from_mp_id(
        material_id, conc_dict)
    plotter = PourbaixPlotter(pourbaix_diagram)
    # Retrieve entry and make into a multientry
    entry = [entry for entry in pbx_entries if entry.entry_id == material_id][0]
    mentry = MultiEntry([entry])

    plt = plotter.plot_entry_stability(mentry, **kwargs)
    return plt


if __name__ == "__main__":
    # TODO: fix options
    mp_id = sys.argv[1]
    stability_plot = get_material_stability_plot(mp_id)
    stability_plot.savefig('{}.png'.format(mp_id))
