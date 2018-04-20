# coding: utf-8
"""
A class for performing analysis of chemical potentials with the grand
canonical linear programming approach
"""
from __future__ import division


from pymatgen.core import Element
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram

def easy_retrieve_chempots_from_MP( defect_entries):
    #make phase diagram out of MaterialsProject data, given elements in the defect_entry list
    elt_list = []
    for defect in defect_entries:
        struct = defect.generate_defect_structure()
        bulk_composition = defect.bulk_structure.composition
        for elt in struct.composition.elements:
            elt_list.append(elt.symbol)

    elt_list = list(set(elt_list))
    with MPRester() as mp:
        entries = mp.get_entries_in_chemsys(elt_list)

    pd = PhaseDiagram(entries)
    return pd.get_all_chempots(bulk_composition)


def easy_retrieve_chempots_from_mpid( mpid, sub_elts=[], full_sub_approach=False):
    """

    :param mpid (str): Mp-id for bulk system you want
    :param sub_elts (list): additional element symbols beyond the bulk composition that you want chemical potentials for
    :param full_sub_approach (Bool): generate chemical potentials by looking at full phase diagram with all sub_elts in it
        (setting to True is really NOT recommended if sub_elts set has more than one element in it...)
            The default approach (full_sub_approach = False) is to only consider facets defined by N-2 phases with
            strictly elements from the BULK composition, and 1 sub_element(+possibly bulk-composition
            element) derived phases (along with the condition for all chemical potentials to be
            defined by the bulk entry, creating N equations to be solved for N atomic chemical potentials).
            This speeds up analysis SIGNFICANTLY when analyzing several substitutional species at once.
            It is essentially the assumption the majority of the elements in the total composition will
            be from the native species present rather than the sub species (a good approximation).
    :return:
    """
    with MPRester() as mp:
        bulk_entry = mp.get_entry_by_material_id(mpid)

    species_symbols = [elt.symbol for elt in bulk_entry.composition.elements]

    if full_sub_approach:
        for elt in sub_elts:
            species_symbols.append(elt)

        with MPRester() as mp:
            entries = mp.get_entries_in_chemsys(species_symbols)

        pd = PhaseDiagram(entries)
        return pd.get_all_chempots(bulk_entry.composition)
    else:
        with MPRester() as mp:
            entries = mp.get_entries_in_chemsys(species_symbols)

        pd = PhaseDiagram(entries)
        out_chempots = pd.get_all_chempots(bulk_entry.composition)
        for bulkfacetname in out_chempots.keys():
            out_chempots[bulkfacetname]['names-to-append'] = []

        for elt in sub_elts:
            this_elt_species_symbols = species_symbols[:]
            this_elt_species_symbols.append(elt)
            with MPRester() as mp:
                entries = mp.get_entries_in_chemsys(this_elt_species_symbols)

            pd = PhaseDiagram(entries)
            this_elt_chempots = pd.get_all_chempots(bulk_entry.composition)
            for facetname, facetcps in this_elt_chempots.items():
                eltcount = facetname.count(elt)
                if eltcount == 1: #only want to include substitution elt-poor regions of phase diagram
                    store_bulk_names = []
                    for possname in facetname.split('-'):
                        if elt in possname:
                            name_to_append = possname
                        else:
                            store_bulk_names.append(possname)
                    for bulkfacetname in out_chempots.keys():
                        if set(store_bulk_names) == set(bulkfacetname.split('-')):
                            out_chempots[bulkfacetname]['names-to-append'].append(name_to_append)
                            out_chempots[bulkfacetname][Element(elt)] = facetcps[Element(elt)]

        final_chempots = {}
        for bulkfacetname, facetcps in out_chempots.items():
            newfacetname = bulkfacetname
            for name_to_append in out_chempots[bulkfacetname]['names-to-append']:
                newfacetname += '-'+str(name_to_append)
            final_chempots[newfacetname] = facetcps
            del final_chempots[newfacetname]['names-to-append']

        return final_chempots

