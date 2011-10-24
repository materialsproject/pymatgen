# -*- coding: iso-8859-15 -*-
"""
Generate a materials collection from the tasks collection.

This is used in the Materials Project.
See http://materials.nersc.gov/ for more information.

Author: Dan Gunter
Created: September 2011
"""
__rcsid__ = '$Id$'


def mapper(doc):
    return doc['crystal_id']

fields = ['elements', 'task_id', 'crystal_id', 
          'space_group', 'output',
	  'reduced_cell_formula', 'pretty_formula',
          'anonymous_formula', 'analysis',
          'oxidation_states','coordination_numbers',
          'is_ordered', 'is_hubbard']

def extractor(doc):
    """Extract only those fields needed for searching and reducing.

    Returns:
      (dict) flattened dictionary with necessary key/value pairs
    """
    r = { 'task_id' : doc['task_id'],
          'crystal_id' : doc['crystal_id'],
          'elements' : doc['elements'],
          'nelements' : len(doc['elements']),
          'reduced_cell_formula' : doc['reduced_cell_formula'],
          'pretty_formula' : doc['pretty_formula'],
          'anonymous_formula' : doc.get('anonymous_formula', ''),
          'crystal_system' : doc['space_group']['crystal_system'],
          'icsd_name' : doc['space_group']['icsd_name'],
          'volume' : doc['output']['crystal']['lattice']['volume'],
          'oxidation_states' : doc['oxidation_states'],
          'coordination_numbers' : doc.get('coordination_numbers',[ ]),
          'is_ordered':doc['is_ordered'],
          'is_hubbard' : doc['is_hubbard'],
          'free_energy_per_atom':doc['output']['free_energy_per_atom'] }
    if doc.has_key('analysis'):
        r.update({'nsites' : doc['analysis'].get('nsites',-1),
                  'density' : doc['analysis'].get('density',-1.0),
                  'e_above_hull' : doc['analysis'].get('e_above_hull',-12345678.0),
                  'formation_energy_per_atom':
                  doc['analysis'].get('formation_energy_per_atom', -1.0)})
    else:
        r.update({'nsites':-1, 'density':-1.0, 'e_above_hull':-12345678.0,
                  'formation_energy_per_atom':-1.0})
    return (r)

def reducer(values):
    """Determine and promote values for 'blessed' task.

    Algorithm for determining the blessed one is:
    does one have a 'U'? is_hubbard -> choose it
    otherwise, lowest energy min(.output.free_energy_per_atom)

    Returns:
      (dict) with keys 'task_id' 'task_ids', and all the other
             keys from the main task identifies with task_id.
    """
    blessed = values[0]
    for v in values:
        if v['is_hubbard']:
            blessed = v
            break
        if v['free_energy_per_atom'] < blessed['free_energy_per_atom']:
            blessed = v # but keep loooking..
    result = blessed
    # add ids of all non-blessed tasks
    result['task_ids'] = [v['task_id'] for v in values]
    result['ntask_ids'] = len(values)
    return result
