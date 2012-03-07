#!/usr/bin/env python

'''
This module performs list of transformations on list of entries, and provides a 
consistent input vs output interface for transformations on db entries.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 2, 2012"

import os
import re
import datetime
import json

from pymatgen.io.cifio import CifParser
from pymatgen.io.vaspio import Poscar

class TransformedStructure(object):
    """
    Container object for new structures that include history of transformations.
    
    Each transformed structure is made up of a sequence of structures with associated
    transformation history.
    """
    
    def __init__(self, structure, transformations, source_info = None):
        self.source_info = source_info
        self.structures = [structure]
        self.transformations = []
        for t in transformations:
            self.append_transformation(t)
    
    def __getitem__(self, index):
        return (self.structures[index], self.transformations[0:index])
    
    def __len__(self):
        return len(self.structures)
    
    def append_transformation(self, trans):
        new_s = trans.apply_transformation(self.structures[-1])
        self.structures.append(new_s)
        self.transformations.append(trans)
    
    def get_vasp_input(self, vasp_input_set, generate_potcar = True):
        """
        Returns VASP input as a dict of vaspio objects.
        
        Args:
            vasp_input_set:
                pymatgen.io.vaspio_set.VaspInputSet like object that creates
                vasp input files from structures
            generate_potcar:
                Set to False to generate a POTCAR.spec file instead of a POTCAR,
                which contains the POTCAR labels but not the actual POTCAR. Defaults
                to True.
        """
        d = vasp_input_set.get_all_vasp_input(self.structures[-1], generate_potcar)
        d['transformations.json'] = json.dumps(self.to_dict)
        return d
    
    def write_vasp_input(self, vasp_input_set, output_dir, create_directory = True):
        """
        Args:
            vasp_input_set:
                pymatgen.io.vaspio_set.VaspInputSet like object that creates
                vasp input files from structures
            output_dir:
                Directory to output files
            create_directory:
                Create the directory if not present. Defaults to True.
        """
        vasp_input_set.write_input(self.structures[-1], output_dir, make_dir_if_not_present = create_directory)
        with open(os.path.join(output_dir, 'transformations.json'), 'w') as fp:
            json.dump(self.to_dict, fp)
        
    @property
    def final_structure(self):
        return self.structures[-1]
    
    @property
    def to_dict(self):
        d = self.structures[-1].to_dict
        history = []
        history.append(self.source_info)
        for i, t in enumerate(self.transformations):
            tdict = t.to_dict
            tdict['input_structure'] = self.structures[i].to_dict
            history.append(tdict)
        d['history'] = history
        return d


class CifTransformedStructure(TransformedStructure):
    """
    Generates NewMaterials from Cifs.
    """
    
    def __init__(self, cif_string, transformations):
        parser = CifParser.from_string(cif_string)
        raw_string = re.sub("'", "\"", cif_string)
        cif_dict = parser.to_dict
        cif_keys = cif_dict.keys()
        s = parser.get_structures()[0]
        partial_cif = cif_dict[cif_keys[0]]
        if '_database_code_ICSD' in partial_cif:
            source = partial_cif['_database_code_ICSD'] + "-ICSD"
        else:
            source = 'uploaded cif'
        source_info = {'source':source,'datetime':str(datetime.datetime.utcnow()), 'original_file':raw_string, 'cif_data':cif_dict[cif_keys[0]]}
        super(CifTransformedStructure, self).__init__(s, transformations, source_info)


class PoscarTransformedStructure(TransformedStructure):
    """
    Generates NewMaterials from Cifs.
    """
    
    def __init__(self, poscar_string, transformations):
        p = Poscar.from_string(poscar_string)
        if not p.true_names:
            raise ValueError("Transformation can be craeted only from POSCAR strings with proper VASP5 element symbols.")
        raw_string = re.sub("'", "\"", poscar_string)
        s = p.struct
        source_info = {'source': "uploaded POSCAR", 'datetime':str(datetime.datetime.utcnow()), 'original_file':raw_string}
        super(PoscarTransformedStructure, self).__init__(s, transformations, source_info)

