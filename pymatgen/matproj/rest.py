#!/usr/bin/env python

'''
Created on Jun 8, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jun 8, 2012"

import urllib
import urllib2

from pymatgen.serializers.json_coders import PMGJSONDecoder
import json


class MPRestAdaptor(object):
    """
    A class to conveniently interface with the Materials Project REST interface.
    
    TODO
    dos - Vasp DOS
    bandstructure - Vasp Bandstructure
    bandgap - Vasp bandgap
    """

    supported_properties = ["structure", "initial_structure", "final_structure", "energy", "energy_per_atom", "formation_energy_per_atom",
         "nsites", "formula", "pretty_formula", "is_hubbard",
         "elements", "nelements", "e_above_hull", "hubbards", "is_compatible", "entry"]

    def __init__(self, api_key):
        self.url = "http://127.0.0.1:8000/rest"
        self.api_key = api_key

    def get_data(self, formula_or_material_id, prop):
        """
        Flexible method to get any data using the Materials Project REST
        interface.
        
        Format of REST return is *always* a list of dict (regardless of the
        number of pieces of data returned. The general format is as follows:
        
        [{
            'material_id': material_id,
            'property_name' : value
        },
        ...
        ]
        """
        url = "{}/{}/vasp/{}?API_KEY={}".format(self.url, formula_or_material_id, prop, self.api_key)
        req = urllib2.Request(url)
        response = urllib2.urlopen(req)
        data = response.read()
        data = json.loads(data, cls=PMGJSONDecoder)
        if data['valid_response']:
            return data['response']
        else:
            print data['error']

    def get_structure_by_material_id(self, material_id, final=True):
        prop = "final_structure" if final else "initial_structure"
        data = self.get_data(material_id, prop=prop)
        return data[0]['structure']

    def get_entry_by_material_id(self, material_id):
        data = self.get_data(material_id, prop="entry")
        return data[0]["entry"]

    def mpquery(self, criteria, properties):
        params = urllib.urlencode({'criteria': criteria, 'properties': properties, 'API_KEY':self.api_key})
        req = urllib2.Request("{}/mpquery".format(self.url), params)
        response = urllib2.urlopen(req)
        data = response.read()
        data = json.loads(data)
        return data


if __name__ == "__main__":

    adaptor = MPRestAdaptor("test_dev")

    for prop in MPRestAdaptor.supported_properties:
        print "{} = {}".format(prop, adaptor.get_data(19017, prop))

    data = adaptor.mpquery(criteria={'elements':{'$in':['Li', 'Na', 'K'], '$all': ['O']}},
                    properties=['formula', 'energy'])

    print data['response']['num_results']
    print len(data['response']['results'])

    print adaptor.get_data(19017, '')
    print adaptor.get_data('Fe2O3', '')

    print adaptor.get_structure_by_material_id(1)
    print adaptor.get_structure_by_material_id(1, False)
    print adaptor.get_entry_by_material_id(1)
