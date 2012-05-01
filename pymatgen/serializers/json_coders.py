#!/usr/bin/env python

'''
Beta version of general JSON encoders and decoders for pymatgen. Will only
support pymatgen objects version > 1.8.3.

Current support for all core objects that obey the to_dict API, including Site,
PeriodicSite, Structure, Specie, Dos, Lattice, etc. and all Entry and  all
Transformations. Note that lists and dicts of these objects are supported as
well.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 30, 2012"

import json

class PMGJSONEncoder(json.JSONEncoder):
    """
    A Pymatgen Json Encoder which supports the to_dict API.
    
    Usage:
        Instead of the typical json.dumps(object), you may use
        enc = PMGJSONEncoder()
        enc.encode(object)
        instead. 
    """

    def default(self, o):
        try:
            return o.to_dict
        except:
            return json.JSONEncoder.default(self, o)


class PMGJSONDecoder(json.JSONDecoder):
    """
    A Pymatgen Json Decoder which supports the from_dict API. By default, the
    decoder attempts to find a module and name associated with a dict. If found,
    the decoder will generate a Pymatgen as a priority.  If that fails, the
    original decoded dictionary from the string is returned. Note that nested
    lists and dicts containing pymatgen object will be decoded correctly as well.
    
    Usage:
        Instead of the typical json.loads(json_string), you may use
        dec = PMGJSONDecoder()
        dec.decode(object)
        instead.
    """

    def _process_decoded(self, d):
        """
        Recursive method to support decoding dicts and lists containing
        pymatgen objects.
        """
        if isinstance(d, dict):
            if 'module' in d and 'class' in d:
                mod = __import__(d['module'], globals(), locals(), [d['class']], -1)
                if hasattr(mod, d['class']):
                    cls = getattr(mod, d['class'])
                    return cls.from_dict(d)
            else:
                return {self._process_decoded(k):self._process_decoded(v) for k, v in d.items()}
        elif isinstance(d, list):
            return [self._process_decoded(x) for x in d]
        return d

    def decode(self, s):
        d = json.JSONDecoder.decode(self, s)
        return self._process_decoded(d)

