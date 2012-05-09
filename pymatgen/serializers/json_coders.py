#!/usr/bin/env python

'''
.. versionadded:: 1.9.0

General JSON encoders and decoders for pymatgen. Only supports pymatgen objects
version >= 1.9.0.

Current support for all core objects that obey the to_dict/from_dict API,
including Site, PeriodicSite, Structure, Specie, Dos, Lattice, etc. and all
Entry and  all Transformations. Note that nested lists and dicts of these
objects are supported as well.

.. note::

    The decoder depends on finding a "module" and "class" key in the dict in
    order to decode the necessary python object. All to_dict properties must
    therefore have the module name and class embedded. In general, the
    PMGJSONEncoder will add these keys if they are not present, but for better
    long term stability, the easiest way is to add the following to any to_dict
    property::
    
        d['module'] = self.__class__.__module__
        d['class'] = self.__class__.__name__
    
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
        Add it as a *cls* keyword when using json.dump
        json.dumps(object, cls=PMGJSONEncoder)
    """

    def default(self, o):
        try:
            d = o.to_dict
            if "module" not in d:
                d["module"] = o.__class__.__module__
            if "class" not in d:
                d['class'] = o.__class__.__name__
            return d
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
        Add it as a *cls* keyword when using json.load
        json.loads(json_string, cls=PMGJSONDecoder)
    """

    def process_decoded(self, d):
        """
        Recursive method to support decoding dicts and lists containing
        pymatgen objects.
        """
        if isinstance(d, dict):
            if 'module' in d and 'class' in d:

                mod = __import__(d['module'], globals(), locals(), [d['class']], -1)
                if hasattr(mod, d['class']):
                    cls = getattr(mod, d['class'])
                    data = {k:v for k, v in d.items() if k not in ["module", "class"]}
                    return cls.from_dict(data)
            else:
                return {self.process_decoded(k):self.process_decoded(v) for k, v in d.items()}
        elif isinstance(d, list):
            return [self.process_decoded(x) for x in d]
        return d

    def decode(self, s):
        d = json.JSONDecoder.decode(self, s)
        return self.process_decoded(d)

