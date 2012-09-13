#!/usr/bin/env python

"""
.. versionadded:: 1.9.0

This module implements the abstract base class for msonable pymatgen objects,
i.e., objects that can be converted to a json representation. MSON stands for
materials json.

It also implements general JSON encoders and decoders for pymatgen. Only
supports pymatgen object version >= 1.9.0.

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

        d["module"] = self.__class__.__module__
        d["class"] = self.__class__.__name__

"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 30, 2012"

import json
import abc

from pymatgen.util.io_utils import zopen


class MSONable(object):
    """
    This is an abstract base class specifying an API for msonable objects. MSON
    is Materials JSON. Essentially, MSONable objects must implement a to_dict
    property and a from_dict static method.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def to_dict(self):
        """
        A JSON serializable dict representation of an object.
        """
        pass

    @staticmethod
    def from_dict(d):
        """
        This simply raises a NotImplementedError to force subclasses to
        implement this static method. Abstract static methods are not
        implemented until Python 3+.
        """
        raise NotImplementedError("MSONable objects must implement a from_dict"
                                  " static method.")

    @property
    def to_json(self):
        """
        Returns a json string representation of the MSONable object.
        """
        return json.dumps(self, cls=PMGJSONEncoder)

    def write_to_json_file(self, filename):
        """
        Writes the mson representation to a file.

        Args:
            filename:
                filename to write to. It is recommended that the file extension
                be ".mson".
        """
        with zopen(filename, "wb") as f:
            json.dump(self, f, cls=PMGJSONEncoder)


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
            if "@module" not in d:
                d["@module"] = o.__class__.__module__
            if "@class" not in d:
                d["@class"] = o.__class__.__name__
            return d
        except:
            return json.JSONEncoder.default(self, o)


class PMGJSONDecoder(json.JSONDecoder):
    """
    A Pymatgen Json Decoder which supports the from_dict API. By default, the
    decoder attempts to find a module and name associated with a dict. If
    found, the decoder will generate a Pymatgen as a priority.  If that fails,
    the original decoded dictionary from the string is returned. Note that
    nested lists and dicts containing pymatgen object will be decoded correctly
    as well.

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
            if "@module" in d and "@class" in d:
                modname = d["@module"]
                classname = d["@class"]
            elif "module" in d and "class" in d:
                modname = d["module"]
                classname = d["class"]
            else:
                modname = None
            if modname:
                mod = __import__(modname, globals(), locals(),
                                 [classname], -1)
                if hasattr(mod, classname):
                    cls = getattr(mod, classname)
                    data = {k: v for k, v in d.items() if k not in ["module",
                                                                    "class",
                                                                    "@module",
                                                                    "@class"]}
                    if hasattr(cls, "from_dict"):
                        return cls.from_dict(data)
            return {self.process_decoded(k): self.process_decoded(v) \
                    for k, v in d.items()}
        elif isinstance(d, list):
            return [self.process_decoded(x) for x in d]
        return d

    def decode(self, s):
        d = json.JSONDecoder.decode(self, s)
        return self.process_decoded(d)
