"""
.. versionadded:: 1.9.0

This module implements the abstract base class for msonable pymatgen objects,
i.e., objects that can be converted to a json representation. MSON stands for
materials json.

It also implements general JSON encoders and decoders for pymatgen. Only
supports pymatgen objects version >= 1.9.0.

Current support for all core objects that obey the to_dict/from_dict API,
including Site, PeriodicSite, Structure, Specie, Dos, Lattice, etc. and all
Entry and  all Transformations. Note that nested lists and dicts of these
objects are supported as well.

.. note::

    The decoder depends on finding a "@module" and "@class" key in the dict in
    order to decode the necessary python object. All to_dict properties must
    therefore have the module name and class embedded. In general, the
    PMGJSONEncoder will add these keys if they are not present, but for better
    long term stability, the easiest way is to add the following to any to_dict
    property::

        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__

"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 30, 2012"

import json
import numpy as np

from abc import ABCMeta, abstractproperty
import datetime

from monty.io import zopen


class MSONable(object):
    """
    This is an abstract base class specifying an API for msonable objects. MSON
    is Materials JSON. Essentially, MSONable objects must implement a to_dict
    property and a from_dict static method.
    """
    __metaclass__ = ABCMeta

    @abstractproperty
    def to_dict(self):
        """
        A JSON serializable dict representation of an object.
        """
        pass

    @classmethod
    def from_dict(cls, d):
        """
        This implements a default from_dict method which supports all
        classes that simply saves all init arguments in a "init_args"
        key. Otherwise, the MSONAble class must override this class method.
        """
        if "init_args" in d:
            return cls(**d['init_args'])
        raise MSONError("Invalid dict for default from_dict. Please "
                        "override from_dict for ".format(cls))

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
        """
        Overriding default method for JSON encoding. This method does two
        things: (a) If an object has a to_dict property, return the to_dict
        output. (b) If the @module and @class keys are not in the to_dict,
        add them to the output automatically. If the object has no to_dict
        property, the default Python json encoder default method is called.

        Args:
            o: Python object.

        Return:
            Python dict representation.
        """
        try:
            if isinstance(o, datetime.datetime):
                return {"@module": "datetime",
                        "@class": "datetime",
                        "string": str(o)}
            elif isinstance(o, np.ndarray):
                return o.tolist()
            elif isinstance(o, np.generic):
                return o.item()

            d = o.to_dict
            if "@module" not in d:
                d["@module"] = o.__class__.__module__
            if "@class" not in d:
                d["@class"] = o.__class__.__name__
            return d
        except AttributeError:
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
                classname = None
            if modname:
                if modname == "datetime" and classname == "datetime":
                    try:
                        dt = datetime.datetime.strptime(d["string"],
                                                        "%Y-%m-%d %H:%M:%S.%f")
                    except ValueError:
                        dt = datetime.datetime.strptime(d["string"],
                                                        "%Y-%m-%d %H:%M:%S")
                    return dt
                mod = __import__(modname, globals(), locals(), [classname], -1)
                if hasattr(mod, classname):
                    cls_ = getattr(mod, classname)
                    data = {k: v for k, v in d.items()
                            if k not in ["module", "class",
                                         "@module", "@class"]}
                    if hasattr(cls_, "from_dict"):
                        return cls_.from_dict(data)
            return {self.process_decoded(k): self.process_decoded(v)
                    for k, v in d.items()}
        elif isinstance(d, list):
            return [self.process_decoded(x) for x in d]

        return d

    def decode(self, s):
        d = json.JSONDecoder.decode(self, s)
        return self.process_decoded(d)


class MSONError(Exception):
    """
    Exception class for serialization errors.
    """
    pass


def json_pretty_dump(obj, filename):
    """
    Serialize obj as a JSON formatted stream to the given filename (
    pretty printing version)
    """
    with open(filename, "w") as fh:
        json.dump(obj, fh, indent=4, sort_keys=4)


def pmg_load(filename, **kwargs):
    """
    Loads a json file and deserialize it with PMGJSONDecoder.

    Args:
        filename (str): Filename of file to open. Can be gzipped or bzipped.
        \*\*kwargs: Any of the keyword arguments supported by the json.load
            method.

    Returns:
        Deserialized pymatgen object. Note that these objects can be lists,
        dicts or otherwise nested pymatgen objects that support the to_dict
        and from_dict MSONAble protocol.
    """
    return json.load(zopen(filename), cls=PMGJSONDecoder, **kwargs)


def pmg_dump(obj, filename, **kwargs):
    """
    Dump an object to a json file using PMGJSONEncoder. Note that these
    objects can be lists, dicts or otherwise nested pymatgen objects that
    support the to_dict and from_dict MSONAble protocol.

    Args:
        obj (object): Object to dump.
        filename (str): Filename of file to open. Can be gzipped or bzipped.
        \*\*kwargs: Any of the keyword arguments supported by the json.load
            method.
    """
    return json.dump(obj, zopen(filename, "w"), cls=PMGJSONEncoder, **kwargs)
