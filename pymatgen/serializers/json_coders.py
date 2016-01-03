# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Most features of this module has been moved to monty. Please refer to
monty.json's documentation.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 30, 2012"

import json
import functools

from monty.serialization import loadfn, dumpfn
from monty.dev import deprecated


def pmg_serialize(method):
    """
    Decorator for methods that add MSON serializations keys 
    to the dictionary. See documentation of MSON for more details
    """
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        self = args[0]
        d = method(*args, **kwargs)
        # Add @module and @class
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    return wrapper


def json_pretty_dump(obj, filename):
    """
    Serialize obj as a JSON formatted stream to the given filename (
    pretty printing version)
    """
    with open(filename, "w") as fh:
        json.dump(obj, fh, indent=4, sort_keys=4)


@deprecated(loadfn, "Will be removed in pmg 4.0.")
def pmg_load(filename, **kwargs):
    """
    Loads a json file and deserialize it with MontyDecoder.

    Args:
        filename (str): Filename of file to open. Can be gzipped or bzipped.
        \*\*kwargs: Any of the keyword arguments supported by the json.load
            method.

    Returns:
        Deserialized pymatgen object. Note that these objects can be lists,
        dicts or otherwise nested pymatgen objects that support the as_dict()
        and from_dict MSONable protocol.
    """
    return loadfn(filename, **kwargs)


@deprecated(dumpfn, "Will be removed in pmg 4.0.")
def pmg_dump(obj, filename, **kwargs):
    """
    Dump an object to a json file using MontyEncoder. Note that these
    objects can be lists, dicts or otherwise nested pymatgen objects that
    support the as_dict() and from_dict MSONable protocol.

    Args:
        obj (object): Object to dump.
        filename (str): Filename of file to open. Can be gzipped or bzipped.
        \*\*kwargs: Any of the keyword arguments supported by the json.dump
            method.
    """
    return dumpfn(obj, filename, **kwargs)
