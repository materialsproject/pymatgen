# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

import json
import functools
import pickle

from pymatgen.core.periodic_table import Element

"""
Most features of this module has been moved to monty. Please refer to
monty.json and monty.serialization documentation.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 30, 2012"


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
    with open(filename, "wt") as fh:
        json.dump(obj, fh, indent=4, sort_keys=4)


class PmgPickler(pickle.Pickler):
    """
    Persistence of External Objects as described in section 12.1.5.1 of
    https://docs.python.org/3/library/pickle.html
    """

    def persistent_id(self, obj):
        """Instead of pickling as a regular class instance, we emit a
        persistent ID."""
        if isinstance(obj, Element):
            # Here, our persistent ID is simply a tuple, containing a tag and
            # a key
            return obj.__class__.__name__, obj.symbol
        else:
            # If obj does not have a persistent ID, return None. This means obj
            # needs to be pickled as usual.
            return None


class PmgUnpickler(pickle.Unpickler):
    """
    Persistence of External Objects as described in section 12.1.5.1 of
    https://docs.python.org/3/library/pickle.html
    """

    def persistent_load(self, pid):
        """
        This method is invoked whenever a persistent ID is encountered.
        Here, pid is the tuple returned by PmgPickler.
        """
        try:
            type_tag, key_id = pid
        except Exception as exc:
            # Sometimes we get a string such as ('Element', u'C') instead
            # of a real tuple. Use ast to evalute the expression (much safer
            # than eval).
            import ast
            type_tag, key_id = ast.literal_eval(pid)

        if type_tag == "Element":
            return Element(key_id)
        else:
            # Always raises an error if you cannot return the correct object.
            # Otherwise, the unpickler will think None is the object referenced
            # by the persistent ID.
            raise pickle.UnpicklingError(
                "unsupported persistent object with pid %s" % pid)


def pmg_pickle_load(filobj, **kwargs):
    """
    Loads a pickle file and deserialize it with PmgUnpickler.

    Args:
        filobj: File-like object
        \\*\\*kwargs: Any of the keyword arguments supported by PmgUnpickler

    Returns:
        Deserialized object.
    """
    return PmgUnpickler(filobj, **kwargs).load()


def pmg_pickle_dump(obj, filobj, **kwargs):
    """
    Dump an object to a pickle file using PmgPickler.

    Args:
        obj (object): Object to dump.
        fileobj: File-like object
        \\*\\*kwargs: Any of the keyword arguments supported by PmgPickler
    """
    return PmgPickler(filobj, **kwargs).dump(obj)


class SlotPickleMixin(object):
    """
    This mixin makes it possible to pickle/unpickle objects with __slots__
    defined.
    """

    def __getstate__(self):
        return dict(
            (slot, getattr(self, slot))
            for slot in self.__slots__ if hasattr(self, slot)
        )

    def __setstate__(self, state):
        for slot, value in state.items():
            setattr(self, slot, value)
