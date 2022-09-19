# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Most features of this module has been moved to monty. Please refer to
monty.json and monty.serialization documentation.
"""

import functools
import json
import pickle

from pymatgen.core.periodic_table import Element


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
        d["@module"] = type(self).__module__
        d["@class"] = type(self).__name__
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
        except Exception:
            # Sometimes we get a string such as ('Element', u'C') instead
            # of a real tuple. Use ast to evaluate the expression (much safer
            # than eval).
            import ast

            type_tag, key_id = ast.literal_eval(pid)

        if type_tag == "Element":
            return Element(key_id)
        # Always raises an error if you cannot return the correct object.
        # Otherwise, the unpickler will think None is the object referenced
        # by the persistent ID.
        raise pickle.UnpicklingError(f"unsupported persistent object with pid {pid}")


def pmg_pickle_load(filobj, **kwargs):
    """
    Loads a pickle file and deserialize it with PmgUnpickler.

    Args:
        filobj: File-like object
        **kwargs: Any of the keyword arguments supported by PmgUnpickler

    Returns:
        Deserialized object.
    """
    return PmgUnpickler(filobj, **kwargs).load()


def pmg_pickle_dump(obj, filobj, **kwargs):
    """
    Dump an object to a pickle file using PmgPickler.

    Args:
        obj : Object to dump.
        fileobj: File-like object
        **kwargs: Any of the keyword arguments supported by PmgPickler
    """
    return PmgPickler(filobj, **kwargs).dump(obj)
