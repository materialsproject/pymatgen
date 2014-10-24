# coding: utf-8
"""
This module implements the pickler objects used in abinitio.
"""

from __future__ import unicode_literals, division, print_function

import pickle

from pymatgen.core.periodic_table import Element


class PmgPickler(pickle.Pickler):
    """
    Persistence of External Objects as described in section 12.1.5.1 of 
    https://docs.python.org/3/library/pickle.html
    """
    def persistent_id(self, obj):
        """Instead of pickling as a regular class instance, we emit a persistent ID."""
        if isinstance(obj, Element):
            # Here, our persistent ID is simply a tuple, containing a tag and a key
            return obj.__class__.__name__, obj._symbol
        else:
            # If obj does not have a persistent ID, return None. This means obj needs to be pickled as usual.
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
        type_tag, key_id = pid

        if type_tag == "Element":
            return Element(key_id)
        else:
            # Always raises an error if you cannot return the correct object.
            # Otherwise, the unpickler will think None is the object referenced by the persistent ID.
            raise pickle.UnpicklingError("unsupported persistent object with pid %s" % pid)


def pmg_pickle_load(filobj, **kwargs):
    """
    Loads a pickle file and deserialize it with PmgUnpickler.

    Args:
        filobj: File-like object
        \*\*kwargs: Any of the keyword arguments supported by PmgUnpickler

    Returns:
        Deserialized object. 
    """
#    return PmgUnpickler(filobj, **kwargs).load()
    return pickle.load(filobj, **kwargs)


def pmg_pickle_dump(obj, filobj, **kwargs):
    """
    Dump an object to a pickle file using PmgPickler.

    Args:
        obj (object): Object to dump.
        fileobj: File-like object
        \*\*kwargs: Any of the keyword arguments supported by PmgPickler
    """
    PmgPickler(filobj, **kwargs).dump(obj)
    return pickle.dump(obj, filobj, **kwargs)
