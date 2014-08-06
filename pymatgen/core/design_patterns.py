"""
This module defines some useful design patterns.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"


class Enum(set):
    """
    Creates an enum out of a set.
    """
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError


class AttrDict(dict):
    """
    Allows to access dict keys as obj.foo in addition to the traditional way
    obj['foo']"
    """
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

    def copy(self):
        newd = super(AttrDict, self).copy()
        return self.__class__(**newd)


class NullFile(object):
    """A file object that is associated to /dev/null."""
    def __init__(self):
        import os
        return open(os.devnull, 'w')


class NullStream(object):
    """A fake stream with a no-op write.."""
    def write(*args):
        """no-op"""


class FrozenDict(dict):
    """A dictionary that does not permit to redefine its keys."""
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, val):
        if key in self:
            raise KeyError("Cannot overwrite existent key: %s" % str(key))

        dict.__setitem__(self, key, val)

    def update(self, *args, **kwargs):
        for k, v in dict(*args, **kwargs).items():
            self[k] = v
