# coding: utf-8
"""
Objects and helper function used to store the results in a MongoDb database
"""
from __future__ import division, print_function, unicode_literals

import collections
import copy

try:
    from matgendb.dbconfig import DBConfig
except ImportError:
    pass


def scan_nestdict(d, key):
    """
    Scan a nested dict d, and return the first value associated
    to the given key. Returns None if key is not found.

    >>> d = {0: 1, 1: {"hello": {"world": {None: [1,2,3]}}}, "foo": [{"bar": 1}, {"color": "red"}]}
    >>> assert scan_nestdict(d, 1) == {"hello": {"world": {None: [1,2,3]}}}
    >>> assert scan_nestdict(d, "hello") == {"world": {None: [1,2,3]}}
    >>> assert scan_nestdict(d, "world") == {None: [1,2,3]}
    >>> assert scan_nestdict(d, None) == [1,2,3]
    >>> assert scan_nestdict(d, "color") == "red"
    """
    if isinstance(d, (list, tuple)):
        #print("got list: ", d)
        for item in d:
            res = scan_nestdict(item, key)
            if res is not None:
                return res
        return None

    if not isinstance(d, collections.Mapping):
        return None

    if key in d:
        return d[key]
    else:
        for v in d.values():
            res = scan_nestdict(v, key)
            if res is not None:
                return res
        return None


class DBConnector(object):
    #@classmethod
    #def from_file(cls, filepath):

    def __init__(self, config_dict=None):
        self.config = DBConfig(config_dict=config_dict)

    def __repr__(self):
        return "<%s object at %s>" % (self.__class__.__name__, id(self))

    def __str__(self):
        return "%s configuration:\n%s" % (self.__class__.__name__, str(self.config))

    def deepcopy(self):
        return copy.deepcopy(self)

    def set_collection_name(self, value):
        """Set the name of the collection."""
        self.config.collection = str(value)

    def get_collection(self):
        """
        Establish a connection with the database and returns the collection
        """
        from pymongo import MongoClient
        config = self.config
        #client = MongoClient(host=config.host, port=config.port)
        client = MongoClient()
        db = client[config.dbname]

        # Authenticate if username is specified.
        #if config.user:
        #    db.autenticate(config.user, password=config.password)

        return db[config.collection]


if __name__ == "__main__":
    connector = DBConnector()
    connector.set_collection_name("foo")
    print(connector)
    print(connector.connect())


