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


def mongo_getattr(rec, key):
    """
    Get value from dict using MongoDB dot-separated path semantics.
    For example:

    >>> assert mongo_getattr({'a': {'b': 1}, 'x': 2}, 'a.b') == 1
    >>> assert mongo_getattr({'a': {'b': 1}, 'x': 2}, 'x') == 2
    >>> assert mongo_getattr({'a': {'b': 1}, 'x': 2}, 'a.b.c') is None

    :param rec: mongodb document
    :param key: path to mongo value
    :param default: default to return if not found
    :return: value, potentially nested, or default if not found
    :raise: AttributeError, if record is not a dict or key is not found.
    """
    if not isinstance(rec, collections.Mapping):
        raise AttributeError('input record must act like a dict')
    if not rec:
        raise AttributeError('Empty dict')

    if not '.' in key:
        return rec.get(key)

    for key_part in key.split('.'):
        if not isinstance(rec, collections.Mapping):
            raise AttributeError('not a mapping for rec_part %s' % key_part)
        if not key_part in rec:
            raise AttributeError('key %s not in dict %s' % key)
        rec = rec[key_part]

    return rec


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



import functools
def try_or_return(exception_list, value_if_exception=None):
    """
    Decorator for functions or methods. Executes the callable in a try block,
    and returns value_if_exception if one of the exceptions listed
    is exception_list is raised. Example

    @try_or_return(ValueError)
    def return_none_if_value_error(self):

    @try_or_return((ValueError, KeyError), "hello")
    def another_method(self):
    """
    # we need a tuple of exceptions.
    if isinstance(exception_list, list): 
        exception_list = tuple(exception_list)
    elif not isinstance(exception_list, tuple): 
        exception_list = (exception_list,)

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except exception_list:
                return value_if_exception
            else:
                raise
        return wrapper
    return decorator

class DBConnector(object):

    def __init__(self, config_dict=None):
        self.config = {}
        if config_dict is not None and config_dict:
            self.config = DBConfig(config_dict=config_dict)

    def __bool__(self):
        return bool(self.config)

    __nonzero__ = __bool__

    def __repr__(self):
        return "<%s object at %s>" % (self.__class__.__name__, id(self))

    def __str__(self):
        return str(self.config)

    def deepcopy(self):
        return copy.deepcopy(self)

    def set_collection_name(self, value):
        """Set the name of the collection, return old value"""
        old = self.config.collection
        self.config.collection = str(value)
        return old

    @try_or_return(Exception, None)
    def get_collection(self, **kwargs):
        """
        Establish a connection with the database. 

        Returns MongoDb collection
        """
        from pymongo import MongoClient
        config = self.config

        # TODO
        #if config.host or config.port:
        #    client = MongoClient(host=config.host, port=config.port)
        #else:
        client = MongoClient()
        db = client[config.dbname]

        # Authenticate if username is specified.
        #if config.user:
        #    db.autenticate(config.user, password=config.password)

        return db[config.collection]


if __name__ == "__main__":
    connector = DBConnector()
    print(connector.get_collection())
    #connector.set_collection_name("foo")
    print(connector)
    print(connector.get_collection())

    import unittest

    class TryOrReturnTest(unittest.TestCase):

        def test_decorator(self):
            class A(object):
                @try_or_return(ValueError, "hello")
                def return_one(self):
                    return 1

                @try_or_return(ValueError, "hello")
                def return_hello(self):
                    raise ValueError()

                @try_or_return(KeyError, "hello")
                def reraise_value_error(self):
                    raise ValueError()

                @try_or_return([KeyError, ValueError], "hello")
                def catch_exc_list(self):
                    import random
                    if random.randint(0, 1) == 0:
                        raise ValueError()
                    else:
                        raise KeyError()

            a = A()
            assert a.return_one() == 1
            assert a.return_hello() == "hello"
            with self.assertRaises(ValueError):
                a.reraise_value_error()
            assert a.catch_exc_list() == "hello"

    unittest.main()
