#!/usr/bin/env python

'''
This module contains useful decorators for a variety of functions.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Dec 31, 2011"

import logging
import datetime
from functools import wraps

def singleton(cls):
    """
    This decorator can be used to create a singleton out of a class.
    """
    
    instances = {}
    
    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance

def cached_class(klass):
    """
    Decorator to cache class instances by constructor arguments.
    This results in a class that behaves like a singleton for each
    set of constructor arguments, ensuring efficiency. 
    
    Note that this should be used for *immutable classes only*.  Having 
    a cached mutable class makes very little sense.  For efficiency,
    avoid using this decorator for situations where there are many 
    constructor arguments permutations.
      
    The keywords argument dictionary is converted to a tuple because
    dicts are mutable; keywords themselves are strings and 
    so are always hashable, but if any arguments (keyword
    or positional) are non-hashable, that set of arguments
    is not cached.
    """
    cache = {}
    
    @wraps(klass, assigned=('__name__', '__module__'), updated=())
    class _decorated(klass):
        # The wraps decorator can't do this because __doc__
        # isn't writable once the class is created
        __doc__ = klass.__doc__
        def __new__(cls, *args, **kwds):
            key = (cls,) + args + tuple(kwds.iteritems())
            if key not in cache:
                o = super(klass, cls).__new__(cls, *args, **kwds)
                cache[key] = o
            return cache[key]

    
    return _decorated

def logged(level = logging.DEBUG):
    def wrap(f):
        logger= logging.getLogger("{}.{}".format(f.__module__, f.__name__))
        def wrapped_f(*args, **kwargs):
            
            logger.log(level, "Called at {} with args = {} and kwargs = {}".format(datetime.datetime.now(), args, kwargs))
            data = f(*args, **kwargs)
            logger.log(level, "Done at {} with args = {} and kwargs = {}".format(datetime.datetime.now(), args, kwargs))
            return data
        return wrapped_f
    return wrap