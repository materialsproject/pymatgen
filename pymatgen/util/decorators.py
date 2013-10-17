#!/usr/bin/env python

"""
This module contains useful decorators for a variety of functions.
"""
from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Dec 31, 2011"

import logging
import datetime
import warnings

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

    @wraps(klass, assigned=("__name__", "__module__"), updated=())
    class _decorated(klass):
        # The wraps decorator can't do this because __doc__
        # isn't writable once the class is created
        __doc__ = klass.__doc__

        def __new__(cls, *args, **kwds):
            #key = (cls,) + args + tuple(kwds.items())
            #if key not in cache:
            #    o = super(klass, cls).__new__(cls, *args, **kwds)
            #    cache[key] = o
            #return cache[key]
            key = (cls,) + args + tuple(kwds.iteritems())
            try:
                inst = cache.get(key, None)
            except TypeError:
                # Can't cache this set of arguments
                inst = key = None
            if inst is None:
                # Technically this is cheating, but it works,
                # and takes care of initializing the instance
                # (so we can override __init__ below safely);
                # calling up to klass.__new__ would be the
                # "official" way to create the instance, but
                # that raises DeprecationWarning if there are
                # args or kwds and klass does not override
                # __new__ (which most classes don't), because
                # object.__new__ takes no parameters (and in
                # Python 3 the warning will become an error)
                inst = klass(*args, **kwds)
                # This makes isinstance and issubclass work
                # properly
                inst.__class__ = cls
                if key is not None:
                    cache[key] = inst
            return inst

        def __init__(self, *args, **kwds):
            # This will be called every time __new__ is
            # called, so we skip initializing here and do
            # it only when the instance is created above
            pass

    return _decorated


def logged(level=logging.DEBUG):
    """
    Useful logging decorator. If a method is logged, the beginning and end of
    the method call will be logged at a pre-specified level.

    Args:
        level:
            Level to log method at. Defaults to DEBUG.
    """
    def wrap(f):
        logger = logging.getLogger("{}.{}".format(f.__module__, f.__name__))

        def wrapped_f(*args, **kwargs):
            logger.log(level, "Called at {} with args = {} and kwargs = {}"
                       .format(datetime.datetime.now(), args, kwargs))
            data = f(*args, **kwargs)
            logger.log(level, "Done at {} with args = {} and kwargs = {}"
                       .format(datetime.datetime.now(), args, kwargs))
            return data

        return wrapped_f
    return wrap


def deprecated(replacement=None):
    """
    Decorator to mark classes or functions as deprecated,
    with a possible replacement.

    Args:
        replacement:
            A replacement class or method.

    Returns:
        Original function, but with a warning to use the updated class.
    """
    def wrap(old):
        def wrapped(*args, **kwargs):
            msg = "{} is deprecated".format(old.__name__)
            if replacement is not None:
                msg += "; use {} in {} instead.".format(
                    replacement.__name__, replacement.__module__)
            warnings.simplefilter('default')
            warnings.warn(msg, DeprecationWarning, stacklevel=2)
            return old(*args, **kwargs)
        return wrapped
    return wrap


class requires(object):
    """
    Decorator to mark classes or functions as requiring a specified condition
    to be true. This can be used to present useful error messages for
    optional dependencies. For example, decorating the following code will
    check if scipy is present and if not, a runtime error will be raised if
    someone attempts to call the use_scipy function::

        try:
            import scipy
        except ImportError:
            scipy = None

        @requires(scipy is not None, "scipy is not present.")
        def use_scipy():
            print scipy.majver

    Args:
        condition:
            Condition necessary to use the class or function.
        message:
            A message to be displayed if the condition is not True.
    """

    def __init__(self, condition, message):
        self.condition = condition
        self.message = message

    def __call__(self, callable):
        @wraps(callable)
        def decorated(*args, **kwargs):
            if not self.condition:
                raise RuntimeError(self.message)
            return callable(*args, **kwargs)
        return decorated


def enable_logging(main):
    """
    This decorator is used to decorate main functions.
    It adds the initialization of the logger and an argument parser that allows one to select the loglevel.
    Useful if we are writing simple main functions that call libraries where the logging module is used

    Args:
        main:
            main function.
    """
    @wraps(main)
    def wrapper(*args, **kwargs):
        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('--loglevel', default="ERROR", type=str,
                            help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

        options = parser.parse_args()

        # loglevel is bound to the string value obtained from the command line argument. 
        # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
        numeric_level = getattr(logging, options.loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % options.loglevel)
        logging.basicConfig(level=numeric_level)

        retcode = main(*args, **kwargs)
        return retcode

    return wrapper
