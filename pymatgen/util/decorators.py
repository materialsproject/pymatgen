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

from collections import namedtuple
from functools import wraps, update_wrapper

try:
    from threading import RLock
except:
    class RLock:
        'Dummy reentrant lock for builds without threads'
        def __enter__(self): pass
        def __exit__(self, exctype, excinst, exctb): pass

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

################################################################################
### LRU Cache function decorator. Backported from python 3 functools
################################################################################

_CacheInfo = namedtuple("CacheInfo", ["hits", "misses", "maxsize", "currsize"])

class _HashedSeq(list):
    """ This class guarantees that hash() will be called no more than once
        per element.  This is important because the lru_cache() will hash
        the key multiple times on a cache miss.

    """

    __slots__ = 'hashvalue'

    def __init__(self, tup, hash=hash):
        self[:] = tup
        self.hashvalue = hash(tup)

    def __hash__(self):
        return self.hashvalue


def _make_key(args, kwds, typed,
             kwd_mark = (object(),),
             fasttypes = {int, str, frozenset, type(None)},
             sorted=sorted, tuple=tuple, type=type, len=len):
    """Make a cache key from optionally typed positional and keyword arguments

    The key is constructed in a way that is flat as possible rather than
    as a nested structure that would take more memory.

    If there is only a single argument and its data type is known to cache
    its hash value, then that argument is returned without a wrapper.  This
    saves space and improves lookup speed.

    """
    key = args
    if kwds:
        sorted_items = sorted(kwds.items())
        key += kwd_mark
        for item in sorted_items:
            key += item
    if typed:
        key += tuple(type(v) for v in args)
        if kwds:
            key += tuple(type(v) for k, v in sorted_items)
    elif len(key) == 1 and type(key[0]) in fasttypes:
        return key[0]
    return _HashedSeq(key)

def lru_cache(maxsize=128, typed=False):
    """Least-recently-used cache decorator.

    If *maxsize* is set to None, the LRU features are disabled and the cache
    can grow without bound.

    If *typed* is True, arguments of different types will be cached separately.
    For example, f(3.0) and f(3) will be treated as distinct calls with
    distinct results.

    Arguments to the cached function must be hashable.

    View the cache statistics named tuple (hits, misses, maxsize, currsize)
    with f.cache_info().  Clear the cache and statistics with f.cache_clear().
    Access the underlying function with f.__wrapped__.

    See:  http://en.wikipedia.org/wiki/Cache_algorithms#Least_Recently_Used

    """

    # Users should only access the lru_cache through its public API:
    #       cache_info, cache_clear, and f.__wrapped__
    # The internals of the lru_cache are encapsulated for thread safety and
    # to allow the implementation to change (including a possible C version).

    # Early detection of an erroneous call to @lru_cache without any arguments
    # resulting in the inner function being passed to maxsize instead of an
    # integer or None.
    if maxsize is not None and not isinstance(maxsize, int):
        raise TypeError('Expected maxsize to be an integer or None')

    # Constants shared by all lru cache instances:
    sentinel = object()          # unique object used to signal cache misses
    make_key = _make_key         # build a key from the function arguments
    PREV, NEXT, KEY, RESULT = 0, 1, 2, 3   # names for the link fields

    def decorating_function(user_function):
        cache = {}
        hits = [0]
        misses = [0]
        full = [False]
        cache_get = cache.get    # bound method to lookup a key or return None
        lock = RLock()           # because linkedlist updates aren't threadsafe
        root = []                # root of the circular doubly linked list
        root[:] = [root, root, None, None]     # initialize by pointing to self
        r = [root]

        if maxsize == 0:

            def wrapper(*args, **kwds):
                # No caching -- just a statistics update after a successful call
                result = user_function(*args, **kwds)
                misses[0] += 1
                return result

        elif maxsize is None:

            def wrapper(*args, **kwds):
                # Simple caching without ordering or size limit
                key = make_key(args, kwds, typed)
                result = cache_get(key, sentinel)
                if result is not sentinel:
                    hits[0] += 1
                    return result
                result = user_function(*args, **kwds)
                cache[key] = result
                misses[0] += 1
                return result

        else:
            def wrapper(*args, **kwds):
                # Size limited caching that tracks accesses by recency
                key = make_key(args, kwds, typed)
                with lock:
                    link = cache_get(key)
                    if link is not None:
                        # Move the link to the front of the circular queue
                        link_prev, link_next, _key, result = link
                        link_prev[NEXT] = link_next
                        link_next[PREV] = link_prev
                        last = r[0][PREV]
                        last[NEXT] = r[0][PREV] = link
                        link[PREV] = last
                        link[NEXT] = r[0]
                        hits[0] += 1
                        return result
                result = user_function(*args, **kwds)
                with lock:
                    if key in cache:
                        # Getting here means that this same key was added to the
                        # cache while the lock was released.  Since the link
                        # update is already done, we need only return the
                        # computed result and update the count of misses.
                        pass
                    elif full[0]:
                        # Use the old root to store the new key and result.
                        oldroot = r[0]
                        oldroot[KEY] = key
                        oldroot[RESULT] = result
                        # Empty the oldest link and make it the new root.
                        # Keep a reference to the old key and old result to
                        # prevent their ref counts from going to zero during the
                        # update. That will prevent potentially arbitrary object
                        # clean-up code (i.e. __del__) from running while we're
                        # still adjusting the links.
                        r[0] = oldroot[NEXT]
                        oldkey = r[0][KEY]
                        oldresult = r[0][RESULT]
                        r[0][KEY] = r[0][RESULT] = None
                        # Now update the cache dictionary.
                        del cache[oldkey]
                        # Save the potentially reentrant cache[key] assignment
                        # for last, after the root and links have been put in
                        # a consistent state.
                        cache[key] = oldroot
                    else:
                        # Put result in a new link at the front of the queue.
                        last = r[0][PREV]
                        link = [last, r[0], key, result]
                        last[NEXT] = r[0][PREV] = cache[key] = link
                        full[0] = (len(cache) >= maxsize)
                    misses[0] += 1
                return result

        def cache_info():
            """Report cache statistics"""
            with lock:
                return _CacheInfo(hits[0], misses[0], maxsize, len(cache))

        def cache_clear():
            """Clear the cache and cache statistics"""
            with lock:
                cache.clear()
                root[:] = [root, root, None, None]
                r[0] = root
                hits[0] = 0
                misses[0] = 0
                full[0] = False

        wrapper.cache_info = cache_info
        wrapper.cache_clear = cache_clear
        return update_wrapper(wrapper, user_function)

    return decorating_function

