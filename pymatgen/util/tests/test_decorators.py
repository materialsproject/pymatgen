#!/usr/bin/env python

"""
TODO: Modify module doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "2/24/13"


import unittest
import warnings

from pymatgen.util.decorators import deprecated, cached_class, singleton, \
    requires


class DecoratorTest(unittest.TestCase):

    def test_deprecated(self):

        def A():
            pass

        @deprecated(A)
        def B():
            pass

        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            B()
            # Verify some things
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))


    def test_cached_class(self):

        @cached_class
        class Cached(object):

            def __init__(self, a):
                self.a = a

            def __str__(self):
                return str(self.a)


        a1 = Cached(1)
        a2 = Cached(2)
        a11 = Cached(1)

        self.assertEqual(id(a1), id(a11))
        self.assertNotEqual(id(a2), id(a11))

    def test_singleton(self):

        @singleton
        class Single(object):

            def __init__(self):
                self.a = 1

        a1 = Single()
        a2 = Single()

        self.assertEqual(id(a1), id(a2))

    def test_requires(self):

        try:
            import fictitious_mod
        except ImportError:
            fictitious_mod = None

        @requires(fictitious_mod is not None, "fictitious_mod is not present.")
        def use_fictitious_mod():
            print "success"

        self.assertRaises(RuntimeError, use_fictitious_mod)

        @requires(unittest is not None, "scipy is not present.")
        def use_unittest():
            return "success"

        self.assertEqual(use_unittest(), "success")

if __name__ == "__main__":
    unittest.main()


