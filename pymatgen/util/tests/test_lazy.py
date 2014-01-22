"""
Test suite for lazy
Taken from https://pypi.python.org/pypi/lazy v1.1
"""
import sys
import unittest

from pymatgen.util.lazy import lazy


class TestCase(unittest.TestCase):

    def assertException(self, exc_cls, pattern, func, *args, **kw):
        """Assert an exception of type 'exc_cls' is raised and
        'pattern' is contained in the exception message.
        """
        try:
            func(*args, **kw)
        except exc_cls, e:
            exc_str = str(e)
        else:
            self.fail('%s not raised' % (exc_cls.__name__,))

        if pattern not in exc_str:
            self.fail('%r not in %r' % (pattern, exc_str))


class LazyTests(TestCase):

    def test_evaluate(self):
        # Lazy attributes should be evaluated when accessed.
        called = []

        class Foo(object):
            @lazy
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertEqual(f.foo, 1)
        self.assertEqual(len(called), 1)

    def test_evaluate_once(self):
        # Lazy attributes should be evaluated only once.
        called = []

        class Foo(object):
            @lazy
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertEqual(f.foo, 1)
        self.assertEqual(f.foo, 1)
        self.assertEqual(f.foo, 1)
        self.assertEqual(len(called), 1)

    def test_private_attribute(self):
        # It should be possible to create private, name-mangled
        # lazy attributes.
        called = []

        class Foo(object):
            @lazy
            def __foo(self):
                called.append('foo')
                return 1
            def get_foo(self):
                return self.__foo

        f = Foo()
        self.assertEqual(f.get_foo(), 1)
        self.assertEqual(f.get_foo(), 1)
        self.assertEqual(f.get_foo(), 1)
        self.assertEqual(len(called), 1)

    def test_reserved_attribute(self):
        # It should be possible to create reserved lazy attributes.
        called = []

        class Foo(object):
            @lazy
            def __foo__(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertEqual(f.__foo__, 1)
        self.assertEqual(f.__foo__, 1)
        self.assertEqual(f.__foo__, 1)
        self.assertEqual(len(called), 1)

    def test_result_shadows_descriptor(self):
        # The result of the function call should be stored in
        # the object __dict__, shadowing the descriptor.
        called = []

        class Foo(object):
            @lazy
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertTrue(isinstance(Foo.foo, lazy))
        self.assertTrue(f.foo is f.foo)
        self.assertTrue(f.foo is f.__dict__['foo']) # !
        self.assertEqual(len(called), 1)

        self.assertEqual(f.foo, 1)
        self.assertEqual(f.foo, 1)
        self.assertEqual(len(called), 1)

        lazy.invalidate(f, 'foo')

        self.assertEqual(f.foo, 1)
        self.assertEqual(len(called), 2)

        self.assertEqual(f.foo, 1)
        self.assertEqual(f.foo, 1)
        self.assertEqual(len(called), 2)

    def test_readonly_object(self):
        # The descriptor should raise an AttributeError when lazy is
        # used on a read-only object (an object with __slots__).
        called = []

        class Foo(object):
            __slots__ = ()
            @lazy
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertEqual(len(called), 0)

        self.assertException(AttributeError,
            "'Foo' object has no attribute '__dict__'",
            getattr, f, 'foo')

        # The function was not called
        self.assertEqual(len(called), 0)

    def test_introspection(self):
        # The lazy decorator should support basic introspection.

        class Foo(object):
            def foo(self):
                """foo func doc"""
            @lazy
            def bar(self):
                """bar func doc"""

        self.assertEqual(Foo.foo.__name__, "foo")
        self.assertEqual(Foo.foo.__doc__, "foo func doc")
        self.assertEqual(Foo.foo.__module__, "pymatgen.util.tests.test_lazy")

        self.assertEqual(Foo.bar.__name__, "bar")
        self.assertEqual(Foo.bar.__doc__, "bar func doc")
        self.assertEqual(Foo.bar.__module__, "pymatgen.util.tests.test_lazy")


class InvalidateTests(TestCase):

    def test_invalidate_attribute(self):
        # It should be possible to invalidate a lazy attribute.
        called = []

        class Foo(object):
            @lazy
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertEqual(f.foo, 1)
        self.assertEqual(len(called), 1)

        lazy.invalidate(f, 'foo')

        self.assertEqual(f.foo, 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_attribute_twice(self):
        # It should be possible to invalidate a lazy attribute
        # twice without causing harm.
        called = []

        class Foo(object):
            @lazy
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertEqual(f.foo, 1)
        self.assertEqual(len(called), 1)

        lazy.invalidate(f, 'foo')
        lazy.invalidate(f, 'foo') # Nothing happens

        self.assertEqual(f.foo, 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_uncalled_attribute(self):
        # It should be possible to invalidate an empty attribute
        # cache without causing harm.
        called = []

        class Foo(object):
            @lazy
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertEqual(len(called), 0)
        lazy.invalidate(f, 'foo') # Nothing happens

    def test_invalidate_private_attribute(self):
        # It should be possible to invalidate a private lazy attribute.
        called = []

        class Foo(object):
            @lazy
            def __foo(self):
                called.append('foo')
                return 1
            def get_foo(self):
                return self.__foo

        f = Foo()
        self.assertEqual(f.get_foo(), 1)
        self.assertEqual(len(called), 1)

        lazy.invalidate(f, '__foo')

        self.assertEqual(f.get_foo(), 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_mangled_attribute(self):
        # It should be possible to invalidate a private lazy attribute
        # by its mangled name.
        called = []

        class Foo(object):
            @lazy
            def __foo(self):
                called.append('foo')
                return 1
            def get_foo(self):
                return self.__foo

        f = Foo()
        self.assertEqual(f.get_foo(), 1)
        self.assertEqual(len(called), 1)

        lazy.invalidate(f, '_Foo__foo')

        self.assertEqual(f.get_foo(), 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_reserved_attribute(self):
        # It should be possible to invalidate a reserved lazy attribute.
        called = []

        class Foo(object):
            @lazy
            def __foo__(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertEqual(f.__foo__, 1)
        self.assertEqual(len(called), 1)

        lazy.invalidate(f, '__foo__')

        self.assertEqual(f.__foo__, 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_nonlazy_attribute(self):
        # Invalidating an attribute that is not lazy should
        # raise an AttributeError.
        called = []

        class Foo(object):
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertException(AttributeError,
            "'Foo.foo' is not a lazy attribute",
            lazy.invalidate, f, 'foo')

    def test_invalidate_nonlazy_private_attribute(self):
        # Invalidating a private attribute that is not lazy should
        # raise an AttributeError.
        called = []

        class Foo(object):
            def __foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertException(AttributeError,
            "'Foo._Foo__foo' is not a lazy attribute",
            lazy.invalidate, f, '__foo')

    def test_invalidate_unknown_attribute(self):
        # Invalidating an unknown attribute should
        # raise an AttributeError.
        called = []

        class Foo(object):
            @lazy
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertException(AttributeError,
            "type object 'Foo' has no attribute 'bar'",
            lazy.invalidate, f, 'bar')

    def test_invalidate_readonly_object(self):
        # Calling invalidate on a read-only object should
        # raise an AttributeError.
        called = []

        class Foo(object):
            __slots__ = ()
            @lazy
            def foo(self):
                called.append('foo')
                return 1

        f = Foo()
        self.assertException(AttributeError,
            "'Foo' object has no attribute '__dict__'",
            lazy.invalidate, f, 'foo')


# A lazy subclass
class cached(lazy):
    pass


class InvalidateSubclassTests(TestCase):

    def test_invalidate_attribute(self):
        # It should be possible to invalidate a cached attribute.
        called = []

        class Bar(object):
            @cached
            def bar(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertEqual(b.bar, 1)
        self.assertEqual(len(called), 1)

        cached.invalidate(b, 'bar')

        self.assertEqual(b.bar, 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_attribute_twice(self):
        # It should be possible to invalidate a cached attribute
        # twice without causing harm.
        called = []

        class Bar(object):
            @cached
            def bar(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertEqual(b.bar, 1)
        self.assertEqual(len(called), 1)

        cached.invalidate(b, 'bar')
        cached.invalidate(b, 'bar') # Nothing happens

        self.assertEqual(b.bar, 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_uncalled_attribute(self):
        # It should be possible to invalidate an empty attribute
        # cache without causing harm.
        called = []

        class Bar(object):
            @cached
            def bar(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertEqual(len(called), 0)
        cached.invalidate(b, 'bar') # Nothing happens

    def test_invalidate_private_attribute(self):
        # It should be possible to invalidate a private cached attribute.
        called = []

        class Bar(object):
            @cached
            def __bar(self):
                called.append('bar')
                return 1
            def get_bar(self):
                return self.__bar

        b = Bar()
        self.assertEqual(b.get_bar(), 1)
        self.assertEqual(len(called), 1)

        cached.invalidate(b, '__bar')

        self.assertEqual(b.get_bar(), 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_mangled_attribute(self):
        # It should be possible to invalidate a private cached attribute
        # by its mangled name.
        called = []

        class Bar(object):
            @cached
            def __bar(self):
                called.append('bar')
                return 1
            def get_bar(self):
                return self.__bar

        b = Bar()
        self.assertEqual(b.get_bar(), 1)
        self.assertEqual(len(called), 1)

        cached.invalidate(b, '_Bar__bar')

        self.assertEqual(b.get_bar(), 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_reserved_attribute(self):
        # It should be possible to invalidate a reserved cached attribute.
        called = []

        class Bar(object):
            @cached
            def __bar__(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertEqual(b.__bar__, 1)
        self.assertEqual(len(called), 1)

        cached.invalidate(b, '__bar__')

        self.assertEqual(b.__bar__, 1)
        self.assertEqual(len(called), 2)

    def test_invalidate_uncached_attribute(self):
        # Invalidating an attribute that is not cached should
        # raise an AttributeError.
        called = []

        class Bar(object):
            def bar(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertException(AttributeError,
            "'Bar.bar' is not a cached attribute",
            cached.invalidate, b, 'bar')

    def test_invalidate_uncached_private_attribute(self):
        # Invalidating a private attribute that is not cached should
        # raise an AttributeError.
        called = []

        class Bar(object):
            def __bar(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertException(AttributeError,
            "'Bar._Bar__bar' is not a cached attribute",
            cached.invalidate, b, '__bar')

    def test_invalidate_unknown_attribute(self):
        # Invalidating an unknown attribute should
        # raise an AttributeError.
        called = []

        class Bar(object):
            @cached
            def bar(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertException(AttributeError,
            "type object 'Bar' has no attribute 'baz'",
            lazy.invalidate, b, 'baz')

    def test_invalidate_readonly_object(self):
        # Calling invalidate on a read-only object should
        # raise an AttributeError.
        called = []

        class Bar(object):
            __slots__ = ()
            @cached
            def bar(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertException(AttributeError,
            "'Bar' object has no attribute '__dict__'",
            cached.invalidate, b, 'bar')

    def test_invalidate_superclass_attribute(self):
        # cached.invalidate CANNOT invalidate a superclass (lazy) attribute.
        called = []

        class Bar(object):
            @lazy
            def bar(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertException(AttributeError,
            "'Bar.bar' is not a cached attribute",
            cached.invalidate, b, 'bar')

    def test_invalidate_subclass_attribute(self):
        # Whereas lazy.invalidate CAN invalidate a subclass (cached) attribute.
        called = []

        class Bar(object):
            @cached
            def bar(self):
                called.append('bar')
                return 1

        b = Bar()
        self.assertEqual(b.bar, 1)
        self.assertEqual(len(called), 1)

        lazy.invalidate(b, 'bar')

        self.assertEqual(b.bar, 1)
        self.assertEqual(len(called), 2)


class AssertExceptionTests(TestCase):

    def test_assert_AttributeError(self):
        self.assertException(AttributeError,
            "'AssertExceptionTests' object has no attribute 'foobar'",
            getattr, self, 'foobar')

    def test_assert_IOError(self):
        self.assertException(IOError,
            "No such file or directory",
            open, './foo/bar/baz/peng/quux', 'rb')

    def test_assert_SystemExit(self):
        self.assertException(SystemExit,
            "",
            sys.exit)

    def test_assert_exception_not_raised(self):
        self.assertRaises(AssertionError,
            self.assertException, AttributeError,
            "'AssertExceptionTests' object has no attribute 'run'",
            getattr, self, 'run')

    def test_assert_pattern_mismatch(self):
        self.assertRaises(AssertionError,
            self.assertException, AttributeError,
            "baz",
            getattr, self, 'foobar')

