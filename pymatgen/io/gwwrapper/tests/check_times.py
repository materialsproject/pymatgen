# coding: utf-8

from __future__ import unicode_literals, print_function, division

#!/usr/bin/env python

__author__ = 'setten'

import timeit
import os


def test_write(m):
    """
    file write test function
    """
    l = []
    for i in range(m):
        l.append(i)
    f = open('test', 'w')
    f.write(str(l))
    f.close()


def test_read():
    """
    file read test function
    """
    f = open('test', mode='r')
    line = f.read()
    f.close()


def test_make_folders():
    """
    creating and removing test function
    """
    for i in range(0, 100, 1):
        os.mkdir('test_folder_%s' % i)
        os.rmdir('test_folder_%s' % i)


def test_cd():
    """
    filesystems navigation test function
    """
    for i in range(0, 100, 1):
        os.mkdir('test_folder_%s' % i)
        os.chdir('test_folder_%s' % i)
        os.chdir('..')
        os.rmdir('test_folder_%s' % i)


if __name__ == '__main__':
    """
    Testing file writing and reading and folder operations, an asserion error is raised when the filesystem seems to be
    too slow
    """

    n = 100

    mk_folders = timeit.timeit("test_make_folders()", setup="from __main__ import test_make_folders", number=n)
    cd_folders = timeit.timeit("test_cd()", setup="from __main__ import test_cd", number=n)

    print 'made and removed 100 folders in', mk_folders, 's'
    print 'made, moved to, returned and removed 100 folders in', cd_folders, 's'

    for my_m in [50000, 100000, 200000]:
        test_write(my_m)
        size = os.path.getsize('test')
        print 'file size', size

        setup = "from __main__ import test_write\nmy_m =" + str(my_m)

        write = timeit.timeit("test_write(my_m)", setup=setup, number=n)
        read = timeit.timeit("test_read()", setup="from __main__ import test_read", number=n*200)

        print n, 'times written in ', write, 's, ', size * n / write / 1000000, 'MB/s'
        print n*200, 'times read in  ', read, 's, ', size * 200 * n / read / 1000000, 'MB/s'

#        assert write < 10 * n / 100 * my_m / 100000
#        assert read < 10 * n / 100 * my_m / 100000

#    assert mk_folders < 1 * n / 100
#    assert cd_folders < 1 * n / 100
