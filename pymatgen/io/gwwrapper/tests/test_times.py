#!/usr/bin/env python

__author__ = 'setten'

import timeit
import os


def test_write():
    """file write test function"""
    l = []
    for i in range(100000):
        l.append(i)
    f = open('test', 'w')
    f.write(str(l))
    f.close()


def test_read():
    f = open('test', mode='r')
    line = f.read()
    f.close()


def test_make_folders():
    for i in range(0, 100, 1):
        os.mkdir('test_folder_%s' % i)
        os.rmdir('test_folder_%s' % i)


def test_cd():
    for i in range(0, 100, 1):
        os.mkdir('test_folder_%s' % i)
        os.chdir('test_folder_%s' % i)
        os.chdir('..')
        os.rmdir('test_folder_%s' % i)


if __name__ == '__main__':
    n = 100

    test_write()
    size = os.path.getsize('test')
    print 'file size', size

    write = timeit.timeit("test_write()", setup="from __main__ import test_write", number=n)
    read = timeit.timeit("test_read()", setup="from __main__ import test_read", number=n*200)
    mk_folders = timeit.timeit("test_make_folders()", setup="from __main__ import test_make_folders", number=n)
    cd_folders = timeit.timeit("test_cd()", setup="from __main__ import test_cd", number=n)

    print n, 'times written: ', write, ', ', size * n / write, 'b/s'
    print n*200, 'times read:  ', read, ', ', size * 200 * n / write, 'b/s'
    print 'made and removed 100 folders', mk_folders
    print 'made, moved to, returned and removed 100 folders', cd_folders

    assert write < 5
    assert read < 5
    assert mk_folders < 1
    assert cd_folders < 1
