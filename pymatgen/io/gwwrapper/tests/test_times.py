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
    print 'write', timeit.timeit("test_write()", setup="from __main__ import test_write", number=n)
    print 'read', timeit.timeit("test_read()", setup="from __main__ import test_read", number=n)
    print 'mk folders', timeit.timeit("test_make_folders()", setup="from __main__ import test_make_folders", number=n)
    print 'cd folders', timeit.timeit("test_cd()", setup="from __main__ import test_cd", number=n)
