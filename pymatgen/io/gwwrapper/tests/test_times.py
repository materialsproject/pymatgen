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


def test_folders():
    for n in range(0, 100, 1):
        os.mkdir('test_folder')
        os.chdir('test_folder')
        os.chdir('..')
        os.rmdir('test_folder')


if __name__ == '__main__':
    n = 100
    print 'write', timeit.timeit("test_write()", setup="from __main__ import test_write", number=n)
    print 'read', timeit.timeit("test_read()", setup="from __main__ import test_read", number=n)
    print 'folders', timeit.timeit("test_folders()", setup="from __main__ import test_folders", number=n)
