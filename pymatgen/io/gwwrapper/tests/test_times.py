#!/usr/bin/env python

__author__ = 'setten'

import timeit

n = 10000


def test_write():
    """Stupid test function"""
    l = []
    for i in range(100000):
        l.append(i)
    f = open('test', 'w')
    f.write(str(l))
    f.close()

if __name__ == '__main__':
    n = 100
    print(timeit.timeit("test_write()", setup="from __main__ import test_write", number=n))
