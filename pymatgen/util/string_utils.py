#!/usr/bin/env python

"""
This module provides utility classes for string operations.
"""

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

def str_delimited(results, header=None, delimiter="\t"):
    """
    Given a tuple of tuples, generate a delimited string form.
    >>> results = [['a','b','c'],['d','e','f'],[1,2,3]]
    >>> print str_delimited(results,delimiter=',')
    a,b,c
    d,e,f
    1,2,3
    
    Args:
        result: 2d sequence of arbitrary types.
        header: optional header
        
    Returns:
        Aligned string output in a table-like format.
    
    """
    returnstr = ''
    if header != None:
        returnstr += delimiter.join(header) + "\n"
    return returnstr + "\n".join([delimiter.join([str(m) for m in result]) for result in results])
        
def str_aligned(results, header=None):
    """
    Given a tuple, generate a nicely aligned string form.
    >>> results = [['a','b','cz'],['d','ez','f'],[1,2,3]]
    >>> print str_aligned(results)
    a    b   cz
    d   ez    f
    1    2    3
    
    Args:
        result: 2d sequence of arbitrary types.
        header: optional header
        
    Returns:
        Aligned string output in a table-like format.
    """
    k = list(zip(*results))
    stringlengths = list()
    count = 0
    for i in k:
        colMaxLength = max([len(str(m)) for m in i])
        if header != None:
            colMaxLength = max([len(str(header[count])), colMaxLength])
        stringlengths.append(colMaxLength)
        count += 1
    formatString = "   ".join(["%" + str(d) + "s" for d in stringlengths])
    returnstr = ''
    if header != None:
        header_str  = formatString % tuple(header)
        returnstr += header_str + "\n"
        returnstr += "-" * len(header_str) + "\n"
        
    return returnstr + "\n".join([formatString % tuple(result) for result in results])

def formula_double_format(afloat, ignore_ones = True, tol=1e-8):
    """
    This function is used to make pretty formulas by formatting the amounts.
    Instead of Li1.0 Fe1.0 P1.0 O4.0, you get LiFePO4.
    
    Arguments:
        afloat - a float
        ignore_ones - if true, floats of 1 are ignored.
        tol - tolerance to round to nearest int. i.e. 2.0000000001 -> 2
    
    Returns:
        A string representation of the float for formulas.
    """
    if ignore_ones and afloat == 1:
        return ""
    elif abs(afloat - int(afloat)) < tol:
        return str(int(afloat))
    else:
        return str(afloat)

if __name__ == "__main__":
    import doctest
    doctest.testmod()