#!/usr/bin/env python

"""
This module provides utility classes for io operations.
"""

__author__="Shyue Ping Ong, Rickard Armiento"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import gzip
import bz2
import re
import numpy

def file_open_zip_aware(filename, *args):
    """
    This wrapper wraps around the bz2, gzip and standard python's open function 
    to deal intelligently with bzipped, gzipped or standard text files.
    
    Arguments:
        filename:
            filename 
        args:
            Standard args for python open(..).  Examples, 'r' for read, 'w' for write,
    
    Returns:
        File handler
    """
    if filename.split(".")[-1].upper() == "BZ2":
        return bz2.BZ2File(filename,*args)
    elif filename.split(".")[-1].upper() == "GZ" or filename.split(".")[-1] == "z":
        return gzip.GzipFile(filename,*args)
    else:
        return open(filename, *args)

def clean_lines(string_list, remove_empty_lines = True):
    """
    Strips whitespace from strings '\n' and '\r' and empty lines from a list (usually obtained with the ".readlines()" method)
    
    Arguments:
        string_list:
            list of strings
        remove_empty_lines:
            Set to True to skip lines which are empty after stripping.
    
    Returns:
        List of clean strings with no whitespaces.
    """
    
    for s in string_list:
        clean_s = s
        if '#' in s:
            ind = s.index('#')
            clean_s = s[:ind]
        clean_s = clean_s.strip()
        if (not remove_empty_lines) or clean_s != '':
            yield clean_s
    
def micro_pyawk(filename, search, results=None, debug=None, postdebug=None):
    """
    Small awk-mimicking search routine.
       
    'file' is file to search through.
    'search' is the "search program", a list of lists/tuples with 3 elements; i.e.
    [[regex,test,run],[regex,test,run],...]
    'results' is a an object that your search program will have access to for storing results.

    Here regex is either as a Regex object, or a string that we compile into a Regex.
    test and run are callable objects.

    This function goes through each line in filename, and if regex matches that line *and* 
    test(results,line)==True (or test == None) we execute run(results,match), 
    where match is the match object from running Regex.match.   
    
    The default results is an empty dictionary. Passing a results object let you interact 
    with it in run() and test(). Hence, in many occasions it is thus clever to use results=self. 
    
    Author: Rickard Armiento
    
    Returns: 
        results
    """
    if results == None: results = {}

    # Compile strings into regexs
    for entry in search:
        if isinstance(entry[0], str):
            entry[0] = re.compile(entry[0])
            
    reader = file_open_zip_aware(filename)
    for line in reader:
        for i in range(len(search)):
            match = search[i][0].search(line)
            if match and (search[i][1] == None or search[i][1](results,line)):
                if debug != None:
                    debug(results,match)
                search[i][2](results,match)
                if postdebug != None:
                    postdebug(results,match)

    reader.close()
    return results

def clean_json(input_json):
    """
    This method cleans an input json-like dict object, either a list or a dictionary, nested or otherwise,
    by converting all non-string dictionary keys (such as int and float) to strings.
    
    Arguments:
        input_dict:
            input dictionary.
    
    Returns:
        Sanitized dict that can be json serialized.
    """
    if isinstance(input_json, (list, numpy.ndarray)):
        return [clean_json(i) for i in input_json]
    elif isinstance(input_json, dict):
        return {str(k): clean_json(v) for k, v in input_json.items()}
    elif isinstance(input_json, (int, float)):
        return input_json
    elif isinstance(input_json, basestring):
        return str(input_json)
    else:
        return clean_json(input_json.to_dict)
