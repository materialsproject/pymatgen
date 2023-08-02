---
layout: default
title: pymatgen.util.io_utils.md
nav_exclude: true
---

# pymatgen.util.io_utils module

This module provides utility classes for io operations.


### pymatgen.util.io_utils.clean_lines(string_list, remove_empty_lines=True)
Strips whitespace, carriage returns and empty lines from a list of strings.


* **Parameters**


    * **string_list** – List of strings


    * **remove_empty_lines** – Set to True to skip lines which are empty after
    stripping.



* **Returns**

    List of clean strings with no whitespaces.



### pymatgen.util.io_utils.micro_pyawk(filename, search, results=None, debug=None, postdebug=None)
Small awk-mimicking search routine.

‘file’ is file to search through.
‘search’ is the “search program”, a list of lists/tuples with 3 elements;
i.e. [[regex,test,run],[regex,test,run],…]
‘results’ is a an object that your search program will have access to for
storing results.

Here regex is either as a Regex object, or a string that we compile into a
Regex. test and run are callable objects.

This function goes through each line in filename, and if regex matches that
line *and* test(results,line)==True (or test is None) we execute
run(results,match),where match is the match object from running
Regex.match.

The default results is an empty dictionary. Passing a results object let
you interact with it in run() and test(). Hence, in many occasions it is
thus clever to use results=self.

Author: Rickard Armiento, Ioannis Petousis


* **Returns**

    results