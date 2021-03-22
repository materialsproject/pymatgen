# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides utility classes for io operations.
"""


import codecs
import errno
import os
import re
import tempfile

from monty.io import zopen

__author__ = "Shyue Ping Ong, Rickard Armiento, Anubhav Jain, G Matteo, Ioannis Petousis"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"


def ask_yesno(question, default=True):
    """
    Args:
        question ():
        default ():

    Returns:

    """
    try:
        answer = input(question)
        return answer.lower().strip() in ["y", "yes"]
    except EOFError:
        return default


def clean_lines(string_list, remove_empty_lines=True):
    """
    Strips whitespace, carriage returns and empty lines from a list of strings.

    Args:
        string_list: List of strings
        remove_empty_lines: Set to True to skip lines which are empty after
            stripping.

    Returns:
        List of clean strings with no whitespaces.
    """

    for s in string_list:
        clean_s = s
        if "#" in s:
            ind = s.index("#")
            clean_s = s[:ind]
        clean_s = clean_s.strip()
        if (not remove_empty_lines) or clean_s != "":
            yield clean_s


def micro_pyawk(filename, search, results=None, debug=None, postdebug=None):
    """
    Small awk-mimicking search routine.

    'file' is file to search through.
    'search' is the "search program", a list of lists/tuples with 3 elements;
    i.e. [[regex,test,run],[regex,test,run],...]
    'results' is a an object that your search program will have access to for
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

    Returns:
        results
    """
    if results is None:
        results = {}

    # Compile strings into regexs
    for entry in search:
        entry[0] = re.compile(entry[0])

    with zopen(filename, "rt") as f:
        for line in f:
            for entry in search:
                match = re.search(entry[0], line)
                if match and (entry[1] is None or entry[1](results, line)):
                    if debug is not None:
                        debug(results, match)
                    entry[2](results, match)
                    if postdebug is not None:
                        postdebug(results, match)

    return results


umask = os.umask(0)
os.umask(umask)


def _maketemp(name, createmode=None):
    """
    Create a temporary file with the filename similar the given ``name``.
    The permission bits are copied from the original file or ``createmode``.
    Returns: the name of the temporary file.
    """
    d, fn = os.path.split(name)
    fd, tempname = tempfile.mkstemp(prefix=".%s-" % fn, dir=d)
    os.close(fd)

    # Temporary files are created with mode 0600, which is usually not
    # what we want. If the original file already exists, just copy its mode.
    # Otherwise, manually obey umask.
    try:
        st_mode = os.lstat(name).st_mode & 0o777
    except OSError as err:
        if err.errno != errno.ENOENT:
            raise
        st_mode = createmode
        if st_mode is None:
            st_mode = ~umask
        st_mode &= 0o666
    os.chmod(tempname, st_mode)

    return tempname


class AtomicFile:
    """
    This is a straight port of Alexander Saltanov's atomicfile package.

    Writeable file object that atomically writes a file.
    All writes will go to a temporary file.
    Call ``close()`` when you are done writing, and AtomicFile will rename
    the temporary copy to the original name, making the changes visible.
    If the object is destroyed without being closed, all your writes are
    discarded.
    If an ``encoding`` argument is specified, codecs.open will be called to open
    the file in the wanted encoding.
    """

    def __init__(self, name, mode="w+b", createmode=None, encoding=None):
        """
        Args:
            name ():
            mode ():
            createmode ():
            encoding ():
        """
        self.__name = name  # permanent name
        self._tempname = _maketemp(name, createmode=createmode)
        if encoding:
            self._fp = codecs.open(self._tempname, mode, encoding)
        else:
            self._fp = open(self._tempname, mode)

        # delegated methods
        self.write = self._fp.write
        self.fileno = self._fp.fileno

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        if exc_type:
            return
        self.close()

    def close(self):
        """
        Close the file.
        """
        if not self._fp.closed:
            self._fp.close()
            # This to avoid:
            #   FileExistsError: [WinError 183] Cannot create a file when that file already exists:
            # On Windows, if dst already exists, OSError will be raised even if it is a file;
            # there may be no way to implement an atomic rename when dst names an existing file.
            if os.name == "nt" and os.path.exists(self.__name):
                os.remove(self.__name)
            os.rename(self._tempname, self.__name)

    def discard(self):
        """
        Discard the file.
        """
        if not self._fp.closed:
            try:
                os.unlink(self._tempname)
            except OSError:
                pass
            self._fp.close()

    def __del__(self):
        if getattr(self, "_fp", None):  # constructor actually did something
            self.discard()
