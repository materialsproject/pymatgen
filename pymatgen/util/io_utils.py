#!/usr/bin/env python

"""
This module provides utility classes for io operations.
"""

__author__ = "Shyue Ping Ong, Rickard Armiento, Anubhav Jain, G Matteo"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import re
import numpy
import os
import time
import errno
from bz2 import BZ2File
from gzip import GzipFile


def zopen(filename, *args, **kwargs):
    """
    This wrapper wraps around the bz2, gzip and standard python's open function
    to deal intelligently with bzipped, gzipped or standard text files.

    Args:
        filename:
            filename
        args:
            Standard args for python open(..). E.g., 'r' for read, 'w' for
            write.
        kwargs:
            Standard kwargs for python open(..).

    Returns:
        File handler
    """
    file_ext = filename.split(".")[-1].upper()
    if file_ext == "BZ2":
        return BZ2File(filename, *args, **kwargs)
    elif file_ext in ("GZ", "Z"):
        return GzipFile(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)


def zpath(filename):
    """
    Returns an existing (zipped or unzipped) file path given the unzipped
    version. If no path exists, returns the filename unmodified

    Args:
        filename:
            filename without zip extension

    Returns:
        filename with a zip extension (unless an unzipped version
        exists)
    """
    for ext in ["", '.gz', '.GZ', '.bz2', '.BZ2', '.z', '.Z']:
        zfilename = "{}{}".format(filename, ext)
        if os.path.exists(zfilename):
            return zfilename
    return filename


def clean_lines(string_list, remove_empty_lines=True):
    """
    Strips whitespace, carriage returns and empty lines from a list of strings.

    Args:
        string_list:
            List of strings
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
    'search' is the "search program", a list of lists/tuples with 3 elements;
    i.e. [[regex,test,run],[regex,test,run],...]
    'results' is a an object that your search program will have access to for
    storing results.

    Here regex is either as a Regex object, or a string that we compile into a
    Regex. test and run are callable objects.

    This function goes through each line in filename, and if regex matches that
    line *and* test(results,line)==True (or test == None) we execute
    run(results,match),where match is the match object from running
    Regex.match.

    The default results is an empty dictionary. Passing a results object let
    you interact with it in run() and test(). Hence, in many occasions it is
    thus clever to use results=self.

    Author: Rickard Armiento

    Returns:
        results
    """
    if results is None:
        results = {}

    # Compile strings into regexs
    for entry in search:
        if isinstance(entry[0], str):
            entry[0] = re.compile(entry[0])

    with zopen(filename) as f:
        for line in f:
            for i in range(len(search)):
                match = search[i][0].search(line)
                if match and (search[i][1] is not None
                              or search[i][1](results, line)):
                    if debug is not None:
                        debug(results, match)
                    search[i][2](results, match)
                    if postdebug is not None:
                        postdebug(results, match)

    return results


def clean_json(input_json, strict=False):
    """
    This method cleans an input json-like dict object, either a list or a dict,
    nested or otherwise, by converting all non-string dictionary keys (such as
    int and float) to strings.

    Args:
        input_dict:
            input dictionary.
        strict:
            This parameters sets the behavior when clean_json encounters an
            object it does not understand. If strict is True, clean_json will
            try to get the to_dict attribute of the object. If no such
            attribute is found, an attribute error will be thrown. If strict is
            False, clean_json will simply call str(object) to convert the
            object to a string representation.

    Returns:
        Sanitized dict that can be json serialized.
    """
    if isinstance(input_json, (list, numpy.ndarray, tuple)):
        return [clean_json(i, strict=strict) for i in input_json]
    elif isinstance(input_json, dict):
        return {str(k): clean_json(v, strict=strict)
                for k, v in input_json.items()}
    elif isinstance(input_json, (int, float)):
        return input_json
    else:
        if not strict:
            return str(input_json)
        else:
            if isinstance(input_json, basestring):
                return str(input_json)
            elif input_json is None:
                return 'None'
            else:
                return clean_json(input_json.to_dict, strict=strict)


def which(program):
    """
    Returns full path to a executable.
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def reverse_readline(m_file, blk_size=4096, max_mem=4000000):
    """
    Generator method to read a file line-by-line, but backwards. This allows
    one to efficiently get data at the end of a file.

    Based on code by Peter Astrand <astrand@cendio.se>, using modifications by
    Raymond Hettinger and Kevin German.
    http://code.activestate.com/recipes/439045-read-a-text-file-backwards
    -yet-another-implementat/

    Reads file forwards and reverses in memory for files smaller than the
    max_mem parameter, or for gzip files where reverse seeks are not supported.

    Files larger than max_mem are dynamically read backwards.

    Args:
        m_file:
            File stream to read (backwards)
        blk_size:
            The buffer size. Defaults to 4096.
        max_mem:
            The maximum amount of memory to involve in this operation. This is
            used to determine when to reverse a file in-memory versus seeking
            portions of a file. For bz2 files, this sets the maximum block
            size.

    Returns:
        Generator that returns lines from the file. Similar behavior to the
        file.readline() method, except the lines are returned from the back
        of the file.
    """

    file_size = os.path.getsize(m_file.name)

    # If the file size is within our desired RAM use, just reverse it in memory
    # GZip files must use this method because there is no way to negative seek
    if file_size < max_mem or isinstance(m_file, GzipFile):
        for line in reversed(m_file.readlines()):
            yield line.rstrip()
    else:
        if isinstance(m_file, BZ2File):
            # for bz2 files, seeks are expensive. It is therefore in our best
            # interest to maximize the blk_size within limits of desired RAM
            # use.
            blk_size = min(max_mem, file_size)

        buf = ""
        m_file.seek(0, 2)
        lastchar = m_file.read(1)
        trailing_newline = (lastchar == "\n")

        while 1:
            newline_pos = buf.rfind("\n")
            pos = m_file.tell()
            if newline_pos != -1:
                # Found a newline
                line = buf[newline_pos + 1:]
                buf = buf[:newline_pos]
                if pos or newline_pos or trailing_newline:
                    line += "\n"
                yield line
            elif pos:
                # Need to fill buffer
                toread = min(blk_size, pos)
                m_file.seek(pos - toread, 0)
                buf = m_file.read(toread) + buf
                m_file.seek(pos - toread, 0)
                if pos == toread:
                    buf = "\n" + buf
            else:
                # Start-of-file
                return


class FileLockException(Exception):
    """Exception raised by FileLock."""


class FileLock(object):
    """
    A file locking mechanism that has context-manager support so you can use
    it in a with statement. This should be relatively cross-compatible as it
    doesn't rely on msvcrt or fcntl for the locking.
    Taken from http://www.evanfosmark.com/2009/01/cross-platform-file-locking
    -support-in-python/
    """
    Error = FileLockException

    def __init__(self, file_name, timeout=10, delay=.05):
        """
        Prepare the file locker. Specify the file to lock and optionally
        the maximum timeout and the delay between each attempt to lock.

        Args:
            file_name:
                Name of file to lock.
            timeout:
                Maximum timeout for locking. Defaults to 10.
            delay:
                Delay between each attempt to lock. Defaults to 0.05.
        """
        self.file_name = os.path.abspath(file_name)
        self.lockfile = os.path.abspath(file_name) + ".lock"
        self.timeout = float(timeout)
        self.delay = float(delay)
        self.is_locked = False

        if self.delay > self.timeout or self.delay <= 0 or self.timeout <= 0:
            raise ValueError("delay and timeout must be positive with delay "
                             "<= timeout")

    def acquire(self):
        """
        Acquire the lock, if possible. If the lock is in use, it check again
        every `delay` seconds. It does this until it either gets the lock or
        exceeds `timeout` number of seconds, in which case it throws
        an exception.
        """
        start_time = time.time()
        while True:
            try:
                self.fd = os.open(self.lockfile,
                                  os.O_CREAT | os.O_EXCL | os.O_RDWR)
                break
            except (OSError,) as e:
                if e.errno != errno.EEXIST:
                    raise
                if (time.time() - start_time) >= self.timeout:
                    raise FileLockException("%s: Timeout occured." %
                                            self.lockfile)
                time.sleep(self.delay)

        self.is_locked = True

    def release(self):
        """ Get rid of the lock by deleting the lockfile.
            When working in a `with` statement, this gets automatically
            called at the end.
        """
        if self.is_locked:
            os.close(self.fd)
            os.unlink(self.lockfile)
            self.is_locked = False

    def __enter__(self):
        """
        Activated when used in the with statement. Should automatically
        acquire a lock to be used in the with block.
        """
        if not self.is_locked:
            self.acquire()
        return self

    def __exit__(self, type, value, traceback):
        """
        Activated at the end of the with statement. It automatically releases
        the lock if it isn't locked.
        """
        if self.is_locked:
            self.release()

    def __del__(self):
        """
        Make sure that the FileLock instance doesn't leave a lockfile
        lying around.
        """
        self.release()
