"""Tools and helper functions for abinit calculations"""
import os.path
import shutil

from pymatgen.util.string_utils import list_strings, StringColorizer


class File(object):
    """
    Very simple class used to store file basenames, absolute paths and directory names.

    Provides wrappers for the most commonly used os.path functions.
    """
    def __init__(self, path):
        self._path = os.path.abspath(path)

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    @property
    def path(self):
        """Absolute path of the file."""
        return self._path

    @property
    def basename(self):
        """File basename."""
        return os.path.basename(self.path)

    @property
    def dirname(self):
        """Absolute path of the directory where the file is located."""
        return os.path.dirname(self.path)

    @property
    def exists(self):
        "True if file exists."
        return os.path.exists(self.path)

    @property
    def isncfile(self):
        "True if self is a NetCDF file"
        return self.basename.endswith(".nc")

    def read(self):
        """Read data from file."""
        with open(self.path, "r") as f:
            return f.read()

    def readlines(self):
        """Read lines from files."""
        with open(self.path, "r") as f:
            return f.readlines()

    def write(self, string):
        """Write string to file."""
        self.make_dir()
        with open(self.path, "w") as f:
            return f.write(string)

    def writelines(self, lines):
        """Write a list of strings to file."""
        self.make_dir()
        with open(self.path, "w") as f:
            return f.writelines()

    def make_dir(self):
        """Make the directory where the file is located."""
        if not os.path.exists(self.dirname):
            os.makedirs(self.dirname)


class Directory(object):
    """
    Very simple class that provides helper functions
    wrapping the most commonly used functions of os.path.
    """
    def __init__(self, path):
        self._path = os.path.abspath(path)

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    @property
    def path(self):
        """Absolute path of the directory."""
        return self._path

    @property
    def basename(self):
        """Directory basename."""
        return os.path.basename(self.path)

    def path_join(self, *p):
        """
        Join two or more pathname components, inserting '/' as needed.
        If any component is an absolute path, all previous path components will be discarded.
        """
        return os.path.join(self.path, *p)

    @property
    def exists(self):
        "True if file exists."
        return os.path.exists(self.path)

    def makedirs(self):
        """
        Super-mkdir; create a leaf directory and all intermediate ones.
        Works like mkdir, except that any intermediate path segment (not
        just the rightmost) will be created if it does not exist.
        """
        if not self.exists:
            os.makedirs(self.path)

    def rmtree(self):
        """Recursively delete the directory tree"""
        shutil.rmtree(self.path, ignore_errors=True)

    def path_in(self, filename):
        """Return the absolute path of filename in the directory."""
        return os.path.join(self.path, filename)

    def list_filepaths(self):
        """Return the list of absolute filepaths in the top-level directory."""
        fnames = [f for f in os.listdir(self.path)]
        return [os.path.join(self.path, f) for f in fnames]

    def has_abiext(self, ext):
        """
        Returns the absolute path of the ABINIT file with extension ext.
        Support both Fortran files and netcdf files. In the later case,
        we check whether a file with extension ext + ".nc" is present 
        in the directory. Returns empty string is file is not present.

        Raises:
            ValueError if multiple files with the given ext are found.
            This implies that this method is not compatible with multiple datasets.
        """
        files = []
        for f in self.list_filepaths():
            if f.endswith(ext) or f.endswith(ext + ".nc"):
                files.append(f)

        if not files:
            return ""

        if len(files) > 1:
            # ABINIT users must learn that multiple datasets are bad!
            err_msg = "Found multiple files with the same extensions\n Please avoid the use of mutiple datasets!"
            raise ValueError(err_msg)

        return files[0]

# This dictionary maps ABINIT file extensions to the 
# variables that must be used to read the file in input.
#
# TODO: It would be nice to pass absolute paths to abinit with getden_path
# so that I can avoid creating symbolic links before running but
# the presence of the C-bindings complicates the implementation
# (gfortran SIGFAULTs if I add strings to dataset_type!
_EXT2VARS = {
    "DEN": {"irdden": 1},
    "WFK": {"irdwfk": 1},
    "SCR": {"irdscr": 1},
    "QPS": {"irdqps": 1},
    "1WFK": {"ird1wf": 1},
    "1DEN": {"ird1den": 1},
    "BSR": {"irdbsreso": 1},
    "BSC": {"irdbscoup": 1},
    "HAYDR_SAVE": {"irdhaydock": 1},
}

def irdvars_for_ext(ext):
    """
    Returns a dictionary with the ABINIT variables 
    that must be used to read the file with extension ext.
    """
    return _EXT2VARS[ext].copy()


def abinit_extensions():
    """List with all the ABINIT extensions that are registered."""
    return list(_EXT2VARS.keys())[:]


def find_file(files, ext, prefix=None, dataset=None, image=None):
    """
    Given a list of file names, return the file with extension "_" + ext, None if not found.

    The prefix, the dataset index and the image index can be specified

    .. warning::

       There are some border cases that will confuse the algorithm
       since the order of dataset and image is not tested.
       Solving this problem requires the knowledge of ndtset and nimages
       This code, however should work in 99.9% of the cases.
    """
    separator = "_"

    for filename in list_strings(files):
        # Remove Netcdf extension (if any)
        f = filename[:-3] if filename.endswith(".nc") else filename
        if separator not in f: 
            continue
        tokens = f.split(separator)

        if tokens[-1] == ext:
            found = True
            if prefix is not None:  
                found = found and filename.startswith(prefix)
            if dataset is not None: 
                found = found and "DS" + str(dataset) in tokens
            if image is not None:   
                found = found and "IMG" + str(image) in tokens
            if found: 
                return filename
    else:
        return None
