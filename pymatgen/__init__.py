__author__ = ", ".join(["Shyue Ping Ong", "Anubhav Jain", "Geoffroy Hautier",
                        "William Davidson Richard", "Stephen Dacek",
                        "Michael Kocher", "Dan Gunter", "Shreyas Cholia",
                        "Vincent L Chevrier", "Rickard Armiento"])
__date__ = "May 12 2013"
__version__ = "2.7.2b"

import json
#Useful aliases for commonly used objects and modules.

from .core import *
from .serializers.json_coders import PMGJSONEncoder, PMGJSONDecoder
from .electronic_structure.core import Spin, Orbital
from .util.io_utils import zopen
from .io.smartio import read_structure, write_structure, read_mol, write_mol
from .matproj.rest import MPRester


def pmg_load(filename, **kwargs):
    """
    Loads a json file and deserialize it with PMGJSONDecoder.

    Args:
        filename:
            Filename of file to open. Can be gzipped or bzipped.
        **kwargs:
            Any of the keyword arguments supported by the json.load method.

    Returns:
        Deserialized pymatgen object. Note that these objects can be lists,
        dicts or otherwise nested pymatgen objects that support the to_dict
        and from_dict MSONAble protocol.
    """
    return json.load(zopen(filename), cls=PMGJSONDecoder, **kwargs)


def pmg_dump(obj, filename, **kwargs):
    """
    Dump an object to a json file using PMGJSONEncoder. Note that these
    objects can be lists, dicts or otherwise nested pymatgen objects that
    support the to_dict and from_dict MSONAble protocol.

    Args:
        obj:
            Object to dump.
        filename:
            Filename of file to open. Can be gzipped or bzipped.
        **kwargs:
            Any of the keyword arguments supported by the json.load method.
    """
    return json.dump(obj, zopen(filename, "w"), cls=PMGJSONEncoder, **kwargs)
