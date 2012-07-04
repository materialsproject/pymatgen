__author__ = "Shyue Ping Ong, Anubhav Jain, Michael Kocher, Geoffroy Hautier, Will Richards, Dan Gunter, Shreyas Cholia, Vincent L Chevrier, Rickard Armiento"
__date__ = "Jun 28, 2012"
__version__ = "2.0.0"

"""
Useful aliases for commonly used objects and modules.
"""
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Structure, Molecule, Composition
from pymatgen.core.lattice import Lattice
from pymatgen.serializers.json_coders import PMGJSONEncoder, PMGJSONDecoder
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.util.io_utils import file_open_zip_aware as openz
