# import pytest

from atomate2.jdftx.io.generic_tags import IntTag, TagContainer
from atomate2.jdftx.io.jdftxinfile_master_format import *

dummy_tagcontainer = TagContainer(
    allow_list_representation=True,
    can_repeat=True,
    subtags={
        "s0": IntTag(write_tagname=False, optional=True),
        "s1": IntTag(write_tagname=False, optional=True),
        "s2": IntTag(write_tagname=False, optional=True),
    },
)
dummy_tagcontainer.validate_value_type("s0", [[1]])


# infile = Path(os.getcwd()) / "tests" / "jdftx" / "io" / "example_files" / "example_sp.in"
# testwrite = Path(os.getcwd()) / "tests" / "jdftx" / "io" / "example_files" / "example_sp_copy.in"
# jif = JDFTXInfile.from_file(infile)
# jif.write_file(testwrite)
# jiflist = jif.get_text_list()
# tag_ex = "fluid-anion"

# get_tag_object(tag_ex)
