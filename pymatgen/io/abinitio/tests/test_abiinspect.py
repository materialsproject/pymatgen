from __future__ import division, print_function

import os
import tempfile

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio.abiinspect import *

#test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')


class YamlTokenizerTest(PymatgenTest):
    """Test YamlTokenizer."""
    def test_base(self):
        string = \
"""---
none: [~, null]
bool: [true, false, on, off]
int: 42
float: 3.14159
list: [LITE, RES_ACID, SUS_DEXT]
dict: {hp: 13, sp: 5}
...

this is not a YAML document!
and the tokenizer will ignore it

--- !Monster
name: Cave spider
hp: [2,6]    # 2d6
ac: 16
attacks: [BITE, HURT]
...

This is not a proper document since it does not start with --- 
the end tag below is ignored 
...
--- !Monster
name: Dragon
hp: [2,6]    # 2d6
ac: 32
attacks: [BITE, HURT]
...
"""
        #for i, line in enumerate(string.splitlines()): print(i, line)
        fd, filename = tempfile.mkstemp(text=True)

        with open(filename, "w") as fh:
            fh.write(string)

        doc_tags = [None, "!Monster", "!Monster"]
        doc_linenos = [1, 13, 23]

        with YamlTokenizer(filename) as r:
            # Iterate the docs
            n = 0
            for i, doc in enumerate(r):
                n += 1
                print("doc", doc)
                self.assertTrue(doc.tag == doc_tags[i])
                self.assertTrue(doc.lineno == doc_linenos[i])

            self.assertTrue(n == len(doc_tags))

            # Read all docs present in the file.
            r.seek(0)
            all_docs = r.all_yaml_docs()
            #print(all_docs)
            self.assertTrue(len(all_docs) == 3)

            # We should be at the begining at the file.
            self.assertTrue(all_docs == r.all_yaml_docs())

            # Find documents by tag.
            r.seek(0)
            monster = r.next_doc_with_tag("!Monster")
            #print("monster",monster)
            self.assertTrue(monster == all_docs[1])

            monster = r.next_doc_with_tag("!Monster")
            self.assertTrue(monster == all_docs[2])

            # this should raise StopIteration
            with self.assertRaises(StopIteration):
                monster = r.next_doc_with_tag("!Monster")

        os.remove(filename)

if __name__ == '__main__':
    import unittest
    unittest.main()
