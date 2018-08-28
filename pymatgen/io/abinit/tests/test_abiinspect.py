# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from __future__ import unicode_literals, division, print_function

import os
import tempfile

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinit.abiinspect import *

_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files', "abinit")

try:
    import matplotlib
    have_matplotlib = "DISPLAY" in os.environ
except ImportError:
    have_matplotlib = False


def ref_file(filename):
    return os.path.join(_test_dir, filename)


def ref_files(*filenames):
    return list(map(ref_file, filenames))


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

        # os.remove(filename)


class AbinitInpectTest(PymatgenTest):

    def test_scfcycle(self):
        """Testing ScfCycle."""
        cycle = GroundStateScfCycle.from_file(ref_file("mgb2_scf.abo"))
        str(cycle)
        cycle.to_string(verbose=2)

        assert cycle.num_iterations == 6
        last = cycle.last_iteration

        assert last["Etot(hartree)"] == -7.1476241568657 and last["vres2"] == 3.879E-08
        assert list(cycle["vres2"]) == [1.769E+02, 7.920E-01, 1.570E-01, 4.259E-03, 4.150E-05, 3.879E-08]

        # TODO: Reactivate
        #if have_matplotlib:
        #    assert cycle.plot(show=False)

        # Testing CyclesPlotter.
        p = CyclesPlotter()
        p.add_label_cycle("mgb2 SCF", cycle)
        p.add_label_cycle("same SCF", cycle)

        # TODO: Reactivate
        #if have_matplotlib:
        #    assert p.combiplot(show=False)
        #    p.slideshow()

    def test_relaxation(self):
        """Testing Relaxation object."""
        relaxation = Relaxation.from_file(ref_file("sic_relax.abo"))
        print(relaxation)
        assert len(relaxation) == 4

        assert relaxation[0]["Etot(hartree)"][-1] == -8.8077409200473
        assert relaxation[-1]["Etot(hartree)"][-1] == -8.8234906607147

        for scf_step in relaxation:
            print(scf_step.num_iterations)

        # TODO: Reactivate
        #if have_matplotlib:
        #    relaxation.plot(show=False)
        #    relaxation.slideshow(show=False)
