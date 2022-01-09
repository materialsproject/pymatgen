from unittest import TestCase

from pymatgen.util.sequence import get_chunks


class SequenceUtilsTest(TestCase):
    def setUp(self):
        self.sequence = list(range(100))

    def test_get_chunks(self):
        lengths = [len(chunk) for chunk in get_chunks(self.sequence, 30)]
        self.assertTrue(all(length == 30 for length in lengths[:-1]))
        self.assertEqual(lengths[-1], 10)
