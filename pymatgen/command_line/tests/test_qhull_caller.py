import unittest

from pymatgen.command_line.qhull_caller import qconvex, qvertex, qvoronoi


class QhullCallerTest(unittest.TestCase):

    def test_qconvex(self):
        data = [[-0.5, -0.5], [-0.5, 0.5], [0.5, -0.5], [0.5, 0.5]]
        facets = qconvex(data)
        self.assertEqual(facets, [[0, 2], [1, 0], [2, 3], [3, 1]],
                         "Qconvex gave wrong answers")

    def test_qvertex(self):
        data = [[-0.5, 0, -0.5], [-0.5, 0.5, 0], [0.5, -0.5, 1], [0.5, 0.5, 2],
                [0.5, 1, 2]]
        facets = qvertex(data)
        self.assertEqual(facets, [[8.0, 3.75, -3.75], [0.5, 0.75, 0.75]],
                         "Qvertex gave wrong answers")

    def test_qvoronoi(self):
        data = [[-0.5, 0, -0.5], [-0.5, 0.5, 0], [0.5, -0.5, 1], [0.5, 0.5, 2],
                [0.5, 1, 2]]
        facets = qvoronoi(data)
        self.assertEqual(facets, [[5, 1, 4, 1, 0, 2], [5, 1, 2, 1, 0, 2],
                                  [5, 2, 4, 0, 1, 2]],
                         "Qvoronoi gave wrong answers")


if __name__ == '__main__':
    unittest.main()
