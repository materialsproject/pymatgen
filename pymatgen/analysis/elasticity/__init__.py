# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import numpy as np

voigt_map = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
reverse_voigt_map = np.array([[0, 5, 4],
                              [5, 1, 3],
                              [4, 3, 2]])
