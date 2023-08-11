"""Development script of the ChemEnv utility to get the equivalent indices of the model coordination environments."""

from __future__ import annotations

import numpy as np

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

if __name__ == "__main__":
    cg_symbol = "O:6"
    equiv_list = []

    # O:6
    if cg_symbol == "O:6":
        opposite_points = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}
        perp_plane = {
            0: [2, 3, 4, 5],
            1: [2, 3, 4, 5],
            2: [0, 1, 4, 5],
            3: [0, 1, 4, 5],
            4: [0, 1, 2, 3],
            5: [0, 1, 2, 3],
        }
        # 0. any point
        for i0 in range(6):
            # 1. point opposite to point 0.
            i1 = opposite_points[i0]
            # 2. one of the 4 points in the perpendicular plane
            for i2 in perp_plane[i0]:
                # 3. point opposite to point 2.
                i3 = opposite_points[i2]
                remaining = list(range(6))
                remaining.remove(i0)
                remaining.remove(i1)
                remaining.remove(i2)
                remaining.remove(i3)
                # 4. one of the 2 remaining points
                for i4 in remaining:
                    # 5. point opposite to point 4.
                    i5 = opposite_points[i4]
                    equiv_list.append([i0, i1, i2, i3, i4, i5])

    # PB:7
    if cg_symbol == "PB:7":
        for i0 in range(5):
            for turn in [1, -1]:
                i1 = np.mod(i0 + turn, 5)
                i2 = np.mod(i1 + turn, 5)
                i3 = np.mod(i2 + turn, 5)
                i4 = np.mod(i3 + turn, 5)
                for i5 in [5, 6]:
                    i6 = 5 if i5 == 6 else 6
                    equiv_list.append([i0, i1, i2, i3, i4, i5, i6])

    # HB:8
    if cg_symbol == "HB:8":
        for i0 in range(6):
            for turn in [1, -1]:
                i1 = np.mod(i0 + turn, 6)
                i2 = np.mod(i1 + turn, 6)
                i3 = np.mod(i2 + turn, 6)
                i4 = np.mod(i3 + turn, 6)
                i5 = np.mod(i4 + turn, 6)
                for i6 in [6, 7]:
                    i7 = 6 if i6 == 7 else 7
                    equiv_list.append([i0, i1, i2, i3, i4, i5, i6, i7])

    # SBT:8
    if cg_symbol == "SBT:8":
        # 0. any point on the square face without cap
        for i0 in [0, 1, 3, 4]:
            # 1. point in this square face but also in the triangular plane of point 0
            # 2. last point in the triangular plane of point 0
            if i0 < 3:
                i1 = 0 if i0 == 1 else 1
                i2 = 2
            else:
                i1 = 3 if i0 == 4 else 4
                i2 = 5
            # 3.4.5. corresponding points in the opposite triangular plane to the one of points 0.1.2.
            i3 = np.mod(i0 + 3, 6)
            i4 = np.mod(i1 + 3, 6)
            i5 = np.mod(i2 + 3, 6)
            # 6. cap point opposite to the first point
            i6 = 7 if i0 in [1, 4] else 6
            # 7. last cap point
            i7 = 6 if i0 in [1, 4] else 7
            equiv_list.append([i0, i1, i2, i3, i4, i5, i6, i7])

    # SA:8
    if cg_symbol == "SA:8":
        sf1 = [0, 2, 1, 3]
        sf2 = [4, 5, 7, 6]
        # 0. any point
        for i0 in range(8):
            # 1. point opposite to point 0. in the square face
            if i0 in [0, 2]:
                i1 = i0 + 1
            elif i0 in [1, 3]:
                i1 = i0 - 1
            elif i0 == 4:
                i1 = 7
            elif i0 == 5:
                i1 = 6
            elif i0 == 6:
                i1 = 5
            elif i0 == 7:
                i1 = 4
            # 2. one of the two last points in the square face
            sfleft = list(sf1) if i0 in sf1 else list(sf2)
            sfleft.remove(i0)
            sfleft.remove(i1)
            for i2 in sfleft:
                sfleft2 = list(sfleft)
                sfleft2.remove(i2)
                # 3. last point in the square face
                i3 = sfleft2[0]
                # 4. point opposite to point 3. and closest to point 0.
                i4 = 0

            # 3.4.5. corresponding points in the opposite triangular plane to the one of points 0.1.2.
            i3 = np.mod(i0 + 3, 6)
            i4 = np.mod(i1 + 3, 6)
            i5 = np.mod(i2 + 3, 6)
            # 6. cap point opposite to the first point
            i6 = 7 if i0 in [1, 4] else 6
            # 7. last cap point
            i7 = 6 if i0 in [1, 4] else 7
            equiv_list.append([i0, i1, i2, i3, i4, i5, i6, i7])

    print(f"Equivalent indices ({len(equiv_list)}) for {cg_symbol} : ")
    print(equiv_list)
