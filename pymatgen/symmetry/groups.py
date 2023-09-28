"""Defines SymmetryGroup parent class and PointGroup and SpaceGroup classes.
Shyue Ping Ong thanks Marc De Graef for his generous sharing of his
SpaceGroup data as published in his textbook "Structure of Materials".
"""

from __future__ import annotations

import os
import re
import warnings
from abc import ABCMeta, abstractmethod
from collections.abc import Sequence
from fractions import Fraction
from itertools import product
from typing import TYPE_CHECKING, ClassVar, Literal, overload

import numpy as np
from monty.design_patterns import cached_class
from monty.serialization import loadfn

from pymatgen.util.string import Stringify

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    # don't import at runtime to avoid circular import
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.operations import SymmOp  # noqa: TCH004


SYMM_DATA = loadfn(os.path.join(os.path.dirname(__file__), "symm_data.json"))
CrystalSystem = Literal["cubic", "hexagonal", "monoclinic", "orthorhombic", "tetragonal", "triclinic", "trigonal"]


class SymmetryGroup(Sequence, Stringify, metaclass=ABCMeta):
    """Abstract class representing a symmetry group."""

    @property
    @abstractmethod
    def symmetry_ops(self) -> set[SymmOp]:
        """
        Returns:
            List of symmetry operations associated with the group.
        """

    def __contains__(self, item: object) -> bool:
        if not isinstance(item, SymmOp):
            return NotImplemented

        return any(np.allclose(i.affine_matrix, item.affine_matrix) for i in self.symmetry_ops)

    def __hash__(self) -> int:
        return len(self)

    @overload
    def __getitem__(self, item: int) -> SymmOp:
        ...

    @overload
    def __getitem__(self, item: slice) -> Sequence[SymmOp]:
        ...

    def __getitem__(self, item: int | slice) -> SymmOp | Sequence[SymmOp]:
        return list(self.symmetry_ops)[item]

    def __len__(self) -> int:
        return len(self.symmetry_ops)

    def is_subgroup(self, supergroup: SymmetryGroup) -> bool:
        """True if this group is a subgroup of the supplied group.

        Args:
            supergroup (SymmetryGroup): Supergroup to test.

        Returns:
            bool: True if this group is a subgroup of the supplied group.
        """
        warnings.warn("This is not fully functional. Only trivial subsets are tested right now. ")
        return set(self.symmetry_ops).issubset(supergroup.symmetry_ops)

    def is_supergroup(self, subgroup: SymmetryGroup) -> bool:
        """True if this group is a supergroup of the supplied group.

        Args:
            subgroup (SymmetryGroup): Subgroup to test.

        Returns:
            bool: True if this group is a supergroup of the supplied group.
        """
        warnings.warn("This is not fully functional. Only trivial subsets are tested right now. ")
        return set(subgroup.symmetry_ops).issubset(self.symmetry_ops)

    def to_latex_string(self) -> str:
        """
        Returns:
            A latex formatted group symbol with proper subscripts and overlines.
        """
        sym = re.sub(r"_(\d+)", r"$_{\1}$", self.to_pretty_string())
        return re.sub(r"-(\d)", r"$\\overline{\1}$", sym)


@cached_class
class PointGroup(SymmetryGroup):
    """Class representing a Point Group, with generators and symmetry operations.

    Attributes:
        symbol (str): Full International or Hermann-Mauguin Symbol.
        generators (list): List of generator matrices. Note that 3x3 matrices are used for Point Groups.
        symmetry_ops (list): Full set of symmetry operations as matrices.
    """

    def __init__(self, int_symbol: str) -> None:
        """Initializes a Point Group from its international symbol.

        Args:
            int_symbol (str): International or Hermann-Mauguin Symbol.
        """
        from pymatgen.core.operations import SymmOp

        self.symbol = int_symbol
        self.generators = [
            SYMM_DATA["generator_matrices"][enc] for enc in SYMM_DATA["point_group_encoding"][int_symbol]
        ]
        self._symmetry_ops = {SymmOp.from_rotation_and_translation(m) for m in self._generate_full_symmetry_ops()}
        self.order = len(self._symmetry_ops)

    @property
    def symmetry_ops(self) -> set[SymmOp]:
        """
        Returns:
            List of symmetry operations associated with the group.
        """
        return self._symmetry_ops

    def _generate_full_symmetry_ops(self) -> list[SymmOp]:
        symm_ops = list(self.generators)
        new_ops = self.generators
        while len(new_ops) > 0:
            gen_ops = []
            for g1, g2 in product(new_ops, symm_ops):
                op = np.dot(g1, g2)
                if not in_array_list(symm_ops, op):
                    gen_ops.append(op)
                    symm_ops.append(op)
            new_ops = gen_ops
        return symm_ops

    def get_orbit(self, p: ArrayLike, tol: float = 1e-5) -> list[np.ndarray]:
        """Returns the orbit for a point.

        Args:
            p: Point as a 3x1 array.
            tol: Tolerance for determining if sites are the same. 1e-5 should
                be sufficient for most purposes. Set to 0 for exact matching
                (and also needed for symbolic orbits).

        Returns:
            ([array]) Orbit for point.
        """
        orbit: list[np.ndarray] = []
        for o in self.symmetry_ops:
            pp = o.operate(p)
            if not in_array_list(orbit, pp, tol=tol):
                orbit.append(pp)
        return orbit


@cached_class
class SpaceGroup(SymmetryGroup):
    """Class representing a SpaceGroup.

    Attributes:
        symbol (str): Full International or Hermann-Mauguin Symbol.
        int_number (int): International number.
        generators (list): List of generator matrices. Note that 4x4 matrices are used for Space Groups.
        order (int): Order of Space Group.
    """

    SYMM_OPS = loadfn(os.path.join(os.path.dirname(__file__), "symm_ops.json"))
    SG_SYMBOLS: ClassVar[set] = set(SYMM_DATA["space_group_encoding"])
    for op in SYMM_OPS:
        op["hermann_mauguin"] = re.sub(r" ", "", op["hermann_mauguin"])
        op["universal_h_m"] = re.sub(r" ", "", op["universal_h_m"])
        SG_SYMBOLS.add(op["hermann_mauguin"])
        SG_SYMBOLS.add(op["universal_h_m"])

    gen_matrices = SYMM_DATA["generator_matrices"]
    # POINT_GROUP_ENC = SYMM_DATA["point_group_encoding"]
    sg_encoding = SYMM_DATA["space_group_encoding"]
    abbrev_sg_mapping = SYMM_DATA["abbreviated_spacegroup_symbols"]
    translations = {k: Fraction(v) for k, v in SYMM_DATA["translations"].items()}
    full_sg_mapping = {v["full_symbol"]: k for k, v in SYMM_DATA["space_group_encoding"].items()}

    def __init__(self, int_symbol: str) -> None:
        """Initializes a Space Group from its full or abbreviated international
        symbol. Only standard settings are supported.

        Args:
            int_symbol (str): Full International (e.g., "P2/m2/m2/m") or
                Hermann-Mauguin Symbol ("Pmmm") or abbreviated symbol. The
                notation is a LaTeX-like string, with screw axes being
                represented by an underscore. For example, "P6_3/mmc".
                Alternative settings can be accessed by adding a ":identifier".
                For example, the hexagonal setting  for rhombohedral cells can be
                accessed by adding a ":H", e.g., "R-3m:H". To find out all
                possible settings for a spacegroup, use the get_settings()
                classmethod. Alternative origin choices can be indicated by a
                translation vector, e.g., 'Fm-3m(a-1/4,b-1/4,c-1/4)'.
        """
        from pymatgen.core.operations import SymmOp

        int_symbol = re.sub(r" ", "", int_symbol)
        if int_symbol in SpaceGroup.abbrev_sg_mapping:
            int_symbol = SpaceGroup.abbrev_sg_mapping[int_symbol]
        elif int_symbol in SpaceGroup.full_sg_mapping:
            int_symbol = SpaceGroup.full_sg_mapping[int_symbol]

        self._symmetry_ops: set[SymmOp] | None

        for spg in SpaceGroup.SYMM_OPS:
            if int_symbol in [spg["hermann_mauguin"], spg["universal_h_m"]]:
                ops = [SymmOp.from_xyz_string(s) for s in spg["symops"]]
                self.symbol = re.sub(r":", "", re.sub(r" ", "", spg["universal_h_m"]))
                if int_symbol in SpaceGroup.sg_encoding:
                    self.full_symbol = SpaceGroup.sg_encoding[int_symbol]["full_symbol"]
                    self.point_group = SpaceGroup.sg_encoding[int_symbol]["point_group"]
                else:
                    self.full_symbol = re.sub(r" ", "", spg["universal_h_m"])
                    self.point_group = spg["schoenflies"]
                self.int_number = spg["number"]
                self.order = len(ops)
                self._symmetry_ops = {*ops}
                break
        else:
            if int_symbol not in SpaceGroup.sg_encoding:
                raise ValueError(f"Bad international symbol {int_symbol!r}")

            data = SpaceGroup.sg_encoding[int_symbol]

            self.symbol = int_symbol
            # TODO: Support different origin choices.
            enc = list(data["enc"])
            inversion = int(enc.pop(0))
            ngen = int(enc.pop(0))
            symm_ops = [np.eye(4)]
            if inversion:
                symm_ops.append(np.array([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]))
            for _ in range(ngen):
                m = np.eye(4)
                m[:3, :3] = SpaceGroup.gen_matrices[enc.pop(0)]
                m[0, 3] = SpaceGroup.translations[enc.pop(0)]
                m[1, 3] = SpaceGroup.translations[enc.pop(0)]
                m[2, 3] = SpaceGroup.translations[enc.pop(0)]
                symm_ops.append(m)
            self.generators = symm_ops
            self.full_symbol = data["full_symbol"]
            self.point_group = data["point_group"]
            self.int_number = data["int_number"]
            self.order = data["order"]

            self._symmetry_ops = None

    def _generate_full_symmetry_ops(self) -> np.ndarray:
        symm_ops = np.array(self.generators)
        for op in symm_ops:
            op[0:3, 3] = np.mod(op[0:3, 3], 1)
        new_ops = symm_ops
        while len(new_ops) > 0 and len(symm_ops) < self.order:
            gen_ops = []
            for g in new_ops:
                temp_ops = np.einsum("ijk,kl", symm_ops, g)
                for op in temp_ops:
                    op[0:3, 3] = np.mod(op[0:3, 3], 1)
                    ind = np.where(np.abs(1 - op[0:3, 3]) < 1e-5)
                    op[ind, 3] = 0
                    if not in_array_list(symm_ops, op):
                        gen_ops.append(op)
                        symm_ops = np.append(symm_ops, [op], axis=0)
            new_ops = gen_ops  # type: ignore[assignment]
        assert len(symm_ops) == self.order
        return symm_ops

    @classmethod
    def get_settings(cls, int_symbol: str) -> set[str]:
        """Returns all the settings for a particular international symbol.

        Args:
            int_symbol (str): Full International (e.g., "P2/m2/m2/m") or
                Hermann-Mauguin Symbol ("Pmmm") or abbreviated symbol. The
                notation is a LaTeX-like string, with screw axes being
                represented by an underscore. For example, "P6_3/mmc".

        Returns:
            set[str]: All possible settings for the given international symbol.
        """
        symbols = []
        if int_symbol in SpaceGroup.abbrev_sg_mapping:
            symbols.append(SpaceGroup.abbrev_sg_mapping[int_symbol])
            int_number = SpaceGroup.sg_encoding[int_symbol]["int_number"]
        elif int_symbol in SpaceGroup.full_sg_mapping:
            symbols.append(SpaceGroup.full_sg_mapping[int_symbol])
            int_number = SpaceGroup.sg_encoding[int_symbol]["int_number"]
        else:
            for spg in SpaceGroup.SYMM_OPS:
                if int_symbol in [
                    re.split(r"\(|:", spg["hermann_mauguin"])[0],
                    re.split(r"\(|:", spg["universal_h_m"])[0],
                ]:
                    int_number = spg["number"]
                    break

        for spg in SpaceGroup.SYMM_OPS:
            if int_number == spg["number"]:
                symbols.append(spg["hermann_mauguin"])
                symbols.append(spg["universal_h_m"])
        return set(symbols)

    @property
    def symmetry_ops(self) -> set[SymmOp]:
        """Full set of symmetry operations as matrices. Lazily initialized as
        generation sometimes takes a bit of time.
        """
        from pymatgen.core.operations import SymmOp

        if self._symmetry_ops is None:
            self._symmetry_ops = {SymmOp(m) for m in self._generate_full_symmetry_ops()}
        return self._symmetry_ops

    def get_orbit(self, p: ArrayLike, tol: float = 1e-5) -> list[np.ndarray]:
        """Returns the orbit for a point.

        Args:
            p: Point as a 3x1 array.
            tol: Tolerance for determining if sites are the same. 1e-5 should
                be sufficient for most purposes. Set to 0 for exact matching
                (and also needed for symbolic orbits).

        Returns:
            ([array]) Orbit for point.
        """
        orbit: list[np.ndarray] = []
        for o in self.symmetry_ops:
            pp = o.operate(p)
            pp = np.mod(np.round(pp, decimals=10), 1)
            if not in_array_list(orbit, pp, tol=tol):
                orbit.append(pp)
        return orbit

    def get_orbit_and_generators(self, p: ArrayLike, tol: float = 1e-5) -> tuple[list, list]:
        """Returns the orbit and its generators for a point.

        Args:
            p: Point as a 3x1 array.
            tol: Tolerance for determining if sites are the same. 1e-5 should
                be sufficient for most purposes. Set to 0 for exact matching
                (and also needed for symbolic orbits).

        Returns:
            ([array], [array]) Orbit and generators for point.
        """
        from pymatgen.core.operations import SymmOp

        orbit: list[np.ndarray] = [np.array(p, dtype=float)]
        identity = SymmOp.from_rotation_and_translation(np.eye(3), np.zeros(3))
        generators: list[np.ndarray] = [identity]
        for o in self.symmetry_ops:
            pp = o.operate(p)
            pp = np.mod(np.round(pp, decimals=10), 1)
            if not in_array_list(orbit, pp, tol=tol):
                orbit.append(pp)
                generators.append(o)
        return orbit, generators

    def is_compatible(self, lattice: Lattice, tol: float = 1e-5, angle_tol: float = 5) -> bool:
        """Checks whether a particular lattice is compatible with the
        *conventional* unit cell.

        Args:
            lattice (Lattice): A Lattice.
            tol (float): The tolerance to check for equality of lengths.
            angle_tol (float): The tolerance to check for equality of angles
                in degrees.
        """
        abc = lattice.lengths
        angles = lattice.angles
        crys_system = self.crystal_system

        def check(param, ref, tolerance):
            return all(abs(i - j) < tolerance for i, j in zip(param, ref) if j is not None)

        if crys_system == "cubic":
            a = abc[0]
            return check(abc, [a, a, a], tol) and check(angles, [90, 90, 90], angle_tol)
        if crys_system == "hexagonal" or (
            crys_system == "trigonal"
            and (
                self.symbol.endswith("H")
                or self.int_number
                in [143, 144, 145, 147, 149, 150, 151, 152, 153, 154, 156, 157, 158, 159, 162, 163, 164, 165]
            )
        ):
            a = abc[0]
            return check(abc, [a, a, None], tol) and check(angles, [90, 90, 120], angle_tol)
        if crys_system == "trigonal":
            a = abc[0]
            alpha = angles[0]
            return check(abc, [a, a, a], tol) and check(angles, [alpha, alpha, alpha], angle_tol)
        if crys_system == "tetragonal":
            a = abc[0]
            return check(abc, [a, a, None], tol) and check(angles, [90, 90, 90], angle_tol)
        if crys_system == "orthorhombic":
            return check(angles, [90, 90, 90], angle_tol)
        if crys_system == "monoclinic":
            return check(angles, [90, None, 90], angle_tol)
        return True

    @property
    def crystal_system(self) -> CrystalSystem:
        """
        Returns:
            str: Crystal system of the space group, e.g., cubic, hexagonal, etc.
        """
        num = self.int_number
        if num <= 2:
            return "triclinic"
        if num <= 15:
            return "monoclinic"
        if num <= 74:
            return "orthorhombic"
        if num <= 142:
            return "tetragonal"
        if num <= 167:
            return "trigonal"
        if num <= 194:
            return "hexagonal"
        return "cubic"

    def is_subgroup(self, supergroup: SymmetryGroup) -> bool:
        """Check if space group is a subgroup of the supplied symmetry group.

        Args:
            supergroup (Spacegroup): Supergroup to test.

        Returns:
            bool: True if this space group is a subgroup of the supplied group.
        """
        if not isinstance(supergroup, SpaceGroup):
            return NotImplemented

        if len(supergroup.symmetry_ops) < len(self.symmetry_ops):
            return False

        groups = [{supergroup.int_number}]
        all_groups = [supergroup.int_number]
        max_subgroups = {int(k): v for k, v in SYMM_DATA["maximal_subgroups"].items()}
        while True:
            new_sub_groups = set()
            for i in groups[-1]:
                new_sub_groups.update([j for j in max_subgroups[i] if j not in all_groups])
            if self.int_number in new_sub_groups:
                return True

            if len(new_sub_groups) == 0:
                break

            groups.append(new_sub_groups)
            all_groups.extend(new_sub_groups)
        return False

    def is_supergroup(self, subgroup: SymmetryGroup) -> bool:
        """True if this space group is a supergroup of the supplied group.

        Args:
            subgroup (Spacegroup): Subgroup to test.

        Returns:
            bool: True if this space group is a supergroup of the supplied group.
        """
        return subgroup.is_subgroup(self)

    @classmethod
    def from_int_number(cls, int_number: int, hexagonal: bool = True) -> SpaceGroup:
        """Obtains a SpaceGroup from its international number.

        Args:
            int_number (int): International number.
            hexagonal (bool): For rhombohedral groups, whether to return the
                hexagonal setting (default) or rhombohedral setting.

        Raises:
            ValueError: If the international number is not valid, i.e. not between 1 and 230 inclusive.

        Returns:
            SpaceGroup: object with the given international number.
        """
        if int_number not in range(1, 231):
            raise ValueError(f"International number must be between 1 and 230, got {int_number}")
        symbol = sg_symbol_from_int_number(int_number, hexagonal=hexagonal)
        if not hexagonal and int_number in (146, 148, 155, 160, 161, 166, 167):
            symbol += ":R"
        return SpaceGroup(symbol)

    def __repr__(self) -> str:
        symbol = self.symbol
        return f"{type(self).__name__}({symbol=})"

    def __str__(self) -> str:
        return (
            f"Spacegroup {self.symbol} with international number {self.int_number} and order {len(self.symmetry_ops)}"
        )

    def to_pretty_string(self) -> str:
        """
        Returns:
            (str): A pretty string representation of the space group.
        """
        return self.symbol


def sg_symbol_from_int_number(int_number: int, hexagonal: bool = True) -> str:
    """Obtains a SpaceGroup name from its international number.

    Args:
        int_number (int): International number.
        hexagonal (bool): For rhombohedral groups, whether to return the
            hexagonal setting (default) or rhombohedral setting.

    Returns:
        str: Spacegroup symbol
    """
    syms = []
    for n, v in SYMM_DATA["space_group_encoding"].items():
        if v["int_number"] == int_number:
            syms.append(n)
    if len(syms) == 0:
        raise ValueError("Invalid international number!")
    if len(syms) == 2:
        for sym in syms:
            if "e" in sym:
                return sym
        if hexagonal:
            syms = list(filter(lambda s: s.endswith("H"), syms))
        else:
            syms = list(filter(lambda s: not s.endswith("H"), syms))
    return syms.pop()


def in_array_list(array_list: list[np.ndarray] | np.ndarray, arr: np.ndarray, tol: float = 1e-5) -> bool:
    """Extremely efficient nd-array comparison using numpy's broadcasting. This
    function checks if a particular array a, is present in a list of arrays.
    It works for arrays of any size, e.g., even matrix searches.

    Args:
        array_list ([array]): A list of arrays to compare to.
        arr (array): The test array for comparison.
        tol (float): The tolerance. Defaults to 1e-5. If 0, an exact match is done.

    Returns:
        (bool)
    """
    if len(array_list) == 0:
        return False
    axes = tuple(range(1, arr.ndim + 1))
    if not tol:
        return any(np.all(array_list == arr[None, :], axes))
    return any(np.sum(np.abs(array_list - arr[None, :]), axes) < tol)
