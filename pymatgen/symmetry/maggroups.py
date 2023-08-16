"""Magnetic space groups."""

from __future__ import annotations

import os
import sqlite3
import textwrap
from array import array
from fractions import Fraction
from typing import TYPE_CHECKING, Sequence

import numpy as np
from monty.design_patterns import cached_class

from pymatgen.core.operations import MagSymmOp
from pymatgen.electronic_structure.core import Magmom
from pymatgen.symmetry.groups import SymmetryGroup, in_array_list
from pymatgen.symmetry.settings import JonesFaithfulTransformation
from pymatgen.util.string import transformation_to_string

if TYPE_CHECKING:
    from pymatgen.core.lattice import Lattice

__author__ = "Matthew Horton, Shyue Ping Ong"

MAGSYMM_DATA = os.path.join(os.path.dirname(__file__), "symm_data_magnetic.sqlite")


@cached_class
class MagneticSpaceGroup(SymmetryGroup):
    """Representation of a magnetic space group."""

    def __init__(self, label, setting_transformation="a,b,c;0,0,0"):
        """
        Initializes a MagneticSpaceGroup from its Belov, Neronova and
        Smirnova (BNS) number supplied as a list or its label supplied
        as a string. To create a magnetic structure in pymatgen, the
        Structure.from_magnetic_spacegroup() method can be used, which
        relies on this class.

        The main difference between magnetic space groups and normal
        crystallographic space groups is the inclusion of a time reversal
        operator that acts on an atom's magnetic moment. This is
        indicated by a prime symbol (') next to the respective symmetry
        operation in its label, e.g. the standard crystallographic
        space group Pnma has magnetic subgroups Pn'ma, Pnm'a, Pnma',
        Pn'm'a, Pnm'a', Pn'ma', Pn'm'a'.

        The magnetic space groups are classified as one of 4 types
        where G = magnetic space group, and F = parent crystallographic
        space group:

        1. G=F no time reversal, i.e. the same as corresponding
            crystallographic group
        2. G=F+F1', "grey" groups, where avg. magnetic moment is zero,
            e.g. a paramagnet in zero ext. mag. field
        3. G=D+(F-D)1', where D is an equi-translation subgroup of F of
            index 2, lattice translations do not include time reversal
        4. G=D+(F-D)1', where D is an equi-class subgroup of F of index 2

        There are two common settings for magnetic space groups, BNS
        and OG. In case 4, the BNS setting != OG setting, and so a
        transformation to go between the two settings is required:
        specifically, the BNS setting is derived from D, and the OG
        setting is derived from F.

        This means that the OG setting refers to the unit cell if magnetic
        order is neglected, and requires multiple unit cells to reproduce
        the full crystal periodicity when magnetic moments are present.
        This does not make the OG setting, in general, useful for
        electronic structure calculations and the BNS setting is preferred.
        However, this class does contain information on the OG setting and
        can be initialized from OG labels or numbers if required.

        Conventions: ITC monoclinic unique axis b, monoclinic cell choice 1,
        hexagonal axis for trigonal groups, origin choice 2 for groups with
        more than one origin choice (ISO-MAG).

        Raw data comes from ISO-MAG, ISOTROPY Software Suite, iso.byu.edu
        http://stokes.byu.edu/iso/magnetic_data.txt
        with kind permission from Professor Branton Campbell, BYU

        Data originally compiled from:
        (1) Daniel B. Litvin, Magnetic Group Tables (International Union
            of Crystallography, 2013) www.iucr.org/publ/978-0-9553602-2-0.
        (2) C. J. Bradley and A. P. Cracknell, The Mathematical Theory of
            Symmetry in Solids (Clarendon Press, Oxford, 1972).

        See http://stokes.byu.edu/iso/magneticspacegroupshelp.php for more
        information on magnetic symmetry.

        :param id: BNS number supplied as list of 2 ints or BNS label as
            str or index as int (1-1651) to iterate over all space groups
        """
        self._data = {}

        # Datafile is stored as sqlite3 database since (a) it can be easily
        # queried for various different indexes (BNS/OG number/labels) and (b)
        # allows binary data to be stored in a compact form similar to that in
        # the source data file, significantly reducing file size.
        # Note that a human-readable JSON format was tested first but was 20x
        # larger and required *much* longer initial loading times.

        # retrieve raw data
        db = sqlite3.connect(MAGSYMM_DATA)
        c = db.cursor()
        if isinstance(label, str):
            label = "".join(label.split())  # remove any white space
            c.execute("SELECT * FROM space_groups WHERE BNS_label=?;", (label,))
        elif isinstance(label, list):
            c.execute("SELECT * FROM space_groups WHERE BNS1=? AND BNS2=?;", (label[0], label[1]))
        elif isinstance(label, int):
            # OG3 index is a 'master' index, going from 1 to 1651
            c.execute("SELECT * FROM space_groups WHERE OG3=?;", (label,))
        raw_data = list(c.fetchone())

        # Jones Faithful transformation
        self.jf = JonesFaithfulTransformation.from_transformation_string("a,b,c;0,0,0")
        if isinstance(setting_transformation, str):
            if setting_transformation != "a,b,c;0,0,0":
                self.jf = JonesFaithfulTransformation.from_transformation_string(setting_transformation)
        elif isinstance(setting_transformation, JonesFaithfulTransformation) and setting_transformation != self.jf:
            self.jf = setting_transformation

        self._data["magtype"] = raw_data[0]  # int from 1 to 4
        self._data["bns_number"] = [raw_data[1], raw_data[2]]
        self._data["bns_label"] = raw_data[3]
        self._data["og_number"] = [raw_data[4], raw_data[5], raw_data[6]]
        self._data["og_label"] = raw_data[7]  # can differ from BNS_label

        def _get_point_operator(idx):
            """Retrieve information on point operator (rotation matrix and Seitz label)."""
            is_hex = self._data["bns_number"][0] >= 143 and self._data["bns_number"][0] <= 194
            c.execute(
                "SELECT symbol, matrix FROM point_operators WHERE idx=? AND hex=?;",
                (idx - 1, is_hex),
            )
            op = c.fetchone()
            return {
                "symbol": op[0],
                "matrix": np.array(op[1].split(","), dtype="f").reshape(3, 3),
            }

        def _parse_operators(b):
            """Parses compact binary representation into list of MagSymmOps."""
            if len(b) == 0:  # e.g. if magtype != 4, OG setting == BNS setting, and b == [] for OG symmops
                return None
            raw_symops = [b[i : i + 6] for i in range(0, len(b), 6)]

            symops = []

            for r in raw_symops:
                point_operator = _get_point_operator(r[0])
                translation_vec = [r[1] / r[4], r[2] / r[4], r[3] / r[4]]
                time_reversal = r[5]
                op = MagSymmOp.from_rotation_and_translation_and_time_reversal(
                    rotation_matrix=point_operator["matrix"],
                    translation_vec=translation_vec,
                    time_reversal=time_reversal,
                )
                # store string representation, e.g. (2x|1/2,1/2,1/2)'
                seitz = (
                    f"({point_operator['symbol']}|"
                    f"{Fraction(translation_vec[0])},{Fraction(translation_vec[1])},{Fraction(translation_vec[2])})"
                )
                if time_reversal == -1:
                    seitz += "'"
                symops.append({"op": op, "str": seitz})

            return symops

        def _parse_wyckoff(b):
            """Parses compact binary representation into list of Wyckoff sites."""
            if len(b) == 0:
                return None

            wyckoff_sites = []

            def get_label(idx):
                if idx <= 25:
                    return chr(97 + idx)  # returns a-z when idx 0-25
                return "alpha"  # when a-z labels exhausted, use alpha, only relevant for a few space groups

            o = 0  # offset
            n = 1  # nth Wyckoff site
            num_wyckoff = b[0]
            while len(wyckoff_sites) < num_wyckoff:
                m = b[1 + o]  # multiplicity
                label = str(b[2 + o] * m) + get_label(num_wyckoff - n)
                sites = []
                for j in range(m):
                    s = b[3 + o + (j * 22) : 3 + o + (j * 22) + 22]  # data corresponding to specific Wyckoff position
                    translation_vec = [s[0] / s[3], s[1] / s[3], s[2] / s[3]]
                    matrix = [
                        [s[4], s[7], s[10]],
                        [s[5], s[8], s[11]],
                        [s[6], s[9], s[12]],
                    ]
                    matrix_magmom = [
                        [s[13], s[16], s[19]],
                        [s[14], s[17], s[20]],
                        [s[15], s[18], s[21]],
                    ]
                    # store string representation, e.g. (x,y,z;mx,my,mz)
                    wyckoff_str = (
                        f"({transformation_to_string(matrix, translation_vec)};"
                        f"{transformation_to_string(matrix_magmom, c='m')})"
                    )
                    sites.append(
                        {
                            "translation_vec": translation_vec,
                            "matrix": matrix,
                            "matrix_magnetic": matrix_magmom,
                            "str": wyckoff_str,
                        }
                    )

                # only keeping string representation of Wyckoff sites for now
                # could do something else with these in future
                wyckoff_sites.append({"label": label, "str": " ".join(s["str"] for s in sites)})
                n += 1
                o += m * 22 + 2

            return wyckoff_sites

        def _parse_lattice(b):
            """Parses compact binary representation into list of lattice vectors/centerings."""
            if len(b) == 0:
                return None
            raw_lattice = [b[i : i + 4] for i in range(0, len(b), 4)]

            lattice = []

            for r in raw_lattice:
                lattice.append(
                    {
                        "vector": [r[0] / r[3], r[1] / r[3], r[2] / r[3]],
                        "str": f"({Fraction(r[0] / r[3]).limit_denominator()},"
                        f"{Fraction(r[1] / r[3]).limit_denominator()},"
                        f"{Fraction(r[2] / r[3]).limit_denominator()})+",
                    }
                )

            return lattice

        def _parse_transformation(b):
            """Parses compact binary representation into transformation between OG and BNS settings."""
            if len(b) == 0:
                return None
            # capital letters used here by convention,
            # IUCr defines P and p specifically
            P = [[b[0], b[3], b[6]], [b[1], b[4], b[7]], [b[2], b[5], b[8]]]
            p = [b[9] / b[12], b[10] / b[12], b[11] / b[12]]
            P = np.array(P).transpose()
            P_string = transformation_to_string(P, components=("a", "b", "c"))
            p_string = (
                f"{Fraction(p[0]).limit_denominator()},"
                f"{Fraction(p[1]).limit_denominator()},"
                f"{Fraction(p[2]).limit_denominator()}"
            )
            return P_string + ";" + p_string

        for i in range(8, 15):
            try:
                raw_data[i] = array("b", raw_data[i])  # construct array from sql binary blobs
            except Exception:
                # array() behavior changed, need to explicitly convert buffer to str in earlier Python
                raw_data[i] = array("b", str(raw_data[i]))

        self._data["og_bns_transform"] = _parse_transformation(raw_data[8])
        self._data["bns_operators"] = _parse_operators(raw_data[9])
        self._data["bns_lattice"] = _parse_lattice(raw_data[10])
        self._data["bns_wyckoff"] = _parse_wyckoff(raw_data[11])
        self._data["og_operators"] = _parse_operators(raw_data[12])
        self._data["og_lattice"] = _parse_lattice(raw_data[13])
        self._data["og_wyckoff"] = _parse_wyckoff(raw_data[14])

        db.close()

    @classmethod
    def from_og(cls, label: Sequence[int] | str) -> MagneticSpaceGroup:
        """
        Initialize from Opechowski and Guccione (OG) label or number.

        :param id: OG number supplied as list of 3 ints or
            or OG label as str
        :return:
        """
        db = sqlite3.connect(MAGSYMM_DATA)
        c = db.cursor()
        if isinstance(label, str):
            c.execute("SELECT BNS_label FROM space_groups WHERE OG_label=?", (label,))
        elif isinstance(label, list):
            c.execute(
                "SELECT BNS_label FROM space_groups WHERE OG1=? and OG2=? and OG3=?",
                (label[0], label[1], label[2]),
            )
        bns_label = c.fetchone()[0]
        db.close()

        return cls(bns_label)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return self._data == other._data

    @property
    def crystal_system(self):
        """Crystal system, e.g., cubic, hexagonal, etc."""
        i = self._data["bns_number"][0]
        if i <= 2:
            return "triclinic"
        if i <= 15:
            return "monoclinic"
        if i <= 74:
            return "orthorhombic"
        if i <= 142:
            return "tetragonal"
        if i <= 167:
            return "trigonal"
        if i <= 194:
            return "hexagonal"
        return "cubic"

    @property
    def sg_symbol(self):
        """Space group symbol."""
        return self._data["bns_label"]

    @property
    def symmetry_ops(self):
        """
        Retrieve magnetic symmetry operations of the space group.
        :return: List of :class:`pymatgen.core.operations.MagSymmOp`.
        """
        ops = [op_data["op"] for op_data in self._data["bns_operators"]]

        # add lattice centerings
        centered_ops = []
        lattice_vectors = [latt["vector"] for latt in self._data["bns_lattice"]]

        for vec in lattice_vectors:
            if not (np.array_equal(vec, [1, 0, 0]) or np.array_equal(vec, [0, 1, 0]) or np.array_equal(vec, [0, 0, 1])):
                for op in ops:
                    new_vec = op.translation_vector + vec
                    new_op = MagSymmOp.from_rotation_and_translation_and_time_reversal(
                        op.rotation_matrix,
                        translation_vec=new_vec,
                        time_reversal=op.time_reversal,
                    )
                    centered_ops.append(new_op)

        ops = ops + centered_ops

        # apply jones faithful transformation
        return [self.jf.transform_symmop(op) for op in ops]

    def get_orbit(self, p, magmom, tol: float = 1e-5):
        """
        Returns the orbit for a point and its associated magnetic moment.

        Args:
            p: Point as a 3x1 array.
            magmom: A magnetic moment, compatible with :class:`pymatgen.electronic_structure.core.Magmom`
            tol: Tolerance for determining if sites are the same. 1e-5 should
                be sufficient for most purposes. Set to 0 for exact matching
                (and also needed for symbolic orbits).

        Returns:
            tuple[list, list]: orbit for point and magnetic moments for orbit.
        """
        orbit: list[np.ndarray] = []
        orbit_magmoms = []
        magmom = Magmom(magmom)
        for sym_op in self.symmetry_ops:
            pp = sym_op.operate(p)
            pp = np.mod(np.round(pp, decimals=10), 1)
            mm = sym_op.operate_magmom(magmom)
            if not in_array_list(orbit, pp, tol=tol):
                orbit.append(pp)
                orbit_magmoms.append(mm)
        return orbit, orbit_magmoms

    def is_compatible(self, lattice: Lattice, tol: float = 1e-5, angle_tol: float = 5) -> bool:
        """
        Checks whether a particular lattice is compatible with the
        *conventional* unit cell.

        Args:
            lattice (Lattice): A Lattice.
            tol (float): The tolerance to check for equality of lengths.
            angle_tol (float): The tolerance to check for equality of angles
                in degrees.

        Returns:
            bool: True if the lattice is compatible with the conventional cell.
        """
        # function from pymatgen.symmetry.groups.SpaceGroup
        abc = lattice.lengths
        angles = lattice.angles
        crys_system = self.crystal_system

        def check(param, ref, tolerance):
            return all(abs(i - j) < tolerance for i, j in zip(param, ref) if j is not None)

        if crys_system == "cubic":
            a = abc[0]
            return check(abc, [a, a, a], tol) and check(angles, [90, 90, 90], angle_tol)
        if crys_system == "hexagonal" or (crys_system == "trigonal" and self.sg_symbol.endswith("H")):
            a = abc[0]
            return check(abc, [a, a, None], tol) and check(angles, [90, 90, 120], angle_tol)
        if crys_system == "trigonal":
            a = abc[0]
            return check(abc, [a, a, a], tol)
        if crys_system == "tetragonal":
            a = abc[0]
            return check(abc, [a, a, None], tol) and check(angles, [90, 90, 90], angle_tol)
        if crys_system == "orthorhombic":
            return check(angles, [90, 90, 90], angle_tol)
        if crys_system == "monoclinic":
            return check(angles, [90, None, 90], angle_tol)
        return True

    def data_str(self, include_og=True):
        """
        Get description of all data, including information for OG setting.
        :return: str.
        """
        # __str__() omits information on OG setting to reduce confusion
        # as to which set of symops are active, this property gives
        # all stored data including OG setting

        desc = {}  # dictionary to hold description strings
        description = ""

        # parse data into strings

        # indicate if non-standard setting specified
        if self.jf != JonesFaithfulTransformation.from_transformation_string("a,b,c;0,0,0"):
            description += "Non-standard setting: .....\n"
            description += repr(self.jf)
            description += "\n\nStandard setting information: \n"

        desc["magtype"] = self._data["magtype"]
        desc["bns_number"] = ".".join(map(str, self._data["bns_number"]))
        desc["bns_label"] = self._data["bns_label"]
        desc["og_id"] = (
            "\t\tOG: " + ".".join(map(str, self._data["og_number"])) + " " + self._data["og_label"]
            if include_og
            else ""
        )
        desc["bns_operators"] = " ".join(op_data["str"] for op_data in self._data["bns_operators"])

        desc["bns_lattice"] = (
            " ".join(lattice_data["str"] for lattice_data in self._data["bns_lattice"][3:])
            if len(self._data["bns_lattice"]) > 3
            else ""
        )  # don't show (1,0,0)+ (0,1,0)+ (0,0,1)+

        desc["bns_wyckoff"] = "\n".join(
            [
                textwrap.fill(
                    wyckoff_data["str"],
                    initial_indent=wyckoff_data["label"] + "  ",
                    subsequent_indent=" " * len(wyckoff_data["label"] + "  "),
                    break_long_words=False,
                    break_on_hyphens=False,
                )
                for wyckoff_data in self._data["bns_wyckoff"]
            ]
        )

        desc["og_bns_transformation"] = (
            f"OG-BNS Transform: ({self._data['og_bns_transform']})\n" if desc["magtype"] == 4 and include_og else ""
        )

        bns_operators_prefix = f"Operators{' (BNS)' if desc['magtype'] == 4 and include_og else ''}: "
        bns_wyckoff_prefix = f"Wyckoff Positions{' (BNS)' if desc['magtype'] == 4 and include_og else ''}: "

        # apply textwrap on long lines
        desc["bns_operators"] = textwrap.fill(
            desc["bns_operators"],
            initial_indent=bns_operators_prefix,
            subsequent_indent=" " * len(bns_operators_prefix),
            break_long_words=False,
            break_on_hyphens=False,
        )

        description += (
            f"BNS: {desc['bns_number']} {desc['bns_label']}{desc['og_id']}\n"
            f"{desc['og_bns_transformation']}"
            f"{desc['bns_operators']}\n"
            f"{bns_wyckoff_prefix}{desc['bns_lattice']}\n"
            f"{desc['bns_wyckoff']}"
        )

        if desc["magtype"] == 4 and include_og:
            desc["og_operators"] = " ".join(op_data["str"] for op_data in self._data["og_operators"])

            # include all lattice vectors because (1,0,0)+ (0,1,0)+ (0,0,1)+
            # not always present in OG setting
            desc["og_lattice"] = " ".join(lattice_data["str"] for lattice_data in self._data["og_lattice"])

            desc["og_wyckoff"] = "\n".join(
                [
                    textwrap.fill(
                        wyckoff_data["str"],
                        initial_indent=wyckoff_data["label"] + "  ",
                        subsequent_indent=" " * len(wyckoff_data["label"] + "  "),
                        break_long_words=False,
                        break_on_hyphens=False,
                    )
                    for wyckoff_data in self._data["og_wyckoff"]
                ]
            )

            og_operators_prefix = "Operators (OG): "

            # apply textwrap on long lines
            desc["og_operators"] = textwrap.fill(
                desc["og_operators"],
                initial_indent=og_operators_prefix,
                subsequent_indent=" " * len(og_operators_prefix),
                break_long_words=False,
                break_on_hyphens=False,
            )

            description += (
                f"\n{desc['og_operators']}\nWyckoff Positions (OG): {desc['og_lattice']}\n{desc['og_wyckoff']}"
            )
        elif desc["magtype"] == 4:
            description += "\nAlternative OG setting exists for this space group."

        return description

    def __str__(self):
        """
        String representation of the space group, specifying the setting
        of the space group, its magnetic symmetry operators and Wyckoff
        positions.
        :return: str.
        """
        return self.data_str(include_og=False)


def _write_all_magnetic_space_groups_to_file(filename):
    """
    Write all magnetic space groups to a human-readable text file.
    Should contain same information as text files provided by ISO-MAG.
    """
    out = (
        "Data parsed from raw data from:\n"
        "ISO-MAG, ISOTROPY Software Suite, iso.byu.edu\n"
        "http://stokes.byu.edu/iso/magnetic_data.txt\n"
        "Used with kind permission from Professor Branton Campbell, BYU\n\n"
    )
    all_msgs = []
    for i in range(1, 1652):
        all_msgs.append(MagneticSpaceGroup(i))
    for msg in all_msgs:
        out += f"\n{msg.data_str()}\n\n--------\n"
    with open(filename, "w") as f:
        f.write(out)
