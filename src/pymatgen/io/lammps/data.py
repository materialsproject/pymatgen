"""
This module implements a core class LammpsData for generating/parsing
LAMMPS data file, and other bridging classes to build LammpsData from
molecules. This module also implements a subclass CombinedData for
merging LammpsData object.

Only point particle styles are supported for now (atom_style in angle,
atomic, bond, charge, full and molecular only). See the pages below for
more info.

    https://docs.lammps.org/atom_style.html
    https://docs.lammps.org/read_data.html
"""

from __future__ import annotations

import itertools
import re
import warnings
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from monty.io import zopen
from monty.json import MSONable
from monty.serialization import loadfn
from ruamel.yaml import YAML

from pymatgen.core import Element, Lattice, Molecule, Structure
from pymatgen.core.operations import SymmOp
from pymatgen.util.io_utils import clean_lines

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, Literal

    from typing_extensions import Self

    from pymatgen.core.sites import Site
    from pymatgen.core.structure import SiteCollection

__author__ = "Kiran Mathew, Zhi Deng, Tingzheng Hou"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "2.0"
__maintainer__ = "Tingzheng Hou"
__email__ = "tingzheng_hou@berkeley.edu"
__date__ = "May 29, 2021"

MODULE_DIR = Path(__file__).resolve().parent

SECTION_KEYWORDS = {
    "atom": [
        "Atoms",
        "Velocities",
        "Masses",
        "Ellipsoids",
        "Lines",
        "Triangles",
        "Bodies",
    ],
    "topology": ["Bonds", "Angles", "Dihedrals", "Impropers"],
    "ff": [
        "Pair Coeffs",
        "PairIJ Coeffs",
        "Bond Coeffs",
        "Angle Coeffs",
        "Dihedral Coeffs",
        "Improper Coeffs",
    ],
    "class2": [
        "BondBond Coeffs",
        "BondAngle Coeffs",
        "MiddleBondTorsion Coeffs",
        "EndBondTorsion Coeffs",
        "AngleTorsion Coeffs",
        "AngleAngleTorsion Coeffs",
        "BondBond13 Coeffs",
        "AngleAngle Coeffs",
    ],
}

CLASS2_KEYWORDS = {
    "Angle Coeffs": ["BondBond Coeffs", "BondAngle Coeffs"],
    "Dihedral Coeffs": [
        "MiddleBondTorsion Coeffs",
        "EndBondTorsion Coeffs",
        "AngleTorsion Coeffs",
        "AngleAngleTorsion Coeffs",
        "BondBond13 Coeffs",
    ],
    "Improper Coeffs": ["AngleAngle Coeffs"],
}

SECTION_HEADERS = {
    "Masses": ["mass"],
    "Velocities": ["vx", "vy", "vz"],
    "Bonds": ["type", "atom1", "atom2"],
    "Angles": ["type", "atom1", "atom2", "atom3"],
    "Dihedrals": ["type", "atom1", "atom2", "atom3", "atom4"],
    "Impropers": ["type", "atom1", "atom2", "atom3", "atom4"],
}

ATOMS_HEADERS = {
    "angle": ["molecule-ID", "type", "x", "y", "z"],
    "atomic": ["type", "x", "y", "z"],
    "bond": ["molecule-ID", "type", "x", "y", "z"],
    "charge": ["type", "q", "x", "y", "z"],
    "full": ["molecule-ID", "type", "q", "x", "y", "z"],
    "molecular": ["molecule-ID", "type", "x", "y", "z"],
}


class LammpsBox(MSONable):
    """Object for representing a simulation box in LAMMPS settings."""

    def __init__(self, bounds: Sequence, tilt: Sequence | None = None) -> None:
        """
        Args:
            bounds: A (3, 2) array/list of floats setting the
                boundaries of simulation box.
            tilt: A (3,) array/list of floats setting the tilt of
                simulation box. Default to None, i.e., use an
                orthogonal box.
        """
        bounds_arr = np.array(bounds)
        if bounds_arr.shape != (3, 2):
            raise ValueError(f"Expecting a (3, 2) array for bounds, got {bounds_arr.shape}")
        self.bounds = bounds_arr.tolist()
        matrix = np.diag(bounds_arr[:, 1] - bounds_arr[:, 0])

        self.tilt = None
        if tilt is not None:
            tilt_arr = np.array(tilt)
            if tilt_arr.shape != (3,):
                raise ValueError(f"Expecting a (3,) array for box_tilt, got {tilt_arr.shape}")
            self.tilt = tilt_arr.tolist()
            matrix[1, 0] = tilt_arr[0]
            matrix[2, 0] = tilt_arr[1]
            matrix[2, 1] = tilt_arr[2]
        self._matrix = matrix

    def __str__(self) -> str:
        return self.get_str()

    def __repr__(self) -> str:
        return self.get_str()

    @property
    def volume(self) -> float:
        """Volume of simulation box."""
        matrix = self._matrix
        return np.dot(np.cross(matrix[0], matrix[1]), matrix[2])

    def get_str(self, significant_figures: int = 6) -> str:
        """Get the string representation of simulation box in LAMMPS data file format.

        Args:
            significant_figures (int): No. of significant figures to
                output for box settings. Default to 6.
        """
        ph = f"{{:.{significant_figures}f}}"
        lines = []
        for bound, d in zip(self.bounds, "xyz", strict=True):
            fillers = bound + [d] * 2
            bound_format = " ".join([ph] * 2 + [" {}lo {}hi"])
            lines.append(bound_format.format(*fillers))
        if self.tilt:
            tilt_format = " ".join([ph] * 3 + [" xy xz yz"])
            lines.append(tilt_format.format(*self.tilt))
        return "\n".join(lines)

    def get_box_shift(self, i: Sequence[int]) -> np.ndarray:
        """
        Calculates the coordinate shift due to PBC.

        Args:
            i: A (n, 3) integer array containing the labels for box
            images of n entries.

        Returns:
            Coordinate shift array with the same shape of i
        """
        return np.inner(i, self._matrix.T)

    def to_lattice(self) -> Lattice:
        """Convert the simulation box to a more powerful Lattice backend.
        Note that Lattice is always periodic in 3D space while a
        simulation box is not necessarily periodic in all dimensions.

        Returns:
            Lattice
        """
        return Lattice(self._matrix)


def lattice_2_lmpbox(lattice: Lattice, origin: Sequence = (0, 0, 0)) -> tuple[LammpsBox, SymmOp]:
    """
    Converts a lattice object to LammpsBox, and calculates the symmetry
    operation used.

    Args:
        lattice (Lattice): Input lattice.
        origin: A (3,) array/list of floats setting lower bounds of
            simulation box. Default to (0, 0, 0).

    Returns:
        tuple[LammpsBox, SymmOp]
    """
    a, b, c = lattice.abc
    xlo, ylo, zlo = origin
    xhi = a + xlo
    matrix = lattice.matrix
    xy = np.dot(matrix[1], matrix[0] / a)
    yhi = np.sqrt(b**2 - xy**2) + ylo
    xz = np.dot(matrix[2], matrix[0] / a)
    yz = (np.dot(matrix[1], matrix[2]) - xy * xz) / (yhi - ylo)
    zhi = np.sqrt(c**2 - xz**2 - yz**2) + zlo
    tilt = None if lattice.is_orthogonal else [xy, xz, yz]
    rot_matrix = np.linalg.solve([[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, zhi - zlo]], matrix)
    bounds = [[xlo, xhi], [ylo, yhi], [zlo, zhi]]
    symm_op = SymmOp.from_rotation_and_translation(rot_matrix, origin)
    return LammpsBox(bounds, tilt), symm_op


class LammpsData(MSONable):
    """Object for representing the data in a LAMMPS data file."""

    def __init__(
        self,
        box: LammpsBox,
        masses: pd.DataFrame,
        atoms: pd.DataFrame,
        velocities: pd.DataFrame = None,
        force_field: dict | None = None,
        topology: dict[str, pd.DataFrame] | None = None,
        atom_style: str = "full",
    ) -> None:
        """Low level constructor designed to work with parsed data or other bridging
        objects (ForceField and Topology). Not recommended to use directly.

        Args:
            box (LammpsBox): Simulation box.
            masses (pandas.DataFrame): DataFrame with one column
                ["mass"] for Masses section.
            atoms (pandas.DataFrame): DataFrame with multiple columns
                for Atoms section. Column names vary with atom_style.
            velocities (pandas.DataFrame): DataFrame with three columns
                ["vx", "vy", "vz"] for Velocities section. Optional
                with default to None. If not None, its index should be
                consistent with atoms.
            force_fieldct (dict): Data for force field sections. Optional
                with default to None. Only keywords in force field and
                class 2 force field are valid keys, and each value is a
                DataFrame.
            topology (dict): Data for topology sections. Optional with
                default to None. Only keywords in topology are valid
                keys, and each value is a DataFrame.
            atom_style (str): Output atom_style. Default to "full".
        """
        if velocities is not None and len(atoms) != len(velocities):
            raise ValueError(f"{len(atoms)=} and {len(velocities)=} mismatch")

        if force_field:
            all_ff_kws = SECTION_KEYWORDS["ff"] + SECTION_KEYWORDS["class2"]
            force_field = {key: values for key, values in force_field.items() if key in all_ff_kws}

        if topology:
            topology = {key: values for key, values in topology.items() if key in SECTION_KEYWORDS["topology"]}

        self.box = box
        self.masses = masses
        self.atoms = atoms
        self.velocities = velocities
        self.force_field = force_field
        self.topology = topology
        self.atom_style = atom_style

    def __str__(self) -> str:
        return self.get_str()

    def __repr__(self) -> str:
        return self.get_str()

    @property
    def structure(self) -> Structure:
        """Exports a periodic structure object representing the simulation
        box.

        Returns:
            Structure
        """
        masses = self.masses
        atoms = self.atoms.copy()
        if "nx" in atoms.columns:
            atoms = atoms.drop(["nx", "ny", "nz"], axis=1)
        atoms["molecule-ID"] = 1
        ld_copy = type(self)(self.box, masses, atoms)
        topologies = ld_copy.disassemble()[-1]
        molecule = topologies[0].sites
        coords = molecule.cart_coords - np.array(self.box.bounds)[:, 0]
        species = molecule.species
        lattice = self.box.to_lattice()
        site_properties = {}
        if "q" in atoms:
            site_properties["charge"] = atoms["q"].to_numpy()
        if self.velocities is not None:
            site_properties["velocities"] = self.velocities.to_numpy()
        return Structure(
            lattice,
            species,
            coords,
            coords_are_cartesian=True,
            site_properties=site_properties,
        )

    def get_str(self, distance: int = 6, velocity: int = 8, charge: int = 4, hybrid: bool = True) -> str:
        """Get the string representation of LammpsData, essentially the string
        to be written to a file. Supports hybrid style coeffs read and write.

        Args:
            distance (int): No. of significant figures to output for
                box settings (bounds and tilt) and atomic coordinates.
                Default to 6.
            velocity (int): No. of significant figures to output for
                velocities. Default to 8.
            charge (int): No. of significant figures to output for
                charges. Default to 4.
            hybrid (bool): Whether to write hybrid coeffs types.
                Default to True. If the data object has no hybrid
                coeffs types and has large coeffs section, one may
                use False to speed up the process. Otherwise, the
                default is recommended.

        Returns:
            str: String representation of LammpsData.
        """
        file_template = """Generated by pymatgen.io.lammps.data.LammpsData

{stats}

{box}

{body}
"""
        box = self.box.get_str(distance)

        body_dict = {}
        body_dict["Masses"] = self.masses
        types = {}
        types["atom"] = len(self.masses)
        if self.force_field:
            all_ff_kws = SECTION_KEYWORDS["ff"] + SECTION_KEYWORDS["class2"]
            ff_kws = [k for k in all_ff_kws if k in self.force_field]
            for kw in ff_kws:
                body_dict[kw] = self.force_field[kw]
                if kw in SECTION_KEYWORDS["ff"][2:]:
                    types[kw.lower()[:-7]] = len(self.force_field[kw])

        body_dict["Atoms"] = self.atoms
        counts = {}
        counts["atoms"] = len(self.atoms)
        if self.velocities is not None:
            body_dict["Velocities"] = self.velocities
        if self.topology:
            for kw in SECTION_KEYWORDS["topology"]:
                if kw in self.topology:
                    body_dict[kw] = self.topology[kw]
                    counts[kw.lower()] = len(self.topology[kw])

        all_stats = list(counts.values()) + list(types.values())
        right_indent = len(str(max(all_stats)))
        count_lines = [f"{v:>{right_indent}}  {k}" for k, v in counts.items()]
        type_lines = [f"{v:>{right_indent}}  {k + ' types'}" for k, v in types.items()]
        stats = "\n".join([*count_lines, "", *type_lines])

        def map_coords(q) -> str:
            return f"{q:.{distance}f}"

        def map_velos(q) -> str:
            return f"{q:.{velocity}f}"

        def map_charges(q) -> str:
            return f"{q:.{charge}f}"

        float_format = "{:.9f}".format
        float_format_2 = "{:.1f}".format
        int_format = "{:.0f}".format
        default_formatters = {
            "x": map_coords,
            "y": map_coords,
            "z": map_coords,
            "vx": map_velos,
            "vy": map_velos,
            "vz": map_velos,
            "q": map_charges,
        }
        coeffs_data_type = loadfn(str(MODULE_DIR / "CoeffsDataType.yaml"))
        coeffs: dict[str, dict] = {}
        for style, types in coeffs_data_type.items():
            coeffs[style] = {}
            for typ, formatter in types.items():
                coeffs[style][typ] = {}
                for coeff, datatype in formatter.items():  # type: ignore[attr-defined]
                    if datatype == "int_format":
                        coeffs[style][typ][coeff] = int_format
                    elif datatype == "float_format_2":
                        coeffs[style][typ][coeff] = float_format_2
                    else:
                        coeffs[style][typ][coeff] = float_format

        section_template = "{kw}\n\n{df}\n"
        parts = []
        for key, val in body_dict.items():
            index = key != "PairIJ Coeffs"
            if hybrid and key in ["Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs"]:
                dfs: list[pd.DataFrame] = np.array_split(val, len(val.index))
                df_string = ""
                for idx, df in enumerate(dfs):
                    if isinstance(df.iloc[0]["coeff1"], str):
                        try:
                            formatters = {
                                **default_formatters,
                                **coeffs[key][df.iloc[0]["coeff1"]],
                            }
                        except KeyError:
                            formatters = default_formatters
                        line_string = df.to_string(
                            header=False,
                            formatters=formatters,
                            index_names=False,
                            index=index,
                            na_rep="",
                        )
                    else:
                        line_string = val.to_string(
                            header=False,
                            formatters=default_formatters,
                            index_names=False,
                            index=index,
                            na_rep="",
                        ).splitlines()[idx]
                    df_string += line_string.replace("nan", "").rstrip() + "\n"
            else:
                df_string = val.to_string(
                    header=False, formatters=default_formatters, index_names=False, index=index, na_rep=""
                )
            parts.append(section_template.format(kw=key, df=df_string))
        body = "\n".join(parts)

        return file_template.format(stats=stats, box=box, body=body)

    def write_file(self, filename: str, distance: int = 6, velocity: int = 8, charge: int = 4) -> None:
        """Write LammpsData to file.

        Args:
            filename (str): Filename.
            distance (int): No. of significant figures to output for
                box settings (bounds and tilt) and atomic coordinates.
                Default to 6.
            velocity (int): No. of significant figures to output for
                velocities. Default to 8.
            charge (int): No. of significant figures to output for
                charges. Default to 4.
        """
        with open(filename, mode="w", encoding="utf-8") as file:
            file.write(self.get_str(distance=distance, velocity=velocity, charge=charge))

    def disassemble(
        self, atom_labels: Sequence[str] | None = None, guess_element: bool = True, ff_label: str = "ff_map"
    ) -> tuple[LammpsBox, ForceField, list[Topology]]:
        """
        Breaks down LammpsData to building blocks
        (LammpsBox, ForceField and a series of Topology).
        RESTRICTIONS APPLIED:

        1. No complex force field defined not just on atom
            types, where the same type or equivalent types of topology
            may have more than one set of coefficients.
        2. No intermolecular topologies (with atoms from different
            molecule-ID) since a Topology object includes data for ONE
            molecule or structure only.

        Args:
            atom_labels ([str]): List of strings (must be different
                from one another) for labelling each atom type found in
                Masses section. Default to None, where the labels are
                automatically added based on either element guess or
                dummy specie assignment.
            guess_element (bool): Whether to guess the element based on
                its atomic mass. Default to True, otherwise dummy
                species "Qa", "Qb", ... will be assigned to various
                atom types. The guessed or assigned elements will be
                reflected on atom labels if atom_labels is None, as
                well as on the species of molecule in each Topology.
            ff_label (str): Site property key for labeling atoms of
                different types. Default to "ff_map".

        Returns:
            LammpsBox, ForceField, [Topology]
        """
        atoms_df = self.atoms.copy()
        if "nx" in atoms_df.columns:
            atoms_df[["x", "y", "z"]] += self.box.get_box_shift(atoms_df[["nx", "ny", "nz"]].values)
        atoms_df = pd.concat([atoms_df, self.velocities], axis=1)

        mids = atoms_df.get("molecule-ID")
        if mids is None:
            unique_mids = [1]
            data_by_mols = {1: {"Atoms": atoms_df}}
        else:
            unique_mids = np.unique(mids)
            data_by_mols = {}
            for k in unique_mids:
                data_by_mols[k] = {"Atoms": atoms_df[atoms_df["molecule-ID"] == k]}

        masses = self.masses.copy()
        masses["label"] = atom_labels
        unique_masses = np.unique(masses["mass"])
        if guess_element:
            ref_masses = {el.name: el.atomic_mass.real for el in Element}
            symbols = list(ref_masses)
            diff = np.abs(np.array(list(ref_masses.values())) - unique_masses[:, None])
            symbols = [Element(symbols[idx]).symbol for idx in np.argmin(diff, axis=1)]
        else:
            symbols = [f"Q{a}" for a in map(chr, range(97, 97 + len(unique_masses)))]
        for um, s in zip(unique_masses, symbols, strict=True):
            masses.loc[masses["mass"] == um, "element"] = s
        if atom_labels is None:  # add unique labels based on elements
            for el, vc in masses["element"].value_counts().items():
                masses.loc[masses["element"] == el, "label"] = [f"{el}{c}" for c in range(1, vc + 1)]
        if masses["label"].nunique(dropna=False) != len(masses):
            raise ValueError("Expecting unique atom label for each type")
        mass_info = [(row.label, row.mass) for row in masses.itertuples()]

        non_bond_coeffs: list = []
        topo_coeffs: dict = {}
        if self.force_field:
            if "PairIJ Coeffs" in self.force_field:
                nbc = self.force_field["PairIJ Coeffs"]
                nbc = nbc.sort_values(["id1", "id2"]).drop(["id1", "id2"], axis=1)
                non_bond_coeffs = [list(t) for t in nbc.itertuples(index=False, name=None)]
            elif "Pair Coeffs" in self.force_field:
                nbc = self.force_field["Pair Coeffs"].sort_index()
                non_bond_coeffs = [list(t) for t in nbc.itertuples(index=False, name=None)]

            topo_coeffs = {k: [] for k in SECTION_KEYWORDS["ff"][2:] if k in self.force_field}
            for kw in topo_coeffs:
                class2_coeffs = {
                    k: list(v.itertuples(index=False, name=None))
                    for k, v in self.force_field.items()
                    if k in CLASS2_KEYWORDS.get(kw, [])
                }
                ff_df = self.force_field[kw]
                for t in ff_df.itertuples(index=True, name=None):
                    coeffs_dict = {"coeffs": list(t[1:]), "types": []}
                    if class2_coeffs:
                        coeffs_dict |= {k: list(v[t[0] - 1]) for k, v in class2_coeffs.items()}
                    topo_coeffs[kw].append(coeffs_dict)

        if self.topology:

            def label_topo(t) -> tuple:
                return tuple(masses.loc[atoms_df.loc[t, "type"], "label"])

            for key, values in self.topology.items():
                ff_kw = key[:-1] + " Coeffs"
                for topo in values.itertuples(index=False, name=None):
                    topo_idx = topo[0] - 1
                    indices = list(topo[1:])
                    mids = atoms_df.loc[indices]["molecule-ID"].unique()
                    if len(mids) != 1:
                        raise RuntimeError(
                            "Do not support intermolecular topology formed by atoms with different molecule-IDs"
                        )
                    label = label_topo(indices)
                    topo_coeffs[ff_kw][topo_idx]["types"].append(label)
                    if data_by_mols[mids[0]].get(key):
                        data_by_mols[mids[0]][key].append(indices)
                    else:
                        data_by_mols[mids[0]][key] = [indices]

        if any(topo_coeffs):
            for v in topo_coeffs.values():
                for coeffs_dict in v:
                    coeffs_dict["types"] = list(set(coeffs_dict["types"]))

        ff = ForceField(
            mass_info=mass_info,
            nonbond_coeffs=non_bond_coeffs if any(non_bond_coeffs) else None,
            topo_coeffs=topo_coeffs if any(topo_coeffs) else None,
        )

        topo_list = []
        for mid in unique_mids:
            data = data_by_mols[mid]
            atoms = data["Atoms"]
            shift = min(atoms.index)
            type_ids = atoms["type"]
            species = masses.loc[type_ids, "element"]
            labels = masses.loc[type_ids, "label"]
            coords = atoms[["x", "y", "z"]]
            mol = Molecule(species.values, coords.values, site_properties={ff_label: labels.to_numpy()})
            charges = atoms.get("q")
            velocities = atoms[["vx", "vy", "vz"]] if "vx" in atoms.columns else None
            topologies = {}
            for kw in SECTION_KEYWORDS["topology"]:
                if data.get(kw):
                    topologies[kw] = (np.array(data[kw]) - shift).tolist()
            topo_list.append(
                Topology(
                    sites=mol,
                    ff_label=ff_label,
                    charges=charges,
                    velocities=velocities,
                    topologies=topologies or None,
                )
            )

        return self.box, ff, topo_list

    @classmethod
    def from_file(cls, filename: str, atom_style: str = "full", sort_id: bool = False) -> Self:
        """
        Constructor that parses a file.

        Args:
            filename (str): Filename to read.
            atom_style (str): Associated atom_style. Default to "full".
            sort_id (bool): Whether sort each section by id. Default to
                True.
        """
        with zopen(filename, mode="rt") as file:
            lines = file.readlines()
        kw_pattern = r"|".join(itertools.chain(*SECTION_KEYWORDS.values()))
        section_marks = [idx for idx, line in enumerate(lines) if re.search(kw_pattern, line)]
        parts = np.split(lines, section_marks)

        float_group = r"([0-9eE.+-]+)"
        header_pattern: dict[str, str] = {}
        header_pattern["counts"] = r"^\s*(\d+)\s+([a-zA-Z]+)$"
        header_pattern["types"] = r"^\s*(\d+)\s+([a-zA-Z]+)\s+types$"
        header_pattern["bounds"] = r"^\s*{}$".format(r"\s+".join([float_group] * 2 + [r"([xyz])lo \3hi"]))
        header_pattern["tilt"] = r"^\s*{}$".format(r"\s+".join([float_group] * 3 + ["xy xz yz"]))

        header: dict[str, Any] = {"counts": {}, "types": {}}
        bounds: dict[str, list[float]] = {}
        for line in clean_lines(parts[0][1:]):  # skip the 1st line
            match = None
            key = None
            for key, val in header_pattern.items():  # noqa: B007
                if match := re.match(val, line):
                    break
            if match and key in {"counts", "types"}:
                header[key][match[2]] = int(match[1])
            elif match and key == "bounds":
                g = match.groups()
                bounds[g[2]] = [float(i) for i in g[:2]]
            elif match and key == "tilt":
                header["tilt"] = [float(i) for i in match.groups()]
        header["bounds"] = [bounds.get(i, [-0.5, 0.5]) for i in "xyz"]
        box = LammpsBox(header["bounds"], header.get("tilt"))

        def parse_section(sec_lines) -> tuple[str, pd.DataFrame]:
            title_info = sec_lines[0].split("#", 1)
            kw = title_info[0].strip()
            str_io = StringIO("".join(sec_lines[2:]))  # skip the 2nd line
            if kw.endswith("Coeffs") and not kw.startswith("PairIJ"):
                dfs = [
                    pd.read_csv(StringIO(line), header=None, comment="#", sep=r"\s+")
                    for line in sec_lines[2:]
                    if line.strip()
                ]
                df_section = pd.concat(dfs, ignore_index=True)
                names = ["id"] + [f"coeff{i}" for i in range(1, df_section.shape[1])]
            else:
                df_section = pd.read_csv(str_io, header=None, comment="#", sep=r"\s+")
                if kw == "PairIJ Coeffs":
                    names = ["id1", "id2"] + [f"coeff{i}" for i in range(1, df_section.shape[1] - 1)]
                    df_section.index.name = None
                elif kw in SECTION_HEADERS:
                    names = ["id"] + SECTION_HEADERS[kw]
                elif kw == "Atoms":
                    names = ["id"] + ATOMS_HEADERS[atom_style]
                    if df_section.shape[1] == len(names):
                        pass
                    elif df_section.shape[1] == len(names) + 3:
                        names += ["nx", "ny", "nz"]
                    else:
                        raise ValueError(f"Format in Atoms section inconsistent with {atom_style=}")
                else:
                    raise NotImplementedError(f"Parser for {kw} section not implemented")
            df_section.columns = names
            if sort_id:
                sort_by = "id" if kw != "PairIJ Coeffs" else ["id1", "id2"]
                df_section = df_section.sort_values(sort_by)
            if "id" in df_section.columns:
                df_section = df_section.set_index("id", drop=True)
                df_section.index.name = None
            return kw, df_section

        err_msg = "Bad LAMMPS data format where "
        body = {}
        seen_atoms = False
        for part in parts[1:]:
            name, df_section = parse_section(part)
            if name == "Atoms":
                seen_atoms = True
            if (
                name in ["Velocities"] + SECTION_KEYWORDS["topology"] and not seen_atoms
            ):  # Atoms must appear earlier than these
                raise RuntimeError(f"{err_msg}{name} section appears before Atoms section")
            body[name] = df_section

        err_msg += "Nos. of {} do not match between header and {} section"
        if len(body["Masses"]) != header["types"]["atom"]:
            raise RuntimeError(err_msg.format("atom types", "Masses"))
        atom_sections = ["Atoms", "Velocities"] if "Velocities" in body else ["Atoms"]
        for atom_sec in atom_sections:
            if len(body[atom_sec]) != header["counts"]["atoms"]:
                raise RuntimeError(err_msg.format("atoms", atom_sec))
        for atom_sec in SECTION_KEYWORDS["topology"]:
            if (
                header["counts"].get(atom_sec.lower(), 0) > 0
                and len(body[atom_sec]) != header["counts"][atom_sec.lower()]
            ):
                raise RuntimeError(err_msg.format(atom_sec.lower(), atom_sec))

        items = {k.lower(): body[k] for k in ["Masses", "Atoms"]}
        items["velocities"] = body.get("Velocities")
        ff_kws = [k for k in body if k in SECTION_KEYWORDS["ff"] + SECTION_KEYWORDS["class2"]]
        items["force_field"] = {k: body[k] for k in ff_kws} if ff_kws else None
        topo_kws = [k for k in body if k in SECTION_KEYWORDS["topology"]]
        items["topology"] = {k: body[k] for k in topo_kws} if topo_kws else None
        items["atom_style"] = atom_style
        items["box"] = box
        return cls(**items)

    @classmethod
    def from_ff_and_topologies(
        cls, box: LammpsBox, ff: ForceField, topologies: Sequence[Topology], atom_style: str = "full"
    ) -> Self:
        """
        Constructor building LammpsData from a ForceField object and a
        list of Topology objects. Do not support intermolecular
        topologies since a Topology object includes data for ONE
        molecule or structure only.

        Args:
            box (LammpsBox): Simulation box.
            ff (ForceField): ForceField object with data for Masses and
                force field sections.
            topologies ([Topology]): List of Topology objects with data
                for Atoms, Velocities and topology sections.
            atom_style (str): Output atom_style. Default to "full".
        """
        atom_types = set.union(*(t.species for t in topologies))
        if not atom_types.issubset(ff.maps["Atoms"]):
            raise ValueError("Unknown atom type found in topologies")

        items = {"box": box, "atom_style": atom_style, "masses": ff.masses, "force_field": ff.force_field}

        mol_ids: list[int] = []
        charges: list[float] = []
        coords: list[np.ndarray] = []
        labels: list[str] = []
        v_collector: list | None = [] if topologies[0].velocities else None
        topo_collector: dict[str, list] = {"Bonds": [], "Angles": [], "Dihedrals": [], "Impropers": []}
        topo_labels: dict[str, list] = {"Bonds": [], "Angles": [], "Dihedrals": [], "Impropers": []}
        for idx, topo in enumerate(topologies):
            if topo.topologies:
                shift = len(labels)
                for k, v in topo.topologies.items():
                    topo_collector[k].append(np.array(v) + shift + 1)
                    topo_labels[k].extend([tuple(topo.type_by_sites[j] for j in t) for t in v])
            if isinstance(v_collector, list):
                v_collector.append(topo.velocities)
            mol_ids.extend([idx + 1] * len(topo.sites))
            labels.extend(topo.type_by_sites)
            coords.append(topo.sites.cart_coords)
            charges.extend(topo.charges or [0.0] * len(topo.sites))

        atoms = pd.DataFrame(np.concatenate(coords), columns=["x", "y", "z"])
        atoms["molecule-ID"] = mol_ids
        atoms["q"] = charges
        atoms["type"] = list(map(ff.maps["Atoms"].get, labels))
        atoms.index += 1
        atoms = atoms[ATOMS_HEADERS[atom_style]]

        velocities = None
        if v_collector:
            velocities = pd.DataFrame(np.concatenate(v_collector), columns=SECTION_HEADERS["Velocities"])
            velocities.index += 1

        topology = {key: pd.DataFrame([]) for key, values in topo_labels.items() if len(values) > 0}
        for key in topology:
            df_topology = pd.DataFrame(np.concatenate(topo_collector[key]), columns=SECTION_HEADERS[key][1:])
            df_topology["type"] = list(map(ff.maps[key].get, topo_labels[key]))
            if any(pd.isna(df_topology["type"])):  # Throw away undefined topologies
                warnings.warn(f"Undefined {key.lower()} detected and removed")
                df_topology = df_topology.dropna(subset=["type"])
                df_topology = df_topology.reset_index(drop=True)
            df_topology.index += 1
            topology[key] = df_topology[SECTION_HEADERS[key]]
        topology = {key: values for key, values in topology.items() if not values.empty}

        items |= {"atoms": atoms, "velocities": velocities, "topology": topology}
        return cls(**items)

    @classmethod
    def from_structure(
        cls,
        structure: Structure,
        ff_elements: Sequence[str] | None = None,
        atom_style: Literal["atomic", "charge"] = "charge",
        is_sort: bool = False,
    ) -> Self:
        """Simple constructor building LammpsData from a structure without
        force field parameters and topologies.

        Args:
            structure (Structure): Input structure.
            ff_elements ([str]): List of strings of elements that must
                be present due to force field settings but not
                necessarily in the structure. Default to None.
            atom_style (str): Choose between "atomic" (neutral) and
                "charge" (charged). Default to "charge".
            is_sort (bool): whether to sort sites
        """
        struct = structure.get_sorted_structure() if is_sort else structure.copy()
        box, symm_op = lattice_2_lmpbox(struct.lattice)
        coords = symm_op.operate_multi(struct.cart_coords)
        site_properties = struct.site_properties
        if "velocities" in site_properties:
            velos = np.array(struct.site_properties["velocities"])
            rot = SymmOp.from_rotation_and_translation(symm_op.rotation_matrix)
            rot_velos = rot.operate_multi(velos)
            site_properties["velocities"] = rot_velos
        boxed_s = Structure(
            box.to_lattice(),
            struct.species,
            coords,
            site_properties=site_properties,
            coords_are_cartesian=True,
        )

        symbols = list(struct.symbol_set)
        if ff_elements:
            symbols.extend(ff_elements)
        elements = sorted(Element(el) for el in set(symbols))
        mass_info = [tuple([i.symbol] * 2) for i in elements]
        ff = ForceField(mass_info)
        topo = Topology(boxed_s)
        return cls.from_ff_and_topologies(box=box, ff=ff, topologies=[topo], atom_style=atom_style)

    def set_charge_atom(self, charges: dict[int, float]) -> None:
        """Set the charges of specific atoms of the data.

        Args:
            charges: A dictionary with atom indexes as keys and
                charges as values, e.g. to set the charge
                of the atom with index 3 to -2, use `{3: -2}`.
        """
        for iat, q in charges.items():
            self.atoms.loc[iat, "q"] = q

    def set_charge_atom_type(self, charges: dict[str | int, float]) -> None:
        """
        Add or modify charges of all atoms of a given type in the data.

        Args:
            charges: Dict containing the charges for the atom types to set.
                The dict should contain atom types as integers or labels and charges.
                Example: change the charge of Li atoms to +3:
                    charges={"Li": 3}
                    charges={1: 3} if Li atoms are of type 1
        """
        for iat, q in charges.items():
            if isinstance(iat, str):
                mass_iat = Element(iat).atomic_mass
                iat = self.masses.loc[self.masses["mass"] == mass_iat].index[0]
            self.atoms.loc[self.atoms["type"] == iat, "q"] = q


class Topology(MSONable):
    """Carry most data in Atoms, Velocities and Molecular Topology sections for
    ONE SINGLE Molecule or Structure object, or a plain list of Sites.
    """

    def __init__(
        self,
        sites: Sequence[Site] | SiteCollection,
        ff_label: str | None = None,
        charges: Sequence | None = None,
        velocities: Sequence[Sequence] | None = None,
        topologies: dict | None = None,
    ) -> None:
        """
        Args:
            sites ([Site] or SiteCollection): A group of sites in a
                list or as a Molecule/Structure.
            ff_label (str): Site property key for labeling atoms of
                different types. Default to None, i.e., use
                site.species_string.
            charges ([q, ...]): Charge of each site in a (n,)
                array/list, where n is the No. of sites. Default to
                None, i.e., search site property for charges.
            velocities ([[vx, vy, vz], ...]): Velocity of each site in
                a (n, 3) array/list, where n is the No. of sites.
                Default to None, i.e., search site property for
                velocities.
            topologies (dict): Bonds, angles, dihedrals and improper
                dihedrals defined by site indices. Default to None,
                i.e., no additional topology. All four valid keys
                listed below are optional.
                {
                    "Bonds": [[i, j], ...],
                    "Angles": [[i, j, k], ...],
                    "Dihedrals": [[i, j, k, l], ...],
                    "Impropers": [[i, j, k, l], ...]
                }.
        """
        if not isinstance(sites, Molecule | Structure):
            sites = Molecule.from_sites(sites)

        type_by_sites = sites.site_properties[ff_label] if ff_label else [site.specie.symbol for site in sites]
        # search for site property if not override
        if charges is None:
            charges = sites.site_properties.get("charge")
        if velocities is None:
            velocities = sites.site_properties.get("velocities")
        # validate shape
        if charges is not None:
            charge_arr = np.array(charges)
            if charge_arr.shape != (len(sites),):
                raise ValueError(f"{charge_arr.shape=} and {(len(sites), )=} mismatch")
            charges = charge_arr.tolist()
        if velocities is not None:
            velocities_arr = np.array(velocities)
            if velocities_arr.shape != (len(sites), 3):
                raise ValueError(f"{velocities_arr.shape=} and {(len(sites), 3)=} mismatch")
            velocities = velocities_arr.tolist()

        if topologies:
            topologies = {k: v for k, v in topologies.items() if k in SECTION_KEYWORDS["topology"]}

        self.sites = sites
        self.ff_label = ff_label
        self.charges = charges
        self.velocities = velocities
        self.topologies = topologies
        self.type_by_sites = type_by_sites
        self.species = set(type_by_sites)

    @classmethod
    def from_bonding(
        cls,
        molecule: Molecule,
        bond: bool = True,
        angle: bool = True,
        dihedral: bool = True,
        tol: float = 0.1,
        **kwargs,
    ) -> Self:
        """
        Another constructor that creates an instance from a molecule.
        Covalent bonds and other bond-based topologies (angles and
        dihedrals) can be automatically determined. Cannot be used for
        non bond-based topologies, e.g. improper dihedrals.

        Args:
            molecule (Molecule): Input molecule.
            bond (bool): Whether find bonds. If set to False, angle and
                dihedral searching will be skipped. Default to True.
            angle (bool): Whether find angles. Default to True.
            dihedral (bool): Whether find dihedrals. Default to True.
            tol (float): Bond distance tolerance. Default to 0.1.
                Not recommended to alter.
            **kwargs: Other kwargs supported by Topology.
        """
        real_bonds = molecule.get_covalent_bonds(tol=tol)
        bond_list = [list(map(molecule.index, [b.site1, b.site2])) for b in real_bonds]
        if not all((bond, bond_list)):
            # do not search for others if not searching for bonds or no bonds
            return cls(sites=molecule, **kwargs)

        angle_list, dihedral_list = [], []
        dests, freq = np.unique(bond_list, return_counts=True)
        hubs = dests[np.where(freq > 1)].tolist()
        bond_arr = np.array(bond_list)
        hub_spokes = {}
        if len(hubs) > 0:
            for hub in hubs:
                ix = np.any(np.isin(bond_arr, hub), axis=1)
                bonds = np.unique(bond_arr[ix]).tolist()
                bonds.remove(hub)
                hub_spokes[hub] = bonds
        # skip angle or dihedral searching if too few bonds or hubs
        dihedral = False if len(bond_list) < 3 or len(hubs) < 2 else dihedral
        angle = False if len(bond_list) < 2 or len(hubs) < 1 else angle

        if angle:
            for k, v in hub_spokes.items():
                angle_list.extend([[i, k, j] for i, j in itertools.combinations(v, 2)])
        if dihedral:
            hub_cons = bond_arr[np.all(np.isin(bond_arr, hubs), axis=1)]
            for ii, jj in hub_cons.tolist():
                ks = [ki for ki in hub_spokes[ii] if ki != jj]
                ls = [li for li in hub_spokes[jj] if li != ii]
                dihedral_list.extend([[ki, ii, jj, li] for ki, li in itertools.product(ks, ls) if ki != li])

        topologies = {
            k: v
            for k, v in zip(SECTION_KEYWORDS["topology"][:3], [bond_list, angle_list, dihedral_list], strict=True)
            if len(v) > 0
        } or None
        return cls(sites=molecule, topologies=topologies, **kwargs)


class ForceField(MSONable):
    """
    Class carrying most data in masses and force field sections.

    Attributes:
        masses (pandas.DataFrame): DataFrame for masses section.
        force_fieldct (dict): Force field section keywords (keys) and
            data (values) as DataFrames.
        maps (dict): Dict for labeling atoms and topologies.
    """

    @staticmethod
    def _is_valid(df: pd.DataFrame):
        return not pd.isna(df).to_numpy().any()

    def __init__(self, mass_info: list, nonbond_coeffs: list | None = None, topo_coeffs: dict | None = None) -> None:
        """
        Args:
            mass_info (list): List of atomic mass info. Elements,
                strings (symbols) and floats are all acceptable for the
                values, with the first two converted to the atomic mass
                of an element. It is recommended to use
                dict.items() to prevent key duplications.
                [("C", 12.01), ("H", Element("H")), ("O", "O"), ...]
            nonbond_coeffs (list): List of Pair or PairIJ
                coefficients, of which the sequence must be sorted
                according to the species in mass_dict. Pair or PairIJ
                determined by the length of list. Optional with default
                to None.
            topo_coeffs (dict): Dict with force field coefficients for
                molecular topologies. Optional with default
                to None. All four valid keys listed below are optional.
                Each value is a list of dicts with non-optional keys
                "coeffs" and "types", and related class2 force field
                keywords as optional keys.
                {
                    "Bond Coeffs":
                        [{"coeffs": [coeff],
                          "types": [("C", "C"), ...]}, ...],
                    "Angle Coeffs":
                        [{"coeffs": [coeff],
                          "BondBond Coeffs": [coeff],
                          "types": [("H", "C", "H"), ...]}, ...],
                    "Dihedral Coeffs":
                        [{"coeffs": [coeff],
                          "BondBond13 Coeffs": [coeff],
                          "types": [("H", "C", "C", "H"), ...]}, ...],
                    "Improper Coeffs":
                        [{"coeffs": [coeff],
                          "AngleAngle Coeffs": [coeff],
                          "types": [("H", "C", "C", "H"), ...]}, ...],
                }
                Topology of same type or equivalent types (e.g.,
                ("C", "H") and ("H", "C") bonds) are NOT ALLOWED to
                be defined MORE THAN ONCE with DIFFERENT coefficients.
        """

        def map_mass(v):
            return (
                v.atomic_mass.real
                if isinstance(v, Element)
                else Element(v).atomic_mass.real
                if isinstance(v, str)
                else v
            )

        index, masses, self.mass_info, atoms_map = [], [], [], {}
        for idx, m in enumerate(mass_info, start=1):
            index.append(idx)
            mass = map_mass(m[1])
            masses.append(mass)
            self.mass_info.append((m[0], mass))
            atoms_map[m[0]] = idx
        self.masses = pd.DataFrame({"mass": masses}, index=index)
        self.maps = {"Atoms": atoms_map}

        ff_dfs = {}

        self.nonbond_coeffs = nonbond_coeffs
        if self.nonbond_coeffs:
            ff_dfs.update(self._process_nonbond())

        self.topo_coeffs = topo_coeffs
        if self.topo_coeffs:
            self.topo_coeffs = {
                key: values for key, values in self.topo_coeffs.items() if key in SECTION_KEYWORDS["ff"][2:]
            }
            for k in self.topo_coeffs:
                coeffs, mapper = self._process_topo(k, self.topo_coeffs)
                ff_dfs.update(coeffs)
                self.maps.update(mapper)

        self.force_field = ff_dfs or None

    def _process_nonbond(self) -> dict:
        pair_df = pd.DataFrame(self.nonbond_coeffs)
        if not self._is_valid(pair_df):
            raise ValueError("Invalid nonbond coefficients with rows varying in length")
        n_pair, n_coeff = pair_df.shape
        pair_df.columns = [f"coeff{i}" for i in range(1, n_coeff + 1)]
        n_mass = len(self.mass_info)
        n_comb = int(n_mass * (n_mass + 1) / 2)
        if n_pair == n_mass:
            kw = "Pair Coeffs"
            pair_df.index = range(1, n_mass + 1)
        elif n_pair == n_comb:
            kw = "PairIJ Coeffs"
            ids = list(itertools.combinations_with_replacement(range(1, n_mass + 1), 2))
            id_df = pd.DataFrame(ids, columns=["id1", "id2"])
            pair_df = pd.concat([id_df, pair_df], axis=1)
        else:
            raise ValueError(
                f"Expecting {n_mass} Pair Coeffs or {n_comb} PairIJ Coeffs for {n_mass} atom types, got {n_pair}"
            )
        return {kw: pair_df}

    def _process_topo(self, kw: str, topo_coeffs: dict) -> tuple[dict, dict]:
        def find_eq_types(label, section) -> list:
            if section.startswith("Improper"):
                label_arr = np.array(label)
                seqs = [[0, 1, 2, 3], [0, 2, 1, 3], [3, 1, 2, 0], [3, 2, 1, 0]]
                return [tuple(label_arr[s]) for s in seqs]
            return [label, label[::-1]]

        main_data, distinct_types = [], []
        class2_data: dict = {key: [] for key in topo_coeffs[kw][0] if key in CLASS2_KEYWORDS.get(kw, [])}
        for dct in topo_coeffs[kw]:
            main_data.append(dct["coeffs"])
            distinct_types.append(dct["types"])
            for key, lst in class2_data.items():
                lst.append(dct[key])
        distinct_types = [set(itertools.chain(*(find_eq_types(t, kw) for t in dt))) for dt in distinct_types]
        type_counts = sum(len(dt) for dt in distinct_types)
        type_union = set.union(*distinct_types)
        if len(type_union) != type_counts:
            raise ValueError(f"Duplicated items found under different coefficients in {kw}")
        atoms = set(np.ravel(list(itertools.chain(*distinct_types))))
        if not atoms.issubset(self.maps["Atoms"]):
            raise ValueError(f"Undefined atom type found in {kw}")
        mapper = {}
        for i, dt in enumerate(distinct_types, start=1):
            for t in dt:
                mapper[t] = i

        def process_data(data) -> pd.DataFrame:
            df_coeffs = pd.DataFrame(data)
            if not self._is_valid(df_coeffs):
                raise ValueError("Invalid coefficients with rows varying in length")
            n, c = df_coeffs.shape
            df_coeffs.columns = [f"coeff{i}" for i in range(1, c + 1)]
            df_coeffs.index = range(1, n + 1)
            return df_coeffs

        all_data = {kw: process_data(main_data)}
        if class2_data:
            all_data |= {k: process_data(v) for k, v in class2_data.items()}
        return all_data, {f"{kw[:-7]}s": mapper}

    def to_file(self, filename: str) -> None:
        """Save force field to a file in YAML format.

        Args:
            filename (str): Filename.
        """
        dct = {
            "mass_info": self.mass_info,
            "nonbond_coeffs": self.nonbond_coeffs,
            "topo_coeffs": self.topo_coeffs,
        }
        with open(filename, mode="w", encoding="utf-8") as file:
            yaml = YAML()
            yaml.dump(dct, file)

    @classmethod
    def from_file(cls, filename: str) -> Self:
        """Constructor that reads in a file in YAML format.

        Args:
            filename (str): Filename.
        """
        with open(filename, encoding="utf-8") as file:
            yaml = YAML()
            d = yaml.load(file)
        return cls.from_dict(d)

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Constructor that reads in a dictionary.

        Args:
            dct (dict): Dictionary to read.
        """
        dct["mass_info"] = [tuple(m) for m in dct["mass_info"]]
        if dct.get("topo_coeffs"):
            for v in dct["topo_coeffs"].values():
                for c in v:
                    c["types"] = [tuple(t) for t in c["types"]]
        return cls(dct["mass_info"], dct["nonbond_coeffs"], dct["topo_coeffs"])


class CombinedData(LammpsData):
    """
    Object for a collective set of data for a series of LAMMPS data file.
    velocities not yet implemented.
    """

    def __init__(
        self,
        list_of_molecules: list,
        list_of_names: list[str],
        list_of_numbers: list[int],
        coordinates: pd.DataFrame,
        atom_style: str = "full",
    ) -> None:
        """
        Args:
            list_of_molecules: A list of LammpsData objects of a chemical cluster.
                Each LammpsData object (cluster) may contain one or more molecule ID.
            list_of_names: A list of name (string) for each cluster. The characters in each name are
                restricted to word characters ([a-zA-Z0-9_]). If names with any non-word characters
                are passed in, the special characters will be substituted by '_'.
            list_of_numbers: A list of Integer for counts of each molecule
            coordinates (pandas.DataFrame): DataFrame at least containing
                columns of ["x", "y", "z"] for coordinates of atoms.
            atom_style (str): Output atom_style. Default to "full".
        """
        self._list_of_molecules = list_of_molecules
        self._list_of_names = list_of_names
        self._list_of_numbers = list_of_numbers
        self._coordinates = coordinates
        self._coordinates.index = self._coordinates.index.map(int)
        max_xyz = self._coordinates[["x", "y", "z"]].max().max()
        min_xyz = self._coordinates[["x", "y", "z"]].min().min()
        self.box = LammpsBox(np.array(3 * [[min_xyz - 0.5, max_xyz + 0.5]]))
        self.atom_style = atom_style
        self.n = sum(self._list_of_numbers)
        self.names = []
        for name in self._list_of_names:
            self.names.append("_".join(re.findall(r"\w+", name)))
        self.mols = self._list_of_molecules
        self.nums = self._list_of_numbers
        self.masses = pd.concat([mol.masses.copy() for mol in self.mols], ignore_index=True)
        self.masses.index += 1
        all_ff_kws = SECTION_KEYWORDS["ff"] + SECTION_KEYWORDS["class2"]
        appeared_kws = {k for mol in self.mols if mol.force_field is not None for k in mol.force_field}
        ff_kws = [k for k in all_ff_kws if k in appeared_kws]
        self.force_field = {}
        for kw in ff_kws:
            self.force_field[kw] = pd.concat(
                [mol.force_field[kw].copy() for mol in self.mols if kw in (mol.force_field or [])],
                ignore_index=True,
            )
            self.force_field[kw].index += 1
        if not bool(self.force_field):
            self.force_field = None

        self.atoms = pd.DataFrame()
        mol_count = type_count = 0
        self.mols_per_data = []
        for idx, mol in enumerate(self.mols):
            atoms_df = mol.atoms.copy()
            atoms_df["molecule-ID"] += mol_count
            atoms_df["type"] += type_count
            mols_in_data = len(atoms_df["molecule-ID"].unique())
            self.mols_per_data.append(mols_in_data)
            for _ in range(self.nums[idx]):
                self.atoms = pd.concat([self.atoms, atoms_df], ignore_index=True)
                atoms_df["molecule-ID"] += mols_in_data
            type_count += len(mol.masses)
            mol_count += self.nums[idx] * mols_in_data
        self.atoms.index += 1
        if len(self.atoms) != len(self._coordinates):
            raise ValueError(f"{len(self.atoms)=} and {len(self._coordinates)=} mismatch")
        self.atoms.update(self._coordinates)

        self.velocities = None
        if self.mols[0].velocities is not None:
            raise RuntimeError("Velocities not supported")

        self.topology = {}
        atom_count = 0
        count = {"Bonds": 0, "Angles": 0, "Dihedrals": 0, "Impropers": 0}
        for idx, mol in enumerate(self.mols):
            for kw in SECTION_KEYWORDS["topology"]:
                if mol.topology and kw in mol.topology:
                    if kw not in self.topology:
                        self.topology[kw] = pd.DataFrame()
                    topo_df = mol.topology[kw].copy()
                    topo_df["type"] += count[kw]
                    for col in topo_df.columns[1:]:
                        topo_df[col] += atom_count
                    for _ in range(self.nums[idx]):
                        self.topology[kw] = pd.concat([self.topology[kw], topo_df], ignore_index=True)
                        for col in topo_df.columns[1:]:
                            topo_df[col] += len(mol.atoms)
                    count[kw] += len(mol.force_field[kw[:-1] + " Coeffs"])
            atom_count += len(mol.atoms) * self.nums[idx]
        for kw in SECTION_KEYWORDS["topology"]:
            if kw in self.topology:
                self.topology[kw].index += 1
        if not self.topology:
            self.topology = None

    @property
    def structure(self) -> Structure:
        """Exports a periodic structure object representing the simulation
        box.

        Returns:
            Structure
        """
        ld_cp = self.as_lammpsdata()
        return ld_cp.structure

    def disassemble(
        self, atom_labels: Sequence[str] | None = None, guess_element: bool = True, ff_label: str = "ff_map"
    ):
        """
        Breaks down each LammpsData in CombinedData to building blocks
        (LammpsBox, ForceField and a series of Topology).
        RESTRICTIONS APPLIED:
        1. No complex force field defined not just on atom
            types, where the same type or equivalent types of topology
            may have more than one set of coefficients.
        2. No intermolecular topologies (with atoms from different
            molecule-ID) since a Topology object includes data for ONE
            molecule or structure only.

        Args:
            atom_labels ([str]): List of strings (must be different
                from one another) for labelling each atom type found in
                Masses section. Default to None, where the labels are
                automatically added based on either element guess or
                dummy specie assignment.
            guess_element (bool): Whether to guess the element based on
                its atomic mass. Default to True, otherwise dummy
                species "Qa", "Qb", ... will be assigned to various
                atom types. The guessed or assigned elements will be
                reflected on atom labels if atom_labels is None, as
                well as on the species of molecule in each Topology.
            ff_label (str): Site property key for labeling atoms of
                different types. Default to "ff_map".

        Returns:
            [(LammpsBox, ForceField, [Topology]), ...]
        """
        return [
            mol.disassemble(atom_labels=atom_labels, guess_element=guess_element, ff_label=ff_label)
            for mol in self.mols
        ]

    # NOTE (@janosh): The following two methods for override parent class LammpsData
    @classmethod
    def from_ff_and_topologies(cls) -> None:
        """Unsupported constructor for CombinedData objects."""
        raise AttributeError("Unsupported constructor for CombinedData objects")

    @classmethod
    def from_structure(cls) -> None:
        """Unsupported constructor for CombinedData objects."""
        raise AttributeError("Unsupported constructor for CombinedData objects")

    @classmethod
    def parse_xyz(cls, filename: str | Path) -> pd.DataFrame:
        """
        Load xyz file generated from packmol (for those who find it hard to install openbabel).

        Returns:
            pandas.DataFrame
        """
        with zopen(filename, mode="rt") as file:
            lines = file.readlines()

        str_io = StringIO("".join(lines[2:]))  # skip the 2nd line
        df_xyz = pd.read_csv(
            str_io,
            header=None,
            comment="#",
            sep=r"\s+",
            names=["atom", "x", "y", "z"],
        )
        df_xyz.index += 1
        return df_xyz

    @classmethod
    def from_files(cls, coordinate_file: str, list_of_numbers: list[int], *filenames: str) -> Self:
        """
        Constructor that parse a series of data file.

        Args:
            coordinate_file (str): The filename of xyz coordinates.
            list_of_numbers (list[int]): A list of numbers specifying counts for each
                clusters parsed from files.
            filenames (str): A series of LAMMPS data filenames in string format.
        """
        names = []
        mols = []
        styles = []
        clusters = []

        for idx, filename in enumerate(filenames, start=1):
            cluster = LammpsData.from_file(filename)
            clusters.append(cluster)
            names.append(f"cluster{idx}")
            mols.append(cluster)
            styles.append(cluster.atom_style)

        if len(set(styles)) != 1:
            raise ValueError("Files have different atom styles.")

        coordinates = cls.parse_xyz(filename=coordinate_file)
        return cls.from_lammpsdata(mols, names, list_of_numbers, coordinates, styles.pop())

    @classmethod
    def from_lammpsdata(
        cls, mols: list, names: list, list_of_numbers: list, coordinates: pd.DataFrame, atom_style: str | None = None
    ) -> Self:
        """
        Constructor that can infer atom_style.
        The input LammpsData objects are used non-destructively.

        Args:
            mols: a list of LammpsData of a chemical cluster.Each LammpsData object (cluster)
                may contain one or more molecule ID.
            names: a list of name for each cluster.
            list_of_numbers: a list of Integer for counts of each molecule
            coordinates (pandas.DataFrame): DataFrame at least containing
                columns of ["x", "y", "z"] for coordinates of atoms.
            atom_style (str): Output atom_style. Default to "full".
        """
        styles = [mol.atom_style for mol in mols]

        if len(set(styles)) != 1:
            raise ValueError("Data have different atom_style.")
        style_return = styles.pop()

        if atom_style and atom_style != style_return:
            raise ValueError("Data have different atom_style as specified.")

        return cls(mols, names, list_of_numbers, coordinates, style_return)

    def get_str(self, distance: int = 6, velocity: int = 8, charge: int = 4, hybrid: bool = True) -> str:
        """Get the string representation of CombinedData, essentially
        the string to be written to a file. Combination info is included
        as a comment. For single molecule ID data, the info format is:
            num name
        For data with multiple molecule ID, the format is:
            num(mols_per_data) name.

        Args:
            distance (int): No. of significant figures to output for
                box settings (bounds and tilt) and atomic coordinates.
                Default to 6.
            velocity (int): No. of significant figures to output for
                velocities. Default to 8.
            charge (int): No. of significant figures to output for
                charges. Default to 4.
            hybrid (bool): Whether to write hybrid coeffs types.
                Default to True. If the data object has no hybrid
                coeffs types and has large coeffs section, one may
                use False to speed up the process. Otherwise, the
                default is recommended.

        Returns:
            str: String representation of CombinedData.
        """
        lines = LammpsData.get_str(self, distance, velocity, charge, hybrid).splitlines()
        info = "# " + " + ".join(
            f"{a} {b}" if c == 1 else f"{a}({c}) {b}"
            for a, b, c in zip(self.nums, self.names, self.mols_per_data, strict=True)
        )
        lines.insert(1, info)
        return "\n".join(lines)

    def as_lammpsdata(self):
        """
        Convert a CombinedData object to a LammpsData object. attributes are deep-copied.

        box (LammpsBox): Simulation box.
        force_fieldct (dict): Data for force field sections. Optional
            with default to None. Only keywords in force field and
            class 2 force field are valid keys, and each value is a
            DataFrame.
        topology (dict): Data for topology sections. Optional with
            default to None. Only keywords in topology are valid
            keys, and each value is a DataFrame.
        """
        items = {}
        items["box"] = LammpsBox(self.box.bounds, self.box.tilt)
        items["masses"] = self.masses.copy()
        items["atoms"] = self.atoms.copy()
        items["atom_style"] = self.atom_style
        items["velocities"] = None  # Velocities not supported
        if self.force_field:
            all_ff_kws = SECTION_KEYWORDS["ff"] + SECTION_KEYWORDS["class2"]
            items["force_field"] = {k: v.copy() for k, v in self.force_field.items() if k in all_ff_kws}

        if self.topology:
            items["topology"] = {k: v.copy() for k, v in self.topology.items() if k in SECTION_KEYWORDS["topology"]}
        return LammpsData(**items)
