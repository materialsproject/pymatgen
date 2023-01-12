# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

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

import numpy as np
import pandas as pd
from monty.io import zopen
from monty.json import MSONable
from monty.serialization import loadfn
from ruamel.yaml import YAML

from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Molecule, Structure
from pymatgen.util.io_utils import clean_lines

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
    """
    Object for representing a simulation box in LAMMPS settings.
    """

    def __init__(self, bounds, tilt=None):
        """
        Args:
            bounds: A (3, 2) array/list of floats setting the
                boundaries of simulation box.
            tilt: A (3,) array/list of floats setting the tilt of
                simulation box. Default to None, i.e., use an
                orthogonal box.
        """
        bounds_arr = np.array(bounds)
        assert bounds_arr.shape == (
            3,
            2,
        ), f"Expecting a (3, 2) array for bounds, got {bounds_arr.shape}"
        self.bounds = bounds_arr.tolist()
        matrix = np.diag(bounds_arr[:, 1] - bounds_arr[:, 0])

        self.tilt = None
        if tilt is not None:
            tilt_arr = np.array(tilt)
            assert tilt_arr.shape == (3,), f"Expecting a (3,) array for box_tilt, got {tilt_arr.shape}"
            self.tilt = tilt_arr.tolist()
            matrix[1, 0] = tilt_arr[0]
            matrix[2, 0] = tilt_arr[1]
            matrix[2, 1] = tilt_arr[2]
        self._matrix = matrix

    def __str__(self):
        return self.get_string()

    def __repr__(self):
        return self.get_string()

    @property
    def volume(self):
        """
        Volume of simulation box.
        """
        m = self._matrix
        return np.dot(np.cross(m[0], m[1]), m[2])

    def get_string(self, significant_figures=6):
        """
        Returns the string representation of simulation box in LAMMPS
        data file format.

        Args:
            significant_figures (int): No. of significant figures to
                output for box settings. Default to 6.

        Returns:
            String representation
        """
        ph = f"{{:.{significant_figures}f}}"
        lines = []
        for bound, d in zip(self.bounds, "xyz"):
            fillers = bound + [d] * 2
            bound_format = " ".join([ph] * 2 + [" {}lo {}hi"])
            lines.append(bound_format.format(*fillers))
        if self.tilt:
            tilt_format = " ".join([ph] * 3 + [" xy xz yz"])
            lines.append(tilt_format.format(*self.tilt))
        return "\n".join(lines)

    def get_box_shift(self, i):
        """
        Calculates the coordinate shift due to PBC.

        Args:
            i: A (n, 3) integer array containing the labels for box
            images of n entries.

        Returns:
            Coorindate shift array with the same shape of i
        """
        return np.inner(i, self._matrix.T)

    def to_lattice(self):
        """
        Converts the simulation box to a more powerful Lattice backend.
        Note that Lattice is always periodic in 3D space while a
        simulation box is not necessarily periodic in all dimensions.

        Returns:
            Lattice
        """
        return Lattice(self._matrix)


def lattice_2_lmpbox(lattice, origin=(0, 0, 0)):
    """
    Converts a lattice object to LammpsBox, and calculates the symmetry
    operation used.

    Args:
        lattice (Lattice): Input lattice.
        origin: A (3,) array/list of floats setting lower bounds of
            simulation box. Default to (0, 0, 0).

    Returns:
        LammpsBox, SymmOp

    """
    a, b, c = lattice.abc
    xlo, ylo, zlo = origin
    xhi = a + xlo
    m = lattice.matrix
    xy = np.dot(m[1], m[0] / a)
    yhi = np.sqrt(b**2 - xy**2) + ylo
    xz = np.dot(m[2], m[0] / a)
    yz = (np.dot(m[1], m[2]) - xy * xz) / (yhi - ylo)
    zhi = np.sqrt(c**2 - xz**2 - yz**2) + zlo
    tilt = None if lattice.is_orthogonal else [xy, xz, yz]
    rot_matrix = np.linalg.solve([[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, zhi - zlo]], m)
    bounds = [[xlo, xhi], [ylo, yhi], [zlo, zhi]]
    symm_op = SymmOp.from_rotation_and_translation(rot_matrix, origin)
    return LammpsBox(bounds, tilt), symm_op


class LammpsData(MSONable):
    """
    Object for representing the data in a LAMMPS data file.

    """

    def __init__(
        self,
        box,
        masses,
        atoms,
        velocities=None,
        force_field=None,
        topology=None,
        atom_style="full",
    ):
        """
        This is a low level constructor designed to work with parsed
        data or other bridging objects (ForceField and Topology). Not
        recommended to use directly.

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
            force_field (dict): Data for force field sections. Optional
                with default to None. Only keywords in force field and
                class 2 force field are valid keys, and each value is a
                DataFrame.
            topology (dict): Data for topology sections. Optional with
                default to None. Only keywords in topology are valid
                keys, and each value is a DataFrame.
            atom_style (str): Output atom_style. Default to "full".
        """
        if velocities is not None:
            assert len(velocities) == len(atoms), "Inconsistency found between atoms and velocities"

        if force_field:
            all_ff_kws = SECTION_KEYWORDS["ff"] + SECTION_KEYWORDS["class2"]
            force_field = {k: v for k, v in force_field.items() if k in all_ff_kws}

        if topology:
            topology = {k: v for k, v in topology.items() if k in SECTION_KEYWORDS["topology"]}

        self.box = box
        self.masses = masses
        self.atoms = atoms
        self.velocities = velocities
        self.force_field = force_field
        self.topology = topology
        self.atom_style = atom_style

    def __str__(self):
        return self.get_string()

    def __repr__(self):
        return self.get_string()

    @property
    def structure(self):
        """
        Exports a periodic structure object representing the simulation
        box.

        Return:
            Structure
        """
        masses = self.masses
        atoms = self.atoms.copy()
        if "nx" in atoms.columns:
            atoms.drop(["nx", "ny", "nz"], axis=1, inplace=True)
        atoms["molecule-ID"] = 1
        ld_copy = self.__class__(self.box, masses, atoms)
        topologies = ld_copy.disassemble()[-1]
        molecule = topologies[0].sites
        coords = molecule.cart_coords - np.array(self.box.bounds)[:, 0]
        species = molecule.species
        latt = self.box.to_lattice()
        site_properties = {}
        if "q" in atoms:
            site_properties["charge"] = atoms["q"].values
        if self.velocities is not None:
            site_properties["velocities"] = self.velocities.values
        return Structure(
            latt,
            species,
            coords,
            coords_are_cartesian=True,
            site_properties=site_properties,
        )

    def get_string(self, distance=6, velocity=8, charge=4, hybrid=True):
        """
        Returns the string representation of LammpsData, essentially
        the string to be written to a file. Support hybrid style
        coeffs read and write.

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
            String representation
        """
        file_template = """Generated by pymatgen.io.lammps.data.LammpsData

{stats}

{box}

{body}
"""
        box = self.box.get_string(distance)

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
        type_lines = [f"{v:>{right_indent}}  {k+ ' types'}" for k, v in types.items()]
        stats = "\n".join(count_lines + [""] + type_lines)

        def map_coords(q):
            return f"{q:.{distance}f}"

        def map_velos(q):
            return f"{q:.{velocity}f}"

        def map_charges(q):
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
        coeffsdatatype = loadfn(str(MODULE_DIR / "CoeffsDataType.yaml"))
        coeffs = {}
        for style, types in coeffsdatatype.items():
            coeffs[style] = {}
            for type, formatter in types.items():
                coeffs[style][type] = {}
                for coeff, datatype in formatter.items():
                    if datatype == "int_format":
                        coeffs[style][type][coeff] = int_format
                    elif datatype == "float_format_2":
                        coeffs[style][type][coeff] = float_format_2
                    else:
                        coeffs[style][type][coeff] = float_format

        section_template = "{kw}\n\n{df}\n"
        parts = []
        for k, v in body_dict.items():
            index = k != "PairIJ Coeffs"
            if (
                k
                in [
                    "Bond Coeffs",
                    "Angle Coeffs",
                    "Dihedral Coeffs",
                    "Improper Coeffs",
                ]
                and hybrid
            ):
                listofdf = np.array_split(v, len(v.index))
                df_string = ""
                for i, df in enumerate(listofdf):
                    if isinstance(df.iloc[0]["coeff1"], str):
                        try:
                            formatters = {
                                **default_formatters,
                                **coeffs[k][df.iloc[0]["coeff1"]],
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
                        line_string = v.to_string(
                            header=False,
                            formatters=default_formatters,
                            index_names=False,
                            index=index,
                            na_rep="",
                        ).splitlines()[i]
                    df_string += line_string.replace("nan", "").rstrip() + "\n"
            else:
                df_string = v.to_string(
                    header=False,
                    formatters=default_formatters,
                    index_names=False,
                    index=index,
                    na_rep="",
                )
            parts.append(section_template.format(kw=k, df=df_string))
        body = "\n".join(parts)

        return file_template.format(stats=stats, box=box, body=body)

    def write_file(self, filename, distance=6, velocity=8, charge=4):
        """
        Writes LammpsData to file.

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
        with open(filename, "w") as f:
            f.write(self.get_string(distance=distance, velocity=velocity, charge=charge))

    def disassemble(self, atom_labels=None, guess_element=True, ff_label="ff_map"):
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
                df = atoms_df[atoms_df["molecule-ID"] == k]
                data_by_mols[k] = {"Atoms": df}

        masses = self.masses.copy()
        masses["label"] = atom_labels
        unique_masses = np.unique(masses["mass"])
        if guess_element:
            ref_masses = [el.atomic_mass.real for el in Element]
            diff = np.abs(np.array(ref_masses) - unique_masses[:, None])
            atomic_numbers = np.argmin(diff, axis=1) + 1
            symbols = [Element.from_Z(an).symbol for an in atomic_numbers]
        else:
            symbols = [f"Q{a}" for a in map(chr, range(97, 97 + len(unique_masses)))]
        for um, s in zip(unique_masses, symbols):
            masses.loc[masses["mass"] == um, "element"] = s
        if atom_labels is None:  # add unique labels based on elements
            for el, vc in masses["element"].value_counts().iteritems():
                masses.loc[masses["element"] == el, "label"] = [f"{el}{c}" for c in range(1, vc + 1)]
        assert masses["label"].nunique(dropna=False) == len(masses), "Expecting unique atom label for each type"
        mass_info = [(row.label, row.mass) for row in masses.itertuples()]

        nonbond_coeffs, topo_coeffs = None, None
        if self.force_field:
            if "PairIJ Coeffs" in self.force_field:
                nbc = self.force_field["PairIJ Coeffs"]
                nbc = nbc.sort_values(["id1", "id2"]).drop(["id1", "id2"], axis=1)
                nonbond_coeffs = [list(t) for t in nbc.itertuples(False, None)]
            elif "Pair Coeffs" in self.force_field:
                nbc = self.force_field["Pair Coeffs"].sort_index()
                nonbond_coeffs = [list(t) for t in nbc.itertuples(False, None)]

            topo_coeffs = {k: [] for k in SECTION_KEYWORDS["ff"][2:] if k in self.force_field}
            for kw in topo_coeffs:
                class2_coeffs = {
                    k: list(v.itertuples(False, None))
                    for k, v in self.force_field.items()
                    if k in CLASS2_KEYWORDS.get(kw, [])
                }
                ff_df = self.force_field[kw]
                for t in ff_df.itertuples(True, None):
                    d = {"coeffs": list(t[1:]), "types": []}
                    if class2_coeffs:
                        d.update({k: list(v[t[0] - 1]) for k, v in class2_coeffs.items()})
                    topo_coeffs[kw].append(d)

        if self.topology:

            def label_topo(t):
                return tuple(masses.loc[atoms_df.loc[t, "type"], "label"])

            for k, v in self.topology.items():
                ff_kw = k[:-1] + " Coeffs"
                for topo in v.itertuples(False, None):
                    topo_idx = topo[0] - 1
                    indices = list(topo[1:])
                    mids = atoms_df.loc[indices]["molecule-ID"].unique()
                    assert (
                        len(mids) == 1
                    ), "Do not support intermolecular topology formed by atoms with different molecule-IDs"
                    label = label_topo(indices)
                    topo_coeffs[ff_kw][topo_idx]["types"].append(label)
                    if data_by_mols[mids[0]].get(k):
                        data_by_mols[mids[0]][k].append(indices)
                    else:
                        data_by_mols[mids[0]][k] = [indices]

        if topo_coeffs:
            for v in topo_coeffs.values():
                for d in v:
                    d["types"] = list(set(d["types"]))

        ff = ForceField(mass_info=mass_info, nonbond_coeffs=nonbond_coeffs, topo_coeffs=topo_coeffs)

        topo_list = []
        for mid in unique_mids:
            data = data_by_mols[mid]
            atoms = data["Atoms"]
            shift = min(atoms.index)
            type_ids = atoms["type"]
            species = masses.loc[type_ids, "element"]
            labels = masses.loc[type_ids, "label"]
            coords = atoms[["x", "y", "z"]]
            m = Molecule(species.values, coords.values, site_properties={ff_label: labels.values})
            charges = atoms.get("q")
            velocities = atoms[["vx", "vy", "vz"]] if "vx" in atoms.columns else None
            topologies = {}
            for kw in SECTION_KEYWORDS["topology"]:
                if data.get(kw):
                    topologies[kw] = (np.array(data[kw]) - shift).tolist()
            topologies = None if not topologies else topologies
            topo_list.append(
                Topology(
                    sites=m,
                    ff_label=ff_label,
                    charges=charges,
                    velocities=velocities,
                    topologies=topologies,
                )
            )

        return self.box, ff, topo_list

    @classmethod
    def from_file(cls, filename, atom_style="full", sort_id=False):
        """
        Constructor that parses a file.

        Args:
            filename (str): Filename to read.
            atom_style (str): Associated atom_style. Default to "full".
            sort_id (bool): Whether sort each section by id. Default to
                True.
        """
        with zopen(filename, "rt") as f:
            lines = f.readlines()
        kw_pattern = r"|".join(itertools.chain(*SECTION_KEYWORDS.values()))
        section_marks = [i for i, l in enumerate(lines) if re.search(kw_pattern, l)]
        parts = np.split(lines, section_marks)

        float_group = r"([0-9eE.+-]+)"
        header_pattern = {}
        header_pattern["counts"] = r"^\s*(\d+)\s+([a-zA-Z]+)$"
        header_pattern["types"] = r"^\s*(\d+)\s+([a-zA-Z]+)\s+types$"
        header_pattern["bounds"] = r"^\s*{}$".format(r"\s+".join([float_group] * 2 + [r"([xyz])lo \3hi"]))
        header_pattern["tilt"] = r"^\s*{}$".format(r"\s+".join([float_group] * 3 + ["xy xz yz"]))

        header = {"counts": {}, "types": {}}
        bounds = {}
        for l in clean_lines(parts[0][1:]):  # skip the 1st line
            match = None
            for k, v in header_pattern.items():  # noqa: B007
                match = re.match(v, l)
                if match:
                    break
            if match and k in ["counts", "types"]:
                header[k][match.group(2)] = int(match.group(1))
            elif match and k == "bounds":
                g = match.groups()
                bounds[g[2]] = [float(i) for i in g[:2]]
            elif match and k == "tilt":
                header["tilt"] = [float(i) for i in match.groups()]
        header["bounds"] = [bounds.get(i, [-0.5, 0.5]) for i in "xyz"]
        box = LammpsBox(header["bounds"], header.get("tilt"))

        def parse_section(sec_lines):
            title_info = sec_lines[0].split("#", 1)
            kw = title_info[0].strip()
            sio = StringIO("".join(sec_lines[2:]))  # skip the 2nd line
            if kw.endswith("Coeffs") and not kw.startswith("PairIJ"):
                df_list = [
                    pd.read_csv(StringIO(line), header=None, comment="#", delim_whitespace=True)
                    for line in sec_lines[2:]
                    if line.strip()
                ]
                df = pd.concat(df_list, ignore_index=True)
                names = ["id"] + [f"coeff{i}" for i in range(1, df.shape[1])]
            else:
                df = pd.read_csv(sio, header=None, comment="#", delim_whitespace=True)
                if kw == "PairIJ Coeffs":
                    names = ["id1", "id2"] + [f"coeff{i}" for i in range(1, df.shape[1] - 1)]
                    df.index.name = None  # pylint: disable=E1101
                elif kw in SECTION_HEADERS:
                    names = ["id"] + SECTION_HEADERS[kw]
                elif kw == "Atoms":
                    names = ["id"] + ATOMS_HEADERS[atom_style]
                    if df.shape[1] == len(names):  # pylint: disable=E1101
                        pass
                    elif df.shape[1] == len(names) + 3:  # pylint: disable=E1101
                        names += ["nx", "ny", "nz"]
                    else:
                        raise ValueError(f"Format in Atoms section inconsistent with atom_style {atom_style}")
                else:
                    raise NotImplementedError(f"Parser for {kw} section not implemented")
            df.columns = names
            if sort_id:
                sort_by = "id" if kw != "PairIJ Coeffs" else ["id1", "id2"]
                df.sort_values(sort_by, inplace=True)
            if "id" in df.columns:
                df.set_index("id", drop=True, inplace=True)
                df.index.name = None
            return kw, df

        err_msg = "Bad LAMMPS data format where "
        body = {}
        seen_atoms = False
        for part in parts[1:]:
            name, section = parse_section(part)
            if name == "Atoms":
                seen_atoms = True
            if (
                name in ["Velocities"] + SECTION_KEYWORDS["topology"] and not seen_atoms
            ):  # Atoms must appear earlier than these
                raise RuntimeError(f"{err_msg}{name} section appears before Atoms section")
            body.update({name: section})

        err_msg += "Nos. of {} do not match between header and {} section"
        assert len(body["Masses"]) == header["types"]["atom"], err_msg.format("atom types", "Masses")
        atom_sections = ["Atoms", "Velocities"] if "Velocities" in body else ["Atoms"]
        for s in atom_sections:
            assert len(body[s]) == header["counts"]["atoms"], err_msg.format("atoms", s)
        for s in SECTION_KEYWORDS["topology"]:
            if header["counts"].get(s.lower(), 0) > 0:
                assert len(body[s]) == header["counts"][s.lower()], err_msg.format(s.lower(), s)

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
    def from_ff_and_topologies(cls, box, ff, topologies, atom_style="full"):
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
        assert atom_types.issubset(ff.maps["Atoms"]), "Unknown atom type found in topologies"

        items = dict(box=box, atom_style=atom_style, masses=ff.masses, force_field=ff.force_field)

        mol_ids, charges, coords, labels = [], [], [], []
        v_collector = [] if topologies[0].velocities else None
        topo_collector = {"Bonds": [], "Angles": [], "Dihedrals": [], "Impropers": []}
        topo_labels = {"Bonds": [], "Angles": [], "Dihedrals": [], "Impropers": []}
        for i, topo in enumerate(topologies):
            if topo.topologies:
                shift = len(labels)
                for k, v in topo.topologies.items():
                    topo_collector[k].append(np.array(v) + shift + 1)
                    topo_labels[k].extend([tuple(topo.type_by_sites[j] for j in t) for t in v])
            if isinstance(v_collector, list):
                v_collector.append(topo.velocities)
            mol_ids.extend([i + 1] * len(topo.sites))
            labels.extend(topo.type_by_sites)
            coords.append(topo.sites.cart_coords)
            q = [0.0] * len(topo.sites) if not topo.charges else topo.charges
            charges.extend(q)

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

        topology = {k: None for k, v in topo_labels.items() if len(v) > 0}
        for k in topology:
            df = pd.DataFrame(np.concatenate(topo_collector[k]), columns=SECTION_HEADERS[k][1:])
            df["type"] = list(map(ff.maps[k].get, topo_labels[k]))
            if any(pd.isnull(df["type"])):  # Throw away undefined topologies
                warnings.warn(f"Undefined {k.lower()} detected and removed")
                df.dropna(subset=["type"], inplace=True)
                df.reset_index(drop=True, inplace=True)
            df.index += 1
            topology[k] = df[SECTION_HEADERS[k]]
        topology = {k: v for k, v in topology.items() if not v.empty}

        items.update({"atoms": atoms, "velocities": velocities, "topology": topology})
        return cls(**items)

    @classmethod
    def from_structure(cls, structure: Structure, ff_elements=None, atom_style="charge", is_sort=False):
        """
        Simple constructor building LammpsData from a structure without
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
        if is_sort:
            s = structure.get_sorted_structure()
        else:
            s = structure.copy()
        box, symm_op = lattice_2_lmpbox(s.lattice)
        coords = symm_op.operate_multi(s.cart_coords)
        site_properties = s.site_properties
        if "velocities" in site_properties:
            velos = np.array(s.site_properties["velocities"])
            rot = SymmOp.from_rotation_and_translation(symm_op.rotation_matrix)
            rot_velos = rot.operate_multi(velos)
            site_properties.update({"velocities": rot_velos})  # type: ignore
        boxed_s = Structure(
            box.to_lattice(),
            s.species,
            coords,
            site_properties=site_properties,
            coords_are_cartesian=True,
        )

        symbols = list(s.symbol_set)
        if ff_elements:
            symbols.extend(ff_elements)
        elements = sorted(Element(el) for el in set(symbols))
        mass_info = [tuple([i.symbol] * 2) for i in elements]
        ff = ForceField(mass_info)
        topo = Topology(boxed_s)
        return cls.from_ff_and_topologies(box=box, ff=ff, topologies=[topo], atom_style=atom_style)


class Topology(MSONable):
    """
    Class carrying most data in Atoms, Velocities and molecular
    topology sections for ONE SINGLE Molecule or Structure
    object, or a plain list of Sites.

    """

    def __init__(self, sites, ff_label=None, charges=None, velocities=None, topologies=None):
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
                }
        """
        if not isinstance(sites, (Molecule, Structure)):
            sites = Molecule.from_sites(sites)

        if ff_label:
            type_by_sites = sites.site_properties.get(ff_label)
        else:
            type_by_sites = [site.specie.symbol for site in sites]
        # search for site property if not override
        if charges is None:
            charges = sites.site_properties.get("charge")
        if velocities is None:
            velocities = sites.site_properties.get("velocities")
        # validate shape
        if charges is not None:
            charge_arr = np.array(charges)
            assert charge_arr.shape == (len(sites),), "Wrong format for charges"
            charges = charge_arr.tolist()
        if velocities is not None:
            velocities_arr = np.array(velocities)
            assert velocities_arr.shape == (
                len(sites),
                3,
            ), "Wrong format for velocities"
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
    def from_bonding(cls, molecule, bond=True, angle=True, dihedral=True, tol: float = 0.1, **kwargs):
        """
        Another constructor that creates an instance from a molecule.
        Covalent bonds and other bond-based topologies (angles and
        dihedrals) can be automatically determined. Cannot be used for
        non bond-based topologies, e.g., improper dihedrals.

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
        if len(hubs) > 0:
            hub_spokes = {}
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
            for i, j in hub_cons.tolist():
                ks = [k for k in hub_spokes[i] if k != j]
                ls = [l for l in hub_spokes[j] if l != i]
                dihedral_list.extend([[k, i, j, l] for k, l in itertools.product(ks, ls) if k != l])

        topologies = {
            k: v for k, v in zip(SECTION_KEYWORDS["topology"][:3], [bond_list, angle_list, dihedral_list]) if len(v) > 0
        } or None
        return cls(sites=molecule, topologies=topologies, **kwargs)


class ForceField(MSONable):
    """
    Class carrying most data in Masses and force field sections.

    Attributes:
        masses (pandas.DataFrame): DataFrame for Masses section.
        force_field (dict): Force field section keywords (keys) and
            data (values) as DataFrames.
        maps (dict): Dict for labeling atoms and topologies.

    """

    @staticmethod
    def _is_valid(df):
        return not pd.isnull(df).values.any()

    def __init__(self, mass_info, nonbond_coeffs=None, topo_coeffs=None):
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
        for i, m in enumerate(mass_info):
            index.append(i + 1)
            mass = map_mass(m[1])
            masses.append(mass)
            self.mass_info.append((m[0], mass))
            atoms_map[m[0]] = i + 1
        self.masses = pd.DataFrame({"mass": masses}, index=index)
        self.maps = {"Atoms": atoms_map}

        ff_dfs = {}

        self.nonbond_coeffs = nonbond_coeffs
        if self.nonbond_coeffs:
            ff_dfs.update(self._process_nonbond())

        self.topo_coeffs = topo_coeffs
        if self.topo_coeffs:
            self.topo_coeffs = {k: v for k, v in self.topo_coeffs.items() if k in SECTION_KEYWORDS["ff"][2:]}
            for k in self.topo_coeffs:
                coeffs, mapper = self._process_topo(k)
                ff_dfs.update(coeffs)
                self.maps.update(mapper)

        self.force_field = None if len(ff_dfs) == 0 else ff_dfs

    def _process_nonbond(self):
        pair_df = pd.DataFrame(self.nonbond_coeffs)
        assert self._is_valid(pair_df), "Invalid nonbond coefficients with rows varying in length"
        npair, ncoeff = pair_df.shape
        pair_df.columns = [f"coeff{i}" for i in range(1, ncoeff + 1)]
        nm = len(self.mass_info)
        ncomb = int(nm * (nm + 1) / 2)
        if npair == nm:
            kw = "Pair Coeffs"
            pair_df.index = range(1, nm + 1)
        elif npair == ncomb:
            kw = "PairIJ Coeffs"
            ids = list(itertools.combinations_with_replacement(range(1, nm + 1), 2))
            id_df = pd.DataFrame(ids, columns=["id1", "id2"])
            pair_df = pd.concat([id_df, pair_df], axis=1)
        else:
            raise ValueError(f"Expecting {nm} Pair Coeffs or {ncomb} PairIJ Coeffs for {nm} atom types, got {npair}")
        return {kw: pair_df}

    def _process_topo(self, kw):
        def find_eq_types(label, section):
            if section.startswith("Improper"):
                label_arr = np.array(label)
                seqs = [[0, 1, 2, 3], [0, 2, 1, 3], [3, 1, 2, 0], [3, 2, 1, 0]]
                return [tuple(label_arr[s]) for s in seqs]
            return [label] + [label[::-1]]

        main_data, distinct_types = [], []
        class2_data = {k: [] for k in self.topo_coeffs[kw][0] if k in CLASS2_KEYWORDS.get(kw, [])}
        for d in self.topo_coeffs[kw]:
            main_data.append(d["coeffs"])
            distinct_types.append(d["types"])
            for k in class2_data:
                class2_data[k].append(d[k])
        distinct_types = [set(itertools.chain(*(find_eq_types(t, kw) for t in dt))) for dt in distinct_types]
        type_counts = sum(len(dt) for dt in distinct_types)
        type_union = set.union(*distinct_types)
        assert len(type_union) == type_counts, f"Duplicated items found under different coefficients in {kw}"
        atoms = set(np.ravel(list(itertools.chain(*distinct_types))))
        assert atoms.issubset(self.maps["Atoms"]), f"Undefined atom type found in {kw}"
        mapper = {}
        for i, dt in enumerate(distinct_types):
            for t in dt:
                mapper[t] = i + 1

        def process_data(data):
            df = pd.DataFrame(data)
            assert self._is_valid(df), "Invalid coefficients with rows varying in length"
            n, c = df.shape
            df.columns = [f"coeff{i}" for i in range(1, c + 1)]
            df.index = range(1, n + 1)
            return df

        all_data = {kw: process_data(main_data)}
        if class2_data:
            all_data.update({k: process_data(v) for k, v in class2_data.items()})
        return all_data, {kw[:-7] + "s": mapper}

    def to_file(self, filename):
        """
        Saves object to a file in YAML format.

        Args:
            filename (str): Filename.
        """
        d = {
            "mass_info": self.mass_info,
            "nonbond_coeffs": self.nonbond_coeffs,
            "topo_coeffs": self.topo_coeffs,
        }
        with open(filename, "w") as f:
            yaml = YAML()
            yaml.dump(d, f)

    @classmethod
    def from_file(cls, filename):
        """
        Constructor that reads in a file in YAML format.

        Args:
            filename (str): Filename.
        """
        with open(filename) as f:
            yaml = YAML()
            d = yaml.load(f)
        return cls.from_dict(d)

    @classmethod
    def from_dict(cls, d):
        """
        Constructor that reads in a dictionary.

        Args:
            d (dict): Dictionary to read.
        """
        d["mass_info"] = [tuple(m) for m in d["mass_info"]]
        if d.get("topo_coeffs"):
            for v in d["topo_coeffs"].values():
                for c in v:
                    c["types"] = [tuple(t) for t in c["types"]]
        return cls(d["mass_info"], d["nonbond_coeffs"], d["topo_coeffs"])


class CombinedData(LammpsData):
    """
    Object for a collective set of data for a series of LAMMPS data file.
    velocities not yet implemented.
    """

    def __init__(
        self,
        list_of_molecules,
        list_of_names,
        list_of_numbers,
        coordinates,
        atom_style="full",
    ):
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
                [mol.force_field[kw].copy() for mol in self.mols if kw in mol.force_field],
                ignore_index=True,
            )
            self.force_field[kw].index += 1
        if not bool(self.force_field):
            self.force_field = None

        self.atoms = pd.DataFrame()
        mol_count = 0
        type_count = 0
        self.mols_per_data = []
        for i, mol in enumerate(self.mols):
            atoms_df = mol.atoms.copy()
            atoms_df["molecule-ID"] += mol_count
            atoms_df["type"] += type_count
            mols_in_data = len(atoms_df["molecule-ID"].unique())
            self.mols_per_data.append(mols_in_data)
            for _ in range(self.nums[i]):
                self.atoms = pd.concat([self.atoms, atoms_df], ignore_index=True)
                atoms_df["molecule-ID"] += mols_in_data
            type_count += len(mol.masses)
            mol_count += self.nums[i] * mols_in_data
        self.atoms.index += 1
        assert len(self.atoms) == len(self._coordinates), "Wrong number of coordinates."
        self.atoms.update(self._coordinates)

        self.velocities = None
        assert self.mols[0].velocities is None, "Velocities not supported"

        self.topology = {}
        atom_count = 0
        count = {"Bonds": 0, "Angles": 0, "Dihedrals": 0, "Impropers": 0}
        for i, mol in enumerate(self.mols):
            for kw in SECTION_KEYWORDS["topology"]:
                if bool(mol.topology) and kw in mol.topology:
                    if kw not in self.topology:
                        self.topology[kw] = pd.DataFrame()
                    topo_df = mol.topology[kw].copy()
                    topo_df["type"] += count[kw]
                    for col in topo_df.columns[1:]:
                        topo_df[col] += atom_count
                    for _ in range(self.nums[i]):
                        self.topology[kw] = pd.concat([self.topology[kw], topo_df], ignore_index=True)
                        for col in topo_df.columns[1:]:
                            topo_df[col] += len(mol.atoms)
                    count[kw] += len(mol.force_field[kw[:-1] + " Coeffs"])
            atom_count += len(mol.atoms) * self.nums[i]
        for kw in SECTION_KEYWORDS["topology"]:
            if kw in self.topology:
                self.topology[kw].index += 1
        if not bool(self.topology):
            self.topology = None

    @property
    def structure(self):
        """
        Exports a periodic structure object representing the simulation
        box.

        Return:
            Structure
        """
        ld_cp = self.as_lammpsdata()
        return ld_cp.structure

    def disassemble(self, atom_labels=None, guess_element=True, ff_label="ff_map"):
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
        disassembles = []
        for mol in self.mols:
            disassembles.append(
                mol.disassemble(atom_labels=atom_labels, guess_element=guess_element, ff_label=ff_label)
            )
        return disassembles

    @classmethod
    def from_ff_and_topologies(cls):
        """
        Unsupported constructor for CombinedData objects.
        """
        raise AttributeError("Unsupported constructor for CombinedData objects.")

    @classmethod
    def from_structure(cls):
        """
        Unsupported constructor for CombinedData objects.
        """
        raise AttributeError("Unsupported constructor for CombinedData objects.")

    @classmethod
    def parse_xyz(cls, filename):
        """
        Load xyz file generated from packmol (for those who find it hard to install openbabel)

        Returns:
            pandas.DataFrame
        """
        with zopen(filename, "rt") as f:
            lines = f.readlines()

        sio = StringIO("".join(lines[2:]))  # skip the 2nd line
        df = pd.read_csv(
            sio,
            header=None,
            comment="#",
            delim_whitespace=True,
            names=["atom", "x", "y", "z"],
        )
        df.index += 1
        return df

    @classmethod
    def from_files(cls, coordinate_file, list_of_numbers, *filenames):
        """
        Constructor that parse a series of data file.

        Args:
            coordinate_file (str): The filename of xyz coordinates.
            list_of_numbers (list): A list of numbers specifying counts for each
                clusters parsed from files.
            filenames (str): A series of LAMMPS data filenames in string format.
        """
        names = []
        mols = []
        styles = []
        coordinates = cls.parse_xyz(filename=coordinate_file)
        for i in range(0, len(filenames)):
            exec(f"cluster{i + 1} = LammpsData.from_file(filenames[i])")
            names.append(f"cluster{i + 1}")
            mols.append(eval(f"cluster{i + 1}"))
            styles.append(eval(f"cluster{i + 1}").atom_style)
        style = set(styles)
        assert len(style) == 1, "Files have different atom styles."
        return cls.from_lammpsdata(mols, names, list_of_numbers, coordinates, style.pop())

    @classmethod
    def from_lammpsdata(cls, mols, names, list_of_numbers, coordinates, atom_style=None):
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
        styles = []
        for mol in mols:
            styles.append(mol.atom_style)
        style = set(styles)
        assert len(style) == 1, "Data have different atom_style."
        style_return = style.pop()
        if atom_style:
            assert atom_style == style_return, "Data have different atom_style as specified."
        return cls(mols, names, list_of_numbers, coordinates, style_return)

    def get_string(self, distance=6, velocity=8, charge=4, hybrid=True):
        """
        Returns the string representation of CombinedData, essentially
        the string to be written to a file. Combination info is included
        as a comment. For single molecule ID data, the info format is:
            num name
        For data with multiple molecule ID, the format is:
            num(mols_per_data) name

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
            String representation
        """
        lines = LammpsData.get_string(self, distance, velocity, charge, hybrid).splitlines()
        info = "# " + " + ".join(
            (str(a) + " " + b) if c == 1 else (str(a) + "(" + str(c) + ") " + b)
            for a, b, c in zip(self.nums, self.names, self.mols_per_data)
        )
        lines.insert(1, info)
        return "\n".join(lines)

    def as_lammpsdata(self):
        """
        Convert a CombinedData object to a LammpsData object. attributes are deepcopied.

        box (LammpsBox): Simulation box.
        force_field (dict): Data for force field sections. Optional
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
