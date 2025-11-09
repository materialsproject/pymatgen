"""
Module for implementing a CTRL file object class for the Stuttgart
LMTO-ASA code. It will primarily be used to generate a pymatgen
Structure object in the pymatgen.electronic_structure.cohp.py module.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np
from monty.io import zopen

from pymatgen.core.structure import Structure
from pymatgen.core.units import Ry_to_eV, bohr_to_angstrom
from pymatgen.electronic_structure.core import Spin
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.num import round_to_sigfigs

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any

    from typing_extensions import Self

    from pymatgen.core.structure import IStructure

__author__ = "Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Marco Esters"
__email__ = "esters@uoregon.edu"
__date__ = "Nov 30, 2017"


class LMTOCtrl:
    """Parse CTRL files from the Stuttgart LMTO-ASA code.
    Currently, only HEADER, VERS and the structure can be used.
    """

    def __init__(self, structure: Structure | IStructure, header: str | None = None, version: str = "LMASA-47") -> None:
        """
        Args:
            structure (Structure): pymatgen object.
            header (str): The header for the CTRL file. Defaults to None.
            version (str): The LMTO version that is used for the VERS category.
                Defaults to version (4.7).
        """
        self.structure = structure
        self.header = header
        self.version = version

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.get_str() == other.get_str()

    def __repr__(self):
        """Representation of the CTRL file is as a string."""
        return self.get_str()

    def __str__(self):
        """String representation of the CTRL file."""
        return self.get_str()

    def get_str(self, sigfigs=8) -> str:
        """Generate the string representation of the CTRL file. This is
        the minimal CTRL file necessary to execute lmhart.run.
        """
        ctrl_dict = self.as_dict()
        lines = [] if "HEADER" not in ctrl_dict else [f"{'HEADER':<10}{self.header}"]
        if "VERS" in ctrl_dict:
            lines.append("VERS".ljust(10) + self.version)

        lines.append(f"{'STRUC'.ljust(10)}ALAT={ctrl_dict['ALAT']:.{sigfigs}f}")
        for idx, latt in enumerate(ctrl_dict["PLAT"]):
            line = "PLAT=".rjust(15) if idx == 0 else " ".ljust(15)
            line += " ".join(str(round(v, sigfigs)) for v in latt)
            lines.append(line)

        for cat in ("CLASS", "SITE"):
            for a, atoms in enumerate(ctrl_dict[cat]):
                lst = [cat.ljust(9)] if a == 0 else [" ".ljust(9)]
                for token, val in sorted(atoms.items()):
                    if token == "POS":  # noqa: S105
                        lst.append("POS=" + " ".join(str(round(p, sigfigs)) for p in val))
                    else:
                        lst.append(f"{token}={val}")
                line = " ".join(lst)
                lines.append(line)

        return "\n".join(lines) + "\n"

    def as_dict(self):
        """Get the CTRL as a dictionary. "SITE" and "CLASS" are of
        the form {'CATEGORY': {'TOKEN': value}}, the rest is of the
        form 'TOKEN'/'CATEGORY': value. It gets the conventional standard
        structure because primitive cells use the conventional
        a-lattice parameter as the scaling factor and not the a-lattice
        parameter of the primitive cell.
        """
        ctrl_dict = {"@module": type(self).__module__, "@class": type(self).__name__}
        if self.header is not None:
            ctrl_dict["HEADER"] = self.header
        if self.version is not None:
            ctrl_dict["VERS"] = self.version
        sga = SpacegroupAnalyzer(self.structure)
        a_len = sga.get_conventional_standard_structure().lattice.a
        plat = self.structure.lattice.matrix / a_len
        # The following is to find the classes (atoms that are not symmetry equivalent,
        # and create labels. Note that LMTO only attaches numbers with the second atom
        # of the same species, e.g. "Bi", "Bi1", "Bi2", etc.
        eq_atoms = sga.get_symmetry_dataset().equivalent_atoms
        ineq_sites_index = list(set(eq_atoms))
        sites = []
        classes = []
        num_atoms = {}
        for idx, site in enumerate(self.structure):
            atom = site.specie
            label_index = ineq_sites_index.index(eq_atoms[idx])
            if atom.symbol in num_atoms:
                if label_index + 1 > sum(num_atoms.values()):
                    num_atoms[atom.symbol] += 1
                    atom_label = f"{atom.symbol}{num_atoms[atom.symbol] - 1}"
                    classes.append({"ATOM": atom_label, "Z": atom.Z})
            else:
                num_atoms[atom.symbol] = 1
                classes.append({"ATOM": atom.symbol, "Z": atom.Z})
            sites.append({"ATOM": classes[label_index]["ATOM"], "POS": site.coords / a_len})

        update = {
            "ALAT": a_len / bohr_to_angstrom,
            "PLAT": plat,
            "CLASS": classes,
            "SITE": sites,
        }
        return {**ctrl_dict, **update}

    def write_file(self, filename="CTRL", **kwargs):
        """Write a CTRL file with structure, HEADER, and VERS that can be
        used as input for lmhart.run.
        """
        with zopen(filename, mode="wt", encoding="utf-8") as file:
            file.write(self.get_str(**kwargs))

    @classmethod
    def from_file(cls, filename: str | Path = "CTRL", **kwargs) -> Self:
        """
        Creates a CTRL file object from an existing file.

        Args:
            filename: The name of the CTRL file. Defaults to 'CTRL'.

        Returns:
            An LMTOCtrl object.
        """
        with zopen(filename, mode="rt", encoding="utf-8") as file:
            contents = file.read()
        return cls.from_str(contents, **kwargs)  # type:ignore[arg-type]

    @classmethod
    def from_str(cls, data: str, sigfigs: int = 8) -> Self:
        """
        Creates a CTRL file object from a string. This will mostly be
        used to read an LMTOCtrl object from a CTRL file. Empty spheres
        are ignored.

        Args:
            data (str): String representation of the CTRL file.

        Returns:
            An LMTOCtrl object.
        """
        lines = data.split("\n")[:-1]
        _struct_lines: dict[str, list] = {
            "HEADER": [],
            "VERS": [],
            "SYMGRP": [],
            "STRUC": [],  # codespell:ignore struc
            "CLASS": [],
            "SITE": [],
        }

        cat = None
        for line in lines:
            if line != "" and not line.isspace():
                if not line[0].isspace():
                    cat = line.split()[0]
                if cat in _struct_lines:
                    _struct_lines[cat].append(line)
                else:
                    pass

        struct_lines: dict[str, str] = {k: " ".join(v).replace("= ", "=") for k, v in _struct_lines.items()}

        structure_tokens: dict[str, Any] = {
            "ALAT": None,
            "PLAT": [],
            "CLASS": [],
            "SITE": [],
        }

        atom = None
        for cat in ("STRUC", "CLASS", "SITE"):  # codespell:ignore struc
            fields = struct_lines[cat].split("=")
            for idx, field in enumerate(fields):
                token = field.split()[-1]
                if token == "ALAT":  # noqa: S105
                    a_lat = round(float(fields[idx + 1].split()[0]), sigfigs)
                    structure_tokens["ALAT"] = a_lat

                elif token == "ATOM":  # noqa: S105
                    atom = fields[idx + 1].split()[0]
                    if not bool(re.match("E[0-9]*$", atom)):
                        if cat == "CLASS":
                            structure_tokens["CLASS"].append(atom)
                        else:
                            structure_tokens["SITE"].append({"ATOM": atom})
                    else:
                        pass

                elif token in {"PLAT", "POS"}:
                    try:
                        arr = np.array([round(float(i), sigfigs) for i in fields[idx + 1].split()])
                    except ValueError:
                        arr = np.array([round(float(i), sigfigs) for i in fields[idx + 1].split()[:-1]])

                    if token == "PLAT":  # noqa: S105
                        structure_tokens["PLAT"] = arr.reshape([3, 3])
                    elif atom is not None and not bool(re.match("E[0-9]*$", atom)):
                        structure_tokens["SITE"][-1]["POS"] = arr
                    else:
                        pass
                else:
                    pass

        try:
            spc_grp_index = struct_lines["SYMGRP"].index("SPCGRP")
            spc_grp = struct_lines["SYMGRP"][spc_grp_index : spc_grp_index + 12]
            structure_tokens["SPCGRP"] = spc_grp.split("=")[1].split()[0]
        except ValueError:
            pass

        for token in ("HEADER", "VERS"):
            try:
                value = re.split(token + r"\s*", struct_lines[token])[1]
                structure_tokens[token] = value.strip()
            except IndexError:
                pass
        return cls.from_dict(structure_tokens)

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Creates a CTRL file object from a dictionary. The dictionary
        must contain the items "ALAT", PLAT" and "SITE".

        Valid dictionary items are:
            ALAT: the a-lattice parameter
            PLAT: (3x3) array for the lattice vectors
            SITE: list of dictionaries: {'ATOM': class label, 'POS': (3x1) array of fractional coordinates}
            CLASS (optional): list of unique atom labels as str
            SPCGRP (optional): space group symbol (str) or number (int)
            HEADER (optional): HEADER text as a str
            VERS (optional): LMTO version as a str

        Args:
            dct: The CTRL file as a dictionary.

        Returns:
            An LMTOCtrl object.
        """
        dct.setdefault("HEADER", None)
        dct.setdefault("VERS", None)
        alat = dct["ALAT"] * bohr_to_angstrom
        plat = dct["PLAT"] * alat
        species = []
        positions = []
        for site in dct["SITE"]:
            species.append(re.split("[0-9*]", site["ATOM"])[0])
            positions.append(site["POS"] * alat)

        # Only check if the structure is to be generated from the space
        # group if the number of sites is the same as the number of classes.
        # If lattice and the spacegroup don't match, assume it's primitive.
        if "CLASS" in dct and "SPCGRP" in dct and len(dct["SITE"]) == len(dct["CLASS"]):
            try:
                structure = Structure.from_spacegroup(
                    dct["SPCGRP"], plat, species, positions, coords_are_cartesian=True
                )
            except ValueError:
                structure = Structure(
                    plat,
                    species,
                    positions,
                    coords_are_cartesian=True,
                    to_unit_cell=True,
                )
        else:
            structure = Structure(plat, species, positions, coords_are_cartesian=True, to_unit_cell=True)

        return cls(structure, header=dct["HEADER"], version=dct["VERS"])


class LMTOCopl:
    """Read COPL files, which contain COHP data.

    Attributes:
        cohp_data (dict): Contains the COHP data of the form:
            {bond: {"COHP": {Spin.up: cohps, Spin.down:cohps},
                "ICOHP": {Spin.up: icohps, Spin.down: icohps},
                "length": bond length}
        efermi (float): The Fermi energy in Ry or eV.
        energies (list): Sequence of energies in Ry or eV.
        is_spin_polarized (bool): True if the calculation is spin-polarized.
    """

    def __init__(self, filename="COPL", to_eV=False):
        """
        Args:
            filename: filename of the COPL file. Defaults to "COPL".
            to_eV: LMTO-ASA gives energies in Ry. To convert energies into
              eV, set to True. Defaults to False for energies in Ry.
        """
        # COPL files have an extra trailing blank line
        with zopen(filename, mode="rt", encoding="utf-8") as file:
            contents = file.read().split("\n")[:-1]
        # The parameters line is the second line in a COPL file. It
        # contains all parameters that are needed to map the file.
        parameters = contents[1].split()
        num_bonds = int(parameters[0])

        if int(parameters[1]) == 2:
            spins = [Spin.up, Spin.down]
            self.is_spin_polarized = True
        else:
            spins = [Spin.up]
            self.is_spin_polarized = False

        # The COHP data start in row num_bonds + 3
        data = np.array([np.array(row.split(), dtype=float) for row in contents[num_bonds + 2 :]]).transpose()
        if to_eV:
            # LMTO energies have 5 sig figs
            self.energies = np.array(
                [round_to_sigfigs(energy, 5) for energy in data[0] * Ry_to_eV],
                dtype=float,
            )
            self.efermi = round_to_sigfigs(float(parameters[-1]) * Ry_to_eV, 5)
        else:
            self.energies = data[0]
            self.efermi = float(parameters[-1])

        cohp_data = {}
        for bond in range(num_bonds):
            label, length, sites = self._get_bond_data(contents[2 + bond])
            cohp = {spin: data[2 * (bond + s * num_bonds) + 1] for s, spin in enumerate(spins)}
            if to_eV:
                icohp = {
                    spin: np.array([round_to_sigfigs(i, 5) for i in data[2 * (bond + s * num_bonds) + 2] * Ry_to_eV])
                    for s, spin in enumerate(spins)
                }
            else:
                icohp = {spin: data[2 * (bond + s * num_bonds) + 2] for s, spin in enumerate(spins)}

            # This takes care of duplicate labels
            if label in cohp_data:
                idx = 1
                lab = f"{label}-{idx}"
                while lab in cohp_data:
                    idx += 1
                    lab = f"{label}-{idx}"
                label = lab

            cohp_data[label] = {
                "COHP": cohp,
                "ICOHP": icohp,
                "length": length,
                "sites": sites,
            }
        self.cohp_data = cohp_data

    @staticmethod
    def _get_bond_data(line):
        """Subroutine to extract bond label, site indices, and length from
        a COPL header line. The site indices are zero-based, so they
        can be easily used with a Structure object.

        Example header line: Fe-1/Fe-1-tr(-1,-1,-1) : 2.482 Ang.

        Args:
            line: line in the COHPCAR header describing the bond.

        Returns:
            The bond label, the bond length and a tuple of the site indices.
        """
        line = line.split()
        length = float(line[2])
        # Replacing "/" with "-" makes splitting easier
        sites = line[0].replace("/", "-").split("-")
        site_indices = tuple(int(ind) - 1 for ind in sites[1:4:2])
        species = tuple(re.split(r"\d+", spec)[0] for spec in sites[:3:2])
        label = f"{species[0]}{site_indices[0] + 1}-{species[1]}{site_indices[1] + 1}"
        return label, length, site_indices
