"""This module implements input and output processing from PWSCF."""

from __future__ import annotations

import re
from collections import defaultdict
from typing import TYPE_CHECKING

from monty.io import zopen
from monty.re import regrep

from pymatgen.core import Element, Lattice, Structure
from pymatgen.util.io_utils import clean_lines

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any, ClassVar

    from typing_extensions import Self


class PWInput:
    """
    Base input file class. Right now, only supports no symmetry and is
    very basic.
    """

    def __init__(
        self,
        structure,
        pseudo=None,
        control=None,
        system=None,
        electrons=None,
        ions=None,
        cell=None,
        kpoints_mode="automatic",
        kpoints_grid=(1, 1, 1),
        kpoints_shift=(0, 0, 0),
        format_options=None,
    ):
        """Initialize a PWSCF input file.

        Args:
            structure (Structure): Input structure. For spin-polarized calculation,
                properties (e.g. {"starting_magnetization": -0.5,
                "pseudo": "Mn.pbe-sp-van.UPF"}) on each site is needed instead of
                pseudo (dict).
            pseudo (dict): A dict of the pseudopotentials to use. Default to None.
            control (dict): Control parameters. Refer to official PWSCF doc
                on supported parameters. Default to {"calculation": "scf"}
            system (dict): System parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            electrons (dict): Electron parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            ions (dict): Ions parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            cell (dict): Cell parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            kpoints_mode (str): Kpoints generation mode. Default to automatic.
            kpoints_grid (sequence): The kpoint grid. Default to (1, 1, 1).
            kpoints_shift (sequence): The shift for the kpoints. Defaults to
                (0, 0, 0).
            format_options (dict): Formatting options when writing into a string.
                Can be used to specify e.g., the number of decimal places
                (including trailing zeros) for real-space coordinate values
                (atomic positions, cell parameters). Defaults to None,
                in which case the following default values are used
                (so as to maintain backwards compatibility):
                {"indent": 2, "kpoints_crystal_b_indent": 1,
                 "coord_decimals": 6, "atomic_mass_decimals": 4,
                 "kpoints_grid_decimals": 4}.
        """
        self.structure = structure
        sections = {}
        sections["control"] = control or {"calculation": "scf"}
        sections["system"] = system or {}
        sections["electrons"] = electrons or {}
        sections["ions"] = ions or {}
        sections["cell"] = cell or {}
        if pseudo is None:
            for site in structure:
                try:
                    site.properties["pseudo"]
                except KeyError:
                    raise PWInputError(f"Missing {site} in pseudo specification!")
        else:
            for species in self.structure.composition:
                if str(species) not in pseudo:
                    raise PWInputError(f"Missing {species} in pseudo specification!")
        self.pseudo = pseudo

        self.sections = sections
        self.kpoints_mode = kpoints_mode
        self.kpoints_grid = kpoints_grid
        self.kpoints_shift = kpoints_shift
        self.format_options = {
            # Default to 2 spaces for indentation
            "indent": 2,
            # Default to 1 space for indent in kpoint grid entries
            # when kpoints_mode == "crystal_b"
            "kpoints_crystal_b_indent": 1,
            # Default to 6 decimal places
            # for atomic position and cell vector coordinates
            "coord_decimals": 6,
            # Default to 4 decimal places for atomic mass values
            "atomic_mass_decimals": 4,
            # Default to 4 decimal places
            # for kpoint grid entries
            # when kpoints_mode == "crystal_b"
            "kpoints_grid_decimals": 4,
        }
        if format_options is None:
            format_options = {}
        self.format_options.update(format_options)

    def __str__(self):
        out = []
        site_descriptions = {}

        if self.pseudo is not None:
            site_descriptions = self.pseudo
        else:
            c = 1
            for site in self.structure:
                name = None
                for k, v in site_descriptions.items():
                    if site.properties == v:
                        name = k

                if name is None:
                    name = f"{site.specie.symbol}{c}"
                    site_descriptions[name] = site.properties
                    c += 1

        def to_str(v):
            if isinstance(v, str):
                return f"{v!r}"
            if isinstance(v, float):
                return f"{str(v).replace('e', 'd')}"
            if isinstance(v, bool):
                if v:
                    return ".TRUE."
                return ".FALSE."
            return v

        indent = " " * self.format_options["indent"]
        for k1 in ["control", "system", "electrons", "ions", "cell"]:
            v1 = self.sections[k1]
            out.append(f"&{k1.upper()}")
            sub = []
            for k2 in sorted(v1):
                if isinstance(v1[k2], list):
                    n = 1
                    for _ in v1[k2][: len(site_descriptions)]:
                        sub.append(f"{indent}{k2}({n}) = {to_str(v1[k2][n - 1])}")
                        n += 1
                else:
                    sub.append(f"{indent}{k2} = {to_str(v1[k2])}")
            if k1 == "system":
                if "ibrav" not in self.sections[k1]:
                    sub.append(f"{indent}ibrav = 0")
                if "nat" not in self.sections[k1]:
                    sub.append(f"{indent}nat = {len(self.structure)}")
                if "ntyp" not in self.sections[k1]:
                    sub.append(f"{indent}ntyp = {len(site_descriptions)}")
            sub.append("/")
            out.append(",\n".join(sub))

        out.append("ATOMIC_SPECIES")
        prec = self.format_options["atomic_mass_decimals"]
        for k, v in sorted(site_descriptions.items(), key=lambda i: i[0]):
            e = re.match(r"[A-Z][a-z]?", k)[0]
            p = v if self.pseudo is not None else v["pseudo"]
            out.append(f"{indent}{k}  {Element(e).atomic_mass:.{prec}f} {p}")

        out.append("ATOMIC_POSITIONS crystal")
        prec = self.format_options["coord_decimals"]
        if self.pseudo is not None:
            for site in self.structure:
                pos_str = [f"{site.specie}"]
                pos_str.extend([f"{v:.{prec}f}" for v in site.frac_coords])
                out.append(f"{indent}{' '.join(pos_str)}")
        else:
            for site in self.structure:
                name = None
                for k, v in sorted(site_descriptions.items(), key=lambda i: i[0]):
                    if v == site.properties:
                        name = k
                pos_str = [f"{name}"]
                pos_str.extend([f"{v:.{prec}f}" for v in site.frac_coords])
                out.append(f"{indent}{' '.join(pos_str)}")

        out.append(f"K_POINTS {self.kpoints_mode}")
        if self.kpoints_mode == "automatic":
            kpt_str = [f"{i}" for i in self.kpoints_grid]
            kpt_str.extend([f"{i}" for i in self.kpoints_shift])
            out.append(f"{indent}{' '.join(kpt_str)}")
        elif self.kpoints_mode == "crystal_b":
            kpt_indent = " " * self.format_options["kpoints_crystal_b_indent"]
            out.append(f"{kpt_indent}{len(self.kpoints_grid)}")
            prec = self.format_options["kpoints_grid_decimals"]
            for i in range(len(self.kpoints_grid)):
                kpt_str = [f"{entry:.{prec}f}" for entry in self.kpoints_grid[i]]
                out.append(f"{kpt_indent}{' '.join(kpt_str)}")
        elif self.kpoints_mode == "gamma":
            pass

        out.append("CELL_PARAMETERS angstrom")
        prec = self.format_options["coord_decimals"]
        for vec in self.structure.lattice.matrix:
            vec_str = [f"{v:.{prec}f}" for v in vec]
            out.append(f"{indent}{' '.join(vec_str)}")
        return "\n".join(out)

    def as_dict(self):
        """
        Create a dictionary representation of a PWInput object.

        Returns:
            dict
        """
        return {
            "structure": self.structure.as_dict(),
            "pseudo": self.pseudo,
            "sections": self.sections,
            "kpoints_mode": self.kpoints_mode,
            "kpoints_grid": self.kpoints_grid,
            "kpoints_shift": self.kpoints_shift,
            "format_options": self.format_options,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Load a PWInput object from a dictionary.

        Args:
            dct (dict): dictionary with PWInput data

        Returns:
            PWInput object
        """
        return cls(
            structure=Structure.from_dict(dct["structure"]),
            pseudo=dct["pseudo"],
            control=dct["sections"]["control"],
            system=dct["sections"]["system"],
            electrons=dct["sections"]["electrons"],
            ions=dct["sections"]["ions"],
            cell=dct["sections"]["cell"],
            kpoints_mode=dct["kpoints_mode"],
            kpoints_grid=dct["kpoints_grid"],
            kpoints_shift=dct["kpoints_shift"],
            format_options=dct["format_options"],
        )

    def write_file(self, filename):
        """Write the PWSCF input file.

        Args:
            filename (str): The string filename to output to.
        """
        with open(filename, mode="w", encoding="utf-8") as file:
            file.write(str(self))

    @classmethod
    def from_file(cls, filename: str | Path) -> Self:
        """
        Reads an PWInput object from a file.

        Args:
            filename (str | Path): Filename for file

        Returns:
            PWInput object
        """
        with zopen(filename, mode="rt") as file:
            return cls.from_str(file.read())

    @classmethod
    def from_str(cls, string: str) -> Self:
        """
        Reads an PWInput object from a string.

        Args:
            string (str): PWInput string

        Returns:
            PWInput object
        """
        lines = list(clean_lines(string.splitlines(), rstrip_only=True))

        def input_mode(line):
            if line[0] == "&":
                return ("sections", line[1:].lower())
            if "ATOMIC_SPECIES" in line:
                return ("pseudo",)
            if "K_POINTS" in line:
                return "kpoints", line.split()[1]
            if "OCCUPATIONS" in line:
                return "occupations"
            if "CELL_PARAMETERS" in line or "ATOMIC_POSITIONS" in line:
                return "structure", line.split()[1]
            if line == "/":
                return None
            return mode

        sections: dict[str, dict] = {
            "control": {},
            "system": {},
            "electrons": {},
            "ions": {},
            "cell": {},
        }
        pseudo = {}
        lattice = []
        species = []
        coords = []
        structure = None
        site_properties: dict[str, list] = {"pseudo": []}
        mode = None
        kpoints_mode = None
        kpoints_grid = (1, 1, 1)
        kpoints_shift = (0, 0, 0)
        coords_are_cartesian = False
        format_options = {}

        for line in lines:
            mode = input_mode(line)
            if mode is None:
                pass
            elif mode[0] == "sections":
                section = mode[1]
                if match := re.match(r"^(\s*)(\w+)\(?(\d*?)\)?\s*=\s*(.*)", line):
                    format_options["indent"] = len(match[1])
                    key = match[2].strip()
                    key_ = match[3].strip()
                    val = match[4].strip().rstrip(",")
                    if key_ != "":
                        if sections[section].get(key) is None:
                            val_ = [0.0] * 20  # MAX NTYP DEFINITION
                            val_[int(key_) - 1] = PWInput.proc_val(key, val)
                            sections[section][key] = val_

                            site_properties[key] = []
                        else:
                            sections[section][key][int(key_) - 1] = PWInput.proc_val(key, val)
                    else:
                        sections[section][key] = PWInput.proc_val(key, val)

            elif mode[0] == "pseudo":
                if match := re.match(r"^(\s*)(\w+\d*[\+-]?)\s+(\d*.\d*)\s+(.*)", line):
                    format_options["indent"] = len(match[1])
                    pseudo[match[2].strip()] = match[4].strip()

            elif mode[0] == "kpoints":
                if match := re.match(r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line):
                    kpoints_grid = (int(match[1]), int(match[2]), int(match[3]))
                    kpoints_shift = (int(match[4]), int(match[5]), int(match[6]))
                else:
                    kpoints_mode = mode[1]

            elif mode[0] == "structure":
                m_l = re.match(r"^(\s*)(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)", line)
                m_p = re.match(r"^(\s*)(\w+\d*[\+-]?)\s+(-?\d+\.\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)", line)
                if m_l:
                    format_options["indent"] = len(m_l[1])
                    lattice += [
                        float(m_l[2]),
                        float(m_l[3]),
                        float(m_l[4]),
                    ]
                    decimals = max(
                        # length of decimal digits; 0 if no decimal digits
                        (len(dec[1]) if len(dec := v.split(".")) == 2 else 0)
                        for v in (m_l[2], m_l[3], m_l[4])
                    )
                    format_options["coord_decimals"] = max(
                        format_options.get("coord_decimals", 0),
                        decimals,
                    )

                elif m_p:
                    format_options["indent"] = len(m_p[1])
                    site_properties["pseudo"].append(pseudo[m_p[2]])
                    species.append(m_p[2])
                    coords += [[float(m_p[3]), float(m_p[4]), float(m_p[5])]]
                    decimals = max(
                        # length of decimal digits; 0 if no decimal digits
                        (len(dec[1]) if len(dec := v.split(".")) == 2 else 0)
                        for v in (m_p[3], m_p[4], m_p[5])
                    )
                    format_options["coord_decimals"] = max(
                        format_options.get("coord_decimals", 0),
                        decimals,
                    )

                    if mode[1] == "angstrom":
                        coords_are_cartesian = True

        structure = Structure(
            Lattice(lattice),
            species,
            coords,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties,
        )
        return cls(
            structure=structure,
            control=sections["control"],
            pseudo=pseudo,
            system=sections["system"],
            electrons=sections["electrons"],
            ions=sections["ions"],
            cell=sections["cell"],
            kpoints_mode=kpoints_mode,
            kpoints_grid=kpoints_grid,
            kpoints_shift=kpoints_shift,
            format_options=format_options,
        )

    @staticmethod
    def proc_val(key, val):
        """Static helper method to convert PWINPUT parameters to proper type, e.g.
        integers, floats, etc.

        Args:
            key: PWINPUT parameter key
            val: Actual value of PWINPUT parameter.
        """
        float_keys = (
            "etot_conv_thr",
            "forc_conv_thr",
            "conv_thr",
            "Hubbard_U",
            "Hubbard_J0",
            "degauss",
            "starting_magnetization",
        )

        int_keys = (
            "nstep",
            "iprint",
            "nberrycyc",
            "gdir",
            "nppstr",
            "ibrav",
            "nat",
            "ntyp",
            "nbnd",
            "nr1",
            "nr2",
            "nr3",
            "nr1s",
            "nr2s",
            "nr3s",
            "nspin",
            "nqx1",
            "nqx2",
            "nqx3",
            "lda_plus_u_kind",
            "edir",
            "report",
            "esm_nfit",
            "space_group",
            "origin_choice",
            "electron_maxstep",
            "mixing_ndim",
            "mixing_fixed_ns",
            "ortho_para",
            "diago_cg_maxiter",
            "diago_david_ndim",
            "nraise",
            "bfgs_ndim",
            "if_pos",
            "nks",
            "nk1",
            "nk2",
            "nk3",
            "sk1",
            "sk2",
            "sk3",
            "nconstr",
        )

        bool_keys = (
            "wf_collect",
            "tstress",
            "tprnfor",
            "lkpoint_dir",
            "tefield",
            "dipfield",
            "lelfield",
            "lorbm",
            "lberry",
            "lfcpopt",
            "monopole",
            "nosym",
            "nosym_evc",
            "noinv",
            "no_t_rev",
            "force_symmorphic",
            "use_all_frac",
            "one_atom_occupations",
            "starting_spin_angle",
            "noncolin",
            "x_gamma_extrapolation",
            "lda_plus_u",
            "lspinorb",
            "london",
            "ts_vdw_isolated",
            "xdm",
            "uniqueb",
            "rhombohedral",
            "realxz",
            "block",
            "scf_must_converge",
            "adaptive_thr",
            "diago_full_acc",
            "tqr",
            "remove_rigid_rot",
            "refold_pos",
        )

        def smart_int_or_float(num_str):
            if num_str.find(".") != -1 or num_str.lower().find("e") != -1:
                return float(num_str)
            return int(num_str)

        try:
            if key in bool_keys:
                if val.lower() == ".true.":
                    return True
                if val.lower() == ".false.":
                    return False
                raise ValueError(f"{key} should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*d?-?\d*", val.lower())[0].replace("d", "e"))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val)[0])

        except ValueError:
            pass

        try:
            return smart_int_or_float(val.replace("d", "e"))
        except ValueError:
            pass

        if "true" in val.lower():
            return True
        if "false" in val.lower():
            return False

        if match := re.match(r"^[\"|'](.+)[\"|']$", val):
            return match[1]
        return None


class PWInputError(BaseException):
    """Error for PWInput."""


class PWOutput:
    """Parser for PWSCF output file."""

    patterns: ClassVar[dict[str, str]] = {
        "energies": r"total energy\s+=\s+([\d\.\-]+)\sRy",
        "ecut": r"kinetic\-energy cutoff\s+=\s+([\d\.\-]+)\s+Ry",
        "lattice_type": r"bravais\-lattice index\s+=\s+(\d+)",
        "celldm1": r"celldm\(1\)=\s+([\d\.]+)\s",
        "celldm2": r"celldm\(2\)=\s+([\d\.]+)\s",
        "celldm3": r"celldm\(3\)=\s+([\d\.]+)\s",
        "celldm4": r"celldm\(4\)=\s+([\d\.]+)\s",
        "celldm5": r"celldm\(5\)=\s+([\d\.]+)\s",
        "celldm6": r"celldm\(6\)=\s+([\d\.]+)\s",
        "nkpts": r"number of k points=\s+([\d]+)",
    }

    def __init__(self, filename):
        """
        Args:
            filename (str): Filename.
        """
        self.filename = filename
        self.data: dict[str, Any] = defaultdict(list)
        self.read_pattern(PWOutput.patterns)
        for k, v in self.data.items():
            if k == "energies":
                self.data[k] = [float(i[0][0]) for i in v]
            elif k in ["lattice_type", "nkpts"]:
                self.data[k] = int(v[0][0][0])
            else:
                self.data[k] = float(v[0][0][0])

    def read_pattern(self, patterns, reverse=False, terminate_on_match=False, postprocess=str):
        r"""
        General pattern reading. Uses monty's regrep method. Takes the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.
                {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"}.
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.

        Renders accessible:
            Any attribute in patterns. For example,
            {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"} will set the
            value of self.data["energy"] = [[-1234], [-3453], ...], to the
            results from regex and postprocess. Note that the returned
            values are lists of lists, because you can grep multiple
            items on one line.
        """
        matches = regrep(
            self.filename,
            patterns,
            reverse=reverse,
            terminate_on_match=terminate_on_match,
            postprocess=postprocess,
        )
        self.data.update(matches)

    def get_celldm(self, idx: int):
        """
        Args:
            idx (int): index.

        Returns:
            Cell dimension along index
        """
        return self.data[f"celldm{idx}"]

    @property
    def final_energy(self) -> float:
        """The final energy from the PW output."""
        return self.data["energies"][-1]

    @property
    def lattice_type(self) -> int:
        """The lattice type."""
        return self.data["lattice_type"]
