# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements input and output processing from PWSCF.
"""

from __future__ import annotations

import re
from collections import defaultdict

from monty.io import zopen
from monty.re import regrep

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.util.io_utils import clean_lines


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
    ):
        """
        Initializes a PWSCF input file.

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
                    name = site.specie.symbol + str(c)
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

        for k1 in ["control", "system", "electrons", "ions", "cell"]:
            v1 = self.sections[k1]
            out.append(f"&{k1.upper()}")
            sub = []
            for k2 in sorted(v1):
                if isinstance(v1[k2], list):
                    n = 1
                    for _ in v1[k2][: len(site_descriptions)]:
                        sub.append(f"  {k2}({n}) = {to_str(v1[k2][n - 1])}")
                        n += 1
                else:
                    sub.append(f"  {k2} = {to_str(v1[k2])}")
            if k1 == "system":
                if "ibrav" not in self.sections[k1]:
                    sub.append("  ibrav = 0")
                if "nat" not in self.sections[k1]:
                    sub.append(f"  nat = {len(self.structure)}")
                if "ntyp" not in self.sections[k1]:
                    sub.append(f"  ntyp = {len(site_descriptions)}")
            sub.append("/")
            out.append(",\n".join(sub))

        out.append("ATOMIC_SPECIES")
        for k, v in sorted(site_descriptions.items(), key=lambda i: i[0]):
            e = re.match(r"[A-Z][a-z]?", k).group(0)
            if self.pseudo is not None:
                p = v
            else:
                p = v["pseudo"]
            out.append(f"  {k}  {Element(e).atomic_mass:.4f} {p}")

        out.append("ATOMIC_POSITIONS crystal")
        if self.pseudo is not None:
            for site in self.structure:
                out.append(f"  {site.specie} {site.a:.6f} {site.b:.6f} {site.c:.6f}")
        else:
            for site in self.structure:
                name = None
                for k, v in sorted(site_descriptions.items(), key=lambda i: i[0]):
                    if v == site.properties:
                        name = k
                out.append(f"  {name} {site.a:.6f} {site.b:.6f} {site.c:.6f}")

        out.append(f"K_POINTS {self.kpoints_mode}")
        if self.kpoints_mode == "automatic":
            kpt_str = [f"{i}" for i in self.kpoints_grid]
            kpt_str.extend([f"{i}" for i in self.kpoints_shift])
            out.append(f"  {' '.join(kpt_str)}")
        elif self.kpoints_mode == "crystal_b":
            out.append(f" {str(len(self.kpoints_grid))}")
            for i in range(len(self.kpoints_grid)):
                kpt_str = [f"{entry:.4f}" for entry in self.kpoints_grid[i]]
                out.append(f" {' '.join(kpt_str)}")
        elif self.kpoints_mode == "gamma":
            pass

        out.append("CELL_PARAMETERS angstrom")
        for vec in self.structure.lattice.matrix:
            out.append(f"  {vec[0]:f} {vec[1]:f} {vec[2]:f}")
        return "\n".join(out)

    def as_dict(self):
        """
        Create a dictionary representation of a PWInput object

        Returns:
            dict
        """
        pwinput_dict = {
            "structure": self.structure.as_dict(),
            "pseudo": self.pseudo,
            "sections": self.sections,
            "kpoints_mode": self.kpoints_mode,
            "kpoints_grid": self.kpoints_grid,
            "kpoints_shift": self.kpoints_shift,
        }
        return pwinput_dict

    @classmethod
    def from_dict(cls, pwinput_dict):
        """
        Load a PWInput object from a dictionary.

        Args:
            pwinput_dict (dict): dictionary with PWInput data

        Returns:
            PWInput object
        """
        pwinput = cls(
            structure=Structure.from_dict(pwinput_dict["structure"]),
            pseudo=pwinput_dict["pseudo"],
            control=pwinput_dict["sections"]["control"],
            system=pwinput_dict["sections"]["system"],
            electrons=pwinput_dict["sections"]["electrons"],
            ions=pwinput_dict["sections"]["ions"],
            cell=pwinput_dict["sections"]["cell"],
            kpoints_mode=pwinput_dict["kpoints_mode"],
            kpoints_grid=pwinput_dict["kpoints_grid"],
            kpoints_shift=pwinput_dict["kpoints_shift"],
        )
        return pwinput

    def write_file(self, filename):
        """
        Write the PWSCF input file.

        Args:
            filename (str): The string filename to output to.
        """
        with open(filename, "w") as f:
            f.write(str(self))

    @staticmethod
    def from_file(filename):
        """
        Reads an PWInput object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            PWInput object
        """
        with zopen(filename, "rt") as f:
            return PWInput.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads an PWInput object from a string.

        Args:
            string (str): PWInput string

        Returns:
            PWInput object
        """
        lines = list(clean_lines(string.splitlines()))

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

        sections = {
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
        site_properties = {"pseudo": []}
        mode = None
        for line in lines:
            mode = input_mode(line)
            if mode is None:
                pass
            elif mode[0] == "sections":
                section = mode[1]
                m = re.match(r"(\w+)\(?(\d*?)\)?\s*=\s*(.*)", line)
                if m:
                    key = m.group(1).strip()
                    key_ = m.group(2).strip()
                    val = m.group(3).strip()
                    if key_ != "":
                        if sections[section].get(key, None) is None:
                            val_ = [0.0] * 20  # MAX NTYP DEFINITION
                            val_[int(key_) - 1] = PWInput.proc_val(key, val)
                            sections[section][key] = val_

                            site_properties[key] = []
                        else:
                            sections[section][key][int(key_) - 1] = PWInput.proc_val(key, val)
                    else:
                        sections[section][key] = PWInput.proc_val(key, val)

            elif mode[0] == "pseudo":
                m = re.match(r"(\w+)\s+(\d*.\d*)\s+(.*)", line)
                if m:
                    pseudo[m.group(1).strip()] = m.group(3).strip()
            elif mode[0] == "kpoints":
                m = re.match(r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
                if m:
                    kpoints_grid = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
                    kpoints_shift = (int(m.group(4)), int(m.group(5)), int(m.group(6)))
                else:
                    kpoints_mode = mode[1]
                    kpoints_grid = (1, 1, 1)
                    kpoints_shift = (0, 0, 0)

            elif mode[0] == "structure":
                m_l = re.match(r"(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)", line)
                m_p = re.match(r"(\w+)\s+(-?\d+\.\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)", line)
                if m_l:
                    lattice += [
                        float(m_l.group(1)),
                        float(m_l.group(2)),
                        float(m_l.group(3)),
                    ]
                elif m_p:
                    site_properties["pseudo"].append(pseudo[m_p.group(1)])
                    species.append(m_p.group(1))
                    coords += [[float(m_p.group(2)), float(m_p.group(3)), float(m_p.group(4))]]

                    if mode[1] == "angstrom":
                        coords_are_cartesian = True
                    elif mode[1] == "crystal":
                        coords_are_cartesian = False
        structure = Structure(
            Lattice(lattice),
            species,
            coords,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties,
        )
        return PWInput(
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
        )

    @staticmethod
    def proc_val(key, val):
        """
        Static helper method to convert PWINPUT parameters to proper type, e.g.,
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
            "defauss",
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

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            return int(numstr)

        try:
            if key in bool_keys:
                if val.lower() == ".true.":
                    return True
                if val.lower() == ".false.":
                    return False
                raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*d?-?\d*", val.lower()).group(0).replace("d", "e"))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))

        except ValueError:
            pass

        try:
            val = val.replace("d", "e")
            return smart_int_or_float(val)
        except ValueError:
            pass

        if "true" in val.lower():
            return True
        if "false" in val.lower():
            return False

        m = re.match(r"^[\"|'](.+)[\"|']$", val)
        if m:
            return m.group(1)


class PWInputError(BaseException):
    """
    Error for PWInput
    """


class PWOutput:
    """
    Parser for PWSCF output file.
    """

    patterns = {
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
            filename (str): Filename
        """
        self.filename = filename
        self.data = defaultdict(list)
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
            patterns (dict): A dict of patterns, e.g.,
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
            idx (int): index

        Returns:
            Cell dimension along index
        """
        return self.data[f"celldm{idx}"]

    @property
    def final_energy(self):
        """
        Returns: Final energy
        """
        return self.data["energies"][-1]

    @property
    def lattice_type(self):
        """
        Returns: Lattice type.
        """
        return self.data["lattice_type"]
