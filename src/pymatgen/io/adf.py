"""IO for ADF files."""

from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING

from monty.io import reverse_readfile
from monty.itertools import chunks
from monty.json import MSONable

from pymatgen.core.structure import Molecule

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import ClassVar

    from typing_extensions import Self

__author__ = "Xin Chen, chenxin13@mails.tsinghua.edu.cn"


class AdfInputError(Exception):
    """The default error class for ADF."""


class AdfOutputError(Exception):
    """The default error class for errors raised by ``AdfOutput``."""


class AdfKey(MSONable):
    """
    The basic input unit for ADF. A key is a string of characters that does not
    contain a delimiter (blank, comma or equal sign). A key may have multiple
    subkeys and a set of options.
    """

    block_keys = (
        "SCF",
        "GEOMETRY",
        "XC",
        "UNITS",
        "ATOMS",
        "CHARGE",
        "BASIS",
        "SYMMETRY",
        "RELATIVISTIC",
        "OCCUPATIONS",
        "SAVE",
        "A1FIT",
        "INTEGRATION",
        "UNRESTRICTED",
        "ZLMFIT",
        "TITLE",
        "EXACTDENSITY",
        "TOTALENERGY",
        "ANALYTICALFREQ",
    )
    sub_keys = ("AtomDepQuality",)

    # Full blocks are blocks that must have an 'END'.
    _full_blocks = ("GEOMETRY", "SCF", "UNITS", "BASIS", "ANALYTICALFREQ")

    def __init__(self, name, options=None, subkeys=None):
        """
        Initialization method.

        Args:
            name (str): The name of this key.
            options : Sized
                The options for this key. Each element can be a primitive object or
                a tuple/list with two elements: the first is the name and the second is a primitive object.
            subkeys (Sized): The subkeys for this key.

        Raises:
            ValueError: If elements in ``subkeys`` are not ``AdfKey`` objects.
        """
        self.name = name
        self.options = options if options is not None else []
        self.subkeys = subkeys if subkeys is not None else []
        if len(self.subkeys) > 0:
            for k in subkeys:
                if not isinstance(k, AdfKey):
                    raise TypeError("Not all subkeys are ``AdfKey`` objects!")
        self._sized_op = None
        if len(self.options) > 0:
            self._sized_op = isinstance(self.options[0], list | tuple)

    def _options_string(self):
        """Return the option string."""
        if len(self.options) > 0:
            opt_str = ""
            for op in self.options:
                if self._sized_op:
                    opt_str += f"{op[0]}={op[1]} "
                else:
                    opt_str += f"{op} "
            return opt_str.strip()
        return ""

    def is_block_key(self) -> bool:
        """Return True if this key is a block key."""
        return self.name.upper() in self.block_keys

    @property
    def key(self) -> str:
        """The name of this key. If this is a block key, the name will be converted to upper cases."""
        if self.is_block_key():
            return self.name.upper()
        return self.name

    def __str__(self):
        """Get the string representation of this ``AdfKey``.

        Notes:
            If this key is 'Atoms' and the coordinates are in Cartesian form,
            a different string format will be used.
        """
        adf_str = f"{self.key}"
        if len(self.options) > 0:
            adf_str += f" {self._options_string()}"
        adf_str += "\n"
        if len(self.subkeys) > 0:
            if self.key.lower() == "atoms":
                for subkey in self.subkeys:
                    adf_str += (
                        f"{subkey.name:2s}  {subkey.options[0]: 14.8f}"
                        f"    {subkey.options[1]: 14.8f}    {subkey.options[2]: 14.8f}\n"
                    )
            else:
                for subkey in self.subkeys:
                    adf_str += str(subkey)
            if self.is_block_key():
                adf_str += "END\n"
            else:
                adf_str += "subend\n"
        elif self.key.upper() in self._full_blocks:
            adf_str += "END\n"
        return adf_str

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AdfKey):
            return False
        return str(self) == str(other)

    def has_subkey(self, subkey: str | AdfKey) -> bool:
        """
        Args:
            subkey (str | AdfKey): A key name or AdfKey object.

        Returns:
            bool: True if this key contains the given subkey.
        """
        if isinstance(subkey, str):
            key = subkey
        elif isinstance(subkey, AdfKey):
            key = subkey.key
        else:
            raise TypeError("The subkey should be an AdfKey or a string!")
        return len(self.subkeys) > 0 and key in (k.key for k in self.subkeys)

    def add_subkey(self, subkey):
        """
        Add a new subkey to this key.

        Args:
            subkey (AdfKey): A new subkey.

        Notes:
            Duplicate check will not be performed if this is an 'Atoms' block.
        """
        if self.key.lower() == "atoms" or not self.has_subkey(subkey):
            self.subkeys.append(subkey)

    def remove_subkey(self, subkey):
        """
        Remove the given subkey, if existed, from this AdfKey.

        Args:
            subkey (str or AdfKey): The subkey to remove.
        """
        if len(self.subkeys) > 0:
            key = subkey if isinstance(subkey, str) else subkey.key
            for idx, sk in enumerate(self.subkeys):
                if sk.key == key:
                    self.subkeys.pop(idx)
                    break

    def add_option(self, option):
        """
        Add a new option to this key.

        Args:
            option : Sized or str or int or float
                A new option to add. This must have the same format
                with existing options.

        Raises:
            TypeError: If the format of the given ``option`` is different.
        """
        if len(self.options) != 0:
            sized_op = isinstance(option, list | tuple)
            if self._sized_op != sized_op:
                raise TypeError("Option type is mismatched!")

        self.options.append(option)

    def remove_option(self, option: str | int) -> None:
        """
        Remove an option.

        Args:
            option (str | int):  The name or index of the option to remove.

        Raises:
            TypeError: If the option has a wrong type.
        """
        if len(self.options) > 0:
            if self._sized_op:
                if not isinstance(option, str):
                    raise TypeError("``option`` should be a name string!")
                for idx, val in enumerate(self.options):
                    if val[0] == option:
                        self.options.pop(idx)
                        break
            else:
                if not isinstance(option, int):
                    raise TypeError("``option`` should be an integer index!")
                self.options.pop(option)

    def has_option(self, option: str) -> bool:
        """
        Args:
            option (str): The option.

        Returns:
            bool: True if this AdfKey has the given option.
        """
        if len(self.options) == 0:
            return False
        return any((self._sized_op and op[0] == option) or op == option for op in self.options)

    def as_dict(self):
        """A JSON-serializable dict representation of self."""
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "name": self.name,
            "options": self.options,
        }
        if len(self.subkeys) > 0:
            subkeys = []
            for subkey in self.subkeys:
                subkeys.append(subkey.as_dict())
            dct["subkeys"] = subkeys
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Construct a MSONable AdfKey object from the JSON dict.

        Args:
            dct (dict): A dict of saved attributes.

        Returns:
            AdfKey: An AdfKey object recovered from the JSON dict.
        """
        key = dct.get("name")
        options = dct.get("options")
        subkey_list = dct.get("subkeys", [])
        subkeys = [AdfKey.from_dict(k) for k in subkey_list] or None
        return cls(key, options, subkeys)

    @classmethod
    def from_str(cls, string: str) -> Self:
        """
        Construct an AdfKey object from the string.

        Args:
            string: str

        Returns:
            An AdfKey object recovered from the string.

        Raises:
            ValueError: Currently nested subkeys are not supported.
                If ``subend`` was found a ValueError would be raised.

        Notes:
            Only the first block key will be returned.
        """

        def is_float(string: str) -> bool:
            return "." in string or "E" in string or "e" in string

        def is_numeric(string: str) -> bool:
            """True if input string is numeric and can be converted to an int or a float."""
            try:
                float(string)
                return True
            except ValueError:
                return False

        if "\n" not in string:
            el: list[str] = string.split()
            if len(el) > 1:
                options: list[list[str]] | list[str] = [s.split("=") for s in el[1:]] if "=" in string else el[1:]
                for idx, op in enumerate(options):
                    if isinstance(op, list):
                        if is_numeric(op[1]):
                            op[1] = float(op[1]) if is_float(op[1]) else int(op[1])
                    elif isinstance(op, str) and is_numeric(op):
                        options[idx] = float(op) if is_float(op) else int(op)  # type: ignore[call-overload]
                return cls(el[0], options)

            return cls(el[0], None)

        if "subend" in string:
            raise ValueError("Nested subkeys are not supported!")

        def iter_lines(s: str) -> Generator[str, None, None]:
            r"""A generator form of s.split('\n') for reducing memory overhead.

            Args:
                s (str): A multi-line string.

            Yields:
                str: line
            """
            prev_nl = -1
            while True:
                next_nl = s.find("\n", prev_nl + 1)
                if next_nl < 0:
                    yield s[(prev_nl + 1) :]
                    break
                yield s[(prev_nl + 1) : next_nl]
                prev_nl = next_nl

        key = None
        for line in iter_lines(string):
            if line == "":
                continue
            el = line.strip().split()
            if len(el) == 0:
                continue
            if el[0].upper() in cls.block_keys:
                if key is None:
                    key = cls.from_str(line)
                else:
                    return key
            elif el[0].upper() == "END":
                return key  # type: ignore[return-value]
            elif key is not None:
                key.add_subkey(cls.from_str(line))

        raise KeyError("Incomplete key: 'END' is missing!")


class AdfTask(MSONable):
    """
    Basic task for ADF. All settings in this class are independent of molecules.

    Notes:
        Unlike other quantum chemistry packages (NWChem, Gaussian, ...),
        ADF does not support calculating force/gradient.
    """

    operations: ClassVar[dict[str, str]] = {
        "energy": "Evaluate the single point energy.",
        "optimize": "Minimize the energy by varying the molecular structure.",
        "frequencies": "Compute second derivatives and print out an analysis of molecular vibrations.",
        "freq": "Same as frequencies.",
        "numerical_frequencies": "Compute molecular frequencies using numerical method.",
    }

    def __init__(
        self,
        operation="energy",
        basis_set=None,
        xc=None,
        title="ADF_RUN",
        units=None,
        geo_subkeys=None,
        scf=None,
        other_directives=None,
    ):
        """
        Initialization method.

        Args:
            operation (str): The target operation.
            basis_set (AdfKey): The basis set definitions for this task. Defaults to 'DZ/Large'.
            xc (AdfKey): The exchange-correlation functionals. Defaults to PBE.
            title (str): The title of this ADF task.
            units (AdfKey): The units. Defaults to Angstroms/Degree.
            geo_subkeys (Sized): The subkeys for the block key 'GEOMETRY'.
            scf (AdfKey): The scf options.
            other_directives (Sized): User-defined directives.
        """
        if operation not in self.operations:
            raise AdfInputError(f"Invalid ADF task {operation}")
        self.operation = operation
        self.title = title
        self.basis_set = basis_set if basis_set is not None else self.get_default_basis_set()
        self.xc = xc if xc is not None else self.get_default_xc()
        self.units = units if units is not None else self.get_default_units()
        self.scf = scf if scf is not None else self.get_default_scf()
        self.other_directives = other_directives if other_directives is not None else []
        self._setup_task(geo_subkeys)

    @staticmethod
    def get_default_basis_set():
        """Get Default basis set."""
        return AdfKey.from_str("Basis\ntype DZ\ncore small\nEND")

    @staticmethod
    def get_default_scf():
        """Get ADF using default SCF."""
        return AdfKey.from_str("SCF\niterations 300\nEND")

    @staticmethod
    def get_default_geo():
        """Get ADFKey using default geometry."""
        return AdfKey.from_str("GEOMETRY SinglePoint\nEND")

    @staticmethod
    def get_default_xc():
        """Get ADFKey using default XC."""
        return AdfKey.from_str("XC\nGGA PBE\nEND")

    @staticmethod
    def get_default_units():
        """Get Default units."""
        return AdfKey.from_str("Units\nlength angstrom\nangle degree\nEnd")

    def _setup_task(self, geo_subkeys):
        """Setup the block 'Geometry' given subkeys and the task.

        Args:
            geo_subkeys (Sized): User-defined subkeys for the block 'Geometry'.

        Notes:
            Most of the run types of ADF are specified in the Geometry
            block except the 'AnalyticFreq'.
        """
        self.geo = AdfKey("Geometry", subkeys=geo_subkeys)
        if self.operation.lower() == "energy":
            self.geo.add_option("SinglePoint")
            if self.geo.has_subkey("Frequencies"):
                self.geo.remove_subkey("Frequencies")
        elif self.operation.lower() == "optimize":
            self.geo.add_option("GeometryOptimization")
            if self.geo.has_subkey("Frequencies"):
                self.geo.remove_subkey("Frequencies")
        elif self.operation.lower() == "numerical_frequencies":
            self.geo.add_subkey(AdfKey("Frequencies"))
        else:
            self.other_directives.append(AdfKey("AnalyticalFreq"))
            if self.geo.has_subkey("Frequencies"):
                self.geo.remove_subkey("Frequencies")

    def __str__(self):
        out = f"""TITLE {self.title}\n
{self.units}
{self.xc}
{self.basis_set}
{self.scf}
{self.geo}\n"""
        for block_key in self.other_directives:
            if not isinstance(block_key, AdfKey):
                raise TypeError(f"{block_key} is not an AdfKey!")
            out += str(block_key) + "\n"
        return out

    def as_dict(self):
        """A JSON-serializable dict representation of self."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "operation": self.operation,
            "title": self.title,
            "xc": self.xc.as_dict(),
            "basis_set": self.basis_set.as_dict(),
            "units": self.units.as_dict(),
            "scf": self.scf.as_dict(),
            "geo": self.geo.as_dict(),
            "others": [k.as_dict() for k in self.other_directives],
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Construct a MSONable AdfTask object from the JSON dict.

        Args:
            dct: A dict of saved attributes.

        Returns:
            An AdfTask object recovered from the JSON dict ``d``.
        """

        def _from_dict(_d):
            return AdfKey.from_dict(_d) if _d is not None else None

        operation = dct.get("operation")
        title = dct.get("title")
        basis_set = _from_dict(dct.get("basis_set"))
        xc = _from_dict(dct.get("xc"))
        units = _from_dict(dct.get("units"))
        scf = _from_dict(dct.get("scf"))
        others = [AdfKey.from_dict(o) for o in dct.get("others", [])]
        geo = _from_dict(dct.get("geo"))

        return cls(operation, basis_set, xc, title, units, geo.subkeys, scf, others)


class AdfInput:
    """A basic ADF input file writer."""

    def __init__(self, task):
        """
        Initialization method.

        Args:
            task (AdfTask): An ADF task.
        """
        self.task = task

    def write_file(self, molecule, inp_file):
        """Write an ADF input file.

        Args:
            molecule (Molecule): The molecule for this task.
            inpfile (str): The name where the input file will be saved.
        """
        mol_blocks = []
        atom_block = AdfKey("Atoms", options=["cartesian"])
        for site in molecule:
            atom_block.add_subkey(AdfKey(str(site.specie), list(site.coords)))
        mol_blocks.append(atom_block)

        if molecule.charge != 0:
            net_q = molecule.charge
            ab = molecule.spin_multiplicity - 1
            charge_block = AdfKey("Charge", [net_q, ab])
            mol_blocks.append(charge_block)
            if ab != 0:
                unres_block = AdfKey("Unrestricted")
                mol_blocks.append(unres_block)

        with open(inp_file, "w+", encoding="utf-8") as file:
            for block in mol_blocks:
                file.write(str(block) + "\n")
            file.write(str(self.task) + "\n")
            file.write("END INPUT")


class AdfOutput:
    """
    A basic ADF output file parser.

    Attributes:
        is_failed (bool): Whether the ADF job is failed.
        is_internal_crash (bool): Whether the job crashed.
            Please read 'TAPE13' of the ADF manual for more detail.
        error (str): The error description.
        run_type (str): The RunType of this ADF job. Possible options are:
            'SinglePoint', 'GeometryOptimization', 'AnalyticalFreq' and 'NUmericalFreq'.
        final_energy (float): The final molecule energy (a.u).
        final_structure (GMolecule): The final structure of the molecule.
        energies (Sized): The energy of each cycle.
        structures (Sized): The structure of each cycle If geometry optimization is performed.
        frequencies (array_like): The frequencies of the molecule.
        normal_modes (array_like): The normal modes of the molecule.
        freq_type (syr): Either 'Analytical' or 'Numerical'.
    """

    def __init__(self, filename):
        """
        Initialization method.

        Args:
            filename (str): The ADF output file to parse.
        """
        self.filename = filename
        self._parse()

    def _parse(self):
        """Parse the ADF outputs. There are two files: one is 'logfile', the other
        is the ADF output file. The final energy and structures are parsed from
        the 'logfile'. Frequencies and normal modes are parsed from the ADF
        output file.
        """
        workdir = os.path.dirname(self.filename)
        logfile = f"{workdir}/logfile"
        for ext in ("", ".gz", ".bz2"):
            if os.path.isfile(f"{logfile}{ext}"):
                logfile = f"{logfile}{ext}"
                break
        else:
            raise FileNotFoundError("The ADF logfile can not be accessed!")

        self.is_failed = False
        self.error = self.final_energy = self.final_structure = None
        self.energies = []
        self.structures = []
        self.frequencies = []
        self.normal_modes = self.freq_type = self.run_type = None
        self.is_internal_crash = False

        self._parse_logfile(logfile)
        if not self.is_failed and self.run_type != "SinglePoint":
            self._parse_adf_output()

    @staticmethod
    def _sites_to_mol(sites):
        """Get a Molecule object given a list of sites.

        Args:
            sites : A list of sites.

        Returns:
            mol (Molecule): A Molecule object.
        """
        return Molecule([site[0] for site in sites], [site[1] for site in sites])

    def _parse_logfile(self, logfile):
        """Parse the formatted logfile."""
        cycle_patt = re.compile(r"Coordinates\sin\sGeometry\sCycle\s(\d+)")
        coord_patt = re.compile(r"\s+([0-9]+)\.([A-Za-z]+)" + 3 * r"\s+([-\.0-9]+)")
        energy_patt = re.compile(r"<.*>\s<.*>\s+current\senergy\s+([-\.0-9]+)\sHartree")
        final_energy_patt = re.compile(r"<.*>\s<.*>\s+Bond\sEnergy\s+([-\.0-9]+)\sa\.u\.")
        error_patt = re.compile(r"<.*>\s<.*>\s+ERROR\sDETECTED:\s(.*)")
        run_type_patt = re.compile(r"<.*>\s<.*>\s+RunType\s+:\s(.*)")
        end_patt = re.compile(r"<.*>\s<.*>\s+END")
        parse_cycle = False
        sites = []
        last_cycle = -1
        parse_final = False

        # Stop parsing the logfile is this job is not terminated successfully.
        # The last non-empty line of the logfile must match the end pattern.
        # Otherwise the job has some internal failure. The TAPE13 part of the
        # ADF manual has a detailed explanation.
        for line in reverse_readfile(logfile):
            if line.strip() == "":
                continue
            if end_patt.search(line) is None:
                self.is_internal_crash = True
                self.error = "Internal crash. TAPE13 is generated!"
                self.is_failed = True
                return
            break

        with open(logfile, encoding="utf-8") as file:
            for line in file:
                if match := error_patt.search(line):
                    self.is_failed = True
                    self.error = match[1]
                    break

                if self.run_type is None:
                    if match := run_type_patt.search(line):
                        if match[1] == "FREQUENCIES":
                            self.freq_type = "Numerical"
                            self.run_type = "NumericalFreq"
                        elif match[1] == "GEOMETRY OPTIMIZATION":
                            self.run_type = "GeometryOptimization"
                        elif match[1] == "CREATE":
                            self.run_type = None
                        elif match[1] == "SINGLE POINT":
                            self.run_type = "SinglePoint"
                        else:
                            raise AdfOutputError("Undefined Runtype!")

                elif self.run_type == "SinglePoint":
                    if match := coord_patt.search(line):
                        sites.append([match.groups()[0], list(map(float, match.groups()[2:]))])
                    elif match := final_energy_patt.search(line):
                        self.final_energy = float(match[1])
                        self.final_structure = self._sites_to_mol(sites)

                elif self.run_type == "GeometryOptimization":
                    if match := cycle_patt.search(line):
                        cycle = int(match[1])
                        if cycle <= 0:
                            raise AdfOutputError(f"Wrong {cycle=}")
                        if cycle > last_cycle:
                            parse_cycle = True
                            last_cycle = cycle
                        else:
                            parse_final = True
                    elif parse_cycle:
                        if match := coord_patt.search(line):
                            sites.append(
                                [
                                    match.groups()[1],
                                    list(map(float, match.groups()[2:])),
                                ]
                            )
                        elif match := energy_patt.search(line):
                            self.energies.append(float(match[1]))
                            mol = self._sites_to_mol(sites)
                            self.structures.append(mol)
                            parse_cycle = False
                            sites = []
                    elif parse_final:
                        if match := final_energy_patt.search(line):
                            self.final_energy = float(match[1])

                elif self.run_type == "NumericalFreq":
                    break

        if not self.is_failed:
            if self.run_type == "GeometryOptimization":
                if len(self.structures) > 0:
                    self.final_structure = self.structures[-1]
                if self.final_energy is None:
                    raise AdfOutputError("The final energy can not be read!")
            elif self.run_type == "SinglePoint":
                if self.final_structure is None:
                    raise AdfOutputError("The final structure is missing!")
                if self.final_energy is None:
                    raise AdfOutputError("The final energy can not be read!")

    def _parse_adf_output(self):
        """Parse the standard ADF output file."""
        numerical_freq_patt = re.compile(r"\s+\*\s+F\sR\sE\sQ\sU\sE\sN\sC\sI\sE\sS\s+\*")
        analytic_freq_patt = re.compile(r"\s+\*\s+F\sR\sE\sQ\sU\sE\sN\sC\sY\s+A\sN\sA\sL\sY\sS\sI\sS\s+\*")
        freq_on_patt = re.compile(r"Vibrations\sand\sNormal\sModes\s+\*+.*\*+")
        freq_off_patt = re.compile(r"List\sof\sAll\sFrequencies:")
        mode_patt = re.compile(r"\s+(\d+)\.([A-Za-z]+)\s+(.*)")
        coord_patt = re.compile(r"\s+(\d+)\s+([A-Za-z]+)" + 6 * r"\s+([0-9\.-]+)")
        coord_on_patt = re.compile(r"\s+\*\s+R\sU\sN\s+T\sY\sP\sE\s:\sFREQUENCIES\s+\*")
        parse_freq = parse_mode = False
        n_next = n_strike = 0
        sites = []

        self.frequencies = []
        self.normal_modes = []

        if self.final_structure is None:
            find_structure = True
            parse_coord = False
            n_atoms = 0
        else:
            find_structure = parse_coord = False
            n_atoms = len(self.final_structure)

        with open(self.filename, encoding="utf-8") as file:
            for line in file:
                if self.run_type == "NumericalFreq" and find_structure:
                    if not parse_coord:
                        if match := coord_on_patt.search(line):
                            parse_coord = True
                    elif match := coord_patt.search(line):
                        sites.append([match[2], list(map(float, match.groups()[2:5]))])
                        n_strike += 1
                    elif n_strike > 0:
                        find_structure = False
                        self.final_structure = self._sites_to_mol(sites)
                        n_atoms = len(self.final_structure)

                elif self.freq_type is None:
                    if numerical_freq_patt.search(line):
                        self.freq_type = "Numerical"
                    elif analytic_freq_patt.search(line):
                        self.freq_type = "Analytical"
                        self.run_type = "AnalyticalFreq"

                elif freq_on_patt.search(line):
                    parse_freq = True

                elif parse_freq:
                    if freq_off_patt.search(line):
                        break
                    el = line.strip().split()
                    if 1 <= len(el) <= 3 and "." in line:
                        n_next = len(el)
                        parse_mode = True
                        parse_freq = False
                        self.frequencies.extend(map(float, el))
                        self.normal_modes.extend([] for _ in range(n_next))

                elif parse_mode and (match := mode_patt.search(line)):
                    v = list(chunks(map(float, match[3].split()), 3))
                    if len(v) != n_next:
                        raise AdfOutputError("Odd Error!")
                    for i, k in enumerate(range(-n_next, 0)):
                        self.normal_modes[k].extend(v[i])
                    if int(match[1]) == n_atoms:
                        parse_freq = True
                        parse_mode = False
        if isinstance(self.final_structure, list):
            self.final_structure = self._sites_to_mol(self.final_structure)

        if self.freq_type is not None:
            if len(self.frequencies) != len(self.normal_modes):
                raise AdfOutputError("The number of normal modes is wrong!")
            if len(self.normal_modes[0]) != n_atoms * 3:
                raise AdfOutputError("The dimensions of the modes are wrong!")
