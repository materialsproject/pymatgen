"""
IO for ADF files.
"""

import os
import re
from typing import Generator

from monty.io import reverse_readline
from monty.itertools import chunks
from monty.json import MSONable

from pymatgen.core.structure import Molecule

__author__ = "Xin Chen, chenxin13@mails.tsinghua.edu.cn"


def is_numeric(s):
    """
    Return True is the string ``s`` is a numeric string.

    Parameters
    ----------
    s : str
        A string.

    Returns
    -------
    res : bool
        If True, ``s`` is a numeric string and can be converted to an int or a
        float. Otherwise False will be returned.

    """
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True


def iterlines(s: str) -> Generator[str, None, None]:
    r"""A generator form of s.split('\n') for reducing memory overhead.

    Args:
        s (str): A multi-line string.

    Yields:
        str: line
    """
    prevnl = -1
    while True:
        nextnl = s.find("\n", prevnl + 1)
        if nextnl < 0:
            yield s[(prevnl + 1) :]
            break
        yield s[(prevnl + 1) : nextnl]
        prevnl = nextnl


class AdfInputError(Exception):
    """
    The default error class for ADF.
    """


class AdfOutputError(Exception):
    """
    The default error class for errors raised by ``AdfOutput``.
    """


class AdfKey(MSONable):
    """
    The basic input unit for ADF. A key is a string of characters that does not
    contain a delimiter (blank, comma or equal sign). A key may have multiple
    subkeys and a set of options.
    """

    block_keys = {
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
    }
    sub_keys = {"AtomDepQuality"}

    # Full blocks are blocks that must have an 'END'.
    _full_blocks = {"GEOMETRY", "SCF", "UNITS", "BASIS", "ANALYTICALFREQ"}

    def __init__(self, name, options=None, subkeys=None):
        """
        Initialization method.

        Parameters
        ----------
        name : str
            The name of this key.
        options : Sized
            The options for this key. Each element can be a primitive object or
            a tuple/list with two elements: the first is the name and the second
            is a primitive object.
        subkeys : Sized
            The subkeys for this key.

        Raises
        ------
        ValueError
            If elements in ``subkeys`` are not ``AdfKey`` objects.

        """
        self.name = name
        self.options = options if options is not None else []
        self.subkeys = subkeys if subkeys is not None else []
        if len(self.subkeys) > 0:
            for k in subkeys:
                if not isinstance(k, AdfKey):
                    raise ValueError("Not all subkeys are ``AdfKey`` objects!")
        self._sized_op = None
        if len(self.options) > 0:
            self._sized_op = isinstance(self.options[0], (list, tuple))

    def _options_string(self):
        """
        Return the option string.
        """
        if len(self.options) > 0:
            s = ""
            for op in self.options:
                if self._sized_op:
                    s += f"{op[0]}={op[1]} "
                else:
                    s += f"{op} "
            return s.strip()
        return ""

    def is_block_key(self):
        """
        Return True if this key is a block key.
        """
        return bool(self.name.upper() in self.block_keys)

    @property
    def key(self):
        """
        Return the name of this key. If this is a block key, the name will be
        converted to upper cases.
        """
        if self.is_block_key():
            return self.name.upper()
        return self.name

    def __str__(self):
        """
        Return the string representation of this ``AdfKey``.

        Notes
        -----
        If this key is 'Atoms' and the coordinates are in Cartesian form, a
        different string format will be used.

        """
        s = f"{self.key}"
        if len(self.options) > 0:
            s += f" {self._options_string()}"
        s += "\n"
        if len(self.subkeys) > 0:
            if self.key.lower() == "atoms":
                for subkey in self.subkeys:
                    s += (
                        f"{subkey.name:2s}  {subkey.options[0]: 14.8f}"
                        f"    {subkey.options[1]: 14.8f}    {subkey.options[2]: 14.8f}\n"
                    )
            else:
                for subkey in self.subkeys:
                    s += str(subkey)
            if self.is_block_key():
                s += "END\n"
            else:
                s += "subend\n"
        elif self.key.upper() in self._full_blocks:
            s += "END\n"
        return s

    def __eq__(self, other):
        if not isinstance(other, AdfKey):
            return False
        return str(self) == str(other)

    def has_subkey(self, subkey):
        """
        Return True if this AdfKey contains the given subkey.

        Parameters
        ----------
        subkey : str or AdfKey
            A key name or an AdfKey object.

        Returns
        -------
        has : bool
            True if this key contains the given key. Otherwise False.

        """
        if isinstance(subkey, str):
            key = subkey
        elif isinstance(subkey, AdfKey):
            key = subkey.key
        else:
            raise ValueError("The subkey should be an AdfKey or a string!")
        if len(self.subkeys) > 0:
            if key in map(lambda k: k.key, self.subkeys):
                return True
        return False

    def add_subkey(self, subkey):
        """
        Add a new subkey to this key.

        Parameters
        ----------
        subkey : AdfKey
            A new subkey.

        Notes
        -----
        Duplicate check will not be performed if this is an 'Atoms' block.

        """
        if self.key.lower() == "atoms" or not self.has_subkey(subkey):
            self.subkeys.append(subkey)

    def remove_subkey(self, subkey):
        """
        Remove the given subkey, if existed, from this AdfKey.

        Parameters
        ----------
        subkey : str or AdfKey
            The subkey to remove.

        """
        if len(self.subkeys) > 0:
            key = subkey if isinstance(subkey, str) else subkey.key
            for i, v in enumerate(self.subkeys):
                if v.key == key:
                    self.subkeys.pop(i)
                    break

    def add_option(self, option):
        """
        Add a new option to this key.

        Parameters
        ----------
        option : Sized or str or int or float
            A new option to add. This must have the same format with existing
            options.

        Raises
        ------
        TypeError
            If the format of the given ``option`` is different.

        """
        if len(self.options) == 0:
            self.options.append(option)
        else:
            sized_op = isinstance(option, (list, tuple))
            if self._sized_op != sized_op:
                raise TypeError("Option type is mismatched!")
            self.options.append(option)

    def remove_option(self, option):
        """
        Remove an option.

        Parameters
        ----------
        option : str or int
            The name (str) or index (int) of the option to remove.

        Raises
        ------
        TypeError
            If the option has a wrong type.

        """
        if len(self.options) > 0:
            if self._sized_op:
                if not isinstance(option, str):
                    raise TypeError("``option`` should be a name string!")
                for i, v in enumerate(self.options):
                    if v[0] == option:
                        self.options.pop(i)
                        break
            else:
                if not isinstance(option, int):
                    raise TypeError("``option`` should be an integer index!")
                self.options.pop(option)

    def has_option(self, option):
        """
        Return True if the option is included in this key.

        Parameters
        ----------
        option : str
            The option.

        Returns
        -------
        has : bool
            True if the option can be found. Otherwise False will be returned.

        """
        if len(self.options) == 0:
            return False
        for op in self.options:
            if (self._sized_op and op[0] == option) or (op == option):
                return True
        return False

    def as_dict(self):
        """
        A JSON-serializable dict representation of self.
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "name": self.name,
            "options": self.options,
        }
        if len(self.subkeys) > 0:
            subkeys = []
            for subkey in self.subkeys:
                subkeys.append(subkey.as_dict())
            d.update({"subkeys": subkeys})
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Construct a MSONable AdfKey object from the JSON dict.

        Parameters
        ----------
        d : dict
            A dict of saved attributes.

        Returns
        -------
        adfkey : AdfKey
            An AdfKey object recovered from the JSON dict ``d``.

        """
        key = d.get("name")
        options = d.get("options", None)
        subkey_list = d.get("subkeys", [])
        if len(subkey_list) > 0:
            subkeys = list(map(lambda k: AdfKey.from_dict(k), subkey_list))
        else:
            subkeys = None
        return cls(key, options, subkeys)

    @staticmethod
    def from_string(string):
        """
        Construct an AdfKey object from the string.

        Parameters
        ----------
        string : str
            A string.

        Returns
        -------
        adfkey : AdfKey
            An AdfKey object recovered from the string.

        Raises
        ------
        ValueError
            Currently nested subkeys are not supported. If ``subend`` was found
            a ValueError would be raised.

        Notes
        -----
        Only the first block key will be returned.

        """

        def is_float(s):
            return "." in s or "E" in s or "e" in s

        if string.find("\n") == -1:
            el = string.split()
            if len(el) > 1:
                if string.find("=") != -1:
                    options = list(map(lambda s: s.split("="), el[1:]))
                else:
                    options = el[1:]
                for i, op in enumerate(options):
                    if isinstance(op, list) and is_numeric(op[1]):
                        op[1] = float(op[1]) if is_float(op[1]) else int(op[1])
                    elif is_numeric(op):
                        options[i] = float(op) if is_float(op) else int(op)
            else:
                options = None
            return AdfKey(el[0], options)

        if string.find("subend") != -1:
            raise ValueError("Nested subkeys are not supported!")

        key = None
        for line in iterlines(string):
            if line == "":
                continue
            el = line.strip().split()
            if len(el) == 0:
                continue
            if el[0].upper() in AdfKey.block_keys:
                if key is None:
                    key = AdfKey.from_string(line)
                else:
                    return key
            elif el[0].upper() == "END":
                return key
            elif key is not None:
                key.add_subkey(AdfKey.from_string(line))
        else:
            raise Exception("IncompleteKey: 'END' is missing!")


class AdfTask(MSONable):
    """
    Basic task for ADF. All settings in this class are independent of molecules.

    Notes
    -----
    Unlike other quantum chemistry packages (NWChem, Gaussian, ...), ADF does
    not support calculating force/gradient.

    """

    operations = {
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

        Parameters
        ----------
        operation : str
            The target operation.
        basis_set : AdfKey
            The basis set definitions for this task. Defaults to 'DZ/Large'.
        xc : AdfKey
            The exchange-correlation functionals. Defaults to PBE.
        title : str
            The title of this ADF task.
        units : AdfKey
            The units. Defaults to Angstroms/Degree.
        geo_subkeys : Sized
            The subkeys for the block key 'GEOMETRY'.
        scf : AdfKey
            The scf options.
        other_directives : Sized
            User-defined directives.

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
        """
        Returns: Default basis set
        """
        return AdfKey.from_string("Basis\ntype DZ\ncore small\nEND")

    @staticmethod
    def get_default_scf():
        """
        Returns: ADF using default SCF.
        """
        return AdfKey.from_string("SCF\niterations 300\nEND")

    @staticmethod
    def get_default_geo():
        """
        Returns: ADFKey using default geometry.
        """
        return AdfKey.from_string("GEOMETRY SinglePoint\nEND")

    @staticmethod
    def get_default_xc():
        """
        Returns: ADFKey using default XC.
        """
        return AdfKey.from_string("XC\nGGA PBE\nEND")

    @staticmethod
    def get_default_units():
        """
        Returns: Default units.
        """
        return AdfKey.from_string("Units\nlength angstrom\nangle degree\nEnd")

    def _setup_task(self, geo_subkeys):
        """
        Setup the block 'Geometry' given subkeys and the task.

        Parameters
        ----------
        geo_subkeys : Sized
            User-defined subkeys for the block 'Geometry'.

        Notes
        -----
        Most of the run types of ADF are specified in the Geometry block except
        the 'AnalyticFreq'.

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
        s = f"""TITLE {self.title}\n
{str(self.units)}
{str(self.xc)}
{str(self.basis_set)}
{str(self.scf)}
{str(self.geo)}"""
        s += "\n"
        for block_key in self.other_directives:
            if not isinstance(block_key, AdfKey):
                raise ValueError(f"{block_key} is not an AdfKey!")
            s += str(block_key) + "\n"
        return s

    def as_dict(self):
        """
        A JSON-serializable dict representation of self.
        """
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
    def from_dict(cls, d):
        """
        Construct a MSONable AdfTask object from the JSON dict.

        Parameters
        ----------
        d : dict
            A dict of saved attributes.

        Returns
        -------
        task : AdfTask
            An AdfTask object recovered from the JSON dict ``d``.

        """

        def _from_dict(_d):
            return AdfKey.from_dict(_d) if _d is not None else None

        operation = d.get("operation")
        title = d.get("title")
        basis_set = _from_dict(d.get("basis_set"))
        xc = _from_dict(d.get("xc"))
        units = _from_dict(d.get("units"))
        scf = _from_dict(d.get("scf"))
        others = [AdfKey.from_dict(o) for o in d.get("others", [])]
        geo = _from_dict(d.get("geo"))

        return cls(operation, basis_set, xc, title, units, geo.subkeys, scf, others)


class AdfInput:
    """
    A basic ADF input file writer.
    """

    def __init__(self, task):
        """
        Initialization method.

        Parameters
        ----------
        task : AdfTask
            An ADF task.

        """
        self.task = task

    def write_file(self, molecule, inpfile):
        """
        Write an ADF input file.

        Parameters
        ----------
        molecule : Molecule
            The molecule for this task.
        inpfile : str
            The name where the input file will be saved.

        """

        mol_blocks = []
        atom_block = AdfKey("Atoms", options=["cartesian"])
        for site in molecule:
            atom_block.add_subkey(AdfKey(str(site.specie), list(site.coords)))
        mol_blocks.append(atom_block)

        if molecule.charge != 0:
            netq = molecule.charge
            ab = molecule.spin_multiplicity - 1
            charge_block = AdfKey("Charge", [netq, ab])
            mol_blocks.append(charge_block)
            if ab != 0:
                unres_block = AdfKey("Unrestricted")
                mol_blocks.append(unres_block)

        with open(inpfile, "w+") as f:
            for block in mol_blocks:
                f.write(str(block) + "\n")
            f.write(str(self.task) + "\n")
            f.write("END INPUT")


class AdfOutput:
    """
    A basic ADF output file parser.

    Attributes
    ----------
    is_failed : bool
        True is the ADF job is terminated without success. Otherwise False.
    is_internal_crash : bool
        True if the job is terminated with internal crash. Please read 'TAPE13'
        of the ADF manual for more detail.
    error : str
        The error description.
    run_type : str
        The RunType of this ADF job. Possible options are: 'SinglePoint',
        'GeometryOptimization', 'AnalyticalFreq' and 'NUmericalFreq'.
    final_energy : float
        The final molecule energy (a.u).
    final_structure : GMolecule
        The final structure of the molecule.
    energies : Sized
        The energy of each cycle.
    structures : Sized
        The structure of each cycle If geometry optimization is performed.
    frequencies : array_like
        The frequencies of the molecule.
    normal_modes : array_like
        The normal modes of the molecule.
    freq_type : str
        Either 'Analytical' or 'Numerical'.

    """

    def __init__(self, filename):
        """
        Initialization method.

        Parameters
        ----------
        filename : str
            The ADF output file to parse.

        """
        self.filename = filename
        self._parse()

    def _parse(self):
        """
        Parse the ADF outputs. There are two files: one is 'logfile', the other
        is the ADF output file. The final energy and structures are parsed from
        the 'logfile'. Frequencies and normal modes are parsed from the ADF
        output file.
        """
        workdir = os.path.dirname(self.filename)
        logfile = os.path.join(workdir, "logfile")
        if not os.path.isfile(logfile):
            raise OSError("The ADF logfile can not be accessed!")

        self.is_failed = False
        self.error = None
        self.final_energy = None
        self.final_structure = None
        self.energies = []
        self.structures = []
        self.frequencies = []
        self.normal_modes = None
        self.freq_type = None
        self.run_type = None
        self.is_internal_crash = False

        self._parse_logfile(logfile)
        if not self.is_failed and self.run_type != "SinglePoint":
            self._parse_adf_output()

    @staticmethod
    def _sites_to_mol(sites):
        """
        Return a ``Molecule`` object given a list of sites.

        Parameters
        ----------
        sites : list
            A list of sites.

        Returns
        -------
        mol : Molecule
            A ``Molecule`` object.

        """
        return Molecule([site[0] for site in sites], [site[1] for site in sites])

    def _parse_logfile(self, logfile):
        """
        Parse the formatted logfile.
        """

        cycle_patt = re.compile(r"Coordinates\sin\sGeometry\sCycle\s(\d+)")
        coord_patt = re.compile(r"\s+([0-9]+)\.([A-Za-z]+)" + 3 * r"\s+([-\.0-9]+)")
        energy_patt = re.compile(r"<.*>\s<.*>\s+current\senergy\s+([-\.0-9]+)\sHartree")
        final_energy_patt = re.compile(r"<.*>\s<.*>\s+Bond\sEnergy\s+([-\.0-9]+)\sa\.u\.")
        error_patt = re.compile(r"<.*>\s<.*>\s+ERROR\sDETECTED:\s(.*)")
        runtype_patt = re.compile(r"<.*>\s<.*>\s+RunType\s+:\s(.*)")
        end_patt = re.compile(r"<.*>\s<.*>\s+END")
        parse_cycle = False
        sites = []
        last_cycle = -1
        parse_final = False

        # Stop parsing the logfile is this job is not terminated successfully.
        # The last non-empty line of the logfile must match the end pattern.
        # Otherwise the job has some internal failure. The TAPE13 part of the
        # ADF manual has a detailed explanantion.
        with open(logfile) as f:
            for line in reverse_readline(f):
                if line == "":
                    continue
                if end_patt.search(line) is None:
                    self.is_internal_crash = True
                    self.error = "Internal crash. TAPE13 is generated!"
                    self.is_failed = True
                    return
                break

        with open(logfile) as f:
            for line in f:
                m = error_patt.search(line)
                if m:
                    self.is_failed = True
                    self.error = m.group(1)
                    break

                if self.run_type is None:
                    m = runtype_patt.search(line)
                    if m:
                        if m.group(1) == "FREQUENCIES":
                            self.freq_type = "Numerical"
                            self.run_type = "NumericalFreq"
                        elif m.group(1) == "GEOMETRY OPTIMIZATION":
                            self.run_type = "GeometryOptimization"
                        elif m.group(1) == "CREATE":
                            self.run_type = None
                        elif m.group(1) == "SINGLE POINT":
                            self.run_type = "SinglePoint"
                        else:
                            raise AdfOutputError("Undefined Runtype!")

                elif self.run_type == "SinglePoint":
                    m = coord_patt.search(line)
                    if m:
                        sites.append([m.groups()[0], list(map(float, m.groups()[2:]))])
                    else:
                        m = final_energy_patt.search(line)
                        if m:
                            self.final_energy = float(m.group(1))
                            self.final_structure = self._sites_to_mol(sites)

                elif self.run_type == "GeometryOptimization":
                    m = cycle_patt.search(line)
                    if m:
                        cycle = int(m.group(1))
                        if cycle <= 0:
                            raise AdfOutputError(f"Wrong cycle {cycle}")
                        if cycle > last_cycle:
                            parse_cycle = True
                            last_cycle = cycle
                        else:
                            parse_final = True
                    elif parse_cycle:
                        m = coord_patt.search(line)
                        if m:
                            sites.append([m.groups()[1], list(map(float, m.groups()[2:]))])
                        else:
                            m = energy_patt.search(line)
                            if m:
                                self.energies.append(float(m.group(1)))
                                mol = self._sites_to_mol(sites)
                                self.structures.append(mol)
                                parse_cycle = False
                                sites = []
                    elif parse_final:
                        m = final_energy_patt.search(line)
                        if m:
                            self.final_energy = float(m.group(1))

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
        """
        Parse the standard ADF output file.
        """
        numerical_freq_patt = re.compile(r"\s+\*\s+F\sR\sE\sQ\sU\sE\sN\sC\sI\sE\sS\s+\*")
        analytic_freq_patt = re.compile(r"\s+\*\s+F\sR\sE\sQ\sU\sE\sN\sC\sY\s+A\sN\sA\sL\sY\sS\sI\sS\s+\*")
        freq_on_patt = re.compile(r"Vibrations\sand\sNormal\sModes\s+\*+.*\*+")
        freq_off_patt = re.compile(r"List\sof\sAll\sFrequencies:")
        mode_patt = re.compile(r"\s+(\d+)\.([A-Za-z]+)\s+(.*)")
        coord_patt = re.compile(r"\s+(\d+)\s+([A-Za-z]+)" + 6 * r"\s+([0-9\.-]+)")
        coord_on_patt = re.compile(r"\s+\*\s+R\sU\sN\s+T\sY\sP\sE\s:\sFREQUENCIES\s+\*")
        parse_freq = False
        parse_mode = False
        nnext = 0
        nstrike = 0
        sites = []

        self.frequencies = []
        self.normal_modes = []

        if self.final_structure is None:
            find_structure = True
            parse_coord = False
            natoms = 0
        else:
            find_structure = False
            parse_coord = False
            natoms = self.final_structure.num_sites

        with open(self.filename) as f:
            for line in f:
                if self.run_type == "NumericalFreq" and find_structure:
                    if not parse_coord:
                        m = coord_on_patt.search(line)
                        if m:
                            parse_coord = True
                    else:
                        m = coord_patt.search(line)
                        if m:
                            sites.append([m.group(2), list(map(float, m.groups()[2:5]))])
                            nstrike += 1
                        elif nstrike > 0:
                            find_structure = False
                            self.final_structure = self._sites_to_mol(sites)
                            natoms = self.final_structure.num_sites

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
                    if 1 <= len(el) <= 3 and line.find(".") != -1:
                        nnext = len(el)
                        parse_mode = True
                        parse_freq = False
                        self.frequencies.extend(map(float, el))
                        for _ in range(nnext):
                            self.normal_modes.append([])

                elif parse_mode:
                    m = mode_patt.search(line)
                    if m:
                        v = list(chunks(map(float, m.group(3).split()), 3))
                        if len(v) != nnext:
                            raise AdfOutputError("Odd Error!")
                        for i, k in enumerate(range(-nnext, 0, 1)):
                            self.normal_modes[k].extend(v[i])
                        if int(m.group(1)) == natoms:
                            parse_freq = True
                            parse_mode = False
        if isinstance(self.final_structure, list):
            self.final_structure = self._sites_to_mol(self.final_structure)

        if self.freq_type is not None:
            if len(self.frequencies) != len(self.normal_modes):
                raise AdfOutputError("The number of normal modes is wrong!")
            if len(self.normal_modes[0]) != natoms * 3:
                raise AdfOutputError("The dimensions of the modes are wrong!")
