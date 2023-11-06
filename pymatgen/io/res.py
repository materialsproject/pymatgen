"""
Provides parsing and read/write support for ShelX .res files as produced by the AIRSS code.

Converting from and back to pymatgen objects is expected to be reversible, i.e. you
should get the same Structure or ComputedStructureEntry back. On the other hand, converting
from and back to a string/file is not guaranteed to be reversible, i.e. a diff on the output
would not be empty. The difference should be limited to whitespace, float precision, and the
REM entries.

"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Literal

import dateutil.parser  # type: ignore[import]
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.core import ParseError

if TYPE_CHECKING:
    from collections.abc import Iterator
    from datetime import date

    from pymatgen.core.trajectory import Vector3D


@dataclass(frozen=True)
class AirssTITL:
    seed: str
    pressure: float
    volume: float
    energy: float
    integrated_spin_density: float
    integrated_absolute_spin_density: float
    spacegroup_label: str
    appearances: int

    def __str__(self) -> str:
        return (
            f"TITL {self.seed:s} {self.pressure:.2f} {self.volume:.4f} {self.energy:.5f} "
            f"{self.integrated_spin_density:f} {self.integrated_absolute_spin_density:f} ({self.spacegroup_label:s}) "
            f"n - {self.appearances}"
        )


@dataclass(frozen=True)
class ResCELL:
    unknown_field_1: float
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float

    def __str__(self) -> str:
        return (
            f"CELL {self.unknown_field_1:.5f} {self.a:.5f} {self.b:.5f} {self.c:.5f} "
            f"{self.alpha:.5f} {self.beta:.5f} {self.gamma:.5f}"
        )


@dataclass(frozen=True)
class Ion:
    specie: str
    specie_num: int
    pos: Vector3D
    occupancy: float
    spin: float | None

    def __str__(self) -> str:
        if self.spin is None:
            ion_fmt = "{:<7s}{:<2d} {:.8f} {:.8f} {:.8f} {:f}"
            return ion_fmt.format(self.specie, self.specie_num, *self.pos, self.occupancy)

        ion_fmt = "{:<7s}{:<2d} {:.8f} {:.8f} {:.8f} {:f} {:5.2f}"
        return ion_fmt.format(self.specie, self.specie_num, *self.pos, self.occupancy, self.spin)


@dataclass(frozen=True)
class ResSFAC:
    species: set[str]
    ions: list[Ion]

    def __str__(self) -> str:
        species = " ".join(f"{specie:<2s}" for specie in self.species)
        ions = "\n".join(map(str, self.ions))
        return f"SFAC {species}\n{ions}\nEND"


@dataclass(frozen=True)
class Res:
    """Representation for the data in a res file."""

    TITL: AirssTITL | None
    REMS: list[str]
    CELL: ResCELL
    SFAC: ResSFAC

    def __str__(self) -> str:
        return "\n".join(
            [
                "TITL" if self.TITL is None else str(self.TITL),
                "\n".join(f"REM {rem}" for rem in self.REMS),
                str(self.CELL),
                "LATT -1",
                str(self.SFAC),
            ]
        )


class ResParseError(ParseError):
    """This exception indicates a problem was encountered during parsing due to unexpected formatting."""


class ResError(ValueError):
    """
    This exception indicates a problem was encountered while trying to retrieve a value or
    perform an action that a provider for the res file does not support.
    """


class ResParser:
    """Parser for the ShelX res file."""

    def __init__(self):
        self.line: int = 0
        self.filename: str | None = None
        self.source: str = ""

    def _parse_titl(self, line: str) -> AirssTITL | None:
        """Parses the TITL entry. Checks for AIRSS values in the entry."""
        fields = line.split(maxsplit=6)
        if len(fields) >= 6:
            # this is probably an AIRSS res file
            seed, pressure, volume, energy, spin, abs_spin = fields[:6]
            spg, nap = "P1", "1"
            if len(fields) == 7:
                rest = fields[6]
                # just extract spacegroup and num appearances from the rest
                lp = rest.find("(")
                rp = rest.find(")")
                spg = rest[lp + 1 : rp]
                nmin = rest.find("n -")
                nap = rest[nmin + 4 :]
            return AirssTITL(
                seed, float(pressure), float(volume), float(energy), float(spin), float(abs_spin), spg, int(nap)
            )
        # there should at least be the first 6 fields if it's an AIRSS res file
        # if it doesn't have them, then just stop looking
        return None

    def _parse_cell(self, line: str) -> ResCELL:
        """Parses the CELL entry."""
        fields = line.split()
        if len(fields) != 7:
            raise ResParseError(f"Failed to parse CELL {line=}, expected 7 fields.")
        field_1, a, b, c, alpha, beta, gamma = map(float, fields)
        return ResCELL(field_1, a, b, c, alpha, beta, gamma)

    def _parse_ion(self, line: str) -> Ion:
        """Parses entries in the SFAC block."""
        fields = line.split()
        if len(fields) == 6:
            spin = None
        elif len(fields) == 7:
            spin = float(fields[-1])
        else:
            raise ResParseError(f"Failed to parse ion entry {line}, expected 6 or 7 fields.")
        specie = fields[0]
        specie_num = int(fields[1])
        x, y, z, occ = map(float, fields[2:6])
        return Ion(specie, specie_num, (x, y, z), occ, spin)

    def _parse_sfac(self, line: str, it: Iterator[str]) -> ResSFAC:
        """Parses the SFAC block."""
        species = set(line.split())
        ions = []
        try:
            while True:
                line = next(it)
                if line == "END":
                    break
                ions.append(self._parse_ion(line))
        except StopIteration:
            raise ResParseError("Encountered end of file before END tag at end of SFAC block.")
        return ResSFAC(species, ions)

    def _parse_txt(self) -> Res:
        """Parses the text of the file."""
        _REMS: list[str] = []
        _TITL: AirssTITL | None = None
        _CELL: ResCELL | None = None
        _SFAC: ResSFAC | None = None

        txt = self.source
        it = iter(txt.splitlines())
        try:
            while True:
                line = next(it)
                self.line += 1
                split = line.split(maxsplit=1)
                splits = len(split)
                if splits == 0:
                    continue
                if splits == 1:
                    first, rest = *split, ""
                else:
                    first, rest = split
                if first == "TITL":
                    _TITL = self._parse_titl(rest)
                elif first == "REM":
                    _REMS.append(rest)
                elif first == "CELL":
                    _CELL = self._parse_cell(rest)
                elif first == "LATT":
                    pass  # ignore
                elif first == "SFAC":
                    _SFAC = self._parse_sfac(rest, it)
                else:
                    raise Warning(f"Skipping {line=}, tag {first} not recognized.")
        except StopIteration:
            pass
        if _CELL is None or _SFAC is None:
            raise ResParseError("Did not encounter CELL or SFAC entry when parsing.")
        return Res(_TITL, _REMS, _CELL, _SFAC)

    @classmethod
    def _parse_str(cls, source: str) -> Res:
        """Parses the res file as a string."""
        self = cls()
        self.source = source
        return self._parse_txt()

    @classmethod
    def _parse_file(cls, filename: str) -> Res:
        """Parses the res file as a file."""
        self = cls()
        self.filename = filename
        with zopen(filename, "r") as file:
            self.source = file.read()
            return self._parse_txt()


class ResWriter:
    """This class provides a means to write a Structure or ComputedStructureEntry to a res file."""

    @classmethod
    def _cell_from_lattice(cls, lattice: Lattice) -> ResCELL:
        """Produce CELL entry from a pymatgen Lattice."""
        return ResCELL(1.0, lattice.a, lattice.b, lattice.c, lattice.alpha, lattice.beta, lattice.gamma)

    @classmethod
    def _ions_from_sites(cls, sites: list[PeriodicSite]) -> list[Ion]:
        """Produce a list of entries for a SFAC block from a list of pymatgen PeriodicSite."""
        ions: list[Ion] = []
        i = 0
        for site in sites:
            for specie, occ in site.species.items():
                i += 1
                x, y, z = map(float, site.frac_coords)
                spin = site.properties.get("magmom")
                spin = spin and float(spin)
                ions.append(Ion(specie, i, (x, y, z), occ, spin))
        return ions

    @classmethod
    def _sfac_from_sites(cls, sites: list[PeriodicSite]) -> ResSFAC:
        """Produce a SFAC block from a list of pymatgen PeriodicSite."""
        ions = cls._ions_from_sites(sites)
        species = {ion.specie for ion in ions}
        return ResSFAC(species, ions)

    @classmethod
    def _res_from_structure(cls, structure: Structure) -> Res:
        """Produce a res file structure from a pymatgen Structure."""
        return Res(None, [], cls._cell_from_lattice(structure.lattice), cls._sfac_from_sites(list(structure)))

    @classmethod
    def _res_from_entry(cls, entry: ComputedStructureEntry) -> Res:
        """Produce a res file structure from a pymatgen ComputedStructureEntry."""
        seed = entry.data.get("seed") or str(hash(entry))
        pres = float(entry.data.get("pressure", 0))
        isd = float(entry.data.get("isd", 0))
        iasd = float(entry.data.get("iasd", 0))
        spg, _ = entry.structure.get_space_group_info()
        rems = [str(x) for x in entry.data.get("rems", [])]
        return Res(
            AirssTITL(seed, pres, entry.structure.volume, entry.energy, isd, iasd, spg, 1),
            rems,
            cls._cell_from_lattice(entry.structure.lattice),
            cls._sfac_from_sites(list(entry.structure)),
        )

    def __init__(self, entry: Structure | ComputedStructureEntry):
        """This class can be constructed from either a pymatgen Structure or ComputedStructureEntry object."""
        func: Callable[[Structure], Res] | Callable[[ComputedStructureEntry], Res]
        func = self._res_from_structure
        if isinstance(entry, ComputedStructureEntry):
            func = self._res_from_entry
        self._res = func(entry)

    def __str__(self):
        return str(self._res)

    @property
    def string(self) -> str:
        """The contents of the res file."""
        return str(self)

    def write(self, filename: str) -> None:
        """Write the res data to a file."""
        with zopen(filename, "w") as file:
            file.write(str(self))


class ResProvider(MSONable):
    """Provides access to elements of the res file in the form of familiar pymatgen objects."""

    def __init__(self, res: Res) -> None:
        """The :func:`from_str` and :func:`from_file` methods should be used instead of constructing this directly."""
        self._res = res

    @classmethod
    def _site_spin(cls, spin: float | None) -> dict[str, float] | None:
        """Check and return a dict with the site spin. Return None if spin is None."""
        if spin is None:
            return None
        return {"magmom": spin}

    @classmethod
    def from_str(cls, string: str) -> ResProvider:
        """Construct a Provider from a string."""
        return cls(ResParser._parse_str(string))

    @classmethod
    def from_file(cls, filename: str) -> ResProvider:
        """Construct a Provider from a file."""
        return cls(ResParser._parse_file(filename))

    @property
    def rems(self) -> list[str]:
        """The full list of REM entries contained within the res file."""
        return self._res.REMS.copy()

    @property
    def lattice(self) -> Lattice:
        """Construct a Lattice from the res file."""
        cell = self._res.CELL
        return Lattice.from_parameters(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)

    @property
    def sites(self) -> list[PeriodicSite]:
        """Construct a list of PeriodicSites from the res file."""
        sfac_tag = self._res.SFAC
        return [
            PeriodicSite(ion.specie, ion.pos, self.lattice, properties=self._site_spin(ion.spin))
            for ion in sfac_tag.ions
        ]

    @property
    def structure(self) -> Structure:
        """Construct a Structure from the res file."""
        return Structure.from_sites(self.sites)


class AirssProvider(ResProvider):
    """
    Provides access to the res file as does ResProvider. This class additionally provides
    access to fields in the TITL entry and various other fields found in the REM entries
    that AIRSS puts in the file. Values in the TITL entry that AIRSS could not get end up as 0.
    If the TITL entry is malformed, empty, or missing then attempting to construct this class
    from a res file will raise a ResError.

    While AIRSS supports a number of geometry and energy solvers, CASTEP is the default. As such,
    fetching the information from the REM entries is only supported if AIRSS was used with CASTEP.
    The other properties that get put in the TITL should still be accessible even if CASTEP was
    not used.

    The :attr:`parse_rems` attribute controls whether functions that fail to retrieve information
    from the REM entries should return ``None``. If this is set to ``"strict"``,
    then a ParseError may be raised, but the return value will not be ``None``.
    If it is set to ``"gentle"``, then ``None`` will be returned instead of raising an
    exception. This setting applies to all methods of this class that are typed to return
    an Optional type. Default is ``"gentle"``.
    """

    _date_fmt = re.compile(r"[MTWFS][a-z]{2}, (\d{2}) ([A-Z][a-z]{2}) (\d{4}) (\d{2}):(\d{2}):(\d{2}) ([+-]?\d{4})")

    def __init__(self, res: Res, parse_rems: Literal["gentle", "strict"] = "gentle"):
        """The :func:`from_str` and :func:`from_file` methods should be used instead of constructing this directly."""
        super().__init__(res)
        if self._res.TITL is None:
            raise ResError(f"{self.__class__} can only be constructed from a res file with a valid TITL entry.")
        if parse_rems not in ["gentle", "strict"]:
            raise ValueError(f"{parse_rems} not valid, must be either 'gentle' or 'strict'.")
        self._TITL = self._res.TITL  # alias for the object so it is guarded by the None check
        self.parse_rems = parse_rems

    @classmethod
    def from_str(cls, string: str, parse_rems: Literal["gentle", "strict"] = "gentle") -> AirssProvider:
        """Construct a Provider from a string."""
        return cls(ResParser._parse_str(string), parse_rems)

    @classmethod
    def from_file(cls, filename: str, parse_rems: Literal["gentle", "strict"] = "gentle") -> AirssProvider:
        """Construct a Provider from a file."""
        return cls(ResParser._parse_file(filename), parse_rems)

    @classmethod
    def _parse_date(cls, string: str) -> date:
        """Parses a date from a string where the date is in the format typically used by CASTEP."""
        match = cls._date_fmt.search(string)
        if match is None:
            raise ResParseError(f"Could not parse the date from {string=}.")
        date_string = match.group(0)
        return dateutil.parser.parse(date_string)

    def _raise_or_none(self, err: ResParseError) -> None:
        if self.parse_rems != "strict":
            return
        raise err

    def get_run_start_info(self) -> tuple[date, str] | None:
        """
        Retrieves the run start date and the path it was started in from the REM entries.

        Returns:
            (date, path)
        """
        for rem in self._res.REMS:
            if rem.strip().startswith("Run started:"):
                date = self._parse_date(rem)
                path = rem.split()[-1]
                return date, path
        self._raise_or_none(ResParseError("Could not find run started information."))
        return None

    def get_castep_version(self) -> str | None:
        """
        Retrieves the version of CASTEP that the res file was computed with from the REM entries.

        Returns:
            version string
        """
        for rem in self._res.REMS:
            if rem.strip().startswith("CASTEP"):
                srem = rem.split()
                return srem[1][:-1]
        self._raise_or_none(ResParseError("No CASTEP version found in REM"))
        return None

    def get_func_rel_disp(self) -> tuple[str, str, str] | None:
        """
        Retrieves the functional, relativity scheme, and dispersion correction from the REM entries.

        Returns:
            (functional, relativity, dispersion)
        """
        for rem in self._res.REMS:
            if rem.strip().startswith("Functional"):
                srem = rem.split()
                return " ".join(srem[1:4]), srem[5], srem[7]
        self._raise_or_none(ResParseError("Could not find functional, relativity, and dispersion."))
        return None

    def get_cut_grid_gmax_fsbc(self) -> tuple[float, float, float, str] | None:
        """
        Retrieves the cut-off energy, grid scale, Gmax, and finite basis set correction setting
        from the REM entries.

        Returns:
            (cut-off, grid scale, Gmax, fsbc)
        """
        for rem in self._res.REMS:
            if rem.strip().startswith("Cut-off"):
                srem = rem.split()
                return float(srem[1]), float(srem[5]), float(srem[7]), srem[10]
        self._raise_or_none(ResParseError("Could not find line with cut-off energy."))
        return None

    def get_mpgrid_offset_nkpts_spacing(self) -> tuple[tuple[int, int, int], Vector3D, int, float] | None:
        """
        Retrieves the MP grid, the grid offsets, number of kpoints, and maximum kpoint spacing.

        Returns:
            (MP grid), (offsets), No. kpts, max spacing)
        """
        for rem in self._res.REMS:
            if rem.strip().startswith("MP grid"):
                srem = rem.split()
                p, q, r = map(int, srem[2:5])
                po, qo, ro = map(float, srem[6:9])
                return (p, q, r), (po, qo, ro), int(srem[11]), float(srem[13])
        self._raise_or_none(ResParseError("Could not find line with MP grid."))
        return None

    def get_airss_version(self) -> tuple[str, date] | None:
        """
        Retrieves the version of AIRSS that was used along with the build date (not compile date).

        Return:
            (version string, date)
        """
        for rem in self._res.REMS:
            if rem.strip().startswith("AIRSS Version"):
                date = self._parse_date(rem)
                v = rem.split()[2]
                return v, date
        self._raise_or_none(ResParseError("Could not find line with AIRSS version."))
        return None

    def _get_compiler(self):
        raise NotImplementedError

    def _get_compile_options(self):
        raise NotImplementedError

    def _get_rng_seeds(self):
        raise NotImplementedError

    def get_pspots(self) -> dict[str, str]:
        """
        Retrieves the OTFG pseudopotential string that can be used to generate the
        pseudopotentials used in the calculation.

        Returns:
            dict[specie, potential]
        """
        pspots: dict[str, str] = {}
        for rem in self._res.REMS:
            srem = rem.split()
            if len(srem) == 2 and Element.is_valid_symbol(srem[0]):
                k, v = srem
                pspots[k] = v
        return pspots

    @property
    def seed(self) -> str:
        """The seed name, typically also the name of the res file."""
        return self._TITL.seed

    @property
    def pressure(self) -> float:
        """Pressure for the run. This is in GPa if CASTEP was used."""
        return self._TITL.pressure

    @property
    def volume(self) -> float:
        """Volume of the structure. This is in cubic Angstroms if CASTEP was used."""
        return self._TITL.volume

    @property
    def energy(self) -> float:
        """Energy of the structure. With CASTEP, this is usually the enthalpy and is in eV."""
        return self._TITL.energy

    @property
    def integrated_spin_density(self) -> float:
        """Corresponds to the last ``Integrated Spin Density`` in the CASTEP file."""
        return self._TITL.integrated_spin_density

    @property
    def integrated_absolute_spin_density(self) -> float:
        """Corresponds to the last ``Integrated |Spin Density|`` in the CASTEP file."""
        return self._TITL.integrated_absolute_spin_density

    @property
    def spacegroup_label(self) -> str:
        """
        The Hermann-Mauguin notation of the spacegroup with ascii characters.
        So no. 225 would be Fm-3m, and no. 194 would be P6_3/mmc.
        """
        return self._TITL.spacegroup_label

    @property
    def appearances(self) -> int:
        """
        This is sometimes the number of times a structure was found in an AIRSS search.
        Using the cryan tool that comes with AIRSS may be a better approach than relying
        on this property.
        """
        return self._TITL.appearances

    @property
    def entry(self) -> ComputedStructureEntry:
        """Get this res file as a ComputedStructureEntry."""
        return ComputedStructureEntry(self.structure, self.energy, data={"rems": self.rems})

    def as_dict(self, verbose: bool = True) -> dict[str, Any]:
        """Get dict with title fields, structure and rems of this AirssProvider."""
        if verbose:
            return super().as_dict()
        return dict(**vars(self._res.TITL), structure=self.structure.as_dict(), rems=self.rems)


class ResIO:
    """
    Class providing convenience methods for converting a Structure or ComputedStructureEntry
    to/from a string or file in the res format as used by AIRSS.

    Note: Converting from and back to pymatgen objects is expected to be reversible, i.e. you
    should get the same Structure or ComputedStructureEntry back. On the other hand, converting
    from and back to a string/file is not guaranteed to be reversible, i.e. a diff on the output
    would not be empty. The difference should be limited to whitespace, float precision, and the
    REM entries.

    If the TITL entry doesn't exist or is malformed or empty, then you can only get
    a Structure. Attempting to get an Entry will raise a ResError.
    """

    @classmethod
    def structure_from_str(cls, string: str) -> Structure:
        """Produces a pymatgen Structure from contents of a res file."""
        return ResProvider.from_str(string).structure

    @classmethod
    def structure_from_file(cls, filename: str) -> Structure:
        """Produces a pymatgen Structure from a res file."""
        return ResProvider.from_file(filename).structure

    @classmethod
    def structure_to_str(cls, structure: Structure) -> str:
        """Produce the contents of a res file from a pymatgen Structure."""
        return str(ResWriter(structure))

    @classmethod
    def structure_to_file(cls, structure: Structure, filename: str) -> None:
        """Write a pymatgen Structure to a res file."""
        return ResWriter(structure).write(filename)

    @classmethod
    def entry_from_str(cls, string: str) -> ComputedStructureEntry:
        """Produce a pymatgen ComputedStructureEntry from contents of a res file."""
        return AirssProvider.from_str(string).entry

    @classmethod
    def entry_from_file(cls, filename: str) -> ComputedStructureEntry:
        """Produce a pymatgen ComputedStructureEntry from a res file."""
        return AirssProvider.from_file(filename).entry

    @classmethod
    def entry_to_str(cls, entry: ComputedStructureEntry) -> str:
        """Produce the contents of a res file from a pymatgen ComputedStructureEntry."""
        return str(ResWriter(entry))

    @classmethod
    def entry_to_file(cls, entry: ComputedStructureEntry, filename: str) -> None:
        """Write a pymatgen ComputedStructureEntry to a res file."""
        return ResWriter(entry).write(filename)
