"""
Provides parsing and read/write support for ShelX .res files as produced by the AIRSS code.
"""

from dataclasses import dataclass
from typing import Iterator, List, Optional, Set, Tuple

from monty.io import zopen

from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry

__all__ = ["ResIO"]


@dataclass
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
        title_fmt = "TITL {:s} {:.2f} {:.4f} {:.5f} {:f} {:f} ({:s}) n - {:d}"
        return title_fmt.format(
            self.seed,
            self.pressure,
            self.volume,
            self.energy,
            self.integrated_spin_density,
            self.integrated_absolute_spin_density,
            self.spacegroup_label,
            self.appearances,
        )


@dataclass
class ResCELL:
    unknown_field_1: float
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float

    def __str__(self) -> str:
        cell_fmt = "CELL {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f}"
        return cell_fmt.format(self.unknown_field_1, self.a, self.b, self.c, self.alpha, self.beta, self.gamma)


@dataclass
class Ion:
    specie: str
    specie_num: int
    pos: Tuple[float, float, float]
    occupancy: float

    def __str__(self) -> str:
        ion_fmt = "{:<7s}{:<2d} {:.8f} {:.8f} {:.8f} {:f}"
        return ion_fmt.format(self.specie, self.specie_num, *self.pos, self.occupancy)


@dataclass
class ResSFAC:
    species: Set[str]
    ions: List[Ion]

    def __str__(self) -> str:
        sfac_fmt = "SFAC {species}\n" "{ions}\n" "END"
        return sfac_fmt.format(
            species=" ".join(f"{specie:<2s}" for specie in self.species), ions="\n".join(map(str, self.ions))
        )


@dataclass
class Res:
    """
    Representation for the data in a res file.
    """

    TITL: Optional[AirssTITL]
    REMS: List[str]
    CELL: ResCELL
    SFAC: ResSFAC

    def __str__(self) -> str:
        return "\n".join(
            [
                "TITL" if self.TITL is None else str(self.TITL),
                "\n".join(self.REMS),
                str(self.CELL),
                "LATT -1",
                str(self.SFAC),
            ]
        )


class ParseError(RuntimeError):
    ...


class ResParser:
    """Parser for the ShelX res file."""

    def __init__(self):
        self.line: int = 0
        self.filename: Optional[str] = None
        self.source: str = ""

    def _parse_titl(self, line: str) -> Optional[AirssTITL]:
        """Parses the TITL entry. Checks for airss values in the entry."""
        fields = line.split(maxsplit=6)
        if len(fields) >= 6:
            # this is probably an airss res file
            seed, pressure, volume, energy, spin, absspin = fields[:6]
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
                seed, float(pressure), float(volume), float(energy), float(spin), float(absspin), spg, int(nap)
            )
        else:
            # there should at least be the first 6 fields if it's an airss res file
            # if it doesn't have them, then just stop looking
            return None

    def _parse_cell(self, line: str) -> ResCELL:
        """Parses the CELL entry."""
        fields = line.split()
        if len(fields) != 7:
            raise ParseError(f"Failed to parse CELL line {line}, expected 7 fields.")
        field_1, a, b, c, alpha, beta, gamma = map(float, fields)
        return ResCELL(field_1, a, b, c, alpha, beta, gamma)

    def _parse_ion(self, line: str) -> Ion:
        """Parses entries in the SFAC block."""
        fields = line.split()
        if len(fields) != 6:
            raise ParseError(f"Failed to parse ion entry {line}, expected 6 fields.")
        specie = fields[0]
        specie_num = int(fields[1])
        x, y, z, occ = map(float, fields[2:])
        return Ion(specie, specie_num, (x, y, z), occ)

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
            raise ParseError("Encountered end of file before END tag at end of SFAC block.")
        return ResSFAC(species, ions)

    def _parse_txt(self) -> Res:
        """Parses the text of the file."""
        _REMS: List[str] = []
        _TITL: Optional[AirssTITL] = None
        _CELL: Optional[ResCELL] = None
        _SFAC: Optional[ResSFAC] = None

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
                elif splits == 1:
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
                    raise Warning(f"Skipping line {line}, tag {first} not recognized.")
        except StopIteration:
            pass
        if _CELL is None or _SFAC is None:
            raise ParseError("Did not encounter CELL or SFAC entry when parsing.")
        return Res(_TITL, _REMS, _CELL, _SFAC)

    @classmethod
    def _parse_str(cls, source: str) -> Res:
        """Parses the res file as a string."""
        self = cls()
        self.source = source
        return self._parse_txt()

    @classmethod
    def _parse_filename(cls, filename: str) -> Res:
        """Parses the res file as a file."""
        self = cls()
        self.filename = filename
        with zopen(filename, "r") as file:
            self.source = file.read()
            return self._parse_txt()


def _lattice(cell: ResCELL) -> Lattice:
    """Produce a pymatgen Lattice from a parsed CELL entry."""
    return Lattice.from_parameters(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)


def _sites(sfactag: ResSFAC, lattice: Lattice) -> List[PeriodicSite]:
    """Produce a list of pymatgen PeriodicSite from a parsed SFAC block and a pymatgen Lattice."""
    return [PeriodicSite(ion.specie, ion.pos, lattice) for ion in sfactag.ions]


def _structure(res: Res) -> Structure:
    """Produce a pymatgen Structure from a parsed Res file."""
    return Structure.from_sites(_sites(res.SFAC, _lattice(res.CELL)))


def _entry(res: Res) -> ComputedStructureEntry:
    """Produce a pymatgen ComputedStructureEntry from a parsed Res file."""
    energy = res.TITL and res.TITL.energy or float("nan")
    return ComputedStructureEntry(_structure(res), energy)


def _cell(lattice: Lattice) -> ResCELL:
    """Produce CELL entry from a pymatgen Lattice."""
    return ResCELL(1.0, lattice.a, lattice.b, lattice.c, lattice.alpha, lattice.beta, lattice.gamma)


def _ions(sites: List[PeriodicSite]) -> List[Ion]:
    """Produce a list of entries for a SFAC block from a list of pymatgen PeriodicSite."""
    ions: List[Ion] = []
    i = 0
    for site in sites:
        for specie, occ in site.species.items():
            i += 1
            x, y, z = map(float, site.frac_coords)
            ions.append(Ion(specie, i, (x, y, z), occ))
    return ions


def _sfac(sites: List[PeriodicSite]) -> ResSFAC:
    """Produce a SFAC block from a list of pymatgen PeriodicSite."""
    ions = _ions(sites)
    species = {ion.specie for ion in ions}
    return ResSFAC(species, ions)


def _structure2res(structure: Structure) -> Res:
    """Produce a res file structure from a pymatgen Structure."""
    return Res(None, [], _cell(structure.lattice), _sfac(list(structure.sites)))


def _entry2res(entry: ComputedStructureEntry) -> Res:
    """Produce a res file structure from a pymatgen ComputedStructureEntry."""
    pres = float(entry.parameters.get("pressure", 0))
    isd = float(entry.parameters.get("isd", 0))
    iasd = float(entry.parameters.get("iasd", 0))
    spg, _ = entry.structure.get_space_group_info()
    rems = [f"{str(k)} : {str(v)}" for k, v in entry.parameters.items() if k not in ["pressure", "isd", "iasd"]]
    return Res(
        AirssTITL(str(hash(entry)), pres, entry.structure.volume, entry.energy, isd, iasd, spg, 1),
        rems,
        _cell(entry.structure.lattice),
        _sfac(list(entry.structure.sites)),
    )


class ResIO:
    """
    Class providing methods for converting a `Structure` or `ComputedStructureEntry`
    to/from a string or file in the res format as used by AIRSS.

    Note: If the TITL entry doesn't exist or is malformed, then the resulting
    `ComputedStructureEntry` will have the energy set to `float('nan')`.
    """

    @classmethod
    def structure_from_txt(cls, txt: str) -> Structure:
        """
        Produces a pymatgen Structure from contents of a res file.
        """
        return _structure(ResParser._parse_str(txt))

    @classmethod
    def structure_from_file(cls, filename: str) -> Structure:
        """
        Produces a pymatgen Structure from a res file.
        """
        return _structure(ResParser._parse_filename(filename))

    @classmethod
    def structure_to_txt(cls, structure: Structure) -> str:
        """
        Produce the contents of a res file from a pymatgen Structure.
        """
        return str(_structure2res(structure))

    @classmethod
    def structure_to_file(cls, structure: Structure, filename: str) -> None:
        """
        Write a pymatgen Structure to a res file.
        """
        with zopen(filename, "w") as file:
            file.write(cls.structure_to_txt(structure))
        return

    @classmethod
    def entry_from_txt(cls, txt: str) -> ComputedStructureEntry:
        """
        Produce a pymatgen ComputedStructureEntry from contents of a res file.
        """
        return _entry(ResParser._parse_str(txt))

    @classmethod
    def entry_from_file(cls, filename: str) -> ComputedStructureEntry:
        """
        Produce a pymatgen ComputedStructureEntry from a res file.
        """
        return _entry(ResParser._parse_filename(filename))

    @classmethod
    def entry_to_txt(cls, entry: ComputedStructureEntry) -> str:
        """
        Produce the contents of a res file from a pymatgen ComputedStructureEntry.
        """
        return str(_entry2res(entry))

    @classmethod
    def entry_to_file(cls, entry: ComputedStructureEntry, filename: str) -> None:
        """
        Write a pymatgen ComputedStructureEntry to a res file.
        """
        with zopen(filename, "w") as file:
            file.write(cls.entry_to_txt(entry))
        return
