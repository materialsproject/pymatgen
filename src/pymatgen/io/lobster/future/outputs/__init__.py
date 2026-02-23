"""Module for reading Lobster output files. For more information on LOBSTER see www.cohp.de.

If you use this module, please cite:

J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier,
"Automated Bonding Analysis with Crystal Orbital Hamilton Populations", ChemPlusChem 2022, e202200123,
DOI: 10.1002/cplu.202200123.
"""

from __future__ import annotations

from pymatgen.io.lobster.future.outputs.bands import BandOverlaps, Fatband, Fatbands
from pymatgen.io.lobster.future.outputs.coxxcar import COBICAR, COBICAR_LCFO, COHPCAR, COHPCAR_LCFO, COOPCAR, COXXCAR
from pymatgen.io.lobster.future.outputs.doscar import DOSCAR, DOSCAR_LCFO
from pymatgen.io.lobster.future.outputs.icoxxlist import (
    ICOBILIST,
    ICOBILIST_LCFO,
    ICOHPLIST,
    ICOHPLIST_LCFO,
    ICOOPLIST,
    ICOXXLIST,
    NcICOBILIST,
)
from pymatgen.io.lobster.future.outputs.lobsterout import LobsterOut
from pymatgen.io.lobster.future.outputs.misc import (
    BWDF,
    BWDFCOHP,
    POLARIZATION,
    LobsterMatrices,
    MadelungEnergies,
    SitePotentials,
    Wavefunction,
)
from pymatgen.io.lobster.future.outputs.populations import CHARGE, CHARGE_LCFO, GROSSPOP, GROSSPOP_LCFO
from pymatgen.util.due import Doi, due  # type: ignore[import]

__author__ = "Tom Demeyere"
__copyright__ = "Copyright 2025, The Materials Project"
__version__ = "0.3"
__maintainer__ = "Tom Demeyere"
__email__ = "tom.demeyere@bam.de"
__date__ = "Sep. 30, 2025"

due.cite(
    Doi("10.1002/cplu.202200123"),
    description="Automated Bonding Analysis with Crystal Orbital Hamilton Populations",
)

__all__ = [
    "BWDF",
    "BWDFCOHP",
    "CHARGE",
    "CHARGE_LCFO",
    "COBICAR",
    "COBICAR_LCFO",
    "COHPCAR",
    "COHPCAR_LCFO",
    "COOPCAR",
    "COXXCAR",
    "DOSCAR",
    "DOSCAR_LCFO",
    "GROSSPOP",
    "GROSSPOP_LCFO",
    "ICOBILIST",
    "ICOBILIST_LCFO",
    "ICOHPLIST",
    "ICOHPLIST_LCFO",
    "ICOOPLIST",
    "ICOXXLIST",
    "POLARIZATION",
    "BandOverlaps",
    "Fatband",
    "Fatbands",
    "LobsterMatrices",
    "LobsterOut",
    "MadelungEnergies",
    "NcICOBILIST",
    "SitePotentials",
    "Wavefunction",
]
