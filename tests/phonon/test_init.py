from __future__ import annotations

import pymatgen.phonon as ph
import pymatgen.phonon.bandstructure as bs
import pymatgen.phonon.gruneisen as gru
from pymatgen.phonon import dos, plotter


def test_convenience_imports():
    assert ph.PhononBandStructure == bs.PhononBandStructure
    assert ph.PhononBandStructureSymmLine == bs.PhononBandStructureSymmLine
    assert ph.PhononDos == dos.PhononDos
    assert ph.CompletePhononDos == dos.CompletePhononDos
    assert ph.GruneisenParameter == gru.GruneisenParameter
    assert ph.GruneisenPhononBandStructure == gru.GruneisenPhononBandStructure
    assert ph.GruneisenPhononBandStructureSymmLine == gru.GruneisenPhononBandStructureSymmLine
    assert ph.PhononDosPlotter == plotter.PhononDosPlotter
    assert ph.PhononBSPlotter == plotter.PhononBSPlotter
    assert ph.GruneisenPlotter == plotter.GruneisenPlotter
    assert ph.GruneisenPhononBSPlotter == plotter.GruneisenPhononBSPlotter
