# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Tests for the Materials Project DFT mixing scheme


"""


__author__ = "Ryan Kingsbury"
__copyright__ = "Copyright 2019-2021, The Materials Project"
__version__ = "0.1"
__email__ = "RKingsbury@lbl.gov"
__date__ = "October 2021"

import unittest
import pytest  # type: ignore

import pandas as pd
import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.entries.mixing_scheme import MaterialsProjectDFTMixingScheme
from pymatgen.entries.compatibility import (
    CompatibilityError,
)
from pymatgen.entries.computed_entries import (
    ComputedEntry,
    ComputedStructureEntry,
)
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.analysis.phase_diagram import PhaseDiagram


class MaterialsProjectDFTMixingSchemeTest(unittest.TestCase):
    def test_no_structure(self):
        # If we try to process a regular ComputedEntry, should get a warning
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedEntry("Sn", 0, parameters={"run_type": "GGA"}),
            ComputedEntry("Br", 0, parameters={"run_type": "GGA"}),
            ComputedEntry("SnBr2", -10, parameters={"run_type": "R2SCAN"}),
            ComputedEntry("SnBr2", -100, parameters={"run_type": "GGA"}),
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                0,
                parameters={"run_type": "GGA"},
                correction=-20,
            ),
        ]

        with pytest.warns(UserWarning, match="not a ComputedStructureEntry"):
            MaterialsProjectDFTMixingScheme(compat_1=None).process_entries(entries)

    def test_empty_entries(self):
        # Test behavior when either gga_entries or scan_entries passed to get_adjustments
        # is empty
        entries = []
        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        compat.process_entries(entries)

    def test_clean(self):
        # make sure the clean=True arg to process_entries works
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
                correction=-20,
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
                correction=-20,
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "GGA"},
                correction=-20,
            ),
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                0,
                parameters={"run_type": "R2SCAN"},
                correction=-20,
            ),
        ]

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        compat.process_entries(entries, clean=False)
        for e in entries:
            assert e.correction == -20

        compat.process_entries(entries, clean=True)
        for e in entries:
            assert e.correction == 0

    def test_no_run_type(self):
        # should raise a ValueError if any of the entries does not have a run_type
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(Structure(lattice, ["Br"], [[0, 0, 0]]), 0, parameters={}),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
        ]

        with pytest.warns(UserWarning, match="missing parameters.run_type"):
            MaterialsProjectDFTMixingScheme(compat_1=None).process_entries(entries)

    def test_no_single_entry(self):
        # Raise CompatibilityError if process_entries is called on a single entry
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "R2SCAN"},
            )
        ]

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        with pytest.warns(UserWarning, match="cannot process single entries"):
            compat.process_entries(entries)

    def test_only_scan_entries(self):
        # If all entries are SCAN, do nothing
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
        ]

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        for e in entries:
            assert (
                compat.get_adjustments(
                    e, EntrySet([]), EntrySet(entries), EntrySet(entries).remove_non_ground_states(), None, [], []
                )
                == []
            )

        compat.process_entries(entries)
        for e in entries:
            assert e.correction == 0

    def test_only_gga_entries(self):
        # If all entries are GGA(+U), do nothing
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                -20,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                -10,
                parameters={"run_type": "GGA"},
            ),
        ]
        compat = MaterialsProjectDFTMixingScheme(compat_1=None)

        for e in entries:
            assert (
                compat.get_adjustments(
                    e,
                    EntrySet(entries),
                    EntrySet([]),
                    EntrySet(entries).remove_non_ground_states(),
                    PhaseDiagram(entries),
                    [False] * len(entries),
                    [False] * len(entries),
                )
                == []
            )
        compat.process_entries(entries)
        for e in entries:
            assert e.correction == 0

    def test_incompatible_run_type(self):
        # If entry.parameters.run_type is not "GGA", "GGA+U", or "R2SCAN", raise
        # a CompatibilityError and ignore that entry
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        lda_entry = ComputedStructureEntry(
            Structure(lattice, ["Sn"], [[0, 0, 0]]),
            0,
            parameters={"run_type": "LDA"},
        )
        gga_entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "GGA"},
            ),
        ]
        scan_entries = []

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        gga_gs = EntrySet(gga_entries)
        gga_gs.remove_non_ground_states()
        with pytest.raises(CompatibilityError, match="Invalid run type LDA"):
            assert (
                compat.get_adjustments(
                    lda_entry,
                    EntrySet(gga_entries),
                    EntrySet(scan_entries),
                    gga_gs,
                    None,
                    [False] * len(gga_entries),
                    [False] * len(gga_entries),
                )
                == []
            )

        entries = compat.process_entries([lda_entry] + gga_entries + scan_entries)
        assert len(entries) == 2

    def test_incomplete_phase_diagram(self):
        # Test behavior when GGA entries don't form a complete phase diagram
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        gga_entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                -10,
                parameters={"run_type": "GGA"},
            ),
        ]
        scan_entries = [
            # SCAN entry with a composition not found in the GGA references
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
        ]

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        with pytest.warns(UserWarning, match="do not form a complete PhaseDiagram!"):
            compat.process_entries(gga_entries + scan_entries)

    def test_majority_scan(self):
        # If there are SCAN entries for all of the GGA structures, SCAN
        # corrections should be zero and the GGA entries should raise
        # CompatibilityError
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        scan_entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
        ]
        gga_entries = [
            # GGA entries with the same structures as SCAN entries
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "GGA"},
            ),
        ]

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        gga_gs = EntrySet(gga_entries)
        gga_gs.remove_non_ground_states()
        for e in scan_entries + gga_entries:
            if e.parameters["run_type"] == "GGA":
                with pytest.raises(CompatibilityError, match="already exists in SCAN"):
                    assert (
                        compat.get_adjustments(
                            e,
                            EntrySet(gga_entries),
                            EntrySet(scan_entries),
                            gga_gs,
                            None,
                            [True] * len(gga_entries),
                            [],
                        )
                        == []
                    )
            elif e.parameters["run_type"] == "R2SCAN":
                assert (
                    compat.get_adjustments(
                        e, EntrySet(gga_entries), EntrySet(scan_entries), gga_gs, None, [True] * len(gga_entries), []
                    )
                    == []
                )

    def test_gga_correction_scan_hull(self):
        # If there are SCAN entries for all of the stable GGA structures, SCAN
        # corrections should be zero, stable GGA entries should raise
        # CompatibilityError, and unstable GGA entries should be corrected
        # to maintain the same energy above the SCAN hull
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        lattice2 = Lattice.from_parameters(a=1, b=1, c=0.8, alpha=120, beta=90, gamma=60)
        scan_entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                -25,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                -25,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                -100,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                -50,
                parameters={"run_type": "R2SCAN"},
            ),
        ]
        gga_entries = [
            # GGA entries with the same structures as SCAN entries
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                -5,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                -10,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                0,
                parameters={"run_type": "GGA"},
            ),
            # unstable GGA entries without corresponding SCAN structures
            # e above hull = 3 eV; final energy should be -25+3 = -22
            ComputedStructureEntry(
                Structure(lattice2, ["Br"], [[0, 0, 0]]),
                -2,
                parameters={"run_type": "GGA"},
            ),
            # e above hull = 6 eV; final energy should be -100+6 = -94 or 31.333 eV/atom
            ComputedStructureEntry(
                Structure(
                    lattice2,
                    ["Sn", "Br", "Br"],
                    [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]],
                ),
                -4,
                parameters={"run_type": "GGA"},
            ),
        ]

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        mixing_state = compat._generate_mixing_scheme_state_data(gga_entries + scan_entries, verbose=True)

        assert (
            compat.get_adjustments(
                gga_entries[-2],
                mixing_state,
            )[0].value
            == pytest.approx(-20)
        )

        assert (
            compat.get_adjustments(
                gga_entries[-1],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                [True] * len(gga_entries),
                [True] * len(gga_entries),
            )[0].value
            == pytest.approx(-90)
        )

        assert (
            compat.get_adjustments(
                scan_entries[2],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                [True] * len(gga_entries),
                [True] * len(gga_entries),
            )
            == []
        )

        with pytest.raises(CompatibilityError, match="already exists in SCAN"):
            assert compat.get_adjustments(
                gga_entries[-3],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                [True] * len(gga_entries),
                [True] * len(gga_entries),
            )[0].value

    def test_missing_comp_gga(self):
        # test the case where the only GGA entry is for a composition not covered
        # by the SCAN calculations
        pass

    def test_missing_comp_scan(self):
        # test the case where the only SCAN entry is for a composition not covered
        # by the GGA calculations (i.e., there is no)
        pass

    def test_no_scan_references(self):
        # test the case where none of the SCAN entries correspond to a GGA
        # reference structure
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        lattice2 = Lattice.from_parameters(a=1, b=1, c=0.8, alpha=120, beta=90, gamma=60)
        gga_entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                -5,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                -10,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(
                    lattice,
                    ["Sn", "Br", "Br", "Br", "Br"],
                    [
                        [0, 0, 0],
                        [0.2, 0.2, 0.2],
                        [0.4, 0.4, 0.4],
                        [0.7, 0.7, 0.7],
                        [1, 1, 1],
                    ],
                ),
                0,
                parameters={"run_type": "GGA"},
            ),
        ]
        scan_entries = [
            # SCAN entries without corresponding SCAN reference structures
            ComputedStructureEntry(
                Structure(lattice2, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            ComputedStructureEntry(
                Structure(
                    lattice2,
                    ["Sn", "Br", "Br"],
                    [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]],
                ),
                -100,
                parameters={"run_type": "R2SCAN"},
            ),
        ]
        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        gga_gs = EntrySet(gga_entries)
        gga_gs.remove_non_ground_states()

        assert (
            compat.get_adjustments(
                gga_entries[0],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                [False] * len(gga_entries),
                [False] * len(gga_entries),
            )
            == []
        )

        with pytest.raises(CompatibilityError, match="there are no SCAN reference"):
            compat.get_adjustments(
                scan_entries[-1],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                [False] * len(gga_entries),
                [False] * len(gga_entries),
            )[0].value == pytest.approx(-20)

            compat.get_adjustments(
                scan_entries[-2],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                [False] * len(gga_entries),
                [False] * len(gga_entries),
            )[0].value == pytest.approx(-20)

    def test_scan_polymorph_pair(self):
        # If we have a pair of SCAN calculations at a single composition, one of
        # which corresponds to the GGA reference structure, we should correct
        # the SCAN entry to GGA and discard the GGA entry
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        lattice2 = Lattice.from_parameters(a=1, b=1, c=0.8, alpha=120, beta=90, gamma=60)
        gga_entries = [
            # GGA reference structures
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                -10,
                parameters={"run_type": "GGA"},
            ),
            # GGA unstable polymorph
            ComputedStructureEntry(
                Structure(
                    lattice2,
                    ["Sn", "Br", "Br"],
                    [[0, 0, 0], [0.7, 0.7, 0.7], [1, 1, 1]],
                ),
                -9,
                parameters={"run_type": "GGA"},
            ),
        ]
        scan_entries = [
            # SCAN reference structure
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                -100,
                parameters={"run_type": "R2SCAN"},
            ),
            # SCAN unstable polymorph structure, e_above_hull = 2 eV
            ComputedStructureEntry(
                Structure(
                    lattice2,
                    ["Sn", "Br", "Br"],
                    [[0, 0, 0], [0.7, 0.7, 0.7], [1, 1, 1]],
                ),
                -98,
                parameters={"run_type": "R2SCAN"},
            ),
        ]

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        gga_gs = EntrySet(gga_entries)
        gga_gs.remove_non_ground_states()
        has_scan_ground_states = []
        has_scan_hull_entries = []
        for e in gga_gs:
            if e.composition.reduced_formula == "SnBr2":
                has_scan_ground_states.append(True)
                has_scan_hull_entries.append(True)
            else:
                has_scan_ground_states.append(False)
                has_scan_hull_entries.append(False)

        assert (
            compat.get_adjustments(
                scan_entries[-1],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                has_scan_ground_states,
                has_scan_hull_entries,
            )[0].value
            == pytest.approx(90)
        )

        assert (
            compat.get_adjustments(
                scan_entries[-2],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                has_scan_ground_states,
                has_scan_hull_entries,
            )[0].value
            == pytest.approx(90)
        )

        assert (
            compat.get_adjustments(
                gga_entries[0],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                has_scan_ground_states,
                has_scan_hull_entries,
            )
            == []
        )

        with pytest.raises(CompatibilityError, match="already exists in SCAN"):
            assert (
                compat.get_adjustments(
                    gga_entries[2],
                    EntrySet(gga_entries),
                    EntrySet(scan_entries),
                    gga_gs,
                    PhaseDiagram(gga_entries),
                    has_scan_ground_states,
                    has_scan_hull_entries,
                )
                == []
            )
            assert (
                compat.get_adjustments(
                    gga_entries[3],
                    EntrySet(gga_entries),
                    EntrySet(scan_entries),
                    gga_gs,
                    PhaseDiagram(gga_entries),
                    has_scan_ground_states,
                    has_scan_hull_entries,
                )
                == []
            )

    def test_wrong_scan_references(self):
        # Raise CompatibilityError if there is no SCAN reference structure
        # for the chosen composition
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        gga_entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                -10,
                parameters={"run_type": "GGA"},
            ),
        ]
        scan_entries = [
            # SCAN reference structure at a different composition
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
            # SCAN entry with a different structure than the GGA reference
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.7, 0.7, 0.7], [1, 1, 1]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
        ]

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        gga_gs = EntrySet(gga_entries)
        gga_gs.remove_non_ground_states()
        has_scan_ground_states = []
        has_scan_hull_entries = []
        for e in gga_gs:
            if e.composition.reduced_formula == "Sn":
                has_scan_ground_states.append(True)
                has_scan_hull_entries.append(True)
            else:
                has_scan_ground_states.append(False)
                has_scan_hull_entries.append(False)

        with pytest.raises(CompatibilityError, match="there is no SCAN reference"):
            compat.get_adjustments(
                scan_entries[-1],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                has_scan_ground_states,
                has_scan_hull_entries,
            )

    def test_none_if_fewer_than_2_scan(self):
        # Raise CompatibilityError if there are fewer than 2 SCAN entries
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        gga_entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                -10,
                parameters={"run_type": "GGA"},
            ),
        ]
        scan_entries = [
            # one SCAN entry with the same structure than the GGA reference
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "R2SCAN"},
            ),
        ]

        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        gga_gs = EntrySet(gga_entries)
        gga_gs.remove_non_ground_states()
        has_scan_ground_states = []
        has_scan_hull_entries = []
        for e in gga_gs:
            if e.composition.reduced_formula == "Sn":
                has_scan_ground_states.append(True)
                has_scan_hull_entries.append(True)
            else:
                has_scan_ground_states.append(False)
                has_scan_hull_entries.append(False)

        with pytest.raises(CompatibilityError, match="fewer than 2 SCAN"):
            compat.get_adjustments(
                scan_entries[0],
                EntrySet(gga_entries),
                EntrySet(scan_entries),
                gga_gs,
                PhaseDiagram(gga_entries),
                has_scan_ground_states,
                has_scan_hull_entries,
            )

    def test_get_adjustments(self):
        """
        Unit tests of the get_adjustments function internal to the mixing scheme
        This function calculates the appropriate energy correction for a single
        entry, given a DataFrame of information about the overall collection
        of entries being mixed.

        The first 3 columns are informational only and not used by get_adjustments,
        so they are populated with dummy values here
        """
        compat = MaterialsProjectDFTMixingScheme(compat_1=None)
        columns = [
            "composition",
            "spacegroup",
            "num_sites",
            "run_type_1",
            "run_type_2",
            "ground_state_energy_1",
            "ground_state_energy_2",
            "is_stable_1",
            "hull_energy_1",
            "hull_energy_2",
        ]
        # Mixing state 1 - we have all run_type_1 ground states in run_type_2
        # The energy hull is built with run_type_2 energies
        # the R2SCAN hull is exactly 10 eV/atom below the GGA in this example
        row_list = [
            [Composition("Ti"), 194, 2, "GGA", "R2SCAN", -7.89, -17.89, True, -7.89, -17.89],
            [Composition("O2"), 12, 8, "GGA", "R2SCAN", -4.95, -14.95, True, -4.95, -14.95],
            [Composition("Ti3O5"), 12, 16, "GGA", "R2SCAN", -9.01, -19.01, True, -9.01, -19.01],
            [Composition("Ti3O4"), 139, 7, "GGA", "R2SCAN", -8.98, -18.98, False, -9, -19],
            [Composition("Ti3O"), 149, 24, "GGA", "R2SCAN", -8.55, -18.55, True, -8.55, -18.55],
            [Composition("TiO"), 189, 6, "GGA", "R2SCAN", -8.99, -18.99, True, -8.99, -18.99],
            [Composition("Ti2O3"), 167, 10, "GGA", "R2SCAN", -9.02, -19.02, True, -9.02, -19.02],
            [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
        ]
        mixing_state = pd.DataFrame(row_list, columns=columns)

        # correction for a run_type_1 entry
        entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
        with pytest.raises(CompatibilityError, match="Discarding GGA entry"):
            adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

        # correction for a run_type_2 entry
        entry = ComputedEntry("TiO2", -15, parameters={"run_type": "R2SCAN"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj == []

        # Mixing state 2 - We have all run_type_1 stable entries (but not all
        # ground states) in run_type_2. The energy hull is built with run_type_2 energies
        # the R2SCAN hull is exactly 10 eV/atom below the GGA in this example
        row_list = [
            [Composition("Ti"), 194, 2, "GGA", "R2SCAN", -7.89, -17.89, True, -7.89, -17.89],
            [Composition("O2"), 12, 8, "GGA", "R2SCAN", -4.95, -14.95, True, -4.95, -14.95],
            [Composition("Ti3O5"), 12, 16, "GGA", "R2SCAN", -9.01, -19.01, True, -9.01, -19.01],
            [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, -19],
            [Composition("Ti3O"), 149, 24, "GGA", "R2SCAN", -8.55, -18.55, True, -8.55, -18.55],
            [Composition("TiO"), 189, 6, "GGA", "R2SCAN", -8.99, -18.99, True, -8.99, -18.99],
            [Composition("Ti2O3"), 167, 10, "GGA", "R2SCAN", -9.02, -19.02, True, -9.02, -19.02],
            [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
        ]
        mixing_state = pd.DataFrame(row_list, columns=columns)

        # correction for a run_type_1 entry
        entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
        with pytest.raises(CompatibilityError, match="Discarding GGA entry"):
            adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

        # correction for a run_type_1 entry where there is no run_type_2 composition
        entry = ComputedEntry("Ti3O4", -15, parameters={"run_type": "GGA"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
        assert np.allclose(adj[0].value, -10 * 7)

        # correction for a run_type_2 entry
        entry = ComputedEntry("TiO2", -15, parameters={"run_type": "R2SCAN"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj == []

        # Mixing state 3 - We have run_type_2 for only a few ground states
        row_list = [
            [Composition("Ti"), 194, 2, "GGA", "R2SCAN", -7.89, -17.89, True, -7.89, -17.89],
            [Composition("O2"), 12, 8, "GGA", "R2SCAN", -4.95, -14.95, True, -4.95, -14.95],
            [Composition("Ti3O5"), 12, 16, "GGA", "R2SCAN", -9.01, -19.01, True, -9.01, -19.01],
            [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, None],
            [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
            [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
            [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
            [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
        ]
        mixing_state = pd.DataFrame(row_list, columns=columns)

        # correction for a run_type_1 entry
        # no correction
        entry = ComputedEntry("Ti3O", -8.55 * 4, parameters={"run_type": "GGA"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj == []

        # correction for a run_type_2 entry that has a run_type_2 ground state
        # this energy should be corrected
        entry = ComputedEntry("Ti3O5", -15 * 8, parameters={"run_type": "R2SCAN"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
        # r2SCAN entry is 4.01 eV/atom above hull. Add to GGA hull energy
        assert np.allclose(adj[0].value, (-9.01 + 4.01 - -15) * 8)

        # correction for a run_type_2 entry that does not have a run_type_2 ground state
        # this entry should be discarded
        entry = ComputedEntry("Ti3O4", -15 * 7, parameters={"run_type": "R2SCAN"})
        with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
            adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

        # Mixing state 4 - We have run_type_2 for only 2 structures, at a single composition
        row_list = [
            [Composition("Ti"), 194, 2, "GGA", None, -7.89, np.nan, True, -7.89, None],
            [Composition("O2"), 12, 8, "GGA", None, -4.95, np.nan, True, -4.95, None],
            [Composition("Ti3O5"), 12, 16, "GGA", None, -9.01, np.nan, True, -9.01, None],
            [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, None],
            [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
            [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
            [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
            [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
            [Composition("TiO2"), 189, 12, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
        ]
        mixing_state = pd.DataFrame(row_list, columns=columns)

        # correction for a run_type_1 entry at this composition
        # no correction
        entry = ComputedEntry("Ti2O4", -8.97 * 6, parameters={"run_type": "GGA"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj == []

        # correction for a run_type_2 corresponding to the reference state
        # should be discarded
        entry = ComputedEntry("Ti2O4", -18.97 * 6, parameters={"run_type": "R2SCAN"})
        with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
            adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

        # correction for a run_type_2 entry that does not have a run_type_2 ground state
        # this entry should be corrected to the GGA hull
        entry = ComputedEntry("Ti4O8", -15.97 * 12, parameters={"run_type": "R2SCAN"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
        # r2SCAN entry is 3 eV/atom above hull. Add to GGA hull energy
        assert np.allclose(adj[0].value, (-8.97 + 3 - entry.energy_per_atom) * 12)

        # Mixing state 5 - We have only 1 run_type_2 material and it is a reference state
        # This is an edge case that should never happen, because _generate_mixing_scheme_state_data
        # Will discard run_type_2 entries if there are fewer than 2 of them, so they'll never be in the DataFrame
        row_list = [
            [Composition("Ti"), 194, 2, "GGA", None, -7.89, np.nan, True, -7.89, None],
            [Composition("O2"), 12, 8, "GGA", None, -4.95, np.nan, True, -4.95, None],
            [Composition("Ti3O5"), 12, 16, "GGA", None, -9.01, np.nan, True, -9.01, None],
            [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, None],
            [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
            [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
            [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
            [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
        ]
        mixing_state = pd.DataFrame(row_list, columns=columns)

        # correction for a run_type_1 entry
        entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj == []

        # correction for the single run_type_2 entry that is the reference state
        entry = ComputedEntry("TiO2", -18.97 * 3, parameters={"run_type": "R2SCAN"})
        with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
            adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

        # correction for a run_type_2 entry at the same composition that was not
        # included in the original list of processed entries (Another edge case that should
        # never happen)
        entry = ComputedEntry("TiO2", -17.97 * 3, parameters={"run_type": "R2SCAN"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
        # r2SCAN entry is 1 eV/atom above hull. Add to GGA hull energy
        assert np.allclose(adj[0].value, (-8.97 + 1 - entry.energy_per_atom) * 3)

        # Mixing state 6 - We have 0 run_type_2 materials
        row_list = [
            [Composition("Ti"), 194, 2, "GGA", None, -7.89, np.nan, True, -7.89, None],
            [Composition("O2"), 12, 8, "GGA", None, -4.95, np.nan, True, -4.95, None],
            [Composition("Ti3O5"), 12, 16, "GGA", None, -9.01, np.nan, True, -9.01, None],
            [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, None],
            [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
            [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
            [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
            [Composition("TiO2"), 141, 6, "GGA", None, -8.97, np.nan, True, -8.97, None],
        ]
        mixing_state = pd.DataFrame(row_list, columns=columns)

        # correction for a run_type_1 entry
        entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj == []

        # correction for a run_type_2 entry
        entry = ComputedEntry("TiO2", -15, parameters={"run_type": "R2SCAN"})
        with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
            adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

        # Mixing state 7 - We have run_type_2_materials for everything, including
        # a composition that is not present in run_type_1
        row_list = [
            [Composition("Ti"), 194, 2, "GGA", "R2SCAN", -7.89, -17.89, True, -7.89, -17.89],
            [Composition("O2"), 12, 8, "GGA", "R2SCAN", -4.95, -14.95, True, -4.95, -14.95],
            [Composition("Ti3O5"), 12, 16, "GGA", "R2SCAN", -9.01, -19.01, True, -9.01, -19.01],
            [Composition("Ti3O4"), 139, 7, None, "R2SCAN", np.nan, -18.98, False, -9, -19],
            [Composition("Ti3O"), 149, 24, "GGA", "R2SCAN", -8.55, -18.55, True, -8.55, -18.55],
            [Composition("TiO"), 189, 6, "GGA", "R2SCAN", -8.99, -18.99, True, -8.99, -18.99],
            [Composition("Ti2O3"), 167, 10, "GGA", "R2SCAN", -9.02, -19.02, True, -9.02, -19.02],
            [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
        ]
        mixing_state = pd.DataFrame(row_list, columns=columns)

        # correction for a run_type_1 entry
        entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
        with pytest.raises(CompatibilityError, match="Discarding GGA entry"):
            adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

        # correction for a run_type_2 entry with no corresponding run_type_1
        entry = ComputedEntry("Ti3O4", -8.97 * 3, parameters={"run_type": "R2SCAN"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj == []

        # correction for a run_type_2 entry
        entry = ComputedEntry("TiO2", -15, parameters={"run_type": "R2SCAN"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj == []

        # Mixing state 8 - We have run_type_2 for only 2 structures, at a single composition
        # that is not on the run_type_1 hull
        row_list = [
            [Composition("Ti"), 194, 2, "GGA", None, -7.89, np.nan, True, -7.89, None],
            [Composition("O2"), 12, 8, "GGA", None, -4.95, np.nan, True, -4.95, None],
            [Composition("Ti3O5"), 12, 16, "GGA", None, -9.01, np.nan, True, -9.01, None],
            [Composition("Ti3O4"), 139, 7, "GGA", "R2SCAN", -8.98, -18.98, False, -9, None],
            [Composition("Ti3O4"), 149, 14, None, "R2SCAN", -8.98, -18.98, False, -9, None],
            [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
            [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
            [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
        ]
        mixing_state = pd.DataFrame(row_list, columns=columns)

        # correction for a run_type_1 entry at this composition
        # no correction
        entry = ComputedEntry("Ti3O4", -8.98 * 7, parameters={"run_type": "GGA"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj == []

        # correction for a run_type_2 corresponding to the reference state
        # should be discarded
        entry = ComputedEntry("Ti3O4", -18.98 * 7, parameters={"run_type": "R2SCAN"})
        with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
            adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

        # correction for a run_type_2 entry that does not have a run_type_2 ground state
        # this entry should be corrected to the GGA hull
        entry = ComputedEntry("Ti6O8", -16.98 * 14, parameters={"run_type": "R2SCAN"})
        adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
        assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
        # r2SCAN entry is 2 eV/atom above hull. Add to GGA hull energy
        assert np.allclose(adj[0].value, (-9 + 2 - entry.energy_per_atom) * 14)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
