# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Tests for the Materials Project DFT mixing scheme

In general there are 3 types of tests

1. Functional tests of the behavior of the process_entries method, which is the
   primary user-facing interface
2. Unit tests that get_adjustments behaves correctly, when provided with a pre-generated
   mixing_scheme_state_data DataFrame
3. Unit tests of _generate_mixing_scheme_state_data, to verify that it generates
   the correct DataFrame when provided with a particular list of entries

The tests are structured around different "mixing states" that are intended to capture
all (or most) likely combinations of run_type_1 and run_type_2 entries. - e.g. all
ground states present for both run_types, only one run_type, etc. Mixing states
are defined as class attributes here as both 1) sets of entries and 2) the
corresponding pandas DataFrame that represents the mixing state. 
"""


__author__ = "Ryan Kingsbury"
__copyright__ = "Copyright 2019-2021, The Materials Project"
__version__ = "0.1"
__email__ = "RKingsbury@lbl.gov"
__date__ = "October 2021"

import pytest
import warnings

import pandas as pd
import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.entries.mixing_scheme import MaterialsProjectDFTMixingScheme
from pymatgen.entries.compatibility import (
    Compatibility,
    CompatibilityError,
)
from pymatgen.entries.computed_entries import (
    ComputedEntry,
    ComputedStructureEntry,
    ConstantEnergyAdjustment,
)
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.analysis.phase_diagram import PhaseDiagram


"""
Define utility classes to make the tests easier to read
"""


class DummyCompatibility(Compatibility):
    """
    Dummy class to test compat1 and compat2 kwargs
    """

    def get_adjustments(self, mixing_scheme_no_compat, entry):
        return [ConstantEnergyAdjustment(-10, name="Dummy adjustment")]


class MixingState:
    """
    Lightweight container class for a mixing state, used by the tests below
    """

    def __init__(self, gga_entries, scan_entries, dataframe):
        """
        Args:
            gga_entries: list of run_type_1 (usually GGA or GGA+U) ComputedEntry
            scan_entries: list of run_type_2 (usually r2SCAN) ComputedEntry
            DataFrame: pandas DataFrame representing the mixing state, of the
                format returned by _generate_mixing_scheme_state_data
        """
        self.gga_entries = gga_entries
        self.scan_entries = scan_entries
        self.state_data = dataframe

    @property
    def all_entries(self):
        return self.gga_entries + self.scan_entries


@pytest.fixture
def mixing_scheme_no_compat():
    """
    Return an instance of MaterialsProjectDFTMixingScheme with no additional
    compatibility schemes (e.g., compat_1=None). Used by most of the tests where
    we are manually supplying energies.
    """
    return MaterialsProjectDFTMixingScheme(compat_1=None)


"""
Define mixing states to test

Note that in the DataFrame, the first 3 columns are informational only and not
used by get_adjustments, so they are populated with dummy values here
"""
columns = [
            "composition",
            "spacegroup",
            "num_sites",
            "entry_id_1",
            "entry_id_2",
            "run_type_1",
            "run_type_2",
            "ground_state_energy_1",
            "ground_state_energy_2",
            "is_stable_1",
            "hull_energy_1",
            "hull_energy_2",
        ]


@pytest.fixture
def ms_complete():
    """
    Mixing state where we have SCAN for all GGA
    SCAN energies are 1 eV/atom below the GGA ones
    """
    # lattices. In general, ground states are all assigned lattice1
    # unstable polymorphs are assigned lattice2 or lattice 3
    lattice1 = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
    lattice2 = Lattice.from_parameters(a=1, b=1, c=1, alpha=120, beta=120, gamma=60)
    lattice3 = Lattice.from_parameters(a=1, b=1, c=1, alpha=120, beta=120, gamma=90)
    lattice_br_gga = Lattice.from_dict({'@module': 'pymatgen.core.lattice',
    '@class': 'Lattice',
    'matrix': [[2.129324, -4.226095, 0.0],
    [2.129324, 4.226095, 0.0],
    [0.0, 0.0, 8.743796]]})
    lattice_br_r2scan = Lattice.from_dict({'@module': 'pymatgen.core.lattice',
    '@class': 'Lattice',
    'matrix': [[0.0, -4.25520892, -0.0],
    [-3.56974866, 2.12760446, 0.0],
    [0.0, 0.0, -8.74536848]]})


    gga_entries = [
        ComputedStructureEntry(
            Structure(lattice1, ["Sn"], [[0, 0, 0]]),
            0,
            parameters={"run_type": "GGA"},
            entry_id = "gga-1"
        ),
        ComputedStructureEntry(
            Structure(lattice1, ["Br"], [[0, 0, 0]]),
            1,
            parameters={"run_type": "GGA"},
            entry_id = "gga-2"
        ),
        ComputedStructureEntry(
            Structure(lattice_br_gga, ["Br","Br","Br","Br"], [[0.642473, 0.642473, 0.117751],
        [0.357527, 0.357527, 0.882249],
        [0.857527, 0.857527, 0.617751],
        [0.142473, 0.142473, 0.382249]]),
            0,
            parameters={"run_type": "GGA"},
            entry_id = "gga-3"
        ),
        ComputedStructureEntry(
            Structure(lattice1, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -18,
            parameters={"run_type": "GGA"},
            entry_id = "gga-4"
        ),
        ComputedStructureEntry(
            Structure(lattice2, ["Sn", "Sn", "Sn", "Sn", "Br", "Br", "Br", "Br", "Br", "Br", "Br", "Br"], [[0.25    , 0.393393, 0.663233],
        [0.75    , 0.606607, 0.336767],
        [0.25    , 0.893393, 0.836767],
        [0.75    , 0.106607, 0.163233],
        [0.25    , 0.662728, 0.548755],
        [0.75    , 0.337272, 0.451245],
        [0.25    , 0.162728, 0.951245],
        [0.75    , 0.837272, 0.048755],
        [0.25    , 0.992552, 0.311846],
        [0.75    , 0.007448, 0.688154],
        [0.25    , 0.492552, 0.188154],
        [0.75    , 0.507448, 0.811846]]),
            -60,
            parameters={"run_type": "GGA"},
            entry_id = "gga-5"
        ),
        ComputedStructureEntry(
            Structure(lattice3, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -12,
            parameters={"run_type": "GGA"},
            entry_id = "gga-6"
        ),
        ComputedStructureEntry(
            Structure(
                lattice1,
                ["Sn", "Br", "Br", "Br", "Br"],
                [
                    [0, 0, 0],
                    [0.2, 0.2, 0.2],
                    [0.4, 0.4, 0.4],
                    [0.7, 0.7, 0.7],
                    [1, 1, 1],
                ],
            ),
            -15,
            parameters={"run_type": "GGA"},
            entry_id = "gga-7"
        ),
    ]
    scan_entries = [
        ComputedStructureEntry(
            Structure(lattice1, ["Sn"], [[0, 0, 0]]),
            -1,
            parameters={"run_type": "R2SCAN"},
            entry_id = "r2scan-1"
        ),
        ComputedStructureEntry(
            Structure(lattice1, ["Br"], [[0, 0, 0]]),
            -1,
            parameters={"run_type": "R2SCAN"},
            entry_id = "r2scan-2"
        ),
        ComputedStructureEntry(
            Structure(lattice_br_r2scan, ["Br","Br","Br","Br"], [[ 0.85985939,  0.        ,  0.38410868],
        [ 0.14014061, -0.        ,  0.61589132],
        [ 0.64014061,  0.        ,  0.88410868],
        [ 0.35985939, -0.        ,  0.11589132]]),
            0,
            parameters={"run_type": "R2SCAN"},
            entry_id = "r2scan-3"
        ),
        ComputedStructureEntry(
            Structure(lattice1, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -21,
            parameters={"run_type": "R2SCAN"},
            entry_id = "r2scan-4"
        ),
        ComputedStructureEntry(
            Structure(lattice2, ["Sn", "Sn", "Sn", "Sn", "Br", "Br", "Br", "Br", "Br", "Br", "Br", "Br"], [[0.25    , 0.393393, 0.663233],
        [0.75    , 0.606607, 0.336767],
        [0.25    , 0.893393, 0.836767],
        [0.75    , 0.106607, 0.163233],
        [0.25    , 0.662728, 0.548755],
        [0.75    , 0.337272, 0.451245],
        [0.25    , 0.162728, 0.951245],
        [0.75    , 0.837272, 0.048755],
        [0.25    , 0.992552, 0.311846],
        [0.75    , 0.007448, 0.688154],
        [0.25    , 0.492552, 0.188154],
        [0.75    , 0.507448, 0.811846]]),
            -96,
            parameters={"run_type": "R2SCAN"},
            entry_id = "r2scan-5"
        ),
        ComputedStructureEntry(
            Structure(lattice3, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -18,
            parameters={"run_type": "R2SCAN"},
            entry_id = "r2scan-6"
        ),
        ComputedStructureEntry(
            Structure(
                lattice1,
                ["Sn", "Br", "Br", "Br", "Br"],
                [
                    [0, 0, 0],
                    [0.2, 0.2, 0.2],
                    [0.4, 0.4, 0.4],
                    [0.7, 0.7, 0.7],
                    [1, 1, 1],
                ],
            ),
            -30,
            parameters={"run_type": "R2SCAN"},
            entry_id = "r2scan-7"
        ),
    ]
    # the fmt command tells the black autoformatter not to mess with this block of code
    # it's easier to edit when all the commas are lined up.
    # fmt: off
    row_list = [
        ["Sn",   191,  1, "gga-1", "r2scan-1", "GGA", "R2SCAN",  0, -1, True,     0, -1],
        ["Br",    65,  1, "gga-2", "r2scan-2", "GGA", "R2SCAN",  1, -1, False,    0, -1],
        ["Br",    65,  1, "gga-3", "r2scan-3", "GGA", "R2SCAN",  0,  0, True,     0, -1],
        ["SnBr2", 65, 12, "gga-4", "r2scan-4", "GGA", "R2SCAN", -6, -7, True,    -6, -8],
        ["SnBr2", 65,  3, "gga-5", "r2scan-5", "GGA", "R2SCAN", -5, -8, False,   -6, -8],
        ["SnBr2", 65,  3, "gga-6", "r2scan-6", "GGA", "R2SCAN", -4, -6, False,   -6, -8],
        ["SnBr4",  8,  5, "gga-7", "r2scan-7", "GGA", "R2SCAN", -3, -6, False, -3.6, -6],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_scan_only():
    """
    Mixing state with only SCAN entries
    """
    # only SCAN, no GGA
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
            -18,
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
            -25,
            parameters={"run_type": "R2SCAN"},
        ),
    ]

    row_list = [
        ["Sn", 194, 1, None, "R2SCAN", np.nan, 0, True, np.nan, 0],
        ["Br", 12, 1, None, "R2SCAN", np.nan, 0, True, np.nan, 0],
        ["SnBr2", 12, 3, None, "R2SCAN", np.nan, -6, True, np.nan, -6],
        ["SnBr4", 139, 5, None, "R2SCAN", np.nan, -5, True, np.nan, -5],
    ]
    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState([], entries, mixing_state)


@pytest.fixture
def ms_gga_only():
    """
    Mixing state with only GGA entries
    """
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
            -18,
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
            -25,
            parameters={"run_type": "GGA"},
        ),
    ]
    row_list = [
        ["Sn", 194, 1, "GGA", None, 0, np.nan, True, 0, np.nan],
        ["Br", 12, 1, "GGA", None, 0, np.nan, True, 0, np.nan],
        ["SnBr2", 12, 3, "GGA", None, -6, np.nan, True, -6, np.nan],
        ["SnBr4", 139, 5, "GGA", None, -5, np.nan, True, -5, np.nan],
    ]
    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState(entries, [], mixing_state)


@pytest.fixture
def ms_invalid_run_type():
    """
    Mixing state with all GGA entries except one with an invalid run_type
    """
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
            Structure(lattice, ["Br"], [[0, 0, 0]]),
            0,
            parameters={"run_type": "LDA"},
        ),
        ComputedStructureEntry(
            Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -18,
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
            -25,
            parameters={"run_type": "GGA"},
        ),
    ]
    row_list = [
        ["Sn", 194, 1, "GGA", None, 0, np.nan, True, 0, np.nan],
        ["Br", 12, 1, "GGA", None, 0, np.nan, True, 0, np.nan],
        ["SnBr2", 12, 3, "GGA", None, -6, np.nan, True, -6, np.nan],
        ["SnBr4", 139, 5, "GGA", None, -5, np.nan, True, -5, np.nan],
    ]
    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState(entries, [], mixing_state)


@pytest.fixture
def ms_incomplete_gga_all_scan():
    """
    Mixing state with an incomplete GGA phase diagram
    """
    lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
    gga_entries = [
        ComputedStructureEntry(
            Structure(lattice, ["Sn"], [[0, 0, 0]]),
            0,
            parameters={"run_type": "GGA"},
        ),
        # ComputedStructureEntry(
        #     Structure(lattice, ["Br"], [[0, 0, 0]]),
        #     0,
        #     parameters={"run_type": "GGA"},
        # ),
        ComputedStructureEntry(
            Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -18,
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
            -25,
            parameters={"run_type": "GGA"},
        ),
    ]
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
            -18,
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
            -25,
            parameters={"run_type": "R2SCAN"},
        ),
    ]
    row_list = [
        ["Sn", 194, 1, "GGA", None, 0, np.nan, False, np.nan, np.nan],
        ["Br", 12, 1, None, "R2SCAN", np.nan, 0, False, np.nan, np.nan],
        ["SnBr2", 12, 3, "GGA", "R2SCAN", -6, -6, False, np.nan, np.nan],
        ["SnBr4", 139, 5, "GGA", "R2SCAN", -5, -5, False, np.nan, np.nan],
    ]
    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState(gga_entries, scan_entries, mixing_state)


class TestMaterialsProjectDFTMixingSchemeArgs:
    """
    Test behavior of args / kwargs passed to the mixing scheme.
    In general, these tests are for the interface rather than the actual
    mixing.
    """

    def test_no_structure(self, mixing_scheme_no_compat):
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
            mixing_scheme_no_compat.process_entries(entries)

    def test_empty_entries(self, mixing_scheme_no_compat):
        # Test behavior when either gga_entries or scan_entries passed to get_adjustments
        # is empty
        entries = []
        mixing_scheme_no_compat.process_entries(entries)

    def test_clean(self, mixing_scheme_no_compat):
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

        mixing_scheme_no_compat.process_entries(entries, clean=False)
        for e in entries:
            assert e.correction == -20

        mixing_scheme_no_compat.process_entries(entries, clean=True)
        for e in entries:
            assert e.correction == 0

    def test_no_run_type(self, mixing_scheme_no_compat):
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
            mixing_scheme_no_compat.process_entries(entries)

    def test_no_single_entry(self, mixing_scheme_no_compat):
        # Raise CompatibilityError if process_entries is called on a single entry
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "R2SCAN"},
            )
        ]

        with pytest.warns(UserWarning, match="cannot process single entries"):
            mixing_scheme_no_compat.process_entries(entries)

    @pytest.mark.skip(reason="Not implemented yet")
    def test_fuzzy_diatomic_matching(self, mixing_scheme_no_compat):
        """
        Test fuzzy diatomic matching
        """
        # TODO - write this test
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_compat_args(self, mixing_scheme_no_compat):
        """
        Test the behavior of compat1 and compat2 kwargs
        """
        # TODO - write this test
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_alternate_structure_matcher(self, mixing_scheme_no_compat):
        """
        Test alternate structure matcher kwargs
        """
        # TODO - write this test
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_chemsys_mismatch(self, mixing_scheme_no_compat):
        """
        Test what happens if the entries aren't in the same chemsys
        """
        # TODO - write this test
        pass


class TestTestMaterialsProjectDFTMixingSchemeStates:
    """
    Test the behavior of the mixing scheme under different mixing states (i.e., different
    combinations of GGA and SCAN entries)
    """

    def test_state_complete_entries(self, mixing_scheme_no_compat, ms_complete):
        """
        Mixing state in which every material is present in both GGA and SCAN
        """
        # If there are SCAN entries for all of the GGA entries, SCAN
        # corrections should be zero and the GGA entries should raise
        # CompatibilityError

        state_data = mixing_scheme_no_compat._generate_mixing_scheme_state_data(ms_complete.all_entries)
        assert isinstance(
            state_data, pd.DataFrame
        ), "_generate_mixing_scheme_state_data failed to generate a DataFrame."
        assert all(state_data["run_type_1"] == "GGA")
        assert all(state_data["run_type_2"] == "R2SCAN")
        assert sum(state_data["is_stable_1"]) == 3
        assert all(state_data["ground_state_energy_1"].notna())
        assert all(state_data["ground_state_energy_2"].notna())
        assert all(state_data["hull_energy_1"].notna())
        assert all(state_data["hull_energy_2"].notna())

        for e in ms_complete.scan_entries:
            assert mixing_scheme_no_compat.get_adjustments(e, ms_complete.state_data) == []

        for e in ms_complete.gga_entries:
            with pytest.raises(CompatibilityError, match="already exists in R2SCAN"):
                mixing_scheme_no_compat.get_adjustments(e, ms_complete.state_data)

        # process_entries should discard all GGA entries and return all R2SCAN
        # with pytest.warns(UserWarning, match="do not form a complete PhaseDiagram"):
        entries = mixing_scheme_no_compat.process_entries(ms_complete.all_entries)
        assert len(entries) == 7

    @pytest.mark.skip(reason="Needs revision")
    def test_state_gga_only(self, mixing_scheme_no_compat, ms_gga_only):
        """
        Mixing state in which we only have GGA entries, forming a complete PhaseDiagram
        """
        # If all entries are GGA(+U), do nothing

        state_data = mixing_scheme_no_compat._generate_mixing_scheme_state_data(ms_gga_only.all_entries)
        assert isinstance(
            state_data, pd.DataFrame
        ), "_generate_mixing_scheme_state_data failed to generate a DataFrame."
        assert all(state_data["run_type_1"] == "GGA")
        assert all(state_data["run_type_2"].isna())
        assert len(state_data["is_stable_1"]) == len(ms_gga_only.all_entries)
        assert all(state_data["ground_state_energy_1"].notna())
        assert all(state_data["ground_state_energy_2"].isna())
        assert all(state_data["hull_energy_1"].notna())
        assert all(state_data["hull_energy_2"].isna())

        for e in ms_gga_only.all_entries:
            assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_only.state_data) == []
        mixing_scheme_no_compat.process_entries(ms_gga_only.all_entries)
        for e in ms_gga_only.all_entries:
            assert e.correction == 0

    @pytest.mark.skip(reason="Needs revision")
    def test_state_scan_only(self, mixing_scheme_no_compat, ms_scan_only):
        """
        Mixing state in which we only have SCAN entries, forming a complete PhaseDiagram

        In this case, the mixing scheme should not do anything
        """
        state_data = mixing_scheme_no_compat._generate_mixing_scheme_state_data(ms_scan_only.all_entries, verbose=True)
        assert isinstance(
            state_data, pd.DataFrame
        ), "_generate_mixing_scheme_state_data failed to generate a DataFrame."
        assert all(state_data["run_type_2"] == "R2SCAN")
        assert all(state_data["run_type_1"].isna())
        assert all(~state_data["is_stable_1"])
        assert all(state_data["ground_state_energy_1"].isna())
        assert all(state_data["ground_state_energy_2"].notna())
        assert all(state_data["hull_energy_1"].isna())
        assert all(state_data["hull_energy_2"].notna())

        for e in ms_scan_only.all_entries:
            assert mixing_scheme_no_compat.get_adjustments(e, ms_scan_only.state_data) == []

        mixing_scheme_no_compat.process_entries(ms_scan_only.all_entries)
        for e in ms_scan_only.all_entries:
            assert e.correction == 0

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_gga_1_scan(self, mixing_scheme_no_compat):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and 1 SCAN entry
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_gga_2_scan_same(self, mixing_scheme_no_compat):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and 2 SCAN entries
        at a single composition
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_gga_2_scan_diff(self, mixing_scheme_no_compat):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and 2 SCAN entries
        at different compositions
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_scan_1_gga(self, mixing_scheme_no_compat):
        """
        Mixing state in which we have a complete SCAN PhaseDiagram and 1 GGA entry
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_scan_2_gga_same(self, mixing_scheme_no_compat):
        """
        Mixing state in which we have a complete SCAN PhaseDiagram and 2 GGA entries
        at a single composition
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_incomplete_gga_2_scan_same(self, mixing_scheme_no_compat):
        """
        Mixing state in which we have an incomplete GGA PhaseDiagram and 2 SCAN
        entries at a single composition
        """
        pass

    @pytest.mark.skip(reason="Needs revision")
    def test_state_incomplete_gga_all_scan(self, mixing_scheme_no_compat, ms_incomplete_gga_all_scan):
        """
        Mixing state in which we have an incomplete GGA PhaseDiagram and all entries
        present in SCAN
        """
        # Test behavior when GGA entries don't form a complete phase diagram
        state_data = mixing_scheme_no_compat._generate_mixing_scheme_state_data(ms_incomplete_gga_all_scan.all_entries)
        assert isinstance(
            state_data, pd.DataFrame
        ), "_generate_mixing_scheme_state_data failed to generate a DataFrame."
        assert sum(state_data["run_type_1"] == "GGA") == len(state_data) - 1
        assert sum(state_data["run_type_2"] == "R2SCAN") == 4
        assert all(~state_data["is_stable_1"])
        assert sum(state_data["ground_state_energy_1"].notna()) == len(state_data) - 1
        assert sum(state_data["ground_state_energy_2"].notna()) == len(state_data)
        assert all(state_data["hull_energy_1"].isna())
        assert all(state_data["hull_energy_2"].notna())

        for e in ms_incomplete_gga_all_scan.all_entries:
            with pytest.raises(CompatibilityError, match="Insufficient combination of entries"):
                mixing_scheme_no_compat.get_adjustments(e, ms_incomplete_gga_all_scan.state_data)

        # process_entries should discard all GGA entries and return all SCAN entries
        with pytest.warns(UserWarning, match="do not form a complete PhaseDiagram"):
            entries = mixing_scheme_no_compat.process_entries(ms_incomplete_gga_all_scan.all_entries)
        assert len(entries) == 4

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_all_gga_scan_gs(self, mixing_scheme_no_compat):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and all GGA
        ground states present as SCAN entries
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_foreign_scan_comp(self, mixing_scheme_no_compat):
        """
        Mixing state in which we try to process a SCAN entry at a composition
        that is not in the GGA PhaseDiagram
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_foreign_gga_comp(self, mixing_scheme_no_compat):
        """
        Mixing state in which we try to process a GGA entry at a composition
        that is not in the SCAN PhaseDiagram
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_foreign_gga_below_hull(self, mixing_scheme_no_compat):
        """
        Mixing state in which we try to process a GGA entry that has an energy
        below the GGA hull computed by _generate_mixing_scheme_state_data
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_state_foreign_scan_below_hull(self, mixing_scheme_no_compat):
        """
        Mixing state in which we try to process a SCAN entry that has an energy
        below the SCAN hull computed by _generate_mixing_scheme_state_data
        """
        pass
    
    @pytest.mark.skip(reason="Needs revision")
    def test_incompatible_run_type(self, mixing_scheme_no_compat, ms_invalid_run_type):
        # If entry.parameters.run_type is not "GGA", "GGA+U", or "R2SCAN", raise
        # a CompatibilityError and ignore that entry

        state_data = mixing_scheme_no_compat._generate_mixing_scheme_state_data(ms_invalid_run_type.all_entries)
        assert isinstance(
            state_data, pd.DataFrame
        ), "_generate_mixing_scheme_state_data failed to generate a DataFrame."
        assert all(state_data["run_type_1"] == "GGA")
        assert all(state_data["run_type_2"].isna())
        assert len(state_data["is_stable_1"]) == len(ms_invalid_run_type.all_entries) - 1
        assert all(state_data["ground_state_energy_1"].notna())
        assert all(state_data["ground_state_energy_2"].isna())
        assert all(state_data["hull_energy_1"].notna())
        assert all(state_data["hull_energy_2"].isna())

        with pytest.raises(CompatibilityError, match="Invalid run type LDA"):
            assert (
                mixing_scheme_no_compat.get_adjustments(
                    [e for e in ms_invalid_run_type.all_entries if e.parameters["run_type"] == "LDA"][0],
                    ms_invalid_run_type.state_data,
                )
                == []
            )

        # process_entries should discard the invalid entry
        entries = mixing_scheme_no_compat.process_entries(ms_invalid_run_type.all_entries)
        assert len(entries) == len(ms_invalid_run_type.all_entries) - 1

    # def test_gga_correction_scan_hull(self, mixing_scheme_no_compat):
    #     # If there are SCAN entries for all of the stable GGA structures, SCAN
    #     # corrections should be zero, stable GGA entries should raise
    #     # CompatibilityError, and unstable GGA entries should be corrected
    #     # to maintain the same energy above the SCAN hull
    #     lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
    #     lattice2 = Lattice.from_parameters(a=1, b=1, c=0.8, alpha=120, beta=90, gamma=60)
    #     scan_entries = [
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn"], [[0, 0, 0]]),
    #             -25,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Br"], [[0, 0, 0]]),
    #             -25,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
    #             -100,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(
    #                 lattice,
    #                 ["Sn", "Br", "Br", "Br", "Br"],
    #                 [
    #                     [0, 0, 0],
    #                     [0.2, 0.2, 0.2],
    #                     [0.4, 0.4, 0.4],
    #                     [0.7, 0.7, 0.7],
    #                     [1, 1, 1],
    #                 ],
    #             ),
    #             -50,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #     ]
    #     gga_entries = [
    #         # GGA entries with the same structures as SCAN entries
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Br"], [[0, 0, 0]]),
    #             -5,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
    #             -10,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(
    #                 lattice,
    #                 ["Sn", "Br", "Br", "Br", "Br"],
    #                 [
    #                     [0, 0, 0],
    #                     [0.2, 0.2, 0.2],
    #                     [0.4, 0.4, 0.4],
    #                     [0.7, 0.7, 0.7],
    #                     [1, 1, 1],
    #                 ],
    #             ),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         # unstable GGA entries without corresponding SCAN structures
    #         # e above hull = 3 eV; final energy should be -25+3 = -22
    #         ComputedStructureEntry(
    #             Structure(lattice2, ["Br"], [[0, 0, 0]]),
    #             -2,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         # e above hull = 6 eV; final energy should be -100+6 = -94 or 31.333 eV/atom
    #         ComputedStructureEntry(
    #             Structure(
    #                 lattice2,
    #                 ["Sn", "Br", "Br"],
    #                 [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]],
    #             ),
    #             -4,
    #             parameters={"run_type": "GGA"},
    #         ),
    #     ]

    #     compat = MaterialsProjectDFTMixingScheme(compat_1=None)
    #     mixing_state = compat._generate_mixing_scheme_state_data(gga_entries + scan_entries, verbose=True)

    #     assert (
    #         compat.get_adjustments(
    #             gga_entries[-2],
    #             mixing_state,
    #         )[0].value
    #         == pytest.approx(-20)
    #     )

    #     assert (
    #         compat.get_adjustments(
    #             gga_entries[-1],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             [True] * len(gga_entries),
    #             [True] * len(gga_entries),
    #         )[0].value
    #         == pytest.approx(-90)
    #     )

    #     assert (
    #         compat.get_adjustments(
    #             scan_entries[2],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             [True] * len(gga_entries),
    #             [True] * len(gga_entries),
    #         )
    #         == []
    #     )

    #     with pytest.raises(CompatibilityError, match="already exists in SCAN"):
    #         assert compat.get_adjustments(
    #             gga_entries[-3],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             [True] * len(gga_entries),
    #             [True] * len(gga_entries),
    #         )[0].value

    # def test_missing_comp_gga(self, mixing_scheme_no_compat):
    #     # test the case where the only GGA entry is for a composition not covered
    #     # by the SCAN calculations
    #     # TODO - write this test
    #     pass

    # def test_missing_comp_scan(self, mixing_scheme_no_compat):
    #     # test the case where the only SCAN entry is for a composition not covered
    #     # by the GGA calculations (i.e., there is no)
    #     # TODO - write this test
    #     pass

    # def test_no_scan_references(self, mixing_scheme_no_compat):
    #     # test the case where none of the SCAN entries correspond to a GGA
    #     # reference structure
    #     lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
    #     lattice2 = Lattice.from_parameters(a=1, b=1, c=0.8, alpha=120, beta=90, gamma=60)
    #     gga_entries = [
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Br"], [[0, 0, 0]]),
    #             -5,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
    #             -10,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(
    #                 lattice,
    #                 ["Sn", "Br", "Br", "Br", "Br"],
    #                 [
    #                     [0, 0, 0],
    #                     [0.2, 0.2, 0.2],
    #                     [0.4, 0.4, 0.4],
    #                     [0.7, 0.7, 0.7],
    #                     [1, 1, 1],
    #                 ],
    #             ),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #     ]
    #     scan_entries = [
    #         # SCAN entries without corresponding SCAN reference structures
    #         ComputedStructureEntry(
    #             Structure(lattice2, ["Sn"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(
    #                 lattice2,
    #                 ["Sn", "Br", "Br"],
    #                 [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]],
    #             ),
    #             -100,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #     ]
    #     compat = MaterialsProjectDFTMixingScheme(compat_1=None)
    #     gga_gs = EntrySet(gga_entries)
    #     gga_gs.remove_non_ground_states()

    #     assert (
    #         compat.get_adjustments(
    #             gga_entries[0],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             [False] * len(gga_entries),
    #             [False] * len(gga_entries),
    #         )
    #         == []
    #     )

    #     with pytest.raises(CompatibilityError, match="there are no SCAN reference"):
    #         compat.get_adjustments(
    #             scan_entries[-1],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             [False] * len(gga_entries),
    #             [False] * len(gga_entries),
    #         )[0].value == pytest.approx(-20)

    #         compat.get_adjustments(
    #             scan_entries[-2],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             [False] * len(gga_entries),
    #             [False] * len(gga_entries),
    #         )[0].value == pytest.approx(-20)

    # def test_scan_polymorph_pair(self, mixing_scheme_no_compat):
    #     # If we have a pair of SCAN calculations at a single composition, one of
    #     # which corresponds to the GGA reference structure, we should correct
    #     # the SCAN entry to GGA and discard the GGA entry
    #     lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
    #     lattice2 = Lattice.from_parameters(a=1, b=1, c=0.8, alpha=120, beta=90, gamma=60)
    #     gga_entries = [
    #         # GGA reference structures
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Br"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
    #             -10,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         # GGA unstable polymorph
    #         ComputedStructureEntry(
    #             Structure(
    #                 lattice2,
    #                 ["Sn", "Br", "Br"],
    #                 [[0, 0, 0], [0.7, 0.7, 0.7], [1, 1, 1]],
    #             ),
    #             -9,
    #             parameters={"run_type": "GGA"},
    #         ),
    #     ]
    #     scan_entries = [
    #         # SCAN reference structure
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
    #             -100,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #         # SCAN unstable polymorph structure, e_above_hull = 2 eV
    #         ComputedStructureEntry(
    #             Structure(
    #                 lattice2,
    #                 ["Sn", "Br", "Br"],
    #                 [[0, 0, 0], [0.7, 0.7, 0.7], [1, 1, 1]],
    #             ),
    #             -98,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #     ]

    #     compat = MaterialsProjectDFTMixingScheme(compat_1=None)
    #     gga_gs = EntrySet(gga_entries)
    #     gga_gs.remove_non_ground_states()
    #     has_scan_ground_states = []
    #     has_scan_hull_entries = []
    #     for e in gga_gs:
    #         if e.composition.reduced_formula == "SnBr2":
    #             has_scan_ground_states.append(True)
    #             has_scan_hull_entries.append(True)
    #         else:
    #             has_scan_ground_states.append(False)
    #             has_scan_hull_entries.append(False)

    #     assert (
    #         compat.get_adjustments(
    #             scan_entries[-1],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             has_scan_ground_states,
    #             has_scan_hull_entries,
    #         )[0].value
    #         == pytest.approx(90)
    #     )

    #     assert (
    #         compat.get_adjustments(
    #             scan_entries[-2],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             has_scan_ground_states,
    #             has_scan_hull_entries,
    #         )[0].value
    #         == pytest.approx(90)
    #     )

    #     assert (
    #         compat.get_adjustments(
    #             gga_entries[0],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             has_scan_ground_states,
    #             has_scan_hull_entries,
    #         )
    #         == []
    #     )

    #     with pytest.raises(CompatibilityError, match="already exists in SCAN"):
    #         assert (
    #             compat.get_adjustments(
    #                 gga_entries[2],
    #                 EntrySet(gga_entries),
    #                 EntrySet(scan_entries),
    #                 gga_gs,
    #                 PhaseDiagram(gga_entries),
    #                 has_scan_ground_states,
    #                 has_scan_hull_entries,
    #             )
    #             == []
    #         )
    #         assert (
    #             compat.get_adjustments(
    #                 gga_entries[3],
    #                 EntrySet(gga_entries),
    #                 EntrySet(scan_entries),
    #                 gga_gs,
    #                 PhaseDiagram(gga_entries),
    #                 has_scan_ground_states,
    #                 has_scan_hull_entries,
    #             )
    #             == []
    #         )

    # def test_wrong_scan_references(self, mixing_scheme_no_compat):
    #     # Raise CompatibilityError if there is no SCAN reference structure
    #     # for the chosen composition
    #     lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
    #     gga_entries = [
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Br"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
    #             -10,
    #             parameters={"run_type": "GGA"},
    #         ),
    #     ]
    #     scan_entries = [
    #         # SCAN reference structure at a different composition
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #         # SCAN entry with a different structure than the GGA reference
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.7, 0.7, 0.7], [1, 1, 1]]),
    #             0,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #     ]

    #     compat = MaterialsProjectDFTMixingScheme(compat_1=None)
    #     gga_gs = EntrySet(gga_entries)
    #     gga_gs.remove_non_ground_states()
    #     has_scan_ground_states = []
    #     has_scan_hull_entries = []
    #     for e in gga_gs:
    #         if e.composition.reduced_formula == "Sn":
    #             has_scan_ground_states.append(True)
    #             has_scan_hull_entries.append(True)
    #         else:
    #             has_scan_ground_states.append(False)
    #             has_scan_hull_entries.append(False)

    #     with pytest.raises(CompatibilityError, match="there is no SCAN reference"):
    #         compat.get_adjustments(
    #             scan_entries[-1],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             has_scan_ground_states,
    #             has_scan_hull_entries,
    #         )

    # def test_none_if_fewer_than_2_scan(self, mixing_scheme_no_compat):
    #     # Raise CompatibilityError if there are fewer than 2 SCAN entries
    #     lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
    #     gga_entries = [
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Br"], [[0, 0, 0]]),
    #             0,
    #             parameters={"run_type": "GGA"},
    #         ),
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
    #             -10,
    #             parameters={"run_type": "GGA"},
    #         ),
    #     ]
    #     scan_entries = [
    #         # one SCAN entry with the same structure than the GGA reference
    #         ComputedStructureEntry(
    #             Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
    #             0,
    #             parameters={"run_type": "R2SCAN"},
    #         ),
    #     ]

    #     compat = MaterialsProjectDFTMixingScheme(compat_1=None)
    #     gga_gs = EntrySet(gga_entries)
    #     gga_gs.remove_non_ground_states()
    #     has_scan_ground_states = []
    #     has_scan_hull_entries = []
    #     for e in gga_gs:
    #         if e.composition.reduced_formula == "Sn":
    #             has_scan_ground_states.append(True)
    #             has_scan_hull_entries.append(True)
    #         else:
    #             has_scan_ground_states.append(False)
    #             has_scan_hull_entries.append(False)

    #     with pytest.raises(CompatibilityError, match="fewer than 2 SCAN"):
    #         compat.get_adjustments(
    #             scan_entries[0],
    #             EntrySet(gga_entries),
    #             EntrySet(scan_entries),
    #             gga_gs,
    #             PhaseDiagram(gga_entries),
    #             has_scan_ground_states,
    #             has_scan_hull_entries,
    #         )

    # def test_get_adjustments(self, mixing_scheme_no_compat):
    #     """
    #     Unit tests of the get_adjustments function internal to the mixing scheme
    #     This function calculates the appropriate energy correction for a single
    #     entry, given a DataFrame of information about the overall collection
    #     of entries being mixed.

    # The first 3 columns are informational only and not used by get_adjustments,
    # so they are populated with dummy values here
    # """
    # compat = MaterialsProjectDFTMixingScheme(compat_1=None)
    # columns = [
    #     "composition",
    #     "spacegroup",
    #     "num_sites",
    #     "run_type_1",
    #     "run_type_2",
    #     "ground_state_energy_1",
    #     "ground_state_energy_2",
    #     "is_stable_1",
    #     "hull_energy_1",
    #     "hull_energy_2",
    # ]
    # # Mixing state 1 - we have all run_type_1 ground states in run_type_2
    # # The energy hull is built with run_type_2 energies
    # # the R2SCAN hull is exactly 10 eV/atom below the GGA in this example
    # row_list = [
    #     [Composition("Ti"), 194, 2, "GGA", "R2SCAN", -7.89, -17.89, True, -7.89, -17.89],
    #     [Composition("O2"), 12, 8, "GGA", "R2SCAN", -4.95, -14.95, True, -4.95, -14.95],
    #     [Composition("Ti3O5"), 12, 16, "GGA", "R2SCAN", -9.01, -19.01, True, -9.01, -19.01],
    #     [Composition("Ti3O4"), 139, 7, "GGA", "R2SCAN", -8.98, -18.98, False, -9, -19],
    #     [Composition("Ti3O"), 149, 24, "GGA", "R2SCAN", -8.55, -18.55, True, -8.55, -18.55],
    #     [Composition("TiO"), 189, 6, "GGA", "R2SCAN", -8.99, -18.99, True, -8.99, -18.99],
    #     [Composition("Ti2O3"), 167, 10, "GGA", "R2SCAN", -9.02, -19.02, True, -9.02, -19.02],
    #     [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
    # ]
    # mixing_state = pd.DataFrame(row_list, columns=columns)

    # # correction for a run_type_1 entry
    # entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
    # with pytest.raises(CompatibilityError, match="Discarding GGA entry"):
    #     adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

    # # correction for a run_type_2 entry
    # entry = ComputedEntry("TiO2", -15, parameters={"run_type": "R2SCAN"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj == []

    # # Mixing state 2 - We have all run_type_1 stable entries (but not all
    # # ground states) in run_type_2. The energy hull is built with run_type_2 energies
    # # the R2SCAN hull is exactly 10 eV/atom below the GGA in this example
    # row_list = [
    #     [Composition("Ti"), 194, 2, "GGA", "R2SCAN", -7.89, -17.89, True, -7.89, -17.89],
    #     [Composition("O2"), 12, 8, "GGA", "R2SCAN", -4.95, -14.95, True, -4.95, -14.95],
    #     [Composition("Ti3O5"), 12, 16, "GGA", "R2SCAN", -9.01, -19.01, True, -9.01, -19.01],
    #     [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, -19],
    #     [Composition("Ti3O"), 149, 24, "GGA", "R2SCAN", -8.55, -18.55, True, -8.55, -18.55],
    #     [Composition("TiO"), 189, 6, "GGA", "R2SCAN", -8.99, -18.99, True, -8.99, -18.99],
    #     [Composition("Ti2O3"), 167, 10, "GGA", "R2SCAN", -9.02, -19.02, True, -9.02, -19.02],
    #     [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
    # ]
    # mixing_state = pd.DataFrame(row_list, columns=columns)

    # # correction for a run_type_1 entry
    # entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
    # with pytest.raises(CompatibilityError, match="Discarding GGA entry"):
    #     adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

    # # correction for a run_type_1 entry where there is no run_type_2 composition
    # entry = ComputedEntry("Ti3O4", -15, parameters={"run_type": "GGA"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
    # assert np.allclose(adj[0].value, -10 * 7)

    # # correction for a run_type_2 entry
    # entry = ComputedEntry("TiO2", -15, parameters={"run_type": "R2SCAN"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj == []

    # # Mixing state 3 - We have run_type_2 for only a few ground states
    # row_list = [
    #     [Composition("Ti"), 194, 2, "GGA", "R2SCAN", -7.89, -17.89, True, -7.89, -17.89],
    #     [Composition("O2"), 12, 8, "GGA", "R2SCAN", -4.95, -14.95, True, -4.95, -14.95],
    #     [Composition("Ti3O5"), 12, 16, "GGA", "R2SCAN", -9.01, -19.01, True, -9.01, -19.01],
    #     [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, None],
    #     [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
    #     [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
    #     [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
    #     [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
    # ]
    # mixing_state = pd.DataFrame(row_list, columns=columns)

    # # correction for a run_type_1 entry
    # # no correction
    # entry = ComputedEntry("Ti3O", -8.55 * 4, parameters={"run_type": "GGA"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj == []

    # # correction for a run_type_2 entry that has a run_type_2 ground state
    # # this energy should be corrected
    # entry = ComputedEntry("Ti3O5", -15 * 8, parameters={"run_type": "R2SCAN"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
    # # r2SCAN entry is 4.01 eV/atom above hull. Add to GGA hull energy
    # assert np.allclose(adj[0].value, (-9.01 + 4.01 - -15) * 8)

    # # correction for a run_type_2 entry that does not have a run_type_2 ground state
    # # this entry should be discarded
    # entry = ComputedEntry("Ti3O4", -15 * 7, parameters={"run_type": "R2SCAN"})
    # with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
    #     adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

    # # Mixing state 4 - We have run_type_2 for only 2 structures, at a single composition
    # row_list = [
    #     [Composition("Ti"), 194, 2, "GGA", None, -7.89, np.nan, True, -7.89, None],
    #     [Composition("O2"), 12, 8, "GGA", None, -4.95, np.nan, True, -4.95, None],
    #     [Composition("Ti3O5"), 12, 16, "GGA", None, -9.01, np.nan, True, -9.01, None],
    #     [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, None],
    #     [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
    #     [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
    #     [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
    #     [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
    #     [Composition("TiO2"), 189, 12, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
    # ]
    # mixing_state = pd.DataFrame(row_list, columns=columns)

    # # correction for a run_type_1 entry at this composition
    # # no correction
    # entry = ComputedEntry("Ti2O4", -8.97 * 6, parameters={"run_type": "GGA"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj == []

    # # correction for a run_type_2 corresponding to the reference state
    # # should be discarded
    # entry = ComputedEntry("Ti2O4", -18.97 * 6, parameters={"run_type": "R2SCAN"})
    # with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
    #     adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

    # # correction for a run_type_2 entry that does not have a run_type_2 ground state
    # # this entry should be corrected to the GGA hull
    # entry = ComputedEntry("Ti4O8", -15.97 * 12, parameters={"run_type": "R2SCAN"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
    # # r2SCAN entry is 3 eV/atom above hull. Add to GGA hull energy
    # assert np.allclose(adj[0].value, (-8.97 + 3 - entry.energy_per_atom) * 12)

    # # Mixing state 5 - We have only 1 run_type_2 material and it is a reference state
    # # This is an edge case that should never happen, because _generate_mixing_scheme_state_data
    # # Will discard run_type_2 entries if there are fewer than 2 of them, so they'll never be in the DataFrame
    # row_list = [
    #     [Composition("Ti"), 194, 2, "GGA", None, -7.89, np.nan, True, -7.89, None],
    #     [Composition("O2"), 12, 8, "GGA", None, -4.95, np.nan, True, -4.95, None],
    #     [Composition("Ti3O5"), 12, 16, "GGA", None, -9.01, np.nan, True, -9.01, None],
    #     [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, None],
    #     [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
    #     [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
    #     [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
    #     [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
    # ]
    # mixing_state = pd.DataFrame(row_list, columns=columns)

    # # correction for a run_type_1 entry
    # entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj == []

    # # correction for the single run_type_2 entry that is the reference state
    # entry = ComputedEntry("TiO2", -18.97 * 3, parameters={"run_type": "R2SCAN"})
    # with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
    #     adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

    # # correction for a run_type_2 entry at the same composition that was not
    # # included in the original list of processed entries (Another edge case that should
    # # never happen)
    # entry = ComputedEntry("TiO2", -17.97 * 3, parameters={"run_type": "R2SCAN"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
    # # r2SCAN entry is 1 eV/atom above hull. Add to GGA hull energy
    # assert np.allclose(adj[0].value, (-8.97 + 1 - entry.energy_per_atom) * 3)

    # # Mixing state 6 - We have 0 run_type_2 materials
    # row_list = [
    #     [Composition("Ti"), 194, 2, "GGA", None, -7.89, np.nan, True, -7.89, None],
    #     [Composition("O2"), 12, 8, "GGA", None, -4.95, np.nan, True, -4.95, None],
    #     [Composition("Ti3O5"), 12, 16, "GGA", None, -9.01, np.nan, True, -9.01, None],
    #     [Composition("Ti3O4"), 139, 7, "GGA", None, -8.98, np.nan, False, -9, None],
    #     [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
    #     [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
    #     [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
    #     [Composition("TiO2"), 141, 6, "GGA", None, -8.97, np.nan, True, -8.97, None],
    # ]
    # mixing_state = pd.DataFrame(row_list, columns=columns)

    # # correction for a run_type_1 entry
    # entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj == []

    # # correction for a run_type_2 entry
    # entry = ComputedEntry("TiO2", -15, parameters={"run_type": "R2SCAN"})
    # with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
    #     adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

    # # Mixing state 7 - We have run_type_2_materials for everything, including
    # # a composition that is not present in run_type_1
    # row_list = [
    #     [Composition("Ti"), 194, 2, "GGA", "R2SCAN", -7.89, -17.89, True, -7.89, -17.89],
    #     [Composition("O2"), 12, 8, "GGA", "R2SCAN", -4.95, -14.95, True, -4.95, -14.95],
    #     [Composition("Ti3O5"), 12, 16, "GGA", "R2SCAN", -9.01, -19.01, True, -9.01, -19.01],
    #     [Composition("Ti3O4"), 139, 7, None, "R2SCAN", np.nan, -18.98, False, -9, -19],
    #     [Composition("Ti3O"), 149, 24, "GGA", "R2SCAN", -8.55, -18.55, True, -8.55, -18.55],
    #     [Composition("TiO"), 189, 6, "GGA", "R2SCAN", -8.99, -18.99, True, -8.99, -18.99],
    #     [Composition("Ti2O3"), 167, 10, "GGA", "R2SCAN", -9.02, -19.02, True, -9.02, -19.02],
    #     [Composition("TiO2"), 141, 6, "GGA", "R2SCAN", -8.97, -18.97, True, -8.97, -18.97],
    # ]
    # mixing_state = pd.DataFrame(row_list, columns=columns)

    # # correction for a run_type_1 entry
    # entry = ComputedEntry("TiO2", -8.97 * 3, parameters={"run_type": "GGA"})
    # with pytest.raises(CompatibilityError, match="Discarding GGA entry"):
    #     adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

    # # correction for a run_type_2 entry with no corresponding run_type_1
    # entry = ComputedEntry("Ti3O4", -8.97 * 3, parameters={"run_type": "R2SCAN"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj == []

    # # correction for a run_type_2 entry
    # entry = ComputedEntry("TiO2", -15, parameters={"run_type": "R2SCAN"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj == []

    # # Mixing state 8 - We have run_type_2 for only 2 structures, at a single composition
    # # that is not on the run_type_1 hull
    # row_list = [
    #     [Composition("Ti"), 194, 2, "GGA", None, -7.89, np.nan, True, -7.89, None],
    #     [Composition("O2"), 12, 8, "GGA", None, -4.95, np.nan, True, -4.95, None],
    #     [Composition("Ti3O5"), 12, 16, "GGA", None, -9.01, np.nan, True, -9.01, None],
    #     [Composition("Ti3O4"), 139, 7, "GGA", "R2SCAN", -8.98, -18.98, False, -9, None],
    #     [Composition("Ti3O4"), 149, 14, None, "R2SCAN", -8.98, -18.98, False, -9, None],
    #     [Composition("Ti3O"), 149, 24, "GGA", None, -8.55, np.nan, True, -8.55, None],
    #     [Composition("TiO"), 189, 6, "GGA", None, -8.99, np.nan, True, -8.99, None],
    #     [Composition("Ti2O3"), 167, 10, "GGA", None, -9.02, np.nan, True, -9.02, None],
    # ]
    # mixing_state = pd.DataFrame(row_list, columns=columns)

    # # correction for a run_type_1 entry at this composition
    # # no correction
    # entry = ComputedEntry("Ti3O4", -8.98 * 7, parameters={"run_type": "GGA"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj == []

    # # correction for a run_type_2 corresponding to the reference state
    # # should be discarded
    # entry = ComputedEntry("Ti3O4", -18.98 * 7, parameters={"run_type": "R2SCAN"})
    # with pytest.raises(CompatibilityError, match="Discarding R2SCAN entry"):
    #     adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)

    # # correction for a run_type_2 entry that does not have a run_type_2 ground state
    # # this entry should be corrected to the GGA hull
    # entry = ComputedEntry("Ti6O8", -16.98 * 14, parameters={"run_type": "R2SCAN"})
    # adj = compat.get_adjustments(entry, mixing_scheme_state_data=mixing_state)
    # assert adj[0].name == "MP GGA(+U)/R2SCAN mixing adjustment"
    # # r2SCAN entry is 2 eV/atom above hull. Add to GGA hull energy
    # assert np.allclose(adj[0].value, (-9 + 2 - entry.energy_per_atom) * 14)


# if __name__ == "__main__":
#     # import sys;sys.argv = ['', 'Test.testName']
#     unittest.main()
