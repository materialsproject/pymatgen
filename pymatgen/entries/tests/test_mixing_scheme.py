# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Tests for the Materials Project DFT mixing scheme

**NOTE FOR FUTURE DEVELOPERS**

PLEASE DO NOT CHANGE THESE TESTS WITHOUT READING THIS ENTIRE DOCUMENTATION!

The MP DFT mixing scheme is conceptually simple but highly intricate to implement.
The tests in this file were crafted carefully over a period of many months
in order to capture as many nuances of the desired behavior as possible, and were
integral to developing the mixing code itself.

Test Structure
--------------

The majority of the tests use "mixing states" to check behavior. Mixing states
are merely combinations of different ComputedStructureEntry with different run_type
- e.g. all ground states present for both run_types, only one run_type, etc. Mixing
states are defined using the `MixingState` utility class, which has attributes that
return 1) the GGA entries, 2) the R2SCAN entries, 3) all entries, and 4) the pandas
DataFrame that represents the mixing state. Note that 'GGA' and 'R2SCAN' are used
throughout this test file to represent run_type_1 and run_type_2, respectively,
but the mixing scheme is design to be able to mix any two functionals.

Most mixing states are subsets of the `ms_complete` mixing state. `ms_complete`
was crafted to capture most of the scenarios that may be encountered when mixing
ComputedStructureEntry from different functionals. It comprises a complete binary
phase diagram, with all entries present as both GGA and R2SCAN calculations.

Note that these entries are inspired by, but NOT equivalent to, the real SnBr2 phase
diagram. Rather than use real energies or structures, arbitrary energies and structures
have been used to keep this test file cleaner and easier to understand. The Bromine
structures are the one exception to this. These structures are taken from real calculations
and will fail to structure match with default arguments (they are included here to
test "fuzzy matching" behavior). Entry-id's are assigned numerically, with the
same number corresponding to equivalent materials. So "gga-1" should have the
same structure as "r2scan-1".

Description of the `ms_complete` mixing state
---------------------------------------------

Future developers are HIGHLY encouraed to plot the PhaseDiagram associated with both
the R2SCAN and the GGA entries in `ms_complete` before attempting to modify any of these
tests.

**GGA entries**

- gga-1: stable ground state of Sn
- gga-2: unstable polymorph of Br
- gga-3: stable polymorph of Br
- gga-4: stable polymorph of SnBr2
- gga-5: unstable polymorph of SnBr2, 1 eV/atom above hull
- gga-6: unstable polymorph of SnBr2, 2 eV/atom above hull
- gga-7: unstable composition (SnBr4), 0.6 eV/atom above hull

**R2SCAN entries**

- All the same entries exist as in GGA
- Entries with corresponding numbers have matching structures
  (i.e. gga-1 and r2scan-1)
- Unless otherwise listed below, energies are 1 eV/atom lower than those in GGA
- for Br, the GGA ground state structure gga-3 does not match r2scan-3 unless you
  use fuzzy matching
- for SnBr2, a different polymorph is stabilized than in GGA (r2scan-5 whereas
  r2scan-4 was the GGA ground state)
- entry r2scan-6 (unstable SnBr2 polymorph) is scaled to 25% the size of gga-6.
  This will match with default StructureMatcher settings but not with customized
  settings.
- SnBr4 (r2scan-7) appears on the hull whereas it is not stable in GGA.
- A few mixing states (but not `ms_complete`) also include an unstable polymorph
  of SnBr4, entry r2scan-8, that does not exist in GGA.

Types of Tests
--------------

In general there are 3 types of tests. Most tests follow this pattern.

1. Unit tests of get_mixing_state_data, to verify that it generates
   the correct DataFrame when provided with a particular list of entries
2. Unit tests that verify get_adjustments behaves correctly, when provided with
   a ComputedEntry and a pre-generated mixing_scheme_state_data DataFrame
3. Functional tests of the behavior of the process_entries method, which is the
   primary user-facing interface

Implementation Notes
--------------------

- Tests are organized into two classes. One class collects tests of the different
  args / kwargs that can be passed to the mixing scheme. The other class collects
  tests of different mixing states.
- Tests are written in pure pytest format, so the class organization is merely
  a convenience to facilitate code folding, etc.
- pytest fixtures are used to define the various mixing states. Using pytest
  fixtures is helpful here since it ensures every test receives fresh, unmodified
  copies of the respective entries. process_entries modifies entry energies
  in place, so tests could cross-contaminate one another if a fixture were not used.
"""

__author__ = "Ryan Kingsbury"
__copyright__ = "Copyright 2019-2021, The Materials Project"
__version__ = "1.0"
__email__ = "RKingsbury@lbl.gov"
__date__ = "October 2021"

import numpy as np
import pandas as pd
import pytest

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.entries.compatibility import Compatibility, CompatibilityError
from pymatgen.entries.computed_entries import (
    CompositionEnergyAdjustment,
    ComputedEntry,
    ComputedStructureEntry,
)
from pymatgen.entries.mixing_scheme import MaterialsProjectDFTMixingScheme

"""
Define utility classes to make the tests easier to read
"""


class DummyCompatibility(Compatibility):
    """
    Dummy class to test compat1 and compat2 kwargs
    """

    def get_adjustments(self, entry):
        return [CompositionEnergyAdjustment(-10, entry.composition.num_atoms, name="Dummy adjustment")]


class MixingState:
    """
    Lightweight container class for a mixing state, used by the tests below
    """

    def __init__(self, gga_entries, scan_entries, dataframe):
        """
        Args:
            gga_entries: list of run_type_1 (usually GGA or GGA+U) ComputedEntry
            scan_entries: list of run_type_2 (usually R2SCAN) ComputedEntry
            DataFrame: pandas DataFrame representing the mixing state, of the
                format returned by get_mixing_state_data
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
Define mixing states to test.
"""
# columns must correspond to the columns in get_mixing_state_data
columns = [
    "formula",
    "spacegroup",
    "num_sites",
    "is_stable_1",
    "entry_id_1",
    "entry_id_2",
    "run_type_1",
    "run_type_2",
    "energy_1",
    "energy_2",
    "hull_energy_1",
    "hull_energy_2",
]

# lattices. In general, ground states are all assigned lattice1
# unstable polymorphs are assigned lattice2 or lattice 3
lattice1 = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
lattice2 = Lattice.from_parameters(a=1, b=1, c=1, alpha=120, beta=120, gamma=60)
lattice3 = Lattice.from_parameters(a=1, b=1, c=1, alpha=120, beta=120, gamma=90)
lattice_br_gga = Lattice.from_dict(
    {
        "@module": "pymatgen.core.lattice",
        "@class": "Lattice",
        "matrix": [[2.129324, -4.226095, 0.0], [2.129324, 4.226095, 0.0], [0.0, 0.0, 8.743796]],
    }
)
lattice_br_r2scan = Lattice.from_dict(
    {
        "@module": "pymatgen.core.lattice",
        "@class": "Lattice",
        "matrix": [[0.0, -4.25520892, -0.0], [-3.56974866, 2.12760446, 0.0], [0.0, 0.0, -8.74536848]],
    }
)


@pytest.fixture
def ms_complete():
    """
    Mixing state where we have R2SCAN for all GGA
    """

    gga_entries = [
        ComputedStructureEntry(
            Structure(lattice1, ["Sn"], [[0, 0, 0]]), 0, parameters={"run_type": "GGA"}, entry_id="gga-1"
        ),
        ComputedStructureEntry(
            Structure(lattice1, ["Br"], [[0, 0, 0]]), 1, parameters={"run_type": "GGA"}, entry_id="gga-2"
        ),
        ComputedStructureEntry(
            Structure(
                lattice_br_gga,
                ["Br", "Br", "Br", "Br"],
                [
                    [0.642473, 0.642473, 0.117751],
                    [0.357527, 0.357527, 0.882249],
                    [0.857527, 0.857527, 0.617751],
                    [0.142473, 0.142473, 0.382249],
                ],
            ),
            0,
            parameters={"run_type": "GGA"},
            entry_id="gga-3",
        ),
        ComputedStructureEntry(
            Structure(lattice1, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -18,
            parameters={"run_type": "GGA"},
            entry_id="gga-4",
        ),
        ComputedStructureEntry(
            Structure(
                lattice2,
                ["Sn", "Sn", "Sn", "Sn", "Br", "Br", "Br", "Br", "Br", "Br", "Br", "Br"],
                [
                    [0.25, 0.393393, 0.663233],
                    [0.75, 0.606607, 0.336767],
                    [0.25, 0.893393, 0.836767],
                    [0.75, 0.106607, 0.163233],
                    [0.25, 0.662728, 0.548755],
                    [0.75, 0.337272, 0.451245],
                    [0.25, 0.162728, 0.951245],
                    [0.75, 0.837272, 0.048755],
                    [0.25, 0.992552, 0.311846],
                    [0.75, 0.007448, 0.688154],
                    [0.25, 0.492552, 0.188154],
                    [0.75, 0.507448, 0.811846],
                ],
            ),
            -60,
            parameters={"run_type": "GGA"},
            entry_id="gga-5",
        ),
        ComputedStructureEntry(
            Structure(lattice3, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -12,
            parameters={"run_type": "GGA"},
            entry_id="gga-6",
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
            entry_id="gga-7",
        ),
    ]
    scan_entries = [
        ComputedStructureEntry(
            Structure(lattice1, ["Sn"], [[0, 0, 0]]), -1, parameters={"run_type": "R2SCAN"}, entry_id="r2scan-1"
        ),
        ComputedStructureEntry(
            Structure(lattice1, ["Br"], [[0, 0, 0]]), -1, parameters={"run_type": "R2SCAN"}, entry_id="r2scan-2"
        ),
        ComputedStructureEntry(
            Structure(
                lattice_br_r2scan,
                ["Br", "Br", "Br", "Br"],
                [
                    [0.85985939, 0.0, 0.38410868],
                    [0.14014061, -0.0, 0.61589132],
                    [0.64014061, 0.0, 0.88410868],
                    [0.35985939, -0.0, 0.11589132],
                ],
            ),
            0,
            parameters={"run_type": "R2SCAN"},
            entry_id="r2scan-3",
        ),
        ComputedStructureEntry(
            Structure(lattice1, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -21,
            parameters={"run_type": "R2SCAN"},
            entry_id="r2scan-4",
        ),
        ComputedStructureEntry(
            Structure(
                lattice2,
                ["Sn", "Sn", "Sn", "Sn", "Br", "Br", "Br", "Br", "Br", "Br", "Br", "Br"],
                [
                    [0.25, 0.393393, 0.663233],
                    [0.75, 0.606607, 0.336767],
                    [0.25, 0.893393, 0.836767],
                    [0.75, 0.106607, 0.163233],
                    [0.25, 0.662728, 0.548755],
                    [0.75, 0.337272, 0.451245],
                    [0.25, 0.162728, 0.951245],
                    [0.75, 0.837272, 0.048755],
                    [0.25, 0.992552, 0.311846],
                    [0.75, 0.007448, 0.688154],
                    [0.25, 0.492552, 0.188154],
                    [0.75, 0.507448, 0.811846],
                ],
            ),
            -96,
            parameters={"run_type": "R2SCAN"},
            entry_id="r2scan-5",
        ),
        ComputedStructureEntry(
            Structure(lattice3.scale(0.25), ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -18,
            parameters={"run_type": "R2SCAN"},
            entry_id="r2scan-6",
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
            entry_id="r2scan-7",
        ),
    ]
    # the fmt command tells the black autoformatter not to mess with this block of code
    # it's easier to edit when all the commas are lined up.
    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3", "r2scan-3", "GGA", "R2SCAN",  0.,  0.,   0., -1.],
        ["Br",   191,  1, False, "gga-2", "r2scan-2", "GGA", "R2SCAN",  1., -1.,   0., -1.],
        ["Sn",   191,  1,  True, "gga-1", "r2scan-1", "GGA", "R2SCAN",  0., -1.,   0., -1.],
        ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN", -6., -7.,  -6., -8.],
        ["SnBr2",  2, 12, False, "gga-5", "r2scan-5", "GGA", "R2SCAN", -5., -8.,  -6., -8.],
        ["SnBr2", 71,  3, False, "gga-6", "r2scan-6", "GGA", "R2SCAN", -4., -6.,  -6., -8.],
        ["SnBr4",  8,  5, False, "gga-7", "r2scan-7", "GGA", "R2SCAN", -3., -6., -3.6, -6.],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_scan_only(ms_complete):
    """
    Mixing state with only R2SCAN entries
    """
    gga_entries = []
    scan_entries = ms_complete.scan_entries

    # fmt: off
    row_list = [
        ["Br",    64,  4, False, None, "r2scan-3", None, "R2SCAN", np.nan,  0., np.nan, -1.],
        ["Br",   191,  1, False, None, "r2scan-2", None, "R2SCAN", np.nan, -1., np.nan, -1.],
        ["Sn",   191,  1, False, None, "r2scan-1", None, "R2SCAN", np.nan, -1., np.nan, -1.],
        ["SnBr2",  2, 12, False, None, "r2scan-5", None, "R2SCAN", np.nan, -8., np.nan, -8.],
        ["SnBr2", 65,  3, False, None, "r2scan-4", None, "R2SCAN", np.nan, -7., np.nan, -8.],
        ["SnBr2", 71,  3, False, None, "r2scan-6", None, "R2SCAN", np.nan, -6., np.nan, -8.],
        ["SnBr4",  8,  5, False, None, "r2scan-7", None, "R2SCAN", np.nan, -6., np.nan, -6.],
    ]
    # fmt: on

    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_gga_only(ms_complete):
    """
    Mixing state with only GGA entries
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = []

    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3", None, "GGA", None,  0., np.nan,   0., np.nan],
        ["Br",   191,  1, False, "gga-2", None, "GGA", None,  1., np.nan,   0., np.nan],
        ["Sn",   191,  1,  True, "gga-1", None, "GGA", None,  0., np.nan,   0., np.nan],
        ["SnBr2", 65,  3,  True, "gga-4", None, "GGA", None, -6., np.nan,  -6., np.nan],
        ["SnBr2",  2, 12, False, "gga-5", None, "GGA", None, -5., np.nan,  -6., np.nan],
        ["SnBr2", 71,  3, False, "gga-6", None, "GGA", None, -4., np.nan,  -6., np.nan],
        ["SnBr4",  8,  5, False, "gga-7", None, "GGA", None, -3., np.nan, -3.6, np.nan],
    ]
    # fmt: on

    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_gga_1_scan(ms_complete):
    """
    Mixing state with all GGA entries and one R2SCAN, corresponding to the GGA
    ground state of SnBr2 (r2scan-4)
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = [e for e in ms_complete.scan_entries if e.entry_id == "r2scan-4"]

    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3",       None, "GGA",     None,  0., np.nan,   0., np.nan],
        ["Br",   191,  1, False, "gga-2",       None, "GGA",     None,  1., np.nan,   0., np.nan],
        ["Sn",   191,  1,  True, "gga-1",       None, "GGA",     None,  0., np.nan,   0., np.nan],
        ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN", -6.,     -7,  -6., np.nan],
        ["SnBr2",  2, 12, False, "gga-5",       None, "GGA",     None, -5., np.nan,  -6., np.nan],
        ["SnBr2", 71,  3, False, "gga-6",       None, "GGA",     None, -4., np.nan,  -6., np.nan],
        ["SnBr4",  8,  5, False, "gga-7",       None, "GGA",     None, -3., np.nan, -3.6, np.nan],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)
    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_gga_1_scan_novel(ms_complete):
    """
    Mixing state with all GGA entries and 1 R2SCAN, corresponding to a composition
    (SnBr) that is not present in the GGA entries.
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = [
        ComputedStructureEntry(
            Structure(
                lattice1,
                ["Sn", "Sn", "Br", "Br"],
                [
                    [0, 0, 0],
                    [0.2, 0.2, 0.2],
                    [0.4, 0.4, 0.4],
                    [0.7, 0.7, 0.7],
                ],
            ),
            -20,
            parameters={"run_type": "R2SCAN"},
            entry_id="r2scan-9",
        ),
    ]

    # fmt: off
    row_list = [
        ["Br",    64,  4,   True, "gga-3",       None, "GGA",     None,      0., np.nan,   0., np.nan],
        ["Br",   191,  1,  False, "gga-2",       None, "GGA",     None,      1., np.nan,   0., np.nan],
        ["Sn",   191,  1,   True, "gga-1",       None, "GGA",     None,      0., np.nan,   0., np.nan],
        ["SnBr",   8,  4,  False,    None, "r2scan-9",  None, "R2SCAN",  np.nan,    -5., -4.5, np.nan],
        ["SnBr2", 65,  3,   True, "gga-4",       None, "GGA",     None,     -6., np.nan,  -6., np.nan],
        ["SnBr2",  2, 12,  False, "gga-5",       None, "GGA",     None,     -5., np.nan,  -6., np.nan],
        ["SnBr2", 71,  3,  False, "gga-6",       None, "GGA",     None,     -4., np.nan,  -6., np.nan],
        ["SnBr4",  8,  5,  False, "gga-7",       None, "GGA",     None,     -3., np.nan, -3.6, np.nan],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)
    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_gga_2_scan_same(ms_complete):
    """
    Mixing state with all GGA entries and 2 R2SCAN, corresponding to the GGA
    ground state and one unstable polymoprh of SnBr2 (r2scan-4 and r2scan-6)
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = [e for e in ms_complete.scan_entries if e.entry_id in ["r2scan-4", "r2scan-6"]]

    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3",       None, "GGA",     None,  0., np.nan,   0., np.nan],
        ["Br",   191,  1, False, "gga-2",       None, "GGA",     None,  1., np.nan,   0., np.nan],
        ["Sn",   191,  1,  True, "gga-1",       None, "GGA",     None,  0., np.nan,   0., np.nan],
        ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN", -6.,     -7,  -6., np.nan],
        ["SnBr2",  2, 12, False, "gga-5",       None, "GGA",     None, -5., np.nan,  -6., np.nan],
        ["SnBr2", 71,  3, False, "gga-6", "r2scan-6", "GGA", "R2SCAN", -4.,     -6,  -6., np.nan],
        ["SnBr4",  8,  5, False, "gga-7",       None, "GGA",     None, -3., np.nan, -3.6, np.nan],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)
    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_gga_2_scan_diff_match(ms_complete):
    """
    Mixing state with all GGA entries and 2 R2SCAN entries corresponding to
    different compositions, where both R2SCAN materials match GGA materials, but
    only one matches a GGA ground state.
    r2scan-4 and r2scan-7
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = [e for e in ms_complete.scan_entries if e.entry_id in ["r2scan-4", "r2scan-7"]]

    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3",       None, "GGA",     None,  0., np.nan,   0., np.nan],
        ["Br",   191,  1, False, "gga-2",       None, "GGA",     None,  1., np.nan,   0., np.nan],
        ["Sn",   191,  1,  True, "gga-1",       None, "GGA",     None,  0., np.nan,   0., np.nan],
        ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN", -6.,     -7,  -6., np.nan],
        ["SnBr2",  2, 12, False, "gga-5",       None, "GGA",     None, -5., np.nan,  -6., np.nan],
        ["SnBr2", 71,  3, False, "gga-6",       None, "GGA",     None, -4., np.nan,  -6., np.nan],
        ["SnBr4",  8,  5, False, "gga-7", "r2scan-7", "GGA", "R2SCAN", -3.,     -6, -3.6, np.nan],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)
    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_gga_2_scan_diff_no_match(ms_complete):
    """
    Mixing state with all GGA entries and 2 R2SCAN, corresponding to the GGA
    ground state of SnBr2 (r2scan-4) and one unstable polymoprh of SnBr4
    that does not match any GGA material (r2scan-8)
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = [e for e in ms_complete.scan_entries if e.entry_id in ["r2scan-4"]]
    scan_entries.append(
        ComputedStructureEntry(
            Structure(
                lattice3,
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
            entry_id="r2scan-8",
        ),
    )

    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3",       None, "GGA",     None,      0., np.nan,   0., np.nan],
        ["Br",   191,  1, False, "gga-2",       None, "GGA",     None,      1., np.nan,   0., np.nan],
        ["Sn",   191,  1,  True, "gga-1",       None, "GGA",     None,      0., np.nan,   0., np.nan],
        ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN",     -6.,    -7.,  -6., np.nan],
        ["SnBr2",  2, 12, False, "gga-5",       None, "GGA",     None,     -5., np.nan,  -6., np.nan],
        ["SnBr2", 71,  3, False, "gga-6",       None, "GGA",     None,     -4., np.nan,  -6., np.nan],
        ["SnBr4",  8,  5, False, "gga-7",       None, "GGA",     None,     -3., np.nan, -3.6, np.nan],
        ["SnBr4", 44,  5, False,    None, "r2scan-8",  None, "R2SCAN",  np.nan,    -5., -3.6, np.nan],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)
    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_all_gga_scan_gs(ms_complete):
    """
    Mixing state with all GGA entries and R2SCAN entries corresponding to all GGA
    ground states, but no others.
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = [e for e in ms_complete.scan_entries if e.entry_id in ["r2scan-1", "r2scan-3", "r2scan-4"]]

    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3", "r2scan-3", "GGA", "R2SCAN",      0.,     0.,   0.,    0.],
        ["Br",   191,  1, False, "gga-2",       None, "GGA",     None,      1., np.nan,   0.,    0.],
        ["Sn",   191,  1,  True, "gga-1", "r2scan-1", "GGA", "R2SCAN",      0.,    -1.,   0.,   -1.],
        ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN",     -6.,    -7.,  -6.,   -7.],
        ["SnBr2",  2, 12, False, "gga-5",       None, "GGA",     None,     -5., np.nan,  -6.,   -7.],
        ["SnBr2", 71,  3, False, "gga-6",       None, "GGA",     None,     -4., np.nan,  -6.,   -7.],
        ["SnBr4",  8,  5, False, "gga-7",       None, "GGA",     None,     -3., np.nan, -3.6,  -4.2],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)
    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_all_gga_scan_gs_plus_novel(ms_all_gga_scan_gs):
    """
    Mixing state with all GGA entries and R2SCAN entries corresponding to all GGA
    ground states, plus one R2SCAN entry at a novel composition not in the GGA
    phase diagram
    """
    gga_entries = ms_all_gga_scan_gs.gga_entries
    scan_entries = ms_all_gga_scan_gs.scan_entries
    scan_entries.append(
        ComputedStructureEntry(
            Structure(
                lattice1,
                ["Sn", "Sn", "Br", "Br"],
                [
                    [0, 0, 0],
                    [0.2, 0.2, 0.2],
                    [0.4, 0.4, 0.4],
                    [0.7, 0.7, 0.7],
                ],
            ),
            -20,
            parameters={"run_type": "R2SCAN"},
            entry_id="r2scan-9",
        ),
    )

    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3", "r2scan-3", "GGA", "R2SCAN",      0.,     0.,   0.,     0.],
        ["Br",   191,  1, False, "gga-2",       None, "GGA",     None,      1., np.nan,   0.,     0.],
        ["Sn",   191,  1,  True, "gga-1", "r2scan-1", "GGA", "R2SCAN",      0.,    -1.,   0.,    -1.],
        ["SnBr",   8,  4, False,    None, "r2scan-9",  None, "R2SCAN",  np.nan,    -5., -4.5,   -5.5],
        ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN",     -6.,    -7.,  -6.,    -7.],
        ["SnBr2",  2, 12, False, "gga-5",       None, "GGA",     None,     -5., np.nan,  -6.,    -7.],
        ["SnBr2", 71,  3, False, "gga-6",       None, "GGA",     None,     -4., np.nan,  -6.,    -7.],
        ["SnBr4",  8,  5, False, "gga-7",       None, "GGA",     None,     -3., np.nan, -3.6,   -4.2],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)
    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_all_scan_novel(ms_complete):
    """
    Mixing state with all GGA entries and all R2SCAN, with an additional unstable
    polymorphs of SnBr4 (r2scan-8) only in R2SCAN.
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = ms_complete.scan_entries
    scan_entries.append(
        ComputedStructureEntry(
            Structure(
                lattice3,
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
            entry_id="r2scan-8",
        ),
    )

    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3", "r2scan-3", "GGA", "R2SCAN",     0.,  0.,   0., -1.],
        ["Br",   191,  1, False, "gga-2", "r2scan-2", "GGA", "R2SCAN",     1., -1.,   0., -1.],
        ["Sn",   191,  1,  True, "gga-1", "r2scan-1", "GGA", "R2SCAN",     0., -1.,   0., -1.],
        ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN",    -6., -7.,  -6., -8.],
        ["SnBr2",  2, 12, False, "gga-5", "r2scan-5", "GGA", "R2SCAN",    -5., -8.,  -6., -8.],
        ["SnBr2", 71,  3, False, "gga-6", "r2scan-6", "GGA", "R2SCAN",    -4., -6.,  -6., -8.],
        ["SnBr4",  8,  5, False, "gga-7", "r2scan-7", "GGA", "R2SCAN",    -3., -6., -3.6, -6.],
        ["SnBr4",  8,  5, False,    None, "r2scan-8",  None, "R2SCAN", np.nan, -6., -3.6, -6.],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)
    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_incomplete_gga_all_scan(ms_complete):
    """
    Mixing state with an incomplete GGA phase diagram
    """
    gga_entries = [e for e in ms_complete.gga_entries if e.composition.reduced_formula != "Sn"]
    scan_entries = ms_complete.scan_entries

    # fmt: off
    row_list = [
        ["Br",    64,  4, False, "gga-3", "r2scan-3", "GGA", "R2SCAN",      0.,  0., np.nan, -1.],
        ["Br",   191,  1, False, "gga-2", "r2scan-2", "GGA", "R2SCAN",      1., -1., np.nan, -1.],
        ["Sn",   191,  1, False,    None, "r2scan-1",  None, "R2SCAN",  np.nan, -1., np.nan, -1.],
        ["SnBr2", 65,  3, False, "gga-4", "r2scan-4", "GGA", "R2SCAN",     -6., -7., np.nan, -8.],
        ["SnBr2",  2, 12, False, "gga-5", "r2scan-5", "GGA", "R2SCAN",     -5., -8., np.nan, -8.],
        ["SnBr2", 71,  3, False, "gga-6", "r2scan-6", "GGA", "R2SCAN",     -4., -6., np.nan, -8.],
        ["SnBr4",  8,  5, False, "gga-7", "r2scan-7", "GGA", "R2SCAN",     -3., -6., np.nan, -6.],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_scan_chemsys_superset(ms_complete):
    """
    Mixing state where we have R2SCAN for all GGA, and there is an additional R2SCAN
    entry outside the GGA chemsys
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = ms_complete.scan_entries
    scan_entries.append(
        ComputedStructureEntry(
            Structure(
                lattice3,
                ["Sn", "Cl", "Cl", "Br", "Br"],
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
            entry_id="r2scan-9",
        ),
    )

    mixing_state = ms_complete.state_data
    return MixingState(gga_entries, scan_entries, mixing_state)


@pytest.fixture
def ms_complete_duplicate_structs(ms_complete):
    """
    Mixing state where we have R2SCAN for all GGA, plus extra entries that duplicate
    the structures of gga-4 and r2scan-4, and have slightly higher energies
    """
    gga_entries = ms_complete.gga_entries
    scan_entries = ms_complete.scan_entries

    gga_entries.append(
        ComputedStructureEntry(
            Structure(lattice1, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -17.9,
            parameters={"run_type": "GGA"},
            entry_id="gga-10",
        ),
    )

    scan_entries.append(
        ComputedStructureEntry(
            Structure(lattice1, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
            -20.9,
            parameters={"run_type": "R2SCAN"},
            entry_id="r2scan-10",
        ),
    )

    # fmt: off
    row_list = [
        ["Br",    64,  4,  True, "gga-3", "r2scan-3", "GGA", "R2SCAN",  0.,  0.,   0., -1.],
        ["Br",   191,  1, False, "gga-2", "r2scan-2", "GGA", "R2SCAN",  1., -1.,   0., -1.],
        ["Sn",   191,  1,  True, "gga-1", "r2scan-1", "GGA", "R2SCAN",  0., -1.,   0., -1.],
        ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN", -6., -7.,  -6., -8.],
        ["SnBr2",  2, 12, False, "gga-5", "r2scan-5", "GGA", "R2SCAN", -5., -8.,  -6., -8.],
        ["SnBr2", 71,  3, False, "gga-6", "r2scan-6", "GGA", "R2SCAN", -4., -6.,  -6., -8.],
        ["SnBr4",  8,  5, False, "gga-7", "r2scan-7", "GGA", "R2SCAN", -3., -6., -3.6, -6.],
    ]
    # fmt: on
    mixing_state = pd.DataFrame(row_list, columns=columns)

    return MixingState(gga_entries, scan_entries, mixing_state)


def test_data_ms_complete(ms_complete):
    """
    Verify that the test chemical system
    ComputedStructureEntry match (or don't match) as intended
    """
    sm = StructureMatcher()
    for g, s in zip(ms_complete.gga_entries, ms_complete.scan_entries):
        if g.entry_id == "gga-3":
            assert not sm.fit(g.structure, s.structure)
        else:
            assert sm.fit(g.structure, s.structure)


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
            ComputedEntry("Sn", 0, parameters={"run_type": "GGA"}, entry_id=1),
            ComputedEntry("Br", 0, parameters={"run_type": "GGA"}, entry_id=2),
            ComputedEntry("SnBr2", -10, parameters={"run_type": "R2SCAN"}, entry_id=3),
            ComputedEntry("SnBr2", -100, parameters={"run_type": "GGA"}, entry_id=4),
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
                entry_id=5,
            ),
        ]

        with pytest.warns(UserWarning, match="not a ComputedStructureEntry"):
            mixing_scheme_no_compat.process_entries(entries)

    def test_empty_entries(self, mixing_scheme_no_compat):
        # Test behavior when either gga_entries or scan_entries passed to get_adjustments
        # is empty
        entries = []
        mixing_scheme_no_compat.process_entries(entries)

    def test_no_entry_ids(self, mixing_scheme_no_compat):
        """
        unique entry_ids are required.
        """
        gga_entries = [
            ComputedStructureEntry(
                Structure(lattice1, ["Sn"], [[0, 0, 0]]), 0, parameters={"run_type": "GGA"}, entry_id="gga"
            ),
            ComputedStructureEntry(
                Structure(lattice1, ["Br"], [[0, 0, 0]]), 1, parameters={"run_type": "GGA"}, entry_id="gga"
            ),
            ComputedStructureEntry(
                Structure(
                    lattice_br_gga,
                    ["Br", "Br", "Br", "Br"],
                    [
                        [0.642473, 0.642473, 0.117751],
                        [0.357527, 0.357527, 0.882249],
                        [0.857527, 0.857527, 0.617751],
                        [0.142473, 0.142473, 0.382249],
                    ],
                ),
                0,
                parameters={"run_type": "GGA"},
                entry_id="gga-3",
            ),
        ]
        scan_entries = [
            ComputedStructureEntry(
                Structure(lattice1, ["Sn"], [[0, 0, 0]]), -1, parameters={"run_type": "R2SCAN"}, entry_id="r2scan-1"
            ),
            ComputedStructureEntry(
                Structure(lattice1, ["Br"], [[0, 0, 0]]), -1, parameters={"run_type": "R2SCAN"}, entry_id="r2scan-2"
            ),
            ComputedStructureEntry(
                Structure(
                    lattice_br_r2scan,
                    ["Br", "Br", "Br", "Br"],
                    [
                        [0.85985939, 0.0, 0.38410868],
                        [0.14014061, -0.0, 0.61589132],
                        [0.64014061, 0.0, 0.88410868],
                        [0.35985939, -0.0, 0.11589132],
                    ],
                ),
                0,
                parameters={"run_type": "R2SCAN"},
                entry_id="r2scan-3",
            ),
        ]
        with pytest.raises(ValueError, match="Unique entry_ids are required"):
            mixing_scheme_no_compat.process_entries(gga_entries + scan_entries)

    def test_clean(self, mixing_scheme_no_compat):
        # make sure the clean=True arg to process_entries works
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
                correction=-20,
                entry_id="entry-1",
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]),
                0,
                parameters={"run_type": "GGA"},
                correction=-20,
                entry_id="entry-2",
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "GGA"},
                correction=-20,
                entry_id="entry-3",
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
                entry_id="entry-4",
            ),
        ]

        mixing_scheme_no_compat.process_entries(entries, clean=False)
        for e in entries:
            assert e.correction == -20

        mixing_scheme_no_compat.process_entries(entries, clean=True)
        for e in entries:
            assert e.correction == 0

    def test_no_run_type(self, mixing_scheme_no_compat):
        """
        If one of the entries doesn't have a run_type attribute, we should get a warning
        from process_entries and a CompatibilityError from get_adjustments
        """
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]), 0, parameters={"run_type": "R2SCAN"}, entry_id="r2scan-1"
            ),
            ComputedStructureEntry(Structure(lattice, ["Br"], [[0, 0, 0]]), 0, parameters={}),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "R2SCAN"},
                entry_id="r2scan-2",
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
                entry_id="r2scan-3",
            ),
        ]

        with pytest.warns(UserWarning, match="missing parameters.run_type"):
            mixing_scheme_no_compat.process_entries(entries)

    def test_incompatible_run_type(self, mixing_scheme_no_compat):
        """
        If entry.parameters.run_type is not "GGA", "GGA+U", or "R2SCAN", we should get a
        warning from process_entries and a CompatibilityError from get_adjustments
        """
        lattice = Lattice.from_parameters(a=1, b=1, c=1, alpha=90, beta=90, gamma=60)
        entries = [
            ComputedStructureEntry(
                Structure(lattice, ["Sn"], [[0, 0, 0]]), 0, parameters={"run_type": "GGA"}, entry_id="gga-1"
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]), 0, parameters={"run_type": "GGA"}, entry_id="gga-2"
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Br"], [[0, 0, 0]]), 0, parameters={"run_type": "LDA"}, entry_id="lda-1"
            ),
            ComputedStructureEntry(
                Structure(lattice, ["Sn", "Br", "Br"], [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]),
                0,
                parameters={"run_type": "GGA"},
                entry_id="gga-3",
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
                entry_id="gga-4",
            ),
        ]

        with pytest.warns(UserWarning, match="Invalid run_type LDA"):
            assert len(mixing_scheme_no_compat.process_entries(entries)) == 4

        state_data = mixing_scheme_no_compat.get_mixing_state_data(entries)
        with pytest.raises(CompatibilityError, match="Invalid run_type LDA"):
            assert (
                mixing_scheme_no_compat.get_adjustments(
                    [e for e in entries if e.parameters["run_type"] == "LDA"][0],
                    state_data,
                )
                == []
            )

    def test_no_single_entry(self, mixing_scheme_no_compat):
        """
        Raise CompatibilityError if process_entries is called on a single Entry
        without any state_data.
        """
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

    def test_no_foreign_entries(self, mixing_scheme_no_compat, ms_complete):
        """
        If process_entries or get_adjustments is called with a populated mixing_state_data
        kwarg and one or more of the entry_ids is not present in the mixing_state_data,
        raise CompatbilityError
        """
        foreign_entry = ComputedStructureEntry(
            Structure(
                lattice3,
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
            entry_id="r2scan-8",
        )

        with pytest.raises(CompatibilityError, match="not found in the mixing state"):
            mixing_scheme_no_compat.get_adjustments(foreign_entry, ms_complete.state_data)

        # process_entries should discard all GGA entries and return all R2SCAN
        entries = mixing_scheme_no_compat.process_entries(
            ms_complete.all_entries + [foreign_entry], mixing_state_data=ms_complete.state_data
        )
        assert len(entries) == 7
        for e in entries:
            assert e.correction == 0
            assert e.parameters["run_type"] == "R2SCAN"

    def test_fuzzy_matching(self, ms_complete):
        """
        Test fuzzy diatomic matching. If fuzzy matching is disabled, then entry r2scan-3 will not
        match GGA ground state gga-3, preventing the mixing scheme from using the R2SCAN energy
        scale.

        In this situation, the mixing scheme should adjust the energies of all R2SCAN entries that
        match GGA ones to be identical to the GGA's, and discard the corresponding GGA entries. Entry
        r2scan-3 should be discarded because there is no reference ground state energy for it.
        """
        # fmt: off
        row_list = [
            ["Br",    64,  4,  True, "gga-3",       None, "GGA",     None,      0.,  np.nan,   0., -1.],
            ["Br",   191,  1, False, "gga-2", "r2scan-2", "GGA", "R2SCAN",      1.,     -1.,   0., -1.],
            ["Br",    64,  4, False,    None, "r2scan-3",  None, "R2SCAN",  np.nan,      0.,   0., -1.],
            ["Sn",   191,  1,  True, "gga-1", "r2scan-1", "GGA", "R2SCAN",      0.,     -1.,   0., -1.],
            ["SnBr2", 65,  3,  True, "gga-4", "r2scan-4", "GGA", "R2SCAN",     -6.,     -7.,  -6., -8.],
            ["SnBr2",  2, 12, False, "gga-5", "r2scan-5", "GGA", "R2SCAN",     -5.,     -8.,  -6., -8.],
            ["SnBr2", 71,  3, False, "gga-6", "r2scan-6", "GGA", "R2SCAN",     -4.,     -6.,  -6., -8.],
            ["SnBr4",  8,  5, False, "gga-7", "r2scan-7", "GGA", "R2SCAN",     -3.,     -6., -3.6, -6.],
        ]
        # fmt: on
        mixing_state = pd.DataFrame(row_list, columns=columns)

        compat = MaterialsProjectDFTMixingScheme(compat_1=None, fuzzy_matching=False)
        state_data = compat.get_mixing_state_data(ms_complete.all_entries)
        pd.testing.assert_frame_equal(state_data, mixing_state)

        for e in ms_complete.all_entries:
            if e.entry_id in ["gga-1", "gga-4"]:
                with pytest.raises(CompatibilityError, match="ground state"):
                    compat.get_adjustments(e, mixing_state)
            elif e.entry_id in ["r2scan-3"]:
                with pytest.raises(CompatibilityError, match="and no R2SCAN ground state"):
                    compat.get_adjustments(e, mixing_state)
            elif e.entry_id in ["gga-2", "gga-5", "gga-6", "gga-7"]:
                with pytest.raises(CompatibilityError, match="there is a matching R2SCAN"):
                    compat.get_adjustments(e, mixing_state)
            elif e.parameters["run_type"] == "GGA":
                assert not compat.get_adjustments(e, mixing_state)
            else:
                assert compat.get_adjustments(e, mixing_state)

        # process_entries should discard all GGA entries and return all R2SCAN
        entries = compat.process_entries(ms_complete.all_entries)
        assert len(entries) == 7
        for e in entries:
            if e.parameters["run_type"] == "GGA":
                assert e.correction == 0
            elif e.entry_id in ["r2scan-2", "r2scan-7"]:
                assert "Replace R2SCAN energy with GGA" in e.energy_adjustments[0].description
            else:
                assert "onto the GGA(+U) hull" in e.energy_adjustments[0].description

    def test_same_run_type(self):
        """
        Test behavior when run_type_1 and run_type_2 are the same
        or overlap
        """
        with pytest.raises(ValueError, match="the same run_type GGA"):
            MaterialsProjectDFTMixingScheme(run_type_1="GGA", run_type_2="GGA")

    def test_alternate_run_types(self):
        """
        Test behavior with alternate run_types
        """
        compat = MaterialsProjectDFTMixingScheme(run_type_1="LDA", run_type_2="GGA(+U)")
        assert compat.valid_rtypes_1 == ["LDA"]
        assert compat.valid_rtypes_2 == ["GGA", "GGA+U"]

    def test_compat_args(self, ms_complete):
        """
        Test the behavior of compat1 and compat2 kwargs
        The DummyCompatibility class defined in this test file should lower the
        energies of all entries by 10 eV/atom
        """
        compat = MaterialsProjectDFTMixingScheme(compat_1=DummyCompatibility(), compat_2=DummyCompatibility())
        rt1_entries, rt2_entries = compat._filter_and_sort_entries(ms_complete.all_entries)
        state_data = compat.get_mixing_state_data(rt1_entries + rt2_entries)
        assert max(state_data["hull_energy_1"]) <= -10  # highest hull energy is -10 eV/atom with the correction
        assert all(state_data["hull_energy_2"] <= -11)  # highest hull energy is -11 eV/atom with the correction

        compat.process_entries(ms_complete.all_entries)
        for e in ms_complete.all_entries:
            assert e.energy_adjustments[0].value == -10 * e.composition.num_atoms

        for e in ms_complete.gga_entries:
            if e.entry_id in ["gga-1", "gga-3", "gga-4"]:
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    compat.get_adjustments(e, state_data)
            else:
                with pytest.raises(CompatibilityError, match="there is a matching R2SCAN"):
                    assert not compat.get_adjustments(e, state_data)

        for e in ms_complete.scan_entries:
            assert not compat.get_adjustments(e, state_data)

    def test_no_mixing_data(self, ms_complete):
        """
        Test the behavior of get_adjustments when mixing_state_data is None
        """
        compat = MaterialsProjectDFTMixingScheme()

        with pytest.raises(CompatibilityError, match="DataFrame is None."):
            compat.get_adjustments(ms_complete.all_entries[0])

    def test_multiple_matching_structures(self, mixing_scheme_no_compat, ms_complete_duplicate_structs):
        """
        Test behavior when the entries contain many structures that match to
        the same material. For this test, entry gga-4 (SnBr2 ground state) and
        its matching r2scan-4 are each duplicated into new entries gga-10 and
        r2scan-10, respectively.

        In this situation, the mixing scheme should only keep the lowest energy entry
        within each run_type, so gga-10 and r2scan-10 should be discarded, and the
        results should otherwise be identical to ms_scan_complete.
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_complete_duplicate_structs.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_complete_duplicate_structs.state_data)

        for e in ms_complete_duplicate_structs.scan_entries:
            if e.entry_id == "r2scan-10":
                with pytest.raises(CompatibilityError, match="not found in the mixing state"):
                    mixing_scheme_no_compat.get_adjustments(e, state_data)
                continue
            assert mixing_scheme_no_compat.get_adjustments(e, state_data) == []

        for e in ms_complete_duplicate_structs.gga_entries:
            if e.entry_id == "gga-10":
                with pytest.raises(CompatibilityError, match="not found in the mixing state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_complete_duplicate_structs.state_data)
                continue
            elif e.entry_id in ["gga-1", "gga-3", "gga-4"]:
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_complete_duplicate_structs.state_data)
            else:
                with pytest.raises(CompatibilityError, match="there is a matching R2SCAN"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_complete_duplicate_structs.state_data)

        # process_entries should discard all GGA entries and return all R2SCAN
        entries = mixing_scheme_no_compat.process_entries(ms_complete_duplicate_structs.all_entries)
        assert len(entries) == 7
        for e in entries:
            assert e.correction == 0
            assert e.parameters["run_type"] == "R2SCAN"

    def test_alternate_structure_matcher(self, ms_complete):
        """
        Test alternate structure matcher kwargs. By setting scale to False, entries
        gga-6 and r2scan-6 will no longer match and will be listed as separate rows
        in the mixing state DataFrame.
        """
        sm = StructureMatcher(scale=False)
        compat = MaterialsProjectDFTMixingScheme(compat_1=None, structure_matcher=sm)
        state_data = compat.get_mixing_state_data(ms_complete.all_entries)
        assert isinstance(state_data, pd.DataFrame), "get_mixing_state_data failed to generate a DataFrame."
        assert len(state_data) == 8
        assert sum(state_data["run_type_1"] == "GGA") == len(state_data) - 1
        assert sum(state_data["run_type_2"] == "R2SCAN") == len(state_data) - 1
        assert sum(state_data["is_stable_1"]) == 3
        assert sum(state_data["energy_1"].isna()) == 1
        assert sum(state_data["energy_2"].isna()) == 1
        assert all(state_data["hull_energy_1"].notna())
        assert all(state_data["hull_energy_2"].notna())

        for e in ms_complete.scan_entries:
            assert not compat.get_adjustments(e, state_data)

        for e in ms_complete.gga_entries:
            if e.entry_id == "gga-6":
                assert compat.get_adjustments(e, state_data)[0].value == -6
                continue

            elif e.entry_id in ["gga-1", "gga-3", "gga-4"]:
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    compat.get_adjustments(e, ms_complete.state_data)
            else:
                with pytest.raises(CompatibilityError, match="there is a matching R2SCAN"):
                    compat.get_adjustments(e, ms_complete.state_data)

        # process_entries should discard all GGA entries except gga-6 and return all R2SCAN
        # entries unmodified. gga-6 should be corrected to the R2SCAN hull
        entries = compat.process_entries(ms_complete.all_entries)
        assert len(entries) == 8


class TestMaterialsProjectDFTMixingSchemeStates:
    """
    Test the behavior of the mixing scheme under different mixing states (i.e., different
    combinations of GGA and R2SCAN entries)
    """

    def test_state_complete_entries(self, mixing_scheme_no_compat, ms_complete):
        """
        Mixing state in which every material is present in both GGA and R2SCAN

        In this state, the mixing scheme should return only the R2SCAN entries, unmodified
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_complete.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_complete.state_data)

        for e in ms_complete.scan_entries:
            assert mixing_scheme_no_compat.get_adjustments(e, ms_complete.state_data) == []

        for e in ms_complete.gga_entries:
            if e.entry_id in ["gga-1", "gga-3", "gga-4"]:
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_complete.state_data)
            else:
                with pytest.raises(CompatibilityError, match="there is a matching R2SCAN"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_complete.state_data)

        # process_entries should discard all GGA entries and return all R2SCAN
        entries = mixing_scheme_no_compat.process_entries(ms_complete.all_entries)
        assert len(entries) == 7
        for e in entries:
            assert e.correction == 0
            assert e.parameters["run_type"] == "R2SCAN"

    def test_state_gga_only(self, mixing_scheme_no_compat, ms_gga_only):
        """
        Mixing state in which we only have GGA entries, forming a complete PhaseDiagram

        In this state, the mixing scheme should not do anything
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_gga_only.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_gga_only.state_data)

        for e in ms_gga_only.all_entries:
            assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_only.state_data) == []
        entries = mixing_scheme_no_compat.process_entries(ms_gga_only.all_entries)
        assert len(entries) == 7
        for e in entries:
            assert e.correction == 0
            assert e.parameters["run_type"] == "GGA"

    def test_state_scan_only(self, mixing_scheme_no_compat, ms_scan_only):
        """
        Mixing state in which we only have R2SCAN entries, forming a complete PhaseDiagram

        In this case, the mixing scheme should not do anything
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_scan_only.all_entries, verbose=True)
        pd.testing.assert_frame_equal(state_data, ms_scan_only.state_data)

        for e in ms_scan_only.all_entries:
            assert mixing_scheme_no_compat.get_adjustments(e, ms_scan_only.state_data) == []

        entries = mixing_scheme_no_compat.process_entries(ms_scan_only.all_entries)
        assert len(entries) == 7
        for e in entries:
            assert e.correction == 0
            assert e.parameters["run_type"] == "R2SCAN"

    def test_state_gga_1_scan(self, mixing_scheme_no_compat, ms_gga_1_scan):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and 1 R2SCAN entry
        The R2SCAN entry chosen is the GGA ground state for SnBr2 (r2scan-4)

        In this state, the mixing scheme should adjust the entry of r2scan-4 to
        match the GGA energy and discard entry gga-4.
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_gga_1_scan.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_gga_1_scan.state_data)

        for e in ms_gga_1_scan.gga_entries:
            if e.entry_id == "gga-4":
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_gga_1_scan.state_data)
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_1_scan.state_data) == []

        for e in ms_gga_1_scan.scan_entries:
            # gga-4 energy is -6 eV/atom, r2scan-4 energy is -7 eV/atom. There are 3 atoms.
            assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_1_scan.state_data)[0].value == 3

        entries = mixing_scheme_no_compat.process_entries(ms_gga_1_scan.all_entries)
        assert len(entries) == 7
        for e in entries:
            if "4" in e.entry_id:
                assert e.correction == 3
                assert e.parameters["run_type"] == "R2SCAN"
            else:
                assert e.correction == 0, f"{e.entry_id}"
                assert e.parameters["run_type"] == "GGA"

    def test_state_gga_1_scan_plus_novel(self, mixing_scheme_no_compat, ms_gga_1_scan_novel):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and 1 R2SCAN entry
        at a composition not in the GGA phase diagram.

        In this state, the mixing scheme should discard the R2SCAN entry
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_gga_1_scan_novel.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_gga_1_scan_novel.state_data)

        for e in ms_gga_1_scan_novel.gga_entries:
            if e.entry_id == "gga-4":
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_1_scan_novel.state_data) == []

        for e in ms_gga_1_scan_novel.scan_entries:
            with pytest.raises(CompatibilityError, match="no R2SCAN ground states at this composition"):
                mixing_scheme_no_compat.get_adjustments(e, ms_gga_1_scan_novel.state_data)

        entries = mixing_scheme_no_compat.process_entries(ms_gga_1_scan_novel.all_entries)
        assert len(entries) == 7

    def test_state_gga_2_scan_same(self, mixing_scheme_no_compat, ms_gga_2_scan_same):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and 2 R2SCAN entries
        at a single composition, one of which is the GGA ground state.

        In this state, the mixing scheme should correct the energy of unstable polymorph
        r2scan-6 to maintain the same e_above_hull (r2scan-6 minus r2scan-4 which is the ground
        state). Entry r2scan-4 (the GGA ground state) should be corrected to the GGA energy.
        Entries gga-4 and gga-6 should be discarded.
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_gga_2_scan_same.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_gga_2_scan_same.state_data)

        for e in ms_gga_2_scan_same.scan_entries:
            if e.entry_id == "r2scan-6":
                # gga-4 energy is -18 eV or -6 eV/atom. r2scan-6 is 1 eV/atom above r2scan-4 (ground state),
                # so r2scan-6 should be adjusted to -5 eV/atom or -15 eV, which is 3 eV higher than its
                # original energy of -18 eV.
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_same.state_data)[0].value == 3
            elif e.entry_id == "r2scan-4":
                # r2scan-4 energy is -7 eV/atom. Needs to be adjusted to -6 eV/atom (3 atoms)
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_same.state_data)[0].value == 3
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_same.state_data) == []

        for e in ms_gga_2_scan_same.gga_entries:
            if e.entry_id == "gga-4":
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_same.state_data)
            elif e.entry_id == "gga-6":
                with pytest.raises(CompatibilityError, match="there is a matching R2SCAN"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_same.state_data)
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_same.state_data) == []

        entries = mixing_scheme_no_compat.process_entries(ms_gga_2_scan_same.all_entries)
        assert len(entries) == 7
        for e in entries:
            if e.entry_id in ["r2scan-4", "r2scan-6"]:
                assert e.correction == 3
                assert e.parameters["run_type"] == "R2SCAN"
            elif e.entry_id == "gga-4":
                raise AssertionError("Entry gga-4 should have been discarded")
            elif e.entry_id == "gga-6":
                raise AssertionError("Entry gga-6 should have been discarded")
            else:
                assert e.correction == 0, f"{e.entry_id}"
                assert e.parameters["run_type"] == "GGA"

    def test_state_gga_2_scan_diff_match(self, mixing_scheme_no_compat, ms_gga_2_scan_diff_match):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and 2 R2SCAN entries
        at different compositions, where both R2SCAN materials match GGA materials but
        only one matches a GGA ground state.

        In this state, the energies of both R2SCAN entries should be set equal to the
        corresponding GGA energies, and the GGA entries discarded.
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_gga_2_scan_diff_match.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_gga_2_scan_diff_match.state_data)

        for e in ms_gga_2_scan_diff_match.scan_entries:
            if e.entry_id == "r2scan-7":
                # r2scan-7 energy is -6 eV/atom, needs to go to -3 eV/atom (5 atoms)
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_match.state_data)[0].value == 15
            elif e.entry_id == "r2scan-4":
                # r2scan-4 energy is -7 eV/atom. Needs to be adjusted to -6 eV/atom (3 atoms)
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_match.state_data)[0].value == 3
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_match.state_data) == []

        for e in ms_gga_2_scan_diff_match.gga_entries:
            if e.entry_id == "gga-4":
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_match.state_data)
            elif e.entry_id == "gga-7":
                with pytest.raises(CompatibilityError, match="there is a matching R2SCAN"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_match.state_data)
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_match.state_data) == []

        entries = mixing_scheme_no_compat.process_entries(ms_gga_2_scan_diff_match.all_entries)
        assert len(entries) == 7
        for e in entries:
            if e.entry_id == "r2scan-4":
                assert e.correction == 3
            elif e.entry_id == "r2scan-7":
                assert e.correction == 15
            elif e.entry_id in ["gga-4"]:
                raise AssertionError(f"Entry {e.entry_id} should have been discarded")
            else:
                assert e.correction == 0, f"{e.entry_id}"
                assert e.parameters["run_type"] == "GGA"

    def test_state_gga_2_scan_diff_nomatch(self, mixing_scheme_no_compat, ms_gga_2_scan_diff_no_match):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and 2 R2SCAN entries
        at different compositions, where one of the R2SCAN materials does not match
        any GGA material.

        In this state, the energy of the matching R2SCAN entry should be adjusted
        to the GGA value, the corresponding GGA entry should be discarded, and the
        novel R2SCAN material that doesn't match anything should be discarded
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_gga_2_scan_diff_no_match.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_gga_2_scan_diff_no_match.state_data)

        for e in ms_gga_2_scan_diff_no_match.scan_entries:
            if e.entry_id == "r2scan-8":
                # there is no matching GGA structure for r2scan-8, so there's no way
                # to adjust its energy onto the GGA hull.
                with pytest.raises(CompatibilityError, match="entry and no R2SCAN ground state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_no_match.state_data)
            elif e.entry_id == "r2scan-4":
                # r2scan-4 energy is -7 eV/atom. Needs to be adjusted to -6 eV/atom (3 atoms)
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_no_match.state_data)[0].value == 3
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_no_match.state_data) == []

        for e in ms_gga_2_scan_diff_no_match.gga_entries:
            if e.entry_id == "gga-4":
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_no_match.state_data)
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_gga_2_scan_diff_no_match.state_data) == []

        entries = mixing_scheme_no_compat.process_entries(ms_gga_2_scan_diff_no_match.all_entries)
        assert len(entries) == 7
        for e in entries:
            if e.entry_id == "r2scan-4":
                assert e.correction == 3
                assert e.parameters["run_type"] == "R2SCAN"
            elif e.entry_id in ["gga-4", "r2scan-8"]:
                raise AssertionError(f"Entry {e.entry_id} should have been discarded")
            else:
                assert e.correction == 0, f"{e.entry_id}"
                assert e.parameters["run_type"] == "GGA"

    def test_state_incomplete_gga_all_scan(self, mixing_scheme_no_compat, ms_incomplete_gga_all_scan):
        """
        Mixing state in which we have an incomplete GGA PhaseDiagram and all entries
        present in R2SCAN.

        This case should fail, because a complete run_type_1 phase diagram is required by the
        mixing scheme
        """
        # Test behavior when GGA entries don't form a complete phase diagram
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_incomplete_gga_all_scan.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_incomplete_gga_all_scan.state_data)

        for e in ms_incomplete_gga_all_scan.all_entries:
            with pytest.raises(CompatibilityError, match="do not form a complete PhaseDiagram"):
                mixing_scheme_no_compat.get_adjustments(e, ms_incomplete_gga_all_scan.state_data)

        # process_entries should discard all entries and issue a warning
        with pytest.warns(UserWarning, match="do not form a complete PhaseDiagram"):
            entries = mixing_scheme_no_compat.process_entries(ms_incomplete_gga_all_scan.all_entries)
        assert len(entries) == 0

    def test_state_all_gga_scan_gs(self, mixing_scheme_no_compat, ms_all_gga_scan_gs):
        """
        Mixing state in which we have a complete GGA PhaseDiagram and all GGA
        ground states present as R2SCAN entries.

        In this situation, we should build the hull using R2SCAN energies, discard
        the GGA ground state materials, and correct the remaining GGA energies onto
        the R2SCAN hull such that they maintain their e_above_hull.
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_all_gga_scan_gs.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_all_gga_scan_gs.state_data)

        for e in ms_all_gga_scan_gs.scan_entries:
            assert mixing_scheme_no_compat.get_adjustments(e, ms_all_gga_scan_gs.state_data) == []

        for e in ms_all_gga_scan_gs.gga_entries:
            if e.entry_id in ["gga-1", "gga-3", "gga-4"]:
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_all_gga_scan_gs.state_data)
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_all_gga_scan_gs.state_data) != []

        # store the energy above hull of each gga entry
        pd_gga = PhaseDiagram(ms_all_gga_scan_gs.gga_entries)
        gga_hull_e = {}
        for e in ms_all_gga_scan_gs.gga_entries:
            gga_hull_e[e.entry_id] = pd_gga.get_e_above_hull(e)

        # process_entries should discard 3 GGA ground state entries
        entries = mixing_scheme_no_compat.process_entries(ms_all_gga_scan_gs.all_entries)
        assert len(entries) == len(ms_all_gga_scan_gs.all_entries) - 3
        pd_mixed = PhaseDiagram(entries)

        for e in entries:
            if e.parameters["run_type"] == "GGA":
                assert "onto the R2SCAN hull" in e.energy_adjustments[0].description
                assert np.allclose(pd_mixed.get_e_above_hull(e), gga_hull_e[e.entry_id])
            else:
                assert e.correction == 0

    def test_state_novel_scan_comp(self, mixing_scheme_no_compat, ms_all_gga_scan_gs_plus_novel):
        """
        Mixing state in which we have all GGA ground states in R2SCAN and then we try to
        process a R2SCAN entry at a composition that is not in the GGA PhaseDiagram.

        In this case, the mixing scheme should preserve the energy of the novel R2SCAN
        entry and discard the 3 GGA ground states
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_all_gga_scan_gs_plus_novel.all_entries)
        pd.testing.assert_frame_equal(state_data, ms_all_gga_scan_gs_plus_novel.state_data)

        for e in ms_all_gga_scan_gs_plus_novel.scan_entries:
            assert mixing_scheme_no_compat.get_adjustments(e, ms_all_gga_scan_gs_plus_novel.state_data) == []

        for e in ms_all_gga_scan_gs_plus_novel.gga_entries:
            if e.entry_id in ["gga-1", "gga-3", "gga-4"]:
                with pytest.raises(CompatibilityError, match="because it is a GGA\\(\\+U\\) ground state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_all_gga_scan_gs_plus_novel.state_data)
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_all_gga_scan_gs_plus_novel.state_data) != []

        # store the energy above hull of each gga entry
        pd_gga = PhaseDiagram(ms_all_gga_scan_gs_plus_novel.gga_entries)
        gga_hull_e = {}
        for e in ms_all_gga_scan_gs_plus_novel.gga_entries:
            gga_hull_e[e.entry_id] = pd_gga.get_e_above_hull(e)

        # process_entries should discard 3 GGA ground state entries
        entries = mixing_scheme_no_compat.process_entries(ms_all_gga_scan_gs_plus_novel.all_entries)
        assert len(entries) == len(ms_all_gga_scan_gs_plus_novel.all_entries) - 3
        pd_mixed = PhaseDiagram(entries)

        for e in entries:
            if e.parameters["run_type"] == "GGA":
                assert "onto the R2SCAN hull" in e.energy_adjustments[0].description
                assert np.allclose(pd_mixed.get_e_above_hull(e), gga_hull_e[e.entry_id])
            else:
                assert e.correction == 0

    def test_state_energy_modified(self, mixing_scheme_no_compat, ms_complete):
        """
        Mixing state in which we try to process an entry whose energy has been
        changed and is no longer consistent with the mixing state data.

        This is a corner case that should never occur, because if the entry_id
        of this entry is not already in the state_data, a CompatibilityError
        should be raised. If the Entry was passed to get_mixing_state_data, then
        its energy should already be there. Hence, this is testing a case
        in which either 1) get_mixing_state_data fails to work properly or 2)
        the energy of the Entry is somehow modified in between calling
        get_mixing_state_data and process_entries. Such a situation could
        potentially arise in e.g. the build pipeline if one is calling
        process_entries with a separately-calculated state_data DataFrame.
        """
        state_data = mixing_scheme_no_compat.get_mixing_state_data(ms_complete.all_entries)
        # lower the energy of the SnBr2 ground state
        e = [e for e in ms_complete.gga_entries if e.entry_id == "gga-4"][0]
        d_compat = DummyCompatibility()
        d_compat.process_entries(e)

        with pytest.raises(CompatibilityError, match="energy has been modified"):
            mixing_scheme_no_compat.get_adjustments(e, state_data)

    def test_chemsys_mismatch(self, mixing_scheme_no_compat, ms_scan_chemsys_superset):
        """
        Test what happens if the entries aren't in the same chemsys. run_type_2 entries
        that are outside the run_type_1 chemsys should be discarded.
        """
        rt1_entries, rt2_entries = mixing_scheme_no_compat._filter_and_sort_entries(
            ms_scan_chemsys_superset.all_entries
        )
        state_data = mixing_scheme_no_compat.get_mixing_state_data(rt1_entries + rt2_entries)
        pd.testing.assert_frame_equal(state_data, ms_scan_chemsys_superset.state_data)

        for e in ms_scan_chemsys_superset.scan_entries:
            if e.entry_id == "r2scan-9":
                with pytest.raises(CompatibilityError, match="not found in the mixing state"):
                    mixing_scheme_no_compat.get_adjustments(e, ms_scan_chemsys_superset.state_data)
            else:
                assert mixing_scheme_no_compat.get_adjustments(e, ms_scan_chemsys_superset.state_data) == []

        # process_entries should discard all GGA entries and return all R2SCAN
        with pytest.warns(UserWarning, match="is larger than GGA\\(\\+U\\) entries chemical system"):
            entries = mixing_scheme_no_compat.process_entries(ms_scan_chemsys_superset.all_entries)
            assert len(entries) == 7
            for e in entries:
                assert e.correction == 0
                assert e.parameters["run_type"] == "R2SCAN"
