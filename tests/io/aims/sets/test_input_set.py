from __future__ import annotations

import copy
import json

import pytest

from pymatgen.core import Structure
from pymatgen.io.aims.sets import AimsInputSet
from pymatgen.util.testing import TEST_FILES_DIR

IN_FILE_DIR = TEST_FILES_DIR / "io/aims/input_files"
SPECIES_DIR = TEST_FILES_DIR / "io/aims/species_directory"


control_in_str = """
#===============================================================================
# FHI-aims geometry file: ./geometry.in
# File generated from pymatgen
# Fri Jan  3 16:35:37 2025
#===============================================================================
xc                                                pbe
k_grid                                            2 2 2
compute_forces                                    .true.
#===============================================================================

################################################################################
#
#  FHI-aims code project
#  VB, Fritz-Haber Institut, 2009
#
#  Suggested "light" defaults for Si atom (to be pasted into control.in file)
#  Be sure to double-check any results obtained with these settings for post-processing,
#  e.g., with the "tight" defaults and larger basis sets.
#
#  2020/09/08 Added f function to "light" after reinspection of Delta test outcomes.
#             This was done for all of Al-Cl and is a tricky decision since it makes
#             "light" calculations measurably more expensive for these elements.
#             Nevertheless, outcomes for P, S, Cl (and to some extent, Si) appear
#             to justify this choice.
#
################################################################################
  species               Si
#     global species definitions
  nucleus               14.0
  mass                  28.0855
#
  l_hartree             4
#
  cut_pot               3.5 1.5 1.0
  basis_dep_cutoff      0.0001
#
  radial_base    42   5.0000
  radial_multiplier        1
  angular_grids  specified
  division       0.5866   50
  division       0.9616  110
  division       1.2249  194
  division       1.3795  302
  outer_grid     302
################################################################################
#
#  Definition of "minimal" basis
#
################################################################################
#     valence basis states
    valence      3  s   2.
    valence      3  p   2.
#     ion occupancy
    ion_occ      3  s   1.
    ion_occ      3  p   1.
################################################################################
#
#  Suggested additional basis functions. For production calculations,
#  uncomment them one after another (the most important basis functions are
#  listed first).
#
#  These were set using pyfhiaims, original may files have additional comments
#
################################################################################
#  First tier
             hydro      3  d  4.2
             hydro      2  p  1.4
             hydro      4  f  6.2
             ionic      3  s  auto
#  Second tier
#            hydro      3  d  9.0
#            hydro      5  g  9.4
#            hydro      4  p  4.0
#            hydro      1  s  0.65
#  Third tier
#            ionic      3  d  auto
#            hydro      3  s  2.6
#            hydro      4  f  8.4
#            hydro      3  d  3.4
#            hydro      3  p  7.8
#  Fourth tier
#            hydro      2  p  1.6
#            hydro      5  g  10.8
#            hydro      5  f  11.2
#            hydro      3  d  1.0
#            hydro      4  s  4.5
#  Further basis functions
#            hydro      4  d  6.6
#            hydro      5  g  16.4
#            hydro      4  d  9.0
"""

control_in_str_rel = """
#===============================================================================
# FHI-aims geometry file: ./geometry.in
# File generated from pymatgen
# Fri Jan  3 16:35:37 2025
#===============================================================================
xc                                                pbe
k_grid                                            2 2 2
relax_geometry                                    trm 1e-3
compute_forces                                    .true.
#===============================================================================

################################################################################
#
#  FHI-aims code project
#  VB, Fritz-Haber Institut, 2009
#
#  Suggested "light" defaults for Si atom (to be pasted into control.in file)
#  Be sure to double-check any results obtained with these settings for post-processing,
#  e.g., with the "tight" defaults and larger basis sets.
#
#  2020/09/08 Added f function to "light" after reinspection of Delta test outcomes.
#             This was done for all of Al-Cl and is a tricky decision since it makes
#             "light" calculations measurably more expensive for these elements.
#             Nevertheless, outcomes for P, S, Cl (and to some extent, Si) appear
#             to justify this choice.
#
################################################################################
  species               Si
#     global species definitions
  nucleus               14.0
  mass                  28.0855
#
  l_hartree             4
#
  cut_pot               3.5 1.5 1.0
  basis_dep_cutoff      0.0001
#
  radial_base    42   5.0000
  radial_multiplier        1
  angular_grids  specified
  division       0.5866   50
  division       0.9616  110
  division       1.2249  194
  division       1.3795  302
  outer_grid     302
################################################################################
#
#  Definition of "minimal" basis
#
################################################################################
#     valence basis states
    valence      3  s   2.
    valence      3  p   2.
#     ion occupancy
    ion_occ      3  s   1.
    ion_occ      3  p   1.
################################################################################
#
#  Suggested additional basis functions. For production calculations,
#  uncomment them one after another (the most important basis functions are
#  listed first).
#
#  These were set using pyfhiaims, original may files have additional comments
#
################################################################################
#  First tier
             hydro      3  d  4.2
             hydro      2  p  1.4
             hydro      4  f  6.2
             ionic      3  s  auto
#  Second tier
#            hydro      3  d  9.0
#            hydro      5  g  9.4
#            hydro      4  p  4.0
#            hydro      1  s  0.65
#  Third tier
#            ionic      3  d  auto
#            hydro      3  s  2.6
#            hydro      4  f  8.4
#            hydro      3  d  3.4
#            hydro      3  p  7.8
#  Fourth tier
#            hydro      2  p  1.6
#            hydro      5  g  10.8
#            hydro      5  f  11.2
#            hydro      3  d  1.0
#            hydro      4  s  4.5
#  Further basis functions
#            hydro      4  d  6.6
#            hydro      5  g  16.4
#            hydro      4  d  9.0
"""

geometry_in_str = """
#===============================================================================
# FHI-aims geometry file: geometry.in
# File generated from pymatgen
# Fri Jan  3 16:35:37 2025
#===============================================================================
lattice_vector 0.000000000000000e+00 2.715000000000000e+00 2.715000000000000e+00
lattice_vector 2.715000000000000e+00 0.000000000000000e+00 2.715000000000000e+00
lattice_vector 2.715000000000000e+00 2.715000000000000e+00 0.000000000000000e+00
atom         0.000000000000e+00   0.000000000000e+00   0.000000000000e+00 Si
atom         1.357500000000e+00   1.357500000000e+00   1.357500000000e+00 Si

"""

geometry_in_str_new = """
#===============================================================================
# FHI-aims geometry file: geometry.in
# File generated from pymatgen
# Fri Jan  3 16:35:37 2025
#===============================================================================
lattice_vector 0.000000000000000e+00 2.715000000000000e+00 2.715000000000000e+00
lattice_vector 2.715000000000000e+00 0.000000000000000e+00 2.715000000000000e+00
lattice_vector 2.715000000000000e+00 2.715000000000000e+00 0.000000000000000e+00
atom         0.000000000000e+00  -2.715000000000e-02  -2.715000000000e-02 Si
atom         1.357500000000e+00   1.357500000000e+00   1.357500000000e+00 Si
"""


def check_file(ref: str, test: str) -> bool:
    ref_lines = [line.strip() for line in ref.split("\n") if len(line.strip()) > 0 and line[0] != "#"]
    test_lines = [line.strip() for line in test.split("\n") if len(line.strip()) > 0 and line[0] != "#"]

    return ref_lines == test_lines


def test_input_set():
    Si = Structure(
        lattice=[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]],
        species=["Si", "Si"],
        coords=[[0.0] * 3, [0.25] * 3],
    )
    params_json = {
        "xc": "pbe",
        "species_dir": f"{SPECIES_DIR}/light",
        "k_grid": [2, 2, 2],
    }
    params_json_rel = {
        "xc": "pbe",
        "species_dir": f"{SPECIES_DIR}/light",
        "k_grid": [2, 2, 2],
        "relax_geometry": "trm 1e-3",
    }

    parameters = {
        "xc": "pbe",
        "species_dir": f"{SPECIES_DIR}/light",
        "k_grid": [2, 2, 2],
    }
    props = ("energy", "free_energy", "forces")

    in_set = AimsInputSet(parameters, Si, props)

    assert check_file(geometry_in_str, in_set.geometry_in)
    assert check_file(control_in_str, in_set.control_in)
    assert params_json == json.loads(in_set.params_json)

    in_set_copy = copy.deepcopy(in_set)
    assert check_file(geometry_in_str, in_set_copy.geometry_in)
    assert check_file(control_in_str, in_set_copy.control_in)
    assert params_json == json.loads(in_set_copy.params_json)

    in_set.set_parameters(**parameters, relax_geometry="trm 1e-3")
    assert check_file(control_in_str_rel, in_set.control_in)
    assert check_file(control_in_str, in_set_copy.control_in)

    assert params_json_rel == json.loads(in_set.params_json)
    assert params_json == json.loads(in_set_copy.params_json)

    in_set.remove_parameters(keys=["relax_geometry"])
    assert check_file(control_in_str, in_set.control_in)
    assert params_json == json.loads(in_set.params_json)

    in_set.remove_parameters(keys=["relax_geometry"], strict=False)
    assert check_file(control_in_str, in_set.control_in)
    assert params_json == json.loads(in_set.params_json)

    with pytest.raises(ValueError, match="key='relax_geometry' not in list"):
        in_set.remove_parameters(keys=["relax_geometry"], strict=True)

    new_struct = Structure(
        lattice=[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]],
        species=["Si", "Si"],
        coords=[[-0.01, 0, 0], [0.25, 0.25, 0.25]],
    )
    in_set.set_structure(new_struct)
    print(in_set.geometry_in)
    assert check_file(geometry_in_str_new, in_set.geometry_in)
    assert check_file(geometry_in_str, in_set_copy.geometry_in)
