from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from pymatgen.core import Structure
from pymatgen.io.aims.sets import AimsInputSet

control_in_str = """
#===============================================================================
# Created using the Atomic Simulation Environment (ASE)
#
# Thu Oct  5 12:27:49 2023
#
#===============================================================================
xc                                 pbe
k_grid                             2 2 2
compute_forces                     .true.
#===============================================================================

################################################################################
#
#  FHI-aims code project
#  VB, Fritz-Haber Institut, 2009
#
#  Suggested "light" defaults for Si atom (to be pasted into control.in file)
#  Be sure to double-check any results obtained with these settings for post-processing,
#  e.g. with the "tight" defaults and larger basis sets.
#
#  2020/09/08 Added f function to "light" after reinspection of Delta test outcomes.
#             This was done for all of Al-Cl and is a tricky decision since it makes
#             "light" calculations measurably more expensive for these elements.
#             Nevertheless, outcomes for P, S, Cl (and to some extent, Si) appear
#             to justify this choice.
#
################################################################################
  species        Si
#     global species definitions
    nucleus             14
    mass                28.0855
#
    l_hartree           4
#
    cut_pot             3.5          1.5  1.0
    basis_dep_cutoff    1e-4
#
    radial_base         42 5.0
    radial_multiplier   1
    angular_grids       specified
      division   0.5866   50
      division   0.9616  110
      division   1.2249  194
      division   1.3795  302
#      division   1.4810  434
#      division   1.5529  590
#      division   1.6284  770
#      division   1.7077  974
#      division   2.4068 1202
#      outer_grid   974
      outer_grid 302
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
#  Constructed for dimers: 1.75 A, 2.0 A, 2.25 A, 2.75 A, 3.75 A
#
################################################################################
#  "First tier" - improvements: -571.96 meV to -37.03 meV
     hydro 3 d 4.2
     hydro 2 p 1.4
     hydro 4 f 6.2
     ionic 3 s auto
#  "Second tier" - improvements: -16.76 meV to -3.03 meV
#     hydro 3 d 9
#     hydro 5 g 9.4
#     hydro 4 p 4
#     hydro 1 s 0.65
#  "Third tier" - improvements: -3.89 meV to -0.60 meV
#     ionic 3 d auto
#     hydro 3 s 2.6
#     hydro 4 f 8.4
#     hydro 3 d 3.4
#     hydro 3 p 7.8
#  "Fourth tier" - improvements: -0.33 meV to -0.11 meV
#     hydro 2 p 1.6
#     hydro 5 g 10.8
#     hydro 5 f 11.2
#     hydro 3 d 1
#     hydro 4 s 4.5
#  Further basis functions that fell out of the optimization - noise
#  level... < -0.08 meV
#     hydro 4 d 6.6
#     hydro 5 g 16.4
#     hydro 4 d 9
"""

control_in_str_rel = """
#===============================================================================
# Created using the Atomic Simulation Environment (ASE)
#
# Thu Oct  5 12:33:50 2023
#
#===============================================================================
xc                                 pbe
k_grid                             2 2 2
relax_geometry                     trm 1e-3
compute_forces                     .true.
#===============================================================================

################################################################################
#
#  FHI-aims code project
#  VB, Fritz-Haber Institut, 2009
#
#  Suggested "light" defaults for Si atom (to be pasted into control.in file)
#  Be sure to double-check any results obtained with these settings for post-processing,
#  e.g. with the "tight" defaults and larger basis sets.
#
#  2020/09/08 Added f function to "light" after reinspection of Delta test outcomes.
#             This was done for all of Al-Cl and is a tricky decision since it makes
#             "light" calculations measurably more expensive for these elements.
#             Nevertheless, outcomes for P, S, Cl (and to some extent, Si) appear
#             to justify this choice.
#
################################################################################
  species        Si
#     global species definitions
    nucleus             14
    mass                28.0855
#
    l_hartree           4
#
    cut_pot             3.5          1.5  1.0
    basis_dep_cutoff    1e-4
#
    radial_base         42 5.0
    radial_multiplier   1
    angular_grids       specified
      division   0.5866   50
      division   0.9616  110
      division   1.2249  194
      division   1.3795  302
#      division   1.4810  434
#      division   1.5529  590
#      division   1.6284  770
#      division   1.7077  974
#      division   2.4068 1202
#      outer_grid   974
      outer_grid 302
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
#  Constructed for dimers: 1.75 A, 2.0 A, 2.25 A, 2.75 A, 3.75 A
#
################################################################################
#  "First tier" - improvements: -571.96 meV to -37.03 meV
     hydro 3 d 4.2
     hydro 2 p 1.4
     hydro 4 f 6.2
     ionic 3 s auto
#  "Second tier" - improvements: -16.76 meV to -3.03 meV
#     hydro 3 d 9
#     hydro 5 g 9.4
#     hydro 4 p 4
#     hydro 1 s 0.65
#  "Third tier" - improvements: -3.89 meV to -0.60 meV
#     ionic 3 d auto
#     hydro 3 s 2.6
#     hydro 4 f 8.4
#     hydro 3 d 3.4
#     hydro 3 p 7.8
#  "Fourth tier" - improvements: -0.33 meV to -0.11 meV
#     hydro 2 p 1.6
#     hydro 5 g 10.8
#     hydro 5 f 11.2
#     hydro 3 d 1
#     hydro 4 s 4.5
#  Further basis functions that fell out of the optimization - noise
#  level... < -0.08 meV
#     hydro 4 d 6.6
#     hydro 5 g 16.4
#     hydro 4 d 9
"""

geometry_in_str = """
#===============================================================================
# Created using the Atomic Simulation Environment (ASE)
#
# Thu Oct  5 12:27:49 2023
#=======================================================
lattice_vector  0.000000000000e+00  2.715000000000e+00  2.715000000000e+00
lattice_vector  2.715000000000e+00  0.000000000000e+00  2.715000000000e+00
lattice_vector  2.715000000000e+00  2.715000000000e+00  0.000000000000e+00
atom  0.000000000000e+00  0.000000000000e+00  0.000000000000e+00 Si
atom  1.357500000000e+00  1.357500000000e+00  1.357500000000e+00 Si
"""
geometry_in_str_new = """
#===============================================================================
# Created using the Atomic Simulation Environment (ASE)
#
# Thu Oct  5 12:27:49 2023
#=======================================================
lattice_vector  0.000000000000e+00  2.715000000000e+00  2.715000000000e+00
lattice_vector  2.715000000000e+00  0.000000000000e+00  2.715000000000e+00
lattice_vector  2.715000000000e+00  2.715000000000e+00  0.000000000000e+00
atom  0.000000000000e+00 -2.715000000000e-02 -2.715000000000e-02 Si
atom  1.357500000000e+00  1.357500000000e+00  1.357500000000e+00 Si
"""

infile_dir = Path(__file__).parent / "input_files"


def check_file(ref: str, test: str) -> bool:
    ref_lines = [line.strip() for line in ref.split("\n") if len(line.strip()) > 0 and line[0] != "#"]
    test_lines = [line.strip() for line in test.split("\n") if len(line.strip()) > 0 and line[0] != "#"]

    return ref_lines == test_lines


def test_input_set():
    species_dir = infile_dir.parents[1] / "species_directory/"
    Si = Structure(
        lattice=[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]],
        species=["Si", "Si"],
        coords=[[0.0] * 3, [0.25] * 3],
    )
    params_json = {
        "xc": "pbe",
        "species_dir": f"{species_dir}/light",
        "k_grid": [2, 2, 2],
    }
    params_json_rel = {
        "xc": "pbe",
        "species_dir": f"{species_dir}/light",
        "k_grid": [2, 2, 2],
        "relax_geometry": "trm 1e-3",
    }

    parameters = {
        "xc": "pbe",
        "species_dir": f"{species_dir}/light",
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
    assert check_file(geometry_in_str_new, in_set.geometry_in)
    assert check_file(geometry_in_str, in_set_copy.geometry_in)
