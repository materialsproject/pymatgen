"""
This module is intended to match crystal structures against known crystallographic "prototype"
structures.

In this module, the AflowPrototypeMatcher uses the AFLOW LIBRARY OF CRYSTALLOGRAPHIC PROTOTYPES.
If using this particular class, please cite their publication appropriately:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
The AFLOW library of crystallographic prototypes: part 1.
Computational Materials Science, 136, S1-S828.
https://doi.org/10.1016/j.commatsci.2017.01.017

The module also contains functions for getting the protostructure labels from a
variety of symmetry detection libraries (spglib, moyopy, aflow-sym). The protostructure
label is defined as the canonicalized aflow label with the alphabetically sorted
chemical system appended - `aflow_sym_label:chemsys`.

The utilities for determining the protostructure label are upstreamed from the
`aviary` package (https://github.com/CompRhys/aviary). If using these functions,
please cite the following publications:

Goodall, R. E., Parackal, A. S., Faber, F. A., Armiento, R., & Lee, A. A. (2022).
Rapid discovery of stable materials by coordinate-free coarse graining.
Science advances, 8(30), eabn4117. https://doi.org/10.1126/sciadv.abn4117

Parackal, A. S., Goodall, R. E., Faber, F. A., & Armiento, R. (2024).
Identifying crystal structures beyond known prototypes from x-ray powder diffraction spectra.
Physical Review Materials, 8(10), 103801. https://doi.org/10.1103/PhysRevMaterials.8.103801
"""

from __future__ import annotations

import os
import re
import subprocess
from collections import Counter, defaultdict
from itertools import chain, groupby, permutations, product
from operator import itemgetter
from shutil import which
from string import ascii_uppercase, digits
from typing import TYPE_CHECKING

import orjson
from monty.fractions import gcd
from monty.serialization import loadfn

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Composition, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing import Literal

    from pymatgen.core.structure import Structure

try:
    from pyxtal import pyxtal

    has_pyxtal = True
except ImportError:
    pyxtal = None
    has_pyxtal = False

try:
    import moyopy
    from moyopy.interface import MoyoAdapter

    has_moyopy = True
except ImportError:
    moyopy = None
    MoyoAdapter = None
    has_moyopy = False

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
AFLOW_PROTOTYPE_LIBRARY = loadfn(f"{MODULE_DIR}/aflow_prototypes.json.gz")
WYCKOFF_MULTIPLICITY_DICT = loadfn(f"{MODULE_DIR}/wyckoff-position-multiplicities.json.gz")
WYCKOFF_POSITION_PARAM_DICT = loadfn(f"{MODULE_DIR}/wyckoff-position-params.json.gz")
WYCKOFF_POSITION_RELAB_DICT = loadfn(f"{MODULE_DIR}/wyckoff-position-relabelings.json.gz")

WYCKOFF_POSITION_RELAB_DICT = {
    spg_num: [{int(key): line for key, line in val.items()} for val in vals]
    for spg_num, vals in WYCKOFF_POSITION_RELAB_DICT.items()
}

CRYSTAL_FAMILY_SYMBOLS = {
    "triclinic": "a",
    "monoclinic": "m",
    "orthorhombic": "o",
    "tetragonal": "t",
    "trigonal": "h",
    "hexagonal": "h",
    "cubic": "c",
}

CRYSTAL_LATTICE_PARAMETERS_COUNTS = {
    "a": 6,
    "m": 4,
    "o": 3,
    "t": 2,
    "h": 2,
    "c": 1,
}

remove_digits = str.maketrans("", "", digits)

# Define regex patterns as constants
RE_WYCKOFF_NO_PREFIX = re.compile(r"((?<![0-9])[A-z])")
RE_ELEMENT_NO_SUFFIX = re.compile(r"([A-z](?![0-9]))")
RE_WYCKOFF = re.compile(r"(?<!\d)([a-zA-Z])")
RE_ANONYMOUS = re.compile(r"([A-Z])(?![0-9])")
RE_SUBST_ONE_PREFIX = r"1\g<1>"
RE_SUBST_ONE_SUFFIX = r"\g<1>1"


@due.dcite(
    Doi("10.1016/j.commatsci.2017.01.017"),
    description="The AFLOW library of crystallographic prototypes: part 1.",
)
class AflowPrototypeMatcher:
    """
    This class will match structures to their crystal prototypes, and will
    attempt to group species together to match structures derived from
    prototypes (e.g. an A_xB_1-x_C from a binary prototype), and will
    give these the names the "-like" suffix.

    This class uses data from the AFLOW LIBRARY OF CRYSTALLOGRAPHIC PROTOTYPES.
    If using this class, please cite their publication appropriately:

    Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
    The AFLOW library of crystallographic prototypes: part 1.
    Computational Materials Science, 136, S1-S828.
    https://doi.org/10.1016/j.commatsci.2017.01.017
    """

    def __init__(
        self,
        initial_ltol: float = 0.2,
        initial_stol: float = 0.3,
        initial_angle_tol: float = 5,
    ) -> None:
        """
        Tolerances as defined in StructureMatcher. Tolerances will be
        gradually decreased until only a single match is found (if possible).

        Args:
            initial_ltol (float): fractional length tolerance.
            initial_stol (float): site tolerance.
            initial_angle_tol (float): angle tolerance.
        """
        self.initial_ltol = initial_ltol
        self.initial_stol = initial_stol
        self.initial_angle_tol = initial_angle_tol

        # Preprocess AFLOW prototypes
        self._aflow_prototype_library: list[tuple[Structure, dict]] = []
        for dct in AFLOW_PROTOTYPE_LIBRARY:
            structure: Structure = dct["snl"].structure
            reduced_structure = self._preprocess_structure(structure)
            self._aflow_prototype_library.append((reduced_structure, dct))

    @staticmethod
    def _preprocess_structure(structure: Structure) -> Structure:
        return structure.get_reduced_structure(reduction_algo="niggli").get_primitive_structure()

    def _match_prototype(
        self,
        structure_matcher: StructureMatcher,
        reduced_structure: Structure,
    ) -> list[dict]:
        tags = []
        for aflow_reduced_structure, dct in self._aflow_prototype_library:
            # Since both structures are already reduced, we can skip the structure reduction step
            match = structure_matcher.fit_anonymous(
                aflow_reduced_structure, reduced_structure, skip_structure_reduction=True
            )
            if match:
                tags.append(dct)
        return tags

    def _match_single_prototype(self, structure: Structure) -> list[dict]:
        sm = StructureMatcher(
            ltol=self.initial_ltol,
            stol=self.initial_stol,
            angle_tol=self.initial_angle_tol,
            primitive_cell=True,
        )
        reduced_structure = self._preprocess_structure(structure)
        tags = self._match_prototype(sm, reduced_structure)
        while len(tags) > 1:
            sm.ltol *= 0.8
            sm.stol *= 0.8
            sm.angle_tol *= 0.8
            tags = self._match_prototype(sm, reduced_structure)
            if sm.ltol < 0.01:
                break
        return tags

    def get_prototypes(self, structure: Structure) -> list[dict] | None:
        """Get prototype(s) structures for a given input structure. If you use this method in
        your work, please cite the appropriate AFLOW publication:

            Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo,
            S. (2017). The AFLOW library of crystallographic prototypes: part 1. Computational
            Materials Science, 136, S1-S828. https://doi.org/10.1016/j.commatsci.2017.01.017

        Args:
            structure (Structure): structure to match

        Returns:
            list[dict] | None: A list of dicts with keys "snl" for the matched prototype and
                "tags", a dict of tags ("mineral", "strukturbericht" and "aflow") of that
                prototype. This should be a list containing just a single entry, but it is
                possible a material can match multiple prototypes.
        """
        tags: list[dict] = self._match_single_prototype(structure)

        return tags or None


def split_alpha_numeric(s: str) -> dict[str, list[str]]:
    """Split a string into separate lists of alpha and numeric groups.

    Args:
        s (str): The input string to split.

    Returns:
        dict[str, list[str]]: A dictionary with keys 'alpha' and 'numeric',
                              each containing a list of the respective groups.
    """
    groups = ["".join(g) for _, g in groupby(s, str.isalpha)]
    return {
        "alpha": [g for g in groups if g.isalpha()],
        "numeric": [g for g in groups if g.isnumeric()],
    }


def count_values_for_wyckoff(
    element_wyckoffs: list[str],
    counts: list[str],
    spg_num: str,
    lookup_dict: dict[str, dict[str, int]],
):
    """Count values from a lookup table and scale by wyckoff multiplicities."""
    return sum(
        int(count) * lookup_dict[spg_num][wyckoff_letter]
        for count, wyckoff_letter in zip(counts, element_wyckoffs, strict=True)
    )


@due.dcite(
    Doi("10.1126/sciadv.abn4117"),
    description="Rapid discovery of stable materials by coordinate-free coarse graining.",
)
def get_protostructure_label(
    struct: Structure,
    method: Literal["aflow", "spglib", "moyopy"],
    raise_errors: bool = False,
    **kwargs,
) -> str | None:
    """Get protostructure label for a pymatgen Structure.

    Args:
        struct (Structure): pymatgen Structure
        method (Literal["aflow", "spglib", "moyopy"]): Method to use for symmetry
        detection
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.
        **kwargs: Additional arguments for the specific method

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    if method == "aflow":
        return get_protostructure_label_from_aflow(struct, raise_errors, **kwargs)
    if method == "spglib":
        return get_protostructure_label_from_spglib(struct, raise_errors, **kwargs)
    if method == "moyopy":
        return get_protostructure_label_from_moyopy(struct, raise_errors, **kwargs)
    raise ValueError(f"Invalid method: {method}")


def get_protostructure_label_from_aflow(
    struct: Structure,
    raise_errors: bool = False,
    aflow_executable: str | None = None,
) -> str:
    """Get protostructure label for a pymatgen Structure. Make sure you're running a
    recent version of the aflow CLI as there's been several breaking changes. This code
    was tested under v3.2.12. The protostructure label is constructed as
    `aflow_label:chemsys`.

    Install guide: https://aflow.org/install-aflow/#install_aflow
        http://aflow.org/install-aflow/install-aflow.sh -o install-aflow.sh
        chmod 555 install-aflow.sh
        ./install-aflow.sh --slim

    Args:
        struct (Structure): pymatgen Structure
        aflow_executable (str): path to aflow executable. Defaults to which("aflow").
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    if aflow_executable is None:
        aflow_executable = which("aflow")

    if which(aflow_executable or "") is None:
        raise FileNotFoundError(
            "AFLOW could not be found, please specify path to its binary with aflow_executable='...'"
        )

    cmd = f"{aflow_executable} --prototype --print=json cat".split()

    output = subprocess.run(
        cmd,
        input=struct.to(fmt="poscar"),
        text=True,
        capture_output=True,
        check=True,
    )

    aflow_proto = orjson.loads(output.stdout)

    aflow_label = aflow_proto["aflow_prototype_label"]
    chemsys = struct.chemical_system
    # check that multiplicities satisfy original composition
    prototype_form, pearson_symbol, spg_num, *element_wyckoffs = aflow_label.split("_")

    element_dict = {}
    for elem, wyk_letters_per_elem in zip(chemsys.split("-"), element_wyckoffs, strict=True):
        # normalize Wyckoff letters to start with 1 if missing digit
        wyk_letters_normalized = re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, wyk_letters_per_elem)
        sep_el_wyks = split_alpha_numeric(wyk_letters_normalized)
        element_dict[elem] = count_values_for_wyckoff(
            sep_el_wyks["alpha"],
            sep_el_wyks["numeric"],
            spg_num,
            WYCKOFF_MULTIPLICITY_DICT,
        )

    element_wyckoffs = "_".join(element_wyckoffs)
    element_wyckoffs = canonicalize_element_wyckoffs(element_wyckoffs, spg_num)

    protostructure_label = f"{prototype_form}_{pearson_symbol}_{spg_num}_{element_wyckoffs}:{chemsys}"

    observed_formula = Composition(element_dict).reduced_formula
    expected_formula = struct.composition.reduced_formula
    if observed_formula != expected_formula:
        err_msg = (
            f"Invalid WP multiplicities - {protostructure_label}, expected {observed_formula} to be {expected_formula}"
        )
        if raise_errors:
            raise ValueError(err_msg)

        return err_msg

    return protostructure_label


def _get_all_wyckoffs_substring_and_element_dict(
    equivalent_wyckoff_labels: list[tuple[int, str, str]],
    spg_num: int | str,
):
    """Get Wyckoff position substring and element dict from equivalent Wyckoff labels.

    Args:
        equivalent_wyckoff_labels (list[tuple[int, str, str]]): List of tuples containing
            (multiplicity, element symbol, Wyckoff letter).
        spg_num (int | str): Space group number.

    Returns:
        tuple[str, dict]: Tuple containing:
            - str: Wyckoff position substring
            - dict: Dictionary mapping element symbols to their multiplicities
    """
    # Pre-sort by element and wyckoff letter to ensure continuous groups in groupby
    equivalent_wyckoff_labels = sorted(equivalent_wyckoff_labels, key=lambda x: (x[1], x[2]))

    # check that multiplicities satisfy original composition
    element_dict = {}
    element_wyckoffs = []
    for el, group in groupby(equivalent_wyckoff_labels, key=lambda x: x[1]):
        # NOTE create a list from the iterator so that we can use it without exhausting
        list_group = list(group)
        element_dict[el] = sum(WYCKOFF_MULTIPLICITY_DICT[str(spg_num)][e[2]] for e in list_group)
        # group by Wyckoff letter to get Wyckoff site multiplicity from len
        element_wyckoffs.append(
            "".join(
                f"{len(list(occurrences))}{wyk_letter}"
                for wyk_letter, occurrences in groupby(list_group, key=lambda x: x[2])
            )
        )
    all_wyckoffs = "_".join(element_wyckoffs)
    all_wyckoffs = canonicalize_element_wyckoffs(all_wyckoffs, spg_num)

    return all_wyckoffs, element_dict


def get_protostructure_label_from_spg_analyzer(
    spg_analyzer: SpacegroupAnalyzer,
    raise_errors: bool = False,
) -> str:
    """Get protostructure label for pymatgen SpacegroupAnalyzer.

    Args:
        spg_analyzer (SpacegroupAnalyzer): pymatgen SpacegroupAnalyzer object.
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    sym_struct = spg_analyzer.get_symmetrized_structure()

    spg_num = spg_analyzer.get_space_group_number()
    pearson_symbol = spg_analyzer.get_pearson_symbol()
    prototype_form = get_prototype_formula_from_composition(sym_struct.composition)
    chemsys = sym_struct.chemical_system

    # get Wyckoff position substring
    equivalent_wyckoff_labels = [
        # tuple of (wp multiplicity, element, wyckoff letter)
        (len(s), s[0].species_string, wyk_letter.translate(remove_digits))
        for s, wyk_letter in zip(sym_struct.equivalent_sites, sym_struct.wyckoff_symbols, strict=True)
    ]

    all_wyckoffs, element_dict = _get_all_wyckoffs_substring_and_element_dict(equivalent_wyckoff_labels, spg_num)

    protostructure_label = f"{prototype_form}_{pearson_symbol}_{spg_num}_{all_wyckoffs}:{chemsys}"

    observed_formula = Composition(element_dict).reduced_formula
    expected_formula = sym_struct.composition.reduced_formula
    if observed_formula != expected_formula:
        err_msg = (
            f"Invalid WP multiplicities - {protostructure_label}, expected {observed_formula} to be {expected_formula}"
        )
        if raise_errors:
            raise ValueError(err_msg)

        return err_msg

    return protostructure_label


def get_protostructure_label_from_spglib(
    struct: Structure,
    raise_errors: bool = False,
    init_symprec: float = 0.1,
    fallback_symprec: float | None = 1e-5,
) -> str:
    """Get AFLOW prototype label for pymatgen Structure.

    Args:
        struct (Structure): pymatgen Structure object.
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.
        init_symprec (float): Initial symmetry precision for spglib. Defaults to 0.1.
        fallback_symprec (float): Fallback symmetry precision for spglib if first
            symmetry detection failed. Defaults to 1e-5.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    attempt_to_recover = False
    try:
        spg_analyzer = SpacegroupAnalyzer(struct, symprec=init_symprec, angle_tolerance=5)
        try:
            aflow_label_with_chemsys = get_protostructure_label_from_spg_analyzer(spg_analyzer, raise_errors)

            if ("Invalid" in aflow_label_with_chemsys) and fallback_symprec is not None:
                attempt_to_recover = True
        except ValueError:
            if fallback_symprec is None:
                raise
            attempt_to_recover = True

        # try again with refined structure if it initially fails
        # NOTE structures with magmoms fail unless all have same magnetic moment
        if attempt_to_recover:
            spg_analyzer = SpacegroupAnalyzer(
                spg_analyzer.get_refined_structure(),
                symprec=fallback_symprec,
                angle_tolerance=-1,
            )
            aflow_label_with_chemsys = get_protostructure_label_from_spg_analyzer(spg_analyzer, raise_errors)
        return aflow_label_with_chemsys  # type: ignore[possibly-undefined]

    except ValueError as exc:
        if not raise_errors:
            return str(exc)
        raise


def get_protostructure_label_from_moyopy(
    struct: Structure,
    raise_errors: bool = False,
    symprec: float = 0.1,
) -> str | None:
    """Get AFLOW prototype label using Moyopy for symmetry detection.

    Args:
        struct (Structure): pymatgen Structure object.
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.
        symprec (float): Initial symmetry precision for Moyopy. Defaults to 0.1.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    if not has_moyopy:
        raise ImportError("moyopy not found, run pip install moyopy")

    # Convert pymatgen Structure to Moyo Cell and get symmetry data
    moyo_cell = MoyoAdapter.from_structure(struct)
    moyo_data = moyopy.MoyoDataset(moyo_cell, symprec=symprec)

    # Get space group number and Pearson symbol
    spg_num = moyo_data.number
    pearson_symbol = moyo_data.pearson_symbol
    prototype_form = get_prototype_formula_from_composition(struct.composition)
    chemsys = struct.chemical_system

    # Group Wyckoff positions by orbit and element
    equivalent_wyckoff_labels = []
    orbit_groups: dict[int, list[int]] = {}

    # Group sites by orbit
    for idx, orbit_id in enumerate(moyo_data.orbits):
        if orbit_id not in orbit_groups:
            orbit_groups[orbit_id] = []
        orbit_groups[orbit_id].append(idx)

    # Create equivalent_wyckoff_labels from orbit groups
    for orbit in orbit_groups.values():
        # All sites in an orbit have the same Wyckoff letter and element
        wyckoff = moyo_data.wyckoffs[orbit[0]]
        element = struct.species[orbit[0]]
        equivalent_wyckoff_labels += [(len(orbit), element.symbol, wyckoff.translate(remove_digits))]

    all_wyckoffs, element_dict = _get_all_wyckoffs_substring_and_element_dict(equivalent_wyckoff_labels, spg_num)

    protostructure_label = f"{prototype_form}_{pearson_symbol}_{spg_num}_{all_wyckoffs}:{chemsys}"

    # Verify multiplicities match composition
    observed_formula = Composition(element_dict).reduced_formula
    expected_formula = struct.composition.reduced_formula
    if observed_formula != expected_formula:
        err_msg = (
            f"Invalid WP multiplicities - {protostructure_label}, expected {observed_formula} to be {expected_formula}"
        )
        if raise_errors:
            raise ValueError(err_msg)
        return err_msg

    return protostructure_label


def canonicalize_element_wyckoffs(element_wyckoffs: str, spg_num: int | str) -> str:
    """Given an element ordering, canonicalize the associated Wyckoff positions
    based on the alphabetical weight of equivalent choices of origin.

    Args:
        element_wyckoffs (str): wyckoff substring section from aflow_label with the
            wyckoff letters for different elements separated by underscores.
        spg_num (int | str): International space group number.

    Returns:
        str: element_wyckoff string with canonical ordering of the wyckoff letters.
    """
    isopointal_element_wyckoffs = list(
        {element_wyckoffs.translate(str.maketrans(trans)) for trans in WYCKOFF_POSITION_RELAB_DICT[str(spg_num)]}
    )

    scored_element_wyckoffs = [
        sort_and_score_element_wyckoffs(element_wyckoffs) for element_wyckoffs in isopointal_element_wyckoffs
    ]

    return min(scored_element_wyckoffs, key=lambda x: (x[1], x[0]))[0]


def sort_and_score_element_wyckoffs(element_wyckoffs: str) -> tuple[str, int]:
    """Determines the order or Wyckoff positions when canonicalizing AFLOW labels.

    Args:
        element_wyckoffs (str): wyckoff substring section from aflow_label with the
            wyckoff letters for different elements separated by underscores.

    Returns:
        tuple: containing
        - str: sorted Wyckoff position substring for AFLOW-style prototype label
        - int: integer score to rank order when canonicalizing
    """
    score = 0
    sorted_element_wyckoffs = []
    for el_wyks in element_wyckoffs.split("_"):
        wp_counts = split_alpha_numeric(el_wyks)
        # TODO: strict=False here because the regex is not perfect and sometimes
        # adds additional "1"s to the numeric list
        sorted_element_wyckoffs.append(
            "".join(
                f"{count}{wyckoff_letter}" if count != "1" else wyckoff_letter
                for count, wyckoff_letter in sorted(
                    zip(wp_counts["numeric"], wp_counts["alpha"], strict=False),
                    key=lambda x: x[1],
                )
            )
        )
        score += sum(0 if wyckoff_letter == "A" else ord(wyckoff_letter) - 96 for wyckoff_letter in wp_counts["alpha"])

    return "_".join(sorted_element_wyckoffs), score


def get_prototype_formula_from_composition(composition: Composition) -> str:
    """An anonymized formula. Unique species are arranged in alphabetical order
    and assigned ascending alphabets. This format is used in the aflow structure
    prototype labelling scheme.

    Args:
        composition (Composition): Pymatgen Composition to process

    Returns:
        str: anonymized formula where the species are in alphabetical order
    """
    reduced = composition.element_composition
    if all(x == int(x) for x in composition.values()):
        reduced /= gcd(*(int(amt) for amt in composition.values()))

    amounts = [reduced[key] for key in sorted(reduced, key=str)]

    anon = ""
    # NOTE:strict=False because ascii_uppercase is full list of uppercase letters
    for elem, amt in zip(ascii_uppercase, amounts, strict=False):
        if amt == 1:
            amt_str = ""
        elif abs(amt % 1) < 1e-8:
            amt_str = str(int(amt))
        else:
            amt_str = str(amt)
        anon += f"{elem}{amt_str}"
    return anon


def get_anonymous_formula_from_prototype_formula(prototype_formula: str) -> str:
    """Get an anonymous formula from a prototype formula."""
    prototype_formula = re.sub(RE_ELEMENT_NO_SUFFIX, RE_SUBST_ONE_SUFFIX, prototype_formula)
    anom_list = split_alpha_numeric(prototype_formula)

    return "".join(
        f"{el}{num}" if num != 1 else el
        for el, num in zip(anom_list["alpha"], sorted(map(int, anom_list["numeric"])), strict=True)
    )


def get_formula_from_protostructure_label(protostructure_label: str) -> str:
    """Get a formula from a protostructure label."""
    aflow_label, chemsys = protostructure_label.split(":")
    prototype_formula = aflow_label.split("_")[0]
    prototype_formula = re.sub(RE_ELEMENT_NO_SUFFIX, RE_SUBST_ONE_SUFFIX, prototype_formula)
    anom_list = split_alpha_numeric(prototype_formula)

    return "".join(
        f"{el}{num}" if num != 1 else el
        for el, num in zip(chemsys.split("-"), map(int, anom_list["numeric"]), strict=True)
    )


def count_distinct_wyckoff_letters(protostructure_label: str) -> int:
    """Count number of distinct Wyckoff letters in protostructure_label.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        int: number of distinct Wyckoff letters in protostructure_label
    """
    aflow_label, _ = protostructure_label.split(":")
    _, _, _, element_wyckoffs = aflow_label.split("_", 3)
    element_wyckoffs = element_wyckoffs.translate(remove_digits).replace("_", "")
    return len(set(element_wyckoffs))  # number of distinct Wyckoff letters


def count_wyckoff_positions(protostructure_label: str) -> int:
    """Count number of Wyckoff positions in protostructure_label.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        int: number of distinct Wyckoff positions in protostructure_label
    """
    aflow_label, _ = protostructure_label.split(":")  # remove chemical system
    # discard prototype formula and spg symbol and spg number
    wyk_letters = aflow_label.split("_", maxsplit=3)[-1]
    # throw Wyckoff positions for all elements together
    wyk_letters = wyk_letters.replace("_", "")
    wyk_list = re.split("[A-z]", wyk_letters)[:-1]  # split on every letter

    # count 1 for letters without prefix
    return sum(1 if len(x) == 0 else int(x) for x in wyk_list)


def count_crystal_dof(protostructure_label: str) -> int:
    """Count number of free parameters in coarse-grained protostructure_label
    representation: how many degrees of freedom would remain to optimize during
    a crystal structure relaxation.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        int: Number of free-parameters in given prototype
    """
    aflow_label, _ = protostructure_label.split(":")  # chop off chemical system
    _, pearson_symbol, spg_num, *element_wyckoffs = aflow_label.split("_")

    return (
        _count_from_dict(element_wyckoffs, WYCKOFF_POSITION_PARAM_DICT, spg_num)
        + CRYSTAL_LATTICE_PARAMETERS_COUNTS[pearson_symbol[0]]
    )


def count_crystal_sites(protostructure_label: str) -> int:
    """Count number of sites from protostructure_label.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        int: Number of free-parameters in given prototype
    """
    aflow_label, _ = protostructure_label.split(":")  # chop off chemical system
    _, _, spg_num, *element_wyckoffs = aflow_label.split("_")

    return _count_from_dict(element_wyckoffs, WYCKOFF_MULTIPLICITY_DICT, spg_num)


def _count_from_dict(element_wyckoffs: list[str], lookup_dict: dict, spg_num: str) -> int:
    """Count number of sites from protostructure_label."""
    n_params = 0

    for wyckoffs in element_wyckoffs:
        # normalize Wyckoff letters to start with 1 if missing digit
        sep_el_wyks = split_alpha_numeric(re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, wyckoffs))
        n_params += count_values_for_wyckoff(
            sep_el_wyks["alpha"],
            sep_el_wyks["numeric"],
            spg_num,
            lookup_dict,
        )

    return int(n_params)


def get_prototype_from_protostructure(protostructure_label: str) -> str:
    """Get a canonicalized string for the prototype. This prototype should be
    the same for all isopointal protostructures.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        str: Canonicalized AFLOW-style prototype label
    """
    aflow_label, _ = protostructure_label.split(":")
    prototype_formula, pearson_symbol, spg_num, *element_wyckoffs = aflow_label.split("_")

    anonymous_formula = get_anonymous_formula_from_prototype_formula(prototype_formula)
    counts = [
        int(x)
        for x in split_alpha_numeric(re.sub(RE_ELEMENT_NO_SUFFIX, RE_SUBST_ONE_SUFFIX, prototype_formula))["numeric"]
    ]

    # map to list to avoid mypy error, zip returns tuples.
    counts, element_wyckoffs = map(list, zip(*sorted(zip(counts, element_wyckoffs, strict=True)), strict=True))
    all_wyckoffs = "_".join(element_wyckoffs)
    all_wyckoffs = re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, all_wyckoffs)
    if len(counts) == len(set(counts)):
        all_wyckoffs = canonicalize_element_wyckoffs(all_wyckoffs, int(spg_num))
        return f"{anonymous_formula}_{pearson_symbol}_{spg_num}_{all_wyckoffs}"

    # credit Stef: https://stackoverflow.com/a/70126643/5517459
    all_wyckoffs_permutations = [
        "_".join(list(map(itemgetter(1), chain.from_iterable(p))))
        for p in product(
            *[
                permutations(g)
                for _, g in groupby(sorted(zip(counts, all_wyckoffs.split("_"), strict=True)), key=lambda x: x[0])
            ]
        )
    ]

    isopointal_all_wyckoffs = list(
        {
            all_wyckoffs.translate(str.maketrans(trans))
            for all_wyckoffs in all_wyckoffs_permutations
            for trans in WYCKOFF_POSITION_RELAB_DICT[spg_num]
        }
    )

    scored_all_wyckoffs = [
        sort_and_score_element_wyckoffs(element_wyckoffs) for element_wyckoffs in isopointal_all_wyckoffs
    ]

    all_wyckoffs = min(scored_all_wyckoffs, key=lambda x: (x[1], x[0]))[0]

    return f"{anonymous_formula}_{pearson_symbol}_{spg_num}_{all_wyckoffs}"


def _get_anonymous_formula_dict(anonymous_formula: str) -> dict:
    """Get a dictionary of element to count from an anonymous formula."""
    result: defaultdict = defaultdict(int)
    element = ""
    count = ""

    for char in anonymous_formula:
        if char.isalpha():
            if element:
                result[element] += int(count) if count else 1
                count = ""
            element = char
        else:
            count += char

    if element:
        result[element] += int(count) if count else 1

    return dict(result)


def _find_translations(dict1: dict[str, int], dict2: dict[str, int]) -> list[dict[str, str]]:
    """Find all possible translations between two dictionaries."""
    if Counter(dict1.values()) != Counter(dict2.values()):
        return []

    keys2 = list(dict2.keys())
    used = set()

    def backtrack(translation, index):
        if index == len(dict1):
            return [translation.copy()]

        key1 = list(dict1.keys())[index]
        value1 = dict1[key1]
        valid_translations = []

        for key2 in keys2:
            if key2 not in used and dict2[key2] == value1:
                used.add(key2)
                translation[key1] = key2
                valid_translations.extend(backtrack(translation, index + 1))
                used.remove(key2)
                del translation[key1]

        return valid_translations

    return backtrack({}, 0)


def get_protostructures_from_aflow_label_and_composition(aflow_label: str, composition: Composition) -> list[str]:
    """Get a canonicalized string for the prototype.

    Args:
        aflow_label (str): AFLOW-style prototype label
        composition (Composition): pymatgen Composition object

    Returns:
        list[str]: List of possible protostructure labels that can be generated
            from combinations of the input aflow_label and composition.
    """
    anonymous_formula, pearson_symbol, spg_num, *element_wyckoffs = aflow_label.split("_")

    ele_amt_dict = composition.get_el_amt_dict()
    proto_formula = get_prototype_formula_from_composition(composition)
    anom_amt_dict = _get_anonymous_formula_dict(anonymous_formula)

    translations = _find_translations(ele_amt_dict, anom_amt_dict)
    anom_ele_to_wyk = dict(zip(anom_amt_dict.keys(), element_wyckoffs, strict=True))
    anonymous_formula = RE_ANONYMOUS.sub(RE_SUBST_ONE_PREFIX, anonymous_formula)

    protostructures = set()
    for t in translations:
        wyckoff_part = "_".join(
            RE_WYCKOFF.sub(RE_SUBST_ONE_PREFIX, anom_ele_to_wyk[t[elem]]) for elem in sorted(t.keys())
        )
        canonicalized_wyckoff = canonicalize_element_wyckoffs(wyckoff_part, spg_num)
        chemical_system = "-".join(sorted(t.keys()))

        protostructures.add(f"{proto_formula}_{pearson_symbol}_{spg_num}_{canonicalized_wyckoff}:{chemical_system}")

    return list(protostructures)


def get_random_structure_for_protostructure(protostructure_label: str, **kwargs) -> Structure:
    """Generate a random structure for a given prototype structure.

    NOTE that due to the random nature of the generation, the output structure
    may be higher symmetry than the requested prototype structure.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.
        **kwargs: Keyword arguments to pass to pyxtal().from_random()
    """
    if not has_pyxtal:
        raise ImportError("pyxtal is required for this function")

    aflow_label, chemsys = protostructure_label.split(":")
    _, _, spg_num, *element_wyckoffs = aflow_label.split("_")

    sep_el_wyks = [split_alpha_numeric(re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, w)) for w in element_wyckoffs]

    species_sites = [
        [
            site
            for count, wyckoff_letter in zip(d["numeric"], d["alpha"], strict=True)
            for site in [f"{WYCKOFF_MULTIPLICITY_DICT[spg_num][wyckoff_letter]}{wyckoff_letter}"] * int(count)
        ]
        for d in sep_el_wyks
    ]

    species_counts = [
        sum(
            WYCKOFF_MULTIPLICITY_DICT[spg_num][wyckoff_letter] * int(count)
            for count, wyckoff_letter in zip(d["numeric"], d["alpha"], strict=True)
        )
        for d in sep_el_wyks
    ]

    p = pyxtal()
    p.from_random(
        dim=3,
        group=int(spg_num),
        species=chemsys.split("-"),
        numIons=species_counts,
        sites=species_sites,
        **kwargs,
    )
    return p.to_pymatgen()
