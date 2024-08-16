"""
This module implements input/output processing from Multiwfn.

Currently, the focus of this module is on processing Quantum Theory of Atoms in Molecules (QTAIM) outputs from
Multiwfn. Additional features may be added over time as needed/desired.
"""

from __future__ import annotations

import copy
import warnings
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

if TYPE_CHECKING:
    from pymatgen.core.structure import Molecule
    from pymatgen.util.typing import PathLike

__author__ = "Santiago Vargas, Evan Spotte-Smith"
__copyright__ = "Copyright 2024, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "santiagovargas921@gmail.com, espottesmith@gmail.com"
__date__ = "July 10, 2024"


QTAIM_CONDITIONALS = {
    "cp_num": ["----------------"],
    "ele_info": ["Corresponding", "nucleus:"],
    "connected_bond_paths": ["Connected", "atoms:"],
    "pos_ang": ["Position", "(Angstrom):"],
    "density_total": ["Density", "of", "all", "electrons:"],
    "density_alpha": ["Density", "of", "Alpha", "electrons:"],
    "density_beta": ["Density", "of", "Beta", "electrons:"],
    "spin_density": ["Spin", "density", "of", "electrons:"],
    "lol": ["Localized", "orbital", "locator"],
    "energy_density": ["Energy", "density", "E(r)"],
    "Lagrangian_K": ["Lagrangian", "kinetic", "energy"],
    "Hamiltonian_K": ["Hamiltonian", "kinetic", "energy"],
    "lap_e_density": ["Laplacian", "electron", "density:"],
    "e_loc_func": ["Electron", "localization", "function"],
    "ave_loc_ion_E": ["Average", "local", "ionization", "energy"],
    "delta_g_promolecular": ["Delta-g", "promolecular"],
    "delta_g_hirsh": ["Delta-g", "Hirshfeld"],
    "esp_nuc": ["ESP", "nuclear", "charges:"],
    "esp_e": ["ESP", "electrons:"],
    "esp_total": ["Total", "ESP:"],
    "grad_norm": ["Components", "gradient", "x/y/z"],
    "lap_norm": ["Components", "Laplacian", "x/y/z"],
    "eig_hess": ["Eigenvalues", "Hessian:"],
    "det_hessian": ["Determinant", "Hessian:"],
    "ellip_e_dens": ["Ellipticity", "electron", "density:"],
    "eta": ["eta", "index:"],
}


def extract_info_from_cp_text(
    lines_split: list[list[str]], cp_type: Literal["atom", "bond", "ring", "cage"], conditionals: dict[str, list[str]]
) -> tuple[str, dict[str, Any]]:
    """
    Extract specific information from a Multiwfn QTAIM output.

    Args:
        lines_split (List[List[str]]): List of lines from a (pre-processed) CP file, containing only information
            regarding one critical point, split by whitespace
        cp_type: Type of critical point. Currently, can be "atom", "bond", "ring", or "cage"
        conditionals (Dict[str, List[str]]): Parameters to extract with strings to search for to see if the
            data is present

    Returns:
        cp_dict: dictionary of critical point information
    """

    cp_dict: dict[str, Any] = dict()

    cp_name = "null"

    this_conditionals = copy.deepcopy(conditionals)

    for ind, i in enumerate(lines_split):
        for k, v in this_conditionals.items():
            if all(x in i for x in v):
                if k == "cp_num":
                    cp_dict[k] = int(i[2][:-1])

                    # Placeholder name
                    # For nuclear critical points, this can be overwritten later (see below)
                    cp_name = f"{cp_dict[k]}_{cp_type}"

                elif k == "ele_info":
                    if i[2] == "Unknown":
                        cp_name = str(cp_dict["cp_num"]) + "_Unknown"
                        cp_dict["number"] = "Unknown"
                        cp_dict["ele"] = "Unknown"
                    else:
                        if len(i) == 3:
                            cp_dict["element"] = i[2].split("(")[1][:-1]
                            cp_dict["number"] = i[2].split("(")[0]
                        else:
                            cp_dict["element"] = i[2].split("(")[1]
                            cp_dict["number"] = i[2].split("(")[0]
                        cp_name = cp_dict["number"] + "_" + cp_dict["element"]

                elif k == "connected_bond_paths":
                    list_raw = list(i[2:])
                    # save only items that contain a number
                    list_raw = [x for x in list_raw if any(char.isdigit() for char in x)]  # type: ignore[misc]
                    list_raw = [int(x.split("(")[0]) for x in list_raw]  # type: ignore[misc]
                    cp_dict[k] = list_raw

                elif k == "pos_ang":
                    cp_dict[k] = [float(x) for x in i[2:]]

                elif k == "esp_total":
                    cp_dict[k] = float(i[2])

                elif k == "eig_hess":
                    cp_dict[k] = np.sum(np.array([float(x) for x in i[-3:]]))

                elif k in ("grad_norm", "lap_norm"):
                    cp_dict[k] = float(lines_split[ind + 2][-1])

                else:
                    cp_dict[k] = float(i[-1])

                this_conditionals.pop(k)
                break

    return cp_name, cp_dict


def parse_cp(lines: list[str]) -> tuple[str | None, dict[str, Any]]:
    """
    Parse information from a single QTAIM critical point.

    Args:
        lines (List[str]): list of lines from a (preparsed) CP file, containing only information regarding one
            critical point
    Returns:
        cp_dict: dictionary of critical point information
    """
    lines_split = [line.split() for line in lines]

    cp_type = None

    # Figure out what kind of critical-point we're dealing with
    if "(3,-3)" in lines_split[0]:
        cp_type = "atom"
        conditionals = {k: v for k, v in QTAIM_CONDITIONALS.items() if k != "connected_bond_paths"}
    elif "(3,-1)" in lines_split[0]:
        cp_type = "bond"
        conditionals = {k: v for k, v in QTAIM_CONDITIONALS.items() if k != "ele_info"}
    elif "(3,+1)" in lines_split[0]:
        cp_type = "ring"
        conditionals = {k: v for k, v in QTAIM_CONDITIONALS.items() if k not in ["connected_bond_paths", "ele_info"]}
    elif "(3,+3)" in lines_split[0]:
        cp_type = "cage"
        conditionals = {k: v for k, v in QTAIM_CONDITIONALS.items() if k not in ["connected_bond_paths", "ele_info"]}
    else:
        return None, dict()

    return extract_info_from_cp_text(lines_split, cp_type, conditionals)  # type: ignore[arg-type]


def get_qtaim_descs(file: PathLike) -> dict[str, dict[str, Any]]:
    """
    Parse CPprop file from multiwfn by parsing each individual critical-point section.

    Args:
        file (PathLike): path to CPprop file

    Returns:
        descriptors (Dict[str, Dict[str, Any]]): Output dictionary of QTAIM descriptors

    """

    cp_sections = list()
    descriptors = dict()

    with open(file) as f:
        lines = f.readlines()
        lines = [line[:-1] for line in lines]

    # section lines into segments on ----------------
    lines_segment: list[str] = list()
    for ind, line in enumerate(lines):
        if "----------------" in line:
            lines_segment = list()
        lines_segment.append(line)
        if ind < len(lines) - 1:
            if "----------------" in lines[ind + 1]:
                cp_sections.append(lines_segment)
        else:
            cp_sections.append(lines_segment)

    for section in cp_sections:
        cp_name, cp_dict = parse_cp(section)

        # Possibility that parsing fails
        if cp_name:
            descriptors[cp_name] = cp_dict

    return descriptors


def separate_cps_by_type(qtaim_descs: dict[Any, dict[str, Any]]) -> dict[str, dict[Any, dict[str, Any]]]:
    """
    Separates QTAIM descriptors by type (atom, bond, ring, or cage)

    Args:
        qtaim_descs (Dict[str, Dict[str, Any]]): Dictionary where keys are CP names and values are dictionaries of
            descriptors obtained from `get_qtaim_descs` and `parse_cp`

    Returns:
        organized_descs (Dict[str, Dict[str, Dict[str, Any]]]): Dictionary of organized QTAIM critical points and their
            descriptors. Keys are "atom", "bond", "ring", and "cage", and values are dicts
            {<CP name>: <QTAIM descriptors>}
    """

    organized_descs: dict[str, dict[str, dict[str, Any]]] = {
        "atom": dict(),
        "bond": dict(),
        "ring": dict(),
        "cage": dict(),
    }

    for k, v in qtaim_descs.items():
        if "bond" in k:
            organized_descs["bond"][k] = v
        elif "ring" in k:
            organized_descs["ring"][k] = v
        elif "cage" in k:
            organized_descs["cage"][k] = v
        elif "Unknown" not in k:
            organized_descs["atom"][k] = v

    return organized_descs


def match_atom_cp(
    molecule: Molecule, index: int, atom_cp_dict: dict[str, dict[str, Any]], max_distance: float = 0.5
) -> tuple[str | None, dict]:
    """
    From a dictionary with an atom's position and element symbol, find the corresponding cp in the atom CP dictionary

    Args:
        molecule (Molecule): structure corresponding to this Multiwfn calculation
        index (int): index of the atom to be matched
        atom_cp_dict (Dict[str, Dict[str, Any]]): Dictionary where keys are critical point names and values are
            descriptors, including element symbol and position
        max_distance (float): Maximum distance (in Angstrom) that a critical point can be away from an atom center
            and be associated with that atom. Default is 0.5 Angstrom

    Returns:
        cp_name (str | None): Key of atom_cp_dict; the name of the atom critical point corresponding to this atom. If
            no match is found, this will be None
        cp_dict (Dict): Dictionary of CP descriptors matching this atom
    """

    atom = molecule.sites[index]
    atom_symbol = atom.species_string

    for cp_name, cp_dict in atom_cp_dict.items():
        if int(cp_name.split("_")[0]) == index + 1 and cp_dict["element"] == atom_symbol:
            return cp_name, cp_dict

        if cp_dict["element"] == atom_symbol:
            distance = np.linalg.norm(np.array(cp_dict["pos_ang"]) - atom.coords)

            if distance < max_distance:
                return cp_name, cp_dict

    # No match
    return None, dict()


def map_atoms_cps(
    molecule: Molecule, atom_cp_dict: dict[str, dict[str, Any]], max_distance: float = 0.5
) -> tuple[
    dict[int, dict[str, Any]],
    list[int],
]:
    """
    Connect atom CPs to atoms by their positions.

    Args:
        molecule (Molecule): structure corresponding to this Multiwfn calculation
        atom_cp_dict (Dict[str, Dict[str, Any]]): Dictionary where keys are critical point names and values are
            descriptors, including element symbol and position
        max_distance (float): Maximum distance (in Angstrom) that a critical point can be away from an atom center
            and be associated with that atom. Default is 0.5 Angstrom
    Returns:
        index_to_cp_desc (Dict[int, Dict[str, Any]]): Dictionary mapping atomic indices to atom critical point
            descriptors
        missing_atoms (List[int]): list of dft atoms that do not have a cp found in qtaim
    """

    index_to_cp_desc = dict()
    missing_atoms = list()

    for index in range(len(molecule)):
        # finds cp by distance and naming scheme from CPprop.txt
        cp_name, this_atom_cp = match_atom_cp(molecule, index, atom_cp_dict, max_distance=max_distance)

        # If this is False, that means no match was found
        if cp_name:
            index_to_cp_desc[index] = this_atom_cp
            index_to_cp_desc[index]["name"] = cp_name
        else:
            index_to_cp_desc[index] = dict()
            missing_atoms.append(index)

    return index_to_cp_desc, missing_atoms


def add_atoms(
    molecule: Molecule,
    organized_cps: dict[str, dict[Any, dict[str, Any]]],
    bond_atom_criterion: Literal["qtaim", "distance", "combined"] = "combined",
    dist_threshold_bond: float = 1.0,
    dist_threshold_ring_cage: float = 3.0,
    distance_margin: float = 0.5,
) -> dict[str, dict[str, dict[str, Any]]]:
    """
    Modify bond, ring, and cage CPs to include information about surrounding critical points. Bonds will include
    information about the atoms that make up the bond. Rings will include the names of neighboring bonds and the
    indices of the atoms involved, and cages will include information on neighboring rings, bonds, and the atoms
    involved.

    Args:
        molecule (Molecule): structure corresponding to this Multiwfn calculation
        organized_cps (Dict[str, Dict[Any, Dict[str, Any]]]): Keys are CP categories ("atom", "bond", "ring", and
            "cage"). Values are themselves dictionaries, where the keys are CP names (or atom indices) and the values
            are CP descriptors
        bond_atom_criterion (Literal["qtaim", "distance", "combined"]): If this is "qtaim", the inherent bonding
            definition obtained from QTAIM/Multiwfn will be used to link bond CPs to atoms involved in those bonds. If
            it is "distance", a distance-based metric will be used instead, where the atoms closest to the bond CP will
            be assumed to be bonded. If this is "combined", then the QTAIM/Multiwfn definition will be used where
            available, and a distance metric will be used for all bond CPs lacking a definition from QTAIM/Multiwfn.
            NOTE: to use "qtaim" or "combined" as `bond_atom_criterion`, you must have used `map_atoms_cps` to link atom
            numbers from Multiwfn to atom indices in `molecule`.
        dist_threshold_bond (float): If the nearest atoms to a bond CP are further from the bond CP than this threshold
            (default 1.0 Angstrom), then a warning will be raised.
        dist_threshold_ring_cage (float): If the nearest bond CPs to a ring CP or the nearest ring CPs to a cage CP are
            further than this threshold (default 3.0 Angstrom), then a warning will be raised.
        distance_margin (float): Rings must be made up of at least three bonds, and cages must be made up of at least
            three rings. We therefore define rings by the three closest bond CPs and cages by the three closest ring
            CPs. We associate additional bond/ring CPs with ring/cage CPs if they are no further than the third-closest
            bond/ring CP, plus this margin. Default is 0.5 Angstrom

    Returns:
        modified_organized_cps (Dict[str, Dict[str, Dict[str, Any]]]): CP dict with additional information added.
    """

    def sort_cps_by_distance(
        position: np.ndarray,
        options: dict[Any, dict[str, Any]],
    ):
        dists = list()
        for oname, oval in options.items():
            opos = np.array(oval["pos_ang"])
            dists.append((np.linalg.norm(position - opos), oname))

        return sorted(dists)

    modified_organized_cps = copy.deepcopy(organized_cps)

    atom_info = {i: {"pos_ang": s.coords} for i, s in enumerate(molecule.sites)}

    # Can only include bonds where the connected atoms are provided
    if bond_atom_criterion == "qtaim":
        modified_organized_cps["bond"] = {
            k: v for k, v in modified_organized_cps["bond"].items() if "connected_bond_paths" in v
        }

    atom_cps = modified_organized_cps["atom"]
    bond_cps = modified_organized_cps["bond"]
    ring_cps = modified_organized_cps["ring"]
    cage_cps = modified_organized_cps["cage"]

    if len(bond_cps) > 0 and len(atom_info) < 2:
        raise ValueError("Cannot have a bond CP with less than two atom CPs!")

    if len(ring_cps) > 0 and len(bond_cps) < 3:
        raise ValueError("Cannot have a ring CP with less than three bond CPs!")

    if len(cage_cps) > 0 and len(ring_cps) < 3:
        raise ValueError("Cannot have a cage CP with less than three ring CPs!")

    for cp_name, cp_desc in bond_cps.items():
        if bond_atom_criterion == "distance":
            # NOTE: for bonds, we associate based on atoms, NOT based on atom CPs
            sorted_atoms = sort_cps_by_distance(np.array(cp_desc["pos_ang"]), atom_info)

            if sorted_atoms[1][0] > dist_threshold_bond:
                warnings.warn("Warning: bond CP is far from bonding atoms")

            # Assume only two atoms involved in bond
            modified_organized_cps["bond"][cp_name]["atom_inds"] = sorted([ca[1] for ca in sorted_atoms[:2]])
        elif bond_atom_criterion == "qtaim":
            bond_atoms_list = list()
            for index in cp_desc["connected_bond_paths"]:
                for true_index, descriptors in atom_cps.items():
                    if int(descriptors["name"].split("_")[0]) == index:
                        bond_atoms_list.append(true_index)
                        break

            modified_organized_cps["bond"][cp_name]["atom_inds"] = sorted(bond_atoms_list)
        else:
            if "connected_bond_paths" in cp_desc:
                bond_atoms_list = list()
                for index in cp_desc["connected_bond_paths"]:
                    for true_index, descriptors in atom_cps.items():
                        if int(descriptors["name"].split("_")[0]) == index:
                            bond_atoms_list.append(true_index)
                            break
                if len(bond_atoms_list) != 2:
                    raise ValueError(f"Could not match all atoms in connected_bond_paths for bond CP {cp_name}")
            else:
                sorted_atoms = sort_cps_by_distance(np.array(cp_desc["pos_ang"]), atom_info)

                if sorted_atoms[1][0] > dist_threshold_bond:
                    warnings.warn("Warning: bond CP is far from bonding atoms")

                bond_atoms_list = sorted([ca[1] for ca in sorted_atoms[:2]])

            modified_organized_cps["bond"][cp_name]["atom_inds"] = sorted(bond_atoms_list)

    for cp_name, cp_desc in ring_cps.items():
        sorted_bonds = sort_cps_by_distance(np.array(cp_desc["pos_ang"]), bond_cps)
        max_close_dist = sorted_bonds[2][0]

        if max_close_dist > dist_threshold_ring_cage:
            warnings.warn("Warning: ring CP is far from closest bond CPs.")

        # Assume that the three closest bonds are all part of the ring
        bond_names = [bcp[1] for bcp in sorted_bonds[:3]]

        # Add additional bonds that are relatively close to this ring
        if len(sorted_bonds) > 3:
            for bond_dist, bond_name in sorted_bonds[3:]:
                if bond_dist < max_close_dist + distance_margin:
                    bond_names.append(bond_name)
                else:
                    break

        # Add all unique atoms involved in this ring
        atom_inds = set()
        for bond_name in bond_names:
            atom_inds.update(bond_cps[bond_name]["atom_inds"])

        modified_organized_cps["ring"][cp_name]["bond_names"] = bond_names
        modified_organized_cps["ring"][cp_name]["atom_inds"] = list(atom_inds)

    for cp_name, cp_desc in cage_cps.items():
        sorted_rings = sort_cps_by_distance(np.array(cp_desc["pos_ang"]), ring_cps)
        max_close_dist = sorted_rings[2][0]

        # Warn if the three closest bonds are further than the max distance
        if max_close_dist > dist_threshold_ring_cage:
            warnings.warn("Warning: cage CP is far from closest ring CPs.")

        # Assume that the three closest rings are all part of the cage
        ring_names = [rcp[1] for rcp in sorted_rings[:3]]

        # Add additional rings that are relatively close to this cage
        if len(sorted_rings) > 3:
            for ring_dist, ring_name in sorted_rings[3:]:
                if ring_dist < max_close_dist + distance_margin:
                    ring_names.append(ring_name)
                else:
                    break

        # Add all unique bonds and atoms involved in this cage
        bond_names_cage = set()
        atom_inds = set()
        for ring_name in ring_names:
            bond_names_cage.update(ring_cps[ring_name]["bond_names"])  # type: ignore[attr-defined]
            atom_inds.update(ring_cps[ring_name]["atom_inds"])  # type: ignore[attr-defined]

        modified_organized_cps["cage"][cp_name]["ring_names"] = ring_names
        modified_organized_cps["cage"][cp_name]["bond_names"] = list(bond_names_cage)
        modified_organized_cps["cage"][cp_name]["atom_inds"] = list(atom_inds)

    return modified_organized_cps


def process_multiwfn_qtaim(
    molecule: Molecule,
    file: PathLike,
    bond_atom_criterion: Literal["qtaim", "distance"] = "distance",
    max_distance_atom: float = 0.5,
    dist_threshold_bond: float = 1.0,
    dist_threshold_ring_cage: float = 3.0,
    distance_margin: float = 0.5,
) -> dict[str, dict[Any, dict[str, Any]]]:
    """
    Process quantum theory of atoms in molecules (QTAIM) outputs from Multiwfn.

    Args:
        molecule (Molecule): structure corresponding to this Multiwfn calculation
        file (PathLike): path to CPprop file containing QTAIM information
        bond_atom_criterion (Literal["qtaim", "distance"]): If "qtaim", the inherent bonding definition obtained from
            QTAIM/Multiwfn will be used to link bond CPs to atoms involved in those bonds; if "distance", a
            distance-based metric will be used instead, where the atoms closest to the bond CP will be assumed to be
            bonded.
        max_distance (float): Maximum distance (in Angstrom) that an atom critical point can be away from an atom
            center and be associated with that atom. Default is 0.5 Angstrom
        dist_threshold_bond (float): If the nearest atoms to a bond CP are further from the bond CP than this threshold
            (default 1.0 Angstrom), then a warning will be raised.
        dist_threshold_ring_cage (float): If the nearest bond CPs to a ring CP or the nearest ring CPs to a cage CP are
            further than this threshold (default 3.0 Angstrom), then a warning will be raised.
        distance_margin (float): Rings must be made up of at least three bonds, and cages must be made up of at least
            three rings. We therefore define rings by the three closest bond CPs and cages by the three closest ring
            CPs. We associate additional bond/ring CPs with ring/cage CPs if they are no further than the third-closest
            bond/ring CP, plus this margin. Default is 0.5 Angstrom

    Returns:
        with_atoms (Dict[str, Dict[str, Dict[str, Any]]]): QTAIM descriptors, organized by type ("atom", "bond",
            "ring", "cage"), with additional metadata added to bond, ring, and cage CPs
    """

    # Initial parsing and organizing
    descriptors_unorganized = get_qtaim_descs(file)
    qtaim_descriptors: dict[str, dict[Any, dict[str, Any]]] = separate_cps_by_type(descriptors_unorganized)

    # Remap atom CPs to atom indices
    remapped_atoms, missing_atoms = map_atoms_cps(molecule, qtaim_descriptors["atom"], max_distance=max_distance_atom)

    if len(missing_atoms) > 0:
        warnings.warn(f"Some atoms not mapped to atom CPs! Indices: {missing_atoms}")

    qtaim_descriptors["atom"] = remapped_atoms

    return add_atoms(
        molecule,
        qtaim_descriptors,
        bond_atom_criterion=bond_atom_criterion,
        dist_threshold_bond=dist_threshold_bond,
        dist_threshold_ring_cage=dist_threshold_ring_cage,
        distance_margin=distance_margin,
    )
