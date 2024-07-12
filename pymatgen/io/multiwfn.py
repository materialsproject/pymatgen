"""
This module implements input/output processing from Multiwfn.

Currently, the focus of this module is on processing Quantum Theory of Atoms in Molecules (QTAIM)
outputs from Multiwfn. Additional features may be added over time as needed/desired.
"""

import copy
import json
import os
import warnings
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import numpy as np


__author__ = "Santiago Vargas, Evan Spotte-Smith"
__copyright__ = "Copyright 2024, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "santiagovargas921@gmail.com, espottesmith@gmail.com"
__date__ = "July 10, 2024"


def extract_info_from_cp_text(
    lines_split: List[str],
    cp_type: str,
    conditionals: Dict[str, List[str]]
) -> Tuple[str, Dict[str, Any]]:
    """
    Extract specific information from a Multiwfn QTAIM output.

    Args:
        lines_split (List[str]): List of lines from a (preparsed) CP file, containing only information regarding one
            critical point, split by whitespace
        cp_type: Type of critical point. Currently, can be "atom", "bond", "ring", or "cage"
        conditionals (Dict[str, List[str]]): Parameters to extract with strings to search for to see if the
            data is present

    Returns:
        cp_dict: dictionary of critical point information
    """

    cp_dict: Dict[str, Any] = dict()
    unknown_id = 0

    cp_name = "null"

    this_conditionals = copy.deepcopy(conditionals)

    for ind, i in enumerate(lines_split):
        for k, v in this_conditionals.items():
            if all(x in i for x in v):
                if k == "cp_num":
                    cp_dict[k] = int(i[2][:-1])

                    # Placeholder name
                    # For nuclear critical points, this can be overwritten later
                    cp_name = f"{cp_dict[k]}_{cp_type}"

                elif k == "ele_info":
                    if i[2] == "Unknown":
                        cp_name = str(unknown_id) + "_Unknown"
                        cp_dict["number"] = "Unknown"
                        cp_dict["ele"] = "Unknown"
                        unknown_id += 1
                    else:
                        if len(i) == 3:
                            cp_dict["element"] = i[2].split("(")[1][:-1]
                            cp_dict["number"] = i[2].split("(")[0]
                        else:
                            cp_dict["element"] = i[2].split("(")[1]
                            cp_dict["number"] = i[2].split("(")[0]
                        cp_name = cp_dict["number"] + "_" + cp_dict["element"]

                elif k == "connected_bond_paths":
                        list_raw = [x for x in i[2:]]
                        # save only items that contain a number
                        list_raw = [
                            x for x in list_raw if any(char.isdigit() for char in x)
                        ]
                        list_raw = [int(x.split("(")[0]) for x in list_raw]
                        cp_dict[k] = list_raw

                elif k == "pos_ang":
                    cp_dict[k] = [float(x) for x in i[2:]]

                elif k == "esp_total":
                    cp_dict[k] = float(i[2])

                elif k == "eig_hess":
                    cp_dict[k] = np.sum(np.array([float(x) for x in i[-3:]]))

                elif k == "grad_norm" or k == "lap_norm":
                    cp_dict[k] = float(lines_split[ind + 2][-1])

                else:
                    cp_dict[k] = float(i[-1])

                this_conditionals.pop(k)
                break

    return cp_name, cp_dict


def parse_cp(lines: List[str]) -> Tuple[str, Dict[str, Any]]:
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
    cp_name = "null"
    
    # Figure out what kind of critical-point we're dealing with
    if "(3,-3)" in lines_split[0]:
        cp_type = "atom"
    elif "(3,-1)" in lines_split[0]:
        cp_type = "bond"
    elif "(3,+1)" in lines_split[0]:
        cp_type = "ring"
    elif "(3,+3)" in lines_split[0]:
        cp_type = "cage"

    else:
        return None, dict()

    conditionals = {
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

    atom_conditions = {k: v for k, v in conditionals.items() if k not in ["connected_bond_paths"]}
    bond_conditionals = {k: v for k, v in conditionals.items() if k not in ["ele_info"]}
    ring_cage_conditionals = {k: v for k, v in conditionals.items() if k not in ["connected_bond_paths", "ele_info"]}

    conditionals = {
        "atom": atom_conditionals,
        "bond": bond_conditionals,
        "ring": ring_cage_conditionals,
        "cage": ring_cage_conditionals
    }

    return extract_info_from_cp_text(lines_split, cp_type, conditionals[cp_type])


def get_qtaim_descs(file: Union[str, Path]) -> Dict[str, Dict[str, Any]]:
    """
    Parse CPprop file from multiwfn by parsing each individual critical-point section.

    Args:
        file (str): path to CPprop file
    
    Returns:
        ret_dict (Dict[str, Dict[str, Any]]): Output dictionary of QTAIM descriptors

    """

    cp_sections = list()
    ret_dict = dict()

    with open(file) as f:
        lines = f.readlines()
        lines = [line[:-1] for line in lines]

    # section lines into segments on ----------------
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
        ret_dict[cp_name] = cp_dict

    return ret_dict


def separate_cps_by_type(qtaim_descs: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[Any, Dict[str, Any]]]:
    """
    Separates QTAIM descriptors by type

    Args:
        qtaim_descs (Dict[str, Dict[str, Any]]): Dictionary where keys are CP names and values are dictionaries of
            descriptors obtained from `get_qtaim_descs` and `parse_cp`

    Returns:
        organized_descs (Dict[str, Dict[str, Dict[str, Any]]]): Dictionary of organized QTAIM critical points and their
            descriptors. Keys are "atom", "bond", "ring", and "cage", and values are dicts 
            {<CP name>: <QTAIM descriptors>}
    """

    organized_descs = {
        "atom": dict(),
        "bond": dict(),
        "ring": dict(),
        "cage": dict()
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
    atom_dict: Dict[str, Any],
    atom_cp_dict: Dict[str, Dict[str, Any]],
    margin: float = 0.5
) -> Tuple[Union[str, None], Dict]:
    """
    From a dictionary of atomic position and element symbol, find the corresponding cp in the atom_cp_dict
    
    Args:
        atom_dict (Dict[str, Any]): Dictionary containing the atom's element symbol and position in Cartesian space
        atom_cp_dict (Dict[str, Dict[str, Any]]): Dictionary where keys are critical point names and values are 
            descriptors, including element symbol and position
        margin (float): Maximum distance (in Angstrom) that a critical point can be away from an atom center
            and be associated with that atom. Default is 0.5 Angstrom

    Returns:
        cp_name (str | None): Key of atom_cp_dict; the name of the atom critical point corresponding to this atom. If no match
            is found, this will be None
        cp_dict (Dict): Dictionary of CP descriptors matching this atom
    """

    for cp_name, cp_dict in atom_cp_dict.items():
        if (
            int(cp_name.split("_")[0]) == atom_dict["ind"] + 1
            and cp_dict["element"] == atom_dict["element"]
        ):
            return cp_name, cp_dict

        else:
            if cp_dict["element"] == atom_dict["element"]:
                distance = np.linalg.norm(
                    np.array(cp_dict["pos_ang"]) - np.array(atom_dict["pos"])
                )

                if distance < margin:
                    return cp_name, cp_dict

    # No match
    return None, dict()


def map_atoms_cps(
    atoms: List[Dict[str, Any]],
    atom_cp_dict: Dict[str, Dict[str, Any]],
    margin: float = 0.5
) -> Tuple[
    Dict[int, Dict[str, Any]],
    Dict[int, Dict[str, Any]],
    List[int],
]:
    """
    Connect atom CPs to atoms by their positions.

    Takes:
        atoms (List[Dict[str, Any]]): List of dictionaries containing the atom's element symbols and positions in
            Cartesian space
        atom_cp_dict (Dict[str, Dict[str, Any]]): Dictionary where keys are critical point names and values are 
            descriptors, including element symbol and position
        margin (float): Maximum distance (in Angstrom) that a critical point can be away from an atom center
            and be associated with that atom. Default is 0.5 Angstrom
    Returns:
        index_to_cp_desc (Dict[int, Dict[str, Any]]): Dictionary mapping atomic indices to atom critical point descriptors
        index_to_cp_key (Dict[int, Dict[str, Any]]): dictionary with qtaim atoms as keys and dft atoms as values
        missing_atoms (List[int]): list of dft atoms that do not have a cp found in qtaim
    """

    index_to_cp_desc, index_to_cp_key = dict(), dict()
    missing_atoms = list()
    
    for atom_info in atoms:
        # finds cp by distance and naming scheme from CPprop.txt
        cp_name, this_atom_cp = match_atom_cp(
            atom_info, atom_cp_dict, margin=margin
        )
        
        # If this is False, that means no match was found
        if cp_name:
            index_to_cp_desc[atom_info["ind"]] = this_atom_cp
            index_to_cp_key[atom_info["ind"]] = {"key": cp_name, "pos": this_atom_cp["pos_ang"]}
        else:
            index_to_cp_desc[atom_info["ind"]] = dict()
            index_to_cp_key[atom_info["ind"]] = {"key": None, "pos": list()}
            missing_atoms.append(atom_info["ind"])

    return index_to_cp_desc, index_to_cp_key, missing_atoms


def match_bond_cp(
    atom_one: int,
    atom_two: int,
    bonds_cps: Dict[str, Dict[str, Any]]
) -> Optional[Dict[str, Any]]:
    """
    From two atom indices, find the corresponding bond cp (if possible)

    Takes:
        atom_one (int): First atom index
        atom_two (int): Second atom index
        bonds_cps: (Dict[str, Dict[str, Any]]): Dictionary where keys are critical point names and values are 
            descriptors
    Returns:
        Dict[str, Any] or None
    """

    for k, v in bonds_cps.items():
        if sorted(v["atom_inds"]) == sorted([atom_one, atom_two]):
            return v

    return None


def add_atoms(
    atom_info: Dict[int, Dict[str, Any]],
    organized_cps: Dict[str, Dict[str, Dict[str, Any]]],
    dist_threshold_bond: float = 1.0,
    dist_threshold_ring_cage: float = 3.0,
    distance_margin: float = 1.0,
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """
    Modify bond, ring, and cage CPs to include information about surrounding critical points. Bonds will include
    information about the atoms that make up the bond. Rings will include the names of neighboring bonds and the
    indices of the atoms involves, and cages will include information on neighboring rings, bonds, and the atoms
    involed.

    Some assumptions:
        - Bonds only involve two atoms
        - The atom CPs closest to a bond correspond to the atoms being bonded. Similarly, the closest bonds to a ring
            and the closest rings to a cage constitute those CPs.

    Args:
        atom_info (Dict[int, Dict[str, Any]]): Mapping between atom indices and atom properties (e.g. atomic positions)
        organized_cps (Dict[str, Dict[str, Dict[str, Any]]]): Keys are CP categories ("atom", "bond", "ring", and
            "cage"). Values are themselves dictionaries, where the keys are CP names and the values are CP descriptors
        dist_threshold_bond (float): If the nearest atoms to a bond CP are further from the bond CP than this threshold
            (default 1.0 Angstrom), then a warning will be raised.
        dist_threshold_ring_cage (float): If the nearest bond CPs to a ring CP or the nearest ring CPs to a cage CP are
            further than this threshold (default 3.0 Angstrom), then a warning will be raised.
        distance_margin (float): Rings must be made up of at least three bonds, and cages must be made up of at least
            three rings. We therefore define rings by the three closest bond CPs and cages by the three closest ring
            CPs. We associate additional bond/ring CPs with ring/cage CPs if they are no further than the third-closest
            bond/ring CP, plus this margin.

    Returns:
        modified_organized_cps (Dict[str, Dict[str, Dict[str, Any]]]): CP dict with additional information added.
    """

    def sort_cps_by_distance(
        position: np.ndarray,
        options: Dict[Any, Dict[str, Any]],
    ):
        dists = list()
        for oname, oval in options.items():
            if "pos_ang" in oval:
                opos = np.array(oval["pos_ang"])
            else:
                # For atom info
                opos = np.array(oval["pos"])

            dists.append((np.linalg.norm(position - opos), oname))

        return sorted(dists)

    if len(bond_cps) > 0 and len(atom_info) < 2:
        raise ValueError("Cannot have a bond CP with less than two atom CPs!")

    if len(ring_cps) > 0 and len(bond_cps) < 3:
        raise ValueError("Cannot have a ring CP with less than three bond CPs!")

    if len(cage_cps) > 0 and len(ring_cps) < 3:
        raise ValueError("Cannot have a cage CP with less than three ring CPs!")

    modified_organized_cps = copy.deepcopy(organizd_cps)

    bond_cps = organized_cps["bond"]
    ring_cps = organized_cps["ring"]
    cage_cps = organized_cps["cage"]

    # NOTE: for bonds, we associate based on atoms, NOT based on atom CPs
    for cp_name, cp_desc in bond_cps.items():
        sorted_atoms = sort_cps_by_distance(np.array(cp_desc["pos_ang"]), atom_info)

        if sorted_atoms[1][0] > dist_threshold_bond:
            warnings.warn("Warning: bond CP is far from bonding atoms")
        
        modified_organized_cps["bond"][cp_name]["atom_inds"] = sorted([ca[1] for ca in closest_atoms[:2]])

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
            for atom_ind in bond_cps[bond_name]["atom_inds"]:
                atom_inds.add(atom_ind)

        modified_organized_cps["ring"][cp_name]["bond_names"] = bond_names
        modified_organized_cps["ring"][cp_name]["atom_inds"] = list(atom_inds)

    for cp_name, cp_desc in cage_cps.items():
        sorted_rings = sort_cps_by_distance(np.array(cp_desc["pos_ang"]), ring_cps)
        max_close_dist = sorted_rings[2][0]

        # yell if the three closest bonds are further than the max distance
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
        bond_names = set()
        atom_inds = set()
        for ring_name in ring_names:
            for bond_name in ring_cps[ring_name]["bond_names"]:
                bond_names.add(bond_name)
            for atom_ind in ring_cps[ring_name]["atom_inds"]:
                atom_inds.add(atom_ind)

        modified_organized_cps["cage"][cp_name]["ring_names"] = ring_names
        modified_organized_cps["cage"][cp_name]["bond_names"] = list(bond_names)
        modified_organized_cps["cage"][cp_name]["atom_inds"] = list(atom_inds)

    return modified_organized_cps


def merge_qtaim_inds(
    qtaim_descs, dft_inp_file, bond_list=None, define_bonds="qtaim", margin=1.0
):
    """
    Gets mapping of qtaim indices to atom indices and remaps atom CP descriptors

    Takes
        qtaim_descs: dict of qtaim descriptors
        dft_inp_file: str input file for dft
    returns:
        dict of qtaim descriptors ordered by atoms in dft_inp_file
    """

    # open dft input file
    dft_dict = dft_inp_to_dict(dft_inp_file)
    # find only atom cps to map
    atom_only_cps, bond_cps = only_atom_cps(qtaim_descs)
    # remap qtaim indices to atom indices
    atom_cps_remapped, qtaim_to_dft, missing_atoms = find_cp_map(
        dft_dict, atom_only_cps, margin=margin
    )
    # remapping bonds
    bond_list_ret = []
    if define_bonds == "qtaim":
        bond_cps_qtaim = {}
        bond_cps = {
            k: v for k, v in bond_cps.items() if "connected_bond_paths" in v.keys()
        }
        #print("bond cps: ", bond_cps)

        for k, v in bond_cps.items():
            bond_list_unsorted = v["connected_bond_paths"]
            #print(bond_list_unsorted)
            bond_list_unsorted = [
                int(qtaim_to_dft[i - 1]["key"].split("_")[0]) - 1
                for i in bond_list_unsorted
            ]
            #print(bond_list_unsorted)
            bond_list_unsorted = sorted(bond_list_unsorted)
            #print(bond_list_unsorted)
            #assert len(bond_list_unsorted) == 2, "bond list not length 2"
            #assert bond_list_unsorted[0] != bond_list_unsorted[1], "bond list same"
            bond_cps_qtaim[tuple(bond_list_unsorted)] = v
            bond_list_ret.append(bond_list_unsorted)
        bond_cps = bond_cps_qtaim

    else:
        bond_cps = bond_cp_distance(bond_cps, bond_list, dft_dict, margin=margin)
    # merge dictionaries
    ret_dict = {**atom_cps_remapped, **bond_cps}
    return ret_dict
