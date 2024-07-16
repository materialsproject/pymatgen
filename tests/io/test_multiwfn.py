from __future__ import annotations

import numpy as np
import pytest
from monty.json import MontyDecoder, jsanitize

from pymatgen.core.structure import Molecule
from pymatgen.io.multiwfn import(
    QTAIM_CONDITIONALS,
    extract_info_from_cp_text,
    parse_cp,
    get_qtaim_descs,
    separate_cps_by_type,
    match_atom_cp,
    map_atoms_cps,
    add_atoms,
    process_multiwfn_qtaim,
)
from pymatgen.util.testing import TEST_FILES_DIR


base_dir = TEST_FILES_DIR / "io" / "multiwfn"


def test_parse_single_cp():
    
    # Test that extract_info_from_cp_text behaves as expected with parse_cp
    # Also tests atom parsing
    with open(base_dir / "cp_just_atom.txt") as file:
        contents = file.readlines()
        name1, desc1 = parse_cp(contents)

        contents_split = [line.split() for line in contents]
        conditionals = {k: v for k, v in QTAIM_CONDITIONALS.items() if k not in ["connected_bond_paths"]}
        name2, desc2 = extract_info_from_cp_text(contents_split, "atom", conditionals)

        assert name1 == name2
        for k, v in desc1.items():
            assert desc2.get(k) == pytest.approx(v)

        assert name1 == "2_N"
        assert desc1["cp_num"] == 102
        assert desc1["element"] == "N"
        # TODO: should we be returning this as an integer?
        assert desc1["number"] == "2"
        assert desc1["pos_ang"] == pytest.approx([1.626104042116, -1.859508691395, -0.405402516863])
        assert desc1["density_total"] == pytest.approx(183.1401128)
        assert "connected_bond_paths" not in desc1

    # Test atom parsing with CP not associated with a known nucleus
    with open(base_dir / "cp_unknown_atom.txt") as file:
        contents = file.readlines()
        name, desc = parse_cp(contents)

        assert name == '142_Unknown'
        assert desc["cp_num"] == 142
        assert desc["number"] == "Unknown"
        assert desc["ele"] == "Unknown"
        assert desc["density_alpha"] == pytest.approx(8.066360869)
        assert desc["density_alpha"] == pytest.approx(desc["density_beta"]) and desc['spin_density'] == pytest.approx(0.0)

    # Test bond parsing
    with open(base_dir / "cp_just_bond.txt") as file:
        contents = file.readlines()
        name, desc = parse_cp(contents)

        assert name == "121_bond"
        assert "ele_info" not in desc
        assert desc["Lagrangian_K"] == pytest.approx(426.4263555)
        assert desc["Hamiltonian_K"] == pytest.approx(106.7023631)
        assert desc["energy_density"] == pytest.approx(-106.7023631)
        assert desc["lap_e_density"] == pytest.approx(1278.89597)

    # Test ring parsing
    with open(base_dir / "cp_just_ring.txt") as file:
        contents = file.readlines()
        name, desc = parse_cp(contents)

        assert name == "123_ring"
        assert "connected_bond_paths" not in desc and "ele_info" not in desc
        assert desc["e_loc_func"] == pytest.approx(0.4201012445)
        assert desc["lol"] == pytest.approx(0.4597922949)
        assert desc["ave_loc_ion_E"] == pytest.approx(6.928709119)
        assert desc["delta_g_promolecular"] == pytest.approx(0.001716120125)
        assert desc["delta_g_hirsh"] == pytest.approx(0.003153281621)
        assert desc["esp_nuc"] == pytest.approx(176.0405167)
        assert desc["esp_e"] == pytest.approx(-79.88321676)
        assert desc["esp_total"] == pytest.approx(96.1572999)

    # Test cage parsing
    with open(base_dir / "cp_just_cage.txt") as file:
        contents = file.readlines()
        name, desc = parse_cp(contents)

        assert name == "56_cage"
        assert "connected_bond_paths" not in desc and "ele_info" not in desc
        assert desc["grad_norm"] == pytest.approx(7.920799975e-18)
        assert desc["lap_norm"] == pytest.approx(0.01583124127)
        assert desc["eig_hess"] == pytest.approx(0.0158312412724)
        assert desc["det_hessian"] == pytest.approx(3.943311116e-08)
        assert desc["ellip_e_dens"] == pytest.approx(-0.759846)
        assert desc["eta"] == pytest.approx(0.083769)

    # Test parsing with unknown/improper CP type
    with open(base_dir / "cp_fake_type.txt") as file:
        contents = file.readlines()
        name, desc = parse_cp(contents)
        assert name is None
        assert len(desc) == 0


def test_parse_cps():
    # Make sure we're not missing any CPs
    all_descs = get_qtaim_descs(base_dir / "CPprop_all.txt")
    assert len(all_descs) == 181

    nums = [int(v["cp_num"]) for v in all_descs.values()]
    for i in range(1, 182):
        assert i in nums

    # Test separation by type
    separated = separate_cps_by_type(all_descs)
    # NOTE: this does not sum to 181, because there are four atom CPs associated with unknown nuclei
    # These "Unknown" CPs are excluded at this point
    assert len(separated["atom"]) == 56
    assert len(separated["bond"]) == 90
    assert len(separated["ring"]) == 28
    assert len(separated["cage"]) == 3


def test_atom_matching():
    mol = Molecule.from_file(base_dir / "mol_all.xyz")

    all_descs = get_qtaim_descs(base_dir / "CPprop_all.txt")
    separated = separate_cps_by_type(all_descs)

    all_descs_fudged = get_qtaim_descs(base_dir / "CPprop_fudged_nuclei.txt")
    separated_fudged = separate_cps_by_type(all_descs_fudged)

    # Test successful single match
    name, desc = match_atom_cp(mol, 0, separated["atom"])
    assert name == "1_Sm"
    assert desc["element"] == "Sm"
    assert desc["number"] == "1"

    # Test successful match by distance
    name, desc = match_atom_cp(mol, 0, separated_fudged["atom"])
    assert name == "78_Sm"

    # Test unsuccessful single match
    name, desc = match_atom_cp(mol, 55, separated_fudged["atom"])
    assert name is None
    assert len(desc) == 0

    # Test overall mapping
    mapping, missing = map_atoms_cps(mol, separated_fudged["atom"])
    assert len(mapping) == 56
    assert mapping[0]["element"] == "Sm"
    assert len(mapping[55]) == 0

    assert len(missing) == 1
    assert missing[0] == 55


def test_add_atoms():
    pass


def test_process_multiwfn_qtaim():
    pass