from __future__ import annotations

import copy

import pytest
from pymatgen.core.structure import Molecule
from pymatgen.io.multiwfn import (
    QTAIM_CONDITIONALS,
    add_atoms,
    extract_info_from_cp_text,
    get_qtaim_descs,
    map_atoms_cps,
    match_atom_cp,
    parse_cp,
    process_multiwfn_qtaim,
    separate_cps_by_type,
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

        assert name == "142_Unknown"
        assert desc["cp_num"] == 142
        assert desc["number"] == "Unknown"
        assert desc["ele"] == "Unknown"
        assert desc["density_alpha"] == pytest.approx(8.066360869)
        assert desc["density_alpha"] == pytest.approx(desc["density_beta"])
        assert desc["spin_density"] == pytest.approx(0.0)

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
        assert "connected_bond_paths" not in desc
        assert "ele_info" not in desc
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
        assert "connected_bond_paths" not in desc
        assert "ele_info" not in desc
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
    mol = Molecule.from_file(base_dir / "mol_all.xyz")

    all_descs = get_qtaim_descs(base_dir / "CPprop_all.txt")
    separated = separate_cps_by_type(all_descs)

    # Test ValueErrors
    mol_minatom = Molecule(["O"], [[0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match=r"bond CP"):
        add_atoms(mol_minatom, separated)

    sep_minbonds = copy.deepcopy(separated)
    sep_minbonds["bond"] = {k: separated["bond"][k] for k in ["1_bond", "2_bond"]}

    with pytest.raises(ValueError, match=r"ring CP"):
        add_atoms(mol, sep_minbonds)

    sep_minrings = copy.deepcopy(separated)
    sep_minrings["ring"] = {k: separated["ring"][k] for k in ["13_ring", "14_ring"]}

    with pytest.raises(ValueError, match=r"cage CP"):
        add_atoms(mol, sep_minrings)

    # Test distance-based metric
    modified = add_atoms(mol, separated, bond_atom_criterion="distance")

    # Test that atom indices are being connected reasonably to bonds
    assert sorted(modified["bond"]["1_bond"]["atom_inds"]) == [3, 14]

    # Test that bonds and atoms are being connected reasonably to rings
    assert sorted(modified["ring"]["13_ring"]["atom_inds"]) == [35, 36, 37, 38, 39, 40, 42, 46]
    assert sorted(modified["ring"]["13_ring"]["bond_names"]) == [
        "11_bond",
        "23_bond",
        "27_bond",
        "28_bond",
        "3_bond",
        "5_bond",
        "8_bond",
    ]

    # Test that rings, bonds, and atoms are being connected reasonably to cages
    assert sorted(modified["cage"]["67_cage"]["atom_inds"]) == [0, 20, 22, 23, 24, 25, 27, 50, 51, 52, 55]
    assert sorted(modified["cage"]["67_cage"]["bond_names"]) == [
        "100_bond",
        "121_bond",
        "134_bond",
        "143_bond",
        "169_bond",
        "171_bond",
        "180_bond",
        "53_bond",
        "55_bond",
        "56_bond",
        "58_bond",
        "60_bond",
        "71_bond",
        "72_bond",
        "74_bond",
        "76_bond",
        "79_bond",
        "81_bond",
        "82_bond",
        "94_bond",
        "95_bond",
    ]
    assert sorted(modified["cage"]["67_cage"]["ring_names"]) == ["62_ring", "66_ring", "70_ring"]

    # Test with QTAIM-defined bonds
    remapped_atoms, _ = map_atoms_cps(mol, separated["atom"])
    separated["atom"] = remapped_atoms
    modified_qtaim = add_atoms(mol, separated, bond_atom_criterion="qtaim")

    assert len(modified_qtaim["bond"]) == 63
    assert modified_qtaim["bond"]["9_bond"]["atom_inds"] == [2, 43]

    # Test with combined QTAIM- + distance-based bonds
    modified_separated = add_atoms(mol, separated)

    assert modified_separated["bond"]["9_bond"]["atom_inds"] == [2, 43]
    assert len(modified_separated["bond"]) == 90


def test_process_multiwfn_qtaim():
    # Don't need to test very thoroughly, since we've already tested everything else
    mol = Molecule.from_file(base_dir / "mol_all.xyz")

    descriptors = process_multiwfn_qtaim(mol, base_dir / "CPprop_all.txt")

    # Checking that everything's been parsed and separated properly
    assert len(descriptors["atom"]) == 56
    assert len(descriptors["bond"]) == 90
    assert len(descriptors["ring"]) == 28
    assert len(descriptors["cage"]) == 3

    # Checking that atoms have been remapped properly
    for i in range(1, 56):
        assert i in descriptors["atom"]

    # Checking that atom info has been added
    assert sorted(descriptors["bond"]["1_bond"]["atom_inds"]) == [3, 14]

    assert sorted(descriptors["ring"]["13_ring"]["atom_inds"]) == [35, 36, 37, 38, 39, 40, 42, 46]
    assert sorted(descriptors["ring"]["13_ring"]["bond_names"]) == [
        "11_bond",
        "23_bond",
        "27_bond",
        "28_bond",
        "3_bond",
        "5_bond",
        "8_bond",
    ]
    assert sorted(descriptors["cage"]["67_cage"]["atom_inds"]) == [0, 20, 22, 23, 24, 25, 27, 50, 51, 52, 55]
    assert sorted(descriptors["cage"]["67_cage"]["bond_names"]) == [
        "100_bond",
        "121_bond",
        "134_bond",
        "143_bond",
        "169_bond",
        "171_bond",
        "180_bond",
        "53_bond",
        "55_bond",
        "56_bond",
        "58_bond",
        "60_bond",
        "71_bond",
        "72_bond",
        "74_bond",
        "76_bond",
        "79_bond",
        "81_bond",
        "82_bond",
        "94_bond",
        "95_bond",
    ]
    assert sorted(descriptors["cage"]["67_cage"]["ring_names"]) == ["62_ring", "66_ring", "70_ring"]
