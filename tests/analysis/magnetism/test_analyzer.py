from __future__ import annotations

import unittest
import warnings
from shutil import which

import numpy as np
import pytest
from monty.serialization import loadfn
from pytest import approx

from pymatgen.analysis.magnetism import (
    CollinearMagneticStructureAnalyzer,
    MagneticStructureEnumerator,
    Ordering,
    magnetic_deformation,
)
from pymatgen.core import Element, Lattice, Species, Structure
from pymatgen.io.cif import CifParser
from pymatgen.util.testing import TEST_FILES_DIR

enum_cmd = which("enum.x") or which("multienum.x")
makestr_cmd = which("makestr.x") or which("makeStr.x") or which("makeStr.py")
enumlib_present = enum_cmd and makestr_cmd


class TestCollinearMagneticStructureAnalyzer(unittest.TestCase):
    def setUp(self):
        parser = CifParser(f"{TEST_FILES_DIR}/Fe.cif")
        self.Fe = parser.get_structures()[0]

        parser = CifParser(f"{TEST_FILES_DIR}/LiFePO4.cif")
        self.LiFePO4 = parser.get_structures()[0]

        parser = CifParser(f"{TEST_FILES_DIR}/Fe3O4.cif")
        self.Fe3O4 = parser.get_structures()[0]

        parser = CifParser(f"{TEST_FILES_DIR}/magnetic.ncl.example.GdB4.mcif")
        self.GdB4 = parser.get_structures()[0]

        parser = CifParser(f"{TEST_FILES_DIR}/magnetic.example.NiO.mcif")
        self.NiO_expt = parser.get_structures()[0]

        latt = Lattice.cubic(4.17)
        species = ["Ni", "O"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        self.NiO = Structure.from_spacegroup(225, latt, species, coords)

        latt = Lattice([[2.085, 2.085, 0.0], [0.0, -2.085, -2.085], [-2.085, 2.085, -4.17]])
        species = ["Ni", "Ni", "O", "O"]
        coords = [[0.5, 0, 0.5], [0, 0, 0], [0.25, 0.5, 0.25], [0.75, 0.5, 0.75]]
        self.NiO_AFM_111 = Structure(latt, species, coords, site_properties={"magmom": [-5, 5, 0, 0]})

        latt = Lattice([[2.085, 2.085, 0], [0, 0, -4.17], [-2.085, 2.085, 0]])
        species = ["Ni", "Ni", "O", "O"]
        coords = [[0.5, 0.5, 0.5], [0, 0, 0], [0, 0.5, 0], [0.5, 0, 0.5]]
        self.NiO_AFM_001 = Structure(latt, species, coords, site_properties={"magmom": [-5, 5, 0, 0]})

        latt = Lattice([[2.085, 2.085, 0], [0, 0, -4.17], [-2.085, 2.085, 0]])
        species = ["Ni", "Ni", "O", "O"]
        coords = [[0.5, 0.5, 0.5], [0, 0, 0], [0, 0.5, 0], [0.5, 0, 0.5]]
        self.NiO_AFM_001_opposite = Structure(latt, species, coords, site_properties={"magmom": [5, -5, 0, 0]})

        latt = Lattice([[2.085, 2.085, 0], [0, 0, -4.17], [-2.085, 2.085, 0]])
        species = ["Ni", "Ni", "O", "O"]
        coords = [[0.5, 0.5, 0.5], [0, 0, 0], [0, 0.5, 0], [0.5, 0, 0.5]]
        self.NiO_unphysical = Structure(latt, species, coords, site_properties={"magmom": [-3, 0, 0, 0]})

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_get_representations(self):
        # tests to convert between storing magnetic moment information
        # on site_properties or on Species 'spin' property

        # test we store magnetic moments on site properties
        self.Fe.add_site_property("magmom", [5])
        msa = CollinearMagneticStructureAnalyzer(self.Fe)
        assert msa.structure.site_properties["magmom"][0] == 5

        # and that we can retrieve a spin representation
        Fe_spin = msa.get_structure_with_spin()
        assert "magmom" not in Fe_spin.site_properties
        assert Fe_spin[0].specie.spin == 5

        # test we can remove magnetic moment information
        msa.get_nonmagnetic_structure()
        assert "magmom" not in Fe_spin.site_properties

        # test with disorder on magnetic site
        self.Fe[0] = {
            Species("Fe", oxidation_state=0, properties={"spin": 5}): 0.5,
            "Ni": 0.5,
        }
        with pytest.raises(
            NotImplementedError,
            match="CollinearMagneticStructureAnalyzer not implemented for disordered structures,"
            " make ordered approximation first.",
        ):
            CollinearMagneticStructureAnalyzer(self.Fe)

    def test_matches(self):
        assert self.NiO.matches(self.NiO_AFM_111)
        assert self.NiO.matches(self.NiO_AFM_001)

        # MSA adds magmoms to Structure, so not equal
        msa = CollinearMagneticStructureAnalyzer(self.NiO, overwrite_magmom_mode="replace_all")
        assert not msa.matches_ordering(self.NiO)
        assert not msa.matches_ordering(self.NiO_AFM_111)
        assert not msa.matches_ordering(self.NiO_AFM_001)

        msa = CollinearMagneticStructureAnalyzer(self.NiO_AFM_001, overwrite_magmom_mode="respect_sign")
        assert not msa.matches_ordering(self.NiO)
        assert not msa.matches_ordering(self.NiO_AFM_111)
        assert msa.matches_ordering(self.NiO_AFM_001)
        assert msa.matches_ordering(self.NiO_AFM_001_opposite)

        msa = CollinearMagneticStructureAnalyzer(self.NiO_AFM_111, overwrite_magmom_mode="respect_sign")
        assert not msa.matches_ordering(self.NiO)
        assert msa.matches_ordering(self.NiO_AFM_111)
        assert not msa.matches_ordering(self.NiO_AFM_001)
        assert not msa.matches_ordering(self.NiO_AFM_001_opposite)

    def test_modes(self):
        mode = "none"
        msa = CollinearMagneticStructureAnalyzer(self.NiO, overwrite_magmom_mode=mode)
        magmoms = msa.structure.site_properties["magmom"]
        assert magmoms == [0, 0]

        mode = "respect_sign"
        msa = CollinearMagneticStructureAnalyzer(self.NiO_unphysical, overwrite_magmom_mode=mode)
        magmoms = msa.structure.site_properties["magmom"]
        assert magmoms == [-5, 0, 0, 0]

        mode = "respect_zeros"
        msa = CollinearMagneticStructureAnalyzer(self.NiO_unphysical, overwrite_magmom_mode=mode)
        magmoms = msa.structure.site_properties["magmom"]
        assert magmoms == [5, 0, 0, 0]

        mode = "replace_all"
        msa = CollinearMagneticStructureAnalyzer(self.NiO_unphysical, overwrite_magmom_mode=mode, make_primitive=False)
        magmoms = msa.structure.site_properties["magmom"]
        assert magmoms == [5, 5, 0, 0]

        mode = "replace_all_if_undefined"
        msa = CollinearMagneticStructureAnalyzer(self.NiO, overwrite_magmom_mode=mode)
        magmoms = msa.structure.site_properties["magmom"]
        assert magmoms == [5, 0]

        mode = "normalize"
        msa = CollinearMagneticStructureAnalyzer(msa.structure, overwrite_magmom_mode="normalize")
        magmoms = msa.structure.site_properties["magmom"]
        assert magmoms == [1, 0]

    def test_net_positive(self):
        msa = CollinearMagneticStructureAnalyzer(self.NiO_unphysical)
        magmoms = msa.structure.site_properties["magmom"]
        assert magmoms == [3, 0, 0, 0]

    def test_get_ferromagnetic_structure(self):
        msa = CollinearMagneticStructureAnalyzer(self.NiO, overwrite_magmom_mode="replace_all_if_undefined")
        s1 = msa.get_ferromagnetic_structure()
        s1_magmoms = [float(m) for m in s1.site_properties["magmom"]]
        s1_magmoms_ref = [5.0, 0.0]
        assert s1_magmoms == s1_magmoms_ref

        _ = CollinearMagneticStructureAnalyzer(self.NiO_AFM_111, overwrite_magmom_mode="replace_all_if_undefined")
        s2 = msa.get_ferromagnetic_structure(make_primitive=False)
        s2_magmoms = [float(m) for m in s2.site_properties["magmom"]]
        s2_magmoms_ref = [5.0, 0.0]
        assert s2_magmoms == s2_magmoms_ref

        s2_prim = msa.get_ferromagnetic_structure(make_primitive=True)
        assert CollinearMagneticStructureAnalyzer(s1).matches_ordering(s2_prim)

    def test_magnetic_properties(self):
        msa = CollinearMagneticStructureAnalyzer(self.GdB4)
        assert not msa.is_collinear

        msa = CollinearMagneticStructureAnalyzer(self.Fe)
        assert not msa.is_magnetic

        self.Fe.add_site_property("magmom", [5])

        msa = CollinearMagneticStructureAnalyzer(self.Fe)
        assert msa.is_magnetic
        assert msa.is_collinear
        assert msa.ordering == Ordering.FM

        msa = CollinearMagneticStructureAnalyzer(
            self.NiO,
            make_primitive=False,
            overwrite_magmom_mode="replace_all_if_undefined",
        )
        assert msa.number_of_magnetic_sites == 4
        assert msa.number_of_unique_magnetic_sites() == 1
        assert msa.types_of_magnetic_species == (Element.Ni,)
        assert msa.get_exchange_group_info() == ("Fm-3m", 225)

    def test_str(self):
        msa = CollinearMagneticStructureAnalyzer(self.NiO_AFM_001)

        ref_msa_str = """Structure Summary
Lattice
    abc : 2.948635277547903 4.17 2.948635277547903
 angles : 90.0 90.0 90.0
 volume : 36.2558565
      A : 2.085 2.085 0.0
      B : 0.0 0.0 -4.17
      C : -2.085 2.085 0.0
Magmoms Sites
+5.00   PeriodicSite: Ni (0.0, 0.0, 0.0) [0.0, 0.0, 0.0]
        PeriodicSite: O (0.0, 0.0, -2.085) [0.0, 0.5, 0.0]
        PeriodicSite: O (0.0, 2.085, 0.0) [0.5, 0.0, 0.5]
-5.00   PeriodicSite: Ni (0.0, 2.085, -2.085) [0.5, 0.5, 0.5]"""

        # just compare lines form 'Magmoms Sites',
        # since lattice param string can vary based on machine precision
        assert "\n".join(str(msa).split("\n")[-5:-1]) == "\n".join(ref_msa_str.split("\n")[-5:-1])

    def test_round_magmoms(self):
        struct = self.NiO_AFM_001.copy()
        struct.add_site_property("magmom", [-5.0143, -5.02, 0.147, 0.146])

        msa = CollinearMagneticStructureAnalyzer(struct, round_magmoms=0.001, make_primitive=False)
        assert np.allclose(msa.magmoms, [5.0171, 5.0171, -0.1465, -0.1465])
        assert msa.magnetic_species_and_magmoms["Ni"] == approx(5.0171)
        assert msa.magnetic_species_and_magmoms["O"] == approx(0.1465)

        struct.add_site_property("magmom", [-5.0143, 4.5, 0.147, 0.146])
        msa = CollinearMagneticStructureAnalyzer(struct, round_magmoms=0.001, make_primitive=False)
        assert np.allclose(msa.magmoms, [5.0143, -4.5, -0.1465, -0.1465])
        assert msa.magnetic_species_and_magmoms["Ni"][0] == approx(4.5)
        assert msa.magnetic_species_and_magmoms["Ni"][1] == approx(5.0143)
        assert msa.magnetic_species_and_magmoms["O"] == approx(0.1465)

    def test_missing_spin(self):
        # This test catches the case where a structure has some species with
        # Species.spin=None. This previously raised an error upon construction
        # of the analyzer).
        latt = Lattice([[2.085, 2.085, 0.0], [0.0, -2.085, -2.085], [-2.085, 2.085, -4.17]])
        species = [Species("Ni", spin=-5), Species("Ni", spin=5), Species("O", spin=None), Species("O", spin=None)]
        coords = [[0.5, 0, 0.5], [0, 0, 0], [0.25, 0.5, 0.25], [0.75, 0.5, 0.75]]
        struct = Structure(latt, species, coords)

        msa = CollinearMagneticStructureAnalyzer(struct, round_magmoms=0.001, make_primitive=False)
        assert msa.structure.site_properties["magmom"] == [-5, 5, 0, 0]


class TestMagneticStructureEnumerator(unittest.TestCase):
    @unittest.skipIf(not enumlib_present, "enumlib not present")
    def test_ordering_enumeration(self):
        # simple afm
        structure = Structure.from_file(f"{TEST_FILES_DIR}/magnetic_orderings/LaMnO3.json")
        enumerator = MagneticStructureEnumerator(structure)
        assert enumerator.input_origin == "afm"

        # ferrimagnetic (Cr produces net spin)
        structure = Structure.from_file(f"{TEST_FILES_DIR}/magnetic_orderings/Cr2NiO4.json")
        enumerator = MagneticStructureEnumerator(structure)
        assert enumerator.input_origin == "ferri_by_Cr"

        # antiferromagnetic on single magnetic site
        structure = Structure.from_file(f"{TEST_FILES_DIR}/magnetic_orderings/Cr2WO6.json")
        enumerator = MagneticStructureEnumerator(structure)
        assert enumerator.input_origin == "afm_by_Cr"

        # afm requiring large cell size
        # (enable for further development of workflow, too slow for CI)

        # structure = Structure.from_file(f"{ref_dir}/CuO.json")
        # enumerator = MagneticOrderingsenumerator(
        #     structure, default_magmoms={"Cu": 1.73}, transformation_kwargs={"max_cell_size": 4}
        # )
        # assert enumerator.input_origin == "afm"

        # antiferromagnetic by structural motif
        structure = Structure.from_file(f"{TEST_FILES_DIR}/magnetic_orderings/Ca3Co2O6.json")
        enumerator = MagneticStructureEnumerator(
            structure,
            strategies=("antiferromagnetic_by_motif",),
            # this example just misses default cut-off, so do not truncate
            truncate_by_symmetry=False,
            transformation_kwargs={"max_cell_size": 2},
        )
        assert enumerator.input_origin == "afm_by_motif_2a"


class TestMagneticDeformation(unittest.TestCase):
    def test_magnetic_deformation(self):
        test_structs = loadfn(f"{TEST_FILES_DIR}/magnetic_deformation.json")
        mag_def = magnetic_deformation(test_structs[0], test_structs[1])

        assert mag_def.type == "NM-FM"
        assert mag_def.deformation == approx(5.0130859485170971)
