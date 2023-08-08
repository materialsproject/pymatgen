from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from pymatgen.core.structure import Structure
from pymatgen.io.abinit.inputs import (
    BasicAbinitInput,
    BasicMultiDataset,
    ShiftMode,
    calc_shiftk,
    ebands_input,
    gs_input,
    ion_ioncell_relax_input,
    num_valence_electrons,
)
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

_test_dir = f"{TEST_FILES_DIR}/abinit"


def abiref_file(filename):
    """Return absolute path to filename in ~pymatgen/tests/files/abinit."""
    return os.path.join(_test_dir, filename)


def abiref_files(*filenames):
    """Return list of absolute paths to filenames in ~pymatgen/tests/files/abinit."""
    return [os.path.join(_test_dir, f) for f in filenames]


class AbinitInputTestCase(PymatgenTest):
    """Unit tests for BasicAbinitInput."""

    def test_api(self):
        """Testing BasicAbinitInput API."""
        # Build simple input with structure and pseudos
        unit_cell = {
            "acell": 3 * [10.217],
            "rprim": [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]],
            "ntypat": 1,
            "znucl": [14],
            "natom": 2,
            "typat": [1, 1],
            "xred": [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
        }

        inp = BasicAbinitInput(structure=unit_cell, pseudos=abiref_file("14si.pspnc"))

        shiftk = [[0.5, 0.5, 0.5], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]]
        assert_array_equal(calc_shiftk(inp.structure), shiftk)
        assert num_valence_electrons(inp.structure, inp.pseudos) == 8

        repr(inp), str(inp)
        assert len(inp) == 0
        assert not inp
        assert inp.get("foo", "bar") == "bar"
        assert inp.pop("foo", "bar") == "bar"
        assert inp.comment is None
        inp.set_comment("This is a comment")
        assert inp.comment == "This is a comment"
        assert inp.isnc
        assert not inp.ispaw

        inp["ecut"] = 1
        assert inp.get("ecut") == 1
        assert len(inp) == 1
        assert "ecut" in inp
        assert "foo" not in inp

        # Test to_string
        assert inp.to_str(with_structure=True, with_pseudos=True)
        assert inp.to_str(with_structure=False, with_pseudos=False)

        inp.set_vars(ecut=5, toldfe=1e-6)
        assert inp["ecut"] == 5
        inp.set_vars_ifnotin(ecut=-10)
        assert inp["ecut"] == 5

        _, tmpname = tempfile.mkstemp(text=True)
        inp.write(filepath=tmpname)

        # Cannot change structure variables directly.
        with pytest.raises(inp.Error):
            inp.set_vars(unit_cell)

        with pytest.raises(TypeError, match="type dict does not have `to_abivars` method"):
            inp.add_abiobjects({})

        with pytest.raises(KeyError) as exc:
            inp.remove_vars("foo", strict=True)
        assert "key='foo' not in self:" in str(exc.value)
        assert not inp.remove_vars("foo", strict=False)

        # Test deepcopy and remove_vars.
        inp["bdgw"] = [1, 2]
        inp_copy = inp.deepcopy()
        inp_copy["bdgw"][1] = 3
        assert inp["bdgw"] == [1, 2]
        assert inp.remove_vars("bdgw")
        assert "bdgw" not in inp

        removed = inp.pop_tolerances()
        assert len(removed) == 1
        assert removed["toldfe"] == 1e-6

        # Test set_spin_mode
        old_vars = inp.set_spin_mode("polarized")
        assert "nsppol" in inp
        assert inp["nspden"] == 2
        assert inp["nspinor"] == 1
        inp.set_vars(old_vars)

        # Test set_structure
        new_structure = inp.structure.copy()
        new_structure.perturb(distance=0.1)
        inp.set_structure(new_structure)
        assert inp.structure == new_structure

        # Compatible with Pickle and MSONable?
        self.serialize_with_pickle(inp, test_eq=False)

    def test_input_errors(self):
        """Testing typical BasicAbinitInput Error."""
        si_structure = Structure.from_file(abiref_file("si.cif"))

        # Ambiguous list of pseudos.
        with pytest.raises(BasicAbinitInput.Error):
            BasicAbinitInput(si_structure, pseudos=abiref_files("14si.pspnc", "14si.4.hgh"))

        # Pseudos do not match structure.
        with pytest.raises(BasicAbinitInput.Error):
            BasicAbinitInput(si_structure, pseudos=abiref_file("H-wdr.oncvpsp"))

        si1_negative_volume = {
            "ntypat": 1,
            "natom": 1,
            "typat": [1],
            "znucl": 14,
            "acell": 3 * [7.60],
            "rprim": [[0.0, 0.5, 0.5], [-0.5, -0.0, -0.5], [0.5, 0.5, 0.0]],
            "xred": [[0.0, 0.0, 0.0]],
        }

        # Negative triple product.
        with pytest.raises(BasicAbinitInput.Error):
            BasicAbinitInput(si1_negative_volume, pseudos=abiref_files("14si.pspnc"))

    def test_helper_functions(self):
        """Testing BasicAbinitInput helper functions."""
        inp = BasicAbinitInput(structure=abiref_file("si.cif"), pseudos="14si.pspnc", pseudo_dir=_test_dir)

        inp.set_kmesh(ngkpt=(1, 2, 3), shiftk=(1, 2, 3, 4, 5, 6))
        assert inp["kptopt"] == 1
        assert inp["nshiftk"] == 2

        inp.set_gamma_sampling()
        assert inp["kptopt"] == 1
        assert inp["nshiftk"] == 1
        assert np.all(inp["shiftk"] == 0)

        inp.set_kpath(ndivsm=3, kptbounds=None)
        assert inp["ndivsm"] == 3
        assert inp["iscf"] == -2
        assert len(inp["kptbounds"]) == 12


class TestMultiDataset(PymatgenTest):
    """Unit tests for BasicMultiDataset."""

    def test_api(self):
        """Testing BasicMultiDataset API."""
        structure = Structure.from_file(abiref_file("si.cif"))
        pseudo = abiref_file("14si.pspnc")
        pseudo_dir = os.path.dirname(pseudo)
        multi = BasicMultiDataset(structure=structure, pseudos=pseudo)
        with pytest.raises(ValueError, match="ndtset=-1 cannot be <=0"):
            BasicMultiDataset(structure=structure, pseudos=pseudo, ndtset=-1)

        multi = BasicMultiDataset(structure=structure, pseudos=pseudo, pseudo_dir=pseudo_dir)

        assert len(multi) == 1
        assert multi.ndtset == 1
        assert multi.isnc
        for i, inp in enumerate(multi):
            assert list(inp) == list(multi[i])

        multi.addnew_from(0)
        assert multi.ndtset == 2
        assert multi[0] is not multi[1]
        assert multi[0].structure == multi[1].structure
        assert multi[0].structure is not multi[1].structure

        multi.set_vars(ecut=2)
        assert all(inp["ecut"] == 2 for inp in multi)
        assert multi.get("ecut") == [2, 2]

        multi[1].set_vars(ecut=1)
        assert multi[0]["ecut"] == 2
        assert multi[1]["ecut"] == 1
        assert multi.get("ecut") == [2, 1]

        assert multi.get("foo", "default") == ["default", "default"]

        multi[1].set_vars(paral_kgb=1)
        assert "paral_kgb" not in multi[0]
        assert multi.get("paral_kgb") == [None, 1]

        pert_structure = structure.copy()
        pert_structure.perturb(distance=0.1)
        assert structure != pert_structure

        assert multi.set_structure(structure) == multi.ndtset * [structure]
        assert all(s == structure for s in multi.structure)
        assert multi.has_same_structures
        multi[1].set_structure(pert_structure)
        assert multi[0].structure != multi[1].structure
        assert multi[1].structure == pert_structure
        assert not multi.has_same_structures

        split = multi.split_datasets()
        assert len(split) == 2
        assert all(split[i] == multi[i] for i in range(multi.ndtset))
        repr(multi)
        str(multi)
        assert multi.to_str(with_pseudos=False)

        tmpdir = tempfile.mkdtemp()
        filepath = f"{tmpdir}/run.abi"
        inp.write(filepath=filepath)
        multi.write(filepath=filepath)

        new_multi = BasicMultiDataset.from_inputs(list(multi))
        assert new_multi.ndtset == multi.ndtset
        assert new_multi.structure == multi.structure

        for old_inp, new_inp in zip(multi, new_multi):
            assert old_inp is not new_inp
            assert old_inp.as_dict() == new_inp.as_dict()

        ref_input = multi[0]
        new_multi = BasicMultiDataset.replicate_input(input=ref_input, ndtset=4)

        assert new_multi.ndtset == 4
        for inp in new_multi:
            assert ref_input is not inp
            assert ref_input.as_dict() == inp.as_dict()

        # Compatible with Pickle and MSONable?
        self.serialize_with_pickle(multi, test_eq=False)


class TestShiftMode(PymatgenTest):
    def test_shiftmode(self):
        """Testing shiftmode."""
        gamma = ShiftMode.GammaCentered
        assert ShiftMode.from_object("G") == gamma
        assert ShiftMode.from_object(gamma) == gamma
        with pytest.raises(TypeError, match="The object provided is not handled: type dict"):
            ShiftMode.from_object({})


class TestFactory(PymatgenTest):
    def setUp(self):
        # Si ebands
        self.si_structure = Structure.from_file(abiref_file("si.cif"))
        self.si_pseudo = abiref_file("14si.pspnc")

    def test_gs_input(self):
        """Testing gs_input factory."""
        inp = gs_input(self.si_structure, self.si_pseudo, kppa=10, ecut=10, spin_mode="polarized")
        str(inp)
        assert inp["nsppol"] == 2
        assert inp["nband"] == 14
        assert_array_equal(inp["ngkpt"], [2, 2, 2])

    def test_ebands_input(self):
        """Testing ebands_input factory."""
        multi = ebands_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)
        str(multi)

        scf_inp, nscf_inp = multi.split_datasets()

        # Test dos_kppa and other options.
        multi_dos = ebands_input(
            self.si_structure,
            self.si_pseudo,
            nscf_nband=10,
            kppa=10,
            ecut=2,
            spin_mode="unpolarized",
            smearing=None,
            charge=2.0,
            dos_kppa=50,
        )
        assert len(multi_dos) == 3
        assert all(i["charge"] == 2 for i in multi_dos)
        assert multi_dos.get("nsppol") == [1, 1, 1]
        assert multi_dos.get("iscf") == [None, -2, -2]

        multi_dos = ebands_input(
            self.si_structure,
            self.si_pseudo,
            nscf_nband=10,
            kppa=10,
            ecut=2,
            spin_mode="unpolarized",
            smearing=None,
            charge=2.0,
            dos_kppa=[50, 100],
        )
        assert len(multi_dos) == 4
        assert multi_dos.get("iscf") == [None, -2, -2, -2]
        str(multi_dos)

    def test_ion_ioncell_relax_input(self):
        """Testing ion_ioncell_relax_input factory."""
        multi = ion_ioncell_relax_input(self.si_structure, self.si_pseudo, kppa=10, ecut=2)

        str(multi)
        ion_inp, ioncell_inp = multi.split_datasets()
        assert ion_inp["chksymbreak"] == 0
        assert ion_inp["ionmov"] == 3
        assert ion_inp["optcell"] == 0
        assert ioncell_inp["ionmov"] == 3
        assert ioncell_inp["optcell"] == 2
