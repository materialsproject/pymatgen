# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest
import numpy as np

from pymatgen import Structure
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinit.inputs import BasicAbinitInput

_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                         'test_files', "abinit")


def abiref_file(filename):
    return os.path.join(_test_dir, filename)


def abiref_files(filenames):
    return [os.path.join(_test_dir, f) for f in filenames]


class AbinitInputTestCase(PymatgenTest):
    """Unit tests for BasicAbinitInput."""

    def test_api(self):
        """Testing BasicAbinitInput API."""
        # Build simple input with structure and pseudos
        unit_cell = {
            "acell": 3*[10.217],
            'rprim': [[.0, .5, .5],
                      [.5, .0, .5],
                      [.5, .5, .0]],
            'ntypat': 1,
            'znucl': [14,],
            'natom': 2,
            'typat': [1, 1],
            'xred': [[.0, .0, .0],
                     [.25,.25,.25]]
        }

        inp = BasicAbinitInput(structure=unit_cell, pseudos=abiref_file("14si.pspnc"))

        repr(inp), str(inp)
        assert len(inp) == 0 and not inp
        assert inp.get("foo", "bar") == "bar" and inp.pop("foo", "bar") == "bar"
        assert inp.comment is None
        inp.set_comment("This is a comment")
        assert inp.comment == "This is a comment"
        assert inp.isnc and not inp.ispaw

        inp["ecut" ] = 1
        assert inp.get("ecut") == 1 and len(inp) == 1 and "ecut" in inp.keys() and "foo" not in inp

        # Test to_string
        assert inp.to_string(with_structure=True, with_pseudos=True)
        assert inp.to_string(with_structure=False, with_pseudos=False)

        inp.set_vars(ecut=5, toldfe=1e-6)
        assert inp["ecut"] == 5
        inp.set_vars_ifnotin(ecut=-10)
        assert inp["ecut"] == 5

        import tempfile
        _, tmpname = tempfile.mkstemp(text=True)
        inp.write(filepath=tmpname)

        # Cannot change structure variables directly.
        with self.assertRaises(inp.Error):
            inp.set_vars(unit_cell)

        with self.assertRaises(TypeError):
            inp.add_abiobjects({})

        with self.assertRaises(KeyError):
            inp.remove_vars("foo", strict=True)
        assert not inp.remove_vars("foo", strict=False)

        # Test deepcopy and remove_vars.
        inp["bdgw"] = [1, 2]
        inp_copy = inp.deepcopy()
        inp_copy["bdgw"][1] = 3
        assert inp["bdgw"] == [1, 2]
        assert inp.remove_vars("bdgw") and "bdgw" not in inp

        removed = inp.pop_tolerances()
        assert len(removed) == 1 and removed["toldfe"] == 1e-6

        # Test set_spin_mode
        old_vars = inp.set_spin_mode("polarized")
        assert "nsppol" in inp and inp["nspden"] == 2 and inp["nspinor"] == 1
        inp.set_vars(old_vars)

        # Test set_structure
        new_structure = inp.structure.copy()
        new_structure.perturb(distance=0.1)
        inp.set_structure(new_structure)
        assert inp.structure == new_structure

        # Compatible with Pickle and MSONable?
        self.serialize_with_pickle(inp, test_eq=False)
        self.assertMSONable(inp)

        #new_inp = BasicAbinitInput(structure=inp.structure, pseudos=inp.pseudos, abi_kwargs=inp.vars)
        #new_inp.pop_par_vars(new_inp)
        #new_inp["npband"] = 2
        #new_inp["npfft"] = 3
        #popped = new_inp.pop_par_vars(new_inp)
        #assert popped["npband"] == 2 and "npband" not in new_inp
        #assert popped["npfft"] == 3 and "npfft" not in new_inp

    def test_input_errors(self):
        """Testing typical BasicAbinitInput Error"""
        #si_structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

        # Ambiguous list of pseudos.
        #with self.assertRaises(BasicAbinitInput.Error):
        #    BasicAbinitInput(si_structure, pseudos=abidata.pseudos("14si.pspnc", "Si.oncvpsp"))

        # Pseudos do not match structure.
        #with self.assertRaises(BasicAbinitInput.Error):
        #    BasicAbinitInput(si_structure, pseudos=abidata.pseudos("13al.981214.fhi"))

        al_negative_volume = dict(
            ntypat=1,
            natom=1,
            typat=[1],
            znucl=13,
            acell=3*[7.60],
            rprim=[[0.0,  0.5,  0.5],
                   [-0.5,  -0.0,  -0.5],
                   [0.5,  0.5,  0.0]],
            xred=[ [0.0 , 0.0 , 0.0]],
        ),

        # Negative triple product.
        #with self.assertRaises(BasicAbinitInput.Error):
        #    s = abidata.structure_from_ucell("Al-negative-volume")
        #    BasicAbinitInput(s, pseudos=abidata.pseudos("13al.981214.fhi"))

    #def test_helper_functions(self):
    #    """Testing BasicAbinitInput helper functions."""
    #    pseudo = abidata.pseudo("14si.pspnc")
    #    pseudo_dir = os.path.dirname(pseudo.filepath)
    #    inp = BasicAbinitInput(structure=abidata.cif_file("si.cif"), pseudos=pseudo.basename, pseudo_dir=pseudo_dir)

    #    inp.set_kmesh(ngkpt=(1, 2, 3), shiftk=(1, 2, 3, 4, 5, 6))
    #    assert inp["kptopt"] == 1 and inp["nshiftk"] == 2
    #    ngkpt, shiftk = inp.get_ngkpt_shiftk()
    #    assert ngkpt.tolist() == [1, 2, 3]
    #    assert len(shiftk) == 2 and shiftk.ravel().tolist() == [1, 2, 3, 4, 5, 6]

    #    inp.pop("ngkpt")
    #    kptrlatt = [1, 0, 0, 0, 4, 0, 0, 0, 8]
    #    shiftk = (0.5, 0.0, 0.0)
    #    inp.set_vars(kptrlatt=kptrlatt, nshiftk=1, shiftk=shiftk)
    #    ngkpt, shiftk = inp.get_ngkpt_shiftk()
    #    assert ngkpt.tolist() == [1, 4, 8]
    #    assert len(shiftk) == 1 and shiftk.ravel().tolist() == [0.5, 0.0, 0.0]

    #    inp.set_vars(kptrlatt = [1, 2, 0, 0, 4, 0, 0, 0, 8], nshiftk=1, shiftk=shiftk)
    #    ngkpt, shiftk = inp.get_ngkpt_shiftk()
    #    assert ngkpt is None

    #    inp.set_gamma_sampling()
    #    assert inp["kptopt"] == 1 and inp["nshiftk"] == 1
    #    assert np.all(inp["shiftk"] == 0)

    #    inp.set_autokmesh(nksmall=2)
    #    assert inp["kptopt"] == 1 and np.all(inp["ngkpt"] == [2, 2, 2]) and inp["nshiftk"] == 4

    #    inp.set_kpath(ndivsm=3, kptbounds=None)
    #    assert inp["ndivsm"] == 3 and inp["iscf"] == -2 and len(inp["kptbounds"]) == 12

    #    inp.set_qpath(ndivsm=3, qptbounds=None)
    #    assert len(inp["ph_qpath"]) == 12 and inp["ph_nqpath"] == 12 and inp["ph_ndivsm"] == 3

    #    inp.set_phdos_qmesh(nqsmall=16, method="tetra")
    #    assert inp["ph_intmeth"] == 2 and np.all(inp["ph_ngqpt"] == 16) and np.all(inp["ph_qshift"] == 0)

    #    inp.set_kptgw(kptgw=(1, 2, 3, 4, 5, 6), bdgw=(1, 2))
    #    assert inp["nkptgw"] == 2 and np.all(inp["bdgw"].ravel() == np.array(len(inp["kptgw"]) * [1,2]).ravel())

    #    linps = inp.linspace("ecut", 2, 6, num=3, endpoint=True)
    #    assert len(linps) == 3 and (linps[0]["ecut"] == 2 and (linps[-1]["ecut"] == 6))

    #    ranginps = inp.arange("ecut", start=3, stop=5, step=1)
    #    assert len(ranginps) == 2 and (ranginps[0]["ecut"] == 3 and (ranginps[-1]["ecut"] == 4))

    #    with self.assertRaises(inp.Error):
    #        inp.product("ngkpt", "tsmear", [[2, 2, 2], [4, 4, 4]])

    #    prod_inps = inp.product("ngkpt", "tsmear", [[2, 2, 2], [4, 4, 4]], [0.1, 0.2, 0.3])
    #    assert len(prod_inps) == 6
    #    assert prod_inps[0]["ngkpt"] == [2, 2, 2] and prod_inps[0]["tsmear"] == 0.1
    #    assert prod_inps[-1]["ngkpt"] ==  [4, 4, 4] and prod_inps[-1]["tsmear"] == 0.3

    #    inp["kptopt"] = 4

    #def test_dict_methods(self):
    #    """ Testing BasicAbinitInput dict methods """
    #    inp = ebands_input(abidata.cif_file("si.cif"), abidata.pseudos("14si.pspnc"), kppa=10, ecut=2)[0]
    #    inp = ideco.SpinDecorator("spinor")(inp)
    #    inp_dict = inp.as_dict()
    #    #self.assertIsInstance(inp_dict['abi_kwargs'], collections.OrderedDict)
    #    assert "abi_args" in inp_dict and len(inp_dict["abi_args"]) == len(inp)
    #    assert all(k in inp for k, _ in inp_dict["abi_args"])
    #    self.assertMSONable(inp)

    #def test_dfpt_methods(self):
    #    """Testing DFPT methods."""
    #    gs_inp = BasicAbinitInput(structure=abidata.structure_from_ucell("AlAs"),
    #                         pseudos=abidata.pseudos("13al.981214.fhi", "33as.pspnc"))

    #    gs_inp.set_vars(
    #        nband=4,
    #        ecut=2,
    #        ngkpt=[4, 4, 4],
    #        nshiftk=4,
    #        shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
    #                0.0, 0.5, 0.0,
    #                0.5, 0.0, 0.0,
    #                0.5, 0.5, 0.5],
    #        #shiftk=[0, 0, 0],
    #        paral_kgb=1,
    #        nstep=25,
    #        tolvrs=1.0e-10,
    #    )

    #    # qpt is not in gs_inp and not passed to method.
    #    with self.assertRaises(ValueError):
    #        gs_inp.abiget_irred_phperts()

    #    ################
    #    # Phonon methods
    #    ################
    #    with self.assertRaises(gs_inp.Error):
    #        try:
    #            ddk_inputs = gs_inp.make_ddk_inputs(tolerance={"tolfoo": 1e10})
    #        except Exception as exc:
    #            print(exc)
    #            raise

    #    phg_inputs = gs_inp.make_ph_inputs_qpoint(qpt=(0, 0, 0), tolerance=None)
    #    #print("phonon inputs at Gamma\n", phg_inputs)
    #    assert len(phg_inputs) == 2
    #    assert np.all(phg_inputs[0]["rfatpol"] == [1, 1])
    #    assert np.all(phg_inputs[1]["rfatpol"] == [2, 2])
    #    assert all(np.all(inp["rfdir"] == [1, 0, 0] for inp in phg_inputs))
    #    assert all(np.all(inp["kptopt"] == 2 for inp in phg_inputs))

    #    # Validate with Abinit
    #    self.abivalidate_multi(phg_inputs)
