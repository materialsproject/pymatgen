# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import copy
import glob
import json
import os
import unittest

import sys
from pymatgen import Molecule
from pymatgen.io.qchem_deprecated import QcTask, QcInput, QcOutput, QcNucVeloc
from pymatgen.util.testing import PymatgenTest

__author__ = 'xiaohuiqu'


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules")


coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]
mol = Molecule(["C", "H", "H", "H", "Cl"], coords)

coords2 = [[0.0, 0.0, -2.4],
           [0.0, 0.0, 0.0],
           [0.0, 0.0, 2.4]]
heavy_mol = Molecule(["Br", "Cd", "Br"], coords2)

coords3 = [[2.632273, -0.313504, -0.750376],
           [3.268182, -0.937310, -0.431464],
           [2.184198, -0.753305, -1.469059]]

water_mol = Molecule(["O", "H", "H"], coords3)


class QcTaskTest(PymatgenTest):

    def elementary_io_verify(self, text, qctask):
        self.to_and_from_dict_verify(qctask)
        self.from_string_verify(contents=text, ref_dict=qctask.as_dict())

    def to_and_from_dict_verify(self, qctask):
        """
        Helper function. This function should be called in each specific test.
        """
        d1 = qctask.as_dict()
        qc2 = QcTask.from_dict(d1)
        d2 = qc2.as_dict()
        self.assertEqual(d1, d2)

    def from_string_verify(self, contents, ref_dict):
        qctask = QcTask.from_string(contents)
        d2 = qctask.as_dict()
        self.assertEqual(ref_dict, d2)

    def test_read_zmatrix(self):
        contents = '''$moLEcule

 1 2

 S
 C  1 1.726563
 H  2 1.085845 1 119.580615
 C  2 1.423404 1 114.230851 3 -180.000000 0
 H  4 1.084884 2 122.286346 1 -180.000000 0
 C  4 1.381259 2 112.717365 1 0.000000 0
 H  6 1.084731 4 127.143779 2 -180.000000 0
 C  6 1.415867 4 110.076147 2 0.000000 0
 F  8 1.292591 6 124.884374 4 -180.000000 0
$end

$reM
   BASIS  =  6-31+G*
   EXCHANGE  =  B3LYP
   jobtype  =  freq
$end

'''
        qctask = QcTask.from_string(contents)
        ans = '''$molecule
 1  2
 S           0.00000000        0.00000000        0.00000000
 C           0.00000000        0.00000000        1.72656300
 H          -0.94431813        0.00000000        2.26258784
 C           1.29800105       -0.00000002        2.31074808
 H           1.45002821       -0.00000002        3.38492732
 C           2.30733813       -0.00000003        1.36781908
 H           3.37622632       -0.00000005        1.55253338
 C           1.75466906       -0.00000003        0.06427152
 F           2.44231414       -0.00000004       -1.03023099
$end


$rem
   jobtype = freq
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        ans_tokens = ans.split('\n')
        ans_text_part = ans_tokens[:2] + ans_tokens[11:]
        ans_coords_part = ans_tokens[2:11]
        converted_tokens = str(qctask).split('\n')
        converted_text_part = converted_tokens[:2] + converted_tokens[11:]
        converted_coords_part = converted_tokens[2:11]
        self.assertEqual(ans_text_part, converted_text_part)
        for ans_coords, converted_coords in zip(ans_coords_part,
                                                converted_coords_part):
            ans_coords_tokens = ans_coords.split()
            converted_coords_tokens = converted_coords.split()
            self.assertEqual(ans_coords_tokens[0], converted_coords_tokens[0])
            xyz1 = ans_coords_tokens[1:]
            xyz2 = converted_coords_tokens[1:]
            for t1, t2 in zip(xyz1, xyz2):
                self.assertTrue(abs(float(t1)-float(t2)) < 0.0001)

    def test_no_mol(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 -1  2
 read
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        qctask = QcTask(molecule="READ", title="Test Methane",
                        exchange="B3LYP", jobtype="SP", charge=-1,
                        spin_multiplicity=2,
                        basis_set="6-31+G*")
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_simple_basis_str(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_fragmented_molecule(self):
        mol1 = copy.deepcopy(mol)
        mol1.set_charge_and_spin(1, 2)
        mol2 = copy.deepcopy(water_mol)
        mol2.set_charge_and_spin(-1, 2)
        qctask = QcTask([mol1, mol2], title="Test Fragments", exchange="B3LYP",
                        jobtype="bsse", charge=0, spin_multiplicity=3, basis_set="6-31++G**")
        ans = """$comment
 Test Fragments
$end


$molecule
 0  3
--
 1  2
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
--
 -1  2
 O           2.63227300       -0.31350400       -0.75037600
 H           3.26818200       -0.93731000       -0.43146400
 H           2.18419800       -0.75330500       -1.46905900
$end


$rem
   jobtype = bsse
  exchange = b3lyp
     basis = 6-31++g**
$end

"""
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_mixed_basis_str(self):
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set=[("C", "6-311G*"), ("H", "6-31g(d,p)"), ("H", "6-31g(d,p)"),
                                   ("H", "6-31g*"), ("cl", "6-31+g*")])
        ans_mixed = """$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = mixed
$end


$basis
 C    1
 6-311g*
 ****
 H    2
 6-31g(d,p)
 ****
 H    3
 6-31g(d,p)
 ****
 H    4
 6-31g*
 ****
 Cl   5
 6-31+g*
 ****
$end

"""
        self.assertEqual(ans_mixed, str(qctask))
        self.elementary_io_verify(ans_mixed, qctask)
        qctask.set_basis_set("6-31+G*")
        ans_simple = """$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
$end

"""
        self.assertEqual(str(qctask), ans_simple)
        qctask.set_basis_set([("C", "6-311G*"), ("H", "6-31g(d,p)"), ("H", "6-31g(d,p)"),
                              ("H", "6-31g*"), ("cl", "6-31+g*")])
        self.assertEqual(str(qctask), ans_mixed)
        self.elementary_io_verify(ans_mixed, qctask)

    def test_velocities(self):
        qctask = QcTask.from_file(
            os.path.join(test_dir, "qc_aimd",
                         "mg2dig_nvt_langevin.inp"))
        qcnv = QcNucVeloc(
            os.path.join(test_dir, "qc_aimd",
                         "NucVeloc.velocities"))
        velocities = qcnv.velocities[-1]
        qctask.set_velocities(velocities)
        qc_text = str(qctask)
        vel_text = qc_text[qc_text.index("$velocity"):]
        self.assertEqual(vel_text.split("\n")[1].strip(),
                         "8.97607E-05    9.45576E-06   -2.39705E-04")
        self.assertEqual(len(vel_text.split("\n")), 66)
        self.assertEqual(vel_text.split("\n")[-4].strip(),
                         "9.05272E-05    1.11329E-03   -9.17663E-04")
        qctask2 = QcTask.from_string(qc_text)
        self.elementary_io_verify(qc_text, qctask2)

    def test_opt_constraint_str(self):
        opt_coords = [[-1.8438708, 1.7639844, 0.0036111],
                      [-0.3186117, 1.7258535, 0.0241264],
                      [0.1990523, 0.2841796, -0.0277432],
                      [1.7243049, 0.2460376, -0.0067397],
                      [-2.1904881, 2.8181992, 0.0419217],
                      [-2.2554858, 1.2221552, 0.8817436],
                      [-2.2293542, 1.2964646, -0.9274861],
                      [0.0400963, 2.2185950, 0.9541706],
                      [0.0663274, 2.2929337, -0.8514870],
                      [-0.1594453, -0.2084377, -0.9579392],
                      [-0.1860888, -0.2830148, 0.8477023],
                      [2.1362687, 0.7881530, -0.8845274],
                      [2.0709344, -0.8081667, -0.0452220],
                      [2.1094213, 0.7132527, 0.9246668]]
        opt_mol = Molecule(["C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"], opt_coords)
        constraint_dict = {'opt': {'CONSTRAINT': [['tors', 1, 2, 3, 4, 180.0]]}}
        ans = """$molecule
 0  1
 C          -1.84387080        1.76398440        0.00361110
 C          -0.31861170        1.72585350        0.02412640
 C           0.19905230        0.28417960       -0.02774320
 C           1.72430490        0.24603760       -0.00673970
 H          -2.19048810        2.81819920        0.04192170
 H          -2.25548580        1.22215520        0.88174360
 H          -2.22935420        1.29646460       -0.92748610
 H           0.04009630        2.21859500        0.95417060
 H           0.06632740        2.29293370       -0.85148700
 H          -0.15944530       -0.20843770       -0.95793920
 H          -0.18608880       -0.28301480        0.84770230
 H           2.13626870        0.78815300       -0.88452740
 H           2.07093440       -0.80816670       -0.04522200
 H           2.10942130        0.71325270        0.92466680
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
$end


$opt
CONSTRAINT
tors 1 2 3 4 180.0
ENDCONSTRAINT
$end

"""
        qctask = QcTask(opt_mol, exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*",
                        optional_params=constraint_dict)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_opt_fixed_atoms(self):
        fixed_dict = {"opt": {"FIXED": {2: "Y", 3: "XYZ"}}}
        qctask1 = QcTask(mol, exchange="B3LYP", jobtype="SP",
                         basis_set="6-31+G*",
                         optional_params=fixed_dict)
        task_text1 = str(qctask1)
        opt_text1 = task_text1[task_text1.index("$opt"):]
        ans1 = """$opt
FIXED
 2 Y
 3 XYZ
ENDFIXED
$end

"""
        self.assertEqual(opt_text1, ans1)
        self.elementary_io_verify(task_text1, qctask1)

        fixed_n_const_dict = {"opt": {"FIXED": {2: "Y", 3: "XYZ"},
                                      "CONSTRAINT": [['tors', 1, 2, 3, 4, 180.0]]}}
        qctask2 = QcTask(mol, exchange="B3LYP", jobtype="SP",
                        basis_set="6-31+G*",
                        optional_params=fixed_n_const_dict)
        task_text2 = str(qctask2)
        opt_text2 = task_text2[task_text2.index("$opt"):]
        ans2 = """$opt
CONSTRAINT
tors 1 2 3 4 180.0
ENDCONSTRAINT

FIXED
 2 Y
 3 XYZ
ENDFIXED
$end

"""
        self.assertEqual(opt_text2, ans2)
        self.elementary_io_verify(task_text2, qctask2)

    def test_method_keyword(self):
        qctask1 = QcTask(mol, method="B3LYP", jobtype="SP",
                         basis_set="6-31+G*")
        task_text = str(qctask1)
        rem_text = task_text[task_text.index("$rem"):]
        ans = """$rem
  jobtype = sp
   method = b3lyp
    basis = 6-31+g*
$end

"""
        self.assertEqual(rem_text, ans)
        self.elementary_io_verify(task_text, qctask1)

    def test_partial_hessian(self):
        qcinp1 = QcInput.from_file(os.path.join(test_dir, "partial_hessian.qcinp"))
        ans = """$molecule
 0  1
 C          -1.76827000        0.46495000        0.28695000
 O           1.78497000       -0.42034000       -0.39845000
 H          -0.77736000        0.78961000        0.66548000
 H          -1.75896000        0.46604000       -0.82239000
 H          -2.54983000        1.16313000        0.65101000
 H          -1.98693000       -0.55892000        0.65381000
 H           2.14698000       -0.07173000        0.45530000
 H           1.25596000       -1.21510000       -0.13726000
$end


$rem
   jobtype = freq
  exchange = b3lyp
     basis = 6-31g*
     n_sol = 3
     phess = true
$end


$alist
 3
 7
 8
$end

"""
        self.assertEqual(ans, str(qcinp1))
        self.elementary_io_verify(ans, qcinp1.jobs[0])
        qcinp1.jobs[0].params["rem"]["jobtype"] = "sp"
        qcinp1.jobs[0].params["rem"]["phess"] = 3
        qcinp1.jobs[0].set_partial_hessian_atoms([2, 3, 4, 5, 6])
        ans = """$molecule
 0  1
 C          -1.76827000        0.46495000        0.28695000
 O           1.78497000       -0.42034000       -0.39845000
 H          -0.77736000        0.78961000        0.66548000
 H          -1.75896000        0.46604000       -0.82239000
 H          -2.54983000        1.16313000        0.65101000
 H          -1.98693000       -0.55892000        0.65381000
 H           2.14698000       -0.07173000        0.45530000
 H           1.25596000       -1.21510000       -0.13726000
$end


$rem
   jobtype = freq
  exchange = b3lyp
     basis = 6-31g*
     n_sol = 5
     phess = True
$end


$alist
 2
 3
 4
 5
 6
$end

"""
        self.assertEqual(ans, str(qcinp1))

    def test_basis2_mixed(self):
        qcinp1 = QcInput.from_file(os.path.join(test_dir, "basis2_mixed.inp"))
        ans = """$molecule
 0  1
 C          -1.76827000        0.46495000        0.28695000
 O           1.78497000       -0.42034000       -0.39845000
 H          -0.77736000        0.78961000        0.66548000
 H          -1.75896000        0.46604000       -0.82239000
 H          -2.54983000        1.16313000        0.65101000
 H          -1.98693000       -0.55892000        0.65381000
 H           2.14698000       -0.07173000        0.45530000
 H           1.25596000       -1.21510000       -0.13726000
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = mixed
    basis2 = basis2_mixed
  purecart = 1111
$end


$basis
 C    1
 6-311+g(3df)
 ****
 O    2
 aug-cc-pvtz
 ****
 H    3
 6-31g*
 ****
 H    4
 6-31g*
 ****
 H    5
 6-31g*
 ****
 H    6
 6-31g*
 ****
 H    7
 cc-pvdz
 ****
 H    8
 cc-pvdz
 ****
$end


$basis2
 C    1
 sto-3g
 ****
 O    2
 sto-3g
 ****
 H    3
 sto-3g
 ****
 H    4
 sto-3g
 ****
 H    5
 sto-3g
 ****
 H    6
 sto-3g
 ****
 H    7
 sto-3g
 ****
 H    8
 sto-3g
 ****
$end

"""
        self.assertEqual(str(qcinp1), ans)
        self.elementary_io_verify(ans, qcinp1.jobs[0])
        basis2 = qcinp1.jobs[0].params["basis2"]
        qcinp2 = copy.deepcopy(qcinp1)
        qcinp2.jobs[0].set_basis2("3-21g")
        self.assertEqual(qcinp2.jobs[0].params["rem"]["basis2"], "3-21g")
        self.assertFalse("basis2" in qcinp2.jobs[0].params)
        qcinp2.jobs[0].set_basis2(basis2)
        self.assertEqual(str(qcinp2), ans)

    def test_aux_basis_str(self):
        ans_gen = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
    jobtype = freq
   exchange = xygjos
      basis = gen
  aux_basis = gen
$end


$aux_basis
 C
 rimp2-cc-pvdz
 ****
 Cl
 rimp2-aug-cc-pvdz
 ****
 H
 rimp2-cc-pvdz
 ****
$end


$basis
 C
 6-31g*
 ****
 Cl
 6-31+g*
 ****
 H
 6-31g*
 ****
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="xygjos",
                        jobtype="Freq",
                        basis_set={"C": "6-31G*", "h": "6-31g*",
                                   "CL": "6-31+g*"},
                        aux_basis_set={"c": "rimp2-cc-pvdz",
                                       "H": "rimp2-cc-pvdz",
                                       "Cl": "rimp2-aug-cc-pvdz"})
        self.assertEqual(str(qctask), ans_gen)
        self.elementary_io_verify(ans_gen, qctask)
        qctask.set_auxiliary_basis_set([("C", "aug-cc-pvdz"), ("H", "cc-pvdz"), ("H", "cc-pvdz"),
                                        ("H", "cc-pvdz"), ("cl", "rimp2-aug-cc-pvdz")])
        ans_mixed_aux = """$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
    jobtype = freq
   exchange = xygjos
      basis = gen
  aux_basis = mixed
$end


$aux_basis
 C    1
 aug-cc-pvdz
 ****
 H    2
 cc-pvdz
 ****
 H    3
 cc-pvdz
 ****
 H    4
 cc-pvdz
 ****
 Cl   5
 rimp2-aug-cc-pvdz
 ****
$end


$basis
 C
 6-31g*
 ****
 Cl
 6-31+g*
 ****
 H
 6-31g*
 ****
$end

"""
        self.assertEqual(ans_mixed_aux, str(qctask))
        self.elementary_io_verify(ans_mixed_aux, qctask)
        qctask.set_basis_set("6-31+G*")
        qctask.set_auxiliary_basis_set("rimp2-cc-pvdz")
        ans_simple = """$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
    jobtype = freq
   exchange = xygjos
      basis = 6-31+g*
  aux_basis = rimp2-cc-pvdz
$end

"""
        self.assertEqual(ans_simple, str(qctask))
        self.elementary_io_verify(ans_simple, qctask)
        qctask.set_basis_set({"C": "6-31G*", "h": "6-31g*",
                              "CL": "6-31+g*"})
        qctask.set_auxiliary_basis_set([("C", "aug-cc-pvdz"), ("H", "cc-pvdz"), ("H", "cc-pvdz"),
                                        ("H", "cc-pvdz"), ("cl", "rimp2-aug-cc-pvdz")])
        self.assertEqual(ans_mixed_aux, str(qctask))
        self.elementary_io_verify(ans_mixed_aux, qctask)

    def test_ecp_str(self):
        ans = '''$comment
 Test ECP
$end


$molecule
 0  1
 Br          0.00000000        0.00000000       -2.40000000
 Cd          0.00000000        0.00000000        0.00000000
 Br          0.00000000        0.00000000        2.40000000
$end


$rem
   jobtype = opt
  exchange = b3lyp
     basis = gen
       ecp = gen
$end


$basis
 Br
 srlc
 ****
 Cd
 srsc
 ****
$end


$ecp
 Br
 srlc
 ****
 Cd
 srsc
 ****
$end

'''
        qctask = QcTask(heavy_mol, title="Test ECP", exchange="B3LYP",
                        jobtype="Opt",
                        basis_set={"Br": "srlc", "Cd": "srsc"},
                        ecp={"Br": "SrlC", "Cd": "srsc"})
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_set_memory(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
     jobtype = sp
    exchange = b3lyp
       basis = 6-31+g*
  mem_static = 500
   mem_total = 18000
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.set_memory(total=18000, static=500)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_qc42_pcm_solvent_format(self):
        text = '''$molecule
 -1  2
 N          -0.00017869        0.00010707        0.20449990
 H           0.89201838        0.20268122       -0.29656572
 H          -0.62191133        0.67135171       -0.29649162
 H          -0.26987729       -0.87406458       -0.29659779
$end


$rem
         jobtype = sp
        exchange = b3lyp
           basis = 6-31+g*
  solvent_method = pcm
$end


$pcm
       theory   ssvpe
     vdwscale   1.1
$end


$pcm_solvent
  dielectric   78.3553
$end

'''
        qctask_qc41 = QcTask.from_string(text)
        qctask_qc42 = copy.deepcopy(qctask_qc41)
        solvent_params = qctask_qc42.params.pop("pcm_solvent")
        qctask_qc42.params["solvent"] = solvent_params
        ans = '''$molecule
 -1  2
 N          -0.00017869        0.00010707        0.20449990
 H           0.89201838        0.20268122       -0.29656572
 H          -0.62191133        0.67135171       -0.29649162
 H          -0.26987729       -0.87406458       -0.29659779
$end


$rem
         jobtype = sp
        exchange = b3lyp
           basis = 6-31+g*
  solvent_method = pcm
$end


$pcm
    theory   ssvpe
  vdwscale   1.1
$end


$solvent
  dielectric   78.3553
$end

'''
        self.assertEqual(str(qctask_qc42), ans)
        self.elementary_io_verify(ans, qctask_qc42)

    def test_set_max_num_of_scratch_files(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
           jobtype = sp
          exchange = b3lyp
             basis = 6-31+g*
  max_sub_file_num = 500
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.set_max_num_of_scratch_files(500)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_set_max_scf_iterations(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
         jobtype = sp
        exchange = b3lyp
           basis = 6-31+g*
  max_scf_cycles = 100
   scf_algorithm = diis_gdm
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.set_scf_algorithm_and_iterations(algorithm="diis_gdm",
                                                iterations=100)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_set_scf_convergence_threshold(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
          jobtype = sp
         exchange = b3lyp
            basis = 6-31+g*
  scf_convergence = 8
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.set_scf_convergence_threshold(exponent=8)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_set_integral_threshold(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
    thresh = 14
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.set_integral_threshold(thresh=14)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_set_dft_grid(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
   xc_grid = 000110000590
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.set_dft_grid(radical_points=110, angular_points=590)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_set_scf_initial_guess(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
    jobtype = sp
   exchange = b3lyp
      basis = 6-31+g*
  scf_guess = gwh
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.set_scf_initial_guess("GWH")
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_geom_opt_max_cycles(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 1  2
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
              jobtype = sp
             exchange = b3lyp
                basis = 6-31+g*
  geom_opt_max_cycles = 100
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP", charge=1, spin_multiplicity=2,
                        basis_set="6-31+G*")
        qctask.set_geom_max_iterations(100)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_set_geom_opt_coords_type(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
          jobtype = sp
         exchange = b3lyp
            basis = 6-31+g*
  geom_opt_coords = 0
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.set_geom_opt_coords_type("cartesian")
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_scale_geom_opt_threshold(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
                    jobtype = sp
                   exchange = b3lyp
                      basis = 6-31+g*
  geom_opt_tol_displacement = 120
        geom_opt_tol_energy = 10
      geom_opt_tol_gradient = 30
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.scale_geom_opt_threshold(gradient=0.1, displacement=0.1,
                                        energy=0.1)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_set_geom_opt_use_gdiis(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
            jobtype = sp
           exchange = b3lyp
              basis = 6-31+g*
  geom_opt_max_diis = -1
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.set_geom_opt_use_gdiis()
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_disable_symmetry(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
     jobtype = sp
    exchange = b3lyp
       basis = 6-31+g*
  sym_ignore = True
    symmetry = False
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.disable_symmetry()
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_use_cosmo(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
             jobtype = sp
            exchange = b3lyp
               basis = 6-31+g*
  solvent_dielectric = 35.0
      solvent_method = cosmo
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.use_cosmo(dielectric_constant=35.0)
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_wrap_comment(self):
        ans = '''$comment
 5_2_2_methoxyethoxy_ethoxy_6_nitro_1_3_dihydro_2_1_3_benzothiadiazole
singlet neutral B3lYP/6-31+G* geometry optimization
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        qctask = QcTask(mol, title="    5_2_2_methoxyethoxy_ethoxy_6_nitro_1_3_dihydro_2_1_3_benzothiadiazole singlet "
                                   "neutral B3lYP/6-31+G* geometry optimization", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)
        title = ''' MgBPh42 singlet neutral PBE-D3/6-31+G* geometry optimization
<SCF Fix Strategy>{
    "current_method_id": 1,
    "methods": [
        "increase_iter",
        "diis_gdm",
        "gwh",
        "rca",
        "gdm",
        "core+gdm"
    ]
}</SCF Fix Strategy>'''
        ans = '''$comment
 MgBPh42 singlet neutral PBE-D3/6-31+G* geometry optimization
<SCF Fix Strategy>{
    "current_method_id": 1,
    "methods": [
        "increase_iter",
        "diis_gdm",
        "gwh",
        "rca",
        "gdm",
        "core+gdm"
    ]
}</SCF Fix Strategy>
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        qctask = QcTask(mol, title=title, exchange="B3LYP", jobtype="SP", basis_set="6-31+G*")
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)
        title = " 5_2_2_methoxyethoxy_ethoxy_6_nitro_1_3_dihydro_2_1_3_benzothiadiazole singlet neutral " \
                "B3lYP/6-31+G* geometry optimization" + \
                '''<SCF Fix Strategy>{
    "current_method_id": 1,
    "methods": [
        "increase_iter",
        "diis_gdm",
        "gwh",
        "rca",
        "gdm",
        "core+gdm"
    ]
}</SCF Fix Strategy>'''
        qctask = QcTask(mol, title=title, exchange="B3LYP", jobtype="SP", basis_set="6-31+G*")
        self.elementary_io_verify(str(qctask), qctask)

    def test_use_pcm_qc41(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
         jobtype = sp
        exchange = b3lyp
           basis = 6-31+g*
  solvent_method = pcm
$end


$pcm
     radii   uff
    theory   ssvpe
  vdwscale   1.1
$end


$pcm_solvent
  dielectric   78.3553
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.use_pcm(solvent_key="pcm_solvent")
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.use_pcm(pcm_params={"Radii": "FF",
                                   "Theory": "CPCM",
                                   "SASrad": 1.5,
                                   "HPoints": 1202},
                       solvent_params={"Dielectric": 20.0,
                                       "Temperature": 300.75,
                                       "NSolventAtoms": 2,
                                       "SolventAtom": [[8, 1, 186, 1.30],
                                                       [1, 2, 187, 1.01]]},
                       radii_force_field="OPLSAA",
                       solvent_key="pcm_solvent")
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
         jobtype = sp
        exchange = b3lyp
           basis = 6-31+g*
      force_fied = oplsaa
  solvent_method = pcm
$end


$pcm
   hpoints   1202
     radii   bondi
    sasrad   1.5
    theory   cpcm
  vdwscale   1.1
$end


$pcm_solvent
     dielectric   20.0
  nsolventatoms   2
    solventatom   8    1    186  1.30
    solventatom   1    2    187  1.01
    temperature   300.75
$end

'''
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_use_pcm_qc42(self):
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
         jobtype = sp
        exchange = b3lyp
           basis = 6-31+g*
  solvent_method = pcm
$end


$pcm
     radii   uff
    theory   ssvpe
  vdwscale   1.1
$end


$solvent
  dielectric   78.3553
$end

'''
        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.use_pcm()
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

        qctask = QcTask(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qctask.use_pcm(pcm_params={"Radii": "FF",
                                   "Theory": "CPCM",
                                   "SASrad": 1.5,
                                   "HPoints": 1202},
                       solvent_params={"Dielectric": 20.0,
                                       "Temperature": 300.75,
                                       "NSolventAtoms": 2,
                                       "SolventAtom": [[8, 1, 186, 1.30],
                                                       [1, 2, 187, 1.01]]},
                       radii_force_field="OPLSAA")
        ans = '''$comment
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
         jobtype = sp
        exchange = b3lyp
           basis = 6-31+g*
      force_fied = oplsaa
  solvent_method = pcm
$end


$pcm
   hpoints   1202
     radii   bondi
    sasrad   1.5
    theory   cpcm
  vdwscale   1.1
$end


$solvent
     dielectric   20.0
  nsolventatoms   2
    solventatom   8    1    186  1.30
    solventatom   1    2    187  1.01
    temperature   300.75
$end

'''
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)

    def test_ghost_atoms(self):
        qctask = QcTask(mol, charge=0, spin_multiplicity=1, exchange="B3LYP", ghost_atoms=[2, 4])
        ans = """$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 @H          1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 @Cl        -0.51336000        0.88916500       -0.36300000
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
$end

"""
        self.assertEqual(str(qctask), ans)
        self.elementary_io_verify(ans, qctask)
        mol1 = copy.deepcopy(mol)
        mol1.set_charge_and_spin(1, 2)
        mol2 = copy.deepcopy(water_mol)
        mol2.set_charge_and_spin(-1, 2)
        qctask = QcTask([mol1, mol2], title="Test Fragments", exchange="B3LYP",
                        jobtype="bsse", charge=0, spin_multiplicity=3, basis_set="6-31++G**",
                        ghost_atoms=[1, 2, 3, 5])
        self.elementary_io_verify(str(qctask), qctask)
        qctask = QcTask(mol, charge=0, spin_multiplicity=2, exchange="B3LYP", ghost_atoms=[2])
        self.assertEqual(qctask.spin_multiplicity, 2)


class TestQcInput(PymatgenTest):
    def test_str_and_from_string(self):
        ans = '''$comment
 Test Methane Opt
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
   jobtype = opt
  exchange = b3lyp
     basis = 6-31+g*
$end


@@@


$comment
 Test Methane Frequency
$end


$molecule
 read
$end


$rem
   jobtype = freq
  exchange = b3lyp
     basis = 6-31+g*
$end


@@@


$comment
 Test Methane Single Point Energy
$end


$molecule
 read
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-311+g(3df,2p)
$end

'''
        qctask1 = QcTask(mol, title="Test Methane Opt", exchange="B3LYP",
                         jobtype="Opt", basis_set="6-31+G*")
        qctask2 = QcTask(molecule="read", title="Test Methane Frequency",
                         exchange="B3LYP", jobtype="Freq", basis_set="6-31+G*")
        qctask3 = QcTask(title="Test Methane Single Point Energy",
                         exchange="B3LYP", jobtype="SP",
                         basis_set="6-311+G(3df,2p)")
        qcinp1 = QcInput(jobs=[qctask1, qctask2, qctask3])
        self.assertEqual(str(qcinp1), ans)
        qcinp2 = QcInput.from_string(ans)
        self.assertEqual(qcinp1.as_dict(), qcinp2.as_dict())

        qcinp_mgbf4 = QcInput.from_file(os.path.join(test_dir, "MgBF4_b_overalpped.qcinp"))
        self.assertEqual(qcinp_mgbf4.jobs[0].ghost_atoms, [0])

    def test_to_and_from_dict(self):
        qctask1 = QcTask(mol, title="Test Methane Opt", exchange="B3LYP",
                         jobtype="Opt", basis_set="6-31+G*")
        qctask2 = QcTask(molecule="read", title="Test Methane Frequency",
                         exchange="B3LYP", jobtype="Freq",
                         basis_set="6-31+G*")
        qctask3 = QcTask(title="Test Methane Single Point Energy",
                         exchange="B3LYP", jobtype="SP",
                         basis_set="6-311+G(3df,2p)")
        qcinp1 = QcInput(jobs=[qctask1, qctask2, qctask3])
        d1 = qcinp1.as_dict()
        qcinp2 = QcInput.from_dict(d1)
        d2 = qcinp2.as_dict()
        self.assertEqual(d1, d2)


class TestQcOutput(PymatgenTest):
    def test_energy(self):
        ref_energies_text = '''
{
    "hf-rimp2.qcout": {
        "RIMP2": -2726.6860779805256,
        "SCF": -2721.541435904716
    },
    "hf_b3lyp.qcout": {
        "SCF": -2733.1747178920828
    },
    "hf_ccsd(t).qcout": {
        "CCSD": -2726.7627121001865,
        "CCSD(T)": -2726.8283514003333,
        "MP2": -2726.685664155242,
        "SCF": -2721.5414360843106
    },
    "hf_cosmo.qcout": {
        "SCF": -2721.1752937496067
    },
    "hf_hf.qcout": {
        "SCF": -2721.541435904716
    },
    "hf_lxygjos.qcout": {
        "SCF": -2724.0769973875713,
        "XYGJ-OS": -2726.3445157759393
    },
    "hf_mosmp2.qcout": {
        "MOS-MP2": -2725.302538779482,
        "SCF": -2721.541435904716
    },
    "hf_mp2.qcout": {
        "MP2": -2726.685661962005,
        "SCF": -2721.541435904716
    },
    "hf_pcm.qcout": {
        "SCF": -2720.703940318968
    },
    "hf_qcisd(t).qcout": {
        "QCISD": -2726.7853751012344,
        "QCISD(T)": -2726.8346541282745,
        "SCF": -2721.5414360843106
    },
    "hf_riccsd(t).qcout": {
        "CCSD": -2726.7641790658904,
        "CCSD(T)": -2726.829853468723,
        "MP2": -2726.6860802173014,
        "SCF": -2721.5414360843106
    },
    "hf_tpssh.qcout": {
        "SCF": -2732.938974944255
    },
    "hf_xyg3.qcout": {
        "SCF": -2728.769906036435,
        "XYG3": -2731.0640917605806
    },
    "hf_xygjos.qcout": {
        "SCF": -2724.0769973875713,
        "XYGJ-OS": -2726.3447230967517
    },
    "hf_wb97xd_gen_scfman.qcout": {
        "GEN_SCFMAN": -30051.134375112342,
        "GEN_SCFMAN": -30051.296918174274,
        "GEN_SCFMAN": -30051.395763612905,
        "GEN_SCFMAN": -30051.458839496852,
        "GEN_SCFMAN": -30051.487970700582,
        "GEN_SCFMAN": -30051.490764186092,
        "GEN_SCFMAN": -30051.491278372443,
        "GEN_SCFMAN": -30051.491359704556,
        "GEN_SCFMAN": -30051.491369799976
    } 
}'''
        ref_energies = json.loads(ref_energies_text)
        parsed_energies = dict()
        # noinspection PyUnresolvedReferences
        for filename in glob.glob(os.path.join(test_dir, "qchem_energies",
                                                         "*.qcout")):
            molname = os.path.basename(filename)
            qcout = QcOutput(filename)
            d = dict(qcout.data[0]["energies"])
            parsed_energies[molname] = d
        self.assertEqual(sorted(ref_energies.keys()),
                         sorted(parsed_energies.keys()))
        mols = sorted(ref_energies.keys())
        for molname in mols:
            self.assertEqual(sorted(ref_energies[molname].keys()),
                             sorted(parsed_energies[molname].keys()))
            methods = sorted(ref_energies[molname].keys())
            for method in methods:
                self.assertAlmostEqual(ref_energies[molname][method],
                                       parsed_energies[molname][method], 2)

    def test_unable_to_determine_lambda_in_geom_opt(self):
        filename = os.path.join(test_dir, "unable_to_determine_lambda_in_geom_opt.qcout")
        qcout = QcOutput(filename)
        self.assertTrue(qcout.data[0]['has_error'])
        self.assertEqual(qcout.data[0]['errors'],
                         ['Lamda Determination Failed',
                          'Geometry optimization failed'])

    def test_geom_opt(self):
        filename = os.path.join(test_dir, "thiophene_wfs_5_carboxyl.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(qcout.data[0]["jobtype"], "opt")
        ans_energies = [(u'SCF', -20179.886441383995),
                        (u'SCF', -20180.12187218424),
                        (u'SCF', -20180.150524404988),
                        (u'SCF', -20180.151628362753),
                        (u'SCF', -20180.151810235497),
                        (u'SCF', -20180.15180854295)]
        self.assertEqual(qcout.data[0]["energies"], ans_energies)
        ans_mol1 = '''Full Formula (H4 C5 S1 O2)
Reduced Formula: H4C5SO2
Charge = -1, Spin Mult = 2
Sites (12)
0 C     0.158839    -0.165379     0.000059
1 C    -0.520531    -1.366720     0.000349
2 C    -1.930811    -1.198460    -0.000041
3 C    -2.297971     0.127429    -0.000691
4 S    -0.938312     1.189630     0.000400
5 H    -0.014720    -2.325340     0.000549
6 H    -2.641720    -2.017721    -0.000161
7 H    -3.301032     0.535659    -0.001261
8 C     1.603079     0.076231    -0.000101
9 O     2.131988     1.173581    -0.000330
10 O     2.322109    -1.079218    -0.000021
11 H     3.262059    -0.820188    -0.000171'''
        ans_mol_last = '''Full Formula (H4 C5 S1 O2)
Reduced Formula: H4C5SO2
Charge = -1, Spin Mult = 2
Sites (12)
0 C     0.194695    -0.158362    -0.001887
1 C    -0.535373    -1.381241    -0.001073
2 C    -1.927071    -1.199274    -0.000052
3 C    -2.332651     0.131916     0.000329
4 S    -0.942111     1.224916    -0.001267
5 H    -0.038260    -2.345185    -0.001256
6 H    -2.636299    -2.025939     0.000620
7 H    -3.339756     0.529895     0.001288
8 C     1.579982     0.071245    -0.002733
9 O     2.196383     1.165675    -0.000178
10 O     2.352341    -1.114671     0.001634
11 H     3.261096    -0.769470     0.003158'''
        self.assertEqual(qcout.data[0]["molecules"][0].__str__(), ans_mol1)
        self.assertEqual(str(qcout.data[0]["molecules"][-1]), ans_mol_last)
        self.assertFalse(qcout.data[0]["has_error"])
        ans_gradient = [{'max_gradient': 0.07996,
                         'gradients': [(-0.0623076, -0.0157774, -2.05e-05),
                                       (0.0260287, 0.0289157, -6e-06),
                                       (-0.015738, 0.0103583, 1.87e-05),
                                       (0.0260219, -0.0028, -1.36e-05),
                                       (-0.0043158, -0.0245896, 2.83e-05),
                                       (4.8e-05, 0.000782, 1.3e-06),
                                       (0.0014679, 0.0020277, 3.9e-06),
                                       (0.0010437, -1.29e-05, -1.04e-05),
                                       (0.0799585, 0.0204159, 1e-06),
                                       (-0.0320357, -0.0421461, 2.1e-06),
                                       (-0.0237691, 0.0247526, -4.6e-06),
                                       (0.0035975, -0.0019264, -3e-07)],
                         'rms_gradient': 0.02244},
                        {'max_gradient': 0.02721,
                         'gradients': [(-0.0195677, -0.0008468, -3.2e-06),
                                       (0.0106798, 0.0039494, 1.11e-05),
                                       (-0.0086473, -0.0012624, -8.1e-06),
                                       (0.0065018, 0.0033749, 5e-07),
                                       (0.0002581, -0.0060831, 7.2e-06),
                                       (-0.0004373, -0.000504, 1.4e-06),
                                       (0.0003216, 0.0001059, -9e-07),
                                       (-0.000814, -5.03e-05, 3e-07),
                                       (0.0272109, 0.001408, -2.06e-05),
                                       (-0.0086971, -0.009251, 8.3e-06),
                                       (-0.0080925, 0.0112191, 2.9e-06),
                                       (0.0012838, -0.0020597, 1.1e-06)],
                         'rms_gradient': 0.007037},
                        {'max_gradient': 0.003444,
                         'gradients': [(0.0021606, 0.0013094, -1.68e-05),
                                       (0.0005757, -0.0002616, -1e-05),
                                       (2.73e-05, -0.0002868, 1.5e-05),
                                       (0.0001088, 0.0006944, -1.23e-05),
                                       (0.0006912, -0.0006523, 6.1e-06),
                                       (-0.0004191, -9.32e-05, -1.3e-06),
                                       (0.0002288, 3.98e-05, 1.8e-06),
                                       (-8.99e-05, -0.0002338, -3.2e-06),
                                       (1.95e-05, -0.0034439, 7.08e-05),
                                       (-0.0008228, -9.18e-05, -2.77e-05),
                                       (-0.0018054, 0.0034031, -2.21e-05),
                                       (-0.0006747, -0.0003834, -3e-07)],
                         'rms_gradient': 0.001008},
                        {'max_gradient': 0.002367,
                         'gradients': [(-0.0001646, 0.0006149, 4.17e-05),
                                       (-0.0004516, -0.0003116, 1.28e-05),
                                       (0.0003366, -3.27e-05, -1.59e-05),
                                       (-0.0003164, 0.0001775, 1.37e-05),
                                       (0.0001399, -0.0001201, -6.9e-06),
                                       (-0.0001374, -1.58e-05, 9e-07),
                                       (-1.19e-05, -3.93e-05, -3.3e-06),
                                       (-1.76e-05, -0.0001233, 5.1e-06),
                                       (9.73e-05, -0.0023668, -0.0001609),
                                       (0.0006998, 0.0009023, 6.31e-05),
                                       (-0.0002169, 0.0014874, 4.95e-05),
                                       (4.28e-05, -0.0001724, 2e-07)],
                         'rms_gradient': 0.0005339},
                        {'max_gradient': 0.001246,
                         'gradients': [(-6.88e-05, 0.0001757, -8.32e-05),
                                       (-0.0002264, -0.0001306, -1.93e-05),
                                       (0.0001526, -1.39e-05, 2.05e-05),
                                       (-0.0001401, 3.8e-06, -2.05e-05),
                                       (1.52e-05, 0.0001152, 8e-06),
                                       (2.01e-05, -3.69e-05, -1e-06),
                                       (-3.62e-05, -3.51e-05, 5.5e-06),
                                       (1.01e-05, -1.23e-05, -6.8e-06),
                                       (9.73e-05, -0.0012462, 0.0003246),
                                       (0.0003926, 0.0008331, -0.0001269),
                                       (-0.0002294, 0.000281, -0.0001009),
                                       (1.3e-05, 6.61e-05, 0.0)],
                         'rms_gradient': 0.0002814},
                        {'max_gradient': 0.0006359,
                         'gradients': [(0.0001036, -0.0001339, 0.0001633),
                                       (0.0001003, 6.98e-05, 3.43e-05),
                                       (-8.28e-05, 1.1e-05, -3.31e-05),
                                       (6.2e-05, -0.0001068, 3.41e-05),
                                       (-5.02e-05, 0.0001346, -1.18e-05),
                                       (8.72e-05, -7.3e-06, 1.5e-06),
                                       (-1.7e-05, 4.9e-06, -1.05e-05),
                                       (1.29e-05, 5.9e-05, 1.26e-05),
                                       (-0.0001059, -5.4e-06, -0.0006359),
                                       (-1.48e-05, 0.0002152, 0.0002469),
                                       (-0.0001335, -0.0003534, 0.0001988),
                                       (3.83e-05, 0.0001124, -1e-07)],
                         'rms_gradient': 0.0001535}]
        self.assertEqual(qcout.data[0]["gradients"], ans_gradient)
        ans_inp = '''$molecule
 -1  2
 C           0.15884000       -0.16538000        0.00006000
 C          -0.52053000       -1.36672000        0.00035000
 C          -1.93081000       -1.19846000       -0.00004000
 C          -2.29797000        0.12743000       -0.00069000
 S          -0.93831000        1.18963000        0.00040000
 H          -0.01472000       -2.32534000        0.00055000
 H          -2.64172000       -2.01772000       -0.00016000
 H          -3.30103000        0.53566000       -0.00126000
 C           1.60308000        0.07623000       -0.00010000
 O           2.13199000        1.17358000       -0.00033000
 O           2.32211000       -1.07922000       -0.00002000
 H           3.26206000       -0.82019000       -0.00017000
$end


$rem
   jobtype = opt
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        self.assertEqual(str(qcout.data[0]['input']), ans_inp)
        self.assertTrue(qcout.data[0]['gracefully_terminated'])
        ans_scf_iter = [[(-743.3130310589, 0.0561),
                         (-741.3557302205, 0.00841),
                         (-740.7031048846, 0.0157),
                         (-741.5589873953, 0.00303),
                         (-741.5918010434, 0.00118),
                         (-741.5966923809, 0.000332),
                         (-741.5970287119, 0.000158),
                         (-741.5971282029, 4.38e-05),
                         (-741.5971448077, 2.17e-05),
                         (-741.5971501973, 7.7e-06),
                         (-741.5971533576, 5.05e-06),
                         (-741.5971541122, 2.7e-06),
                         (-741.5971544119, 9.48e-07),
                         (-741.5971544408, 2.61e-07),
                         (-741.5971544436, 1.21e-07),
                         (-741.5971544441, 5.45e-08),
                         (-741.5971544442, 1.77e-08),
                         (-741.5971544442, 7.79e-09)],
                        [(-741.5552794274, 0.00265),
                         (-741.6048574279, 0.000515),
                         (-741.6037290502, 0.000807),
                         (-741.6056978336, 0.000188),
                         (-741.6057976553, 4.78e-05),
                         (-741.6058045572, 1.54e-05),
                         (-741.6058057373, 4.51e-06),
                         (-741.6058061671, 2.91e-06),
                         (-741.6058062822, 8.32e-07),
                         (-741.6058063435, 7.17e-07),
                         (-741.6058063636, 1.97e-07),
                         (-741.6058063662, 5.03e-08),
                         (-741.6058063666, 3.35e-08),
                         (-741.6058063666, 1.24e-08),
                         (-741.6058063666, 5.25e-09)],
                        [(-741.6023833754, 0.0013),
                         (-741.6065067966, 0.000305),
                         (-741.6057886337, 0.000559),
                         (-741.6068434004, 7.61e-05),
                         (-741.6068555361, 3.4e-05),
                         (-741.6068589376, 5.66e-06),
                         (-741.6068591778, 2.95e-06),
                         (-741.60685927, 1.27e-06),
                         (-741.6068592962, 4.82e-07),
                         (-741.6068593106, 3.84e-07),
                         (-741.6068593157, 9.23e-08),
                         (-741.6068593162, 2.49e-08),
                         (-741.6068593163, 1.52e-08),
                         (-741.6068593163, 5.71e-09)],
                        [(-741.6012175391, 0.000209),
                         (-741.6068794773, 7.2e-05),
                         (-741.606851035, 0.000117),
                         (-741.606899078, 1.53e-05),
                         (-741.6068997567, 6.01e-06),
                         (-741.6068998747, 1.68e-06),
                         (-741.6068998849, 5.32e-07),
                         (-741.6068998857, 2.76e-07),
                         (-741.606899886, 6.41e-08),
                         (-741.606899886, 3.08e-08),
                         (-741.606899886, 9.5e-09)],
                        [(-741.6067290885, 0.0001),
                         (-741.6069044268, 2.64e-05),
                         (-741.6068991026, 5.29e-05),
                         (-741.6069065234, 3.51e-06),
                         (-741.6069065452, 2.49e-06),
                         (-741.6069065686, 3.57e-07),
                         (-741.6069065693, 2.59e-07),
                         (-741.6069065696, 7.05e-08),
                         (-741.6069065696, 4.44e-08),
                         (-741.6069065697, 1.52e-08),
                         (-741.6069065697, 8.17e-09)],
                        [(-741.6074251344, 0.000129),
                         (-741.6069044127, 2.43e-05),
                         (-741.6068998551, 4.95e-05),
                         (-741.6069064294, 4.49e-06),
                         (-741.606906478, 2.77e-06),
                         (-741.6069065049, 5.85e-07),
                         (-741.6069065068, 2.74e-07),
                         (-741.6069065073, 6.99e-08),
                         (-741.6069065074, 3.37e-08),
                         (-741.6069065075, 1.89e-08),
                         (-741.6069065075, 7.38e-09)]]
        self.assertEqual(qcout.data[0]['scf_iteration_energies'], ans_scf_iter)

    def test_multiple_step_job(self):
        filename = os.path.join(test_dir, "CdBr2.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(len(qcout.data), 3)
        self.assertEqual(qcout.data[0]['jobtype'], 'opt')
        self.assertEqual(qcout.data[1]['jobtype'], 'freq')
        ans_thermo_corr_text = '''
{
    "Rotational Enthalpy": 0.025714259,
    "Rotational Entropy": 0.000833523586,
    "Total Enthalpy": 0.199729978,
    "Total Entropy": 0.003218965579,
    "Translational Enthalpy": 0.038549707,
    "Translational Entropy": 0.001851513374,
    "Vibrational Enthalpy": 0.109795116,
    "Vibrational Entropy": 0.000533928619,
    "ZPE": 0.039330241,
    "Zero point vibrational energy": 0.039330241,
    "gas constant (RT)": 0.025714259
}'''
        ans_thermo_corr = json.loads(ans_thermo_corr_text)
        self.assertEqual(sorted(qcout.data[1]['corrections'].keys()),
                         sorted(ans_thermo_corr.keys()))
        for k, ref in ans_thermo_corr.items():
            self.assertAlmostEqual(qcout.data[1]['corrections'][k], ref)
        self.assertEqual(len(qcout.data[1]['molecules']), 1)
        ans_mol1 = '''Full Formula (Cd1 Br2)
Reduced Formula: CdBr2
Charge = 0, Spin Mult = 1
Sites (3)
0 Br     0.000000     0.000000    -2.453720
1 Cd     0.000000     0.000000     0.000000
2 Br     0.000000     0.000000     2.453720'''
        self.assertEqual(str(qcout.data[1]['molecules'][0]), ans_mol1)
        self.assertFalse(qcout.data[1]['has_error'])
        self.assertEqual(qcout.data[1]['gradients'], [])
        ans_inp = '''$molecule
 read
$end


$rem
         jobtype = freq
        exchange = b3lyp
           basis = gen
             ecp = gen
  max_scf_cycles = 100
       scf_guess = gwh
$end


$basis
 Br
 srlc
 ****
 Cd
 srsc
 ****
$end


$ecp
 Br
 srlc
 ****
 Cd
 srsc
 ****
$end

'''
        self.assertEqual(str(qcout.data[1]['input']), ans_inp)
        ans_freq = [{'vib_mode': ((0.17, -0.475, 0.0),
                                  (-0.236, 0.659, 0.0),
                                  (0.17, -0.475, 0.0)),
                     'frequency': 61.36},
                    {'vib_mode': ((-0.475, -0.17, 0.0),
                                  (0.659, 0.236, 0.0),
                                  (-0.475, -0.17, 0.0)),
                     'frequency': 61.36},
                    {'vib_mode': ((0.0, 0.0, 0.707),
                                  (0.0, 0.0, 0.0),
                                  (0.0, 0.0, -0.707)),
                     'frequency': 199.94},
                    {'vib_mode': ((0.0, 0.0, -0.505),
                                  (0.0, 0.0, 0.7),
                                  (0.0, 0.0, -0.505)),
                     'frequency': 311.74}]
        self.assertEqual(qcout.data[1]['frequencies'], ans_freq)
        self.assertAlmostEqual(qcout.data[2]['energies'][0][1],
                               -5296.720741780598, 5)
        ans_scf_iter_ene = [[(-176.9147092199, 0.779),
                             (-156.8236033975, 0.115),
                             (-152.9396694452, 0.157),
                             (-183.2743425778, 0.138),
                             (-182.2994943574, 0.142),
                             (-181.990425533, 0.143),
                             (-182.1690180647, 0.142),
                             (-106.6454708618, 0.239),
                             (-193.8056267625, 0.0432),
                             (-193.0854096948, 0.0455),
                             (-194.6340538334, 0.0062),
                             (-194.6495072245, 0.00205),
                             (-194.6508787796, 0.000189),
                             (-194.6508984743, 2.18e-05),
                             (-194.6508986262, 2.17e-06)]]
        self.assertEqual(qcout.data[2]['scf_iteration_energies'],
                         ans_scf_iter_ene)

    def test_solvent_method(self):
        filename = os.path.join(test_dir, "thiophene_wfs_5_carboxyl.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(qcout.data[0]["solvent_method"], "NA")

        filename = os.path.join(test_dir, "qchem_energies", "hf_cosmo.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(qcout.data[0]["solvent_method"], "cosmo")

        filename = os.path.join(test_dir, "qchem_energies", "hf_pcm.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(qcout.data[0]["solvent_method"], "pcm")

    def test_failed_message(self):
        scf_file = os.path.join(test_dir, "hf.qcout")
        scf_qcout = QcOutput(scf_file)
        self.assertTrue(scf_qcout.data[0]['has_error'])
        self.assertEqual(scf_qcout.data[0]['errors'],
                         ['Bad SCF convergence',
                          'Molecular charge is not found',
                          'Geometry optimization failed'])
        geom_file = os.path.join(test_dir, "hf_opt_failed.qcout")
        geom_qcout = QcOutput(geom_file)
        self.assertTrue(geom_qcout.data[0]['has_error'])
        self.assertEqual(geom_qcout.data[0]['errors'],
                         ['Geometry optimization failed'])

    def test_abnormal_exit(self):
        no_reading_file = os.path.join(test_dir, "no_reading.qcout")
        no_reading_qcout = QcOutput(no_reading_file)
        self.assertTrue(no_reading_qcout.data[0]['has_error'])
        self.assertEqual(no_reading_qcout.data[0]['errors'],
                         ['Exit Code 134',
                          'Molecular charge is not found',
                          'No input text',
                          'Bad SCF convergence'])
        exit_code_134_file = os.path.join(test_dir, "exit_code_134.qcout")
        ec134_qcout = QcOutput(exit_code_134_file)
        self.assertTrue(ec134_qcout.data[0]['has_error'])
        self.assertEqual(ec134_qcout.data[0]['errors'],
                         ['Exit Code 134',
                          'Molecular charge is not found',
                          'Bad SCF convergence'])

    def test_chelp_and_mulliken_charges(self):
        filename = os.path.join(test_dir, 'chelpg_charges.qcout')
        qcout = QcOutput(filename)
        mulliken_charges = [0.393961, -0.281545, 0.066432, 0.019364, -0.186041,
                            -0.16007, 0.315659, 0.30631, 0.064257, 0.056438,
                            -0.17695, 0.16976, -0.13326, -0.131853, -0.178711,
                            0.163697, 0.170148, 0.143329, 0.152702, 0.152929,
                            0.170475, -0.451542, -0.441554, -0.709834,
                            -0.592718, 0.20506, 0.211043, 0.204389, 0.546173,
                            -0.414558, 0.346511]
        self.assertEqual(qcout.data[0]['charges']['mulliken'],
                         mulliken_charges)
        chelpg_charges = [0.399404, -0.277179, -0.057502, -0.110085, -0.07107,
                          -0.274987, 0.475781, 0.423117, -0.054079, -0.101424,
                          -0.05793, 0.115179, -0.116069, -0.10949, -0.06664,
                          0.161442, 0.135438, 0.158081, 0.125881, 0.125324,
                          0.115863, -0.425251, -0.42309, -0.602375, -0.458844,
                          0.140267, 0.139084, 0.139995, 0.698011, -0.487911,
                          0.341061]
        self.assertEqual(qcout.data[0]['charges']['chelpg'], chelpg_charges)

    @unittest.skipIf(sys.platform not in ["linux", 'darwin'],
                     "Skip unix file path test on Windows")
    def test_scr_dir(self):
        filename = os.path.join(test_dir, 'chelpg_charges.qcout')
        qcout = QcOutput(filename)
        self.assertEqual(qcout.data[0]['scratch_dir'],
                         "/Users/xiaohuiqu/scratch/qchem7101")

    def test_no_message_scf_opt_fail(self):
        so_failfile = os.path.join(test_dir, 'scf_opt_no_message_fail.qcout')
        so_failqcout = QcOutput(so_failfile)
        self.assertTrue(so_failqcout.data[0]['has_error'])
        self.assertEqual(so_failqcout.data[0]['errors'],
                         ['Exit Code 134',
                          'Molecular charge is not found',
                          'Bad SCF convergence',
                          'Geometry optimization failed'])
        o_failfile = os.path.join(test_dir, 'opt_fail_no_message.qcout')
        o_failqcout = QcOutput(o_failfile)
        self.assertEqual(o_failqcout.data[0]['errors'],
                         ['Geometry optimization failed'])
        s_failfile = os.path.join(test_dir, 'scf_no_message_fail.qcout')
        s_failqcout = QcOutput(s_failfile)
        self.assertEqual(s_failqcout.data[0]['errors'],
                         ['Exit Code 134',
                          'Molecular charge is not found',
                          'Bad SCF convergence'])
        so_successfile = os.path.join(test_dir,
                                      'thiophene_wfs_5_carboxyl.qcout')
        so_successqcout = QcOutput(so_successfile)
        self.assertFalse(so_successqcout.data[0]['has_error'])

    def test_negative_eigen(self):
        filename = os.path.join(test_dir, "negative_eigen.qcout")
        qcout = QcOutput(filename)
        self.assertTrue(qcout.data[0]['has_error'])
        self.assertEqual(qcout.data[0]["errors"],
                         ['Negative Eigen',
                          'Molecular charge is not found',
                          'Bad SCF convergence',
                          'Geometry optimization failed'])

    def test_insufficient_memory(self):
        filename = os.path.join(test_dir, "insufficient_memory.qcout")
        qcout = QcOutput(filename)
        self.assertTrue(qcout.data[0]['has_error'])
        self.assertEqual(qcout.data[0]['errors'],
                         ['Insufficient static memory',
                          'Molecular charge is not found',
                          'Bad SCF convergence',
                          'Geometry optimization failed'])

    def test_freq_seg_too_small(self):
        filename = os.path.join(test_dir, "freq_seg_too_small.qcout")
        qcout = QcOutput(filename)
        self.assertTrue(qcout.data[0]['has_error'])
        self.assertEqual(qcout.data[0]['errors'],
                         ['Freq Job Too Small',
                          'Exit Code 134'])

    def test_not_enough_total_memory(self):
        filename = os.path.join(test_dir, "not_enough_total_memory.qcout")
        qcout = QcOutput(filename)
        self.assertTrue(qcout.data[1]['has_error'])
        self.assertEqual(qcout.data[1]["errors"],
                         ['Not Enough Total Memory',
                          'Exit Code 134'])

    def test_killed(self):
        filename = os.path.join(test_dir, "killed.qcout")
        qcout = QcOutput(filename)
        self.assertFalse(qcout.data[0]["has_error"])
        self.assertTrue(qcout.data[1]["has_error"])
        self.assertEqual(qcout.data[1]["errors"],
                         ['Killed',
                          'Molecular charge is not found',
                          'Bad SCF convergence'])

    def test_gdm_scf(self):
        filename = os.path.join(test_dir, "gmd_scf.qcout")
        qcout = QcOutput(filename)
        self.assertTrue(qcout.data[0]['has_error'])
        self.assertEqual(qcout.data[0]['errors'],
                         ['Exit Code 134',
                          'Bad SCF convergence',
                          'Geometry optimization failed'])
        self.assertEqual(len(qcout.data[0]['scf_iteration_energies']), 2)
        self.assertEqual(len(qcout.data[0]['scf_iteration_energies'][-1]), 192)
        self.assertAlmostEqual(qcout.data[0]['scf_iteration_energies'][-1][-1][0],
                               -1944.945908459, 5)

    def test_crazy_scf_values(self):
        filename = os.path.join(test_dir, "crazy_scf_values.qcout")
        qcout = QcOutput(filename)
        ans = [(-28556254.06737586, 6.49e-06),
               (-28556254.067382727, 9.45e-06),
               (-28556254.067382865, 6.14e-06)]
        self.assertEqual(qcout.data[0]["scf_iteration_energies"][-1][-3:], ans)

    def test_crowd_gradient_number(self):
        filename = os.path.join(test_dir, "crowd_gradient_number.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(qcout.data[0]['gradients'][0]['gradients'],
                         [(-0.0307525, 0.0206536, -0.0396255),
                          (0.0008938, -0.000609, 0.0082746),
                          (0.042143, -0.0240514, 0.0380298),
                          (-0.0843578, 0.0002757, 0.0884924),
                          (0.0356689, -0.0444656, -0.0710646),
                          (-0.0190554, -0.0308886, -0.0297994),
                          (0.0470543, -0.0263915, -0.0690973),
                          (-0.0297801, 0.0296872, -0.0104344),
                          (0.0504581, -0.0014272, 0.0262245),
                          (-0.0927323, 0.0750046, 0.0128003),
                          (0.0183242, -0.0084638, 0.0127388),
                          (-0.0083989, 0.0111579, -0.0002461),
                          (-0.0316941, 267.34455, 878.3493251),
                          (0.017459, 0.0487124, -0.0276365),
                          (-0.3699134, 0.0110442, 0.0260809),
                          (0.363931, 0.24044, 0.5192852),
                          (0.026669, -0.0284192, -0.0347528),
                          (0.0047475, 0.0049706, 0.0148794),
                          (-0.077804, 0.003402, 0.000852),
                          (-6772.1697035, -267.4471902, -878.585931),
                          (-0.0029556, -0.0616073, -0.0180577),
                          (-0.0001915, 0.0021213, 0.0006193),
                          (0.0320436, -0.0073456, -0.01509),
                          (0.0155112, -0.0035725, 0.0015675),
                          (-0.0034309, 0.0170739, 0.0074455),
                          (-0.0088735, -0.0129874, 0.0092329),
                          (-0.0271963, -0.0258714, 0.0246954),
                          (0.0025065, 0.0062934, 0.0209733),
                          (0.0152829, -0.0080239, -0.018902),
                          (0.0461304, 0.0071952, 0.0012227),
                          (-0.0272755, -0.0280053, 0.0325455),
                          (0.0122118, 0.027816, -0.0167773),
                          (0.0168893, -0.0014211, 0.0039917),
                          (-0.0048723, 0.0026667, -0.0159952),
                          (-0.1840467, -0.1425887, -0.3235801),
                          (0.015975, -0.0922797, 0.0640925),
                          (0.0267234, 0.1031154, -0.0299014),
                          (-0.0175591, 0.0081813, -0.0165425),
                          (0.0119225, 0.0113174, 0.0154056),
                          (0.0138491, 0.0083436, 0.0188022),
                          (-0.0151146, -0.0015971, -0.0054462)])

    def test_nbo_charges(self):
        filename = os.path.join(test_dir, "quinoxaline_anion.qcout")
        qcout = QcOutput(filename)
        ans = [-0.29291, -0.29807, 0.12715, 0.12715, -0.29807, -0.29291,
               0.21284, 0.22287, 0.22287, 0.21284, -0.10866, -0.10866,
               0.19699, -0.5602, -0.5602, 0.19699]
        self.assertEqual(qcout.data[0]["charges"]["nbo"], ans)
        filename = os.path.join(test_dir, "tfsi_nbo.qcout")
        qcout = QcOutput(filename)
        ans = [2.2274, 2.23584, -0.94183, -0.94575, -0.94719, -0.9423,
               0.86201, 0.85672, -0.35698, -0.35373, -0.35782, -0.35647,
               -0.35646, -0.35787, -1.26555]
        self.assertEqual(qcout.data[0]["charges"]["nbo"], ans)
        filename = os.path.join(test_dir, "crowd_nbo_charges.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(
            qcout.data[0]["charges"]["nbo"],
            [-0.33917, -0.6104, -0.15912, -0.17751, -0.61817, -0.3357, 0.24671,
             0.19942, 0.19325, 0.2362, 0.23982, 0.21985, 0.2305, 0.20444,
             0.23179, 0.20491, 0.85965, -0.59655, -0.59561, -0.14789, -0.13859,
             -0.32712, -0.33359, 0.21602, 0.22383, 0.2123, 0.22759, 0.2507,
             0.20098, 0.18631, 0.24945, 0.19709, 0.20274, -0.34831, -0.56307,
             -0.14572, -0.1431, -0.55866, -0.3572, 0.22695, 0.21983, 0.1963,
             0.20977, 0.22298, 0.20875, 0.21081, 0.19586, 0.24708, 0.20067,
             -0.34288, -0.55793, -0.16806, -0.15609, -0.56464, -0.34695,
             0.22555, 0.20417, 0.206, 0.20825, 0.22409, 0.25415, 0.20977,
             0.18976, 0.24647, 0.1993, -0.33605, -0.59395, -0.15985, -0.18024,
             -0.60646, -0.32742, 0.22909, 0.19347, 0.21872, 0.2203, 0.23518,
             0.25185, 0.23523, 0.18666, 0.22737, 0.2205, -0.35902, -0.56138,
             -0.14552, -0.14903, -0.55491, -0.3493, 0.22826, 0.21789, 0.19075,
             0.20898, 0.21343, 0.21715, 0.20794, 0.19695, 0.2429, 0.18482,
             -0.33943, -0.55659, -0.16437, -0.14503, -0.56155, -0.34131,
             0.22339, 0.20483, 0.19376, 0.23395, 0.20784, 0.2096, 0.21945,
             0.19192, 0.23089, 0.20493, -0.32963, -0.56949, -0.1446, -0.15244,
             -0.55482, -0.34848, 0.22802, 0.20471, 0.19704, 0.20744, 0.22332,
             0.2206, 0.20734, 0.18871, 0.22907, 0.20741, -0.33856, -0.564,
             -0.16575, -0.17422, -0.56032, -0.3426, 0.22585, 0.20169, 0.20529,
             0.20836, 0.21329, 0.25353, 0.23374, 0.19306, 0.23582, 0.20196,
             -0.34069, -0.56522, -0.17228, -0.17503, -0.55505, -0.34264,
             0.22696, 0.19604, 0.20515, 0.23964, 0.2437, 0.2111, 0.21204,
             0.19975, 0.2347, 0.18835, -0.34324, -0.55184, -0.16086, -0.15907,
             -0.56319, -0.3384, 0.23866, 0.19808, 0.19728, 0.20205, 0.24698,
             0.21416, 0.20398, 0.20475, 0.2265, 0.20141, -0.34339, -0.56344,
             -0.14955, -0.14878, -0.55906, -0.34506, 0.23937, 0.20027, 0.19671,
             0.2085, 0.21693, 0.22164, 0.20863, 0.20703, 0.22889, 0.1916])

    def test_simple_aimd(self):
        filename = os.path.join(test_dir, "h2o_aimd.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(len(qcout.data[0]["molecules"]), 11)

    def test_homo_lumo(self):
        filename = os.path.join(test_dir, "quinoxaline_anion.qcout")
        qcout = QcOutput(filename)
        for a, b in zip(qcout.data[0]["HOMO/LUMOs"][-1],
                        [1.00682120282, 2.80277253758]):
            self.assertAlmostEqual(a, b, 5)
        filename = os.path.join(test_dir, "qchem_energies", "hf_ccsd(t).qcout")
        qcout = QcOutput(filename)
        self.assertArrayAlmostEqual(qcout.data[0]["HOMO/LUMOs"],
                              [[-17.741823053011334, 5.224585929721129],
                               [-17.741823053011334, 5.224585929721129]], 4)
        filename = os.path.join(test_dir, "crowd_gradient_number.qcout")
        qcout = QcOutput(filename)
        self.assertArrayAlmostEqual(
            qcout.data[0]["HOMO/LUMOs"], [[-5.741602245683116,
                                                        -4.544301303455358],
                                                       [-4.9796834642654515,
                                                        -4.2993988379996795], [-4.761992383860404, -3.8095939070883236]], 4)

    def test_bsse(self):
        filename = os.path.join(test_dir, "bsse.qcout")
        qcout = QcOutput(filename)
        self.assertAlmostEqual(qcout.data[0]["bsse"], -0.164210762949, 5)
        self.assertEqual(qcout.data[0]["jobtype"], "bsse")

    def test_hirshfeld_charge(self):
        filename = os.path.join(test_dir, "hirshfeld_population.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(qcout.data[0]["charges"]["hirshfeld"],
                         [-0.286309, 0.143134, 0.143176])
        self.assertFalse(qcout.data[0]["has_error"])

    def test_ghost_atoms(self):
        filename = os.path.join(test_dir, "ghost_atoms.qcout")
        qcout = QcOutput(filename)
        elements = [a.specie.symbol for a in qcout.data[-1]["molecules"][-1].sites]
        self.assertEqual(elements, ['O', 'H', 'H', 'C', 'H', 'H', 'H', 'H'])
        filename = os.path.join(test_dir, "MgBF4_b_overalpped.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(qcout.data[0]["input"].ghost_atoms, [0])

    def test_final_energy(self):
        filename = os.path.join(test_dir, "thiophene_wfs_5_carboxyl.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(qcout.final_energy, -20180.15180854295)

    def test_final_structure(self):
        filename = os.path.join(test_dir, "thiophene_wfs_5_carboxyl.qcout")
        qcout = QcOutput(filename)
        ans = '''Full Formula (H4 C5 S1 O2)
Reduced Formula: H4C5SO2
Charge = -1, Spin Mult = 2
Sites (12)
0 C     0.194695    -0.158362    -0.001887
1 C    -0.535373    -1.381241    -0.001073
2 C    -1.927071    -1.199274    -0.000052
3 C    -2.332651     0.131916     0.000329
4 S    -0.942111     1.224916    -0.001267
5 H    -0.038260    -2.345185    -0.001256
6 H    -2.636299    -2.025939     0.000620
7 H    -3.339756     0.529895     0.001288
8 C     1.579982     0.071245    -0.002733
9 O     2.196383     1.165675    -0.000178
10 O     2.352341    -1.114671     0.001634
11 H     3.261096    -0.769470     0.003158'''
        self.assertEqual(qcout.final_structure.__str__(), ans)

    def test_time_nan_values(self):
        filename = os.path.join(test_dir, "time_nan_values.qcout")
        qcout = QcOutput(filename)
        self.assertFalse(qcout.data[0]["has_error"])

    def test_pcm_solvent_deprecated(self):
        filename = os.path.join(test_dir, "pcm_solvent_deprecated.qcout")
        qcout = QcOutput(filename)
        self.assertTrue(qcout.data[-1]["has_error"])
        ans = ['pcm_solvent deprecated',
               'Molecular charge is not found',
               'No input text',
               'Bad SCF convergence']
        self.assertEqual(qcout.data[-1]["errors"], ans)

    def test_qc43_batch_job(self):
        filename = os.path.join(test_dir, "qchem43_batch_job.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(len(qcout.data), 2)
        self.assertEqual(len(qcout.data[0]["scf_iteration_energies"][0]), 22)
        self.assertTrue("pcm_solvent deprecated" in qcout.data[1]["errors"])

    def test_output_file_wierd_encoding(self):
        filename = os.path.join(test_dir, "ferrocenium_1pos.qcout")
        qcout = QcOutput(filename)
        self.assertFalse(qcout.data[1]["has_error"])
        self.assertEqual(qcout.data[1]["frequencies"][0]["frequency"], -157.11)

    def test_homo_lumo_nan_values(self):
        filename = os.path.join(test_dir, "homo_lumo_nan_values.qcout")
        qcout = QcOutput(filename)
        self.assertTrue(qcout.data[0]["has_error"])

    def test_ordinal_not_in_range(self):
        filename = os.path.join(test_dir, "ordinal_not_in_range.qcout.gz")
        qcout = QcOutput(filename)
        self.assertEqual(len(qcout.data), 1)

    def test_aux_mpi_time_in_the_end_of_job(self):
        filename = os.path.join(test_dir, "aux_mpi_time_mol.qcout")
        qcout = QcOutput(filename)
        self.assertEqual(len(qcout.data), 2)

    def test_parse_opt_contraint(self):
        filename = os.path.join(test_dir, "pt_dft_180.0.qcout")
        qcout = QcOutput(filename)
        qcin = qcout.data[-1]['input']
        qcin_ans = '''$molecule
 0  1
 S           1.82267924       -1.19997629        0.28714109
 C           3.20006180       -0.17260711        0.06528466
 C           2.82980603        1.10216298       -0.25610036
 C           1.41909100        1.26345446       -0.34254814
 C           0.71738150        0.10901545       -0.08456145
 H           0.93627498        2.19419272       -0.61095402
 C          -0.71741859       -0.10899254       -0.08455524
 S          -1.82328469        1.20374179       -0.44105740
 C          -1.41912820       -1.26343144        0.17343142
 C          -3.19922829        0.16690023       -0.25767458
 C          -2.82941826       -1.10493701        0.07562280
 H          -3.53750269       -1.90709774        0.23645949
 H           4.19429620       -0.57452886        0.18632814
 H           3.53860725        1.89960515       -0.43610218
 H          -4.19239866        0.56181917       -0.40716131
 H          -0.93481970       -2.20399421        0.40193462
$end


$rem
         jobtype = opt
        exchange = b3lyp
           basis = 6-31++g**
  max_scf_cycles = 75
      mem_static = 100
       mem_total = 1500
$end


$opt
CONSTRAINT
tors 4 5 7 9 180.0
ENDCONSTRAINT
$end

'''
        self.assertEqual(str(qcin), qcin_ans)
        constraint = qcin.params['opt']['CONSTRAINT']
        constraint_ans = [['tors', 4, 5, 7, 9, 180.0]]
        self.assertEqual(constraint, constraint_ans)

        stre_text = """CONSTRAINT
stre 70 9 3.795
stre 13 44 3.656
ENDCONSTRAINT"""
        stre_d = qcin._parse_opt(stre_text.split('\n'))
        qctask = QcTask(mol, exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*",
                        optional_params={"opt": stre_d})
        stre_text_2 = "\n".join(qctask._format_opt())
        self.assertEqual(stre_text_2, stre_text)

class TestQcNucVeloc(PymatgenTest):

    def test_parse(self):
        qcnv = QcNucVeloc(os.path.join(test_dir, "qc_aimd", "NucVeloc.velocities"))
        self.assertEqual(len(qcnv.step_times), 302)
        self.assertEqual(qcnv.step_times[-1], 14.56168)
        self.assertEqual(len(qcnv.velocities), 302)
        self.assertEqual(qcnv.velocities[0][0],
                         (1.42192e-05, 6.65659e-05, 5.22453e-05))
        self.assertEqual(qcnv.velocities[0][-1],
                         (-6.52039e-06, -0.000213074, -0.000769596))
        self.assertEqual(qcnv.velocities[-1][0],
                         (8.976072e-05, 9.455759e-06, -0.0002397046))
        self.assertEqual(qcnv.velocities[-1][-1],
                         (9.052722e-05, 0.001113288, -0.0009176628))


if __name__ == "__main__":
    unittest.main()
