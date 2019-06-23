# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
Helper methods for generating gw input / and work flows.
"""

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "May 2014"

import time
import os
import ast
import copy
import math
import shutil
import numpy as np
from pymatgen.core.units import eV_to_Ha
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import OxidationStateRemovalTransformation, \
    PrimitiveCellTransformation, SupercellTransformation


def now():
    """
    helper to return a time string
    """
    return time.strftime("%H:%M:%S %d/%m/%Y")


def read_extra_abivars():
    ea = {}
    if os.path.isfile('extra_abivars'):
        f = open('extra_abivars')
        ea = ast.literal_eval(f.read())
        if not isinstance(ea, dict):
            raise RuntimeError
    return ea


def refine_structure(structure, symprec=1e-3):
    remove_ox = OxidationStateRemovalTransformation()
    structure = remove_ox.apply_transformation(structure)
    sym_finder = SpacegroupAnalyzer(structure=structure, symprec=symprec)
    structure = sym_finder.get_refined_structure()
    get_prim = PrimitiveCellTransformation()
    structure = get_prim.apply_transformation(structure)
    m = structure.lattice.matrix
    x_prod = np.dot(np.cross(m[0], m[1]), m[2])
    if x_prod < 0:
        print(x_prod)
        trans = SupercellTransformation(((1, 0, 0), (0, 0, 1), (0, 1, 0)))
        structure = trans.apply_transformation(structure)
        m = structure.lattice.matrix
        x_prod = np.dot(np.cross(m[0], m[1]), m[2])
        print(x_prod)
        if x_prod < 0:
            raise RuntimeError
    return structure


def s_name(structure):
    if os.path.isfile('old'):
        name_ = str(structure.composition.reduced_formula)
    else:
        name_ = str(structure.composition.reduced_formula) + '_' + str(structure.item)
    return name_


def clean(some_string, uppercase=False):
    """
    helper to clean up an input string
    """
    if uppercase:
        return some_string.strip().upper()
    else:
        return some_string.strip().lower()


def expand(tests, level):
    from abipy.gw.codeinterfaces import get_all_ecuteps, get_all_nbands
    new_tests = copy.deepcopy(tests)
    for test in tests.keys():
        if test in get_all_ecuteps():
            ec = str(test)
            ec_range = tests[ec]['test_range']
            ec_step = ec_range[-1] - ec_range[-2]
            if int(level / 2) == level / 2:
                print('new ec wedge')
                # even level of grid extension > new ec wedge
                new_ec_range = (ec_range[-1] + int(level / 2 * ec_step),)
            else:
                print('new nb wedge')
                # odd level of grid extension > new nb wedge
                extension = tuple(range(ec_range[-1] + ec_step, ec_range[-1] + (1 + int((level - 1) / 2)) * ec_step, ec_step))
                new_ec_range = ec_range + extension
            new_tests[ec].update({'test_range': new_ec_range})
        if test in get_all_nbands():
            nb = str(test)
            nb_range = tests[nb]['test_range']
            nb_step = nb_range[-1] - nb_range[-2]
            print(nb_step)
            if int(level / 2) == level / 2:
                # even level of grid extension > new ec wedge
                extension = tuple(range(nb_range[-1] + nb_step, nb_range[-1] + (1 + int(level / 2)) * nb_step, nb_step))
                new_nb_range = nb_range + extension
            else:
                # odd level of grid extension > new nb wedge
                new_nb_range = (nb_range[-1] + int((level + 1) / 2 * nb_step),)
            new_tests[nb].update({'test_range': new_nb_range})
    print(new_tests)
    return new_tests


def print_gnuplot_header(filename, title='', mode='convplot', filetype='jpeg'):
    xl = 'set xlabel "nbands"\n'
    yl = 'set ylabel "ecuteps (eV)"\n'
    zl = 'set zlabel "gap (eV)"\n'
    if mode == 'convplot':
        f = open(filename, mode='a')
        if filetype is not None:
            f.write('set terminal '+filetype+'\n')
        f.write('set title "'+title+'"\n')
        f.write(xl)
        f.write(yl)
        f.write(zl)
        f.close()


def read_grid_from_file(filename):
    """
    Read the results of a full set of calculations from file
    """
    try:
        f = open(filename, mode='r')
        full_res = ast.literal_eval(f.read())
        f.close()
    except SyntaxError:
        print('Problems reading ', filename)
        full_res = {'grid': 0, 'all_done': False}
    except (OSError, IOError):
        full_res = {'grid': 0, 'all_done': False}
    return full_res


def is_converged(hartree_parameters, structure, return_values=False):
    filename = s_name(structure) + ".conv_res"
    to_return = {}
    try:
        f = open(filename, mode='r')
        conv_res = ast.literal_eval(f.read())
        f.close()
        converged = True if True in conv_res['control'].values() else False
    except (IOError, OSError, ValueError):
        if return_values:
            print('Inputfile ', filename, ' not found, the convergence calculation did not finish properly' \
                                          ' or was not parsed ...')
        converged = False
        return converged
    if return_values and converged:
        if hartree_parameters:
            try:
                conv_res['values']['ecut'] = 4 * math.ceil(conv_res['values']['ecut'] * eV_to_Ha / 4)
            except (KeyError, ArithmeticError, FloatingPointError, SyntaxError) as ex:
                print('exception in is_converged %s', ex.message)
            try:
                conv_res['values']['ecuteps'] = 4 * math.ceil(conv_res['values']['ecuteps'] * eV_to_Ha / 4)
            except (KeyError, ArithmeticError, FloatingPointError, SyntaxError) as ex:
                print('exception in is_converged %s', ex.message)
        for k in conv_res['values'].keys():
            if conv_res['values'][k] != 0 and conv_res['values'][k] != np.inf:
                to_return.update({k: conv_res['values'][k]})
        return to_return
    else:
        return converged


def store_conv_results(name, folder):
    print("| Storing results for %s" % name)
    if not os.path.isdir(folder):
        os.mkdir(folder)
    shutil.copy(name+'.full_res', os.path.join(folder, name+'.full_res'))
    for data_file in ['conv_res', 'log', 'conv.log', 'str', 'fitdat', 'convdat', 'data']:
        try:
            os.rename(name+'.'+data_file, os.path.join(folder, name+'.'+data_file))
        except OSError:
            pass


def add_gg_gap(structure):
    structure.vbm_l = "G"
    structure.cbm_l = "G"
    structure.cbm = (0.0, 0.0, 0.0)
    structure.vbm = (0.0, 0.0, 0.0)
    return structure
