"""
Helper methods for generating gw input / and work flows.
"""

from __future__ import division

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
from pymatgen.core.units import Ha_to_eV, eV_to_Ha


class SplineInputError(Exception):
    def __init__(self, msg):
        self.msg = msg


def now():
    """
    helper to return a time string
    """
    return time.strftime("%H:%M:%S %d/%m/%Y")

import string
import random


def id_generator(size=8, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def s_name(structure):
    name_ = str(structure.composition.reduced_formula) #+ '_' + str(structure.item)
    return name_


def clean(string, uppercase=False):
    """
    helper to clean up an input string
    """
    if uppercase:
        return string.strip().upper()
    else:
        return string.strip().lower()


def get_derivatives(xs, ys):
    """
    return the derivatives of y(x) at the points x
    if scipy is available a spline is generated to calculate the derivatives
    if scipy is not available the left and right slopes are calculated, if both exist the average is returned
    """
    try:
        if len(xs) < 4:
            er = SplineInputError('too few data points')
            raise er
        from scipy.interpolate import UnivariateSpline
        spline = UnivariateSpline(xs, ys)
        d = spline.derivative(1)(xs)
    except (ImportError, SplineInputError):
        d = []
        m, left, right = 0, 0, 0
        for n in range(0, len(xs), 1):
            try:
                left = (ys[n] - ys[n-1]) / (xs[n] - xs[n-1])
                m += 1
            except IndexError:
                pass
            try:
                right = (ys[n+1] - ys[n]) / (xs[n+1] - xs[n])
                m += 1
            except IndexError:
                pass
            d.append(left + right / m)
    return d


def reciprocal(x, a, b, n):
    """
    reciprocal function to the power n to fit convergence data
    """
    import numpy as np
    if n < 1:
        n = 1
    elif n > 5:
        n = 5
    #print a, b, x
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / x_v ** n)
        y = np.array(y_l)
    else:
        y = a + b / x ** n
    #print y
    #print type(y)
    return y


def exponential(x, a, b, n):
    """
    exponential function base n to fit convergence data
    """
    import numpy as np
    if n < 1.001:
        n = 1.001
    elif n > 1.2:
        n = 1.2
    #print a, b, x
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b * n ** -x_v)
        y = np.array(y_l)
    else:
        y = a + b * n ** -x
    #print y
    #print type(y)
    return y


def p0exp(xs, ys):
    n0 = 1.1
    b0 = (n0 ** -xs[-1] - n0 ** -xs[0]) / (ys[-1] - ys[0])
    a0 = ys[0] - b0 * n0 ** -xs[0]
    return [a0, b0, n0]


def double_reciprocal(x, a, b, c):
    """
    reciprocal function to the power n to fit convergence data
    """
    import numpy as np
    #print a, b, x
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / x_v + c / x_v ** 2)
        y = np.array(y_l)
    else:
        y = a + b / x + c / x ** 2
    #print y
    #print type(y)
    return y


def p0reci(xs, ys):
    """
    predictor for first guess for reciprocal
    """
    a0 = ys[len(ys) - 1]
    b0 = ys[0]*xs[0] - a0*xs[0]
    return [a0, b0, 1]


def test_conv(xs, ys, name, tol=0.0001):
    """
    test it and at which x_value dy(x)/dx < tol for all x >= x_value, conv is true is such a x_value exists.
    """
    conv = False
    x_value = float('inf')
    y_value = None
    n_value = None
    popt = [None, None, None]
    if len(xs) > 2:
        ds = get_derivatives(xs[0:len(ys)], ys)
        try:
            import numpy as np
            from scipy.optimize import curve_fit
            if None not in ys:
                popt, pcov = curve_fit(exponential, xs, ys, p0exp(xs, ys), maxfev=8000)
                # todo print this to file via a method in helper, as dict
                f = open(name+'.fitdat', mode='a')
                f.write('{')
                f.write('"popt": ' + str(popt) + ', ')
                f.write('"pcov": ' + str(pcov) + ', ')
                f.write('"data": [')
                for n in range(0, len(ys), 1):
                    f.write('[' + str(xs[n]) + ' ' + str(ys[n]) + ']')
                f.write(']}\n')
                f.close()
               # print 'plot ', popt[0], ' + ', popt[1], "/x**", popt[2], ', "'+name+'.convdat"'
              #  print 'plot ', popt[0], ' + ', popt[1], "/x", popt[2], '/x**2, "'+name+'.convdat"'
                id = id_generator()
                print 'plot ', popt[0], ' + ', popt[1], "* ", popt[2], " ** -x, '", 'convdat.'+id, "'"
                f = open('convdat.'+str(id), mode='w')
                for n in range(0, len(ys), 1):
                    f.write(str(xs[n]) + ' ' + str(ys[n]) + '\n')
                f.write('\n')
                f.close()
        except ImportError:
            popt, pcov = None, None
        for n in range(0, len(ds), 1):

            if tol < 0:
                if popt[0] is not None:
                    test = abs(popt[0] - ys[n])
                else:
                    test = float('inf')
            else:
                test = abs(ds[n])

            if test < abs(tol):
                conv = True
                if xs[n] < x_value:
                    x_value = xs[n]
                    y_value = ys[n]
                    n_value = n
            else:
                conv = False
                x_value = float('inf')
        if n_value is None:
            return [conv, x_value, y_value, n_value, popt[0], None]
        else:
            return [conv, x_value, y_value, n_value, popt[0], ds[n_value]]
    else:
        return [conv, x_value, y_value, n_value, popt[0], None]


def expand_tests(tests, level):
    from pymatgen.io.gwwrapper.codeinterfaces import get_all_ecuteps, get_all_nbands
    new_tests = copy.deepcopy(tests)
    for test in tests.keys():
        if test in get_all_ecuteps():
            ec = str(test)
            ec_range = tests[ec]['test_range']
            ec_step = ec_range[-1] - ec_range[-2]
            if int(level / 2) == level / 2:
                print 'new ec wedge'
                # even level of grid extension > new ec wedge
                new_ec_range = (ec_range[-1] + int(level / 2 * ec_step),)
            else:
                print 'new nb wedge'
                # odd level of grid extension > new nb wedge
                extension = tuple(range(ec_range[-1] + ec_step, ec_range[-1] + (1 + int((level - 1) / 2)) * ec_step, ec_step))
                new_ec_range = ec_range + extension
            new_tests[ec].update({'test_range': new_ec_range})
        if test in get_all_nbands():
            nb = str(test)
            nb_range = tests[nb]['test_range']
            nb_step = nb_range[-1] - nb_range[-2]
            print nb_step
            if int(level / 2) == level / 2:
                # even level of grid extension > new ec wedge
                extension = tuple(range(nb_range[-1] + nb_step, nb_range[-1] + (1 + int(level / 2)) * nb_step, nb_step))
                new_nb_range = nb_range + extension
            else:
                # odd level of grid extension > new nb wedge
                new_nb_range = (nb_range[-1] + int((level + 1) / 2 * nb_step),)
            new_tests[nb].update({'test_range': new_nb_range})
    print new_tests
    return new_tests


def print_gnuplot_header(filename, title='', mode='convplot', filetype='jpeg'):
    xl = 'set xlabel "nbands"\n'
    yl = 'set ylabel "ecuteps (eV)"\n'
    zl = 'set zlabel "gap (eV)"\n'
    if mode == 'convplot':
        f = open(filename, mode='a')
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
        print 'Problems reading ', filename
        full_res = {'grid': 0, 'all_done': False}
    except (OSError, IOError):
        full_res = {'grid': 0, 'all_done': False}
    return full_res


def is_converged(hartree_parameters, structure, return_values=False):
    filename = s_name(structure) + ".conv_res"
    try:
        f = open(filename, mode='r')
        conv_res = ast.literal_eval(f.read())
        f.close()
        converged = conv_res['control']['nbands']
    except (IOError, OSError):
        if return_values:
            print 'Inputfile ', filename, ' not found, the convergence calculation did not finish properly' \
                                          ' or was not parsed ...'
        converged = False
        return converged
    if return_values and converged:
        if hartree_parameters:
            conv_res['values']['ecuteps'] = 4 * math.ceil(conv_res['values']['ecuteps'] * eV_to_Ha / 4)
        return conv_res['values']
    else:
        return converged


def store_conv_results(name, folder):
    print "| Storing results for %s" % name
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
