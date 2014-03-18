"""
Helper modules for generating gw input / workflows.
"""

from __future__ import division

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "Oct 23, 2013"

import time


class SplineInputError(Exception):
    def __init__(self, msg):
        self.msg = msg


def now():
    """
    helper to return a time string
    """
    return time.strftime("%H:%M:%S %d/%m/%Y")


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
            er = SplineInputError('test')
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


def reciprocal(x, a, b):
    n = 2
    print a, b, x
    if isinstance(x, list):
        y = []
        for x_v in x:
            y.append(a + b / x_v ** n)
    else:
        y = a + b / x ** n
    print y
    return y


#def reciprocal(x, a, b):
#    print a, b, x
#    y = a + b / x
#    print y
#    return y


def p0reci(xs, ys):
    """
    predictor for first gues
    """
    a0 = ys[len(ys) - 1]
    b0 = ys[0]*xs[0] - a0*xs[0]
    return [a0, b0]


def test_conv(xs, ys, tol=0.0001):
    """
    test it and at which x_value dy(x)/dx < tol for all x >= x_value, conv is true is such a x_value exists.
    """
    conv = False
    x_value = float('inf')
    y_value = None
    n_value = None
    popt = None
    if len(xs) > 1:
        ds = get_derivatives(xs[0:len(ys)], ys)
        try:
            import numpy as np
            from scipy.optimize import curve_fit
            #print 'xs    ', xs
            #print 'ys    ', ys
            if None not in ys:
                popt, pcov = curve_fit(reciprocal, xs, ys, p0reci(xs, ys))
                print 'plot ', popt[0], ' + ', popt[1], "/x**2, '-' w p"
                for n in range(0, len(ys), 1):
                    print xs[n], ys[n]
                print 'e'
        except ImportError:
            popt, pcov = None, None
        for n in range(0, len(ds), 1):

            if tol < 0:
                if popt is not None:
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

    return [conv, x_value, y_value, n_value, popt[0]]


def print_gnuplot_header(filename, title='', mode='convplot', filetype='jpeg'):
    xl = 'set xlabel "nbands"\n'
    yl = 'set ylabel "encutgw (eV)"\n'
    zl = 'set zlabel "gap (eV)"\n'
    if mode == 'convplot':
        f = open(filename, mode='a')
        f.write('set terminal '+filetype+'\n')
        f.write('set title "'+title+'"\n')
        f.write(xl)
        f.write(yl)
        f.write(zl)
        f.close()