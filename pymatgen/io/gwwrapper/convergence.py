"""
function for calculating the convergence of an x, y data set
main function:
test_conv(xs, ys, name, tol)
tries to fit multiple functions to the x, y data
calculates which function fits best
for tol < 0
returns the x value for which y is converged within tol of the assymtotic value
for tol > 0
returns the x_value for which dy(x)/dx < tol for all x >= x_value, conv is true is such a x_value exists
for the best fit a gnuplot line is printed plotting the data, the function and the assymthotic value
"""

from __future__ import division

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "May 2014"

import string
import random


def id_generator(size=8, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


class SplineInputError(Exception):
    def __init__(self, msg):
        self.msg = msg


def get_derivatives(xs, ys, fd=False):
    """
    return the derivatives of y(x) at the points x
    if scipy is available a spline is generated to calculate the derivatives
    if scipy is not available the left and right slopes are calculated, if both exist the average is returned
    putting fd to zero always returns the finite difference slopes
    """
    try:
        if fd:
            raise SplineInputError('no spline wanted')
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


"""
functions used in the fitting procedure, with initial guesses
"""


def reciprocal(x, a, b, n):
    """
    reciprocal function to the power n to fit convergence data
    """
    import numpy as np
    if n < 1:
        n = 1
    elif n > 5:
        n = 5
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / x_v ** n)
        y = np.array(y_l)
    else:
        y = a + b / x ** n
    return y


def p0_reciprocal(xs, ys):
    """
    predictor for first guess for reciprocal
    """
    a0 = ys[len(ys) - 1]
    b0 = ys[0]*xs[0] - a0*xs[0]
    return [a0, b0, 1]


def exponential(x, a, b, n):
    """
    exponential function base n to fit convergence data
    """
    import numpy as np
    if n < 1.000001:
        n = 1.000001
        #print n
    elif n > 1.2:
        n = 1.2
        #print n
    if b < -10:
        b = -10
        #print b
    elif b > 10:
        b = 10
        #print b
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


def p0_exponential(xs, ys):
    n0 = 1.005
    b0 = (n0 ** -xs[-1] - n0 ** -xs[1]) / (ys[-1] - ys[1])
    a0 = ys[1] - b0 * n0 ** -xs[1]
    #a0 = ys[-1]
    #b0 = (ys[0] - a0) / n0 ** xs[0]
    return [a0, b0, n0]


def single_reciprocal(x, a, b, c):
    """
    reciprocal function to fit convergence data
    """
    import numpy as np
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / (x_v - c))
        y = np.array(y_l)
    else:
        y = a + b / (x - c)
    return y


def p0_single_reciprocal(xs, ys):
    c = 1
    b = (1/(xs[-1] - c)-1/(xs[1] - c)) / (ys[-1] - ys[1])
    a = ys[1] - b / (xs[1] - c)
    return [a, b, c]


def simple_reciprocal(x, a, b):
    """
    reciprocal function to fit convergence data
    """
    import numpy as np
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / x_v)
        y = np.array(y_l)
    else:
        y = a + b / x
    return y


def p0_simple_reciprocal(xs, ys):
    c = 0
    b = (1/(xs[-1] - c) - 1/(xs[1] - c)) / (ys[-1] - ys[1])
    a = ys[1] - b / (xs[1] - c)
    return [a, b]


def simple_2reciprocal(x, a, b):
    """
    reciprocal function to fit convergence data
    """
    import numpy as np
    c = 2
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / x_v ** c)
        y = np.array(y_l)
    else:
        y = a + b / x ** c
    return y


def p0_simple_2reciprocal(xs, ys):
    c = 2
    b = (1/xs[-1]**c - 1/xs[1]**c) / (ys[-1] - ys[1])
    a = ys[1] - b / xs[1]**c
    return [a, b]


def simple_4reciprocal(x, a, b):
    """
    reciprocal function to fit convergence data
    """
    import numpy as np
    c = 4
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / x_v ** c)
        y = np.array(y_l)
    else:
        y = a + b / x ** c
    return y


def p0_simple_4reciprocal(xs, ys):
    c = 4
    b = (1/xs[-1]**c - 1/xs[1]**c) / (ys[-1] - ys[1])
    a = ys[1] - b / xs[1]**c
    return [a, b]


def simple_5reciprocal(x, a, b):
    """
    reciprocal function to fit convergence data
    """
    import numpy as np
    c = 0.5
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / x_v ** c)
        y = np.array(y_l)
    else:
        y = a + b / x ** c
    return y


def p0_simple_5reciprocal(xs, ys):
    c = 0.5
    b = (1/xs[-1]**c - 1/xs[1]**c) / (ys[-1] - ys[1])
    a = ys[1] - b / xs[1]**c
    return [a, b]


def measure(function, xs, ys, popt, weights):
    """
    measure the quality of the fit
    """
    m = 0
    n = 0
    for x in xs:
        if len(popt) == 2:
            m += abs(ys[n] - function(x, popt[0], popt[1]))  # * weights[n]
        elif len(popt) == 3:
            m += abs(ys[n] - function(x, popt[0], popt[1], popt[2]))   # * weights[n]
        else:
            raise NotImplementedError
        n += 1
    return m


def get_weigts(xs, ys, mode=1):
    if mode == 1:
        ds = get_derivatives(xs, ys, fd=True)
        import numpy as np
        mind = np.inf
        for d in ds:
            mind = min(abs(d), mind)
        weights = []
        for d in ds:
            weights.append(abs((mind / d)))
    else:
        weights = [1] * len(xs)
    return weights


def multy_curve_fit(xs, ys, verbose):
    """
    fit multiple functions to the x, y data, return the best fit
    """
    #functions = {exponential: p0_exponential, reciprocal: p0_reciprocal, single_reciprocal: p0_single_reciprocal}
    functions = {
        exponential: p0_exponential,
        reciprocal: p0_reciprocal,
        single_reciprocal: p0_single_reciprocal,
        simple_reciprocal: p0_simple_reciprocal,
        simple_2reciprocal: p0_simple_2reciprocal,
        simple_4reciprocal: p0_simple_4reciprocal,
        simple_5reciprocal: p0_simple_5reciprocal
    }
    import numpy as np
    from scipy.optimize import curve_fit
    fit_results = {}
    best = ['', np.inf]
    for function in functions:
        try:
            weights = get_weigts(xs, ys)
            popt, pcov = curve_fit(function, xs, ys, functions[function](xs, ys), maxfev=8000, sigma=weights)
            m = measure(function, xs, ys, popt, weights)
            perr = max(np.sqrt(np.diag(pcov)))
            #print 'pcov:\n', pcov
            #print 'diag:\n', np.sqrt(np.diag(pcov))
            #print 'function:\n', function, perr, m
            fit_results.update({function: {'measure': m, 'perr': perr, 'popt': popt, 'pcov': pcov}})
            for f in fit_results:
                if fit_results[f]['measure'] <= best[1]:
                    best = f, fit_results[f]['measure']
            print function, m
        except RuntimeError:
            if True:
                print 'no fit found for ', function

    return fit_results[best[0]]['popt'], fit_results[best[0]]['pcov'], best


def print_plot_line(function, popt, xs, ys, name, extra=''):
    """
    print the gnuplot command line to plot the x, y data with the fitted function using the popt parameters
    """
    idp = id_generator()
    f = open('convdat.'+str(idp), mode='w')
    for n in range(0, len(ys), 1):
        f.write(str(xs[n]) + ' ' + str(ys[n]) + '\n')
        f.write('\n')
    f.close()
    line = ''
    if function is exponential:
        line = "plot %s + %s * %s ** -x, 'convdat.%s' pointsize 4 lt 0, %s" % \
               (popt[0], popt[1], min(max(1.00001, popt[2]), 1.2), idp, popt[0])
    elif function is reciprocal:
        line = "plot %s + %s / x**%s, 'convdat.%s' pointsize 4 lt 0, %s" % \
               (popt[0], popt[1], min(max(1, popt[2]), 5), idp, popt[0])
    elif function is single_reciprocal:
        line = "plot %s + %s / (x - %s), 'convdat.%s' pointsize 4 lt 0, %s" % \
               (popt[0], popt[1], popt[2], idp, popt[0])
    elif function is simple_reciprocal:
        line = "plot %s + %s / x, 'convdat.%s' pointsize 4 lt 0, %s" % \
               (popt[0], popt[1], idp, popt[0])
    elif function is simple_2reciprocal:
        line = "plot %s + %s / x**2, 'convdat.%s' pointsize 4 lt 0, %s" % \
               (popt[0], popt[1], idp, popt[0])
    elif function is simple_4reciprocal:
        line = "plot %s + %s / x**4, 'convdat.%s' pointsize 4 lt 0, %s" % \
               (popt[0], popt[1], idp, popt[0])
    elif function is simple_5reciprocal:
        line = "plot %s + %s / x**0.5, 'convdat.%s' pointsize 4 lt 0, %s" % \
               (popt[0], popt[1], idp, popt[0])
    else:
        print function, ' no plot '
    f = open('plot-fits', mode='a')
    f.write('pause -1 \n')
    f.write('set title "' + name + ' - ' + extra + '"\n')
    f.write("set output '" + name + '-' + idp + ".gif'" + '\n')
    f.write(line + '\n')
    f.close()


def test_conv(xs, ys, name, tol=0.0001, extra='', verbose=False):
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
                #popt, pcov = curve_fit(exponential, xs, ys, p0_exponential(xs, ys), maxfev=8000)
                #perr = np.sqrt(np.diag(pcov))
                #print perr
                popt, pcov, func = multy_curve_fit(xs, ys, verbose)
                #print popt
                #print pcov
                #print func
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
                print_plot_line(func[0], popt, xs, ys, name, extra=extra)
              # print 'plot ', popt[0], ' + ', popt[1], "/x**", popt[2], ', "'+name+'.convdat"'
              #  print 'plot ', popt[0], ' + ', popt[1], "/x", popt[2], '/x**2, "'+name+'.convdat"'
              #  id = id_generator()
              #  print 'plot ', popt[0], ' + ', popt[1], "* ", popt[2], " ** -x," "'"+'convdat.'+id+"'"

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

            if verbose:
                print test

            if test < abs(tol):
                if verbose:
                    print 'converged'
                conv = True
                if xs[n] < x_value:
                    x_value = xs[n]
                    y_value = ys[n]
                    n_value = n
            else:
                if verbose:
                    print 'not converged'
                conv = False
                x_value = float('inf')
        if n_value is None:
            return [conv, x_value, y_value, n_value, popt[0], None]
        else:
            return [conv, x_value, y_value, n_value, popt[0], ds[n_value]]
    else:
        return [conv, x_value, y_value, n_value, popt[0], None]


