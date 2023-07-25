"""
Functions for calculating the convergence of an x, y data set.

Main API:

test_conv(xs, ys, name, tol)

tries to fit multiple functions to the x, y data

calculates which function fits best
for tol < 0
returns the x value for which y is converged within tol of the asymptotic value
for tol > 0
returns the x_value for which dy(x)/dx < tol for all x >= x_value, conv is true is such a x_value exists
for the best fit a gnuplot line is printed plotting the data, the function and the asymptotic value
"""

from __future__ import annotations

import random
import string

import numpy as np
from monty.dev import deprecated

deprecated(message="pymatgen.util.convergence is deprecated and will be removed in a late 2023 release.")

__author__ = "Michiel van Setten"
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "June 2014"


def id_generator(size: int = 8, chars: str = string.ascii_uppercase + string.digits) -> str:
    """
    Generate a random string of specified size and characters.

    Args:
        size (int): The length of the generated string.
        chars (str): The characters to use for generating the string.

    Returns:
        str: The generated random string.
    """
    return "".join(random.choice(chars) for _ in range(size))


class SplineInputError(Exception):
    """Error for Spline input."""

    def __init__(self, msg):
        """
        Args:
            msg (str): Message.
        """
        self.msg = msg


def get_derivatives(xs, ys, fd=False):
    """
    return the derivatives of y(x) at the points x
    if scipy is available a spline is generated to calculate the derivatives
    if scipy is not available the left and right slopes are calculated, if both exist the average is returned
    putting fd to zero always returns the finite difference slopes.
    """
    try:
        if fd:
            raise SplineInputError("no spline wanted")
        if len(xs) < 4:
            er = SplineInputError("too few data points")
            raise er
        from scipy.interpolate import UnivariateSpline

        spline = UnivariateSpline(xs, ys)
        d = spline.derivative(1)(xs)
    except (ImportError, SplineInputError):
        d = []
        m, left, right = 0, 0, 0
        for n in range(0, len(xs), 1):
            try:
                left = (ys[n] - ys[n - 1]) / (xs[n] - xs[n - 1])
                m += 1
            except IndexError:
                pass
            try:
                right = (ys[n + 1] - ys[n]) / (xs[n + 1] - xs[n])
                m += 1
            except IndexError:
                pass
            d.append(left + right / m)
    return d


"""
functions used in the fitting procedure, with initial guesses
"""


def print_and_raise_error(xs, ys, name):
    """
    Print error message and raise a RuntimeError.

    Args:
        xs (list): List of x values.
        ys (list): List of y values.
        name (str): Name of the function where the error occurred.
    """
    print("Index error in", name)
    print("ys: ", ys)
    print("xs: ", xs)
    raise RuntimeError


def reciprocal(x, a, b, n):
    """Reciprocal function to the power n to fit convergence data."""
    if n < 1:
        n = 1
    elif n > 5:
        n = 5
    if isinstance(x, list):
        return np.array([a + b / x_v**n for x_v in x])
    return a + b / x**n


def p0_reciprocal(xs, ys):
    """Predictor for first guess for reciprocal."""
    a0 = ys[len(ys) - 1]
    b0 = ys[0] * xs[0] - a0 * xs[0]
    return [a0, b0, 1]


def exponential(x, a, b, n):
    """Exponential function base n to fit convergence data."""
    if n < 1.000001:
        n = 1.000001
    elif n > 1.2:
        n = 1.2
    if b < -10:
        b = -10
    elif b > 10:
        b = 10
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b * n**-x_v)
        y = np.array(y_l)
    else:
        y = a + b * n**-x
    return y


def p0_exponential(xs, ys):
    """
    Calculate the initial guess parameters for the exponential function.

    Args:
        xs (list): List of x values.
        ys (list): List of y values.

    Returns:
        list: List of initial guess parameters [a, b, n].
    """
    n0 = 1.005
    b0 = (n0 ** -xs[-1] - n0 ** -xs[1]) / (ys[-1] - ys[1])
    a0 = ys[1] - b0 * n0 ** -xs[1]
    # a0 = ys[-1]
    # b0 = (ys[0] - a0) / n0 ** xs[0]
    return [a0, b0, n0]


def single_reciprocal(x, a, b, c):
    """Reciprocal function to fit convergence data."""
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / (x_v - c))
        y = np.array(y_l)
    else:
        y = a + b / (x - c)
    return y


def p0_single_reciprocal(xs, ys):
    """
    Calculate the initial guess parameters for the single reciprocal function.

    Args:
        xs (list): List of x values.
        ys (list): List of y values.

    Returns:
        list: List of initial guess parameters [a, b, c].
    """
    c = 1
    b = (1 / (xs[-1] - c) - 1 / (xs[1] - c)) / (ys[-1] - ys[1])
    a = ys[1] - b / (xs[1] - c)
    return [a, b, c]


def simple_reciprocal(x, a, b):
    """Reciprocal function to fit convergence data."""
    y_l = [a + b / x_v for x_v in x]
    return np.array(y_l) if isinstance(x, list) else a + b / x


def p0_simple_reciprocal(xs, ys):
    """
    Calculate the initial guess parameters for the simple reciprocal function.

    Args:
        xs (list): List of x values.
        ys (list): List of y values.

    Returns:
        list: List of initial guess parameters [a, b].
    """
    b = (ys[-1] - ys[-2]) / (1 / (xs[-1]) - 1 / (xs[-2]))
    a = ys[-2] - b / (xs[-2])
    return [a, b]


def simple_2reciprocal(x, a, b):
    """Reciprocal function to fit convergence data."""
    c = 2
    if isinstance(x, list):
        y_l = []
        for x_v in x:
            y_l.append(a + b / x_v**c)
        y = np.array(y_l)
    else:
        y = a + b / x**c
    return y


def p0_simple_2reciprocal(xs, ys):
    """
    Calculate the initial guess parameters for the simple reciprocal function with a power of 2.

    Args:
        xs (list): List of x values.
        ys (list): List of y values.

    Returns:
        list: List of initial guess parameters [a, b].
    """
    c = 2
    b = (ys[-1] - ys[1]) / (1 / xs[-1] ** c - 1 / xs[1] ** c)
    a = ys[1] - b / xs[1] ** c
    return [a, b]


def simple_4reciprocal(x, a, b):
    """Reciprocal function to fit convergence data."""
    c = 4
    return [a + b / x_v**c for x_v in x] if isinstance(x, list) else a + b / x**c


def p0_simple_4reciprocal(xs, ys):
    """
    Calculate the initial guess parameters for the simple reciprocal function with a power of 4.

    Args:
        xs (list): List of x values.
        ys (list): List of y values.

    Returns:
        list: The initial guess parameters [a, b].
    """
    c = 4
    b = (ys[-1] - ys[1]) / (1 / xs[-1] ** c - 1 / xs[1] ** c)
    a = ys[1] - b / xs[1] ** c
    return [a, b]


def simple_5reciprocal(x, a, b):
    """Reciprocal function to fit convergence data."""
    c = 0.5
    return [a + b / x_v**c for x_v in x] if isinstance(x, list) else a + b / x**c


def p0_simple_5reciprocal(xs, ys):
    """
    Calculate the initial guess parameters for the simple reciprocal function with a power of 0.5.

    Args:
        xs (list): List of x values.
        ys (list): List of y values.

    Returns:
        list: List of parameters [a, b].
    """
    c = 0.5
    b = (ys[-1] - ys[1]) / (1 / xs[-1] ** c - 1 / xs[1] ** c)
    a = ys[1] - b / xs[1] ** c
    return [a, b]


def extrapolate_simple_reciprocal(xs, ys):
    """
    Extrapolate simple reciprocal function to fit convergence data.

    Args:
        xs: List of x values.
        ys: List of y values.

    Returns:
        List of parameters [a, b].
    """
    b = (ys[-2] - ys[-1]) / (1 / (xs[-2]) - 1 / (xs[-1]))
    a = ys[-1] - b / (xs[-1])
    return [a, b]


def extrapolate_reciprocal(xs, ys, n, noise):
    """Return the parameters such that a + b / x^n hits the last two data points."""
    if len(xs) > 4 and noise:
        y1 = (ys[-3] + ys[-4]) / 2
        y2 = (ys[-1] + ys[-2]) / 2
        x1 = (xs[-3] + xs[-4]) / 2
        x2 = (xs[-1] + xs[-2]) / 2
        try:
            b = (y1 - y2) / (1 / x1**n - 1 / x2**n)
            a = y2 - b / x2**n
        except IndexError:
            print_and_raise_error(xs, ys, "extrapolate_reciprocal")
    else:
        try:
            b = (ys[-2] - ys[-1]) / (1 / (xs[-2]) ** n - 1 / (xs[-1]) ** n)
            a = ys[-1] - b / (xs[-1]) ** n
        except IndexError:
            print_and_raise_error(xs, ys, "extrapolate_reciprocal")
    return [a, b, n]


def measure(function, xs, ys, popt, weights):
    """Measure the quality of a fit."""
    m = 0
    n = 0
    for x in xs:
        try:
            if len(popt) == 2:
                m += (ys[n] - function(x, popt[0], popt[1])) ** 2 * weights[n]
            elif len(popt) == 3:
                m += (ys[n] - function(x, popt[0], popt[1], popt[2])) ** 2 * weights[n]
            else:
                raise NotImplementedError
            n += 1
        except IndexError as exc:
            raise IndexError(f"y does not exist for {x = }, this should not happen") from exc

    return m


def get_weights(xs, ys, mode=2):
    """
    Args:
        xs (list): List of x values.
        ys (list): List of y values.
        mode (int): Mode for calculating weights.

    Returns:
        list: List of weights.
    """
    ds = get_derivatives(xs, ys, fd=True)
    if mode == 1:
        mind = np.inf
        for d in ds:
            mind = min(abs(d), mind)
        weights = [abs(mind / d) for d in ds]
    elif mode == 2:
        x_max = max(xs) ** 2
        weights = [x**2 / x_max for x in xs]
    else:
        weights = [1] * len(xs)
    return weights


def multi_curve_fit(xs, ys, verbose):
    """Fit multiple functions to the x, y data, return the best fit."""
    # functions = {exponential: p0_exponential, reciprocal: p0_reciprocal, single_reciprocal: p0_single_reciprocal}
    functions = {
        exponential: p0_exponential,
        reciprocal: p0_reciprocal,
        # single_reciprocal: p0_single_reciprocal,
        simple_reciprocal: p0_simple_reciprocal,
        simple_2reciprocal: p0_simple_2reciprocal,
        simple_4reciprocal: p0_simple_4reciprocal,
        simple_5reciprocal: p0_simple_5reciprocal,
    }
    from scipy.optimize import curve_fit

    fit_results = {}
    best = ["", np.inf]
    for k, v in functions.items():
        try:
            weights = get_weights(xs, ys)
            popt, pcov = curve_fit(
                k,
                xs,
                ys,
                v(xs, ys),
                maxfev=8000,
                sigma=weights,
            )
            pcov = []
            m = measure(k, xs, ys, popt, weights)
            fit_results.update({k: {"measure": m, "popt": popt, "pcov": pcov}})
            for f, v in fit_results.items():
                if v["measure"] <= best[1]:
                    best = f, v["measure"]
            if verbose:
                print(str(k), m)
        except RuntimeError:
            print("no fit found for ", k)

    return fit_results[best[0]]["popt"], fit_results[best[0]]["pcov"], best


def multi_reciprocal_extra(xs, ys, noise=False):
    """
    Calculates for a series of powers ns the parameters for which the last two points are at the curve.
    With these parameters measure how well the other data points fit.
    return the best fit.
    """
    ns = np.linspace(0.5, 6.0, num=56)
    best = ["", np.inf]
    fit_results = {}
    weights = get_weights(xs, ys)
    for n in ns:
        popt = extrapolate_reciprocal(xs, ys, n, noise)
        m = measure(reciprocal, xs, ys, popt, weights)
        pcov = []
        fit_results.update({n: {"measure": m, "popt": popt, "pcov": pcov}})
    for n, v in fit_results.items():
        if v["measure"] <= best[1]:
            best = reciprocal, v["measure"], n
    return fit_results[best[2]]["popt"], fit_results[best[2]]["pcov"], best


def print_plot_line(function, popt, xs, ys, name, tol: float = 0.05, extra=""):
    """Print the gnuplot command line to plot the x, y data with the fitted function using the popt parameters."""
    idp = id_generator()
    with open(f"convdat.{idp}", mode="w") as f:
        for n in range(0, len(ys), 1):
            f.write(f"{xs[n]} {ys[n]}\n")
    tol = abs(tol)
    line = f"plot 'convdat.{idp}' pointsize 4 lt 0, "
    line += f"{popt[0]} lt 3, {popt[0] - tol} lt 4, {popt[0] + tol} lt 4, "
    if function is exponential:
        line += f"{popt[0]} + {popt[1]} * {min(max(1.00001, popt[2]), 1.2)} ** -x"
    elif function is reciprocal:
        line += f"{popt[0]} + {popt[1]} / x**{min(max(0.5, popt[2]), 6)}"
    elif function is single_reciprocal:
        line += f"{popt[0]} + {popt[1]} / (x - {popt[2]})"
    elif function is simple_reciprocal:
        line += f"{popt[0]} + {popt[1]} / x"
    elif function is simple_2reciprocal:
        line += f"{popt[0]} + {popt[1]} / x**2"
    elif function is simple_4reciprocal:
        line += f"{popt[0]} + {popt[1]} / x**4"
    elif function is simple_5reciprocal:
        line += f"{popt[0]} + {popt[1]} / x**0.5"
    else:
        print(function, " no plot ")

    with open("plot-fits", mode="a") as f:
        f.write('set title "' + name + " - " + extra + '"\n')
        f.write("set output '" + name + "-" + idp + ".gif'\n")
        f.write(f"set yrange [{popt[0] - 5 * tol}:{popt[0] + 5 * tol}]\n")
        f.write(line + "\n")
        f.write("pause -1 \n")


def determine_convergence(xs, ys, name, tol: float = 0.0001, extra="", verbose=False, mode="extra", plots=True):
    """Test it and at which x_value dy(x)/dx < tol for all x >= x_value, conv is true is such a x_value exists."""
    if len(xs) != len(ys):
        raise RuntimeError("the range of x and y are not equal")
    conv = False
    x_value = float("inf")
    y_value = n_value = None
    popt = [None, None, None]
    if len(xs) > 2:
        ds = get_derivatives(xs[0 : len(ys)], ys)
        try:
            if None not in ys:
                if mode == "fit":
                    popt, pcov, func = multi_curve_fit(xs, ys, verbose)
                elif mode == "extra":
                    res = multi_reciprocal_extra(xs, ys)
                    if res is not None:
                        popt, pcov, func = multi_reciprocal_extra(xs, ys)
                    else:
                        print(xs, ys)
                        popt, pcov = None, None  # type: ignore[assignment]
                elif mode == "extra_noise":
                    popt, pcov, func = multi_reciprocal_extra(xs, ys, noise=True)
                else:
                    raise NotImplementedError("unknown mode for test conv")
                if func[1] > abs(tol):
                    print(f"warning function {func[0]} as the best fit but not a good fit: {func[1]}")
                # TODO print this to file via a method in helper, as dict
                if plots:
                    with open(f"{name}.fitdat", mode="a") as f:
                        f.write("{")
                        f.write(f'"popt": {popt}, ')
                        f.write(f'"pcov": {pcov}, ')
                        f.write('"data": [')
                        for n in range(0, len(ys), 1):
                            f.write(f"[{xs[n]} {ys[n]}]")
                        f.write("]}\n")

                    print_plot_line(func[0], popt, xs, ys, name, tol=tol, extra=extra)

        except ImportError:
            popt, pcov = None, None  # type: ignore[assignment]
        for n in range(0, len(ds), 1):
            if verbose:
                print(n, ys[n])
                print(ys)
            test = (abs(popt[0] - ys[n]) if popt[0] is not None else float("inf")) if tol < 0 else abs(ds[n])
            if verbose:
                print(test)
            if test < abs(tol):
                if verbose:
                    print("converged")
                conv = True
                if xs[n] < x_value:
                    x_value = xs[n]
                    y_value = ys[n]
                    n_value = n
            else:
                if verbose:
                    print("not converged")
                conv = False
                x_value = float("inf")
        if n_value is None:
            return [conv, x_value, y_value, n_value, popt[0], None]

        return [conv, x_value, y_value, n_value, popt[0], ds[n_value]]

    return [conv, x_value, y_value, n_value, popt[0], None]
