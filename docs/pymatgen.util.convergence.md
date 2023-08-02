---
layout: default
title: pymatgen.util.convergence.md
nav_exclude: true
---

# pymatgen.util.convergence module

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


### _exception_ pymatgen.util.convergence.SplineInputError(msg)
Bases: `Exception`

Error for Spline input.


* **Parameters**

    **msg** (*str*) – Message.



### pymatgen.util.convergence.determine_convergence(xs, ys, name, tol: float = 0.0001, extra='', verbose=False, mode='extra', plots=True)
Test it and at which x_value dy(x)/dx < tol for all x >= x_value, conv is true is such a x_value exists.


### pymatgen.util.convergence.exponential(x, a, b, n)
Exponential function base n to fit convergence data.


### pymatgen.util.convergence.extrapolate_reciprocal(xs, ys, n, noise)
Return the parameters such that a + b / x^n hits the last two data points.


### pymatgen.util.convergence.extrapolate_simple_reciprocal(xs, ys)
Extrapolate simple reciprocal function to fit convergence data.


* **Parameters**


    * **xs** – List of x values.


    * **ys** – List of y values.



* **Returns**

    List of parameters [a, b].



### pymatgen.util.convergence.get_derivatives(xs, ys, fd=False)
return the derivatives of y(x) at the points x
if scipy is available a spline is generated to calculate the derivatives
if scipy is not available the left and right slopes are calculated, if both exist the average is returned
putting fd to zero always returns the finite difference slopes.


### pymatgen.util.convergence.get_weights(xs, ys, mode=2)

* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.


    * **mode** (*int*) – Mode for calculating weights.



* **Returns**

    List of weights.



* **Return type**

    list



### pymatgen.util.convergence.id_generator(size: int = 8, chars: str = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
Generate a random string of specified size and characters.


* **Parameters**


    * **size** (*int*) – The length of the generated string.


    * **chars** (*str*) – The characters to use for generating the string.



* **Returns**

    The generated random string.



* **Return type**

    str



### pymatgen.util.convergence.measure(function, xs, ys, popt, weights)
Measure the quality of a fit.


### pymatgen.util.convergence.multi_curve_fit(xs, ys, verbose)
Fit multiple functions to the x, y data, return the best fit.


### pymatgen.util.convergence.multi_reciprocal_extra(xs, ys, noise=False)
Calculates for a series of powers ns the parameters for which the last two points are at the curve.
With these parameters measure how well the other data points fit.
return the best fit.


### pymatgen.util.convergence.p0_exponential(xs, ys)
Calculate the initial guess parameters for the exponential function.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b, n].



* **Return type**

    list



### pymatgen.util.convergence.p0_reciprocal(xs, ys)
Predictor for first guess for reciprocal.


### pymatgen.util.convergence.p0_simple_2reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function with a power of 2.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b].



* **Return type**

    list



### pymatgen.util.convergence.p0_simple_4reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function with a power of 4.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    The initial guess parameters [a, b].



* **Return type**

    list



### pymatgen.util.convergence.p0_simple_5reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function with a power of 0.5.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of parameters [a, b].



* **Return type**

    list



### pymatgen.util.convergence.p0_simple_reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b].



* **Return type**

    list



### pymatgen.util.convergence.p0_single_reciprocal(xs, ys)
Calculate the initial guess parameters for the single reciprocal function.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b, c].



* **Return type**

    list



### pymatgen.util.convergence.print_and_raise_error(xs, ys, name)
Print error message and raise a RuntimeError.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.


    * **name** (*str*) – Name of the function where the error occurred.



### pymatgen.util.convergence.print_plot_line(function, popt, xs, ys, name, tol: float = 0.05, extra='')
Print the gnuplot command line to plot the x, y data with the fitted function using the popt parameters.


### pymatgen.util.convergence.reciprocal(x, a, b, n)
Reciprocal function to the power n to fit convergence data.


### pymatgen.util.convergence.simple_2reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### pymatgen.util.convergence.simple_4reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### pymatgen.util.convergence.simple_5reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### pymatgen.util.convergence.simple_reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### pymatgen.util.convergence.single_reciprocal(x, a, b, c)
Reciprocal function to fit convergence data.