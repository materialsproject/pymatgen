---
layout: default
title: pymatgen.core.spectrum.md
nav_exclude: true
---

# pymatgen.core.spectrum module

This module defines classes to represent any type of spectrum, essentially any
x y value pairs.


### _class_ pymatgen.core.spectrum.Spectrum(x: ArrayLike, y: ArrayLike, \*args, \*\*kwargs)
Bases: `MSONable`

Base class for any type of xas, essentially just x, y values. Examples
include XRD patterns, XANES, EXAFS, NMR, DOS, etc.

Implements basic tools like application of smearing, normalization, addition
multiplication, etc.

Subclasses should extend this object and ensure that super is called with
ALL args and kwargs. That ensures subsequent things like add and mult work
properly.


* **Parameters**


    * **x** (*ndarray*) – A ndarray of N values.


    * **y** (*ndarray*) – A ndarray of N x k values. The first dimension must be
    the same as that of x. Each of the k values are interpreted as separate.


    * **\*args** – All subclasses should provide args other than x and y
    when calling super, e.g., super().__init__(
    x, y, arg1, arg2, kwarg1=val1, ..). This guarantees the +, -,

    ```
    *
    ```

    ,
    etc. operators work properly.


    * **\*\*kwargs** – Same as that for

    ```
    *
    ```

    args.



#### XLABEL(_ = 'x_ )

#### YLABEL(_ = 'y_ )

#### copy()

* **Returns**

    Copy of Spectrum object.



#### get_interpolated_value(x: float)
Returns an interpolated y value for a particular x value.


* **Parameters**

    **x** – x value to return the y value for



* **Returns**

    Value of y at x



#### normalize(mode: Literal['max', 'sum'] = 'max', value: float = 1.0)
Normalize the spectrum with respect to the sum of intensity.


* **Parameters**


    * **mode** (*"max"** | **"sum"*) – Normalization mode. “max” sets the max y value to value,
    e.g., in XRD patterns. “sum” sets the sum of y to a value, i.e., like a
    probability density.


    * **value** (*float*) – Value to normalize to. Defaults to 1.



#### smear(sigma: float = 0.0, func: str | Callable = 'gaussian')
Apply Gaussian/Lorentzian smearing to spectrum y value.


* **Parameters**


    * **sigma** – Std dev for Gaussian smear function


    * **func** – “gaussian” or “lorentzian” or a callable. If this is a callable, the sigma value is ignored. The
    callable should only take a single argument (a numpy array) and return a set of weights.



### pymatgen.core.spectrum.lorentzian(x, x_0: float = 0, sigma: float = 1.0)

* **Parameters**


    * **x** – x values


    * **x_0** – Center


    * **sigma** – FWHM



* **Returns**

    Value of lorentzian at x.