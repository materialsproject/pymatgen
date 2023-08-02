---
layout: default
title: pymatgen.analysis.xps.md
nav_exclude: true
---

# pymatgen.analysis.xps module

This is a module for XPS analysis. It is modelled after the Galore package ([https://github.com/SMTG-UCL/galore](https://github.com/SMTG-UCL/galore)), but
with some modifications for easier analysis from pymatgen itself. Please cite the following original work if you use
this:

```default
Adam J. Jackson, Alex M. Ganose, Anna Regoutz, Russell G. Egdell, David O. Scanlon (2018). Galore: Broadening and
weighting for simulation of photoelectron spectroscopy. Journal of Open Source Software, 3(26), 773,
doi: 10.21105/joss.007733
```

You may wish to look at the optional dependency galore for more functionality such as plotting and other cross-sections.
Note that the atomic_subshell_photoionization_cross_sections.csv has been reparsed from the original compilation:

```default
Yeh, J. J.; Lindau, I. Atomic Subshell Photoionization Cross Sections and Asymmetry Parameters: 1 ⩽ Z ⩽ 103.
Atomic Data and Nuclear Data Tables 1985, 32 (1), 1-155. https://doi.org/10.1016/0092-640X(85)90016-6.
```

This version contains all detailed information for all orbitals.


### _class_ pymatgen.analysis.xps.XPS(x: ArrayLike, y: ArrayLike, \*args, \*\*kwargs)
Bases: [`Spectrum`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum)

Class representing an X-ray photoelectron spectra.


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



#### XLABEL(_ = 'Binding Energy (eV)_ )

#### YLABEL(_ = 'Intensity_ )

#### _classmethod_ from_dos(dos: [CompleteDos](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos))

* **Parameters**


    * **dos** – CompleteDos object with project element-orbital DOS. Can be obtained from Vasprun.get_complete_dos.


    * **sigma** – Smearing for Gaussian.



* **Returns**

    XPS