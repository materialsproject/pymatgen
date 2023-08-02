---
layout: default
title: pymatgen.analysis.eos.md
nav_exclude: true
---

# pymatgen.analysis.eos module

This module implements various equation of states.

Note: Most of the code were initially adapted from ASE and deltafactor by
@gmatteo but has since undergone major refactoring.


### _class_ pymatgen.analysis.eos.Birch(volumes, energies)
Bases: `EOSBase`

Birch EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



### _class_ pymatgen.analysis.eos.BirchMurnaghan(volumes, energies)
Bases: `EOSBase`

BirchMurnaghan EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



### _class_ pymatgen.analysis.eos.DeltaFactor(volumes, energies)
Bases: `PolynomialEOS`

Fitting a polynomial EOS using delta factor.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### fit(order=3)
Overridden since this eos works with volume\*\*(2/3) instead of volume.


### _class_ pymatgen.analysis.eos.EOS(eos_name='murnaghan')
Bases: `object`

Convenient wrapper. Retained in its original state to ensure backward
compatibility.

Fit equation of state for bulk systems.

The following equations are supported:

```default
murnaghan: PRB 28, 5480 (1983)

birch: Intermetallic compounds: Principles and Practice, Vol I:
    Principles. pages 195-210

birch_murnaghan: PRB 70, 224107

pourier_tarantola: PRB 70, 224107

vinet: PRB 70, 224107

deltafactor

numerical_eos: 10.1103/PhysRevB.90.174107.
```

Usage:

```default
eos = EOS(eos_name='murnaghan')
eos_fit = eos.fit(volumes, energies)
eos_fit.plot()
```


* **Parameters**

    **eos_name** (*str*) – Type of EOS to fit.



#### MODELS(_ = {'birch': <class 'pymatgen.analysis.eos.Birch'>, 'birch_murnaghan': <class 'pymatgen.analysis.eos.BirchMurnaghan'>, 'deltafactor': <class 'pymatgen.analysis.eos.DeltaFactor'>, 'murnaghan': <class 'pymatgen.analysis.eos.Murnaghan'>, 'numerical_eos': <class 'pymatgen.analysis.eos.NumericalEOS'>, 'pourier_tarantola': <class 'pymatgen.analysis.eos.PourierTarantola'>, 'vinet': <class 'pymatgen.analysis.eos.Vinet'>_ )

#### fit(volumes, energies)
Fit energies as function of volumes.


* **Parameters**


    * **volumes** (*list/np.array*) –


    * **energies** (*list/np.array*) –



* **Returns**

    EOSBase object



* **Return type**

    EOSBase



### _class_ pymatgen.analysis.eos.EOSBase(volumes, energies)
Bases: `object`

Abstract class that must be subclassed by all equation of state
implementations.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _property_ b0()
Returns the bulk modulus.
Note: the units for the bulk modulus: unit of energy/unit of volume^3.


#### _property_ b0_GPa()
Returns the bulk modulus in GPa.
Note: This assumes that the energy and volumes are in eV and Ang^3

> respectively.


#### _property_ b1()
Returns the derivative of bulk modulus wrt pressure(dimensionless).


#### _property_ e0()
Returns the min energy.


#### fit()
Do the fitting. Does least square fitting. If you want to use custom
fitting, must override this.


#### func(volume)
The equation of state function with the parameters other than volume set
to the ones obtained from fitting.


* **Parameters**

    **volume** (*list/numpy.array*) –



* **Returns**

    numpy.array



#### plot(width=8, height=None, plt=None, dpi=None, \*\*kwargs)
Plot the equation of state.


* **Parameters**


    * **width** (*float*) – Width of plot in inches. Defaults to 8in.


    * **height** (*float*) – Height of plot in inches. Defaults to width \*
    golden ratio.


    * **plt** (*matplotlib.pyplot*) – If plt is supplied, changes will be made
    to an existing plot. Otherwise, a new plot will be created.


    * **dpi** –


    * **kwargs** (*dict*) – additional args fed to pyplot.plot.
    supported keys: style, color, text, label



* **Returns**

    Matplotlib plot object.



#### plot_ax(ax=None, fontsize=12, \*\*kwargs)
Plot the equation of state on axis ax.


* **Parameters**


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **fontsize** – Legend fontsize.


    * **color** (*str*) – plot color.


    * **label** (*str*) – Plot label


    * **text** (*str*) – Legend text (options)



* **Returns**

    Matplotlib figure object.


Keyword arguments controlling the display of the figure:

| kwargs

 | Meaning

 |
| ------------ | ------------------------------------------------------------------------------------------------- |  |  |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### _property_ results()
Returns a summary dict.


* **Returns**

    dict



#### _property_ v0()
Returns the minimum or the reference volume in Ang^3.


### _exception_ pymatgen.analysis.eos.EOSError()
Bases: `Exception`

Error class for EOS fitting.


### _class_ pymatgen.analysis.eos.Murnaghan(volumes, energies)
Bases: `EOSBase`

Murnaghan EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



### _class_ pymatgen.analysis.eos.NumericalEOS(volumes, energies)
Bases: `PolynomialEOS`

A numerical EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### fit(min_ndata_factor=3, max_poly_order_factor=5, min_poly_order=2)
Fit the input data to the ‘numerical eos’, the equation of state employed
in the quasiharmonic Debye model described in the paper:
10.1103/PhysRevB.90.174107.

credits: Cormac Toher


* **Parameters**


    * **min_ndata_factor** (*int*) – parameter that controls the minimum number
    of data points that will be used for fitting.
    minimum number of data points =

    > total data points-2\*min_ndata_factor



    * **max_poly_order_factor** (*int*) – parameter that limits the max order
    of the polynomial used for fitting.
    max_poly_order = number of data points used for fitting -

    > max_poly_order_factor



    * **min_poly_order** (*int*) – minimum order of the polynomial to be
    considered for fitting.



### _class_ pymatgen.analysis.eos.PolynomialEOS(volumes, energies)
Bases: `EOSBase`

Derives from EOSBase. Polynomial based equations of states must subclass
this.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### fit(order)
Do polynomial fitting and set the parameters. Uses numpy polyfit.


* **Parameters**

    **order** (*int*) – order of the fit polynomial



### _class_ pymatgen.analysis.eos.PourierTarantola(volumes, energies)
Bases: `EOSBase`

PourierTarantola EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



### _class_ pymatgen.analysis.eos.Vinet(volumes, energies)
Bases: `EOSBase`

Vinet EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.