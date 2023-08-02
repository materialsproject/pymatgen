---
layout: default
title: pymatgen.phonon.ir_spectra.md
nav_exclude: true
---

# pymatgen.phonon.ir_spectra module

This module provides classes to handle the calculation of the IR spectra
This implementation is adapted from Abipy
[https://github.com/abinit/abipy](https://github.com/abinit/abipy)
where it was originally done by Guido Petretto and Matteo Giantomassi.


### _class_ pymatgen.phonon.ir_spectra.IRDielectricTensor(oscillator_strength, ph_freqs_gamma, epsilon_infinity, structure)
Bases: `MSONable`

Class to handle the Ionic Dielectric Tensor
The implementation is adapted from Abipy
See the definitions Eq.(53-54) in

```
:cite:`Gonze1997`
```

 PRB55, 10355 (1997).


* **Parameters**


    * **oscillator_strength** – IR oscillator strengths as defined in Eq. 54 in


    ```
    :cite:`Gonze1997`
    ```

     PRB55, 10355 (1997).


    * **ph_freqs_gamma** – Phonon frequencies at the Gamma point


    * **epsilon_infinity** – electronic susceptibility as defined in Eq. 29.


    * **structure** – A Structure object corresponding to the structure used for the calculation.



#### as_dict()
JSON-serializable dict representation of IRDielectricTensor.


#### _classmethod_ from_dict(d)
Returns IRDielectricTensor from dict representation.


#### get_ir_spectra(broad=5e-05, emin=0, emax=None, divs=500)
The IR spectra is obtained for the different directions.


* **Parameters**


    * **broad** – a list of broadenings or a single broadening for the phonon peaks


    * **emin** (*float*) – minimum energy in which to obtain the spectra. Defaults to 0.


    * **emax** (*float*) – maximum energy in which to obtain the spectra. Defaults to None.


    * **divs** – number of frequency samples between emin and emax



* **Returns**

    divs array with the frequencies at which the

        dielectric tensor is calculated

    dielectric_tensor: divsx3x3 numpy array with the dielectric tensor

        for the range of frequencies




* **Return type**

    frequencies



#### get_plotter(components=('xx',), reim='reim', broad=5e-05, emin=0, emax=None, divs=500, \*\*kwargs)
Return an instance of the Spectrum plotter containing the different requested components.


* **Parameters**


    * **components** – A list with the components of the dielectric tensor to plot.
    Can be either two indexes or a string like ‘xx’ to plot the (0,0) component


    * **reim** – If ‘re’ (im) is present in the string plots the real (imaginary) part of the dielectric tensor


    * **broad** (*float*) – a list of broadenings or a single broadening for the phonon peaks. Defaults to 0.00005.


    * **emin** (*float*) – minimum energy in which to obtain the spectra. Defaults to 0.


    * **emax** (*float*) – maximum energy in which to obtain the spectra. Defaults to None.


    * **divs** – number of frequency samples between emin and emax


    * **\*\*kwargs** – Passed to IRDielectricTensor.get_spectrum()



#### get_spectrum(component, reim, broad=5e-05, emin=0, emax=None, divs=500, label=None)
component: either two indexes or a string like ‘xx’ to plot the (0,0) component
reim: only “re” or “im”
broad: a list of broadenings or a single broadening for the phonon peaks.


#### _property_ max_phfreq()
Maximum phonon frequency.


#### _property_ nph_freqs()
Number of phonon frequencies.


#### plot(components=('xx',), reim='reim', show_phonon_frequencies=True, xlim=None, ylim=None, \*\*kwargs)
Helper function to generate the Spectrum plotter and directly plot the results.


* **Parameters**


    * **components** – A list with the components of the dielectric tensor to plot.
    Can be either two indexes or a string like ‘xx’ to plot the (0,0) component


    * **reim** – If ‘re’ (im) is present in the string plots the real (imaginary) part of the dielectric tensor


    * **show_phonon_frequencies** – plot a dot where the phonon frequencies are to help identify IR inactive modes


    * **xlim** – x-limits of the plot. Defaults to None for automatic determination.


    * **ylim** – y-limits of the plot. Defaults to None for automatic determination.


    * **kwargs** – keyword arguments passed to the plotter


Keyword arguments controlling the display of the figure:

| kwargs

 | Meaning

 |
| ------------ | ------------------------------------------------------------------------------------------------- |  |  |  |  |  |  |  |  |  |  |
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

#### write_json(filename)
Save a json file with this data.