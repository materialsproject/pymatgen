---
layout: default
title: pymatgen.io.abinit.pseudos.md
nav_exclude: true
---

# pymatgen.io.abinit.pseudos module

This module provides objects describing the basic parameters of the
pseudopotentials used in Abinit, and a parser to instantiate pseudopotential objects..


### _class_ pymatgen.io.abinit.pseudos.AbinitHeader()
Bases: `dict`

Dictionary whose keys can be also accessed as attributes.


### _class_ pymatgen.io.abinit.pseudos.AbinitPseudo(path, header)
Bases: `Pseudo`

An AbinitPseudo is a pseudopotential whose file contains an abinit header.


* **Parameters**


    * **path** – Filename.


    * **header** – `AbinitHeader` instance.



#### _property_ Z()
The atomic number of the atom.


#### _property_ Z_val()
Valence charge.


#### _property_ l_local()
Angular momentum used for the local part.


#### _property_ l_max()
Maximum angular momentum.


#### _property_ summary()
Summary line reported in the ABINIT header.


#### _property_ supports_soc()
True if the pseudo can be used in a calculation with spin-orbit coupling.
Base classes should provide a concrete implementation that computes this value.


### _class_ pymatgen.io.abinit.pseudos.Hint(ecut, pawecutdg=None)
Bases: `object`

Suggested value for the cutoff energy [Hartree units]
and the cutoff energy for the dense grid (only for PAW pseudos).


#### as_dict()
Return dictionary for MSONable protocol.


#### _classmethod_ from_dict(d)
Build instance from dictionary (MSONable protocol).


### _class_ pymatgen.io.abinit.pseudos.NcAbinitHeader(summary, \*\*kwargs)
Bases: `AbinitHeader`

The abinit header found in the NC pseudopotential files.


#### _static_ fhi_header(filename, ppdesc)
Parse the FHI abinit header. Example:

Troullier-Martins psp for element  Sc        Thu Oct 27 17:33:22 EDT 1994

    21.00000   3.00000    940714                zatom, zion, pspdat

        1    1    2    0      2001    .00000      pspcod,pspxc,lmax,lloc,mmax,r2well

1.80626423934776     .22824404341771    1.17378968127746   rchrg,fchrg,qchrg


#### _static_ gth_header(filename, ppdesc)
Parse the GTH abinit header. Example:

Goedecker-Teter-Hutter  Wed May  8 14:27:44 EDT 1996
1   1   960508                     zatom,zion,pspdat
2   1   0    0    2001    0.       pspcod,pspxc,lmax,lloc,mmax,r2well
0.2000000 -4.0663326  0.6778322 0 0     rloc, c1, c2, c3, c4
0 0 0                              rs, h1s, h2s
0 0                                rp, h1p

> 1.36 .2   0.6                    rcutoff, rloc


#### _static_ hgh_header(filename, ppdesc)
Parse the HGH abinit header. Example:

Hartwigsen-Goedecker-Hutter psp for Ne,  from PRB58, 3641 (1998)

    > 10   8  010605 zatom,zion,pspdat

    3 1   1 0 2001 0  pspcod,pspxc,lmax,lloc,mmax,r2well


#### _static_ oncvpsp_header(filename, ppdesc)
Parse the ONCVPSP abinit header. Example:

Li    ONCVPSP  r_core=  2.01  3.02

    > > 3.0000      3.0000      140504    zatom,zion,pspd

    > 8     2     1     4   600     0    pspcod,pspxc,lmax,lloc,mmax,r2well

    5.99000000  0.00000000  0.00000000    rchrg fchrg qchrg

        > 2     2     0     0     0    nproj
        > 0                 extension_switch

        0                        -2.5000025868368D+00 -1.2006906995331D+00

            1  0.0000000000000D+00  0.0000000000000D+00  0.0000000000000D+00
            2  1.0000000000000D-02  4.4140499497377D-02  1.9909081701712D-02


#### _static_ tm_header(filename, ppdesc)
Parse the TM abinit header. Example:

Troullier-Martins psp for element Fm         Thu Oct 27 17:28:39 EDT 1994
100.00000  14.00000    940714                zatom, zion, pspdat

> 1    1    3    0      2001    .00000      pspcod,pspxc,lmax,lloc,mmax,r2well
> 0   4.085   6.246    0   2.8786493        l,e99.0,e99.9,nproj,rcpsp
> .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
> 1   3.116   4.632    1   3.4291849        l,e99.0,e99.9,nproj,rcpsp
> .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
> 2   4.557   6.308    1   2.1865358        l,e99.0,e99.9,nproj,rcpsp
> .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
> 3  23.251  29.387    1   2.4776730        l,e99.0,e99.9,nproj,rcpsp
> .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
> 3.62474762267880     .07409391739104    3.07937699839200   rchrg,fchrg,qchrg


### _class_ pymatgen.io.abinit.pseudos.NcAbinitPseudo(path, header)
Bases: `NcPseudo`, `AbinitPseudo`

Norm-conserving pseudopotential in the Abinit format.


* **Parameters**


    * **path** – Filename.


    * **header** – `AbinitHeader` instance.



#### _property_ Z()
The atomic number of the atom.


#### _property_ Z_val()
Number of valence electrons.


#### _property_ l_local()
Angular momentum used for the local part.


#### _property_ l_max()
Maximum angular momentum.


#### _property_ nlcc_radius()
Radius at which the core charge vanish (i.e. cut-off in a.u.).
Returns 0.0 if nlcc is not used.


#### _property_ summary()
Summary line reported in the ABINIT header.


### _class_ pymatgen.io.abinit.pseudos.NcPseudo()
Bases: `object`

Abstract class defining the methods that must be implemented
by the concrete classes representing norm-conserving pseudopotentials.


#### _property_ has_nlcc()
True if the pseudo is generated with non-linear core correction.


#### _abstract property_ nlcc_radius()
Radius at which the core charge vanish (i.e. cut-off in a.u.).
Returns 0.0 if nlcc is not used.


#### _property_ rcore()
Radius of the pseudization sphere in a.u.


### _class_ pymatgen.io.abinit.pseudos.PawAbinitHeader(summary, \*\*kwargs)
Bases: `AbinitHeader`

The abinit header found in the PAW pseudopotential files.


#### _static_ paw_header(filename, ppdesc)
Parse the PAW abinit header. Examples:

Paw atomic data for element Ni - Generated by AtomPAW (N. Holzwarth) + AtomPAW2Abinit v3.0.5

    > 28.000  18.000 20061204               : zatom,zion,pspdat
    > 7  7  2 0   350 0.                    : pspcod,pspxc,lmax,lloc,mmax,r2well

    paw3 1305

        5 13                                  : basis_size,lmn_size

    0 0 1 1 2                              : orbitals
    3                                      : number_of_meshes
    1 3  350 1.1803778368E-05 3.5000000000E-02 : mesh 1, type,size,rad_step[,log_step]
    2 1  921 2.500000000000E-03                : mesh 2, type,size,rad_step[,log_step]
    3 3  391 1.1803778368E-05 3.5000000000E-02 : mesh 3, type,size,rad_step[,log_step]

    > 2.3000000000                          : r_cut(SPH)

    2 0.

Another format:

C  (US d-loc) - PAW data extracted from US-psp (D.Vanderbilt) - generated by USpp2Abinit v2.3.0

    > > 6.000   4.000 20090106               : zatom,zion,pspdat

    > 7 11  1 0   560 0.                    : pspcod,pspxc,lmax,lloc,mmax,r2well

    paw4 2230

        4  8                                  : basis_size,lmn_size

    0 0 1 1                                : orbitals
    5                                      : number_of_meshes
    1 2  560 1.5198032759E-04 1.6666666667E-02 : mesh 1, type,size,rad_step[,log_step]
    2 2  556 1.5198032759E-04 1.6666666667E-02 : mesh 2, type,size,rad_step[,log_step]
    3 2  576 1.5198032759E-04 1.6666666667E-02 : mesh 3, type,size,rad_step[,log_step]
    4 2  666 1.5198032759E-04 1.6666666667E-02 : mesh 4, type,size,rad_step[,log_step]
    5 2  673 1.5198032759E-04 1.6666666667E-02 : mesh 5, type,size,rad_step[,log_step]

    > 1.5550009124                          : r_cut(PAW)

    3 0.                                   : shape_type,rshape

Yet nnother one:

Paw atomic data for element Si - Generated by atompaw v3.0.1.3 & AtomPAW2Abinit v3.3.1

    > 14.000   4.000 20120814               : zatom,zion,pspdat
    > 7      11  1 0   663 0.               : pspcod,pspxc,lmax,lloc,mmax,r2well

    paw5 1331

        4  8                                  : basis_size,lmn_size

    0 0 1 1                                : orbitals
    5                                      : number_of_meshes
    1 2  663 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 1, type,size,rad_step[,log_step]
    2 2  658 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 2, type,size,rad_step[,log_step]
    3 2  740 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 3, type,size,rad_step[,log_step]
    4 2  819 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 4, type,size,rad_step[,log_step]
    5 2  870 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 5, type,size,rad_step[,log_step]

    > 1.5669671236                          : r_cut(PAW)

    2 0.                                   : shape_type,rshape


### _class_ pymatgen.io.abinit.pseudos.PawAbinitPseudo(path, header)
Bases: `PawPseudo`, `AbinitPseudo`

Paw pseudopotential in the Abinit format.


* **Parameters**


    * **path** – Filename.


    * **header** – `AbinitHeader` instance.



#### _property_ paw_radius()
Radius of the PAW sphere in a.u.


#### _property_ supports_soc()
True if the pseudo can be used in a calculation with spin-orbit coupling.
Base classes should provide a concrete implementation that computes this value.


### _class_ pymatgen.io.abinit.pseudos.PawPseudo()
Bases: `object`

Abstract class that defines the methods that must be implemented
by the concrete classes representing PAW pseudopotentials.


#### _abstract property_ paw_radius()
Radius of the PAW sphere in a.u.


#### _property_ rcore()
Alias of paw_radius.


### _class_ pymatgen.io.abinit.pseudos.PawXmlSetup(filepath)
Bases: `Pseudo`, `PawPseudo`

Setup class for PawXml.


* **Parameters**

    **filepath** (*str*) – Path to the XML file.



#### _property_ Z()
The atomic number of the atom.


#### _property_ Z_val()
Number of valence electrons.


#### ae_core_density()
The all-electron radial density.


#### ae_partial_waves()
Dictionary with the AE partial waves indexed by state.


#### _property_ l_local()
Angular momentum used for the local part.


#### _property_ l_max()
Maximum angular momentum.


#### _property_ paw_radius()
Radius of the PAW sphere in a.u.


#### plot_densities(ax=None, \*\*kwargs)
Plot the PAW densities.


* **Parameters**

    **ax** – matplotlib `Axes` or None if a new figure should be created.



* **Returns**

    matplotlib figure


Keyword arguments controlling the display of the figure:

| kwargs

 | Meaning

 |
| ------------ | ------------------------------------------------------------------------------------------------- |  |  |  |  |  |  |  |  |
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

#### plot_projectors(ax=None, fontsize=12, \*\*kwargs)
Plot the PAW projectors.


* **Parameters**

    **ax** – matplotlib `Axes` or None if a new figure should be created.


Returns: matplotlib figure

Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
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

#### plot_waves(ax=None, fontsize=12, \*\*kwargs)
Plot the AE and the pseudo partial waves.


* **Parameters**


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **fontsize** – fontsize for legends and titles


Returns: matplotlib figure

Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
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

#### projector_functions()
Dictionary with the PAW projectors indexed by state.


#### pseudo_core_density()
The pseudized radial density.


#### _property_ pseudo_partial_waves()
Dictionary with the pseudo partial waves indexed by state.


#### root()
Root tree of XML.


#### _property_ summary()
String summarizing the most important properties.


#### _property_ supports_soc()
Here I assume that the ab-initio code can treat the SOC within the on-site approximation.


#### yield_figs(\*\*kwargs)
This function *generates* a predefined list of matplotlib figures with minimal input from the user.


### _class_ pymatgen.io.abinit.pseudos.Pseudo()
Bases: `MSONable`

Abstract base class defining the methods that must be
implemented by the concrete pseudo-potential sub-classes.


#### _abstract property_ Z(_: in_ )
The atomic number of the atom.


#### _abstract property_ Z_val(_: in_ )
Valence charge.


#### as_dict(\*\*kwargs)
Return dictionary for MSONable protocol.


#### _classmethod_ as_pseudo(obj)
Convert obj into a pseudo. Accepts:

>
> * Pseudo object.


> * string defining a valid path.


#### as_tmpfile(tmpdir=None)
Copy the pseudopotential to a temporary a file and returns a new pseudopotential object.
Useful for unit tests in which we have to change the content of the file.


* **Parameters**

    **tmpdir** – If None, a new temporary directory is created and files are copied here
    else tmpdir is used.



#### _property_ basename(_: st_ )
File basename.


#### compute_md5()
Compute and return MD5 hash value.


#### _property_ djrepo_path()
The path of the djrepo file. None if file does not exist.


#### _property_ element(_: [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element_ )
Pymatgen `Element`.


#### _property_ filepath(_: st_ )
Absolute path to pseudopotential file.


#### _classmethod_ from_dict(dct)
Build instance from dictionary (MSONable protocol).


#### _static_ from_file(filename)
Build an instance of a concrete Pseudo subclass from filename.
Note: the parser knows the concrete class that should be instantiated
Client code should rely on the abstract interface provided by Pseudo.


#### _property_ has_dojo_report()
True if the pseudo has an associated DOJO_REPORT section.


#### _property_ has_hints()
True if self provides hints on the cutoff energy.


#### hint_for_accuracy(accuracy='normal')
Returns a `Hint` object with the suggested value of ecut [Ha] and
pawecutdg [Ha] for the given accuracy.
ecut and pawecutdg are set to zero if no hint is available.


* **Parameters**

    **accuracy** – [“low”, “normal”, “high”]



#### _property_ isnc(_: boo_ )
True if norm-conserving pseudopotential.


#### _property_ ispaw(_: boo_ )
True if PAW pseudopotential.


#### _abstract property_ l_local(_: in_ )
Angular momentum used for the local part.


#### _abstract property_ l_max(_: in_ )
Maximum angular momentum.


#### md5()
MD5 hash value.


#### open_pspsfile(ecut=20, pawecutdg=None)
Calls Abinit to compute the internal tables for the application of the
pseudopotential part. Returns `PspsFile` object providing methods
to plot and analyze the data or None if file is not found or it’s not readable.


* **Parameters**


    * **ecut** – Cutoff energy in Hartree.


    * **pawecutdg** – Cutoff energy for the PAW double grid.



#### _abstract property_ summary(_: st_ )
String summarizing the most important properties.


#### _abstract property_ supports_soc()
True if the pseudo can be used in a calculation with spin-orbit coupling.
Base classes should provide a concrete implementation that computes this value.


#### _property_ symbol(_: st_ )
Element symbol.


#### to_str(verbose=0)
String representation.


#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead


#### _property_ type(_: st_ )
Type of pseudo.


### _exception_ pymatgen.io.abinit.pseudos.PseudoParseError()
Bases: [`ParseError`](pymatgen.io.core.md#pymatgen.io.core.ParseError)

Base Error class for the exceptions raised by `PseudoParser`.


### _class_ pymatgen.io.abinit.pseudos.PseudoParser()
Bases: `object`

Responsible for parsing pseudopotential files and returning pseudopotential objects.

Usage:

```default
pseudo = PseudoParser().parse("filename")
```


#### Error()
alias of `PseudoParseError`


#### parse(filename)
Read and parse a pseudopotential file. Main entry point for client code.


* **Returns**

    pseudopotential object or None if filename is not a valid pseudopotential file.



#### read_ppdesc(filename)
Read the pseudopotential descriptor from filename.


* **Returns**

    Pseudopotential descriptor. None if filename is not a valid pseudopotential file.



* **Raises**

    **PseudoParseError** –



#### scan_directory(dirname, exclude_exts=(), exclude_fnames=())
Analyze the files contained in directory dirname.


* **Parameters**


    * **dirname** – directory path


    * **exclude_exts** – list of file extensions that should be skipped.


    * **exclude_fnames** – list of file names that should be skipped.



* **Returns**

    List of pseudopotential objects.



### _class_ pymatgen.io.abinit.pseudos.PseudoTable(pseudos)
Bases: `Sequence`, `MSONable`

Define the pseudopotentials from the element table.
Individidual elements are accessed by name, symbol or atomic number.

For example, the following all retrieve iron:

print elements[26]
Fe
print elements.Fe
Fe
print elements.symbol(‘Fe’)
Fe
print elements.name(‘iron’)
Fe
print elements.isotope(‘Fe’)
Fe


* **Parameters**

    **pseudos** – List of pseudopotentials or filepaths.



#### all_combinations_for_elements(element_symbols)
Return a list with all the possible combination of pseudos
for the given list of element_symbols.
Each item is a list of pseudopotential objects.

Example:

```default
table.all_combinations_for_elements(["Li", "F"])
```


#### _property_ allnc()
True if all pseudos are norm-conserving.


#### _property_ allpaw()
True if all pseudos are PAW.


#### as_dict(\*\*kwargs)
Return dictionary for MSONable protocol.


#### _classmethod_ as_table(items)
Return an instance of `PseudoTable` from the iterable items.


#### _classmethod_ from_dict(d)
Build instance from dictionary (MSONable protocol).


#### _classmethod_ from_dir(top, exts=None, exclude_dirs='_\*')
Find all pseudos in the directory tree starting from top.


* **Parameters**


    * **top** – Top of the directory tree


    * **exts** – List of files extensions. if exts == “all_files”
    we try to open all files in top


    * **exclude_dirs** – Wildcard used to exclude directories.


return: `PseudoTable` sorted by atomic number Z.


#### get_pseudos_for_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Return the list of `Pseudo` objects to be used for this `Structure`.


* **Parameters**

    **structure** – pymatgen `Structure`.



* **Raises**


    * **ValueError** –


    * **multiple occurrences are present in the table.** –



#### is_complete(zmax=118)
True if table is complete i.e. all elements with Z < zmax have at least on pseudopotential.


#### print_table(stream=<_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>, filter_function=None)
A pretty ASCII printer for the periodic table, based on some filter_function.


* **Parameters**


    * **stream** – file-like object


    * **filter_function** – A filtering function that take a Pseudo as input and returns a boolean.
    For example, setting filter_function = lambda p: p.Z_val > 2 will print
    a periodic table containing only pseudos with Z_val > 2.



#### pseudo_with_symbol(symbol, allow_multi=False)
Return the pseudo with the given chemical symbol.


* **Parameters**


    * **symbols** – String with the chemical symbol of the element


    * **allow_multi** – By default, the method raises ValueError
    if multiple occurrences are found. Use allow_multi to prevent this.



* **Raises**

    **ValueError if symbol is not found**** or ****multiple occurrences are present and not allow_multi** –



#### pseudos_with_symbols(symbols)
Return the pseudos with the given chemical symbols.


* **Raises**

    **ValueError if one**** of ****the symbols is not found**** or ****multiple occurrences are present.** –



#### select(condition)
Select only those pseudopotentials for which condition is True.


* **Parameters**

    **condition** – Function that accepts a `Pseudo` object and returns True or False.



* **Returns**

    New PseudoTable instance with pseudos for which condition is True.



* **Return type**

    PseudoTable



#### select_family(family)
Return PseudoTable with element belonging to the specified family, e.g. family=”alkaline”.


#### select_rows(rows)
Return new class:PseudoTable object with pseudos in the given rows of the periodic table.
rows can be either a int or a list of integers.


#### select_symbols(symbols, ret_list=False)
Return a `PseudoTable` with the pseudopotentials with the given list of chemical symbols.


* **Parameters**


    * **symbols** – str or list of symbols
    Prepend the symbol string with “-”, to exclude pseudos.


    * **ret_list** – if True a list of pseudos is returned instead of a `PseudoTable`



#### sort_by_z()
Return a new `PseudoTable` with pseudos sorted by Z.


#### sorted(attrname, reverse=False)
Sort the table according to the value of attribute attrname.


* **Returns**

    PseudoTable object



* **Return type**

    New class



#### to_table(filter_function=None)
Return string with data in tabular form.


#### with_dojo_report()
Select pseudos containing the DOJO_REPORT section. Return new class:PseudoTable object.


#### _property_ zlist()
Ordered list with the atomic numbers available in the table.


### _class_ pymatgen.io.abinit.pseudos.RadialFunction(mesh, values)
Bases: `RadialFunction`

Radial Function class.

Create new instance of RadialFunction(mesh, values)


### pymatgen.io.abinit.pseudos.l2str(l_ang_mom)
Convert the angular momentum l (int) to string.


### pymatgen.io.abinit.pseudos.str2l(s)
Convert a string to the angular momentum l (int).


### pymatgen.io.abinit.pseudos.straceback()
Returns a string with the traceback.