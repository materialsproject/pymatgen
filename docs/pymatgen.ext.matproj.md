---
layout: default
title: pymatgen.ext.matproj.md
nav_exclude: true
---

# pymatgen.ext.matproj module

This module provides classes to interface with the Materials Project REST
API v2 to enable the creation of data structures and pymatgen objects using
Materials Project data.

To make use of the Materials API, you need to be a registered user of the
Materials Project, and obtain an API key by going to your dashboard at
[https://materialsproject.org/dashboard](https://materialsproject.org/dashboard).


### _exception_ pymatgen.ext.matproj.MPRestError()
Bases: `Exception`

Exception class for legacy MPRestAdaptor. Raised when query is malformed.


### _class_ pymatgen.ext.matproj.MPRester(\*args, \*\*kwargs)
Bases: `object`

A class to conveniently interface with the new and legacy Materials Project REST interface.

The recommended way to use MPRester is as a context manager to ensure
that sessions are properly closed after usage:

> with MPRester(“API_KEY”) as mpr:

>     docs = mpr.call_some_method()

MPRester uses the “requests” package, which provides HTTP connection
pooling. All connections are made via https for security.

For more advanced uses of the Materials API, please consult the API
documentation at [https://materialsproject.org/api](https://materialsproject.org/api) and [https://docs.materialsproject.org](https://docs.materialsproject.org).

This class handles the transition between old and new MP API, making it easy to switch between them
by passing a new (length 32) or old (15 <= length <= 17) API key. See [https://docs.materialsproject.org](https://docs.materialsproject.org)
for which API to use.


* **Parameters**


    * **\*args** – Pass through to either legacy or new MPRester.


    * **\*\*kwargs** – Pass through to either legacy or new MPRester.



### _class_ pymatgen.ext.matproj.TaskType(value)
Bases: `Enum`

task types available in legacy MP data.


#### GGAU_DEF(_ = 'GGA+U Deformation_ )

#### GGAU_LINE(_ = 'GGA+U NSCF Line_ )

#### GGAU_OPT(_ = 'GGA+U Structure Optimization_ )

#### GGAU_STATIC(_ = 'GGA+U Static_ )

#### GGAU_STATIC_DIEL(_ = 'GGA+U Static Dielectric_ )

#### GGAU_UNIFORM(_ = 'GGA+U NSCF Uniform_ )

#### GGA_DEF(_ = 'GGA Deformation_ )

#### GGA_LINE(_ = 'GGA NSCF Line_ )

#### GGA_OPT(_ = 'GGA Structure Optimization_ )

#### GGA_STATIC(_ = 'GGA Static_ )

#### GGA_STATIC_DIEL(_ = 'GGA Static Dielectric_ )

#### GGA_UNIFORM(_ = 'GGA NSCF Uniform_ )

#### LDA_STATIC_DIEL(_ = 'LDA Static Dielectric_ )

#### SCAN_OPT(_ = 'SCAN Structure Optimization_ )

### pymatgen.ext.matproj.get_chunks(sequence: Sequence[Any], size=1)

* **Parameters**


    * **sequence** (*Sequence**[**Any**]*) – Any sequence.


    * **size** (*int*) – Chunk length. Defaults to 1.



* **Returns**

    input sequence in chunks of length size.



* **Return type**

    list[Sequence[Any]]