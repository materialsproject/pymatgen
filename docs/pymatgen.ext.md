---
layout: default
title: pymatgen.ext.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.ext namespace


## pymatgen.ext.cod module

This module provides classes to interface with the Crystallography Open
Database. If you use data from the COD, please cite the following works (as
stipulated by the COD developers).

> Merkys, A., Vaitkus, A., Butkus, J., Okulič-Kazarinas, M., Kairys, V. &
> Gražulis, S. (2016) “COD::CIF::Parser: an error-correcting CIF parser for
> the Perl language”. Journal of Applied Crystallography 49.

> Gražulis, S., Merkys, A., Vaitkus, A. & Okulič-Kazarinas, M. (2015)
> “Computing stoichiometric molecular composition from crystal structures”.
> Journal of Applied Crystallography 48, 85-91.

> Gražulis, S., Daškevič, A., Merkys, A., Chateigner, D., Lutterotti, L.,
> Quirós, M., Serebryanaya, N. R., Moeck, P., Downs, R. T. & LeBail, A.
> (2012) “Crystallography Open Database (COD): an open-access collection of
> crystal structures and platform for world-wide collaboration”. Nucleic
> Acids Research 40, D420-D427.

> Grazulis, S., Chateigner, D., Downs, R. T., Yokochi, A. T., Quiros, M.,
> Lutterotti, L., Manakova, E., Butkus, J., Moeck, P. & Le Bail, A. (2009)
> “Crystallography Open Database - an open-access collection of crystal
> structures”. J. Appl. Cryst. 42, 726-729.

> Downs, R. T. & Hall-Wallace, M. (2003) “The American Mineralogist Crystal
> Structure Database”. American Mineralogist 88, 247-250.


### _class_ pymatgen.ext.cod.COD()
Bases: `object`

An interface to the Crystallography Open Database.

Blank __init__. No args required.


#### get_cod_ids(formula)
Queries the COD for all cod ids associated with a formula. Requires
mysql executable to be in the path.


* **Parameters**

    **formula** (*str*) – Formula.



* **Returns**

    List of cod ids.



#### get_structure_by_formula(formula: str, \*\*kwargs)
Queries the COD for structures by formula. Requires mysql executable to
be in the path.


* **Parameters**


    * **formula** (*str*) – Chemical formula.


    * **kwargs** – All kwargs supported by
    `pymatgen.core.structure.Structure.from_str()`.



* **Returns**

    Structure, “cod_id”: int, “sg”: “P n m a”}]



* **Return type**

    A list of dict of the format [{“structure”



#### get_structure_by_id(cod_id, \*\*kwargs)
Queries the COD for a structure by id.


* **Parameters**


    * **cod_id** (*int*) – COD id.


    * **kwargs** – All kwargs supported by
    `pymatgen.core.structure.Structure.from_str()`.



* **Returns**

    A Structure.



#### query(sql: str)
Perform a query.


* **Parameters**

    **sql** – SQL string



* **Returns**

    Response from SQL query.


## pymatgen.ext.matproj module

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


## pymatgen.ext.optimade module

Optimade support.


### _class_ pymatgen.ext.optimade.OptimadeRester(aliases_or_resource_urls: str | list[str] | None = None, refresh_aliases: bool = False, timeout: int = 5)
Bases: `object`

Class to call OPTIMADE-compliant APIs, see [https://optimade.org](https://optimade.org) and [1].

This class is ready to use but considered in-development and subject to change.

[1] Andersen, C.W., *et al*.

    OPTIMADE, an API for exchanging materials data.
    Sci Data 8, 217 (2021). [https://doi.org/10.1038/s41597-021-00974-z](https://doi.org/10.1038/s41597-021-00974-z)

OPTIMADE is an effort to provide a standardized interface to retrieve information
from many different materials science databases.

This is a client to retrieve structures from OPTIMADE v1 compliant endpoints. It
does not yet support all features of the OPTIMADE v1 specification but is intended
as a way to quickly search an endpoint in a way familiar to users of pymatgen without
needing to know the full OPTIMADE specification.

For advanced usage, please see the OPTIMADE documentation at optimade.org and
consider calling the APIs directly.

For convenience, known OPTIMADE endpoints have been given aliases in pymatgen to save
typing the full URL.

To get an up-to-date list aliases, generated from the current list of OPTIMADE providers
at optimade.org, call the refresh_aliases() method or pass refresh_aliases=True when
creating instances of this class.


* **Parameters**


    * **aliases_or_resource_urls** – the alias or structure resource URL or a list of


    * **URLs** (*aliases** or **resource*) –


    * **not** (*if providing the resource URL directly it should*) –


    * **index** (*be an*) –


    * **"v1/structures"** (*this interface can only currently access the*) –


    * **URL** (*information from the specified resource*) –


    * **refresh_aliases** – if True, use an up-to-date list of providers/aliases from the live


    * **https** (*list** of **OPTIMADE providers hosted at*) – //providers.optimade.org.


    * **timeout** – number of seconds before an attempted request is abandoned, a good


    * **providers** (*timeout is useful when querying many*) –


    * **offline** (*some** of **which may be*) –



#### aliases(_ = {'aflow': 'http://aflow.org/API/optimade/', 'cod': 'https://www.crystallography.net/cod/optimade', 'jarvis': 'https://jarvis.nist.gov/optimade/jarvisdft', 'mcloud.2dtopo': 'https://aiida.materialscloud.org/2dtopo/optimade', 'mcloud.curated-cofs': 'https://aiida.materialscloud.org/curated-cofs/optimade', 'mcloud.mc2d': 'https://aiida.materialscloud.org/mc2d/optimade', 'mcloud.mc3d': 'https://aiida.materialscloud.org/mc3d/optimade', 'mcloud.optimade-sample': 'https://aiida.materialscloud.org/optimade-sample/optimade', 'mcloud.pyrene-mofs': 'https://aiida.materialscloud.org/pyrene-mofs/optimade', 'mcloud.scdm': 'https://aiida.materialscloud.org/autowannier/optimade', 'mcloud.stoceriaitf': 'https://aiida.materialscloud.org/stoceriaitf/optimade', 'mcloud.tc-applicability': 'https://aiida.materialscloud.org/tc-applicability/optimade', 'mcloud.tin-antimony-sulfoiodide': 'https://aiida.materialscloud.org/tin-antimony-sulfoiodide/optimade', 'mp': 'https://optimade.materialsproject.org', 'mpds': 'https://api.mpds.io', 'nmd': 'https://nomad-lab.eu/prod/rae/optimade/', 'odbx': 'https://optimade.odbx.science', 'odbx.odbx_misc': 'https://optimade-misc.odbx.science', 'omdb.omdb_production': 'http://optimade.openmaterialsdb.se', 'oqmd': 'http://oqmd.org/optimade/', 'tcod': 'https://www.crystallography.net/tcod/optimade', 'twodmatpedia': 'http://optimade.2dmatpedia.org'_ )

#### describe()
Provides human-readable information about the resources being searched by the OptimadeRester.


#### get_snls(elements: list[str] | str | None = None, nelements: int | None = None, nsites: int | None = None, chemical_formula_anonymous: str | None = None, chemical_formula_hill: str | None = None, additional_response_fields: str | list[str] | set[str] | None = None)
Retrieve StructureNL from OPTIMADE providers.

A StructureNL is an object provided by pymatgen which combines Structure with
associated metadata, such as the URL is was downloaded from and any additional namespaced
data.

Not all functionality of OPTIMADE is currently exposed in this convenience method. To
use a custom filter, call get_structures_with_filter().


* **Parameters**


    * **elements** – List of elements


    * **nelements** – Number of elements, e.g. 4 or [2, 5] for the range >=2 and <=5


    * **nsites** – Number of sites, e.g. 4 or [2, 5] for the range >=2 and <=5


    * **chemical_formula_anonymous** – Anonymous chemical formula


    * **chemical_formula_hill** – Chemical formula following Hill convention


    * **additional_response_fields** – Any additional fields desired from the OPTIMADE API,


    * **dictionary.** (*these will be stored under the '_optimade' key in each StructureNL.data*) –


Returns: Dict of (Dict of StructureNLs keyed by that database’s id system) keyed by provider


#### get_snls_with_filter(optimade_filter: str, additional_response_fields: str | list[str] | set[str] | None = None)
Get structures satisfying a given OPTIMADE filter.


* **Parameters**


    * **optimade_filter** – An OPTIMADE-compliant filter


    * **additional_response_fields** – Any additional fields desired from the OPTIMADE API,


Returns: Dict of Structures keyed by that database’s id system


#### get_structures(elements: list[str] | str | None = None, nelements: int | None = None, nsites: int | None = None, chemical_formula_anonymous: str | None = None, chemical_formula_hill: str | None = None)
Retrieve Structures from OPTIMADE providers.

Not all functionality of OPTIMADE is currently exposed in this convenience method. To
use a custom filter, call get_structures_with_filter().


* **Parameters**


    * **elements** – List of elements


    * **nelements** – Number of elements, e.g. 4 or [2, 5] for the range >=2 and <=5


    * **nsites** – Number of sites, e.g. 4 or [2, 5] for the range >=2 and <=5


    * **chemical_formula_anonymous** – Anonymous chemical formula


    * **chemical_formula_hill** – Chemical formula following Hill convention


Returns: Dict of (Dict Structures keyed by that database’s id system) keyed by provider


#### get_structures_with_filter(optimade_filter: str)
Get structures satisfying a given OPTIMADE filter.


* **Parameters**

    **optimade_filter** – An OPTIMADE-compliant filter


Returns: Dict of Structures keyed by that database’s id system


#### mandatory_response_fields(_ = ('lattice_vectors', 'cartesian_site_positions', 'species', 'species_at_sites'_ )

#### refresh_aliases(providers_url='https://providers.optimade.org/providers.json')
Updates available OPTIMADE structure resources based on the current list of OPTIMADE
providers.


### _class_ pymatgen.ext.optimade.Provider(name, base_url, description, homepage, prefix)
Bases: `tuple`

Create new instance of Provider(name, base_url, description, homepage, prefix)


#### base_url()
Alias for field number 1


#### description()
Alias for field number 2


#### homepage()
Alias for field number 3


#### name()
Alias for field number 0


#### prefix()
Alias for field number 4