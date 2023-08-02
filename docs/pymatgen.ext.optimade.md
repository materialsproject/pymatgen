---
layout: default
title: pymatgen.ext.optimade.md
nav_exclude: true
---

# pymatgen.ext.optimade module

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