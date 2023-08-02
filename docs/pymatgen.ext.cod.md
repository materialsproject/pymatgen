---
layout: default
title: pymatgen.ext.cod.md
nav_exclude: true
---

# pymatgen.ext.cod module

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