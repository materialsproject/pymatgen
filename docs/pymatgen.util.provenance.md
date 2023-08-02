---
layout: default
title: pymatgen.util.provenance.md
nav_exclude: true
---

# pymatgen.util.provenance module

Classes and methods related to the Structure Notation Language (SNL).


### _class_ pymatgen.util.provenance.Author(name, email)
Bases: `Author`

An Author contains two fields: name and email. It is meant to represent
the author of a Structure or the author of a code that was applied to a Structure.

Create new instance of Author(name, email)


#### as_dict()
Returns: MSONable dict.


#### _static_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    Author



#### _static_ parse_author(author)
Parses an Author object from either a String, dict, or tuple.


* **Parameters**

    **author** – A String formatted as “NAME <[email@domain.com](mailto:email@domain.com)>”,
    (name, email) tuple, or a dict with name and email keys.



* **Returns**

    An Author object.



### _class_ pymatgen.util.provenance.HistoryNode(name, url, description)
Bases: `HistoryNode`

A HistoryNode represents a step in the chain of events that lead to a
Structure. HistoryNodes leave ‘breadcrumbs’ so that you can trace back how
a Structure was created. For example, a HistoryNode might represent pulling
a Structure from an external database such as the ICSD or CSD. Or, it might
represent the application of a code (e.g. pymatgen) to the Structure, with
a custom description of how that code was applied (e.g. a site removal
Transformation was applied).

A HistoryNode contains three fields:


#### name()
The name of a code or resource that this Structure encountered in
its history (String)


#### url()
The URL of that code/resource (String)


#### description()
A free-form description of how the code/resource is related to the
Structure (dict).

Create new instance of HistoryNode(name, url, description)


#### as_dict()
Returns: Dict.


#### _static_ from_dict(dct: dict[str, str])

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    HistoryNode



#### _static_ parse_history_node(h_node)
Parses a History Node object from either a dict or a tuple.


* **Parameters**

    **h_node** – A dict with name/url/description fields or a 3-element
    tuple.



* **Returns**

    History node.



### _class_ pymatgen.util.provenance.StructureNL(struct_or_mol, authors, projects=None, references='', remarks=None, data=None, history=None, created_at=None)
Bases: `object`

The Structure Notation Language (SNL, pronounced ‘snail’) is a container for a pymatgen
Structure/Molecule object with some additional fields for enhanced provenance.

It is meant to be imported/exported in a JSON file format with the following structure:


    * sites


    * lattice (optional)


    * about


            * created_at


            * authors


            * projects


            * references


            * remarks


            * data


            * history


* **Parameters**


    * **struct_or_mol** – A pymatgen.core.structure Structure/Molecule object


    * **authors** – *List* of {“name”:’’, “email”:’’} dicts,
    *list* of Strings as ‘John Doe <[johndoe@gmail.com](mailto:johndoe@gmail.com)>’,
    or a single String with commas separating authors


    * **projects** – List of Strings [‘Project A’, ‘Project B’]


    * **references** – A String in BibTeX format


    * **remarks** – List of Strings [‘Remark A’, ‘Remark B’]


    * **data** – A free form dict. Namespaced at the root level with an
    underscore, e.g. {“_materialsproject”: <custom data>}


    * **history** – List of dicts - [{‘name’:’’, ‘url’:’’, ‘description’:{}}]


    * **created_at** – A datetime object.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    Class



#### _classmethod_ from_structures(structures: Sequence[[Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure)], authors: Sequence[dict[str, str]], projects=None, references='', remarks=None, data=None, histories=None, created_at=None)
A convenience method for getting a list of StructureNL objects by
specifying structures and metadata separately. Some of the metadata
is applied to all of the structures for ease of use.


* **Parameters**


    * **structures** – A list of Structure objects


    * **authors** – *List* of {“name”:’’, “email”:’’} dicts,
    *list* of Strings as ‘John Doe <[johndoe@gmail.com](mailto:johndoe@gmail.com)>’,
    or a single String with commas separating authors


    * **projects** – List of Strings [‘Project A’, ‘Project B’]. This
    applies to all structures.


    * **references** – A String in BibTeX format. Again, this applies to all
    structures.


    * **remarks** – List of Strings [‘Remark A’, ‘Remark B’]


    * **data** – A list of free form dict. Namespaced at the root level
    with an underscore, e.g. {“_materialsproject”:<custom data>}
    . The length of data should be the same as the list of
    structures if not None.


    * **histories** – List of list of dicts - [[{‘name’:’’, ‘url’:’’,
    ‘description’:{}}], …] The length of histories should be the
    same as the list of structures if not None.


    * **created_at** – A datetime object



### pymatgen.util.provenance.is_valid_bibtex(reference: str)
Use pybtex to validate that a reference is in proper BibTeX format.


* **Parameters**

    **reference** – A String reference in BibTeX format.



* **Returns**

    Boolean indicating if reference is valid bibtex.