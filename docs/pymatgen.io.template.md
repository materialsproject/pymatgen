---
layout: default
title: pymatgen.io.template.md
nav_exclude: true
---

# pymatgen.io.template module

This module defines a simple concrete implementation of the InputGenerator class that can be
used to facilitate writing large numbers of input files based on a template.


### _class_ pymatgen.io.template.TemplateInputGen()
Bases: [`InputGenerator`](pymatgen.io.core.md#pymatgen.io.core.InputGenerator)

Concrete implementation of InputGenerator that is based on a single template input
file with variables.

This class is provided as a low-barrier way to support new codes and to provide
an intuitive way for users to transition from manual scripts to pymatgen I/O
classes.


#### get_input_set(template: str | Path, variables: dict | None = None, filename: str = 'input.txt')

* **Parameters**


    * **template** – the input file template containing variable strings to be
    replaced.


    * **variables** – dict of variables to replace in the template. Keys are the
    text to replaced with the values, e.g. {“TEMPERATURE”: 298} will
    replace the text $TEMPERATURE in the template. See Python’s
    Template.safe_substitute() method documentation for more details.


    * **filename** – name of the file to be written.