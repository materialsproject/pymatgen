---
layout: default
title: pymatgen.io.abinit.variable.md
nav_exclude: true
---

# pymatgen.io.abinit.variable module

Support for Abinit input variables.


### _class_ pymatgen.io.abinit.variable.InputVariable(name, value, units='', valperline=3)
Bases: `object`

An Abinit input variable.


* **Parameters**


    * **name** – Name of the variable.


    * **value** – Value of the variable.


    * **units** – String specifying one of the units supported by Abinit. Default: atomic units.


    * **valperline** – Number of items printed per line.



#### _property_ basename()
Return the name trimmed of any dataset index.


#### _property_ dataset()
Return the dataset index in string form.


#### format_list(values, floatdecimal=0)
Format a list of values into a string.
The result might be spread among several lines.


#### _static_ format_list2d(values, floatdecimal=0)
Format a list of lists.


#### _static_ format_scalar(val, floatdecimal=0)
Format a single numerical value into a string
with the appropriate number of decimal.


#### get_value()
Return the value.


#### _property_ name()
Name of the variable.


#### _property_ units()
Return the units.