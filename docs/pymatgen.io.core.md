---
layout: default
title: pymatgen.io.core.md
nav_exclude: true
---

# pymatgen.io.core module

This module defines the abstract interface for reading and writing calculation
inputs in pymatgen. The interface comprises a 3-tiered hierarchy of classes.


1. An InputFile object represents the contents of a single input file, e.g.
the INCAR. This class standardizes file read and write operations.


2. An InputSet is a dict-like container that maps filenames (keys) to file
contents (either strings or InputFile objects). This class provides a standard
write_input() method.


3. InputGenerator classes implement a get_input_set method that, when provided
with a structure, return an InputSet object with all parameters set correctly.
Calculation input files can be written to disk with the write_inputs method.

If you want to implement a new InputGenerator, please take note of the following:


1. You must implement a get_input_set method that returns an InputSet


2. All customization of calculation parameters should be done in the __init__
method of the InputGenerator. The idea is that the generator contains
the “recipe”, but nothing that is specific to a particular system. get_input_set
takes system-specific information (such as structure) and applies the recipe.


3. All InputGenerator must save all supplied args and kwargs as instance variables.
E.g., self.my_arg = my_arg and self.kwargs = kwargs in the __init__. This
ensures the as_dict and from_dict work correctly.


### _class_ pymatgen.io.core.InputFile()
Bases: `MSONable`

Abstract base class to represent a single input file. Note that use of this class
is optional; it is possible create an InputSet that does not rely on underlying
InputFile objects.

All InputFile classes must implement a get_string method, which is called by
write_file.

If InputFile classes implement an __init__ method, they must assign all arguments
to __init__ as attributes.


#### _classmethod_ from_file(path: str | Path)
Creates an InputFile object from a file.


* **Parameters**

    **path** – Filename to read, including path.



* **Returns**

    InputFile



#### _abstract classmethod_ from_str(contents: str)
Create an InputFile object from a string.


* **Parameters**

    **contents** – The contents of the file as a single string



* **Returns**

    InputFile



#### _abstract_ get_string()
Return a string representation of an entire input file.


#### write_file(filename: str | Path)
Write the input file.


* **Parameters**

    **filename** – The filename to output to, including path.



### _class_ pymatgen.io.core.InputGenerator()
Bases: `MSONable`

InputGenerator classes serve as generators for Input objects. They contain
settings or sets of instructions for how to create Input from a set of
coordinates or a previous calculation directory.


#### _abstract_ get_input_set()
Generate an InputSet object. Typically the first argument to this method
will be a Structure or other form of atomic coordinates.


### _class_ pymatgen.io.core.InputSet(inputs: dict[str | Path, str | InputFile] | None = None, \*\*kwargs)
Bases: `MSONable`, `MutableMapping`

Abstract base class for all InputSet classes. InputSet are dict-like
containers for all calculation input data.

Since InputSet inherits dict, it can be instantiated in the same manner,
or a custom __init__ can be provided. Either way, self should be
populated with keys that are filenames to be written, and values that are
InputFile objects or strings representing the entire contents of the file.

All InputSet must implement from_directory. Implementing the validate method
is optional.

Instantiate an InputSet.


* **Parameters**


    * **inputs** – The core mapping of filename: file contents that defines the InputSet data.
    This should be a dict where keys are filenames and values are InputFile objects
    or strings representing the entire contents of the file. If a value is not an
    InputFile object nor a str, but has a __str__ method, this str representation
    of the object will be written to the corresponding file. This mapping will
    become the .inputs attribute of the InputSet.


    * **\*\*kwargs** – Any kwargs passed will be set as class attributes e.g.
    InputSet(inputs={}, foo=’bar’) will make InputSet.foo == ‘bar’.



#### _classmethod_ from_directory(directory: str | Path)
Construct an InputSet from a directory of one or more files.


* **Parameters**

    **directory** – Directory to read input files from



#### validate()
A place to implement basic checks to verify the validity of an
input set. Can be as simple or as complex as desired.

Will raise a NotImplementedError unless overloaded by the inheriting class.


#### write_input(directory: str | Path, make_dir: bool = True, overwrite: bool = True, zip_inputs: bool = False)
Write Inputs to one or more files.


* **Parameters**


    * **directory** – Directory to write input files to


    * **make_dir** – Whether to create the directory if it does not already exist.


    * **overwrite** – Whether to overwrite an input file if it already exists.


    * **generate_inputs** (*Additional kwargs are passed to*) –


    * **zip_inputs** – If True, inputs will be zipped into a file with the
    same name as the InputSet (e.g., InputSet.zip)



### _exception_ pymatgen.io.core.ParseError()
Bases: `SyntaxError`

This exception indicates a problem was encountered during parsing due to unexpected formatting.