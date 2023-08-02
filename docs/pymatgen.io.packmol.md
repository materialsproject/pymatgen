---
layout: default
title: pymatgen.io.packmol.md
nav_exclude: true
---

# pymatgen.io.packmol module

This module provides a pymatgen I/O interface to packmol.

This adopts the minimal core I/O interface (see pymatgen/io/core).
In this case, only a two classes are used. PackmolSet(InputSet) is the container
class that provides a run() method for running packmol locally.

PackmolBoxGen(InputGenerator) provides a recipe for packing molecules into a
box, which returns a PackmolSet object.

For the run() method to work, you need to install the packmol package
See [http://m3g.iqm.unicamp.br/packmol](http://m3g.iqm.unicamp.br/packmol) or
[http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml](http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml)
for download and setup instructions. Note that packmol versions prior to 20.3.0
do not support paths with spaces.
After installation, you may need to manually add the path of the packmol
executable to the PATH environment variable.


### _class_ pymatgen.io.packmol.PackmolBoxGen(tolerance: float = 2.0, seed: int = 1, control_params: dict | None = None, inputfile: str | Path = 'packmol.inp', outputfile: str | Path = 'packmol_out.xyz', stdoutfile: str | Path = 'packmol.stdout')
Bases: [`InputGenerator`](pymatgen.io.core.md#pymatgen.io.core.InputGenerator)

Generator for a Packmol InputSet that packs one or more molecules into a rectangular
simulation box.

Instantiate a PackmolBoxGen class. The init method defines simulations parameters
like filenames, random seed, tolerance, etc.


* **Parameters**


    * **tolerance** – Tolerance for packmol, in Å.


    * **seed** – Random seed for packmol. Use a value of 1 (default) for deterministic
    output, or -1 to generate a new random seed from the current time.


    * **inputfile** – Path to the input file. Default to ‘packmol.inp’.


    * **outputfile** – Path to the output file. Default to ‘output.xyz’.


    * **stdoutfile** – Path to the file where stdout will be recorded. Default to ‘packmol.stdout’



#### get_input_set(molecules: list[dict], box: list[float] | None = None)
Generate a Packmol InputSet for a set of molecules.


* **Parameters**

    **molecules** – A list of dict containing information about molecules to pack
    into the box. Each dict requires three keys:

    >
    > 1. ”name” - the structure name


    > 2. ”number” - the number of that molecule to pack into the box


    > 3. ”coords” - Coordinates in the form of either a Molecule object or

    >     a path to a file.



### Example

{“name”: “water”,

    “number”: 500,
    “coords”: “/path/to/input/file.xyz”}

    > box: A list of box dimensions xlo, ylo, zlo, xhi, yhi, zhi, in Å. If set to None

    >     (default), pymatgen will estimate the required box size based on the volumes of
    >     the provided molecules.


### _class_ pymatgen.io.packmol.PackmolSet(inputs: dict[str | Path, str | [InputFile](pymatgen.io.core.md#pymatgen.io.core.InputFile)] | None = None, \*\*kwargs)
Bases: [`InputSet`](pymatgen.io.core.md#pymatgen.io.core.InputSet)

InputSet for the Packmol software. This class defines several attributes related
to.

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

    **directory** (*str** | **Path*) – Directory to read input files from.



#### run(path: str | Path, timeout=30)
Run packmol and write out the packed structure.


* **Parameters**


    * **path** – The path in which packmol input files are located.


    * **timeout** – Timeout in seconds.



* **Raises**


    * **ValueError if packmol does not succeed in packing the box.** –


    * **TimeoutExpiredError if packmold does not finish within the timeout.** –