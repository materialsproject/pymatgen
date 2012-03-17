Overview
========

This page provides new users of the pymatgen code base with a quick overview of 
the pymatgen code base.

Pymatgen is structured in a highly object-oriented manner. Almost everything
(Element, Site, Structure, etc.) is an object.  Currently, the code is heavily
biased towards the representation and manipulation of crystals with periodic 
boundary conditions, though flexibility has been built in for molecules and other
materials.

The core modules are in the (yes, you guess it) pymatgen.core package. Given the 
importance of this package for the overall functioning of the code, we have 
provided a quick summary of the various modules here:

1. pymatgen.core.periodic_table : Everything begins here, where the Element and 
   Specie (Element with an oxidation state) objects are defined.  Unlike typical 
   implementations, pymatgen's Element object is rich, which means that each 
   Element contains many useful properties associated with it, including atomic 
   numbers, atomic masses, melting points, boiling points, just to name a few. 

2. pymatgen.core.lattice : This module defines a Lattice object, which 
   essentially defines the lattice vectors in three dimensions. The Lattice 
   object provides convenience methods for performing fractional to cartesian 
   coordinates and vice version, lattice parameter and angles computations, etc.
 
3. pymatgen.core.structure : Defines the Site, PeriodicSite, Structure and 
   Composition objects. A Site is essentially a coordinate point containing an 
   Element or Specie. A PeriodicSite contains a Lattice as well. A Structure is 
   simply a list of Sites having the same Lattice. Finally, a Composition is 
   mapping of Element/Specie to amounts.

Side-note : to_dict/from_dict
=============================

As you explore the code, you may notice that many of the objects have a to_dict 
property and a from_dict static method implemented.  For most of the non-basic
objects, we have designed pymatgen such that it is easy to save objects for 
subsequent use. While python does provide pickling functionality, pickle tends to
be extremely fragile with respect to code changes. Pymatgen's to_dict provide a
means to save your work in a more robust manner, which also has the added benefit
of being more readable. The dictionary representation is also particularly useful
for entering such objects into certain databases, such as MongoDb.

The output from a to_dict method is always json/yaml serializable. So if you 
want to save a structure, you may do the following:

::

   with open('structure.json','w') as f:
      json.dump(structure.to_dict, f)

Similarly, to get the structure back from a json, you can do the following to
restore the structure (or any object with a to_dict method) from the json as
follows:

::

   with open('structure.json', 'r') as f:
      d = json.load(f)
      structure = Structure.from_dict(d)

You may replace any of the above json commands with yaml in the PyYAML package
to create a yaml file instead. There are certain tradeoffs between the two 
choices. JSON is much more efficient as a format, with extremely fast read/write
speed, but is much less readable. YAML is much slower in terms of io, but is 
human readable.

Creating Structures
===================

For most applications, you will be creating and manipulating Structure objects. 
There are several ways to create these objects:

Creating a Structure manually
-----------------------------

This is generally the most painful method. Though sometimes necessary, it is 
seldom the method you would use.  An example of creating the basic silicon 
crystal is provided below:

::

   si = Element("Si")
   coords = list()
   coords.append([0,0,0])
   coords.append([0.75,0.5,0.75])
   lattice = Lattice.from_parameters(a = 3.84, b = 3.84, c = 3.84, alpha = 120, beta = 90, gamma = 60)
   struct = Structure(lattice, [si, si], coords)


Creating Structures using the pymatgen.io packages
--------------------------------------------------

More often, you would already have the Structure that you want in a 
Crystallographic Information Format (CIF) file or from VASP input and output 
files. 

Pymatgen provides convenient packages to parse such files to obtain a Structure 
as well as other information associated with these files.

For example, to create a Structure from a cif,

::

   from pymatgen.io.cifio import CifParser
   parser = CifParser("mycif.cif")
   structure = parser.get_structures()[0]

Another example, creating a Structure from a VASP POSCAR/CONTCAR file.

::

   from pymatgen.io.vaspio import Poscar
   poscar = Poscar.from_file("POSCAR")
   struct = poscar.struct

Many of these io packages also provide the means to write a Structure to various 
output formats, e.g. the CifWriter in pymatgen.io.cifio. In particular, the
pymatgen.io.vaspio_set provides a powerful way to generate complete sets of VASP 
input files from a Structure.

Things you can do with Structures
=================================

This section is a work in progress.  But just to give an overview of the kind of 
analysis you can do:

1. Modify Structures using either pymatgen.core.structure_modifier, or even 
   better, using the pymatgen.transformations and pymatgen.alchemy packages.
2. Analyse Structures. E.g., compute the Ewald sum using the 
   pymatgen.analysis.ewald package, compare two structures for similarity using 
   pymatgen.analysis.structure_fitter.


Example scripts
===============

Some example scripts have been provided in the scripts directory. In general, 
most file format conversions, manipulations and io can be done with a few quick 
lines of code. For example, to read a POSCAR and write a cif:

::

   from pymatgen.io.vaspio import Poscar
   from pymatgen.io.cifio import CifWriter

   p = Poscar.from_file('POSCAR')
   w = CifWriter(p.struct)
   w.write_file('mystructure.cif')

More examples will be added soon.

