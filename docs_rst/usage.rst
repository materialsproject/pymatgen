=====
Usage
=====

This page provides new users of the pymatgen code base with a quick overview of
the pymatgen code base. It should also be pointed out that there is an
:doc:`examples page </examples>` with many ipython notebook examples with
actual code demonstrating the use of the code. Learning from those examples
is the fastest way to get started.

Pymatgen is structured in a highly object-oriented manner. Almost everything
(Element, Site, Structure, etc.) is an object.  Currently, the code is heavily
biased towards the representation and manipulation of crystals with periodic
boundary conditions, though flexibility has been built in for molecules.

The core modules are in the (yes, you guess it) pymatgen.core package. Given the
importance of this package for the overall functioning of the code, we have
provided a quick summary of the various modules here:

1. :mod:`pymatgen.core.periodic_table`: Everything begins here, where the
   Element and Specie (Element with an oxidation state) objects are defined.
   Unlike typical implementations, pymatgen's Element object is rich,
   which means that each Element contains many useful properties associated
   with it, including atomic numbers, atomic masses, melting points,
   boiling points, just to name a few.

2. :mod:`pymatgen.core.lattice`: This module defines a Lattice object, which
   essentially defines the lattice vectors in three dimensions. The Lattice
   object provides convenience methods for performing fractional to cartesian
   coordinates and vice versa, lattice parameter and angles computations, etc.

3. :mod:`pymatgen.core.sites`: Defines the Site and PeriodicSite objects. A
   Site is essentially a coordinate point containing an Element or Specie. A
   PeriodicSite contains a Lattice as well.

4. :mod:`pymatgen.core.structure`: Defines the Structure and Molecule objects.
   A Structure and Molecule are simply a list of PeriodicSites and Site
   respectively.

5. :mod:`pymatgen.core.composition`: A Composition is simply a mapping of
   Element/Specie to amounts.

All units in pymatgen are typically assumed to be in atomic units, i.e.,
angstroms for lengths, eV for energies, etc. However, most objects do not
assume any units per se and it should be perfectly fine for the most part no
matter what units are being used, as long as they are used consistently.

Side-note : as_dict / from_dict
===============================

As you explore the code, you may notice that many of the objects have an as_dict
method and a from_dict static method implemented.  For most of the non-basic
objects, we have designed pymatgen such that it is easy to save objects for
subsequent use. While python does provide pickling functionality, pickle tends
to be extremely fragile with respect to code changes. Pymatgen's as_dict
provide a means to save your work in a more robust manner, which also has the
added benefit of being more readable. The dict representation is also
particularly useful for entering such objects into certain databases,
such as MongoDb. This as_dict specification is provided in the monty library,
which is a general python supplementary library arising from pymatgen.

The output from an as_dict method is always json/yaml serializable. So if you
want to save a structure, you may do the following::

    with open('structure.json','w') as f:
        json.dump(structure.as_dict(), f)

Similarly, to get the structure back from a json, you can do the following to
restore the structure (or any object with a as_dict method) from the json as
follows::

    with open('structure.json', 'r') as f:
        d = json.load(f)
        structure = Structure.from_dict(d)

You may replace any of the above json commands with yaml in the PyYAML package
to create a yaml file instead. There are certain tradeoffs between the two
choices. JSON is much more efficient as a format, with extremely fast
read/write speed, but is much less readable. YAML is an order of magnitude
or more slower in terms of parsing, but is more human readable.

MontyEncoder/Decoder
--------------------

Extensions of the standard Python JSONEncoder and JSONDecoder has been
implemented to support pymatgen objects. The MontyEncoder uses the as_dict
API of pymatgen to generate the necessary dict for converting into json. To
use the MontyEncoder, simply add it as the *cls* kwarg when using json.
For example,::

    json.dumps(object, cls=MontyEncoder)

The MontyDecoder depends on finding a "@module" and "@class" key in the dict
to decode the necessary python object. In general, the MontyEncoder will
add these keys if they are not present, but for better long term stability
(e.g., there may be situations where to_dict is called directly rather than
through the encoder), the easiest way is to add the following to any to_dict
property::

    d["@module"] = self.__class__.__module__
    d["@class"] = self.__class__.__name__

To use the MontyDecoder, simply specify it as the *cls* kwarg when using json
load, e.g.::

    json.loads(json_string, cls=MontyDecoder)

The decoder is written in such a way that it supports nested list and dict of
pymatgen objects. When going through the nesting hirerachy, the decoder will
look for the highest level module/class names specified and convert those to
pymatgen objects.

The MontyEncoder/Decoder also supports datetime and numpy arrays out of box.

Structures and Molecules
========================

For most applications, you will be creating and manipulating
Structure/Molecule objects. There are several ways to create these objects:

Creating a Structure manually
-----------------------------

This is generally the most painful method. Though sometimes necessary, it is
seldom the method you would use.  An example of creating the basic silicon
crystal is provided below::

    from pymatgen.core import Lattice, Structure, Molecule

    coords = [[0, 0, 0], [0.75,0.5,0.75]]
    lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120,
                                      beta=90, gamma=60)
    struct = Structure(lattice, ["Si", "Si"], coords)

    coords = [[0.000000, 0.000000, 0.000000],
              [0.000000, 0.000000, 1.089000],
              [1.026719, 0.000000, -0.363000],
              [-0.513360, -0.889165, -0.363000],
              [-0.513360, 0.889165, -0.363000]]
    methane = Molecule(["C", "H", "H", "H", "H"], coords)

Note that both elements and species (elements with oxidation states) are
supported. So both "Fe" and "Fe2+" are valid specifications.

Reading and writing Structures/Molecules
----------------------------------------

More often, you would already have the Structure/Molecule in one of many
typical formats used (e.g., the Cystallographic Information Format (CIF),
electronic structure code input / output, xyz, mol, etc.).

Pymatgen provides a convenient way to read structures and molecules via the
from_file and to methods::

    # Read a POSCAR and write to a CIF.
    structure = Structure.from_file("POSCAR")
    structure.to(filename="CsCl.cif")

    # Read an xyz file and write to a Gaussian Input file.
    methane = Molecule.from_file("methane.xyz")
    methane.to(filename="methane.gjf")

The format is automatically guessed from the filename.

For more fine-grained control over which parsed to use, you can specify
specific io packages. For example, to create a Structure from a cif::

    from pymatgen.io.cif import CifParser
    parser = CifParser("mycif.cif")
    structure = parser.get_structures()[0]

Another example, creating a Structure from a VASP POSCAR/CONTCAR file::

    from pymatgen.io.vasp import Poscar
    poscar = Poscar.from_file("POSCAR")
    structure = poscar.structure

Many of these io packages also provide the means to write a Structure to
various output formats, e.g. the CifWriter in :mod:`pymatgen.io.cif`. In
particular, the :mod:`pymatgen.io.vasp.sets` provides a powerful way to
generate complete sets of VASP input files from a Structure. In general,
most file format conversions can be done with a few quick lines of code. For
example, to read a POSCAR and write a cif::

    from pymatgen.io.vasp import Poscar
    from pymatgen.io.cif import CifWriter

    p = Poscar.from_file('POSCAR')
    w = CifWriter(p.structure)
    w.write_file('mystructure.cif')

For molecules, pymatgen has in-built support for XYZ and Gaussian input and
output files via the :mod:`pymatgen.io.xyz` and
:mod:`pymatgen.io.gaussian` respectively::

    from pymatgen.io.xyz import XYZ
    from pymatgen.io.gaussian import GaussianInput

    xyz = XYZ.from_file('methane.xyz')
    gau = GaussianInput(xyz.molecule,
                        route_parameters={'SP': "", "SCF": "Tight"})
    gau.write_file('methane.inp')

There is also support for more than 100 file types via the OpenBabel
interface. But that requires you to install openbabel with Python bindings.
Please see the :doc:`installation guide </installation>`.

Things you can do with Structures
---------------------------------

This section is a work in progress.  But just to give an overview of the kind of
analysis you can do:

1. Modify Structures directly or even better, using the :mod:`pymatgen
   .transformations` and :mod:`pymatgen.alchemy` packages.
2. Analyse Structures. E.g., compute the Ewald sum using the
   :mod:`pymatgen.analysis.ewald` package, compare two structures for
   similarity using :mod:`pymatgen.analysis.structure_matcher`.

It should be noted that Structure and Molecule are designed to be mutable. In
fact, they are the most basic mutable units (everything below in the class
hierarchy such as Element, Specie, Site, PeriodicSite, Lattice are immutable).
If you need guarantees of immutability for Structure/Molecule,
you should use the IStructure and IMolecule classes instead.

Modifying Structures or Molecules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pymatgen supports a highly Pythonic interface for modifying Structures and
Molecules. For example, you can change any site simply with::

    # Change the specie at site position 1 to a fluorine atom.
    structure[1] = "F"
    molecule[1] = "F"

    # Change species and coordinates (fractional assumed for Structures,
    # cartesian for Molecules)
    structure[1] = "Cl", [0.51, 0.51, 0.51]
    molecule[1] = "F", [1.34, 2, 3]

    # Structure/Molecule also supports typical list-like operators,
    # such as reverse, extend, pop, index, count.
    structure.reverse()
    molecule.reverse()

    structure.append("F", [0.9, 0.9, 0.9])
    molecule.append("F", [2.1, 3,.2 4.3])

There are also many typical transforms you can do on Structures. Here are
some examples::

    # Make a supercell
    structure.make_supercell([2, 2, 2])

    # Get a primitive version of the Structure
    structure.get_primitive_structure()

    # Interpolate between two structures to get 10 structures (typically for
    # NEB calculations.
    structure.interpolate(another_structure, nimages=10)

The above is just some examples of typical use cases. A lot more is possible
and you may explore the actual API doc for the structure and molecule classes.

.. _entries:

Entries - Basic analysis unit
=============================

Beyond the core Element, Site and Structure objects, most analyses within in
pymatgen (e.g., creating a PhaseDiagram) are performed using Entry objects. An
Entry in its most basic form contains a calculated energy and a composition,
and may optionally contain other input or calculated data. In most instances,
you will use the ComputedEntry or ComputedStructureEntry objects defined in
:mod:`pymatgen.entries.computed_entries`. ComputedEntry objects can be created
by either manually parsing calculated data calculations, or by using the
:mod:`pymatgen.apps.borg` package.

.. _compatibility:

Compatibility - Mixing GGA and GGA+U runs
-----------------------------------------

The Ceder group has developed a scheme where by GGA and GGA+U calculations can
be "mixed" such that analyses may be performed using the type of calculation
most appropriate for each entry. For instance, to generate a Fe-P-O phase
diagram, metallic phases such as Fe and FexPy are most appropriately modelled
using standard GGA, while a hubbard U should be applied for the oxides such
as FexOy and FexPyOz.

In the :mod:`pymatgen.io.vasp.sets` module, pre-defined parameter sets have
been coded to allow users to generate VASP input files that are consistent
with input parameters that are compatible with the Materials Project data.
Users who wish to perform analysis using runs calculated using these
parameters should post-process entries generated from these runs using the
appropriate compatibility. For example, if a user wants to generate a phase
diagram from a list of entries generated from Fe-P-O vasp runs,
he should use the following procedure::

   from pymatgen.entries.compatibility import MaterialsProjectCompatibility
   from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

   # Get unprocessed_entries using pymatgen.borg or other means.

   # Process the entries for compatibility
   compat = MaterialsProjectCompatibility()
   processed_entries = compat.process_entries(unprocessed_entries)

   # These few lines generates the phase diagram using the ComputedEntries.
   pd = PhaseDiagram(processed_entries)
   plotter = PDPlotter(pd)
   plotter.show()

pymatgen.io - Managing calculation inputs and outputs
=====================================================

The :mod:`pymatgen.io` module contains classes to facilitate writing input files
and parsing output files from a variety of computational codes, including VASP,
Q-Chem, LAMMPS, CP2K, AbInit, and many more.

The core class for managing inputs is the :class:`InputSet`. An :class:`InputSet` object contains
all the data necessary to write one or more input files for a calculation.
Specifically, every :class:`InputSet` has a `write_input()` method that writes all the
necessary files to a location you specify. There are also :class:`InputGenerator` classes
that yield :class:`InputSet` with settings tailored to specific calculation types (for example,
a structure relaxation). You can think of :class:`InputGenerator` classes as "recipes" for
accomplishing specific computational tasks, while :class:`InputSet` contain those recipes
applied to a specific system or structure.

Custom settings can be provided to :class:`InputGenerator` on instantiation. For example,
to construct an :class:`InputSet` for a VASP structure relaxation using default Materials
Project parameters, but change the `NSW` parameter from the default (99) to 500::

    from pymatgen.io.generators import MPRelaxGen

    input_gen = MPRelaxGen(user_incar_settings={"NSW": 500})
    vasp_input_set = input_gen.get_input_set(structure)
    vasp_input_set.write_input('/path/to/calc/directory')

You can also use `InputSet.from_directory()` to construct a pymatgen :class:`InputSet`
from a directory containing calculation inputs.

Many codes also contain classes for parsing output files into pymatgen objects that
inherit from :class:`InputFile`, which provides a standard interface for reading and
writing individual files.

Use of :class:`InputFile`, :class:`InputSet`, and :class:`InputGenerator` classes is
not yet fully implemented by all codes supported by pymatgen, so please refer to the
respective module documentation for each code for more details.

pymatgen.borg - High-throughput data assimilation
=================================================

The borg package is still a work in progress, but a lot can already be done with
it. The basic concept is to provide a convenient means to
assimilate large quantities of data in a directory structure. For now, the main
application is the assimilation of entire directory structures of VASP
calculations into usable pymatgen entries, which can then be used for phase
diagram and other analyses.  The outline of how it works is as follows:

1. Drones are defined in the :mod:`pymatgen.apps.borg.hive` module. A Drone
   is essentially an object which defines how a directory is parsed into a
   pymatgen object. For example, the VaspToComputedEntryDrone defines how a
   directory containing a vasp run (with a vasprun.xml file) is converted
   into ComputedEntry.
2. The BorgQueen object in :mod:`pymatgen.apps.borg.queen` module uses Drones
   to assimilate an entire subdirectory structure. Parallel processing is
   used where possible to speed up the process.

Simple example - Making a phase diagram
---------------------------------------

Let's say you want to make the Li-O phase diagram. You have calculated all
Li, O, and Li-O compounds you are interested in and the runs are in the
directory "Li-O_runs". You can then generate the phase diagram using the
following few lines of code::

   from pymatgen.borg.hive import VaspToComputedEntryDrone
   from pymatgen.borg.queen import BorgQueen
   from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

   # These three lines assimilate the data into ComputedEntries.
   drone = VaspToComputedEntryDrone()
   queen = BorgQueen(drone, "Li-O_runs", 2)
   entries = queen.get_data()

   # It's a good idea to perform a save_data, especially if you just assimilated
   # a large quantity of data which took some time. This allows you to reload
   # the data using a BorgQueen initialized with only the drone argument and
   # calling queen.load_data("Li-O_entries.json")
   queen.save_data("Li-O_entries.json")

   # These few lines generates the phase diagram using the ComputedEntries.
   pd = PhaseDiagram(entries)
   plotter = PDPlotter(pd)
   plotter.show()

In this example, neither Li nor O requires a Hubbard U. However, if you are
making a phase diagram from a mix of GGA and GGA+U entries, you may need to
post-process the assimilated entries with a Compatibility object before
running the phase diagram code. See earlier section on entries_ and
compatibility_.

Another example - Calculating reaction energies
-----------------------------------------------

Another example of a cool thing you can do with the loaded entries is to
calculate reaction energies. For example, reusing the Li-O data we have saved
in the above step::

   from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
   from pymatgen.apps.borg.queen import BorgQueen
   from pymatgen.analysis.reaction_calculator import ComputedReaction

   # These three lines assimilate the data into ComputedEntries.
   drone = VaspToComputedEntryDrone()
   queen = BorgQueen(drone)
   queen.load_data("Li-O_entries.json")
   entries = queen.get_data()

   #Extract the correct entries and compute the reaction.
   rcts = filter(lambda e: e.composition.reduced_formula in ["Li", "O2"], entries)
   prods = filter(lambda e: e.composition.reduced_formula == "Li2O", entries)
   rxn = ComputedReaction(rcts, prods)
   print rxn
   print rxn.calculated_reaction_energy

pymatgen.transformations
========================

The :mod:`pymatgen.transformations` package is the standard package for
performing transformations on structures. Many transformations are already
supported today, from simple transformations such as adding and removing
sites, and replacing species in a structure to more advanced one-to-many
transformations such as partially removing a fraction of a certain species
from a structure using an electrostatic energy criterion. The Transformation
classes follow a strict API. A typical usage is as follows::

   from pymatgen.io.cif import CifParser
   from pymatgen.transformations.standard_transformations import RemoveSpecieTransformations

   # Read in a LiFePO4 structure from a cif.
   parser = CifParser('LiFePO4.cif')
   struct = parser.get_structures()[0]

   t = RemoveSpeciesTransformation(["Li"])
   modified_structure = t.apply_transformation(struct)

pymatgen.alchemy - High-throughput transformations
==================================================

The :mod:`pymatgen.alchemy` package is a framework for performing
high-throughput (HT) structure transformations. For example, it allows a user
to define a series of transformations to be applied to a set of structures,
generating new structures in the process. The framework is also designed to
provide proper logging of all changes performed on structures,
with infinite undo. The main classes are:

1. :class:`pymatgen.alchemy.materials.TransformedStructure` - Standard object
   representing a TransformedStructure. Takes in an input structure and a list
   of transformations as an input. Can also be generated from cifs and POSCARs.
2. :class:`pymatgen.alchemy.transmuters.StandardTransmuter` - An example of
   a Transmuter class, which takes a list of structures, and apply a sequence
   of transformations on all of them.

Usage example - replace Fe with Mn and remove all Li in all structures::

   from pymatgen.alchemy.transmuters import CifTransmuter
   from pymatgen.transformations.standard_transformations import SubstitutionTransformation, RemoveSpeciesTransformation

   trans = []
   trans.append(SubstitutionTransformation({"Fe":"Mn"}))
   trans.append(RemoveSpecieTransformation(["Lu"]))
   transmuter = CifTransmuter.from_filenames(["MultiStructure.cif"], trans)
   structures = transmuter.transformed_structures

pymatgen.matproj.rest - Integration with the Materials Project REST API
=======================================================================

In version 2.0.0 of pymatgen, we introduced one of the most powerful and useful
tools yet - an adaptor to the Materials Project REST API. The Materials Project
REST API (simply Materials API) was introduced to provide a means for
users to programmatically query for materials data. This allows users to
efficiently perform structure manipulation and analyses without going through
the web interface.

In parallel, we have coded in the :mod:`pymatgen.ext.matproj` module a
MPRester, a user-friendly high-level interface to the Materials API to obtain
useful pymatgen objects for further analyses.  To use the Materials API,
your need to first register with the Materials Project and generate your API
key in your dashboard at https://www.materialsproject.org/dashboard. In the
examples below, the user's Materials API key is designated as "USER_API_KEY".

The MPRester provides many convenience methods, but we will just highlight
a few key methods here.

To obtain information on a material with Materials Project Id "mp-1234",
one can use the following::

    from pymatgen.ext.matproj import MPRester
    with MPRester("USER_API_KEY") as m:

        # Structure for material id
        structure = m.get_structure_by_material_id("mp-1234")

        # Dos for material id
        dos = m.get_dos_by_material_id("mp-1234")

        # Bandstructure for material id
        bandstructure = m.get_bandstructure_by_material_id("mp-1234")

The Materials API also allows for query of data by formulas::

    # To get a list of data for all entries having formula Fe2O3
    data = m.get_data("Fe2O3")

    # To get the energies of all entries having formula Fe2O3
    energies = m.get_data("Fe2O3", "energy")

Finally, the MPRester provides methods to obtain all entries in a
chemical system. Combined with the borg framework, this provides a
particularly powerful way to combine one's own calculations with Materials
Project data for analysis. The code below demonstrates the phase stability of
a new calculated material can be determined::

   from pymatgen.ext.matproj import MPRester
   from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
   from pymatgen.apps.borg.queen import BorgQueen
   from pymatgen.entries.compatibility import MaterialsProjectCompatibility
   from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

   # Assimilate VASP calculations into ComputedEntry object. Let's assume that
   # the calculations are for a series of new LixFeyOz phases that we want to
   # know the phase stability.
   drone = VaspToComputedEntryDrone()
   queen = BorgQueen(drone, rootpath=".")
   entries = queen.get_data()

   # Obtain all existing Li-Fe-O phases using the Materials Project REST API
   with MPRester("USER_API_KEY") as m:
       mp_entries = m.get_entries_in_chemsys(["Li", "Fe", "O"])

   # Combined entry from calculated run with Materials Project entries
   entries.extend(mp_entries)

   # Process entries using the MaterialsProjectCompatibility
   compat = MaterialsProjectCompatibility()
   entries = compat.process_entries(entries)

   # Generate and plot Li-Fe-O phase diagram
   pd = PhaseDiagram(entries)
   plotter = PDPlotter(pd)
   plotter.show()

The query method
----------------

For the most flexibility, you can also use the query method of the MPRester.
This method allows any kind of mongo query to be performed on the Materials
Project database. It also supports a simple string syntax with wild cards.
Examples are given below::

   from pymatgen.ext.matproj import MPRester

   with MPRester("USER_API_KEY") as m:

       # Get all energies of materials with formula "*2O".
       results = m.query("*2O", ['energy'])

       # Get the formulas and energies of materials with materials_id mp-1234
       # or with formula FeO.
       results = m.query("FeO mp-1234", ['pretty_formula', 'energy'])

       # Get all compounds of the form ABO3
       results = m.query("**O3", ['pretty_formula', 'energy'])

It is highly recommended that you consult the Materials API documentation at
http://bit.ly/materialsapi, which provides a comprehensive explanation of the
document schema used in the Materials Project and how best to query for the
relevant information you need.

Setting the PMG_MAPI_KEY in the config file
-------------------------------------------

MPRester can also read the API key via the pymatgen config file. Simply run::

    pmg config --add PMG_MAPI_KEY <USER_API_KEY>

to add this to the `.pmgrc.yaml`, and you can now call MPRester without any
arguments. This makes it much easier for heavy users of the Materials API to
use MPRester without having to constantly insert their API key in the scripts.
