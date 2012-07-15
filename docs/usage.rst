Overview
========

This page provides new users of the pymatgen code base with a quick overview of 
the pymatgen code base.

Pymatgen is structured in a highly object-oriented manner. Almost everything
(Element, Site, Structure, etc.) is an object.  Currently, the code is heavily
biased towards the representation and manipulation of crystals with periodic 
boundary conditions, though flexibility has been built in for molecules.

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
 
3. pymatgen.core.structure : Defines the Site, PeriodicSite, Structure,
   Molecule, and Composition objects. A Site is essentially a coordinate point
   containing an Element or Specie. A PeriodicSite contains a Lattice as well.
   A Structure and Molecule are simply a list of PeriodicSites and Site
   respectively. Finally, a Composition is mapping of Element/Specie to amounts.

All units in pymatgen are typically assumed to be in atomic units, i.e.,
angstroms for lengths, eV for energies, etc. However, most objects do not assume
any units per se and it should be perfectly fine for the most part no matter
what units are being used, as long as they are used consistently.

Side-note : to_dict / from_dict
===============================

As you explore the code, you may notice that many of the objects have a to_dict 
property and a from_dict static method implemented.  For most of the non-basic
objects, we have designed pymatgen such that it is easy to save objects for 
subsequent use. While python does provide pickling functionality, pickle tends
to be extremely fragile with respect to code changes. Pymatgen's to_dict provide
a means to save your work in a more robust manner, which also has the added
benefit of being more readable. The dict representation is also particularly
useful for entering such objects into certain databases, such as MongoDb.

The output from a to_dict method is always json/yaml serializable. So if you 
want to save a structure, you may do the following::

   with open('structure.json','w') as f:
      json.dump(structure.to_dict, f)

Similarly, to get the structure back from a json, you can do the following to
restore the structure (or any object with a to_dict method) from the json as
follows::

   with open('structure.json', 'r') as f:
      d = json.load(f)
      structure = Structure.from_dict(d)

You may replace any of the above json commands with yaml in the PyYAML package
to create a yaml file instead. There are certain tradeoffs between the two 
choices. JSON is much more efficient as a format, with extremely fast read/write
speed, but is much less readable. YAML is an order of magnitude or more slower
in terms of parsing, but is more human readable.

PMG JSON encoder/decoder
------------------------

.. versionadded:: 1.9.0

In version 1.9.0 of pymatgen, a brand new serialization framework is
implemented. Extensions of the standard Python JSONEncoder and JSONDecoder is
introduced to support pymatgen objects. Given that these coders depend on
certain new dict keys, they will only support pymatgen objects coming from
version >= 1.9.0.

The PMGJSONEncoder uses the to_dict API of pymatgen to generate the necessary
dict for converting into json. To use the PMGJSONEncoder, simply add it as the
*cls* kwarg when using json. For example,::

   json.dumps(object, cls=PMGJSONEncoder)

The PMGJSONDecoder depends on finding a "module" and "class" key in the dict in
order to decode the necessary python object. In general, the PMGJSONEncoder will
add these keys if they are not present, but for better long term stability
(e.g., there may be situations where to_dict is called directly rather than
through the encoder), the easiest way is to add the following to any to_dict
property::
    
   d['module'] = self.__class__.__module__
   d['class'] = self.__class__.__name__
        
To use the PMGJSONDecoder, simply specify it as the *cls* kwarg when using json
load, e.g.,

::

   json.loads(json_string, cls=PMGJSONDecoder)

The decoder is written in such a way that it supports nested list and dict of
pymatgen objects. When going through the nesting hirerachy, the decoder will
look for the highest level module/class names specified and convert those to
pymatgen objects.

Structures
==========

For most applications, you will be creating and manipulating Structure objects.
The construction of Molecule follows a similar API. There are several ways to
create these objects:

Creating a Structure manually
-----------------------------

This is generally the most painful method. Though sometimes necessary, it is 
seldom the method you would use.  An example of creating the basic silicon 
crystal is provided below::

   from pymatgen import Lattice, Structure
   
   coords = list()
   coords.append([0,0,0])
   coords.append([0.75,0.5,0.75])
   lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, 
                                     beta=90, gamma=60)
   struct = Structure(lattice, ["Si", "Si"], coords)

Note that both elements and species (elements with oxidation states) are
supported. So both "Fe" and "Fe2+" are valid specifications.

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
input files from a Structure. In general, most file format conversions can be
done with a few quick lines of code. For example, to read a POSCAR and write a
cif::

   from pymatgen.io.vaspio import Poscar
   from pymatgen.io.cifio import CifWriter

   p = Poscar.from_file('POSCAR')
   w = CifWriter(p.struct)
   w.write_file('mystructure.cif')


Things you can do with Structures
---------------------------------

This section is a work in progress.  But just to give an overview of the kind of 
analysis you can do:

1. Modify Structures using either pymatgen.core.structure_modifier, or even 
   better, using the pymatgen.transformations and pymatgen.alchemy packages.
2. Analyse Structures. E.g., compute the Ewald sum using the 
   pymatgen.analysis.ewald package, compare two structures for similarity using 
   pymatgen.analysis.structure_fitter.

.. _entries:

Entries - Basic analysis unit
=============================

Beyond the core Element, Site and Structure objects, most analyses within in
pymatgen (e.g., creating a PhaseDiagram) is performed using Entry objects. An 
Entry in its most basic form contains a calculated energy and a composition, 
and may optionally contain other input or calculated data. In most instances, 
you will use the ComputedEntry or ComputedStructureEntry objects defined in the 
pymatgen.entries.computed_entries module. ComputedEntry objects can be created 
by either manually parsing calculated data calculations, or by using the 
pymatgen.borg package.

.. _compatibility:

Compatibility - Mixing GGA and GGA+U runs
-----------------------------------------

The Ceder group has developed a scheme where by GGA and GGA+U calculations can
be "mixed" such that analyses may be performed using the type of calculation
most appropriate for each entry. For instance, to generate a Fe-P-O phase diagram,
metallic phases such as Fe and FexPy are most appropriately modelled using 
standard GGA, while a hubbard U should be applied for the oxides such as FexOy 
and FexPyOz.

In the pymatgen.io.vaspio_set module, pre-defined parameter sets have been coded
to allow users to generate VASP input files that are consistent with input 
parameters that are compatible with the Materials Project data. Users who wish to 
perform analysis using runs calculated using these parameters should post-process 
entries generated from these runs using the appropriate compatibility. For 
example, if a user wants to generate a phase diagram from a list of entries 
generated from Fe-P-O vasp runs, he should use the following procedure:

::

   from pymatgen.entries.compatibility import MaterialsProjectCompatibility
   from pymatgen.phasediagram.pdmaker import PhaseDiagram
   from pymatgen.phasediagram.plotter import PDPlotter
   
   # Get unprocessed_entries using pymatgen.borg or other means.
   
   # Process the entries for compatibility
   compat = MaterialsProjectCompatibility()
   processed_entries = compat.process_entries(unprocessed_entries)
     
   # These few lines generates the phase diagram using the ComputedEntries. 
   pd = PhaseDiagram(processed_entries)
   plotter = PDPlotter(pd)
   plotter.show()

pymatgen.borg - High-throughput data assimilation
=================================================

The borg package is still a work in progress, but a lot can already be done with
it. The basic concept is to provide a convenient means to
assimilate large quantities of data in a directory structure. For now, the main
application is the assimilation of entire directory structures of VASP 
calculations into usable pymatgen entries, which can then be used for phase 
diagram and other analyses.  The outline of how it works is as follows:

1. Drones are defined in the pymatgen.borg.hive module. A Drone is essentially
   an object which defines how a directory is parsed into a pymatgen object. For
   example, the VaspToComputedEntryDrone defines how a directory containing a 
   vasp run (with a vasprun.xml file) is converted into ComputedEntry.
2. The BorgQueen object in pymatgen.borg.queen module uses Drones to assimilate
   an entire subdirectory structure. Parallel processing is used where possible
   to speed up the process.

Simple example - Making a phase diagram
---------------------------------------

Let's say you want to make the Li-O phase diagram. You have calculated all
Li, O, and Li-O compounds you are interested in and the runs are in the directory
"Li-O_runs". You can then generate the phase diagram using the following few lines
of code:

::
   
   from pymatgen.borg.hive import VaspToComputedEntryDrone
   from pymatgen.borg.queen import BorgQueen
   from pymatgen.phasediagram.pdmaker import PhaseDiagram
   from pymatgen.phasediagram.plotter import PDPlotter
   
   # These three lines assimilate the data into ComputedEntries.
   drone = VaspToComputedEntryDrone()
   queen = BorgQueen(drone, "Li-O_runs", 2)   
   entries = queen.get_data()
   
   # It's a good idea to perform a save_data, especially if you just assimilated
   # a large quantity of data which took some time. This allows you to reload the
   # data using a BorgQueen initialized with only the drone argument and calling
   # queen.load_data("Li-O_entries.json")
   queen.save_data("Li-O_entries.json")
   
   # These few lines generates the phase diagram using the ComputedEntries. 
   pd = PhaseDiagram(entries)
   plotter = PDPlotter(pd)
   plotter.show()

In this example, neither Li nor O requires a Hubbard U. However, if you are making
a phase diagram from a mix of GGA and GGA+U entries, you may need to post-process
the assimilated entries with a Compatibility object before running the phase
diagram code. See earlier section on entries_ and compatibility_.

Another example - Calculating reaction energies
-----------------------------------------------

Another example of a cool thing you can do with the loaded entries is to calculate
reaction energies. For example, reusing the Li-O data we have saved in the above
step,

::
   
   from pymatgen.borg.hive import VaspToComputedEntryDrone
   from pymatgen.borg.queen import BorgQueen
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

The pymatgen.transformations package is the standard package for performing
transformations on structures. Many transformations are already supported today,
from simple transformations such as adding and removing sites, and replacing
species in a structure to more advanced one-to-many transformations such as
partially removing a fraction of a certain species from a structure using an
electrostatic energy criterion. The Transformation classes follow a strict API.
A typical usage is as follows:

::

   from pymatgen.io.cifio import CifParser
   from pymatgen.transformations.standard_transformations import RemoveSpecieTransformations
   
   # Read in a LiFePO4 structure from a cif.
   parser = CifParser('LiFePO4.cif')
   struct = parser.get_structures()[0]
   
   t = RemoveSpeciesTransformation(["Li"])
   modified_structure = t.apply_transformation(struct)

pymatgen.alchemy - High-throughput transformations
==================================================

The pymatgen.alchemy package is a framework for performing high-throughput (HT)
structure transformations. For example, it allows a user to define a series of
transformations to be applied to a set of structures, generating new structures
in the process. The framework is also designed to provide proper logging of all
changes performed on structures, with infinite undo. The main classes are:

1. pymatgen.alchemy.materials.TransformedStructure - Standard object
   representing a TransformedStructure. Takes in an input structure and a list
   of transformations as an input. Can also be generated from cifs and POSCARs.
2. pymatgen.alchemy.transmuters.TransformedStructureTransmuter - An example of
   a Transmuter class, which takes a list of structures, and apply a sequence
   of transformations on all of them.
   
Usage example - replace Fe with Mn and remove all Li in all structures:

::

   from pymatgen.alchemy.transmuters import TransformedStructureTransmuter
   from pymatgen.transformations.standard_transformations import SubstitutionTransformation, RemoveSpeciesTransformation

   trans = []
   trans.append(SubstitutionTransformation({"Fe":"Mn"}))
   trans.append(RemoveSpecieTransformation(["Lu"]))
   transmuter = TransformedStructureTransmuter.from_cifs(["MultiStructure.cif"], trans)
   structures = transmuter.get_transformed_structures()

pymatgen.matproj.rest - Integration with the Materials Project REST API
=======================================================================

.. versionadded:: 2.0.0

In version 2.0.0 of pymatgen, we introduced one of the most powerful and useful
tools yet - an adaptor to the Materials Project REST API. The Materials Project
REST API (currently in a limited beta) was introduced to provide a means for
users to programmatically query for materials data. This allows users to
efficiently perform structure manipulatino and analyses without going through
the web interface.

In parallel, we have coded in the pymatgen.matproj.rest module a MPRestAdaptor,
a user-friendly adaptor to interface with the MP REST API to obtain useful
pymatgen objects for further analyses.  To use the MP REST API, a user first
needs to apply for an api key at the Materials Project website. In the examples
below, the user's Materials Project API key is designated as "USER_API_KEY".

The MPRestAdaptor provides many convenience methods, but we will just highlight
a few key methods here.

To obtain information on a material with Materials Project Id 1234, one can use
the following::

   adaptor = MPRestAdaptor("USER_API_KEY")
   
   #Structure for material id
   structure = adaptor.get_structure_by_material_id(1234) 

   #Dos for material id
   dos = adaptor.get_dos_by_material_id(1234) 

   #Bandstructure for material id
   bandstructure = adaptor.get_bandstructure_by_material_id(1234) 

The MP REST interface also allows for query of data by formulas::

   #To get a list of data for all entries having formula Fe2O3   
   data = adaptor.get_data("Fe2O3")
   
   #To get the energies of all entries having formula Fe2O3   
   energies = adaptor.get_data("Fe2O3", "energy")

Finally, the MPRestAdaptor provides methods to obtain all entries in a
chemical system. Combined with the borg framework, this provides a
particularly powerful way to combine one's own calculations with Materials
Project data for analysis. The code below demonstrates the phase stability of
a new calculated material can be determined::

   from pymatgen.matproj.rest import MPRestAdaptor
   from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
   from pymatgen.apps.borg.queen import BorgQueen
   from pymatgen.entries.compatibility import MaterialsProjectCompatibility
   from pymatgen.phasediagram.pdmaker import PhaseDiagram
   from pymatgen.phasediagram.plotter import PDPlotter

   # Assimilate VASP calculations into ComputedEntry object. Let's assume that
   # the calculations are for a series of new LixFeyOz phases that we want to
   # know the phase stability.
   drone = VaspToComputedEntryDrone()
   queen = BorgQueen(drone, rootpath=".")
   entries = queen.get_data()
   
   # Obtain all existing Li-Fe-O phases using the Materials Project REST API
   adaptor = MPRestAdaptor("USER_API_KEY")
   mp_entries = adaptor.get_entries_in_chemsys(["Li", "Fe", "O"])
   
   # Combined entry from calculated run with Materials Project entries
   entries.extend(mp_entries)
   
   # Process entries using the MaterialsProjectCompatibility
   compat = MaterialsProjectCompatibility()
   entries = compat.process_entries(entries)
   
   # Generate and plot Li-Fe-O phase diagram
   pd = PhaseDiagram(entries)
   plotter = PDPlotter(pd)
   plotter.show()

Example scripts
===============

A good way to explore the functionality of pymatgen is to look at examples. We
have written some example scripts to perform some commonly desired
functionality, e.g., file format conversion, determining the spacegroup of a
structure, plotting the DOS of a VASP run, visualizing a structure using VTK,
etc. These example scripts can be found in the `scripts directory of pymatgen's
github repo <https://github.com/materialsproject/pymatgen/tree/master/scripts>`_
or the `downloaded source from PyPI <http://pypi.python.org/pypi/pymatgen>`_. 

More examples will be added to the scripts directory in future.

