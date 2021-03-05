Change log
==========

v2021.3.4
---------
* **Backwards incompatible**: Pymatgen root imports have been removed from
  v2021.3.4 in preparation for a change to a more modular, extensible
  architecture that will allow more developers to contribute.

  If your existing code uses `from pymatgen import <something>`, you will need to make
  modifications. The easiest way is to use an IDE to run a Search and Replace.
  First, replace any `from pymatgen import MPRester` with
  `from pymatgen.ext.matproj import MPRester`. Then, replace
  `from pymatgen import` with `from pymatgen.core import`. Alternatively, if you
  are using a Mac command line, you can do::

    find . -name '*.py' | xargs sed -i "" 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '*.py' | xargs sed -i "" 's/from pymatgen import/from pymatgen.core import/g'

  From a Linux command line, you can do::

    find . -name '*.py' | xargs sed -i 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '*.py' | xargs sed -i 's/from pymatgen import/from pymatgen.core import/g'

  This should resolve most import errors and only a few more modifications may
  need to be done by hand.

  Specifically, the following "convenience imports" have been removed in favor of
  their canonical import::

    from pymatgen import Composition  # now "from pymatgen.core.composition import Composition"
    from pymatgen import Lattice  # now "from pymatgen.core.lattice import Lattice"
    from pymatgen import SymmOp  # now "from pymatgen.core.operations import SymmOp"
    from pymatgen import DummySpecie, DummySpecies, Element, Specie, Species  # now "from pymatgen.core.periodic_table ..."
    from pymatgen import PeriodicSite, Site  # now "from pymatgen.core.sites ..."
    from pymatgen import IMolecule, IStructure, Molecule, Structure  # now "from pymatgen.core.structure ..."
    from pymatgen import ArrayWithUnit, FloatWithUnit, Unit  # now "from pymatgen.core.units ..."
    from pymatgen import Orbital, Spin  # now "from pymatgen.electronic_structure.core ..."
    from pymatgen import MPRester  # now "from pymatgen.ext.matproj ..."

