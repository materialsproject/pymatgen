Change log
==========

v2021.3.4
---------
* **Backwards incompatible** Pymatgen root imports have been removed from 
  v2021.3.4 in preparation for a change to a more modular, extensible 
  architecture that will allow more developers to contribute. If your existing
  code uses `from pymatgen import <something>`, you will need to make
  modifications. The easiest way is to use an IDE to run a Search and Replace.
  First, replace any `from pymatgen import MPRester` with 
  `from pymatgen.ext.matproj import MPRester`. Then, replace 
  `from pymatgen import` with `from pymatgen.core import`. Alternative, if you
  are using a Mac command line, you can do::

    find . -name '*.py' | xargs sed -i "" 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '*.py' | xargs sed -i "" 's/from pymatgen import/from pymatgen.core import/g'

  From a Linux command line, you can do::

    find . -name '*.py' | xargs sed -i 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '*.py' | xargs sed -i 's/from pymatgen import/from pymatgen.core import/g'

  This should resolve most import errors and only a few more modifications may
  need to be done by hand. 
