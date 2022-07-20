Change log
==========

v2022.7.19
----------
This will be the final release with the pymatgen.analysis.defects 
module included in the standard pymatgen package. This release will 
include the older defects code by default, but can also be replaced with 
the newer defects code through installation of ``pymatgen-analysis-defects``.

Subsequent versions of pymatgen will require 
the additional installation of `pymatgen-analysis-defects <https://github.com/materialsproject/pymatgen-analysis-defects>`_ for all defect-related 
functionality via ``pip install pymatgen-analysis-defects``.

Relevant imports will still be from the ``pymatgen.analysis.defects`` namespace but the code will now be maintained and developed in this separate repository.

There will be significant changes to the defects code to support new functionality.Existing PyCDT users should use this version of pymatgen or older. Any questions 
about this change should be directed to Jimmy-Xuan Shen, @jmmshn.

For more information about other pymatgen "add-on" packages, please see 
`this page in our documentation <https://pymatgen.org/addons.html>`_.

* Preparation for the removal of the defects module, PR #2582 by @jmmshn 
