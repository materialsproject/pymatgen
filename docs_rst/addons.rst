Pymatgen Add-ons and External Tools
===================================

Add-ons
-------

With effect from v2022.0.3, pymatgen, pymatgen.analysis, pymatgen.ext and and pymatgen.io are now
`namespace packages <http://packaging.python.org/guides/packaging-namespace-packages>`_. You may refer to the
:doc:`contributing page </contributing>`. for details on how to write such packages. This page serves as a universal
resource page to list known pymatgen add-ons.

It should be noted that the **pymatgen maintainers provide no guarantees whatsoever on the quality or reliability of
any of the add-ons** listed here. End users should make their
own assessment of the functionality and quality.

Please submit a pull request to update this page when if release a new add-on package.

Add-ons for Analysis
~~~~~~~~~~~~~~~~~~~~

* `pymatgen-analysis-diffusion <http://pypi.org/project/pymatgen-analysis-diffusion>`_: Provides modules for diffusion
  analysis, including path determination for NEB calculations, analysis of MD trajectories (RDF, van Hove, Arrhenius
  plots, etc.). This package is maintained by the Materials Virtual Lab.

* `pymatgen-analysis-defects <https://pypi.org/project/pymatgen-analysis-defects>`_: Provides functionality related to
  defect analysis. This package is maintained by Jimmy-Xuan Shen, and officially supported by the Materials Project.

Add-ons for Input/Output
~~~~~~~~~~~~~~~~~~~~~~~~

* `pymatgen-io-fleur <http://pypi.org/project/pymatgen-io-fleur>`_: Provides modules for reading and writing
  files used by the `fleur <https://www.flapw.de/rel>`_ DFT code. This package is maintained by the juDFT team.
  
Add-ons for Externals Services
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* None at present.

External Tools
--------------

If you would like your own tool to be listed here, please `submit a PR <https://github.com/materialsproject/pymatgen/edit/master/docs_rst/addons.rst>`_! For a more complete but less curated list, have a
look at `pymatgen dependents <https://github.com/materialsproject/pymatgen/network/dependents>`_.

* `QuAcc <https://github.com/arosen93/quacc>`_: A platform to enable high-throughput, database-driven quantum
  chemistry and computational materials science.

* `LobsterPy <https://github.com/JaGeo/LobsterPy>`_: Automatically analyze Lobster (http://www.cohp.de) runs.

* `pymatviz <https://github.com/janosh/pymatviz>`_: Complements ``pymatgen`` with additional plotting
  functionality for larger datasets common in materials informatics.

* `M3GNet <https://github.com/materialsvirtuallab/m3gnet>`_: Materials graph network with 3-body interactions featuring
  a DFT surrogate crystal relaxer and property predictor.

* `DiSCoVeR <https://github.com/sparks-baird/mat_discover>`_: A materials discovery algorithm geared towards exploring
  high-performance candidates in new chemical spaces.

* `rxn-network <https://github.com/GENESIS-EFRC/reaction-network>`_: Reaction Network is a Python package for predicting likely
  inorganic chemical reaction pathways using graph theory.

* `Matbench <https://github.com/materialsproject/matbench>`_: Benchmarks for materials science property prediction.
