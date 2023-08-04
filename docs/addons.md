---
layout: default
title: Add-ons
nav_order: 7
---

# Add-ons

With effect from v2022.0.3, pymatgen, pymatgen.analysis, pymatgen.ext and pymatgen.io are now
[namespace packages](http://packaging.python.org/guides/packaging-namespace-packages). You may refer to the
[contributing page](/contributing). for details on how to write such packages. This page serves as a universal
resource page to list known pymatgen add-ons.

It should be noted that the **pymatgen maintainers provide no guarantees whatsoever on the quality or reliability of
any of the add-ons** listed here. End users should make their
own assessment of the functionality and quality.

Please submit a pull request to update this page when if release a new add-on package.

## Add-ons for Analysis

* [pymatgen-analysis-diffusion](http://pypi.org/project/pymatgen-analysis-diffusion): Provides modules for diffusion analysis, including path determination for NEB calculations, analysis of MD trajectories (RDF, van Hove, Arrhenius plots, etc.). This package is maintained by the Materials Virtual Lab.
* [pymatgen-analysis-defects](https://pypi.org/project/pymatgen-analysis-defects): Provides functionality related to defect analysis. This package is maintained by Jimmy-Xuan Shen, and officially supported by the Materials Project.

## Add-ons for Input/Output

* [pymatgen-io-fleur](http://pypi.org/project/pymatgen-io-fleur): Provides modules for reading and writing files used by the [`fleur`](https://www.flapw.de/rel) DFT code. This package is maintained by the juDFT team.
* [pymatgen-io-openmm](https://github.com/orionarcher/pymatgen-io-openmm): Provides easy IO for performing molecular dynamics on solutions with OpenMM. This package is maintained by Orion Archer Cohen.

# External Tools

If you would like your own tool to be listed here, please submit a PR. For a more complete but less curated list, have a
look at [pymatgen dependents](https://github.com/materialsproject/pymatgen/network/dependents).

* [Atomate2](https://github.com/materialsproject/atomate2): atomate2 is a library of computational materials science workflows.
* [LobsterPy](https://github.com/JaGeo/LobsterPy): Automatically analyze [Lobster runs](https://cohp.de).
* [pymatviz](https://github.com/janosh/pymatviz): Complements `pymatgen` with additional plotting functionality for larger datasets common in materials informatics.
* [DiSCoVeR](https://github.com/sparks-baird/mat_discover): A materials discovery algorithm geared towards exploring
  high-performance candidates in new chemical spaces.
* [rxn-network](https://github.com/GENESIS-EFRC/reaction-network): Reaction Network is a Python package for predicting likely
  inorganic chemical reaction pathways using graph theory.
* [Matbench](https://github.com/materialsproject/matbench): Benchmarks for machine learning property prediction.
* [Matbench Discovery](https://github.com/janosh/matbench-discovery): Benchmark for machine learning crystal stability prediction.
* [matgl](https://github.com/materialsvirtuallab/matgl): Graph deep learning library for materials. Implements M3GNet and MEGNet in DGL and Pytorch with more to come.
* [chgnet](https://github.com/CederGroupHub/chgnet): Pretrained universal neural network potential for charge-informed atomistic modeling.
