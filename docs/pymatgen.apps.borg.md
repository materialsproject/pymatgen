---
layout: default
title: pymatgen.apps.borg.md
nav_exclude: true
---

# pymatgen.apps.borg package

The borg package contains modules that assimilate large quantities of data into
pymatgen objects for analysis.



* [pymatgen.apps.borg.hive module](pymatgen.apps.borg.hive.md)


    * [`AbstractDrone`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.AbstractDrone)


        * [`AbstractDrone.assimilate()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.AbstractDrone.assimilate)


        * [`AbstractDrone.get_valid_paths()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.AbstractDrone.get_valid_paths)


    * [`GaussianToComputedEntryDrone`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone)


        * [`GaussianToComputedEntryDrone.as_dict()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone.as_dict)


        * [`GaussianToComputedEntryDrone.assimilate()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone.assimilate)


        * [`GaussianToComputedEntryDrone.from_dict()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone.from_dict)


        * [`GaussianToComputedEntryDrone.get_valid_paths()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone.get_valid_paths)


    * [`SimpleVaspToComputedEntryDrone`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.SimpleVaspToComputedEntryDrone)


        * [`SimpleVaspToComputedEntryDrone.as_dict()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.SimpleVaspToComputedEntryDrone.as_dict)


        * [`SimpleVaspToComputedEntryDrone.assimilate()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.SimpleVaspToComputedEntryDrone.assimilate)


        * [`SimpleVaspToComputedEntryDrone.from_dict()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.SimpleVaspToComputedEntryDrone.from_dict)


    * [`VaspToComputedEntryDrone`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone)


        * [`VaspToComputedEntryDrone.as_dict()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone.as_dict)


        * [`VaspToComputedEntryDrone.assimilate()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone.assimilate)


        * [`VaspToComputedEntryDrone.from_dict()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone.from_dict)


        * [`VaspToComputedEntryDrone.get_valid_paths()`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone.get_valid_paths)


* [pymatgen.apps.borg.queen module](pymatgen.apps.borg.queen.md)


    * [`BorgQueen`](pymatgen.apps.borg.queen.md#pymatgen.apps.borg.queen.BorgQueen)


        * [`BorgQueen.get_data()`](pymatgen.apps.borg.queen.md#pymatgen.apps.borg.queen.BorgQueen.get_data)


        * [`BorgQueen.load_data()`](pymatgen.apps.borg.queen.md#pymatgen.apps.borg.queen.BorgQueen.load_data)


        * [`BorgQueen.parallel_assimilate()`](pymatgen.apps.borg.queen.md#pymatgen.apps.borg.queen.BorgQueen.parallel_assimilate)


        * [`BorgQueen.save_data()`](pymatgen.apps.borg.queen.md#pymatgen.apps.borg.queen.BorgQueen.save_data)


        * [`BorgQueen.serial_assimilate()`](pymatgen.apps.borg.queen.md#pymatgen.apps.borg.queen.BorgQueen.serial_assimilate)


    * [`order_assimilation()`](pymatgen.apps.borg.queen.md#pymatgen.apps.borg.queen.order_assimilation)