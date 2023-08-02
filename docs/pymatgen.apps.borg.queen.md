---
layout: default
title: pymatgen.apps.borg.queen.md
nav_exclude: true
---

# pymatgen.apps.borg.queen module

This module defines the BorgQueen class, which manages drones to assimilate
data using Python’s multiprocessing.


### _class_ pymatgen.apps.borg.queen.BorgQueen(drone, rootpath=None, number_of_drones=1)
Bases: `object`

The Borg Queen controls the drones to assimilate data in an entire
directory tree. Uses multiprocessing to speed up things considerably. It
also contains convenience methods to save and load data between sessions.


* **Parameters**


    * **drone** (*Drone*) – An implementation of
    [`pymatgen.apps.borg.hive.AbstractDrone`](pymatgen.apps.borg.hive.md#pymatgen.apps.borg.hive.AbstractDrone) to use for
    assimilation.


    * **rootpath** (*str*) – The root directory to start assimilation. Leave it
    as None if you want to do assimilation later, or is using the
    BorgQueen to load previously assimilated data.


    * **number_of_drones** (*int*) – Number of drones to parallelize over.
    Typical machines today have up to four processors. Note that you
    won’t see a 100% improvement with two drones over one, but you
    will definitely see a significant speedup of at least 50% or so.
    If you are running this over a server with far more processors,
    the speedup will be even greater.



#### get_data()
Returns an list of assimilated objects.


#### load_data(filename)
Load assimilated data from a file.


#### parallel_assimilate(rootpath)
Assimilate the entire subdirectory structure in rootpath.


#### save_data(filename)
Save the assimilated data to a file.


* **Parameters**

    **filename** (*str*) – filename to save the assimilated data to. Note
    that if the filename ends with gz or bz2, the relevant gzip
    or bz2 compression will be applied.



#### serial_assimilate(rootpath)
Assimilate the entire subdirectory structure in rootpath serially.


### pymatgen.apps.borg.queen.order_assimilation(args)
Internal helper method for BorgQueen to process assimilation.