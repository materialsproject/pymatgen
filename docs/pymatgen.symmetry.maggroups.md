---
layout: default
title: pymatgen.symmetry.maggroups.md
nav_exclude: true
---

# pymatgen.symmetry.maggroups module

Magnetic space groups.


### _class_ pymatgen.symmetry.maggroups.MagneticSpaceGroup(\*args, \*\*kwargs)
Bases: `MagneticSpaceGroup`

Representation of a magnetic space group.

Initializes a MagneticSpaceGroup from its Belov, Neronova and
Smirnova (BNS) number supplied as a list or its label supplied
as a string. To create a magnetic structure in pymatgen, the
Structure.from_magnetic_spacegroup() method can be used, which
relies on this class.

The main difference between magnetic space groups and normal
crystallographic space groups is the inclusion of a time reversal
operator that acts on an atom’s magnetic moment. This is
indicated by a prime symbol (’) next to the respective symmetry
operation in its label, e.g. the standard crystallographic
space group Pnma has magnetic subgroups Pn’ma, Pnm’a, Pnma’,
Pn’m’a, Pnm’a’, Pn’ma’, Pn’m’a’.

The magnetic space groups are classified as one of 4 types
where G = magnetic space group, and F = parent crystallographic
space group:


1. G=F no time reversal, i.e. the same as corresponding

    crystallographic group


2. G=F+F1’, “grey” groups, where avg. magnetic moment is zero,

    e.g. a paramagnet in zero ext. mag. field


3. G=D+(F-D)1’, where D is an equi-translation subgroup of F of

    index 2, lattice translations do not include time reversal


4. G=D+(F-D)1’, where D is an equi-class subgroup of F of index 2

There are two common settings for magnetic space groups, BNS
and OG. In case 4, the BNS setting != OG setting, and so a
transformation to go between the two settings is required:
specifically, the BNS setting is derived from D, and the OG
setting is derived from F.

This means that the OG setting refers to the unit cell if magnetic
order is neglected, and requires multiple unit cells to reproduce
the full crystal periodicity when magnetic moments are present.
This does not make the OG setting, in general, useful for
electronic structure calculations and the BNS setting is preferred.
However, this class does contain information on the OG setting and
can be initialized from OG labels or numbers if required.

Conventions: ITC monoclinic unique axis b, monoclinic cell choice 1,
hexagonal axis for trigonal groups, origin choice 2 for groups with
more than one origin choice (ISO-MAG).

Raw data comes from ISO-MAG, ISOTROPY Software Suite, iso.byu.edu
[http://stokes.byu.edu/iso/magnetic_data.txt](http://stokes.byu.edu/iso/magnetic_data.txt)
with kind permission from Professor Branton Campbell, BYU

Data originally compiled from:
(1) Daniel B. Litvin, Magnetic Group Tables (International Union

> of Crystallography, 2013) www.iucr.org/publ/978-0-9553602-2-0.


1. C. J. Bradley and A. P. Cracknell, The Mathematical Theory of
Symmetry in Solids (Clarendon Press, Oxford, 1972).

See [http://stokes.byu.edu/iso/magneticspacegroupshelp.php](http://stokes.byu.edu/iso/magneticspacegroupshelp.php) for more
information on magnetic symmetry.


* **Parameters**

    **id** – BNS number supplied as list of 2 ints or BNS label as
    str or index as int (1-1651) to iterate over all space groups