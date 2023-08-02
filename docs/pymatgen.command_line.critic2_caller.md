---
layout: default
title: pymatgen.command_line.critic2_caller.md
nav_exclude: true
---

# pymatgen.command_line.critic2_caller module

This module implements an interface to the critic2 Bader analysis code.

For most Bader analysis purposes, users are referred to
pymatgen.command_line.bader_caller instead, this module is for advanced
usage requiring identification of critical points in the charge density.

This module depends on a compiled critic2 executable available in the path.
Please follow the instructions at [https://github.com/aoterodelaroza/critic2](https://github.com/aoterodelaroza/critic2)
to compile.

New users are *strongly* encouraged to read the critic2 manual first.

In brief,
\* critic2 searches for critical points in charge density
\* a critical point can be one of four types: nucleus, bond, ring
or cage
\* it does this by seeding locations for likely critical points
and then searching in these regions
\* there are two lists of critical points in the output, a list
of non-equivalent points (with in-depth information about the
field at those points), and a full list of points generated
by the appropriate symmetry operations
\* connectivity between these points is also provided when
appropriate (e.g. the two nucleus critical points linked to

> a bond critical point)


* critic2 can do many other things besides

If you use this module, please cite:

A. Otero-de-la-Roza, E. R. Johnson and V. Luaña,
Comput. Phys. Communications 185, 1007-1018 (2014)
([https://doi.org/10.1016/j.cpc.2013.10.026](https://doi.org/10.1016/j.cpc.2013.10.026))

A. Otero-de-la-Roza, M. A. Blanco, A. Martín Pendás and
V. Luaña, Comput. Phys. Communications 180, 157-166 (2009)
([https://doi.org/10.1016/j.cpc.2008.07.018](https://doi.org/10.1016/j.cpc.2008.07.018))


### _class_ pymatgen.command_line.critic2_caller.Critic2Analysis(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), stdout=None, stderr=None, cpreport=None, yt=None, zpsp=None)
Bases: `MSONable`

Class to process the standard output from critic2 into pymatgen-compatible objects.

This class is used to store results from the Critic2Caller.

To explore the bond graph, use the “structure_graph”
method, which returns a user-friendly StructureGraph
class with bonding information. By default, this returns
a StructureGraph with edge weights as bond lengths, but
can optionally return a graph with edge weights as any
property supported by the CriticalPoint class, such as
bond ellipticity.

This class also provides an interface to explore just the
non-symmetrically-equivalent critical points via the
critical_points attribute, and also all critical
points (via nodes dict) and connections between them
(via edges dict). The user should be familiar with critic2
before trying to understand these.

Indexes of nucleus critical points in the nodes dict are the
same as the corresponding sites in structure, with indices of
other critical points arbitrarily assigned.

Only one of (stdout, cpreport) required, with cpreport preferred
since this is a new, native JSON output from critic2.


* **Parameters**


    * **structure** – associated Structure


    * **stdout** – stdout from running critic2 in automatic
    mode


    * **stderr** – stderr from running critic2 in automatic
    mode


    * **cpreport** – json output from CPREPORT command


    * **yt** – json output from YT command


    * **(****dict****)** (*zpsp*) – Dict of element/symbol name to number of electrons


(ZVAL in VASP pseudopotential), with which to calculate charge transfer.
Optional.


#### get_critical_point_for_site(n: int)

* **Parameters**

    **n** (*int*) – Site index.


Returns: A CriticalPoint instance


#### get_volume_and_charge_for_site(n)

* **Parameters**

    **n** – Site index n.


Returns: A dict containing “volume” and “charge” keys,
or None if YT integration not performed


#### structure_graph(include_critical_points=('bond', 'ring', 'cage'))
A StructureGraph object describing bonding information
in the crystal.


* **Parameters**


    * **include_critical_points** – add DummySpecies for


    * **themselves** (*the critical points*) –


    * **of** (*a list*) –


    * **"nucleus"** –


    * **"bond"** –


    * **"ring"** –


    * **"cage"** –


    * **None** (*set to*) –


    * **disable** (*to*) –


Returns: a StructureGraph


### _class_ pymatgen.command_line.critic2_caller.Critic2Caller(input_script)
Bases: `object`

Class to call critic2 and store standard output for further processing.

Run Critic2 on a given input script.


* **Parameters**

    **input_script** – string defining the critic2 input



#### _classmethod_ from_chgcar(structure, chgcar=None, chgcar_ref=None, user_input_settings=None, write_cml=False, write_json=True, zpsp=None)
Run Critic2 in automatic mode on a supplied structure, charge
density (chgcar) and reference charge density (chgcar_ref).

The reason for a separate reference field is that in
VASP, the CHGCAR charge density only contains valence
electrons and may be missing substantial charge at
nuclei leading to misleading results. Thus, a reference
field is commonly constructed from the sum of AECCAR0
and AECCAR2 which is the total charge density, but then
the valence charge density is used for the final analysis.

If chgcar_ref is not supplied, chgcar will be used as the
reference field. If chgcar is not supplied, the promolecular
charge density will be used as the reference field – this can
often still give useful results if only topological information
is wanted.

User settings is a dictionary that can contain:
\* GRADEPS, float (field units), gradient norm threshold
\* CPEPS, float (Bohr units in crystals), minimum distance between

> critical points for them to be equivalent


* NUCEPS, same as CPEPS but specifically for nucleus critical
points (critic2 default is dependent on grid dimensions)


* NUCEPSH, same as NUCEPS but specifically for hydrogen nuclei
since associated charge density can be significantly displaced
from hydrogen nucleus


* EPSDEGEN, float (field units), discard critical point if any
element of the diagonal of the Hessian is below this value,
useful for discarding points in vacuum regions


* DISCARD, float (field units), discard critical points with field
value below this value, useful for discarding points in vacuum
regions


* SEED, list of strings, strategies for seeding points, default
is [‘WS 1’, ‘PAIR 10’] which seeds critical points by
sub-dividing the Wigner-Seitz cell and between every atom pair
closer than 10 Bohr, see critic2 manual for more options


* **Parameters**


    * **structure** – Structure to analyze


    * **chgcar** – Charge density to use for analysis. If None, will
    use promolecular density. Should be a Chgcar object or path (string).


    * **chgcar_ref** – Reference charge density. If None, will use
    chgcar as reference. Should be a Chgcar object or path (string).


    * **(****dict****)** (*user_input_settings*) – as explained above


    * **(****bool****)** (*write_json*) – Useful for debug, if True will write all
    critical points to a file ‘table.cml’ in the working directory
    useful for visualization


    * **(****bool****)** – Whether to write out critical points


and YT json. YT integration will be performed with this setting.
:param zpsp (dict): Dict of element/symbol name to number of electrons
(ZVAL in VASP pseudopotential), with which to properly augment core regions
and calculate charge transfer. Optional.


#### _classmethod_ from_path(path, suffix='', zpsp=None)
Convenience method to run critic2 analysis on a folder with typical VASP output files.

This method will:

1. Look for files CHGCAR, AECAR0, AECAR2, POTCAR or their gzipped
counterparts.

2. If AECCAR\* files are present, constructs a temporary reference
file as AECCAR0 + AECCAR2.

3. Runs critic2 analysis twice: once for charge, and a second time
for the charge difference (magnetization density).


* **Parameters**


    * **path** – path to folder to search in


    * **suffix** – specific suffix to look for (e.g. ‘.relax1’ for
    ‘CHGCAR.relax1.gz’)


    * **zpsp** – manually specify ZPSP if POTCAR not present



* **Returns**




### _class_ pymatgen.command_line.critic2_caller.CriticalPoint(index, type, frac_coords, point_group, multiplicity, field, field_gradient, coords=None, field_hessian=None)
Bases: `MSONable`

Access information about a critical point and the field values at that point.

Class to characterise a critical point from a topological
analysis of electron charge density.

Note this class is usually associated with a Structure, so
has information on multiplicity/point group symmetry.


* **Parameters**


    * **index** – index of point


    * **type** – type of point, given as a string


    * **coords** – Cartesian coordinates in Angstroms


    * **frac_coords** – fractional coordinates


    * **point_group** – point group associated with critical point


    * **multiplicity** – number of equivalent critical points


    * **field** – value of field at point (f)


    * **field_gradient** – gradient of field at point (grad f)


    * **field_hessian** – hessian of field at point (del^2 f)



#### _property_ ellipticity()
Most meaningful for bond critical points,
can be physically interpreted as e.g. degree
of pi-bonding in organic molecules. Consult
literature for more information.
Returns: The ellpiticity of the field at the critical point.


#### _property_ laplacian()
The Laplacian of the field at the critical point.


* **Type**

    Returns



#### _property_ type()
Instance of CriticalPointType.


* **Type**

    Returns



### _class_ pymatgen.command_line.critic2_caller.CriticalPointType(value)
Bases: `Enum`

Enum type for the different varieties of critical point.


#### bond(_ = 'bond_ )

#### cage(_ = 'cage_ )

#### nnattr(_ = 'nnattr_ )

#### nucleus(_ = 'nucleus_ )

#### ring(_ = 'ring_ )

### pymatgen.command_line.critic2_caller.get_filepath(filename, warning, path, suffix)

* **Parameters**


    * **filename** – Filename


    * **warning** – Warning message


    * **path** – Path to search


    * **suffix** – Suffixes to search.