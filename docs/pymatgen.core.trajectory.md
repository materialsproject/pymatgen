---
layout: default
title: pymatgen.core.trajectory.md
nav_exclude: true
---

# pymatgen.core.trajectory module

This module provides classes to define a simulation trajectory, which could come from
either relaxation or molecular dynamics.


### _class_ pymatgen.core.trajectory.Trajectory(species: list[str | [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies) | [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition)], coords: list[list[Vector3D]] | np.ndarray | list[np.ndarray], charge: int | float | None = None, spin_multiplicity: int | float | None = None, lattice: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice) | Matrix3D | list[[Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice)] | list[Matrix3D] | np.ndarray | None = None, \*, site_properties: SitePropsType | None = None, frame_properties: list[dict] | None = None, constant_lattice: bool | None = True, time_step: int | float | None = None, coords_are_displacement: bool = False, base_positions: list[list[Vector3D]] | np.ndarray | None = None)
Bases: `MSONable`

Trajectory of a geometry optimization or molecular dynamics simulation.

Provides basic functions such as slicing trajectory, combining trajectories, and
obtaining displacements.

In below, N denotes the number of sites in the structure, and M denotes the
number of frames in the trajectory.


* **Parameters**


    * **species** – shape (N,). List of species on each site. Can take in flexible
    input, including:
    i.  A sequence of element / species specified either as string

    > symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    > e.g., (3, 56, …) or actual Element or Species objects.


        1. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** – shape (M, N, 3). fractional coordinates of the sites.


    * **charge** – int or float. Charge of the system. This is only used for Molecule-based
    trajectories.


    * **spin_multiplicity** – int or float. Spin multiplicity of the system. This is only
    used for Molecule-based trajectories.


    * **lattice** – shape (3, 3) or (M, 3, 3). Lattice of the structures in the
    trajectory; should be used together with constant_lattice.
    If constant_lattice=True, this should be a single lattice that is
    common for all structures in the trajectory (e.g. in an NVT run).
    If constant_lattice=False, this should be a list of lattices,
    each for one structure in the trajectory (e.g. in an NPT run or a
    relaxation that allows changing the cell size). This is only used for
    Structure-based trajectories.


    * **site_properties** – Properties associated with the sites. This should be a
    list of M dicts for a single dict. If a list of dicts, each provides
    the site properties for a frame. Each value in a dict should be a
    sequence of length N, giving the properties of the N sites.
    For example, for a trajectory with M=2 and N=4, the
    site_properties can be: [{“magmom”:[5,5,5,5]}, {“magmom”:[5,5,5,5]}].
    If a single dict, the site properties in the dict apply to all frames
    in the trajectory. For example, for a trajectory with M=2 and N=4,
    {“magmom”:[2,2,2,2]} means that, through the entire trajectory,
    the magmom are kept constant at 2 for all four atoms.


    * **frame_properties** – Properties associated with the structure (e.g. total
    energy). This should be a sequence of M dicts, with each dict
    providing the properties for a frame. For example, for a trajectory with
    M=2, the frame_properties can be [{‘energy’:1.0}, {‘energy’:2.0}].


    * **constant_lattice** – Whether the lattice changes during the simulation.
    Should be used together with lattice. See usage there. This is only
    used for Structure-based trajectories.


    * **time_step** – Time step of MD simulation in femto-seconds. Should be None
    for a trajectory representing a geometry optimization.


    * **coords_are_displacement** – Whether coords are given in displacements
    (True) or positions (False). Note, if this is True, coords
    of a frame (say i) should be relative to the previous frame (i.e.
    i-1), but not relative to the base_position.


    * **base_positions** – shape (N, 3). The starting positions of all atoms in the
    trajectory. Used to reconstruct positions when converting from
    displacements to positions. Only needs to be specified if
    coords_are_displacement=True. Defaults to the first index of
    coords when coords_are_displacement=False.



#### as_dict()
Return the trajectory as a MSONable dict.


#### extend(trajectory: Trajectory)
Append a trajectory to the current one.

The lattice, coords, and all other properties are combined.


* **Parameters**

    **trajectory** – Trajectory to append.



#### _classmethod_ from_file(filename: str | Path, constant_lattice: bool = True, \*\*kwargs)
Create trajectory from XDATCAR or vasprun.xml file.


* **Parameters**


    * **filename** – Path to the file to read from.


    * **constant_lattice** – Whether the lattice changes during the simulation,
    such as in an NPT MD simulation.


    * **\*\*kwargs** – Additional kwargs passed to Trajectory constructor.



* **Returns**

    A trajectory from the file.



#### _classmethod_ from_molecules(molecules: list[[pymatgen.core.structure.Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)], \*\*kwargs)
Create trajectory from a list of molecules.

Note: Assumes no atoms removed during simulation.


* **Parameters**


    * **molecules** – pymatgen Molecule objects.


    * **\*\*kwargs** – Additional kwargs passed to Trajectory constructor.



* **Returns**

    A trajectory from the structures.



#### _classmethod_ from_structures(structures: list[[pymatgen.core.structure.Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure)], constant_lattice: bool = True, \*\*kwargs)
Create trajectory from a list of structures.

Note: Assumes no atoms removed during simulation.


* **Parameters**


    * **structures** – pymatgen Structure objects.


    * **constant_lattice** – Whether the lattice changes during the simulation,
    such as in an NPT MD simulation.


    * **\*\*kwargs** – Additional kwargs passed to Trajectory constructor.



* **Returns**

    A trajectory from the structures.



#### get_molecule(idx: int)
Get molecule at specified index.


* **Parameters**

    **idx** – Index of molecule.



* **Returns**

    A pymatgen Molecule object.



#### get_structure(idx: int)
Get structure at specified index.


* **Parameters**

    **idx** – Index of structure.



* **Returns**

    A pymatgen Structure object.



#### to_displacements()
Converts positions of trajectory into displacements between consecutive frames.

base_positions and coords should both be in fractional coords. Does
not work for absolute coords because the atoms are to be wrapped into the
simulation box.

This is the opposite operation of to_positions().


#### to_positions()
Convert displacements between consecutive frames into positions.

base_positions and coords should both be in fractional coords or
absolute coords.

This is the opposite operation of to_displacements().


#### write_Xdatcar(filename: str | Path = 'XDATCAR', system: str | None = None, significant_figures: int = 6)
Writes to Xdatcar file.

The supported kwargs are the same as those for the
Xdatcar_from_structs.get_string method and are passed through directly.


* **Parameters**


    * **filename** – Name of file to write.  It’s prudent to end the filename with
    ‘XDATCAR’, as most visualization and analysis software require this
    for autodetection.


    * **system** – Description of system (e.g. 2D MoS2).


    * **significant_figures** – Significant figures in the output file.