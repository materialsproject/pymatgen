---
layout: default
title: pymatgen.io.common.md
nav_exclude: true
---

# pymatgen.io.common module

Module for defining common data used and produced by atomistic simulation packages.


### _class_ pymatgen.io.common.VolumetricData(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), data, distance_matrix=None, data_aug=None)
Bases: `MSONable`

Simple volumetric object. Used to read LOCPOT/CHGCAR files produced by
vasp as well as cube files produced by other codes.


#### structure()
Structure associated with the Volumetric Data object

..attribute:: is_spin_polarized

> True if run is spin polarized

..attribute:: dim

> Tuple of dimensions of volumetric grid in each direction (nx, ny, nz).

..attribute:: data

> Actual data as a dict of {string: np.array}. The string are “total”
> and “diff”, in accordance to the output format of Vasp LOCPOT and
> CHGCAR files where the total spin density is written first, followed
> by the difference spin density.


#### ngridpts()
Total number of grid points in volumetric data.

Typically, this constructor is not used directly and the static
from_file constructor is used. This constructor is designed to allow
summation and other operations between VolumetricData objects.


* **Parameters**


    * **structure** – Structure associated with the volumetric data


    * **data** – Actual volumetric data. If the data is provided as in list format,
    it will be converted into an np.array automatically


    * **data_aug** – Any extra information associated with volumetric data
    (typically augmentation charges)


    * **distance_matrix** – A pre-computed distance matrix if available.
    Useful so pass distance_matrices between sums,
    short-circuiting an otherwise expensive operation.



#### copy()

* **Returns**

    Copy of Volumetric object



#### _classmethod_ from_cube(filename)
Initialize the cube object and store the data as data.


* **Parameters**

    **filename** (*str*) – of the cube to read



#### _classmethod_ from_hdf5(filename, \*\*kwargs)
Reads VolumetricData from HDF5 file.


* **Parameters**

    **filename** – Filename



* **Returns**

    VolumetricData



#### get_average_along_axis(ind)
Get the averaged total of the volumetric data a certain axis direction.
For example, useful for visualizing Hartree Potentials from a LOCPOT
file.


* **Parameters**

    **ind** (*int*) – Index of axis.



* **Returns**

    Average total along axis



#### get_axis_grid(ind)
Returns the grid for a particular axis.


* **Parameters**

    **ind** (*int*) – Axis index.



#### get_integrated_diff(ind, radius, nbins=1)
Get integrated difference of atom index ind up to radius. This can be
an extremely computationally intensive process, depending on how many
grid points are in the VolumetricData.


* **Parameters**


    * **ind** (*int*) – Index of atom.


    * **radius** (*float*) – Radius of integration.


    * **nbins** (*int*) – Number of bins. Defaults to 1. This allows one to
    obtain the charge integration up to a list of the cumulative
    charge integration values for radii for [radius/nbins,
    2 \* radius/nbins, ….].



* **Returns**

    Differential integrated charge as a np array of [[radius, value],
    …]. Format is for ease of plotting. E.g., plt.plot(data[:,0],
    data[:,1])



#### linear_add(other, scale_factor=1.0)
Method to do a linear sum of volumetric objects. Used by + and -
operators as well. Returns a VolumetricData object containing the
linear sum.


* **Parameters**


    * **other** (*VolumetricData*) – Another VolumetricData object


    * **scale_factor** (*float*) – Factor to scale the other data by.



* **Returns**

    VolumetricData corresponding to self + scale_factor \* other.



#### linear_slice(p1, p2, n=100)
Get a linear slice of the volumetric data with n data points from
point p1 to point p2, in the form of a list.


* **Parameters**


    * **p1** (*list*) – 3-element list containing fractional coordinates of the first point.


    * **p2** (*list*) – 3-element list containing fractional coordinates of the second point.


    * **n** (*int*) – Number of data points to collect, defaults to 100.



* **Returns**

    List of n data points (mostly interpolated) representing a linear slice of the
    data from point p1 to point p2.



#### scale(factor)
Scale the data in place by a factor.


#### _property_ spin_data()
data}.
Essentially, this provides the actual Spin.up and Spin.down data
instead of the total and diff. Note that by definition, a
non-spin-polarized run would have Spin.up data == Spin.down data.


* **Type**

    The data decomposed into actual spin data as {spin



#### to_cube(filename, comment=None)
Write the total volumetric data to a cube file format, which consists of two comment lines,
a header section defining the structure IN BOHR, and the data.


* **Parameters**


    * **filename** (*str*) – Name of the cube file to be written.


    * **comment** (*str*) – If provided, this will be added to the second comment line



#### to_hdf5(filename)
Writes the VolumetricData to a HDF5 format, which is a highly optimized
format for reading storing large data. The mapping of the VolumetricData
to this file format is as follows:

VolumetricData.data -> f[“vdata”]
VolumetricData.structure ->

> f[“Z”]: Sequence of atomic numbers
> f[“fcoords”]: Fractional coords
> f[“lattice”]: Lattice in the pymatgen.core.lattice.Lattice matrix

> > format

> f.attrs[“structure_json”]: String of json representation


* **Parameters**

    **filename** (*str*) – Filename to output to.



#### value_at(x, y, z)
Get a data value from self.data at a given point (x, y, z) in terms
of fractional lattice parameters. Will be interpolated using a
RegularGridInterpolator on self.data if (x, y, z) is not in the original
set of data points.


* **Parameters**


    * **x** (*float*) – Fraction of lattice vector a.


    * **y** (*float*) – Fraction of lattice vector b.


    * **z** (*float*) – Fraction of lattice vector c.



* **Returns**

    Value from self.data (potentially interpolated) correspondisng to
    the point (x, y, z).