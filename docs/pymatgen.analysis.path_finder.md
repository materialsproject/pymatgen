---
layout: default
title: pymatgen.analysis.path_finder.md
nav_exclude: true
---

# pymatgen.analysis.path_finder module

This module finds diffusion paths through a structure based on a given
potential field.

If you use PathFinder algorithm for your research, please consider citing the
following work:

```default
Ziqin Rong, Daniil Kitchaev, Pieremanuele Canepa, Wenxuan Huang, Gerbrand
Ceder, The Journal of Chemical Physics 145 (7), 074112
```


### _class_ pymatgen.analysis.path_finder.ChgcarPotential(chgcar, smear=False, normalize=True)
Bases: `StaticPotential`

Implements a potential field based on the charge density output from VASP.


* **Parameters**


    * **chgcar** – Chgcar object based on a VASP run of the structure of
    interest (Chgcar.from_file(“CHGCAR”))


    * **smear** – Whether or not to apply a Gaussian smearing to the
    potential


    * **normalize** – Whether or not to normalize the potential to range
    from 0 to 1



### _class_ pymatgen.analysis.path_finder.FreeVolumePotential(struct, dim, smear=False, normalize=True)
Bases: `StaticPotential`

Implements a potential field based on geometric distances from atoms in the
structure - basically, the potential
is lower at points farther away from any atoms in the structure.


* **Parameters**


    * **struct** – Unit cell on which to base the potential


    * **dim** – Grid size for the potential


    * **smear** – Whether or not to apply a Gaussian smearing to the
    potential


    * **normalize** – Whether or not to normalize the potential to range
    from 0 to 1



### _class_ pymatgen.analysis.path_finder.MixedPotential(potentials, coefficients, smear=False, normalize=True)
Bases: `StaticPotential`

Implements a potential that is a weighted sum of some other potentials.


* **Parameters**


    * **potentials** – List of objects extending the StaticPotential superclass


    * **coefficients** – Mixing weights for the elements of the potentials list


    * **smear** – Whether or not to apply a Gaussian smearing to the potential


    * **normalize** – Whether or not to normalize the potential to range from
    0 to 1.



### _class_ pymatgen.analysis.path_finder.NEBPathfinder(start_struct, end_struct, relax_sites, v, n_images=20, mid_struct=None)
Bases: `object`

General pathfinder for interpolating between two structures, where the
interpolating path is calculated with the elastic band method with
respect to the given static potential for sites whose indices are given
in relax_sites, and is linear otherwise.


* **Parameters**


    * **start_struct** – Endpoint structures to interpolate


    * **end_struct** – Endpoint structures to interpolate


    * **relax_sites** – List of site indices whose interpolation paths should
    be relaxed


    * **v** – Static potential field to use for the elastic band relaxation


    * **n_images** – Number of interpolation images to generate


    * **mid_struct** – (optional) additional structure between the start and end structures to help.



#### _property_ images()
Returns a list of structures interpolating between the start and
endpoint structures.


#### interpolate()
Finds a set of n_images from self.s1 to self.s2, where all sites except
for the ones given in relax_sites, the interpolation is linear (as in
pymatgen.core.structure.interpolate), and for the site indices given
in relax_sites, the path is relaxed by the elastic band method within
the static potential V.

If a mid point is defined we will interpolate from s1–> mid –>s2
The final number of images will still be n_images.


#### plot_images(outfile)
Generates a POSCAR with the calculated diffusion path with respect to the first endpoint.
:param outfile: Output file for the POSCAR.


#### _static_ string_relax(start, end, V, n_images=25, dr=None, h=3.0, k=0.17, min_iter=100, max_iter=10000, max_tol=5e-06)
Implements path relaxation via the elastic band method. In general, the
method is to define a path by a set of points (images) connected with
bands with some elasticity constant k. The images then relax along the
forces found in the potential field V, counterbalanced by the elastic
response of the elastic band. In general the endpoints of the band can
be allowed to relax also to their local minima, but in this calculation
they are kept fixed.


* **Parameters**


    * **start** – Starting point of the path calculation given in discrete
    coordinates with respect to the grid in V.


    * **end** – Endpoints of the path calculation.


    * **V** – potential field through which to calculate the path


    * **n_images** – number of images used to define the path. In general
    anywhere from 20 to 40 seems to be good.


    * **dr** – Conversion ratio from discrete coordinates to real coordinates
    for each of the three coordinate vectors


    * **h** – Step size for the relaxation. h = 0.1 works reliably, but is
    slow. h=10 diverges with large gradients but for the types of
    gradients seen in CHGCARs, works pretty reliably


    * **k** – Elastic constant for the band (in real units, not discrete)


    * **min_iter** – Minimum number of iterations to perform. Defaults to 100.


    * **max_iter** – Number of optimization steps the string will
    take before exiting (even if unconverged). Defaults to 10000.


    * **max_tol** – Convergence threshold such that if the string moves by
    less than max_tol in a step, and at least min_iter steps have
    passed, the algorithm will terminate. Depends strongly on the
    size of the gradients in V, but 5e-6 works reasonably well for
    CHGCARs.



### _class_ pymatgen.analysis.path_finder.StaticPotential(struct, pot)
Bases: `object`

Defines a general static potential for diffusion calculations. Implements
grid-rescaling and smearing for the potential grid. Also provides a
function to normalize the potential from 0 to 1 (recommended).


* **Parameters**


    * **struct** – atomic structure of the potential


    * **pot** – volumentric data to be used as a potential



#### gaussian_smear(r)
Applies an isotropic Gaussian smear of width (standard deviation) r to
the potential field. This is necessary to avoid finding paths through
narrow minima or nodes that may exist in the field (although any
potential or charge distribution generated from GGA should be
relatively smooth anyway). The smearing obeys periodic
boundary conditions at the edges of the cell.

:param r - Smearing width in Cartesian coordinates, in the same units

    as the structure lattice vectors


#### get_v()
Returns the potential.


#### normalize()
Sets the potential range 0 to 1.


#### rescale_field(new_dim)
Changes the discretization of the potential field by linear
interpolation. This is necessary if the potential field
obtained from DFT is strangely skewed, or is too fine or coarse. Obeys
periodic boundary conditions at the edges of
the cell. Alternatively useful for mixing potentials that originally
are on different grids.


* **Parameters**

    **new_dim** – tuple giving the numpy shape of the new grid