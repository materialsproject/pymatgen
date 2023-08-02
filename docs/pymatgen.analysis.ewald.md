---
layout: default
title: pymatgen.analysis.ewald.md
nav_exclude: true
---

# pymatgen.analysis.ewald module

This module provides classes for calculating the Ewald sum of a structure.


### _class_ pymatgen.analysis.ewald.EwaldMinimizer(matrix, m_list, num_to_return=1, algo=0)
Bases: `object`

This class determines the manipulations that will minimize an Ewald matrix,
given a list of possible manipulations. This class does not perform the
manipulations on a structure, but will return the list of manipulations
that should be done on one to produce the minimal structure. It returns the
manipulations for the n lowest energy orderings. This class should be used
to perform fractional species substitution or fractional species removal to
produce a new structure. These manipulations create large numbers of
candidate structures, and this class can be used to pick out those with the
lowest Ewald sum.

An alternative (possibly more intuitive) interface to this class is the
order disordered structure transformation.

Author - Will Richards


* **Parameters**


    * **matrix** – A matrix of the Ewald sum interaction energies. This is stored
    in the class as a diagonally symmetric array and so
    self._matrix will not be the same as the input matrix.


    * **m_list** – list of manipulations. each item is of the form
    (multiplication fraction, number_of_indices, indices, species)
    These are sorted such that the first manipulation contains the
    most permutations. this is actually evaluated last in the
    recursion since I’m using pop.


    * **num_to_return** – The minimizer will find the number_returned lowest
    energy structures. This is likely to return a number of duplicate
    structures so it may be necessary to overestimate and then
    remove the duplicates later. (duplicate checking in this
    process is extremely expensive).



#### ALGO_BEST_FIRST(_ = _ )
Slowly increases the speed (with the cost of decreasing
accuracy) as the minimizer runs. Attempts to limit the run time to
approximately 30 minutes.


* **Type**

    ALGO_TIME_LIMIT



#### ALGO_COMPLETE(_ = _ )

#### ALGO_FAST(_ = _ )

#### ALGO_TIME_LIMIT(_ = _ )

#### add_m_list(matrix_sum, m_list)
This adds an m_list to the output_lists and updates the current
minimum if the list is full.


#### best_case(matrix, m_list, indices_left)
Computes a best case given a matrix and manipulation list.


* **Parameters**


    * **matrix** – the current matrix (with some permutations already
    performed)


    * **m_list** – [(multiplication fraction, number_of_indices, indices,
    species)] describing the manipulation


    * **indices** – Set of indices which haven’t had a permutation
    performed on them.



#### _property_ best_m_list()
Best m_list found.


* **Type**

    Returns



#### _classmethod_ get_next_index(matrix, manipulation, indices_left)
Returns an index that should have the most negative effect on the
matrix sum.


#### minimize_matrix()
This method finds and returns the permutations that produce the lowest
Ewald sum calls recursive function to iterate through permutations.


#### _property_ minimized_sum()
Minimized sum.


* **Type**

    Returns



#### _property_ output_lists()
output lists.


* **Type**

    Returns



### _class_ pymatgen.analysis.ewald.EwaldSummation(structure, real_space_cut=None, recip_space_cut=None, eta=None, acc_factor=12.0, w=0.7071067811865475, compute_forces=False)
Bases: `MSONable`

Calculates the electrostatic energy of a periodic array of charges using
the Ewald technique.

Ref:

    Ewald summation techniques in perspective: a survey
    Abdulnour Y. Toukmaji and John A. Board Jr.
    DOI: 10.1016/0010-4655(96)00016-1
    URL: [http://www.ee.duke.edu/~ayt/ewaldpaper/ewaldpaper.html](http://www.ee.duke.edu/~ayt/ewaldpaper/ewaldpaper.html)

This matrix can be used to do fast calculations of Ewald sums after species
removal.

E = E_recip + E_real + E_point

Atomic units used in the code, then converted to eV.

Initializes and calculates the Ewald sum. Default convergence
parameters have been specified, but you can override them if you wish.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure that must have proper
    Species on all sites, i.e. Element with oxidation state. Use
    Structure.add_oxidation_state… for example.


    * **real_space_cut** (*float*) – Real space cutoff radius dictating how
    many terms are used in the real space sum. Defaults to None,
    which means determine automagically using the formula given
    in gulp 3.1 documentation.


    * **recip_space_cut** (*float*) – Reciprocal space cutoff radius.
    Defaults to None, which means determine automagically using
    the formula given in gulp 3.1 documentation.


    * **eta** (*float*) – The screening parameter. Defaults to None, which means
    determine automatically.


    * **acc_factor** (*float*) – No. of significant figures each sum is
    converged to.


    * **w** (*float*) – Weight parameter, w, has been included that represents
    the relative computational expense of calculating a term in
    real and reciprocal space. Default of 0.7 reproduces result
    similar to GULP 4.2. This has little effect on the total
    energy, but may influence speed of computation in large
    systems. Note that this parameter is used only when the
    cutoffs are set to None.


    * **compute_forces** (*bool*) – Whether to compute forces. False by
    default since it is usually not needed.



#### CONV_FACT(_ = 14.3996454784256_ )

#### as_dict(verbosity: int = 0)
Json-serialization dict representation of EwaldSummation.


* **Parameters**

    **verbosity** (*int*) – Verbosity level. Default of 0 only includes the
    matrix representation. Set to 1 for more details.



#### compute_partial_energy(removed_indices)
Gives total Ewald energy for certain sites being removed, i.e. zeroed
out.


#### compute_sub_structure(sub_structure, tol: float = 0.001)
Gives total Ewald energy for an sub structure in the same
lattice. The sub_structure must be a subset of the original
structure, with possible different charges.


* **Parameters**


    * **substructure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Substructure to compute Ewald sum for.


    * **tol** (*float*) – Tolerance for site matching in fractional coordinates.



* **Returns**

    Ewald sum of substructure.



#### _property_ eta()
eta value used in Ewald summation.


* **Type**

    Returns



#### _property_ forces()
The forces on each site as a Nx3 matrix. Each row corresponds to a
site.


#### _classmethod_ from_dict(d: dict[str, Any], fmt: str | None = None, \*\*kwargs)
Create an EwaldSummation instance from JSON-serialized dictionary.


* **Parameters**


    * **d** (*dict*) – Dictionary representation


    * **fmt** (*str**, **optional*) – Unused. Defaults to None.



* **Returns**

    class instance



* **Return type**

    EwaldSummation



#### get_site_energy(site_index)
Compute the energy for a single site in the structure.


* **Parameters**

    **site_index** (*int*) – Index of site


ReturnS:
(float) - Energy of that site


#### _property_ point_energy()
The point energy.


#### _property_ point_energy_matrix()
The point space matrix. A diagonal matrix with the point terms for each
site in the diagonal elements.


#### _property_ real_space_energy()
The real space energy.


#### _property_ real_space_energy_matrix()
The real space energy matrix. Each matrix element (i, j) corresponds to
the interaction energy between site i and site j in real space.


#### _property_ reciprocal_space_energy()
The reciprocal space energy.


#### _property_ reciprocal_space_energy_matrix()
The reciprocal space energy matrix. Each matrix element (i, j)
corresponds to the interaction energy between site i and site j in
reciprocal space.


#### _property_ total_energy()
The total energy.


#### _property_ total_energy_matrix()
The total energy matrix. Each matrix element (i, j) corresponds to the
total interaction energy between site i and site j.

Note that this does not include the charged-cell energy, which is only important
when the simulation cell is not charge balanced.


### pymatgen.analysis.ewald.compute_average_oxidation_state(site)
Calculates the average oxidation state of a site.


* **Parameters**

    **site** – Site to compute average oxidation state



* **Returns**

    Average oxidation state of site.