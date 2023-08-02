---
layout: default
title: pymatgen.analysis.quasiharmonic.md
nav_exclude: true
---

# pymatgen.analysis.quasiharmonic module

This module implements the Quasi-harmonic Debye approximation that can
be used to compute thermal properties.

See the following papers for more info:

> [https://doi.org/10.1016/j.comphy.2003.12.001](https://doi.org/10.1016/j.comphy.2003.12.001) (2004)
> [https://doi.org/10.1103/PhysRevB.90.174107](https://doi.org/10.1103/PhysRevB.90.174107) (2014)


### _class_ pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox(energies, volumes, structure, t_min=300.0, t_step=100, t_max=300.0, eos='vinet', pressure=0.0, poisson=0.25, use_mie_gruneisen=False, anharmonic_contribution=False)
Bases: `object`

Quasiharmonic approximation.


* **Parameters**


    * **energies** (*list*) – list of DFT energies in eV


    * **volumes** (*list*) – list of volumes in Ang^3


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – pymatgen structure object


    * **t_min** (*float*) – min temperature


    * **t_step** (*float*) – temperature step


    * **t_max** (*float*) – max temperature


    * **eos** (*str*) – equation of state used for fitting the energies and the
    volumes.
    options supported by pymatgen: “quadratic”, “murnaghan”, “birch”,

    > ”birch_murnaghan”, “pourier_tarantola”, “vinet”,
    > “deltafactor”, “numerical_eos”



    * **pressure** (*float*) – in GPa, optional.


    * **poisson** (*float*) – poisson ratio.


    * **use_mie_gruneisen** (*bool*) – whether or not to use the mie-gruneisen
    formulation to compute the gruneisen parameter.
    The default is the slater-gamma formulation.


    * **anharmonic_contribution** (*bool*) – whether or not to consider the anharmonic
    contribution to the Debye temperature. Cannot be used with
    use_mie_gruneisen. Defaults to False.



#### _static_ debye_integral(y)
Debye integral. Eq(5) in  doi.org/10.1016/j.comphy.2003.12.001.


* **Parameters**

    **y** (*float*) – debye temperature/T, upper limit



* **Returns**

    unitless



* **Return type**

    float



#### debye_temperature(volume)
Calculates the debye temperature.
Eq(6) in doi.org/10.1016/j.comphy.2003.12.001. Thanks to Joey.

Eq(6) above is equivalent to Eq(3) in doi.org/10.1103/PhysRevB.37.790
which does not consider anharmonic effects. Eq(20) in the same paper
and Eq(18) in doi.org/10.1016/j.commatsci.2009.12.006 both consider
anharmonic contributions to the Debye temperature through the Gruneisen
parameter at 0K (Gruneisen constant).

The anharmonic contribution is toggled by setting the anharmonic_contribution
to True or False in the QuasiharmonicDebyeApprox constructor.


* **Parameters**

    **volume** (*float*) – in Ang^3



* **Returns**

    debye temperature in K



* **Return type**

    float



#### get_summary_dict()
Returns a dict with a summary of the computed properties.


#### gruneisen_parameter(temperature, volume)
Slater-gamma formulation(the default):

    gruneisen parameter = - d log(theta)/ d log(V) = - (1/6 + 0.5 d log(B)/ d log(V))

        = - (1/6 + 0.5 V/B dB/dV), where dB/dV = d^2E/dV^2 + V \* d^3E/dV^3.

Mie-gruneisen formulation:

    Eq(31) in doi.org/10.1016/j.comphy.2003.12.001
    Eq(7) in Blanco et. al. Joumal of Molecular Structure (Theochem)

    > 368 (1996) 245-255

    Also se J.P. Poirier, Introduction to the Physics of the Earth’s

        Interior, 2nd ed. (Cambridge University Press, Cambridge,
        2000) Eq(3.53)


* **Parameters**


    * **temperature** (*float*) – temperature in K


    * **volume** (*float*) – in Ang^3



* **Returns**

    unitless



* **Return type**

    float



#### optimize_gibbs_free_energy()
Evaluate the Gibbs free energy as a function of V, T and P i.e
G(V, T, P), minimize G(V, T, P) wrt V for each T and store the
optimum values.

Note: The data points for which the equation of state fitting fails

    are skipped.


#### optimizer(temperature)
Evaluate G(V, T, P) at the given temperature(and pressure) and
minimize it wrt V.


1. Compute the  vibrational Helmholtz free energy, A_vib.


2. Compute the Gibbs free energy as a function of volume, temperature

    and pressure, G(V,T,P).


3. Perform an equation of state fit to get the functional form of

    Gibbs free energy:G(V, T, P).


4. Finally G(V, P, T) is minimized with respect to V.


* **Parameters**

    **temperature** (*float*) – temperature in K



* **Returns**

    G_opt(V_opt, T, P) in eV and V_opt in Ang^3.



* **Return type**

    float, float



#### thermal_conductivity(temperature, volume)
Eq(17) in 10.1103/PhysRevB.90.174107.


* **Parameters**


    * **temperature** (*float*) – temperature in K


    * **volume** (*float*) – in Ang^3



* **Returns**

    thermal conductivity in W/K/m



* **Return type**

    float



#### vibrational_free_energy(temperature, volume)
Vibrational Helmholtz free energy, A_vib(V, T).
Eq(4) in doi.org/10.1016/j.comphy.2003.12.001.


* **Parameters**


    * **temperature** (*float*) – temperature in K


    * **volume** (*float*) –



* **Returns**

    vibrational free energy in eV



* **Return type**

    float



#### vibrational_internal_energy(temperature, volume)
Vibrational internal energy, U_vib(V, T).
Eq(4) in doi.org/10.1016/j.comphy.2003.12.001.


* **Parameters**


    * **temperature** (*float*) – temperature in K


    * **volume** (*float*) – in Ang^3



* **Returns**

    vibrational internal energy in eV



* **Return type**

    float