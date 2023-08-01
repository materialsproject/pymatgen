---
layout: default
title: pymatgen.analysis.diffraction.tests.md
nav_exclude: true
---

# pymatgen.analysis.diffraction.tests package

TODO: Modify module doc.


## pymatgen.analysis.diffraction.tests.test_neutron module


### _class_ pymatgen.analysis.diffraction.tests.test_neutron.NDCalculatorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_pattern()

#### test_get_plot()
## pymatgen.analysis.diffraction.tests.test_tem module

Unit tests for TEM calculator.


### _class_ pymatgen.analysis.diffraction.tests.test_tem.TEMCalculatorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_TEM_dots()

#### test_bragg_angles()

#### test_cell_intensity()

#### test_cell_scattering_factors()

#### test_electron_scattering_factors()

#### test_generate_points()

#### test_get_first_point()

#### test_get_interplanar_spacings()

#### test_get_pattern()

#### test_get_plot_2d()

#### test_get_plot_2d_concise()

#### test_get_plot_coeffs()

#### test_get_positions()

#### test_get_s2()

#### test_interplanar_angle()

#### test_is_parallel()

#### test_normalized_cell_intensity()

#### test_wavelength_rel()

#### test_x_ray_factors()

#### test_zone_axis_filter()
## pymatgen.analysis.diffraction.tests.test_xrd module


### _class_ pymatgen.analysis.diffraction.tests.test_xrd.XRDCalculatorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_pattern()

#### test_type_wavelength()
Test TypeError is raised if wavelength is unaccepted type