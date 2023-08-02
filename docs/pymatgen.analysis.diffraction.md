---
layout: default
title: pymatgen.analysis.diffraction.md
nav_exclude: true
---

# pymatgen.analysis.diffraction package

This package implements various diffraction analyses.



* [pymatgen.analysis.diffraction.core module](pymatgen.analysis.diffraction.core.md)


    * [`AbstractDiffractionPatternCalculator`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator)


        * [`AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL)


        * [`AbstractDiffractionPatternCalculator.TWO_THETA_TOL`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.TWO_THETA_TOL)


        * [`AbstractDiffractionPatternCalculator.get_pattern()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.get_pattern)


        * [`AbstractDiffractionPatternCalculator.get_plot()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.get_plot)


        * [`AbstractDiffractionPatternCalculator.plot_structures()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.plot_structures)


        * [`AbstractDiffractionPatternCalculator.show_plot()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.show_plot)


    * [`DiffractionPattern`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.DiffractionPattern)


        * [`DiffractionPattern.XLABEL`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.DiffractionPattern.XLABEL)


        * [`DiffractionPattern.YLABEL`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.DiffractionPattern.YLABEL)


    * [`get_unique_families()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.get_unique_families)


* [pymatgen.analysis.diffraction.neutron module](pymatgen.analysis.diffraction.neutron.md)


    * [`NDCalculator`](pymatgen.analysis.diffraction.neutron.md#pymatgen.analysis.diffraction.neutron.NDCalculator)


        * [`NDCalculator.get_pattern()`](pymatgen.analysis.diffraction.neutron.md#pymatgen.analysis.diffraction.neutron.NDCalculator.get_pattern)


* [pymatgen.analysis.diffraction.tem module](pymatgen.analysis.diffraction.tem.md)


    * [`TEMCalculator`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator)


        * [`TEMCalculator.bragg_angles()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.bragg_angles)


        * [`TEMCalculator.cell_intensity()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.cell_intensity)


        * [`TEMCalculator.cell_scattering_factors()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.cell_scattering_factors)


        * [`TEMCalculator.electron_scattering_factors()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.electron_scattering_factors)


        * [`TEMCalculator.generate_points()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.generate_points)


        * [`TEMCalculator.get_first_point()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_first_point)


        * [`TEMCalculator.get_interplanar_angle()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_interplanar_angle)


        * [`TEMCalculator.get_interplanar_spacings()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_interplanar_spacings)


        * [`TEMCalculator.get_pattern()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_pattern)


        * [`TEMCalculator.get_plot_2d()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_plot_2d)


        * [`TEMCalculator.get_plot_2d_concise()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_plot_2d_concise)


        * [`TEMCalculator.get_plot_coeffs()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_plot_coeffs)


        * [`TEMCalculator.get_positions()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_positions)


        * [`TEMCalculator.get_s2()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_s2)


        * [`TEMCalculator.is_parallel()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.is_parallel)


        * [`TEMCalculator.normalized_cell_intensity()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.normalized_cell_intensity)


        * [`TEMCalculator.tem_dots()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.tem_dots)


        * [`TEMCalculator.wavelength_rel()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.wavelength_rel)


        * [`TEMCalculator.x_ray_factors()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.x_ray_factors)


        * [`TEMCalculator.zone_axis_filter()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.zone_axis_filter)


* [pymatgen.analysis.diffraction.xrd module](pymatgen.analysis.diffraction.xrd.md)


    * [`XRDCalculator`](pymatgen.analysis.diffraction.xrd.md#pymatgen.analysis.diffraction.xrd.XRDCalculator)


        * [`XRDCalculator.AVAILABLE_RADIATION`](pymatgen.analysis.diffraction.xrd.md#pymatgen.analysis.diffraction.xrd.XRDCalculator.AVAILABLE_RADIATION)


        * [`XRDCalculator.get_pattern()`](pymatgen.analysis.diffraction.xrd.md#pymatgen.analysis.diffraction.xrd.XRDCalculator.get_pattern)