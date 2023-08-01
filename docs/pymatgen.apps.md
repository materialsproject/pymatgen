---
layout: default
title: pymatgen.apps.md
nav_exclude: true
---

# pymatgen.apps package

This package is an overarching package for subpackages for specific materials
applications.

## Subpackages


* [pymatgen.apps.battery package](pymatgen.apps.battery.md)


    * [Subpackages](pymatgen.apps.battery.md#subpackages)


        * [pymatgen.apps.battery.tests package](pymatgen.apps.battery.tests.md)




            * [pymatgen.apps.battery.tests.test_analyzer module](pymatgen.apps.battery.tests.md#module-pymatgen.apps.battery.tests.test_analyzer)


                * [`BatteryAnalyzerTest`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_analyzer.BatteryAnalyzerTest)


                    * [`BatteryAnalyzerTest.load_from_cif()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_analyzer.BatteryAnalyzerTest.load_from_cif)


                    * [`BatteryAnalyzerTest.load_from_internal()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_analyzer.BatteryAnalyzerTest.load_from_internal)


                    * [`BatteryAnalyzerTest.setUp()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_analyzer.BatteryAnalyzerTest.setUp)


                    * [`BatteryAnalyzerTest.test_capacitygrav_calculations()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_analyzer.BatteryAnalyzerTest.test_capacitygrav_calculations)


                    * [`BatteryAnalyzerTest.test_capacityvol_calculations()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_analyzer.BatteryAnalyzerTest.test_capacityvol_calculations)


                    * [`BatteryAnalyzerTest.test_ion_removal()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_analyzer.BatteryAnalyzerTest.test_ion_removal)


                    * [`BatteryAnalyzerTest.test_oxide_check()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_analyzer.BatteryAnalyzerTest.test_oxide_check)


            * [pymatgen.apps.battery.tests.test_conversion_battery module](pymatgen.apps.battery.tests.md#module-pymatgen.apps.battery.tests.test_conversion_battery)


                * [`ConversionElectrodeTest`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_conversion_battery.ConversionElectrodeTest)


                    * [`ConversionElectrodeTest.setUp()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_conversion_battery.ConversionElectrodeTest.setUp)


                    * [`ConversionElectrodeTest.test_composite()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_conversion_battery.ConversionElectrodeTest.test_composite)


                    * [`ConversionElectrodeTest.test_init()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_conversion_battery.ConversionElectrodeTest.test_init)


                    * [`ConversionElectrodeTest.test_summary()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_conversion_battery.ConversionElectrodeTest.test_summary)


            * [pymatgen.apps.battery.tests.test_insertion_battery module](pymatgen.apps.battery.tests.md#module-pymatgen.apps.battery.tests.test_insertion_battery)


                * [`InsertionElectrodeTest`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest)


                    * [`InsertionElectrodeTest.setUp()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.setUp)


                    * [`InsertionElectrodeTest.test_capacities()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_capacities)


                    * [`InsertionElectrodeTest.test_entries()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_entries)


                    * [`InsertionElectrodeTest.test_get_all_entries()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_get_all_entries)


                    * [`InsertionElectrodeTest.test_get_instability()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_get_instability)


                    * [`InsertionElectrodeTest.test_get_muO2()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_get_muO2)


                    * [`InsertionElectrodeTest.test_get_summary_dict()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_get_summary_dict)


                    * [`InsertionElectrodeTest.test_init_no_structure()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_init_no_structure)


                    * [`InsertionElectrodeTest.test_to_from_dict()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_to_from_dict)


                    * [`InsertionElectrodeTest.test_voltage()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_voltage)


                    * [`InsertionElectrodeTest.test_voltage_pair()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest.test_voltage_pair)


            * [pymatgen.apps.battery.tests.test_plotter module](pymatgen.apps.battery.tests.md#module-pymatgen.apps.battery.tests.test_plotter)


                * [`VoltageProfilePlotterTest`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_plotter.VoltageProfilePlotterTest)


                    * [`VoltageProfilePlotterTest.setUp()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_plotter.VoltageProfilePlotterTest.setUp)


                    * [`VoltageProfilePlotterTest.testName()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_plotter.VoltageProfilePlotterTest.testName)


                    * [`VoltageProfilePlotterTest.testPlotly()`](pymatgen.apps.battery.tests.md#pymatgen.apps.battery.tests.test_plotter.VoltageProfilePlotterTest.testPlotly)




    * [pymatgen.apps.battery.analyzer module](pymatgen.apps.battery.md#module-pymatgen.apps.battery.analyzer)


        * [`BatteryAnalyzer`](pymatgen.apps.battery.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer)


            * [`BatteryAnalyzer.get_max_capgrav()`](pymatgen.apps.battery.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.get_max_capgrav)


            * [`BatteryAnalyzer.get_max_capvol()`](pymatgen.apps.battery.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.get_max_capvol)


            * [`BatteryAnalyzer.get_removals_int_oxid()`](pymatgen.apps.battery.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.get_removals_int_oxid)


            * [`BatteryAnalyzer.max_ion_insertion`](pymatgen.apps.battery.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.max_ion_insertion)


            * [`BatteryAnalyzer.max_ion_removal`](pymatgen.apps.battery.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.max_ion_removal)


        * [`is_redox_active_intercalation()`](pymatgen.apps.battery.md#pymatgen.apps.battery.analyzer.is_redox_active_intercalation)


    * [pymatgen.apps.battery.battery_abc module](pymatgen.apps.battery.md#module-pymatgen.apps.battery.battery_abc)


        * [`AbstractElectrode`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode)


            * [`AbstractElectrode.voltage_pairs`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.voltage_pairs)


            * [`AbstractElectrode.working_ion`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.working_ion)


            * [`AbstractElectrode.working_ion_entry`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.working_ion_entry)


            * [`AbstractElectrode.framework_formula`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.framework_formula)


            * [`AbstractElectrode.framework`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.framework)


            * [`AbstractElectrode.framework_formula`](pymatgen.apps.battery.md#id0)


            * [`AbstractElectrode.get_average_voltage()`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_average_voltage)


            * [`AbstractElectrode.get_capacity_grav()`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_capacity_grav)


            * [`AbstractElectrode.get_capacity_vol()`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_capacity_vol)


            * [`AbstractElectrode.get_energy_density()`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_energy_density)


            * [`AbstractElectrode.get_specific_energy()`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_specific_energy)


            * [`AbstractElectrode.get_sub_electrodes()`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_sub_electrodes)


            * [`AbstractElectrode.get_summary_dict()`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_summary_dict)


            * [`AbstractElectrode.max_delta_volume`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.max_delta_volume)


            * [`AbstractElectrode.max_voltage`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.max_voltage)


            * [`AbstractElectrode.max_voltage_step`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.max_voltage_step)


            * [`AbstractElectrode.min_voltage`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.min_voltage)


            * [`AbstractElectrode.normalization_mass`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.normalization_mass)


            * [`AbstractElectrode.normalization_volume`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.normalization_volume)


            * [`AbstractElectrode.num_steps`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.num_steps)


            * [`AbstractElectrode.voltage_pairs`](pymatgen.apps.battery.md#id1)


            * [`AbstractElectrode.working_ion`](pymatgen.apps.battery.md#id2)


            * [`AbstractElectrode.working_ion_entry`](pymatgen.apps.battery.md#id3)


            * [`AbstractElectrode.x_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.x_charge)


            * [`AbstractElectrode.x_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.x_discharge)


        * [`AbstractVoltagePair`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair)


            * [`AbstractVoltagePair.voltage`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.voltage)


            * [`AbstractVoltagePair.mAh`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.mAh)


            * [`AbstractVoltagePair.mass_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.mass_charge)


            * [`AbstractVoltagePair.mass_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.mass_discharge)


            * [`AbstractVoltagePair.vol_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.vol_charge)


            * [`AbstractVoltagePair.vol_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.vol_discharge)


            * [`AbstractVoltagePair.frac_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.frac_charge)


            * [`AbstractVoltagePair.frac_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.frac_discharge)


            * [`AbstractVoltagePair.working_ion_entry`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.working_ion_entry)


            * [`AbstractVoltagePair.framework_formula`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.framework_formula)


            * [`AbstractVoltagePair.frac_charge`](pymatgen.apps.battery.md#id4)


            * [`AbstractVoltagePair.frac_discharge`](pymatgen.apps.battery.md#id5)


            * [`AbstractVoltagePair.framework`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.framework)


            * [`AbstractVoltagePair.framework_formula`](pymatgen.apps.battery.md#id6)


            * [`AbstractVoltagePair.mAh`](pymatgen.apps.battery.md#id7)


            * [`AbstractVoltagePair.mass_charge`](pymatgen.apps.battery.md#id8)


            * [`AbstractVoltagePair.mass_discharge`](pymatgen.apps.battery.md#id9)


            * [`AbstractVoltagePair.vol_charge`](pymatgen.apps.battery.md#id10)


            * [`AbstractVoltagePair.vol_discharge`](pymatgen.apps.battery.md#id11)


            * [`AbstractVoltagePair.voltage`](pymatgen.apps.battery.md#id12)


            * [`AbstractVoltagePair.working_ion`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.working_ion)


            * [`AbstractVoltagePair.working_ion_entry`](pymatgen.apps.battery.md#id13)


            * [`AbstractVoltagePair.x_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.x_charge)


            * [`AbstractVoltagePair.x_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.x_discharge)


    * [pymatgen.apps.battery.conversion_battery module](pymatgen.apps.battery.md#module-pymatgen.apps.battery.conversion_battery)


        * [`ConversionElectrode`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode)


            * [`ConversionElectrode.from_composition_and_entries()`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.from_composition_and_entries)


            * [`ConversionElectrode.from_composition_and_pd()`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.from_composition_and_pd)


            * [`ConversionElectrode.get_sub_electrodes()`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.get_sub_electrodes)


            * [`ConversionElectrode.get_summary_dict()`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.get_summary_dict)


            * [`ConversionElectrode.initial_comp`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.initial_comp)


            * [`ConversionElectrode.initial_comp_formula`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.initial_comp_formula)


            * [`ConversionElectrode.is_super_electrode()`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.is_super_electrode)


        * [`ConversionVoltagePair`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair)


            * [`ConversionVoltagePair.rxn`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.rxn)


            * [`ConversionVoltagePair.voltage`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.voltage)


            * [`ConversionVoltagePair.mAh`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.mAh)


            * [`ConversionVoltagePair.vol_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.vol_charge)


            * [`ConversionVoltagePair.vol_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.vol_discharge)


            * [`ConversionVoltagePair.mass_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.mass_charge)


            * [`ConversionVoltagePair.mass_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.mass_discharge)


            * [`ConversionVoltagePair.frac_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.frac_charge)


            * [`ConversionVoltagePair.frac_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.frac_discharge)


            * [`ConversionVoltagePair.entries_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.entries_charge)


            * [`ConversionVoltagePair.entries_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.entries_discharge)


            * [`ConversionVoltagePair.working_ion_entry`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.working_ion_entry)


            * [`ConversionVoltagePair.entries_charge`](pymatgen.apps.battery.md#id14)


            * [`ConversionVoltagePair.entries_discharge`](pymatgen.apps.battery.md#id15)


            * [`ConversionVoltagePair.from_steps()`](pymatgen.apps.battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.from_steps)


            * [`ConversionVoltagePair.rxn`](pymatgen.apps.battery.md#id16)


    * [pymatgen.apps.battery.insertion_battery module](pymatgen.apps.battery.md#module-pymatgen.apps.battery.insertion_battery)


        * [`InsertionElectrode`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode)


            * [`InsertionElectrode.as_dict_legacy()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.as_dict_legacy)


            * [`InsertionElectrode.from_dict_legacy()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.from_dict_legacy)


            * [`InsertionElectrode.from_entries()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.from_entries)


            * [`InsertionElectrode.fully_charged_entry`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.fully_charged_entry)


            * [`InsertionElectrode.fully_discharged_entry`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.fully_discharged_entry)


            * [`InsertionElectrode.get_all_entries()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_all_entries)


            * [`InsertionElectrode.get_max_instability()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_max_instability)


            * [`InsertionElectrode.get_max_muO2()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_max_muO2)


            * [`InsertionElectrode.get_min_instability()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_min_instability)


            * [`InsertionElectrode.get_min_muO2()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_min_muO2)


            * [`InsertionElectrode.get_stable_entries()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_stable_entries)


            * [`InsertionElectrode.get_sub_electrodes()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_sub_electrodes)


            * [`InsertionElectrode.get_summary_dict()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_summary_dict)


            * [`InsertionElectrode.get_unstable_entries()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_unstable_entries)


            * [`InsertionElectrode.stable_entries`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.stable_entries)


            * [`InsertionElectrode.unstable_entries`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.unstable_entries)


        * [`InsertionVoltagePair`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionVoltagePair)


            * [`InsertionVoltagePair.entry_charge`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionVoltagePair.entry_charge)


            * [`InsertionVoltagePair.entry_discharge`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionVoltagePair.entry_discharge)


            * [`InsertionVoltagePair.from_entries()`](pymatgen.apps.battery.md#pymatgen.apps.battery.insertion_battery.InsertionVoltagePair.from_entries)


    * [pymatgen.apps.battery.plotter module](pymatgen.apps.battery.md#module-pymatgen.apps.battery.plotter)


        * [`VoltageProfilePlotter`](pymatgen.apps.battery.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter)


            * [`VoltageProfilePlotter.add_electrode()`](pymatgen.apps.battery.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.add_electrode)


            * [`VoltageProfilePlotter.get_plot()`](pymatgen.apps.battery.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.get_plot)


            * [`VoltageProfilePlotter.get_plot_data()`](pymatgen.apps.battery.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.get_plot_data)


            * [`VoltageProfilePlotter.get_plotly_figure()`](pymatgen.apps.battery.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.get_plotly_figure)


            * [`VoltageProfilePlotter.save()`](pymatgen.apps.battery.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.save)


            * [`VoltageProfilePlotter.show()`](pymatgen.apps.battery.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.show)


* [pymatgen.apps.borg package](pymatgen.apps.borg.md)


    * [Subpackages](pymatgen.apps.borg.md#subpackages)


        * [pymatgen.apps.borg.tests package](pymatgen.apps.borg.tests.md)




            * [pymatgen.apps.borg.tests.test_hive module](pymatgen.apps.borg.tests.md#module-pymatgen.apps.borg.tests.test_hive)


                * [`GaussianToComputedEntryDroneTest`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.GaussianToComputedEntryDroneTest)


                    * [`GaussianToComputedEntryDroneTest.setUp()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.GaussianToComputedEntryDroneTest.setUp)


                    * [`GaussianToComputedEntryDroneTest.tearDown()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.GaussianToComputedEntryDroneTest.tearDown)


                    * [`GaussianToComputedEntryDroneTest.test_assimilate()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.GaussianToComputedEntryDroneTest.test_assimilate)


                    * [`GaussianToComputedEntryDroneTest.test_get_valid_paths()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.GaussianToComputedEntryDroneTest.test_get_valid_paths)


                    * [`GaussianToComputedEntryDroneTest.test_to_from_dict()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.GaussianToComputedEntryDroneTest.test_to_from_dict)


                * [`SimpleVaspToComputedEntryDroneTest`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.SimpleVaspToComputedEntryDroneTest)


                    * [`SimpleVaspToComputedEntryDroneTest.setUp()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.SimpleVaspToComputedEntryDroneTest.setUp)


                    * [`SimpleVaspToComputedEntryDroneTest.tearDown()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.SimpleVaspToComputedEntryDroneTest.tearDown)


                    * [`SimpleVaspToComputedEntryDroneTest.test_get_valid_paths()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.SimpleVaspToComputedEntryDroneTest.test_get_valid_paths)


                    * [`SimpleVaspToComputedEntryDroneTest.test_to_from_dict()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.SimpleVaspToComputedEntryDroneTest.test_to_from_dict)


                * [`VaspToComputedEntryDroneTest`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.VaspToComputedEntryDroneTest)


                    * [`VaspToComputedEntryDroneTest.setUp()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.VaspToComputedEntryDroneTest.setUp)


                    * [`VaspToComputedEntryDroneTest.tearDown()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.VaspToComputedEntryDroneTest.tearDown)


                    * [`VaspToComputedEntryDroneTest.test_assimilate()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.VaspToComputedEntryDroneTest.test_assimilate)


                    * [`VaspToComputedEntryDroneTest.test_get_valid_paths()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.VaspToComputedEntryDroneTest.test_get_valid_paths)


                    * [`VaspToComputedEntryDroneTest.test_to_from_dict()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_hive.VaspToComputedEntryDroneTest.test_to_from_dict)


            * [pymatgen.apps.borg.tests.test_queen module](pymatgen.apps.borg.tests.md#module-pymatgen.apps.borg.tests.test_queen)


                * [`BorgQueenTest`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_queen.BorgQueenTest)


                    * [`BorgQueenTest.setUp()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_queen.BorgQueenTest.setUp)


                    * [`BorgQueenTest.tearDown()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_queen.BorgQueenTest.tearDown)


                    * [`BorgQueenTest.test_get_data()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_queen.BorgQueenTest.test_get_data)


                    * [`BorgQueenTest.test_load_data()`](pymatgen.apps.borg.tests.md#pymatgen.apps.borg.tests.test_queen.BorgQueenTest.test_load_data)




    * [pymatgen.apps.borg.hive module](pymatgen.apps.borg.md#module-pymatgen.apps.borg.hive)


        * [`AbstractDrone`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.AbstractDrone)


            * [`AbstractDrone.assimilate()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.AbstractDrone.assimilate)


            * [`AbstractDrone.get_valid_paths()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.AbstractDrone.get_valid_paths)


        * [`GaussianToComputedEntryDrone`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone)


            * [`GaussianToComputedEntryDrone.as_dict()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone.as_dict)


            * [`GaussianToComputedEntryDrone.assimilate()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone.assimilate)


            * [`GaussianToComputedEntryDrone.from_dict()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone.from_dict)


            * [`GaussianToComputedEntryDrone.get_valid_paths()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.GaussianToComputedEntryDrone.get_valid_paths)


        * [`SimpleVaspToComputedEntryDrone`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.SimpleVaspToComputedEntryDrone)


            * [`SimpleVaspToComputedEntryDrone.as_dict()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.SimpleVaspToComputedEntryDrone.as_dict)


            * [`SimpleVaspToComputedEntryDrone.assimilate()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.SimpleVaspToComputedEntryDrone.assimilate)


            * [`SimpleVaspToComputedEntryDrone.from_dict()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.SimpleVaspToComputedEntryDrone.from_dict)


        * [`VaspToComputedEntryDrone`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone)


            * [`VaspToComputedEntryDrone.as_dict()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone.as_dict)


            * [`VaspToComputedEntryDrone.assimilate()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone.assimilate)


            * [`VaspToComputedEntryDrone.from_dict()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone.from_dict)


            * [`VaspToComputedEntryDrone.get_valid_paths()`](pymatgen.apps.borg.md#pymatgen.apps.borg.hive.VaspToComputedEntryDrone.get_valid_paths)


    * [pymatgen.apps.borg.queen module](pymatgen.apps.borg.md#module-pymatgen.apps.borg.queen)


        * [`BorgQueen`](pymatgen.apps.borg.md#pymatgen.apps.borg.queen.BorgQueen)


            * [`BorgQueen.get_data()`](pymatgen.apps.borg.md#pymatgen.apps.borg.queen.BorgQueen.get_data)


            * [`BorgQueen.load_data()`](pymatgen.apps.borg.md#pymatgen.apps.borg.queen.BorgQueen.load_data)


            * [`BorgQueen.parallel_assimilate()`](pymatgen.apps.borg.md#pymatgen.apps.borg.queen.BorgQueen.parallel_assimilate)


            * [`BorgQueen.save_data()`](pymatgen.apps.borg.md#pymatgen.apps.borg.queen.BorgQueen.save_data)


            * [`BorgQueen.serial_assimilate()`](pymatgen.apps.borg.md#pymatgen.apps.borg.queen.BorgQueen.serial_assimilate)


        * [`order_assimilation()`](pymatgen.apps.borg.md#pymatgen.apps.borg.queen.order_assimilation)