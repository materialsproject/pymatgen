---
layout: default
title: pymatgen.apps.battery.md
nav_exclude: true
---

# pymatgen.apps.battery package

This package contains all battery-related application classes, including
representations of InsertionElectrodes and ConversionElectrodes.



* [pymatgen.apps.battery.analyzer module](pymatgen.apps.battery.analyzer.md)


    * [`BatteryAnalyzer`](pymatgen.apps.battery.analyzer.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer)


        * [`BatteryAnalyzer.get_max_capgrav()`](pymatgen.apps.battery.analyzer.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.get_max_capgrav)


        * [`BatteryAnalyzer.get_max_capvol()`](pymatgen.apps.battery.analyzer.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.get_max_capvol)


        * [`BatteryAnalyzer.get_removals_int_oxid()`](pymatgen.apps.battery.analyzer.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.get_removals_int_oxid)


        * [`BatteryAnalyzer.max_ion_insertion`](pymatgen.apps.battery.analyzer.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.max_ion_insertion)


        * [`BatteryAnalyzer.max_ion_removal`](pymatgen.apps.battery.analyzer.md#pymatgen.apps.battery.analyzer.BatteryAnalyzer.max_ion_removal)


    * [`is_redox_active_intercalation()`](pymatgen.apps.battery.analyzer.md#pymatgen.apps.battery.analyzer.is_redox_active_intercalation)


* [pymatgen.apps.battery.battery_abc module](pymatgen.apps.battery.battery_abc.md)


    * [`AbstractElectrode`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode)


        * [`AbstractElectrode.voltage_pairs`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.voltage_pairs)


        * [`AbstractElectrode.working_ion`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.working_ion)


        * [`AbstractElectrode.working_ion_entry`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.working_ion_entry)


        * [`AbstractElectrode.framework_formula`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.framework_formula)


        * [`AbstractElectrode.framework`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.framework)


        * [`AbstractElectrode.framework_formula`](pymatgen.apps.battery.battery_abc.md#id0)


        * [`AbstractElectrode.get_average_voltage()`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_average_voltage)


        * [`AbstractElectrode.get_capacity_grav()`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_capacity_grav)


        * [`AbstractElectrode.get_capacity_vol()`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_capacity_vol)


        * [`AbstractElectrode.get_energy_density()`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_energy_density)


        * [`AbstractElectrode.get_specific_energy()`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_specific_energy)


        * [`AbstractElectrode.get_sub_electrodes()`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_sub_electrodes)


        * [`AbstractElectrode.get_summary_dict()`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.get_summary_dict)


        * [`AbstractElectrode.max_delta_volume`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.max_delta_volume)


        * [`AbstractElectrode.max_voltage`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.max_voltage)


        * [`AbstractElectrode.max_voltage_step`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.max_voltage_step)


        * [`AbstractElectrode.min_voltage`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.min_voltage)


        * [`AbstractElectrode.normalization_mass`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.normalization_mass)


        * [`AbstractElectrode.normalization_volume`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.normalization_volume)


        * [`AbstractElectrode.num_steps`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.num_steps)


        * [`AbstractElectrode.voltage_pairs`](pymatgen.apps.battery.battery_abc.md#id1)


        * [`AbstractElectrode.working_ion`](pymatgen.apps.battery.battery_abc.md#id2)


        * [`AbstractElectrode.working_ion_entry`](pymatgen.apps.battery.battery_abc.md#id3)


        * [`AbstractElectrode.x_charge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.x_charge)


        * [`AbstractElectrode.x_discharge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode.x_discharge)


    * [`AbstractVoltagePair`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair)


        * [`AbstractVoltagePair.voltage`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.voltage)


        * [`AbstractVoltagePair.mAh`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.mAh)


        * [`AbstractVoltagePair.mass_charge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.mass_charge)


        * [`AbstractVoltagePair.mass_discharge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.mass_discharge)


        * [`AbstractVoltagePair.vol_charge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.vol_charge)


        * [`AbstractVoltagePair.vol_discharge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.vol_discharge)


        * [`AbstractVoltagePair.frac_charge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.frac_charge)


        * [`AbstractVoltagePair.frac_discharge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.frac_discharge)


        * [`AbstractVoltagePair.working_ion_entry`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.working_ion_entry)


        * [`AbstractVoltagePair.framework_formula`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.framework_formula)


        * [`AbstractVoltagePair.frac_charge`](pymatgen.apps.battery.battery_abc.md#id4)


        * [`AbstractVoltagePair.frac_discharge`](pymatgen.apps.battery.battery_abc.md#id5)


        * [`AbstractVoltagePair.framework`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.framework)


        * [`AbstractVoltagePair.framework_formula`](pymatgen.apps.battery.battery_abc.md#id6)


        * [`AbstractVoltagePair.mAh`](pymatgen.apps.battery.battery_abc.md#id7)


        * [`AbstractVoltagePair.mass_charge`](pymatgen.apps.battery.battery_abc.md#id8)


        * [`AbstractVoltagePair.mass_discharge`](pymatgen.apps.battery.battery_abc.md#id9)


        * [`AbstractVoltagePair.vol_charge`](pymatgen.apps.battery.battery_abc.md#id10)


        * [`AbstractVoltagePair.vol_discharge`](pymatgen.apps.battery.battery_abc.md#id11)


        * [`AbstractVoltagePair.voltage`](pymatgen.apps.battery.battery_abc.md#id12)


        * [`AbstractVoltagePair.working_ion`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.working_ion)


        * [`AbstractVoltagePair.working_ion_entry`](pymatgen.apps.battery.battery_abc.md#id13)


        * [`AbstractVoltagePair.x_charge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.x_charge)


        * [`AbstractVoltagePair.x_discharge`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair.x_discharge)


* [pymatgen.apps.battery.conversion_battery module](pymatgen.apps.battery.conversion_battery.md)


    * [`ConversionElectrode`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode)


        * [`ConversionElectrode.from_composition_and_entries()`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.from_composition_and_entries)


        * [`ConversionElectrode.from_composition_and_pd()`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.from_composition_and_pd)


        * [`ConversionElectrode.get_sub_electrodes()`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.get_sub_electrodes)


        * [`ConversionElectrode.get_summary_dict()`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.get_summary_dict)


        * [`ConversionElectrode.initial_comp`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.initial_comp)


        * [`ConversionElectrode.initial_comp_formula`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.initial_comp_formula)


        * [`ConversionElectrode.is_super_electrode()`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionElectrode.is_super_electrode)


    * [`ConversionVoltagePair`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair)


        * [`ConversionVoltagePair.rxn`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.rxn)


        * [`ConversionVoltagePair.voltage`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.voltage)


        * [`ConversionVoltagePair.mAh`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.mAh)


        * [`ConversionVoltagePair.vol_charge`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.vol_charge)


        * [`ConversionVoltagePair.vol_discharge`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.vol_discharge)


        * [`ConversionVoltagePair.mass_charge`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.mass_charge)


        * [`ConversionVoltagePair.mass_discharge`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.mass_discharge)


        * [`ConversionVoltagePair.frac_charge`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.frac_charge)


        * [`ConversionVoltagePair.frac_discharge`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.frac_discharge)


        * [`ConversionVoltagePair.entries_charge`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.entries_charge)


        * [`ConversionVoltagePair.entries_discharge`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.entries_discharge)


        * [`ConversionVoltagePair.working_ion_entry`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.working_ion_entry)


        * [`ConversionVoltagePair.entries_charge`](pymatgen.apps.battery.conversion_battery.md#id0)


        * [`ConversionVoltagePair.entries_discharge`](pymatgen.apps.battery.conversion_battery.md#id1)


        * [`ConversionVoltagePair.from_steps()`](pymatgen.apps.battery.conversion_battery.md#pymatgen.apps.battery.conversion_battery.ConversionVoltagePair.from_steps)


        * [`ConversionVoltagePair.rxn`](pymatgen.apps.battery.conversion_battery.md#id2)


* [pymatgen.apps.battery.insertion_battery module](pymatgen.apps.battery.insertion_battery.md)


    * [`InsertionElectrode`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode)


        * [`InsertionElectrode.as_dict_legacy()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.as_dict_legacy)


        * [`InsertionElectrode.from_dict_legacy()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.from_dict_legacy)


        * [`InsertionElectrode.from_entries()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.from_entries)


        * [`InsertionElectrode.fully_charged_entry`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.fully_charged_entry)


        * [`InsertionElectrode.fully_discharged_entry`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.fully_discharged_entry)


        * [`InsertionElectrode.get_all_entries()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_all_entries)


        * [`InsertionElectrode.get_max_instability()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_max_instability)


        * [`InsertionElectrode.get_max_muO2()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_max_muO2)


        * [`InsertionElectrode.get_min_instability()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_min_instability)


        * [`InsertionElectrode.get_min_muO2()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_min_muO2)


        * [`InsertionElectrode.get_stable_entries()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_stable_entries)


        * [`InsertionElectrode.get_sub_electrodes()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_sub_electrodes)


        * [`InsertionElectrode.get_summary_dict()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_summary_dict)


        * [`InsertionElectrode.get_unstable_entries()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.get_unstable_entries)


        * [`InsertionElectrode.stable_entries`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.stable_entries)


        * [`InsertionElectrode.unstable_entries`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionElectrode.unstable_entries)


    * [`InsertionVoltagePair`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionVoltagePair)


        * [`InsertionVoltagePair.entry_charge`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionVoltagePair.entry_charge)


        * [`InsertionVoltagePair.entry_discharge`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionVoltagePair.entry_discharge)


        * [`InsertionVoltagePair.from_entries()`](pymatgen.apps.battery.insertion_battery.md#pymatgen.apps.battery.insertion_battery.InsertionVoltagePair.from_entries)


* [pymatgen.apps.battery.plotter module](pymatgen.apps.battery.plotter.md)


    * [`VoltageProfilePlotter`](pymatgen.apps.battery.plotter.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter)


        * [`VoltageProfilePlotter.add_electrode()`](pymatgen.apps.battery.plotter.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.add_electrode)


        * [`VoltageProfilePlotter.get_plot()`](pymatgen.apps.battery.plotter.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.get_plot)


        * [`VoltageProfilePlotter.get_plot_data()`](pymatgen.apps.battery.plotter.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.get_plot_data)


        * [`VoltageProfilePlotter.get_plotly_figure()`](pymatgen.apps.battery.plotter.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.get_plotly_figure)


        * [`VoltageProfilePlotter.save()`](pymatgen.apps.battery.plotter.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.save)


        * [`VoltageProfilePlotter.show()`](pymatgen.apps.battery.plotter.md#pymatgen.apps.battery.plotter.VoltageProfilePlotter.show)