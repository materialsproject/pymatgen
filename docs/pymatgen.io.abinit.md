---
layout: default
title: pymatgen.io.abinit.md
nav_exclude: true
---

# pymatgen.io.abinit package

This package implements basic input and output capabilities for Abinit.



* [pymatgen.io.abinit.abiobjects module](pymatgen.io.abinit.abiobjects.md)


    * [`AbivarAble`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.AbivarAble)


        * [`AbivarAble.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.AbivarAble.to_abivars)


    * [`Constraints`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Constraints)


        * [`Constraints.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Constraints.to_abivars)


    * [`Electrons`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Electrons)


        * [`Electrons.as_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Electrons.as_dict)


        * [`Electrons.from_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Electrons.from_dict)


        * [`Electrons.nspden`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Electrons.nspden)


        * [`Electrons.nspinor`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Electrons.nspinor)


        * [`Electrons.nsppol`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Electrons.nsppol)


        * [`Electrons.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Electrons.to_abivars)


    * [`ElectronsAlgorithm`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm)


        * [`ElectronsAlgorithm.as_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm.as_dict)


        * [`ElectronsAlgorithm.from_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm.from_dict)


        * [`ElectronsAlgorithm.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm.to_abivars)


    * [`ExcHamiltonian`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian)


        * [`ExcHamiltonian.inclvkb`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.inclvkb)


        * [`ExcHamiltonian.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.to_abivars)


        * [`ExcHamiltonian.use_cg`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.use_cg)


        * [`ExcHamiltonian.use_direct_diago`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.use_direct_diago)


        * [`ExcHamiltonian.use_haydock`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.use_haydock)


    * [`HilbertTransform`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.HilbertTransform)


        * [`HilbertTransform.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.HilbertTransform.to_abivars)


    * [`KSampling`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling)


        * [`KSampling.as_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.as_dict)


        * [`KSampling.automatic_density()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.automatic_density)


        * [`KSampling.explicit_path()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.explicit_path)


        * [`KSampling.from_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.from_dict)


        * [`KSampling.gamma_centered()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.gamma_centered)


        * [`KSampling.gamma_only()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.gamma_only)


        * [`KSampling.is_homogeneous`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.is_homogeneous)


        * [`KSampling.monkhorst()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.monkhorst)


        * [`KSampling.monkhorst_automatic()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.monkhorst_automatic)


        * [`KSampling.path_from_structure()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.path_from_structure)


        * [`KSampling.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSampling.to_abivars)


    * [`KSamplingModes`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSamplingModes)


        * [`KSamplingModes.automatic`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSamplingModes.automatic)


        * [`KSamplingModes.monkhorst`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSamplingModes.monkhorst)


        * [`KSamplingModes.path`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.KSamplingModes.path)


    * [`ModelDielectricFunction`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ModelDielectricFunction)


        * [`ModelDielectricFunction.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.ModelDielectricFunction.to_abivars)


    * [`PPModel`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModel)


        * [`PPModel.as_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModel.as_dict)


        * [`PPModel.as_ppmodel()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModel.as_ppmodel)


        * [`PPModel.from_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModel.from_dict)


        * [`PPModel.get_noppmodel()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModel.get_noppmodel)


        * [`PPModel.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModel.to_abivars)


    * [`PPModelModes`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModelModes)


        * [`PPModelModes.farid`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModelModes.farid)


        * [`PPModelModes.godby`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModelModes.godby)


        * [`PPModelModes.hybersten`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModelModes.hybersten)


        * [`PPModelModes.linden`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModelModes.linden)


        * [`PPModelModes.noppmodel`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.PPModelModes.noppmodel)


    * [`RelaxationMethod`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod)


        * [`RelaxationMethod.IONMOV_DEFAULT`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.IONMOV_DEFAULT)


        * [`RelaxationMethod.OPTCELL_DEFAULT`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.OPTCELL_DEFAULT)


        * [`RelaxationMethod.as_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.as_dict)


        * [`RelaxationMethod.atoms_and_cell()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.atoms_and_cell)


        * [`RelaxationMethod.atoms_only()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.atoms_only)


        * [`RelaxationMethod.from_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.from_dict)


        * [`RelaxationMethod.move_atoms`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.move_atoms)


        * [`RelaxationMethod.move_cell`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.move_cell)


        * [`RelaxationMethod.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.to_abivars)


    * [`Screening`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Screening)


        * [`Screening.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Screening.to_abivars)


        * [`Screening.use_hilbert`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Screening.use_hilbert)


    * [`SelfEnergy`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SelfEnergy)


        * [`SelfEnergy.gwcalctyp`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SelfEnergy.gwcalctyp)


        * [`SelfEnergy.symsigma`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SelfEnergy.symsigma)


        * [`SelfEnergy.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SelfEnergy.to_abivars)


        * [`SelfEnergy.use_ppmodel`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SelfEnergy.use_ppmodel)


    * [`Smearing`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Smearing)


        * [`Smearing.as_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Smearing.as_dict)


        * [`Smearing.as_smearing()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Smearing.as_smearing)


        * [`Smearing.from_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Smearing.from_dict)


        * [`Smearing.mode`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Smearing.mode)


        * [`Smearing.nosmearing()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Smearing.nosmearing)


        * [`Smearing.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.Smearing.to_abivars)


    * [`SpinMode`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SpinMode)


        * [`SpinMode.as_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SpinMode.as_dict)


        * [`SpinMode.as_spinmode()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SpinMode.as_spinmode)


        * [`SpinMode.from_dict()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SpinMode.from_dict)


        * [`SpinMode.to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.SpinMode.to_abivars)


    * [`contract()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.contract)


    * [`lattice_from_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.lattice_from_abivars)


    * [`species_by_znucl()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.species_by_znucl)


    * [`structure_from_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.structure_from_abivars)


    * [`structure_to_abivars()`](pymatgen.io.abinit.abiobjects.md#pymatgen.io.abinit.abiobjects.structure_to_abivars)


* [pymatgen.io.abinit.abitimer module](pymatgen.io.abinit.abitimer.md)


    * [`AbinitTimer`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer)


        * [`AbinitTimer.cpuwall_histogram()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.cpuwall_histogram)


        * [`AbinitTimer.get_dataframe()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.get_dataframe)


        * [`AbinitTimer.get_section()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.get_section)


        * [`AbinitTimer.get_values()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.get_values)


        * [`AbinitTimer.names_and_values()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.names_and_values)


        * [`AbinitTimer.ncpus`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.ncpus)


        * [`AbinitTimer.order_sections()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.order_sections)


        * [`AbinitTimer.pie()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.pie)


        * [`AbinitTimer.scatter_hist()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.scatter_hist)


        * [`AbinitTimer.sum_sections()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.sum_sections)


        * [`AbinitTimer.to_csv()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.to_csv)


        * [`AbinitTimer.to_table()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.to_table)


        * [`AbinitTimer.totable()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimer.totable)


    * [`AbinitTimerParseError`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParseError)


    * [`AbinitTimerParser`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser)


        * [`AbinitTimerParser.BEGIN_TAG`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.BEGIN_TAG)


        * [`AbinitTimerParser.END_TAG`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.END_TAG)


        * [`AbinitTimerParser.Error`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.Error)


        * [`AbinitTimerParser.filenames`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.filenames)


        * [`AbinitTimerParser.get_sections()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.get_sections)


        * [`AbinitTimerParser.parse()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.parse)


        * [`AbinitTimerParser.pefficiency()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.pefficiency)


        * [`AbinitTimerParser.plot_all()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.plot_all)


        * [`AbinitTimerParser.plot_efficiency()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.plot_efficiency)


        * [`AbinitTimerParser.plot_pie()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.plot_pie)


        * [`AbinitTimerParser.plot_stacked_hist()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.plot_stacked_hist)


        * [`AbinitTimerParser.section_names()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.section_names)


        * [`AbinitTimerParser.summarize()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.summarize)


        * [`AbinitTimerParser.timers()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.timers)


        * [`AbinitTimerParser.walk()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.walk)


    * [`AbinitTimerSection`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerSection)


        * [`AbinitTimerSection.FIELDS`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.FIELDS)


        * [`AbinitTimerSection.NUMERIC_FIELDS`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.NUMERIC_FIELDS)


        * [`AbinitTimerSection.STR_FIELDS`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.STR_FIELDS)


        * [`AbinitTimerSection.fake()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.fake)


        * [`AbinitTimerSection.to_csvline()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.to_csvline)


        * [`AbinitTimerSection.to_dict()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.to_dict)


        * [`AbinitTimerSection.to_tuple()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.to_tuple)


    * [`ParallelEfficiency`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.ParallelEfficiency)


        * [`ParallelEfficiency.bad_sections()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.ParallelEfficiency.bad_sections)


        * [`ParallelEfficiency.good_sections()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.ParallelEfficiency.good_sections)


        * [`ParallelEfficiency.totable()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.ParallelEfficiency.totable)


    * [`alternate()`](pymatgen.io.abinit.abitimer.md#pymatgen.io.abinit.abitimer.alternate)


* [pymatgen.io.abinit.inputs module](pymatgen.io.abinit.inputs.md)


    * [`AbstractInput`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.AbstractInput)


        * [`AbstractInput.deepcopy()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.AbstractInput.deepcopy)


        * [`AbstractInput.pop_vars()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.AbstractInput.pop_vars)


        * [`AbstractInput.remove_vars()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.AbstractInput.remove_vars)


        * [`AbstractInput.set_vars()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.AbstractInput.set_vars)


        * [`AbstractInput.set_vars_ifnotin()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.AbstractInput.set_vars_ifnotin)


        * [`AbstractInput.to_str()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.AbstractInput.to_str)


        * [`AbstractInput.vars`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.AbstractInput.vars)


        * [`AbstractInput.write()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.AbstractInput.write)


    * [`BasicAbinitInput`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput)


        * [`BasicAbinitInput.Error`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.Error)


        * [`BasicAbinitInput.add_abiobjects()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.add_abiobjects)


        * [`BasicAbinitInput.as_dict()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.as_dict)


        * [`BasicAbinitInput.comment`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.comment)


        * [`BasicAbinitInput.from_dict()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.from_dict)


        * [`BasicAbinitInput.isnc`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.isnc)


        * [`BasicAbinitInput.ispaw`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.ispaw)


        * [`BasicAbinitInput.new_with_vars()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.new_with_vars)


        * [`BasicAbinitInput.pop_irdvars()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.pop_irdvars)


        * [`BasicAbinitInput.pop_tolerances()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.pop_tolerances)


        * [`BasicAbinitInput.pseudos`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.pseudos)


        * [`BasicAbinitInput.set_comment()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_comment)


        * [`BasicAbinitInput.set_gamma_sampling()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_gamma_sampling)


        * [`BasicAbinitInput.set_kmesh()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_kmesh)


        * [`BasicAbinitInput.set_kpath()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_kpath)


        * [`BasicAbinitInput.set_spin_mode()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_spin_mode)


        * [`BasicAbinitInput.set_structure()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_structure)


        * [`BasicAbinitInput.structure`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.structure)


        * [`BasicAbinitInput.to_str()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.to_str)


        * [`BasicAbinitInput.to_string()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.to_string)


        * [`BasicAbinitInput.vars`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInput.vars)


    * [`BasicAbinitInputError`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicAbinitInputError)


    * [`BasicMultiDataset`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset)


        * [`BasicMultiDataset.Error`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.Error)


        * [`BasicMultiDataset.addnew_from()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.addnew_from)


        * [`BasicMultiDataset.append()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.append)


        * [`BasicMultiDataset.deepcopy()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.deepcopy)


        * [`BasicMultiDataset.extend()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.extend)


        * [`BasicMultiDataset.from_inputs()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.from_inputs)


        * [`BasicMultiDataset.has_same_structures`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.has_same_structures)


        * [`BasicMultiDataset.isnc`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.isnc)


        * [`BasicMultiDataset.ispaw`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.ispaw)


        * [`BasicMultiDataset.ndtset`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.ndtset)


        * [`BasicMultiDataset.pseudos`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.pseudos)


        * [`BasicMultiDataset.replicate_input()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.replicate_input)


        * [`BasicMultiDataset.split_datasets()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.split_datasets)


        * [`BasicMultiDataset.to_str()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.to_str)


        * [`BasicMultiDataset.to_string()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.to_string)


        * [`BasicMultiDataset.write()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.BasicMultiDataset.write)


    * [`ShiftMode`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.ShiftMode)


        * [`ShiftMode.GammaCentered`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.ShiftMode.GammaCentered)


        * [`ShiftMode.MonkhorstPack`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.ShiftMode.MonkhorstPack)


        * [`ShiftMode.OneSymmetric`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.ShiftMode.OneSymmetric)


        * [`ShiftMode.Symmetric`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.ShiftMode.Symmetric)


        * [`ShiftMode.from_object()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.ShiftMode.from_object)


    * [`as_structure()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.as_structure)


    * [`calc_shiftk()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.calc_shiftk)


    * [`ebands_input()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.ebands_input)


    * [`gs_input()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.gs_input)


    * [`ion_ioncell_relax_input()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.ion_ioncell_relax_input)


    * [`num_valence_electrons()`](pymatgen.io.abinit.inputs.md#pymatgen.io.abinit.inputs.num_valence_electrons)


* [pymatgen.io.abinit.netcdf module](pymatgen.io.abinit.netcdf.md)


    * [`ETSF_Reader`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.ETSF_Reader)


        * [`ETSF_Reader.chemical_symbols()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.ETSF_Reader.chemical_symbols)


        * [`ETSF_Reader.read_abinit_hdr()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.ETSF_Reader.read_abinit_hdr)


        * [`ETSF_Reader.read_abinit_xcfunc()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.ETSF_Reader.read_abinit_xcfunc)


        * [`ETSF_Reader.read_structure()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.ETSF_Reader.read_structure)


        * [`ETSF_Reader.typeidx_from_symbol()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.ETSF_Reader.typeidx_from_symbol)


    * [`NO_DEFAULT`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NO_DEFAULT)


    * [`NetcdfReader`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader)


        * [`NetcdfReader.Error`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader.Error)


        * [`NetcdfReader.close()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader.close)


        * [`NetcdfReader.print_tree()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader.print_tree)


        * [`NetcdfReader.read_dimvalue()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_dimvalue)


        * [`NetcdfReader.read_keys()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_keys)


        * [`NetcdfReader.read_value()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_value)


        * [`NetcdfReader.read_variable()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_variable)


        * [`NetcdfReader.read_varnames()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_varnames)


        * [`NetcdfReader.walk_tree()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.NetcdfReader.walk_tree)


    * [`as_etsfreader()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.as_etsfreader)


    * [`as_ncreader()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.as_ncreader)


    * [`structure_from_ncdata()`](pymatgen.io.abinit.netcdf.md#pymatgen.io.abinit.netcdf.structure_from_ncdata)


* [pymatgen.io.abinit.pseudos module](pymatgen.io.abinit.pseudos.md)


    * [`AbinitHeader`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.AbinitHeader)


    * [`AbinitPseudo`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.AbinitPseudo)


        * [`AbinitPseudo.Z`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.AbinitPseudo.Z)


        * [`AbinitPseudo.Z_val`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.AbinitPseudo.Z_val)


        * [`AbinitPseudo.l_local`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.AbinitPseudo.l_local)


        * [`AbinitPseudo.l_max`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.AbinitPseudo.l_max)


        * [`AbinitPseudo.summary`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.AbinitPseudo.summary)


        * [`AbinitPseudo.supports_soc`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.AbinitPseudo.supports_soc)


    * [`Hint`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Hint)


        * [`Hint.as_dict()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Hint.as_dict)


        * [`Hint.from_dict()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Hint.from_dict)


    * [`NcAbinitHeader`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitHeader)


        * [`NcAbinitHeader.fhi_header()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.fhi_header)


        * [`NcAbinitHeader.gth_header()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.gth_header)


        * [`NcAbinitHeader.hgh_header()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.hgh_header)


        * [`NcAbinitHeader.oncvpsp_header()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.oncvpsp_header)


        * [`NcAbinitHeader.tm_header()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.tm_header)


    * [`NcAbinitPseudo`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo)


        * [`NcAbinitPseudo.Z`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.Z)


        * [`NcAbinitPseudo.Z_val`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.Z_val)


        * [`NcAbinitPseudo.l_local`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.l_local)


        * [`NcAbinitPseudo.l_max`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.l_max)


        * [`NcAbinitPseudo.nlcc_radius`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.nlcc_radius)


        * [`NcAbinitPseudo.summary`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.summary)


    * [`NcPseudo`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcPseudo)


        * [`NcPseudo.has_nlcc`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcPseudo.has_nlcc)


        * [`NcPseudo.nlcc_radius`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcPseudo.nlcc_radius)


        * [`NcPseudo.rcore`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.NcPseudo.rcore)


    * [`PawAbinitHeader`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawAbinitHeader)


        * [`PawAbinitHeader.paw_header()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawAbinitHeader.paw_header)


    * [`PawAbinitPseudo`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawAbinitPseudo)


        * [`PawAbinitPseudo.paw_radius`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawAbinitPseudo.paw_radius)


        * [`PawAbinitPseudo.supports_soc`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawAbinitPseudo.supports_soc)


    * [`PawPseudo`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawPseudo)


        * [`PawPseudo.paw_radius`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawPseudo.paw_radius)


        * [`PawPseudo.rcore`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawPseudo.rcore)


    * [`PawXmlSetup`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup)


        * [`PawXmlSetup.Z`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.Z)


        * [`PawXmlSetup.Z_val`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.Z_val)


        * [`PawXmlSetup.ae_core_density()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.ae_core_density)


        * [`PawXmlSetup.ae_partial_waves()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.ae_partial_waves)


        * [`PawXmlSetup.l_local`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.l_local)


        * [`PawXmlSetup.l_max`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.l_max)


        * [`PawXmlSetup.paw_radius`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.paw_radius)


        * [`PawXmlSetup.plot_densities()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.plot_densities)


        * [`PawXmlSetup.plot_projectors()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.plot_projectors)


        * [`PawXmlSetup.plot_waves()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.plot_waves)


        * [`PawXmlSetup.projector_functions()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.projector_functions)


        * [`PawXmlSetup.pseudo_core_density()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.pseudo_core_density)


        * [`PawXmlSetup.pseudo_partial_waves`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.pseudo_partial_waves)


        * [`PawXmlSetup.root()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.root)


        * [`PawXmlSetup.summary`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.summary)


        * [`PawXmlSetup.supports_soc`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.supports_soc)


        * [`PawXmlSetup.yield_figs()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PawXmlSetup.yield_figs)


    * [`Pseudo`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo)


        * [`Pseudo.Z`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.Z)


        * [`Pseudo.Z_val`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.Z_val)


        * [`Pseudo.as_dict()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.as_dict)


        * [`Pseudo.as_pseudo()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.as_pseudo)


        * [`Pseudo.as_tmpfile()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.as_tmpfile)


        * [`Pseudo.basename`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.basename)


        * [`Pseudo.compute_md5()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.compute_md5)


        * [`Pseudo.djrepo_path`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.djrepo_path)


        * [`Pseudo.element`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.element)


        * [`Pseudo.filepath`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.filepath)


        * [`Pseudo.from_dict()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.from_dict)


        * [`Pseudo.from_file()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.from_file)


        * [`Pseudo.has_dojo_report`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.has_dojo_report)


        * [`Pseudo.has_hints`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.has_hints)


        * [`Pseudo.hint_for_accuracy()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.hint_for_accuracy)


        * [`Pseudo.isnc`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.isnc)


        * [`Pseudo.ispaw`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.ispaw)


        * [`Pseudo.l_local`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.l_local)


        * [`Pseudo.l_max`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.l_max)


        * [`Pseudo.md5()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.md5)


        * [`Pseudo.open_pspsfile()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.open_pspsfile)


        * [`Pseudo.summary`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.summary)


        * [`Pseudo.supports_soc`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.supports_soc)


        * [`Pseudo.symbol`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.symbol)


        * [`Pseudo.to_str()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.to_str)


        * [`Pseudo.to_string()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.to_string)


        * [`Pseudo.type`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.Pseudo.type)


    * [`PseudoParseError`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoParseError)


    * [`PseudoParser`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoParser)


        * [`PseudoParser.Error`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoParser.Error)


        * [`PseudoParser.parse()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoParser.parse)


        * [`PseudoParser.read_ppdesc()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoParser.read_ppdesc)


        * [`PseudoParser.scan_directory()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoParser.scan_directory)


    * [`PseudoTable`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable)


        * [`PseudoTable.all_combinations_for_elements()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.all_combinations_for_elements)


        * [`PseudoTable.allnc`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.allnc)


        * [`PseudoTable.allpaw`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.allpaw)


        * [`PseudoTable.as_dict()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.as_dict)


        * [`PseudoTable.as_table()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.as_table)


        * [`PseudoTable.from_dict()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.from_dict)


        * [`PseudoTable.from_dir()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.from_dir)


        * [`PseudoTable.get_pseudos_for_structure()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.get_pseudos_for_structure)


        * [`PseudoTable.is_complete()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.is_complete)


        * [`PseudoTable.print_table()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.print_table)


        * [`PseudoTable.pseudo_with_symbol()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.pseudo_with_symbol)


        * [`PseudoTable.pseudos_with_symbols()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.pseudos_with_symbols)


        * [`PseudoTable.select()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.select)


        * [`PseudoTable.select_family()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.select_family)


        * [`PseudoTable.select_rows()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.select_rows)


        * [`PseudoTable.select_symbols()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.select_symbols)


        * [`PseudoTable.sort_by_z()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.sort_by_z)


        * [`PseudoTable.sorted()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.sorted)


        * [`PseudoTable.to_table()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.to_table)


        * [`PseudoTable.with_dojo_report()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.with_dojo_report)


        * [`PseudoTable.zlist`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.PseudoTable.zlist)


    * [`RadialFunction`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.RadialFunction)


    * [`l2str()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.l2str)


    * [`str2l()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.str2l)


    * [`straceback()`](pymatgen.io.abinit.pseudos.md#pymatgen.io.abinit.pseudos.straceback)


* [pymatgen.io.abinit.variable module](pymatgen.io.abinit.variable.md)


    * [`InputVariable`](pymatgen.io.abinit.variable.md#pymatgen.io.abinit.variable.InputVariable)


        * [`InputVariable.basename`](pymatgen.io.abinit.variable.md#pymatgen.io.abinit.variable.InputVariable.basename)


        * [`InputVariable.dataset`](pymatgen.io.abinit.variable.md#pymatgen.io.abinit.variable.InputVariable.dataset)


        * [`InputVariable.format_list()`](pymatgen.io.abinit.variable.md#pymatgen.io.abinit.variable.InputVariable.format_list)


        * [`InputVariable.format_list2d()`](pymatgen.io.abinit.variable.md#pymatgen.io.abinit.variable.InputVariable.format_list2d)


        * [`InputVariable.format_scalar()`](pymatgen.io.abinit.variable.md#pymatgen.io.abinit.variable.InputVariable.format_scalar)


        * [`InputVariable.get_value()`](pymatgen.io.abinit.variable.md#pymatgen.io.abinit.variable.InputVariable.get_value)


        * [`InputVariable.name`](pymatgen.io.abinit.variable.md#pymatgen.io.abinit.variable.InputVariable.name)


        * [`InputVariable.units`](pymatgen.io.abinit.variable.md#pymatgen.io.abinit.variable.InputVariable.units)