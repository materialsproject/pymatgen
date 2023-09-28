---
layout: default
title: pymatgen.io.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.io namespace

## Subpackages


* [pymatgen.io.abinit package](pymatgen.io.abinit.md)




    * [pymatgen.io.abinit.abiobjects module](pymatgen.io.abinit.md#module-pymatgen.io.abinit.abiobjects)


        * [`AbivarAble`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.AbivarAble)


            * [`AbivarAble._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.AbivarAble._abc_impl)


            * [`AbivarAble.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.AbivarAble.to_abivars)


        * [`Constraints`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Constraints)


            * [`Constraints._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Constraints._abc_impl)


            * [`Constraints.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Constraints.to_abivars)


        * [`Electrons`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Electrons)


            * [`Electrons._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Electrons._abc_impl)


            * [`Electrons.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Electrons.as_dict)


            * [`Electrons.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Electrons.from_dict)


            * [`Electrons.nspden`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Electrons.nspden)


            * [`Electrons.nspinor`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Electrons.nspinor)


            * [`Electrons.nsppol`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Electrons.nsppol)


            * [`Electrons.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Electrons.to_abivars)


        * [`ElectronsAlgorithm`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm)


            * [`ElectronsAlgorithm._DEFAULT`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm._DEFAULT)


            * [`ElectronsAlgorithm._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm._abc_impl)


            * [`ElectronsAlgorithm.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm.as_dict)


            * [`ElectronsAlgorithm.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm.from_dict)


            * [`ElectronsAlgorithm.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ElectronsAlgorithm.to_abivars)


        * [`ExcHamiltonian`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian)


            * [`ExcHamiltonian._ALGO2VAR`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian._ALGO2VAR)


            * [`ExcHamiltonian._COULOMB_MODES`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian._COULOMB_MODES)


            * [`ExcHamiltonian._EXC_TYPES`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian._EXC_TYPES)


            * [`ExcHamiltonian._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian._abc_impl)


            * [`ExcHamiltonian.inclvkb`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.inclvkb)


            * [`ExcHamiltonian.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.to_abivars)


            * [`ExcHamiltonian.use_cg`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.use_cg)


            * [`ExcHamiltonian.use_direct_diago`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.use_direct_diago)


            * [`ExcHamiltonian.use_haydock`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ExcHamiltonian.use_haydock)


        * [`HilbertTransform`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.HilbertTransform)


            * [`HilbertTransform._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.HilbertTransform._abc_impl)


            * [`HilbertTransform.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.HilbertTransform.to_abivars)


        * [`KSampling`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling)


            * [`KSampling._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling._abc_impl)


            * [`KSampling._path()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling._path)


            * [`KSampling.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.as_dict)


            * [`KSampling.automatic_density()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.automatic_density)


            * [`KSampling.explicit_path()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.explicit_path)


            * [`KSampling.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.from_dict)


            * [`KSampling.gamma_centered()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.gamma_centered)


            * [`KSampling.gamma_only()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.gamma_only)


            * [`KSampling.is_homogeneous`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.is_homogeneous)


            * [`KSampling.monkhorst()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.monkhorst)


            * [`KSampling.monkhorst_automatic()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.monkhorst_automatic)


            * [`KSampling.path_from_structure()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.path_from_structure)


            * [`KSampling.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSampling.to_abivars)


        * [`KSamplingModes`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSamplingModes)


            * [`KSamplingModes.automatic`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSamplingModes.automatic)


            * [`KSamplingModes.monkhorst`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSamplingModes.monkhorst)


            * [`KSamplingModes.path`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.KSamplingModes.path)


        * [`ModelDielectricFunction`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ModelDielectricFunction)


            * [`ModelDielectricFunction._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ModelDielectricFunction._abc_impl)


            * [`ModelDielectricFunction.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.ModelDielectricFunction.to_abivars)


        * [`PPModel`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModel)


            * [`PPModel._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModel._abc_impl)


            * [`PPModel.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModel.as_dict)


            * [`PPModel.as_ppmodel()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModel.as_ppmodel)


            * [`PPModel.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModel.from_dict)


            * [`PPModel.get_noppmodel()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModel.get_noppmodel)


            * [`PPModel.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModel.to_abivars)


        * [`PPModelModes`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModelModes)


            * [`PPModelModes.farid`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModelModes.farid)


            * [`PPModelModes.godby`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModelModes.godby)


            * [`PPModelModes.hybersten`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModelModes.hybersten)


            * [`PPModelModes.linden`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModelModes.linden)


            * [`PPModelModes.noppmodel`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.PPModelModes.noppmodel)


        * [`RelaxationMethod`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod)


            * [`RelaxationMethod.IONMOV_DEFAULT`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.IONMOV_DEFAULT)


            * [`RelaxationMethod.OPTCELL_DEFAULT`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.OPTCELL_DEFAULT)


            * [`RelaxationMethod._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod._abc_impl)


            * [`RelaxationMethod._default_vars`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod._default_vars)


            * [`RelaxationMethod.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.as_dict)


            * [`RelaxationMethod.atoms_and_cell()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.atoms_and_cell)


            * [`RelaxationMethod.atoms_only()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.atoms_only)


            * [`RelaxationMethod.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.from_dict)


            * [`RelaxationMethod.move_atoms`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.move_atoms)


            * [`RelaxationMethod.move_cell`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.move_cell)


            * [`RelaxationMethod.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.RelaxationMethod.to_abivars)


        * [`Screening`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Screening)


            * [`Screening._SC_MODES`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Screening._SC_MODES)


            * [`Screening._WTYPES`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Screening._WTYPES)


            * [`Screening._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Screening._abc_impl)


            * [`Screening.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Screening.to_abivars)


            * [`Screening.use_hilbert`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Screening.use_hilbert)


        * [`SelfEnergy`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SelfEnergy)


            * [`SelfEnergy._SC_MODES`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SelfEnergy._SC_MODES)


            * [`SelfEnergy._SIGMA_TYPES`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SelfEnergy._SIGMA_TYPES)


            * [`SelfEnergy._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SelfEnergy._abc_impl)


            * [`SelfEnergy.gwcalctyp`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SelfEnergy.gwcalctyp)


            * [`SelfEnergy.symsigma`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SelfEnergy.symsigma)


            * [`SelfEnergy.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SelfEnergy.to_abivars)


            * [`SelfEnergy.use_ppmodel`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SelfEnergy.use_ppmodel)


        * [`Smearing`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Smearing)


            * [`Smearing._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Smearing._abc_impl)


            * [`Smearing._mode2occopt`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Smearing._mode2occopt)


            * [`Smearing.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Smearing.as_dict)


            * [`Smearing.as_smearing()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Smearing.as_smearing)


            * [`Smearing.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Smearing.from_dict)


            * [`Smearing.mode`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Smearing.mode)


            * [`Smearing.nosmearing()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Smearing.nosmearing)


            * [`Smearing.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.Smearing.to_abivars)


        * [`SpinMode`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SpinMode)


            * [`SpinMode._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SpinMode._abc_impl)


            * [`SpinMode.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SpinMode.as_dict)


            * [`SpinMode.as_spinmode()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SpinMode.as_spinmode)


            * [`SpinMode.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SpinMode.from_dict)


            * [`SpinMode.to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.SpinMode.to_abivars)


        * [`contract()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.contract)


        * [`lattice_from_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.lattice_from_abivars)


        * [`species_by_znucl()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.species_by_znucl)


        * [`structure_from_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.structure_from_abivars)


        * [`structure_to_abivars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abiobjects.structure_to_abivars)


    * [pymatgen.io.abinit.abitimer module](pymatgen.io.abinit.md#module-pymatgen.io.abinit.abitimer)


        * [`AbinitTimer`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer)


            * [`AbinitTimer._reduce_sections()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer._reduce_sections)


            * [`AbinitTimer.cpuwall_histogram()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.cpuwall_histogram)


            * [`AbinitTimer.get_dataframe()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.get_dataframe)


            * [`AbinitTimer.get_section()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.get_section)


            * [`AbinitTimer.get_values()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.get_values)


            * [`AbinitTimer.names_and_values()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.names_and_values)


            * [`AbinitTimer.ncpus`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.ncpus)


            * [`AbinitTimer.order_sections()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.order_sections)


            * [`AbinitTimer.pie()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.pie)


            * [`AbinitTimer.scatter_hist()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.scatter_hist)


            * [`AbinitTimer.sum_sections()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.sum_sections)


            * [`AbinitTimer.to_csv()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.to_csv)


            * [`AbinitTimer.to_table()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.to_table)


            * [`AbinitTimer.totable()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimer.totable)


        * [`AbinitTimerParseError`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParseError)


        * [`AbinitTimerParser`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser)


            * [`AbinitTimerParser.BEGIN_TAG`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.BEGIN_TAG)


            * [`AbinitTimerParser.END_TAG`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.END_TAG)


            * [`AbinitTimerParser.Error`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.Error)


            * [`AbinitTimerParser._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser._abc_impl)


            * [`AbinitTimerParser._read()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser._read)


            * [`AbinitTimerParser.filenames`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.filenames)


            * [`AbinitTimerParser.get_sections()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.get_sections)


            * [`AbinitTimerParser.parse()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.parse)


            * [`AbinitTimerParser.pefficiency()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.pefficiency)


            * [`AbinitTimerParser.plot_all()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.plot_all)


            * [`AbinitTimerParser.plot_efficiency()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.plot_efficiency)


            * [`AbinitTimerParser.plot_pie()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.plot_pie)


            * [`AbinitTimerParser.plot_stacked_hist()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.plot_stacked_hist)


            * [`AbinitTimerParser.section_names()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.section_names)


            * [`AbinitTimerParser.summarize()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.summarize)


            * [`AbinitTimerParser.timers()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.timers)


            * [`AbinitTimerParser.walk()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerParser.walk)


        * [`AbinitTimerSection`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerSection)


            * [`AbinitTimerSection.FIELDS`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.FIELDS)


            * [`AbinitTimerSection.NUMERIC_FIELDS`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.NUMERIC_FIELDS)


            * [`AbinitTimerSection.STR_FIELDS`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.STR_FIELDS)


            * [`AbinitTimerSection.fake()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.fake)


            * [`AbinitTimerSection.to_csvline()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.to_csvline)


            * [`AbinitTimerSection.to_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.to_dict)


            * [`AbinitTimerSection.to_tuple()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.AbinitTimerSection.to_tuple)


        * [`ParallelEfficiency`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.ParallelEfficiency)


            * [`ParallelEfficiency._order_by_peff()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.ParallelEfficiency._order_by_peff)


            * [`ParallelEfficiency.bad_sections()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.ParallelEfficiency.bad_sections)


            * [`ParallelEfficiency.good_sections()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.ParallelEfficiency.good_sections)


            * [`ParallelEfficiency.totable()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.ParallelEfficiency.totable)


        * [`alternate()`](pymatgen.io.abinit.md#pymatgen.io.abinit.abitimer.alternate)


    * [pymatgen.io.abinit.inputs module](pymatgen.io.abinit.md#module-pymatgen.io.abinit.inputs)


        * [`AbstractInput`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput)


            * [`AbstractInput._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput._abc_impl)


            * [`AbstractInput._check_varname()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput._check_varname)


            * [`AbstractInput.deepcopy()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput.deepcopy)


            * [`AbstractInput.pop_vars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput.pop_vars)


            * [`AbstractInput.remove_vars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput.remove_vars)


            * [`AbstractInput.set_vars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput.set_vars)


            * [`AbstractInput.set_vars_ifnotin()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput.set_vars_ifnotin)


            * [`AbstractInput.to_str()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput.to_str)


            * [`AbstractInput.vars`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput.vars)


            * [`AbstractInput.write()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.AbstractInput.write)


        * [`BasicAbinitInput`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput)


            * [`BasicAbinitInput.Error`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.Error)


            * [`BasicAbinitInput._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput._abc_impl)


            * [`BasicAbinitInput._check_varname()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput._check_varname)


            * [`BasicAbinitInput.add_abiobjects()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.add_abiobjects)


            * [`BasicAbinitInput.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.as_dict)


            * [`BasicAbinitInput.comment`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.comment)


            * [`BasicAbinitInput.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.from_dict)


            * [`BasicAbinitInput.isnc`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.isnc)


            * [`BasicAbinitInput.ispaw`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.ispaw)


            * [`BasicAbinitInput.new_with_vars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.new_with_vars)


            * [`BasicAbinitInput.pop_irdvars()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.pop_irdvars)


            * [`BasicAbinitInput.pop_tolerances()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.pop_tolerances)


            * [`BasicAbinitInput.pseudos`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.pseudos)


            * [`BasicAbinitInput.set_comment()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_comment)


            * [`BasicAbinitInput.set_gamma_sampling()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_gamma_sampling)


            * [`BasicAbinitInput.set_kmesh()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_kmesh)


            * [`BasicAbinitInput.set_kpath()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_kpath)


            * [`BasicAbinitInput.set_spin_mode()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_spin_mode)


            * [`BasicAbinitInput.set_structure()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.set_structure)


            * [`BasicAbinitInput.structure`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.structure)


            * [`BasicAbinitInput.to_str()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.to_str)


            * [`BasicAbinitInput.to_string()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.to_string)


            * [`BasicAbinitInput.vars`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInput.vars)


        * [`BasicAbinitInputError`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicAbinitInputError)


        * [`BasicMultiDataset`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset)


            * [`BasicMultiDataset.Error`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.Error)


            * [`BasicMultiDataset.addnew_from()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.addnew_from)


            * [`BasicMultiDataset.append()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.append)


            * [`BasicMultiDataset.deepcopy()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.deepcopy)


            * [`BasicMultiDataset.extend()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.extend)


            * [`BasicMultiDataset.from_inputs()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.from_inputs)


            * [`BasicMultiDataset.has_same_structures`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.has_same_structures)


            * [`BasicMultiDataset.isnc`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.isnc)


            * [`BasicMultiDataset.ispaw`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.ispaw)


            * [`BasicMultiDataset.ndtset`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.ndtset)


            * [`BasicMultiDataset.pseudos`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.pseudos)


            * [`BasicMultiDataset.replicate_input()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.replicate_input)


            * [`BasicMultiDataset.split_datasets()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.split_datasets)


            * [`BasicMultiDataset.to_str()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.to_str)


            * [`BasicMultiDataset.to_string()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.to_string)


            * [`BasicMultiDataset.write()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.BasicMultiDataset.write)


        * [`ShiftMode`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.ShiftMode)


            * [`ShiftMode.GammaCentered`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.ShiftMode.GammaCentered)


            * [`ShiftMode.MonkhorstPack`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.ShiftMode.MonkhorstPack)


            * [`ShiftMode.OneSymmetric`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.ShiftMode.OneSymmetric)


            * [`ShiftMode.Symmetric`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.ShiftMode.Symmetric)


            * [`ShiftMode.from_object()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.ShiftMode.from_object)


        * [`_find_ecut_pawecutdg()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs._find_ecut_pawecutdg)


        * [`_find_scf_nband()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs._find_scf_nband)


        * [`_get_shifts()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs._get_shifts)


        * [`_stopping_criterion()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs._stopping_criterion)


        * [`as_structure()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.as_structure)


        * [`calc_shiftk()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.calc_shiftk)


        * [`ebands_input()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.ebands_input)


        * [`gs_input()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.gs_input)


        * [`ion_ioncell_relax_input()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.ion_ioncell_relax_input)


        * [`num_valence_electrons()`](pymatgen.io.abinit.md#pymatgen.io.abinit.inputs.num_valence_electrons)


    * [pymatgen.io.abinit.netcdf module](pymatgen.io.abinit.md#module-pymatgen.io.abinit.netcdf)


        * [`AbinitHeader`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.AbinitHeader)


            * [`AbinitHeader.to_str()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.AbinitHeader.to_str)


            * [`AbinitHeader.to_string()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.AbinitHeader.to_string)


        * [`ETSF_Reader`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.ETSF_Reader)


            * [`ETSF_Reader.chemical_symbols()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.ETSF_Reader.chemical_symbols)


            * [`ETSF_Reader.read_abinit_hdr()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.ETSF_Reader.read_abinit_hdr)


            * [`ETSF_Reader.read_abinit_xcfunc()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.ETSF_Reader.read_abinit_xcfunc)


            * [`ETSF_Reader.read_structure()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.ETSF_Reader.read_structure)


            * [`ETSF_Reader.typeidx_from_symbol()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.ETSF_Reader.typeidx_from_symbol)


        * [`NO_DEFAULT`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NO_DEFAULT)


        * [`NetcdfReader`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader)


            * [`NetcdfReader.Error`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader.Error)


            * [`NetcdfReader._read_dimensions()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader._read_dimensions)


            * [`NetcdfReader._read_variables()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader._read_variables)


            * [`NetcdfReader.close()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader.close)


            * [`NetcdfReader.print_tree()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader.print_tree)


            * [`NetcdfReader.read_dimvalue()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_dimvalue)


            * [`NetcdfReader.read_keys()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_keys)


            * [`NetcdfReader.read_value()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_value)


            * [`NetcdfReader.read_variable()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_variable)


            * [`NetcdfReader.read_varnames()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader.read_varnames)


            * [`NetcdfReader.walk_tree()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReader.walk_tree)


        * [`NetcdfReaderError`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.NetcdfReaderError)


        * [`_H`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf._H)


            * [`_H.doc`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf._H.doc)


            * [`_H.etsf_name`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf._H.etsf_name)


            * [`_H.name`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf._H.name)


        * [`_asreader()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf._asreader)


        * [`as_etsfreader()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.as_etsfreader)


        * [`as_ncreader()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.as_ncreader)


        * [`structure_from_ncdata()`](pymatgen.io.abinit.md#pymatgen.io.abinit.netcdf.structure_from_ncdata)


    * [pymatgen.io.abinit.pseudos module](pymatgen.io.abinit.md#module-pymatgen.io.abinit.pseudos)


        * [`AbinitHeader`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.AbinitHeader)


        * [`AbinitPseudo`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.AbinitPseudo)


            * [`AbinitPseudo.Z`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.AbinitPseudo.Z)


            * [`AbinitPseudo.Z_val`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.AbinitPseudo.Z_val)


            * [`AbinitPseudo._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.AbinitPseudo._abc_impl)


            * [`AbinitPseudo.l_local`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.AbinitPseudo.l_local)


            * [`AbinitPseudo.l_max`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.AbinitPseudo.l_max)


            * [`AbinitPseudo.summary`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.AbinitPseudo.summary)


            * [`AbinitPseudo.supports_soc`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.AbinitPseudo.supports_soc)


        * [`Hint`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Hint)


            * [`Hint.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Hint.as_dict)


            * [`Hint.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Hint.from_dict)


        * [`NcAbinitHeader`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitHeader)


            * [`NcAbinitHeader._VARS`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitHeader._VARS)


            * [`NcAbinitHeader.fhi_header()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.fhi_header)


            * [`NcAbinitHeader.gth_header()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.gth_header)


            * [`NcAbinitHeader.hgh_header()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.hgh_header)


            * [`NcAbinitHeader.oncvpsp_header()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.oncvpsp_header)


            * [`NcAbinitHeader.tm_header()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitHeader.tm_header)


        * [`NcAbinitPseudo`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo)


            * [`NcAbinitPseudo.Z`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.Z)


            * [`NcAbinitPseudo.Z_val`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.Z_val)


            * [`NcAbinitPseudo._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo._abc_impl)


            * [`NcAbinitPseudo.l_local`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.l_local)


            * [`NcAbinitPseudo.l_max`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.l_max)


            * [`NcAbinitPseudo.nlcc_radius`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.nlcc_radius)


            * [`NcAbinitPseudo.summary`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcAbinitPseudo.summary)


        * [`NcPseudo`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcPseudo)


            * [`NcPseudo._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcPseudo._abc_impl)


            * [`NcPseudo.has_nlcc`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcPseudo.has_nlcc)


            * [`NcPseudo.nlcc_radius`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcPseudo.nlcc_radius)


            * [`NcPseudo.rcore`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.NcPseudo.rcore)


        * [`PawAbinitHeader`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawAbinitHeader)


            * [`PawAbinitHeader._VARS`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawAbinitHeader._VARS)


            * [`PawAbinitHeader.paw_header()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawAbinitHeader.paw_header)


        * [`PawAbinitPseudo`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawAbinitPseudo)


            * [`PawAbinitPseudo._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawAbinitPseudo._abc_impl)


            * [`PawAbinitPseudo.paw_radius`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawAbinitPseudo.paw_radius)


            * [`PawAbinitPseudo.supports_soc`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawAbinitPseudo.supports_soc)


        * [`PawPseudo`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawPseudo)


            * [`PawPseudo._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawPseudo._abc_impl)


            * [`PawPseudo.paw_radius`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawPseudo.paw_radius)


            * [`PawPseudo.rcore`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawPseudo.rcore)


        * [`PawXmlSetup`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup)


            * [`PawXmlSetup.Z`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.Z)


            * [`PawXmlSetup.Z_val`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.Z_val)


            * [`PawXmlSetup._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup._abc_impl)


            * [`PawXmlSetup._eval_grid()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup._eval_grid)


            * [`PawXmlSetup._parse_all_radfuncs()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup._parse_all_radfuncs)


            * [`PawXmlSetup._parse_radfunc()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup._parse_radfunc)


            * [`PawXmlSetup.ae_core_density()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.ae_core_density)


            * [`PawXmlSetup.ae_partial_waves()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.ae_partial_waves)


            * [`PawXmlSetup.l_local`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.l_local)


            * [`PawXmlSetup.l_max`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.l_max)


            * [`PawXmlSetup.paw_radius`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.paw_radius)


            * [`PawXmlSetup.plot_densities()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.plot_densities)


            * [`PawXmlSetup.plot_projectors()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.plot_projectors)


            * [`PawXmlSetup.plot_waves()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.plot_waves)


            * [`PawXmlSetup.projector_functions()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.projector_functions)


            * [`PawXmlSetup.pseudo_core_density()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.pseudo_core_density)


            * [`PawXmlSetup.pseudo_partial_waves`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.pseudo_partial_waves)


            * [`PawXmlSetup.root()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.root)


            * [`PawXmlSetup.summary`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.summary)


            * [`PawXmlSetup.supports_soc`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.supports_soc)


            * [`PawXmlSetup.yield_figs()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PawXmlSetup.yield_figs)


        * [`Pseudo`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo)


            * [`Pseudo.Z`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.Z)


            * [`Pseudo.Z_val`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.Z_val)


            * [`Pseudo._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo._abc_impl)


            * [`Pseudo.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.as_dict)


            * [`Pseudo.as_pseudo()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.as_pseudo)


            * [`Pseudo.as_tmpfile()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.as_tmpfile)


            * [`Pseudo.basename`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.basename)


            * [`Pseudo.compute_md5()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.compute_md5)


            * [`Pseudo.djrepo_path`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.djrepo_path)


            * [`Pseudo.element`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.element)


            * [`Pseudo.filepath`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.filepath)


            * [`Pseudo.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.from_dict)


            * [`Pseudo.from_file()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.from_file)


            * [`Pseudo.has_dojo_report`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.has_dojo_report)


            * [`Pseudo.has_hints`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.has_hints)


            * [`Pseudo.hint_for_accuracy()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.hint_for_accuracy)


            * [`Pseudo.isnc`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.isnc)


            * [`Pseudo.ispaw`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.ispaw)


            * [`Pseudo.l_local`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.l_local)


            * [`Pseudo.l_max`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.l_max)


            * [`Pseudo.md5()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.md5)


            * [`Pseudo.open_pspsfile()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.open_pspsfile)


            * [`Pseudo.summary`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.summary)


            * [`Pseudo.supports_soc`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.supports_soc)


            * [`Pseudo.symbol`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.symbol)


            * [`Pseudo.to_str()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.to_str)


            * [`Pseudo.to_string()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.to_string)


            * [`Pseudo.type`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.Pseudo.type)


        * [`PseudoParseError`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoParseError)


        * [`PseudoParser`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoParser)


            * [`PseudoParser.Error`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoParser.Error)


            * [`PseudoParser._PSPCODES`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoParser._PSPCODES)


            * [`PseudoParser.parse()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoParser.parse)


            * [`PseudoParser.read_ppdesc()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoParser.read_ppdesc)


            * [`PseudoParser.scan_directory()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoParser.scan_directory)


        * [`PseudoTable`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable)


            * [`PseudoTable._abc_impl`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable._abc_impl)


            * [`PseudoTable.all_combinations_for_elements()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.all_combinations_for_elements)


            * [`PseudoTable.allnc`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.allnc)


            * [`PseudoTable.allpaw`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.allpaw)


            * [`PseudoTable.as_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.as_dict)


            * [`PseudoTable.as_table()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.as_table)


            * [`PseudoTable.from_dict()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.from_dict)


            * [`PseudoTable.from_dir()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.from_dir)


            * [`PseudoTable.get_pseudos_for_structure()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.get_pseudos_for_structure)


            * [`PseudoTable.is_complete()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.is_complete)


            * [`PseudoTable.print_table()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.print_table)


            * [`PseudoTable.pseudo_with_symbol()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.pseudo_with_symbol)


            * [`PseudoTable.pseudos_with_symbols()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.pseudos_with_symbols)


            * [`PseudoTable.select()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.select)


            * [`PseudoTable.select_family()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.select_family)


            * [`PseudoTable.select_rows()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.select_rows)


            * [`PseudoTable.select_symbols()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.select_symbols)


            * [`PseudoTable.sort_by_z()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.sort_by_z)


            * [`PseudoTable.sorted()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.sorted)


            * [`PseudoTable.to_table()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.to_table)


            * [`PseudoTable.with_dojo_report()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.with_dojo_report)


            * [`PseudoTable.zlist`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.PseudoTable.zlist)


        * [`RadialFunction`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.RadialFunction)


        * [`_dict_from_lines()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos._dict_from_lines)


        * [`_int_from_str()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos._int_from_str)


        * [`_read_nlines()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos._read_nlines)


        * [`l2str()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.l2str)


        * [`str2l()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.str2l)


        * [`straceback()`](pymatgen.io.abinit.md#pymatgen.io.abinit.pseudos.straceback)


    * [pymatgen.io.abinit.variable module](pymatgen.io.abinit.md#module-pymatgen.io.abinit.variable)


        * [`InputVariable`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.InputVariable)


            * [`InputVariable.basename`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.InputVariable.basename)


            * [`InputVariable.dataset`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.InputVariable.dataset)


            * [`InputVariable.format_list()`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.InputVariable.format_list)


            * [`InputVariable.format_list2d()`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.InputVariable.format_list2d)


            * [`InputVariable.format_scalar()`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.InputVariable.format_scalar)


            * [`InputVariable.get_value()`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.InputVariable.get_value)


            * [`InputVariable.name`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.InputVariable.name)


            * [`InputVariable.units`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.InputVariable.units)


        * [`flatten()`](pymatgen.io.abinit.md#pymatgen.io.abinit.variable.flatten)


* [pymatgen.io.cp2k package](pymatgen.io.cp2k.md)




    * [pymatgen.io.cp2k.inputs module](pymatgen.io.cp2k.md#module-pymatgen.io.cp2k.inputs)


        * [`AtomicMetadata`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata)


            * [`AtomicMetadata.info`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.info)


            * [`AtomicMetadata.element`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.element)


            * [`AtomicMetadata.potential`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.potential)


            * [`AtomicMetadata.name`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.name)


            * [`AtomicMetadata.alias_names`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.alias_names)


            * [`AtomicMetadata.filename`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.filename)


            * [`AtomicMetadata.version`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.version)


            * [`AtomicMetadata.alias_names`](pymatgen.io.cp2k.md#id0)


            * [`AtomicMetadata.element`](pymatgen.io.cp2k.md#id1)


            * [`AtomicMetadata.filename`](pymatgen.io.cp2k.md#id2)


            * [`AtomicMetadata.get_hash()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.get_hash)


            * [`AtomicMetadata.get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.get_str)


            * [`AtomicMetadata.get_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.get_string)


            * [`AtomicMetadata.info`](pymatgen.io.cp2k.md#id3)


            * [`AtomicMetadata.name`](pymatgen.io.cp2k.md#id4)


            * [`AtomicMetadata.potential`](pymatgen.io.cp2k.md#id5)


            * [`AtomicMetadata.softmatch()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.AtomicMetadata.softmatch)


            * [`AtomicMetadata.version`](pymatgen.io.cp2k.md#id6)


        * [`Band_Structure`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Band_Structure)


            * [`Band_Structure.from_kpoints()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Band_Structure.from_kpoints)


        * [`BasisFile`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisFile)


            * [`BasisFile.from_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisFile.from_str)


        * [`BasisInfo`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo)


            * [`BasisInfo.electrons`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.electrons)


            * [`BasisInfo.core`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.core)


            * [`BasisInfo.valence`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.valence)


            * [`BasisInfo.polarization`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.polarization)


            * [`BasisInfo.diffuse`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.diffuse)


            * [`BasisInfo.cc`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.cc)


            * [`BasisInfo.pc`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.pc)


            * [`BasisInfo.sr`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.sr)


            * [`BasisInfo.molopt`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.molopt)


            * [`BasisInfo.admm`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.admm)


            * [`BasisInfo.lri`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.lri)


            * [`BasisInfo.contracted`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.contracted)


            * [`BasisInfo.xc`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.xc)


            * [`BasisInfo.admm`](pymatgen.io.cp2k.md#id7)


            * [`BasisInfo.cc`](pymatgen.io.cp2k.md#id8)


            * [`BasisInfo.contracted`](pymatgen.io.cp2k.md#id9)


            * [`BasisInfo.core`](pymatgen.io.cp2k.md#id10)


            * [`BasisInfo.diffuse`](pymatgen.io.cp2k.md#id11)


            * [`BasisInfo.electrons`](pymatgen.io.cp2k.md#id12)


            * [`BasisInfo.from_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.from_str)


            * [`BasisInfo.from_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.from_string)


            * [`BasisInfo.lri`](pymatgen.io.cp2k.md#id13)


            * [`BasisInfo.molopt`](pymatgen.io.cp2k.md#id14)


            * [`BasisInfo.pc`](pymatgen.io.cp2k.md#id15)


            * [`BasisInfo.polarization`](pymatgen.io.cp2k.md#id16)


            * [`BasisInfo.softmatch()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BasisInfo.softmatch)


            * [`BasisInfo.sr`](pymatgen.io.cp2k.md#id17)


            * [`BasisInfo.valence`](pymatgen.io.cp2k.md#id18)


            * [`BasisInfo.xc`](pymatgen.io.cp2k.md#id19)


        * [`BrokenSymmetry`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BrokenSymmetry)


            * [`BrokenSymmetry.from_el()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.BrokenSymmetry.from_el)


        * [`Cell`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cell)


        * [`Coord`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Coord)


        * [`Cp2kInput`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cp2kInput)


            * [`Cp2kInput._from_dict()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cp2kInput._from_dict)


            * [`Cp2kInput._from_lines()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cp2kInput._from_lines)


            * [`Cp2kInput.from_file()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cp2kInput.from_file)


            * [`Cp2kInput.from_lines()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cp2kInput.from_lines)


            * [`Cp2kInput.from_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cp2kInput.from_str)


            * [`Cp2kInput.from_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cp2kInput.from_string)


            * [`Cp2kInput.get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cp2kInput.get_str)


            * [`Cp2kInput.write_file()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Cp2kInput.write_file)


        * [`DOS`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DOS)


        * [`DataFile`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DataFile)


            * [`DataFile.from_file()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DataFile.from_file)


            * [`DataFile.from_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DataFile.from_str)


            * [`DataFile.from_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DataFile.from_string)


            * [`DataFile.get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DataFile.get_str)


            * [`DataFile.get_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DataFile.get_string)


            * [`DataFile.objects`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DataFile.objects)


            * [`DataFile.write_file()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DataFile.write_file)


        * [`Davidson`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Davidson)


        * [`Dft`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Dft)


        * [`DftPlusU`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.DftPlusU)


        * [`Diagonalization`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Diagonalization)


        * [`E_Density_Cube`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.E_Density_Cube)


        * [`ForceEval`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.ForceEval)


        * [`GaussianTypeOrbitalBasisSet`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet)


            * [`GaussianTypeOrbitalBasisSet.info`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.info)


            * [`GaussianTypeOrbitalBasisSet.nset`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.nset)


            * [`GaussianTypeOrbitalBasisSet.n`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.n)


            * [`GaussianTypeOrbitalBasisSet.lmax`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.lmax)


            * [`GaussianTypeOrbitalBasisSet.lmin`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.lmin)


            * [`GaussianTypeOrbitalBasisSet.nshell`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.nshell)


            * [`GaussianTypeOrbitalBasisSet.exponents`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.exponents)


            * [`GaussianTypeOrbitalBasisSet.coefficients`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.coefficients)


            * [`GaussianTypeOrbitalBasisSet.coefficients`](pymatgen.io.cp2k.md#id20)


            * [`GaussianTypeOrbitalBasisSet.exponents`](pymatgen.io.cp2k.md#id21)


            * [`GaussianTypeOrbitalBasisSet.from_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.from_str)


            * [`GaussianTypeOrbitalBasisSet.from_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.from_string)


            * [`GaussianTypeOrbitalBasisSet.get_keyword()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.get_keyword)


            * [`GaussianTypeOrbitalBasisSet.get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.get_str)


            * [`GaussianTypeOrbitalBasisSet.get_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.get_string)


            * [`GaussianTypeOrbitalBasisSet.info`](pymatgen.io.cp2k.md#id22)


            * [`GaussianTypeOrbitalBasisSet.lmax`](pymatgen.io.cp2k.md#id23)


            * [`GaussianTypeOrbitalBasisSet.lmin`](pymatgen.io.cp2k.md#id24)


            * [`GaussianTypeOrbitalBasisSet.n`](pymatgen.io.cp2k.md#id25)


            * [`GaussianTypeOrbitalBasisSet.nexp`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.nexp)


            * [`GaussianTypeOrbitalBasisSet.nset`](pymatgen.io.cp2k.md#id26)


            * [`GaussianTypeOrbitalBasisSet.nshell`](pymatgen.io.cp2k.md#id27)


        * [`Global`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Global)


        * [`GthPotential`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential)


            * [`GthPotential.info`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.info)


            * [`GthPotential.n_elecs`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.n_elecs)


            * [`GthPotential.r_loc`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.r_loc)


            * [`GthPotential.nexp_ppl`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.nexp_ppl)


            * [`GthPotential.c_exp_ppl`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.c_exp_ppl)


            * [`GthPotential.radii`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.radii)


            * [`GthPotential.nprj`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.nprj)


            * [`GthPotential.nprj_ppnl`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.nprj_ppnl)


            * [`GthPotential.hprj_ppnl`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.hprj_ppnl)


            * [`GthPotential.c_exp_ppl`](pymatgen.io.cp2k.md#id28)


            * [`GthPotential.from_section()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.from_section)


            * [`GthPotential.from_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.from_str)


            * [`GthPotential.from_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.from_string)


            * [`GthPotential.get_keyword()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.get_keyword)


            * [`GthPotential.get_section()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.get_section)


            * [`GthPotential.get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.get_str)


            * [`GthPotential.get_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.GthPotential.get_string)


            * [`GthPotential.hprj_ppnl`](pymatgen.io.cp2k.md#id29)


            * [`GthPotential.n_elecs`](pymatgen.io.cp2k.md#id30)


            * [`GthPotential.nexp_ppl`](pymatgen.io.cp2k.md#id31)


            * [`GthPotential.nprj`](pymatgen.io.cp2k.md#id32)


            * [`GthPotential.nprj_ppnl`](pymatgen.io.cp2k.md#id33)


            * [`GthPotential.r_loc`](pymatgen.io.cp2k.md#id34)


            * [`GthPotential.radii`](pymatgen.io.cp2k.md#id35)


        * [`Keyword`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Keyword)


            * [`Keyword.as_dict()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Keyword.as_dict)


            * [`Keyword.from_dict()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Keyword.from_dict)


            * [`Keyword.from_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Keyword.from_str)


            * [`Keyword.from_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Keyword.from_string)


            * [`Keyword.get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Keyword.get_str)


            * [`Keyword.get_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Keyword.get_string)


            * [`Keyword.verbosity()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Keyword.verbosity)


        * [`KeywordList`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.KeywordList)


            * [`KeywordList.append()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.KeywordList.append)


            * [`KeywordList.extend()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.KeywordList.extend)


            * [`KeywordList.get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.KeywordList.get_str)


            * [`KeywordList.get_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.KeywordList.get_string)


            * [`KeywordList.verbosity()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.KeywordList.verbosity)


        * [`Kind`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kind)


        * [`Kpoint_Set`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kpoint_Set)


        * [`Kpoints`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kpoints)


            * [`Kpoints.from_kpoints()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kpoints.from_kpoints)


        * [`LDOS`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.LDOS)


        * [`MO_Cubes`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.MO_Cubes)


        * [`Mgrid`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Mgrid)


        * [`OrbitalTransformation`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.OrbitalTransformation)


        * [`PBE`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PBE)


        * [`PDOS`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PDOS)


        * [`PotentialFile`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialFile)


            * [`PotentialFile.from_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialFile.from_str)


        * [`PotentialInfo`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialInfo)


            * [`PotentialInfo.electrons`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialInfo.electrons)


            * [`PotentialInfo.potential_type`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialInfo.potential_type)


            * [`PotentialInfo.nlcc`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialInfo.nlcc)


            * [`PotentialInfo.xc`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialInfo.xc)


            * [`PotentialInfo.electrons`](pymatgen.io.cp2k.md#id36)


            * [`PotentialInfo.from_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialInfo.from_str)


            * [`PotentialInfo.from_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialInfo.from_string)


            * [`PotentialInfo.nlcc`](pymatgen.io.cp2k.md#id37)


            * [`PotentialInfo.potential_type`](pymatgen.io.cp2k.md#id38)


            * [`PotentialInfo.softmatch()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.PotentialInfo.softmatch)


            * [`PotentialInfo.xc`](pymatgen.io.cp2k.md#id39)


        * [`QS`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.QS)


        * [`Scf`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Scf)


        * [`Section`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section)


            * [`Section._get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section._get_str)


            * [`Section._update()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section._update)


            * [`Section.add()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.add)


            * [`Section.by_path()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.by_path)


            * [`Section.check()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.check)


            * [`Section.get()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.get)


            * [`Section.get_keyword()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.get_keyword)


            * [`Section.get_section()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.get_section)


            * [`Section.get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.get_str)


            * [`Section.get_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.get_string)


            * [`Section.inc()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.inc)


            * [`Section.insert()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.insert)


            * [`Section.safeset()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.safeset)


            * [`Section.set()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.set)


            * [`Section.setitem()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.setitem)


            * [`Section.silence()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.silence)


            * [`Section.unset()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.unset)


            * [`Section.update()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.update)


            * [`Section.verbosity()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Section.verbosity)


        * [`SectionList`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.SectionList)


            * [`SectionList._get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.SectionList._get_str)


            * [`SectionList.append()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.SectionList.append)


            * [`SectionList.extend()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.SectionList.extend)


            * [`SectionList.get()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.SectionList.get)


            * [`SectionList.get_str()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.SectionList.get_str)


            * [`SectionList.get_string()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.SectionList.get_string)


            * [`SectionList.verbosity()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.SectionList.verbosity)


        * [`Smear`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Smear)


        * [`Subsys`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Subsys)


        * [`V_Hartree_Cube`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.V_Hartree_Cube)


        * [`Xc_Functional`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Xc_Functional)


    * [pymatgen.io.cp2k.outputs module](pymatgen.io.cp2k.md#module-pymatgen.io.cp2k.outputs)


        * [`Cp2kOutput`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput)


            * [`Cp2kOutput._gauss_smear()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput._gauss_smear)


            * [`Cp2kOutput.as_dict()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.as_dict)


            * [`Cp2kOutput.band_structure`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.band_structure)


            * [`Cp2kOutput.calculation_type`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.calculation_type)


            * [`Cp2kOutput.charge`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.charge)


            * [`Cp2kOutput.complete_dos`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.complete_dos)


            * [`Cp2kOutput.completed`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.completed)


            * [`Cp2kOutput.convergence()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.convergence)


            * [`Cp2kOutput.cp2k_version`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.cp2k_version)


            * [`Cp2kOutput.is_hubbard`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.is_hubbard)


            * [`Cp2kOutput.is_metal`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.is_metal)


            * [`Cp2kOutput.is_molecule`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.is_molecule)


            * [`Cp2kOutput.multiplicity`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.multiplicity)


            * [`Cp2kOutput.num_warnings`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.num_warnings)


            * [`Cp2kOutput.parse_atomic_kind_info()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_atomic_kind_info)


            * [`Cp2kOutput.parse_bandstructure()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_bandstructure)


            * [`Cp2kOutput.parse_cell_params()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_cell_params)


            * [`Cp2kOutput.parse_chi_tensor()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_chi_tensor)


            * [`Cp2kOutput.parse_cp2k_params()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_cp2k_params)


            * [`Cp2kOutput.parse_dft_params()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_dft_params)


            * [`Cp2kOutput.parse_dos()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_dos)


            * [`Cp2kOutput.parse_energies()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_energies)


            * [`Cp2kOutput.parse_files()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_files)


            * [`Cp2kOutput.parse_forces()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_forces)


            * [`Cp2kOutput.parse_global_params()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_global_params)


            * [`Cp2kOutput.parse_gtensor()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_gtensor)


            * [`Cp2kOutput.parse_hirshfeld()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_hirshfeld)


            * [`Cp2kOutput.parse_homo_lumo()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_homo_lumo)


            * [`Cp2kOutput.parse_hyperfine()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_hyperfine)


            * [`Cp2kOutput.parse_initial_structure()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_initial_structure)


            * [`Cp2kOutput.parse_input()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_input)


            * [`Cp2kOutput.parse_ionic_steps()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_ionic_steps)


            * [`Cp2kOutput.parse_mo_eigenvalues()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_mo_eigenvalues)


            * [`Cp2kOutput.parse_mulliken()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_mulliken)


            * [`Cp2kOutput.parse_nmr_shift()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_nmr_shift)


            * [`Cp2kOutput.parse_opt_steps()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_opt_steps)


            * [`Cp2kOutput.parse_overlap_condition()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_overlap_condition)


            * [`Cp2kOutput.parse_plus_u_params()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_plus_u_params)


            * [`Cp2kOutput.parse_qs_params()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_qs_params)


            * [`Cp2kOutput.parse_raman()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_raman)


            * [`Cp2kOutput.parse_scf_opt()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_scf_opt)


            * [`Cp2kOutput.parse_scf_params()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_scf_params)


            * [`Cp2kOutput.parse_stresses()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_stresses)


            * [`Cp2kOutput.parse_structures()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_structures)


            * [`Cp2kOutput.parse_tddfpt()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_tddfpt)


            * [`Cp2kOutput.parse_timing()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_timing)


            * [`Cp2kOutput.parse_total_numbers()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_total_numbers)


            * [`Cp2kOutput.project_name`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.project_name)


            * [`Cp2kOutput.ran_successfully()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.ran_successfully)


            * [`Cp2kOutput.read_pattern()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.read_pattern)


            * [`Cp2kOutput.read_table_pattern()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.read_table_pattern)


            * [`Cp2kOutput.run_type`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.run_type)


            * [`Cp2kOutput.spin_polarized`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.Cp2kOutput.spin_polarized)


        * [`parse_dos()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.parse_dos)


        * [`parse_energy_file()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.parse_energy_file)


        * [`parse_pdos()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.outputs.parse_pdos)


    * [pymatgen.io.cp2k.sets module](pymatgen.io.cp2k.md#module-pymatgen.io.cp2k.sets)


        * [`CellOptSet`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.CellOptSet)


        * [`Cp2kValidationError`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.Cp2kValidationError)


            * [`Cp2kValidationError.CP2K_VERSION`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.Cp2kValidationError.CP2K_VERSION)


        * [`DftSet`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet)


            * [`DftSet.activate_epr()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_epr)


            * [`DftSet.activate_fast_minimization()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_fast_minimization)


            * [`DftSet.activate_hybrid()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_hybrid)


            * [`DftSet.activate_hyperfine()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_hyperfine)


            * [`DftSet.activate_localize()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_localize)


            * [`DftSet.activate_motion()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_motion)


            * [`DftSet.activate_nmr()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_nmr)


            * [`DftSet.activate_nonperiodic()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_nonperiodic)


            * [`DftSet.activate_polar()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_polar)


            * [`DftSet.activate_robust_minimization()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_robust_minimization)


            * [`DftSet.activate_spinspin()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_spinspin)


            * [`DftSet.activate_tddfpt()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_tddfpt)


            * [`DftSet.activate_vdw_potential()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_vdw_potential)


            * [`DftSet.activate_very_strict_minimization()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.activate_very_strict_minimization)


            * [`DftSet.create_subsys()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.create_subsys)


            * [`DftSet.get_basis_and_potential()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.get_basis_and_potential)


            * [`DftSet.get_cutoff_from_basis()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.get_cutoff_from_basis)


            * [`DftSet.get_xc_functionals()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.get_xc_functionals)


            * [`DftSet.modify_dft_print_iters()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.modify_dft_print_iters)


            * [`DftSet.print_bandstructure()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_bandstructure)


            * [`DftSet.print_dos()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_dos)


            * [`DftSet.print_e_density()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_e_density)


            * [`DftSet.print_forces()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_forces)


            * [`DftSet.print_hirshfeld()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_hirshfeld)


            * [`DftSet.print_ldos()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_ldos)


            * [`DftSet.print_mo()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_mo)


            * [`DftSet.print_mo_cubes()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_mo_cubes)


            * [`DftSet.print_mulliken()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_mulliken)


            * [`DftSet.print_pdos()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_pdos)


            * [`DftSet.print_v_hartree()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.print_v_hartree)


            * [`DftSet.set_charge()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.set_charge)


            * [`DftSet.validate()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.validate)


            * [`DftSet.write_basis_set_file()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.write_basis_set_file)


            * [`DftSet.write_potential_file()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.DftSet.write_potential_file)


        * [`HybridCellOptSet`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.HybridCellOptSet)


        * [`HybridRelaxSet`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.HybridRelaxSet)


        * [`HybridStaticSet`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.HybridStaticSet)


        * [`RelaxSet`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.RelaxSet)


        * [`StaticSet`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.sets.StaticSet)


    * [pymatgen.io.cp2k.utils module](pymatgen.io.cp2k.md#module-pymatgen.io.cp2k.utils)


        * [`chunk()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.utils.chunk)


        * [`get_truncated_coulomb_cutoff()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.utils.get_truncated_coulomb_cutoff)


        * [`get_unique_site_indices()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.utils.get_unique_site_indices)


        * [`natural_keys()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.utils.natural_keys)


        * [`postprocessor()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.utils.postprocessor)


        * [`preprocessor()`](pymatgen.io.cp2k.md#pymatgen.io.cp2k.utils.preprocessor)


* [pymatgen.io.exciting package](pymatgen.io.exciting.md)




    * [pymatgen.io.exciting.inputs module](pymatgen.io.exciting.md#module-pymatgen.io.exciting.inputs)


        * [`ExcitingInput`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput)


            * [`ExcitingInput.structure`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.structure)


            * [`ExcitingInput.title`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.title)


            * [`ExcitingInput.lockxyz`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.lockxyz)


            * [`ExcitingInput._dicttoxml()`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput._dicttoxml)


            * [`ExcitingInput._indent()`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput._indent)


            * [`ExcitingInput.bohr2ang`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.bohr2ang)


            * [`ExcitingInput.from_file()`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.from_file)


            * [`ExcitingInput.from_str()`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.from_str)


            * [`ExcitingInput.from_string()`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.from_string)


            * [`ExcitingInput.lockxyz`](pymatgen.io.exciting.md#id0)


            * [`ExcitingInput.write_etree()`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.write_etree)


            * [`ExcitingInput.write_file()`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.write_file)


            * [`ExcitingInput.write_string()`](pymatgen.io.exciting.md#pymatgen.io.exciting.inputs.ExcitingInput.write_string)


* [pymatgen.io.feff package](pymatgen.io.feff.md)




    * [pymatgen.io.feff.inputs module](pymatgen.io.feff.md#module-pymatgen.io.feff.inputs)


        * [`Atoms`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Atoms)


            * [`Atoms._set_cluster()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Atoms._set_cluster)


            * [`Atoms.atoms_string_from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Atoms.atoms_string_from_file)


            * [`Atoms.cluster`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Atoms.cluster)


            * [`Atoms.cluster_from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Atoms.cluster_from_file)


            * [`Atoms.get_lines()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Atoms.get_lines)


            * [`Atoms.write_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Atoms.write_file)


        * [`FeffParseError`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.FeffParseError)


        * [`Header`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Header)


            * [`Header.formula`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Header.formula)


            * [`Header.from_cif_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Header.from_cif_file)


            * [`Header.from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Header.from_file)


            * [`Header.from_str()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Header.from_str)


            * [`Header.from_string()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Header.from_string)


            * [`Header.header_string_from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Header.header_string_from_file)


            * [`Header.structure_symmetry`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Header.structure_symmetry)


            * [`Header.write_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Header.write_file)


        * [`Paths`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Paths)


            * [`Paths.write_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Paths.write_file)


        * [`Potential`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Potential)


            * [`Potential.pot_dict_from_string()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Potential.pot_dict_from_string)


            * [`Potential.pot_string_from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Potential.pot_string_from_file)


            * [`Potential.write_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Potential.write_file)


        * [`Tags`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags)


            * [`Tags._stringify_val()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags._stringify_val)


            * [`Tags.as_dict()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags.as_dict)


            * [`Tags.diff()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags.diff)


            * [`Tags.from_dict()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags.from_dict)


            * [`Tags.from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags.from_file)


            * [`Tags.get_str()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags.get_str)


            * [`Tags.get_string()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags.get_string)


            * [`Tags.proc_val()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags.proc_val)


            * [`Tags.write_file()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Tags.write_file)


        * [`get_absorbing_atom_symbol_index()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.get_absorbing_atom_symbol_index)


        * [`get_atom_map()`](pymatgen.io.feff.md#pymatgen.io.feff.inputs.get_atom_map)


    * [pymatgen.io.feff.outputs module](pymatgen.io.feff.md#module-pymatgen.io.feff.outputs)


        * [`Eels`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Eels)


            * [`Eels.as_dict()`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Eels.as_dict)


            * [`Eels.atomic_background`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Eels.atomic_background)


            * [`Eels.energies`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Eels.energies)


            * [`Eels.fine_structure`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Eels.fine_structure)


            * [`Eels.from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Eels.from_file)


            * [`Eels.total_spectrum`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Eels.total_spectrum)


        * [`LDos`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.LDos)


            * [`LDos.charge_transfer_from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.LDos.charge_transfer_from_file)


            * [`LDos.charge_transfer_to_string()`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.LDos.charge_transfer_to_string)


            * [`LDos.from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.LDos.from_file)


        * [`Xmu`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu)


            * [`Xmu.as_dict()`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.as_dict)


            * [`Xmu.calc`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.calc)


            * [`Xmu.chi`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.chi)


            * [`Xmu.e_fermi`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.e_fermi)


            * [`Xmu.edge`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.edge)


            * [`Xmu.energies`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.energies)


            * [`Xmu.from_file()`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.from_file)


            * [`Xmu.material_formula`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.material_formula)


            * [`Xmu.mu`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.mu)


            * [`Xmu.mu0`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.mu0)


            * [`Xmu.relative_energies`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.relative_energies)


            * [`Xmu.source`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.source)


            * [`Xmu.wavenumber`](pymatgen.io.feff.md#pymatgen.io.feff.outputs.Xmu.wavenumber)


    * [pymatgen.io.feff.sets module](pymatgen.io.feff.md#module-pymatgen.io.feff.sets)


        * [`AbstractFeffInputSet`](pymatgen.io.feff.md#pymatgen.io.feff.sets.AbstractFeffInputSet)


            * [`AbstractFeffInputSet._abc_impl`](pymatgen.io.feff.md#pymatgen.io.feff.sets.AbstractFeffInputSet._abc_impl)


            * [`AbstractFeffInputSet.all_input()`](pymatgen.io.feff.md#pymatgen.io.feff.sets.AbstractFeffInputSet.all_input)


            * [`AbstractFeffInputSet.atoms`](pymatgen.io.feff.md#pymatgen.io.feff.sets.AbstractFeffInputSet.atoms)


            * [`AbstractFeffInputSet.header()`](pymatgen.io.feff.md#pymatgen.io.feff.sets.AbstractFeffInputSet.header)


            * [`AbstractFeffInputSet.potential`](pymatgen.io.feff.md#pymatgen.io.feff.sets.AbstractFeffInputSet.potential)


            * [`AbstractFeffInputSet.tags`](pymatgen.io.feff.md#pymatgen.io.feff.sets.AbstractFeffInputSet.tags)


            * [`AbstractFeffInputSet.write_input()`](pymatgen.io.feff.md#pymatgen.io.feff.sets.AbstractFeffInputSet.write_input)


        * [`FEFFDictSet`](pymatgen.io.feff.md#pymatgen.io.feff.sets.FEFFDictSet)


            * [`FEFFDictSet._abc_impl`](pymatgen.io.feff.md#pymatgen.io.feff.sets.FEFFDictSet._abc_impl)


            * [`FEFFDictSet.atoms`](pymatgen.io.feff.md#pymatgen.io.feff.sets.FEFFDictSet.atoms)


            * [`FEFFDictSet.from_directory()`](pymatgen.io.feff.md#pymatgen.io.feff.sets.FEFFDictSet.from_directory)


            * [`FEFFDictSet.header()`](pymatgen.io.feff.md#pymatgen.io.feff.sets.FEFFDictSet.header)


            * [`FEFFDictSet.potential`](pymatgen.io.feff.md#pymatgen.io.feff.sets.FEFFDictSet.potential)


            * [`FEFFDictSet.tags`](pymatgen.io.feff.md#pymatgen.io.feff.sets.FEFFDictSet.tags)


        * [`MPEELSDictSet`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPEELSDictSet)


            * [`MPEELSDictSet._abc_impl`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPEELSDictSet._abc_impl)


        * [`MPELNESSet`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPELNESSet)


            * [`MPELNESSet.CONFIG`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPELNESSet.CONFIG)


            * [`MPELNESSet._abc_impl`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPELNESSet._abc_impl)


        * [`MPEXAFSSet`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPEXAFSSet)


            * [`MPEXAFSSet.CONFIG`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPEXAFSSet.CONFIG)


            * [`MPEXAFSSet._abc_impl`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPEXAFSSet._abc_impl)


        * [`MPEXELFSSet`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPEXELFSSet)


            * [`MPEXELFSSet.CONFIG`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPEXELFSSet.CONFIG)


            * [`MPEXELFSSet._abc_impl`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPEXELFSSet._abc_impl)


        * [`MPXANESSet`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPXANESSet)


            * [`MPXANESSet.CONFIG`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPXANESSet.CONFIG)


            * [`MPXANESSet._abc_impl`](pymatgen.io.feff.md#pymatgen.io.feff.sets.MPXANESSet._abc_impl)


* [pymatgen.io.lammps package](pymatgen.io.lammps.md)




    * [pymatgen.io.lammps.data module](pymatgen.io.lammps.md#module-pymatgen.io.lammps.data)


        * [`CombinedData`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData)


            * [`CombinedData.as_lammpsdata()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.as_lammpsdata)


            * [`CombinedData.disassemble()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.disassemble)


            * [`CombinedData.from_ff_and_topologies()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.from_ff_and_topologies)


            * [`CombinedData.from_files()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.from_files)


            * [`CombinedData.from_lammpsdata()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.from_lammpsdata)


            * [`CombinedData.from_structure()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.from_structure)


            * [`CombinedData.get_str()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.get_str)


            * [`CombinedData.get_string()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.get_string)


            * [`CombinedData.parse_xyz()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.parse_xyz)


            * [`CombinedData.structure`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.CombinedData.structure)


        * [`ForceField`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField)


            * [`ForceField.masses`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField.masses)


            * [`ForceField.force_field`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField.force_field)


            * [`ForceField.maps`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField.maps)


            * [`ForceField._is_valid()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField._is_valid)


            * [`ForceField._process_nonbond()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField._process_nonbond)


            * [`ForceField._process_topo()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField._process_topo)


            * [`ForceField.from_dict()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField.from_dict)


            * [`ForceField.from_file()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField.from_file)


            * [`ForceField.to_file()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.ForceField.to_file)


        * [`LammpsBox`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsBox)


            * [`LammpsBox.get_box_shift()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsBox.get_box_shift)


            * [`LammpsBox.get_str()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsBox.get_str)


            * [`LammpsBox.get_string()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsBox.get_string)


            * [`LammpsBox.to_lattice()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsBox.to_lattice)


            * [`LammpsBox.volume`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsBox.volume)


        * [`LammpsData`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData)


            * [`LammpsData.disassemble()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.disassemble)


            * [`LammpsData.from_ff_and_topologies()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.from_ff_and_topologies)


            * [`LammpsData.from_file()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.from_file)


            * [`LammpsData.from_structure()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.from_structure)


            * [`LammpsData.get_str()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.get_str)


            * [`LammpsData.get_string()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.get_string)


            * [`LammpsData.set_charge_atom()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.set_charge_atom)


            * [`LammpsData.set_charge_atom_type()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.set_charge_atom_type)


            * [`LammpsData.structure`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.structure)


            * [`LammpsData.write_file()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.LammpsData.write_file)


        * [`Topology`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.Topology)


            * [`Topology.from_bonding()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.Topology.from_bonding)


        * [`lattice_2_lmpbox()`](pymatgen.io.lammps.md#pymatgen.io.lammps.data.lattice_2_lmpbox)


    * [pymatgen.io.lammps.generators module](pymatgen.io.lammps.md#module-pymatgen.io.lammps.generators)


        * [`BaseLammpsGenerator`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.BaseLammpsGenerator)


            * [`BaseLammpsGenerator.calc_type`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.calc_type)


            * [`BaseLammpsGenerator.get_input_set()`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.get_input_set)


            * [`BaseLammpsGenerator.keep_stages`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.keep_stages)


            * [`BaseLammpsGenerator.settings`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.settings)


            * [`BaseLammpsGenerator.template`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.template)


        * [`LammpsMinimization`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.LammpsMinimization)


            * [`LammpsMinimization.atom_style`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.LammpsMinimization.atom_style)


            * [`LammpsMinimization.boundary`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.LammpsMinimization.boundary)


            * [`LammpsMinimization.dimension`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.LammpsMinimization.dimension)


            * [`LammpsMinimization.force_field`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.LammpsMinimization.force_field)


            * [`LammpsMinimization.read_data`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.LammpsMinimization.read_data)


            * [`LammpsMinimization.settings`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.LammpsMinimization.settings)


            * [`LammpsMinimization.template`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.LammpsMinimization.template)


            * [`LammpsMinimization.units`](pymatgen.io.lammps.md#pymatgen.io.lammps.generators.LammpsMinimization.units)


    * [pymatgen.io.lammps.inputs module](pymatgen.io.lammps.md#module-pymatgen.io.lammps.inputs)


        * [`LammpsInputFile`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile)


            * [`LammpsInputFile._add_command()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile._add_command)


            * [`LammpsInputFile._add_comment()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile._add_comment)


            * [`LammpsInputFile._check_stage_format()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile._check_stage_format)


            * [`LammpsInputFile._clean_lines()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile._clean_lines)


            * [`LammpsInputFile._get_blocks()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile._get_blocks)


            * [`LammpsInputFile._initialize_stage()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile._initialize_stage)


            * [`LammpsInputFile.add_commands()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.add_commands)


            * [`LammpsInputFile.add_stage()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.add_stage)


            * [`LammpsInputFile.append()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.append)


            * [`LammpsInputFile.contains_command()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.contains_command)


            * [`LammpsInputFile.from_file()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.from_file)


            * [`LammpsInputFile.from_str()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.from_str)


            * [`LammpsInputFile.from_string()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.from_string)


            * [`LammpsInputFile.get_args()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.get_args)


            * [`LammpsInputFile.get_str()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.get_str)


            * [`LammpsInputFile.get_string()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.get_string)


            * [`LammpsInputFile.merge_stages()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.merge_stages)


            * [`LammpsInputFile.ncomments`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.ncomments)


            * [`LammpsInputFile.nstages`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.nstages)


            * [`LammpsInputFile.remove_command()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.remove_command)


            * [`LammpsInputFile.remove_stage()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.remove_stage)


            * [`LammpsInputFile.rename_stage()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.rename_stage)


            * [`LammpsInputFile.set_args()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.set_args)


            * [`LammpsInputFile.stages_names`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.stages_names)


            * [`LammpsInputFile.write_file()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsInputFile.write_file)


        * [`LammpsRun`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsRun)


            * [`LammpsRun.md()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsRun.md)


            * [`LammpsRun.template_dir`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsRun.template_dir)


            * [`LammpsRun.write_inputs()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsRun.write_inputs)


        * [`LammpsTemplateGen`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsTemplateGen)


            * [`LammpsTemplateGen.get_input_set()`](pymatgen.io.lammps.md#pymatgen.io.lammps.inputs.LammpsTemplateGen.get_input_set)


    * [pymatgen.io.lammps.outputs module](pymatgen.io.lammps.md#module-pymatgen.io.lammps.outputs)


        * [`LammpsDump`](pymatgen.io.lammps.md#pymatgen.io.lammps.outputs.LammpsDump)


            * [`LammpsDump.as_dict()`](pymatgen.io.lammps.md#pymatgen.io.lammps.outputs.LammpsDump.as_dict)


            * [`LammpsDump.from_dict()`](pymatgen.io.lammps.md#pymatgen.io.lammps.outputs.LammpsDump.from_dict)


            * [`LammpsDump.from_str()`](pymatgen.io.lammps.md#pymatgen.io.lammps.outputs.LammpsDump.from_str)


            * [`LammpsDump.from_string()`](pymatgen.io.lammps.md#pymatgen.io.lammps.outputs.LammpsDump.from_string)


        * [`parse_lammps_dumps()`](pymatgen.io.lammps.md#pymatgen.io.lammps.outputs.parse_lammps_dumps)


        * [`parse_lammps_log()`](pymatgen.io.lammps.md#pymatgen.io.lammps.outputs.parse_lammps_log)


    * [pymatgen.io.lammps.sets module](pymatgen.io.lammps.md#module-pymatgen.io.lammps.sets)


        * [`LammpsInputSet`](pymatgen.io.lammps.md#pymatgen.io.lammps.sets.LammpsInputSet)


            * [`LammpsInputSet._abc_impl`](pymatgen.io.lammps.md#pymatgen.io.lammps.sets.LammpsInputSet._abc_impl)


            * [`LammpsInputSet.from_directory()`](pymatgen.io.lammps.md#pymatgen.io.lammps.sets.LammpsInputSet.from_directory)


            * [`LammpsInputSet.validate()`](pymatgen.io.lammps.md#pymatgen.io.lammps.sets.LammpsInputSet.validate)


    * [pymatgen.io.lammps.utils module](pymatgen.io.lammps.md#module-pymatgen.io.lammps.utils)


        * [`LammpsRunner`](pymatgen.io.lammps.md#pymatgen.io.lammps.utils.LammpsRunner)


            * [`LammpsRunner.run()`](pymatgen.io.lammps.md#pymatgen.io.lammps.utils.LammpsRunner.run)


        * [`Polymer`](pymatgen.io.lammps.md#pymatgen.io.lammps.utils.Polymer)


            * [`Polymer._add_monomer()`](pymatgen.io.lammps.md#pymatgen.io.lammps.utils.Polymer._add_monomer)


            * [`Polymer._align_monomer()`](pymatgen.io.lammps.md#pymatgen.io.lammps.utils.Polymer._align_monomer)


            * [`Polymer._create()`](pymatgen.io.lammps.md#pymatgen.io.lammps.utils.Polymer._create)


            * [`Polymer._next_move_direction()`](pymatgen.io.lammps.md#pymatgen.io.lammps.utils.Polymer._next_move_direction)


* [pymatgen.io.lobster package](pymatgen.io.lobster.md)




    * [pymatgen.io.lobster.inputs module](pymatgen.io.lobster.md#module-pymatgen.io.lobster.inputs)


        * [`Lobsterin`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin)


            * [`Lobsterin.AVAILABLE_KEYWORDS`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.AVAILABLE_KEYWORDS)


            * [`Lobsterin.BOOLEAN_KEYWORDS`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.BOOLEAN_KEYWORDS)


            * [`Lobsterin.FLOAT_KEYWORDS`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.FLOAT_KEYWORDS)


            * [`Lobsterin.LISTKEYWORDS`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.LISTKEYWORDS)


            * [`Lobsterin.STRING_KEYWORDS`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.STRING_KEYWORDS)


            * [`Lobsterin._get_nbands()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin._get_nbands)


            * [`Lobsterin._get_potcar_symbols()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin._get_potcar_symbols)


            * [`Lobsterin.as_dict()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.as_dict)


            * [`Lobsterin.diff()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.diff)


            * [`Lobsterin.from_dict()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.from_dict)


            * [`Lobsterin.from_file()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.from_file)


            * [`Lobsterin.get_all_possible_basis_functions()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.get_all_possible_basis_functions)


            * [`Lobsterin.get_basis()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.get_basis)


            * [`Lobsterin.standard_calculations_from_vasp_files()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.standard_calculations_from_vasp_files)


            * [`Lobsterin.write_INCAR()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.write_INCAR)


            * [`Lobsterin.write_KPOINTS()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.write_KPOINTS)


            * [`Lobsterin.write_POSCAR_with_standard_primitive()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.write_POSCAR_with_standard_primitive)


            * [`Lobsterin.write_lobsterin()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.Lobsterin.write_lobsterin)


        * [`get_all_possible_basis_combinations()`](pymatgen.io.lobster.md#pymatgen.io.lobster.inputs.get_all_possible_basis_combinations)


    * [pymatgen.io.lobster.lobsterenv module](pymatgen.io.lobster.md#module-pymatgen.io.lobster.lobsterenv)


        * [`ICOHPNeighborsInfo`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo)


            * [`ICOHPNeighborsInfo._asdict()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo._asdict)


            * [`ICOHPNeighborsInfo._field_defaults`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo._field_defaults)


            * [`ICOHPNeighborsInfo._fields`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo._fields)


            * [`ICOHPNeighborsInfo._make()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo._make)


            * [`ICOHPNeighborsInfo._replace()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo._replace)


            * [`ICOHPNeighborsInfo.atoms`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.atoms)


            * [`ICOHPNeighborsInfo.central_isites`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.central_isites)


            * [`ICOHPNeighborsInfo.labels`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.labels)


            * [`ICOHPNeighborsInfo.list_icohps`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.list_icohps)


            * [`ICOHPNeighborsInfo.n_bonds`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.n_bonds)


            * [`ICOHPNeighborsInfo.total_icohp`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.total_icohp)


        * [`LobsterLightStructureEnvironments`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterLightStructureEnvironments)


            * [`LobsterLightStructureEnvironments.as_dict()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterLightStructureEnvironments.as_dict)


            * [`LobsterLightStructureEnvironments.from_Lobster()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterLightStructureEnvironments.from_Lobster)


            * [`LobsterLightStructureEnvironments.uniquely_determines_coordination_environments`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterLightStructureEnvironments.uniquely_determines_coordination_environments)


        * [`LobsterNeighbors`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors)


            * [`LobsterNeighbors._adapt_extremum_to_add_cond()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._adapt_extremum_to_add_cond)


            * [`LobsterNeighbors._determine_unit_cell()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._determine_unit_cell)


            * [`LobsterNeighbors._evaluate_ce()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._evaluate_ce)


            * [`LobsterNeighbors._find_environments()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._find_environments)


            * [`LobsterNeighbors._find_relevant_atoms_additional_condition()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._find_relevant_atoms_additional_condition)


            * [`LobsterNeighbors._get_atomnumber()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._get_atomnumber)


            * [`LobsterNeighbors._get_icohps()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._get_icohps)


            * [`LobsterNeighbors._get_limit_from_extremum()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._get_limit_from_extremum)


            * [`LobsterNeighbors._get_plot_label()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._get_plot_label)


            * [`LobsterNeighbors._split_string()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors._split_string)


            * [`LobsterNeighbors.anion_types`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.anion_types)


            * [`LobsterNeighbors.get_anion_types()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_anion_types)


            * [`LobsterNeighbors.get_info_cohps_to_neighbors()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_info_cohps_to_neighbors)


            * [`LobsterNeighbors.get_info_icohps_between_neighbors()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_info_icohps_between_neighbors)


            * [`LobsterNeighbors.get_info_icohps_to_neighbors()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_info_icohps_to_neighbors)


            * [`LobsterNeighbors.get_light_structure_environment()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_light_structure_environment)


            * [`LobsterNeighbors.get_nn_info()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_nn_info)


            * [`LobsterNeighbors.molecules_allowed`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.molecules_allowed)


            * [`LobsterNeighbors.plot_cohps_of_neighbors()`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.plot_cohps_of_neighbors)


            * [`LobsterNeighbors.structures_allowed`](pymatgen.io.lobster.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.structures_allowed)


    * [pymatgen.io.lobster.outputs module](pymatgen.io.lobster.md#module-pymatgen.io.lobster.outputs)


        * [`Bandoverlaps`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Bandoverlaps)


            * [`Bandoverlaps.maxDeviation`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Bandoverlaps.maxDeviation)


            * [`Bandoverlaps._read()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Bandoverlaps._read)


            * [`Bandoverlaps.has_good_quality_check_occupied_bands()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Bandoverlaps.has_good_quality_check_occupied_bands)


            * [`Bandoverlaps.has_good_quality_maxDeviation()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Bandoverlaps.has_good_quality_maxDeviation)


        * [`Charge`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Charge)


            * [`Charge.atomlist`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Charge.atomlist)


            * [`Charge.types`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Charge.types)


            * [`Charge.Mulliken`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Charge.Mulliken)


            * [`Charge.Loewdin`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Charge.Loewdin)


            * [`Charge.num_atoms`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Charge.num_atoms)


            * [`Charge.get_structure_with_charges()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Charge.get_structure_with_charges)


        * [`Cohpcar`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Cohpcar)


            * [`Cohpcar.cohp_data`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Cohpcar.cohp_data)


            * [`Cohpcar.efermi`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Cohpcar.efermi)


            * [`Cohpcar.energies`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Cohpcar.energies)


            * [`Cohpcar.is_spin_polarized`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Cohpcar.is_spin_polarized)


            * [`Cohpcar.orb_cohp`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Cohpcar.orb_cohp)


            * [`Cohpcar._get_bond_data()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Cohpcar._get_bond_data)


        * [`Doscar`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Doscar)


            * [`Doscar.completedos`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Doscar.completedos)


            * [`Doscar.pdos`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Doscar.pdos)


            * [`Doscar.tdos`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Doscar.tdos)


            * [`Doscar.energies`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Doscar.energies)


            * [`Doscar.tdensities`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Doscar.tdensities)


            * [`Doscar.itdensities`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Doscar.itdensities)


            * [`Doscar.is_spin_polarized`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Doscar.is_spin_polarized)


            * [`Doscar._parse_doscar()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Doscar._parse_doscar)


            * [`Doscar.completedos`](pymatgen.io.lobster.md#id0)


            * [`Doscar.energies`](pymatgen.io.lobster.md#id1)


            * [`Doscar.is_spin_polarized`](pymatgen.io.lobster.md#id2)


            * [`Doscar.itdensities`](pymatgen.io.lobster.md#id3)


            * [`Doscar.pdos`](pymatgen.io.lobster.md#id4)


            * [`Doscar.tdensities`](pymatgen.io.lobster.md#id5)


            * [`Doscar.tdos`](pymatgen.io.lobster.md#id6)


        * [`Fatband`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband)


            * [`Fatband.efermi`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.efermi)


            * [`Fatband.eigenvals`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.eigenvals)


            * [`Fatband.is_spin_polarized`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.is_spin_polarized)


            * [`Fatband.kpoints_array`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.kpoints_array)


            * [`Fatband.label_dict`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.label_dict)


            * [`Fatband.lattice`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.lattice)


            * [`Fatband.nbands`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.nbands)


            * [`Fatband.p_eigenvals`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.p_eigenvals)


            * [`Fatband.structure`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.structure)


            * [`Fatband.get_bandstructure()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Fatband.get_bandstructure)


        * [`Grosspop`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Grosspop)


            * [`Grosspop.list_dict_grosspop`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Grosspop.list_dict_grosspop)


            * [`Grosspop.get_structure_with_total_grosspop()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Grosspop.get_structure_with_total_grosspop)


        * [`Icohplist`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Icohplist)


            * [`Icohplist.are_coops`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Icohplist.are_coops)


            * [`Icohplist.is_spin_polarized`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Icohplist.is_spin_polarized)


            * [`Icohplist.Icohplist`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Icohplist.Icohplist)


            * [`Icohplist.IcohpCollection`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Icohplist.IcohpCollection)


            * [`Icohplist.icohpcollection`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Icohplist.icohpcollection)


            * [`Icohplist.icohplist`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Icohplist.icohplist)


        * [`Lobsterout`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout)


            * [`Lobsterout.basis_functions`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.basis_functions)


            * [`Lobsterout.basis_type`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.basis_type)


            * [`Lobsterout.charge_spilling`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.charge_spilling)


            * [`Lobsterout.dft_program`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.dft_program)


            * [`Lobsterout.elements`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.elements)


            * [`Lobsterout.has_charge`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_charge)


            * [`Lobsterout.has_cohpcar`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_cohpcar)


            * [`Lobsterout.has_madelung`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_madelung)


            * [`Lobsterout.has_coopcar`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_coopcar)


            * [`Lobsterout.has_cobicar`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_cobicar)


            * [`Lobsterout.has_doscar`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_doscar)


            * [`Lobsterout.has_doscar_lso`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_doscar_lso)


            * [`Lobsterout.has_projection`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_projection)


            * [`Lobsterout.has_bandoverlaps`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_bandoverlaps)


            * [`Lobsterout.has_density_of_energies`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_density_of_energies)


            * [`Lobsterout.has_fatbands`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_fatbands)


            * [`Lobsterout.has_grosspopulation`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.has_grosspopulation)


            * [`Lobsterout.info_lines`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.info_lines)


            * [`Lobsterout.info_orthonormalization`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.info_orthonormalization)


            * [`Lobsterout.is_restart_from_projection`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.is_restart_from_projection)


            * [`Lobsterout.lobster_version`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.lobster_version)


            * [`Lobsterout.number_of_spins`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.number_of_spins)


            * [`Lobsterout.number_of_threads`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.number_of_threads)


            * [`Lobsterout.timing`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.timing)


            * [`Lobsterout.total_spilling`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.total_spilling)


            * [`Lobsterout.warning_lines`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.warning_lines)


            * [`Lobsterout._get_all_info_lines()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_all_info_lines)


            * [`Lobsterout._get_all_warning_lines()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_all_warning_lines)


            * [`Lobsterout._get_dft_program()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_dft_program)


            * [`Lobsterout._get_elements_basistype_basisfunctions()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_elements_basistype_basisfunctions)


            * [`Lobsterout._get_lobster_version()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_lobster_version)


            * [`Lobsterout._get_number_of_spins()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_number_of_spins)


            * [`Lobsterout._get_spillings()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_spillings)


            * [`Lobsterout._get_threads()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_threads)


            * [`Lobsterout._get_timing()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_timing)


            * [`Lobsterout._get_warning_orthonormalization()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._get_warning_orthonormalization)


            * [`Lobsterout._has_fatband()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout._has_fatband)


            * [`Lobsterout.get_doc()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Lobsterout.get_doc)


        * [`MadelungEnergies`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.MadelungEnergies)


            * [`MadelungEnergies.madelungenergies_Mulliken`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.MadelungEnergies.madelungenergies_Mulliken)


            * [`MadelungEnergies.madelungenergies_Loewdin`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.MadelungEnergies.madelungenergies_Loewdin)


            * [`MadelungEnergies.ewald_splitting`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.MadelungEnergies.ewald_splitting)


        * [`SitePotential`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential)


            * [`SitePotential.atomlist`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential.atomlist)


            * [`SitePotential.types`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential.types)


            * [`SitePotential.num_atoms`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential.num_atoms)


            * [`SitePotential.sitepotentials_Mulliken`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential.sitepotentials_Mulliken)


            * [`SitePotential.sitepotentials_Loewdin`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential.sitepotentials_Loewdin)


            * [`SitePotential.madelung_Mulliken`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential.madelung_Mulliken)


            * [`SitePotential.madelung_Loewdin`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential.madelung_Loewdin)


            * [`SitePotential.ewald_splitting`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential.ewald_splitting)


            * [`SitePotential.get_structure_with_site_potentials()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.SitePotential.get_structure_with_site_potentials)


        * [`Wavefunction`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction)


            * [`Wavefunction.grid`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.grid)


            * [`Wavefunction.points`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.points)


            * [`Wavefunction.real`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.real)


            * [`Wavefunction.imaginary`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.imaginary)


            * [`Wavefunction.distance`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.distance)


            * [`Wavefunction._parse_file()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction._parse_file)


            * [`Wavefunction.get_volumetricdata_density()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.get_volumetricdata_density)


            * [`Wavefunction.get_volumetricdata_imaginary()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.get_volumetricdata_imaginary)


            * [`Wavefunction.get_volumetricdata_real()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.get_volumetricdata_real)


            * [`Wavefunction.set_volumetric_data()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.set_volumetric_data)


            * [`Wavefunction.write_file()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.Wavefunction.write_file)


        * [`get_orb_from_str()`](pymatgen.io.lobster.md#pymatgen.io.lobster.outputs.get_orb_from_str)


* [pymatgen.io.qchem package](pymatgen.io.qchem.md)




    * [pymatgen.io.qchem.inputs module](pymatgen.io.qchem.md#module-pymatgen.io.qchem.inputs)


        * [`QCInput`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput)


            * [`QCInput.almo_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.almo_template)


            * [`QCInput.cdft_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.cdft_template)


            * [`QCInput.find_sections()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.find_sections)


            * [`QCInput.from_file()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.from_file)


            * [`QCInput.from_multi_jobs_file()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.from_multi_jobs_file)


            * [`QCInput.from_str()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.from_str)


            * [`QCInput.geom_opt_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.geom_opt_template)


            * [`QCInput.get_str()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.get_str)


            * [`QCInput.get_string()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.get_string)


            * [`QCInput.molecule_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.molecule_template)


            * [`QCInput.multi_job_string()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.multi_job_string)


            * [`QCInput.nbo_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.nbo_template)


            * [`QCInput.opt_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.opt_template)


            * [`QCInput.pcm_nonels_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.pcm_nonels_template)


            * [`QCInput.pcm_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.pcm_template)


            * [`QCInput.plots_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.plots_template)


            * [`QCInput.read_almo()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_almo)


            * [`QCInput.read_cdft()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_cdft)


            * [`QCInput.read_geom_opt()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_geom_opt)


            * [`QCInput.read_molecule()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_molecule)


            * [`QCInput.read_nbo()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_nbo)


            * [`QCInput.read_opt()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_opt)


            * [`QCInput.read_pcm()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_pcm)


            * [`QCInput.read_pcm_nonels()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_pcm_nonels)


            * [`QCInput.read_plots()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_plots)


            * [`QCInput.read_rem()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_rem)


            * [`QCInput.read_scan()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_scan)


            * [`QCInput.read_smx()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_smx)


            * [`QCInput.read_solvent()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_solvent)


            * [`QCInput.read_svp()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_svp)


            * [`QCInput.read_vdw()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.read_vdw)


            * [`QCInput.rem_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.rem_template)


            * [`QCInput.scan_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.scan_template)


            * [`QCInput.smx_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.smx_template)


            * [`QCInput.solvent_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.solvent_template)


            * [`QCInput.svp_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.svp_template)


            * [`QCInput.van_der_waals_template()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.van_der_waals_template)


            * [`QCInput.write_multi_job_file()`](pymatgen.io.qchem.md#pymatgen.io.qchem.inputs.QCInput.write_multi_job_file)


    * [pymatgen.io.qchem.outputs module](pymatgen.io.qchem.md#module-pymatgen.io.qchem.outputs)


        * [`QCOutput`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput)


            * [`QCOutput._check_completion_errors()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._check_completion_errors)


            * [`QCOutput._detect_general_warnings()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._detect_general_warnings)


            * [`QCOutput._get_grad_format_length()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._get_grad_format_length)


            * [`QCOutput._read_SCF()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_SCF)


            * [`QCOutput._read_almo_msdft()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_almo_msdft)


            * [`QCOutput._read_cdft()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_cdft)


            * [`QCOutput._read_charge_and_multiplicity()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_charge_and_multiplicity)


            * [`QCOutput._read_charges_and_dipoles()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_charges_and_dipoles)


            * [`QCOutput._read_cmirs_information()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_cmirs_information)


            * [`QCOutput._read_coefficient_matrix()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_coefficient_matrix)


            * [`QCOutput._read_eigenvalues()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_eigenvalues)


            * [`QCOutput._read_fock_matrix()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_fock_matrix)


            * [`QCOutput._read_force_data()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_force_data)


            * [`QCOutput._read_frequency_data()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_frequency_data)


            * [`QCOutput._read_geometries()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_geometries)


            * [`QCOutput._read_gradients()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_gradients)


            * [`QCOutput._read_isosvp_information()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_isosvp_information)


            * [`QCOutput._read_nbo_data()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_nbo_data)


            * [`QCOutput._read_optimization_data()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_optimization_data)


            * [`QCOutput._read_pcm_information()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_pcm_information)


            * [`QCOutput._read_scan_data()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_scan_data)


            * [`QCOutput._read_smd_information()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_smd_information)


            * [`QCOutput._read_species_and_inital_geometry()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput._read_species_and_inital_geometry)


            * [`QCOutput.as_dict()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput.as_dict)


            * [`QCOutput.multiple_outputs_from_file()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput.multiple_outputs_from_file)


        * [`check_for_structure_changes()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.check_for_structure_changes)


        * [`get_percentage()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.get_percentage)


        * [`jump_to_header()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.jump_to_header)


        * [`nbo_parser()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.nbo_parser)


        * [`parse_hybridization_character()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.parse_hybridization_character)


        * [`parse_hyperbonds()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.parse_hyperbonds)


        * [`parse_natural_populations()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.parse_natural_populations)


        * [`parse_perturbation_energy()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.parse_perturbation_energy)


        * [`z_int()`](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.z_int)


    * [pymatgen.io.qchem.sets module](pymatgen.io.qchem.md#module-pymatgen.io.qchem.sets)


        * [`ForceSet`](pymatgen.io.qchem.md#pymatgen.io.qchem.sets.ForceSet)


        * [`FreqSet`](pymatgen.io.qchem.md#pymatgen.io.qchem.sets.FreqSet)


        * [`OptSet`](pymatgen.io.qchem.md#pymatgen.io.qchem.sets.OptSet)


        * [`PESScanSet`](pymatgen.io.qchem.md#pymatgen.io.qchem.sets.PESScanSet)


        * [`QChemDictSet`](pymatgen.io.qchem.md#pymatgen.io.qchem.sets.QChemDictSet)


            * [`QChemDictSet.write()`](pymatgen.io.qchem.md#pymatgen.io.qchem.sets.QChemDictSet.write)


        * [`SinglePointSet`](pymatgen.io.qchem.md#pymatgen.io.qchem.sets.SinglePointSet)


        * [`TransitionStateSet`](pymatgen.io.qchem.md#pymatgen.io.qchem.sets.TransitionStateSet)


    * [pymatgen.io.qchem.utils module](pymatgen.io.qchem.md#module-pymatgen.io.qchem.utils)


        * [`lower_and_check_unique()`](pymatgen.io.qchem.md#pymatgen.io.qchem.utils.lower_and_check_unique)


        * [`process_parsed_HESS()`](pymatgen.io.qchem.md#pymatgen.io.qchem.utils.process_parsed_HESS)


        * [`process_parsed_coords()`](pymatgen.io.qchem.md#pymatgen.io.qchem.utils.process_parsed_coords)


        * [`process_parsed_fock_matrix()`](pymatgen.io.qchem.md#pymatgen.io.qchem.utils.process_parsed_fock_matrix)


        * [`read_matrix_pattern()`](pymatgen.io.qchem.md#pymatgen.io.qchem.utils.read_matrix_pattern)


        * [`read_pattern()`](pymatgen.io.qchem.md#pymatgen.io.qchem.utils.read_pattern)


        * [`read_table_pattern()`](pymatgen.io.qchem.md#pymatgen.io.qchem.utils.read_table_pattern)


* [pymatgen.io.vasp package](pymatgen.io.vasp.md)




    * [pymatgen.io.vasp.help module](pymatgen.io.vasp.md#module-pymatgen.io.vasp.help)


        * [`VaspDoc`](pymatgen.io.vasp.md#pymatgen.io.vasp.help.VaspDoc)


            * [`VaspDoc.get_help()`](pymatgen.io.vasp.md#pymatgen.io.vasp.help.VaspDoc.get_help)


            * [`VaspDoc.get_incar_tags()`](pymatgen.io.vasp.md#pymatgen.io.vasp.help.VaspDoc.get_incar_tags)


            * [`VaspDoc.print_help()`](pymatgen.io.vasp.md#pymatgen.io.vasp.help.VaspDoc.print_help)


            * [`VaspDoc.print_jupyter_help()`](pymatgen.io.vasp.md#pymatgen.io.vasp.help.VaspDoc.print_jupyter_help)


    * [pymatgen.io.vasp.inputs module](pymatgen.io.vasp.md#module-pymatgen.io.vasp.inputs)


        * [`BadIncarWarning`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.BadIncarWarning)


        * [`Incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar)


            * [`Incar.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.as_dict)


            * [`Incar.check_params()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.check_params)


            * [`Incar.diff()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.diff)


            * [`Incar.from_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.from_dict)


            * [`Incar.from_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.from_file)


            * [`Incar.from_str()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.from_str)


            * [`Incar.from_string()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.from_string)


            * [`Incar.get_str()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.get_str)


            * [`Incar.get_string()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.get_string)


            * [`Incar.proc_val()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.proc_val)


            * [`Incar.write_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Incar.write_file)


        * [`Kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints)


            * [`Kpoints.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.as_dict)


            * [`Kpoints.automatic()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.automatic)


            * [`Kpoints.automatic_density()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.automatic_density)


            * [`Kpoints.automatic_density_by_lengths()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_lengths)


            * [`Kpoints.automatic_density_by_vol()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_vol)


            * [`Kpoints.automatic_gamma_density()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.automatic_gamma_density)


            * [`Kpoints.automatic_linemode()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.automatic_linemode)


            * [`Kpoints.from_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.from_dict)


            * [`Kpoints.from_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.from_file)


            * [`Kpoints.from_str()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.from_str)


            * [`Kpoints.from_string()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.from_string)


            * [`Kpoints.gamma_automatic()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.gamma_automatic)


            * [`Kpoints.monkhorst_automatic()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.monkhorst_automatic)


            * [`Kpoints.style`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.style)


            * [`Kpoints.supported_modes`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.supported_modes)


            * [`Kpoints.write_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints.write_file)


        * [`KpointsSupportedModes`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.KpointsSupportedModes)


            * [`KpointsSupportedModes.Automatic`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Automatic)


            * [`KpointsSupportedModes.Cartesian`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Cartesian)


            * [`KpointsSupportedModes.Gamma`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Gamma)


            * [`KpointsSupportedModes.Line_mode`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Line_mode)


            * [`KpointsSupportedModes.Monkhorst`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Monkhorst)


            * [`KpointsSupportedModes.Reciprocal`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Reciprocal)


            * [`KpointsSupportedModes.from_str()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.from_str)


            * [`KpointsSupportedModes.from_string()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.from_string)


        * [`Orbital`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital)


            * [`Orbital.E`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital.E)


            * [`Orbital._asdict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital._asdict)


            * [`Orbital._field_defaults`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital._field_defaults)


            * [`Orbital._fields`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital._fields)


            * [`Orbital._make()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital._make)


            * [`Orbital._replace()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital._replace)


            * [`Orbital.j`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital.j)


            * [`Orbital.l`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital.l)


            * [`Orbital.n`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital.n)


            * [`Orbital.occ`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Orbital.occ)


        * [`OrbitalDescription`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription)


            * [`OrbitalDescription.E`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription.E)


            * [`OrbitalDescription.Rcut`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription.Rcut)


            * [`OrbitalDescription.Rcut2`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription.Rcut2)


            * [`OrbitalDescription.Type`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription.Type)


            * [`OrbitalDescription.Type2`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription.Type2)


            * [`OrbitalDescription._asdict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription._asdict)


            * [`OrbitalDescription._field_defaults`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription._field_defaults)


            * [`OrbitalDescription._fields`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription._fields)


            * [`OrbitalDescription._make()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription._make)


            * [`OrbitalDescription._replace()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription._replace)


            * [`OrbitalDescription.l`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.OrbitalDescription.l)


        * [`Poscar`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar)


            * [`Poscar.structure`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.structure)


            * [`Poscar.comment`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.comment)


            * [`Poscar.true_names`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.true_names)


            * [`Poscar.selective_dynamics`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.selective_dynamics)


            * [`Poscar.velocities`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.velocities)


            * [`Poscar.predictor_corrector`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.predictor_corrector)


            * [`Poscar.predictor_corrector_preamble`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.predictor_corrector_preamble)


            * [`Poscar.temperature`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.temperature)


            * [`Poscar.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.as_dict)


            * [`Poscar.from_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.from_dict)


            * [`Poscar.from_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.from_file)


            * [`Poscar.from_str()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.from_str)


            * [`Poscar.from_string()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.from_string)


            * [`Poscar.get_str()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.get_str)


            * [`Poscar.get_string()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.get_string)


            * [`Poscar.natoms`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.natoms)


            * [`Poscar.predictor_corrector`](pymatgen.io.vasp.md#id0)


            * [`Poscar.predictor_corrector_preamble`](pymatgen.io.vasp.md#id1)


            * [`Poscar.selective_dynamics`](pymatgen.io.vasp.md#id2)


            * [`Poscar.set_temperature()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.set_temperature)


            * [`Poscar.site_symbols`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.site_symbols)


            * [`Poscar.velocities`](pymatgen.io.vasp.md#id3)


            * [`Poscar.write_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Poscar.write_file)


        * [`Potcar`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Potcar)


            * [`Potcar.FUNCTIONAL_CHOICES`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Potcar.FUNCTIONAL_CHOICES)


            * [`Potcar.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Potcar.as_dict)


            * [`Potcar.from_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Potcar.from_dict)


            * [`Potcar.from_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Potcar.from_file)


            * [`Potcar.set_symbols()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Potcar.set_symbols)


            * [`Potcar.spec`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Potcar.spec)


            * [`Potcar.symbols`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Potcar.symbols)


            * [`Potcar.write_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Potcar.write_file)


        * [`PotcarSingle`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle)


            * [`PotcarSingle.data`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.data)


            * [`PotcarSingle.keywords`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.keywords)


            * [`PotcarSingle.atomic_no`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.atomic_no)


            * [`PotcarSingle.electron_configuration`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.electron_configuration)


            * [`PotcarSingle.element`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.element)


            * [`PotcarSingle.from_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.from_file)


            * [`PotcarSingle.from_symbol_and_functional()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.from_symbol_and_functional)


            * [`PotcarSingle.functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.functional)


            * [`PotcarSingle.functional_class`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.functional_class)


            * [`PotcarSingle.functional_dir`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.functional_dir)


            * [`PotcarSingle.functional_tags`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.functional_tags)


            * [`PotcarSingle.get_potcar_file_hash()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.get_potcar_file_hash)


            * [`PotcarSingle.get_potcar_hash()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.get_potcar_hash)


            * [`PotcarSingle.get_sha256_file_hash()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.get_sha256_file_hash)


            * [`PotcarSingle.identify_potcar()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.identify_potcar)


            * [`PotcarSingle.nelectrons`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.nelectrons)


            * [`PotcarSingle.parse_functions`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.parse_functions)


            * [`PotcarSingle.potential_type`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.potential_type)


            * [`PotcarSingle.symbol`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.symbol)


            * [`PotcarSingle.verify_potcar()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.verify_potcar)


            * [`PotcarSingle.write_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.PotcarSingle.write_file)


        * [`UnknownPotcarWarning`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.UnknownPotcarWarning)


        * [`VaspInput`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.VaspInput)


            * [`VaspInput.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.VaspInput.as_dict)


            * [`VaspInput.from_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.VaspInput.from_dict)


            * [`VaspInput.from_directory()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.VaspInput.from_directory)


            * [`VaspInput.run_vasp()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.VaspInput.run_vasp)


            * [`VaspInput.write_input()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.VaspInput.write_input)


        * [`_parse_bool()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs._parse_bool)


        * [`_parse_float()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs._parse_float)


        * [`_parse_int()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs._parse_int)


        * [`_parse_list()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs._parse_list)


        * [`_parse_string()`](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs._parse_string)


    * [pymatgen.io.vasp.optics module](pymatgen.io.vasp.md#module-pymatgen.io.vasp.optics)


        * [`DielectricFunctionCalculator`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator)


            * [`DielectricFunctionCalculator.cder`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.cder)


            * [`DielectricFunctionCalculator.cder_imag`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.cder_imag)


            * [`DielectricFunctionCalculator.cder_real`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.cder_real)


            * [`DielectricFunctionCalculator.cshift`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.cshift)


            * [`DielectricFunctionCalculator.deltae`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.deltae)


            * [`DielectricFunctionCalculator.efermi`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.efermi)


            * [`DielectricFunctionCalculator.eigs`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.eigs)


            * [`DielectricFunctionCalculator.from_directory()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.from_directory)


            * [`DielectricFunctionCalculator.from_vasp_objects()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.from_vasp_objects)


            * [`DielectricFunctionCalculator.get_epsilon()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.get_epsilon)


            * [`DielectricFunctionCalculator.ismear`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.ismear)


            * [`DielectricFunctionCalculator.ispin`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.ispin)


            * [`DielectricFunctionCalculator.kweights`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.kweights)


            * [`DielectricFunctionCalculator.nedos`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.nedos)


            * [`DielectricFunctionCalculator.plot_weighted_transition_data()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.plot_weighted_transition_data)


            * [`DielectricFunctionCalculator.sigma`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.sigma)


            * [`DielectricFunctionCalculator.volume`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.volume)


        * [`delta_func()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.delta_func)


        * [`delta_methfessel_paxton()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.delta_methfessel_paxton)


        * [`epsilon_imag()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.epsilon_imag)


        * [`get_delta()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.get_delta)


        * [`get_step()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.get_step)


        * [`kramers_kronig()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.kramers_kronig)


        * [`step_func()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.step_func)


        * [`step_methfessel_paxton()`](pymatgen.io.vasp.md#pymatgen.io.vasp.optics.step_methfessel_paxton)


    * [pymatgen.io.vasp.outputs module](pymatgen.io.vasp.md#module-pymatgen.io.vasp.outputs)


        * [`BSVasprun`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.BSVasprun)


            * [`BSVasprun.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.BSVasprun.as_dict)


        * [`Chgcar`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Chgcar)


            * [`Chgcar.from_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Chgcar.from_file)


            * [`Chgcar.net_magnetization`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Chgcar.net_magnetization)


        * [`Dynmat`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Dynmat)


            * [`Dynmat.data`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Dynmat.data)


            * [`Dynmat.get_phonon_frequencies()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Dynmat.get_phonon_frequencies)


            * [`Dynmat.masses`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Dynmat.masses)


            * [`Dynmat.natoms`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Dynmat.natoms)


            * [`Dynmat.ndisps`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Dynmat.ndisps)


            * [`Dynmat.nspecs`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Dynmat.nspecs)


        * [`Eigenval`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval)


            * [`Eigenval.filename`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.filename)


            * [`Eigenval.occu_tol`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.occu_tol)


            * [`Eigenval.ispin`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.ispin)


            * [`Eigenval.nelect`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.nelect)


            * [`Eigenval.nkpt`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.nkpt)


            * [`Eigenval.nbands`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.nbands)


            * [`Eigenval.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.kpoints)


            * [`Eigenval.kpoints_weights`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.kpoints_weights)


            * [`Eigenval.eigenvalues`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.eigenvalues)


            * [`Eigenval.eigenvalue_band_properties`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Eigenval.eigenvalue_band_properties)


        * [`Elfcar`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Elfcar)


            * [`Elfcar.from_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Elfcar.from_file)


            * [`Elfcar.get_alpha()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Elfcar.get_alpha)


        * [`Locpot`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Locpot)


            * [`Locpot.from_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Locpot.from_file)


        * [`Oszicar`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Oszicar)


            * [`Oszicar.electronic_steps`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Oszicar.electronic_steps)


            * [`Oszicar.ionic_steps`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Oszicar.ionic_steps)


            * [`Oszicar.all_energies`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Oszicar.all_energies)


            * [`Oszicar.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Oszicar.as_dict)


            * [`Oszicar.final_energy`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Oszicar.final_energy)


        * [`Outcar`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar)


            * [`Outcar.magnetization`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.magnetization)


            * [`Outcar.chemical_shielding`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.chemical_shielding)


            * [`Outcar.unsym_cs_tensor`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.unsym_cs_tensor)


            * [`Outcar.cs_g0_contribution`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.cs_g0_contribution)


            * [`Outcar.cs_core_contribution`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.cs_core_contribution)


            * [`Outcar.efg`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.efg)


            * [`Outcar.charge`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.charge)


            * [`Outcar.is_stopped`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.is_stopped)


            * [`Outcar.run_stats`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.run_stats)


            * [`Outcar.elastic_tensor`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.elastic_tensor)


            * [`Outcar.drift`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.drift)


            * [`Outcar.ngf`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.ngf)


            * [`Outcar.sampling_radii`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.sampling_radii)


            * [`Outcar.electrostatic_potential`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.electrostatic_potential)


            * [`Outcar.final_energy_contribs`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.final_energy_contribs)


            * [`Outcar.efermi`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.efermi)


            * [`Outcar.filename`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.filename)


            * [`Outcar.final_energy`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.final_energy)


            * [`Outcar.final_energy_wo_entrp`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.final_energy_wo_entrp)


            * [`Outcar.final_fr_energy`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.final_fr_energy)


            * [`Outcar.has_onsite_density_matrices`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.has_onsite_density_matrices)


            * [`Outcar.lcalcpol`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.lcalcpol)


            * [`Outcar.lepsilon`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.lepsilon)


            * [`Outcar.nelect`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.nelect)


            * [`Outcar.spin`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.spin)


            * [`Outcar.total_mag`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.total_mag)


            * [`Outcar._parse_sci_notation()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar._parse_sci_notation)


            * [`Outcar.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.as_dict)


            * [`Outcar.read_avg_core_poten()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_avg_core_poten)


            * [`Outcar.read_chemical_shielding()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_chemical_shielding)


            * [`Outcar.read_core_state_eigen()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_core_state_eigen)


            * [`Outcar.read_corrections()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_corrections)


            * [`Outcar.read_cs_core_contribution()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_cs_core_contribution)


            * [`Outcar.read_cs_g0_contribution()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_cs_g0_contribution)


            * [`Outcar.read_cs_raw_symmetrized_tensors()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_cs_raw_symmetrized_tensors)


            * [`Outcar.read_elastic_tensor()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_elastic_tensor)


            * [`Outcar.read_electrostatic_potential()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_electrostatic_potential)


            * [`Outcar.read_fermi_contact_shift()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_fermi_contact_shift)


            * [`Outcar.read_freq_dielectric()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_freq_dielectric)


            * [`Outcar.read_igpar()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_igpar)


            * [`Outcar.read_internal_strain_tensor()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_internal_strain_tensor)


            * [`Outcar.read_lcalcpol()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_lcalcpol)


            * [`Outcar.read_lepsilon()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_lepsilon)


            * [`Outcar.read_lepsilon_ionic()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_lepsilon_ionic)


            * [`Outcar.read_neb()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_neb)


            * [`Outcar.read_nmr_efg()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_nmr_efg)


            * [`Outcar.read_nmr_efg_tensor()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_nmr_efg_tensor)


            * [`Outcar.read_onsite_density_matrices()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_onsite_density_matrices)


            * [`Outcar.read_pattern()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_pattern)


            * [`Outcar.read_piezo_tensor()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_piezo_tensor)


            * [`Outcar.read_pseudo_zval()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_pseudo_zval)


            * [`Outcar.read_table_pattern()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar.read_table_pattern)


        * [`Procar`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Procar)


            * [`Procar.data`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Procar.data)


            * [`Procar.weights`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Procar.weights)


            * [`Procar.phase_factors`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Procar.phase_factors)


            * [`Procar.nbands`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Procar.nbands)


            * [`Procar.nkpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Procar.nkpoints)


            * [`Procar.nions`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Procar.nions)


            * [`Procar.get_occupation()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Procar.get_occupation)


            * [`Procar.get_projection_on_elements()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Procar.get_projection_on_elements)


        * [`UnconvergedVASPWarning`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.UnconvergedVASPWarning)


        * [`VaspParseError`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.VaspParseError)


        * [`Vasprun`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun)


            * [`Vasprun.ionic_steps`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.ionic_steps)


            * [`Vasprun.tdos`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.tdos)


            * [`Vasprun.idos`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.idos)


            * [`Vasprun.pdos`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.pdos)


            * [`Vasprun.efermi`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.efermi)


            * [`Vasprun.eigenvalues`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.eigenvalues)


            * [`Vasprun.projected_eigenvalues`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.projected_eigenvalues)


            * [`Vasprun.projected_magnetisation`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.projected_magnetisation)


            * [`Vasprun.other_dielectric`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.other_dielectric)


            * [`Vasprun.nionic_steps`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.nionic_steps)


            * [`Vasprun.force_constants`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.force_constants)


            * [`Vasprun.normalmode_eigenvals`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.normalmode_eigenvals)


            * [`Vasprun.normalmode_eigenvecs`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.normalmode_eigenvecs)


            * [`Vasprun.md_data`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.md_data)


            * [`Vasprun.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.incar)


            * [`Vasprun.parameters`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.parameters)


            * [`Vasprun.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.kpoints)


            * [`Vasprun.actual_kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.actual_kpoints)


            * [`Vasprun.actual_kpoints_weights`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.actual_kpoints_weights)


            * [`Vasprun.atomic_symbols`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.atomic_symbols)


            * [`Vasprun.potcar_symbols`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.potcar_symbols)


            * [`Vasprun._parse()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse)


            * [`Vasprun._parse_atominfo()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_atominfo)


            * [`Vasprun._parse_calculation()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_calculation)


            * [`Vasprun._parse_chemical_shielding_calculation()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_chemical_shielding_calculation)


            * [`Vasprun._parse_diel()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_diel)


            * [`Vasprun._parse_dos()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_dos)


            * [`Vasprun._parse_dynmat()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_dynmat)


            * [`Vasprun._parse_eigen()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_eigen)


            * [`Vasprun._parse_kpoints()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_kpoints)


            * [`Vasprun._parse_optical_transition()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_optical_transition)


            * [`Vasprun._parse_params()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_params)


            * [`Vasprun._parse_projected_eigen()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_projected_eigen)


            * [`Vasprun._parse_structure()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun._parse_structure)


            * [`Vasprun.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.as_dict)


            * [`Vasprun.calculate_efermi()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.calculate_efermi)


            * [`Vasprun.complete_dos`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.complete_dos)


            * [`Vasprun.complete_dos_normalized`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.complete_dos_normalized)


            * [`Vasprun.converged`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.converged)


            * [`Vasprun.converged_electronic`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.converged_electronic)


            * [`Vasprun.converged_ionic`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.converged_ionic)


            * [`Vasprun.dielectric`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.dielectric)


            * [`Vasprun.eigenvalue_band_properties`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.eigenvalue_band_properties)


            * [`Vasprun.epsilon_ionic`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.epsilon_ionic)


            * [`Vasprun.epsilon_static`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.epsilon_static)


            * [`Vasprun.epsilon_static_wolfe`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.epsilon_static_wolfe)


            * [`Vasprun.final_energy`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.final_energy)


            * [`Vasprun.get_band_structure()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.get_band_structure)


            * [`Vasprun.get_computed_entry()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.get_computed_entry)


            * [`Vasprun.get_potcars()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.get_potcars)


            * [`Vasprun.get_trajectory()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.get_trajectory)


            * [`Vasprun.hubbards`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.hubbards)


            * [`Vasprun.is_hubbard`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.is_hubbard)


            * [`Vasprun.is_spin`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.is_spin)


            * [`Vasprun.optical_absorption_coeff`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.optical_absorption_coeff)


            * [`Vasprun.run_type`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.run_type)


            * [`Vasprun.structures`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.structures)


            * [`Vasprun.update_charge_from_potcar()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.update_charge_from_potcar)


            * [`Vasprun.update_potcar_spec()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Vasprun.update_potcar_spec)


        * [`VolumetricData`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.VolumetricData)


            * [`VolumetricData.parse_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.VolumetricData.parse_file)


            * [`VolumetricData.write_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.VolumetricData.write_file)


        * [`WSWQ`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.WSWQ)


            * [`WSWQ.nspin`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.WSWQ.nspin)


            * [`WSWQ.nkpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.WSWQ.nkpoints)


            * [`WSWQ.nbands`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.WSWQ.nbands)


            * [`WSWQ.me_real`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.WSWQ.me_real)


            * [`WSWQ.me_imag`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.WSWQ.me_imag)


            * [`WSWQ.data`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.WSWQ.data)


            * [`WSWQ.from_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.WSWQ.from_file)


            * [`WSWQ.me_imag`](pymatgen.io.vasp.md#id4)


            * [`WSWQ.me_real`](pymatgen.io.vasp.md#id5)


            * [`WSWQ.nbands`](pymatgen.io.vasp.md#id6)


            * [`WSWQ.nkpoints`](pymatgen.io.vasp.md#id7)


            * [`WSWQ.nspin`](pymatgen.io.vasp.md#id8)


        * [`Wavecar`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar)


            * [`Wavecar.filename`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.filename)


            * [`Wavecar.vasp_type`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.vasp_type)


            * [`Wavecar.nk`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.nk)


            * [`Wavecar.nb`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.nb)


            * [`Wavecar.encut`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.encut)


            * [`Wavecar.efermi`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.efermi)


            * [`Wavecar.a`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.a)


            * [`Wavecar.b`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.b)


            * [`Wavecar.vol`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.vol)


            * [`Wavecar.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.kpoints)


            * [`Wavecar.band_energy`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.band_energy)


            * [`Wavecar.Gpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.Gpoints)


            * [`Wavecar.coeffs`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.coeffs)


            * [`Wavecar._generate_G_points()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar._generate_G_points)


            * [`Wavecar._generate_nbmax()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar._generate_nbmax)


            * [`Wavecar.evaluate_wavefunc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.evaluate_wavefunc)


            * [`Wavecar.fft_mesh()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.fft_mesh)


            * [`Wavecar.get_parchg()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.get_parchg)


            * [`Wavecar.write_unks()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Wavecar.write_unks)


        * [`Waveder`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder)


            * [`Waveder.cder_real`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder.cder_real)


            * [`Waveder.cder_imag`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder.cder_imag)


            * [`Waveder.cder`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder.cder)


            * [`Waveder.cder_imag`](pymatgen.io.vasp.md#id9)


            * [`Waveder.cder_real`](pymatgen.io.vasp.md#id10)


            * [`Waveder.from_binary()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder.from_binary)


            * [`Waveder.from_formatted()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder.from_formatted)


            * [`Waveder.get_orbital_derivative_between_states()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder.get_orbital_derivative_between_states)


            * [`Waveder.nbands`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder.nbands)


            * [`Waveder.nkpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder.nkpoints)


            * [`Waveder.nspin`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Waveder.nspin)


        * [`Xdatcar`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Xdatcar)


            * [`Xdatcar.structures`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Xdatcar.structures)


            * [`Xdatcar.comment`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Xdatcar.comment)


            * [`Xdatcar.concatenate()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Xdatcar.concatenate)


            * [`Xdatcar.get_str()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Xdatcar.get_str)


            * [`Xdatcar.get_string()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Xdatcar.get_string)


            * [`Xdatcar.natoms`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Xdatcar.natoms)


            * [`Xdatcar.site_symbols`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Xdatcar.site_symbols)


            * [`Xdatcar.write_file()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Xdatcar.write_file)


        * [`_parse_from_incar()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs._parse_from_incar)


        * [`_parse_parameters()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs._parse_parameters)


        * [`_parse_v_parameters()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs._parse_v_parameters)


        * [`_parse_varray()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs._parse_varray)


        * [`_vasprun_float()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs._vasprun_float)


        * [`get_adjusted_fermi_level()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.get_adjusted_fermi_level)


        * [`get_band_structure_from_vasp_multiple_branches()`](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.get_band_structure_from_vasp_multiple_branches)


    * [pymatgen.io.vasp.sets module](pymatgen.io.vasp.md#module-pymatgen.io.vasp.sets)


        * [`BadInputSetWarning`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.BadInputSetWarning)


        * [`DictSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet)


            * [`DictSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet._abc_impl)


            * [`DictSet.calculate_ng()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet.calculate_ng)


            * [`DictSet.estimate_nbands()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet.estimate_nbands)


            * [`DictSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet.incar)


            * [`DictSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet.kpoints)


            * [`DictSet.nelect`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet.nelect)


            * [`DictSet.poscar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet.poscar)


            * [`DictSet.potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet.potcar_functional)


            * [`DictSet.structure`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet.structure)


            * [`DictSet.write_input()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.DictSet.write_input)


        * [`LobsterSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.LobsterSet)


            * [`LobsterSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.LobsterSet._abc_impl)


            * [`LobsterSet._valid_potcars`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.LobsterSet._valid_potcars)


            * [`LobsterSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.LobsterSet.incar)


            * [`LobsterSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.LobsterSet.user_potcar_functional)


        * [`MITMDSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITMDSet)


            * [`MITMDSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITMDSet._abc_impl)


            * [`MITMDSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITMDSet.kpoints)


            * [`MITMDSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITMDSet.user_potcar_functional)


        * [`MITNEBSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITNEBSet)


            * [`MITNEBSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITNEBSet._abc_impl)


            * [`MITNEBSet._process_structures()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITNEBSet._process_structures)


            * [`MITNEBSet.poscar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITNEBSet.poscar)


            * [`MITNEBSet.poscars`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITNEBSet.poscars)


            * [`MITNEBSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITNEBSet.user_potcar_functional)


            * [`MITNEBSet.write_input()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITNEBSet.write_input)


        * [`MITRelaxSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITRelaxSet)


            * [`MITRelaxSet.CONFIG`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITRelaxSet.CONFIG)


            * [`MITRelaxSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITRelaxSet._abc_impl)


            * [`MITRelaxSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MITRelaxSet.user_potcar_functional)


        * [`MPAbsorptionSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPAbsorptionSet)


            * [`MPAbsorptionSet.SUPPORTED_MODES`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPAbsorptionSet.SUPPORTED_MODES)


            * [`MPAbsorptionSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPAbsorptionSet._abc_impl)


            * [`MPAbsorptionSet.from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPAbsorptionSet.from_prev_calc)


            * [`MPAbsorptionSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPAbsorptionSet.incar)


            * [`MPAbsorptionSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPAbsorptionSet.kpoints)


            * [`MPAbsorptionSet.override_from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPAbsorptionSet.override_from_prev_calc)


            * [`MPAbsorptionSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPAbsorptionSet.user_potcar_functional)


        * [`MPHSEBSSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSEBSSet)


            * [`MPHSEBSSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSEBSSet._abc_impl)


            * [`MPHSEBSSet.from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSEBSSet.from_prev_calc)


            * [`MPHSEBSSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSEBSSet.kpoints)


            * [`MPHSEBSSet.override_from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSEBSSet.override_from_prev_calc)


            * [`MPHSEBSSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSEBSSet.user_potcar_functional)


        * [`MPHSERelaxSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSERelaxSet)


            * [`MPHSERelaxSet.CONFIG`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSERelaxSet.CONFIG)


            * [`MPHSERelaxSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSERelaxSet._abc_impl)


            * [`MPHSERelaxSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPHSERelaxSet.user_potcar_functional)


        * [`MPMDSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPMDSet)


            * [`MPMDSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPMDSet._abc_impl)


            * [`MPMDSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPMDSet.incar)


            * [`MPMDSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPMDSet.kpoints)


            * [`MPMDSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPMDSet.user_potcar_functional)


        * [`MPMetalRelaxSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPMetalRelaxSet)


            * [`MPMetalRelaxSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPMetalRelaxSet._abc_impl)


            * [`MPMetalRelaxSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPMetalRelaxSet.user_potcar_functional)


        * [`MPNMRSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNMRSet)


            * [`MPNMRSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNMRSet._abc_impl)


            * [`MPNMRSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNMRSet.incar)


            * [`MPNMRSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNMRSet.user_potcar_functional)


        * [`MPNonSCFSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNonSCFSet)


            * [`MPNonSCFSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNonSCFSet._abc_impl)


            * [`MPNonSCFSet.from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNonSCFSet.from_prev_calc)


            * [`MPNonSCFSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNonSCFSet.incar)


            * [`MPNonSCFSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNonSCFSet.kpoints)


            * [`MPNonSCFSet.override_from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNonSCFSet.override_from_prev_calc)


            * [`MPNonSCFSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPNonSCFSet.user_potcar_functional)


        * [`MPRelaxSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPRelaxSet)


            * [`MPRelaxSet.CONFIG`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPRelaxSet.CONFIG)


            * [`MPRelaxSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPRelaxSet._abc_impl)


            * [`MPRelaxSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPRelaxSet.user_potcar_functional)


        * [`MPSOCSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPSOCSet)


            * [`MPSOCSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPSOCSet._abc_impl)


            * [`MPSOCSet.from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPSOCSet.from_prev_calc)


            * [`MPSOCSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPSOCSet.incar)


            * [`MPSOCSet.override_from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPSOCSet.override_from_prev_calc)


            * [`MPSOCSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPSOCSet.user_potcar_functional)


        * [`MPScanRelaxSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanRelaxSet)


            * [`MPScanRelaxSet.CONFIG`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanRelaxSet.CONFIG)


            * [`MPScanRelaxSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanRelaxSet._abc_impl)


            * [`MPScanRelaxSet._valid_potcars`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanRelaxSet._valid_potcars)


            * [`MPScanRelaxSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanRelaxSet.user_potcar_functional)


        * [`MPScanStaticSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanStaticSet)


            * [`MPScanStaticSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanStaticSet._abc_impl)


            * [`MPScanStaticSet.from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanStaticSet.from_prev_calc)


            * [`MPScanStaticSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanStaticSet.incar)


            * [`MPScanStaticSet.override_from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanStaticSet.override_from_prev_calc)


            * [`MPScanStaticSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPScanStaticSet.user_potcar_functional)


        * [`MPStaticSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPStaticSet)


            * [`MPStaticSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPStaticSet._abc_impl)


            * [`MPStaticSet.from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPStaticSet.from_prev_calc)


            * [`MPStaticSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPStaticSet.incar)


            * [`MPStaticSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPStaticSet.kpoints)


            * [`MPStaticSet.override_from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPStaticSet.override_from_prev_calc)


            * [`MPStaticSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MPStaticSet.user_potcar_functional)


        * [`MVLElasticSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLElasticSet)


            * [`MVLElasticSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLElasticSet._abc_impl)


            * [`MVLElasticSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLElasticSet.user_potcar_functional)


        * [`MVLGBSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGBSet)


            * [`MVLGBSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGBSet._abc_impl)


            * [`MVLGBSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGBSet.incar)


            * [`MVLGBSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGBSet.kpoints)


            * [`MVLGBSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGBSet.user_potcar_functional)


        * [`MVLGWSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGWSet)


            * [`MVLGWSet.CONFIG`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGWSet.CONFIG)


            * [`MVLGWSet.SUPPORTED_MODES`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGWSet.SUPPORTED_MODES)


            * [`MVLGWSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGWSet._abc_impl)


            * [`MVLGWSet.from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGWSet.from_prev_calc)


            * [`MVLGWSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGWSet.incar)


            * [`MVLGWSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGWSet.kpoints)


            * [`MVLGWSet.override_from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGWSet.override_from_prev_calc)


            * [`MVLGWSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLGWSet.user_potcar_functional)


        * [`MVLNPTMDSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLNPTMDSet)


            * [`MVLNPTMDSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLNPTMDSet._abc_impl)


            * [`MVLNPTMDSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLNPTMDSet.incar)


            * [`MVLNPTMDSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLNPTMDSet.user_potcar_functional)


        * [`MVLRelax52Set`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLRelax52Set)


            * [`MVLRelax52Set.CONFIG`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLRelax52Set.CONFIG)


            * [`MVLRelax52Set._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLRelax52Set._abc_impl)


            * [`MVLRelax52Set._valid_potcars`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLRelax52Set._valid_potcars)


            * [`MVLRelax52Set.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLRelax52Set.user_potcar_functional)


        * [`MVLScanRelaxSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLScanRelaxSet)


            * [`MVLScanRelaxSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLScanRelaxSet._abc_impl)


            * [`MVLScanRelaxSet._valid_potcars`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLScanRelaxSet._valid_potcars)


            * [`MVLScanRelaxSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLScanRelaxSet.user_potcar_functional)


        * [`MVLSlabSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLSlabSet)


            * [`MVLSlabSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLSlabSet._abc_impl)


            * [`MVLSlabSet.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLSlabSet.as_dict)


            * [`MVLSlabSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLSlabSet.incar)


            * [`MVLSlabSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLSlabSet.kpoints)


            * [`MVLSlabSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MVLSlabSet.user_potcar_functional)


        * [`MatPESStaticSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MatPESStaticSet)


            * [`MatPESStaticSet.CONFIG`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MatPESStaticSet.CONFIG)


            * [`MatPESStaticSet.INHERITED_INCAR_PARAMS`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MatPESStaticSet.INHERITED_INCAR_PARAMS)


            * [`MatPESStaticSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MatPESStaticSet._abc_impl)


            * [`MatPESStaticSet.from_prev_calc()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MatPESStaticSet.from_prev_calc)


            * [`MatPESStaticSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MatPESStaticSet.incar)


            * [`MatPESStaticSet.user_potcar_functional`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.MatPESStaticSet.user_potcar_functional)


        * [`VaspInputSet`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet)


            * [`VaspInputSet._abc_impl`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet._abc_impl)


            * [`VaspInputSet._valid_potcars`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet._valid_potcars)


            * [`VaspInputSet.as_dict()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet.as_dict)


            * [`VaspInputSet.get_vasp_input()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet.get_vasp_input)


            * [`VaspInputSet.incar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet.incar)


            * [`VaspInputSet.kpoints`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet.kpoints)


            * [`VaspInputSet.poscar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet.poscar)


            * [`VaspInputSet.potcar`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet.potcar)


            * [`VaspInputSet.potcar_symbols`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet.potcar_symbols)


            * [`VaspInputSet.write_input()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet.write_input)


        * [`_load_yaml_config()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets._load_yaml_config)


        * [`batch_write_input()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.batch_write_input)


        * [`get_structure_from_prev_run()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.get_structure_from_prev_run)


        * [`get_valid_magmom_struct()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.get_valid_magmom_struct)


        * [`get_vasprun_outcar()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.get_vasprun_outcar)


        * [`next_num_with_prime_factors()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.next_num_with_prime_factors)


        * [`primes_less_than()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.primes_less_than)


        * [`standardize_structure()`](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.standardize_structure)


* [pymatgen.io.xtb package](pymatgen.io.xtb.md)




    * [pymatgen.io.xtb.inputs module](pymatgen.io.xtb.md#module-pymatgen.io.xtb.inputs)


        * [`CRESTInput`](pymatgen.io.xtb.md#pymatgen.io.xtb.inputs.CRESTInput)


            * [`CRESTInput.constrains_template()`](pymatgen.io.xtb.md#pymatgen.io.xtb.inputs.CRESTInput.constrains_template)


            * [`CRESTInput.write_input_files()`](pymatgen.io.xtb.md#pymatgen.io.xtb.inputs.CRESTInput.write_input_files)


    * [pymatgen.io.xtb.outputs module](pymatgen.io.xtb.md#module-pymatgen.io.xtb.outputs)


        * [`CRESTOutput`](pymatgen.io.xtb.md#pymatgen.io.xtb.outputs.CRESTOutput)


            * [`CRESTOutput._parse_crest_output()`](pymatgen.io.xtb.md#pymatgen.io.xtb.outputs.CRESTOutput._parse_crest_output)



## pymatgen.io.adf module

IO for ADF files.


### _class_ AdfInput(task)
Bases: `object`

A basic ADF input file writer.

Initialization method.


* **Parameters**

    **task** (*AdfTask*) – An ADF task.



#### write_file(molecule, inpfile)
Write an ADF input file.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – The molecule for this task.


    * **inpfile** (*str*) – The name where the input file will be saved.



### _exception_ AdfInputError()
Bases: `Exception`

The default error class for ADF.


### _class_ AdfKey(name, options=None, subkeys=None)
Bases: `MSONable`

The basic input unit for ADF. A key is a string of characters that does not
contain a delimiter (blank, comma or equal sign). A key may have multiple
subkeys and a set of options.

Initialization method.


* **Parameters**


    * **name** (*str*) – The name of this key.


    * **options** (*Sized*) – The options for this key. Each element can be a primitive object or
    a tuple/list with two elements: the first is the name and the second
    is a primitive object.


    * **subkeys** (*Sized*) – The subkeys for this key.


    * **Raises** –


    * **------** –


    * **ValueError** – If elements in `subkeys` are not `AdfKey` objects.



#### _full_blocks(_ = ('GEOMETRY', 'SCF', 'UNITS', 'BASIS', 'ANALYTICALFREQ'_ )

#### _options_string()
Return the option string.


#### add_option(option)
Add a new option to this key.


* **Parameters**


    * **option** (*Sized** or **str** or **int** or **float*) – A new option to add. This must have the same format with existing
    options.


    * **Raises** –


    * **------** –


    * **TypeError** – If the format of the given `option` is different.



#### add_subkey(subkey)
Add a new subkey to this key.


* **Parameters**


    * **subkey** (*AdfKey*) – A new subkey.


    * **Notes** –


    * **-----** –


    * **block.** (*Duplicate check will not be performed if this is an 'Atoms'*) –



#### as_dict()
A JSON-serializable dict representation of self.


#### block_keys(_ = ('SCF', 'GEOMETRY', 'XC', 'UNITS', 'ATOMS', 'CHARGE', 'BASIS', 'SYMMETRY', 'RELATIVISTIC', 'OCCUPATIONS', 'SAVE', 'A1FIT', 'INTEGRATION', 'UNRESTRICTED', 'ZLMFIT', 'TITLE', 'EXACTDENSITY', 'TOTALENERGY', 'ANALYTICALFREQ'_ )

#### _classmethod_ from_dict(d)
Construct a MSONable AdfKey object from the JSON dict.


* **Parameters**


    * **d** (*dict*) – A dict of saved attributes.


    * **Returns** –


    * **-------** –


    * **adfkey** (*AdfKey*) – An AdfKey object recovered from the JSON dict `d`.



#### _static_ from_str(string)
Construct an AdfKey object from the string.


* **Parameters**


    * **string** (*str*) – A string.


    * **Returns** –


    * **-------** –


    * **adfkey** (*AdfKey*) – An AdfKey object recovered from the string.


    * **Raises** –


    * **------** –


    * **ValueError** – Currently nested subkeys are not supported. If `subend` was found
    a ValueError would be raised.


    * **Notes** –


    * **-----** –


    * **returned.** (*Only the first block key will be*) –



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### has_option(option)
Return True if the option is included in this key.


* **Parameters**


    * **option** (*str*) – The option.


    * **Returns** –


    * **-------** –


    * **has** (*bool*) – True if the option can be found. Otherwise False will be returned.



#### has_subkey(subkey)
Return True if this AdfKey contains the given subkey.


* **Parameters**


    * **subkey** (*str** or **AdfKey*) – A key name or an AdfKey object.


    * **Returns** –


    * **-------** –


    * **has** (*bool*) – True if this key contains the given key. Otherwise False.



#### is_block_key()
Return True if this key is a block key.


#### _property_ key()
Return the name of this key. If this is a block key, the name will be
converted to upper cases.


#### remove_option(option)
Remove an option.


* **Parameters**


    * **option** (*str** or **int*) – The name (str) or index (int) of the option to remove.


    * **Raises** –


    * **------** –


    * **TypeError** – If the option has a wrong type.



#### remove_subkey(subkey)
Remove the given subkey, if existed, from this AdfKey.


* **Parameters**

    **subkey** (*str** or **AdfKey*) – The subkey to remove.



#### sub_keys(_ = ('AtomDepQuality',_ )

### _class_ AdfOutput(filename)
Bases: `object`

A basic ADF output file parser.

### Attributes:

is_failed

    True is the ADF job is terminated without success. Otherwise False.

is_internal_crash

    True if the job is terminated with internal crash. Please read ‘TAPE13’
    of the ADF manual for more detail.

error

    The error description.

run_type

    The RunType of this ADF job. Possible options are: ‘SinglePoint’,
    ‘GeometryOptimization’, ‘AnalyticalFreq’ and ‘NUmericalFreq’.

final_energy

    The final molecule energy (a.u).

final_structure

    The final structure of the molecule.

energies

    The energy of each cycle.

structures

    The structure of each cycle If geometry optimization is performed.

frequencies

    The frequencies of the molecule.

normal_modes

    The normal modes of the molecule.

freq_type

    Either ‘Analytical’ or ‘Numerical’.

Initialization method.


* **param filename**

    The ADF output file to parse.



* **type filename**

    str



#### _parse()
Parse the ADF outputs. There are two files: one is ‘logfile’, the other
is the ADF output file. The final energy and structures are parsed from
the ‘logfile’. Frequencies and normal modes are parsed from the ADF
output file.


#### _parse_adf_output()
Parse the standard ADF output file.


#### _parse_logfile(logfile)
Parse the formatted logfile.


#### _static_ _sites_to_mol(sites)
Return a `Molecule` object given a list of sites.


* **Parameters**


    * **sites** (*list*) – A list of sites.


    * **Returns** –


    * **-------** –


    * **mol** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – A `Molecule` object.



### _exception_ AdfOutputError()
Bases: `Exception`

The default error class for errors raised by `AdfOutput`.


### _class_ AdfTask(operation='energy', basis_set=None, xc=None, title='ADF_RUN', units=None, geo_subkeys=None, scf=None, other_directives=None)
Bases: `MSONable`

Basic task for ADF. All settings in this class are independent of molecules.

### Notes:

Unlike other quantum chemistry packages (NWChem, Gaussian, …), ADF does
not support calculating force/gradient.

Initialization method.


* **param operation**

    The target operation.



* **type operation**

    str



* **param basis_set**

    The basis set definitions for this task. Defaults to ‘DZ/Large’.



* **type basis_set**

    AdfKey



* **param xc**

    The exchange-correlation functionals. Defaults to PBE.



* **type xc**

    AdfKey



* **param title**

    The title of this ADF task.



* **type title**

    str



* **param units**

    The units. Defaults to Angstroms/Degree.



* **type units**

    AdfKey



* **param geo_subkeys**

    The subkeys for the block key ‘GEOMETRY’.



* **type geo_subkeys**

    Sized



* **param scf**

    The scf options.



* **type scf**

    AdfKey



* **param other_directives**

    User-defined directives.



* **type other_directives**

    Sized



#### _setup_task(geo_subkeys)
Setup the block ‘Geometry’ given subkeys and the task.


* **Parameters**


    * **geo_subkeys** (*Sized*) – User-defined subkeys for the block ‘Geometry’.


    * **Notes** –


    * **-----** –


    * **except** (*Most** of **the run types** of **ADF are specified in the Geometry block*) –


    * **'AnalyticFreq'.** (*the*) –



#### as_dict()
A JSON-serializable dict representation of self.


#### _classmethod_ from_dict(d)
Construct a MSONable AdfTask object from the JSON dict.


* **Parameters**


    * **d** (*dict*) – A dict of saved attributes.


    * **Returns** –


    * **-------** –


    * **task** (*AdfTask*) – An AdfTask object recovered from the JSON dict `d`.



#### _static_ get_default_basis_set()
Returns: Default basis set.


#### _static_ get_default_geo()
Returns: ADFKey using default geometry.


#### _static_ get_default_scf()
Returns: ADF using default SCF.


#### _static_ get_default_units()
Returns: Default units.


#### _static_ get_default_xc()
Returns: ADFKey using default XC.


#### operations(_ = {'energy': 'Evaluate the single point energy.', 'freq': 'Same as frequencies.', 'frequencies': 'Compute second derivatives and print out an analysis of molecular vibrations.', 'numerical_frequencies': 'Compute molecular frequencies using numerical method.', 'optimize': 'Minimize the energy by varying the molecular structure.'_ )

### is_numeric(s)
Return True is the string `s` is a numeric string.


* **Parameters**


    * **s** (*str*) – A string.


    * **Returns** –


    * **-------** –


    * **res** (*bool*) – If True, `s` is a numeric string and can be converted to an int or a
    float. Otherwise False will be returned.



### iterlines(s: str)
A generator form of s.split(’n’) for reducing memory overhead.


* **Parameters**

    **s** (*str*) – A multi-line string.



* **Yields**

    *str* – line


## pymatgen.io.ase module

This module provides conversion between the Atomic Simulation Environment
Atoms object and pymatgen Structure objects.


### _class_ AseAtomsAdaptor()
Bases: `object`

Adaptor serves as a bridge between ASE Atoms and pymatgen objects.


#### _static_ get_atoms(structure: [SiteCollection](pymatgen.core.md#pymatgen.core.structure.SiteCollection), \*\*kwargs)
Returns ASE Atoms object from pymatgen structure or molecule.


* **Parameters**


    * **structure** ([*SiteCollection*](pymatgen.core.md#pymatgen.core.structure.SiteCollection)) – pymatgen Structure or Molecule


    * **\*\*kwargs** – passed to the ASE Atoms constructor



* **Returns**

    ASE Atoms object



* **Return type**

    [Atoms](pymatgen.io.feff.md#pymatgen.io.feff.inputs.Atoms)



#### _static_ get_molecule(atoms: ~ase.atoms.Atoms, cls: type[~pymatgen.core.structure.Molecule] = <class 'pymatgen.core.structure.Molecule'>, \*\*cls_kwargs)
Returns pymatgen molecule from ASE Atoms.


* **Parameters**


    * **atoms** – ASE Atoms object


    * **cls** – The Molecule class to instantiate (defaults to pymatgen molecule)


    * **\*\*cls_kwargs** – Any additional kwargs to pass to the cls



* **Returns**

    Equivalent pymatgen.core.structure.Molecule



* **Return type**

    [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule)



#### _static_ get_structure(atoms: ~ase.atoms.Atoms, cls: type[~pymatgen.core.structure.Structure] = <class 'pymatgen.core.structure.Structure'>, \*\*cls_kwargs)
Returns pymatgen structure from ASE Atoms.


* **Parameters**


    * **atoms** – ASE Atoms object


    * **cls** – The Structure class to instantiate (defaults to pymatgen Structure)


    * **\*\*cls_kwargs** – Any additional kwargs to pass to the cls



* **Returns**

    Equivalent pymatgen.core.structure.Structure


## pymatgen.io.atat module

Classes for reading/writing mcsqs files following the rndstr.in format.


### _class_ Mcsqs(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Bases: `object`

Handle input/output for the crystal definition format
used by mcsqs and other ATAT codes.


* **Parameters**

    **Structure** – input Structure.



#### _static_ structure_from_str(data)
Parses a rndstr.in, lat.in or bestsqs.out file into pymatgen’s
Structure format.


* **Parameters**

    **data** – contents of a rndstr.in, lat.in or bestsqs.out file



* **Returns**

    Structure object



#### structure_from_string(\*\*kwargs)

#### to_str()

* **Returns**

    a structure in mcsqs rndstr.in format.



* **Return type**

    str



#### to_string(\*\*kwargs)
## pymatgen.io.babel module

OpenBabel interface module, which opens up access to the hundreds of file
formats supported by OpenBabel. Requires openbabel with python bindings to be
installed. Please consult the openbabel docs [https://openbabel.org](https://openbabel.org).


### _class_ BabelMolAdaptor(mol: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule) | openbabel.OBMol | pybel.Molecule)
Bases: `object`

Adaptor serves as a bridge between OpenBabel’s Molecule and pymatgen’s
Molecule.

Initializes with pymatgen Molecule or OpenBabel’s OBMol.


* **Parameters**

    **mol** – pymatgen’s Molecule/IMolecule or OpenBabel OBMol



#### add_hydrogen()
Add hydrogens (make all hydrogen explicit).


#### confab_conformers(forcefield='mmff94', freeze_atoms=None, rmsd_cutoff=0.5, energy_cutoff=50.0, conf_cutoff=100000, verbose=False)
Conformer generation based on Confab to generate all diverse low-energy
conformers for molecules. This is different from rotor_conformer or
gen3d_conformer as it aims to not simply to find a low energy
conformation but to generate several different conformations.


* **Parameters**


    * **forcefield** (*str*) – Default is mmff94. Options are ‘gaff’, ‘ghemical’,
    ‘mmff94’, ‘mmff94s’, and ‘uff’.


    * **freeze_atoms** (*[**int**]*) – index of atoms to be freezed when performing
    conformer search, default is None.


    * **rmsd_cutoff** (*float*) – rmsd_cufoff, default is 0.5 Angstrom.


    * **energy_cutoff** (*float*) – energy_cutoff, default is 50.0 kcal/mol.


    * **conf_cutoff** (*float*) – max number of conformers to test,
    default is 1 million.


    * **verbose** (*bool*) – whether to display information on torsions found,
    default is False.



* **Returns**

    Molecule objects for generated conformers.



* **Return type**

    list[[Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule)]



#### _static_ from_file(filename, file_format='xyz', return_all_molecules=False)
Uses OpenBabel to read a molecule from a file in all supported formats.


* **Parameters**


    * **filename** – Filename of input file


    * **file_format** – String specifying any OpenBabel supported formats.


    * **return_all_molecules** – If `True`, will return a list of
    `BabelMolAdaptor` instances, one for each molecule found in
    the file. If `False`, will return only the first molecule.



* **Returns**

    BabelMolAdaptor object or list thereof



#### _static_ from_molecule_graph(mol)
Read a molecule from a pymatgen MoleculeGraph object.


* **Parameters**

    **mol** – pymatgen MoleculeGraph object.



* **Returns**

    BabelMolAdaptor object



#### from_str()
staticmethod(function) -> method

Convert a function to be a static method.

A static method does not receive an implicit first argument.
To declare a static method, use this idiom:

> class C:

>     @staticmethod
>     def f(arg1, arg2, …):

>     > …

It can be called either on the class (e.g. C.f()) or on an instance
(e.g. C().f()). Both the class and the instance are ignored, and
neither is passed implicitly as the first argument to the method.

Static methods in Python are similar to those found in Java or C++.
For a more advanced concept, see the classmethod builtin.


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### gen3d_conformer()
A combined method to first generate 3D structures from 0D or 2D
structures and then find the minimum energy conformer:


1. Use OBBuilder to create a 3D structure using rules and ring templates


2. Do 250 steps of a steepest descent geometry optimization with the
MMFF94 forcefield


3. Do 200 iterations of a Weighted Rotor conformational search
(optimizing each conformer with 25 steps of a steepest descent)


4. Do 250 steps of a conjugate gradient geometry optimization.

Warning from openbabel docs:
For many applications where 100s if not 1000s of molecules need to be
processed, gen3d is rather SLOW. Sometimes this function can cause a
segmentation fault.
A future version of Open Babel will provide options for slow/medium/fast
3D structure generation which will involve different compromises
between speed and finding the global energy minimum.


#### localopt(forcefield='mmff94', steps=500)
A wrapper to pybel’s localopt method to optimize a Molecule.


* **Parameters**


    * **forcefield** – Default is mmff94. Options are ‘gaff’, ‘ghemical’,
    ‘mmff94’, ‘mmff94s’, and ‘uff’.


    * **steps** – Default is 500.



#### make3d(forcefield='mmff94', steps=50)
A wrapper to pybel’s make3D method generate a 3D structure from a
2D or 0D structure.
The 3D structure is made very quickly using a combination of rules
(e.g. sp3 atoms should have four bonds arranged in a tetrahedron) and
ring templates (e.g. cyclohexane is shaped like a chair). Once 3D
coordinates are generated, hydrogens are added and a quick local
optimization is carried out as default.

The generated 3D structure can have clashes or have high energy
structures due to some strain. Please consider to use the conformer
search or geometry optimization to further optimize the structure.


* **Parameters**


    * **forcefield** – Default is mmff94. Options are ‘gaff’, ‘ghemical’,
    ‘mmff94’, ‘mmff94s’, and ‘uff’.


    * **steps** – Default is 50.



#### _property_ openbabel_mol()
Returns OpenBabel’s OBMol.


#### _property_ pybel_mol()
Returns Pybel’s Molecule object.


#### _property_ pymatgen_mol()
Returns pymatgen Molecule object.


#### remove_bond(idx1, idx2)
Remove a bond from an openbabel molecule.


* **Parameters**


    * **idx1** – The atom index of one of the atoms participating the in bond


    * **idx2** – The atom index of the other atom participating in the bond



#### rotor_conformer(\*rotor_args, algo='WeightedRotorSearch', forcefield='mmff94')
Conformer search based on several Rotor Search algorithms of openbabel.
If the input molecule is not 3D, make3d will be called (generate 3D
structure, add hydrogen, a quick localopt). All hydrogen atoms need
to be made explicit.


* **Parameters**


    * **rotor_args** – pass args to Rotor Search in openbabel.
    for “WeightedRotorSearch”: (conformers, geomSteps,
    sampleRingBonds-default False)
    for “SystematicRotorSearch”: (geomSteps-default 2500,
    sampleRingBonds-default False)
    for “RandomRotorSearch”: (conformers, geomSteps-default 2500,
    sampleRingBonds-default False)


    * **algo** (*str*) – Default is “WeightedRotorSearch”. Options are
    “SystematicRotorSearch”, “RandomRotorSearch”, and
    “WeightedRotorSearch”.


    * **forcefield** (*str*) – Default is mmff94. Options are ‘gaff’, ‘ghemical’,
    ‘mmff94’, ‘mmff94s’, and ‘uff’.



#### write_file(filename, file_format='xyz')
Uses OpenBabel to output all supported formats.


* **Parameters**


    * **filename** – Filename of file to output


    * **file_format** – String specifying any OpenBabel supported formats.


## pymatgen.io.cif module

Wrapper classes for Cif input and output from Structures.


### _class_ CifBlock(data, loops, header)
Bases: `object`

Object for storing cif data. All data is stored in a single dictionary.
Data inside loops are stored in lists in the data dictionary, and
information on which keys are grouped together are stored in the loops
attribute.


* **Parameters**


    * **data** – dict of data to go into the cif. Values should be convertible to string,
    or lists of these if the key is in a loop


    * **loops** – list of lists of keys, grouped by which loop they should appear in


    * **header** – name of the block (appears after the

    ```
    data_
    ```

     on the first line).



#### _format_field(v)

#### _loop_to_string(loop)

#### _classmethod_ _process_string(string)

#### _classmethod_ from_str(string)
Reads CifBlock from string.


* **Parameters**

    **string** – String representation.



* **Returns**

    CifBlock



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### maxlen(_ = 7_ )

### _class_ CifFile(data, orig_string=None, comment=None)
Bases: `object`

Reads and parses CifBlocks from a .cif file or string.


* **Parameters**


    * **data** (*dict*) – Of CifBlock objects.


    * **orig_string** (*str*) – The original cif string.


    * **comment** (*str*) – Comment string.



#### _classmethod_ from_file(filename)
Reads CifFile from a filename.


* **Parameters**

    **filename** – Filename



* **Returns**

    CifFile



#### _classmethod_ from_str(string)
Reads CifFile from a string.


* **Parameters**

    **string** – String representation.



* **Returns**

    CifFile



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


### _class_ CifParser(filename: str | StringIO, occupancy_tolerance: float = 1.0, site_tolerance: float = 0.0001, frac_tolerance: float = 0.0001)
Bases: `object`

Parses a CIF file. Attempts to fix CIFs that are out-of-spec, but will
issue warnings if corrections applied. These are also stored in the
CifParser’s errors attribute.


* **Parameters**


    * **filename** (*str*) – CIF filename, gzipped or bzipped CIF files are fine too.


    * **occupancy_tolerance** (*float*) – If total occupancy of a site is between 1 and occupancy_tolerance, the
    occupancies will be scaled down to 1.


    * **site_tolerance** (*float*) – This tolerance is used to determine if two sites are sitting in the same position,
    in which case they will be combined to a single disordered site. Defaults to 1e-4.


    * **frac_tolerance** (*float*) – This tolerance is used to determine is a coordinate should be rounded to an ideal
    value. E.g., 0.6667 is rounded to 2/3. This is desired if symmetry operations are going to be applied.
    However, for very large CIF files, this may need to be set to 0.



#### _get_structure(data: dict[str, Any], primitive: bool, symmetrized: bool, check_occu: bool = False)
Generate structure from part of the cif.


#### _parse_symbol(sym)
Parse a string with a symbol to extract a string representing an element.


* **Parameters**

    **sym** (*str*) – A symbol to be parsed.



* **Returns**

    A string with the parsed symbol. None if no parsing was possible.



#### _sanitize_data(data)
Some CIF files do not conform to spec. This function corrects
known issues, particular in regards to Springer materials/
Pauling files.

This function is here so that CifParser can assume its
input conforms to spec, simplifying its implementation.
:param data: CifBlock


* **Returns**

    data CifBlock



#### _unique_coords(coords: list[Vector3D], magmoms: list[[Magmom](pymatgen.electronic_structure.md#pymatgen.electronic_structure.core.Magmom)] | None = None, lattice: [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice) | None = None, labels: dict[Vector3D, str] | None = None)
Generate unique coordinates using coord and symmetry positions
and also their corresponding magnetic moments, if supplied.


#### as_dict()
MSONable dict


#### _static_ from_str(cif_string: str, \*\*kwargs)
Creates a CifParser from a string.


* **Parameters**


    * **cif_string** (*str*) – String representation of a CIF.


    * **\*\*kwargs** – Passthrough of all kwargs supported by CifParser.



* **Returns**

    CifParser



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_bibtex_string()
Get BibTeX reference from CIF file.
:param data:


* **Returns**

    BibTeX string.



#### get_lattice(data, length_strings=('a', 'b', 'c'), angle_strings=('alpha', 'beta', 'gamma'), lattice_type=None)
Generate the lattice from the provided lattice parameters. In
the absence of all six lattice parameters, the crystal system
and necessary parameters are parsed.


#### _static_ get_lattice_no_exception(data, length_strings=('a', 'b', 'c'), angle_strings=('alpha', 'beta', 'gamma'), lattice_type=None)
Take a dictionary of CIF data and returns a pymatgen Lattice object.


* **Parameters**


    * **data** – a dictionary of the CIF file


    * **length_strings** – The strings that are used to identify the length parameters in the CIF file.


    * **angle_strings** – The strings that are used to identify the angles in the CIF file.


    * **lattice_type** – The type of lattice.  This is a string, and can be any of the following:



* **Returns**

    Lattice object



#### get_magsymops(data)
Equivalent to get_symops except for magnetic symmetry groups.
Separate function since additional operation for time reversal symmetry
(which changes magnetic moments on sites) needs to be returned.


#### get_structures(primitive: bool = True, symmetrized: bool = False, check_occu: bool = True, on_error: Literal['ignore', 'warn', 'raise'] = 'warn')
Return list of structures in CIF file.


* **Parameters**


    * **primitive** (*bool*) – Set to False to return conventional unit cells.
    Defaults to True. With magnetic CIF files, will return primitive
    magnetic cell which may be larger than nuclear primitive cell.


    * **symmetrized** (*bool*) – If True, return a SymmetrizedStructure which will
    include the equivalent indices and symmetry operations used to
    create the Structure as provided by the CIF (if explicit symmetry
    operations are included in the CIF) or generated from information
    in the CIF (if only space group labels are provided). Note that
    currently Wyckoff labels and space group labels or numbers are
    not included in the generated SymmetrizedStructure, these will be
    notated as “Not Parsed” or -1 respectively.


    * **check_occu** (*bool*) – If False, site occupancy will not be checked, allowing unphysical
    occupancy != 1. Useful for experimental results in which occupancy was allowed
    to refine to unphysical values. Warning: unphysical site occupancies are incompatible
    with many pymatgen features. Defaults to True.


    * **on_error** (*'ignore'** | **'warn'** | **'raise'*) – What to do in case of KeyError or ValueError
    while parsing CIF file. Defaults to ‘warn’.



* **Returns**

    All structures in CIF file.



* **Return type**

    list[[Structure](pymatgen.core.md#pymatgen.core.structure.Structure)]



#### get_symops(data)
In order to generate symmetry equivalent positions, the symmetry
operations are parsed. If the symops are not present, the space
group symbol is parsed, and symops are generated.


#### _property_ has_errors()
Whether there are errors/warnings detected in CIF parsing.


#### _static_ parse_magmoms(data, lattice=None)
Parse atomic magnetic moments from data dictionary.


#### _static_ parse_oxi_states(data)
Parse oxidation states from data dictionary.


### _class_ CifWriter(struct, symprec=None, write_magmoms=False, significant_figures=8, angle_tolerance=5.0, refine_struct=True)
Bases: `object`

A wrapper around CifFile to write CIF files from pymatgen structures.


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure to write


    * **symprec** (*float*) – If not none, finds the symmetry of the structure
    and writes the cif with symmetry information. Passes symprec
    to the SpacegroupAnalyzer. See also refine_struct.


    * **write_magmoms** (*bool*) – If True, will write magCIF file. Incompatible
    with symprec


    * **significant_figures** (*int*) – Specifies precision for formatting of floats.
    Defaults to 8.


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding. Passes
    angle_tolerance to the SpacegroupAnalyzer. Used only if symprec
    is not None.


    * **refine_struct** – Used only if symprec is not None. If True, get_refined_structure
    is invoked to convert input structure from primitive to conventional.



#### _property_ ciffile()
CifFile associated with the CifWriter.


* **Type**

    Returns



#### write_file(filename)
Write the cif file.


### str2float(text)
Remove uncertainty brackets from strings and return the float.

## pymatgen.io.common module

Module for defining common data used and produced by atomistic simulation packages.


### _class_ VolumetricData(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), data, distance_matrix=None, data_aug=None)
Bases: `MSONable`

Simple volumetric object. Used to read LOCPOT/CHGCAR files produced by
vasp as well as cube files produced by other codes.


#### structure()
Structure associated with the Volumetric Data object.


* **Type**

    [Structure](pymatgen.core.md#pymatgen.core.structure.Structure)



#### is_spin_polarized()
True if run is spin polarized.


* **Type**

    bool



#### dim()
Tuple of dimensions of volumetric grid in each direction (nx, ny, nz).


* **Type**

    tuple



#### data()
Actual data as a dict of {string: np.array}. The string are “total”
and “diff”, in accordance to the output format of Vasp LOCPOT and
CHGCAR files where the total spin density is written first, followed
by the difference spin density.


* **Type**

    dict



#### ngridpts()
Total number of grid points in volumetric data.


* **Type**

    int


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
> f[“lattice”]: Lattice in the pymatgen.core.Lattice matrix

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


## pymatgen.io.core module

This module defines the abstract interface for reading and writing calculation
inputs in pymatgen. The interface comprises a 3-tiered hierarchy of classes.


1. An InputFile object represents the contents of a single input file, e.g.
the INCAR. This class standardizes file read and write operations.


2. An InputSet is a dict-like container that maps filenames (keys) to file
contents (either strings or InputFile objects). This class provides a standard
write_input() method.


3. InputGenerator classes implement a get_input_set method that, when provided
with a structure, return an InputSet object with all parameters set correctly.
Calculation input files can be written to disk with the write_inputs method.

If you want to implement a new InputGenerator, please take note of the following:


1. You must implement a get_input_set method that returns an InputSet


2. All customization of calculation parameters should be done in the __init__
method of the InputGenerator. The idea is that the generator contains
the “recipe”, but nothing that is specific to a particular system. get_input_set
takes system-specific information (such as structure) and applies the recipe.


3. All InputGenerator must save all supplied args and kwargs as instance variables.
E.g., self.my_arg = my_arg and self.kwargs = kwargs in the __init__. This
ensures the as_dict and from_dict work correctly.


### _class_ InputFile()
Bases: `MSONable`

Abstract base class to represent a single input file. Note that use of this class
is optional; it is possible create an InputSet that does not rely on underlying
InputFile objects.

All InputFile classes must implement a get_string method, which is called by
write_file.

If InputFile classes implement an __init__ method, they must assign all arguments
to __init__ as attributes.


#### _classmethod_ from_file(path: str | Path)
Creates an InputFile object from a file.


* **Parameters**

    **path** – Filename to read, including path.



* **Returns**

    InputFile



#### _abstract classmethod_ from_str(contents: str)
Create an InputFile object from a string.


* **Parameters**

    **contents** – The contents of the file as a single string



* **Returns**

    InputFile



#### _abstract classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead

Create an InputFile object from a string.


* **Parameters**

    **contents** – The contents of the file as a single string



* **Returns**

    InputFile



#### _abstract_ get_str()
Return a string representation of an entire input file.


#### _abstract_ get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead

Return a string representation of an entire input file.


#### write_file(filename: str | Path)
Write the input file.


* **Parameters**

    **filename** – The filename to output to, including path.



### _class_ InputGenerator()
Bases: `MSONable`

InputGenerator classes serve as generators for Input objects. They contain
settings or sets of instructions for how to create Input from a set of
coordinates or a previous calculation directory.


#### _abstract_ get_input_set()
Generate an InputSet object. Typically the first argument to this method
will be a Structure or other form of atomic coordinates.


### _class_ InputSet(inputs: dict[str | Path, str | InputFile] | None = None, \*\*kwargs)
Bases: `MSONable`, `MutableMapping`

Abstract base class for all InputSet classes. InputSet are dict-like
containers for all calculation input data.

Since InputSet inherits dict, it can be instantiated in the same manner,
or a custom __init__ can be provided. Either way, self should be
populated with keys that are filenames to be written, and values that are
InputFile objects or strings representing the entire contents of the file.

All InputSet must implement from_directory. Implementing the validate method
is optional.

Instantiate an InputSet.


* **Parameters**


    * **inputs** – The core mapping of filename: file contents that defines the InputSet data.
    This should be a dict where keys are filenames and values are InputFile objects
    or strings representing the entire contents of the file. If a value is not an
    InputFile object nor a str, but has a __str__ method, this str representation
    of the object will be written to the corresponding file. This mapping will
    become the .inputs attribute of the InputSet.


    * **\*\*kwargs** – Any kwargs passed will be set as class attributes e.g.
    InputSet(inputs={}, foo=’bar’) will make InputSet.foo == ‘bar’.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_directory(directory: str | Path)
Construct an InputSet from a directory of one or more files.


* **Parameters**

    **directory** – Directory to read input files from



#### validate()
A place to implement basic checks to verify the validity of an
input set. Can be as simple or as complex as desired.

Will raise a NotImplementedError unless overloaded by the inheriting class.


#### write_input(directory: str | Path, make_dir: bool = True, overwrite: bool = True, zip_inputs: bool = False)
Write Inputs to one or more files.


* **Parameters**


    * **directory** – Directory to write input files to


    * **make_dir** – Whether to create the directory if it does not already exist.


    * **overwrite** – Whether to overwrite an input file if it already exists.


    * **generate_inputs** (*Additional kwargs are passed to*) –


    * **zip_inputs** – If True, inputs will be zipped into a file with the
    same name as the InputSet (e.g., InputSet.zip)



### _exception_ ParseError()
Bases: `SyntaxError`

This exception indicates a problem was encountered during parsing due to unexpected formatting.

## pymatgen.io.cssr module

This module provides input and output from the CSSR file format.


### _class_ Cssr(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Bases: `object`

Basic object for working with Cssr file. Right now, only conversion from
a Structure to a Cssr file is supported.


* **Parameters**

    **structure** (*Structure/IStructure*) – A structure to create the Cssr object.



#### _static_ from_file(filename)
Reads a CSSR file to a Cssr object.


* **Parameters**

    **filename** (*str*) – Filename to read from.



* **Returns**

    Cssr object.



#### _static_ from_str(string)
Reads a string representation to a Cssr object.


* **Parameters**

    **string** (*str*) – A string representation of a CSSR.



* **Returns**

    Cssr object.



#### write_file(filename)
Write out a CSSR file.


* **Parameters**

    **filename** (*str*) – Filename to write to.


## pymatgen.io.fiesta module

This module implements input and output for Fiesta ([http://perso.neel.cnrs.fr/xavier.blase/fiesta/index.html](http://perso.neel.cnrs.fr/xavier.blase/fiesta/index.html)).

and

-Nwchem2Fiesta class: to create the input files needed for a Fiesta run
-Fiesta_run: run gw_fiesta and bse_fiesta
-Localised Basis set reader


### _class_ BSEOutput(filename)
Bases: `object`

A bse output file parser. The start…

All energies are in eV.


* **Parameters**

    **filename** – Filename to read.



#### _static_ _parse_job(output)

### _class_ BasisSetReader(filename)
Bases: `object`

A basis set reader.
Basis set are stored in data as a dict:
:key l_zeta_ng for each nl orbitals which contain list of tuple (alpha, coef) for each of the ng gaussians
in l_zeta orbital.


* **Parameters**

    **filename** – Filename to read.



#### _static_ _parse_file(input)

#### infos_on_basis_set()
Infos on the basis set as in Fiesta log.


#### set_n_nlmo()
the number of nlm orbitals for the basis set


### _class_ FiestaInput(mol, correlation_grid: dict[str, str] | None = None, Exc_DFT_option: dict[str, str] | None = None, COHSEX_options: dict[str, str] | None = None, GW_options: dict[str, str] | None = None, BSE_TDDFT_options: dict[str, str] | None = None)
Bases: `MSONable`

Input File for Fiesta called “cell.in” by default (mandatory in Fiesta for now).


* **Parameters**


    * **mol** – pymatgen mol


    * **correlation_grid** – dict


    * **Exc_DFT_option** – dict


    * **COHSEX_options** – dict


    * **GW_options** – dict


    * **BSE_TDDFT_options** – dict



#### as_dict()
MSONable dict


#### dump_BSE_data_in_GW_run(BSE_dump=True)

* **Parameters**

    **BSE_dump** – boolean



* **Returns**

    set the “do_bse” variable to one in cell.in



#### dump_TDDFT_data_in_GW_run(TDDFT_dump=True)

* **Parameters**

    **TDDFT_dump** – boolean



* **Returns**

    set the do_tddft variable to one in cell.in



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    FiestaInput



#### _classmethod_ from_file(filename)
Read an Fiesta input from a file. Currently tested to work with
files generated from this class itself.


* **Parameters**

    **filename** – Filename to parse.



* **Returns**

    FiestaInput object



#### _classmethod_ from_str(string_input)
Read an FiestaInput from a string. Currently tested to work with
files generated from this class itself.


* **Parameters**

    **string_input** – string_input to parse.



* **Returns**

    FiestaInput object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### _property_ infos_on_system()
Returns infos on initial parameters as in the log file of Fiesta.


#### _static_ make_FULL_BSE_Densities_folder(folder)
Mkdir “FULL_BSE_Densities” folder (needed for bse run) in the desired folder.


#### _property_ molecule()
Returns molecule associated with this FiestaInput.


#### set_BSE_options(n_excitations=10, nit_bse=200)
Set parameters in cell.in for a BSE computation
:param nv_bse: number of valence bands
:param nc_bse: number of conduction bands
:param n_excitations: number of excitations
:param nit_bse: number of iterations.


#### set_GW_options(nv_band=10, nc_band=10, n_iteration=5, n_grid=6, dE_grid=0.5)
Set parameters in cell.in for a GW computation
:param nv__band: number of valence bands to correct with GW
:param nc_band: number of conduction bands to correct with GW
:param n_iteration: number of iteration
:param n_grid and dE_grid:: number of points and spacing in eV for correlation grid.


#### set_auxiliary_basis_set(folder, auxiliary_folder, auxiliary_basis_set_type='aug_cc_pvtz')
copy in the desired folder the needed auxiliary basis set “X2.ion” where X is a specie.
:param auxiliary_folder: folder where the auxiliary basis sets are stored
:param auxiliary_basis_set_type: type of basis set (string to be found in the extension of the file name; must

> be in lower case). ex: C2.ion_aug_cc_pvtz_RI_Weigend find “aug_cc_pvtz”.


#### write_file(filename)
Write FiestaInput to a file
:param filename: Filename.


### _class_ FiestaOutput(filename)
Bases: `object`

A Fiesta output file parser.

All energies are in eV.


* **Parameters**

    **filename** – Filename to read.



#### _static_ _parse_job(output)

### _class_ FiestaRun(folder: str | None = None, grid: tuple[int, int, int] = (2, 2, 2), log_file: str = 'log')
Bases: `MSONable`

To run FIESTA inside python:

    if grid is [x,x] then bse runs
    if grid is [x,x,y] the fiesta(gw) runs
    otherwise it breaks.


* **Parameters**


    * **folder** – Folder to look for runs.


    * **grid** –


    * **log_file** – logfile of Fiesta.



#### _gw_run()
Performs FIESTA (gw) run.


#### as_dict()
MSONable dict


#### bse_run()
Performs BSE run.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    FiestaRun



#### run()
Performs FIESTA (gw) run.


### _class_ Nwchem2Fiesta(folder, filename='nwchem', log_file='log_n2f')
Bases: `MSONable`

To run NWCHEM2FIESTA inside python:

If nwchem.nw is the input, nwchem.out the output, and structure.movecs the
“movecs” file, the syntax to run NWCHEM2FIESTA is: NWCHEM2FIESTA
nwchem.nw  nwchem.nwout  structure.movecs > log_n2f

folder: where are stored the nwchem
filename: name of nwchem files read by NWCHEM2FIESTA (filename.nw, filename.nwout and filename.movecs)
logfile: logfile of NWCHEM2FIESTA.

the run method launches NWCHEM2FIESTA


#### as_dict()
MSONable dict


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    Nwchem2Fiesta



#### run()
Performs actual NWCHEM2FIESTA run.

## pymatgen.io.gaussian module

This module implements input and output processing from Gaussian.


### _class_ GaussianInput(mol, charge=None, spin_multiplicity=None, title=None, functional='HF', basis_set='6-31G(d)', route_parameters=None, input_parameters=None, link0_parameters=None, dieze_tag='#P', gen_basis=None)
Bases: `object`

An object representing a Gaussian input file.


* **Parameters**


    * **mol** – Input molecule. It can either be a Molecule object,
    a string giving the geometry in a format supported by Gaussian,
    or `None`. If the molecule is `None`, you will need to use
    read it in from a checkpoint. Consider adding `CHK` to the
    `link0_parameters`.


    * **charge** – Charge of the molecule. If None, charge on molecule is used.
    Defaults to None. This allows the input file to be set a
    charge independently from the molecule itself.
    If `mol` is not a Molecule object, then you must specify a charge.


    * **spin_multiplicity** – Spin multiplicity of molecule. Defaults to None,
    which means that the spin multiplicity is set to 1 if the
    molecule has no unpaired electrons and to 2 if there are
    unpaired electrons. If `mol` is not a Molecule object, then you

    > must specify the multiplicity



    * **title** – Title for run. Defaults to formula of molecule if None.


    * **functional** – Functional for run.


    * **basis_set** – Basis set for run.


    * **route_parameters** – Additional route parameters as a dict. For example,
    {‘SP’:””, “SCF”:”Tight”}


    * **input_parameters** – Additional input parameters for run as a dict. Used
    for example, in PCM calculations. E.g., {“EPS”:12}


    * **link0_parameters** – Link0 parameters as a dict. E.g., {“%mem”: “1000MW”}


    * **dieze_tag** – # preceding the route line. E.g. “#p”


    * **gen_basis** – allows a user-specified basis set to be used in a Gaussian
    calculation. If this is not None, the attribute `basis_set` will
    be set to “Gen”.



#### _static_ _parse_coords(coord_lines)
Helper method to parse coordinates.


#### _xyz_patt(_ = re.compile('^(\\\\w+)[\\\\s,]+([\\\\d\\\\.eE\\\\-]+)[\\\\s,]+([\\\\d\\\\.eE\\\\-]+)[\\\\s,]+([\\\\d\\\\.eE\\\\-]+)[\\\\-\\\\.\\\\s,\\\\w.]\*$'_ )

#### _zmat_patt(_ = re.compile('^(\\\\w+)\*([\\\\s,]+(\\\\w+)[\\\\s,]+(\\\\w+))\*[\\\\-\\\\.\\\\s,\\\\w]\*$'_ )

#### as_dict()
MSONable dict


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – dict



* **Returns**

    GaussianInput



#### _static_ from_file(filename)
Creates GaussianInput from a file.


* **Parameters**

    **filename** – Gaussian input filename



* **Returns**

    GaussianInput object



#### _static_ from_str(contents)
Creates GaussianInput from a string.


* **Parameters**

    **contents** – String representing an Gaussian input file.



* **Returns**

    GaussianInput object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_cart_coords()
Return the Cartesian coordinates of the molecule.


#### get_zmatrix()
Returns a z-matrix representation of the molecule.


#### _property_ molecule()
Returns molecule associated with this GaussianInput.


#### to_str(cart_coords=False)
Return GaussianInput string.


* **Parameters**

    **cart_coords** (*bool*) – If True, return Cartesian coordinates instead of z-matrix.
    Defaults to False.



#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead


#### write_file(filename, cart_coords=False)
Write the input string into a file.

Option: see __str__ method


### _class_ GaussianOutput(filename)
Bases: `object`

Parser for Gaussian output files.

Note: Still in early beta.


#### structures()
All structures from the calculation in the standard orientation. If the
symmetry is not considered, the standard orientation is not printed out
and the input orientation is used instead. Check the standard_orientation
attribute.


* **Type**

    list[[Structure](pymatgen.core.md#pymatgen.core.structure.Structure)]



#### structures_input_orientation()
All structures from the calculation in the input
orientation or the Z-matrix orientation (if an opt=z-matrix was requested).


* **Type**

    list



#### opt_structures()
All optimized structures from the calculation in the standard
orientation, if the attribute ‘standard_orientation’ is True, otherwise in the input
or the Z-matrix orientation.


* **Type**

    list



#### energies()
All energies from the calculation.


* **Type**

    list



#### eigenvalues()
List of eigenvalues for the last geometry.


* **Type**

    list



#### MO_coefficients()
Matrix of MO coefficients for the last geometry.


* **Type**

    list



#### cart_forces()
All Cartesian forces from the calculation.


* **Type**

    list



#### frequencies()
A list for each freq calculation and for each mode of a dict with
{

> > “frequency”: freq in cm-1,
> > “symmetry”: symmetry tag
> > “r_mass”: Reduce mass,
> > “f_constant”: force constant,
> > “IR_intensity”: IR Intensity,
> > “mode”: normal mode

> }

The normal mode is a 1D vector of dx, dy dz of each atom.


* **Type**

    list



#### hessian()
Matrix of second derivatives of the energy with respect to cartesian
coordinates in the input orientation frame. Need #P in the route section in order to
be in the output.


* **Type**

    ndarray



#### properly_terminated()
True if run has properly terminated.


* **Type**

    bool



#### is_pcm()
True if run is a PCM run.


* **Type**

    bool



#### is_spin()
True if it is an unrestricted run.


* **Type**

    bool



#### stationary_type()
If it is a relaxation run, indicates whether it is a minimum
(Minimum) or a saddle point (“Saddle”).


* **Type**

    str



#### corrections()
Thermochemical corrections if this run is a Freq run as a dict. Keys
are “Zero-point”, “Thermal”, “Enthalpy” and “Gibbs Free Energy”.


* **Type**

    dict



#### functional()
Functional used in the run.


* **Type**

    str



#### basis_set()
Basis set used in the run.


* **Type**

    str



#### route()
Additional route parameters as a dict. For example,
{‘SP’:””, “SCF”:”Tight”}.


* **Type**

    dict



#### dieze_tag()
# preceding the route line, e.g. “#P”.


* **Type**

    str



#### link0()
Link0 parameters as a dict. E.g., {“%mem”: “1000MW”}.


* **Type**

    dict



#### charge()
Charge for structure.


* **Type**

    int



#### spin_multiplicity()
Spin multiplicity for structure.


* **Type**

    int



#### num_basis_func()
Number of basis functions in the run.


* **Type**

    int



#### electrons()
Number of alpha and beta electrons as (N alpha, N beta).


* **Type**

    tuple



#### pcm()
PCM parameters and output if available.


* **Type**

    dict



#### errors()
Error if not properly terminated (list to be completed in error_defs).


* **Type**

    list



#### Mulliken_charges()
Mulliken atomic charges.


* **Type**

    list



#### eigenvectors()
Matrix of shape (num_basis_func, num_basis_func). Each column is an
eigenvectors and contains AO coefficients of an MO.
eigenvectors[Spin] = mat(num_basis_func, num_basis_func).


* **Type**

    dict



#### molecular_orbital()
MO development coefficients on AO in a more convenient array dict
for each atom and basis set label.
mo[Spin][OM j][atom i] = {AO_k: coeff, AO_k: coeff … }.


* **Type**

    dict



#### atom_basis_labels()
Labels of AO for each atoms. These labels are those used in the
output of molecular orbital coefficients (POP=Full) and in the molecular_orbital array
dict. atom_basis_labels[iatom] = [AO_k, AO_k, …].


* **Type**

    list



#### resumes()
List of gaussian data resume given at the end of the output file before
the quotation. The resumes are given as string.


* **Type**

    list



#### title()
Title of the gaussian run.


* **Type**

    str



#### standard_orientation()
If True, the geometries stored in the structures are in the
standard orientation. Else, the geometries are in the input orientation.


* **Type**

    bool



#### bond_orders()
Dict of bond order values read in the output file such as:
{(0, 1): 0.8709, (1, 6): 1.234, …}.
The keys are the atom indexes and the values are the Wiberg bond indexes that are
printed using pop=NBOREAD and $nbo bndidx $end.


* **Type**

    dict


Methods:
.. method:: to_input()

> Return a GaussianInput object using the last geometry and the same
> calculation parameters.


#### read_scan()
Read a potential energy surface from a gaussian scan calculation.


#### get_scan_plot()
Get a matplotlib plot of the potential energy surface


#### save_scan_plot()
Save a matplotlib plot of the potential energy surface to a file


* **Parameters**

    **filename** – Filename of Gaussian output file.



#### _check_pcm(line)

#### _parse(filename)

#### _parse_hessian(file, structure)
Parse the hessian matrix in the output file.


* **Parameters**


    * **file** – file object


    * **structure** – structure in the output file



#### as_dict()
JSON-serializable dict representation.


#### _property_ final_energy()
Final energy in Gaussian output.


#### _property_ final_structure()
Final structure in Gaussian output.


#### get_scan_plot(coords=None)
Get a matplotlib plot of the potential energy surface.


* **Parameters**

    **coords** – internal coordinate name to use as abscissa.



#### get_spectre_plot(sigma=0.05, step=0.01)
Get a matplotlib plot of the UV-visible xas. Transitions are plotted
as vertical lines and as a sum of normal functions with sigma with. The
broadening is applied in energy and the xas is plotted as a function
of the wavelength.


* **Parameters**


    * **sigma** – Full width at half maximum in eV for normal functions.


    * **step** – bin interval in eV



* **Returns**

    {“energies”: values, “lambda”: values, “xas”: values}

        where values are lists of abscissa (energies, lamba) and
        the sum of gaussian functions (xas).

    A matplotlib plot.




* **Return type**

    A dict



#### read_excitation_energies()
Read a excitation energies after a TD-DFT calculation.


* **Returns**

    A list of tuple for each transition such as

        [(energie (eV), lambda (nm), oscillatory strength), … ]




* **Return type**

    A list



#### read_scan()
Read a potential energy surface from a gaussian scan calculation.


* **Returns**

    {“energies”: [ values ],

        ”coords”: {“d1”: [ values ], “A2”, [ values ], … }}

    ”energies” are the energies of all points of the potential energy
    surface. “coords” are the internal coordinates used to compute the
    potential energy surface and the internal coordinates optimized,
    labelled by their name as defined in the calculation.




* **Return type**

    A dict



#### save_scan_plot(filename='scan.pdf', img_format='pdf', coords=None)
Save matplotlib plot of the potential energy surface to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **coords** – internal coordinate name to use as abcissa.



#### save_spectre_plot(filename='spectre.pdf', img_format='pdf', sigma=0.05, step=0.01)
Save matplotlib plot of the spectre to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **sigma** – Full width at half maximum in eV for normal functions.


    * **step** – bin interval in eV



#### to_input(mol=None, charge=None, spin_multiplicity=None, title=None, functional=None, basis_set=None, route_parameters=None, input_parameters=None, link0_parameters=None, dieze_tag=None, cart_coords=False)
Create a new input object using by default the last geometry read in
the output file and with the same calculation parameters. Arguments
are the same as GaussianInput class.


* **Returns**

    the gaussian input object



* **Return type**

    gaunip (GaussianInput)



### read_route_line(route)
read route line in gaussian input/output and return functional basis_set
and a dictionary of other route parameters.


* **Parameters**

    **route** (*str*) – the route line



* **Returns**

    the method (HF, PBE …)
    basis_set (str) : the basis set
    route (dict) : dictionary of parameters



* **Return type**

    functional (str)


## pymatgen.io.jarvis module

This module provides conversion between the JARVIS
Atoms object and pymatgen Structure objects.


### _class_ JarvisAtomsAdaptor()
Bases: `object`

Adaptor serves as a bridge between JARVIS Atoms and pymatgen objects.


#### _static_ get_atoms(structure)
Returns JARVIS Atoms object from pymatgen structure.


* **Parameters**

    **structure** – pymatgen.core.structure.Structure



* **Returns**

    JARVIS Atoms object



#### _static_ get_structure(atoms)
Returns pymatgen structure from JARVIS Atoms.


* **Parameters**

    **atoms** – JARVIS Atoms object



* **Returns**

    Equivalent pymatgen.core.structure.Structure


## pymatgen.io.lmto module

Module for implementing a CTRL file object class for the Stuttgart
LMTO-ASA code. It will primarily be used to generate a pymatgen
Structure object in the pymatgen.electronic_structure.cohp.py module.


### _class_ LMTOCopl(filename='COPL', to_eV=False)
Bases: `object`

Class for reading COPL files, which contain COHP data.


#### cohp_data()
Contains the COHP data of the form:
{bond: {“COHP”: {Spin.up: cohps, Spin.down:cohps},

> “ICOHP”: {Spin.up: icohps, Spin.down: icohps},
> “length”: bond length}


* **Type**

    dict



#### efermi()
The Fermi energy in Ry or eV.


* **Type**

    float



#### energies()
Sequence of energies in Ry or eV.


* **Type**

    list



#### is_spin_polarized()
Boolean to indicate if the calculation is spin polarized.


* **Type**

    bool



* **Parameters**


    * **filename** – filename of the COPL file. Defaults to “COPL”.


    * **to_eV** – LMTO-ASA gives energies in Ry. To convert energies into
    eV, set to True. Defaults to False for energies in Ry.



#### _static_ _get_bond_data(line)
Subroutine to extract bond label, site indices, and length from
a COPL header line. The site indices are zero-based, so they
can be easily used with a Structure object.

Example header line: Fe-1/Fe-1-tr(-1,-1,-1) : 2.482 Ang.


* **Parameters**

    **line** – line in the COHPCAR header describing the bond.



* **Returns**

    The bond label, the bond length and a tuple of the site indices.



### _class_ LMTOCtrl(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), header: str | None = None, version: str = 'LMASA-47')
Bases: `object`

Class for parsing CTRL files from the Stuttgart LMTO-ASA code.
Currently, only HEADER, VERS and the structure can be used.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – pymatgen object.


    * **header** (*str*) – The header for the CTRL file. Defaults to None.


    * **version** (*str*) – The LMTO version that is used for the VERS category.
    Defaults to version (4.7).



#### as_dict()
Returns the CTRL as a dictionary. “SITE” and “CLASS” are of
the form {‘CATEGORY’: {‘TOKEN’: value}}, the rest is of the
form ‘TOKEN’/’CATEGORY’: value. It gets the conventional standard
structure because primitive cells use the conventional
a-lattice parameter as the scaling factor and not the a-lattice
parameter of the primitive cell.


#### _classmethod_ from_dict(dct)
Creates a CTRL file object from a dictionary. The dictionary
must contain the items “ALAT”, PLAT” and “SITE”.

Valid dictionary items are:

    ALAT: the a-lattice parameter
    PLAT: (3x3) array for the lattice vectors
    SITE: list of dictionaries: {‘ATOM’: class label, ‘POS’: (3x1) array of fractional coordinates}
    CLASS (optional): list of unique atom labels as str
    SPCGRP (optional): space group symbol (str) or number (int)
    HEADER (optional): HEADER text as a str
    VERS (optional): LMTO version as a str


* **Parameters**

    **dct** – The CTRL file as a dictionary.



* **Returns**

    An LMTOCtrl object.



#### _classmethod_ from_file(filename='CTRL', \*\*kwargs)
Creates a CTRL file object from an existing file.


* **Parameters**

    **filename** – The name of the CTRL file. Defaults to ‘CTRL’.



* **Returns**

    An LMTOCtrl object.



#### _classmethod_ from_str(data: str, sigfigs: int = 8)
Creates a CTRL file object from a string. This will mostly be
used to read an LMTOCtrl object from a CTRL file. Empty spheres
are ignored.


* **Parameters**

    **data** (*str*) – String representation of the CTRL file.



* **Returns**

    An LMTOCtrl object.



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_str(sigfigs=8)
Generates the string representation of the CTRL file. This is
the minimal CTRL file necessary to execute lmhart.run.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### write_file(filename='CTRL', \*\*kwargs)
Writes a CTRL file with structure, HEADER, and VERS that can be
used as input for lmhart.run.

## pymatgen.io.nwchem module

This module implements input and output processing from Nwchem.

2015/09/21 - Xin Chen ([chenxin13@mails.tsinghua.edu.cn](mailto:chenxin13@mails.tsinghua.edu.cn)):

> NwOutput will read new kinds of data:

> >
> > 1. normal hessian matrix.       [“hessian”]


> > 2. projected hessian matrix.    [“projected_hessian”]


> > 3. normal frequencies.          [“normal_frequencies”]

> For backward compatibility, the key for accessing the projected frequencies
> is still ‘frequencies’.

2015/10/12 - Xin Chen

    NwOutput will read new kinds of data:

    >
    > 1. forces.                      [“forces”]


### _class_ NwInput(mol, tasks, directives=None, geometry_options=('units', 'angstroms'), symmetry_options=None, memory_options=None)
Bases: `MSONable`

An object representing a Nwchem input file, which is essentially a list
of tasks on a particular molecule.


* **Parameters**


    * **mol** – Input molecule. If molecule is a single string, it is used as a
    direct input to the geometry section of the Gaussian input
    file.


    * **tasks** – List of NwTasks.


    * **directives** – List of root level directives as tuple. E.g.,
    [(“start”, “water”), (“print”, “high”)]


    * **geometry_options** – Additional list of options to be supplied to the
    geometry. E.g., [“units”, “angstroms”, “noautoz”]. Defaults to
    (“units”, “angstroms”).


    * **symmetry_options** – Addition list of option to be supplied to the
    symmetry. E.g. [“c1”] to turn off the symmetry


    * **memory_options** – Memory controlling options. str.
    E.g “total 1000 mb stack 400 mb”.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    NwInput



#### _classmethod_ from_file(filename)
Read an NwInput from a file. Currently tested to work with
files generated from this class itself.


* **Parameters**

    **filename** – Filename to parse.



* **Returns**

    NwInput object



#### _classmethod_ from_str(string_input)
Read an NwInput from a string. Currently tested to work with
files generated from this class itself.


* **Parameters**

    **string_input** – string_input to parse.



* **Returns**

    NwInput object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### _property_ molecule()
Returns molecule associated with this GaussianInput.


#### write_file(filename)

* **Parameters**

    **filename** (*str*) – Filename.



### _exception_ NwInputError()
Bases: `Exception`

Error class for NwInput.


### _class_ NwOutput(filename)
Bases: `object`

A Nwchem output file parser. Very basic for now - supports only dft and
only parses energies and geometries. Please note that Nwchem typically
outputs energies in either au or kJ/mol. All energies are converted to
eV in the parser.


* **Parameters**

    **filename** – Filename to read.



#### _static_ _parse_job(output)

#### _static_ _parse_preamble(preamble)

#### get_excitation_spectrum(width=0.1, npoints=2000)
Generate an excitation spectra from the singlet roots of TDDFT
calculations.


* **Parameters**


    * **width** (*float*) – Width for Gaussian smearing.


    * **npoints** (*int*) – Number of energy points. More points => smoother
    curve.



* **Returns**

    (ExcitationSpectrum) which can be plotted using

        pymatgen.vis.plotters.SpectrumPlotter.




#### parse_tddft()
Parses TDDFT roots. Adapted from nw_spectrum.py script.


* **Returns**

    {

        “singlet”: [

            {

                “energy”: float,
                “osc_strength: float

            }

        ],
        “triplet”: [

        > {

        >     “energy”: float

        > }

        ]

    }




### _class_ NwTask(charge, spin_multiplicity, basis_set, basis_set_option='cartesian', title=None, theory='dft', operation='optimize', theory_directives=None, alternate_directives=None)
Bases: `MSONable`

Base task for Nwchem.

Very flexible arguments to support many types of potential setups.
Users should use more friendly static methods unless they need the
flexibility.


* **Parameters**


    * **charge** – Charge of the molecule. If None, charge on molecule is
    used. Defaults to None. This allows the input file to be set a
    charge independently from the molecule itself.


    * **spin_multiplicity** – Spin multiplicity of molecule. Defaults to None,
    which means that the spin multiplicity is set to 1 if the
    molecule has no unpaired electrons and to 2 if there are
    unpaired electrons.


    * **basis_set** – The basis set used for the task as a dict. E.g.,
    {“C”: “6-311++G\*\*”, “H”: “6-31++G\*\*”}.


    * **basis_set_option** – cartesian (default) | spherical,


    * **title** – Title for the task. Defaults to None, which means a title
    based on the theory and operation of the task is
    autogenerated.


    * **theory** – The theory used for the task. Defaults to “dft”.


    * **operation** – The operation for the task. Defaults to “optimize”.


    * **theory_directives** – A dict of theory directives. For example,
    if you are running dft calculations, you may specify the
    exchange correlation functional using {“xc”: “b3lyp”}.


    * **alternate_directives** – A dict of alternate directives. For
    example, to perform cosmo calculations and dielectric
    constant of 78, you’d supply {‘cosmo’: {“dielectric”: 78}}.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ dft_task(mol, xc='b3lyp', \*\*kwargs)
A class method for quickly creating DFT tasks with optional
cosmo parameter .


* **Parameters**


    * **mol** – Input molecule


    * **xc** – Exchange correlation to use.


    * **kwargs** – Any of the other kwargs supported by NwTask. Note the
    theory is always “dft” for a dft task.



#### _classmethod_ esp_task(mol, \*\*kwargs)
A class method for quickly creating ESP tasks with RESP
charge fitting.


* **Parameters**


    * **mol** – Input molecule


    * **kwargs** – Any of the other kwargs supported by NwTask. Note the
    theory is always “dft” for a dft task.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    NwTask



#### _classmethod_ from_molecule(mol, theory, charge=None, spin_multiplicity=None, basis_set='6-31g', basis_set_option='cartesian', title=None, operation='optimize', theory_directives=None, alternate_directives=None)
Very flexible arguments to support many types of potential setups.
Users should use more friendly static methods unless they need the
flexibility.


* **Parameters**


    * **mol** – Input molecule


    * **charge** – Charge of the molecule. If None, charge on molecule is
    used. Defaults to None. This allows the input file to be set a
    charge independently from the molecule itself.


    * **spin_multiplicity** – Spin multiplicity of molecule. Defaults to None,
    which means that the spin multiplicity is set to 1 if the
    molecule has no unpaired electrons and to 2 if there are
    unpaired electrons.


    * **basis_set** – The basis set to be used as string or a dict. E.g.,
    {“C”: “6-311++G\*\*”, “H”: “6-31++G\*\*”} or “6-31G”. If string,
    same basis set is used for all elements.


    * **basis_set_option** – cartesian (default) | spherical,


    * **title** – Title for the task. Defaults to None, which means a title
    based on the theory and operation of the task is
    autogenerated.


    * **theory** – The theory used for the task. Defaults to “dft”.


    * **operation** – The operation for the task. Defaults to “optimize”.


    * **theory_directives** – A dict of theory directives. For example,
    if you are running dft calculations, you may specify the
    exchange correlation functional using {“xc”: “b3lyp”}.


    * **alternate_directives** – A dict of alternate directives. For
    example, to perform cosmo calculations with DFT, you’d supply
    {‘cosmo’: “cosmo”}.



#### operations(_ = {'': 'dummy', 'dynamics': 'Perform classical molecular dynamics.', 'energy': 'Evaluate the single point energy.', 'freq': 'Same as frequencies.', 'frequencies': 'Compute second derivatives and print out an analysis of molecular vibrations.', 'gradient': 'Evaluate the derivative of the energy with respect to nuclear coordinates.', 'hessian': 'Compute second derivatives.', 'optimize': 'Minimize the energy by varying the molecular structure.', 'property': 'Calculate the properties for the wave function.', 'saddle': 'Conduct a search for a transition state (or saddle point).', 'thermodynamics': 'Perform multi-configuration thermodynamic integration using classical MD.', 'vscf': 'Compute anharmonic contributions to the vibrational modes.'_ )

#### theories(_ = {'band': 'Pseudopotential plane-wave DFT for solids using NWPW', 'ccsd': 'Coupled-cluster single and double excitations', 'ccsd(t)': 'Coupled-cluster linearized triples approximation', 'ccsd+t(ccsd)': 'Fourth order triples contribution', 'dft': 'DFT', 'direct_mp2': 'MP2 using a full-direct algorithm', 'esp': 'ESP', 'g3gn': 'some description', 'mcscf': 'Multiconfiguration SCF', 'md': 'Classical molecular dynamics simulation', 'mp2': 'MP2 using a semi-direct algorithm', 'pspw': 'Pseudopotential plane-wave DFT for molecules and insulating solids using NWPW', 'rimp2': 'MP2 using the RI approximation', 'scf': 'Hartree-Fock', 'selci': 'Selected CI with perturbation correction', 'sodft': 'Spin-Orbit DFT', 'tce': 'Tensor Contraction Engine', 'tddft': 'Time Dependent DFT'_ )
## pymatgen.io.packmol module

This module provides a pymatgen I/O interface to packmol.

This adopts the minimal core I/O interface (see pymatgen/io/core).
In this case, only a two classes are used. PackmolSet(InputSet) is the container
class that provides a run() method for running packmol locally.

PackmolBoxGen(InputGenerator) provides a recipe for packing molecules into a
box, which returns a PackmolSet object.

For the run() method to work, you need to install the packmol package
See [http://m3g.iqm.unicamp.br/packmol](http://m3g.iqm.unicamp.br/packmol) or
[http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml](http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml)
for download and setup instructions. Note that packmol versions prior to 20.3.0
do not support paths with spaces.
After installation, you may need to manually add the path of the packmol
executable to the PATH environment variable.


### _class_ PackmolBoxGen(tolerance: float = 2.0, seed: int = 1, control_params: dict | None = None, inputfile: str | Path = 'packmol.inp', outputfile: str | Path = 'packmol_out.xyz', stdoutfile: str | Path = 'packmol.stdout')
Bases: `InputGenerator`

Generator for a Packmol InputSet that packs one or more molecules into a rectangular
simulation box.

Instantiate a PackmolBoxGen class. The init method defines simulations parameters
like filenames, random seed, tolerance, etc.


* **Parameters**


    * **tolerance** – Tolerance for packmol, in Å.


    * **seed** – Random seed for packmol. Use a value of 1 (default) for deterministic
    output, or -1 to generate a new random seed from the current time.


    * **inputfile** – Path to the input file. Default to ‘packmol.inp’.


    * **outputfile** – Path to the output file. Default to ‘output.xyz’.


    * **stdoutfile** – Path to the file where stdout will be recorded. Default to ‘packmol.stdout’



#### get_input_set(molecules: list[dict], box: list[float] | None = None)
Generate a Packmol InputSet for a set of molecules.


* **Parameters**

    **molecules** – A list of dict containing information about molecules to pack
    into the box. Each dict requires three keys:

    >
    > 1. ”name” - the structure name


    > 2. ”number” - the number of that molecule to pack into the box


    > 3. ”coords” - Coordinates in the form of either a Molecule object or

    >     a path to a file.



### Example

{“name”: “water”,

    “number”: 500,
    “coords”: “/path/to/input/file.xyz”}

    > box: A list of box dimensions xlo, ylo, zlo, xhi, yhi, zhi, in Å. If set to None

    >     (default), pymatgen will estimate the required box size based on the volumes of
    >     the provided molecules.


### _class_ PackmolSet(inputs: dict[str | Path, str | InputFile] | None = None, \*\*kwargs)
Bases: `InputSet`

InputSet for the Packmol software. This class defines several attributes related to.

Instantiate an InputSet.


* **Parameters**


    * **inputs** – The core mapping of filename: file contents that defines the InputSet data.
    This should be a dict where keys are filenames and values are InputFile objects
    or strings representing the entire contents of the file. If a value is not an
    InputFile object nor a str, but has a __str__ method, this str representation
    of the object will be written to the corresponding file. This mapping will
    become the .inputs attribute of the InputSet.


    * **\*\*kwargs** – Any kwargs passed will be set as class attributes e.g.
    InputSet(inputs={}, foo=’bar’) will make InputSet.foo == ‘bar’.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_directory(directory: str | Path)
Construct an InputSet from a directory of one or more files.


* **Parameters**

    **directory** (*str** | **Path*) – Directory to read input files from.



#### run(path: str | Path, timeout=30)
Run packmol and write out the packed structure.


* **Parameters**


    * **path** – The path in which packmol input files are located.


    * **timeout** – Timeout in seconds.



* **Raises**


    * **ValueError if packmol does not succeed in packing the box.** –


    * **TimeoutExpiredError if packmold does not finish within the timeout.** –


## pymatgen.io.phonopy module

Module for interfacing with phonopy, see [https://atztogo.github.io/phonopy/](https://atztogo.github.io/phonopy/).


### _extrapolate_grun(b, distance, gruneisenparameter, gruneisenband, i, pa)

### eigvec_to_eigdispl(v, q, frac_coords, mass)
Converts a single eigenvector to an eigendisplacement in the primitive cell
according to the formula:

```default
exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
```

Compared to the modulation option in phonopy, here all the additional
multiplicative and phase factors are set to 1.


* **Parameters**


    * **v** – the vector that should be converted. A 3D complex numpy array.


    * **q** – the q point in fractional coordinates


    * **frac_coords** – the fractional coordinates of the atom


    * **mass** – the mass of the atom



### get_complete_ph_dos(partial_dos_path, phonopy_yaml_path)
Creates a pymatgen CompletePhononDos from a partial_dos.dat and
phonopy.yaml files.
The second is produced when generating a Dos and is needed to extract
the structure.


* **Parameters**


    * **partial_dos_path** – path to the partial_dos.dat file.


    * **phonopy_yaml_path** – path to the phonopy.yaml file.



### get_displaced_structures(pmg_structure, atom_disp=0.01, supercell_matrix=None, yaml_fname=None, \*\*kwargs)
Generate a set of symmetrically inequivalent displaced structures for
phonon calculations.


* **Parameters**


    * **pmg_structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – A pymatgen structure object.


    * **atom_disp** (*float*) – Atomic displacement. Default is 0.01 $\\AA$.


    * **supercell_matrix** (*3x3 array*) – Scaling matrix for supercell.


    * **yaml_fname** (*str*) – If not None, it represents the full path to
    the outputting displacement yaml file, e.g. disp.yaml.


    * **\*\*kwargs** – Parameters used in Phonopy.generate_displacement method.



* **Returns**

    A list of symmetrically inequivalent structures with displacements, in
    which the first element is the perfect supercell structure.



### get_gruneisen_ph_bs_symm_line(gruneisen_path, structure=None, structure_path=None, labels_dict=None, fit=False)
Creates a pymatgen GruneisenPhononBandStructure from a band.yaml file.
The labels will be extracted from the dictionary, if present.
If the ‘eigenvector’ key is found the eigendisplacements will be
calculated according to the formula:
\\exp(2\*pi\*i\*(frac_coords \\dot q) / sqrt(mass) \* v

> and added to the object.


* **Parameters**


    * **gruneisen_path** – path to the band.yaml file


    * **structure** – pymaten Structure object


    * **structure_path** – path to a structure file (e.g., POSCAR)


    * **labels_dict** – dict that links a qpoint in frac coords to a label.


    * **fit** – Substitute Grueneisen parameters close to the gamma point
    with points obtained from a fit to a spline if the derivate from
    a smooth curve (i.e. if the slope changes by more than 200% in the
    range of 10% around the gamma point).
    These derivations occur because of very small frequencies
    (and therefore numerical inaccuracies) close to gamma.



### get_gruneisenparameter(gruneisen_path, structure=None, structure_path=None)
Get Gruneisen object from gruneisen.yaml file, as obtained from phonopy (Frequencies in THz!).
The order is structure > structure path > structure from gruneisen dict.
Newer versions of phonopy include the structure in the yaml file,
the structure/structure_path is kept for compatibility.


* **Parameters**


    * **gruneisen_path** – Path to gruneisen.yaml file (frequencies have to be in THz!)


    * **structure** – pymatgen Structure object


    * **structure_path** – path to structure in a file (e.g., POSCAR)



* **Returns**

    GruneisenParameter



### get_gs_ph_bs_symm_line_from_dict(gruneisen_dict, structure=None, structure_path=None, labels_dict=None, fit=False)
Creates a pymatgen GruneisenPhononBandStructure object from the dictionary
extracted by the gruneisen.yaml file produced by phonopy. The labels
will be extracted from the dictionary, if present. If the ‘eigenvector’
key is found the eigendisplacements will be calculated according to the
formula:

```default
exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
```

and added to the object. A fit algorithm can be used to replace diverging
Gruneisen values close to gamma.


* **Parameters**


    * **gruneisen_dict** (*dict*) – the dictionary extracted from the gruneisen.yaml file


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – pymatgen structure object


    * **structure_path** – path to structure file


    * **labels_dict** (*dict*) – dict that links a qpoint in frac coords to a label.
    Its value will replace the data contained in the band.yaml.


    * **fit** (*bool*) – Substitute Grueneisen parameters close to the gamma point
    with points obtained from a fit to a spline if the derivate from
    a smooth curve (i.e. if the slope changes by more than 200% in the
    range of 10% around the gamma point).
    These derivations occur because of very small frequencies
    (and therefore numerical inaccuracies) close to gamma.



### get_ph_bs_symm_line(bands_path, has_nac=False, labels_dict=None)
Creates a pymatgen PhononBandStructure from a band.yaml file.
The labels will be extracted from the dictionary, if present.
If the ‘eigenvector’  key is found the eigendisplacements will be
calculated according to the formula:
\\exp(2\*pi\*i\*(frac_coords \\dot q) / sqrt(mass) \* v

> and added to the object.


* **Parameters**


    * **bands_path** – path to the band.yaml file


    * **has_nac** – True if the data have been obtained with the option
    –nac option. Default False.


    * **labels_dict** – dict that links a qpoint in frac coords to a label.



### get_ph_bs_symm_line_from_dict(bands_dict, has_nac=False, labels_dict=None)
Creates a pymatgen PhononBandStructure object from the dictionary
extracted by the band.yaml file produced by phonopy. The labels
will be extracted from the dictionary, if present. If the ‘eigenvector’
key is found the eigendisplacements will be calculated according to the
formula:

```default
exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
```

and added to the object.


* **Parameters**


    * **bands_dict** – the dictionary extracted from the band.yaml file


    * **has_nac** – True if the data have been obtained with the option
    –nac option. Default False.


    * **labels_dict** – dict that links a qpoint in frac coords to a label.
    Its value will replace the data contained in the band.yaml.



### get_ph_dos(total_dos_path)
Creates a pymatgen PhononDos from a total_dos.dat file.


* **Parameters**

    **total_dos_path** – path to the total_dos.dat file.



### get_phonon_band_structure_from_fc(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), supercell_matrix: ndarray, force_constants: ndarray, mesh_density: float = 100.0, \*\*kwargs)
Get a uniform phonon band structure from phonopy force constants.


* **Parameters**


    * **structure** – A structure.


    * **supercell_matrix** – The supercell matrix used to generate the force
    constants.


    * **force_constants** – The force constants in phonopy format.


    * **mesh_density** – The density of the q-point mesh. See the docstring
    for the `mesh` argument in Phonopy.init_mesh() for more details.


    * **\*\*kwargs** – Additional kwargs passed to the Phonopy constructor.



* **Returns**

    The uniform phonon band structure.



### get_phonon_band_structure_symm_line_from_fc(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), supercell_matrix: ndarray, force_constants: ndarray, line_density: float = 20.0, symprec: float = 0.01, \*\*kwargs)
Get a phonon band structure along a high symmetry path from phonopy force
constants.


* **Parameters**


    * **structure** – A structure.


    * **supercell_matrix** – The supercell matrix used to generate the force
    constants.


    * **force_constants** – The force constants in phonopy format.


    * **line_density** – The density along the high symmetry path.


    * **symprec** – Symmetry precision passed to phonopy and used for determining
    the band structure path.


    * **\*\*kwargs** – Additional kwargs passed to the Phonopy constructor.



* **Returns**

    The line mode band structure.



### get_phonon_dos_from_fc(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), supercell_matrix: ndarray, force_constants: ndarray, mesh_density: float = 100.0, num_dos_steps: int = 200, \*\*kwargs)
Get a projected phonon density of states from phonopy force constants.


* **Parameters**


    * **structure** – A structure.


    * **supercell_matrix** – The supercell matrix used to generate the force
    constants.


    * **force_constants** – The force constants in phonopy format.


    * **mesh_density** – The density of the q-point mesh. See the docstring
    for the `mesh` argument in Phonopy.init_mesh() for more details.


    * **num_dos_steps** – Number of frequency steps in the energy grid.


    * **\*\*kwargs** – Additional kwargs passed to the Phonopy constructor.



* **Returns**

    The density of states.



### get_phonopy_structure(pmg_structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Convert a pymatgen Structure object to a PhonopyAtoms object.


* **Parameters**

    **pmg_structure** (*pymatgen Structure*) – A Pymatgen structure object.



### get_pmg_structure(phonopy_structure: PhonopyAtoms)
Convert a PhonopyAtoms object to pymatgen Structure object.


* **Parameters**

    **phonopy_structure** (*PhonopyAtoms*) – A phonopy structure object.



### get_structure_from_dict(d)
Extracts a structure from the dictionary extracted from the output
files of phonopy like phonopy.yaml or band.yaml.
Adds “phonopy_masses” in the site_properties of the structures.
Compatible with older phonopy versions.


### get_thermal_displacement_matrices(thermal_displacements_yaml='thermal_displacement_matrices.yaml', structure_path='POSCAR')
Function to read “thermal_displacement_matrices.yaml” from phonopy and return a list of
ThermalDisplacementMatrices objects
:param thermal_displacements_yaml: path to thermal_displacement_matrices.yaml
:param structure_path: path to POSCAR.

Returns:

## pymatgen.io.prismatic module

Write Prismatic ([http://prism-em.com](http://prism-em.com)) input files.


### _class_ Prismatic(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), comment: str = 'Generated by pymatgen')
Bases: `object`

Class to write Prismatic  ([http://prism-em.com/](http://prism-em.com/)) input files.
This is designed for STEM image simulation.


* **Parameters**


    * **structure** – pymatgen Structure


    * **comment** (*str*) – comment.



#### to_str()

* **Returns**

    Prismatic XYZ file. This is similar to XYZ format

        but has specific requirements for extra fields, headers, etc.




* **Return type**

    str



#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead

## pymatgen.io.pwscf module

This module implements input and output processing from PWSCF.


### _class_ PWInput(structure, pseudo=None, control=None, system=None, electrons=None, ions=None, cell=None, kpoints_mode='automatic', kpoints_grid=(1, 1, 1), kpoints_shift=(0, 0, 0))
Bases: `object`

Base input file class. Right now, only supports no symmetry and is
very basic.

Initializes a PWSCF input file.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure. For spin-polarized calculation,
    properties (e.g. {“starting_magnetization”: -0.5,
    “pseudo”: “Mn.pbe-sp-van.UPF”}) on each site is needed instead of
    pseudo (dict).


    * **pseudo** (*dict*) – A dict of the pseudopotentials to use. Default to None.


    * **control** (*dict*) – Control parameters. Refer to official PWSCF doc
    on supported parameters. Default to {“calculation”: “scf”}


    * **system** (*dict*) – System parameters. Refer to official PWSCF doc
    on supported parameters. Default to None, which means {}.


    * **electrons** (*dict*) – Electron parameters. Refer to official PWSCF doc
    on supported parameters. Default to None, which means {}.


    * **ions** (*dict*) – Ions parameters. Refer to official PWSCF doc
    on supported parameters. Default to None, which means {}.


    * **cell** (*dict*) – Cell parameters. Refer to official PWSCF doc
    on supported parameters. Default to None, which means {}.


    * **kpoints_mode** (*str*) – Kpoints generation mode. Default to automatic.


    * **kpoints_grid** (*sequence*) – The kpoint grid. Default to (1, 1, 1).


    * **kpoints_shift** (*sequence*) – The shift for the kpoints. Defaults to
    (0, 0, 0).



#### as_dict()
Create a dictionary representation of a PWInput object.


* **Returns**

    dict



#### _classmethod_ from_dict(pwinput_dict)
Load a PWInput object from a dictionary.


* **Parameters**

    **pwinput_dict** (*dict*) – dictionary with PWInput data



* **Returns**

    PWInput object



#### _static_ from_file(filename)
Reads an PWInput object from a file.


* **Parameters**

    **filename** (*str*) – Filename for file



* **Returns**

    PWInput object



#### _static_ from_str(string)
Reads an PWInput object from a string.


* **Parameters**

    **string** (*str*) – PWInput string



* **Returns**

    PWInput object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### _static_ proc_val(key, val)
Static helper method to convert PWINPUT parameters to proper type, e.g.,
integers, floats, etc.


* **Parameters**


    * **key** – PWINPUT parameter key


    * **val** – Actual value of PWINPUT parameter.



#### write_file(filename)
Write the PWSCF input file.


* **Parameters**

    **filename** (*str*) – The string filename to output to.



### _exception_ PWInputError()
Bases: `BaseException`

Error for PWInput.


### _class_ PWOutput(filename)
Bases: `object`

Parser for PWSCF output file.


* **Parameters**

    **filename** (*str*) – Filename.



#### _property_ final_energy()
Final energy.


* **Type**

    Returns



#### get_celldm(idx: int)

* **Parameters**

    **idx** (*int*) – index.



* **Returns**

    Cell dimension along index



#### _property_ lattice_type()
Lattice type.


* **Type**

    Returns



#### patterns(_ = {'celldm1': 'celldm\\\\(1\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm2': 'celldm\\\\(2\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm3': 'celldm\\\\(3\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm4': 'celldm\\\\(4\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm5': 'celldm\\\\(5\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm6': 'celldm\\\\(6\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'ecut': 'kinetic\\\\-energy cutoff\\\\s+=\\\\s+([\\\\d\\\\.\\\\-]+)\\\\s+Ry', 'energies': 'total energy\\\\s+=\\\\s+([\\\\d\\\\.\\\\-]+)\\\\sRy', 'lattice_type': 'bravais\\\\-lattice index\\\\s+=\\\\s+(\\\\d+)', 'nkpts': 'number of k points=\\\\s+([\\\\d]+)'_ )

#### read_pattern(patterns, reverse=False, terminate_on_match=False, postprocess=<class 'str'>)
General pattern reading. Uses monty’s regrep method. Takes the same
arguments.


* **Parameters**


    * **patterns** (*dict*) – A dict of patterns, e.g.,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”}.


    * **reverse** (*bool*) – Read files in reverse. Defaults to false. Useful for
    large files, esp OUTCARs, especially when used with
    terminate_on_match.


    * **terminate_on_match** (*bool*) – Whether to terminate when there is at
    least one match in each key in pattern.


    * **postprocess** (*callable*) – A post processing function to convert all
    matches. Defaults to str, i.e., no change.


Renders accessible:

    Any attribute in patterns. For example,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”} will set the
    value of self.data[“energy”] = [[-1234], [-3453], …], to the
    results from regex and postprocess. Note that the returned
    values are lists of lists, because you can grep multiple
    items on one line.

## pymatgen.io.res module

Provides parsing and read/write support for ShelX .res files as produced by the AIRSS code.

Converting from and back to pymatgen objects is expected to be reversible, i.e. you
should get the same Structure or ComputedStructureEntry back. On the other hand, converting
from and back to a string/file is not guaranteed to be reversible, i.e. a diff on the output
would not be empty. The difference should be limited to whitespace, float precision, and the
REM entries.


### _class_ AirssProvider(res: Res, parse_rems: Literal['gentle', 'strict'] = 'gentle')
Bases: `ResProvider`

Provides access to the res file as does ResProvider. This class additionally provides
access to fields in the TITL entry and various other fields found in the REM entries
that AIRSS puts in the file. Values in the TITL entry that AIRSS could not get end up as 0.
If the TITL entry is malformed, empty, or missing then attempting to construct this class
from a res file will raise a ResError.

While AIRSS supports a number of geometry and energy solvers, CASTEP is the default. As such,
fetching the information from the REM entries is only supported if AIRSS was used with CASTEP.
The other properties that get put in the TITL should still be accessible even if CASTEP was
not used.

The `parse_rems` attribute controls whether functions that fail to retrieve information
from the REM entries should return `None`. If this is set to `"strict"`,
then a ParseError may be raised, but the return value will not be `None`.
If it is set to `"gentle"`, then `None` will be returned instead of raising an
exception. This setting applies to all methods of this class that are typed to return
an Optional type. Default is `"gentle"`.

The `from_str()` and `from_file()` methods should be used instead of constructing this directly.


#### _date_fmt(_ = re.compile('[MTWFS][a-z]{2}, (\\\\d{2}) ([A-Z][a-z]{2}) (\\\\d{4}) (\\\\d{2}):(\\\\d{2}):(\\\\d{2}) ([+-]?\\\\d{4})'_ )

#### _get_compile_options()

#### _get_compiler()

#### _get_rng_seeds()

#### _classmethod_ _parse_date(string: str)
Parses a date from a string where the date is in the format typically used by CASTEP.


#### _raise_or_none(err: ResParseError)

#### _property_ appearances(_: in_ )
This is sometimes the number of times a structure was found in an AIRSS search.
Using the cryan tool that comes with AIRSS may be a better approach than relying
on this property.


#### as_dict(verbose: bool = True)
Get dict with title fields, structure and rems of this AirssProvider.


#### _property_ energy(_: floa_ )
Energy of the structure. With CASTEP, this is usually the enthalpy and is in eV.


#### _property_ entry(_: [ComputedStructureEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry_ )
Get this res file as a ComputedStructureEntry.


#### _classmethod_ from_file(filename: str, parse_rems: Literal['gentle', 'strict'] = 'gentle')
Construct a Provider from a file.


#### _classmethod_ from_str(string: str, parse_rems: Literal['gentle', 'strict'] = 'gentle')
Construct a Provider from a string.


#### get_airss_version()
Retrieves the version of AIRSS that was used along with the build date (not compile date).


* **Returns**

    (version string, date)



#### get_castep_version()
Retrieves the version of CASTEP that the res file was computed with from the REM entries.


* **Returns**

    version string



#### get_cut_grid_gmax_fsbc()
Retrieves the cut-off energy, grid scale, Gmax, and finite basis set correction setting
from the REM entries.


* **Returns**

    (cut-off, grid scale, Gmax, fsbc)



#### get_func_rel_disp()
Retrieves the functional, relativity scheme, and dispersion correction from the REM entries.


* **Returns**

    (functional, relativity, dispersion)



#### get_mpgrid_offset_nkpts_spacing()
Retrieves the MP grid, the grid offsets, number of kpoints, and maximum kpoint spacing.


* **Returns**

    (MP grid), (offsets), No. kpts, max spacing)



#### get_pspots()
Retrieves the OTFG pseudopotential string that can be used to generate the
pseudopotentials used in the calculation.


* **Returns**

    dict[specie, potential]



#### get_run_start_info()
Retrieves the run start date and the path it was started in from the REM entries.


* **Returns**

    (date, path)



#### _property_ integrated_absolute_spin_density(_: floa_ )
Corresponds to the last `Integrated |Spin Density|` in the CASTEP file.


#### _property_ integrated_spin_density(_: floa_ )
Corresponds to the last `Integrated Spin Density` in the CASTEP file.


#### _property_ pressure(_: floa_ )
Pressure for the run. This is in GPa if CASTEP was used.


#### _property_ seed(_: st_ )
The seed name, typically also the name of the res file.


#### _property_ spacegroup_label(_: st_ )
The Hermann-Mauguin notation of the spacegroup with ascii characters.
So no. 225 would be Fm-3m, and no. 194 would be P6_3/mmc.


#### _property_ volume(_: floa_ )
Volume of the structure. This is in cubic Angstroms if CASTEP was used.


### _class_ AirssTITL(seed: 'str', pressure: 'float', volume: 'float', energy: 'float', integrated_spin_density: 'float', integrated_absolute_spin_density: 'float', spacegroup_label: 'str', appearances: 'int')
Bases: `object`


#### appearances(_: in_ )

#### energy(_: floa_ )

#### integrated_absolute_spin_density(_: floa_ )

#### integrated_spin_density(_: floa_ )

#### pressure(_: floa_ )

#### seed(_: st_ )

#### spacegroup_label(_: st_ )

#### volume(_: floa_ )

### _class_ Ion(specie: 'str', specie_num: 'int', pos: 'Vector3D', occupancy: 'float', spin: 'float | None')
Bases: `object`


#### occupancy(_: floa_ )

#### pos(_: Vector3_ )

#### specie(_: st_ )

#### specie_num(_: in_ )

#### spin(_: float | Non_ )

### _class_ Res(TITL: AirssTITL | None, REMS: list[str], CELL: ResCELL, SFAC: ResSFAC)
Bases: `object`

Representation for the data in a res file.


#### CELL(_: ResCEL_ )

#### REMS(_: list[str_ )

#### SFAC(_: ResSFA_ )

#### TITL(_: AirssTITL | Non_ )

### _class_ ResCELL(unknown_field_1: 'float', a: 'float', b: 'float', c: 'float', alpha: 'float', beta: 'float', gamma: 'float')
Bases: `object`


#### a(_: floa_ )

#### alpha(_: floa_ )

#### b(_: floa_ )

#### beta(_: floa_ )

#### c(_: floa_ )

#### gamma(_: floa_ )

#### unknown_field_1(_: floa_ )

### _exception_ ResError()
Bases: `ValueError`

This exception indicates a problem was encountered while trying to retrieve a value or
perform an action that a provider for the res file does not support.


### _class_ ResIO()
Bases: `object`

Class providing convenience methods for converting a Structure or ComputedStructureEntry
to/from a string or file in the res format as used by AIRSS.

Note: Converting from and back to pymatgen objects is expected to be reversible, i.e. you
should get the same Structure or ComputedStructureEntry back. On the other hand, converting
from and back to a string/file is not guaranteed to be reversible, i.e. a diff on the output
would not be empty. The difference should be limited to whitespace, float precision, and the
REM entries.

If the TITL entry doesn’t exist or is malformed or empty, then you can only get
a Structure. Attempting to get an Entry will raise a ResError.


#### _classmethod_ entry_from_file(filename: str)
Produce a pymatgen ComputedStructureEntry from a res file.


#### _classmethod_ entry_from_str(string: str)
Produce a pymatgen ComputedStructureEntry from contents of a res file.


#### _classmethod_ entry_to_file(entry: [ComputedStructureEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry), filename: str)
Write a pymatgen ComputedStructureEntry to a res file.


#### _classmethod_ entry_to_str(entry: [ComputedStructureEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))
Produce the contents of a res file from a pymatgen ComputedStructureEntry.


#### _classmethod_ structure_from_file(filename: str)
Produces a pymatgen Structure from a res file.


#### _classmethod_ structure_from_str(string: str)
Produces a pymatgen Structure from contents of a res file.


#### _classmethod_ structure_to_file(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), filename: str)
Write a pymatgen Structure to a res file.


#### _classmethod_ structure_to_str(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Produce the contents of a res file from a pymatgen Structure.


### _exception_ ResParseError()
Bases: `ParseError`

This exception indicates a problem was encountered during parsing due to unexpected formatting.


### _class_ ResParser()
Bases: `object`

Parser for the ShelX res file.


#### _parse_cell(line: str)
Parses the CELL entry.


#### _classmethod_ _parse_file(filename: str)
Parses the res file as a file.


#### _parse_ion(line: str)
Parses entries in the SFAC block.


#### _parse_sfac(line: str, it: Iterator[str])
Parses the SFAC block.


#### _classmethod_ _parse_str(source: str)
Parses the res file as a string.


#### _parse_titl(line: str)
Parses the TITL entry. Checks for AIRSS values in the entry.


#### _parse_txt()
Parses the text of the file.


### _class_ ResProvider(res: Res)
Bases: `MSONable`

Provides access to elements of the res file in the form of familiar pymatgen objects.

The `from_str()` and `from_file()` methods should be used instead of constructing this directly.


#### _classmethod_ _site_spin(spin: float | None)
Check and return a dict with the site spin. Return None if spin is None.


#### _classmethod_ from_file(filename: str)
Construct a Provider from a file.


#### _classmethod_ from_str(string: str)
Construct a Provider from a string.


#### _property_ lattice(_: [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice_ )
Construct a Lattice from the res file.


#### _property_ rems(_: list[str_ )
The full list of REM entries contained within the res file.


#### _property_ sites(_: list[[PeriodicSite](pymatgen.core.md#pymatgen.core.sites.PeriodicSite)_ )
Construct a list of PeriodicSites from the res file.


#### _property_ structure(_: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure_ )
Construct a Structure from the res file.


### _class_ ResSFAC(species: 'set[str]', ions: 'list[Ion]')
Bases: `object`


#### ions(_: list[Ion_ )

#### species(_: set[str_ )

### _class_ ResWriter(entry: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | [ComputedStructureEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))
Bases: `object`

This class provides a means to write a Structure or ComputedStructureEntry to a res file.

This class can be constructed from either a pymatgen Structure or ComputedStructureEntry object.


#### _classmethod_ _cell_from_lattice(lattice: [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice))
Produce CELL entry from a pymatgen Lattice.


#### _classmethod_ _ions_from_sites(sites: list[[PeriodicSite](pymatgen.core.md#pymatgen.core.sites.PeriodicSite)])
Produce a list of entries for a SFAC block from a list of pymatgen PeriodicSite.


#### _classmethod_ _res_from_entry(entry: [ComputedStructureEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))
Produce a res file structure from a pymatgen ComputedStructureEntry.


#### _classmethod_ _res_from_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Produce a res file structure from a pymatgen Structure.


#### _classmethod_ _sfac_from_sites(sites: list[[PeriodicSite](pymatgen.core.md#pymatgen.core.sites.PeriodicSite)])
Produce a SFAC block from a list of pymatgen PeriodicSite.


#### _property_ string(_: st_ )
The contents of the res file.


#### write(filename: str)
Write the res data to a file.

## pymatgen.io.shengbte module

This module implements reading and writing of ShengBTE CONTROL files.


### _class_ Control(ngrid: list[int] | None = None, temperature: float | dict[str, float] = 300, \*\*kwargs)
Bases: `MSONable`, `dict`

Class for reading, updating, and writing ShengBTE CONTROL files.
See  [https://bitbucket.org/sousaw/shengbte/src/master/](https://bitbucket.org/sousaw/shengbte/src/master/) for more
detailed description and default values of CONTROL arguments.


* **Parameters**


    * **ngrid** – Reciprocal space grid density as a list of 3 ints.


    * **temperature** – The temperature to calculate the lattice thermal
    conductivity for. Can be given as a single float, or a dictionary
    with the keys “min”, “max”, “step”.


    * **\*\*kwargs** – Other ShengBTE parameters. Several parameters are required
    for ShengBTE to run - we have listed these parameters below:
    - nelements (int): number of different elements in the compound
    - natoms (int): number of atoms in the unit cell
    - lattvec (size 3x3 array): real-space lattice vectors, in units

    > of lfactor


        * lfactor (float): unit of measurement for lattice vectors (nm).

        I.e., set to 0.1 if lattvec given in Angstrom.


        * types (size natom list): a vector of natom integers, ranging
    from 1 to nelements, assigning an element to each atom in the
    system


        * elements (size natom list): a vector of element names


        * positions (size natomx3 array): atomic positions in lattice
    coordinates


        * scell (size 3 list): supercell sizes along each crystal axis
    used for the 2nd-order force constant calculation




#### allocations_keys(_ = ('nelements', 'natoms', 'ngrid', 'norientations'_ )

#### as_dict()
Returns: MSONable dict.


#### crystal_keys(_ = ('lfactor', 'lattvec', 'types', 'elements', 'positions', 'masses', 'gfactors', 'epsilon', 'born', 'scell', 'orientations'_ )

#### data_keys(_ = ('nelements', 'natoms', 'ngrid', 'lattvec', 'types', 'elements', 'positions', 'scell'_ )

#### flags_keys(_ = ('nonanalytic', 'convergence', 'isotopes', 'autoisotopes', 'nanowires', 'onlyharmonic', 'espresso'_ )

#### _classmethod_ from_dict(control_dict: dict)
Write a CONTROL file from a Python dictionary. Description and default
parameters can be found at
[https://bitbucket.org/sousaw/shengbte/src/master/](https://bitbucket.org/sousaw/shengbte/src/master/).
Note some parameters are mandatory. Optional parameters default here to
None and will not be written to file.


* **Parameters**

    **control_dict** – A Python dictionary of ShengBTE input parameters.



#### _classmethod_ from_file(filepath: str)
Read a CONTROL namelist file and output a ‘Control’ object.


* **Parameters**

    **filepath** – Path of the CONTROL file.



* **Returns**

    ‘Control’ object with parameters instantiated.



#### _classmethod_ from_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), reciprocal_density: int | None = 50000, \*\*kwargs)
Get a ShengBTE control object from a structure.


* **Parameters**


    * **structure** – A structure object.


    * **reciprocal_density** – If not None, the q-point grid (“ngrid”) will be
    set using this density.


    * **kwargs** – Additional options to be passed to the Control constructor.
    See the docstring of the __init__ method for more details



* **Returns**

    A ShengBTE control object.



#### get_structure()
Get a pymatgen Structure from a ShengBTE control object.

The control object must have the “lattvec”, “types”, “elements”, and
“positions” settings otherwise an error will be thrown.


* **Returns**

    The structure.



#### params_keys(_ = ('t', 't_min', 't_max', 't_step', 'omega_max', 'scalebroad', 'rmin', 'rmax', 'dr', 'maxiter', 'nticks', 'eps'_ )

#### required_params(_ = ('nelements', 'natoms', 'ngrid', 'lattvec', 'types', 'elements', 'positions', 'scell'_ )

#### to_file(filename: str = 'CONTROL')
Writes ShengBTE CONTROL file from ‘Control’ object.


* **Parameters**

    **filename** – A file name.



### _get_subdict(master_dict, subkeys)
Helper method to get a set of keys from a larger dictionary.

## pymatgen.io.template module

This module defines a simple concrete implementation of the InputGenerator class that can be
used to facilitate writing large numbers of input files based on a template.


### _class_ TemplateInputGen()
Bases: `InputGenerator`

Concrete implementation of InputGenerator that is based on a single template input
file with variables.

This class is provided as a low-barrier way to support new codes and to provide
an intuitive way for users to transition from manual scripts to pymatgen I/O
classes.


#### get_input_set(template: str | Path, variables: dict | None = None, filename: str = 'input.txt')

* **Parameters**


    * **template** – the input file template containing variable strings to be
    replaced.


    * **variables** – dict of variables to replace in the template. Keys are the
    text to replaced with the values, e.g. {“TEMPERATURE”: 298} will
    replace the text $TEMPERATURE in the template. See Python’s
    Template.safe_substitute() method documentation for more details.


    * **filename** – name of the file to be written.


## pymatgen.io.wannier90 module

Modules for working with wannier90 input and output.


### _class_ Unk(ik: int, data: ndarray)
Bases: `object`

Object representing the data in a UNK file.


#### ik()
Index of kpoint for this file.


* **Type**

    int



#### data()
Numpy array that contains the wavefunction data for in the UNK file.
The shape should be (nbnd, ngx, ngy, ngz) for regular calculations and (nbnd, 2, ngx, ngy, ngz)
for noncollinear calculations.


* **Type**

    numpy.ndarray



#### is_noncollinear()
Boolean that specifies if data is from a noncollinear calculation.


* **Type**

    bool



#### nbnd()
Number of bands in data.


* **Type**

    int



#### ng()
Sequence of three integers that correspond to the grid size of the given data.
The definition is ng = (ngx, ngy, ngz).


* **Type**

    tuple


Initialize Unk class.


* **Parameters**


    * **ik** (*int*) – index of the kpoint UNK file is for


    * **data** (*np.ndarray*) – data from the UNK file that has shape (nbnd,
    ngx, ngy, ngz) or (nbnd, 2, ngx, ngy, ngz) if noncollinear



#### _property_ data(_: ndarra_ )
contains the wavefunction data for in the UNK file.
The shape should be (nbnd, ngx, ngy, ngz) for regular calculations and
(nbnd, 2, ngx, ngy, ngz) for noncollinear calculations.


* **Type**

    np.ndarray



#### _static_ from_file(filename: str)
Reads the UNK data from file.


* **Parameters**

    **filename** (*str*) – path to UNK file to read



* **Returns**

    Unk object



#### ik(_: in_ )

#### is_noncollinear(_: boo_ )

#### nbnd(_: in_ )

#### ng(_: Sequence[int_ )

#### write_file(filename: str)
Write the UNK file.


* **Parameters**

    **filename** (*str*) – path to UNK file to write, the name should have the
    form ‘UNKXXXXX.YY’ where XXXXX is the kpoint index (Unk.ik) and
    YY is 1 or 2 for the spin index or NC if noncollinear


## pymatgen.io.xcrysden module

Support for reading XCrysDen files.


### _class_ XSF(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Bases: `object`

Class for parsing XCrysden files.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object.



#### _classmethod_ from_str(input_string, cls_=None)
Initialize a Structure object from a string with data in XSF format.


* **Parameters**


    * **input_string** – String with the structure in XSF format.
    See [http://www.xcrysden.org/doc/XSF.html](http://www.xcrysden.org/doc/XSF.html)


    * **cls** – Structure class to be created. default: pymatgen structure



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### to_str(atom_symbol=True)
Returns a string with the structure in XSF format
See [http://www.xcrysden.org/doc/XSF.html](http://www.xcrysden.org/doc/XSF.html).


* **Parameters**

    **atom_symbol** (*bool*) – Uses atom symbol instead of atomic number. Defaults to True.



#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead

## pymatgen.io.xr module

This module provides input and output mechanisms
for the xr file format, which is a modified CSSR
file format and, for example, used in GULP.
In particular, the module makes it easy
to remove shell positions from relaxations
that employed core-shell models.


### _class_ Xr(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Bases: `object`

Basic object for working with xr files.


* **Parameters**

    **structure** (*Structure/IStructure*) – Structure object to create the
    Xr object.



#### _static_ from_file(filename, use_cores=True, thresh=0.0001)
Reads an xr-formatted file to create an Xr object.


* **Parameters**


    * **filename** (*str*) – name of file to read from.


    * **use_cores** (*bool*) – use core positions and discard shell
    positions if set to True (default). Otherwise,
    use shell positions and discard core positions.


    * **thresh** (*float*) – relative threshold for consistency check
    between cell parameters (lengths and angles) from
    header information and cell vectors, respectively.



* **Returns**

    Xr object corresponding to the input

        file.




* **Return type**

    xr (Xr)



#### _static_ from_str(string, use_cores=True, thresh=0.0001)
Creates an Xr object from a string representation.


* **Parameters**


    * **string** (*str*) – string representation of an Xr object.


    * **use_cores** (*bool*) – use core positions and discard shell
    positions if set to True (default). Otherwise,
    use shell positions and discard core positions.


    * **thresh** (*float*) – relative threshold for consistency check
    between cell parameters (lengths and angles) from
    header information and cell vectors, respectively.



* **Returns**

    Xr object corresponding to the input

        string representation.




* **Return type**

    xr (Xr)



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### write_file(filename)
Write out an xr file.


* **Parameters**

    **filename** (*str*) – name of the file to write to.


## pymatgen.io.xyz module

Module implementing an XYZ file object class.


### _class_ XYZ(mol: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule) | [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | Sequence[[Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule) | [Structure](pymatgen.core.md#pymatgen.core.structure.Structure)], coord_precision: int = 6)
Bases: `object`

Basic class for importing and exporting Molecules or Structures in XYZ
format.

**NOTE**: Exporting periodic structures in the XYZ format will lose information
about the periodicity. Essentially, only Cartesian coordinates are
written in this format and no information is retained about the
lattice.


* **Parameters**


    * **mol** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)* | *[*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input molecule or structure or list thereof.


    * **coord_precision** – Precision to be used for coordinates.



#### _frame_str(frame_mol)

#### _static_ _from_frame_string(contents)
Convert a single frame XYZ string to a molecule.


#### _property_ all_molecules(_: list[[Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule)_ )
Returns all the frames of molecule associated with this XYZ.


#### as_dataframe()
Generates a coordinates data frame with columns: atom, x, y, and z
In case of multiple frame XYZ, returns the last frame.


* **Returns**

    pandas.DataFrame



#### _static_ from_file(filename)
Creates XYZ object from a file.


* **Parameters**

    **filename** – XYZ filename



* **Returns**

    XYZ object



#### _static_ from_str(contents)
Creates XYZ object from a string.


* **Parameters**

    **contents** – String representing an XYZ file.



* **Returns**

    XYZ object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### _property_ molecule(_: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule_ )
Returns molecule associated with this XYZ. In case of multi-frame
XYZ, returns the last frame.


#### write_file(filename: str)
Writes XYZ to file.


* **Parameters**

    **filename** (*str*) – File name of output file.


## pymatgen.io.zeopp module

Module implementing classes and functions to use Zeo++
by Maciej Haranczyk.

If using this module, cite the following paper on Zeo++:
T.F. Willems, C.H. Rycroft, M. Kazi, J.C. Meza, and M. Haranczyk,
Algorithms and tools for high-throughput geometry-based analysis of crystalline porous materials,
Microporous and Mesoporous Materials, 149 (2012) 134-141.

### Zeo++ Installation Steps:

A stable version of Zeo++ can be obtained from [http://zeoplusplus.org](http://zeoplusplus.org).
Instructions can be found at [http://www.zeoplusplus.org/download.html](http://www.zeoplusplus.org/download.html)

### Zeo++ Post-Installation Checking:


1. Go to pymatgen/io/tests and run “python test_zeoio.py”
If Zeo++ python bindings are properly installed, the tests should
pass. One or two tests will be skipped.


1. Go to pymatgen/analysis/defects/tests and run
“python test_point_defects.py”. Lots of tests will be skipped if GULP
is not installed. But there should be no errors.


### _class_ ZeoCssr(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Bases: `Cssr`

ZeoCssr adds extra fields to CSSR sites to conform with Zeo++
input CSSR format. The coordinate system is rotated from xyz to zyx.
This change aligns the pivot axis of pymatgen (z-axis) to pivot axis
of Zeo++ (x-axis) for structural modifications.


* **Parameters**

    **structure** – A structure to create ZeoCssr object.



#### _static_ from_file(filename)
Reads a CSSR file to a ZeoCssr object.


* **Parameters**

    **filename** – Filename to read from.



* **Returns**

    ZeoCssr object.



#### _static_ from_str(string)
Reads a string representation to a ZeoCssr object.


* **Parameters**

    **string** – A string representation of a ZeoCSSR.



* **Returns**

    ZeoCssr object.



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


### _class_ ZeoVoronoiXYZ(mol)
Bases: `XYZ`

Class to read Voronoi Nodes from XYZ file written by Zeo++.
The sites have an additional column representing the voronoi node radius.
The voronoi node radius is represented by the site property voronoi_radius.


* **Parameters**

    **mol** – Input molecule holding the voronoi node information.



#### _static_ from_file(filename)
Creates XYZ object from a file.


* **Parameters**

    **filename** – XYZ filename



* **Returns**

    XYZ object



#### _static_ from_str(contents)
Creates Zeo++ Voronoi XYZ object from a string.
from_string method of XYZ class is being redefined.


* **Parameters**

    **contents** – String representing Zeo++ Voronoi XYZ file.



* **Returns**

    ZeoVoronoiXYZ object



### get_free_sphere_params(structure, rad_dict=None, probe_rad=0.1)
Analyze the void space in the input structure using voronoi decomposition
Calls Zeo++ for Voronoi decomposition.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **rad_dict** (*optional*) – Dictionary of radii of elements in structure.
    If not given, Zeo++ default values are used.
    Note: Zeo++ uses atomic radii of elements.
    For ionic structures, pass rad_dict with ionic radii


    * **probe_rad** (*optional*) – Sampling probe radius in Angstroms. Default is
    0.1 A



* **Returns**

    voronoi nodes as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure
    voronoi face centers as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure



### get_high_accuracy_voronoi_nodes(structure, rad_dict, probe_rad=0.1)
Analyze the void space in the input structure using high accuracy
voronoi decomposition.
Calls Zeo++ for Voronoi decomposition.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **rad_dict** (*optional*) – Dictionary of radii of elements in structure.
    If not given, Zeo++ default values are used.
    Note: Zeo++ uses atomic radii of elements.
    For ionic structures, pass rad_dict with ionic radii


    * **probe_rad** (*optional*) – Sampling probe radius in Angstroms.
    Default is 0.1 A



* **Returns**

    voronoi nodes as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure
    voronoi face centers as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure



### get_voronoi_nodes(structure, rad_dict=None, probe_rad=0.1)
Analyze the void space in the input structure using voronoi decomposition
Calls Zeo++ for Voronoi decomposition.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **rad_dict** (*optional*) – Dictionary of radii of elements in structure.
    If not given, Zeo++ default values are used.
    Note: Zeo++ uses atomic radii of elements.
    For ionic structures, pass rad_dict with ionic radii


    * **probe_rad** (*optional*) – Sampling probe radius in Angstroms. Default is
    0.1 A



* **Returns**

    voronoi nodes as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure
    voronoi face centers as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure