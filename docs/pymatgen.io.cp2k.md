---
layout: default
title: pymatgen.io.cp2k.md
nav_exclude: true
---

# pymatgen.io.cp2k package

Module for CP2K input/output parsing as well as sets for standard calculations.



* [pymatgen.io.cp2k.inputs module](pymatgen.io.cp2k.inputs.md)


    * [`AtomicMetadata`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata)


        * [`AtomicMetadata.info`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.info)


        * [`AtomicMetadata.element`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.element)


        * [`AtomicMetadata.potential`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.potential)


        * [`AtomicMetadata.name`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.name)


        * [`AtomicMetadata.alias_names`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.alias_names)


        * [`AtomicMetadata.filename`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.filename)


        * [`AtomicMetadata.version`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.version)


        * [`AtomicMetadata.alias_names`](pymatgen.io.cp2k.inputs.md#id0)


        * [`AtomicMetadata.element`](pymatgen.io.cp2k.inputs.md#id1)


        * [`AtomicMetadata.filename`](pymatgen.io.cp2k.inputs.md#id2)


        * [`AtomicMetadata.get_hash()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.get_hash)


        * [`AtomicMetadata.get_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.get_string)


        * [`AtomicMetadata.info`](pymatgen.io.cp2k.inputs.md#id3)


        * [`AtomicMetadata.name`](pymatgen.io.cp2k.inputs.md#id4)


        * [`AtomicMetadata.potential`](pymatgen.io.cp2k.inputs.md#id5)


        * [`AtomicMetadata.softmatch()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.AtomicMetadata.softmatch)


        * [`AtomicMetadata.version`](pymatgen.io.cp2k.inputs.md#id6)


    * [`Band_Structure`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Band_Structure)


        * [`Band_Structure.from_kpoints()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Band_Structure.from_kpoints)


    * [`BasisFile`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisFile)


        * [`BasisFile.from_str()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisFile.from_str)


    * [`BasisInfo`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo)


        * [`BasisInfo.electrons`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.electrons)


        * [`BasisInfo.core`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.core)


        * [`BasisInfo.valence`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.valence)


        * [`BasisInfo.polarization`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.polarization)


        * [`BasisInfo.diffuse`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.diffuse)


        * [`BasisInfo.cc`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.cc)


        * [`BasisInfo.pc`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.pc)


        * [`BasisInfo.sr`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.sr)


        * [`BasisInfo.molopt`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.molopt)


        * [`BasisInfo.admm`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.admm)


        * [`BasisInfo.lri`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.lri)


        * [`BasisInfo.contracted`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.contracted)


        * [`BasisInfo.xc`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.xc)


        * [`BasisInfo.admm`](pymatgen.io.cp2k.inputs.md#id7)


        * [`BasisInfo.cc`](pymatgen.io.cp2k.inputs.md#id8)


        * [`BasisInfo.contracted`](pymatgen.io.cp2k.inputs.md#id9)


        * [`BasisInfo.core`](pymatgen.io.cp2k.inputs.md#id10)


        * [`BasisInfo.diffuse`](pymatgen.io.cp2k.inputs.md#id11)


        * [`BasisInfo.electrons`](pymatgen.io.cp2k.inputs.md#id12)


        * [`BasisInfo.from_str()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.from_str)


        * [`BasisInfo.from_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.from_string)


        * [`BasisInfo.lri`](pymatgen.io.cp2k.inputs.md#id13)


        * [`BasisInfo.molopt`](pymatgen.io.cp2k.inputs.md#id14)


        * [`BasisInfo.pc`](pymatgen.io.cp2k.inputs.md#id15)


        * [`BasisInfo.polarization`](pymatgen.io.cp2k.inputs.md#id16)


        * [`BasisInfo.softmatch()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BasisInfo.softmatch)


        * [`BasisInfo.sr`](pymatgen.io.cp2k.inputs.md#id17)


        * [`BasisInfo.valence`](pymatgen.io.cp2k.inputs.md#id18)


        * [`BasisInfo.xc`](pymatgen.io.cp2k.inputs.md#id19)


    * [`BrokenSymmetry`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BrokenSymmetry)


        * [`BrokenSymmetry.from_el()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.BrokenSymmetry.from_el)


    * [`Cell`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Cell)


    * [`Coord`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Coord)


    * [`Cp2kInput`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Cp2kInput)


        * [`Cp2kInput.from_file()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Cp2kInput.from_file)


        * [`Cp2kInput.from_lines()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Cp2kInput.from_lines)


        * [`Cp2kInput.from_str()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Cp2kInput.from_str)


        * [`Cp2kInput.from_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Cp2kInput.from_string)


        * [`Cp2kInput.get_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Cp2kInput.get_string)


        * [`Cp2kInput.write_file()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Cp2kInput.write_file)


    * [`DOS`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.DOS)


    * [`DataFile`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.DataFile)


        * [`DataFile.from_file()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.DataFile.from_file)


        * [`DataFile.from_str()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.DataFile.from_str)


        * [`DataFile.from_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.DataFile.from_string)


        * [`DataFile.get_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.DataFile.get_string)


        * [`DataFile.objects`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.DataFile.objects)


        * [`DataFile.write_file()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.DataFile.write_file)


    * [`Davidson`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Davidson)


    * [`Dft`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Dft)


    * [`DftPlusU`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.DftPlusU)


    * [`Diagonalization`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Diagonalization)


    * [`E_Density_Cube`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.E_Density_Cube)


    * [`ForceEval`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.ForceEval)


    * [`GaussianTypeOrbitalBasisSet`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet)


        * [`GaussianTypeOrbitalBasisSet.info`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.info)


        * [`GaussianTypeOrbitalBasisSet.nset`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.nset)


        * [`GaussianTypeOrbitalBasisSet.n`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.n)


        * [`GaussianTypeOrbitalBasisSet.lmax`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.lmax)


        * [`GaussianTypeOrbitalBasisSet.lmin`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.lmin)


        * [`GaussianTypeOrbitalBasisSet.nshell`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.nshell)


        * [`GaussianTypeOrbitalBasisSet.exponents`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.exponents)


        * [`GaussianTypeOrbitalBasisSet.coefficients`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.coefficients)


        * [`GaussianTypeOrbitalBasisSet.coefficients`](pymatgen.io.cp2k.inputs.md#id20)


        * [`GaussianTypeOrbitalBasisSet.exponents`](pymatgen.io.cp2k.inputs.md#id21)


        * [`GaussianTypeOrbitalBasisSet.from_str()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.from_str)


        * [`GaussianTypeOrbitalBasisSet.from_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.from_string)


        * [`GaussianTypeOrbitalBasisSet.get_keyword()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.get_keyword)


        * [`GaussianTypeOrbitalBasisSet.get_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.get_string)


        * [`GaussianTypeOrbitalBasisSet.info`](pymatgen.io.cp2k.inputs.md#id22)


        * [`GaussianTypeOrbitalBasisSet.lmax`](pymatgen.io.cp2k.inputs.md#id23)


        * [`GaussianTypeOrbitalBasisSet.lmin`](pymatgen.io.cp2k.inputs.md#id24)


        * [`GaussianTypeOrbitalBasisSet.n`](pymatgen.io.cp2k.inputs.md#id25)


        * [`GaussianTypeOrbitalBasisSet.nexp`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet.nexp)


        * [`GaussianTypeOrbitalBasisSet.nset`](pymatgen.io.cp2k.inputs.md#id26)


        * [`GaussianTypeOrbitalBasisSet.nshell`](pymatgen.io.cp2k.inputs.md#id27)


    * [`Global`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Global)


    * [`GthPotential`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential)


        * [`GthPotential.info`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.info)


        * [`GthPotential.n_elecs`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.n_elecs)


        * [`GthPotential.r_loc`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.r_loc)


        * [`GthPotential.nexp_ppl`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.nexp_ppl)


        * [`GthPotential.c_exp_ppl`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.c_exp_ppl)


        * [`GthPotential.radii`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.radii)


        * [`GthPotential.nprj`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.nprj)


        * [`GthPotential.nprj_ppnl`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.nprj_ppnl)


        * [`GthPotential.hprj_ppnl`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.hprj_ppnl)


        * [`GthPotential.c_exp_ppl`](pymatgen.io.cp2k.inputs.md#id28)


        * [`GthPotential.from_section()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.from_section)


        * [`GthPotential.from_str()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.from_str)


        * [`GthPotential.from_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.from_string)


        * [`GthPotential.get_keyword()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.get_keyword)


        * [`GthPotential.get_section()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.get_section)


        * [`GthPotential.get_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.GthPotential.get_string)


        * [`GthPotential.hprj_ppnl`](pymatgen.io.cp2k.inputs.md#id29)


        * [`GthPotential.n_elecs`](pymatgen.io.cp2k.inputs.md#id30)


        * [`GthPotential.nexp_ppl`](pymatgen.io.cp2k.inputs.md#id31)


        * [`GthPotential.nprj`](pymatgen.io.cp2k.inputs.md#id32)


        * [`GthPotential.nprj_ppnl`](pymatgen.io.cp2k.inputs.md#id33)


        * [`GthPotential.r_loc`](pymatgen.io.cp2k.inputs.md#id34)


        * [`GthPotential.radii`](pymatgen.io.cp2k.inputs.md#id35)


    * [`Keyword`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Keyword)


        * [`Keyword.as_dict()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Keyword.as_dict)


        * [`Keyword.from_dict()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Keyword.from_dict)


        * [`Keyword.from_str()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Keyword.from_str)


        * [`Keyword.from_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Keyword.from_string)


        * [`Keyword.get_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Keyword.get_string)


        * [`Keyword.verbosity()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Keyword.verbosity)


    * [`KeywordList`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.KeywordList)


        * [`KeywordList.append()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.KeywordList.append)


        * [`KeywordList.extend()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.KeywordList.extend)


        * [`KeywordList.get_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.KeywordList.get_string)


        * [`KeywordList.verbosity()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.KeywordList.verbosity)


    * [`Kind`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kind)


    * [`Kpoint_Set`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoint_Set)


    * [`Kpoints`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)


        * [`Kpoints.from_kpoints()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints.from_kpoints)


    * [`LDOS`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.LDOS)


    * [`MO_Cubes`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.MO_Cubes)


    * [`Mgrid`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Mgrid)


    * [`OrbitalTransformation`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.OrbitalTransformation)


    * [`PBE`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PBE)


    * [`PDOS`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PDOS)


    * [`PotentialFile`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialFile)


        * [`PotentialFile.from_str()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialFile.from_str)


    * [`PotentialInfo`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialInfo)


        * [`PotentialInfo.electrons`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialInfo.electrons)


        * [`PotentialInfo.potential_type`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialInfo.potential_type)


        * [`PotentialInfo.nlcc`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialInfo.nlcc)


        * [`PotentialInfo.xc`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialInfo.xc)


        * [`PotentialInfo.electrons`](pymatgen.io.cp2k.inputs.md#id36)


        * [`PotentialInfo.from_str()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialInfo.from_str)


        * [`PotentialInfo.from_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialInfo.from_string)


        * [`PotentialInfo.nlcc`](pymatgen.io.cp2k.inputs.md#id37)


        * [`PotentialInfo.potential_type`](pymatgen.io.cp2k.inputs.md#id38)


        * [`PotentialInfo.softmatch()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.PotentialInfo.softmatch)


        * [`PotentialInfo.xc`](pymatgen.io.cp2k.inputs.md#id39)


    * [`QS`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.QS)


    * [`Scf`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Scf)


    * [`Section`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section)


        * [`Section.add()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.add)


        * [`Section.by_path()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.by_path)


        * [`Section.check()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.check)


        * [`Section.get()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.get)


        * [`Section.get_keyword()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.get_keyword)


        * [`Section.get_section()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.get_section)


        * [`Section.get_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.get_string)


        * [`Section.inc()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.inc)


        * [`Section.insert()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.insert)


        * [`Section.safeset()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.safeset)


        * [`Section.set()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.set)


        * [`Section.setitem()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.setitem)


        * [`Section.silence()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.silence)


        * [`Section.unset()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.unset)


        * [`Section.update()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.update)


        * [`Section.verbosity()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Section.verbosity)


    * [`SectionList`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.SectionList)


        * [`SectionList.append()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.SectionList.append)


        * [`SectionList.extend()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.SectionList.extend)


        * [`SectionList.get()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.SectionList.get)


        * [`SectionList.get_string()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.SectionList.get_string)


        * [`SectionList.verbosity()`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.SectionList.verbosity)


    * [`Smear`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Smear)


    * [`Subsys`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Subsys)


    * [`V_Hartree_Cube`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.V_Hartree_Cube)


    * [`Xc_Functional`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Xc_Functional)


* [pymatgen.io.cp2k.outputs module](pymatgen.io.cp2k.outputs.md)


    * [`Cp2kOutput`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput)


        * [`Cp2kOutput.as_dict()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.as_dict)


        * [`Cp2kOutput.band_structure`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.band_structure)


        * [`Cp2kOutput.calculation_type`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.calculation_type)


        * [`Cp2kOutput.charge`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.charge)


        * [`Cp2kOutput.complete_dos`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.complete_dos)


        * [`Cp2kOutput.completed`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.completed)


        * [`Cp2kOutput.convergence()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.convergence)


        * [`Cp2kOutput.cp2k_version`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.cp2k_version)


        * [`Cp2kOutput.is_hubbard`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.is_hubbard)


        * [`Cp2kOutput.is_metal`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.is_metal)


        * [`Cp2kOutput.is_molecule`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.is_molecule)


        * [`Cp2kOutput.multiplicity`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.multiplicity)


        * [`Cp2kOutput.num_warnings`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.num_warnings)


        * [`Cp2kOutput.parse_atomic_kind_info()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_atomic_kind_info)


        * [`Cp2kOutput.parse_bandstructure()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_bandstructure)


        * [`Cp2kOutput.parse_cell_params()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_cell_params)


        * [`Cp2kOutput.parse_chi_tensor()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_chi_tensor)


        * [`Cp2kOutput.parse_cp2k_params()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_cp2k_params)


        * [`Cp2kOutput.parse_dft_params()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_dft_params)


        * [`Cp2kOutput.parse_dos()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_dos)


        * [`Cp2kOutput.parse_energies()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_energies)


        * [`Cp2kOutput.parse_files()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_files)


        * [`Cp2kOutput.parse_forces()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_forces)


        * [`Cp2kOutput.parse_global_params()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_global_params)


        * [`Cp2kOutput.parse_gtensor()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_gtensor)


        * [`Cp2kOutput.parse_hirshfeld()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_hirshfeld)


        * [`Cp2kOutput.parse_homo_lumo()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_homo_lumo)


        * [`Cp2kOutput.parse_hyperfine()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_hyperfine)


        * [`Cp2kOutput.parse_initial_structure()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_initial_structure)


        * [`Cp2kOutput.parse_input()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_input)


        * [`Cp2kOutput.parse_ionic_steps()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_ionic_steps)


        * [`Cp2kOutput.parse_mo_eigenvalues()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_mo_eigenvalues)


        * [`Cp2kOutput.parse_mulliken()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_mulliken)


        * [`Cp2kOutput.parse_nmr_shift()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_nmr_shift)


        * [`Cp2kOutput.parse_opt_steps()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_opt_steps)


        * [`Cp2kOutput.parse_overlap_condition()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_overlap_condition)


        * [`Cp2kOutput.parse_plus_u_params()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_plus_u_params)


        * [`Cp2kOutput.parse_qs_params()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_qs_params)


        * [`Cp2kOutput.parse_raman()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_raman)


        * [`Cp2kOutput.parse_scf_opt()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_scf_opt)


        * [`Cp2kOutput.parse_scf_params()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_scf_params)


        * [`Cp2kOutput.parse_stresses()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_stresses)


        * [`Cp2kOutput.parse_structures()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_structures)


        * [`Cp2kOutput.parse_tddfpt()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_tddfpt)


        * [`Cp2kOutput.parse_timing()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_timing)


        * [`Cp2kOutput.parse_total_numbers()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.parse_total_numbers)


        * [`Cp2kOutput.project_name`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.project_name)


        * [`Cp2kOutput.ran_successfully()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.ran_successfully)


        * [`Cp2kOutput.read_pattern()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.read_pattern)


        * [`Cp2kOutput.read_table_pattern()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.read_table_pattern)


        * [`Cp2kOutput.run_type`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.run_type)


        * [`Cp2kOutput.spin_polarized`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.Cp2kOutput.spin_polarized)


    * [`parse_dos()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.parse_dos)


    * [`parse_energy_file()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.parse_energy_file)


    * [`parse_pdos()`](pymatgen.io.cp2k.outputs.md#pymatgen.io.cp2k.outputs.parse_pdos)


* [pymatgen.io.cp2k.sets module](pymatgen.io.cp2k.sets.md)


    * [`CellOptSet`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.CellOptSet)


    * [`Cp2kValidationError`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.Cp2kValidationError)


        * [`Cp2kValidationError.CP2K_VERSION`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.Cp2kValidationError.CP2K_VERSION)


    * [`DftSet`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet)


        * [`DftSet.activate_epr()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_epr)


        * [`DftSet.activate_fast_minimization()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_fast_minimization)


        * [`DftSet.activate_hybrid()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_hybrid)


        * [`DftSet.activate_hyperfine()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_hyperfine)


        * [`DftSet.activate_localize()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_localize)


        * [`DftSet.activate_motion()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_motion)


        * [`DftSet.activate_nmr()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_nmr)


        * [`DftSet.activate_nonperiodic()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_nonperiodic)


        * [`DftSet.activate_polar()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_polar)


        * [`DftSet.activate_robust_minimization()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_robust_minimization)


        * [`DftSet.activate_spinspin()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_spinspin)


        * [`DftSet.activate_tddfpt()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_tddfpt)


        * [`DftSet.activate_vdw_potential()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_vdw_potential)


        * [`DftSet.activate_very_strict_minimization()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.activate_very_strict_minimization)


        * [`DftSet.create_subsys()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.create_subsys)


        * [`DftSet.get_basis_and_potential()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.get_basis_and_potential)


        * [`DftSet.get_cutoff_from_basis()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.get_cutoff_from_basis)


        * [`DftSet.get_xc_functionals()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.get_xc_functionals)


        * [`DftSet.modify_dft_print_iters()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.modify_dft_print_iters)


        * [`DftSet.print_bandstructure()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_bandstructure)


        * [`DftSet.print_dos()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_dos)


        * [`DftSet.print_e_density()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_e_density)


        * [`DftSet.print_forces()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_forces)


        * [`DftSet.print_hirshfeld()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_hirshfeld)


        * [`DftSet.print_ldos()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_ldos)


        * [`DftSet.print_mo()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_mo)


        * [`DftSet.print_mo_cubes()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_mo_cubes)


        * [`DftSet.print_mulliken()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_mulliken)


        * [`DftSet.print_pdos()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_pdos)


        * [`DftSet.print_v_hartree()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.print_v_hartree)


        * [`DftSet.set_charge()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.set_charge)


        * [`DftSet.validate()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.validate)


        * [`DftSet.write_basis_set_file()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.write_basis_set_file)


        * [`DftSet.write_potential_file()`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.DftSet.write_potential_file)


    * [`HybridCellOptSet`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.HybridCellOptSet)


    * [`HybridRelaxSet`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.HybridRelaxSet)


    * [`HybridStaticSet`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.HybridStaticSet)


    * [`RelaxSet`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.RelaxSet)


    * [`StaticSet`](pymatgen.io.cp2k.sets.md#pymatgen.io.cp2k.sets.StaticSet)


* [pymatgen.io.cp2k.utils module](pymatgen.io.cp2k.utils.md)


    * [`chunk()`](pymatgen.io.cp2k.utils.md#pymatgen.io.cp2k.utils.chunk)


    * [`get_truncated_coulomb_cutoff()`](pymatgen.io.cp2k.utils.md#pymatgen.io.cp2k.utils.get_truncated_coulomb_cutoff)


    * [`get_unique_site_indices()`](pymatgen.io.cp2k.utils.md#pymatgen.io.cp2k.utils.get_unique_site_indices)


    * [`natural_keys()`](pymatgen.io.cp2k.utils.md#pymatgen.io.cp2k.utils.natural_keys)


    * [`postprocessor()`](pymatgen.io.cp2k.utils.md#pymatgen.io.cp2k.utils.postprocessor)


    * [`preprocessor()`](pymatgen.io.cp2k.utils.md#pymatgen.io.cp2k.utils.preprocessor)