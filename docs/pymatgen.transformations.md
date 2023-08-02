---
layout: default
title: pymatgen.transformations.md
nav_exclude: true
---

# pymatgen.transformations package

The transformations package defines various transformations that can be applied
on structures, i.e., converting one structure to another.



* [pymatgen.transformations.advanced_transformations module](pymatgen.transformations.advanced_transformations.md)


    * [`AddAdsorbateTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.AddAdsorbateTransformation)


        * [`AddAdsorbateTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.AddAdsorbateTransformation.apply_transformation)


        * [`AddAdsorbateTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.AddAdsorbateTransformation.inverse)


        * [`AddAdsorbateTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.AddAdsorbateTransformation.is_one_to_many)


    * [`ChargeBalanceTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.ChargeBalanceTransformation)


        * [`ChargeBalanceTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.ChargeBalanceTransformation.apply_transformation)


        * [`ChargeBalanceTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.ChargeBalanceTransformation.inverse)


        * [`ChargeBalanceTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.ChargeBalanceTransformation.is_one_to_many)


    * [`CubicSupercellTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.CubicSupercellTransformation)


        * [`CubicSupercellTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.CubicSupercellTransformation.apply_transformation)


        * [`CubicSupercellTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.CubicSupercellTransformation.inverse)


        * [`CubicSupercellTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.CubicSupercellTransformation.is_one_to_many)


    * [`DisorderOrderedTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.DisorderOrderedTransformation)


        * [`DisorderOrderedTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.DisorderOrderedTransformation.apply_transformation)


        * [`DisorderOrderedTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.DisorderOrderedTransformation.inverse)


        * [`DisorderOrderedTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.DisorderOrderedTransformation.is_one_to_many)


    * [`DopingTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.DopingTransformation)


        * [`DopingTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.DopingTransformation.apply_transformation)


        * [`DopingTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.DopingTransformation.inverse)


        * [`DopingTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.DopingTransformation.is_one_to_many)


    * [`EnumerateStructureTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.EnumerateStructureTransformation)


        * [`EnumerateStructureTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.EnumerateStructureTransformation.apply_transformation)


        * [`EnumerateStructureTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.EnumerateStructureTransformation.inverse)


        * [`EnumerateStructureTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.EnumerateStructureTransformation.is_one_to_many)


    * [`GrainBoundaryTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.GrainBoundaryTransformation)


        * [`GrainBoundaryTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.GrainBoundaryTransformation.apply_transformation)


        * [`GrainBoundaryTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.GrainBoundaryTransformation.inverse)


        * [`GrainBoundaryTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.GrainBoundaryTransformation.is_one_to_many)


    * [`MagOrderParameterConstraint`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MagOrderParameterConstraint)


        * [`MagOrderParameterConstraint.satisfies_constraint()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MagOrderParameterConstraint.satisfies_constraint)


    * [`MagOrderingTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MagOrderingTransformation)


        * [`MagOrderingTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MagOrderingTransformation.apply_transformation)


        * [`MagOrderingTransformation.determine_min_cell()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MagOrderingTransformation.determine_min_cell)


        * [`MagOrderingTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MagOrderingTransformation.inverse)


        * [`MagOrderingTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MagOrderingTransformation.is_one_to_many)


    * [`MonteCarloRattleTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MonteCarloRattleTransformation)


        * [`MonteCarloRattleTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MonteCarloRattleTransformation.apply_transformation)


        * [`MonteCarloRattleTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MonteCarloRattleTransformation.inverse)


        * [`MonteCarloRattleTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MonteCarloRattleTransformation.is_one_to_many)


    * [`MultipleSubstitutionTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MultipleSubstitutionTransformation)


        * [`MultipleSubstitutionTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MultipleSubstitutionTransformation.apply_transformation)


        * [`MultipleSubstitutionTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MultipleSubstitutionTransformation.inverse)


        * [`MultipleSubstitutionTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.MultipleSubstitutionTransformation.is_one_to_many)


    * [`SQSTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SQSTransformation)


        * [`SQSTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SQSTransformation.apply_transformation)


        * [`SQSTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SQSTransformation.inverse)


        * [`SQSTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SQSTransformation.is_one_to_many)


    * [`SlabTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SlabTransformation)


        * [`SlabTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SlabTransformation.apply_transformation)


        * [`SlabTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SlabTransformation.inverse)


        * [`SlabTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SlabTransformation.is_one_to_many)


    * [`SubstituteSurfaceSiteTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SubstituteSurfaceSiteTransformation)


        * [`SubstituteSurfaceSiteTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SubstituteSurfaceSiteTransformation.apply_transformation)


        * [`SubstituteSurfaceSiteTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SubstituteSurfaceSiteTransformation.inverse)


        * [`SubstituteSurfaceSiteTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SubstituteSurfaceSiteTransformation.is_one_to_many)


    * [`SubstitutionPredictorTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SubstitutionPredictorTransformation)


        * [`SubstitutionPredictorTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SubstitutionPredictorTransformation.apply_transformation)


        * [`SubstitutionPredictorTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SubstitutionPredictorTransformation.inverse)


        * [`SubstitutionPredictorTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SubstitutionPredictorTransformation.is_one_to_many)


    * [`SuperTransformation`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SuperTransformation)


        * [`SuperTransformation.apply_transformation()`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SuperTransformation.apply_transformation)


        * [`SuperTransformation.inverse`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SuperTransformation.inverse)


        * [`SuperTransformation.is_one_to_many`](pymatgen.transformations.advanced_transformations.md#pymatgen.transformations.advanced_transformations.SuperTransformation.is_one_to_many)


* [pymatgen.transformations.site_transformations module](pymatgen.transformations.site_transformations.md)


    * [`AddSitePropertyTransformation`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.AddSitePropertyTransformation)


        * [`AddSitePropertyTransformation.apply_transformation()`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.AddSitePropertyTransformation.apply_transformation)


        * [`AddSitePropertyTransformation.inverse`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.AddSitePropertyTransformation.inverse)


        * [`AddSitePropertyTransformation.is_one_to_many`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.AddSitePropertyTransformation.is_one_to_many)


    * [`InsertSitesTransformation`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.InsertSitesTransformation)


        * [`InsertSitesTransformation.apply_transformation()`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.InsertSitesTransformation.apply_transformation)


        * [`InsertSitesTransformation.inverse`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.InsertSitesTransformation.inverse)


        * [`InsertSitesTransformation.is_one_to_many`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.InsertSitesTransformation.is_one_to_many)


    * [`PartialRemoveSitesTransformation`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation)


        * [`PartialRemoveSitesTransformation.ALGO_BEST_FIRST`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation.ALGO_BEST_FIRST)


        * [`PartialRemoveSitesTransformation.ALGO_COMPLETE`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation.ALGO_COMPLETE)


        * [`PartialRemoveSitesTransformation.ALGO_ENUMERATE`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation.ALGO_ENUMERATE)


        * [`PartialRemoveSitesTransformation.ALGO_FAST`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation.ALGO_FAST)


        * [`PartialRemoveSitesTransformation.apply_transformation()`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation.apply_transformation)


        * [`PartialRemoveSitesTransformation.inverse`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation.inverse)


        * [`PartialRemoveSitesTransformation.is_one_to_many`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation.is_one_to_many)


    * [`RadialSiteDistortionTransformation`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.RadialSiteDistortionTransformation)


        * [`RadialSiteDistortionTransformation.apply_transformation()`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.RadialSiteDistortionTransformation.apply_transformation)


        * [`RadialSiteDistortionTransformation.inverse`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.RadialSiteDistortionTransformation.inverse)


        * [`RadialSiteDistortionTransformation.is_one_to_many`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.RadialSiteDistortionTransformation.is_one_to_many)


        * [`RadialSiteDistortionTransformation.use_multiprocessing`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.RadialSiteDistortionTransformation.use_multiprocessing)


    * [`RemoveSitesTransformation`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.RemoveSitesTransformation)


        * [`RemoveSitesTransformation.apply_transformation()`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.RemoveSitesTransformation.apply_transformation)


        * [`RemoveSitesTransformation.inverse`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.RemoveSitesTransformation.inverse)


        * [`RemoveSitesTransformation.is_one_to_many`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.RemoveSitesTransformation.is_one_to_many)


    * [`ReplaceSiteSpeciesTransformation`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.ReplaceSiteSpeciesTransformation)


        * [`ReplaceSiteSpeciesTransformation.apply_transformation()`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.ReplaceSiteSpeciesTransformation.apply_transformation)


        * [`ReplaceSiteSpeciesTransformation.inverse`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.ReplaceSiteSpeciesTransformation.inverse)


        * [`ReplaceSiteSpeciesTransformation.is_one_to_many`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.ReplaceSiteSpeciesTransformation.is_one_to_many)


    * [`TranslateSitesTransformation`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.TranslateSitesTransformation)


        * [`TranslateSitesTransformation.apply_transformation()`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.TranslateSitesTransformation.apply_transformation)


        * [`TranslateSitesTransformation.as_dict()`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.TranslateSitesTransformation.as_dict)


        * [`TranslateSitesTransformation.inverse`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.TranslateSitesTransformation.inverse)


        * [`TranslateSitesTransformation.is_one_to_many`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.TranslateSitesTransformation.is_one_to_many)


* [pymatgen.transformations.standard_transformations module](pymatgen.transformations.standard_transformations.md)


    * [`AutoOxiStateDecorationTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.AutoOxiStateDecorationTransformation)


        * [`AutoOxiStateDecorationTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.AutoOxiStateDecorationTransformation.apply_transformation)


        * [`AutoOxiStateDecorationTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.AutoOxiStateDecorationTransformation.inverse)


        * [`AutoOxiStateDecorationTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.AutoOxiStateDecorationTransformation.is_one_to_many)


    * [`ChargedCellTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ChargedCellTransformation)


        * [`ChargedCellTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ChargedCellTransformation.apply_transformation)


        * [`ChargedCellTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ChargedCellTransformation.inverse)


        * [`ChargedCellTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ChargedCellTransformation.is_one_to_many)


    * [`ConventionalCellTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ConventionalCellTransformation)


        * [`ConventionalCellTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ConventionalCellTransformation.apply_transformation)


        * [`ConventionalCellTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ConventionalCellTransformation.inverse)


        * [`ConventionalCellTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ConventionalCellTransformation.is_one_to_many)


    * [`DeformStructureTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.DeformStructureTransformation)


        * [`DeformStructureTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.DeformStructureTransformation.apply_transformation)


        * [`DeformStructureTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.DeformStructureTransformation.inverse)


        * [`DeformStructureTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.DeformStructureTransformation.is_one_to_many)


    * [`DiscretizeOccupanciesTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.DiscretizeOccupanciesTransformation)


        * [`DiscretizeOccupanciesTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.DiscretizeOccupanciesTransformation.apply_transformation)


        * [`DiscretizeOccupanciesTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.DiscretizeOccupanciesTransformation.inverse)


        * [`DiscretizeOccupanciesTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.DiscretizeOccupanciesTransformation.is_one_to_many)


    * [`OrderDisorderedStructureTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OrderDisorderedStructureTransformation)


        * [`OrderDisorderedStructureTransformation.ALGO_BEST_FIRST`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OrderDisorderedStructureTransformation.ALGO_BEST_FIRST)


        * [`OrderDisorderedStructureTransformation.ALGO_COMPLETE`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OrderDisorderedStructureTransformation.ALGO_COMPLETE)


        * [`OrderDisorderedStructureTransformation.ALGO_FAST`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OrderDisorderedStructureTransformation.ALGO_FAST)


        * [`OrderDisorderedStructureTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OrderDisorderedStructureTransformation.apply_transformation)


        * [`OrderDisorderedStructureTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OrderDisorderedStructureTransformation.inverse)


        * [`OrderDisorderedStructureTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OrderDisorderedStructureTransformation.is_one_to_many)


        * [`OrderDisorderedStructureTransformation.lowest_energy_structure`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OrderDisorderedStructureTransformation.lowest_energy_structure)


    * [`OxidationStateDecorationTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OxidationStateDecorationTransformation)


        * [`OxidationStateDecorationTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OxidationStateDecorationTransformation.apply_transformation)


        * [`OxidationStateDecorationTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OxidationStateDecorationTransformation.inverse)


        * [`OxidationStateDecorationTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OxidationStateDecorationTransformation.is_one_to_many)


    * [`OxidationStateRemovalTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OxidationStateRemovalTransformation)


        * [`OxidationStateRemovalTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OxidationStateRemovalTransformation.apply_transformation)


        * [`OxidationStateRemovalTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OxidationStateRemovalTransformation.inverse)


        * [`OxidationStateRemovalTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.OxidationStateRemovalTransformation.is_one_to_many)


    * [`PartialRemoveSpecieTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PartialRemoveSpecieTransformation)


        * [`PartialRemoveSpecieTransformation.ALGO_BEST_FIRST`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PartialRemoveSpecieTransformation.ALGO_BEST_FIRST)


        * [`PartialRemoveSpecieTransformation.ALGO_COMPLETE`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PartialRemoveSpecieTransformation.ALGO_COMPLETE)


        * [`PartialRemoveSpecieTransformation.ALGO_ENUMERATE`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PartialRemoveSpecieTransformation.ALGO_ENUMERATE)


        * [`PartialRemoveSpecieTransformation.ALGO_FAST`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PartialRemoveSpecieTransformation.ALGO_FAST)


        * [`PartialRemoveSpecieTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PartialRemoveSpecieTransformation.apply_transformation)


        * [`PartialRemoveSpecieTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PartialRemoveSpecieTransformation.inverse)


        * [`PartialRemoveSpecieTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PartialRemoveSpecieTransformation.is_one_to_many)


    * [`PerturbStructureTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PerturbStructureTransformation)


        * [`PerturbStructureTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PerturbStructureTransformation.apply_transformation)


        * [`PerturbStructureTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PerturbStructureTransformation.inverse)


        * [`PerturbStructureTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PerturbStructureTransformation.is_one_to_many)


    * [`PrimitiveCellTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PrimitiveCellTransformation)


        * [`PrimitiveCellTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PrimitiveCellTransformation.apply_transformation)


        * [`PrimitiveCellTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PrimitiveCellTransformation.inverse)


        * [`PrimitiveCellTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.PrimitiveCellTransformation.is_one_to_many)


    * [`RemoveSpeciesTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.RemoveSpeciesTransformation)


        * [`RemoveSpeciesTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.RemoveSpeciesTransformation.apply_transformation)


        * [`RemoveSpeciesTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.RemoveSpeciesTransformation.inverse)


        * [`RemoveSpeciesTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.RemoveSpeciesTransformation.is_one_to_many)


    * [`RotationTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.RotationTransformation)


        * [`RotationTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.RotationTransformation.apply_transformation)


        * [`RotationTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.RotationTransformation.inverse)


        * [`RotationTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.RotationTransformation.is_one_to_many)


    * [`ScaleToRelaxedTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ScaleToRelaxedTransformation)


        * [`ScaleToRelaxedTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ScaleToRelaxedTransformation.apply_transformation)


        * [`ScaleToRelaxedTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ScaleToRelaxedTransformation.inverse)


        * [`ScaleToRelaxedTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.ScaleToRelaxedTransformation.is_one_to_many)


    * [`SubstitutionTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.SubstitutionTransformation)


        * [`SubstitutionTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.SubstitutionTransformation.apply_transformation)


        * [`SubstitutionTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.SubstitutionTransformation.inverse)


        * [`SubstitutionTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.SubstitutionTransformation.is_one_to_many)


    * [`SupercellTransformation`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.SupercellTransformation)


        * [`SupercellTransformation.apply_transformation()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.SupercellTransformation.apply_transformation)


        * [`SupercellTransformation.from_scaling_factors()`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.SupercellTransformation.from_scaling_factors)


        * [`SupercellTransformation.inverse`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.SupercellTransformation.inverse)


        * [`SupercellTransformation.is_one_to_many`](pymatgen.transformations.standard_transformations.md#pymatgen.transformations.standard_transformations.SupercellTransformation.is_one_to_many)


* [pymatgen.transformations.transformation_abc module](pymatgen.transformations.transformation_abc.md)


    * [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)


        * [`AbstractTransformation.apply_transformation()`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation.apply_transformation)


        * [`AbstractTransformation.inverse`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation.inverse)


        * [`AbstractTransformation.is_one_to_many`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation.is_one_to_many)


        * [`AbstractTransformation.use_multiprocessing`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation.use_multiprocessing)