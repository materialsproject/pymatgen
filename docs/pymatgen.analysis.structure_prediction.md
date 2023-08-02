---
layout: default
title: pymatgen.analysis.structure_prediction.md
nav_exclude: true
---

# pymatgen.analysis.structure_prediction package

Utilities to predict new structures.



* [pymatgen.analysis.structure_prediction.dopant_predictor module](pymatgen.analysis.structure_prediction.dopant_predictor.md)


    * [`get_dopants_from_shannon_radii()`](pymatgen.analysis.structure_prediction.dopant_predictor.md#pymatgen.analysis.structure_prediction.dopant_predictor.get_dopants_from_shannon_radii)


    * [`get_dopants_from_substitution_probabilities()`](pymatgen.analysis.structure_prediction.dopant_predictor.md#pymatgen.analysis.structure_prediction.dopant_predictor.get_dopants_from_substitution_probabilities)


* [pymatgen.analysis.structure_prediction.substitution_probability module](pymatgen.analysis.structure_prediction.substitution_probability.md)


    * [`SubstitutionPredictor`](pymatgen.analysis.structure_prediction.substitution_probability.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor)


        * [`SubstitutionPredictor.composition_prediction()`](pymatgen.analysis.structure_prediction.substitution_probability.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor.composition_prediction)


        * [`SubstitutionPredictor.list_prediction()`](pymatgen.analysis.structure_prediction.substitution_probability.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor.list_prediction)


    * [`SubstitutionProbability`](pymatgen.analysis.structure_prediction.substitution_probability.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionProbability)


* [pymatgen.analysis.structure_prediction.substitutor module](pymatgen.analysis.structure_prediction.substitutor.md)


    * [`Substitutor`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor)


        * [`Substitutor.as_dict()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.as_dict)


        * [`Substitutor.from_dict()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.from_dict)


        * [`Substitutor.get_allowed_species()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.get_allowed_species)


        * [`Substitutor.pred_from_comp()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.pred_from_comp)


        * [`Substitutor.pred_from_list()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.pred_from_list)


        * [`Substitutor.pred_from_structures()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.pred_from_structures)


* [pymatgen.analysis.structure_prediction.volume_predictor module](pymatgen.analysis.structure_prediction.volume_predictor.md)


    * [`DLSVolumePredictor`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor)


        * [`DLSVolumePredictor.get_predicted_structure()`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor.get_predicted_structure)


        * [`DLSVolumePredictor.predict()`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor.predict)


    * [`RLSVolumePredictor`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor)


        * [`RLSVolumePredictor.get_predicted_structure()`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor.get_predicted_structure)


        * [`RLSVolumePredictor.predict()`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor.predict)