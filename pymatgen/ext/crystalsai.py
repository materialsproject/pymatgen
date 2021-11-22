# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module provides an interface to the MaterialsVirtualLab's MEGNet REST API
for rapid property prediction.
"""

import requests


class CrystalAIRester:
    """
    This is a high-level interface to the REST API for http://megnet.crystals.ai
    for property prediction. Using this API, you can use MatErials Graph
    Networks (MEGNet) developed by the Materials Virtual Lab to predict the
    properties of any given crystal. These models are trained on the latest
    versions of the Materials Project. The open-source code
    implementing these MEGNet models are available at
    https://github.com/materialsvirtuallab/megnet.

    For the details of MEGNet and benchmarks, please refer to the following work:
    Chen, C.; Ye, W.; Zuo, Y.; Zheng, C.; Ong, S. P. <i>Graph Networks as a Universal Machine Learning Framework
        for Molecules and Crystals.</i> Chemistry of Materials 2019, acs.chemmater.9b01294.
     DOI: <a href="https://doi.org/10.1021/acs.chemmater.9b01294">10.1021/acs.chemmater.9b01294</a>.</p>
    """

    def __init__(self):
        """
        Init for Rester.
        """
        self.session = requests.Session()
        self.url = "http://megnet.crystals.ai"

    def __enter__(self):
        """
        Support for "with" context.
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Support for "with" context.
        """
        self.session.close()

    def get_available_models(self):
        """
        Returns:
            Available model names. It should be noted that model names starting
            with log10 are for the log10 of that quantity.
        """
        response = self.session.get(self.url + "/models")
        return response.json()

    def predict_mp(self, model_name, mp_id):
        """
        Predict using the http://megnet.crystals.ai API.

        :param model_name: An available model in the
            http://megnet.crystals.ai API. Use get_available_models to find
            the model names.
        :param mp_id: A Materials Project id.
        :return: Predicted value. It should be noted that model names starting
            with log10 are for the log10 of that quantity and you should apply
            10 ** prediction to get the actual value.
        """
        response = None
        url = self.url + "/predict_mp/%s/%s" % (model_name, mp_id)
        try:
            response = self.session.get(url)
            if response.status_code in [200, 400]:
                return response.json()
            raise ValueError("REST query returned with error status code {}".format(response.status_code))
        except Exception as ex:
            msg = "{}. Content: {}".format(str(ex), response.content) if hasattr(response, "content") else str(ex)
            raise ValueError(msg)

    def predict_structure(self, model_name, structure):
        """
        Predict using the http://megnet.crystals.ai API.

        :param model_name: An available model in the
            http://megnet.crystals.ai API. Use get_available_models to find
            the model names.
        :param structure: A Pymatgen Structure.
        :return: Predicted value. It should be noted that model names starting
            with log10 are for the log10 of that quantity and you should apply
            10 ** prediction to get the actual value.
        """
        response = None
        url = self.url + "/predict_structure/%s" % model_name
        try:
            data = {"structure": structure.to(fmt="POSCAR"), "fmt": "POSCAR"}
            response = self.session.post(url, data=data)
            if response.status_code in [200, 400]:
                return response.json()
            raise ValueError("REST query returned with error status code {}".format(response.status_code))
        except Exception as ex:
            msg = "{}. Content: {}".format(str(ex), response.content) if hasattr(response, "content") else str(ex)
            raise ValueError(msg)
