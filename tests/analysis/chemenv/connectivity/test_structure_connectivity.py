from __future__ import annotations

import json
import os

from pymatgen.analysis.chemenv.connectivity.connectivity_finder import ConnectivityFinder
from pymatgen.analysis.chemenv.connectivity.structure_connectivity import StructureConnectivity
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import (
    LightStructureEnvironments,
    StructureEnvironments,
)
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

try:
    import bson
except ModuleNotFoundError:
    bson = None  # type: ignore

__author__ = "waroquiers"


class TestStructureConnectivity(PymatgenTest):
    def test_serialization(self):
        BaTiO3_se_fpath = os.path.join(
            TEST_FILES_DIR,
            "chemenv",
            "structure_environments_files",
            "se_mp-5020.json",
        )
        with open(BaTiO3_se_fpath) as file:
            dd = json.load(file)
        struct_envs = StructureEnvironments.from_dict(dd)
        lse = LightStructureEnvironments.from_structure_environments(
            strategy=SimplestChemenvStrategy(), structure_environments=struct_envs
        )
        cf = ConnectivityFinder()
        sc = cf.get_structure_connectivity(light_structure_environments=lse)
        sc_from_dict = StructureConnectivity.from_dict(sc.as_dict())
        assert sc.light_structure_environments == sc_from_dict.light_structure_environments
        assert set(sc._graph.nodes()) == set(sc_from_dict._graph.nodes())
        assert set(sc._graph.edges()) == set(sc_from_dict._graph.edges())

        sc_from_json = StructureConnectivity.from_dict(json.loads(json.dumps(sc.as_dict())))
        assert sc.light_structure_environments == sc_from_json.light_structure_environments
        assert set(sc._graph.nodes()) == set(sc_from_json._graph.nodes())
        assert set(sc._graph.edges()) == set(sc_from_json._graph.edges())

        if bson is not None:
            bson_data = bson.BSON.encode(sc.as_dict())
            sc_from_bson = StructureConnectivity.from_dict(bson_data.decode())
            assert sc.light_structure_environments == sc_from_bson.light_structure_environments
            assert set(sc._graph.nodes()) == set(sc_from_bson._graph.nodes())
            assert set(sc._graph.edges()) == set(sc_from_bson._graph.edges())
