from __future__ import annotations

import json

from pymatgen.analysis.chemenv.connectivity.connectivity_finder import ConnectivityFinder
from pymatgen.analysis.chemenv.connectivity.structure_connectivity import StructureConnectivity
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import (
    LightStructureEnvironments,
    StructureEnvironments,
)
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "waroquiers"


class TestStructureConnectivity(PymatgenTest):
    def test_serialization(self):
        BaTiO3_se_fpath = f"{TEST_FILES_DIR}/analysis/chemenv/structure_environments/se_mp-5020.json"
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

        json_str = self.assert_msonable(sc)
        sc_from_json = StructureConnectivity.from_dict(json.loads(json_str))
        assert sc.light_structure_environments == sc_from_json.light_structure_environments
        assert set(sc._graph.nodes()) == set(sc_from_json._graph.nodes())
        assert set(sc._graph.edges()) == set(sc_from_json._graph.edges())
