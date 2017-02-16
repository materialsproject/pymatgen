from pymatgen.analysis.chemenv.connectivity.structure_connectivity import StructureConnectivity
import logging

__author__ = 'waroquiers'

class ConnectivityFinder(object):

    def __init__(self, parameters=None):
        self.parameters = parameters
        #self.setup_graph()

    def setup_structure_environments(self, light_structure_environments=None):
        if light_structure_environments is not None:
            self.light_structure_environments = light_structure_environments
        else:
            raise NotImplementedError("light_structure_environments should be provided")

    def setup_graph(self):
        logging.info('Setup of structure connectivity graph')
        self.structure_connectivity = StructureConnectivity(self.light_structure_environments)
        self.structure_connectivity.add_sites()
        for isite, site in enumerate(self.light_structure_environments.structure):
            site_neighbors_sets = self.light_structure_environments.neighbors_sets[isite]
            if site_neighbors_sets is None:
                continue
            if len(site_neighbors_sets) > 1:
                raise NotImplementedError('Mix of environments (with mix of neighbors sets) are not implemented')
            site_neighbors_set = site_neighbors_sets[0]
            self.structure_connectivity.add_bonds(isite, site_neighbors_set)
    #
    # def setup_environment_subgraph(self, myelement, myenvironment):
    #     self.structure_connectivity.setup_environment_subgraph(myelement, myenvironment)
    #     self.structure_connectivity.setup_boundary_graph_and_cycles()
    #     self.structure_connectivity.get_fiber_cycles(start_cycle=-1, end_cycle=1, to_periodic_site=False)

