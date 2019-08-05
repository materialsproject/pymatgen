__author__ = 'waroquiers'

import abc
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from monty.json import MSONable
from six import with_metaclass
import operator


class AbstractEnvironmentNode(MSONable):

    COORDINATION_ENVIRONMENT = 0
    NUMBER_OF_NEIGHBORING_COORDINATION_ENVIRONMENTS = 1
    NUMBER_OF_NEIGHBORING_CES = NUMBER_OF_NEIGHBORING_COORDINATION_ENVIRONMENTS
    NEIGHBORING_COORDINATION_ENVIRONMENTS = 2
    NEIGHBORING_CES = NEIGHBORING_COORDINATION_ENVIRONMENTS
    NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_COORDINATION_ENVIRONMENT = 3
    NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_CE = NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_COORDINATION_ENVIRONMENT
    LIGANDS_ARRANGEMENT = 4
    NEIGHBORS_LIGANDS_ARRANGEMENT = 5
    ATOM = 6
    CE_NNBCES_NBCES_LIGANDS = -1
    DEFAULT_EXTENSIONS = [ATOM, COORDINATION_ENVIRONMENT]
    # DEFAULT_EXTENSIONS = [COORDINATION_ENVIRONMENT,
    #                       NUMBER_OF_NEIGHBOURS,
    #                       NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBOUR]

    def __init__(self, central_site, i_central_site):
        """

        :param central_site: pymatgen Site or subclass of Site (e.g. PeriodicSite, ...)
        """
        #if not isinstance(central_site, Site):
        #    raise ClassTypeChemenvError(central_site, Site)
        self.central_site = central_site
        self.i_central_site = i_central_site

    @property
    def isite(self):
        return self.i_central_site

    def __hash__(self):
        return self.central_site.__hash__()

    def __eq__(self, other):
        return self.isite == other.isite

    def everything_equal(self, other):
        return self.__eq__(other) and self.central_site == other.central_site

    @abc.abstractproperty
    def coordination_environment(self):
        return

    def number_of_neighboring_coordination_environments(self, environments_subgraph):
        # One cannot use the MultiGraph.neighbors(self) method because for self-loops, it yields the neighbor only once
        incident_edges = environments_subgraph.edges(self)
        count = 0
        for edge in incident_edges:
            count += 1
            if edge[0] == edge[1]:
                count += 1
        return str(count)

    def neighboring_coordination_environments(self, environments_subgraph):
        pass
        # neighboring_environments_nodes = environments_subgraph.neighbors(self)
        # return str(len(neighboring_environments_nodes))

    def descriptor(self, extensions=DEFAULT_EXTENSIONS, **kwargs):
        return self.ce_extended(extensions, **kwargs)

    def ce_extended(self, extensions=DEFAULT_EXTENSIONS, **kwargs):
        npoints = max(extensions)
        res = ['']*(npoints + 1)
        for extension in extensions:
            res[extension] = self.get_descriptor(extension, **kwargs)
        return '.'.join(res)

    def get_descriptor(self, extension, **kwargs):
        if extension == AbstractEnvironmentNode.COORDINATION_ENVIRONMENT:
            return self.coordination_environment
        elif extension == AbstractEnvironmentNode.NUMBER_OF_NEIGHBORING_CES:
            env_subgraph = kwargs['environments_subgraph']
            return self.number_of_neighboring_coordination_environments(environments_subgraph=env_subgraph)
        elif extension == AbstractEnvironmentNode.CE_NNBCES_NBCES_LIGANDS:
            descr = str(self.coordination_environment)
            descr += '.'
            env_subgraph = kwargs['environments_subgraph']
            neighboring_ces_nodes = env_subgraph.neighbors(self)
            if len(list(neighboring_ces_nodes)) == 0:
                descr += '0-..'
            else:
                my_neighboring_nodes = []
                for ces_node in neighboring_ces_nodes:
                    for iedge in env_subgraph[self][ces_node]:
                        my_neighboring_nodes.append({'neighbor': ces_node,
                                                     'edge_data': env_subgraph[self][ces_node][iedge]})
                        # Special case for self-loops :
                        if ces_node == self:
                            opposed_data = dict(env_subgraph[self][ces_node][iedge])
                            opposed_data['delta'] = tuple([-ii for ii in opposed_data['delta']])
                            my_neighboring_nodes.append({'neighbor': ces_node,
                                                         'edge_data': opposed_data})
                my_neighboring_nodes.sort(key=lambda x: len(x['edge_data']['ligands']))
                descr += str(len(my_neighboring_nodes))
                descr += '-'
                descr += ','.join([str(len(nn['edge_data']['ligands'])) for nn in my_neighboring_nodes])
                descr += '.'
                descr += ','.join([str(nn['neighbor'].coordination_environment) for nn in my_neighboring_nodes])
                #TODO: check sorting according to the nomenclature !
            return descr
        elif extension == AbstractEnvironmentNode.ATOM:
            return self.atom_symbol
        else:
            return 'NULL'

    @property
    def ce(self):
        return self.coordination_environment

    @property
    def mp_symbol(self):
        return self.coordination_environment

    @property
    def atom_symbol(self):
        return self.central_site.specie.symbol

    def __str__(self):
        return 'Node #{:d} ({})'.format(self.isite, self.mp_symbol)


class EnvironmentNode(AbstractEnvironmentNode):

    def __init__(self, central_site, i_central_site, ce_symbol):
        AbstractEnvironmentNode.__init__(self, central_site, i_central_site)
        self.ce_symbol = ce_symbol

    @property
    def coordination_environment(self):
        return self.ce_symbol

    def everything_equal(self, other):
        return (super().everything_equal(other) and
                self.coordination_environment == other.coordination_environment)


class OctahedralEnvironmentNode(AbstractEnvironmentNode):

    CG = AllCoordinationGeometries().get_geometry_from_mp_symbol('O:6')

    def __init__(self, central_site, i_central_site):
        AbstractEnvironmentNode.__init__(self, central_site, i_central_site)

    @property
    def coordination_environment(self):
        return self.CG.mp_symbol


class TetrahedralEnvironmentNode(AbstractEnvironmentNode):

    CG = AllCoordinationGeometries().get_geometry_from_mp_symbol('T:4')

    def __init__(self, central_site, i_central_site):
        AbstractEnvironmentNode.__init__(self, central_site, i_central_site)

    @property
    def coordination_environment(self):
        return self.CG.mp_symbol


allowed_environment_nodes = {'O:6': OctahedralEnvironmentNode,
                             'T:4': TetrahedralEnvironmentNode}


def get_environment_node(central_site, i_central_site, ce_symbol):
    """
    Get the EnvironmentNode class or subclass for the given site and symbol.

    :param central_site: Central site of the environment
    :param i_central_site: Index of the central site in the structure
    :param ce_symbol: Environment symbol
    :return: An EnvironmentNode object
    """
    return EnvironmentNode(central_site=central_site, i_central_site=i_central_site, ce_symbol=ce_symbol)
