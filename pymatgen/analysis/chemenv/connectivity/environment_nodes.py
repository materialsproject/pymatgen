__author__ = 'waroquiers'

import abc
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from six import with_metaclass


class AbstractEnvironmentNode(with_metaclass(abc.ABCMeta)):

    COORDINATION_ENVIRONMENT = 0
    NUMBER_OF_NEIGHBOURS = 1
    NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBOUR = 2
    LIGANDS_ARRANGEMENT = 3
    NEIGHBOURS_LIGANDS_ARRANGEMENT = 4
    ATOM = 5
    DEFAULT_EXTENSIONS = [COORDINATION_ENVIRONMENT]
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

    @abc.abstractproperty
    def coordination_environment(self):
        return

    def descriptor(self, extensions=DEFAULT_EXTENSIONS):
        return self.ce_extended(extensions)

    def ce_extended(self, extensions=DEFAULT_EXTENSIONS):
        npoints = max(extensions)
        res = ['']*(npoints + 1)
        for extension in extensions:
            res[extension] = self.get_descriptor(extension)
        return '.'.join(res)

    @abc.abstractmethod
    def get_descriptor(self, extension):
        return

    @property
    def ce(self):
        return self.coordination_environment

    @property
    def mp_symbol(self):
        return self.coordination_environment

    def __str__(self):
        return 'Node #{:d} ({})'.format(self.isite, self.mp_symbol)


class OctahedralEnvironmentNode(AbstractEnvironmentNode):

    CG = AllCoordinationGeometries().get_geometry_from_mp_symbol('O:6')

    def __init__(self, central_site, i_central_site):
        AbstractEnvironmentNode.__init__(self, central_site, i_central_site)

    @property
    def coordination_environment(self):
        return self.CG.mp_symbol

    def get_descriptor(self, extension):
        if extension == AbstractEnvironmentNode.COORDINATION_ENVIRONMENT:
            return self.coordination_environment
        else:
            return 'NULL'


class TetrahedralEnvironmentNode(AbstractEnvironmentNode):

    CG = AllCoordinationGeometries().get_geometry_from_mp_symbol('T:4')

    def __init__(self, central_site, i_central_site):
        AbstractEnvironmentNode.__init__(self, central_site, i_central_site)

    @property
    def coordination_environment(self):
        return self.CG.mp_symbol

    def get_descriptor(self, extension):
        if extension == AbstractEnvironmentNode.COORDINATION_ENVIRONMENT:
            return self.coordination_environment
        else:
            return 'NULL'


allowed_environment_nodes = {'O:6': OctahedralEnvironmentNode,
                             'T:4': TetrahedralEnvironmentNode}


def get_environment_node(central_site, i_central_site, mp_symbol):
    """

    :param central_site:
    :param mp_symbol:
    :raise NotImplementedError:
    """
    if mp_symbol not in allowed_environment_nodes:
        raise NotImplementedError('Coordination environment "{}" is not yet allowed for connectivity description')
    return allowed_environment_nodes[mp_symbol](central_site, i_central_site)