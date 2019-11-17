from monty.json import MSONable
from monty.io import zopen


class Section(dict, MSONable):

    required_params = []

    def __init__(self, section_name, section_keyword=None, section_params=None):
        self.section_name = section_name
        self.section_keyword = section_keyword
        self.section_params = section_params

        super(Section).__init__()
        if section_keyword:
            section_head = self.section_name+" "+self.section_keyword
        else:
            section_head = self.section_name

        self.update({section_head: self.section_params})

        for k in self.required_params:
            if k not in section_params:
                raise ValueError("{}: Required parameter {} not specified!".format(section_params, k))

    def __str__(self):
        return self.get_string()

    def __dict__(self):
        return self.as_dict()

    def get_string(self):
        return Section._get_string(self)

    @staticmethod
    def _get_string(d, indent=0):
        string = ''
        if isinstance(d, dict):
            for key, value in d.items():
                if isinstance(value, dict) or isinstance(value, list):
                    string += '\t' * indent + '&' + key.__str__() + '\n'
                    string += Section._get_string(value, indent + 1)
                    string += '\t' * indent + '&END ' + key.__str__().split()[0] + '\n'
                else:
                    string += '\t' * indent + key.__str__()+' '+value.__str__() + '\n'
        elif isinstance(d, list):
            for i in d:
                string += '\t' * indent + i.__str__() + '\n'
        return string

    def as_dict(self):
        d = dict(self)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        for k,v in d.items():
            h = k.split()
            n = h[0]
            if len(h) > 1: kw = h[1]
            else: kw=None
            return Section(section_name=n, section_keyword=kw,
                           section_params={k: v for k, v in d.items() if k not in ('@module',
                                                                               '@class')})

    def write_file(self, filename):
        with zopen(filename, 'wt') as f:
            f.write(self.__str__())


class Input(Section):

    """
    Input file object, slight extension of Sections object

    header
    comments
    sections
    """

    def get_string(self):
        return Section._get_string(self[self.section_name])


class Global(Section):

    required_params = []

    def __init__(self, project='CP2K_PROJECT', run_type='GEO_OPT',
                 print_level='MEDIUM', params={}):

        d = {}

        d['PROJECT']  = project
        d['RUN_TYPE'] = run_type
        d['PRINT_LEVEL'] = print_level

        super(Global, self).__init__(section_name='GLOBAL',
                                     section_keyword=None,
                                     section_params=d)


class ForceEval(Section):

    required_params = []

    def __init__(self, params={}):

        d = {}

        d['METHOD'] = params.get('METHOD', 'QUICKSTEP')
        d['DFT'] = params.get('DFT', Dft())
        d['SUBSYS'] = params.get('SUBSYS', Subsys())

        d.update(params)

        super(ForceEval, self).__init__(section_name='FORCE_EVAL',
                                        section_keyword=None,
                                        section_params=d)


class Dft(Section):

    required_params = []

    def __init__(self, params={}):

        d = {}

        d['BASIS_SET_FILENAME'] = params.get('BASIS_SET_FILE_NAME', 'GTH_BASIS_SETS')
        d['POTENTIAL_FILE_NAME'] = params.get('POTENTIAL_FILE_NAME', 'GTH_POTENTIALS')
        d.update(params.get('QS', QS()))

        params.update(d)

        super(Dft, self).__init__(section_name='DFT',
                                  section_keyword=None,
                                  section_params=params)


class Subsys(Section):

    def __init__(self, params={}):

        super(Subsys, self).__init__(section_name='SUBSYS',
                                     section_keyword=None,
                                     section_params=params)


class Mgrid(Section):

    def __init__(self, params={}):

        d = {}

        d['NGRIDS'] = params.get('NGRIDS', 4)
        d['CUTOFF'] = params.get('CUTOFF', 600)
        d['REL_CUTOFF'] = params.get('REL_CUTOFF', 60)

        params.update(d)

        super(Mgrid, self).__init__(section_name='MGRID',
                                    section_keyword=None,
                                    section_params=params)


class QS(Section):

    required_params = []

    def __init__(self, params={}):

        d = {}

        d['METHOD'] = params.get('METHOD', 'GPW')
        d['EPS_DEFAULT'] = params.get('EPS_DEFAULT', 1.0E-7)
        d['EXTRAPOLATION'] = params.get('EXTRAPOLATION', 'ASPC')

        params.update(d)

        super(QS, self).__init__(section_name='QS',
                                 section_keyword=None,
                                 section_params=params)


class Poisson(Section):

    required_params = []

    def __init__(self, params={}):

        d = {}

        d['PERIODIC'] = params.get('PERIODIC', 'XYZ')

        params.update(d)

        super(Poisson, self).__init__(section_name='Poisson',
                                      section_keyword=None,
                                      section_params=params)


class Scf(Section):

    def __init__(self, params={}):

        d = {}

        d['SCF_GUESS'] = params.get('SCF_GUESS', 'ATOMIC')
        d['EPS_SCF'] = params.get('EPS_SCF', 1.0E-5)
        d['MAX_SCF'] = params.get('MAX_SCF', 50)

        params.update(d)

        super(Scf, self).__init__(section_name='SCF',
                                  section_keyword=None,
                                  section_params=params)


class OrbitalTransformation(Section):

    def __init__(self, params={}):

        d = {}

        d['PRECONDITIONER'] = params.get('PRECONDITIONER', 'FULL_ALL')
        d['MINIMIZER'] = params.get('MINIMIZER', 'DIIS')

        params.update(d)

        super(OrbitalTransformation, self).__init__(section_name='OT',
                                                    section_keyword=None,
                                                    section_params=params)


class Cell(Section):

    def __init__(self, lattice):

        A = "".join([" ".join(["%.6f" % i for i in lattice.matrix[0]])]) + "\n"
        B = "".join([" ".join(["%.6f" % i for i in lattice.matrix[1]])]) + "\n"
        C = "".join([" ".join(["%.6f" % i for i in lattice.matrix[2]])]) + "\n"

        d = {'A': A, 'B': B, 'C': C}

        super(Cell, self).__init__(section_name='CELL',
                                   section_keyword=None,
                                   section_params=d)


class Kind(Section):

    def __init__(self, params=None):

        specie    = params.get('specie')
        basis_set = params.get('BASIS_SET')
        potential = params.get('POTENTIAL')

        d = {'Element': specie,
             'BASIS_SET': basis_set,
             'POTENTIAL': potential
             }

        super(Kind, self).__init__(section_name='KIND',
                                   section_keyword=specie,
                                   section_params=d)


class Coord(Section):

    def __init__(self, structure):

        coords = [s.specie.symbol + " " + " ".join([c.__str__() for c in s.coords]) for s in structure.sites]

        super(Coord, self).__init__(section_name='COORDS',
                                     section_keyword=None,
                                     section_params=coords)
