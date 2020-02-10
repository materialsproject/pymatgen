import os
import copy
import textwrap
from monty.json import MSONable

# TODO: Use multisets for keywords?
# Todo: update method for section?
# TODO: Location/Dependencies
# TODO can either use kwargs accessed with kwargs.get(,default) or actual default keyword arguments. Which one?


class Section(MSONable):

    # TODO: Not sure if I want to use this required_sections/keywords idea
    required_sections = []
    required_keywords = []
    subsections = {}  # Must exist before __init__ is even called

    def __init__(self, name, subsections, repeats=False, description=None, keywords=[],
                 section_parameters=[], location=None, verbose=True, **kwargs):
        """
        Basic object representing a CP2K Section. Sections activate different parts of the calculation.
        For example, FORCE_EVAL section will activate CP2K's ability to calculate forces. Sections are
        described with the following:

        Args:
            name (str): The name of the section (must match name in CP2K)
            subsections (dict): A dictionary of subsections that are nested in this section. Format is
                {'NAME': Section(*args, **kwargs). The name you chose for 'NAME' to index that subsection
                does not *have* to be the same as the section's true name, but we recommend matching
                them. You can specify a blank dictionary if there are no subsections, or if you want to
                insert the subsections later.
            repeats (bool): Whether or not this section can be repeated. Most sections cannot. Default=False.
            description (str): Description of this section for easier readability/more descriptiveness
            keywords (list): the keywords to be set for this section. Each element should be a Keyword
                object. This can be more cumbersome than simply using kwargs for building a class in a script,
                but is more convenient for the class instantiations of CP2K sections (see below).
            section_parameters (list): the section parameters for this section. Section parameters are
                specialized keywords that modify the behavior of the section overall. Most
                sections do not have section parameters, but some do. Unlike normal Keywords, these are
                specified as strings and not as Keyword objects.
            location (str): the path to the section in the form 'SECTION/SUBSECTION1/SUBSECTION3', example for
                QS module: 'FORCE_EVAL/DFT/QS'. This location is used to automatically determine if a subsection
                requires a supersection to be activated.
            verbose (str): Controls how much is printed to Cp2k input files (Also see Keyword). If True, then
                a description of the section will be printed with it as a comment (if description is set).
                Default=True.

            kwargs are interpreted as keyword, value pairs and added to the keywords array as Keyword objects
        """

        self.name = name
        self.repeats = repeats
        self.description = description
        self.keywords = keywords
        for k, v in kwargs.items():
            self.keywords.append(Keyword(k, v))
        self.section_parameters = section_parameters
        self.subsections = subsections
        self.location = location
        self.verbose = verbose

        for k in self.required_sections:
            if not self.check(k):
                raise UserWarning("WARNING: REQUIRED SECTION {} HAS NOT BEEN INITIALIZED".format(k))
        for k in self.required_keywords:
            if k not in self.keywords:
                raise UserWarning("WARNING: REQUIRED KEYWORD {} HAS NOT BEEN PROVIDED".format(k))

    def __str__(self):
        return self.get_string()

    def __dict__(self):
        return self.as_dict()

    def __deepcopy__(self, memodict={}):
        return Section.from_dict(copy.deepcopy(self.as_dict()))

    def __getitem__(self, d):
        return self.subsections[d]

    def __setitem__(self, key, value):
        if isinstance(value, Section):
            self.subsections[key] = value.__deepcopy__()
        else:
            self.subsections[key] = Section(value, subsections={})

    def __delitem__(self, key):
        if key in self.subsections.keys():
            del self.subsections[key]
        else:
            for i, k in enumerate(self.keywords):
                if k.name is key:
                    del self.keywords[i]

    def get_keyword(self, kwd):
        for i, k in enumerate(self.keywords):
            if k.name == kwd:
                return i, k

    def set(self, d):
        return self.update(d)

    def pop(self, d):
        for k, v in d.items():
            if isinstance(v, dict):
                self[k].pop(v)
            else:
                self[k].__delitem__(v)

    def as_dict(self):
        d = {}
        d["@module"] = self.__class__.__module__
        d["@module"] = self.__class__.__name__
        d['name'] = self.name
        d['repeats'] = self.repeats
        d['description'] = self.description
        d['keywords'] = []
        for k in self.keywords:
            d['keywords'].append(k.as_dict())
        d['section_parameters'] = self.section_parameters
        d['subsections'] = {}
        for k, v in self.subsections.items():
            d['subsections'][v.name] = v.as_dict()
        return d

    @classmethod
    def from_dict(cls, d):
        subsections = {}
        keywords = []
        for v in d['subsections'].values():
            subsections[v['name']] = Section.from_dict(v)
        for k in d['keywords']:
            keywords.append(Keyword.from_dict(copy.deepcopy(k)))
        return Section(d['name'], repeats=d['repeats'], description=d['description'],
                       keywords=keywords, section_parameters=d['section_parameters'],
                       subsections=subsections)

    def update(self, d):
        """
        Update the Section according to a dictionary argument. This is most useful
        for providing user-override settings to default parameters. As you pass a
        dictionary the class variables like "description", "location", or "repeats"
        are not included. Therefore, it is recommended that this be used to modify
        existing Section objects to a user's needs, but not used for the creation
        of new Section child-classes.

        Args:
            d (dict): A dictionary containing the update information. Should use nested dictionaries to
                specify the full path of the update. If a section or keyword does not exist, it will be created,
                but only with the values that are provided in "d", not using default values from a Section object.

                {'SUBSECTION1': {'SUBSEC2': {'NEW_KEYWORD', 'NEW_VAL'},{'NEW_SUBSEC': {'NEW_KWD': 'NEW_VAL'}}}
        """
        Section._update(self, d)

    @staticmethod
    def _update(d1, d2):
        """
        Helper method for self.update(d) method (see above).
        """
        for k, v in d2.items():
            if isinstance(v, (str, float, bool)):
                kwds = [kwd.name for kwd in d1.keywords]
                if k in kwds:
                    i = kwds.index(k)
                    d1.keywords[i] = Keyword(k, v)
                else:
                    d1.add(Keyword(k, v))
            elif isinstance(v, dict):
                if k not in list(d1.subsections.keys()):
                    d1.insert(Section(k, subsections={}))
                return Section._update(d1.subsections[k], v)

    def insert(self, d):
        """
        Insert a new section as a subsection of the current one
        """
        self.subsections[d.name] = d.__deepcopy__()

    def add(self, d):
        """
        Add a keyword to the section keywords.
        """
        self.keywords.append(d)

    def check(self, path):
        """
        Check if section exists within the current. Can be useful for cross-checking whether or not
        required dependencies have been satisfied, which CP2K does not enforce.

        Args:
            path (str): Path to section of form 'SUBSECTION1/SUBSECTION2/SUBSECTION_OF_INTEREST'
        """
        path = path.split('/')
        s = self.subsections
        for p in path:
            if p in s.keys():
                s = s[p].subsections
            else:
                return False
        return True

    def by_path(self, path):
        """
        Access a sub-section using a path. Used by the file parser.

        Args:
            path (str): Path to section of form 'SUBSECTION1/SUBSECTION2/SUBSECTION_OF_INTEREST'

        """
        path = path.split('/')
        if path[0] == self.name:
            path = path[1:]
        s = self
        for p in path:
            s = s[p]
        return s

    def get_string(self):
        return Section._get_string(self)

    @staticmethod
    def _get_string(d, indent=0):
        """
        Helper function to return a pretty string of the section. Includes indentation and
        descriptions (if present).
        """
        string = ''
        if d.description and d.verbose:
            string += '\n' + textwrap.fill(d.description, initial_indent='\t' * indent + '! ',
                                           subsequent_indent='\t' * indent + '! ', width=50) + '\n'
        string += '\t' * indent + '&' + d.name
        string += ' ' + ' '.join(map(str, d.section_parameters)) + '\n'

        for s in d.keywords:
            string += '\t' * (indent+1) + s.__str__() + '\n'
        for k,v in d.subsections.items():
            string += v._get_string(v, indent + 1)
        string += '\t' * indent + '&END ' + '\n'

        return string

    def silence(self):
        """
        Silence the section recursively by turning off the printing of descriptions. This may reduce
        the appealing documentation of input files, but also helps de-clutter them.
        """
        self.verbose = False
        for kwd in self.keywords:
            kwd.silence()
        for k, v in self.subsections.items():
            v.silence()


class Keyword:

    def __init__(self, name, *args, description=None, units=None,
                 verbose=True, repeats=False):
        """
        Class representing a keyword argument in CP2K. Within CP2K Secitons, which activate features
        of the CP2K code, the keywords are arguments that control the functionality of that feature.
        For example, the section "FORCE_EVAL" activates the evaluation of forces/energies, but within
        "FORCE_EVAL" the keyword "METHOD" controls whether or not this will be done with, say,
        "Quickstep" (DFT) or "EIP" (empirical interatomic potential). These Keywords and the value passed
        to them are sometimes as simple as KEYWORD VALUE, but can also be more elaborate such as
        KEYWORD [UNITS] VALUE1 VALUE2, which is why this class exists: to handle many values and control
        easy printing to an input file.

        Args:
            name (str): The name of this keyword. Must match an acceptable keyword from CP2K
            args: All non-keyword arguments after 'name' are interpreted as the values to set for
                this keyword. i.e: KEYWORD ARG1 ARG2 would provide to values to the keyword.
            description (str): The description for this keyword. This can make readability of
                input files easier for some. Default=None.
            units (str): The units for this keyword. If not specified, CP2K default units will be
                used. Consult manual for default units. Default=None.
            repeats (bool): Whether or not this keyword may be repeated. Default=False.
        """
        self.name = name
        self.values = [a for a in args]
        self.description = description
        self.repeats = repeats
        self.units = units
        self.verbose = verbose

    def __str__(self):
        return self.name.__str__() + ' ' + \
               ('[{}] '.format(self.units) if self.units else "") + \
               ' '.join(map(str, self.values)) + \
               (' ! '+self.description if (self.description and self.verbose) else "")

    def __eq__(self, other):
        if self.name == other.name:
            if self.values == other.values:
                if self.units == self.units:
                    return True
        return False

    def as_dict(self):
        d = {}
        d['name'] = self.name
        d['values'] = self.values
        d['description'] = self.description
        d['repeats'] = self.repeats
        d['units'] = self.units
        d['verbose'] = self.verbose
        return d

    @classmethod
    def from_dict(cls, d):
        return Keyword(d['name'], *d['values'], description=d['description'],
                       repeats=d['repeats'], units=d['units'], verbose=d['verbose'])

    @classmethod
    def from_string(self, s):
        return Keyword(*s.split())

    def silence(self):
        """
        Deactivate the printing of this keyword's description.
        """
        self.verbose = False


class Cp2kInput(Section):

    def __init__(self, name='CP2K_INPUT', subsections={}, **kwargs):
        """
        Special instance of 'Section' class that is meant to represent the overall cp2k input.
        Distinguishes itself from Section by overriding get_string() to not print this section's
        title and by implementing the file i/o
        """

        description = "CP2K Input"
        super(Cp2kInput, self).__init__(name, repeats=False, description=description, keywords=[],
                                        section_parameters=[], subsections=subsections,
                                        **kwargs)

    def get_string(self):
        s = ''
        for k in self.subsections.keys():
            s += self.subsections[k].get_string()
        return s

    @classmethod
    def from_dict(cls, d):
        return Cp2kInput('CP2K_INPUT', subsections=Section.from_dict(d).subsections)

    @staticmethod
    def from_file(file):
        with open(file) as f:
            lines = f.read().splitlines()
            lines = [line.replace('\t', '') for line in lines]
            lines = [line for line in lines if line]
            return Cp2kInput.from_lines(lines)

    @classmethod
    def from_lines(cls, lines):
        cp2k_input = Cp2kInput('CP2K_INPUT', subsections={})
        Cp2kInput._from_lines(cp2k_input, lines)
        return cp2k_input

    def _from_lines(self, lines):
        current = self.name
        for line in lines:
            if line.startswith('!') or line.startswith('#'):
                continue
            elif line.startswith('&END'):
                current = '/'.join(current.split('/')[:-1])
            elif line.startswith('&'):
                name, subsection_params = line.split()[0][1:], line.split()[1:]
                s = Section(name, section_parameters=subsection_params,
                            subsections={})
                self.by_path(current).insert(s)
                current = current+'/'+name
            else:
                args = line.split()
                self.by_path(current).keywords.append(Keyword(*args))

    def write_file(self, input_filename='cp2k.inp', output_dir='.', make_dir_if_not_present=True):
        if not os.path.isdir(output_dir) and make_dir_if_not_present:
            os.mkdir(output_dir)
        filepath = os.path.join(output_dir, input_filename)
        with open(filepath, 'w') as f:
            f.write(self.get_string())


class Global(Section):

    def __init__(self, project_name='CP2K', run_type='ENERGY_FORCE', subsections={}, **kwargs):
        description = 'Section with general information regarding which kind of simulation' + \
                      'to perform an parameters for the whole PROGRAM'

        keywords = [
            Keyword('PROJECT_NAME', project_name),
            Keyword('RUN_TYPE', run_type)
        ]

        super().__init__('GLOBAL', description=description,
                         keywords=keywords, subsections=subsections, **kwargs)


class ForceEval(Section):

    def __init__(self, subsections={}, **kwargs):
        description = 'parameters needed to calculate energy and forces and describe the system you want to analyze.'

        keywords = [
            Keyword('METHOD', kwargs.get('METHOD', 'QS'))
        ]

        super(ForceEval, self).__init__('FORCE_EVAL', repeats=True, description=description,
                                        keywords=keywords, subsections=subsections, **kwargs)


class Dft(Section):

    def __init__(self, subsections={}, **kwargs):

        description = 'parameter needed by dft programs'

        keywords = [
            Keyword('BASIS_SET_FILE_NAME', kwargs.get('BASIS_SET_FILE_NAME', 'BASIS_MOLOPT')),
            Keyword('POTENTIAL_FILE_NAME', kwargs.get('POTENTIAL_FILE_NAME', 'GTH_POTENTIALS')),
            Keyword('UKS', kwargs.get('UKS', True))
        ]

        super().__init__('DFT', description=description, keywords=keywords,
                         subsections=subsections, **kwargs)


class Subsys(Section):

    def __init__(self, subsections={}, **kwargs):
        description = 'a subsystem: coordinates, topology, molecules and cell'
        super(Subsys, self).__init__('SUBSYS', description=description,
                                     subsections=subsections, **kwargs)


class QS(Section):

    def __init__(self, method='GPW', eps_default=1e-7, extrapolation='ASPC',
                 subsections={}, **kwargs):

        description = 'parameters needed to set up the Quickstep framework'

        keywords = [
            Keyword('METHOD', method),
            Keyword('EPS_DEFAULT', eps_default),
            Keyword('EXTRAPOLATION', extrapolation)
        ]

        super(QS, self).__init__('QS', description=description, keywords=keywords,
                                 subsections=subsections, **kwargs)


class Scf(Section):

    def __init__(self, max_scf=50, eps_scf=1e-6, scf_guess='RESTART',
                 subsections={}, **kwargs):

        description = 'Parameters needed to perform an SCF run.'

        keywords = [
            Keyword('MAX_SCF', max_scf),
            Keyword('EPS_SCF', eps_scf),
            Keyword('SCF_GUESS', scf_guess)  # Uses Restart file if present, and ATOMIC if not present
        ]

        super(Scf, self).__init__('SCF', description=description, keywords=keywords,
                                  subsections=subsections, **kwargs)


class Mgrid(Section):

    def __init__(self, cutoff=1200, rel_cutoff=80, ngrids=5,
                 progression_factor=3, multigrid_set=False,
                 subsections={}, **kwargs):

        description = 'Multigrid information. Multigrid allows for sharp gaussians and diffuse ' + \
                        'gaussians to be treated on different grids, where the spacing of FFT integration ' + \
                        'points can be tailored to the degree of sharpness/diffusiveness of the gaussians.'

        keywords = [
            Keyword('CUTOFF', cutoff, description='Cutoff in [Ry] for finest level of the MG.'),
            Keyword('REL_CUTOFF', rel_cutoff,
                    description='Controls which gaussians are mapped to which level of the MG'),
            Keyword('NGRIDS', ngrids),
            Keyword('PROGRESSION_FACTOR', progression_factor),
            Keyword('MULTIGRID_SET', multigrid_set)
        ]

        super(Mgrid, self).__init__('MGRID', description=description, keywords=keywords,
                                    subsections=subsections, **kwargs)


class Diagonalization(Section):
    
    def __init__(self, eps_adapt=0, eps_iter=1e-8, eps_jacobi=0,
                 jacobi_threshold=1e-7, max_iter=2, subsections={}):

        location = 'CP2K_INPUT/FORCE_EVAL/DFT/SCF/DIAGONALIZATION'

        keywords = [
            Keyword('EPS_ADAPT', eps_adapt),
            Keyword('EPS_ITER', eps_iter),
            Keyword('EPS_JACOBI', eps_jacobi),
            Keyword('JACOBI_THRESHOLD', jacobi_threshold)
        ]

        super(Diagonalization, self).__init__('DIAGONALIZATION', keywords=keywords,
                                              subsections=subsections, repeats=False,
                                              location=location)


class OrbitalTransformation(Section):

    def __init__(self, minimizer='CG', preconditioner='FULL_ALL', energy_gap=0.01, **kwargs):
        """
        Turns on the Orbital Transformation scheme for diagonalizing the Hamiltonian. Much faster and with
        guaranteed convergence compared to normal diagonalization, but requires the system to have a band
        gap.

        NOTE: OT has poor convergence for metallic systems and cannot use SCF mixing or smearing. Therefore,
        you should not use it for metals or systems with 'small' band gaps. In that case, use normal
        diagonalization, which will be slower, but will converge properly.

        Args:
            minimizer (str): The minimizer to use with the OT method. Default is conjugate gradient method,
                which is more robust, but more well-behaved systems should use DIIS, which can be as much
                as 50% faster.
            preconditioner (str): Preconditionar to use for OT, FULL_ALL tends to be most robust, but is
                not always most efficient. Consult the manual.
            energy_gap (float): Guess for the band gap. Should be smaller than the actual band gap, so simply
                using 0.01 is a robust value. Choosing a larger value will help if you start with a bad
                initial guess though.
        """

        description = 'Sets the various options for the orbital transformation (OT) method. ' + \
                      'Default settings already provide an efficient, yet robust method. Most ' + \
                      'systems benefit from using the FULL_ALL preconditioner combined with a small ' + \
                      'value (0.001) of ENERGY_GAP. Well-behaved systems might benefit from using ' + \
                      'a DIIS minimizer. Advantages: It\'s fast, because no expensive diagonalisation' + \
                      'is performed. If preconditioned correctly, method guaranteed to find minimum. ' + \
                      'Disadvantages: Sensitive to preconditioning. A good preconditioner can be ' + \
                      'expensive. No smearing, or advanced SCF mixing possible: POOR convergence for ' + \
                      'metalic systems.'

        keywords = [
            Keyword('MINIMIZER', minimizer),
            Keyword('PRECONDITIONER', preconditioner),
            Keyword('ENERGY_GAP', energy_gap)
        ]
        
        super(OrbitalTransformation, self).__init__('OT', description=description, keywords=keywords,
                                                    subsections={}, **kwargs)


class Cell(Section):

    def __init__(self, lattice, subsections={}, **kwargs):

        description = 'Input parameters needed to set up the CELL.'

        keywords = [
            Keyword('A', *lattice.matrix[0]),
            Keyword('B', *lattice.matrix[1]),
            Keyword('C', *lattice.matrix[2])
        ]

        super(Cell, self).__init__('CELL', description=description, keywords=keywords,
                                   subsections=subsections, **kwargs)


class Kind(Section):

    def __init__(self, specie, alias=None, magnetization=0.0,
                 basis_set='GTH_BASIS', potential='GTH_POTENTIALS',
                 **kwargs):
        """
        Specifies the information for the different atom types being simulated.

        Args:
            specie (Species or Element): Object representing the atom.
            alias (str): Alias for the atom, can be used for specifying modifcations
                to certain atoms but not all, e.g. Mg_1 and Mg_2 to force difference
                oxidation states on the two atoms.
            magnetization (float): From the CP2K Manual: The magnetization used
                in the atomic initial guess. Adds magnetization/2 spin-alpha
                electrons and removes magnetization/2 spin-beta electrons.
            basis_set (str): Basis set for this atom, accessible from the
                basis set file specified
            potential (str): Pseudopotential for this atom, accessible from the
                potential file
            kwargs: Additional kwargs to pass to Section()
        """

        description = 'The description of the kind of the atoms (mostly for QM)'

        keywords = [
            Keyword('ELEMENT', specie.__str__()),
            Keyword('MAGNETIZATION', magnetization),
            Keyword('BASIS_SET', basis_set),
            Keyword('POTENTIAL', potential)
        ]

        section_param = alias if alias else specie.__str__()

        super(Kind, self).__init__('KIND', subsections={}, description=description, keywords=keywords,
                                   section_parameters=[section_param], **kwargs)


class Coord(Section):

    def __init__(self, structure, **kwargs):
        """
        Specifies the coordinates of the atoms using a pymatgen structure object.

        Args:
            structure: Pymatgen structure object
            kwargs: additional kwargs to pass to Section()
        """

        description = 'The coordinates for simple systems (like small QM cells) are specified ' + \
                      'here by default using explicit XYZ coordinates. More complex systems ' + \
                      'should be given via an external coordinate file in the SUBSYS%TOPOLOGY section.'

        keywords = [Keyword(s.specie.symbol, *s.coords) for s in structure.sites]

        super(Coord, self).__init__('COORD', description=description, keywords=keywords,
                                    subsections={}, **kwargs)


class PDOS(Section):

    def __init__(self, nlumo=-1):

        description = "Controls printing of the projected density of states"

        keywords = [
            Keyword('NLUMO', nlumo),
            Keyword('COMPONENTS')
        ]

        super(PDOS, self).__init__('PDOS', description=description,
                                   keywords=keywords, subsections={})


class LDOS(Section):

    def __init__(self):

        description = "Controls printing of the projected density of states decomposed by atom type"

        keywords = [
            Keyword('COMPONENTS')
        ]

        super(LDOS, self).__init__('LDOS', description=description,
                                   keywords=keywords, subsections={})


class V_Hartree_Cube(Section):

    def __init__(self):

        description = "Controls the printing of a cube file with eletrostatic potential generated by " + \
                      "the total density (electrons+ions). It is valid only for QS with GPW formalism. " + \
                      "Note that by convention the potential has opposite sign than the expected physical one."

        keywords = []

        super(V_Hartree_Cube, self).__init__('V_HARTREE_CUBE', description=description,
                                             keywords=keywords, subsections={})


class MO_Cubes(Section):

    def __init__(self, write_cube=False, nhomo=1, nlumo=1):

        description = "Controls the printing of a cube file with eletrostatic potential generated by " + \
                      "the total density (electrons+ions). It is valid only for QS with GPW formalism. " + \
                      "Note that by convention the potential has opposite sign than the expected physical one."

        keywords = [
            Keyword('WRITE_CUBE', write_cube),
            Keyword("NHOMO", nhomo),
            Keyword("NLUMO", nlumo)
        ]

        super(MO_Cubes, self).__init__('MO_CUBES', description=description,
                                       keywords=keywords, subsections={})


class E_Density_Cube(Section):

    def __init__(self):

        description = "Controls the printing of cube files with the electronic density and, for LSD " + \
                      "calculations, the spin density."

        keywords = []

        super(E_Density_Cube, self).__init__('E_DENSITY_CUBE', description=description,
                                             keywords=keywords, subsections={})

class Smear(Section):

    def __init__(self, elec_temp=300, method='FERMI_DIRAC', fixed_magnetic_moment=-1e2):

        description = "Activates smearing of electron occupations"

        keywords = [
            Keyword('ELEC_TEMP', elec_temp),
            Keyword('METHOD', method),
            Keyword('FIXED_MAGNETIC_MOMENT', fixed_magnetic_moment)
        ]

        super(Smear, self).__init__('SMEAR', description=description,
                                    subsections={}, keywords=keywords)


class BrokenSymmetry(Section):

    def __init__(self, L_alpha = -1, N_alpha = 0, NEL_alpha = -1,
                 L_beta = -1, N_beta = 0, NEL_beta = -1):
        """
        Define the required atomic orbital occupation assigned in initialization
        of the density matrix, by adding or subtracting electrons from specific
        angular momentum channels. It works only with GUESS ATOMIC

        Args:
            L_alpha: Angular momentum quantum number of the orbitals whose occupation is changed
            N_alpha: Principal quantum number of the orbitals whose occupation is changed.
                Default is the first not occupied
            NEL_alpha: Orbital occupation change per angular momentum quantum number. In
                unrestricted calculations applied to spin alpha
            L_beta: Same as L_alpha for beta channel
            N_beta: Same as N_alpha for beta channel
            NEL_beta: Same as NEL_alpha for beta channel
        """

        description = 'Define the required atomic orbital occupation assigned in initialization ' + \
                      'of the density matrix, by adding or subtracting electrons from specific ' + \
                      'angular momentum channels. It works only with GUESS ATOMIC'

        keywords_alpha = [
            Keyword('L', L_alpha),
            Keyword('N', N_alpha),
            Keyword('NEL', NEL_alpha)
        ]
        alpha = Section('ALPHA', keywords=keywords_alpha, subsections={}, repeats=False)

        keywords_beta = [
            Keyword('L', L_beta),
            Keyword('N', N_beta),
            Keyword('NEL', NEL_beta)
        ]
        beta = Section('BETA', keywords=keywords_beta, subsections={}, repeats=False)

        super(BrokenSymmetry, self).__init__('BS', description=description,
                                             subsections={'ALPHA': alpha, 'BETA': beta},
                                             keywords=[], repeats=False)


class XC_FUNCTIONAL(Section):

    def __init__(self, functional, subsections={}):

        location = 'CP2K_INPUT/FORCE_EVAL/DFT/XC/XC_FUNCTIONAL'

        built_in = [
            'BL3YLP',
            'BEEFVDW',
            'BLYP',
            'BP',
            'HCTH120',
            'LDA',
            'NONE',
            'NO_SHORTCUT',
            'OLYP',
            'PADE',
            'PBE',
            'PBE0',
            'TPSS'
        ]

        if functional not in built_in:
            raise Warning("The selected functional does not exist in CP2K.")

        super(XC_FUNCTIONAL, self).__init__('XC_FUNCTIONAL', keywords=[], subsections=subsections,
                                            location=location, repeats=False,
                                            section_parameters=[functional])


class PBE(Section):

    def __init__(self, parameterization='ORIG', scale_c=1, scale_x=1):
        """
        Info about the PBE functional.
        
        Args:
            parameterization (str):
                ORIG: original PBE
                PBESOL: PBE for solids/surfaces
                REVPBE: revised PBE
            scale_c (float): scales the correlation part of the functional.
            scale_x (float): scales the exchange part of the functional.
        """

        location = 'CP2K_INPUT/FORCE_EVAL/DFT/XC/XC_FUNCTIONAL/PBE'

        keywords = [
            Keyword('PARAMETRIZATION', parameterization),
            Keyword('SCALE_C', scale_c),
            Keyword('SCALE_X', scale_x)
        ]

        super(PBE, self).__init__('PBE', subsections={}, repeats=False,
                                  location=location, section_parameters=[],
                                  keywords=keywords)



