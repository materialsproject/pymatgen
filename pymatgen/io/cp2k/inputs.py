"""
This module defines the building blocks of a CP2K input file. The cp2k input structure is essentially a collection
of "sections" which are similar to dictionary objects that activate modules of the cp2k executable, and then
"keywords" which adjust variables inside of those modules. For example, FORCE_EVAL section will activate CP2K's ability
to calculate forces, and inside FORCE_EVAL, the Keyword "METHOD can be set to "QS" to set the method of force evaluation
to be the quickstep (DFT) module.

A quick overview of the module:

-- Section class defines the basis of Cp2k input and contains methods for manipulating these objects similarly to Dicts.
-- Keyword class defines the keywords used inside of Section objects that changes variables in Cp2k program
-- Cp2kInput class is special instantiation of Section that is used to represent the full cp2k calculation input.
-- The rest of the classes are children of Section intended to make initialization of common sections easier.
"""

import os
import copy
import textwrap
from monty.json import MSONable
from monty.io import zopen
from pymatgen.io.cp2k.utils import _postprocessor

__author__ = "Nicholas Winner"
__version__ = "0.2"
__email__ = "nwinner@berkeley.edu"
__date__ = "January 2019"


class Section(MSONable):

    """
    Basic input representation of input to Cp2k. Activates functionality inside of the Cp2k executable.
    """

    # TODO: Not sure if I want to use this required_sections/keywords idea
    required_sections: tuple = ()
    required_keywords: tuple = ()

    def __init__(
        self,
        name,
        subsections=None,
        repeats=False,
        description=None,
        keywords=[],
        section_parameters=[],
        location=None,
        verbose=True,
        alias=None,
        **kwargs
    ):
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
        self.subsections = subsections if subsections is not None else {}
        self.repeats = repeats
        self.description = description
        self.keywords = keywords
        self.section_parameters = section_parameters
        self.location = location
        self.verbose = verbose
        self.alias = alias
        self.kwargs = kwargs
        for k, v in self.kwargs.items():
            self.keywords.append(Keyword(k, v))

        for k in self.required_sections:
            if not self.check(k):
                raise UserWarning(
                    "WARNING: REQUIRED SECTION {} HAS NOT BEEN INITIALIZED".format(
                        k
                    )
                )
        for k in self.required_keywords:
            if k not in self.keywords:
                raise UserWarning(
                    "WARNING: REQUIRED KEYWORD {} HAS NOT BEEN PROVIDED".format(
                        k
                    )
                )

    def __str__(self):
        return self.get_string()

    def __dict__(self):
        return self.as_dict()

    def __deepcopy__(self, memodict={}):
        c = copy.deepcopy(self.as_dict())
        return getattr(
            __import__(c['@module'], globals(),
                       locals(), c['@class'], 0),
            c['@class']
        ).from_dict(copy.deepcopy(self.as_dict()))

    def __getitem__(self, d):
        for k in self.keywords:
            if k.name.upper() == d.upper():
                return k.values[0] if len(k.values) == 1 else k.values
        return self.subsections[d]

    def __setitem__(self, key, value):
        if isinstance(value, Section):
            if key in self.subsections:
                self.subsections[key] = value.__deepcopy__()
            else:
                self.insert(value)
        if isinstance(value, str):
            value = Keyword(key, value)
        if isinstance(value, Keyword):
            for i, k in enumerate(self.keywords):
                if k.name == value.name:
                    self.keywords[i] = value
                    return
            self.keywords.append(value)

    def __delitem__(self, key):
        if key in self.subsections.keys():
            del self.subsections[key]
        else:
            for i, k in enumerate(self.keywords):
                if k.name is key:
                    del self.keywords[i]

    def __eq__(self, other):
        return self.as_dict().__eq__(other.as_dict())

    def get_keyword(self, kwd):
        """
        Retrieve a keyword given its name
        """
        for i, k in enumerate(self.keywords):
            if k.name.lower() == kwd.lower():
                return k

    def set_keyword(self, kwd):
        """
        Set a keyword given another keyword with the same name
        """
        for i, k in enumerate(self.keywords):
            if k.name == kwd.name:
                self.keywords[i] = kwd
                return
        self.keywords.append(kwd)

    def set(self, d):
        """
        Alias for update
        """
        return self.update(d)

    def pop(self, d):
        """
        pop function (dictionary) that deletes from a Section according to a dictionary of what to delete
        """
        for k, v in d.items():
            if isinstance(v, dict):
                self[k].pop(v)
            else:
                self[k].__delitem__(v)

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
                kwds = [kwd.name.upper() for kwd in d1.keywords]
                if k.upper() in kwds:
                    i = kwds.index(k)
                    d1.keywords[i] = Keyword(k, v)
                else:
                    d1.add(Keyword(k, v))
            elif isinstance(v, dict):
                if k not in list(d1.subsections.keys()):
                    d1.insert(Section(k, subsections={}))
                return Section._update(d1.subsections[k], v)
            elif isinstance(v, Section):
                d1.insert(v)

    def insert(self, d):
        """
        Insert a new section as a subsection of the current one
        """
        self.subsections[d.alias or d.name] = d.__deepcopy__()

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
        path = path.split("/")
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
        path = path.split("/")
        if path[0] == self.name:
            path = path[1:]
        s = self
        for p in path:
            s = s[p]
        return s

    def get_string(self):
        """
        Get string representation of Section
        """
        return Section._get_string(self)

    @staticmethod
    def _get_string(d, indent=0):
        """
        Helper function to return a pretty string of the section. Includes indentation and
        descriptions (if present).
        """
        string = ""
        if d.description and d.verbose:
            string += (
                "\n"
                + textwrap.fill(
                    d.description,
                    initial_indent="\t" * indent + "! ",
                    subsequent_indent="\t" * indent + "! ",
                    width=50,
                )
                + "\n"
            )
        string += "\t" * indent + "&" + d.name
        string += " " + " ".join(map(str, d.section_parameters)) + "\n"

        for s in d.keywords:
            string += "\t" * (indent + 1) + s.__str__() + "\n"
        for k, v in d.subsections.items():
            string += v._get_string(v, indent + 1)
        string += "\t" * indent + "&END " + d.name + "\n"

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


class Keyword(MSONable):

    """
    Class representing a keyword argument in CP2K. Within CP2K Secitons, which activate features
    of the CP2K code, the keywords are arguments that control the functionality of that feature.
    For example, the section "FORCE_EVAL" activates the evaluation of forces/energies, but within
    "FORCE_EVAL" the keyword "METHOD" controls whether or not this will be done with, say,
    "Quickstep" (DFT) or "EIP" (empirical interatomic potential).
    """

    def __init__(
        self,
        name,
        *values,
        description=None,
        units=None,
        verbose=True,
        repeats=False,
    ):
        """
        Initializes a keyword. These Keywords and the value passed to them are sometimes as simple as KEYWORD VALUE,
        but can also be more elaborate such as KEYWORD [UNITS] VALUE1 VALUE2, which is why this class exists:
        to handle many values and control easy printing to an input file.

        Args:
            name (str): The name of this keyword. Must match an acceptable keyword from CP2K
            args: All non-keyword arguments after 'name' are interpreted as the values to set for
                this keyword. i.e: KEYWORD ARG1 ARG2 would provide two values to the keyword.
            description (str): The description for this keyword. This can make readability of
                input files easier for some. Default=None.
            units (str): The units for this keyword. If not specified, CP2K default units will be
                used. Consult manual for default units. Default=None.
            repeats (bool): Whether or not this keyword may be repeated. Default=False.
        """
        self.name = name
        self.values = values
        self.description = description
        self.repeats = repeats
        self.units = units
        self.verbose = verbose

    def __str__(self):
        return (
            self.name.__str__()
            + " "
            + ("[{}] ".format(self.units) if self.units else "")
            + " ".join(map(str, self.values))
            + (
                " ! " + self.description
                if (self.description and self.verbose)
                else ""
            )
        )

    def __eq__(self, other):
        if self.name == other.name:
            if self.values == other.values:
                if self.units == self.units:
                    return True
        return False

    def as_dict(self):
        """
        Get a dictionary representation of the Keyword
        """
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["name"] = self.name
        d["values"] = self.values
        d["description"] = self.description
        d["repeats"] = self.repeats
        d["units"] = self.units
        d["verbose"] = self.verbose
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Initialise from dictonary
        """
        return Keyword(
            d["name"],
            *d["values"],
            description=d["description"],
            repeats=d["repeats"],
            units=d["units"],
            verbose=d["verbose"]
        )

    @classmethod
    def from_string(self, s):
        """
        Initialise from a string
        """
        return Keyword(*s.split())

    def silence(self):
        """
        Deactivate the printing of this keyword's description.
        """
        self.verbose = False


class Cp2kInput(Section):

    """
    Special instance of 'Section' class that is meant to represent the overall cp2k input.
    Distinguishes itself from Section by overriding get_string() to not print this section's
    title and by implementing the file i/o
    """

    def __init__(self, name="CP2K_INPUT", subsections={}, **kwargs):
        """
        Initialize Cp2kInput by calling the super
        """
        self.name = name
        self.subsections = subsections
        self.kwargs = kwargs

        description = "CP2K Input"
        super(Cp2kInput, self).__init__(
            name,
            repeats=False,
            description=description,
            keywords=[],
            section_parameters=[],
            subsections=subsections,
            **kwargs
        )

    def get_string(self):
        """
        Get string representation of the Cp2kInput
        """
        s = ""
        for k in self.subsections.keys():
            s += self.subsections[k].get_string()
        return s

    @classmethod
    def _from_dict(cls, d):
        """
        Initialize from a dictionary
        """

        return Cp2kInput(
            "CP2K_INPUT", subsections=getattr(
                __import__(d['@module'], globals(),
                           locals(), d['@class'], 0),
                d['@class']
                ).from_dict(d).subsections
        )

    @staticmethod
    def from_file(file):
        """
        Initialize from a file
        """
        with zopen(file, "rt") as f:
            lines = f.read().splitlines()
            lines = [line.replace("\t", "") for line in lines]
            lines = [line for line in lines if line]
            return Cp2kInput.from_lines(lines)

    @classmethod
    def from_lines(cls, lines):
        """
        Helper method to read lines of file
        """
        cp2k_input = Cp2kInput("CP2K_INPUT", subsections={})
        Cp2kInput._from_lines(cp2k_input, lines)
        return cp2k_input

    def _from_lines(self, lines):
        """
        Helper method, reads lines of text to get a Cp2kInput
        """
        current = self.name
        for line in lines:
            if line.startswith("!") or line.startswith("#"):
                continue
            elif line.startswith("&END"):
                current = "/".join(current.split("/")[:-1])
            elif line.startswith("&"):
                name, subsection_params = line.split()[0][1:], line.split()[1:]
                s = Section(
                    name, section_parameters=subsection_params, subsections={}
                )
                self.by_path(current).insert(s)
                current = current + "/" + name
            else:
                args = map(_postprocessor, line.split())
                self.by_path(current).keywords.append(Keyword(*args))

    def write_file(
        self,
        input_filename="cp2k.inp",
        output_dir=".",
        make_dir_if_not_present=True,
    ):
        """
        Write input to a file.
        """
        if not os.path.isdir(output_dir) and make_dir_if_not_present:
            os.mkdir(output_dir)
        filepath = os.path.join(output_dir, input_filename)
        with open(filepath, "w") as f:
            f.write(self.get_string())


class Global(Section):

    """
    Controls 'global' settings for cp2k execution such as RUN_TYPE and PROJECT_NAME
    """

    def __init__(
        self,
        project_name="CP2K",
        run_type="ENERGY_FORCE",
        **kwargs
    ):
        """
        Initialize the global section
        """

        self.project_name = project_name
        self.run_type = run_type
        self.kwargs = kwargs

        description = (
            "Section with general information regarding which kind of simulation"
            + "to perform an parameters for the whole PROGRAM"
        )

        keywords = [
            Keyword("PROJECT_NAME", project_name),
            Keyword("RUN_TYPE", run_type),
            Keyword("EXTENDED_FFT_LENGTHS", True)
        ]

        super().__init__(
            "GLOBAL",
            description=description,
            keywords=keywords,
            subsections={},
            **kwargs
        )


class ForceEval(Section):

    """
    Controls the calculation of energy and forces in Cp2k
    """

    def __init__(self, subsections={}, **kwargs):
        """
        Initialize the ForceEval section
        """

        self.subsections = subsections
        self.kwargs = kwargs

        description = "parameters needed to calculate energy and forces and describe the system you want to analyze."

        keywords = [
            Keyword("METHOD", kwargs.get("METHOD", "QS")),
            Keyword("STRESS_TENSOR", kwargs.get("STRESS_TENSOR", "ANALYTICAL")),
        ]

        super(ForceEval, self).__init__(
            "FORCE_EVAL",
            repeats=True,
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs
        )


class Dft(Section):

    """
    Controls the DFT parameters in Cp2k
    """

    def __init__(
        self,
        basis_set_filenames="BASIS_MOLOPT",
        potential_filename="GTH_POTENTIALS",
        uks=True,
        wfn_restart_file_name=None,
        subsections={},
        **kwargs
    ):
        """
        Initialize the DFT section

        Args:
            subsections: Any subsections to initialize with
            basis_set_filename: Name of the file that contains the basis set information
            potential_filename: Name of the file that contains the pseudopotential information
            uks: Whether to run unrestricted Kohn Sham (spin polarized)
        """

        self.basis_set_filenames = basis_set_filenames
        self.potential_filename = potential_filename
        self.uks = uks
        self.wfn_restart_file_name = wfn_restart_file_name
        self.subsections = subsections
        self.kwargs = kwargs

        description = "parameter needed by dft programs"

        keywords = [Keyword('BASIS_SET_FILE_NAME', k) for k in basis_set_filenames]
        keywords.append(Keyword("POTENTIAL_FILE_NAME", potential_filename))
        keywords.append(Keyword("UKS", uks))

        if wfn_restart_file_name:
            keywords.append(Keyword('WFN_RESTART_FILE_NAME', wfn_restart_file_name))

        super().__init__(
            "DFT",
            description=description,
            keywords=keywords,
            subsections=self.subsections,
            **kwargs
        )


class Subsys(Section):

    """
    Controls the definition of the system to be simulated
    """

    def __init__(self, subsections={}, **kwargs):
        """
        Initialize the subsys section
        """
        self.subsections = subsections
        self.kwargs = kwargs
        description = "a subsystem: coordinates, topology, molecules and cell"
        super(Subsys, self).__init__(
            "SUBSYS", description=description, subsections=subsections, **kwargs
        )


class QS(Section):

    """
    Controls the quickstep settings (DFT driver)
    """

    def __init__(
        self,
        method="GPW",
        eps_default=1e-7,
        extrapolation="PS",
        subsections={},
        **kwargs
    ):
        """
        Initialize the QS Section

        Args:
            method: What dft methodology to use. Can be GPW (Gaussian Plane Waves) for DFT with pseudopotentials
                or GAPW (Gaussian Augmented Plane Waves) for all electron calculations
            eps_default: The default level of convergence accuracy. NOTE: This is a global value for all the numerical
                value of all EPS_* values in QS module. It is not the same as EPS_SCF, which sets convergence accuracy
                of the SCF cycle alone.
            extrapolation: Method use for extrapolation. If using gamma-point-only calculation, then one should use
                PS for relaxations and ASPC for MD. See the manual for other options.
            subsections: Subsections to initialize with
        """

        self.method = method
        self.eps_default = eps_default
        self.extrapolation = extrapolation
        self.subsections = subsections
        self.kwargs = kwargs

        description = "parameters needed to set up the Quickstep framework"

        keywords = [
            Keyword("METHOD", method),
            Keyword("EPS_DEFAULT", eps_default),
            Keyword("EXTRAPOLATION", extrapolation),
        ]

        super(QS, self).__init__(
            "QS",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs
        )


class Scf(Section):

    """
    Controls the self consistent field loop
    """

    def __init__(
        self,
        max_scf=50,
        eps_scf=1e-6,
        scf_guess="RESTART",
        subsections={},
        **kwargs
    ):
        """
        Initialize the Scf section

        Args:
            max_scf: maximum number of SCF loops before terminating
            eps_scf: convergence criteria for SCF loop
            scf_guess: Initial guess for SCF loop (RESTART will switch to ATOMIC when no restart file
                is present)
        """

        self.max_scf = max_scf
        self.eps_scf = eps_scf
        self.scf_guess = scf_guess
        self.subsections = subsections
        self.kwargs = kwargs

        description = "Parameters needed to perform an SCF run."

        keywords = [
            Keyword("MAX_SCF", max_scf),
            Keyword("EPS_SCF", eps_scf),
            Keyword(
                "SCF_GUESS", scf_guess
            ),  # Uses Restart file if present, and ATOMIC if not present
            Keyword("MAX_ITER_LUMO", kwargs.get('max_iter_lumo', 400))
        ]

        super(Scf, self).__init__(
            "SCF",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs
        )


class Mgrid(Section):

    """
    Controls the multigrid for numerical integration
    """

    def __init__(
        self,
        cutoff=1200,
        rel_cutoff=80,
        ngrids=5,
        progression_factor=3,
        subsections={},
        **kwargs
    ):
        """
        Initialize the MGRID section

        Args:
            cutoff: Cutoff energy (in Rydbergs for historical reasons) defining how find of Gaussians will be used
            rel_cutoff: The relative cutoff energy, which defines how to map the Gaussians onto the multigrid. If the
                the value is too low then, even if you have a high cutoff with sharp Gaussians, they will be mapped
                to the course part of the multigrid
            ngrids: number of grids to use
            progression_factor: divisor that decides how to map Gaussians the multigrid after the highest mapping is
                decided by rel_cutoff
        """

        self.cutoff = cutoff
        self.rel_cutoff = rel_cutoff
        self.ngrids = ngrids
        self.progression_factor = progression_factor
        self.subsections = subsections
        self.kwargs = kwargs

        description = (
            "Multigrid information. Multigrid allows for sharp gaussians and diffuse "
            + "gaussians to be treated on different grids, where the spacing of FFT integration "
            + "points can be tailored to the degree of sharpness/diffusiveness of the gaussians."
        )

        keywords = [
            Keyword(
                "CUTOFF",
                cutoff,
                description="Cutoff in [Ry] for finest level of the MG.",
            ),
            Keyword(
                "REL_CUTOFF",
                rel_cutoff,
                description="Controls which gaussians are mapped to which level of the MG",
            ),
            Keyword("NGRIDS", ngrids),
            Keyword("PROGRESSION_FACTOR", progression_factor),
        ]

        super(Mgrid, self).__init__(
            "MGRID",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs
        )


class Diagonalization(Section):

    """
    Controls diagonalization settings (if using traditional diagonalization).
    """

    def __init__(
        self,
        eps_adapt=0,
        eps_iter=1e-8,
        eps_jacobi=0,
        jacobi_threshold=1e-7,
        subsections={},
        **kwargs
    ):
        """
        Initialize the diagronalization section
        """

        self.eps_adapt = eps_adapt
        self.eps_iter = eps_iter
        self.eps_jacobi = eps_jacobi
        self.jacobi_threshold = jacobi_threshold
        self.subsections = subsections
        self.kwargs = kwargs

        location = "CP2K_INPUT/FORCE_EVAL/DFT/SCF/DIAGONALIZATION"

        keywords = [
            Keyword("EPS_ADAPT", eps_adapt),
            Keyword("EPS_ITER", eps_iter),
            Keyword("EPS_JACOBI", eps_jacobi),
            Keyword("JACOBI_THRESHOLD", jacobi_threshold),
        ]

        super(Diagonalization, self).__init__(
            "DIAGONALIZATION",
            keywords=keywords,
            repeats=False,
            location=location,
            subsections=self.subsections,
            **kwargs
        )


class OrbitalTransformation(Section):

    """
    Turns on the Orbital Transformation scheme for diagonalizing the Hamiltonian. Much faster and with
    guaranteed convergence compared to normal diagonalization, but requires the system to have a band
    gap.

    NOTE: OT has poor convergence for metallic systems and cannot use SCF mixing or smearing. Therefore,
    you should not use it for metals or systems with 'small' band gaps. In that case, use normal
    diagonalization, which will be slower, but will converge properly.
    """

    def __init__(
        self,
        minimizer="CG",
        preconditioner="FULL_ALL",
        algorithm="STRICT",
        energy_gap=0.01,
        linesearch="2PNT",
        subsections={},
        **kwargs
    ):
        """
        Initialize the OT section

        Args:
            minimizer (str): The minimizer to use with the OT method. Default is conjugate gradient method,
                which is more robust, but more well-behaved systems should use DIIS, which can be as much
                as 50% faster.
            preconditioner (str): Preconditionar to use for OT, FULL_ALL tends to be most robust, but is
                not always most efficient. For difficult systems, FULL_SINGLE_INVERSE can be more robust,
                and is reasonably efficient with large systems. For huge, but well behaved, systems,
                where construction of the preconditioner can take a very long time, FULL_KINETIC can be a good
                choice.
            energy_gap (float): Guess for the band gap. For FULL_ALL, should be smaller than the actual band gap,
                so simply using 0.01 is a robust value. Choosing a larger value will help if you start with a bad
                initial guess though. For FULL_SINGLE_INVERSE, energy_gap is treated as a lower bound. Values lower
                than 0.05 in this case can lead to stability issues.
            algorithm (str): What algorithm to use for OT. 'Strict': Taylor or diagonalization based algorithm.
                IRAC: Orbital Transformation based Iterative Refinement of the Approximative Congruence
                transformation (OT/IR).
            linesearch (str): From the manual: 1D line search algorithm to be used with the OT minimizer,
                in increasing order of robustness and cost. MINIMIZER CG combined with LINESEARCH
                GOLD should always find an electronic minimum. Whereas the 2PNT minimizer is almost always OK,
                3PNT might be needed for systems in which successive OT CG steps do not decrease the total energy.
        """

        self.minimizer = minimizer
        self.preconditioner = preconditioner
        self.algorithm = algorithm
        self.energy_gap = energy_gap
        self.linesearch = linesearch
        self.subsections = subsections
        self.kwargs = kwargs

        description = (
            "Sets the various options for the orbital transformation (OT) method. "
            + "Default settings already provide an efficient, yet robust method. Most "
            + "systems benefit from using the FULL_ALL preconditioner combined with a small "
            + "value (0.001) of ENERGY_GAP. Well-behaved systems might benefit from using "
            + "a DIIS minimizer. Advantages: It's fast, because no expensive diagonalisation"
            + "is performed. If preconditioned correctly, method guaranteed to find minimum. "
            + "Disadvantages: Sensitive to preconditioning. A good preconditioner can be "
            + "expensive. No smearing, or advanced SCF mixing possible: POOR convergence for "
            + "metalic systems."
        )

        keywords = [
            Keyword("MINIMIZER", minimizer),
            Keyword("PRECONDITIONER", preconditioner),
            Keyword("ENERGY_GAP", energy_gap),
            Keyword("ALGORITHM", algorithm),
            Keyword("LINESEARCH", linesearch),
        ]

        super(OrbitalTransformation, self).__init__(
            "OT",
            description=description,
            keywords=keywords,
            subsections=self.subsections,
            **kwargs
        )


class Cell(Section):

    """
    Defines the simulation cell (lattice)
    """

    def __init__(self, lattice, **kwargs):
        """
        Initialize the cell section.

        Args:
            lattice: pymatgen lattice object
        """

        self.lattice = lattice
        self.kwargs = kwargs

        description = "Input parameters needed to set up the CELL."

        keywords = [
            Keyword("A", *lattice.matrix[0]),
            Keyword("B", *lattice.matrix[1]),
            Keyword("C", *lattice.matrix[2]),
        ]

        super(Cell, self).__init__(
            "CELL",
            description=description,
            keywords=keywords,
            subsections={},
            **kwargs
        )


class Kind(Section):

    """
    Specifies the information for the different atom types being simulated.
    """

    def __init__(
        self,
        specie,
        alias=None,
        magnetization=0.0,
        subsections={},
        basis_set="GTH_BASIS",
        potential="GTH_POTENTIALS",
        ghost=False,
        **kwargs
    ):
        """
        Initialize a KIND section

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

        self.name = "KIND"
        self.specie = specie
        self.alias = alias
        self.magnetization = magnetization
        self.subsections = subsections
        self.basis_set = basis_set
        self.potential = potential
        self.ghost = ghost
        self.kwargs = kwargs

        self.description = "The description of the kind of the atoms (mostly for QM)"

        keywords = [
            Keyword("ELEMENT", specie.__str__()),
            Keyword("MAGNETIZATION", magnetization),
            Keyword("BASIS_SET", basis_set),
            Keyword("POTENTIAL", potential),
        ]
        if ghost:
            keywords.append(Keyword('GHOST', ghost))

        kind_name = alias if alias else specie.__str__()
        self.alias = kind_name

        self.section_parameters = [kind_name]
        self.location = None
        self.verbose = True
        self.repeats = False

        super(Kind, self).__init__(
            name=self.name,
            subsections=self.subsections,
            description=self.description,
            keywords=keywords,
            section_parameters=self.section_parameters,
            alias=self.alias,
            location=self.location,
            verbose=self.verbose,
            **self.kwargs
        )


class Coord(Section):

    """
    Specifies the coordinates of the atoms using a pymatgen structure object.
    """

    def __init__(self, structure, alias=False, subsections={}, **kwargs):
        """
        Args:
            structure: Pymatgen structure object
            alias (bool): whether or not to identify the sites by Element + number so you can do things like
                assign unique magnetization do different elements.
        """

        self.structure = structure
        self.alias = alias
        self.subsections = subsections
        self.kwargs = kwargs

        description = (
            "The coordinates for simple systems (like small QM cells) are specified "
            + "here by default using explicit XYZ coordinates. More complex systems "
            + "should be given via an external coordinate file in the SUBSYS%TOPOLOGY section."
        )
        if alias:
            keywords = []
            for k, v in alias.items():
                keywords.extend([Keyword(k, *structure[i].coords) for i in v])
        else:
            keywords = [
                Keyword(s.specie.symbol, *s.coords) for s in structure.sites
            ]

        super(Coord, self).__init__(
            name="COORD", description=description, keywords=keywords,
            subsections=self.subsections, **kwargs
        )


class PDOS(Section):

    """
    Controls printing of projected density of states onto the different atom KINDS
    (elemental decomposed DOS).
    """

    def __init__(self, nlumo=-1, **kwargs):
        """
        Initialize the PDOS section

        Args:
            nlumo: how many unoccupied orbitals to include (-1==ALL)
        """

        self.nlumo = nlumo
        self.kwargs = kwargs

        description = "Controls printing of the projected density of states"

        keywords = [Keyword("NLUMO", nlumo), Keyword("COMPONENTS")]

        super(PDOS, self).__init__(
            "PDOS", description=description, keywords=keywords,
            subsections={}, **kwargs
        )


class LDOS(Section):

    """
    Controls printing of the LDOS (List-Density of states). i.e. projects onto specific atoms.
    """

    def __init__(self, index, alias=None, **kwargs):
        """
        Initialize the LDOS section

        Args:
            index: Index of the atom to project onto
        """

        self.index = index
        self.alias = alias
        self.kwargs = kwargs

        description = "Controls printing of the projected density of states decomposed by atom type"

        keywords = [Keyword("COMPONENTS"), Keyword("LIST", index)]

        super(LDOS, self).__init__(
            "LDOS",
            subsections={},
            alias=alias,
            description=description,
            keywords=keywords,
            **kwargs
        )


class V_Hartree_Cube(Section):

    """
    Controls printing of the hartree potential as a cube file.
    """

    def __init__(self, **kwargs):
        """
        Initialize the V_HARTREE_CUBE section
        """

        self.kwargs = kwargs

        description = (
            "Controls the printing of a cube file with eletrostatic potential generated by "
            + "the total density (electrons+ions). It is valid only for QS with GPW formalism. "
            + "Note that by convention the potential has opposite sign than the expected physical one."
        )

        keywords = [
            Keyword('stride', *kwargs.get('STRIDE', [1, 1, 1]))
        ]

        super(V_Hartree_Cube, self).__init__(
            "V_HARTREE_CUBE",
            subsections={},
            description=description,
            keywords=keywords,
            **kwargs
        )


class MO_Cubes(Section):

    """
    Controls printing of the molecular orbital eigenvalues
    """

    def __init__(self, write_cube=False, nhomo=1, nlumo=1, **kwargs):
        """
            Initialize the MO_CUBES section
        """

        self.write_cube = write_cube
        self.nhomo = nhomo
        self.nlumo = nlumo
        self.kwargs = kwargs

        description = (
            "Controls the printing of a cube file with eletrostatic potential generated by "
            + "the total density (electrons+ions). It is valid only for QS with GPW formalism. "
            + "Note that by convention the potential has opposite sign than the expected physical one."
        )

        keywords = [
            Keyword("WRITE_CUBE", write_cube),
            Keyword("NHOMO", nhomo),
            Keyword("NLUMO", nlumo),
        ]

        super(MO_Cubes, self).__init__(
            "MO_CUBES",
            subsections={},
            description=description,
            keywords=keywords,
            **kwargs
        )


class E_Density_Cube(Section):

    """
    Controls printing of the electron density cube file
    """

    def __init__(self, **kwargs):
        """
        Initialize the E_DENSITY_CUBE Section
        """

        self.kwargs = kwargs

        description = (
            "Controls the printing of cube files with the electronic density and, for LSD "
            + "calculations, the spin density."
        )

        keywords = []

        super(E_Density_Cube, self).__init__(
            "E_DENSITY_CUBE",
            subsections={},
            description=description,
            keywords=keywords,
            **kwargs
        )


class Smear(Section):

    """
    Control electron smearing
    """

    def __init__(
        self, elec_temp=300, method="FERMI_DIRAC", fixed_magnetic_moment=-1e2, **kwargs
    ):
        """
        Initialize the SMEAR section
        """

        self.elec_temp = elec_temp
        self.method = method
        self.fixed_magnetic_moment = fixed_magnetic_moment
        self.kwargs = kwargs

        description = "Activates smearing of electron occupations"

        keywords = [
            Keyword("ELEC_TEMP", elec_temp),
            Keyword("METHOD", method),
            Keyword("FIXED_MAGNETIC_MOMENT", fixed_magnetic_moment),
        ]

        super(Smear, self).__init__(
            "SMEAR", description=description, keywords=keywords,
            subsections={}, **kwargs
        )


class BrokenSymmetry(Section):

    """
    Define the required atomic orbital occupation assigned in initialization
    of the density matrix, by adding or subtracting electrons from specific
    angular momentum channels. It works only with GUESS ATOMIC
    """

    def __init__(
        self,
        l_alpha=-1,
        n_alpha=0,
        nel_alpha=-1,
        l_beta=-1,
        n_beta=0,
        nel_beta=-1,
    ):
        """
        Initialize the broken symmetry section

        Args:
            l_alpha: Angular momentum quantum number of the orbitals whose occupation is changed
            n_alpha: Principal quantum number of the orbitals whose occupation is changed.
                Default is the first not occupied
            nel_alpha: Orbital occupation change per angular momentum quantum number. In
                unrestricted calculations applied to spin alpha
            l_beta: Same as L_alpha for beta channel
            n_beta: Same as N_alpha for beta channel
            nel_beta: Same as NEL_alpha for beta channel
        """

        self.l_alpha = l_alpha
        self.n_alpha = n_alpha
        self.nel_alpha = nel_alpha
        self.l_beta = l_beta
        self.n_beta = n_beta
        self.nel_beta = nel_beta

        description = (
            "Define the required atomic orbital occupation assigned in initialization "
            + "of the density matrix, by adding or subtracting electrons from specific "
            + "angular momentum channels. It works only with GUESS ATOMIC"
        )

        keywords_alpha = [
            Keyword("L", l_alpha),
            Keyword("N", n_alpha),
            Keyword("NEL", nel_alpha),
        ]
        alpha = Section(
            "ALPHA", keywords=keywords_alpha, subsections={}, repeats=False
        )

        keywords_beta = [
            Keyword("L", l_beta),
            Keyword("N", n_beta),
            Keyword("NEL", nel_beta),
        ]
        beta = Section(
            "BETA", keywords=keywords_beta, subsections={}, repeats=False
        )

        super(BrokenSymmetry, self).__init__(
            "BS",
            description=description,
            subsections={"ALPHA": alpha, "BETA": beta},
            keywords=[],
            repeats=False,
        )


class XC_FUNCTIONAL(Section):

    """
    Defines the XC functional to use
    """

    def __init__(self, functional, subsections={}, **kwargs):
        """
        Initialize the XC_FUNCTIONAL class
        """

        self.functional = functional
        self.subsections = subsections
        self.kwargs = kwargs

        location = "CP2K_INPUT/FORCE_EVAL/DFT/XC/XC_FUNCTIONAL"

        built_in = [
            "BL3YLP",
            "BEEFVDW",
            "BLYP",
            "BP",
            "HCTH120",
            "LDA",
            "NONE",
            "NO_SHORTCUT",
            "OLYP",
            "PADE",
            "PBE",
            "PBE0",
            "TPSS",
        ]

        if functional not in built_in:
            raise Warning("The selected functional does not exist in CP2K.")

        super(XC_FUNCTIONAL, self).__init__(
            "XC_FUNCTIONAL",
            keywords=[],
            subsections=self.subsections,
            location=location,
            repeats=False,
            section_parameters=[functional],
            **kwargs
        )


class PBE(Section):

    """
        Info about the PBE functional.
    """

    def __init__(self, parameterization="ORIG", scale_c=1, scale_x=1):
        """
        Args:
            parameterization (str):
                ORIG: original PBE
                PBESOL: PBE for solids/surfaces
                REVPBE: revised PBE
            scale_c (float): scales the correlation part of the functional.
            scale_x (float): scales the exchange part of the functional.
        """

        self.parameterization = parameterization
        self.scale_c = scale_c
        self.scale_x = scale_x

        location = "CP2K_INPUT/FORCE_EVAL/DFT/XC/XC_FUNCTIONAL/PBE"

        keywords = [
            Keyword("PARAMETRIZATION", parameterization),
            Keyword("SCALE_C", scale_c),
            Keyword("SCALE_X", scale_x),
        ]

        super(PBE, self).__init__(
            "PBE",
            subsections={},
            repeats=False,
            location=location,
            section_parameters=[],
            keywords=keywords,
        )


class Kpoints(Section):

    """
    Description of the k-points to use for the calculation.
    """

    def __init__(self, kpts, weights=None, eps_geo=1e-6, full_grid=False,
                 parallel_group_size=-1,
                 scheme='MONKHORST-PACK', symmetry=False,
                 units='B_VECTOR', verbose=False,
                 wavefunctions='COMPLEX'):
        """
        Args:
            kpts (list, tuple): a 2D array for the kpoints of the form
                [(1,1,1),]. If len(kpts) == 1. Then it is taken as subdivisions
                for automatic kpoint scheme. If it has more entries, it is
                taken as manual entries for kpoints.
            eps_geo (float): tolerance for symmetry. Default=1e-6
            full_grid (bool): use full (not reduced) kpoint grid. Default=False.
            parallel_group_size (int): from cp2k manual: Number of processors
                to be used for a single kpoint. This number must divide the
                total number of processes. The number of groups must divide
                the total number of kpoints. Value=-1 (smallest possible
                number of processes per group, satisfying the constraints).
                Value=0 (all processes). Value=n (exactly n processes).
                Default=-1.
            scheme (str): kpoint generation scheme. Default='Monkhorst-Pack'
            symmetry (bool): Use symmetry to reduce the number of kpoints.
                Default=False.
            units (str): Units for the kpoint coordinates (reciprocal coordinates
                or cartesian). Default='B_VECTOR' (reciprocal)
            verbose (bool): verbose output for kpoints. Default=False
            wavefunctions (str): Whether to use complex or real valued wavefunctions
                (if available). Default='complex'
        """
        description = 'Sets up the kpoints'
        keywords = []

        self.kpts = kpts
        self.weights = weights if weights else [1]*len(kpts)
        self.eps_geo = eps_geo
        self.full_grid = full_grid
        self.parallel_group_size = parallel_group_size
        self.scheme = scheme
        self.symmetry = symmetry
        self.units = units
        self.verbose = verbose
        self.wavefunctions = wavefunctions

        if len(kpts) == 1:
            keywords.append(
                Keyword('SCHEME', scheme, *kpts[0])
            )
        elif len(kpts) > 1:
            for k, w in zip(kpts, weights):
                keywords.append(
                    Keyword('KPOINT', *k, w)
                )
        else:
            raise ValueError("No k-points provided!")

        keywords.extend(
            [
                Keyword('EPS_GEO', eps_geo),
                Keyword('FULL_GRID', full_grid),
                Keyword('PARALLEL_GROUP_SIZE', parallel_group_size),
                Keyword('SYMMETRY', symmetry),
                Keyword('UNITS', units),
                Keyword('VERBOSE', verbose),
                Keyword('WAVEFUNCTIONS', wavefunctions)
            ]
        )

        super(Kpoints, self).__init__(name='KPOINTS', subsections=None, repeats=False,
                                      description=description, keywords=keywords)

    @classmethod
    def from_kpoints(cls, kpoints):
        """
        Initialize the section from a Kpoints object (pymatgen.io.vasp.inputs).

        Args:
            kpoints: A pymatgen kpoints object.
        """
        k = kpoints.as_dict()
        kpoints = k['kpoints']
        weights = k['kpts_weights']
        scheme = k['generation_style']
        if scheme.lower() == 'Monkhorst':
            scheme = 'MONKHORST-PACK'
        units = k['coord_type']
        if k['coord_type']:
            if k['coord_type'].lower() == 'reciprocal':
                units = 'B_VECTOR'
            elif k['coord_type'].lower() == 'cartesian':
                units = 'CART_ANGSTROM'
        else:
            units = 'B_VECTOR'

        return Kpoints(kpts=kpoints, weights=weights, scheme=scheme, units=units)
