"""
This module defines input sets for CP2K and is a work in progress. The structure/philosophy
of this module is based on the Vasp input sets in Pymatgen. These sets are meant to contain
tested parameters that will result in successful, reproducible, consistent calculations without
need for intervention 99% of the time. 99% of the time, you only need to provide a pymatgen
structure object and let the defaults take over from there.

The sets are intended to be very general, e.g. a set for geometry relaxation, and so most of the
time, if you have specific needs, you can simply specify them via the keyword argument
override_default_params (see Section.update() method). If you have the need to create a new input
set (say for a standardized high throughput calculation) then you can create a new child of the
Cp2kInputSet class.

In order to implement a new Set within the current code structure, follow this 3 step flow:
    (1) Inherit from Cp2kInputSet or one of its children and call the super() constructor
    (2) Create the new sections and insert them into self and its subsections as needed
    (3) Call self.update(override_default_params) in order to allow user settings.
"""

from pathlib import Path
from pymatgen.io.cp2k.inputs import (
    Cp2kInput,
    Section,
    Keyword,
    Global,
    ForceEval,
    Mgrid,
    MO_Cubes,
    OrbitalTransformation,
    XC_FUNCTIONAL,
    V_Hartree_Cube,
    Dft,
    E_Density_Cube,
    LDOS,
    PDOS,
    Coord,
    Kind,
    QS,
    PBE,
    Cell,
    Subsys,
    Scf,
    Kpoints,
    Smear
)
from pymatgen.io.cp2k.utils import (
    get_basis_and_potential,
    get_aux_basis,
    get_unique_site_indices,
)

__author__ = "Nicholas Winner"
__version__ = "0.2"
__email__ = "nwinner@berkeley.edu"
__date__ = "January 2019"

MODULE_DIR = Path(__file__).resolve().parent


class Cp2kInputSet(Cp2kInput):

    """
    The basic representation of a CP2K input set as a collection of "sections" defining the simulation
    connected to a structure object. At the most basis level, CP2K requires a &GLOBAL section and
    &FORCE_EVAL section. Global sets parameters like "RUN_TYPE" or the overall verbosity. FORCE_EVAL is
    the largest section usually, containing the cell and coordinates of atoms, the DFT settings, and more.
    This top level input set is meant to initialize GLOBAL and FORCE_EVAL based on a structure object and
    and sections that the user provides.

    Like everything that goes into a cp2k input file, this base input set is essentially a section object.
    These sets are distinguished by saving default settings for easy implementation of calculations such
    as relaxation and static calculations. This base set is here to transfer a pymatgen structure object
    into the input format for cp2k and associate the basis set and pseudopotential to use with each
    element in the structure.

    Generally, this class will not be used directly, and instead one of
    its child-classes will be used, which contain more predefined initializations of various sections, and,
    if modifications are required, the user can specify override_default_settings.
    """

    def __init__(
        self,
        structure,
        potential_and_basis={},
        multiplicity=0,
        project_name="CP2K",
        override_default_params={},
        **kwargs
    ):
        """
        Args:
            structure: (Structure) pymatgen structure object used to define the lattice,
                        coordinates, and elements

            sections: (dict) the sections that define the simulation parameters.
                Most sections are optional, but two are required: (1) GLOBAL and (2) FORCE_EVAL.
                GLOBAL does not, itself, have required values, but the section must exist. FORCE_EVAL,
                technically, only needs SUBSYS to exist, which is written so long as a structure argument
                passed, but the code won't do anything. Sections that are defined by the structure object
                automatically are:

                    -SUBSYS Requires a list of "KIND" section objects, which describes the elements in
                    the calculation.
                    -CELL defines the lattice vectors as A, B, C
                    -COORDS defines the coordinates of the atoms within the CELL in cartesian coordinates

            potential_and_basis: (dict) Specifies what basis set and potential to use. Specify each as "basis" and
                "potential". Note: These must be consistent with one another, e.g. a GTH basis
                should go with a GTH potential. They potentials of all the elements should also
                match, e.g. if you use a double zeta potential for one element, you should
                use it for all.

                    Ex: {'potential': 'GTH', 'functional': 'PBE', 'basis': 'DZVP'} for a
                        GTH type potential with a double-zeta valence GTH basis and PBE
                        XC functional.

            override_default_params: (dict) Specifies user-defined settings to override the settings of any
                input set (See Section.update())
        """
        super(Cp2kInputSet, self).__init__(name="CP2K_INPUT", subsections={})

        # Important CP2K set parameters
        self.structure = structure
        self.charge = structure.charge
        self.potential_and_basis = potential_and_basis
        self.multiplicity = multiplicity  # spin multiplicity = 2s+1
        self.override_default_params = override_default_params
        self.project_name = project_name
        self.kwargs = kwargs

        cell = Cell(structure.lattice)
        subsys = Subsys(subsections={"CELL": cell})

        # Decide what basis sets/pseudopotentials to use
        basis_and_potential = get_basis_and_potential(
            [str(s) for s in structure.species],
            potential_and_basis)
        self.basis_set_file_names = basis_and_potential["basis_filenames"]
        self.potential_file_name = basis_and_potential["potential_filename"]

        # Insert atom kinds by identifying the unique sites (unique element and site properties)
        unique_kinds = get_unique_site_indices(structure)
        for k in unique_kinds.keys():
            kind = k.split("_")[0]
            if "magmom" in self.structure.site_properties.keys():
                magnetization = self.structure.site_properties["magmom"][
                    unique_kinds[k][0]
                ]
            else:
                magnetization = 0
            if 'ghost' in self.structure.site_properties:
                ghost = self.structure.site_properties.get("ghost")[unique_kinds[k][0]]
            else:
                ghost = False
            subsys.insert(
                Kind(
                    kind,
                    alias=k,
                    basis_set=basis_and_potential[kind]["basis"],
                    magnetization=magnetization,
                    potential=basis_and_potential[kind]["potential"],
                    ghost=ghost
                )
            )
        coord = Coord(structure, aliases=unique_kinds)
        subsys.insert(coord)

        if not self.check("FORCE_EVAL"):
            self.insert(ForceEval())

        self["FORCE_EVAL"].insert(subsys)
        self.print_forces()
        self.print_motion()

        self.update(override_default_params)

    def create_structure(self, structure=None):
        """
        Create the structure for the input
        """
        pass

    def activate_hybrid(
        self,
        hybrid_functional="PBE0",
        hf_fraction=0.25,
        gga_x_fraction=0.75,
        gga_c_fraction=1,
        max_memory=2000,
        cutoff_radius=8.,
        omega=.2,
        aux_basis=None,
        **kwargs
    ):

        """
        Basic set for activating hybrid DFT calculation using Auxiliary Density Matrix Method.

        Note 1: When running ADMM with cp2k, memory is very important. If the memory requirements exceed
        what is available (see max_memory), then CP2K will have to calculate the 4-electron integrals
        for HFX during each step of the SCF cycle. ADMM provides a huge speed up by making the memory
        requirements *feasible* to fit into RAM, which means you only need to calculate the integrals
        once each SCF cycle. But, this only works if it fits into memory. When setting up ADMM
        calculations, we recommend doing whatever is possible to fit all the 4EI into memory.

        Note 2: This set is designed for reliable high-throughput calculations, NOT for extreme
        accuracy. Please review the in-line comments in this method if you want more control.

        Args:
            method (str): Type of hybrid functional. This set supports HSE (screened) and PBE0
                (truncated). Default is PBE0, which converges easier in the GPW basis used by
                cp2k.
            hf_fraction (float): fraction of exact HF exchange energy to mix. Default: 0.25
            gga_x_fraction (float): fraction of gga exchange energy to retain. Default: 0.75
            gga_c_fraction (float): fraction of gga correlation energy to retain. Default: 1.0
            max_memory (int): Maximum memory available to each MPI process (in Mb) in the calculation.
                Most modern computing nodes will have ~2Gb per core, or 2048 Mb, but check for
                your specific system. This value should be as large as possible while still leaving
                some memory for the other parts of cp2k. Important: If this value is set larger
                than the memory limits, CP2K will likely seg-fault.
                Default: 2000
            cutoff_radius (float): for truncated hybrid functional (i.e. PBE0), this is the cutoff
                radius. The default is selected as that which generally gives convergence, but
                maybe too low (if you want very high accuracy) or too high (if you want a quick
                screening). Default: 8 angstroms
        """
        if aux_basis is None:
            basis = get_aux_basis(self.structure.symbol_set)
        else:
            basis = aux_basis
        self["FORCE_EVAL"]["DFT"].keywords.extend(
            [Keyword("BASIS_SET_FILE_NAME", k) for k in ['BASIS_ADMM', 'BASIS_ADMM_MOLOPT']],
        )

        for k, v in self["FORCE_EVAL"]["SUBSYS"].subsections.items():
            if "KIND" == v.name.upper():
                kind = v.get_keyword("ELEMENT").values[0]
                v.keywords.append(
                    Keyword("BASIS_SET", "AUX_FIT", basis[kind])
                )

        aux_matrix_params = [
            Keyword("ADMM_PURIFICATION_METHOD", "NONE"),  # Use NONE for accurate eigenvalues (static calcs)
            Keyword("METHOD", "BASIS_PROJECTION"),  # Don't change unless you know what you're doing
        ]
        aux_matrix = Section(
            "AUXILIARY_DENSITY_MATRIX_METHOD",
            keywords=aux_matrix_params,
            subsections={},
        )

        # Define the GGA functional as PBE
        pbe = PBE("ORIG", scale_c=gga_c_fraction, scale_x=gga_x_fraction)
        xc_functional = XC_FUNCTIONAL("PBE", subsections={"PBE": pbe})

        # Define screening of the HF term
        screening = Section(
            "SCREENING",
            subsections={},
            keywords=[
                Keyword("EPS_SCHWARZ", 1e-6),  # Aggressive screening for high throughput
                Keyword("EPS_SCHWARZ_FORCES", 1e-6),  # Aggressive screening for high throughput
                Keyword("SCREEN_ON_INITIAL_P", True),  # Better performance. Requires decent starting wfn
                Keyword("SCREEN_P_FORCES", True),  # Better performance. Requires decent starting wfn
            ],
        )

        if hybrid_functional == "HSE06":
            xc_functional.insert(
                Section(
                    "XWPBE",
                    subsections={},
                    keywords=[
                        Keyword("SCALE_X0", 1),
                        Keyword("SCALE_X", -hf_fraction),
                        Keyword("OMEGA", omega),
                    ],
                )
            )
            ip_keywords = [
                Keyword("POTENTIAL_TYPE", "SHORTRANGE"),
                Keyword(
                    "OMEGA", omega  # HSE06 screening parameter
                ),
            ]
        elif hybrid_functional == "PBE0":
            ip_keywords = [
                Keyword("POTENTIAL_TYPE", "TRUNCATED"),
                Keyword(
                    "CUTOFF_RADIUS", cutoff_radius
                ),
                Keyword("T_C_G_DATA", "t_c_g.dat"),
            ]
        interaction_potential = Section(
            "INTERACTION_POTENTIAL", subsections={}, keywords=ip_keywords
        )

        load_balance = Section(
            "LOAD_BALANCE",
            keywords=[Keyword("RANDOMIZE", True)],
            subsections={},
        )
        memory = Section(
            "MEMORY",
            subsections={},
            keywords=[
                Keyword("EPS_STORAGE_SCALING", 0.1),  # squashes the integrals for efficient storage
                Keyword("MAX_MEMORY", max_memory),
            ],
        )
        hf = Section(
            "HF",
            keywords=[Keyword("FRACTION", hf_fraction)],
            subsections={
                "SCREENING": screening,
                "INTERACTION_POTENTIAL": interaction_potential,
                "LOAD_BALANCE": load_balance,
                "MEMORY": memory,
            },
        )
        xc = Section(
            "XC", subsections={"XC_FUNCTIONAL": xc_functional, "HF": hf}
        )

        self.subsections["FORCE_EVAL"]["DFT"].insert(aux_matrix)
        self.subsections["FORCE_EVAL"]["DFT"].insert(xc)

    def print_forces(self):
        """
        Print out the forces and stress during calculation
        """
        self["FORCE_EVAL"].insert(Section("PRINT", subsections={}))
        self["FORCE_EVAL"]["PRINT"].insert(Section("FORCES", subsections={}))
        self["FORCE_EVAL"]["PRINT"].insert(
            Section("STRESS_TENSOR", subsections={})
        )

    def print_motion(self):
        """
        Print the motion info (trajectory, cell, forces, stress
        """
        if not self.check("MOTION"):
            self.insert(Section("MOTION", subsections={}))
        self["MOTION"].insert(Section("PRINT", subsections={}))
        self["MOTION"]["PRINT"].insert(
            Section("TRAJECTORY", section_parameters=["ON"], subsections={})
        )
        self["MOTION"]["PRINT"].insert(Section("CELL", subsections={}))
        self["MOTION"]["PRINT"].insert(Section("FORCES", subsections={}))
        self["MOTION"]["PRINT"].insert(Section("STRESS", subsections={}))


class DftSet(Cp2kInputSet):
    """
    Base for an input set using the Quickstep module (i.e. a DFT calculation). The DFT section is pretty vast
    in CP2K, so this set hopes to make the DFT setup fairly simple. The provided parameters are pretty conservative,
    and so they should not need to be changed very often.
    """

    def __init__(
        self,
        structure,
        ot=True,
        band_gap=0.01,
        eps_default=1e-12,
        eps_scf=1e-7,
        max_scf=None,
        minimizer="DIIS",
        preconditioner="FULL_ALL",
        algorithm="STRICT",
        linesearch="2PNT",
        cutoff=1200,
        rel_cutoff=80,
        ngrids=5,
        progression_factor=3,
        override_default_params={},
        wfn_restart_file_name=None,
        kpoints=None,
        smearing=None,
        **kwargs
    ):
        """
        Args:
            structure: Pymatgen structure object
            ot (bool): Whether or not to use orbital transformation method for matrix diagonalization. OT is the
                flagship matrix diagonalizer of CP2K, and will provide huge speed-ups for this part of the calculation,
                but the system must have a band gap for OT to be used (higher band-gap --> faster convergence).
                Band gap is also used by the preconditioner for OT, and should be set as a value SMALLER than the true
                band gap to get good efficiency. Generally, this parameter does not need to be changed from
                default of 0.01
            band_gap (float): The band gap can also be specified in order to determine if ot should be turned on.
            eps_default (float): Replaces all EPS_XX Keywords in the DFT section (NOT its subsections!) to have this
                value, ensuring an overall accuracy of at least this much.
            minimizer (str): The minimization scheme. DIIS can be as much as 50% faster than the more robust conjugate
                gradient method, and so it is chosen as default. Switch to CG if dealing with a difficult system.
            preconditioner (str): Preconditioner for the OT method. The FULL_ALL preconditioner with a small estimate
                of the band gap is usually very robust. Should only change when simulation cell gets to be very large.
            cutoff (int): Cutoff energy (in Ry) for the finest level of the multigrid. A high cutoff will allow you to
                have very accurate calculations PROVIDED that REL_CUTOFF is appropriate.
            rel_cutoff (int): This cutoff decides how the Guassians are mapped onto the different levels of the
                multigrid. From CP2K: A Gaussian is mapped onto the coarsest level of the multi-grid, on which the
                    function will cover number of grid points greater than or equal to the number of grid points
                    will cover on a reference grid defined by REL_CUTOFF.
            progression_factor (int): Divisor of CUTOFF to get the cutoff for the next level of the multigrid.

            Takeaway for the cutoffs: https://www.cp2k.org/howto:converging_cutoff
            If CUTOFF is too low, then all grids will be coarse and the calculation may become inaccurate; and if
            REL_CUTOFF is too low, then even if you have a high CUTOFF, all Gaussians will be mapped onto the coarsest
            level of the multi-grid, and thus the effective integration grid for the calculation may still be too
            coarse.
        """

        super(DftSet, self).__init__(structure, **kwargs)

        self.structure = structure
        self.ot = ot
        self.band_gap = band_gap
        self.eps_default = eps_default
        self.eps_scf = eps_scf
        self.max_scf = max_scf
        self.minimizer = minimizer
        self.preconditioner = preconditioner
        self.algorithm = algorithm
        self.linesearch = linesearch
        self.cutoff = cutoff
        self.rel_cutoff = rel_cutoff
        self.ngrids = ngrids
        self.progression_factor = progression_factor
        self.override_default_params = override_default_params
        self.wfn_restart_file_name = wfn_restart_file_name
        self.kpoints = kpoints
        self.smearing = smearing
        self.kwargs = kwargs

        # Build the QS Section
        qs = QS(eps_default=eps_default)
        max_scf = max_scf if max_scf else 50 if ot else 400  # If ot, max_scf is for inner loop
        scf = Scf(eps_scf=eps_scf, max_scf=max_scf, subsections={})

        # If there's a band gap, always use OT, else use Davidson
        if ot:
            if band_gap <= 0:
                raise UserWarning('Orbital Transformation method is being used for'
                                  'a system without a bandgap. OT can have very poor'
                                  'convergence for metallic systems, proceed with caution.')
            scf.insert(
                OrbitalTransformation(
                    minimizer=minimizer,
                    preconditioner=preconditioner,
                    energy_gap=band_gap,
                    algorithm=algorithm,
                    linesearch=linesearch,
                )
            )
            scf.insert(
                Section(
                    "OUTER_SCF",
                    subsections={},
                    keywords=[
                        Keyword("MAX_SCF", kwargs.get("outer_max_scf", 20)),
                        Keyword(
                            "EPS_SCF", kwargs.get("outer_eps_scf", eps_scf)
                        ),
                    ],
                )
            )
        else:
            scf.insert(Section("DIAGONALIZATION", subsections={}))
            mixing_kwds = [
                Keyword('METHOD', 'BROYDEN_MIXING'),
                Keyword('ALPHA', 0.2),
                Keyword('NBUFFER', 5)
            ]
            mixing = Section('MIXING', keywords=mixing_kwds, subsections=None)
            scf.insert(mixing)
            davidson_kwds = [
                Keyword('PRECONDITIONER', 'FULL_SINGLE_INVERSE')
            ]
            davidson = Section('DAVIDSON', keywords=davidson_kwds, subsections=None)
            scf["DIAGONALIZATION"].insert(davidson)

        # Create the multigrid for FFTs
        mgrid = Mgrid(
            cutoff=cutoff,
            rel_cutoff=rel_cutoff,
            ngrids=ngrids,
            progression_factor=progression_factor,
        )

        # Set the DFT calculation with global parameters
        dft = Dft(
            MULTIPLICITY=self.multiplicity,
            CHARGE=self.charge,
            basis_set_filenames=self.basis_set_file_names,
            potential_filename=self.potential_file_name,
            subsections={"QS": qs, "SCF": scf, "MGRID": mgrid},
            wfn_restart_file_name=wfn_restart_file_name
        )

        if kpoints:
            dft.insert(Kpoints.from_kpoints(kpoints))
        if smearing or (band_gap <= 0.0):
            scf['ADDED_MOS'] = 100  # TODO: how to grab the appropriate number?
            scf.insert(Smear())

        # Create subsections and insert into them
        self["FORCE_EVAL"].insert(dft)
        xc_functional = XC_FUNCTIONAL(functional="PBE")
        xc = Section("XC", subsections={"XC_FUNCTIONAL": xc_functional})
        self["FORCE_EVAL"]["DFT"].insert(xc)
        self["FORCE_EVAL"]["DFT"].insert(Section("PRINT", subsections={}))

        # TODO: Ugly workaround...
        if kwargs.get('print_pdos', True):
            self.print_pdos()
        if kwargs.get('print_ldos', False):
            self.print_ldos()
        if kwargs.get('print_mo_cubes', True):
            self.print_mo_cubes()
        if kwargs.get('print_hartree_potential', False):
            self.print_hartree_potential(stride=kwargs.get('stride', [1, 1, 1]))
        if kwargs.get('print_e_density', False):
            self.print_e_density()
        if kwargs.get('activate_fast_minimization', False):
            self.activate_fast_minimization()
        if kwargs.get('activate_robust_minimization', False):
            self.activate_robust_minimization()
        if kwargs.get('activate_very_strict_minimization', False):
            self.activate_very_strict_minimization()
        self.update(override_default_params)

    def print_pdos(self, nlumo=-1):
        """
        Activate creation of the PDOS file.

        Args:
            nlumo (int): Number of virtual orbitals to be added to the MO set (-1=all).
                CAUTION: Setting this value to be higher than the number of states present may cause a Cholesky error.
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/PDOS"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(PDOS(nlumo=nlumo))

    def print_ldos(self, nlumo=-1):
        """
        Activate the printing of LDOS files, printing one for each atom kind by default

        Args:
            nlumo (int): Number of virtual orbitals to be added to the MO set (-1=all).
                CAUTION: Setting this value to be higher than the number of states present may cause a Cholesky error.
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/PDOS"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(PDOS(nlumo=nlumo))
        for i in range(self.structure.num_sites):
            self["FORCE_EVAL"]["DFT"]["PRINT"]["PDOS"].insert(
                LDOS(i+1, alias="LDOS {}".format(i+1), verbose=False)
            )

    def print_mo_cubes(self, write_cube=False, nlumo=1, nhomo=1):
        """
        Activate printing of molecular orbitals.

        Args:
            write_cube (bool): whether to write cube file for the MOs (setting false will just print levels in out file)
            nlumo (int): Controls the number of lumos that are printed and dumped as a cube (-1=all)
            nhomo (int): Controls the number of homos that are printed and dumped as a cube (-1=all)
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/MO_CUBES"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(
                MO_Cubes(write_cube=write_cube, nlumo=nlumo, nhomo=nhomo)
            )

    def print_mo(self):
        """
        Print molecular orbitals when running non-OT diagonalization
        """
        raise NotImplementedError
        pass

    def print_hartree_potential(self, stride=[1, 1, 1]):
        """
        Controls the printing of a cube file with eletrostatic potential generated by the total density
        (electrons+ions). It is valid only for QS with GPW formalism.
        Note that by convention the potential has opposite sign than the expected physical one.
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(V_Hartree_Cube(stride=stride))

    def print_e_density(self):
        """
        Controls the printing of cube files with the electronic density and, for LSD calculations, the spin density
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/E_DENSITY_CUBE"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(E_Density_Cube())

    def set_charge(self, charge):
        """
        Set the overall charge of the simulation cell
        """
        self["FORCE_EVAL"]["DFT"].set_keyword(Keyword("CHARGE", charge))

    def activate_fast_minimization(self):
        """
        Method to modify the set to use fast SCF minimization.
        """
        ot = OrbitalTransformation(
            minimizer="DIIS",
            preconditioner="FULL_ALL",
            algorithm="IRAC",
            energy_gap=0.01,
            linesearch="2PNT",
        )
        self.update({"FORCE_EVAL": {"DFT": {"SCF": {"OT": ot}}}})

    def activate_robust_minimization(self):
        """
        Method to modify the set to use more robust SCF minimization technique
        """
        ot = OrbitalTransformation(
            minimizer="CG",
            preconditioner="FULL_SINGLE_INVERSE",
            algorithm="STRICT",
            energy_gap=0.05,
            linesearch="3PNT",
        )
        self.update({"FORCE_EVAL": {"DFT": {"SCF": {"OT": ot}}}})

    def activate_very_strict_minimization(self):
        """
        Method to modify the set to use very strict SCF minimization scheme
        :return:
        """
        ot = OrbitalTransformation(
            minimizer="CG",
            preconditioner="FULL_SINGLE_INVERSE",
            algorithm="STRICT",
            energy_gap=0.05,
            linesearch="GOLD",
        )
        self.update({"FORCE_EVAL": {"DFT": {"SCF": {"OT": ot}}}})

    def activate_nonperiodic(self):
        """
        Activates a calculation with non-periodic calculations by turning of PBC and
        changing the poisson solver.
        """
        self['FORCE_EVAL']['SUBSYS']['CELL'].keywords.apped(Keyword('PERIODIC', 'NONE'))
        kwds = [Keyword('POISSON_SOLVER', 'MT'), Keyword('PERIODIC', 'NONE')]
        self['FORCE_EVAL']['DFT'].insert(Section('POISSON', subsections={}, keywords=kwds))


class StaticSet(DftSet):

    """
    Basic static energy calculation. Turns on Quickstep module, sets the run_type in global,
    and uses structure object to build the subsystem.
    """

    def __init__(
        self,
        structure,
        project_name="Static",
        run_type="ENERGY_FORCE",
        override_default_params={},
        **kwargs
    ):
        """
        Args:
            structure: Pymatgen structure object
            project_name (str): What to name this cp2k project (controls naming of files printed out)
            run_type (str): Run type. As a static set it should be one of the static aliases, like 'ENERGY_FORCE'
        """
        super(StaticSet, self).__init__(structure, **kwargs)
        global_section = Global(project_name=project_name, run_type=run_type)
        self.structure = structure
        self.project_name = project_name
        self.run_type = run_type
        self.override_default_params = override_default_params
        self.insert(global_section)
        self.update(override_default_params)
        self.kwargs = kwargs


class RelaxSet(DftSet):

    """
    CP2K input set containing the basic settings for performing geometry optimization. Values are all cp2k
    defaults, and should be good for most systems of interest.
    """

    def __init__(
        self,
        structure,
        max_drift=3e-3,
        max_force=4.5e-3,
        max_iter=200,
        project_name="Relax",
        optimizer="BFGS",
        override_default_params={},
        **kwargs
    ):

        """
        Args:
            structure:
            max_drift: Convergence criterion for the maximum geometry change between the current and the
                last optimizer iteration. This keyword cannot be repeated and it expects precisely one real.
                Default value: 3.00000000E-003
                Default unit: [bohr]
            max_force (float): Convergence criterion for the maximum force component of the current configuration.
                This keyword cannot be repeated and it expects precisely one real.
                Default value: 4.50000000E-004
                Default unit: [bohr^-1*hartree]
            max_iter (int): Specifies the maximum number of geometry optimization steps.
                One step might imply several force evaluations for the CG and LBFGS optimizers.
                This keyword cannot be repeated and it expects precisely one integer.
                Default value: 200
            optimizer (str): Specify which method to use to perform a geometry optimization.
                This keyword cannot be repeated and it expects precisely one keyword. BFGS is a
                quasi-newtonian method, and will best for "small" systems near the minimum. LBFGS
                is a limited memory version that can be used for "large" (>1000 atom) systems when
                efficiency outweights robustness. CG is more robust, especially when you are far from
                the minimum, but it slower.
                Default value: BFGS
        """
        super(RelaxSet, self).__init__(structure, **kwargs)

        global_section = Global(project_name=project_name, run_type="GEO_OPT")

        geo_opt_params = [
            Keyword("TYPE", "MINIMIZATION"),
            Keyword("MAX_DR", max_drift),
            Keyword("MAX_FORCE", max_force),
            Keyword("RMS_DR", 1.5e-3),
            Keyword("MAX_ITER", max_iter),
            Keyword("OPTIMIZER", optimizer),
        ]
        geo_opt = Section("GEO_OPT", subsections={}, keywords=geo_opt_params)

        if not self.check("MOTION"):
            self.insert(Section("MOTION", subsections={}))
        self["MOTION"].insert(geo_opt)

        self.insert(global_section)
        self.update(override_default_params)


# TODO: add defaults
class CellOptSet(DftSet):

    """
    CP2K input set containing the basic settings for performing geometry optimization. Values are all cp2k
    defaults, and should be good for most systems of interest.
    """

    def __init__(
        self,
        structure,
        project_name="CellOpt",
        override_default_params={},
        **kwargs
    ):

        """
        Args:
            structure:
            max_drift: Convergence criterion for the maximum geometry change between the current and the
                last optimizer iteration. This keyword cannot be repeated and it expects precisely one real.
                Default value: 3.00000000E-003
                Default unit: [bohr]
            max_force (float): Convergence criterion for the maximum force component of the current configuration.
                This keyword cannot be repeated and it expects precisely one real.
                Default value: 4.50000000E-004
                Default unit: [bohr^-1*hartree]
            max_iter (int): Specifies the maximum number of geometry optimization steps.
                One step might imply several force evaluations for the CG and LBFGS optimizers.
                This keyword cannot be repeated and it expects precisely one integer.
                Default value: 200
            optimizer (str): Specify which method to use to perform a geometry optimization.
                This keyword cannot be repeated and it expects precisely one keyword. BFGS is a
                quasi-newtonian method, and will best for "small" systems near the minimum. LBFGS
                is a limited memory version that can be used for "large" (>1000 atom) systems when
                efficiency outweights robustness. CG is more robust, especially when you are far from
                the minimum, but it slower.
                Default value: BFGS
        """
        super(CellOptSet, self).__init__(structure, **kwargs)
        global_section = Global(project_name=project_name, run_type="CELL_OPT")
        self.insert(global_section)
        self.update(override_default_params)


class HybridStaticSet(StaticSet):

    """
    Static calculation using hybrid DFT with the ADMM formalism in Cp2k.
    """

    def __init__(
        self,
        structure,
        hybrid_functional="PBE0",
        hf_fraction=0.25,
        project_name="Hybrid-Static",
        gga_x_fraction=0.75,
        gga_c_fraction=1,
        override_default_params={},
        **kwargs
    ):
        """
        Args:
            structure: pymatgen structure object
            method: hybrid dft method to use (currently select between HSE06 and PBE0)
            hf_fraction: percentage of exact HF to mix-in
            project_name: what to call this project
            gga_x_fraction: percentage of gga exchange to use
            gga_c_fraction: percentage of gga correlation to use
            override_default_params: override settings (see above).
        """
        super(HybridStaticSet, self).__init__(
            structure, project_name=project_name, **kwargs
        )
        self.activate_hybrid(
            hybrid_functional=hybrid_functional,
            hf_fraction=hf_fraction,
            gga_x_fraction=gga_x_fraction,
            gga_c_fraction=gga_c_fraction,
            **kwargs
        )
        self.update(override_default_params)


class HybridRelaxSet(RelaxSet):

    """
    Static calculation using hybrid DFT with the ADMM formalism in Cp2k.
    """

    def __init__(
        self,
        structure,
        hybrid_functional="PBE0",
        hf_fraction=0.25,
        project_name="Hybrid-Relax",
        gga_x_fraction=0.75,
        gga_c_fraction=1,
        override_default_params={},
        **kwargs
    ):
        """
        Args:
            structure: pymatgen structure object
            method: hybrid dft method to use (currently select between HSE06 and PBE0)
            hf_fraction: percentage of exact HF to mix-in
            project_name: what to call this project
            gga_x_fraction: percentage of gga exchange to use
            gga_c_fraction: percentage of gga correlation to use
            override_default_params: override settings (see above).
        """
        super(HybridRelaxSet, self).__init__(
            structure, project_name=project_name, **kwargs
        )
        self.activate_hybrid(
            hybrid_functional=hybrid_functional,
            hf_fraction=hf_fraction,
            gga_x_fraction=gga_x_fraction,
            gga_c_fraction=gga_c_fraction,
            **kwargs
        )
        self.update(override_default_params)


class HybridCellOptSet(CellOptSet):

    """
    Static calculation using hybrid DFT with the ADMM formalism in Cp2k.
    """

    def __init__(
        self,
        structure,
        hybrid_functional="PBE0",
        hf_fraction=0.25,
        project_name="Hybrid-Relax",
        gga_x_fraction=0.75,
        gga_c_fraction=1,
        override_default_params={},
        **kwargs
    ):
        """
        Args:
            structure: pymatgen structure object
            method: hybrid dft method to use (currently select between HSE06 and PBE0)
            hf_fraction: percentage of exact HF to mix-in
            project_name: what to call this project
            gga_x_fraction: percentage of gga exchange to use
            gga_c_fraction: percentage of gga correlation to use
            override_default_params: override settings (see above).
        """
        super(HybridCellOptSet, self).__init__(
            structure, project_name=project_name, **kwargs
        )
        self.activate_hybrid(
            hybrid_functional=hybrid_functional,
            hf_fraction=hf_fraction,
            gga_x_fraction=gga_x_fraction,
            gga_c_fraction=gga_c_fraction,
            **kwargs
        )
        self.update(override_default_params)
