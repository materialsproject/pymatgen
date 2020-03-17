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

import math
from pathlib import Path
from pymatgen.io.cp2k.inputs import Cp2kInput, Section, Keyword, Global, ForceEval, Mgrid, MO_Cubes, \
    OrbitalTransformation, XC_FUNCTIONAL, V_Hartree_Cube, Dft, E_Density_Cube, LDOS, PDOS, \
    Coord, Kind, QS, PBE, Cell, Subsys, Scf
from pymatgen.io.cp2k.utils import get_basis_and_potential, get_aux_basis, get_unique_site_indices

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

    def __init__(self, structure, potential_and_basis={},
                 kppa=1000, break_symmetry=False,
                 multiplicity=0, project_name='CP2K',
                 override_default_params={}, **kwargs):
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
        super(Cp2kInputSet, self).__init__(name='CP2K_INPUT', subsections={})

        # Important CP2K set parameters
        self.structure = structure
        self.charge = structure.charge
        self.multiplicity = multiplicity  # spin multiplicity = 2s+1
        self.override_default_params = override_default_params
        self.project_name = project_name
        self.kppa = kppa

        # This is a simple way to make a supercell according to the number of k-points / atom
        # TODO: should be anisotropic, scaling the directions as needed
        if kppa > 0:
            # Only sample gamma point, need k-points / num_atoms to be < 8.
            structure.make_supercell(math.ceil(math.pow(kppa / structure.num_sites / 8, 1/3)))

        cell = Cell(structure.lattice)
        subsys = Subsys(subsections={'CELL': cell})

        # Decide what basis sets/pseudopotentials to use
        basis_type = potential_and_basis.get("basis", 'MOLOPT')
        potential = potential_and_basis.get("potential", 'GTH')
        cardinality = potential_and_basis.get('cardinality', 'DZVP')
        basis_and_potential = get_basis_and_potential(structure.symbol_set, functional=potential, basis_type=basis_type,
                                                      cardinality=cardinality)
        self.basis_set_file_name = basis_and_potential['basis_filename']
        self.potential_file_name = basis_and_potential['potential_filename']

        # Insert atom kinds by identifying the unique sites (unique element and site properties)
        unique_kinds = get_unique_site_indices(structure)
        for k in unique_kinds.keys():
            kind = k.split('_')[0]
            subsys.insert(Kind(kind, alias=k, basis_set=basis_and_potential[kind]['basis'],
                          potential=basis_and_potential[kind]['potential']))
        coord = Coord(structure, alias=unique_kinds)
        subsys.insert(coord)

        if not self.check('FORCE_EVAL'):
            self.insert(ForceEval())

        self['FORCE_EVAL'].insert(subsys)
        self.print_forces()
        self.print_motion()

        self.update(override_default_params)

    def create_structure(self, structure=None):
        """
        Create the structure for the input
        """
        pass

    def activate_hybrid(self, structure, method='HSE06', hf_fraction=0.25, gga_x_fraction=0.75, gga_c_fraction=1):

        """
        Basic set for activating hybrid DFT calculation using Auxiliary Density Matrix Method.
        """
        basis = get_aux_basis(self.structure.symbol_set)
        self['FORCE_EVAL']['DFT'].keywords.extend([
            Keyword('BASIS_SET_FILE_NAME', b) for b in basis['basis_filename']])

        for k, v in self['FORCE_EVAL']['SUBSYS'].subsections.items():
            if 'KIND' in k:
                kind = v.get_keyword('ELEMENT').values[0]
                v.keywords.append(Keyword('BASIS_SET', 'AUX_FIT', basis[kind]['basis']))

        aux_matrix_params = [
            Keyword('ADMM_PURIFICATION_METHOD', 'MO_DIAG'),
            Keyword('METHOD', 'BASIS_PROJECTION')
        ]
        aux_matrix = Section("AUXILIARY_DENSITY_MATRIX_METHOD", keywords=aux_matrix_params,
                             subsections={})

        # Define the GGA functional as PBE
        pbe = PBE('ORIG', scale_c=gga_c_fraction, scale_x=gga_x_fraction)
        xc = XC_FUNCTIONAL('PBE', subsections={'PBE': pbe})
        xc_functional = Section('XC_FUNCTIONAL', subsections={'PBE': pbe})

        # Define screening of the HF term
        screening = Section('SCREENING', subsections={},
                            keywords=[Keyword('EPS_SCHWARZ', 1e-7),
                                      Keyword('EPS_SCHWARZ_FORCES', 1e-7),
                                      Keyword('SCREEN_ON_INITIAL_P', True),
                                      Keyword('SCREEN_P_FORCES', True)])

        if method == 'HSE06':
            xc_functional.insert(Section('XWPBE', subsections={},
                                         keywords=[Keyword('SCALE_X0', 1),
                                                   Keyword('SCALE_X', -hf_fraction),
                                                   Keyword('OMEGA', 0.11)]))
            ip_keywords = [
                Keyword('POTENTIAL_TYPE', 'SHORTRANGE'),
                Keyword('OMEGA', 0.11)  # TODO This value should be tested, can we afford to increase it to an opt val?
            ]
        elif method == 'PBE0':
            ip_keywords = [
                Keyword('POTENTIAL_TYPE', 'TRUNCATED'),
                Keyword('CUTOFF_RADIUS', structure.lattice.matrix[structure.lattice.matrix.nonzero()].min() / 2),
                Keyword('T_C_G_DATA', 't_c_g.dat')
            ]
        interaction_potential = Section('INTERACTION_POTENTIAL', subsections={},
                                        keywords=ip_keywords)

        load_balance = Section('LOAD_BALANCE', keywords=[Keyword('RANDOMIZE', True)], subsections={})
        memory = Section('MEMORY', subsections={},
                         keywords=[Keyword('EPS_STORAGE_SCALING', 0.1),
                                   Keyword('MAX_MEMORY', 5000)])
        hf = Section('HF', keywords=[Keyword('FRACTION', hf_fraction)],
                     subsections={'SCREENING': screening, 'INTERACTION_POTENTIAL': interaction_potential,
                                  'LOAD_BALANCE': load_balance, 'MEMORY': memory})

        xc = Section('XC', subsections={'XC_FUNCTIONAL': xc_functional, 'HF': hf})

        self.subsections['FORCE_EVAL']['DFT'].insert(aux_matrix)
        self.subsections['FORCE_EVAL']['DFT'].insert(xc)

    def print_forces(self):
        """
        Print out the forces and stress during calculation
        """
        self['FORCE_EVAL'].insert(Section('PRINT', subsections={}))
        self['FORCE_EVAL']['PRINT'].insert(Section('FORCES', subsections={}))
        self['FORCE_EVAL']['PRINT'].insert(Section('STRESS_TENSOR', subsections={}))

    def print_motion(self):
        """
        Print the motion info (trajectory, cell, forces, stress
        """
        if not self.check('MOTION'):
            self.insert(Section('MOTION', subsections={}))
        self['MOTION'].insert(Section('PRINT', subsections={}))
        self['MOTION']['PRINT'].insert(Section('TRAJECTORY',
                                               section_parameters=['HIGH'],
                                               subsections={}))
        self['MOTION']['PRINT'].insert(Section('CELL', subsections={}))
        self['MOTION']['PRINT'].insert(Section('FORCES', subsections={}))
        self['MOTION']['PRINT'].insert(Section('STRESS', subsections={}))


class DftSet(Cp2kInputSet):
    """
    Base for an input set using the Quickstep module (i.e. a DFT calculation). The DFT section is pretty vast
    in CP2K, so this set hopes to make the DFT setup fairly simple. The provided parameters are pretty conservative,
    and so they should not need to be changed very often.
    """

    def __init__(self, structure, ot=True, band_gap=0.01, eps_default=1e-12,
                 eps_scf=1e-7, max_scf=50, minimizer='DIIS', preconditioner='FULL_ALL',
                 cutoff=1200, rel_cutoff=80, ngrids=5, progression_factor=3,
                 override_default_params={}, **kwargs):
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

        # Build the QS Section
        qs = QS(eps_default=eps_default)
        scf = Scf(eps_scf=eps_scf, max_scf=max_scf, subsections={})

        # If there's a band gap, always use OT, else use Davidson
        if band_gap or ot:
            scf.insert(OrbitalTransformation(minimizer=minimizer,
                                             preconditioner=preconditioner,
                                             energy_gap=band_gap))
            scf.insert(Section('OUTER_SCF', subsections={},
                               keywords=[
                                   Keyword('MAX_SCF', kwargs.get('outer_max_scf', 20)),
                                   Keyword('EPS_SCF', kwargs.get('outer_eps_scf', eps_scf))
                               ]))
        else:
            scf.insert(Section('DIAGONALIZATION', subsections={}))
            scf['DIAGONALIZATION'].insert(Section('DAVIDSON', subsections={}))

        # Create the multigrid for FFTs
        mgrid = Mgrid(cutoff=cutoff, rel_cutoff=rel_cutoff, ngrids=ngrids, progression_factor=progression_factor)

        # Set the DFT calculation with global parameters
        dft = Dft(MULTIPLICITY=self.multiplicity, CHARGE=self.charge,
                  basis_set_filename=self.basis_set_file_name,
                  potential_filename=self.potential_file_name,
                  subsections={'QS': qs, 'SCF': scf, 'MGRID': mgrid})

        # Create subsections and insert into them
        self['FORCE_EVAL'].insert(dft)
        xc_functional = XC_FUNCTIONAL(functional='PBE')
        xc = Section('XC', subsections={'XC_FUNCTIONAL': xc_functional})
        self['FORCE_EVAL']['DFT'].insert(xc)
        self['FORCE_EVAL']['DFT'].insert(Section('PRINT', subsections={}))

        self.print_pdos()
        self.print_ldos()
        self.print_mo_cubes()

        self.update(override_default_params)

    def print_pdos(self, nlumo=-1):
        """
        Activate creation of the PDOS file.

        Args:
            nlumo (int): Number of virtual orbitals to be added to the MO set (-1=all).
                CAUTION: Setting this value to be higher than the number of states present may cause a Cholesky error.
        """
        if not self.check('FORCE_EVAL/DFT/PRINT/PDOS'):
            self['FORCE_EVAL']['DFT']['PRINT'].insert(PDOS(nlumo=nlumo))

    def print_ldos(self, nlumo=-1):
        """
        Activate the printing of LDOS files, printing one for each atom kind by default

        Args:
            nlumo (int): Number of virtual orbitals to be added to the MO set (-1=all).
                CAUTION: Setting this value to be higher than the number of states present may cause a Cholesky error.
        """
        if not self.check('FORCE_EVAL/DFT/PRINT/PDOS'):
            self['FORCE_EVAL']['DFT']['PRINT'].insert(PDOS(nlumo=nlumo))
        for i in range(len(self.structure)):
            self['FORCE_EVAL']['DFT']['PRINT']['PDOS'].insert(LDOS(i))

    def print_mo_cubes(self, write_cube=False, nlumo=1, nhomo=1):
        """
        Activate printing of molecular orbitals.

        Args:
            write_cube (bool): whether to write cube file for the MOs (setting false will just print levels in out file)
            nlumo (int): Controls the number of lumos that are printed and dumped as a cube (-1=all)
            nhomo (int): Controls the number of homos that are printed and dumped as a cube (-1=all)
        """
        if not self.check('FORCE_EVAL/DFT/PRINT/MO_CUBES'):
            self['FORCE_EVAL']['DFT']['PRINT'].insert(MO_Cubes(write_cube=write_cube, nlumo=nlumo, nhomo=nhomo))

    def print_hartree_potential(self):
        """
        Controls the printing of a cube file with eletrostatic potential generated by the total density
        (electrons+ions). It is valid only for QS with GPW formalism.
        Note that by convention the potential has opposite sign than the expected physical one.
        """
        if not self.check('FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE'):
            self['FORCE_EVAL']['DFT']['PRINT'].insert(V_Hartree_Cube())

    def print_e_density(self):
        """
        Controls the printing of cube files with the electronic density and, for LSD calculations, the spin density
        """
        if not self.check('FORCE_EVAL/DFT/PRINT/E_DENSITY_CUBE'):
            self['FORCE_EVAL']['DFT']['PRINT'].insert(E_Density_Cube())

    def set_charge(self, charge):
        """
        Set the overall charge of the simulation cell
        """
        self['FORCE_EVAL']['DFT'].set_keyword(Keyword('CHARGE', charge))


class StaticSet(DftSet):

    """
    Basic static energy calculation. Turns on Quickstep module, sets the run_type in global,
    and uses structure object to build the subsystem.
    """

    def __init__(self, structure, project_name='Static', run_type='ENERGY_FORCE',
                 override_default_params={}, **kwargs):
        """
        Args:
            structure: Pymatgen structure object
            project_name (str): What to name this cp2k project (controls naming of files printed out)
            run_type (str): Run type. As a static set it should be one of the static aliases, like 'ENERGY_FORCE'
        """
        super(StaticSet, self).__init__(structure, **kwargs)
        global_section = Global(project_name=project_name, run_type=run_type)
        self.insert(global_section)
        self.update(override_default_params)


class RelaxSet(DftSet):

    """
    CP2K input set containing the basic settings for performing geometry optimization. Description
    of args for geo optimization are from CP2K Manual
    """

    def __init__(self, structure, max_drift=1e-3, max_force=1e-3, max_iter=200,
                 project_name='Relax', optimizer='CG', override_default_params={}, **kwargs):

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
                This keyword cannot be repeated and it expects precisely one keyword. BFGS is
                the CP2K default as it is more efficient for small systems, but CG is more robust
                overall, and thus a better choice if you don't know what type of system will be
                calculated.
                Default value: CG
        """
        super(RelaxSet, self).__init__(structure, **kwargs)

        global_section = Global(project_name=project_name, run_type='GEO_OPT')

        geo_opt_params = [
            Keyword('TYPE', 'MINIMIZATION'),
            Keyword('MAX_DR', max_drift),
            Keyword('MAX_FORCE', max_force),
            Keyword('RMS_DR', 1.0e-3),
            Keyword('MAX_ITER', max_iter),
            Keyword('OPTIMIZER', optimizer)
        ]
        cg = Section('CG', subsections={}, keywords=[Keyword('MAX_STEEP_STEPS', 0),
                                                     Keyword('RESTART_LIMIT', 9.0e-1)])
        geo_opt = Section('GEO_OPT', subsections={}, keywords=geo_opt_params)
        geo_opt.insert(cg)

        if not self.check('MOTION'):
            self.insert(Section('MOTION', subsections={}))
        self['MOTION'].insert(geo_opt)

        self.insert(global_section)
        self.update(override_default_params)


class HybridStaticSet(StaticSet):

    """
    Static calculation using hybrid DFT with the ADMM formalism in Cp2k.
    """

    def __init__(self, structure, method='HSE06', hf_fraction=0.25, project_name='Hybrid-Static',
                 gga_x_fraction=0.75, gga_c_fraction=1, override_default_params={}, **kwargs):
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
        super(HybridStaticSet, self).__init__(structure, project_name=project_name, **kwargs)
        self.activate_hybrid(structure, method=method, hf_fraction=hf_fraction,
                             gga_x_fraction=gga_x_fraction, gga_c_fraction=gga_c_fraction)
        self.update(override_default_params)


class HybridRelaxSet(RelaxSet):

    """
    Static calculation using hybrid DFT with the ADMM formalism in Cp2k.
    """

    def __init__(self, structure, method='HSE06', hf_fraction=0.25, project_name='Hybrid-Relax',
                 gga_x_fraction=0.75, gga_c_fraction=1, override_default_params={}, **kwargs):
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
        super(HybridRelaxSet, self).__init__(structure, project_name=project_name, **kwargs)
        self.activate_hybrid(structure, method=method, hf_fraction=hf_fraction,
                             gga_x_fraction=gga_x_fraction, gga_c_fraction=gga_c_fraction)
        self.update(override_default_params)
