import math
from pathlib import Path
from pymatgen.io.cp2k.inputs import *
from pymatgen.io.cp2k.utils import get_basis_and_potential

import warnings

__author__ = "Nicholas Winner"
__version__ = "0.2"
__email__ = "nwinner@berkeley.edu"
__date__ = "January 2019"

MODULE_DIR = Path(__file__).resolve().parent

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

from pymatgen.core.structure import Structure
class Cp2kInputSet(Cp2kInput):

    def __init__(self, structure, potential_and_basis={}, kppa=1000,
                 auto_supercell=True,
                 override_default_params={}, **kwargs):
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

        if auto_supercell:
            # Only sample gamma point, need k-points / num_atoms to be < 8
            structure.make_supercell(math.ceil(math.pow(kppa / structure.num_sites / 8, 1/3)))

        coord = Coord(structure)
        cell = Cell(structure.lattice)
        subsys = Subsys(subsections={'CELL': cell, 'COORD': coord})

        basis_type = potential_and_basis.get("basis", 'MOLOPT')
        potential = potential_and_basis.get("potential", 'GTH')
        functional = potential_and_basis.get('functional', 'PBE')
        cardinality = potential_and_basis.get('cardinality', 'DZVP')
        basis_and_potential = get_basis_and_potential(structure.symbol_set, potential_type=potential,
                                                      functional=functional, basis_type=basis_type,
                                                      cardinality=cardinality)

        for s in structure.symbol_set:
            subsys.insert(Kind(s, basis_set=basis_and_potential[s]['basis'],
                          potential=basis_and_potential[s]['potential']))

        if 'FORCE_EVAL' not in self.subsections.keys():
            self.insert(ForceEval())

        self['FORCE_EVAL'].insert(subsys)

        self.update(override_default_params)

    def activate_hybrid(self, structure, method='HSE06', hf_fraction=0.25, gga_x_fraction=0.75, gga_c_fraction=1):

        """
        Basic set for activating hybrid DFT calculation using Auxiliary Density Matrix Method.
        """

        self['FORCE_EVAL']['DFT'].keywords.append(Keyword('BASIS_SET_FILE_NAME', 'BASIS_ADMM'))

        for k, v in self['FORCE_EVAL']['SUBSYS'].subsections.items():
            if 'KIND' in k:
                v.keywords.append(Keyword('BASIS_SET', 'AUX_FIT', 'cFIT3'))

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

        if method is 'HSE06':
            ip_keywords = [
                Keyword('POTENTIAL_TYPE', 'SHORTRANGE'),
                Keyword('OMEGA', 0.11)  # TODO This value should be tested, can we afford to increase it to an opt val?
            ]
        elif method is 'PBE0':
            ip_keywords = [
                Keyword('POTENTIAL_TYPE', 'TRUNCATED'),
                Keyword('CUTOFF_RADIUS', min(structure.lattice.abc) / 2),
                Keyword('T_C_G_DATA', 't_c_g.data')
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


class QsSet(Cp2kInputSet):
    """
    Base for an input set using the Quickstep module (i.e. a DFT calculation).

    Optional Params:
        BAND_GAP:
        MINIMIZER
        PRECONDITIONER

        CUTOFF
        REL_CUTOFF
        NGRIDS
        PROGRESSION_FACTOR
    """

    def __init__(self, structure, ot=True, band_gap=0.01, eps_default=1e-6,
                 minimizer='CG', preconditioner='FULL_ALL',
                 cutoff=1200, rel_cutoff=80, ngrids=5, progression_factor=3, multigrid_set=False,
                 override_default_params={}, **kwargs):

        super(QsSet, self).__init__(structure, **kwargs)

        # Build the QS Section
        qs = QS(eps_default=eps_default)
        scf = Scf(subsections={})

        # If there's a band gap, always use OT, else use Davidson
        if band_gap or ot:
            scf.insert(OrbitalTransformation(minimizer=minimizer,
                                             preconditioner=preconditioner,
                                             energy_gap=band_gap))
            scf.insert(Section('OUTER_SCF', subsections={},
                               keywords=[
                                   Keyword('MAX_SCF', kwargs.get('outer_max_scf', 20)),
                                   Keyword('EPS_SCF', kwargs.get('outer_eps_scf', 1e-6))
                               ]))
        else:
            scf.insert(Section('DIAGONALIZATION', subsections={}))
            scf['DIAGONALIZATION'].insert(Section('DAVIDSON', subsections={}))

        mgrid = Mgrid(cutoff=cutoff, rel_cutoff=rel_cutoff, ngrids=ngrids, progression_factor=progression_factor,
                      multigrid_set=multigrid_set)
        dft = Dft(subsections={'QS': qs, 'SCF': scf, 'MGRID': mgrid})

        # Create subsections and insert into them
        self['FORCE_EVAL'].insert(dft)
        xc_functional = XC_FUNCTIONAL(functional='PBE')
        xc = Section('XC', subsections={'XC_FUNCTIONAL': xc_functional})
        self['FORCE_EVAL']['DFT'].insert(xc)
        self['FORCE_EVAL']['DFT'].insert(Section('PRINT', subsections={}))
        self['FORCE_EVAL']['DFT']['PRINT'].insert(PDOS())
        self['FORCE_EVAL']['DFT']['PRINT']['PDOS'].insert(LDOS())

        self.update(override_default_params)

    def print_pdos(self, nlumo=-1):
        if not self.check('FORCE_EVAL/DFT/PRINT'):
            self['FORCE_EVAL']['DFT'].insert(Section('PRINT', subsections={}))
        self['FORCE_EVAL']['DFT']['PRINT'].insert(PDOS(nlumo=nlumo))

    def print_mo_cubes(self, write_cube=False, nlumo=1, nhomo=1):
        if not self.check('FORCE_EVAL/DFT/PRINT'):
            self['FORCE_EVAL']['DFT'].insert(Section('PRINT', subsections={}))
        self['FORCE_EVAL']['DFT']['PRINT'].insert(MO_Cubes(write_cube=write_cube, nlumo=nlumo, nhomo=nhomo))

    def print_hartree_potential(self):
        if not self.check('FORCE_EVAL/DFT/PRINT'):
            self['FORCE_EVAL']['DFT'].insert(Section('PRINT', subsections={}))
        self['FORCE_EVAL']['DFT']['PRINT'].insert(V_Hartree_Cube())

    def print_e_density(self):
        if not self.check('FORCE_EVAL/DFT/PRINT'):
            self['FORCE_EVAL']['DFT'].insert(Section('PRINT', subsections={}))
        self['FORCE_EVAL']['DFT']['PRINT'].insert(E_Density_Cube())


class StaticSet(QsSet):
    """
    Basic static energy calculation. Turns on Quickstep module, sets the run_type in global,
    and uses structure object to build the subsystem.
    """

    def __init__(self, structure, project_name='CP2K', run_type='ENERGY_FORCE',
                 override_default_params={}, **kwargs):
        super(StaticSet, self).__init__(structure, **kwargs)
        global_section = Global(project_name=project_name, run_type=run_type)
        self.insert(global_section)
        self.update(override_default_params)


class RelaxSet(QsSet):

    def __init__(self, structure, max_drift=1e-3, max_force=1e-3, max_iter=200,
                 optimizer='CG', override_default_params={}, **kwargs):

        """
        CP2K input set containing the basic settings for performing geometry optimization. Description
        of args for geo optimization are from CP2K Manual

        Args:
            structure:
            max_drift: Convergence criterion for the maximum geometry change between the current and the last optimizer iteration.
                This keyword cannot be repeated and it expects precisely one real.
                Default value: 3.00000000E-003
                Default unit: [bohr]
            max_force (float): Convergence criterion for the maximum force component of the current configuration.
                This keyword cannot be repeated and it expects precisely one real.
                Default value: 4.50000000E-004
                Default unit: [bohr^-1*hartree]
            max_iter (int): Specifies the maximum number of geometry optimization steps. One step might imply several force evaluations for the CG and LBFGS optimizers.
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

        s = structure.composition.formula.replace(' ', '-')
        global_section = Global(project_name='{}-GEO_OPT'.format(s), run_type='GEO_OPT')

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

        motion = Section('MOTION', subsections={'GEO_OPT': geo_opt})

        self.insert(global_section)
        self.insert(motion)
        self.update(override_default_params)


class HybridStaticSet(StaticSet):

    def __init__(self, structure, override_default_params={}, **kwargs):
        super(HybridStaticSet, self).__init__(structure, **kwargs)
        self.activate_hybrid(structure, method='HSE06', hf_fraction=0.25, gga_x_fraction=0.75, gga_c_fraction=1)
        self.update(override_default_params)


class HybridRelaxSet(RelaxSet):

    def __init__(self, structure, override_default_params={}, **kwargs):
        super(HybridRelaxSet, self).__init__(structure, **kwargs)
        self.activate_hybrid(structure, method='HSE06', hf_fraction=0.25, gga_x_fraction=0.75, gga_c_fraction=1)
        self.update(override_default_params)


