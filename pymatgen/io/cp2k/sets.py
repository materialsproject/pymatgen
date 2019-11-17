import os
from pathlib import Path
from pymatgen.io.cp2k.inputs import *
from pymatgen.io.cp2k.utils import get_basis_and_potential

import warnings

__author__ = "Nicholas Winner"

MODULE_DIR = Path(__file__).resolve().parent


class Cp2kInputSet(MSONable):

    def __init__(self, structure, sections={}, potential_and_basis={}, override_default_settings={}):
        """
        The basic representation of a CP2K input set as a collection of "sections" defining the simulation
        connected to a structure object. At the most basis level, CP2K requires a &GLOBAL section and
        &FORCE_EVAL section. Global sets parameters like "RUN_TYPE" or the overall verbosity. FORCE_EVAL is
        the largest section usually, containing the cell and coordinates of atoms, the DFT settings, and more.
        This top level input set is meant to initialize GLOBAL and FORCE_EVAL based on a structure object and
        and sections that the user provides. Generally, this class will not be used directly, and instead one of
        its child-classes will be used, which contain more predefined initializations of various sections, and,
        if modifications are required, the user can specify override_default_settings.

        Args:
            structure: (Structure) pymatgen structure object used to define the lattice, coordinates, and elements

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

            override_default_settings: (dict) Specifies user-defined settings to override the settings of any
                input set.
        """
        if "GLOBAL" not in sections.keys():
            warnings.warn("No GLOBAL section found in arguments. A DEFAULT GLOBAL section will be used.")
            sections.update(Global())
        if "FORCE_EVAL" not in sections.keys():
            warnings.warn("No FORCE_EVAL section found in arguments. This section is required, be careful.")
            sections.update(ForceEval())

        coords = Coord(structure)
        cell = Cell(structure.lattice)

        basis_set = potential_and_basis.get("basis", 'DZVP')
        potential = potential_and_basis.get("potential", 'GTH')
        functional = potential_and_basis.get('functional', 'PBE')
        basis_and_potential = get_basis_and_potential(structure.symbol_set, potential_type=potential,
                                                      functional=functional, basis_type=basis_set)

        kind = [Kind(params={'specie': s, 'BASIS_SET': basis_and_potential[s]['basis'],
                             "POTENTIAL": basis_and_potential[s]['potential']}) for s in structure.symbol_set]

        subsys_params = {}
        subsys_params.update(coords)
        subsys_params.update(cell)
        for k in kind: subsys_params.update(k)
        subsys = Subsys(params=subsys_params)

        sections['FORCE_EVAL'].update(subsys)

        sections.update(override_default_settings)
        self.sections = sections

    def write_input(self, input_filename='cp2k.inp', output_dir='.', make_dir_if_not_present=True):
        if not os.path.isdir(output_dir) and make_dir_if_not_present:
            os.mkdir(output_dir)
        filepath = os.path.join(output_dir, input_filename)
        Input(section_name='CP2K_INPUT', section_keyword=None, section_params=self.sections).write_file(filepath)


class QsSet(Cp2kInputSet):
    """
    Base for an input set using the Quickstep module (i.e. a DFT calculation).
    """

    def __init__(self, structure, sections={}, **kwargs):

        mgrid = Mgrid()
        qs  = QS()
        scf = Scf()

        dft = Dft(params={**mgrid, **qs, **scf})
        force_eval = ForceEval(params=dft)

        sections.update(force_eval)

        super(QsSet, self).__init__(structure, sections=sections, **kwargs)


class RelaxSet(QsSet):

    def __init__(self, structure, **kwargs):

        sections = {}

        sections.update(Global(params={'RUN_TYPE': 'GEO_OPT'}))

        # TODO: This is hardcoded while we are testing. In the future, EVERYTHING in in the "sets" module
        # TODO: Should only use defaults in the initialization, and the keyword "params" will let users overwrite
        geo_opt_params = {
            'TYPE': 'MINIMIZATION',
            'MAX_DR': 1.0e-3,
            'MAX_FORCE': 1.0e-3,
            'RMS_DR': 1.0e-3,
            'MAX_ITER': 200,
            'OPTIMIZER': 'CG',
            'CG': {
                'MAX_STEEP_STEPS': 0,
                'RESTART_LIMIT': 9.0e-1
            }
        }
        geo_opt= Section(section_name='GEO_OPT', section_keyword=None,
                         section_params=geo_opt_params)

        motion_params = {}
        motion_params.update(geo_opt)
        motion = Section(section_name='MOTION', section_keyword=None,
                         section_params=motion_params)

        sections.update(motion)

        super(RelaxSet, self).__init__(structure, sections=sections, **kwargs)