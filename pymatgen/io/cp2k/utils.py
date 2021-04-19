"""
Utility functions for assisting with cp2k IO
"""

import os
import re
import numpy as np
from ruamel import yaml
from monty.serialization import loadfn
from monty.io import zopen
from pathlib import Path

from pymatgen.core import SETTINGS
from pymatgen.core.periodic_table import Element

MODULE_DIR = Path(__file__).resolve().parent


def _postprocessor(s):
    """
    Helper function to post process the results of the pattern matching functions in Cp2kOutput and turn them to
    python types.
    """
    s = s.rstrip()  # Remove leading/trailing whitespace
    s = s.replace(" ", "_")  # Remove whitespaces

    if s.lower() == "no":
        return False
    elif s.lower() == "none":
        return None
    elif s.lower() == "yes" or s.lower() == 'true':
        return True
    elif re.match(r"^-?\d+$", s):
        try:
            return int(s)
        except ValueError:
            raise IOError("Error in parsing CP2K file.")
    elif re.match(r"^[+\-]?(?=.)(?:0|[1-9]\d*)?(?:\.\d*)?(?:\d[eE][+\-]?\d+)?$", s):
        try:
            return float(s)
        except ValueError:
            raise IOError("Error in parsing CP2K file.")
    elif re.match(r"\*+", s):
        try:
            return np.NaN
        except ValueError:
            raise IOError("Error in parsing CP2K file.")
    else:
        return s


def _preprocessor(s, d='.'):
    """
    Cp2k contains internal preprocessor flags that are evaluated before
    execution. This helper function recognizes those preprocessor flags
    and replacees them with an equivalent cp2k input (this way everything
    is contained neatly in the cp2k input structure, even if the user
    preferred to use the flags.

    CP2K preprocessor flags (with arguments) are:

        @INCLUDE FILENAME: Insert the contents of FILENAME into the file at
            this location.
        @SET VAR VALUE: set a variable, VAR, to have the value, VALUE.
        $VAR or ${VAR}: replace these with the value of the variable, as set
            by the @SET flag.
        @IF/@ELIF: Not implemented yet.

    Args:
        s (str): string representation of cp2k input to preprocess
    """
    includes = re.findall(r"(@include.+)", s, re.IGNORECASE)
    for incl in includes:
        inc = incl.split()
        assert len(inc) == 2  # @include filename
        inc = inc[1].strip('\'')
        inc = inc.strip('\"')
        with zopen(os.path.join(d, inc)) as f:
            s = re.sub(r"{}".format(incl), f.read(), s)
    variable_sets = re.findall(r"(@SET.+)", s, re.IGNORECASE)
    for match in variable_sets:
        v = match.split()
        assert len(v) == 3  # @SET VAR value
        var, value = v[1:]
        s = re.sub(r"{}".format(match), "", s)
        s = re.sub(r"\${?"+var+"}?", value, s)

    c1 = re.findall(r"@IF", s, re.IGNORECASE)
    c2 = re.findall(r"@ELIF", s, re.IGNORECASE)
    if len(c1) > 0 or len(c2) > 0:
        raise NotImplementedError("This cp2k input processer does not currently "
                                  "support conditional blocks.")
    return s


def natural_keys(text):
    """
    Sort text by numbers coming after an underscore with natural number
    convention,
    Ex: [file_1, file_12, file_2] becomes [file_1, file_2, file_12]
    """
    def atoi(t):
        return int(t) if t.isdigit() else t
    return [atoi(c) for c in re.split(r'_(\d+)', text)]


def get_basis_and_potential(species, basis_and_potential_map, cardinality='DZVP', functional='PBE'):
    """
    Given a specie and a potential/basis type, this function accesses the available basis sets and potentials.
    Generally, the GTH potentials are used with the GTH basis sets.

    Args:
        species: (list) list of species for which to get the potential/basis strings
        d: (dict) a dictionary specifying how bases and/or potentials should be assigned to species
            E.g. {'Si': {'cardinality': 'DZVP', 'sr': True}, 'O': {'cardinality': 'TZVP'}}
        functional: (str) functional type. Default: 'PBE'
        basis_type: (str) the basis set type. Default: 'MOLOPT'
        cardinality: (str) basis set cardinality. Default: 'DZVP'

    Returns:
        (dict) of the form {'specie': {'potential': potential, 'basis': basis}...}
    """

    potential_filename = SETTINGS.get(
        "PMG_DEFAULT_CP2K_POTENTIAL_FILE", "GTH_POTENTIALS"
    )
    basis_filenames = ['BASIS_MOLOPT', 'BASIS_MOLOPT_UCL']

    functional = functional or SETTINGS.get(
        "PMG_DEFAULT_FUNCTIONAL", "PBE"
    )
    cardinality = cardinality or SETTINGS.get(
        "PMG_DEFAULT_BASIS_CARDINALITY", "DZVP"
    )
    basis_and_potential = {
        "basis_filenames": basis_filenames,
        "potential_filename": potential_filename,
    }

    if basis_and_potential_map == 'best':
        with open(os.path.join(MODULE_DIR, 'best_basis.yaml'), 'rt') as f:
            d = yaml.load(f, Loader=yaml.Loader)
    else:
        d = basis_and_potential_map
        for s in species:
            if s not in d:
                d[s] = {}
            if 'sr' not in d[s]:
                d[s]['sr'] = True
            if 'cardinality' not in d[s]:
                if s == 'Mg':
                    d[s]['cardinality'] = 'DZVPd'
                else:
                    d[s]['cardinality'] = cardinality

    with open(os.path.join(MODULE_DIR, 'basis_molopt.yaml'), 'rt') as f:
        data_b = yaml.load(f, Loader=yaml.Loader)
    with open(os.path.join(MODULE_DIR, 'gth_potentials.yaml'), 'rt') as f:
        data_p = yaml.load(f, Loader=yaml.Loader)

    for s in species:
        basis_and_potential[s] = {}
        if 'basis' in d[s]:
            b = [d[s]['basis']] if d[s]['basis'].upper() in [_.upper() for _ in data_b[s]] else []
        else:
            b = [_ for _ in data_b[s] if d[s]['cardinality'] in _.split('-')]
            if d[s]['sr'] and any(['SR' in _ for _ in b]):
                b = [_ for _ in b if 'SR' in _]
            else:
                b = [_ for _ in b if 'SR' not in _]
            if 'q' in d[s]:
                b = [_ for _ in b if d[s]['q'] in _]
            else:
                def srt(x):
                    return int(x.split('q')[-1])
                b = sorted(b, key=srt)[-1:]
        if len(b) == 0:
            raise LookupError('NO BASIS OF THAT TYPE AVAILABLE')
        elif len(b) > 1:
            raise LookupError('AMBIGUITY IN BASIS. PLEASE SPECIFY FURTHER')

        basis_and_potential[s]['basis'] = b[0]
        p = [_ for _ in data_p[s] if functional in _.split('-')]
        if len(p) == 0:
            raise LookupError('NO PSEUDOPOTENTIAL OF THAT TYPE AVAILABLE')
        if len(p) > 1:
            raise LookupError('AMBIGUITY IN POTENTIAL. PLEASE SPECIFY FURTHER')

        basis_and_potential[s]['potential'] = p[0]

    return basis_and_potential


def get_aux_basis(basis_type, default_basis_type='cpFIT'):
    """
    Get auxiliary basis info for a list of species.

    Args:
        basis_type (dict): dict of auxiliary basis sets to use. i.e:
            basis_type = {'Si': 'cFIT', 'O': 'cpFIT'}. Basis type needs to
            exist for that species.

            Basis types:
                FIT
                cFIT
                pFIT
                cpFIT
                GTH-def2
                aug-{FIT,cFIT,pFIT,cpFIT, GTH-def2}

        default_basis_type (str) default basis type if n

    """
    default_basis_type = default_basis_type or SETTINGS.get("PMG_CP2K_DEFAULT_AUX_BASIS_TYPE")
    basis_type = {k: basis_type[k] if basis_type[k] else default_basis_type for k in basis_type}
    basis = {k: {} for k in basis_type}
    aux_bases = loadfn(os.path.join(MODULE_DIR, 'aux_basis.yaml'))
    for k in basis_type:
        if aux_bases.get(k) is None:
            basis[k] = None
            continue
        for i in aux_bases[k]:
            if i.startswith(basis_type[k]):
                basis[k] = i
                break
    for k in basis:
        if not basis[k]:
            if aux_bases[k]:
                basis[k] = aux_bases[k][0]
            else:
                raise LookupError("NO BASIS OF THAT TYPE")
    return basis


def get_unique_site_indices(structure):
    """
    Get unique site indices for a structure according to site properties. Whatever site-property has the most
    unique values is used for indexing.

    For example, if you have magnetic CoO with half Co atoms having a positive moment, and the other
    half having a negative moment. Then this function will create a dict of sites for Co_1, Co_2, O. This
    function also deals with "Species" properties like oxi_state and spin by pushing them to site
    properties.

    This creates unique sites, based on site properties, but does not have anything to do with turning
    those site properties into CP2K input parameters. This will only be done for properties which can be
    turned into CP2K input parameters, which are stored in parsable_site_properties.
    """
    spins = []
    oxi_states = []
    parsable_site_properties = {
        'magmom', 'oxi_state', 'spin',
        'u_minus_j', 'basis', 'potential',
        'ghost'
    }

    for site in structure:
        for sp, occu in site.species.items():
            oxi_states.append(getattr(sp, "oxi_state", 0))
            spins.append(getattr(sp, "_properties", {}).get('spin', 0))

    structure.add_site_property('oxi_state', oxi_states)
    structure.add_site_property('spin', spins)
    structure.remove_oxidation_states()
    items = [
        (
            site.species_string,
            *[structure.site_properties[k][i]
              for k in structure.site_properties if k.lower() in parsable_site_properties]
        )
        for i, site in enumerate(structure)
    ]
    unique_itms = list(set(items))
    _sites = {u: [] for u in unique_itms}
    for i, itm in enumerate(items):
        _sites[itm].append(i)
    sites = {}
    nums = {s: 1 for s in structure.symbol_set}
    for s in _sites:
        sites['{}_{}'.format(s[0], nums[s[0]])] = _sites[s]
        nums[s[0]] += 1

    return sites


def get_cutoff_from_basis(els, bases, rel_cutoff=50):
    """
    Gets the appropriate cutoff for the calculation given the elements/basis sets being used
    and the desired relative cutoff.

    Args:
        els: list of element symbols
        bases: corresponding basis set names for the elements
        rel_cutoff: The desired relative cutoff
    Returns:
        Ideal cutoff for calculation.
    """
    with open(os.path.join(MODULE_DIR, 'basis_largest_exponents.yaml'), 'rt') as f:
        _exponents = yaml.load(f, Loader=yaml.Loader)
        exponents = {el.upper(): {b.upper(): v for b, v in basis.items()} for el, basis in _exponents.items()}
        return max(
            [
                np.ceil(exponents[el.upper()][basis.upper()])*rel_cutoff for el, basis in zip(els, bases)
            ]
        )


def generate_max_exponents(files):
    """
    Helper function to get the max exponent for each basis set and write to a yaml file.
    """
    exponents = {}
    els = [str(e) for e in Element]
    for file in files:
        with open(file) as f:
            lines = f.readlines()
            i = 0
            while i < len(lines):
                if lines[i].split():
                    el = lines[i].split()[0]
                    basis = lines[i].split()[-1]
                    if el in els:
                        i += 3
                        if el not in exponents:
                            exponents[el] = {}
                        exponents[el][basis] = float(lines[i].split()[0])
                i += 1
    with open(os.path.join(MODULE_DIR, 'basis_largest_exponents.yaml'), 'w') as f:
        yaml.dump(exponents, f, Dumper=yaml.Dumper, allow_unicode=True, default_flow_style=False)
