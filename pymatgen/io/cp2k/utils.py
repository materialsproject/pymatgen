import os
import re
from monty.io import zopen
from ruamel import yaml
from pathlib import Path

from pymatgen import SETTINGS

MODULE_DIR = Path(__file__).resolve().parent


# TODO: only loads the GTH type potentials. (Usually this is what is used with CP2K, but more exist.)
def update_potentials(directory=None, potential_files=[]):
    """
    Updates the POTENTIALS.yaml files in the cp2k module.

    Can give this function a different data directory from which to read the files.
    This can be useful if you have custom potential files.
    """
    if directory:
        cp2k_data = directory
    else:
        cp2k_data = SETTINGS.get("PMG_CP2K_DATA_DIR", '')
    if not potential_files:
        potential_files = [
            'GTH_POTENTIALS',
            'NLCC_POTENTIALS'
        ]

    class CustomDumper(yaml.Dumper):
        # Super neat hack to preserve the mapping key order. See https://stackoverflow.com/a/52621703/1497385
        # Preserves ordering of the elements in the potential.yaml file
        def represent_dict_preserve_order(self, data):
            return self.represent_dict(data.items())
    CustomDumper.add_representer(dict, CustomDumper.represent_dict_preserve_order)

    for potential_file in potential_files:
        potentials = read_potentials(filename=os.path.join(cp2k_data, potential_file))
        with open('{}.yaml'.format(potential_file), 'w') as outfile:
            yaml.dump(potentials, outfile, default_flow_style=False, Dumper=CustomDumper)


def read_potentials(filename):
    """
    Reads in the pseudopotentials information from GTH formatted potential file as
    is found in cp2k/data directory

    Args:
        filename: (str) filename to be parsed.

    Returns:
         (dict) representation of the pseudopotential file
    """
    with zopen(filename) as f:
        lines = f.readlines()
    potentials = {}
    current_functional = ''
    for i in range(len(lines)):
        if lines[i].__contains__('functional'):
            match = re.search(r'(\w+) functional', lines[i])
            current_functional = match.groups(0)[0]
            potentials[current_functional] = {}
        if current_functional:
            if lines[i].__contains__('GTH'):
                l = lines[i].split()
                element = l[0]
                potential = l[1]
                potentials[current_functional][element] = {}
                potentials[current_functional][element][potential] = {}
                potentials[current_functional][element][potential]['alias'] = [s for s in l[2:]]

                i += 1
                potentials[current_functional][element][potential]['nelect'] = \
                    [int(s) for s in lines[i].split()]
                i += 1

                l3 = lines[i].split()  # r_loc nexp_ppl cexp_ppl(1) ... cexp_ppl(nexp_ppl)
                potentials[current_functional][element][potential]['r_loc'] = \
                    float(l3[0])
                potentials[current_functional][element][potential]['nexp_ppl'] = \
                    int(l3[1])
                potentials[current_functional][element][potential]['cexp_ppl'] = \
                    [float(s) for s in l3[2:]]
                i += 1

                if lines[i].split()[0] == 'NLCC':
                    potentials[current_functional][element][potential]['NLCC'] = {
                        'n_nlcc': int(lines[i].split()[-1]),
                        'r_core': float(lines[i+1].split()[0]),
                        'n_core': float(lines[i+1].split()[1]),
                        'c_core': float(lines[i+1].split()[2])
                    }
                    i += 2

                potentials[current_functional][element][potential]['nprj'] = int(lines[i].split()[0])
                i += 1

                potentials[current_functional][element][potential][
                    'r'] = []  # Radius of non-local part for ang. mom. quantum number l
                potentials[current_functional][element][potential][
                    'nprj_ppnl'] = []  # number of nonlocal projectors for ang mom = l
                potentials[current_functional][element][potential]['hprj_ppnl'] = []  # coeff of nonlocal projectors funcs
                for j in range(potentials[current_functional][element][potential]['nprj']):
                    l = lines[i].split()
                    potentials[current_functional][element][potential]['r'].append(float(l[0]))
                    potentials[current_functional][element][potential]['nprj_ppnl'].append(int(l[1]))
                    for k in range(1, potentials[current_functional][element][potential]['nprj_ppnl'][-1] + 1).__reversed__():
                        l = lines[i].split()
                        potentials[current_functional][element][potential]['hprj_ppnl'].extend(
                            [float(s) for s in l[-k:]]
                        )
                        i += 1
                i += 1
    return potentials


# TODO: Setting the default basis set to triple zeta double valence potential (highest accuracy). Check this.
def get_basis_and_potential(species, potential_type='GTH', functional='PBE', basis_type='TZV2P'):

    """
    Given a specie and a potential/basis type, this function accesses the available basis sets and potentials in
    the available "*_POTENTIAL.yaml" files. These files should come with this module distribution, but can be
    updated or remade if needed (see utils.py). Generally, the GTH potentials are used with the GTH basis sets.

    Note: as with most cp2k inputs, the convention is to use all caps, so use type="GTH" instead of "gth"

    Args:
        specie: (list) list of species for which to get the potential/basis strings
        potential_type: (str) the potential type. Default: 'GTH'
        basis_type: (str) the basis set type. Default: 'TZV2P'
        functional: (str) functional type. Default: 'PBE'

            functionals available in CP2K:
                - BLYP
                - BP
                - HCTH120
                - HCTH407
                - PADE
                - PBE
                - PBEsol
                - OLYP

    Returns:
        (dict) of the form {'specie': {'potential': potential, 'basis': basis}...}
    """

    with zopen(os.path.join(MODULE_DIR, '{}_POTENTIALS.yaml'.format(potential_type))) as f:
        potentials = yaml.safe_load(f)

    d = {}
    for specie in species:
        l = list(potentials[functional][specie].keys())
        if len(l) == 1:
            s = l[0].split('-')
            d[specie] = {'potential': l[0],
                         'basis': "{}-GTH-{}".format(basis_type, s[-1])}
        else:
            raise AttributeError('FOUND MORE THAN ONE FUNCTIONAL FOR {} WITH TYPE {}'.format(specie, type),
                                 'AMBIGUITY CANNOT BE HANDLED. MUST MANUALLY SET THE POTENTIAL')
    return d
