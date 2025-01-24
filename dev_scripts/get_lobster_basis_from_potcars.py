"""Script to generate potential basis functions for the LOBSTER input files. """


from monty.io import zopen
from pymatgen.core import SETTINGS
from pymatgen.io.vasp.inputs import Potcar
from pathlib import Path
import yaml
import os

def get_valence(filename: str = "POTCAR", available_basis_functions: list[str] | None = None, occupation_cutoff: float = 0.01) -> str:
    """
    Extracts the valence configuration from a POTCAR file.

    Args:
        filename (str): Path to the POTCAR file.
        available_basis_functions (list[str] | None): List of available basis functions.
        occupation_cutoff (float): Cutoff below which the orbital is considered to be empty.

    Returns:
        str: A string of basis functions.
    """
    nextline = False
    ready = False
    counter = 0
    orbitals: list[list[int]] = []
    occupations: list[float] = []
    linenumber = 0

    with zopen(filename, 'rt') as f:
        for line in f:
            if linenumber == 1:
                valence = float(line.split()[0])
            if not ready:
                if not nextline:
                    if "Atomic configuration" in line:
                        nextline = True
                else:
                    number_atomic_orbitals = int(line.split()[0])
                    ready = True
            else:
                if counter <= number_atomic_orbitals:
                    if counter != 0:
                        orbitals.append([int(line.split()[0]), int(line.split()[1])])
                        occupations.append(float(line.split()[4]))
                    counter += 1
                else:
                    break
            linenumber += 1

    # Reverse the lists for processing
    rev_orb = list(reversed(orbitals))

    rev_occ = list(reversed(occupations))
    basis: list[list[int]] = []

    valence_checker = valence
    for orbital in range(len(rev_orb)):
        if not abs(valence_checker - 0.0) < 0.01:
            if abs(rev_occ[orbital] - 0.0) >= occupation_cutoff:
                valence_checker -= rev_occ[orbital]
                basis.append(rev_orb[orbital])

    basis = list(reversed(basis))

    # Convert orbital numbers to strings
    basis_result: list[str] = []
    for orbital in basis:
        match orbital[1]:
            case 0:
                basis_result.append(f"{orbital[0]}s")
            case 1:
                basis_result.append(f"{orbital[0]}p")
            case 2:
                basis_result.append(f"{orbital[0]}d")
            case 3:
                basis_result.append(f"{orbital[0]}f")

    if available_basis_functions is not None:
        basis_string = " ".join(e for e in sorted(set(basis_result)) if e in available_basis_functions)
    else:
        basis_string = " ".join(sorted(set(basis_result)))

    return basis_string

# Set occupation cutoff
occupation_cutoff = 0.2  # if set to 0, also empty basis functions will be included

# Path to the POTCAR directory
source_potcar = Path(SETTINGS["PMG_VASP_PSP_DIR"]) / "POT_GGA_PAW_PBE_54"

# Get all POTCAR subfolders, excluding GW potcars
potcar_names: list[str] = [o for o in os.listdir(source_potcar)]

# Read in all available basis functions
available_basis_functions: dict[str, list[str]] = {}
with zopen("all_basis_functions_pbeVaspFit2015_lobster.txt") as f:
    data = f.read().split("\n")

for datum in data:
    datum_raw = datum.split(" ")
    if len(datum_raw) > 2:
        available_basis_functions[datum_raw[1]] = datum_raw[2:]

basis_dict: dict[str, dict[str, str]] = {'BASIS': {}}
for potcar in potcar_names:
    potential_type = potcar.split('.')[1]
    element = potential_type.split('_')[0]

    if element in available_basis_functions:
        potcar_object = Potcar.from_file(source_potcar / potcar)
        basis_dict['BASIS'][potcar_object.spec[0]['symbol']] = get_valence(
            filename=source_potcar / potcar,
            available_basis_functions=available_basis_functions[element],
            occupation_cutoff=occupation_cutoff
        )

# Save the basis dictionary to a YAML file
with open('BASIS_PBE_54.yaml', 'w') as outfile:
    yaml.dump(basis_dict, outfile, default_flow_style=False)
