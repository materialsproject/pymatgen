"""Module for JDFTx IO module output utils.

Module for JDFTx IO module output utils. Functions kept in this module are here if they are
used by multiple submodules, or if they are anticipated to be used by multiple
submodules in the future.
"""

from __future__ import annotations

from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeAlias

import numpy as np

from pymatgen.electronic_structure.core import Orbital

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray


def check_file_exists(func: Callable) -> Any:
    """Check if file exists.

    Check if file exists (and continue normally) or raise an exception if
    it does not.
    """

    @wraps(func)
    def wrapper(filename: str) -> Any:
        filepath = Path(filename)
        if not filepath.is_file():
            raise OSError(f"'{filename}' file doesn't exist!")
        return func(filename)

    return wrapper


@check_file_exists
def read_file(file_name: str) -> list[str]:
    """
    Read file into a list of str.

    Parameters
    ----------
    filename: Path or str
        name of file to read

    Returns
    -------
    text: list[str]
        list of strings from file
    """
    with open(file_name) as f:
        text = f.readlines()
    f.close()
    return text


def read_outfile_slices(file_name: str) -> list[list[str]]:
    """
    Read slice of out file into a list of str.

    Parameters
    ----------
    filename: Path or str
        name of file to read
    out_slice_idx: int
        index of slice to read from file

    Returns
    -------
    texts: list[list[str]]
        list of out file slices (individual calls of JDFTx)
    """
    _text = read_file(file_name)
    start_lines = get_start_lines(_text, add_end=True)
    texts = []
    for i in range(len(start_lines) - 1):
        text = _text[start_lines[i] : start_lines[i + 1]]
        texts.append(text)
    return texts


def _brkt_list_of_3_to_nparray(line: str) -> NDArray[np.float64]:
    """Return 3x1 numpy array.

    Convert a string of the form "[ x y z ]" to a 3x1 numpy array

    Parameters
    ----------
    line: str
        A string of the form "[ x y z ]"
    """
    return np.array([float(x) for x in line.split()[1:-1]])


def _brkt_list_of_3x3_to_nparray(lines: list[str], i_start: int = 0) -> NDArray[np.float64]:
    """Return 3x3 numpy array.

    Convert a list of strings of the form "[ x y z ]" to a 3x3 numpy array

    Parameters
    ----------
    lines: list[str]
        A list of strings of the form "[ x y z ]"
    i_start: int
        The index of the first line in lines

    Returns
    -------
    out: np.ndarray
        A 3x3 numpy array
    """
    out = np.zeros([3, 3])
    for i in range(3):
        out[i, :] += _brkt_list_of_3_to_nparray(lines[i + i_start])
    return out


# Named "t1" in unmet anticipation of multiple ways that a float would be needed
# to be read following the variable string with a colon.
def get_colon_val(linetext: str, lkey: str) -> float | np.float64 | None:
    """Return float val from '...lkey: val...' in linetext.

    Read a float from an elec minimization line assuming value appears as
    "... lkey value ...".

    Parameters
    ----------
    linetext: str
        A line of text from a JDFTx out file
    lkey: str
        A string that appears before the float value in linetext. Must include
        the colon.
    """
    colon_var = None
    if lkey in linetext:
        val = linetext.split(lkey)[1].strip().split(" ")[0]
        colon_var = np.nan if val == "nan" else float(linetext.split(lkey)[1].strip().split(" ")[0])
    return colon_var


# Temporary alias until the outside modules are merged with the renaming
get_colon_var_t1 = get_colon_val


# This function matches the format of the generic "is_<x>_start_line" functions specific to
# the JOutStructure object initialization, but is not moved to the joutstructure module
# as it is also used in methods for other JDFTx IO modules (e.g. JOutStructures) so it
# is kept here to avoid circular imports.
def is_lowdin_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of Lowdin log message.

    Return True if the line_text is the start of a Lowdin population analysis
    in a JDFTx out file.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file

    Returns
    -------
    is_line: bool
        True if the line_text is the start of a Lowdin population analysis in a
        JDFTx out file
    """
    return "#--- Lowdin population analysis ---" in line_text


# TODO: Figure out if this ever actually gets called, (I think I added it to make JOutStructure user friendly
# but JOutStructure is not intended to be a user-initialized class)
def correct_geom_opt_type(opt_type: str | None) -> str | None:
    """Return recognizable opt_type string.

    Correct the opt_type string to match the JDFTx convention.

    Parameters
    ----------
    opt_type:
        The type of optimization step

    Returns
    -------
    opt_type: str | None
        The corrected type of optimization step
    """
    if opt_type is not None:
        if "lattice" in opt_type.lower():
            opt_type = "LatticeMinimize"
        elif "ionic" in opt_type.lower():
            opt_type = "IonicDynamics" if "dyn" in opt_type.lower() else "IonicMinimize"
        else:
            opt_type = None
    return opt_type


def get_start_lines(
    text: list[str],
    start_key: str = "*************** JDFTx",
    add_end: bool = False,
) -> list[int]:
    """Get start line numbers for JDFTx calculations.

    Get the line numbers corresponding to the beginning of separate JDFTx calculations
    (in case of multiple calculations appending the same out file).

    Args:
        text: output of read_file for out file
    """
    start_lines = []
    i = None
    for i, line in enumerate(text):
        if start_key in line:
            start_lines.append(i)
    if add_end and i is not None:
        start_lines.append(i)
    if i is None:
        raise ValueError("Outfile parser fed an empty file.")
    if not start_lines:
        raise ValueError("No JDFTx calculations found in file.")
    return start_lines


def find_key_first(key_input: str, tempfile: list[str]) -> int | None:
    """Find the first instance of a key in the output file.

    Args:
        key_input (str): Key string to match.
        tempfile (list[str]): Output from readlines() function in read_file method.

    Returns:
        int | None: The index of the first occurrence of the key in the tempfile list, or None if the key is not found.
    """
    key_input = str(key_input)
    line = None
    for i in range(len(tempfile)):
        if key_input in tempfile[i]:
            line = i
            break
    return line


def find_key(key_input: str, tempfile: list[str]) -> int | None:
    """
    Find the last instance of a key in the output file.

    Args:
        key_input (str): Key string to match.
        tempfile (list[str]): Output from readlines() function in read_file method.

    Returns:
        int | None: The index of the last occurrence of the key in the tempfile list, or None if the key is not found.
    """
    key_input = str(key_input)
    line = None
    lines = find_all_key(key_input, tempfile)
    if len(lines):
        line = lines[-1]
    return line


def find_first_range_key(
    key_input: str,
    tempfile: list[str],
    startline: int = 0,
    endline: int = -1,
    skip_pound: bool = False,
) -> list[int]:
    """
    Find all lines that exactly begin with key_input in a range of lines.

    Args:
        key_input (str): Key string to match.
        tempfile (list[str]): Output from readlines() function in read_file method.
        startline (int): Line to start searching from.
        endline (int): Line to stop searching at.
        skip_pound (bool): Whether to skip lines that begin with a pound sign.

    Returns:
        list[int]: List of line numbers where key_input occurs.
    """
    key_input = str(key_input)
    startlen = len(key_input)
    line_list = []

    if endline == -1:
        endline = len(tempfile)
    for i in range(startline, endline):
        line = tempfile[i]
        if skip_pound:
            for _ in range(10):  # repeat to make sure no really weird formatting
                line = line.lstrip()
                line = line.lstrip("#")
        line = line[0:startlen]
        if line == key_input:
            line_list.append(i)
    return line_list


def key_exists(key_input: str, tempfile: list[str]) -> bool:
    """
    Check if key_input exists in tempfile.

    Search through tempfile for key_input. Return True if found, False otherwise.

    Args:
        key_input (str): Key string to match.
        tempfile (list[str]): Output from readlines() function in read_file method.

    Returns:
        bool: True if key_input exists in tempfile, False otherwise.
    """
    line = find_key(key_input, tempfile)
    return line is not None


def find_all_key(key_input: str, tempfile: list[str], startline: int = 0) -> list[int]:
    """Find all lines containing key_input.

    Search through tempfile for all lines containing key_input. Returns a list of line numbers.

    Args:
        key_input (str): Key string to match.
        tempfile (list[str]): Output from readlines() function in read_file method.
        startline (int): Line to start searching from.

    Returns:
        list[int]: List of line numbers where key_input occurs.
    """
    return [i for i in range(startline, len(tempfile)) if key_input in tempfile[i]]


def _init_dict_from_colon_dump_lines(lines: list[str]):
    varsdict = {}
    for line in lines:
        if ":" in line:
            lsplit = line.split(":")
            key = lsplit[0].strip()
            val = lsplit[1].split()[0].strip()
            varsdict[key] = val
    return varsdict


def _parse_bandfile_complex(bandfile_filepath: str | Path) -> NDArray[np.complex64]:
    Dtype: TypeAlias = np.complex64
    token_parser = _complex_token_parser
    return _parse_bandfile_reader(bandfile_filepath, Dtype, token_parser)


def _parse_bandfile_normalized(bandfile_filepath: str | Path) -> NDArray[np.float32]:
    Dtype: TypeAlias = np.float32
    token_parser = _normalized_token_parser
    return _parse_bandfile_reader(bandfile_filepath, Dtype, token_parser)


def _get__from_bandfile_filepath(bandfile_filepath: Path | str, tok_idx: int) -> int:
    """
    Get arbitrary integer from header of bandprojections file.

    Args:
        bandfile_filepath (Path | str): Path to bandprojections file.
        tok_idx (int): Index of token to return.

    Returns:
        int: Integer from header of bandprojections file.
    """
    ret_data = None
    bandfile = read_file(bandfile_filepath)
    for iLine, line in enumerate(bandfile):
        tokens = line.split()
        if iLine == 0:
            ret_data = int(tokens[tok_idx])
            break
    if ret_data is None:
        raise ValueError("Provided an empty file")
    return ret_data


def _get_nstates_from_bandfile_filepath(bandfile_filepath: Path | str) -> int:
    """Get number of states from bandprojections file.

    Get the number of states from the bandprojections file.

    Args:
        bandfile_filepath (Path | str): Path to bandprojections file.

    Returns:
        int: Number of states.
    """
    return _get__from_bandfile_filepath(bandfile_filepath, 0)


def _get_nbands_from_bandfile_filepath(bandfile_filepath: Path | str) -> int:
    """Get number of bands from bandprojections file.

    Get the number of bands from the bandprojections file. The output here should match up with `nbands` from the
    JDFTXOutfile object.

    Args:
        bandfile_filepath (Path | str): Path to bandprojections file.

    Returns:
        int: Number of bands.
    """
    return _get__from_bandfile_filepath(bandfile_filepath, 2)


def _get_nproj_from_bandfile_filepath(bandfile_filepath: Path | str) -> int:
    """Get number of projections from bandprojections file.

    Get the number of projections from the bandprojections file.

    Args:
        bandfile_filepath (Path | str): Path to bandprojections file.

    Returns:
        int: Number of projections.
    """
    return _get__from_bandfile_filepath(bandfile_filepath, 4)


def _get_nspecies_from_bandfile_filepath(bandfile_filepath: Path | str) -> int:
    """Get number of species (ion types) from bandprojections file.

    Get the number of species (ion types) from the bandprojections file.

    Args:
        bandfile_filepath (Path | str): Path to bandprojections file.

    Returns:
        int: Number of species.
    """
    return _get__from_bandfile_filepath(bandfile_filepath, 6)


def _parse_bandfile_reader(bandfile_filepath: str | Path, arr_dtype: TypeAlias, token_parser: Callable) -> NDArray:
    nstates = _get_nstates_from_bandfile_filepath(bandfile_filepath)
    nbands = _get_nbands_from_bandfile_filepath(bandfile_filepath)
    nproj = _get_nproj_from_bandfile_filepath(bandfile_filepath)
    nspecies = _get_nspecies_from_bandfile_filepath(bandfile_filepath)
    # Header of length 3, and then each states occupies 1 (header) + nbands lineas
    bandfile = read_file(bandfile_filepath)
    expected_length = 2 + nspecies + (nstates * (1 + nbands))
    if not expected_length == len(bandfile):
        raise RuntimeError("Bandprojections file does not match expected length - ensure no edits have been made.")
    proj_tju: NDArray[arr_dtype] = np.zeros((nstates, nbands, nproj), dtype=arr_dtype)
    for line, text in enumerate(bandfile):
        tokens = text.split()
        if line >= nspecies + 2:
            istate = (line - (nspecies + 2)) // (nbands + 1)
            iband = (line - (nspecies + 2)) - istate * (nbands + 1) - 1
            if iband >= 0 and istate < nstates:
                proj_tju[istate, iband] = np.array(token_parser(tokens))
    return proj_tju


def _complex_token_parser(tokens: list[str]) -> NDArray[np.complex64]:
    out = np.zeros(int(len(tokens) / 2), dtype=np.complex64)
    ftokens = np.array(tokens, dtype=np.float32)
    out += 1j * ftokens[1::2]
    out += ftokens[::2]
    return out


def _normalized_token_parser(tokens: list[str]) -> NDArray[np.float32]:
    normalized_tokens: NDArray[np.float32] = np.array(tokens, dtype=np.float32)
    return normalized_tokens


def get_proj_tju_from_file(bandfile_filepath: Path | str) -> NDArray[np.float32 | np.complex64]:
    """Return projections from file in tju shape.

    Return projections from file in (state, band, proj) shape. Collected in this shape before sabcju shape due to ready
    availability of this shape in the file.

    Args:
        bandfile_filepath (Path | str): Path to bandprojections file.

    Returns:
        np.ndarray: Projections array in shape (state, band, proj).
    """
    is_complex = _is_complex_bandfile_filepath(bandfile_filepath)
    return _parse_bandfile_complex(bandfile_filepath) if is_complex else _parse_bandfile_normalized(bandfile_filepath)


def _parse_kptsfrom_bandprojections_file(bandfile_filepath: str | Path) -> tuple[list[float], list[NDArray]]:
    """Parse kpts from bandprojections file.

    Parse kpts from bandprojections file.

    Args:
        bandfile_filepath (Path | str): Path to bandprojections file.

    Returns:
        tuple[list[float], list[np.ndarray[float]]]: Tuple of k-point weights and k-points
    """
    wk_list: list[float] = []
    k_points_list: list[NDArray] = []
    kpt_lines = []
    with open(bandfile_filepath) as f:
        for line in f:
            if line.startswith("#") and ";" in line:
                _line = line.split(";")[0].lstrip("#")
                kpt_lines.append(_line)
    for line in kpt_lines:
        k_points = line.split("[")[1].split("]")[0].strip().split()
        _k_points_floats: list[float] = [float(v) for v in k_points]
        k_points_list.append(np.array(_k_points_floats))
        wk = float(line.split("]")[1].strip().split()[0])
        wk_list.append(wk)
    return wk_list, k_points_list


def _is_complex_bandfile_filepath(bandfile_filepath: str | Path) -> bool:
    """Determine if bandprojections file is complex.

    Determine if the bandprojections file is complex. Needed before attempting pCOHP analysis.

    Args:
        bandfile_filepath (Path | str): Path to bandprojections file.

    Returns:
        bool: True if the bandprojections file is complex, False otherwise.
    """
    hash_lines = 0
    val = True
    with open(bandfile_filepath) as f:
        for _i, line in enumerate(f):
            if "#" in line:
                hash_lines += 1
                if hash_lines == 2:
                    val = "|projection|^2" not in line
                    break
    f.close()
    return val


# TODO: This is very likely redundant to something in pymatgen - replace with that if possible.
orb_ref_list = [
    ["s"],
    ["py", "pz", "px"],
    ["dxy", "dyz", "dz2", "dxz", "dx2-y2"],
    ["fy(3x2-y2)", "fxyz", "fyz2", "fz3", "fxz2", "fz(x2-y2)", "fx(x2-3y2)"],
]
orb_ref_to_o_dict = {
    "s": int(Orbital.s),
    "py": int(Orbital.py),
    "pz": int(Orbital.pz),
    "px": int(Orbital.px),
    "dxy": int(Orbital.dxy),
    "dyz": int(Orbital.dyz),
    "dz2": int(Orbital.dz2),
    "dxz": int(Orbital.dxz),
    "dx2-y2": int(Orbital.dx2),
    # Keep the f-orbitals arbitrary-ish until they get designated names in pymatgen.
    orb_ref_list[-1][0]: int(Orbital.f_3),
    orb_ref_list[-1][1]: int(Orbital.f_2),
    orb_ref_list[-1][2]: int(Orbital.f_1),
    orb_ref_list[-1][3]: int(Orbital.f0),
    orb_ref_list[-1][4]: int(Orbital.f1),
    orb_ref_list[-1][5]: int(Orbital.f2),
}


def _get_atom_orb_labels_map_dict(bandfile_filepath: Path) -> dict[str, list[str]]:
    """
    Return a dictionary mapping each atom symbol to pymatgen-compatible orbital projection string representations.

    Identical to _get_atom_orb_labels_ref_dict, but doesn't include the numbers in the labels.



    Args:
        bandfile_filepath (str | Path): The path to the bandfile.

    Returns:
        dict[str, list[str]]: A dictionary mapping each atom symbol to all atomic orbital projection string
        representations.
    """
    bandfile = read_file(bandfile_filepath)
    labels_dict: dict[str, list[str]] = {}

    for i, line in enumerate(bandfile):
        if i > 1:
            if "#" in line:
                break
            lsplit = line.strip().split()
            sym = lsplit[0]
            labels_dict[sym] = []
            lmax = int(lsplit[3])
            # Would prefer to use "l" rather than "L" here (as uppercase "L" means something else entirely) but
            # pr*-c*mm*t thinks "l" is an ambiguous variable name.
            for L in range(lmax + 1):
                mls = orb_ref_list[L]
                nshells = int(lsplit[4 + L])
                for _n in range(nshells):
                    if nshells > 1:
                        for ml in mls:
                            labels_dict[sym].append(f"{ml}")
                    else:
                        labels_dict[sym] += mls
    return labels_dict


def _get_atom_orb_labels_ref_dict(bandfile_filepath: Path) -> dict[str, list[str]]:
    """
    Return a dictionary mapping each atom symbol to all atomic orbital projection string representations.

    Return a dictionary mapping each atom symbol to all atomic orbital projection string representations.
    For example:
    {
        "H": ["s"],
        "O": ["s", "px", "py", "pz", "dxy", "dxz", "dyz", "dx2y2", "dz2"],
        "Pt": ["0s", "1s", "0px", "0py", "0pz", "1px", "1py", "1pz", "dxy", "dxz", "dyz", "dx2y2", "dz2",
               "fx3-3xy2", "fyx2-yz2", "fxz2", "fz3", "fyz2", "fxyz", "f3yx2-y3"]
    }
    where the numbers are needed when using pseudopotentials with multiple valence shells of the same angular momentum
    and are NOT REPRESENTATIVE OF THE TRUE PRINCIPAL QUANTUM NUMBER.

    Args:
        bandfile_filepath (str | Path): The path to the bandfile.

    Returns:
        dict[str, list[str]]: A dictionary mapping each atom symbol to all atomic orbital projection string
        representations.
    """
    bandfile = read_file(bandfile_filepath)
    labels_dict: dict[str, list[str]] = {}

    for i, line in enumerate(bandfile):
        if i > 1:
            if "#" in line:
                break
            lsplit = line.strip().split()
            sym = lsplit[0]
            labels_dict[sym] = []
            lmax = int(lsplit[3])
            # Would prefer to use "l" rather than "L" here (as uppercase "L" means something else entirely) but
            # pr*-c*mm*t thinks "l" is an ambiguous variable name.
            for L in range(lmax + 1):
                mls = orb_ref_list[L]
                nshells = int(lsplit[4 + L])
                for n in range(nshells):
                    if nshells > 1:
                        for ml in mls:
                            labels_dict[sym].append(f"{n}{ml}")
                    else:
                        labels_dict[sym] += mls
    return labels_dict


def _get_atom_count_list(bandfile_filepath: Path) -> list[tuple[str, int]]:
    """
    Return a list of tuples of atom symbols and counts.

    Return a list of tuples of atom symbols and counts. This is superior to a dictionary as it maintains the order of
    the atoms in the bandfile.

    Args:
        bandfile_filepath (str | Path): The path to the bandfile.

    Returns:
        list[tuple[str, int]]: A list of tuples of atom symbols and counts.
    """
    bandfile = read_file(bandfile_filepath)
    atom_count_list = []

    for i, line in enumerate(bandfile):
        if i > 1:
            if "#" in line:
                break
            lsplit = line.strip().split()
            sym = lsplit[0].strip()
            count = int(lsplit[1].strip())
            atom_count_list.append((sym, count))
    return atom_count_list


def _get_orb_label_list_expected_len(labels_dict: dict[str, list[str]], atom_count_list: list[tuple[str, int]]) -> int:
    """
    Return the expected length of the atomic orbital projection string representation list.

    Return the expected length of the atomic orbital projection string representation list.

    Args:
        labels_dict (dict[str, list[str]]): A dictionary mapping each atom symbol to all atomic orbital projection
        string representations.
        atom_count_list (list[tuple[str, int]]): A list of tuples of atom symbols and counts.

    Returns:
        int: The expected length of the atomic orbital projection string representation list.
    """
    expected_len = 0
    for ion_tuple in atom_count_list:
        ion = ion_tuple[0]
        count = ion_tuple[1]
        orbs = labels_dict[ion]
        expected_len += count * len(orbs)
    return expected_len


def _get_orb_label(ion: str, idx: int, orb: str) -> str:
    """
    Return the string representation for an orbital projection.

    Return the string representation for an orbital projection.

    Args:
        ion (str): The symbol of the atom.
        idx (int): The index of the atom.
        orb (str): The atomic orbital projection string representation.

    Returns:
        str: The atomic orbital projection string representation for the atom.
    """
    return f"{ion}#{idx + 1}({orb})"


def _get_u_to_oa_map(bandfile_filepath: Path) -> list[tuple[int, int]]:
    """
    Return a list, where the u'th element is a tuple of the atomic orbital index and the ion index.

    Args:
        bandfile_filepath (str | Path): The path to the bandfile.

    Returns:
        list[tuple[int, int]]: A list, where the u'th element is a tuple of the atomic orbital index and the ion index.
    """
    map_labels_dict = _get_atom_orb_labels_map_dict(bandfile_filepath)
    atom_count_list = _get_atom_count_list(bandfile_filepath)
    u_to_oa_map = []
    a = 0
    for ion, ion_count in atom_count_list:
        for _i in range(ion_count):
            for orb in map_labels_dict[ion]:
                u_to_oa_map.append((orb_ref_to_o_dict[orb], a))
            a += 1
    return u_to_oa_map


def _get_orb_label_list(bandfile_filepath: Path) -> tuple[str, ...]:
    """
    Return a tuple of all atomic orbital projection string representations.

    Return a tuple of all atomic orbital projection string representations.

    Args:
        bandfile_filepath (str | Path): The path to the bandfile.

    Returns:
        tuple[str]: A list of all atomic orbital projection string representations.
    """
    labels_dict = _get_atom_orb_labels_ref_dict(bandfile_filepath)
    atom_count_list = _get_atom_count_list(bandfile_filepath)
    read_file(bandfile_filepath)
    labels_list: list[str] = []
    for ion_tuple in atom_count_list:
        ion = ion_tuple[0]
        orbs = labels_dict[ion]
        count = ion_tuple[1]
        for i in range(count):
            for orb in orbs:
                labels_list.append(_get_orb_label(ion, i, orb))
    # This is most likely unnecessary, but it is a good check to have.
    if len(labels_list) != _get_orb_label_list_expected_len(labels_dict, atom_count_list):
        raise RuntimeError("Number of atomic orbital projections does not match expected length.")
    return tuple(labels_list)
