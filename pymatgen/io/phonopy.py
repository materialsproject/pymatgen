"""Module for interfacing with phonopy, see https://atztogo.github.io/phonopy/."""

from __future__ import annotations

import numpy as np
from monty.dev import requires
from monty.serialization import loadfn
from scipy.interpolate import InterpolatedUnivariateSpline

from pymatgen.core import Lattice, Structure
from pymatgen.phonon.bandstructure import PhononBandStructure, PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos, PhononDos
from pymatgen.phonon.gruneisen import GruneisenParameter, GruneisenPhononBandStructureSymmLine
from pymatgen.phonon.thermal_displacements import ThermalDisplacementMatrices
from pymatgen.symmetry.bandstructure import HighSymmKpath

try:
    from phonopy import Phonopy
    from phonopy.file_IO import write_disp_yaml
    from phonopy.structure.atoms import PhonopyAtoms
except ImportError:
    Phonopy = write_disp_yaml = PhonopyAtoms = None


@requires(Phonopy, "phonopy not installed!")
def get_pmg_structure(phonopy_structure: PhonopyAtoms) -> Structure:
    """
    Convert a PhonopyAtoms object to pymatgen Structure object.

    Args:
        phonopy_structure (PhonopyAtoms): A phonopy structure object.
    """
    lattice = phonopy_structure.cell
    frac_coords = phonopy_structure.scaled_positions
    symbols = phonopy_structure.symbols
    masses = phonopy_structure.masses
    mms = getattr(phonopy_structure, "magnetic_moments", None) or [0] * len(symbols)

    return Structure(
        lattice,
        symbols,
        frac_coords,
        site_properties={"phonopy_masses": masses, "magnetic_moments": mms},
    )


@requires(Phonopy, "phonopy not installed!")
def get_phonopy_structure(pmg_structure: Structure) -> PhonopyAtoms:
    """
    Convert a pymatgen Structure object to a PhonopyAtoms object.

    Args:
        pmg_structure (pymatgen Structure): A Pymatgen structure object.
    """
    symbols = [site.specie.symbol for site in pmg_structure]

    return PhonopyAtoms(
        symbols=symbols,
        cell=pmg_structure.lattice.matrix,
        scaled_positions=pmg_structure.frac_coords,
    )


def get_structure_from_dict(d):
    """
    Extracts a structure from the dictionary extracted from the output
    files of phonopy like phonopy.yaml or band.yaml.
    Adds "phonopy_masses" in the site_properties of the structures.
    Compatible with older phonopy versions.
    """
    species = []
    frac_coords = []
    masses = []
    if "points" in d:
        for p in d["points"]:
            species.append(p["symbol"])
            frac_coords.append(p["coordinates"])
            masses.append(p["mass"])
    elif "atoms" in d:
        for p in d["atoms"]:
            species.append(p["symbol"])
            frac_coords.append(p["position"])
            masses.append(p["mass"])
    else:
        raise ValueError("The dict does not contain structural information")

    return Structure(d["lattice"], species, frac_coords, site_properties={"phonopy_masses": masses})


def eigvec_to_eigdispl(v, q, frac_coords, mass):
    r"""
    Converts a single eigenvector to an eigendisplacement in the primitive cell
    according to the formula::

        exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v

    Compared to the modulation option in phonopy, here all the additional
    multiplicative and phase factors are set to 1.

    Args:
        v: the vector that should be converted. A 3D complex numpy array.
        q: the q point in fractional coordinates
        frac_coords: the fractional coordinates of the atom
        mass: the mass of the atom
    """
    c = np.exp(2j * np.pi * np.dot(frac_coords, q)) / np.sqrt(mass)

    return c * v


def get_ph_bs_symm_line_from_dict(bands_dict, has_nac=False, labels_dict=None):
    r"""
    Creates a pymatgen PhononBandStructure object from the dictionary
    extracted by the band.yaml file produced by phonopy. The labels
    will be extracted from the dictionary, if present. If the 'eigenvector'
    key is found the eigendisplacements will be calculated according to the
    formula::

        exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v

    and added to the object.

    Args:
        bands_dict: the dictionary extracted from the band.yaml file
        has_nac: True if the data have been obtained with the option
            --nac option. Default False.
        labels_dict: dict that links a qpoint in frac coords to a label.
            Its value will replace the data contained in the band.yaml.
    """
    structure = get_structure_from_dict(bands_dict)

    qpts = []
    frequencies = []
    eigendisplacements = []
    phonopy_labels_dict = {}
    for p in bands_dict["phonon"]:
        q = p["q-position"]
        qpts.append(q)
        bands = []
        eig_q = []
        for b in p["band"]:
            bands.append(b["frequency"])
            if "eigenvector" in b:
                eig_b = []
                for i, eig_a in enumerate(b["eigenvector"]):
                    v = np.zeros(3, complex)
                    for x in range(3):
                        v[x] = eig_a[x][0] + eig_a[x][1] * 1j
                    eig_b.append(
                        eigvec_to_eigdispl(
                            v,
                            q,
                            structure[i].frac_coords,
                            structure.site_properties["phonopy_masses"][i],
                        )
                    )
                eig_q.append(eig_b)
        frequencies.append(bands)
        if "label" in p:
            phonopy_labels_dict[p["label"]] = p["q-position"]
        if eig_q:
            eigendisplacements.append(eig_q)

    qpts = np.array(qpts)
    # transpose to match the convention in PhononBandStructure
    frequencies = np.transpose(frequencies)
    if eigendisplacements:
        eigendisplacements = np.transpose(eigendisplacements, (1, 0, 2, 3))

    rec_latt = Lattice(bands_dict["reciprocal_lattice"])

    labels_dict = labels_dict or phonopy_labels_dict

    return PhononBandStructureSymmLine(
        qpts,
        frequencies,
        rec_latt,
        has_nac=has_nac,
        labels_dict=labels_dict,
        structure=structure,
        eigendisplacements=eigendisplacements,
    )


def get_ph_bs_symm_line(bands_path, has_nac=False, labels_dict=None):
    r"""
    Creates a pymatgen PhononBandStructure from a band.yaml file.
    The labels will be extracted from the dictionary, if present.
    If the 'eigenvector'  key is found the eigendisplacements will be
    calculated according to the formula:
    \\exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
     and added to the object.

    Args:
        bands_path: path to the band.yaml file
        has_nac: True if the data have been obtained with the option
            --nac option. Default False.
        labels_dict: dict that links a qpoint in frac coords to a label.
    """
    return get_ph_bs_symm_line_from_dict(loadfn(bands_path), has_nac, labels_dict)


def get_ph_dos(total_dos_path):
    """
    Creates a pymatgen PhononDos from a total_dos.dat file.

    Args:
        total_dos_path: path to the total_dos.dat file.
    """
    a = np.loadtxt(total_dos_path)
    return PhononDos(a[:, 0], a[:, 1])


def get_complete_ph_dos(partial_dos_path, phonopy_yaml_path):
    """
    Creates a pymatgen CompletePhononDos from a partial_dos.dat and
    phonopy.yaml files.
    The second is produced when generating a Dos and is needed to extract
    the structure.

    Args:
        partial_dos_path: path to the partial_dos.dat file.
        phonopy_yaml_path: path to the phonopy.yaml file.
    """
    a = np.loadtxt(partial_dos_path).transpose()
    d = loadfn(phonopy_yaml_path)

    structure = get_structure_from_dict(d["primitive_cell"])

    total_dos = PhononDos(a[0], a[1:].sum(axis=0))

    pdoss = {}
    for site, pdos in zip(structure, a[1:]):
        pdoss[site] = pdos.tolist()

    return CompletePhononDos(structure, total_dos, pdoss)


@requires(Phonopy, "phonopy not installed!")
def get_displaced_structures(pmg_structure, atom_disp=0.01, supercell_matrix=None, yaml_fname=None, **kwargs):
    r"""
    Generate a set of symmetrically inequivalent displaced structures for
    phonon calculations.

    Args:
        pmg_structure (Structure): A pymatgen structure object.
        atom_disp (float): Atomic displacement. Default is 0.01 $\\AA$.
        supercell_matrix (3x3 array): Scaling matrix for supercell.
        yaml_fname (str): If not None, it represents the full path to
            the outputting displacement yaml file, e.g. disp.yaml.
        **kwargs: Parameters used in Phonopy.generate_displacement method.

    Return:
        A list of symmetrically inequivalent structures with displacements, in
        which the first element is the perfect supercell structure.
    """
    is_plusminus = kwargs.get("is_plusminus", "auto")
    is_diagonal = kwargs.get("is_diagonal", True)
    is_trigonal = kwargs.get("is_trigonal", False)

    ph_structure = get_phonopy_structure(pmg_structure)

    if supercell_matrix is None:
        supercell_matrix = np.eye(3) * np.array((1, 1, 1))

    phonon = Phonopy(unitcell=ph_structure, supercell_matrix=supercell_matrix)
    phonon.generate_displacements(
        distance=atom_disp,
        is_plusminus=is_plusminus,
        is_diagonal=is_diagonal,
        is_trigonal=is_trigonal,
    )

    if yaml_fname is not None:
        displacements = phonon.get_displacements()
        write_disp_yaml(
            displacements=displacements,
            supercell=phonon.get_supercell(),
            filename=yaml_fname,
        )

    # Supercell structures with displacement
    disp_supercells = phonon.get_supercells_with_displacements()
    # Perfect supercell structure
    init_supercell = phonon.get_supercell()
    # Structure list to be returned
    structure_list = [get_pmg_structure(init_supercell)]

    for c in disp_supercells:
        if c is not None:
            structure_list.append(get_pmg_structure(c))

    return structure_list


@requires(Phonopy, "phonopy is required to calculate phonon density of states")
def get_phonon_dos_from_fc(
    structure: Structure,
    supercell_matrix: np.ndarray,
    force_constants: np.ndarray,
    mesh_density: float = 100.0,
    num_dos_steps: int = 200,
    **kwargs,
) -> CompletePhononDos:
    """
    Get a projected phonon density of states from phonopy force constants.

    Args:
        structure: A structure.
        supercell_matrix: The supercell matrix used to generate the force
            constants.
        force_constants: The force constants in phonopy format.
        mesh_density: The density of the q-point mesh. See the docstring
            for the ``mesh`` argument in Phonopy.init_mesh() for more details.
        num_dos_steps: Number of frequency steps in the energy grid.
        **kwargs: Additional kwargs passed to the Phonopy constructor.

    Returns:
        The density of states.
    """
    structure_phonopy = get_phonopy_structure(structure)
    phonon = Phonopy(structure_phonopy, supercell_matrix=supercell_matrix, **kwargs)
    phonon.set_force_constants(force_constants)
    phonon.run_mesh(
        mesh_density,
        is_mesh_symmetry=False,
        with_eigenvectors=True,
        is_gamma_center=True,
    )

    # get min, max, step frequency
    frequencies = phonon.get_mesh_dict()["frequencies"]
    freq_min = frequencies.min()
    freq_max = frequencies.max()
    freq_pitch = (freq_max - freq_min) / num_dos_steps

    phonon.run_projected_dos(freq_min=freq_min, freq_max=freq_max, freq_pitch=freq_pitch)

    dos_raw = phonon.projected_dos.get_partial_dos()
    pdoss = dict(zip(structure, dos_raw[1]))

    total_dos = PhononDos(dos_raw[0], dos_raw[1].sum(axis=0))
    return CompletePhononDos(structure, total_dos, pdoss)


@requires(Phonopy, "phonopy is required to calculate phonon band structures")
def get_phonon_band_structure_from_fc(
    structure: Structure,
    supercell_matrix: np.ndarray,
    force_constants: np.ndarray,
    mesh_density: float = 100.0,
    **kwargs,
) -> PhononBandStructure:
    """
    Get a uniform phonon band structure from phonopy force constants.

    Args:
        structure: A structure.
        supercell_matrix: The supercell matrix used to generate the force
            constants.
        force_constants: The force constants in phonopy format.
        mesh_density: The density of the q-point mesh. See the docstring
            for the ``mesh`` argument in Phonopy.init_mesh() for more details.
        **kwargs: Additional kwargs passed to the Phonopy constructor.

    Returns:
        The uniform phonon band structure.
    """
    structure_phonopy = get_phonopy_structure(structure)
    phonon = Phonopy(structure_phonopy, supercell_matrix=supercell_matrix, **kwargs)
    phonon.set_force_constants(force_constants)
    phonon.run_mesh(mesh_density, is_mesh_symmetry=False, is_gamma_center=True)
    mesh = phonon.get_mesh_dict()

    return PhononBandStructure(mesh["qpoints"], mesh["frequencies"], structure.lattice)


@requires(Phonopy, "phonopy is required to calculate phonon band structures")
def get_phonon_band_structure_symm_line_from_fc(
    structure: Structure,
    supercell_matrix: np.ndarray,
    force_constants: np.ndarray,
    line_density: float = 20.0,
    symprec: float = 0.01,
    **kwargs,
) -> PhononBandStructureSymmLine:
    """
    Get a phonon band structure along a high symmetry path from phonopy force
    constants.

    Args:
        structure: A structure.
        supercell_matrix: The supercell matrix used to generate the force
            constants.
        force_constants: The force constants in phonopy format.
        line_density: The density along the high symmetry path.
        symprec: Symmetry precision passed to phonopy and used for determining
            the band structure path.
        **kwargs: Additional kwargs passed to the Phonopy constructor.

    Returns:
        The line mode band structure.
    """
    structure_phonopy = get_phonopy_structure(structure)
    phonon = Phonopy(structure_phonopy, supercell_matrix=supercell_matrix, symprec=symprec, **kwargs)
    phonon.set_force_constants(force_constants)

    kpath = HighSymmKpath(structure, symprec=symprec)

    kpoints, labels = kpath.get_kpoints(line_density=line_density, coords_are_cartesian=False)

    phonon.run_qpoints(kpoints)
    frequencies = phonon.qpoints.get_frequencies().T

    labels_dict = {a: k for a, k in zip(labels, kpoints) if a != ""}

    return PhononBandStructureSymmLine(kpoints, frequencies, structure.lattice, labels_dict=labels_dict)


def get_gruneisenparameter(gruneisen_path, structure=None, structure_path=None) -> GruneisenParameter:
    """
    Get Gruneisen object from gruneisen.yaml file, as obtained from phonopy (Frequencies in THz!).
    The order is structure > structure path > structure from gruneisen dict.
    Newer versions of phonopy include the structure in the yaml file,
    the structure/structure_path is kept for compatibility.

    Args:
        gruneisen_path: Path to gruneisen.yaml file (frequencies have to be in THz!)
        structure: pymatgen Structure object
        structure_path: path to structure in a file (e.g., POSCAR)

    Returns: GruneisenParameter object

    """
    gruneisen_dict = loadfn(gruneisen_path)

    if structure_path and structure is None:
        structure = Structure.from_file(structure_path)
    else:
        try:
            structure = get_structure_from_dict(gruneisen_dict)
        except ValueError as exc:
            raise ValueError("Please provide a structure or structure path") from exc

    qpts, multiplicities, frequencies, gruneisen = ([] for _ in range(4))
    phonopy_labels_dict = {}

    for p in gruneisen_dict["phonon"]:
        q = p["q-position"]
        qpts.append(q)
        m = p["multiplicity"] if "multiplicity" in p else 1
        multiplicities.append(m)
        bands, gruneisenband = [], []
        for b in p["band"]:
            bands.append(b["frequency"])
            if "gruneisen" in b:
                gruneisenband.append(b["gruneisen"])
        frequencies.append(bands)
        gruneisen.append(gruneisenband)
        if "label" in p:
            phonopy_labels_dict[p["label"]] = p["q-position"]

    qpts_np = np.array(qpts)
    multiplicities_np = np.array(multiplicities)
    # transpose to match the convention in PhononBandStructure
    frequencies_np = np.transpose(frequencies)
    gruneisen_np = np.transpose(gruneisen)

    return GruneisenParameter(
        gruneisen=gruneisen_np,
        qpoints=qpts_np,
        multiplicities=multiplicities_np,
        frequencies=frequencies_np,
        structure=structure,
    )


def get_gs_ph_bs_symm_line_from_dict(
    gruneisen_dict, structure=None, structure_path=None, labels_dict=None, fit=False
) -> GruneisenPhononBandStructureSymmLine:
    r"""
    Creates a pymatgen GruneisenPhononBandStructure object from the dictionary
    extracted by the gruneisen.yaml file produced by phonopy. The labels
    will be extracted from the dictionary, if present. If the 'eigenvector'
    key is found the eigendisplacements will be calculated according to the
    formula::

        exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v

    and added to the object. A fit algorithm can be used to replace diverging
    Gruneisen values close to gamma.

    Args:
        gruneisen_dict (dict): the dictionary extracted from the gruneisen.yaml file
        structure (Structure): pymatgen structure object
        structure_path: path to structure file
        labels_dict (dict): dict that links a qpoint in frac coords to a label.
            Its value will replace the data contained in the band.yaml.
        fit (bool): Substitute Grueneisen parameters close to the gamma point
            with points obtained from a fit to a spline if the derivate from
            a smooth curve (i.e. if the slope changes by more than 200% in the
            range of 10% around the gamma point).
            These derivations occur because of very small frequencies
            (and therefore numerical inaccuracies) close to gamma.
    """
    if structure_path and structure is None:
        structure = Structure.from_file(structure_path)
    else:
        try:
            structure = get_structure_from_dict(gruneisen_dict)
        except ValueError as exc:
            raise ValueError("Please provide a structure or structure path") from exc

    q_points, frequencies, gruneisen_params = [], [], []
    phonopy_labels_dict: dict[str, dict[str, str]] = {}

    if fit:
        for pa in gruneisen_dict["path"]:
            phonon = pa["phonon"]  # This is a list
            start = pa["phonon"][0]
            end = pa["phonon"][-1]

            if start["q-position"] == [0, 0, 0]:  # Gamma at start of band
                qpts_temp, frequencies_temp = [], []
                gruneisen_temp: list[list[float]] = []
                distance: list[float] = []
                for i in range(pa["nqpoint"]):
                    bands = []
                    gruneisen_band: list[float] = []
                    for b in phonon[pa["nqpoint"] - i - 1]["band"]:
                        bands.append(b["frequency"])
                        # Fraction of leftover points in current band
                        gruen = _extrapolate_grun(b, distance, gruneisen_temp, gruneisen_band, i, pa)
                        gruneisen_band.append(gruen)
                    q = phonon[pa["nqpoint"] - i - 1]["q-position"]
                    qpts_temp.append(q)
                    d = phonon[pa["nqpoint"] - i - 1]["distance"]
                    distance.append(d)
                    frequencies_temp.append(bands)
                    gruneisen_temp.append(gruneisen_band)
                    if "label" in phonon[pa["nqpoint"] - i - 1]:
                        phonopy_labels_dict[phonon[pa["nqpoint"] - i - 1]]["label"] = phonon[pa["nqpoint"] - i - 1][
                            "q-position"
                        ]

                q_points.extend(list(reversed(qpts_temp)))
                frequencies.extend(list(reversed(frequencies_temp)))
                gruneisen_params.extend(list(reversed(gruneisen_temp)))

            elif end["q-position"] == [0, 0, 0]:  # Gamma at end of band
                distance = []
                for i in range(pa["nqpoint"]):
                    bands, gruneisen_band = [], []
                    for b in phonon[i]["band"]:
                        bands.append(b["frequency"])
                        gruen = _extrapolate_grun(b, distance, gruneisen_params, gruneisen_band, i, pa)
                        gruneisen_band.append(gruen)
                    q = phonon[i]["q-position"]
                    q_points.append(q)
                    d = phonon[i]["distance"]
                    distance.append(d)
                    frequencies.append(bands)
                    gruneisen_params.append(gruneisen_band)
                    if "label" in phonon[i]:
                        phonopy_labels_dict[phonon[i]["label"]] = phonon[i]["q-position"]

            else:  # No Gamma in band
                distance = []
                for i in range(pa["nqpoint"]):
                    bands, gruneisen_band = [], []
                    for b in phonon[i]["band"]:
                        bands.append(b["frequency"])
                        gruneisen_band.append(b["gruneisen"])
                    q = phonon[i]["q-position"]
                    q_points.append(q)
                    d = phonon[i]["distance"]
                    distance.append(d)
                    frequencies.append(bands)
                    gruneisen_params.append(gruneisen_band)
                    if "label" in phonon[i]:
                        phonopy_labels_dict[phonon[i]["label"]] = phonon[i]["q-position"]

    else:
        for pa in gruneisen_dict["path"]:
            for p in pa["phonon"]:
                q = p["q-position"]
                q_points.append(q)
                bands, gruneisen_bands = [], []
                for b in p["band"]:
                    bands.append(b["frequency"])
                    gruneisen_bands.append(b["gruneisen"])
                frequencies.append(bands)
                gruneisen_params.append(gruneisen_bands)
                if "label" in p:
                    phonopy_labels_dict[p["label"]] = p["q-position"]

    rec_latt = structure.lattice.reciprocal_lattice
    labels_dict = labels_dict or phonopy_labels_dict
    return GruneisenPhononBandStructureSymmLine(
        qpoints=np.array(q_points),
        # transpose to match the convention in PhononBandStructure
        frequencies=np.transpose(frequencies),
        gruneisenparameters=np.transpose(gruneisen_params),
        lattice=rec_latt,
        labels_dict=labels_dict,
        structure=structure,
        eigendisplacements=None,
    )


def _extrapolate_grun(b, distance, gruneisenparameter, gruneisenband, i, pa):
    leftover_fraction = (pa["nqpoint"] - i - 1) / pa["nqpoint"]
    if leftover_fraction < 0.1:
        diff = abs(b["gruneisen"] - gruneisenparameter[-1][len(gruneisenband)]) / abs(
            gruneisenparameter[-2][len(gruneisenband)] - gruneisenparameter[-1][len(gruneisenband)]
        )
        if diff > 2:
            x = list(range(len(distance)))
            y = [i[len(gruneisenband)] for i in gruneisenparameter]
            y = y[-len(x) :]  # Only elements of current band
            extrapolator = InterpolatedUnivariateSpline(x, y, k=5)
            g_extrapolated = extrapolator(len(distance))
            gruen = float(g_extrapolated)
        else:
            gruen = b["gruneisen"]
    else:
        gruen = b["gruneisen"]
    return gruen


def get_gruneisen_ph_bs_symm_line(gruneisen_path, structure=None, structure_path=None, labels_dict=None, fit=False):
    r"""
    Creates a pymatgen GruneisenPhononBandStructure from a band.yaml file.
    The labels will be extracted from the dictionary, if present.
    If the 'eigenvector' key is found the eigendisplacements will be
    calculated according to the formula:
    \\exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
     and added to the object.

    Args:
        gruneisen_path: path to the band.yaml file
        structure: pymaten Structure object
        structure_path: path to a structure file (e.g., POSCAR)
        labels_dict: dict that links a qpoint in frac coords to a label.
        fit: Substitute Grueneisen parameters close to the gamma point
            with points obtained from a fit to a spline if the derivate from
            a smooth curve (i.e. if the slope changes by more than 200% in the
            range of 10% around the gamma point).
            These derivations occur because of very small frequencies
            (and therefore numerical inaccuracies) close to gamma.
    """
    return get_gs_ph_bs_symm_line_from_dict(loadfn(gruneisen_path), structure, structure_path, labels_dict, fit)


def get_thermal_displacement_matrices(
    thermal_displacements_yaml="thermal_displacement_matrices.yaml", structure_path="POSCAR"
):
    """
    Function to read "thermal_displacement_matrices.yaml" from phonopy and return a list of
    ThermalDisplacementMatrices objects
    Args:
        thermal_displacements_yaml: path to thermal_displacement_matrices.yaml
        structure_path: path to POSCAR.

    Returns:
    """
    thermal_displacements_dict = loadfn(thermal_displacements_yaml)

    if not structure_path:
        raise ValueError("Please provide a structure_path")
    structure = Structure.from_file(structure_path)

    thermal_displacement_objects_list = []
    for matrix in thermal_displacements_dict["thermal_displacement_matrices"]:
        thermal_displacement_objects_list.append(
            ThermalDisplacementMatrices(
                thermal_displacement_matrix_cart=matrix["displacement_matrices"],
                temperature=matrix["temperature"],
                structure=structure,
                thermal_displacement_matrix_cif=matrix["displacement_matrices_cif"],
            )
        )

    return thermal_displacement_objects_list
