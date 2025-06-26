"""Module for interfacing with phonopy, see https://atztogo.github.io/phonopy/."""

from __future__ import annotations

import typing

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

if typing.TYPE_CHECKING:
    from pymatgen.core.structure import IStructure

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
    magmoms = getattr(phonopy_structure, "magnetic_moments", [0] * len(symbols))
    site_props = {"phonopy_masses": phonopy_structure.masses, "magmom": magmoms}

    return Structure(lattice, symbols, frac_coords, site_properties=site_props)


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
        magnetic_moments=pmg_structure.site_properties.get("magmom"),
    )


def get_structure_from_dict(dct) -> Structure:
    """Extracts a structure from the dictionary extracted from the output
    files of phonopy like phonopy.yaml or band.yaml.
    Adds "phonopy_masses" in the site_properties of the structures.
    Compatible with older phonopy versions.
    """
    species = []
    frac_coords = []
    masses = []
    if "points" in dct:
        for pt in dct["points"]:
            species.append(pt["symbol"])
            frac_coords.append(pt["coordinates"])
            masses.append(pt["mass"])
    elif "atoms" in dct:
        for pt in dct["atoms"]:
            species.append(pt["symbol"])
            frac_coords.append(pt["position"])
            masses.append(pt["mass"])
    else:
        raise ValueError("The dict does not contain structural information")

    return Structure(dct["lattice"], species, frac_coords, site_properties={"phonopy_masses": masses})


def eigvec_to_eigdispl(eig_vec, q, frac_coords, mass):
    """
    Converts a single eigenvector to an eigendisplacement in the primitive cell
    according to the formula:

        exp(2*pi*i*(frac_coords dot q) / sqrt(mass) * v

    Compared to the modulation option in phonopy, here all the additional
    multiplicative and phase factors are set to 1.

    Args:
        v: the vector that should be converted. A 3D complex numpy array.
        q: the q point in fractional coordinates
        frac_coords: the fractional coordinates of the atom
        mass: the mass of the atom
    """
    c = np.exp(2j * np.pi * np.dot(frac_coords, q)) / np.sqrt(mass)

    return c * eig_vec


def get_ph_bs_symm_line_from_dict(bands_dict, has_nac=False, labels_dict=None) -> PhononBandStructureSymmLine:
    r"""Create a pymatgen PhononBandStructure object from the dictionary
    extracted by the band.yaml file produced by phonopy. The labels
    will be extracted from the dictionary, if present. If the 'eigenvector'
    key is found the eigendisplacements will be calculated according to the
    formula:

        exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v

    and added to the object.

    Args:
        bands_dict: the dictionary extracted from the band.yaml file
        has_nac: True if the data have been obtained with the option
            --nac option. Default False.
        labels_dict: dict that links a qpoint in frac coords to a label.
            Its value will replace the data contained in the band.yaml.

    Returns:
        PhononBandStructure: the phonon band structure
    """
    structure = get_structure_from_dict(bands_dict)

    q_pts = []
    frequencies = []
    eigen_displacements = []
    phonopy_labels_dict = {}
    for phonon in bands_dict["phonon"]:
        q_pos = phonon["q-position"]
        q_pts.append(q_pos)
        bands = []
        eig_q = []
        for band in phonon["band"]:
            bands.append(band["frequency"])
            if "eigenvector" in band:
                eig_b = []
                for idx, eig_a in enumerate(band["eigenvector"]):
                    eig_vec = np.zeros(3, complex)
                    for x in range(3):
                        eig_vec[x] = eig_a[x][0] + eig_a[x][1] * 1j
                    eig_b.append(
                        eigvec_to_eigdispl(
                            eig_vec,
                            q_pos,
                            structure[idx].frac_coords,
                            structure.site_properties["phonopy_masses"][idx],
                        )
                    )
                eig_q.append(eig_b)
        frequencies.append(bands)
        if "label" in phonon:
            phonopy_labels_dict[phonon["label"]] = phonon["q-position"]
        if eig_q:
            eigen_displacements.append(eig_q)

    q_pts = np.array(q_pts)  # type:ignore[assignment]
    # transpose to match the convention in PhononBandStructure
    frequencies = np.transpose(frequencies)  # type:ignore[assignment]
    if eigen_displacements:
        eigen_displacements = np.transpose(eigen_displacements, (1, 0, 2, 3))  # type:ignore[assignment]

    rec_lattice = Lattice(bands_dict["reciprocal_lattice"])

    labels_dict = labels_dict or phonopy_labels_dict

    return PhononBandStructureSymmLine(
        q_pts,
        frequencies,
        rec_lattice,
        has_nac=has_nac,
        labels_dict=labels_dict,
        structure=structure,
        eigendisplacements=eigen_displacements,
    )


def get_ph_bs_symm_line(bands_path, has_nac=False, labels_dict=None):
    r"""Create a pymatgen PhononBandStructure from a band.yaml file.
    The labels will be extracted from the dictionary, if present.
    If the 'eigenvector'  key is found the eigendisplacements will be
    calculated according to the formula:
    \\exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
     and added to the object.

    Args:
        bands_path: path to the band.yaml file
        has_nac: True if the data have been obtained with the option
            --nac option. Default False.
        labels_dict: dict that links a q-point in frac coords to a label.
    """
    return get_ph_bs_symm_line_from_dict(loadfn(bands_path), has_nac, labels_dict)


def get_ph_dos(total_dos_path):
    """Create a pymatgen PhononDos from a total_dos.dat file.

    Args:
        total_dos_path: path to the total_dos.dat file.
    """
    arr = np.loadtxt(total_dos_path)
    return PhononDos(arr[:, 0], arr[:, 1])


def get_complete_ph_dos(partial_dos_path, phonopy_yaml_path):
    """Create a pymatgen CompletePhononDos from a partial_dos.dat and
    phonopy.yaml files.
    The second is produced when generating a Dos and is needed to extract
    the structure.

    Args:
        partial_dos_path: path to the partial_dos.dat file.
        phonopy_yaml_path: path to the phonopy.yaml file.
    """
    arr = np.loadtxt(partial_dos_path).transpose()
    dct = loadfn(phonopy_yaml_path)

    structure = get_structure_from_dict(dct["primitive_cell"])

    total_dos = PhononDos(arr[0], arr[1:].sum(axis=0))

    partial_doses = {}
    for site, p_dos in zip(structure, arr[1:], strict=True):
        partial_doses[site] = p_dos.tolist()

    return CompletePhononDos(structure, total_dos, partial_doses)


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
            the outputting displacement YAML file, e.g. disp.yaml.
        **kwargs: Parameters used in Phonopy.generate_displacement method.

    Returns:
        A list of symmetrically inequivalent structures with displacements, in
        which the first element is the perfect supercell structure.
    """
    is_plus_minus = kwargs.get("is_plusminus", "auto")
    is_diagonal = kwargs.get("is_diagonal", True)
    is_trigonal = kwargs.get("is_trigonal", False)

    ph_structure = get_phonopy_structure(pmg_structure)

    if supercell_matrix is None:
        supercell_matrix = np.eye(3) * np.array((1, 1, 1))

    phonon = Phonopy(unitcell=ph_structure, supercell_matrix=supercell_matrix)
    phonon.generate_displacements(
        distance=atom_disp,
        is_plusminus=is_plus_minus,
        is_diagonal=is_diagonal,
        is_trigonal=is_trigonal,
    )

    if yaml_fname is not None:
        displacements = phonon.displacements
        write_disp_yaml(
            displacements=displacements,
            supercell=phonon.supercell,
            filename=yaml_fname,
        )

    # Supercell structures with displacement
    disp_supercells = phonon.supercells_with_displacements
    # Perfect supercell structure
    init_supercell = phonon.supercell
    # Structure list to be returned
    structure_list = [get_pmg_structure(init_supercell)]

    for cell in disp_supercells:
        if cell is not None:
            structure_list.append(get_pmg_structure(cell))

    return structure_list


@requires(Phonopy, "phonopy is required to calculate phonon density of states")
def get_phonon_dos_from_fc(
    structure: Structure | IStructure,
    supercell_matrix: np.ndarray,
    force_constants: np.ndarray,
    mesh_density: float = 100.0,
    num_dos_steps: int = 200,
    **kwargs,
) -> CompletePhononDos:
    """Get a projected phonon density of states from phonopy force constants.

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
    phonon.force_constants = force_constants
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

    dos_raw = (phonon.projected_dos.frequency_points, phonon.projected_dos.projected_dos)
    p_doses = dict(zip(structure, dos_raw[1], strict=True))

    total_dos = PhononDos(dos_raw[0], dos_raw[1].sum(axis=0))
    return CompletePhononDos(structure, total_dos, p_doses)  # type:ignore[arg-type]


@requires(Phonopy, "phonopy is required to calculate phonon band structures")
def get_phonon_band_structure_from_fc(
    structure: Structure | IStructure,
    supercell_matrix: np.ndarray,
    force_constants: np.ndarray,
    mesh_density: float = 100.0,
    **kwargs,
) -> PhononBandStructure:
    """Get a uniform phonon band structure from phonopy force constants.

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
    phonon.force_constants = force_constants
    phonon.run_mesh(mesh_density, is_mesh_symmetry=False, is_gamma_center=True)
    mesh = phonon.get_mesh_dict()

    return PhononBandStructure(mesh["qpoints"], mesh["frequencies"], structure.lattice)


@requires(Phonopy, "phonopy is required to calculate phonon band structures")
def get_phonon_band_structure_symm_line_from_fc(
    structure: Structure | IStructure,
    supercell_matrix: np.ndarray,
    force_constants: np.ndarray,
    line_density: float = 20.0,
    symprec: float = 0.01,
    **kwargs,
) -> PhononBandStructureSymmLine:
    """Get a phonon band structure along a high symmetry path from phonopy force
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
    phonon.force_constants = force_constants

    k_path = HighSymmKpath(structure, symprec=symprec)

    kpoints, labels = k_path.get_kpoints(line_density=line_density, coords_are_cartesian=False)

    phonon.run_qpoints(kpoints)
    frequencies = phonon.qpoints.frequencies.T

    labels_dict = {a: k for a, k in zip(labels, kpoints, strict=True) if a != ""}

    return PhononBandStructureSymmLine(kpoints, frequencies, structure.lattice, labels_dict=labels_dict)


def get_gruneisenparameter(gruneisen_path, structure=None, structure_path=None) -> GruneisenParameter:
    """Get Gruneisen object from gruneisen.yaml file, as obtained from phonopy (Frequencies in THz!).
    The order is structure > structure path > structure from gruneisen dict.
    Newer versions of phonopy include the structure in the YAML file,
    the structure/structure_path is kept for compatibility.

    Args:
        gruneisen_path: Path to gruneisen.yaml file (frequencies have to be in THz!)
        structure: pymatgen Structure object
        structure_path: path to structure in a file (e.g., POSCAR)

    Returns:
        GruneisenParameter
    """
    gruneisen_dict = loadfn(gruneisen_path)

    if structure_path and structure is None:
        structure = Structure.from_file(structure_path)
    else:
        try:
            structure = get_structure_from_dict(gruneisen_dict)
        except ValueError as exc:
            raise ValueError("Please provide a structure or structure path") from exc

    q_pts, multiplicities, frequencies, gruneisen = ([] for _ in range(4))
    phonopy_labels_dict = {}

    for p in gruneisen_dict["phonon"]:
        q_pos = p["q-position"]
        q_pts.append(q_pos)
        m = p.get("multiplicity", 1)
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

    q_pts_np = np.array(q_pts)
    multiplicities_np = np.array(multiplicities)
    # transpose to match the convention in PhononBandStructure
    frequencies_np = np.transpose(frequencies)
    gruneisen_np = np.transpose(gruneisen)

    return GruneisenParameter(
        gruneisen=gruneisen_np,
        qpoints=q_pts_np,
        multiplicities=multiplicities_np,  # type:ignore[arg-type]
        frequencies=frequencies_np,
        structure=structure,
    )


def get_gs_ph_bs_symm_line_from_dict(
    gruneisen_dict, structure=None, structure_path=None, labels_dict=None, fit=False
) -> GruneisenPhononBandStructureSymmLine:
    r"""Create a pymatgen GruneisenPhononBandStructure object from the dictionary
    extracted by the gruneisen.yaml file produced by phonopy. The labels
    will be extracted from the dictionary, if present. If the 'eigenvector'
    key is found the eigendisplacements will be calculated according to the
    formula:

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
                q_pts_temp, frequencies_temp = [], []
                gruneisen_temp: list[list[float]] = []
                distance: list[float] = []
                for idx in range(pa["nqpoint"]):
                    bands = []
                    gruneisen_band: list[float] = []
                    for b in phonon[pa["nqpoint"] - idx - 1]["band"]:
                        bands.append(b["frequency"])
                        # Fraction of leftover points in current band
                        gruen = _extrapolate_grun(b, distance, gruneisen_temp, gruneisen_band, idx, pa)
                        gruneisen_band.append(gruen)
                    q_pos = phonon[pa["nqpoint"] - idx - 1]["q-position"]
                    q_pts_temp.append(q_pos)
                    d = phonon[pa["nqpoint"] - idx - 1]["distance"]
                    distance.append(d)
                    frequencies_temp.append(bands)
                    gruneisen_temp.append(gruneisen_band)
                    if "label" in phonon[pa["nqpoint"] - idx - 1]:
                        phonopy_labels_dict[phonon[pa["nqpoint"] - idx - 1]]["label"] = phonon[pa["nqpoint"] - idx - 1][
                            "q-position"
                        ]

                q_points.extend(list(reversed(q_pts_temp)))
                frequencies.extend(list(reversed(frequencies_temp)))
                gruneisen_params.extend(list(reversed(gruneisen_temp)))

            elif end["q-position"] == [0, 0, 0]:  # Gamma at end of band
                distance = []
                for idx in range(pa["nqpoint"]):
                    bands, gruneisen_band = [], []
                    for b in phonon[idx]["band"]:
                        bands.append(b["frequency"])
                        gruen = _extrapolate_grun(b, distance, gruneisen_params, gruneisen_band, idx, pa)
                        gruneisen_band.append(gruen)
                    q_pos = phonon[idx]["q-position"]
                    q_points.append(q_pos)
                    d = phonon[idx]["distance"]
                    distance.append(d)
                    frequencies.append(bands)
                    gruneisen_params.append(gruneisen_band)
                    if "label" in phonon[idx]:
                        phonopy_labels_dict[phonon[idx]["label"]] = phonon[idx]["q-position"]

            else:  # No Gamma in band
                distance = []
                for idx in range(pa["nqpoint"]):
                    bands, gruneisen_band = [], []
                    for b in phonon[idx]["band"]:
                        bands.append(b["frequency"])
                        gruneisen_band.append(b["gruneisen"])
                    q_pos = phonon[idx]["q-position"]
                    q_points.append(q_pos)
                    d = phonon[idx]["distance"]
                    distance.append(d)
                    frequencies.append(bands)
                    gruneisen_params.append(gruneisen_band)
                    if "label" in phonon[idx]:
                        phonopy_labels_dict[phonon[idx]["label"]] = phonon[idx]["q-position"]

    else:
        for pa in gruneisen_dict["path"]:
            for p in pa["phonon"]:
                q_pos = p["q-position"]
                q_points.append(q_pos)
                bands, gruneisen_bands = [], []
                for b in p["band"]:
                    bands.append(b["frequency"])
                    gruneisen_bands.append(b["gruneisen"])
                frequencies.append(bands)
                gruneisen_params.append(gruneisen_bands)
                if "label" in p:
                    phonopy_labels_dict[p["label"]] = p["q-position"]

    rec_lattice = structure.lattice.reciprocal_lattice
    labels_dict = labels_dict or phonopy_labels_dict
    return GruneisenPhononBandStructureSymmLine(
        qpoints=np.array(q_points),
        # transpose to match the convention in PhononBandStructure
        frequencies=np.transpose(frequencies),
        gruneisenparameters=np.transpose(gruneisen_params),
        lattice=rec_lattice,
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
    r"""Create a pymatgen GruneisenPhononBandStructure from a band.yaml file.
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
    thermal_displacements_yaml="thermal_displacement_matrices.yaml",
    structure_path="POSCAR",
):
    """Read "thermal_displacement_matrices.yaml" from phonopy and return a list of
    ThermalDisplacementMatrices objects.

    Args:
        thermal_displacements_yaml: path to thermal_displacement_matrices.yaml
        structure_path: path to POSCAR.

    Returns:
        list[ThermalDisplacementMatrices]
    """
    thermal_displacements_dict = loadfn(thermal_displacements_yaml)

    structure = Structure.from_file(structure_path)

    thermal_displacement_objects = []
    for matrix in thermal_displacements_dict["thermal_displacement_matrices"]:
        thermal_displacement_objects.append(
            ThermalDisplacementMatrices(
                thermal_displacement_matrix_cart=matrix["displacement_matrices"],
                temperature=matrix["temperature"],
                structure=structure,
                thermal_displacement_matrix_cif=matrix["displacement_matrices_cif"],
            )
        )

    return thermal_displacement_objects
