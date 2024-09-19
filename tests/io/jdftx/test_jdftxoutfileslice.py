from pathlib import Path

from pymatgen.core.units import Ha_to_eV

ex_files_dir = Path(__file__).parents[0] / "example_files"
ex_outslice_fname1 = ex_files_dir / "ex_out_slice_latmin"
ex_outslice1 = []
with open(ex_outslice_fname1) as f:
    for line in f:
        ex_outslice1.append(line)
ex_outslice1_known = {
    "Nbands": 42,
    "broadening_type": "Fermi",
    "broadening_value": 0.001,
    "truncation_type": "Periodic",
    "truncation_radius": None,
    "fluid": None,
    "prefix": None,
    "kgrid": [6, 2, 7],
    "latnIter": 100,
    "ionnIter": 0,
    "pptype": "GBRV",
    "is_gc": False,
    "geom_opt": True,
    "geom_opt_type": "lattice",
    "fftgrid": [28, 80, 28],
    "pwcut": 20 * Ha_to_eV,
    "rhocut": 100 * Ha_to_eV,
    "Emin": -1.780949 * Ha_to_eV,
    "HOMO": 0.704289 * Ha_to_eV,
    "EFermi": 0.704399 * Ha_to_eV,
    "LUMO": 0.704651 * Ha_to_eV,
    "Emax": 0.949497 * Ha_to_eV,
    "Egap": 0.000362 * Ha_to_eV,
    "traj_len": 7,
    "total_electrons": 64.000000,
}
