from __future__ import annotations

import json

from pymatgen.core.periodic_table import Element

element_properties = [
    "X",
    "atomic_mass",
    "atomic_mass_number",
    "atomic_orbitals_eV",
    "atomic_radius",
    "average_anionic_radius",
    "average_cationic_radius",
    "average_ionic_radius",
    "block",
    "common_oxidation_states",
    "electron_affinity",
    "electronic_structure",
    "full_electronic_structure",
    "ground_state_term_symbol",
    "group",
    "icsd_oxidation_states",
    "ionic_radii",
    "ionization_energy",
    "is_actinoid",
    "is_alkali",
    "is_alkaline",
    "is_chalcogen",
    "is_halogen",
    "is_lanthanoid",
    "is_metal",
    "is_metalloid",
    "is_noble_gas",
    "is_post_transition_metal",
    "is_quadrupolar",
    "is_radioactive",
    "is_rare_earth",
    "is_transition_metal",
    "iupac_ordering",
    "max_oxidation_state",
    "min_oxidation_state",
    "n_electrons",
    "nmr_quadrupole_moment",
    "number",
    "oxidation_states",
    "row",
    "term_symbols",
    "valence",
]

extra_properties = [
    "mendeleev_no",
    "electrical_resistivity",
    "velocity_of_sound",
    "reflectivity",
    "refractive_index",
    "poissons_ratio",
    "molar_volume",
    "thermal_conductivity",
    "boiling_point",
    "melting_point",
    "critical_temperature",
    "superconduction_temperature",
    "liquid_range",
    "bulk_modulus",
    "youngs_modulus",
    "brinell_hardness",
    "rigidity_modulus",
    "mineral_hardness",
    "vickers_hardness",
    "density_of_solid",
    "atomic_radius_calculated",
    "van_der_waals_radius",
    "atomic_orbitals",
    "coefficient_of_linear_thermal_expansion",
    "ground_state_term_symbol",
    "valence",
    "ground_level",
    "ionization_energies",
    "metallic_radius",
]

all_properties = sorted(set(element_properties + extra_properties))


def dump_properties(path: str):
    data: dict = {prop: {} for prop in all_properties}
    for el in Element:
        for prop in all_properties:
            try:
                val = getattr(el, prop)
                if isinstance(val, set):
                    val = list(val)
                data[prop][el.symbol] = val
            except Exception:
                data[prop][el.symbol] = None

    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)


if __name__ == "__main__":
    dump_properties("element_properties.json")
