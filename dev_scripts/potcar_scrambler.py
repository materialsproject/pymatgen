from __future__ import annotations

import os
import shutil
import warnings

import numpy as np
from monty.os.path import zpath
from monty.serialization import zopen

from pymatgen.core import SETTINGS
from pymatgen.io.vasp import Potcar, PotcarSingle
from pymatgen.io.vasp.sets import _load_yaml_config


class PotcarScrambler:

    """
    Takes a POTCAR and replaces its values with completely random values
    Does type matching and attempts precision matching on floats to ensure
    file is read correctly by Potcar and PotcarSingle classes.

    Used to generate copyright-compliant POTCARs for PMG tests.

    In case of questions, contact Aaron Kaplan <adkaplan@lbl.gov>.

    Recommended use:
        PotcarScrambler.from_file(
            input_filename = <input POTCAR name as str>,
            output_filename = <name of desired randomized POTCAR as str>
        )
    to generate a POTCAR with name `output_filename` with completely random values
    from existing POTCAR `input_filename`
    """

    def __init__(self, potcars: Potcar | PotcarSingle):
        if isinstance(potcars, PotcarSingle):
            self.PSP_list = [potcars]
        else:
            self.PSP_list = potcars
        self.scrambled_potcars_str = ""
        for psp in self.PSP_list:
            scrambled_potcar_str = self.scramble_single_potcar(psp)
            self.scrambled_potcars_str += scrambled_potcar_str
        return

    def _rand_float_from_str_with_prec(self, input_str: str, bloat: float = 1.5):
        n_prec = len(input_str.split(".")[1])
        bd = max(1, bloat * abs(float(input_str)))
        return round(bd * np.random.rand(1)[0], n_prec)

    def _read_fortran_str_and_scramble(self, input_str: str, bloat: float = 1.5):
        input_str = input_str.strip()

        if input_str.lower() in ("t", "f", "true", "false"):
            return bool(np.random.randint(2))

        if input_str.upper() == input_str.lower() and input_str[0].isnumeric():
            if "." in input_str:
                return self._rand_float_from_str_with_prec(input_str, bloat=bloat)
            integer = int(input_str)
            fac = int(np.sign(integer))  # return int of same sign
            return fac * np.random.randint(abs(max(1, int(np.ceil(bloat * integer)))))
        try:
            float(input_str)
            return self._rand_float_from_str_with_prec(input_str, bloat=bloat)
        except ValueError:
            return input_str

    def scramble_single_potcar(self, potcar: PotcarSingle):
        scrambled_potcar_str = ""
        for line in potcar.data.split("\n")[:-1]:
            single_line_rows = line.split(";")
            if "SHA256" in line or "COPYR" in line:
                # files not copyrighted, remove copyright statement
                # sha256 no longer applicable
                continue

            cline = ""
            for idx, row in enumerate(single_line_rows):
                split_row = row.split()
                for itmp, tmp in enumerate(split_row):
                    cline += f"{self._read_fortran_str_and_scramble(tmp)}"
                    if itmp < len(split_row) - 1:
                        cline += " "
                if len(single_line_rows) > 1 and idx == 0:
                    cline += "; "

            aux_str = ""
            if "TITEL" in line:
                aux_str = " FAKE"
            scrambled_potcar_str += f"{cline}{aux_str}\n"
        return scrambled_potcar_str

    def to_file(self, filename: str):
        with zopen(filename, "wt") as f:
            f.write(self.scrambled_potcars_str)

    @classmethod
    def from_file(cls, input_filename: str, output_filename: str | None = None):
        psp = Potcar.from_file(input_filename)
        psp_scrambled = cls(psp)
        if output_filename:
            psp_scrambled.to_file(output_filename)
        return psp_scrambled


def generate_fake_potcar_libraries():
    """
    To test the `_gen_potcar_summary_stats` function in `pymatgen.io.vasp.inputs`,
    need a library of fake POTCARs which do not violate copyright
    """
    mp_relax_set = _load_yaml_config("MPRelaxSet")
    psp_variants = [mp_relax_set["POTCAR"][element] for element in mp_relax_set["POTCAR"]]

    output_dir = "./fake_potcar_library/"
    shutil.rmtree(output_dir, ignore_errors=True)

    vasp_psp_dir = SETTINGS.get("PMG_VASP_PSP_DIR")
    src_dirs = [f"{vasp_psp_dir}/{func_dir}" for func_dir in PotcarSingle.functional_dir.values()]

    if not any(map(os.path.isdir, src_dirs)):
        raise RuntimeError(f"No input POTCAR library found, tried {src_dirs}")

    for func_dir in src_dirs:
        if not os.path.isdir(func_dir):
            continue

        for psp_name in psp_variants:
            rebase_dir = f"{output_dir}/{func_dir}/{psp_name}/"
            paths_to_try = [
                zpath(f"{func_dir}/POTCAR.{psp_name}"),
                zpath(f"{func_dir}/{psp_name}/POTCAR"),
            ]
            if not any(map(os.path.isfile, paths_to_try)):
                warnings.warn(f"Could not find {psp_name} in {paths_to_try}")
            for potcar_path in paths_to_try:
                if os.path.isfile(potcar_path):
                    os.makedirs(rebase_dir, exist_ok=True)
                    PotcarScrambler.from_file(input_filename=potcar_path, output_filename=f"{rebase_dir}/POTCAR.gz")
                    break


if __name__ == "__main__":
    generate_fake_potcar_libraries()
