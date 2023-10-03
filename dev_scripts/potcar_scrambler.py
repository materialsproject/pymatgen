from __future__ import annotations

import numpy as np
from monty.serialization import zopen

from pymatgen.io.vasp import Potcar, PotcarSingle


class POTCAR_SCRAMBLER:

    """
    Takes a POTCAR and replaces its values with completely random values
    Does type matching and attempts prec matching on floats to ensure
    file is read correctly by Potcar and PotcarSingle classes

    Used to generate copyright-compliant POTCARs for PMG tests

    Recommended use:
        POTCAR_SCRAMBLER.from_file(
            input_filename = <input POTCAR name as str>,
            output_filename = <name of desired randomized POTCAR as str>
        )
    to generate a POTCAR with name `output_filename` with completely random values
    from existing POTCAR `input_filename`
    """

    def __init__(self, potcars: Potcar | PotcarSingle):
        if isinstance(potcars, PotcarSingle):
            potcars = [potcars]
        self.PSP_list = potcars
        # self.scrambled_potcars = []
        self.scrambled_potcars_str = ""
        for psp in self.PSP_list:
            scrambled_potcar_str = self.scramble_single_potcar(psp)
            self.scrambled_potcars_str += scrambled_potcar_str
        #    self.scrambled_potcars.append(PotcarSingle(scrambled_potcar_str))
        return

    def _rand_float_from_str_with_prec(self, istr: str, bloat: float = 1.5):
        nprec = len(istr.split(".")[1])
        bd = max(1, bloat * abs(float(istr)))
        return round(bd * np.random.rand(1)[0], nprec)

    def _read_fortran_str_and_scramble(self, istr: str, bloat: float = 1.5):
        istr = istr.strip()

        if istr.lower() in ["t", "f"] or istr.lower() in ["true", "false"]:
            return bool(np.random.randint(2))

        if (istr.upper() == istr.lower()) and istr[0].isnumeric():
            if "." in istr:
                return self._rand_float_from_str_with_prec(istr, bloat=bloat)
            else:
                intval = int(istr)
                fac = int(np.sign(intval))  # return int of same sign
                return fac * np.random.randint(abs(max(1, int(np.ceil(bloat * intval)))))
        else:
            try:
                float(istr)
                return self._rand_float_from_str_with_prec(istr, bloat=bloat)
            except ValueError:
                return istr

    def scramble_single_potcar(self, potcar: PotcarSingle):
        scrambled_potcar_str = ""
        for aline in potcar.data.split("\n")[:-1]:
            single_line_rows = aline.split(";")
            if "SHA256" in aline or "COPYR" in aline:
                # files not copyrighted, remove copyright statement
                # sha256 no longer applicable
                continue

            cline = ""
            for ibrow, brow in enumerate(single_line_rows):
                split_row = brow.split()
                for itmp, tmp in enumerate(split_row):
                    cline += f"{self._read_fortran_str_and_scramble(tmp)}"
                    if itmp < len(split_row) - 1:
                        cline += " "
                if len(single_line_rows) > 1 and ibrow == 0:
                    cline += "; "

            aux_str = ""
            if "TITEL" in aline:
                aux_str = " FAKE"
            scrambled_potcar_str += f"{cline}{aux_str}\n"
        return scrambled_potcar_str

    def to_file(self, filename: str):
        with zopen(filename, "wt") as f:
            f.write(self.scrambled_potcars_str)

    @staticmethod
    def from_file(input_filename: str, output_filename: str | None = None):
        psp = Potcar.from_file(input_filename)
        psp_scrambled = POTCAR_SCRAMBLER(psp)
        if output_filename:
            psp_scrambled.to_file(output_filename)
        return psp_scrambled


def generate_fake_potcar_libraries():
    """
    To test the `_gen_potcar_summary_stats` function in `pymatgen.io.vasp.inputs`,
    need a library of fake POTCARs which do not violate copyright
    """
    from os import path, system

    from pymatgen.io.vasp.sets import _load_yaml_config

    MPset = _load_yaml_config("MPRelaxSet")
    psp_variants = [MPset["POTCAR"][element] for element in MPset["POTCAR"]]

    output_base_dir = "./fake_POTCAR_library/"
    if path.isdir(output_base_dir):
        system(f"rm -rf {output_base_dir}")

    for func in PotcarSingle.functional_dir:
        func_dir = PotcarSingle.functional_dir[func]
        if not path.isdir(func_dir):
            continue

        for psp_name in psp_variants:
            rebase_dir = f"{output_base_dir}/{func_dir}/{psp_name}/"
            paths_to_try = [
                f"{func_dir}/POTCAR.{psp_name}",
                f"{func_dir}/POTCAR.{psp_name}.gz",
                f"{func_dir}/{psp_name}/POTCAR",
                f"{func_dir}/{psp_name}/POTCAR.gz",
            ]
            for potcar_path in paths_to_try:
                if path.isfile(potcar_path):
                    system(f"mkdir -p {rebase_dir}")
                    POTCAR_SCRAMBLER.from_file(input_filename=potcar_path, output_filename=f"{rebase_dir}/POTCAR.gz")
                    break


if __name__ == "__main__":
    generate_fake_potcar_libraries()
