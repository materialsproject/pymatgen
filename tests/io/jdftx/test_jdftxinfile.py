from pathlib import Path
from atomate2.jdftx.io.jdftxinfile import JDFTXInfile
from pytest import approx
import pytest
from pymatgen.util.typing import PathLike
import os

ex_files_dir = Path(__file__).parents[0] / "example_files"

ex_infile1_fname = ex_files_dir / "CO.in"
ex_infile1_knowns = {
    "dump-name": "$VAR",
    "initial-state": "$VAR",
    "elec-ex-corr": "gga",
    "van-der-waals": "D3",
    "elec-cutoff": {"Ecut": 20.0,"EcutRho":100.0},
    "elec-n-bands": 15,
    # "kpoint-folding": [1,1,1],
    "kpoint-folding": {"n0": 1, "n1": 1, "n2": 1},
    "spintype": "z-spin",
    "core-overlap-check": "none",
    "converge-empty-states": True,
    "latt-move-scale": {"s0": 0., "s1": 0., "s2": 0.},
    # "latt-move-scale": [0,0,0],
    "symmetries": "none",
    "fluid": {"type": "LinearPCM"},
    "pcm-variant": "CANDLE",
    "fluid-solvent": [{"name": "H2O"}],
    "fluid-cation": {"name": "Na+", "concentration": 0.5},
    "fluid-anion": {"name": "F-", "concentration": 0.5},
}
# jif = JDFTXInfile.from_file(ex_infile1_fname)
# out = jif.get_list_representation(jif)
# jif.get_text_list()

@pytest.mark.parametrize("infile_fname,knowns", [(ex_infile1_fname, ex_infile1_knowns)])
def test_JDFTXInfile_knowns_simple(infile_fname: PathLike, knowns: dict):
    jif = JDFTXInfile.from_file(infile_fname)
    for key in knowns:
        assert is_identical_jif_val(jif[key], knowns[key])


@pytest.mark.parametrize("infile_fname", [ex_infile1_fname])
def test_JDFTXInfile_self_consistency(infile_fname: PathLike):
    jif = JDFTXInfile.from_file(infile_fname)
    dict_jif = jif.as_dict()
    # # Commenting out tests with jif2 due to the list representation asserted
    jif2 = JDFTXInfile.get_dict_representation(JDFTXInfile.from_dict(dict_jif))
    jif3 = JDFTXInfile.from_str(str(jif))
    tmp_fname = ex_files_dir / "tmp.in"
    jif.write_file(tmp_fname)
    jif4 = JDFTXInfile.from_file(tmp_fname)
    jifs = [jif, jif2, jif3, jif4]
    for i in range(len(jifs)):
        for j in range(i+1, len(jifs)):
            print(f"{i}, {j}")
            assert is_identical_jif(jifs[i], jifs[j])
    os.remove(tmp_fname)


def is_identical_jif(jif1: JDFTXInfile | dict, jif2: JDFTXInfile | dict):
    for key in jif1:
        if key not in jif2:
            return False
        else:
            v1 = jif1[key]
            v2 = jif2[key]
            assert is_identical_jif_val(v1, v2)
    return True


def is_identical_jif_val(v1, v2):
    if type(v1) != type(v2):
        return False
    elif isinstance(v1, float):
        return v1 == approx(v2)
    elif True in [isinstance(v1, str), isinstance(v1, int)]:
        return v1 == v2
    elif True in [isinstance(v1, list)]:
        if len(v1) != len(v2):
            return False
        for i, v in enumerate(v1):
            if not is_identical_jif_val(v, v2[i]):
                return False
        return True
    elif True in [isinstance(v1, dict)]:
        return is_identical_jif(v1, v2)



# test_JDFTXInfile_self_consistency(ex_infile1_fname)
