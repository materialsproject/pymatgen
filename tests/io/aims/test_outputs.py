from __future__ import annotations

import gzip
import json

from monty.json import MontyDecoder, MontyEncoder
from numpy.testing import assert_allclose

from pymatgen.core import Structure
from pymatgen.io.aims.outputs import AimsOutput
from pymatgen.util.testing import TEST_FILES_DIR

OUT_FILE_DIR = TEST_FILES_DIR / "io/aims/output_files"


def comp_images(test, ref):
    assert test.species == ref.species
    if isinstance(test, Structure):
        assert_allclose(test.lattice.matrix, ref.lattice.matrix, rtol=1e-7, atol=1e-7)
        test.translate_sites(range(len(test)), [0.0555, 0.0555, 0.0555])
        ref.translate_sites(range(len(ref)), [0.0555, 0.0555, 0.0555])

    assert_allclose(test.cart_coords, ref.cart_coords, rtol=1e-7, atol=1e-7)

    for key, val in ref.site_properties.items():
        assert_allclose(val, ref.site_properties[key])

    for key, val in ref.properties.items():
        if val is None and ref.properties[key] is None:
            continue
        assert_allclose(val, ref.properties[key])


def test_aims_output_si():
    si = AimsOutput.from_outfile(f"{OUT_FILE_DIR}/si.out")
    with gzip.open(f"{OUT_FILE_DIR}/si_ref.json.gz") as ref_file:
        si_ref = json.load(ref_file, cls=MontyDecoder)

    assert si_ref.metadata == si.metadata
    assert si_ref.structure_summary == si.structure_summary

    assert si_ref.n_images == si.n_images
    for ii in range(si.n_images):
        comp_images(si.get_results_for_image(ii), si_ref.get_results_for_image(ii))


def test_aims_output_h2o():
    h2o = AimsOutput.from_outfile(f"{OUT_FILE_DIR}/h2o.out")
    with gzip.open(f"{OUT_FILE_DIR}/h2o_ref.json.gz", mode="rt") as ref_file:
        h2o_ref = json.load(ref_file, cls=MontyDecoder)

    assert h2o_ref.metadata == h2o.metadata
    assert h2o_ref.structure_summary == h2o.structure_summary

    assert h2o_ref.n_images == h2o.n_images
    for ii in range(h2o.n_images):
        comp_images(h2o.get_results_for_image(ii), h2o_ref.get_results_for_image(ii))


def test_aims_output_si_dict():
    si = AimsOutput.from_outfile(f"{OUT_FILE_DIR}/si.out")
    si = json.loads(json.dumps(si.as_dict(), cls=MontyEncoder), cls=MontyDecoder)

    with gzip.open(f"{OUT_FILE_DIR}/si_ref.json.gz") as ref_file:
        si_ref = json.load(ref_file, cls=MontyDecoder)

    assert si_ref.metadata == si.metadata
    assert si_ref.structure_summary == si.structure_summary

    assert si_ref.n_images == si.n_images
    for ii in range(si.n_images):
        comp_images(si.get_results_for_image(ii), si_ref.get_results_for_image(ii))


def test_aims_output_h2o_dict():
    h2o = AimsOutput.from_outfile(f"{OUT_FILE_DIR}/h2o.out")
    h2o = json.loads(json.dumps(h2o.as_dict(), cls=MontyEncoder), cls=MontyDecoder)

    with gzip.open(f"{OUT_FILE_DIR}/h2o_ref.json.gz", mode="rt") as ref_file:
        h2o_ref = json.load(ref_file, cls=MontyDecoder)

    assert h2o_ref.metadata == h2o.metadata
    assert h2o_ref.structure_summary == h2o.structure_summary

    assert h2o_ref.n_images == h2o.n_images
    for ii in range(h2o.n_images):
        comp_images(h2o.get_results_for_image(ii), h2o_ref.get_results_for_image(ii))
