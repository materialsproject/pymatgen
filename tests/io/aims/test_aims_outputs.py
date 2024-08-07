from __future__ import annotations

import gzip
import json
from pathlib import Path

from monty.json import MontyDecoder, MontyEncoder
from numpy.testing import assert_allclose

from pymatgen.core import Structure
from pymatgen.io.aims.outputs import AimsOutput

outfile_dir = Path(__file__).parent / "output_files"


def comp_images(test, ref):
    assert test.species == ref.species
    if isinstance(test, Structure):
        assert_allclose(test.lattice.matrix, ref.lattice.matrix)
    assert_allclose(test.cart_coords, ref.cart_coords, atol=1e-12)

    for key, val in test.site_properties.items():
        assert_allclose(val, ref.site_properties[key])

    for key, val in test.properties.items():
        assert_allclose(val, ref.properties[key])


def test_aims_output_si():
    si = AimsOutput.from_outfile(f"{outfile_dir}/si.out.gz")
    with gzip.open(f"{outfile_dir}/si_ref.json.gz") as ref_file:
        si_ref = json.load(ref_file, cls=MontyDecoder)

    assert si_ref.metadata == si.metadata
    assert si_ref.structure_summary == si.structure_summary

    assert si_ref.n_images == si.n_images
    for ii in range(si.n_images):
        comp_images(si.get_results_for_image(ii), si_ref.get_results_for_image(ii))


def test_aims_output_h2o():
    h2o = AimsOutput.from_outfile(f"{outfile_dir}/h2o.out.gz")
    with gzip.open(f"{outfile_dir}/h2o_ref.json.gz", mode="rt") as ref_file:
        h2o_ref = json.load(ref_file, cls=MontyDecoder)

    assert h2o_ref.metadata == h2o.metadata
    assert h2o_ref.structure_summary == h2o.structure_summary

    assert h2o_ref.n_images == h2o.n_images
    for ii in range(h2o.n_images):
        comp_images(h2o.get_results_for_image(ii), h2o_ref.get_results_for_image(ii))


def test_aims_output_si_str():
    with gzip.open(f"{outfile_dir}/si.out.gz", mode="rt") as si_out:
        si = AimsOutput.from_str(si_out.read())

    with gzip.open(f"{outfile_dir}/si_ref.json.gz", mode="rt") as ref_file:
        si_ref = json.load(ref_file, cls=MontyDecoder)

    assert si_ref.metadata == si.metadata
    assert si_ref.structure_summary == si.structure_summary

    assert si_ref.n_images == si.n_images
    for ii in range(si.n_images):
        comp_images(si.get_results_for_image(ii), si_ref.get_results_for_image(ii))


def test_aims_output_h2o_str():
    with gzip.open(f"{outfile_dir}/h2o.out.gz", mode="rt") as h2o_out:
        h2o = AimsOutput.from_str(h2o_out.read())

    with gzip.open(f"{outfile_dir}/h2o_ref.json.gz", mode="rt") as ref_file:
        h2o_ref = json.load(ref_file, cls=MontyDecoder)

    assert h2o_ref.metadata == h2o.metadata
    assert h2o_ref.structure_summary == h2o.structure_summary

    assert h2o_ref.n_images == h2o.n_images
    for ii in range(h2o.n_images):
        comp_images(h2o.get_results_for_image(ii), h2o_ref.get_results_for_image(ii))


def test_aims_output_si_dict():
    si = AimsOutput.from_outfile(f"{outfile_dir}/si.out.gz")
    si = json.loads(json.dumps(si.as_dict(), cls=MontyEncoder), cls=MontyDecoder)

    with gzip.open(f"{outfile_dir}/si_ref.json.gz") as ref_file:
        si_ref = json.load(ref_file, cls=MontyDecoder)

    assert si_ref.metadata == si.metadata
    assert si_ref.structure_summary == si.structure_summary

    assert si_ref.n_images == si.n_images
    for ii in range(si.n_images):
        comp_images(si.get_results_for_image(ii), si_ref.get_results_for_image(ii))


def test_aims_output_h2o_dict():
    h2o = AimsOutput.from_outfile(f"{outfile_dir}/h2o.out.gz")
    h2o = json.loads(json.dumps(h2o.as_dict(), cls=MontyEncoder), cls=MontyDecoder)

    with gzip.open(f"{outfile_dir}/h2o_ref.json.gz", mode="rt") as ref_file:
        h2o_ref = json.load(ref_file, cls=MontyDecoder)

    assert h2o_ref.metadata == h2o.metadata
    assert h2o_ref.structure_summary == h2o.structure_summary

    assert h2o_ref.n_images == h2o.n_images
    for ii in range(h2o.n_images):
        comp_images(h2o.get_results_for_image(ii), h2o_ref.get_results_for_image(ii))
