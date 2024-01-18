from __future__ import annotations

import gzip


def compare_files(ref_file, test_file):
    with open(test_file) as tf:
        test_lines = tf.readlines()[5:]

    with gzip.open(f"{ref_file}.gz", mode="rt") as rf:
        ref_lines = rf.readlines()[5:]

    for test_line, ref_line in zip(test_lines, ref_lines):
        if "species_dir" in ref_line:
            continue
        assert test_line.strip() == ref_line.strip()
