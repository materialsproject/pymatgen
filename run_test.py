from __future__ import division

"""
Script to wrap around nosetests to only run on files modified in the past 10
commits as well as 10% of all files.
"""

import os
import sys
import subprocess
import random


run_ratio = 1/10

output = subprocess.check_output(["git", "diff", "--name-only", "HEAD~10"])
files_changed = [f for f in output.split("\n") if f.startswith("pymatgen")]

must_run = []
for f in files_changed:
    d = os.path.dirname(f)
    b = os.path.basename(f)
    testname = os.path.join(d, "tests", "test_" + b)
    if os.path.exists(testname):
        must_run.append(testname)
        
can_run = []
for parent, subdir, files in os.walk("pymatgen"):
    for f in files:
        if parent.endswith("tests") and f.startswith("test_") and f.endswith(".py") and f not in must_run:
            can_run.append(os.path.join(parent, f))

nrun = int(run_ratio * len(can_run))

if random.randint(1, 50) % 50 == 0:
    #One in fifty times, we will run a full test.
    to_run = must_run + can_run
else:
    to_run = random.sample(can_run, nrun) + must_run

status = subprocess.call(["nosetests", "-v"] + to_run)
sys.exit(status)
