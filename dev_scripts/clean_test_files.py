#!/usr/bin/env python


import os

to_delete = []
for parent, subdir, files in os.walk("../test_files"):
    for fn in files:
        found = None
        for p, s, t in os.walk("../pymatgen"):
            for fn2 in t:
                if p.endswith("tests"):
                    with open(os.path.join(p, fn2)) as f:
                        contents = f.read()
                    if fn in contents:
                        found = os.path.join(p, fn2)
                        break
        if found is None:
            to_delete.append(os.path.join(parent, fn))
            os.remove(os.path.join(parent, fn))

print to_delete