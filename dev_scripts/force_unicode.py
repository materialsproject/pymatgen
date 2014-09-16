# coding: utf-8
import os
for parent, subdir, files in os.walk("../pymatgen"):
    for fn in files:
        if fn.endswith(".py"):
            fp = os.path.join(parent, fn)
            with open(fp) as f:
                contents = f.read()
            if "unicode_literals" not in contents:
                print fp
                contents = "from __future__ import unicode_literals\n" + contents
                with open(fp, "w") as f:
                    f.write(contents)

