import os

for parent, subdir, files in os.walk("../pymatgen"):
    for fn in files:
        if fn.endswith(".py") and fn != "__init__.py":
            fp = os.path.join(parent, fn)
            with open(fp) as f:
                contents = f.read()
            if "unicode_literals" not in contents:
                contents = "from __future__ import unicode_literals\n" + contents
            lines = contents.split("\n")
            future = []
            clean = []
            while len(lines) > 0:
                l = lines.pop(0)
                if l.strip().startswith("from __future__"):
                    future.append(l.strip().split("import")[-1].strip())
                elif not l.strip().startswith("# coding"):
                    clean.append(l)

            clean = (
                "# coding: utf-8\n\nfrom __future__ import "
                + ", ".join(future)
                + "\n\n"
                + "\n".join(clean).strip()
                + "\n"
            )
            with open(fp, "w") as f:
                f.write(clean)
