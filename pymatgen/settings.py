"""
Global settings for pymatgen.
"""

import os

try:
    import ruamel.yaml as yaml
except ImportError:
    try:
        import ruamel_yaml as yaml  # type: ignore  # noqa
    except ImportError:
        import yaml  # type: ignore # noqa

SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings():
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
    except IOError:
        # If there are any errors, default to using environment variables
        # if present.
        d = {}
        for k, v in os.environ.items():
            if k.startswith("PMG_"):
                d[k] = v
            elif k in ["VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"]:
                d["PMG_" + k] = v
    d = d or {}
    return dict(d)


SETTINGS = _load_pmg_settings()
locals().update(SETTINGS)
