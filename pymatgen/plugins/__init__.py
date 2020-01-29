"""
Plugin architecture for pymatgen.
"""

import pkg_resources

discovered_plugins = {
    entry_point.name: entry_point.load()
    for entry_point
    in pkg_resources.iter_entry_points('pymatgen.plugins')
}

locals().update(discovered_plugins)
