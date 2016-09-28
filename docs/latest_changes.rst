Change log
==========

v4.4.0
------
* Much more Pythonic API for modifying Structure/Molecule species. Now,
  strings, slices, and sequences should magically work, in addition to the
  previous API of simple int indices. Examples::

    s[0] = "Fe"
    s[0] = "Fe", [0.5, 0.5, 0.5]  # Replaces site and fractional coordinates.
    s[0] = "Fe", [0.5, 0.5, 0.5], {"spin": 2}  # Replaces site and fractional coordinates and properties.
    s[(0, 2, 3)] = "Fe"  # Replaces sites 0, 2 and 3 with Fe.
    s[0::2] = "Fe"  # Replaces all even index sites with Fe.
    s["Mn"] = "Fe"  # Replaces all Mn in the structure with Fe.
    s["Mn"] = "Fe0.5Co0.5"  # Replaces all Mn in the structure with Fe: 0.5, Co: 0.5, i.e.,creates a disordered structure!

* Massive update to internal representation of Bandstructure objects for
  memory and computational efficiency.
* Bug fixes to CIF parsing in some edge cases. (Will Richards).
