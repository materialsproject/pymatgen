

def bulk_coordination(slab, bondlength, bond):

    # -If bond is list of 2 species, we find coordination for
    #   the first specie relative to the second only [s1, s2]
    # -If list of one species, find coordination for all neighbors [s]
    # -If first item is a list, we are looking for
    #   the coordination of a polyhedron [[s1a, s1b], s2]
    # IMPORTANT NOTE, cannot specify the specific bondlength of your polyhedron
    # if you are looking for the coordination of a polyhedron. The bondlength
    # will be the same as that of the polyhedron vertex and the next species

    center_ion = bond[0] if type(bond[0]).__name__ == 'str' else bond[0][0]
    mean_cn = []
    ucell = slab.oriented_unit_cell
    for i, site in enumerate(ucell):
        cn = 0
        if site.species_string == center_ion:
            nn = ucell.get_neighbors(site, bondlength,
                                     include_index=True)

            for n in nn:
                # If we're dealing with the coordination of a single atoms
                if type(bond[0]).__name__ == 'str':
                    if len(bond) == 2:
                        if n[0].species_string == bond[1]:
                            cn += 1
                    else:
                        cn += 1
                # If we're dealing with the coordination of a polyhedron
                else:
                    # Check if the vertex of the polyhedron is the correct species
                    if n[0].species_string == bond[0][1]:
                        # Get coordinated sites with that vertex
                        vert_n = ucell.get_neighbors(ucell[n[2]], bondlength,
                                                     include_index=True)
                        for nnn in vert_n:
                            # Check the coordinated site of vertex is
                            # not the center of the polyhedron
                            if nnn[2] == i:
                                continue
                            if len(bond) == 2:
                                if nnn[0].species_string == bond[1]:
                                    cn += 1
                            else:
                                cn += 1

            mean_cn.append(cn)

    return min(mean_cn)


def surface_coordination(slab, bonds, top=True):

    """
    A function that analyzes the coordination environment of bulk atoms, surface atoms
        and broken bonds. Returns a dictionary describing each type of environment for
        each type of bond: eg. {"bulk": {(species1, species2): 12}, "surface": {(species1,
        species2): 6}, "broken": {(species1, species2): 3}}

    Args:
        slab (Slab): Initial input slab.
        bonds ({(specie1, specie2): max_bond_dist}: bonds are
            specified as a dict of tuples: float of specie1, specie2
            and the max bonding distance. For example, PO4 groups may be
            defined as {("P", "O"): 3}.
    """

    bulk_bonds, surface_bonds, broken_bonds = {}, {}, {}
    for bond in bonds.keys():

        # First we check the cn of the bulk for each type of bond
        cn = bulk_coordination(slab, bonds[bond], bond)

        # Next we use the cn of the bulk as
        # reference to find the number of broken bonds
        center_ion = bond[0] if type(bond[0]).__name__ == 'str' else bond[0][0]
        cnb, tot_surf_cn, bb_cn = 0, 0, 0
        for i, site in enumerate(slab):
            if str(site.specie) == center_ion:
                nn = slab.get_neighbors(site, bonds[bond],
                                        include_index=True)

                def count_cnb(cnb, tot_surf_cn, bb_cn):
                    slab_cn = 0
                    for n in nn:

                        # If we're dealing with the coordination of a single atoms
                        if type(bond[0]).__name__ == 'str':
                            if len(bond) == 2:
                                if n[0].species_string == bond[1]:
                                    slab_cn += 1
                            else:
                                slab_cn += 1

                        # If we're dealing with the coordination of a polyhedron
                        else:
                            # Check if the vertex of the polyhedron is the correct species
                            if n[0].species_string == bond[0][1]:
                                # Get coordinated sites with that vertex
                                vert_n = slab.get_neighbors(slab[n[2]], bonds[bond],
                                                            include_index=True)
                                for nnn in vert_n:
                                    # Check the coordinated site of vertex is
                                    # not the center of the polyhedron
                                    if nnn[2] == i:
                                        continue
                                    if len(bond) == 2:
                                        if nnn[0].species_string == bond[1]:
                                            slab_cn += 1
                                    else:
                                        slab_cn += 1

                    cnb += cn
                    tot_surf_cn += slab_cn
                    bb_cn += cn - slab_cn
                    return cnb, tot_surf_cn, bb_cn

                if top and site.frac_coords[2] > slab.center_of_mass:
                    cnb, tot_surf_cn, bb_cn = count_cnb(cnb, tot_surf_cn, bb_cn)
                if not top and site.frac_coords[2] < slab.center_of_mass:
                    cnb, tot_surf_cn, bb_cn = count_cnb(cnb, tot_surf_cn, bb_cn)

        bulk_bonds[bond] = cnb / slab.surface_area
        surface_bonds[bond] = tot_surf_cn / slab.surface_area
        broken_bonds[bond] = bb_cn / slab.surface_area

    return {"bulk": bulk_bonds, "surface": surface_bonds, "broken": broken_bonds}