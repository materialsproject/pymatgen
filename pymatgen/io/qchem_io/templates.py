"""Templates for each section of a qchem input file. These functions return a qchem formatted string for each section"""


def molecule_template(molecule):
    mol_list = []
    mol_list.append("$molecule")
    mol_list.append(" {charge} {spin_mult}".format(charge=molecule.charge, spin_mult=molecule.spin_multiplicity))
    for site in molecule.sites:
        mol_list.append(" {atom}     {x: .10f}     {y: .10f}     {z: .10f}".format(atom=site.species_string,
                                                                                   x=site.x, y=site.y, z=site.z))
    mol_list.append("$end")
    return '\n'.join(mol_list)


def rem_template(rem):
    rem_list = []
    rem_list.append("$rem")
    for key, value in rem.items():
        rem_list.append("   {key} = {value}".format(key=key, value=value))
    rem_list.append("$end")
    return '\n'.join(rem_list)


def opt_template(opt):
    opt_list = []
    opt_list.append("$opt")
    # loops over all opt sections
    for key, value in opt.items():
        opt_list.append("   {section}".format(section=key))
        # loops over all values within the section
        for i in value:
            opt_list.append("   {val}".format(val=i))
        opt_list.append("   END{section}".format(section=key))
        opt_list.append("")
    del opt_list[-1]
    opt_list.append("$end")
    return '\n'.join(opt_list)