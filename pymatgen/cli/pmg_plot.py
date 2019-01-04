#!/usr/bin/env python
# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from collections import OrderedDict
from pymatgen import Structure
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp import Vasprun, Chgcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.diffraction.xrd import XRDCalculator


def get_dos_plot(args):
    v = Vasprun(args.dos_file)
    dos = v.complete_dos

    all_dos = OrderedDict()
    all_dos["Total"] = dos

    structure = v.final_structure

    if args.site:
        for i in range(len(structure)):
            site = structure[i]
            all_dos["Site " + str(i) + " " + site.specie.symbol] = \
                dos.get_site_dos(site)

    if args.element:
        syms = [tok.strip() for tok in args.element[0].split(",")]
        all_dos = {}
        for el, dos in dos.get_element_dos().items():
            if el.symbol in syms:
                all_dos[el] = dos
    if args.orbital:
        all_dos = dos.get_spd_dos()

    plotter = DosPlotter()
    plotter.add_dos_dict(all_dos)
    return plotter.get_plot()


def get_chgint_plot(args):
    chgcar = Chgcar.from_file(args.chgcar_file)
    s = chgcar.structure

    if args.inds:
        atom_ind = [int(i) for i in args.inds[0].split(",")]
    else:
        finder = SpacegroupAnalyzer(s, symprec=0.1)
        sites = [sites[0] for sites in
                 finder.get_symmetrized_structure().equivalent_sites]
        atom_ind = [s.sites.index(site) for site in sites]

    from pymatgen.util.plotting import pretty_plot
    plt = pretty_plot(12, 8)
    for i in atom_ind:
        d = chgcar.get_integrated_diff(i, args.radius, 30)
        plt.plot(d[:, 0], d[:, 1],
                 label="Atom {} - {}".format(i, s[i].species_string))
    plt.legend(loc="upper left")
    plt.xlabel("Radius (A)")
    plt.ylabel("Integrated charge (e)")
    plt.tight_layout()
    return plt


def get_xrd_plot(args):
    s = Structure.from_file(args.xrd_structure_file)
    c = XRDCalculator()
    return c.get_plot(s)


def plot(args):
    plt = None
    if args.chgcar_file:
        plt = get_chgint_plot(args)
    elif args.xrd_structure_file:
        plt = get_xrd_plot(args)
    elif args.dos_file:
        plt = get_dos_plot(args)

    if plt:
        if args.out_file:
            plt.savefig(args.out_file)
        else:
            plt.show()
