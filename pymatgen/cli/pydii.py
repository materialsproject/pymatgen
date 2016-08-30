#!/usr/bin/env python

from __future__ import division, unicode_literals

"""
A script with tools for computing point defect concentrations.
Manual and citation for the script, DOI: 10.1016/j.cpc.2015.03.015
"""

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "3.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__date__ = "March 16, 2015"

import argparse
import os
import glob

from monty.serialization import loadfn, dumpfn
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.analysis.defects.point_defects import Vacancy
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp import Kpoints
from pymatgen.io.vasp import Vasprun
from pymatgen.matproj.rest import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.dilute_solution_model import \
            compute_defect_density, solute_defect_density


def get_sc_scale(inp_struct, final_site_no):
    lengths = inp_struct.lattice.abc
    no_sites = inp_struct.num_sites
    mult = (final_site_no/no_sites*lengths[0]*lengths[1]*lengths[2]) ** (1/3)
    num_mult = [int(round(mult/l)) for l in lengths]
    num_mult = [i if i > 0 else 1 for i in num_mult]
    return num_mult


def vac_antisite_def_struct_gen(args):
    mpid = args.mpid
    mapi_key = args.mapi_key
    cellmax = args.cellmax

    if not mpid:
        print ("============\nERROR: Provide an mpid\n============")
        return

    if not mapi_key:
        with MPRester() as mp:
            struct = mp.get_structure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            struct = mp.get_structure_by_material_id(mpid)

    prim_struct_sites = len(struct.sites)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    conv_struct_sites = len(struct.sites)
    conv_prim_rat = int(conv_struct_sites/prim_struct_sites)
    sc_scale = get_sc_scale(struct,cellmax)

    mpvis = MPRelaxSet(struct, user_incar_settings={"LDAU": False})

    # Begin defaults: All default settings.
    blk_vasp_incar_param = {'IBRION':-1,'EDIFF':1e-4,'EDIFFG':0.001,'NSW':0,}
    def_vasp_incar_param = {'ISIF':2,'NELM':99,'IBRION':2,'EDIFF':1e-6, 
                            'EDIFFG':0.001,'NSW':40,}
    kpoint_den = 6000
    # End defaults
    
    ptcr_flag = True
    try:
        potcar = mpvis.potcar
    except:
        print ("VASP POTCAR folder not detected.\n" \
              "Only INCAR, POSCAR, KPOINTS are generated.\n" \
              "If you have VASP installed on this system, \n" \
              "refer to pymatgen documentation for configuring the settings.")
        ptcr_flag = False


    vac = Vacancy(struct, {}, {})
    scs = vac.make_supercells_with_defects(sc_scale)
    site_no = scs[0].num_sites
    if site_no > cellmax:
        max_sc_dim = max(sc_scale)
        i = sc_scale.index(max_sc_dim)
        sc_scale[i] -= 1
        scs = vac.make_supercells_with_defects(sc_scale)

    for i in range(len(scs)):
        sc = scs[i]
        mpvis = MPRelaxSet(sc, user_incar_settings={"LDAU": False})
        poscar = mpvis.poscar
        kpoints = Kpoints.automatic_density(sc,kpoint_den)
        incar = mpvis.incar
        if ptcr_flag:
            potcar = mpvis.potcar

        interdir = mpid
        if not i:
            fin_dir = os.path.join(interdir,'bulk')
            try:
                os.makedirs(fin_dir)
            except:
                pass
            incar.update(blk_vasp_incar_param)
            incar.write_file(os.path.join(fin_dir,'INCAR'))
            poscar.write_file(os.path.join(fin_dir,'POSCAR'))
            if ptcr_flag:
                potcar.write_file(os.path.join(fin_dir,'POTCAR'))
            kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))
        else:
            blk_str_sites = set(scs[0].sites)
            vac_str_sites = set(sc.sites)
            vac_sites = blk_str_sites - vac_str_sites
            vac_site = list(vac_sites)[0]
            site_mult = int(vac.get_defectsite_multiplicity(i-1)/conv_prim_rat)
            vac_site_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol

            vac_dir ='vacancy_{}_mult-{}_sitespecie-{}'.format(str(i),
                    site_mult, vac_symbol)
            fin_dir = os.path.join(interdir,vac_dir)
            try:
                os.makedirs(fin_dir)
            except:
                pass
            incar.update(def_vasp_incar_param)
            poscar.write_file(os.path.join(fin_dir,'POSCAR'))
            incar.write_file(os.path.join(fin_dir,'INCAR'))
            if ptcr_flag:
                potcar.write_file(os.path.join(fin_dir,'POTCAR'))
            kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))

            # Antisite generation at all vacancy sites
            struct_species = scs[0].types_of_specie
            for specie in set(struct_species)-set([vac_site_specie]):
                subspecie_symbol = specie.symbol
                anti_struct = sc.copy()
                anti_struct.append(specie, vac_site.frac_coords)
                mpvis = MPRelaxSet(anti_struct, user_incar_settings={"LDAU": False})
                poscar = mpvis.poscar
                incar = mpvis.incar
                incar.update(def_vasp_incar_param)
                as_dir ='antisite_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                        str(i), site_mult, vac_symbol, subspecie_symbol)
                fin_dir = os.path.join(interdir,as_dir)
                try:
                    os.makedirs(fin_dir)
                except:
                    pass
                poscar.write_file(os.path.join(fin_dir,'POSCAR'))
                incar.write_file(os.path.join(fin_dir,'INCAR'))
                if ptcr_flag:
                    potcar.write_file(os.path.join(fin_dir,'POTCAR'))
                kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))


def substitute_def_struct_gen(args):
    mpid = args.mpid 
    solute = args.solute 
    mapi_key = args.mapi_key 
    cellmax = args.cellmax 

    if not mpid:
        print ("============\nERROR: Provide an mpid\n============")
        return
    if not solute:
        print ("============\nERROR: Provide solute atom\n============")
        return

    if not mapi_key:
        with MPRester() as mp:
            struct = mp.get_structure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            struct = mp.get_structure_by_material_id(mpid)
    prim_struct_sites = len(struct.sites)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    conv_struct_sites = len(struct.sites)
    conv_prim_rat = int(conv_struct_sites/prim_struct_sites)

    mpvis = MPRelaxSet(struct, user_incar_settings={"LDAU": False})

    # Begin defaults: All default settings.
    blk_vasp_incar_param = {'IBRION':-1,'EDIFF':1e-4,'EDIFFG':0.001,'NSW':0,}
    def_vasp_incar_param = {'ISIF':2,'NELM':99,'IBRION':2,'EDIFF':1e-6, 
                            'EDIFFG':0.001,'NSW':40,}
    kpoint_den = 6000
    # End defaults
    
    # Check if POTCAR file can be geneated
    ptcr_flag = True
    try:
        potcar = mpvis.potcar
    except:
        print ("VASP POTCAR folder not detected.\n" \
              "Only INCAR, POSCAR, KPOINTS are generated.\n" \
              "If you have VASP installed on this system, \n" \
              "refer to pymatgen documentation for configuring the settings.")
        ptcr_flag = False

    vac = Vacancy(struct, {}, {})
    sc_scale = get_sc_scale(struct,cellmax)
    scs = vac.make_supercells_with_defects(sc_scale)
    site_no = scs[0].num_sites
    if site_no > cellmax:
            max_sc_dim = max(sc_scale)
            i = sc_scale.index(max_sc_dim)
            sc_scale[i] -= 1
            scs = vac.make_supercells_with_defects(sc_scale)

    interdir = mpid
    blk_str_sites = set(scs[0].sites)
    for i in range(1,len(scs)):
        sc = scs[i]
        vac_str_sites = set(sc.sites)
        vac_sites = blk_str_sites - vac_str_sites
        vac_site = list(vac_sites)[0]
        site_mult = int(vac.get_defectsite_multiplicity(i-1)/conv_prim_rat)
        vac_site_specie = vac_site.specie
        vac_specie = vac_site.specie.symbol

        # Solute substitution defect generation at all vacancy sites
        struct_species = scs[0].types_of_specie
        solute_struct = sc.copy()
        solute_struct.append(solute, vac_site.frac_coords)

        mpvis = MPRelaxSet(solute_struct, user_incar_settings={"LDAU": False})

        incar = mpvis.incar
        incar.update(def_vasp_incar_param)
        poscar = mpvis.poscar
        kpoints = Kpoints.automatic_density(solute_struct,kpoint_den)
        if ptcr_flag:
            potcar = mpvis.potcar

        sub_def_dir ='solute_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                str(i), site_mult, vac_specie, solute)
        fin_dir = os.path.join(interdir,sub_def_dir)
        try:
            os.makedirs(fin_dir)
        except:
            pass
        poscar.write_file(os.path.join(fin_dir,'POSCAR'))
        incar.write_file(os.path.join(fin_dir,'INCAR'))
        kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))
        if ptcr_flag:
            potcar.write_file(os.path.join(fin_dir,'POTCAR'))


def solute_def_parse_energy(args):
    mpid = args.mpid 
    solute = args.solute 
    mapi_key = args.mapi_key 

    if not mpid:
        print ("============\nERROR: Provide an mpid\n============")
        return 
    if not solute:
        print ("============\nERROR: Provide solute element\n============")
        return 

    if not mapi_key:
        with MPRester() as mp:
            structure = mp.get_structure_by_material_id(mpid)      
    else:
        with MPRester(mapi_key) as mp:
            structure = mp.get_structure_by_material_id(mpid)      

    energy_dict = {}

    solutes = []
    def_folders = glob.glob(os.path.join(
        mpid,"solute*subspecie-{}".format(solute)))
    def_folders += glob.glob(os.path.join(mpid,"bulk"))
    for defdir in def_folders:
        fldr_name = os.path.split(defdir)[1]
        vr_file = os.path.join(defdir,'vasprun.xml') 
        if not os.path.exists(vr_file):
            print (fldr_name, ": vasprun.xml doesn't exist in the folder. " \
                   "Abandoning parsing of energies for {}".format(mpid))
            break       # Further processing for the mpid is not useful

        try:
            vr = Vasprun(vr_file)
        except:
            print (fldr_name, ":Failure, couldn't parse vaprun.xml file. "
                   "Abandoning parsing of energies for {}".format(mpid))
            break

        if not vr.converged:
            print (fldr_name, ": Vasp calculation not converged. "
                   "Abandoning parsing of energies for {}".format(mpid))
            break       # Further processing for the mpid is not useful

        fldr_fields = fldr_name.split("_")
        if 'bulk' in fldr_fields:
            bulk_energy = vr.final_energy
            bulk_sites = vr.structures[-1].num_sites
        elif 'solute' in fldr_fields:
            site_index = int(fldr_fields[1])
            site_multiplicity = int(fldr_fields[2].split("-")[1])
            site_specie = fldr_fields[3].split("-")[1]
            substitution_specie = fldr_fields[4].split("-")[1]
            energy = vr.final_energy
            solutes.append({'site_index':site_index,
                'site_specie':site_specie,'energy':energy,
                'substitution_specie':substitution_specie,
                'site_multiplicity':site_multiplicity
                })
    else:
        if not solutes:
            print("Solute folders do not exist")
            return {}

        print("Solute {} calculations successful for {}".format(solute,mpid))
        for solute in solutes:
            solute_flip_energy = solute['energy']-bulk_energy
            solute['energy'] = solute_flip_energy
        solutes.sort(key=lambda entry: entry['site_index'])
        energy_dict[mpid] = {'solutes':solutes}
        fl_nm = mpid+'_solute-'+args.solute+'_raw_defect_energy.json'
        dumpfn(energy_dict, fl_nm, indent=2, cls=MontyEncoder)


def vac_antisite_def_parse_energy(args):
    mpid = args.mpid
    mapi_key = args.mapi_key 

    if not mpid:
        print("============\nERROR: Provide an mpid\n============")
        return 

    if not mapi_key:
        with MPRester() as mp:
            structure = mp.get_structure_by_material_id(mpid)      
    else:
        with MPRester(mapi_key) as mp:
            structure = mp.get_structure_by_material_id(mpid)      

    energy_dict = {}

    antisites = []
    vacancies = []
    def_folders = glob.glob(os.path.join(mpid,"vacancy*"))
    def_folders += glob.glob(os.path.join(mpid,"antisite*"))
    def_folders += glob.glob(os.path.join(mpid,"bulk"))
    for defdir in def_folders:
        fldr_name = os.path.split(defdir)[1]
        vr_file = os.path.join(defdir,'vasprun.xml') 
        if not os.path.exists(vr_file):
            print (fldr_name, ": vasprun.xml doesn't exist in the folder. " \
                   "Abandoning parsing of energies for {}".format(mpid))
            break       # Further processing for the mpid is not useful

        try:
            vr = Vasprun(vr_file)
        except:
            print (fldr_name, ":Failure, couldn't parse vaprun.xml file. "
                   "Abandoning parsing of energies for {}".format(mpid))
            break

        if not vr.converged:
            print (fldr_name, ": Vasp calculation not converged. "
                   "Abandoning parsing of energies for {}".format(mpid))
            break       # Further processing for the mpid is not useful

        fldr_fields = fldr_name.split("_")
        if 'bulk' in fldr_fields:
            bulk_energy = vr.final_energy
            bulk_sites = vr.structures[-1].num_sites
        elif 'vacancy' in fldr_fields:
            site_index = int(fldr_fields[1])
            site_multiplicity = int(fldr_fields[2].split("-")[1])
            site_specie = fldr_fields[3].split("-")[1]
            energy = vr.final_energy
            vacancies.append({'site_index':site_index,
                'site_specie':site_specie,'energy':energy,
                'site_multiplicity':site_multiplicity
                })
        elif 'antisite' in fldr_fields:
            site_index = int(fldr_fields[1])
            site_multiplicity = int(fldr_fields[2].split("-")[1])
            site_specie = fldr_fields[3].split("-")[1]
            substitution_specie = fldr_fields[4].split("-")[1]
            energy = vr.final_energy
            antisites.append({'site_index':site_index,
                'site_specie':site_specie,'energy':energy,
                'substitution_specie':substitution_specie,
                'site_multiplicity':site_multiplicity
                })
    else:
        print("All calculations successful for ", mpid)
        e0 = bulk_energy/bulk_sites*structure.num_sites
        for vac in vacancies:
            vac_flip_energy = vac['energy']-bulk_energy
            vac['energy'] = vac_flip_energy
        vacancies.sort(key=lambda entry: entry['site_index'])
        for antisite in antisites:
            as_flip_energy = antisite['energy']-bulk_energy
            antisite['energy'] = as_flip_energy
        antisites.sort(key=lambda entry: entry['site_index'])
        energy_dict[str(mpid)] = {u"structure":structure,
                'e0':e0,'vacancies':vacancies,'antisites':antisites}

        fl_nm = args.mpid+'_raw_defect_energy.json'
        dumpfn(energy_dict, fl_nm, cls=MontyEncoder, indent=2)


def get_def_profile(args):
    if not args.mpid and not args.file:
        print ("------------\nERROR: mpid is not given.\n========")
        return
    mpid = args.mpid 
    T = args.T
    if args.file:
        file = args.file 
    else:
        file = mpid+'_raw_defect_energy.json'

    raw_energy_dict = loadfn(file,cls=MontyDecoder)

    e0 = raw_energy_dict[mpid]['e0']
    struct = raw_energy_dict[mpid]['structure']
    vacs = raw_energy_dict[mpid]['vacancies']
    antisites = raw_energy_dict[mpid]['antisites']
    vacs.sort(key=lambda entry: entry['site_index'])
    antisites.sort(key=lambda entry: entry['site_index'])
    for vac_def in vacs:
        if not vac_def:
            print('All vacancy defect energies not present')
            continue
    for antisite_def in antisites:
        if not antisite_def:
            print('All antisite defect energies not preset')
            continue

    try:
        def_conc, def_en, mu = compute_defect_density(struct, e0, vacs, antisites, T,
                plot_style='gnuplot')
        fl_nm = mpid+'_def_concentration.dat'
        with open(fl_nm,'w') as fp:
            for row in def_conc:
                print >> fp, row
        fl_nm = mpid+'_def_energy.dat'
        with open(fl_nm,'w') as fp:
            for row in def_en:
                print >> fp, row
        fl_nm = mpid+'_chem_pot.dat'
        with open(fl_nm,'w') as fp:
            for row in mu:
                print >> fp, row
    except:
        raise


def get_solute_def_profile(args):
    if not args.mpid:
        print ('===========\nERROR: mpid is not given.\n===========')
        return
    if not args.solute:
        print ('===========\nERROR: Solute atom is not given.\n===========')
        return

    mpid = args.mpid 
    solute = args.solute 
    solute_conc = args.solute_conc/100.0
    T = args.T 

    def_file = mpid + '_raw_defect_energy.json'
    raw_energy_dict = loadfn(def_file,cls=MontyDecoder)
    sol_file = mpid+'_solute-'+solute+'_raw_defect_energy.json'
    sol_raw_energy_dict = loadfn(sol_file,cls=MontyDecoder)

    #try:
    e0 = raw_energy_dict[mpid]['e0']
    struct = raw_energy_dict[mpid]['structure']
    vacs = raw_energy_dict[mpid]['vacancies']
    antisites = raw_energy_dict[mpid]['antisites']
    solutes = sol_raw_energy_dict[mpid]['solutes']

    for vac_def in vacs:
        if not vac_def:
            print('All vacancy defect energies not present')
            continue
    for antisite_def in antisites:
        if not antisite_def:
            print('All antisite defect energies not preset')
            continue
    for solute_def in solutes:
        if not solute_def:
            print('All solute defect energies not preset')
            continue

    try:
        def_conc = solute_defect_density(struct, e0, vacs, 
                antisites, solutes, solute_concen=solute_conc, T=T, 
                plot_style="gnuplot")
        fl_nm = args.mpid+'_solute-'+args.solute+'_def_concentration.dat'
        with open(fl_nm,'w') as fp: 
            for row in def_conc:
                print >> fp, row 
    except:
        raise


def main():
    parser = argparse.ArgumentParser(description="""
    pydii is a script that generates vasp inputs, parses vasp output files
    and computes the point defect concentrations. This script works based 
    on several sub-commands with their own options. To see the options for 
    the sub-commands, type "pydii sub-command -h".""",
                                     epilog="""
    Author: Bharat Medasani
    Version: {}
    Last updated: {}""".format(__version__, __date__))

    subparsers = parser.add_subparsers()
    MP_string = "Materials Project id of the intermetallic structure.\n" \
                "For more info on Materials Project, please refer to " \
                "www.materialsproject.org"
    MAPI_string = "Your Materials Project REST API key.\nFor more info, " \
                  "please refer to www.materialsproject.org/opne"
    cell_string = "Maximum number of atoms in supercell.\n" \
                  "The default is 128.\n" \
                  "Keep in mind the number of atoms in the supercell " \
                  "may vary from the provided number including the default."

    parser_vasp_inp = subparsers.add_parser("gen_def_structure", 
            help="Vasp input files for intrinsic defects.")
    parser_vasp_inp.add_argument("--mpid", type=str.lower, dest="mpid",
            help=MP_string)
    parser_vasp_inp.add_argument("--mapi_key", default = None, 
            dest="mapi_key", help=MAPI_string)
    parser_vasp_inp.add_argument("--cellmax", type=int, default=128, 
            dest="cellmax", help=cell_string)
    parser_vasp_inp.set_defaults(func=vac_antisite_def_struct_gen)


    parser_sol_vasp_inp = subparsers.add_parser("gen_sol_pref_structure", 
            help="Vasp input files for extrinsic substitutional defects.")
    parser_sol_vasp_inp.add_argument("--mpid", type=str.lower, dest="mpid",
            help=MP_string)
    parser_sol_vasp_inp.add_argument("--solute", dest="solute", 
            help="Solute Element")
    parser_sol_vasp_inp.add_argument("--mapi_key", default = None, 
            dest="mapi_key", help=MAPI_string)
    parser_sol_vasp_inp.add_argument("--cellmax", type=int, default=128, 
            dest="cellmax", help=cell_string)
    parser_sol_vasp_inp.set_defaults(func=substitute_def_struct_gen)


    parser_vasp_out = subparsers.add_parser("gen_def_energy", 
            help = 'Command to parse vacancy and antisite defect ' \
                   'energies for intermetallics from the VASP DFT ' \
                   'calculations.')
    parser_vasp_out.add_argument("--mpid", type=str.lower, dest="mpid",
            help=MP_string)
    parser_vasp_out.add_argument("--mapi_key", default = None, dest="mapi_key",
            help=MAPI_string)
    parser_vasp_out.set_defaults(func=vac_antisite_def_parse_energy)


    parser_sol_vasp_out = subparsers.add_parser("gen_sol_def_energy", 
            help = 'Command to parse solute substitution defect ' \
                   'energies for intermetallics from the VASP DFT ' \
                   'calculations.')
    parser_sol_vasp_out.add_argument("--mpid", type=str.lower, dest="mpid",
            help=MP_string)
    parser_sol_vasp_out.add_argument("--solute", dest="solute", 
            help="Solute Element")
    parser_sol_vasp_out.add_argument("--mapi_key", default = None, 
            dest="mapi_key", help=MAPI_string)
    parser_sol_vasp_out.set_defaults(func=solute_def_parse_energy)


    parser_conc = subparsers.add_parser("gen_def_profile", 
            help = 'Command to generate vacancy and antisite defect ' \
                   'concentration for intermetallics from the raw defect '\
                   'energies.' )
    parser_conc.add_argument("--mpid", type=str.lower, dest="mpid",
            help=MP_string)
    parser_conc.add_argument('-T', "--temp", type=float, default=1000,
            dest="T", help="Temperature in Kelvin")
    parser_conc.add_argument("--file", default = None, dest="file",
            help = "The default file is 'mpid'+'_raw_defect_energy.json'.\n" \
                   "If the file is named differently supply it.")
    parser_conc.set_defaults(func=get_def_profile)


    parser_sol_conc = subparsers.add_parser("gen_sol_site_pref", 
            help = 'Command to generate solute defect site preference ' \
                   'for intermetallics from the raw defect energies.' )
    parser_sol_conc.add_argument("--mpid", type=str.lower, dest="mpid",
            help=MAPI_string)
    parser_sol_conc.add_argument("--solute", dest="solute", 
            help="Solute Element")
    parser_sol_conc.add_argument("--sol_conc", type=float, default=1.0,
            dest="solute_conc", 
            help="Solute Concentration in %. Default is 1%")
    parser_sol_conc.add_argument('-T', "--temp", type=float, default=1000.0,
            dest="T", help="Temperature in Kelvin")
    parser_sol_conc.set_defaults(func=get_solute_def_profile)


    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
