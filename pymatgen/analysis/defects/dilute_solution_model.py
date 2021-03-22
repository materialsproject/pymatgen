# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Evaluate the defect concentration based on composition, temperature,
and defect energies using "Dilute Solution Model"
Reference: Phys Rev B, 63, 094103, 2001,
"Density of constitutional and thermal point defects in L12 Al3Sc",
C. Woodward, M. Asta, G. Kresse and J. Hafner.
Manual and citation for the code, DOI: 10.1016/j.cpc.2015.03.015
"""

import copy
import math

import numpy as np
from monty.dev import deprecated
from monty.fractions import gcd
from sympy import Float, Integer, Matrix, Symbol, exp, nsolve, solve

__author__ = "Bharat Medasani"
__version__ = "0.2"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Alpha"
__date__ = "6/4/14"

# physical consts
k_B = 8.6173324e-5  # eV/K


# Check the inputs
def _check_input(def_list):
    flag = True
    for defect in def_list:
        if not defect:
            flag = False
            break
    return flag


@deprecated(message="Refactoring of the defects module will eventualy remove this function")
def dilute_solution_model(structure, e0, vac_defs, antisite_defs, T, trial_chem_pot=None, generate="plot"):
    """
    Compute the defect densities using dilute solution model.

    Args:
        structure: pymatgen.core.structure.Structure object representing the
            primitive or unitcell of the crystal.
        e0: The total energy of the undefected system.
            This is E0 from VASP calculation.
        vac_defs: List of vacancy defect parameters in the dictionary format.
            The keys of the dict associated with each vacancy defect are
            1) site_index, 2) site_specie, 3) site_multiplicity, and
            4) energy. 1-3 can be obtained from
            pymatgen.analysis.defects.point_defects.Vacancy class.
            Site index is expected to start with 1 (fortran index).
        antisite_defs: List of antisite defect parameters in the dictionary
            format. The keys of the dict associated with each antisite defect
            are 1) site_index, 2) site_specie, 3) site_multiplicity,
            4) substitution_specie, and 5) energy. 1-3 can be obtained
            from pymatgen.analysis.defects.point_defects.Vacancy class.
        T: Temperature in Kelvin
        trial_chem_pot (optional): Trial chemical potentials to speedup
            the plot generation. Format is {el1:mu1,...}
        generate (string): Options are plot or energy
            Chemical potentials are also returned with energy option.
            If energy option is not chosen, plot is generated.

    Returns:
        If generate=plot, the plot data is generated and returned in
        HighCharts format.
        If generate=energy, defect formation enthalpies and chemical
        potentials are returned.
    """

    if not _check_input(vac_defs):
        raise ValueError("Vacancy energy is not defined")
    if not _check_input(antisite_defs):
        raise ValueError("Antisite energy is not defined")

    formation_energies = {}
    formation_energies["vacancies"] = copy.deepcopy(vac_defs)
    formation_energies["antisites"] = copy.deepcopy(antisite_defs)
    for vac in formation_energies["vacancies"]:
        del vac["energy"]
    for asite in formation_energies["antisites"]:
        del asite["energy"]
    # Setup the system
    site_species = [vac_def["site_specie"] for vac_def in vac_defs]
    multiplicity = [vac_def["site_multiplicity"] for vac_def in vac_defs]
    m = len(set(site_species))  # distinct species
    n = len(vac_defs)  # inequivalent sites

    # Reduce the system and associated parameters such that only distinctive
    # atoms are retained
    comm_div = gcd(*tuple(multiplicity))
    multiplicity = [val / comm_div for val in multiplicity]
    e0 = e0 / comm_div
    T = Float(T)

    # c0 = np.diag(multiplicity)
    c0 = np.diag(np.ones(n))
    mu = [Symbol("mu" + i.__str__()) for i in range(m)]

    # Generate maps for hashing
    # Generate specie->mu map and use it for site->mu map
    specie_order = []  # Contains hash for site->mu map    Eg: [Al, Ni]
    site_specie_set = set()  # Eg: {Ni, Al}
    for i in range(n):
        site_specie = site_species[i]
        if site_specie not in site_specie_set:
            site_specie_set.add(site_specie)
            specie_order.append(site_specie)
    site_mu_map = []  # Eg: [mu0,mu0,mu0,mu1] where mu0->Al, and mu1->Ni
    for i in range(n):
        site_specie = site_species[i]
        j = specie_order.index(site_specie)
        site_mu_map.append(j)
    specie_site_index_map = []  # Eg: [(0,3),(3,4)] for Al & Ni
    for i in range(m):
        low_ind = site_species.index(specie_order[i])
        if i < m - 1:
            hgh_ind = site_species.index(specie_order[i + 1])
        else:
            hgh_ind = n
        specie_site_index_map.append((low_ind, hgh_ind))
    """
    dC: delta concentration matrix:
    dC[i,j,k]: Concentration change of atom i, due to presence of atom
    j on lattice site k
    Special case is [i,i,i] which is considered as vacancy
    Few cases: dC[i,i,i] = -1 due to being vacancy special case
                dC[k,k,i] = +1 due to increment in k at i lattice if i
                               lattice type is of different element
                dC[i,k,i] = -1 due to decrement of ith type atom due to
                presence of kth type atom on ith sublattice and kth type
                atom specie is different from ith sublattice atom specie
                dC[i,k,k] = 0 due to no effect on ith type atom
                dC[i,j,k] = 0 if i!=j!=k
    """
    dC = np.zeros((n, n, n), dtype=np.int_)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if i == j and site_species[j] != site_species[k] and site_species[i] != site_species[k]:
                    dC[i, j, k] = 1
        for j in range(n):
            for k in range(n):
                if i == k:
                    dC[i, j, k] = -1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                if i != j:
                    if site_species[j] == site_species[k]:
                        dC[i, j, k] = 0

    for ind_map in specie_site_index_map:
        if ind_map[1] - ind_map[0] > 1:
            for index1 in range(ind_map[0] + 1, ind_map[1]):
                for index2 in range(ind_map[0]):
                    for i in range(n):
                        dC[i, index1, index2] = 0
                for index2 in range(ind_map[1], n):
                    for i in range(n):
                        dC[i, index1, index2] = 0

    # dE matrix: Flip energies (or raw defect energies)
    els = [vac_def["site_specie"] for vac_def in vac_defs]
    dE = []
    for i in range(n):
        dE.append([])
    for i in range(n):
        for j in range(n):
            dE[i].append(0)

    for j in range(n):
        for i in range(n):
            if i == j:
                dE[i][j] = vac_defs[i]["energy"]
            else:
                sub_specie = vac_defs[i]["site_specie"]
                site_specie = vac_defs[j]["site_specie"]
                if site_specie == sub_specie:
                    dE[i][j] = 0
                else:
                    for as_def in antisite_defs:
                        if int(as_def["site_index"]) == j + 1 and sub_specie == as_def["substitution_specie"]:
                            dE[i][j] = as_def["energy"]
                            break
    dE = np.array(dE)

    # Initialization for concentrations
    # c(i,p) == presence of ith type atom on pth type site
    c = Matrix(n, n, [0] * n ** 2)
    for i in range(n):
        for p in range(n):
            c[i, p] = Integer(c0[i, p])
            site_flip_contribs = []
            for epi in range(n):
                sum_mu = sum([mu[site_mu_map[j]] * Integer(dC[j, epi, p]) for j in range(n)])
                flip = Integer(dC[i, epi, p]) * exp(-(dE[epi, p] - sum_mu) / (k_B * T))
                if flip not in site_flip_contribs:
                    site_flip_contribs.append(flip)
                    c[i, p] += flip

    total_c = []
    for ind in specie_site_index_map:
        val = 0
        for i in range(*ind):
            sum_i = sum([c[i, j] * multiplicity[j] for j in range(n)])
            val += sum_i
        total_c.append(val)

    c_ratio = [total_c[-1] / total_c[i] for i in range(m)]

    # Expression for Omega, the Grand Potential
    omega1 = e0 - sum([mu[site_mu_map[i]] * sum(c0[i, :]) * multiplicity[i] for i in range(n)])
    omega2 = []
    fm_en_eff = []
    used_dEs = []
    for p_r in range(n):
        for epi in range(n):
            sum_mu = sum([mu[site_mu_map[j]] * dC[j, epi, p_r] for j in range(n)])
            if p_r != epi and site_mu_map[p_r] == site_mu_map[epi]:
                continue
            if dE[epi, p_r] not in used_dEs:
                omega2.append(k_B * T * multiplicity[p_r] * exp(-(dE[epi, p_r] - sum_mu) / (k_B * T)))
                fm_en_eff.append(dE[epi, p_r] - sum_mu)
                used_dEs.append(dE[epi, p_r])
    omega = omega1 - sum(omega2)

    # Compute composition range
    li = specie_site_index_map[0][0]
    hi = specie_site_index_map[0][1]
    comp1_min = sum(multiplicity[li:hi]) / sum(multiplicity) * 100 - 1
    comp1_max = sum(multiplicity[li:hi]) / sum(multiplicity) * 100 + 1
    delta = float(comp1_max - comp1_min) / 120.0
    yvals = []
    for comp1 in np.arange(comp1_min, comp1_max + delta, delta):
        comp2 = 100 - comp1
        y = comp2 / comp1
        yvals.append(y)

    def reduce_mu():
        omega = [e0 - sum([mu[site_mu_map[i]] * sum(c0[i, :]) for i in range(n)])]
        x = solve(omega)
        return x

    def compute_mus_by_search():
        # Compute trial mu
        mu_red = reduce_mu()

        mult = multiplicity
        specie_concen = [sum(mult[ind[0] : ind[1]]) for ind in specie_site_index_map]
        y_vect = [specie_concen[-1] / specie_concen[i] for i in range(m)]
        vector_func = [y_vect[i] - c_ratio[i] for i in range(m - 1)]
        vector_func.append(omega)
        min_diff = 1e10
        mu_vals = None
        c_val = None
        m1_min = -20.0
        if e0 > 0:
            m1_max = 10  # Search space needs to be modified
        else:
            m1_max = 0
        for m1 in np.arange(m1_min, m1_max, 0.01):
            m0 = mu_red[mu[0]].subs(mu[-1], m1)

            try:
                x = nsolve(vector_func, mu, [m0, m1], module="numpy")
            except Exception:
                continue

            c_val = c.subs(dict(zip(mu, x)))
            # if all(x >= 0 for x in c_val):
            specie_concen = []
            for ind in specie_site_index_map:
                specie_concen.append(sum([sum(c_val[i, :]) for i in range(*ind)]))
            y_comp = [specie_concen[-1] / specie_concen[i] for i in range(m)]
            diff = math.sqrt(sum([pow(abs(y_comp[i] - y_vect[i]), 2) for i in range(m)]))
            if diff < min_diff:
                min_diff = diff
                mu_vals = x

        if mu_vals:
            mu_vals = [float(mu_val) for mu_val in mu_vals]
        else:
            raise ValueError()
        return mu_vals

    def compute_def_formation_energies():
        i = 0
        for vac_def in vac_defs:
            site_specie = vac_def["site_specie"]
            ind = specie_order.index(site_specie)
            uncor_energy = vac_def["energy"]
            formation_energy = uncor_energy + mu_vals[ind]
            formation_energies["vacancies"][i]["formation_energy"] = formation_energy
            specie_ind = site_mu_map[i]
            indices = specie_site_index_map[specie_ind]
            specie_ind_del = indices[1] - indices[0]
            cur_ind = i - indices[0] + 1
            if not specie_ind_del - 1:
                label = "$V_{" + site_specie + "}$"
            else:
                label = "$V_{" + site_specie + "_" + str(cur_ind) + "}$"
            formation_energies["vacancies"][i]["label"] = label
            i += 1
        i = 0
        for as_def in antisite_defs:
            site_specie = as_def["site_specie"]
            sub_specie = as_def["substitution_specie"]
            ind1 = specie_order.index(site_specie)
            ind2 = specie_order.index(sub_specie)
            uncor_energy = as_def["energy"]
            formation_energy = uncor_energy + mu_vals[ind1] - mu_vals[ind2]
            formation_energies["antisites"][i]["formation_energy"] = formation_energy
            specie_ind = site_mu_map[i]
            indices = specie_site_index_map[specie_ind]
            specie_ind_del = indices[1] - indices[0]
            cur_ind = i - indices[0] + 1
            if not specie_ind_del - 1:
                label = "$" + sub_specie + "_{" + site_specie + "}$"
            else:
                label = "$" + sub_specie + "_{" + site_specie + "_" + str(cur_ind) + "}$"
            formation_energies["antisites"][i]["label"] = label
            i += 1
        return formation_energies

    # If generate option is energy compute effective formation energies
    # at ideal stoichiometry and return the formation energies and chem pot.
    if generate == "energy":
        if not trial_chem_pot:
            mu_vals = compute_mus_by_search()
        else:
            mu_vals = [trial_chem_pot[element] for element in specie_order]

        formation_energies = compute_def_formation_energies()
        mu_dict = dict(zip(specie_order, mu_vals))
        return formation_energies, mu_dict

    if not trial_chem_pot:
        # Try computing mus by assuming one of the defects is dominant at 0.01
        # concen.  First vacancy is tried and then antisite

        # Generate trial mus assuming vacancy as dominant defect
        # for specie-0 at lower yval
        li = specie_site_index_map[0][0]
        hi = specie_site_index_map[0][1]
        li1 = specie_site_index_map[1][0]
        hi1 = specie_site_index_map[1][1]
        spec_mult = [sum(multiplicity[li:hi]), sum(multiplicity[li1:hi1])]
        ln_def_conc = 4.60517
        for i in range(li, hi):
            vac_flip_en = vac_defs[i]["energy"]
            mu_vals = [ln_def_conc * k_B * T - vac_flip_en]
            mu_vals.append((e0 - spec_mult[0] * mu_vals[0]) / spec_mult[1])
            comp_ratio = yvals[0]

            # Test if the trial mus are good
            vector_func = [comp_ratio - c_ratio[0]]
            vector_func.append(omega)
            try:
                mu_vals = nsolve(vector_func, mu, mu_vals)
                if mu_vals:
                    mu_vals = [float(mu_val) for mu_val in mu_vals]
                break
            except Exception:  # Go for antisite as dominant defect
                mu_gs = [Symbol("mu_gs" + j.__str__()) for j in range(m)]

                eqs = [mu_gs[0] - mu_gs[1] - (ln_def_conc * k_B * T - antisite_defs[i]["energy"])]
                eqs.append(spec_mult[0] * mu_gs[0] + spec_mult[1] * mu_gs[1] - e0)
                x = solve(eqs, mu_gs)
                # mu_names = sorted([key.name for key in x.keys()])
                mu_vals = []
                for key in sorted(x.keys(), key=lambda inp: inp.name):
                    mu_vals.append(x[key])
                vector_func = [comp_ratio - c_ratio[0]]
                vector_func.append(omega)

                try:
                    mu_vals = nsolve(vector_func, mu, mu_vals)
                    if mu_vals:
                        mu_vals = [float(mu_val) for mu_val in mu_vals]
                    break
                except Exception:  # Go to the default option (search the space)
                    pass
        else:
            mu_vals = compute_mus_by_search()

    else:
        try:
            mu_vals = [trial_chem_pot[element] for element in specie_order]
        except Exception:
            mu_vals = compute_mus_by_search()

    # Compile mu's for all composition ratios in the range
    # +/- 1% from the stoichiometry
    result = {}
    i = 0
    failed_y, failed_i = [], []
    for y in yvals:
        vector_func = [y - c_ratio[0]]
        vector_func.append(omega)
        try:
            x = nsolve(vector_func, mu, mu_vals, module="numpy")
            if x:
                mu_vals = [float(mu_val) for mu_val in x]
        except Exception:
            failed_y.append(y)
            failed_i.append(i)
            continue
        result[y] = list(mu_vals)
        x = None
        i += 1

    def get_next_mu_val(i):
        if i >= len(yvals):
            return None

        y = yvals[i + 1]
        x = result.get(y, None)
        if x:
            mu_vals = [float(mu_val) for mu_val in x]
            return mu_vals
        return get_next_mu_val(i + 1)

    def get_prev_mu_val(i):
        if i <= 0:
            return None

        y = yvals[i - 1]
        x = result.get(y, None)
        if x:
            mu_vals = [float(mu_val) for mu_val in x]
            return mu_vals
        return get_next_mu_val(i - 1)

    # Try to get better trial mus for failed cases
    for j, y in enumerate(failed_y):
        i = failed_i[j]

        prev_mu_val = get_prev_mu_val(i)
        if not prev_mu_val:
            continue
        next_mu_val = get_next_mu_val(i)
        if not next_mu_val:
            continue

        vector_func = [y - c_ratio[0]]
        vector_func.append(omega)
        trial_mu = list(map(lambda x: float(sum(x)) / len(x), zip(prev_mu_val, next_mu_val)))
        try:
            x = nsolve(vector_func, mu, trial_mu, module="numpy")
            if x:
                mu_vals = [float(mu_val) for mu_val in x]
        except Exception:
            continue
        result[y] = mu_vals
        x = None

    # Alternate way of calculating trial mus for failed cases
    # by taking average of trial mus at extremes.
    # for j in range(len(failed_y)):
    #    y = yvals[0]
    #    prev_mu_val = result[y]
    #    y = yvals[-1]
    #    next_mu_val = result[y]
    #
    #    trial_mu = list(map(lambda x: float(sum(x))/len(x),  \
    #            zip(prev_mu_val,next_mu_val)))
    #    y = failed_y[j]
    #    vector_func = [y-c_ratio[0]]
    #    vector_func.append(omega)
    #    try:
    #        x = nsolve(vector_func,mu,trial_mu,module="numpy")
    #        if x:
    #            mu_vals = [float(mu_val) for mu_val in x]
    #    except Exception:
    #        continue
    #    result[y] = list(mu_vals)

    if len(result.keys()) < len(yvals) / 2:
        raise ValueError("Not sufficient data")

    res = []
    new_mu_dict = {}
    # Compute the concentrations for all the compositions
    for key in sorted(result.keys()):
        mu_val = result[key]
        total_c_val = [total_c[i].subs(dict(zip(mu, mu_val))) for i in range(len(total_c))]
        c_val = c.subs(dict(zip(mu, mu_val)))
        res1 = []
        # Concentration of first element/over total concen
        res1.append(float(total_c_val[0] / sum(total_c_val)))
        new_mu_dict[res1[0]] = mu_val
        for i in range(n):
            for j in range(n):
                if i == j:  # Vacancy
                    vac_conc = float(exp(-(mu_val[site_mu_map[i]] + dE[i, i]) / (k_B * T)))
                    res1.append(vac_conc)
                else:  # Antisite
                    res1.append(float(c_val[i, j] / c0[j, j]))
        res.append(res1)

    res = np.array(res)
    dtype = [(str("x"), np.float64)] + [(str("y%d%d" % (i, j)), np.float64) for i in range(n) for j in range(n)]
    res1 = np.sort(res.view(dtype), order=[str("x")], axis=0)

    conc_data = {}
    """Because all the plots have identical x-points storing it in a
    single array"""
    conc_data["x"] = [dat[0][0] for dat in res1]  # x-axis data
    # Element whose composition is varied. For x-label
    conc_data["x_label"] = els[0] + " mole fraction"
    conc_data["y_label"] = "Point defect concentration"
    conc = []
    for i in range(n):
        conc.append([])
        for j in range(n):
            conc[i].append([])
    for i in range(n):
        for j in range(n):
            y1 = [dat[0][i * n + j + 1] for dat in res1]
            conc[i][j] = y1

    y_data = []
    for i in range(n):
        data = conc[i][i]
        specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        vac_string = "$Vac_{"
        if not specie_ind_del - 1:
            label = vac_string + specie + "}$"
        else:
            label = vac_string + specie + "_" + str(cur_ind) + "}$"
        # Plot data and legend info
        y_data.append({"data": data, "name": label})

    for i in range(n):
        site_specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        for j in range(m):  # Antisite plot dat
            sub_specie = specie_order[j]
            if sub_specie == site_specie:
                continue
            if not specie_ind_del - 1:
                label = "$" + sub_specie + "_{" + site_specie + "}$"
            else:
                label = "$" + sub_specie + "_{" + site_specie + "_" + str(cur_ind) + "}$"
            inds = specie_site_index_map[j]
            # TODO: Investigate the value below
            data = np.sum([conc[ind][i] for ind in range(*inds)], axis=0)
            data = data.tolist()
            y_data.append({"data": data, "name": label})

    conc_data["y"] = y_data

    # Compute the  formation energies
    def compute_vac_formation_energies(mu_vals):
        en = []
        for vac_def in vac_defs:
            site_specie = vac_def["site_specie"]
            ind = specie_order.index(site_specie)
            uncor_energy = vac_def["energy"]
            formation_energy = uncor_energy + mu_vals[ind]
            en.append(float(formation_energy))
        return en

    en_res = []
    for key in sorted(new_mu_dict.keys()):
        mu_val = new_mu_dict[key]
        en_res.append(compute_vac_formation_energies(mu_val))

    en_data = {"x_label": els[0] + " mole fraction", "x": []}
    en_data["x"] = [dat[0][0] for dat in res1]  # x-axis data

    i = 0
    y_data = []
    for vac_def in vac_defs:
        data = [data[i] for data in en_res]
        site_specie = vac_def["site_specie"]
        ind = specie_order.index(site_specie)
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        vac_string = "$Vac_{"
        if not specie_ind_del - 1:
            label = vac_string + site_specie + "}$"
        else:
            label = vac_string + site_specie + "_" + str(cur_ind) + "}$"
        y_data.append({"data": data, "name": label})
        i += 1

    def compute_as_formation_energies(mu_vals):
        en = []
        for as_def in antisite_defs:
            site_specie = as_def["site_specie"]
            sub_specie = as_def["substitution_specie"]
            ind1 = specie_order.index(site_specie)
            ind2 = specie_order.index(sub_specie)
            uncor_energy = as_def["energy"]
            form_en = uncor_energy + mu_vals[ind1] - mu_vals[ind2]
            en.append(form_en)
        return en

    en_res = []
    for key in sorted(new_mu_dict.keys()):
        mu_val = new_mu_dict[key]
        en_res.append(compute_as_formation_energies(mu_val))
    i = 0
    for as_def in antisite_defs:
        data = [data[i] for data in en_res]
        site_specie = as_def["site_specie"]
        sub_specie = as_def["substitution_specie"]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        if not specie_ind_del - 1:
            label = "$" + sub_specie + "_{" + site_specie + "}$"
        else:
            label = "$" + sub_specie + "_{" + site_specie + "_" + str(cur_ind) + "}$"
        y_data.append({"data": data, "name": label})
        i += 1

    en_data["y"] = y_data

    # Return chem potential as well
    mu_data = {"x_label": els[0] + " mole fraction", "x": []}
    mu_data["x"] = [dat[0][0] for dat in res1]  # x-axis data

    y_data = []
    for j in range(m):
        specie = specie_order[j]
        mus = [new_mu_dict[key][j] for key in sorted(new_mu_dict.keys())]
        y_data.append({"data": mus, "name": specie})
    mu_data["y"] = y_data

    return conc_data, en_data, mu_data


def compute_defect_density(
    structure,
    e0,
    vac_defs,
    antisite_defs,
    T=800,
    trial_chem_pot=None,
    plot_style="highcharts",
):
    """
    Wrapper for the dilute_solution_model.
    The computed plot data is prepared based on plot_style.

    Args:
        structure: pymatgen.core.structure.Structure object representing the
            primitive or unitcell of the crystal.
        e0: The total energy of the undefected system.
            This is E0 from VASP calculation.
        vac_defs: List of vacancy defect parameters in the dictionary format.
            The keys of the dict associated with each vacancy defect are
            1) site_index, 2) site_specie, 3) site_multiplicity, and
            4) energy. 1-3 can be obtained from
            pymatgen.analysis.defects.point_defects.Vacancy class.
            Site index is expected to start with 1 (fortran index).
        antisite_defs: List of antisite defect parameters in the dictionary
            format. The keys of the dict associated with each antisite defect
            are 1) site_index, 2) site_specie, 3) site_multiplicity,
            4) substitution_specie, and 5) energy. 1-3 can be obtained
            from pymatgen.analysis.defects.point_defects.Vacancy class.
        T: Temperature in Kelvin
        trial_chem_pot (optional): Trial chemical potentials to speedup
            the plot generation. Format is {el1:mu1,...}
        plot_style (string): Allowed options are
            1) highcharts (default)
            2) gnuplot

    Returns:
        The plot data is generated and returned in asked format.
    """
    conc_data, en_data, mu_data = dilute_solution_model(
        structure, e0, vac_defs, antisite_defs, T, trial_chem_pot=trial_chem_pot
    )

    if plot_style == "highcharts":
        "Energy data is ignored in this mode"
        hgh_chrt_data = {}
        hgh_chrt_data["xAxis"] = conc_data["x_label"]
        hgh_chrt_data["yAxis"] = conc_data["y_label"]

        series = []
        x = conc_data["x"]
        for y_data in conc_data["y"]:
            y = y_data["data"]
            xy = zip(x, y)
            xy = [list(el) for el in xy]
            name = y_data["name"].strip("$")
            flds = name.split("_")
            def_string = flds[0]
            site_string = flds[1].strip("{}")
            name = def_string + "<sub>" + site_string + "</sub>"
            # series.append({'data':xy, 'name':y_data['name']})
            series.append({"data": xy, "name": name})
        hgh_chrt_data["series"] = series
        return hgh_chrt_data
    if plot_style == "gnuplot":

        def data_to_rows(inp_data):
            rows = []
            labels = []
            labels.append(inp_data["x_label"])
            labels += [y["name"] for y in inp_data["y"]]
            # labels.sort()
            rows.append("#" + "\t".join(labels))
            m = len(inp_data["x"])
            for i in range(m):
                data = []
                data.append(inp_data["x"][i])
                data += [y["data"][i] for y in inp_data["y"]]
                data = [float(x) for x in data]
                rows.append("\t".join(list(map(str, data))))
            return rows

        conc_rows = data_to_rows(conc_data)
        en_rows = data_to_rows(en_data)
        mu_rows = data_to_rows(mu_data)

        return conc_rows, en_rows, mu_rows
    raise ValueError("Invalid plot_style")


# solute_site_preference_finder is based on dilute_solution_model and so most
# of the code is same. However differences exist in setting up and processing
# hence new function
def solute_site_preference_finder(
    structure,
    e0,
    T,
    vac_defs,
    antisite_defs,
    solute_defs,
    solute_concen=0.01,
    trial_chem_pot=None,
):
    """
    Compute the solute defect densities using dilute solution model.
    Args:
        structure: pymatgen.core.structure.Structure object representing the
            primitive or unitcell of the crystal.
        e0: The total energy of the undefected system.
            This is E0 from VASP calculation.
        T: Temperature in Kelvin
        vac_defs: List of vacancy defect parameters in the dictionary format.
            The keys of the dict associated with each vacancy defect are
            1) site_index, 2) site_specie, 3) site_multiplicity, and
            4) energy. 1-3 can be obtained from
            pymatgen.analysis.defects.point_defects.Vacancy class.
            Site index is expected to start with 1 (fortran index).
        antisite_defs: List of antisite defect parameters in the dictionary
            format. The keys of the dict associated with each antisite
            defect are 1) site_index, 2) site_specie, 3) site_multiplicity,
            4) substitution_specie, and 5) energy. 1-3 can be obtained
            from pymatgen.analysis.defects.point_defects.Vacancy class.
        solute_defs: List of solute defect parameters in the dictionary
            format. Similary to that of antisite defs, wtih solute specie
            specified in substitution_specie
        solute_concen: Solute concentration (in fractional value)
        trial_chem_pot: Trial chemical potentials to speedup the plot
            generation. Format is {el1:mu1,...}

    Returns:
        plot_data: The data for plotting the solute defect concentration.
    """

    if not _check_input(vac_defs):
        raise ValueError("Vacancy energy is not defined")
    if not _check_input(antisite_defs):
        raise ValueError("Antisite energy is not defined")

    formation_energies = {}
    formation_energies["vacancies"] = copy.deepcopy(vac_defs)
    formation_energies["antisites"] = copy.deepcopy(antisite_defs)
    formation_energies["solute"] = copy.deepcopy(solute_defs)
    for vac in formation_energies["vacancies"]:
        del vac["energy"]
    for asite in formation_energies["antisites"]:
        del asite["energy"]
    for solute in formation_energies["solute"]:
        del solute["energy"]
    # Setup the system
    site_species = [vac_def["site_specie"] for vac_def in vac_defs]
    solute_specie = solute_defs[0]["substitution_specie"]
    site_species.append(solute_specie)
    multiplicity = [vac_def["site_multiplicity"] for vac_def in vac_defs]
    m = len(set(site_species))  # distinct species
    n = len(vac_defs)  # inequivalent sites

    # Reduce the system and associated parameters such that only distinctive
    # atoms are retained
    comm_div = gcd(*tuple(multiplicity))
    multiplicity = [val / comm_div for val in multiplicity]
    multiplicity.append(0)
    e0 = e0 / comm_div
    T = Float(T)

    # c0 = np.diag(multiplicity)
    c0 = np.diag(np.ones(n + 1))
    c0[n, n] = 0
    mu = [Symbol("mu" + str(i)) for i in range(m)]

    # Generate maps for hashing
    # Generate specie->mu map and use it for site->mu map
    specie_order = []  # Contains hash for site->mu map    Eg: [Al, Ni]
    site_specie_set = set()  # Eg: {Ni, Al}
    for site_specie in site_species:
        if site_specie not in site_specie_set:
            site_specie_set.add(site_specie)
            specie_order.append(site_specie)
    site_mu_map = []  # Eg: [mu0,mu0,mu0,mu1] where mu0->Al, and mu1->Ni
    for site_specie in site_species:
        j = specie_order.index(site_specie)
        site_mu_map.append(j)
    specie_site_index_map = []  # Eg: [(0,3),(3,4)] for Al & Ni
    for i in range(m):
        low_ind = site_species.index(specie_order[i])
        if i < m - 1:
            hgh_ind = site_species.index(specie_order[i + 1])
        else:
            hgh_ind = len(site_species)
        specie_site_index_map.append((low_ind, hgh_ind))
    """
    dC: delta concentration matrix:
    dC[i,j,k]: Concentration change of atom i, due to presence of atom
    j on lattice site k
    Special case is [i,i,i] which is considered as vacancy
    Few cases: dC[i,i,i] = -1 due to being vacancy special case
                dC[k,k,i] = +1 due to increment in k at i lattice if i
                               lattice type is of different element
                dC[i,k,i] = -1 due to decrement of ith type atom due to
                presence of kth type atom on ith sublattice and kth type
                atom specie is different from ith sublattice atom specie
                dC[i,k,k] = 0 due to no effect on ith type atom
                dC[i,j,k] = 0 if i!=j!=k
    """
    dC = np.zeros((n + 1, n + 1, n), dtype=np.int_)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if i == j and site_species[j] != site_species[k] and site_species[i] != site_species:
                    dC[i, j, k] = 1
        for j in range(n + 1):
            for k in range(n):
                if i == k:
                    dC[i, j, k] = -1
    for k in range(n):
        dC[n, n, k] = 1
    for k in range(n):
        for j in range(n):
            if i != j:
                if site_species[i] == site_species[k]:
                    dC[i, j, k] = 0

    for ind_map in specie_site_index_map:
        if ind_map[1] - ind_map[0] > 1:
            for index1 in range(ind_map[0] + 1, ind_map[1]):
                for index2 in range(ind_map[0]):
                    for i in range(n):
                        dC[i, index1, index2] = 0
                for index2 in range(ind_map[1], n):
                    for i in range(n):
                        dC[i, index1, index2] = 0

    # dE matrix: Flip energies (or raw defect energies)
    els = [vac_def["site_specie"] for vac_def in vac_defs]
    dE = []
    for i in range(n + 1):
        dE.append([])
    for i in range(n + 1):
        for j in range(n):
            dE[i].append(0)

    for j in range(n):
        for i in range(n):
            if i == j:
                dE[i][j] = vac_defs[i]["energy"]
            else:
                sub_specie = vac_defs[i]["site_specie"]
                site_specie = vac_defs[j]["site_specie"]
                if site_specie == sub_specie:
                    dE[i][j] = 0
                else:
                    for as_def in antisite_defs:
                        if int(as_def["site_index"]) == j + 1 and sub_specie == as_def["substitution_specie"]:
                            dE[i][j] = as_def["energy"]
                            break
        # Solute
        site_specie = vac_defs[j]["site_specie"]
        for solute_def in solute_defs:
            def_site_ind = int(solute_def["site_index"])
            def_site_specie = solute_def["site_specie"]
            if def_site_specie == site_specie and def_site_ind == j + 1:
                dE[n][j] = solute_def["energy"]
                break

    dE = np.array(dE)
    # np.where(dE == np.array(None), 0, dE)

    # Initialization for concentrations
    # c(i,p) == presence of ith type atom on pth type site
    c = Matrix(n + 1, n, [0] * n * (n + 1))
    for i in range(n + 1):
        for p in range(n):
            c[i, p] = Integer(c0[i, p])
            site_flip_contribs = []
            for epi in range(n + 1):
                sum_mu = sum([mu[site_mu_map[j]] * Integer(dC[j, epi, p]) for j in range(n + 1)])
                flip = dC[i, epi, p] * exp(-(dE[epi, p] - sum_mu) / (k_B * T))
                if flip not in site_flip_contribs:
                    site_flip_contribs.append(flip)
                    c[i, p] += flip
    host_c = Matrix(n, n, [0] * n * n)
    for i in range(n):
        for p in range(n):
            host_c[i, p] = Integer(c0[i, p])
            site_flip_contribs = []
            for epi in range(n):
                sum_mu = sum([mu[site_mu_map[j]] * Integer(dC[j, epi, p]) for j in range(n)])
                flip = dC[i, epi, p] * exp(-(dE[epi, p] - sum_mu) / (k_B * T))
                if flip not in site_flip_contribs:
                    site_flip_contribs.append(flip)
                    host_c[i, p] += flip

    # specie_concen = [sum(mult[ind[0]:ind[1]]) for ind in specie_site_index_map]
    # total_c = [sum(c[ind[0]:ind[1]]) for ind in specie_site_index_map]
    total_c = []
    for ind in specie_site_index_map:
        val = 0
        for i in range(*ind):
            sum_i = sum([c[i, j] * multiplicity[j] for j in range(n)])
            val += sum_i
        total_c.append(val)

    c_ratio = [total_c[i] / sum(total_c) for i in range(m)]

    host_total_c = []
    for ind in specie_site_index_map[:-1]:
        val = 0
        for i in range(*ind):
            sum_i = sum([host_c[i, j] * multiplicity[j] for j in range(n)])
            val += sum_i
        host_total_c.append(val)

    host_c_ratio = [host_total_c[i] / sum(host_total_c) for i in range(m - 1)]

    # Expression for Omega, the Grand Potential
    omega1 = e0 - sum([mu[site_mu_map[i]] * sum(c0[i, :]) * multiplicity[i] for i in range(n)])
    omega = omega1

    used_dEs = []
    for p_r in range(n):
        for epi in range(n):
            sum_mu1 = sum([mu[site_mu_map[j]] * Integer(dC[j, epi, p_r]) for j in range(n)])
            sum_mu = sum_mu1 - mu[site_mu_map[n]] * dC[n, epi, p_r]
            if p_r != epi and site_mu_map[p_r] == site_mu_map[epi]:
                continue
            if dE[epi, p_r] not in used_dEs:
                omega1 -= k_B * T * multiplicity[p_r] * exp(-(dE[epi, p_r] - sum_mu1) / (k_B * T))
                omega -= k_B * T * multiplicity[p_r] * exp(-(dE[epi, p_r] - sum_mu) / (k_B * T))
                used_dEs.append(dE[epi, p_r])

    # Compute composition ranges
    max_host_specie_concen = 1 - solute_concen
    mult = multiplicity
    specie_concen = [sum(mult[ind[0] : ind[1]]) for ind in specie_site_index_map]
    host_specie_concen_ratio = [specie_concen[i] / sum(specie_concen) * max_host_specie_concen for i in range(m)]
    host_specie_concen_ratio[-1] = solute_concen
    li = specie_site_index_map[0][0]
    hi = specie_site_index_map[0][1]
    comp1_min = sum(multiplicity[li:hi]) / sum(multiplicity) * max_host_specie_concen - 0.01
    comp1_max = sum(multiplicity[li:hi]) / sum(multiplicity) * max_host_specie_concen + 0.01
    delta = (comp1_max - comp1_min) / 50.0

    # def reduce_mu():
    #    omega = [e0 - sum([mu[site_mu_map[i]]*sum(c0[i,:]) for i in range(n)])]
    #    x = solve(omega)
    #    return x
    def reduce_mu():
        host_concen = 1 - solute_concen
        new_c0 = c0.astype(float)
        for i in range(n):
            new_c0[i, i] = host_concen * c0[i, i]
        new_c0[n, n] = 2 * solute_concen
        omega = [e0 - sum([mu[site_mu_map[i]] * sum(new_c0[i, :]) for i in range(n + 1)])]
        x = solve(omega)
        return x

    def compute_solute_mu_by_lin_search(host_mu_vals):
        # Compute trial mu
        mult = multiplicity
        specie_concen = [sum(mult[ind[0] : ind[1]]) for ind in specie_site_index_map]
        max_host_specie_concen = 1 - solute_concen
        host_specie_concen_ratio = [specie_concen[i] / sum(specie_concen) * max_host_specie_concen for i in range(m)]
        host_specie_concen_ratio[-1] = solute_concen
        y_vect = host_specie_concen_ratio
        vector_func = [y_vect[i] - c_ratio[i] for i in range(m)]
        vector_func.append(omega)
        mu_vals = None
        m1_min = -20.0
        if e0 > 0:
            m1_max = 10  # Search space needs to be modified
        else:
            m1_max = 0
        for m1 in np.arange(m1_min, m1_max, 0.1):
            trial_mus = host_mu_vals + [m1]
            try:
                x = nsolve(vector_func, mu, trial_mus, module="numpy")
                if x:
                    mu_vals = [float(mu_val) for mu_val in x]
                break
            except Exception:
                continue
        else:
            raise ValueError()
        return mu_vals

    def compute_mus():

        # Compute trial mu
        mu_red = reduce_mu()
        mult = multiplicity
        specie_concen = [sum(mult[ind[0] : ind[1]]) for ind in specie_site_index_map]
        max_host_specie_concen = 1 - solute_concen
        host_specie_concen_ratio = [specie_concen[i] / sum(specie_concen) * max_host_specie_concen for i in range(m)]
        host_specie_concen_ratio[-1] = solute_concen

        y_vect = host_specie_concen_ratio
        vector_func = [y_vect[i] - c_ratio[i] for i in range(m)]
        vector_func.append(omega)
        mu_vals = None
        m_min = -15.0
        if e0 > 0:
            m_max = 10  # Search space needs to be modified
        else:
            m_max = 0
        for m1 in np.arange(m_min, m_max, 0.3):
            for m2 in np.arange(m_min, m_max, 0.3):
                m0 = mu_red[mu[0]].subs([(mu[1], m1), (mu[2], m2)])
                try:
                    mu_vals = nsolve(vector_func, mu, [m0, m1, m2], module="numpy")
                    # Line needs to be modified to include all mus when n > 2
                except Exception:
                    continue
                break
            if mu_vals:
                mu_vals = [float(mu_val) for mu_val in mu_vals]
                break
        else:
            raise ValueError("Couldn't find mus")
        return mu_vals

    if not trial_chem_pot:
        # Try computing mus by assuming one of the defects is dominant at 0.01
        # concen.  First vacancy is tried and then antisite

        # Generate trial mus assuming vacancy as dominant defect
        # for specie-0 at lower yval
        li = specie_site_index_map[0][0]
        hi = specie_site_index_map[0][1]
        li1 = specie_site_index_map[1][0]
        hi1 = specie_site_index_map[1][1]
        spec_mult = [sum(multiplicity[li:hi]), sum(multiplicity[li1:hi1])]
        ln_def_conc = 4.60517
        for i in range(li, hi):
            vac_flip_en = vac_defs[i]["energy"]
            mu_vals = [ln_def_conc * k_B * T - vac_flip_en]
            mu_vals.append((e0 - spec_mult[0] * mu_vals[0]) / spec_mult[1])
            comp_ratio = comp1_min

            # Test if the trial mus are good
            vector_func = [comp_ratio - host_c_ratio[0]]
            vector_func.append(omega1)
            try:
                host_mu_vals = nsolve(vector_func, mu[:-1], mu_vals)
                if host_mu_vals:
                    host_mu_vals = [float(mu_val) for mu_val in host_mu_vals]
                compute_solute_mu_by_lin_search(host_mu_vals)
                break
            except Exception:  # Go for antisite as dominant defect
                mu_gs = [Symbol("mu_gs" + j.__str__()) for j in range(m - 1)]

                eqs = [mu_gs[0] - mu_gs[1] - (ln_def_conc * k_B * T - antisite_defs[i]["energy"])]
                eqs.append(spec_mult[0] * mu_gs[0] + spec_mult[1] * mu_gs[1] - e0)
                x = solve(eqs, mu_gs)
                host_mu_vals = []
                for key in sorted(x.keys(), key=lambda inp: inp.name):
                    host_mu_vals.append(x[key])
                vector_func = [comp_ratio - host_c_ratio[0]]
                vector_func.append(omega1)

                try:
                    host_mu_vals = nsolve(vector_func, mu[:-1], host_mu_vals)
                    if host_mu_vals:
                        host_mu_vals = [float(mu_val) for mu_val in host_mu_vals]
                    mu_vals = compute_solute_mu_by_lin_search(host_mu_vals)
                    break
                except Exception:  # Go to the default option (search the space)
                    pass
        else:
            mu_vals = compute_mus()

    else:
        try:
            mu_vals = [trial_chem_pot[element] for element in specie_order]
        except Exception:
            mu_vals = compute_mus()

    # Compile mu's for all composition ratios in the range
    # +/- 1% from the stoichiometry
    result = {}
    for y in np.arange(comp1_min, comp1_max + delta, delta):
        y_vect = []
        y_vect.append(y)
        y2 = max_host_specie_concen - y
        y_vect.append(y2)
        y_vect.append(solute_concen)
        vector_func = [y_vect[i] - c_ratio[i] for i in range(1, m)]
        vector_func.append(omega)

        try:
            x = nsolve(vector_func, mu, mu_vals)
            if x:
                mu_vals = [float(mu_val) for mu_val in x]
        except Exception:
            continue
        result[y] = mu_vals

    res = []

    # Compute the concentrations for all the compositions
    for key in sorted(result.keys()):
        mu_val = result[key]
        total_c_val = [total_c[i].subs(dict(zip(mu, mu_val))) for i in range(len(total_c))]
        c_val = c.subs(dict(zip(mu, mu_val)))
        # Concentration of first element/over total concen
        res1 = []
        res1.append(float(total_c_val[0] / sum(total_c_val)))

        for i in range(n + 1):
            for j in range(n):
                if i == j:  # Vacancy
                    vac_conc = float(exp(-(mu_val[site_mu_map[i]] + dE[i, i]) / (k_B * T)))
                    res1.append(vac_conc)
                else:  # Antisite
                    res1.append(float(c_val[i, j] / c0[j, j]))
        res.append(res1)

    res = np.array(res)
    dtype = [(str("x"), np.float64)] + [(str("y%d%d" % (i, j)), np.float64) for i in range(n + 1) for j in range(n)]
    res1 = np.sort(res.view(dtype), order=[str("x")], axis=0)

    conc = []
    for i in range(n + 1):
        conc.append([])
        for j in range(n):
            conc[i].append([])
    for i in range(n + 1):  # Append vacancies
        for j in range(n):
            y1 = [dat[0][i * n + j + 1] for dat in res1]
            conc[i][j] = y1

    # Compute solute site preference
    # Removing the functionality
    # site_pref_data = {}
    """Because all the plots have identical x-points storing it in a
    single array"""
    # site_pref_data['x'] = [dat[0][0] for dat in res1]         # x-axis data
    # Element whose composition is varied. For x-label
    # site_pref_data['x_label'] = els[0]+ "_mole_fraction"
    # site_pref_data['y_label'] = "$"+solute_specie+"_{"+els[0]+"}/("+\
    #    solute_specie+"_{"+els[0]+"}+"+solute_specie+"_{"+els[1]+"})$"

    # y_data = []
    # inds = specie_site_index_map[m-1]
    # data1 = np.sum([multiplicity[0]*conc[ind][0] for ind in range(*inds)],axis=0)
    # data2 = np.sum([multiplicity[1]*conc[ind][1] for ind in range(*inds)],axis=0)
    # frac_data = data1/(data1+data2)
    # frac_data = frac_data.tolist()
    # y_data.append({'data':frac_data})

    # site_pref_data['y'] = y_data

    #  Return all defect concentrations
    conc_data = {}
    """Because all the plots have identical x-points storing it in a
    single array"""
    conc_data["x"] = [dat[0][0] for dat in res1]  # x-axis data
    # Element whose composition is varied. For x-label
    conc_data["x_label"] = els[0] + " mole fraction"
    conc_data["y_label"] = "Point defect concentration"

    y_data = []
    # Vacancy
    for i in range(n):
        data = conc[i][i]
        specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        vac_string = "$Vac_{"
        if not specie_ind_del - 1:
            label = vac_string + specie + "}$"
        else:
            label = vac_string + specie + "_" + str(cur_ind) + "}$"
        # Plot data and legend info
        y_data.append({"data": data, "name": label})

    # Antisites and solute
    for i in range(n):
        site_specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        for j in range(m):
            sub_specie = specie_order[j]
            if sub_specie == site_specie:
                continue
            if not specie_ind_del - 1:
                label = "$" + sub_specie + "_{" + site_specie + "}$"
            else:
                label = "$" + sub_specie + "_{" + site_specie + "_" + str(cur_ind) + "}$"
            inds = specie_site_index_map[j]
            # TODO: Investigate the value below
            data = np.sum([conc[ind][i] for ind in range(*inds)], axis=0)
            data = data.tolist()
            y_data.append({"data": data, "name": label})

    conc_data["y"] = y_data
    # return site_pref_data, conc_data
    return conc_data


def solute_defect_density(
    structure,
    e0,
    vac_defs,
    antisite_defs,
    solute_defs,
    solute_concen=0.01,
    T=800,
    trial_chem_pot=None,
    plot_style="highchargs",
):
    """
    Wrapper for the solute_site_preference_finder.
    The computed plot data is prepared based on plot_style.

    Args:
        structure: pymatgen.core.structure.Structure object representing the
            primitive or unitcell of the crystal.
        e0: The total energy of the undefected system.
            This is E0 from VASP calculation.
        vac_defs: List of vacancy defect parameters in the dictionary format.
            The keys of the dict associated with each vacancy defect are
            1) site_index, 2) site_specie, 3) site_multiplicity, and
            4) energy. 1-3 can be obtained from
            pymatgen.analysis.defects.point_defects.Vacancy class.
            Site index is expected to start with 1 (fortran index).
        antisite_defs: List of antisite defect parameters in the dictionary
            format. The keys of the dict associated with each antisite defect
            are 1) site_index, 2) site_specie, 3) site_multiplicity,
            4) substitution_specie, and 5) energy. 1-3 can be obtained
            from pymatgen.analysis.defects.point_defects.Vacancy class.
        solute_defs: List of solute defect parameters in the dictionary
            format. Similary to that of antisite defs, wtih solute specie
            specified in substitution_specie
        solute_concen: Solute concentration (in fractional value)
        T: Temperature in Kelvin
        trial_chem_pot (optional): Trial chemical potentials to speedup
            the plot generation. Format is {el1:mu1,...}
        plot_style (string): Allowed options are
            1) highcharts (default)
            2) gnuplot

    Returns:
        The plot data is generated and returned in asked format.
    """

    # solute_site_pref_data, def_conc_data = solute_site_preference_finder(
    def_conc_data = solute_site_preference_finder(
        structure,
        e0,
        T,
        vac_defs,
        antisite_defs,
        solute_defs,
        solute_concen=solute_concen,
        trial_chem_pot=trial_chem_pot,
    )

    if plot_style == "highcharts":
        "Energy data is ignored in this mode"
        hgh_chrt_data = {}
        hgh_chrt_data["xAxis"] = def_conc_data["x_label"]
        hgh_chrt_data["yAxis"] = def_conc_data["y_label"]

        series = []
        x = def_conc_data["x"]
        for y_data in def_conc_data["y"]:
            y = y_data["data"]
            xy = zip(x, y)
            xy = [list(el) for el in xy]
            name = y_data["name"].strip("$")
            flds = name.split("_")
            def_string = flds[0]
            site_string = flds[1].strip("{}")
            name = def_string + "<sub>" + site_string + "</sub>"
            # series.append({'data':xy, 'name':y_data['name']})
            series.append({"data": xy, "name": name})
        hgh_chrt_data["series"] = series
        return hgh_chrt_data
    if plot_style == "gnuplot":

        def data_to_rows(inp_data, y_lbl_flg):
            rows = []
            labels = []
            labels.append(inp_data["x_label"])
            if y_lbl_flg:
                labels.append(inp_data["y_label"])
            else:
                labels += [y["name"] for y in inp_data["y"]]
            rows.append("#" + "\t".join(labels))
            m = len(inp_data["x"])
            for i in range(m):
                data = []
                data.append(inp_data["x"][i])
                data += [y["data"][i] for y in inp_data["y"]]
                data = [float(x) for x in data]
                rows.append("\t".join(list(map(str, data))))
            return rows

        # solute_site_pref_rows = data_to_rows(solute_site_pref_data, True)
        pt_def_conc_rows = data_to_rows(def_conc_data, False)
        # return solute_site_pref_rows, pt_def_conc_rows
        return pt_def_conc_rows

    raise ValueError("Invalid plot_style.")
