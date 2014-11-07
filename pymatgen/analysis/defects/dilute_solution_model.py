# coding: utf-8

from __future__ import unicode_literals
from __future__ import division

"""
Evaluate the defect concentration based on composition, temperature,
and defect energies using "Dilute Solution Model"
Reference: Phys Rev B, 63, 094103, 2001,
"Density of constitutional and thermal point defects in L12 Al3Sc",
C. Woodward, M. Asta, G. Kresse and J. Hafner.
"""


__author__ = 'Bharat Medasani'
__version__ = "0.2"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Alpha"
__date__ = "6/4/14"

import math
import copy
import numpy as np
from six.moves import zip

from monty.dev import requires
from monty.fractions import gcd


try:
    from sympy import Symbol, nsolve, Integer, Float, Matrix, exp, solve, Eq
    sympy_found = True
except ImportError:
    sympy_found = False

# physical consts
k_B=8.6173324e-5                # eV/K

# Check the inputs
def check_input(def_list):
    flag = True
    for defect in def_list:
        if not defect:
            flag = False
            break
    return flag

# TODO: Cleanup the initialization

@requires(sympy_found,
            "dilute_solution_model requires Sympy module. Please install it.")
def dilute_solution_model(structure, e0, vac_defs, antisite_defs, T,
        trial_chem_pot = None, generate='plot'):

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

    if not check_input(vac_defs):
        raise ValueError('Vacancy energy is not defined')
    if not check_input(antisite_defs):
        raise ValueError('Antisite energy is not defined')

    formation_energies = {}
    formation_energies['vacancies'] = copy.deepcopy(vac_defs)
    formation_energies['antisites'] = copy.deepcopy(antisite_defs)
    for vac in formation_energies['vacancies']:
        del vac['energy']
    for asite in formation_energies['antisites']:
        del asite['energy']
    # Setup the system
    site_species = [vac_def['site_specie'] for vac_def in vac_defs]
    multiplicity = [vac_def['site_multiplicity'] for vac_def in vac_defs]
    m = len(set(site_species))      # distinct species
    n = len(vac_defs)           # inequivalent sites

    # Reduce the system and associated parameters such that only distinctive
    # atoms are retained
    comm_div = gcd(*tuple(multiplicity))
    multiplicity = [val/comm_div for val in multiplicity]
    e0 = e0/comm_div
    T = Integer(T)

    c0 = np.diag(multiplicity)
    #print ('c0',c0)
    mu = [Symbol('mu'+i.__str__()) for i in range(m)]

    # Generate maps for hashing
    # Generate specie->mu map and use it for site->mu map
    specie_order = []       # Contains hash for site->mu map    Eg: [Al, Ni]
    site_specie_set = set()             # Eg: {Ni, Al}
    for i in range(n):
        site_specie  = site_species[i]
        if site_specie not in site_specie_set:
            site_specie_set.add(site_specie)
            specie_order.append(site_specie)
    site_mu_map = []     # Eg: [mu0,mu0,mu0,mu1] where mu0->Al, and mu1->Ni
    for i in range(n):
        site_specie  = site_species[i]
        j = specie_order.index(site_specie)
        site_mu_map.append(j)
    specie_site_index_map = []      # Eg: [(0,3),(3,4)] for Al & Ni
    for i in range(m):
        low_ind = site_species.index(specie_order[i])
        if i < m-1:
            hgh_ind = site_species.index(specie_order[i+1])
        else:
            hgh_ind = n
        specie_site_index_map.append((low_ind,hgh_ind))


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
    dC = np.zeros((n,n,n), dtype=np.int)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if i == j and site_species[j] != site_species[k] and \
                                site_species[i] != site_species[k]:
                    dC[i,j,k] = 1
        for j in range(n):
            for k in range(n):
                if i == k:
                    #if j == k or site_species[j] != site_species[k]:
                        dC[i,j,k] = -1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                if i != j:
                    if site_species[j] == site_species[k]:
                        dC[i,j,k] = 0
        #for j in range(n):
        #    if k != j:
        #        for i in range(n):
        #            if i == j:
        #                if abs(i-j) <= 1 and abs(j-k) <= 1 and abs(i-k) <= 1:
        #                    dC[i,j,k] = 0

    for ind_map in specie_site_index_map:
        if ind_map[1]-ind_map[0] > 1:
            for index1 in range(ind_map[0]+1,ind_map[1]):
                for index2 in range(ind_map[0]):
                    for i in range(n):
                        #print (i, index1, index2)
                        dC[i,index1,index2] = 0
                for index2 in range(ind_map[1],n):
                    for i in range(n):
                        #print (i, index1, index2)
                        dC[i,index1,index2] = 0


    #print ('dC', dC)
    # dE matrix: Flip energies (or raw defect energies)
    els = [vac_def['site_specie'] for vac_def in vac_defs]
    dE = []
    for i in range(n):
        dE.append([])
    for i in range(n):
        for j in range(n):
            dE[i].append(0)

    for j in range(n):
        for i in range(n):
            if i == j:
                dE[i][j] = vac_defs[i]['energy']
            else:
                sub_specie = vac_defs[i]['site_specie']
                site_specie = vac_defs[j]['site_specie']
                if site_specie == sub_specie:
                    dE[i][j] = 0
                else:
                    for as_def in antisite_defs:
                        if int(as_def['site_index']) == j+1 and \
                                sub_specie == as_def['substitution_specie']:
                            dE[i][j] = as_def['energy']
                            break
    dE = np.array(dE)
    #np.where(dE is None, dE, 0)
    #print ('dE', dE)

    # Initialization for concentrations
    # c(i,p) == presence of ith type atom on pth type site
    c = Matrix(n,n,[0]*n**2)
    for i in range(n):
        for p in range(n):
            c[i,p] = Integer(c0[i,p])
            site_flip_contribs = []
            for epi in range(n):
                sum_mu = sum([mu[site_mu_map[j]]*Integer(dC[j,epi,p]) \
                        for j in range(n)])
                #print (i, epi, p, dC[i,epi, p], sum_mu)
                #c[i,p] += Integer(multiplicity[p]*dC[i,epi,p]) * \
                #        exp(-(dE[epi,p]-sum_mu)/(k_B*T))
                flip = Integer(multiplicity[p]*dC[i,epi,p]) * \
                        exp(-(dE[epi,p]-sum_mu)/(k_B*T))
                if flip not in site_flip_contribs:
                    site_flip_contribs.append(flip)
                    c[i,p] += flip
                #else:
                    #print (i, epi, p)
                    #print ('flip already present in site_flips')

    #print ('c', c)
    total_c = []
    for ind in specie_site_index_map:
        total_c.append(sum([sum(c[i,:]) for i in range(*ind)]))
    c_ratio = [total_c[-1]/total_c[i] for i in range(m)]
    #print ('c_ratio')
    #for i in range(len(c_ratio)):
    #    print(c_ratio[i])

    # Expression for Omega, the Grand Potential
    omega = e0 - sum([mu[site_mu_map[i]]*sum(c0[i,:]) for i in range(n)])
    for p_r in range(n):
        for epi in range(n):
            sum_mu = sum([mu[site_mu_map[j]]*Float(
                    dC[j,epi,p_r]) for j in range(n)])
            omega -= k_B*T*multiplicity[p_r]*exp(-(dE[epi,p_r]-sum_mu)/(k_B*T))

    def compute_mus():

        def reduce_mu():
            omega = [e0 - sum([mu[site_mu_map[i]]*sum(c0[i,:]) for i in range(n)])]
            x = solve(omega)
            return x

        # Compute trial mu
        mu_red = reduce_mu()

        mult = multiplicity
        specie_concen = [sum(mult[ind[0]:ind[1]]) for ind in specie_site_index_map]
        y_vect = [specie_concen[-1]/specie_concen[i] for i in range(m)]
        vector_func = [y_vect[i]-c_ratio[i] for i in range(m-1)]
        vector_func.append(omega)
        min_diff = 1e10
        mu_vals = None
        c_val = None
        m1_min = -20.0
        if e0 > 0:
            m1_max = 10            # Search space needs to be modified
        else:
            m1_max = 0
        for m1 in np.arange(m1_min,m1_max,0.1):
            m0 = mu_red[mu[0]].subs(mu[-1],m1)

            try:
                x = nsolve(vector_func,mu,[m0,m1],module="numpy")
            except:
                continue

            c_val = c.subs(dict(zip(mu,x)))
            #if all(x >= 0 for x in c_val):
            specie_concen = []
            for ind in specie_site_index_map:
                specie_concen.append(sum([sum(c_val[i,:]) for i in range(*ind)]))
            y_comp = [specie_concen[-1]/specie_concen[i] for i in range(m)]
            diff = math.sqrt(sum([pow(abs(y_comp[i]-y_vect[i]),2) for i in range(m)]))
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
            site_specie = vac_def['site_specie']
            ind = specie_order.index(site_specie)
            uncor_energy = vac_def['energy']
            formation_energy = uncor_energy + mu_vals[ind]
            formation_energies['vacancies'][i]['formation_energy'] = formation_energy
            specie_ind = site_mu_map[i]
            indices = specie_site_index_map[specie_ind]
            specie_ind_del = indices[1]-indices[0]
            cur_ind = i - indices[0] + 1
            if not specie_ind_del-1:
                label = '$V_{'+site_specie+'}$'
            else:
                label = '$V_{'+site_specie+'_'+str(cur_ind)+'}$'
            formation_energies['vacancies'][i]['label'] = label
            i += 1
        i = 0
        for as_def in antisite_defs:
            site_specie = as_def['site_specie']
            sub_specie = as_def['substitution_specie']
            ind1 = specie_order.index(site_specie)
            ind2 = specie_order.index(sub_specie)
            uncor_energy = as_def['energy']
            formation_energy = uncor_energy + mu_vals[ind1] - mu_vals[ind2]
            formation_energies['antisites'][i]['formation_energy'] = formation_energy
            specie_ind = site_mu_map[i]
            indices = specie_site_index_map[specie_ind]
            specie_ind_del = indices[1]-indices[0]
            cur_ind = i - indices[0] + 1
            if not specie_ind_del-1:
                label = '$'+sub_specie+'_{'+site_specie+'}$'
            else:
                label = '$'+sub_specie+'_{'+site_specie+'_'+str(cur_ind)+'}$'
            formation_energies['antisites'][i]['label'] = label
            i += 1
        return formation_energies

    if not trial_chem_pot:
        mu_vals = compute_mus()
    else:
        try:
            mu_vals = [trial_chem_pot[element] for element in specie_order]
        except:
            mu_vals = compute_mus()


    if generate == 'energy':
        formation_energies = compute_def_formation_energies()
        mu_dict = dict(zip(specie_order,mu_vals))
        return formation_energies, mu_dict

    #sys.exit()

    # Compute ymax
    li = specie_site_index_map[0][0]
    hi = specie_site_index_map[0][1]
    comp1_min = int(sum(multiplicity[li:hi])/sum(multiplicity)*100)-1
    comp1_max = int(sum(multiplicity[li:hi])/sum(multiplicity)*100)+1
    delta = float(comp1_max-comp1_min)/120.0
    yvals = []
    for comp1 in np.arange(comp1_min,comp1_max+delta,delta):
        comp2 = 100-comp1
        y = comp2/comp1
        yvals.append(y)
    #ymin = comp2_min/comp1_max
    #print comp1_min, comp1_max
    #print comp2_min, comp2_max
    #print ymax, ymin

    # Compile mu's for all composition ratios in the range
    #+/- 1% from the stoichiometry
    result = {}
    #for y in np.arange(ymin,ymax+delta,delta):
    for y in yvals:
        result[y] = []
        vector_func = [y-c_ratio[0]]
        vector_func.append(omega)
        x = nsolve(vector_func,mu,mu_vals,module="numpy")
        result[y].append(x[0])
        result[y].append(x[1])

    res = []
    new_mu_dict = {}
    # Compute the concentrations for all the compositions
    for key in sorted(result.keys()):
        mu_val = result[key]
        total_c_val = [total_c[i].subs(dict(zip(mu,mu_val))) \
                for i in range(len(total_c))]
        c_val = c.subs(dict(zip(mu,mu_val)))
        res1 = []
        # Concentration of first element/over total concen
        res1.append(float(total_c_val[0]/sum(total_c_val)))
        new_mu_dict[res1[0]] = mu_val
        sum_c0 = sum([c0[i,i] for i in range(n)])
        #print res1[0]
        for i in range(n):
            for j in range(n):
                if i == j:              # Vacancy
                    #print (i,  c_val[:,i])
                    # Consider numerical accuracy
                    #res1.append(float((c0[i,i]-sum(c_val[:,i]))/c0[i,i]))
                    #print ((mu_val[site_mu_map[i]]-dE[i,i])/(k_B*T))
                    vac_conc = float(exp(-(mu_val[site_mu_map[i]]+dE[i,i])/(k_B*T)))
                    #print vac_conc
                    res1.append(vac_conc)
                else:                   # Antisite
                    res1.append(float(c_val[i,j]/c0[j,j]))
        res.append(res1)

    res = np.array(res)
    dtype = [(str('x'),np.float64)]+[(str('y%d%d' % (i, j)), np.float64) \
            for i in range(n) for j in range(n)]
    res1 = np.sort(res.view(dtype), order=[str('x')],axis=0)


    conc_data = {}
    """Because all the plots have identical x-points storing it in a
    single array"""
    conc_data['x'] = [dat[0][0] for dat in res1]         # x-axis data
    # Element whose composition is varied. For x-label
    conc_data['x_label'] = els[0]+ " mole fraction"
    conc_data['y_label'] = "Point defect concentration"
    conc = []
    for i in range(n):
        conc.append([])
        for j in range(n):
            conc[i].append([])
    for i in range(n): 
        for j in range(n):
            y1 = [dat[0][i*n+j+1] for dat in res1]
            conc[i][j] = y1

    y_data = []
    for i in range(n):      
        data = conc[i][i]
        specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1]-indices[0]
        cur_ind = i - indices[0] + 1
        #if 'V' not in els:
        #    vac_string = "$V_{"
        #else:
        vac_string = "$Vac_{"
        if not specie_ind_del-1:
            label = vac_string+specie+'}$'
        else:
            label = vac_string+specie+'_'+str(cur_ind)+'}$'
        # Plot data and legend info
        y_data.append({'data':data,'name':label})

    for i in range(n):      
        site_specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1]-indices[0]
        cur_ind = i - indices[0] + 1
        for j in range(m):          # Antisite plot dat
            sub_specie = specie_order[j]
            if sub_specie == site_specie:
                continue
            if not specie_ind_del-1:
                label = '$'+sub_specie+'_{'+site_specie+'}$'
            else:
                label = '$'+sub_specie+'_{'+site_specie+'_'+str(cur_ind)+'}$'
            inds = specie_site_index_map[j]
            data = np.sum([conc[ind][i] for ind in range(*inds)],axis=0)
            data = data.tolist()
            y_data.append({'data':data,'name':label})

    conc_data['y'] = y_data

    # Compute the  formation energies
    def compute_vac_formation_energies(mu_vals):
        en = []
        for vac_def in vac_defs:
            site_specie = vac_def['site_specie']
            ind = specie_order.index(site_specie)
            uncor_energy = vac_def['energy']
            formation_energy = uncor_energy + mu_vals[ind]
            en.append(float(formation_energy))
        return en
    en_res = []
    for key in sorted(new_mu_dict.keys()):
        mu_val = new_mu_dict[key]
        en_res.append(compute_vac_formation_energies(mu_val))

    en_data = {'x_label':els[0]+' mole fraction', 'x':[]}
    en_data['x'] = [dat[0][0] for dat in res1]         # x-axis data

    i = 0
    y_data = []
    for vac_def in vac_defs:
        data = [data[i] for data in en_res]
        site_specie = vac_def['site_specie']
        ind = specie_order.index(site_specie)
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1]-indices[0]
        cur_ind = i - indices[0] + 1
        vac_string = "$Vac_{"
        if not specie_ind_del-1:
            label = vac_string+site_specie+'}$'
        else:
            label = vac_string+site_specie+'_'+str(cur_ind)+'}$'
        y_data.append({'data':data,'name':label})
        i += 1

    def compute_as_formation_energies(mu_vals):
        en = []
        for as_def in antisite_defs:
            site_specie = as_def['site_specie']
            sub_specie = as_def['substitution_specie']
            ind1 = specie_order.index(site_specie)
            ind2 = specie_order.index(sub_specie)
            uncor_energy = as_def['energy']
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
        site_specie = as_def['site_specie']
        sub_specie = as_def['substitution_specie']
        ind1 = specie_order.index(site_specie)
        ind2 = specie_order.index(sub_specie)
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1]-indices[0]
        cur_ind = i - indices[0] + 1
        if not specie_ind_del-1:
            label = '$'+sub_specie+'_{'+site_specie+'}$'
        else:
            label = '$'+sub_specie+'_{'+site_specie+'_'+str(cur_ind)+'}$'
        y_data.append({'data':data,'name':label})
        i += 1

    en_data['y'] = y_data

    # Return chem potential as well
    mu_data = {'x_label':els[0]+' mole fraction', 'x':[]}
    mu_data['x'] = [dat[0][0] for dat in res1]         # x-axis data

    y_data = []
    for j in range(m):
        specie = specie_order[j]
        mus = [new_mu_dict[key][j] for key in sorted(new_mu_dict.keys())]
        y_data.append({'data':mus, 'name':specie})
    mu_data['y'] = y_data

    return conc_data, en_data, mu_data


@requires(sympy_found,
          "comute_defect_density requires Sympy module. Please install it.")
def compute_defect_density(structure, e0, vac_defs, antisite_defs, T=800,
        trial_chem_pot=None, plot_style="highcharts"):
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
            structure,e0,vac_defs,antisite_defs,T,
            trial_chem_pot=trial_chem_pot)

    if plot_style == 'highcharts':
        "Energy data is ignored in this mode"
        hgh_chrt_data = {}
        hgh_chrt_data['xAxis'] = conc_data['x_label']
        hgh_chrt_data['yAxis'] = conc_data['y_label']

        series = []
        x = conc_data['x']
        for y_data in conc_data['y']:
            y = y_data['data']
            xy = zip(x,y)
            xy = [list(el) for el in xy]
            name = y_data['name'].strip('$')
            flds= name.split('_')
            def_string = flds[0]
            site_string = flds[1].strip('{}')
            name = def_string+"<sub>"+site_string+"</sub>"
            #series.append({'data':xy, 'name':y_data['name']})
            series.append({'data':xy, 'name':name})
        hgh_chrt_data['series'] = series
        return hgh_chrt_data
    elif plot_style == 'gnuplot':
        def data_to_rows(inp_data):
            rows = []
            labels = []
            labels.append(inp_data['x_label'])
            labels += [y['name'] for y in inp_data['y']]
            #labels.sort()
            rows.append('#'+'\t'.join(labels))
            m = len(inp_data['x'])
            for i in range(m):
                data = []
                data.append(inp_data['x'][i])
                data += [y['data'][i] for y in inp_data['y']]
                data = [float(x) for x in data]
                rows.append('\t'.join(list(map(str,data))))
            return rows
        conc_rows = data_to_rows(conc_data)
        en_rows = data_to_rows(en_data)
        mu_rows = data_to_rows(mu_data)

        return conc_rows, en_rows, mu_rows


#solute_site_preference_finder is based on dilute_solution_model and so most
#of the code is same. However differences exist in setting up and processing
#hence new function
@requires(sympy_found, "solute_site_preference_finder requires Sympy module. "\
        "Please install it.")
def solute_site_preference_finder(
        structure, e0, T, vac_defs, antisite_defs,  solute_defs,
        solute_concen=0.01, trial_chem_pot = None):

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

    if not check_input(vac_defs):
        raise ValueError('Vacancy energy is not defined')
    if not check_input(antisite_defs):
        raise ValueError('Antisite energy is not defined')

    formation_energies = {}
    formation_energies['vacancies'] = copy.deepcopy(vac_defs)
    formation_energies['antisites'] = copy.deepcopy(antisite_defs)
    formation_energies['solute'] = copy.deepcopy(solute_defs)
    for vac in formation_energies['vacancies']:
        del vac['energy']
    for asite in formation_energies['antisites']:
        del asite['energy']
    for solute in formation_energies['solute']:
        del solute['energy']
    # Setup the system
    site_species = [vac_def['site_specie'] for vac_def in vac_defs]
    solute_specie = solute_defs[0]['substitution_specie']
    site_species.append(solute_specie)
    multiplicity = [vac_def['site_multiplicity'] for vac_def in vac_defs]
    m = len(set(site_species))      # distinct species
    n = len(vac_defs)           # inequivalent sites

    # Reduce the system and associated parameters such that only distinctive
    # atoms are retained
    comm_div = gcd(*tuple(multiplicity))
    multiplicity = [val/comm_div for val in multiplicity]
    multiplicity.append(0)
    e0 = e0/comm_div
    T = Integer(T)

    c0 = np.diag(multiplicity)
    #print(('c0', c0))
    mu = [Symbol('mu'+str(i)) for i in range(m)]

    # Generate maps for hashing
    # Generate specie->mu map and use it for site->mu map
    specie_order = []       # Contains hash for site->mu map    Eg: [Al, Ni]
    site_specie_set = set()             # Eg: {Ni, Al}
    for i in range(len(site_species)):
        site_specie  = site_species[i]
        if site_specie not in site_specie_set:
            site_specie_set.add(site_specie)
            specie_order.append(site_specie)
    site_mu_map = []     # Eg: [mu0,mu0,mu0,mu1] where mu0->Al, and mu1->Ni
    for i in range(len(site_species)):
        site_specie  = site_species[i]
        j = specie_order.index(site_specie)
        site_mu_map.append(j)
    specie_site_index_map = []      # Eg: [(0,3),(3,4)] for Al & Ni
    for i in range(m):
        low_ind = site_species.index(specie_order[i])
        if i < m-1:
            hgh_ind = site_species.index(specie_order[i+1])
        else:
            hgh_ind = len(site_species)
        specie_site_index_map.append((low_ind,hgh_ind))

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
    dC = np.zeros((n+1,n+1,n), dtype=np.int)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if i == j and site_species[j] != site_species[k] and \
                        site_species[i] != site_species:
                    dC[i,j,k] = 1
        for j in range(n+1):
            for k in range(n):
                if i == k:
                    #if j == k or site_species[j] != site_species[k]:
                        dC[i,j,k] = -1
    for k in range(n):
        dC[n,n,k] = 1
    for k in range(n):
        for j in range(n):
            if i != j:
                if site_species[i] == site_species[k]:
                    dC[i,j,k] = 0

    for ind_map in specie_site_index_map:
        if ind_map[1]-ind_map[0] > 1:
            for index1 in range(ind_map[0]+1,ind_map[1]):
                for index2 in range(ind_map[0]):
                    for i in range(n):
                        print (i, index1, index2)
                        dC[i,index1,index2] = 0
                for index2 in range(ind_map[1],n):
                    for i in range(n):
                        print (i, index1, index2)
                        dC[i,index1,index2] = 0

    print ('dC', dC)
    # dE matrix: Flip energies (or raw defect energies)
    els = [vac_def['site_specie'] for vac_def in vac_defs]
    dE = []
    for i in range(n+1):
        dE.append([])
    for i in range(n+1):
        for j in range(n):
            dE[i].append(0)

    for j in range(n):
        for i in range(n):
            if i == j:
                dE[i][j] = vac_defs[i]['energy']
            else:
                sub_specie = vac_defs[i]['site_specie']
                site_specie = vac_defs[j]['site_specie']
                if site_specie == sub_specie:
                    dE[i][j] = 0
                else:
                    for as_def in antisite_defs:
                        if int(as_def['site_index']) == j+1 and \
                                sub_specie == as_def['substitution_specie']:
                            dE[i][j] = as_def['energy']
                            break
        # Solute
        site_specie = vac_defs[j]['site_specie']
        for solute_def in solute_defs:
            def_site_ind = int(solute_def['site_index'])
            def_site_specie = solute_def['site_specie']
            #print((def_site_specie, site_specie))
            #print((def_site_ind, j+1))
            if def_site_specie == site_specie and def_site_ind == j+1:
                #print(('se', solute_def['energy']))
                dE[n][j] = solute_def['energy']
                break

    dE = np.array(dE)
    #np.where(dE == np.array(None), 0, dE)

    # Initialization for concentrations
    # c(i,p) == presence of ith type atom on pth type site
    c = Matrix(n+1,n,[0]*n*(n+1))
    for i in range(n+1):
        for p in range(n):
            c[i,p] = Integer(c0[i,p])
            #print((c[i,p]))
            for epi in range(n+1):
                sum_mu = sum([mu[site_mu_map[j]]*Integer(
                        dC[j,epi,p]) for j in range(n+1)])
                #print((sum_mu))
                #print((multiplicity[p], dC[i,epi,p], dE[epi,p]))
                c[i,p] += Integer(multiplicity[p]*dC[i,epi,p]) * \
                        exp(-(dE[epi,p]-sum_mu)/(k_B*T))
    #print(("--------c---------"))
    #for i in range(n+1):
    #    print((c[i,:]))
    #print(("--------c---------"))

    #specie_concen = [sum(mult[ind[0]:ind[1]]) for ind in specie_site_index_map]
    #total_c = [sum(c[ind[0]:ind[1]]) for ind in specie_site_index_map]
    total_c = []
    for ind in specie_site_index_map:
        total_c.append(sum([sum(c[i,:]) for i in range(*ind)]))
    #total_c = [sum(c[i,:]) for i in range(n)]
    c_ratio = [total_c[i]/sum(total_c) for i in range(m)]
    print(('-------c_ratio-------------'))
    for i in range(m):
        print((c_ratio[i]))

    # Expression for Omega, the Grand Potential
    omega = e0 - sum([mu[site_mu_map[i]]*sum(c0[i,:]) for i in range(n+1)])
    for p_r in range(n):
        for epi in range(n):
            sum_mu = sum([mu[site_mu_map[j]]*Integer(
                    dC[j,epi,p_r]) for j in range(n+1)])
            omega -= k_B*T*multiplicity[p_r]*exp(-(dE[epi,p_r]-sum_mu)/(k_B*T))

    print ('omega')
    print (omega)
    def compute_mus():

        def reduce_mu():
            host_concen = 1-solute_concen
            new_c0 = c0.astype(float)
            for i in range(n):
                new_c0[i,i] = host_concen*c0[i,i]
            new_c0[n,n] = 2*solute_concen
            omega = [
                e0-sum([mu[site_mu_map[i]]*sum(new_c0[i,:])
                    for i in range(n+1)])]
            x = solve(omega)
            return x

        # Compute trial mu
        mu_red = reduce_mu()
        #print(('mu_red', mu_red))

        mult = multiplicity
        #for ind in specie_site_index_map:
        #    print(ind[0], ind[1])
        specie_concen = [
            sum(mult[ind[0]:ind[1]]) for ind in specie_site_index_map]
        #print('specie_concen', specie_concen)
        max_host_specie_concen = 1-solute_concen
        host_specie_concen_ratio = [specie_concen[i]/sum(specie_concen)* \
                                    max_host_specie_concen for i in range(m)]
        host_specie_concen_ratio[-1] = solute_concen
        #print('hostspecie_concen_rat', host_specie_concen_ratio)


        y_vect = host_specie_concen_ratio
        #print(('y_vect', y_vect))
        vector_func = [y_vect[i]-c_ratio[i] for i in range(m)]
        vector_func.append(omega)
        #print((vector_func))
        mu_vals = None
        c_val = None
        m_min = -15.0
        if e0 > 0:
            m_max = 10            # Search space needs to be modified
        else:
            m_max = 0
        for m1 in np.arange(m_min,m_max,0.3):
            for m2 in np.arange(m_min,m_max,0.3):
                m0 = mu_red[mu[0]].subs([(mu[1],m1),(mu[2],m2)])
                try:
                    #print(m1,m2)
                    mu_vals = nsolve(vector_func,mu,[m0,m1,m2],module="numpy")
                    #mu_vals = nsolve(vector_func,mu,[m0,m1,m2])
                    # Line needs to be modified to include all mus when n > 2
                except:
                    continue
                break
            if mu_vals:
                mu_vals = [float(mu_val) for mu_val in mu_vals]
                break
        else:
            raise ValueError("Couldn't find mus")
        #print (('mu_vals', mu_vals))
        return mu_vals

    if not trial_chem_pot:
        mu_vals = compute_mus()
    else:
        try:
            mu_vals = [trial_chem_pot[element] for element in specie_order]
        except:
            mu_vals = compute_mus()
    #print(('mu_vals', mu_vals))


    # Compute ymax
    max_host_specie_concen = 1-solute_concen
    mult = multiplicity
    specie_concen = [
            sum(mult[ind[0]:ind[1]]) for ind in specie_site_index_map]
    host_specie_concen_ratio = [specie_concen[i]/sum(specie_concen)* \
                                max_host_specie_concen for i in range(m)]
    host_specie_concen_ratio[-1] = solute_concen
    #print('hostspecie_concen_rat', host_specie_concen_ratio)
    li = specie_site_index_map[0][0]
    hi = specie_site_index_map[0][1]
    comp1_min = sum(multiplicity[li:hi])/sum(multiplicity)* \
                max_host_specie_concen - 0.01
    comp1_max = sum(multiplicity[li:hi])/sum(multiplicity)* \
                max_host_specie_concen + 0.01
    #print(ymin, ymax)
    delta = (comp1_max - comp1_min)/50.0

    #for i in range(len(mu)):
    #    print(mu[i], mu_vals[i])

    # Compile mu's for all composition ratios in the range
    #+/- 1% from the stoichiometry
    result = {}
    for y in np.arange(comp1_min,comp1_max+delta,delta):
        result[y] = []
        y_vect = []
        y_vect.append(y)
        y2 = max_host_specie_concen - y
        y_vect.append(y2)
        y_vect.append(solute_concen)
        #print ('y_vect', y_vect)
        vector_func = [y_vect[i]-c_ratio[i] for i in range(1,m)]
        vector_func.append(omega)
        #try:
        #print (vector_func)
        #print (mu)
        #print (mu_vals)

        x = nsolve(vector_func,mu,mu_vals,module="numpy")
        #x = nsolve(vector_func,mu,mu_vals)
        #except:
        #    del result[y]
        #    continue

        result[y].append(x[0])
        result[y].append(x[1])
        result[y].append(x[2])


    res = []

    #print ('result', result.keys())
    # Compute the concentrations for all the compositions
    for key in sorted(result.keys()):
        mu_val = result[key]
        total_c_val = [total_c[i].subs(dict(zip(mu,mu_val))) \
                for i in range(len(total_c))]
        c_val = c.subs(dict(zip(mu,mu_val)))
        #print ('c_val', c_val)
        # Concentration of first element/over total concen
        res1 = []
        res1.append(float(total_c_val[0]/sum(total_c_val)))

        sum_c0 = sum([c0[i,i] for i in range(n)])
        #print ('c_val_n_0', c_val[n,0])
        #print ('c_val_n_1', c_val[n,1])
        for i in range(n+1):
            for j in range(n):
                if i == j:              # Vacancy
                    #res1.append(float((c0[i,i]-sum(c_val[:,i]))/c0[i,i]))
                    vac_conc = float(exp(-(mu_val[site_mu_map[i]]+dE[i,i])/(k_B*T)))
                    res1.append(vac_conc)
                else:                   # Antisite
                    res1.append(float(c_val[i,j]/c0[j,j]))
        res.append(res1)

    #print ('c0', c0)
    #print ('res', res)
    res = np.array(res)
    #print ('res', res)
    dtype = [(str('x'),np.float64)]+[(str('y%d%d' % (i, j)), np.float64) \
            for i in range(n+1) for j in range(n)]
    res1 = np.sort(res.view(dtype),order=[str('x')],axis=0)

    conc = []
    for i in range(n+1):
        conc.append([])
        for j in range(n):
            conc[i].append([])
    #print(conc)
    for i in range(n+1): # Append vacancies
        for j in range(n):
            y1 = [dat[0][i*n+j+1] for dat in res1]
            conc[i][j] = y1
    #print(type(conc[i][j]))

    plot_data = {}
    """Because all the plots have identical x-points storing it in a
    single array"""
    #for dat in res1:
    #    print dat
    plot_data['x'] = [dat[0][0] for dat in res1]         # x-axis data
    # Element whose composition is varied. For x-label
    plot_data['x_label'] = els[0]+ "_mole_fraction"
    plot_data['y_label'] = "Fraction of {} at {} sites".format(
        solute_specie,els[0])

    y_data = []
    #for i in range(n):      # Vacancy plots
    inds = specie_site_index_map[m-1]
    #print ('inds', inds)
    data1 = np.sum([conc[ind][0] for ind in range(*inds)],axis=0)
    #print ('data1', data1)
    data2 = np.sum([conc[ind][1] for ind in range(*inds)],axis=0)
    #print ('data2', data2)
    frac_data = data1/(data1+data2)
    #print ('frac_data', frac_data)
    frac_data = frac_data.tolist()
    y_data.append({'data':frac_data})

    plot_data['y'] = y_data

    return plot_data


@requires(sympy_found,
          "solute_defect_density requires Sympy module. Please install it.")
def solute_defect_density(structure, e0, vac_defs, antisite_defs, solute_defs,
        solute_concen=0.01, T=800, trial_chem_pot = None, 
        plot_style="highchargs"):
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

    solute_conc_data = solute_site_preference_finder(
        structure, e0, T, vac_defs, antisite_defs, solute_defs,
        solute_concen=solute_concen, trial_chem_pot=trial_chem_pot)

    if plot_style == 'highcharts':
        "Energy data is ignored in this mode"
        hgh_chrt_data = {}
        hgh_chrt_data['xAxis'] = conc_data['x_label']
        hgh_chrt_data['yAxis'] = conc_data['y_label']

        series = []
        x = conc_data['x']
        for y_data in solute_conc_data['y']:
            y = y_data['data']
            xy = zip(x,y)
            xy = [list(el) for el in xy]
            name = y_data['name'].strip('$')
            flds= name.split('_')
            def_string = flds[0]
            site_string = flds[1].strip('{}')
            name = def_string+"<sub>"+site_string+"</sub>"
            #series.append({'data':xy, 'name':y_data['name']})
            series.append({'data':xy, 'name':name})
        hgh_chrt_data['series'] = series
        return hgh_chrt_data
    elif plot_style == 'gnuplot':
        def data_to_rows(inp_data):
            rows = []
            labels = []
            labels.append(inp_data['x_label'])
            labels.append(inp_data['y_label'])
            rows.append('#'+'\t'.join(labels))
            m = len(inp_data['x'])
            for i in range(m):
                data = []
                data.append(inp_data['x'][i])
                data += [y['data'][i] for y in inp_data['y']]
                data = [float(x) for x in data]
                rows.append('\t'.join(list(map(str,data))))
            return rows
        conc_rows = data_to_rows(solute_conc_data)
        return conc_rows

