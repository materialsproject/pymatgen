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

from monty.dev import requires
from monty.fractions import gcd


try:
    from sympy import Symbol, nsolve, Integer, Rational, Matrix, exp, solve, Eq
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
            "comute_defect_density requires Sympy module. Please install it.")
def dilute_solution_model(structure, e0, vac_defs, antisite_defs, T, 
        trial_chem_pot = None, generate='plot'):

    """
    Compute the defect densities for a structure based on the input parameters
    using dilute solution model.
    Args:
        structure:
            pymatgen.core.structure.Structure object representing the
            primitive or unitcell of the crystal.
        e0:
            The total energy of the undefected system.
            This is E0 from VASP calculation.
        vac_defs:
            List of vacancy defect parameters in the dictionary format.
            The keys of the dict associated with each vacancy defect are
            1) site_index, 2) site_specie, 3) site_multiplicity, and
            4) energy. 1-3 can be obtained from
            pymatgen.analysis.defects.point_defects.Vacancy class.
            Site index is expected to start with 1 (fortran index).
        antisite_defs:
            List of antisite defect parameters in the dictionary format.
            The keys of the dict associated with each antisite defect are
            1) site_index, 2) site_specie, 3) site_multiplicity,
            4) substitution_specie, and 5) energy. 1-3 can be obtained
            from pymatgen.analysis.defects.point_defects.Vacancy class.
        T:
            Temperature in Kelvin
        trial_chem_pot:
            Trial chemical potentials to speedup the plot generation
            Format is {el1:mu1,...}
        generate (string): Options are plot or energy
            Chemical potentials are also returned with energy option.
            If energy option is not chosen, plot is generated.
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
    #print multiplicity
    m = len(set(site_species))      # distinct species
    n = len(vac_defs)           # inequivalent sites

    # Reduce the system and associated parameters such that only distinctive
    # atoms are retained
    comm_div = gcd(*tuple(multiplicity))
    multiplicity = [val/comm_div for val in multiplicity]
    e0 = e0/comm_div
    T = Integer(T)

    c0 = np.diag(multiplicity)
    #print 'c0', c0
    mu = [Symbol('mu'+str(i)) for i in range(m)]

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

    #print 'specie_site_index_map', specie_site_index_map
    #for el in specie_site_index_map:
    #    print range(*el)
    #print site_species
    #print 'site_mu_map', site_mu_map
    #print specie_order


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
                if i == j:# and site_species[j] != site_species[k]:
                    dC[i,j,k] = 1
        for j in range(n):
            for k in range(n):
                if i == k:
                    #if j == k or site_species[j] != site_species[k]:
                        dC[i,j,k] = -1

    #print dC
    # dE matrix: Flip energies (or raw defect energies)
    els = [vac_def['site_specie'] for vac_def in vac_defs]
    dE = []
    for i in range(n):
        dE.append([])
    for i in range(n):
        for j in range(n):
            dE[i].append(None)

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
                        if as_def['site_index'] == j+1 and \
                                sub_specie == as_def['substitution_specie']:
                            dE[i][j] = as_def['energy']
                            break
    dE = np.array(dE)
    np.where(dE == None, dE, 0)
    #print dE

    # Initialization for concentrations
    # c(i,p) == presence of ith type atom on pth type site
    c = Matrix(n,n,[0]*n**2)
    for i in range(n):
        for p in range(n):
            c[i,p] = Integer(c0[i,p])
            for epi in range(n):
                sum_mu = sum([mu[site_mu_map[j]]*Integer(
                        dC[j,epi,p]) for j in range(n)])
                c[i,p] += Integer(multiplicity[p]*dC[i,epi,p]) * \
                        exp(-(dE[epi,p]-sum_mu)/(k_B*T))

    #specie_concen = [sum(mult[ind[0]:ind[1]]) for ind in specie_site_index_map]
    #total_c = [sum(c[ind[0]:ind[1]]) for ind in specie_site_index_map]
    total_c = []
    for ind in specie_site_index_map:
        total_c.append(sum([sum(c[i,:]) for i in range(*ind)]))
    #total_c = [sum(c[i,:]) for i in range(n)]
    c_ratio = [total_c[-1]/total_c[i] for i in range(m)]
    #print 'c_ratio'
    #for i in range(len(c_ratio)):
        #print c_ratio[i]

    # Expression for Omega, the Grand Potential
    omega = e0 - sum([mu[site_mu_map[i]]*sum(c0[i,:]) for i in range(n)])
    for p_r in range(n):
        for epi in range(n):
            sum_mu = sum([mu[site_mu_map[j]]*Integer(
                    dC[j,epi,p_r]) for j in range(n)])
            omega -= k_B*T*multiplicity[p_r]*exp(-(dE[epi,p_r]-sum_mu)/(k_B*T))

    def compute_mus():

        def reduce_mu():
            omega = [e0 - sum([mu[site_mu_map[i]]*sum(c0[i,:]) for i in range(n)])]
            x = solve(omega)
            return x

        # Compute trial mu
        mu_red = reduce_mu()
        #print mu_red

        mult = multiplicity
        #for ind in specie_site_index_map:
        #    print ind[0], ind[1]
        specie_concen = [sum(mult[ind[0]:ind[1]]) for ind in specie_site_index_map]
        #print 'specie_concent', specie_concen
        y_vect = [specie_concen[-1]/specie_concen[i] for i in range(m)]
        #print 'y_vect', y_vect
        vector_func = [y_vect[i]-c_ratio[i] for i in range(m-1)]
        vector_func.append(omega)
        #print vector_func
        #vector_func.append(mu_equalities)
        #print 'y0', y0
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
                # Line needs to be modified to include all mus when n > 2
            except:
                continue

            c_val = c.subs(dict(zip(mu,x)))
            #print c_val
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
        print mu_vals
        return mu_vals
        #print els

    def compute_def_formation_energies():
        i = 0
        for vac_def in vac_defs:
            site_specie = vac_def['site_specie']
            ind = specie_order.index(site_specie)
            uncor_energy = vac_def['energy']
            formation_energy = uncor_energy + mu_vals[ind]
            print site_specie, 'vancancy formation_energy', formation_energy
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
            print site_specie, sub_specie, 'antisite ', formation_energy
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
            mu_vals = [trail_chem_pot[element] for element in specie_ordger]
        except:
            mu_vals = compute_mus()

    if generate == 'energy':
        formation_energies = compute_def_formation_energies()
        mu_dict = dict(zip(specie_order,mu_vals)) 
        return formation_energies, mu_dict


    # Compute ymax
    li = specie_site_index_map[0][0]
    hi = specie_site_index_map[0][1]
    comp1_min = int(sum(multiplicity[li:hi])/sum(multiplicity)*100)-1
    comp2_max = 100-comp1_min
    ymax = comp2_max/comp1_min
    comp1_max = int(sum(multiplicity[li:hi])/sum(multiplicity)*100)+1
    comp2_min = 100-comp1_max
    ymin = comp2_min/comp1_max
    #print ymin, ymax
    delta = (ymax-ymin)/40.0

    #for i in range(len(mu)):
    #    print mu[i], mu_vals[i]

    # Compile mu's for all composition ratios in the range 
    #+/- 1% from the stoichiometry
    result = {}
    for y in np.arange(ymin,ymax,delta):
        result[y] = []
        vector_func = [y-c_ratio[0]]
        vector_func.append(omega)
        x = nsolve(vector_func,mu,mu_vals,module="numpy")
        result[y].append(x[0])
        result[y].append(x[1])


    res = []

    # Compute the concentrations for all the compositions
    for key in result:
        mu_val = result[key]
        total_c_val = [total_c[i].subs(dict(zip(mu,mu_val))) \
                for i in range(len(total_c))]
        c_val = c.subs(dict(zip(mu,mu_val)))
        res1 = []
        # Concentration of first element/over total concen
        res1.append(float(total_c_val[0]/sum(total_c_val)))    
        sum_c0 = sum([c0[i,i] for i in range(n)])
        for i in range(n):
            for j in range(n):
                if i == j:              # Vacancy
                    res1.append(float((c0[i,i]-sum(c_val[:,i]))/c0[i,i]))     
                else:                   # Antisite
                    res1.append(float(c_val[i,j]/c0[i,i]))                    
        res.append(res1)

    res = np.array(res)
    dtype = [('x',np.float64)]+[('y'+str(i)+str(j),np.float64) \
            for i in range(n) for j in range(n)]
    res1 = np.sort(res.view(dtype),order=['x'],axis=0)


    plot_data = {}
    """Because all the plots have identical x-points storing it in a 
    single array"""
    plot_data['x'] = [dat[0][0] for dat in res1]         # x-axis data
    # Element whose composition is varied. For x-label
    plot_data['x_label'] = els[0]+ " mole fraction" 
    plot_data['y_label'] = "Point defect concentration"
    conc = []
    for i in range(n):
        conc.append([])
        for j in range(n):
            conc[i].append([])
    #print conc
    for i in range(n): # Append vacancies
        for j in range(n):
            y1 = [dat[0][i*n+j+1] for dat in res1]
            conc[i][j] = y1
    #print type(conc[i][j])

    y_data = []
    for i in range(n):      # Vacancy plots
        data = conc[i][i]
        specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1]-indices[0]
        cur_ind = i - indices[0] + 1
        if 'V' not in els:
            vac_string = "$V_{"
        else:
            vac_string = "$Vac_{"
        if not specie_ind_del-1:
            label = vac_string+specie+'}$'
        else:
            label = vac_string+specie+'_'+str(cur_ind)+'}$'
        # Plot data and legend info
        y_data.append({'data':data,'name':label})       

        site_specie = els[i]
        for j in range(m):          # Antisite plot dat
            sub_specie = specie_order[j]
            if sub_specie == site_specie:
                continue
            if not specie_ind_del-1:
                label = '$'+sub_specie+'_{'+specie+'}$'
            else:
                label = '$'+sub_specie+'_{'+specie+'_'+str(cur_ind)+'}$'
            inds = specie_site_index_map[j]
            data = np.sum([conc[ind][i] for ind in range(*inds)],axis=0)
            data = data.tolist()
            y_data.append({'data':data,'name':label})

    plot_data['y'] = y_data

    return plot_data


def compute_defect_density(structure, e0, vac_defs, antisite_defs, T=800, 
        trial_chem_pot=None, plot_style="HighCharts"):
    """
    Wrapper for the dilute_solution_model where the computed plot data is 
    prepared based on plot_style. Only "HighCharts" is supported at this point

    :param structure:
    :param e0:
    :param vac_defs:
    :param antisite_defs:
    :param T:
    :param plot_style:
    :return:
    """
    plot_data = dilute_solution_model(structure,e0,vac_defs,antisite_defs,T,
            trial_chem_pot=trial_chem_pot)

    if plot_style == 'HighCharts':
        hgh_chrt_data = {}
        hgh_chrt_data['xAxis'] = plot_data['x_label']
        hgh_chrt_data['yAxis'] = plot_data['y_label']

        series = []
        x = plot_data['x']
        for y_data in plot_data['y']:
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
