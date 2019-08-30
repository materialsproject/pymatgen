from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import ruamel.yaml

from monty.serialization import loadfn
from adjustText import adjust_text

from pymatgen import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.reaction_calculator import ComputedReaction


class CorrectionCalculator:
    
    species = ['oxide', 'peroxide', 'superoxide', 'F', 'Cl', 'Br', 'I', 'N', 'S', 'Se',\
               'Si', 'Sb', 'Te', 'V', 'Cr', 'Mn', 'Fe', 'Co','Ni', 'Cu', 'Mo'] #species that we're fitting corrections for

    def __init__(self, exp_json, comp_json):
        self.exp_compounds = loadfn(exp_json) #experimental data
        self.calc_compounds = loadfn(comp_json) #computed entries
        self.corrections = []
        self.corrections_std_error = []
        self.corrections_dict = {} #{'compound': (value, error)}

        #these three lists are just to help the graph_residual_error_per_species() method
        self.oxides = [] 
        self.peroxides = []
        self.superoxides = []

    def compute_corrections(self, allow_polyanions=False, allow_large_errors=False, allow_unstable=False):
        self.names = []
        self.diffs = []
        self.coeff_mat = []
        self.exp_uncer = []

        self.mpids = []
        for cmpd_info in self.exp_compounds:
            name = cmpd_info['formula']
            warnings = cmpd_info['warnings']

            if allow_polyanions:
                warnings.pop('polyanion', None)
            if allow_large_errors:
                warnings.pop('large_uncertainty', None)
            if allow_unstable:
                warnings.pop('unstable', None)

            if name in self.calc_compounds and not warnings:
                
                comp = Composition(name)
                elems = list(comp.as_dict())
                
                compound = self.calc_compounds[name]
                
                reactants = []
                for elem in elems:
                    try:
                        reactants.append(self.calc_compounds[elem])
                    except KeyError:
                        raise ValueError('Computed entries missing ' + elem)
                
                rxn = ComputedReaction(reactants, [compound])
                rxn.normalize_to(comp)
                energy = rxn.calculated_reaction_energy
                
                if compound.data['oxide_type'] == 'oxide':
                    coeff = [comp['O'], 0, 0]
                    self.oxides.append(name)
                elif compound.data['oxide_type'] == 'peroxide':
                    coeff = [0, comp['O'], 0]
                    self.peroxides.append(name)
                elif compound.data['oxide_type'] == 'superoxide':
                    coeff = [0, 0, comp['O']]
                    self.superoxides.append(name)
                else:
                    coeff = [0, 0, 0]
                coeff += [comp[elem] for elem in self.species[3:]]
                
                self.names.append(name)
                self.diffs.append((cmpd_info['exp energy'] - energy)/comp.num_atoms)
                self.coeff_mat.append([i/comp.num_atoms for i in coeff])
                self.exp_uncer.append((cmpd_info['uncertainty'])/comp.num_atoms)
                            
                self.mpids.append(compound.entry_id)   

        #for any exp entries with no uncertainty value, assign average uncertainty value
        sigma = np.array(self.exp_uncer)
        sigma[sigma == 0] = np.nan
        mean_uncer = np.nanmean(sigma)
        sigma = np.where(np.isnan(sigma), mean_uncer, sigma)

        f = lambda x, *m: np.dot(x, m)

        popt, pcov = curve_fit(f, self.coeff_mat, self.diffs, p0=np.ones(21), sigma=sigma, absolute_sigma=True)
        self.corrections = popt.tolist()
        self.corrections_std_error = np.sqrt(np.diag(pcov)).tolist()
        for i in range(len(self.species)):
            self.corrections_dict[self.species[i]] = (round(self.corrections[i], 3), round(self.corrections_std_error[i], 4))
        return self.corrections_dict


    def graph_residual_error(self):
        if len(self.corrections) == 0:
            self.compute_corrections()
        
        indices = [i for i in range(len(self.diffs))]
        abs_errors = [abs(i) for i in (self.diffs - np.dot(self.coeff_mat, self.corrections))]
        labels_graph = self.names.copy()
        uncertainty_graph = self.exp_uncer.copy()
        abs_errors, labels_graph, uncertainty_graph = (list(t) for t in zip(*sorted(zip(abs_errors, labels_graph, uncertainty_graph)))) #sort by error
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 20))
        num = len(indices)//2
        ax2.scatter(indices[num:], abs_errors[num:])
        ax2.set_ylim(bottom=ax1.get_ylim()[0])
        texts = []
        for i, txt in enumerate(labels_graph[num:]):
            texts.append(ax2.text(indices[i+num], abs_errors[i+num], txt, fontsize = 12))
        adjust_text(texts, ax = ax2)
        ax1.scatter(indices[:num], abs_errors[:num])
        ax1.set_ylim(top=ax2.get_ylim()[1])
        texts = []
        for i, txt in enumerate(labels_graph[:num]):
            texts.append(ax1.text(indices[i], abs_errors[i], txt, fontsize = 12))
        adjust_text(texts, ax = ax1)
        plt.show()
        print('Residual Error:')
        print('Median = ' + str(np.median(np.array(abs_errors))))
        print('Mean = ' + str(np.mean(np.array(abs_errors))))
        print('Std Dev = ' + str(np.std(np.array(abs_errors))))
        print('Original Error:')
        print('Median = ' + str(abs(np.median(np.array(self.diffs)))))
        print('Mean = ' + str(abs(np.mean(np.array(self.diffs)))))
        print('Std Dev = ' + str(np.std(np.array(self.diffs))))   
        
    def graph_residual_error_per_species(self, specie):
        if specie not in self.species:
            raise Exception('not a valid specie')

        if len(self.corrections) == 0:
            self.compute_corrections()

        abs_errors = [abs(i) for i in (self.diffs - np.dot(self.coeff_mat, self.corrections))]
        labels_species = self.names.copy()
        diffs_cpy = self.diffs.copy()
        num = len(labels_species)
        
        if specie == 'oxide' or specie == 'peroxide' or specie == 'superoxide':
            if specie == 'oxide':
                compounds = self.oxides
            elif specie == 'peroxide':
                compounds = self.peroxides
            else:
                compounds = self.superoxides
            for i in range(num):
                if labels_species[num-i-1] not in compounds:
                    del labels_species[num-i-1]
                    del abs_errors[num-i-1]
                    del diffs_cpy[num-i-1]  
        else:
            for i in range(num):
                if not Composition(labels_species[num-i-1])[specie]:
                    del labels_species[num-i-1]
                    del abs_errors[num-i-1]
                    del diffs_cpy[num-i-1]        
        abs_errors, labels_species, diffs_cpy = (list(t) for t in zip(*sorted(zip(abs_errors, labels_species, diffs_cpy)))) #sort by error
        indices = [i for i in range(len(diffs_cpy))]
        if len(indices) > 20:
            plt.figure(figsize=(20, 10))
        plt.scatter(indices, abs_errors)
        texts = []
        for i, txt in enumerate(labels_species):
            texts.append(plt.text(indices[i], abs_errors[i], txt, fontsize = 12))
        adjust_text(texts)
        plt.show()
        print('Residual Error:')
        print('Median = ' + str(np.median(np.array(abs_errors))))
        print('Mean = ' + str(np.mean(np.array(abs_errors))))
        print('Std Dev = ' + str(np.std(np.array(abs_errors))))
        print('Original Error:')
        print('Median = ' + str(abs(np.median(np.array(diffs_cpy)))))
        print('Mean = ' + str(abs(np.mean(np.array(diffs_cpy)))))
        print('Std Dev = ' + str(np.std(np.array(diffs_cpy)))) 

        print(labels_species) 

    def make_yaml(self, name='MP'):
        """
        create the _name_Compatibility.yaml that stores corrections as well as _name_CompatibilityErrors.yaml for correction errors
        """

        if len(self.corrections) == 0:
            self.compute_corrections()

        #from old mpcompatibility
        aqueous = OrderedDict()
        aqueous['O2'] = -0.316731
        aqueous['N2'] = -0.295729
        aqueous['F2'] = -0.313025
        aqueous['Cl2'] = -0.344373
        aqueous['Br'] = -0.235039
        aqueous['Hg'] = -0.234421
        aqueous['H2'] = -3.6018845
        aqueous['H2O'] = -4.972

        compatibility = OrderedDict()
        anion_corr = OrderedDict()
        advanced = OrderedDict()
        gas_corr = OrderedDict()
        u_corr = OrderedDict()
        o = OrderedDict()
        f = OrderedDict()

        compatibility_error = OrderedDict()
        anion_corr_error = OrderedDict()
        advanced_error = OrderedDict()
        gas_corr_error = OrderedDict()
        u_corr_error = OrderedDict()
        o_error = OrderedDict()
        f_error = OrderedDict()

        anion_corr['oxide'] = self.corrections_dict['oxide'][0]
        anion_corr['peroxide'] = self.corrections_dict['peroxide'][0]
        anion_corr['superoxide'] = self.corrections_dict['superoxide'][0]
        anion_corr['ozonide'] = 0 #do i need this??

        anion_corr_error['oxide'] = self.corrections_dict['oxide'][1]
        anion_corr_error['peroxide'] = self.corrections_dict['peroxide'][1]
        anion_corr_error['superoxide'] = self.corrections_dict['superoxide'][1]
        anion_corr_error['ozonide'] = 0 #do i need this??

        anion_corr['sulfide'] = self.corrections_dict['S'][0]
        anion_corr_error['sulfide'] = self.corrections_dict['S'][1]


        for elem in ['Br', 'I', 'Se', 'Si', 'Sb', 'Te']:
            anion_corr[elem] = self.corrections_dict[elem][0]
            anion_corr_error[elem] = self.corrections_dict[elem][1]


        for elem in ['F', 'Cl', 'N']:
            entry = self.calc_compounds[elem]
            key = entry.composition.reduced_formula
            val = entry.energy_per_atom - self.corrections_dict[elem][0]
            gas_corr[key] = val
            gas_corr_error[key] = self.corrections_dict[elem][1]

        gas_corr['H2'] =  -3.23973666138 #from old mpcompatibility
        gas_corr_error['H2'] = 0

        for elem in ['V', 'Cr', 'Mn', 'Fe', 'Co','Ni', 'Cu', 'Mo']:
            o[elem] = self.corrections_dict[elem][0]
            f[elem] = self.corrections_dict[elem][0]

            o_error[elem] = self.corrections_dict[elem][1]
            f_error[elem] = self.corrections_dict[elem][1]
    
        u_corr['O'] = o
        u_corr['F'] = f
        advanced['UCorrections'] = u_corr
        compatibility['Name'] = name
        compatibility['Advanced'] = advanced
        compatibility['GasCorrections'] = gas_corr
        compatibility['AnionCorrections'] = anion_corr
        compatibility['AqueuousCompoundEnergies'] = aqueous

        u_corr_error['O'] = o_error
        u_corr_error['F'] = f_error
        advanced_error['UCorrections'] = u_corr_error
        compatibility_error['Name'] = name
        compatibility_error['Advanced'] = advanced_error
        compatibility_error['GasCorrections'] = gas_corr_error
        compatibility_error['AnionCorrections'] = anion_corr_error
        

        fn = name + 'Compatibility.yaml'

        f = open(fn, 'w')
        yaml = ruamel.yaml.YAML()
        yaml.Representer.add_representer(OrderedDict, yaml.Representer.represent_dict)
        yaml.default_flow_style = False
        yaml.dump(compatibility, f)
        f.close()

        fn = name + 'CompatibilityErrors.yaml'
        f = open(fn, 'w')
        yaml = ruamel.yaml.YAML()
        yaml.Representer.add_representer(OrderedDict, yaml.Representer.represent_dict)
        yaml.default_flow_style = False
        yaml.dump(compatibility_error, f)
        f.close()