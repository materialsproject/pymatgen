'''
This module performs common batches of transformations on a transmuter object.
Basically just wrappers to perform common sets of operations
'''
from pymatgen.transformations.standard_transformations import SubstitutionTransformation, ChargeBalanceTransformation, OrderDisorderedStructureTransformation

class MultipleSubstitutionTransformation(object):
    '''
    Performs multiple substitutions on a transmuter. For example, can do a fractional replacement of Ge in LiGePS with a list of species,
    creating one structure for each substitution. Ordering is done using a dummy element so only one ordering must be done per substitution
    oxidation state. Charge balancing of the structure is optionally performed.
    
    Note:
    There are no checks to make sure that removal fractions are possible and rounding may occur
    Currently charge balancing only works for removal of species
    '''
    
    def __init__(self, sp_to_replace, r_fraction, substitution_dict, charge_balance_species = None):
        '''
        Performs multiple fractional substitutions on a transmuter
        
        Args:
            sp_to_replace
                species to be replaced
            r_fraction
                fraction of that specie to replace
            substitution_dict
                dictionary of the format
                {2: ["Mg", "Ti", "V", "As", "Cr", "Ta", "N", "Nb"], 
                3: ["Ru", "Fe", "Co", "Ce", "As", "Cr", "Ta", "N", "Nb"],
                4: ["Ru", "V", "Cr", "Ta", "N", "Nb"], 
                5: ["Ru", "W", "Mn"]
                }
                The number is the charge used for each of the list of elements (an element can be present in multiple lists)
            charge_balance_species:
                if specified, will balance the charge on the structure using that specie
        '''
        self._sp_to_replace = sp_to_replace
        self._r_fraction = r_fraction
        self._substitution_dict = substitution_dict
        self._charge_balance_species = charge_balance_species
        
    def apply_transformation(self, transmuter, order = True):
        substitutions = []
        for x in self._substitution_dict.keys():
            mapping = {}
            if x>0:
                sign = '+'
            else:
                sign = '-'
            mapping[self._sp_to_replace] ={self._sp_to_replace:1-self._r_fraction, 'X{}{}'.format(str(x),sign):self._r_fraction}
            substitutions.append(SubstitutionTransformation(mapping))
        transmuter.branch_collection(substitutions)
        if self._charge_balance_species is not None:
            cbt = ChargeBalanceTransformation(self._charge_balance_species)
            transmuter.append_transformation(cbt)
        if order:
            odst = OrderDisorderedStructureTransformation()
            transmuter.append_transformation(odst)
        
        substitutions = []
        for oxistate, element_list in self._substitution_dict.items():
            for element in element_list:
                if x>0:
                    sign = '+'
                else:
                    sign = '-'
                substitutions.append(SubstitutionTransformation( {'X{}+'.format(str(oxistate)): '{}{}{}'.format(element,str(oxistate),sign)}) )
        transmuter.branch_collection(substitutions)
        
            