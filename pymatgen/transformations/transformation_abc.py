#!/usr/bin/env python

'''
Defines an abstract base class contract for Transformation object.
'''

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"


import abc
import json


class AbstractTransformation(object):
    """
    Abstract transformation class.    
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def apply_transformation(self, structure):
        '''
        Applies the transformation to a structure.
        
        Args:
            structure:
                input structure
                
        Returns:
            Transformed structure
        '''
        return

    @abc.abstractproperty
    def inverse(self):
        '''
        Returns the inverse transformation if available.
        Otherwise, should return None.
        '''
        return

    @abc.abstractproperty
    def is_one_to_many(self):
        '''
        Determines if a Transformation is a one-to-many transformation. If a
        Transformation is a one-to-many transformation, the apply_transformation
        method should have a keyword arg "return_ranked_list" which allows for
        the transformed structures to be returned as a ranked list.
        '''
        return False

    def to_json(self):
        return json.dumps(self.to_dict)


    @abc.abstractproperty
    def to_dict(self):
        '''
        Creates a json representation of the transformation, which must contain
        all necessary information to allow a reconstruction from a from_json static method.
        Format should follow:
        {'name' : transformation_class_name, 'init_arguments' : (init arguments)}
        '''
        return
