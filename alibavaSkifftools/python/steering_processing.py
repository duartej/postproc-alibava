#!/usr/bin/env python
"""Processing the steering files
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

# The path to the steering files
import os
the_steering_file=os.path.join(os.path.dirname(__file__),'steering_files','01-ab_converter.xml')
# -- FIXME 


# Dummy class to define what steering file is associated
# with a given step and what arguments should be provided
# to fill the template steering file
class marlin_step(object):
    def __init__(self,step_name):
        self.step_name = step_name
        self.token = '@'

# Concrete implementations
class pedestal_conversion(marlin_step):
    def __init__(self,step_name):
        super(pedestal_conversion).__init__(step_name)
        self.steering_file = '01-ab_converter.xml'
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAME', 'OUTPUT_FILENAME')
        self.required_arguments_description = (\
                'Name of the output root file created by the AIDA processor',\
                'The run number of the input file', 'The input file name (ALIBAVA RAW data)',\
                'Name of the output LCIO file created by the LCIOOutputProcessor')

# Static dict mapping the step with the template steering 
# files and their available word-substitutions
_wordsubs_list = {
        'pedestal-conversion': '01-ab_converter.xml',
        }
