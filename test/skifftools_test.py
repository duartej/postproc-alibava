#!/usr/bin/env python
"""
Suite of tests for the alibavaSkifftools package
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "v0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

import unittest
from alibavaSkifftools import steering_processing

class TestProcessingModule(unittest.TestCase):
    """Test the steering_processing module
    """
    def test_correct_template_path(self):
        import os
        self.assertTrue(os.path.isdir(steering_processing.get_template_path()))
        
    def test_available_marlin_steps(self):
        for _step in steering_processing.available_steps:
            self.assertTrue(issubclass(_step,steering_processing.marlin_step))

    def test_required_arguments(self):
        for _step in steering_processing.available_steps:
            a = _step()
            for req_arg in a.required_arguments:
                self.assertTrue(steering_processing._ARGUMENTS.has_key(req_arg),\
                        "Argument '{0}' required in class '{1}' but not present"\
                        " in _ARGUMENTS static dict".format(req_arg,a.step_name))

    def test_required_arguments_present_in_template(self):
        import os
        for _step in steering_processing.available_steps:
            a = _step()
            for req_arg in a.required_arguments:
                self.assertTrue(a.steering_file_content.find(req_arg)!=-1,\
                        "Argument '{0}' required in class '{1}' but not present"\
                        " in the associated template file '{2}'".format(req_arg,a.step_name,\
                        os.path.basename(a.steering_file_template)))
    
if __name__ == "__main__":
    unittest.main()

