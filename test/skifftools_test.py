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
    
    def test_required_arguments_from_template_file(self):
        import os
        import re
        for _step in steering_processing.available_steps:
            a = _step()
            if isinstance(a,steering_processing.alibava_full_reco) \
                or isinstance(a,steering_processing.telescope_full_reco):
                # No sense for the metaclass alibava/telescope full reco
                continue
            # Extract the arguments to be substituted present in the templates
            # Note we need to split first potential double arguments 
            # (as @ARG1@_@ARG2@ -> @ARG1@ @ARG2@)
            arglist=re.findall("@(.*)@",a.steering_file_content.replace("@_@","@\n@").replace("@ @","@\n@"))
            for arg in arglist:
                self.assertTrue(arg in a.required_arguments,\
                        "Argument '{0}' present in the template file '{1}'"\
                        " is not explicitely required in the class '{2}'"\
                        .format(arg,os.path.basename(a.steering_file_template),a.step_name))
if __name__ == "__main__":
    unittest.main()

