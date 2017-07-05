#!/usr/bin/env python
"""Processing the steering files
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

# The path to the steering files, after installation
def get_template_path():
    """Get the path where the steering files templates are
    installed. Be careful, this functions highly depends of the
    setup.py entry `package_data`

    Return
    ------
    str: the local path folder where the steering files live
    """
    import os
    return os.path.join(os.path.dirname(__file__),'steering_files')

# Available arguments present in the steering file templates.
# These arguments are dynamically changed on run time
_ARGUMENTS = { 'ROOT_FILENAME': 'Name of the output root file created by the AIDA processor',
        'RUN_NUMBER': 'The run number of the input file',
        'ALIBAVA_INPUT_FILENAME': 'The input file name (ALIBAVA RAW data)',
        'INPUT_FILENAMES': 'The list of input file names (LCIO DATA)',
        'OUTPUT_FILENAME': 'Name of the output LCIO file created by the LCIOOutputProcessor',
        'GEAR_FILE': 'The name of the gear file to be used',
        'PEDESTAL_OUTPUT_FILENAME': 'Name of the output LCIO file created by the AlibavaPedestalNoiseProcessor',
        }

# Marlin step class definition
class marlin_step(object):
    """Abstract class to define what steering file is associated
    with a given step and what arguments should be provided
    to fill the template steering file. 

    Attributes
    ----------
    step_name: str
        The name of the step, it corresponds to the name of the 
        concrete class
    token: str
        The token to be used to define the substitutible arguments 
        in the template steering file
    steering_file: str
        The absolute path name of the steering file created by 
        the class
    steering_file_template: str
        The absolute path name of the template associated
    steering_file_content: str (property)
        The proper content of the steering file
    required_arguments: tuple(str)
        The substitutible arguments of the template steering file
    argument_values: dict(str,GenericType)
        The map bewteen the argument template name and the concrete
        value used
    """
    def __init__(self,step_name):
        """
        Parameters
        ----------
        step_name: str
            name of the step, which corresponds a concrete class
            name 
        """
        import os
        self.step_name = step_name
        self.token = '@'
        self.steering_file = os.path.join(os.getcwd(),"{0}.xml".format(self.step_name))
        # To be fill with the concrete class
        # The complete path of the template file
        self.steering_file_template = None
        # A list of the required arguments and its default values
        self.required_arguments     = None
        # Dict mapping the required argument names with the particular
        # and concrete values used
        self.argument_values        = {}

    @property
    def steering_file_content(self):
        """Getter for the steering file template content.
        If the data member still does not exist, if not, it will
        create it by reading the content of the steering_file_template
        file, otherwise it will return just the content of 
        the data member

        Return
        ------
        steering_file_content: str
            the content of the steering file template
        """
        try:
            return self._steering_file_content
        except AttributeError:
            with open(self.steering_file_template) as f:
                self._steering_file_content =  f.read()
            return self._steering_file_content
    @steering_file_content.setter
    def steering_file_content(self,newcontent):
        """Update the content of the steering file
        
        Parameters
        ----------
        newcontent: str
            the new content
        """
        self._steering_file_content = newcontent
    
    def required_arguments_description(self,argument):
        """The description of the argument

        Parameters
        ----------
        argument: str|int
            the argument name or the argument index 

        Raises
        ------
        KeyError
            If the argument does not present in _ARGUMENTS 
            static dict
        IndexError:
            If the argument index is higher than the number
            of arguments present
        """
        if type(argument) == str:
            return _ARGUMENTS[argument]
        elif type(argument) == int:
            return _ARGUMENTS[self.required_arguments(argument)]

    def set_argument_value(self,argument,value):
        """Substitute the argument by the introduced value
        in the steering file and fill the argument dict, mapping
        the argument name with its concrete value

        Parameters
        ----------
        argument: str
            the name of the argument
        value: Any python type
            the value to be substituted

        Raises
        ------
        KeyError
            If the argument is not defined in the _ARGUMENTS static
            dict
        """
        self.steering_file_content = self.steering_file_content.replace("{0}{1}{0}".format(self.token,argument),str(value))
        self.argument_values[argument] = value

    def get_default_argument_value(self,argument):
        """Get the argument value per default

        Parameters
        ----------
        argument: str
            the name of the argument

        Return
        ------
        str: the default value

        Raises
        ------
        RuntimeError
            If the argument must be set by the user
        """
        import os

        if argument == 'ROOT_FILENAME':
            return self.step_name
        elif argument == 'RUN_NUMBER':
            return -1
        elif argument == 'GEAR_FILE':
            return 'gear_dummy.xml'
        elif argument == 'ALIBAVA_INPUT_FILENAME':
            pass
        elif argument == 'INPUT_FILENAMES':
            pass
        elif argument == 'OUTPUT_FILENAME':
            if self.argument_values.has_key('ALIBAVA_INPUT_FILENAME'):
                return os.path.basename(self.argument_values['ALIBAVA_INPUT_FILENAME'].replace('.dat','.slcio'))
            elif self.argument.values.has_key('INPUT_FILENAMES'):
                return os.path.basename(self.argument_values['INPUT_FILENAMES'].replace('.slcio','_{0}.slcio'.format(self.step_name)))
        elif argument == 'PEDESTAL_OUTPUT_FILENAME':
            return os.path.basename(self.argument_values['INPUT_FILENAMES'].replace('.slcio','_{0}_PEDESTALFILE.slcio'.format(self.step_name)))
               
        raise RuntimeError('Argument "{0}" must be explicitely set'.format(argument))

    def publish_steering_file(self,**kwd):
        """Creates a copy of the steering file with the particular
        values introduced substituted. If a given argument is not
        provided, the default value will be used

        Parameters
        ----------

        Raises
        ------
        NotImplementedError
            If any of the introduced arguments is not defined in 
            the _ARGUMENTS static dictionary
        """
        # Check for inconsistencies
        for key in kwd.keys():
            if key not in _ARGUMENTS.keys():
                raise NotImplementedError("The argument '{0}' is not implemented".format(key))
        # Setting the values    
        for arg,value in kwd.iteritems():
            self.set_argument_value(arg,value)
        # Those arguments not set by the user, set using the
        # default values
        for argument in filter(lambda _a: _a not in kwd.keys(),self.required_arguments):
            val = self.get_default_argument_value(argument)
            self.set_argument_value(argument,val)

        # Create the steering file
        with open(self.steering_file,"w") as f:
            f.write(self.steering_file_content)


# Marlin step concrete implementations
# ------------------------------------

class pedestal_conversion(marlin_step):
    def __init__(self):
        import os
        super(pedestal_conversion,self).__init__('pedestal_conversion')

        self.steering_file_template = os.path.join(get_template_path(),'01-ab_converter.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'ALIBAVA_INPUT_FILENAME', 'OUTPUT_FILENAME')
    
    @staticmethod
    def get_description():
        return 'Convert RAW binary data into LCIO, pedestal runs'  

class pedestal_preevaluation(marlin_step):
    def __init__(self):
        import os
        super(pedestal_preevaluation,self).__init__('pedestal_preevaluation')

        self.steering_file_template = os.path.join(get_template_path(),'02-ped_preevaluation.xml')
        self.required_arguments = ('ROOT_FILENAME','INPUT_FILENAMES', 'PEDESTAL_OUTPUT_FILENAME')
    
    @staticmethod
    def get_description():
        return 'Estimate pedestal and noise from a gaussian distribution'


# The available marlin_steps classes (ordered)
available_steps = (pedestal_conversion,pedestal_preevaluation)

# Factory to create the concrete marlin step given its name of class
def create_marlin_step(step_name):
    """Creates a marlin_step concrete class instance given the name 
    of the class. 
    
    Parameters
    ----------
    step_name: str
        Name of the concrete marlin_step class 

    Return
    ------
    marlin_step concrete class instance

    Raises
    ------
    NotImplementedError
        If the class name does not implemented
    """
    try:
        theclass = filter(lambda x: x.__name__ == step_name,available_steps)[0]
    except IndexError:
        raise NotImplementedError("Not Implemented class '{0}'".format(step_name))
    return theclass()

