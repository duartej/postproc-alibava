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
        'ALIBAVA_INPUT_FILENAME': 'The input file name (ALIBAVA RAW data or slcio for the merger)',
        'TELESCOPE_INPUT_FILENAME': 'The input file name (ACONITE Telescope RAW data or slcio for the merger)',
        'INPUT_FILENAMES': 'The list of input file names (LCIO DATA)',
        'OUTPUT_FILENAME': 'Name of the output LCIO file created by the LCIOOutputProcessor',
        'GEAR_FILE': 'The name of the gear file to be used',
        'PEDESTAL_OUTPUT_FILENAME': 'Name of the output LCIO file created by the AlibavaPedestalNoiseProcessor',
        'PEDESTAL_INPUT_FILENAME': 'Name of the input LCIO file created by the AlibavaPedestalNoiseProcessor,'\
                ' containing the pedestals',
        'CALIBRATION_OUTPUT_FILENAME': 'Name of the output LCIO file created by the AlibavaCalibrateProcessor',
        'MAXADC':'Max ADCs counts for the common mode histograms (600 per default)',
        'MINADC':'Min ADCs counts for the common mode histograms (400 per default)',
        'NBINS': 'Number of bins to be used in the histograms (200 default)',
        'TIMECUT_MIN': 'The minimum TDC time that is acceptable to use an event',
        'TIMECUT_MAX': 'The maximum TDC time that is acceptable to use an event',
        'MAXCMMDERR': 'Maximum value for the common mode error histogram',
        'MINCMMDERR': 'Minimum value for the common mode error histogram',
        'CMMDCUT_MIN': 'The minimum common mode noise ADC counts acceptable to use an event',
        'CMMDCUT_MAX': 'The maximum common mode noise ADC counts acceptable to use an event',
        'SNRCUT_SEED': 'The minimum signal to noise ratio that a channel has to pass to be used as cluster seed', 
        'SNRCUT_NGB': 'The minimum signal to noise ratio that a neighbour channel has to pass to be added to a cluster', 
        'SIGNAL_POLARITY': 'The polarity of the signal (-1 for negative signals)',
        'SENSORID_STARTS': 'The sensor ID for the alibava data which will be stored as SENSORID_STARTS+chip_number',
        'MAX_FIRING_FREQ_PIXEL': 'The maximum allowed firing frequency to consider a pixel as hot',
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
        # List of auxiliary files (basename) to be copied from
        # get_template_path() to the cwd
        self.auxiliary_files = []

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
        #import shutil
        from .SPS2017TB_metadata import filename_parser

        if argument == 'ROOT_FILENAME':
            return self.step_name
        elif argument == 'RUN_NUMBER':
            # if RS or beam, the run number can be extracted from
            # the datafile 
            for key in ['ALIBAVA_INPUT_FILENAME', 'TELESCOPE_INPUT_FILENAME','INPUT_FILENAMES' ]:
                if self.argument_values.has_key(key):
                    # First, check if is a telescope file: run000XXX.suffix
                    if os.path.basename(self.argument_values[key]).find('run000') == 0:
                        rawnumber = os.path.basename(self.argument_values[key]).replace('run000','').split('.')[0]
                        return int(rawnumber)
                    # if not telescope data, then check the standard alibava 
                    # filenames with the parser
                    fnp = filename_parser(self.argument_values[key])
                    if fnp.is_beam:
                        return fnp.run_number
            return -1
        elif argument == 'GEAR_FILE':
            # First copy the file to the cwd
            self.auxiliary_files.append('dummy_gear.xml')
            return 'dummy_gear.xml'
        elif argument == 'ALIBAVA_INPUT_FILENAME':
            pass
        elif argument == 'TELESCOPE_INPUT_FILENAME':
            pass
        elif argument == 'INPUT_FILENAMES':
            pass
        elif argument == 'PEDESTAL_INPUT_FILENAME':
            pass
        elif argument == 'OUTPUT_FILENAME':
            # The merger use of the the Alibava and Telescope input filenames
            if self.argument_values.has_key('ALIBAVA_INPUT_FILENAME') and \
                    self.argument_values.has_key('TELESCOPE_INPUT_FILENAME'):
                # Get the parser and use
                fnp = filename_parser(self.argument_values['ALIBAVA_INPUT_FILENAME'])
                return "{0}_{1}_{2}_{3}_{4}_DATAMERGED.slcio".format(fnp.run_number,fnp.self.date,self.sensor_name,\
                        self.voltage_bias,self.temperature)
            # The raw binary use of the arguments
            if self.argument_values.has_key('ALIBAVA_INPUT_FILENAME'):
                return os.path.basename(self.argument_values['ALIBAVA_INPUT_FILENAME'].replace('.dat','.slcio'))
            elif self.argument_values.has_key('TELESCOPE_INPUT_FILENAME'):
                return os.path.basename(self.argument_values['TELESCOPE_INPUT_FILENAME'].replace('.raw','.slcio'))
            elif self.argument_values.has_key('INPUT_FILENAMES'):
                return os.path.basename(self.argument_values['INPUT_FILENAMES'].replace('.slcio','.{0}.slcio'.format(self.step_name.replace("_","."))))
        elif argument == 'PEDESTAL_OUTPUT_FILENAME':
            return os.path.basename(self.argument_values['INPUT_FILENAMES'].replace('.slcio','_{0}_PEDESTALFILE.slcio'.format(self.step_name)))
        elif argument == 'CALIBRATION_OUTPUT_FILENAME':
            return os.path.basename(self.argument_values['INPUT_FILENAMES'].replace('.slcio','_{0}_CALIBRATIONFILE.slcio'.format(self.step_name)))
        elif argument == 'MAXADC':
            return 600.0
        elif argument == 'MINADC':
            return 400.0
        elif argument == 'NBINS':
            return 200
        elif argument == 'TIMECUT_MIN':
            return 3.0
        elif argument == 'TIMECUT_MAX':
            return 30.0
        elif argument == 'MAXCMMDERR':
            return 20.0
        elif argument == 'MINCMMDERR':
            return 0.0
        elif argument == 'CMMDCUT_MIN':
            return -10.0
        elif argument == 'CMMDCUT_MAX':
            return 10.0
        elif argument == 'SNRCUT_SEED':
            return 5
        elif argument == 'SNRCUT_NGB':
            return 3
        elif argument == 'SIGNAL_POLARITY':
            return -1
        elif argument == 'SENSORID_STARTS':
            return 5
        elif argument == 'MAX_FIRING_FREQ_PIXEL':
            return 0.001
               
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
        import shutil
        import os

        # Check for inconsistencies
        for key in kwd.keys():
            if key not in _ARGUMENTS.keys():
                raise NotImplementedError("The argument '{0}' is not implemented".format(key))
        # Setting the values provided by the user
        for arg,value in kwd.iteritems():
            self.set_argument_value(arg,value)
        # Those arguments not set by the user (i.e. not present in the
        # argument_values dictionary -- this allows to define tuned defaults
        # depending of the step), set using the default values
        #for argument in filter(lambda _a: _a not in kwd.keys(),self.required_arguments):
        for argument in filter(lambda _a: _a not in self.argument_values.keys(),self.required_arguments):
            val = self.get_default_argument_value(argument)
            self.set_argument_value(argument,val)
        # And the tuned default values (note that if the user already set it,
        # the next lines do not affect anything as the @KEY@ is not present
        # anymore in steering_file_content
        for argdef,val in self.argument_values.iteritems():
            self.set_argument_value(argdef,val)

        # Create the steering file
        with open(self.steering_file,"w") as f:
            f.write(self.steering_file_content)

        # And copy the needed auxiliary file, if any
        for _f in self.auxiliary_files:
            shutil.copyfile(os.path.join(get_template_path(),_f), os.path.join(os.getcwd(),_f))



# Marlin step concrete implementations
# ------------------------------------

class pedestal_conversion(marlin_step):
    def __init__(self):
        import os
        super(pedestal_conversion,self).__init__('pedestal_conversion')

        self.steering_file_template = os.path.join(get_template_path(),'01-ab_converter.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'ALIBAVA_INPUT_FILENAME', 'OUTPUT_FILENAME','GEAR_FILE')
    
    @staticmethod
    def get_description():
        return 'Convert RAW binary data into LCIO, pedestal runs'  

class pedestal_preevaluation(marlin_step):
    def __init__(self):
        import os
        super(pedestal_preevaluation,self).__init__('pedestal_preevaluation')

        self.steering_file_template = os.path.join(get_template_path(),'02-ped_preevaluation.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_OUTPUT_FILENAME','GEAR_FILE')
    
    @staticmethod
    def get_description():
        return 'Estimate pedestal and noise from a gaussian distribution'

class cmmd_calculation(marlin_step):
    def __init__(self):
        import os
        super(cmmd_calculation,self).__init__('cmmd_calculation')

        self.steering_file_template = os.path.join(get_template_path(),'03-ped_cmmd_calculation.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_INPUT_FILENAME',\
                'OUTPUT_FILENAME','MAXADC','MINADC','NBINS','MAXCMMDERR','MINCMMDERR','GEAR_FILE')
    
    @staticmethod
    def get_description():
        return 'Common mode noise calculation'

class pedestal_evaluation(marlin_step):
    def __init__(self):
        import os
        super(pedestal_evaluation,self).__init__('pedestal_evaluation')

        self.steering_file_template = os.path.join(get_template_path(),'04-ped_evaluation.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_OUTPUT_FILENAME',\
                'MAXADC','MINADC','NBINS','GEAR_FILE')
        # Define a tuned default for the histogram bin and ranges
        self.argument_values['MAXADC']=800.0
        self.argument_values['MINADC']=200.0
        self.argument_values['NBINS']=600
    
    @staticmethod
    def get_description():
        return 'Estimate pedestal and noise from a gaussian distribution (common mode subtracted)'

class calibration_conversion(pedestal_conversion):
    def __init__(self):
        import os
        super(calibration_conversion,self).__init__()
        # Change the step name
        self.step_name='calibration_conversion'
        self.steering_file = self.steering_file.replace('pedestal_conversion',self.step_name)
    
    @staticmethod
    def get_description():
        return 'Convert RAW binary data into LCIO, calibration runs'  

class calibration_extraction(marlin_step):
    def __init__(self):
        import os
        super(calibration_extraction,self).__init__('calibration_extraction')

        self.steering_file_template = os.path.join(get_template_path(),'02-cal_extraction.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'CALIBRATION_OUTPUT_FILENAME',\
                'GEAR_FILE')
    
    @staticmethod
    def get_description():
        return 'Extract calibration constant per channel'

class rs_conversion(marlin_step):
    def __init__(self):
        import os
        super(rs_conversion,self).__init__('rs_conversion')

        self.steering_file_template = os.path.join(get_template_path(),'01-ab_converter_rs.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'ALIBAVA_INPUT_FILENAME', 'OUTPUT_FILENAME','GEAR_FILE',\
                "TIMECUT_MIN","TIMECUT_MAX")
    
    @staticmethod
    def get_description():
        return 'Convert RAW binary data into LCIO, beam or RS runs'  

class signal_reconstruction(marlin_step):
    def __init__(self):
        import os
        super(signal_reconstruction,self).__init__('signal_reconstruction')

        self.steering_file_template = os.path.join(get_template_path(),'02-signal_reconstruction.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_INPUT_FILENAME',\
                'MAXADC','MINADC','NBINS','MAXCMMDERR','MINCMMDERR','CMMDCUT_MIN','CMMDCUT_MAX','GEAR_FILE',\
                'OUTPUT_FILENAME')
        # Define a tuned default for the histogram bin and ranges
        self.argument_values['MAXADC']=1000.0
        self.argument_values['MINADC']=-1000.0
        self.argument_values['NBINS']=2000
        self.argument_values['MAXCMMDERR']=30.0
        self.argument_values['MINCMMDERR']=0.0
    
    @staticmethod
    def get_description():
        return 'Signal reconstruction: pedestal and common mode subtraction, plus common mode cut'

class alibava_clustering(marlin_step):
    def __init__(self):
        import os
        super(alibava_clustering,self).__init__('alibava_clustering')

        self.steering_file_template = os.path.join(get_template_path(),'03-ab_clustering.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_INPUT_FILENAME',\
                'SNRCUT_NGB','SNRCUT_SEED','SIGNAL_POLARITY','GEAR_FILE','SENSORID_STARTS','OUTPUT_FILENAME')
    
    @staticmethod
    def get_description():
        return 'Cluster finding algorithm and conversion from Alibava clusters to EUTelSparseCluster'

# Metaclass to deal with the full reconstruction for ALIBAVA
class alibava_full_reco(marlin_step):
    def __init__(self):
        import os
        super(alibava_full_reco,self).__init__('alibava_full_reco')
        
        # -- Dummy 
        self.required_arguments = () 

        # The list of steps with their needed arguments
        self.step_chain = ( 
                (pedestal_conversion(), { 'ALIBAVA_INPUT_FILENAME': self.pedestal_raw_file},self.update_output),
                (pedestal_preevaluation(), {'INPUT_FILENAMES':self.last_output_filename},self.update_pedestal),
                (cmmd_calculation(), { 'INPUT_FILENAMES': self.last_output_filename, 'PEDESTAL_INPUT_FILENAME': self.pedestal_file },self.update_output),
                (pedestal_evaluation(),{'INPUT_FILENAMES': self.last_output_filename},self.update_pedestal),
                (rs_conversion(), { 'ALIBAVA_INPUT_FILENAME': self.beam_raw_file},self.update_output),
                (signal_reconstruction(), {'INPUT_FILENAMES': self.last_output_filename, 'PEDESTAL_INPUT_FILENAME': self.pedestal_file},self.update_output),
                (alibava_clustering(), {'INPUT_FILENAMES': self.last_output_filename, 'PEDESTAL_INPUT_FILENAME': self.pedestal_file},self.update_output),
                )

    # Some datamembers used in the step_chain are not going to be 
    # populated until the call to the step is performed, using them
    # as properties. Note the use of the updaters (like the setters of
    # a property) methods (defined below, are the mechanism to update
    # this elements)
    def pedestal_raw_file(self):
        return self._pedestal_raw_file
    
    def beam_raw_file(self):
        return self._beam_raw_file

    def last_output_filename(self):
        return self._last_output_filename
    
    def pedestal_file(self):
        return self._pedestal_file
    
    # Updaters
    def update_output(self,step_inst):
        """Update the last_output_filename data member using the value
        from the last step
        """
        self._last_output_filename = step_inst.argument_values['OUTPUT_FILENAME']

    def update_pedestal(self,step_inst):
        self._pedestal_file = step_inst.argument_values['PEDESTAL_OUTPUT_FILENAME']

    @staticmethod
    def get_description():
        return 'Metaclass to perform the whole bunch of steps to process the ALiBaVa data'

    def publish_steering_file(self,**kwd):
        """Creates every steering file needed for reconstruct
        the ALiBaVa data using only the RAW alibava beam and 
        pedestal input filenames. Note that the other values should
        be change manually later in the created steering files

        Parameters
        ----------
        ALIBAVA_INPUT_FILENAME: str
            the name of the raw alibava data beam file
        PEDESTAL_INPUT_FILENAME: str
            the name of the raw alibava data pedestal file

        Raises
        ------
        NotImplementedError
            If any of the introduced arguments is not defined in 
            the _ARGUMENTS static dictionary
        """
        import time
        import datetime
        import os
        import stat
    
        # Check for inconsistencies
        for key in ['ALIBAVA_INPUT_FILENAME','PEDESTAL_INPUT_FILENAME']:
            if key not in kwd.keys():
                raise NotImplementedError("Relevant argument '{0}' not present!".format(key))
        # The input files
        self._pedestal_raw_file = kwd['PEDESTAL_INPUT_FILENAME']
        self._beam_raw_file     = kwd['ALIBAVA_INPUT_FILENAME']
        
        # remove all the lcio files except the last one...
        toremove = set([])
        # The bash file to concatenate all the process
        thebash = '#!/bin/bash\n\n'
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        thebash += '# Automaticaly created by open_sesame script at {0}\n'.format(st)
        thebash += '# ALIBAVA data chain reconstruction\n\n'.format(st)
        # Setting the values provided by the user
        for (step,args,action) in self.step_chain:
            # Redefine the args dict, activating the values of the dict
            newargs = dict(map(lambda (x,y): (x,y()), args.iteritems()))
            # Create the steering file for this step
            step.publish_steering_file(**newargs)
            # The particular action defined: it will update file names...
            # See the properties and updaters defined above
            action(step)
            thebash += 'echo "\033[1;34mRUNNING\033[1;m: \033[1;29mMarlin {0}\033[1;m"\n'.format(step.steering_file)
            thebash += 'time Marlin {0}\n'.format(step.steering_file)
            # remove the intermediate created files
            toremove.add(self.last_output_filename())
            try:
                toremove.add(self.pedestal_file())
            except AttributeError:
                # still not available
                pass
        thebash += "\nrm "
        for i in filter(lambda x: x not in [self.pedestal_file(),self.last_output_filename()],toremove):
            thebash += i+" "
        thebash+='\n\necho "ALiBaVa marlin data reconstruction chain DONE"\n'
        bashname = "alibava_full_reconstruction.sh"
        print "Created '{0}' script".format(bashname)
        # create the file
        with open(bashname,"w") as f:
            f.write(thebash)
        # get the mode
        bash_st = os.stat(bashname)
        os.chmod(bashname, bash_st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

# Telescope related
class telescope_conversion(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(telescope_conversion,self).__init__('telescope_conversion')

        self.steering_file_template = os.path.join(get_template_path(),'01-telescope_converter.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'TELESCOPE_INPUT_FILENAME', 'OUTPUT_FILENAME','GEAR_FILE',\
                "MAX_FIRING_FREQ_PIXEL")
        # Define a tuned default for the gear file, describes
        # telescope with no DUTs at all
        self.argument_values['GEAR_FILE']='gear_TB2017_CERNSPS_SETUP00_TELESCOPE_noDUTs.xml'
        # And copy the gear file to the relevant place
        self.auxiliary_files.append(self.argument_values['GEAR_FILE'])
    
    @staticmethod
    def get_description():
        return 'Convert Telescope ACONITE RAW binary data into LCIO, masking hot pixels'  

class telescope_clustering(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(telescope_clustering,self).__init__('telescope_clustering')

        self.steering_file_template = os.path.join(get_template_path(),'02-telescope_clustering.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAMES', 'OUTPUT_FILENAME','GEAR_FILE')
        # Define a tuned default for the gear file, describes
        # telescope with no DUTs at all
        self.argument_values['GEAR_FILE']='gear_TB2017_CERNSPS_SETUP00_TELESCOPE_noDUTs.xml'
        # And copy the gear file to the relevant place
        self.auxiliary_files.append(self.argument_values['GEAR_FILE'])
        self.auxiliary_files.append('histoinfo_telescope.xml')
    
    @staticmethod
    def get_description():
        return 'Find telescope cluster patterns, removing clusters with hot pixels'  

class telescope_filter(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(telescope_filter,self).__init__('telescope_filter')

        self.steering_file_template = os.path.join(get_template_path(),'03-telescope_filter.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAMES', 'OUTPUT_FILENAME','GEAR_FILE')
        # Define a tuned default for the gear file, describes
        # telescope with no DUTs at all
        self.argument_values['GEAR_FILE']='gear_TB2017_CERNSPS_SETUP00_TELESCOPE_noDUTs.xml'
        # And copy the gear file to the relevant place
        self.auxiliary_files.append(self.argument_values['GEAR_FILE'])
        self.auxiliary_files.append('histoinfo_telescope.xml')
    
    @staticmethod
    def get_description():
        return 'Filters the telescope clusters (very slow process!)'  

# Merge telescope and alibava data
class merger(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(merger,self).__init__('merger')

        self.steering_file_template = os.path.join(get_template_path(),'10-merger.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'TELESCOPE_INPUT_FILENAME',\
                'ALIBAVA_INPUT_FILENAME', 'OUTPUT_FILENAME','GEAR_FILE')
        # Define a tuned default for the gear file, describes
        # telescope with no DUTs at all
        #self.argument_values['GEAR_FILE']='gear_TB2017_CERNSPS_SETUP00_TELESCOPE_noDUTs.xml'
        # And copy the gear file to the relevant place
        #self.auxiliary_files.append(self.argument_values['GEAR_FILE'])
    
    @staticmethod
    def get_description():
        return 'Merge the alibava and telescope data'  

# ==================================================================================================
# The available marlin_steps classes (ordered)
available_steps = (pedestal_conversion,pedestal_preevaluation,cmmd_calculation,pedestal_evaluation,\
        calibration_conversion,calibration_extraction,\
        rs_conversion,signal_reconstruction,alibava_clustering,
        alibava_full_reco,
        # Telescope related
        telescope_conversion,telescope_clustering,telescope_filter,
        # Join both 
        merger
        )
# ==================================================================================================

# END -- Marlin step concrete implementations
# -------------------------------------------


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

