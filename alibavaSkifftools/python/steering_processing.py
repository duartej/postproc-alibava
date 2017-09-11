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

# The common name for the gear files 
def get_gear_filename_common():
    """The per default name of a gear file, centralized in order
    to be common along the whole package

    Return
    ------
    gear_filename: str
        a generic filename for the gear files, without extension
    """
    return 'gear_file'

# Available arguments present in the steering file templates.
# These arguments are dynamically changed on run time
_ARGUMENTS = { 'ROOT_FILENAME': 'Name of the output root file created by the AIDA processor',
        'CURRENT_WORKING_DIR': 'Working directory',
        'RUN_NUMBER': 'The run number of the input file',
        'ALIBAVA_INPUT_FILENAME': 'The input file name (ALIBAVA RAW data or slcio for the merger)',
        'ALIBAVA_REF_INPUT_FILENAME': 'The input file name for the reference ALIBAVA (the slcio for the merger)',
        'ACTIVE_CHIP': 'The beetle chip used, automaticaly defined given the run number and the sensor name',
        'ACTIVE_CHANNELS': 'The list of the active channels in ranges comma-separated (range edges are'\
                ' included as actives), for instance, A:B,C:D ...',
        'ENABLE_AUTOMASKING': 'This option activates the noisy channel auto masking, by masking'\
                ' those channels with a |noise_ch - <noise>| > [C] sigma_noise. See also '\
                '"CRITERIUM_AUTOMASKING" option',
        'CRITERIUM_AUTOMASKING': 'The [C]-value to be used when "ENABLE_AUTOMASKING" is ON [Default: 2.5]',
        'GEO_ID': 'The geometrical identificator, automaticaly defined given the run number and the sensor name',
        'TELESCOPE_INPUT_FILENAME': 'The input file name (ACONITE Telescope RAW data or slcio for the merger)',
        'INPUT_FILENAMES': 'The list of input file names (LCIO DATA)',
        'OUTPUT_FILENAME': 'Name of the output LCIO file created by the LCIOOutputProcessor',
        'GEAR_FILE': 'The name of the gear file to be used',
        'PEDESTAL_OUTPUT_FILENAME': 'Name of the output LCIO file created by the AlibavaPedestalNoiseProcessor',
        'PEDESTAL_INPUT_FILENAME': 'Name of the input LCIO file created by the AlibavaPedestalNoiseProcessor,'\
                ' containing the pedestals',
        'CALIBRATION_OUTPUT_FILENAME': 'Name of the output LCIO file created by the AlibavaCalibrateProcessor',
        'CALIBRATION_INPUT_FILENAME': 'Name of the input alibava file for calibration (Only for alibava_full_reco)',
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
        'PREALIGN_DUMP_GEAR': 'Whether to use or not a gear file to store the alignment constants, or use a LCIO file',
        'DUT_PLANES': 'The list of planes which need to find the missing coordinate',
        'MAX_RESIDUAL': 'The Maximum distance to determine if a hit is correlated [mm]',
        'REF_PLANE_LEFT': 'The telescope planes nearest to the DUT from the left, to extrapolate the hit',
        'REF_PLANE_RIGHT': 'The telescope planes nearest to the DUT from the right, to extrapolate the hit',
        'ITERATION': 'The number of iteration of the alignment job, mandatory argument of the alignment step',
        'ALIGNMENT_PROCESSOR_LOAD': 'Internal use: The placeholder to put the '\
                'processor load for the alignment',
        'ALIGNMENT_PROCESSOR_DESCRIPTION': 'Internal use: The placeholder to '\
                'put the processor description for the constant alignment load',
        'PREITERATION': 'Internal use: The previous number of iteration',
        'ALIGN_CTE_NAME': 'The name of the alignment constant, automaticaly set from the ITERATION argument',
        'ALIGN_CTE_LIST': 'Internal use: The list of the alignment constant to be used',
        'ALIGNED_HIT_LIST': 'Internal use: The list of hits after the alignment is applied',
        'DUMMY_REF_LIST': 'Internal use: Just a dummy list with the proper name of elements',
        'RESIDUAL_XMAX_U': 'The residual cut in X for the upstream planes (0,1)',
        'RESIDUAL_XMIN_U': 'The residual cut in X for the upstream planes (0,1)',
        'RESIDUAL_YMAX_U': 'The residual cut in Y for the upstream planes (0,1)',
        'RESIDUAL_YMIN_U': 'The residual cut in Y for the upstream planes (0,1)',
        'RESIDUAL_XMAX_D': 'The residual cut in X for the downstream planes (2,3,4)',
        'RESIDUAL_XMIN_D': 'The residual cut in X for the downstream planes (2,3,4)',
        'RESIDUAL_YMAX_D': 'The residual cut in Y for the downstream planes (2,3,4)',
        'RESIDUAL_YMIN_D': 'The residual cut in Y for the downstream planes (2,3,4)',
        'RESOLUTION_X_U':  'The telescope resolution (in X) in the upstream planes (0,1)',
        'RESOLUTION_X_D':  'The telescope resolution (in X) in the downstream planes (2,3,4)',
        'RESOLUTION_Y_U':  'The telescope resolution (in Y) in the upstream planes (0,1)',
        'RESOLUTION_Y_D':  'The telescope resolution (in Y) in the downstream planes (2,3,4)',
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
        The map between the argument template name and the concrete
        value useda
    auxiliary_files: list(str)
        Files needed to be copied from the get_template_path() to
        the current working directory
    DUT: FIXME DOC
    parser_instance: FIXME DOC
    devices: list(str)
        The list of devices (telescope or DUT) which this step 
        can be applied
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
        # The DUT if any, and if the SPS2017TB.filename_parser 
        # instance used to obtained
        self.DUT = None
        self.parser_instance=None
        # The device to be applied this step
        self.devices = []

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
        from .SPS2017TB_metadata import filename_parser

        self.steering_file_content = self.steering_file_content.replace("{0}{1}{0}".format(self.token,argument),str(value))
        self.argument_values[argument] = value
        # An special case: just create the parser instance if not there
        if argument == "ALIBAVA_INPUT_FILENAME" and self.parser_instance is None:
            fnp = filename_parser(value)
            self.parser_instance = fnp

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
        elif argument == 'CURRENT_WORKING_DIR':
            return os.getcwd()
        elif argument == 'ACTIVE_CHANNELS':
            # As the chip number is not available until after this function
            # is called, just put in front of the ranges the @ACTIVE_CHIP@
            # in order to be substituted afterwards.
            # Per default, all of them are active
            return "$@ACTIVE_CHIP@:0-127$"
        elif argument == 'ENABLE_AUTOMASKING':
            return 1
        elif argument == 'CRITERIUM_AUTOMASKING':
            return 2.5
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
                    self.parser_instance = fnp
                    if fnp.is_beam:
                        return fnp.run_number
            return -1
        elif argument == 'GEAR_FILE':
            # First copy the file to the cwd
            #self.auxiliary_files.append('dummy_gear.xml')
            return get_gear_filename_common()
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
                # The same name than the alibava data, but added the DATAMERGED
                return "{0}_beam.DATAMERGED.slcio".format("_".join(fnp.filename.split("_")[:-1]))
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
            return 0.01
        elif argument == 'DUT_PLANES':
            return 5
        elif argument == 'MAX_RESIDUAL':
            return 0.005
        elif argument == 'REF_PLANE_LEFT':
            return 1
        elif argument == 'REF_PLANE_RIGHT':
            return 2
        elif argument == 'PREALIGN_DUMP_GEAR':    
            return 'false'
               
        raise RuntimeError('Argument "{0}" must be explicitely set'.format(argument))

    def special_preprocessing(self,**kwd):
        """Performs any special preprocessing step needed
        in the 'publish_steering_file' function, just before
        of filling the arguments.
        Note this function must be implemented in the concrete steps, 
        except for the 'ACTIVE_CHANNELS' argument which is present 
        in several alibava classes, and therefore, useful to code it
        here. 
        In that case, the user could enter an string like: A:B,C:D,...
        where A B C D are the edges of the channel ranges considered
        actives (edges included). As the active chip is given by 
        the name of the sensor, the template ACTIVE_CHIP is set and
        be filled afterwards

        Parameters
        ----------
        kwd: dict
            the dictionary of arguments, which must be defined
            at _ARGUMENTS
        
        Return
        ------
        kwd: the updated dictionary
        
        Raises
        ------
        NotImplementedError
            If any of the introduced arguments is not defined in 
            the _ARGUMENTS static dictionary
        """
        if kwd.has_key('ACTIVE_CHANNELS'):
            active_channels = ''
            for channels in kwd['ACTIVE_CHANNELS'].split(","):
                ch_min,ch_max = channels.split(":")
                active_channels+= "$@ACTIVE_CHIP@:{0}-{1}$ ".format(ch_min,ch_max)
            # Re-write the channels with the proper formatted string
            kwd['ACTIVE_CHANNELS'] = active_channels
        return kwd

    def publish_steering_file(self,**kwd):
        """Creates a copy of the steering file with the particular
        values introduced substituted. If a given argument is not
        provided, the default value will be used (see the __init__
        methods of the concrete marlin_step classes definitions, 
        and the `self.get_argument_default` function)

        Parameters
        ----------
        kwd: dict
            the dictionary of arguments, which must be defined 
            at _ARGUMENTS            

        Raises
        ------
        NotImplementedError
            If any of the introduced arguments is not defined in 
            the _ARGUMENTS static dictionary
        """
        import shutil
        import os
        from .SPS2017TB_metadata import get_beetle
        from .SPS2017TB_metadata import get_gear_content
        from .SPS2017TB_metadata import get_geo_id
        from .SPS2017TB_metadata import get_standard_sensor_name as ssnm

        # Check for inconsistencies
        for key in kwd.keys():
            if key not in _ARGUMENTS.keys():
                raise NotImplementedError("The argument '{0}' is not implemented".format(key))
        
        # Pre-formatted special cases
        kwd = self.special_preprocessing(**kwd)
        
        # Setting the values provided by the user
        for arg,value in kwd.iteritems():
            self.set_argument_value(arg,value)
        
        # The default values (note that if the user already set it,
        # the next lines do not affect anything as the @KEY@ is not present
        # anymore in steering_file_content
        for argument in filter(lambda _a: _a not in self.argument_values.keys(),self.required_arguments):
            val = self.get_default_argument_value(argument)
            self.set_argument_value(argument,val)
        
        # Apply automaticaly the chip number and gear file content,
        # once the run number is available.
        # Check for the sensor to be process (if any)
        if self.parser_instance is not None:
            sensor_name = ssnm(self.parser_instance.sensor_name)
            self.DUT = get_beetle(sensor_name)
            # Obtain the proper GEAR file and create in the working dir
            # Also include the REFerence detector
            gear_content = get_gear_content(self.argument_values['RUN_NUMBER'],\
                    sensor_name=sensor_name,include_ref=True)
            geoid = get_geo_id(self.argument_values['RUN_NUMBER'],sensor_name)
            # XXX -- Be CAREFUL WITH this file name (to be centralized)
            gear_filename = 'gear_file'
            with open(gear_filename+".xml",'w') as f:
                f.write(gear_content)
            # More setters: 
            # First the active channels
            self.set_argument_value('ACTIVE_CHIP',self.DUT)
            self.set_argument_value('GEO_ID',geoid)
            # INCLUDE all this info int the gear_file 
            self.set_argument_value('GEAR_FILE',gear_filename)
        # The telescope case needs also to create the gear file
        if self.argument_values.has_key('TELESCOPE_INPUT_FILENAME') \
                and not self.argument_values.has_key('ALIBAVA_INPUT_FILENAME'):
            gear_content = get_gear_content(self.argument_values['RUN_NUMBER'])
            geoid = get_geo_id(self.argument_values['RUN_NUMBER'])
            # XXX -- Be CAREFUL WITH this file name (to be centralized)
            gear_filename = 'gear_file.xml'
            with open(gear_filename,'w') as f:
                f.write(gear_content)
        
        # Those arguments not set by the user but set by the class,
        # i.e. not present in the argument_values dictionary. This 
        # allows to define tuned defaults when defining the concrete
        # marlin_reco daughters classes inside the __init__ function
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

# Alibava Reconstruction 
# ======================
class pedestal_conversion(marlin_step):
    def __init__(self):
        import os
        super(pedestal_conversion,self).__init__('pedestal_conversion')
        self.devices = ['DUT']

        self.steering_file_template = os.path.join(get_template_path(),'01-ab_converter.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'ALIBAVA_INPUT_FILENAME',\
                'OUTPUT_FILENAME','GEAR_FILE','ACTIVE_CHIP','GEO_ID' )
        # Dummy definitions to let the class populate the needed
        # data member (parser_instance) when reading the file name
        self.argument_values['ACTIVE_CHIP'] = -1
        self.argument_values['GEO_ID'] = -1
    
    @staticmethod
    def get_description():
        return 'Convert RAW binary data into LCIO, pedestal runs'  

class pedestal_preevaluation(marlin_step):
    def __init__(self):
        import os
        super(pedestal_preevaluation,self).__init__('pedestal_preevaluation')
        self.devices = ['DUT']

        self.steering_file_template = os.path.join(get_template_path(),'02-ped_preevaluation.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', \
                'PEDESTAL_OUTPUT_FILENAME','GEAR_FILE', 'ACTIVE_CHANNELS')
    
    @staticmethod
    def get_description():
        return 'Estimate pedestal and noise from a gaussian distribution'

class cmmd_calculation(marlin_step):
    def __init__(self):
        import os
        super(cmmd_calculation,self).__init__('cmmd_calculation')
        self.devices = ['DUT']

        self.steering_file_template = os.path.join(get_template_path(),'03-ped_cmmd_calculation.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_INPUT_FILENAME',\
                'OUTPUT_FILENAME','MAXADC','MINADC','NBINS','MAXCMMDERR','MINCMMDERR','GEAR_FILE',\
                'ACTIVE_CHANNELS')
    
    @staticmethod
    def get_description():
        return 'Common mode noise calculation'

class pedestal_evaluation(marlin_step):
    def __init__(self):
        import os
        super(pedestal_evaluation,self).__init__('pedestal_evaluation')
        self.devices = ['DUT']

        self.steering_file_template = os.path.join(get_template_path(),'04-ped_evaluation.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_OUTPUT_FILENAME',\
                'MAXADC','MINADC','NBINS','GEAR_FILE','ACTIVE_CHANNELS')
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
        self.devices = ['DUT']
        # Change the step name
        self.step_name='calibration_conversion'
        self.steering_file = self.steering_file.replace('pedestal_conversion',self.step_name)
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER',\
                'ALIBAVA_INPUT_FILENAME','GEAR_FILE','OUTPUT_FILENAME','ACTIVE_CHIP','GEO_ID')
        # Dummy definitions to let the class populate the needed
        # data member (parser_instance) when reading the file name
        self.argument_values['ACTIVE_CHIP'] = -1
        self.argument_values['GEO_ID'] = -1
    
    @staticmethod
    def get_description():
        return 'Convert RAW binary data into LCIO, calibration runs'  

class calibration_extraction(marlin_step):
    def __init__(self):
        import os
        super(calibration_extraction,self).__init__('calibration_extraction')
        self.devices = ['DUT']

        self.steering_file_template = os.path.join(get_template_path(),'02-cal_extraction.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'CALIBRATION_OUTPUT_FILENAME',\
                'GEAR_FILE','ACTIVE_CHANNELS')
    
    @staticmethod
    def get_description():
        return 'Extract calibration constant per channel'

class rs_conversion(marlin_step):
    def __init__(self):
        import os
        super(rs_conversion,self).__init__('rs_conversion')
        self.devices = ['DUT']

        self.steering_file_template = os.path.join(get_template_path(),'01-ab_converter_rs.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'ALIBAVA_INPUT_FILENAME', \
                'ACTIVE_CHIP','GEO_ID','OUTPUT_FILENAME','GEAR_FILE',"TIMECUT_MIN","TIMECUT_MAX")
        # Dummy definitions to let the class populate the needed
        # data member (parser_instance) when reading the file name
        self.argument_values['ACTIVE_CHIP'] = -1
        self.argument_values['GEO_ID'] = -1
    
    @staticmethod
    def get_description():
        return 'Convert RAW binary data into LCIO, beam or RS runs'  

class signal_reconstruction(marlin_step):
    def __init__(self):
        import os
        super(signal_reconstruction,self).__init__('signal_reconstruction')
        self.devices = ['DUT']

        self.steering_file_template = os.path.join(get_template_path(),'02-signal_reconstruction.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_INPUT_FILENAME',\
                'MAXADC','MINADC','NBINS','MAXCMMDERR','MINCMMDERR','CMMDCUT_MIN','CMMDCUT_MAX','GEAR_FILE',\
                'OUTPUT_FILENAME','ACTIVE_CHANNELS',"ENABLE_AUTOMASKING","CRITERIUM_AUTOMASKING")
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
        self.devices = ['DUT']

        self.steering_file_template = os.path.join(get_template_path(),'03-ab_clustering.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_INPUT_FILENAME',\
                'SNRCUT_NGB','SNRCUT_SEED','SIGNAL_POLARITY','GEAR_FILE','SENSORID_STARTS','OUTPUT_FILENAME',\
                'ACTIVE_CHANNELS',"ENABLE_AUTOMASKING","CRITERIUM_AUTOMASKING")
    
    @staticmethod
    def get_description():
        return 'Cluster finding algorithm and conversion from Alibava clusters to EUTelSparseCluster'

class cluster_histograms(marlin_step):
    def __init__(self):
        import os
        super(cluster_histograms,self).__init__('cluster_histograms')
        self.devices = ['DUT']

        self.steering_file_template = os.path.join(get_template_path(),'04-cluster_histograms.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER','INPUT_FILENAMES', 'PEDESTAL_INPUT_FILENAME',\
                'CALIBRATION_INPUT_FILENAME','SIGNAL_POLARITY','GEAR_FILE','ACTIVE_CHANNELS',\
                'ENABLE_AUTOMASKING','CRITERIUM_AUTOMASKING','CURRENT_WORKING_DIR')
        # Needed files
        self.auxiliary_files.append('histoinfo_alibava.xml')
    
    @staticmethod
    def get_description():
        return 'Cluster plotter'
    
# Metaclass to deal with the full reconstruction for ALIBAVA
class alibava_full_reco(marlin_step):
    def __init__(self):
        import os
        super(alibava_full_reco,self).__init__('alibava_full_reco')
        self.devices = ['DUT']
        
        # -- Dummy 
        self.required_arguments = ()

        # The list of steps with their needed arguments
        self.step_chain = ( 
                (pedestal_conversion(), { 'ALIBAVA_INPUT_FILENAME': self.pedestal_raw_file},self.update_output),
                (pedestal_preevaluation(), {'INPUT_FILENAMES':self.last_output_filename,'ACTIVE_CHANNELS':self.active_channels},self.update_pedestal),
                (cmmd_calculation(), { 'INPUT_FILENAMES': self.last_output_filename, 'PEDESTAL_INPUT_FILENAME': self.pedestal_file, \
                        'ACTIVE_CHANNELS': self.active_channels },self.update_output),
                (pedestal_evaluation(),{'INPUT_FILENAMES': self.last_output_filename, 'ACTIVE_CHANNELS':self.active_channels},self.update_pedestal),
                (calibration_conversion(), {'ALIBAVA_INPUT_FILENAME': self.calibration_raw_file},self.update_output), 
                (calibration_extraction(), {'INPUT_FILENAMES': self.last_output_filename,'ACTIVE_CHANNELS':self.active_channels},self.update_calibration),
                (rs_conversion(), { 'ALIBAVA_INPUT_FILENAME': self.beam_raw_file},self.update_output),
                (signal_reconstruction(), {'INPUT_FILENAMES': self.last_output_filename, 'PEDESTAL_INPUT_FILENAME': self.pedestal_file,\
                        'ACTIVE_CHANNELS': self.active_channels },self.update_output),
                (alibava_clustering(), {'INPUT_FILENAMES': self.last_output_filename, 'PEDESTAL_INPUT_FILENAME': self.pedestal_file,\
                        'ACTIVE_CHANNELS': self.active_channels },self.update_output),
                # Note that the cluster_histograms
                (cluster_histograms(), { 'INPUT_FILENAMES': self.last_output_filename, \
                        'PEDESTAL_INPUT_FILENAME': self.pedestal_file, 'CALIBRATION_INPUT_FILENAME': self.calibration_file,
                        'ACTIVE_CHANNELS': self.active_channels}, self.dummy )
                )

    # Some datamembers used in the step_chain are not going to be 
    # populated until the call to the step is performed, using them
    # as properties. Note the use of the updaters (like the setters of
    # a property) methods (defined below, are the mechanism to update
    # this elements)
    def pedestal_raw_file(self):
        return self._pedestal_raw_file
    
    def calibration_raw_file(self):
        return self._calibration_raw_file
    
    def beam_raw_file(self):
        return self._beam_raw_file

    def last_output_filename(self):
        return self._last_output_filename
    
    def pedestal_file(self):
        return self._pedestal_file
    
    def calibration_file(self):
        return self._calibration_file

    def active_channels(self):
        return self._active_channels
    
    # Updaters
    def update_output(self,step_inst):
        """Update the last_output_filename data member using the value
        from the last step
        """
        self._last_output_filename = step_inst.argument_values['OUTPUT_FILENAME']

    def update_pedestal(self,step_inst):
        self._pedestal_file = step_inst.argument_values['PEDESTAL_OUTPUT_FILENAME']
    
    def update_calibration(self,step_inst):
        self._calibration_file = step_inst.argument_values['CALIBRATION_OUTPUT_FILENAME']

    def dummy(self,step_inst):
        pass

    @staticmethod
    def get_description():
        return '\033[1;32mMetaclass\033[1;m to perform the whole bunch of steps to process the ALiBaVa data'

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
        CALIBRATION_INPUT_FILENAME: str
            the name of the raw alibava data calibration file

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
        for key in ['ALIBAVA_INPUT_FILENAME','PEDESTAL_INPUT_FILENAME','CALIBRATION_INPUT_FILENAME']:
            if key not in kwd.keys():
                raise NotImplementedError("Relevant argument '{0}' not present!".format(key))
        # The input files
        self._pedestal_raw_file   = kwd['PEDESTAL_INPUT_FILENAME']
        self._calibration_raw_file= kwd['CALIBRATION_INPUT_FILENAME']
        self._beam_raw_file       = kwd['ALIBAVA_INPUT_FILENAME']

        # The active channels, if any
        try:
            self._active_channels = kwd['ACTIVE_CHANNELS']
        except KeyError:
            # use the default active channels
            from .SPS2017TB_metadata import active_channels,filename_parser,get_standard_sensor_name
            sn = get_standard_sensor_name(filename_parser(self._beam_raw_file).sensor_name)
            ac_list = active_channels[sn]
            self._active_channels = ''
            # Set the expected format 'ch1:ch2,ch3:ch4,...'
            for i in xrange(0,len(ac_list),2):
                self._active_channels += '{0}:{1},'.format(ac_list[i],ac_list[i+1])
            self._active_channels = self._active_channels[:-1]
        
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
        for i in filter(lambda x: x not in [self.pedestal_file(),self.calibration_file(),self.last_output_filename()],toremove):
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

# Telescope Reconstruction 
# ========================
class telescope_conversion(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(telescope_conversion,self).__init__('telescope_conversion')
        self.devices = ['Telescope']

        self.steering_file_template = os.path.join(get_template_path(),'01-telescope_converter.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'TELESCOPE_INPUT_FILENAME', 'OUTPUT_FILENAME','GEAR_FILE',\
                "MAX_FIRING_FREQ_PIXEL")
        # Define a tuned default for the gear file, describes
        # telescope with no DUTs at all
        #self.argument_values['GEAR_FILE']='gear_TB2017_CERNSPS_SETUP00_TELESCOPE_noDUTs.xml'
        # And copy the gear file to the relevant place
        #self.auxiliary_files.append(self.argument_values['GEAR_FILE'])
    
    @staticmethod
    def get_description():
        return 'Convert Telescope ACONITE RAW binary data into LCIO, masking hot pixels'  

class telescope_clustering(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(telescope_clustering,self).__init__('telescope_clustering')
        self.devices = ['Telescope']

        self.steering_file_template = os.path.join(get_template_path(),'02-telescope_clustering.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAMES', \
                'OUTPUT_FILENAME','GEAR_FILE','CURRENT_WORKING_DIR')
        # Define a tuned default for the gear file, describes
        # telescope with no DUTs at all
        #self.argument_values['GEAR_FILE']='gear_TB2017_CERNSPS_SETUP00_TELESCOPE_noDUTs.xml'
        # And copy the gear file to the relevant place
        #self.auxiliary_files.append(self.argument_values['GEAR_FILE'])
        self.auxiliary_files.append('histoinfo_telescope.xml')
    
    @staticmethod
    def get_description():
        return 'Find telescope cluster patterns, removing clusters with hot pixels'  

class telescope_filter(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(telescope_filter,self).__init__('telescope_filter')
        self.devices = ['Telescope']

        self.steering_file_template = os.path.join(get_template_path(),'03-telescope_filter.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAMES', \
                'OUTPUT_FILENAME','GEAR_FILE','CURRENT_WORKING_DIR')
        # Define a tuned default for the gear file, describes
        # telescope with no DUTs at all
        #self.argument_values['GEAR_FILE']='gear_TB2017_CERNSPS_SETUP00_TELESCOPE_noDUTs.xml'
        # And copy the gear file to the relevant place
        #self.auxiliary_files.append(self.argument_values['GEAR_FILE'])
        self.auxiliary_files.append('histoinfo_telescope.xml')
    
    @staticmethod
    def get_description():
        return 'Filters the telescope clusters (very slow process!)'  

class telescope_alignment(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(telescope_alignment,self).__init__('telescope_alignment')
        self.devices = ['Telescope']

        self.steering_file_template = os.path.join(get_template_path(),'04-telescope_alignment.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAMES', \
                'GEAR_FILE',\
                'ITERATION','ALIGN_CTE_NAME','ALIGN_CTE_LIST',\
                'ALIGNED_HIT_LIST','DUMMY_REF_LIST',\
                'RESIDUAL_XMAX_U','RESIDUAL_XMIN_U','RESIDUAL_XMAX_D','RESIDUAL_XMIN_D',
                'RESIDUAL_YMAX_U','RESIDUAL_YMIN_U','RESIDUAL_YMAX_D','RESIDUAL_YMIN_D',
                'RESOLUTION_X_U','RESOLUTION_X_D','RESOLUTION_Y_U','RESOLUTION_Y_D',
                'ALIGNMENT_PROCESSOR_LOAD','ALIGNMENT_PROCESSOR_DESCRIPTION')
        # Define iteration-dependent cuts
        self._dependent_cuts = { 
                'RESIDUAL_XMAX_U': { 0: 300, 1: 100, 2: 50, 3: 50, 4: 25},
                'RESIDUAL_XMIN_U': { 0: -300, 1: -100, 2: -50, 3: -50, 4: -25},
                'RESIDUAL_XMAX_D': { 0: 500, 1: 200, 2: 100, 3: 50, 4: 25},
                'RESIDUAL_XMIN_D': { 0: -500, 1: -200, 2: -100, 3: -50, 4: -25},
                'RESIDUAL_YMAX_U': { 0: 300, 1: 100, 2: 50, 3: 50, 4: 25},
                'RESIDUAL_YMIN_U': { 0: -300, 1: -100, 2: -50, 3: -50, 4: -25},
                'RESIDUAL_YMAX_D': { 0: 500, 1: 200, 2: 100, 3: 50, 4: 25},
                'RESIDUAL_YMIN_D': { 0: -500, 1: -200, 2: -100, 3: -50, 4: -25},
                'RESOLUTION_X_U': { 0: 10, 1: 8, 2: 5, 3: 3, 4: 3},
                'RESOLUTION_X_D': { 0: 10, 1: 8, 2: 5, 3: 3, 4: 3},
                'RESOLUTION_Y_U': { 0: 10, 1: 8, 2: 5, 3: 3, 4: 3},
                'RESOLUTION_Y_D': { 0: 10, 1: 8, 2: 5, 3: 3, 4: 3},
                }
    
    @staticmethod
    def get_description():
        return 'Telescope alignment iteratively'  
    
    def special_preprocessing(self,**kwd):
        """Concrete implementation of the virtual function.
        Set the alignment constant name, the residuals and
        resolutions corresponding to the given iteration.

        Parameters
        ----------
        kwd: dict
            the dictionary of arguments, which must be defined
            at _ARGUMENTS. Must contain
             - ITERATION: the alignment iteration
            Optionally, could contain (if not, default values
            are going to be filled):
             - RESIDUAL_ZZZZ_Z: the (max/min) diference between hit measured
                   and hit predicted (straight line) in the (X/Y) axis to 
                   be accepted as valid hit belonging to a track
             - RESOLUTION_Z_Z: the telescope resolution in the X/Y axis

        Return
        ------
        kwd: the updated dictionary
        
        Raises
        ------
        TypeError
            If the 'ITERATION' is not an integer
        
        RuntimeError
            If ITERATION is not present

        NotImplementedError
            If any of the introduced arguments is not defined in 
            the _ARGUMENTS static dictionary
        """
        import os

        #   - ITERATION
        if not kwd.has_key('ITERATION'):
            raise RuntimeError("Needed 'ITERATION' argument not provided")
        if int(kwd['ITERATION']) == 0:
            kwd['PREITERATION'] = 'pre'
        elif kwd['ITERATION'].isdigit():
            kwd['PREITERATION'] = int(kwd['ITERATION'])-1
        else:
            raise TypeError("Not valid 'ITERATION' argument passed: "+kwd['ITERATION'])
        kwd['ALIGN_CTE_NAME'] = 'alignment{0}'.format(kwd['ITERATION'])

        # The run number is already needed, so the input filename as well
        # in order to find the run number -- CAREFUL, breaking the logic flow
        # contained in the 'publish_steering_file' function
        self.set_argument_value('INPUT_FILENAMES',kwd['INPUT_FILENAMES'])
        kwd['RUN_NUMBER'] = self.get_default_argument_value('RUN_NUMBER')

        # Prepare a processor for each alignment already performed (the prealigment is the 
        # minimum performed)
        kwd['ALIGNMENT_PROCESSOR_LOAD']='<processor name="LoadPreAlignment"/>'
        kwd['ALIGNMENT_PROCESSOR_DESCRIPTION']='<processor name="LoadPreAlignment"'\
                'type="ConditionsProcessor">\n        <parameter name="DBInit" type="string" value="localhost:lccd_test:align:tel"/>\n'\
                '        <parameter name="SimpleFileHandler" type="StringVec">'\
                'prealign {0}/{1}-telescopeOnly-prealign-db.slcio alignment </parameter>'\
                '\n    </processor>'.format(os.getcwd(),kwd['RUN_NUMBER'])
        kwd['ALIGN_CTE_LIST'] = 'prealign'
        kwd['ALIGNED_HIT_LIST'] = ''
        kwd['DUMMY_REF_LIST'] = 'dummy0'
        for i in xrange(1,int(kwd['ITERATION'])+1):
            kwd['ALIGNMENT_PROCESSOR_LOAD']+='\n\t<processor name="Load{0}Alignment"/>'.format(i-1)
            kwd['ALIGNMENT_PROCESSOR_DESCRIPTION']+='\n    <processor name="Load{0}Alignment"'\
                    'type="ConditionsProcessor">\n        <parameter name="DBInit" type="string" value="localhost:lccd_test:align:tel"/>\n'\
                    '        <parameter name="SimpleFileHandler" type="StringVec">'\
                    'alignment{0} {1}/{2}-telescopeOnly-{0}align-db.slcio alignment{0}</parameter>'\
                    '\n    </processor>'.format(i-1,os.getcwd(),kwd['RUN_NUMBER'])
            kwd['ALIGN_CTE_LIST']   = 'alignment{0} {1}'.format(i-1,kwd['ALIGN_CTE_LIST'])
            kwd['ALIGNED_HIT_LIST'] = 'alignedHit_{0} {1}'.format(i-1,kwd['ALIGNED_HIT_LIST'])
            kwd['DUMMY_REF_LIST'] = 'dummy{0} {1}'.format(i,kwd['DUMMY_REF_LIST'])

        # Set the per default values if the user did not provide them
        # Note the dependenced of the iteration, therefore cannot 
        # be set in the __init__
        args_to_set = ['RESIDUAL_XMAX_U','RESIDUAL_XMIN_U',\
                'RESIDUAL_XMAX_D','RESIDUAL_XMIN_D',
                'RESIDUAL_YMAX_U','RESIDUAL_YMIN_U',\
                'RESIDUAL_YMAX_D','RESIDUAL_YMIN_D',\
                'RESOLUTION_X_U','RESOLUTION_X_D',\
                'RESOLUTION_Y_U','RESOLUTION_Y_D']
        for thearg in args_to_set:
            if not kwd.has_key(thearg):
                k = int(kwd['ITERATION'])
                # -- Find the last ITERATION value in the 
                #    _dependent_cuts dict. Higher iterations
                #    are equivalent than the last available one 
                while k >= 0:
                    try:
                        kwd[thearg] = self._dependent_cuts[thearg][k]
                        break
                    except KeyError:
                        k-=1
        return kwd

class telescope_update_gear(marlin_step):
    def __init__(self,in_full_reco_mode=False):
        import os
        import shutil
        super(telescope_update_gear,self).__init__('telescope_update_gear')
        self.devices = ['Telescope']

        self.steering_file_template = os.path.join(get_template_path(),'041-telescope_update_gear.xml')
        self.required_arguments = ('GEAR_FILE','ALIGN_CTE_LIST',\
                'ALIGNMENT_PROCESSOR_LOAD','ALIGNMENT_PROCESSOR_DESCRIPTION')
        # Auxiliary member to take into account if the step is run
        # within the telescope_full_reco metaclass, in that case 
        # the extraction of the iteration in the special_preprocesing
        # method is performed in a different approach
        self._embebbed_metaclass = in_full_reco_mode
        
        # A dummy lcio with one event to enter in the processEvent
        self.auxiliary_files.append('dummy_lcio.slcio')
    
    @staticmethod
    def get_description():
        return 'Dump the alignment constants to a gear file'  
    
    def special_preprocessing(self,**kwd):
        """Concrete implementation of the virtual function.
        Set the alignment constant name, the residuals and
        resolutions corresponding to the given iteration.

        Parameters
        ----------
        kwd: dict
            the dictionary of arguments, which must be defined
            at _ARGUMENTS. Must contain
             - ITERATION: the alignment iteration

        Return
        ------
        kwd: the updated dictionary
        
        Raises
        ------
        TypeError
            If the 'ITERATION' is not an integer
        
        RuntimeError
            If ITERATION is not present

        IOError
            If not found any of the slcio alignment DB files
        """
        import os
        import glob

        # The run number is already needed, so we need to obtain it from the 
        # calculated align db present in the working directory
        if not self._embebbed_metaclass:
            # WARNING: this is assuming the telescope_alignment step has been
            #          performed previously (no sense otherwise) and we are
            #          running in the working directory where it was run the
            #          telescope_alingment (this is more a defect of this code)
            #          Providing  a directory by the user can fix this issue
            dbfiles = glob.glob(os.path.join(os.getcwd(),"*-telescopeOnly-*align-db.slcio"))
            try:
                aligdb_file = os.path.basename(dbfiles[0])
            except IndexError:
                raise IOError("Not found in the current working directory the alignment DB files."\
                        " The 'telescope_alignment' step must be performed previously. "\
                        " If it was, then lauch the script in the same directory were"\
                        " the 'telescope_alignment' was run")
                # Very specific format
                run_number = int(aligdb_file.split("-")[0])
                # Extract the maximum number of iterations available 
                iter_max = len(dbfiles)
        else:
            # Just using the steering files created by the metaclass for the telescope_alignment
            # steps to obtain the run_number and the number of iterations. Also it is needed to
            # be added plus 1 in order to take into account the prealignment (not present 
            # as part of the telescope_alignment steps but inside the filter)
            iter_max = len(glob.glob(os.path.join(os.getcwd(),"telescope_alignment*xml")))+1
            run_number = kwd['RUN_NUMBER']


        # Prepare a processor for each alignment already performed (the prealigment is the 
        # minimum performed)
        kwd['ALIGNMENT_PROCESSOR_LOAD']='<processor name="LoadPreAlignment"/>'
        kwd['ALIGNMENT_PROCESSOR_DESCRIPTION']='<processor name="LoadPreAlignment"'\
                'type="ConditionsProcessor">\n        <parameter name="DBInit" type="string" value="localhost:lccd_test:align:tel"/>\n'\
                '        <parameter name="SimpleFileHandler" type="StringVec">'\
                'prealign {0}/{1}-telescopeOnly-prealign-db.slcio alignment </parameter>'\
                '\n    </processor>'.format(os.getcwd(),run_number)
        kwd['ALIGN_CTE_LIST'] = 'prealign'
        for i in xrange(1,iter_max):
            kwd['ALIGNMENT_PROCESSOR_LOAD']+='\n\t<processor name="Load{0}Alignment"/>'.format(i-1)
            kwd['ALIGNMENT_PROCESSOR_DESCRIPTION']+='\n    <processor name="Load{0}Alignment"'\
                    'type="ConditionsProcessor">\n        <parameter name="DBInit" type="string" value="localhost:lccd_test:align:tel"/>\n'\
                    '        <parameter name="SimpleFileHandler" type="StringVec">'\
                    'alignment{0} {1}/{2}-telescopeOnly-{0}align-db.slcio alignment{0}</parameter>'\
                    '\n    </processor>'.format(i-1,os.getcwd(),run_number)
            kwd['ALIGN_CTE_LIST']   = '{0} alignment{1}'.format(kwd['ALIGN_CTE_LIST'],i-1)

        return kwd

class telescope_fitter(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(telescope_fitter,self).__init__('telescope_fitter')
        self.devices = ['Telescope']

        self.steering_file_template = os.path.join(get_template_path(),'05-telescope_fitter.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAMES', \
                'OUTPUT_FILENAME','GEAR_FILE')
        # The new @GEAR_FILE@_aligned.xml is needed to be copied in the
        # working directory --> NOT ACTUALLY As the file will be created after 
        # the alignement process... 
        # self.auxiliary_files.append('@GEAR_FILE@_aligned.xml')
    
    @staticmethod
    def get_description():
        return 'Track fitter using a Deterministic annealing filter (DAF)'  
    
    #def special_preprocessing(self,**kwd):
    #    """Concrete implementation of the virtual function.
    #    Just fix the name of the GEAR file in order to be able to
    #    do the copy. The dictionary arguments is not modified, just
    #    the auxiliary_files content @GEAR_FILE@

    #    Parameters
    #    ----------
    #    kwd: dict
    #        the dictionary of arguments, which must be defined
    #        at _ARGUMENTS.

    #    Return
    #    ------
    #    kwd: the dictionary
    #    
    #    """
    #    gear_filename=''
    #    if kwd.has_key('GEAR_FILE'):
    #        gear_filename = kwd['GEAR_FILE']
    #    else:
    #        gear_filename = get_gear_filename_common()
    #    # Do the thing
    #    index = self.auxiliary_files.index('@GEAR_FILE@_aligned.xml')
    #    self.auxiliary_files[index] = self.auxiliary_files[index].replace('@GEAR_FILE@',gear_filename)

    #    return kwd
        
# Metaclass to deal with the full reconstruction for ALIBAVA
class telescope_full_reco(marlin_step):
    """Class to gather a set of marlin_step instances to run serialized, 
    where the inputs and outputs are related between them. 
    The 'step_chain' data-member (A,B,C) is a 3-tuple where A is the 
    marlin_step instance, B is a dictionary with the arguments and values
    to be used by the instance, and C is the needed functor to be called 
    when use the 3-tuple (at the 'publish_steering_file' method)
    
    """
    def __init__(self):
        import os
        super(telescope_full_reco,self).__init__('telescope_full_reco')
        self.devices = ['Telescope']
        
        # -- Dummy 
        self.required_arguments = () 

        # Iteration for the alignment
        self._iteration = 0
        self._iter_init = False
        
        # XXX FIXME: Include the alignment steps
        # The list of steps with their needed arguments
        self.step_chain = ( 
                (telescope_conversion(), { 'TELESCOPE_INPUT_FILENAME': self.raw_file},self.update_output),
                (telescope_clustering(), {'INPUT_FILENAMES':self.last_output_filename},self.update_output),
                (telescope_filter(), {'INPUT_FILENAMES': self.last_output_filename},self.update_output),
                (telescope_alignment(),{'INPUT_FILENAMES': self.last_output_filename, 'ITERATION': self.iteration},self.dummy),
                (telescope_alignment(),{'INPUT_FILENAMES': self.last_output_filename, 'ITERATION': self.iteration},self.dummy),
                (telescope_alignment(),{'INPUT_FILENAMES': self.last_output_filename, 'ITERATION': self.iteration},self.dummy),
                (telescope_alignment(),{'INPUT_FILENAMES': self.last_output_filename, 'ITERATION': self.iteration},self.dummy),
                (telescope_alignment(),{'INPUT_FILENAMES': self.last_output_filename, 'ITERATION': self.iteration},self.dummy),
                (telescope_update_gear(True),{'RUN_NUMBER': self.get_run_number},self.dummy),
                # Remove that, because the fitter does not store events, therefore the ALIBAVA merging 
                # MUST be done before the fitter
                # (telescope_fitter(),{'INPUT_FILENAMES': self.last_output_filename},self.update_output)
                )

    # Some datamembers used in the step_chain are not going to be 
    # populated until the call to the step is performed, using them
    # as properties. Note the use of the updaters (like the setters of
    # a property) methods (defined below, are the mechanism to update
    # this elements)
    def raw_file(self):
        return self._raw_file

    def last_output_filename(self):
        return self._last_output_filename
    
    def get_iteration_without_modification(self):
        if self._iter_init:
            return self._iteration-1
        else:
            return self._iteration
    
    # Updaters
    def update_output(self,step_inst):
        """Update the last_output_filename data member using the value
        from the last step
        """
        self._last_output_filename = step_inst.argument_values['OUTPUT_FILENAME']

    def iteration(self):
        """Return and update iteration
        """
        self._iter_init = True
        self._iteration += 1
        # Need it as string for posterior processing
        return "{0}".format(self._iteration-1)

    def get_run_number(self):
        """Obtain the run number from the current output filename (which should be
        the corresponding to the telescope_filter step)
        """
        import os

        if os.path.basename(self._last_output_filename).find('run000') != 0:
            raise RuntimeError("The file created by the 'telescope_filter'"\
                    " step does not follow the name convection: '{0}'".format(self._last_output_filename))
        rawnumber = os.path.basename(self._last_output_filename).replace('run000','').split('.')[0]
        return int(rawnumber)

    def dummy(self,step_inst):
        pass

    @staticmethod
    def get_description():
        return '\033[1;32mMetaclass\033[1;m to perform the whole bunch of steps to process the Telescope data'

    def publish_steering_file(self,**kwd):
        """Creates every steering file needed for reconstruct
        the TELESCOPE data using only the RAW telescope data filenames. 
        Note that the other values should be change manually later in 
        the created steering files

        Parameters
        ----------
        TELESCOPE_INPUT_FILENAME: str
            the name of the raw alibava data beam file

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
        import shutil
    
        # Check for inconsistencies
        for key in ['TELESCOPE_INPUT_FILENAME']:
            if key not in kwd.keys():
                raise NotImplementedError("Relevant argument '{0}' not present!".format(key))
        # The input files, initialization
        self._raw_file       = kwd['TELESCOPE_INPUT_FILENAME']

        # remove all the lcio files except the last one...
        toremove = set([])
        # The bash file to concatenate all the process
        thebash = '#!/bin/bash\n\n'
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        thebash += '# Automaticaly created by open_sesame script at {0}\n'.format(st)
        thebash += '# TELESCOPE data chain reconstruction\n\n'.format(st)
        # Setting the values provided by the user
        for (step,args,action) in self.step_chain:
            # Redefine the args dict, activating the values of the dict
            newargs = dict(map(lambda (x,y): (x,y()), args.iteritems()))
            # Create the steering file for this step
            step.publish_steering_file(**newargs)
            # In case of alignment, change the file name to include the alignment
            # iteration
            if isinstance(step,telescope_alignment):
                old_steering_filename = step.steering_file
                step.steering_file = old_steering_filename.replace(".xml","_{0}.xml".format(self.get_iteration_without_modification()))
                shutil.move(old_steering_filename,step.steering_file)
            # The particular action defined: it will update file names...
            # See the properties and updaters defined above
            action(step)
            thebash += 'echo "\033[1;34mRUNNING\033[1;m: \033[1;29mMarlin {0}\033[1;m"\n'.format(step.steering_file)
            thebash += 'time Marlin {0}\n'.format(step.steering_file)
            # remove the intermediate created files
            toremove.add(self.last_output_filename())
        thebash += "\nrm "
        for i in filter(lambda x: x not in [self.last_output_filename()],toremove):
            thebash += i+" "
        thebash+='\n\necho "TELESCOPE marlin data reconstruction chain DONE"\n'
        thebash+='\n\necho "\033[1;31mSOME STEPS ARE VERY DEMANDING, CONSIDER TO SEND THEM THE CLUSTER\033[1;m"\n'
        bashname = "telescope_full_reconstruction.sh"
        print "Created '{0}' script".format(bashname)
        # create the file
        with open(bashname,"w") as f:
            f.write(thebash)
        # get the mode
        bash_st = os.stat(bashname)
        os.chmod(bashname, bash_st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

# Merge telescope and alibava data
# ================================
class merger(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(merger,self).__init__('merger')
        self.devices = ['Telescope','DUT']

        self.steering_file_template = os.path.join(get_template_path(),'10-merger.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'TELESCOPE_INPUT_FILENAME',\
                'ALIBAVA_INPUT_FILENAME', 'ALIBAVA_REF_INPUT_FILENAME', 'OUTPUT_FILENAME','GEAR_FILE')
        # Define a tuned default for the gear file, describes
        # telescope with no DUTs at all
        #self.argument_values['GEAR_FILE']='gear_TB2017_CERNSPS_SETUP00_TELESCOPE_noDUTs.xml'
        # And copy the gear file to the relevant place
        #self.auxiliary_files.append(self.argument_values['GEAR_FILE'])
    
    @staticmethod
    def get_description():
        return 'Merge the alibava and telescope data'  

class hitmaker(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(hitmaker,self).__init__('hitmaker')
        self.devices = ['Telescope','DUT']

        self.steering_file_template = os.path.join(get_template_path(),'11-hitmaker.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAMES',\
                 'OUTPUT_FILENAME','GEAR_FILE')
    
    @staticmethod
    def get_description():
        return 'Hit local position'  

class prealignment(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(prealignment,self).__init__('prealignment')
        self.devices = ['Telescope','DUT']

        self.steering_file_template = os.path.join(get_template_path(),'12-prealignment.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAMES',\
                 'OUTPUT_FILENAME','GEAR_FILE','PREALIGN_DUMP_GEAR')
    
    @staticmethod
    def get_description():
        return 'Pre-alignment using the distance of the hits to the first telescope plane '  
        
    def special_preprocessing(self,**kwd):
        """Concrete implementation of the virtual function.
        Just change the boolean content to the equivalent string

        Parameters
        ----------
        kwd: dict
            the dictionary of arguments, which must be defined
            at _ARGUMENTS. Must contain
             - PREALIGN_DUMP_GEAR

        Return
        ------
        kwd: the updated dictionary
        
        Raises
        ------
        RuntimeError
            If PRE_ALIGNED_DUMP_GEAR is not present

        NotImplementedError
            If any of the introduced arguments is not defined in 
            the _ARGUMENTS static dictionary
        """
        if not kwd.has_key('PREALIGN_DUMP_GEAR'):
            raise RuntimeError("Needed 'PREALIGNED_DUMP_GEAR', argument not present")
        if kwd['PREALIGN_DUMP_GEAR'] == True:
            kwd['PREALIGN_DUMP_GEAR'] ='true'
        else:
            kwd['PREALIGN_DUMP_GEAR'] ='false'

        return kwd

class simple_coordinate_finder_DUT(marlin_step):
    def __init__(self):
        import os
        import shutil
        super(simple_coordinate_finder_DUT,self).__init__('simple_coordinate_finder_DUT')
        self.devices = ['Telescope','DUT']

        self.steering_file_template = os.path.join(get_template_path(),'121-simple_coordinate_finder_DUT.xml')
        self.required_arguments = ('ROOT_FILENAME','RUN_NUMBER', 'INPUT_FILENAMES',\
                 'OUTPUT_FILENAME','GEAR_FILE', 'DUT_PLANES', 'MAX_RESIDUAL',
                 'REF_PLANE_LEFT','REF_PLANE_RIGHT','CURRENT_WORKING_DIR')
        # NOTE: Defined a tuned default for the gear file, provided by the hitmaker
        # steps (in the PreAligner, see 11-hitmaker.xml). 
        # Copy the histogram file
        self.auxiliary_files.append('histoinfo_alibava.xml')
    
    @staticmethod
    def get_description():
        return 'Simple coordinate finder for the DUTs: a line extrapolation from Telescope planes'  

# ==================================================================================================
# The available marlin_steps classes (ordered)
available_steps = (pedestal_conversion,pedestal_preevaluation,cmmd_calculation,pedestal_evaluation,\
        calibration_conversion,calibration_extraction,\
        rs_conversion,signal_reconstruction,alibava_clustering,
        cluster_histograms,
        alibava_full_reco,
        # Telescope related
        telescope_conversion,telescope_clustering,telescope_filter,telescope_alignment, \
                telescope_update_gear,telescope_fitter,
        telescope_full_reco,
        # Join both 
        merger, hitmaker,prealignment,
        simple_coordinate_finder_DUT
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

