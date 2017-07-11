#!/usr/bin/env python
"""Several dictionaries with static information from the SPS
Test Beam performed at May, 2017 at CERN 

Information extracted partially from 
https://docs.google.com/spreadsheets/d/1Z4nlyHUdAhCy-oNC-472c0Sydg54GraveDxROsnZKrM/edit#gid=0
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

__all__ = [
        "sensor_names",
        "filename_parser",
        "associated_filenames",
        "active_sensors"
        ]

# the EOS path where to find the alibava data
eospath="/eos/user/d/duarte/alibava_data"

# the list of used sensors
sensor_names = [ 'LGAD7859W1H6_0_b1', 'iLGAD8533W1K05T_0_b2', 'REF_0_b1', 'M1-5_0_b2', 
        'M1-8_7e15_b2', 'M2-3_1e16_b2', 'N1-3_0_b1', 'N1-7_7e15_b2', 'N1-8_1e16_b1' ]

# A simple parser class for the filenames
class filename_parser(object):
    """A simple parser for the file names. 
    Check the list of fields at the __init__ func
    """
    def __init__(self,filename):
        """XXX: DOC
        """
        import os

        self.filename = filename
        self.is_beam = False
        self.is_calibration = False
        self.is_pedestal = False

        naming_len=len(os.path.basename(filename).split("_"))
        if naming_len == 10:
            self.date,self.hour,self.pcname,self.motherboard,self.sensor_name,\
                self.voltage_bias,self.current_leak,self.temperature,self.latency,\
                self.run_type = os.path.basename(filename).split("_")
        elif naming_len == 11:
            self.run_number,self.date,self.hour,self.pcname,self.motherboard,self.sensor_name,\
                self.voltage_bias,self.current_leak,self.temperature,self.latency,\
                self.run_type = os.path.basename(filename).split("_")
        else:
            raise NameError("Invalid naming convection for the file:"\
                    "'{0}'".format(os.path.basename(filename)))

    def __str__(self):
        if self.is_beam:
            runnumber_str = "[RUN:        {0}] ".format(self.run_number)
        elif self.is_pedestal:
            runnumber_str = "[PEDESTAL    RUN] "
        elif self.is_calibration:
            runnumber_str = "[CALIBRATION RUN] "
        message = "{0} Sensor:{1}, bias voltage:{2:0.1f} V, "\
            "leak current:{3:0.1f} uA".format(runnumber_str,\
            self.sensor_name,self.voltage_bias,self.current_leak)
        return message
    
    @property
    def hour(self):
        return self._hour
    @hour.setter
    def hour(self,hour):
        # Convert into seconds (from 00:00)
        h,m=hour.split("-")
        self._hour = (float(h)*3600.+float(m)*60.)

    @property
    def run_number(self):
        return self._run_number
    @run_number.setter
    def run_number(self,run_number):
        if "-" in run_number or 'RunN' in run_number:
            raise RuntimeError("Invalid format for the run number. This"\
                " file should be masked")
        self._run_number = int(run_number)
    
    @property
    def run_type(self):
        return self._run_type
    @run_type.setter
    def run_type(self,run_type):
        self._run_type = run_type.split(".")[0]
        if self.run_type == "beam":
            self.is_beam = True
        elif self.run_type == "cal":
            self.is_calibration = True
        elif self.run_type == "ped":
            self.is_pedestal = True
        else:
            raise RuntimeError("Not a valid run mode: {0}".format(run_type))
    
    @property
    def voltage_bias(self):
        return self._voltage_bias
    @voltage_bias.setter
    def voltage_bias(self,voltage):
        try:
            self._voltage_bias = float(voltage)
        except ValueError:
            self._voltage_bias = float(voltage.upper().replace("V",""))

    @property
    def current_leak(self):
        return self._current_leak
    @current_leak.setter
    def current_leak(self,current):
        try:
            self._current_leak = float(current)
        except ValueError:
            self._current_leak = float(current.upper().replace("D",".").replace("UA",""))
        
    @property
    def temperature(self):
        return self._temperature
    @temperature.setter
    def temperature(self,temperature):
        try:
            self._temperature = float(temperature)
        except ValueError:
            self._temperature = float(temperature.upper().replace("C",""))

    @property
    def latency(self):
        return self._latency
    @latency.setter
    def latency(self,latency):
        try:
            self._latency = int(latency)
        except ValueError:
            self._latency = int(latency.upper().replace("LAT",""))

    def get_epoch_time(self):
        """Return the date and hour in epoch time

        Returns
        -------
        epoch time
        """
        try:
            return self._epoch_time
        except AttributeError:
            import time
            pattern='%Y-%m-%d'
            self._epoch_time = int(time.mktime(time.strptime(self.date,pattern)))+self.hour
            return self._epoch_time

    def __eq__(self,other):
        """Pseudo equality to match the beam runs with the pedestal ones
        """
        if self.sensor_name != other.sensor_name or \
                self.motherboard != other.motherboard:
            return False
        # And also not so far in time (5 hours)
        # First convert to an Epoch time to fairly compare
        our_hour  = self.get_epoch_time()
        other_hour= other.get_epoch_time() 
        if abs(our_hour-other_hour) > 5.0*3600.0:
            return False
        # The propose float equalizers :
        equalizers = [ ('voltage_bias',1e-19), ('current_leak',10.0), 
                ('temperature',1e-19) ]
        for (equalizer,max_diff) in equalizers:
            if abs(getattr(other,equalizer)-getattr(self,equalizer)) > max_diff:
                return False
        # we are here, then everything is equal
        return True

    def closest(self,list_of_others):
        """Returning the instance (pedestal and/or calibration instances)
        closest to this beam instance.

        Note
        ----
        This functions is only useful (and make sense) for a beam instance
        and a list of pedestal/calibration instances

        Parameters
        ----------
        list_of_others: list(filename_parser instances)
            the list of pedestal/calibration files

        Returns
        -------
        the filename_parser closest to this beam instance
        """
        closest=list_of_others[0]
        for other in list_of_others[1:]:
            # Just the time taken as the nearby handle
            delta_hour_other   = self.get_epoch_time()-other.get_epoch_time()
            delta_hour_closest = self.get_epoch_time()-closest.get_epoch_time()
            # A positive distance is preferable (the beam run has been taken after
            # the calibration/pedestal run) 
            if delta_hour_other > 0 and delta_hour_closest < 0:
                closest = other
                continue
            if abs(delta_hour_other) < abs(delta_hour_closest):
                closest = other
        return closest

# A class to associate a beam file with the pedestal and calibration
class associated_filenames(object):
    """Create a 3-tuple of filename_parsers instances which links each
    beam run with its corresponding pedestal and calibration run
    """
    def __init__(self,fn_instance,list_available_pedcal):
        """XXX: DOC. MISSING
        """
        # Matching the beam file with the pedestal and calibrate
        pre_associated_files = filter(lambda _fcp: _fcp == fn_instance and \
                (_fcp.is_calibration or _fcp.is_pedestal), list_available_pedcal)
        # Just cleaning a bit if there is too close 
        if len(pre_associated_files) > 2:
            pedfile = fn_instance.closest(filter(lambda x: x.is_pedestal,pre_associated_files))
            calfile = fn_instance.closest(filter(lambda x: x.is_calibration,pre_associated_files))
            #pedfile = sorted(filter(lambda x: x.is_pedestal,pre_associated_files),key=lambda x: x.get_epoch_time())[0]
            #calfile = sorted(filter(lambda x: x.is_calibration,pre_associated_files),key=lambda x: x.get_epoch_time())[0]
            associated_files = [pedfile,calfile]
        else:
            associated_files = pre_associated_files

        if len(associated_files) == 0:
            raise IOError("No pedestal nor calibration runs found for the "\
                    " current 'filename_parser' instance:\n {0}".format(fn_instance))
        elif len(associated_files) < 2:
            # WARNING maybe?
            raise RuntimeError("Just found 1 non-beam run type for the"\
                    " current 'filename_parser' instance:\n{0}\n{1}".format(fn_instance,associated_files[0]))
        elif len(associated_files) > 2:
            # Never it could be here!!!
            message = ""
            for af in associated_files:
                message += " - "+af.filename+"\n"
            raise RuntimeError("Oversize number of ped/cal files."\
                " Current 'filename_parser' instance: \n{0}\n"\
                "List of cal/ped files:\n{1}".format(fn_instance,message))

        self.beam_instance        = fn_instance
        self.pedestal_instance    = filter(lambda x: x.is_pedestal,associated_files)[0]
        self.calibration_instance = filter(lambda x: x.is_calibration,associated_files)[0]

    def __str__(self):
        message = "Beam instance: {0}\n".format(self.beam_instance)
        message += "  - Beam file: {0}\n".format(self.beam_instance.filename)
        message += "  - Pedestal file:{0}\n".format(self.pedestal_instance.filename)
        message += "  - Calibration file:{0}".format(self.calibration_instance.filename)
        return message
# -----------------------------------------------------------------------------

# list of active sensors per run number, # they are ordered with the
# beam hitting in increasing z-plane (Motherboard). The first 
# element of the 2-tuple is the name of the sensor while the second 
# is the beetle number (0,1). Note that the name of the sensor contains
# the beetle number b1-b2
active_sensors = {
    274: [ ('LGAD7859W1H6_0_b1',0), ('REF_0_b1',0), ('M2-3_1e16_b2',1), ('N1-7_7e15_b2',1) ],
    345: [ ('iLGAD8533W1K05T_0_b2',1),('REF_0_b1',0),('N1-8_1e16_b1',0),('N1-7_7e15_b2',1) ],
    362: [ ('iLGAD8533W1K05T_0_b2',1),('REF_0_b1',0),('M2-3_1e16_b2',1),('N1-7_7e15_b2',1) ],
    378: [ ('REF_0_b1',0),('M2-3_1e16_b2',1), ('M1-5_0_b2',1), ('M1-8_7e15_b2',1) ],
    391: [ ('REF_0_b1',0),('N1-8_1e16_b1',0), ('N1-3_0_b1',0), ('N1-7_7e15_b2',1) ], 
    }
# append those with the same configuration
_samecfg = { 
    274: range(275,285)+[286]+range(299,303)+range(304,311)+range(313,316),
    345: range(346,353),
    362: range(363,373),
    378: range(379,386),
    391: range(392,296)+range(399,410) 
    }
# Build the final map
for k,runlist in _samecfg.iteritems():
    for irun in runlist:
        active_sensors[irun] = active_sensors[k]


