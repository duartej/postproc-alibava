#!/usr/bin/env python
"""Several dictionaries and functions with static information from 
the SPS Test Beam performed at May, 2017 at CERN 

Information extracted partially from 
https://docs.google.com/spreadsheets/d/1Z4nlyHUdAhCy-oNC-472c0Sydg54GraveDxROsnZKrM/edit#gid=0
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "1.0"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

__all__ = [
        "eospath",
        "run_numbers",
        "filename_parser",
        "associated_filenames",
        "sensor_names",
        "sensor_ids",
        "active_channels",
        "equivalent_run_number",
        "get_active_sensor_list",
        "get_beetle",
        "get_setup",
        "get_gear_content"
        ]

# the EOS path where to find the alibava data
eospath="/eos/user/d/duarte/alibava_data"

# Run numbers defined at https://docs.google.com/spreadsheets/d/1Z4nlyHUdAhCy-oNC-472c0Sydg54GraveDxROsnZKrM/
run_numbers = range(275,284+1)+[286]+range(299,302+1)+range(304,315+1)+\
        range(345,352+1)+range(362,372+1)+range(378,385+1)+\
        range(391,395+1)+range(399,410+1)
# Remove some bad file not able to be processed
for badrun in [ 286,311,312 ]:
    run_numbers.remove(badrun)

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
        if abs(our_hour-other_hour) > 15.0*3600.0:
            return False
        # The propose float equalizers :
        equalizers = [ ('voltage_bias',1e-19), ('current_leak',20.0), 
                ('temperature',1e-19) ]
        # [HARDCODED CASES]
        ###################
        if self.sensor_name == "M2-3" and other.run_number in [ 282,283,284]:
            equalizers[0] = ('voltage_bias',5.0)
            equalizers[1] = ('current_leak', 400.0)
        # [END HARDCODED CASES]
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
        self.run_number           = fn_instance.run_number

    def __str__(self):
        message = "Beam instance: {0}\n".format(self.beam_instance)
        message += "  - Beam file: {0}\n".format(self.beam_instance.filename)
        message += "  - Pedestal file:{0}\n".format(self.pedestal_instance.filename)
        message += "  - Calibration file:{0}".format(self.calibration_instance.filename)
        return message

# the list of used sensors and an integer identifying them
sensor_names = [ 'LGAD7859W1H6_0_b1', 'iLGAD8533W1K05T_0_b2', 'REF_0_b1', 'M1-5_0_b2', 
        'M1-8_7e15_b2', 'M2-3_1e16_b2', 'N1-3_0_b1', 'N1-7_7e15_b2', 'N1-8_1e16_b1' ]
sensor_ids = dict(map(lambda (x,name): (name,x),enumerate(sensor_names)))
# A map to convert the names without the fluence nor the beetle
standard_sensor_name_map = { 'LGAD7859W1H6': 'LGAD7859W1H6_0_b1', 
        'iLGAD8533W1K05T': 'iLGAD8533W1K05T_0_b2',
        'REF': 'REF_0_b1', 
        'M1-5':'M1-5_0_b2', 'M1-8': 'M1-8_7e15_b2', 'M2-3': 'M2-3_1e16_b2',
        'N1-3':'N1-3_0_b1', 'N1-7':'N1-7_7e15_b2','N1-8':'N1-8_1e16_b1'
        }

# The list of active channels 
active_channels = dict(map(lambda x: (x,[0,127]),sensor_names))
# --- LGAD and iLGAD are not using all of them
active_channels['iLGAD8533W1K05T_0_b2'] = [28,72]
active_channels['LGAD7859W1H6_0_b1'] = [41,70]

# -----------------------------------------------------------------------------
# Some characteristics of the sensors
def _binary_resolution(pitch):
    """Extract the binary resolution
    """
    from math import sqrt
    return float(pitch)/sqrt(12.0)

class specs_sensor():
    def __init__(self,sizeX,sizeY,pitchX,pitchY,thickness,polarity=-1.0):
        # All units are mm
        # Note that the size is standarize to the 128 sensor 
        # strip length in order to keep the same indices along 
        # the processors
        force_strips = 128
        self.sizeX = force_strips*pitchX #sizeX
        self.sizeY = sizeY
        self.pitchX = pitchX
        self.pitchY = pitchY
        self.thickness = thickness
        self.polarity = polarity 
        self.resolution = _binary_resolution(self.pitchX)
# Instances
mtype = specs_sensor(6.4,7.5,0.05,0.05,0.23)
ntype = specs_sensor(3.2,7.5,0.0250,0.100,0.23)
lgad  = specs_sensor(5.12,5.12,0.160,5.12,0.3)
ilgad = specs_sensor(7.2,7.2,0.160,7.2,0.3,polarity=1.0)
ref   = specs_sensor(10.24,10.24,0.08,10.24,0.3)

# Maps the name of the sensor with the proper specs_sensor instance
sensor_name_spec_map = { 'REF_0_b1': ref, 'LGAD7859W1H6_0_b1': lgad,\
        'iLGAD8533W1K05T_0_b2': ilgad,\
        'M2-3_1e16_b2' : mtype, 'M1-5_0_b2': mtype, 'M1-8_7e15_b2': mtype, \
        'N1-3_0_b1': ntype, 'N1-7_7e15_b2': ntype, 'N1-8_1e16_b1': ntype
        }

# list of active sensors per run number. They are ordered with the
# beam hitting in increasing z-plane (Motherboard). The first 
# element of the 3-tuple is the name of the sensor while the second 
# is the beetle number (0,1); the third element is the z-position in mm taking
# as z=0 the first plane of the telescope.
# Note that the name of the sensor contains the beetle number b1-b2
active_sensors = {
    315: [ ('LGAD7859W1H6_0_b1',0,-158.0), ('REF_0_b1',0,-106.0), ('M2-3_1e16_b2',1,218.0), ('N1-7_7e15_b2',1,339.0) ],
    352: [ ('iLGAD8533W1K05T_0_b2',1,-158.0),('REF_0_b1',0,-106.0),('N1-8_1e16_b1',0,218.0) ],
    372: [ ('iLGAD8533W1K05T_0_b2',1,-158.0),('REF_0_b1',0,-106.0),('M2-3_1e16_b2',1,218.0),('N1-7_7e15_b2',1,339.0) ],
    385: [ ('REF_0_b1',0,-106.0),('M2-3_1e16_b2',1,218.0), ('M1-5_0_b2',1,278.0), ('M1-8_7e15_b2',1,339.0) ],
    410: [ ('REF_0_b1',0,-106.0),('N1-8_1e16_b1',0,218.0), ('N1-3_0_b1',0,278.0), ('N1-7_7e15_b2',1,339.0) ], 
    }

# The gear file content template. The fields to be filled are defined as
# ---------------------------------------------------------------
# 0: "WITHOUT DUT"|"WITH DUT <dut_name>
# 1: "XX"|<dut_name>
# 2: "XX"|<dut z-position>
# 3: 1XXXY (the GeoID, XXX: equivalent_run_number, Y: sensor_id), 
#    no DUTs means Y=0
# 4: 5|6 (with DUT)
# 5: ""|"Mimosa26.so" (with DUT)
# 6: ""|gear_dut_template(sensor_name)
# ---------------------------------------------------------------
# Note that the first column is referring to a gear with Telescope
# only, and the second one, with a given DUT
gear_content_template = """
<gear>
  <!--
      GEAR file for CMSPixel CERN SPS July 2017 testbeam setup
      https://docs.google.com/spreadsheets/d/1Z4nlyHUdAhCy-oNC-472c0Sydg54GraveDxROsnZKrM/edit#gid=0

      TELESCOPE {0}

                                                       |<< DUT >>|
      SETUP:                             REF  T0   T1     {1}      T2   T3   T4
      The positions (mm) z   :          -106  0   101     {2}    493  594   694
  -->

  <global detectorName="EUTelescope"/>
  <BField type="ConstantBField" x="0.0" y="0.0" z="0.0"/>
  <detectors>
    <detector name="SiPlanes" geartype="SiPlanesParameters">
      <siplanesID ID="{3}"/>
      <siplanesType type="TelescopeWithoutDUT"/>
      <siplanesNumber number="{4}"/>
      <parameter name="Geometry" 
                 type="StringVec" 
                 value="{5} Mimosa26.so Mimosa26.so Mimosa26.so Mimosa26.so Mimosa26.so"/>
      <layers>
	<!--DATURA-Plane 0 - EUD0 -->
	<layer>
	  <ladder 	ID="0"
			positionX="0.00"	positionY="0.00"	positionZ="0.00"
			rotationZY="0.00"	rotationZX="0.0"	rotationXY="0.0" 
			sizeX="21.2"		sizeY="10.6"		thickness="0.050"
			radLength="93.660734"
			/>
	  <sensitive 	ID="0"
			positionX="0.00"	positionY="0.00"	positionZ="0.00" 
			sizeX="21.2"		sizeY="10.6"		thickness="0.070"
			npixelX="1152"		npixelY="576" 
			pitchX="0.018402778"	pitchY="0.018402778" 	resolution="0.0045" 
			rotation1="-1.0" 	rotation2="0.0" 
			rotation3="0.0"		rotation4="-1.0" 
			radLength="93.660734"
			/>
	</layer>
	<!--DATURA-Plane 1 - EUD1 -->
	<layer>
	  <ladder 	ID="1" 
			positionX="0.00"	positionY="0.00"	positionZ="101.0" 
			rotationZY="0.0"	rotationZX="0.0"	rotationXY="0.0" 
			sizeX="21.2"		sizeY="10.6"		thickness="0.050" 
			radLength="93.660734"
			/>
	  <sensitive 	ID="1" 
			positionX="0.00"	positionY="0.00"	positionZ="101.0" 
			sizeX="21.2"		sizeY="10.6"		thickness="0.070"
			npixelX="1152"		npixelY="576" 
			pitchX="0.018402778"	pitchY="0.018402778"	resolution="0.0045" 
			rotation1="-1.0"	rotation2="0.0" 
			rotation3="0.0"		rotation4="-1.0" 
			radLength="93.660734"
			/>
	</layer>
	<!--DATURA-Plane 2 - EUD2 -->
	<layer>
	  <ladder 	ID="2" 
			positionX="0.00"	positionY="0.00"	positionZ="493.0" 
			rotationZY="0.0"	rotationZX="0.0"	rotationXY="0.0" 
			sizeX="21.2"		sizeY="10.6"		thickness="0.050" 
			radLength="93.660734"
			/>
	  <sensitive 	ID="2" 
			positionX="0.00"	positionY="0.00"	positionZ="493.0" 
			sizeX="21.2"		sizeY="10.6"		thickness="0.070"
			npixelX="1152"		npixelY="576"
			pitchX="0.018402778" 	pitchY="0.018402778" 	resolution="0.0045"
			rotation1="-1.0"	rotation2="0.0"
			rotation3="0.0" 	rotation4="-1.0" 
			radLength="93.660734"
			/>
	</layer>
	
        <!--DATURA-Plane 3 - EUD0 -->
        <layer>
          <ladder         ID="3"
                          positionX="0.00"        positionY="0.00"        positionZ="594.0"
                          rotationZY="0.00"       rotationZX="0.0"        rotationXY="0.0" 
                          sizeX="21.2"            sizeY="10.6"            thickness="0.050"
                          radLength="93.660734"
                          />
          <sensitive      ID="3"
                          positionX="0.00"        positionY="0.00"        positionZ="594.0" 
                          sizeX="21.2"            sizeY="10.6"            thickness="0.070"
                          npixelX="1152"          npixelY="576" 
                          pitchX="0.018402778"    pitchY="0.018402778"    resolution="0.0045" 
                          rotation1="-1.0"        rotation2="0.0" 
                          rotation3="0.0"         rotation4="-1.0" 
                          radLength="93.660734"
                          />
        </layer>
        <!--DATURA-Plane 4 - EUD0 -->
        <layer>
          <ladder         ID="4"
                          positionX="0.00"        positionY="0.00"        positionZ="694.00"
                          rotationZY="0.00"       rotationZX="0.0"        rotationXY="0.0" 
                          sizeX="21.2"            sizeY="10.6"            thickness="0.050"
                          radLength="93.660734"
                          />
          <sensitive      ID="4"
                          positionX="0.00"        positionY="0.00"        positionZ="694.00" 
                          sizeX="21.2"            sizeY="10.6"            thickness="0.070"
                          npixelX="1152"          npixelY="576" 
                          pitchX="0.018402778"    pitchY="0.018402778"    resolution="0.0045" 
                          rotation1="-1.0"        rotation2="0.0" 
                          rotation3="0.0"         rotation4="-1.0" 
                          radLength="93.660734"
                          />

        </layer>
        {6}
      </layers>
    </detector>
  </detectors>
</gear>
"""
# The DUT layer template to be place in the gear template.
# 0: dut_name
# 1: chip number
# 2: size in X-direction (sensitive direction)
# 3: size in Y-direction
# 4: pitch in X-direction (sensitive direction)
# 5: pitch in Y-direction
# 6: spatial resolution (just the binary resolution)
# 7: thickness
gear_dut_template="""<!--{0} - chip {1} -->
        <!-- WARNING, not sure about specs, first tentative values collected from several 
             sources. Resolution: (binary resolution, i.e. p/sqrt(12))-->
        <layer> 
          <ladder         ID="{9}"
                          positionX="0.00"        positionY="0.0"       positionZ="{8}"
                          rotationZY="0.00"       rotationZX="0.0"      rotationXY="0.0" 
                          sizeX="{2}"           sizeY="{3}"           thickness="{7}"
                          radLength="93.660734"
                          />
          <sensitive      ID="{9}"
                          positionX="0.00"        positionY="0.00"      positionZ="{8}" 
                          sizeX="{2}"             sizeY="{3}"           thickness="0.230"
                          npixelX="128"           npixelY="1" 
                          pitchX="{4}"           pitchY="{5}"        resolution="{6}" 
                          rotation1="-1.0"        rotation2="0.0" 
                          rotation3="0.0"         rotation4="-1.0" 
                          radLength="93.660734"
                          />
        </layer>"""

# ------------------------------------------------------------------------------
# Useful functions
def get_standard_sensor_name(sensor_name):
    """Given the name of a sensor, returns the standarized name,
    which includes the fluence and the beetle where it was bonded

    Paramaters
    ----------
    sensor_name: str
        the (non-standard) sensor name
    
    Returns
    -------
    str: the standarized sensor name
    """
    try:
        return standard_sensor_name_map[sensor_name]
    except IndexError:
        raise RuntimeError("Invalid sensor name '{0}'".format(sensor_name))
    

# Mapping a run number with its setup
# Available info at https://docs.google.com/spreadsheets/d/1Z4nlyHUdAhCy-oNC-472c0Sydg54GraveDxROsnZKrM/edit#gid=343850123
# The setups are linked to a range of run numbers. The equivalent_run_number
# is defined to be as the maximum run number of each range
#  * run number <= 315    --> Setup-0
#  * run number (315,352] --> Setup-1
#  * ...
setups = [ 315, 352, 372, 385, 410 ]
def equivalent_run_number(run_number):
    """Helper function to obtain the equivalent run number given a run 
    number. Setups are linked to a range of run numbers. The equivalent
    run number is defined to be as the maximum run number of each range

    Parameters
    ----------
    run_number: int
        the run  number 

    Returns
    -------
    int: the equivalent run number

    Raises
    ------
    RuntimeError
        if the run number is higher than the maximum recorded
    """
    for (k,maxrunnumber) in enumerate(setups):
        if int(run_number) <= maxrunnumber:
            return maxrunnumber
    raise RuntimeError("Invalid run number: {0} (> {1}, max. stored)".format(run_number,setups[-1]))

def get_active_sensor_list(run_number=-1):
    """Helper function to obtain the active sensor given a run number
    If no run number is provided, then all the available sensors
    are returned

    Note
    ----
    When returning the list without a run number assumes
    that all the sensors do not change between runs, which
    it is the case so far, but it could change in the future

    Parameters
    ----------
    run_number: int, optional
        the run number

    Returns
    -------
    list(ntuple(str,int,float)): the name, the chip and the z-position
    """
    if run_number == -1:
        sensor_list = []
        for sensor_ntuple_list in active_sensors.values():
            for sensor_ntuple in sensor_ntuple_list:
                sensor_list.append( sensor_ntuple )
    else:
        ern = equivalent_run_number(run_number)
        sensor_list = active_sensors[ern]
    return sensor_list


def get_beetle(sensor_name):
    """Helper function to obtain the beetle given the sensor name

    Parameters
    ----------
    sensor_name: str

    Returns
    -------
    int: the beetle id

    Raises
    ------
    RuntimeError
        if the run number is higher than the maximum recorded
    """
    sensorlist = get_active_sensor_list()
    try:
        return filter(lambda (sname,beetle,zpos): sname == sensor_name,sensorlist)[0][1]
    except IndexError:
        raise RuntimeError("Invalid sensor name '{0}' ".format(sensor_name))

def get_z(run_number,sensor_name):
    """Helper function to obtain z-position of a sensor

    Parameters
    ----------
    run_number: int
        the run number 
    sensor_name: str
        

    Returns
    -------
    int: the z position [mm]

    Raises
    ------
    RuntimeError
        if the run number is higher than the maximum recorded
    """
    sensorlist = get_active_sensor_list(run_number)
    try:
        return filter(lambda (sname,beetle,zpos): sname == sensor_name,sensorlist)[0][2]
    except IndexError:
        raise RuntimeError("Invalid sensor name '{0}' OR not present in the run '{1}' ".format(sensor_name,run_number))

def get_setup(run_number):
    """Helper function to obtain the setup given a run number

    Parameters
    ----------
    run_number: int
        the run  number 

    Returns
    -------
    int: the setup index
    """
    for (k,maxrunnumber) in enumerate(setups):
        if int(run_number) <= maxrunnumber:
            return k
    raise RuntimeError("Invalid run number '{0}' > {0} (max. stored)".format(run_number,setups[-1]))

def get_geo_id(run_number,sensor_name=""):
    """Get the geoID to be related with a particular setup: 
        ZXXXY ->   Z: 1-telescope only, 2-telescope+DUT
                 XXX: equivalent_run_number 
                   Y: sensor_id 

    Parameters
    ----------
    run_number: int
        the run  number 
    sensor_name: str
        the name of the DUT to be included in the gear file. 
        Empty string assumes only Telescope with no DUT

    Returns
    -------
    int: the geoId (ZXXXY)
    """
    if sensor_name == "":
        return int("1{0}0".format(equivalent_run_number(run_number)))
    else:
        return int("2{0}{1}".format(equivalent_run_number(run_number),sensor_ids[sensor_name]))


def get_gear_content(run_number,sensor_name="",include_ref=False):
    """Get the gear file corresponding to a run number

    Parameters
    ----------
    run_number: int
        the run  number 
    sensor_name: str
        the name of the DUT to be included in the gear file. 
        Empty string assumes only Telescope with no DUT

    Returns
    -------
    str: the gear template filled

    Raises
    ------
    RuntimeError
        if the run number is higher than the maximum recorded
    """
    ern = equivalent_run_number(run_number)
    if sensor_name == "":
        # Only telescope (see gear_content_template)
        geoid = get_geo_id(equivalent_run_number(run_number))
        filler_gear = ("WITHOUT DUT","XXX","XXX",geoid,5,"","")
    else:
        # Sensor id
        sensorID = 5
        if sensor_name.find("REF") == 0:
            sensorID = 6
        # A DUT is defined, see gear_dut_template
        # Find which instance we need
        specs=sensor_name_spec_map[sensor_name]
        # Fill the need input of the template for the dut layer
        filler_dut = (sensor_name,get_beetle(sensor_name),\
                specs.sizeX,specs.sizeY,specs.pitchX,specs.pitchY,specs.resolution,\
                specs.thickness,get_z(run_number,sensor_name),sensorID)
        dut_layer = gear_dut_template.format(*filler_dut)
        # And for the Gear file
        geoid = get_geo_id(equivalent_run_number(run_number),sensor_name)
        #filler_gear = ("WITH DUT "+sensor_name, sensor_name,get_z(run_number,sensor_name),\
        #        geoid,6,"Mimosa26.so",dut_layer)
        
        # And if the REF should be include
        if include_ref:
            sensorID_ref = 6
            specs_ref=sensor_name_spec_map['REF_0_b1']
            filler_ref = ("REF_0_b1", get_beetle("REF_0_b1"),\
                    specs_ref.sizeX,specs_ref.sizeY,specs_ref.pitchX,specs_ref.pitchY,specs_ref.resolution,\
                    specs_ref.thickness,get_z(run_number,"REF_0_b1"),sensorID_ref)
            ref_layer = gear_dut_template.format(*filler_ref)
            # And for the Gear file
            filler_gear = ("WITH DUT "+sensor_name, sensor_name,get_z(run_number,sensor_name),\
                    geoid,7,"Mimosa26.so Mimosa26.so",dut_layer+"\n"+ref_layer)
        else:
            filler_gear = ("WITH DUT "+sensor_name, sensor_name,get_z(run_number,sensor_name),\
                    geoid,6,"Mimosa26.so",dut_layer)


    return gear_content_template.format(*filler_gear)
