#!/usr/bin/env python
"""Several dictionaries with static information from the SPS
Test Beam performed at May, 2017 at CERN 

Dictionary with the pedestal and calibration file corresponding to
a given run number. Information extracted partially from 
https://docs.google.com/spreadsheets/d/1Z4nlyHUdAhCy-oNC-472c0Sydg54GraveDxROsnZKrM/edit#gid=0
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

# the EOS path where to find the alibava data
eospath="/eos/user/d/duarte/alibava_data"

# list of active sensors per run number,
# they are ordered with the beam hitting in 
# increasing z-plane (Motherboard). The first
# element of the 2-tuple is the name of the sensor
# while the second is the beetle number (0,1) 
active_sensors = {
        274: [ ('LGAD7859W1H6_0_b1',0), ('REF_0_b1',0), ('M2-3_1e16_b2',0), ('N1-7_7e15_b2',0) ],
        345: [ ('iLGAD8533W1K05T_0_b2',1),('REF_0_b1',0),('N1-8_1e16_b1',0),('N1-7_7e15_b2',0) ],
        362: [ ('iLGAD8533W1K05T_0_b2',1),('REF_0_b1',0),('M2-3_1e16_b2',0),('N1-7_7e15_b2',0) ],
        378: [ ('REF_0_b1',0),('M2-3_1e16_b2',0), ('M1-5_0_b2',0), ('M1-8_7e15_b2',0) ],
        391: [ ('REF_0_b1',0),('N1-8_1e16_b1',0), ('N1-3_0_b1',0), ('N1-7_7e15_b2',0) ], 
}
# append those with the same configuration
samecfg = { 
        274: range(275,285)+[286]+range(299,303)+range(304,311)+range(313,316),
        345: range(346,353),
        362: range(363,373),
        378: range(379,386),
        391: range(392,296)+range(399,410) 
}
# Build the final map
for k,runlist in samecfg.iteritems():
    for irun in runlist:
        active_sensors[irun] = active_sensors[k]

# 
#pedcalfiles = {
#       274: (,),
#        410: (,)
#}

