#!/usr/bin/env python
"""Python wrapper to the AlibavaSensorAnalysis C-code. 
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

def callit_once(f):
    """Decorator  to assure that the function argument is called 
    only once, otherwise raises an exception
    """
    def wrapper(*args, **kwargs):
        if not wrapper.has_run:
            wrapper.has_run = True
            return f(*args, **kwargs)
        wrapper.has_run = False
        raise RuntimeError("Unexpected call to the method '{0}'"\
                " which can be only called ONCE".format(f.__name__))

class alibava_analysis(object):
    """Wrapper to the C++ class IOFortythieves and AlibavaSensorAnalysis 
    which implement the created ROOT file by the `fortythieves` exec. 
    centralized access, and the algorithms to perform cluster reconstruction,
    finding, positioning and analysis (eta, charge, ...)
    """
    def __init__(self,filename,chip):
        """Initialization
        """
        import ctypes

        # some useful attributes
        # -- list of (python9 functions which can only 
        #    be called once
        self._called_functions = []

        # Library
        self._lib = ctypes.cdll.LoadLibrary("libAlibavaSensorAnalysis.so")
        # All involved types for the IOFortythieves
        self._lib.ioft_new.restype = ctypes.c_void_p
        self._lib.ioft_new.argtypes = [ ctypes.c_char_p, ctypes.c_int ]
        self._lib.ioft_delete.argtypes = [ ctypes.c_void_p ]
        self._lib.ioft_initialize.argtypes = [ ctypes.c_void_p ]
        self._lib.ioft_process.argtypes = [ ctypes.c_void_p, ctypes.c_int ]
        self._lib.ioft_get_entries.argtypes = [ ctypes.c_void_p ]
        self._lib.ioft_get_entries.restype = ctypes.c_int 
        # All involved types for the AlibavaSensorAnalysis
        self._lib.aa_new.restype = ctypes.c_void_p
        self._lib.aa_new.argtypes = [ ctypes.c_int ]
        self._lib.aa_delete.argtypes = [ ctypes.c_void_p ]
        self._lib.aa_configure_time_cut.argtypes = [ ctypes.c_void_p, ctypes.c_float, ctypes.c_float ]
        self._lib.aa_configure_masked_channels.argtypes = [ ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_int ]

        # The objects
        self._ioft              = self._lib.ioft_new(filename,chip)
        self._sensor_analysis   = self._lib.aa_new(self._ioft)

    def __del__(self):
        """Freeing memory by calling the destructor function
        """
        self._lib.aa_delete(self._sensor_analysis)
        self._lib.ioft_delete(self._ioft)

    def initialize(self):
        """Wrapper to the initialization:
        The noise, pedestal and calibration data are filled
        """
        self._lib.ioft_initialize(self._ioft)

    def get_entries(self):
        """
        """
        return self._lib.ioft_get_entries(self._ioft)
    
    def process(self,i):
        """
        """
        return self._lib.ioft_process(self._ioft,i)
    

    #@callit_once
    def set_time_cut(self,t0,t1):
        """Set the TDC time cut 

        Parameters
        ----------
        t0: float
            The minimum TDC time
        t1: float 
            The maximum TDC time
        """
        self._lib.aa_configure_time_cut(self._sensor_analysis,t0,t1)

    def set_masked_channels(self,masked_chan):
        """User defined masked channels

        Parameters
        ----------
        masked_chan: list(int)
            The list of masked channels 
        """
        import ctypes
        # build the string to create the array of ints
        mch_arr_eval = "(ctypes.c_int*len(masked_chan))("
        for ch in masked_chan:
            mch_arr_eval += "{0},".format(ch)
        mch_arr_eval = mch_arr_eval[:-1]+")"
        mch_arr = eval(mch_arr_eval)
        # Actually call the function 
        self._lib.aa_configure_masked_channels(self._sensor_analysis,mch_arr,len(masked_chan))

