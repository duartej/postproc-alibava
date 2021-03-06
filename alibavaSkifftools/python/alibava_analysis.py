#!/usr/bin/env python
"""Alibava cluster analysis from a ROOT file created using
the `fortythieves`. 
Python wrapper to the AlibavaSensorAnalysis C-code. 
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

import ctypes

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

# Helper classes to avoid the segfault due to the returned
# datatypes in the AlibavaAnalysisSensor. Not needed in principle
# for the other two (IOFortythieves and IOASAResults), but added
# and separated just in case. (See https://stackoverflow.com/\
#        questions/20806070/ctypes-correctly-sublcass-c-void-\
#        p-for-passing-and-returning-custom-data-types )
# XXX: The issue is not fully understood (specially the different
#      behaviours depending the machine/compiler)
class c_ioft(ctypes.c_void_p):
    pass

class c_aa(ctypes.c_void_p):
    pass

class c_iores(ctypes.c_void_p):
    pass

class alibava_analysis(object):
    """Wrapper to the C++ class IOFortythieves and AlibavaSensorAnalysis 
    which implement the created ROOT file by the `fortythieves` exec. 
    centralized access, and the algorithms to perform cluster reconstruction,
    finding, positioning and analysis (eta, charge, ...)
    """
    def __init__(self,filename,chip):
        """Initialization
        """
        # some useful attributes
        # -- whether or not the results should be stored in a file
        #    using the IOASAResults class
        self._store_results = False
        # -- the iterations to be applied for the xtcorrection
        self.xtcorrection   = None

        # Library
        self._lib = ctypes.cdll.LoadLibrary("libAlibavaSensorAnalysis.so")
        # All involved types for the IOFortythieves
        self._lib.ioft_new.restype = c_ioft 
        self._lib.ioft_new.argtypes = [ ctypes.c_char_p, ctypes.c_int ]
        self._lib.ioft_delete.restype = None
        self._lib.ioft_delete.argtypes = [ c_ioft ]
        self._lib.ioft_initialize.restype = None
        self._lib.ioft_initialize.argtypes = [ c_ioft ]
        self._lib.ioft_process.restype = None
        self._lib.ioft_process.argtypes = [ c_ioft, ctypes.c_int ]
        self._lib.ioft_get_entries.restype = ctypes.c_int 
        self._lib.ioft_get_entries.argtypes = [ c_ioft ]
        # All involved types for the AlibavaSensorAnalysis
        # -- constructor and destructor
        self._lib.aa_new.restype = c_aa
        self._lib.aa_new.argtypes = [ ]
        self._lib.aa_delete.argtypes = [ c_aa ]
        # All involved types for the IOASAResults
        self._lib.ioresults_new.restype   = c_iores
        self._lib.ioresults_new.argtypes = [ ctypes.c_char_p ]
        self._lib.ioresults_delete.restype = None
        self._lib.ioresults_delete.argtypes = [ c_iores ]
        self._lib.ioresults_book_tree.restype = None
        self._lib.ioresults_book_tree.argtypes = [ c_iores ]
        self._lib.ioresults_fill_header.restype = None
        self._lib.ioresults_fill_header.argtypes = [ c_iores, c_ioft ]
        self._lib.ioresults_fill_tree.restype = None
        self._lib.ioresults_fill_tree.argtypes = [ c_iores, c_ioft, c_aa ]
        # -- some setters
        self._lib.aa_configure_polarity.argtypes = [ c_aa, ctypes.c_int ]
        self._lib.aa_configure_time_cut.argtypes = [ c_aa, ctypes.c_float, ctypes.c_float ]
        self._lib.aa_configure_masked_channels.argtypes = [ c_aa, ctypes.POINTER(ctypes.c_int), ctypes.c_int ]
        self._lib.aa_configure_mask_criterium.argtypes = [ c_aa, ctypes.c_float ]
        self._lib.aa_configure_snr_seed.argtypes = [ c_aa, ctypes.c_float ]
        self._lib.aa_configure_snr_neighbour.argtypes = [ c_aa, ctypes.c_float ]
        # -- getters of configuration member
        self._lib.aa_polarity_getter.argtype = [ c_aa ]
        self._lib.aa_polarity_getter.restype = ctypes.c_int 
        self._lib.aa_time_cut_getter.argtype = [ c_aa ]
        self._lib.aa_time_cut_getter.restype = ctypes.POINTER(ctypes.c_float) 
        self._lib.aa_number_masked_channels.argtype = [ c_aa ]
        self._lib.aa_number_masked_channels.restype = ctypes.c_int
        self._lib.aa_masked_channels_getter.argtype = [ c_aa ]
        self._lib.aa_masked_channels_getter.restype = ctypes.POINTER(ctypes.c_int)
        self._lib.aa_mask_criterium_getter.argtype = [ c_aa ]
        self._lib.aa_mask_criterium_getter.restype = ctypes.c_float
        self._lib.aa_snr_seed_getter.argtypes = [ c_aa ]
        self._lib.aa_snr_seed_getter.restype = ctypes.c_float
        self._lib.aa_snr_neighbour_getter.argtypes = [ c_aa ]
        self._lib.aa_snr_neighbour_getter.restype = ctypes.c_float
        # print configuration
        self._lib.aa_print_configuration.argtypes = [ c_aa, c_ioft ]
        self._lib.aa_print_configuration.restype  = None

        # -- modifiers
        self._lib.aa_mask_channels.argtypes = [ c_aa, c_ioft ]
        self._lib.aa_mask_channels.restype  = None
        self._lib.aa_update_crosstalk_factors.argtypes = [ c_aa ]
        self._lib.aa_update_crosstalk_factors.restype = ctypes.c_int 

        # The objects
        self._ioft            = self._lib.ioft_new(filename,chip)
        self._sensor_analysis = self._lib.aa_new()
        self._results         = None 

    def __del__(self):
        """Freeing memory by calling the destructor function
        """
        self._lib.aa_delete(self._sensor_analysis)
        self._lib.ioft_delete(self._ioft)
        if self._store_results:
            self._lib.ioresults_delete(self._results)
    
    # -- Functions related with IOFortythieves
    def get_entries(self):
        """
        """
        return self._lib.ioft_get_entries(self._ioft)
    
    # XXX: Should be not exposed? Given that it is used by
    # the analysis_sensor class? (see find_clusters)
    def process_event(self,i):
        """
        """
        return self._lib.ioft_process(self._ioft,i)
    
    # -- Functions related with AlibavaSensorAnalysis
    # ==============================================
    # Configurables
    @property
    def polarity(self):
        """Get the signal polarity, using the C-function 
        (and therefore using the value contained in the 
        data member of the C++ class)

        Return
        ------
        int: the signal polarity (a value 0 implies no initialization
            or configuration of the class)
        """
        # Get the pointer object and convert it to a 2-tuple
        return self._lib.aa_polarity_getter(self._sensor_analysis)
    @polarity.setter
    def polarity(self,p):
        """Set the signal polarity of this sensor

        Parameters
        ----------
        t1: int
            The signal polarity, -1|+1
        """
        self._lib.aa_configure_polarity(self._sensor_analysis,p)
    

    @property
    def time_cut(self):
        """Get the TDC time cuts, using the C-function 
        (and therefore using the value contained in the 
        data member of the C++ class)

        Return
        ------
        tuple(float,float): the low and up allowed times
        """
        # Get the pointer object and convert it to a 2-tuple
        fpointer = self._lib.aa_time_cut_getter(self._sensor_analysis)
        return (fpointer[0],fpointer[1])

    @time_cut.setter
    def time_cut(self,t0,t1):
        """Set the TDC time cut 

        Parameters
        ----------
        t0: float
            The minimum TDC time
        t1: float 
            The maximum TDC time
        """
        self._lib.aa_configure_time_cut(self._sensor_analysis,t0,t1)
    
    @property
    def masked_channels(self):
        """Get the masked channels, using the C-function (and 
        therefore using the value contained in the data member
        of the C++ class)

        Return
        ------
        masked_channels: n-tuple(int) | None
        """
        # First obtain the number of channels, if any
        n = self._lib.aa_number_masked_channels(self._sensor_analysis)
        if n == 0 :
            return None
        # Obtain the pointer
        mch_pointer = self._lib.aa_masked_channels_getter(self._sensor_analysis)
        # And now convert it to a list (before return the n-tuple
        return tuple(map(lambda i: mch_pointer[i] ,xrange(n)))
    
    @masked_channels.setter
    def masked_channels(self,masked_chan):
        """User defined masked channels

        Parameters
        ----------
        masked_chan: list(int)
            The list of masked channels 
        """
        # A check before
        if len(masked_chan) == 0:
            raise RuntimeError("Reset to no-masked channels is"\
                    " Not implemented yet") 
        # build the string to create the array of ints
        # (avoiding to deal with std::vectors)
        mch_arr_eval = "(ctypes.c_int*len(masked_chan))("
        for ch in masked_chan:
            mch_arr_eval += "{0},".format(ch)
        mch_arr_eval = mch_arr_eval[:-1]+")"
        mch_arr = eval(mch_arr_eval)
        # Actually call the function 
        self._lib.aa_configure_masked_channels(self._sensor_analysis,mch_arr,len(masked_chan))
    
    @property
    def automask_criterium(self):
        """Get the masked channels criterium used, i.e. how many sigmas
        away from the mean noise, a channel is to be considered a noisy 
        channel. Note that this setter is using the C-function (and 
        therefore using the value contained in the data member of the 
        C++ class):
        The noisy channel condition is given by:

        .. math:: i-strip is noisy \Leftrightarrow |N_{i} - \langle Noise\\rangle| > X\sigma

        Return
        ------
        float: The X 
        """
        # Get the pointer object and convert it to a 2-tuple
        return self._lib.aa_mask_criterium_getter(self._sensor_analysis)

    @automask_criterium.setter
    def automask_criterium(self,sigma):
        """Set the masked channels criterium used, i.e. how many sigmas
        away from the mean noise, a channel is to be considered a noisy 
        channel. Note that this setter is using the C-function (and 
        therefore using the value contained in the data member of the 
        C++ class):
        The noisy channel condition is given by:

        .. math:: i-strip is noisy \Leftrightarrow |N_{i} - \langle Noise\\rangle| > X\sigma

        Parameters
        ----------
        sigma: float
        """
        self._lib.aa_configure_mask_criterium(self._sensor_analysis,sigma)

    @property
    def snr_seed(self):
        """Get the Signal-to-Noise Ratio (SNR) for a seed-cluster
        candidate, using the C-function (and therefore using the
        value contained in the data member of the C++ class)
        """
        return self._lib.aa_snr_seed_getter(self._sensor_analysis)

    @snr_seed.setter
    def snr_seed(self,snr):
        """Set the Signal-to-Noise Ratio (SNR) for a seed-cluster
        candidate

        Parameters
        ----------
        snr: float
            The SNR
        """
        self._lib.aa_configure_snr_seed(self._sensor_analysis,snr)

    @property
    def snr_neighbour(self):
        """Get the Signal-to-Noise Ratio (SNR) for a strip
        candidate, using the C-function (and therefore using the
        value contained in the data member of the C++ class)
        """
        return self._lib.aa_snr_neighbour_getter(self._sensor_analysis)
    @snr_neighbour.setter
    def snr_neighbour(self,snr):
        """Set the Signal-to-Noise Ratio (SNR) for a strip candidate
        to be included in a cluster

        Parameters
        ----------
        snr: float
            The SNR
        """
        self._lib.aa_configure_snr_neighbour(self._sensor_analysis,snr)
    
    # ==============================================
    # Main methods
    def print_configuration(self):
        """Print the current configuration of the AlibavaSensorAnalysis
        class
        """
        self._lib.aa_print_configuration(self._sensor_analysis,self._ioft)

    def configure(self,**cfg):
        """Configure the ...
        
        Parameters
        ----------
        """
        for key,val in cfg.iteritems():
            if hasattr(self,key):
                setattr(self,key,val)
        #self.print_configuration()
    
    def initialize(self):
        """Wrapper to the initialization of the data and
        the algorithms. The noise, pedestal and calibration
        data are filled and the masked channels are set.
        """
        # Initialize the data
        self._lib.ioft_initialize(self._ioft)
        # And mask the noisy channels 
        self._lib.aa_mask_channels(self._sensor_analysis,self._ioft)
        # Print the configuration now
        self.print_configuration()

    def book_results(self,filename):
        """Wrapper to the storage manager to store the analized
        data: book and prepare the output ROOT file

        Parameters
        ----------
        filename: str
            The name of the output ROOT file
        """
        # instantiate the results data member and book the tree
        # XXX need status flag??
        self._results = self._lib.ioresults_new(filename)
        self._lib.ioresults_book_tree(self._results)
        self._store_results = True
        self._results_filename = filename

    def update_crosstalk_factors(self):
        """Wrapper to the update of the cross-talk factors
        This function should be called after the event loop
        XXX: check that this is true? In the 

        Return
        ------
        int: Whether the update was performed without errors
            (0) or not (-1)
        """
        # Update the factors
        status = self._lib.aa_update_crosstalk_factors(self._sensor_analysis)
        # Reset the ioresult 
        # XXX is there a better approach than deleting and reintializating
        if self._store_results:
            self._lib.ioresults_delete(self._results)
            self.book_results(self._results_filename)
            self._lib.ioresults_fill_header(self._results,self._ioft)

    def process(self,i=-1):
        """Processing the event, i.e finding cluster algorithm 
        in action.
        XXX DOC

        Parameters
        ----------
        i: int, optional
            The entry to process, if -1 means all events will
            be processed
        """
        import sys

        # If there is a results instance booked
        if self._store_results:
            # the header
            self._lib.ioresults_fill_header(self._results,self._ioft)
            analyze_event=lambda: self._lib.ioresults_fill_tree(self._results,self._ioft,self._sensor_analysis)
        else:
            analyze_event=lambda: self._lib.aa_find_clusters(self._sensor_analysis,self._ioft)
        # Whether or not process all the events
        if i==-1:
            evts = xrange(self.get_entries())
        else:
            evts = [ i ]
            # Force not to do the cross-talk correction
            self.xtcorrection = 0
        
        for xt_iter in xrange(self.xtcorrection+1):
            print "\r\033[1;34mINFO::::::::\033[1;m Cross-talk correction ITERATION:{0}".format(xt_iter)
            point = float(self.get_entries())/100.0
            for k in evts:
                # Progress bar 
                sys.stdout.write("\r\033[1;34mINFO\033[1;m Alibava cluster analysis "+\
                        "[ "+"\b"+str(int(float(k+1)/point)).rjust(3)+"%]")
                sys.stdout.flush()
                self.process_event(k)
                analyze_event()
            if xt_iter < self.xtcorrection:
                self.update_crosstalk_factors()
            print
    
