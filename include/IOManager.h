// ***********************************************
// 
// Output manager to persistify ROOT files
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#ifndef IOMANAGER_H
#define IOMANAGER_H

#include "input_options.h"

#include <string>
#include <vector>
#include <map>

// forward declarations
class TFile;
class TTree;
class TObject;

class AlibavaRunHeader;
class AlibavaEvent;
class CalibrationParameters;
class AlibavaDiagnosis;

using PedestalNoiseBeetleMap = std::map<int,std::pair<std::vector<float>,std::vector<float>> >;
using CalibrateBeetleMap = std::map<int,std::vector<int>>;


class IOManager
{
    public:
        IOManager(const std::string & rootfilename);
        ~IOManager();

        // Set the calibration parameters (if needed)
        void set_calibration_parameters(const std::vector<int> & param_v);
        // And get them
        inline CalibrationParameters * get_calibration_parameters() const { return _cal_parameters; }

        void book_tree_header();
        void book_tree();
        void book_monitor_plots();
        
        void fill_header(const AlibavaRunHeader * aheader) const;
        void fill_event(const AlibavaEvent * anAlibavaEvent);

        // Diagnostic plots filling functions. A set of predefined 
        // histograms are filled by using the next functions
        // - using only beam data, 
        void set_diagnostic_plots(bool is_pedestal_file_present);
        // - using beam and pedestal data
        void set_diagnostic_plots(const PedestalNoiseBeetleMap & pednoise_m);
        // - histograms getter
        //template<typename ROOTTYPE> 
        //    ROOTTYPE* get_diagnostic_plot(const std::string & plotname, const int & chip);
        // - book histograms on runtime, created and fill externally
        void book_monitor_plot(const std::string & name, const TObject * theplot, const int & chip);
        // - Fill the predefined monitor histograms 
        template <class T1,class T2>
            void update_diagnostic_plot(const int & chip, const std::string & plotname, const T1 & x, const T2 & y);
        // - Fill the predefined monitor histograms (chip-independent version)
        template <class T1,class T2>
            void update_diagnostic_plot(const std::string & plotname, const T1 & x, const T2 & y);
        // - set the 3dimensional histogram for the calibration curves
        //   using the information on the calibration manager
        void set_calibration_plot(const IOManager & calibration_manager);
        // - get all the objects needed to get the 3dim calibration curves
        const std::vector<TObject*> get_calibration_objects(const int & chipnumber) const;
        // - get the pedestal and noise vectors (per beetle) by using the header
        //   note if the header is corrupted, it will use a default values
        //   XXX -- PRIVATE?
        //PedestalNoiseBeetleMap get_pednoise_from_header(const float & ped_default = 512.0,const float & noise_default=3.0);
        
        // Update the runheader tree with the electron per ADC conversion
        void update(const CalibrateBeetleMap & calib_m);

        // Update the Tree (Runheader one) with pedestal/noise info
        // and the Events tree with the signal pedestal-noise subtracted
        void update(const PedestalNoiseBeetleMap & pednoise_m);
        
        void close();
    
        // Event Tree related functions
        // Events tree getter
        TTree * get_events_tree() const;
        // Speed up access to the tree
        void set_events_tree_access(const std::vector<std::string> & branch_list) const;
        // Re-activate everything (and reset branch addresses)
        void reset_events_tree() const;
        // Probably a template, but just vector of floats so far
        void set_events_tree_branch_address(const std::string & branch, std::vector<float> ** v) const;
        void set_events_tree_branch_address(const std::string & branch, float * v ) const;
        // Number of entries and Get the entry
        int get_events_number_entries() const;
        void get_events_entry(const int & i) const;
    
    private:
        // datamembers
        std::string _rootfilename;
        TFile * _file;
        TTree * _tree_header;
        TTree * _tree_events;
        int     _eventsProcessed;

        // The auxiliary data
        // Whether or not the monitor has been booked
        bool _monitor_plots_booked;
        CalibrationParameters * _cal_parameters;
        AlibavaRunHeader * _runheader;
        AlibavaEvent * _events;

        // The monitor/diagnostic plots
        std::map<int,AlibavaDiagnosis*> _monitor_plots;
        // Related methods
        // - get the pedestal and noise vectors (per beetle) by using the header
        //   note if the header is corrupted, it will use a default values
        PedestalNoiseBeetleMap get_pednoise_from_header(const float & ped_default = 512.0,const float & noise_default=3.0);
        // - get the map of signal free 
        
        // store the friends
        void aux_store_friends(TTree * tree);
        // Resurrecting the events tree
        void resurrect_events_tree();
};

#endif
