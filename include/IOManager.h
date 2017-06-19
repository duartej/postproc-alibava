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

class AlibavaRunHeader;
class AlibavaEvent;
class CalibrationParameters;

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
        
        void fill_header(const AlibavaRunHeader * aheader) const;
        void fill_event(const AlibavaEvent * anAlibavaEvent);
        
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
        CalibrationParameters * _cal_parameters;
        AlibavaRunHeader * _runheader;
        AlibavaEvent * _events;
        
        // store the friends
        void aux_store_friends(TTree * tree);
        // Resurrecting the events tree
        void resurrect_events_tree();
};

#endif