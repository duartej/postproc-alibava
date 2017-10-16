// ***************************************************************
// 
// Manager to read the ROOT output for the fortythieve
// exec. https://github.com/duartej/postproc-alibava/blob/master/bin/fortythieves.cc
//
// Entry point for the post-processing and cluster 
// analysis of the ROOT created by fortythieves. 
//
// This class should be updated every time the ROOT 
// trees structure or branch names from the fortythieves
// are modified.
//
// ***************************************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***************************************************************
#ifndef IOFORTYTHIEVES_H
#define IOFORTYTHIEVES_H

#include <string>
#include <map>
#include <list>
#include <vector>

// forward declaraations
class TFile;
class TTree;

class IOFortythieves
{
    public:
        IOFortythieves() = delete;
        IOFortythieves(const std::string & filename,const int & chip);
        ~IOFortythieves();
    
        // Prepare all the variables and fill the run header related vectors
        // (noise, calibration, ...) ready to be used.
        void initialize();

        // Some informational getters
        const char * get_filename() const;
        inline int get_chip_number() const { return _chip; }

        // Process the event entry
        void process(const int & i);
        // Get the total entries for the Events tree
        inline int get_entries() const { return _entries; }

        // Getter for the noise vector 
        inline const std::vector<float> & noise() const { return _noise; }
        inline const std::vector<float> & pedestal() const { return _pedestal; }
        inline const std::vector<float> & calibration() const { return _calibration; }

        // Getters for event-related 
        inline int event_number() const { return _event_number; }
        inline int run_number() const { return _run_number; }
        inline float event_time() const { return _event_time; }
        inline float temperature() const { return _temperature; }
        inline float common_mode() const { return _common_mode; }
        inline float event_noise() const { return _event_noise; }
        inline const std::vector<float> & adc_data() const { return *_adc_data; }

    private:
        // The current beetle
        int _chip;

        // The total number of entries in the Events tree
        int _entries; 

        // Data members
        TFile * _file;
        // A helper list to 
        static std::list<std::string> _TREENAMES;
        // The trees 
        std::map<std::string,TTree*> _trees;

        // The accessors attached to the trees
        // -- the event time, temperature, run and event numbers
        int _run_number;  // Check if it's worth it to use it instead of the one in the runHeader
        int _event_number;
        float _event_time;
        float _temperature;
        // -- Common mode, i.e. the ADC mean in an event (not included those
        //    strips which could be signal)
        float _common_mode;
        // -- Event noise, i.e the standard deviation of the ADC counts in an
        //    event (not included those which could be signal
        float _event_noise;
        // -- the data (in ADC counts, pedestal and common mode extracted)
        std::vector<float> * _adc_data;
        // The noise, pedestal and calibration vectors
        std::vector<float> _noise;
        std::vector<float> _pedestal;
        // Note that this vector provides the number of electrons
        // which corresponds to a ADC (per channel)
        std::vector<float> _calibration;
};


#endif
