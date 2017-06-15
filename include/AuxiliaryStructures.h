// ***********************************************
// 
// Auxiliary data structures
//
// **********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#ifndef AUXILIARYSTRUCTURES_H
#define AUXILIARYSTRUCTURES_H

#include "ALIBAVA.h"

// Define alibava header and events (not using ROOT dictionaries
// although it can be think to use it later)
struct AlibavaRunHeader
{
    AlibavaRunHeader() 
    {
        version = -1;
        data_type = -1;
        run_number = -1;
        header = "";
        date_time = "";
        chipSelection.push_back(0);
        chipSelection.push_back(1);
    }
    int version;
    int data_type;
    int run_number;
    std::string header;
    std::string date_time;
    std::vector<int> chipSelection;
    std::vector<float> header_pedestal;
    std::vector<float> header_noise;
};

struct AlibavaEvent
{
    AlibavaEvent() 
    { 
        runNumber = -1;
        eventNumber = -1;
        eventType   = 9999;
        eventSize   = 9999;
        eventClock  = 9999;
        eventTime = -1.0;
        eventTemp = -1.0;
        calCharge = -9999;
        calDelay  = -9999;

        beetle1_chipheader.reserve(ALIBAVA::NOOFCHANNELS);
        beetle2_chipheader.reserve(ALIBAVA::NOOFCHANNELS);

        beetle1_data.reserve(ALIBAVA::NOOFCHANNELS);
        beetle2_data.reserve(ALIBAVA::NOOFCHANNELS);
    }
    unsigned int eventType ;
    unsigned int eventSize ;
    unsigned int eventClock;
    int runNumber;
    int eventNumber;
    float eventTime ;
    float eventTemp ;
    float calCharge ;
    float calDelay  ;

    std::vector<float> beetle1_chipheader;
    std::vector<float> beetle2_chipheader;
    
    std::vector<float> beetle1_data;
    std::vector<float> beetle2_data;
};

// Auxiliary structure with the needed parameters to calibrate
struct CalibrationParameters
{
    CalibrationParameters()
    {
        nPulses = -1;
        initialCharge = -99999;
        finalCharge   = -99999;
        deltaCharge   = -1;
        nSamplesPerPulse = -1;
    }

    int get_injected_charge(const int & evtNumber)
    {
        // WARNING: XXX ---> maybe use std::modf function??
        return initialCharge+deltaCharge*((evtNumber)/(nSamplesPerPulse));
    }

    int nPulses;
    int initialCharge;
    int finalCharge;
    int deltaCharge;
    int nSamplesPerPulse;
};

#endif
