// ***********************************************
// 
// Helper structure useful for gather all input 
// options to the executables. 
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************

#ifndef INPUT_OPTIONS_H
#define INPUT_OPTIONS_H

#include <string>

struct input_options
{
    input_options() 
    {
        // Defaults
        outputFilename = "fortythieves.root";
        pedestal_file = "";
        calibration_file = "";
        storeHeaderPedestalNoise = false;
        runNumber = -1;
        cmndfile  = "";
        startEventNum = -1;
        stopEventNum  = -1;
    }
    std::string cmndfile;
    std::string outputFilename;
    std::string pedestal_file;
    std::string calibration_file;
    int runNumber;
    int startEventNum;
    int stopEventNum;
    bool storeHeaderPedestalNoise;
};

#endif
