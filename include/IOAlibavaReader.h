// ************************************************
// 
// Manager to read binary raw data from the ALIBAVA
// DAQs
//
// ************************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ************************************************
#ifndef IOALIBAVAREADER_H
#define IOALIBAVAREADER_H

#include "input_options.h"
#include "IOManager.h"

#include <string>

class IOAlibavaReader
{
    public:
        // FIXME: Possibly delete
        IOAlibavaReader();
        
        // Read the file defined in the opt, and store it using the IOManager
        static int read_data(const input_options & opt,const IOManager & );

        // Get TDC time
        static double tdc_time(const unsigned int & tdcTime);
        static double get_temperature(const unsigned short & temp);
};

#endif
