// ***********************************************
// 
// Pool of functions to post-process and analyse 
// the ROOT data created by the fortythieves exec
// (https://github.com/duartej/postproc-alibava/blob/master/bin/fortythieves.cc)
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#ifndef ALIBAVASENSORANALYSIS_H
#define ALIBAVASENSORANALYSIS_H

#include "StripCluster.h"

// ROOT

// System headers
#include <string>
#include <map>
#include <vector>

// forward
class IOFortythieves;

class AlibavaSensorAnalysis
{
    public:
        AlibavaSensorAnalysis(IOFortythieves * io_ft);
        AlibavaSensorAnalysis() = delete;
        ~AlibavaSensorAnalysis() { ; }

        // Configuration functions
        inline void configure_time_cut(const float & t0,const float & t1) { _tdc_t0 = t0; _tdc_t1 = t1; }
        inline void configure_masked_channels(const std::vector<int> & masked_channels) 
                                { _masked_channels = new std::vector<int>(masked_channels); }
        inline void configure_snr_seed(const float & snr_min) { _snr_seed = snr_min; }
        inline void configure_snr_neighbour(const float & snr_min) { _snr_neighbour = snr_min; }
        
        // Finding algorithm
        StripCluster find_clusters(const std::vector<float> & adc_corrected);


    private:
        // The input data
        IOFortythieves * _ft;

        // TDC time cut
        float _tdc_t0;
        float _tdc_t1;
        // Masked channels
        std::vector<int> * _masked_channels;
        // Seed and neighbour cuts for the cluster finding
        float _snr_seed;
        float _snr_neighbour;
};


#endif
