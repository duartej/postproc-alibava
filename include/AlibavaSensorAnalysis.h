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

// XXX
#include "IOFortythieves.h"

// ROOT

// System headers
#include <string>
#include <map>
#include <vector>
#include <memory>

// forward
class IOFortythieves;

class AlibavaSensorAnalysis
{
    public:
        AlibavaSensorAnalysis();
        ~AlibavaSensorAnalysis();

        // Configuration functions
        inline void configure_polarity(const int & polarity) { _polarity = polarity; }
        inline void configure_time_cut(const float & t0,const float & t1) { _tdc_cut[0] = t0; _tdc_cut[1] = t1; }
        inline void configure_masked_channels(const std::vector<int> & masked_channels) 
        { 
            if(_masked_channels != nullptr)
            {
                _masked_channels->resize(0);
                _masked_channels->assign(masked_channels.begin(),masked_channels.end());
            }
            else
            {
                _masked_channels = new std::vector<int>(masked_channels); 
            }
        }
        inline void configure_mask_criterium(const float & sigma) { _mask_criterium = sigma; }
        inline void configure_snr_seed(const float & snr_min) { _snr_seed = snr_min; }
        inline void configure_snr_neighbour(const float & snr_min) { _snr_neighbour = snr_min; }

        // Getters 
        inline int get_polarity() const { return _polarity; }
        inline const std::vector<float> & get_time_cut() const { return _tdc_cut; }
        inline const std::vector<int> * get_masked_channels() const { return _masked_channels; }
        inline float get_mask_criterium() const { return _mask_criterium; }
        inline float get_snr_seed() const { return _snr_seed; }
        inline float get_snr_neighbour() const { return _snr_neighbour; }
        // Per event base getters, wrappers to the IOFortythieves instance
        /*inline int get_event_number() { return _ft->event_number(); }
        inline int get_run_number() { return _ft->run_number(); }
        inline int event_time() { return _ft->event_time(); }
        inline int get_temperature() { return _ft->temperature(); }*/
        
        // Automatic noisy channel masking
        void mask_channels(const IOFortythieves * ioft_inst);

        // Finding algorithm
        std::vector<std::unique_ptr<StripCluster> > find_clusters(const IOFortythieves * ioft_inst);

        // Obtain whetere a channel number is masked or not
        bool is_channel_masked(const int & ich);


    private:
        // Signal polarity of the sensor
        int _polarity;

        // TDC time cut
        std::vector<float> _tdc_cut;

        // Masked channels
        std::vector<int> * _masked_channels;
        // Masked criterium (if automatic masking is enabled)
        // how many sigmas away a channel has to be considered noisy
        float _mask_criterium; 
        
        // Seed and neighbour cuts for the cluster finding
        float _snr_seed;
        float _snr_neighbour;
};


#endif
