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

// System headers
#include <string>
#include <map>
#include <vector>
#include <memory>

// forward
class IOFortythieves;
class TH2F;

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
        inline void configure_not_automasking() { _automasking = false; }
        inline void configure_mask_criterium(const float & sigma) { _mask_criterium = sigma; }
        inline void configure_snr_seed(const float & snr_min) { _snr_seed = snr_min; }
        inline void configure_snr_neighbour(const float & snr_min) { _snr_neighbour = snr_min; }

        // print configuration 
        void print_configuration(IOFortythieves * ioft = nullptr);

        // Getters 
        inline int get_polarity() const { return _polarity; }
        inline const std::vector<float> & get_time_cut() const { return _tdc_cut; }
        inline const std::vector<int> * get_masked_channels() const { return _masked_channels; }
        inline float get_mask_criterium() const { return _mask_criterium; }
        inline float get_snr_seed() const { return _snr_seed; }
        inline float get_snr_neighbour() const { return _snr_neighbour; }
        
        // Automatic noisy channel masking
        void mask_channels(const IOFortythieves * ioft_inst);

        // Masking event whenever a cut is not fulfill
        bool check_analysis_cuts(const IOFortythieves * ioft_inst) const;

        // Finding algorithm
        std::vector<std::unique_ptr<StripCluster> > find_clusters(const IOFortythieves * ioft_inst);

        // Obtain whetere a channel number is masked or not
        bool is_channel_masked(const int & ich);

        // Calculate the eta of a seed (only using their two neighbours)
        float calculate_seed_eta(const int & seed_channel, const std::vector<float> & adc_data);

        // Asymmetric cross-talk correction
        int update_crosstalk_factors();


    private:
        // Signal polarity of the sensor
        int _polarity;

        // TDC time cut
        std::vector<float> _tdc_cut;

        // Allow automasking
        bool _automasking;
        // Masked channels
        std::vector<int> * _masked_channels;
        // Masked criterium (if automatic masking is enabled)
        // how many sigmas away a channel has to be considered noisy
        float _mask_criterium; 
        
        // Seed and neighbour cuts for the cluster finding
        float _snr_seed;
        float _snr_neighbour;
        
        // Asymmetric cross-talk related data-members
        // Cross-talk counte, the number of times the xt-correction 
        // method has been called
        int _current_xt_iteration;
        // Neighbours to consider in the correction
        int _nNeighbours;
        // The neighbour factors: element=neighbour order+1
        std::vector<float> _xtfactors;
        // The histogram to calculate the cross-talk factors
        // The histogram is storing the charge difference between
        // the i(-bin) neighbour of the cluster seed
        TH2F * _nb_charge_diff_h;

        // Method to correct the asymmetric cross-talk
        std::vector<float> get_data_xt_corrected(const std::vector<float> & original_data);
};


#endif
