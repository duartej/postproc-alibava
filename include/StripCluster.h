// ***********************************************
// 
// Cluster definition  (heavly inspired from the 
// AlibavaCluster class at EUTelescope)
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#ifndef STRIPCLUSTER_H
#define STRIPCLUSTER_H

// ROOT

// System headers
#include <map>
#include <vector>

class StripCluster
{
    public:
        StripCluster();
        ~StripCluster();

        inline void set_polarity(const int & polarity) { _signal_polarity = polarity; }
        inline void set_eta_seed(const float & etaseed) { _eta_seed = etaseed; }
        inline void set_sensitive_direction(const int & direction) { _sensitive_direction = direction; }
        void add(const int & channel, const float & adc);
        // Getters
        // -- Number of strips belonging to the cluster
        inline unsigned int size() const { return _channels.size(); }
        // -- The number of channels given the element index 
        //    Note that the seed is always the element 0
        inline int channels(const int & element_index) { return _channels[element_index]; }
        // -- Get the leftest channel in the cluster
        int left_channel() const;
        // -- Get the rightest channel in the cluster
        int right_channel() const;
        // -- The charge of the channel 
        inline float charge(const int & element_index) { return _signal_polarity*_signal_map[_element_index[element_index]]; }
        // -- The total charge of the cluster
        float charge() const;
        // The calibrated charge
        float calibrated_charge(const std::vector<float> & calibration); 
        // The cluster signal to noise
        float snr(const std::vector<float> & noise) ;
        // -- The eta of the cluster (see details in the implementation)
        float eta();
        // -- The eta of the seed (see details in the implementation)
        inline float eta_seed() const { return _eta_seed; }
        // -- The position (in channel coordinates) estimated with 
        //    the Center of Gravity algorithm
        float position_cog();

    private:
        int _signal_polarity;
        // The channel number ordered by charge,
        // therefore the seed is always at zero position
        std::vector<int> _channels;
        // Map relating the channel number with its charge
        std::map<int,float> _signal_map;
        // Map to related the index in the _channel vector 
        // with the channel number
        std::map<unsigned int, int> _element_index;
        
        // 0: X, 1:Y and 2:Z
        int _sensitive_direction;

        bool _eta_calculated;
        float _eta_seed;
        float _eta_cluster;
};

#endif
