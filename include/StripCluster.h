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
        inline void set_eta_seed(const float & eta_seed) { _eta_seed = eta_seed; }
        inline void set_sensitive_direction(const int & direction) { _sensitive_direction = direction; }
        void add(const int & channel, const float & adc);

    private:
        int _signal_polarity;
        // The channel number ordered by charge,
        // therefore the seed is always at zero position
        std::vector<int> _channels;
        std::map<int,float> _signal_map;
        
        // 0: X, 1:Y and 2:Z
        int _sensitive_direction;

        float _eta_seed;
        float _eta_cluster;
};

#endif
