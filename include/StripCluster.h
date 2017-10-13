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

    private:
        int _signal_polarity;
        // The channel number ordered by charge,
        // therefore the seed is always at zero position
        std::vector<int> _channels;
        std::map<int,float> _signal_map;
};

#endif
